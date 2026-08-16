#include <cstring>
#include <string>
#include <vector>
#include <lammpsio/frame.hpp>

namespace lammpsio {
namespace lammps_data {

constexpr size_t HEADER_SCAN_LIMIT = 8192;

struct AtomStyleLayout {
    const char* name;
    int idxId;
    int idxType;
    int idxX;
};

constexpr AtomStyleLayout ATOM_STYLES[] = {
    { "atomic",       0,   1,   2 },
    { "charge",       0,   1,   3 },
    { "full",         0,   2,   4 },
    { "molecular",    0,   2,   3 },
    { "bond",         0,   2,   3 },
    { "angle",        0,   2,   3 },
    { "dipole",       0,   1,   3 },
    { "sphere",       0,   1,   4 },
    { "peri",         0,   1,   4 },
    { "electron",     0,   1,   5 },
    { "meso",         0,   1,   5 },
    { "spin",         0,   1,   2 },
    { "template",     0,   4,   5 },
    { "body",         0,   1,   5 },
    { "ellipsoid",    0,   1,   5 },
    { "line",         0,   3,   6 },
    { "tri",          0,   3,   6 }
};

struct HeaderScan {
    FrameHeader header;
    bool valid = false;
};

bool isSectionHeader(const char* content, const char* lineEnd, const char* word, size_t length) {
    if ((size_t)(lineEnd - content) < length || strncmp(content, word, length) != 0) return false;
    const char* after = skipWhitespace(content + length, lineEnd);
    return after >= lineEnd || *after == '#';
}

HeaderScan parseHeader(const char* RESTRICT data, size_t size) {
    HeaderScan scan;
    const char* limit = data + (size < HEADER_SCAN_LIMIT ? size : HEADER_SCAN_LIMIT);
    const char* p = data;
    uint8_t found = 0;

    while (p < limit) {
        const char* lineEnd = findLineEnd(p, limit);
        const char* content = skipWhitespace(p, lineEnd);

        if (UNLIKELY(content >= lineEnd || *content == '#')) {
            p = lineEnd + 1;
            continue;
        }

        if (!(found & 1)) {
            const char* tokenEnd = findTokenEnd(content, lineEnd);
            const char* keyword = skipWhitespace(tokenEnd, lineEnd);
            if (keyword + 5 <= lineEnd && strncmp(keyword, "atoms", 5) == 0) {
                scan.header.atomCount = fastAtoi(content, tokenEnd);
                found |= 1;
                p = lineEnd + 1;
                continue;
            }
        }

        {
            double values[3] = { 0.0, 0.0, 0.0 };
            const char* labels[3] = { nullptr, nullptr, nullptr };
            const char* token = content;
            int index = 0;

            for (; index < 6 && token < lineEnd; index++) {
                const char* tokenEnd = findTokenEnd(token, lineEnd);
                if (index < 3) values[index] = fastAtof(token, tokenEnd);
                else labels[index - 3] = token;
                token = skipWhitespace(tokenEnd, lineEnd);
            }

            if (index == 6 && token >= lineEnd &&
                labels[0][0] == 'x' && labels[0][1] == 'y' &&
                labels[1][0] == 'x' && labels[1][1] == 'z' &&
                labels[2][0] == 'y' && labels[2][1] == 'z') {
                scan.header.box.xy = values[0];
                scan.header.box.xz = values[1];
                scan.header.box.yz = values[2];
                p = lineEnd + 1;
                continue;
            }
        }

        if (lineEnd - content > 4) {
            const char* loPos = (const char*)memmem(content, lineEnd - content, "lo", 2);
            if (loPos && loPos > content) {
                const char axis = *(loPos - 1);

                const char* token = content;
                const char* tokenEnd = findTokenEnd(token, lineEnd);
                const double lo = fastAtof(token, tokenEnd);

                token = skipWhitespace(tokenEnd, lineEnd);
                tokenEnd = findTokenEnd(token, lineEnd);
                const double hi = fastAtof(token, tokenEnd);

                if (axis == 'x' && !(found & 2)) {
                    scan.header.box.xlo = lo;
                    scan.header.box.xhi = hi;
                    found |= 2;
                } else if (axis == 'y' && !(found & 4)) {
                    scan.header.box.ylo = lo;
                    scan.header.box.yhi = hi;
                    found |= 4;
                } else if (axis == 'z' && !(found & 8)) {
                    scan.header.box.zlo = lo;
                    scan.header.box.zhi = hi;
                    found |= 8;
                }
            }
        }

        p = lineEnd + 1;
    }

    scan.valid = found == 15;
    if (scan.valid) deriveCellFromBox(scan.header);
    return scan;
}

void parseMasses(const char* RESTRICT data, size_t size, ParsedFrame& frame) {
    const char* end = data + size;
    const char* p = data;
    const char* section = nullptr;

    while (p < end) {
        const char* lineEnd = findLineEnd(p, end);
        const char* content = skipWhitespace(p, lineEnd);
        if (isSectionHeader(content, lineEnd, "Masses", 6)) {
            section = lineEnd + 1;
            break;
        }
        p = lineEnd + 1;
    }
    if (!section) return;

    p = section;
    while (p < end) {
        const char* lineEnd = findLineEnd(p, end);
        const char* content = skipWhitespace(p, lineEnd);
        if (content < lineEnd && *content != '#') break;
        p = lineEnd + 1;
    }

    struct Row {
        int type;
        double mass;
        std::string hint;
    };
    std::vector<Row> rows;
    int maxType = 0;

    while (p < end) {
        const char* lineEnd = findLineEnd(p, end);
        const char* content = skipWhitespace(p, lineEnd);

        if (UNLIKELY(content >= lineEnd)) break;
        if (*content == '#') { p = lineEnd + 1; continue; }
        if (*content >= 'A' && *content <= 'Z') break;

        const char* tokenEnd = findTokenEnd(content, lineEnd);
        const int type = fastAtoi(content, tokenEnd);

        const char* massToken = skipWhitespace(tokenEnd, lineEnd);
        const char* massEnd = findTokenEnd(massToken, lineEnd);
        const double mass = fastAtof(massToken, massEnd);

        std::string hint;
        const char* hash = (const char*)memchr(massEnd, '#', lineEnd - massEnd);
        if (hash) {
            const char* hintStart = skipWhitespace(hash + 1, lineEnd);
            const char* hintEnd = findTokenEnd(hintStart, lineEnd);
            if (hintEnd > hintStart) hint.assign(hintStart, hintEnd - hintStart);
        }

        if (type > 0) {
            rows.push_back({ type, mass, std::move(hint) });
            if (type > maxType) maxType = type;
        }
        p = lineEnd + 1;
    }

    if (maxType == 0) return;

    frame.massesByType.assign((size_t)maxType, 0.0);
    frame.elementHintsByType.assign((size_t)maxType, std::string());
    for (const auto& row : rows) {
        frame.massesByType[(size_t)row.type - 1] = row.mass;
        frame.elementHintsByType[(size_t)row.type - 1] = row.hint;
    }
}

const char* findAtomsSection(const char* data, size_t size, std::string& style) {
    const char* end = data + size;
    const char* p = data;

    while (p < end) {
        const char* lineEnd = findLineEnd(p, end);
        const char* content = skipWhitespace(p, lineEnd);

        if (isSectionHeader(content, lineEnd, "Atoms", 5)) {
            const char* hash = (const char*)memchr(content, '#', lineEnd - content);
            if (hash) {
                const char* styleStart = skipWhitespace(hash + 1, lineEnd);
                const char* styleEnd = findTokenEnd(styleStart, lineEnd);
                if (styleEnd > styleStart) style.assign(styleStart, styleEnd - styleStart);
            }
            return lineEnd + 1;
        }

        p = lineEnd + 1;
    }

    return nullptr;
}

bool applyAtomStyle(const std::string& style, ColumnMapping& cols) {
    if (style.empty()) return false;

    for (const auto& layout : ATOM_STYLES) {
        if (style != layout.name) continue;
        cols.idxId = layout.idxId;
        cols.idxType = layout.idxType;
        cols.idxX = layout.idxX;
        cols.idxY = layout.idxX + 1;
        cols.idxZ = layout.idxX + 2;
        cols.computeMaxIdx();
        return true;
    }

    return false;
}

int countColumns(const char* p, const char* end) {
    const char* lineEnd = findLineEnd(p, end);
    const char* token = skipWhitespace(p, lineEnd);
    int columns = 0;

    while (token < lineEnd) {
        columns++;
        token = skipWhitespace(findTokenEnd(token, lineEnd), lineEnd);
    }

    return columns;
}

bool sniff(const MappedFile& file) {
    const HeaderScan scan = parseHeader(file.data, file.size);
    return scan.valid;
}

bool scan(const MappedFile& file, std::vector<FrameIndexEntry>& frames, std::string& error) {
    const HeaderScan scanned = parseHeader(file.data, file.size);
    if (!scanned.valid) {
        error = "Invalid LAMMPS data format";
        return false;
    }

    FrameIndexEntry entry;
    entry.index = 0;
    entry.byteOffset = 0;
    entry.byteLength = file.size;
    entry.timestep = 0;
    entry.atomCount = scanned.header.atomCount;
    frames.push_back(entry);
    return true;
}

bool readHeader(const MappedFile& file, const FrameIndexEntry&,
                FrameHeader& header, std::string& error) {
    HeaderScan scanned = parseHeader(file.data, file.size);
    if (!scanned.valid) {
        error = "Invalid LAMMPS data format";
        return false;
    }

    header = std::move(scanned.header);
    return true;
}

bool readFrame(const MappedFile& file, const FrameIndexEntry&, const ReadOptions& options,
               FrameAllocator& allocator, ParsedFrame& frame, std::string& error) {
    HeaderScan scanned = parseHeader(file.data, file.size);
    if (!scanned.valid) {
        error = "Invalid LAMMPS data format";
        return false;
    }

    std::string style;
    const char* p = findAtomsSection(file.data, file.size, style);
    if (!p) {
        error = "LAMMPS data file has no Atoms section";
        return false;
    }

    const char* end = file.data + file.size;
    while (p < end) {
        const char* content = skipWhitespace(p, end);
        if (content < end && *content != '\n' && *content != '\r' && *content != '#') break;
        p = jumpToNextLine(p, end);
    }

    ColumnMapping cols;
    if (!applyAtomStyle(style, cols)) {
        detectDataColumnStyle(countColumns(p, end), cols);
    }

    frame.header = std::move(scanned.header);
    frame.hasIds = options.includeIds && cols.idxId >= 0;
    parseMasses(file.data, file.size, frame);

    const int expectedAtoms = frame.header.atomCount;
    FrameBuffers buffers = allocator.allocate(expectedAtoms, frame.hasIds);

    frame.bbox.init();
    const int maxColumn = cols.maxIdx;
    int atomIndex = 0;

    while (p < end && atomIndex < expectedAtoms) {
        const char* lineEnd = findLineEnd(p, end);
        const char* content = skipWhitespace(p, lineEnd);

        if (UNLIKELY(content >= lineEnd || *content == '#')) {
            p = lineEnd + 1;
            continue;
        }
        if (UNLIKELY(*content >= 'A' && *content <= 'Z')) break;

        double x = 0, y = 0, z = 0;
        int type = 0;
        uint32_t id = 0;

        const char* token = content;
        int column = 0;

        while (token < lineEnd && column <= maxColumn) {
            const char* tokenEnd = findTokenEnd(token, lineEnd);

            if (column == cols.idxX) {
                x = fastAtof(token, tokenEnd);
            } else if (column == cols.idxY) {
                y = fastAtof(token, tokenEnd);
            } else if (column == cols.idxZ) {
                z = fastAtof(token, tokenEnd);
            } else if (column == cols.idxType) {
                type = fastAtoi(token, tokenEnd);
            } else if (buffers.ids && column == cols.idxId) {
                id = (uint32_t)fastAtoi(token, tokenEnd);
            }

            token = skipWhitespace(tokenEnd, lineEnd);
            column++;
        }

        buffers.setPosition(atomIndex, x, y, z);
        buffers.types[atomIndex] = (uint16_t)type;
        if (buffers.ids) buffers.ids[atomIndex] = id;

        frame.bbox.update((float)x, (float)y, (float)z);
        atomIndex++;
        p = lineEnd + 1;
    }

    return true;
}

extern const FormatReader reader = {
    format_id::LammpsData,
    sniff,
    scan,
    readHeader,
    readFrame
};

}
}
