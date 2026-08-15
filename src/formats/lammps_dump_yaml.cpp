// Ported from OVITO's LAMMPSDumpYAMLImporter (MIT option of its dual license).

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <lammpsio/frame.hpp>

namespace lammpsio {
namespace lammps_dump_yaml {

namespace {

struct HeaderScan {
    FrameHeader header;
    ColumnMapping cols;
    bool scaledCoords = false;
    const char* dataSection = nullptr;
    bool valid = false;
};

bool startsWith(const char* content, const char* lineEnd, const char* prefix) {
    const size_t length = std::strlen(prefix);
    return (size_t)(lineEnd - content) >= length && std::strncmp(content, prefix, length) == 0;
}

bool isDocumentStart(const char* content, const char* lineEnd) {
    if (!startsWith(content, lineEnd, "---")) return false;
    return skipWhitespace(content + 3, lineEnd) >= lineEnd;
}

void splitFlowSequence(const char* content, const char* lineEnd, std::vector<std::string>& items) {
    const char* open = (const char*)std::memchr(content, '[', lineEnd - content);
    if (!open) return;

    const char* close = lineEnd;
    while (close > open && *(close - 1) != ']') close--;
    if (close <= open) close = lineEnd;
    else close--;

    const char* p = open + 1;
    while (p < close) {
        const char* comma = p;
        while (comma < close && *comma != ',') comma++;

        const char* itemStart = skipWhitespace(p, comma);
        const char* itemEnd = comma;
        while (itemEnd > itemStart && (unsigned char)*(itemEnd - 1) <= ' ') itemEnd--;
        if (itemEnd > itemStart) items.emplace_back(itemStart, itemEnd - itemStart);

        p = comma + 1;
    }
}

const char* valueAfterKey(const char* content, const char* lineEnd, size_t keyLength) {
    return skipWhitespace(content + keyLength, lineEnd);
}

void applyKeywords(const std::vector<std::string>& keywords, HeaderScan& scan) {
    for (size_t index = 0; index < keywords.size(); index++) {
        std::string name = keywords[index];
        for (char& c : name) if (c >= 'A' && c <= 'Z') c += 32;

        if (name == "type") {
            scan.cols.idxType = (int)index;
        } else if (name == "id") {
            scan.cols.idxId = (int)index;
        } else if (!name.empty() && (name[0] == 'x' || name[0] == 'y' || name[0] == 'z') &&
                   isPositionColumn(name.c_str(), name.size())) {
            if (name.size() >= 2 && name[1] == 's') scan.scaledCoords = true;
            if (name[0] == 'x') scan.cols.idxX = (int)index;
            else if (name[0] == 'y') scan.cols.idxY = (int)index;
            else scan.cols.idxZ = (int)index;
        }

        scan.header.headers.push_back(std::move(name));
    }

    scan.cols.computeMaxIdx();
}

HeaderScan parseHeader(const char* RESTRICT data, const char* RESTRICT end) {
    HeaderScan scan;
    const char* p = data;
    bool boxSeen = false;
    double tilt[3] = { 0.0, 0.0, 0.0 };
    double bounds[3][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } };
    int boxRow = 0;
    bool inBox = false;

    while (p < end) {
        const char* lineEnd = findLineEnd(p, end);
        const char* content = skipWhitespace(p, lineEnd);

        if (content >= lineEnd) {
            p = lineEnd + 1;
            continue;
        }

        if (p != data && isDocumentStart(content, lineEnd)) break;
        if (startsWith(content, lineEnd, "...")) break;

        if (inBox) {
            if (*content == '-') {
                std::vector<std::string> row;
                splitFlowSequence(content, lineEnd, row);
                if (boxRow < 3 && row.size() >= 2) {
                    bounds[boxRow][0] = fastAtof(row[0].c_str(), row[0].c_str() + row[0].size());
                    bounds[boxRow][1] = fastAtof(row[1].c_str(), row[1].c_str() + row[1].size());
                } else if (boxRow == 3 && row.size() >= 3) {
                    for (int index = 0; index < 3; index++) {
                        tilt[index] = fastAtof(row[index].c_str(), row[index].c_str() + row[index].size());
                    }
                }
                boxRow++;
                p = lineEnd + 1;
                continue;
            }
            inBox = false;
        }

        if (startsWith(content, lineEnd, "timestep:")) {
            const char* value = valueAfterKey(content, lineEnd, 9);
            scan.header.timestep = fastAtoi(value, findTokenEnd(value, lineEnd));
        } else if (startsWith(content, lineEnd, "natoms:")) {
            const char* value = valueAfterKey(content, lineEnd, 7);
            scan.header.atomCount = fastAtoi(value, findTokenEnd(value, lineEnd));
        } else if (startsWith(content, lineEnd, "boundary:")) {
            std::vector<std::string> flags;
            splitFlowSequence(content, lineEnd, flags);
            if (flags.size() == 6) {
                for (int axis = 0; axis < 3; axis++) {
                    scan.header.periodic[axis] = flags[axis * 2] == "p" && flags[axis * 2 + 1] == "p";
                }
            }
        } else if (startsWith(content, lineEnd, "keywords:")) {
            std::vector<std::string> keywords;
            splitFlowSequence(content, lineEnd, keywords);
            applyKeywords(keywords, scan);
        } else if (startsWith(content, lineEnd, "box:")) {
            inBox = true;
            boxSeen = true;
            boxRow = 0;
        } else if (startsWith(content, lineEnd, "data:")) {
            scan.dataSection = lineEnd + 1;
            break;
        }

        p = lineEnd + 1;
    }

    if (boxSeen) {
        scan.header.box.xlo = bounds[0][0]; scan.header.box.xhi = bounds[0][1];
        scan.header.box.ylo = bounds[1][0]; scan.header.box.yhi = bounds[1][1];
        scan.header.box.zlo = bounds[2][0]; scan.header.box.zhi = bounds[2][1];
        scan.header.box.xy = tilt[0];
        scan.header.box.xz = tilt[1];
        scan.header.box.yz = tilt[2];

        const double minX = std::min(std::min(0.0, tilt[0]), std::min(tilt[1], tilt[0] + tilt[1]));
        const double maxX = std::max(std::max(0.0, tilt[0]), std::max(tilt[1], tilt[0] + tilt[1]));
        scan.header.box.xlo -= minX;
        scan.header.box.xhi -= maxX;
        scan.header.box.ylo -= std::min(0.0, tilt[2]);
        scan.header.box.yhi -= std::max(0.0, tilt[2]);
        deriveCellFromBox(scan.header);
    }

    scan.valid = scan.dataSection != nullptr && boxSeen &&
                 scan.cols.idxX >= 0 && scan.cols.idxY >= 0 && scan.cols.idxZ >= 0;
    return scan;
}

const char* frameEndPointer(const MappedFile& file, const FrameIndexEntry& entry) {
    const size_t end = entry.byteOffset + entry.byteLength;
    return file.data + (end > file.size ? file.size : end);
}

}

bool sniff(const MappedFile& file) {
    const char* end = file.data + file.size;
    const char* p = skipWhitespace(file.data, end);
    if (p >= end) return false;

    const char* lineEnd = findLineEnd(p, end);
    if (!isDocumentStart(p, lineEnd)) return false;
    if (lineEnd >= end) return false;

    p = skipWhitespace(lineEnd + 1, end);
    lineEnd = findLineEnd(p, end);
    return startsWith(p, lineEnd, "creator: LAMMPS");
}

bool scan(const MappedFile& file, std::vector<FrameIndexEntry>& frames, std::string& error) {
    const char* end = file.data + file.size;
    const char* p = file.data;

    while (p < end) {
        const char* content = skipWhitespace(p, end);
        if (content >= end) break;

        const char* lineEnd = findLineEnd(content, end);
        if (!isDocumentStart(content, lineEnd)) break;

        HeaderScan header = parseHeader(content, end);
        if (!header.valid) {
            if (frames.empty()) {
                error = "Invalid LAMMPS YAML dump format";
                return false;
            }
            break;
        }

        FrameIndexEntry entry;
        entry.index = (int)frames.size();
        entry.byteOffset = (size_t)(content - file.data);
        entry.timestep = header.header.timestep;
        entry.atomCount = header.header.atomCount;

        const char* cursor = header.dataSection;
        for (int atom = 0; atom < header.header.atomCount && cursor < end; atom++) {
            cursor = jumpToNextLine(cursor, end);
        }
        while (cursor < end) {
            const char* rowEnd = findLineEnd(cursor, end);
            const char* rowStart = skipWhitespace(cursor, rowEnd);
            if (rowStart < rowEnd && !startsWith(rowStart, rowEnd, "...")) break;
            cursor = rowEnd < end ? rowEnd + 1 : end;
            if (rowStart < rowEnd) break;
        }

        entry.byteLength = (size_t)(cursor - file.data) - entry.byteOffset;
        frames.push_back(entry);
        p = cursor;
    }

    if (frames.empty()) {
        error = "Invalid LAMMPS YAML dump format";
        return false;
    }

    return true;
}

bool readHeader(const MappedFile& file, const FrameIndexEntry& entry,
                FrameHeader& header, std::string& error) {
    HeaderScan scanned = parseHeader(file.data + entry.byteOffset, frameEndPointer(file, entry));
    if (!scanned.valid) {
        error = "Invalid LAMMPS YAML dump format";
        return false;
    }

    header = std::move(scanned.header);
    return true;
}

bool readFrame(const MappedFile& file, const FrameIndexEntry& entry, const ReadOptions& options,
               FrameAllocator& allocator, ParsedFrame& frame, std::string& error) {
    const char* frameEnd = frameEndPointer(file, entry);
    HeaderScan header = parseHeader(file.data + entry.byteOffset, frameEnd);
    if (!header.valid) {
        error = "Invalid LAMMPS YAML dump format";
        return false;
    }

    ColumnMapping cols = header.cols;
    resolveExtraColumns(options.properties, header.header.headers, cols);

    frame.header = header.header;
    frame.hasIds = options.includeIds && cols.idxId >= 0;

    const int atomCount = header.header.atomCount;
    FrameBuffers buffers = allocator.allocate(atomCount, frame.hasIds);

    const int extraCount = (int)cols.extraPropIndices.size();
    frame.extras.resize(extraCount);
    for (int extra = 0; extra < extraCount; extra++) {
        frame.extras[extra].name = cols.extraPropNames[extra];
        frame.extras[extra].values.resize((size_t)atomCount);
        frame.extras[extra].dtype = ColumnDtype::Int32;
    }

    const double lx = header.header.box.xhi - header.header.box.xlo;
    const double ly = header.header.box.yhi - header.header.box.ylo;
    const double lz = header.header.box.zhi - header.header.box.zlo;

    frame.bbox.init();
    const char* p = header.dataSection;
    int atomIndex = 0;

    while (p < frameEnd && atomIndex < atomCount) {
        const char* lineEnd = findLineEnd(p, frameEnd);
        const char* content = skipWhitespace(p, lineEnd);

        if (content >= lineEnd) {
            p = lineEnd + 1;
            continue;
        }
        if (*content != '-') break;

        std::vector<std::string> row;
        splitFlowSequence(content, lineEnd, row);
        if (row.empty()) {
            p = lineEnd + 1;
            continue;
        }

        const auto field = [&row](int index) -> const std::string* {
            return index >= 0 && index < (int)row.size() ? &row[(size_t)index] : nullptr;
        };
        const auto asDouble = [](const std::string* value) {
            return value ? fastAtof(value->c_str(), value->c_str() + value->size()) : 0.0;
        };

        double x = asDouble(field(cols.idxX));
        double y = asDouble(field(cols.idxY));
        double z = asDouble(field(cols.idxZ));

        if (header.scaledCoords) {
            const double fx = x, fy = y, fz = z;
            x = header.header.box.xlo + fx * lx + fy * header.header.box.xy + fz * header.header.box.xz;
            y = header.header.box.ylo + fy * ly + fz * header.header.box.yz;
            z = header.header.box.zlo + fz * lz;
        }

        buffers.setPosition(atomIndex, x, y, z);
        buffers.types[atomIndex] = (uint16_t)asDouble(field(cols.idxType));
        if (buffers.ids) buffers.ids[atomIndex] = (uint32_t)asDouble(field(cols.idxId));

        for (int extra = 0; extra < extraCount; extra++) {
            const std::string* value = field(cols.extraPropIndices[extra]);
            if (!value) continue;
            const char* start = value->c_str();
            const char* stop = start + value->size();
            frame.extras[extra].values[(size_t)atomIndex] = fastAtof(start, stop);
            if (!isIntegerToken(start, stop)) frame.extras[extra].dtype = ColumnDtype::Float32;
        }

        frame.bbox.update((float)x, (float)y, (float)z);
        atomIndex++;
        p = lineEnd + 1;
    }

    return true;
}

extern const FormatReader reader = {
    format_id::LammpsDumpYaml,
    sniff,
    scan,
    readHeader,
    readFrame
};

}
}
