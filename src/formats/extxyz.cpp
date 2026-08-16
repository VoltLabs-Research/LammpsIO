// Ported from OVITO's XYZImporter (MIT option of its dual license); extended-XYZ support

#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>
#include <lammpsio/frame.hpp>

namespace lammpsio {
namespace extxyz {

constexpr const char* SPECIES_NAMES[] = { "species", "type", "element", "atom_types" };

struct ColumnSpec {
    std::string name;
    int index = 0;
};

struct Layout {
    ColumnMapping cols;
    std::vector<std::string> headers;
    bool symbolicSpecies = true;
};

struct HeaderScan {
    FrameHeader header;
    Layout layout;
    const char* atomsSection = nullptr;
    bool valid = false;
};

bool isSpeciesName(const std::string& name) {
    for (const char* candidate : SPECIES_NAMES) {
        if (name == candidate) return true;
    }
    return false;
}

const char* findKey(const char* content, const char* lineEnd, const char* key, size_t keyLength) {
    for (const char* p = content; p + keyLength <= lineEnd; p++) {
        size_t matched = 0;
        while (matched < keyLength) {
            char left = p[matched];
            char right = key[matched];
            if (left >= 'A' && left <= 'Z') left += 32;
            if (right >= 'A' && right <= 'Z') right += 32;
            if (left != right) break;
            matched++;
        }
        if (matched == keyLength) return p + keyLength;
    }
    return nullptr;
}

void readKeyValue(const char* start, const char* lineEnd, const char*& valueStart, const char*& valueEnd) {
    if (start < lineEnd && (*start == '"' || *start == '\'')) {
        const char quote = *start;
        valueStart = start + 1;
        valueEnd = valueStart;
        while (valueEnd < lineEnd && *valueEnd != quote) valueEnd++;
    } else {
        valueStart = start;
        valueEnd = findTokenEnd(start, lineEnd);
    }
}

bool parseLattice(const char* content, const char* lineEnd, FrameHeader& header) {
    const char* after = findKey(content, lineEnd, "lattice=", 8);
    if (!after) return false;

    const char* valueStart;
    const char* valueEnd;
    readKeyValue(after, lineEnd, valueStart, valueEnd);

    double values[9] = { 0 };
    const char* token = skipWhitespace(valueStart, valueEnd);
    int count = 0;

    while (token < valueEnd && count < 9) {
        const char* tokenEnd = findTokenEnd(token, valueEnd);
        values[count++] = fastAtof(token, tokenEnd);
        token = skipWhitespace(tokenEnd, valueEnd);
    }

    if (count != 9) return false;

    for (int vector = 0; vector < 3; vector++) {
        for (int axis = 0; axis < 3; axis++) {
            header.cell[vector][axis] = values[vector * 3 + axis];
        }
    }
    deriveBoxFromCell(header);
    return true;
}

void parsePbc(const char* content, const char* lineEnd, FrameHeader& header, bool hasLattice) {
    for (int axis = 0; axis < 3; axis++) header.periodic[axis] = hasLattice;

    const char* after = findKey(content, lineEnd, "pbc=", 4);
    if (!after) return;

    const char* valueStart;
    const char* valueEnd;
    readKeyValue(after, lineEnd, valueStart, valueEnd);

    const char* token = skipWhitespace(valueStart, valueEnd);
    for (int axis = 0; axis < 3 && token < valueEnd; axis++) {
        header.periodic[axis] = *token == 'T' || *token == 't' || *token == '1';
        token = skipWhitespace(findTokenEnd(token, valueEnd), valueEnd);
    }
}

bool parseProperties(const char* content, const char* lineEnd, Layout& layout) {
    const char* after = findKey(content, lineEnd, "properties=", 11);
    if (!after) return false;

    const char* valueStart;
    const char* valueEnd;
    readKeyValue(after, lineEnd, valueStart, valueEnd);

    std::vector<std::string> fields;
    const char* p = valueStart;
    while (p < valueEnd) {
        const char* separator = p;
        while (separator < valueEnd && *separator != ':') separator++;
        fields.emplace_back(p, separator - p);
        p = separator + 1;
    }

    if (fields.size() < 3) return false;

    int column = 0;
    for (size_t triple = 0; triple + 2 < fields.size(); triple += 3) {
        std::string name = fields[triple];
        for (char& c : name) if (c >= 'A' && c <= 'Z') c += 32;
        const char type = fields[triple + 1].empty() ? 'R' : fields[triple + 1][0];
        const int columns = fastAtoi(fields[triple + 2].c_str(),
                                     fields[triple + 2].c_str() + fields[triple + 2].size());

        for (int component = 0; component < (columns > 0 ? columns : 1); component++) {
            if (name == "pos" && component < 3) {
                if (component == 0) layout.cols.idxX = column;
                else if (component == 1) layout.cols.idxY = column;
                else layout.cols.idxZ = column;
                layout.headers.push_back(component == 0 ? "x" : (component == 1 ? "y" : "z"));
            } else if (isSpeciesName(name) && columns == 1) {
                layout.cols.idxType = column;
                layout.symbolicSpecies = type == 'S' || type == 's';
                layout.headers.push_back("type");
            } else if (name == "id" && columns == 1) {
                layout.cols.idxId = column;
                layout.headers.push_back("id");
            } else {
                layout.headers.push_back(columns > 1 ? name + "_" + std::to_string(component) : name);
            }
            column++;
        }
    }

    layout.cols.computeMaxIdx();
    return layout.cols.idxX >= 0;
}

void assumePlainLayout(const char* firstRow, const char* end, Layout& layout) {
    const char* lineEnd = findLineEnd(firstRow, end);
    const char* token = skipWhitespace(firstRow, lineEnd);
    int columns = 0;
    while (token < lineEnd) {
        columns++;
        token = skipWhitespace(findTokenEnd(token, lineEnd), lineEnd);
    }

    layout.cols.idxType = 0;
    layout.cols.idxX = 1;
    layout.cols.idxY = 2;
    layout.cols.idxZ = 3;
    layout.symbolicSpecies = true;
    layout.headers = { "type", "x", "y", "z" };
    for (int column = 4; column < columns; column++) {
        layout.headers.push_back("column_" + std::to_string(column));
    }
    layout.cols.computeMaxIdx();
}

HeaderScan parseHeader(const char* RESTRICT data, const char* RESTRICT end) {
    HeaderScan scan;

    const char* countLine = skipWhitespace(data, end);
    const char* countEnd = findLineEnd(countLine, end);
    const char* countTokenEnd = findTokenEnd(countLine, countEnd);
    if (countTokenEnd == countLine || !isIntegerToken(countLine, countTokenEnd)) return scan;

    scan.header.atomCount = fastAtoi(countLine, countTokenEnd);
    if (scan.header.atomCount < 0) return scan;

    const char* commentLine = countEnd + 1;
    if (commentLine > end) return scan;
    const char* commentEnd = findLineEnd(commentLine, end);

    const bool hasLattice = parseLattice(commentLine, commentEnd, scan.header);
    parsePbc(commentLine, commentEnd, scan.header, hasLattice);

    scan.atomsSection = commentEnd < end ? commentEnd + 1 : end;

    if (!parseProperties(commentLine, commentEnd, scan.layout)) {
        assumePlainLayout(scan.atomsSection, end, scan.layout);
    }

    scan.header.headers = scan.layout.headers;
    scan.valid = scan.layout.cols.idxX >= 0 && scan.layout.cols.idxY >= 0 && scan.layout.cols.idxZ >= 0;
    return scan;
}

bool sniff(const MappedFile& file) {
    const char* end = file.data + file.size;
    const char* p = skipWhitespace(file.data, end);
    if (p >= end) return false;

    const char* lineEnd = findLineEnd(p, end);
    const char* tokenEnd = findTokenEnd(p, lineEnd);
    if (tokenEnd == p || !isIntegerToken(p, tokenEnd)) return false;
    if (skipWhitespace(tokenEnd, lineEnd) < lineEnd) return false;

    return lineEnd < end && parseHeader(file.data, end).valid;
}

bool scan(const MappedFile& file, std::vector<FrameIndexEntry>& frames, std::string& error) {
    const char* end = file.data + file.size;
    const char* p = file.data;

    while (p < end) {
        const char* content = skipWhitespace(p, end);
        if (content >= end) break;

        HeaderScan header = parseHeader(content, end);
        if (!header.valid) {
            if (frames.empty()) {
                error = "Invalid XYZ format";
                return false;
            }
            break;
        }

        FrameIndexEntry entry;
        entry.index = (int)frames.size();
        entry.byteOffset = (size_t)(content - file.data);
        entry.timestep = entry.index;
        entry.atomCount = header.header.atomCount;

        const char* cursor = header.atomsSection;
        for (int atom = 0; atom < header.header.atomCount && cursor < end; atom++) {
            cursor = jumpToNextLine(cursor, end);
        }

        entry.byteLength = (size_t)(cursor - file.data) - entry.byteOffset;
        frames.push_back(entry);
        p = cursor;
    }

    if (frames.empty()) {
        error = "Invalid XYZ format";
        return false;
    }

    return true;
}

bool readHeader(const MappedFile& file, const FrameIndexEntry& entry,
                FrameHeader& header, std::string& error) {
    HeaderScan scanned = parseHeader(file.data + entry.byteOffset, frameEndPointer(file, entry));
    if (!scanned.valid) {
        error = "Invalid XYZ format";
        return false;
    }

    header = std::move(scanned.header);
    header.timestep = entry.timestep;
    return true;
}

bool readFrame(const MappedFile& file, const FrameIndexEntry& entry, const ReadOptions& options,
               FrameAllocator& allocator, ParsedFrame& frame, std::string& error) {
    const char* frameEnd = frameEndPointer(file, entry);
    HeaderScan header = parseHeader(file.data + entry.byteOffset, frameEnd);
    if (!header.valid) {
        error = "Invalid XYZ format";
        return false;
    }

    ColumnMapping cols = header.layout.cols;
    resolveExtraColumns(options.properties, header.layout.headers, cols);

    frame.header = header.header;
    frame.header.timestep = entry.timestep;
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

    std::unordered_map<std::string, int> speciesIds;

    frame.bbox.init();
    const char* p = header.atomsSection;
    const int maxColumn = cols.maxIdx;
    int atomIndex = 0;

    while (p < frameEnd && atomIndex < atomCount) {
        const char* lineEnd = findLineEnd(p, frameEnd);
        const char* content = skipWhitespace(p, lineEnd);

        if (UNLIKELY(content >= lineEnd)) {
            p = lineEnd + 1;
            continue;
        }

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
                if (header.layout.symbolicSpecies && !isIntegerToken(token, tokenEnd)) {
                    const std::string symbol(token, tokenEnd - token);
                    auto existing = speciesIds.find(symbol);
                    if (existing != speciesIds.end()) {
                        type = existing->second;
                    } else {
                        type = (int)speciesIds.size() + 1;
                        speciesIds.emplace(symbol, type);
                        frame.elementHintsByType.push_back(symbol);
                    }
                } else {
                    type = fastAtoi(token, tokenEnd);
                }
            } else if (buffers.ids && column == cols.idxId) {
                id = (uint32_t)fastAtoi(token, tokenEnd);
            }

            for (int extra = 0; extra < extraCount; extra++) {
                if (column != cols.extraPropIndices[extra]) continue;
                const double value = fastAtof(token, tokenEnd);
                frame.extras[extra].values[(size_t)atomIndex] = value;
                if (!isIntegerToken(token, tokenEnd) ||
                    value < -2147483648.0 || value > 2147483647.0) {
                    frame.extras[extra].dtype = ColumnDtype::Float32;
                }
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

    if (frame.header.cell[0][0] == 0.0 && frame.header.cell[1][1] == 0.0 &&
        frame.header.cell[2][2] == 0.0 && atomIndex > 0) {
        frame.header.cell[0][0] = frame.bbox.maxX - frame.bbox.minX;
        frame.header.cell[1][1] = frame.bbox.maxY - frame.bbox.minY;
        frame.header.cell[2][2] = frame.bbox.maxZ - frame.bbox.minZ;
        frame.header.origin[0] = frame.bbox.minX;
        frame.header.origin[1] = frame.bbox.minY;
        frame.header.origin[2] = frame.bbox.minZ;
        deriveBoxFromCell(frame.header);
    }

    return true;
}

extern const FormatReader reader = {
    format_id::ExtXyz,
    sniff,
    scan,
    readHeader,
    readFrame
};

}
}
