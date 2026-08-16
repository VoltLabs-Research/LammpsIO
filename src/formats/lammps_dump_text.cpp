#include <algorithm>
#include <cstring>
#include <future>
#include <string>
#include <thread>
#include <vector>
#include <lammpsio/frame.hpp>

namespace lammpsio {
namespace lammps_dump_text {

constexpr int MULTITHREAD_ATOM_THRESHOLD = 50000;

ALWAYS_INLINE bool isItemLine(const char* content, const char* lineEnd) {
    return lineEnd - content >= 5 && content[0] == 'I' && content[4] == ':';
}

ALWAYS_INLINE const char* itemLabel(const char* content) {
    return content + 6;
}

struct HeaderScan {
    FrameHeader header;
    const char* atomsSection = nullptr;
    bool valid = false;
};

HeaderScan parseHeader(const char* RESTRICT data, const char* RESTRICT end, ColumnMapping& cols) {
    HeaderScan scan;
    const char* p = data;
    uint8_t found = 0;

    while (p < end && found != 15) {
        const char* lineEnd = findLineEnd(p, end);
        const char* content = skipWhitespace(p, lineEnd);

        if (UNLIKELY(content >= lineEnd)) {
            p = lineEnd + 1;
            continue;
        }

        if (isItemLine(content, lineEnd)) {
            const char* label = itemLabel(content);

            if (!(found & 1) && strncmp(label, "TIMESTEP", 8) == 0) {
                p = lineEnd + 1;
                const char* valueEnd = findLineEnd(p, end);
                scan.header.timestep = fastAtoi(skipWhitespace(p, valueEnd), valueEnd);
                found |= 1;
                p = valueEnd;
            } else if (!(found & 2) && strncmp(label, "NUMBER OF ATOMS", 15) == 0) {
                p = lineEnd + 1;
                const char* valueEnd = findLineEnd(p, end);
                scan.header.atomCount = fastAtoi(skipWhitespace(p, valueEnd), valueEnd);
                found |= 2;
                p = valueEnd;
            } else if (!(found & 4) && strncmp(label, "BOX BOUNDS", 10) == 0) {
                const char* flag = skipWhitespace(label + 10, lineEnd);
                if (flag + 1 < lineEnd && flag[0] == 'x' && flag[1] == 'y') {
                    for (int skipped = 0; skipped < 3 && flag < lineEnd; skipped++) {
                        flag = skipWhitespace(findTokenEnd(flag, lineEnd), lineEnd);
                    }
                }
                for (int axis = 0; axis < 3 && flag < lineEnd; axis++) {
                    scan.header.periodic[axis] = flag[0] == 'p';
                    flag = skipWhitespace(findTokenEnd(flag, lineEnd), lineEnd);
                }

                double tiltXY = 0.0, tiltXZ = 0.0, tiltYZ = 0.0;
                for (int axis = 0; axis < 3 && p < end; axis++) {
                    p = lineEnd + 1;
                    lineEnd = findLineEnd(p, end);

                    const char* token = skipWhitespace(p, lineEnd);
                    const char* tokenEnd = findTokenEnd(token, lineEnd);
                    const double lo = fastAtof(token, tokenEnd);

                    token = skipWhitespace(tokenEnd, lineEnd);
                    tokenEnd = findTokenEnd(token, lineEnd);
                    const double hi = fastAtof(token, tokenEnd);

                    token = skipWhitespace(tokenEnd, lineEnd);
                    double tilt = 0.0;
                    if (token < lineEnd) {
                        tokenEnd = findTokenEnd(token, lineEnd);
                        tilt = fastAtof(token, tokenEnd);
                    }

                    if (axis == 0) { scan.header.box.xlo = lo; scan.header.box.xhi = hi; tiltXY = tilt; }
                    else if (axis == 1) { scan.header.box.ylo = lo; scan.header.box.yhi = hi; tiltXZ = tilt; }
                    else { scan.header.box.zlo = lo; scan.header.box.zhi = hi; tiltYZ = tilt; }
                }

                scan.header.box.xy = tiltXY;
                scan.header.box.xz = tiltXZ;
                scan.header.box.yz = tiltYZ;

                const double minX = std::min(std::min(0.0, tiltXY), std::min(tiltXZ, tiltXY + tiltXZ));
                const double maxX = std::max(std::max(0.0, tiltXY), std::max(tiltXZ, tiltXY + tiltXZ));
                scan.header.box.xlo -= minX;
                scan.header.box.xhi -= maxX;
                scan.header.box.ylo -= std::min(0.0, tiltYZ);
                scan.header.box.yhi -= std::max(0.0, tiltYZ);
                deriveCellFromBox(scan.header);
                found |= 4;
            } else if (!(found & 8) && strncmp(label, "ATOMS", 5) == 0) {
                const char* head = skipWhitespace(label + 5, lineEnd);
                int columnIndex = 0;

                while (head < lineEnd) {
                    const char* tokenEnd = findTokenEnd(head, lineEnd);
                    const size_t length = tokenEnd - head;
                    const char first = (head[0] >= 'A' && head[0] <= 'Z') ? head[0] + 32 : head[0];

                    if (length == 4 && first == 't' && head[1] == 'y') {
                        cols.idxType = columnIndex;
                    } else if (length == 2 && first == 'i' && head[1] == 'd') {
                        cols.idxId = columnIndex;
                    } else if ((first == 'x' || first == 'y' || first == 'z') && isPositionColumn(head, length)) {
                        if (length >= 2 && head[1] == 's') cols.scaledCoords = true;
                        if (first == 'x') cols.idxX = columnIndex;
                        else if (first == 'y') cols.idxY = columnIndex;
                        else cols.idxZ = columnIndex;
                    }

                    std::string name(head, length);
                    for (char& c : name) if (c >= 'A' && c <= 'Z') c += 32;
                    scan.header.headers.push_back(std::move(name));

                    head = skipWhitespace(tokenEnd, lineEnd);
                    columnIndex++;
                }

                cols.computeMaxIdx();
                scan.header.positionsWereScaled = cols.scaledCoords;
                scan.atomsSection = lineEnd + 1;
                found |= 8;
            } else {
                const char* nameEnd = lineEnd;
                while (nameEnd > label && (unsigned char)*(nameEnd - 1) <= ' ') nameEnd--;
                std::string name(label, nameEnd - label);

                std::string value;
                if (lineEnd < end) {
                    const char* valueLine = lineEnd + 1;
                    const char* valueEnd = findLineEnd(valueLine, end);
                    const char* valueStart = skipWhitespace(valueLine, valueEnd);
                    const char* trimmed = valueEnd;
                    while (trimmed > valueStart && (unsigned char)*(trimmed - 1) <= ' ') trimmed--;
                    value.assign(valueStart, trimmed - valueStart);
                    lineEnd = valueEnd;
                }

                scan.header.extraSections.emplace_back(std::move(name), std::move(value));
            }
        }

        p = lineEnd + 1;
    }

    scan.valid = found == 15 && cols.idxType >= 0 &&
                 cols.idxX >= 0 && cols.idxY >= 0 && cols.idxZ >= 0;
    return scan;
}

struct ChunkResult {
    BoundingBox bbox;
    int count = 0;
    std::vector<uint8_t> nonInteger;
};

HOT void parseChunk(const char* RESTRICT chunkStart, const char* RESTRICT chunkEnd,
                    const char* RESTRICT globalEnd, FrameBuffers buffers, int startIndex,
                    const ColumnMapping& cols, const SimulationBox& box,
                    ChunkResult* result, double** extras, int extraCount) {
    const char* p = chunkStart;
    int atomIndex = startIndex;
    BoundingBox bbox;
    bbox.init();

    result->nonInteger.assign(extraCount, 0);
    uint8_t* RESTRICT nonInteger = extraCount > 0 ? result->nonInteger.data() : nullptr;

    const int maxColumn = cols.maxIdx;
    const bool scaled = cols.scaledCoords;
    const double lx = box.xhi - box.xlo;
    const double ly = box.yhi - box.ylo;
    const double lz = box.zhi - box.zlo;

    while (p < chunkEnd) {
        const char* lineEnd = findLineEnd(p, globalEnd);
        const char* content = skipWhitespace(p, lineEnd);

        if (UNLIKELY(content >= lineEnd)) {
            p = lineEnd + 1;
            continue;
        }

        if (UNLIKELY(isItemLine(content, lineEnd))) break;

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

            for (int extra = 0; extra < extraCount; extra++) {
                if (column != cols.extraPropIndices[extra]) continue;
                const double value = fastAtof(token, tokenEnd);
                extras[extra][atomIndex] = value;
                if (!isIntegerToken(token, tokenEnd) ||
                    value < -2147483648.0 || value > 2147483647.0) {
                    nonInteger[extra] = 1;
                }
            }

            token = skipWhitespace(tokenEnd, lineEnd);
            column++;
        }

        if (scaled) {
            const double fx = x, fy = y, fz = z;
            x = box.xlo + fx * lx + fy * box.xy + fz * box.xz;
            y = box.ylo + fy * ly + fz * box.yz;
            z = box.zlo + fz * lz;
        }

        buffers.setPosition(atomIndex, x, y, z);
        buffers.types[atomIndex] = (uint16_t)type;
        if (buffers.ids) buffers.ids[atomIndex] = id;

        bbox.update((float)x, (float)y, (float)z);
        atomIndex++;
        p = lineEnd + 1;
    }

    result->bbox = bbox;
    result->count = atomIndex - startIndex;
}

int countAtomsInChunk(const char* start, const char* end, const char* globalEnd) {
    int count = 0;
    const char* p = start;

    while (p < end) {
        const char* lineEnd = findLineEnd(p, globalEnd);
        const char* content = skipWhitespace(p, lineEnd);

        if (content < lineEnd) {
            if (UNLIKELY(isItemLine(content, lineEnd))) break;
            count++;
        }
        p = lineEnd + 1;
    }

    return count;
}

bool sniff(const MappedFile& file) {
    const char* end = file.data + file.size;
    const char* p = file.data;

    while (p < end) {
        const char* lineEnd = findLineEnd(p, end);
        const char* content = skipWhitespace(p, lineEnd);
        if (content < lineEnd) {
            return isItemLine(content, lineEnd) && strncmp(itemLabel(content), "TIMESTEP", 8) == 0;
        }
        p = lineEnd + 1;
    }

    return false;
}

bool scan(const MappedFile& file, std::vector<FrameIndexEntry>& frames, std::string& error) {
    const char* end = file.data + file.size;
    const char* p = file.data;

    while (p < end) {
        const char* lineEnd = findLineEnd(p, end);
        const char* content = skipWhitespace(p, lineEnd);

        if (content >= lineEnd) {
            p = lineEnd + 1;
            continue;
        }

        if (!isItemLine(content, lineEnd) || strncmp(itemLabel(content), "TIMESTEP", 8) != 0) {
            break;
        }

        ColumnMapping cols;
        HeaderScan header = parseHeader(p, end, cols);
        if (!header.valid || !header.atomsSection) {
            if (frames.empty()) {
                error = "Invalid LAMMPS dump format";
                return false;
            }
            break;
        }

        FrameIndexEntry entry;
        entry.index = (int)frames.size();
        entry.byteOffset = (size_t)(p - file.data);
        entry.timestep = header.header.timestep;
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
        error = "Invalid LAMMPS dump format";
        return false;
    }

    return true;
}

bool readHeader(const MappedFile& file, const FrameIndexEntry& entry,
                FrameHeader& header, std::string& error) {
    ColumnMapping cols;
    HeaderScan scanned = parseHeader(file.data + entry.byteOffset, frameEndPointer(file, entry), cols);
    if (!scanned.valid) {
        error = "Invalid LAMMPS dump format";
        return false;
    }

    header = std::move(scanned.header);
    return true;
}

bool readFrame(const MappedFile& file, const FrameIndexEntry& entry, const ReadOptions& options,
               FrameAllocator& allocator, ParsedFrame& frame, std::string& error) {
    const char* frameStart = file.data + entry.byteOffset;
    const char* frameEnd = frameEndPointer(file, entry);

    ColumnMapping cols;
    HeaderScan header = parseHeader(frameStart, frameEnd, cols);
    if (!header.valid || !header.atomsSection) {
        error = "Invalid LAMMPS dump format";
        return false;
    }

    resolveExtraColumns(options.properties, header.header.headers, cols);

    frame.header = header.header;
    frame.hasIds = options.includeIds && cols.idxId >= 0;

    const int atomCount = header.header.atomCount;
    FrameBuffers buffers = allocator.allocate(atomCount, frame.hasIds);

    const int extraCount = (int)cols.extraPropIndices.size();
    frame.extras.resize(extraCount);
    std::vector<double*> extraPointers(extraCount);
    for (int extra = 0; extra < extraCount; extra++) {
        frame.extras[extra].name = cols.extraPropNames[extra];
        frame.extras[extra].values.resize((size_t)atomCount);
        extraPointers[extra] = frame.extras[extra].values.data();
    }
    double** extras = extraCount > 0 ? extraPointers.data() : nullptr;

    const char* dataStart = header.atomsSection;
    const unsigned int threadCount =
        resolveThreadCount(options.maxThreads, atomCount, MULTITHREAD_ATOM_THRESHOLD);

    std::vector<ChunkResult> results(threadCount);

    if (threadCount == 1) {
        parseChunk(dataStart, frameEnd, frameEnd, buffers, 0,
                   cols, header.header.box, &results[0], extras, extraCount);
    } else {
        const size_t chunkSize = (size_t)(frameEnd - dataStart) / threadCount;
        std::vector<const char*> bounds(threadCount + 1);
        bounds[0] = dataStart;
        bounds[threadCount] = frameEnd;
        for (unsigned int i = 1; i < threadCount; i++) {
            bounds[i] = jumpToNextLine(dataStart + i * chunkSize, frameEnd);
        }

        std::vector<std::future<int>> counts;
        for (unsigned int i = 0; i < threadCount; i++) {
            counts.push_back(std::async(std::launch::async,
                countAtomsInChunk, bounds[i], bounds[i + 1], frameEnd));
        }

        std::vector<int> offsets(threadCount, 0);
        int running = 0;
        for (unsigned int i = 0; i < threadCount; i++) {
            offsets[i] = running;
            running += counts[i].get();
        }

        std::vector<std::thread> threads;
        for (unsigned int i = 0; i < threadCount; i++) {
            threads.emplace_back(parseChunk, bounds[i], bounds[i + 1], frameEnd, buffers,
                                 offsets[i], std::cref(cols), std::cref(header.header.box),
                                 &results[i], extras, extraCount);
        }
        for (auto& thread : threads) thread.join();
    }

    frame.bbox.init();
    for (const auto& result : results) {
        if (result.count > 0) frame.bbox.merge(result.bbox);
    }

    for (int extra = 0; extra < extraCount; extra++) {
        frame.extras[extra].dtype = ColumnDtype::Int32;
        for (const auto& result : results) {
            if (extra < (int)result.nonInteger.size() && result.nonInteger[extra]) {
                frame.extras[extra].dtype = ColumnDtype::Float32;
                break;
            }
        }
    }

    return true;
}

extern const FormatReader reader = {
    format_id::LammpsDumpText,
    sniff,
    scan,
    readHeader,
    readFrame
};

}
}
