// LAMMPS text dump reader (`ITEM:`-delimited), the format `dump ... atom/custom` writes.
//
// Multi-frame aware: a dump file is a concatenation of frames, and scan() reports the
// byte range of each so a caller can split one file into single-frame files without
// reparsing. The previous version of this parser read the first frame and silently
// ignored the rest.

#include <algorithm>
#include <cstring>
#include <future>
#include <string>
#include <thread>
#include <vector>
#include "../frame.hpp"

namespace lammps_dump_text {

namespace {

constexpr int MULTITHREAD_ATOM_THRESHOLD = 50000;

/** True for a line that opens an `ITEM:` section. */
ALWAYS_INLINE bool isItemLine(const char* content, const char* lineEnd) {
    return lineEnd - content >= 5 && content[0] == 'I' && content[4] == ':';
}

/** Points past `ITEM: ` on a line already known to be a section header. */
ALWAYS_INLINE const char* itemLabel(const char* content) {
    return content + 6;
}

struct HeaderScan {
    FrameHeader header;
    const char* atomsSection = nullptr;
    bool valid = false;
};

/**
 * Reads the four `ITEM:` sections that precede the atom rows.
 *
 * Triclinic dumps print a third value per bounds line (the tilt factors) and their
 * lo/hi are the *bounding box* inflated by those tilts, so the true edges are recovered
 * below. For orthogonal cells the tilts are zero and the recovery is a no-op.
 */
HeaderScan parseHeader(const char* RESTRICT data, const char* RESTRICT end, ColumnMapping& cols) {
    HeaderScan scan;
    const char* p = data;
    uint8_t found = 0; // 1=timestep, 2=natoms, 4=bounds, 8=atoms

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
                // The rest of this line carries the boundary styles, optionally preceded
                // by the literal `xy xz yz` that marks a triclinic dump.
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
                        // x/xu are Cartesian; xs/xsu are fractions of the box and need
                        // mapping through it. The 's' right after the axis letter marks
                        // the scaled styles.
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
                scan.atomsSection = lineEnd + 1;
                found |= 8;
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
    /**
     * Per-extra-column: set when this chunk saw a token that is not integer-formatted.
     * Merged across chunks, because a column is i32 only if no chunk disagreed.
     */
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

        // The next frame's header ends this one.
        if (UNLIKELY(isItemLine(content, lineEnd))) break;

        float x = 0, y = 0, z = 0;
        int type = 0;
        uint32_t id = 0;

        const char* token = content;
        int column = 0;

        while (token < lineEnd && column <= maxColumn) {
            const char* tokenEnd = findTokenEnd(token, lineEnd);

            if (column == cols.idxX) {
                x = (float)fastAtof(token, tokenEnd);
            } else if (column == cols.idxY) {
                y = (float)fastAtof(token, tokenEnd);
            } else if (column == cols.idxZ) {
                z = (float)fastAtof(token, tokenEnd);
            } else if (column == cols.idxType) {
                type = fastAtoi(token, tokenEnd);
            } else if (buffers.ids && column == cols.idxId) {
                id = (uint32_t)fastAtoi(token, tokenEnd);
            }

            for (int extra = 0; extra < extraCount; extra++) {
                if (column != cols.extraPropIndices[extra]) continue;
                const double value = fastAtof(token, tokenEnd);
                extras[extra][atomIndex] = value;
                // A '.'/exponent, or an integer past int32 range (LAMMPS ids and mol
                // numbers can exceed it), downgrades the whole column to f32.
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
            x = (float)(box.xlo + fx * lx + fy * box.xy + fz * box.xz);
            y = (float)(box.ylo + fy * ly + fz * box.yz);
            z = (float)(box.zlo + fz * lz);
        }

        const int base = atomIndex * 3;
        buffers.positions[base] = x;
        buffers.positions[base + 1] = y;
        buffers.positions[base + 2] = z;
        buffers.types[atomIndex] = (uint16_t)type;
        if (buffers.ids) buffers.ids[atomIndex] = id;

        bbox.update(x, y, z);
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

const char* frameEndPointer(const MappedFile& file, const FrameIndexEntry& entry) {
    const size_t end = entry.byteOffset + entry.byteLength;
    return file.data + (end > file.size ? file.size : end);
}

} // namespace

bool sniff(const MappedFile& file) {
    // The first non-blank line of a dump is always `ITEM: TIMESTEP`.
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
            // Not a frame boundary: a truncated or malformed tail, nothing more to index.
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

        // Skipping exactly natoms lines is much cheaper than testing every line for an
        // `ITEM:` prefix, and it lands on the next frame's first line.
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
    unsigned int threadCount = std::thread::hardware_concurrency();
    if (threadCount == 0) threadCount = 1;
    if (atomCount < MULTITHREAD_ATOM_THRESHOLD) threadCount = 1;

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

        // Each chunk needs to know where its atoms start in the output buffers, which
        // means counting the earlier chunks' rows before any of them can write.
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

// `extern` is load-bearing: a const object at namespace scope defaults to internal
// linkage, so without it the registry's declaration finds no definition to link to.
extern const FormatReader reader = {
    format_id::LammpsDumpText,
    sniff,
    scan,
    readHeader,
    readFrame
};

} // namespace lammps_dump_text
