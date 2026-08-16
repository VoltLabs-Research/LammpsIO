// reader is a port of its LAMMPSBinaryDumpImporter (MIT option of its dual license).

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <lammpsio/frame.hpp>

namespace lammpsio {
namespace lammps_dump_binary {

constexpr int32_t ENDIAN_MARKER = 0x0001;
constexpr int32_t FORMAT_REVISION = 0x0002;
constexpr int32_t MAX_SIZE_ONE = 40;
constexpr int64_t MAX_ATOMS = 100000000000LL;

constexpr const char* ATOM_STYLE_COLUMNS = "id type xs ys zs";

struct Cursor {
    const char* p;
    const char* end;
    bool bigEndian;
    bool wideBigInt;

    bool remaining(size_t bytes) const { return (size_t)(end - p) >= bytes; }

    template <typename T>
    bool readRaw(T& out) {
        if (!remaining(sizeof(T))) return false;
        char bytes[sizeof(T)];
        std::memcpy(bytes, p, sizeof(T));
        if (bigEndian) std::reverse(bytes, bytes + sizeof(T));
        std::memcpy(&out, bytes, sizeof(T));
        p += sizeof(T);
        return true;
    }

    bool readInt(int32_t& out) { return readRaw(out); }
    bool readDouble(double& out) { return readRaw(out); }

    bool readBigInt(int64_t& out) {
        if (wideBigInt) return readRaw(out);
        int32_t narrow = 0;
        if (!readRaw(narrow)) return false;
        out = narrow;
        return true;
    }

    bool readBytes(size_t length, std::string& out) {
        if (!remaining(length)) return false;
        out.assign(p, length);
        p += length;
        return true;
    }

    bool skip(size_t length) {
        if (!remaining(length)) return false;
        p += length;
        return true;
    }
};

struct BinaryHeader {
    FrameHeader header;
    int32_t sizeOne = 0;
    int32_t chunkCount = 0;
    std::string columns;
    const char* body = nullptr;
    bool scaledCoords = false;
};

void recoverBoxEdges(SimulationBox& box) {
    const double xy = box.xy, xz = box.xz, yz = box.yz;
    const double minX = std::min(std::min(0.0, xy), std::min(xz, xy + xz));
    const double maxX = std::max(std::max(0.0, xy), std::max(xz, xy + xz));
    box.xlo -= minX;
    box.xhi -= maxX;
    box.ylo -= std::min(0.0, yz);
    box.yhi -= std::max(0.0, yz);
}

void applyColumns(const std::string& columns, BinaryHeader& parsed, ColumnMapping& cols) {
    const char* p = columns.c_str();
    const char* end = p + columns.size();
    int index = 0;

    while (p < end) {
        p = skipWhitespace(p, end);
        if (p >= end) break;
        const char* tokenEnd = findTokenEnd(p, end);

        std::string name(p, tokenEnd - p);
        for (char& c : name) if (c >= 'A' && c <= 'Z') c += 32;

        if (name == "type") {
            cols.idxType = index;
        } else if (name == "id") {
            cols.idxId = index;
        } else if (!name.empty() && (name[0] == 'x' || name[0] == 'y' || name[0] == 'z') &&
                   isPositionColumn(name.c_str(), name.size())) {
            if (name.size() >= 2 && name[1] == 's') parsed.scaledCoords = true;
            if (name[0] == 'x') cols.idxX = index;
            else if (name[0] == 'y') cols.idxY = index;
            else cols.idxZ = index;
        }

        parsed.header.headers.push_back(std::move(name));
        p = tokenEnd;
        index++;
    }

    cols.computeMaxIdx();
}

bool parseFrameHeader(Cursor cursor, BinaryHeader& parsed, ColumnMapping& cols) {
    int64_t timestep = 0;
    if (!cursor.readBigInt(timestep)) return false;

    int32_t revision = 0;
    bool atomStyle = false;

    if (timestep < 0) {
        const int64_t magicLength = -timestep;
        if (magicLength != 8 && magicLength != 10) return false;

        std::string magic;
        if (!cursor.readBytes((size_t)magicLength, magic)) return false;
        if (magic != "DUMPATOM" && magic != "DUMPCUSTOM") return false;
        atomStyle = magic == "DUMPATOM";

        int32_t endianMarker = 0;
        if (!cursor.readInt(endianMarker) || endianMarker != ENDIAN_MARKER) return false;
        if (!cursor.readInt(revision) || revision != FORMAT_REVISION) return false;
        if (!cursor.readBigInt(timestep) || timestep < 0) return false;
    }

    int64_t atomCount = 0;
    if (!cursor.readBigInt(atomCount) || atomCount < 0 || atomCount > MAX_ATOMS) return false;

    const char* beforeBoundary = cursor.p;
    int32_t triclinic = -1;
    if (!cursor.readInt(triclinic)) return false;

    int32_t boundary[3][2] = { { 0, 0 }, { 0, 0 }, { 0, 0 } };
    bool boundaryLooksValid = true;
    for (int axis = 0; axis < 3; axis++) {
        for (int face = 0; face < 2; face++) {
            if (!cursor.readInt(boundary[axis][face])) return false;
            if (boundary[axis][face] < 0 || boundary[axis][face] > 3) boundaryLooksValid = false;
        }
    }

    if (!boundaryLooksValid) {
        cursor.p = beforeBoundary;
        triclinic = -1;
        for (int axis = 0; axis < 3; axis++) {
            boundary[axis][0] = boundary[axis][1] = 0;
        }
    }

    SimulationBox box = {};
    double bounds[3][2];
    for (int axis = 0; axis < 3; axis++) {
        if (!cursor.readDouble(bounds[axis][0]) || !cursor.readDouble(bounds[axis][1])) return false;
        if (!(bounds[axis][1] > bounds[axis][0])) return false;
        for (int face = 0; face < 2; face++) {
            if (!std::isfinite(bounds[axis][face]) ||
                bounds[axis][face] < -1e9 || bounds[axis][face] > 1e9) return false;
        }
    }
    if (triclinic < -1 || triclinic > 1) return false;

    box.xlo = bounds[0][0]; box.xhi = bounds[0][1];
    box.ylo = bounds[1][0]; box.yhi = bounds[1][1];
    box.zlo = bounds[2][0]; box.zhi = bounds[2][1];

    if (triclinic != 0) {
        const char* beforeTilt = cursor.p;
        double tilt[3];
        bool tiltLooksValid = true;
        for (int index = 0; index < 3; index++) {
            if (!cursor.readDouble(tilt[index])) return false;
            if (!std::isfinite(tilt[index])) tiltLooksValid = false;
        }
        if (tiltLooksValid) {
            for (int axis = 0; axis < 3; axis++) {
                const double span = bounds[axis][1] - bounds[axis][0];
                if (tilt[axis] < -span || tilt[axis] > span) tiltLooksValid = false;
            }
        }

        if (tiltLooksValid) {
            box.xy = tilt[0];
            box.xz = tilt[1];
            box.yz = tilt[2];
        } else {
            cursor.p = beforeTilt;
        }
    }

    recoverBoxEdges(box);

    int32_t sizeOne = 0;
    if (!cursor.readInt(sizeOne) || sizeOne <= 0 || sizeOne > MAX_SIZE_ONE) return false;

    std::string columns;
    if (revision >= FORMAT_REVISION) {
        int32_t unitStyleLength = 0;
        if (!cursor.readInt(unitStyleLength) || unitStyleLength < 0) return false;
        if (!cursor.skip((size_t)unitStyleLength)) return false;

        int32_t timeFlag = 0;
        std::string flagByte;
        if (!cursor.readBytes(1, flagByte)) return false;
        timeFlag = (unsigned char)flagByte[0];
        if (timeFlag) {
            double simulationTime = 0;
            if (!cursor.readDouble(simulationTime)) return false;
        }

        int32_t columnsLength = 0;
        if (!cursor.readInt(columnsLength) || columnsLength < 0) return false;
        if (!cursor.readBytes((size_t)columnsLength, columns)) return false;
    }

    int32_t chunkCount = 0;
    if (!cursor.readInt(chunkCount) || chunkCount <= 0 || (int64_t)chunkCount > atomCount) return false;

    parsed.header.timestep = (int)timestep;
    parsed.header.atomCount = (int)atomCount;
    parsed.header.box = box;
    for (int axis = 0; axis < 3; axis++) {
        parsed.header.periodic[axis] = boundary[axis][0] == 0 && boundary[axis][1] == 0;
    }
    deriveCellFromBox(parsed.header);

    parsed.sizeOne = sizeOne;
    parsed.chunkCount = chunkCount;
    parsed.body = cursor.p;

    if (columns.empty()) {
        columns = atomStyle ? ATOM_STYLE_COLUMNS : "id type x y z";
        for (int32_t extra = 5; extra < sizeOne; extra++) {
            columns += " column_" + std::to_string(extra);
        }
    }

    parsed.columns = columns;
    applyColumns(parsed.columns, parsed, cols);

    return cols.idxX >= 0 && cols.idxY >= 0 && cols.idxZ >= 0;
}

bool measureFrameBody(const BinaryHeader& parsed, const char* end, const char*& frameEnd) {
    const char* p = parsed.body;

    for (int32_t chunk = 0; chunk < parsed.chunkCount; chunk++) {
        if ((size_t)(end - p) < sizeof(int32_t)) return false;
        int32_t values = 0;
        std::memcpy(&values, p, sizeof(int32_t));
        p += sizeof(int32_t);
        if (values < 0) return false;
        const size_t bytes = (size_t)values * sizeof(double);
        if ((size_t)(end - p) < bytes) return false;
        p += bytes;
    }

    frameEnd = p;
    return true;
}

struct Convention {
    bool bigEndian;
    bool wideBigInt;
};

constexpr Convention CONVENTIONS[] = {
    { false, true },
    { false, false },
    { true, true },
    { true, false }
};

bool detectConvention(const MappedFile& file, size_t offset, Convention& convention,
                      BinaryHeader& parsed, ColumnMapping& cols) {
    for (const Convention candidate : CONVENTIONS) {
        BinaryHeader attempt;
        ColumnMapping attemptCols;
        Cursor cursor{ file.data + offset, file.data + file.size, candidate.bigEndian, candidate.wideBigInt };

        if (!parseFrameHeader(cursor, attempt, attemptCols)) continue;

        const char* frameEnd = nullptr;
        if (!measureFrameBody(attempt, file.data + file.size, frameEnd)) continue;

        convention = candidate;
        parsed = std::move(attempt);
        cols = attemptCols;
        return true;
    }

    return false;
}

bool sniff(const MappedFile& file) {
    Convention convention{};
    BinaryHeader parsed;
    ColumnMapping cols;
    return detectConvention(file, 0, convention, parsed, cols);
}

bool scan(const MappedFile& file, std::vector<FrameIndexEntry>& frames, std::string& error) {
    const char* end = file.data + file.size;
    size_t offset = 0;

    while (offset < file.size) {
        Convention convention{};
        BinaryHeader parsed;
        ColumnMapping cols;

        if (!detectConvention(file, offset, convention, parsed, cols)) break;

        const char* frameEnd = nullptr;
        if (!measureFrameBody(parsed, end, frameEnd)) break;

        FrameIndexEntry entry;
        entry.index = (int)frames.size();
        entry.byteOffset = offset;
        entry.byteLength = (size_t)(frameEnd - file.data) - offset;
        entry.timestep = parsed.header.timestep;
        entry.atomCount = parsed.header.atomCount;
        frames.push_back(entry);

        offset = (size_t)(frameEnd - file.data);
    }

    if (frames.empty()) {
        error = "Invalid LAMMPS binary dump format";
        return false;
    }

    return true;
}

bool readHeader(const MappedFile& file, const FrameIndexEntry& entry,
                FrameHeader& header, std::string& error) {
    Convention convention{};
    BinaryHeader parsed;
    ColumnMapping cols;

    if (!detectConvention(file, entry.byteOffset, convention, parsed, cols)) {
        error = "Invalid LAMMPS binary dump format";
        return false;
    }

    header = std::move(parsed.header);
    return true;
}

bool readFrame(const MappedFile& file, const FrameIndexEntry& entry, const ReadOptions& options,
               FrameAllocator& allocator, ParsedFrame& frame, std::string& error) {
    Convention convention{};
    BinaryHeader parsed;
    ColumnMapping cols;

    if (!detectConvention(file, entry.byteOffset, convention, parsed, cols)) {
        error = "Invalid LAMMPS binary dump format";
        return false;
    }

    resolveExtraColumns(options.properties, parsed.header.headers, cols);

    frame.header = parsed.header;
    frame.hasIds = options.includeIds && cols.idxId >= 0;

    const int atomCount = parsed.header.atomCount;
    FrameBuffers buffers = allocator.allocate(atomCount, frame.hasIds);

    const int extraCount = (int)cols.extraPropIndices.size();
    frame.extras.resize(extraCount);
    for (int extra = 0; extra < extraCount; extra++) {
        frame.extras[extra].name = cols.extraPropNames[extra];
        frame.extras[extra].values.resize((size_t)atomCount);
        frame.extras[extra].dtype = ColumnDtype::Int32;
    }

    const double lx = parsed.header.box.xhi - parsed.header.box.xlo;
    const double ly = parsed.header.box.yhi - parsed.header.box.ylo;
    const double lz = parsed.header.box.zhi - parsed.header.box.zlo;

    frame.bbox.init();

    Cursor cursor{ parsed.body, file.data + file.size, convention.bigEndian, convention.wideBigInt };
    int atomIndex = 0;

    for (int32_t chunk = 0; chunk < parsed.chunkCount && atomIndex < atomCount; chunk++) {
        int32_t values = 0;
        if (!cursor.readInt(values) || values < 0) {
            error = "Truncated LAMMPS binary dump chunk";
            return false;
        }

        const int rows = values / parsed.sizeOne;
        for (int row = 0; row < rows && atomIndex < atomCount; row++) {
            double fields[MAX_SIZE_ONE];
            for (int32_t field = 0; field < parsed.sizeOne; field++) {
                if (!cursor.readDouble(fields[field])) {
                    error = "Truncated LAMMPS binary dump chunk";
                    return false;
                }
            }

            double x = fields[cols.idxX];
            double y = fields[cols.idxY];
            double z = fields[cols.idxZ];

            if (parsed.scaledCoords) {
                const double fx = x, fy = y, fz = z;
                x = parsed.header.box.xlo + fx * lx + fy * parsed.header.box.xy + fz * parsed.header.box.xz;
                y = parsed.header.box.ylo + fy * ly + fz * parsed.header.box.yz;
                z = parsed.header.box.zlo + fz * lz;
            }

            buffers.setPosition(atomIndex, x, y, z);
            buffers.types[atomIndex] = cols.idxType >= 0 ? (uint16_t)fields[cols.idxType] : 0;
            if (buffers.ids) buffers.ids[atomIndex] = (uint32_t)fields[cols.idxId];

            for (int extra = 0; extra < extraCount; extra++) {
                const double value = fields[cols.extraPropIndices[extra]];
                frame.extras[extra].values[(size_t)atomIndex] = value;
                if (value != (double)(int64_t)value ||
                    value < -2147483648.0 || value > 2147483647.0) {
                    frame.extras[extra].dtype = ColumnDtype::Float32;
                }
            }

            frame.bbox.update((float)x, (float)y, (float)z);
            atomIndex++;
        }
    }

    return true;
}

extern const FormatReader reader = {
    format_id::LammpsDumpBinary,
    sniff,
    scan,
    readHeader,
    readFrame
};

}
}
