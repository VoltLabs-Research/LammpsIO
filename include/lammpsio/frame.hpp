#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <thread>
#include <utility>
#include <vector>
#include <lammpsio/common.hpp>

namespace lammpsio {

namespace format_id {
    constexpr const char* LammpsDumpText = "lammps-dump";
    constexpr const char* LammpsDumpBinary = "lammps-dump-binary";
    constexpr const char* LammpsDumpYaml = "lammps-dump-yaml";
    constexpr const char* LammpsData = "lammps-data";
    constexpr const char* ExtXyz = "extxyz";
}

struct ReadOptions {
    bool includeIds = false;
    std::vector<std::string> properties;
    int frame = 0;
    int maxThreads = 0;
};

struct ExtraColumn {
    std::string name;
    std::vector<double> values;
    ColumnDtype dtype = ColumnDtype::Int32;
};

struct FrameHeader {
    const char* format = nullptr;
    int timestep = 0;
    int atomCount = 0;
    SimulationBox box = {};
    double cell[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    double origin[3] = { 0, 0, 0 };
    bool periodic[3] = { true, true, true };
    std::vector<std::string> headers;
    bool positionsWereScaled = false;
    std::vector<std::pair<std::string, std::string>> extraSections;
};

inline void deriveCellFromBox(FrameHeader& header) {
    const SimulationBox& box = header.box;
    header.cell[0][0] = box.xhi - box.xlo;
    header.cell[0][1] = 0.0;
    header.cell[0][2] = 0.0;
    header.cell[1][0] = box.xy;
    header.cell[1][1] = box.yhi - box.ylo;
    header.cell[1][2] = 0.0;
    header.cell[2][0] = box.xz;
    header.cell[2][1] = box.yz;
    header.cell[2][2] = box.zhi - box.zlo;
    header.origin[0] = box.xlo;
    header.origin[1] = box.ylo;
    header.origin[2] = box.zlo;
}

inline void deriveBoxFromCell(FrameHeader& header) {
    header.box.xlo = header.origin[0];
    header.box.ylo = header.origin[1];
    header.box.zlo = header.origin[2];
    header.box.xhi = header.origin[0] + header.cell[0][0];
    header.box.yhi = header.origin[1] + header.cell[1][1];
    header.box.zhi = header.origin[2] + header.cell[2][2];
    header.box.xy = header.cell[1][0];
    header.box.xz = header.cell[2][0];
    header.box.yz = header.cell[2][1];
}

struct FrameIndexEntry {
    int index = 0;
    size_t byteOffset = 0;
    size_t byteLength = 0;
    int timestep = 0;
    int atomCount = 0;
};

struct FrameBuffers {
    float* positions32 = nullptr;
    double* positions64 = nullptr;
    uint16_t* types = nullptr;
    uint32_t* ids = nullptr;

    HOT ALWAYS_INLINE void setPosition(int atomIndex, double x, double y, double z) {
        const int base = atomIndex * 3;
        if (positions32) {
            positions32[base] = (float)x;
            positions32[base + 1] = (float)y;
            positions32[base + 2] = (float)z;
        } else {
            positions64[base] = x;
            positions64[base + 1] = y;
            positions64[base + 2] = z;
        }
    }
};

enum class PositionPrecision { Float32, Float64 };

struct FrameAllocator {
    virtual FrameBuffers allocate(int atomCount, bool withIds) = 0;
    virtual PositionPrecision positionPrecision() const { return PositionPrecision::Float32; }
    virtual ~FrameAllocator() = default;
};

inline unsigned int resolveThreadCount(int maxThreads, int atomCount, int multithreadThreshold) {
    if (atomCount < multithreadThreshold) return 1;

    unsigned int available = maxThreads > 0
        ? (unsigned int)maxThreads
        : std::thread::hardware_concurrency();

    return available == 0 ? 1 : available;
}

struct ParsedFrame {
    FrameHeader header;
    BoundingBox bbox = {};
    std::vector<ExtraColumn> extras;
    bool hasIds = false;
    std::vector<double> massesByType;
    std::vector<std::string> elementHintsByType;
};

struct FormatReader {
    const char* id;
    bool (*sniff)(const MappedFile& file);
    bool (*scan)(const MappedFile& file, std::vector<FrameIndexEntry>& frames, std::string& error);
    bool (*readHeader)(const MappedFile& file, const FrameIndexEntry& entry,
                       FrameHeader& header, std::string& error);
    bool (*readFrame)(const MappedFile& file, const FrameIndexEntry& entry,
                      const ReadOptions& options, FrameAllocator& allocator,
                      ParsedFrame& frame, std::string& error);
};

inline void resolveExtraColumns(const std::vector<std::string>& requested,
                                const std::vector<std::string>& headers,
                                ColumnMapping& cols) {
    if (requested.empty()) return;

    const auto isConsumed = [&cols](int index) {
        return index == cols.idxId || index == cols.idxType ||
               index == cols.idxX || index == cols.idxY || index == cols.idxZ;
    };

    bool wildcard = false;
    for (const auto& name : requested) {
        if (name == "*") { wildcard = true; break; }
    }

    if (wildcard) {
        for (size_t index = 0; index < headers.size(); index++) {
            if (isConsumed((int)index)) continue;
            cols.extraPropIndices.push_back((int)index);
            cols.extraPropNames.push_back(headers[index]);
        }
    } else {
        for (const auto& name : requested) {
            for (size_t index = 0; index < headers.size(); index++) {
                if (headers[index] != name || isConsumed((int)index)) continue;
                cols.extraPropIndices.push_back((int)index);
                cols.extraPropNames.push_back(name);
                break;
            }
        }
    }

    if (!cols.extraPropIndices.empty()) cols.computeMaxIdx();
}

}
