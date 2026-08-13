#pragma once

// The format-agnostic model every reader produces, plus the interface a reader
// implements. Nothing here knows about N-API: the bridge in napi_bridge.hpp is the
// only place that does, so a new format is a self-contained translation unit that
// never repeats the result-building boilerplate.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>
#include "common.hpp"

/** Wire identifiers, also what detectFormat() returns. */
namespace format_id {
    constexpr const char* LammpsDumpText = "lammps-dump";
    constexpr const char* LammpsDumpBinary = "lammps-dump-binary";
    constexpr const char* LammpsDumpYaml = "lammps-dump-yaml";
    constexpr const char* LammpsData = "lammps-data";
    constexpr const char* ExtXyz = "extxyz";
}

struct ReadOptions {
    bool includeIds = false;
    /** Extra per-atom columns to extract. A single "*" means every non-base column. */
    std::vector<std::string> properties;
    /** Which frame of a multi-frame file to read. */
    int frame = 0;
};

/**
 * One extra per-atom column, staged as double while it is being read.
 *
 * Double is exact for integers up to 2^53, so both the dtype decision and the final
 * narrowing cast are lossless — a categorical column written as "2" stays i32, and one
 * written as "2.0" becomes f32.
 */
struct ExtraColumn {
    std::string name;
    std::vector<double> values;
    ColumnDtype dtype = ColumnDtype::Int32;
};

struct FrameHeader {
    int timestep = 0;
    int atomCount = 0;
    /**
     * LAMMPS-shaped cell: lo/hi per axis plus tilt factors. Cannot express a general
     * lattice, which is why `cell` below is the primary representation.
     */
    SimulationBox box = {};
    /**
     * The three cell vectors, `cell[0]` = a, and the cell origin. A general lattice —
     * an extended-XYZ file can carry one that is not upper triangular, so lo/hi/tilt
     * alone would silently reshape it.
     */
    double cell[3][3] = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };
    double origin[3] = { 0, 0, 0 };
    bool periodic[3] = { true, true, true };
    /** Column names as they appeared in the file, lowercased. */
    std::vector<std::string> headers;
};

/**
 * Fills the cell vectors from a LAMMPS box description, which is upper triangular by
 * construction: a along x, b tilted by xy, c tilted by xz and yz.
 */
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

/**
 * Fills the LAMMPS box description from cell vectors, for formats that carry a lattice
 * directly. Exact when the lattice is upper triangular (what LAMMPS and ASE write); for
 * anything else the tilts still describe the shear, but lo/hi describe the bounding
 * extent rather than an exact LAMMPS cell — `cell` remains the faithful representation.
 */
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

/**
 * Where a frame lives inside its file. `byteOffset`/`byteLength` are what let a caller
 * split a multi-frame file into single-frame files by copying bytes, with no reparse.
 */
struct FrameIndexEntry {
    int index = 0;
    size_t byteOffset = 0;
    size_t byteLength = 0;
    int timestep = 0;
    int atomCount = 0;
};

/** Raw destination buffers for the bulk per-atom data. */
struct FrameBuffers {
    float* positions = nullptr;
    uint16_t* types = nullptr;
    /** Null unless ids were requested and the format carries them. */
    uint32_t* ids = nullptr;
};

/**
 * Hands out the buffers a reader fills.
 *
 * The N-API implementation allocates them straight inside V8-visible ArrayBuffers, so
 * atom data is written exactly once and never copied on its way to JavaScript. A
 * reader must call this exactly once, after it knows the atom count.
 */
struct FrameAllocator {
    virtual FrameBuffers allocate(int atomCount, bool withIds) = 0;
    virtual ~FrameAllocator() = default;
};

struct ParsedFrame {
    FrameHeader header;
    BoundingBox bbox = {};
    std::vector<ExtraColumn> extras;
    bool hasIds = false;
    /** Data files only: 1-indexed by LAMMPS type, so index 0 holds type 1. */
    std::vector<double> massesByType;
    /** Data files only: the trailing `# <symbol>` comment on each Masses row. */
    std::vector<std::string> elementHintsByType;
};

/**
 * A format implementation. Plain function pointers rather than a class hierarchy —
 * there is one static instance per format and no state to carry.
 *
 * `sniff` gets the whole mapped file and decides cheaply whether this is its format.
 * `scan` lists the frames; single-frame formats return one entry spanning the file.
 * The two read functions take an entry from `scan`, so they never search for a frame.
 * All three report failure by returning false and filling `error`.
 */
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

/**
 * Resolves requested column names to indices into `headers`.
 *
 * "*" expands to every column that was not already consumed as id/type/position.
 * Exclusion is by *index*, not by name, because a scaled dump's position columns are
 * called xs/ys/zs — excluding a fixed name list reported those as extra columns
 * alongside the Cartesian positions they had already been turned into.
 */
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
