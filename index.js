// Native readers for LAMMPS and LAMMPS-adjacent trajectory formats.
//
// One entry point per question, and every one of them identifies the format itself —
// callers never have to guess, and never have to try one reader and fall back to
// another. Frame indices come from scanFrames().
//
//   detectFormat(path)
//     → 'lammps-dump' | 'lammps-data' | null
//     null means "nothing here recognizes it". A missing or unreadable file throws.
//
//   scanFrames(path)
//     → { format, frames: [{ index, byteOffset, byteLength, timestep, natoms }] }
//     Frame boundaries as byte ranges, so a multi-frame file can be split into
//     single-frame files by copying bytes, with no reparse. Single-frame formats
//     report one entry spanning the file.
//
//   readHeader(path, { frame = 0 })
//     → { format, timestep, natoms, boxBounds, pbc, headers }
//     Header only, no atom data: what a caller needs to accept an upload and describe
//     its cell without paying for the positions.
//
//   readFrame(path, { frame = 0, includeIds = false, properties = [] })
//     → { positions: Float32Array, types: Uint16Array, ids?: Uint32Array,
//         properties?: Record<string, Int32Array | Float32Array>,
//         propertyDtypes?: Record<string, 'i32' | 'f32'>,
//         metadata, min, max,
//         massesByType?: number[],              // data files, 1-indexed by type
//         elementHintsByType?: (string|null)[] } // data files, trailing "# <symbol>"
//     `properties: ['*']` requests every per-atom column the file carries beyond
//     id/type/position.
//
// boxBounds holds the true (untilted) cell edges plus the triclinic tilt factors
// xy/xz/yz, already recovered from the bounding box a dump prints.
module.exports = require('./build/Release/lammps_io.node');
