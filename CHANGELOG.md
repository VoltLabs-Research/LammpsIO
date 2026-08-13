# Changelog

## 2.1.0

The reader core becomes consumable from C++ as well, so the analysis plugins read
trajectories through the same code the daemon does instead of their own parser. Additive
for the npm package: nothing in the JavaScript API changed shape.

### Added

- A CMake target and a Conan recipe (`lammpsio/2.1.0`) that build the readers as a plain
  C++ static library. `include/lammpsio/reader_registry.hpp` is the public surface:
  `detectFormat`, `scanFrames`, `readHeader`, `readFrame`, plus a `VectorFrameAllocator`
  for consumers with no special memory requirements.
- Positions are written in the precision the consumer asks for. `FrameAllocator` declares
  float32 or float64 and readers keep values in double until the write, so a renderer and a
  float64 geometry pipeline are both served without a conversion pass. The Node addon still
  gets float32 straight into its V8 buffers.
- `maxThreads` in the read options, for callers already inside a parallel region. A plugin
  running under a TBB task would otherwise start a nested pool on every core it was given.
- `positionsWereScaled` on the header: whether the file held fractional coordinates.
  Positions come back Cartesian either way, but a consumer re-emitting the file needs to
  know how it arrived.
- `extraSections` on the header: `ITEM:` sections the reader has no meaning for, in order,
  so `ITEM: TIME` and friends survive a round trip instead of being dropped.

### Changed

- Everything lives in a `lammpsio` namespace. Names like `BoundingBox` and `FrameHeader`
  at global scope would collide in any sizeable C++ project that links this.
- Public headers moved to `include/lammpsio/`, implementations stay in `src/`, so
  `#include <lammpsio/...>` resolves identically in the source tree and against the
  installed package.
- Format dispatch moved out of the N-API layer into `reader_registry.cpp`. Both builds
  compile the same readers; `registry.cpp` is now only argument reading and result
  conversion.

### Fixed

- `findLineEnd` and `jumpToNextLine` clamp when the cursor has already passed the end of
  the mapping, which a truncated frame can cause. Previously that handed `memchr` a huge
  unsigned length.

## 2.0.0

One addon with one API, replacing the two separate parsers and their per-parser entry
points. Breaking: `dumpParser.parseDump` and `dataParser.parseData` are gone.

### Added

- `detectFormat(path)` — names the format, or returns `null` when nothing recognizes the
  file. Callers no longer try one parser and fall back to another.
- `scanFrames(path)` — the byte range, timestep and atom count of every frame. A
  multi-frame file can be split into single-frame files by copying bytes, no reparse.
- `readHeader(path, { frame })` — cell, atom count and column names without reading any
  atom data.
- `readFrame(path, { frame, includeIds, properties })` — replaces both old parsers, and
  takes a frame index.
- Extended XYZ and plain XYZ (`extxyz`): `Lattice=`, `Properties=`, `pbc=`, element
  symbols mapped to numeric types with the symbols returned as element hints.
- LAMMPS binary dumps (`lammps-dump-binary`), with the endianness and integer width of
  the writing build recovered by probing, both the pre- and post-2018 header layouts, and
  the chunked body.
- LAMMPS YAML dumps (`lammps-dump-yaml`). Reads the restricted subset LAMMPS emits rather
  than depending on a YAML parser, and fails instead of guessing on anything else.
- `cellVectors` and `cellOrigin` in every header: the three cell vectors as read. An
  extended-XYZ lattice need not be upper triangular, which `boxBounds` cannot express.
- Per-axis periodic boundary flags (`pbc`), read from a dump's `BOX BOUNDS` line.
- LAMMPS data files honor the `atom_style` declared on the `Atoms` line instead of
  guessing the layout from the column count. Several styles share a count, so guessing
  read `full`'s charge column as a coordinate.
- LAMMPS data files parse the `xy xz yz` tilt line.
- A differential test suite that checks every fixture against OVITO 3.15.5, with the
  ground truth committed so it runs without OVITO installed.

### Fixed

- **Multi-frame files no longer read as their first frame.** `parseDump` stopped after
  frame 0 and silently dropped the rest, which a caller could only notice by counting.
- **The `'*'` property wildcard no longer duplicates position columns.** Base columns
  were excluded by name (`id`, `type`, `x`, `y`, `z`), so a scaled dump reported its
  `xs`/`ys`/`zs` columns as data on top of the Cartesian positions they had become.
  Exclusion is now by resolved column index.
- **A triclinic data file no longer reads as orthogonal.** The header scan stopped as
  soon as the required values were in, and the optional tilt line comes after them.
- Element hints are returned for any format that carries them, not only alongside
  masses — an extended-XYZ file has species symbols and no masses.
- An unreadable file and an unrecognized one are now distinguishable: the first throws,
  the second returns `null` from `detectFormat`.

### Changed

- Built as a single `lammps_io.node` target. The N-API result building that was
  duplicated across the two parsers now lives in one place (`src/napi_bridge.hpp`), and
  a format is a self-contained translation unit registered in one line.
- Dropped `-ffast-math`. The hot path is byte scanning and `fast_float`, neither of which
  it helps, and it makes floating-point behavior around non-finite values unpredictable.
- Added a `LICENSE` (MIT), with the attribution the ported OVITO readers require.

### Retained

- The memory-mapped, `fast_float`-based scanning and the multi-threaded chunked parse for
  text formats (single-threaded below 50k atoms, where the threads cost more than they
  save).
- Zero-copy delivery of atom data: readers write positions, types and ids straight into
  the V8-visible buffers they are handed.

## 1.0.2

- Per-column dtypes (`propertyDtypes`), `massesByType` and `elementHintsByType`.
- `properties: ['*']` wildcard.
- Removed the unused `statsParser` module.

## 1.0.1

- Scaled-coordinate dumps (`xs`/`ys`/`zs`) map to Cartesian through the box and its tilt
  factors, instead of being read as if already Cartesian.
