// Native LAMMPS parsers (N-API addons).
//
// dumpParser.parseDump(path, { includeIds, properties }) → {
//   positions: Float32Array, types: Uint16Array, ids?: Uint32Array,
//   properties?: Record<string, Int32Array | Float32Array>,   // typed per column
//   propertyDtypes?: Record<string, 'i32' | 'f32'>,           // parallel dtype map
//   metadata: { timestep, natoms, boxBounds, headers }, min, max
// }
//
// dataParser.parseData(path, { includeIds }) → {
//   positions, types, ids?, metadata, min, max,
//   massesByType?: number[],         // 1-indexed by LAMMPS type (index 0 = type 1)
//   elementHintsByType?: (string|null)[]  // trailing "# <symbol>" comment per type
// }
module.exports = {
  dataParser: require('./build/Release/data_parser.node'),
  dumpParser: require('./build/Release/dump_parser.node')
};
