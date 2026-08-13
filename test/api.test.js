// The four entry points: format detection, frame indexing, header-only reads, and full
// frame reads. The multi-frame cases are the ones that matter most — the previous
// version of this addon read frame 0 of a multi-frame file and silently dropped the
// rest, which is invisible to a caller that never counts.

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');

const { detectFormat, scanFrames, readHeader, readFrame } = require('..');
const { ORTHO_BOUNDS, write, writeDump } = require('./fixtures');

const threeFrameDump = () => writeDump('three-frames.dump', [0, 100, 200].map((timestep, index) => ({
    timestep,
    bounds: ORTHO_BOUNDS,
    columns: 'id type x y z c_pe',
    rows: [
        `1 1 ${index}.0 0.0 0.0 -3.${index}`,
        `2 2 ${index}.5 1.0 1.0 -2.${index}`
    ]
})));

const DATA_FILE = [
    'LAMMPS data file',
    '',
    '4 atoms',
    '2 atom types',
    '',
    '0.0 10.0 xlo xhi',
    '0.0 20.0 ylo yhi',
    '0.0 30.0 zlo zhi',
    '',
    'Masses',
    '',
    '1 55.845   # Fe',
    '2 12.011   # C',
    '',
    'Atoms # atomic',
    '',
    '1 1 0.0 0.0 0.0',
    '2 1 2.0 0.0 0.0',
    '3 2 0.0 2.0 0.0',
    '4 2 0.0 0.0 2.0',
    ''
].join('\n');

test('detectFormat names the format, or returns null for anything unrecognized', () => {
    assert.equal(detectFormat(threeFrameDump()), 'lammps-dump');
    assert.equal(detectFormat(write('atoms.data', DATA_FILE)), 'lammps-data');
    assert.equal(detectFormat(write('notes.txt', 'hello world\nnot a trajectory\n')), null);
});

test('detectFormat throws on a missing file rather than calling it unrecognized', () => {
    assert.throws(() => detectFormat('/nonexistent/path.dump'), /Failed to open file/);
});

test('scanFrames indexes every frame of a multi-frame dump', () => {
    const { format, frames } = scanFrames(threeFrameDump());

    assert.equal(format, 'lammps-dump');
    assert.equal(frames.length, 3);
    assert.deepEqual(frames.map((frame) => frame.timestep), [0, 100, 200]);
    assert.deepEqual(frames.map((frame) => frame.index), [0, 1, 2]);
    assert.deepEqual(frames.map((frame) => frame.natoms), [2, 2, 2]);
});

test('frame byte ranges tile the file exactly and each one is a valid dump on its own', () => {
    // This is the contract a caller relies on to split a multi-frame upload into
    // per-frame files by copying bytes, with no reparse and no reserialization.
    const file = threeFrameDump();
    const bytes = fs.readFileSync(file);
    const { frames } = scanFrames(file);

    assert.equal(frames[0].byteOffset, 0);
    assert.equal(
        frames.reduce((total, frame) => total + frame.byteLength, 0),
        bytes.length
    );

    frames.forEach((frame, index) => {
        if (index > 0) {
            const previous = frames[index - 1];
            assert.equal(frame.byteOffset, previous.byteOffset + previous.byteLength);
        }

        const slice = write(`slice-${frame.timestep}.dump`,
            bytes.subarray(frame.byteOffset, frame.byteOffset + frame.byteLength));

        assert.equal(scanFrames(slice).frames.length, 1);
        assert.equal(readFrame(slice).metadata.timestep, frame.timestep);
    });
});

test('readFrame reads the requested frame, not always the first', () => {
    const file = threeFrameDump();

    assert.equal(readFrame(file).metadata.timestep, 0);
    assert.equal(readFrame(file, { frame: 1 }).metadata.timestep, 100);
    assert.equal(readFrame(file, { frame: 2 }).metadata.timestep, 200);
    assert.deepEqual(Array.from(readFrame(file, { frame: 2 }).positions.slice(0, 3)), [2, 0, 0]);
});

test('an out-of-range frame is a range error naming what the file actually has', () => {
    const file = threeFrameDump();
    assert.throws(() => readFrame(file, { frame: 3 }), /Frame 3 is out of range.*3 frame/s);
    assert.throws(() => readFrame(file, { frame: -1 }), RangeError);
});

test('readHeader describes a frame without reading its atoms', () => {
    const header = readHeader(threeFrameDump(), { frame: 1 });

    assert.equal(header.format, 'lammps-dump');
    assert.equal(header.timestep, 100);
    assert.equal(header.natoms, 2);
    assert.deepEqual(header.headers, ['id', 'type', 'x', 'y', 'z', 'c_pe']);
    assert.deepEqual(header.pbc, [true, true, true]);
    assert.equal(header.boxBounds.xlo, 0);
    assert.equal(header.boxBounds.xhi, 10);
    assert.equal(header.boxBounds.xy, 0);
    // Header reads carry no atom data at all.
    assert.equal(header.positions, undefined);
});

test('boundary flags are reported per axis', () => {
    const file = writeDump('mixed-boundaries.dump', [{
        timestep: 0,
        bounds: ORTHO_BOUNDS,
        boundaryFlags: 'pp ff pp',
        columns: 'id type x y z',
        rows: ['1 1 1.0 2.0 3.0']
    }]);

    assert.deepEqual(readHeader(file).pbc, [true, false, true]);
});

test('property dtypes distinguish a categorical integer column from a continuous one', () => {
    const file = writeDump('dtypes.dump', [{
        timestep: 0,
        bounds: ORTHO_BOUNDS,
        columns: 'id type x y z mol c_pe',
        rows: ['1 1 0.0 0.0 0.0 3 -1.5', '2 1 1.0 1.0 1.0 4 -2.5']
    }]);

    const frame = readFrame(file, { includeIds: true, properties: ['*'] });

    assert.deepEqual(frame.propertyDtypes, { mol: 'i32', c_pe: 'f32' });
    assert.ok(frame.properties.mol instanceof Int32Array);
    assert.ok(frame.properties.c_pe instanceof Float32Array);
});

test('named properties are honored and unknown names are simply absent', () => {
    const file = writeDump('named.dump', [{
        timestep: 0,
        bounds: ORTHO_BOUNDS,
        columns: 'id type x y z mol c_pe',
        rows: ['1 1 0.0 0.0 0.0 3 -1.5']
    }]);

    const frame = readFrame(file, { properties: ['c_pe', 'nonexistent'] });
    assert.deepEqual(Object.keys(frame.properties), ['c_pe']);
});

test('ids are only materialized when asked for', () => {
    const file = threeFrameDump();

    assert.equal(readFrame(file).ids, undefined);
    assert.deepEqual(Array.from(readFrame(file, { includeIds: true }).ids), [1, 2]);
});

test('a data file is a single frame and carries its Masses section', () => {
    const file = write('masses.data', DATA_FILE);
    const { format, frames } = scanFrames(file);

    assert.equal(format, 'lammps-data');
    assert.equal(frames.length, 1);
    assert.equal(frames[0].natoms, 4);

    const frame = readFrame(file, { includeIds: true });
    assert.equal(frame.metadata.natoms, 4);
    assert.deepEqual(Array.from(frame.types), [1, 1, 2, 2]);
    assert.deepEqual(Array.from(frame.ids), [1, 2, 3, 4]);
    // 1-indexed by type, so index 0 is type 1.
    assert.deepEqual(frame.massesByType, [55.845, 12.011]);
    assert.deepEqual(frame.elementHintsByType, ['Fe', 'C']);
});

test('the declared atom_style picks the column layout over the column count', () => {
    // `full` is id mol type q x y z — seven columns, the same count as several other
    // styles. Guessing from the count alone would read the charge as a coordinate.
    const file = write('full.data', [
        'LAMMPS data file', '',
        '1 atoms', '1 atom types', '',
        '0.0 10.0 xlo xhi', '0.0 10.0 ylo yhi', '0.0 10.0 zlo zhi', '',
        'Atoms # full', '',
        '1 7 1 -0.5 1.0 2.0 3.0', ''
    ].join('\n'));

    const frame = readFrame(file, { includeIds: true });
    assert.deepEqual(Array.from(frame.positions), [1, 2, 3]);
    assert.deepEqual(Array.from(frame.types), [1]);
});

test('a triclinic data file reports its tilt factors', () => {
    const file = write('tilted.data', [
        'LAMMPS data file', '',
        '1 atoms', '1 atom types', '',
        '0.0 10.0 xlo xhi', '0.0 10.0 ylo yhi', '0.0 10.0 zlo zhi',
        '2.0 0.5 0.25 xy xz yz', '',
        'Atoms # atomic', '',
        '1 1 1.0 2.0 3.0', ''
    ].join('\n'));

    const { boxBounds } = readHeader(file);
    assert.equal(boxBounds.xy, 2);
    assert.equal(boxBounds.xz, 0.5);
    assert.equal(boxBounds.yz, 0.25);
});

test('an unrecognized file is an error naming the path, not a silent empty read', () => {
    const file = write('unsupported.txt', 'just some text\n');
    assert.throws(() => readFrame(file), /Unsupported trajectory format.*unsupported\.txt/);
    assert.throws(() => scanFrames(file), /Unsupported trajectory format/);
});

test('each frame of a multi-frame XYZ gets a distinct timestep', () => {
    // An XYZ frame carries no timestep, so its index stands in for one. Callers key
    // frames by timestep; identical values would collapse the trajectory to one frame.
    const file = write('two-frames.xyz', [
        '1', 'Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3', 'Fe 1.0 1.0 1.0',
        '1', 'Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3', 'Fe 2.0 2.0 2.0',
        ''
    ].join('\n'));

    const { frames } = scanFrames(file);
    assert.equal(frames.length, 2);
    assert.deepEqual(frames.map((frame) => frame.timestep), [0, 1]);
    assert.deepEqual([0, 1].map((index) => readHeader(file, { frame: index }).timestep), [0, 1]);
    assert.deepEqual([0, 1].map((index) => readFrame(file, { frame: index }).metadata.timestep), [0, 1]);
});

test('XYZ species become numeric types with their symbols returned as hints', () => {
    const file = write('species.xyz', [
        '4', 'Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3',
        'Fe 0.0 0.0 0.0', 'C 1.0 0.0 0.0', 'Fe 2.0 0.0 0.0', 'O 3.0 0.0 0.0', ''
    ].join('\n'));

    const frame = readFrame(file);
    // Numbered in order of first appearance, 1-based.
    assert.deepEqual(Array.from(frame.types), [1, 2, 1, 3]);
    assert.deepEqual(frame.elementHintsByType, ['Fe', 'C', 'O']);
});

test('a dump section the reader has no meaning for is carried through, not dropped', () => {
    // `ITEM: TIME` is the common one. A consumer re-emitting the frame needs it back, so
    // unknown sections travel with the header instead of being skipped.
    const file = write('with-time.dump', [
        'ITEM: TIMESTEP', '400',
        'ITEM: TIME', '0.0004',
        'ITEM: NUMBER OF ATOMS', '1',
        'ITEM: BOX BOUNDS pp pp pp', '0.0 10.0', '0.0 10.0', '0.0 10.0',
        'ITEM: ATOMS id type x y z', '1 1 1.0 2.0 3.0', ''
    ].join('\n'));

    const header = readHeader(file);
    assert.equal(header.timestep, 400);
    assert.deepEqual(header.extraSections, [['TIME', '0.0004']]);
});

test('positionsWereScaled records what the file held, not what was returned', () => {
    const scaled = writeDump('was-scaled.dump', [{
        timestep: 0, bounds: ORTHO_BOUNDS, columns: 'id type xs ys zs',
        rows: ['1 1 0.5 0.5 0.5']
    }]);
    const cartesian = writeDump('was-cartesian.dump', [{
        timestep: 0, bounds: ORTHO_BOUNDS, columns: 'id type x y z',
        rows: ['1 1 5.0 10.0 15.0']
    }]);

    // Both return Cartesian positions; only one came that way.
    assert.equal(readHeader(scaled).positionsWereScaled, true);
    assert.equal(readHeader(cartesian).positionsWereScaled, false);
    assert.deepEqual(Array.from(readFrame(scaled).positions), [5, 10, 15]);
    assert.deepEqual(Array.from(readFrame(cartesian).positions), [5, 10, 15]);
});
