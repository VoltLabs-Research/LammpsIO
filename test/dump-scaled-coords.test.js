// Regression coverage for LAMMPS scaled-coordinate dumps (xs/ys/zs).
//
// Before the fix, the dump parser matched any 1–2 char column starting with
// x/y/z, so it read fractional `xs ys zs` as if they were Cartesian `x y z`
// and emitted positions in [0,1]. Downstream that produced a unit-cube GLB
// that rendered tiny inside a full-size simulation cell. These tests pin the
// fractional→Cartesian mapping (including triclinic tilt) and guard the
// unscaled path against regressing.

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');

const { dumpParser } = require('..');

const writeDump = (name, body) => {
    const file = path.join(os.tmpdir(), `lammpsio-${name}-${process.pid}.dump`);
    fs.writeFileSync(file, body);
    return file;
};

const round = (values, digits = 2) =>
    Array.from(values, (value) => Number(value.toFixed(digits)));

test('unscaled x/y/z coordinates pass through unchanged', () => {
    const file = writeDump('ortho', [
        'ITEM: TIMESTEP', '0',
        'ITEM: NUMBER OF ATOMS', '3',
        'ITEM: BOX BOUNDS pp pp pp',
        '0.0 10.0', '0.0 20.0', '0.0 30.0',
        'ITEM: ATOMS id type x y z',
        '1 1 1.0 2.0 3.0',
        '2 1 5.0 10.0 15.0',
        '3 2 9.0 18.0 27.0', ''
    ].join('\n'));

    const { positions } = dumpParser.parseDump(file, { includeIds: true });
    assert.deepEqual(round(positions), [1, 2, 3, 5, 10, 15, 9, 18, 27]);
});

test('scaled xs/ys/zs map to Cartesian via the orthogonal box', () => {
    const file = writeDump('scaled-ortho', [
        'ITEM: TIMESTEP', '0',
        'ITEM: NUMBER OF ATOMS', '3',
        'ITEM: BOX BOUNDS pp pp pp',
        '0.0 10.0', '0.0 20.0', '0.0 30.0',
        'ITEM: ATOMS id type xs ys zs',
        '1 1 0.0 0.0 0.0',
        '2 1 0.5 0.5 0.5',
        '3 1 1.0 1.0 1.0', ''
    ].join('\n'));

    const { positions } = dumpParser.parseDump(file, { includeIds: true });
    // x = xs*lx, y = ys*ly, z = zs*lz  (lx,ly,lz = 10,20,30)
    assert.deepEqual(round(positions), [0, 0, 0, 5, 10, 15, 10, 20, 30]);
});

test('scaled xs/ys/zs honor triclinic tilt and recovered box edges', () => {
    // Printed bbox x:[0,12] with xy tilt 2 → true edge xlo=0, xhi=10 (lx=10).
    const file = writeDump('scaled-tri', [
        'ITEM: TIMESTEP', '0',
        'ITEM: NUMBER OF ATOMS', '3',
        'ITEM: BOX BOUNDS xy xz yz pp pp pp',
        '0.0 12.0 2.0', '0.0 20.0 0.0', '0.0 30.0 0.0',
        'ITEM: ATOMS id type xs ys zs',
        '1 1 0.0 0.0 0.0',
        '2 1 0.5 0.5 0.5',
        '3 1 1.0 1.0 1.0', ''
    ].join('\n'));

    const { positions } = dumpParser.parseDump(file, { includeIds: true });
    // x = xlo + xs*lx + ys*xy + zs*xz ; y = ylo + ys*ly + zs*yz ; z = zlo + zs*lz
    assert.deepEqual(round(positions), [0, 0, 0, 6, 10, 15, 12, 20, 30]);
});
