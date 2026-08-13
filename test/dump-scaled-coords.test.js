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

const { readFrame } = require('..');
const { ORTHO_BOUNDS, writeDump } = require('./fixtures');

const round = (values, digits = 2) =>
    Array.from(values, (value) => Number(value.toFixed(digits)));

test('unscaled x/y/z coordinates pass through unchanged', () => {
    const file = writeDump('ortho.dump', [{
        timestep: 0,
        bounds: ORTHO_BOUNDS,
        columns: 'id type x y z',
        rows: ['1 1 1.0 2.0 3.0', '2 1 5.0 10.0 15.0', '3 2 9.0 18.0 27.0']
    }]);

    const { positions } = readFrame(file, { includeIds: true });
    assert.deepEqual(round(positions), [1, 2, 3, 5, 10, 15, 9, 18, 27]);
});

test('scaled xs/ys/zs map to Cartesian via the orthogonal box', () => {
    const file = writeDump('scaled-ortho.dump', [{
        timestep: 0,
        bounds: ORTHO_BOUNDS,
        columns: 'id type xs ys zs',
        rows: ['1 1 0.0 0.0 0.0', '2 1 0.5 0.5 0.5', '3 1 1.0 1.0 1.0']
    }]);

    // x = xs*lx, y = ys*ly, z = zs*lz  (lx,ly,lz = 10,20,30)
    const { positions } = readFrame(file, { includeIds: true });
    assert.deepEqual(round(positions), [0, 0, 0, 5, 10, 15, 10, 20, 30]);
});

test('scaled xs/ys/zs honor triclinic tilt and recovered box edges', () => {
    // Printed bbox x:[0,12] with xy tilt 2 → true edge xlo=0, xhi=10 (lx=10).
    const file = writeDump('scaled-tri.dump', [{
        timestep: 0,
        bounds: ['0.0 12.0 2.0', '0.0 20.0 0.0', '0.0 30.0 0.0'],
        boundaryFlags: 'xy xz yz pp pp pp',
        columns: 'id type xs ys zs',
        rows: ['1 1 0.0 0.0 0.0', '2 1 0.5 0.5 0.5', '3 1 1.0 1.0 1.0']
    }]);

    // x = xlo + xs*lx + ys*xy + zs*xz ; y = ylo + ys*ly + zs*yz ; z = zlo + zs*lz
    const { positions } = readFrame(file, { includeIds: true });
    assert.deepEqual(round(positions), [0, 0, 0, 6, 10, 15, 12, 20, 30]);
});

test('a scaled dump does not report its xs/ys/zs columns as extra properties', () => {
    // The wildcard used to exclude base columns by name, so a scaled dump reported
    // xs/ys/zs as data columns on top of the positions they had become.
    const file = writeDump('scaled-wildcard.dump', [{
        timestep: 0,
        bounds: ORTHO_BOUNDS,
        columns: 'id type xs ys zs c_pe',
        rows: ['1 1 0.0 0.0 0.0 -1.0', '2 1 0.5 0.5 0.5 -2.0']
    }]);

    const frame = readFrame(file, { includeIds: true, properties: ['*'] });
    assert.deepEqual(Object.keys(frame.properties), ['c_pe']);
    assert.deepEqual(round(frame.positions), [0, 0, 0, 5, 10, 15]);
    // The raw column names stay visible for callers that want them.
    assert.deepEqual(frame.metadata.headers, ['id', 'type', 'xs', 'ys', 'zs', 'c_pe']);
});
