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

    const { positions } = readFrame(file, { includeIds: true });
    assert.deepEqual(round(positions), [0, 0, 0, 5, 10, 15, 10, 20, 30]);
});

test('scaled xs/ys/zs honor triclinic tilt and recovered box edges', () => {
    const file = writeDump('scaled-tri.dump', [{
        timestep: 0,
        bounds: ['0.0 12.0 2.0', '0.0 20.0 0.0', '0.0 30.0 0.0'],
        boundaryFlags: 'xy xz yz pp pp pp',
        columns: 'id type xs ys zs',
        rows: ['1 1 0.0 0.0 0.0', '2 1 0.5 0.5 0.5', '3 1 1.0 1.0 1.0']
    }]);

    const { positions } = readFrame(file, { includeIds: true });
    assert.deepEqual(round(positions), [0, 0, 0, 6, 10, 15, 12, 20, 30]);
});

test('a scaled dump does not report its xs/ys/zs columns as extra properties', () => {
    const file = writeDump('scaled-wildcard.dump', [{
        timestep: 0,
        bounds: ORTHO_BOUNDS,
        columns: 'id type xs ys zs c_pe',
        rows: ['1 1 0.0 0.0 0.0 -1.0', '2 1 0.5 0.5 0.5 -2.0']
    }]);

    const frame = readFrame(file, { includeIds: true, properties: ['*'] });
    assert.deepEqual(Object.keys(frame.properties), ['c_pe']);
    assert.deepEqual(round(frame.positions), [0, 0, 0, 5, 10, 15]);
    assert.deepEqual(frame.metadata.headers, ['id', 'type', 'xs', 'ys', 'zs', 'c_pe']);
});
