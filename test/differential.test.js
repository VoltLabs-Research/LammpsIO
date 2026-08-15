const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');

const { scanFrames, readFrame, readHeader } = require('..');

const FIXTURES = path.join(__dirname, 'fixtures');
const oracle = JSON.parse(fs.readFileSync(path.join(__dirname, 'oracle.json'), 'utf8'));

const TOLERANCE = 1e-4;

const OVITO_PROPERTY_ALIASES = {
    'Molecule Identifier': 'mol',
    Charge: 'q',
    'Potential Energy': 'local_energy',
    'Kinetic Energy': 'kinetic_energy',
    'Total Energy': 'total_energy',
    Coordination: 'n_neighb',
    Cluster: 'cluster',
    Radius: 'radius'
};

const NOT_COMPARED = {
    Mass: 'returned as massesByType instead of per particle',
    'Velocity Magnitude': 'a magnitude OVITO derives, not a column the file carries',
    Charge: 'data files carry no column names to key extra columns by',
    'Molecule Identifier': 'data files carry no column names to key extra columns by'
};

const closeTo = (actual, expected, what) => {
    assert.ok(
        Math.abs(actual - expected) <= TOLERANCE,
        `${what}: ${actual} is not within ${TOLERANCE} of ${expected}`
    );
};

const cellFromBox = (box) => ({
    a: [box.xhi - box.xlo, 0, 0],
    b: [box.xy, box.yhi - box.ylo, 0],
    c: [box.xz, box.yz, box.zhi - box.zlo],
    origin: [box.xlo, box.ylo, box.zlo]
});

const isDataFile = (name) => name.endsWith('.data');

for (const [name, expected] of Object.entries(oracle.fixtures)) {
    const file = path.join(FIXTURES, name);

    test(`${name}: frame count matches OVITO`, () => {
        assert.equal(scanFrames(file).frames.length, expected.frameCount);
    });

    expected.frames.forEach((reference, index) => {
        test(`${name}: frame ${index} matches OVITO`, () => {
            const frame = readFrame(file, {
                frame: index,
                includeIds: true,
                properties: ['*']
            });

            assert.equal(frame.metadata.natoms, reference.natoms, 'atom count');
            assert.equal(frame.positions.length, reference.natoms * 3, 'position buffer length');

            reference.positions.forEach((position, atom) => {
                position.forEach((value, axis) => {
                    closeTo(frame.positions[atom * 3 + axis], value, `atom ${atom} axis ${axis}`);
                });
            });

            assert.deepEqual(Array.from(frame.types), reference.types, 'types');
            if (reference.ids) {
                assert.deepEqual(Array.from(frame.ids), reference.ids, 'ids');
            }

            if (reference.cell) {
                reference.cell.forEach((row, axis) => {
                    frame.metadata.cellVectors.forEach((vector, index) => {
                        closeTo(vector[axis], row[index], `cell vector ${index}, axis ${axis}`);
                    });
                    closeTo(frame.metadata.cellOrigin[axis], row[3], `cell origin axis ${axis}`);
                });

                const derived = cellFromBox(frame.metadata.boxBounds);
                ['a', 'b', 'c'].forEach((label, index) => {
                    derived[label].forEach((value, axis) => {
                        closeTo(value, frame.metadata.cellVectors[index][axis], `boxBounds ${label}[${axis}]`);
                    });
                });

                assert.deepEqual(frame.metadata.pbc, reference.pbc, 'periodic boundary flags');
            } else {
                assert.deepEqual(frame.metadata.pbc, [false, false, false],
                    'a file with no declared cell is not periodic');
            }

            for (const [ovitoName, values] of Object.entries(reference.properties)) {
                if (NOT_COMPARED[ovitoName]) continue;

                const ours = OVITO_PROPERTY_ALIASES[ovitoName] ?? ovitoName;
                const column = frame.properties?.[ours];
                assert.ok(column, `expected column ${ours} (OVITO's "${ovitoName}") to be present`);
                values.forEach((value, atom) => {
                    closeTo(column[atom], value, `${ours}[${atom}]`);
                });
            }
        });
    });

    if (isDataFile(name)) {
        test(`${name}: per-type masses cover what OVITO reports per particle`, () => {
            const frame = readFrame(file, { includeIds: true });
            const reference = expected.frames[0];
            const perParticleMass = reference.properties.Mass;

            if (!perParticleMass) return;

            assert.ok(frame.massesByType, 'massesByType');
            reference.types.forEach((type, atom) => {
                closeTo(frame.massesByType[type - 1], perParticleMass[atom], `mass of atom ${atom}`);
            });
        });
    }

    test(`${name}: readHeader agrees with readFrame about the same frame`, () => {
        const header = readHeader(file);
        const frame = readFrame(file);

        assert.equal(header.format, frame.metadata.format);
        assert.equal(header.timestep, frame.metadata.timestep);
        assert.equal(header.natoms, frame.metadata.natoms);
        assert.deepEqual(header.pbc, frame.metadata.pbc);
        assert.deepEqual(header.headers, frame.metadata.headers);

        if (expected.frames[0].cell) {
            assert.deepEqual(header.boxBounds, frame.metadata.boxBounds);
            assert.deepEqual(header.cellVectors, frame.metadata.cellVectors);
        } else {
            assert.deepEqual(header.cellVectors, [[0, 0, 0], [0, 0, 0], [0, 0, 0]]);
        }
    });
}

test('the oracle covers every fixture on disk', () => {
    const onDisk = fs.readdirSync(FIXTURES).filter((name) => !name.startsWith('.')).sort();
    assert.deepEqual(onDisk, Object.keys(oracle.fixtures).sort());
});
