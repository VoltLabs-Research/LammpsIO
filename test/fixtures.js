'use strict';

const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');

const DIRECTORY = fs.mkdtempSync(path.join(os.tmpdir(), 'lammpsio-test-'));

const write = (name, body) => {
    const file = path.join(DIRECTORY, name);
    fs.writeFileSync(file, body);
    return file;
};

const dumpFrame = ({ timestep, bounds, boundaryFlags = 'pp pp pp', columns, rows }) => [
    'ITEM: TIMESTEP', String(timestep),
    'ITEM: NUMBER OF ATOMS', String(rows.length),
    `ITEM: BOX BOUNDS ${boundaryFlags}`,
    ...bounds,
    `ITEM: ATOMS ${columns}`,
    ...rows,
    ''
].join('\n');

const ORTHO_BOUNDS = ['0.0 10.0', '0.0 20.0', '0.0 30.0'];

const writeDump = (name, frames) => write(name, frames.map(dumpFrame).join(''));

module.exports = {
    DIRECTORY,
    ORTHO_BOUNDS,
    dumpFrame,
    write,
    writeDump
};
