#!/usr/bin/env python3
"""Regenerates test/oracle.json: ground truth for every fixture, read by OVITO itself.

Run under an OVITO Python interpreter, which is the only thing here that needs OVITO:

    OVITO_BIN=/path/to/ovitos npm run oracle

The output is committed so the differential tests run anywhere — CI has no OVITO. Rerun
this only when a fixture changes or a new one is added, and review the diff: a change in
oracle.json means the reference moved, which is a much bigger claim than a test edit.

STDOUT is reserved for progress messages; the data goes to the JSON file.
"""
import json
import os
import sys

from ovito.io import import_file

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(HERE, 'fixtures')
OUTPUT = os.path.join(HERE, 'oracle.json')

RENAMED = {
    'Position': None,
    'Particle Type': None,
    'Particle Identifier': None,
}


def frame_record(data):
    """The parts of a frame both implementations can be held to."""
    particles = data.particles
    record = {
        'natoms': int(particles.count),
        'positions': [[round(float(c), 5) for c in p] for p in particles.positions],
        'types': [int(t) for t in particles.particle_types],
        'cell': None,
        'pbc': None,
        'properties': {},
    }

    if data.cell is not None:
        record['cell'] = [[round(float(c), 5) for c in row] for row in data.cell[...]]
        record['pbc'] = [bool(flag) for flag in data.cell.pbc]

    if 'Particle Identifier' in particles:
        record['ids'] = [int(i) for i in particles['Particle Identifier']]

    for name in particles.keys():
        if name in RENAMED:
            continue
        values = particles[name]
        if values.ndim != 1:
            continue
        record['properties'][name] = [round(float(v), 5) for v in values]

    return record


def fixture_record(path):
    pipeline = import_file(path)
    frames = []
    for index in range(pipeline.source.num_frames):
        frames.append(frame_record(pipeline.compute(index)))
    return {
        'frameCount': int(pipeline.source.num_frames),
        'frames': frames,
    }


def main():
    if not os.path.isdir(FIXTURES):
        sys.exit(f'no fixtures directory at {FIXTURES}')

    names = sorted(n for n in os.listdir(FIXTURES) if not n.startswith('.'))
    if not names:
        sys.exit(f'no fixtures in {FIXTURES}')

    oracle = {'ovito': None, 'fixtures': {}}

    import ovito
    oracle['ovito'] = ovito.version_string

    for name in names:
        path = os.path.join(FIXTURES, name)
        print(f'  reading {name}', flush=True)
        try:
            oracle['fixtures'][name] = fixture_record(path)
        except Exception as exc:  # noqa: BLE001 - a bad fixture must name itself
            sys.exit(f'{name}: {exc}')

    with open(OUTPUT, 'w') as handle:
        json.dump(oracle, handle, indent=2, sort_keys=True)
        handle.write('\n')

    print(f'wrote {OUTPUT} from OVITO {oracle["ovito"]} ({len(names)} fixtures)')


if __name__ == '__main__':
    main()
