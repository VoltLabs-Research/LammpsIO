#!/usr/bin/env python3
"""Writes a LAMMPS binary dump fixture, in the layout LAMMPS itself writes.

Plain python3 — no OVITO needed. Run it when the binary fixture needs regenerating:

    python3 test/make-binary-dump.py

The fixture is committed, so this exists to document the byte layout the reader is
written against as much as to produce the file. Little-endian, 64-bit bigint, format
revision 2, which is what a current LAMMPS build produces on x86-64.

Per frame:
    bigint  -len(magic)          negative, marks the post-2018 format
    bytes   magic                "DUMPCUSTOM"
    int32   endian               0x0001
    int32   format revision      0x0002
    bigint  ntimestep
    bigint  natoms
    int32   triclinic
    int32   boundary[3][2]
    double  bbox[3][2]
    double  tilt[3]              only when triclinic
    int32   size_one             values per atom
    int32   len(unit style); bytes unit style
    char    time flag; double time  (only when the flag is set)
    int32   len(columns); bytes columns
    int32   nchunk
    per chunk: int32 n; double[n]
"""
import os
import struct

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT = os.path.join(HERE, 'fixtures', 'dump-binary-2frames.bin')

MAGIC = b'DUMPCUSTOM'
ENDIAN = 0x0001
FORMAT_REVISION = 0x0002
COLUMNS = 'id type x y z c_pe'

# Two frames, each with the same four atoms moved along x. Column order matches COLUMNS.
FRAMES = [
    (0, [
        (1, 1, 1.0, 2.0, 3.0, -3.15),
        (2, 1, 5.0, 10.0, 15.0, -3.25),
        (3, 2, 9.0, 18.0, 27.0, -2.05),
        (4, 2, 0.5, 0.5, 0.5, -1.95),
    ]),
    (750, [
        (1, 1, 1.5, 2.5, 3.5, -3.10),
        (2, 1, 5.5, 10.5, 15.5, -3.20),
        (3, 2, 9.5, 18.5, 27.5, -2.00),
        (4, 2, 1.0, 1.0, 1.0, -1.90),
    ]),
]

BBOX = [(0.0, 10.0), (0.0, 20.0), (0.0, 30.0)]
# Triclinic, so the tilt block is present and the reader has to handle it.
TILT = (2.0, 0.5, 0.25)
TRICLINIC = 1
# 0 = periodic on both faces of every axis, which is LAMMPS's `p p p`.
BOUNDARY = [[0, 0], [0, 0], [0, 0]]


def frame_bytes(timestep, atoms):
    out = bytearray()
    out += struct.pack('<q', -len(MAGIC))
    out += MAGIC
    out += struct.pack('<i', ENDIAN)
    out += struct.pack('<i', FORMAT_REVISION)
    out += struct.pack('<q', timestep)
    out += struct.pack('<q', len(atoms))
    out += struct.pack('<i', TRICLINIC)
    for axis in BOUNDARY:
        out += struct.pack('<ii', *axis)
    for lo, hi in BBOX:
        out += struct.pack('<dd', lo, hi)
    out += struct.pack('<ddd', *TILT)

    size_one = len(COLUMNS.split())
    out += struct.pack('<i', size_one)

    unit_style = b'metal'
    out += struct.pack('<i', len(unit_style)) + unit_style
    out += struct.pack('<b', 0)  # no simulation time recorded
    columns = COLUMNS.encode()
    out += struct.pack('<i', len(columns)) + columns

    # Two chunks, to exercise the chunked body rather than assuming a single block.
    split = len(atoms) // 2
    chunks = [atoms[:split], atoms[split:]]
    out += struct.pack('<i', len(chunks))
    for chunk in chunks:
        values = [float(value) for atom in chunk for value in atom]
        out += struct.pack('<i', len(values))
        out += struct.pack(f'<{len(values)}d', *values)

    return bytes(out)


def main():
    payload = b''.join(frame_bytes(timestep, atoms) for timestep, atoms in FRAMES)
    with open(OUTPUT, 'wb') as handle:
        handle.write(payload)
    print(f'wrote {OUTPUT} ({len(payload)} bytes, {len(FRAMES)} frames)')


if __name__ == '__main__':
    main()
