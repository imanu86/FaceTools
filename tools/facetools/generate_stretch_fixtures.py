# Copyright (C) 2026 Emmanuele Semprini
# SPDX-License-Identifier: AGPL-3.0-only
# Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

"""Create deterministic, closed L/tube STL fixtures under ignored build/ only.

Refuses to overwrite existing files. No user model, profile or app is accessed.
The L is 40 x 10 x 24 mm, with a 4 mm arm and an underside plane at Z=20.
The tube is 20 x 20 x 30 mm, with a 10 x 10 mm through-hole.
"""
from collections import Counter
from pathlib import Path
import argparse
import hashlib
import json
import math
import struct


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def l_bracket():
    section = [(0, 0), (10, 0), (10, 20), (40, 20), (40, 24), (0, 24)]
    vertices = [(x, y, z) for y in (0, 10) for x, z in section]
    cap = [(0, 1, 2), (0, 2, 5), (2, 3, 4), (2, 4, 5)]
    faces = list(cap) + [(a + 6, c + 6, b + 6) for a, b, c in cap]
    for i in range(6):
        j = (i + 1) % 6
        faces.extend([(i, i + 6, j + 6), (i, j + 6, j)])
    return vertices, faces, 3600.0


def hollow_tube():
    section = [(-10, -10), (10, -10), (10, 10), (-10, 10),
               (-5, -5), (5, -5), (5, 5), (-5, 5)]
    vertices = [(x, y, z) for z in (0, 30) for x, y in section]
    faces = []
    for i in range(4):
        j = (i + 1) % 4
        cap = [(i, j, j + 4), (i, j + 4, i + 4)]
        faces.extend((a + 8, b + 8, c + 8) for a, b, c in cap)
        faces.extend((a, c, b) for a, b, c in cap)
        faces.extend([(i, j, j + 8), (i, j + 8, i + 8)])
        a, b = i + 4, j + 4
        faces.extend([(a, a + 8, b + 8), (a, b + 8, b)])
    return vertices, faces, 9000.0


def encode_checked(vertices, faces, expected_volume):
    edges = Counter()
    volume = 0.0
    out = bytearray(b'FaceTools deterministic geometry fixture'.ljust(80, b'\0'))
    out += struct.pack('<I', len(faces))
    for triangle in faces:
        assert len(set(triangle)) == 3
        a, b, c = (vertices[i] for i in triangle)
        normal = cross(sub(b, a), sub(c, a))
        norm = math.sqrt(dot(normal, normal))
        assert math.isfinite(norm) and norm > 0
        volume += dot(a, cross(b, c)) / 6.0
        for i in range(3):
            edges[(triangle[i], triangle[(i + 1) % 3])] += 1
        values = tuple(x / norm for x in normal) + a + b + c
        assert all(math.isfinite(x) for x in values)
        out += struct.pack('<12fH', *values, 0)
    assert all(count == 1 and edges[(b, a)] == 1 for (a, b), count in edges.items())
    assert math.isclose(volume, expected_volume, rel_tol=1e-12, abs_tol=1e-8)
    assert len(out) == 84 + 50 * len(faces)
    return bytes(out), {'vertices': len(vertices), 'triangles': len(faces),
                        'signed_volume_mm3': volume, 'closed_and_consistently_wound': True,
                        'sha256': hashlib.sha256(out).hexdigest()}


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    root = Path(__file__).resolve().parents[2]
    output = root / 'build/facetools-v4-demo'
    assert output.resolve() == output and not output.is_symlink()
    fixtures = {'l-bracket.stl': l_bracket(), 'hollow-tube.stl': hollow_tube()}
    paths = [output / name for name in fixtures] + [output / 'fixture-provenance.json']
    if any(path.exists() or path.is_symlink() for path in paths):
        raise SystemExit('Refusing to overwrite an existing fixture or provenance file.')
    encoded = {name: encode_checked(*data) for name, data in fixtures.items()}
    output.mkdir(parents=True, exist_ok=True)
    provenance = {}
    for name, (data, details) in encoded.items():
        with (output / name).open('xb') as stream:
            stream.write(data)
        provenance[name] = details
        print(f'{name}: {details["triangles"]} triangles, {details["signed_volume_mm3"]:g} mm3, closed/winding PASS')
    with (output / 'fixture-provenance.json').open('x', encoding='utf-8') as stream:
        json.dump(provenance, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
