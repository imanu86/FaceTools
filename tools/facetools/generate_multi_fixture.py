# Copyright (C) 2026 Emmanuele Semprini
# SPDX-License-Identifier: AGPL-3.0-only
# Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

"""Create a connected, closed two-hole fixture for same-facing multiselection.
Python standard library only; no app, profiles, accounts or external libraries.
Writes only new files under repository build/face-extrusion-demo, never tools.
The original local-host path assertion is retained intentionally.
Existing fixtures are refused. --help does not generate files.
"""
from collections import Counter, defaultdict, deque
from pathlib import Path
import argparse
import hashlib
import json
import math
import struct

argparse.ArgumentParser(description=__doc__).parse_args()
root = Path(__file__).resolve().parents[2]
base = root / "build/face-extrusion-demo"
assert base.resolve() == base and not base.is_symlink()
output = base / 'two-level-rings-demo.stl'
report = base / 'two-level-rings-info.json'
if output.exists() or report.exists():
    raise SystemExit('Refusing to overwrite an existing fixture or report.')
base.mkdir(parents=True, exist_ok=True)

xs = [0, 5, 15, 20, 35, 40, 50, 55]
ys = [0, 5, 15, 20]
zs = [0, 1, 2, 4]
filled = set()
for i in range(len(xs) - 1):
    for j in range(len(ys) - 1):
        x = (xs[i] + xs[i + 1]) / 2
        y = (ys[j] + ys[j + 1]) / 2
        if x < 20:
            height = 0 if 5 < x < 15 and 5 < y < 15 else 2
        elif x >= 35:
            height = 0 if 40 < x < 50 and 5 < y < 15 else 4
        else:
            height = 1 if y < 5 else 0
        for k in range(len(zs) - 1):
            if (zs[k] + zs[k + 1]) / 2 < height:
                filled.add((i, j, k))

# Outward winding; common lattice subdivisions prevent hanging/T-junction edges.
faces = [
    ((-1, 0, 0), ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0))),
    ((1, 0, 0), ((1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1))),
    ((0, -1, 0), ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1))),
    ((0, 1, 0), ((0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0))),
    ((0, 0, -1), ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0))),
    ((0, 0, 1), ((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))),
]
vertices = []
indices = {}
triangles = []
normals = []
for i, j, k in sorted(filled):
    for normal, corners in faces:
        neighbor = (i + normal[0], j + normal[1], k + normal[2])
        if neighbor in filled:
            continue
        quad = []
        for dx, dy, dz in corners:
            point = (xs[i + dx], ys[j + dy], zs[k + dz])
            if point not in indices:
                indices[point] = len(vertices)
                vertices.append(point)
            quad.append(indices[point])
        for triangle in ((quad[0], quad[1], quad[2]), (quad[0], quad[2], quad[3])):
            triangles.append(triangle)
            normals.append(normal)

edges = defaultdict(list)
directed = Counter()
for face_id, triangle in enumerate(triangles):
    for a, b in zip(triangle, triangle[1:] + triangle[:1]):
        edges[tuple(sorted((a, b)))].append(face_id)
        directed[(a, b)] += 1
assert all(len(owners) == 2 for owners in edges.values())
assert all(directed[(a, b)] == 1 and directed[(b, a)] == 1 for a, b in edges)
adjacency = defaultdict(set)
for owners in edges.values():
    a, b = owners
    adjacency[a].add(b)
    adjacency[b].add(a)
visited = {0}
queue = deque([0])
while queue:
    for neighbor in adjacency[queue.popleft()]:
        if neighbor not in visited:
            visited.add(neighbor)
            queue.append(neighbor)
assert len(visited) == len(triangles)
euler = len(vertices) - len(edges) + len(triangles)
assert euler == -2  # One connected closed surface, genus two.

def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])

volumes = []
area_by_top = defaultdict(float)
for tri, normal in zip(triangles, normals):
    a, b, c = (vertices[index] for index in tri)
    bc = cross(b, c)
    volumes.append(sum(a[axis] * bc[axis] for axis in range(3)) / 6)
    ab, ac = tuple(b[q] - a[q] for q in range(3)), tuple(c[q] - a[q] for q in range(3))
    geometric_normal = cross(ab, ac)
    assert sum(geometric_normal[q] * normal[q] for q in range(3)) > 0
    if normal == (0, 0, 1):
        area_by_top[a[2]] += math.sqrt(sum(v * v for v in geometric_normal)) / 2
volume = math.fsum(volumes)
assert abs(volume - 1875.0) < 1e-9
assert dict(area_by_top) == {1: 75.0, 2: 300.0, 4: 300.0}

payload = bytearray(b'FaceTools: connected rings, top faces z=2 and z=4, two holes'.ljust(80, b' '))
payload.extend(struct.pack('<I', len(triangles)))
for triangle, normal in zip(triangles, normals):
    values = list(normal) + [coordinate for index in triangle for coordinate in vertices[index]]
    payload.extend(struct.pack('<12fH', *values, 0))
with output.open('xb') as stream:
    stream.write(payload)
metadata = {
    'file': output.name,
    'sha256': hashlib.sha256(payload).hexdigest(),
    'bounds_mm': [[0, 0, 0], [55, 20, 4]],
    'vertices_welded': len(vertices),
    'triangles': len(triangles),
    'closed_edge_manifold': True,
    'consistent_winding': True,
    'surface_components': 1,
    'euler_characteristic': euler,
    'genus': 2,
    'volume_mm3': volume,
    'upward_face_areas_by_height_mm': dict(area_by_top),
    'instructions': 'Select the ring at z=2, then Shift-click the ring at z=4. Both normals point up; the two holes stay open. The lower bridge top at z=1 is a third optional face.',
}
with report.open('x', encoding='utf-8') as stream:
    json.dump(metadata, stream, indent=2)
    stream.write('\n')
print(json.dumps(metadata, indent=2))
