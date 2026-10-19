# Copyright (C) 2026 Emmanuele Semprini
# SPDX-License-Identifier: AGPL-3.0-only
# Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

"""Create the existing self-validated planar-frame STL for FaceTools acceptance.
Python standard library only; no app, profiles, accounts or external libraries.
Writes new files only under repository build/face-extrusion-demo, never tools.
Existing fixtures are refused. --help does not generate files.
"""
from collections import defaultdict
from pathlib import Path
import argparse
import json
import math

argparse.ArgumentParser(description=__doc__).parse_args()
root = Path(__file__).resolve().parents[2]
folder = root / "build/face-extrusion-demo"
assert folder.resolve() == folder and not folder.is_symlink()
stl_path = folder / "planar-frame-demo.stl"
info_path = folder / "fixture-info.json"
if stl_path.exists() or info_path.exists():
    raise SystemExit("Refusing to overwrite an existing demo fixture.")
folder.mkdir(parents=True, exist_ok=True)

# 60 x 40 mm rectangular frame, 40 x 20 mm through-hole, 4 mm thick.
base = [(-30., -20., 0.), (30., -20., 0.), (30., 20., 0.), (-30., 20., 0.),
        (-20., -10., 0.), (20., -10., 0.), (20., 10., 0.), (-20., 10., 0.)]
cap = [(0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
       (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7)]
vertices = base + [(x, y, 4.) for x, y, z in base]
triangles = [(a, c, b) for a, b, c in cap]
triangles += [(a + 8, b + 8, c + 8) for a, b, c in cap]
patch_edges = defaultdict(list)
for face in cap:
    for a, b in zip(face, face[1:] + face[:1]):
        patch_edges[tuple(sorted((a, b)))].append((a, b))
for incidences in patch_edges.values():
    if len(incidences) == 1:
        a, b = incidences[0]
        triangles.extend(((a, b, b + 8), (a, b + 8, a + 8)))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


closed_edges = defaultdict(list)
volume = 0.
for face in triangles:
    a, b, c = (vertices[i] for i in face)
    n = cross(subtract(b, a), subtract(c, a))
    assert dot(n, n) > 0
    volume += dot(a, cross(b, c)) / 6.
    for p, q in zip(face, face[1:] + face[:1]):
        closed_edges[tuple(sorted((p, q)))].append((p, q))
assert all(len(e) == 2 and e[0] == e[1][::-1] for e in closed_edges.values())
assert len(vertices) - len(closed_edges) + len(triangles) == 0  # one through-hole
assert math.isclose(volume, 6400., rel_tol=1e-12)
assert len(triangles) == 32

lines = ["solid FaceTools_demo_frame"]
for face in triangles:
    a, b, c = (vertices[i] for i in face)
    normal = cross(subtract(b, a), subtract(c, a))
    length = math.sqrt(dot(normal, normal))
    normal = tuple(value / length for value in normal)
    lines.append("  facet normal " + " ".join(f"{value:.9g}" for value in normal))
    lines.append("    outer loop")
    for vertex in (a, b, c):
        lines.append("      vertex " + " ".join(f"{value:.9g}" for value in vertex))
    lines.extend(("    endloop", "  endfacet"))
lines.append("endsolid FaceTools_demo_frame")
with stl_path.open("x", encoding="ascii") as stream:
    stream.write("\n".join(lines) + "\n")
info = {
    "name": "FaceTools planar-frame demo",
    "units": "mm",
    "outer_size": [60, 40, 4],
    "through_hole_size": [40, 20],
    "triangles": len(triangles),
    "signed_volume_mm3": volume,
    "closed_oriented_edges": True,
    "euler_characteristic": 0,
    "manual_check": "Select the top ring at z=4; try +5 mm and -2 mm. The source must remain unchanged and the center hole must remain open."
}
with info_path.open("x", encoding="utf-8") as stream:
    stream.write(json.dumps(info, indent=2) + "\n")
print(f"Created {stl_path.name}: {len(triangles)} triangles; {volume:g} mm^3; one open through-hole.")
