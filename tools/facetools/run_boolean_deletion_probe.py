#!/usr/bin/env python3
# Copyright (C) 2026 Emmanuele Semprini
# SPDX-License-Identifier: AGPL-3.0-only
# Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

"""Build an isolated deletion probe using an existing, completed CMake build.

Never starts CMake, rebuilds/copies archives, or touches application profiles.
Do not run while the application/test build is writing those archives.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import struct
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(__file__).with_name("boolean_deletion_probe.cpp").resolve()
CALLER = ROOT / "src/slic3r/GUI/Gizmos/GLGizmoMeshBoolean.cpp"


def fixture_bytes():
    """Versioned input definition; byte-identical to the original 228-face STL."""
    xs, ys, zs = [0, 5, 15, 20, 35, 40, 50, 55], [0, 5, 15, 20], [0, 1, 2, 4]
    filled = set()
    for i in range(len(xs) - 1):
        for j in range(len(ys) - 1):
            x, y = (xs[i] + xs[i + 1]) / 2, (ys[j] + ys[j + 1]) / 2
            if x < 20:
                height = 0 if 5 < x < 15 and 5 < y < 15 else 2
            elif x >= 35:
                height = 0 if 40 < x < 50 and 5 < y < 15 else 4
            else:
                height = 1 if y < 5 else 0
            for k in range(len(zs) - 1):
                if (zs[k] + zs[k + 1]) / 2 < height:
                    filled.add((i, j, k))
    faces = [
        ((-1, 0, 0), ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0))),
        ((1, 0, 0), ((1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1))),
        ((0, -1, 0), ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1))),
        ((0, 1, 0), ((0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0))),
        ((0, 0, -1), ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0))),
        ((0, 0, 1), ((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))),
    ]
    records = []
    for i, j, k in sorted(filled):
        for normal, corners in faces:
            if (i + normal[0], j + normal[1], k + normal[2]) in filled:
                continue
            points = [(xs[i + a], ys[j + b], zs[k + c]) for a, b, c in corners]
            for triangle in ((0, 1, 2), (0, 2, 3)):
                coords = [value for index in triangle for value in points[index]]
                records.append(struct.pack("<12fH", *normal, *coords, 0))
    header = b"FaceTools: connected rings, top faces z=2 and z=4, two holes".ljust(80, b" ")
    data = header + struct.pack("<I", len(records)) + b"".join(records)
    expected = "50455ad3c4345a1d52cb3fb0a94e12546474589acd8c1d875fe01b26a077c313"
    if hashlib.sha256(data).hexdigest() != expected:
        raise RuntimeError("The versioned fixture no longer matches the original audit input.")
    return data


def without_outputs(arguments):
    result, skip = [], False
    for arg in arguments:
        if skip:
            skip = False
        elif arg in ("-o", "-MF", "-MT", "-MQ", "-MJ", "--serialize-diagnostics"):
            skip = True
        elif arg not in ("-MD", "-MMD", "-MP"):
            result.append(arg)
    return result


def run_logged(command, cwd, log, env=None):
    with log.open("x", encoding="utf-8") as stream:
        completed = subprocess.run(command, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT)
    print(f"{log.name}: exit {completed.returncode}", flush=True)
    return completed.returncode


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-dir", type=Path, default=Path("build/facetools-arm64"),
                        help="Completed Makefiles build, relative to repository root or absolute")
    parser.add_argument("--baseline", action="store_true",
                        help="Also run the unsorted baseline; completed failing scenarios are expected")
    args = parser.parse_args()
    build = (ROOT / args.build_dir).resolve()
    entries = json.loads((build / "compile_commands.json").read_text(encoding="utf-8"))
    entry = next(e for e in entries if e["file"].endswith("/src/slic3r/GUI/Gizmos/GLGizmoMeshBoolean.cpp"))
    if (Path(entry["directory"]) / entry["file"]).resolve() != CALLER.resolve():
        raise RuntimeError("Compile commands belong to another checkout; configure/build this checkout first.")
    production = CALLER.read_text(encoding="utf-8")
    function = production.index("void BooleanOperationEngine::delete_volumes_from_model(")
    start = production.index("        // ObjectList consumes these indices in reverse order.", function)
    end = production.index("        wxGetApp().obj_list()->delete_from_model_and_list(items);", start)
    snippet = production[start:end]
    if snippet.count("std::sort(") != 1 or snippet.count("std::unique(") != 1:
        raise RuntimeError("Production deletion code changed; review the verbatim extraction before running.")
    linkdir = build / "tests/libslic3r"
    link_args = shlex.split((linkdir / "CMakeFiles/face_extrusion_tests.dir/link.txt").read_text(encoding="utf-8"))
    output_root = ROOT / "build/facetools-probes"
    subprocess.run(["git", "check-ignore", "--quiet", "build/facetools-probes/ignore-check"], cwd=ROOT, check=True)
    output_root.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix="boolean-deletion-", dir=output_root))
    print(f"Artifacts: {output}", flush=True)
    (output / "production_delete_order.inc").write_text(snippet, encoding="utf-8")
    fixture = output / "two-level-rings-demo.stl"
    fixture.write_bytes(fixture_bytes())
    obj, executable = output / "boolean_deletion_probe.o", output / "boolean_deletion_probe"
    compile_args = without_outputs(entry.get("arguments") or shlex.split(entry["command"]))
    if compile_args.count(entry["file"]) != 1:
        raise RuntimeError("Expected exactly one source argument in the cached compiler command.")
    compile_args = [str(SOURCE) if arg == entry["file"] else arg for arg in compile_args]
    compile_args += ["-I" + str(ROOT / "tests"), "-I" + str(output), "-DCATCH_CONFIG_FAST_COMPILE",
                     "-DCATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS", "-o", str(obj)]
    if run_logged(compile_args, entry["directory"], output / "compile.log"):
        raise RuntimeError("Probe compilation failed; see compile.log in the artifact directory.")
    link = []
    replacements = 0
    for arg in without_outputs(link_args):
        if arg.endswith("/test_face_extrusion.cpp.o"):
            link.append(str(obj))
            replacements += 1
        elif not arg.endswith("/test_face_extrusion_model.cpp.o"):
            link.append(arg)
    if replacements != 1:
        raise RuntimeError("Unexpected test link layout; no link was attempted.")
    link += ["-o", str(executable)]
    if run_logged(link, linkdir, output / "link.log"):
        raise RuntimeError("Probe linking failed; see link.log in the artifact directory.")
    provenance = {
        "caller_sha256": hashlib.sha256(production.encode()).hexdigest(),
        "snippet_sha256": hashlib.sha256(snippet.encode()).hexdigest(),
        "fixture_sha256": hashlib.sha256(fixture.read_bytes()).hexdigest(),
        "probe_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "build_dir": str(build),
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    for apply_sort in ([False, True] if args.baseline else [True]):
        log = output / ("fixed.log" if apply_sort else "baseline.log")
        env = dict(os.environ, MP_APPLY_SORT=str(int(apply_sort)), FACETOOLS_PROBE_FIXTURE=str(fixture))
        code = run_logged([str(executable), "[multipart_delete_audit]", "-r", "compact"], output, log, env)
        text = log.read_text(encoding="utf-8")
        summaries = re.findall(r"(?:DELETE|MULTIOBJECT)_SUMMARY fixed=([01]) cases=(\d+) failed=(\d+)", text)
        complete = len(summaries) == 2 and [int(row[1]) for row in summaries] == [900, 240]
        failures = sum(int(row[2]) for row in summaries)
        if not complete or any(row[0] != str(int(apply_sort)) for row in summaries):
            raise RuntimeError(f"Incomplete scenario matrix in {log.name}; not a valid baseline/pass.")
        if apply_sort and (code != 0 or failures != 0):
            raise RuntimeError("Corrected deletion regression failed.")
        if not apply_sort and (code <= 0 or failures == 0):
            raise RuntimeError("Baseline did not demonstrate the expected assertion failures.")
        print(f"{'Fixed' if apply_sort else 'Baseline'}: 1140 scenarios, {failures} failures", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, RuntimeError, StopIteration, subprocess.CalledProcessError) as error:
        raise SystemExit(f"Probe stopped: {error}")
