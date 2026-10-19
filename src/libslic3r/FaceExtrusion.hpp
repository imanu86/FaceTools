// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_FaceExtrusion_hpp_
#define slic3r_FaceExtrusion_hpp_

#include "TriangleMesh.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace Slic3r {

enum class FaceExtrusionError {
    None,
    EmptySelection,
    InvalidTriangleIndex,
    InvalidVertex,
    DuplicateTriangle,
    DegenerateTriangle,
    InvalidTolerance,
    InvalidDisplacement,
    NonPlanar,
    Disconnected,
    NonManifold,
    InconsistentWinding,
    SelfIntersection,
    PrecisionLoss,
    TooLarge
};

struct FaceExtrusionOptions {
    double planarity_tolerance_mm = 1e-4;
    double length_tolerance_mm = 1e-6;
    double angular_tolerance_rad = 1e-3;
};

struct FaceExtrusionResult {
    indexed_triangle_set mesh;
    FaceExtrusionError error = FaceExtrusionError::None;
    std::string message;
    size_t boundary_edges = 0;
    double max_planarity_error_mm = 0.0;

    explicit operator bool() const noexcept { return error == FaceExtrusionError::None; }
};

// Create a separate, closed solid from ONE edge-connected planar patch.
// Both source vertices and displacement must already be in the SAME world-mm
// coordinate system. The caller applies instance/volume transforms (including
// handedness correction), and converts the result back to its storage frame.
// Displacement must be perpendicular to the patch; either sign is supported.
//
// Only the selected faces are inspected. Exact coincident vertices are welded;
// nearby vertices are never merged. Source cap positions/triangulation are
// preserved (planar within the given tolerance, not silently flattened).
// Boundary edges alone receive walls, including the boundaries of holes.
// Malformed, disconnected, non-manifold or overlapping patches are rejected.
// Input and its face properties are never modified. Output properties are reset.
// On failure mesh is empty; message is an English diagnostic for the caller.
FaceExtrusionResult extrude_planar_face(
    const indexed_triangle_set &source,
    const std::vector<int> &triangle_ids,
    const Vec3d &displacement,
    const FaceExtrusionOptions &options = {});

} // namespace Slic3r

#endif // slic3r_FaceExtrusion_hpp_
