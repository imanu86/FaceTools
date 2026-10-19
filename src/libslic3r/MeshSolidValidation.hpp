// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_MeshSolidValidation_hpp_
#define slic3r_MeshSolidValidation_hpp_

#include "TriangleMesh.hpp"

#include <functional>

namespace Slic3r {

// Validate a non-empty, finite, representable closed solid without changing it.
// Rejects degenerate/duplicate faces, unused vertices, non-manifold topology,
// self-intersections and inconsistent shell nesting/orientation. Properly inward
// cavity shells and disjoint outward solids are allowed. No welding, snapping,
// repair or orientation changes are performed, including on failure.
//
// Invalid geometry throws std::runtime_error. check_cancel is called at bounded
// loop checkpoints and between the non-interruptible diagnostics/CGAL calls;
// any exception from it propagates unchanged to the caller. Empty Boolean-set
// identities must be handled by their adapter before calling this non-empty
// solid validator. Surface validity is not a proof of Boolean set semantics.
void validate_solid_mesh(const indexed_triangle_set &mesh,
                         const std::function<void()> &check_cancel = {});

} // namespace Slic3r

#endif // slic3r_MeshSolidValidation_hpp_
