// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_CutGuideSelection_hpp_
#define slic3r_CutGuideSelection_hpp_

#include "Point.hpp"

#include <cstddef>
#include <memory>
#include <vector>

namespace Slic3r {

class TriangleMesh;

struct CutGuidePatch {
    int plane_id = -1;
    std::vector<int> triangles;
    size_t seed_facet = 0;
    Vec3d hit_local = Vec3d::Zero();
    Vec3d normal_local = Vec3d::Zero();
};

struct CutGuideSelectionLimits {
    size_t max_source_facets = 2000000;
    size_t max_guides = 4096;
};

// Owns an immutable source snapshot and a private, unpainted native selector.
// Single-threaded: angular_region updates only its private seed-fill/cache state.
// Invalid input throws std::invalid_argument; exceeded limits throw
// std::length_error. Neither operation modifies source geometry or caller guides.
class CutGuideSelection {
public:
    CutGuideSelection(std::shared_ptr<const TriangleMesh> mesh,
                      std::vector<int> facet_plane_ids,
                      CutGuideSelectionLimits limits = {});
    ~CutGuideSelection();

    // Run before constructing native Measuring/TriangleSelector. All original
    // vertices/facets must fit native int indices; vertices are additionally
    // bounded by 3 * max_source_facets (including unused vertices). Reject float
    // normal arithmetic that underflows, overflows, or becomes nonfinite/zero.
    // Open meshes are allowed: this does not impose solid/manifold constraints.
    static void validate_source(const TriangleMesh &mesh, CutGuideSelectionLimits limits = {});

    CutGuideSelection(const CutGuideSelection &) = delete;
    CutGuideSelection &operator=(const CutGuideSelection &) = delete;

    // Uses the native local-angle/EPSILON/clamped-dot behavior, without painting,
    // clipping, triangle splitting, or topology reconstruction. The immutable
    // sorted mask is reused while angle is unchanged and seed remains inside it.
    std::shared_ptr<const std::vector<int>> angular_region(size_t seed,
                                                         const Vec3d &hit,
                                                         float angle_degrees);

    // Exact masks, not entire measurement planes. If any requested facet is new,
    // add the entire region; otherwise subtract it. Repartition the exact result
    // by plane ID and native adjacency. Input masks may be unordered but must be
    // nonempty, unique, in range, and mutually disjoint across current patches.
    // Seed/hit provenance is retained where possible; output normals always come
    // from original triangle vertices promoted to double before the cross product.
    std::vector<CutGuidePatch> toggle(const std::vector<CutGuidePatch> &current,
                                     const std::vector<int> &region,
                                     size_t clicked_seed,
                                     const Vec3d &clicked_hit) const;

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace Slic3r

#endif // slic3r_CutGuideSelection_hpp_
