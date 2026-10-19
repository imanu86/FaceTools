// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include "CutGuideSelection.hpp"

#include "TriangleMesh.hpp"
#include "TriangleSelector.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>

namespace Slic3r {
namespace {

bool finite_local_point(const Vec3d &point)
{
    return point.allFinite() && point.cast<float>().allFinite();
}

Vec3d source_normal(const indexed_triangle_set &its, size_t facet)
{
    const Vec3i &triangle = its.indices[facet];
    const Vec3d a = its.vertices[triangle[0]].cast<double>();
    const Vec3d b = its.vertices[triangle[1]].cast<double>();
    const Vec3d c = its.vertices[triangle[2]].cast<double>();
    const Vec3d cross = (b - a).cross(c - a);
    const Vec3d scaled = cross / cross.cwiseAbs().maxCoeff();
    return Vec3d(scaled / scaled.norm());
}

} // namespace

void CutGuideSelection::validate_source(const TriangleMesh &mesh, CutGuideSelectionLimits limits)
{
    if (limits.max_source_facets == 0 || limits.max_guides == 0)
        throw std::invalid_argument("Cut guide selection limits must be positive.");
    const indexed_triangle_set &its = mesh.its;
    const size_t native_max = size_t(std::numeric_limits<int>::max());
    const uint64_t vertex_limit = 3ULL * std::min(limits.max_source_facets, native_max);
    if (its.indices.size() > limits.max_source_facets || its.indices.size() > native_max ||
        its.vertices.size() > native_max || uint64_t(its.vertices.size()) > vertex_limit)
        throw std::length_error("The cut guide source exceeds the native selector limits.");
    if (its.indices.empty() || its.vertices.empty())
        throw std::invalid_argument("The cut guide source is empty.");
    for (const Vec3f &vertex : its.vertices)
        if (!vertex.allFinite())
            throw std::invalid_argument("The cut guide source contains a nonfinite vertex.");
    for (const Vec3i &triangle : its.indices) {
        for (int index : triangle)
            if (index < 0 || size_t(index) >= its.vertices.size())
                throw std::invalid_argument("The cut guide source contains an invalid vertex index.");
        const Vec3f vertices[3] = {its.vertices[triangle[0]], its.vertices[triangle[1]], its.vertices[triangle[2]]};
        const Vec3d a = vertices[0].cast<double>();
        const Vec3d cross = (vertices[1].cast<double>() - a).cross(vertices[2].cast<double>() - a);
        if (!cross.allFinite() || !(cross.cwiseAbs().maxCoeff() > 0.0))
            throw std::invalid_argument("The cut guide source contains a degenerate triangle.");
        // Use the actual arithmetic consumed by native Measuring/TriangleSelector.
        // Reject unsafe representation, without changing their angle/EPSILON rules.
        const Vec3f native_cross = face_normal(vertices);
        const float squared_length = native_cross.squaredNorm();
        if (!native_cross.allFinite() || !std::isnormal(squared_length) || !(squared_length > 0.0f))
            throw std::invalid_argument("A source triangle has an unrepresentable native float normal.");
        const Vec3f native_normal = native_cross.normalized();
        if (!native_normal.allFinite() || !(native_normal.squaredNorm() > 0.0f))
            throw std::invalid_argument("A source triangle has an invalid native float normal.");
    }
}

struct CutGuideSelection::Impl {
    std::shared_ptr<const TriangleMesh> mesh;
    std::vector<int> plane_ids;
    CutGuideSelectionLimits limits;
    std::unique_ptr<TriangleSelector> selector;
    std::shared_ptr<const std::vector<int>> cached_region;
    float cached_angle = 0.0f;

    Impl(std::shared_ptr<const TriangleMesh> source, std::vector<int> ids, CutGuideSelectionLimits bounds)
        : mesh(std::move(source)), plane_ids(std::move(ids)), limits(bounds)
    {
        if (!mesh)
            throw std::invalid_argument("The cut guide source is unavailable.");
        CutGuideSelection::validate_source(*mesh, limits);
        const size_t count = mesh->its.indices.size();
        if (plane_ids.size() != count)
            throw std::invalid_argument("The cut guide plane map does not match the source.");
        for (int plane : plane_ids)
            if (plane < 0 || size_t(plane) >= count)
                throw std::invalid_argument("The cut guide plane map contains an invalid plane ID.");
        // No selector allocation, neighbor walk or native normal calculation precedes
        // source validation. This selector is never painted, split or deserialized.
        selector = std::make_unique<TriangleSelector>(*mesh);
        const auto &neighbors = selector->get_neighbors();
        if (neighbors.size() != count || selector->get_triangles().size() != count)
            throw std::logic_error("The native cut guide selector lost original facet identity.");
        for (const Vec3i &adjacent : neighbors)
            for (int facet : adjacent)
                if (facet < -1 || (facet >= 0 && size_t(facet) >= count))
                    throw std::logic_error("The native cut guide selector returned an invalid neighbor.");
    }
};

CutGuideSelection::CutGuideSelection(std::shared_ptr<const TriangleMesh> mesh,
                                   std::vector<int> facet_plane_ids,
                                   CutGuideSelectionLimits limits)
    : m_impl(std::make_unique<Impl>(std::move(mesh), std::move(facet_plane_ids), limits))
{}

CutGuideSelection::~CutGuideSelection() = default;

std::shared_ptr<const std::vector<int>> CutGuideSelection::angular_region(size_t seed,
                                                                        const Vec3d &hit,
                                                                        float angle_degrees)
{
    Impl &impl = *m_impl;
    const size_t count = impl.mesh->its.indices.size();
    if (seed >= count || !finite_local_point(hit) || !std::isfinite(angle_degrees) ||
        angle_degrees < 0.0f || angle_degrees > 90.0f)
        throw std::invalid_argument("The cut guide angular selection request is invalid.");
    if (impl.cached_region && impl.cached_angle == angle_degrees &&
        std::binary_search(impl.cached_region->begin(), impl.cached_region->end(), int(seed)))
        return impl.cached_region;

    // Force recomputation on a cache miss: native seed-fill otherwise reuses its
    // selected seed even when the requested angle changed. Do not project a rounded
    // GUI hit, change the native threshold, or reconstruct a selected mesh.
    impl.selector->seed_fill_select_triangles(hit.cast<float>(), int(seed), Transform3d::Identity(),
                                              TriangleSelector::ClippingPlane{}, angle_degrees, 0.0f, true);
    auto region = std::make_shared<std::vector<int>>();
    const auto &triangles = impl.selector->get_triangles();
    if (triangles.size() != count)
        throw std::logic_error("The native cut guide selector no longer contains only original facets.");
    for (size_t i = 0; i < count; ++i) {
        const auto &triangle = triangles[i];
        if (triangle.valid() && !triangle.is_split() && triangle.is_selected_by_seed_fill()) {
            if (triangle.source_triangle != int(i))
                throw std::logic_error("The native cut guide selector returned a nonoriginal facet.");
            region->push_back(int(i));
        }
    }
    if (!std::binary_search(region->begin(), region->end(), int(seed)))
        throw std::logic_error("The native cut guide region does not contain its seed.");
    impl.cached_region = std::move(region);
    impl.cached_angle = angle_degrees;
    return impl.cached_region;
}

std::vector<CutGuidePatch> CutGuideSelection::toggle(const std::vector<CutGuidePatch> &current,
                                                   const std::vector<int> &region,
                                                   size_t clicked_seed,
                                                   const Vec3d &clicked_hit) const
{
    const Impl &impl = *m_impl;
    const size_t count = impl.mesh->its.indices.size();
    if (current.size() > impl.limits.max_guides || current.size() > count)
        throw std::length_error("The current cut guide list exceeds the selection limits.");
    if (region.empty() || region.size() > count || clicked_seed >= count || !finite_local_point(clicked_hit))
        throw std::invalid_argument("The cut guide toggle request is invalid.");

    // Validate the entire transaction before computing any replacement. Owner IDs
    // retain prior guide order/metadata while membership remains an exact facet set.
    std::vector<int> owner(count, -1);
    size_t total = 0;
    for (size_t i = 0; i < current.size(); ++i) {
        const CutGuidePatch &patch = current[i];
        if (patch.triangles.empty() || patch.triangles.size() > count - total || patch.plane_id < 0 ||
            size_t(patch.plane_id) >= count || patch.seed_facet >= count || !finite_local_point(patch.hit_local) ||
            !patch.normal_local.allFinite() || !(patch.normal_local.cwiseAbs().maxCoeff() > 0.0))
            throw std::invalid_argument("A current cut guide has invalid mask or metadata.");
        total += patch.triangles.size();
        bool has_seed = false;
        for (int facet : patch.triangles) {
            if (facet < 0 || size_t(facet) >= count || owner[facet] != -1 || impl.plane_ids[facet] != patch.plane_id)
                throw std::invalid_argument("Current cut guide masks overlap or contain invalid facets.");
            owner[facet] = int(i);
            has_seed = has_seed || size_t(facet) == patch.seed_facet;
        }
        if (!has_seed)
            throw std::invalid_argument("A current cut guide seed is outside its mask.");
    }

    std::vector<unsigned char> selected(count, 0);
    bool add = false;
    for (int facet : region) {
        if (facet < 0 || size_t(facet) >= count || selected[facet])
            throw std::invalid_argument("The requested cut guide mask contains invalid or duplicate facets.");
        selected[facet] = 1;
        add = add || owner[facet] == -1;
    }
    if (!selected[clicked_seed])
        throw std::invalid_argument("The clicked cut guide seed is outside the requested mask.");
    for (size_t i = 0; i < count; ++i)
        selected[i] = add ? (owner[i] >= 0 || selected[i]) : (owner[i] >= 0 && !selected[i]);

    struct Component {
        int plane_id;
        size_t first_old;
        size_t seed_old;
        std::vector<int> triangles;
    };
    std::vector<Component> components;
    std::vector<int> component_of(count, -1);
    std::vector<int> queue;
    const auto &neighbors = impl.selector->get_neighbors();
    for (size_t root = 0; root < count; ++root) {
        if (!selected[root] || component_of[root] >= 0) continue;
        if (components.size() >= impl.limits.max_guides)
            throw std::length_error("The exact cut guide result exceeds the guide limit.");
        const int component = int(components.size());
        components.push_back({impl.plane_ids[root], current.size(), current.size(), {}});
        Component &out = components.back();
        queue.clear();
        queue.push_back(int(root));
        component_of[root] = component;
        for (size_t head = 0; head < queue.size(); ++head) {
            const int facet = queue[head];
            if (owner[facet] >= 0) out.first_old = std::min(out.first_old, size_t(owner[facet]));
            for (int adjacent : neighbors[facet]) {
                if (adjacent >= 0 && selected[adjacent] && component_of[adjacent] < 0 &&
                    impl.plane_ids[adjacent] == out.plane_id) {
                    component_of[adjacent] = component;
                    queue.push_back(adjacent);
                }
            }
        }
    }
    // Enumerate original indices in order rather than sorting/copying a native
    // selected mesh. Unselected facets of the same measured plane stay absent.
    for (size_t i = 0; i < count; ++i)
        if (component_of[i] >= 0) components[component_of[i]].triangles.push_back(int(i));
    for (size_t i = 0; i < current.size(); ++i) {
        const int component = component_of[current[i].seed_facet];
        if (component >= 0 && components[component].seed_old == current.size())
            components[component].seed_old = i;
    }
    std::vector<size_t> order(components.size());
    std::iota(order.begin(), order.end(), size_t(0));
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        const Component &left = components[a], &right = components[b];
        return std::make_tuple(left.first_old, left.seed_old == current.size(), left.triangles.front()) <
               std::make_tuple(right.first_old, right.seed_old == current.size(), right.triangles.front());
    });

    std::vector<CutGuidePatch> result;
    result.reserve(components.size());
    for (size_t index : order) {
        Component &component = components[index];
        CutGuidePatch patch;
        patch.plane_id = component.plane_id;
        patch.triangles = std::move(component.triangles);
        if (component.seed_old < current.size()) {
            const CutGuidePatch &previous = current[component.seed_old];
            patch.seed_facet = previous.seed_facet;
            patch.hit_local = previous.hit_local;
        } else if (component_of[clicked_seed] == int(index)) {
            patch.seed_facet = clicked_seed;
            patch.hit_local = clicked_hit;
        } else {
            patch.seed_facet = size_t(patch.triangles.front());
            patch.hit_local = impl.mesh->its.vertices[impl.mesh->its.indices[patch.seed_facet][0]].cast<double>();
        }
        patch.normal_local = source_normal(impl.mesh->its, patch.seed_facet);
        result.push_back(std::move(patch));
    }
    return result;
}

} // namespace Slic3r
