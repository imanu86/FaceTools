#include <catch2/catch.hpp>

#include "libslic3r/FaceExtrusion.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/MeshDiagnostics.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "slic3r/GUI/Gizmos/FaceExtrusionModelPolicy.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace Slic3r;
using Slic3r::GUI::prepare_face_extrusion_part_storage;

namespace {

void require_exact_mesh(const indexed_triangle_set &actual, const indexed_triangle_set &expected)
{
    REQUIRE(actual.vertices.size() == expected.vertices.size());
    REQUIRE(actual.indices.size() == expected.indices.size());
    REQUIRE(actual.properties.size() == expected.properties.size());
    for (size_t i = 0; i < actual.vertices.size(); ++i)
        CHECK((actual.vertices[i].array() == expected.vertices[i].array()).all());
    for (size_t i = 0; i < actual.indices.size(); ++i)
        CHECK((actual.indices[i].array() == expected.indices[i].array()).all());
    for (size_t i = 0; i < actual.properties.size(); ++i) {
        CHECK(actual.properties[i].type == expected.properties[i].type);
        CHECK(actual.properties[i].area == expected.properties[i].area);
    }
}

void require_closed_solid(const indexed_triangle_set &mesh)
{
    REQUIRE_FALSE(mesh.vertices.empty());
    REQUIRE_FALSE(mesh.indices.empty());
    std::set<std::array<int, 3>> faces;
    for (const Vec3f &vertex : mesh.vertices) REQUIRE(vertex.allFinite());
    for (const Vec3i &face : mesh.indices) {
        for (int i = 0; i < 3; ++i) {
            REQUIRE(face[i] >= 0);
            REQUIRE(size_t(face[i]) < mesh.vertices.size());
        }
        const Vec3d a = mesh.vertices[face[0]].cast<double>();
        const Vec3d b = mesh.vertices[face[1]].cast<double>();
        const Vec3d c = mesh.vertices[face[2]].cast<double>();
        REQUIRE((b - a).cross(c - a).squaredNorm() > 0.0);
        std::array<int, 3> key{face[0], face[1], face[2]};
        std::sort(key.begin(), key.end());
        REQUIRE(faces.insert(key).second);
    }
    const auto stats = its_mesh_diagnostics(mesh);
    REQUIRE(stats.open_edges == 0);
    REQUIRE(stats.non_manifold_edges == 0);
    REQUIRE(stats.non_manifold_vertices == 0);
    REQUIRE(stats.same_direction_edges == 0);
    REQUIRE_FALSE(stats.has_reversed_faces);
    REQUIRE(std::isfinite(its_volume(mesh)));
    REQUIRE(its_volume(mesh) > 0.f);
    auto cgal = MeshBoolean::cgal::triangle_mesh_to_cgal(mesh);
    REQUIRE_FALSE(MeshBoolean::cgal::does_self_intersect(*cgal));
    REQUIRE(MeshBoolean::cgal::does_bound_a_volume(*cgal));
}

long euler_characteristic(const indexed_triangle_set &mesh)
{
    std::set<int> vertices;
    std::set<std::pair<int, int>> edges;
    for (const Vec3i &face : mesh.indices)
        for (int i = 0; i < 3; ++i) {
            vertices.insert(face[i]);
            edges.insert(std::minmax(face[i], face[(i + 1) % 3]));
        }
    return long(vertices.size()) - long(edges.size()) + long(mesh.indices.size());
}

indexed_triangle_set square_patch()
{
    indexed_triangle_set mesh;
    mesh.vertices = {{-2.f, -2.f, 0.f}, {2.f, -2.f, 0.f}, {2.f, 2.f, 0.f}, {-2.f, 2.f, 0.f}};
    mesh.indices = {{0, 1, 2}, {0, 2, 3}};
    mesh.properties.assign(mesh.indices.size(), FaceProperty{eNormal, 0.25});
    return mesh;
}

indexed_triangle_set narrow_u_patch()
{
    // Synthetic fixture only: no user STL or repaired geometry belongs in this test.
    // D and G are distinct in world float storage but both round to x=64 after
    // the real world-to-object transform. No individual triangle loses area.
    const float g = 1.5e-6f;
    indexed_triangle_set mesh;
    mesh.vertices = {{-2.f, -2.f, 0.f}, {2.f, -2.f, 0.f}, {2.f, 2.f, 0.f}, {g, 1.f, 0.f},
                     {1.f, 0.f, 0.f}, {-1.f, 0.f, 0.f}, {-g, 1.f, 0.f}, {-2.f, 2.f, 0.f}};
    mesh.indices = {{0, 1, 4}, {0, 4, 5}, {1, 2, 4}, {2, 3, 4}, {0, 5, 7}, {7, 5, 6}};
    return mesh;
}

indexed_triangle_set ring_patch()
{
    indexed_triangle_set mesh;
    mesh.vertices = {{-2.f, -2.f, 0.f}, {2.f, -2.f, 0.f}, {2.f, 2.f, 0.f}, {-2.f, 2.f, 0.f},
                     {-1.f, -1.f, 0.f}, {1.f, -1.f, 0.f}, {1.f, 1.f, 0.f}, {-1.f, 1.f, 0.f}};
    mesh.indices = {{0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                    {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    return mesh;
}

indexed_triangle_set cavity_solid(bool orient_cavity_inward, bool include_inner_island = false)
{
    auto mesh = make_cube(10, 10, 10).its;
    TriangleMesh inner = make_cube(6, 6, 6);
    inner.translate(2, 2, 2);
    if (orient_cavity_inward) inner.flip_triangles();
    its_merge(mesh, inner.its);
    if (include_inner_island) {
        TriangleMesh island = make_cube(2, 2, 2);
        island.translate(4, 4, 4);
        its_merge(mesh, island.its);
    }
    return mesh;
}

indexed_triangle_set vertex_pinched_solid()
{
    // Two outward tetrahedra share ONLY vertex 0. Every edge has two opposite
    // halfedges, but the incident face link at vertex 0 has two components.
    indexed_triangle_set mesh;
    mesh.vertices = {{0.f, 0.f, 0.f}, {1.f, 0.f, 0.f}, {0.f, 1.f, 0.f}, {0.f, 0.f, 1.f},
                     {-1.f, 0.f, 0.f}, {0.f, -1.f, 0.f}, {0.f, 0.f, -1.f}};
    mesh.indices = {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3},
                    {0, 4, 5}, {0, 6, 4}, {0, 5, 6}, {4, 6, 5}};
    return mesh;
}

void require_boolean_rejected_unchanged(const indexed_triangle_set &invalid)
{
    // These regressions exercise the enabled production wrapper, never legacy
    // MCUT or a silently skipped backend. No invalid TriangleMesh import/repair.
    REQUIRE(MeshBoolean::openmeshcraft::is_enabled());
    TriangleMesh source;
    source.its = invalid;
    const TriangleMesh empty;
    const TriangleMesh sentinel = make_cube(10, 10, 10);
    const auto source_before = source.its;
    const auto empty_before = empty.its;
    const auto sentinel_before = sentinel.its;
    std::vector<TriangleMesh> output{sentinel};
    int failures = 0;
    float progress = 0.f;
    REQUIRE_NOTHROW(MeshBoolean::mcut::make_boolean(source, empty, output, "UNION", [] { return false; },
        [&](float p) { progress = p; }, [&] { ++failures; }));
    CHECK(failures == 1);
    CHECK(progress < 100.f);
    REQUIRE(output.size() == 1);
    require_exact_mesh(output.front().its, sentinel_before);
    require_exact_mesh(source.its, source_before);
    require_exact_mesh(empty.its, empty_before);
}

void require_boolean_orphan_rejection(const Vec3f &unused)
{
    auto mesh = make_cube(10, 10, 10).its;
    mesh.properties.assign(mesh.indices.size(), FaceProperty{eNormal, 0.25});
    mesh.vertices.push_back(unused);
    require_boolean_rejected_unchanged(mesh);
}

void require_cgal_orphan_rejection(const Vec3f &unused)
{
    auto mesh = make_cube(10, 10, 10).its;
    mesh.vertices.push_back(unused);
    const auto before = mesh;
    auto cgal = MeshBoolean::cgal::triangle_mesh_to_cgal(mesh);
    const auto cgal_before = MeshBoolean::cgal::cgal_to_indexed_triangle_set(*cgal);
    REQUIRE(cgal_before.vertices.size() == mesh.vertices.size());
    bool bounds = true;
    REQUIRE_NOTHROW(bounds = MeshBoolean::cgal::does_bound_a_volume(*cgal));
    CHECK_FALSE(bounds);
    require_exact_mesh(mesh, before);
    // The safety guard must reject, not silently delete the isolated vertex.
    require_exact_mesh(MeshBoolean::cgal::cgal_to_indexed_triangle_set(*cgal), cgal_before);
}

} // namespace

TEST_CASE("Mesh quality rejects storage collapse below the positional tolerance",
          "[mesh_quality][mesh_quality_storage_reject]")
{
    const auto patch = narrow_u_patch();
    const auto patch_before = patch;
    const auto extrusion = extrude_planar_face(patch, {0, 1, 2, 3, 4, 5}, Vec3d(0, 0, 1));
    REQUIRE(extrusion);
    const TriangleMesh world{extrusion.mesh};
    const auto world_before = world.its;
    require_closed_solid(world.its);
    Transform3d instance = Transform3d::Identity();
    instance.translation().x() = -64.0;
    const auto prepared = prepare_face_extrusion_part_storage(world, instance);
    CHECK_FALSE(prepared.valid);
    CHECK_FALSE(bool(prepared));
    CHECK_FALSE(prepared.message.empty());
    CHECK(prepared.mesh.its.vertices.empty());
    CHECK(prepared.mesh.its.indices.empty());
    require_exact_mesh(world.its, world_before);
    require_exact_mesh(patch, patch_before);
}

TEST_CASE("Mesh quality preserves signed and mirrored valid NewPart storage",
          "[mesh_quality][mesh_quality_control]")
{
    const auto patch = square_patch();
    const auto patch_before = patch;
    for (double depth : {1.0, -1.0, 0.25, -0.25})
        for (double mirror : {1.0, -1.0}) {
            INFO("Signed depth=" << depth << ", mirror=" << mirror);
            const auto extrusion = extrude_planar_face(patch, {0, 1}, Vec3d(0, 0, depth));
            REQUIRE(extrusion);
            const TriangleMesh world{extrusion.mesh};
            const auto world_before = world.its;
            Transform3d instance = Transform3d::Identity();
            instance.linear() = Eigen::AngleAxisd(0.35, Vec3d::UnitZ()).toRotationMatrix() *
                                Vec3d(mirror * 1.25, 0.75, 2.0).asDiagonal();
            instance.translation() = Vec3d(8.0, -4.0, 2.0);
            const auto prepared = prepare_face_extrusion_part_storage(world, instance);
            REQUIRE(prepared.valid);
            REQUIRE(bool(prepared));
            CHECK(prepared.message.empty());
            require_closed_solid(prepared.mesh.its);
            TriangleMesh roundtrip = prepared.mesh;
            roundtrip.transform(instance, true);
            require_closed_solid(roundtrip.its);
            REQUIRE(roundtrip.its.vertices.size() == world.its.vertices.size());
            REQUIRE(roundtrip.its.indices.size() == world.its.indices.size());
            for (size_t i = 0; i < world.its.vertices.size(); ++i)
                CHECK((roundtrip.its.vertices[i].cast<double>() - world.its.vertices[i].cast<double>()).norm() <= 1e-4);
            for (size_t i = 0; i < world.its.indices.size(); ++i)
                CHECK((roundtrip.its.indices[i].array() == world.its.indices[i].array()).all());
            CHECK(its_volume(roundtrip.its) == Approx(16.0 * std::abs(depth)).epsilon(1e-5));
            require_exact_mesh(world.its, world_before);
            require_exact_mesh(patch, patch_before);
        }
}

TEST_CASE("Mesh quality preserves a through hole in NewPart storage",
          "[mesh_quality][mesh_quality_control]")
{
    const auto patch = ring_patch();
    const auto patch_before = patch;
    const auto extrusion = extrude_planar_face(patch, {0, 1, 2, 3, 4, 5, 6, 7}, Vec3d(0, 0, 1.5));
    REQUIRE(extrusion);
    const TriangleMesh world{extrusion.mesh};
    const auto world_before = world.its;
    const auto prepared = prepare_face_extrusion_part_storage(world, Transform3d::Identity());
    REQUIRE(prepared.valid);
    require_closed_solid(prepared.mesh.its);
    CHECK(euler_characteristic(prepared.mesh.its) == 0);
    CHECK(its_split(prepared.mesh.its).size() == 1);
    CHECK(its_volume(prepared.mesh.its) == Approx(18.0).epsilon(1e-5));
    require_exact_mesh(prepared.mesh.its, world_before);
    require_exact_mesh(world.its, world_before);
    require_exact_mesh(patch, patch_before);
}

TEST_CASE("Mesh quality preserves cavity orientation and nested solid islands",
          "[mesh_quality][mesh_quality_control]")
{
    for (bool island : {false, true}) {
        INFO("Nested positive island=" << island);
        const TriangleMesh world{cavity_solid(true, island)};
        const auto before = world.its;
        require_closed_solid(world.its);
        const auto prepared = prepare_face_extrusion_part_storage(world, Transform3d::Identity());
        REQUIRE(prepared.valid);
        require_closed_solid(prepared.mesh.its);
        std::vector<float> shell_volumes;
        for (const auto &shell : its_split(prepared.mesh.its)) shell_volumes.push_back(its_volume(shell));
        std::sort(shell_volumes.begin(), shell_volumes.end());
        REQUIRE(shell_volumes.size() == (island ? 3 : 2));
        CHECK(shell_volumes.front() == Approx(-216.0).margin(1e-3));
        CHECK(shell_volumes.back() == Approx(1000.0).margin(1e-3));
        if (island) CHECK(shell_volumes[1] == Approx(8.0).margin(1e-3));
        CHECK(its_volume(prepared.mesh.its) == Approx(island ? 792.0 : 784.0).margin(1e-3));
        require_exact_mesh(prepared.mesh.its, before);
        require_exact_mesh(world.its, before);
    }
}

TEST_CASE("Mesh quality permits disconnected positive solid shells",
          "[mesh_quality][mesh_quality_control]")
{
    auto mesh = make_cube(2, 2, 2).its;
    TriangleMesh second = make_cube(2, 2, 2);
    second.translate(4, 0, 0);
    its_merge(mesh, second.its);
    const TriangleMesh world{mesh};
    const auto prepared = prepare_face_extrusion_part_storage(world, Transform3d::Identity());
    REQUIRE(prepared.valid);
    require_closed_solid(prepared.mesh.its);
    CHECK(its_split(prepared.mesh.its).size() == 2);
    CHECK(its_volume(prepared.mesh.its) == Approx(16.0).margin(1e-4));
    require_exact_mesh(prepared.mesh.its, mesh);
    require_exact_mesh(world.its, mesh);
}

TEST_CASE("Mesh quality refuses empty NewPart storage without fabricating geometry",
          "[mesh_quality][mesh_quality_control]")
{
    const TriangleMesh empty;
    const auto before = empty.its;
    const auto prepared = prepare_face_extrusion_part_storage(empty, Transform3d::Identity());
    CHECK_FALSE(prepared.valid);
    CHECK_FALSE(prepared.message.empty());
    CHECK(prepared.mesh.its.vertices.empty());
    CHECK(prepared.mesh.its.indices.empty());
    require_exact_mesh(empty.its, before);
}

TEST_CASE("Mesh quality Boolean rejects invalid indexed solids transactionally",
          "[mesh_quality][mesh_quality_control]")
{
    const auto cube = make_cube(10, 10, 10).its;
    using Corrupt = std::function<void(indexed_triangle_set &)>;
    const std::vector<std::pair<std::string, Corrupt>> corruptions{
        {"naked edge", [](auto &mesh) { mesh.indices.pop_back(); }},
        {"same-winding duplicate", [](auto &mesh) { mesh.indices.push_back(mesh.indices.front()); }},
        {"opposite-winding duplicate", [](auto &mesh) {
            const auto face = mesh.indices.front(); mesh.indices.emplace_back(face[0], face[2], face[1]);
        }},
        {"zero area with distinct indices", [](auto &mesh) {
            mesh.vertices[mesh.indices[0][0]] = mesh.vertices[mesh.indices[0][1]];
        }},
        {"single inverted face", [](auto &mesh) { std::swap(mesh.indices[0][0], mesh.indices[0][1]); }},
        {"inward solid", [](auto &mesh) { its_flip_triangles(mesh); }},
        {"self-intersecting shells", [](auto &mesh) {
            TriangleMesh overlap = make_cube(10, 10, 10); overlap.translate(5, 5, 5); its_merge(mesh, overlap.its);
        }}
    };
    for (const auto &corrupt : corruptions) {
        INFO(corrupt.first);
        auto mesh = cube;
        corrupt.second(mesh);
        require_boolean_rejected_unchanged(mesh);
    }
    SECTION("nested shells with the same orientation") {
        require_boolean_rejected_unchanged(cavity_solid(false));
    }
    SECTION("vertex fan pinched without a non-manifold edge") {
        const auto mesh = vertex_pinched_solid();
        const auto stats = its_mesh_diagnostics(mesh);
        REQUIRE(stats.open_edges == 0);
        REQUIRE(stats.non_manifold_edges == 0);
        REQUIRE(stats.same_direction_edges == 0);
        REQUIRE(stats.non_manifold_vertices == 1);
        require_boolean_rejected_unchanged(mesh);
    }
    SECTION("isolated coincident opposite-winding triangle pair") {
        auto mesh = cube;
        const int a = int(mesh.vertices.size());
        mesh.vertices.emplace_back(20.f, 0.f, 0.f);
        mesh.vertices.emplace_back(21.f, 0.f, 0.f);
        mesh.vertices.emplace_back(20.f, 1.f, 0.f);
        mesh.indices.emplace_back(a, a + 1, a + 2);
        mesh.indices.emplace_back(a, a + 2, a + 1);
        const auto stats = its_mesh_diagnostics(mesh);
        REQUIRE(stats.open_edges == 0);
        REQUIRE(stats.non_manifold_edges == 0);
        REQUIRE(stats.non_manifold_vertices == 0);
        REQUIRE(stats.same_direction_edges == 0);
        require_boolean_rejected_unchanged(mesh);
    }
}

TEST_CASE("Mesh quality Boolean preserves empty set identities",
          "[mesh_quality][mesh_quality_control]")
{
    REQUIRE(MeshBoolean::openmeshcraft::is_enabled());
    const TriangleMesh cube = make_cube(2, 2, 2);
    const TriangleMesh empty;
    const auto cube_before = cube.its;
    const auto empty_before = empty.its;
    for (const std::string operation : {"UNION", "INTERSECTION", "A_NOT_B", "B_NOT_A"})
        for (bool a_empty : {false, true})
            for (bool b_empty : {false, true}) {
                if (!a_empty && !b_empty) continue;
                INFO(operation << ", empty A=" << a_empty << ", empty B=" << b_empty);
                const TriangleMesh &a = a_empty ? empty : cube;
                const TriangleMesh &b = b_empty ? empty : cube;
                const bool nonempty = operation == "UNION" ? (!a_empty || !b_empty) :
                    operation == "A_NOT_B" ? !a_empty && b_empty :
                    operation == "B_NOT_A" ? a_empty && !b_empty : false;
                std::vector<TriangleMesh> output{cube};
                int failures = 0;
                float progress = 0.f;
                REQUIRE_NOTHROW(MeshBoolean::mcut::make_boolean(a, b, output, operation, [] { return false; },
                    [&](float p) { progress = p; }, [&] { ++failures; }));
                CHECK(failures == 0);
                CHECK(progress == 100.f);
                REQUIRE(output.size() == (nonempty ? 2 : 1));
                require_exact_mesh(output.front().its, cube_before);
                if (nonempty) require_exact_mesh(output.back().its, cube_before);
                require_exact_mesh(cube.its, cube_before);
                require_exact_mesh(empty.its, empty_before);
            }
}

// Run each orphan test in its own subprocess for the pre-fix RED baseline:
// the old production CGAL predicate can fault before Catch2 gets control back.
TEST_CASE("Mesh quality Boolean safely rejects an unused origin vertex",
          "[mesh_quality][mesh_quality_orphan_boolean][mesh_quality_orphan_origin]")
{
    require_boolean_orphan_rejection(Vec3f(0, 0, 0));
}

TEST_CASE("Mesh quality Boolean safely rejects an unused inside vertex",
          "[mesh_quality][mesh_quality_orphan_boolean][mesh_quality_orphan_inside]")
{
    require_boolean_orphan_rejection(Vec3f(5, 5, 5));
}

TEST_CASE("Mesh quality Boolean safely rejects an unused outside vertex",
          "[mesh_quality][mesh_quality_orphan_boolean][mesh_quality_orphan_outside]")
{
    require_boolean_orphan_rejection(Vec3f(20, 20, 20));
}

TEST_CASE("Mesh quality CGAL volume guard safely rejects an unused origin vertex",
          "[mesh_quality][mesh_quality_orphan_cgal][mesh_quality_orphan_origin]")
{
    require_cgal_orphan_rejection(Vec3f(0, 0, 0));
}

TEST_CASE("Mesh quality CGAL volume guard safely rejects an unused inside vertex",
          "[mesh_quality][mesh_quality_orphan_cgal][mesh_quality_orphan_inside]")
{
    require_cgal_orphan_rejection(Vec3f(5, 5, 5));
}

TEST_CASE("Mesh quality CGAL volume guard safely rejects an unused outside vertex",
          "[mesh_quality][mesh_quality_orphan_cgal][mesh_quality_orphan_outside]")
{
    require_cgal_orphan_rejection(Vec3f(20, 20, 20));
}

namespace Slic3r { std::vector<Vec3d> its_reversed_face_ray_directions(); }

TEST_CASE("Reversed-face ray set avoids axis, diagonal and symmetry planes", "[mesh_quality][reversed_faces]")
{
    // Regression: axis-aligned rays on mirror-symmetric parts could run inside a
    // symmetry plane / along triangle edges, drop a crossing and make the
    // outermost hit a back face, flagging a valid solid as reversed (Cut and
    // extend then failed at specific distances with a Boolean-validity error).
    const std::vector<Vec3d> dirs = Slic3r::its_reversed_face_ray_directions();
    REQUIRE(dirs.size() == 46);
    for (size_t i = 0; i < dirs.size(); ++i) {
        CHECK(dirs[i].norm() == Approx(1.0));
        for (int a = 0; a < 3; ++a) {
            CHECK(std::abs(dirs[i][a]) > 1e-3);                   // not inside any axis plane
            for (int b = a + 1; b < 3; ++b)
                CHECK(std::abs(std::abs(dirs[i][a]) - std::abs(dirs[i][b])) > 1e-3); // not on a diagonal plane
        }
        for (size_t j = i + 1; j < dirs.size(); ++j)
            CHECK(dirs[i].dot(dirs[j]) < 1.0 - 1e-9);             // no duplicate rays (opposite rays are intended)
    }
}

TEST_CASE("Mirror-symmetric solids with raised ridges are never flagged as reversed", "[mesh_quality][reversed_faces]")
{
    const auto box = [](Vec3d lo, Vec3d hi) {
        indexed_triangle_set m = its_make_cube(hi.x() - lo.x(), hi.y() - lo.y(), hi.z() - lo.z());
        its_translate(m, lo.cast<float>());
        return TriangleMesh(std::move(m));
    };
    int checked = 0;
    for (double h = 8.0; h <= 30.0; h += 0.5)
        for (double hb = 8.0; hb <= 30.0; hb += 3.5) {
            TriangleMesh part = box({-44.15, -115.8, 0.0}, {44.15, 115.8, 3.0});
            MeshBoolean::cgal::plus(part, box({-25.0, -112.6, 2.9}, {25.0, -108.8, h}));
            MeshBoolean::cgal::plus(part, box({-25.0, 3.2, 2.9}, {25.0, 7.0, hb}));
            const MeshDiagnosticStats stats = its_mesh_diagnostics(part.its);
            INFO("h=" << h << " hb=" << hb);
            CHECK(!stats.has_reversed_faces);
            ++checked;
        }
    REQUIRE(checked > 100);

    // The check must still catch a really inside-out solid.
    TriangleMesh flipped = box({-44.15, -115.8, 0.0}, {44.15, 115.8, 21.0});
    its_flip_triangles(flipped.its);
    CHECK(its_mesh_diagnostics(flipped.its).has_reversed_faces);
}
