// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

// ONE immutable test source for NEUTRAL RED -> STRICT GREEN.
// Parent must wire the actual GUI export producer to this production helper BEFORE
// running RED. This file does not copy Plater orchestration or emulate a backend.
#include <catch2/catch.hpp>
#include "libslic3r/FFFStlExportSafety.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/MeshSolidValidation.hpp"

#include <algorithm>
#include <utility>

using namespace Slic3r;
namespace {
constexpr auto raw = FffStlExportMode::RawPositives;
constexpr auto boolean_mode = FffStlExportMode::BooleanPositivesMinusNegatives;

ModelVolume *add(ModelObject &object, TriangleMesh mesh, ModelVolumeType type = ModelVolumeType::MODEL_PART)
{
    return object.add_volume(std::move(mesh), type, false);
}
TriangleMesh translated_cube(float side, float x, float y = 0.f, float z = 0.f)
{
    TriangleMesh mesh = make_cube(side, side, side);
    mesh.translate(x, y, z);
    return mesh;
}
void same_mesh(const indexed_triangle_set &a, const indexed_triangle_set &b)
{
    REQUIRE(a.vertices.size() == b.vertices.size());
    REQUIRE(a.indices.size() == b.indices.size());
    REQUIRE(a.properties.size() == b.properties.size());
    for (size_t i = 0; i < a.vertices.size(); ++i) CHECK((a.vertices[i].array() == b.vertices[i].array()).all());
    for (size_t i = 0; i < a.indices.size(); ++i) CHECK((a.indices[i].array() == b.indices[i].array()).all());
    for (size_t i = 0; i < a.properties.size(); ++i) {
        CHECK(a.properties[i].type == b.properties[i].type);
        CHECK(a.properties[i].area == b.properties[i].area);
    }
}
struct Stamp {
    ObjectID id;
    ModelVolumeType type;
    const TriangleMesh *identity;
    Transform3d transform;
    indexed_triangle_set mesh;
};
std::vector<Stamp> capture(const ModelObject &object)
{
    std::vector<Stamp> result;
    for (const auto *volume : object.volumes)
        result.push_back({volume->id(), volume->type(), volume->mesh_ptr(), volume->get_matrix(), volume->mesh().its});
    return result;
}
void unchanged(const ModelObject &object, const std::vector<Stamp> &before)
{
    REQUIRE(object.volumes.size() == before.size());
    for (size_t i = 0; i < before.size(); ++i) {
        const auto *volume = object.volumes[i];
        CHECK(volume->id() == before[i].id);
        CHECK(volume->type() == before[i].type);
        CHECK(volume->mesh_ptr() == before[i].identity);
        CHECK((volume->get_matrix().matrix().array() == before[i].transform.matrix().array()).all());
        same_mesh(volume->mesh().its, before[i].mesh);
    }
}
void require_failure(const StlGeometryResult &result)
{
    CHECK(result.status == StlGeometryStatus::Failure);
    CHECK_FALSE(result.message.empty());
    CHECK(result.mesh.its.vertices.empty());
    CHECK(result.mesh.its.indices.empty());
}
void require_solid(const StlGeometryResult &result)
{
    REQUIRE(result.status == StlGeometryStatus::Success);
    REQUIRE_NOTHROW(validate_solid_mesh(result.mesh.its));
}
} // namespace

TEST_CASE("FFF strict export preserves a successful empty subtraction", "[strict_fff_export][strict_fff_export_red]")
{
    REQUIRE(MeshBoolean::openmeshcraft::is_enabled());
    for (bool containing : {false, true}) {
        INFO("Containing cutter=" << containing);
        Model model;
        REQUIRE_FALSE(model.is_need_backup());
        ModelObject *object = model.add_object();
        object->add_instance();
        add(*object, make_cube(10, 10, 10));
        add(*object, containing ? translated_cube(12, -1, -1, -1) : make_cube(10, 10, 10), ModelVolumeType::NEGATIVE_VOLUME);
        const auto before = capture(*object);
        const auto result = prepare_fff_stl_object(*object, -1, boolean_mode);
        CHECK(result.status == StlGeometryStatus::Empty);
        CHECK(result.mesh.its.vertices.empty());
        CHECK(result.mesh.its.indices.empty());
        const auto batch = finalize_stl_export_batch({result});
        CHECK_FALSE(batch.ready);
        CHECK(batch.outputs.empty());
        unchanged(*object, before);
    }
}

TEST_CASE("FFF Boolean export unions positives before subtracting cutters", "[strict_fff_export][strict_fff_export_red]")
{
    REQUIRE(MeshBoolean::openmeshcraft::is_enabled());
    SECTION("A minus A plus a disjoint B exports B, not the old accumulator A") {
        Model model;
        auto *object = model.add_object(); object->add_instance();
        add(*object, make_cube(10, 10, 10));
        add(*object, make_cube(10, 10, 10), ModelVolumeType::NEGATIVE_VOLUME);
        add(*object, translated_cube(10, 20));
        const auto before = capture(*object);
        const auto result = prepare_fff_stl_object(*object, -1, boolean_mode);
        require_solid(result);
        CHECK(its_volume(result.mesh.its) == Approx(1000.f));
        CHECK(result.mesh.bounding_box().min.x() == Approx(20.0));
        unchanged(*object, before);
    }
    SECTION("A negative listed first still cuts all positive copies") {
        Model model;
        auto *object = model.add_object(); object->add_instance();
        add(*object, translated_cube(6, 2, 2, 2), ModelVolumeType::NEGATIVE_VOLUME);
        add(*object, make_cube(10, 10, 10));
        const auto before = capture(*object);
        const auto result = prepare_fff_stl_object(*object, -1, boolean_mode);
        require_solid(result);
        CHECK(its_volume(result.mesh.its) == Approx(784.f));
        unchanged(*object, before);
    }
}

TEST_CASE("FFF raw export rejects touching or overlapping parts without repair", "[strict_fff_export][strict_fff_export_red]")
{
    for (float x : {10.f, 5.f, 0.f}) {
        INFO("Second cube X=" << x);
        Model model;
        auto *object = model.add_object(); object->add_instance();
        add(*object, make_cube(10, 10, 10));
        add(*object, translated_cube(10, x));
        const auto before = capture(*object);
        require_failure(prepare_fff_stl_object(*object, -1, raw));
        unchanged(*object, before);
    }
    SECTION("Explicit volume Boolean mode may unite overlapping operands") {
        REQUIRE(MeshBoolean::openmeshcraft::is_enabled());
        Model model;
        auto *object = model.add_object(); object->add_instance();
        add(*object, make_cube(10, 10, 10));
        add(*object, translated_cube(10, 5));
        const auto result = prepare_fff_stl_object(*object, -1, boolean_mode);
        require_solid(result);
        CHECK(its_volume(result.mesh.its) == Approx(1500.f));
    }
}

TEST_CASE("FFF operand two failure cannot publish operand one or a partial batch", "[strict_fff_export][strict_fff_export_red]")
{
    REQUIRE(MeshBoolean::openmeshcraft::is_enabled());
    Model model;
    auto *good = model.add_object(); good->add_instance();
    add(*good, make_cube(10, 10, 10));
    auto *bad = model.add_object(); bad->add_instance();
    add(*bad, translated_cube(10, 30));
    TriangleMesh open = translated_cube(2, 31, 1, 1);
    open.its.indices.pop_back(); // Valid indices: safe to characterize the real legacy path.
    add(*bad, std::move(open), ModelVolumeType::NEGATIVE_VOLUME);
    const auto before = capture(*bad);
    // NO means raw positives: the unused negative must not be silently consumed.
    require_solid(prepare_fff_stl_object(*bad, -1, raw));
    const auto first = prepare_fff_stl_object(*good, -1, raw);
    const auto second = prepare_fff_stl_object(*bad, -1, boolean_mode);
    require_solid(first);
    require_failure(second);
    require_failure(combine_stl_export_results({first, second}));
    const auto batch = finalize_stl_export_batch({first, second});
    CHECK_FALSE(batch.ready);
    CHECK(batch.outputs.empty());
    unchanged(*bad, before);
}

TEST_CASE("FFF valid cavities disconnected solids and mirrors are unchanged", "[strict_fff_export][strict_fff_export_control]")
{
    TriangleMesh cavity = make_cube(10, 10, 10);
    TriangleMesh inner = translated_cube(6, 2, 2, 2);
    its_flip_triangles(inner.its);
    cavity.merge(inner);
    Model model;
    auto *object = model.add_object();
    auto *instance = object->add_instance();
    add(*object, std::move(cavity));
    add(*object, translated_cube(2, 20));
    Transform3d mirror = Transform3d::Identity();
    mirror.linear() = Vec3d(-1.25, 0.75, 2.0).asDiagonal();
    mirror.translation() = Vec3d(4, 8, 12);
    instance->set_transformation(Geometry::Transformation(mirror));
    const auto before = capture(*object);
    const auto object_id = object->id();
    const auto instance_id = instance->id();
    const auto result = prepare_fff_stl_object(*object, 0, raw);
    require_solid(result);
    CHECK(its_volume(result.mesh.its) == Approx(float((784.0 + 8.0) * 1.25 * 0.75 * 2.0)));
    CHECK(object->id() == object_id);
    CHECK(instance->id() == instance_id);
    CHECK((instance->get_matrix().matrix().array() == mirror.matrix().array()).all());
    unchanged(*object, before);
}

TEST_CASE("FFF instance and object composition stays unfused", "[strict_fff_export][strict_fff_export_red]")
{
    Model model;
    auto *object = model.add_object();
    object->add_instance(); object->add_instance();
    add(*object, make_cube(10, 10, 10));
    require_failure(prepare_fff_stl_object(*object, -1, raw));
    require_failure(prepare_fff_stl_object(*object, -1, boolean_mode));
    require_solid(prepare_fff_stl_object(*object, 0, raw));
    require_failure(prepare_fff_stl_object(*object, 9, raw));
    object->instances[1]->set_offset(Vec3d(20, 0, 0));
    const auto disjoint = prepare_fff_stl_object(*object, -1, raw);
    require_solid(disjoint);
    const auto first = prepare_fff_stl_object(*object, 0, raw);
    require_failure(combine_stl_export_results({first, first}));
}

TEST_CASE("FFF final gate sees actual float collapse after every consumed transform", "[strict_fff_export][strict_fff_export_red]")
{
    const auto source = make_cube(1, 1, 1).its;
    const auto before = source;
    Transform3d far = Transform3d::Identity();
    far.translation().x() = 1e8;
    require_failure(prepare_stl_volume(source, {far}));
    Transform3d singular = Transform3d::Identity(); singular.linear()(2, 2) = 0.0;
    require_failure(prepare_stl_volume(source, {singular}));
    Model model;
    auto *object = model.add_object();
    auto *instance = object->add_instance();
    auto *volume = add(*object, make_cube(1, 1, 1));
    volume->set_transformation(far);
    require_failure(prepare_fff_stl_object(*object, 0, raw));
    volume->set_transformation(Transform3d::Identity());
    instance->set_offset(Vec3d(1e8, 0, 0));
    require_failure(prepare_fff_stl_object(*object, 0, raw));
    auto prepared = prepare_stl_volume(source, {});
    require_solid(prepared);
    // Simulates the caller's REAL final origin translation, not a copied export loop.
    prepared.mesh.translate(1e8f, 0, 0);
    const auto batch = finalize_stl_export_batch({prepared});
    CHECK_FALSE(batch.ready);
    CHECK(batch.outputs.empty());
    same_mesh(source, before);
}

TEST_CASE("FFF physical precision budget rejects manifold distortion and amplified earlier rounding", "[strict_fff_export][strict_fff_export_red]")
{
    const auto source = make_cube(10, 10, 10).its;
    Transform3d far = Transform3d::Identity(); far.translation().x() = 1e8;
    TriangleMesh distorted(source);
    distorted.transform(far, true);
    // Stored cube remains a closed manifold, but its X extent is now 8 rather
    // than 10 mm. Topological validity alone is therefore not the contract.
    REQUIRE_NOTHROW(validate_solid_mesh(distorted.its));
    CHECK(its_volume(distorted.its) != Approx(1000.f));
    require_failure(prepare_stl_volume(source, {far}));

    Transform3d small_shift = Transform3d::Identity(); small_shift.translation().x() = 1e-5;
    Transform3d amplify = Transform3d::Identity(); amplify.linear()(0, 0) = 1024.0;
    require_solid(prepare_stl_volume(source, {small_shift}));
    // Multiplication by 1024 is exact in float: step-local comparison would miss
    // the first step's sub-budget error amplified beyond 1e-4 mm.
    require_failure(prepare_stl_volume(source, {small_shift, amplify}));

    Model model;
    auto *object = model.add_object(); auto *instance = object->add_instance();
    auto *volume = add(*object, make_cube(10, 10, 10));
    volume->set_transformation(far);
    require_failure(prepare_fff_stl_object(*object, 0, raw));
    volume->set_transformation(Transform3d::Identity());
    instance->set_offset(Vec3d(1e8, 0, 0));
    require_failure(prepare_fff_stl_object(*object, 0, raw));
    volume->set_transformation(small_shift);
    instance->set_transformation(Geometry::Transformation(amplify));
    const auto before = capture(*object);
    require_failure(prepare_fff_stl_object(*object, 0, raw));
    require_failure(prepare_fff_stl_object(*object, 0, boolean_mode));
    unchanged(*object, before);
}

TEST_CASE("FFF legitimate Empty contributes nothing but preserves batch name mapping", "[strict_fff_export][strict_fff_export_red]")
{
    Model model;
    auto *object = model.add_object(); object->add_instance();
    add(*object, make_cube(10, 10, 10));
    add(*object, make_cube(10, 10, 10), ModelVolumeType::NEGATIVE_VOLUME);
    const auto empty = prepare_fff_stl_object(*object, -1, boolean_mode);
    const auto solid = prepare_stl_volume(translated_cube(10, 20).its, {});
    const auto combined = combine_stl_export_results({empty, solid});
    require_solid(combined);
    CHECK(its_volume(combined.mesh.its) == Approx(1000.f));
    const auto batch = finalize_stl_export_batch({empty, solid});
    REQUIRE(batch.ready);
    REQUIRE(batch.outputs.size() == 1);
    CHECK(batch.outputs.front().request_index == 1);
}

TEST_CASE("FFF cancellation exceptions propagate without publishing or changing source", "[strict_fff_export][strict_fff_export_cancel]")
{
    struct Canceled {};
    Model model;
    auto *object = model.add_object(); object->add_instance();
    add(*object, make_cube(10, 10, 10));
    const auto before = capture(*object);
    CHECK_THROWS_AS(prepare_fff_stl_object(*object, -1, boolean_mode, [] { throw Canceled{}; }), Canceled);
    int checks = 0;
    CHECK_THROWS_AS(prepare_stl_volume(object->volumes.front()->mesh().its, {}, [&] {
        if (++checks == 3) throw Canceled{};
    }), Canceled);
    unchanged(*object, before);
}
