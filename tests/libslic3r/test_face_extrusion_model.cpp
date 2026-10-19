// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>

#include "libslic3r/FaceExtrusion.hpp"
#include "libslic3r/FaceExtrusionAnnotations.hpp"
#include "libslic3r/Format/3mf.hpp"
#include "libslic3r/Format/bbs_3mf.hpp"
#include "libslic3r/MeshDiagnostics.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Preset.hpp"

#include <boost/filesystem/operations.hpp>

#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Slic3r;

namespace {
namespace fs = boost::filesystem;

// Independent closed fixture: a 20x20x2 square ring with a 10x10 through hole.
// The eight +Z cap triangles have IDs [8,16), before any volume transforms.
indexed_triangle_set ring_body()
{
    indexed_triangle_set mesh;
    mesh.vertices = {{-10.f, -10.f, 0.f}, {10.f, -10.f, 0.f}, {10.f, 10.f, 0.f}, {-10.f, 10.f, 0.f},
                     {-5.f, -5.f, 0.f}, {5.f, -5.f, 0.f}, {5.f, 5.f, 0.f}, {-5.f, 5.f, 0.f}};
    mesh.vertices.reserve(16);
    for (int i = 0; i < 8; ++i)
        mesh.vertices.emplace_back(mesh.vertices[i].x(), mesh.vertices[i].y(), 2.f);
    const std::vector<Vec3i> cap{{0, 1, 5}, {0, 5, 4}, {1, 2, 6}, {1, 6, 5},
                               {2, 3, 7}, {2, 7, 6}, {3, 0, 4}, {3, 4, 7}};
    for (const Vec3i &f : cap) mesh.indices.emplace_back(f[0], f[2], f[1]);
    for (const Vec3i &f : cap) mesh.indices.emplace_back(f[0] + 8, f[1] + 8, f[2] + 8);
    auto wall = [&mesh](int a, int b) {
        mesh.indices.emplace_back(a, b, b + 8);
        mesh.indices.emplace_back(a, b + 8, a + 8);
    };
    for (int i = 0; i < 4; ++i) {
        wall(i, (i + 1) % 4);
        wall(4 + (i + 1) % 4, 4 + i);
    }
    mesh.properties.assign(mesh.indices.size(), FaceProperty{eNormal, 0.0});
    return mesh;
}

std::map<std::string, std::string> config_values(const ModelConfig &config)
{
    std::map<std::string, std::string> values;
    for (const std::string &key : config.keys()) values.emplace(key, config.opt_serialize(key));
    return values;
}

void require_exact_mesh(const indexed_triangle_set &actual, const indexed_triangle_set &expected)
{
    REQUIRE(actual.vertices.size() == expected.vertices.size());
    REQUIRE(actual.indices.size() == expected.indices.size());
    REQUIRE(actual.properties.size() == expected.properties.size());
    for (size_t i = 0; i < actual.vertices.size(); ++i) CHECK(actual.vertices[i] == expected.vertices[i]);
    for (size_t i = 0; i < actual.indices.size(); ++i) CHECK(actual.indices[i] == expected.indices[i]);
    for (size_t i = 0; i < actual.properties.size(); ++i) {
        CHECK(actual.properties[i].type == expected.properties[i].type);
        CHECK(actual.properties[i].area == expected.properties[i].area);
    }
}

indexed_triangle_set world_mesh(const ModelVolume &volume, const ModelInstance &instance)
{
    auto mesh = volume.mesh().its;
    const Transform3d transform = instance.get_matrix() * volume.get_matrix();
    its_transform(mesh, transform, true);
    return mesh;
}

void require_closed_ring(const indexed_triangle_set &mesh)
{
    const auto stats = its_mesh_diagnostics(mesh);
    CHECK(stats.open_edges == 0);
    CHECK(stats.non_manifold_edges == 0);
    CHECK(stats.non_manifold_vertices == 0);
    CHECK(stats.same_direction_edges == 0);
    CHECK_FALSE(stats.has_reversed_faces);
    std::set<std::pair<int, int>> edges;
    for (const Vec3i &f : mesh.indices)
        for (int side = 0; side < 3; ++side)
            edges.insert(std::minmax(f[side], f[(side + 1) % 3]));
    CHECK(int(mesh.vertices.size()) - int(edges.size()) + int(mesh.indices.size()) == 0); // genus 1
    CHECK(its_volume(mesh) > 0.f);
}

// Compare triangle geometry, not only a bounding box or an unordered point cloud.
// Importers may recenter volumes, change vertex order or correct handedness.
// Return actual vertices ordered like the reference to measure cap separation.
std::vector<Vec3d> require_same_world_geometry(const indexed_triangle_set &actual,
                                             const indexed_triangle_set &expected,
                                             double tolerance_mm)
{
    REQUIRE(actual.vertices.size() == expected.vertices.size());
    REQUIRE(actual.indices.size() == expected.indices.size());
    std::vector<int> remap(actual.vertices.size(), -1);
    std::vector<bool> used(expected.vertices.size(), false);
    std::vector<Vec3d> ordered(expected.vertices.size());
    for (size_t i = 0; i < actual.vertices.size(); ++i) {
        double best_distance = std::numeric_limits<double>::max();
        size_t best = expected.vertices.size();
        for (size_t j = 0; j < expected.vertices.size(); ++j) {
            const double distance = (actual.vertices[i].cast<double>() - expected.vertices[j].cast<double>()).norm();
            if (!used[j] && distance < best_distance) { best_distance = distance; best = j; }
        }
        INFO("world vertex " << i << ", error mm = " << best_distance);
        REQUIRE(best < expected.vertices.size());
        REQUIRE(best_distance <= tolerance_mm);
        remap[i] = int(best);
        used[best] = true;
        ordered[best] = actual.vertices[i].cast<double>();
    }
    std::multiset<std::array<int, 3>> expected_faces, actual_faces;
    for (const Vec3i &f : expected.indices) {
        std::array<int, 3> face{{f[0], f[1], f[2]}};
        std::sort(face.begin(), face.end());
        expected_faces.insert(face);
    }
    for (const Vec3i &f : actual.indices) {
        for (int i = 0; i < 3; ++i) {
            REQUIRE(f[i] >= 0);
            REQUIRE(size_t(f[i]) < remap.size());
        }
        std::array<int, 3> face{{remap[f[0]], remap[f[1]], remap[f[2]]}};
        std::sort(face.begin(), face.end());
        actual_faces.insert(face);
    }
    CHECK(actual_faces == expected_faces);
    require_closed_ring(actual);
    return ordered;
}

void require_thickness(const std::vector<Vec3d> &vertices, const Vec3d &displacement, double tolerance_mm)
{
    REQUIRE(vertices.size() == 16);
    for (size_t i = 0; i < 8; ++i) {
        const Vec3d actual = vertices[i + 8] - vertices[i];
        CHECK((actual - displacement).norm() <= tolerance_mm);
        CHECK(actual.norm() == Approx(displacement.norm()).margin(tolerance_mm));
    }
}

struct ModelFixture {
    Model model;
    ModelObject *object = nullptr;
    ModelVolume *source = nullptr;
    ModelInstance *instance = nullptr;

    explicit ModelFixture(int mirror_mask)
    {
        object = model.add_object();
        object->name = "Face extrusion model fixture";
        object->config.set_key_value("extruder", new ConfigOptionInt(1));
        object->config.set_key_value("layer_height", new ConfigOptionFloat(0.24));
        source = object->add_volume(TriangleMesh(ring_body()), ModelVolumeType::MODEL_PART, false);
        source->name = "Original ring";
        source->config.set_key_value("extruder", new ConfigOptionInt(2));
        source->config.set_key_value("wall_loops", new ConfigOptionInt(4));
        source->source.input_file = "synthetic-ring.stl";
        source->source.object_idx = 3;
        source->source.volume_idx = 7;
        source->ensure_part_guid();
        Geometry::Transformation volume_transform;
        volume_transform.set_offset(Vec3d(8.0, -4.0, 25.0));
        volume_transform.set_rotation(Vec3d(0.3, 0.4, 0.2));
        volume_transform.set_scaling_factor(Vec3d(1.2, 0.7, 1.8));
        volume_transform.set_mirror(Vec3d((mirror_mask & 1) ? -1.0 : 1.0, 1.0, 1.0));
        source->set_transformation(volume_transform);
        instance = object->add_instance();
        Geometry::Transformation instance_transform;
        instance_transform.set_offset(Vec3d(75.0, 90.0, 35.0));
        instance_transform.set_rotation(Vec3d(-0.2, 0.1, 0.5));
        instance_transform.set_scaling_factor(Vec3d(0.8, 1.3, 0.9));
        instance_transform.set_mirror(Vec3d(1.0, (mirror_mask & 2) ? -1.0 : 1.0, 1.0));
        instance->set_transformation(instance_transform);
        object->invalidate_bounding_box();
    }
};

struct InsertedPart {
    ModelVolume *volume = nullptr;
    indexed_triangle_set preview_world;
    Vec3d displacement;
};

// The same storage contract as the gizmo, without wxWidgets/undo/scene calls:
// kernel in world mm -> inverse INSTANCE only -> identity MODEL_PART volume.
InsertedPart insert_part(ModelFixture &fixture, double signed_depth)
{
    ModelVolume &source = *fixture.source;
    ModelObject &object = *fixture.object;
    const TriangleMesh *original_mesh_ptr = source.mesh_ptr();
    const auto original_mesh = source.mesh().its;
    const Transform3d original_volume_matrix = source.get_matrix();
    const Transform3d original_instance_matrix = fixture.instance->get_matrix();
    const auto original_source_config = config_values(source.config);
    const auto original_object_config = config_values(object.config);
    const ObjectID original_id = source.id();
    const std::string original_name = source.name;
    const std::string original_guid = source.part_guid();

    const Transform3d to_world = fixture.instance->get_matrix() * source.get_matrix();
    const Vec3d normal_world = (to_world.linear().inverse().transpose() * Vec3d::UnitZ()).normalized();
    const Vec3d displacement = signed_depth * normal_world;
    auto input = world_mesh(source, *fixture.instance);
    std::vector<int> selected(8);
    std::iota(selected.begin(), selected.end(), 8);
    auto result = extrude_planar_face(input, selected, displacement);
    INFO(result.message);
    REQUIRE(bool(result));
    CHECK(result.boundary_edges == 8);
    InsertedPart inserted;
    inserted.preview_world = result.mesh;
    inserted.displacement = displacement;
    TriangleMesh stored(std::move(result.mesh));
    stored.transform(fixture.instance->get_matrix().inverse(), true);
    inserted.volume = object.add_volume(std::move(stored), ModelVolumeType::MODEL_PART, false);
    inserted.volume->name = "Extruded ring part";
    inserted.volume->source.is_from_builtin_objects = true;
    inserted.volume->config.set_key_value("extruder", new ConfigOptionInt(source.extruder_id()));
    object.invalidate_bounding_box();

    REQUIRE(object.volumes.size() == 2);
    CHECK(object.volumes.front() == &source);
    CHECK(inserted.volume != &source);
    CHECK(inserted.volume->id() != source.id());
    CHECK(inserted.volume->type() == ModelVolumeType::MODEL_PART);
    CHECK(inserted.volume->get_matrix().matrix().isIdentity(1e-12));
    CHECK(inserted.volume->source.input_file.empty());
    CHECK(inserted.volume->source.object_idx == -1);
    CHECK(inserted.volume->source.volume_idx == -1);
    CHECK(inserted.volume->config.keys().size() == 1);
    CHECK(inserted.volume->config.opt_int("extruder") == 2);
    CHECK_FALSE(inserted.volume->config.has("wall_loops"));
    CHECK(source.mesh_ptr() == original_mesh_ptr);
    CHECK(source.id() == original_id);
    CHECK(source.name == original_name);
    CHECK(source.part_guid() == original_guid);
    CHECK(source.source.input_file == "synthetic-ring.stl");
    CHECK(source.source.object_idx == 3);
    CHECK(source.source.volume_idx == 7);
    CHECK(source.get_matrix().matrix() == original_volume_matrix.matrix());
    CHECK(fixture.instance->get_matrix().matrix() == original_instance_matrix.matrix());
    CHECK(config_values(source.config) == original_source_config);
    CHECK(config_values(object.config) == original_object_config);
    require_exact_mesh(source.mesh().its, original_mesh);
    const auto actual_world = world_mesh(*inserted.volume, *fixture.instance);
    require_thickness(require_same_world_geometry(actual_world, inserted.preview_world, 1e-4), displacement, 1e-4);
    return inserted;
}

// Never touch TEST_DATA_DIR. A directory is acquired atomically before assigning
// any filename. Cleanup accepts only known children, checks canonical absolute
// paths and rejects symlinks; no recursive deletion is used.
class Temporary3mfDirectory {
public:
    Temporary3mfDirectory()
    {
        const fs::path parent = fs::canonical(fs::temp_directory_path());
        for (int attempt = 0; attempt < 16; ++attempt) {
            const fs::path candidate = parent / fs::unique_path("bambu-face-extrusion-%%%%-%%%%-%%%%");
            if (fs::create_directory(candidate)) {
                m_directory = fs::canonical(candidate);
                m_file = (m_directory / "roundtrip.3mf").string();
                return;
            }
        }
        throw std::runtime_error("Cannot create a unique face extrusion test directory");
    }
    Temporary3mfDirectory(const Temporary3mfDirectory &) = delete;
    Temporary3mfDirectory &operator=(const Temporary3mfDirectory &) = delete;
    ~Temporary3mfDirectory()
    {
        try {
            if (!verified_directory()) return;
            for (const char *name : {"roundtrip.3mf", "roundtrip.3mf.tmp", "_temp_1.config", "_temp_2.config", "_temp_3.config"}) {
                const fs::path file = m_directory / name;
                if (fs::is_regular_file(fs::symlink_status(file)) && fs::canonical(file) == file) {
                    boost::system::error_code ec;
                    fs::remove(file, ec);
                }
            }
            const fs::path auxiliary = m_directory / "Auxiliaries";
            if (fs::is_directory(fs::symlink_status(auxiliary)) && fs::canonical(auxiliary) == auxiliary && fs::is_empty(auxiliary)) {
                boost::system::error_code ec;
                fs::remove(auxiliary, ec);
            }
            if (verified_directory() && fs::is_empty(m_directory)) {
                boost::system::error_code ec;
                fs::remove(m_directory, ec);
            }
        } catch (...) {
            // A failed safety check leaves artifacts for diagnosis, never removes
            // an unexpected path. Tests also assert that normal cleanup succeeded.
        }
    }
    const std::string &filename() const { return m_file; }
    const fs::path &directory() const { return m_directory; }
private:
    bool verified_directory() const
    {
        return m_directory.is_absolute() && fs::is_directory(fs::symlink_status(m_directory)) &&
               fs::canonical(m_directory) == m_directory;
    }
    fs::path m_directory;
    std::string m_file;
};

// BBS import/export asks Model for a backup path even in a non-GUI test.
// Constrain it to our existing private directory, then detach before Model's
// destructor can enqueue asynchronous removal through the GUI backup manager.
class ScopedModelScratch {
public:
    ScopedModelScratch(Model &model, const fs::path &directory) : m_model(model)
    {
        m_model.set_backup_path(directory.string());
    }
    ~ScopedModelScratch() { m_model.set_backup_path("detach"); }
    ScopedModelScratch(const ScopedModelScratch &) = delete;
    ScopedModelScratch &operator=(const ScopedModelScratch &) = delete;
private:
    Model &m_model;
};

struct ImportedBbsData {
    PlateDataPtrs plates;
    std::vector<Preset *> presets;
    ~ImportedBbsData()
    {
        release_PlateData_list(plates);
        for (Preset *preset : presets) delete preset;
    }
};

const ModelVolume &volume_named(const ModelObject &object, const std::string &name)
{
    const auto found = std::find_if(object.volumes.begin(), object.volumes.end(), [&name](const ModelVolume *volume) {
        return volume->name == name;
    });
    REQUIRE(found != object.volumes.end());
    return **found;
}

void require_roundtrip_model(const Model &loaded, const ModelFixture &fixture, const InsertedPart &inserted,
                             const indexed_triangle_set &source_world, const indexed_triangle_set &part_world,
                             bool bbs_format)
{
    REQUIRE(loaded.objects.size() == 1);
    const ModelObject &object = *loaded.objects.front();
    REQUIRE(object.instances.size() == 1);
    REQUIRE(object.volumes.size() == 2);
    CHECK(object.name == fixture.object->name);
    const ModelVolume &source = volume_named(object, fixture.source->name);
    const ModelVolume &part = volume_named(object, inserted.volume->name);
    CHECK(&source != &part);
    CHECK(source.type() == ModelVolumeType::MODEL_PART);
    CHECK(part.type() == ModelVolumeType::MODEL_PART);
    CHECK(source.config.has("extruder"));
    CHECK(part.config.has("extruder"));
    CHECK(source.config.opt_int("extruder") == 2);
    CHECK(part.config.opt_int("extruder") == 2);
    CHECK(source.source.input_file == "synthetic-ring.stl");
    // Both importers fill previously empty provenance with the loaded archive.
    CHECK(fs::path(part.source.input_file).extension() == ".3mf");
    CHECK(part.source.input_file != source.source.input_file);
    // The legacy Prusa-compatible loader intentionally filters other config
    // keys. Native BBS is the GUI's format and must also preserve these overrides.
    if (bbs_format) {
        REQUIRE(source.config.has("wall_loops"));
        CHECK(source.config.opt_int("wall_loops") == 4);
        REQUIRE(object.config.has("layer_height"));
        CHECK(object.config.opt_float("layer_height") == Approx(0.24));
        CHECK(source.part_guid() == fixture.source->part_guid());
        CHECK(part.part_guid() == inserted.volume->part_guid());
    }
    require_same_world_geometry(world_mesh(source, *object.instances.front()), source_world, 1e-3);
    const auto ordered = require_same_world_geometry(world_mesh(part, *object.instances.front()), part_world, 1e-3);
    require_thickness(ordered, inserted.displacement, 2e-3);
    // No assertion on loaded volume identity: import may recenter or bake a
    // transform. The physical world geometry and separate parts are the contract.
}

} // namespace

TEST_CASE("Face extrusion inserts an identity Model part without changing the source", "[face_extrusion_model]")
{
    for (int mirrors = 0; mirrors < 4; ++mirrors)
        for (double signed_depth : {3.0, -3.0}) {
            INFO("mirror mask " << mirrors << ", signed depth " << signed_depth);
            ModelFixture fixture(mirrors);
            const auto inserted = insert_part(fixture, signed_depth);
            require_closed_ring(inserted.volume->mesh().its);
        }
}

TEST_CASE("Face extrusion parts survive a legacy Prusa-compatible 3MF roundtrip", "[face_extrusion_model][3mf]")
{
    for (double signed_depth : {3.0, -3.0}) {
        fs::path temporary_directory;
        {
            Temporary3mfDirectory file;
            temporary_directory = file.directory();
            INFO("3MF fixture " << file.filename());
            ModelFixture fixture(1); // mirrored volume, rotated/nonuniform instance
            const auto inserted = insert_part(fixture, signed_depth);
            const auto source_world = world_mesh(*fixture.source, *fixture.instance);
            const auto part_world = world_mesh(*inserted.volume, *fixture.instance);
            REQUIRE(store_3mf(file.filename().c_str(), &fixture.model, nullptr, false));
            REQUIRE(fs::is_regular_file(file.filename()));
            REQUIRE(fs::file_size(file.filename()) > 0);
            Model loaded;
            DynamicPrintConfig config;
            ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Disable};
            REQUIRE(load_3mf(file.filename().c_str(), config, substitutions, &loaded, false));
            require_roundtrip_model(loaded, fixture, inserted, source_world, part_world, false);
        }
        CHECK_FALSE(fs::exists(temporary_directory));
    }
}

TEST_CASE("Face extrusion parts survive the native Bambu 3MF model roundtrip", "[face_extrusion_model][bbs_3mf]")
{
    for (double signed_depth : {3.0, -3.0}) {
        fs::path temporary_directory;
        {
            Temporary3mfDirectory file;
            temporary_directory = file.directory();
            INFO("3MF fixture " << file.filename());
            ModelFixture fixture(2); // mirrored instance, rotated/nonuniform volume
            Model loaded;
            ScopedModelScratch source_scratch(fixture.model, file.directory());
            ScopedModelScratch loaded_scratch(loaded, file.directory());
            const auto inserted = insert_part(fixture, signed_depth);
            const auto source_world = world_mesh(*fixture.source, *fixture.instance);
            const auto part_world = world_mesh(*inserted.volume, *fixture.instance);
            // full_print_config() clones static defaults, whose generic enum
            // vectors may lack keys_map. Apply through the dynamic option
            // factory instead: it binds enum definitions while copying the same
            // values, as required by the native JSON project-settings exporter.
            DynamicPrintConfig config;
            config.apply(DynamicPrintConfig::full_print_config());
            for (const std::string &key : config.keys()) {
                INFO("Serializable default enum: " << key);
                const ConfigOption *option = config.option(key);
                if (const auto *value = dynamic_cast<const ConfigOptionEnumsGeneric *>(option))
                    REQUIRE(value->keys_map != nullptr);
                if (const auto *value = dynamic_cast<const ConfigOptionEnumsGenericNullable *>(option))
                    REQUIRE(value->keys_map != nullptr);
                if (const auto *value = dynamic_cast<const ConfigOptionEnumGeneric *>(option))
                    REQUIRE(value->keys_map != nullptr);
            }
            // BBS clamps part extruders to filament_settings_id.size() on load.
            // Represent the two real filament slots used by this model fixture.
            config.set_key_value("filament_settings_id", new ConfigOptionStrings({"Test filament 1", "Test filament 2"}));
            config.set_key_value("filament_colour", new ConfigOptionStrings({"#FFFFFF", "#000000"}));
            config.set_key_value("filament_type", new ConfigOptionStrings({"PLA", "PLA"}));
            PlateData plate;
            plate.plate_index = 0;
            plate.plate_name = "Face extrusion test plate";
            plate.objects_and_instances.emplace_back(0, 0);
            StoreParams params;
            params.path = file.filename().c_str();
            params.model = &fixture.model;
            params.config = &config;
            params.plate_data_list.push_back(&plate);
            params.strategy = SaveStrategy::Zip64 | SaveStrategy::Silence | SaveStrategy::SkipAuxiliary;
            REQUIRE(store_bbs_3mf(params));
            REQUIRE(fs::is_regular_file(file.filename()));
            REQUIRE(fs::file_size(file.filename()) > 0);
            DynamicPrintConfig loaded_config;
            ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Disable};
            ImportedBbsData imported;
            bool is_bambu = false;
            Semver version;
            const LoadStrategy strategy = LoadStrategy::LoadModel | LoadStrategy::LoadConfig | LoadStrategy::Silence;
            REQUIRE(load_bbs_3mf(file.filename().c_str(), &loaded_config, &substitutions, &loaded,
                                 &imported.plates, &imported.presets, &is_bambu, &version, nullptr, strategy));
            CHECK(is_bambu);
            CHECK(imported.plates.size() == 1);
            CHECK(imported.presets.empty());
            require_roundtrip_model(loaded, fixture, inserted, source_world, part_world, true);
        }
        CHECK_FALSE(fs::exists(temporary_directory));
    }
}

namespace {

bool openmeshcraft_wrapper_regressions_enabled()
{
    if (MeshBoolean::openmeshcraft::is_enabled())
        return true;
    // Catch2 v2 has no native SKIP. Make non-execution visible instead of
    // pretending to validate the legacy MCUT backend with this wrapper suite.
    WARN("SKIPPED: OpenMeshCraft is disabled; wrapper validation/fallback assertions "
         "were not executed and do not provide MCUT coverage.");
    return false;
}

indexed_triangle_set boolean_ring(float bottom = 0.f, float top = 4.f, bool alternate = false)
{
    auto mesh = ring_body();
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
        mesh.vertices[i].x() *= i % 8 < 4 ? 3.f : 4.f;
        mesh.vertices[i].y() *= 2.f;
        mesh.vertices[i].z() = i < 8 ? bottom : top;
    }
    if (alternate) {
        for (int i = 0; i < 4; ++i) {
            const int a = i, b = (i + 1) % 4, c = b + 4, d = a + 4;
            mesh.indices[2 * i] = Vec3i(a, d, b);
            mesh.indices[2 * i + 1] = Vec3i(b, d, c);
            mesh.indices[8 + 2 * i] = Vec3i(a + 8, b + 8, d + 8);
            mesh.indices[9 + 2 * i] = Vec3i(b + 8, c + 8, d + 8);
        }
    }
    return mesh;
}

struct SolidBooleanResult {
    std::vector<TriangleMesh> meshes;
    std::vector<float> progress;
    int failures = 0;
};

SolidBooleanResult run_solid_boolean(const TriangleMesh &a, const TriangleMesh &b, const std::string &operation)
{
    const auto a_before = a.its;
    const auto b_before = b.its;
    SolidBooleanResult result;
    MeshBoolean::mcut::make_boolean(a, b, result.meshes, operation, [] { return false; },
        [&](float p) { result.progress.push_back(p); }, [&] { ++result.failures; });
    require_exact_mesh(a.its, a_before);
    require_exact_mesh(b.its, b_before);
    return result;
}

void require_boolean_result(const SolidBooleanResult &result, double volume, size_t shells = 1, long euler = 0)
{
    REQUIRE(result.failures == 0);
    REQUIRE_FALSE(result.progress.empty());
    CHECK(std::is_sorted(result.progress.begin(), result.progress.end()));
    CHECK(result.progress.back() == 100.f);
    if (volume == 0.0) {
        CHECK(result.meshes.empty());
        return;
    }
    REQUIRE(result.meshes.size() == 1);
    const auto &mesh = result.meshes.front().its;
    const auto stats = its_mesh_diagnostics(mesh);
    CHECK(stats.open_edges == 0);
    CHECK(stats.non_manifold_edges == 0);
    CHECK(stats.non_manifold_vertices == 0);
    CHECK(stats.same_direction_edges == 0);
    CHECK_FALSE(stats.has_reversed_faces);
    CHECK(its_volume(mesh) == Approx(volume).epsilon(2e-5).margin(1e-3));
    // TriangleMesh::split() flips cavity shells and discards small components.
    // its_split() is the non-repairing topology inspection needed here.
    CHECK(its_split(mesh).size() == shells);
    std::set<int> vertices;
    std::set<std::pair<int, int>> edges;
    for (const Vec3i &face : mesh.indices)
        for (int i = 0; i < 3; ++i) {
            vertices.insert(face[i]);
            edges.emplace(std::min(face[i], face[(i + 1) % 3]), std::max(face[i], face[(i + 1) % 3]));
        }
    CHECK(long(vertices.size()) - long(edges.size()) + long(mesh.indices.size()) == euler);
    auto cgal_mesh = MeshBoolean::cgal::triangle_mesh_to_cgal(mesh);
    CHECK_FALSE(MeshBoolean::cgal::does_self_intersect(*cgal_mesh));
    CHECK(MeshBoolean::cgal::does_bound_a_volume(*cgal_mesh));
}

std::vector<int> ring_top_faces()
{
    std::vector<int> faces(8);
    std::iota(faces.begin(), faces.end(), 8);
    return faces;
}

} // namespace

TEST_CASE("Solid boolean validates coincident primitive and extrusion unions", "[face_extrusion_model][mesh_boolean]")
{
    if (!openmeshcraft_wrapper_regressions_enabled()) return;
    const TriangleMesh source{boolean_ring()};
    auto extrusion = extrude_planar_face(source.its, ring_top_faces(), Vec3d(0, 0, 5));
    REQUIRE(extrusion);
    const TriangleMesh part{std::move(extrusion.mesh)};
    require_boolean_result(run_solid_boolean(source, part, "UNION"), 14400.0);
    for (bool alternate : {false, true}) {
        INFO("Independent annular prism, alternate diagonal=" << alternate);
        const TriangleMesh primitive{boolean_ring(4, 9, alternate)};
        require_boolean_result(run_solid_boolean(source, primitive, "UNION"), 14400.0);
    }
    TriangleMesh a = make_cube(60, 40, 4);
    TriangleMesh b = make_cube(60, 40, 5);
    b.translate(0, 0, 4);
    require_boolean_result(run_solid_boolean(a, b, "UNION"), 21600.0, 1, 2);
}

TEST_CASE("Solid boolean preserves signed extrusion set operations and holes", "[face_extrusion_model][mesh_boolean]")
{
    if (!openmeshcraft_wrapper_regressions_enabled()) return;
    const TriangleMesh source{boolean_ring()};
    for (double depth : {5.0, -2.0, -5.0}) {
        INFO("Signed depth=" << depth);
        auto extrusion = extrude_planar_face(source.its, ring_top_faces(), Vec3d(0, 0, depth));
        REQUIRE(extrusion);
        const TriangleMesh part{std::move(extrusion.mesh)};
        require_boolean_result(run_solid_boolean(source, part, "UNION"), depth > 0 ? 14400.0 : depth == -2 ? 6400.0 : 8000.0);
        require_boolean_result(run_solid_boolean(source, part, "A_NOT_B"), depth > 0 ? 6400.0 : depth == -2 ? 3200.0 : 0.0);
        require_boolean_result(run_solid_boolean(source, part, "B_NOT_A"), depth > 0 ? 8000.0 : depth == -2 ? 0.0 : 1600.0);
        require_boolean_result(run_solid_boolean(source, part, "INTERSECTION"), depth > 0 ? 0.0 : depth == -2 ? 3200.0 : 6400.0);
    }
}

TEST_CASE("Solid boolean validates transformed reflected Model extrusion unions", "[face_extrusion_model][mesh_boolean]")
{
    if (!openmeshcraft_wrapper_regressions_enabled()) return;
    for (int mirror_mask = 0; mirror_mask < 4; ++mirror_mask) {
        INFO("Mirror mask=" << mirror_mask);
        ModelFixture fixture(mirror_mask);
        const auto inserted = insert_part(fixture, 3.0);
        const auto source_world = world_mesh(*fixture.source, *fixture.instance);
        const auto part_world = world_mesh(*inserted.volume, *fixture.instance);
        double area = 0.0;
        for (int id : ring_top_faces()) {
            const Vec3i &face = source_world.indices[id];
            const Vec3d a = source_world.vertices[face[0]].cast<double>();
            const Vec3d b = source_world.vertices[face[1]].cast<double>();
            const Vec3d c = source_world.vertices[face[2]].cast<double>();
            area += (b - a).cross(c - a).norm() / 2.0;
        }
        const double expected = its_volume(source_world) + area * 3.0;
        require_boolean_result(run_solid_boolean(TriangleMesh{source_world}, TriangleMesh{part_world}, "UNION"), expected);
    }
}

TEST_CASE("Solid boolean subtraction preserves negatively oriented cavity shells", "[face_extrusion_model][mesh_boolean]")
{
    if (!openmeshcraft_wrapper_regressions_enabled()) return;
    TriangleMesh outer = make_cube(10, 10, 10);
    TriangleMesh inner = make_cube(6, 6, 6);
    inner.translate(2, 2, 2);
    const auto result = run_solid_boolean(outer, inner, "A_NOT_B");
    require_boolean_result(result, 784.0, 2, 4);
    std::vector<float> volumes;
    for (const auto &shell : its_split(result.meshes.front().its)) volumes.push_back(its_volume(shell));
    std::sort(volumes.begin(), volumes.end());
    REQUIRE(volumes.size() == 2);
    CHECK(volumes[0] == Approx(-216.0).margin(1e-3));
    CHECK(volumes[1] == Approx(1000.0).margin(1e-3));
}

TEST_CASE("Solid boolean distinguishes valid empty sets from failures", "[face_extrusion_model][mesh_boolean]")
{
    if (!openmeshcraft_wrapper_regressions_enabled()) return;
    const TriangleMesh empty;
    const TriangleMesh cube = make_cube(2, 2, 2);
    TriangleMesh disjoint = cube;
    disjoint.translate(4, 0, 0);
    require_boolean_result(run_solid_boolean(cube, cube, "A_NOT_B"), 0.0);
    require_boolean_result(run_solid_boolean(cube, cube, "B_NOT_A"), 0.0);
    require_boolean_result(run_solid_boolean(cube, disjoint, "INTERSECTION"), 0.0);
    require_boolean_result(run_solid_boolean(cube, disjoint, "UNION"), 16.0, 2, 4);
    for (const std::string operation : {"UNION", "A_NOT_B", "B_NOT_A", "INTERSECTION"}) {
        INFO(operation);
        require_boolean_result(run_solid_boolean(empty, empty, operation), 0.0);
        require_boolean_result(run_solid_boolean(cube, empty, operation), operation == "UNION" || operation == "A_NOT_B" ? 8.0 : 0.0, 1, 2);
        require_boolean_result(run_solid_boolean(empty, cube, operation), operation == "UNION" || operation == "B_NOT_A" ? 8.0 : 0.0, 1, 2);
    }
}

TEST_CASE("Solid boolean rejects invalid operands without touching existing output", "[face_extrusion_model][mesh_boolean]")
{
    if (!openmeshcraft_wrapper_regressions_enabled()) return;
    const TriangleMesh good = make_cube(10, 10, 10);
    const auto sentinel = good.its;
    using Corrupt = std::function<void(TriangleMesh &)>;
    const std::vector<std::pair<std::string, Corrupt>> corruptions{
        {"open", [](TriangleMesh &m) { m.its.indices.pop_back(); }},
        {"negative index", [](TriangleMesh &m) { m.its.indices[0][0] = -1; }},
        {"out of range", [](TriangleMesh &m) { m.its.indices[0][0] = int(m.its.vertices.size()); }},
        {"NaN", [](TriangleMesh &m) { m.its.vertices[0].x() = std::numeric_limits<float>::quiet_NaN(); }},
        {"infinity", [](TriangleMesh &m) { m.its.vertices[0].x() = std::numeric_limits<float>::infinity(); }},
        {"degenerate", [](TriangleMesh &m) { m.its.indices[0][0] = m.its.indices[0][1]; }},
        {"zero area", [](TriangleMesh &m) { m.its.vertices[m.its.indices[0][0]] = m.its.vertices[m.its.indices[0][1]]; }},
        {"unrepresentable float volume", [](TriangleMesh &m) {
            m.its.vertices.insert(m.its.vertices.begin(), Vec3f(std::numeric_limits<float>::max(), 0, 0));
            for (Vec3i &face : m.its.indices) face.array() += 1;
        }},
        {"duplicate face", [](TriangleMesh &m) { m.its.indices.push_back(m.its.indices[0]); }},
        {"one reversed face", [](TriangleMesh &m) { std::swap(m.its.indices[0][0], m.its.indices[0][1]); }},
        {"inward solid", [](TriangleMesh &m) { m.flip_triangles(); }},
        {"self-intersection", [](TriangleMesh &m) { TriangleMesh overlap = m; overlap.translate(5, 5, 5); m.merge(overlap); }}
    };
    for (const auto &corrupt : corruptions) {
        INFO(corrupt.first);
        TriangleMesh bad = good;
        corrupt.second(bad); // Mutate only after safe construction; cached stats must not be trusted.
        std::vector<TriangleMesh> output{good};
        int failures = 0;
        float progress = 0;
        MeshBoolean::mcut::make_boolean(bad, good, output, "UNION", [] { return false; },
            [&](float p) { progress = p; }, [&] { ++failures; });
        CHECK(failures == 1);
        CHECK(progress < 100.f);
        REQUIRE(output.size() == 1);
        require_exact_mesh(output.front().its, sentinel);
        require_exact_mesh(good.its, sentinel);
    }
    std::vector<TriangleMesh> output{good};
    int failures = 0;
    MeshBoolean::mcut::make_boolean(good, good, output, "invalid-operation", nullptr, nullptr, [&] { ++failures; });
    CHECK(failures == 1);
    REQUIRE(output.size() == 1);
    require_exact_mesh(output.front().its, sentinel);
}

TEST_CASE("Solid boolean cancellation and callback errors publish no partial result", "[face_extrusion_model][mesh_boolean]")
{
    if (!openmeshcraft_wrapper_regressions_enabled()) return;
    const TriangleMesh source{boolean_ring()};
    const TriangleMesh part{boolean_ring(4, 9)};
    const auto sentinel = source.its;
    for (float stop_at : {0.f, 60.f, 90.f, 100.f}) {
        INFO("Cancel stage=" << stop_at);
        std::vector<TriangleMesh> output{source};
        bool canceled = stop_at == 0.f;
        int failures = 0;
        MeshBoolean::mcut::make_boolean(source, part, output, "UNION", [&] { return canceled; },
            [&](float p) { if (p >= stop_at) canceled = true; }, [&] { ++failures; });
        CHECK(failures == 0);
        REQUIRE(output.size() == 1);
        require_exact_mesh(output.front().its, sentinel);
        require_exact_mesh(source.its, sentinel);
    }
    for (float throw_at : {90.f, 100.f}) {
        INFO("Throwing progress stage=" << throw_at);
        std::vector<TriangleMesh> output{source};
        int failures = 0;
        MeshBoolean::mcut::make_boolean(source, part, output, "UNION", [] { return false; },
            [&](float p) { if (p >= throw_at) throw std::runtime_error("test callback failure"); }, [&] { ++failures; });
        CHECK(failures == 1);
        REQUIRE(output.size() == 1);
        require_exact_mesh(output.front().its, sentinel);
    }
}

TEST_CASE("Solid boolean publishes nothing when valid solids have a non-manifold union", "[face_extrusion_model][mesh_boolean]")
{
    if (!openmeshcraft_wrapper_regressions_enabled()) return;
    const TriangleMesh source = make_cube(10, 10, 10);
    TriangleMesh other = source;
    other.translate(10, 10, 0); // Only one vertical edge is shared, not a face.
    const auto source_before = source.its;
    const auto other_before = other.its;
    std::vector<TriangleMesh> output{source};
    int failures = 0;
    float progress = 0;
    MeshBoolean::mcut::make_boolean(source, other, output, "UNION", [] { return false; },
        [&](float p) { progress = p; }, [&] { ++failures; });
    CHECK(failures == 1);
    CHECK(progress < 100.f);
    REQUIRE(output.size() == 1);
    require_exact_mesh(output.front().its, source_before);
    require_exact_mesh(source.its, source_before);
    require_exact_mesh(other.its, other_before);
}

TEST_CASE("Multi-face extrusion preserves independent holed parts through native BBS 3MF", "[face_extrusion_model][multiface]")
{
    // Model-layer batch contract only: GUI face picking and the single undo
    // snapshot are deliberately not simulated by this integration test.
    for (int mirror_mask : {1, 2})
        for (double signed_depth : {1.25, -1.25}) {
            INFO("Batch mirror mask=" << mirror_mask << ", signed world depth=" << signed_depth);
            fs::path temporary_directory;
            {
                Temporary3mfDirectory file;
                temporary_directory = file.directory();
                INFO("Multi-face 3MF fixture " << file.filename());
                ModelFixture fixture(mirror_mask);
                Model loaded;
                ScopedModelScratch source_scratch(fixture.model, file.directory());
                ScopedModelScratch loaded_scratch(loaded, file.directory());

                // Two disconnected annular bodies in ONE source volume. Their
                // +Z patches are parallel but not coplanar: local z=2 and z=4.
                auto combined = ring_body();
                auto second = ring_body();
                for (Vec3f &vertex : second.vertices) {
                    vertex.x() += 35.f;
                    vertex.z() *= 2.f;
                }
                const int second_face_offset = int(combined.indices.size());
                its_merge(combined, second);
                combined.properties.assign(combined.indices.size(), FaceProperty{eNormal, 0.0});
                fixture.source->set_mesh(TriangleMesh{std::move(combined)});
                fixture.source->calculate_convex_hull();
                fixture.source->name = "Original two-ring source";
                fixture.source->source.input_file = "synthetic-two-rings.stl";
                // Cover both a volume override and a filament inherited from
                // its parent object, matching the GUI's effective extruder rule.
                if (mirror_mask == 2) fixture.source->config.erase("extruder");
                fixture.object->invalidate_bounding_box();

                ModelVolume &source = *fixture.source;
                const TriangleMesh *const original_mesh_ptr = source.mesh_ptr();
                const auto original_mesh = source.mesh().its;
                const auto original_source_config = config_values(source.config);
                const auto original_object_config = config_values(fixture.object->config);
                const Transform3d original_volume_matrix = source.get_matrix();
                const Transform3d original_instance_matrix = fixture.instance->get_matrix();
                const ObjectID original_id = source.id();
                const std::string original_name = source.name;
                const std::string original_guid = source.part_guid();
                const auto original_material_id = source.material_id();
                const int extruder = source.extruder_id();
                REQUIRE(extruder == (mirror_mask == 1 ? 2 : 1));
                const auto source_world = world_mesh(source, *fixture.instance);
                const auto source_shells = its_split(source_world);
                REQUIRE(source_shells.size() == 2);
                for (const auto &shell : source_shells) require_closed_ring(shell);

                auto require_source_unchanged = [&] {
                    REQUIRE_FALSE(fixture.object->volumes.empty());
                    CHECK(fixture.object->volumes.front() == &source);
                    CHECK(source.mesh_ptr() == original_mesh_ptr);
                    CHECK(source.id() == original_id);
                    CHECK(source.name == original_name);
                    CHECK(source.part_guid() == original_guid);
                    CHECK(source.material_id() == original_material_id);
                    CHECK(source.source.input_file == "synthetic-two-rings.stl");
                    CHECK(source.source.object_idx == 3);
                    CHECK(source.source.volume_idx == 7);
                    CHECK(source.get_matrix().matrix() == original_volume_matrix.matrix());
                    REQUIRE(fixture.object->instances.size() == 1);
                    CHECK(fixture.object->instances.front() == fixture.instance);
                    CHECK(fixture.instance->get_matrix().matrix() == original_instance_matrix.matrix());
                    CHECK(config_values(source.config) == original_source_config);
                    CHECK(config_values(fixture.object->config) == original_object_config);
                    require_exact_mesh(source.mesh().its, original_mesh);
                };

                const Transform3d to_world = original_instance_matrix * original_volume_matrix;
                const Vec3d normal = (to_world.linear().inverse().transpose() * Vec3d::UnitZ()).normalized();
                const Vec3d displacement = signed_depth * normal;
                std::array<indexed_triangle_set, 2> previews;
                std::array<double, 2> plane_height{};
                std::vector<TriangleMesh> prepared;
                prepared.reserve(2);
                for (size_t patch = 0; patch < 2; ++patch) {
                    std::vector<int> triangles(8);
                    std::iota(triangles.begin(), triangles.end(), 8 + (patch ? second_face_offset : 0));
                    const Vec3i &face = source_world.indices[triangles.front()];
                    const Vec3d a = source_world.vertices[face[0]].cast<double>();
                    const Vec3d b = source_world.vertices[face[1]].cast<double>();
                    const Vec3d c = source_world.vertices[face[2]].cast<double>();
                    CHECK(((b - a).cross(c - a).normalized() - normal).norm() < 1e-5);
                    plane_height[patch] = normal.dot(a);
                    auto result = extrude_planar_face(source_world, triangles, displacement);
                    INFO(result.message);
                    REQUIRE(bool(result));
                    CHECK(result.boundary_edges == 8);
                    REQUIRE(result.mesh.vertices.size() == 16);
                    std::set<std::array<float, 3>> expected_base, actual_base;
                    for (int id : triangles)
                        for (int corner = 0; corner < 3; ++corner) {
                            const Vec3f &vertex = source_world.vertices[source_world.indices[id][corner]];
                            expected_base.insert({{vertex.x(), vertex.y(), vertex.z()}});
                        }
                    for (size_t i = 0; i < 8; ++i) {
                        const Vec3f &vertex = result.mesh.vertices[i];
                        actual_base.insert({{vertex.x(), vertex.y(), vertex.z()}});
                    }
                    REQUIRE(expected_base.size() == 8);
                    CHECK(actual_base == expected_base); // Independent XY/height oracle, not only the preview.
                    previews[patch] = result.mesh;
                    TriangleMesh stored(std::move(result.mesh));
                    stored.transform(original_instance_matrix.inverse(), true);
                    prepared.push_back(std::move(stored));
                }
                CHECK(plane_height[1] - plane_height[0] > 0.1);
                REQUIRE(prepared.size() == 2);
                REQUIRE(fixture.object->volumes.size() == 1); // Entire batch prepared before any insertion.
                require_source_unchanged();

                std::array<ModelVolume *, 2> parts{};
                for (size_t patch = 0; patch < parts.size(); ++patch) {
                    auto *part = fixture.object->add_volume(std::move(prepared[patch]), ModelVolumeType::MODEL_PART, false);
                    parts[patch] = part;
                    part->name = "Extruded parallel patch " + std::to_string(patch + 1);
                    part->source.is_from_builtin_objects = true;
                    part->config.set_key_value("extruder", new ConfigOptionInt(extruder));
                    CHECK(part->type() == ModelVolumeType::MODEL_PART);
                    CHECK(part->get_matrix().matrix().isIdentity(1e-12));
                    CHECK(part->id() != source.id());
                    CHECK(part->source.input_file.empty());
                    CHECK(part->source.object_idx == -1);
                    CHECK(part->source.volume_idx == -1);
                    REQUIRE(part->config.keys().size() == 1);
                    CHECK(part->config.opt_int("extruder") == extruder);
                    CHECK_FALSE(part->config.has("wall_loops"));
                    require_source_unchanged();
                }
                fixture.object->invalidate_bounding_box();
                REQUIRE(fixture.object->volumes.size() == 3);
                CHECK(parts[0] != parts[1]);
                CHECK(parts[0]->id() != parts[1]->id());

                auto require_part_geometry = [&](const ModelVolume &part, const ModelInstance &instance,
                                                 size_t patch, double tolerance) {
                    const auto actual_world = world_mesh(part, instance);
                    REQUIRE(its_split(actual_world).size() == 1);
                    const auto ordered = require_same_world_geometry(actual_world, previews[patch], tolerance);
                    require_thickness(ordered, displacement, tolerance * 2.0);
                    for (size_t i = 0; i < 8; ++i)
                        CHECK(normal.dot(ordered[i]) == Approx(plane_height[patch]).margin(tolerance));
                };
                for (size_t patch = 0; patch < parts.size(); ++patch)
                    require_part_geometry(*parts[patch], *fixture.instance, patch, 1e-4);

                // Use the real dynamic-option factory, not direct cloning of
                // static enum defaults (which may lack their serialization map).
                DynamicPrintConfig config;
                config.apply(DynamicPrintConfig::full_print_config());
                config.set_key_value("filament_settings_id", new ConfigOptionStrings({"Test filament 1", "Test filament 2"}));
                config.set_key_value("filament_colour", new ConfigOptionStrings({"#FFFFFF", "#000000"}));
                config.set_key_value("filament_type", new ConfigOptionStrings({"PLA", "PLA"}));
                PlateData plate;
                plate.plate_index = 0;
                plate.plate_name = "Multi-face extrusion plate";
                plate.objects_and_instances.emplace_back(0, 0);
                StoreParams params;
                params.path = file.filename().c_str();
                params.model = &fixture.model;
                params.config = &config;
                params.plate_data_list.push_back(&plate);
                params.strategy = SaveStrategy::Zip64 | SaveStrategy::Silence | SaveStrategy::SkipAuxiliary;
                REQUIRE(store_bbs_3mf(params));
                REQUIRE(fs::is_regular_file(file.filename()));
                REQUIRE(fs::file_size(file.filename()) > 0);
                require_source_unchanged();
                std::set<std::string> guids{source.part_guid()};
                for (const auto *part : parts) {
                    REQUIRE_FALSE(part->part_guid().empty());
                    CHECK(guids.insert(part->part_guid()).second);
                }

                DynamicPrintConfig loaded_config;
                ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Disable};
                ImportedBbsData imported;
                bool is_bambu = false;
                Semver version;
                const LoadStrategy strategy = LoadStrategy::LoadModel | LoadStrategy::LoadConfig | LoadStrategy::Silence;
                REQUIRE(load_bbs_3mf(file.filename().c_str(), &loaded_config, &substitutions, &loaded,
                                     &imported.plates, &imported.presets, &is_bambu, &version, nullptr, strategy));
                CHECK(is_bambu);
                CHECK(imported.plates.size() == 1);
                CHECK(imported.presets.empty());
                REQUIRE(loaded.objects.size() == 1);
                const ModelObject &reopened = *loaded.objects.front();
                REQUIRE(reopened.instances.size() == 1);
                REQUIRE(reopened.volumes.size() == 3);
                CHECK(reopened.name == fixture.object->name);
                const ModelVolume &reopened_source = volume_named(reopened, original_name);
                CHECK(reopened_source.type() == ModelVolumeType::MODEL_PART);
                CHECK(reopened_source.part_guid() == original_guid);
                CHECK(reopened_source.extruder_id() == extruder);
                CHECK(reopened_source.source.input_file == "synthetic-two-rings.stl");
                REQUIRE(reopened_source.config.has("wall_loops"));
                CHECK(reopened_source.config.opt_int("wall_loops") == 4);
                REQUIRE(reopened.config.has("layer_height"));
                CHECK(reopened.config.opt_float("layer_height") == Approx(0.24));
                const auto reopened_world = world_mesh(reopened_source, *reopened.instances.front());
                require_same_world_geometry(reopened_world, source_world, 1e-3);
                const auto reopened_shells = its_split(reopened_world);
                REQUIRE(reopened_shells.size() == 2);
                for (const auto &shell : reopened_shells) require_closed_ring(shell);
                for (size_t patch = 0; patch < parts.size(); ++patch) {
                    const ModelVolume &part = volume_named(reopened, parts[patch]->name);
                    CHECK(&part != &reopened_source);
                    CHECK(part.type() == ModelVolumeType::MODEL_PART);
                    CHECK(part.part_guid() == parts[patch]->part_guid());
                    REQUIRE(part.config.has("extruder"));
                    CHECK(part.config.opt_int("extruder") == extruder);
                    CHECK_FALSE(part.config.has("wall_loops"));
                    CHECK(fs::path(part.source.input_file).extension() == ".3mf");
                    require_part_geometry(part, *reopened.instances.front(), patch, 1e-3);
                }
                require_source_unchanged();
            }
            CHECK_FALSE(fs::exists(temporary_directory));
        }
}

namespace {

using AnnotationMember = FacetsAnnotation ModelVolume::*;
const std::array<AnnotationMember, 5> annotation_members{{
    &ModelVolume::supported_facets, &ModelVolume::seam_facets, &ModelVolume::mmu_segmentation_facets,
    &ModelVolume::fuzzy_skin_facets, &ModelVolume::exterior_facets
}};

void seed_annotation_layers(ModelVolume &volume, bool exterior = true)
{
    for (size_t i = 0; i < (exterior ? 5u : 4u); ++i) {
        auto &layer = volume.*annotation_members[i];
        REQUIRE(layer.empty());
        // Valid unsplit leaf states: 4 is state 1, 8 is state 2.
        layer.set_triangle_from_string(int(i), i == 1 ? "8" : "4");
    }
}

void seed_legacy_annotations(ModelVolume &volume)
{
    TriangleMesh mesh = volume.mesh();
    mesh.its.properties.assign(mesh.its.indices.size(), FaceProperty{eNormal, 0.0});
    mesh.its.properties[0] = {eSmallOverhang, 4.0};
    mesh.its.properties[1] = {eSmallHole, 8.0};
    mesh.its.properties[2] = {eExteriorAppearance, 16.0};
    mesh.its.properties[3] = {eMaxNumFaceTypes, 32.0}; // Unknown/sentinel is not silently ignored.
    mesh.its.properties[7] = {eNormal, 17.5}; // Normal metadata must remain untouched.
    volume.set_mesh(std::move(mesh));
}

std::array<bool, 5> annotation_flags(const FaceExtrusionAnnotationInfo &info)
{
    return {{info.support, info.seam, info.mmu, info.fuzzy_skin, info.exterior}};
}

} // namespace

TEST_CASE("Face annotation inspection reports only target payloads and ignores normal defaults", "[face_extrusion_model][face_extrusion_annotations]")
{
    ModelFixture fixture(0);
    auto &source = *fixture.source;
    const TriangleMesh *const mesh_ptr = source.mesh_ptr();
    const auto initial_mesh = source.mesh().its;
    const auto initial_stamp = capture_face_extrusion_annotation_stamp(source);
    REQUIRE_FALSE(initial_mesh.properties.empty());
    CHECK_FALSE(inspect_face_extrusion_annotations(source).any());
    std::array<ObjectBase::Timestamp, 5> timestamps{};
    for (size_t i = 0; i < annotation_members.size(); ++i)
        timestamps[i] = (source.*annotation_members[i]).timestamp();
    clear_face_extrusion_annotations(source);
    CHECK(source.mesh_ptr() == mesh_ptr);
    CHECK(capture_face_extrusion_annotation_stamp(source) == initial_stamp);
    require_exact_mesh(source.mesh().its, initial_mesh);
    for (size_t i = 0; i < annotation_members.size(); ++i)
        CHECK((source.*annotation_members[i]).timestamp() == timestamps[i]);

    auto *other = fixture.object->add_volume_with_shared_mesh(source, ModelVolumeType::PARAMETER_MODIFIER);
    seed_annotation_layers(*other);
    const auto other_stamp = capture_face_extrusion_annotation_stamp(*other);
    CHECK_FALSE(inspect_face_extrusion_annotations(source).any());
    clear_face_extrusion_annotations(source);
    CHECK(capture_face_extrusion_annotation_stamp(*other) == other_stamp);
    CHECK(other->is_modifier());

    for (size_t active = 0; active < annotation_members.size(); ++active) {
        auto &layer = source.*annotation_members[active];
        const ObjectID id = layer.id();
        layer.set_triangle_from_string(0, "4");
        std::array<ObjectBase::Timestamp, 5> before_clear{};
        for (size_t i = 0; i < annotation_members.size(); ++i)
            before_clear[i] = (source.*annotation_members[i]).timestamp();
        const auto info = inspect_face_extrusion_annotations(source);
        std::array<bool, 5> expected{};
        expected[active] = true;
        CHECK(annotation_flags(info) == expected);
        CHECK(info.legacy == 0);
        CHECK(info.any());
        clear_face_extrusion_annotations(source);
        CHECK_FALSE(inspect_face_extrusion_annotations(source).any());
        CHECK(layer.id() == id);
        CHECK(layer.empty());
        for (size_t i = 0; i < annotation_members.size(); ++i) {
            const bool changed = (source.*annotation_members[i]).timestamp() != before_clear[i];
            CHECK(changed == (i == active)); // Empty unrelated layers must not be touched.
        }
        CHECK(source.mesh_ptr() == mesh_ptr); // No legacy data: no geometry copy.
        CHECK(capture_face_extrusion_annotation_stamp(*other) == other_stamp);
    }
    seed_legacy_annotations(source);
    const auto legacy = inspect_face_extrusion_annotations(source);
    CHECK(legacy.legacy == 4);
    CHECK(legacy.legacy_small_overhang == 1);
    CHECK(legacy.legacy_small_hole == 1);
    CHECK(legacy.legacy_exterior == 1);
    CHECK(legacy.legacy_other == 1);
    CHECK(legacy.any());
    const std::array<bool, 5> none{};
    CHECK(annotation_flags(legacy) == none);
}

TEST_CASE("Face annotation stamps detect equal-count equal-timestamp content changes", "[face_extrusion_model][face_extrusion_annotations]")
{
    for (AnnotationMember member : annotation_members) {
        ModelFixture fixture(0);
        auto *a = fixture.source;
        auto *b = fixture.object->add_volume_with_shared_mesh(*a);
        auto *c = fixture.object->add_volume_with_shared_mesh(*a);
        (a->*member).set_triangle_from_string(0, "4");
        (b->*member).set_triangle_from_string(0, "8");
        (c->*member).set_triangle_from_string(1, "4");
        REQUIRE((a->*member).timestamp() == (b->*member).timestamp());
        REQUIRE((a->*member).timestamp() == (c->*member).timestamp());
        REQUIRE((a->*member).get_data().first.size() == (b->*member).get_data().first.size());
        REQUIRE((a->*member).get_data().second.size() == (b->*member).get_data().second.size());
        CHECK(inspect_face_extrusion_annotations(*a) == inspect_face_extrusion_annotations(*b));
        CHECK(inspect_face_extrusion_annotations(*a) == inspect_face_extrusion_annotations(*c));
        const auto stamp = capture_face_extrusion_annotation_stamp(*a);
        CHECK(stamp == capture_face_extrusion_annotation_stamp(*a));
        CHECK(stamp != capture_face_extrusion_annotation_stamp(*b)); // Different encoded leaf state.
        CHECK(stamp != capture_face_extrusion_annotation_stamp(*c)); // Different triangle ID.
        const auto copied_stamp = stamp;
        CHECK(copied_stamp == stamp);
    }
    ModelFixture fixture(0);
    auto &source = *fixture.source;
    auto mesh = source.mesh().its;
    mesh.properties[0] = {eSmallHole, 4.0};
    source.set_mesh(mesh);
    const auto before = capture_face_extrusion_annotation_stamp(source);
    const auto info = inspect_face_extrusion_annotations(source);
    mesh.properties[0].area = 5.0;
    source.set_mesh(mesh);
    CHECK(inspect_face_extrusion_annotations(source) == info);
    CHECK(capture_face_extrusion_annotation_stamp(source) != before);
    mesh.properties[0] = {eNormal, 0.0};
    mesh.properties[1] = {eSmallHole, 4.0};
    source.set_mesh(mesh);
    CHECK(inspect_face_extrusion_annotations(source) == info);
    CHECK(capture_face_extrusion_annotation_stamp(source) != before);
    mesh.properties[1].area = std::numeric_limits<double>::quiet_NaN();
    source.set_mesh(mesh);
    const auto nan_stamp = capture_face_extrusion_annotation_stamp(source);
    CHECK(nan_stamp == capture_face_extrusion_annotation_stamp(source));
    mesh.properties.assign(mesh.indices.size(), FaceProperty{eNormal, 0.0});
    source.set_mesh(mesh);
    const auto normal_stamp = capture_face_extrusion_annotation_stamp(source);
    mesh.properties[0].area = -0.0;
    source.set_mesh(mesh);
    CHECK_FALSE(inspect_face_extrusion_annotations(source).any());
    CHECK(capture_face_extrusion_annotation_stamp(source) != normal_stamp); // Exact area representation.
    const TriangleMesh *const ptr = source.mesh_ptr();
    const auto signed_zero_stamp = capture_face_extrusion_annotation_stamp(source);
    clear_face_extrusion_annotations(source);
    CHECK(source.mesh_ptr() == ptr);
    CHECK(capture_face_extrusion_annotation_stamp(source) == signed_zero_stamp);
}

TEST_CASE("Face annotation clear is copy-on-write and preserves model identity and unrelated metadata", "[face_extrusion_model][face_extrusion_annotations]")
{
    ModelFixture fixture(3);
    auto &source = *fixture.source;
    seed_annotation_layers(source);
    seed_legacy_annotations(source);
    source.set_material_id("annotation-material");
    TextInfo text;
    text.m_text = "Keep procedural text";
    text.m_font_name = "Test font";
    text.m_thickness = 3.5f;
    source.set_text_info(text);
    source.emboss_shape.emplace();
    source.emboss_shape->scale = 0.002;
    source.emboss_shape->svg_file.emplace();
    source.emboss_shape->svg_file->path = "metadata-only.svg";
    source.cut_info.is_connector = true;
    source.cut_info.radius = 3.f;
    source.cut_info.height = 6.f;
    source.cut_info.source_part_idx = 17;
    source.source.mesh_offset = Vec3d(1, 2, 3);
    source.source.is_converted_from_inches = true;
    source.mmuseg_ts = 0;
    CHECK(source.get_extruders() == std::vector<int>{1, 2});
    auto *modifier = fixture.object->add_volume_with_shared_mesh(source, ModelVolumeType::PARAMETER_MODIFIER);
    modifier->name = "Untouched modifier";
    modifier->config.set_key_value("wall_loops", new ConfigOptionInt(7));
    seed_annotation_layers(*modifier);
    const auto source_stamp = capture_face_extrusion_annotation_stamp(source);
    const auto modifier_stamp = capture_face_extrusion_annotation_stamp(*modifier);
    const auto source_mesh = source.mesh().its;
    const auto source_config = config_values(source.config);
    const auto object_config = config_values(fixture.object->config);
    const auto modifier_config = config_values(modifier->config);
    Model copied = fixture.model;
    auto &target = *copied.objects.front()->volumes.front();
    auto &other = *copied.objects.front()->volumes[1];
    REQUIRE(target.mesh_ptr() == source.mesh_ptr());
    REQUIRE(other.mesh_ptr() == source.mesh_ptr());
    REQUIRE(capture_face_extrusion_annotation_stamp(target) == source_stamp);
    std::array<ObjectID, 5> ids{};
    std::array<ObjectBase::Timestamp, 5> timestamps{};
    for (size_t i = 0; i < annotation_members.size(); ++i) {
        ids[i] = (target.*annotation_members[i]).id();
        timestamps[i] = (target.*annotation_members[i]).timestamp();
    }
    const ObjectID target_id = target.id();
    const ObjectID config_id = target.config.id();
    const auto config_timestamp = static_cast<const ModelConfig &>(target.config).timestamp();
    const TriangleMesh *const old_mesh_ptr = target.mesh_ptr();
    clear_face_extrusion_annotations(target);
    CHECK_FALSE(inspect_face_extrusion_annotations(target).any());
    CHECK(target.mesh_ptr() != old_mesh_ptr);
    auto expected_mesh = source_mesh;
    for (auto &property : expected_mesh.properties)
        if (property.type != eNormal) property = {eNormal, 0.0};
    require_exact_mesh(target.mesh().its, expected_mesh);
    require_exact_mesh(source.mesh().its, source_mesh);
    require_exact_mesh(other.mesh().its, source_mesh);
    CHECK(source.mesh_ptr() == old_mesh_ptr);
    CHECK(other.mesh_ptr() == old_mesh_ptr);
    CHECK(capture_face_extrusion_annotation_stamp(source) == source_stamp);
    CHECK(capture_face_extrusion_annotation_stamp(*modifier) == modifier_stamp);
    CHECK(capture_face_extrusion_annotation_stamp(other) == modifier_stamp);
    CHECK(target.id() == target_id);
    CHECK(target.config.id() == config_id);
    CHECK(static_cast<const ModelConfig &>(target.config).timestamp() == config_timestamp);
    CHECK(config_values(target.config) == source_config);
    CHECK(config_values(copied.objects.front()->config) == object_config);
    CHECK(config_values(other.config) == modifier_config);
    CHECK(other.is_modifier());
    CHECK(other.name == modifier->name);
    CHECK(other.id() == modifier->id());
    CHECK(target.type() == source.type());
    CHECK(target.name == source.name);
    CHECK(target.material_id() == source.material_id());
    CHECK(target.part_guid() == source.part_guid());
    CHECK(target.get_matrix().matrix() == source.get_matrix().matrix());
    CHECK(copied.objects.front()->instances.front()->get_matrix().matrix() == fixture.instance->get_matrix().matrix());
    CHECK(target.source.input_file == source.source.input_file);
    CHECK(target.source.object_idx == source.source.object_idx);
    CHECK(target.source.volume_idx == source.source.volume_idx);
    CHECK(target.source.mesh_offset == source.source.mesh_offset);
    CHECK(target.source.is_converted_from_inches == source.source.is_converted_from_inches);
    CHECK(target.get_text_info().m_text == text.m_text);
    CHECK(target.get_text_info().m_font_name == text.m_font_name);
    CHECK(target.get_text_info().m_thickness == text.m_thickness);
    REQUIRE(target.emboss_shape.has_value());
    CHECK(target.emboss_shape->scale == 0.002);
    REQUIRE(target.emboss_shape->svg_file.has_value());
    CHECK(target.emboss_shape->svg_file->path == "metadata-only.svg");
    CHECK(target.cut_info.is_connector);
    CHECK(target.cut_info.radius == 3.f);
    CHECK(target.cut_info.height == 6.f);
    CHECK(target.cut_info.source_part_idx == 17);
    CHECK(target.mmuseg_extruders.empty());
    CHECK(target.mmuseg_ts == 0);
    CHECK(target.get_extruders() == std::vector<int>{2});
    for (size_t i = 0; i < annotation_members.size(); ++i) {
        CHECK((target.*annotation_members[i]).id() == ids[i]);
        CHECK((target.*annotation_members[i]).timestamp() != timestamps[i]);
        timestamps[i] = (target.*annotation_members[i]).timestamp();
    }
    const auto cleared_stamp = capture_face_extrusion_annotation_stamp(target);
    const TriangleMesh *const cleared_ptr = target.mesh_ptr();
    clear_face_extrusion_annotations(target);
    CHECK(target.mesh_ptr() == cleared_ptr);
    CHECK(capture_face_extrusion_annotation_stamp(target) == cleared_stamp);
    for (size_t i = 0; i < annotation_members.size(); ++i)
        CHECK((target.*annotation_members[i]).timestamp() == timestamps[i]);
}

TEST_CASE("Face annotation consent gate clears only a matching successful model commit", "[face_extrusion_model][face_extrusion_annotations]")
{
    // Exercise the real production clear/stamp helpers behind a caller-side
    // commit guard. This is not a GUI, selection-list or Undo simulation.
    enum class Outcome { NoConsent, Failure, Canceled, Stale, Commit };
    for (Outcome outcome : {Outcome::NoConsent, Outcome::Failure, Outcome::Canceled, Outcome::Stale, Outcome::Commit}) {
        ModelFixture fixture(0);
        auto &target = *fixture.source;
        seed_annotation_layers(target);
        seed_legacy_annotations(target);
        const auto consent_stamp = capture_face_extrusion_annotation_stamp(target);
        const auto original_info = inspect_face_extrusion_annotations(target);
        if (outcome == Outcome::Stale) {
            target.seam_facets.reset();
            target.seam_facets.set_triangle_from_string(1, "4"); // Same flags/count, different state.
            REQUIRE(inspect_face_extrusion_annotations(target) == original_info);
            REQUIRE(capture_face_extrusion_annotation_stamp(target) != consent_stamp);
        }
        const auto current_stamp = capture_face_extrusion_annotation_stamp(target);
        const TriangleMesh *const original_ptr = target.mesh_ptr();
        const auto original_config = config_values(target.config);
        const auto original_guid = target.part_guid();
        const auto original_pose = target.get_matrix();
        const bool consent = outcome != Outcome::NoConsent;
        const bool success = outcome != Outcome::Failure;
        const bool canceled = outcome == Outcome::Canceled;
        bool committed = false;
        if (success && consent && !canceled && capture_face_extrusion_annotation_stamp(target) == consent_stamp) {
            clear_face_extrusion_annotations(target);
            committed = true;
        }
        CHECK(committed == (outcome == Outcome::Commit));
        if (!committed) {
            CHECK(target.mesh_ptr() == original_ptr);
            CHECK(capture_face_extrusion_annotation_stamp(target) == current_stamp);
        } else {
            CHECK_FALSE(inspect_face_extrusion_annotations(target).any());
        }
        CHECK(config_values(target.config) == original_config);
        CHECK(target.part_guid() == original_guid);
        CHECK(target.get_matrix().matrix() == original_pose.matrix());
    }
}

TEST_CASE("Cleared target and untouched painted neighbor survive native BBS 3MF", "[face_extrusion_model][face_extrusion_annotations]")
{
    fs::path temporary_directory;
    {
        Temporary3mfDirectory file;
        temporary_directory = file.directory();
        ModelFixture fixture(2);
        Model loaded;
        ScopedModelScratch source_scratch(fixture.model, file.directory());
        ScopedModelScratch loaded_scratch(loaded, file.directory());
        auto &target = *fixture.source;
        target.name = "Explicitly cleared target";
        auto mesh = target.mesh().its;
        mesh.properties[2] = {eSmallHole, 9.0};
        target.set_mesh(mesh);
        seed_annotation_layers(target);
        auto *other = fixture.object->add_volume_with_shared_mesh(target);
        other->name = "Untouched painted neighbor";
        other->config.set_key_value("extruder", new ConfigOptionInt(2));
        other->set_transformation(target.get_transformation());
        other->set_offset(other->get_offset() + Vec3d(35, 0, 0));
        // Native BBS serializes support/seam/MMU/fuzzy and legacy properties.
        // Exterior FacetsAnnotation is not serialized by this existing format;
        // its in-memory preservation is covered above, not falsely claimed here.
        seed_annotation_layers(*other, false);
        other->mmuseg_ts = 0; // Explicitly initialize the fixture's derived cache.
        const auto other_stamp = capture_face_extrusion_annotation_stamp(*other);
        const auto target_world = world_mesh(target, *fixture.instance);
        const auto other_world = world_mesh(*other, *fixture.instance);
        const auto target_config = config_values(target.config);
        const auto target_guid = target.part_guid();
        clear_face_extrusion_annotations(target); // Explicit successful-commit consent in this fixture.
        REQUIRE_FALSE(inspect_face_extrusion_annotations(target).any());
        REQUIRE(capture_face_extrusion_annotation_stamp(*other) == other_stamp);
        DynamicPrintConfig config;
        config.apply(DynamicPrintConfig::full_print_config());
        config.set_key_value("filament_settings_id", new ConfigOptionStrings({"Test filament 1", "Test filament 2"}));
        config.set_key_value("filament_colour", new ConfigOptionStrings({"#FFFFFF", "#000000"}));
        config.set_key_value("filament_type", new ConfigOptionStrings({"PLA", "PLA"}));
        PlateData plate;
        plate.plate_index = 0;
        plate.plate_name = "Annotation consent test";
        plate.objects_and_instances.emplace_back(0, 0);
        StoreParams params;
        params.path = file.filename().c_str();
        params.model = &fixture.model;
        params.config = &config;
        params.plate_data_list.push_back(&plate);
        params.strategy = SaveStrategy::Zip64 | SaveStrategy::Silence | SaveStrategy::SkipAuxiliary;
        REQUIRE(store_bbs_3mf(params));
        CHECK(config_values(target.config) == target_config);
        CHECK(target.part_guid() == target_guid);
        CHECK(capture_face_extrusion_annotation_stamp(*other) == other_stamp);
        DynamicPrintConfig loaded_config;
        ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Disable};
        ImportedBbsData imported;
        bool is_bambu = false;
        Semver version;
        const LoadStrategy strategy = LoadStrategy::LoadModel | LoadStrategy::LoadConfig | LoadStrategy::Silence;
        REQUIRE(load_bbs_3mf(file.filename().c_str(), &loaded_config, &substitutions, &loaded,
                             &imported.plates, &imported.presets, &is_bambu, &version, nullptr, strategy));
        REQUIRE(is_bambu);
        REQUIRE(loaded.objects.size() == 1);
        const auto &object = *loaded.objects.front();
        REQUIRE(object.volumes.size() == 2);
        REQUIRE(object.instances.size() == 1);
        const auto &clean = volume_named(object, target.name);
        const auto &painted = volume_named(object, other->name);
        CHECK_FALSE(inspect_face_extrusion_annotations(clean).any());
        CHECK(capture_face_extrusion_annotation_stamp(painted) == other_stamp);
        CHECK(clean.part_guid() == target_guid);
        CHECK(clean.config.opt_int("extruder") == 2);
        CHECK(clean.config.opt_int("wall_loops") == 4);
        require_same_world_geometry(world_mesh(clean, *object.instances.front()), target_world, 1e-3);
        require_same_world_geometry(world_mesh(painted, *object.instances.front()), other_world, 1e-3);
    }
    CHECK_FALSE(fs::exists(temporary_directory));
}
