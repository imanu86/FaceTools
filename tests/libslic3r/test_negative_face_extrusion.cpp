// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>

#include "libslic3r/FaceExtrusion.hpp"
#include "libslic3r/FaceExtrusionAnnotations.hpp"
#include "libslic3r/Format/bbs_3mf.hpp"
#include "libslic3r/MeshBoolean.hpp"
#include "libslic3r/MeshDiagnostics.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Preset.hpp"
#include "slic3r/GUI/Gizmos/FaceExtrusionModelPolicy.hpp"

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
using namespace Slic3r::GUI;

namespace {
namespace fs = boost::filesystem;
using Mode = FaceExtrusionToolMode;
using Role = ModelVolumeType;
using Layer = FacetsAnnotation ModelVolume::*;
const std::array<Layer, 5> layers{{&ModelVolume::supported_facets, &ModelVolume::seam_facets,
    &ModelVolume::mmu_segmentation_facets, &ModelVolume::fuzzy_skin_facets, &ModelVolume::exterior_facets}};

// Independent outward 20 x 20 x 2 annular prism; cap area 300, volume 600.
// Original top-cap facets are [8,16). No user mesh or project is loaded.
indexed_triangle_set ring()
{
    indexed_triangle_set mesh;
    mesh.vertices = {{-10,-10,0},{10,-10,0},{10,10,0},{-10,10,0},
                     {-5,-5,0},{5,-5,0},{5,5,0},{-5,5,0}};
    mesh.vertices.reserve(16);
    for (int i = 0; i < 8; ++i) mesh.vertices.emplace_back(mesh.vertices[i].x(), mesh.vertices[i].y(), 2.f);
    const std::vector<Vec3i> cap{{0,1,5},{0,5,4},{1,2,6},{1,6,5},{2,3,7},{2,7,6},{3,0,4},{3,4,7}};
    for (const auto &f : cap) mesh.indices.emplace_back(f[0],f[2],f[1]);
    for (const auto &f : cap) mesh.indices.emplace_back(f[0]+8,f[1]+8,f[2]+8);
    const auto wall = [&](int a, int b) {
        mesh.indices.emplace_back(a,b,b+8);
        mesh.indices.emplace_back(a,b+8,a+8);
    };
    for (int i = 0; i < 4; ++i) { wall(i,(i+1)%4); wall(4+(i+1)%4,4+i); }
    mesh.properties.assign(mesh.indices.size(), FaceProperty{eNormal,0.});
    return mesh;
}

TriangleMesh cube(double x, double y, double z)
{
    TriangleMesh mesh = make_cube(x,y,z);
    mesh.its.properties.assign(mesh.its.indices.size(), FaceProperty{eNormal,0.});
    return mesh;
}

std::vector<int> top_faces(int offset = 0)
{
    std::vector<int> faces(8);
    std::iota(faces.begin(), faces.end(), offset + 8);
    return faces;
}

std::map<std::string,std::string> config_values(const ModelConfig &config)
{
    std::map<std::string,std::string> result;
    for (const auto &key : config.keys()) result.emplace(key,config.opt_serialize(key));
    return result;
}

bool same_matrix(const Transform3d &a, const Transform3d &b)
{ return (a.matrix().array() == b.matrix().array()).all(); }

void check_mesh(const indexed_triangle_set &actual, const indexed_triangle_set &expected)
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

void require_solid(const indexed_triangle_set &mesh, bool annular = false)
{
    const auto diagnostics = its_mesh_diagnostics(mesh);
    CHECK(diagnostics.open_edges == 0);
    CHECK(diagnostics.non_manifold_edges == 0);
    CHECK(diagnostics.non_manifold_vertices == 0);
    CHECK(diagnostics.same_direction_edges == 0);
    CHECK_FALSE(diagnostics.has_reversed_faces);
    CHECK(its_volume(mesh) > 0.f); // A NEGATIVE_VOLUME is still an outward solid.
    if (annular) {
        std::set<std::pair<int,int>> edges;
        for (const auto &face : mesh.indices)
            for (int j = 0; j < 3; ++j) edges.insert(std::minmax(face[j],face[(j+1)%3]));
        CHECK(long(mesh.vertices.size()) - long(edges.size()) + long(mesh.indices.size()) == 0);
        CHECK(its_split(mesh).size() == 1);
    }
}

indexed_triangle_set world_mesh(const ModelVolume &volume, const ModelInstance &instance)
{
    auto mesh = volume.mesh().its;
    its_transform(mesh, instance.get_matrix() * volume.get_matrix(), true);
    return mesh;
}

// Native import may reorder/recenter vertices; compare all world triangles, not
// merely volume/bounds. The existing Faces budgets stay 1e-4 (storage), 1e-3 (3MF).
void require_world_geometry(const indexed_triangle_set &actual, const indexed_triangle_set &expected, double tolerance)
{
    REQUIRE(actual.vertices.size() == expected.vertices.size());
    REQUIRE(actual.indices.size() == expected.indices.size());
    std::vector<int> remap(actual.vertices.size(), -1);
    std::vector<bool> used(expected.vertices.size(), false);
    for (size_t i = 0; i < actual.vertices.size(); ++i) {
        double best_distance = std::numeric_limits<double>::max();
        size_t best = expected.vertices.size();
        for (size_t j = 0; j < expected.vertices.size(); ++j) {
            const double distance = (actual.vertices[i].cast<double>() - expected.vertices[j].cast<double>()).norm();
            if (!used[j] && distance < best_distance) { best = j; best_distance = distance; }
        }
        INFO("world vertex " << i << ", error " << best_distance);
        REQUIRE(best < expected.vertices.size());
        REQUIRE(best_distance <= tolerance);
        remap[i] = int(best);
        used[best] = true;
    }
    std::multiset<std::array<int,3>> a, b;
    for (const auto &face : actual.indices) {
        std::array<int,3> indices{{remap.at(face[0]),remap.at(face[1]),remap.at(face[2])}};
        std::sort(indices.begin(),indices.end()); a.insert(indices);
    }
    for (const auto &face : expected.indices) {
        std::array<int,3> indices{{face[0],face[1],face[2]}};
        std::sort(indices.begin(),indices.end()); b.insert(indices);
    }
    CHECK(a == b);
    require_solid(actual);
}

void paint_four_layers(ModelVolume &volume)
{
    for (size_t i = 0; i < 4; ++i) (volume.*layers[i]).set_triangle_from_string(int(i), i == 1 ? "8" : "4");
    volume.mmuseg_ts = 0;
}

struct SavedVolume {
    const ModelVolume *volume;
    ObjectID id, config_id;
    Role role;
    std::shared_ptr<const TriangleMesh> mesh, hull;
    indexed_triangle_set geometry;
    Transform3d matrix;
    std::map<std::string,std::string> config;
    std::string name, material, guid;
    ModelVolume::Source source;
    FaceExtrusionAnnotationStamp annotations;
    std::array<ObjectID,5> layer_ids;
    std::array<ObjectBase::Timestamp,5> timestamps;

    explicit SavedVolume(const ModelVolume &v)
        : volume(&v), id(v.id()), config_id(v.config.id()), role(v.type()),
          mesh(v.get_mesh_shared_ptr()), hull(v.get_convex_hull_shared_ptr()), geometry(v.mesh().its),
          matrix(v.get_matrix()), config(config_values(v.config)), name(v.name), material(v.material_id()),
          guid(v.part_guid()), source(v.source), annotations(capture_face_extrusion_annotation_stamp(v))
    {
        for (size_t i = 0; i < layers.size(); ++i) {
            layer_ids[i] = (v.*layers[i]).id(); timestamps[i] = (v.*layers[i]).timestamp();
        }
    }

    void metadata() const
    {
        CHECK(volume->type() == role); CHECK(volume->name == name);
        CHECK(config_values(volume->config) == config); CHECK(volume->material_id() == material);
        CHECK(volume->part_guid() == guid); CHECK(same_matrix(volume->get_matrix(),matrix));
    }

    void unchanged() const
    {
        metadata();
        CHECK(volume->id() == id); CHECK(volume->config.id() == config_id);
        CHECK(volume->get_mesh_shared_ptr() == mesh); CHECK(volume->get_convex_hull_shared_ptr() == hull);
        check_mesh(volume->mesh().its,geometry);
        CHECK(capture_face_extrusion_annotation_stamp(*volume) == annotations);
        for (size_t i = 0; i < layers.size(); ++i) {
            CHECK((volume->*layers[i]).id() == layer_ids[i]);
            CHECK((volume->*layers[i]).timestamp() == timestamps[i]);
        }
        CHECK(volume->source.input_file == source.input_file);
        CHECK(volume->source.object_idx == source.object_idx); CHECK(volume->source.volume_idx == source.volume_idx);
        CHECK(volume->source.mesh_offset == source.mesh_offset);
        CHECK(volume->source.is_from_builtin_objects == source.is_from_builtin_objects);
        CHECK(volume->source.is_converted_from_inches == source.is_converted_from_inches);
        CHECK(volume->source.is_converted_from_meters == source.is_converted_from_meters);
        CHECK(same_matrix(volume->source.transform.get_matrix(),source.transform.get_matrix()));
    }
};

struct Fixture {
    Model model;
    ModelObject *object;
    ModelInstance *instance;
    ModelVolume *source, *carrier, *neighbor;
    std::vector<std::vector<int>> patches;

    Fixture(Role type = Role::NEGATIVE_VOLUME, bool multiple = false, int mirror_mask = -1)
    {
        object = model.add_object(); object->name = "Synthetic negative face fixture";
        object->config.set_key_value("extruder",new ConfigOptionInt(1));
        object->config.set_key_value("layer_height",new ConfigOptionFloat(0.24));
        carrier = object->add_volume(cube(80,40,10),Role::MODEL_PART,false);
        carrier->name = "Untouched carrier"; carrier->set_offset(Vec3d(-15,-20,0)); carrier->ensure_part_guid();
        auto source_mesh = ring(); patches.push_back(top_faces());
        if (multiple) {
            auto second = ring();
            for (auto &p : second.vertices) { p.x() += 35.f; p.z() *= 2.f; }
            patches.push_back(top_faces(int(source_mesh.indices.size())));
            its_merge(source_mesh,second);
            source_mesh.properties.assign(source_mesh.indices.size(),FaceProperty{eNormal,0.});
        }
        source = object->add_volume(TriangleMesh(std::move(source_mesh)),type,false);
        source->name = "Selected cutter"; source->ensure_part_guid();
        source->config.set_key_value("extruder",new ConfigOptionInt(2));
        source->config.set_key_value("wall_loops",new ConfigOptionInt(4));
        source->set_material_id("synthetic-source-material");
        source->source.input_file = "synthetic-negative-ring.stl";
        source->source.object_idx = 3; source->source.volume_idx = 7;
        source->source.is_converted_from_inches = true;
        neighbor = object->add_volume(cube(3,4,5),Role::MODEL_PART,false);
        neighbor->name = "Untouched painted neighbor"; neighbor->set_offset(Vec3d(70,0,0));
        neighbor->config.set_key_value("extruder",new ConfigOptionInt(2)); neighbor->ensure_part_guid();
        paint_four_layers(*neighbor);
        for (Role sentinel_type : {Role::PARAMETER_MODIFIER,Role::SUPPORT_BLOCKER}) {
            auto *sentinel = object->add_volume(cube(2,2,2),sentinel_type,false);
            sentinel->name = sentinel_type == Role::PARAMETER_MODIFIER ? "Untouched modifier" : "Untouched support blocker";
            sentinel->set_offset(Vec3d(-30,0,-4));
            sentinel->config.set_key_value("wall_loops",new ConfigOptionInt(7));
        }
        instance = object->add_instance();
        if (mirror_mask >= 0) {
            Geometry::Transformation volume_pose;
            volume_pose.set_offset(Vec3d(8,-4,15)); volume_pose.set_rotation(Vec3d(0.3,0.4,0.2));
            volume_pose.set_scaling_factor(Vec3d(1.2,0.7,1.8));
            volume_pose.set_mirror(Vec3d(mirror_mask & 1 ? -1. : 1.,1.,1.));
            source->set_transformation(volume_pose);
            Geometry::Transformation instance_pose;
            instance_pose.set_offset(Vec3d(75,90,0)); instance_pose.set_rotation(Vec3d(-0.2,0.1,0.5));
            instance_pose.set_scaling_factor(Vec3d(0.8,1.3,0.9));
            instance_pose.set_mirror(Vec3d(1.,mirror_mask & 2 ? -1. : 1.,1.));
            instance->set_transformation(instance_pose);
            object->invalidate_bounding_box();
            Vec3d offset = instance->get_offset(); offset.z() -= object->get_min_z(); instance->set_offset(offset);
        }
        object->invalidate_bounding_box();
    }

    std::vector<SavedVolume> saved(bool include_source = true) const
    {
        std::vector<SavedVolume> result;
        for (const auto *v : object->volumes) if (include_source || v != source) result.emplace_back(*v);
        return result;
    }
};

struct Patch {
    TriangleMesh world, object;
    Vec3d displacement;
};
Patch prepare_patch(const Fixture &fixture, size_t patch, double distance)
{
    const Transform3d world = fixture.instance->get_matrix() * fixture.source->get_matrix();
    const Vec3d displacement = distance * (world.linear().inverse().transpose() * Vec3d::UnitZ()).normalized();
    const auto input = world_mesh(*fixture.source,*fixture.instance);
    auto extrusion = extrude_planar_face(input,fixture.patches.at(patch),displacement);
    INFO(extrusion.message); REQUIRE(extrusion);
    CHECK(extrusion.boundary_edges == 8); require_solid(extrusion.mesh,true);
    const size_t cap_size = extrusion.mesh.vertices.size()/2;
    REQUIRE(cap_size == 8);
    for (size_t i = 0; i < cap_size; ++i) {
        const Vec3d actual = extrusion.mesh.vertices[i+cap_size].cast<double>() - extrusion.mesh.vertices[i].cast<double>();
        CHECK((actual-displacement).norm() <= 1e-4);
    }
    Patch result;
    result.world = TriangleMesh(std::move(extrusion.mesh)); result.object = result.world; result.displacement = displacement;
    result.object.transform(fixture.instance->get_matrix().inverse(),true);
    for (size_t i = 0; i < result.object.its.vertices.size(); ++i)
        REQUIRE((fixture.instance->get_matrix()*result.object.its.vertices[i].cast<double>() - result.world.its.vertices[i].cast<double>()).norm() <= 1e-4);
    return result;
}

std::shared_ptr<const TriangleMesh> boolean_result(const Fixture &fixture, bool subtract)
{
    const auto source = world_mesh(*fixture.source,*fixture.instance);
    auto patch = prepare_patch(fixture,0,subtract ? -1. : 2.);
    std::vector<TriangleMesh> outputs; int failures = 0;
    MeshBoolean::mcut::make_boolean(TriangleMesh(source),patch.world,outputs,subtract ? "A_NOT_B" : "UNION",
        [] { return false; },nullptr,[&] { ++failures; });
    REQUIRE(failures == 0); REQUIRE(outputs.size() == 1);
    require_solid(outputs.front().its,true);
    // These fixtures are untransformed for the independent exact volume oracle.
    CHECK(its_volume(outputs.front().its) == Approx(subtract ? 300. : 1200.).epsilon(2e-5).margin(1e-3));
    auto solid = MeshBoolean::cgal::triangle_mesh_to_cgal(outputs.front().its);
    CHECK_FALSE(MeshBoolean::cgal::does_self_intersect(*solid)); CHECK(MeshBoolean::cgal::does_bound_a_volume(*solid));
    outputs.front().transform((fixture.instance->get_matrix()*fixture.source->get_matrix()).inverse(),true);
    return std::make_shared<const TriangleMesh>(std::move(outputs.front()));
}

void check_replacement_metadata(const SavedVolume &before)
{
    before.metadata();
    const auto &source = *before.volume;
    CHECK(source.id() != before.id); CHECK(source.config.id() != before.config_id);
    for (size_t i = 0; i < 4; ++i) CHECK((source.*layers[i]).id() != before.layer_ids[i]);
    CHECK(source.source.input_file.empty()); CHECK(source.source.object_idx == -1); CHECK(source.source.volume_idx == -1);
    CHECK(source.source.is_from_builtin_objects); CHECK_FALSE(source.source.is_converted_from_inches);
    CHECK_FALSE(source.source.is_converted_from_meters);
    CHECK(source.source.mesh_offset == Vec3d::Zero());
}

// Same private-scratch safety contract as the existing native BBS tests. Never
// reads TEST_DATA_DIR or an app profile; removes only canonical known children.
class Temporary3mfDirectory {
public:
    Temporary3mfDirectory()
    {
        const fs::path parent = fs::canonical(fs::temp_directory_path());
        for (int i = 0; i < 16; ++i) {
            const fs::path candidate = parent/fs::unique_path("bambu-negative-face-%%%%-%%%%-%%%%");
            if (fs::create_directory(candidate)) { directory = fs::canonical(candidate); file = (directory/"roundtrip.3mf").string(); return; }
        }
        throw std::runtime_error("Cannot create private negative face test directory");
    }
    ~Temporary3mfDirectory()
    {
        try {
            if (!verified()) return;
            for (const char *name : {"roundtrip.3mf","roundtrip.3mf.tmp","_temp_1.config","_temp_2.config","_temp_3.config"}) {
                const fs::path child = directory/name;
                if (fs::is_regular_file(fs::symlink_status(child)) && fs::canonical(child) == child) {
                    boost::system::error_code ec; fs::remove(child,ec);
                }
            }
            const fs::path auxiliary = directory/"Auxiliaries";
            if (fs::is_directory(fs::symlink_status(auxiliary)) && fs::canonical(auxiliary) == auxiliary && fs::is_empty(auxiliary)) {
                boost::system::error_code ec; fs::remove(auxiliary,ec);
            }
            if (verified() && fs::is_empty(directory)) { boost::system::error_code ec; fs::remove(directory,ec); }
        } catch (...) {} // Unexpected artifacts remain for diagnosis, never recursive deletion.
    }
    Temporary3mfDirectory(const Temporary3mfDirectory &) = delete;
    Temporary3mfDirectory &operator=(const Temporary3mfDirectory &) = delete;
    fs::path directory;
    std::string file;
private:
    bool verified() const { return directory.is_absolute() && fs::is_directory(fs::symlink_status(directory)) && fs::canonical(directory) == directory; }
};

class ScopedModelScratch {
public:
    ScopedModelScratch(Model &m, const fs::path &path) : model(m) { model.set_backup_path(path.string()); }
    ~ScopedModelScratch() { model.set_backup_path("detach"); }
    ScopedModelScratch(const ScopedModelScratch &) = delete;
    ScopedModelScratch &operator=(const ScopedModelScratch &) = delete;
private:
    Model &model;
};

struct ImportedBbsData {
    PlateDataPtrs plates;
    std::vector<Preset *> presets;
    ~ImportedBbsData() { release_PlateData_list(plates); for (auto *preset : presets) delete preset; }
};

const ModelVolume &volume_named(const ModelObject &object, const std::string &name)
{
    const auto it = std::find_if(object.volumes.begin(),object.volumes.end(),[&](const ModelVolume *v) { return v->name == name; });
    REQUIRE(it != object.volumes.end()); return **it;
}

void require_bbs_roundtrip(Fixture &fixture)
{
    fs::path directory;
    {
        Temporary3mfDirectory temporary; directory = temporary.directory;
        Model loaded;
        ScopedModelScratch source_scratch(fixture.model,temporary.directory), loaded_scratch(loaded,temporary.directory);
        const auto unchanged = fixture.saved();
        DynamicPrintConfig config; config.apply(DynamicPrintConfig::full_print_config());
        config.set_key_value("filament_settings_id",new ConfigOptionStrings({"Synthetic filament 1","Synthetic filament 2"}));
        config.set_key_value("filament_colour",new ConfigOptionStrings({"#FFFFFF","#000000"}));
        config.set_key_value("filament_type",new ConfigOptionStrings({"PLA","PLA"}));
        PlateData plate; plate.plate_index = 0; plate.plate_name = "Synthetic negative face plate";
        plate.objects_and_instances.emplace_back(0,0);
        StoreParams params; params.path = temporary.file.c_str(); params.model = &fixture.model; params.config = &config;
        params.plate_data_list.push_back(&plate); params.strategy = SaveStrategy::Zip64 | SaveStrategy::Silence | SaveStrategy::SkipAuxiliary;
        REQUIRE(store_bbs_3mf(params));
        for (const auto &saved : unchanged) saved.unchanged();
        DynamicPrintConfig loaded_config;
        ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Disable};
        ImportedBbsData imported; bool is_bambu = false; Semver version;
        REQUIRE(load_bbs_3mf(temporary.file.c_str(),&loaded_config,&substitutions,&loaded,&imported.plates,&imported.presets,
                            &is_bambu,&version,nullptr,LoadStrategy::LoadModel | LoadStrategy::LoadConfig | LoadStrategy::Silence));
        CHECK(is_bambu); CHECK(imported.plates.size() == 1); CHECK(imported.presets.empty());
        REQUIRE(loaded.objects.size() == 1);
        const ModelObject &object = *loaded.objects.front(); REQUIRE(object.instances.size() == 1);
        REQUIRE(object.volumes.size() == fixture.object->volumes.size());
        CHECK(object.parts_count() == fixture.object->parts_count());
        for (const auto *original : fixture.object->volumes) {
            const auto &reopened = volume_named(object,original->name);
            CHECK(reopened.type() == original->type());
            require_world_geometry(world_mesh(reopened,*object.instances.front()),world_mesh(*original,*fixture.instance),1e-3);
            if (!original->part_guid().empty()) CHECK(reopened.part_guid() == original->part_guid());
            // Imported cache IDs and storage transforms need not equal the originals.
        }
        const auto &source = volume_named(object,fixture.source->name);
        REQUIRE(source.config.has("wall_loops")); CHECK(source.config.opt_int("wall_loops") == 4);
        CHECK(source.config.opt_int("extruder") == 2);
        if (!fixture.source->source.input_file.empty()) CHECK(source.source.input_file == fixture.source->source.input_file);
        else if (fixture.source->is_model_part()) CHECK(fs::path(source.source.input_file).extension() == ".3mf");
        else CHECK(source.source.input_file.empty());
        CHECK(capture_face_extrusion_annotation_stamp(volume_named(object,fixture.neighbor->name)) ==
              capture_face_extrusion_annotation_stamp(*fixture.neighbor));
        CHECK(inspect_face_extrusion_annotations(source) == inspect_face_extrusion_annotations(*fixture.source));
    }
    CHECK_FALSE(fs::exists(directory));
}
} // namespace

TEST_CASE("Face model extraction preserves positive factory and replacement contracts", "[face_extrusion_model][negative_face_extrusion][negative_face_control]")
{
    for (Mode mode : {Mode::Faces,Mode::PointProfile}) {
        Fixture fixture(Role::MODEL_PART);
        const auto originals = fixture.saved(); const auto instance_matrix = fixture.instance->get_matrix();
        const auto object_config = config_values(fixture.object->config);
        auto patch = prepare_patch(fixture,0,2.); const auto preview = patch.world.its;
        auto *created = insert_face_extrusion_part(*fixture.object,std::move(patch.object),"Positive new part",mode,
                                                  Role::MODEL_PART,fixture.source->extruder_id());
        REQUIRE(created); CHECK(created->is_model_part()); CHECK(created->get_matrix().matrix().isIdentity(0.));
        REQUIRE(created->config.keys().size() == 1); CHECK(created->config.opt_int("extruder") == 2);
        CHECK(created->source.is_from_builtin_objects); CHECK(created->source.input_file.empty());
        CHECK_FALSE(inspect_face_extrusion_annotations(*created).any());
        require_world_geometry(world_mesh(*created,*fixture.instance),preview,1e-4);
        for (const auto &saved : originals) saved.unchanged();
        CHECK(same_matrix(fixture.instance->get_matrix(),instance_matrix)); CHECK(config_values(fixture.object->config) == object_config);
    }
    for (Mode mode : {Mode::Faces,Mode::PointProfile,Mode::CutAndStretch}) {
        Fixture fixture(Role::MODEL_PART); const SavedVolume source(*fixture.source); const auto neighbors = fixture.saved(false);
        const auto result = std::make_shared<const TriangleMesh>(ring());
        const auto hull = std::make_shared<const TriangleMesh>(result->convex_hull_3d());
        replace_face_extrusion_source(*fixture.source,mode,Role::MODEL_PART,result,hull,false);
        check_replacement_metadata(source); CHECK(fixture.source->get_mesh_shared_ptr() == result);
        for (const auto &saved : neighbors) saved.unchanged();
    }
    CHECK(face_extrusion_new_part_above_bed(Role::MODEL_PART,0.));
    CHECK_FALSE(face_extrusion_new_part_above_bed(Role::MODEL_PART,-0.1));
}

TEST_CASE("Face production role policy admits negative faces only", "[face_extrusion_model][negative_face_extrusion][negative_face_desired]")
{
    const std::array<Role,6> roles{{Role::MODEL_PART,Role::NEGATIVE_VOLUME,Role::PARAMETER_MODIFIER,
        Role::SUPPORT_BLOCKER,Role::SUPPORT_ENFORCER,Role::INVALID}};
    for (Mode mode : {Mode::Faces,Mode::PointProfile,Mode::CutAndStretch}) {
        CHECK_FALSE(face_extrusion_selection_roles_allowed(mode,{}));
        CHECK_FALSE(face_extrusion_selection_roles_allowed(mode,{Role::NEGATIVE_VOLUME,Role::NEGATIVE_VOLUME}));
        for (Role role : roles) {
            INFO("mode=" << int(mode) << ", role=" << int(role));
            const bool expected = role == Role::MODEL_PART || (mode == Mode::Faces && role == Role::NEGATIVE_VOLUME);
            CHECK(face_extrusion_source_role_allowed(mode,role) == expected);
            CHECK(face_extrusion_selection_roles_allowed(mode,{role}) == expected);
            CHECK(face_extrusion_selection_roles_allowed(mode,{Role::MODEL_PART,role}) == (role == Role::MODEL_PART));
            Fixture fixture(role);
            CHECK(face_extrusion_source_type_matches(mode,role,*fixture.source) == expected);
        }
    }
}

TEST_CASE("Face NewParts keep negative role and preserve every existing volume", "[face_extrusion_model][negative_face_extrusion][negative_face_desired]")
{
    for (bool multiple : {false,true}) for (int mirror_mask = 0; mirror_mask < 4; ++mirror_mask) for (double depth : {1.25,-1.25}) {
        INFO("multiple=" << multiple << ", mirrors=" << mirror_mask << ", depth=" << depth);
        Fixture fixture(Role::NEGATIVE_VOLUME,multiple,mirror_mask); paint_four_layers(*fixture.source);
        const auto originals = fixture.saved(); const auto instance_matrix = fixture.instance->get_matrix();
        const auto object_config = config_values(fixture.object->config); const size_t count = fixture.object->volumes.size();
        const size_t positive_count = fixture.object->parts_count();
        std::vector<Patch> prepared;
        for (size_t i = 0; i < fixture.patches.size(); ++i) prepared.push_back(prepare_patch(fixture,i,depth));
        CHECK(fixture.object->volumes.size() == count); // No part is inserted during preparation.
        std::set<size_t> created_ids;
        for (size_t i = 0; i < prepared.size(); ++i) {
            auto *created = insert_face_extrusion_part(*fixture.object,std::move(prepared[i].object),"Negative extrusion " + std::to_string(i),
                Mode::Faces,Role::NEGATIVE_VOLUME,fixture.source->extruder_id());
            REQUIRE(created); CHECK(created->is_negative_volume());
            CHECK(created->get_matrix().matrix().isIdentity(0.)); CHECK(created->id() != fixture.source->id());
            CHECK(created_ids.insert(created->id().id).second);
            CHECK(created->config.keys().empty()); CHECK(created->material_id().empty()); CHECK(created->part_guid().empty());
            CHECK_FALSE(inspect_face_extrusion_annotations(*created).any());
            CHECK(created->source.input_file.empty()); CHECK(created->source.object_idx == -1); CHECK(created->source.volume_idx == -1);
            CHECK(created->source.is_from_builtin_objects); require_solid(created->mesh().its,true);
            require_world_geometry(world_mesh(*created,*fixture.instance),prepared[i].world.its,1e-4);
        }
        fixture.object->invalidate_bounding_box();
        CHECK(fixture.object->volumes.size() == count+prepared.size()); CHECK(fixture.object->parts_count() == positive_count);
        for (const auto &saved : originals) saved.unchanged();
        CHECK(same_matrix(fixture.instance->get_matrix(),instance_matrix)); CHECK(config_values(fixture.object->config) == object_config);
    }
}

TEST_CASE("Face Boolean production replacement changes the negative cutter only", "[face_extrusion_model][negative_face_extrusion][negative_face_desired]")
{
    for (bool subtract : {false,true}) {
        INFO("subtract=" << subtract);
        Fixture fixture; paint_four_layers(*fixture.source);
        auto annotated = fixture.source->mesh(); annotated.its.properties[0] = FaceProperty{eSmallHole,9.};
        fixture.source->set_mesh(std::move(annotated));
        const SavedVolume source(*fixture.source); const auto neighbors = fixture.saved(false);
        const auto instance_matrix = fixture.instance->get_matrix(); const auto object_config = config_values(fixture.object->config);
        const auto result = boolean_result(fixture,subtract);
        const auto hull = std::make_shared<const TriangleMesh>(result->convex_hull_3d());
        source.unchanged(); for (const auto &saved : neighbors) saved.unchanged();
        // This is the real mutation seam with an explicitly prepared clear flag,
        // NOT a simulation/proof of GUI annotation consent or its Undo transaction.
        replace_face_extrusion_source(*fixture.source,Mode::Faces,Role::NEGATIVE_VOLUME,result,hull,true);
        check_replacement_metadata(source); CHECK(fixture.source->is_negative_volume());
        CHECK_FALSE(inspect_face_extrusion_annotations(*fixture.source).any());
        CHECK(fixture.source->mmuseg_extruders.empty()); CHECK(fixture.source->mmuseg_ts == 0);
        for (size_t i = 0; i < 4; ++i) CHECK((fixture.source->*layers[i]).timestamp() != source.timestamps[i]);
        require_solid(fixture.source->mesh().its,true);
        CHECK(its_volume(fixture.source->mesh().its) == Approx(subtract ? 300. : 1200.).epsilon(2e-5).margin(1e-3));
        CHECK(same_matrix(fixture.instance->get_matrix(),instance_matrix)); CHECK(config_values(fixture.object->config) == object_config);
        for (const auto &saved : neighbors) saved.unchanged();
    }
}

TEST_CASE("Face role changes and unsupported factories reject without model mutation", "[face_extrusion_model][negative_face_extrusion][negative_face_desired]")
{
    for (const auto &change : std::vector<std::pair<Role,Role>>{{Role::MODEL_PART,Role::NEGATIVE_VOLUME},
             {Role::NEGATIVE_VOLUME,Role::MODEL_PART},{Role::NEGATIVE_VOLUME,Role::PARAMETER_MODIFIER},
             {Role::NEGATIVE_VOLUME,Role::SUPPORT_ENFORCER}}) {
        Fixture fixture(change.first);
        const ObjectID id = fixture.source->id(); const auto mesh = fixture.source->get_mesh_shared_ptr();
        const auto matrix = fixture.source->get_matrix();
        fixture.source->set_type(change.second); // The actual native setter changes no ID/mesh/pose.
        REQUIRE(fixture.source->id() == id); REQUIRE(fixture.source->get_mesh_shared_ptr() == mesh);
        REQUIRE(same_matrix(fixture.source->get_matrix(),matrix));
        const auto after_external_change = fixture.saved();
        CHECK_FALSE(face_extrusion_source_type_matches(Mode::Faces,change.first,*fixture.source));
        REQUIRE_THROWS_AS(replace_face_extrusion_source(*fixture.source,Mode::Faces,change.first,mesh,
            fixture.source->get_convex_hull_shared_ptr(),false),std::invalid_argument);
        for (const auto &saved : after_external_change) saved.unchanged();
    }
    // The factory guards its DECLARED role. It has no live-source parameter;
    // stale-live-source rejection above belongs to the production match/seam.
    for (Role type : {Role::PARAMETER_MODIFIER,Role::SUPPORT_BLOCKER,Role::SUPPORT_ENFORCER,Role::INVALID}) {
        Fixture fixture; const auto before = fixture.saved(); const size_t count = fixture.object->volumes.size();
        REQUIRE_THROWS_AS(insert_face_extrusion_part(*fixture.object,TriangleMesh(ring()),"Unsupported",Mode::Faces,type,2),std::invalid_argument);
        CHECK(fixture.object->volumes.size() == count); for (const auto &saved : before) saved.unchanged();
    }
    for (Mode mode : {Mode::PointProfile,Mode::CutAndStretch}) {
        Fixture fixture; const auto before = fixture.saved();
        REQUIRE_THROWS_AS(insert_face_extrusion_part(*fixture.object,TriangleMesh(ring()),"Unsupported",mode,Role::NEGATIVE_VOLUME,2),std::invalid_argument);
        for (const auto &saved : before) saved.unchanged();
    }
    Fixture unchanged_negative;
    CHECK(face_extrusion_source_type_matches(Mode::Faces,Role::NEGATIVE_VOLUME,*unchanged_negative.source));
}

TEST_CASE("Face negative bed preflight matches native positive-only placement exactly", "[face_extrusion_model][negative_face_extrusion][negative_face_desired]")
{
    struct Case { size_t positives; double bottom, height; };
    const std::array<Case,10> cases{{{1,0.,1.},{1,2.,1.},{1,-0.5,1.},{1,-2.,1.},{1,1e-9,1.},
        {1,-0.0005,1.},{1,0.,0.02},{2,0.,1.},{2,-0.01,0.02},{2,2.,1.}}};
    for (const auto &test : cases) {
        INFO("positive_count=" << test.positives << ", bottom=" << test.bottom << ", height=" << test.height);
        Model model; auto *object = model.add_object();
        for (size_t i = 0; i < test.positives; ++i) {
            auto *carrier = object->add_volume(cube(2,2,test.height),Role::MODEL_PART,false);
            carrier->set_offset(Vec3d(4.*i,0.,test.bottom));
            carrier->set_convex_hull_shared_ptr(carrier->get_mesh_shared_ptr());
        }
        auto *source = object->add_volume(cube(2,2,2),Role::NEGATIVE_VOLUME,false);
        source->set_offset(Vec3d(0,0,-4));
        auto *instance = object->add_instance(); instance->set_offset(Vec3d(75,90,0));
        object->invalidate_bounding_box();
        std::vector<SavedVolume> before; for (const auto *v : object->volumes) before.emplace_back(*v);
        const auto instance_before = instance->get_matrix();
        const auto context = capture_face_extrusion_placement(*object,*instance,*source,instance->get_matrix()*source->get_matrix());
        CHECK(context.source_type == Role::NEGATIVE_VOLUME); CHECK(context.model_parts == test.positives);
        REQUIRE(context.volumes.size() == object->volumes.size());
        for (size_t i = 0; i < context.volumes.size(); ++i) {
            CHECK(context.volumes[i].id == object->volumes[i]->id());
            CHECK(context.volumes[i].type == object->volumes[i]->type());
            CHECK(context.volumes[i].identity == object->volumes[i]->get_mesh_shared_ptr());
        }
        TriangleMesh result = cube(8,8,20); result.translate(-20,-20,-10);
        const auto placement = face_extrusion_replacement_placement(context,result,result);
        CHECK(placement.min_z == object->get_instance_min_z(0));
        CHECK(placement.max_z == object->get_instance_max_z(0));
        const auto expected_bounds = object->instance_bounding_box(0);
        CHECK(placement.printable_bounds.min == expected_bounds.min);
        CHECK(placement.printable_bounds.max == expected_bounds.max);
        // Independent native oracle, on a PRIVATE Model copy: install the known
        // negative candidate through native APIs, not a test reimplementation of
        // ensure_on_bed. This does not assert the production factory's role for it.
        Model replacement_oracle = model;
        auto *replacement_object = replacement_oracle.objects.front();
        auto *replacement_source = replacement_object->volumes.back();
        replacement_source->set_mesh(result);
        replacement_source->set_convex_hull_shared_ptr(replacement_source->get_mesh_shared_ptr());
        replacement_object->invalidate_bounding_box(); replacement_object->ensure_on_bed(true);
        const bool replacement_stays = same_matrix(replacement_object->instances.front()->get_matrix(),instance_before);
        CHECK(face_extrusion_bed_safe(context.model_parts,placement.min_z,placement.max_z) == replacement_stays);
        Model insertion_oracle = model; auto *insertion_object = insertion_oracle.objects.front();
        insertion_object->add_volume(TriangleMesh(result),Role::NEGATIVE_VOLUME,false);
        insertion_object->invalidate_bounding_box(); insertion_object->ensure_on_bed(true);
        const bool insertion_stays = same_matrix(insertion_object->instances.front()->get_matrix(),instance_before);
        CHECK(face_extrusion_new_part_above_bed(Role::NEGATIVE_VOLUME,-10.));
        CHECK(face_extrusion_new_parts_bed_safe(Role::NEGATIVE_VOLUME,*object,*instance,
            std::max(object->get_max_z(),10.)) == insertion_stays);
        CHECK(same_matrix(instance->get_matrix(),instance_before)); for (const auto &saved : before) saved.unchanged();
    }
    CHECK_FALSE(face_extrusion_bed_safe(0,0.,1.));
    CHECK_FALSE(face_extrusion_new_part_above_bed(Role::MODEL_PART,-10.));
}

TEST_CASE("Face negative NewParts retain roles and geometry through native BBS 3MF", "[face_extrusion_model][negative_face_extrusion][negative_face_desired][bbs_3mf]")
{
    for (double depth : {1.25,-1.25}) {
        Fixture fixture(Role::NEGATIVE_VOLUME,true,1); paint_four_layers(*fixture.source);
        for (size_t i = 0; i < fixture.patches.size(); ++i) {
            auto patch = prepare_patch(fixture,i,depth);
            auto *created = insert_face_extrusion_part(*fixture.object,std::move(patch.object),"Saved negative extrusion " + std::to_string(i),
                Mode::Faces,Role::NEGATIVE_VOLUME,fixture.source->extruder_id());
            REQUIRE(created->is_negative_volume()); CHECK(created->config.keys().empty());
        }
        fixture.object->invalidate_bounding_box(); require_bbs_roundtrip(fixture);
    }
}

TEST_CASE("Face negative Boolean results retain roles through native BBS 3MF", "[face_extrusion_model][negative_face_extrusion][negative_face_desired][bbs_3mf]")
{
    for (bool subtract : {false,true}) {
        Fixture fixture; paint_four_layers(*fixture.source);
        const auto result = boolean_result(fixture,subtract);
        const auto hull = std::make_shared<const TriangleMesh>(result->convex_hull_3d());
        replace_face_extrusion_source(*fixture.source,Mode::Faces,Role::NEGATIVE_VOLUME,result,hull,true);
        REQUIRE(fixture.source->is_negative_volume()); require_bbs_roundtrip(fixture);
    }
}
