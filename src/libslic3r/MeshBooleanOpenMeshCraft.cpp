#include "Exception.hpp"
#include "MeshBoolean.hpp"
#include "MeshDiagnostics.hpp"
#include "libslic3r/TriangleMesh.hpp"
#include "libslic3r/TryCatchSignal.hpp"

#include <OpenMeshCraft/OpenMeshCraftBoolean.h>

#include "boost/log/trivial.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace Slic3r {
namespace MeshBoolean {
namespace openmeshcraft {

namespace {

enum class Operation { Union, Intersection, Difference };
struct Canceled {};

void check_canceled(const BooleanCancelCB &cancel_cb)
{
    if (cancel_cb && cancel_cb())
        throw Canceled{};
}

struct SolidInfo
{
    bool empty = true;
    double volume = 0.0;
    double area = 0.0;
    Vec3d min = Vec3d::Zero();
    Vec3d max = Vec3d::Zero();
};

// Validate indices and arithmetic BEFORE constructing TriangleMesh or invoking
// diagnostics/CGAL, all of which assume valid indexed triangles. No welding,
// snapping, orientation repair or geometry-changing tolerance is applied here.
SolidInfo validate_solid(const indexed_triangle_set &mesh, const BooleanCancelCB &cancel_cb)
{
    check_canceled(cancel_cb);
    if (mesh.vertices.size() > size_t(std::numeric_limits<int>::max()) ||
        mesh.indices.size() > size_t(std::numeric_limits<int>::max()))
        throw RuntimeError("Boolean mesh exceeds the indexed mesh size limit.");
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
        if ((i & 1023) == 0) check_canceled(cancel_cb);
        if (!mesh.vertices[i].allFinite())
            throw RuntimeError("Boolean mesh has a non-finite vertex.");
    }

    SolidInfo info;
    if (mesh.indices.empty())
        return info; // The empty set is a valid operand/result.
    info.empty = false;
    info.min.setConstant(std::numeric_limits<double>::infinity());
    info.max.setConstant(-std::numeric_limits<double>::infinity());
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        if ((i & 1023) == 0) check_canceled(cancel_cb);
        const Vec3i &face = mesh.indices[i];
        for (int index : {face[0], face[1], face[2]}) {
            if (index < 0 || size_t(index) >= mesh.vertices.size())
                throw RuntimeError("Boolean mesh has an invalid triangle index.");
            const Vec3d vertex = mesh.vertices[index].cast<double>();
            info.min = info.min.cwiseMin(vertex);
            info.max = info.max.cwiseMax(vertex);
        }
        if (face[0] == face[1] || face[1] == face[2] || face[2] == face[0])
            throw RuntimeError("Boolean mesh has a degenerate triangle.");
        const Vec3d a = mesh.vertices[face[0]].cast<double>();
        const Vec3d b = mesh.vertices[face[1]].cast<double>();
        const Vec3d c = mesh.vertices[face[2]].cast<double>();
        const double area = (b - a).cross(c - a).norm() * 0.5;
        const Vec3f float_cross = (mesh.vertices[face[1]] - mesh.vertices[face[0]]).cross(
                                  mesh.vertices[face[2]] - mesh.vertices[face[0]]);
        if (!std::isfinite(area) || area <= 0.0 || !std::isfinite(float_cross.squaredNorm()) ||
            float_cross.squaredNorm() <= 0.f)
            throw RuntimeError("Boolean mesh has a degenerate or unrepresentable triangle.");
        info.area += area;
    }

    // Translation-invariant double accumulation, not TriangleMesh's cached float
    // volume. Compensated summation also preserves small inner-cavity volumes.
    const Vec3d origin = 0.5 * (info.min + info.max);
    double compensation = 0.0;
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        if ((i & 1023) == 0) check_canceled(cancel_cb);
        const Vec3i &face = mesh.indices[i];
        const Vec3d a = mesh.vertices[face[0]].cast<double>() - origin;
        const Vec3d b = mesh.vertices[face[1]].cast<double>() - origin;
        const Vec3d c = mesh.vertices[face[2]].cast<double>() - origin;
        const double term = a.dot(b.cross(c)) / 6.0 - compensation;
        const double sum = info.volume + term;
        compensation = (sum - info.volume) - term;
        info.volume = sum;
    }
    if (!std::isfinite(info.area) || !std::isfinite(info.volume) || info.volume <= 0.0 ||
        info.volume > std::numeric_limits<float>::max())
        throw RuntimeError("Boolean mesh does not have a finite positive solid volume.");
    // TriangleMesh still computes cached statistics in float with its first
    // vertex as origin, which may even be an unreferenced far-away vertex.
    const float stored_volume = its_volume(mesh);
    if (!std::isfinite(stored_volume) || stored_volume <= 0.f)
        throw RuntimeError("Boolean mesh volume is not representable by TriangleMesh.");

    const MeshDiagnosticStats stats = its_mesh_diagnostics(mesh);
    check_canceled(cancel_cb);
    if (stats.open_edges || stats.non_manifold_edges || stats.non_manifold_vertices ||
        stats.same_direction_edges || stats.has_reversed_faces)
        throw RuntimeError("Boolean mesh is not closed, manifold and consistently oriented.");
    auto cgal_mesh = cgal::triangle_mesh_to_cgal(mesh);
    check_canceled(cancel_cb);
    if (cgal::does_self_intersect(*cgal_mesh))
        throw RuntimeError("Boolean mesh has self-intersections.");
    check_canceled(cancel_cb);
    if (!cgal::does_bound_a_volume(*cgal_mesh))
        throw RuntimeError("Boolean mesh shells do not bound a consistently oriented volume.");
    check_canceled(cancel_cb);
    return info;
}

void validate_result(const SolidInfo &result, const SolidInfo &a, const SolidInfo &b, Operation operation)
{
    const double magnitude = std::max({a.min.cwiseAbs().maxCoeff(), a.max.cwiseAbs().maxCoeff(),
                                       b.min.cwiseAbs().maxCoeff(), b.max.cwiseAbs().maxCoeff(),
                                       result.min.cwiseAbs().maxCoeff(), result.max.cwiseAbs().maxCoeff()});
    // These bounds accommodate float storage rounding only. They never modify
    // coordinates, create overlap, or permit a topological defect.
    const double length_error = 8.0 * std::numeric_limits<float>::epsilon() * magnitude;
    const double volume_error = length_error * (a.area + b.area + result.area) +
                                32.0 * std::numeric_limits<double>::epsilon() * (a.volume + b.volume);
    double lower = 0.0;
    double upper = 0.0;
    if (operation == Operation::Union) {
        lower = std::max(a.volume, b.volume);
        upper = a.volume + b.volume;
        if (result.empty && (!a.empty || !b.empty))
            throw RuntimeError("Boolean union unexpectedly removed all geometry.");
    } else if (operation == Operation::Intersection) {
        upper = std::min(a.volume, b.volume);
    } else {
        lower = std::max(0.0, a.volume - b.volume);
        upper = a.volume;
    }
    if (result.volume < lower - volume_error || result.volume > upper + volume_error)
        throw RuntimeError("Boolean result violates solid-volume bounds.");
    if (result.empty)
        return;

    auto inside = [length_error](const SolidInfo &inner, const SolidInfo &outer) {
        return inner.empty || (!outer.empty &&
               (inner.min.array() >= outer.min.array() - length_error).all() &&
               (inner.max.array() <= outer.max.array() + length_error).all());
    };
    if (operation == Operation::Union) {
        SolidInfo combined = a.empty ? b : a;
        if (!a.empty && !b.empty) {
            combined.min = a.min.cwiseMin(b.min);
            combined.max = a.max.cwiseMax(b.max);
        }
        if (!inside(a, result) || !inside(b, result) || !inside(result, combined))
            throw RuntimeError("Boolean union violates operand bounds.");
    } else if (!inside(result, a) || (operation == Operation::Intersection && !inside(result, b))) {
        throw RuntimeError("Boolean result extends outside its operands.");
    }
}

struct OMCInput
{
    std::vector<BambuOMCPoint> points;
    std::vector<BambuOMCTriangle> triangles;
    BambuOMCMesh view() const { return {points.data(), points.size(), triangles.data(), triangles.size()}; }
};
struct OMCOutput
{
    std::vector<BambuOMCPoint> points;
    std::vector<BambuOMCTriangle> triangles;
};

OMCInput triangle_mesh_to_openmeshcraft(const indexed_triangle_set &mesh, const BooleanCancelCB &cancel_cb)
{
    OMCInput input;
    input.points.reserve(mesh.vertices.size());
    input.triangles.reserve(mesh.indices.size());
    for (size_t i = 0; i < mesh.vertices.size(); ++i) {
        if ((i & 1023) == 0) check_canceled(cancel_cb);
        const Vec3f &v = mesh.vertices[i];
        input.points.push_back({double(v.x()), double(v.y()), double(v.z())});
    }
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        if ((i & 1023) == 0) check_canceled(cancel_cb);
        const Vec3i &f = mesh.indices[i];
        input.triangles.push_back({{size_t(f[0]), size_t(f[1]), size_t(f[2])}});
    }
    return input;
}

void BAMBU_OMC_CALL receive_openmeshcraft_result(const BambuOMCPoint *points, size_t point_count,
                                                const BambuOMCTriangle *triangles, size_t triangle_count, void *user_data)
{
    if ((point_count && !points) || (triangle_count && !triangles) ||
        point_count > size_t(std::numeric_limits<int>::max()) || triangle_count > size_t(std::numeric_limits<int>::max()))
        throw RuntimeError("OpenMeshCraft returned invalid mesh buffers.");
    auto &output = *static_cast<OMCOutput *>(user_data);
    output.points.clear();
    output.triangles.clear();
    if (point_count) output.points.assign(points, points + point_count);
    if (triangle_count) output.triangles.assign(triangles, triangles + triangle_count);
}

indexed_triangle_set openmeshcraft_to_mesh(const OMCOutput &output, const BooleanCancelCB &cancel_cb)
{
    indexed_triangle_set mesh;
    mesh.vertices.reserve(output.points.size());
    mesh.indices.reserve(output.triangles.size());
    for (size_t i = 0; i < output.points.size(); ++i) {
        if ((i & 1023) == 0) check_canceled(cancel_cb);
        const auto &p = output.points[i];
        for (double value : {p.x, p.y, p.z})
            if (!std::isfinite(value) || std::abs(value) > std::numeric_limits<float>::max())
                throw RuntimeError("OpenMeshCraft returned an unrepresentable vertex.");
        mesh.vertices.emplace_back(float(p.x), float(p.y), float(p.z));
    }
    for (size_t i = 0; i < output.triangles.size(); ++i) {
        if ((i & 1023) == 0) check_canceled(cancel_cb);
        const auto &f = output.triangles[i];
        for (size_t index : f.vertices)
            if (index >= output.points.size() || index > size_t(std::numeric_limits<int>::max()))
                throw RuntimeError("OpenMeshCraft returned an invalid triangle index.");
        mesh.indices.emplace_back(int(f.vertices[0]), int(f.vertices[1]), int(f.vertices[2]));
    }
    return mesh;
}

} // namespace

void make_boolean(const TriangleMesh &src_mesh, const TriangleMesh &cut_mesh, std::vector<TriangleMesh> &dst_mesh,
                  const std::string &boolean_opts, const BooleanCancelCB &cancel_cb,
                  const BooleanProgressCB &progress_cb, const BooleanFailedCB &failed_cb)
{
    try {
        check_canceled(cancel_cb);
        if (progress_cb) progress_cb(5.0f);
        const Operation operation = boolean_opts == "UNION" ? Operation::Union :
                                    boolean_opts == "INTERSECTION" ? Operation::Intersection : Operation::Difference;
        if (boolean_opts != "UNION" && boolean_opts != "INTERSECTION" &&
            boolean_opts != "A_NOT_B" && boolean_opts != "B_NOT_A")
            throw RuntimeError("Unsupported mesh boolean operation: " + boolean_opts);
        const auto &a = boolean_opts == "B_NOT_A" ? cut_mesh.its : src_mesh.its;
        const auto &b = boolean_opts == "B_NOT_A" ? src_mesh.its : cut_mesh.its;
        const SolidInfo a_info = validate_solid(a, cancel_cb);
        const SolidInfo b_info = validate_solid(b, cancel_cb);
        if (progress_cb) progress_cb(15.0f);
        check_canceled(cancel_cb);
        indexed_triangle_set candidate;
        bool candidate_validated = false;

        if (a_info.empty || b_info.empty) {
            // Exact set identities; don't send empty soups to either backend.
            if (operation == Operation::Union) candidate = a_info.empty ? b : a;
            else if (operation == Operation::Difference) candidate = a;
        } else {
            BOOST_LOG_TRIVIAL(info) << "Mesh boolean backend: OpenMeshCraft, op='" << boolean_opts << "'";
            OMCInput src_input = triangle_mesh_to_openmeshcraft(a, cancel_cb);
            OMCInput cut_input = triangle_mesh_to_openmeshcraft(b, cancel_cb);
            const BambuOMCMesh src_view = src_input.view();
            const BambuOMCMesh cut_view = cut_input.view();
            OMCOutput output;
            bool boolean_ok = false;
            bool hw_fault = false;
            std::string omc_error;
            // Catch C++ exceptions INSIDE the signal/SEH lambda. Never unwind a
            // C++ exception through MSVC __try/__except. TBB policy is unchanged.
            try_catch_signal({SIGSEGV, SIGFPE}, [&]() {
                try {
                    const BambuOMCOperation op = operation == Operation::Union ? BAMBU_OMC_UNION :
                                                operation == Operation::Intersection ? BAMBU_OMC_INTERSECTION : BAMBU_OMC_SUBTRACTION;
                    std::array<char, 1024> message{};
                    boolean_ok = bambu_omc_boolean(&src_view, &cut_view, op, receive_openmeshcraft_result,
                                                   &output, message.data(), message.size());
                    if (!boolean_ok) omc_error = message.data();
                } catch (const std::exception &e) {
                    omc_error = e.what();
                }
            }, [&]() { hw_fault = true; });
            if (progress_cb) progress_cb(60.0f);
            check_canceled(cancel_cb);
            if (hw_fault) omc_error = "hardware fault";
            else if (!boolean_ok && omc_error.empty()) omc_error = "operation did not complete";
            if (boolean_ok && !hw_fault && omc_error.empty()) {
                try {
                    candidate = openmeshcraft_to_mesh(output, cancel_cb);
                    validate_result(validate_solid(candidate, cancel_cb), a_info, b_info, operation);
                    candidate_validated = true;
                    // Zero is mathematically valid for subtraction/intersection,
                    // but volume bounds alone cannot confirm an empty answer.
                    if (candidate.indices.empty())
                        omc_error = "confirming empty solid result with CGAL";
                } catch (const Canceled &) { throw; }
                catch (const std::exception &e) { omc_error = e.what(); }
            }
            if (!omc_error.empty()) {
                BOOST_LOG_TRIVIAL(warning) << "OpenMeshCraft result rejected (" << omc_error
                                           << "); trying CGAL for '" << boolean_opts << "'.";
                if (progress_cb) progress_cb(70.0f);
                check_canceled(cancel_cb);
                // Operate only on a local copy. Neither backend may publish an
                // intermediate result or mutate either caller-owned operand.
                candidate = a;
                if (operation == Operation::Union) cgal::plus(candidate, b);
                else if (operation == Operation::Intersection) cgal::intersect(candidate, b);
                else cgal::minus(candidate, b);
                check_canceled(cancel_cb);
                if (progress_cb) progress_cb(85.0f);
                validate_result(validate_solid(candidate, cancel_cb), a_info, b_info, operation);
                candidate_validated = true;
            }
        }

        if (!candidate_validated)
            validate_result(validate_solid(candidate, cancel_cb), a_info, b_info, operation);
        if (progress_cb) progress_cb(90.0f);
        check_canceled(cancel_cb);
        // Keep inner cavity shells negative. Flipping every disconnected shell
        // to positive volume would turn a valid subtraction cavity into solid.
        const bool append = !candidate.indices.empty();
        if (append) {
            TriangleMesh final_mesh{std::move(candidate)};
            check_canceled(cancel_cb);
            dst_mesh.push_back(std::move(final_mesh));
        }
        try {
            if (progress_cb) progress_cb(100.0f);
            check_canceled(cancel_cb);
        } catch (...) {
            // Even a terminal callback exception/cancellation must not leave a
            // partly committed answer in the caller's existing output vector.
            if (append) dst_mesh.pop_back();
            throw;
        }
    } catch (const Canceled &) {
        // Cancellation is not a backend error and never publishes a result.
    } catch (const std::exception &e) {
        if (cancel_cb && cancel_cb()) return;
        BOOST_LOG_TRIVIAL(error) << "Mesh boolean failed: " << e.what();
        if (failed_cb) failed_cb();
    }
}

} // namespace openmeshcraft
} // namespace MeshBoolean
} // namespace Slic3r
