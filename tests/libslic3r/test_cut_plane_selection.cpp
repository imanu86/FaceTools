// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include <catch2/catch.hpp>
#include "libslic3r/CutPlaneSelection.hpp"
#include <Eigen/LU>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>

using namespace Slic3r;

namespace {
void invalid(const CutPlaneSelectionResult &result)
{
    CHECK_FALSE(result);
    CHECK(result.side == 0);
    CHECK(result.moving_normal == Vec3d::Zero());
    CHECK_FALSE(result.message.empty());
}
void selected(const CutPlaneSelectionResult &result, int side, const Vec3d &expected)
{
    REQUIRE(result);
    CHECK(result.side == side);
    CHECK(result.message.empty());
    CHECK(result.moving_normal.norm() == Approx(1.0).margin(1e-14));
    CHECK((result.moving_normal - expected).norm() < 1e-13);
}
std::array<uint64_t, 3> bits(const Vec3d &value)
{
    std::array<uint64_t, 3> out;
    static_assert(sizeof(double) == sizeof(uint64_t), "double representation required");
    for (size_t i = 0; i < 3; ++i) std::memcpy(&out[i], &value[i], sizeof(double));
    return out;
}
}

TEST_CASE("Cut guide L top and underside select material relative to an independent plane", "[cut_plane_selection]")
{
    const Vec3d origin(0, 0, 20), normal = Vec3d::UnitZ();
    selected(select_cut_moving_side({{Vec3d(18, 5, 24), Vec3d::UnitZ()}}, origin, normal), 1, normal);
    // Underside of the upper arm, exactly on the cut: material is ABOVE -Z.
    selected(select_cut_moving_side({{Vec3d(18, 5, 20), -Vec3d::UnitZ()}}, origin, normal), 1, normal);
    selected(select_cut_moving_side({{Vec3d(18, 5, 20), Vec3d::UnitZ()}}, origin, normal), -1, -normal);
    selected(select_cut_moving_side({{Vec3d(18, 5, 16), -Vec3d::UnitZ()}}, origin, normal), -1, -normal);
    // An off-plane point chooses its location, NEVER -outward unconditionally.
    selected(select_cut_moving_side({{Vec3d(18, 5, 24), Vec3d::UnitZ()}}, origin, -normal), -1, normal);
}

TEST_CASE("Cut guides with different outward normals select one rigid side without averaging", "[cut_plane_selection]")
{
    const Vec3d origin(0, 0, 20);
    const std::vector<CutPlaneFaceGuide> faces = {
        {Vec3d(25, 2, 22), Vec3d::UnitX()}, // first face is lateral, not Z
        {Vec3d(10, 8, 24), Vec3d::UnitZ()},
        {Vec3d(18, -2, 23), -Vec3d::UnitY()},
    };
    std::array<int, 3> order{{0, 1, 2}};
    do {
        std::vector<CutPlaneFaceGuide> guides;
        for (int i : order) guides.push_back(faces[i]);
        selected(select_cut_moving_side(guides, origin, Vec3d::UnitZ()), 1, Vec3d::UnitZ());
        selected(select_cut_moving_side(guides, origin, -Vec3d::UnitZ()), -1, Vec3d::UnitZ());
    } while (std::next_permutation(order.begin(), order.end()));
}

TEST_CASE("Cut tube outer inner and cap guides do not confuse outward direction with moving direction", "[cut_plane_selection]")
{
    // Pure guides only: the parent's workflow tests cover the actual tube kernel.
    const std::vector<CutPlaneFaceGuide> tube = {
        {Vec3d(5, 0, 12), Vec3d::UnitX()},
        {Vec3d(3, 0, 13), -Vec3d::UnitX()},
        {Vec3d(0, 5, 14), Vec3d::UnitY()},
        {Vec3d(0, 3, 15), -Vec3d::UnitY()},
    }; // mean outward normal is zero; all points still select +Z
    selected(select_cut_moving_side(tube, Vec3d(0, 0, 10), Vec3d::UnitZ()), 1, Vec3d::UnitZ());
    selected(select_cut_moving_side({{Vec3d(4, 0, 25), Vec3d::UnitZ()}}, Vec3d(0, 0, 25), Vec3d::UnitZ()), -1, -Vec3d::UnitZ());
    selected(select_cut_moving_side({{Vec3d(4, 0, 0), -Vec3d::UnitZ()}}, Vec3d::Zero(), Vec3d::UnitZ()), 1, Vec3d::UnitZ());
}

TEST_CASE("Cut mixed guides reject even with a large majority and never leak partial output", "[cut_plane_selection]")
{
    std::vector<CutPlaneFaceGuide> guides(32, {Vec3d(1, 1, 1000), Vec3d::UnitX()});
    guides.push_back({Vec3d(1, 1, -1), Vec3d::UnitZ()});
    invalid(select_cut_moving_side(guides, Vec3d::Zero(), Vec3d::UnitZ()));
    std::reverse(guides.begin(), guides.end());
    invalid(select_cut_moving_side(guides, Vec3d::Zero(), Vec3d::UnitZ()));
    invalid(select_cut_moving_side({{Vec3d(1, 0, 2), Vec3d::UnitX()}, {Vec3d::Zero(), Vec3d::UnitZ()}}, Vec3d::Zero(), Vec3d::UnitZ()));
    invalid(select_cut_moving_side({{Vec3d::Zero(), -Vec3d::UnitZ()}, {Vec3d(1, 0, 0), Vec3d::UnitZ()}}, Vec3d::Zero(), Vec3d::UnitZ()));
    // A definite guide must not override another ambiguous on-plane guide.
    invalid(select_cut_moving_side({{Vec3d(0, 0, 2), Vec3d::UnitZ()}, {Vec3d::Zero(), Vec3d::UnitX()}}, Vec3d::Zero(), Vec3d::UnitZ()));
}

TEST_CASE("Cut near-plane band never becomes the coplanar material-side fallback", "[cut_plane_selection]")
{
    for (double distance : {8e-5, -8e-5, 1e-4, -1e-4, 1e-8, -1e-8, 1e-12, -1e-12,
                            std::numeric_limits<double>::denorm_min(), -std::numeric_limits<double>::denorm_min()}) {
        for (double outward_sign : {-1.0, 1.0}) {
            INFO("distance=" << distance << " outward=" << outward_sign);
            invalid(select_cut_moving_side({{Vec3d(12, 18, distance), outward_sign * Vec3d::UnitZ()}}, Vec3d::Zero(), Vec3d::UnitZ()));
        }
    }
    selected(select_cut_moving_side({{Vec3d(0, 0, 2e-4), Vec3d::UnitZ()}}, Vec3d::Zero(), Vec3d::UnitZ()), 1, Vec3d::UnitZ());
    selected(select_cut_moving_side({{Vec3d(0, 0, -2e-4), -Vec3d::UnitZ()}}, Vec3d::Zero(), Vec3d::UnitZ()), -1, -Vec3d::UnitZ());
    // Even an exactly cancelling dot at huge tangent coordinates is rejected if
    // the arithmetic uncertainty would intrude into the geometric tolerance band.
    invalid(select_cut_moving_side({{Vec3d(1e14, -1e14, 0), Vec3d(-1, -1, 0)}}, Vec3d::Zero(), Vec3d(1, 1, 0)));
    const Vec3d far_origin(1e8, -1e8, 1e8);
    invalid(select_cut_moving_side({{far_origin + Vec3d(3, 4, 8e-5), -Vec3d::UnitZ()}}, far_origin, Vec3d::UnitZ()));
    selected(select_cut_moving_side({{far_origin + Vec3d(3, 4, 0), -Vec3d::UnitZ()}}, far_origin, Vec3d::UnitZ()), 1, Vec3d::UnitZ());
}

TEST_CASE("Cut coplanar nonparallel guides reject rather than using a rounded dot-product parallel test", "[cut_plane_selection]")
{
    invalid(select_cut_moving_side({{Vec3d::Zero(), Vec3d::UnitX()}}, Vec3d::Zero(), Vec3d::UnitZ()));
    invalid(select_cut_moving_side({{Vec3d::Zero(), Vec3d(1, 0, 1)}}, Vec3d::Zero(), Vec3d::UnitZ()));
    // 1 - abs(dot(normalized(...), Z)) may round to zero for this real tilt.
    invalid(select_cut_moving_side({{Vec3d::Zero(), Vec3d(1e-8, 0, 1)}}, Vec3d::Zero(), Vec3d::UnitZ()));
    selected(select_cut_moving_side({{Vec3d(3, 4, 0), Vec3d(0, 0, -7)}}, Vec3d::Zero(), Vec3d(0, 0, 3)), 1, Vec3d::UnitZ());
}

TEST_CASE("Cut side and moving normal are equivariant under world frames and normal reversal", "[cut_plane_selection]")
{
    for (int sample = 0; sample < 48; ++sample) {
        Transform3d frame = Transform3d::Identity();
        frame.linear() = Eigen::AngleAxisd(0.137 * sample, Vec3d(1, 3, 2).normalized()).toRotationMatrix() *
                         Vec3d(sample & 1 ? -2.0 : 2.0, 0.7, sample & 2 ? -3.0 : 3.0).asDiagonal();
        frame.translation() = Vec3d(47.25, -101.75, 33.5);
        const Matrix3d normal_matrix = frame.linear().inverse().transpose();
        const Vec3d plane_origin = frame * Vec3d(0, 0, 20);
        const Vec3d plane_normal = normal_matrix * Vec3d::UnitZ();
        const Vec3d expected = plane_normal.normalized();
        INFO("world frame=" << sample);
        const std::vector<CutPlaneFaceGuide> far = {
            {frame * Vec3d(2, 3, 24), normal_matrix * Vec3d::UnitX()},
            {frame * Vec3d(4, 2, 23), normal_matrix * Vec3d::UnitY()},
        };
        const std::vector<CutPlaneFaceGuide> underside = {
            {frame * Vec3d(3, 4, 20), normal_matrix * -Vec3d::UnitZ()},
        };
        selected(select_cut_moving_side(far, plane_origin, plane_normal), 1, expected);
        selected(select_cut_moving_side(far, plane_origin, -plane_normal), -1, expected);
        selected(select_cut_moving_side(underside, plane_origin, plane_normal), 1, expected);
        selected(select_cut_moving_side(underside, plane_origin, -plane_normal), -1, expected);
    }
}

TEST_CASE("Cut normals normalize robustly without overflow or underflow", "[cut_plane_selection]")
{
    for (double scale : {1e300, 1e-300, std::numeric_limits<double>::denorm_min()}) {
        selected(select_cut_moving_side({{Vec3d(2, 3, 4), Vec3d(0, 0, -scale)}}, Vec3d::Zero(), Vec3d(0, 0, scale)), 1, Vec3d::UnitZ());
        selected(select_cut_moving_side({{Vec3d::Zero(), Vec3d(0, 0, -scale)}}, Vec3d::Zero(), Vec3d(0, 0, scale)), 1, Vec3d::UnitZ());
        selected(select_cut_moving_side({{Vec3d::Zero(), Vec3d(0, 0, -scale)}}, Vec3d::Zero(), Vec3d(0, 0, -scale)), -1, Vec3d::UnitZ());
    }
    selected(select_cut_moving_side({{Vec3d::Zero(), Vec3d(-1e300, -1e300, -1e300)}}, Vec3d::Zero(), Vec3d(1e300, 1e300, 1e300)),
             1, Vec3d(1, 1, 1).normalized());
}

TEST_CASE("Cut invalid inputs and guide count limits fail closed with immutable inputs", "[cut_plane_selection]")
{
    const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
    const std::vector<CutPlaneFaceGuide> good{{Vec3d(-0.0, 2, 3), Vec3d(0, 0, 1)}};
    invalid(select_cut_moving_side({}, Vec3d::Zero(), Vec3d::UnitZ()));
    for (const Vec3d &normal : {Vec3d(0, 0, 0), Vec3d(nan, 0, 1), Vec3d(0, inf, 1)})
        invalid(select_cut_moving_side(good, Vec3d::Zero(), normal));
    for (double tolerance : {0.0, -1.0, nan, inf})
        invalid(select_cut_moving_side(good, Vec3d::Zero(), Vec3d::UnitZ(), tolerance));
    invalid(select_cut_moving_side(good, Vec3d(0, nan, 0), Vec3d::UnitZ()));
    invalid(select_cut_moving_side(good, Vec3d(inf, 0, 0), Vec3d::UnitZ()));
    for (const CutPlaneFaceGuide &bad : {
        CutPlaneFaceGuide{Vec3d(0, 0, nan), Vec3d::UnitZ()},
        CutPlaneFaceGuide{Vec3d(inf, 0, 1), Vec3d::UnitZ()},
        CutPlaneFaceGuide{Vec3d(0, 0, 1), Vec3d::Zero()},
        CutPlaneFaceGuide{Vec3d(0, 0, 1), Vec3d(nan, 0, 1)},
        CutPlaneFaceGuide{Vec3d(0, 0, 1), Vec3d(0, inf, 1)},
    }) {
        std::vector<CutPlaneFaceGuide> guides = good;
        guides.push_back(bad); // failure after an already-valid guide must not leak its side
        const auto before = guides;
        const Vec3d origin(-0.0, 0.0, 0.0), normal(0.0, -0.0, 1.0);
        const auto origin_bits = bits(origin), normal_bits = bits(normal);
        invalid(select_cut_moving_side(guides, origin, normal));
        for (size_t i = 0; i < guides.size(); ++i) {
            CHECK(bits(guides[i].point_world) == bits(before[i].point_world));
            CHECK(bits(guides[i].outward_normal_world) == bits(before[i].outward_normal_world));
        }
        CHECK(bits(origin) == origin_bits);
        CHECK(bits(normal) == normal_bits);
    }
    auto maximum = std::vector<CutPlaneFaceGuide>(4096, good.front());
    selected(select_cut_moving_side(maximum, Vec3d::Zero(), Vec3d::UnitZ()), 1, Vec3d::UnitZ());
    maximum.push_back(good.front());
    invalid(select_cut_moving_side(maximum, Vec3d::Zero(), Vec3d::UnitZ()));
    const auto good_before = good;
    selected(select_cut_moving_side(good, Vec3d::Zero(), Vec3d::UnitZ()), 1, Vec3d::UnitZ());
    CHECK(bits(good[0].point_world) == bits(good_before[0].point_world));
    CHECK(bits(good[0].outward_normal_world) == bits(good_before[0].outward_normal_world));
    const double maximum_coordinate = std::numeric_limits<double>::max();
    invalid(select_cut_moving_side({{Vec3d(maximum_coordinate, 0, 0), Vec3d::UnitX()}}, Vec3d(-maximum_coordinate, 0, 0), Vec3d::UnitX()));
    invalid(select_cut_moving_side({{Vec3d::Constant(maximum_coordinate), Vec3d::UnitX()}}, Vec3d::Zero(), Vec3d(1, 1, 1)));
}
