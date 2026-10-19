// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

// Production encoder and POSIX writer tests. Synthetic data only; never user STL.
// IO tests are HIDDEN by default and must be explicitly selected later by parent.
#include <catch2/catch.hpp>
#include "libslic3r/CheckedStlExport.hpp"
#include "libslic3r/TriangleMesh.hpp"

#include <boost/filesystem.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace Slic3r;
namespace {
class RecordingSink final : public CheckedStlSink {
public:
    std::vector<unsigned char> bytes;
    size_t byte_limit = std::numeric_limits<size_t>::max();
    bool flush_ok = true, close_ok = true, reported_error = false;
    size_t closes = 0, flushes = 0;
    size_t write(const unsigned char *data, size_t size) override
    {
        const size_t count = std::min(size, byte_limit - std::min(byte_limit, bytes.size()));
        bytes.insert(bytes.end(), data, data + count);
        return count;
    }
    bool error() const override { return reported_error; }
    bool flush() override { ++flushes; return flush_ok; }
    bool close() override { ++closes; return close_ok; }
};
uint32_t read_u32(const unsigned char *data)
{
    uint32_t value = 0;
    for (unsigned i = 0; i < 4; ++i) value |= uint32_t(data[i]) << (8 * i);
    return value;
}
float read_float(const unsigned char *data)
{
    const uint32_t bits = read_u32(data);
    float result;
    std::memcpy(&result, &bits, sizeof(result));
    return result;
}
#ifndef _WIN32
namespace fs = boost::filesystem;
class Scratch {
public:
    Scratch() : root(fs::canonical(fs::temp_directory_path()))
    {
        for (int attempt = 0; attempt < 16; ++attempt) {
            const fs::path candidate = root / fs::unique_path("bambu-stl-export-test-%%%%-%%%%-%%%%");
            if (fs::create_directory(candidate)) { directory = fs::canonical(candidate); return; }
        }
        throw std::runtime_error("Cannot allocate isolated STL test directory.");
    }
    ~Scratch()
    {
        try {
            // Recursive cleanup is confined to this exclusively CREATED scratch,
            // never an env-supplied/user directory. Refuse a replaced symlink/path.
            if (directory.empty() || directory.parent_path() != root ||
                directory.filename().string().find("bambu-stl-export-test-") != 0 || fs::is_symlink(directory)) return;
            if (fs::canonical(directory) != directory) return;
            boost::system::error_code error;
            fs::remove_all(directory, error);
        } catch (...) {}
    }
    fs::path file(const char *name) const { return directory / name; }
    fs::path root, directory;
};
std::vector<unsigned char> own_bytes(const fs::path &path)
{
    std::ifstream input(path.string(), std::ios::binary);
    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}
size_t temporary_count(const fs::path &directory)
{
    size_t count = 0;
    for (fs::directory_iterator it(directory), end; it != end; ++it)
        if (it->path().filename().string().find(".bambu-stl-export-") == 0) ++count;
    return count;
}
#endif
} // namespace

TEST_CASE("Checked STL encoder writes actual little-endian finite oriented facets", "[strict_stl_writer]")
{
    const auto mesh = make_cube(10, 10, 10).its;
    RecordingSink sink;
    const auto result = write_binary_stl_checked(mesh, sink);
    REQUIRE(result.complete);
    REQUIRE(sink.bytes.size() == 84 + 50 * mesh.indices.size());
    CHECK(result.bytes_written == sink.bytes.size());
    CHECK(read_u32(sink.bytes.data() + 80) == mesh.indices.size());
    CHECK(sink.flushes == 1);
    CHECK(sink.closes == 1);
    for (size_t i = 0; i < mesh.indices.size(); ++i) {
        const unsigned char *facet = sink.bytes.data() + 84 + 50 * i;
        Vec3d normal;
        for (int axis = 0; axis < 3; ++axis) normal[axis] = read_float(facet + 4 * axis);
        REQUIRE(normal.allFinite());
        CHECK(normal.norm() == Approx(1.0).margin(1e-6));
        const auto &face = mesh.indices[i];
        const Vec3d cross = (mesh.vertices[face[1]].cast<double>() - mesh.vertices[face[0]].cast<double>()).cross(
                            mesh.vertices[face[2]].cast<double>() - mesh.vertices[face[0]].cast<double>());
        CHECK(normal.dot(cross) > 0.0);
        for (int vertex = 0; vertex < 3; ++vertex)
            for (int axis = 0; axis < 3; ++axis)
                CHECK(read_float(facet + 12 + 12 * vertex + 4 * axis) == mesh.vertices[face[vertex]][axis]);
        CHECK(facet[48] == 0);
        CHECK(facet[49] == 0);
    }
}

TEST_CASE("Checked STL encoder reports short-write flush error and close failure", "[strict_stl_writer]")
{
    const auto mesh = make_cube(10, 10, 10).its;
    for (size_t limit : {size_t(1), size_t(82), size_t(100)}) {
        INFO("Short write after " << limit << " bytes");
        RecordingSink sink; sink.byte_limit = limit;
        const auto result = write_binary_stl_checked(mesh, sink);
        CHECK_FALSE(result.complete);
        CHECK_FALSE(result.message.empty());
        CHECK(sink.closes == 1);
    }
    SECTION("flush failure") {
        RecordingSink sink; sink.flush_ok = false;
        CHECK_FALSE(write_binary_stl_checked(mesh, sink).complete);
        CHECK(sink.closes == 1);
    }
    SECTION("close failure") {
        RecordingSink sink; sink.close_ok = false;
        CHECK_FALSE(write_binary_stl_checked(mesh, sink).complete);
        CHECK(sink.closes == 1);
    }
    SECTION("ferror despite a full write") {
        RecordingSink sink; sink.reported_error = true;
        CHECK_FALSE(write_binary_stl_checked(mesh, sink).complete);
        CHECK(sink.closes == 1);
    }
    SECTION("invalid index checked before serialization") {
        auto invalid = mesh; invalid.indices.back()[0] = -1;
        RecordingSink sink;
        CHECK_FALSE(write_binary_stl_checked(invalid, sink).complete);
        CHECK(sink.bytes.empty());
        CHECK(sink.closes == 1);
    }
}

TEST_CASE("Atomic STL export rejects Cancel before writing", "[strict_stl_writer]")
{
    const auto mesh = make_cube(10, 10, 10).its;
    CHECK_FALSE(write_stl_export_atomic("", mesh, StlDestinationPolicy::NoReplace).committed);
    CHECK_FALSE(write_stl_export_atomic("<cancel>", mesh, StlDestinationPolicy::NoReplace).committed);
}

#ifndef _WIN32
TEST_CASE("Atomic STL writes only complete files inside isolated scratch", "[.strict_stl_export_io][strict_stl_writer_io]")
{
    Scratch scratch;
    const auto target = scratch.file("mesh.stl");
    const auto mesh = make_cube(10, 10, 10).its;
    const auto first = write_stl_export_atomic(target.string(), mesh, StlDestinationPolicy::NoReplace);
    REQUIRE(first.committed);
    CHECK(first.cleanup_warning.empty());
    const auto before = own_bytes(target);
    CHECK(before.size() == 84 + 50 * mesh.indices.size());
    CHECK(temporary_count(scratch.directory) == 0);

    const auto collision = write_stl_export_atomic(target.string(), mesh, StlDestinationPolicy::NoReplace);
    CHECK_FALSE(collision.committed);
    CHECK(own_bytes(target) == before);

    auto invalid = mesh; invalid.indices.pop_back();
    CHECK_FALSE(write_stl_export_atomic(target.string(), invalid, StlDestinationPolicy::ReplaceConfirmedRegularFile).committed);
    CHECK(own_bytes(target) == before);
    CHECK(temporary_count(scratch.directory) == 0);

    const auto replacement = make_cube(20, 20, 20).its;
    CHECK(write_stl_export_atomic(target.string(), replacement, StlDestinationPolicy::ReplaceConfirmedRegularFile).committed);
    CHECK(own_bytes(target) != before);
    CHECK(temporary_count(scratch.directory) == 0);

    const auto alias = scratch.file("alias.stl");
    fs::create_symlink(target, alias);
    const auto replaced_bytes = own_bytes(target);
    CHECK_FALSE(write_stl_export_atomic(alias.string(), mesh, StlDestinationPolicy::ReplaceConfirmedRegularFile).committed);
    CHECK(fs::is_symlink(alias));
    CHECK(own_bytes(target) == replaced_bytes);

    const auto missing = scratch.directory / "missing" / "mesh.stl";
    CHECK_FALSE(write_stl_export_atomic(missing.string(), mesh, StlDestinationPolicy::NoReplace).committed);
    CHECK_FALSE(fs::exists(missing.parent_path()));
    CHECK_FALSE(write_stl_export_atomic(scratch.directory.string(), mesh, StlDestinationPolicy::ReplaceConfirmedRegularFile).committed);
    CHECK(temporary_count(scratch.directory) == 0);
}
TEST_CASE("Atomic STL NoReplace obeys real filesystem case and Unicode equivalence", "[.strict_stl_export_io][strict_stl_writer_io]")
{
    Scratch scratch;
    const auto mesh = make_cube(10, 10, 10).its;
    const std::vector<std::pair<std::string, std::string>> names = {
        {"Case.stl", "case.stl"}, {"caf\xc3\xa9.stl", "cafe\xcc\x81.stl"}};
    for (const auto &pair : names) {
        const auto first_path = scratch.directory / pair.first;
        const auto second_path = scratch.directory / pair.second;
        REQUIRE(write_stl_export_atomic(first_path.string(), mesh, StlDestinationPolicy::NoReplace).committed);
        const auto first_bytes = own_bytes(first_path);
        // No naive case-sensitive reservation set: ask THIS filesystem.
        const bool already_exists = fs::exists(second_path);
        const auto second = write_stl_export_atomic(second_path.string(), mesh, StlDestinationPolicy::NoReplace);
        CHECK(second.committed == !already_exists);
        CHECK(own_bytes(first_path) == first_bytes);
    }
    CHECK(temporary_count(scratch.directory) == 0);
}
#else
TEST_CASE("Atomic STL Windows draft refuses unsupported writes", "[strict_stl_writer]")
{
    const auto result = write_stl_export_atomic("must-not-be-created.stl", make_cube(10, 10, 10).its,
                                               StlDestinationPolicy::NoReplace);
    CHECK_FALSE(result.committed);
    CHECK(result.message.find("Unsupported") != std::string::npos);
}
#endif
