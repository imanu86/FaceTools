// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#include "CheckedStlExport.hpp"
#include "libslic3r/MeshSolidValidation.hpp"

#include <boost/filesystem.hpp>
#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <system_error>
#include <utility>
#include <vector>
#ifndef _WIN32
#include <unistd.h>
#endif

namespace Slic3r {
namespace {
void put_u32(unsigned char *output, uint32_t value)
{
    for (unsigned i = 0; i < 4; ++i) output[i] = static_cast<unsigned char>((value >> (8 * i)) & 255u);
}
void put_float(unsigned char *output, float value)
{
    static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559,
                  "Binary STL requires IEEE-754 binary32.");
    uint32_t bits;
    std::memcpy(&bits, &value, sizeof(bits));
    put_u32(output, bits);
}
Vec3f stored_normal(const indexed_triangle_set &mesh, const Vec3i &face)
{
    const Vec3d a = mesh.vertices[size_t(face[0])].cast<double>();
    const Vec3d b = mesh.vertices[size_t(face[1])].cast<double>();
    const Vec3d c = mesh.vertices[size_t(face[2])].cast<double>();
    const Vec3d cross = (b - a).cross(c - a);
    const double length = cross.norm();
    if (!cross.allFinite() || !std::isfinite(length) || length <= 0.0)
        throw std::runtime_error("STL facet has no finite normal.");
    const Vec3f normal = (cross / length).cast<float>();
    if (!normal.allFinite() || !std::isfinite(normal.squaredNorm()) || normal.squaredNorm() <= 0.f ||
        normal.cast<double>().dot(cross) <= 0.0)
        throw std::runtime_error("STL normal cannot be represented without changing orientation.");
    return normal;
}
void validate_serializable(const indexed_triangle_set &mesh)
{
    validate_solid_mesh(mesh);
    if (mesh.indices.size() > size_t(std::numeric_limits<uint32_t>::max()))
        throw std::runtime_error("STL triangle count exceeds the binary format.");
    // Check every normal BEFORE opening/creating a destination temporary file.
    for (const Vec3i &face : mesh.indices) (void)stored_normal(mesh, face);
}
#ifndef _WIN32
std::string os_error(const char *operation, int value)
{
    return std::string(operation) + ": " + std::error_code(value, std::generic_category()).message();
}
class FileSink final : public CheckedStlSink {
public:
    explicit FileSink(FILE *file) : m_file(file) {}
    ~FileSink() override { if (m_file) std::fclose(m_file); }
    size_t write(const unsigned char *bytes, size_t size) override { return std::fwrite(bytes, 1, size, m_file); }
    bool error() const override { return m_file && std::ferror(m_file) != 0; }
    bool flush() override { return std::fflush(m_file) == 0; }
    bool close() override
    {
        FILE *file = m_file;
        m_file = nullptr;
        return file && std::fclose(file) == 0;
    }
private:
    FILE *m_file;
};
class OwnedTemporary {
public:
    OwnedTemporary(std::string path, std::string parent) : m_path(std::move(path)), m_parent(std::move(parent)) {}
    ~OwnedTemporary() { try { if (m_owned) (void)cleanup(); } catch (...) {} }
    void created_as(const char *path) noexcept
    {
        // mkstemp changes only the final six characters; storage is allocated
        // before creation, so establishing ownership cannot allocate or throw.
        std::copy_n(path, m_path.size(), m_path.begin());
        m_owned = true;
    }
    const std::string &path() const { return m_path; }
    void released_by_rename() noexcept { m_owned = false; }
    std::string cleanup()
    {
        if (!m_owned) return {};
        m_owned = false;
        // Delete ONLY the exact absolute filename created exclusively by mkstemp.
        // Never delete a destination or recursively remove its directory.
        const boost::filesystem::path path(m_path);
        if (!path.is_absolute() || path.parent_path().string() != m_parent ||
            path.filename().string().compare(0, 18, ".bambu-stl-export-") != 0)
            return "Refused cleanup of an unexpected temporary path: " + m_path;
        if (::unlink(m_path.c_str()) != 0 && errno != ENOENT)
            return os_error("Could not remove owned STL temporary", errno) + " (" + m_path + ")";
        return {};
    }
private:
    std::string m_path, m_parent;
    bool m_owned = false;
};
void check_destination(const boost::filesystem::path &destination, StlDestinationPolicy policy)
{
    namespace fs = boost::filesystem;
    boost::system::error_code error;
    const auto status = fs::symlink_status(destination, error);
    if (error && error != boost::system::errc::make_error_condition(boost::system::errc::no_such_file_or_directory))
        throw std::runtime_error("Cannot inspect STL destination: " + error.message());
    if (status.type() == fs::file_not_found) return;
    if (status.type() != fs::regular_file)
        throw std::runtime_error("STL destination is not a regular file (symlinks are unsupported).");
    if (policy == StlDestinationPolicy::NoReplace)
        throw std::runtime_error("STL destination already exists.");
}
#endif
} // namespace

StlStreamWriteResult write_binary_stl_checked(const indexed_triangle_set &mesh, CheckedStlSink &sink)
{
    StlStreamWriteResult result;
    bool closed = false;
    try {
        validate_serializable(mesh);
        const auto emit = [&](const unsigned char *bytes, size_t size) {
            const size_t written = sink.write(bytes, size);
            result.bytes_written += written;
            if (written != size || sink.error()) throw std::runtime_error("STL write failed or was short.");
        };
        std::array<unsigned char, 80> header{};
        constexpr char label[] = "BambuStudio checked solid STL";
        std::memcpy(header.data(), label, sizeof(label) - 1);
        emit(header.data(), header.size());
        std::array<unsigned char, 4> count{};
        put_u32(count.data(), uint32_t(mesh.indices.size()));
        emit(count.data(), count.size());
        for (const Vec3i &face : mesh.indices) {
            std::array<unsigned char, 50> facet{};
            const Vec3f normal = stored_normal(mesh, face);
            for (size_t axis = 0; axis < 3; ++axis) put_float(facet.data() + 4 * axis, normal[int(axis)]);
            for (size_t corner = 0; corner < 3; ++corner)
                for (size_t axis = 0; axis < 3; ++axis)
                    put_float(facet.data() + 12 + 12 * corner + 4 * axis,
                              mesh.vertices[size_t(face[int(corner)])][int(axis)]);
            emit(facet.data(), facet.size());
        }
        if (!sink.flush() || sink.error()) throw std::runtime_error("STL flush failed.");
        closed = true; // Even a failed/throwing close must not be attempted twice.
        if (!sink.close()) throw std::runtime_error("STL close failed.");
        result.complete = true;
    } catch (const std::exception &error) {
        result.message = error.what();
    } catch (...) {
        result.message = "STL serialization failed.";
    }
    if (!closed) {
        try {
            if (!sink.close()) result.message += " STL cleanup close also failed.";
        } catch (...) { result.message += " STL cleanup close threw."; }
    }
    return result;
}

StlFileWriteResult write_stl_export_atomic(const std::string &destination_name,
    const indexed_triangle_set &mesh, StlDestinationPolicy policy)
{
    StlFileWriteResult result;
    if (destination_name.empty() || destination_name == "<cancel>") {
        result.message = "No STL destination was selected.";
        return result;
    }
#ifdef _WIN32
    (void)mesh;
    (void)policy;
    result.message = "Unsupported: checked atomic STL export is not implemented for Windows in this draft.";
    return result; // BEFORE any write/temp creation.
#else
    try {
        if (policy != StlDestinationPolicy::NoReplace && policy != StlDestinationPolicy::ReplaceConfirmedRegularFile)
            throw std::runtime_error("Invalid STL destination policy.");
        validate_serializable(mesh); // Invalid input must not create even a temporary file.
        namespace fs = boost::filesystem;
        const fs::path supplied(destination_name);
        if (supplied.filename().empty() || supplied.filename() == "." || supplied.filename() == "..")
            throw std::runtime_error("Invalid STL destination filename.");
        const fs::path parent = fs::canonical(supplied.has_parent_path() ? supplied.parent_path() : fs::path("."));
        const fs::path destination = parent / supplied.filename();
        check_destination(destination, policy);
        std::string pattern = (parent / ".bambu-stl-export-XXXXXX").string();
        std::vector<char> buffer(pattern.begin(), pattern.end());
        buffer.push_back('\0');
        OwnedTemporary temporary(pattern, parent.string());
        const int descriptor = ::mkstemp(buffer.data());
        if (descriptor < 0) throw std::runtime_error(os_error("Cannot create STL temporary", errno));
        temporary.created_as(buffer.data());
        try {
            FILE *file = ::fdopen(descriptor, "wb");
            if (!file) {
                const int error = errno;
                const int close_error = ::close(descriptor) == 0 ? 0 : errno;
                result.message = os_error("Cannot open STL temporary stream", error);
                if (close_error) result.message += " " + os_error("Temporary descriptor close failed", close_error);
                result.cleanup_warning = temporary.cleanup();
                return result;
            }
            FileSink sink(file);
            const auto encoded = write_binary_stl_checked(mesh, sink);
            result.bytes_written = encoded.bytes_written;
            if (!encoded.complete) {
                result.message = encoded.message;
                result.cleanup_warning = temporary.cleanup();
                return result;
            }
            // Recheck known destination types immediately before commit. ReplaceConfirmed
            // is not a concurrency-security primitive; NoReplace also uses atomic link.
            check_destination(destination, policy);
            if (policy == StlDestinationPolicy::NoReplace) {
                if (::link(temporary.path().c_str(), destination.string().c_str()) != 0) {
                    result.message = os_error("Atomic STL publication failed (no overwrite)", errno);
                    result.cleanup_warning = temporary.cleanup();
                    return result;
                }
                result.committed = true;
                result.cleanup_warning = temporary.cleanup();
            } else {
                if (::rename(temporary.path().c_str(), destination.string().c_str()) != 0) {
                    result.message = os_error("Atomic STL replacement failed", errno);
                    result.cleanup_warning = temporary.cleanup();
                    return result;
                }
                temporary.released_by_rename();
                result.committed = true;
            }
        } catch (...) {
            result.cleanup_warning = temporary.cleanup();
            throw;
        }
    } catch (const std::exception &error) {
        result.message = error.what();
    } catch (...) {
        result.message = "Checked STL export failed.";
    }
    return result;
#endif
}
} // namespace Slic3r
