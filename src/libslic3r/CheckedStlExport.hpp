// Copyright (C) 2026 Emmanuele Semprini
// SPDX-License-Identifier: AGPL-3.0-only
// Part of FaceTools, an unofficial modification of Bambu Studio (AGPL-3.0).

#ifndef slic3r_CheckedStlExport_hpp_
#define slic3r_CheckedStlExport_hpp_

#include "libslic3r/TriangleMesh.hpp"

#include <cstddef>
#include <string>

namespace Slic3r {

// Narrow production seam: both the real FILE sink and deterministic short-write /
// flush / close tests use the SAME encoder. No alternate test serialization.
class CheckedStlSink {
public:
    virtual ~CheckedStlSink() = default;
    virtual size_t write(const unsigned char *bytes, size_t size) = 0;
    virtual bool error() const = 0;
    virtual bool flush() = 0;
    virtual bool close() = 0;
};
struct StlStreamWriteResult {
    bool complete = false;
    size_t bytes_written = 0;
    std::string message;
};
// Validates before emitting bytes; closes sink exactly once, including on error.
// No cancellation is offered while writing. The GUI must resolve Cancel beforehand.
StlStreamWriteResult write_binary_stl_checked(const indexed_triangle_set &mesh, CheckedStlSink &sink);

enum class StlDestinationPolicy { NoReplace, ReplaceConfirmedRegularFile };
struct StlFileWriteResult {
    bool committed = false;
    size_t bytes_written = 0;
    std::string message;
    std::string cleanup_warning; // A committed file remains a SUCCESS if temp cleanup fails.
};

// Export-only POSIX implementation (Mac/Linux). Windows returns Unsupported BEFORE
// any file write. Temporary file is exclusive and in DESTDIR. Complete encoding,
// ferror, fflush and fclose must succeed before atomic destination publication.
// NoReplace: link(temp,destination), failing closed on unsupported filesystems.
// ReplaceConfirmedRegularFile: rename(temp,destination), WITHOUT pre-unlink.
// Refuses existing symlinks/non-regular destinations. Caller chooses each name
// just before this call. No global-batch atomicity or power-loss durability claim.
// ReplaceConfirmed has no protection against hostile concurrent path replacement.
StlFileWriteResult write_stl_export_atomic(const std::string &destination,
    const indexed_triangle_set &mesh, StlDestinationPolicy policy);

} // namespace Slic3r
#endif
