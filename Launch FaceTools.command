#!/bin/zsh
# Local development launcher. Never reuse the installed Bambu Studio profile.
set -eu

ROOT="${0:A:h}"
BINARY="${ROOT}/build/facetools-arm64/src/BambuStudio.app/Contents/MacOS/BambuStudio"
PROFILE="${ROOT}/build/facetools-profile"
OUTPUT="${ROOT}/build/facetools-launch"

if [[ ! -x "${BINARY}" ]]; then
    print -u2 -- "The local FaceTools build is unavailable. Installed Bambu Studio will not be launched."
    exit 1
fi

mkdir -p "${PROFILE}" "${OUTPUT}"
LOG="${OUTPUT}/launch-$(date +%Y%m%d-%H%M%S)-$$.log"
print -- "Starting local FaceTools with a separate profile."
print -- "Profile: ${PROFILE}"
print -- "Log: ${LOG}"
print -- "Tool: Extrude face as new part (Prepare tab)."

# A normal double-click opens our synthetic frame, not a user's existing project.
DEMO="${ROOT}/build/face-extrusion-demo/planar-frame-demo.stl"
if (( $# == 0 )) && [[ -f "${DEMO}" ]]; then
    set -- "${DEMO}"
    print -- "Opening the sample frame for testing."
fi

exec "${BINARY}" --datadir "${PROFILE}" --outputdir "${OUTPUT}" "$@" >"${LOG}" 2>&1
