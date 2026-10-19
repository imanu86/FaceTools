# FaceTools

*Unofficial planar-face tools for Bambu Studio — AGPL-3.0 fork.*

> **Unofficial.** An independent modification of [Bambu Studio](https://github.com/bambulab/BambuStudio) (AGPL-3.0). Not affiliated with, endorsed or supported by Bambu Lab. "Bambu Studio" and "Bambu Lab" are trademarks of their owners.

FaceTools adds direct, face-based mesh editing tools to the Prepare view, aimed at quick fixes of already-modelled STL/3MF parts without going back to CAD.

## Features
- **Extrude / move planar faces** as a new part, or union/subtract with the source, with live preview and numeric or drag distance.
- **Point-defined profiles** for extrusion.
- **Position by point coincidence** — snap one part's point onto another's.
- **Cut and extend** — cut a part at a face-defined plane and move/stretch one side (guide faces select the connected region to move).
- **Native Cut from a face**, **Smart fill by angle**, and support for **faces of negative parts** (cutters).
- Undo/Redo integration and unit tests for the geometry core (`tests/libslic3r/test_face_extrusion*.cpp`).

## Status and honest limits
- **Prototype.** It is a mesh-based tool, **not a parametric CAD kernel**: no BRep, no curved-surface editing, no shortening, no general "Move Face".
- **GUI acceptance is only partial.** Automated tests pass, but the manual GUI scenarios (see [docs/facetools/TEST_CHECKLIST.md](docs/facetools/TEST_CHECKLIST.md)) are not all verified. **"Cut and extend" is known to need more work** on some shapes — expect bugs and always keep a backup of your project.
- **Printing to a cloud printer is not authorised for forks.** Bambu restricts "Print" for unsigned third-party builds; sending files may work, but direct Print can fail. No attempt is made to bypass this. The proprietary Bambu network plugin is **not** included.
- macOS arm64 only so far; binaries (if provided) are ad-hoc signed, not notarized.

## Build / install
Follow the upstream build guide (https://github.com/bambulab/BambuStudio/wiki) and `BuildMac.sh`. Needs Homebrew `zstd` and `imath` among others. Releases, when available, list the exact source commit.

## License
AGPL-3.0 — see `LICENSE`. This is a derivative work of Bambu Studio, PrusaSlicer and Slic3r; third-party licences stay with their code. Changes are listed in `FACETOOLS_CHANGES.md`. Source for any binary: the tagged commit of this repository.


---

# Upstream README (Bambu Studio)

![image](https://user-images.githubusercontent.com/106916061/179006347-497d24c0-9bd6-45b7-8c49-d5cc8ecfe5d7.png)
# BambuStudio
Bambu Studio is a cutting-edge, feature-rich slicing software.  
It contains project-based workflows, systematically optimized slicing algorithms, and an easy-to-use graphic interface, bringing users an incredibly smooth printing experience.

Prebuilt Windows, macOS 64-bit and Linux releases are available through the [github releases page](https://github.com/bambulab/BambuStudio/releases/).

Bambu Studio is based on [PrusaSlicer](https://github.com/prusa3d/PrusaSlicer) by Prusa Research, which is from [Slic3r](https://github.com/Slic3r/Slic3r) by Alessandro Ranellucci and the RepRap community.

See the [wiki](https://github.com/bambulab/BambuStudio/wiki) and the [documentation directory](https://github.com/bambulab/BambuStudio/tree/master/doc) for more information.

# What are Bambu Studio's main features?
Key features are:
- Basic slicing features & GCode viewer
- Multiple plates management
- Remote control & monitoring
- Auto-arrange objects
- Auto-orient objects
- Hybrid/Tree/Normal support types, Customized support
- multi-material printing and rich painting tools
- multi-platform (Win/Mac/Linux) support
- Global/Object/Part level slicing parameters

Other major features are:
- Advanced cooling logic controlling fan speed and dynamic print speed
- Auto brim according to mechanical analysis
- Support arc path(G2/G3)
- Support STEP format
- Assembly & explosion view
- Flushing transition-filament into infill/object during filament change

# How to compile
Following platforms are currently supported to compile:
- Windows 64-bit, [Compile Guide](https://github.com/bambulab/BambuStudio/wiki/Windows-Compile-Guide)
- Mac 64-bit, [Compile Guide](https://github.com/bambulab/BambuStudio/wiki/Mac-Compile-Guide)
- Linux, [Compile Guide](https://github.com/bambulab/BambuStudio/wiki/Linux-Compile-Guide)
  - currently we only provide linux appimages on [github releases](https://github.com/bambulab/BambuStudio/releases) for Ubuntu/Fedora, and a [flathub version](https://flathub.org/apps/com.bambulab.BambuStudio) can be used for all the linux platforms

# Report issue
You can add an issue to the [github tracker](https://github.com/bambulab/BambuStudio/issues) if **it isn't already present.**

# License
Bambu Studio is licensed under the GNU Affero General Public License, version 3. Bambu Studio is based on PrusaSlicer by PrusaResearch.

PrusaSlicer is licensed under the GNU Affero General Public License, version 3. PrusaSlicer is owned by Prusa Research. PrusaSlicer is originally based on Slic3r by Alessandro Ranellucci.

Slic3r is licensed under the GNU Affero General Public License, version 3. Slic3r was created by Alessandro Ranellucci with the help of many other contributors.

The GNU Affero General Public License, version 3 ensures that if you use any part of this software in any way (even behind a web server), your software must be released under the same license.

The bambu networking plugin is based on non-free libraries. It is optional to the Bambu Studio and provides extended networking functionalities for users.
By default, after installing Bambu Studio without the networking plugin, you can initiate printing through the SD card after slicing is completed.

