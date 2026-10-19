# Experimental planar-face extrusion

FaceTools adds local mesh-editing tools to Bambu Studio, not a parametric CAD
kernel. Face and point-profile extrusions offer **New part**, **Union** and
**Subtract**. New part preserves the source; Boolean results replace only the
selected source part after a successful calculation. V4 also adds point alignment
in Prepare and rigid **Cut and stretch** with explicit plane and moving side.
Original licenses/notices remain unchanged. Development is versioned in a private
repository; no public FaceTools release or notarized distribution is provided.

## Delivery and verification status

**V11 is installed and verified after NEW explicit saved/closed consent, consumed
once.** The normal app was compiled from
`c2cfe1b564042474ee47062cf882320f1f8dbcdb` at `20261016-222520`; later docs-only
commits/tags are not a recompilation. Actual NORMAL geometry **158 / 115,438**, wx
**7 / 112** and separate hidden IO **2 / 26** pass. Do not add repeated CTest, XML
or private runs. The delivery ledger
keeps CTest907 native PASS distinct from its immutable wrapper failure, and records
normal build, package, five isolated CLI checks, signing and installed verification.
V10 app/guide backups were verified; user-data preservation is limited to the
observed installation checkpoint, not a byte proof of opaque profiles/projects.

See the [installed v11 contracts](<#v11-contracts--installed-build>),
historical private RED/GREEN
and source/tag publication procedure.
Live mouse/GL, Undo/Redo, save/reopen and slicing acceptance remain **NOT RUN**;
wx console tests and UNREPAIRED CLI statistics do not establish those outcomes.

### Installed v10 baseline — historical verification

**Current delivery: v10 installed after NEW explicit save/closure consent.** Faces
supports one explicitly selected negative part: NewParts stay negative, and
Union/Subtract modify only the cutter, without changing positive-body placement.
The tool isolates the selected cutter for picking; see the
[negative-source interaction and limits](<#negative-source-parts--v10>).

Normal arm64 app and integrated tests pass **125 cases / 93,593 assertions**;
strict/deep local signing and five isolated CLI checks pass. Compiled source is
`9d4f4d249d5d4cc9edfb925dedd8f2b3b27a2cc1`. Installed hash/signature are verified;
v9 app/guide backups and all user profiles, plugins, projects and fixture paths are
preserved. Native/static evidence is not live mouse/GL, Undo or slicing acceptance.
See current evidence and artifact hashes.

### Preserved v9 standard-Cut delivery

**V9 behavior, retained in v10:** Standard Cut starts with **Keep both
pieces in place** enabled, retaining both results as child parts of one object
without automatic separation, reorientation or loader bed placement. Reset restores
that default; temporary context changes do not erase the chosen output, and all
seven output choices are archived for Undo. Unsupported combinations are explained
and blocked, not silently converted into a moving/separate-object operation.
Select an individual child part in the Objects list before moving only that piece.
See the v9 standard-Cut workflow.

The user's failed v8 workflow was reproduced with actual production decisions and
the real native cutter: default output moved a test half 22 mm before loading.
The same desired regressions pass after the fix; the final focused suite has
**10 cases / 3,885 assertions**. Earlier forced-Parts tests were not acceptance of
that default GUI path. V7's face-derived plane/preview fix and V8's optional
**Smart fill by angle** remain intact; see the
smart-guide scope and native-angle limits.

The normal arm64 app build and integrated CTest passed **117 cases / 69,171
assertions**; strict/deep local signing and five isolated native CLI checks passed.
Compiled source: `43f0673e28fd4b81de94c4c07fd2277f43cf22e2`. After NEW explicit
save/closure confirmation, v9 was installed and its signature/hash verified again.
V8 and its guide are backed up alongside all earlier backups. Profiles, plugins,
projects and historical fixture paths remain preserved. See
current status, artifact hashes and evidence.

**P2S cloud Print authorization is separate:** the user confirms successful Send
and screenshots show the unsigned-software warning for Print. Local ad-hoc Mac
signing does not confer Bambu caller authorization. No public self-service procedure
was found to authorize this fork's cloud Print; vendor-approved integration is not
ruled out. No alternative mode/handoff was selected and direct Print remains
unresolved; see verified facts and official options.
Real mouse/GL, Undo/Redo, save/reopen and slicing acceptance remain separate.
The older delivery entries below are historical evidence, not proof that the new
GUI workflow has been accepted.

- The original application and the first custom application built successfully.
- The user opened the separate Mac app and confirmed extrusion and slicing.
- The delivered v2 refinement adds hover, multiselection, a distance arrow,
  camera gestures and explicit Boolean modes. The user confirmed hover/highlights
  and rotation, but reported release-only preview and an annotation refusal.
- Corrective v3 adds live preview scheduling, a non-interactive distance badge,
  relative 1 mm drag increments, and explicit annotation-loss consent that keeps
  the face selection. Numeric entry remains unrestricted.
- The full v3 `BambuStudio` and `face_extrusion_tests` targets rebuilt successfully.
  Normal CTest passed **27 cases / 30,459 assertions**. The staged app passed
  strict signature verification and native CLI startup/import of both fixtures.
- After the user confirmed saving/quitting, the actual v3 app was installed.
  The sibling profile/projects were preserved; the previous bundle was retained
  as a backup. Strict/deep signature verification passed on the installed bundle.
  Interactive v3 acceptance is still required; static/helper tests are not GUI tests.
- V4 implements ordered point profiles with holes, point coincidence in Prepare,
  and **Cut and extend**. The normal full build and dedicated
  integrated CTest passed **57 cases / 41,545 assertions**. Staging passed strict/deep
  signing, CLI startup and native import of all four fixtures.
- After fresh user save/quit consent and an owned-process check, v4 was installed
  and its signature/hash verified again. V3 and its guide are retained as backups;
  the existing profile/projects were not read or copied and sibling metadata was
  verified unchanged. The two new examples were added without overwriting user files.
  See status and v4 scope/limits.
  **No interactive acceptance is implied.**

- V5 was installed after fresh save/quit consent. Cut/extend follows **guide faces
  → independent yellow plane → normal translation**; the first guide no longer
  dictates the plane and there is no manual Along/Opposite choice. The print-send
  error UI is dispatched through the owning Job's native wx main-thread event queue.
- Normal app and both regression targets built successfully: **69 geometry cases /
  44,139 assertions + 7 wx cases / 112 assertions**, total **76 / 44,251**. The wx
  regression also produced 44 expected failures against the original v4 implementation.
  Native CLI help and all four fixture imports passed with a fresh isolated data dir.
- V5 strict/deep ad-hoc signing passed before and after installation. V4 and its guide
  are backed up; profile, plugin binary metadata and user sibling metadata are unchanged.
  A first attempt safely rolled back because Finder updated its cache; the successful
  retry excludes only that cache from sibling metadata equality. No project/profile
  contents were read/copied, no plugin reinstalled, and no print was sent.
- V5 evidence: compiled commit
  `be08d6bcd3575fd7bc5ad342ceb0ba114d31d721`, installed binary SHA-256
  `0cc4d0d5f8bb724d9b43e5b4a10999682e3434d7930aebb890f76e27107cf49e`.
  Real GUI acceptance and the original send-failure cause remain unverified.

- **V6 is installed** after new closure-for-update confirmation and process checks.
  Later-plate cuts use a nearby physical-mm computation frame without repositioning
  the model; any selected G1/G2/etc. guide can explicitly define the cut plane.
  Exact source-facet provenance and double-derived seed normals avoid approximate
  world-plane reconstruction; nonzero nearby offsets are never silently snapped.
- Normal arm64 app and integrated tests passed **85 cases / 48,435 assertions**;
  strict/deep signing and five isolated native CLI checks also passed. V5 and its
  guide were preserved, with no profile, plugin, project or example replacement.
  Compiled source: `6ed1ab01ff2c34c46f8122c2cc2b5111c691ab17`.
  See v6 artifact hashes and evidence.
  Interactive acceptance on the user's later-plate model remains separate.

## Interaction

1. In Prepare, select an FFF object with a single instance and open **Extrude faces**.
2. Hover a planar surface to highlight it. Click to select it.
3. **Shift+click** toggles additional parallel, concordantly oriented faces of the
   **same source part**. Different plane heights are allowed. Opposite normals
   and another source part are rejected without discarding the valid selection.
4. Set a shared, signed distance in physical millimeters, or drag the arrow along
   the world normal. The arrow is anchored to a real picked point on a face, not
   a bounding-box center that might lie inside a hole.
5. Choose **New part**, **Union** or **Subtract**, then explicitly apply.

**New part** is the default and creates one separate, closed volume per selected
patch, in one model operation: MODEL_PART for a positive source; the v10 negative
source workflow below creates NEGATIVE_VOLUME instead. There is no attempt to concatenate
intersecting shells and call them a valid single solid. Each patch retains its
cap triangulation and has boundary walls, including walls around holes.

Union/Subtract operate on the source plus all selected extrusion operands.
A red subtraction preview shows the **cutter**, not a precomputed Boolean result.
Negative distance does not implicitly mean subtraction: the operation is an
explicit choice. An empty final Boolean result is not installed as an empty part.

The arrow changes only the preview, never the model or undo history. Drag updates
are rendered live, with a non-interactive distance badge. A drag uses 1 mm steps
relative to its starting value; typed fractional values are not quantized.
Escape during an arrow drag restores its starting distance; otherwise Escape first
clears the face selection. Near-parallel viewing rays do not produce unstable
distance jumps: rotate the view slightly or use the numeric value.

Click and drag are different gestures. Dragging away from the arrow orbits the
camera without moving the object or changing the pinned faces. Cmd navigation,
right/middle pan and wheel navigation remain available. Button-up handling must
not fall through to ordinary object/plate selection. A shared Canvas cleanup
resets the orbit pivot even when the gizmo consumes the release.

## V4 interaction and geometry

### Point profiles

Use the separate four-point toolbar icon in Prepare. Pick at least three coplanar
**mesh vertices in perimeter order**, all on one positive source part. The first
pick pins the Boolean target. Vertex picking uses the hit triangle, not a
planar-feature or normal agreement requirement; the normal comes from the validated
polygon. Contours close implicitly; do not repeat the first vertex. Inner contours
explicitly define holes, which are not automatically inherited from the source.
Markers, perimeter, undo-last and clear controls expose the selection.

The initial extrusion is 10 mm. Signed world-mm distance, relative 1 mm drag steps
and New part / Union / Subtract reuse the protected face-tool application pipeline.
The [profile kernel](<src/libslic3r/PlanarProfile.cpp>) validates original and stored
float contours, planarity, boundary direction, area, manifoldness and connectivity.
It retains all input vertices; repairing an internal triangulation diagonal does
not move, merge or flatten them. Extrusion delegates to the existing kernel.
The 4096-point limit bounds validation work.

### Point coincidence in Prepare

Select parts of the same object, use **Assemble → Point and point assembly** in
Prepare, and pick **Point 1 (Fixed)** followed by **Point 2 (Moving)**. Do not switch
to the separate Assemble tab. Coincide points applies the complete translation in
one action; world XYZ target offsets mean Point 2 minus Point 1 and remain editable
at tiny/zero separation. Enter or Apply commits all components together.

The [translation plan](<src/libslic3r/PointAlignment.hpp>) changes only the moving
volume's offset and leaves its linear transform bit-exact. The object must have one
instance: volume poses are shared by its copies. Live mesh/ID/pose, whole-scene GL
state and native bed-placement checks precede the single snapshot. No-op creates
no history entry. Shift/Plane casts against the real mesh, not the raised highlight;
explicit circle-center Points remain valid even inside holes. Native Face/Face,
Circle/Circle distance controls and Assemble-tab poses are not reimplemented.

### Cut and stretch

The following v4–v6 account describes the historical whole-half operation. For
the installed v11 guided default, full-region confirmation and narrower scope
rules, use the [v11 contract](<#v11-contracts--installed-build>).
Standard native Cut and its v9 preserving default are separate and unchanged.

Use the separate split-L toolbar icon. The v4 interaction was rejected by the user:
using the first picked face as the plane and manually choosing Along/Opposite did
not match the requested face-first operation. **Installed v5 correction** (see
current delivery status) separates three phases:

1. **Guides:** click to add/remove planar faces of the same source. Different
   outward normals are allowed here; the first face does not orient the cut plane.
2. **Plane:** choose an independent yellow XY/XZ/YZ plane, or a plane based on a
   face picked in this second phase. Click, enter its absolute world coordinate or
   drag the plane handle. This handle moves the plane, not the extrusion distance.
3. **Distance:** confirm the plane, prepare the split and move the orange portion
   normally to the plane with the arrow/numeric distance. Back preserves guides;
   only Apply changes the model. The gray portion remains fixed.

The clicked guide positions determine a consistent moving side, not a majority
vote or averaged normal. A truly coplanar parallel guide identifies its material
side (opposite its outward normal); near-plane uncertainty or conflicting guides
is rejected. This identifies the entire connected half, **not arbitrary individual
CAD face offsets**. Both cut halves must still satisfy the unchanged kernel checks.
Plane and distance drags capture their parameter separately; Esc restores the
correct one. Changing steps/plane/guides invalidates stale jobs and restores source
visibility. Heavy preparation starts only after plane confirmation.

The [kernel](<src/libslic3r/CutAndStretch.cpp>) keeps the stationary half fixed,
translates the moving half rigidly and creates a bridge only over actual contact.
It never unions the unchanged original L back into the result. Native cutting is
uncapped; constrained triangulation caps only the generated cut boundary, preserving
its 3D vertices and collinear seam points. Exact projected predicates prove that
the stationary cap is included in the moving cap; unsupported partial contacts are
rejected rather than filled. Through-holes are retained. Apply re-extrudes the
original bridge base at the actual distance; it does not scale a unit prism.

Preparation is asynchronous and cached; distance dragging does not run CSG. The
preview distinguishes stationary gray, moving orange and bridge blue. The original
source is hidden only while a valid replacement preview is available. Applying
replaces only that source, using the shared stale/cancel/annotation/snapshot guards.

**Positive extension only**; zero is a non-applicable no-op and shortening is not
implemented. One connected surface component is required: through-holes are allowed,
but separate cavity shells/disconnected surfaces and unsupported shelf contacts are
rejected. Source/cap limits are 2,000,000/20,000 triangles. Native cutting's approximately
1e-3 mm band is not permission to snap geometry; ambiguous near-band cases and float
precision beyond the 1e-4 mm budget are refused.

V5 example: select the L arm's top (optionally also its lateral face), then choose
an independent **XY plane at absolute Z=20**, confirm and move +10 mm. The chosen
faces identify the +Z portion without an Along/Opposite control. Total height
becomes 34 mm while the arm stays 4 mm thick. For the tube, select its upper portion,
then set XY at **Z=15** and move +10: height becomes 40 mm and the 10×10 mm hole stays
open. The workflow regression deliberately selects a lateral guide FIRST, so its
normal cannot accidentally orient the plane or the motion.

#### V6 installed correction: later plates and guide-face planes

The v5 world-float budget incorrectly rejected small parts placed far from the
world origin, including ordinary later print plates. V6 prepares directly from
the original local mesh and its double placement into a nearby physical-mm cut
frame. Preview still uses the actual world placement. The result returns directly
to the original source-local frame, with physical-mm float-storage and round-trip
checks; neither the instance nor the source part is repositioned.

In the Plane step, choose **G1, G2, etc.** and press **Use guide face as plane**.
Any selected guide can define the plane, not just the first. The guide remains
selected for movement. Clicking that same guide in **From a face** also reuses it.
For example, select the L arm's **underside**, use it as the plane, confirm, and
move +10 mm: the arm remains 4 mm thick and the total height becomes 34 mm, on the
first or a later plate.

An unoffset face-defined plane carries its original source triangle into the
kernel. Exact source-plane predicates, rather than inference from nearby world
coordinates, preserve genuine coplanarity through rotated/scaled placements.
A numeric or dragged offset uses the separate general-plane path; nonzero nearby
offsets are not silently snapped back to a face. The existing genuine-partition,
contact, hole, topology, cancellation, annotation and bed-placement checks remain.
The Cut picker derives the original seed triangle normal in double precision;
it does not reuse the measurement cache's float-derived normal. See
installed artifact and validation evidence and
scope/limits. Live GUI acceptance remains separate.

## Negative source parts — v10

**Delivered in installed v10**, with native RED→GREEN, normal build, local signing,
isolated CLI and fresh-consent installation verified in STATUS.
This does not claim live mouse/GL, Undo/Redo or user-model slicing acceptance.

In **Objects**, explicitly select **one negative part**, then open **Extrude faces**
in FFF Prepare. The containing object must have a positive body and one instance.
Parent/mixed selections, multiple negative parts, parameter modifiers and support
volumes are not eligible negative sources. Point profiles and Cut and extend
remain positive-part-only.

The selected cutter is temporarily isolated so that a positive body does not hide
its faces. Picking targets that negative mesh, not the final sliced cavity surface.
The lower bed clipping cutoff is temporarily removed for that isolated cutter;
other camera/mesh validity restrictions remain. Leaving the tool restores the
scene visibility and clipping state. Clear selection/cancel calculation stay in
the tool, so isolation continues while editing that cutter. The physical plate
may still occlude below-bed geometry depending on the view. Actual interaction
and restoration still need GUI acceptance.

- **New part:** one new **negative volume** per selected patch; source and existing
  volumes are unchanged. New negative parts inherit no extruder/configuration,
  material, painting or source GUID.
- **Union:** enlarge the negative cutter and therefore the region it removes.
- **Subtract:** remove geometry from the cutter and therefore reduce that region.

The Boolean operands are the negative source and extrusion geometry, **not the
positive host**. Neither the geometric normal, distance sign nor CSG operation is
automatically inverted. Preview geometry is the extrusion operand, not a promised
preview of the final sliced cavity. Existing annotation consent and procedural,
Cut and exterior-data safeguards still apply.

Source roles are pinned at selection and rechecked before preview/Apply, because
changing a volume's type need not change its mesh or IDs. Negative geometry does
not contribute to printable bounds or the positive-part count. It may extend
below the bed, but an operation is still rejected if the normal model refresh
would translate the positive body, even by a tiny nonzero amount.

Persistence regressions pass through actual **native Bambu project 3MF** store/load
for negative NewParts and Boolean results, not the generic legacy 3MF path. A separate legacy negative-role token mismatch was found statically;
it is not reproduced/fixed by this feature, and legacy negative roundtrip is not
claimed. GUI Undo/Redo, native project save/reopen, actual cavity slicing and
visibility/clipping restoration remain explicit acceptance checks.

<a id="v11-candidate-contracts--not-installed"></a>
## V11 contracts — installed build

These contracts describe the installed v11 build from compiled c2.
STATUS separates actual normal
build/runtime, package and installation evidence from historical private runs and
still-unexecuted live GUI acceptance.

### Guide-associated Cut and stretch

- The GUI defaults to **Guide-associated piece**, locked to **one ModelVolume**.
  ALL original facet IDs in all guide masks must associate with ONE connected
  upper cut component, not just the clicked seeds. Wholly fixed, unassociated,
  ambiguous, duplicate or out-of-range IDs refuse; there is no nearest, majority
  or automatic whole-half fallback. Coplanar association uses exact source-plane,
  adjacency and material-side checks, not a near-plane guess.
- Contact uses whole lower-cap triangles contained in the chosen footprint;
  partial overlap refuses instead of clipping/filling it. The stationary remainder
  is a **real Boolean union of lower material and excluded upper components**, not
  concatenated closed shells. Moving/bridge contact, overlap or containment with
  an excluded upper region refuses. Apply replaces only the locked source.
- With a suitable separating plane, a U can move only its selected left OR right
  leg. The native synthetic U case has source/stationary/moving volumes
  **7,000 / 5,500 / 1,500 mm³**, output **8,000 mm³**, selected top **40 mm** and
  other top **30 mm**. A T's still-connected bar is ONE region: both arms move;
  guides do not invent a CAD branch boundary. **Whole moving half-space** is an
  explicit alternative, never a fallback. Choose another plane if orange is too broad.
- Positive extension only; zero is a no-op and shortening remains unsupported.
  The source still requires one connected surface boundary: separate closed cavity
  shells/multiple source boundaries are unsupported. Only the upper cut collection
  may split into separately validated components (at most 64). Through-holes remain
  supported where the contact checks allow them. General solid validation/export
  support for cavities or disjoint solids is NOT a global Cut-and-stretch promise.

Inspect and explicitly confirm the **ENTIRE orange moving region**, not only the
guide highlights. Confirmation binds immutable plan identity AND scope. Changing
source, pose, plane, masks, scope, plan or clearing invalidates it; a distance-only
change on the same region may retain it. Current GUI completion guards also reject
stale/canceled/reordered source state before mutation. Weak-identity/ABA tests pass,
but a matching GL cache identity is not proof of successful draw/presentation.
**Live preview, GUI Undo/Redo and slicing have not been run for v11.** See the
[kernel API](<src/libslic3r/CutAndStretch.hpp>) and
[confirmation policy](<src/slic3r/GUI/Gizmos/CutAndStretchScopeConfirmation.hpp>).
The working **Smart fill by angle** filter and standard native Cut are unchanged.

### Solid storage and physical precision

The committed [solid validator](<src/libslic3r/MeshSolidValidation.cpp>) checks finite
coordinates, indices, degenerate/duplicate faces, unused vertices, topology including
vertex fans and edge winding, positive oriented volume/bounds, CGAL self-intersections
and nested cavity orientation. The public CGAL volume guard rejects null-halfedge
orphan vertices. [New-part storage](<src/slic3r/GUI/Gizmos/FaceExtrusionModelPolicy.hpp>)
validates the ACTUAL inverse-stored F32 mesh and forward-materialized F32 world mesh
before mutation, with no recentering on insertion. Positional error alone is not
proof of valid topology. Valid through-holes, cavities and disjoint solids are kept
where the particular API supports them; unsupported geometry refuses without
silent welding, snapping, filling, remeshing, flipping or repair.

The **1e-4 mm** float-storage budget is per API/transform sequence, with a double
reference from original vertices; raw FFF placement uses `(instance * volume) *
original`. Boolean operand-error amplification plus final-rounding bounds are NOT
a global CSG-fidelity proof. Separate API calls reset the reference, and the GUI's
final translation is another stage: no end-to-end workflow bound is claimed.

### Strict FFF STL export is a separate path

The [strict export API](<src/libslic3r/FFFStlExportSafety.hpp>) and its
[Plater consumer](<src/slic3r/GUI/Plater.cpp>) do not change legacy combine,
ObjectList, SLA or global Boolean adapters. This checked path refuses SLA.
Results are typed **Success / Empty / Failure**. A legitimate Empty never falls
back to raw positives; it contributes nothing, and an all-empty batch is not
writable. Geometry failure publishes no partial output vector.

Explicit Boolean export unions positive parts FIRST, then subtracts negative
cutters **within each object**, independently of volume ordering. Different objects
and instances remain UNFUSED: concatenation is followed by solid/contact/overlap
validation, not an implicit union. The actual enabled OpenMeshCraft path through
`mcut::make_boolean` checks failure, completion and cancellation. Tests cover initial
and checkpoint-3 cancellation, NOT terminal OMC cancellation. Export's modal stale
check covers IDs/indices/selection mode, not a deep mesh/instance/pose/technology
snapshot; it is not the stronger Cut/face-tool source guard.

The final geometry batch gate runs after the LAST transform/merge and before
generated output-name discovery or writes (the initial destination dialog is
separate). Empty outputs preserve request/name mapping. No silent repair is added.

### Checked file publication and remaining delivery gates

The [POSIX writer](<src/libslic3r/CheckedStlExport.hpp>) creates an exclusive owned
temporary in the destination directory; checked writes, stream error, flush and
close must all succeed. Atomic rename is only for an explicitly confirmed regular
replacement; **NoReplace uses a hard link**, so new destinations require filesystem
hard-link support. Symlinks/non-regular destinations refuse; Windows refuses before
writing. This is **per-file atomicity**, not batch atomicity, power-loss durability
or hostile-concurrency protection. Earlier completed files can remain if a later
file fails. Two isolated IO cases now pass **26 assertions**, including replacement,
NoReplace, invalid-input preservation and case/Unicode aliases on this filesystem;
this is not universal filesystem compatibility.

Normal build/runtime, source-to-bundle no-follow preflight, package/signature,
five isolated CLI checks and installation verification have passed their separate
recorded gates; the fresh one-use saved/closed consent is consumed. These results
do not establish GUI acceptance or remote publication. See the
delivery ledger and source/tag procedure.
No printer/cloud settings change: P2S Print trust remains unresolved; historical
Send success is separate. No bypass, LAN/developer-mode, plugin/account/host edits.
The local arm64/internal1/public0/deployment11 build uses cached dependencies;
Homebrew zstd/Imath target newer macOS. No older-macOS runtime compatibility,
portability or notarization is promised.

## Scope and safety

- One object instance in FFF Prepare. Faces supports positive parts and the single
  explicitly selected negative source described above; the other extrusion modes
  stay positive-only. No assembly/exploded transforms, SLA, arbitrary curved-surface
  editing or automatic mesh repair.
- Each selected patch must separately be edge-connected, planar, consistently
  oriented and valid. Multiselection does not relax the kernel's topology checks.
- Instance/volume transforms are baked with handedness correction. World normals
  use the inverse transpose, so scaling/reflection do not alter the requested mm.
- All New part outputs are prepared before the one snapshot/model mutation.
  Original mesh/configuration/transforms/annotations remain unchanged. New parts
  inherit the effective extruder only for positive outputs, not unrelated settings
  or face painting. New negative outputs inherit no configuration.
- Boolean calculations run on copies in the existing UI worker. They never
  manipulate live GUI/Model objects on the worker thread or cancel unrelated jobs.
  Finalization checks cancellation, lifetime/epoch, IDs, source/all-volume roles,
  mesh identities, transformations and positive-part count before taking a
  snapshot or applying a result.
- Boolean mode lists detected support/seam/MMU/fuzzy painting and non-default
  legacy mesh tags. Explicit consent authorizes source-only removal at a successful
  commit; clicking the checkbox never clears annotations or the face selection.
  Complete content stamps are checked at consent, enqueue and finalization.
- Procedural text/emboss, connectors and structural Cut relationships remain
  non-bypassable for Boolean mode. The dedicated `exterior_facets` field lacks
  native Undo serialization: its presence anywhere in the Model blocks all three
  tool applications, including New part. Legacy `eExteriorAppearance` mesh tags
  are different and remain consentible. This does not fix other native tools' Undo.
- Boolean replacement retains source name/configuration values/material/transform/
  part GUID; cache IDs are renewed as in Simplify and stale reload provenance is
  discarded for that modified part only. Unrelated parts are not cleared.
- Post-operation bed bounds are checked using the actual model/hull arithmetic.
  The single-model-part and multiple-part `ensure_on_bed` branches differ.
  A result that would silently translate the object is refused, not clamped or
  shifted. An extrusion Boolean operand may extend below the bed. Negative-source
  results may also extend there, but are excluded from printable bounds/count:
  the positive body's exact native placement check still applies.
- Cancellation leaves the model unchanged but cannot interrupt a single running
  CGAL/OpenMeshCraft call internally. Large solids can therefore take time to stop.

## Geometry and Boolean architecture

- [FaceExtrusion.hpp](<src/libslic3r/FaceExtrusion.hpp>) and
  [FaceExtrusion.cpp](<src/libslic3r/FaceExtrusion.cpp>): non-mutating world-mm
  extrusion kernel. Exact-position welding, patch validation, preserved caps and
  closed boundary walls. No epsilon overlap, plane projection, filling or BRep.
- [GUI tool](<src/slic3r/GUI/Gizmos/GLGizmoFaceExtrude.cpp>): measurement-based
  positive picking and a private owning raycaster for the selected negative,
  independently cached hover, pinned patch/role selection, preview/arrow,
  explicit operation modes and main-thread model application.
- [Shared Model policy](<src/slic3r/GUI/Gizmos/FaceExtrusionModelPolicy.hpp>): non-wx
  role/selection eligibility, immutable placement capture, exact bed decisions and
  actual insert/replace mutations, consumed by the GUI and native regressions.
  Snapshot/Undo, consent, lifetime/cancellation and full Apply remain GUI-owned;
  helper tests alone do not establish interactive acceptance.
- [OpenMeshCraft adapter](<src/libslic3r/MeshBooleanOpenMeshCraft.cpp>): input and
  result validation, with CGAL fallback on copies. Checks include finite and
  representable geometry, indices, degeneracy, manifold/winding, self-intersection,
  solid nesting and operation-specific volume/bounds sanity. Shells are not split
  and individually flipped: an enclosed cavity must retain its negative shell.

A real probe reproduced OpenMeshCraft returning success with invalid geometry for
exactly touching faces, including independently constructed primitive cubes. The
original Union matrix passed 4/13 through that backend versus 13/13 through CGAL;
the corrected adapter passed all 13 cases. This is not claimed to be the exact
popup once observed by the user. See diagnosis and
fix/test report.

Topology and volume/bounds checks are necessary safeguards, **not a mathematical
proof that every nonempty output is exactly the intended set operation**.
The existing non-MSVC signal guard does not catch SIGSEGV/SIGFPE on macOS; no
hardware-fault recovery guarantee is made.

## Tests and local build

```sh
cmake --build build/facetools-arm64 --target face_extrusion_tests --parallel 3
ctest --test-dir build/facetools-arm64 -R '^face_extrusion_tests$' --output-on-failure
```

The dedicated target links real `libslic3r` and `libslic3r_gui`, with the vendored
Catch2 runner, not stubs or a wxWidgets event-loop fixture. The broader historical
upstream suite has unrelated obsolete `load_obj` and `set_rotation` calls and is
not counted as passing.

- V4 isolated new production/test objects linked to the stable normal v3 archives:
  [profiles](<tests/libslic3r/test_planar_profile.cpp>) **12 cases / 6,669 assertions**;
  [cut/stretch](<tests/libslic3r/test_cut_and_stretch.cpp>) **11 / 2,722**;
  [point alignment](<tests/libslic3r/test_point_alignment.cpp>) **7 / 1,695**.
  These preliminary runs were followed by the normal integrated v4 CTest:
  **57 cases / 41,545 assertions PASS**. Full build,
  staging/signature, isolated CLI
  and installed-bundle ledger are local-only evidence.
  The compiled code revision is `ea4714cd7c4a6a757214a882f85e8e2714a87a2f`;
  the installed binary SHA-256 is `71cf7546adeedfdefd478cca3c1bd16c44ebbaa6c5e34b04ce29eac1c2d4bee3`.
  Independent read-only kernel reviews found no concrete blockers. The final GUI
  syntax pass and 59,155 extracted helper/Eigen/shim checks are **not runtime GUI QA**.
- [Kernel tests](<tests/libslic3r/test_face_extrusion.cpp>): 10 cases / 1,196
  assertions, also passed under AddressSanitizer/UndefinedBehaviorSanitizer.
- [Model/Boolean/3MF tests](<tests/libslic3r/test_face_extrusion_model.cpp>): real
  model insertion, immutable source, transformed/signed geometry, cavity and
  invalid/cancel/failure cases, and native Bambu plus legacy 3MF save/reopen.
- The new multiface case creates two perforated outputs from parallel faces at
  different heights on one part, checks effective material inheritance, and
  compares world triangles after native 3MF reload. It covers 11,284 assertions.
- Final normal v3 CTest passed
  **27 cases / 30,459 assertions**, including five new annotation cases.
  Full app build,
  staging/signature checks and
  installed-bundle ledger record delivery.
  These build/log links are **local-only evidence**, deliberately excluded from Git.
- The [versioned native deletion probe](<tools/facetools/README.md>) compiled and
  linked against the completed v3 build: **1,140 scenarios / 5,956 assertions PASS**.
  Its unsorted baseline fails 1,005 scenarios. The original 96-case geometric audit
  also passed; it does not prove correctness for the still-missing user 3MF.
- Eight adapter-specific cases explicitly skip when OpenMeshCraft is disabled.
  Catch2 v2 uses a visible warning/return rather than native SKIP. This is **not**
  MCUT coverage; the backend-OFF configuration was not built/run here.
- V3 extracted helper/shim checks passed 14,457 drag/snap checks, 310 numerical
  ray/camera checks and 2,116 Cut/project-Undo policy checks: 16,883 in total.
  They are not mouse/capture/runtime GUI tests.

The native 3MF fixture initializes dynamic defaults through the normal option
factory, preserving enum serialization maps. Dependency installation must use a
real path without spaces (MPFR resolves symlinks). Dependencies stay project-local.

## Local app and example models

The delivered [Mac app](<../FaceTools Mac/Bambu Studio FaceTools.app>) has physically
copied resources, a separate bundle identifier and a local ad-hoc signature.
Keep it beside its `configuration` folder. Do not replace the official application
or reuse its profile. Updating the bundle must wait until the user saves/closes it
and must leave that configuration and user projects in place.

This is a local Apple Silicon build, not a notarized/portable distribution. It
still uses this Mac's Homebrew zstd/Imath libraries; a deployment target of macOS
11 is not proof that those dependencies support it. The experimental app does not
claim official URL callbacks or default STL/3MF associations.

- Single-ring fixture:
  60×40×4 mm, one through-hole. +5 mm reaches z=9 without changing the source.
- Multiselection fixture:
  two perforated rings at z=2/z=4 joined by a low bridge; one connected closed
  surface, genus 2, 228 triangles, volume 1,875 mm³. Select both ring tops.
- L fixture: 40×10×24 mm,
  20 triangles, volume 3,600 mm³. Underside at Z=20; +10 mm upward extension
  gives 34 mm total height, 4 mm arm thickness and 4,600 mm³ volume.
- Hollow tube: 20×20×30 mm,
  10×10 mm through-hole, 32 triangles, volume 9,000 mm³. A +10 mm stretch
  gives 40 mm height and 12,000 mm³ volume. Both examples come from the
  [deterministic generator](<tools/facetools/generate_stretch_fixtures.py>), which
  checks closure, winding and volume and refuses to overwrite existing files.
- [Developer launcher](<Launch FaceTools.command>) uses its own build-tree profile;
  it is not required for the packaged app's Finder launch.

## Networking and outstanding view issue

The user subsequently reactivated networking work. Live PRE403/public200 checks
confirmed the routing issue; after fresh save/quit consent, only the prototype's
`iot_environment` was changed from 2 to 3. All other configuration bytes and file
permissions were verified unchanged. The user initially deferred installer testing,
then confirmed the plugin was working and supplied a crash report from print sending.
The report lists the plugin as loaded and identifies an off-main-thread dialog resize
while reporting a failed send. V5 queues the error through the owning Job's wx event
handler and joins/cancels that job before dialog reset/teardown. The original send
error code and actual printer outcome are still unknown. No plugin/profile changes
or real print tests were performed for this correction.
The later P2S report confirms cloud mode and successful Send; screenshots show
unsigned-software refusal for direct Print. That caller-authorization issue remains
unresolved and is separate from the historical V5 crash. See the current
verified network/signing facts; no alternative mode or
handoff was selected, and no printer/account/plugin setting was changed for v10.

The earlier read-only diagnosis found the prototype's development profile targeting PRE:
the catalog returned HTTP 403, while the public catalog returned a nominally
compatible Mac plugin. The minimal existing UI switch is **Preferences → Developer
Tools → Host Setting → Product host → debug save button → confirm**. This may log
out an account in the prototype. Keep SSL enabled and use the normal official
installer; no plugin archive, account or profile has been copied by this work.
Those catalog probes alone did not establish installation, ABI or loading; the later
user report establishes loading, not a successful print. See the historical
sanitized routing evidence.

The user's perspective-view defect also occurs outside this tool. No cause or
projection fix has been demonstrated. Restoring orbit input in this gizmo is a
separate correction, not a claim to solve the general perspective issue.

## Manual acceptance still required

### V11 Cut and stretch — seven live scenarios (NOT RUN)

These are the seven original CUT acceptance categories, not seven new native tests.
Installation and wx console tests do not complete them. Preserve the angular filter,
standard Cut and v10 negative-part behavior throughout.

- [ ] **U/Y:** inspect the entire guide-associated branch/leg; only the intended
      connected region moves and excluded regions stay fixed, with no fallback.
- [ ] **T:** inspect and confirm the entire connected bar when both arms are one
      upper component; guides must not invent a CAD branch boundary.
- [ ] **Ambiguous guides:** fixed/unassociated/ambiguous masks refuse without
      silently choosing a majority, nearest region or whole-half scope.
- [ ] **Invalidation:** source/pose/plane/masks/scope/new-plan/clear invalidate
      old region consent and stale work; cancellation/reordering cannot commit it.
- [ ] **Distance:** same-region distance-only edits may retain identity-bound
      consent; zero is no-op, negative extension refuses, and plane/distance stay separate.
- [ ] **First frame:** the FULL orange region is actually drawn/presented and
      inspectable before confirmation; GL cache identity alone is insufficient.
- [ ] **Undo/Redo:** successful Apply and one Undo/Redo preserve the intended
      source pose, fixed remainder and unrelated parts; rejected work changes nothing.

### Additional v11 GUI checks — separate from the seven CUT scenarios

- [ ] Strict-export modal stale state, destination selection, cancellation,
      failure/partial-file reporting and checked output require live GUI acceptance.
- [ ] Native project save/reopen and regenerated slicing require separate checks.

### Earlier acceptance evidence and still-open checks

- [x] User confirmed first-iteration app opening, extrusion and slicing.
- [x] User confirmed hover/highlights and rotation in v2.
- [ ] V3 preview and distance badge update during drag, not just on release.
- [ ] Relative 1 mm drag increments preserve a fractional starting value;
      arbitrary numeric input remains available.
- [ ] Annotation consent retains the selected faces; cancel/failure changes nothing.
- [ ] Hover does not change pinned faces/distance; camera movement refreshes it.
- [ ] Face extrusion: Click/Shift-click toggles the correct same-part faces;
      opposite normals reject here (unlike cut/stretch guide selection).
- [ ] Left orbit, Cmd navigation and right/middle pan preserve face selection.
- [ ] Arrow works at zero and ±distance; Esc/focus/capture loss clean up safely.
- [ ] New part creates N closed parts with one Undo/Redo and no source mutation.
- [ ] Union/Subtract preserve unrelated parts and holes, and handle failure/cancel.
- [ ] Tool close/reload/undo during a calculation cannot apply a stale result.
- [ ] V10: explicitly selected negative singleton activates Faces, isolates the
      cutter and permits correct Click/Shift-click; parent/mixed/multiple negatives,
      parameter modifiers/support and negative Profile/Cut are rejected.
- [ ] V10: signed NewParts remain negative; Union/Subtract change only the cutter;
      positive placement, unrelated parts and source-only annotation consent remain
      correct through one Undo/Redo, native project save/reopen and cavity slicing.
- [ ] V10: switching source/type blocks stale preview/Apply; clear/cancel calculation
      stays isolated, while leaving/reloading restores visibility and lower clipping.
      Check below-bed picking from an unoccluded camera view.
- [ ] Point profile: ordered 3+ vertices, explicit hole, invalid/unfinished loops,
      source lock, undo-last/clear and all three output modes with both signs.
- [ ] Prepare coincidence: fixed/moving order, circle-center Point, real Shift/Plane,
      tiny/zero XYZ, one Undo and unrelated parts/instances unchanged.
- [ ] V5 face-first interaction: select top/lateral guides with different normals,
      then independently set the yellow XY plane at Z=20 (L) / Z=15 (tube).
      No Along/Opposite radio or plane inherited from the first guide.
- [ ] V6: on a later print plate, use the selected L underside itself as the plane
      via **G1/G2/etc. → Use guide face as plane**; the source stays on its plate.
      Repeat with a non-first guide, an oblique face and a rotated/scaled instance.
- [ ] V6: **From a face** can click another triangle of the same selected guide
      without removing it. A typed/dragged offset remains distinct from that face;
      an outer tangent plane rejects rather than inventing a CAD face offset.
- [ ] Plane drag changes only the plane; distance drag changes only translation.
      Esc/focus/capture loss restores the correct parameter; Back keeps guides.
- [ ] After plane confirmation and +10: L height34/arm4/old arm absent;
      tube height40 and hole open; orange preview matches the portion actually moved.
- [ ] Stretch zero/cancel/negative/rejected shelf create no model/history change;
      source visibility returns on clear/off/stale and dragging runs no CSG.
- [ ] Save/reopen 3MF and regenerate slicing for each new operation mode.
- [ ] General perspective issue is retested when a reproducer is available.
- [x] User reports the plugin working; the crash report also lists `libBambuSource`.
- [ ] V5 send-error message appears without crashing; cancel/close/reopen are safe.
      Check the printer/Handy for an existing job before any retry. No automated
      test sends a print, and the original send-failure cause remains unverified.

No completed refinement GUI validation or full CAD support is implied.
