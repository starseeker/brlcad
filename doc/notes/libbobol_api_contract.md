# libBObol API, ownership, and ABI contract

This document is the reviewed contract for libBObol 1.x.  The header lists in
`include/BObol/api_tiers.cmake` are the machine-readable source of truth used
by installation and header-compilation tests.

## API tiers

The **stable source tier** is the host-facing surface included by `BObol.h`:

`BDatabaseSource.h`, `BDefines.h`, `BDisplayEndpoint.h`,
`BDisplaySession.h`, `BExportAction.h`, `BFramebuffer.h`,
`BHeadlessWindowHost.h`, `BHostFactory.h`, `BInit.h`, `BMeasureAction.h`,
`BInput.h`, `BPickDetail.h`, `BRtRender.h`, `BSceneController.h`, `BSceneGroup.h`,
`BSnapAction.h`, `BSourceMeshRequest.h`, `BViewController.h`,
`BViewQuery.h`, `BViewStore.h`, and `BWindowHost.h`.

The **advanced source tier** is installed for BRL-CAD drawing integrations but
is not pulled in by the umbrella header:

`BADC.h`, `BAxes.h`, `BDrawCache.h`, `BEditPreview.h`,
`BEvaluatedPoints.h`, `BGrid.h`, `BHUDLabelOverlay.h`, `BImagePlane.h`,
`BImageSource.h`, `BLineLayerOverlay.h`, `BLodMeshShape.h`,
`BLodRealization.h`, `BLodService.h`, `BLodUpdateAction.h`,
`BMaterialObject.h`, `BMeshLodCache.h`, `BMeshLodSubmitAction.h`,
`BMeshResidencyAction.h`, `BMeshShape.h`, `BPerformance.h`,
`BRealizeAction.h`, `BViewAttachment.h`, `BViewLod.h`,
`BViewportImage.h`, and `BVListShape.h`.

Headers under `src/libBObol`, including realization repositories, performance
implementation helpers, and CAD-assembly build state, are private.  They are
not installed and their symbols must remain hidden.

## ABI policy

The stable C entry points and opaque C handles are the binary compatibility
boundary for libBObol 1.x.  Exported host-facing C++ facades use pImpl where
their Coin runtime contract does not require public fields.  Coin node/action
classes and all advanced C++ APIs are source-supported but their C++ ABI is
experimental: consumers must rebuild them with the matching BRL-CAD release.
An incompatible change to the stable C ABI requires a libBObol SOVERSION bump.
The current `SOVERSION 1` therefore promises the stable C ABI, not a frozen ABI
for every installed C++ class.

The machine-readable `LIBBOBOL_STABLE_C_ABI_HEADERS` list identifies the
headers parsed by the zero-warning documentation gate.  Full reference output
also includes the stable-source and advanced C++ tiers, but undocumented
experimental C++ declarations do not silently expand the SOVERSION promise.

## Owners and lifetime

The display endpoint is the sole attachment boundary.  A host owns or borrows
an endpoint according to the endpoint creation flags; the endpoint owns its
per-view `BObolViewController`.  A GED owns one shared scene controller and
may attach that scene to several endpoints.  Replacing or destroying an
endpoint never transfers scene ownership.

Feature, polygon, compact-instance, and GED scene references are values.  They
contain an owner identity, object identity, and generation or revision.  A
value does not retain its owner.  Every operation resolves all identity fields
before using backing storage, and a removed object, rebuilt store, replaced
controller, or changed scene rejects the stale value.

Pointers returned by `get*`, `find*`, or controller/store accessors are
borrowed.  Unless a declaration explicitly says otherwise, they remain valid
only until the next mutation of that owner and must not be deleted or retained
past owner teardown.

## Supported mutation inventory

This is the finite libBObol 1.x application-mutation inventory. A row names one
state family rather than repeating every overload. Adding a new live mutation
family requires adding a row, naming one owner and adding an acceptance check.
Detached construction, owner-thread publication, worker-result delivery,
derived output, borrowed inspection and observer notification are distinct
lifetime classes.

| State family | Supported entry and owner | Lifetime and revision effect | Failure, callback and caller audit | Status |
|---|---|---|---|---|
| Detached Coin construction | Configure fields and call the type's explicit builder before attachment; examples are `SoBRLViewportImage::rebuildGeometry()` and source realization templates | The caller owns the unattached node. No scene or frame revision is promised | Preparation may leave caller-assigned fields but preserves predecessor generated geometry. Production source constructors are internal to libBObol | Closed |
| Endpoint, session and framebuffer | `BObolDisplaySession`, endpoint/window host methods and `BObolFramebufferStream` configure, ensure, write, present, flush, compose, cursor, reset and close operations | The endpoint owns each per-view attachment; raw stream/framebuffer arguments remain caller-owned and outlive their attachment; the host/stream owns image presentation and delivery revisions | Host publication commits storage, root membership and render work before observers. GED framebuffer and Qt host callers use this boundary | Closed |
| View, camera and render control | `BObolViewController` camera, viewport, interaction, render-request, LoD-policy and lighting entries | One endpoint controller owns view-local state and host-work revisions; shared source state is not rewritten | Owner-thread transitions request later work; callbacks observe complete state. Controller and host regressions cover replacement and teardown | Closed |
| Scene topology and groups | `BObolSceneController` root, group, move, rename, removal and transaction entries | The scene controller owns indexes, membership, structural revision and the corresponding frame revision | Complete successor topology and indexes commit before observers; preparation failure preserves the graph and callback failure leaves it committed. GED draw reducers are the production caller | Closed |
| Database-source identity and configuration | Scene-controller publish, replace, rename, representation, draw-mode, hierarchy and material-policy entries | The controller owns each attached instance; source, input and structural revisions invalidate stamped results | Detached workers receive templates. Attached source fields and direct node configuration are not application entries. `draw_obol*.cpp` callers use scene-controller methods | Closed |
| Database-source display, placement, bounds and metadata | Scene-controller state, display-patch, placement, bounds, draw-metadata and material-refresh entries | Frame revision and every attached peer controller advance at the source commit; geometry remains source-local | Preparation failure preserves the predecessor. Observers see source and owned records together; callback failure leaves the complete change committed | Closed |
| Compact occurrence presentation | `applyPresentationTransaction` and the scene-controller compact display, visibility-frontier, visibility-override, selected-path and selection-delta entries | The source owns retained sparse rules and occurrence revisions; the scene controller owns peer frame effects | Targets are resolved by instance key and retained through notification. GED selection/frontier and `ged_draw_perf` callers no longer mutate borrowed source nodes | Closed 2026-09-20 |
| Source invalidation and realization delivery | Scene-controller stale, stamped merge/certify/adopt/snapshot, realization-state and role/view-policy entries | Invalidation publishes source revision, stale state, owned records and frame revision together. A captured realization stamp gates every late result | Rejected stamps have no live effect. Successful commits precede observers; completed-prefix rules apply to multi-source sweeps | Closed 2026-09-20 |
| External and auxiliary source geometry | Scene-controller external line, point, triangle, annotation, primitive and auxiliary publication entries | The attached source owns immutable geometry and child identity; placement remains a separate source record | Preparation failure preserves preceding geometry. Screen annotations reject the model-coordinate legacy fallback and use native compact publication | Closed |
| LoD service work, capacity and display state | The view controller owns generation/admission, result application, convergence and its user-facing progress classification; `BObolLodService` owns queues, shared producers, reservations, residency, subscriptions, work status and capacity publication | Tasks carry values, generation and source/view identity. Database access is held by explicit leases; result application occurs on the view owner thread. `workStatus()` owns the complete-idle predicate; `generationWorkStatus()` owns coherent generation readiness and includes shared-producer leases. `residentCapacityStatus()` combines the service limit and growth reservation with one directly published stable-byte scalar. Exact stable growth is published before its reservation is released, permitting conservative double counting but no capacity gap. `BObolLodConvergenceStatus::progressDisplayStatus()` derives visibility, terminal readiness and a stable display class from one immutable convergence snapshot. The controller retains notification age, not a second result-availability flag | Batch submission is all-or-none. Cancellation retires a generation without cancelling another generation's shared producer. Generation zero is an empty scope. Unsubscribe is the callback barrier. Phase diagnostics may inspect individual counters, but complete readiness, resident headroom/pressure and progress-display scheduling use the owner status contracts. Libged owns faceplate content and styling; Qt owns publication cadence and event-loop scheduling | Closed for API ownership; lifecycle and load qualification remain in LIFE-01 and CONTROL-01 |
| Retained features | `BObolViewController::features()` returns its `BObolFeatureStore`; store publication, edit and removal methods are the only live feature entries | Store handles contain identity/revision and do not retain the store. The store owns records, nodes, root order and presentation revision | Single and coordinated publication prepare complete successors before commit. GED view, vdraw, annotation/faceplate and Qt edit callers use the store | Closed |
| Retained polygons | `BObolViewController::polygons()` returns its `BObolPolygonStore`; create/update/style/geometry/import/export/remove entries own polygon changes | Store handles are generation-scoped; borrowed geometry/name pointers expire on the documented next mutation | Invalid handles or inputs preserve the record. GED polygon commands use the store; editing qualification remains in EDIT-01 | Closed for API ownership |
| Semantic selection | The view selection store owns path and pick records; GED's selection reducer publishes deltas and then routes compact presentation through scene controllers | Selection records are values. Feature/source handles are resolved at use; presentation has its own revisions | Path-delta publication is bounded and later source streams inherit retained selection. GED selection is the production caller | Closed for API ownership |
| Query, pick, measure, snap and export | Controller/store summaries, value handles, viewport-scoped actions and export helpers | Returned values are snapshots. Returned nodes or pointers are borrowed until the owner's next mutation | Queries do not mutate canonical source state. Viewport actions install and restore camera state for the traversal | Closed |
| GED orchestration | `ged_scene_*` reduction and the Obol backend own command semantics; scene/view/store/controller entries own effects | Reducer values cross the backend boundary; libged does not own a second renderer scene model | Toolkit/C boundaries contain exceptions. Draw, annotate, selection, view-object and deferred-realization callers were included in the source/store audit | Closed for API ownership; FLOW-01 demonstrates the complete workflow on the September 20 Linux software-rendered checkout |
| Application observers | Registered frame, graph, result and store callbacks inspect committed state and queue later work | Callback contexts are borrowed for the call. Removal/unsubscribe is the documented lifetime barrier | Nested mutation is outside the application contract. Throwing observers cannot roll back a committed publication | Closed |

Public Inventor fields remain readable and mutable for Coin integration, and
advanced node methods remain source-supported for detached builders and
library internals. Once a node is attached, raw fields and direct
`SoBRLDatabaseSource` mutators are outside the application contract. Borrowed
query pointers are read-only to applications even when their C++ type is not
`const`.

The September 20 source audit searched production libged, libqtcad, qged,
mged and gtools callers. It found a split `sourceRevision`/`markStale`
transition, GED compact selection/frontier updates, the draw performance
tool's style update and a wireframe fallback invoked on an attached source.
Those paths now publish through the scene controller, including peer frame
revisions before source notification and exact-stamp rejection for the
fallback. Two unused GED inspection branches that could attach raw shape
children were deleted. No production caller in those trees directly writes
an attached source field, attaches a child or invokes the audited source
display, compact-presentation, configuration, realization or geometry
mutators. Tests may still use the node API to construct isolated fixtures.

Internal prepared publication composes every required participant before
notification. Its candidate is not another public live mutation surface.
This keeps the withdrawal of the unsealed September 13 viewport field watcher:
host/stream methods own live viewport updates, and raw writes or connections
do not rebuild host-owned geometry. No stable C ABI changed. Advanced C++
consumers must rebuild with the matching release.

## Thread and callback rules

Scene nodes, controllers, stores, endpoints, and rendering providers are
owner-thread objects.  Their public mutating and query methods run on the host
owner thread.  LoD workers receive database snapshots or plain request values
and return plain result data; they do not mutate live Coin nodes, global type
state, or an OpenGL context.  Source realization has one narrow exception: a
detached, untraversed database-source template with all field sensors disabled
is transferred exclusively to one worker.  The GUI reads an immutable launch
stamp and the occurrence stream while that worker runs; it may inspect or
adopt the detached node only after the job publishes a terminal state with
release semantics.  The owner thread drains and publishes all live scene
changes.

The source registry merges streamed drawing data by occurrence path. A
whole-target overview and a bare primitive can share that path. Once a leaf
record is published, a later priority overview must not replace it, even if
both currently contain bounding-box geometry. Geometry tier alone does not
establish authority: the leaf owns its source request and interaction state.
Later leaf geometry may refine that same record normally.

Per-view CAD presentation derives its selected part, placement and cut from
the view's staged presentation. It must not write those derived records into
the shared source index, which also supplies callback/export compilation.
Complete resets feed an indexed record reader into Obol's prepared replacement;
sparse changes construct only the existing bounded publication batches. The
reader is a synchronous internal conversion, not an application observer, and
does not retain its source or mutate the destination. This preserves one source
baseline across views without a second full occurrence vector.

Source-realization submission is transactional.  The coordinator validates
and prepares the complete request batch before its queue lock commits any
ownership transfer.  Failure leaves every source reference, database handle,
stream, and callback context with the caller; success consumes all source and
database handles in the batch atomically.  A shared callback context remains
alive through the last worker callback even when the client cancels and drops
its job handle.

Frame-request and presentation callbacks may be invoked while the controller
is active, but never after callback removal returns.  A callback must not
destroy the controller from within itself.  It may queue work or request a
subsequent frame.  Publication callbacks receive stack-owned context and must
not retain it.

Application observers inspect committed state and defer subsequent mutations
to a later owner-thread event-loop turn. They must contain their own exceptions
at toolkit/C boundaries. Legacy paths may defensively tolerate nested mutation
or throwing observers; this is not a promise that arbitrary combinations are
supported. No universal deferral queue or reentry guard is implied here.

Prepared publication has two distinct failure phases: preparation failure
preserves the predecessor; observer failure after commit preserves the
successor. Some C++ methods still propagate an observer exception after
finishing publication. Callers must not interpret that as rollback or blindly
retry a non-idempotent operation. The remaining API audit must document each
entry's result and exception boundary before changing its callers.

## Naming and result glossary

`BObol*` names are services and value types.  `SoBRL*` names are Inventor
nodes, actions, details, or elements.  `bobol_*` names are C ABI functions and
opaque C values.  Variables use their role—`endpoint`, `scene`, `source`,
`feature`, and `handle`—rather than a historical renderer layer.

Return values follow these rules:

| Form | Meaning |
|---|---|
| `SbBool` / `bool` | success or a documented predicate; no count is encoded |
| pointer | borrowed/created object as documented; null means unavailable or failure |
| `size_t` getter | a count; it never encodes an error |
| `int` mutator returning `-1/0/1` | error / no change / changed |
| `int` collection operation | negative error, otherwise an explicit count |
| `BRLCAD_OK` / `BRLCAD_ERROR` C adapter | command-style success / failure |

Declarations whose historical integer convention differs must state it at the
declaration.  New stable APIs use a typed status when more than the states in
the table are needed; they do not overload one integer as an undocumented
boolean, count, and error.

## Coordinates and units

Database-source geometry is source-local.  `drawMatrix` and compact occurrence
transforms carry placement separately.  View-controller camera and clipping
methods use BRL-CAD model units; screen coordinates are pixels in the bound
viewport; colors use normalized Inventor channels unless a C declaration says
RGB bytes.  LoD time budgets are microseconds and render timing values are
nanoseconds, as named in their declarations.

`sourceBounds` is also source-local.  Before realization certifies complete
coverage it may be a conservative union derived from the current compact
registry.  Once a realization stream publishes its complete coverage bound,
that value is authoritative across incremental geometry replacement and final
adoption; consumers apply the source placement only when requesting effective
scene bounds.

### Annotation path: owners and representation boundaries

`annotate` is a database authoring command. It creates or replaces
`rt_annot_internal` records; GED events reconcile existing draw intents. The
command does not own a second display list or a client-specific renderer.
The following inventory is part of API-01, with remaining visual acceptance
under QUALITY-01/EDIT-01:

| Operation | Owner and contract | Current evidence or gap |
|---|---|---|
| Creation and update | libged `annotate` writes geometry, styles and update attributes; successful GED events invalidate displayed occurrences | Command tests cover hidden, nested, erased and failed updates |
| Show/hide | GED draw/erase owns intent; source-stream admission cancels only obsolete owners | Delivered-stroke regression covers unrelated pending draws and a subsequent view update |
| Database color | libBObol material sweep resolves the full database path for direct primitives and combination leaves; standard `color` takes precedence over legacy `rgb` | Color and annotation tests cover aliases, inheritance, draw overrides and update round trips; effective color resolution preserves ordinary primitive region identity |
| Model-space plotting | librt `rt_annot_plot` applies the anchor and plane basis once; source compilation retains those coordinates | Duplicate anchor translation is repaired; radial/diameter/angular export assertions check authored spans |
| Autoview | Source bounds supply tight geometry extents; libbv fits their bounding sphere with viewport aspect | Radial scene export checks complete coverage and centered geometry; fitted old cameras align all seven former controls, isolating the intentional framing change before the controls were updated |
| Display-space plotting | Immutable native `PartGeometry::displayPlane` retains anchor and pixel density; camera-specific assembly presentation derives placement without changing source records | OSMesa and System GL tests cover independent cameras, zoom, rotation, resize and nonuniform instance placement |
| Stroke style | The annotation record owns segment width, pattern, color, font and decoration; retained geometry preserves them independently of occurrence selection and draw overrides | Font outlines/italic geometry and compact width/color/pattern runs survive fixed, GLSL, instanced, flat, immediate and direct-software rendering, exact library/Qt export and declared PostScript/PNG/plot mappings |
| Filled areas | librt triangulates annotation loops; immutable native triangle storage has an explicit unlit fill role plus compact authored-color and background-mask runs, visible before strokes in every draw mode | Area/hole export, fill color/replacement and mask pixels pass; PostScript and PNG retain fills while plot documents its area-fill limitation |
| Bounds, picking, export and queries | Native `cadDisplayPlaneTransform` supplies camera-specific bounds, ray picking, rectangle selection, scene collection, viewport export, snapping and measurement | Native pixel/bounds/picking checks, GED rectangle/export checks and viewport snap/measurement checks pass; ordinary action traversal remains deliberately viewless |
| Overview replacement during an edit | Source publication installs complete leaf identity and interaction state; geometry-only refresh applies to established leaves | A temporary root overview cannot be cached as primitive geometry; a focused refresh regression covers its transition to a selectable leaf |

Display-plane offsets use the annotation command's nominal density of
96 pixels per inch divided by 25.4 millimeters per inch; explicit `--dpi` is
already reflected in stored text height. Do not apply it twice. Endpoint
viewport size, camera orientation and projection determine display placement;
they must not be baked into a shared source's model coordinates. Qualification
must include two views with different cameras, resize, zoom, rotation, nested
transforms, and both rendering backends. A passing command or an existing
stroke count does not close these representation gaps.

Authored line presentation uses immutable `WireRep::styleRuns`, indexed by
explicit segment ordinal. A run contains the width multiplier, optional RGBA
and optional 16-bit pattern/factor; the implicit prefix uses occurrence style
and unit width. Widths stay fixed in pixels through camera changes, and the
existing presentation API changes the base width. Admission requires finite
positive scales, valid colors and factors, and strictly increasing in-range
starts. Runs currently apply to explicit non-progressive segments; polylines
and triangle-derived edges retain their occurrence style. No per-segment
occurrence or extra scene owner is created. Ordinary geometry avoids a
per-segment array; authored runs participate in identity and memory accounting.

`SoBRLExportAction::LineRecord` carries a floating-point effective width plus
the exact effective 16-bit pattern, factor and authored-alpha-composed
transparency; `lineStyle` remains the
solid/patterned compatibility classification. Consumers must rebuild against
the additive record change. Source/export regressions cover mixed widths,
colors and patterns, default resets, selection replacement, wire/shaded
realization and live base-width updates. Native pixel checks cover fixed VBO,
GLSL VBO, instanced, flat-batch and direct software execution, including style
replacement under existing part IDs in a reused rendering context. A
compatibility-profile diagnostic override supplies a required immediate-GL
tier-zero run on OSMesa and System GL; it passes the same width, color, pattern,
zoom and replacement assertions while remaining distinct from direct software
wire. Native hardware qualification remains open.

Output consumers use these declared mappings:

- PostScript multiplies effective pixel width by the command's `-l` output
  unit. It converts the exact cyclic 16-bit mask and factor to a PostScript dash
  array, composites transparency over the white page and emits annotation fill
  triangles. Semantic masks fill with page white.
- PNG samples the exact 16-bit mask and factor along each rasterized segment,
  rounds effective width to pixels and composites alpha over the current pixel.
  It rasterizes annotation triangles with one owner for shared edges; semantic
  masks use the command's `-c` background. Ordinary lines, annotation fills and
  annotation strokes are emitted in that order.
- Plot has no arbitrary line-width or area-fill representation. It composites
  color over white and maps continuous, dotted, dashed, center and phantom
  masks to `solid`, `dotted`, `shortdashed`, `dotdashed` and `longdashed`;
  other visible masks use `dotdashed` and a zero mask is omitted.
- Qt object queries return exact effective line and triangle styles in vectors
  parallel to their geometry. Global triangle records also retain their exact
  effective style and semantic-mask flag.

The command regression uses an isolated styled/fill fixture to check
PostScript width, mask, color and fill output, the declared plot mode, PNG line
style and fill effects, and valid cleanup to the selected background. The Qt
export regression checks the parallel style vectors. These mappings close the
known consumer representation gap; visual-format fidelity across clients
remains acceptance work.

Filled annotation areas use the existing `PartGeometry::shaded` storage with
`shadedIsFill`, and share the unlit rendering pass with authored points. They
retain librt's triangulated holes; conversion does not emit compatibility
outlines or triangulation edges as strokes. Fill geometry is non-progressive,
two-sided and ineligible for subpixel replacement. It draws in every mode,
before authored strokes, using a compact per-triangle authored color where
present and the effective occurrence color otherwise. Selection, highlighting
and explicit draw color suppress authored geometry colors. Automatic
wire-mode picking includes fill interiors but excludes hidden ordinary model
triangles. Direct software wire rendering falls back to GL when fills or
authored points require the unlit pass.

`RT_ANNOT_ROLE_MASK` remains a semantic blanking fill rather than receiving an
occurrence or authored color. `FillStyle::backgroundMask` is mutually exclusive
with an authored fill color. The render action publishes the solid or gradient
background already used for the frame, so fixed and GLSL rendering reproduce
that background without storing view state in geometry. Selection does not
replace mask semantics. Exact library and Qt geometry export retain an explicit
mask flag and zero effective transparency; a downstream output supplies its own
background.

The existing geometry store, frame plan and triangle GPU buffers own this
representation. A change to the fill role invalidates the cached plan;
screen-space fills use the same display-plane transform as strokes. Source
summaries and realization reuse retain the annotation's plotted identity.
Library export includes both fill triangles and authored lines with effective
colors and composed transparency. Native tests cover holes, authored fill color and replacement, all four
native draw modes, same-ID role replacement, mixed point/fill/stroke geometry
and camera changes on both software GL backends. Mask tests compare covered and
exposed solid/gradient background pixels in model and display coordinates,
including selection and occurrence-color replacement. Source tests check
fill-only, mixed, colored and mask annotations. Focused consumer tests establish
the mappings above; they do not establish client-wide visual parity.

`SoBRLExportAction::applyViewport`, `SoBRLSnapAction::applyViewport` and
`SoBRLMeasureAction::applyViewport` capture the viewport camera and pixel size
for one traversal. GED export/query/snap consumers and Qt export/measurement
use these entry points. Ordinary `apply(node/path)` keeps the deterministic
viewless XY layout at the stored anchor; reusing an action after viewport
traversal does not retain a camera. Path-local measurement always exposes the
stored display offsets. World-space viewport measurement reports the displayed
geometry, while viewport snapping targets that same geometry.
Mixed fill/stroke annotations still follow the existing measurement/snapping
channel policy: a wire channel takes precedence over triangles. Area and fill
boundary semantics need an explicit decision; triangulation edges are not
authored strokes. The viewport-query repair does not add new fill-boundary
semantics.

The older `SoBRLVListShape` and `BObolExternalAnnotation` contracts use
source-local model coordinates and have no camera-facing display plane.
`publishPrimitiveWireframe` therefore returns `-1` for a screen-space
`rt_annot_internal`; it cannot silently reinterpret display offsets as model
units. Database screen annotations use native compact realization and its
`CadDisplayPlane`. A regression covers successful native realization plus
failure without legacy geometry publication.

The preceding projection candidate is recorded under
`.build-main/obol-qualification/20260919-annotation-projection`. Camera and
viewport changes patch derived annotation placements using the existing
presentation journal; they do not rebuild ordinary geometry or write into
shared instance records. Source bounds contain the model anchor, while
camera-aware bounds include the projected display offsets. These checks do
not close stroke/fill fidelity, native hardware or the full editor matrix.
