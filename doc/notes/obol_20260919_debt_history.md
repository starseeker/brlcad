> Historical snapshot preserved on 2026-09-19. Current instructions and status
> live in [libbobol_active_debt.md](libbobol_active_debt.md). Dated passes and percentage estimates
> below do not qualify the current checkout.

# libBObol active debt

Last reviewed: 2026-09-13

This is the sole remaining-work list for the Obol drawing stack.  It records
work which is not complete; resolved failure analysis belongs in
`libbobol_engineering_lessons.md`, and executable release rows belong in
`obol_production_readiness.md`.  Historical timing directories and
chronological debugging logs are not design authority.

The [simplification guide](obol_simplification_guide.md) defines the bounded
migration and completion gates.  This document owns the tasks which close them.

Any TLC state count retained below is a dated historical observation.  The
sole current formal-suite result and count source is
`tla/baselines/tlc-2.19.json`.

## Resume priorities (2026-09-12)

Start with the [session handoff](obol_session_handoff.md) for checkout and
reproduction facts, and the [simplification guide](obol_simplification_guide.md)
for migration gates S1--S6.  Gate identifiers organize this backlog; they do not
replace its detailed requirements.

**Continue the remaining custom-node/traversal cleanup and lifecycle writer inventory.**
[ROOT-001](libbobol_tla_conformance.md#cxx-root-001--root-and-repository-ownership-retires-before-complete-scene-state)
is qualified and no longer owns active work. Use the
[direct realization writer inventory](libbobol_tla_conformance.md#direct-realization-writer-inventory)
and the classified field/engine, traversal and cleanup work below for the next bounded
source/service operations. GED-005--020 close base promotion,
deferred-proxy retirement,
retained-current marking, complete source-record application and resident-occurrence
adoption, streamed batch publication, terminal deferred adoption, compact proxy bootstrap,
the enclosing retained-current source sequence, cross-scene source copying, path/global
visibility, highlight and transparency publication, and public shape/group-reference
display, deferred-selection target acceptance, intrinsic compact subtract style and composed
database-object rename, redraw/erase, exact mesh-LoD resource attachment and immutable stream
discovery certification/reservation. The split public
adoption operation and zero-caller local GED mesh/vlist
factories were removed instead of acquiring new acceptance surfaces.
The finite GED/source presentation writer inventory has no open rows.

[CONFIG-008](libbobol_tla_conformance.md#cxx-config-008--direct-source-fields-notify-before-dependent-invalidation)
closes all twelve formerly sensor-owned source fields. Direct assignment now reaches the
same prepared invalidation boundary before immediate or delayed observers; connected fields,
engines, quiet notification, reentry, exceptions and allocation denial have permanent
coverage. The private no-op auditors exist only to keep Obol notification live for a quiet
container. Evidence is retained under
`.build/obol-qualification/20260910-direct-field-publication`; the full 80-check state suite
passes. [CXX-SOURCE-026](libbobol_tla_conformance.md#cxx-source-026--throwing-delete-callbacks-interrupt-graph-retirement)
also contains throwing delete callbacks at the native sensor API and proves node, path,
field and parent-owned child retirement.
[CXX-SOURCE-027](libbobol_tla_conformance.md#cxx-source-027--source-clear-exposes-partially-retired-realized-children)
now publishes source realized-geometry clearing as one prepared child/path transaction.
[CXX-SOURCE-028](libbobol_tla_conformance.md#cxx-source-028--axes-and-adc-rebuild-expose-empty-or-partial-overlay-graphs)
closes axes/ADC replacement and visible-to-hidden clearing, including the native
zero-reference-parent lifetime exposed by prepared child edits.
[CXX-SOURCE-029](libbobol_tla_conformance.md#cxx-source-029--grid-rebuild-separates-derived-fields-from-its-hud-graph)
closes grid derived-field/HUD replacement and clearing. Its enclosing
view-input composition is now closed by
[CXX-SOURCE-078](libbobol_tla_conformance.md#cxx-source-078--grid-view-input-and-derived-geometry-publish-one-generation).
Continue the custom-node rebuild family after
[CXX-SOURCE-030](libbobol_tla_conformance.md#cxx-source-030--image-rebuild-clears-presentation-before-stream-preparation),
which closes image-plane/viewport-image replacement, cached pointers and revisions plus the
shared quad builder. Image-source mutation remains separate.
[CXX-SOURCE-031](libbobol_tla_conformance.md#cxx-source-031--hud-and-line-layer-rebuilds-expose-partial-child-graphs)
closes HUD-label and line-layer replacement and clearing.
[CXX-SOURCE-032](libbobol_tla_conformance.md#cxx-source-032--edit-manipulator-rebuilds-expose-partial-recursive-graphs)
closes the axis and indexed manipulator `rebuildGeometry` child transactions;
their public-field and topology/axis setter composition remains open.
[CXX-SOURCE-033](libbobol_tla_conformance.md#cxx-source-033--edit-preview-publishes-children-and-status-in-separate-steps)
closes direct, transformed, failed and cleared preview output; its three input
writers remain separate.
[CXX-SOURCE-034](libbobol_tla_conformance.md#cxx-source-034--navigation-rebuild-clears-cached-hud-state-before-construction)
closes direct navigation-gizmo rebuild publication; camera synchronization and
hover/active input setters remain separate.
[CXX-SOURCE-035](libbobol_tla_conformance.md#cxx-source-035--cutting-affordance-and-scene-lights-rebuild-incrementally)
closes cutting-affordance and scene-light output rebuilds plus stable camera-rig/
scene-light root ordering; their input and enabled setters remain separate.
[CXX-SOURCE-036](libbobol_tla_conformance.md#cxx-source-036--render-environment-repair-publishes-children-and-root-order-separately)
closes live render-environment creation/repair and its coordinated viewport-root
order.
[CXX-SOURCE-037](libbobol_tla_conformance.md#cxx-source-037--feature-overlay-reordering-publishes-one-complete-root-order)
closes pure feature-overlay root reordering, including aliased custom nodes and
retained paths.
[CXX-SOURCE-038](libbobol_tla_conformance.md#cxx-source-038--framebuffer-composition-publishes-viewport-and-layer-roots-together)
closes window-host framebuffer composition across viewport state, realized
children and the three layer roots.
[CXX-SOURCE-039](libbobol_tla_conformance.md#cxx-source-039--raytrace-layer-migration-publishes-one-retained-composition)
closes retained-raytrace display-endpoint layer migration, shares the three-root
transaction with the window host, preserves realized image identity, and avoids
an unrelated raytrace restart.
[CXX-SOURCE-040](libbobol_tla_conformance.md#cxx-source-040--framebuffer-close-publishes-attachment-and-root-retirement-together)
closes explicit framebuffer retirement across the host attachment inventory,
all three layer roots and retained paths for all four composition modes.
[CXX-SOURCE-041](libbobol_tla_conformance.md#cxx-source-041--framebuffer-open-publishes-attachment-and-root-insertion-together)
closes explicit framebuffer image construction and open publication across
scoped node ownership, the host attachment inventory and all layer roots.
[CXX-SOURCE-042](libbobol_tla_conformance.md#cxx-source-042--framebuffer-preparation-precedes-neutral-host-policy-commit)
closes descriptor/open-state rollback for neutral hosts by preparing every
fallible framebuffer resource before the host-policy commit.
[CXX-SOURCE-043](libbobol_tla_conformance.md#cxx-source-043--framebuffer-host-destruction-does-not-run-live-publication)
closes framebuffer-host destruction under allocation denial and callback
failure. Owned controllers retire the graph; borrowed controllers retain a
complete root-owned graph when fallible root preparation cannot run.
[CXX-SOURCE-044](libbobol_tla_conformance.md#cxx-source-044--viewport-size-publishes-one-complete-controller-state)
closes full-window viewport-size publication for ordinary and automatic-LoD
controllers. All viewport copies and the LoD view revision commit before one
host wake; allocation denial preserves the complete old state.
[CXX-SOURCE-045](libbobol_tla_conformance.md#cxx-source-045--both-viewport-setters-share-one-publication-transaction)
closes arbitrary sub-viewport publication and unifies both public viewport
setters behind the same complete-state transaction.
[CXX-SOURCE-046](libbobol_tla_conformance.md#cxx-source-046--established-root-neutral-host-open-publishes-window-policy-and-controller-viewport-together)
closes established-root neutral base-host open composition across the
descriptor, open flag, controller viewport, immediate LoD revision,
progressive-work state and final render request.
[CXX-SOURCE-047](libbobol_tla_conformance.md#cxx-source-047--host-owned-root-exists-before-first-open-can-publish)
closes host-owned first-root construction and first-open publication. A new
host is exposed with one group root, three matching viewport copies and no
pending construction work.
[CXX-SOURCE-048](libbobol_tla_conformance.md#cxx-source-048--default-and-borrowed-controllers-satisfy-one-host-attachment-invariant)
closes default-controller construction, borrowed-controller root validation,
allocation-free identical attachment, C endpoint construction containment and
allocation-free LoD-service detachment during controller destruction.
[CXX-SOURCE-049](libbobol_tla_conformance.md#cxx-source-049--headless-camera-state-exists-before-host-lifecycle-publication)
closes headless camera creation and adoption. The default controller owns one
private viewport camera during construction; headless open validates that
invariant and never searches or changes modeled scene content.
[CXX-SOURCE-050](libbobol_tla_conformance.md#cxx-source-050--headless-provider-and-host-policy-publish-as-one-lifecycle-transition)
closes headless context-provider publication for open, live replacement and
close. Provider registration, renderer invalidation, descriptor, viewport and
host state prepare before one final callback; allocation failure retains the
preceding complete policy. Native Obol provider and cache-key registry changes
now preserve their old entry until replacement registration succeeds.
[CXX-SOURCE-051](libbobol_tla_conformance.md#cxx-source-051--feature-record-node-and-root-publish-as-one-replacement)
closes existing-feature node replacement and reparenting. A copied record,
replacement node and every affected root order prepare before one commit;
aliased custom nodes remain attached until their final owner moves or retires.
All existing-record feature publishers which rebuild retained nodes now use
that publication primitive.
[CXX-SOURCE-052](libbobol_tla_conformance.md#cxx-source-052--feature-controller-migration-publishes-all-roots-and-ownership-together)
closes multi-record feature-controller replacement, detachment and adoption.
All record placements and affected root orders prepare before controller,
attachment and presentation ownership publish together; old and new
controllers are then both awakened without stopping after the first callback
exception.
[CXX-SOURCE-053](libbobol_tla_conformance.md#cxx-source-053--qt-canvas-validation-precedes-host-identity-and-ownership-changes)
closes invalid Qt canvas/controller rebinding. The adapter validates the
candidate controller before changing or destroying its canvas, reports
acceptance, preserves borrowed and owned predecessors on rejection, and uses a
scope-bound factory-bind preservation flag. Valid live rebind under
framebuffer-close or canvas-resource failure remains a distinct S4 lifecycle
boundary.
[CXX-SOURCE-054](libbobol_tla_conformance.md#cxx-source-054--qt-fallback-faceplate-nodes-prepare-before-live-publication)
closes per-node grid, axes and ADC publication in the direct no-GED Qt
fallback. Each complete successor is built off graph and installed through one
prepared child-list replacement; invalidation-free view/frame no-ops retain
node identity. The public input setters remain separate.
[CXX-SOURCE-055](libbobol_tla_conformance.md#cxx-source-055--initial-feature-publication-commits-identity-graph-and-frame-obligation-together)
closes first-record feature publication. Candidate payload, both indexes,
retained node/root order, presentation identity and controller request prepare
before one allocation-free commit. Failure preserves the empty predecessor and
its next identity; callbacks see the complete successor.
[CXX-SOURCE-056](libbobol_tla_conformance.md#cxx-source-056--ged-faceplate-records-publish-as-one-logical-composition)
closes whole-faceplate GED feature composition. Mixed line, HUD-label and
custom-node replacement/removal prepares every record, index, retained node,
root order and one render request before one allocation-free commit. The GED
adapter now submits center dot, interactive rectangle, grid, ADC, parameters,
LoD progress, scale and axes as one logical publication. Framebuffer compose/
present still follows that transaction, and the general public
`ged_view_feature_batch` path remains sequential.
[CXX-SOURCE-057](libbobol_tla_conformance.md#cxx-source-057--feature-removal-and-clear-commit-before-observation)
closes public handle/name/owned-name, prefix, scope and clear retirement.
Matching sets, final root orders, result notifications and one controller
request prepare before an allocation-free commit; clear also resets command-
owner generations inside that commit. Callback failure leaves the complete
retirement visible. Destructor-time feature teardown remains separate.
[CXX-SOURCE-058](libbobol_tla_conformance.md#cxx-source-058--feature-controller-migration-includes-both-render-obligations)
closes the render-obligation remainder of feature-controller migration. Old
detach and new attach requests and their LoD transition scopes now prepare
with all affected roots and commit with controller and attachment ownership;
graph observers see both standing obligations and callback failure does not
stop later notifications.
[CXX-SOURCE-059](libbobol_tla_conformance.md#cxx-source-059--direct-feature-records-and-mesh-fields-prepare-before-publication)
closes the remaining direct feature-record and retained-mesh field writers.
Complete record, point/normal or selection arrays, source revision, render
request and result payload prepare before one allocation-free commit. The new
native `SoMField` value-replacement primitive preserves field and node identity,
restores the preceding notification state and reports the completed successor
once. Metadata-only edits prepare their record and result without a synthetic
frame effect.
[CXX-SOURCE-060](libbobol_tla_conformance.md#cxx-source-060--ged-edit-transactions-publish-one-complete-retained-preview)
closes GED edit-preview transaction composition. Primitive plotting finishes
before one complete edit-preview publication installs geometry, style,
identity, intent, revisions, overlay metadata and a standing render request.
Identical updates are no-ops, terminal events retire directly, and the
uncomposed feature-store touch API is removed.
[CXX-SOURCE-061](libbobol_tla_conformance.md#cxx-source-061--initial-custom-nodes-publish-configured-overlay-state-in-one-commit)
closes initial custom-node composition in the qged ARB, ellipsoid, sketch and
BoT presenters and the libged navigation gizmo. Each detached, fully configured
node now publishes with its typed overlay, root order and frame obligation in
one feature transaction. The focused allocation and observer test passes 75
positions; exact predecessor/current evidence is retained under
`.build/obol-qualification/20260912-custom-node-overlay-publication`.
[CXX-SOURCE-062](libbobol_tla_conformance.md#cxx-source-062--edit-manipulator-setters-publish-one-complete-node-state)
closes node-local edit-manipulator setter atomicity. Axis and indexed topology,
selection and visibility changes prepare owned successor children before one
private/scalar/child commit and notify only complete state. The focused sweep
passes 767 allocation positions and rejects the exact CXX-SOURCE-061 library.
[CXX-SOURCE-063](libbobol_tla_conformance.md#cxx-source-063--qged-manipulator-updates-publish-one-store-visible-successor)
closes that enclosing feature-store composition for the small qged ARB,
ellipsoid and sketch manipulators. Each live geometry, selection or pointer
change now replaces one complete detached node while advancing the record and
publishing its frame obligation. Exact preceding plugins fail the new retained-
revision assertion; current single/quad replays pass. Plugin state retains only
stable feature identity and resolves the current node on use.
[CXX-SOURCE-064](libbobol_tla_conformance.md#cxx-source-064--bot-shared-mesh-edits-publish-coordinated-surface-successors)
closes the BoT surface half of shared-mesh point moves. All per-view wrapper,
record, ordered-root and frame successors prepare before an allocation-free
one-to-three-point shared-geometry patch; every store commits before callbacks.
The 128-position coordinated sweep preserves both stores and the heavy geometry
on preparation failure. Exact preceding qged keeps the old wrapper; current
single/quad replays replace it while retaining the same shared mesh. The GED
point-handle batch and cross-controller feature composition remain separate.
[CXX-SOURCE-065](libbobol_tla_conformance.md#cxx-source-065--navigation-gizmo-state-publishes-as-retained-successors)
closes semantic navigation-gizmo mutation. Style, visibility, interaction and
camera orientation now form one private value snapshot and one detached custom-
node successor; record, ordered root and frame obligation commit before
callbacks. The fixed-camera node API prevents later sensor/traversal mutation,
and the focused sweep passes 1,941 allocation positions. Viewport anchoring
remains render-local layout.
[CXX-SOURCE-066](libbobol_tla_conformance.md#cxx-source-066--live-camera-replacement-publishes-one-controller-state)
closes direct live camera replacement. Native viewport child/path state and
the viewport, controller and render-manager camera pointers now commit with the
LoD revision and frame obligation before callbacks. Replacement, removal and
addition pass 181 fresh-process allocation positions. The navigation plugin
retargets its camera field observer at the completed viewport-root edge, so a
retired camera cannot publish another gizmo generation.
[CXX-SOURCE-067](libbobol_tla_conformance.md#cxx-source-067--view-camera-input-publishes-one-derived-controller-state)
closes the enclosing `syncCameraFromViewContext()` input transaction. Initial
viewport adoption, same-generation camera fields, projection replacement,
tracked lighting, clipping, section-aid geometry, LoD state and one frame
obligation now prepare and publish together. Its ten ordinary/automatic
operation cases pass 1,578 fresh allocation positions; exact equal input
allocates and publishes nothing. The navigation plugin now also publishes one
retained successor for a camera-generation change whose orientation is equal.
[CXX-SOURCE-068](libbobol_tla_conformance.md#cxx-source-068--cutting-plane-input-publishes-intent-clip-and-aid-together)
closes the public cutting-plane setter transaction. Private intent, retained
clip fields, visible section aid and one capacity render request now prepare
and publish together. Eight ordinary/automatic operation cases pass 711 fresh
allocation positions; invalid and exact inputs allocate and publish nothing,
and callback reentry retains its later successor.
[CXX-SOURCE-069](libbobol_tla_conformance.md#cxx-source-069--scene-light-input-publishes-values-nodes-and-policy-together)
closes scene-light content and enablement publication. Copied input, retained
children or enablement fields and one capacity request now commit together;
enablement retains node identity and exact input does no work. A combined
public operation removes the preceding two-call pair from the GED lighting
synchronizer. Twelve operations in both ordinary and automatic modes pass 443
fresh allocation positions.
[CXX-SOURCE-070](libbobol_tla_conformance.md#cxx-source-070--complete-lighting-input-publishes-one-controller-generation)
closes the remaining `ged_view_lighting_sync()` sequence. Profile, normalized
offset, camera tracking, headlight policy, environment and camera-light
fields, scene-light content/enablement and one capacity request now prepare
and commit as one controller generation. Existing leaf setters share the same
owner; exact and invalid input does no work. The fresh complete-state sweep
passes 109 allocation positions across ordinary and automatic LoD, and the
exact full compact matrix passes 165 direct checks and 1,075 PASS rows.
[CXX-SOURCE-071](libbobol_tla_conformance.md#cxx-source-071--master-lighting-publishes-shading-rig-and-renderer-epoch-together)
closes the master shading input. Flat/PHONG state, the profile-dependent camera
rig, renderer timing/capacity invalidation and one typed frame now prepare and
commit together. The Studio/MGED ordinary/automatic sweep passes 120 fresh
allocation positions; exact input preserves renderer evidence and does no
work, and the immediate predecessor fails the behavioral contract.
[CXX-SOURCE-072](libbobol_tla_conformance.md#cxx-source-072--scalar-appearance-inputs-publish-retained-and-frame-state-together)
closes background, depth-test, headlight color/intensity and depth-cue input
publication. One shared prepared transaction commits private and Coin fields,
derived fog state, renderer invalidation where applicable and one typed frame
before callbacks. The seven-operation ordinary/automatic sweep passes 128
fresh allocation positions; exact and invalid input does no work, stale fog
state is repaired, and the immediate predecessor fails the behavioral
contract.
[CXX-SOURCE-073](libbobol_tla_conformance.md#cxx-source-073--render-action-appearance-has-one-redraw-and-renderer-epoch)
closes transparency and antialiasing publication. Controller intent, onscreen
and cached-offscreen actions, renderer invalidation and one reason-matched
frame now commit before the callback; exact input also repairs stale backend
state. All eight change/repair and ordinary/automatic operation states are
allocation-free, and the immediate predecessor fails the behavioral contract.
[CXX-SOURCE-074](libbobol_tla_conformance.md#cxx-source-074--clip-bounds-publish-derived-planes-and-renderer-epoch-together)
closes camera-relative clip-bound publication. Private bounds, both derived
retained planes, renderer invalidation and one reason-matched frame now commit
before callbacks. Active, disabled, one-sided and stale-repair cases pass 104
fresh allocation positions across ordinary and automatic LoD; exact and
invalid input does no work, and the immediate predecessor fails the behavioral
contract.
[CXX-SOURCE-075](libbobol_tla_conformance.md#cxx-source-075--software-wire-policy-has-one-complete-renderer-generation)
closes the renderer-appearance family. Controller, LoD-wrapper and compact-
batch software-wire policy, renderer evidence and one reason-matched frame now
commit before callbacks. Exact checks repair stale render paths; absent-LoD,
invalid-input, callback-failure and reentry cases are covered. Seven operations
in ordinary and automatic LoD pass 56 fresh allocation positions, and the
immediate predecessor fails the complete-state contract.
[CXX-SOURCE-076](libbobol_tla_conformance.md#cxx-source-076--controller-edit-preview-input-publishes-one-retained-generation)
closes the direct controller edit-preview input. A complete detached candidate
prepares before the existing preview's scalar fields and child are committed
with its capacity-classified frame; insertion and removal use prepared root
edits. Existing preview identity is preserved, internal stale sensors are
handled without suppressing external observers, and graph/field/frame failures
drain. Insert/replace/remove in both LoD modes pass 364 stable allocation
positions, and the immediate predecessor fails the complete-state contract.
[CXX-SOURCE-077](libbobol_tla_conformance.md#cxx-source-077--controller-overlay-inputs-publish-one-complete-graph-and-frame)
closes the direct controller line-layer/HUD-label input family. Complete
detached candidates, prepared root edits and one capacity-classified frame now
commit before observers. HUD replacement preserves its existing node identity
while scalar fields, child geometry and cached label pointer publish together;
line replacement preserves its established new-node semantics. Line operations
pass 348 stable allocation positions and HUD operations pass 538 across both
LoD modes. The immediate predecessor fails both complete-state selectors.
[CXX-SOURCE-078](libbobol_tla_conformance.md#cxx-source-078--grid-view-input-and-derived-geometry-publish-one-generation)
closes `bobol_grid_configure_from_view()` and its context overload. A detached
grid candidate now derives the whole successor before prepared scalar fields
and child geometry preserve the attached grid identity in one commit. Visible,
hidden and snap-only configurations pass 9,422 stable allocation positions;
the immediate predecessor fails the complete-state selector.
[CXX-SOURCE-079](libbobol_tla_conformance.md#cxx-source-079--axes-and-adc-view-inputs-publish-one-retained-generation)
closes composition for the view fields represented by `SoBRLAxes` and
`SoBRLADC`. Two explicit state APIs derive detached candidates and use one
shared scalar/child publisher while preserving the retained node identity.
Axes and ADC visible/hidden operations pass 222 stable allocation positions;
the immediate predecessor exposes a partial input generation. Rich axes
labels, ticks and style remain capability scope rather than publication debt.
[CXX-SOURCE-080](libbobol_tla_conformance.md#cxx-source-080--navigation-gizmo-inputs-publish-one-retained-generation)
closes the supported navigation gizmo method boundary. Hover, active, camera
add/replace/remove, fixed snapshots and traversal-driven rotation now stage a
complete detached HUD and commit scalar fields, child geometry, private
pointers and camera ownership before callbacks. Realized and unrealized camera
operations cover 5,172 stable allocation positions; the immediate predecessor
exposes a partial generation. Direct style-field writes and renderer-derived
anchor placement remain separately classified.
[CXX-SOURCE-081](libbobol_tla_conformance.md#cxx-source-081--image-source-inputs-publish-one-retained-generation)
closes `SoBRLImageSource` stream, image, clear and refresh publication. Stream
subscription and owned-image conversion now prepare off target; one quiet
commit installs the complete metadata and callback ownership before retiring
the predecessor and notifying observers. Eight operations cover 174 stable
allocation positions, including clean replacement after dirty metadata; owned
self-reattachment is an exact no-op. The immediate predecessor exposes a
partial generation.
The private ownership refactor preserves the preceding 928-byte C++ object
layout. Image display rebuilding remains a distinct consumer step.
[CXX-SOURCE-082](libbobol_tla_conformance.md#cxx-source-082--framebuffer-cursor-input-publishes-one-presentation-generation)
closes the three framebuffer cursor-state entry points. A fixed three-field
publisher commits visibility, image position, shape and the standing
presentation request before callbacks. Four ordinary/automatic operations pass
44 stable allocation positions; the immediate predecessor exposes a partial
cursor generation. Exact input is allocation- and notification-free.
[CXX-SOURCE-083](libbobol_tla_conformance.md#cxx-source-083--framebuffer-view-input-publishes-one-retained-presentation-generation)
closes the framebuffer view method. A staged viewport candidate rebuilds the
HUD from accepted texture bytes and revisions, then commits center, zoom,
geometry, cache pointers and the presentation request before callbacks. Pending
stream data remains owned by flush. Pan/zoom operations cover 1,322 stable
allocation positions; the immediate predecessor exposes a partial generation.
Exact input is allocation- and notification-free.
[CXX-SOURCE-084](libbobol_tla_conformance.md#cxx-source-084--framebuffer-reset-publishes-one-retained-presentation-generation)
closes framebuffer reset. The retained viewport transaction commits center,
zoom, cursor visibility, HUD geometry, cache pointers and the presentation
request together without consuming pending source pixels. Ordinary/automatic
reset covers 442 stable allocation positions; the immediate predecessor exposes
a partial generation. Exact repeated reset does no work.
[CXX-SOURCE-085](libbobol_tla_conformance.md#cxx-source-085--framebuffer-viewport-input-composes-placement-controller-size-and-frame-obligation)
closes framebuffer viewport placement. Retained geometry prepares before the
controller size transaction; placement, geometry, controller region/LoD state
and one capacity-aware request commit before callbacks. Four ordinary/automatic
operation shapes cover 1,776 stable allocation positions; the immediate
predecessor exposes a partial composed generation. Exact input does no work.
[CXX-SOURCE-086](libbobol_tla_conformance.md#cxx-source-086--framebuffer-flush-publishes-one-source-to-display-generation)
closes framebuffer source-to-display flush. The existing private source
publisher, stable stream payload, detached viewport successor and presentation
request now form one transaction. Dirty-pixel ordinary/automatic coverage
passes 486 stable allocation positions; the immediate predecessor exposes a
partial source/display generation. Settled flush does no work.
[CXX-SOURCE-087](libbobol_tla_conformance.md#cxx-source-087--retained-rt-image-publication-no-longer-exposes-an-empty-or-split-generation)
closes retained RT first-image and restart publication. Scoped off-scene
resources and one source/payload/viewport/root/request transaction prevent an
empty attachment or split source/display generation. Exact staged restart
pixels survive failed preparation, and restart generations make callback-time
successors retire older worker launches. Ordinary/automatic restart passes 538
stable allocation positions; two observed activations are complete, and the
immediate predecessor exposes an empty first viewport. Exact RT policy input
does no work.
[CXX-SOURCE-088](libbobol_tla_conformance.md#cxx-source-088--renderer-selection-publishes-one-complete-endpoint-generation)
closes renderer-engine selection. Controller invalidation and outgoing root
replacement prepare before mutation; RT activation commits through the image
publisher's precommit hook. Retiring RT state remains alive until all callbacks
drain, and failed postcommit render-manager/request effects remain retryable.
Graphical frame notification waits until render-manager synchronization has
succeeded.
Twelve transitions in ordinary/automatic LoD cover 2,887 stable allocation
positions; the exact predecessor exposes a partial engine/root generation.
Settled repeated policy does no work.
[CXX-SOURCE-089](libbobol_tla_conformance.md#cxx-source-089--endpoint-destruction-is-terminal-no-throw-and-callback-safe)
closes endpoint terminal teardown. Destruction uses a dedicated no-throw path
that retires every controller work owner even when presentation observers or
automatic-LoD service cancellation fail. Callback-time destroy defers final
deletion until the synchronous endpoint operation drains, while further
mutation rejects the terminal endpoint. Host detachment, RT worker/source
ownership, borrowed-controller root restoration and controller-owned
render-manager retry state are closed by the same boundary. The focused matrix
covers 30 owned/borrowed engine allocation positions and 9 safe borrowed-root
fallbacks; the exact predecessor terminates with status 42. Source and full
compact suites, focused ASan, adjacent publication/CTest rows, independent
compilation and public contracts pass.
Continue with the remaining classified traversal/lifecycle writers. The next
focused boundary is direct public viewport field writes. Classify Qt widget
policy, borrowed-controller construction failures and other constructor-only
cleanup separately. Continue the remaining service-lifecycle and
resource-magnitude rows without reopening qualified field, delete-callback or
source-clear behavior.

[GED-020](libbobol_tla_conformance.md#cxx-ged-020--stream-discovery-facts-and-reservation-bypass-source-acceptance)
makes expected count and optional source profile one exact immutable discovery contract.
The producer publishes its final validated facts once, and the scene controller accepts them
only for the current source key and realization stamp. Capacity reservation occurs inside
that same stamped batch delivery and preserves preceding geometry on failure. Resource
identity invalidation revokes the contract; view-only invalidation and failed delivery retain
valid discovery evidence. The three split public source writers are removed.
Nine focused checks pass in 26.22 seconds under
`.build/obol-qualification/20260910-ged-stream-contract`, including allocation denial,
stream cancellation, LoD append traversal, full GED draw sync and its CTest row. The saved
GED-019 library supplies the split-symbol baseline. This closes the row without changing the
dated stage estimates.

[GED-019](libbobol_tla_conformance.md#cxx-ged-019--mesh-lod-handle-and-bounds-attach-through-split-path-writers)
keeps a prepared mesh-LoD reader in caller ownership until an exact source key and
realization stamp accept the complete reader-and-bounds record. Non-finite or inverted bounds,
stale inputs and replaced sources reject without changing source ownership. Source, input and
database invalidation retire the record; view-only invalidation retains it. Mesh attachment is
private service state and does not publish a frame. Split source/path setters and dead runtime
fields are removed.

[GED-018](libbobol_tla_conformance.md#cxx-ged-018--redraw-and-erase-use-live-order-global-realization-and-duplicate-visibility-writers)
captures all redraw/erase targets before the first callback. Exact and prefix erase retire
source/group edges through one scene transaction. Scoped redraw realizes only its exact
stamped source; descendant cleanup stays beneath that owner. Retained compact visibility is
a durable presentation override with an optional pre-removal stamp, so later callback edits
and replacements win. Pure removals no longer start unrelated pending realization.

[GED-017](libbobol_tla_conformance.md#cxx-ged-017--database-rename-publishes-source-compact-group-and-repository-state-in-pieces)
prepares repository residency, source identity and hierarchy, compact semantics and
renderer hierarchy, group paths/names/intents, indexes and peer revisions before one
database-rename commit. Path-derived GED keys retarget while opaque keys stay stable.
Callbacks see the complete scene, and later edits, removal and exceptions remain current.

[GED-016](libbobol_tla_conformance.md#cxx-ged-016--streamed-subtract-style-is-missing-from-realization-publication)
derives the default dashed style from the durable compact boolean operation, so direct,
streamed and reconstructed occurrences publish identity and style together. The later
GED redraw setter was redundant for eager work and ineffective for streaming; it is
removed. The compact setter remains available for explicit later overrides.

[GED-015](libbobol_tla_conformance.md#cxx-ged-015--deferred-selection-callback-can-supersede-its-source-target)
retains the exact source while compact selection initialization notifies observers, then
requires its key, node, routing identity, path and representation to remain current before
capturing worker state. Callback retarget, same-key replacement and removal survive and
prevent work from starting for the superseded target.

[GED-014](libbobol_tla_conformance.md#cxx-ged-014--public-source-display-publishes-configuration-and-appearance-in-pieces)
snapshots all same-path sources before publication and applies draw mode, representation,
display and material changes through one complete source state. The public shape-reference
path canonicalizes its database-path alias before that publication, so it cannot publish
the same source twice and overwrite a callback edit. Pending source edits, same-key
replacement and removal supersede the saved work. The public group-reference branch keeps
its already prepared single group publication.

[GED-013](libbobol_tla_conformance.md#cxx-ged-013--presentation-commands-publish-compact-and-aggregate-state-separately)
routes path/global visibility, highlight and transparency through one scene transaction.
Each source commits aggregate display, retained compact rules, effective occurrence and
renderer state, and its frame effect before observers. Multiple sources and groups form a
consistent prefix, and callback edits, removals or replacements supersede pending work.
The old direct GED writer loops are removed.

[GED-012](libbobol_tla_conformance.md#cxx-ged-012--cross-scene-source-copy-applies-material-policy-later)
copies material policy with the rest of a primary source summary in one secondary-scene
publication. Root callbacks see the complete inserted source, and edits, removal or a
same-key replacement cannot be changed by a later material-policy lookup.

[GED-011](libbobol_tla_conformance.md#cxx-ged-011--retained-current-source-policy-publishes-before-current-state)
includes current external realization state in the source publisher which changes its
policy. The live GED regression now exercises that enclosing path: callbacks see the new
policy and current geometry together, and edits, removal or a deliberately stale same-key
replacement survive. The obsolete follow-up helper and its private declaration are removed.

[GED-010](libbobol_tla_conformance.md#cxx-ged-010--compact-proxy-bootstrap-publishes-one-source-snapshot)
routes cached-manifest, structural and view-envelope compact bootstrap through one stamped
scene-controller operation. It constructs the registry off-scene and commits occurrence
identity, certified bounds/profile, terminal realization state, primary-child retirement
and frame effects before observers. A commit witness lets GED retain a complete snapshot
when an observer throws `std::bad_alloc`. The live adapter regression observes the compact
replacement and retains a later callback edit; the redundant post-line-set AABB bounds
publication is removed because external shape publication already owns that exact fact.

[GED-009](libbobol_tla_conformance.md#cxx-ged-009--deferred-result-publication-uses-stamped-scene-acceptance)
routes every live streamed occurrence batch and the terminal detached result through
the stamped scene controller. Source state, source-local child/index effects and frame
revisions commit before callbacks. Completion witnesses distinguish an internal
allocation failure from an observer throwing after a complete batch, bound or terminal
commit. GED re-resolves the exact routing identity after callbacks and captures proxy
descendant identities before terminal notification, so later callback replacement,
removal or reparenting is not overwritten by cleanup.

[GED-008](libbobol_tla_conformance.md#cxx-ged-008--resident-occurrence-adoption-notifies-before-donor-retirement)
now prepares a detached target compact-index candidate together with donor edge
retirement, source indexes and controller revisions. Those effects commit before
source or hierarchy observers. GED groups donors by stable instance key and resolves
the retained objects inside the controller transaction, so callbacks may edit,
replace or remove the target without invalidating a cached pointer. Allocation
failure preserves the preceding scene, and the old payload-only public adoption
entry point was removed after its callers moved to the composed owner.

[GED-007](libbobol_tla_conformance.md#cxx-ged-007--source-record-application-publishes-seven-field-families-separately)
now translates the GED source record once and publishes configuration, display,
material, realization, ownership and view policy through the complete source
transaction. The explicit realization record is applied after configuration/view
invalidation on the private candidate, so an authoritative current record remains
current. Live callbacks see one frame and may edit, replace or remove the current
source. The saved preceding libraries reproduce the partial callback state.

[GED-006](libbobol_tla_conformance.md#cxx-ged-006--deferred-proxy-state-and-ownership-publish-separately)
now publishes external-to-native invalidation and stale-to-current realization with
their role ownership in one source transaction. Retained proxy geometry remains as
continuity data while marked stale, and current callbacks may edit, replace or remove
the source without a later write reaching the new target. The preceding libraries
reproduce both partial-state failures; current controller and live GED regressions pass.

[GED-005](libbobol_tla_conformance.md#cxx-ged-005--base-source-promotion-publishes-identity-before-representation)
publishes current-to-next source identity, representation and source/group state as
one composed operation. The same callback sees the new key and shaded representation,
and its later display edit survives. The adjacent private local shape audit found no
callers or dynamic exports and removed about 320 lines of unsafe factory code.

[GROUP-008](libbobol_tla_conformance.md#cxx-group-008--shape-state-notifies-before-complete-fields-and-frame-effects)
now stages the four shape-state setters' requested fields and frame effects before
observers. Shared geometry, existing tolerances and original-target lifetime remain
intact; later callback edits survive. The duplicate live writers and dispatch wrappers
are removed. Root/repository ownership is qualified; the remaining enclosing GED
occurrence/adoption callers stay open.

[GROUP-007](libbobol_tla_conformance.md#cxx-group-007--shape-membership-notifies-before-complete-scene-effects)
now composes shape movement/removal through prepared child publication. Both parent edits,
child paths and scene effects publish before observers, retaining the original shape and
shared geometry through callback removal, replacement or movement. Lookup still excludes
database-source subtrees.

[GROUP-006](libbobol_tla_conformance.md#cxx-group-006--source-removal-notifies-before-complete-scene-effects)
now composes source removal and source-wide clearing through the same child publisher.
Shared source edges, auxiliary subtrees, complete descendant indexes and repository
ownership survive correctly. Path removal retains its selected identity when instance
keys collide; duplicate child paths retain the selected edge position. Later callbacks
remain current, and the obsolete live removal/index writers are removed.

[GROUP-005](libbobol_tla_conformance.md#cxx-group-005--recursive-group-removal-notifies-before-complete-scene-effects)
now composes subtree removal and whole-group clearing through the existing child
publisher. Shared cache ownership and remaining graph edges survive; path/node observers
see complete scene effects. The duplicate live removal/clearing sequences are removed.
Repository attachment/retirement remains a separate task.

[GROUP-004](libbobol_tla_conformance.md#cxx-group-004--child-membership-notifies-before-scene-and-repository-effects)
now prepares child append/removal with complete subtree indexes, repository ownership and
one structural/frame effect before observers. Shared nodes retain their surviving edges
and lookups; repository seeding preserves published geometry metadata and prepares cache
ownership before commit. Later callback changes survive. Root replacement and other
enclosing scene writers still need their own publication audits.

[GROUP-003](libbobol_tla_conformance.md#cxx-group-003--group-rename-notifies-before-index-and-revision-effects)
now prepares complete name/path/index/revision publication through the existing group
owner. It preserves source state and later callback changes, removes the recursive live
path writer and repairs native name-registration allocation failure. Its production,
native and graphical checks qualify this boundary.

[GROUP-002](libbobol_tla_conformance.md#cxx-group-002--source-movement-publishes-intermediate-membership)
now composes source movement through the existing prepared publisher, preserving source
configuration/geometry and later callback edits. Both movement and composed source publication
reject cycles through newly prepared descendants. GROUP-001's complete group creation remains
qualified. These repairs remove the separate live creation/movement sequences.

[GED-004](libbobol_tla_conformance.md#cxx-ged-004--group-synchronization-overwrites-a-source-callbacks-edit)
now composes group metadata with source/scene publication before observers. The standalone
group intent/display setters reuse the same prepared scalar helpers and publish their
frame effect before callbacks. New/existing groups, both callback directions, allocation,
quiet fields, no-ops, group-only effects, batches and lifetime/exception behavior have
production regressions. GED's later group-sync/regroup writer is removed (50 net lines).

[GED-003](libbobol_tla_conformance.md#cxx-ged-003--redraw-resets-and-restores-current-appearance)
removes the root redraw's default-style reset and saved-style restoration. Current source
appearance now goes directly into the existing publisher; callbacks can edit, replace or
remove that source without a later stale appearance write. The repair removes 118 net lines.
CONFIG-007's prepared source and scene/index publication remains qualified.

[GED-002](libbobol_tla_conformance.md#cxx-ged-002--unreachable-leaf-redraw-machinery)
removes the unreachable leaf metadata/style adapter and its traversal, along with the
same helper's unused deferred branch and options. All six accepted modes use root sources.
GED-001's metadata/frame and retained material-refresh target repairs remain in place;
CONFIG-007 now qualifies the live publisher's composition separately.

CONFIG-001--006 and earlier realization/prototype/policy repairs remain in place.
SPARSE-001 covers seven entry-local setters; PICK-001 registers retained-assembly ray
traversal, exposed by their renderer regression. SPARSE-002 adds retained intent and removes
three live reapply writers. Dated qualification belongs to
[readiness](obol_production_readiness.md#september-9-complete-shape-membership).

**Evidence recovery:** the cited historical `/tmp/obol*` and `/tmp/qged*` raw
reports, scripts, captures and caches are unavailable as of September 8. Their
dated summaries and formal baseline catalog remain, but cannot substitute for
raw current qualification. Reconstruct the open Lucy/planning/frame reproductions
and rerun the complete formal suite before release acceptance. New evidence is
retained under `.build/obol-qualification`, with hashes and matching runtime copies.

1. **S1--S2: inventory ownership, then prove one production path.** Complete
   the mutable-state and policy-writer audit below and map the remaining work
   to acceptance gates.  The source observation/consumption boundary from
   `CXX-SOURCE-001` now has focused closure evidence in the conformance record.
   Static-quality completed-frame selection is now in the existing trial
   owner, with the `CXX-STATIC-001` sentinel regression and
   [writer inventory](libbobol_tla_conformance.md#static-quality-writer-inventory).
   The [CXX-CURSOR-001 retarget regression](libbobol_tla_conformance.md#cxx-cursor-001--an-untouched-pass-acquires-unnecessary-rescan-debt)
   now distinguishes untouched passes from consumed prefixes.  The
   [submission writer inventory](libbobol_tla_conformance.md#submission-writer-and-successor-inventory)
   classifies all 51 direct mutations and the indirect capacity/repair entries.
   `CXX-DEMAND-001` fixes presentation repair consuming ordinary demand.
   [CXX-DEMAND-002](libbobol_tla_conformance.md#cxx-demand-002--source-reconciliation-discards-pending-demand)
   preserves ordinary demand through source coverage and normals repair while
   retaining minimum-coverage priority.  The
   [visible-demand follow-up](obol_production_readiness.md#visible-demand-replacement-follow-up)
   proves current projection, actual presentation and finite retirement.
   The [constrained-source follow-up](obol_production_readiness.md#september-5-shared-minimum-coverage-under-source-replacement)
   now adds visible entries beyond the consumed census prefix, constrained
   redistribution and the exact frame preceding its final allocation.
   `CXX-COVERAGE-001` fixes structural repair rejecting a rich shared prefix
   before trying its affordable minimum.  `CXX-TIMING-001` supplies the
   completed-duration boundary used for controlled frame evidence.
   The static lifecycle audit now enumerates all 19 direct trial mutations;
   `CXX-STATIC-002/003` close duplicate-configuration work and missing deadline
   policy/frame invalidation.  Start and interruption policy remains an S3
   extraction boundary.  The [detached-source inventory](libbobol_tla_conformance.md#detached-source-submission-and-release-inventory)
   now classifies the production job/queue/reservation boundary and closes
   `CXX-SOURCE-002`: failed submission preserves caller ownership and the
   pre-existing queue.  `CXX-SOURCE-003/004` add partial-pool cleanup, GED startup
   error containment and cancellation through the existing worker slots.
   The [compact-adoption audit](libbobol_tla_conformance.md#compact-source-adoption-boundary)
   closes `CXX-SOURCE-005`, which inferred payload delivery from matching
   metadata. The [delivery audit](libbobol_tla_conformance.md#cxx-source-006--producer-completion-is-not-successful-delivery)
   adds memory-denied delivery, failure retention across view policy changes,
   and retirement of an obsolete presentation binding. Its checkpoint records
   passed C++/formal checks and remaining graphical qualification.
   The [aggregation follow-up](obol_production_readiness.md#september-5-exact-zero-work-assembly-aggregation)
   fixes exact zero-work contributions stranding a denied source outside the
   view. The [publication-boundary follow-up](obol_production_readiness.md#september-6-policy-and-erase-publication-boundaries)
   closes the small replay's unnamed policy/erase transitions through the
   existing notification path. The [overview-extent repair](obol_production_readiness.md#september-6-retained-overview-extent-and-delivery)
   now preserves the geometry/transform pair, publishes generic exact coverage
   through the priority lane, and advances placement when an overview changes
   coordinate systems. Actual retained vertices and the two-backend small
   reserve/merge-denial replays qualify that boundary. Coverage before an exact
   bound is available and multiple-root failure remain unqualified. Internal
   merge allocation consistency now has the focused coverage linked below. The [source-routing repair](obol_production_readiness.md#september-6-source-routing-and-scoped-cancellation)
   now rejects queued output for both a replacement owner and an in-place path
   change with reused numeric stamps. Existing stream cancellation retires stale
   delivery, and worker tests preserve a healthy sibling. The
   [source-identity repair](obol_production_readiness.md#september-7-source-realization-identity)
   centralizes admission in the source and rejects same-owner population,
   database, representation-key and tolerance changes, including database
   rebinding A to B to A. The [stream cancellation repair](obol_production_readiness.md#september-7-cancelled-stream-ownership)
   now releases unconsumed stream geometry, staged imports and persistence
   journals while a healthy sibling remains active, rejects late publication,
   and preserves completed staging transferred to a consumer. Complete the
   remaining traversal/custom cleanup, service-lifecycle and journal-retention
   boundaries. CXX-CONFIG-008 closes the twelve direct source-field/engine paths.
   The [terminal-item repair](obol_production_readiness.md#september-7-terminal-item-resource-retirement)
   now releases unsuccessful sources/database handles and each finished item's
   callback context before terminal item publication, including admission denial
   and queued shutdown. GED borrows only completed results and no longer keeps a
   raw worker-source alias. Qualify the magnitude/lifetime of necessary completed
   results separately from worker reservations and real snapshot-file cleanup.
   The [queued-cancellation repair](obol_production_readiness.md#september-7-queued-cancellation-and-admission-wakeup)
   releases unstarted items without worker or memory admission, preserves healthy
   siblings and wakes newly admissible work before cleanup finishes. Its private
   weak endpoint and shutdown retirement counter have controlled lifecycle tests.
   Qualify queue/data magnitude and cancellation latency under real pressure.
   Cancellation cleanup is synchronous and scales with backlog (about 86 ms
   for the measured 150k queued metadata records); broader geometry, memory
   pressure and native cancellation limits remain unqualified. The
   [compact-publication repair](obol_production_readiness.md#september-7-compact-merge-publication)
   prepares one occurrence before commit, publishes a coherent prefix on failure,
   releases unpublished parts, and falls back to authoritative refresh when sparse
   journals cannot allocate. Its allocator sweeps and GED fault after the first
   committed leaf cover the prior untested internal-merge boundary. The
   [owned-index installation repair](obol_production_readiness.md#september-7-owned-source-index-installation)
   now prepares the complete candidate, including indexed presentation state,
   before transferring ownership and handle identity. Its allocator regression
   qualifies old-or-complete registry consistency and candidate release. The
   [terminal-adoption repair](obol_production_readiness.md#september-7-terminal-adoption-and-notification-recovery)
   now extends preparation/commit to streamed overview retirement, final metadata
   and staged ownership, with ordinary child/path retirement before field/node
   observers. It also repairs dependency notification cleanup after allocation or
   callback failure and contains denied terminal adoption in GED. The
   [child/path retirement repair](obol_production_readiness.md#september-7-child-path-and-auditor-retirement)
   now prepares removal before adoption, commits child/path state before direct
   path observers, and retires default graph/sensor references without allocating
   an auditor snapshot or dispatching unrelated callbacks. The
   [realization-state repair](obol_production_readiness.md#september-7-realization-state-publication-and-auditor-allocation)
   now prepares the status setter's source/owned-shape records together, preserves
   nested owners and removes the unrelated compact-display rebuild. It also
   closes allocation failures in the dependency's field-auditor lock and
   pointer-list capacity publication. The
   [direct-writer inventory](libbobol_tla_conformance.md#direct-realization-writer-inventory)
   now distinguishes private worker construction from live cached realization,
   external publication, compact edits, auxiliary/prototype publication and
   invalidation. The [qualified append repair](obol_production_readiness.md#september-7-child-append-ownership)
   closes `CXX-SOURCE-017` reference/parent-link ownership failures. The
   [external publication repair](obol_production_readiness.md#september-7-external-geometry-publication)
   closes `CXX-SOURCE-018` for line/point/mesh/annotation, clear/empty geometry
   and primitive adapters, with prepared child replacement and coherent
   compact/compiled retirement. It also repairs GED's noncompact occurrence
   lookup exposed by that retirement. The
   [live cached realization repair](obol_production_readiness.md#september-7-live-cached-realization)
   closes `CXX-SOURCE-019` for public/action cached publication, including
   existing compiled drawings, sparse overrides and provider failure. Workers
   keep explicit private construction; the obsolete preservation flag is removed.
   It also repairs enum mapping replacement during failed field copies.
   The [combination edit repair](obol_production_readiness.md#september-7-combination-edit-publication) closes
   `CXX-SOURCE-020`, including obsolete exact bounds and LoD asset cache lock
   recovery. The [leaf edit repair](obol_production_readiness.md#september-7-incremental-leaf-edit-publication)
   closes `CXX-SOURCE-021`: affected occurrences and owner metadata publish
   together, indexed paths resolve, exact bounds grow/shrink, and obsolete
   source profiles are invalidated. The
   [auxiliary line repair](obol_production_readiness.md#september-7-auxiliary-line-publication-and-retirement)
   closes `CXX-SOURCE-022` for named line replacement and prepared auxiliary
   removal/clearing; it also repairs doubled nested placement and ordinary
   controller structural revisions. The [auxiliary source repair](obol_production_readiness.md#september-8-auxiliary-source-publication)
   closes `CXX-SOURCE-023` for enclosing configuration/terminal publication,
   old-worker rejection and ordinary controller revision recovery after observer
   exceptions. It also repairs allocation-failed Obol field-sensor attachment.
   `CXX-SOURCE-024` now closes prototype source publication through the shared
   prepared boundary. `CXX-REALIZE-001` now closes source/action/scene progress
   publication and unwind retention. `CXX-REALIZE-002` now closes view successor
   publication, ordinary request preparation and trace-unwind failure handling.
   `CXX-SOURCE-025` now closes the inner invalidation record and producer-data
   retirement. `CXX-CONFIG-001` now prepares batched configuration and all three
   source setters while retaining sensors and complete owned metadata. The enclosing
   scene draw/representation/rename effects now close through `CXX-CONFIG-002`.
   Display-name, hierarchy and material-policy effects now close through `CXX-CONFIG-003`.
   Display and placement publication now close through `CXX-CONFIG-004`.
   Remaining scene state setters now close through `CXX-CONFIG-005`.
   Explicit database metadata and source/scene material refresh now close through
   [CXX-CONFIG-006](obol_production_readiness.md#september-9-database-metadata-and-material-refresh),
   including full-path records, compiler effects and bulk committed-prefix revisions.
   SPARSE-001 closes the seven entry-local display/metadata setters and retains sparse
   appearance and journal effects. SPARSE-002 closes retained rule/frontier/selection
   publication with the same owner. REGION-001 closes evaluated-region source and bulk
   scene publication. GED-001 closes draw-metadata scene effects and material-refresh
   target acceptance. Resume enclosing direct/deferred draw transactions,
   then the remaining classified success/failure writers,
   callers of legacy individual removal/truncation, traversal and custom-node cleanup,
   connected-field/engine and delayed-sensor
   exceptions, sparse setters and journal retention. The journal fallback cannot repair an inconsistent
   upstream writer. Continue the LoD service,
   GED, host and presentation writer/release-edge audit, acceptance-threshold closure and applicable
   graphical qualification before declaring S1/S2 complete. The SOURCE-018
   single/quad edit replays pass, but their small-pane viewport text overlaps
   and clips; retain those captures for S5 visual/layout acceptance. The
   cost-controlled test is numeric/effect evidence, not performance or full
   resource/platform qualification.
2. **S3: finish responsibility extraction.** Close the inventoried boundaries
   in the control-state and physical-extraction sections below.  Preserve
   compact data, immutable sharing, byte bounds and existing formal
   compositions.  The untracked `brlobol_rework.txt` is a proposal, not an
   adopted replacement specification.
3. **S4: close asynchronous/resource failures.** Retain the distinct quiet
   planning-cycle and desktop lost-frame reproductions.  Qualify interrupted
   preparation/cancellation, cold giant and multiple-root admission, shared
   cache contention, typed bounded BREP production and teardown.  At the
   source coordinator, the focused startup/cooperative-shutdown boundaries now
   pass (`CXX-SOURCE-003/004`); qualify actual geometry cancellation latency and
   resource magnitude under pressure, plus required native and shared-stack
   sanitizer runs.  The new Linux thread and controlled-reservation tests do
   not close those wider conditions.
4. **S5: finish rendering, scale and feature qualification.** The
   [handoff](obol_session_handoff.md) identifies current binaries and their
   checks. The [framebuffer feature checkpoint](obol_production_readiness.md#september-7-framebuffer-feature-evidence)
   resolves the nonterminal HUD discrepancy: an interrupted traversal retained
   older pixels despite an empty render latch. Both hosts now carry the traversed
   feature revision with their pixels. Eight current Generic rows, all strict
   traces and 38 matching HUD images pass; one older frame is correctly excluded.
   The checker still rejects a hidden current fill or an unpresented terminal
   error. The [passive-checkpoint entry](obol_production_readiness.md#september-7-passive-checkpoints-and-hud-frame-correspondence)
   retains the earlier Generic wire timing and software Lucy quality failures;
   the current pass does not close a different adaptive timing history.
   The [terminal-error publication repair](obol_production_readiness.md#september-7-terminal-error-hud-publication)
   now publishes the final red bar and retires its frame on both backends;
   its causal regression disables periodic samples. The
   [pose observation repair](obol_production_readiness.md#september-7-pose-deadline-observation)
   identifies a new hard-deadline interruption omitted by the earlier System
   warm-wire validator, and uses exact presented populations for retention.
   Preserve original failed reports and separate corrected verdicts. This does
   not qualify broader pose/zoom latency or retire other timing histories.
   The last warm software Lucy run fails the close floor at cut 28 despite
   equal cut-24 normals populations. Establish repeatable zoom/floor evidence.
   Preserve the earlier normals mismatch, planning-cycle trace, System GL
   retained-cut and default close-floor failures, and software cut-25-start
   failure; isolated warm passes do not retire them.  Extend same-cut
   normals/images, spatial continuity, cold/warm turnover,
   repeated/distinct assets, 50k/150k and real geometry evidence.  Keep the
   interaction/editing and shared-client requirements below.
5. **S6: qualify the final candidate.** Required native hosts, shared-stack
   sanitizers and release rows must qualify the final binaries.  Renderer
   defaults, visual floors, deadlines and memory bounds remain unchanged.
   Make dependency staging reproducible: the overview checkpoint reproduced a
   Qt crash after CMake replaced the qualified OSMesa library from an updated
   external bundle. Restoring the matching library repairs this checkout;
   qualification must pin a coherent dependency set and verify loader/hash
   provenance after configuration, not rely on that manual restoration.

The [readiness matrix](obol_production_readiness.md) owns detailed run evidence.
The [frozen resume history](obol_20260905_resume_history.md) preserves the
previous chronological priorities and reproductions.  Passing historical rows
do not qualify later binaries or different adaptive conditions.

## Maturity estimates — September 10

The [dated readiness estimates](obol_production_readiness.md#september-12-stage-completion-estimates)
are the sole percentage baseline. They are engineering judgments against the guide's exit
criteria, not effort or time forecasts. All six gates remain incomplete. S1–S3 define the
structural simplification finish line; S6 defines the production finish line. Update the
baseline when a gate's evidence changes materially, not after each individual repair.

## Gate assignment

This assigns the existing backlog to the guide's finish line.  It does not add
requirements or mark a gate complete.  Exact tests and limits remain owned by
their linked acceptance documents.

| Existing work | Gate and closure evidence |
|---|---|
| Remaining mutable-state/writer inventory and acceptance gaps | S1: [conformance inventory](libbobol_tla_conformance.md#s1-ownership-inventory-september-5-first-pass), complete writer/lifecycle classification and declared row thresholds |
| First complete production path | S2: source evidence now has deterministic decision/effect tests and removed caller policy; applicable graphical failures and the remaining inventory prevent declaring the gate complete |
| Remaining controller, source/service, GED and presentation responsibility extraction | S3: close the inventoried ownership violations with independent compilation and regression evidence |
| Planning/frame-delivery stalls, cancellation, byte bounds, cache failure, compaction and teardown | S4: finite progress/retirement, released reservations, applicable models and measured lifecycle/resource tests |
| Large-BREP provider | S4: bounded validated typed production; S5: original BREP image/growth/reclamation qualification |
| Small/giant/shared/distinct geometry, discovery/storage costs, 50k/150k, vehicle importance, normals and view turnover | S5: existing real-model and pressure rows, declared latency/memory limits and inspected images |
| Interaction/editing and shared-client feature completion | S5: existing editing and GUI matrix, including actual input and command readback |
| Required native hosts and full shared-stack sanitizers | S4/S5 for implementation/lifecycle closure; S6 for final-candidate evidence |
| Release qualification | S6: required rows qualify the exact final binaries with no unresolved production requirement |
| Cross-process generation lease | Conditional S4: required if measured cold contention violates the working-set contract; record the contention result either way |
| First-cold manifest OBB enrichment | Conditional S5: decide from coverage/importance qualification; no enrichment mechanism is required when the existing cue meets the contract |
| Terminal aggregate proxy for multi-page assets | Conditional S4/S5: required if measured pressure needs that representation; then qualify one atomic occurrence owner for all pages |

S1 remains incomplete until the partial writer audit and any unspecified
required acceptance thresholds are closed.  Existing passing evidence is
retained; it does not substitute for the untested compositions named above.

## Current implementation baseline

The following architecture is implemented and must not be reintroduced as
debt under another name:

- libged owns renderer-neutral draw roots, exact/subtree erase, selection,
  edit scopes, annotations, polygons, faceplate state, and typed scene deltas.
- libBObol owns database realization, compact occurrence inventories, cache
  and worker services, view demand, resource policy, presentation planning,
  and convergence evidence.
- Obol owns retained part/instance storage, renderer planning/execution,
  view-local preparation, picking primitives, and completed-frame evidence.
- Each view owns a distinct `SoCADAssembly` and `CadViewState`.  Immutable,
  admitted `PartGeometry` may be shared; camera-local plans and reports may
  not.
- Part and occurrence identities are distinct fixed-width strong types.
  Persistent mesh topology uses fixed-width values and is validated before it
  becomes renderer-visible.
- `PartGeometryBuilder` is the sole mutable producer form.  Its admission
  validates once and constructs a private, const `PartGeometry` snapshot;
  cached snapshots recover tokens in O(1), and renderers cannot receive a
  mutable alias.  Scene mutation batches have pure preflight functions and are
  atomic with respect to validation failure.  libBObol stages geometry,
  instances, and semantics without touching the live assembly, then commits
  through checked mutation calls.
- Complete-scene replacement constructs its candidate before the no-throw
  publication swap and reports resource denial without changing the live
  scene.  Sparse publication remains proportional to its bounded journal:
  validation and libBObol staging are atomic, but Obol does not clone a 150k
  retained scene to recover from process-wide allocator exhaustion during the
  mechanical commit.  Capacity is reserved from the manifest and sparse
  journal size remains bounded.
- Deterministic precommit denial covers retained-scene, presentation, durable
  hierarchy-record, and final name-mapping publication.  A denied scene commit
  preserves the prior node notification identity; a denied cache replacement
  preserves the prior discoverable immutable hierarchy.  Filesystem crash
  consistency remains a separate operational gate.
- Obol update windows are move-only, nest-safe RAII values.  Manual
  `beginUpdate()`/`endUpdate()` pairing is no longer public.
- One retained assembly presentation never chooses CAD LoD.  It executes the
  producer-certified cut and reports exact cost/resource evidence.
- Progressive mesh policy reads one shared immutable generation snapshot per
  allocation candidate; enrichment or compaction publishes a new generation
  without invalidating an in-flight transaction.  Renderer timing evidence is
  point-threshold stamped, and hosts must classify capacity relevance and
  actual CAD execution with distinct required types.
- Progressive spatial pages are storage/render partitions of one logical CAD
  occurrence.  They inherit its transform, style, selection, picking, and
  semantic path and are not promoted to independent CAD objects.
- Optional PCA-oriented bounds are immutable cache/part metadata, not a LoD
  state.  Monolithic retained parts use them for screen-significant batched
  aggregate boxes in shaded and wire modes; AABB remains the validated
  fallback.  Exact Obol admission checks real renderer positions rather than
  requiring a rotated OBB to contain the artificial corners of its own AABB;
  only an invalid optional OBB is discarded.  The optional projection pass
  occurs after cold coverage can be published, so it cannot delay the first
  useful preview.
- Workload names and renderer names are qualification profiles, not control
  modes.  Lucy, Hubble, 50k, 150k, System GL, and OSMesa all use the same
  state transition relation.
- Authorization evidence is exact.  Renderer cache keys retain complete typed
  tuples, while capacity populations and bounded structural frontiers receive
  non-reused tokens only after exact comparison.  Compact-source and
  cross-source presentation reuse compares exact revision/membership vectors.
  Digests remain diagnostic, content-addressed, or bucket-selection values and
  never independently certify current mutable state.  Multi-assembly draw,
  timing, and resource observations use exact canonical source tuples keyed
  by non-reused Obol assembly identities; aggregate tokens never depend on an
  object address, hash, or saturating sum.

Focused Obol CAD tests, libBObol rendering/realization/update tests,
`ged_test_obol_draw_sync`, the pivot guard, and the installed-package consumer
pass after the 2026-08-28 control-contract cleanup.  This is a development
gate, not production clearance.

## P0: finish the control-state reduction

The canonical contract is `libbobol_progressive_pipeline_contract.md`: one
six-domain evidence stamp (inventory, availability, visibility, view, policy,
capacity),
at most one bounded plan cursor, at most one presentation transaction, and a
finite work ledger.  HUD outcome is a projection of those values, never a
second phase machine.

Much of this contract is implemented: typed revision advancement, stale-plan
rejection, finite interaction/static-quality/point-quality states, one
presentation transaction, bounded preparation evidence, and the control
refinement map all have focused tests.  The former 6.1k-line coordinator
header is now a small umbrella over independently compile-checked admission,
capacity, delivery, scene-evidence, presentation, and view-policy boundaries.
Submission source repositioning now atomically resets its source-local plan,
allocation cursor, and visibility census.  Fresh/retired submission passes
atomically consume predecessor rescan state, while deliberate inventory pauses
retain it.  Deadline-safe cuts, retained visual-error bounds, and renderer
timing/upload observations are keyed allocation-free evidence values rather
than independently writable companion scalars.
Terminal renderer-capacity evidence now remains authoritative through both
the final mechanical cut application and subsequent ordinary admission
pumps.  The completed allocation is one typed policy transition: it either
retires only in-flight measurement while preserving a current certificate, or
atomically invalidates the certificate and its budget-limited witness.  A
current terminal certificate caps ordinary planning at its safe budget, so a
throughput estimate cannot immediately reselect the rejected discrete
population and recreate the same search.
Constrained terminal outcomes now expose a typed evidence mask and the runtime
refinement checker rejects an unwitnessed constraint.  This distinguishes an
honest deadline/memory/presentation endpoint from a controller which merely
stops while visual debt remains.  Keep this evidence derived from the owning
policy values; do not add another writable outcome-reason field.
`ObolLodComposition.tla` now closes the formal seam between admission,
retained growth, exact presentation, capacity, structural repair, point
quality, static quality, and terminal publication.
`ObolTerminalConvergenceComposition.tla` closes the previously atomic tail of
that seam.  It permits the production relationships which matter near
termination—an exact visibility census followed by its exact framebuffer
classification before reallocation, a distinct static population reopening a
bounded capacity search,
each exact presentation consuming monotone renderer-preparation units, an
over-budget allocation consuming a strictly coarser local representation, and
quiet compaction after foreground convergence—but requires every re-entry to
consume a finite semantic member.  Its current 2026-08-31 TLC run checked 639
  distinct states to depth 81.  The corresponding executable handoff now
removes the temporary global ceiling after a successful occurrence-local
allocation; only a revision-bound protected-minimum constraint may retain it.
This eliminates the former path which re-walked global PoP ordinals after the
local allocation had already succeeded.
`ObolControlLifecycleComposition.tla` closes the orthogonal policy/provider/
host seam.  It proves that provider registration is not provider work,
policy-off camera bookkeeping cannot arm automatic demand, exact-frame debt
distinguishes work awaiting a request from work already attached to a frame,
the first capacity-relevant hard-deadline miss owns quality recovery and a
fresh exact-frame obligation, and terminal/HUD outcome is derived rather than
stored.  The 2026-09-01 TLC pass checked 1,092,377 generated / 227,787 distinct
lifecycle states to depth 32.  It includes two bounded interruptions which
retire a runnable demand cursor without discharging its level-triggered
importance-census obligation, and semantic-only selection/style mutations
before and during exact presentation.  Those mutations advance their own
presentation revision but preserve every LoD control fact; even an interrupted
semantic-only frame may request only an exact repaint, never capacity recovery.
The then-current
28-fact refinement map checked 57,344 states: every concrete field
independently and every combination of its ten distinct owner classes.  The
offline checker now requires that complete fact mask, independently projects
its obligations and fixed-precedence owner, and includes it in A/B/A cycle
identity.  Its adversarial gate rejects alias-only cycles, owner/obligation
mismatches, a presentation owner without a typed finite successor, and a
missing concrete mask.  Dense 50k OSMesa tracing then found two refinement
defects which sparse event samples had missed: a consumed capacity candidate
could publish a replacement plan without advancing the capacity domain, and
the runtime validator did not recognize the controller-scoped PUMP-to-RENDER
transition already present in both `ObolControlLifecycleComposition` and
`ObolHostWork`.  Capacity replacement now advances its semantic revision, and
presentation diagnostics export the exact witness source; an unrelated shared
host pump cannot satisfy that check.  The resulting dense trace passes the
refinement checker, although the 30-second 50k OSMesa gate can still expire
while monotone command preparation and worker realization are active.  This
was distinguished from an end-state cycle by an extended full workflow: the
warm 50k OSMesa selection/erase/redraw replay completed in 31.6 seconds, and
its three post-selection convergence waits returned to ownerless readiness in
0.55--0.92 seconds with no obligation or violation.
The controller transition journal now closes the remaining observation gap.
It records typed before/after endpoints at every registered effect boundary,
includes renderer-preparation target and remaining-unit ranks, and emits an
explicit `unnamed` failure record if observable state changes outside those
boundaries.  The qged checker rejects missing fields, discontinuity, event/
owner mismatch, truncation, and every prior sampled-state violation.  Its first
complete trace found a queued-render-to-in-flight-frame seam: exact-frame debt
temporarily had no witness after host claim.  Production now exposes the
`renderInFlight` state already required by `ObolHostWork` as a typed claimed-
frame witness and retires it on frame completion or interruption.  A 512-step
deterministic C++ trace and a 57-transition OSMesa draw/zoom trace pass with no
unnamed transition or drop.  Keep this journal diagnostic-only and bounded;
do not turn it into another scheduler or a monolithic policy facade.
One September 10 `libBObol_lod_update_action` invocation run concurrently with
two other test executables reported a missing current exact-frame witness. The
saved predecessor and four isolated current invocations pass. Reproduce this
under controlled contention as S4 load-sensitive work before deciding whether
the test timing contract or production witness ordering needs repair.
The 2026-08-30 implementation refinement now gives exact-frame debt distinct
request-required and frame-awaited states, including reattachment to a
coalesced host request after a newer semantic mutation supersedes its target.
Provider registration and provider pending remain separate at the host/status
boundary.  One derived controller-local pump projection covers all reducer
obligations, while service queues and provider work are composed at the host
boundary.  The policy-off transaction asserts the complete modeled retirement
postcondition and also retires the inventory-coalescing deadline.  The
controller-to-ledger projection is now canonical as well: inventory
coalescing, visibility deferral, source deltas, quality probes, resumable
retained allocation, importance census, and resident-admission retry can no
longer keep the pump active while disappearing from convergence diagnostics.
At that stage the public diagnostic snapshot also began carrying the complete
concrete fact mask; the current 29-fact inventory is defined by
`ObolControlRefinement` and the exhaustive C++ map, and its proof boundary is
cataloged in `tla/models.json`.
This distinguishes those aliases without exposing private reducer objects.
Focused off/on/off, idle-provider, exact-frame, and host-wakeup tests exercise these
bridges.  Requesting an exact or batched-publication frame now transfers the
runnable host level from PUMP to RENDER instead of polling an action which
cannot advance before that frame completes; `BObolProgressiveStatus::hasMore`
still reports the unfinished transaction without becoming a second scheduler.
The HUD audit found no independently writable terminal outcome or controller
phase: its phase and completion are projections of current evidence, while Qt
retains only publication-change observations.
The 2026-08-29 ownership audit retired the last detached lifecycle combinations
relevant to proxy admission: structural relaxation is one four-state value,
the recovery-plan witness lives with point-quality state, renderer-feedback
consumption lives with the interaction session, and the host render/capacity
request is one three-way value.
`ObolAssetPublicationComposition.tla` and `ObolCadFrameComposition.tla` also
close the adjacent producer/cache and retained-scene/frame seams.  The first
composition run found a real ownerless successor: a producer constrained under
an old view demand did not resume for a newer demand.  Production now retries
the current demand, and the focused service test proves a new-view result is
actually published rather than merely clearing the old failure query.  The CAD
frame seam rejects stale target reports across concurrent scene/view changes.
The production shape is nevertheless still concentrated in two places:

- `lod_admission_policy_private.h` now retains the allocation-free evidence
  values and planner declarations; the 71-method numeric planner is compiled
  independently in `lod_admission_policy.cpp`.  Point and structural evidence
  remain in the same 1.6k-line header because their small methods preserve the
  trivially-copyable value contract, but any further extraction must keep the
  header-alone compile and exhaustive policy tests.
- `view_lod_coordinator_state_private.h` is about 1.6k lines.  Submission,
  repair, capacity, presentation, point-quality, and interaction ownership are
  cohesive typed values.  Further changes should remove a measured effect-
  shell responsibility, not add accessors or a generic lifecycle wrapper.
- `view_controller.cpp` is about 8.3k lines after lifetime/camera/host-request,
  exact-picking, and residency extraction.  It remains both reducer caller and
  the LoD policy effect executor.  Exact-view quality history is now an
  independently compiled allocation-free policy rather than another inline
  controller implementation.

Required completion work (`tla/RISK_COVERAGE.md` supplies the detailed
formal/executable evidence mapping for the audit-derived items, and
`libbobol_tla_conformance.md` records the model-to-production conformance
findings).  Result authentication and shared-producer lifetime now have
production gates, focused formulas, and complete asynchronous lifecycle
matrices; they are no longer debt:

1. Finish inventorying the remaining mutable controller fields by sole owner,
   revision domain, progress witness, and terminal transition.  Delete
   write-only or derivable fields; do not split a keyed certificate back into
   convenient writable scalars.
   The [September 5 inventory](libbobol_tla_conformance.md#s1-ownership-inventory-september-5-first-pass)
   classifies the first state families.  Preserve the now-closed
   `CXX-SOURCE-001` source-order/lifetime regression; the remaining host,
   source/service, GED and submission/presentation audit still needs closure.
2. Continue moving nontrivial pure-policy method bodies behind the new
   compile-checked private boundaries.  Keep hot occurrence storage dense and
   allocation-free; this is a responsibility split, not an object hierarchy.
3. Preserve the complete transition-journal gate whenever a controller effect
   boundary changes.  An `unnamed` event, dropped record, discontinuous
   endpoint, or owner without a finite witness is a contract failure, not a
   diagnostic warning.
4. Re-run the applicable TLC models after ownership/liveness changes and the
   graphical matrix after numeric policy changes.  Formal models do not judge
   visual quality or wall-clock performance.
5. Preserve the completed identity gates: evidence stamps require all six
   typed domains, inactive holders use the administrative sentinel, and every
   authentication, ordering, or ownership counter uses the checked fail-stop
   successor.  Diagnostic counters may saturate but must never authorize work.

Acceptance: unchanged evidence cannot reopen planning; an invalid/stale plan
cannot commit; one event cannot select two successor owners; no terminal HUD
state has foreground work; no nonterminal state is ownerless; asynchronous
results authenticate their complete identity; shared cancellation has no live
consumer lease; and constraint evidence is never relabeled as failure.

Two 2026-08-28 150k System-GL traces exposed the same missing refinement at
different effect boundaries.  A point-calibration request pauses the active
submission cursor, so that cursor cannot also be counted as the frame producer
which will present the calibration.  Producer classification is now one typed
input value, all 256 input combinations have an executable truth-table test,
and the controller recomputes ownership after deadline effects mutate the
calibration state.  Runtime refinement validation additionally rejects a
`PRESENTATION` owner without a requested frame, independent producer, or
finite publication timer.  Retain those checks while completing the reducer;
an owner label, progressive-pump level, or paused cursor is not a progress
witness.

After immutable generation snapshots removed repeated shared-generation loads
and CAD timing became threshold stamped, the 2026-08-28 warm 150k shaded
replays terminated box-free in about 40.4 seconds on System GL and 36.8 seconds
on OSMesa.  The prior 32/64-pixel balancing cycle was gone: terminal samples
have no owner, obligation, pending foreground/background work, or control
violation.  This resolves the known same-population liveness defect, but not
all scale debt.  Profile true-cold useful-preview latency and cache I/O using
common policy events rather than an object-count regime.

The stricter 2026-08-31 matched-quality workflow selects materially richer
populations.  Its initial 150k System-GL view reaches 20.49M faces through
78,341 meshed occurrences in 220 seconds; OSMesa reaches its constrained
1.17M-face endpoint in 51 seconds.  Both are ownerless, box-free, proxy-free,
and have no prominent-floor violation.  A current 50k `perf` replay confirms
that the former 13.6% Obol scene-mutation snapshot hotspot has been retired;
no current exclusive symbol exceeds 3%.  Treat the remaining 150k latency as
distributed realization/publication throughput debt.  Instrument phase-level
work and test batching or parallel publication before proposing another local
micro-optimization.  `obol_lod_visual_quality.md` owns the exact metrics and
the distinction between matched controls and unsafe whole-scene controls.

The 2026-08-29 draw-metadata migration replay exposed a separate cross-owner
gap: a policy revision could create a coverage census after the exact
camera/source proof was already current, retire its cursor through a stronger
capacity path, and leave compaction waiting on the orphaned census.  Coverage
is now an explicit inventory obligation, its producer is level-triggered, and
policy-only quality changes preserve the existing proof.  The corrected warm
System-GL 50k/150k matrices both pass; the 50k endpoint is fully ownerless,
while 150k presents a terminal box-free framebuffer and truthfully continues
bounded resident-prefix compaction in the background.  Do not weaken that
distinction by treating background reclamation as unfinished visual work.

A subsequent 50k exact-subpath erase exposed a different revision-boundary
defect: compact retained visibility changed without changing immutable mesh
inventory, so the presentation was correct while the convergence denominator
remained at 50,000.  Mesh inventory and effective presentation visibility now
have independent monotonic revisions and bounded journals.  The controller
unions their changed-entry sets, while the host performs a level-triggered
revision check after presentation synchronization.  Exact visibility is also
a first-class revision in retained-allocation plans.  It requires one
successor allocation after its bounded source delta is applied and the
resulting framebuffer has been classified exactly.  The source census alone
is insufficient: a restored occurrence may have lost its retained payload
while hidden, and only the successor frame exposes the structural replacement
work which must precede allocation.  This edge does not
invalidate the renderer-capacity certificate: visibility changes allocation,
not the cost model.  A newer exact edit supersedes an unpublished predecessor
allocation so the latest delta cannot be stranded behind stale work.

The final warm OSMesa 50k workflow updates 50,000 to 49,984 and back in
1.36/1.40 seconds while retaining 1,206,806 faces, zero boxes, and the same
5,130,680-unit certified budget.  The corresponding 150k workflow updates
150,000 to 149,984 and back in 6.05/3.96 seconds while retaining 835,318 faces,
zero boxes, and the same 1,446,651-unit budget within one rounding unit.  Both
terminate ownerless and the six-domain runtime trace passes.  Preserve this
distinction: visibility
is a planning input, not a geometry-inventory or renderer-capacity mutation,
and presentation-only edits still require a progress witness.

A later timing-sensitive 50k replay made this ordering failure deterministic:
erase-time allocation retired three visible mesh payloads, redraw updated the
50,000-occurrence census, and immediate allocation inspected the predecessor
frame.  Sixteen restored occurrences were absent from both the terminal ledger
and every producer.  The controller now requires an exact presentation frame
on every exact visibility delta, and
`BObolLodPlanningObligations::exactVisibilityReallocationReady` is the sole
gate from that frame to reallocation.  The strengthened submission and
terminal-composition models prove the ordering and quiescent liveness.  The
post-fix warm 50k OSMesa workflow at
`/tmp/qged-visibility-prerequisite-50k-osmesa-20260901` restores all 50,000
occurrences in 0.83 seconds, box-free and ownerless.

The pre-typed-host 150k OSMesa report had three quality-floor misses, including
one prominent synthetic occurrence.  The rebuilt host rerun has zero total or
prominent floor misses and zero visual-importance debt, confirming that the
old faceplate/CAD timing mix affected the selected capacity rather than proving
an allocator-ordering failure.  Realistic wheels/blades/hulls must still
demonstrate that the scene-wide importance ordering preserves genuinely
prominent forms.  Treat any future unaffordable residual as numeric
allocation/qualification debt; do not reopen a terminal control search merely
because proven-constrained quality demand remains.

## P0: complete physical responsibility extraction

The initial extraction and the Obol renderer split are real improvements, but
several units still obstruct review:

- `database_source.cpp` is about 18.3k lines.  Compact CAD presentation is now
  an independently compiled 2.1k-line owner, and its bounded copy-on-write
  staging journal remains side-effect-free until Obol accepts the complete
  renderer transaction.  BoT geometry conversion and normal construction are
  also independently compiled worker-safe units.  The next genuine boundary
  is discovery/traversal versus worker realization orchestration; do not move
  hot compact records into per-occurrence objects while extracting it.
- `draw_obol.cpp` is about 12.9k lines after endpoint, geometry, overlay, and
  scene-record extraction.  The libged private mutation vocabulary is now
  consistently named `ged_scene_reducer_*`; separate reducer orchestration
  from source realization/backend effects without reintroducing an adapter
  transaction type.
- Obol's former 5.8k-line `SoCADAssembly.cpp` is now a roughly 2.7k-line Coin
  node/mutation/action surface.  Retained plan maintenance is compiled in
  `CadAssemblyPlan.cpp`, camera-local subpixel classification in
  `CadAssemblyClassification.cpp`, and their state-only private contract in
  `CadAssemblyImpl.h`.  Renderer indirect, instanced, flattened, and
  retained/direct executors remain separate.  Further extraction should occur
  only at a measured retained-cache, picking, or Coin-action boundary.

The September 5 audit also measures `lod_service.cpp` at 7.5k lines,
`mesh_lod_cache.cpp` at 6.9k, `view_lod.cpp` at 6.3k, and
`mesh_lod_submit_action.cpp` at 5.9k.  Line count is a review-cost signal, not
proof of a faulty boundary.  In particular, separate service task/lease,
residency/reservation, and persistence responsibilities without losing the
shared producer lifetime or adding per-occurrence ownership overhead.

Do not replace large files with textual `.inc` fragments, mirrored state, or
thin forwarding abstractions.  Independent compilation must enforce one-way
dependencies.  Preserve sparse deltas, bounded 512-record publication,
immutable geometry ownership, and the no-second-150k-vector rule.

Remove superseded APIs and code in the same change which replaces them.
Compatibility with intermediate branch APIs/cache formats is not required.

## P0: graphical release qualification

Run the exact final binaries using true-cold isolated caches and certified
warm caches.  Required cases and controls are enumerated in
`obol_production_readiness.md`; the minimum model set is Generic Twin, Lucy,
multi-Lucy and xpushed multi-Lucy, Hubble, Havoc, NIST BREPs, heterogeneous
50k/150k scenes, Stanford meshes, and independent multi-gigabyte vehicle
models.

For both System GL and OSMesa, cover shaded, wire, hidden-line/evaluated modes
where applicable, LoD on/off, resize and fractional DPR, zoom, rotation,
translation, selection, exact/subpath erase-redraw, and memory turnover.
Compare System GL and OSMesa semantics/images within declared tolerances.
Use APNG and apitrace when diagnosing corruption, flashes, or camera jumps.
Use the matched full-detail methodology and provisional corpus targets in
`obol_lod_visual_quality.md`; libicv SSIM/PHASH and silhouette disagreement are
required evidence, not substitutes for named-feature inspection.

Release evidence must prove:

- a useful bounded initial presentation and truthful cold-work HUD;
- monotone useful convergence or an explicit typed resource constraint;
- no box/mesh cycling, mesh-to-box regression, holes from invalid PoP
  topology, lingering fallback boxes, or premature `View ready`;
- protected prominent shapes meet the visual-significance floor;
- rotation/translation restore proven affordable detail, while zoom starts
  from retained detail and changes demand incrementally;
- invisible/noncritical residency is reclaimed under pressure without OOM;
- input remains responsive and expensive rendering is interruptible; and
- repeated view turnover settles without leaked work or zombie GUI processes.

Run ASan/UBSan across the shared dynamic stack with worker teardown, corrupt
cache input, endpoint replacement, plugin reload, edit cancellation, and rapid
view close.  Run TSan/LSan where the native runtime supports them.  Static
linkage is reserved for explicit static-link tests.

## P0: interaction and editing completion

`qged_editing.md` owns the detailed editing plan.  The remaining production
gate explicitly includes:

- preserve the completed descriptor-classification gate and audit every
  advertised operation's runtime legality, rejection, constraint preservation,
  error reporting, and command readback;
- complete the reusable manipulator vocabulary (including ARB vertex/edge/
  face interaction) and richer primitive-specific manipulators;
- implement specialized sketch curve/profile interaction using the current
  libqtcad widget only where it remains useful; delete obsolete demos/widgets;
- validate qged and MGED/gsh sessions against the same edit command state,
  including final MGED `sed`/`oed` qualification;
- drive real mouse paths for polygon creation, resize, move, and boolean
  operations, not only equivalent GED commands;
- cover point/rectangle selection and modifiers, hierarchy-scale tree row
  styling, selection/highlight persistence across draw/erase/edit,
  lighting/navigation/axes/grid/framebuffer faceplate control, single/quad
  layouts, resize, and fractional DPR; and
- prove GUI widgets, command readback, retained manipulators, and database
  state agree regardless of which surface made the last change.

The deterministic measurement gate now covers 2D and exact-hit 3D mouse
gestures, distance/angle and degree/radian readback, cancellation, tool
replacement, resize, single/quad layouts, and System GL/OSMesa presentation.
It remains part of fractional-DPR hardware repetition, but is no longer an
unimplemented interaction mode.

The deterministic view-settings gate now covers bidirectional GED/widget
state for ADC, center dot, grid, model/view axes, scale, parameter/FPS text,
framebuffer composition, and the cutting plane in single/quad layouts on both
renderers.  It verifies retained Obol feature presence, command readback,
visible clipping, and exact restoration.  Fractional-DPR repetition and the
interactive in-scene cutting-plane affordance remain; ordinary settings
readback is no longer an unqualified control path.

The existing edit runtime and deterministic GUI replays are regression
substrate, not a substitute for physical-pointer, plugin-lifecycle,
hierarchy-scale, and renderer-backed qualification.

## P1: scale, cache, and geometry completion

- Measure discovery on cold local storage, warm page cache, and slower storage.
  Parallel read-only discovery must be bounded by I/O and memory and publish
  ordinary librt/GED hierarchy data without a parallel object ecosystem.
- Complete xpushed multi-Lucy cold preparation and visibility-turnover tests.
  Shared-instance reuse and thousands of genuinely distinct assets are
  separate workloads and both are required.
- Qualify cold spatial previews as globally representative, budget-aware
  presentations.  A local source page may never replace whole-object coverage.
- Preserve the cross-renderer NIST adaptive-wire regression and the passing
  shaded NIST growth/zoom-out/reclamation/cache-restore lifecycle.  Repeat the
  same constrained-residency qualification on at least one large real BREP;
  the partial Big Boy BoT is not evidence for its original BREP hierarchy.
  The current indexed-face guard rejects the known partial Big Boy tire and
  cache version 3 prevents its replay, while fresh-cache NIST remains green.
  This is only fail-closed containment.  Replace the legacy aggregate-success
  CDT call with one bounded provider contract carrying deadline,
  memory/result limits, cancellation, per-face completion, and typed outcome.
  Both LoD-on and LoD-off shaded drawing must consume that same validated mesh
  before the original Big Boy hierarchy can become a release gate.
- Verify constrained-memory cache reload, corrupt/incompatible cache
  invalidation, background compaction, and cancellation during persistence.
- Qualify simultaneous true-cold qged processes against the same large-asset
  cache.  LMDB transaction publication and the final completeness witness have
  passed the Generic Twin cross-process correctness gate, but independent
  processes still duplicate source classification until one publishes that
  witness.  Measure the memory/CPU cost on a large unique asset and add a
  crash-recoverable per-content generation lease if the stampede can violate
  the working-set contract.
- Measure frame, publication, cache, and resident-memory costs at 50k/150k and
  multi-gigabyte scale.  Preserve essential summaries and shared fixtures;
  remove disposable captures and duplicate large directories after each run.
- Extend the passing heterogeneous 10k OBB cold/warm matrix to no-LoD
  ground-truth comparisons on realistic wheels, blades, booms, and hulls and
  to the 50k/150k pressure range.  Persistence and cross-renderer selection are
  already qualified: the constrained OSMesa pair deterministically retains
  693 OBBs across cold-manifest and warm-manifest runs.  Record the perceptual
  improvement over AABB and the cache-generation cost on real geometry.
- Decide whether a first-cold discovery manifest should be enriched after
  detached PoP characterization.  The live payload already carries its OBB;
  the sealed discovery journal and its cache lifecycle lock deliberately make
  the initial structural cue an AABB.  Any enrichment must remain bounded and
  must not introduce nested draw-cache access or a second occurrence registry.
- If multi-page assets need terminal proxying under measured pressure, add one
  occurrence-level aggregate owner which suppresses/restores all pages
  atomically; never attach the whole-source OBB independently to each page.

## Exit condition

The new stack is production-ready only when the exact final binaries pass the
release matrix on required renderers/platforms, the remaining controller and
source boundaries are understandable and independently testable, real large
models meet the fidelity/responsiveness/memory requirements, and qged, MGED,
gsh, Archer, and rtwizard retain all user-facing drawing and editing behavior.
