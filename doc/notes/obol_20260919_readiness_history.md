> Historical snapshot preserved on 2026-09-19. Current instructions and status
> live in [obol_production_readiness.md](obol_production_readiness.md). Dated passes and percentage estimates
> below do not qualify the current checkout.

# Obol production-readiness matrix

Last reviewed: 2026-09-12

This is the release qualification plan for the BRL-CAD Obol drawing stack.  A
row is green only when it was run against the final binaries and records its
model, renderer, cache state, event script, viewport/DPR, report, images, and
performance evidence.  Dated evidence is retained below; the order of the
narratives does not establish which binaries are current.  Resolved failure
mechanisms belong in `libbobol_engineering_lessons.md`.

TLC counts retained in evidence narratives are dated snapshots.  The sole
current formal-suite result and count source is
`tla/baselines/tlc-2.19.json`.

The [simplification guide](obol_simplification_guide.md) defines migration and
completion gates; this matrix owns release acceptance and run evidence.
For current blockers and reproducible resume commands, start with the
[session handoff](obol_session_handoff.md) and its linked latest evidence.
The historical matrix is evidence of exercised behavior, not an all-green
declaration for the latest shared binaries.  In particular, Lucy quality is still unqualified
and the large-scene warm caches cited by older rows have been removed.

For acceptance requirements, go directly to [common criteria](#common-acceptance-criteria),
[required models](#required-models), the [release matrix](#matrix), and
[evidence retention](#evidence-retention).  The dated narratives preserve
reproductions; their historical pass labels do not close current backlog items.

The current stack publishes cached, external and prototype geometry through
prepared source boundaries. Realization commits action progress, scene revisions
and the view's standing render request before observers, preserving completed
prefixes and attempting the host wakeup despite observer exceptions. The
reproduced policy-enable failures remain repaired. Invalidation now prepares its
realization records and producer-data retirement before observers. Batched
configuration and the separate draw, representation and retarget setters now prepare
complete metadata and retain sensors. Scene draw/representation/rename now prepare
paired channels, indexes, conflict retirement and revisions before observers.
Display-name, hierarchy and material policy now publish complete metadata and scene
revisions before observers. Compact database colors survive policy changes, and hierarchy
updates retained renderer parents and compiler evidence. Placement now publishes fields,
child paths, compact matrices and scene effects before observers. Display-state and
display-patch publication now prepare sparse appearance and revision effects together.
Bounds, explicit realization, role-flags and view-policy setters now include scene
revisions and compiler effects in their commit. Database metadata and material refresh
now prepare full-path records and commit source/scene effects before observers. Seven
entry-local display/metadata setters now prepare changed records and retain sparse
appearance; native retained-assembly ray traversal now dispatches its picker.
Retained rules/frontiers/selection now publish through that same prepared owner. Evaluated
region edits now prepare all owned shape markers and publish each source's scene effect
before observers; completed prefixes and later callback edits survive. GED metadata now
includes scene effects; material refresh retains original targets and commits each changed
source's frame revision before observers. The obsolete leaf redraw branch is removed.
Scene source publication now composes configuration, display, policy, roles, view and
placement with source insertion/moves, indexes and revisions before observers (CONFIG-007).
Root redraw preserves current source appearance directly (GED-003). Group metadata now
commits with source/scene publication before observers, and standalone group intent/display
setters share prepared scalar publication (GED-004). Group creation now shares the source
publisher's hierarchy preparation and commits the entire path, indexes and revision effects
before observers (GROUP-001). Source movement now uses the same prepared source/scene
publisher, and both membership-only and composed publication reject cycles through newly
prepared descendant groups (GROUP-002). The group writer inventory bounds the remaining
removal and root replacement audit. Group rename now composes names, descendant
paths, indexes and revisions before observers, preserving source state and accepting
later callback edits (GROUP-003). The native name registry also preserves its previous
registration on allocation failure. Child append/removal now composes complete subtree
membership, indexes, repository ownership and revisions before observers (GROUP-004).
Shared nodes retain surviving edges and lookups; prepared cache updates preserve published
geometry identity. Recursive group removal and whole-group clearing now use that same
publisher, including complete path effects and surviving shared ownership (GROUP-005).
Source removal and bulk clearing now share complete child/scene publication, retain the
source selected by a colliding path lookup and preserve duplicate-edge path positions
(GROUP-006). The obsolete live source removal/index writers are removed.
Shape movement/removal now uses the same prepared child publisher, committing both
parent edits and path/scene effects before observers while preserving shape data and
shared geometry (GROUP-007). The four scalar shape-state setters now share preparation
of requested fields and frame effects, retaining original targets and accepting later
callback edits (GROUP-008). Repository clear, view invalidation and rename now use prepared
cache retirement with retained scene invocation targets. Root replacement, source indexes,
repository ownership, shared-controller propagation and realization child publication are
qualified together at the ROOT-001 boundary. GED local shape caller acceptance remains.
Other direct writers, Lucy quality and broader
resource/platform qualification remain open.
The formal baseline catalog is unchanged; its earlier raw `/tmp` output is unavailable.

## September 10: complete root and repository ownership

`CXX-ROOT-001` is repaired and qualified at its production boundary. The retained
GROUP-008 library rejects `root-publication-transaction`, `root-publication-edges` and
`repository-share-transaction`; the current library passes all three. Root replacement
now commits the root, complete indexes, repository membership and revisions before
retirement callbacks. Shared graph edits prepare every affected controller's lookup and
repository ownership. Repository clear, invalidation and rename retain their invocation
target and complete cache bookkeeping before releasing values. Last-controller teardown
does not allocate.

The integration run exposed two contract violations outside the original reproducers.
Ordinary realization replaced source children and advanced the structural revision even
when database-source membership did not change. It now requests a frame effect, while the
prepared membership diff promotes the transaction to structural when a nested database
source is actually added or removed. A view-controller regression also edited an attached
root with raw `SoGroup` calls, leaving the controller's authoritative source index pointed
at a retired node. It now uses the scene-controller removal/append transaction. The public
controller contract states that source and named-group membership changes after attachment
must use this boundary. Fixture construction may still finish the raw graph before the
controller is attached.

The complete source suite passes **45 direct checks plus 21 scenarios in 114.20 seconds**;
the state suite passes **66 checks in 179.20 seconds**. The identical 22-test matrix which
first reported two failures now passes **22/22 in 187.97 seconds**. An isolated full build
passes in **76.25 seconds**, and the five-client target set is current. Six translation
units compile independently, GED synchronization passes in **8.62 seconds**, the renderer
smoke test passes, and **23 freshly linked probes provide 29 passing rows**. Native Obol
passes **1,577 unit tests with one existing skip** and **63 integration tests**. OSMesa
remains `097aac27...`.

Fresh graphical evidence is retained under
`.build/obol-qualification/20260910-root-integration-2`. All **24 rows pass in 129.01
seconds** of summed stage execution: two eager rows, ten delivery/resource-denial rows,
eight Generic Twin cold/warm shaded/wire rows and four 207-event primitive-edit rows.
Twenty-four selected original captures were inspected. The first attempt failed before
qged because `/tmp/.X11-unix` had nonstandard host ownership; the retained replacement
used a private user/mount namespace and fresh artifacts. Private-Xvfb System GL is CPU
llvmpipe. Small-pane HUD overlap, offscreen edit framing, BoT surface-quality acceptance
and native GPU hosts remain S5 obligations.

Focused evidence, identities and verification live under
`.build/obol-qualification/20260909-root-repository`. ROOT-001 closure does not qualify
the downstream GED, source/service and presentation boundaries, known planning/frame
stalls, Lucy quality, broad resource limits or an exact release candidate. GED-005--013
below close the audited live GED target/return, compact-bootstrap, source-copy and core presentation
boundaries; the other obligations remain
in S1 and S3--S6.

## September 10: live GED source transitions

`CXX-GED-005--008` close five adjacent live GED publication defects. Base-source
promotion used to publish the new instance identity before its representation; its
callback observed `REPRESENTATION_DEFAULT`. External proxy replacement cleared and
invalidated the live source before the complete native source publication. Retained-
current marking added external ownership and realization status in two separately
observable controller calls. The source-record bridge applied one record through seven
separately observable controller calls, and its final view-policy update could invalidate
a record which had just been marked current. Broad-source redraw adopted compact
occurrences before retiring their narrow-source owners, so observers could see new payload
with old hierarchy, indexes and revisions while the GED caller retained a raw target
pointer across callbacks.

Base promotion now uses the source publisher's current-to-next identity overload, so
identity, representation, group state, indexes and frame effects commit together. The
external-to-native source publication marks the retained proxy stale while preserving it
as continuity geometry. A composed realization-state overload publishes status, revisions,
diagnostic, staleness, role ownership and owned-shape state under one notification envelope.
The retained-current helper calls that overload once. Callback edits remain current, and
callback replacement/removal cannot receive a later write intended for the old source.
The complete source publisher now accepts an optional explicit realization record. It
applies that record after configuration/view invalidation on the private candidate. The
GED record bridge translates configuration, display, material, realization, ownership
and view policy once, preserves existing appearance, and publishes current-to-current in
one scene transaction.

Resident-occurrence adoption now builds a detached target compact-index candidate while
the controller prepares all selected donor edges, the source index and affected revisions.
The target commits and donor owners retire before source or hierarchy observers, including
zero-change adoption. GED groups stable donor keys by target key and resolves retained
objects inside the transaction. The old payload-only public adoption operation was removed
once the live caller moved to the composed boundary.

The neighboring private local mesh/vlist factories had no source callers and no exported
dynamic symbols. Removing their common helper and declarations deletes about **320 lines**
and removes an otherwise unnecessary borrowed-target acceptance boundary.

Evidence is retained under `.build/obol-qualification/20260910-ged-live-callers`.
Four saved boundaries reject the new contracts: base promotion exits 1 with the exact
default representation, pre-deferred libraries expose partial external retirement and
lose the retained proxy transition, and pre-retained-current libraries expose split
status/ownership. The pre-source-record libraries expose the new revisions while
realization and roles remain old, then exit on that first partial callback. The resident-
adoption predecessor exports and calls the split operation but lacks the composed
controller method; the current contract fails to link against it on that missing symbol.
The current GED regression passes base promotion, retained proxy invalidation, retained-
current, source-record callback change, same-key replacement and removal, and the natural
narrow-to-broad compact occurrence and selection transition.
The source-level edge test covers external retirement; the scene-state edge test covers
combined realization ownership, a throwing observer and a later reentrant edit.

The final artifact runs **seven focused checks in 24.41 seconds**: source and source/group
publication allocation sweeps, source-publish edges, scene-state edges, resident-adoption
allocation and callback edges, and the complete GED draw-sync regression. Source
transition sweeps cover **465 positions** (456 preserved, eight committed) and source/group
transition sweeps cover **495 positions** (486 preserved, eight committed). Resident
adoption covers **239 positions** (236 preserved, two
committed), plus complete-state, edit, replacement, removal, throwing-observer and invalid-
target edges. The final compact-test hash is `a6d07a4e...`; GED draw-sync is
`2bebb3aa...`. CTest's GED contract row also passes in 7.95 seconds. This focused
checkpoint does not rerun the full client, native, graphical,
resource or platform matrix. At this GED-008 checkpoint, the live deferred adapters were
the next acceptance boundary. S1--S6 and the stage estimates remained unchanged.

## September 10: deferred realization acceptance

GED-009 closes the live deferred-result target/reference boundary. Incremental occurrence
batches now enter through `mergeDatabaseSourceInstanceCompactOccurrences`; terminal
results enter through `adoptDatabaseSourceInstanceRealization`. Both authenticate the
captured source stamp, retain the exact source through callbacks and publish affected frame
effects before observers. Coverage, batch and terminal completion witnesses distinguish
precommit allocation failure from callback failure after a complete commit. The GED pump
stores keys across callbacks and re-resolves the exact routing identity before continuing.
Terminal proxy cleanup captures parent/ancestor/routing identity before the callback and
revalidates each descendant before removal.

Evidence is retained under
`.build/obol-qualification/20260910-ged-deferred-acceptance`. Seven checks pass in **16.13
seconds**: direct streamed terminal publication, scene terminal allocation/callback checks,
scene stream allocation/callback checks, the complete GED draw-sync regression and its
CTest row. Terminal scene publication covers **135 allocation positions** (84 preserved,
50 complete commits); streamed scene publication covers **88** (39 preserved, 39
consistent prefixes, nine completed batches). The live regression proves complete
coverage and terminal callbacks retain later edits. Compact-test SHA-256 is `97e877e5...`;
GED draw-sync is `e3e3fe66...`; libBObol is `8f0c5a67...`; libged is `85e14c1c...`.
The saved pre-GED-008 library lacks both controller APIs, which is useful symbol evidence
but not an immediate-predecessor behavioral baseline. This focused boundary does not rerun
the full client, native, graphical, resource or platform matrix. S1--S6 and the estimates
below remain unchanged.

## September 10: compact proxy snapshot publication

GED-010 closes the split compact-bootstrap writer. Cached overview, decoded manifest,
synchronous structural and view-envelope paths now call
`publishDatabaseSourceInstanceCompactSnapshot` with an exact source stamp. The operation
builds its compact index off-scene, then publishes the registry, optional certified
bounds/profile, realized/current record, external ownership, primary-child retirement and
scene frame effect as one prepared source transaction. Auxiliary children remain. Its
commit witness distinguishes preparation failure from observer allocation failure after
the snapshot is complete. GED retains the latter result and does not write through a
borrowed source afterward. The redundant AABB bounds setter was removed because external
line publication already commits the same exact derived bounds. Deferred roots consume a
cached AABB directly in the compact transaction instead of publishing an external node
which the next call would immediately retire.

Evidence is retained under
`.build/obol-qualification/20260910-ged-compact-snapshot`. Eight checks pass in **16.50
seconds**. Compact snapshot publication covers **263 allocation positions** (201 preserved,
61 complete commits), with later edit, same-key replacement, removal, throwing observer,
invalid target and stale-stamp cases. Terminal and stream scene checks remain green. The
full GED draw-sync regression and its CTest row pass; the live deferred draw changes line
width from the compact-replacement callback, throws `std::bad_alloc`, and retains both the
snapshot and later edit. Compact-test SHA-256 is `e4394ab5...`; GED draw-sync is
`562020e2...`; libBObol is `b490b7fd...`; libged is `dd56cc6c...`. The retained GED-009
symbol log lacks the new API while current libBObol exports it and current libged calls it.
This focused checkpoint does not rerun the full client, native, graphical, resource or
platform matrix. S1--S6 and the stage estimates remain unchanged.

## September 10: retained-current source sequence composed

GED-011 closes the enclosing sequence left outside GED-006. A live source-policy
publication could mark retained external geometry stale and notify before a second helper
restored its current realization record. That helper also looked up the instance again,
allowing its later write to reach callback-time same-key replacement geometry.

The GED publisher now carries the retained external roles and final current realization
record in the same `BObolDatabaseSourcePublishState` as configuration, view policy,
appearance, group state and placement. Candidate invalidation precedes the authoritative
realization record off-scene; one commit then exposes the complete source. The follow-up
lookup, setter, helper implementation and private declaration are removed.

Evidence is retained under
`.build/obol-qualification/20260910-ged-retained-sequence`. The current regression passes
display edit, removal and stale same-key replacement cases through the real GED source
ensure path. Linked against the saved GED-010 archive, it exits on the first partial
callback. That archive has the follow-up helper symbol; the current libged does not. Four
focused checks pass in **23.60 seconds**, including all seven source-publication allocation
scenarios, source edges, full GED draw sync and its CTest row. GED draw-sync SHA-256 is
`0137776b...`; libged is `93f39cac...`; the saved libged baseline is `dd56cc6c...`.
This focused checkpoint leaves the stage estimates unchanged.

## September 10: cross-scene source copy composed

GED-012 closes a two-step secondary-scene insertion. The copy helper published identity,
representation, revisions, display, roles, view policy and placement, then applied material
policy through a separate instance lookup. Root callbacks could see the default inherited
policy and replace the source before the saved policy write.

The source summary translator now includes `materialPolicyValid` and `materialPolicy` in
the existing complete publication. The second setter and result merge are gone. A live
root-observer regression invokes the reducer draw path with a separate scene controller;
the callback sees the complete source and its edit, removal or same-key replacement remains
current.

Evidence is retained under `.build/obol-qualification/20260910-ged-scene-copy`. The current
test passes all three callback cases; the same test object linked with the saved GED-011
library exits on the partial material-policy callback. Four focused checks pass in **23.90
seconds**, including the full seven-case source publication allocation sweep and the GED
CTest row. GED draw-sync SHA-256 is `90acdd8e...`; libged is `2be67cea...`; the saved libged
baseline is `93f39cac...`. The stage estimates remain unchanged.

## September 10: GED presentation publication composed

GED-013 closes the path/global visibility, highlight and transparency writer family. The
old reducer paths published compact retained rules, source aggregate state and groups in
separate calls. A source callback could observe a split presentation, and callback changes
to later sources or groups could be overwritten. Global highlight clear also removed
retained compact rules outside the aggregate publication.

The scene controller now accepts one `BObolScenePresentationTransaction`. Each source
prepares and commits its display fields, retained compact rules, effective occurrence
visibility/highlight/transparency, renderer records, sparse visibility effects and frame
effect before notifying. Groups use the existing prepared scalar publisher. All targets
are captured before publication; invalid/duplicate targets reject before the first commit,
and callback edits, removal or replacement supersede pending work. Multi-target failure
therefore leaves a consistent prefix of complete source/group records.

Evidence is retained under `.build/obol-qualification/20260910-ged-presentation`. The
isolated current GED contract linked with saved GED-011 libraries exits 1 on the split
highlight callback; GED-012 did not change this presentation path. Current highlight and
global clear callbacks see aggregate state, every compact occurrence and frame effects
together, and retain a callback line-width edit. The scene test additionally verifies
retained leaf visibility/highlight/transparency and renderer opacity/visibility, pending
source/group edits, prevalidation and **106 allocation positions**: 71 preserve the old
scene and 36 leave complete prefixes. Five focused checks pass in **16.25 seconds**,
including adjacent retained-presentation/display edges, the full GED draw-sync executable
and its CTest row. GED draw-sync SHA-256 is `6593e5a6...`; libBObol is `9072d9ae...`;
libged is `862f0ab1...`; saved libged is `93f39cac...`. This focused checkpoint leaves the
stage estimates unchanged.

## September 10: public source-display publication composed

GED-014 closes the public shape/group-reference display row. The source bridge previously
published draw mode, representation and display/material fields separately. Public shape
references then repeated the update through a database-path alias. A callback could see a
partial source record, and either the remaining field writes or the alias replay could
overwrite its edit.

The bridge now translates a retained source summary into one complete publication state.
All same-path targets are retained and snapshotted before the first commit, and each pending
target must still have the captured node identity. The public shape path canonicalizes its
database path before publishing once. The existing public group branch remains one prepared
group display publication. The shared summary translator also removes the duplicated
cross-scene copy mapping.

Evidence is retained under `.build/obol-qualification/20260910-ged-source-display`. The
saved GED-013 library exits 1 because the public shape alias replay overwrites a callback
visibility edit; current code passes public shape/group callbacks, a simultaneous
draw/representation/display/material update, and pending-target edit, same-key replacement
and removal across two same-path sources. Five focused checks pass in **21.73 seconds**,
including the complete scene source-publication sweep, adjacent display edges, full GED
draw sync and its CTest row. GED draw-sync SHA-256 is `dea9967d...`; libBObol is
`9072d9ae...`; libged is `fb45acd0...`; saved libged is `862f0ab1...`. This focused
checkpoint leaves the stage estimates unchanged.

## September 10: deferred-selection target acceptance composed

GED-015 closes the deferred-selection initialization row. Deferred startup previously
validated a source, synchronized compact selection and then continued using the borrowed
source pointer. Because selection synchronization notifies, a callback could retarget,
replace or remove that source before its worker state was captured.

Startup now retains the source across selection callbacks and repeats one shared target
predicate afterward. The instance key must still resolve to the retained node, and its
routing identity, path and representation must still match the captured target before the
realization stamp, database snapshot and detached template are created.

Evidence is retained under
`.build/obol-qualification/20260910-ged-deferred-selection`. The saved GED-014 library
exits 1 with `sourcePreparationPending=1` after the callback retarget; current code
preserves retarget, same-key replacement and removal and launches no superseded work.
Five focused checks pass in **21.62 seconds**, including the complete scene
source-publication sweep, adjacent display edges, full GED draw sync and its CTest row.
GED draw-sync SHA-256 is `db51df3d...`; libBObol is `9072d9ae...`; libged is
`fa8d45d0...`; saved libged is `fb45acd0...`. This focused checkpoint leaves the stage
estimates unchanged.

## September 10: compact subtract style composed with realization

GED-016 closes the compact subtract-style row. Direct realization carried a transient
dashed flag into compact construction, while streamed occurrences retained only their
boolean operation. Stream reconstruction therefore published subtractive wires with a
solid style. A later eager-redraw setter was a no-op for direct realization and could not
repair the streamed path.

Compact construction now derives the default dashed style from either the direct flag or
the durable subtract operation. Direct, streamed and reconstructed occurrences publish
their boolean identity, styles and renderer records together. The redundant GED
post-realization writer is removed; explicit later style overrides retain their existing
compact API.

Evidence is retained under `.build/obol-qualification/20260910-ged-subtract-style`. The
saved GED-015 libraries exit 1 with `/progressive_root.c/box.s` at subtract operation and
`style=0`; current code publishes dashed style. Six focused checks pass in **18.54
seconds**, covering the real deferred producer, direct compact callback, sparse style
override/allocation behavior, GED-015, full GED draw sync and its CTest row. GED draw-sync
SHA-256 is `163ef1df...`; compact publication is `3a49734c...`; libBObol is
`81e3b4e3...`; libged is `376a782f...`; saved libraries are `9072d9ae...` and
`fa8d45d0...`. This focused checkpoint leaves the stage estimates unchanged.

## September 10: database rename composed across the live scene

GED-017 closes the rename and compact-path-retargeting row. The reducer previously renamed
repository keys, compact paths, individual sources and groups in separate notification
windows, then restored saved source state. Observers could see mixed old/new identities and
indexes, callback edits could be overwritten, and a nested object rename could collapse the
retained source path to the renamed leaf.

The scene controller now prepares every affected repository, source identity and owned
shape, descendant parent key, compact semantic and renderer record, group path/name/intent,
index and peer revision as one operation. Path-derived GED keys retarget explicitly; opaque
keys stay stable. The committed scene is realized only after observer delivery. This also
stops active realization in every affected controller and includes controllers that share a
source beneath otherwise independent roots.

Evidence is retained under `.build/obol-qualification/20260910-ged-rename`. The saved
GED-016 libraries exit 1 because the current callback contract observes incomplete state;
the current binaries pass and preserve a callback display edit and nested component path.
The direct sweep covers **658 allocation positions**: 649 preserve the old scene and eight
fail after a complete commit. It covers same-root and split-root peers, three independent
repositories, opaque source identity, descendant hierarchy, stable compact handles and
renderer-parent rebuild. Callback edit, removal, exception and invalid-input cases pass.
Ten focused checks pass in **36.84 seconds**, including scene configuration, group and
repository rename regressions, GED-015/016, full GED draw sync and its CTest row. GED
draw-sync SHA-256 is `03e41974...`; compact publication is `1a246212...`; libBObol is
`2bac97bf...`; libged is `dc76e322...`; the saved libraries are `81e3b4e3...` and
`376a782f...`. This focused checkpoint leaves the stage estimates unchanged.

## September 10: redraw and erase target acceptance composed

GED-018 closes the redraw/erase inventory row. Redraw previously walked a live source array
which callbacks could reorder, then used global realization for a scoped request. Erase
wrote immediate authored compact visibility, removed sources and groups in separate
notifications, and later rescanned retained sources after removal callbacks. Those paths
could skip an accepted source, realize unrelated work, expose partial hierarchy, lose
visibility during compact population replacement, or overwrite a newer callback edit.

Redraw now captures stable source keys and exact realization stamps before its first
callback, including wildcard-mode requests. Scoped realization accepts only the current
key/stamp, limits compact descendant cleanup to that owner and leaves sibling representations
untouched. Exact and prefix erase publish all selected source/group
edges through one removal transaction. Retained compact visibility uses the durable
presentation override owner. A presentation stamp captured before removal rejects later
callback edits and same-key replacements, and only compact sources which cover the erased
path receive the rule. Pure removal no longer launches global realization.

Evidence is retained under
`.build/obol-qualification/20260910-ged-redraw-erase`. Nine focused checks pass in **17.22
seconds**. The scene presentation sweep covers **106 allocation positions** (71 preserved,
36 complete prefixes), and removal covers **206** (197 preserved, eight complete commits).
Realization coverage includes 88 callback histories, 296 descendant and 248 cleanup
allocation positions. The real GED regression covers compact-root retirement, callback
order shifts, scoped realization, same-key replacement, mode-scoped erase, newer retained
visibility and unrelated pending work after prefix erase. Full GED draw-sync and its CTest
row pass. GED draw-sync SHA-256 is `b83e6041...`; compact publication is `7e26506f...`;
libBObol is `6a348118...`; libged is `77226d29...`. No immediate predecessor was sealed;
the saved GED-016 library lacks all four new source/controller entry points and supplies API
evidence only. This focused checkpoint leaves the stage estimates unchanged.

## September 10: mesh-LoD attachment and bounds ownership composed

GED-019 closes the mesh-LoD attachment inventory row. The BoT path previously gave its raw
reader to the first semantic-path source before verifying the initial cut, then wrote bounds
through a separate all-matching-sources loop. Failure could leave a half-prepared owned
reader, sibling representations could receive bounds without its handle, and source identity
invalidation retained stale cache authority.

Preparation now owns the reader locally until hierarchy, cut and finite ordered bounds are
ready. One controller operation accepts the exact source instance and pre-work realization
stamp, then transfers reader and bounds as one record. Rejected work stays caller-owned. The
new record is installed before the old reader is destroyed. Source, input and database
identity changes retire the complete record; view-only changes retain it. Because this is
private service state, attachment does not advance a frame. The following mesh publication
owns the visible scene change. The former split public source setters and path-based GED
setters are gone.

Evidence is retained under `.build/obol-qualification/20260910-ged-mesh-lod`. Six focused
checks pass in **22.26 seconds**. The live GED check covers invalid bounds, atomic transfer,
stale rejection, retirement/rebuild and exact wire/shaded ownership. Direct source
invalidation/configuration sweeps, mesh-cache CTest, full GED draw sync and its CTest row pass.
GED draw-sync SHA-256 is `67ccd04f...`; compact publication is `11087cf9...`; libBObol is
`96b26047...`; libged is `bc203a86...`. The saved GED-016 libraries expose the removed split
symbols and lack the exact transfer, providing API evidence. The AddressSanitizer tree is
currently blocked by stale staged Obol headers; sanitizer qualification remains in S5. This
focused checkpoint leaves the stage estimates unchanged.

## September 10: immutable stream discovery and stamped reservation

GED-020 closes the last row in the finite GED/source presentation writer inventory. Compact
discovery previously used independent public source writers for expected count, profile and
capacity. Count publication silently selected a maximum, profiles were not checked for
internal consistency, parallel discovery published before terminal-empty leaves were
resolved, and capacity reservation was outside exact stamped delivery. Resource identity
changes could therefore retain stale discovery authority, while denied allocation could
leave a partial source-side record.

The stream now publishes one exact immutable count and optional validated profile after final
discovery. Conflicting and post-cancellation writes are rejected. A stamped scene-controller
operation certifies those facts for one current source, and capacity reservation occurs
inside the exact batch-delivery operation with strong failure preservation. Source, input and
database identity invalidation revoke the contract. View-only invalidation and failed
delivery preserve valid discovery evidence for the retained preview. The three split public
source setters are removed.

Evidence is retained under
`.build/obol-qualification/20260910-ged-stream-contract`. Nine focused checks pass in
**26.22 seconds**. They cover malformed/conflicting/stale facts, extreme-size certification,
capacity denial, retained geometry, stream concurrency/cancellation, incremental and
terminal publication, source invalidation/configuration allocation sweeps, LoD append
traversal, the live GED denial path, full GED draw sync and its CTest row. GED draw-sync
SHA-256 is `6f04833e...`; compact publication is `08b9e6c4...`; libBObol is `d7d092a5...`;
libged is `08f3f06a...`. The saved GED-019 library exposes the removed setters and lacks the
new stamped controller operation. This focused checkpoint leaves the stage estimates
unchanged.

## September 10: direct source-field publication

`CXX-CONFIG-008` closes the twelve watched source fields which previously
published a new field value before stale state, owned-shape metadata and source
resource retirement. The retained baseline under
`.build/obol-qualification/20260910-direct-field-sensors` exits 1 after observing
that partial state on a raw `path` assignment.

The source now intercepts those field notifications synchronously, establishes
an allocation-free fail-safe invalidation, and then uses the existing prepared
configuration publisher before inherited field/node delivery. Internal prepared
notifications carry a one-shot marker, so they do not recursively invalidate the
source and callback-time raw writes remain independent events. Twelve private
no-op auditors keep notification live when the source container is quiet; they
perform no state mutation. Direct fields disabled individually retain ordinary
Obol semantics and publish after re-enable plus `touch()`.

Evidence is retained under
`.build/obol-qualification/20260910-direct-field-publication`. The focused row
passes all twelve immediate/delayed fields, quiet source/field edges, field and
engine connections, reentry, observer exception recovery and 87 allocation
positions (one preserved, 85 explicit fail-safe, one committed). The complete
state suite passes 80 checks and 805 scenario/sweep rows in **172.74 seconds**.
The compatible external-draw matrix passes 220 allocation positions across four
mesh/quiet cases, GED draw-sync CTest passes in 8.21 seconds, and the static-link
probe passes. Current compact-test, static libBObol and shared libBObol SHA-256
values are `fe4a92a2...`, `3ca4cd0a...` and `651e2d39...`; the saved failing
libBObol is `d7d092a5...`.

This closes the direct field/engine inventory row at its C++ publication
boundary. It does not qualify graphical behavior, native hosts, shared-stack
sanitizers, resource pressure or a final release candidate.

## September 10: throwing deletion callbacks

`CXX-SOURCE-026` closes the native delete-callback exception boundary left by
the earlier child/path retirement work. A standalone probe attaches a throwing
delete callback to a node, path, node-owned field and parent-owned child. The
preceding Obol library propagates node/path deletion before the object is
retired (exit 1), while field and child destruction cross a `noexcept`
destructor and reach the terminate handler (exit 42).

`SoDataSensor::invokeDeleteCallback` now contains exceptions from the user
hook, allowing every sensor to detach and irreversible object retirement to
finish. The public callback documentation makes this behavior explicit.
Ordinary sensor change callbacks retain their existing exception behavior. A
permanent native test verifies one callback, detachment and actual destructor
completion for all four object families.

Evidence is retained under
`.build/obol-qualification/20260910-delete-callbacks`. All four current probe
modes exit 0. Native Obol passes 1,578 of 1,579 unit tests with one established
skip, all 63 integration tests and all three lifecycle tests. The complete
libBObol compact-publication executable passes 126 direct checks, 21 scenarios
and 1,021 pass rows in **309.79 seconds**. Eight focused child, sensor, path,
graph, auditor, direct-field and custom-traversal rows pass as well. Current and
baseline Obol SHA-256 values are `6d194213...` and `acbe7137...`.

Custom-node cleanup and legacy individual remove/truncate callers remain in S1
and S3. This focused closure does not change the dated completion ranges below.

## September 10: source realized-geometry clear

`CXX-SOURCE-027` closes `SoBRLDatabaseSource::clearRealizedGeometry`. The
preceding implementation retired compact/compiled state and then removed
primary children individually. An immediate source or path observer therefore
saw compact state absent with only a prefix of the old child graph removed; an
observer exception could stop the rest of the clear.

The source now prepares one complete `SoChildList::Removal`, commits children
and paths quietly, retires compact and compiled state, and only then attempts
path and source notifications. Preparation failure preserves the old source;
callback failure leaves the complete clear committed and does not suppress the
remaining notification attempts. Placement, auxiliary policy, terminal fields
and the established return behavior have explicit regression coverage.

Evidence is retained under
`.build/obol-qualification/20260910-clear-realized-geometry`. The baseline fails
the new immediate-observer row. Current coverage passes both auxiliary modes and
24 allocation positions (21 preserved, three committed), throwing observers,
exact path truncation/reindexing, retry, and an allocation-free already-clear
child graph. Eight neighboring publication/lifetime rows pass. The full compact
executable passes 127 direct checks, 21 scenarios and 1,022 pass rows in
**310.16 seconds**. Prototype and GED draw-sync CTest rows pass in 0.80 and 8.08
seconds; GED covers redraw through this API in all six display modes. Current
compact-test/shared/static libBObol hashes are `0414f3ea...`, `0384541c...` and
`5ab79953...`; the corresponding baseline values are `d02f6606...`,
`651e2d39...` and `3ca4cd0a...`.

This closes one concrete legacy removal caller. Custom-node rebuilds and the
remaining live detach/reparent callers stay in S1/S3, so the completion ranges
below remain unchanged.

## September 10: axes and ADC overlay rebuild publication

`CXX-SOURCE-028` closes the axes/ADC part of the custom-node rebuild family.
Both nodes previously removed their children before allocating and configuring
the replacement shapes. Immediate node/path observers could see an empty or
partial overlay, and allocation failure destroyed the preceding presentation.

The rebuilds now retain newly built shapes locally, prepare one complete child
replacement, commit children and audited paths quietly, and notify only after
the complete graph is visible. Hidden overlays publish one complete clear; an
already hidden and empty overlay is allocation- and notification-free.
Preparation failure preserves the old graph, while callback failure leaves the
new graph committed and does not suppress later notifications.

This repair also found a native lifetime defect: prepared child edits deleted a
valid zero-reference parent when releasing their temporary reference. Obol now
preserves the caller-owned zero-reference state while retaining normal final
release for a parent that began owned. A native regression covers committed
replacement, committed removal and abandoned replacement.

Evidence is retained under
`.build/obol-qualification/20260910-overlay-rebuild-publication`. The baseline
libBObol fails because an observer sees a partial graph; repaired libBObol with
the preceding Obol separately fails the prototype axes-segment check. Current
focused coverage passes axes/ADC replacement and clearing across 67 allocation
positions (63 preserved and four committed), throwing observers, exact path
retirement, retry and no-op clearing. Seven neighboring rows pass. The full
compact executable passes 128 direct checks, 21 scenarios and 1,023 pass rows
in **308.17 seconds**. Prototype and Qt faceplate CTest rows pass in 0.80 and
0.04 seconds. Native Obol passes 1,579 of 1,580 unit tests with one established
skip, all 63 integration tests and all three lifecycle tests. Current compact,
shared/static libBObol and Obol hashes are `1878b08b...`, `63565e32...`,
`b4d180c8...` and `21872f7f...`.

Grid and the other custom-node builders remain open, followed by live
detach/reparent callers. This bounded subfamily does not materially change the
completion ranges below.

## September 10: grid rebuild publication

`CXX-SOURCE-029` closes `SoBRLGrid::rebuildGeometry`. The preceding method
removed its HUD child and published seven derived spacing/count values in
separate steps before the replacement graph existed. Immediate node, field and
path observers could see mixed generations, while allocation failure destroyed
the preceding grid.

The rebuild now calculates all spacing and line results locally, owns every new
shape and HUD node through preparation, and publishes one child replacement
together with one prepared seven-field notification set. The live graph and
fields commit with notifications disabled; callbacks run only after the whole
result is visible. Failure before commit preserves the old graph and fields.
Failure during notification leaves the new state complete, drains the remaining
callbacks and then rethrows. Hidden clearing uses the same path, with a true
already-empty no-op.

Evidence is retained under
`.build/obol-qualification/20260910-grid-rebuild-publication`. The baseline
fails the immediate-observer regression. Current focused coverage passes 156
allocation positions (144 preserved and 12 committed), visible replacement,
clearing, all seven fields, exact old-child/path ownership, retry and throwing
observers. Eight neighboring rows pass. The full compact executable passes 129
direct checks, 21 scenarios and 1,024 pass rows in **307.59 seconds**. Prototype
and Qt faceplate CTest rows pass in 0.80 and 0.04 seconds. Current compact,
shared/static libBObol and Obol hashes are `1c62689b...`, `0e6ea03d...`,
`aadd76d9...` and `21872f7f...`.

The bulk `bobol_grid_configure_from_view` input writer remains a separate
composition boundary. The other custom-node and live detach/reparent families
remain open, so the completion ranges below are unchanged.

## September 10: image display rebuild publication

`CXX-SOURCE-030` closes image-plane and viewport-image rebuild publication. The
preceding methods cleared children and cached texture/face pointers before
stream loading and replacement construction, exposing an empty or partial
display and losing the old image on failure.

Both nodes now retain their preceding presentation until the payload, fit,
textured quad and complete child replacement are ready. The common quad builder
owns its root and every intermediate node during configuration. A shared
publication helper commits the child graph, cached pointers and two realized
revision fields quietly, restores notification policy and then drains all
callbacks. Precommit failure preserves the complete old display. Callback
failure leaves the complete new display installed. Hidden clearing retains the
established revisions and has an already-empty no-op.

The referenced image source may notify while payload loading refreshes it. Those
preflight callbacks see a complete old display; the regression accepts that
state or the complete new display and rejects intermediate combinations.

Evidence is retained under
`.build/obol-qualification/20260910-image-rebuild-publication`. The baseline
fails the immediate-observer row. Current focused coverage passes four
plane/viewport replacement/clear cases across 289 allocation positions (280
preserved and nine committed), cached pointer and revision checks, dependency
preflight, retry and throwing observers. Eight neighboring rows pass. The full
compact executable passes 130 direct checks, 21 scenarios and 1,025 pass rows
in **308.81 seconds**. Image-display, prototype and Qt faceplate CTest rows pass
in 0.01, 0.82 and 0.04 seconds. Current compact, shared/static libBObol,
libimgstream, Obol and OSMesa hashes are `5106a469...`, `3f7e63dd...`,
`683e6bb0...`, `e2019786...`, `21872f7f...` and `097aac27...`.

CMake restaged a different OSMesa binary while adding the compact test's
explicit libimgstream dependency. The hash change was detected before testing,
the coherent qualified dependency was restored, and the retained loader trace
shows all four relevant libraries resolve from the evidence directory.

Image-source mutation and the remaining custom-node and detach/reparent families
stay open, so the completion ranges below are unchanged.

## September 10: HUD and line-layer rebuild publication

`CXX-SOURCE-031` closes HUD-label and line-layer overlay rebuild publication.
The preceding HUD-label method cleared its cached pointer and live graph before
replacement allocation. The preceding line-layer method removed the old graph,
attached a depth node and then added replacement shapes individually. Immediate
node and path observers could see empty or partial results, and allocation
failure lost the complete preceding overlay.

Both rebuilds now construct every replacement node under reference ownership
outside the live child list and prepare one complete replacement before commit.
The HUD-label method publishes its cached pointer before callbacks. The line-
layer method installs the depth state and every usable shape together. Hidden
clearing uses the same prepared operation, and an already empty overlay performs
no allocation or notification. Failure before commit preserves the old graph
and audited paths; callback failure leaves the complete new graph installed and
drains the remaining notifications before rethrowing.

Evidence is retained under
`.build/obol-qualification/20260910-hud-line-rebuild-publication`. The exact
retained baseline fails the immediate-observer row. Current focused coverage
passes both overlay types under replacement and clearing across 142 allocation
positions (138 preserved and four committed), cached-pointer/content checks,
retry and throwing observers. Eight neighboring rows pass. The full compact
executable passes 131 direct checks, 21 scenarios and 1,026 pass rows in
**308.60 seconds**. Its source-only subset passes 51 direct checks and 221 pass
rows in 132.38 seconds. Image-display, prototype and Qt faceplate CTest rows
pass in 0.01, 0.79 and 0.04 seconds. Current compact, shared/static libBObol,
libbg, libimgstream, Obol and OSMesa hashes are `c1eea7b8...`, `76bce770...`,
`54352ccc...`, `87994660...`, `e2019786...`, `21872f7f...` and `097aac27...`.
The retained loader trace resolves all five relevant shared libraries from the
evidence directory.

Image-source mutation, edit preview/manipulator, navigation, cutting-plane and
scene-light builders, the enclosing grid writer and live detach/reparent family
stay open. This bounded closure does not change the completion ranges below.

## September 10: edit manipulator rebuild publication

`CXX-SOURCE-032` closes the child-publication boundary in the axis and indexed
edit-manipulator `rebuildGeometry` methods. Both previously removed their live
children and installed replacement light, axis, edge, point and emphasis
subtrees incrementally. Immediate node and path observers could see empty or
partial graphs, while allocation failure discarded the preceding manipulator.

Both methods now construct every top-level and nested node under reference
ownership outside the live graph. One shared helper prepares and commits the
complete top-level order, then notifies only after every recursive subtree is
complete. Hidden or empty rebuilds use the same prepared replacement and have
an already-empty no-op. Preparation failure preserves the old graph and paths;
throwing callbacks leave the complete new graph installed and do not prevent
later notifications.

Evidence is retained under
`.build/obol-qualification/20260910-edit-manipulator-rebuild-publication`.
The exact retained baseline fails the immediate-observer row. In an isolated
process, focused coverage passes both manipulator types under replacement and
clearing across 449 allocation positions (444 preserved and five committed).
After the preceding compact rows initialize process-global notification state,
the same check has 448 positions (444 preserved and four committed). Six
neighboring rows pass. The full compact executable passes 132 direct checks,
21 scenarios and 1,027 pass rows in **307.53 seconds**. Edit-manipulator,
image-display, prototype and Qt faceplate CTest rows pass in 0.01, 0.01, 0.79
and 0.04 seconds. Current compact, shared/static libBObol, Obol and OSMesa
hashes are `59af6aad...`, `d9ce205b...`, `897b62d3...`, `21872f7f...` and
`097aac27...`. The retained loader trace resolves the preserved dependency
stack.

The public-field and `setEllipsoidAxes`/`setTopology` composition boundaries
remain open, as do edit preview, image-source mutation, navigation, cutting-
plane, scene-light and live detach/reparent work. This bounded closure does not
change the completion ranges below.

## September 10: edit preview publication

`CXX-SOURCE-033` closes the output boundary in `setLineSet`,
`setTransformedLineSet` and `clearPreview`. Those methods previously removed
the live graph before replacement allocation and then wrote realized revisions,
preview status and stale state sequentially. Immediate node, field and path
observers could see mixed generations, and failure discarded the old preview.

A shared helper now prepares the complete child replacement and a four-field
notification set before live mutation. It commits the graph, realized revision
pair, status and stale state quietly, restores notification policy and drains
all callbacks after the complete result is visible. Direct and transformed
graphs own every node during construction. Invalid input preserves the last
realized pair while publishing empty FAILED/stale output; explicit clear
advances the pair and publishes EMPTY/non-stale. Repeated empty outcomes do no
work.

Evidence is retained under
`.build/obol-qualification/20260910-edit-preview-publication`. The exact
retained baseline fails the immediate-observer row. Isolated focused coverage
passes direct, transformed, clear and invalid outcomes across 96 allocation
positions (76 preserved and 20 committed). After preceding compact checks
initialize process-global notification state, the row has 93 positions (76
preserved and 17 committed). Six neighboring rows pass. The full compact
executable passes 133 direct checks, 21 scenarios and 1,028 pass rows in
**306.86 seconds**. Prototype, view-store, compact-edit, edit-manipulator and
Qt faceplate CTest rows pass in 0.83, 0.02, 0.01, 0.01 and 0.04 seconds.
Current compact, shared/static libBObol, Obol and OSMesa hashes are
`3db95b8f...`, `3396f728...`, `441f7d2d...`, `21872f7f...` and
`097aac27...`. The retained loader trace resolves the preserved dependency
stack.

Input publication through `setEditIntent`, `markSourceRevision` and
`markInputsRevision` remains separate. Manipulator setters, image-source
mutation, navigation, cutting-plane, scene-light and live detach/reparent work
also remain open. This bounded closure does not change the completion ranges
below.

## September 10: navigation-gizmo rebuild publication

`CXX-SOURCE-034` closes direct `SoBRLNavigationGizmo::rebuildGeometry`
publication. The preceding implementation cleared its cached HUD and anchor
pointers and removed the live graph before allocating the replacement. It then
built and attached every nested node through raw pointers. Immediate observers
could see empty state, while allocation failure discarded the preceding gizmo
and could leak an incomplete private subtree.

The complete HUD now builds offline under reference ownership. One prepared
child replacement commits the new graph and both cached pointers before path
and node callbacks. Hidden clearing uses the same operation and an already empty
hidden graph does no work. Failure before commit preserves the exact old HUD,
cached HUD identity and audited path; callback failure leaves the complete new
graph visible and allows later callbacks to run.

Evidence is retained under
`.build/obol-qualification/20260910-navigation-gizmo-rebuild-publication`. The
current test executable forced to load the exact saved predecessor fails the
immediate-observer row. Isolated focused coverage passes replacement and clear
across 949 allocation positions (945 preserved and four committed), including
the recursive 20-child circles widget, cached HUD identity, retry, throwing
observers and repeated hidden no-op behavior. In the full process the same row
has 948 positions (945 preserved and three committed). Five neighboring graph
rows and the dedicated navigation CTest pass. The full compact executable passes
134 direct checks, 21 scenarios and 1,029 pass rows in **304.57 seconds**.
Current compact, navigation-test, shared/static libBObol, Obol and OSMesa hashes
are `15269f7c...`, `3263867d...`, `d286ae7a...`, `fbf312a4...`, `21872f7f...`
and `097aac27...`. The retained loader trace resolves libBObol, libbg,
libimgstream, Obol and OSMesa from the evidence directory.

Camera synchronization and hover/active setters remain input-plus-rebuild
composition boundaries. Image-source mutation, grid configuration, preview and
manipulator input setters, cutting-plane, scene-light and live detach/reparent
work also remain open. This bounded closure does not change the completion
ranges below.

## September 10: cutting-affordance and scene-light rebuild publication

`CXX-SOURCE-035` closes two controller-owned output builders. The cutting-plane
affordance previously removed its HUD before plane/camera/viewport validation,
projection and replacement allocation. The scene-light path removed all current
lights and added directional, spot and point nodes one at a time. Its group and
camera-rig helpers also used observable detach/reinsert sequences on the
viewport root.

The affordance now computes its full projected line set and builds its HUD and
six-node widget offline under reference ownership. Existing output publishes
through one prepared child replacement; a new affordance attaches only when its
HUD is complete. Disabled/invalid output clears atomically and repeated empty
output does no work. Scene lights similarly build every typed node offline and
publish one complete group replacement. Camera-rig and scene-light group order
changes use prepared root orders and unchanged orders are no-ops.

Evidence is retained under
`.build/obol-qualification/20260910-cutting-scene-light-rebuild-publication`.
The current executable forced to load the exact saved predecessor fails both
immediate-observer rows. Cutting coverage passes 147 allocation positions (145
preserved and two committed), including 36 projected points, 18 line segments,
retained HUD/path, clearing, retry, throwing observers and disabled no-op
behavior. Scene-light coverage passes 80 positions (69 preserved and 11
committed), including all three light types, exact fields, retained children/
paths, clearing, retry, throwing observers and empty no-op behavior. Five
neighboring graph rows pass. Window-host and retained-raytrace CTests pass in
0.47 and 0.75 seconds. The full compact executable passes 136 direct checks, 21
scenarios and 1,031 pass rows in **305.13 seconds**. Current compact, window,
raytrace, shared/static libBObol, Obol and OSMesa hashes are `a9795135...`,
`f647efaa...`, `a68aca60...`, `e9b3fca3...`, `f348d6dc...`, `21872f7f...`
and `097aac27...`. The retained loader trace resolves libBObol, libbg,
libimgstream, Obol and OSMesa from the evidence directory.

Cutting-plane and scene-light input/enabled setters remain composition
boundaries. Image-source mutation, grid configuration, preview, manipulator and
navigation input setters, render-environment repair and live detach/reparent
work also remain open. This bounded closure does not change the completion
ranges below.

## September 10: render-environment repair publication

`CXX-SOURCE-036` closes live render-environment creation/repair and the paired
viewport-root order. The preceding repair inserted missing environment children
one at a time, detached/reinserted a misplaced environment, and repaired the
camera rig afterward. Immediate observers and allocation failure could retain
those partial states. Absent-environment creation could also abandon an
unowned migration headlight when a named root headlight already existed.

All missing nodes are now reference-owned before publication. Existing
extension children remain in order; both the environment and viewport-root
replacement operations prepare before either commits; both complete lists are
visible before notification starts. A new environment is assembled offline and
published once. Root order is deterministic with and without a camera, repeated
owned occurrences are normalized, notification exceptions do not suppress the
other prepared notification, and an already complete graph is an allocation-
free no-op.

Evidence is retained under
`.build/obol-qualification/20260910-render-environment-repair-publication`.
The current executable forced to load the exact saved predecessor fails with
the expected partial-children/root-order diagnostic. Focused coverage passes
147 allocation positions (131 exact preceding states and 16 complete commits)
across absent creation and in-place repair. It checks all required nodes and
defaults, unrelated extensions, every retained path, retry, throwing observers
and complete no-op behavior. Six neighboring publication rows pass. OSMesa-
render, window-host and retained-raytrace CTests pass in 0.20, 0.41 and 0.63
seconds. The retained loader trace resolves libBObol, libbg, libimgstream, Obol
and OSMesa from the evidence directory. In the full process the focused row
passes 145 positions (131 preserved and 14 committed), and the exact executable
passes 137 direct checks, 21 scenarios and 1,032 pass rows in **305.88
seconds**. Current compact, window, render, raytrace, shared/static libBObol,
Obol and OSMesa hashes are
`b97e2f64...`, `f647efaa...`, `3d0f9bcc...`, `a68aca60...`, `00375975...`,
`db31baa5...`, `21872f7f...` and `097aac27...`.

One `libBObol_lod_update_action` invocation run concurrently with two other
executables reported a missing current exact-frame witness. The saved
predecessor and four isolated current invocations pass. Treat that as a
load-sensitive S4 reproduction lead; it does not qualify or reject this
publication boundary.

This closes the lighting-side live detach/reparent inventory. Framebuffer host
and endpoint migration, feature-store attachment/reordering and Qt faceplate
adaptation remain open. The input-setter composition boundaries and the S4--S6
resource, platform and graphical gates are unchanged, so the conservative
completion ranges below remain current.

## September 10: feature-overlay order publication

`CXX-SOURCE-037` closes the pure feature-overlay root reorder step. The previous
store detached all feature overlays and reattached sorted records one by one,
exposing truncated roots and paths to immediate observers and failure. It also
depended on repeated attachment checks to keep a custom node shared by multiple
records unique.

The store now builds the complete successor order offline, preserves unrelated
endpoint children in relative order, reduces shared-node records to one root
occurrence, and publishes one prepared child-list replacement. Every retained
path remaps directly to its final index. An exact repeated overlay-setting call
allocates and notifies nothing.

Evidence is retained under
`.build/obol-qualification/20260910-feature-overlay-order-publication`. The
current executable forced to load the exact saved predecessor fails with the
expected partial-root diagnostic. Isolated coverage passes 129 allocation
positions (125 exact preceding states and four complete commits) across typed
HUD lines, an aliased custom node, unrelated extensions, every retained path,
retry, a throwing observer and the exact no-op. Six neighboring publication
rows pass. Feature-store and Qt faceplate CTests pass in 0.02 and 0.04 seconds.
The compact loader resolves five relevant libraries from the evidence directory;
the faceplate loader resolves seven. In the full process the focused row passes
128 positions (125 preserved and three committed), and the exact executable
passes 138 direct checks, 21 scenarios and 1,033 pass rows in **305.40 seconds**.
Current compact, feature-store, Qt faceplate, shared/static libBObol, Obol and
OSMesa hashes are `9d403eee...`, `b3617626...`, `c49c1629...`, `f6a048a4...`,
`61e51c22...`, `21872f7f...` and `097aac27...`.

Feature replacement/reparenting, multi-record controller migration, and the
enclosing metadata/node/order setter remain open. Framebuffer host, display-
endpoint and Qt faceplate live-child work also remains. This bounded closure
does not change the conservative completion ranges below.

## September 10: framebuffer composition publication

`CXX-SOURCE-038` closes composition changes for an already open window-host
framebuffer. The prior setter removed the live viewport from its layer roots,
changed fields, rebuilt realized children, and attached the destination in
separate observable steps. Failure could retain the detached viewport, and a
path through the realized image could retire before the intended layer existed.

The host now builds a private viewport candidate without refreshing or notifying
the shared live image source. Viewport scalar fields, realized children, cached
render-node pointers, and all three layer-root orders prepare first. They commit
as one complete state before notifications begin. Unrelated children keep their
relative order, owned viewport occurrences are normalized, every prepared
notification is attempted after callback failure, and an exact repeated mode is
an allocation-free no-op.

Evidence is retained under
`.build/obol-qualification/20260910-framebuffer-composition-publication`. The
current executable forced to load the exact saved predecessor fails with the
expected partial-viewport/root diagnostic. Isolated coverage passes 257
allocation positions (254 exact preceding states and three complete commits).
It checks all three roots, viewport and image-source identity, scalar/realized
state, cached texture/face pointers, direct paths, a path through a retired
realized child, retry, throwing observers and the exact no-op. Six neighboring
publication rows pass. Image-source, image-display, window-host, headless-host
and Qt window-host CTests pass in 0.01, 0.01, 0.42, 0.04 and 0.54 seconds. The
compact loader resolves five relevant libraries from the evidence directory;
the Qt loader resolves seven. In the full process the focused row passes 256
positions (254 preserved and two committed), and the exact executable passes
139 direct checks, 21 scenarios and 1,034 pass rows in **305.23 seconds**.
Current compact, window-host, Qt window-host, shared/static libBObol, Obol and
OSMesa hashes are `a7aa017e...`, `b35f3d74...`, `8ef3825b...`, `f70823e3...`,
`0bd218f9...`, `21872f7f...` and `097aac27...`.

Retained-raytrace display-endpoint migration is closed separately by
CXX-SOURCE-039 below. Framebuffer open/close cleanup, viewport/cursor setters,
feature-store reparenting and Qt faceplate composition remain open. This bounded
closure does not change the conservative completion ranges below.

## September 10: raytrace layer composition publication

`CXX-SOURCE-039` closes composition-layer migration for the retained-raytrace
display endpoint. The prior path removed the viewport, changed its layer,
rebuilt its image graph, and attached the destination in separate observable
steps. The property setter then restarted ray work even though the render
content had not changed.

The endpoint now shares the window host's complete three-root preparation. It
prepares the scalar layer change and every affected root, commits them together
with endpoint policy, and only then notifies. The retained RT image remains the
first child in its destination; unrelated layer children retain order. Its
source, realized children and cached texture/face nodes retain identity, and a
composition-only change requests a frame without restarting the renderer.
Every prepared notification is attempted after callback failure. An exact
repeat allocates and notifies nothing.

Evidence is retained under
`.build/obol-qualification/20260910-rt-layer-composition-publication`. The
current executable forced to load the exact saved predecessor fails because an
observer sees partial cross-root state. Isolated coverage passes 80 allocation
positions (78 exact preceding states and two complete commits). It checks all
three roots, endpoint policy, retained viewport/source/geometry identity,
direct and nested paths, retry, throwing observers and the exact no-op. Six
neighboring publication rows, retained-raytrace, window-host, headless-host and
Qt window-host tests, and the static link check pass. The compact loader resolves
five relevant libraries from the evidence directory; the Qt loader resolves
seven. The exact full executable passes 140 direct checks, 21 scenarios and
1,035 pass rows in **308.97 seconds**. Current compact, retained-raytrace,
shared/static libBObol, Obol and OSMesa hashes are `a896c5f9...`, `d084adb1...`,
`8e5e26b6...`, `e462056f...`, `21872f7f...` and `097aac27...`.

Endpoint presentation cleanup, framebuffer open cleanup, viewport/cursor
setters, feature-store reparenting and Qt faceplate composition remain open.
Explicit framebuffer close is closed separately by CXX-SOURCE-040 below. This
bounded closure does not change the conservative completion ranges below.

## September 10: framebuffer close publication

`CXX-SOURCE-040` closes explicit retirement of a window-host framebuffer. The
prior close path removed the live viewport and released viewport/source
ownership before erasing the host attachment. Immediate root and path observers
could see a detached graph while host count and lookup still reported the
framebuffer, and callback failure could retain that mixed state.

Close now prepares the complete three-root successor transaction before any
live change. It erases the attachment and commits all root/path changes before
notification, then releases retained object ownership and requests
presentation. Every cleanup and notification step preserves the first callback
failure while attempting the rest. Allocation failure during preparation keeps
the exact preceding attachment, and a repeated close allocates and notifies
nothing.

Evidence is retained under
`.build/obol-qualification/20260910-framebuffer-close-publication`. The current
executable forced to load the exact saved predecessor fails because an observer
sees partial attachment/root state. Isolated coverage passes 54 allocation
positions (50 exact preceding states and four complete commits) across off,
underlay, interlay and overlay. It checks all roots, attachment count/lookups,
viewport/source lifetime, direct and nested paths, retry, throwing observers
and the exact no-op. Six neighboring publication rows, image-source, image-
display, window-host, headless-host and Qt window-host tests, and the static
link check pass. The compact loader resolves five relevant libraries from the
evidence directory; the Qt loader resolves seven. The exact full executable
passes 141 direct checks, 21 scenarios and 1,036 pass rows in **308.76
seconds**. Current compact, shared/static libBObol, Obol and OSMesa hashes are
`d532d0c4...`, `d5a6d699...`, `168f5d8e...`, `21872f7f...` and `097aac27...`.

Framebuffer construction and attachment/root insertion are closed separately
by CXX-SOURCE-041 below. Host-policy rollback, window-host destructor-only
cleanup, endpoint presentation cleanup, viewport/cursor setters, feature-store
reparenting and Qt faceplate composition remain separate at this checkpoint.
Host-policy rollback and framebuffer-host destruction are closed by
CXX-SOURCE-042 and CXX-SOURCE-043 below. This bounded closure does not change
the conservative completion ranges below.

## September 10: framebuffer open publication

`CXX-SOURCE-041` closes explicit framebuffer image construction and open
publication. The prior path held its new source and viewport through raw
references, inserted the attachment, and then performed ordinary root
insertion. Allocation could leak the private nodes or retain an attachment
whose viewport was absent from every layer root.

Construction now uses scoped node references. The complete realized viewport
and all three successor roots prepare while private. Only after all fallible
preparation succeeds does the host insert the attachment, commit the roots and
transfer its node ownership. Prepared notification then exposes one complete
state. A notification or presentation-request failure occurs after that commit;
precommit failure releases the private graph and preserves the exact absent
state. Repeated open for the same framebuffer allocates and notifies nothing.

Evidence is retained under
`.build/obol-qualification/20260910-framebuffer-open-publication`. The current
executable forced to load the exact saved predecessor fails an allocation row
which retains the attachment without its root insertion. Isolated coverage
passes 458 allocation positions (266 exact absent states and 192 complete
commits). It checks all three roots, attachment count/lookups, viewport/source
realization, retained extension paths, retry, throwing observers and the exact
no-op. Six neighboring publication rows, image-source, image-display, window-
host, headless-host and Qt window-host tests, and the static link check pass.
The compact loader resolves five relevant libraries from the evidence
directory; the Qt loader resolves seven. The exact full executable passes 142
direct checks, 21 scenarios and 1,037 pass rows in **308.04 seconds**. Current
compact, shared/static libBObol, Obol and OSMesa hashes are `0269e48f...`,
`e6ba539e...`, `9f6dc606...`, `21872f7f...` and `097aac27...`.

Neutral descriptor/open-state rollback is closed separately by CXX-SOURCE-042
below. Controller viewport and platform-host policy, window-host destructor-only
cleanup, endpoint presentation cleanup, viewport/cursor setters, feature-store
reparenting and Qt faceplate composition remain separate at this checkpoint.
Framebuffer-host destruction is closed by CXX-SOURCE-043 below. This bounded
closure does not change the conservative completion ranges below.

## September 10: framebuffer open policy publication

`CXX-SOURCE-042` closes neutral-host descriptor/open-state leakage during
framebuffer open when the controller viewport size is unchanged. The prior path
changed host policy before validating and constructing the stream-backed image,
attachment capacity and root successors. A later failure could therefore leave
the requested descriptor and open flag without the framebuffer attachment.

The base host now prepares its copied descriptor privately and publishes it by
pointer swap after its prerequisites succeed. Framebuffer resources, attachment
capacity and every root replacement prepare before host open. The successful
host call is followed only by no-allocation attachment insertion, root commit
and ownership transfer; callbacks therefore see the complete committed state.
The headless host also establishes its scene root and camera before the base
descriptor commit.

Evidence is retained under
`.build/obol-qualification/20260910-framebuffer-open-policy-publication`. The
current executable forced to load the exact saved CXX-SOURCE-041 library fails
because an allocation row retains the new host policy without the attachment.
Both initially closed and already-open neutral hosts pass 603 allocation
positions (439 exact preceding policies and 164 complete commits). Coverage
checks every descriptor field, open state, unchanged controller viewport,
render request, attachment, all roots and retained paths, retry, throwing
observers and the exact no-op. Six neighboring publication rows, image-source,
image-display, window-host, headless-host and Qt window-host tests, and the
static link check pass. The compact loader resolves five relevant libraries
from the evidence directory; the Qt loader resolves seven. The exact full
executable passes 143 direct checks, 21 scenarios and 1,038 pass rows in
**308.71 seconds**. Current compact, shared/static libBObol, Obol and OSMesa
hashes are `ed914796...`, `2eba79a8...`, `0b00def1...`, `21872f7f...` and
`097aac27...`.

At this checkpoint controller viewport-size publication, headless context-
manager lifecycle, Qt-owned widget policy, endpoint cleanup, viewport/cursor
setters, feature-store reparenting and Qt faceplate composition remain open.
Framebuffer-host destruction is closed separately by CXX-SOURCE-043 below, and
full-window viewport-size publication by CXX-SOURCE-044.
This bounded closure does not change the conservative completion ranges below.

## September 10: framebuffer host destruction

`CXX-SOURCE-043` closes framebuffer cleanup reached by base, headless and Qt
host destruction. The prior destructors reused live close, provider, camera and
scene-root publication paths. Those paths allocate and dispatch callbacks; a
failure therefore escaped an implicit-noexcept destructor and terminated the
process.

The base host now uses a private no-throw cleanup before controller retirement,
and derived hosts run it before destroying platform resources. Completed close
publication still detaches each attachment coherently. If root preparation is
denied, host-owned source and viewport references are released without raw
root mutation. An owned controller retires the remaining graph in the same
destruction. A borrowed controller keeps a complete valid root-owned graph
until its owner can retire it. View-controller destruction directly releases
providers, camera and render-graph ownership without constructing and
publishing replacement live state.

Evidence is retained under
`.build/obol-qualification/20260910-framebuffer-host-destruction`. The current
focused executable forced to load the exact CXX-SOURCE-042 library exits 42 on
the first denied allocation. Current owned and borrowed-controller coverage
passes 18 allocation positions and throwing overlay-root observers. Eight
borrowed-controller denied-preparation cases preserve the complete root-owned
fallback, while every other case fully detaches. Seven neighboring publication
rows, five image/window/Qt tests, the static-link check and all seven model-
conformance tests pass; the catalog validates 34 scenarios and 110 action
mappings. The compact loader resolves five relevant libraries from the
evidence directory and the Qt loader resolves seven. The exact full executable
passes 144 direct checks, 21 scenarios and 1,039 pass rows in **310.168
seconds**. Current compact, shared/static libBObol, Obol and OSMesa hashes are
`02904d05...`, `4fd8c93c...`, `321acd20...`, `21872f7f...` and `097aac27...`.

Exact borrowed-controller root detachment under allocation denial, arbitrary
camera/root observer failure, other endpoint/destructor owners, general
sub-viewport and platform-host policy, cursor setters, feature-store
reparenting and Qt faceplate composition remain open. This bounded closure does
not change the conservative completion ranges below.

## September 10: viewport-size publication

`CXX-SOURCE-044` closes `BObolViewController::setViewportSize()` for ordinary
and automatic-LoD controllers. The preceding implementation changed its stored
region and `SoViewport`, redundantly reattached the render manager's unchanged
scene graph, and let nested LoD/progressive requests wake the endpoint before
the final `viewport-size` request existed. Allocation denial could leave the
three public regions changed or split under the preceding LoD revision.

The resize transaction now prepares its final request first, updates only the
render manager's viewport, and lets LoD synchronization advance its revision
and progressive-work level without an intermediate endpoint wake. A failed
automatic convergence snapshot restores the controller, `SoViewport` and
render-manager regions. Success commits the exact reason and wakes once;
callback failure propagates only after the complete state exists. Existing
private helper entry points remain available through compatibility overloads.

Evidence is retained under
`.build/obol-qualification/20260910-viewport-size-publication`. The current
focused executable forced to load the exact CXX-SOURCE-043 library exits 1 on
the intermediate success callback. Current ordinary and automatic-LoD coverage
passes six allocation outcomes (four exact old states and two complete
commits), retry, exact allocation-free no-op and throwing callbacks. Seven
neighboring publication rows, eight prototype/LoD/image/host/Qt tests, the
static-link check and all seven model-conformance tests pass; the catalog
validates 34 scenarios and 110 action mappings. The compact loader resolves
five relevant libraries from the evidence directory and the Qt loader resolves
seven. The exact full executable passes 145 direct checks, 21 scenarios and
1,040 pass rows in **367.41 seconds**. Current compact, shared/static libBObol,
Obol and OSMesa hashes are `d9ab4dad...`, `fc6c558d...`, `a85333a8...`,
`21872f7f...` and `097aac27...`.

General sub-viewport publication is closed separately by CXX-SOURCE-045 below.
Established-root neutral host composition is closed separately by
CXX-SOURCE-046 below. Arbitrary Coin field observers, first-root and derived-
host policy, camera/cursor setters, feature-store reparenting and Qt faceplate
composition remain open. This bounded closure does not change the conservative
completion ranges below.

## September 10: viewport-region publication

`CXX-SOURCE-045` closes `BObolViewController::setViewportRegion()` and
consolidates both public viewport setters behind the CXX-SOURCE-044 publication
transaction. The preceding general setter still changed visible region copies,
reattached the render manager's scene root, and dispatched nested LoD/
progressive wakes before its final `viewport` request.

The shared transaction handles exact no-op detection, prepared request state,
the three viewport copies, allocation-failure rollback, immediate LoD view
revision, progressive-work state and one final endpoint wake. Full-window size
normalization and the two public request reasons remain distinct. The general
setter preserves the caller's window dimensions, sub-viewport origin and
extent, and pixel density.

Evidence is retained under
`.build/obol-qualification/20260910-viewport-region-publication`. The current
focused executable forced to load the exact CXX-SOURCE-044 library exits 1 on
the intermediate success callback. Current ordinary and automatic-LoD coverage
passes six allocation outcomes (four exact old states and two complete
commits), retry, exact allocation-free no-op and throwing callbacks. The prior
viewport-size row plus seven neighboring publication rows, eight prototype/
LoD/image/host/Qt tests, the static-link check and all seven model-conformance
tests pass; the catalog validates 34 scenarios and 110 action mappings. The
compact loader resolves five relevant libraries from the evidence directory
and the Qt loader resolves seven. The exact full executable passes 146 direct
checks, 21 scenarios and 1,041 pass rows in **341.95 seconds**. Current compact,
shared/static libBObol, Obol and OSMesa hashes are `82c61f31...`, `77dad89d...`,
`7f0abf24...`, `21872f7f...` and `097aac27...`.

Arbitrary Coin field observers, platform-host sizing policy, camera/cursor
setters, feature-store reparenting and Qt faceplate composition remain open.
This bounded closure does not change the conservative completion ranges below.

## September 10: window-open viewport publication

`CXX-SOURCE-046` closes established-root neutral
`BObolWindowHost::open()` composition between host policy and the controller
viewport transaction. The preceding implementation changed the viewport and
dispatched its endpoint callback before swapping the descriptor and setting the
open flag. That callback could therefore observe a mixed
old-host/new-controller policy.

Open now rejects an exact already-open no-op before allocation. After the
already-established group root passes validation, it privately prepares the
normalized descriptor and controller viewport publication, then installs the
descriptor and open flag before finishing the controller transaction.
Allocation failure preserves the complete preceding host/controller state.
Success publishes the three viewport copies, immediate LoD revision,
progressive-work state and exact request reason with one callback; callback
failure propagates after the complete commit.

Evidence is retained under
`.build/obol-qualification/20260910-window-open-viewport-publication`. The
current focused executable forced to load the exact saved CXX-SOURCE-045
library exits 1 when the callback sees the preceding host policy. Current
closed/open and ordinary/automatic-LoD coverage after root establishment
passes 40 allocation outcomes (36 exact preceding policies and four complete
commits), retry, exact allocation-free repeated open and throwing callbacks.
The two viewport rows plus seven neighboring publication checks, eight
prototype/LoD/image/host/Qt
tests, the static-link check and all seven model-conformance tests pass; the
catalog validates 34 scenarios and 110 action mappings. The compact loader
resolves five relevant libraries from the evidence directory and the Qt loader
resolves seven. The exact full executable passes 147 direct checks, 21
scenarios and 1,042 pass rows in **311.13 seconds**. Current compact,
shared/static libBObol, Obol and OSMesa hashes are `0deaa573...`,
`5179915b...`, `3169d391...`, `21872f7f...` and `097aac27...`.

Host-owned first-root construction is closed separately by CXX-SOURCE-047
below. Borrowed-controller root establishment, headless context/camera setup,
Qt widget policy, arbitrary camera/cursor setters, feature-store reparenting
and Qt faceplate composition remain open. This bounded closure does not change
the conservative completion ranges below.

## September 10: window first-open publication

`CXX-SOURCE-047` closes construction and first open for the neutral
host-owned controller. The preceding host exposed a rootless controller and
created its empty group through the live scene-root setter during
`BObolWindowHost::open()`. Its `scene-root` callback could see the new root
with the default descriptor, closed flag and old viewport; the final viewport
request could then coalesce without another wake. The controller's stored and
render-manager regions also began at 1x1 while `SoViewport` retained a
different library default.

Host private construction now keeps the controller under temporary unique
ownership, establishes its empty group before the host is observable, retires
construction-only host work, and releases ownership only after success.
Controller private construction initializes both viewport consumers from the
stored region. Construction failure therefore exposes no object; successful
construction exposes one group root, three matching viewport copies and no
pending work. The first visible open proceeds through the CXX-SOURCE-046
descriptor/viewport transaction and emits one complete callback.

Evidence is retained under
`.build/obol-qualification/20260910-window-first-open-publication`. The
current focused executable forced to load the exact saved CXX-SOURCE-046
library exits 1 on the mixed first-root callback. Current ordinary and
automatic-LoD first opens pass 20 allocation outcomes (18 exact closed states
and two complete commits), retry, allocation-free exact repeat and throwing
callbacks. A 220-position construction sweep rejects 205 incomplete
constructions and produces 15 complete, quiet hosts. Ten neighboring
publication checks, eight prototype/LoD/image/host/Qt tests, the static-link
check and all seven model-conformance tests pass; the catalog validates 34
scenarios and 110 action mappings. The compact loader resolves five relevant
libraries from the evidence directory and the Qt loader resolves seven. The
exact full executable passes 148 direct checks, 21 scenarios and 1,043 pass
rows in **349.88 seconds**. Current compact, shared/static libBObol, Obol and
OSMesa hashes are `5dd6fcdd...`, `0d7988fe...`, `8a27709c...`,
`21872f7f...` and `097aac27...`.

Borrowed-controller attachment/root establishment, headless context/camera
setup, Qt widget policy, arbitrary camera/cursor setters, feature-store
reparenting and Qt faceplate composition remain open. This bounded closure
does not change the conservative completion ranges below.

## September 10: window controller attachment

`CXX-SOURCE-048` closes default-controller construction, borrowed-controller
attachment validation and the first host publication after attachment. The
preceding default `BObolViewController` was rootless. A host could borrow that
incomplete controller and create its root through a separate live setter during
the first open, reproducing the split root/descriptor/viewport publication
which CXX-SOURCE-047 had closed only for host-owned controllers. The display
endpoint C constructor also assembled a default controller and root through a
separate expression which could let a C++ exception cross the C ABI.

Both controller constructors now share one graph-initialization path. The
default form owns an empty `SoBRLSceneGroup`, initializes all three viewport
copies to the same 1x1 region and clears construction-only work before it
returns. Local unique owners retain every graph component until construction
commits. The base host constructs its descriptor first and retains the new
controller in a local `unique_ptr` until both private members are complete, so
failure cannot leak the raw controller.

`BObolWindowHost::attachController()` now reports success, rejects null-root
and non-group-root controllers without changing the existing attachment, and
treats an identical controller as an allocation-free operation. Headless,
display-session, endpoint, GED framebuffer and Qt factory callers propagate
rejection. Default C endpoint creation catches every exception and transfers
ownership only after the endpoint is complete. Controller destruction uses
the LoD service unsubscribe operation directly as its callback-quiescence
barrier, then retires generation and resident-consumer state without invoking
the fallible live service-policy transition.

Evidence is retained under
`.build/obol-qualification/20260910-window-controller-attachment`. The focused
executable forced to load the exact CXX-SOURCE-047 library exits 1 because the
default controller is not hostable. Current construction coverage passes 326
allocation positions: 230 reject without exposing an object and 96 return a
complete quiet controller. Endpoint creation covers 233 positions, containing
232 failures at the C boundary and returning one complete measured endpoint.
Both invalid-root forms reject without allocation or state change. Borrowed
ordinary and automatic-LoD controllers pass 20 first-open positions,
allocation-free identical attachment and throwing callbacks after commit.
Attachment and traced LoD-service teardown each perform zero allocations.

Twelve neighboring publication checks, eleven prototype/LoD/image/RT/host/GED/
Qt/qged tests, the static-link check and all seven model-conformance tests pass;
the catalog validates 34 scenarios and 110 action mappings. The compact loader
resolves five relevant libraries from the evidence directory and the Qt loader
resolves seven. The exact full executable passes 149 direct checks, 21
scenarios and 1,044 pass rows in **360.27 seconds**. Current compact,
shared/static libBObol, Obol and OSMesa hashes are `1888ab0e...`,
`922f2260...`, `60a7349c...`, `21872f7f...` and `097aac27...`.

Headless context/camera setup, Qt widget lifetime policy, arbitrary camera/
cursor setters, feature replacement/reparenting, multi-record controller
migration, forced-failure borrowed framebuffer retirement and Qt faceplate
composition remain open. This bounded closure does not change the conservative
completion ranges below.

## September 10: headless camera invariant

`CXX-SOURCE-049` removes camera creation and scene mutation from the headless
host lifecycle. The preceding headless `open()` searched the modeled scene for
a camera and, when none existed, added a new orthographic camera to that scene
before also installing it in `SoViewport`. This contradicted `SoViewport`'s
contract that the scene exclude cameras, caused the same camera to be reachable
through both the private render root and application content, and exposed a
host-specific scene write during first open.

The default controller now constructs its preserved orthographic camera with
its empty `SoBRLSceneGroup`, installs it only in the private viewport graph and
synchronizes it with the render manager before returning. Construction remains
quiet. Headless open and render only validate the existing group root and
active camera; they never search, adopt or parent scene cameras. Camera-less
borrowed controllers reject before descriptor allocation or provider binding,
and the built-in headless factory validates before attachment and contains
binding exceptions.

Evidence is retained under
`.build/obol-qualification/20260910-headless-camera-invariant`. The focused
executable forced to load the exact CXX-SOURCE-048 library exits 1 because its
default controller has no active camera. Current coverage proves one camera is
shared by the controller, viewport and render manager, is present in the
private viewport root and absent from modeled content. Successful open retains
camera identity and scene children. Camera-less empty-scene and misplaced-
scene-camera controllers reject without allocation or mutation. The
strengthened default-construction sweep passes 380 positions (259 rejected,
121 complete); endpoint construction passes 262 positions (261 contained
rejections, one complete).

Thirteen neighboring publication checks, eleven prototype/LoD/image/RT/host/
GED/Qt/qged tests, the static-link check and all seven model-conformance tests
pass; the catalog validates 34 scenarios and 110 action mappings. The compact
loader resolves five relevant libraries from the evidence directory and the
Qt loader resolves seven. The exact full executable passes 150 direct checks,
21 scenarios and 1,045 pass rows in **332.96 seconds**. Current compact,
shared/static libBObol, Obol and OSMesa hashes are `da0b8783...`,
`9e40d9db...`, `e6ec6b5c...`, `21872f7f...` and `097aac27...`.

This closes headless camera creation/adoption. Headless context-provider
publication, general live camera replacement, Qt widget lifetime policy and
the other listed lifecycle and resource rows remain open. This bounded closure
does not change the conservative completion ranges below.

## September 10: headless context-provider publication

`CXX-SOURCE-050` composes provider binding with the headless host lifecycle.
The preceding `open()` invalidated renderer evidence and requested a frame
before assigning the new provider or publishing the host descriptor, viewport
and open flag. The callback therefore observed a closed host and the preceding
provider. Close had the symmetric defect, and live provider replacement first
changed the host pointer before the controller operation could fail.

Controller provider changes now have explicit prepare, bind, renderer-
invalidation commit and notification phases. The window host similarly
separates prepared descriptor/viewport state from commit and notification.
Headless open combines those phases, suppresses the redundant renderer wake
when the viewport request subsumes it, and issues one callback after provider,
renderer evidence, descriptor, viewport and open state agree. Replacement and
close use the same provider transaction. Initial factory binding waits for the
open transaction; rebinding an open instance goes through the host publisher.

Native Obol now registers a replacement provider before changing
`SoGLRenderAction`, so first registration preserves the preceding action when
map growth fails and replacement updates the existing registry entry without
allocation. Cache-context replacement similarly registers the new key before
retiring the old key. The retained predecessor stacks independently reject the
premature host callback and the allocating registered-provider replacement.

The focused selector proves allocation-free registered-provider replacement,
two cache-key allocation outcomes (one preserved and one committed), 14
ordinary/automatic first-open outcomes (12 preserved and two committed), and
16 live replacement/close outcomes (eight preserved and eight committed).
Exact no-ops allocate and notify zero times; callback exceptions occur after
the complete state. Fourteen neighboring publication selectors, twelve
prototype/render/LoD/image/RT/host/GED/Qt/qged tests, three public/install
contract tests, the static consumer, 27 native Obol render-action tests and all
seven model-conformance tests pass. The exact full executable passes 151 direct
checks, 21 scenarios and 1,046 rows in **356.62 seconds**.

Evidence is retained under
`.build/obol-qualification/20260910-headless-context-publication`. Current
compact/shared/static libBObol, Obol and OSMesa hashes are `34a85acc...`,
`ad9d0a4f...`, `b4644c79...`, `4e36cf66...` and `097aac27...`.

This closes headless provider publication. General live camera replacement,
Qt widget lifetime policy and the other lifecycle, resource, graphical and
platform rows remain open. The conservative completion ranges below are
unchanged.

## September 10: feature replacement and reparenting publication

`CXX-SOURCE-051` composes an existing feature's record, retained node,
attachment root and overlay order. The preceding implementation changed its
live record before node construction, detached the old node before adding its
successor, and reordered overlays separately. Root, path and result callbacks
could therefore observe incomplete state; allocation failure could preserve
the old graph with new metadata; and replacing one of two records which shared
a custom node could detach the other record's presentation.

The store now copies and configures an existing record off-live, retains the
candidate node, computes each affected root from the effective successor
record set, and prepares every child-list replacement before mutation. It
commits root edits and the record pointer without allocation, releases the old
node after the successor is installed, and only then notifies graph observers
and requests presentation. Shared custom nodes retain one root occurrence
until the final attachment moves or is removed. All existing-record feature
publishers which rebuild retained nodes use this common path; ordinary overlay
reorder and transactional replacement use the same ordering rule.

The exact CXX-SOURCE-050 shared library rejects the focused selector with
`feature node replacement or reparent exposed partial state`. Current typed,
custom-alias and screen-reparent cases pass 108 forced-allocation outcomes: 98
retain the complete predecessor and ten reach the complete successor. They
also cover retry, an allocation-free exact reparent no-op, immediate root/path
and result observers, and exceptions from graph and result callbacks. Nine
neighboring publication selectors, ten feature/GED/Qt integration tests,
three public/install contract checks, static link and all seven conformance
checks pass. The exact full executable passes 152 direct checks, 21 scenarios
and 1,050 rows in **356.73 seconds**. Current compact, shared/static libBObol,
Obol and OSMesa hashes are `c838a8d8...`, `8f94eaba...`, `60b37629...`,
`4e36cf66...` and `097aac27...`. Exact baseline/current binaries,
dependencies, source snapshots, loader traces, logs and hashes are retained under
`.build/obol-qualification/20260910-feature-node-publication`.

This closes existing feature replacement/reparenting and
`setOverlayInfo()` record/node/root/order composition. Multi-record controller
migration, initial feature insertion, in-place metadata and mesh-field edits,
Qt/lifecycle/resource gates, graphical qualification and platform rows remain
open. This bounded closure does not change the conservative completion ranges
below.

## September 12: feature controller migration publication

`CXX-SOURCE-052` composes all feature records and roots when a standalone
feature store changes controllers. The preceding `setController()` detached
and attached records one at a time, wrote each attachment pointer during that
loop, changed the public controller afterward, and then reordered overlays.
Immediate root and path observers could see a prefix paired with the old
controller. Allocation or callback failure could strand records across old and
new roots without a presentation-revision witness.

The store now projects all retained nodes into their desired roots, preserving
unrelated children, shared custom-node occurrences and deterministic overlay
order. It prepares every affected old and new child-list replacement before
mutation. Root orders, controller identity, per-record attachment roots and one
presentation revision commit together. All root/path notifications are then
attempted, followed by presentation requests to both controllers, while the
first callback exception is retained. Rebinding the current controller is an
allocation-free, notification-free no-op. Single-record replacement, ordinary
overlay reorder and whole-store migration share the same root-order composer.

The exact CXX-SOURCE-051 library rejects the focused selector with `feature
controller migration exposed partial roots or ownership`. Current coverage
uses five model/screen records, a two-record custom-node alias, four unrelated
root extensions and every retained feature path. Replacement, detachment and
adoption pass 170 forced-allocation outcomes: 159 retain the complete
predecessor and 11 reach the complete successor. Retry, exact no-op, all four
old/new roots, both controller frame callbacks and throwing graph/frame
callbacks also pass. Ten neighboring publication selectors, ten
feature/GED/Qt integration tests, three public/install contract checks, static
link and all seven conformance checks pass. The exact full executable passes
153 direct checks, 21 scenarios and 1,054 rows in **360.67 seconds**. Current
compact, shared/static libBObol, Obol and OSMesa hashes are `d933b9a1...`,
`d6c1a1b5...`, `140a4b1f...`, `4e36cf66...` and `097aac27...`. Exact
binaries, selected dependencies, source snapshots, loader traces, logs and hashes are retained
under `.build/obol-qualification/20260912-feature-controller-migration`.

This closes multi-record feature-controller migration. Initial feature
insertion, in-place metadata and mesh-field edits, Qt faceplate/lifetime work,
other lifecycle/resource gates, graphical qualification and platform rows
remain open. This bounded closure does not change the conservative completion
ranges below.

## September 12: Qt host rebinding validation

`CXX-SOURCE-053` closes the destructive invalid-input path in
`QgObolWindowHost::setCanvas()`. The old implementation destroyed an owned
canvas and published the replacement pointer before the base host validated
the replacement controller. A controller with a missing or non-group scene
root was then rejected after the adapter had already split its public canvas
and controller identities.

The Qt host now uses the same hostability predicate for canvas replacement and
direct controller binding before changing live state. `setCanvas()` returns
its acceptance result, invalid input preserves the complete borrowed or owned
predecessor, and a same-canvas request is an exact no-op. The factory-bind
canvas-preservation flag also has a scope-bound lifetime so an exception cannot
silently change later close ownership.

The exact CXX-SOURCE-052 libqtcad fails the focused executable at `rejected
borrowed-canvas replacement preserves host identity`. Current coverage checks
direct controller rejection, borrowed and owned canvas rejection, owned-canvas
lifetime, exact no-op and accepted replacement. The focused executable,
twelve Qt functional tests, Qt plugin reuse, qged framebuffer-host integration
and staged public header copy pass. The System-GL-only host row records its
defined environment skip on this seat. Exact focused executable,
predecessor/current libqtcad, libBObol, Obol and OSMesa hashes are
`d9d6b275...`, `ef3a5c4e...`, `6709cefc...`, `d6c1a1b5...`,
`4e36cf66...` and `097aac27...`. Binaries, source snapshots, loader traces,
logs and hashes are retained under
`.build/obol-qualification/20260912-qtcad-host-rebinding`.

This closes rejected Qt canvas/controller input and its ownership consequence.
Forced failure during a valid live rebind, framebuffer retirement, remaining
faceplate synchronization and Qt widget policy stay separate. This bounded
closure does not change the conservative completion ranges below.

## September 12: Qt fallback faceplate node publication

`CXX-SOURCE-054` closes direct live-node mutation in the no-GED Qt faceplate
fallback. The previous grid path assigned center, spacing, divisions, view
inputs and derived geometry to the attached node one field at a time. Axes and
ADC followed the same live-input pattern. Immediate sensors could therefore
observe a mixture of predecessor and successor inputs before the rebuilt child
graph was ready.

The fallback now builds each grid, axes or ADC successor off graph. A shared
publication helper prepares one complete root child list, commits it once, and
then attempts graph notification and the frame request while preserving the
first exception. It covers add, replace and remove; a typed lookup replaces
the three former child scans. The view/frame-hash no-op continues to preserve
existing nodes and geometry.

The exact CXX-SOURCE-053 libqtcad fails the focused executable at `Qt
faceplate grid exposed partially updated input fields`. Current focused
coverage observes only the complete predecessor or successor, verifies the
new adaptive grid, preserves unchanged-refresh identity, and retains existing
axes/ADC/removal checks. Twelve Qt functional tests, Qt plugin reuse, qged
framebuffer-host integration and the graphical GED faceplate row pass. The
System-GL-only host row records its defined environment skip. Exact focused
executable, predecessor/current libqtcad, libBObol, Obol and OSMesa hashes are
`0d00e4c0...`, `6709cefc...`, `9ca01e87...`, `d6c1a1b5...`,
`4e36cf66...` and `097aac27...`. Exact artifacts are retained under
`.build/obol-qualification/20260912-qtcad-faceplate-node-publication`.

This closes per-node publication in the direct Qt fallback. Whole-faceplate
GED composition, direct mutation through `bobol_grid_configure_from_view`,
axes/ADC public input setters and Qt widget policy remain separate. This
bounded closure does not change the conservative completion ranges below.

## September 12: initial feature publication

`CXX-SOURCE-055` closes first-record feature publication. The former path
inserted the record and name indexes before payload copying, retained-node
construction and root attachment. Failures could consume an identity or leave
a partial live record. It also notified the root before its presentation
revision and controller request existed, and a request-preparation failure
after commit could leave an identical retry with no repaint to schedule.

The feature store now constructs the candidate and advances a local next-ID
copy off-store. It preallocates both map entries, prepares the complete child
replacement and prepares the controller request before mutation. Commit then
publishes record identity, both indexes, node ownership, root order,
presentation revision and the standing request without allocation. Graph and
frame notifications are both attempted after that state is visible, retaining
the first exception and the exact prepared controller target. Existing-record
node publication shares the request preparation and ordering.

The exact CXX-SOURCE-054 libBObol fails the current contract at `feature
insertion failure left a partial record, graph or request`. The current focused
selector passes 138 positions: 48 preserve the exact empty predecessor and its
next identity; 90 expose only the complete successor. Graph, frame and
command-result callback failures retain the committed successor. Ten adjacent
feature checks, ten production integration tests, three public/install tests
and the static consumer pass. The complete compact executable passes 154
direct checks, 21 scenarios and 1,055 PASS rows. Exact executable,
predecessor/current shared library, current static library, Obol and OSMesa
hashes are `c933c619...`, `d6c1a1b5...`, `7dbacbea...`, `dc42e322...`,
`4e36cf66...` and `097aac27...`. Exact artifacts are retained under
`.build/obol-qualification/20260912-feature-initial-publication`.

This bounded closure does not qualify whole-faceplate composition, direct
in-place metadata/mesh fields, removal/clear, controller-migration request
preparation or the wider S4--S6 matrix. The completion ranges below remain
conservative.

## September 12: whole-faceplate GED composition publication

`CXX-SOURCE-056` closes multi-record faceplate composition in the GED adapter.
The former refresh published or removed center dot, interactive rectangle,
grid, ADC, parameters, LoD progress, scale and axes records one at a time.
Immediate graph or frame observers could therefore see a mixture of old and
new faceplate records before the logical refresh finished.

The feature store now accepts one mixed replace/remove publication for line
sets, HUD labels and retained custom nodes. It validates the complete input,
stages all records and nodes, preallocates new index entries, prepares every
affected final root order and prepares one controller presentation request
before mutation. One allocation-free commit installs the complete successor,
advances the presentation revision once and makes the standing request visible
before callbacks. Exact no-ops preserve identities and issue no callbacks;
invalid or duplicate entries reject without changing the predecessor. Graph,
frame and command-result callbacks are all attempted after commit while the
first exception is retained.

The GED faceplate adapter now accumulates every enabled replacement and
disabled removal and submits them through that primitive. Its exact
CXX-SOURCE-055 predecessor exposes the defect in a clean saved-library probe:
18 immediate graph callbacks include a partial faceplate state before the final
complete state. The current real-GED selector observes only complete insertion
and replacement/removal refreshes, including overlay metadata and the single
feature-store frame request. The direct failure sweep covers 373 allocation
positions: 231 preserve the exact predecessor and 142 expose the complete
successor. Neighboring feature checks, eleven production integration tests,
three public/install tests, the static consumer and seven conformance tests
pass. The complete compact executable passes 155 direct checks, 21 scenarios
and 1,056 PASS rows. Exact direct/GED executables, predecessor/current
libBObol and predecessor/current libged hashes are `42bfc460...`,
`f1941c87...`, `7dbacbea...`, `01c58527...`, `5abb7712...` and
`d70e5358...`. Exact artifacts are retained under
`.build/obol-qualification/20260912-ged-faceplate-composition`.

Framebuffer composition and presentation still follow the feature commit, so
they are not part of this store transaction. The public
`ged_view_feature_batch` adapter also continues to replay general operations
sequentially. Direct in-place metadata/mesh fields, controller-migration
request preparation and the wider S4--S6 matrix remain separate. This bounded
closure does not change the conservative completion ranges below.

## September 12: feature removal and clear publication

`CXX-SOURCE-057` closes the public feature removal and clear paths. Handle
removal formerly erased the name index, invoked its result callback, detached
the node, erased the record and only then requested a frame. A callback could
observe or strand that partial state. Prefix and scope removal repeated the
sequence for each matching handle, exposing partial sets and advancing the
presentation once per record. Clear likewise notified records before retiring
them and reset command-owner generations only after those callbacks.

Handle, name, owned-name, prefix, scope and public clear now use the prepared
feature publication transaction. Each wrapper builds the complete removal set
before entering the store. The transaction stages result notifications,
prepares all final root orders and one reason-preserving controller request,
then commits indexes, records, roots, retained ownership, presentation revision
and, for clear, command-owner generations without allocation. Graph, frame and
all result callbacks run against the complete retired state; the first callback
exception is rethrown only after the remaining notifications are attempted.
Unmatched removals and an empty clear remain presentation no-ops.

A focused executable linked to the exact CXX-SOURCE-056 libBObol rejects the
old path at `feature removal measurement exposed a partial publication`.
Current handle, prefix, scope and clear sweeps cover 439 forced-allocation
positions: 228 preserve the exact predecessor and 211 expose only the complete
successor. They also cover throwing graph, frame and result callbacks, filtered
survivors, clear-time generation reset, generation-only clear and unmatched
no-ops. Four adjacent feature selectors, thirteen production integration tests,
three public/install tests, the static consumer and seven conformance tests
pass. The complete compact executable passes 156 direct checks, 21 scenarios
and 1,061 PASS rows. Exact focused executable and predecessor/current shared
libBObol hashes are `b6e41cf6...`, `01c58527...` and `edb42f3e...`.
Exact artifacts are retained under
`.build/obol-qualification/20260912-feature-removal-publication`.

This boundary covers the public store APIs; destructor-time record teardown is
still a separate lifecycle concern. Direct in-place metadata and mesh fields,
generic GED batch composition and the wider S4--S6 matrix also remain
separate. The completion ranges below stay conservative.

## September 12: feature-controller migration request publication

`CXX-SOURCE-058` completes the render-obligation side of feature-controller
migration. CXX-SOURCE-052 already prepared and committed all affected old/new
root orders, controller ownership and record attachment roots together. It
still notified those graph changes before calling
`requestPresentationRender` on the previous and next controllers. Those calls
prepared requests after commit, so an allocation failure could leave migrated
content without a standing repaint, and graph observers saw the complete graph
before either render obligation existed.

The controller-independent request preparation used by feature publication is
now shared with migration. Previous-controller detach and next-controller
attach requests, including their LoD transition scopes, prepare after all root
replacements and before any live change. The allocation-free commit installs
the roots, controller and attachment ownership, presentation revision and both
standing requests. Root observers then see the complete obligations; request
notifications for both controllers are attempted while retaining the first
callback exception. A same-controller migration remains an allocation-free
no-op.

The enhanced migration selector rejects the exact CXX-SOURCE-057 libBObol at
`feature controller migration exposed partial roots or ownership`. Current
replace, detach and attach sweeps pass all 170 allocation positions: 167 retain
the exact predecessor and the three allocation-free commits expose complete
successors. Immediate root and frame observers verify standing old/new
requests and exact detach/attach reasons. Throwing graph and frame callbacks do
not stop later notifications. Five neighboring feature selectors, thirteen
production integration tests, three public/install tests, the static consumer
and seven conformance tests pass. The complete compact executable passes 156
direct checks, 21 scenarios and 1,061 PASS rows. Exact focused executable and
predecessor/current shared libBObol hashes are `c5aace3a...`, `edb42f3e...`
and `c3d93019...`. Exact artifacts are retained under
`.build/obol-qualification/20260912-feature-controller-request-publication`.

No formal baseline changes. At this checkpoint, destructor-time feature
teardown, direct in-place metadata and mesh fields, generic GED batch
composition and the wider S4--S6 matrix remained separate. This closure did
not change the conservative stage ranges.

## September 12: direct feature record and mesh-field publication

`CXX-SOURCE-059` closes five direct public feature-store edits which still
changed published state in pieces. Indexed-face point edits patched the live
record before replacing point and normal arrays; mesh selection/highlight edits
changed record and `SoMFInt32` state separately. Their source revision and
controller request followed those mutations. Metadata and primitive metadata
also assigned into the live record before their result callback payload had
been prepared.

Obol now prepares a same-type `SoMField` successor off graph. Its allocation-
free commit swaps owned storage and dimensions without notification, and its
one-shot notification restores the target's previous notification state. The
feature store composes these replacements with a copied record, mesh source
revision, presentation revision, standing request and prebuilt result. All
fallible preparation precedes the first live change; field, node, frame and
result observers then see the complete successor, with every notification
attempted if an earlier one throws. The retained mesh node and field objects do
not change identity. Metadata-only edits use the candidate-record/result half
of that boundary and do not request a frame.

The exact CXX-SOURCE-058 libBObol fails the enhanced indexed-point case at
`direct feature mutation exposed partial record, field or request`. The current
five-operation sweep passes **529 allocation positions**: 144 preserve the
exact predecessor and 385 expose the complete successor. It includes point and
selection array growth/shrinkage, feature/field/node identity, exact request
reasons, immediate field/node/frame/result observation and graph, frame and
result callback failures. Six neighboring feature selectors pass, including
the adjusted preparation counts for feature batch publication. Thirteen
feature/GED/Qt integration tests, three public/install tests, the static
consumer and seven conformance tests pass.

Native Obol's focused field replacement test passes, as do **1,580 unit tests
with one existing profiler skip** and all **63 integration tests**. The complete
compact executable passes **157 direct checks, 21 scenarios and 1,067 PASS
rows in 357.81 seconds**. OSMesa was restored after CMake staging and remains
`097aac27...`. Exact binaries, dependency-local replays, source snapshots,
logs, hashes and the verifier are retained under
`.build/obol-qualification/20260912-feature-direct-publication`.

No formal baseline changes. Externally owned custom-node mutation,
destructor-time feature teardown, generic GED batch composition and the wider
S4--S6 matrix remain separate. The bounded closure improves the S1/S3 inventory
but is too small relative to the remaining gate uncertainty to move the
conservative ranges below.

## September 12: atomic GED edit-transaction publication

`CXX-SOURCE-060` closes the retained edit transaction which MGED and Qt use to
publish a live primitive preview. The old implementation called overlay ensure,
overlay event, color, geometry/semantic replacement and feature `touch()` in
sequence. A first preview update exposed three retained graph generations; an
immediate observer saw an incomplete preview during one of them. The direct
touch writer also had no presentation revision or frame obligation of its own.

Primitive transformation, plotting and VLIST conversion now finish before any
retained mutation. The GED adapter builds one complete edit-preview publication
containing geometry, commands, color, semantic identity and intent, source and
input revisions, typed overlay metadata and owner. The feature store prepares
its retained node, root order, presentation revision and standing render
request, then commits the publication once. Identical live events publish
nothing; terminal events retire the preview without doing obsolete geometry
work. The feature-store and GED touch APIs were removed.

The exact CXX-SOURCE-059 libged/libBObol pair fails the dynamic causal probe
with three graph callbacks, one partial observation and revision 0 to 3. Current
records one complete graph callback, one frame request and revision 0 to 1.
Focused GED coverage exercises point-array insertion and replacement,
transformed internal-primitive insertion, exact no-op behavior and terminal
retirement. The mixed feature batch passes **570 allocation positions** and
graph, frame and result callback failure with its edit-preview entry.

No formal baseline changes. Externally owned custom-node mutation and its
revision/notification contract, destructor-time feature teardown, the general
public `ged_view_feature_batch` composition and the wider S4--S6 matrix remain
separate. This bounded S1/S3 repair does not move the conservative stage ranges.

## September 12: atomic initial custom-node overlay publication

`CXX-SOURCE-061` closes the initial publication sequence used by qged's ARB,
ellipsoid, sketch and BoT edit presenters and by the libged navigation gizmo.
Those callers previously attached a custom node before publishing its overlay
classification. ARB, ellipsoid and sketch also populated manipulator geometry
and selection state after attachment. An immediate observer could therefore
see an unclassified or default custom node even though the eventual state was
correct.

The feature store now accepts a complete overlay value when publishing a
custom node and stages it in the same candidate record as node ownership,
ordered-root placement, presentation identity and the frame request. The
production callers construct this value first; ARB, ellipsoid and sketch also
finish initial node configuration while detached. The redundant initial frame
requests are gone. The compatible original overload remains available and
preserves an existing overlay on replacement.

Against the exact CXX-SOURCE-060 libBObol, the causal probe records two graph
callbacks, one partial observation, one coalesced frame callback and a
presentation-revision delta of two. Current records one complete graph
callback, one frame callback and a delta of one. The permanent focused test
also proves complete graph, frame and result observation for the typed X-ray
edit-handle overlay across **75 forced-allocation positions**: 58 preserve the
empty predecessor and 17 expose the complete successor. The complete compact
executable passes **158 direct checks, 21 scenarios and 1,068 PASS rows in
316.99 seconds**. Fourteen integration/API tests, all seven conformance tests,
public-header self-containment and the static consumer pass. Exact artifacts
and a verifier are retained under
`.build/obol-qualification/20260912-custom-node-overlay-publication`.

No formal baseline changes. Mutable updates to an already published external
node remain separate because the store cannot prepare changes made directly by
the owner. Destructor-time feature teardown, general GED feature-batch
composition and wider S4--S6 qualification also remain open. This bounded
S1/S3 closure does not move the conservative ranges below.

## September 12: atomic edit-manipulator setter publication

`CXX-SOURCE-062` closes the node-local half of later custom-node mutation for
the two retained edit manipulators. Their setters previously wrote private
axes/topology or public selection and visibility fields before rebuilding
children. An allocation failure during that rebuild could preserve old
children with new input state. Field and node observers could see the same
partial prefix; indexed topology could be cleared before even copying its
successor.

Axis and indexed geometry now build into retained off-node owners. The child-
list replacement prepares before live mutation, then commits with the private
candidate and quiet scalar fields. Notification policy is restored before
path, node and field callbacks run, and all callbacks are attempted before the
first exception is rethrown. Indexed point/edge/face queries and rendering use
one candidate representation. Unchanged scalar setters remain allocation- and
callback-free.

The new permanent selector exercises seven setter paths, retained paths,
immediate observers, retry, exact no-ops and callback exceptions. It passes
**767 allocation positions**: 759 preserve the complete predecessor and eight
publish the complete successor. Forced loading of the exact CXX-SOURCE-061
library fails with a partial-state diagnostic. The adjacent rebuild sweep
passes 449 positions. The complete compact executable passes **159 direct
checks, 21 scenarios and 1,069 PASS rows in 315.45 seconds**. Fourteen integration/API tests,
seven conformance tests, public-header self-containment and the static consumer
pass.

No public API or formal baseline changed. Feature-store revision/frame
composition for live external-node updates, the BoT shared-mesh hot path and
navigation-gizmo mutation remain open. This bounded S1/S3 closure does not move
the conservative ranges below.

## September 12: atomic qged manipulator replacement publication

`CXX-SOURCE-063` closes the store-facing half of later updates for qged's small
ARB, ellipsoid and sketch manipulators. Those adapters previously invoked
several now-atomic node setters and then requested a frame. The custom feature's
record revision did not change, and one input or session update could expose
several node generations before its frame obligation existed.

Each adapter now constructs the complete successor off graph and replaces the
same custom-node feature through one store transaction. Geometry, session
revision, selection domain/index, hover/active state, typed X-ray overlay,
ordered attachment and render request therefore become visible together. ARB
and sketch fold session selection into the geometry successor. Plugin state
keeps the stable feature handle rather than a borrowed node pointer and resolves
the current record at use, including after commit callbacks and reentrant graph
or frame observation.

The permanent primitive-edit replay asserts that all three feature records
advance during live interaction. Exact preceding plugin binaries fail the first
ellipsoid assertion at revision one. Current passes the complete single-view
and quad-view replays; five adjacent store, manipulator, API, GED and Qt tests
pass. Exact artifacts and a verifier are retained under
`.build/obol-qualification/20260912-qged-manipulator-publication`.

No public libBObol API or formal baseline changed. BoT continues to share a
large retained mesh and is not a candidate for whole-node replacement. Its hot
path and the navigation gizmo retain separate qualification. This bounded
S1/S3 closure does not move the conservative ranges below.

## September 12: coordinated BoT shared-mesh publication

`CXX-SOURCE-064` closes the retained-surface half of interactive BoT point
moves. The adapter formerly changed the shared mesh and manually cleared the
published wrapper caches without replacing their feature records. The first
repair shape also exposed a preparation-failure prefix by changing the shared
points before allocating the successor records.

The feature store now has a coordinated multi-store operation. It prepares all
records, custom nodes, root replacements and frame requests, runs one declared
allocation-free shared-state commit, commits every participating store and
then delivers callbacks. Custom mesh publication records snapshot selected and
highlighted primitive state and reject records which disagree with the node.
QBot prepares one small wrapper per view and retains the single heavy mesh and
topology. Its commit changes only one to three points, clears stale normals and
advances the shared source identity. All per-view wrappers then change identity
and rebuild their private render cache over that same heavy-geometry pointer.

The new two-store sweep passes **128 allocation positions**: 79 preserve the
complete two-store/shared-geometry predecessor and 49 expose the complete
successor. Immediate graph and frame observers see both stores committed. The
neighboring mixed feature batch remains green across **602 positions**. The
complete compact executable passes **160 direct checks, 21 scenarios and 1,070
PASS rows in 315.28 seconds**. The exact preceding BoT plugin fails the focused
GUI successor check; current passes the focused case and full single-view and
quad-view primitive-edit replays. Four adjacent store/API/GED/Qt tests, public
header self-containment and the static consumer pass. Evidence is retained at
`.build/obol-qualification/20260912-qged-bot-shared-mesh-publication`.

The GED point handle still publishes through its separate batch, and a general
cross-controller composition remains open. The GUI fixture proves that the
heavy mesh is retained; it does not qualify large-BoT timing or memory. The
navigation gizmo remains the next live custom-node boundary. This bounded
S1/S3 closure does not move the conservative ranges below.

## September 12: retained navigation-gizmo successor publication

`CXX-SOURCE-065` closes semantic live mutation in the libged navigation gizmo.
The plugin formerly cached a published node pointer, changed its style,
visibility and interaction fields in place, rebuilt its HUD directly and
requested frames separately. Its live camera sensor and render traversal could
also replace orientation-dependent children without advancing the feature
record.

The plugin now owns one value snapshot and resolves the current node by stable
feature identity. Style, visibility, hover/active state and camera orientation
are used to build a complete detached fixed-camera node. One coordinated
feature publication commits the snapshot, record, node, typed faceplate
overlay, ordered root and frame obligation before callbacks. Equal
presentations are exact no-ops. Click and drag updates suppress the camera
sensor while updating the BRL-CAD view, synchronize the display endpoint even
when the GED host callback is absent, and publish the resulting interaction and
camera state once.

The new fixed-orientation node API replaces existing orientation-dependent HUD
children through the failure-atomic Coin child transaction. The expanded
focused selector passes **1,941 allocation positions**: 1,936 complete
predecessors and five complete successors. The exact preceding plugin fails all
thirteen new retained-generation assertions; current passes. The exact full
compact executable passes **160 direct checks, 21 scenarios and 1,070 PASS rows
in 335.65 seconds**. Eight adjacent tests, header self-containment and the
static consumer pass. The full primitive-edit replay and a focused gizmo
enable/orient/disable replay pass with OSMesa in single and quad layouts.
System GL was attempted but the configured X display was unavailable, so it is
not claimed by this boundary. Exact artifacts and a verifier are retained at
`.build/obol-qualification/20260912-navigation-gizmo-publication`.

Viewport-anchor movement remains renderer-derived layout state. Input-layer
installation/removal is qualified only under the endpoint's owner-thread
ordering, and idle replacement of the controller camera object remains in the
general live-camera boundary. The separate GED point-handle batch,
cross-controller composition and wider S4--S6 matrix remain open. This bounded
S1/S3 closure does not materially change the conservative ranges below.

## September 12: live camera replacement publication

`CXX-SOURCE-066` closes direct replacement, removal and addition of the live
controller camera. The preceding native viewport removed and notified the old
camera before inserting its successor, and the controller then changed its own
pointer, render-manager pointer, LoD view state and frame request in separate
steps. Root and path observers could see mixed camera owners, and failure could
interrupt the transition.

Native Obol now prepares the complete root child order while retaining both
camera generations. Its commit changes the root and viewport camera without
callbacks; notification follows the complete state. The controller prepares
that transaction and its capacity-level render request, computes the candidate
LoD signature without an endpoint wake, and commits all three camera pointers,
the LoD revision and the standing request before native and frame observers.
It attempts every notification and rethrows the first observer exception after
the successor is complete. Identical camera pointers are allocation-free exact
no-ops.

The navigation-gizmo plugin observes the completed viewport-root publication
to retarget its camera-orientation field sensor. Direct replacement publishes
one new fixed-orientation gizmo; later mutation of the retired camera publishes
nothing, while mutation of the replacement publishes one successor. The exact
preceding plugin fails those two new integration assertions.

The fresh focused allocation sweep covers **181 positions**: 174 complete
predecessors and seven complete successors across three operations and both
ordinary and automatic LoD. The warmed full matrix covers 180 positions, 174
predecessors and six successors. The exact preceding controller/native pair
fails the new complete-state assertion. Current Obol passes all 13 viewport
tests; the affected BObol, GED and Qt tests, public umbrella, static consumer
and exported-symbol checks pass. The exact full compact executable passes
**161 direct checks, 21 scenarios and 1,071 PASS rows in 357.19 seconds**.
Evidence is retained at
`.build/obol-qualification/20260912-live-camera-replacement`.

The enclosing `syncCameraFromViewContext()` calculation still composes camera
fields with tracked lighting, clipping and the cutting-plane affordance. It
remains a separate input-field boundary. Qt widget lifetime, cross-controller
feature composition and the wider S4--S6 matrix also remain open. This bounded
S1/S3 closure does not materially change the conservative ranges below.

## September 12: complete view-camera input publication

`CXX-SOURCE-067` closes `syncCameraFromViewContext()` as the enclosing
controller input boundary. Previously, projection replacement exposed a
default camera before its fields were set; viewport adoption, three camera
lights, three clip nodes, the section aid, camera fields, LoD state and the
frame request then changed in separate observable steps. Failure could retain
any prefix. The section aid was built against the preceding camera, equal
inputs could rebuild it, and clip-only input was not reported or rendered.

The controller now computes one complete candidate and prepares every fallible
scalar, graph, LoD and render-request resource before an allocation-free
commit. All public camera and viewport copies, lights, clips, section aid,
private depth/affordance caches, LoD revision and standing request are visible
before root, field or endpoint callbacks. Notification state is restored
before callback dispatch, all callbacks are attempted, and observer failure is
reported after the committed successor remains visible. Same-projection
navigation retains camera identity; projection changes replace it; exact
input allocates and notifies nothing.

The dynamic navigation-gizmo test now drives the complete input boundary.
It proves one successor for a same-projection orientation change, no successor
for an exact repeat, one retargeted successor for projection replacement, and
no response from the retired camera. That last projection case exposed and
repaired a plugin equality shortcut which suppressed the successor whenever
the new camera happened to have the same orientation.

The fresh allocation sweep covers **1,578 positions** across five operations
and ordinary/automatic LoD: 1,536 complete predecessors and 42 complete
successors. The saved CXX-SOURCE-066 library rejects the new complete-state
test, and its plugin rejects the same-orientation camera-generation test.
The focused selector and dynamic plugin pass under ASan. The direct camera,
viewport, cutting-affordance, scene-light, render-environment and source-field
publication checks pass, as do the affected LoD, host, public API, GED and Qt
tests. The exact full compact executable passes 162 direct checks, 21
scenarios and 1,072 PASS rows in 363.92 seconds. Evidence is retained at
`.build/obol-qualification/20260912-view-camera-sync`.

Public cutting-plane setters, Qt widget lifetime, cross-controller feature
composition and the wider S4--S6 matrix remain open. This bounded S1/S3 closure
does not materially change the conservative ranges below.

## September 12: complete cutting-plane input publication

`CXX-SOURCE-068` closes the public cutting-plane setters as one controller
input boundary. Previously, private plane intent, live clip fields, visible
section-aid geometry and the render request changed in sequence. Immediate
observers could see a mixed generation, and allocation or callback failure
could stop the remaining effects. Equal planes still rebuilt HUD geometry;
non-finite input could enter the controller.

The two setters now share one prepared transaction. They construct detached
clip scalar state, an optional affordance successor and one capacity-level
render request before committing the clip fields, affordance, private values
and standing request without allocation. Notification state is restored before
callbacks, all committed notifications are attempted, and the first callback
exception is propagated after the complete state remains visible. Exact input
does no work, while zero-normal and non-finite planes are rejected before
allocation.

The focused matrix covers plane changes while disabled and enabled plus enable
and disable operations, each under ordinary and automatic LoD. Its 711 forced-
allocation positions retain 636 complete predecessors and expose 75 complete
successors. Immediate root, clip-field, aid and frame observers see only a
complete state; failed predecessors retry; callback exceptions drain; reentry
retains the later publication. The exact CXX-SOURCE-067 library fails the new
complete-state assertion. Focused ASan and adjacent camera, affordance,
environment, host, GED and Qt checks pass. The exact full compact executable
passes 163 direct checks, 21 scenarios and 1,073 PASS rows in 356.77 seconds.
Evidence is retained at
`.build/obol-qualification/20260912-cutting-plane-input`.

The scene-light setters, Qt widget lifetime, cross-controller feature
composition and the wider S4--S6 matrix remain open. This bounded S1/S3 closure
does not materially change the conservative ranges below.

## September 12: complete scene-light input publication

`CXX-SOURCE-069` closes controller scene-light content and enablement as one
input family. Previously, content entered the private vector before retained
children were rebuilt and without a frame request. Enablement changed private
intent and each live light field in sequence before requesting a frame. GED
then invoked both setters consecutively, exposing another intermediate
generation. Equal input still rebuilt nodes or requested a redundant frame.

The controller now prepares a copied value snapshot, every detached child or
retained-node scalar edit, and one capacity-level render request before an
allocation-free commit. Private state, child or field state and the standing
request are complete before callbacks. Enablement preserves node identity;
notification state is restored before delivery; callback failures drain; and
reentry retains the later successor. Exact content and enablement allocate and
notify nothing. A combined public overload lets GED publish discovered lights
and their policy bit in one operation.

The focused matrix covers installation, replacement, clearing, nonempty and
empty enablement, combined changes and ordinary/automatic LoD. Its 443 forced-
allocation positions retain 392 complete predecessors and expose 51 complete
successors. Immediate graph, path, light-field and frame observers see only a
complete state; failed predecessors retry; callback exceptions drain; both
replacement and enablement reentry preserve their later successor. The saved
CXX-SOURCE-068 library lacks the combined public entry, while its source
snapshot records the preceding sequential implementations. Focused ASan and
adjacent scene-light, camera, environment, host, view-store, GED and Qt checks
pass. The exact full compact executable passes 164 direct checks, 21 scenarios
and 1,074 PASS rows in 374.97 seconds. Evidence is retained at
`.build/obol-qualification/20260912-scene-light-input`.

The enclosing camera-rig/profile lighting synchronizer, Qt widget lifetime,
cross-controller feature composition and the wider S4--S6 matrix remain open.
This bounded S1/S3 closure does not materially change the conservative ranges
below.

## September 12: complete lighting-state input publication

`CXX-SOURCE-070` closes the remaining multi-call lighting synchronizer.
Previously, one `bv_lighting_state` passed through profile, offset, tracking
and headlight setters before GED published scene lights. Observers and failures
could therefore separate one source record into several controller
generations.

One additive controller operation now prepares environment and camera-light
fields, copied scene-light intent, replacement children or enablement fields,
and one capacity render request. Its allocation-free commit publishes the
complete private, retained and frame state before callbacks. Existing leaf
setters use the same owner. Stable direction canonicalization prevents an
already normalized offset from drifting during an unrelated setter. Exact
input is allocation-free; invalid profiles and degenerate or non-finite
directions do no work; callback failures drain; and reentry retains its later
successor.

The fresh complete-state matrix covers ordinary and automatic LoD. Its 109
forced-allocation positions retain 98 complete predecessors and expose 11
complete successors. Immediate environment, camera-light, scene, path and
frame observers see only complete state; retries pass. Direct checks exercise
the profile, offset, tracking and headlight leaf APIs and their exact no-ops.
The saved CXX-SOURCE-069 library lacks the additive complete-state symbol.
Focused ASan, adjacent controller, renderer, LoD, host, GED and Qt checks pass.
The exact full compact executable passes 165 direct checks, 21 scenarios and
1,075 PASS rows in 347.54 seconds. Evidence is retained at
`.build/obol-qualification/20260912-lighting-state-input`.

Master shading and the remaining renderer-appearance inputs, Qt widget
lifetime, cross-controller feature composition and the wider S4--S6 matrix
remain open. This bounded S1/S3 closure does not materially change the
conservative ranges below.

## September 12: master-lighting input publication

`CXX-SOURCE-071` closes the master shading setter. Previously, one call wrote
the retained light model and three camera-light fields sequentially, then
published renderer invalidation and a second frame request. Field observers
could see a partial rig, and the callback reason did not match the final
standing request.

The setter now prepares detached scalar candidates and a reason-aware renderer
invalidation before mutation. Its allocation-free commit publishes flat/PHONG
shading, the complete profile-dependent rig, cleared renderer timing/capacity
evidence, automatic-LoD work and one typed `lighting` frame before callbacks.
Exact input allocates nothing and preserves renderer evidence; unchanged scene
lights retain their node identity; callback failures drain; and reentry retains
the later successor. The prior private invalidation entry remains available to
its existing callers.

The Studio/MGED matrix covers ordinary and automatic LoD. Its 120 forced-
allocation positions retain 104 complete predecessors and expose 16 complete
successors. Immediate light-model, camera-light and frame observers see only a
complete state. The saved CXX-SOURCE-070 library fails this selector. Focused
ASan, adjacent publication, LoD, renderer, host, GED and Qt checks, the public
umbrella and static link pass. The exact full compact executable passes 166
direct checks, 21 scenarios and 1,076 PASS rows in 371.20 seconds. Evidence is
retained at
`.build/obol-qualification/20260912-master-lighting-input`.

The remaining renderer-appearance inputs, Qt widget lifetime,
cross-controller feature composition and the wider S4--S6 matrix remain open.
This bounded S1/S3 closure does not materially change the conservative ranges
below.

## September 12: scalar appearance input publication

`CXX-SOURCE-072` closes background, depth-test, headlight color/intensity and
depth-cue setters. Previously, they wrote private or retained fields before a
frame request; depth-related setters also published a separately named
renderer invalidation, headlight values were misclassified as capacity work,
and an equal fog mode concealed stale derived fog fields.

The setters now share one prepared scalar transaction. Each detached candidate
and typed frame or renderer invalidation is ready before one allocation-free
commit. Background private colors and fog color, depth test/write, headlight
fields, fog mode/color/visibility, renderer evidence and the frame level are
complete before callbacks. Background and headlight values preserve renderer
timing as presentation-only changes. Depth-test and depth-cue clear timing
evidence and publish automatic-LoD capacity work only under that policy. Exact
input allocates nothing; non-finite input does no work; stale derived fog state
is repaired; callback failures drain; and reentry retains the later successor.

The seven-operation ordinary/automatic matrix covers 128 forced-allocation
positions: 96 retain complete predecessors and 32 expose complete successors.
Immediate environment, depth-buffer, headlight and frame observers see only a
complete state. The saved CXX-SOURCE-071 library fails this selector. Focused
ASan, nine adjacent publication selectors, ten adjacent CTests, the public
umbrella and static link pass. The exact full compact executable passes 167
direct checks, 21 scenarios and 1,077 PASS rows in 338.31 seconds. Evidence is
retained at
`.build/obol-qualification/20260912-scalar-appearance-input`.

Transparency, antialiasing, clip bounds and software-wire mode remain in the
renderer-appearance family. Qt widget lifetime, cross-controller feature
composition and the wider S4--S6 matrix also remain open. This bounded S1/S3
closure does not materially change the conservative ranges below.

## September 12: render-action appearance input publication

`CXX-SOURCE-073` closes transparency and antialiasing setters. Previously,
they changed controller and render-action state before a separately named
renderer invalidation and final request. Antialiasing also scheduled a toolkit
redraw outside the controller frame owner, while exact private state concealed
stale action state.

Both setters now prepare one reason-matched renderer invalidation before an
allocation-free commit of controller intent, onscreen and cached-offscreen
actions, renderer evidence and the frame level. Direct nonallocating action
setters enforce the named single-pass antialiasing policy and leave the typed
controller request as the sole redraw owner. Exact checks include all available
actions and repair stale backend state; callback failures retain the committed
state, and reentry retains the later successor.

The four change/repair operations in ordinary and automatic LoD cover eight
allocation-free operation states. The saved CXX-SOURCE-072 library fails this
selector. Focused ASan, ten adjacent publication selectors, ten adjacent
CTests, the public umbrella and static link pass. The exact full compact
executable passes 168 direct checks, 21 scenarios and 1,078 PASS rows in 356.16
seconds. Evidence is retained at
`.build/obol-qualification/20260912-render-action-appearance-input`.

Clip bounds and software-wire mode remain in the renderer-appearance family.
Qt widget lifetime, cross-controller feature composition and the wider S4--S6
matrix also remain open. This bounded S1/S3 closure does not materially change
the conservative ranges below.

## September 12: clip-bound input publication

`CXX-SOURCE-074` closes camera-relative clip bounds. Previously, the setter
changed only its private minimum and maximum, published a separately named
renderer invalidation and then requested a `clip-bounds` frame. The retained
clip nodes continued to hold their preceding world-space planes until a later
camera synchronization, and equal private values concealed stale planes.

Camera synchronization and the setter now use one finite plane derivation. The
setter prepares detached successors for the changed clip nodes and one reason-
matched renderer invalidation, then commits both derived planes, private bounds,
renderer evidence and the frame level before callbacks. Clip enablement remains
unchanged in active and disabled views, and a one-sided change preserves the
other plane. Exact input allocates nothing; invalid and unrepresentable input
does no work; stale retained state is repaired; callback failures drain; and
reentry retains the later successor. The shared multi-node scalar helper
removes the duplicate camera-local implementation.

The five active/disabled/change/repair operations in ordinary and automatic
LoD cover 104 forced-allocation positions: 80 retain complete predecessors and
24 expose complete successors. Immediate clip-node, clip-field and frame
observers see only complete state. The saved CXX-SOURCE-073 library fails this
selector. Focused ASan, twelve adjacent publication selectors, ten adjacent
CTests, the public umbrella and static link pass. The exact full compact
executable passes 169 direct checks, 21 scenarios and 1,079 PASS rows in 323.52
seconds. Evidence is retained at
`.build/obol-qualification/20260912-clip-bounds-input`.

Software-wire mode is the remaining renderer-appearance input. Qt widget
lifetime, cross-controller feature composition and the wider S4--S6 matrix
also remain open. This bounded S1/S3 closure does not materially change the
conservative ranges below.

## September 12: software-wire input publication

`CXX-SOURCE-075` closes the renderer-appearance input family. Previously,
software-wire mode was copied into controller intent, the LoD wrapper and the
compact render batch, but equality consulted only controller intent. A stale
render path could therefore survive an exact setter call. A change also
published a `renderer-performance` callback before renaming the standing frame
request to `software-wire-mode`.

The setter now checks every available runtime copy and prepares one reason-
matched renderer invalidation before committing controller, LoD-wrapper and
compact-batch policy, renderer evidence and the frame level. Exact input does
no work while equal controller intent repairs stale render paths. Invalid
input selects AUTO, an absent LoD wrapper retains a working compact path,
callback failure preserves the committed generation, and reentry retains the
later successor.

Seven mode/change/repair/render-path operations in ordinary and automatic LoD
cover 56 forced-allocation positions: 28 retain complete predecessors and 28
expose complete successors. The saved CXX-SOURCE-074 library fails this
selector. Focused ASan, thirteen adjacent publication selectors, ten adjacent
CTests, the public umbrella and static link pass. The source suite passes 275
PASS rows in 139.29 seconds. The exact full compact executable passes 170
direct checks, 21 scenarios and 1,080 PASS rows in 320.72 seconds. Evidence is
retained at
`.build/obol-qualification/20260912-software-wire-input`.

The renderer-appearance family is closed. Qt widget lifetime, cross-controller
feature composition, remaining traversal/lifecycle inputs and the wider
S4--S6 matrix remain open. This bounded S1/S3 closure does not materially
change the conservative ranges below.

## September 12: controller edit-preview input publication

`CXX-SOURCE-076` closes the direct controller edit-preview boundary. Existing
previews previously received identity, intent and revision fields before
geometry preparation could fail, and successful replacement, insertion and
removal notified the scene graph before creating the frame obligation.

The controller now builds a complete detached preview candidate first. For an
existing preview, prepared scalar and child-list operations preserve node
identity while committing all requested/realized fields, status, geometry and
the capacity-classified frame before observation. Initial insertion and
removal use prepared root child edits. The preview's own stale-marking sensors
are handled as part of the transaction while external field, node and path
observers remain active. Callback failures drain and reentry retains the later
successor; invalid input and missing removal do no work.

Insert, replace and remove in ordinary and automatic LoD cover 364 stable
forced-allocation positions: 350 retain complete predecessors and 14 expose
complete successors. The saved CXX-SOURCE-075 library fails this selector.
Focused ASan, fourteen adjacent publication selectors, ten adjacent CTests,
the public umbrella and static link pass. The source suite passes 276 PASS rows
in 140.46 seconds. The exact full compact executable passes 171 direct checks,
21 scenarios and 1,081 PASS rows in 320.73 seconds. Evidence is retained at
`.build/obol-qualification/20260912-controller-edit-preview-input`.

Direct line-layer/HUD-label controller inputs, Qt widget lifetime, cross-
controller feature composition, remaining traversal/lifecycle inputs and the
wider S4--S6 matrix remain open. This bounded S1/S3 closure does not materially
change the conservative ranges below.

## September 12: controller overlay input publication

`CXX-SOURCE-077` closes the direct controller line-layer and HUD-label input
family. Line-layer operations previously held their detached candidate without
reference ownership and changed/notified the root before creating the frame
obligation. HUD replacement wrote live fields before allocating replacement
geometry, and HUD insert/remove also notified the graph before the frame.

Both paths now construct complete reference-owned candidates first. Line-layer
insert/replace prepares the whole root order and removal prepares retirement.
HUD insertion uses the same root transaction; replacement preserves the live
HUD overlay identity while prepared scalar fields, child geometry and the
cached label pointer commit together. Removal is prepared. Every operation
commits its capacity-classified frame before graph or field callbacks. The
shared callback drain attempts all notifications while retaining the first
exception; invalid input does no work and reentry retains the later successor.

Line insert/replace/remove in ordinary and automatic LoD cover 348 stable
forced-allocation positions: 322 complete predecessors and 26 complete
successors. HUD operations cover 538 positions: 490 predecessors and 48
successors. The saved CXX-SOURCE-076 library fails both selectors. Focused ASan,
fourteen adjacent publication selectors, ten adjacent CTests, the public
umbrella and static link pass. The source suite passes 278 PASS rows in 140.33
seconds. The exact full compact executable passes 173 direct checks, 21
scenarios and 1,083 PASS rows in 320.43 seconds. Evidence is retained at
`.build/obol-qualification/20260912-controller-overlay-inputs`.

Remaining manipulator/navigation field inputs, Qt widget lifetime, cross-
controller feature composition, traversal/lifecycle inputs and the wider
S4--S6 matrix remain open. This bounded S1/S3 closure does not materially
change the conservative ranges below.

## September 12: grid view-input publication

`CXX-SOURCE-078` closes `bobol_grid_configure_from_view()` and its context
overload. The helper previously wrote roughly twenty live grid input fields
before an allocating rebuild published derived spacing, segment counts and HUD
geometry. Observers and allocation failure could retain arbitrary new-input/
old-output combinations.

The helper now copies retained grid policy into a reference-owned detached
candidate, applies the complete sanitized view input and rebuilds the candidate
off graph. Prepared scalar fields and a prepared child replacement then commit
all input, derived and geometry state while preserving the attached grid node.
Child and field callbacks drain while retaining the first exception. Invalid
input does no work and reentry retains the later successor.

Visible, hidden and snap-only configurations cover 9,422 stable forced-
allocation positions: 229 complete predecessors and 9,193 complete successors.
The high successor count sweeps every allocating post-commit field callback;
each retains the completed generation and remaining notifications. The saved
CXX-SOURCE-077 library fails the selector. Focused ASan, seven adjacent
publication selectors, ten adjacent CTests, the public umbrella and static link
pass. The source suite passes 279 PASS rows in 145.02 seconds. The exact full
compact executable passes 174 direct checks, 21 scenarios and 1,084 PASS rows
in 323.87 seconds. Evidence is retained at
`.build/obol-qualification/20260912-grid-input-publication`.

Axes/ADC and other manipulator/navigation input composition, image-source
mutation, Qt widget lifetime, remaining traversal/lifecycle work and the wider
S4--S6 matrix remain open. This bounded S1/S3 closure does not materially
change the conservative ranges below.

## September 12: axes and ADC view-input publication

`CXX-SOURCE-079` closes the retained input boundary for the fields represented
by `SoBRLAxes` and `SoBRLADC`. Their geometry rebuilds were already atomic, but
callers could only write several attached public fields and then rebuild.
Observers or allocation failure could leave a new input prefix paired with the
preceding geometry.

The new `bobol_axes_configure_from_view()` and
`bobol_adc_configure_from_view()` APIs build the supported state mapping and
geometry on a detached candidate, then commit prepared scalar fields and the
complete child list while preserving the target node. A shared private helper
now owns this commit and callback-drain protocol for axes, ADC and grid. The Qt
fallback and GED ADC producer use the canonical mapping. Invalid pointers do
no work, callback failures leave the complete successor and drain remaining
notifications, and reentry retains the later successor.

Axes and ADC visible/hidden configurations cover 222 stable forced-allocation
positions: 200 complete predecessors and 22 complete successors. The exact
CXX-SOURCE-078 library fails the selector on a partial input observation.
Focused ASan, eight adjacent publication selectors, ten adjacent CTests, the
two exported symbols, public umbrella and static link pass. The source suite
passes 95 direct checks, 21 scenarios and 280 PASS rows in 144.70 seconds. The
exact full compact executable passes 175 direct checks, 21 scenarios and 1,085
PASS rows in 323.13 seconds. Evidence is retained at
`.build/obol-qualification/20260912-overlay-input-publication`.

This closes publication for axes draw/position/size and ADC draw/center/angle/
distance/colors/width while retaining ADC crosshair/tick sizing. Rich axes
labels, ticks and style remain separate capability scope. Other manipulator/
navigation input composition, image-source mutation, Qt widget lifetime,
remaining traversal/lifecycle work and the wider S4--S6 matrix remain open.
This bounded S1/S3 closure does not materially change the conservative ranges
below.

## September 13: navigation gizmo input publication

`CXX-SOURCE-080` closes the supported input methods on
`SoBRLNavigationGizmo`. Hover, active and camera changes previously mutated live
state before allocating the replacement HUD. Camera replacement also detached
the preceding sensor before the new sensor, camera ownership and geometry were
known to be ready. Allocation failure and immediate observers could therefore
leave or see mixed generations.

The methods now build a reference-owned detached candidate and prepare its
scalar fields and complete child list before changing the attached gizmo.
Private rotation, HUD and anchor pointers commit in the same quiet phase.
Camera changes prepare the new retained camera and sensor before swapping
ownership; an unrealized gizmo commits only that private state. Fixed snapshots
and traversal-driven rotation use the same publication path. Exact no-ops do
no work, callback failures drain, and reentry retains the later successor.

Realized hover, active and camera add/replace/remove plus unrealized camera
add/replace/remove cover 5,172 stable forced-allocation positions: 5,162
complete predecessors and 10 complete successors. The exact CXX-SOURCE-079
library fails the selector on a partial input observation. Focused ASan, nine
adjacent publication selectors, ten adjacent CTests, four public method
symbols, the public umbrella and static link pass. The source suite passes 96
direct checks, 21 scenarios and 281 PASS rows in 169.71 seconds. The exact full
compact executable passes 176 direct checks, 21 scenarios and 1,086 PASS rows
in 374.37 seconds. Evidence is retained at
`.build/obol-qualification/20260913-navigation-gizmo-input-publication`.

This closes supported method composition and the fixed-snapshot production
producer. Direct public style-field writes, renderer-derived anchor placement,
image-source mutation, Qt widget lifetime, remaining traversal/lifecycle work
and the wider S4--S6 matrix remain open. This bounded S1/S3 closure does not
materially change the conservative ranges below.

## September 13: image-source input publication

`CXX-SOURCE-081` closes the supported source methods on `SoBRLImageSource`.
Stream replacement, image conversion, clear and refresh previously released
ownership or changed live scalar fields before the complete successor was known
to be publishable. Allocation failure could lose the predecessor, immediate
observers could see mixed metadata, and reattaching an owned stream through
`setStream(getStream())` could reuse a pointer after destroying it.

The source now prepares scalar metadata and a subscription owner off target.
One quiet commit installs stream ownership, subscription identity, pending and
realized generations and all public fields before retiring the predecessor or
delivering callbacks. Candidate dirty callbacks remain isolated until commit.
Settled refreshes, same-stream attachment and repeated clear/failure states are
exact no-ops. Callback failures drain, and reentry retains the later successor.
The ownership pointer occupies the preceding stream-pointer slot; a compiled
before/current probe confirms the class remains 928 bytes.

Stream add, dirty/clean stream replacement, static-image replacement, clear,
invalid stream, invalid image and refresh cover 174 stable forced-allocation
positions: 143 complete predecessors and 31 complete successors. Owned-stream
self-reattachment is also an exact no-op. The exact CXX-SOURCE-080 library fails
the selector on a partial input observation. Focused ASan, nine adjacent
publication selectors, ten adjacent CTests, the public umbrella and static link
pass. The source suite passes 97 direct checks, 21 scenarios and 282 PASS rows.
The source run completes in 171.80 seconds. The exact full compact executable
passes 177 direct checks, 21 scenarios and 1,087 PASS rows in 373.51 seconds.
Evidence is retained at
`.build/obol-qualification/20260913-image-source-input-publication`.

This closes the supported image-source method boundary. Direct public field
writes, composition with a display consumer, Qt widget lifetime, remaining
traversal/lifecycle work and the wider S4--S6 matrix remain open. This bounded
S1/S3 closure does not materially change the conservative ranges below.

## September 13: framebuffer cursor input publication

`CXX-SOURCE-082` closes the three framebuffer cursor-state methods on
`BObolWindowHost`. They previously changed cursor visibility, image position
and cursor shape in separate observable writes, followed by a presentation
request. A field or node observer could see a mixed state without a standing
frame, and a throwing callback could stop the remaining publication.

The methods now share one fixed three-field transaction. It prepares the owned
presentation reason before changing live state, commits all cursor fields and
the standing presentation request quietly, then drains field and endpoint
callbacks while preserving the first exception. Reentry retains the later
successor and exact input does no work. Preparation is limited to the three
changed fields and one render-request value.

Cursor enable, screen-cursor disable, custom-shape installation and shape reset
under ordinary and automatic LoD cover 44 stable forced-allocation positions:
4 complete predecessors and 40 complete successors. The exact CXX-SOURCE-081
library fails on its first partial cursor observation. Focused ASan, ten
adjacent publication selectors, ten adjacent CTests, public method symbols,
the public umbrella and static link pass. The source suite passes 98 direct
checks, 21 scenarios and 283 PASS rows in 155.23 seconds. The exact full
compact executable passes 178 direct checks, 21 scenarios and 1,088 PASS rows
in 335.61 seconds. Evidence is retained at
`.build/obol-qualification/20260913-framebuffer-cursor-input-publication`.

This closes cursor-state publication. Custom bitmap rendering, framebuffer
reset, viewport/view input, source-to-display refresh, direct viewport fields,
endpoint cleanup, Qt widget lifetime and the wider S4--S6 matrix remain open.
This bounded S1/S3 closure does not materially change the conservative ranges
below.

## September 13: framebuffer view input publication

`CXX-SOURCE-083` closes `BObolWindowHost::setFramebufferView()`. The method
previously published source center, zoom, refreshed source/geometry and the
presentation request in sequence. Observers and allocation failure could retain
a mixed transform/geometry/request state, and a view-only operation could
consume a pending image generation before its display publication was ready.

The host now edits a detached viewport candidate and prepares the presentation
request before commit. A retained-payload path copies the accepted texture bytes
and realized revisions, and one shared internal builder derives the successor
HUD and texture coordinates. Scalar fields, child geometry, cache pointers and
the standing request commit before callbacks. Pending pixels and image-source
fields remain unchanged for flush. Exact input is an allocation-free no-op,
callback failures drain, and reentry retains the later view and its valid
in-memory texture. The implementation adds no public layout or exported helper
surface.

Combined pan/zoom, pan-only and zoom-only operations under ordinary and automatic
LoD cover 1,322 stable forced-allocation positions: 1,308 complete predecessors
and 14 complete successors. The exact CXX-SOURCE-082 library fails the focused
contract. Focused ASan, ten adjacent publication selectors, ten adjacent CTests,
the public symbol/umbrella checks and static link pass. The source suite passes
99 direct checks, 21 scenarios and 284 PASS rows in 154.94 seconds. The exact
full compact executable passes 179 direct checks, 21 scenarios and 1,089 PASS
rows in 337.67 seconds. OSMesa remains `097aac27...`. Evidence is retained at
`.build/obol-qualification/20260913-framebuffer-view-input-publication`.

This closes framebuffer view publication. Framebuffer reset and viewport input,
source-to-display flush publication, direct viewport fields, endpoint cleanup,
Qt widget lifetime and the wider S4--S6 matrix remain open. This bounded S1/S3
closure does not materially change the conservative ranges below.

## September 13: framebuffer reset input publication

`CXX-SOURCE-084` closes `BObolWindowHost::resetFramebuffer()`. The method
previously changed center, zoom and cursor visibility one field at a time,
refreshed source and viewport state, then requested presentation. A failure or
observer could expose a partial reset, and reset could consume pending pixels
outside the display refresh boundary.

Reset now checks the exact target state before allocation, edits the detached
viewport candidate and reuses the retained viewport transaction from the view
boundary. Reset fields, complete HUD, cached texture/face pointers and the
standing request commit before callbacks. Currently displayed bytes and
revisions remain unchanged, while pending stream data stays pending for flush.
The reset and view methods share one commit and callback-drain implementation.

Ordinary and automatic LoD cover 442 stable forced-allocation positions: 436
complete predecessors and 6 complete successors. The exact CXX-SOURCE-083
library fails on a partial reset observation. Focused ASan, ten adjacent
publication selectors and ten adjacent CTests pass. Exact repeated reset is
allocation- and notification-free; callback failure and reentry pass. The
source suite passes 100 direct checks, 21 scenarios and 285 PASS rows in
155.17 seconds. The exact full compact executable passes 180 direct checks,
21 scenarios and 1,090 PASS rows in 337.34 seconds. Public symbols, the public
umbrella and static link pass; OSMesa remains `097aac27...`. Evidence is retained at
`.build/obol-qualification/20260913-framebuffer-reset-input-publication`.

This closes framebuffer reset publication. Framebuffer viewport input,
source-to-display flush publication, direct viewport fields, endpoint cleanup,
Qt widget lifetime and the wider S4--S6 matrix remain open. This bounded S1/S3
closure does not materially change the conservative ranges below.

## September 13: framebuffer viewport input publication

`CXX-SOURCE-085` closes `BObolWindowHost::setFramebufferViewport()`. The old
method published controller size before framebuffer position/size and HUD
geometry, then requested another frame. Endpoint and scene observers could see
different generations, failure could retain a split state, and placement could
consume pending pixels outside flush.

The host now prepares retained HUD geometry before composing the controller
size transition. One quiet commit installs framebuffer position/size, complete
geometry and cache pointers, all controller region/LoD state and one
`fb-viewport` request before callbacks. Size changes retain the controller's
capacity semantics; position-only changes request presentation. Exact input
does no work, and pending source data remains pending for flush.

Move/resize, move-only, resize-only and invalid-extent fallback in ordinary and
automatic LoD cover 1,776 stable forced-allocation positions: 1,756 complete
predecessors and 20 complete successors. The exact CXX-SOURCE-084 library fails
on a partial viewport observation. Focused ASan, ten adjacent publication
selectors and ten adjacent CTests pass. Failure retry, callback drain and
reentry preserve one complete composed generation. The source suite passes 101
direct checks, 21 scenarios and 286 PASS rows in 157.54 seconds. The exact full
compact executable passes 181 direct checks, 21 scenarios and 1,091 PASS rows
in 338.82 seconds. Public symbols, the public umbrella and static link pass;
OSMesa remains `097aac27...`. Evidence is retained at
`.build/obol-qualification/20260913-framebuffer-viewport-input-publication`.

This closes framebuffer viewport publication. Source-to-display flush
publication, direct viewport fields, endpoint cleanup, Qt widget lifetime and
the wider S4--S6 matrix remain open. This bounded S1/S3 closure does not
materially change the conservative ranges below.

## September 13: framebuffer source-to-display flush publication

`CXX-SOURCE-086` closes `BObolWindowHost::flushFramebuffer()`. The old method
published source refresh, viewport synchronization and its frame request as
three separate observable operations. A callback or preparation failure could
leave one accepted source generation paired with old displayed pixels, old
realized revisions or no standing frame obligation.

Flush now reuses the private image-source publication and the detached
viewport builder. It prepares source metadata, a stream payload whose metadata
is stable across the read, complete HUD geometry/cache state and one
`fb-flush` request before committing any live field. Source and viewport
notification gates are both restored before the first callback, and all
source, viewport/path and endpoint notifications drain after the complete
commit. Exact settled flush is a nonallocating no-op.

Dirty-pixel flush in ordinary and automatic LoD covers 486 stable
forced-allocation positions: 476 complete pending predecessors and 10 complete
displayed successors. The exact CXX-SOURCE-085 library fails on a partial
source observation. Focused ASan, eleven adjacent publication selectors and
ten adjacent CTests pass. Failure retry, callback drain and callback-time later
pixel publication preserve one complete generation. The source suite passes
102 direct checks, 21 scenarios and 287 PASS rows in 156.55 seconds. The exact
full compact executable passes 182 direct checks, 21 scenarios and 1,092 PASS
rows in 338.58 seconds. Separate implementation-unit compilation, public
symbols, the public umbrella and static link pass; OSMesa remains
`097aac27...`. Evidence is retained at
`.build/obol-qualification/20260913-framebuffer-flush-input-publication`.

This closes framebuffer source-to-display refresh. Direct viewport fields,
endpoint cleanup, Qt widget lifetime and the wider S4--S6 matrix remain open.
This bounded S1/S3 closure does not materially change the conservative ranges
below.

## September 13: retained RT image publication

`CXX-SOURCE-087` closes retained RT first-image and restart publication. The
old endpoint attached an empty viewport, then published source refresh and
viewport synchronization separately. Source/root observers could therefore
see mismatched generations. A nested restart from an immediate callback could
also leave a newer worker joinable when the older invocation tried to assign
its worker, terminating the process.

RT stream, source and viewport resources now prepare under scoped ownership.
One internal transaction stages the exact stream generation, source metadata,
stable pixels, complete viewport geometry/cache state, optional first root
attachment and the `rt-restart` request before a quiet commit. Both node
notification gates restore before callbacks, and all source, viewport/path,
root and request notifications drain. Failed restart preparation retains the
exact staged payload for retry only while the private stream generation still
matches. Dequeued worker frames are restored when publication has not
committed, and a restart generation makes callback-time successors supersede
older stack frames. RT numeric/quality properties roll back after a rejected
restart and exact values do no work.

Restart in ordinary and automatic LoD covers 538 stable forced-allocation
positions: 474 complete predecessors and 64 complete successors. Two observed
first activations attach only complete images. The exact CXX-SOURCE-086 library
fails on the empty first attachment. Focused ASan, ten adjacent publication
selectors and ten adjacent CTests pass. Failure retry, callback drain and
callback-time nested restart preserve one complete generation. The source
suite passes 103 direct checks, 21 scenarios and 288 PASS rows in 165.80
seconds. The exact full compact executable passes 183 direct checks, 21
scenarios and 1,093 PASS rows in 351.12 seconds. Independent implementation-
unit compilation, public symbols, the public umbrella and static link pass;
OSMesa remains `097aac27...`. Evidence is retained at
`.build/obol-qualification/20260913-rt-image-publication`.

This closes the retained RT image source-to-display boundary. Renderer-engine
selection, endpoint teardown, direct public viewport fields, Qt widget lifetime
and the wider S4--S6 matrix remain open. This bounded S1/S3 closure does not
materially change the conservative ranges below.

## September 13: renderer-engine publication

`CXX-SOURCE-088` closes the renderer-policy transition boundary. The old
setter changed the public engine, destroyed outgoing RT state and only then
performed allocation-capable controller invalidation and render-manager work.
Failures could expose a target policy with the preceding root, or the old
policy after losing its RT presentation.

Controller invalidation and outgoing RT root removal now prepare first. RT
activation joins the retained image transaction through a precommit hook, so
engine policy, graphical eligibility, performance invalidation, source,
viewport, root attachment and the capacity request are all live before any
observer runs. An outgoing RT state is retired from the endpoint at commit but
held until notification stacks finish, allowing callback-time replacement to
install a new RT successor safely. Allocation-capable Coin synchronization and
non-graphical request clearing retain explicit postcommit retry obligations. A
graphical frame notification remains deferred until synchronization succeeds.
Once settled, an exact repeated selection does no work.

Twelve transitions across all six engine policies in ordinary and automatic
LoD cover 2,887 forced-allocation positions: 1,932 complete predecessors and
955 complete successors. Root/path observers, retry, callback exceptions and
both directions of callback-time RT replacement pass. The exact
CXX-SOURCE-087 library exposes a partial engine/root generation. Focused ASan
and adjacent endpoint/host checks pass. The source suite passes 104 direct
checks, 21 scenarios and 289 PASS rows in 168.70 seconds. The exact full
compact executable passes 184 direct checks, 21 scenarios and 1,094 PASS rows
in 349.59 seconds. Independent implementation-unit compilation, public
symbols, the public umbrella and static link pass; OSMesa remains
`097aac27...`. Evidence is retained at
`.build/obol-qualification/20260913-renderer-engine-publication`.

This closes renderer-engine selection. Endpoint teardown, direct public
viewport fields, Qt widget lifetime and the wider S4--S6 matrix remain open.
This bounded S1/S3 closure does not materially change the conservative ranges
below.

## September 13: endpoint terminal destruction

`CXX-SOURCE-089` closes endpoint teardown. The old destructor depended on
live publication and cancellation paths: an observer exception or allocation
failure could skip later cleanup, cross the C or C++ destruction boundary, or
leave controller work live. A destroy request made by a synchronous endpoint
callback could delete state still in use by the outer public operation.

The terminal path now contains every exception and independently retires
presentation state, automatic-LoD generations, render requests, capacity
claims and progressive work. An owner-thread operation scope defers final
deletion until a callback stack drains, while terminal endpoints reject new
mutation. Host detach is no-throw, adopted RT streams belong to their image
source, and borrowed controllers receive a complete safe root before RT state
is released. Pending render-manager synchronization belongs to the controller
and can be settled by a replacement endpoint.

Owned and borrowed endpoints across `AUTO`, `NONE`, `DIAGNOSTIC` and `RT`
cover 30 forced-allocation positions and 9 self-contained borrowed-root
fallbacks. Throwing observers, callback-time destruction and reentry, an
active automatic-LoD generation, RT ownership, engine activation callbacks
and replacement-endpoint retry pass. The exact CXX-SOURCE-088 library exits 42
after terminal cleanup failure. Focused ASan, ten adjacent publication
selectors and ten adjacent CTests pass. The source suite passes 105 direct
checks, 21 scenarios and 290 PASS rows in 170.39 seconds. The exact full
compact executable passes 185 direct checks, 21 scenarios and 1,095 PASS rows
in 384.03 seconds. Four implementation units compile independently; public
symbols, the public umbrella and static linking pass. Evidence is retained at
`.build/obol-qualification/20260913-endpoint-destruction`.

This bounded S3/S4 closure removes endpoint teardown from the open lifecycle
list. Direct public viewport fields, Qt widget lifetime and the wider S4--S6
matrix remain open. The conservative ranges below are unchanged.

## September 13: stage completion estimates

These ranges apply the simplification guide's exit criteria to retained current evidence.
They estimate completed scope, not remaining time, and should not be averaged.

| Stage | Estimated completion | Main remaining obligation |
|---|---:|---|
| S1: ownership and acceptance inventory | 75--82% | Finish remaining traversal, lifecycle and caller-alias acceptance mapping |
| S2: one complete production path | 85--92% | Reconcile the completed source-evidence path with the final S1 inventory and declare its superseded state closed |
| S3: physical ownership boundaries | 55--68% | Close remaining traversal and service-lifecycle ownership and remove superseded writers |
| S4: resources, async work and lifecycle | 45--65% | Close planning/frame stalls and qualify memory, cancellation, BREP and shared-cache behavior |
| S5: capability qualification | 35--50% | Close large/cold scene, visual, native-host and shared-client rows |
| S6: release candidate qualification | 5--10% | Freeze one coherent stack and run the complete required matrix against its exact identities |

The structural simplification finish line, S1--S3, is roughly **68--77%** complete by
exit-criterion coverage. The production finish line remains roughly **40--55%** complete
because S4--S6 contain the largest timing, resource, visual and platform uncertainty.

## September 9: root and repository counterexamples

`CXX-ROOT-001` was reproduced against GROUP-008; those baseline binaries remain retained.
A working repair now passes the four rebuilt reproducers, but is not qualified. The working evidence is
`.build/obol-qualification/20260909-root-repository`; its README and audit status distinguish
this baseline audit from the preceding sealed shape-state checkpoint.

All four standalone probe builds pass. `root_repository_probe` exposes obsolete root and
revision state in both replacement/clear deletion callbacks, and lost cache retention in
three shared-controller cases. `root_repository_index_probe` additionally primes the
sibling index and reproduces a removed source remaining in the shared-root lookup.
`root_failure_probe` injects 38 allocation positions: 20 preserve the preceding state and
18 leave the next root installed with old revisions and incomplete cache ownership.
`root_destruction_probe` denies the first allocation during destruction; its terminate
handler confirms `std::bad_alloc` and exits 42 without a core dump. The other three
reproducers exit 1, as expected for the recorded baseline failures.

The working implementation prepares root, indexes and repository membership together,
separates controller ownership from source cache refresh, and uses native parent auditors
to propagate shared child edits. Source conflict retirement uses the same publication path.
The publication and GED draw-sync test targets compile (`build-root-3`, 43.54 s). All four
rebuilt reproducers exit 0: both retirement callbacks observe complete state, shared cache
and primed lookup cases pass, all **73 allocation positions preserve the old scene**, and
destruction retires its cache ownership with **zero allocations attempted**. All **62 selected
scene state checks pass in 177.14 s**, within the unchanged 180-second allowance but with
less margin than the preceding 148.65-second run. Investigate that cost before qualification.
The working audit status and `working-check-identities.json` record the checks and binary
identities; all tool handles are terminal and the final host process check is empty.
These focused checks do not yet supply
production regression integration, complete shared-writer/callback coverage or graphical
qualification. Live libBObol and static GED archives differ from GROUP-008.

The shared-controller follow-up adds a reproducible 18-row probe: 15 rename/metadata
rows fail on the working root repair, while three movement rows already pass. All 18
pass after composing peer revisions with field and hierarchy publication. Scalar source
setters use this same prepared peer context. Bulk material refresh prepares effects at
each target's existing acceptance boundary, so unchanged or retired peers do not advance.
Three new production checks pass: **513 scene allocation positions** (509 preserved,
four complete commits), **176 layout/batch/quiet/no-op/callback/conflict cases**, and
**122 material-sweep allocation positions plus four callback histories**. The working
audit status records subsequent broader regressions and baseline comparisons. These
checks do not close ROOT-001 or qualify new graphical/native behavior.

The first broader run exposed an older fixture's expectation that shared peers never
advance. Its material/geometry assertions are retained, with explicit revision checks
for every fixture controller. A subsequent run hit the unchanged 180-second limit.
Short graph traversals now use native inline storage, retain hashed membership for larger
graphs and cycle checks, and visit peer owners without a temporary result set. The
controlled source/group publication workload drops from 3,115 to 2,918 allocation positions
and from 7.06 to 6.56 s. The first compile of this optimization caught protected access to
the native list's capacity; the helper now extends that storage through private inheritance.

Final working-step evidence: **65 scene checks pass in 160.89 s**, GED draw synchronization
passes, **21 freshly linked probes** pass and **four translation units** compile independently.
The three strengthened production checks reject GROUP-008 with loader resolution retained;
the initial missing-library launch is recorded separately. Final shared sweeps cover **444
scene allocation positions** (441 preserved, three complete commits) and **106 material
allocation positions**, with the 176 edge cases and four material callback histories retained.
The root sweep now covers **61 allocation positions**, all preserving the old scene;
destruction still attempts no allocation. `shared-check-identities.json` records tested
inputs/binaries. All tool handles are terminal and the final host process check is empty.
Remaining repository operations, broader root/repository regressions and
full client/graphical qualification remain open. Stage estimates remain unchanged.

### Realization callback acceptance follow-up

The corrected standalone reproducer exposes four failures in the preceding working
implementation: missing peer frame publication, an obsolete tail after root detachment,
and crashes after current-source removal or repository replacement. The first probe
fixture realized an unretained source; those four initial crashes are retained as fixture
errors and are excluded from implementation evidence. The corrected inputs and runs are
under `realization-acceptance/`, with the distinction in `scope-and-attempts.json`.

Realization now retains the captured repository, root and current source. Each source
prepares peer effects with the existing scene/view commit. Scene hierarchy or repository
replacement stops every enclosing native traversal, including suspended nested calls and
explicit mutation batches. This prevents a child cursor or cache tail from continuing
against retired inputs. Committed progress survives interruption, and a subsequent call
visits the current scene. Source-owned descendant effects and compact-hierarchy cleanup
were retained for the separate follow-up recorded below.

Expanded callback coverage found a native `SoBase::getAuditors` defect: clearing a cached
list by increasing its index skipped entries as the list shrank, retaining a destroyed
peer's sensor. Clearing from the end fixes it. The new native test rejects the preceding
Obol library. Final native validation passes **1,577 unit tests** with one existing skip
and **63 integration tests**. The explicit workspace install preserves OSMesa's
`097aac27...` identity.

Two new production checks pass **64 callback histories**, **eight source-root cases**,
sole-owner root retirement and **189 allocation positions** with complete prefixes
`[103, 73, 13]`. Both checks reject the preceding working libBObol with current Obol.
The full suites pass **102 direct checks and 21 source scenarios**: source publication
in **89.20 s**, scene state in **161.52 s**, and GED synchronization in **8.11 s**.
Both publication limits remain 180 seconds. Five translation units compile independently, and all 22 freshly linked scene/GED/root/realization probes pass.
An early overlapping regression launch was cancelled; it is excluded, its process audit
is empty, and the sequential replacement completed successfully. The working status and
`realization-check-identities.json` retain final inputs and checks. Root/repository,
five-client/graphical and S1--S6 qualification remain open; stage estimates are unchanged.

### Source-descendant and cleanup acceptance follow-up

A controller rooted at a nested source previously retained stale source indexes after its
parent realized, and exact last-owner cache retirement failed. The cleanup loop also stored
only instance-key strings, so a callback replacement or reparenting with the same key could
be removed as if it were the original target. The corrected standalone probe retains these
four preceding failures under `source-descendants/`; both production selectors reject the
preceding working libBObol with loader resolution recorded.

Successful realization now prepares its exact child order together with descendant indexes,
repository ownership and all affected scene revisions. Changed retained auxiliary metadata
joins that same peer-effect collector. The source and scene effects commit before callbacks,
while the expected child edit does not stop its own traversal. Cleanup retains source and
ancestor nodes and revalidates their identity, parent chain and compact state immediately
before removal, so later callback edits remain current.

The new production coverage passes **12 descendant histories**, **eight changed/unchanged
auxiliary cases**, **one mixed success/failure sequence**, **six cleanup callback histories**,
**296 descendant allocation positions** (256 preceding states and 40 complete commits), and **248 cleanup allocation positions**
with complete prefixes `[178, 69, 1]`. The final suite contains **38 source checks plus 21
source scenarios** and **66 scene checks**; the combined run took **102.05 s**, **172.48 s**
and **8.26 s** for source, scene and GED synchronization. Six translation units compile
independently, and 23 rebuilt probe binaries provide 29 passing rows. These are working-step
results: ROOT-001 integration, fresh graphical
coverage and S1--S6 qualification remain open.

The mixed sequence caught and rejected an intermediate repair that kept the preceding
successful source's prepared hierarchy effects when a later source failed. The final
adapter clears its mutually exclusive child/peer state at every source boundary.

### Auxiliary and external late-source publication follow-up

Publication after initial scene construction previously bypassed the prepared scene-effect
boundary for auxiliary sources and external source data. A shared controller could therefore
observe new source children without the matching source-root, descendant-index, repository
or peer-revision state. Auxiliary removal had the same gap when callbacks changed a retained
source or its enclosing scene.

The controller now supplies one internal child-effect adapter to auxiliary insertion,
replacement and removal and to external line, point, mesh, annotation, clear, primitive-wire
and submodel publication. Each operation prepares its exact child order and all affected
indexes, repository membership and controller revisions before source mutation. Scene
effects commit before callbacks. A newly inserted source may establish its private identity
before preparation so traversal can discover it; it remains retained and unattached until
the complete publication commits.

The external selector passes **seven callback histories** and **494 allocation positions**
(492 preceding states and two complete commits). The auxiliary selector passes **eight
removal histories** and **281 allocation positions** (268 preceding states and 13 complete
commits). Both reject the saved immediately preceding libBObol. The final combined run
passes **40 source checks plus 21 scenarios in 116.59 s**, **66 scene checks in 180.84 s**
and GED synchronization in **8.00 s**. The exhaustive publication tests now have a documented
300-second scheduler allowance. Six translation units compile independently, and 23 rebuilt
probes provide all 29 passing rows.

CMake regeneration during this run restaged OSMesa `75dfd5f3...`. The independent renderer
and GED test both crashed in `_mesa_GetError`; GDB showed the public lookup resolving to an
incompatible internal entry point. Restoring the documented local OSMesa `097aac27...`
made the independent renderer and final GED row pass. The failed attempts are retained as
dependency-staging evidence and excluded from source-publication results. This recurrence
belongs to S6 coherent dependency staging; no source-publication defect is attributed to it.

The late-source step remains working evidence rather than a sealed checkpoint. Its next
repository operations are qualified below; ROOT-001 production integration, fresh graphical
coverage and S1--S6 qualification remain open. The stage estimates remain unchanged.

### Repository operation acceptance follow-up

The preceding repository cache released values while its maps were only partly updated.
A destruction callback could publish a later cache entry which the original loop then
removed. Rename mutated seven cache maps before every allocation and residency change had
succeeded. It also added old and new reference totals, overcounting a source which owned
both names and leaking the final cache lease. Scene wrappers invoked a raw pointer through
their member `shared_ptr`; a destruction callback which replaced that repository could
destroy the active callee and crash.

Clear and object/view invalidation now extract all affected entries before releasing them.
Rename prepares only the exact and variant keys plus affected source-owner records, retains
the staged geometry and commits with node handles after all allocation succeeds. Its
destination reference count is the union of owning sources. Retired values remain held
until residency commits and then invoke callbacks against complete state. Scene wrappers
retain the invoked repository and stop an active realization traversal before mutation,
preserving callback repository replacement.

`repository-cache-operations` passes **four direct cache callback histories**, allocation-
free clear/view invalidation, **eight rename allocation positions** and exact last-owner
retirement. `scene-repository-operations` passes **three repository-replacement callback
histories**. `realization-acceptance` passes **88 callback histories**, eight source-root
cases and sole-owner root retirement. The saved preceding libBObol rejects the first and
third selectors and segfaults on the isolated scene selector; all three current selectors
pass.

The final source suite passes **42 direct checks plus 21 scenarios in 113.23 seconds**.
The unchanged **66 scene checks pass in 179.04 seconds**. Static libBObol and GED were then
rebuilt; GED synchronization passes in **7.89 seconds**, six translation units compile
independently, and 23 freshly linked probes provide **29 passing rows**. The independent
renderer smoke test passes and OSMesa remains `097aac27...`. The earlier probe launch made
before rebuilding static libBObol is retained but superseded; it is excluded from final
probe evidence.

This closes the three repository operations at their focused boundary. ROOT-001 still
requires composed root/repository production integration, wider lifetime/failure coverage,
five-client and fresh graphical qualification. S1--S6 remain open and the stage estimates
below are unchanged.

The existing full builds and graphical passes remain dated GROUP-008 evidence. The new
counterexamples require combined root/index/repository publication, separate controller
ownership from source-data refresh, shared lookup propagation and allocation-free final
retirement. The [conformance entry](libbobol_tla_conformance.md#cxx-root-001--root-and-repository-ownership-retires-before-complete-scene-state)
owns these boundary obligations; active debt retains the wider release scope. The
[stage estimates](#september-9-stage-completion-estimates) remain unchanged and S1--S6 open.

## September 9: complete shape state

`CXX-GROUP-008` repairs `setShapeDrawState`, `setShapeDisplayState`,
`setShapePlacementState` and `setShapeSourceState`. All sixteen preceding-library probe
cases expose partial updates; all eight editing callbacks also lose their later edit.

One preparation path now stages only requested fields, retaining the original shape.
The existing scalar publication machinery commits values and one frame effect before
observers. Shared geometry, arrays and unrelated metadata keep their owners. Later
callback edits survive, including reentrant calls and removal/replacement. Structural
revisions, nullable strings, aliased inputs, existing tolerances and batch semantics
remain covered. The four live typed writers and duplicate mesh/wire dispatch wrappers
are removed; the scene controller is **287 net lines smaller**.

Evidence is retained in `.build/obol-qualification/20260909-shape-state`. Both new production
checks reject the preceding shared library with verified loader resolution. Sixty-four
allocation scenarios cover **1,472 injected positions**: 1,440 preserve the preceding scene
and 32 preserve the complete commit when notification fails. Sixty-four unconstrained
runs also pass. Coverage includes shared shape/geometry/cache ownership, node/field/path
observers, quiet fields/nodes, batches, throwing observers and tested allocation-free
no-ops. **72 callback cases** cover reentrant updates, removal/movement/replacement, root
replacement, quiet policy, aliased requests and unrelated metadata edits. Eight lifetime
cases remove targets during callbacks without an external reference. Separate cases
exercise aliased colors/strings/matrices, nullable strings and existing tolerances.

All five clients rebuild in **65.86 s**. All **22 focused checks** pass in
**148.70 s**, including **97 direct publication checks and 21 scenarios**, each
run once. All preceding 95 direct checks and 30 GED appearance cases remain present.
The controller compiles independently with production flags. All **16 rebuilt scene/GED
probes** pass, including complete shape-state observation, retained callback edits and the
preceding cached-bounds diagnostic. The first two build attempts caught a floating-point
comparison and a shadowed test variable; both were corrected before qualification.

All **24 graphical rows pass**, with **131.82 s** of summed driver execution:
two eager rows, ten delivery/rejection rows, eight cold/warm Generic Twin rows and four
primitive-edit rows with 207 events each. Eight camera contracts and 44 visible-HUD checks
pass. All **24 selected original captures were inspected**. Denied partial merges retain
a filled triangle, bounds and an explicit error cue; denied adoption retains four filled
triangles. Aircraft geometry and primitive-edit handles remain visible on both backends.
Text clipping/overlap and offscreen edit geometry remain S5 work; BoT edit captures do
not independently qualify surface quality. Private-Xvfb System GL uses CPU llvmpipe and
supplies no native GPU qualification.

`current-hashes.json` and `runtime.json` retain **58 runtime artifacts**,
**84 resolved shared libraries** and **141 physical files**, with loader resolution.
libBObol is `bee5a583...` and libged `4e28cd31...`.
Obol and OSMesa are unchanged; their native suites were not rerun. The formal baseline
and user's proposal are unchanged. Full TLC, shared-stack sanitizers, coherent dependency
staging and broad resource/native-host qualification remain release obligations.

Next audit root replacement and repository sharing/clearing, view-variant invalidation and
object rename, including cross-controller source ownership. The enclosing GED local shape
factories retain separate target/current-return acceptance work. Source/service and
presentation ownership follow. The [stage estimates](#september-9-stage-completion-estimates)
remain unchanged. **GROUP-008 closes the four shape-state setters; S1--S6 remain incomplete.**

## September 9: complete shape membership

`CXX-GROUP-007` repairs `moveShapeToGroup` and `removeShape`. The preceding library
exposes partial scene state in all four shape-probe cases: mesh/wire movement and
removal. Movement notifies after detaching the original edge and before adding the
destination edge; both operations notify before their structural/frame effect.

Both entry points now use `ChildPublication`. Movement stages both parent edits before
commit, and removal uses the existing selected-edge transaction. Membership, child paths
and one structural/frame effect commit before observers. Shape fields, shared geometry
and unrelated source/cache ownership remain intact. The original target stays alive
through notification, including callback removal during movement. Later style, membership,
parent/root and notification-policy changes remain current. The direct membership/late-
revision sequences are removed without another state owner. Existing lookup exclusion
for database-source subtrees, first-match path behavior and no-op/error semantics remain.

Evidence is retained in `.build/obol-qualification/20260909-shape-membership`. Both new
production checks reject the preceding shared library, with loader resolution recorded.
Sixty-four allocation scenarios cover **1,632 injected positions**, all preserving the
preceding scene, plus sixty-four unconstrained reference runs. They cover mesh/wire
shapes, shared geometry, shared nodes and parent groups, duplicate child positions,
node/path observation, batches and quiet fields. Throwing observers preserve the complete
publication and remaining notification. **176 callback cases** cover observation, style,
drain, movement, reinsertion/replacement, parent/root removal, quiet policy and aliased
requests. Four lifetime cases remove targets without an external reference, including
callback removal of a newly moved target. Separate cases cover colliding paths, invalid
requests, source-subtree exclusion and existing destination edges. The source-removal
and shape fixtures now share one independent path-ownership helper.

All five clients rebuild in **67.60 s**. All **22 focused checks** pass in
**143.97 s**, including **95 direct publication checks and 21 scenarios**, each
run once. All preceding 93 direct checks and 30 GED appearance cases remain present.
The scene controller compiles independently with production flags. Fourteen rebuilt
scene/GED probes pass; the original shape probe now observes complete state in all four
cases. An additional two-case bounds probe primes the origin's bounding-box query and
observes an empty origin from path callbacks after movement and removal. The first probe
compile omitted the concrete shape headers; that failed compile remains recorded and
was corrected before baseline execution.

The first combined GUI attempt ended with tool status **143**, without a runner
completion record. A surviving qged child subsequently exited; its stderr records a
broken X11 connection and an idle timeout at event 65. The cause of the runner termination
remains undetermined. `gui-interrupted-1` retains that attempt, the original-path mapping,
its unsuccessful final report and the eight initially inspected captures. It supplies
no complete graphical qualification. Host-visible checks confirmed that all its processes
were terminal before the fresh run.

The fresh **24 graphical rows pass**: two eager LoD-off/on rows, ten delivery/rejection
rows, eight cold/warm Generic Twin rows and four primitive-edit rows with 207 events each.
The stages run separately with their own completion records; `gui-matrix.json` records
**131.59 s** of summed driver execution, excluding gaps between invocations. Eight camera
contracts and 43 visible-HUD checks pass. All **24 selected fresh original captures were
inspected**. Partial-merge denial retains a filled triangle, bounds and an explicit error
cue; adoption denial retains four filled triangles. Eager/aircraft geometry and edit
handles remain visible on both backends. Text clipping/overlap and offscreen edit geometry
remain S5 work; BoT edit captures do not independently qualify surface quality. Private-
Xvfb System GL uses CPU llvmpipe and supplies no native GPU qualification.

`current-hashes.json` and `runtime.json` retain **57 runtime artifacts**, **84 resolved
shared libraries** and **140 physical files**, with binary identities and loader resolution.
qged remains `3eca3e98...`; libBObol is `ee379836...` and libged `b2b2e344...`.
Obol remains `5d28342f...` and OSMesa `097aac27...`; their native suites were not rerun
for unchanged binaries. The formal baseline and user's proposal are unchanged. Full TLC,
shared-stack sanitizers, coherent dependency staging and broad resource/native-host
qualification remain release obligations.

Next audit the four scalar shape-state setters, then root replacement and repository
sharing/clearing, view-variant invalidation and object rename. The adjacent GED local
mesh/vlist factories return borrowed targets after append/move callbacks and now have
explicit caller-acceptance work. Source/service and presentation ownership follow.
The [stage estimates](#september-9-stage-completion-estimates) remain unchanged.
**GROUP-007 closes shape movement/removal; S1--S6 remain incomplete.**

## September 9: complete source removal

`CXX-GROUP-006` repairs `removeDatabaseSource`, `removeDatabaseSourceInstance` and
`clearDatabaseSources`. The preceding library exposes partial scene state in all six
source-removal probe cases. Single-source removal also drops still-reachable lookups
and fails to retire the complete descendant ownership. Bulk clear notifies between
parent edits. A separate collision reproducer shows path removal deleting a different
source because it retargets the request through a shared instance key.

The existing `ChildPublication` now prepares a selected source edge or all source edges
to clear. Membership, child paths, descendant indexes, repository ownership and one
structural/frame effect commit before observers. Shared sources keep surviving edges
and lookups; bulk clear preserves ordinary nodes and auxiliary source subtrees. Path
removal keeps the selected source identity. Single-child edits use the native prepared
removal's exact position so duplicate node references update the correct paths. Later
callback changes remain current. The old recursive clear helper and incremental source
unindex writers are deleted; the scene controller is **65 net lines smaller**.

Evidence is retained in `.build/obol-qualification/20260909-source-removal`. Both new
production checks reject the preceding shared library, with loader resolution recorded.
Forty-eight allocation scenarios cover **4,948 injected positions**, all preserving the
preceding scene, plus forty-eight unconstrained reference runs. They include shared
cached geometry, duplicate child edges, shared parent groups, nested/auxiliary sources,
node/path observation, batches and quiet fields. Throwing observers preserve complete
publication and remaining notification. Ninety-two callback cases cover observation,
drain, reinsertion/replacement, parent/root removal, quiet policy and aliased requests.
Separate cases cover colliding path/instance targets, invalid/missing requests, explicit
auxiliary removal and original-target lifetime without an external source reference.

All five clients rebuild in **46.76 s**. All **22 focused checks** pass in
**139.37 s**, including **93 direct publication checks and 21 scenarios**, each
run once. All preceding 91 direct checks and 30 GED appearance cases remain present.
The scene controller compiles independently with production flags. Thirteen rebuilt
scene/GED probes pass; the original six source-removal cases now expose complete scene
state and the collision case removes only the selected source.

All **24 graphical rows** pass in **131.79 s**: two eager LoD-off/on rows, ten delivery
and rejection rows, eight cold/warm Generic Twin rows and four primitive-edit rows with
207 events each. Eight camera contracts and 43 visible-HUD checks pass. All **24 selected
original captures were inspected**. Both partial-merge-denied captures show a filled
triangle, retained bounds and an explicit error cue; adoption denial retains four filled
triangles. Eager and aircraft geometry and primitive-edit handles remain visible on both
backends. Text clipping/overlap and offscreen edit geometry remain S5 work; BoT edit
captures do not independently qualify surface quality. Private-Xvfb System GL uses CPU
llvmpipe and supplies no native GPU qualification.

`current-hashes.json` and `runtime.json` retain **55 runtime artifacts**, **84 resolved
shared libraries** and **138 physical files**, with binary identities and loader resolution.
qged remains `3eca3e98...`; libBObol is `5ab7573f...` and libged `8d9d462c...`.
Obol remains `5d28342f...` and OSMesa `097aac27...`; their native suites were not rerun
for unchanged binaries. The formal baseline and user's proposal are unchanged. Full TLC,
shared-stack sanitizers, coherent dependency staging and broad resource/native-host
qualification remain release obligations.

Next audit shape movement/removal, root replacement and repository sharing/clearing,
view-variant invalidation and object rename before the remaining GED/source/service/
presentation work. Shared-source cache ownership across controllers remains explicit
acceptance work. The [stage estimates](#september-9-stage-completion-estimates) remain
unchanged. **GROUP-006 closes source removal and clearing; S1--S6 remain incomplete.**

## September 9: complete group removal

`CXX-GROUP-005` repairs `eraseGroupSubpath`, `removeGroup` and `clearGroup`. The preceding
library exposes partial scene state to observers in all six removal-probe cases, covering
all three entry points with shared and unshared subtrees. Each writer changed children
before advancing revisions and separately released repository entries.

The existing `ChildPublication` now prepares a whole-group clear as well as a single
child edit. It retains all affected children, commits membership, path changes, indexes,
repository ownership and one structural/frame effect before observers, and preserves
sources still reachable through another edge. Old cache values remain retained through
notification. The two hierarchy-path entry points share one resolver and publisher;
their old live helper and the separate clear/remove/index/revision sequences are removed.
The scene controller is **17 net lines smaller**. No new publication framework is added.

Evidence is retained in `.build/obol-qualification/20260909-group-removal`. Both new
production checks reject the preceding shared library, with loader resolution retained.
Twenty-four allocation scenarios cover **2,304 injected positions**, all preserving the
preceding scene, plus twenty-four unconstrained reference runs. They include shared cached
geometry, duplicate child edges, nested/auxiliary sources, node/path observation, batches
and quiet fields. Throwing observers preserve complete publication and remaining notification.
Forty-eight callback cases cover no-ops, refill, movement, parent removal/replacement,
root replacement, quiet policy and aliased request strings. Root clearing and invalid,
empty and missing paths have separate cases. Hierarchy removal retains component parsing;
indexed clearing retains its canonical-key lookup contract.

Early test runs exposed fixture assumptions: field snapshots retain geometry references,
cross-fixture shared geometry must compare contents while per-fixture identity stays fixed,
and indexed group lookup does not accept a trailing slash on a non-root key. Two shadowed
test variables also failed compilation. Those issues were corrected before final
qualification; failed runs remain retained, including an early run named
`removal-transaction-qualified` which did not pass and supplies no qualification.

All five clients rebuild in **42.92 s**. All **22 focused checks** pass in **124.57 s**,
including **91 direct publication checks and 21 scenarios**, each run once. All preceding
89 direct checks and 30 GED appearance cases remain present. The changed scene-controller
translation unit compiles independently with production flags. Twelve rebuilt scene/GED
probes pass; the original removal probe now observes complete state in all six cases.

All **24 graphical rows** pass in **134.53 s**: two eager LoD-off/on rows, ten delivery
and rejection rows, eight cold/warm Generic Twin rows and four primitive-edit rows with
207 events each. Eight camera contracts and 43 visible-HUD checks pass. All **24 selected
original captures were inspected**. Both partial-merge-denied captures retain a filled
triangle, bounds and an explicit failure cue; adoption denial retains four triangles.
Eager and aircraft geometry and edit handles remain visible on both backends. HUD clipping,
annotation overlap and offscreen edit geometry remain S5 work; BoT edit captures do not
independently qualify surface quality. Private-Xvfb System GL uses CPU llvmpipe.

`current-hashes.json` and `runtime.json` retain **54 runtime artifacts**, **84 resolved
shared libraries** and **137 physical files**, with final binary identities and loader resolution.
qged is `3eca3e98...`, libBObol `415c7ee2...` and libged `281853c0...`.
Obol remains `5d28342f...` and OSMesa `097aac27...`; native dependency suites were not rerun
for unchanged binaries. The formal baseline and user's proposal are unchanged. Full TLC,
shared-stack sanitizers, coherent dependency staging and broad resource/native-host
qualification remain release obligations.

The inventory now explicitly includes source-instance removal, source-wide clearing and
repository sharing/clearing/view-variant invalidation/object rename. These existing scene
writers retain their own acceptance obligations. Next audit source removal/clearing, then
shape, root and repository ownership before the remaining GED/source/service/presentation
work. The [stage estimates](#september-9-stage-completion-estimates) remain unchanged.
**GROUP-005 closes recursive group removal and clearing; S1--S6 remain incomplete.**

## September 9: complete child membership

`CXX-GROUP-004` repairs `appendChildToGroup` and `removeChildFromGroup`. The preceding
library exposes partial membership/index/revision state in all six source, group and
ordinary-node probe cases. Direct-source append/removal also leaves incomplete lookups.
The new child publisher composes existing hierarchy and index preparation with a prepared
repository update, then commits one complete structural/frame effect before observers.
Nested and auxiliary sources are included. Shared nodes remain indexed and owned while
another edge reaches them; colliding path, instance and group keys resolve against the
completed graph. Ordinary new insertion prepares the incoming subtree, while removals and
collisions inspect reachability without making another full scene copy.

Repository seeding previously rewrote published geometry identity and could partially
change cache ownership on allocation failure. Seeding and release now share a private
prepared update for affected object counts and cache entries. Commit transfers map nodes;
retired values remain owned until after the enclosing observers. Cache-reference insertion
returns its acquired reference on failure. The child publisher retains original targets
and performs no later writes after callbacks. The old live child/index sequence, direct
source-index helper and repository object-release writer are removed.

Evidence is retained in `.build/obol-qualification/20260909-group-membership`. All three
new production checks reject the preceding shared library, with loader resolution recorded.
Forty child allocation scenarios cover **3,708 injected positions**, all preserving the
preceding scene, plus forty unconstrained reference runs. Explicit throwing-observer cases
preserve complete publication. Eight wire/mesh repository scenarios cover **92 injected
positions** plus eight reference runs. Eighty callback cases cover no-ops, inverse changes,
movement, parent removal/replacement, root replacement, quiet policy and aliased arguments.
Additional cases cover shared-parent lookup fallback, colliding keys and two child cycles;
the preceding six source-movement cycle cases remain passing.

All five clients rebuild in **36.16 s**. All **22 focused checks** pass in **122.95 s**,
including **89 direct publication checks and 21 scenarios**, each run once. All preceding
86 direct checks and the 30 GED appearance cases remain present. All four changed production
translation units compile independently with production flags. Eleven rebuilt scene/GED
probes pass; the original membership probe now sees complete state in all six cases.
Early probe/test compilation mistakes and one stale-binary selector failure were corrected
before final qualification; their failed runs remain retained.

All **24 graphical rows** pass in **131.40 s**: two eager LoD-off/on rows, ten delivery
and rejection rows, eight cold/warm Generic Twin rows and four primitive-edit rows with
207 events each. Eight camera contracts and 42 visible-HUD checks pass. All **24 selected
original captures were inspected**. Both partial-merge-denied captures retain a filled
triangle, bounds and an explicit failure cue; adoption denial retains the preceding four
triangles. Eager and aircraft geometry and edit handles remain visible on both backends.
HUD clipping, annotation overlap and offscreen edit geometry remain S5 work; BoT edit
captures do not independently qualify surface quality. Private-Xvfb System GL uses CPU llvmpipe.

`current-hashes.json` and `runtime.json` identify **53 runtime artifacts**, **84 resolved
shared libraries** and **136 retained physical files**, with final loader resolution.
qged is `3eca3e98...`, libBObol `646bf2a4...` and libged `d8c4e2df...`.
Obol remains `5d28342f...` and OSMesa `097aac27...`; native dependency suites were not rerun
for unchanged binaries. The formal baseline and user's proposal are unchanged. Full TLC,
shared-stack sanitizers, coherent staging and broad resource/native-host qualification
remain release obligations.

Next audit `eraseGroupSubpath`, `removeGroup` and `clearGroup`, then shape and root operations
and the remaining GED/source/service/presentation owners. The
[stage estimates](#september-9-stage-completion-estimates) remain unchanged.
**GROUP-004 closes child membership and its prepared repository update; S1--S6 remain incomplete.**

## September 9: complete group rename

`CXX-GROUP-003` repairs group rename publication. The preceding library completes the
rename correctly, but both root callbacks observe incomplete paths, indexes and revisions.
Rename now prepares descendant paths through the existing `GroupPublication` owner and
prepares the group index before changing the name. All scene effects commit before
notification. Original participants stay retained, all notification policies are restored
before the first callback, and later callback changes survive. The separate recursive
path writer and postnotification index invalidation are removed. Shared subtrees retain
the preceding traversal's final canonical path, with each group prepared once; traversal
through ordinary scene groups remains supported. Rename targets retain the existing
BRL scene-group lookup contract.

Allocation injection also reproduced an underlying Obol defect: `SoBase::setName` removed
the old registration before allocating its replacement. The two private name maps now
use standard containers with owned name lists. Normalization, destination capacity and
map entries prepare before removal; name clearing and object destruction retire empty
registration buckets. Same-name registration still selects the last registered object.
No public API or separate name-publication framework is added.

Evidence is retained in `.build/obol-qualification/20260909-group-rename`. The three new
production checks reject the preceding shared libraries, with loader resolution retained.
All **2,040 injected positions** in the twelve rename allocation sweeps
preserve the preceding scene. Twelve unconstrained runs establish the bounds; explicit
throwing-observer cases also preserve complete publication. Six native-name scenarios sweep **10 positions**
plus their reference runs, covering first registration, rename, shared/same names,
normalization and clearing. Twenty-four callback cases cover later/descendant renames,
removal/replacement, root replacement, quiet policy, no-ops and aliased arguments.
Source fields, geometry, realization stamps and source lookups remain unchanged.

The five-client build passes in **64.34 s**. All **22 focused checks** pass in
**113.44 s**, including 86 direct publication checks and 21 scenarios, each run once.
All preceding 83 direct checks remain, as do the 30 GED redraw appearance cases.
The scene controller and both changed native implementation units compile independently
with production flags. Ten rebuilt group/source/GED probes pass; the original rename
probe now observes complete state in both callbacks and one frame increment.
Native Obol runs pass **1,576 unit tests** (one existing profiler concurrency skip) and
**63 integration tests**, including name ordering, normalization and lifetime coverage.
The earlier name-test failure from a stale exception flag and the unsupported plain-group
rename fixture were corrected before final qualification; their failed runs remain retained.

All **24 graphical rows** pass in **131.57 s**: two eager LoD-off/on rows, ten delivery
and rejection rows, eight cold/warm Generic Twin rows and four primitive-edit rows with
207 events each. Eight Generic camera contracts pass, with 44 visible-HUD checks.
Twenty-four original captures were inspected. Eager geometry, intentional rejection cues,
aircraft geometry and edit handles remain present on both backends. The injected partial
merge failure shows a retained triangle plus bounds on OSMesa and nested bounds on System
GL; both retain two compact entries, exact source bounds, an explicit failure cue and no
pending work. This row qualifies those failure outcomes. Small-pane text
clipping/overlap and offscreen edit geometry remain S5 work; the BoT edit captures do not
independently qualify surface quality. Private-Xvfb System GL uses CPU llvmpipe.

`current-hashes.json` and `runtime.json` identify **52 runtime artifacts**, **84 resolved
shared libraries** and **135 retained physical files**. qged is `3eca3e98...`,
libBObol `14de95dd...`, libged `d95c0228...` and Obol `5d28342f...`.
The explicit Obol install prefix targets `.build`; OSMesa remains `097aac27...` without
restaging. The formal baseline and user's proposal are unchanged. Full TLC and shared-stack
sanitizers were not run; coherent staging, complete formal evidence and the existing
large/cold/resource/native-platform qualification remain release obligations.

The remaining group append/removal/root and GED writers retain their inventory rows.
The [stage estimates](#september-9-stage-completion-estimates) remain unchanged.
**GROUP-003 closes this rename/name-registration boundary; S1--S6 remain incomplete.**

## September 9: complete source movement

`CXX-GROUP-002` repairs the standalone source move boundary. Previously, creating a
destination could notify before the source was retained, then live remove/add/reindex
operations exposed partial membership and stale scene revisions. A new destination advanced
the frame twice. Movement now uses an explicit membership-only constructor on the existing
`SourcePublication` owner. It shares hierarchy/index preparation and commit with composed
source publication, while preserving source configuration, geometry, realization identity
and existing lookup/fallback semantics. The separate live mutation sequence is removed.

The prepared-hierarchy cycle concern is also reproduced and repaired. The preceding
composed publisher rejected an existing descendant but accepted a new destination beneath
it. Reachability now reads prepared child orders where present and live children otherwise.
Both publication APIs reject existing descendants and one/two-level new suffixes without
changing the scene or its revisions.

Evidence is retained in `.build/obol-qualification/20260909-source-movement`.
Three new production regression groups cover six cycle cases, existing/new/root destinations,
the path wrapper, no-ops, batches, quiet/throwing observers, aliased requests and callback
edits, moves, removal/replacement, group rename and root replacement. All three reject the
preceding shared library, with loader resolution retained. Twelve allocation sweeps exercise
**892 injected positions**, all preserving the preceding scene, plus twelve unconstrained
runs. Explicit throwing-observer cases preserve a complete commit. Source field/geometry
snapshots and realization stamps remain unchanged by membership-only movement.

The five-client build passes in **64.23 s**. All **22 focused checks** pass in **112.52 s**,
including 83 direct publication checks and 21 scenarios, each run once, and 30 GED redraw
appearance/group cases. All preceding 80 direct checks remain. The changed scene controller
also compiles independently with production warning flags. Nine rebuilt group/source/GED
probes pass, including the original source-move reproducer: new and existing destinations
now expose complete state and exactly one frame increment.

All **24 graphical rows** pass in **132.06 s**: two eager LoD-off/on rows, ten delivery
and rejection rows, eight cold/warm Generic Twin rows and four primitive-edit rows with
207 events each. All eight Generic camera contracts pass, with 43 visible-HUD checks.
Twenty-four original captures were inspected, preserving the eager triangles, intentional
rejection cues, aircraft geometry and edit handles on both backends. Small-pane text
overlap/clipping and offscreen edit geometry remain S5 work; the BoT edit captures do not
independently qualify surface quality. Private-Xvfb System GL uses CPU llvmpipe.

`current-hashes.json` and `runtime.json` identify **52 runtime artifacts**, **84 resolved
shared libraries** and **135 retained physical files**. qged remains `3eca3e98...`;
libBObol is `2e47aca8...` and libged is `c7fca5b8...`. Obol (`194d89fe...`), OSMesa
(`097aac27...`), the formal baseline and the user's proposal are unchanged. No dependency
restaging occurred. Native dependency suites and full TLC were not rerun for this C++
repair. Shared-stack sanitizers, coherent staging, full formal evidence and the existing
large/cold/resource/native-platform qualification remain release obligations.

`CXX-GROUP-003` remains open. The retained `group_rename_probe` finishes the rename
correctly, but both observer calls see incomplete path/index/revision state. This is the
next repair in the finite group writer inventory. Passing a callback rename inside the
move test does not qualify the rename operation's own notification boundary.

The [stage estimates](#september-9-stage-completion-estimates) remain unchanged.
**GROUP-002 closes this membership/cycle boundary; S1--S6 remain incomplete.**

## September 9: complete group creation

`CXX-GROUP-001` repairs `ensureGroup`: observers previously ran after each individual
append, before its index entry and before the enclosing structural/frame revision.
The minimal probe reproduces partial observation with both a new hierarchy and an existing
prefix. Both cases now observe the complete hierarchy and revisions. After notification,
the API resolves the canonical requested path again, accepting callback removal, rename,
replacement or root replacement without returning a detached prepared target.

The existing source publisher's hierarchy preparation is shared through the private
`HierarchyPublication` owner. It retains affected nodes, prepares child-list replacements
and index entries, commits them without notification, then delivers observers after scene
effects. This removes the separate live creation loop and its unused `indexSceneGroup`
writer. No source configuration, geometry, representation policy or renderer behavior is
added. The public header documents the current-target return contract.

Evidence is retained under `.build/obol-qualification/20260909-group-creation`.
Two new production regression groups cover eight allocation scenarios, new/existing
prefixes, batches, quiet/throwing observers, recursive no-ops/extensions, removal,
replacement, rename, root replacement and aliased request strings. All **844 injected
allocation positions** preserve either the preceding scene (**812**) or a complete commit
(**32**); eight unconstrained runs establish the sweep bounds. Postcommit allocation
failures include final path lookup, so committed outcomes are not all notification faults.
The same two production checks reject the preceding shared libBObol; retained loader
resolution confirms that baseline. Standalone GED-style probes instead link archives.

The five-client build passes in 34.86 s. All **22 focused checks** pass in **114.12 s**,
including 80 direct publication checks and 21 scenarios, each run once, and 30 GED redraw
appearance/group callback cases. All preceding 78 direct checks remain. The changed
`scene_controller.cpp` also compiles independently with production warning flags. Seven
rebuilt prior source/GED probes pass, preserving the earlier composed publication repairs.

All **24 graphical rows** pass in **131.92 s**: two eager LoD-off/on rows, ten delivery
and rejection rows, eight cold/warm Generic Twin rows and four primitive-edit rows with
207 events each. Generic camera contracts pass on all eight rows, with 43 visible-HUD
checks. Twenty-four original captures were inspected. They retain the eager triangles,
intentional rejection cues, aircraft geometry and edit handles on both backends.
Text overlap/clipping and offscreen edit geometry remain S5 work; the BoT edit capture
does not independently qualify surface quality. Private-Xvfb System GL uses CPU llvmpipe.
The report/image checks do not retire the retained Lucy, desktop stall, resource or native
platform failures.

`current-hashes.json` and `runtime.json` identify 51 runtime artifacts, 84 resolved shared
libraries and 134 retained physical files. qged remains `3eca3e98...`; libBObol is `8756c6d3...` and libged is
`560341d6...`. Obol (`194d89fe...`), OSMesa (`097aac27...`), the formal baseline and the
user's proposal are unchanged. No dependency restaging occurred. Native dependency suites
and full TLC were not rerun for this C++ hierarchy change; shared-stack sanitizers, coherent
dependency staging and final formal evidence remain release obligations.

The first runtime-retention tool session terminated with status 143 before its wrapper
wrote a completion record; the cause is not established. A host check found no surviving
job. Retention resumed successfully with streaming hashes cached by physical file; checkpoint
verification checks the existing content-addressed copies. This recording interruption is separate from the
completed build, focused and graphical checks.

`CXX-GROUP-002` remains open: the new `group_move_probe` sees partial source movement with
both destination conditions. A new destination produces four callbacks and two frame
increments; an existing one produces two callbacks and one increment. Both finish with
the source moved, but their intermediate observations violate the publication boundary.
This retained failure is the next repair, followed by the enumerated group/GED writers.
Cycle rejection through newly prepared descendant paths is a separate code-review concern
to reproduce during that repair; no runtime result for that case is claimed here.

The [stage completion estimates](#september-9-stage-completion-estimates) remain unchanged.
**GROUP-001 closes creation only; GROUP-002 and S1--S6 remain incomplete.**

## September 9: source and group publication composed

`CXX-GED-004` repairs the enclosing source/group metadata transaction. GED previously
published the source, then wrote group intent and display from an earlier snapshot. Its
source callback could set group line width to 7, only to have the later group sync restore
3. The scene now prepares group metadata with the source, membership, indexes and revisions,
commits the full update, and restores notification policy before observers. The later GED
group-sync/regroup helper is deleted, removing **50 net lines** from `draw_obol.cpp`.

The standalone group intent/display setters use the same prepared group owner and publish
their frame effect before notification. Existing scalar copy, commit and notification
helpers are shared through `scalar_publication_private.h`; geometry and source-specific
ownership remain in the source implementation. The existing source-only API remains.
Group-only updates preserve source state and revisions; invalid metadata requests publish
nothing. No separate notification framework or callback prohibition is introduced.

Evidence is retained in `.build/obol-qualification/20260909-source-group-publication`.
The original six-mode probe fails on the preceding GED implementation and passes after
the repair. Thirty production GED cases include six group-edit callbacks; the same test
object linked against the preceding GED archive passes the first four source-appearance
cases, then rejects the new group case. Four additional production regression groups cover
both callback directions, new/existing groups, quiet fields, no-ops, group-only effects,
batches, removal/replacement, exceptions, aliased intent and invalid metadata targets.
Eleven allocation sweeps exercise **2,766 injected positions plus 11 unconstrained runs**:
2,722 preserve the preceding state and 44 leave a complete commit when notification fails.

All **22 focused checks pass in 115.16 seconds**, exercising all **78 direct publication
checks and 21 scenarios once**, plus the 30 GED redraw cases. All five clients rebuild;
`database_source.cpp`, `scene_controller.cpp` and `draw_obol.cpp` compile independently.
The existing source publication, six-mode redraw and GED metadata/material-target probes
remain green. An initial test build rejected an exact float comparison; the corrected
comparison passes. An invalid-target fixture was corrected to use the ordinary scene root,
since an ordinary named group has no indexed metadata path.

The standalone group probe rejects the retained preceding libBObol archive for both
intent and display: observers see partial metadata/frame state and lose their newer edit.
The repaired binary passes both. Its initial fixture omitted a scene root and exited
before exercising a setter; an explicit retained root corrected the fixture. Loader guards
also exposed that this probe links libBObol statically. The baseline comparison therefore
relinks the same object with the retained archive and records the exact link command.

All **24 fresh graphical rows pass in 126.53 seconds**: two eager draws, ten delivery/failure
rows, eight Generic Twin cold/warm shaded/wire rows and four primitive-edit rows. Editing
completes 207 events per row; Generic adds eight camera contracts and 41 HUD pixel checks.
All **24 selected original captures were inspected**. Eager geometry remains filled,
constrained delivery retains its valid prefix/preceding geometry and diagnostic, and
ordinary-model and editing geometry remains visible. HUD clipping, annotation overlap
and offscreen edit geometry remain S5 acceptance work.

The checkpoint retains **49 runtime artifacts, 84 distinct resolved libraries and
132 physical files**, with loader resolution and hashes. The final stack is
qged `3eca3e98...`, libBObol `e6ef2e46...`, libged `92a07fda...` and
libqtcad `372682d2...`. Obol, the restored OSMesa, the formal baseline and the user's proposal
are unchanged. Native dependency suites were not rerun for unchanged binaries. This Debug
build has sanitizers off; private-Xvfb System GL uses CPU llvmpipe. Required native-host,
shared-stack sanitizer and broad resource/capability qualification remain open.

CMake configuration again restaged the external bundle's OSMesa `75dfd5f3...`. After the
build ended, the matching local OSMesa `097aac27...` was restored with an explicit `.build`
prefix and verified before tests. Coherent dependency staging remains S6 work.

The separate group hierarchy entry points, base promotion, deferred-proxy clearing,
retained-current marking, broad/narrow occurrence adoption and live deferred target/reference
acceptance still require their own audits. Source/service and presentation ownership follow.
**GED-004 qualifies the composed metadata boundary and two standalone metadata setters;
S1--S6 and their estimates remain unchanged.**

## September 9: redraw appearance preserved directly

`CXX-GED-003` removes **118 net lines** from `draw_obol.cpp`. Root redraw previously saved
matching source appearance, published default appearance, and restored the saved values.
Both field and node probes reproduced every accepted mode observing line width 1 instead
of the preceding 3, then losing a newer callback edit to 7. The redraw now supplies no
appearance override: the existing publisher obtains the exact current source's appearance
while preparing its candidate. The saved-record type, collection/restoration writers,
default settings and an inert root metadata record are removed.

Evidence is retained in `.build/obol-qualification/20260909-redraw-appearance`.
The before probes fail in all six modes. The final field probe has no temporary width
notification; the node probe observes 3 and retains its newer 7 in every mode. Twenty-four
production GED cases cover unchanged appearance and callback edit, replacement and removal
across modes 0--5, including current identity/revision and realized geometry. Linking the
same regression object against the retained preceding GED archive rejects the transient
appearance. The complete focused run passes **22 checks in 114.77 seconds**, including
all **74 direct publication checks and 21 scenarios once**, plus the 24 new GED cases.
All five clients rebuild; `draw_obol.cpp` compiles independently. The six-mode redraw,
source publication and GED metadata/material-target probes remain green.

The first probe driver received SIGTERM (tool exit 143) after completing the publication
probe. Host inspection found no surviving related children. Its logs and successful probe
were retained; only the unfinished probes were resumed, and their final checks completed.

All **24 fresh graphical rows pass in 128.38 seconds**: two eager draws, ten delivery/failure
rows, eight Generic Twin cold/warm shaded/wire rows and four primitive-edit rows. Editing
completes 207 events per row; Generic adds eight camera contracts and 39 HUD pixel checks.
All **24 selected original captures were inspected**. Eager geometry remains filled,
constrained delivery retains its valid prefix/preceding geometry and diagnostic, and
ordinary-model and editing geometry remains visible. HUD clipping, annotation overlap
and offscreen edit geometry remain S5 acceptance work.

The checkpoint retains **48 runtime artifacts, 84 distinct resolved libraries and
131 physical files**, with loader resolution and hashes. The final stack is
qged `3eca3e98...`, libBObol `d57d9a37...`, libged `8340979a...` and
libqtcad `372682d2...`. libBObol, Obol, OSMesa, the formal baseline and the user's proposal
are unchanged. Native dependency suites were not rerun for unchanged binaries. This Debug
build has sanitizers off; private-Xvfb System GL uses CPU llvmpipe. Required native-host,
shared-stack sanitizer and broad resource/capability qualification remain open.

The separate `group_probe` still fails in all six modes: a source callback sets its group's
line width to 7, then `ged_obol_sync_group_state` overwrites it with 3. Source appearance
and geometry remain correct. This is **GED-004**, the next source/group composition repair.
Base promotion, proxy clearing, retained-current marking, broad/narrow occurrence adoption
and live deferred target/reference acceptance still require their own audits.
**GED-003 is closed at the source-appearance boundary; S1--S6 and their estimates remain
unchanged.**

## September 9: stage completion estimates

These planning ranges estimate progress against the simplification guide's exit
criteria. They are engineering judgments, not measured fractions of effort, time
forecasts, or release acceptance. The stages differ in size and overlap, so their
percentages should not be averaged. None is complete.

| Stage | Estimated completion | Main remaining obligation |
|---|---:|---|
| S1: ownership and acceptance inventory | 65--80% | Finish writer classification and explicit acceptance coverage |
| S2: one complete production path | 80--90% | Consolidate demonstrated boundary repairs into stage closure |
| S3: physical ownership boundaries | 40--60% | Finish source/service, GED and presentation ownership and remove superseded writers |
| S4: resources, async work and lifecycle | 45--65% | Close stalls and qualify memory, cancellation, BREP and shared-cache behavior |
| S5: capability qualification | 30--50% | Large/cold scenes, visual failures, editing/shared clients and required native hosts |
| S6: release candidate qualification | 0--10% | Freeze and qualify the exact coherent stack against the required matrix |

Focused repairs support the method and reduce known defects. Broader resource
and platform qualification remains the largest source of uncertainty. Revise
these estimates when a gate's evidence changes materially, rather than after
each individual passing regression.

## September 9: composed scene source publication

`CXX-CONFIG-007` repairs the enclosing scene source publisher. Its revision callback
previously saw old line style and frame state, and a newer callback edit was overwritten.
The source now prepares a complete scalar candidate using existing normalization, owned
metadata, compact style and placement helpers. The scene prepares insertion, moves,
missing parent groups and indexes without changing the live graph. All effects commit
before notification, and prepared references retain callback targets. The old sequence
of separately notifying live setters is removed. Existing rename index preparation is
shared; immutable geometry remains shared.

Evidence is retained in `.build/obol-qualification/20260909-scene-source-publication`.
The original probe fails on the preceding stack and passes on the final stack, including
the later callback edit. New transaction and callback tests also reject the retained
preceding libBObol, with verified loader resolution. Six allocation sweeps cover 2,379
injected positions plus six unconstrained runs: 2,338 preserve the preceding scene and
41 complete the commit before notification fails. Three production regression groups
cover complete source/scene observation, quiet fields, batches, nested/auxiliary owners,
insertion, missing groups, shared paths, retargeting, moves, no-ops, cycle rejection and
callback style/removal/replacement/move. Actual compact renderer records preserve geometry
and local bounds while applying placement, line style, opacity and full-path materials.

All **22 focused checks pass in 112.59 seconds**, exercising all **74 direct publication
checks and 21 scenarios once**. Both changed production translation units compile
independently. All five application clients rebuild. The original publication probe,
six-mode GED redraw probe and GED metadata/material-target probe pass on final libraries.
Early build failures were corrected shadowing/unused-helper and test API mistakes; an
initial edge invocation used the preceding test binary after a failed test build and
reported an unknown selector. The corrected build and all final checks pass.

All **24 fresh graphical rows pass in 131.56 seconds**: two eager draws, ten delivery/failure
rows, eight Generic Twin cold/warm shaded/wire rows and four primitive-edit rows. Editing
completes 207 events per row; Generic adds eight camera contracts and 39 HUD pixel checks.
All **24 selected original captures were inspected**. Eager drawings remain filled,
constrained delivery retains its valid prefix/preceding geometry and diagnostic, and
ordinary-model and editing geometry remains visible. Annotation overlap, offscreen edit
geometry and clipped HUD text remain S5 acceptance work.

The checkpoint retains **45 runtime artifacts, 84 distinct resolved libraries and
128 physical files**, with loader resolution and hashes. The final stack is
qged `3eca3e98...`, libBObol `d57d9a37...`, libged `ea218288...` and
libqtcad `372682d2...`. Obol, OSMesa, the formal baseline and the user's proposal are unchanged.
Native dependency suites were not rerun for unchanged binaries. This Debug build has
sanitizers off; private-Xvfb System GL uses CPU llvmpipe. Required native-host,
shared-stack sanitizer and broad resource/capability qualification remain open.
Continue enclosing GED replace/redraw, resident-occurrence adoption and live deferred
target/reference acceptance, then source/service and presentation ownership.
**CONFIG-007 is closed at this boundary; S1--S6 and their estimates remain unchanged.**

## September 9: unreachable GED leaf redraw removed

`CXX-GED-002` removes **438 net lines** from `draw_obol.cpp`. The previously queued
metadata/line-style adapter was unreachable: the entry gate accepts exactly modes 0--5,
and all six return through the root-source branch before the leaf walk. Its sole caller
also always selects synchronous redraw, mixed-mode preservation and display preservation,
with no result output. The dead leaf traversal, instance-key generation, separate metadata
and style adapter, deferred branch, result counters, redundant mode predicate and unused
options are removed. `ged_obol_redraw_source_paths` now takes the existing path vector and
mode directly. Actual deferred drawing retains its separate live adapter.

Evidence is retained in `.build/obol-qualification/20260909-ged-root-draw-audit`.
The source audit records both predicates and the complete caller chain. The preceding
executable's debugger trace enters the root helper once for each mode 0--5 and never enters
the leaf callback. The actual six-mode redraw probe passes before and after cleanup,
preserving line width, mode and realized wire/mesh geometry. The final debugger trace
confirms all six entries in the simplified helper; obsolete symbols are absent from the
rebuilt archive.

Initial probe setup needed the private endpoint declaration, a scene root, the actual
source-update request rather than invalidation, and an explicitly stale source after
clearing geometry. Its initial geometry oracle counted legacy nodes but omitted compact
renderer geometry. These tool/fixture failures remain alongside the corrected passing
probe; they are not implementation-build failures.

All five clients rebuild, the edited translation unit compiles independently, and
**all 22 focused checks pass in 106.50 seconds**. Existing publication groups still
execute all 71 direct checks and 21 installation/merge scenarios once. The preceding GED
metadata/material-refresh probe passes. **All 24 fresh graphical rows pass**: two eager,
ten delivery/failure, eight Generic Twin and four primitive-edit rows, with 207 edit events
per row, eight camera contracts and 41 HUD pixel checks. All 24 selected original
captures were inspected. Annotation overlap, offscreen edit geometry and HUD clipping
remain S5 work.

The checkpoint retains **45 runtime artifacts, 84 resolved libraries and 128
physical files**. Current libged is `1e0249c5...`; libBObol, Obol,
OSMesa, the formal baseline and the user's proposal are unchanged. Native dependency
suites were not rerun for unchanged binaries. This Debug stack has sanitizers off and
private-Xvfb System GL uses CPU llvmpipe. All jobs are terminal; the final host check is empty.

The audit also reproduces **open `CXX-CONFIG-007`** in the live scene source publisher:
its source-revision callback sees old line style and an old scene frame, and a newer
callback style edit is overwritten by the remaining setter. The retained
`publication_probe` intentionally fails before and after cleanup. The next repair must
compose the existing source publish state and its scene/index effects before observers,
including insertion and movement; then continue enclosing GED target/reference acceptance.
This is a remaining defect, not closure of the whole draw transaction. **S1--S6 remain open.**

## September 9: GED metadata and material-refresh effects

`CXX-GED-001` closes metadata's missing scene effects and material refresh's acceptance
of later edited/replaced targets. The scene controller normalizes the existing typed
metadata record and contributes frame effects to the aggregate and compact publishers.
GED's duplicate converter and key-by-key refresh loop are removed. Targeted and full-scene
refresh share the existing database sweep, retaining original owners and checking material,
CAD and node revisions plus scene membership before each source. Targeted refresh keeps
per-source database fallback; an explicit override applies to all targets. A single cache
is reused while the database is unchanged and replaced when switching databases. Each
changed source now advances the frame before observers; the old first-source-only
suppression is removed.
Explicit mutation batches retain coalescing, and compact stamp-only refresh remains quiet
with respect to frame effects. Unchanged valid metadata reports successful application to
GED without advancing the frame.

Evidence is under `.build/obol-qualification/20260909-ged-metadata-effects`. The retained
preceding probe reproduces missing frame effects for both metadata representations and
overwrites of a later material revision and replacement owner. The current probe passes
all four cases. The new GED regression rejects verified preceding static `libged.a` and
`libBObol.a`; current GED passes ordinary/throwing observers and unchanged application too.
Review added a mixed-database regression that rejected the first implementation's use of
the first source's database for every target. Per-source fallback, explicit override and
skipping a missing database without suppressing later valid targets now pass. Full-scene
refresh and the explicit-database free helper retain their established database semantics.
Initial probe compilation errors concerned reference-wrapper construction and indentation.
The first full-scene sweep fixture reused a nested-source identity; distinct owner keys
correct the fixture without changing the production acceptance check. These failures are
retained separately from passing implementation and qualification builds. The combined
publication test exceeded its existing 180-second CTest timeout. Source and state sweeps
now run as two explicit groups, each retaining that timeout; no case or assertion was
removed. Execution markers verify the full partition, and individual selectors remain
available. The focused matrix consequently has 22 rows.

New scene/batch aggregate and compact metadata coverage comprises **8
scenarios, 386 C++ allocation positions: 372 preserved
and 14 committed**. Independent metadata/RGB expectations and actual compact
renderer readback verify scope, retained geometry/style, quiet flags, no-ops and retries.
The separate targeted/full-scene sweeps cover **4 scenarios and
866 positions: 378 rejected, 406
completed prefixes and 82 with both primary sources committed**. Observers
check each source/frame pair. Later revisions, quiet metadata edits without a material
revision change, replacement/removal, reentrant current edits, throwing observers and
explicit batches retain direct coverage. These tests do not establish large-workload
memory, cache contention or latency bounds.

The five clients rebuild and **all 22 focused checks pass in 106.34 seconds**. Four
production translation units compile independently; prior controller, metadata, state,
remaining-scene, database-material, sparse, retained and region probes pass. **All 24 fresh
graphical rows pass**: two eager draws, ten delivery/failure rows, eight Generic Twin
cold/warm shaded/wire rows and four primitive-edit rows. Editing completes 207 events per
row; Generic adds eight camera contracts and 42 HUD pixel checks. All
**24 selected originals** were inspected. Geometry and explicit constrained outcomes remain
visible; annotation overlap, offscreen edit geometry and HUD clipping remain S5 work.

The checkpoint retains **51 runtime artifacts, 84 resolved libraries and 134
physical files**. qged remains `f14ea62d...`; current libBObol is
`8f7035d1...` and libged `e6006976...`. Obol,
OSMesa, the formal baseline and the user's proposal are unchanged. Native dependency suites
were not rerun for unchanged binaries; their preceding sealed results identify those
binaries. CMake restaged the known mismatched OSMesa bundle while registering the split
tests. Host-process checks confirmed no active test/GUI before the matching local OSMesa
was restored with an explicit workspace install prefix; hashes were verified before any
subsequent test or GUI. Coherent dependency staging remains S6 work. This Debug build has
sanitizers off, and private-Xvfb System GL uses CPU llvmpipe.
All jobs are terminal and the final host-process check is empty.

The direct-leaf adapter now pins its owner and copies its key across metadata callbacks,
skipping a replaced owner. Its following line-style step and the enclosing direct/deferred
draw transaction still require combined acceptance and reference-lifetime qualification.
Continue that audit, then source/service and presentation ownership. Retained stalls,
Lucy quality, measured resource limits and required native-host qualification remain open.
**S1--S6 remain incomplete; maturity estimates are unchanged by this individual repair.**

## September 9: evaluated-region source and scene publication

`CXX-REGION-001` repairs GED's partial evaluated-region updates: individual wire callbacks
could observe incomplete state, only the first mesh changed, and the scene frame effect
was absent. The source now prepares just the changed marker fields on all matching owned
shapes. It preserves unrelated region/material fields, geometry, producer stamps, auxiliary
shapes and nested owners. Shared path matching serves both representations; compact
sources retain the existing sparse publisher. Compiler/batch effects and each changed
source's scene revision commit before observers. The bulk operation preserves completed
prefixes and skips targets superseded by callback edits or replacement, including quiet
edits. Two unused private GED helpers are removed.

Evidence is retained under `.build/obol-qualification/20260909-evaluated-region-publication`.
Three production groups cover **8 source/publication scenarios and 144 C++ allocation
positions: 132 preserved and 12 committed**. The separate two-source sweeps cover **70
positions: 38 rejected, 30 completed prefixes and 2 complete operations**. Independent
expected fields and geometry checks cover ordinary/quiet/throwing observers, retry/no-op,
path boundaries and instance
suffixes, owner replacement, reentrant later edits and explicit mutation batches.
Noncompact picking checks semantic readback through source traversal; compiler retirement
and rebuild are checked separately. Compact checks inspect actual renderer records.
These allocation probes do not establish large-population memory or latency bounds.

The retained preceding production probe reports incomplete shapes, a partial observer and
no frame advance; the current probe reports complete shapes and a committed frame. The
new real GED regression likewise rejects preceding static `libged.a`/`libBObol.a` and
passes current libraries. Archive paths and hashes are retained; private GED entries are
not exported by the shared library. Early probe-link and test-fixture failures remain
recorded separately: companion shapes prevented compilation, noncompact compiled state
is not an attached child, and a NULL database requests removal rather than replacement.
The final fixtures exercise these distinctions without weakening production assertions.

All **21 focused checks pass in 177.11 seconds**. The five application clients rebuild;
four affected production translation units compile independently. Existing controller,
metadata, state, remaining-scene, database-material, sparse and retained probes pass.
All **24 fresh graphical rows** pass: two eager draws, ten delivery/failure rows, eight
Generic Twin cold/warm shaded/wire rows and four primitive-edit rows. Editing completes
207 events per row; Generic adds eight camera contracts and 42 HUD pixel checks.
All **24 selected original captures** were inspected. Filled eager drawings
and constrained delivery outcomes remain visible. Annotation overlap, offscreen edit
geometry and clipped HUD text remain S5 acceptance work.

The checkpoint retains **50 runtime artifacts, 84 distinct resolved libraries
and 133 physical files**, reusing unchanged retained content by hash. Current qged is
`f14ea62d...`, libBObol `ee4f84dc...` and libged `3168d64c...`; Obol, OSMesa, the formal
baseline and the user's proposal are unchanged. Native dependency suites were not rerun for unchanged
binaries; their preceding sealed results remain evidence for those identities. This Debug
build has sanitizers off. System GL uses private-Xvfb CPU llvmpipe; required native-host
and shared-stack sanitizer qualification remain open. Jobs are terminal and the final
host-process check is empty.

Continue GED draw-metadata application's enclosing scene effects and material-refresh
target/reentrant-owner acceptance, then source/service and presentation ownership.
Retained stalls, Lucy quality and the broader resource/capability matrix remain open.
**S1--S6 remain incomplete; the stage estimates are unchanged by this individual repair.**

## September 9: retained presentation and selection publication

`CXX-SPARSE-002` extends the existing prepared sparse publisher to nine retained operations:
visibility/highlight/transparency rules, highlight clearing, visibility frontier/override
replacement and clearing, selected-path replacement and selection deltas. Intent, affected
records, memberships, renderer state and compiler/journal effects commit before observers.
The three separate live reapply functions are removed. Shared matching/style helpers serve
current records and later arrivals; geometry, indexes and producer stamps retain ownership.
Meaningful retained-intent changes notify even without current matching occurrences.
Removing a selected child path now preserves coverage from a retained parent, and delta
matching agrees with replacement and streamed arrivals. Unrelated immediate selection survives.

The preceding library exposes partial publication in eight operation probes and incorrect
overlap removal in the ninth. The final probe passes all nine. Three regression groups cover
**36 scenarios and 342 C++ allocation positions: 282 preserved and 60 committed**. Independent
expected results cover current and future records, empty sources, quiet and throwing
observers, independently enabled selected-field notification, reentrant later rules,
authored appearance, exact visibility journals, normalized object/instance paths, overlaps,
retry and no-op behavior. All three groups reject the retained preceding library with
verified loader resolution. These are publication tests, not large-population performance
or measured memory-pressure qualification.

The five shared clients rebuild; all **21 focused checks** pass in 178.13 seconds.
Four affected production translation units compile independently with production flags;
controller, metadata, state, remaining-scene, database-material, sparse-entry and retained
probes pass. The initial implementation build's parameter-shadow warnings were corrected.
The initial sweep counted an externally equivalent empty-source overlap removal ambiguously;
its retry result now distinguishes preserved intent from committed removal. Both earlier
logs remain retained, separately from the final checks.

Fresh graphical regression rows pass: two eager LoD-off/on, ten delivery/adoption,
eight Generic cold/warm shaded/wire across both backends, and four primitive-editing rows
with 207 events each. All eight Generic camera contracts and 42 HUD pixel checks pass;
24 original images were inspected. Eager images retain four filled objects. Both partial-merge
denial images retain one filled leaf plus aggregate bounds and the red error HUD; adoption
denial retains four filled objects. Generic aircraft and editing geometry/manipulators
remain present. Warm Generic text overlaps the nose, and small editing panes retain
viewport text overlap/clipping, including the System GL ellipsoid status bar. Those S5
visual/layout issues remain open.
These rows qualify the changed boundary; large/cold/multiple-giant workloads, required native
hosts, full editing parity, latency/memory limits and known Lucy/planning/frame failures
retain their existing S4/S5 requirements.

Evidence is retained at
`.build/obol-qualification/20260909-retained-presentation-publication`, including before/after
files, diffs, negative and passing probes, scripts, reports, inspected originals, hashes and
runtime copies. The runtime manifest covers 49 artifacts, 84
resolved libraries and 132 retained physical files. `libBObol` is `a4860c61...`.
Obol `194d89fe...`, OSMesa `097aac27...`, the formal baseline and the user's proposal are
unchanged. Native dependency suites were not rerun for this BRL-CAD-only repair; their
preceding sealed results identify those unchanged binaries. Shared-stack sanitizers and
native-GPU qualification remain open. System GL here uses private-Xvfb CPU llvmpipe. All jobs are terminal and the final host-process
check is empty. Continue the noncompact GED region/metadata writers, enclosing operations
and remaining ownership inventory. **S1--S6 remain incomplete.**

## September 9: sparse entry publication and picking dispatch

`CXX-SPARSE-001` replaces partial mutation and full display rebuilds in seven entry-local
setters with one prepared publication: display state, transparency, selectability, region
ID, region metadata, full metadata and subtraction line style. It prepares changed records,
shader storage and membership capacity before committing numeric/semantic state, effective
styles, renderer records and compiler/journal effects. Shared intrinsic metadata and RGB
helpers also serve source-wide material refresh. Geometry, placements, retained overlays
and producer stamps retain their owners. Changed-field counts and notification behavior
are preserved; later reentrant edits survive outer callback failure.

The preceding library reproduces partial updates in five setters and resets transparency
and subtraction styling during an ordinary material edit. The final probe passes all seven.
Three production regression groups cover **28 scenarios and 382 C++ allocation positions:
322 preserved and 60 committed**. They also check exact/subtree/object and legacy matching,
instance suffixes and prefix boundaries, null/empty/missing queries, unchanged entries,
opacity clamping, invalid-color normalization, visibility journals, compiler retirement,
numeric revisions, effective summaries, immutable renderer geometry and actual picking.
All three groups reject the physical preceding libBObol with verified preload resolution.

The picking check exposed `CXX-PICK-001`: normal action traversal never called the retained
assembly's ray picker because its `SoNode` base supplies generic pick dispatch. BRL-CAD's
source wrapper explicitly calls the picker and had masked the missing registration.
`SoCADAssembly::initClass` now registers the existing ray dispatcher. The native regression
requires a translated wire hit, unpickable/hidden suppression and restoration; it rejects
the verified old Obol library. Derived BRL-CAD assembly traversal also passes. Numeric
picking algorithms and the source wrapper are unchanged.

Evidence is retained under `.build/obol-qualification/20260909-sparse-entry-publication`.
All **21 focused checks pass in 178.36 s**; all five application clients were rebuilt.
Four implementation translation units compile independently with their production flags,
including the scene controller's flags extracted from its actual unity compilation.
All five preceding public probes and the new sparse probe pass. Obol's full unit suite
passes **1,575 tests**, with one skip because profiling support is disabled; all **63
integration tests** pass. These native tests qualify the changed dependency rather than
relying on the earlier selected notification-test run.

All **24 graphical rows pass**: two eager, ten delivery, eight Generic and four editing rows
across OSMesa and System GL. Editing executes 207 events per row; eight Generic camera
contracts and **41 visible-HUD pixel checks** pass. Twenty-four original captures were
inspected. The partial-merge denial captures differ in prepared tier: OSMesa displays one
filled object and aggregate bounds; System GL retains leaf and aggregate boxes. Both meet
that row's retained-population, exact-bound, terminal-error and HUD criteria; this is not
a surface-quality pass. Text overlap/clipping and offscreen edit geometry remain S5 work.
Private-Xvfb System GL uses CPU llvmpipe, not native GPU qualification.

Initial test/tool failures are retained: the new fixture first assumed a pointer return
from an optional-valued renderer query and used exact scalar float comparisons. Its
unchanged positive picking control then exposed the real native dispatch defect; adding
a view-policy node alone did not repair it. The independent compilation script initially
requested a nonexistent standalone scene-controller object and now resolves its actual
unity object before compiling the source separately. No acceptance assertion was weakened.

Current libBObol is `d56e8010...`, Obol `194d89fe...` and OSMesa `097aac27...`.
Loader/hash evidence retains **48 runtime artifacts, 84 distinct resolved libraries and
131 physical runtime files**. Before/after files include both repositories; baseline
libraries, negative regressions, source diffs and the manifest retain the checkpoint.
OSMesa, the formal baseline and the user proposal are unchanged. Sanitizers remain off;
historical raw `/tmp` evidence is still unavailable.

**Next:** retained visibility/highlight/transparency rules, visibility frontiers,
selected-path replacement/deltas and their reapply functions, then GED's noncompact region
helpers and enclosing multi-source operations. Complete source/service, presentation,
resource and platform obligations remain open. **S1--S6 are incomplete.**

## September 9: database metadata and material refresh

`CXX-CONFIG-006` prepares explicit source metadata and owned records before notification.
Aggregate metadata applies to the source path, including slash/instance normalization;
distinct descendants and authored auxiliaries retain their own records. Material refresh
resolves each occurrence's full-path region, color and shader, prepares string storage,
then commits compact semantics, effective RGB, renderer records and compiler/batch effects.
Geometry, selection, visibility, line styling and presentation opacity retain their owners.
The existing compact stamp-only return/revision behavior is preserved.

Single-source, free bulk and scene callers share the same publication path. Bulk refresh
uses one material resolution sweep, pins the source set and advances the scene at its first
committed source. Exceptions retain completed prefixes; an enclosing sweep cannot overwrite
a later reentrant material revision. GED's duplicate post-publication shape writes are
removed. Across the four changed production files, the repair removes a net 15 lines.

Evidence is retained under `.build/obol-qualification/20260909-database-material-publication`.
**22 scenarios sweep 2,477 C++ allocation positions: 2,178 preserved and 299 committed.**
They cover explicit metadata, legacy/compact refresh, source/scene/free/bulk/batched callers,
quiet notification, aliases, retries and allocation-free no-ops. Throwing observers require
complete owned and compact records. Separate checks cover aggregate boundaries, borrowed
shaders, missing targets, stamp-only updates, committed prefixes and reentrant revisions.
A mixed-region production scene checks shaders and colors against librt's full-path oracle,
retained immutable geometry and actual compiled renderer records.

All **21 focused checks pass in 179.27 s**. All five application clients and affected
consumers were rebuilt; the three changed implementation translation units pass independent
compilation with production warning flags. All five retained public probes pass, including
both previously failing database metadata rows. All three new regression groups fail against
the preceding physical library; explicit preload and loader/hash verification identify it.

Initial failures remain retained separately: a parameter-shadow warning, a temporary shader
in the no-op fixture, and a copied string in the implementation's generic comparison. The
comparison now takes references. The region fixture initially retained its old RGB attribute
and now synchronizes combination attributes before writing. The renderer fixture initially
derived an ID from a source key and now uses the actual occurrence handle. These corrections
do not weaken the material, appearance or publication requirements.

All **24 graphical rows pass**: two eager, ten source-delivery, eight Generic and four editing
rows across OSMesa and System GL. Each editing row executes 207 events. Strict traces,
eight Generic camera contracts and **42 visible-HUD pixel checks** pass. Twenty-four original
captures were inspected. Text overlap/clipping and offscreen edit geometry remain S5 work;
private-Xvfb System GL uses CPU llvmpipe, not native GPU qualification.

Current libBObol is `d0b84e52...` and libged `e0d282d0...`; qged remains `282b9b6d...`,
Obol `adc2a691...` and OSMesa `097aac27...`. Loader/hash evidence covers **45 runtime
artifacts, 84 distinct resolved libraries and 128 retained runtime files**. Before/after files,
diffs and the manifest retain the checkpoint. Obol/OSMesa, the formal baseline and user
proposal are unchanged. The earlier 321 selected Obol notification tests were not rerun;
sanitizers remain off and historical raw `/tmp` evidence remains unavailable.

**Next:** sparse metadata/display setters and their GED callers, then the remaining
source/service and presentation inventory. This checkpoint does not qualify direct field/
engine callbacks, enclosing multi-step GED operations, complete resource bounds, large/cold
workloads or native hosts. **S1--S6 remain open.**

## September 9: scene state publication

`CXX-CONFIG-005` closes scene bounds, explicit realization, role-flags and view-policy
publication. Bounds/realization reuse their qualified inner commit with a scene hook.
Roles publish one scalar with compiler/batch retirement. View policy prepares its
fields and owner metadata with the existing configuration publication and retains the
threshold sensor. Unrelated geometry, compact appearance, authored auxiliary state,
nested ownership, source-production stamps and staged data survive. External failure
and exact-bound retention preserve view-only invalidation semantics. The old repeated
view-policy mutation and role/view-policy appearance rebuild paths are removed.
The [conformance record](libbobol_tla_conformance.md#cxx-config-005--remaining-scene-state-setters-lose-committed-revisions)
owns scope and the remaining inventory.

Evidence is retained under `.build/obol-qualification/20260909-scene-state-publication`.
**72 scenarios sweep 3,072 allocation positions:
2,892 preserved and 180 committed outcomes.**
They cover four setters, direct/scene/batched callers, current/internal-failed/external-
failed sources and ordinary/quiet notification. Immediate observers require complete
state and scene revisions. Throwing observers cannot skip later field notifications.
Separate compiler tests use eligible primary geometry; auxiliary sources exercise
retained editing state through their own rendering path. Reentrant commits, missing
targets, normalized no-ops, numeric tolerance, borrowed bounds/diagnostics and subsequent
threshold invalidation have explicit checks.

All **21 focused checks pass in 172.52 s**. All five application clients and affected
focused consumers were rebuilt; the two changed library translation units pass
independent syntax checks with production warning flags. The controller, metadata,
display/placement and four remaining-scene probes all pass. All three new regression
groups fail against the preceding physical library, with its loader resolution retained.
An initial negative attempt resolved the current library because the retained basename
did not match the versioned SONAME; that passing run is not negative evidence. Accepted
negative runs explicitly preload the retained library. Initial fixture builds accessed
a private batch revision, used prohibited float comparisons and omitted an aggregate
initializer. The first runtime fixture incorrectly required compiled rendering for a
source with auxiliary children. These fixture failures and corrections remain retained;
compiler and auxiliary behavior now have separate valid tests. The next metadata probe's
initial missing include/build flags are also retained separately from its reproduction.

All **24 graphical rows pass**: two eager LoD off/on/erase, ten healthy/injected source-
delivery failure, eight Generic cold/warm shaded/wire and four single/quad editing rows
across OSMesa and System GL. Each editing row executes 207 events. Strict traces,
eight Generic camera contracts and **41 Generic visible-HUD pixel checks** pass,
together with eager/delivery HUD checks. Twenty-four original captures were inspected.
Annotation overlap, offscreen edit geometry and HUD clipping remain S5 work. Private-Xvfb
System GL is CPU llvmpipe; native GPU and broader capability qualification remain open.

Current qged is `282b9b6d...`, libBObol `9e4d2be9...`,
libged `a5ab140c...`, libqtcad `45465d09...`,
Obol `adc2a691...` and OSMesa `097aac27...`.
Loader resolution and full hashes cover **45 runtime artifacts, 84 distinct resolved
libraries and 128 retained runtime files**. Before/after files and diffs preserve
the checkpoint. Obol/OSMesa, the formal baseline catalog and user proposal are unchanged.
The earlier 321 selected Obol notification tests were not rerun. Sanitizers remain off;
historical raw `/tmp` evidence remains unavailable.

**Next: `CXX-CONFIG-006`, database metadata and material refresh.** The retained probe
exits 1 with both rows incomplete after observer exceptions. These failures are outside
the passing suite. Close those setters and their enclosing callers, then continue the
source/service, sparse-journal and presentation inventory. S1--S6 remain open.

## September 9: display publication

`CXX-CONFIG-004` now closes full display-state and display-patch publication with
scene effects, retaining the preceding placement repair. Both setters use one
prepared field mapping. Source/owned properties, compact membership and appearance,
source invalidation and scene revisions commit before observers. Revision sensors
remain attached. Unrelated selection, opacity, names and CSG line patterns survive;
intrinsic database colors and sparse visibility/highlight/transparency precedence
remain intact. Registry selection stays occurrence-owned, and retired overviews
cannot be resurrected. Compact preparation allocates membership/mask buffers without
cloning geometry or the occurrence registry. The
[conformance record](libbobol_tla_conformance.md#cxx-config-004--display-and-placement-writers-can-interrupt-their-effects)
owns detailed scope and the next inventoried boundary.

Evidence is retained under `.build/obol-qualification/20260909-display-publication`.
**72 scenarios sweep 6,948 allocation positions:
6,768 preserved and 180 committed outcomes.**
They cover direct/scene/batched calls, ordinary occurrences/registries, realized/failed
sources, ordinary/quiet notification, throwing display/revision observers, sparse
overlays and source/input invalidation. Journal allocation failure preserves the
complete mutation and requires authoritative rescan. Source no-ops allocate nothing;
independent subsequent source/input changes verify both sensors. Production database
combinations exercise CSG line patterns, legacy/compact material precedence, renderer
color round trips, unrelated styling and geometry identity. Reentrant patches,
ignored/borrowed arguments, overlay removal and retired overviews have explicit checks.
All three new regression groups fail against the preceding physical library.

All **21 focused checks pass in 163.95 s**. The five application clients and
affected focused-test consumers were rebuilt. All three changed library translation
units pass independent syntax checks with production warning flags. The initial
test build accessed private realization-stamp fields; the fixture now captures
public revisions. A subsequent assertion treated the streamed whole-target overview
as a leaf and placed a visibility rule on the whole root. The fixture now targets
explicit leaf paths while retaining the overview. Both failed runs and their test
inputs are retained separately; these are fixture corrections, not library repairs.
Controller, metadata and original state probes all pass. Obol/OSMesa are unchanged;
the earlier 321 selected Obol notification tests were not rerun. Sanitizers remain off.

All **24 graphical rows pass**: two eager LoD off/on/erase, ten healthy/injected
source-delivery-failure, eight Generic cold/warm shaded/wire and four single/quad
editing rows across OSMesa and System GL. Each editing row executes 207 events.
Strict traces, eight Generic camera contracts and **43 Generic visible-HUD pixel
checks** pass, together with eager/delivery HUD checks. Twenty-four original captures
were inspected. Annotation overlap, offscreen edit geometry and HUD clipping remain
S5 work. Private-Xvfb System GL is CPU llvmpipe; native GPU and broader capability
qualification remain open.

Current qged is `282b9b6d...`, libBObol `c72952c7...`,
libged `a5ab140c...`, libqtcad `45465d09...`,
Obol `adc2a691...` and OSMesa `097aac27...`.
Full hashes and loader resolution cover **44 runtime artifacts, 84 distinct resolved
libraries and 127 retained runtime files**. Before/after files and diffs preserve
the checkpoint. The formal baseline catalog and user proposal remain unchanged;
historical raw `/tmp` evidence remains unavailable.

**Next: `CXX-CONFIG-005`, remaining scene state publication.** The retained
`remaining_scene_probe.cpp` exits 1: bounds and explicit realization commit their
inner source state but lose the scene revision after observer exceptions. Role-flags
and view-policy changes also leave incomplete effects. These failures are outside
the passing suite. Continue those effects and the remaining database metadata/material
refresh, source/service, sparse-journal and presentation inventory. S1--S6 remain open.

## September 8: placement publication

The placement part of `CXX-CONFIG-004` now commits source and owned metadata,
transform nodes and child paths, compact matrices, compiler/batch/LoD invalidation
and scene revisions before observers. Ordinary nested sources keep their own state;
auxiliaries recursively inherit metadata without another transform. Shared owned
shapes publish once. Geometry, appearance, source bounds, realization stamps and
staged producer results remain retained. Graph-only normalization repairs missing,
reordered and duplicate transforms. Borrowed arguments, numeric no-ops and reentrant
callbacks retain their semantics. The
[conformance record](libbobol_tla_conformance.md#cxx-config-004--display-and-placement-writers-can-interrupt-their-effects)
owns scope and remaining display work.

Evidence is retained under `.build/obol-qualification/20260908-placement-publication`.
**39 scenarios sweep 8,852 allocation positions:
8,795 preserved and 57 committed outcomes.**
They include direct/scene/batched calls, ordinary/quiet notification, realized/failed
sources, insertion/movement/removal and path-audited graph repair. Independent checks
exercise scaled/translated production compact renderer matrices, unchanged local
geometry and auxiliary traversal bounds. The three new regression groups fail
against the preceding physical library. The retained state probe now reports
`frame_advanced=1 state_complete=1` for placement despite its throwing observer;
its two display rows remain broken, so the whole probe still exits 1.

All **21 focused checks pass in 153.61 s**. The five application clients and
affected focused-test consumers were rebuilt. Both changed library translation
units pass independent syntax checks with production warning flags. The initial
auxiliary bounding-box test used a default-constructed, uninitialized vector as
zero. Its failed output and initial test fragment are retained; the corrected test
uses explicit coordinates. This is a fixture correction, not another library fix.
Controller and metadata probes remain passing. Obol/OSMesa are unchanged; the
earlier 321 selected Obol notification tests were not rerun. Sanitizers remain off.

All **24 graphical rows pass**: two eager LoD off/on/erase, ten healthy/injected
source-delivery-failure, eight Generic cold/warm shaded/wire and four single/quad
editing rows across OSMesa and System GL. Each editing row executes 207 events.
Strict traces, eight Generic camera contracts and **41 Generic visible-HUD pixel
checks** pass, together with eager/delivery HUD checks. Twenty-four original
captures were inspected. Existing annotation overlap, offscreen edit geometry and
HUD clipping remain S5 work. Private-Xvfb System GL is CPU llvmpipe; native GPU and
the broader capability envelope remain unqualified.

Current qged is `282b9b6d...`, libBObol `823f63ef...`,
libged `a5ab140c...`, libqtcad `45465d09...`,
Obol `adc2a691...` and OSMesa `097aac27...`.
Full hashes and loader resolution cover **43 runtime artifacts, 84 distinct resolved
libraries and 126 retained runtime files**. Before/after files and diffs preserve
this checkpoint. The formal baseline catalog and user proposal remain unchanged;
historical raw `/tmp` evidence remains unavailable.

**Next: finish `CXX-CONFIG-004` display-state/display-patch publication and scene
effects**, then continue the remaining revision/material, source/service,
sparse-journal and presentation inventory. S1--S6 remain incomplete.

## September 8: source metadata publication

`CXX-CONFIG-003` now closes display-name, hierarchy and material-policy publication
and scene revision effects. Source and changed owned metadata publish before observers,
with compact appearance/parent records and compiler/batch invalidation. Preparation
failure preserves preceding state; notification failure retains a complete commit.
Authored non-database names/materials, nested ownership, per-occurrence overrides,
geometry, realization stamps and staged producer data remain retained. Compact
occurrences keep their intrinsic full-path database color independently of inherited
appearance. Material-only changes preserve line patterns, width and opacity; hierarchy
preserves leaf occurrence indices/operations and retires obsolete parent cache evidence.
The [conformance record](libbobol_tla_conformance.md#cxx-config-003--metadata-setters-can-interrupt-propagation-and-scene-revisions)
owns detailed semantics and remaining boundaries.

Evidence is retained under `.build/obol-qualification/20260908-source-metadata`.
**36 scenarios sweep 2,940 allocation positions: 2,886 preserved and 54 committed
outcomes.** They cover direct/scene/batched setters, realized/failed sources, enabled/
quiet notifications, complete observer state, aliases, normalized no-ops, reentrant
callbacks and subsequent input changes. A production database combination with
subtraction checks legacy/compact color round trips, unchanged geometry and styling,
renderer records, hierarchy cache retirement and source-name/path fallbacks.
`metadata-counts.json` owns individual counts. All three new regression groups fail
against the preceding physical library. The corrected metadata probe exits 0, with
complete metadata and advanced scene frames for all three rows; the preceding library
fails it. The original name probe incorrectly expected an authored diagnostic prototype
to be renamed. Its source/output remain retained alongside the corrected primary-
geometry fixture, so the correction is not counted as a library repair.

All **21 focused checks pass in 122.86 s**. The five application clients and affected
focused-test consumers were rebuilt. Both touched library translation units pass
independent syntax checks with production warning flags. The initial test build used
an unavailable public revision member and a forbidden floating-point equality; the
first sweep also incorrectly required allocation-free scene lookup. These fixture
failures were corrected. The initial focused run then caught a real missing compiler
invalidation for hierarchy, which is repaired and covered by the renderer regression.
All failed runs remain separate from final qualification. Obol/OSMesa are unchanged;
the preceding 321 selected Obol notification tests were not rerun. This Debug BRL-CAD/
Release Obol stack still has sanitizers off.

All **24 graphical rows pass**: two eager LoD off/on/erase, ten healthy/injected
source-delivery-failure, eight Generic cold/warm shaded/wire and four single/quad
editing rows across OSMesa and System GL. Each editing row executes 207 events.
Strict traces, eight Generic camera contracts and **41 Generic visible-HUD pixel
checks** pass, with eager/delivery HUD checks. Twenty-four original captures were
inspected. Existing annotation overlap, offscreen edit geometry and HUD clipping
remain S5 acceptance work. Private-Xvfb System GL is CPU llvmpipe; native GPU and
broader capability qualification remain open.

Current qged is `282b9b6d...`, libBObol `ec7e4621...`, libged `a5ab140c...`,
libqtcad `45465d09...`, Obol `adc2a691...` and OSMesa `097aac27...`.
Full hashes and loader resolution cover **43 runtime artifacts, 84 distinct resolved
libraries and 126 retained runtime files**. Before/after files and diffs preserve the
checkpoint. The formal baseline catalog and user proposal remain unchanged;
historical raw `/tmp` evidence remains unavailable.

**Next: `CXX-CONFIG-004`, display/placement publication and scene effects.** The
retained state probe exits 1: full display updates, display patches and placement
leave incomplete state and skip their scene frame revision after observer exceptions.
These failures remain outside the passing suite. Continue those writers and the
remaining revision/material, source/service, sparse-journal and presentation inventory.
S1--S6 remain incomplete.

## September 8: scene configuration publication

`CXX-CONFIG-002` now closes the scene draw/representation/rename boundary. Source
configuration and scene revisions publish before observers. Rename also prepares
only affected indexes, conflict child/path removal and repository retirement;
allocation rejection preserves the preceding whole state, while notification
failure retains a complete commit. Removed nested owners lose their indexes and
cache leases; shared descendants retain their surviving parent and one path
membership. Compatible external draw geometry stays realized across the controller's
paired representation update. Explicit revision-only rename is honored, and graph
conflicts which cannot be replaced coherently are rejected before mutation. The
[conformance record](libbobol_tla_conformance.md#cxx-config-002--scene-configuration-effects-can-be-interrupted)
owns detailed semantics and remaining boundaries.

Evidence is retained under `.build/obol-qualification/20260908-scene-configuration`.
The original controller probe now exits 0: draw, representation and rename all report
advanced frames, retained source indexes and complete draw channels. Both new
regression groups fail against the preceding physical libBObol. **36 scenarios
sweep 4,868 allocation positions: 4,764 preserved and 104 committed outcomes.**
They cover ordinary/batched revisions, normal/quiet notifications, shared-path
rename, nested conflict retirement, cached geometry references, path truncation,
source stamps, lookup/order/routing consistency, compatible external mesh/wire
geometry, throwing observers and retry/no-op behavior. Edge checks cover aliases,
normalized no-ops, rejected ancestor/shared conflicts and reentrant changes.
`scene-configuration-counts.json` owns individual counts.

All **21 focused checks pass in 120.83 s**. The five application clients and affected
focused-test consumers were rebuilt. The three touched library translation units
also pass independent syntax checks with the production warning flags. Initial
include-path, shadowed-name and standalone include failures were corrected, as was
an initially unreleased raw-pointer cache reference caught by the lease regression.
Their logs, including a repeat against the preceding executable after a failed
build, remain separate from final qualification. Obol and OSMesa are unchanged;
the preceding checkpoint's 321 selected Obol notification tests were not rerun.
This Debug BRL-CAD/Release Obol stack still has sanitizers off.

All **24 graphical rows pass**: two eager LoD off/on/erase, ten healthy/injected
source-delivery-failure, eight Generic cold/warm shaded/wire and four single/quad
editing rows across OSMesa and System GL. Each editing row executes 207 events.
Strict traces, eight Generic camera contracts and **39 Generic visible-HUD pixel
checks** pass, together with eager/delivery HUD checks. Twenty-four original
captures were inspected. Existing annotation overlap, offscreen edit geometry and
HUD clipping remain S5 acceptance work. Private-Xvfb System GL is CPU llvmpipe;
native GPU and broader capability qualification remain open.

Current qged is `3be2bee2...`, libBObol `90de6fb3...`, libged `0fb3ddd9...`,
libqtcad `45465d09...`, Obol `adc2a691...` and OSMesa `097aac27...`.
Full hashes and loader resolution cover **42 runtime artifacts, 84 distinct
resolved libraries and 125 retained runtime files**. Before/after files and diffs
preserve the checkpoint. The formal baseline catalog and user proposal remain
unchanged; historical raw `/tmp` evidence remains unavailable.

**Next: `CXX-CONFIG-003`, metadata setters and scene effects.** The retained probe
exits 1: display-name and hierarchy observers leave incomplete metadata; display-name,
hierarchy and material-policy wrappers all skip their scene frame revision. The
material row proves the committed policy field, not wider propagation. These
failures are outside the passing suite. Continue prepared metadata publication and
its scene effects, then the existing display/revision, source/service, sparse-journal
and presentation inventory. S1--S6 remain incomplete.

## September 8: source configuration setters

`CXX-CONFIG-001` now covers batched configuration and all three separate source
setters through one prepared publication boundary. Each setter retains its own
defaulting, external-realization and revision semantics. Source and owned metadata,
producer-data retirement and compact draw channels commit before observers while
all field sensors remain attached. Representation refresh preserves descendant
paths, and retargeting maps their suffixes without matching unrelated prefixes.
Materials retain multi-value properties; shared children publish once and nested
sources keep their ownership. The
[conformance record](libbobol_tla_conformance.md#cxx-config-001--interrupted-configuration-leaves-mixed-fields-and-lost-sensors)
owns the detailed scope.

Evidence is retained under `.build/obol-qualification/20260908-configuration-setters`.
The configuration, setter and path probes now exit 0. Each of the three setter
rows reports `threw=1 invalidated_at_exception=1 later_sensor_recovered=1`.
The preceding physical libBObol fails the three new regression groups. Across
**28 new scenarios, 3,512 allocation positions produce 3,438 preserved and 74
committed outcomes**. These cover draw, representation, path/instance/revision
changes, failed/realized state, quiet/normal notification and compatible external
wire/mesh geometry. Further checks cover throwing observers, subsequent unchanged
and changed field touches, aliased inputs, path boundaries, null keys and invalid
input. `setter-counts.json` retains individual sweep counts. An initial missing
constructor initializer and a path-normalization regression were fixed; their
failed build/test logs remain separate from the passing final results.

All **21 focused checks pass in 105.79 s**. The five application clients and
affected test consumers were rebuilt. Obol and OSMesa are unchanged from the
preceding checkpoint; its 321 selected Obol notification tests were not rerun.
The Debug BRL-CAD/Release Obol stack still has sanitizers off.

All **24 graphical rows pass**: two eager LoD off/on/erase, ten healthy/injected
delivery-failure, eight Generic cold/warm shaded/wire and four single/quad editing
rows across OSMesa and System GL. Each editing row executes 207 events. Strict
traces, eight Generic camera contracts and **44 Generic visible-HUD pixel checks**
pass, along with eager/delivery HUD checks. Twenty-four original captures were
inspected. Existing annotation overlap, offscreen edit geometry and HUD clipping
remain S5 acceptance work. Private-Xvfb System GL is CPU llvmpipe; native GPU and
broader capability qualification remain open.

Current qged is `3be2bee2...`, libBObol `f045e1ce...`, libged `e6117229...`,
libqtcad `45465d09...`, Obol `adc2a691...` and OSMesa `097aac27...`.
Full hashes and loader resolution cover **44 runtime artifacts, 84 distinct
resolved libraries and 127 retained runtime files**. Before/after sources and
diffs preserve the checkpoint. The formal baseline catalog and user proposal are
unchanged; historical raw `/tmp` evidence remains unavailable.

**Next: `CXX-CONFIG-002`, enclosing scene configuration.** The retained controller
probe exits 1: observer exceptions skip frame revisions for draw/representation,
interrupt the paired draw-channel update, and leave a renamed source absent from
the index. These failures are outside the passing suite. Source/index/conflict and
revision effects need one prepared boundary. Remaining display/revision writers,
source/service, sparse journals, presentation and all S1--S6 gates remain open.

## September 8: batched source configuration

The batched portion of `CXX-CONFIG-001` now prepares source configuration and
owned legacy metadata before publication. Existing field sensors stay attached;
Obol's field notification can omit the sensor whose effects have already committed
without suppressing a reentrant change. Database binding, exact-bound/source-data
retirement and compact draw channels commit before observers. Shared children
publish once and nested sources retain ownership; geometry, occurrence identities
and sparse overrides remain retained. Unchanged configuration allocates and notifies
nothing. The [conformance record](libbobol_tla_conformance.md#cxx-config-001--interrupted-configuration-leaves-mixed-fields-and-lost-sensors)
owns the detailed boundary and remaining setter obligations.

Evidence is under `.build/obol-qualification/20260908-source-configuration`.
The original probe now exits 0 with `threw=1 partial_configuration=0
later_input_sensor_recovered=1`. The preceding physical libBObol fails the new
regression. Eight configuration scenarios cover identity/draw changes, realized/
failed initial states and enabled/disabled notifications. Their **778 allocation
positions produce 740 preserved and 38 committed outcomes**. Additional tests cover
throwing observers, all twelve retained sensors, reentrant writes, database A/B/A,
aliased strings, default identities, auxiliary ownership and no-op calls.
`configuration-counts.json` owns exact per-case counts. The first draw-only fixture
also changed its empty representation identity; its failure is retained separately
from the corrected, normalized fixture. An initial launcher addressed the test
executable under the wrong build directory and launched no test; that error is
likewise separate from test results.

All **21 focused checks pass in 99.90 s**. The five application clients and affected
test consumers were rebuilt. Obol's existing field/sensor/notification selection
also passes **321 tests**, with zero disabled tests, failures or errors. Obol was
rebuilt and installed into the explicit workspace prefix; the matching OSMesa
remained unchanged. This Debug BRL-CAD/Release Obol stack has sanitizers off.

All **24 graphical rows pass**: two eager LoD off/on/erase, ten healthy/injected
source-delivery-failure, eight Generic cold/warm shaded/wire and four single/quad
editing rows across OSMesa and System GL. Each editing row executes 207 events.
Strict traces, eight Generic camera contracts and **38 Generic visible-HUD pixel
checks** pass, along with eager/delivery HUD checks. Twenty-four original captures
were inspected. Existing small-pane annotation overlap, offscreen edit geometry and
HUD clipping remain S5 acceptance work. Private-Xvfb System GL is CPU llvmpipe;
these results do not qualify a native GPU or the broader capability envelope.

Current qged is `3be2bee2...`, libBObol `d64fea0c...`, libged `e6117229...`,
libqtcad `45465d09...`, Obol `adc2a691...` and OSMesa `097aac27...`.
Full hashes and loader resolution cover 43 runtime artifacts, 84 distinct resolved
libraries and 126 retained runtime files. Before/after sources preserve this
checkpoint. The formal baseline catalog and user proposal are
unchanged; historical raw `/tmp` evidence remains unavailable.

**`CXX-CONFIG-001` remains open.** The separate draw, representation and retarget
setter probe exits 1; each row reports `threw=1 invalidated_at_exception=0
later_sensor_recovered=0`. Those setters still detach sensors and publish fields
before their remaining work. These failing rows are outside the passing suite.
Extend prepared publication to them, then continue the existing writer inventory.
All S1--S6 gates remain incomplete; retained planning/frame stalls, Lucy quality,
large/cold/shared resource bounds, bounded BREP, full editing/native coverage,
shared-stack sanitizers and coherent final dependency qualification remain required.

## September 8: source invalidation

`CXX-SOURCE-025` publishes the source's stale reason, realization record, owned
legacy records and applicable exact-bound/source-data retirement before observers.
Rejected preparation preserves the preceding state and staged data; throwing
observers retain the complete invalidation. Compact geometry and display records
remain untouched. Configuration callers now explicitly own their metadata/display
propagation. The [conformance record](libbobol_tla_conformance.md#cxx-source-025--invalidation-publishes-before-resource-and-owner-retirement)
owns the contract and remaining writer limits.

Evidence is under `.build/obol-qualification/20260908-source-invalidation`.
The original probe previously observed partial records with and without an
exception. It now exits 0 and reports `partial_observation=0`, `stale=1`,
`reason=1`, `status=0`, `exact=0`, `owner_stale=1` in both cases. The preceding
physical library fails the new permanent regression. Its **42 scenarios** cover
seven reason combinations, current/internal-failed/external-failed states and
enabled/disabled notifications. The final suite exercises **398 allocation
positions: 336 preserved and 62 committed outcomes**, plus throwing/removing
observers and zero-reason calls. `invalidation-counts.json` owns exact counts;
this is logical publication/lifetime evidence, not large-population memory or
cancellation qualification.

All **21 focused checks pass in 96.98 s**. All five application clients and
affected test consumers were rebuilt. An initial full run exposed draw-mode
propagation which had depended on `markStale`; the representation and retarget
callers now retain that responsibility explicitly. A subsequent wrapper exited
143 before a complete CTest verdict; no surviving test process remained, and its
log is retained separately from the final passing rerun. Test-fixture initialization
and build corrections likewise retain their earlier failures.

All **24 graphical rows pass** across OSMesa and System GL: two eager LoD
off/on/erase, ten healthy/injected-delivery-failure, eight Generic cold/warm
shaded/wire, and four single/quad primitive-edit rows. Strict traces, eight Generic
camera contracts and **42 Generic visible-HUD pixel checks** pass; eager/delivery
HUD checks also pass. Each editing row executes 207 events. Twenty-four original
representative captures were inspected. Annotation overlap, offscreen edit geometry
and HUD clipping remain S5 work. Private-Xvfb System GL is CPU llvmpipe; this Debug
stack has sanitizers off.

Current qged is `9d25f692...`, libBObol `0d7c40f3...` and libged `33b2ab28...`.
Obol `54ee73b2...` and OSMesa `097aac27...` are unchanged. Full hashes and loader
resolution are retained with 42 runtime artifacts, 84 distinct resolved libraries
and 125 retained runtime files. Unchanged retained content is shared by hard link;
live build outputs are never linked into evidence. Before/after files, diffs,
scripts, reports, images and final process/manifest audits preserve this checkpoint.

**Open `CXX-CONFIG-001`:** the configuration probe still exits 1 with `threw=1
partial_configuration=1 later_input_sensor_recovered=0`. A new instance key and
old path remain together after an observer exception, and subsequent realization
does not recover the detached input-revision sensor. This failing probe is outside
the passing suite. Prepare the enclosing configuration/owned-metadata boundary
while preserving sensor ownership, then resume the remaining source/service,
sparse-journal and presentation audit. All S1--S6 gates remain incomplete; retained
planning/frame stalls, large/cold/shared resource limits, BREP production, Lucy
quality, full editing/native coverage, shared-stack sanitizers, evidence recovery
and coherent final dependency staging remain mandatory acceptance work.

## September 8: realization render requests

`CXX-REALIZE-002` extends the existing source/action/scene commit to the view's
standing render request. It prepares reasons and policy classification before
source writes, merges through the existing request owner, and notifies the host
even when source observers throw. Completed prefixes retain their work; a request
consumed by the host is not recreated at return. Ordinary render requests share
the same prepared boundary. Diagnostic trace cleanup now records missing evidence
without throwing during operation exit. The
[conformance record](libbobol_tla_conformance.md#cxx-realize-002--observer-failure-bypasses-the-view-render-request)
owns the contract and remaining writer limits.

Evidence is under
`.build/obol-qualification/20260908-realization-render-request`. The original
scene and view probes both exit 0; the view reports `threw=1 source_committed=1
frame_before=1 frame_after=2 render_requested=1`. With the preceding physical
library preloaded, the new view and ordinary request tests both exit 1, and the
trace-unwind test aborts on an allocating scope destructor. Those negative runs
retain their library identity and exact exits; they are expected failure evidence.

All **21 focused checks pass in 96.95 s**. The new permanent tests cover 32 view
scenarios across prototype, cached wire/mesh and failed realization, tracing
off/on, existing/empty requests and explicit LoD off/auto. Their **8,542 allocation
positions** produce 6,488 preserved and 2,054 committed outcomes. Four ordinary
request sweeps cover **40 positions**, and two three-source prefix sweeps cover
**538 positions** with every complete prefix represented. These **9,120 new
positions**, first-exception checks, host consumption and targeted trace-unwind
checks supplement the preceding source/action/scene coverage. `request-counts.json`
owns per-case counts from the final suite. Incomplete traces under deliberately
denied allocation are explicit diagnostic gaps, not valid formal witnesses.

An initial test incorrectly assumed the default internal view attachment had LoD
disabled. Its failure and diagnostic runs are retained separately. The corrected
fixture sets the policy explicitly and tests both settings; the final production
regression suite passes. All five application clients and affected realization
action consumers were rebuilt. This Debug stack has sanitizers off.

All **24 graphical rows pass**: two eager LoD off/on/erase rows, ten healthy and
injected-delivery-failure rows, eight Generic cold/warm shaded/wire rows, and four
single/quad primitive-edit rows across OSMesa and System GL. All strict traces
pass, as do eight Generic camera contracts and **43 Generic visible-HUD pixel
checks**; eager and delivery HUD checks also pass. Each editing row executes
207 events with command/readback and feature assertions. Twenty-four original
representative images were inspected. Existing viewport annotation overlap,
offscreen edit geometry and clipped HUD text remain S5 visual/layout acceptance
work. Private-Xvfb System GL is CPU llvmpipe, not native GPU qualification.

Current qged is `dcddc951...`, libBObol `b0293b4e...`, libged `3f01b818...` and
libqtcad `025fbc61...`. Obol `54ee73b2...` and OSMesa `097aac27...` are unchanged.
`current-hashes.json` retains full hashes; `runtime.json` records 42 artifacts,
84 distinct resolved libraries and 125 physical runtime copies. Before/after
files, diffs, scripts, reports, images and process/manifest audits preserve the
checkpoint.

All S1--S6 gates remain incomplete. Resume configuration/invalidation, then the
remaining source/service/sparse and presentation boundaries. This repair does
not qualify the retained planning/frame stalls, Lucy quality, large/cold/shared
resource bounds, bounded BREP production, full editing/platform coverage,
shared-stack sanitizers or coherent dependency staging. Historical evidence
recovery and the full formal rerun remain required before release acceptance.

## September 8: realization effects

`CXX-REALIZE-001` extends the existing source commit to action progress and scene
frame publication. Failed-source diagnostics are prepared before publication;
committed counts, diagnostics and frame revisions precede source observers.
Each changed source has its own frame revision; explicit mutation batches still
coalesce. Scoped scene progress survives exceptions without allocating during
unwind. The [conformance record](libbobol_tla_conformance.md#cxx-realize-001--observer-failure-bypasses-scene-frame-publication)
owns the boundary and remaining writer limits.

Evidence is under `.build/obol-qualification/20260908-realization-effects`.
The preceding prototype-qualified library fails the new observer/progress test.
The original probe now prints `threw=1 source_committed=1 frame_before=1
frame_after=2` and exits 0. Twelve action/scene scenarios cover prototype, compact
wire/mesh, evaluated wire/points and failed database realization. Full allocation
sweeps cover all except evaluated points, which retains first/last denial and
success samples because its ray-sampling provider is expensive. The two three-source
prefix sweeps include rejection before any commit and interruption after each
completed prefix. Additional checks cover throwing observers, accumulated failure
diagnostics, custom traversal failure, retry, disabled notifications and explicit
mutation batches. The final suite exercises **2,678 single-source allocation positions** (2,142
old and 536 committed outcomes), plus **519 three-source prefix positions**.
Exact per-case counts are retained in `realization-counts.json`.

All **21 focused checks pass in 78.61 s**, including OSMesa prototype/HUD rendering
and the rebuilt NanoRT action client. All five application clients and affected
action consumers were rebuilt. Initial build warnings and successful follow-up
builds are retained separately. The `bobol_sanitizer` test label does not imply
this Debug build enabled sanitizers.

All **ten graphical rows pass**: two eager LoD off/on/erase rows and eight Generic
cold/warm, shaded/wire rows across OSMesa and System GL. Their strict traces,
eight camera contracts and **42 visible-HUD pixel checks** pass. Twelve original
representative captures were inspected. Existing viewport annotation overlap and
clipping remain S5 debt. Private-Xvfb System GL is CPU llvmpipe, not native GPU
qualification. The previous delivery-fault and primitive-edit rows were not
rerun here; they retain their dated binary provenance.

Current libBObol is `776794bc...` and libged `700ef390...`; qged remains
`be39cdbf...`. Obol `54ee73b2...` and OSMesa `097aac27...` are unchanged.
`current-hashes.json` retains full hashes. `runtime.json` records 42 artifacts,
84 distinct resolved libraries and 125 physical copies, including both probes
and the rebuilt NanoRT consumer. The checkpoint retains before/after files,
diffs, commands, reports, images and final process/manifest audits.

**`CXX-REALIZE-002` found here, subsequently
[repaired](#september-8-realization-render-requests):** the separate
`view_realization_probe` exits 1 with
`frame_before=1 frame_after=2 render_requested=0` after an observer exception.
The scene commit survives, but `BObolViewController::realizePending` bypasses its
render request. The probe remains outside the passing suite and is not a pass.

All S1--S6 gates remain incomplete. This checkpoint does not qualify the remaining
view/host successor, configuration/invalidation, source/service/sparse writers,
Lucy quality, large/cold/shared-cache resource bounds, full editing/platform
coverage, shared-stack sanitizers or coherent dependency staging. Historical
`/tmp` evidence recovery and the full formal rerun remain required.

## September 8: prototype publication

`CXX-SOURCE-024` replaces live prototype edits with private construction followed
by the existing prepared source publication. The diagnostic square retains its
revision-dependent geometry and intent. Its source-local bounds now describe
that geometry, and auxiliary lines are preserved even when they precede the
primary drawing. Source/owner records, paths and obsolete compact/compiled
retirement publish coherently. The [conformance record](libbobol_tla_conformance.md#cxx-source-024--prototype-realization-mutates-a-live-drawing)
owns the exact source-level boundary.

Evidence is under `.build/obol-qualification/20260908-prototype-publication`.
The new regression rejects the preceding policy-qualified library: 308 immediate
callbacks include partial scene state. With the repair, all eight direct/action
cases pass; the standalone sweep exercises **1,528 allocation positions, 1,305
preserved and 223 committed outcomes**. The default suite repeat has its own
allocation counts in `prototype-counts.json`. Cases cover compact retirement,
auxiliary and nested ownership, placement, disabled notifications, reference/path
recovery, retry and throwing observers. Existing prototype tests also exercise
stale-field refresh, picking, snapping, measurement and export.

All **20 focused checks pass in 58.65 s**, including the actual prototype and
color/HUD rendering checks through OSMesa. All five clients and the affected
executables were rebuilt. The evidence root retains the successful build and
preceding attempts; interrupted or failed attempts are not qualification.

All **24 graphical rows pass** their applicable functional, strict-trace and HUD
checks: two eager policy rows, ten deferred-delivery rows, eight Generic
shaded/wire cold/warm rows and four single/quad primitive-edit rows. The Generic
rows also pass eight camera contracts and **42 visible-HUD pixel checks**; the
editing rows each execute 207 events. Twenty-four original representative
captures were inspected. The existing viewport text overlap/clipping remains
S5 debt. The number of visible-HUD checks varies with the captured state; this
run's per-row counts and the unchanged predicate are retained.

**Controller defect found at this checkpoint (`CXX-REALIZE-001`, subsequently
repaired [above](#september-8-realization-effects)):** the retained
`scene_revision_probe.cpp` and executable call `BObolSceneController::realizePending`
with an immediate source observer which throws. The source commits revision 2,
but the scene frame revision remains 1. The probe exits 1 and is not part of the
passing default suite. This demonstrates missing controller bookkeeping after
publication, independently of the repaired source transaction. Repair the common
realization effect boundary, including observer ordering, cached success/failure
and partial traversal, before claiming controller-level exception recovery.

Current qged is `be39cdbf...`, libBObol `1f0dfce1...`, libged `3c6da76a...`,
libbv `7120bb05...` and libqtcad `95758c5e...`; Obol `54ee73b2...` and OSMesa
`097aac27...` are unchanged. `current-hashes.json` retains full hashes.
`runtime.json` records 40 artifacts, 84 distinct resolved libraries and 123
physical copies, including the prototype renderer test and failing controller
probe. Before/after files, diffs, negative evidence, commands, reports and images
are retained with a final process audit and manifest.

All S1--S6 gates remain incomplete. The new controller defect, configuration and
invalidation, remaining source/service/sparse writers, Lucy quality, large/cold
and multiple-giant resource behavior, complete client workflows and native hosts
remain open. System GL in the private-Xvfb rows is llvmpipe. Formal relations and
strict graphical predicates are unchanged, and the user's proposal is preserved.

## September 8: live policy publication

`CXX-POLICY-001` repairs the publication order exposed by the preceding graphical
failures. GED previously changed the attachment before the controller retired
work and advanced its policy revision. The controller now owns those effects as
one external-input operation. Equivalent sanitized policies preserve current
work, and disabling LoD preserves the retained drawing. The shared comparison
helper removes duplicate policy tests; the adapter's clear/revision sequence is
removed. The [conformance record](libbobol_tla_conformance.md#cxx-policy-001--live-policy-published-before-its-control-epoch)
owns the exact boundary.

Evidence is retained under
`.build/obol-qualification/20260908-policy-publication`:

- All **19 focused checks pass in 53.09 s**, including BV, the controller policy
  regression, source-publication allocation sweeps and GED/Qt integration.
- Both eager-draw rows pass the unchanged strict trace, command, bounds, payload
  and HUD checks. The event order matches the preceding failing replay; command
  actions additionally require success, and capture paths identify this run.
- All **ten deferred-delivery rows pass** functional and strict checks on OSMesa
  and System GL. This includes partial-merge and terminal-adoption denial, which
  failed four preceding rows. Failure status, exact bounds, retained leaves,
  policy disable/re-enable and erase are checked separately from the trace.
- Eight Generic shaded/wire cold/warm rows, their eight camera contracts and
  40 HUD pixel checks pass with complete strict traces.
- Four OSMesa/System GL single/quad primitive-edit rows pass 207 events each,
  including ARB, ellipsoid, sketch and BoT checks, with complete strict traces.
- Twenty-four representative original images were inspected: two eager, six
  delivery, eight Generic and eight ARB/ellipsoid editing captures. Existing
  viewport text overlap/clipping remains S5 debt.

The retained preceding traces explain the apparent partial-delivery difference:
SOURCE-022 had zero retained leaf payloads and remained nonterminal before
re-enable, while SOURCE-023 had one leaf and a terminal presentation. The former
did not exercise terminal reopening. Both exposed the attachment change as
`producer_progress` under the old policy revision. This comparison does not
establish why the older run lacked its leaf presentation. The old eager and
adoption failures remain retained negative evidence.

All five clients were rebuilt. Current qged is `9004b9dc...`, libBObol
`8f17b7b8...`, libbv `7120bb05...`, libged `309b5283...` and libqtcad
`95758c5e...`; Obol `54ee73b2...` and OSMesa `097aac27...` are unchanged.
`current-hashes.json` retains full hashes. `runtime.json` records 38 artifacts,
84 distinct resolved libraries and 121 retained physical copies. The evidence
root includes source/header/test diffs, preceding binaries, commands, reports,
images and before/after policy-transition comparisons. The four fixed denial
rows explicitly reach a terminal state with their retained leaf payloads before
re-enable. All tool handles are terminal and the final host-process audit is
empty.

All S1--S6 gates remain incomplete. Prototype/configuration publication, other
writers, Lucy quality, large/cold/multiple-giant resource behavior, complete
client workflows, coherent dependency staging and native/platform qualification
remain open. System GL here uses private-Xvfb llvmpipe. The formal baseline and
strict predicates are unchanged; this repair does not restore the unavailable
historical `/tmp` qualification output.

## September 8: auxiliary source publication

`CXX-SOURCE-023` completes the enclosing nonempty auxiliary-source transaction.
Configuration, line geometry, retained owner metadata, terminal state and
source/parent child paths now publish together. Named source routing identity,
input sensors, semantic overrides and independent nested owners survive. Obsolete
bounds/profile/staging and older worker stamps retire, including replacements
with unchanged configuration. Ordinary controller structural/frame revisions
advance at commit before observers; throwing observers leave a complete edit.
The [conformance record](libbobol_tla_conformance.md#cxx-source-023--auxiliary-source-configuration-and-publication)
owns the API boundaries and remaining custom-callback limits.

New-source allocation denial exposed an Obol field-sensor constructor-cleanup
crash. Attachment now follows successful registration, and registration rolls
back if a custom connection hook rejects it. Regressions cover initial storage,
auditor-list growth, hook rejection, retry and constructor unwinding. The final
test rejects the preceding SOURCE-022 stack for incomplete source publication
and late scene revisions, and rejects the intermediate attachment-only fix for
a dangling auditor during field destruction. Crash stacks and resolved negative
libraries are retained.

The 20 auxiliary scenario definitions run directly and through the controller:
**40 cases, 3,980 allocation positions, 3,782 preserved and 198 committed
outcomes**. Source publication accounts for 18 cases/3,039 positions; line and
removal cover 22/941. All **18 focused checks pass in 53.11 s**. The existing Obol
unit-test executable reports **1,575 passes and one skip** against the final
installed dependency. All five clients were rebuilt. These checks do not qualify
arbitrary reentrant hooks, physical OOM throughout GED, memory/latency magnitude,
full client workflows or native hosts.

Final graphical evidence is under
`.build/obol-qualification/20260908-auxiliary-source/final`:

- Four single/quad primitive-edit rows pass 207 events and complete strict traces
  each, on OSMesa and System GL.
- Ten small deferred-delivery rows pass command, retained-geometry/bounds, policy,
  erase and applicable HUD checks. **Only six strict traces pass.** Partial-merge
  denial and terminal-adoption denial each fail `terminal-reopened-without-revision`
  when LoD is re-enabled, on both backends. The matrix correctly returns failure.
- Eight Generic rows and their eight camera contracts pass, with complete strict
  traces and 40 HUD pixel checks.
- Twenty-two original representative images were inspected: six delivery, eight
  Generic initial-ready and eight ARB/ellipsoid editing captures.

These reconstructed delivery rows start with automatic LoD, explicitly check
command success and retain four actual leaf payloads separately from the
inactive overview entry. Failed delivery preserves one overview, one leaf plus
its overview, or four leaves plus the overview, depending on the failed stage.
Source configuration may become stale after a display-policy change; framebuffer
readiness and preservation of a FAILED source remain separate assertions.

An additional eager-start diagnostic draws with LoD initially off, then enables
LoD. Its settled images are functional, but transitions 20/21 record mask 128
(`NONTERMINAL_WITHOUT_PROGRESS`) before the policy revision advances. The retained
SOURCE-022 stack reproduces those same violations. The preceding adoption-denial
row also reproduces terminal reopening; its partial-delivery row passes, so the
cause of that difference still needs a transition/timing comparison. These are
open acceptance defects at the SOURCE-023 checkpoint. The subsequent
[policy repair](#september-8-live-policy-publication) closes those reproductions
with the predicates unchanged.

The evidence root retains scripts, source/header/test diffs, negative probes,
crash stacks, reports/images and loader resolution. Final qged is `5a09de31...`,
libBObol `e477beba...`, libged `27b88cb4...`, libqtcad `ec6055f5...`, Obol
`54ee73b2...` and OSMesa `097aac27...`; `current-hashes.json` has full hashes.
`runtime.json` identifies 37 runtime artifacts and their 84 distinct resolved
libraries, with retained physical copies. Earlier sibling graphical directories
used the intermediate Obol `14ef8b38...`; use `final/` for SOURCE-023 qualification.
`artifact-relocations.json` resolves the early eager-start reports' moved images.
Private-Xvfb System GL remains CPU llvmpipe. Viewport text overlap/clipping remains
an S5 acceptance limit.

**Evidence availability changed:** prior `/tmp/obol*` and `/tmp/qged*` artifacts
are unavailable. The historical narratives below describe then-observed results;
their claims of retained reports/scripts/images/caches no longer establish
current availability. The formal catalog/baseline and source survive, but raw
TLC output and the old Lucy/planning/frame reproductions must be reconstructed
before release qualification. New evidence stays under the workspace build tree.
Formal relations, quality floors, renderer defaults and the user's proposal are
unchanged. **S1--S6 remain incomplete.**

## September 7: auxiliary line publication and retirement

`CXX-SOURCE-022` repairs named auxiliary line publication, empty auxiliary removal
and complete auxiliary clearing. Geometry, placement and child paths publish
together through child machinery shared with primary publication. Node names,
primitive selection, edit intent and notification settings survive replacement;
primary geometry and compact handles remain retained. Ordinary controller calls
advance scene structure. Nested auxiliaries inherit the parent transform once:
the bounding-box regression now reports minimum `(10, 20, 30)`, replacing the
preceding doubled `(20, 40, 60)`, for both insertion and subsequent replacement.
The [conformance record](libbobol_tla_conformance.md#cxx-source-022--auxiliary-line-publication-and-retirement)
owns the boundary and API semantics.

Eleven scenarios pass **453 allocation positions** (441 preserved, 12 committed).
The unconstrained insertion probe observes complete geometry/metadata throughout
its 12 callbacks; the preceding library exposes incomplete state in a 66-callback
run. Tests cover default commands, display overrides, disabled fields, missing
placement, repeated children, nested-source removal, retry and throwing observers.
All **18 focused checks pass in 38.99 s**, and all five clients were rebuilt.
The final test rejects SOURCE-021 independently for incomplete publication,
missing structural revision advancement and doubled placement.

Ten current delivery rows, eight current Generic rows/camera contracts and four
current single/quad primitive-edit rows pass on both backends. All 22 strict
traces pass; Generic has 40 HUD pixel checks, and each editing row has 207 events.
Twenty-two original representative images were inspected: six delivery, eight
Generic initial-ready and eight ARB/ellipsoid edits. Geometry and handles remain
visible; small-viewport text overlap/clipping remains an S5 acceptance limit.
Private-Xvfb System GL remains CPU llvmpipe, not native GPU qualification.

`/tmp/obol-auxiliary-publication-20260907` retains source/header diffs, tests,
negative probes, reports/images and loader resolution. Current qged is
`310279d29afd...`, libBObol `3a5b1a1b3e94...`, libged `b622f7c26b7c...`,
libqtcad `ec6055f5bbe6...`, Obol `5ddf8e89e12b...` and OSMesa
`097aac27db41...`. Eighteen runtime artifacts and 78 resolved libraries have
retained physical copies. All build/test/GUI handles are terminal and the final
host process audit is empty.

The enclosing nonempty auxiliary source configuration/terminal transaction and
controller revision recovery after observer exceptions remain open. Continue
prototype publication and invalidation afterward. Child/path preparation is not
a complete byte/latency bound; broader provider/traversal ownership, Lucy/large
scenes, full editing and native-host qualification remain open. SOURCE-019 owns
the preceding exhaustive sampled-point sweeps; current routine checks retain
their three selected boundaries. Formal relations, quality floors and renderer
defaults are unchanged. **S1--S6 remain incomplete.**

## September 7: incremental leaf edit publication

`CXX-SOURCE-021` publishes affected leaf geometry, part ownership, revision and
owner metadata together. It fixes indexed selection-path lookup, recomputes exact
bounds when needed, invalidates obsolete source profiles and retains the existing
compiled object and sensors. The
[conformance record](libbobol_tla_conformance.md#cxx-source-021--incremental-leaf-edit-publication)
owns the boundary. No full source/index copy or persistent mode was added.

Twelve scenarios pass **444 allocation positions** (428 preserved, 16 committed).
All **18 focused checks pass in 33.91 s**; all five clients were rebuilt.
The scaling check records **56 preparation allocations** with both 32 and 4,096
unrelated occurrences. Indexed-path editing updates both repeated occurrences,
and exact bounds grow and shrink. Four scenarios named provider failure exercise
missing-object rejection specifically; they do not cover every provider failure.
Allocation counts do not establish total byte or latency bounds, and replacing an
extremal occurrence can require scanning the index. SOURCE-019 retains exhaustive
sampled-point sweeps; current routine cases cover their three selected boundaries.

Ten current delivery rows pass functional, strict and applicable HUD checks.
Four current primitive-edit rows pass 207 events and strict traces each on OSMesa
and System GL, in single/quad layouts. Fourteen original captures were inspected:
eight ARB/ellipsoid edits and six healthy/denied-adoption delivery images. Geometry
and handles remain visible, and denied adoption retains the drawing with a
terminal error indicator. Viewport text overlap/clipping remains an S5 acceptance
limit. Generic was not rerun; SOURCE-019 owns its preceding eight-row matrix.

`/tmp/obol-leaf-edit-20260907` retains source diffs, tests, negative probes, runtime
hashes, loader resolution and graphical reports/images. Current libBObol is
`2a2a1117bfe3...`; qged remains `dd29cfbd6fda...`, Obol `5ddf8e89e12b...`,
and OSMesa `097aac27db41...`. Eighteen runtime artifacts and all 78 resolved
libraries have retained physical copies. Final tests reject SOURCE-020 for stale
source profiles and indexed-path lookup; the earlier observer counterexample is
also retained. All build/test/GUI handles are terminal and the host process audit
is empty.

Auxiliary/prototype publication and invalidation are next. Other cache operations,
broader provider/traversal cleanup, resource bounds, Lucy/large scenes, full editing
and native-host qualification remain open. Formal relations, quality floors and
renderer defaults are unchanged. **S1--S6 remain incomplete.**

## September 7: combination edit publication

`CXX-SOURCE-020` repairs combination revision/geometry publication and stale exact
bounds. It also fixes a draw-cache deadlock exposed by allocation failure during
BoT reuse. The [conformance record](libbobol_tla_conformance.md#cxx-source-020--combination-edit-publication-and-cache-lock-recovery)
owns the implementation boundary. No controller, scene-sized mirror or workload
mode was added; the existing prepared source publication is shared.

Eight production scenarios pass **985 allocation positions** (977 preserved,
8 committed). These are measured positions, dependent on cache/test ordering.
All **18 focused checks pass in 34.51 s**, including the existing cached/external
publication, source/service, prototype, GED and qtcad checks. The sampled-point
cases retain their routine three-boundary coverage; SOURCE-019 owns their earlier
exhaustive sweeps. All five clients were rebuilt.

Four current primitive-edit replays pass on OSMesa and System GL in single/quad
layouts, with 207 events and complete strict traces each. Eight representative
ARB/ellipsoid captures were inspected. Geometry and handles remain visible;
viewport text crowding/clipping remains an S5 acceptance limit. These replays
exercise integration; the direct production test owns the combination mutation
and failure-path assertions. SOURCE-019 retains the preceding Generic and ten-row
delivery matrix; those rows were not rerun for this checkpoint.

`/tmp/obol-combination-edit-20260907` retains patches, current/preceding libraries,
tests, timeout/stack evidence, GUI scripts/reports/images and loader resolution.
Current libBObol is `30ebfc0d981e...`; qged is
`dd29cfbd6fda...`, Obol `5ddf8e89e12b...` and OSMesa
`097aac27db41...`. Eighteen runtime artifacts and 78 resolved libraries
have retained physical copies. The final test rejects SOURCE-019's stale bounds.
The intermediate BoT timeout is a diagnosed semaphore leak, distinct from the
preceding sampled-point sweep's allocator contention. Initial fixture and compile
corrections remain in their preliminary logs.

Leaf edit publication, other cache operations, broader provider/traversal cleanup,
resource bounds, Lucy/large scenes, full editing and native-host qualification
remain open. Formal relations, quality floors and renderer defaults are unchanged.
**S1--S6 remain incomplete.**

## September 7: live cached realization

`CXX-SOURCE-019` closes the live wire/mesh and evaluated representation publication
boundary. Public APIs and the realization action construct an owned private
candidate and publish it with the existing prepared source/child machinery.
Workers explicitly retain private construction and incremental stream handoff.
The repair removes the cache preservation flag and repeated private terminal
assignments, scopes template/evaluated temporary ownership, and publishes provider
failure coherently while retaining the last drawing. The
[conformance entry](libbobol_tla_conformance.md#cxx-source-019--live-cached-realization-exposes-private-construction)
owns the operation boundary and its limits.

The retained public probes reproduce incomplete observations against SOURCE-018:
12/14 wire callbacks, 12/14 shaded-mesh callbacks, 11/48 evaluated-wire callbacks
and 14/48 evaluated-points callbacks. Current probes observe **zero incomplete
callbacks**, with ten ordinary wire/mesh notifications and eleven for each
evaluated representation. The final production test independently rejects the
preceding libBObol's terminal/child state. Allocation injection also reproduced
`SoSFEnum::setEnums` freeing its current mappings before replacement allocation;
the repaired dependency owns both arrays before replacing either. The final enum
test rejects the preceding Obol, tests both allocation denials, and verifies
complete replacement and self-assignment.

The final cached qualification covers **18 scenarios and 3,691 allocation
positions** (3,469 preserved, 222 committed). It includes public/action APIs,
explicit NULL/non-NULL streams, repeated assemblies, compact and compiled prior
drawings, selection/appearance overrides, auxiliary/nested owners, placement,
disabled fields, provider errors, child/path retirement, retry and throwing
callbacks. The compiled cases account for 250 positions. The separate enum sweep
passes three positions. `publication-metrics.json` preserves the per-case counts.
The two named evaluated-points cases run their full sweeps separately; routine
CTest checks three allocation boundaries for each while fully sweeping the
other cases. Pseudo-random sampled coordinates need not equal another call;
observers check complete triangles/normals and exact terminal/owner state.

All **18 focused checks pass in 44.32 s** and **304 existing Obol checks pass**
against the installed library. Ten small delivery replays, eight Generic Twin
rows/camera contracts and four single/quad primitive-edit replays pass on both
backends. Generic has complete strict traces and **42 HUD pixel checks**; all
four edit traces pass, with 207 events each. Thirty original images were inspected.
Small-pane viewport text retains its S5 visual/layout acceptance limit.

`/tmp/obol-live-cached-realization-20260907` retains source patches, counterexamples,
allocator/integration logs, event scripts, reports, images and loader resolution.
SHA-256 prefixes: qged `dd29cfbd6fda`, libBObol
`db5a61d2c398`, libged `420d171b79ab`, libqtcad
`ec6055f5bbe6`, Obol `5ddf8e89e12b`, OSMesa
`097aac27db41`. Eighteen runtime artifacts and all 78 resolved
libraries have retained physical copies. Obol was installed with the explicit
workspace prefix; OSMesa was restored to the matching artifact after configuration
restaged the external bundle.

Preliminary test issues are retained separately: the initial mesh fixture still
selected wire display, exact comparison was inappropriate for pseudo-random
sampled coordinates, and the action fixture had outstanding work on its nested
source. Stack samples confirmed the long point sweep advanced through allocator
contention; they do not establish a deadlock. The enum double-free and live
publication counterexamples are repaired production defects.

This qualifies owning-thread publication under the tested C++ allocation and
observer failures. It does not establish complete C allocator/provider/traversal
cleanup or geometry/snapshot memory bounds. Generic's pose checks use the
existing soft-deadline contract; private-Xvfb System GL is CPU rendering.
Renderer defaults, quality floors and formal relations remain unchanged.
Compact edit refresh is next. Other writers, Lucy/large scenes, full editing
acceptance, native hosts and shared-stack sanitizers remain required.
**S1--S6 remain incomplete.**

## September 7: external geometry publication

`CXX-SOURCE-018` repairs external line, point, triangle, annotation and primitive
publication. Geometry is prepared privately, and source/owned auxiliary metadata,
placement, bounds, child/path replacement and compact/compiled retirement commit
before callbacks. Clear retains existing primary shape objects. Nested sources
keep their own state. The submodel adapter stages its replacement leaves,
aggregates transformed bounds and preserves repeated occurrence identity.
The [conformance boundary](libbobol_tla_conformance.md#cxx-source-018--external-replacement-discards-the-preceding-drawing)
defines the qualified operations and remaining limits.

The production line probe now retains the old drawing on allocation denial;
ordinary replacement has **nine callbacks and zero incomplete observations**.
The final test fails against the preceding retained libBObol. The expanded
sweeps pass **25 external scenarios, 2,591 allocation positions** (2,542 preserved,
49 committed), plus **21 prepared-child replacement positions** (20 preserved,
one committed). They cover immediate source/shape/path observers, throwing
callbacks and retry, shared geometry, compact handles, compiled children,
placement, annotation text/precision, empty input and actual primitive adapters.
Invalid inputs do not allocate or notify. Temporary plot blocks are counted
across conversion denial. Existing installation/terminal sweeps pass 1,205
positions; merge passes 618, realization state 13, child retirement nine and
append five. Counts describe this final test process and its measured allocation
order, not universal constants.

This also fixes float-only annotation points being cleared when optional double
coordinates are absent. Repeated submodel instances remain distinct with
combination instance IDs enabled or disabled. The broader GED test exposed an
existing fallback occurrence visitor that stopped at the first unrelated record;
its corrected continuation convention keeps a cleared source addressable and
allows subsequent highlighting/readback.

All **18 focused checks pass in 20.02 s**; **304 existing Obol checks pass** against
the installed library. All **ten small delivery replays**, **eight Generic Twin
rows/camera contracts** and **four primitive-edit replays** pass on OSMesa and
System GL. The edit rows preserve the existing 199 functional events and add
four frame waits/captures, for 207 events each, in single and quad layouts.
Generic strict traces are complete and **35 HUD pixel checks** pass. Thirty
original images were inspected: fourteen drawing views and sixteen editing
captures. Editing handles and geometry are present; small-pane text still
crowds/clips, and these captures are not pixel-equivalence or full visual-layout
qualification.

`/tmp/obol-external-publication-20260907` retains patches, allocator and integration
logs, negative evidence, event scripts, reports, images, executables and loader
resolution. SHA-256 prefixes: qged `dd29cfbd6fda`, libBObol `dd504a447206`,
libged `ab0093c68883`, libqtcad `ec6055f5bbe6`, Obol `464c91b70e6a`, OSMesa
`097aac27db41`. Eighteen runtime artifacts and all 78 resolved libraries have
retained physical copies. Matching Obol provides the new prepared replacement
API without changing child-list layout; it was installed with the explicit
workspace prefix.

Retained preliminary failures are classified separately: a faulty path target
in the first compact fixture, an auxiliary-bearing fixture ineligible for compiled
presentation, and an edit screenshot requested before its frame existed are test
setup failures. The initial GED fallback failure and annotation/submodel failures
are repaired production defects. Generic wire rows still exercise the existing
soft-deadline contract, including completed peaks of 57.31/56.74 ms and explicit
interruptions around 50 ms; no hard-deadline claim follows. Renderer defaults,
quality floors and the formal baseline are unchanged. Private-Xvfb System GL
remains CPU rendering. Live cached realization is next; other writers, resource
bounds, actual snapshot cleanup, Lucy/large scenes, full editing acceptance,
native hosts and shared-stack sanitizers remain open. **S1--S6 are incomplete.**

## September 7: child append ownership

`CXX-SOURCE-017` repairs two allocation-ownership defects in the Obol dependency.
Base-list append/insertion acquired references before storage growth; child
append also installed its parent auditor before growth. Failure could leak a
reference and route notifications through a parent with no corresponding child.
The exact final tests reject both defects with the preceding retained Obol and
pass with the repaired library. `/tmp/obol-child-append-20260907` retains the
source changes, negative probes, tests, binaries and loader resolution.

Child append now reserves storage before connecting the parent auditor, then
publishes the slot and reference without allocating. Base-list append and
insertion acquire references only after storage succeeds. `SbPList::reserve`
changes capacity alone and retains the existing growth policy. List layout and
ordinary append notification order are unchanged. The
[publication inventory](libbobol_tla_conformance.md#cxx-source-017--failed-append-leaves-unowned-references-or-parent-links)
defines the qualified boundary.

The append sweep passes **five allocation positions**, three preserving old
ownership and two leaving complete commits. It covers new and shared children,
reference counts, parent notification routing, existing paths, retry and throwing
observers. The separate base-list check covers both append and insertion.
Existing sweeps also pass: 1,201 installation/terminal positions, 618 merge
positions, 13 realization-state positions and nine child-retirement positions.
All **18 focused checks pass in 15.82 s**, and the retained existing Obol test
executable passes **304 checks** against the installed library.

All **ten small delivery replays** pass functional, strict and applicable HUD
checks. All **eight Generic Twin rows and camera contracts** pass on OSMesa and
System GL, with complete strict traces and **37 HUD pixel checks**. Fourteen
original images were inspected: six healthy/denied-adoption views and eight
Generic initial-ready views. `generic/contract-audit.json` retains the actual
pose costs, including interrupted warm OSMesa wire frames and System wire
completed peaks of 55.59/64.54 ms. These are the existing soft-deadline contracts;
this checkpoint does not establish a hard frame deadline or close Lucy failures.

Final SHA-256 prefixes: qged `dd29cfbd6fda`, libBObol `7e8f7bcaf707`, libged
`a9dd4bcba3eb`, libqtcad `ec6055f5bbe6`, Obol `768872853d63`, OSMesa
`097aac27db41`. The drawing/client targets were rebuilt against matching Obol,
installed with the explicit workspace prefix. Eighteen runtime artifacts and
all 78 libraries resolved for the qualified commands have retained physical
copies. The Obol core-test command explicitly selects `.build/lib`; its loader
record uses the same environment. The initial preload attempt lacked this
search path and is retained as `negative-preload-invalid-loader.json`, not defect
evidence. Renderer defaults, quality floors, formal relations and the accepted
baseline are unchanged. Private-Xvfb System GL remains CPU rendering.

The [direct-writer inventory](libbobol_tla_conformance.md#direct-realization-writer-inventory)
now distinguishes private workers from live cached realization, external
publication, compact edits, auxiliary/prototype publication and invalidation.
The same production helpers serve both worker and live sources. At that checkpoint, an external
line-set probe lost the old geometry under allocation denial, and 53 of
91 immediate callbacks saw incomplete replacement during an ordinary successful
call. The subsequent [external publication repair](#september-7-external-geometry-publication)
closes that `CXX-SOURCE-018` boundary; the original counterexample is retained. Individual child insertion/removal/replacement, traversal cleanup,
other source/service/presentation owners, resource/geometry bounds, editing,
large scenes, native hosts and full shared-stack sanitizers remain required.
**S1--S6 remain incomplete.**

## September 7: realization-state publication and auditor allocation

`CXX-SOURCE-016` repairs the status setter used by the scene controller and GED
failed stream delivery. The initial regression fails against the preceding
candidate without allocation denial: an immediate observer sees a partially
updated source or owned shape. The expanded fixture also rejects a parent
status change overwriting a nested source's shape metadata.
`/tmp/obol-realization-state-20260907` retains the negative probes, source
patches, final binaries and reports. Its
[publication inventory](libbobol_tla_conformance.md#cxx-source-016--realization-status-publishes-before-owned-shape-metadata)
defines the changed boundary.

The setter prepares the eight realization fields and changed diagnostic storage
for the source and each owned wire/mesh shape. It visits shared legacy nodes
once, stops at nested source owners and retains prepared shapes across callbacks.
Commit changes all records quietly and invalidates compiled/batch eligibility
before notification. It no longer rebuilds compact display records or applies
placement/material/display state during a status update. The fixed-size field
notification guard is shared with terminal adoption; all participating flags
are restored before callbacks, and the first callback exception propagates
after the remaining notification attempts. An unchanged request allocates
nothing and sends no notification. Temporary storage scales with legacy nodes
and changed diagnostics, independently of the compact occurrence population.

The full regression exposed a dependency failure hidden by the isolated case:
allocation denial while snapshotting a field-auditor list leaked its recursive
notification lock. The SoDB counter returned to zero, but the next owner thread
blocked. The retained symbolized GDB stack identifies that lock; the stalled
test was terminated before rebuilding. Obol now scopes auditor-list locking
through mutation and snapshot construction. A second probe showed `SbPList`
publishing larger capacity before allocating its replacement buffer, allowing
a retry to trust nonexistent storage. Growth now commits capacity after allocation
and copying succeed. The final isolated probes fail or time out with the
preceding Obol and pass against the repaired library. Existing list layout,
growth policy and notification ordering are unchanged.

The status sweep passes **13 allocation positions**, 11 preserving the old
records and two leaving complete commits. It covers source/shape observers,
shared children, nested owners, compact-record preservation, notification flags
and retry. Additional cases pass for explicit revisions, invalid/unrealized
requests, disabled notifications, throwing observers and callback-driven child
removal. The field-auditor regression proves a later thread can notify; the
pointer-list regression proves denied growth preserves capacity and retry.
Existing coverage also passes: 1,201 installation/terminal positions, 618 merge
positions and nine child-retirement positions. All **18 focused checks pass
in 16.38 s**. The retained existing Obol test executable passes **304 string,
field, sensor, path, list and node checks** against the current library.

All **ten small delivery rows** pass functional, strict-trace and applicable HUD
checks on OSMesa and System GL. All **eight Generic Twin rows and camera
contracts** pass with complete strict traces and **39 HUD pixel checks**.
Fourteen original images were inspected: healthy delivery with LoD off/on,
denied adoption on both backends and eight Generic initial-ready views. Geometry
and placement agree across the paired runs; denied adoption retains four leaves
and the marker with the terminal error display. `generic/contract-audit.json`
records soft-deadline peaks and interruptions. These passes do not establish
hard frame deadlines or resolve the retained Lucy quality failures.

Final SHA-256 prefixes are qged `d04c1dc3aa45`, libBObol `df11a220e5f9`, libged
`88f419c0878e`, libqtcad `941db3374b80`, Obol `66ec2c649eea` and OSMesa
`097aac27db41`. The dependency and drawing/client targets were rebuilt; Obol
was installed with the explicit workspace prefix. Eighteen runtime artifacts
and all 78 `ldd`-resolved libraries have retained physical copies. The final audit
removed one trailing blank line from the dependency source; the qualified file
and exact difference are retained, with no executable token changes. Renderer defaults, quality floors,
formal relations and the accepted baseline are unchanged. Private-Xvfb System
GL is CPU rendering, not native GPU qualification.

`mark_source_realized_current`, other direct realization success/failure writers,
legacy child mutation, traversal cleanup, sparse setters/journal retention,
field connection/engine and delayed-sensor lifetimes, and throwing deletion or
custom cleanup remain open. This setter repair does not qualify the full GED
diagnostic path under persistent physical OOM. Completed-result/geometry bounds,
actual snapshot-file cleanup, Lucy/timing reproductions, large scenes, editing,
native hosts and shared-stack sanitizers remain required. **S1--S6 remain
incomplete.**

## September 7: child, path and auditor retirement

`CXX-SOURCE-015` closes a gap in the preceding terminal-adoption checkpoint:
direct path observers could see a source with only some obsolete children
removed. All four strengthened terminal scenarios fail against the preceding
libBObol, including ordinary publication without allocation denial.
`/tmp/obol-child-retirement-20260907` retains those negative results, the
independent child-removal probe, source patches, binaries and loader resolution.
The [publication inventory](libbobol_tla_conformance.md#cxx-source-015--path-callbacks-observe-partially-retired-source-children)
defines this boundary and its remaining callers.

Obol's new prepared `SoChildList::Removal` owns removal indices and references
to retired nodes and affected paths. Preparation validates indices and identifies
shared children whose parent link must survive. Commit compacts the child list
once and updates all affected paths before notification. Temporary storage scales
with removed children and registered path auditors; it does not copy the scene
graph or path chains. Terminal adoption prepares this operation before mutating
the source and commits it alongside the existing prepared realization. Path and
parent notifications each get an attempt even when another observer throws; the
first exception propagates afterward. Existing individual `remove`/`truncate`
callers keep their previous notification order and require separate review.

The lifetime probe also reproduced an allocation-denial abort during node
release, a crash when a deletion callback destroyed another pending sensor,
and unrelated immediate callbacks dispatched by path/group destruction.
Node destruction now retires sensors directly from the live auditor tree without
allocating a snapshot or retaining pointers to sensors another callback can
delete. Node-sensor detachment cancels its queued callback. Path and child-list
destructors retire links quietly, leaving unrelated queued work for normal
dispatch. Immediate-queue processing attempts the remaining callbacks within
its existing rescheduling limit and restores its locks before rethrowing the
first exception. Negative libraries and crash/probe records are retained.

The final four terminal sweeps pass **217 allocation positions**: 205 preserve
the old source and 12 leave a complete commit. They check terminal fields,
staging ownership, auxiliary-child retention, direct path observers and retry.
All eleven installation/terminal sweeps pass **1,201 positions**, with 1,182
preserved and 19 committed outcomes; ten merge sweeps pass 618 positions.
The separate child-retirement sweep passes **nine positions**, seven preserved
and two committed. Its additional checks cover abandoned/invalid preparation,
duplicate shared children, multiple affected paths, throwing path observers,
eight attached node sensors, retired subtrees and recovery. Queue detachment,
quiet path/group destruction and deletion of another pending auditor also pass.
All **18 focused checks pass in 16.10 s**. The retained existing Obol test
executable passes **304 string, field, sensor, path, list and node checks**
against the final workspace library; this is not a rebuild of that test suite.

Final graphical evidence is under **`qualification-final/`**. All **ten small
delivery rows** pass functional, strict-trace and applicable HUD checks on OSMesa
and System GL. All **eight Generic Twin rows and camera contracts** pass, with
complete strict traces and **37 HUD pixel checks**. Fourteen original images
were inspected: healthy delivery with LoD off/on and denied adoption on both
backends, plus eight Generic initial-ready views. The expected four leaves and
marker survive denied adoption with the terminal error display; Generic
silhouettes and placement agree across cold/warm and backend pairs. The pose
audit records completed peaks above the soft target in the System wire cases;
passing the retained-cut contract is not a hard frame-time guarantee.
The artifact-root `generic/` and `replays-final/` directories retain an earlier
passing run before the final child-list destructor repair. Their names remain
unchanged because reports embed absolute image paths.

Final SHA-256 prefixes are qged `b31a7b10022f`, libBObol `05e67011a480`, libged
`05d0beffb8a8`, libqtcad `941db3374b80`, Obol `412284ca34a6` and OSMesa
`097aac27db41`. Eighteen runtime artifacts and all 78 `ldd`-resolved libraries
have retained physical copies. The new child-removal API requires the matching
Obol dependency; existing child-list layout is unchanged. Obol was installed
with the explicit workspace prefix and the drawing targets were rebuilt.
One obsolete comment describing the deleted auditor snapshot helper was removed
after qualification; its exact difference and qualified source are retained.
No executable code changed after the final tests.

Other realization writers, legacy child-mutation callers, traversal cleanup,
sparse setters/journal retention, connected-field/engine and delayed-sensor
exceptions remain open. Throwing deletion callbacks and custom node cleanup
are outside this tested destruction boundary. The existing GED diagnostic path
is not a persistent physical-OOM guarantee. Completed-result/geometry bounds,
actual snapshot-file cleanup, Lucy/timing reproductions, large scenes, editing,
native hosts and full shared-stack sanitizers remain required. Renderer defaults,
quality floors and formal relations are unchanged. Private-Xvfb System GL is
CPU rendering. **S1--S6 remain incomplete.**

## September 7: terminal adoption and notification recovery

`CXX-SOURCE-014` strengthens the preceding registry-only adoption boundary. The
corrected persistent allocation probe finds six mixed replacements and nine
mixed streamed adoptions against `CXX-SOURCE-013`; all 81 positions pass with the
repair. The four strengthened terminal scenarios also fail against that preceding
library, including hidden overview state which a later visibility change could
revive. `/tmp/obol-terminal-adoption-20260907` retains the probes, negative results,
source patches, runtime artifacts and loader resolution. The corrected probe
accepts the intentionally revoked LoD requests of a preserved stale source; the
preliminary positive probe's request-aggregate failures were an oracle error.

Adoption prepares the identity string, bounds and field notification state before
changing the live source. The streamed branch reserves sparse retirement storage
before retiring visible or temporarily hidden authored overviews; replacement uses
the owned-index installer. Staging, current request stamps and terminal metadata
commit before ordinary child/path retirement and final field/node notification.
The fixed fourteen-field guard restores all notification settings even on failure.
The [publication inventory](libbobol_tla_conformance.md#cxx-source-014--terminal-adoption-exposes-mixed-realization-state)
and `publication-ownership.json` record the exact owners and limits.

Immediate-observer injection exposed a dependency defect: failed notification
left a global transaction open or a field recursion flag set, suppressing later
callbacks. Obol now owns notification cleanup through scope guards, releases the
notification lock after throwing callbacks, and releases immediate-queue insertion
locks after allocation failure. Cleanup preserves the original exception without
allocating diagnostics. `SoSFString` gains an rvalue setter so a prepared identity
can transfer into a quiet field without copying; field layout is unchanged.
Both notification regressions fail with the preceding Obol library and pass with
the repair, including a subsequent field/node writer on another thread.

The final four terminal sweeps pass **185 allocation positions**: 173 preserve the
old source and 12 leave a complete committed source. They include metadata/bounds,
renderer records and memberships, publication identities, immediate observers,
restored notification flags, ordinary child-path retirement, staged import
ownership/claim/release, hidden overview non-resurrection and retry. All eleven
installation/terminal sweeps pass **1,169 positions**, with 1,150 preserved and 19
committed outcomes. The ten existing merge sweeps pass all 618 positions.
All **18 focused checks pass in 15.89 s**. A retained Obol unit-test executable
also passes 154 existing string, field and sensor checks against the repaired
workspace library; its executable and loader resolution are retained.

GED now catches denied terminal adoption through its existing failed-delivery
path. All **ten** small delivery rows pass functional, strict-trace and applicable
HUD checks on OSMesa and System GL. The two new rows retain all four delivered
leaves and certified whole-root overview, report FAILED/terminal error, settle
across LoD policies and permit erase. All **eight** cold/warm shaded/wire Generic
Twin rows and camera contracts pass, with complete strict traces and **36 HUD
pixel checks**. Fourteen original images were inspected: eight Generic
initial-ready views, healthy delivery with LoD off/on, and denied-adoption
terminal displays on both backends. The latter show four leaves and the marker,
a full red bar and the geometry-error annotation. Pose durations and interruptions
are recorded in `generic/contract-audit.json`; these are soft deadline contracts.

Final SHA-256 prefixes are qged `332495d5a4ff`, libBObol `05e06cc9ebcf`, libged
`e2d6bc61797a`, libqtcad `965f59158f0a`, Obol `246c8e0360b1` and OSMesa
`097aac27db41`. Obol and drawing targets were rebuilt; the dependency was installed
with an explicit workspace prefix. Eighteen runtime artifacts and all 78
`ldd`-resolved libraries have retained physical copies. Renderer defaults, quality
floors, formal relations and the accepted baseline did not change. Private-Xvfb
System GL remains CPU software rendering, not native GPU qualification.

Direct path-auditor exceptions during child-list mutation, other realization
writers (`setRealizationState`, `mark_source_realized_current`), traversal cleanup,
connected-field/engine and delayed-sensor exception paths remain outside this
qualified boundary. GED's diagnostic publication is not a persistent physical-OOM
survival guarantee. Real completed-result/geometry bounds and snapshot cleanup,
Lucy/timing reproductions, large scenes, editing, native clients and full
shared-stack sanitizers remain open. **S1--S6 are incomplete.**

## September 7: owned source index installation

`CXX-SOURCE-013` repairs full-index candidate ownership and installation. The
initial persistent allocation probe finds **140 inconsistent positions out of
168**: leaked unpublished geometry, obsolete source bounds and incomplete
presentation. The expanded regression fails against the preceding library in
all nine installation scenarios; detached adoption invalidates an old live
handle at allocation zero. `/tmp/obol-source-installation-20260907` retains the
probe, negative results, source patches, binaries and loader resolution.

Setters, tree-walk data and detached transfer now own their candidate explicitly.
Runtime-state transfer, lookups, indexed visibility/selection and renderer
memberships are prepared before installation. Commit transfers the registry and
handle-source identity together, publishes its revisions and bounds, and releases
old history before field notification. The existing exact whole-root bounds are
preserved. Full replacement uses the existing path indexes and two temporary
scalar masks; it does not scan every selection path for every occurrence.
The [ownership inventory](libbobol_tla_conformance.md#cxx-source-013--failed-full-index-preparation-leaks-geometry-or-invalidates-live-handles)
and `publication-ownership.json` identify the qualified boundary and its limits.

The final nine installation sweeps pass **1,147 allocation positions**: 1,126
preserve the old registry and 21 leave a complete replacement. Checks cover
public lookups, old handles and publication identity, whole-population state,
source contracts/bounds, exact-bound preservation, renderer geometry/placement/
opacity/memberships, uncommitted-geometry release and retry. The 21 committed
cases include failures in later detached realization-field writes; the registry
oracle does not qualify those fields as a terminal transaction. The ten existing
merge sweeps also pass their 618 positions, 311 partial prefixes and 331 full
refreshes. All **18 focused checks pass in 15.79 s**.

All eight small delivery rows pass functional, strict-trace and applicable HUD
checks on OSMesa and System GL, including full/reserve/partial denial, healthy
delivery, LoD changes and erase. All eight cold/warm shaded/wire Generic Twin rows
and camera contracts pass. Their strict traces are complete, without dropped or
truncated transitions, and **41 HUD pixel checks** pass. Fourteen original images
were inspected: eight Generic initial-ready views, healthy delivery with LoD
off/on on both backends, and both partial-denial terminal displays. The latter
retain one leaf, exact whole-root overview, full red bar and error annotation.
`generic/contract-audit.json` records durations and deadline interruptions;
passing pose contracts do not claim every completed frame met its soft target.

The final libBObol SHA prefix is `ced3af429fc2`; qged remains `7307d903b15d`,
Obol `626f851905dc` and OSMesa `097aac27db41`. All drawing targets were rebuilt
and final checks used this stack. Seventeen binaries and 78 `ldd`-resolved
libraries have retained physical copies. This turn changed no upstream Obol
source, renderer default, quality floor, formal relation or accepted baseline.
Private-Xvfb System GL is CPU software rendering, not native GPU qualification.

Candidate ownership and registry consistency are qualified by these cases.
Stream-preserving overview retirement, child retirement, final realization
fields/staging transfer and traversal resource cleanup remain outside this
transaction. Other source setters/sparse writers, real completed-result and
geometry-memory bounds, snapshot cleanup, Lucy/timing failures, large-scene and
editing cases, native clients and full shared-stack sanitizers remain open.
S1--S6 are incomplete.

## September 7: compact merge publication

`CXX-SOURCE-012` repairs a source-publication failure below the previous
pre-merge GED fault. A three-leaf persistent allocation sweep finds 111 of 143
positions leaving incomplete lookup state or missing publication evidence.
The expanded regression fails against the predecessor libBObol in all ten
scenarios, including retained unpublished parts, partially replaced requests,
and root-leaf publication which ignores an active visibility frontier.
`/tmp/obol-compact-merge-allocation-20260907` retains those negative results.

Each occurrence now prepares its entry, part identity, lookup slots, memberships
and presentation rules before committing. It stages one occurrence and shares
immutable arrays. Append's deque insertion is the last allocating operation;
replacement uses a checked nonthrowing move at the stable entry address. A
later failure publishes the complete prefix, including source-contract stamps,
overview state and bounds. Allocation-failed CAD/LoD/visibility journals clear
sparse history and advance their floor. Consumers explicitly refresh from the
source. Bounds fields become coherent before observers are notified. The
[ownership inventory](libbobol_tla_conformance.md#cxx-source-012--allocation-failure-exposes-an-incomplete-compact-merge)
and `publication-ownership.json` classify these writers and limits.

The final production regression passes **618 allocation positions across ten
scenarios**, including 311 partial prefixes and 331 authoritative-refresh cases.
It covers empty/existing sources, short paths and shared leaf keys, active
presentation rules, geometry replacement, same-tier source enrichment,
evolving overviews and root-leaf replacement. Checks include public lookup,
old-or-complete entry state, exact renderer records/geometry/memberships,
source aggregates and bounds, release of unpublished geometry and successful
retry. Allocation failure remains enabled during unwinding. The three-leaf
initial repaired probe separately passes all 122 positions; it precedes the
expanded bounds/ownership regression and is not the final qualification scope.

GED's new controlled fault waits for its small worker to complete, then throws
inside a multi-leaf merge after the first leaf commits. The source remains FAILED
with one usable leaf and exact whole-root overview; cancelled delivery cannot
adopt the completed worker's full registry. Eighteen final focused checks pass
in **16.05 s**, including this integration and the allocator regression.

Obol was rebuilt from the sibling checkout with the defaulted `SbString` move
constructor and assignment, installed with an explicit workspace prefix, and the
drawing binaries relinked. Configuration again restaged the incompatible OSMesa
bundle; the matching library was restored after an empty host process audit.
Current SHA-256 prefixes are qged `7307d903b15d`, libBObol `18b7b82e7ed6`, libged
`5e36e6a0cb2f`, libqtcad `1c0c46c5e26c`, Obol `626f851905dc` and OSMesa
`097aac27db41`. Seventeen binaries and all 78 `ldd`-resolved libraries have
retained physical copies. The earlier small replays and focused run before the
Obol rebuild are retained separately and do not qualify these final binaries.

All eight final small 25-event rows pass functional, strict-trace and applicable
HUD checks: merge denial, reserve denial, partial merge denial and healthy
delivery on each backend. OSMesa durations are 1,114/1,122/1,127/1,137 ms and
System durations 1,048/1,055/1,056/1,070 ms. Partial-denial images retain one leaf
inside the whole-target overview; with LoD enabled, the full red bar and
`View incomplete: 1 geometry error` annotation report the failure. Healthy
images show four leaves. Policy changes and erase pass.

The final Generic Twin matrix passes all eight cold/warm shaded/wire rows and
camera contracts on both backends. Eight strict traces are complete without
dropped/truncated transitions, and **40 HUD pixel checks** pass. Durations are
10,044/9,764/10,031/9,657/9,939/10,038/10,555/10,304 ms; transition counts are
444/411/473/462/306/276/360/333. System cold wire exercises deadline recovery at
50.23 and 50.16 ms against 50 ms and passes the existing pose contract. Eighteen
original images were inspected: ten small-delivery views and eight Generic
initial-ready views. `generic/contract-audit.json`, `replays-final/assertions.json`
and `visual-inspection.json` retain those observations.

This closes the tested compact-merge consistency boundary. It does not establish
physical OOM survival, real-geometry cancellation/memory limits, full source
installation/setter exception safety, every sparse writer, or native/client/
shared-stack sanitizer qualification. Private Xvfb System GL remains CPU
llvmpipe. Lucy and the retained planning/frame/quality failures were not rerun.
No rendering quality floor, control trace predicate or formal relation changed.
**S1--S6 remain incomplete.**

## September 7: pose deadline observation

`CXX-TIMING-002` explains two retained System warm-wire assertion failures. In
the queued-cancellation run, completed renders peak at 48.465424 ms, below the
52.5 ms tolerant target, but an intervening frame interrupts at 50.140745 ms
against the 50 ms hard deadline. Its counter advances from 3 to 4 and transition
277 records `render-deadline`. The earlier passive-checkpoint run also advances
from 3 to 4 at 50.140426 ms. The former validator ignored interrupted frames,
so it incorrectly rejected the recovery ceiling on a later fast frame.

The matrix now calls `lod_pose_contract.jq`. A responsive System GL pose keeps
its ceiling absent and at least 95% of its presented population. The existing
5% completed-frame timing tolerance and 1.01 pixel limit are unchanged. A new
deadline interruption during the gesture is pressure evidence even when the
following completed render is fast. Historical interruption counts and evidence
after the held checkpoint cannot justify an earlier ceiling. Counters must be
nonnegative integers and monotonic; new interruptions require a positive hard
deadline and a measured duration at or beyond it.

Population retention now uses exact presented faces in shaded mode and lines
in wire mode. The old deep diagnostic face counter was zero in both tested modes
when deep collection was disabled, making its comparison vacuous. Missing,
inexact or empty baseline population cannot certify responsive retention.
OSMesa retains the previous predicate's scope; other visual and trace gates
continue to apply to both backends.

Evidence is at `/tmp/obol-pose-deadline-contract-20260907`. The extracted old
predicate fails both retained warm-wire reports; the new filter passes those
reports and all eight reports from the queued-cancellation matrix. The original
reports, validation logs and failed matrix summaries are preserved. These are
separate offline verdicts in `offline/audit.json`, not relabeled full-row passes.
Twenty-eight synthetic contract cases pass, including stale/later interruptions,
missing/decreasing/fractional counters, inconsistent deadlines, uncorrected
pressure, shaded/wire population loss and exact threshold boundaries.

Normal CMake registration/build passes; 17 focused checks pass in 14.14 s,
including the new `qged_lod_pose_contract`. Configuration again restages the
incompatible OSMesa library. After an empty host process audit, the matching
library is restored into the workspace before graphical qualification. Sixteen
binaries and all 78 `ldd`-resolved libraries are physically retained. SHA-256
prefixes are qged `c2e92c0100bd`, libBObol `2758c68c77a5`, libged `69f2749c8e88`,
libqtcad `aecad5fcc488`, Obol `f281b59bdbd2` and OSMesa `097aac27db41`.

The fresh Generic Twin matrix passes all eight cold/warm shaded/wire rows and
all eight camera contracts on OSMesa and System GL. Strict traces are complete,
without dropped or truncated transitions, and 38 applicable HUD pixel checks
pass. Replay durations are 10,054/9,730/10,280/9,825/10,009/9,863/10,310/10,086 ms
in matrix order; transition counts are 432/410/487/468/313/283/353/329.
System warm wire again exercises the repaired observation: its completed peak
is 50.050568 ms and it interrupts at 50.057722 ms against 50 ms. The current
validator accepts this measured recovery. Ten inspected images include all
eight initial-ready views and the warm System wire pre-rotation/held pair.
Aircraft geometry and cold/warm placement agree; the held image retains the
aircraft and displays the interactive-detail bar. `generic/SUMMARY.md`,
`generic/contract-audit.json` and `visual-inspection.json` own these observations.

This changes observation and validation, not rendering policy or production
control ownership. The accepted formal baseline is unchanged. The preceding
queued-cancellation checkpoint retains its own 16 focused passes, six small
replays and original 7/8 Generic result. Full platform, client and package
qualification remains open, including coherent dependency staging. Lucy was
not rerun, and its prominent close-floor failure is not explained by this
validator repair. **S1--S6 remain incomplete.**

## September 7: queued cancellation and admission wakeup

`CXX-SOURCE-011` closes cancellation which waited for unrelated work to finish.
The predecessor leaves all three queued entries and their resources retained in
ten regression cases with every worker occupied or the source memory budget
fully reserved. Individual stream cancellation, whole-job cancellation, interest
handle destruction, cancellation before submission and streamless jobs reproduce
the defect. `/tmp/obol-queued-cancellation-20260907/tests-negative.log` and the
retained initial test source/executable and predecessor library own this evidence.

The existing coordinator now detaches cancelled entries with list splicing under
its queue mutex, then finishes items and retires their remaining counts outside
the lock. Streams receive a weak endpoint when submission commits; cancellation
before binding is detected at that boundary. Job cancellation also notifies after
closing streams, covering streamless jobs and shared-stream peers. Item cleanup
closes publication without recursively scanning the queue for each sibling.
Admission denial closes its batch and then sends one queue notification.

Two additional failures establish the cleanup ordering. An isolated library
without the new shutdown wait fails the actual process-shutdown observer
(`without-shutdown-wait/test.log`); its caller destructor still owns resources
after workers join. The coordinator now waits for in-flight queue retirement
before cache teardown. A candidate which notifies only after cleanup fails when
a destructor waits for newly admissible healthy peers
(`tests-frontier-negative.log`, 2.47 s). Queue removal now wakes workers before
cleanup and retirement completion wakes the shutdown barrier. These variants
preserve their own intermediate sources and binaries; they are not both variants
of the final source snapshot.

The final source test covers sixteen queued cases, adding concurrent individual
cancellation and shared streams under job cancellation and admission denial. It
checks source/context release, no cancelled callback execution, unchanged active
reservations, healthy sibling completion and admission wakeup. The shutdown
fixture retains a caller destructor across actual pool teardown and cancels a
surviving completed stream after coordinator destruction. Completed source
results keep the borrowed lifetime established by `CXX-SOURCE-010`.

Independent compilation exposed a missing `bu/str.h` include in the stream file,
previously supplied by a unity neighbor. Both source realization and occurrence
stream files now compile independently in the normal CMake build. The failed
standalone compile and final normal compile commands are retained. Configuration
again restaged the incompatible OSMesa library; after an empty host process audit,
the matching library was restored into this workspace with the explicit install
prefix recorded in the handoff. Coherent dependency staging remains S6 debt.

The final build passes. `tests-candidate.log` records sixteen focused
source/pool/stream/LoD/GED/Qt/contract checks passing in 13.20 s. Six final small
25-event rows pass functional, strict and applicable HUD checks in
`replays-final`: OSMesa merge-denied/reserve-denied/healthy take
1,123/1,150/1,099 ms with 80/80/75 transitions; System takes
1,040/1,040/1,063 ms with 61/61/57. Six inspected images preserve whole-root
coverage and full error bars on denial, and show four leaves on success.
The earlier `replays` directory is preliminary evidence from before independent
compilation and configuration; final qualification uses `replays-final`.

The final Generic Twin matrix has **seven PASS rows and one FAIL**: System warm
wire fails the responsive retained-cut assertion. Before rotation the render
time is 48.465424 ms, below the 52.5 ms threshold; the observed rotation peak is
also 48.465424 ms. The held frame has pixel error 1 and a ceiling of 0 instead
of -1. Presented lines change from 470,979 to 470,931, with no pending source
preparation throughout rotation. This repeats the earlier System warm-wire
failure; no causal link to queued cancellation has been established.
`generic/failed-row-audit.json` preserves the sample indices and predicate facts.

Follow-up inspection identifies omitted timing evidence: the first rotation wait
increments `presentation_interrupted_frames` from 3 to 4, with a 50.140745 ms
interrupted frame against the 50 ms hard deadline. Transition 277 records
`render-deadline`. The earlier retained System warm-wire failure also increments
that counter from 3 to 4 (50.140426 ms). The predicate considers only completed
render times and therefore misses this pressure. The original rows remain FAIL
pending a separately tested validator correction; `generic/pose-interruption-audit.json`
owns the cross-row comparison. Do not change rendering thresholds on this evidence.

All eight replay scripts succeed and all eight strict traces are complete,
without dropped or truncated records. All 41 applicable HUD pixel checks pass.
Seven camera contracts pass; the failed row stops before its camera contract.
Eight inspected initial-ready images show aircraft geometry with cold/warm
agreement. The failed row's pre-rotation and held images were also inspected:
the aircraft remains visible, and the held image shows the interactive-detail
bar. These image observations do not override the failed policy assertion.
`generic/SUMMARY.md`, `generic/contract-audit.json` and `visual-inspection.json`
own the row and image details. Lucy was not rerun.

The artifact retains scoped patches, source snapshots, writer/call-site deltas,
16 binaries and physical copies of all 78 libraries resolved by `ldd`. Library
copies are indexed by `retained-resolved-libraries.json`; this does not inventory
every dynamically loaded plugin. SHA-256 prefixes are qged `c2e92c0100bd`,
libBObol `2758c68c77a5`, libged `73b6278bfead`, libqtcad `aecad5fcc488`,
Obol `f281b59bdbd2` and OSMesa `097aac27db41`. Public function signatures are
unchanged, but the coordinator's private member changes its class layout;
full client/package qualification remains required. No formal relation, accepted
TLC baseline, scheduler or thread was added.

**S1--S6 remain incomplete.** Direct queue retirement closes an ownership edge;
its scan and payload destruction still scale with retained work. Real geometry
latency and byte bounds, necessary completed-result retention, snapshot-file
cleanup, remaining source/sparse writers, partial-merge allocation failure,
wire/Lucy quality, large/cold/shared geometry, editing and native/shared-stack
sanitizer qualification remain open. Private-display System GL is CPU llvmpipe
evidence and does not qualify a native GPU.

## September 7: terminal item resource retirement

`CXX-SOURCE-010` closes ownership retained by a terminal item while its sibling
continues. A six-case production regression holds a sibling callback active
while the first item completes, cancels or throws from a warm probe or cold
manifest callback. The old implementation retains every callback context and
unsuccessful source. Linux file-descriptor observation also finds retained
database handles after cold cancellation/failure and warm-probe failure.
`/tmp/obol-item-retirement-20260907/tests-negative.log` records the failures;
the initial test source, executable and matching old library are retained.

The item's `finish` method now owns terminal publication. Once callbacks return
or unwind, it releases unsuccessful streams and source/database/snapshot
resources, drops callback context and publishes item state with release ordering.
Source destruction precedes database close. The worker still releases its
reservation before aggregate terminal publication. Admission denial and queued
shutdown use the same cleanup outside the coordinator lock. The callback-owner
probe queries accounting during destruction, exercising that lock boundary.

Completed sources and their borrowed database handles remain valid until the
job handle is released, including after explicit cancellation. Completed
callback scratch retires immediately. `itemResult` acquires state before reading
the source pointer, and returns a source only for COMPLETE. GED clears its
transferred source pointer, removes the ownership flag and obtains the completed
result at adoption. Stream cancellation and the source-owned launch stamp still
control whether that result can be adopted. Public signatures/layout are
unchanged; exposing a source in non-COMPLETE states is deliberately removed.

The five source/pool/stream/GED checks pass in 4.20 s; eleven additional
LoD/Qt/contract checks pass in 10.58 s. The final three source/pool checks also
pass after adding a defensive null check to the fixture. Tests verify source
and context destruction while a sibling runs, live callback retention,
completed cold database readback after cancellation, and release after dropping
the result. Admission denial checks source/context retirement without a worker.
The process-shutdown observer verifies active and queued source/context release
while client handles are still retained, plus no queued callback execution.

Six small 25-event denial/success/policy/erase rows pass functional, strict and
applicable HUD checks. OSMesa merge-denied/reserve-denied/healthy take
1,127/1,125/1,104 ms, with 80/80/75 transitions; System takes
1,046/1,042/1,051 ms with 61/61/57. Six inspected images preserve whole-root
coverage and full error bars on denial, and four leaves on success.

Eight Generic Twin cold/warm shaded/wire rows and their camera contracts pass
on OSMesa and System GL. All strict traces are complete, without dropped or
truncated records; 37 applicable HUD pixel checks pass. Transition counts are
452/422/500/475/302/271/359/318 in matrix order. Eight inspected initial-ready
images show aircraft source geometry, with cold/warm agreement in both modes.
`generic/SUMMARY.md`, `generic/contract-audit.json` and `visual-inspection.json`
own the row and image details.

The artifact directory retains source snapshots, scoped patches, the complete
current writer delta, 16 binaries and hashes for 78 resolved libraries. SHA-256
prefixes are qged `4c34315fe051`, libBObol `5eb39b8025a9`, libged `02d38762c13e`,
libqtcad `4293b8ebaaf6`, Obol `f281b59bdbd2` and OSMesa `097aac27db41`.
No formal relation, accepted TLC baseline, scheduler or thread was added.

**S1--S6 remain incomplete.** These small fixtures prove release edges, not
large-geometry memory magnitude or cancellation latency. Ordinary cancellation
while an item waits for queue admission, necessary completed-result retention,
real snapshot-file cleanup and shared-stack sanitizers remain separate resource
qualification. The earlier synchronous stream-backlog measurements remain dated
evidence. Remaining source/sparse writers, allocation failure inside merge,
wire/Lucy failures, large/cold/shared geometry, native hosts, editing and coherent
dependency staging remain open. Lucy was not rerun; the passing Generic rows do
not close another adaptive timing history. Private-display System GL remains
CPU llvmpipe evidence, not native GPU qualification.

## September 7: cancelled stream ownership

`CXX-SOURCE-009` makes cancellation an ownership boundary in the existing
occurrence stream. The old implementation only set a flag: a cancelled item
retained its queued occurrence and 88-byte source import while a healthy sibling
continued, and all publication lanes accepted late output. The stream now
detaches queued geometry, priority coverage, staging and persistence journals
without allocating. Cleanup runs outside the stream mutex; subsequent writers
are rejected. Fixed progress facts remain observable. Claimed imports and
drained geometry retain their independent owners.

Successful job retirement preserves a completed stream transferred to a
consumer. Both the coordinator interest handle and GED adapter previously
requested cancellation unconditionally in their destructors. They now cancel
unfinished work while preserving successful consumer ownership. Explicit
cancellation still releases unconsumed data.

Evidence is at `/tmp/obol-stream-cancellation-20260907`. The two initial
regressions fail with the old code (`tests-negative.log`). The source fixture
now waits for only the healthy worker to remain before checking import expiry,
empty storage and rejection of a late valid import. The positive test drops a
completed producer handle, claims its surviving import, then verifies that
stream cancellation cannot revoke the claimed owner. The stream test covers
all late publication lanes, discarded journals and reentrant storage cleanup.

That cleanup probe exposed a shutdown consequence: cancelling under the
coordinator lock deadlocks when a storage owner queries accounting. The
intermediate candidate times out at 15.02 s (`shutdown-negative.log`). Shutdown
now closes admission, cancels active jobs outside the lock and retires queued
items through the same final-count helper as workers. Admission-constrained
cleanup also runs outside the queue lock. No additional scheduler or thread
was introduced.

Sixteen focused source/stream/LoD/GED/Qt/contract checks pass in 14.00 s. After
strengthening the source fixture's reservation-retirement barrier, its three
source/pool checks pass again. The actual pool-shutdown observer confirms
cancelled active/queued jobs, no queued callback execution, storage release and
eventual callback-context destruction. These checks retain the earlier partial
construction, submission exception and healthy-sibling contracts.

The final isolated `measure_cancellation` probe queues and seals 5k/50k/150k
metadata records, then cancels. Cleanup takes 2.24/24.05/86.33 ms, with zero
queued records remaining. At 150k, RSS is 417,935,360 bytes before cancellation
and 333,934,592 afterward; ownership release does not imply allocator pages
immediately return to the OS. These are measured costs, not a universal latency
bound or full distinct-geometry qualification. Synchronous cleanup scales with
backlog, and broader memory pressure/geometry/native limits remain open.

All six small 25-event delivery/policy/erase rows pass functional, strict and
applicable HUD checks. OSMesa merge-denied/reserve-denied/healthy take
1,123/1,121/1,098 ms with 80/80/75 transitions; System takes
1,056/1,027/1,062 ms with 61/61/57. Six inspected images retain whole-root
coverage and full terminal bars on denial, and four leaves on success.
The initial sandbox attempt could not connect Qt to Xvfb and produced no
report; `display-unavailable/` preserves it. The qualified private-display runs
use host display access.

All eight Generic Twin shaded/wire cold/warm rows and camera contracts pass.
Strict traces contain 472/451/491/459/305/277/351/339 transitions in matrix
order, without dropped or truncated records. Forty HUD pixel checks pass.
Eight initial-ready images were inspected and show source geometry on both
backends. `generic/SUMMARY.md` and `generic/contract-audit.json` own row details.

Current SHA-256 prefixes are qged `f152a9ffe690`, libBObol `6b2c24055b72`,
libged `fd8b000997f2`, libqtcad `4293b8ebaaf6`, Obol `f281b59bdbd2`, and
OSMesa `097aac27db41`. Retained binaries, loader/library hashes, source
snapshots and scoped patches identify the candidate. Formal cancellation
outcomes and the accepted TLC baseline are unchanged; allocation, container
ownership and lock mechanics remain executable qualification boundaries.

**S1--S6 remain incomplete.** Batch-held detached sources, database handles and
callback contexts have not acquired per-item retirement. Remaining source
writers/journals, allocation failure inside merge, earlier wire/Lucy failures,
large/cold/shared geometry, native hosts, editing, sanitizers and coherent
dependency staging retain their independent requirements. Lucy was not rerun;
this Generic pass does not retire another adaptive timing history.

## September 7: framebuffer feature evidence

`CXX-PRESENTATION-006` authenticates the HUD content in a captured image using
the feature-store revision traversed by the displayed frame. Exact CAD work
and an empty render-request latch are insufficient: an interrupted traversal
can execute CAD, retire its latch and retain older complete pixels while
recovery is pending. The source-identity checkpoint's System wire cold images
at `history-return-stable` and `zoom-in-motion` are byte-identical, both at
presentation/completion serial 42. The intervening traversal is interrupted
at 50.01291 ms. `saved-frame-audit.json` retains this causal sequence.

Evidence is at `/tmp/obol-frame-feature-evidence-20260907`. Rendering captures
the existing feature revision before completion feedback can publish another
overlay. Both Qt hosts retain that revision with completed pixels and expose
a passive observation; provisional images have no certificate. qged reports
current and presented revisions as strings. The HUD checker requires them to
match before comparing nonterminal fill geometry with pixels. Missing or
malformed identity fails, and a terminal error with no pending render must
have reached a matching framebuffer. The 100-pixel floor, fill color and
terminal bar/label requirements remain unchanged. No scheduling or LoD policy
decision changed.

The old checker fails the new interrupted-frame synthetic case. The fixed
checker passes all 37 cases, including rejection of a hidden fill after
presentation resumes and of an unpresented terminal error. The real Qt
deadline fixtures publish a new HUD label, abort a traversal, retire the
request, and require unchanged pixels carrying the old revision. Recovery
must paint new pixels carrying the new revision. Passive capture preserves
pending work; resized software pixels cannot certify the old viewport.
Both software and System GL paths pass under private Xvfb at DPR 1.5.

All nine focused Qt/prototype/contract tests pass in 1.90 s. Rebuilding the
prototype exposed a stale direct-BoT stream assertion: the earlier overview
repair legitimately emits a priority coverage record beside one authoritative
leaf. The test now validates that optional coverage record separately while
preserving the leaf, manifest and source-profile requirements. Its negative
diagnostic and fixed pass are retained; this required no source-producer change.
The controller's optional render output and Qt's new virtual observation
change the C++ ABI; affected repository callers were rebuilt. Full shared-client
and package qualification remains S6 work.

All six small 25-event replays pass functional, strict and applicable HUD
checks. OSMesa merge-denied/reserve-denied/healthy take 1,111/1,122/1,113 ms
with 80/80/75 transitions; System takes 1,048/1,050/1,072 ms with 61/61/57.
All six inspected `on.png` images carry matching current/presented revisions:
denied delivery retains whole-root coverage and full red terminal bars;
healthy delivery shows four leaves. Policy changes and erase pass.

All eight Generic Twin shaded/wire cold/warm rows and their applicable camera
contracts pass. `generic/contract-audit.json` records each result. Strict
transitions in matrix order are 449/401/482/467/305/275/377/333. The HUD audit
checks 38 images and excludes one older frame in System warm wire at
`zoom-in-motion`: current feature revision 119, presented revision 112,
presentation serial 36, one interrupted traversal. Its pixels equal the
preceding stable image. Against this same report the old checker exits 1 with
zero fill pixels; the fixed checker exits 0 and checks five current fills.
Eight representative HUD images and that retained-image pair were inspected.
This paired evidence preserves the historical failed report and establishes
why its image assertion was invalid.

Current SHA-256 prefixes are qged `377ea00ea466`, libBObol `4b3ad0abe1e9`,
libged `97daa97a5211`, libqtcad `f900369188e5`, Obol `f281b59bdbd2`, and
OSMesa `097aac27db41`. Eleven binaries and 78 resolved libraries have retained
hashes; scoped implementation/documentation patches and source snapshots
identify the candidate. No formal relation or accepted TLC baseline changed.

**S1--S6 remain incomplete.** Current Generic passes do not retire the earlier
System warm-wire retained-cut failure under a different timing history. Lucy
was not rerun; its software close-floor failure remains open. Remaining source
writers, journal/cancelled-payload retention, partial-merge allocation, large
and cold workloads, native hosts, coherent dependency staging, editing and
shared-stack sanitizers retain their qualification requirements.

## September 7: source realization identity

`CXX-SOURCE-008` rejects old output after same-owner population replacement or
source reconfiguration, using one immutable stamp captured and authenticated by
`SoBRLDatabaseSource`. The GED adapter's separate field list and predicate are
removed. Stream appends and camera/view policy changes remain independent of
source content admission. Database rebinding uses a non-reused epoch, including
A to B to A, without retaining an old database solely for address comparison.
The public detached-adoption API now requires the prelaunch stamp; all repository
call sites are updated.

Evidence is at `/tmp/obol-source-identity-20260907`. Before the repair, the real
GED pump adds an old overview to a replacement population, growing it from four
entries to five (`test-negative.log`). The negative executable, test source and
preceding production source are retained. The fixed test covers nine replacement
cases while preserving numeric revisions and modes: owner, path, population,
database, database round trip, representation key and each tessellation tolerance.
It verifies retained geometry, placement, bounds, population and realization
status until stale preparation retires without error. A separate update-action
regression rejects a superseded population through direct unstreamed adoption;
existing positive tests preserve valid appends, selection and camera independence.

The first focused pair passes in 8.37 s. After final cleanup and rebuilding the
affected test targets, all 26 focused source/Qt/contract checks pass in 13.77 s,
including seven conformance checks. No formal relation or accepted TLC baseline
changed; this repair strengthens implementation of existing result authentication.

All six small 25-event GUI replays pass functional, strict-trace and applicable
HUD checks. OSMesa merge-denied/reserve-denied/healthy take
1,125/1,135/1,086 ms with 80/80/75 strict transitions; System takes
1,052/1,054/1,069 ms with 61/61/57 transitions. Six images were inspected: failed
delivery retains whole-root coverage and a full red terminal bar, healthy
delivery shows four leaves, and policy changes and erase remain correct.

The Generic Twin matrix in `generic/` has six passing rows: four OSMesa
shaded/wire cold/warm and two System shaded cold/warm, including their camera
contracts. System wire cold fails HUD validation; its dependent warm row is
skipped. All seven completed replays report success and pass strict traces
(432/405/493/478/313/277/363 transitions in matrix order). Seven initial-ready
images and the failed zoom image were inspected. Actual source geometry is
present on both backends.

The failed `zoom-in-motion.png` has **zero HUD fill pixels**. Its report records
an exact CAD frame, no pending render, 466,021 presented lines, all 2,249 source
leaves available, and no source preparation pending. Retained features describe
a yellow seven-pixel-wide fill from y=115.8928 to y=168.3248 and the label
`Interactive detail`, but neither appears in the image. This is nonterminal
(phase 2, no terminal error), distinct from the terminal-error publication
regression. `generic/wire-cold-hud-failure.json` preserves the facts. Its causal
relationship to this source repair is not established. Keep the HUD predicate
and 100-pixel floor unchanged; this is a failed release row.

Current SHA-256 prefixes are qged `afe41c249564`, libBObol `4882deed9a37`, libged
`84be374351b2`, libqtcad `e1664b7e5f94`, Obol `f281b59bdbd2`, and OSMesa
`097aac27db41`. `runtime-hashes.json`, `resolved-library-hashes.json`, retained
binaries and scoped source/documentation diffs establish provenance. The OSMesa
staging hazard remains open. Lucy was not rerun; the previous software close-floor
failure and System warm-wire retained-cut failure remain open. **S1--S6 remain
incomplete.** Remaining source writers, journal retention, cancelled storage,
partial-merge allocation, large/cold workloads, native hosts and sanitizers
retain their independent qualification requirements.

## September 7: terminal-error HUD publication

`CXX-PRESENTATION-005` distinguishes active and terminal errors in the private
Qt HUD publication grouping. The controller already reports the terminal
outcome; this repair ensures the final retained bar is published when another
provider retires, even if no periodic sample or geometry frame remains.

Evidence is at `/tmp/obol-terminal-error-hud-20260907`. The finalized C++
regression fails with the old grouping (`final-test-negative.log`) because no
terminal sample is published. It uses a real GED/controller scene, a controlled
pending provider and the production publication helper with periodic sampling
disabled. Real rendering retires the remaining frame obligations. The fixed
test requires a full red bar and error label, and verifies unchanged feature
revisions and no new render/progressive work after the final HUD frame.

The HUD oracle now rejects an absent or unfinished terminal-error fill instead
of skipping it. The retained System GL report fails the new assertion in
`hud-oracle-saved-negative.log`; 29 synthetic cases cover normal eligibility,
terminal appearance, missing/hidden fill, pending frames and disabled LoD.
Eight focused Qt/contract checks pass. These include the existing passive
capture, window-host and faceplate checks as well as the new transition.

All six small 25-event replays pass functional, strict-trace and applicable
HUD checks. OSMesa merge-denied/reserve-denied/healthy take
1,135/1,115/1,088 ms with 80/80/75 strict transitions; System GL takes
1,047/1,053/1,072 ms with 61/61/57 transitions. Both denied backends now
retain full red bars with the correct error label and no pending render at
`on.png`. Six current images were inspected. Policy off/on/off, whole-root
coverage, healthy delivery and erase remain correct.

Current hashes: qged `96a30741...`, libqtcad `ccf1b1f1...`, libBObol
`87f278f4...`, libged `384d1443...`, Obol `f281b59b...`, OSMesa
`097aac27...`. Reconfiguration reproduced the dependency restaging failure:
OSMesa `75dfd5f3...` caused a test SIGSEGV; restoring the matching workspace
library removed that failure. Coherent dependency staging remains S6 debt.

Generic Twin and Lucy were not rerun for this terminal-error-only change.
Their previous warm-wire retained-cut and software close-floor failures remain
open predecessor evidence in the following entry. Native GPU, resource,
remaining ownership/extraction and shared-client qualification remain open.
No formal relation or accepted TLC baseline changed. S1--S6 remain incomplete.

## September 7: passive checkpoints and HUD frame correspondence

`CXX-PRESENTATION-004` closes an observation defect. In the retained routing
matrix, System GL's `subpath-erased` and `subpath-redraw-return` images are
pixel-identical at presentation serial 67. The latter sample describes newly
published HUD geometry with a pending render. Asserting that unpresented fill
against the preceding image produced the reported zero-pixel failure.
`/tmp/obol-hud-frame-evidence-20260906/saved-hud-frame-mismatch.json` records
the independent image comparison and both snapshots.

Software checkpoints also used image export, executing a renderer callback
without completing a presentation. `negative-passive-capture.log` observes
one traversal; its pixels happen to remain unchanged in this fixture.
`QgSW::get_presented_frame_image` now samples painted pixels, and qged shares
one passive capture helper across backends. The software canvas keeps a
reference to the displayed image separately from its exact abort fallback;
ordinary completed presentations share the same storage. The getter rejects
viewport size or pixel-ratio mismatch. The first GUI run exposed a startup
barrier accepting a pre-resize paint; the barrier now requires an available
image for the requested viewport. That negative candidate and its binaries
remain in `passive-before-viewport-barrier/` and
`/tmp/qged-obol-passive-checkpoint-generic-20260907`.

The HUD checker now inspects every sufficiently large non-idle fill with an
exact CAD frame and no outstanding render request, for both faces and wire
lines. It retains the 100-pixel floor and color tolerance. Invalid presentation
state, malformed geometry, missing images and any hidden eligible fill fail;
pending/incomplete frames cannot authorize an image assertion. The synthetic
negative reproduces the stale-frame selection, and 21 cases include rejection
of a later hidden fill after an earlier visible one. This is a corrected
observation prerequisite and broader pixel coverage, not a relaxed floor.

Seven focused Qt and contract checks pass in 1.45 s. Both passive-capture
branches also pass under private Xvfb at display scale 1.5: no rendering
callback, unchanged frame/request serials, pending work preserved, and software
readback equal to QPainter's actual orientation and logical pixel scale.
All six small delivery/policy/erase replays pass functional and strict checks;
six current `on.png` images were inspected. OSMesa denied/reserve-denied/healthy
have 76/76/75 transitions; System GL has 58/58/57.

The final Generic Twin matrix at
`/tmp/qged-obol-passive-checkpoint-final-generic-20260907` passes seven of eight
rows and their applicable camera contracts. All eight strict traces pass.
The final HUD checker passes all 39 eligible images, with no skipped row;
eight representative images were inspected. System GL warm wire still fails
the retained-cut policy check: pre-rotation cost is 47.67 ms, rotation peaks at
50.74 ms below the 52.5 ms threshold, but the held ceiling is 0 instead of -1.
The failed warm row and earlier cold failures remain open.

Current warm OSMesa Lucy at
`/tmp/qged-obol-passive-checkpoint-lucy-osmesa-20260907` completes all 217 events
in 106,975 ms with 2,351 strict transitions, but fails its close-view quality
floor at events 106/108/110: cut 28, 2,619,533 presented faces, normalized error
3.5813 and one prominent-floor violation. All four normals checkpoints use
cut 24 and 2,101,208 faces; the close image, four normals images and final image
were inspected. This is a failed qualification row, not a current Lucy pass.
The causal relationship between the sampling change and the adaptive failure
has not been established. System GL Lucy remains predecessor evidence.

The small System GL error replay exposed a separate publication
defect: its estimate is 1 with no pending render, but the retained fill is
absent; OSMesa has a full red fill. Both display the correct error label.
`terminal-error-hud-difference.json` and `error-publication-gdb/` retain the
snapshots and published error stages. Its active/terminal grouping omission is
addressed by the later [terminal-error checkpoint](#september-7-terminal-error-hud-publication).

Checkpoint hashes are qged `96a30741...`, libqtcad `347d68e1...`, libBObol
`87f278f4...`, libged `d7437100...`, Obol `f281b59b...` and OSMesa
`097aac27...`. No formal relation, strict transition checker or accepted TLC
baseline changed. Native GPU/platform, resource, remaining writer/extraction,
editing/shared-client and coherent dependency qualification remain open.
S1--S6 are incomplete.

## September 6: source routing and scoped cancellation

`CXX-SOURCE-007` authenticates the captured source route, path, source/input
revisions and modes before GED consumes streamed metadata or geometry and
before terminal adoption. Rejected delivery cancels its existing stream;
cancelled queues no longer strand provider retirement. Adoption evaluates each
item independently. Workers observe individual stream cancellation at stage
and completion boundaries, allowing healthy siblings to finish. This removes
the redundant delivery-denial flag and unused view-revision snapshot; it adds
no controller, persistent state or public API. Job-wide failure rules remain
unchanged.

`/tmp/obol-source-routing-20260906` retains negative/fixed tests, source diffs,
executables and hashes. Against the predecessor libraries,
`negative-routing.log` shows old queued output growing an in-place replacement
from one occurrence to two. `negative-owner-and-cancellation.log` shows a new
owner with the same key/path/stamps growing from four to five, and a separately
cancelled worker item reporting COMPLETE. Its healthy peer was not cancelled
in that negative run. The predecessor shared libraries remain archived in
`/tmp/obol-overview-extent-20260906/final`.

The fixed GED regression preserves replacement bounds, geometry, placement,
population and REALIZED state through every pump, and retires preparation
without error. The coordinator test reports CANCELLED for the interrupted item,
keeps a held sibling active, then completes that sibling and releases its worker
reservation. It tests two items in the real coordinator; the GED replacement
cases use one-item jobs and synchronous replacement geometry.
`focused-final.log` passes 23 checks in 15.03 s, including GED, Qt, source/pool,
streaming, rendering and controller/service contracts. The updated regression
catalog and all seven conformance checks pass in 0.56 s.

All six small 25-event GUI rows pass functional assertions and strict tracing:
OSMesa merge-denied/reserve-denied/healthy have 71/71/70 transitions; System GL
has 58/58/57. Six current `on.png` images were inspected. All 24 checkpoint scene
regions are pixel-identical to the inspected overview predecessor; the comparison
excludes the lower timing text, as recorded in `scene-pixel-comparison.json`.
This preserves the bounded coverage/policy/erase evidence; it is not inspection
of all 24 full current images or a repair of the System GL HUD.

The full Generic Twin matrix at
`/tmp/qged-obol-source-routing-generic-20260906` passes four OSMesa rows and their
camera contracts. System GL cold shaded again fails zero HUD fill at
`subpath-redraw-return`. System GL cold wire fails the retained-cut policy check:
pre-rotation cost is 49.66 ms and rotation peaks at 50.65 ms, both below the
52.5 ms threshold, yet the held checkpoint sets the progressive ceiling to 0
instead of -1. `retained-cut-diagnostic.json` preserves the exact values. Both
System warm rows are skipped. All six executed traces pass the strict checker;
seven representative images covering those rows and the held rotation were
inspected. The earlier wire pass does not close this recurring failure, and
this run does not establish that the routing repair caused it. Lucy was not
rerun; its older binaries remain predecessor evidence.

Current hashes are qged `d72c0061...`, libBObol `87f278f4...`, libged
`d7437100...`, libqtcad `6a1a9193...`, Obol `f281b59b...` and OSMesa
`097aac27...`. No formal relation or accepted baseline changed. The GED scenario
maps to result authentication as a regression, not a stepwise trace; the
one-provider lifecycle model does not prove multi-item worker cancellation.
Cancelled payload storage can remain until the whole job retires, even after
its worker reservation is released. Same-lifetime source reconfiguration and
population identity, storage retention beside a running sibling, large-geometry
cancellation latency, remaining writers, merge exception safety and coherent
dependency staging stay open. The next bounded graphical investigation is HUD
pixel/report frame correspondence. S1--S6 remain incomplete.

## September 6: retained overview extent and delivery

`CXX-PRESENTATION-003` preserves the initial view-envelope geometry/transform
pair, publishes generic exact coverage through the existing priority lane, and
advances placement when replacement geometry changes coordinate systems. GED
records the expected population without allocating its storage and commits the
priority overview before attempting leaf reservation or merge. This adds one
count-only source API and reuses the internal coverage builder; it adds no
controller, persistent control state or stream lane.

`/tmp/obol-overview-extent-20260906` retains the sources, implementation diff,
executables, hashes and loader evidence. `negative-transform.log` measures a
unit box where the intended preview spans approximately
(-1.597,-1.842,-1.573) to (1.597,1.842,1.573).
`fixed-transform.log` passes that boundary. `negative-coverage.log` and
`negative-generic-publication.log` expose missing whole-root delivery and the
generic producer's absent priority overview. With that producer repaired,
`negative-retained-placement.log` measures the old transform applied to new
source-local vertices: approximately (-4.790,-5.526,-4.720) to
(295.361,82.896,73.942), instead of (-1,-1,-1) to (93,23,24).
Advancing placement closes this final counterexample.

`fixed-placement.log` passes GED, Qt, source realization, compact streaming and
the private-symbol contract (5/5, 4.42 s). `focused-final.log` passes 19 related
checks, including real rendering, controller/service and the Qt recheck
(14.96 s). All seven conformance checks pass after the regression catalog update
(0.50 s). The GED fixture measures actual wire vertices under both the initial
and updated retained instance transforms, captures a successful overview commit
before denying the next leaf merge, and retains the expected population, error,
policy off/on/off rendering, erase and fresh-draw recovery. The count-only test
records `SIZE_MAX` without creating an instance index. That is allocation
independence evidence, not a physical OOM experiment.

All six 25-event small replays in `replays/` pass functional assertions and the
unchanged strict trace checker. OSMesa merge-denied/reserve-denied/healthy have
71/71/70 transitions; System GL has 58/58/57. All 24 images were inspected.
Denied delivery retains the whole-root box, healthy delivery shows four leaves,
and erase returns to the ready marker. Reported planning bounds equal the exact
root bounds (-4,-4,-1) to (4,4,1); the C++ test separately proves actual vertices.
The original camera sequence is unchanged, including the initially off-screen
root. The private reserve fault precedes instance reservation; the merge fault
now precedes leaf merge after overview commit. Neither tests allocation failure
inside mutation. These completed-producer faults are only for streams small
enough to complete without consumer drainage.

The full Generic Twin matrix at
`/tmp/qged-obol-overview-extent-generic-20260906` passes four OSMesa and two
System GL wire rows, including their camera contracts. System GL cold shaded
still fails zero HUD fill at `subpath-redraw-return`; its warm row is skipped.
All seven executed rows pass strict tracing. A representative image from every
row was inspected; this does not close the older adaptive retained-cut failure.
Lucy was not rerun, and its predecessor evidence does not qualify these binaries.

Final hashes are qged `3e749079...`, libBObol `2913f99f...`, libged
`31b39b00...`, libqtcad `6a1a9193...`, Obol `f281b59b...` and OSMesa
`097aac27...`. CMake initially restaged an updated external bundle and replaced
OSMesa with `75dfd5f3...`; Qt then crashed in `_mesa_GetError`. The offscreen gdb
stack and `negative-osmesa-only.log` reproduce the failure, including selection
of that library with unchanged final C++ binaries. Reinstalling the matching
local OSMesa artifact restores the qualified hash and the Qt pass. This repairs
the checkout, not dependency reproducibility: S6 still requires a coherent,
pinned dependency set that survives configuration.

No formal relation or accepted baseline changed; the existing complete TLC
result remains predecessor model evidence. Whole-root coverage before exact
bounds exist, multiple-root failure, strong exception safety within merge,
physical memory pressure and full platform/resource qualification remain open.
System GL under private Xvfb is CPU llvmpipe, not native GPU qualification.
S1--S6 remain incomplete.

## September 6: policy and erase publication boundaries

`CXX-TRACE-001` uses the existing source-input notification and diagnostic
scope. `clearViewLodState` now records cancellation and binding retirement in
one external-input transition. Notifications record committed scene changes
with automatic LoD off while leaving automatic planning disabled. GED notifies
attached views after shared policy/source publication and an independent view
after its local source publication, before clearing its old view state. The
existing post-frontier notification remains. No event kind, control owner,
public API, or checker predicate was added.

`/tmp/obol-policy-trace-20260905` retains the predecessor, gdb stacks,
`implementation.patch`, tested sources/executables, full hashes and loader
records. `negative.log` reproduces the GED policy/erase failure. The corrected
controller fixture fails specifically at the clear boundary against predecessor
libBObol `c11db473...` in `negative-final.log`; it also requires a no-op
notification to leave a later unannounced write unnamed without scheduling
automatic work. The early fixture left a simulated frame claimed; its separate
setup failure and correction are retained.

The controller/GED/Qt checks pass in 13.63 s. Nineteen additional focused checks,
including actual rendering, pass in 5.63 s; seven conformance checks pass in
0.32 s. qged and affected consumers were rebuilt. Current hashes are qged
`4659c4ec...`, libBObol `a570f237...`, libged `8d785e20...` and libqtcad
`6a1a9193...`; Obol `f281b59b...` and OSMesa `097aac27...` are unchanged.
No formal relation or accepted baseline changed; the full TLC suite was not
repeated for this publication-boundary repair.

All four 25-event small replays pass their functional assertions and the
unchanged strict trace checker: OSMesa denied/healthy record 71/70 transitions,
System GL denied/healthy 60/57. All 16 images were inspected. Denied delivery
retains a visible error preview; healthy delivery shows four authoritative
boxes; erase returns to the ready marker. **Whole-root preview extent remains
unqualified:** the retained overview is still a unit box. These runs do not
prove physical OOM handling, latency limits, multiple-view refinement or native
GPU behavior.

The full Generic Twin matrix at
`/tmp/qged-obol-policy-trace-generic-20260906` passes all four OSMesa rows and
both System GL wire rows, including their camera contracts. System GL cold
shaded still fails the HUD fill check at `subpath-redraw-return` (zero fill
pixels); its warm row is skipped. All seven executed rows pass strict tracing,
and representative images from each were inspected. The current wire passes
do not causally retire the preceding timing-dependent retained-cut failure.
Lucy was not rerun at these library hashes; preserve its preceding evidence
and known adaptive failures. S1--S6 remain incomplete.

## September 5: exact zero-work assembly aggregation

`CXX-PRESENTATION-002` fixes a completed-frame aggregation defect. An exact
zero contribution from a fully culled assembly invalidated the total primitive
and render-cost observations. With a denied source outside the current view,
the visible marker rendered 48 lines, but the controller repeatedly requested
`lod-exact-payload-replay`. The preceding System GL and OSMesa runs both timed
out after 10 seconds at event 7; they are retained under
`/tmp/obol-ged-stream-denial-20260905/gui-delivery/{system,osmesa}-denied-v2`.
Framing the entire scene avoided that failure, which is why the earlier
headless off/on/off pass did not establish this boundary.

`refreshCadPresentationFrameStatus` now sums exact zero contributions. Empty
assemblies are already excluded; execution, requested-control and exactness
checks remain intact. Policy continues to decide whether a measured aggregate
is useful capacity evidence. No new state or model relation is introduced.

`/tmp/obol-zero-work-aggregation-20260905` retains predecessor sources/library,
the final renderer executable, loader records and `implementation.patch`.
`negative-final.log` fails the final regression against libBObol `28ac9da7...`;
`fixed-quality.log` passes in 0.21 s. The real quality renderer checks mixed
visible/culled assemblies, exact zero for an entirely culled frame and rejection
of an unexecuted frame. The initial direct-wire fixture measured submitted lines
before pixel clipping, so its setup failures are retained separately.

All 19 focused CTests pass in 26.54 s, GED/Qt integration in 4.71 s and seven
conformance CTests in 0.38 s. qged and affected consumers were rebuilt. Current
hashes are qged `4659c4ec...`, libBObol `c11db473...`, libged `eca61503...`,
libqtcad `6a1a9193...`, Obol `f281b59b...` and OSMesa `097aac27...`;
`hashes.json` and `ldd.txt` record full identities. The renderer test is
`dd20942a...` and static-linked GED test `75118705...`. The build remains Debug
without sanitizers. Catalog lint passes 45 pairs; the unchanged CAD-frame
composition matches the accepted baseline. Concrete work magnitude is below
that abstraction; no further formal baseline is accepted.

All four small System GL/OSMesa healthy/denied replays complete 25 events in
about 1.0--1.2 seconds under `gui`. Their functional checks require retained
failure through off/on/off, current bounds, no source-preparation debt and ready
recovery after erase. Healthy images show four boxes and the marker. **The full
rows fail:** the unchanged strict checker rejects unnamed transitions at policy
changes and erase, also present in the predecessor's all-visible replay.
`predecessor-unnamed-transitions.json` identifies CAD/resident-demand, source
population and outcome changes outside named boundaries. Denied preview images
also show a small unit box; compact planning bounds are `[0,1]` while certified
source bounds are `[-4,-4,-1]` to `[4,4,1]`. Retained visibility is demonstrated;
whole-root preview coverage still requires investigation and qualification.

The affected Generic/Lucy rows are terminal, with unchanged validators:

| Artifact under `/tmp` | Current result |
|---|---|
| `qged-obol-zero-work-generic-20260905` | Four OSMesa rows pass. System cold shaded fails HUD fill (zero pixels); System cold wire fails responsive retained-cut policy, and both System warm rows are skipped. The wire rotation peaks at 47.17 ms against 52.5 ms, but applies ceiling 0. |
| `qged-obol-zero-work-osmesa-20260905` | Full warm row passes: 217 events, 128,047 ms, 2,443 strict transitions. Smooth zoom starts at cut 24, progresses through effective cuts 24/24/26 and settles close at effective cut 29 with 3,680,298 presented faces and no prominent-floor violations. |
| `qged-obol-zero-work-glsl-20260905` | Full opt-in warm row passes: 217 events, 84,039 ms, 980 strict transitions. Smooth zoom starts at cut 24 and settles close at cut 29 with 3,445,998 presented faces and no prominent-floor violations. This is private-Xvfb CPU llvmpipe, not native hardware-GPU qualification. |

All four normals checkpoints in each Lucy row use cut 24 and 2,101,208 presented
faces; all eight images were inspected, along with close/zoom-out and
representative Generic checkpoints. OSMesa reports 6,303,624 normals in each;
programmable flat/authored report zero and smooth reports 4,356,455. The sibling
`generic-summary.json` and `lucy-{osmesa,glsl}-summary.json` retain the relevant
measurements. These cut-24 passes do not close the older cut-25 software
failure, nor establish a cause for the prior programmable close-floor failure.
Preserve those adaptive reproductions; no graphical timing improvement is
causally attributed to the aggregation repair.

Native hosts, physical OOM, multiple failed roots and shared-stack sanitizers
remain open, as do S1--S6. All build, test, formal and GUI handles are terminal.

## September 5: failed source delivery and display transitions

`CXX-SOURCE-006` prevents a completed worker from bypassing a memory-denied
owner-thread merge through detached adoption. The source retains preview
geometry and expected inventory, publishes its existing FAILED outcome, and
stops counting absent, undeliverable leaves as pending work. Current-source error counts replace a
historical synchronous-action count. View-only invalidation preserves a failed
external producer's outcome, and policy-off convergence honors explicit source
work and errors. No persistent C++ control state is added.

The expanded rendering regression found `CXX-PRESENTATION-001`: switching from
LoD back to an ordinary source representation left an unused cached assembly
in the required frame-execution set. The representation boundary now retires
that binding through its existing setter. Exact-frame checks remain unchanged.

`/tmp/obol-ged-stream-denial-20260905` retains originals, the instrumented negative
executable, final sources/executable, `implementation.patch` and logs.
`negative.log` demonstrates erroneous adoption after actual worker completion;
`diagnostic.log` demonstrates why skipping adoption alone leaves permanent
inventory work. `display-diagnostic.log` records failure erased by policy
invalidation. `display-frame.log` and `display-assembly.log` identify the unused
ordinary source's LoD assembly holding exact presentation open. With both
boundaries repaired, `display-retirement.log` passes in 4.09 s: the real GED
provider preserves its preview, renders through off/on/off, reaches terminal
error, and a fresh draw produces four authoritative occurrences.

All 19 focused CTests pass in 16.03 s; GED and Qt drawing integration pass in
4.17 s total; all seven `bobol_conformance` CTests pass in 0.37 s. qged and
affected consumers were rebuilt. `hashes.json` and `ldd.txt` identify qged
`4659c4ec...`, libBObol `28ac9da7...`, libged `eca61503...`, libqtcad
`6a1a9193...`, Obol `f281b59b...` and OSMesa `097aac27...`. The static-linked GED
test is `e6046df0...`. External renderers and public C++ layouts are unchanged;
the build is Debug without sanitizers.

The lifecycle composition now models provider failure and error preservation
across policy changes, with success/failure coverage required by the catalog.
The initial focused relation passes; a mutant hiding policy-off error fails
`ErrorHasWitness` in three states. The final wrappers expose separate action
coverage without changing that relation. Catalog lint and the complete suite
pass all 45 pairs at `tla-full`; all 61 required actions have nonzero coverage.
`CompleteProvider` and `FailProvider` each execute 58,779 times. The review in
`baseline-audit.json` verifies the unchanged tool pin and lifecycle-only
count/depth change; the complete result is explicitly accepted as the current
baseline. The canonical check took 29 min 7 s; model-process time totals
34 min 14 s. The new failure slice covers one provider obligation, with the
independently live-provider projection checked by the C++ coordinator test.

The headless test checks real rendering and retained geometry, not pixel
fidelity or performance. Physical OOM, real geometry under pressure, multiple
failed roots, native hosts and shared-stack sanitizers remain unqualified.
The full Generic/Lucy matrix has not been rerun at this library identity; its
earlier graphical failures remain open. These C++ repairs do not close S1--S6.

## September 5: authoritative source adoption

`CXX-SOURCE-005` removes the shortcut that treated matching occurrence metadata
as proof that a detached payload was already installed. Only the existing
explicit drained-stream certificate now permits reuse; other callers install
the detached registry. GED's normal authoritative streaming path retains its
existing fast adoption. No new state, registry mirror or rendering policy is
introduced.

`/tmp/obol-source-adoption-20260905` retains original sources/library, the final
test executable with its negative loader record, patches and validation.
`negative-final.log` fails the new regression in 7.00 s against predecessor
libBObol `b97409d8...`. The final test passes with the fix and requires the new
immutable geometry in the retained assembly, retained selection, stable source
routing, advanced population identity and invalidated inventory deltas. Existing
streamed overview-retirement and staged-provider lifetime tests also pass.

All 19 focused CTests pass in 13.99 s; GED and Qt drawing integration pass in
4.18 s total. qged and affected consumers were rebuilt. `hashes.json` and
`ldd.txt` identify qged `4659c4ec...`, libBObol `4798589c...`, libged
`6395306b...`, libqtcad `6a1a9193...`, Obol `f281b59b...` and OSMesa
`097aac27...`. These are Debug tests without sanitizers. Public layouts and
external renderer libraries are unchanged.

Source-evidence and CAD-frame composition checks preserve their baseline
counts/depths, and catalog lint passes all 45 pairs. Concrete payload comparison
is below those abstractions; the negative/fixed executable supplies its proof.
No model, conformance mapping or baseline was changed or accepted. The prior
source-demand checkpoint retains the latest complete formal-suite run.

The full Generic/Lucy matrix was not rerun. Its preceding zoom, prominent-floor,
normals, HUD and frame-delivery failures remain open, and historical passes do
not qualify these binaries. This bounded adoption repair does not qualify all
source mutation/pressure paths or close S1--S6.

## September 5: source worker-pool construction and shutdown

`CXX-SOURCE-003/004` close two pool-lifetime defects and the associated GED
startup-error escape.  Construction now owns its candidate through RAII and
joins partially started workers before propagating an error.  GED contains
startup errors while it still owns the requests.  Shutdown uses the same stop
path, reaching active jobs through their fixed worker slots as well as queued
jobs; it cancels their streams before joining.

`/tmp/obol-source-pool-lifetime-20260905` retains original sources/library,
instrumented negative executables, the final patch, writer inventory and logs.
`negative.log` shows both new source tests failing against the instrumented
predecessor `719a1549...` (original library `43e78f04...`): partial construction
leaves two threads against a baseline of one, and shutdown fails its callback
cancellation/retirement check after 2.06 s.  The separate `ged-negative.log`
fails in 0.07 s because startup failure escapes the draw transaction.

The final construction test checks three failed starts independently return to
the baseline thread count, then creates exactly the requested two-worker pool.
Its independent `/proc/self/task` oracle makes this CTest Linux-specific.
The shutdown test keeps client handles alive across actual coordinator
destruction.  One cooperative callback holds the full source allowance while
a separate two-item job remains queued.  The post-destruction observer requires
cancellation, no queued callback execution and released callback contexts.
The strengthened GED fixture requires retained nonempty proxy bounds, no
pending source producer, and successful erase/redraw into its four real
occurrences after removing the injected fault.

Both new source tests have headless and sanitizer labels; the latter names a
required future sanitizer run, not a sanitizer result.  All 19 focused CTests
pass in 13.56 s (`tests-focused.log`).  GED and Qt drawing integration tests
pass in 4.45 s total (`tests-integration.log`).  These are Debug tests with
sanitizers off.  The original source submission, callback lifetime, admission,
aging and cache-shutdown regressions remain passing.

qged and affected consumers were rebuilt.  `hashes.json`, `ldd.txt` and
`ged-ldd.txt` retain full identity and loader evidence: qged `ccf76635...`,
libBObol `b97409d8...`, libged `bd917708...`, libqtcad `6a1a9193...`, Obol
`f281b59b...` and OSMesa `097aac27...`.  The GED regression links the production
static library objects; its negative and final executable hashes are retained.
External Obol/OSMesa and public C++ layouts are unchanged.

The control-lifecycle composition and deferred-autoview models pass with
unchanged tool identity, state counts and depths (`focused-model-audit.json`).
Catalog lint passes all 45 pairs.  Thread creation/join mechanics remain
outside those control abstractions; the new executable regressions supply the
pool evidence.  No model, conformance mapping or baseline was changed or
accepted.  The source-demand checkpoint retains the latest complete formal run.

Worker slots add one current-job reference per configured worker, bounded by
the existing worker count.  Cancellation visits each job's immutable item
vector once.  No new scheduler or lifecycle phase was introduced.  This closes
the tested startup/cooperative-shutdown boundaries; physical OOM, actual
geometry cancellation latency, native hosts and shared-stack sanitizers remain
unqualified.  No numeric rendering policy changed and the full Generic/Lucy
matrix was not rerun.  Its earlier failures remain open, and historical passes
do not qualify these binaries.  S1--S6 remain incomplete.

## September 5: atomic source submission under allocation failure

`CXX-SOURCE-002` fixes allocation after caller resources had been consumed.
Queue growth could publish a partial batch; allocating the returned job handle
could fail after publishing the whole batch or consuming a constrained batch.
Submission now prepares an unattached handle and queue storage first.  Queue
failure rolls back only the new suffix while holding the publication lock.
An allocation failure returns an empty handle with every caller resource and
stream preserved; successful admission denial still returns CONSTRAINED.

`/tmp/obol-source-submit-transaction-20260905` retains original sources/library,
the instrumented negative library and executable, the final patch, writer
inventory and validation.  `negative.log` fails all three injected cases:
partial queue append, normal handle allocation and constrained handle allocation.
Each loses caller ownership; normal cases leave queue sizes 2/1 and 3/1 relative
to pre-existing work.  The negative library is `e6b73d11...`, containing only
test instrumentation on predecessor `1b8fc9ad...`.

The final regression preserves source/database/stream/callback identity,
uncancelled streams, pre-existing queue entries and active reservations.
It rejects unintended callbacks, retries the preserved requests, checks the
normal or constrained result and waits for all reservations to retire.
`fixed.log` passes in 0.34 s.  The final 17 focused CTests pass in 13.49 s,
including the source regression, existing drawing/cache fault checks and
conformance catalog (`tests-focused.log`).  GED and Qt drawing integration
tests both pass in 4.56 s total (`tests-integration.log`).

Affected consumers and qged were rebuilt.  `hashes.json` and `ldd.txt` identify
qged `5ffd9e63...`, libBObol `43e78f04...`, libqtcad `6a1a9193...`, Obol
`f281b59b...` and OSMesa `097aac27...`.  External renderer libraries and C++
layouts are unchanged.  These are Debug tests without sanitizers; deterministic
exceptions exercise allocation-failure handling, not physical OOM magnitude.

The existing control-lifecycle composition and deferred-autoview models pass;
`focused-model-audit.json` confirms unchanged tool identity, state counts and
depths.  Catalog lint passes all 45 pairs.  These models check downstream control
contracts; allocation mechanics are outside their abstraction and are verified
by the C++ fault matrix.  No model, conformance mapping or baseline was changed
or accepted.  The source-demand checkpoint remains the latest complete
unfiltered 45-model run.

The source inventory records 23 explicit lifecycle/queue/reservation/stop
mutations, with initialization, copied inputs and bypass aging classified
separately.  Full source registry/adoption, LoD service, GED, host and
presentation audits remain open.  Pool-construction failure, shutdown under
pressure and shared-stack sanitizers remain S4/S6 qualification work.
No numeric rendering policy changed, and the full Generic/Lucy matrix was not
rerun.  Its earlier zoom, prominent-floor and HUD failures remain open; no
historical graphical pass qualifies these binaries.  S1--S6 are incomplete.

## September 5: static-quality configuration lifecycle

`/tmp/obol-static-lifecycle-audit-20260905` retains the two negative builds,
their sources and libraries, the final implementation patch, writer inventory
and validation.  These are focused configuration fixes; the system remains
unqualified for production.

| Counterexample | Before | Fixed behavior |
|---|---|---|
| `CXX-STATIC-002`: repeat automatic submission on a settled constrained view | `autosubmit-negative.log`, predecessor libBObol `67f1a5c6...`: render request appears and work mask becomes 32; constraint mask 1 and allocation serial 17 stay unchanged | Identical normalized intent returns without resetting or resynchronizing automatic work. The regression preserves quiescence, constraint evidence and the allocation serial. |
| `CXX-STATIC-003`: expand the stable frame allowance | `deadline-negative.log`, intermediate library `24dab975...` containing only the first fix: policy stays 3 and no frame is requested | Changed deadlines advance policy, mark work pending and request a capacity frame. Normal settling produces an exact three-occurrence frame with more actual triangles than before the change. |

The expanded deadline is 500 ms in the controlled-cost fixture.  It is a test
input, not a changed default or latency qualification.  The same-value deadline
guard remains.  The existing policy owner retires stale static/retained evidence;
no new field, scheduler or trial state was introduced.  The audit enumerates
19 direct static-trial mutations in `static-writers.json`; the deadline setter
uses the existing policy-revision writer indirectly.

All 14 focused CTests pass in 10.34 s (`tests-focused.log`).  After merging two
adjacent test guards and rebuilding qged, the affected LoD test passes again in
9.14 s (`test-final.log`).  Qt controller and window-host offscreen tests pass
in 0.74 s total.  The private-Xvfb System GL window-host CTest passes in 0.36 s,
without a skip.  A direct private-Xvfb controller invocation also exits zero;
its internally conditional GL framebuffer branch was not separately attested.
These are Debug tests with sanitizers off, not full qged workflows or hardware
GPU qualification.

`hashes.json` and `ldd.txt` identify qged `d4e3345d...`, libBObol `1b8fc9ad...`,
libqtcad `6a1a9193...`, Obol `f281b59b...` and OSMesa `097aac27...`, as well as
the saved negative libraries and test executables.  External Obol/OSMesa are
unchanged.  Affected consumers were rebuilt; no shared-library layout changed.

Static-quality, capacity-search and LoD composition models pass.
`focused-model-audit.json` confirms unchanged tool identity, state counts and
depths.  Catalog lint passes all 45 model/config pairs.  Only a production
mapping comment and regression metadata changed; no transition or baseline
was altered or accepted.  The source-demand checkpoint retains the latest
complete unfiltered 45-model run.

The full Generic/Lucy matrix was not rerun for these setter changes.  The
[preceding shared-minimum rows](#september-5-shared-minimum-coverage-under-source-replacement)
retain the software zoom, System GL prominent-floor and HUD failures; their
historical passes do not qualify these binaries.  Static start/interruption
policy is now inventoried but remains an S3 extraction boundary.  The other
writer/release-edge audits, acceptance thresholds and S1--S6 gates remain open.

## September 5: shared minimum coverage under source replacement

`CXX-COVERAGE-001` fixes a shared-asset repair path which rejected a rich resident
prefix before trying its minimum.  The controller then settled the appended
visible occurrence as a box.  In the controlled reproduction, repair has 449
cost units remaining and the shared minimum costs 74.  Adding structural repair
to the existing minimum-prefix selection rule removes that false rejection.
No new control state, scheduler or mode was added.

`/tmp/obol-constrained-source-allocation-20260905` retains predecessor sources
and library, the negative test source/executable, loader verification and logs.
The strengthened test passes with the fix and fails with the same executable
loading predecessor libBObol `30142145...`: the appended occurrence still has a
12-face box and no progressive mesh.  `qualified-negative.log` and the final
`final-negative.log` retain that replay and its clarified diagnostic.

The shared fixture covers ample and constrained allocation.  The constrained
case uses the real 19,602-face progressive asset at scales 1/2/3, 3,000 offscreen
entries followed by the retained visible entries, actual OSMesa traversal
reports and a fixed 100,000 ns per reported cost unit.  Actual render-start
identity is preserved.  This is controlled numeric/effect evidence, not a
wall-clock performance measurement.  Assertions cover:

- Demand and importance are pending while the visible suffix still has stale
  projection; source input and normals repair overlap that census.
- Current projections preserve uniform scale; current allocated cuts increase
  with importance and at least one retained occurrence gives up detail.
- Three actual mesh occurrences and their selected triangle counts reach the
  renderer, with flat presentation normals and no terminal/structural boxes.
- Selected cost is positive and within its certificate, below pixel demand.
  The journal is continuous, states conform, an exact current frame follows
  demand retirement and precedes the final allocation, and all work retires.

The diagnostic positive run records cuts 8/10 before replacement and 5/7/8
afterward, with 1,550 selected/certified units against 3,148 demanded units.
These numbers describe that run; the assertions check their relationships.
`positive-trace-1.log` preserves the full development trace.  The final journal
starts at source replacement; its checks do not claim the separate asynchronous
result-publication transaction, which shared reuse need not create.

qged and affected consumers were rebuilt.  All 14 focused CTests pass in 9.54 s
(`tests-final.log`), Debug without sanitizers.  Current hashes are qged
`25751c71...`, libBObol `67f1a5c6...`, Obol `f281b59b...` and OSMesa `097aac27...`;
`hashes.json` and `ldd.txt` retain their full identity and resolution.  The
external renderer libraries are unchanged.

The structural-frontier, terminal-quality ordering and LoD composition models
pass; `focused-model-audit.json` confirms unchanged tool identity, state counts
and depths.  Catalog lint passes all 45 model/config pairs.  Only conformance
mapping/comment metadata changed; no model transition or baseline was altered.
An initial sandboxed TLC attempt could not open its local RMI socket; the
successful checks are under `tla-checked`.  The previous source-demand run
remains the complete unfiltered 45-model evidence.

| Shared-minimum graphical artifact | Result |
|---|---|
| `/tmp/qged-obol-shared-minimum-generic-20260905` | Six full rows and camera checks pass. Cold System GL shaded fails HUD fill; warm shaded is skipped. Warm-wire strict trace passes, but its 79.449 ms rotation peak exceeds the 52.5 ms responsive threshold, leaving the earlier responsive retained-cut failure open. |
| `/tmp/qged-obol-shared-minimum-osmesa-20260905` | All 217 events complete in 133.297 s; 2,447 strict control transitions pass. The full warm row fails continuous zoom from effective cut 25; in-4/in-8/in-12 present effective cuts 25/25/23, and active refinement remains at 25. |
| `/tmp/qged-obol-shared-minimum-glsl-20260905` | All 217 events complete in 98.482 s; 1,056 strict control transitions pass. The full opt-in System GL row fails the close-view prominent floor at cut 28, 2,619,533 faces and normalized error 3.5813, target 0.25. |

Software normals checks present equal cut 24, 2,101,208 faces and 6,303,624 normals
at all four checkpoints.  Its inspected close view reaches cut 29, 3,445,998
faces, no explicit normals and normalized error 2.8587 at target 0.25, without a
prominent-floor violation.  The software smooth-normal image and failed Generic
HUD image were also inspected.  `checkpoint-diagnosis.json` records exact
populations and limits.  System GL's inspected close image still fails its floor;
its flat/authored-reference/smooth normals checkpoints use cut 24 and 2,101,208
faces, but authored-return drops to cut 22 and 861,234 faces.  It supplies no
equal-cut four-way normals comparison.  Runtime provenance matches the hashes
above for all three graphical roots; validators and thresholds are unchanged.
Adaptive differences have not been attributed to this fix.  S1/S2 remain incomplete, and the failed full rows remain open S4/S5
qualification.  The private-Xvfb System GL context is CPU llvmpipe, not hardware
GPU evidence.  No native-platform, sanitizer or main-baseline claim is made.

## September 5: completed traversal measurement

`CXX-TIMING-001` preserves the measured duration supplied with an exact
traversal.  The old path sampled the clock again at delivery: the regression's
41 ms traversal plus 25 ms delivery delay became 66.001017 ms.  The fixed
completion consumer records exactly 41 ms.  On-time completion and the existing
incomplete-frame behavior also pass.  `renderToImage` now sends the measurement
already taken at traversal completion to the same consumer.

`/tmp/obol-completed-frame-duration-20260905` retains the predecessor sources
and library, failing test source/executable, `negative.log`, and
`implementation.patch`.  qged and the affected build dependencies were rebuilt;
all 14 focused CTests pass in 9.76 s (`tests-fixed.log`).  These Debug tests do
not use sanitizers.  `hashes.json` identifies qged `25751c71...`, libBObol
`30142145...`, Obol `f281b59b...` and OSMesa `097aac27...`.

The timing-evidence, CAD-frame composition and LoD composition models pass;
`focused-model-audit.json` confirms their tool identity, status, state counts
and depths match the accepted baseline.  Catalog lint passes for all 45 pairs.
The timing model adds a mapping comment and production test metadata; its
transition relation is unchanged.  No baseline was accepted here.  The preceding
source-demand run remains the complete unfiltered 45-model evidence.

| Completed-frame timing artifact | Result |
|---|---|
| `/tmp/qged-obol-frame-duration-generic-20260905` | Six complete rows and their camera checks pass: four OSMesa shaded/wire cold/warm rows and two System GL wire rows. Cold System GL shaded still fails HUD fill; warm shaded is skipped. |
| `/tmp/qged-obol-frame-duration-osmesa-20260905` | Full warm row and camera contract pass. All 217 events complete in 83.345 s with 1,364 valid control transitions. Continuous zoom starts at effective cut 23 and reaches 24 during active refinement. |
| `/tmp/qged-obol-frame-duration-glsl-20260905` | Opt-in System GL completes 217 events in 85.104 s with 943 valid control transitions. The full row fails the close-view prominent floor: cut 28, 2,619,533 faces, normalized error 3.5813. |

Both Lucy traces and Generic warm-wire pass independent strict validation.
Generic warm-wire's 55.251 ms rotation peak exceeds its 52.5 ms responsive
threshold; the allowed ceiling remains zero.  This pass does not exercise the
earlier failure with a responsive peak and an unjustified ceiling.  The failed
HUD image still has no progress fill and was inspected.

Software Lucy's four normals checkpoints present the same cut 21, 559,494 faces
and 1,678,482 normals.  Its inspected close view presents 2,619,533 faces at cut
28 with no floor violation, target error 0.5 and normalized error 1.7906.  The
System GL close view has the same face count but target error 0.25; retain those
different policy conditions when interpreting its floor failure.  Software
normals and both close images were inspected; `checkpoint-diagnosis.json`
records exact populations and active/requested presentation limits.

Runtime provenance matches the hashes above.  The renderer libraries,
validators, deadlines, quality floors and warm cache are unchanged.  These
private-Xvfb System GL results use CPU llvmpipe, not a hardware GPU; no main
comparison or sanitizer qualification was performed.  All process handles are
terminal.  The software pass from cut 23 does not close the earlier cut-25 zoom
failure, and adaptive differences have not been causally attributed to this
fix.  The retained-detail and HUD failures remain open as well.

The constrained source-allocation proof remains open.  Its next test can use
the existing completed-frame delivery API with controlled duration while
checking actual projected geometry and presentation prerequisites.  S1--S6
remain incomplete.

## September 5: demand preservation through source coverage

`CXX-DEMAND-002` preserves pending full-view demand across source input and
normals repair.  The existing coverage owner distinguishes source coverage
from a view census, so preserving demand does not permit optional refinement
ahead of new geometry's minimum coverage.  The conformance finding owns the
implementation boundary and remaining importance-allocation proof.

`/tmp/obol-source-demand-supersession-20260905` retains source baselines,
`implementation.patch`, and separate failing inventory/visibility probes.
`before-fix` contains the failing test source, executable and matching
libBObol.  The first visibility probe had an incorrect expected hidden-entry
count; `visibility-demand-counted-probe.log` is the corrected negative result.
Both corrected probes complete the expected repair population but lose demand.
The fixed automatic cases preserve it and retire after one ordinary successor.

All 14 focused CTests pass in 8.76 s after rebuilding the affected consumers;
`build-final.log` and `tests-final.log` retain the commands' results.  The
one-mesh fixture now finishes its pending demand and source-coverage setup
before checking unchanged requests and scale-only retargeting.  Its original
behavior assertions remain.  These Debug runs do not use sanitizers.

The expanded `ObolSubmissionPass` component passes with 135 generated and
66 distinct states at depth eight.  The full 45-pair catalog lint passes.
`model-negative` rejects the old coverage-completion wakeup with
`DemandHasProducer`; `model-consumed-negative` rejects premature demand
consumption with `DemandRetiresOnlyAfterDemandPass`.  The full unfiltered
45-model suite passes in `tla-full` (2,237 summed model seconds).  Its accepted
baseline was audited against the saved predecessor: tool identity, model
catalog and all other status/state-count/depth results match.  Only the
expanded submission model changes, from 90/45/depth 11 to 135/66/depth 8.
`baseline-audit.json` retains the comparison.  All formal processes exited;
this formal result does not qualify the graphical failures below.

`hashes.txt`, `ldd.txt` and GUI runtime provenance identify qged `2b28af79...`,
libBObol `7f48a11a...`, Obol `f281b59b...`, and OSMesa `097aac27...`.
The external renderer libraries are unchanged.

| Source-demand checkpoint artifact | Result |
|---|---|
| `/tmp/qged-obol-source-demand-generic-20260905` | Four OSMesa cold/warm shaded/wire rows and cold System GL wire pass, including camera checks. Cold System GL shaded fails HUD fill; warm shaded is skipped. Warm System GL wire fails responsive retained-cut preservation. |
| `/tmp/qged-obol-source-demand-osmesa-20260905` | All 217 events complete in 130.070 s with 1,915 valid control transitions. Coverage and equal-cut normals pass; the full warm row fails continuous-zoom responsiveness. |
| `/tmp/qged-obol-source-demand-glsl-20260905` | Opt-in System GL completes 217 events in 83.930 s with 946 valid control transitions. The full warm row fails the close-view prominent floor. |

Both Lucy traces and the failed Generic warm-wire trace pass independent
strict validation.  Generic warm-wire rotation peaks at 48.465 ms, within its
52.5 ms acceptance limit, but the held state installs ceiling zero instead of
retaining an unrestricted cut.  `retained-cut-diagnosis.json` records this
predicate failure and its image was inspected.  The failed shaded
`subpath-redraw-return.png` again has no HUD track/fill and was inspected.

Software Lucy's flat/authored/smooth/authored-return checkpoints all present
cut 24, 2,101,208 faces and 6,303,624 normals.  Authored and smooth images were
inspected.  Continuous zoom starts at effective cut 25; the in-4/in-8/in-12
checkpoints have effective cuts 25, 25 and 23, with no improvement over the
starting cut at those checkpoints.  Its close view has no prominent-floor
violation and normalized error 2.8587.  System GL's close view presents cut 28,
2,619,533 faces, normalized error 3.5813 and one floor violation.  Both close
images were inspected; `checkpoint-diagnosis.json` in each case retains exact
populations and separate active/requested/presented facts.

Five complete Generic rows pass at this checkpoint.  The preceding checkpoint's
six Generic passes and full software Lucy pass remain dated evidence; they do
not qualify these binaries or prove repeatability.  Preserve all earlier
normals, zoom, retained-cut, planning-cycle and frame-delivery reproductions.
No causal attribution of the adaptive graphical differences to this fix is
established.  Validators, deadlines, quality floors and the shared warm cache
are preserved.  Private-Xvfb System GL is CPU llvmpipe, not hardware GPU
qualification; no main-baseline comparison was run.  S1--S6 remain incomplete.

### Visible demand replacement follow-up

`/tmp/obol-source-importance-20260905` adds a production integration test on the
same libBObol `7f48a11a...` stack.  Two retained tetrahedra have different
projected sizes; 3,000 off-screen entries make the census interruptible.
After rotation, the test verifies pending demand and importance, appends a
third visible size, changes normals, and pumps the actual software renderer
to an exact frame with no remaining control obligation.

The test checks current view/policy stamps, changed projection after rotation,
linear projected extent under uniform orthographic scaling, four faces per
occurrence, and flat normals in the geometry supplied to the renderer.  A
payload's cached producer normal style is not the presentation-style authority.
All 14 focused CTests pass in 8.82 s (`tests-qualified.log`); the update-action
test takes 4.23 s.  `projection-test.patch` and `before` retain the change and
its predecessor.  No production library or graphical validator changed.

This ample-capacity fixture settles without a retained allocation certificate.
It proves projection/presentation replacement and finite retirement, **not**
redistribution of a constrained budget or the exact frame that authorizes that
allocation.  That effect proof remains open in active debt.  Earlier probe
logs record corrected fixture assumptions (camera observation before debounce,
retained meshes before the large census, optional allocation, and cached versus
presented normals); they are not confirmed production counterexamples.

## September 5: demand preservation through presentation repair

The submission audit classifies all 51 direct pass mutations and the indirect
capacity/repair entries.  It exposed `CXX-DEMAND-001`: completing a normals
repair could consume an ordinary full-view demand obligation, although repair
can rebind an existing cut without evaluating richer demand.  The caller now
includes presentation repair in the existing demand-completion guard.  No
state, revision, pass mode or scheduler is added.

`/tmp/obol-submission-owner-audit-20260905` retains the initial sources,
`submission-writers.json`, isolated `implementation.patch`, and
`repair-demand-probe.log`.  The failed automatic-controller probe visits all
3,014 repair entries and loses demand; `before-fix` retains its test source,
executable and matching old libBObol.  The fixed regression preserves demand
through repair, visits one complete ordinary successor and proves terminal
demand/cursor retirement.  It uses the existing compact fixture with no worker
tasks or cache writes.

`build-fix.log` and `build-final.log` record successful rebuilds, including qged
and the affected tests.  All 14 focused CTests pass in 8.82 s in
`tests-final.log`; these Debug checks run without sanitizers.  Full catalog
lint passes for all 45 model/config pairs.  Focused TLC checks of
`ObolSubmissionPass`, `ObolTerminalConvergenceComposition`, and
`ObolControlLifecycleComposition` pass with baseline comparison enabled.
Their status, tool identity, state counts and depths match the accepted
baseline.  The model transition relation is unchanged; the new conformance
entry maps a production regression to its existing demand-preservation rule.
The [source-evidence checkpoint](#september-5-source-observation-and-consumption)
remains the latest unfiltered 45-model pass.  No baseline was accepted here.

`hashes.txt`, `ldd.txt` and GUI runtime provenance identify qged `2b28af79...`,
libBObol `0835efd7...`, Obol `f281b59b...`, and OSMesa `097aac27...`.
The external renderer libraries are unchanged.

| Demand-repair checkpoint artifact | Result |
|---|---|
| `/tmp/qged-obol-demand-repair-generic-20260905` | Four default OSMesa cold/warm shaded/wire rows and both System GL wire rows pass, including camera checks. System GL cold shaded still fails HUD fill; warm shaded is skipped. |
| `/tmp/qged-obol-demand-repair-osmesa-20260905` | The full default OSMesa warm row and camera contract pass: 217 events in 120.576 s, 2,244 valid control transitions. |
| `/tmp/qged-obol-demand-repair-glsl-20260905` | Opt-in System GL completes 217 events in 82.376 s with 983 valid control transitions. The full warm row fails the close-view prominent floor. |

The failed Generic shaded row also passes independent strict control-trace
validation.  Its inspected `subpath-redraw-return.png` still has no HUD
track/fill despite reported geometry.  Software Lucy passes the same-cut
normals comparison: flat/authored/smooth/authored-return all present cut 24,
2,101,208 faces and 6,303,624 normals.  Its continuous-zoom test starts at
effective cut 25 and reaches 26 during input.  The close checkpoint presents
3,874,829 faces at ceiling 29, normalized error 2.8587 and no prominent-floor
violation.  Authored/smooth and close captures were inspected;
`checkpoint-diagnosis.json` retains the measured populations.

System GL Lucy's close checkpoint presents cut 27 with 1,806,720 faces,
normalized error 4.3543 and one prominent-floor violation.  Zoom-out presents
cut 25 with no floor violation on this run.  Its independent strict control
trace passes, and the close image was inspected.  This row used
`OBOL_CAD_SOFTWARE_GLSL=1`; its report and `checkpoint-diagnosis.json` retain
the exact measurements.  All build, model and GUI handles are terminal.

The software pass qualifies this run.  It does not establish that the repair
fix caused every adaptive difference or retire earlier normals, cut-25-start,
planning-cycle and retained-cut reproductions without repeatability evidence.
Defaults, validators, deadlines, quality floors and the shared Lucy warm cache
are preserved.  Private-Xvfb System GL uses the existing CPU llvmpipe
environment; it is not hardware GPU qualification.  No main-baseline comparison
was run.  S1 still requires the identified source-supersession effect proof
and remaining state-family audit; S1--S6 are not declared complete.

## September 5: submission retarget prefix

`submitLodRequestsIfNeeded` now distinguishes an untouched submission cursor
from a consumed prefix when inputs change.  It rebuilds the untouched
source-local plan without inventing a second whole-scene pass.  A consumed
prefix retains its cursor and owes one full rescan; independently owed
selective-delta rescans remain intact.  The fix adds no state or scheduler.

`/tmp/obol-submission-rescan-20260905` retains the initial files, isolated
`implementation.patch`, and `auto-retarget-probe.log`.  That automatic
production-controller probe fails before the fix; `before-auto-fix` retains
its source, test executable and old libBObol.  The earlier manual-controller
probes do not exercise the automatic lifecycle and are not negative evidence
for this finding.  The final test covers policy and scale changes before a
pass starts, and a policy change after a real prefix.  It checks exact visits
and completion over the existing 3,014-occurrence fixture without provider work.

`build-fix.log` records the successful 42-target rebuild.  `tests-final.log`
records all 14 focused CTests passing in 8.80 s.  Final review moved a fixture
comment and corrected a failure message; the rebuilt affected test passes
again in `test-review.log`.  These Debug checks run without sanitizers.

Full catalog lint passes for all 45 model/config pairs.  Focused TLC checks of
`ObolSubmissionPass`, `ObolSourceEvidence`, `ObolStaticQuality`,
`ObolTerminalConvergenceComposition`, and `ObolControlLifecycleComposition`
all pass with baseline comparison enabled.  Per-model `results.json` files
and `tla-focused.log` retain the results.  Only the submission model's boundary
comment and conformance/catalog mapping changed; its transition relation is
unchanged.  The source-evidence checkpoint remains the latest unfiltered
45-model pass, and no baseline was accepted by this checkpoint.

`hashes.txt`, `ldd.txt` and each GUI runner's runtime provenance identify qged
`2b28af79...`, libBObol `dd20abf3...`, Obol `f281b59b...`, and OSMesa
`097aac27...`.  External renderer libraries are unchanged.

| Submission-retarget checkpoint artifact | Result |
|---|---|
| `/tmp/qged-obol-submission-rescan-generic-20260905` | All four default OSMesa cold/warm shaded/wire rows and both System GL wire rows pass, including camera checks. System GL cold shaded still fails HUD fill; warm shaded is skipped. |
| `/tmp/qged-obol-submission-rescan-osmesa-20260905` | Default OSMesa completes 217 events in 134.420 s with 2,353 valid control transitions. The full row fails the same-cut normals comparison. |
| `/tmp/qged-obol-submission-rescan-glsl-20260905` | Opt-in System GL completes 217 events in 100.968 s with 1,034 valid control transitions. The full row fails the close-view prominent floor. |

The failed shaded Generic Twin and both Lucy rows pass independent strict
control-trace checks.  Current wire passes support the retarget regression;
they do not prove every historical planning-cycle or retained-cut timing
condition is resolved.  The shaded `subpath-redraw-return.png` has no HUD
track/fill on inspection despite reported fill geometry.  Its next diagnostic
boundary is correspondence between the captured presented frame and the
reported faceplate state; keep the validator unchanged while establishing it.

Software Lucy's authored reference and return present cut 24 with 2,101,208
faces; smooth presents cut 23 with 1,284,926 faces.  All are terminal and ready
at policy revision 39, inventory 2 and availability 2.  Deadline interruptions
increase across these events, but causality is not established.  The captures
were inspected and `checkpoint-diagnosis.json` retains the measurements.
The same-cut requirement remains open; unequal populations do not establish
normal-policy image correctness.  This failure prevents the runner from
reaching its later continuous-zoom check, so the preceding cut-25-start
failure remains open.

System GL Lucy's close view reports normalized error 3.5813 and one prominent
floor violation.  Its zoom-out checkpoint has no floor violation on this run;
earlier zoom-out failures remain open.  The close image was inspected and
`checkpoint-diagnosis.json` retains the measurements.  These private-Xvfb
System GL runs use the existing CPU llvmpipe environment and do not qualify
hardware GPU rendering.  Lucy used `OBOL_CAD_SOFTWARE_GLSL=1`; defaults,
validators, deadlines, floors and the shared warm cache were preserved.
No main-baseline comparison was run.  All model, build and GUI processes are
terminal.  The bounded retarget defect is closed at its tested boundary;
the remaining S1 writer/successor audit and S4/S5 failures are still required.

## September 5: static completed-frame ownership

The existing static-quality trial now selects preparation/steady replay,
fractional acceptance, predicted rejection, cut advancement and population
handoff in independently compiled `lod_static_quality.cpp`.  The controller's
competing completed-frame branches are removed.  No persistent state or
geometry storage is added.  `CXX-STATIC-001` additionally corrects an
unavailable single-occurrence prediction being published as an unlimited
renderer ceiling: it transfers the measured budget to the existing occurrence
allocator, retaining the guard and avoiding a fabricated rejection.

`/tmp/obol-static-frame-owner-20260905` retains the initial source files,
isolated `implementation.patch`, failing `prediction-probe.log`, and pre-fix
policy/test/executable.  `build-final.log` records the successful 84-target
completion after fixing the integration test's private policy link inputs.
`tests-final.log` records all 14 focused CTests passing in 8.97 s.  The new
cases drive the compiled production decision, coordinator admission effects,
and real two-progressive-occurrence cost provider; existing lifecycle and
handoff regressions also pass.

Full catalog lint passes for all 45 model/config pairs.  Focused TLC checks of
`ObolStaticQuality`, `ObolCapacityPresentationHandoff`,
`ObolTerminalQualityOrdering`, `ObolTerminalConvergenceComposition`, and
`ObolControlLifecycleComposition` all pass with `TLA_COMPARE_BASELINE=ON`.
Their per-model `results.json` files and `tla-focused.log` retain the results.
No TLA transition relation changed in this extraction.  These checks do not
model the numeric sentinel or replace the complete formal-suite gate; the
source-evidence run remains the latest unfiltered 45-model pass, and the
accepted baseline is unchanged by this checkpoint.

`hashes.txt` and `ldd.txt` identify qged `2b28af79...`, libBObol
`73227704...`, Obol `f281b59b...`, and OSMesa `097aac27...`.
Every graphical runner records matching runtime provenance.

| Static-frame checkpoint artifact | Result |
|---|---|
| `/tmp/qged-obol-static-frame-generic-20260905` | All four default OSMesa cold/warm shaded/wire rows and their camera checks pass. System GL cold shaded fails HUD fill; cold wire fails the planning-cycle trace check. Both System GL warm rows are skipped. |
| `/tmp/qged-obol-static-frame-osmesa-20260905` | Default OSMesa completes 217 events in 125.962 s with 2,264 valid control transitions, but fails continuous-zoom refinement. |
| `/tmp/qged-obol-static-frame-glsl-20260905` | Opt-in System GL completes 217 events in 87.845 s with 948 valid control transitions, but fails the close-view prominent floor. |

The failed shaded Generic Twin and both Lucy rows pass independent strict
control-trace checks.  The wire failure is a distinct current reproduction:
planning observations 75--77 switch SUBMISSION -> SUBMISSION + SUBMISSION_RESCAN ->
SUBMISSION, moving the cursor 0 -> 2048 -> 0 at unchanged revisions,
transaction 14 and completed frame 19.  The script itself completes 72 events;
`control-cycle-diagnosis.json` retains the surrounding transitions.  This trace
prompted the [submission-retarget audit](#september-5-submission-retarget-prefix)
recorded above.  The trace alone does not prove an infinite loop or explain
every successor.  The prior wire retained-cut failure remains unresolved as
well; this different run does not retire it.

The shaded `subpath-redraw-return.png` again lacks the reported HUD fill,
confirmed by inspection.  Software Lucy starts at cut 25 and exposes effective
in-gesture cuts 25, 23 and 25.  Its resident prefix grows and cache loads
increase, but active/requested cut both equal 27, so the existing discrete-cut
exception does not apply.  `continuous-zoom-diagnosis.json` retains the actual
values.  System GL Lucy's close frame is cut 27 with 1,806,720 presented faces
and normalized error 4.3543, one prominent-floor violation.  Zoom-out has no
such violation on this run; earlier failed zoom-out evidence remains open.

Both Lucy close captures were inspected.  All build, model and GUI handles
are terminal.  `glxinfo.txt` records Mesa 25.2.8 llvmpipe on the CPU, so these
runs do not qualify hardware GPU rendering.  The System GL Lucy run used
`OBOL_CAD_SOFTWARE_GLSL=1`; defaults, validators, deadlines, quality floors and
the shared warm cache were preserved.  No main-baseline comparison was run.
These results close the tested static decision defect, while the graphical
failures and S1--S6 gates remain open.  Their causal relationship to the
extraction is not established.

## September 5: source observation and consumption

Notifications and submission now use one `BObolLodSourceEvidence` owner.
Distinct pending observations publish their exact admission domain; duplicates
and consumption do not publish another revision.  A stale submission cannot
consume a newer observation.  The former unkeyed flag and competing controller
publication branches are removed.  Immutable source signatures preserve the
sparse submitted baseline without copying occurrence or geometry storage.

Final evidence is under `/tmp/obol-source-evidence-20260905`.  Its isolated
`implementation.patch` separates this change from earlier uncommitted work.
`build-final.log` records the successful 82-target rebuild, and
`tests-final.log` records all 14 affected CTests passing in 8.70 s.  The
production-owner and real-controller regressions cover pending change order,
duplicate delivery, stale consumption, replacement/removal/reappearance,
notification-free submission and empty-source retirement.

The complete pinned TLC suite passes all 45 models; `tla-full/results.json`
and the accepted baseline record the results.  `ObolSourceEvidence` reaches
1,408 distinct states at depth nine.  Restoring the old pending-publication
behavior in the isolated `model-mutant` copy violates
`ObservedInputIsDelivered` after two pending inventory observations (exit 12,
counterexample retained in `model-mutant/tlc.log`).

`hashes.txt` records final qged `2b28af79...`, libBObol `a98f3c9a...`, Obol
`f281b59b...` and OSMesa `097aac27...`; `ldd.txt` records resolution.

| Source-evidence checkpoint artifact | Result |
|---|---|
| `/tmp/qged-obol-source-evidence-generic-20260905` | All four default OSMesa cold/warm shaded/wire rows and their camera checks pass. System GL cold shaded fails HUD fill; cold wire fails responsive retained-cut acceptance. Both System GL warm rows are skipped because cold qualification failed. |
| `/tmp/qged-obol-source-evidence-osmesa-20260905` | Default OSMesa completes 217 events in 141.763 s with 2,459 valid control transitions, but fails continuous-zoom acceptance. |
| `/tmp/qged-obol-source-evidence-glsl-20260905` | Opt-in System GL completes 217 events in 76.088 s with 911 valid control transitions, but fails the prominent quality floor at close and zoom-out checkpoints. |

The failed System GL Generic Twin rows and the System GL Lucy row also pass
independent strict control-trace checks, retained beside their reports.
Passing control traces do not change their full-row failures:

- The Generic Twin shaded `subpath-redraw-return.png` lacks the progress fill
  which its sample reports present; inspection confirms the discrepancy.
- Generic Twin wire retains ceiling 0 during rotation despite meeting the
  validator's responsive-cut condition.  The pre-rotation frame is 48.383 ms
  and the rotation peak 50.192 ms, near the 50 ms presentation deadline and
  below the validator's 52.5 ms comparison.  Reconcile the evidence and policy
  at their owning boundary; the validator and deadline remain unchanged.
- Lucy OSMesa starts at cut 25.  Its sampled effective in-gesture cuts are
  25, 23 and 25, although the resident prefix grows to cut 28.  The current
  requested cut also equals 28, so the existing discrete-limit exception does
  not apply.  `continuous-zoom-diagnosis.json` retains these values.  All 42
  wheel callbacks meet 250 ms (peak 2.120 ms); callback latency is not the
  failed condition.  This preserves the earlier cut-25-start reproduction.
- Opt-in System GL Lucy has one prominent-floor violation at both close and
  zoom-out: normalized errors 3.5813 and 3.2554 respectively.  A terminal
  ownerless state with performance-limit evidence does not excuse that floor.

Representative OSMesa Generic Twin final shaded/wire images and both Lucy
close-view captures were inspected.  All jobs exited.  `glxinfo.txt` records
the same private-X-server setup's Mesa llvmpipe renderer; these are CPU runs,
not native hardware-GPU qualification.  System GL Lucy used
`OBOL_CAD_SOFTWARE_GLSL=1`; defaults, deadlines, quality floors and reusable
caches were preserved.  No main-baseline comparison was available.

This closes the tested source identity defect, with the graphical failures
retained under S5.  Their causal relationship to the source change is not
established.  S1/S2 and full release qualification remain incomplete; the
remaining lifecycle, scale, rendering, editing and platform requirements are
still required.  Earlier passing Lucy rows do not qualify these new binaries.

## Dated capability evidence

| Area | Dated evidence | Recorded status for the cited runs |
|---|---|---|
| focused CTest gate | After a complete current-tree Ninja build on 2026-09-01, all 28 `bobol_headless` tests pass in 9.0 seconds, including renderer, LoD service/coordinator/update, cache, compact ownership, retained allocation, edit manipulator, host, API/symbol, GED draw-sync, and view-command contracts.  The policy-disable regression additionally proves that off/on/off transitions preserve the current presentation, retire every automatic owner, keep explicit manual generations usable, and cannot be rearmed by a renderer-style change or capacity-labelled repaint.  The service test proves that a coalesced asset producer may complete against its retained latest demand; the view-state test proves that a superseded result cannot create a terminal occurrence failure, while genuinely stale source/cache data retains its failure semantics.  Independent draw scope is exercised with a real local endpoint across draw/erase/redraw/zap.  The linked Obol CAD suite includes its 131,072-occurrence classifier-reservation regression.  The current OSMesa/offscreen interaction gate passes all 13 qged event, measurement, settings, polygon, polygon/sketch, primitive-edit, and framebuffer-host tests in 14.7 seconds.  After explicitly relinking the two model-test executables against the current controller ABI, all 19 qtcad Pinewood, Havoc, M35, NIST, and Generic Twin real-model/progressive tests pass in 46 seconds. | broader graphical production suite still required |
| shared clients | Current gsh, MGED, Archer/TkObol, and rtwizard binaries were explicitly rebuilt against the current shared stack.  A fresh private-X-server gate passes all 15 thread-affinity, Tk widget/attach, Archer, gsh draw/ERT/rt-routing/progressive-LoD, MGED host/draw/framebuffer/ERT/progressive-LoD/edit-restore, and rtwizard smoke tests in 135 seconds. | shared-client smoke and core framebuffer/edit routes green; comprehensive interactive client qualification remains |
| qged retained interaction | The final 2026-08-29 binaries pass the focused 12-test qtcad/qged interaction gate and the dual-backend graphical matrices.  `/tmp/qged-selection-presentation-revision-20260829` passes System GL and OSMesa point add/remove, fractional-DPR rectangle placement/commit, selected styling, erase-selected, and redraw-selected.  The 2026-08-31 follow-up at `/tmp/qged-selection-ui-delivery-final-20260831` additionally proves every observed selection frame drains to idle without changing the LoD policy or requesting renderer-capacity evidence.  Single-object ARB/ELL/sketch manipulators and the navigation gizmo now use the typed presentation-only endpoint API, and direct GL explicitly presents the sole exact-frame barrier when Qt would otherwise coalesce it indefinitely.  Full Hubble hierarchy selection/erase/redraw passes on System GL and OSMesa in 9/26 seconds at `/tmp/qged-selection-stability-hubble-{system-final2,final}-20260831`, with the selection frame idle before the next checkpoint and no LoD service work.  A later intermittent quiet-handoff stall was traced to a level-triggered demand refresh whose bounded cursor had been retired by a stronger presentation owner.  The guarded cursor restart now passes a fresh System GL/OSMesa Hubble lifecycle in 11/14 seconds plus two independent OSMesa repetitions at `/tmp/qged-selection-demand-restart*-20260831`; every selected observation preserves the policy revision and returns to IDLE/ready with zero obligations.  After terminal capacity certificates were made to retire predecessor handoff debt atomically, fresh point-selection replays at `/tmp/qged-hubble-point-selection-post-handoff-{system,osmesa}-report.json` pass in 5.2/9.9 seconds: selection creates one exact-presentation obligation, preserves the view/policy revisions, and the following observation is ready with zero control or service debt.  `/tmp/qged-polygon-visual-presentation-revision-20260829` passes retained polygon styling, selection, movement, resize, zoom, cleanup, and cross-renderer placement; immediate selection/move deltas are 1,583/4,894 pixels on System GL and 1,586/4,921 on OSMesa.  `/tmp/qged-framebuffer-presentation-revision-20260829` passes in 37/37 seconds and combines progressive Generic Twin drawing, hierarchy selection, an ellipsoid edit manipulator, faceplate center marker, external raytrace framebuffer underlay/overlay/interlay, resize, rerender, and framebuffer disable.  The current-binary reruns at `/tmp/qged-selection-ui-contract-final-20260901`, `/tmp/qged-polygon-visual-contract-final4-20260901`, and `/tmp/qged-framebuffer-contract-final2-20260901` pass on both renderers; the matching System-GL single/quad widget replays are under `/tmp/qged-system-ui-final-20260901`.  Shared feature/polygon stores now publish content revisions independently of their coalesced controller render latch, so every actual shared mutation wakes the attached views while queries and equal publication remain passive. | focused cross-renderer selection/polygon/framebuffer/edit gate green; full control, physical-pointer, and model matrix remains |
| qged measurement | `qged_measure_ui_replay` and its quad-layout counterpart drive the actual palette and canvas through 2D and exact-hit 3D distance/angle gestures, degree/radian readback, cancellation, resize, and tool replacement.  They pass on OSMesa in the ordinary CTest gate and on final System OpenGL and OSMesa binaries in single and quad layouts; the current System-GL evidence is under `/tmp/qged-system-ui-final-20260901/measure-*`.  The replay asserts line-layer point counts and field values, requires visible framebuffer deltas for both measurement modes, and requires cancellation to reproduce the exact baseline (zero changed pixels).  It exposed and fixed three production defects: no endpoint binding for right-click cancel, palette deactivation that failed to detach its semantic filter, and measurement lines hidden by shaded geometry.  Line-layer depth behavior is now explicit, with measurement guides depth-independent and other diagnostic overlays depth-tested by default. | dual-backend, single/quad, 2D/3D measurement gate green; fractional-DPR physical-device repetition remains |
| qged view settings | `qged_view_settings_ui_replay` and its quad-layout counterpart pass on final System OpenGL and OSMesa binaries in single and quad layouts; the current System-GL evidence is under `/tmp/qged-system-ui-final-20260901/settings-*`.  The replay drives the real settings palette, toolbar, GED commands, endpoint properties, retained Obol features, and framebuffer.  It proves bidirectional ADC/center-dot/grid/model-axes/scale/view-axes/parameter/FPS/framebuffer state, exact cutting-plane command readback, a visible shaded clipping delta, and pixel-exact restoration when clipping is disabled.  Stable replay IDs cover every exercised control.  System GL exposed a host race in which a fast renderer consumed the transient refresh latch before qged's post-command comparison; canvas checkpoints now compare libbv's monotonic frame revision as the durable semantic witness. | dual-backend, single/quad settings and cutting-plane control gate green; fractional-DPR and in-scene plane-affordance qualification remain |
| qged draw modes and resize | `/tmp/qged-resize-moss-all-modes-contract-final2-20260901` passes all 24 moss workflows: modes 0--5, managed LoD and full-detail policy, on System GL and OSMesa.  The apparent mode-5 camera mismatch was a harness race: a delayed native X11 configure changed the window aspect after the resize wait but before draw/autoview.  The runner now requires one second of initial geometry stability and proves the semantic draw batch executes at the requested geometry.  The isolated replacement at `/tmp/qged-resize-moss-m5-framing-final-20260901` passes all four mode-5 rows with identical cross-policy camera size and height.  The current mode-3 rook replacement at `/tmp/qged-resize-rook-m3-contract-final-20260901` passes its four rows in 5--13 seconds.  Each run resizes during realization and after settling, exercises minimize/maximize/fullscreen restore plus a resize storm and policy round trip, and terminates with zero structural boxes and no progressive work.  GUI replay uses an application-owned temporary QSettings scope and never reads or writes the operator's geometry or plugin choices.  Top-level resize requests carry generations and an optional native stability barrier, so a delayed window-manager acknowledgement cannot change autoview or let an old retry resurrect a superseded size; the event-player regression injects that race directly. | baseline all-mode, dual-backend resize/policy matrix green; repeat on larger models and fractional-DPR hardware |
| control-plane models | The focused publication, host-work, convergence, admission, arbitration, canonical-pipeline, cold-preview, static-quality, renderer-preparation, interaction-session, deadline-ownership, capacity/handoff, timing-evidence, resident-growth, point-recovery, structural-frontier, and producer-demand models pass.  The canonical pipeline explored 2,358,764 generated / 1,095,220 distinct states.  The current 29-fact control refinement explored 475,078 / 237,568 states and distinguishes exclusive interrupted replay from ordinary exact-presentation debt, while the lifecycle composition explored 1,092,377 / 227,787 states to depth 32 and proves recovery after two bounded demand-cursor interruptions plus isolation of semantic-only exact frames from capacity recovery.  The cross-boundary composition model explored 35,600 / 18,136 states to depth 15 and proves that admission, independent/capacity-owned growth, exact visibility presentation, ceiling reconciliation, capacity sampling, structural repair, point quality, and terminal publication share one compatible owner/terminal contract.  The terminal-convergence composition explores 1,520 / 1,187 states to depth 81 and requires an exact visibility census and its exact framebuffer classification before reallocation.  The bounded capacity search includes the protected minimum in its immutable numeric domain, freezes capacity-owned allowance inputs, and exports an exact allocation/presentation/three-sample progress rank; its current exploration count is recorded in `libbobol_formal_models.md`.  `ObolCadFrameComposition` now treats camera and effective renderer controls as one exact-frame revision (4,921 / 1,568 states, depth 17).  `ObolHostWork` includes unsubmitted source revisions, timing-driven work opened by frame completion, and provider retirement while an independent exact frame remains; it explores 78,500 / 12,902 states to depth 17.  The 2026-09-01 `ObolSubmissionPass` explores 90 / 45 states to depth 11 and proves that scene-wide demand debt survives a selective source pass until a complete successor consumes it or a newer semantic revision supersedes it.  `ObolLodConvergence` explores 742 / 423 states to depth 25 with the matching quiet-quality owner rule.  On 2026-08-29 `ObolPointTerminalEvidence` added the missing idempotence contract for a constrained point cut: an unchanged exact structural census cannot erase the only terminal witness after its finer-preload decision was consumed.  TLC explored 243 generated / 135 distinct states to depth 15.  The revised admission model explored 321 / 151 states to depth 6 and proves late safe-scene certification may atomically promote an existing PoP payload by marginal cost.  The 2026-09-01 capacity handoff explores 159,114 / 17,598 states to depth 19, requires a resident retry to resolve as drawable or explicitly constrained, rejects a zero-work sample unless an exact census proves the scene is empty, and requires a terminal certificate to retire older reconciliation debt atomically.  Corresponding coordinator, update-action, and dual-renderer Generic Twin regressions pass.  Runtime refinement checks sampled regression, duplicate plans, spontaneous reopen, unwitnessed presentation/constraints, invalid readiness, and six-domain revision monotonicity. | formal composition and sampled-runtime boundaries green; complete event/effect reducer coverage and per-transition records remain |
| Lucy OSMesa | The exact-current 2026-08-26 true-cold shaded lifecycle at `/tmp/qged-lucy-latest-demand-contract-20260826` passes in 76 seconds and separately records globally representative coverage and the first real CAD mesh.  Camera changes made while the immutable hierarchy is being built no longer discard its live pages or final result: the producer resolves page selection against the service-owned latest demand, and service validation uses that same demand.  Two consecutive certified-warm replays after splitting superseded work from stale source failures pass in 63/62 seconds at `/tmp/qged-lucy-superseded-fix-{a,b}-20260826`, each with a terminal ready, box-free mesh and zero occurrence failures.  The typed submission-owner replay at `/tmp/qged-submission-owner-lucy-20260826` passes in 82 seconds with 1.95M presented faces and 0.632-pixel maximum certified error.  After the exact source-delta owner replaced its independent active/source/plan fields, `/tmp/qged-submission-delta-lucy-20260826` also passes in 94 seconds with a box-free 1.265-pixel terminal view.  These lifecycles exercise continuous zoom/refinement, zoom-out recovery, rotation, lighting, hierarchy selection, subpath erase/redraw, camera ownership, and the strict HUD contract.  The 2026-08-28 exact-current warm wire lifecycle at `/tmp/qged-lucy-wire-capacity-transfer-revision-20260828` passes in 53 seconds after older frame-barrier precedence and capacity-transfer revision ownership were fixed; it completes all 191 events with 861,234 final mesh faces, zero boxes, and a clean external control trace.  Two canonical-handoff replays at `/tmp/qged-lucy-wire-canonical-handoff-{a,b}-20260828.json` first proved request-owned reconciliation budget and terminal-measurement retirement.  A later timing trace showed that ordinary planning could still ignore the retained terminal certificate and reselect the adjacent known-slow PoP population.  The terminal budget is now authoritative for ordinary admission, and two independent full replays at `/tmp/qged-lucy-wire-terminal-planner-clamp-{a,b}-20260828.json` pass all 191 events in 51.4/52.6 seconds.  Events 33 and 147 terminate ready instead of repeating the search; both final frames are box-free with no owner, obligation, pending calibration, or runtime-contract violation.  The final-binary true-cold/certified-warm wire pair at `/tmp/qged-production-lucy-wire-osmesa-20260828` passes in 51/58 seconds.  Both 191-event replays finish box-free, ownerless, and quality-floor-clean; cold presents 559,494 faces within the measured software budget and warm safely admits 1,284,926 faces. | cold/warm shaded and wire OSMesa interaction rows green |
| 50k scale | Cold shaded scale checks pass on System GL and OSMesa in 8.7/10.0 s to terminal CONSTRAINED.  A later exact-current run exposed interactive deadline recovery walking down one PoP ordinal per missed frame even when the measured cost ratio already proved that insufficient.  Recovery now operates in render-cost space during motion as well as at rest; the prior-pose deadline floor remains quiet-only.  The warm OSMesa shaded matrix at `/tmp/qged-50k-superseded-fix-20260826` passes in 78 seconds: held-motion recovery reaches its responsive ceiling directly, and the terminal view has 7,566 mesh payloads, 1.52M faces, zero boxes/failures, and no control owner.  After replacing the raw source-plan/census latches with typed owner-thread values, `/tmp/qged-submission-owner-50k-20260826` passes in 72 seconds and terminates ready with 6,633 progressive payloads, 1.32M faces, zero boxes/failures, and no pending work.  The subsequent exact source-delta ownership replay at `/tmp/qged-submission-delta-50k-20260826` passes in 95 seconds with 7,448 progressive payloads, 1.50M faces, 3.000-pixel certified error, and the same box-free terminal contract.  The exact-current warm OSMesa wire lifecycle remains green in 74 s.  The post-terminal-authority endpoint at `/tmp/qged-50k-terminal-planner-clamp-20260828.json` settles terminal and box-free in 7.0 seconds with 255 meshes, 49,745 subpixel occurrences, no owner/obligation/violation, and matching requested, active, and certified budgets.  Independent final-binary true-cold shaded OSMesa runs at `/tmp/qged-50k-constraint-witness-{final,repeat}-20260828` pass the complete interaction/camera matrix in 48.8/46.0 seconds.  Both select the same 6,935 mesh payloads and 43,065 subpixel occurrences, present 1.57--1.63M faces, have zero quality-floor misses/control violations, and explicitly record stable-budget, subpixel-aggregation, and static-deadline witnesses (`constraintEvidenceMask=13`) rather than relying on an inferred `performanceLimited` label.  The 2026-08-29 final-binary warm shaded OSMesa replay at `/tmp/qged-production-50k-osmesa-shaded-idempotent-20260829` passes the complete scripted interaction in 20 seconds after point-terminal census idempotence was enforced.  It ends box-free with 753 meshes, 49,282 aggregated occurrences, 363,752 faces, and a 126 ms software frame. | warm shaded/wire and repeated true-cold shaded OSMesa correctness, liveness, and interaction gates green; remaining System-GL cold/wire rows remain |
| 150k scale | The bounded cold crash gate passed under a 16 GiB address-space cap on System GL and OSMesa with all 150,001 leaves discovered and no terminal failures.  With the default scheduler (eight active realization workers on this host), the exact-current System-GL replay at `/tmp/qged_150k_timing_fix_report.json` reaches a stable terminal endpoint in 40.4 seconds with 78,341 mesh occurrences, zero visible boxes/failures, and no owner, obligation, background work, or runtime-contract violation.  The rebuilt typed-host OSMesa replay at `/tmp/qged_150k_typed_timing_osmesa_report_retry.json` terminates in 36.8 seconds with 4,097 meshes, 384,755 faces, 145,903 structural occurrences represented by 1,426 point draw records, and the same zero-box/control contract.  It is responsiveness constrained but has zero total/prominent quality-floor violations and zero visual-importance debt.  Threshold-stamped CAD timing prevents the prior 32/64-pixel balancing loop and stops faceplate-only frames from inflating CAD capacity.  The post-terminal-authority endpoint at `/tmp/qged-150k-terminal-planner-clamp-20260828.json` settles terminal and box-free in 29.3 seconds with 686 meshes, 149,344 subpixel occurrences, no owner/obligation/violation, and matching requested, active, and certified budgets.  A later trace found an unchanged structural-distribution seed reopening the already constrained point cut after its finer preload was rejected.  The idempotent reducer fix and formal contract are qualified by `/tmp/qged-production-150k-idempotent-seed-20260829`: the final-binary warm System-GL wire interaction passes in 27 seconds, remains box-free, and each constrained finer candidate is consumed once with no identical terminal replay.  This is liveness/performance evidence rather than final realistic visual-significance clearance. | current shaded System-GL/OSMesa liveness green; true-cold/wire rows and realistic prominent-object quality remain |
| System GL smoke | The current shaded Generic Twin cold/warm pair at `/tmp/qged-optional-policy-signature-generic-20260826` passes in 12/11 seconds after camera identity became an optional snapshot.  Both exact-view returns recall history and every terminal checkpoint is ready with zero structural boxes and occurrence failures; the final certified errors are 0.236/0.241 pixels.  The 2026-08-27 post-ownership full System GL interaction replay at `/tmp/qged-generic-growth-handoff-20260827/formal-contract-system-report.json` passed in 9.7 s: all 12 waits were terminal/ready, every wait was box-free with at least 673 meshes, and the final frame had 709 meshes and 135k faces.  The post-terminal-authority replay at `/tmp/qged-generic-system-terminal-planner-clamp-20260828.json` passes in 9.5 seconds and ends ready with 709 meshes, 135,073 faces, 57,102 lines, zero boxes, and no owner/obligation/violation.  The prior wire evidence is `/tmp/qged-generic-wire-system-debug-20260824`; both terminal wire frames contain 709 mesh payloads, zero boxes/failures/pending work, and matching cold/warm camera state.  The exact-current Lucy shaded cold/warm interaction replay at `/tmp/qged-lucy-system-20260824` also passed: both final frames are exact/ready, box-free, and quality-floor-clean with one resident source payload, 7.41M displayed PoP faces, a 230 MB resident mesh set, and one quality-history recall.  The final-binary Lucy true-cold/certified-warm wire pair at `/tmp/qged-production-lucy-wire-system-20260828` passes in 39/16 seconds.  Both runs finish with the same 3,183,110-face cut, 9,549,330 wire lines, no boxes, pending work, control violation, or quality-floor miss, and matching camera state.  Hubble shaded/wire cold/warm pass after the overview lifecycle repair. | Lucy shaded/wire green; complete the remaining System-GL large-model matrix |
| OSMesa Generic Twin | The current shaded cold/warm pair at `/tmp/qged-optional-policy-signature-generic-20260826` passes in 15/14 seconds.  Both exact-view returns recall history and terminate ready with zero boxes/failures; the final certified errors are 0.621/0.676 pixels.  The 2026-08-27 post-ownership full replay at `/tmp/qged-generic-growth-handoff-20260827/formal-contract-osmesa-report.json` passed in 14.3 s: all 12 waits were terminal/ready, every wait was box-free with at least 673 meshes, and the final frame had 709 meshes and 135k faces.  The final scene content and camera match the paired System-GL images.  The 2026-08-26 cross-renderer wire replay at `/tmp/qged-generic-wire-hud-final-20260826` passes cold and warm on System GL and OSMesa after the availability/HUD changes, with no invalid ready-label sample or terminal box.  The exact 2026-08-29 cold shaded cross-renderer replay at `/tmp/qged-production-generic-shaded-idempotent-20260829` passes in 11 seconds on System GL and 14 seconds on OSMesa; both end with all 709 meshes, approximately 135k faces, a one-pixel point threshold, and zero structural fallback.  The 2026-08-30 isolated-cold and same-cache hidden-line runs at `/tmp/qged-generic-osmesa-m4-{cold-r12,warm-r13}` pass in 11/13 seconds after cost-matched structural timing and exact-presentation recertification.  Every checkpoint contains all expected meshes with zero boxes/proxies/control violations; libicv SSIM is 0.998978--0.999833 with 0.0997%--0.2804% silhouette disagreement.  The differing cold/warm renderer costs are expected flat-batch versus retained-VBO submitted-work currencies and are consumed only with their paired durations. | compact direct path green across shaded, wire, and hidden-line cold/warm; continue larger real-model matrix |
| spatial Lucy, System GL | Exact-current 2026-08-25 warm shaded retained-page replay passes.  It refines during continuous zoom, reaches its quiet pixel target, compacts on zoom-out, restores the prior view through quality history, and has no boxes or page-level subpixel proxies. | focused retained-page regression green; cold and wire rows remain |
| direct-mesh Generic Twin | The 2026-08-29 cold replays at `/tmp/qged-generic-promotion-{cold,osmesa}` pass on System GL and OSMesa.  Late safe-scene certification now promotes an already visible PoP payload atomically by marginal replacement cost.  Each terminal frame has all 709 BoT occurrences in direct full detail, 185,388 faces, zero progressive/proxy/structural-box payloads, and no pending work.  The overall extent preview remains a discovery-only presentation and is not counted as a semantic leaf. | focused admission regression green; keep discovery-preview latency distinct from terminal-mesh admission |
| matched terminal image quality | The matched full-detail/LoD shaded matrix uses the same qged binary, per-orientation camera, physical canvas, lighting, and renderer, with libicv SSIM/PHASH plus exact and one-pixel-tolerant silhouette disagreement.  Generic Twin is effectively identical to full detail; Hubble, certified-warm Lucy, and NIST BREP satisfy their pixel-target or typed-constrained contracts without missing topology.  The current 4.96M-face heterogeneous 5k control compares at 0.978799--0.990360 SSIM on System GL while managed LoD presents 1.56M--2.26M faces; its tolerant silhouette disagreement is 0.004%--0.031%.  OSMesa presents 988k--1.66M faces at 0.943401--0.981091 SSIM under an explicit 169--333 ms responsiveness constraint.  The partial real-world Big Boy BoT conversion retains 14%--23% of 6.93M instance-expanded faces at 0.968001--0.982677 SSIM on System GL.  Its corrected OSMesa path reaches 0.946401--0.973746 SSIM without boxes, terminal proxies, or protected-floor debt.  A cross-renderer comparison found that non-exact PoP cuts disabled culling but fixed-function OSMesa still used one-sided lighting.  Culling and face orientation now follow the displayed-cut contract, while missing normals consistently preserve the source surface as a stable LoD appearance attribute.  The split is covered by Obol CPU, shader-source, fixed/GLSL two-sided, and real Lucy regressions.  Named wheel, running-gear, and boiler/cab crops show that the remaining high exact side-view disagreement is predominantly one-pixel boundary movement rather than a missing prominent component.  The 2026-08-31 wire and shaded-with-edges comparison adds Generic Twin and Hubble on both renderers.  Generic Twin reaches the complete mesh at 0.998978--0.999986 SSIM.  Hubble wireframe has below 0.008% one-pixel-tolerant silhouette disagreement; shaded-with-edges remains below 0.65%, including explicitly responsiveness-limited OSMesa cuts.  Exact metrics, content-digest controls, safe large-control rules, and provisional targets are in `obol_lod_visual_quality.md`. | simple/moderate shaded, wire, and shaded-with-edges; heterogeneous 5k; and partial real-world Big Boy cross-renderer baselines green; full BREP train, multi-large-mesh image comparison, and independent production-vehicle qualification remain |
| 50k/150k visibility mutation | The final warm 50k OSMesa workflow at `/tmp/qged-selection-capacity-50k-osmesa-final-20260831` selects in 0.93 s, erases 16 exact occurrences in 1.36 s, restores them in 1.40 s, and deselects in 0.90 s.  Visibility changes 50,000/49,984/50,000 while the terminal retained scene keeps 1,206,806 faces, zero boxes, and the same 5,130,680-unit certified budget; the complete matrix and camera contract pass.  A later timing-sensitive replay retired three mesh payloads while the subpath was hidden and exposed a missing prerequisite: source visibility became current before the successor framebuffer had classified restored occurrences.  The controller now enforces exact census, then exact presentation/classification, then reallocation.  `/tmp/qged-visibility-prerequisite-50k-osmesa-20260901` restores all 50,000 occurrences in 0.83 s and terminates with 1,343 meshes, 48,657 subpixel occurrences, zero boxes, and no owner or obligation.  The 150k OSMesa qualification at `/tmp/qged-selection-capacity-150k-osmesa-fixed-20260831` selects in 0.65 s, erases in 6.05 s, restores in 3.96 s, and deselects in 0.63 s.  It retains 835,318 faces, zero boxes, and the same 1,446,651-unit budget within one rounding unit across 150,000/149,984/150,000 visible occurrences.  All terminate ownerless with monotonic visibility revisions and pass the six-domain runtime trace.  Exact visibility deltas avoid a full inventory rescan and preserve renderer-capacity evidence. | shaded OSMesa hierarchy selection and exact mutation green at 50k/150k; wire and nested/deep hierarchy rows remain |
| ordinary-model terminal latency | The final Generic Twin cold/warm interaction repeats reach the first terminal mesh presentation in 0.71--0.93 seconds on System GL and 1.27--1.55 seconds on OSMesa.  Each endpoint has 673 active CAD payloads, zero structural boxes, no background/progressive obligation, and a ready convergence witness.  A later intermittent OSMesa run reached 99% after frame timing opened capacity/handoff work without arming another Qt timer.  Frame completion now re-evaluates the level-triggered host-work predicate; `ObolHostWork` checks that transition, and three consecutive untraced OSMesa quality runs terminate in 8--9 seconds with 134k--135k faces and no surviving qged process.  A separate Lucy race showed provider teardown clearing this shared pump after its final mutation left exact-frame debt; teardown now synchronizes every independent obligation, the executable regression passes, and three repeated System-GL Lucy lifecycles settle identically in 4.42--4.50 seconds.  The gate measures draw-to-idle rather than a fixed-delay screenshot.  Its nested-loop presenter forces the real QOpenGLWidget paint path, and a one-command test batch no longer disables its only canvas update.  Structural timing evidence is now an inseparable `(presented cost, duration)` pair, so a cheap box-only frame cannot authorize a richer current scene.  Compact direct scenes use a 400 ms prominent replacement deadline while large censuses retain their 200 ms finite deadline, and an exact presentation above a stale allocation certificate re-enters capacity reconciliation instead of becoming an ownerless 672/673 endpoint. | dual-backend cold/warm latency, exact presentation controls, cost-matched structural admission, and completed-frame wakeup green |
| concurrent cold cache | On 2026-08-28 two System-GL qged processes opened the same `Generic_Twin.g` and one initially empty `BU_DIR_CACHE` simultaneously.  Both exited successfully with terminal ready, box-free mesh presentations and no occurrence failure.  `mdb_stat` validated all three resulting LMDB environments (117 LoD, 83 draw, and 4 name-map entries), and an independent third process reopened the shared result, loaded from cache, and again converged terminal/ready without a cache write. | cross-process cache correctness green; large-asset duplicate-generation resource pressure remains to qualify |
| Hubble OSMesa | Exact-current 2026-08-25 shaded cold/warm lifecycles and camera contracts pass after the static-quality restoration change, with 1,804 warm shaded payloads and zero structural boxes.  The final-binary cold/warm wire pair at `/tmp/qged-production-hubble-wire-osmesa-20260828` passes in 16/16 seconds.  Both finish terminal, box-free, ownerless, and quality-floor-clean with 1,965 mesh payloads plus 553 subpixel occurrences.  The hierarchy replay retains one selected CAD instance across exact subpath erase/redraw, loads 2,530 tree items without a fallback scan, and applies its path update in at most 436 microseconds.  `/tmp/qged-resize-hubble-body-selection-final-20260829` adds four shaded managed/full-detail resize runs on System GL and OSMesa using the substantial `all.g/BODY` region rather than the old small probe.  Selection survives every window transition and the resize storm; the terminal selected frame is box-free and idle, clearing it changes 4,763--4,925 scene pixels, and all four final cameras match exactly.  The visually inspected frames show the same body panels receiving selected styling on both renderers.  The 2026-08-30 selected-subpath erase replay exposed allocator output feeding back through its external-cost input.  `/tmp/qged-hubble-selection-fixed-report.json` now settles after three bounded OSMesa allocation steps, and `/tmp/qged-hubble-selection-fixed-system-report.json` is visually stable from the first 250 ms checkpoint; both end terminal with no owner, obligation, or render request. | shaded selection/resize and shaded/wire cold/warm focused gates green; complete wire resize and broader physical-pointer modes |
| NIST BREP | The final-binary cold/warm wire lifecycle at `/tmp/qged-nist-resident-wire-final-20260828` passes on System GL and OSMesa.  All four runs select cut 10 with 19,457 lines, refine to cut 11 with 28,892 lines on zoom-in, and return to cut 10/19,457 lines on zoom-out.  Each endpoint is terminal with one direct progressive occurrence, zero structural boxes, and zero LoD-service tasks; cold HUD state remains `Discovering model` until that occurrence arrives.  The focused source-mutation test replaces the resident progressive part through the sparse delta journal and proves its cut and aggregate metrics retire.  The paired shaded OSMesa regression at `/tmp/qged-nist-shaded-regression-20260828` remains box-free and view-responsive.  The current resize/policy matrix at `/tmp/qged-resize-nist-contract-final-20260901` passes all eight System-GL/OSMesa, managed/off, shaded/edge workflows in 5--12 seconds; inspection confirms matching shaded silhouettes and the expected richer adaptive edge tessellation when LoD is enabled. | adaptive wire cold/warm and dual-backend shaded/edge resize gates green; shaded zoom-memory reclamation and larger real-BREP qualification remain |
| multi-Lucy/xpush capacity | Exact-current 2026-08-26 warm initial-view runs pass on both renderers after capacity samples were ordered behind ceiling-free occurrence-plan handoff.  The pre-fix timing replay split 2/4 between pixel demand and a premature ordinary-deadline endpoint; after the applied-but-hidden candidate fix, 6/6 early-checkpoint OSMesa runs reached the identical approximately 984.1k-face, 0.457-pixel pixel-demand endpoint in 4.59--7.04 s with no boxes or control work.  Post-recovery-ownership OSMesa replays remain green: multi-Lucy reaches 985,820 faces before a performance-limited zoom endpoint, while xpush reaches 985,808 faces initially and 1.53M after zoom, all box-free and terminal.  System GL reached 3.81M faces at 0.228 pixels in 2.45 s.  The 2026-08-27 full 75-event OSMesa turnover/zoom/rotation replay at `/tmp/qged-multi-lucy-growth-handoff-20260827/formal-contract-report.json` completed in 29.2 s after completed-pass ownership and annotation lifetime were made atomic.  All waits settled terminal/ready with zero failures or control violations; the final view had eight meshes, 319k faces, zero boxes, and no remaining owner or obligation. | warm full-interaction regression green; cold eight-distinct-asset preparation, compaction, and wire rows remain |

The 2026-09-01 post-visibility-contract image qualification supersedes the
older shaded point measurements in the table without weakening the remaining
release rows.  Fresh byte-matched full-detail/LoD comparisons cover Generic
Twin, Hubble, Lucy, heterogeneous 5k, partial Big Boy BoT, and NIST BREP on
System GL and OSMesa; all 48 checkpoints satisfy the exact presentation,
geometric error, topology, proxy, prominent-floor, and constraint-evidence
contracts.  Managed-only multi-Lucy, 50k, and 150k tiers also terminate
box-free and proxy-free with every classified prominent occurrence meshed.
The current 150k System-GL/OSMesa endpoints take 218/84 seconds, present
20.49M/1.12M faces, and peak at 9.27/7.52 GB process RSS while retaining
688/79 MB of controller-accounted mesh data.  Complete reports and realistic
expectations are in `obol_lod_visual_quality.md`; artifacts are under
`/tmp/qged-lod-quality-post-visibility-*-20260901`.

The subsequent terminal-handoff audit found two concrete refinement defects,
not failures of the abstract convergence rule.  A stronger selective pass
could retire the local cursor while leaving scene-wide quality debt, and a
point-threshold change could create a new allocation without advancing its
capacity revision.  The demand obligation is now level-triggered, all
classifier setters publish their semantic revision, and the runtime producer
validator recognizes complete inventory/submission coverage without accepting
an arbitrary selective submission.  The matched Hubble rerun at
`/tmp/qged-lod-quality-hubble-contract-final-20260901` passes modes 0 and 4 on
System GL in 8/9 seconds and OSMesa in 18/27 seconds.  Every endpoint is
terminal and ready with all 2,030 payloads satisfied, zero terminal proxies,
and a clean external control trace; the OSMesa hidden-line return is explicitly
responsiveness constrained rather than silently declared pixel-exact.
Focused point-selection replays at
`/tmp/qged-hubble-point-selection-contract-current-{system,osmesa}-20260901.json`
also pass the external trace.  Selection preserves the view and policy
revisions and the same 1,782/1,782 satisfied population; both renderers are
already idle and ready with zero boxes, proxies, obligations, service work, or
control violations at the first GUI readback and remain so after two seconds.

The 2026-09-01 capacity/HUD audit fixed two further refinement defects before
changing presentation text.  The protected visual floor is now the capacity
search's exact lower bound and capacity-owned allowances remain frozen for the
search key; this removes the out-of-domain candidate/floor alternation.  Exact
frame debt is distinct from exclusive interrupted replay, so a capacity
candidate allocates before its downstream frame.  Finally, a selective source
scope now retires with its completed cursor instead of blocking the broader
result-demand successor.  The final warm 5k OSMesa gate at
`/tmp/qged-lod-hud-unique5k-osmesa-final-20260901` passes in 32 seconds with
987,826 faces, zero boxes/proxies, a terminal ownerless endpoint, and a clean
325-transition external trace.  Its bounded capacity rank reaches 40/40.

The subsequent Hubble OSMesa mode-4 audit found two concrete-to-contract
mapping defects.  A complete pixel-demand endpoint was being assigned a new
plan serial when only its now-irrelevant protected-floor trial changed, and a
missing coverage pass could clear the demand-refresh producer of a still-live
importance census.  Allocation now canonicalizes that complete endpoint, and
only a completed dense ordinary pass may retire physical-demand refresh.  The
exact reproducer at
`/tmp/qged-contract-fixed2-hubble-osmesa-m4-20260901` reaches terminal/ready in
26.3 seconds with all 130 transitions clean.  The focused HUD replay records
capacity search as an exact local rank (for example, `search 13%`, `probe 2/8`)
alongside whole-view progress and names final allocation, publication, and
presentation owners instead of describing every last-percent handoff as
`Improving view 99%`.

The final wire/edge replay at
`/tmp/qged-lod-quality-wire-edge-contract-final-20260901` then passes Generic
Twin and Hubble modes 0 and 4 on System GL and OSMesa: all eight rows terminate
ready in 4--22 seconds with zero structural boxes, terminal proxies, prominent
floor violations, or runtime-contract violations.  All sixteen Generic Twin
views meet their strict safe-direct image targets.  Hubble's software
hidden-line endpoint is explicitly responsiveness constrained; paired-image
inspection retains every major component and the long cable rather than
hiding a topology loss behind aggregate image metrics.

The complete post-contract visual ladder is retained under
`/tmp/qged-lod-quality-post-contract-*-20260901`.  Fresh matched shaded runs
cover Generic Twin, Hubble, NIST BREP, Lucy, 5k distinct meshes, and the
partial Big Boy BoT on both renderers; all 48 views pass their image,
presentation, and external control-trace contracts.  Managed two-L Lucy,
50k, and 150k tiers terminate box/proxy-free with no prominent-floor miss.
The current 150k System-GL/OSMesa initial views finish in 76/67 seconds under a
12 GiB process limit, reach their exact 40/40 capacity rank, and retain all
current prominence candidates.  System GL reduces the preceding terminal
population from 20.49M to 4.51M faces while changing only 0.037 percent of the
one-pixel-tolerant silhouette; the OSMesa change is 0.574 percent.  This is
evidence that the screen-space allocator can remove substantial invisible
work without losing the whole-model envelope.  It does not replace the open
independent production-vehicle significance gate.

The quality harness now emits executable matched and managed assessments.
Safe-direct SSIM is a hard contract.  For larger shaded scenes, SSIM remains
an advisory photometric target while pixel error and one-pixel-tolerant
silhouette are independent geometry gates; coarse normals may change lighting
without moving a silhouette.  Managed-only telemetry can establish a valid
bounded endpoint but is always marked for visual review because it has no safe
whole-scene oracle.

The 2026-08-31 50k selection-stabilization replay exposed a missing render
classification.  Capacity relevance alone could not distinguish an LoD-owned
structural classifier from a selection/style presentation; both are
non-capacity frames when no measurable mesh population is present.  The frame
transaction now carries independent capacity and LoD-planning relevance.
Selection is excluded from both, while the classifier remains planning
relevant.  `ObolSemanticPresentation` proves that its isolated exact successor
can manufacture neither kind of evidence.  The strengthened graphical matrix
waits for selection readiness before erase.  Warm shaded 50k replays pass on
OSMesa and System GL under
`/tmp/qged-50k-selection-planning-{final,system}-20260831`.  The final OSMesa
binary retires the exact selection frame in 0.69 seconds to owner/obligation
zero without changing its 1,195,385-face population.  System GL returns in
0.19 seconds to the same
pre-selection background-compaction owner and unchanged 14,493,432-face
population.  Neither opens a capacity sample or headroom-refinement frame.

The 2026-08-31 selection follow-up also covers a cold shared-mesh publication
edge which the earlier single-object UI probe did not exercise.  A provisional
prepared preview was incorrectly treated as a resident progressive binding,
so one occurrence could keep projection refresh alive after the shared asset
was ready.  The capability-based classifier and direct shared-asset
replacement now pass the focused update-action test, a true-cold eight-instance
Lucy lifecycle, the dual-renderer selection UI matrix, and cold/warm Hubble
hierarchy selection/erase/redraw.  Every terminal observation has zero owner,
obligation, queued task, and in-flight task; all eight Lucy occurrences own
progressive payloads.  Reports are retained under
`/tmp/qged-selection-stabilize-*`; regenerated cache payloads were removed.

The subsequent full Hubble OSMesa hierarchy replay showed that the selection
frame was healthy but selected subpath erase/redraw could still expose a
capacity/triangle-recovery ownership violation.  A stronger `ALLOCATING`
capacity candidate had been represented by the generic preserve-budget request
and was clamped to the recovery ceiling; every exact frame then rejected that
different population as stale and restarted the search.  Capacity candidates
now have a distinct retained-allocation request kind which weaker recovery,
deadline, and reconciliation requests cannot replace.  The focused policy
test exercises the simultaneous recovery ceiling, and
`ObolPointQualityOwnership` passes 66 generated / 36 distinct states to depth
8.  Fresh exact-current warm Hubble workflows pass on System GL and OSMesa at
`/tmp/qged-hubble-post-selection-capacity-owner-final-{system,osmesa}-20260831`.
Selection returns to idle without changing its 620,449-face population, and
erase/redraw/clear-selection all terminate without an owner, obligation, or
capacity-sample loop.

The latest high-cardinality System-GL evidence supersedes the older terminal
populations in the cumulative scale rows above.  The bounded-batch 150k shaded
replay at `/tmp/qged-production-current-150k-batched-20260830` reaches a
whole-model preview during discovery, then terminates pixel-target and box-free
with 91,710 mesh payloads, 61,771 subpixel occurrences, 9.54M triangles, about
642 MB resident mesh data, and no prominent-floor or visual-importance debt.
Its first pose takes 179 seconds and the complete interaction workflow 229
seconds, so tens-of-thousands-of-distinct-mesh throughput remains an
optimization target even though the view is useful much earlier.  The focused
structural-frontier model now covers the bounded batches and their strictly
decreasing exact remainder over 60,535 generated / 29,266 distinct states to
depth 386.

The 2026-08-29 aggregate-proxy update keeps a genuinely tiny occurrence as one
point but preserves any footprint larger than five pixels as a batched box.
Cache-backed monolithic meshes use a validated PCA OBB when it is at least five
percent tighter by surface measure; absent, invalid, or insignificant metadata
uses the AABB.  Wire mode draws its 12 edges, shaded mode draws 12 lit
triangles, and shaded-with-edges draws both without creating per-occurrence
scene draw calls.
The same classification now applies to indirect-atlas admission pressure;
that path formerly emitted an unconditional point and discarded the source
bounds.  A forced one-MiB-atlas regression proves 128 screen-significant
pressure replacements produce 1,536 shaded triangles (and 1,536 edges when
requested), then collapse to exactly 128 points only after zooming below the
five-pixel ceiling.  Cache round-trip/validation and realization tests preserve
the OBB metadata, and the focused Obol CAD suite passes all 19 contracts on
System GL and OSMesa.  `/tmp/qged-obb-generic-20260829` additionally passes
the final-binary Generic Twin shaded cold/warm matrix on both renderers: all
four terminal frames have 709 direct full-detail meshes and zero AABB, OBB,
generic proxy, progressive, or structural-fallback payloads.  This proves the
OBB path does not defeat direct-to-mesh admission for a modest scene.

The heterogeneous 10k fixture at
`/tmp/obol-current-cache-matrix/unique-obb-10k` mixes mesh size, topology,
color, hierarchy, and baked non-axis-aligned local geometry.  A fresh two-pass
prewarm generated all 10,000 PoP hierarchies (8,349 then 1,651) and retained
1,960 materially tighter PCA boxes in both the asset and draw caches.  Under a
deliberately constrained OSMesa policy, the cold-manifest and warm-manifest
runs select the same 693 OBBs, 2,856 AABBs, 6,442 aggregate points, and nine
mesh payloads, with 29,571 presented faces, no structural fallback, and a
terminal ready state in 1.78/1.83 seconds.  The paired warm System-GL run
instead admits 6,962 mesh payloads and represents only 3,587 subpixel
occurrences as points; it terminates box-free with 920,108 presented faces.
The GUI matrix now derives its smooth-zoom target from database bounds, so
these images cannot pass by examining an empty fixed coordinate.  Repeated
explicit prewarm also skips valid current cache entries instead of resubmitting
the same prefix indefinitely.

The independently versioned draw metadata was then migrated over the preserved
current-format PoP caches instead of erasing 8.8 GiB of valid hierarchy data.
The final System-GL warm interaction replays pass at both pressure scales:
`/tmp/obol-current-cache-matrix/unique_mesh-50k-coverage-contract-warm`
finishes ownerless and box-free in 29 seconds with 34,746 active mesh payloads,
17,906 subpixel occurrences, and 4.59 million presented faces;
`/tmp/obol-current-cache-matrix/unique_mesh-150k-coverage-contract-warm-v2`
passes in 62 seconds with 13,937 resident mesh payloads, 97,087 aggregate
points, 8.86 million presented faces, and zero structural boxes or control
violations.  The 150k final framebuffer is already terminal and usable while
the HUD correctly reports its remaining bounded resident-prefix compaction as
background memory/cache work.  A policy-only revision no longer manufactures
a redundant coverage census; an active current-view census is now explicit
foreground inventory work and is reattached to a bounded cursor if a stronger
owner had paused it.

This clears persistence, cross-renderer semantics, and heterogeneous synthetic
pressure as implementation risks.  Real wheels/blades/booms/hulls still need
perceptual comparison against no-LoD ground truth.  Multi-page terminal proxy
ownership and optional enrichment of a manifest sealed before first-cold PoP
characterization also remain open; the active cold payload already carries
the OBB, but the initial structural cue deliberately remains an AABB.
The final Generic Twin resize/policy matrix at
`/tmp/qged-production-generic-resize-policy-disable-final-20260829` passes
auto and initially-off shaded cases on System GL in 6/10 seconds and OSMesa in
9/10 seconds.  Every endpoint has zero control owner, obligation, violation,
pending render/progressive work, and structural fallback boxes.  The
off/on/off checkpoints retain meshes rather than restarting at boxes, and
viewport/style repaints while disabled remain presentation-only.
The final point-bracket policy combines its timing-ratio correction with the
known safe/unsafe bounds, so it may skip redundant thresholds but can never
cross the proven coarse side.  `ObolCadTimingEvidence` exhaustively checks all
bounded jumps across the full 64-cut production domain (7,726,131 generated /
1,222,796 distinct states, depth 12), and the focused coordinator regression
covers the measured correction numerically.  The current 150k shaded OSMesa
replay at `/tmp/qged-production-150k-osmesa-threshold-acceleration-20260829`
passes in 49 seconds: its former twelve-frame `32 -> ... -> 64` staircase is
now `32 -> 55.3 -> 64`, and it terminates with 2,075 meshes, 147,964 aggregated
occurrences in 274 renderer records, zero structural fallbacks, and no control
owner or obligation.  Current follow-up replays pass for 50k shaded OSMesa
(`/tmp/qged-production-50k-osmesa-threshold-acceleration-20260829`, 20 s),
Lucy shaded OSMesa
(`/tmp/qged-production-lucy-osmesa-threshold-acceleration-20260829`, 45 s,
340,082 faces, 1.265-pixel maximum error, no quality-floor violation), and
cold shaded Generic Twin on both renderers
(`/tmp/qged-production-generic-shaded-threshold-acceleration-20260829`,
12/14 s).  These runs validate the complete client path in addition to the
renderer-level primitive/accounting checks.

The 2026-08-28 exact-current shaded smokes cover both renderers without
replacing the broader interaction rows above.  Generic Twin reaches a terminal
box-free endpoint in 1.85 seconds on System GL and 4.34 seconds on the rebuilt
OSMesa host,
with 673 meshes.  Lucy reaches a terminal box-free endpoint in 1.75 seconds on
System GL and 2.70 seconds on OSMesa; both present 559,494 faces at a
0.928-pixel maximum certified error.  The 50k fixture reaches a terminal
box-free endpoint in 11.1 seconds on System GL with 27,968 mesh occurrences,
23,006 point occurrences, and 2.27M faces.  OSMesa reaches a constrained
box-free endpoint in 14.1 seconds with 1,094 meshes, 1,431 point draw records,
and 378k faces.  All report zero runtime-contract violations.  The next scale
gates are true-cold large-asset latency and realistic visual-significance
qualification rather than another object-count-specific controller change.

The final-binary Generic Twin wire matrix at
`/tmp/qged-generic-wire-cross-renderer-fix-20260828` passes cold and warm on
System GL in 10/10 seconds and OSMesa in 15/19 seconds.  All four endpoints
contain 709 mesh occurrences, about 135k scene faces, zero structural
fallbacks, zero control work or violations, and exact camera-contract matches.
This row also regresses the completed selective source-delta to scene-wide
handoff: the former OSMesa warm loop repeatedly visited 14 satisfied entries
and could not certify its 709-occurrence allocation.

The 150k terminal view may validly be performance-limited.  It is not visually
qualified until prominent-object quality floors are demonstrated on realistic
large models, not merely synthetic coverage fixtures.

The 2026-09-01 BREP fallback and progress pass adds a current-tree 64/64
`drawing_baseline` run, a 28/28 `bobol_headless` run, and focused NIST PMI7-10
shaded and wire real-model passes.  The shaded real-model test finishes from
an isolated writable cache in about twelve seconds and the wire case in about
six seconds.  A true-cold OSMesa GUI replay exposes the exact four-item source
preparation rank, replaces each temporary box with its completed mesh, and
terminates with no boxes.  The matched managed/control image reaches SSIM
0.998831 and zero one-pixel-tolerant silhouette disagreement while using fewer
faces than the control.  Cache-context failure is now separately diagnosed
and preserves valid tessellation as terminal non-LoD geometry.  This closes
the misleading PoP-rejection and indefinite-box failure modes.

The subsequent constrained-residency lifecycle at
`/tmp/qged-lod-brep-pmi7-10-reclaim-trace15-20260901` (System GL) and
`/tmp/qged-lod-brep-pmi7-10-reclaim-osmesa-20260901` (OSMesa) qualifies the
NIST shaded BREP zoom/reclaim/return/restore path against matched LoD-off
controls.  Both backends present all four targets with zero structural boxes,
terminal proxies, or prominent-floor violations.  The close view and its
cache-backed restore both present 220,637 faces.  System GL records SSIM
0.993568 and OSMesa 0.994099--0.994100; every checkpoint has zero silhouette
disagreement after a one-pixel tolerance.  The deliberately reduced resident
limit produces an honest memory- or responsiveness-constrained terminal
status without losing the pixel-target presentation.  A changed resident
capacity epoch now permits one byte-governed physical-demand retry past its
stale allocation; previously that retry scheduled zero work and could publish
a false coarse terminal state.  NIST adaptive shaded BREP residency is green;
the larger original Big Boy BREP hierarchy remains an open release row.

The HUD now reports finite source preparation and capacity-search ranks rather
than parking at a generic `Improving view 99%`.  A warm OSMesa NIST probe
(`/tmp/qged-lod-progress-probe-report2.json`) showed the source rank advancing
0/4, 1/4, 2/4, and 3/4 before terminal readiness.  The fourth BREP dominates
the remaining wall time, so this is deliberately a finite work rank rather
than a fabricated ETA; the label continues to identify the active producer.

The 2026-09-03 HUD follow-up preserves those exact ranks and adds a separate
observed-rate estimate.  The estimate is available only when every unfinished
foreground rank has both a finite denominator and a learned rate; otherwise
the bar is explicitly indeterminate and the label retains the exact stage
counts.  Terminal and same-episode reopen unit tests prove that only a terminal
view may publish 100 percent.  Fresh cold System-GL and OSMesa Generic Twin
runs end ready with all 709 CAD payloads and zero structural boxes.  The final
OSMesa evidence is retained at `/tmp/qged-hud-estimator-final`.  A fresh cold
OSMesa Lucy lifecycle at `/tmp/qged-hud-estimator-lucy-final` ends ready with
2,101,208 shaded faces and zero boxes; no nonterminal sample reports 100
percent and its bounded 20-unit capacity search exposes a changing ETA.
All 36 `bobol_headless` tests pass after a complete relink.  The current full
GUI reports still fail the newly expanded transition-journal audit on unnamed
internal transitions and transient owner/revision observations even though
their framebuffer, interaction, and terminal predicates pass; resolving that
instrumentation-to-atomic-reducer mapping remains a qualification item and
must not be hidden by weakening the checker.

The final current-tree visual spot audit after the background-result handoff
repair is retained under `/tmp/qged-lod-quality-*-20260901` and summarized in
`obol_lod_visual_quality.md`.  Generic Twin and Hubble pass matched control/LoD
wire and shaded comparisons on both backends; NIST BREP and the partial Big Boy
BoT pass their matched System-GL rows; and single/multi-Lucy plus shaded
50k/150k managed endpoints are coherent, terminal, box-free, proxy-free, and
prominent-floor-clean.  A queued cache result can no longer demote an already
resolved semantic revision from terminal to foreground refinement.  The 50k
wire endpoints are correct but explicitly constrained at 105--145 ms per
frame and approximately 2--3 pixels maximum projected error.  Improving
high-cardinality wire submission/raster throughput without changing visible
edge semantics is therefore an optimization gate, not a prerequisite for the
shaded quality result.  The HUD reports `render primitives` from the latest
exact completed-frame record, so direct wire and full-detail channels are no
longer mislabeled or reported as zero.  It falls back to the active PoP face
estimate only while no exact frame record is available.

The September 3 tiled-lighting audit found valid spatial-page positions and
indices but inconsistent fallback-normal semantics between Obol renderer
tiers.  All tiers now retain the quantized position for drawing and lighting
location while deriving a missing face normal from the exact source triangle;
the displayed face is only a degenerate-source fallback.  Focused CPU, shader,
two-sided, and complete CAD renderer contracts pass, and short Lucy captures
on System GL and OSMesa show coherent folds without page-local
lighting changes.  Immutable spatial pages now honor all view-normal modes.
The cache retains a fixed-width source-vertex ID beside each page-local point;
workers use it to synthesize crease-aware smooth normals across the complete
admitted page set and canonicalize one bounded page at a time.  Flat mode keeps
the compact normal-free representation, and compaction cannot replace either
mode with an authored-normal variant.  The GUI thread performs no normal
expansion or large mesh copy.

The next BREP production-quality gate is the larger original Big Boy hierarchy:
complete its shaded source conversion, repeat the now-passing NIST constrained
close/zoom-out/cache-restore lifecycle, and use named train features for an
independent real-model significance review.  Cache-unavailable, invalid-input,
and PoP-construction failures still share a low-level integer return API;
replacing that ambiguity with typed diagnostics is cleanup debt, although valid
BREP tessellation is now preserved correctly in every case.

The September 1 single-wheel probe makes the source-conversion prerequisite
concrete.  The legacy CDT silently aggregated a partial result after five
constrained-edge face failures: only 32 of 368 allocated points were
referenced, and the resulting 128 indices covered one upper plane rather than
the tire's full thickness.  PoP classified that defective source faithfully;
it did not reject or repair it.  The indexed-face provider now fails closed
when referenced output omits the exact BREP boundary extent or escapes the
surface envelope, and shaded BREP cache version 3 invalidates the old partial
payload.  The enhanced `brepdraw` audit shades the same primitive with all
eight faces represented by 1,336 triangles in 0.009 seconds.  A full-train
comparison remains blocked until both display policies share that bounded,
validated source contract, including deadline, memory/result limits,
cancellation, per-face completion, and typed failure.

The guard is selective rather than a global BREP fallback.  All four NIST
PMI7-10 solids pass the provider diagnostic, and the fresh-cache OSMesa
comparison at `/tmp/qged-lod-nist-brep-v3-20260901` completes in 18/24 seconds
for control/LoD.  It presents 178,561 LoD faces with SSIM 0.980645, zero
one-pixel-tolerant silhouette disagreement, and no boxes, proxies,
quality-floor violations, or terminal error.

## Common acceptance criteria

Every graphical row must show:

- correct final camera, extent, draw mode, colors/materials, lighting, and
  requested normal style across ordinary and spatial-page meshes, and retained
  scene composition;
- exact visible occurrence classification, no unexpected boxes, no empty
  frame, no invalid geometry, and coherent shaded/wire cut semantics;
- no owner-thread pending-without-witness loop, stale result, coordinator
  invariant violation, or unbounded worker/resident allocation;
- semantic selection, tree state, highlights, erase/redraw, and edit
  promotion consistent with command and widget state;
- stable-frame, first-useful-image, interaction-latency, resident/peak-memory,
  cache, and representation telemetry within the row's declared threshold;
- visually inspected checkpoint images.  APNG is required for suspected
  flicker; apitrace is required for System-GL-only corruption or state leaks.

The GUI runner may allow cold asset construction time, but a timeout, crash,
sanitizer finding, empty capture, unexplained terminal box, pending idle loop,
or a missing report fails the row.

## Required models

| Model | Primary purpose |
|---|---|
| Generic Twin | production visual baseline; tail/hull/engine, autoview, lighting |
| Lucy | spatial-page demand, close zoom, smooth refinement, compaction |
| multi-Lucy/xpush | shared-asset reuse and visibility turnover |
| Hubble | deep hierarchy, selection/tree scale, small-component importance |
| Havoc | mixed hierarchy and real-model interaction |
| NIST BREP corpus | wire/shaded BREP and adaptive tessellation behavior |
| Stanford meshes | varied manifold/non-manifold mesh behavior |
| 5k/50k/150k varied fixtures | distinct assets, mixed sizes/colors/regions/hierarchy, aggregate budgets |
| independent multi-gigabyte vehicle | real visual significance, I/O, memory, and unique-mesh behavior |

Synthetic fixtures must include mostly manifold inputs, mixed mesh sizes,
visually prominent components, color/region variation, and deep hierarchy.
Repeated-instance and distinct-asset cases are separate obligations.

The explicit `test_bobol_database_draw_corpus` qualifier exercises every
named solid and combination, including nested roots, rather than only the
usual database tops.  Its normal view-managed pass currently has no contract
failure across the installed corpus; the deliberately repeated `geometry`
aggregate in `primitives.g` and large volumetric roots remain separately
time-bounded performance qualifications.  The exact parallel-axis ETO from
`photo_holder.g` has a dedicated finite-output regression, and bare Buddha
wire/shaded GUI checks must end with nonzero mesh data, zero terminal boxes,
and zero occurrence failures on both renderers.  The 2026-09-04 nested-root
pass covers all 105 `primitives.g` objects in wire and shaded mode with zero
failures (22 deliberate representation fallbacks and four valid empty roots).
Its tolerance-sampled VOL provider completes the formerly unbounded
`vol.r` rows in 0.1--0.3 seconds.  Separate Happy Buddha headless checks pass
both its bare BoT and region roots in both modes.  Current bare-BoT shaded and
wire GUI replays under `/tmp/qged-buddha-bare[-wire]-{system,osmesa}` terminate
ready on both renderers with one mesh, zero boxes, and zero occurrence
failures.  System GL retains a 952,476-face/2,857,428-line cut through motion.
OSMesa retains its resident 730,656-face/2,191,968-line cut, temporarily
submits a roughly 53,500-face visible working set during the measured shaded
rotation, and restores the same cut from view history; that transient
reduction is explicit interaction-budget behavior rather than source
rebuilding or eviction.

The final 2026-09-04 current-tree corpus rerun extends that evidence to every
installed database rather than only `primitives.g`.  With one isolated cache,
both shaded and wire passes covered 64 databases and 16,917 named solids,
combinations, and nested paths with zero failures.  Shaded reported 24
deliberate representation fallbacks, wire reported 12, and both classified the
same 43 semantically empty roots.  The focused mesh-cache, LoD-coordinator,
and update-action tests passed immediately before the corpus run.  This is a
provider/routing qualification; framebuffer composition and interactive
behavior remain the responsibility of the graphical rows above.

The same complete 33,834-realization corpus was repeated after the final
direct-leaf and presentation-liveness fixes with automatic LoD both enabled and
disabled.  Both configurations passed all 64 databases and 16,917 independent
roots with zero failures; automatic LoD reported 36 deliberate representation
fallbacks and the full-detail route reported 15, while both classified 86
mode-counted empty roots.  This closes the earlier LoD-off manifest-identity and
missing-bound observations rather than treating them as expected differences.

A final shaded Buddha motion replay on the 2026-09-04 binaries records the
renderer-dependent quality policy explicitly.  System GL presents all 952,476
resident faces before, during, and after motion, with measured scene time below
0.6 ms.  OSMesa presents the 597,464-face resident cut at rest in 131--165 ms,
temporarily caps the exact motion frame at cut 16 / 78,516 faces in 25.6 ms to
meet the 20 FPS interaction target, and restores the retained cut 445 ms after
release without provider rebuilding or boxes.  Both terminate ready with no
ceiling or constraint.  The same audit repaired a timing-dependent suffix-
publication seam: a ceiling hiding a richer retained population is now
foreground convergence debt even if compaction is the only other work.  The
System GL and OSMesa reports are respectively
`/tmp/qged-buddha-current-shaded-system/fixed-report.json` and
`/tmp/qged-buddha-current-shaded-osmesa/fixed-report.json`.

The matching bare-Buddha wire replay closes the second renderer-ceiling edge.
A static trial which was already probing when a richer cumulative suffix became
resident now requests one exact successor, and an inert ceiling equal to the
active maximum is removed through one ceiling-free commit frame.  System GL
terminates at the full 952,476-face / 2,857,428-line population with a 0.30 ms
scene measurement.  OSMesa terminates at cut 23, 730,656 faces, and 2,191,968
lines with a 190 ms stable measurement; its measured motion cut is 16 / 78,516
faces at 30 ms and the retained stable cut is restored without boxes.  Both
end ownerless, unconstrained, and with no convergence obligation.  Their reports
are `/tmp/qged-buddha-wire-system-current/fixed2-warm-report.json` and
`/tmp/qged-buddha-wire-osmesa-current/fixed2-warm-report.json`.

The hidden Qt progressive tests previously consumed a request after reporting
frame completion, which could swallow the successor frame published during
completion and strand a final structural box.  They now share the production
claim-before-traversal protocol.  Cold OSMesa Generic Twin settles in eight
feedback frames with 532 mesh payloads, 177 subpixel occurrences, and zero
boxes; both shaded/startup and wire CTests complete in under two seconds.

After enforcing the allocator-certificate prerequisite on spatial-page cut
downgrades, both warm Buddha wire replays again terminate with zero boxes,
owners, obligations, or violations.  System GL retains the full cut 28 /
952,476 faces through interaction and measures 0.22--0.69 ms.  OSMesa
terminates at cut 25 / 876,752 faces; its stable and post-motion images are
byte-identical, while the 32.8 ms interaction frame uses the cut-16 renderer
ceiling without mutating the retained cut.  The current reports are
`/tmp/qged-buddha-wire-system-current/certified-limit-report.json` and
`/tmp/qged-buddha-wire-osmesa-current/certified-limit-report.json`.

The corresponding drawing-baseline qualification completed with no failure on
2026-09-04: 63 tests passed and the Xvfb System-GL context probe reported its
defined environment skip.  That run includes qged polygon and primitive-edit UI traces,
selection/picking/snap/measurement behavior, faceplate and framebuffer
composition, MGED and rtwizard views, NIST BREP cases, and the real-model
progressive tests.  The two legacy libged image comparisons initially rejected
the retained renderer's deliberately richer database wireframe (SSIM 0.974
against a coarse display-manager control); the captures were complete and
coherent, so their structural thresholds now admit that background
tessellation difference while retaining the semantic geometry checks.

### 2026-09-04 preparation-cost and frame-delivery follow-up

The Buddha intermediate-quality reduction above is renderer-budget behavior,
not missing residency.  It is not evidence that every reduction is optimal.
The warm Lucy audit at `/tmp/qged-lucy-capacity-audit-20260904` exposed a
separate software issue: cut 27 rendered repeatedly in about 276 ms, but its
first restored frame took 691 ms and was classified as draw overload rather
than preparation.  Successive fallback first frames then rejected affordable
cuts.  Do not relax the close-zoom prominent-quality floor to hide this.

Obol's fixed-function shaded, explicit-wire, and indexed triangle-wire cut
uploads now publish finite, exact-target preparation certificates rather than
only activity serials.  Their occurrence-request frontier advances only after
a successful retained upload; a replay or re-upload behind the frontier cannot
create another progress witness.  The target includes view, context, plan,
geometry, ceiling, fraction, and projection identity.  The immutable work domain
does not depend on whether another executor handles one draw channel.  This
implements the existing preparation-progress contract; it does not increase
memory/FPS limits or prove bounded cancellation latency inside one large upload.
That latter implementation obligation remains explicit.

The new fixed-cut finer/coarser/replay regression passes.  The software CAD
suite passes 28/28; System GL passes 27/27 applicable cases, with the software-
specific regression deliberately skipped.  The rebuilt BRL-CAD coordinator and
update-action CTests also pass.  Logs are
`/tmp/obol-fixed-preparation-{software,system}.log`.

Full Lucy qualification is **not green**.  The desktop replay at
`/tmp/qged-lucy-fixed-preparation-20260904` stopped at event 28 after requesting
`render-preparation-replay`, with host flags 7, owner 6, no workers, and no later
render before the 60-second timeout.  A similar requested-frame stall preceded
this change in `/tmp/qged-lucy-capacity-trace-20260904`.  Diagnose canvas exposure,
Qt delivery, and endpoint claiming before attributing this to capacity policy.
Timeouts now include canvas visibility, update enablement, and visible-region
state.  Normal-style same-cut comparisons and the close-zoom quality floor also
remain mandatory; earlier warm runs had timing-dependent failures in these
checks.  Retain the certified shared Lucy cache at
`/tmp/qged-lucy-control-contract-cold-20260904/caches/lucy-osmesa-shaded-swapdefault/cache`;
reuse it rather than generating duplicate large caches for each replay.

The isolated-display follow-up at
`/tmp/qged-lucy-fixed-preparation-isolated-20260904` completes the full 217-event
script in 147 seconds without the desktop's missing-frame timeout.  This
narrows, but does not prove, the exposure/delivery hypothesis.  Its log confirms
the new fixed-function preparation kind reaches the controller and unchanged
replays stop advancing its certificate.  Qualification still fails: the initial
terminal cut is one boundary below the existing recognizable-prefix threshold,
the close view retains a prominent-floor violation despite 2.62M presented
faces, and smooth/authored normal transitions still require separate same-cut
analysis.  Do not report this as a passing Lucy matrix or assume the fixed-cut
reporting gap explains every earlier flat-atlas preparation/capacity failure.
Both graphical test processes exited; the shared cache was reused and no large
geometry duplicates were created.

### Normal-layer readiness follow-up

The 2026-09-04 replay at `/tmp/qged-lucy-normal-readiness-20260904`
exposed and verifies a narrower implementation correction: drawable old page
layers do not prove that a newly requested normal policy has been delivered.
The incremental presentable-payload count now includes immutable page normal
policy/crease identity.  An explicit policy change reclassifies the retained
payloads once; ordinary frame readiness remains O(1).  Wire, culled spatial
occurrences, and admitted terminal coverage proxies do not acquire fictitious
normal-layer debt.  The submit, cut-retarget, and readiness paths share the
same normal-policy comparison.  Existing meshes/assemblies stay resident until
the replacement arrives.  The coordinator/update-action tests pass, including
pending-worker/result, cancel/readback, publication, and crease-change checks.

The full 217-event OSMesa Lucy replay finished in 132 seconds.  The smooth
command's ready wait now takes about 16 seconds and returns with no worker
remaining; it previously returned while the required worker was still active.
This is **not** a passing full qualification.  Smooth cut 24 incurred repeated
400 ms non-preparation deadline interruptions and settled at cut 23 (1,284,926
faces versus the authored reference's 2,101,208).  The trace therefore gives
actual renderer-deadline evidence for this reduction, not just a suspected
normal-count cost-unit change.  Profile the renderer at controlled same-cut
normal styles before changing cost weights or relaxing the same-cut visual
check.  The smooth close view still presents 2,619,533 faces at cut 28 with one
prominent-quality-floor violation; the earlier discrete close also remains
too coarse.  Normal and close-view images were inspected and retained.

The first automated validation failure in this run was a separate unnamed
trace edge at observation 1135/event 141, before the normal commands.  The
reported phase changes from IDLE to BACKGROUND without changing the owner,
obligations, revisions, or terminal outcome in the exported trace.  This needed
a concrete producer witness, not suppressed unnamed-edge validation.
The subsequent admission-trace correction below resolves this missing edge;
it does not make the whole Lucy replay green.
All four Generic Twin shaded cold/warm replays (System GL and OSMesa) pass,
including their camera contracts, at
`/tmp/qged-generic-normal-readiness-20260904`.  Both suites ran sequentially
and exited; the Lucy run reused the single shared cache.

### Admission publication trace and software profile

The admission epoch which a worker publishes before its callback reaches the
controller was already included in convergence readiness, but absent from its
trace snapshot.  A reclamation could therefore change IDLE to BACKGROUND with
no exported producer fact to explain the edge.  The status snapshot now records
the service's published epoch and the controller's observed epoch; readiness
uses those same captured values.  Trace equality, producer-edge classification,
and complete-endpoint continuity use them too.  A deterministic no-worker test
raises the service limit without delivering a controller callback: it failed
before the fix, passes afterwards, and verifies that unchanged observations do
not generate repeated events.  Checker fixtures accept the corresponding
producer edge and reject a discontinuous admission history.  Unnamed-edge and
no-progress-cycle checks remain strict.

The 199 Hz OSMesa perf replay at `/tmp/qged-lucy-normal-perf-20260904`
completed all 217 GUI events in 95 seconds with no unnamed transitions; the
updated complete-trace checker passes.  It records the actual zoom-out edge at
event 139: service admission 1 -> 2, observed epoch still 1, named
`producer_progress`.  All four normal checkpoints in this run keep cut 23 and
1,284,926 faces (roughly 250-263 ms).  This is useful same-cut timing evidence,
not a passing Lucy quality result: the composite recognizability gate still
fails, and close-view qualification remains open.  The September 5 audit below
identifies representation-specific assertions in that composite gate; its
failure text alone did not establish a coarse initial prefix.  Profile
mode omits expensive deep diagnostics, so its overall elapsed time must not
be presented as an improvement over the previous fully instrumented replay.
The coordinator, update-action, and trace-checker CTests pass, as do the four
Generic Twin shaded cold/warm System GL/OSMesa GUI runs at
`/tmp/qged-generic-admission-trace-20260904`.  The superseded Generic Twin
normal-readiness run retains its evidence but no longer its duplicate caches.

The recording has 19,385 samples and no lost samples (1.84 MB).  Self samples:
OSMesa 65.64%, libBObol 13.70%, Obol 12.16%.  The largest symbols are OSMesa
`light_fast_rgba_twoside` (21.34%), `smooth_rgba_z32_triangle` (14.34%), and
`triangle_offset_twoside_rgba` (8.97%).  Obol flat-shaded packing is about
5.57%, and progressive cut preparation about 2.69%.  Kernel symbols were not
available; do not assign the unidentified samples to a specific kernel cause.
BRL-CAD is a Debug build, whereas Obol and its bundled OSMesa are optimized.

### September 5 resumed build and observation qualification

`/tmp/qged-obol-resume-20260905` retains this resume's reports, images,
reproductions, and `provenance` directory.  Obol/OSMesa were rebuilt, explicitly
installed into BRL-CAD's `.build`, and BRL-CAD consumers rebuilt.  Installed and
build-tree SHA-256 identities match: Obol `dd3ac7b1...`, OSMesa `097aac27...`.
The retained `qged-ldd.txt` resolves both through `.build/lib`.  These are
incremental rebuilds, not a clean-build or platform qualification.

The observation repair serializes uncollected compact-registry diagnostics as
JSON null, retaining the existing `deep_lod_diagnostics` availability flag.
Measured registries still require valid counts and entry-backed payloads.
`lod_lucy_contract.jq` independently names every failed coverage, presentation,
and prominent-floor condition before the other gates can return early.  It
also diagnoses historical reports whose disabled counts were serialized as
zeros; those reports and their original validation logs are preserved.
The old startup proxy was at 43 ms, before autoview; actual mesh coverage was
present by 1.715 s.  The first-useful-image and cold/warm mesh deadlines remain
unchanged.  Exact frames after observed mesh coverage, or claiming readiness,
may not regress to points.  Initial coverage publication can precede System GL
frame certification; the first-mesh checkpoint separately requires exact
rendered geometry.

Both fresh warm Lucy runs disable deep diagnostics and use the one retained
Lucy cache, one GUI process at a time on a private Xvfb display:

| Backend | Replay/control result | Remaining image/quality evidence |
|---|---|---|
| OSMesa | All 217 events finish in 123.478 s; all 2,505 control transitions pass. Base report and diagnostics-availability checks pass. | Discrete close: 680,242 faces, normalized error 3.705. Smooth close: 2,619,533 faces, error 3.581. Both miss the prominent floor despite a performance constraint. Zoom-out and final view meet the floor in this run; final presentation is 2,101,208 faces. |
| System GL | All 217 events finish in 105.435 s; base report and diagnostics-availability checks pass. The 1,209-transition trace fails the A/B/A checker at indices 714--716, event 139. | Smooth close: 510,112 faces, normalized error 8.709, with a prominent-floor miss. Final presentation is 2,101,208 faces and floor-clean. |

**Both Lucy qualifications remain FAIL.** The System GL control failure
surrounds `lod-static-overscan-resident-growth`; request/completion returns to
the same semantic tuple without a recorded finite-work decrease.  A following
planning observation advances the presented-frame serial.  Serial growth
alone is intentionally insufficient and the checker was not relaxed.
`provenance/system-control-cycle.json` retains the surrounding transitions.
This completed isolated-display replay does not resolve the earlier desktop
frame-delivery timeout.  The software close endpoints also retain
`lod_allocation_protected_floor_allowed=false` and zero maximum protected
budget; their numeric admission/rejection boundary needs investigation.

Software flat/authored/smooth/return checkpoints hold cut 24 and 2,101,208
faces at approximately 326/317/369/308 ms.  System GL uses cut 22 and 861,234
faces for flat/authored, drops to cut 21 and 559,494 faces for smooth, and
returns to cut 22.  Its timings are therefore not a same-cut normal comparison.
Neither run establishes a speedup, cross-page continuity, or hardware-GL
qualification merely by selecting a rendering backend under Xvfb.

All four Generic Twin wire cold/warm System GL/OSMesa matrix rows pass, including
camera contracts, in approximately 12 seconds each.  Focused OSMesa probes
of Lucy's combination root and bare Buddha each reach a ready, exact, nonzero
mesh with diagnostics on/off.  Measured final registries contain one entry and
one payload; disabled fields are null.  Buddha's retained historical cache
path was absent, so the authoritative repeats use the existing private
`buddha-cache` directory and `buddha-bare-private-deep{0,1}` reports; earlier
`buddha-bare-deep{0,1}` probes do not establish cache conditions.

A distinct bare-giant failure is retained at `bare-deep0/report.json`:
`draw -m1 lucy.s` is source-working-set constrained and the mesh-ready gate
expires at 60 seconds with no CAD payload.  The source coordinator estimates
serialized-leaf import bytes before manifest probing; combinations reserve
the finite coordinator capacity instead.  Do not fix this by disabling the
pre-import bound.  A bounded cache/metadata admission path is still needed.

The five focused CTests pass in 2.45 s.  The new adversarial Lucy validator
test also passes after the System GL coverage/publication distinction was
added.  No TLA model or controller transition implementation changed.

### September 5 follow-up: replay ownership and bare-root source admission

The retained A/B/A occurred while a submission pass was still active.  The
static-quality requester queued a frame which its completion consumer was
required to ignore until submission/capacity/publication finished.  Both now
use the same readiness predicate.  The focused executable regression exercises
submission, presentation-barrier and publication ownership, then confirms the
trial can request its successor after those producers complete.  Existing
terminal-composition phases already order allocation/presentation before static
quality; no new model state or generic frame-serial progress witness was added.

`/tmp/qged-obol-static-replay-20260905` retains the follow-up System GL warm
matrix.  All 217 events complete in 131.530 s and all 1,184 control transitions
pass the unchanged checker.  **Lucy quality still fails:** smooth close has
1,110,715 faces and normalized error 5.717 with one prominent-floor miss.
Flat/authored/smooth/return now hold cut 24 and 2,101,208 faces, with observed
scene times of 289/262/270/277 ms.  These are same-cut observations, not a
controlled end-to-end speedup result.

Bare streamed v5 BoTs now reserve the finite source allowance for the existing
serialized census, matching the combination route.  Coverage and detail keep
their own working-set admission.  The route skips direct primitive import and
rejects incomplete coverage before the unrestricted serial fallback.  Other
primitive/LoD-off imports and explicit oversized requests retain their original
admission checks.  `libBObol_source_realization` now exercises a real detached
source, automatic admission with LoD on/off, and a cold serialized lazy contract.

`/tmp/qged-obol-bare-coverage-20260905` uses the shared warm Lucy mesh cache and
private settings/display.  All four bare-root probes pass strict trace and
semantic final-image checks:

| Backend/mode | Elapsed | Control transitions | Exact presentation |
|---|---:|---:|---|
| System GL shaded | 5.694 s | 194 | 1,284,926 faces |
| System GL wire | 4.662 s | 84 | 3,854,778 lines |
| OSMesa shaded | 4.091 s | 124 | 1,284,926 faces |
| OSMesa wire | 5.118 s | 143 | 3,854,778 lines |

Each final view has one logical occurrence/entry/payload, exact ready geometry,
zero source/occurrence failures, boxes, terminal proxies and prominent-floor
misses.  Peak governed work is 256--643 MiB under the unchanged 1 GiB limit;
this is ledger evidence, not a measured process-RSS bound.  Cold giant and
multi-root contention remain unqualified.  Seven focused CTests pass in 2.62 s;
the `moss.g`, `cube.g`, `prim.g` corpus passes 156 shaded/wire checks over 78
roots.  All 44 model/config pairs pass lint.  Focused TLC explores terminal
convergence (1,187 distinct states), terminal quality ordering (10,232) and
static quality (128), all passing.  Results are in
`/tmp/brlcad-obol-terminal-composition-tla-20260905` and
`/tmp/brlcad-obol-static-quality-tla-20260905`.  The initial sandbox attempt
could not open Java's local listener; the authorized rerun uses local socket
access.  The broad canonical run was interrupted before completion, so this
is not a new full-suite TLC result.
The earlier desktop frame-delivery failure and full Lucy quality gate remain
open.

The subsequent shortened budget diagnostic is retained at
`/tmp/qged-obol-close-budget-20260905` (111 events, 47.332 s, 589 passing
control transitions).  Its final exact ready close view has 1,806,720 faces and
normalized error 4.354, with one prominent-floor miss.  Existing trace knobs
write to **stdout**, captured by qged's command logger.  At view 28/policy 26,
the trace selects a 4,333,525-cost protected allocation, later admits the
6,942,511-cost pixel-demand allocation, then reconciles a rejected static
ordinal trial to 2,516,547 with protected allocation disabled.  The preceding
deadline traces still report zero floor misses/rejection.  This identifies a
boundary to investigate, not proof that the earlier floor remained affordable:
page publication, canonical allocation cost and actual presented work must be
matched before retaining a stronger numeric claim.  Follow-up inspection
confirms one allocation candidate per logical progressive payload and one
common active cut for its presentation layers.  Pages have independent
residency, not independently allocated presentation cuts.  Layer count alone
therefore does not justify changing single-occurrence terminal routing.

### September 5 follow-up: resident page presentation repair

The shortened software profiler at `/tmp/qged-obol-close-perf-20260905`
failed at event 25 after a 60-second ready timeout and 21,244 control
transitions.  `/tmp/qged-obol-planning-cycle-20260905` reproduced the planning
loop with a 15-second timeout: the same selected allocation remained unapplied
while capacity revisions advanced, with no worker or requested frame.  These
are planning failures, not useful close-zoom renderer profiles.  Other repeats
completed; synchronous diagnostic logging changes the timing.  Compact request
names are materialized lazily, so use `BOBOL_LOD_TRACE_OBJECT='*'` to trace
retained admission before provider submission supplies the object name.
The trace filter accepts the failed 15-second report as revisions keep changing;
the script's ready timeout remains a failure.  Audit those revision witnesses
before treating a passing trace alone as evidence of convergence.

A deterministic regression identified a missing presentation guard in the
retained submission fast path.  Updating demand at an unchanged cut can succeed
without installing the required page geometry.  The path nevertheless skipped
the provider whenever its mesh bytes were already resident.  It now requires
prepared presentation at the admitted cut as well as any requested resident
prefetch.  Missing geometry reaches the existing bounded worker repair even
with prefetch disabled and zero additional draw allowance.  No validator,
quality floor, memory bound or logical occurrence count was relaxed.

The regression failed before the fix with zero tasks and zero cut updates.
Afterward it prepares the missing pages exactly once, preserves the allocated
cut and settles; shaded, wire and hidden-line variants pass.  Final focused
validation passes nine CTests in 2.75 s.  The existing
`ObolRetainedAllocationPresentation` model passes TLC with 57 distinct states
at `/tmp/brlcad-obol-spatial-repair-tla-20260905`; it already separates cut
application from spatial presentation readiness.  No model change was needed.

The warm 111-event zoom repeats are at
`/tmp/qged-obol-spatial-repair-20260905/{osmesa,system}`:

| Backend | Script time | Strict control trace | Final presented faces | Final normalized error |
|---|---:|---|---:|---:|
| OSMesa, verbose admission diagnostics | 56.349 s | 1,472 transitions pass | 2,619,533 | 3.581 |
| System GL, ordinary diagnostics | 42.194 s | 466 transitions pass | 510,112 | 8.709 |

Both finish exactly presented and view-ready, with one prominent quality-floor
miss.  These repeats support the preparation fix but do not qualify visual
quality, establish a timing comparison, or retire every historical planning
failure.  The separate desktop frame-delivery stall remains open.

The final-binary full software lifecycle is retained at
`/tmp/qged-obol-spatial-repair-full-20260905`.  All 217 events complete in
111.814 s, all 1,567 control transitions pass, and the endpoint is ownerless,
has no host work, and presents 2,101,208 faces.  The matrix still **fails**:
`quality.prominent_floor` rejects the smooth close stable checkpoints at
events 106/108/110, normalized error 3.581.  Coverage/presentation and the other
validators pass.  This is a full lifecycle control result, not production
qualification or an end-to-end performance comparison.

### September 5 follow-up: constant face lighting in software batches

Obol's `CadRendererGLFlat.cpp` now selects `GL_FLAT` only when every batch mesh
has no explicit normals, all supplied lights are directional, and the fixed
lighting model uses an infinite viewer.  The batch supplies one face normal
and one material per triangle, so its lit vertex colors are constant.  This
avoids redundant pixel color interpolation and restores the caller's shade
model after drawing, including interrupted draws.  Authored/smooth normals,
positional lights and local-viewer lighting keep the existing shading model.
No LoD thresholds, frame deadlines, allocation policy or TLA+ ownership changed.

`FlatFaceLightingMatchesExplicitNormals` compares entire nonempty images with
explicit equivalent normals, including opposite windings, directional/point
lighting, local viewer, mixed smooth normals and an incoming flat shade model.
It checks GL state restoration and requires the flat atlas preparation/tier;
128 distinct parts are necessary to select this route in the ordinary-mesh
fixture.  The corrected reference also passes against the preserved original
library.  Explicit software selection (`OBOL_TEST_RENDER_BACKEND=swrast`) is
required in this dual-backend build; the executable defaults to native GL.
Seven focused software tests pass, and six System GL tests pass with the
software-only fixed-cut test skipped.  Nine BRL-CAD focused CTests pass in
3.01 s.  No new full-suite TLC, sanitizers or platform qualification is claimed.

The reviewed Obol build and BRL-CAD installed library both hash to
`e02411a84f3c6a98e5ba4ad0dbcc992437922991bb25916afa0bfa539a031a87`;
OSMesa remains `097aac27...`.  Obol's two source/test files are uncommitted.
The original `dd3ac7b1...` library, build/install logs and focused test results
are under `/tmp/obol-flat-face-lighting-20260905`.  qged was rebuilt after the
explicit install into BRL-CAD's `.build`.

| Warm replay | Result |
|---|---|
| Original close profile, `/tmp/qged-obol-close-render-perf-20260905` | 111 events, 54.281 s, control pass. Protected cut 29/cost 4,333,525 draws in 419.763 ms with preparation=0 and is rejected by the 400 ms deadline; settles cut 28, 2,619,533 faces, error 3.581 |
| Optimized close profile, `/tmp/qged-obol-flat-lighting-perf-20260905` | 111 events, 62.113 s, all 840 transitions pass. Close cut 29, 3,445,998 faces, error 2.859, no floor miss, 372.390 ms at the stable checkpoint |
| Full software, `/tmp/qged-obol-flat-lighting-full-20260905` | 217 events, 66.049 s, all 1,241 transitions and Lucy geometry/quality checks pass. Close cut 31, 3,445,998 faces, error 1.429, 385.501 ms. **Matrix FAIL** on separate HUD pixel check |
| Full System GL, `/tmp/qged-obol-flat-lighting-system-20260905` | 217 events, 94.470 s, all 860 transitions pass. **Matrix FAIL** on close prominent floor: cut 27, 1,806,720 faces, error 4.354. Final frame ready and ownerless, 3,550,842 faces |

The close profiles use the same 199 Hz perf/deadline/budget instrumentation,
private display/settings and shared warm cache.  The optimized profile records
`flat_rgba_z32_triangle`, confirming the intended raster path executes.  These
adaptive replays do different work and take different total times; the result
records a close-floor improvement, not a general percentage speedup.  A
subsequent route audit shows both close checkpoints use tier 1 (retained cuts),
so this batch-only change does not prove a direct close-up rendering gain.
Preparation at protected cut 29 still includes a 1,159.929 ms interrupted frame
and needs separate cancellation-latency work.

The full software row's failing event 8 (`ae90-0200ms.png`) visibly contains
its short amber fill.  The captured geometry spans 12.743 pixels at width 7,
about 89.2 pixel area; rasterization produces the 91 pixels reported by the
validator.  Its newer estimated fraction is 0.1533, while the captured fill
occupies about 0.0486 of the track.  The fixed 100-pixel minimum consequently
rejects a correctly visible short bar.  This observation contract remains to
be repaired with captured-geometry and hidden-fill regression evidence; the
failed artifact and validator remain unchanged.  Software normal-policy
checkpoints all retain cut 23 and 1,284,926 faces, with approximately
188/189/243/195 ms for flat/authored/smooth/return.  They pass the existing
image reversibility check, but are not the older cut-24 workload.

### September 5 follow-up: captured HUD geometry and retained-cut lighting

The standalone `lod_hud_contract.sh` chooses the first refining checkpoint
whose recorded fill contains at least 100 interior pixels.  It then retains
the original 100 phase-colored pixel minimum and crop/tolerance.  A newer
controller estimate cannot select a too-short rendered fill.  No eligible
capture explicitly reports a skipped visibility check; missing images,
malformed fill geometry, gray-only tracks and insufficient actual colored
pixels fail.  A pixel failure cannot fall through to a later successful image.
`qged_lod_hud_contract` exercises these cases, including the observed 91-pixel
bar followed by a larger visible fill.  All three offline qged contracts pass.
The archived failing HUD artifact passes this repaired check with 1,358 pixels
in its later first-coverage checkpoint; the original report/log is preserved.

A clean full software replay at `/tmp/qged-obol-hud-contract-final-20260905`
completes all 217 events in 102.957 s and passes all 1,998 control transitions
and Lucy quality checks.  The new HUD helper independently counts 1,365 pixels.
The matrix **fails before the HUD gate** on the normal-policy comparison:
flat and authored report active cut 25/3,183,110 faces, but authored actually
presents 2,101,208 faces; smooth and returned authored settle at cut 24 with
that smaller face count.  All four checkpoints report ready.  This needs a
controlled same-cut qualification and an allocation/presentation audit, not a
looser pixel comparison.  The earlier `hud-contract-full` run was terminated
because it overlapped a CMake refresh/rebuild and is excluded from evidence.

The System GL profile identifies the isolated display as Mesa 25.2.8 llvmpipe,
LLVM 20.1.2.  Its close view uses the fixed-function retained-cut executor
(tier 1), even though earlier replay frames spend time in the flat atlas.
The same tier occurs in the prior software close profiles.  Obol now applies
the constant-face lighting optimization to retained cuts as well, factoring
its directional-light/infinite-viewer eligibility into one private helper.
Authored normal items retain the caller's interpolation independently of
normal-free items.  State restoration covers interrupted loop exits.

The exact-image fixture now separately requires the batch atlas and retained
cut preparation kinds.  A finite quantization cut forces both normal-free and
explicit-normal references through retained cut preparation.  Both references
pass against the preserved pre-change library, and eight focused software /
seven System GL tests pass against the new implementation (one software-only
cache test skips on System GL).  Ten BRL-CAD CTests pass in 2.97 s.  The build
and explicitly installed Obol library hash is
`5400bd138f8e0a076c9deff6193edaf113f161989ab2d79dd07023b52eb0fc84`;
OSMesa remains `097aac27...`.  Logs and the preserved `e02411a8...` library are
in `/tmp/obol-fixed-face-lighting-20260905`.  Five Obol files are uncommitted;
no public layout, controller/TLA+ model, deadline, or allocation limit changed.

| System GL close profile | Result |
|---|---|
| Before retained-cut change, `/tmp/qged-obol-system-close-perf-20260905` | 111 events, 37.904 s, all 472 transitions pass. Cut 28 rejected at 439.327 ms with preparation=0; close cut 27, 1,806,720 faces, error 4.354, 354.316 ms |
| After, `/tmp/qged-obol-fixed-face-system-perf-20260905` | 111 events, 38.790 s, all 482 transitions pass. Cut 28 rejected at 442.349 ms with preparation=0; close cut 27, same faces/error, 351.844 ms |

This pair proves neither a meaningful System GL speedup nor close-floor
qualification.  Retain the unchanged quality/deadline failures.  Profiling
still finds substantial atlas construction and driver work; the close snapshot
also retains a roughly 512 MiB batch allocation alongside fixed cut buffers.
Investigate actual executor selection and preparation/reuse costs without
lifting the renderer's memory bound.

### September 5 follow-up: fractional preparation retry identity

The full software replay with the retained-cut lighting change, before the
controller retry correction, is preserved at
`/tmp/qged-obol-fixed-face-full-osmesa-20260905`.  It times out at event 105
after the unchanged 60-second readiness wait, with 106 samples in 99.337 s.
No workers remain; the owner is presentation (6), host flags are 7, and the
requested successor is `render-preparation-replay`.  This is active target
churn: event 105 records 487 transitions and 79 preparation signatures while
the view/policy remain fixed and frames continue to arrive.  The final
snapshot reports 124/124 prepared units, 141 interrupted frames, and a
626.225 ms last interruption.  All 1,809 strict control transitions pass,
exposing another limit of revision/target changes as liveness evidence.

`notePresentationRenderInterrupted()` retained the integer ceiling but
republished it with the default zero next-cut fraction.  That discarded the
page subset being prepared; completing the cheaper whole cut could start
another fractional candidate.  The handler now preserves both coordinates
when the ceiling is unchanged.  A real integer-ceiling correction still
starts at zero fraction.  No deadline, allocation budget, renderer memory
bound, quality floor, or trace checker was relaxed.

The corrected-binary full software matrix at
`/tmp/qged-obol-fractional-retry-20260905` **passes** all 217 events in 98.545 s,
all 1,960 control transitions, Lucy quality, normal-policy coherence/images,
and the HUD visibility gate (1,365 pixels).  The close checkpoint presents
3,863,263 faces, normalized error 2.859, no prominent-floor miss, and a
397.978 ms last frame.  The final view is ready, ownerless, and has no host
work flags.  These adaptive populations/timings are not a same-work speedup
measurement, and this pass does not retire the earlier normal-policy failure.

The actual interruption sequence supplies the regression evidence: at event
105, transitions 1312--1315 request/claim the fractional candidate; 1316
publishes its preparation target; 1317 requests the interrupted preparation
retry; 1320 completes it and requests one reusable timing frame; and 1323
releases the presentation owner.  The preparation target stays unchanged
through those retries and acceptance.  This is an exercised controller path,
not merely a replay which avoided the fractional branch.  No isolated C++
fixture for arming that private static-quality state was added.

Ten focused BRL-CAD CTests pass in 2.98 s.  The existing
`ObolPresentationPreparation` model passes TLC (526 distinct / 1,143 generated
states), and all 44 model/config pairs pass lint.  That finite model already
requires immutable preparation identity; its success does not independently
prove the C++ fraction field.  Build/test logs are
`/tmp/brlcad-obol-fractional-retry-{build,tests}-20260905.log`, with TLC results
under `/tmp/brlcad-obol-fractional-retry-tla-20260905`.  Obol/OSMesa identities
remain `5400bd13...` / `097aac27...`.

The full System GL replay before the controller retry correction,
`/tmp/qged-obol-fixed-face-full-system-20260905`, completes 217 events in
93.654 s and passes 854 transitions, but **fails** the close prominent floor
(cut 27, 1,806,720 faces, error 4.354).  This remains a separate quality issue;
the desktop lost-frame and earlier planning-revision reproducers likewise
retain their own unresolved evidence.

The same final binaries were then replayed on System GL at
`/tmp/qged-obol-fractional-retry-system-20260905`: all 217 events finish in
102.925 s, and all 916 strict transitions pass with final readiness, no owner,
and no host work.  The matrix still **fails** the close prominent floor at
events 106/108/110 (cut 27, 1,806,720 faces, error 4.354, 353.815 ms last frame).
The early quality rejection prevents the later normal-image/HUD gates from
running.  Recorded normal checkpoints use cut 26 for flat/authored/returned
and cut 24 for smooth, so this is not a controlled same-cut normal comparison.
The software pass and System GL failure share the same source/library state;
broader backend quality and repeatability remain unqualified.

### September 5 follow-up: static capacity-search deadline identity

An opt-in System GL programmable-renderer experiment exposed a controller
regression before any renderer default was changed.  Both cases below use
`OBOL_CAD_SOFTWARE_GLSL=1` on the private llvmpipe display, the same certified
warm Lucy cache, the 111-event close-zoom script, and unchanged deadlines.

| Controller state | Result |
|---|---|
| Before, `/tmp/qged-obol-system-glsl-close-20260905` | 111 events in 22.782 s; 440 strict transitions pass. Final cut 21, 134,950 faces, error 17.417, prominent-floor violation, 18.479 ms frame |
| After, `/tmp/qged-obol-static-domain-glsl-close-20260905` | 111 events in 30.832 s; 430 strict transitions pass. Final cut 29, 3,445,998 faces, error 2.859, no floor violation, 376.402 ms frame; ready, ownerless, no host work |

The failed trace first completes a richer static population, then starts a
new capacity search under the STEADY goal with preferred and maximum targets
both set to 400 ms.  That applies the earlier 258,437-unit cadence ceiling
and eventually discards the rich population.  The coordinator now constructs
one complete search key for both completed-pass and interrupted-frame entry:
it retains the ordinary preferred duration and explicitly starts STATIC for
an active static-quality trial or retained ready-view handoff.  Existing
active searches retain their frozen key.  Negative ceilings remain partitioned
by deadline; neither deadline nor budget limits were increased.

The deterministic `test_static_quality_capacity_search_domain` exercises the
production coordinator key and bounded capacity certificate.  A preferred
deadline miss cannot coarsen a proven static population, and a genuine static
miss still bounds the next search.  Ten focused CTests pass in 3.01 s.
`ObolCapacitySearch` now explicitly includes a `staticTrial` handoff origin;
TLC passes 357,336 distinct / 639,650 generated states, and its individual
baseline entry is updated.  `ObolStaticQuality` passes its existing 128 states;
all 44 model/config pairs pass lint.  These are focused checks, not a new
full-suite result.  Logs are `/tmp/brlcad-obol-static-search-domain-*20260905*`
and `/tmp/brlcad-obol-static-search-entry-tla-20260905`.

Default System GL at `/tmp/qged-obol-static-search-domain-system-20260905`
completes all 217 events in 91.338 s and passes all 832 transitions, with final
readiness and no owner/host work.  The matrix still **fails** close quality:
cut 27, 1,806,720 faces, error 4.354, 351.133 ms last frame.  Its later
normal-image/HUD gates do not run after that quality rejection.  The opt-in
programmable result is useful cross-path evidence for the controller fix;
it does not establish a generally faster or fully qualified renderer default.
Obol/OSMesa remain `5400bd13...` / `097aac27...` with no renderer source changes
in this follow-up.

The default OSMesa row
`/tmp/qged-obol-static-search-domain-osmesa-20260905` **passes** the complete
matrix: 217 events in 103.681 s, all 1,481 control transitions, Lucy quality,
normal-policy coherence/images, and HUD visibility (1,365 pixels).  Close cut
29 presents 3,445,998 faces, error 2.859, no floor violation, and a 351.262 ms
last frame.  Flat/authored/smooth/returned checkpoints all retain cut 24 and
2,101,208 faces, with approximately 269/293/336/274 ms frames.  The final view
is ready with no owner or host work.  This validates the default software path
after the shared controller fix without retiring the separate System GL miss.

The full opt-in programmable System GL row
`/tmp/qged-obol-static-search-domain-glsl-full-20260905` completes 217 events
in 76.798 s and passes all 914 transitions, Lucy settled quality, normal-policy
coherence/images, and HUD visibility (364 pixels).  Close cut 29/error 2.859
passes at 374.834 ms.  Normal checkpoints all use cut 24/2,101,208 faces.
The matrix nevertheless **fails** the later smooth-zoom responsiveness gate:
none of the four in-gesture captures satisfies `spatially_advanced`.  The
starting active cut is 26 with renderer ceiling 25; the captures have rendered
cut bounds 25/23/23/23, and all retain the starting 262,094,432 resident bytes.
Only the later quiet close view reaches the richer population.  The retained
`gesture-failure.json` names the failed witnesses.  This is the next
programmable-path refinement reproducer, not permission to waive the gate or
change renderer defaults.  Final readiness is ownerless with no host work.

### September 5 follow-up: motion recovery and failed-cut evidence

The programmable continuous-zoom reproducer exposed another numeric-policy
boundary.  Cost-based deadline recovery can skip several cut ordinals, but
`noteQualityMiss()` treated the cheap fallback as the maximum safe probe.
Intermediate cuts were rejected without measurement.  Completed and aborted
frame paths now supply the actually attempted cut separately; only that cut
and richer cuts are excluded.  The recovery hint, deadlines, memory limits,
renderer defaults, and existing cheaper-view/epoch witness retirement remain
unchanged.  `BOBOL_LOD_TRACE_HEADROOM` now records the probe decision and its
floor, retry limit, and work witnesses.

The coordinator regression recovers from failed cut 21 to cut 16, then requires
completed-frame witnesses for probes 17--20 while rejecting 21.  Ten focused
CTests pass in 3.39 s; logs are
`/tmp/brlcad-obol-motion-miss-witness-{build,test-build,tests}-20260905.log`.
After comment cleanup and removing a redundant fixture assignment, the final
58-target rebuild and all ten CTests pass again (9.17 s); logs are
`/tmp/brlcad-obol-motion-miss-witness-final-{build,tests}-20260905.log`.
This changes a numeric retry boundary within the existing frame barrier;
no new formal model or full-suite TLC result is claimed.

The instrumented opt-in 111-event replay
`/tmp/qged-obol-motion-miss-witness-glsl-close-20260905` completes in 23.984 s
and passes all 426 strict transitions.  Its trace actually advances ceilings
20/21/22 before vetoing 23, then permits further work after the existing
cheaper-view witness changes.  Final cut 29 presents 3,445,998 faces, error
2.859, no floor miss, at 383.789 ms.

Full warm rows use the same binaries and cache, one private display at a time:

| Renderer and artifact under `/tmp` | Full result |
|---|---|
| Default OSMesa, `qged-obol-motion-miss-witness-osmesa-20260905` | **PASS**: 217 events, 69.244 s, 1,270 strict transitions, quality, normal-policy images, smooth zoom and HUD (1,358 pixels). Close cut 29/3,445,998 faces, error 1.429, 341.719 ms. Normal checkpoints all use cut 23/1,284,926 faces, approximately 189/213/244/187 ms. |
| Default System GL, `qged-obol-motion-miss-witness-system-20260905` | **FAIL**: 217 events, 92.909 s, 778 strict transitions pass. Close cut 27/1,806,720 faces still misses the floor (error 4.354, 329.981 ms). Event 32's first-mesh checkpoint also reports zero presented geometry and `presented_cad_work_exact=false` despite an active mesh. Later image/HUD gates are not reached. |
| Opt-in programmable System GL, `qged-obol-motion-miss-witness-glsl-full-20260905` | **FAIL**: 217 events, 73.107 s, 900 strict transitions, settled quality, normal-policy images and HUD (364 pixels) pass. Close cut 29/3,445,998 faces has error 2.859 at 376.916 ms. Continuous refinement advances, but the atlas-generation reuse condition still rejects the row. |

Every final view is ready with no owner or host work.  The opt-in row now
passes `spatially_advanced` and `spatially_realized` at events 100/102: actual
cut bound 27 exceeds starting cut 26, resident bytes grow from 262,094,432 to
279,277,676, cache loads increase from six to seven, and the frame takes
39.479 ms against 100 ms.  All 13 other outer zoom conditions pass, including
quiet restoration, memory turnover and wheel dispatch.

The remaining nested condition requires either atlas generation reuse,
ordinary generation reuse, or an ordinary-part first upload.  Atlas reuse
stays at 1,160, ordinary reuse stays zero, and upload counters remain unchanged.
Both captures use indirect renderer tier 6; the atlas counter in
`CadGpuResources.cpp` increments on a generation change, not every successful
redraw.  This identifies an observation boundary to investigate, not permission
to bypass its incremental-realization proof.  Retain `gesture-conditions.json`
and `gesture-failure.json` in the opt-in artifact.  The validator is unchanged.
Also retain the new default System GL first-mesh observation and the older
normal-policy, planning and desktop failures; this motion fix does not close
those mechanisms or qualify programmable rendering across scales.

### September 5 follow-up: indirect atlas-ceiling admission

The preceding zoom failure had a real renderer component.  The indirect
occurrence-cut patch admitted missing GPU suffixes, but the renderer-ceiling
patch only clamped commands to the already uploaded prefix.  Raising the
ceiling without changing the occurrence cut could leave valid CPU detail
unused and report atlas pressure without attempting admission.

`IndirectCeilingGrowsResidentAtlas` reproduces the defect on the old library:
128 triangles remain drawn instead of the required 38,400.  The two patch
paths now share `prepareIndirectAtlasPrefix()`.  It retains the existing
memory governor, coherent-prefix fallback under real pressure, and exact
preparation fallback when relocation changes retained command offsets.
Both occurrence-cut and ceiling-only regressions require the expected counts,
a nonempty rich image, a different coarse image, and exact rich-image
restoration without another upload or generation change.

The old failing test is `/tmp/obol-atlas-ceiling-before-system-20260905.log`.
Ten focused System GL tests pass, including actual indirect execution,
pressure/proxy bounds, abort validation and lighting; six software tests pass.
Build/install/test logs are in `/tmp/obol-atlas-ceiling-20260905`.  Installed
Obol and its build tree both hash `f281b59b...`; the preceding `5400bd13...`
library is retained under `baseline`.  OSMesa remains `097aac27...`.

The first full programmable Lucy replay after this renderer fix,
`/tmp/qged-obol-atlas-ceiling-glsl-full-20260905`, completes 217 events in
74.630 s with 909 transitions and passing settled quality, normal images,
and HUD (364 pixels).  During zoom, atlas suffix uploads advance from
59,448,624 to 83,509,080 bytes while generation reuse stays at 1,158.
Its original validator still rejects that unchanged generation counter.
`lod_atlas_reuse.jq` now recognizes measured same-generation suffix uploads;
13 adversarial cases reject unchanged data, full replacement alone, reservation
growth alone, and missing/malformed counters.  The complete revised zoom
filter accepts this report and still rejects the prior no-upload reproducer.

A separate first-mesh observation bug was exposed by default System GL.
`wait_progressive_cad_mesh_ready` checked CPU adoption before delivering a
queued frame and could return after that frame aborted.  It now inspects
exact completed geometry after frame delivery, maintaining the existing
quiet interval and deadlines.  All eleven focused BRL-CAD CTests pass (3.00 s).
Two accidentally overlapping rebuilds were stopped/reconciled before GUI
qualification; the sequential 94-target rebuild succeeds.  No graphical
timing row overlapped that build work.

Final all-changes warm rows use the same qged (`2b28af79...`), installed
libraries and cache, with isolated displays and no concurrent builds:

| Artifact under `/tmp` | Result |
|---|---|
| `qged-obol-atlas-ceiling-glsl-final-20260905` | **FAIL**: 217 events in 69.874 s, 860 strict transitions pass. Close cut 29/3,445,998 faces/error 2.859 passes at 383.222 ms. Quiet zoom-out instead settles at cut 20/337,610 faces/error 4.274 and 41.700 ms. |
| `qged-obol-atlas-ceiling-system-final-20260905` | **FAIL**: 217 events in 100.532 s, 929 strict transitions pass. The first-mesh checkpoint is exact with 340,082 faces. Close cut 27/1,806,720 faces/error 4.354 still misses the floor at 346.632 ms. Later normal-image/HUD gates are not reached. |
| `qged-obol-atlas-ceiling-osmesa-final-20260905` | **FAIL**: 217 events in 104.708 s, 1,460 strict transitions, settled quality, normal images and HUD (1,414 pixels) pass. Close cut 29/error 2.859 takes 340.344 ms; all normal snapshots use cut 24/2,101,208 faces. Continuous zoom fails its rendered-cut advancement condition. |

Every final view is ready with no owner or host work.  The programmable row's
revised in-gesture branch passes, including its actual suffix upload.  Its only
failed outer zoom conditions are the quiet prominent floor and zoom-out
residency settlement.  At that checkpoint active cost is 1,067,000 against a
10,270,264-unit budget, with no renderer ceiling or GPU pressure; resident
bytes are 328,746,776 and the compaction plan is not current (revision zero).
The payload has a memory-denial witness at admission revision 1 and resident
cut 26 behind displayed cut 20.  That resource fact was missed in the initial
budget-only diagnosis.  Retain `gesture-conditions.json` and the follow-up
below; the original row does not justify a longer render deadline.
Later normal/HUD gates do not run after the quality rejection.

The software row starts at cut 25 and captures in-gesture rendered bounds
25/23/25/25.  Resident bytes grow from 262,094,432 to 279,277,676 and cache
loads from five to six, but none of those captured bounds exceeds the starting
cut.  The bounded-discrete exception is also false because requested and
active cuts both equal 28.  Its other 13 outer zoom conditions pass; quiet
zoom-out compacts to 193,515,360 bytes with a current plan.  This differs from
the programmable terminal failure.  All three first-mesh checkpoints now
report exact nonzero geometry.  Keep these failures and renderer defaults;
no full-stack production or cross-scale programmable qualification is claimed.

### September 5 follow-up: prepared residency after a memory denial

The quiet zoom-out investigation found that a current memory-denial witness
restricted the occurrence allocator to `activeCut`, even when a worker had
already loaded and prepared a richer prefix behind that displayed cut.  The
capacity search consequently classified unchanged coarse geometry throughout
its numeric budget domain.  A denial still constrains new source growth; it
must not erase usable immutable data.

The allocator now considers the richest requested cut supported by resident
data and prepared renderer geometry.  Spatial demand requires every populated
page, and the existing fixed bounded-coverage exception remains unchanged.
The submission action permits that certified prepared allocation through the
memory-denial gate without enabling provider retries.  No memory limits,
render deadlines, visual floors, or renderer defaults change.

The expanded `libBObol_lod_update_action` regression fails before the fix on an
already prepared richer suffix.  All three cases pass afterward: unavailable
suffix, resident-but-unprepared suffix, and prepared suffix.  The first two
retain the denial; the last advances without a provider task.  Twelve focused
CTests, including the service, allocation oracle and coordinator, pass in
8.23 s.  Logs use `/tmp/obol-resident-denial-{before,after,final-tests}-20260905.log`;
sequential build logs use the same prefix.  No control-state transition or
formal model changed; no new TLC run is claimed.  The allocation oracle and
GUI control ledger cover this application change.

The new full opt-in System GL warm row,
`/tmp/qged-obol-resident-denial-glsl-20260905`, **passes every matrix gate**:
217 events in 78.677 s and 928 strict transitions.  Close detail uses renderer
ceiling 29/error 2.859 at 390.679 ms.  Quiet zoom-out reaches cut 25/3,162,945
faces/error 1.339 at 347.882 ms, with one compaction and resident bytes falling
from 327,791,888 to 275,887,428.  Its current-plan flag is false at that
checkpoint, but measured reclamation satisfies the unchanged turnover gate;
the return-view plan is current.  Normal snapshots all use cut 24/2,101,208
faces, and the normal, HUD, interaction, first-mesh and camera gates pass.
Final readiness has no owner or host work.

The companion default OSMesa row,
`/tmp/qged-obol-resident-denial-osmesa-20260905`, also **passes every matrix
gate**: 217 events in 69.911 s and 1,263 strict transitions.  Its starting cut
is 23; the active zoom reaches occurrence cut 26 behind renderer ceiling 25.
Close cut 29/error 1.429 and quiet zoom-out cut 24/error 0.814 pass.  One
compaction reduces resident bytes from 222,923,888 to 95,123,820, with a current
plan.  Normal snapshots all use cut 23/1,284,926 faces.  Final readiness again
has no owner or host work.  This does not retire the earlier software
cut-25-start input reproducer, which follows a different adaptive path.

These are two full passing rows, not proof of schedule-independent graphical
behavior.  The preceding diagnostic replays also reached cut 25 before the
fix; retain the original failing report and deterministic regression rather
than attributing all timing differences to the change.  The short diagnostic
artifacts are `/tmp/qged-obol-zoomout-headroom-trace{-b,-c}-20260905` and
`/tmp/qged-obol-zoomout-memory{,-shallow}-trace-20260905`.  Current qged remains
`2b28af79...`, Obol `f281b59b...`; libBObol is `ce8d39fa...` (full hashes in
`/tmp/obol-resident-denial-20260905/hashes.txt`).

### September 5 handoff audit: software optimization and remaining failures

The repeated-normal optimization is implemented in OSMesa commit `bfa90de`.
Only the constant-material, multi-directional-light fast path reuses the
preceding vertex's result, after a bitwise-equal normal comparison.  Reuse is
local to one lighting call; per-vertex material updates and general/positional
lighting retain their independent paths.  The new `osmesa_lighting_reuse`
test compares whole nonempty RGBA images against the material-processing path
with equal constant per-vertex material values.  Its 16 cases cover repeated
and varying normals, one/three lights, one/two-sided lighting, normal
transforms, and material changes between draws.  It is not a test of arbitrary
varying per-vertex materials or a timing-threshold CI test.

Three alternating baseline/optimized microbenchmark pairs preserve every image
hash.  Normalizing against the reference loop within each case suggests about
15--17% improvement for repeated-normal multi-light batches, with varying
normals essentially unchanged.  Keep the raw paired timings and baseline
library in `.build/testing-artifacts/osmesa-lighting`; absolute timings showed
enough machine noise that they must not be sold as an end-to-end speedup.

The retained warm GUI/perf run at
`/tmp/qged-lucy-lighting-reuse-perf-20260904` completes all 217 events in
114.466 seconds (117 seconds in the harness summary).  All 1,504 transitions
pass the current complete-trace checker, with no dropped/truncated records.
The final frame is ready, ownerless, obligation-free, and presents 2,101,208
faces at cut 24 with no prominent-floor miss.  Flat/authored/smooth/returned
normal checkpoints all retain that same cut and face population; recorded
frame times are approximately 316/335/365/305 ms.  This is useful same-cut
evidence, not yet repeated cross-backend image qualification.

The overall result is still **FAIL**, but its first diagnostic needs correction:

- The composite `Lucy did not converge to a recognizable resident PoP mesh`
  check requires final compact-registry entries, but perf mode disables deep
  diagnostics with `QGED_TEST_DEEP_LOD_REPORT=0`.  In
  `QgProgressiveDiagnostics.cpp` those counters are incremented only inside
  `collectDeepLodDiagnostics`; the serializer still emits their initialized
  zeros.  The validator mistakes unmeasured fields for measured empty
  inventories despite one active CAD payload, one resident service asset, and
  a real rendered mesh.  Fix the diagnostic-availability contract, and retain
  semantic coverage checks for both bare-object and combination routes.
- The same check rejects a proxy in the very first draw-return frame (event
  4).  Determine whether this is an unacceptable initial cue or an incorrectly
  scoped settled-frame assertion; do not silently exempt later regressions.
- Its background-cache checkpoint actually has zero prominent-floor misses.
  The preceding `/tmp/qged-lucy-normal-perf-20260904` report also fails the
  unmeasured-registry/initial-proxy clauses while its background-cache checkpoint
  has zero floor misses and normalized error 0.846.  The old generic failure
  text was not evidence that either initial settled population was too coarse.
- Independent real deficits remain: discrete zoom settles at 264,723 faces
  with normalized error 5.91; smooth close zoom at 1,806,720 faces/error 4.35;
  zoom-out at 555,864 faces/error 3.26.  Each reports one prominent-floor miss
  and a performance constraint.  These are numeric/visual qualification
  failures even though the control trace drains correctly.  Preserve their
  hard quality checks while repairing the observation contract.

At the handoff review, the user updated Obol to `a544ae38`, recording the
same `bfa90de` OSMesa source as the sibling checkout.  Source integration is
therefore no longer open.  Binary provenance still needs reconciliation:
BRL-CAD's installed `libosmesa.so` hashes to `097aac27...`, while Obol's
build-tree copy still hashes to the preserved baseline `0fc61246...`.
Rebuild/install the matching dependency chain before new qualification; do
not infer a common renderer from clean git status.  Obol builds its real
`external/osmesa` submodule, not a symlink to the sibling checkout.

Fresh checks during this documentation review pass the four convergence,
update-action, progress-estimator, and trace-checker CTests (2.47 s), the two
sibling OSMesa lighting/thread-context CTests (0.33 s), and TLA lint for all
44 cataloged model/config pairs.  All 44 retained baseline entries report
pass; full TLC exploration, a clean dependency rebuild, new GUI runs, and
cross-platform qualification were not performed during this review.

## Matrix

### Lower-level/shared stack

- Clean Ninja build of qged, MGED, gsh, Archer/TkObol, rtwizard, and plugins.
- `bobol_contract`, `bobol_headless`, `drawing_baseline`, libbg polygon,
  librt discovery/edit, libbu cache, installed-package consumer, public-header,
  and symbol-manifest tests.
- TLC configurations for occurrence publication, host work, bounded
  Lucy/large-scene convergence, and deferred-autoview ownership.  The latter
  proves control ownership only; real qged cold/warm camera replays remain
  required for Qt/renderer timing and image behavior.
- Shared-stack ASan/UBSan: workers active during teardown, cache corruption,
  endpoint replacement, plugin reload, edit cancellation, and rapid close.
- Native-worker TSan/LSan; do not count a container runtime limitation as a
  successful dynamic analysis run.

### qged controls and interaction

Run each applicable row in single and quad layouts, DPR 1 and fractional DPR,
before/after resize, on System GL and OSMesa.

- Camera: MGED orbit, rotate, translate, center, zoom, smooth round trip,
  autoview, aspect changes, navigation gizmo, and camera/history readback.
  Close orthographic zoom must retain its scene-depth range and never behave
  as an implicit cut; verify this at extreme Hubble zoom on both renderers.
- Selection: point and rectangle, append/toggle/subtract, tree-to-scene and
  scene-to-tree propagation, selection under erase/redraw, and Hubble latency.
- Faceplate: axes/ticks, grid, ADC, HUD progress, lighting profiles,
  framebuffer modes, raytrace overlay/underlay/interlay, snapping, and view
  settings.  The `view cutting` controller is disabled by default and is
  independent of camera clipping; qualify its visible plane affordance,
  faceplate controls, persistence, exact picking, and both renderers.
- Polygon/measurement: direct mouse creation and editing, booleans, snapping,
  persistence, cancellation, command readback, all measurement modes, labels,
  units, and resize.
- Primitive/sketch editing: actual widgets and mouse gestures, CLI/widget/
  manipulator round trips, commit/cancel/revert, invalid-input immutability,
  plugin reload, multiple views, and database mutation invalidation.

### Drawing and LoD

- Modes 0--5 on appropriate moss/rook models, with LoD auto and off.
- Generic Twin shaded/wire cold/warm and LoD auto/off; compare `ae 90 0` and
  oblique views to the main baseline for skin seams, engine, tail, underside,
  lighting, and boxes.
- Lucy cold/warm zoom in/out and rotation; check spatial page coverage,
  crack-free wings, retained history, in-motion refinement, uniform zoom-out
  quality, and resident-memory recovery.
- Hubble shaded/wire selection, resize, erase/redraw, and importance floors.
- NIST BREP shaded/wire, LoD auto/off, adaptive tessellation growth, and
  zoom-out reclamation.
- Multi-Lucy/xpush turnover; bring old and new instances into view while
  retaining shared assets and interaction responsiveness.
- 5k/50k/150k cold and certified-warm shaded/wire runs, including rotation,
  zoom, selection, exact subpath erase/redraw, and retained terminal frames.

### Discovery, performance, and visual significance

- Measure read-only discovery on cold local storage, warm page cache, and
  slower storage.  Record elapsed time, worker concurrency, I/O behavior,
  peak bytes, and tree/GED parity.
- For 50k and a real model, collect perf.  Separate discovery, compact
  planning, worker/cache, upload/preparation, and raster work from diagnostic
  JSON/image compression.
- Capture APNG for any visual cycle; compare terminal and representative active
  frames to the main baseline where available.
- Record time to overview, leaf coverage, first useful mesh, view ready, stable
  return, p50/p95/max input latency, completed-frame time/FPS, resident/peak
  memory, cache growth, and quality-floor debt.
- Evaluate visual-importance results using wheels, blades, tails, hulls, and
  fasteners.  Synthetic terminal coverage alone is not evidence of adequate
  visual quality.

## Evidence retention

Temporary-artifact cleanup on 2026-09-04 reclaimed approximately 40 GiB.
Completed Lucy replay caches were removed, retaining the format-24 shared
cache identified above and all replay reports, scripts, logs, and images.
The caches in `/tmp/obol-current-cache-matrix` were still format 23 and are
not reusable with the current format-24 implementation; they were removed,
but the historical comparison evidence remains at its original paths (39 MiB).
Future stress qualification must rebuild these caches from the retained shared
model inputs; the earlier warm-run evidence is not current-format qualification.
The abandoned canonical TLC state store and two temporary documentation
configure builds were also removed.  Configure diagnostics are archived in
`.build/testing-artifacts/tmp-cleanup-20260904/tla-configure-diagnostics.tar.gz`;
other TLC logs and counterexample evidence were retained.  No shared geometry,
crash records, or other sessions' temporary artifacts were removed.

Keep one shared copy of reusable large inputs and canonical warm caches.  Keep
the event script, final report, representative checkpoint images, and perf/APNG
only when they explain a result.  Remove redundant transient screenshots,
trace logs, and duplicate generated models after recording the needed evidence.
