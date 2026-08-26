# Live spatial-page publication

## Purpose

A cold large serialized BoT must become useful before its durable PoP cache is
complete.  The bounded coverage preview guarantees a complete extent and an
immediate visual, but it is intentionally too coarse to be the only result
while spatial pages are being built.  A completed spatial page must therefore
be publishable before the final cache marker exists.

## Ownership and representation

Use the existing `SoCADAssembly` multi-part capability.  Do not add mutable
geometry to Coin nodes or combine all pages into one growing GUI-thread mesh.

Each logical CAD occurrence has one retained presentation containing:

1. one immutable coverage part, present until every visible spatial page has a
   richer replacement;
2. zero or more immutable detail-page parts, each backed by a completed
   producer or cache page; and
3. after the hierarchy is complete, the exact resident page set needed by the
   current view at one common active cut.

There is deliberately no final monolithic copy of a multi-page hierarchy.
Warm cache reads and live cold publication use the same retained layer form.
Packing every resident page into one cumulative array made a local page change
copy the rest of a large mesh and retained replaced pages as tombstones.

The coverage part and detail parts share the occurrence transform, selection
identity, color/style, clipping state, and semantic path.  Detail pages have
their own assembly part identifiers derived from `(source content identity,
spatial partition version, page id)`.  They are not new BRL-CAD objects and
must not create rows in the tree, duplicate selection entries, or alter draw /
erase path semantics.

The source owner thread updates only the affected shared part set.  A producer
worker builds `PartGeometry` off thread and publishes an immutable result; it
never touches a Coin node.

## Producer contract

After the serialized stream has been partitioned, a page becomes live only
after all of the following are true:

- its local source-face stream is complete;
- every index is validated and the local cumulative cut metadata is complete;
- its geometry has been copied into an immutable worker-owned page record;
- its page key and partition version are known; and
- cancellation has not been requested.

The producer may leave a durable cache record absent or incomplete while live
pages exist.  The final hierarchy marker remains the only cross-process cache
discovery witness.  A cancelled producer discards unpublished pages and leaves
no discoverable hierarchy marker.

The first implementation uses the existing global quantization schedule.  A
page at cut `n` uses source-coordinate snapping for that schedule, so adjacent
pages rendered at the same cut agree at boundaries.  Unequal-cut stitching is
not introduced by live publication.

## Presentation contract

Coverage is removed only where detail pages provide a replacement, never merely
because a page is being prepared.  The initial implementation may retain
coverage behind detail pages, but it must use a non-occluding coverage draw
mode so detail is visible and must not produce z-fighting.  A follow-up compact
coverage mask may remove covered cells after a bounded batch of page arrivals.

The retained page set is bounded by the same view demand and resident-memory
policy as cache-complete pages.  Page publication is batched by a fixed cost
allowance; it may not enqueue one GUI transaction per source face or bypass
the scene-wide render budget.  Invisible pages may be persisted without being
presented.  On a memory edge, detail pages retire first while coverage remains
available; a later view demand reloads them from cache if durable.

## Implementation sequence

The existing compact adapter deliberately has one `activePart` per logical
occurrence.  Do not bypass that model by adding untracked scene nodes or extra
compact-registry entries: either would duplicate selection/tree identity and
make erase/redraw races inevitable.  Extend the retained presentation record
to own a bounded layer set instead:

1. Add a worker-owned immutable layer value to `BObolLodResult` and
   `CadPayload`: `(layer key, geometry, geometry revision, cut, channel
   metadata)`.  A result with layers replaces its old layer set atomically;
   the legacy single geometry is the one-layer form until it is removed.
2. In the cold mesh preview context, copy each callback page into a
   `PartGeometry` off the GUI thread, retain it by deterministic page key,
   and publish bounded snapshots containing the standing coverage layer plus
   currently visible detail pages.  It must not publish a page as the
   occurrence's only geometry.
3. Let `SoBRLCadAssembly::CompactInstancePresentation` map the one logical
   occurrence to deterministic auxiliary renderer instance IDs derived from
   `(logical instance id, layer key)`.  Those instances are an assembly
   detail, not source occurrences: each maps back to the base
   `InstanceSemantic`.
4. Propagate base transform, effective style, hidden/unpickable/selected
   state, LoD cut, and semantic pick detail to every auxiliary instance in the
   same update transaction.  Tree and selection sets remain keyed only by the
   base occurrence.  A pick of a page reports the base semantic path.
5. Add/remove auxiliary parts only after their reference count reaches zero;
   update page batches through `upsertSharedParts`/`upsertInstances`, never by
   rebuilding the entire compact registry or combining page arrays on the GUI
   thread.  Coverage is a retained non-occluding layer until page replacement
   logic can remove its covered cells safely.

The current page representation uses one active cut for every presented page.
Pages may have different resident capacities, but a page whose capacity is
below the common cut cannot be published in that presentation.  Unequal active
cuts at neighboring page boundaries are not permitted.

## Required API seams

- The mesh-cache producer now exposes `BObolMeshLodSpatialPageCallback` in
  `BMeshLodCache.h`.  It carries a complete deterministic page id, cumulative
  page data at the common terminal cut, and a page-local hierarchy snapshot.
  The record is borrowed only for the callback; the receiver must copy it into
  an immutable worker-owned layer.  It is distinct from a global PoP prefix
  and from coverage points.  The cache regression verifies completeness,
  local index bounds, deterministic page order, whole-source face coverage,
  and cancellation after one delivered page without a durable hierarchy
  marker.
- `BObolLodProgressiveMesh` retains decoded cache pages independently and
  prepares immutable page-local `PartGeometry` lazily.  Unchanged pages are
  weakly cached and adopted directly by the scene; a view-local update neither
  rebases their indices nor copies their arrays into an aggregate.
- `BObolLodResult` and `CadPayload` carry an ordered vector of immutable
  `BObolLodPresentationLayer` values.  Result canonicalization rejects invalid
  or duplicate layer keys.  The single-geometry member is the unpaged one-layer
  form, not an alternative multi-page protocol.
- `SoBRLDatabaseSource` maps layers to existing shared-part updates and
  deterministic auxiliary renderer instances under one compact source
  occurrence.  Scene selection and tree state remain keyed by the logical
  source binding, not by layer part or auxiliary instance id.

## Acceptance tests

1. A synthetic serialized mesh with at least three spatial pages publishes
   coverage first, then a detail page, without visible structural boxes or
   source-order fragments.
2. Replacing a visible page preserves the GPU prefix/lineage for unaffected
   pages and does not copy their CPU arrays on the owner thread.
3. Cancelling after one live-page publication preserves that valid coverage /
   detail presentation and publishes no durable hierarchy marker.
4. Cold Lucy shows coverage promptly and at least one richer visible page
   before durable cache completion; warm Lucy uses cache-selective pages.
5. Selection, subpath erase/redraw, color override, wire/shaded/hidden-line,
   resize, and view replacement remain semantic per logical occurrence.
6. System GL and OSMesa screenshots contain the same visible page set within
   their declared cut/performance tolerance, with no z-fighting or occlusion
   of detail by coverage.

## Control-plane model

`ObolLiveSpatialPublication.tla` is the bounded TLA+ model for this seam.  It
does not claim to model triangle topology, frame duration, or visual error.
It establishes the ownership conditions which must remain true independent of
those data-plane details:

- an emitted page was first complete and validated;
- coverage remains present throughout incomplete publication;
- only a complete final hierarchy may retire coverage and mark the durable
  cache; and
- cancellation freezes publication and cannot create a discoverable cache
  marker.

Run it when changing the producer callback, payload layering, or cache-marker
ordering:

```
java -XX:+UseParallelGC -jar /home/cyapp/tla+/tla2tools.jar -workers 1 \
  -config doc/notes/ObolLiveSpatialPublication.cfg \
  doc/notes/ObolLiveSpatialPublication.tla
```

On 2026-08-25 TLC explored 343 generated / 212 distinct states to depth 11
with no invariant, deadlock, or liveness error.  The retained-page cache and
LoD-service tests also cover page replacement and page-selective compaction.
A warm System GL Lucy interaction replay passed with in-gesture refinement,
box-free terminal presentation, zoom-out compaction, and exact-view history
restoration.  The corresponding OSMesa replay remained correct and box-free,
but exposed a separate software-throughput quality-floor limit; it is not
release evidence for the OSMesa visual-quality row.
