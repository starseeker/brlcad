# GED draw/view C API contract

This document defines the public vocabulary and behavioral rules for the
renderer-neutral GED draw API.  It is the design contract for the installed
headers; Obol adapters and database-source realization services are private
libged implementation details.

## Vocabulary and headers

| Header | Vocabulary | Responsibility |
| --- | --- | --- |
| `ged/view.h` | `ged_view_*` | Opaque view and view-set lifecycle, host policy, and display endpoint |
| `ged/scene.h` | `ged_scene_*` | Semantic draw intent, scene queries, occurrences, deltas, and edit scopes |
| `ged/view_feature.h` | `ged_view_feature_*` | Read-only queries for managed per-view presentation features |
| `ged/view_feature_batch.h` | `ged_view_feature_batch_*` | Atomic typed feature publication by in-tree producers |
| `ged/view_polygon.h` | `ged_view_polygon_*` | View-owned polygon editing, presentation, and sketch conversion |
| `ged/view_edit.h` | `ged_view_edit_*` | Retained edit-preview publication and lifecycle |
| `ged/view_export.h` | `ged_view_database_export_*` | Renderer-neutral database presentation export |
| `ged/selection.h` | `ged_selection_*`, `ged_pick_*` | Per-view selection and immutable pick results |
| `ged/draw.h` | none | Convenience umbrella for the drawing-domain headers |

Names do not repeat `draw_view_context`: the opaque first argument already
identifies a GED view.  Renderer names do not appear in renderer-neutral
headers.  Realized shape/group/source records, database traversal, LoD
realization, framebuffer internals, and direct controller operations are
private services and are not stable libged API.

The view-polygon domain owns retained interaction state, not polygon
algorithms.  `bg/polygon.h` and `bg/polygon_types.h` are the canonical lower
level representation and topology API for area, boolean operations, overlap,
and triangulation.  A generally useful polygon algorithm belongs in libbg;
libged should only adapt it to a view, a database sketch, or a presentation
feature.

## Ownership and lifetime

`ged_view_context` and `ged_view_set` are incomplete C types.  A view is owned
by its creating application unless it is registered with
`ged_view_context_owned_add`; a GED-owned view is destroyed by `ged_destroy`.
The view does not expose or transfer ownership of its libbv or renderer
implementation.  Its display endpoint is the only public renderer attachment
boundary and follows the ownership rules in
`obol_endpoint_lifecycle.md`.

Scene occurrence, edit, feature, and polygon references are copyable values.
They do not retain their owner.  Occurrence references remain valid across
unrelated scene commits and become stale when that occurrence is retired;
retired identifiers are not reused.  Record pointers and strings returned in
callbacks are borrowed and valid only for the callback, or until the next
documented owner mutation.  Feature batches own deep-copied staged geometry
until commit or abort; either operation consumes the batch.

## Threading

View lifecycle, endpoint, feature, polygon, edit-preview, selection, and pick
calls run on the view owner's thread.  Value-reference resolution asserts this
rule in assertion builds and rejects cross-thread use otherwise.  GED scene
transactions and observer registration are likewise serialized by the GED
owner.  Worker services may compute source geometry, but publication back to a
view is marshalled to the owner thread.  No installed call grants concurrent
mutation merely because its handle is opaque.

## Return and error conventions

Unless a declaration states a count:

- `int` returns `1` for success and `0` for absence, stale input, unsupported
  operation, or validation failure.
- Transactional functions that need to distinguish a hard failure return a
  negative value on failure, zero for no change, and a positive affected count.
- Pointer and value-handle results return their documented null value when not
  found or on failure.
- `size_t` query functions return zero for both an empty result and invalid
  input; callers that need to distinguish those cases use the corresponding
  validation/query operation.
- Functions accepting `struct bu_vls *result` append actionable diagnostics;
  other validation failures do not mutate the GED command result string.

Stale scene, feature, and polygon handles are normal validation failures: the
operation returns failure without dereferencing retired storage.  Mutations
are atomic at their documented transaction/commit boundary.  Observer
callbacks run only after the resulting state and revision are stable.

## Stable surface

`src/libged/ged_draw_api.symbols` is the reviewed draw/view export manifest.
The `ged_draw_api_check` build target extracts exported declarations from the
installed domain headers and fails if a symbol is added or removed without an
intentional manifest update.  Private Obol and librt adapters must never appear
in that manifest.  The current reviewed manifest contains 233 symbols.
Dynamic-symbol checks also reject the former realized
shape/group/source/scene-record families even though their definitions still
exist as private implementation seams.
