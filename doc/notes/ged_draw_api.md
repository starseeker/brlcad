# GED draw/view C API contract

This document defines the public vocabulary and behavioral rules for the
renderer-neutral GED draw API.  It is the design contract for the installed
headers; Obol adapters and database-source realization services are private
libged implementation details.

## Vocabulary and headers

| Header | Vocabulary | Responsibility |
| --- | --- | --- |
| `ged/view.h` | `ged_view_*` | Opaque view and view-set lifecycle, host policy, and display endpoint |
| `ged/scene.h` | `ged_scene_*`, retained `ged_draw_*` records | Semantic draw transactions, groups, shapes, and scene queries |
| `ged/draw_source.h` | `ged_draw_source_*` | Renderer-neutral source snapshots and source-state queries |
| `ged/view_feature.h` | `ged_view_feature_*`, `ged_annotation_*` | Managed per-view features and annotations |
| `ged/polygon.h` | `ged_polygon_*` | Managed view polygons and sketch conversion |
| `ged/selection.h` | `ged_selection_*`, `ged_pick_*` | Per-view selection and immutable pick results |
| `ged/result.h` | `ged_result_*` | Atomic command-result scene publication |
| `ged/draw.h` | none | Convenience umbrella for the six draw-domain headers |

Names do not repeat `draw_view_context`: the opaque first argument already
identifies a GED view.  Renderer names do not appear in renderer-neutral
headers.  Database traversal, publication, LoD realization, framebuffer, and
direct controller operations are private services and are not stable libged
API.

## Ownership and lifetime

`ged_view_context` and `ged_view_set` are incomplete C types.  A view is owned
by its creating application unless it is registered with
`ged_view_context_owned_add`; a GED-owned view is destroyed by `ged_destroy`.
The view does not expose or transfer ownership of its libbv or renderer
implementation.  Its display endpoint is the only public renderer attachment
boundary and follows the ownership rules in
`obol_endpoint_lifecycle.md`.

Scene, feature, and polygon references are copyable values.  They do not retain
their owner.  Record pointers and strings returned in callbacks are borrowed
and valid only for the callback, or until the documented scene revision
changes.  Functions returning an allocated result or array identify the
matching free operation in their declaration.  A command-result scene owns its
staged records until commit or abort; either operation consumes that staging
object.

## Threading

View lifecycle, endpoint, feature, annotation, polygon, selection, pick, and
result-scene calls run on the view owner's thread.  Value-reference resolution
asserts this rule in assertion builds and rejects cross-thread use otherwise.
GED scene transactions and observer registration are likewise serialized by
the GED owner.  Worker services may compute source geometry, but publication
back to a view is marshalled to the owner thread.  No installed call grants
concurrent mutation merely because its handle is opaque.

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
in that manifest.

