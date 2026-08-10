# Direct libged/libBObol drawing architecture

This is the entry point for the current drawing architecture.  Historical
migration and realization notes are evidence records; they do not define a
second compatibility design.

The drawing-aware libged commands are C++ and call libBObol's typed services
directly.  libged supplies command parsing, GED transaction policy, and small
non-owning value conversions.  It does not own a C++ draw session or mirror a
controller registry.  C hosts use the renderer-neutral GED C surface, which
reaches the same state owners.

The display endpoint is the sole public attachment boundary.  One GED owns the
shared database scene, and every hosted view owns or borrows one endpoint with
one view controller.  Endpoint replacement recomposes the existing shared
scene; it never transfers or replays scene ownership.  See
`obol_endpoint_lifecycle.md` for the ownership diagram and teardown order.

libBObol owns scene sources, compact instances, realization, view-local LoD,
features, selection/query state, framebuffer composition, and rendering.
Value handles include owner identity, object identity, and generation or
revision, and resolution rejects stale values before accessing storage.  The
stable/advanced header tiers, ABI boundary, threading rules, naming glossary,
and result conventions are defined by `libbobol_api_contract.md` and the
machine-readable `include/BObol/api_tiers.cmake` manifest.

Supported hosts and worker/owner-thread boundaries are recorded in
`libbobol_platform_threading_matrix.md`.  The deliberately stable GED draw
surface is documented in `ged_draw_api.md` and enforced by
`src/libged/ged_draw_api.symbols`.  Current enabling work is limited to the
items in `libbobol_active_debt.md`.

The older `brl_obol_*` and `obol_*_coverage` notes are archived milestone and
behavior evidence.  When they describe a retired renderer or migration step,
the contracts above take precedence.

The progressive drawing safety, liveness, budget, and complexity properties
are specified in `libbobol_lod_state_contract.md`.  That state contract is the
review boundary for LoD scheduler changes; renderer and GUI tests map back to
its explicit proof obligations.
