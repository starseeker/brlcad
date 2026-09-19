# Obol drawing: session handoff

Updated 2026-09-19. Start with the [simplification guide](obol_simplification_guide.md)
for direction, [active debt](libbobol_active_debt.md) for tasks, and
[production readiness](obol_production_readiness.md) for acceptance. This is a
replaceable checkout snapshot. Previous details are preserved in the
[September 19 handoff archive](obol_20260919_handoff_history.md).

## Current work

The user has authorized source/documentation rework toward a smaller supported
mutation surface, and then prioritized rebasing on current `origin/main`.
Squashing branch history is authorized if useful. Preserve the full old branch
and dirty work before integration. Review upstream behavior and port it through
the new drawing owners; resolving text conflicts alone is insufficient.

The complete pre-integration state is preserved at `123010cebd` on
`brlobol7-before-main-20260919`. Integration onto fetched main `0d745fca35`
uses one combined drawing commit from the common ancestor `6c94ae71d2`.
There are 553 upstream commits and 261 files changed by both sides.
[The integration record](obol_main_integration.md) maps upstream behavior to
its new owner and records measured results and remaining acceptance checks.
The original branch now uses the integrated baseline; the safety branch keeps
all old history and working changes. No remote changes are made.

The libraries, plugins, MGED, qged and launcher build from the original
checkout in `.build-main`; repository/license checks and default loader
resolution pass. The 34 selected behavioral rows pass after separating Qt's
software/offscreen and window-system GL test paths; both additional Xvfb GL
rows also pass. Targeted annotation,
color, framebuffer, geometry, interpreter, MGED and qged behavioral tests pass.
Launcher resize and pointer quit were exercised and visually inspected.
Both publication sweeps pass in isolated runs; an earlier concurrent run's
scene-light assertion remains an intermittent CONTROL-01 witness. All seven
upstream annotation image comparisons still differ: no controls were changed.
Full visual parity remains open under QUALITY-01/EDIT-01.

Use `.build-main` and the preserved `.build-main-deps` prefix for the rebased
checkout. Keep the older `.build` unchanged as historical evidence. Current
logs, images, the upstream/conflict inventory and exact dependency manifest
are under `.build-main/obol-qualification/20260919-main-integration`.
The next work is to establish BUILD-01's repeatable native-source baseline,
then complete the bounded cold drawing workflow. Annotation visual differences
and the intermittent publication result have concrete retained reproductions.

## Transition prepared before integration

- The unsealed direct viewport field watcher is removed, including its private
  accepted-state snapshot, 13 auditors and extra source reference. Existing
  host/stream APIs and prepared viewport publication remain the live path.
- Raw fields configure detached viewports; explicit rebuild publishes their
  geometry. The API/header contract now distinguishes this from live host
  mutation and documents postcommit observer failure.
- Ordinary framebuffer hide/show, pan/zoom, resize, close/reopen tests check
  independent quad/texture coordinates and exact pixel bytes. The texture
  rectangle helper now keeps image pan/zoom separate from screen fitting.
- The raw-field experiment's tests are removed; preceding allocation-failure
  and lifetime regression coverage is retained with the original source-ref
  expectations. The rebased image-display/window-host tests and isolated
  publication sweeps now pass; the intermittent sweep result remains open.
- Guide, current backlog and readiness are compact current authorities; old
  narratives and the original text rewrite proposal are archived. Completion
  percentages are retired in favor of explicit S0--S6 acceptance criteria.

The first transition build exposed stale Obol and zlib staging in the old build
directory. The separate integration build avoids that stale state. Historical
attempts remain under `.build/obol-qualification/20260919-simplification-transition`.

## Dependencies and evidence limits

Native Obol is the sibling checkout, with local renderer/publication/notification
repairs; matching OSMesa is built there. Both dirty source state and runtime
libraries matter. Do not rely on branch names or old binaries alone. Local
installation must use an explicit build-local prefix: the Obol build's default
install prefix points to the external bundle. Do not replace libraries while
an application is using them. Recheck loader resolution after configuration.
BUILD-01 now treats coherent repeatable staging as an early prerequisite.

The last sealed focused repair is CXX-SOURCE-089 endpoint destruction, under
`.build/obol-qualification/20260913-endpoint-destruction`. Later raw-field work
was never sealed. September 19's assessment ran 14 selected existing CTest rows
and independently reproduced the raw visibility round-trip defect; these were
pre-transition binaries, not current qualification. Its evidence is under
`.build/obol-qualification/20260919-assessment`.

Old `/tmp` GUI reports and warm caches cited in historical notes are unavailable.
Required current qualification includes native graphics, Lucy failing quality
conditions, bounded BREP, shared-cache/large-scale resource behavior, and full
client/edit/platform rows. Software System GL on Xvfb is not native GPU evidence.
The [GUI runner instructions](../../src/qged/tests/README.md) own replay options;
use fresh persistent artifact directories and serialize heavy qualification jobs.
