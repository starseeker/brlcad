# Obol LoD visual-quality qualification

Last reviewed: 2026-09-01

This document defines how BRL-CAD qualifies the image produced by managed LoD
against a full-detail control.  It complements the state and performance
contracts: a view can be live, bounded, and terminal while still choosing a
visually poor allocation.  Conversely, a constrained terminal view is not
expected to contain every source triangle when those triangles cannot affect
the current image.

## Comparison method

`src/qged/tests/qged_lod_quality_matrix.sh` renders full-detail and managed-LoD
views with the same qged binary, renderer, draw mode, lighting profile, canvas,
device-pixel ratio, camera, and event sequence.  Each checkpoint first proves
that the camera and physical canvas match.  Each process also verifies that the
database content is unchanged while it runs, and the pair is rejected unless
both policies used identical content.  The runner hashes each distinct input
once for provenance, records device, inode, size, mtime, and ctime around every
process, and rehashes only when that metadata changes.  This avoids repeatedly
reading a multi-gigabyte input merely to prove identity while still handling
BRL-CAD's harmless rewrite of an identical writable database header.  The
runner records the fitted
camera for each orientation independently and replays its center and size in
BRL-CAD base units (`units mm`); database-local units must never be applied a
second time to report values.  A reusable control without explicit per-view
camera fields is rejected.  The geometry comparison crops the documented
progress-HUD bands, then records:

- libicv `imgdiff -A` structural similarity (SSIM) and perceptual-hash Hamming
  distance;
- the fraction of pixels differing by more than two percent;
- exact foreground-mask disagreement and disagreement after allowing a
  one-physical-pixel boundary displacement, both relative to the union of the
  silhouettes;
- presented faces, render cost and budget, maximum normalized visual error of
  the exact presented cuts, prominent-floor debt, point threshold, structural
  boxes, exact-work status, render time, and the typed performance/memory
  constraint witnesses.  The allocator candidate's maximum normalized error
  is retained in a separate column.  A pose-only view change may preserve a
  richer resident cut than the allocator's currently affordable candidate;
  image acceptance measures the framebuffer, not that unused coarser plan.

The original uncropped images and both foreground masks remain with the test
artifact.  SSIM is the primary whole-image signal, but it is never a topology
oracle: a narrow missing blade, cable, wheel, or skin panel can be important
despite changing few pixels.  Exact silhouette mismatch intentionally flags
subpixel rasterization changes; the one-pixel-tolerant result distinguishes
those from missing geometry.  Qualification therefore combines image metrics,
both silhouette measures, named-feature crops for realistic models, and the
runtime per-occurrence quality-floor diagnostics.  PHASH is useful for gross
regression triage but is too coarse to be a release threshold by itself.

Full-detail controls may require much more time and memory than managed LoD.
They must not be run unbounded when the control is likely to exhaust the
machine.  For 50k/150k and multi-gigabyte scenes, use representative full-detail
subsets and named-feature crops, source-semantic checks, and a bounded LoD
whole-scene run.  An OOM-prone whole-scene control is not useful evidence.

## Current matched evidence

The following shaded measurements use 1100 by 800 qged windows and the current
2026-09-01 binaries.  The ranges cover `ae 90 0`, `ae 35 25`, a two-times close
zoom, and the exact return view.

| Model/backend | LoD faces | SSIM range | Silhouette disagreement | Interpretation |
|---|---:|---:|---:|---|
| Generic Twin/System GL | 134k--135k | 0.999929--0.999982 | 0 | Direct/safe-scene drawing is effectively indistinguishable from full detail. |
| Generic Twin/OSMesa | 134k--135k | 0.999896--0.999994 | at most 0.0037% exact; 0 at 1 px | Software rendering reaches the same full-mesh endpoint; the small residual is raster variation rather than proxy geometry. |
| Hubble/System GL | 342k--708k | 0.973634--0.987800 | 0.34%--1.05% exact; at most 0.031% at 1 px | The current matched run preserves the body, solar arrays, appendages, and formerly truncated cable/tail-scale features.  The close-view SSIM difference is an expected PoP shading advisory rather than missing topology. |
| Hubble/OSMesa | 342k--708k | 0.970219--0.986824 | 0.32%--0.88% exact; at most 0.003% at 1 px | The current software run reaches the same pixel target with no typed resource constraint or prominent-floor violation; renderer time is higher, but topology remains within the one-pixel contract. |
| Lucy/System GL, certified warm | 4.63M--6.12M | 0.984892--0.994758 | 0.31%--0.53% exact; 0 at 1 px | The post-selection-liveness hardware run confirms that a 28M-face source retains high silhouette fidelity while shading and fine detail vary with the performance and memory budgets. |
| Lucy/OSMesa, certified warm | 1.28M--2.10M | 0.971613--0.992283 | 0.41%--1.10% exact; 0 at 1 px | Static-start capacity recovery preserves a rich retained image.  The close zoom records explicit performance and memory constraints rather than stopping at an unowned coarse cut; the exact return population may differ while remaining silhouette-equivalent at one pixel. |
| NIST BREP/System GL | 38.6k--53.2k | 0.992120--0.996329 | 0.09%--0.58% exact; 0 at 1 px | No missing silhouette was observed; the close view uses a memory-constrained pixel tessellation/cut rather than all triangles. |
| NIST BREP/OSMesa | 38.6k--53.2k | 0.992623--0.998831 | 0.09%--0.22% exact; 0 at 1 px | Cross-renderer geometry behavior agrees within the expected raster and lighting differences. |
| Generic Twin wire/edges, both backends | 134k--135k | 0.998978--0.999986 | 0.004%--0.28% exact; 0--0.017% at 1 px | Wireframe and shaded-with-edges reach the same complete mesh endpoint.  Thin-line rasterization makes the exact mask more sensitive than shaded silhouettes. |
| Hubble wireframe, both backends | 342k--708k | 0.940862--0.979983 | 2.24%--4.72% exact; below 0.008% at 1 px | The visible wire topology is retained; almost all exact disagreement is a one-pixel line displacement. |
| Hubble shaded-with-edges, both backends | 179k--708k | 0.908252--0.983435 | 6.17%--17.55% exact; 0.005%--0.65% at 1 px | System GL reaches the pixel target.  OSMesa records responsiveness limits and lower triangle cuts, but visual inspection retains the major panels, appendages, and cable rather than dropping components. |
| Heterogeneous 5k/System GL | 1.56M--2.26M | 0.978818--0.990477 | 0.16%--0.49% exact; 0.004%--0.031% at 1 px | The allocator preserves the recognizable mixed-size scene while drawing 31%--46% of its 4.96M control faces.  All 5,000 occurrences remain represented without structural boxes or terminal proxies. |
| Heterogeneous 5k/OSMesa, capacity-dependent | 988k--1.67M | 0.956942--0.981945 | 0.32%--0.80% exact; 0.008%--0.050% at 1 px | The current run is explicitly responsiveness-limited where needed at 213--398 ms per measured frame.  It retains every occurrence with zero structural boxes, terminal proxies, or prominent-floor violations; host throughput changes shading detail but not topology. |
| Two transformed Lucys/System GL, managed | 6.37M--12.96M | control unsafe | inspected, no holes or missing instance | Both large occurrences refine independently.  A typed performance/memory limit leaves a 1.50 normalized-error close endpoint without starving either instance. |
| Two transformed Lucys/OSMesa, managed | 1.72M--2.57M | control unsafe | inspected, no holes or missing instance | Software rendering assigns cuts 22 and 23 according to their slightly different projected importance.  Both remain solid and recognizable with zero boxes, proxies, or prominent-floor debt. |
| Heterogeneous 50k, initial managed view | System 9.26M; OSMesa 701k | control unsafe | envelope and large-feature inspection | System GL meshes 26,994 occurrences and aggregates 23,006 subpixel parts.  OSMesa meshes all 838 prominent occurrences and aggregates the other 49,162; both are terminal, coherent, box-free, and proxy-free. |
| Heterogeneous 50k wire, initial managed view | System 1.88M faces/5.65M segments; OSMesa 461k faces/1.40M segments | control unsafe | inspected, complete occurrence coverage | Both backends preserve the crowded multicolor envelope and every classified prominent occurrence with zero boxes or terminal proxies.  They are explicitly responsiveness constrained at 105/145 ms measured frames and may stop near 2--3 pixels of projected error.  Wire throughput is a separate production optimization gate, not a shaded-quality regression. |
| Heterogeneous 150k, initial managed view | System 4.51M; OSMesa 1.17M | control unsafe | envelope and large-feature inspection | System GL meshes 62,752 occurrences and aggregates 87,248 below a two-pixel classifier.  OSMesa meshes all 1,865 current prominence candidates and aggregates 148,135 below 32 pixels.  Both retain a one-pixel maximum projected-error certificate, preserve the scene envelope and largest forms, and terminate box/proxy-free.  Software fine detail is intentionally sparse under its measured renderer constraint. |
| Big Boy partial BoT/System GL | 949k--1.60M | 0.968011--0.982675 | 0.62%--4.37% | The partial real-world train conversion retains 14%--23% of its instance-expanded faces.  Named-feature and one-pixel-tolerant comparisons confirm that the high exact side-view disagreement is mostly matching wheel/rod boundaries, not missing prominent components. |
| Big Boy partial BoT/OSMesa | 949k--1.60M | 0.946403--0.973746 | 0.68%--4.02% exact; 0--0.65% at 1 px | The responsiveness-limited software path preserves the wheels, running gear, boiler, cab, and roof without boxes or terminal proxies.  The cold-preview run may begin at 497k faces, but the current certified-warm comparison reaches the same cuts as System GL. |

### Final current-tree spot audit

The following runs were repeated after the terminal background-result handoff
repair.  They are the shortest authoritative current-tree set for checking
whether a later change has altered the quality regimes above:

- `/tmp/qged-lod-quality-basic-{system,osmesa}-20260901` passes all sixteen
  Generic Twin and Hubble control/LoD processes.  Generic Twin SSIM is
  0.999779--0.999987 on System GL and 0.999861--0.999980 on OSMesa; its
  one-pixel silhouette disagreement is zero.  Hubble SSIM is
  0.979149--0.987799 and its one-pixel disagreement is below 0.009 percent.
- `/tmp/qged-lod-quality-nist-system-20260901` presents 178,561 managed BREP
  faces at 0.981867 SSIM, 0.9785 normalized error, and zero one-pixel
  silhouette disagreement.  Ordinary LoD-off BREP tessellation is a fixed-
  tolerance control, not necessarily the finest possible adaptive reference.
- `/tmp/qged-lod-quality-bigboy-system-20260901` presents 988,654 partial-train
  BoT faces at 0.980164 SSIM and 0.67 percent one-pixel silhouette
  disagreement.  Wheels, rods, boiler, cab, and roof remain present.
- `/tmp/qged-lod-quality-lucy-{system,osmesa}-20260901` and
  `/tmp/qged-lod-quality-lucy-terminal-fix-cold-20260901` verify coherent
  single-large-mesh endpoints without boxes or proxies.  System GL retains
  4.63M faces; the cold software run retains 2.10M under its measured
  responsiveness limit.  A whole 28M-face LoD-off control exceeded both the
  6.5 GiB and 9 GiB bounded process limits and is not a safe routine oracle on
  this host.
- `/tmp/qged-lod-quality-multi-lucy-{system,osmesa}-final-20260901` presents
  both transformed Lucy instances without starving either: 9.26M faces on
  System GL and 2.57M on OSMesa, with zero boxes, proxies, or prominent-floor
  violations.
- `/tmp/qged-lod-quality-{50k,150k}-current-20260901` covers the managed shaded
  high-cardinality endpoints.  The 50k pair terminates in 36/19 seconds with
  9.26M/701k faces; the 150k pair terminates in 76/66 seconds with
  4.51M/1.17M faces.  Every visible occurrence is either meshed or represented
  by its subpixel aggregate, and all current prominence candidates are meshed.
- `/tmp/qged-lod-quality-50k-wire-current-20260901` is the current wire stress
  result.  It terminates without cycling or fallback geometry, but its
  105--145 ms frame times and 2.15--2.60 pixel maximum error identify wire
  batching/raster throughput as a remaining optimization target.

These runs also prove that a queued result belonging to already-useful
background cache work cannot revoke terminal readiness under the same semantic
revision.  Applying a representation-changing result still advances the
revision and must satisfy the ordinary foreground contract before the new
endpoint becomes terminal.

A focused current-binary NIST PMI7-10 OSMesa rerun after the BREP terminal-
fallback repair presents 38,629 managed faces versus 46,249 control faces.  It
records SSIM 0.998831, 0.5299 percent changed pixels, 0.2155 percent exact
silhouette disagreement, and zero disagreement with a one-pixel tolerance.
It terminates with exact presented work, no resource constraint, and no
structural box or proxy.  This is the intended pixel-target result: fewer
triangles than the LoD-off control, but no missing image-scale topology.

The matching true-cold GUI timeline contains four source work items.  It moves
from unknown preparation work to 0/4, 1/4, and 3/4 while completed meshes
replace four temporary boxes, then terminates with four meshes, 159,894 faces,
and no boxes.  The public phase remains `PREPARING` throughout that source
work, and its fraction advances from 40 to 49 to 66 percent rather than
appearing to restart discovery after provisional autoview.  These intermediate
boxes are truthful unavailable-source coverage, not PoP cuts or terminal
proxies.

The matching reports and images are in
`/tmp/qged-lod-quality-current-{simple-large,hubble-system,hubble-osmesa,nist}`,
`/tmp/qged-lod-quality-static-start-regression-20260831`,
`/tmp/qged-lod-quality-post-selection-lucy-hw-system-20260831`,
`/tmp/qged-lod-quality-current-unique5k-hw-system-20260831`,
`/tmp/qged-lod-quality-current-unique5k-osmesa-20260831`,
`/tmp/qged-lod-quality-current-bigboy-hw-system-20260831`, and
`/tmp/qged-lod-quality-bigboy-osmesa-two-sided-fix-20260831`.  The current
dual-backend Big Boy rerun is in
`/tmp/qged-lod-quality-terminal-fix-bigboy-20260831`.  The repaired
certified-warm Lucy OSMesa comparison is in
`/tmp/qged-lod-quality-lucy-osmesa-static-start-20260831`.
The post-terminal-composition reruns are in
`/tmp/qged-lod-quality-terminal-fix-{simple-hubble,lucy-system,lucy-osmesa,5k,5k-single,5k-system,multi-lucy-pair,50k-initial,150k-initial}-20260831`.
The `5k-single` case fixes `LP_NUM_THREADS=1` to reduce available throughput;
it is a deliberate capacity perturbation, not a preferred configuration.
The cold-preview-to-spatial-generation repair is qualified by
`/tmp/qged-lod-quality-bigboy-preview-transition-osmesa-20260831`: it reaches
four terminal views in 17 seconds with zero terminal occurrence or renderer-
page failures.  Its lower-capacity 497k--949k-face allocation broadens the
OSMesa range above and records side-view feature SSIM of 0.791 for the front
wheels, 0.781 for the running gear, and 0.879 for the boiler/cab.
The current wire/edge matrix is in
`/tmp/qged-lod-quality-wire-final-20260831`; its one failed pre-fix compact
OSMesa replay is replaced by the passing exact rerun at
`/tmp/qged-lod-quality-generic-m4-osmesa-debounce-fix-20260831`.
The 2026-09-01 post-contract wire/edge matrix is in
`/tmp/qged-lod-quality-post-contract-wire-edge-20260901`.  Its Hubble OSMesa
mode-4 image completed successfully but the strict trace exposed a transient
producerless importance census; the contract-clean replacement at
`/tmp/qged-contract-fixed2-hubble-osmesa-m4-20260901` terminates ready in 26.3
seconds with the same box-free visual class and zero violation across all 130
recorded transitions.  The final current-binary rerun at
`/tmp/qged-lod-quality-wire-edge-contract-final-20260901` passes all eight LoD
rows in 4--22 seconds.  Generic Twin satisfies every safe-direct SSIM target.
Hubble remains box/proxy-free with no prominent-floor debt; direct inspection
of the advisory OSMesa mode-4 pairs confirms the body, panels, appendages, and
cable remain present at the responsiveness-constrained endpoint.
These paths are development evidence, not permanent release artifacts; a
release run must archive its reports and images outside `/tmp`.

The 2026-09-01 post-visibility-contract qualification is the current coherent
complexity ladder:

- `/tmp/qged-lod-quality-post-visibility-simple-hubble-20260901` contains
  byte-matched Generic Twin and Hubble full-detail/LoD pairs on both renderers;
- `/tmp/qged-lod-quality-post-visibility-nist-20260901` covers adaptive shaded
  BREP presentation;
- `/tmp/qged-lod-quality-post-visibility-lucy-20260901` compares the 28M-face
  Lucy source, and records 103/143-second full-detail controls versus
  12/32-second managed System-GL/OSMesa runs;
- `/tmp/qged-lod-quality-post-visibility-unique5k-20260901` covers 5,000
  distinct, mixed-size meshes;
- `/tmp/qged-lod-quality-post-visibility-bigboy-20260901` contains the partial
  real-world train comparison and named wheel, running-gear, and boiler/cab
  crops; and
- `/tmp/qged-lod-quality-post-visibility-{multi-lucy-pair,50k-initial,
  150k-initial}-20260901` contains the managed-only tiers for which a complete
  control is unsafe or misleading.

The post-capacity-contract rerun supersedes those terminal measurements while
retaining their controls where the database digest, camera, renderer, canvas,
and lighting are identical:

- `/tmp/qged-lod-quality-post-contract-{generic,hubble,nist,unique5k,
  bigboy}-20260901` contains fresh matched System-GL and OSMesa runs;
- `/tmp/qged-lod-quality-post-contract-lucy-cold-20260901` uses a fresh LoD
  cache against the byte-matched full-detail controls;
- `/tmp/qged-lod-quality-post-contract-multi-lucy-pair-20260901` covers two
  independently allocated Lucy-scale occurrences; and
- `/tmp/qged-lod-quality-post-contract-{50k-initial,150k-initial}-20260901`
  contains the current managed-only high-cardinality endpoints.

Every one of the 48 matched checkpoints passes its executable quality and
presentation contract.  Generic Twin remains indistinguishable from the
full-detail control, Hubble and NIST retain their one-pixel topology, Lucy
remains crack-free through close zoom and return, and the 5k and partial Big
Boy comparisons preserve the mixed-size and named prominent features.  Both
managed high-cardinality searches reach their exact finite 40/40 capacity
rank and terminate with owner, obligation, violation, structural-box,
terminal-proxy, and prominent-floor counts all zero.

The current-binary Lucy audit at
`/tmp/qged-lod-quality-current-lucy-20260901` and
`/tmp/qged-lod-quality-current-lucy-osmesa-presented-20260901` records
103/140-second full-detail controls versus 36/25-second managed System-GL and
warm OSMesa rows.  Across both backends, exact silhouette disagreement is
0.31--1.10 percent and becomes zero at a one-pixel tolerance; inspection shows
no wing, drapery, hand, base, or spatial-page hole.  The OSMesa `ae35` view
also exercises pose continuity: cut 24 remains presented at 0.712 normalized
error even though the new allocator candidate is the coarser cut 23 at 1.103.
Qualification now records both values and correctly applies the presented
value to the captured-image contract.  The close and return views remain
explicitly constrained when their presented error exceeds one.

`/tmp/qged-lod-quality-hubble-contract-final-20260901` supersedes the Hubble
rows after the terminal-handoff audit.  Modes 0 and 4 pass on both renderers;
all four scripted views in every run have exact presentation evidence, all
2,030 payloads satisfied, and zero structural boxes, terminal proxies, or
prominent-floor violations.  The external trace also proves that the temporary
one-to-64-pixel point-classifier change and its one-pixel restoration occupy
distinct capacity revisions rather than publishing competing allocation plans
under one revision.

All 48 matched checkpoints satisfy the geometric/presentation contract, and
all managed checkpoints terminate box-free and proxy-free with zero prominent
floor violations.  Generic Twin remains a strict safe-direct case with SSIM
above 0.999 and zero silhouette disagreement.  For the larger shaded cases,
one-pixel-tolerant silhouette disagreement is the topology gate; SSIM is also
reported, but simplified PoP normals can change illumination without moving a
silhouette.  The lowest fresh matched SSIM is 0.946358 for the close OSMesa
Big Boy view, whose tolerant silhouette disagreement is zero and whose named
wheels/running gear remain present.  The 5k OSMesa close view similarly scores
0.956942 SSIM but only 0.0078 percent tolerant silhouette disagreement.  These
are photometric review signals, not missing-geometry evidence.

The current 150k System-GL/OSMesa runs finish in 76/66 seconds.  The current
harness hashes the 3.4 GB database once per matrix and uses metadata to avoid
re-reading it around every process; an identity-changing process still forces
a validating rehash.  The runs peak at 8.58/6.89 GB process RSS.  Their
controller-accounted resident mesh sets are only 297/29 MB
respectively; the much larger RSS includes reclaimable database and cache
mappings and must not be reported as LoD mesh residency.  The hardware
endpoint draws 4.51M faces in a 30 ms measured frame.  The software endpoint
draws 1.17M faces in 267 ms with an explicit responsiveness witness.  Both
certify at most one pixel of projected error.

The lower hardware face count is an allocation improvement, not evidence of a
missing population.  Against the preceding 20.49M-face terminal image at the
identical camera, the new image has SSIM 0.972266 and only 15 one-pixel-
tolerant silhouette samples differ, 0.037 percent of the silhouette union.
The OSMesa endpoint has SSIM 0.948675 and 0.574 percent tolerant silhouette
disagreement against its preceding richer endpoint, remaining comfortably
inside the provisional 2.5 percent topology threshold.  Visual inspection
retains the same envelope and large forms.  Thus the production expectation
at 150k is a coherent useful preview followed by roughly minute-scale warm-
cache terminal convergence, not small-model latency, exhaustive mesh loading,
or hardware-equivalent software shading density.

`qged_lod_quality_matrix.sh` now writes `matched_assessment.csv` and
`managed_assessment.csv`.  The former makes safe-direct image identity,
pixel-error/topology, exact presentation, and typed constraint evidence
executable while retaining shaded SSIM as an advisory target outside the safe
class.  The latter validates managed telemetry but always requires visual
review because a run without a safe full-detail control cannot prove image
quality.

The post-structural-repair OSMesa hidden-line qualification is in
`/tmp/qged-generic-osmesa-m4-{cold-r12,warm-r13}`.  Both the isolated-cold and
same-cache warm runs pass all four views in 11--13 seconds.  Every nonempty
checkpoint contains all expected 673, 680, or 709 terminal mesh occurrences,
uses a one-pixel point threshold, and has zero structural boxes, terminal
proxies, pending exact work, or control violations.  Libicv reports SSIM
0.998978--0.999833 and silhouette disagreement 0.0997%--0.2804%.  The two
runs currently select different retained renderer representations (flat
batching versus ordinary retained buffers), so their render-cost currencies
are 2.56M--2.58M and 1.87M--1.88M units despite matching faces and images.
That representation-dependent accounting is intentional: it measures work
submitted by the backend and must not be read as a topology count or compared
across renderer routes without the paired duration.  Allocator decisions use
their separate canonical retained-scene currency.

A true-cold Lucy run reached a useful spatial preview substantially before its
28M-face hierarchy finished.  Its first comparison was consequently much
coarser than the certified-warm terminal result.  That is intended progressive
behavior, but the HUD must distinguish **view usable** from **asset preparation
complete**.  Cold-preview quality and warm terminal quality are separate
measurements and must not be combined into one completion percentage.

OSMesa quality is expected to be lower, but not arbitrary.  Earlier 5k runs
made under different host load selected very different measured capacities:
one presented 262k--370k faces while a heavily loaded run fell as low as 8k.
Those runs exposed controller defects and are not the current quality
baseline.  After the terminal-composition repair, normal and deliberately
single-threaded runs retain 421k--1.67M faces with SSIM
0.942882--0.981945, zero terminal proxies, and no prominent-floor violations.
The remaining range is legitimate adaptation to measured throughput, not a
stable visual-quality guarantee.
Release reports must therefore record system load, renderer throughput, the
typed constraint, and prominent-floor debt rather than publishing one OSMesa
SSIM number as a backend property.  A production operator who selects a lower
minimum FPS or a larger static-frame allowance should obtain a richer retained
image without changing topology or restarting discovery.

The current warm System-GL 50k interaction run at
`/tmp/qged-production-current-50k-batched-20260830` terminates with 34,760
mesh occurrences, 15,240 subpixel occurrences, 14.49M triangles, and no
structural boxes or terminal proxies.  A formerly hard-rejected relaxation
from the four-pixel projection bucket now preloads its 10,849-occurrence
frontier in deterministic batches of 8,192 and 2,657 while the preceding
coherent point image remains public.  The exact remaining rank decreases
10,849 to 2,657 to zero, permitting the view to reach the one-pixel point
threshold instead of stopping at 16 pixels and 6.32M triangles.  The complete
draw, zoom, rotation, selection, erase, and redraw workflow takes 83 seconds;
its first pose is terminal in 39 seconds.  This is a quality improvement, not
an unbounded first wave: the complete candidate is first projected against the
static frame deadline, every provider batch remains bounded, and failure to
reduce the exact rank rejects the private candidate.

The post-terminal-composition initial-view replay retains 9.26M faces through
26,994 meshed occurrences and aggregates the other 23,006 subpixel
occurrences.  It terminates in 46 seconds with zero boxes, terminal proxies,
or quality-floor violations.  OSMesa terminates in 24 seconds with 617k faces:
all 838 prominent occurrences are meshed and the other 49,162 are aggregated.
Its 149 ms terminal frame has an explicit responsiveness constraint.  This is
a recognizable whole-model result, but its component-level shading is visibly
coarser than System GL and should be presented to users as constrained quality.

The corresponding 150k System-GL run at
`/tmp/qged-production-current-150k-batched-20260830` uses the same 6.6 GB
certified cache without copying it.  It provides a whole-model structural
preview while discovery is active, reaches a useful shaded scene before its
pixel target, and eventually presents 9.54M triangles through 91,710 mesh
payloads plus 61,771 aggregated subpixel occurrences.  It retains zero
structural boxes, terminal proxies, prominent-floor violations, or visual-
importance debt while keeping resident mesh data to about 642 MB.  The first
pose takes 179 seconds to reach its pixel-target endpoint and the complete GUI
workflow takes 229 seconds.  That latency is not a small-scene expectation:
it is the current cost of reading, binding, and allocating tens of thousands
of distinct cached meshes.  The early preview is therefore the usability
contract at this tier; final pixel-target convergence may take minutes, and
further throughput work must reduce that time without weakening the bounded
memory, coherent-frame, or visual-priority contracts.  A full-detail 150k
control is deliberately not attempted on this 15 GB test machine.

The post-terminal-composition initial-view replay terminates in 220 seconds on
System GL with 20.49M faces, 78,341 meshed occurrences, 71,659 subpixel
aggregates, and 688 MB of resident mesh data.  Its measured 51 ms frame and
687 MB resident population carry explicit performance and memory witnesses.
OSMesa reaches its more aggressively constrained endpoint in 51 seconds with
1.17M faces, all 1,865 prominent occurrences meshed, and 148,135 subpixel
aggregates.  Both endpoints are ownerless, box-free, proxy-free, and preserve
the same large-feature envelope.  The large System-GL latency remains a
production throughput issue even though terminal correctness now holds.

A current 199 Hz `perf` sample of the warm 50k System-GL initial view is in
`/tmp/qged-lod-50k-terminal-fix-20260831.perf`.  It captures 13,139 samples
over a successful 30.5-second replay.  The former 13.6% sparse scene-mutation
snapshot hotspot does not appear; the upstream Obol data-structure work has
retired it.  No current symbol exceeds 3% exclusive cost.  The leading work is
distributed across memory movement, instanced rendering, matrix and normal
transforms, allocation, mesh validation, and parallel cache decoding.  The
remaining high-cardinality latency is therefore accumulated realization and
publication throughput.  Do not optimize the obsolete rollback snapshot or a
single sub-percent helper without new evidence; the next performance pass
needs phase-level timing and batching/parallel-publication experiments on the
150k path.

## User-visible quality contract

The production contract is based on projected error and resource evidence,
not on a promise to load all geometry:

1. A safe, readily drawable scene goes directly to mesh data.  A terminal
   Generic-Twin-scale view must contain no persistent aggregate or structural
   boxes and should be visually indistinguishable from full detail.
2. An unconstrained view refines until every visible occurrence satisfies its
   projected pixel target.  Different valid PoP cuts or adaptive BREP
   tessellations may be selected after a round trip; equal triangle counts are
   not required if the projected-error and topology contracts remain true.
3. A resource-constrained view presents the best scene-wide allocation proven
   to fit the configured frame and memory budgets.  Prominent, selected, and
   large-screen-error features receive priority; subpixel occurrences may be
   aggregated.  The HUD must identify the limiting resource rather than claim
   unconstrained completion.
4. OSMesa follows the same geometry and state contracts as System GL.  Its
   software renderer may select fewer triangles because measured render cost is
   higher, but renderer choice cannot justify holes, invalid topology, stale
   boxes, or cycling.
5. Structural AABBs are temporary coverage.  A persistent AABB or cached OBB
   is acceptable only as an explicit terminal proxy selected under measured
   scene pressure.  It must retain the occurrence's style, selection,
   transform, picking semantics, and conservative extent.
6. Camera motion starts from retained cuts.  Rotation and translation restore
   the last proven affordable stationary allocation; orthographic zoom changes
   demand incrementally.  Neither transition may restart an unchanged scene at
   structural boxes.
7. A terminal frame is stable for its exact inventory, visibility, view,
   policy, and capacity revisions.  Background cache creation or memory
   compaction may continue, but it cannot repeatedly reopen visible balancing
   or make represented geometry disappear.
8. Renderer presentation controls are part of the exact-frame identity.  A
   completed frame is not current after its progressive ceiling, fractional
   cut, point threshold, or motion-reuse mode changes.  Frame completion may
   itself open capacity, handoff, or demand work and must arm the host pump;
   progress cannot depend on a later pointer event or repaint.
9. A renderer duration and its exact presented render cost are one evidence
   item.  Structural admission, capacity search, and deadline recovery may not
   combine the cost of the current scene with timing measured for a different
   representation or population.  A direct, compact scene may spend the
   400-millisecond prominent one-shot deadline to replace temporary coverage
   with useful meshes; an unbounded scene retains the finite 200-millisecond
   structural-census deadline.
10. A persistent terminal proxy requires a current capacity-rejection witness.
    Exact geometry presented above an older allocation certificate is not an
    ownerless endpoint: capacity reconciliation must recertify it or select a
    bounded successor.
11. The capacity search's protected quality floor is part of its immutable
    numeric domain.  Applying the floor after selecting a lower candidate can
    alternate forever between the numeric candidate and the richer rendered
    population.  Capacity-owned timing allowances are frozen for the same
    reason.
12. Interrupted presentation replay is the exclusive completion of an
    existing transaction; ordinary exact-presentation debt is not.  When a
    capacity candidate needs both allocation and a subsequent exact frame,
    allocation owns the next transition.  The trace checker recomputes that
    owner from the individual exported facts.
13. Selective source-delta scope retires with its completed cursor.  It cannot
    survive a later presentation barrier as ownerless planning state or block
    a complete demand pass.  Runtime and offline refinement checks both reject
    a selective scope without a submission cursor.
14. A valid source tessellation is terminal drawable geometry even when PoP or
    persistent cache acceleration is unavailable.  Such a fallback loses
    view-aware refinement and memory trimming for that asset, but it cannot be
    represented as a PoP rejection or a persistent structural box.

## Provisional qualification targets

These are corpus-building targets, not permanent global constants:

- Safe/direct scenes: SSIM at least 0.999, zero silhouette disagreement, zero
  persistent proxies, and no topology-significant named-feature difference.
- Unconstrained pixel-target scenes: one-pixel-tolerant silhouette
  disagreement below 2.5%, maximum normalized error at most 1.0, exact
  presentation, and no topology-significant named-feature difference.  Shaded
  SSIM of 0.97 remains the photometric target, but is advisory outside the
  safe-direct class: a lower value requires image review and must not override
  stronger topology evidence merely because simplified normals changed
  illumination.
- Constrained scenes: a valid typed constraint witness, no structural fallback
  debt, and no cheaper allocation with a better measured importance/error
  reduction.  Prominent quality-floor debt must be zero whenever the protected
  population fits the hard static-frame deadline.  If it cannot fit, the exact
  count and named features are release evidence, not something hidden behind a
  generic `ready` label.  Whole-image SSIM is reported but is not by itself a
  pass/fail criterion.

Thresholds must be revised using a broader realistic corpus rather than tuned
to Lucy or Hubble.  Required additions are vehicle wheels and blades, mixed
large/tiny components, multi-Lucy and distinct-large-mesh scenes, shaded and
wire modes, and the heterogeneous 50k/150k fixtures.  Named-feature crops need
stable semantic paths so a globally high SSIM cannot hide a conspicuous local
failure.

The practical expectation for a large model is therefore bounded rather than
uniform: every visible occurrence has a conservative representation quickly;
recognizable large features remain preferentially meshed; fine topology
approaches the pixel target as capacity permits; and a constrained terminal
view is retained without cycling.  System GL should normally satisfy the
prominent floor on the current 5k tier.  OSMesa may need coarse triangles and
subpixel aggregates to keep input responsive, but persistent structural boxes,
holes, invalid vertices, missing named features without an explicit quality
constraint, and false `View ready` reports are failures on either backend.

The current allocator should not be retuned from SSIM alone.  Its
recognizability floor is minimax over projected error, followed by a
rate-distortion queue using affected screen span, error reduction, and added
render cost.  The fresh Hubble, Big Boy, multi-Lucy, and 50k/150k evidence
shows that this ordering preserves the tested large and thin features.  The
remaining shaded disagreement is predominantly normal/illumination change at
otherwise valid cuts.  Computing SSIM or scanning triangles in the view hot
path would be both renderer-dependent and too expensive.  If a broader
vehicle corpus demonstrates prominent interior shading failures despite good
silhouettes, the appropriate experiment is fixed-size builder-time per-cut
metadata such as a conservative normal-cone or projected surface-change
bound, incorporated once into the marginal benefit.  Such a cache-format
change requires A/B evidence that it improves named-feature quality without
reducing convergence or 50k/150k throughput; it is not justified by the
current corpus.

## Correctness fixes found by the matrix

The comparison work found control defects rather than allocator tuning
issues.  First, exact presentation reports were keyed only by view identity;
changing the renderer's point threshold or progressive ceiling could leave an
old framebuffer certified as current.  Exactness is now checked against all
effective presentation controls, and every such mutation opens an exact-frame
obligation.  Focused renderer tests prove that an old report is rejected and a
successor report commits the new threshold.

Second, a completed Qt frame can open new timing-driven controller work after
the previous progressive timer has retired.  A level which was already true
does not emit another producer edge, so an intermittent OSMesa Generic Twin
run stopped at 99% until unrelated input arrived.  Frame completion now always
re-evaluates the level-triggered host-work predicate.  The focused
`ObolHostWork` model includes this transition, and three consecutive
untraced OSMesa GUI runs terminate in 8--9 seconds with full meshes and no
surviving qged process.

Third, structural repair used to pair the current scene cost with the last
duration recorded at the same point threshold, even when that duration
belonged to a box-only or otherwise cheaper frame.  Renderer evidence now
stores cost and duration together and an admission decision consumes only an
exact cost match.  Compact direct scenes reserve source-estimated work from a
single scene budget and receive the prominent one-shot deadline; large or
unknown inventories retain equal-share bounded census admission.  Finally,
terminal-quality ordering treats an exact presentation above a stale
certificate as capacity-owned work.  This prevents both premature persistent
OBBs and the observed 672-of-673 ownerless terminal tail.

Fourth, provider teardown cleared the shared host pump when the provider list
became empty even if the provider's last geometry mutation still owed an exact
CAD frame.  Teardown and result draining now call one race-safe synchronization
predicate covering provider, submission, result, publication, residency, and
presentation debt.  The integrated `ObolHostWork` model checks provider
retirement, unsubmitted source revisions, and pump-to-exact-render handoff
over 78,500 generated / 12,902 distinct states to depth 17.  A focused
executable regression reproduces the
old lost wakeup, and three consecutive System-GL Lucy lifecycles settle in
4.42--4.50 seconds at the same 4,628,212-face cut with no exact-frame debt or
control violation.  The matched legacy full-detail images have identical
camera and canvas records; current libicv SSIM is 0.980581--0.991645 and
silhouette disagreement is 0.28%--1.45% across the four views.

Fifth, a private finer-point preload invalidated itself when its own resident
results advanced the source inventory revision.  Hubble's first System-GL
view consequently stopped at a 32-pixel point cutoff with only 57,049 faces,
despite a sub-millisecond frame and ample scene budget.  Candidate-domain
invalidation now distinguishes expected immutable-result publication from
source/coverage, visibility, view, and policy changes.  Shared resident records
may require another occurrence-binding pass, so finalization admits only a
strictly smaller exact remainder.  The updated structural-frontier model proves
that this rank terminates.  The matched replay reaches the one-pixel cutoff
with 341,875 faces, SSIM 0.987807, 1.048% silhouette disagreement, and no boxes,
terminal proxies, constraint witness, or control debt.

Sixth, the same private preload treated the structural first-wave batch size
as a semantic ceiling on the complete candidate.  A 50k frame which could
afford the projected finer population was consequently stopped because 10,849
occurrences exceeded the 8,192-request service quantum.  Candidate pricing and
provider batching are now separate decisions.  The controller prices the
complete finer population, admits deterministic bounded subsets, and requires
an exact strictly decreasing remaining-occurrence rank after every subset.
The public point threshold and capacity certificate remain unchanged until the
rank reaches zero.  `ObolStructuralFrontierOwnership` now models the batch as
well as the rank; TLC exhaustively checked 60,535 generated / 29,266 distinct
states to depth 386 with no safety, deadlock, or liveness failure.

Seventh, exact hierarchy erase/redraw was incorrectly classified as immutable
inventory and capacity invalidation.  That discarded a still-safe renderer
budget, selected a much coarser 150k successor, and could leave newly visible
entries behind an allocation launched before their sparse delta was applied.
Visibility now owns the sixth planning revision.  Its bounded journal is
applied first, followed by one exact framebuffer classification and then one
allocation successor that preserves the capacity certificate; a newer edit
supersedes an unpublished predecessor allocation.  The frame prerequisite is
essential because a restored occurrence may have lost its mesh while hidden,
and the predecessor framebuffer cannot report that structural replacement
work.  The final OSMesa 50k/150k workflows retain 1,206,806/835,318
faces and zero boxes across selection, 16-occurrence erase, redraw, and
deselection.  The 150k erase/redraw settles in 6.05/3.96 seconds rather than
the prior 28.4-second quality collapse.  Both runtime traces pass the
six-domain checker.  `ObolSubmissionPass` now distinguishes census-current from
presentation-current and explores 55 generated / 28 distinct states to depth
8.  The terminal composition separately proves the same ordering before
reallocation over 1,520 generated / 1,187 distinct states to depth 81.

GUI checkpoint waits also consume one already-requested endpoint presentation
before returning `viewReady`.  This does not wait for cache compaction, but it
prevents a framebuffer capture from retaining the preceding balancing label
while controller telemetry already describes the terminal or background-
optimization HUD state.

Eighth, a retained spatial result could be prepared under a valid current
allocation and then arrive after its immutable mesh publication invalidated
the allocator's census record.  Result adoption consequently clamped the
prepared cut back to the older active cut and immediately requested the same
suffix again.  On Lucy OSMesa this looked like high CPU use with no visible
progress.  Provider results now carry a view/policy/cut admission certificate;
the controller accepts a richer result only when either the current allocation
or that exact certificate authorizes it.  Request coalescing cannot rewrite
the certificate's provenance.  The repaired run reaches cuts 11, 14, 17, 19,
and 21/23 in ten service tasks, terminates with no boxes or proxies, and
produces the matched image-quality row above.

Ninth, single-object selection can activate an edit manipulator.  Those
retained style/manipulator changes were using the capacity-relevant endpoint
request, so selection reopened LoD calibration; System GL could additionally
leave the exact presentation behind a coalesced `QOpenGLWidget` update.  The
endpoint now exposes a distinct presentation-only frame request used by ARB,
ellipsoid, sketch, and navigation-gizmo updates.  Direct GL explicitly renders
the sole exact-presentation barrier when it has no other wake.  The
dual-backend selection UI matrix and full Hubble hierarchy selection,
erase/redraw, and deselection workflows now return to idle without a capacity
sample, policy revision, or service task.

A focused full-Hubble point-selection replay after the cold-preview fixes
provides the narrower timing check.  System GL and OSMesa each request one
semantic presentation frame on mouse release, return to IDLE/ready in 101 ms
and 121 ms respectively, preserve the same 481,208-face allocation, and have
no queued or in-flight service work.  Selected frames sampled immediately
after that handoff and again through two seconds of quiet are pixel-identical.
Selection may therefore repaint retained style once; it must not reopen
balancing or refinement.

A later timing-sensitive OSMesa run exposed a separate post-selection/quiet-
handoff liveness hole: the retained importance census remained pending after a
stronger presentation owner retired its bounded submission cursor.  The
unchanged-signature fast path now recreates one ordinary dense demand pass from
that level-triggered obligation.  The repaired full Hubble lifecycle passes in
11 seconds on System GL and 14 seconds on OSMesa; two additional independent
OSMesa repetitions also pass.  At the selected checkpoint every run preserves
the LoD policy revision and returns to IDLE/ready with no controller obligation
or service work.  A fresh matched Generic Twin matrix then completes all four
full-detail/LoD renderer rows in 4--5 seconds.  Its LoD silhouette disagreement
is zero at exact and one-pixel tolerance in every tested view, with SSIM
0.999860--0.999992 and no boxes, terminal proxies, prominent-floor violations,
or foreground control work.

Tenth, an unbracketed zoom could complete its first motion frame, remain in
the intended quiet debounce, and then miss an OSMesa presentation deadline.
The interruption armed a replay after the motion-frame gate was already
cleared; replay priority prevented the quiet reducer from running, so the
compact Generic Twin shaded-with-edges case rendered for two minutes without
settling.  Deadline observation now retires replay whenever the interaction
session is quiet-ready, whether the motion frame completed normally or its
gate expired at the deadline.  The strengthened interaction model checks this
composition, and the exact failed replay now completes in 11 seconds with four
idle, box-free, full-mesh checkpoints.

Eleventh, non-exact PoP cuts correctly disabled back-face culling because
quantization can expose a triangle's reverse side, but the fixed-function
software renderer continued to choose one-sided lighting from the exact source
mesh's oriented flag.  Large surfaces consequently appeared nearly black in
OSMesa even though System GL selected the same cuts and its GLSL shader flipped
the normal for a back-facing fragment.  Flat-occurrence planning now derives
two-sided lighting from the displayed cut's culling decision.  The Obol
renderer regression exercises an oriented source at a non-exact cut on both
fixed-function and GLSL routes.  In the end-to-end Big Boy replay, the black
boiler/cab/roof patches disappear and the four matched-view SSIM range improves
from 0.898631--0.955979 to 0.946401--0.973746 without weakening the silhouette,
box, proxy, or quality-floor contracts.

Twelfth, a cold spatial producer temporarily published a valid whole-prefix
preview before its durable private-page hierarchy completed.  The page
residency query incorrectly treated that unpaged generation as proof that any
requested private page set was drawable.  Final delivery therefore skipped
the representation transition and classified the later missing page geometry
as a terminal failure; depending on task order, Big Boy either stopped near
98 percent or recovered only after a later view change.  Nonempty page demand
now requires an actual spatial generation.  The provider also prices the
first view-local drawable population rather than the hierarchy-wide minimum,
admits an underestimated indivisible minimum as the bounded measured
correction, and treats a legitimately empty page population at a coarse cut as
a successful zero-draw state.  Focused realization/service tests cover all
three boundaries.  `ObolLodColdPreview` now forbids a resident binding before
the spatial generation exists; TLC exhaustively checks 1,902 generated / 918
distinct states to depth 13.  The fresh cold OSMesa run cited above has zero
terminal failures, zero page-preparation failures, zero boxes, and zero
terminal proxies.

Thirteenth, occurrence-local static-quality reconciliation could finish while
leaving its temporary renderer-wide ceiling installed.  A richer local cut
would then restart the same global-ordinal walk, alternating balancing and
refinement near the terminal endpoint.  The handoff now has one typed commit
edge: an ordinary successful allocation releases the ceiling and requests one
deadline-bounded ceiling-free measurement; only a revision-bound protected-
minimum constraint may retain it.  `ObolTerminalConvergenceComposition`
models exact-target renderer preparation, static-to-capacity re-entry,
strictly coarser local reconciliation, and subordinate background compaction
as finite semantic ranks.  TLC checked 1,520 generated / 1,187 distinct states
to depth 81.  Post-fix Generic Twin, Hubble,
Lucy, multi-Lucy, 5k, 50k, and 150k runs all terminate ownerless with no
structural boxes, terminal proxies, or prominent-floor violations.

## Big Boy mesh scope and shaded BREP qualification

`.build/bigboy.g` now contains a partial `big_boy.bot` conversion.  The current
database reports 263 BoT leaves in that hierarchy: 172 nonempty leaves with
4,397,952 faces and 91 zero-vertex, zero-face leaves.  (Repeated references to
empty leaves produce 206 terminal-empty occurrences among 576 hierarchy
occurrences, leaving 370 drawable occurrences.)  The
nonempty portion is a useful real-world test of visual-importance allocation:
compare managed LoD against LoD-off drawing of the same `big_boy.bot`
hierarchy.  It is not a completeness oracle for the original train and the
empty leaves must not be counted as successful geometry realization.

The August 30 System-GL qualification uses the byte-identical control and LoD
pair in `/tmp/qged-lod-quality-bigboy-stable-hash-20260830`; both runs record source
digest `48723c5f229c02bb9c1bc763c4972540c0d6816701418dec2d80fa3465adbb45`.
The database mtime changed during the control despite that identical digest,
which prompted the harness's content-identity check.  Including two 543 MB
digest reads per process, it draws a nonblank LoD-off control in eleven seconds
and the managed view in ten seconds.  Across the standard
oblique, orthographic, close-zoom, and return views, LoD presents
948,928--1,596,578 faces with no structural boxes or terminal proxies.  Against
the 6,932,072 instance-expanded faces in the LoD-off control, this is a
14%--23% retained population.  All four checkpoints satisfy their quarter-pixel
projected-error target without a performance or memory constraint, and the
137--257 prominent visible occurrences have zero protected-floor debt.  SSIM
is 0.968001--0.982677 and silhouette disagreement is 0.62%--4.37%.
Inspection of the worst, side-on mask shows that its 4.37% value is dominated
by one-pixel raster-edge disagreement around matching wheels, rods, boiler,
cab, and rail-bed silhouettes rather than a missing prominent component.  This
is confirmed by its 0.67% one-pixel-tolerant whole-scene disagreement.  Named
side-view crops record exact/tolerant disagreements of 6.66%/1.30% for the
front wheels, 6.00%/1.05% for the running gear, and 1.52%/0.01% for the
boiler/cab.  This is why named-feature inspection remains part of the contract
and exact silhouette disagreement is not used as a lone pass/fail oracle.
These figures are
visibility-allocation evidence only; they do not qualify the missing
conversion results.

The NIST PMI7-10 shaded BREP hierarchy now qualifies the adaptive lifecycle on
both renderers.  The matched System-GL run at
`/tmp/qged-lod-brep-pmi7-10-reclaim-trace15-20260901` and OSMesa run at
`/tmp/qged-lod-brep-pmi7-10-reclaim-osmesa-20260901` use a deliberately reduced
resident limit after the first close view, then zoom out and return.  All five
views present all four targets with no structural boxes, terminal proxies, or
prominent-floor violations.  The close and restored views select the same
220,637-face pixel-target population.  Their matched-control SSIM is 0.993568
on System GL and 0.994099--0.994100 on OSMesa, and their one-pixel-tolerant
silhouette disagreement is zero.  This demonstrates reclamation and
cache-backed restoration without requiring every tessellation triangle to
remain resident.

That result does not qualify the original full Big Boy BREP hierarchy.  The
partial `big_boy.bot` control remains useful for scene-wide visual-importance
testing, but the source BREP still needs complete shaded conversion and the
same constrained close/return/restore comparison before it can serve as the
large real-BREP production gate.

A pre-validation probe of
`Default/68_inch_tire_1/68_inch_tire_1_item/68_inch_tire_1_item.s` records the
reason this gate remains open.  The LoD path publishes one 64-face adaptive
payload with a 0.201-pixel reported projected error, but the LoD-off mode-2
control publishes no shaded CAD payload; its visible ring is the conventional
BREP edge presentation.  Comparing those images would therefore score two
different representation contracts, not LoD fidelity.  The retained evidence
is `/tmp/qged-lod-bigboy-brep-wheel-current-20260901`.

That 64-face payload was not a valid coarse representation.  The legacy CDT
reported constrained-edge failures for five faces, retained the successful
face fragments, and returned aggregate success.  Its arrays contained 368
allocated points but only 32 referenced points and 128 indices.  The
referenced result occupied only the upper plane near `z=65.109`, while the
source tire extends to `z=56.637`.  PoP then classified exactly that incomplete
input; it neither rejected nor repaired it.  The indexed-face provider now
rejects non-finite output, output outside the BREP surface envelope, and output
which omits the exact BREP boundary-vertex extent.  The shaded BREP cache key
is `brep-indexed-face-set-v3`, so the old partial payload cannot be replayed.
The exact tire diagnostic now returns an explicit provider failure in 0.55
seconds instead of publishing crossed-ring geometry.

The same boundary validation accepts all four NIST PMI7-10 solids.  A
fresh-v3-cache OSMesa matched run at
`/tmp/qged-lod-nist-brep-v3-20260901` completes the control in 18 seconds and
LoD view in 24 seconds.  It presents 178,561 faces with no boxes, proxies,
quality-floor violations, or terminal error; SSIM is 0.980645 and the
one-pixel-tolerant silhouette disagreement is zero.  This is a regression
against broad BREP suppression as well as a matched visual-quality result.

The `brepdraw` branch's bounded audit generates a finite shaded mesh for the
same source primitive (718 vertices and 1,336 triangles with one worker in
0.009 seconds), confirming that the source BREP is usable and that the missing
control is a lower-level shaded-display limitation in this branch.  Do not
expand this wheel result into a whole-train LoD claim.  First integrate or
supersede the applicable bounded BREP display work, then require the ordinary
and LoD paths to consume the same validated source tessellation contract before
running the full hierarchy comparison.  Post-hoc coverage validation cannot
bound a tessellator which fails to return or report which faces completed.  A
production provider therefore needs the bounded API's deadline, memory/result
limits, cancellation, per-face completion, and typed outcome in addition to
the output validation retained here.
