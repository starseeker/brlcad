# Ayam NURBS editing survey and BRL-CAD implementation notes

This note records a source-level survey of the Ayam tree in `ayam/`, compares
its NURBS editing facilities with openNURBS/libbrep, and identifies suitable
paths for qged integration.  It also records the topology and licensing
boundaries that must govern any implementation.

## Executive conclusions

- Ayam's useful NURBS algorithms do not fundamentally depend on the Ayam
  `NPatch`, `ConcatNP`, or other scene-object structures.  Those structures
  supply storage, construction history, notification, property UI, and
  conversion callbacks.  Curve and surface algorithms can be re-expressed
  against `ON_NurbsCurve`, `ON_NurbsSurface`, and `ON_Brep`.
- Reproducing Ayam's construction history is a separate product decision.
  One-shot operations in libbrep/libged do not need an Ayam-style scene DAG.
  A persistent parametric `Skin`, `Sweep`, or `ConcatNP` equivalent would need
  a BRL-CAD dependency/history representation in addition to the geometry
  algorithm.
- Ayam is a surface modeler, not a manifold B-rep editor.  Its patch
  notification code repairs patch-local properties and regenerates dependent
  construction objects, but it does not reconcile a shared `ON_BrepEdge`,
  both faces, their 2-D trims, and their 3-D edge curve after a boundary edit.
- Direct CV editing must consequently be conservative in BRL-CAD.  An edit is
  safe without topology repair when the edited basis function has zero
  support on every trim of the face.  A new exact backend also admits a
  strict class of two-face isoparametric C0 edits; other trim-influencing CVs
  remain visible and selectable but locked.
- Eigen is sufficient for the linear systems in interpolation,
  approximation, compatibility, and first-generation fairness tools.  No
  additional solver dependency is currently justified.  Robust B-rep
  boundary work is primarily a geometry, constraint, and validation problem,
  not a missing linear solver problem.

## How Ayam represents and updates NURBS

`NPatch` is Ayam's directly editable tensor-product NURBS patch object.  It
stores dimensions, orders, closure flags, knot types/vectors, rational
Euclidean CVs, display settings, and child trim objects.  Its notification
callback in `ayam/src/objects/npatch.c` checks or regenerates knots and closure
state, rebuilds multiple-point data, caps and bevels, and invalidates display
and tessellation caches.

Construction objects such as `ConcatNP`, `Skin`, `Sweep`, `Swing`, `Birail1`,
`Birail2`, `Gordon`, `Revolve`, `Extrude`, `Cap`, and `Bevel` own parameters
and children.  Ayam's notify/provide mechanism in
`ayam/src/aycore/notify.c` regenerates their resulting `NPatch` objects after
a child or parameter change.  The resulting patch is therefore an output
cache of a scene-graph operation.

This architecture divides cleanly into three layers:

| Ayam responsibility | BRL-CAD expression |
|---|---|
| NURBS curve/surface math | libbrep, implemented on openNURBS geometry |
| Valid B-rep mutation and validation | libbrep transaction operations |
| User-visible command and database transaction | librt edit descriptors and libged |
| Picking, handles, control net, preview, property UI | reusable libqtcad filters/widgets |
| Optional persistent construction history | future database/graph design, not required for one-shot tools |

Ayam's object types therefore do not preclude implementation in BRL-CAD.
They should not be mechanically transplanted: the geometry algorithms and
their behavior are the reusable concepts, while ownership and notification
must follow BRL-CAD's database and edit-transaction model.

## Ayam capability inventory

### Direct point editing and inspection

Ayam displays control polygons/nets and control points, supports single and
multi-point selection, screen/region-based selection, point transforms,
rational weight editing and weight visualization, and coincident
"multiple-point" handling.  It also supplies curvature plots, boundary
display, U-direction display, and surface/curve property panels.

These are good libqtcad candidates.  Selection and drag behavior belongs in a
reusable view filter; indices, coordinates, weights, lock state, and numeric
entry belong in a small editor panel.  All mutation should be routed through
the same librt/libbrep operation used by command-line editing.

### Shape-preserving representation edits

Ayam supports the following curve operations:

- reverse, open/close, and shift the seam of closed curves;
- insert/refine knots, refine using another curve, and remove knots;
- remove superfluous knots, clamp/unclamp, and rescale knot domains;
- elevate and reduce degree/order;
- make curves compatible;
- split, trim, concatenate, extend, and extract;
- reparameterize, interpolate, approximate, fair, collapse, and explode
  points.

Its corresponding surface tools include:

- reverse U/V and swap U/V;
- open/close U/V;
- insert/refine/remove knots and remove superfluous knots;
- clamp/unclamp, elevate/reduce degree, and rescale knot domains;
- make surfaces compatible;
- split and extract subpatches;
- interpolate, approximate, fair, cap, and build/break surfaces from/to
  curves.

Knot insertion/refinement, clamping, degree elevation, parameter reversal,
and transposition are the strongest initial candidates.  openNURBS already
implements their geometry kernels, and they preserve the represented locus.
For a face, parameter reversal or transposition must also transform its 2-D
trims; `ON_BrepFace` supplies that face-aware operation.

Knot removal, degree reduction, approximation, and fairness are good later
candidates but require an explicit error tolerance, rank/conditioning policy,
and regression tests.  A surface-level operation is not automatically a
valid face-level operation: any operation that changes the surface locus at a
trim requires coupled topology handling.

### Curve and surface construction

Ayam provides:

- circles, circular B-splines, rectangles, trim rectangles, and tween curves;
- ruled/build-from-curves, revolve, extrude, swing, sweep, birail, skin,
  Gordon, tween-surface, cap, bevel, concatenate, extract, and offset
  operations;
- projection of curves to surfaces and creation/use of trim curves;
- conversion among direct NURBS objects and construction objects.

Most are reasonable one-shot libbrep operations.  Revolve, ruled/tensor
product, extraction, and basic interpolation already have partial libbrep
support.  Skin, sweep, birail, Gordon, compatibility, and concatenate are
medium-sized additions: the core math is conventional, but tolerances,
orientation, degenerate inputs, rational coordinates, knot compatibility,
and the topology of the result need BRL-CAD-specific treatment.

Offset surfaces, automatic bevel/cap construction, projection-to-trim, and
arbitrary trimmed-surface concatenation are poor early candidates.  A useful
industrial implementation must handle singularities, self-intersections,
multiple loops, curve pullbacks, edge tolerances, and failure without
partially corrupting the input B-rep.

## Topology preservation

Ayam's closure routines in `ayam/src/nurbs/npt.c` copy or average seam rows
and depend on a compatible knot vector.  This maintains the local tensor
product patch convention.  It is not a global manifold edit.  Similarly,
the `NPatch` notification callback rebuilds patch-local cached and generated
objects; it does not discover neighboring B-rep faces and solve a shared
boundary constraint.

For an openNURBS B-rep, an arbitrary trimmed boundary has at least four
coupled representations:

1. the first face's surface and its 2-D trim;
2. the shared 3-D edge curve and vertices;
3. the second face's 2-D trim;
4. the second face's surface.

Moving a face CV can change the surface locus under its trim.  The new locus
need not lie on the adjacent surface, so merely changing the edge curve or
pulling back a new trim is insufficient.  The adjacent surface may need a
constrained deformation or the faces may need re-intersection.  Continuity
beyond C0 adds derivative constraints.

The implemented policy classifies both edges and selected CVs.  Interior CVs
use the direct face-aware edit path.  When a CV influences exactly one mated
edge, the constrained path may:

1. verify that both trims are monotone isoparametric curves on distinct,
   unshared, non-transposed NURBS surfaces;
2. extract each surface's exact NURBS isocurve over the trim interval and
   orient it to the shared edge;
3. require compatible trace degree, knots, CV count, rational state, and
   weights;
4. perturb candidate surface CVs to construct their exact linear influence
   on the trace CVs;
5. use Eigen `CompleteOrthogonalDecomposition` to find adjacent-face CV
   translations reproducing the requested trace displacement;
6. reject endpoint movement, secondary-trim influence, excessive
   coefficients, or an unacceptable residual;
7. update both surfaces and the shared 3-D edge in a trial B-rep; and
8. commit only if the traces agree and the complete B-rep remains valid.

This does not assume corresponding surface rows.  It supports compatible
natural clamped boundaries and compatible interior isoparametric trims
through one backend.  It currently keeps pcurves and rational weights fixed,
and rejects periodic seams, vertices/junctions, multiple affected trims,
one-sided isoparametric edges, arbitrary trims, shared underlying surfaces,
and incompatible trace spaces.

A local census of shared edges in `db/nist/NIST_MBE_PMI_{1,2,4}.stp`
illustrates why this exact case should be a backend rather than the whole
design:

| Model | Shared edges | Natural on both | Isoparametric on both | Isoparametric on at least one |
|---|---:|---:|---:|---:|
| NIST 1 | 310 | 24 (7.7%) | 188 (60.6%) | 296 (95.5%) |
| NIST 2 | 1,475 | 292 (19.8%) | 672 (45.6%) | 1,295 (87.8%) |
| NIST 4 | 1,156 | 241 (20.8%) | 767 (66.3%) | 1,140 (98.6%) |
| Combined | 2,941 | 557 (18.9%) | 1,627 (55.3%) | 2,731 (92.9%) |

These figures classify trim geometry only; they do not say that every edge
has compatible trace spaces or satisfies the current edit acceptance rules.
They do show that a natural-only implementation would be too narrow, while
exact isoparametric handling has substantial value inside a general
affected-topology constraint graph.

Recommended capability stages are:

1. Permit a CV edit only when its tensor-product basis support does not
   intersect any face trim.  Implemented.
2. Add exact C0 coupling for compatible isoparametric traces, including
   natural boundaries.  A strict transactional subset is implemented.
3. Generalize C0 propagation with an affected-topology graph, exact
   constraints where available, adaptive geometric constraints for arbitrary
   trims, and face re-intersection/trim reconstruction when deformation is
   insufficient.
4. Add optional G1/G2 derivative constraints and coupled adjacent rows.

This policy is intentionally stricter than Ayam.  It supports ordinary
interior and eligible exact-seam editing without implying that arbitrary
manifold boundary editing has been solved.

## Solver assessment

Ayam's approximation code forms normal equations in places and calls its own
solver.  Its `ayam/src/nurbs/act.c` also contains an SVD described as an ANSI C
rendering of Golub-Reinsch ALGOL code "believed to be in the public domain."
That is the ambiguous academic-code reference and should not be copied.

Eigen covers the anticipated BRL-CAD needs:

- `HouseholderQR` or `ColPivHouseholderQR` for ordinary dense least squares;
- `CompleteOrthogonalDecomposition` for rank-deficient systems;
- `JacobiSVD` for small ill-conditioned fitting systems and diagnostics;
- `SparseQR`, `SimplicialLDLT`, or an iterative sparse solver for larger
  fairness/constraint grids.

New code should solve the least-squares system directly rather than form
`N^T N`, which squares the condition number.  Exact interpolation and
tridiagonal systems do not require an external package.  Boundary
deformation can initially use linearized equality-constrained least squares;
the hard part remains defining geometrically correct constraints and
validating the result.

No third-party solver should be added now.  A nonlinear least-squares or
large sparse direct dependency should be considered only after a concrete
operation and benchmark demonstrate that Eigen is inadequate, followed by a
license review of the exact package and all optional backends.

The implemented exact C0 trace coupling uses
`CompleteOrthogonalDecomposition`, including residual and coefficient
checks.  Its successful natural and interior-isoparametric fixtures confirm
that Eigen is adequate for this stage; geometry coverage, not solver
capability, is the present limitation.

## Licensing cautions

`ayam/License.txt` is a permissive notice allowing use, modification,
distribution, and relicensing, provided existing copyright notices are
retained and the notice is included verbatim in distributions.  Ayam-authored
code is broadly compatible with BRL-CAD's use, but copied or closely derived
code must carry the required attribution and notice.

The top-level notice does not eliminate provenance concerns in individual
files:

- `ayam/src/nurbs/act.c` says its Golub-Reinsch SVD is only *believed* to be
  public domain.
- `ayam/src/nurbs/nb.c` identifies code adapted from *The NURBS Book*, code
  derived from NURBS++, and LU code originating in Skip Carter's Kalman
  library.  Those chains need independent license verification before any
  copying.
- `ayam/src/nurbs/rtess.c` cites example code associated with an academic
  tessellation paper/project, and `stess.c` says an intersection routine was
  taken from the `c.g.algorithms` FAQ.
- Other parts of the Ayam distribution include their own notices, including
  Affine/RenderMan and platform/library code.  They are not made safe for
  BRL-CAD merely by residing in the Ayam tree.

The preferred policy is:

1. use openNURBS or existing BRL-CAD implementations where available;
2. implement published mathematical algorithms independently from their
   descriptions, using BRL-CAD naming, data structures, tests, and error
   policy;
3. use Eigen for linear algebra instead of transplanting Ayam solvers;
4. copy Ayam implementation code only after file-level provenance review and
   with the complete required attribution.

## Current implementation tranches

The initial implementation establishes:

- reference-safe curve/surface removal and shape-preserving replacement in
  libbrep;
- curve/surface rational weight, clamp, degree-elevation, reverse,
  transpose, and knot-insertion operations;
- face-aware CV position/translation/weight edits that isolate shared
  surfaces and reject unsupported trim constraints;
- edge/CV constraint classification and command-line summaries;
- exact transactional C0 coupling for eligible natural and interior
  isoparametric two-face seams, including adjacent-face solve and shared-edge
  regeneration;
- face-aware reverse/transpose operations;
- libged `brep ... geo` commands plus whole-B-rep validation;
- librt edit descriptors for CV selection, move, absolute position, and
  rational weight;
- reusable libqtcad CV pick/drag filters with interior, exact-C0-coupled, and
  locked status;
- a qged edit panel with an in-memory control-net/wireframe preview and
  explicit Apply/Reset transaction.

Regression tests exercise reference remapping, shape invariance after knot
insertion/degree elevation, rational Euclidean CV behavior, shared-surface
isolation, rejection and rollback of unsupported trim-affecting edits, exact
natural and interior-isoparametric seam coupling, preserved endpoints and
topology-array counts, complete B-rep validity, and the coupled librt path.

The tranche deliberately does not claim arbitrary boundary propagation,
trace compatibility refinement, rational boundary-weight or junction edits,
G1/G2 coupled editing, construction history, knot removal, degree reduction,
fairing, fitting, or robust offset/bevel/cap generation.

## Ayam source map

The primary files used for this survey are:

- `ayam/License.txt`
- `ayam/doc/html/ayam-4.html` (object types and properties)
- `ayam/doc/html/ayam-5.html` (NURBS modeling tools)
- `ayam/src/objects/npatch.c`
- `ayam/src/objects/concatnp.c`
- `ayam/src/objects/{skin,sweep,swing,birail1,birail2,gordon}.c`
- `ayam/src/aycore/notify.c`
- `ayam/src/nurbs/{nct,npt,knots,act,apt,nb}.c`
