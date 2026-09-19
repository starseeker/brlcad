# BObol formal model suite

This directory is the authoritative entry point for BObol's TLA+ models.  TLC
exhaustively checks the configured finite state spaces.  Those checks are not
proofs of the unbounded C++ program: executable refinement tests, renderer
tests, image comparisons, resource measurements, and production traces remain
required evidence.

## Proof hierarchy

The suite has three tiers:

1. `canonical/` defines the finite progressive-control contract and the
   concrete-fact refinement projection.
2. `compositions/` checks that adjacent component contracts have compatible
   ownership, progress, and terminal assumptions.
3. `components/` isolates one focused protocol boundary.

Physical placement records a model's tier and primary domain.  It cannot
express the complete many-to-many proof relationship.  `models.json` is the
sole machine-readable catalog of parent flows, component membership,
production owners, executable tests, and verification cost.

The component domains are:

- `lifecycle/`: host scheduling, user interaction, policy lifecycle, and
  bounded submission ownership;
- `planning-quality/`: admission, allocation, capacity, convergence, and
  terminal-quality policy;
- `retained-publication/`: retained CAD mutation and exact presentation; and
- `asset-production-residency/`: immutable producer, cache, live-page, and
  resident-memory protocols.

`AUDIT.md` records the semantic audit and its deliberate proof boundary.
`FAIRNESS.md`, `TERMINOLOGY.md`, and `RISK_COVERAGE.md` define the progress
contract, local vocabulary, and formal-to-implementation evidence roadmap.
They are intentionally separate from historical engineering lessons and from
generated TLC results.

## Required model header

Every model must state, in its leading comment or catalog record:

- the behavior it models and deliberately excludes;
- its tier, parent composition flows, and production owner;
- the concrete state or C++ transition which refines each abstract fact;
- its finite bounds and what they abstract;
- the input-closure and fairness assumptions used by liveness properties;
- its finite progress measure and terminal outcomes; and
- the executable regression tests which protect the refinement boundary.

Focused models may add detail to a canonical fact.  They may not create an
additional production owner, revision meaning, control mode, or terminal
outcome without first changing the canonical and applicable composition
models.

## Running checks

The suite is pinned initially to TLA+ tools 1.7.4 / TLC 2.19.  Supply the JAR
explicitly; state files and logs are always written outside the source tree.

```sh
cmake \
  -DTLA2TOOLS_JAR=/path/to/tla2tools.jar \
  -DTLA_MODE=lint \
  -P misc/CMake/RunTLA.cmake

cmake \
  -DTLA2TOOLS_JAR=/path/to/tla2tools.jar \
  -DTLA_MODE=check \
  -DTLA_FLOW=cad-frame \
  -DTLA_WORK_DIR=/tmp/brlcad-tla \
  -P misc/CMake/RunTLA.cmake
```

Selectors are `TLA_MODEL`, `TLA_FLOW`, and `TLA_TIER`.  With no selector,
`check` runs the complete catalog.  `lint` always validates and parses the
complete suite.  A flow check includes the canonical/refinement pair unless a
tier or model selector explicitly excludes them.  `TLA_COMPARE_BASELINE=ON`
rejects state-count or search-depth drift.  After reviewing a successful,
unfiltered tool-upgrade or semantic change run, `TLA_ACCEPT_BASELINE=ON`
explicitly replaces the checked baseline.  See `verification.md` for the
review checklist and baseline policy.

## Change discipline

- Move or rename models without changing their transition relation; verify
  identical state counts before semantic work.
- Change one ownership/refinement boundary at a time and update its parent
  compositions and C++ regression in the same series.
- Turn every meaningful TLC counterexample into an executable regression.
- Keep numeric quality, geometry, elapsed time, and memory-size claims out of
  control models; test those properties in the data plane.
- Never copy TLC state directories or generated reports into this tree.
