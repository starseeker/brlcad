# TLA+ verification and audit policy

## Pinned baseline

The initial reproducible toolchain is:

- TLA+ tools 1.7.4, TLC 2.19, revision `5a47802`;
- `tla2tools.jar` SHA-256
  `936a262061c914694dfd669a543be24573c45d5aa0ff20a8b96b23d01e050e88`;
- Java 11 or newer; and
- one TLC worker for baseline comparison and liveness checking.

The checked baseline is the sole source of recorded generated/distinct state
counts, search depth, and runtime observations.  Architecture documents may
summarize a successful gate but must not copy those numbers.

## Model review checklist

For every model and configuration, review and record:

1. Every variable is covered by `TypeOK`, and every action preserves it.
2. Success, rejection, supersession, constraint, and failure paths are
   reachable when they belong to the abstraction.
3. Required actions have nonzero TLC coverage, or their configuration-specific
   absence is documented.
4. `Init`, constants, constraints, and symmetry do not remove the interleaving
   the model claims to check.  Symmetry-sensitive models receive a smaller
   confirmation run with symmetry disabled.
5. `CHECK_DEADLOCK FALSE` has a documented terminal-stuttering reason; it is
   not a blanket waiver for an unintended deadlock.
6. Every temporal property names the input-closure condition under which it is
   promised and is checked for vacuity.
7. Every weak or strong fairness assumption maps to a concrete progress
   witness which production can schedule.
8. Every nonterminal state has a finite rank or an explicit external witness.
9. Ready, constrained, failed, and subordinate-background meanings agree with
   the canonical and parent composition models.
10. Every abstract production action has a typed reducer transition and an
    adversarial executable regression, or is recorded as active refinement
    debt.

## Composition review

Each composition model has a dependency table in `models.json`: its
`components` identify the focused contracts supplying premises, while each
component's `parents` identifies the flows which check its seams.  The
`compositionContracts` records reviewed, qualified assumption/guarantee
operators.  Lint verifies their existence, ownership, and component
membership.  Review must reject:

- a component with no parent flow unless it is explicitly an isolated
  data-plane boundary;
- a focused owner or terminal outcome absent from the canonical relation;
- a composition premise supplied only by prose;
- two components which can simultaneously own one semantic obligation; or
- a liveness cycle in which owners recreate one another's finite work without
  consuming a semantic member.

Required high-risk actions and the deadlock policy are likewise cataloged
under `semanticAudit`.  Those fields are gate inputs, not advisory prose: a
missing required coverage line, fair environment action, broad fairness on
`Next`, or unclassified deadlock configuration fails the run.

`conformance.json` is the implementation side of that contract.  Every
`semanticAudit.requiredActions` entry must map to an existing C++ regression;
stepwise mappings must execute the formal action name explicitly and check the
shared observation invariants after the action.  Run its Java-independent
drift check with:

```sh
cmake -P misc/CMake/ValidateBObolConformance.cmake
```

## Verification tiers

- **Lint:** catalog/schema/link checks and SANY parsing for all models.
- **Affected flow:** a changed model, all parent compositions named by the
  catalog, and the canonical/refinement pair.
- **Full:** all configurations with coverage and baseline comparison.
- **Tool upgrade:** complete old/new tool runs before changing the pin.
- **C++ conformance:** run `ctest -L bobol_conformance`; use the separately
  labeled sanitizer, graphics, and performance gates for behavior outside the
  bounded formal abstraction.

TLC metadata and logs belong under a caller-selected build or temporary
directory.  The runner also assigns a private Java temporary directory so
independent SANY/TLC invocations cannot race over extracted standard modules.
Updating the checked baseline is an explicit review action, never an automatic
consequence of an ordinary run.
