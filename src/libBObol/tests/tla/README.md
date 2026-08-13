# BObol LoD TLA+ models

These modules model the control-plane contracts of BObol's progressive LoD
implementation. They intentionally do not model geometry generation,
floating-point projection, perceptual quality, or real renderer timing.

`LodCoordinator.tla` is a conformance model for
`BObolLodStateMachine`. It exhausts the ten Boolean observations, every public
event, and the same three counter boundary cases used by
`test_lod_coordinator.cpp`. Each arbitrary observation is followed by a valid
quiet observation to check bounded discharge from contradictory input.

`LodScene.tla` composes three occurrences, three cuts, asynchronous requests,
epoch admission, coverage, allocation, publication, frame completion,
calibration, memory pressure, and compaction. The deliberately tiny domains
make all orderings checkable by TLC. All three occurrences have the same cost
curve; their differing importance makes `ImportanceOrdering` a direct check
that an equivalent lower-importance candidate cannot win because of scan
order. Equal-score allocations may differ, but their checked utility and
protected-floor behavior must not.

## Running TLC

A Java runtime and an external `tla2tools.jar` are required. The jar is not
committed to BRL-CAD. A direct coordinator run is:

```sh
java -Xmx2g -cp /path/to/tla2tools.jar tlc2.TLC \
  -cleanup -workers 1 -config LodCoordinator.cfg LodCoordinator
```

Run that command from this directory. Scene configurations use the same form,
for example `-config LodScene_Safety.cfg LodScene`.

For CTest integration, configure with either the cache value or environment
variable `TLA2TOOLS_JAR`. An optional `TLA2TOOLS_SHA256` cache value or
environment variable pins the exact jar:

```sh
cmake -S . -B build \
  -DTLA2TOOLS_JAR=/path/to/tla2tools.jar \
  -DTLA2TOOLS_SHA256=<sha256>
ctest --test-dir build -L bobol_tla --output-on-failure
```

When Java or the jar is absent, the TLA+ tests are not registered and ordinary
builds are unaffected. TLC prints its release in each test log, while the test
runner always prints the jar SHA-256. Record both when publishing a
counterexample or updating the expected results.

## Configurations

- `LodCoordinator.cfg`: exhaustive current coordinator contract.
- `LodScene_Safety.cfg`: disruptive input, stale completion, allocation,
  publication, budget, and compaction safety.
- `LodScene_Liveness.cfg`: weakly fair internal work after a fair one-way
  environment quiescence action.
- `LodScene_Pressure.cfg`: pressure is present before minimum coverage.
- `LodScene_ResidentFloor.cfg`: the three minimum payloads require six memory
  units while the soft resident limit is five.
- `LodScene_Cancelled.cfg`: a cancelled source begins with a late result from
  the retired generation.
- `LodCoordinator_PrematureStable.cfg`, `LodScene_StaleMutation.cfg`, and
  `LodScene_PrematureStable.cfg`: expected-failure mutation checks proving
  that the principal invariants are capable of rejecting the named bugs.

## Code correspondence

| Model concept | Current implementation contract |
| --- | --- |
| Coordinator phases/events and canonical input | `BObolLodStateMachine` in `lod_coordinator_private.h` |
| Coverage cursor and completed proof | `BObolLodCoveragePolicy` |
| Quiet readiness and terminal limitation | `BObolLodConvergencePolicy` |
| Scene allowance and calibration | `BObolLodBudgetPolicy` |
| Pressure edge and one-shot recovery | `BObolLodResourcePolicy` |
| Resident maintenance admission | `BObolLodCompactionPolicy` |
| Request/result identity | `BLodIdentifiers.h`, `lod_service.cpp`, and result admission in `view_lod.cpp` |
| Importance and protected floor | retained allocator in `view_controller.cpp` |
| Publication/frame barrier | `BObolLodPublicationPolicy` and controller frame completion |
| Executable comparison tests | `test_lod_coordinator.cpp` |

The repository-level `tla_overview.txt` records the rationale, abstraction
boundary, known omissions, and the procedure for reconciling these initial
models with ongoing LoD development.
