# Observational refinement checker for qged progressive-control reports.
# This filter must never be consulted by production scheduling.  It verifies
# that sampled concrete executions refine the finite revision/work contract;
# it deliberately does not mistake an unchanged, long-running worker for a
# cycle.

def required_control_fields:
  [
    "lod_control_obligation_mask",
    "lod_control_owner",
    "lod_control_violation_mask",
    "lod_inventory_revision",
    "lod_availability_revision",
    "lod_control_view_revision",
    "lod_control_policy_revision",
    "lod_capacity_revision",
    "lod_allocation_plan_serial",
    "lod_presentation_transaction_serial",
    "lod_presentation_required_render_serial",
    "lod_control_presented_frame_serial",
    "lod_view_revision",
    "lod_policy_revision",
    "lod_convergence_phase",
    "lod_convergence_outcome",
    "lod_convergence_terminal",
    "lod_convergence_view_ready"
  ];

def revision_tuple:
  [
    (.lod_inventory_revision // 0),
    (.lod_availability_revision // 0),
    (.lod_control_view_revision // 0),
    (.lod_control_policy_revision // 0),
    (.lod_capacity_revision // 0)
  ];

def control_state:
  [
    revision_tuple,
    (.lod_control_obligation_mask // 0),
    (.lod_control_owner // 0),
    (.lod_allocation_plan_serial // 0),
    (.lod_presentation_transaction_serial // 0),
    (.lod_presentation_required_render_serial // 0),
    (.lod_control_presented_frame_serial // 0),
    (.lod_convergence_phase // 0),
    (.lod_convergence_outcome // 0),
    (.lod_convergence_terminal // false),
    (.lod_convergence_view_ready // false)
  ];

def controller_samples:
  [.samples[] | select(.controller_available == true)];

# Commands, camera gestures, resize, and UI state changes are explicit input
# edges.  Waits and checkpoints only observe the controller.  A host may
# publish the input before the owner-thread synchronization which advances its
# semantic revision, so that one observation boundary is not a spontaneous
# terminal reopen.
def explicit_input:
  (.action // "") as $action |
  ($action != "") and
  ($action != "checkpoint") and
  (($action | startswith("wait")) | not);

def missing_field_violations($samples):
  [
    $samples | to_entries[] as $entry |
    required_control_fields[] as $field |
    select(($entry.value | has($field)) | not) |
    {
      kind: "missing-control-field",
      sample: $entry.key,
      field: $field
    }
  ];

def revision_regressions($samples):
  [
    range(1; $samples | length) as $index |
    ($samples[$index - 1] | revision_tuple) as $before |
    ($samples[$index] | revision_tuple) as $after |
    range(0; 5) as $domain |
    select($after[$domain] < $before[$domain]) |
    {
      kind: "revision-regression",
      sample: $index,
      domain: $domain,
      before: $before[$domain],
      after: $after[$domain]
    }
  ];

def serial_regressions($samples; $field; $kind):
  [
    range(1; $samples | length) as $index |
    ($samples[$index - 1][$field] // 0) as $before |
    ($samples[$index][$field] // 0) as $after |
    select($after < $before) |
    {
      kind: $kind,
      sample: $index,
      before: $before,
      after: $after
    }
  ];

# A semantic revision tuple certifies one retained allocation plan.  A zero
# serial means there is no committed plan and is intentionally ignored.
def duplicate_plan_violations($samples):
  [
    $samples |
    sort_by(revision_tuple) |
    group_by(revision_tuple)[] |
    . as $group |
    ([$group[].lod_allocation_plan_serial |
      select(. != null and . > 0)] | unique) as $plans |
    select($plans | length > 1) |
    {
      kind: "multiple-plans-for-revision",
      revision: ($group[0] | revision_tuple),
      plans: $plans
    }
  ];

# Once a sampled tuple is terminal it may remain terminal while background
# work drains.  It may not manufacture foreground work without a semantic
# revision edge.  Group entries retain original sample indices so sorting by
# the tuple cannot hide temporal order.
def terminal_reopen_violations($samples):
  [
    $samples |
    to_entries |
    sort_by(.value | revision_tuple) |
    group_by(.value | revision_tuple)[] |
    (sort_by(.key)) as $group |
    ([range(0; $group | length) |
      select($group[.].value.lod_convergence_terminal == true)] |
      first // null) as $terminal_offset |
    select($terminal_offset != null) |
    range($terminal_offset + 1; $group | length) as $offset |
    select($group[$offset].value.lod_convergence_terminal != true) |
    select((any($group[$terminal_offset + 1:$offset + 1][];
      (.value | explicit_input))) | not) |
    {
      kind: "terminal-reopened-without-revision",
      revision: ($group[0].value | revision_tuple),
      terminal_sample: $group[$terminal_offset].key,
      reopened_sample: $group[$offset].key
    }
  ];

# An A/B/A observation is a nonprogress cycle only when the complete semantic
# tuple, allocation plan, presentation transaction, and committed-frame
# witness are unchanged.  Identical repeated states are allowed: a bounded
# worker or cache write can legitimately outlive the diagnostic sampling
# interval.
def nonprogress_cycle_violations($samples):
  [
    range(2; $samples | length) as $index |
    $samples[$index - 2] as $first |
    $samples[$index - 1] as $middle |
    $samples[$index] as $last |
    select(($first | control_state) == ($last | control_state)) |
    select(($first | control_state) != ($middle | control_state)) |
    select(($first | revision_tuple) == ($middle | revision_tuple)) |
    select(($first.lod_allocation_plan_serial // 0) ==
      ($middle.lod_allocation_plan_serial // 0)) |
    select(($first.lod_presentation_transaction_serial // 0) ==
      ($middle.lod_presentation_transaction_serial // 0)) |
    select(($first.lod_control_presented_frame_serial // 0) ==
      ($middle.lod_control_presented_frame_serial // 0)) |
    select(all([$first, $middle, $last][];
      (.lod_interactive // false) == false)) |
    select(all([$first, $middle, $last][];
      (explicit_input | not))) |
    {
      kind: "nonprogress-control-cycle",
      samples: [$index - 2, $index - 1, $index],
      revision: ($first | revision_tuple),
      states: [($first | control_state), ($middle | control_state)]
    }
  ];

controller_samples as $samples |
(
  missing_field_violations($samples) +
  revision_regressions($samples) +
  serial_regressions($samples;
    "lod_presentation_transaction_serial";
    "presentation-transaction-regression") +
  serial_regressions($samples;
    "lod_control_presented_frame_serial";
    "presented-frame-regression") +
  duplicate_plan_violations($samples) +
  terminal_reopen_violations($samples) +
  nonprogress_cycle_violations($samples) +
  [
    $samples | to_entries[] as $entry |
    select(($entry.value.lod_control_violation_mask // 0) != 0) |
    {
      kind: "local-control-invariant",
      sample: $entry.key,
      mask: $entry.value.lod_control_violation_mask
    }
  ] +
  [
    $samples | to_entries[] as $entry |
    select(($entry.value.lod_control_view_revision // 0) !=
      ($entry.value.lod_view_revision // 0)) |
    {
      kind: "view-revision-observer-mismatch",
      sample: $entry.key
    }
  ] +
  [
    $samples | to_entries[] as $entry |
    select(($entry.value.lod_control_policy_revision // 0) !=
      ($entry.value.lod_policy_revision // 0)) |
    {
      kind: "policy-revision-observer-mismatch",
      sample: $entry.key
    }
  ]
) as $violations |
if ($violations | length) == 0 then
  true
else
  error($violations | tojson)
end
