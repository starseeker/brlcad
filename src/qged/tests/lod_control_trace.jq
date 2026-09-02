# Observational refinement checker for qged progressive-control reports.
# This filter must never be consulted by production scheduling.  It verifies
# that sampled concrete executions refine the finite revision/work contract;
# it deliberately does not mistake an unchanged, long-running worker for a
# cycle.

def required_control_fields:
  [
    "lod_control_fact_mask",
    "lod_control_obligation_mask",
    "lod_control_owner",
    "lod_control_violation_mask",
    "lod_control_presentation_witness_mask",
    "host_work_flags",
    "lod_inventory_revision",
    "lod_availability_revision",
    "lod_visibility_revision",
    "lod_control_view_revision",
    "lod_control_policy_revision",
    "lod_capacity_revision",
    "lod_cad_revision",
    "lod_resident_demand_revision",
    "lod_allocation_plan_serial",
    "lod_presentation_transaction_serial",
    "lod_presentation_required_render_serial",
    "lod_control_presented_frame_serial",
    "lod_capacity_search_phase",
    "lod_capacity_search_goal",
    "lod_capacity_search_samples_remaining",
    "lod_capacity_search_measured_candidates",
    "lod_capacity_search_total_measured_candidates",
    "lod_capacity_search_candidate_limit",
    "lod_capacity_search_maximum_candidates",
    "lod_capacity_search_sample_limit",
    "lod_capacity_search_completed_units",
    "lod_capacity_search_total_units",
    "lod_renderer_preparation_target_signature",
    "lod_renderer_preparation_total_units",
    "lod_renderer_preparation_completed_units",
    "lod_renderer_preparation_remaining_units",
    "lod_renderer_preparation_reserved_bytes",
    "lod_renderer_preparation_target_count",
    "lod_renderer_preparation_preparing_targets",
    "lod_renderer_preparation_constrained_targets",
    "lod_renderer_preparation_failed_targets",
    "lod_renderer_preparation_invalid_targets",
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
    (.lod_visibility_revision // 0),
    (.lod_control_view_revision // 0),
    (.lod_control_policy_revision // 0),
    (.lod_capacity_revision // 0),
    (.lod_cad_revision // 0),
    (.lod_resident_demand_revision // 0)
  ];

# Allocation authority is owned by the six admission domains.  CAD binding and
# resident-demand revisions validate an in-flight transaction, but applying a
# plan advances them and must not turn that plan's output into a new allocation
# problem.
def allocation_revision_tuple:
  [
    (.lod_inventory_revision // 0),
    (.lod_availability_revision // 0),
    (.lod_visibility_revision // 0),
    (.lod_control_view_revision // 0),
    (.lod_control_policy_revision // 0),
    (.lod_capacity_revision // 0)
  ];

def control_state:
  [
    revision_tuple,
    (.lod_control_fact_mask // 0),
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

# Transaction and frame serials are clocks, not finite convergence ranks.  A
# defective controller can alternate equivalent semantic states while issuing
# a fresh transaction and frame each time.  Cycle identity therefore excludes
# those serials; strict_progress_between names the bounded work which may
# legitimately justify returning to an earlier control state.
def cycle_state:
  [
    revision_tuple,
    (.lod_control_fact_mask // 0),
    (.lod_control_obligation_mask // 0),
    (.lod_control_owner // 0),
    (.lod_allocation_plan_serial // 0),
    (.lod_renderer_preparation_target_signature // 0),
    (.lod_convergence_phase // 0),
    (.lod_convergence_outcome // 0),
    (.lod_convergence_terminal // false),
    (.lod_convergence_view_ready // false)
  ];

def unfinished_worker_units:
  (.lod_convergence_pending_tasks // 0) +
  (.lod_convergence_in_flight // 0) +
  (.lod_convergence_queued_results // 0);

def submission_cursor:
  [
    (.lod_submission_source_index // 0),
    (.lod_submission_entry_offset // 0)
  ];

def strict_progress_between($before; $after):
  (($after.lod_capacity_search_completed_units // 0) >
    ($before.lod_capacity_search_completed_units // 0)) or
  (($after | submission_cursor) > ($before | submission_cursor)) or
  (($after.lod_convergence_available_leaves // 0) >
    ($before.lod_convergence_available_leaves // 0)) or
  (($after.lod_convergence_satisfied_payloads // 0) >
    ($before.lod_convergence_satisfied_payloads // 0)) or
  (($after.lod_convergence_source_preparation_providers // 0) <
    ($before.lod_convergence_source_preparation_providers // 0)) or
  ((($before.lod_renderer_preparation_target_signature // 0) != 0) and
   (($after.lod_renderer_preparation_target_signature // 0) ==
    ($before.lod_renderer_preparation_target_signature // 0)) and
   (($after.lod_renderer_preparation_remaining_units // 0) <
    ($before.lod_renderer_preparation_remaining_units // 0))) or
  (($after.lod_renderer_preparation_target_count // 0) <
    ($before.lod_renderer_preparation_target_count // 0)) or
  (($after | unfinished_worker_units) <
    ($before | unfinished_worker_units)) or
  (($after.lod_convergence_compactions // 0) >
    ($before.lod_convergence_compactions // 0)) or
  (($after.lod_convergence_compaction_candidates // 0) <
    ($before.lod_convergence_compaction_candidates // 0)) or
  (($after.lod_convergence_fraction // 0) >
    (($before.lod_convergence_fraction // 0) + 0.000001));

def fact_active($mask; $bit):
  ((($mask / $bit) | floor) % 2) == 1;

# Keep this projection numerically identical to
# BObolLodControlRefinement::evaluate.  jq has no portable bitwise operators,
# so masks are assembled from exact powers of two.  The maximum fact bit is
# 2^28, which remains exactly representable by jq's numeric type.
def projected_obligation_mask($facts):
  (if fact_active($facts; 1) then 1 else 0 end) +
  (if fact_active($facts; 2) then 2 else 0 end) +
  (if fact_active($facts; 4) then 4 else 0 end) +
  (if fact_active($facts; 8) or fact_active($facts; 16)
   then 8 else 0 end) +
  (if any([32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
           16384, 32768, 65536, 134217728][];
          fact_active($facts; .))
   then 16 else 0 end) +
  (if any([131072, 262144, 524288, 1048576, 2097152, 4194304,
           8388608, 268435456][];
          fact_active($facts; .))
   then 32 else 0 end) +
  (if fact_active($facts; 16777216) then 64 else 0 end) +
  (if fact_active($facts; 33554432) then 128 else 0 end) +
  (if fact_active($facts; 67108864) then 256 else 0 end);

def projected_owner($facts; $obligations):
  if fact_active($facts; 131072) then 6
  elif fact_active($facts; 1) then 1
  elif fact_active($facts; 8) or fact_active($facts; 16) then 4
  elif fact_active($facts; 2) then 2
  elif fact_active($facts; 4) then 3
  elif fact_active($facts; 8192) then 5
  elif fact_active($obligations; 32) then 6
  elif fact_active($obligations; 16) then 5
  elif fact_active($facts; 16777216) then 7
  elif fact_active($facts; 33554432) then 8
  elif fact_active($facts; 67108864) then 9
  else 0
  end;

def concrete_projection_violations($samples):
  [
    $samples | to_entries[] as $entry |
    ($entry.value.lod_control_fact_mask // 0) as $facts |
    (projected_obligation_mask($facts)) as $projected_obligations |
    ($entry.value.lod_control_obligation_mask // 0) as $actual_obligations |
    select($actual_obligations != $projected_obligations) |
    {
      kind: "concrete-obligation-mismatch",
      sample: $entry.key,
      facts: $facts,
      expected: $projected_obligations,
      actual: $actual_obligations
    }
  ] +
  [
    $samples | to_entries[] as $entry |
    ($entry.value.lod_control_fact_mask // 0) as $facts |
    (projected_obligation_mask($facts)) as $projected_obligations |
    (projected_owner($facts; $projected_obligations)) as $projected_owner |
    ($entry.value.lod_control_owner // 0) as $actual_owner |
    select($actual_owner != $projected_owner) |
    {
      kind: "concrete-owner-mismatch",
      sample: $entry.key,
      facts: $facts,
      obligations: $projected_obligations,
      expected: $projected_owner,
      actual: $actual_owner
    }
  ];

# A selective source-delta scope is a mode of the bounded submission cursor.
# It cannot survive cursor completion as independent planning debt: doing so
# blocks the scene-wide demand pass which must subsume the completed delta.
def planning_producer_violations($samples):
  [
    $samples | to_entries[] as $entry |
    ($entry.value.lod_control_fact_mask // 0) as $facts |
    select(fact_active($facts; 128)) |
    select(fact_active($facts; 32) | not) |
    {
      kind: "selective-scope-without-submission-cursor",
      sample: $entry.key,
      facts: $facts
    }
  ];

def controller_samples:
  (if ((.control_transitions // []) | length) > 0
   then .control_transitions else (.samples // []) end) |
  map(select(.controller_available == true));

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
    range(0; 8) as $domain |
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

# Each retained target has an immutable denominator and monotone completion
# rank.  Check that concrete aggregate independently of the controller's
# categorical STARTED/ADVANCED result so a malformed producer cannot validate
# its own retry.
def preparation_certificate_violations($samples):
  [
    $samples | to_entries[] as $entry |
    ($entry.value.lod_renderer_preparation_target_signature // 0) as $signature |
    ($entry.value.lod_renderer_preparation_target_count // 0) as $targets |
    ($entry.value.lod_renderer_preparation_total_units // 0) as $total |
    ($entry.value.lod_renderer_preparation_completed_units // 0) as $completed |
    ($entry.value.lod_renderer_preparation_remaining_units // 0) as $remaining |
    ($entry.value.lod_renderer_preparation_preparing_targets // 0) as $preparing |
    ($entry.value.lod_renderer_preparation_constrained_targets // 0) as $constrained |
    ($entry.value.lod_renderer_preparation_failed_targets // 0) as $failed |
    ($entry.value.lod_renderer_preparation_invalid_targets // 0) as $invalid |
    select(
      (($signature == 0) != ($targets == 0)) or
      ($completed > $total) or
      ($remaining != ($total - $completed)) or
      (($preparing + $constrained + $failed) > $targets) or
      ($invalid != 0)) |
    {
      kind: "invalid-renderer-preparation-certificate",
      sample: $entry.key,
      signature: $signature,
      targets: $targets,
      total: $total,
      completed: $completed,
      remaining: $remaining,
      preparing: $preparing,
      constrained: $constrained,
      failed: $failed,
      invalid: $invalid
    }
  ] +
  [
    range(1; $samples | length) as $index |
    $samples[$index - 1] as $before |
    $samples[$index] as $after |
    ($before.lod_renderer_preparation_target_signature // 0) as $signature |
    select($signature != 0) |
    select(($after.lod_renderer_preparation_target_signature // 0) ==
      $signature) |
    select(
      (($after.lod_renderer_preparation_total_units // 0) !=
       ($before.lod_renderer_preparation_total_units // 0)) or
      (($after.lod_renderer_preparation_completed_units // 0) <
       ($before.lod_renderer_preparation_completed_units // 0)) or
      (($after.lod_renderer_preparation_remaining_units // 0) >
       ($before.lod_renderer_preparation_remaining_units // 0))) |
    {
      kind: "renderer-preparation-rank-regression",
      sample: $index,
      signature: $signature,
      before: [
        ($before.lod_renderer_preparation_total_units // 0),
        ($before.lod_renderer_preparation_completed_units // 0),
        ($before.lod_renderer_preparation_remaining_units // 0)
      ],
      after: [
        ($after.lod_renderer_preparation_total_units // 0),
        ($after.lod_renderer_preparation_completed_units // 0),
        ($after.lod_renderer_preparation_remaining_units // 0)
      ]
    }
  ];

# The capacity HUD rank is derived from the same finite candidate protocol as
# the production search.  Validate that projection independently so a broken
# progress display cannot conceal an unbounded or regressing implementation.
def capacity_search_certificate_violations($samples):
  [
    $samples | to_entries[] as $entry |
    ($entry.value.lod_capacity_search_phase // -1) as $phase |
    ($entry.value.lod_capacity_search_goal // -1) as $goal |
    ($entry.value.lod_capacity_search_samples_remaining // -1) as $remaining |
    ($entry.value.lod_capacity_search_measured_candidates // -1) as $measured |
    ($entry.value.lod_capacity_search_total_measured_candidates // -1) as $total_measured |
    ($entry.value.lod_capacity_search_candidate_limit // -1) as $candidate_limit |
    ($entry.value.lod_capacity_search_maximum_candidates // -1) as $maximum_candidates |
    ($entry.value.lod_capacity_search_sample_limit // -1) as $sample_limit |
    ($entry.value.lod_capacity_search_completed_units // -1) as $completed |
    ($entry.value.lod_capacity_search_total_units // -1) as $total |
    ($sample_limit + 2) as $candidate_units |
    (if $phase == 0 then 0
     elif $phase == 1 then $total_measured * $candidate_units
     elif $phase == 2 then $total_measured * $candidate_units + 1
     elif $phase == 3 then
       $total_measured * $candidate_units + 2 +
         ($sample_limit - $remaining)
     elif $phase == 4 then $total
     else -1
     end) as $expected_completed |
    select(
      ($phase < 0 or $phase > 4) or
      ($goal < 0 or $goal > 1) or
      ($candidate_limit <= 0) or
      ($sample_limit <= 0) or
      ($measured < 0 or $measured > $candidate_limit) or
      ($total_measured < 0 or $total_measured > $maximum_candidates) or
      ($remaining < 0 or $remaining > $sample_limit) or
      ($completed < 0 or $completed > $total) or
      (if $phase == 0 then
         ($maximum_candidates != 0 or $completed != 0 or $total != 0 or
          $remaining != 0 or $measured != 0 or $total_measured != 0)
       else
         ($maximum_candidates <= 0 or
          ($maximum_candidates != $candidate_limit and
           $maximum_candidates != (2 * $candidate_limit)) or
          $total != ($maximum_candidates * $candidate_units) or
          $completed != $expected_completed or
          (($phase != 3) and $remaining != 0) or
          (($phase >= 1 and $phase <= 3) and $completed >= $total))
       end)) |
    {
      kind: "invalid-capacity-search-certificate",
      sample: $entry.key,
      phase: $phase,
      goal: $goal,
      completed: $completed,
      expected_completed: $expected_completed,
      total: $total,
      measured: $measured,
      total_measured: $total_measured,
      samples_remaining: $remaining
    }
  ];

# A semantic revision tuple certifies one retained allocation plan.  A zero
# serial means there is no committed plan and is intentionally ignored.
def duplicate_plan_violations($samples):
  [
    $samples |
    sort_by(allocation_revision_tuple) |
    group_by(allocation_revision_tuple)[] |
    . as $group |
    ([$group[].lod_allocation_plan_serial |
      select(. != null and . > 0)] | unique) as $plans |
    select($plans | length > 1) |
    {
      kind: "multiple-plans-for-revision",
      revision: ($group[0] | allocation_revision_tuple),
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

# An A/B/A observation is a nonprogress cycle when the semantic state returns
# without explicit input or a strict decrease of finite work.  Transaction and
# frame serial growth is deliberately insufficient: repeatedly drawing the
# same two states is the failure this check is intended to expose.  Identical
# repeated states are allowed because bounded worker or cache work can outlive
# the diagnostic sampling interval.
def nonprogress_cycle_violations($samples):
  [
    range(2; $samples | length) as $index |
    $samples[$index - 2] as $first |
    $samples[$index - 1] as $middle |
    $samples[$index] as $last |
    select(($first | cycle_state) == ($last | cycle_state)) |
    select(($first | cycle_state) != ($middle | cycle_state)) |
    select(($first | revision_tuple) == ($middle | revision_tuple)) |
    select(strict_progress_between($first; $last) | not) |
    select(all([$first, $middle, $last][];
      (.lod_interactive // false) == false)) |
    select(all([$first, $middle, $last][];
      (explicit_input | not))) |
    {
      kind: "nonprogress-control-cycle",
      samples: [$index - 2, $index - 1, $index],
      revision: ($first | revision_tuple),
      states: [($first | cycle_state), ($middle | cycle_state)],
      transaction_serials: [
        ($first.lod_presentation_transaction_serial // 0),
        ($middle.lod_presentation_transaction_serial // 0),
        ($last.lod_presentation_transaction_serial // 0)
      ],
      frame_serials: [
        ($first.lod_control_presented_frame_serial // 0),
        ($middle.lod_control_presented_frame_serial // 0),
        ($last.lod_control_presented_frame_serial // 0)
      ]
    }
  ];

# PRESENTATION is an executable owner, so one typed finite successor must be
# visible at the same observation boundary.  Recheck the exported source mask
# independently of the controller's aggregate violation mask: otherwise a
# defect in the runtime validator could validate itself.
def presentation_witness_violations($samples):
  [
    $samples | to_entries[] as $entry |
    select(($entry.value.lod_control_owner // 0) == 6) |
    select(($entry.value.lod_control_presentation_witness_mask // 0) == 0) |
    {
      kind: "presentation-owner-without-progress-witness",
      sample: $entry.key,
      facts: ($entry.value.lod_control_fact_mask // 0),
      host_work_flags: ($entry.value.host_work_flags // 0)
    }
  ];

# A quiet, nonterminal view must have one finite successor owner or concrete
# work which can produce one.  This is the observational counterpart of the
# terminal-composition leads-to property: the categorical phase label alone
# cannot make progress.  Exclude explicit input boundaries because the host
# may publish the input sample immediately before it installs the successor.
def ownerless_nonterminal_violations($samples):
  [
    $samples | to_entries[] as $entry |
    select(($entry.value | explicit_input) | not) |
    select(($entry.value.lod_convergence_has_state // false) == true) |
    select(($entry.value.lod_convergence_terminal // false) == false) |
    select(($entry.value.lod_interactive // false) == false) |
    select(($entry.value.lod_convergence_phase // 0) != 0) |
    select(($entry.value.lod_control_fact_mask // 0) == 0) |
    select(($entry.value.lod_control_obligation_mask // 0) == 0) |
    select(($entry.value.lod_control_owner // 0) == 0) |
    select(($entry.value | unfinished_worker_units) == 0) |
    select(($entry.value.lod_renderer_preparation_remaining_units // 0) == 0) |
    select(($entry.value.lod_renderer_preparation_preparing_targets // 0) == 0) |
    select(($entry.value.lod_convergence_source_preparation_pending // false) == false) |
    select(($entry.value.render_requested // false) == false) |
    select(($entry.value.progressive_pending // false) == false) |
    select(($entry.value.lod_results_pending // false) == false) |
    select(($entry.value.lod_submissions_pending // false) == false) |
    select(($entry.value.lod_refinement_frame_pending // false) == false) |
    {
      kind: "ownerless-nonterminal-state",
      sample: $entry.key,
      revision: ($entry.value | revision_tuple),
      phase: ($entry.value.lod_convergence_phase // 0),
      fraction: ($entry.value.lod_convergence_fraction // 0)
    }
  ];

((.samples // []) | map(select(.controller_available == true))) as $observations |
controller_samples as $samples |
(
  (if (.control_transition_trace_truncated // false) then
    [{
      kind: "truncated-control-transition-trace",
      dropped: (.control_transition_trace_dropped // 0)
    }]
   else [] end) +
  missing_field_violations($samples) +
  concrete_projection_violations($samples) +
  planning_producer_violations($samples) +
  preparation_certificate_violations($samples) +
  capacity_search_certificate_violations($samples) +
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
  presentation_witness_violations($samples) +
  ownerless_nonterminal_violations($observations) +
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
