#!/bin/sh

set -eu

if test "$#" -ne 2; then
    echo "Usage: $0 JQ FILTER" >&2
    exit 2
fi

jq_executable="$1"
trace_filter="$2"

make_report()
{
    "$jq_executable" -n --arg scenario "$1" \
	--argjson facts "${2:-0}" \
	--argjson obligations "${3:-0}" \
	--argjson owner "${4:-0}" '
        def sample($phase; $owner; $plan; $terminal):
          {
            controller_available: true,
            lod_control_fact_mask:
              (if $terminal then 0 else 32 end),
            lod_control_obligation_mask:
              (if $terminal then 0 else 16 end),
            lod_control_owner: $owner,
            lod_control_violation_mask: 0,
            lod_control_presentation_witness_mask:
              (if $owner == 6 then 1 else 0 end),
            host_work_flags: 0,
            lod_inventory_revision: 1,
            lod_availability_revision: 1,
            lod_visibility_revision: 1,
            lod_control_view_revision: 1,
            lod_view_revision: 1,
            lod_control_policy_revision: 1,
            lod_policy_revision: 1,
            lod_capacity_revision: 1,
            lod_cad_revision: 1,
            lod_resident_demand_revision: 1,
            lod_allocation_plan_serial: $plan,
            lod_presentation_transaction_serial: 3,
            lod_presentation_required_render_serial: 0,
            lod_control_presented_frame_serial: 4,
            lod_renderer_preparation_target_signature: 0,
            lod_renderer_preparation_total_units: 0,
            lod_renderer_preparation_completed_units: 0,
            lod_renderer_preparation_remaining_units: 0,
            lod_renderer_preparation_reserved_bytes: 0,
            lod_renderer_preparation_target_count: 0,
            lod_renderer_preparation_preparing_targets: 0,
            lod_renderer_preparation_constrained_targets: 0,
            lod_renderer_preparation_failed_targets: 0,
            lod_renderer_preparation_invalid_targets: 0,
            lod_convergence_phase: $phase,
            lod_convergence_outcome: (if $terminal then 1 else 0 end),
            lod_convergence_terminal: $terminal,
            lod_convergence_view_ready: $terminal,
            lod_convergence_has_state: true,
            lod_convergence_background_pending: false,
            lod_capacity_search_phase: 0,
            lod_capacity_search_goal: 0,
            lod_capacity_search_samples_remaining: 0,
            lod_capacity_search_measured_candidates: 0,
            lod_capacity_search_total_measured_candidates: 0,
            lod_capacity_search_candidate_limit: 4,
            lod_capacity_search_maximum_candidates: 0,
            lod_capacity_search_sample_limit: 3,
            lod_capacity_search_completed_units: 0,
            lod_capacity_search_total_units: 0,
            lod_submission_source_index: 0,
            lod_submission_entry_offset: 0,
            lod_convergence_expected_leaves: 1,
            lod_convergence_available_leaves: 1,
            lod_convergence_satisfied_payloads: 1,
            lod_convergence_pending_tasks: 0,
            lod_convergence_in_flight: 0,
            lod_convergence_queued_results: 0,
            lod_convergence_source_preparation_providers: 0,
            lod_convergence_source_preparation_pending: false,
            lod_convergence_compactions: 0,
            lod_convergence_compaction_candidates: 0,
            lod_convergence_fraction: (if $terminal then 1 else 0.5 end),
            render_requested: false,
            progressive_pending: (if $terminal then false else true end),
            lod_results_pending: false,
            lod_submissions_pending: false,
            lod_refinement_frame_pending: false,
            lod_interactive: false,
            action: "wait"
          };
        if $scenario == "valid" then
          {samples: [sample(1; 5; 7; false), sample(0; 0; 7; true)]}
        elif $scenario == "cycle" then
          {samples: [
            sample(1; 5; 7; false),
            (sample(5; 8; 7; false) |
              .lod_control_fact_mask = 33554432 |
              .lod_control_obligation_mask = 128),
            sample(1; 5; 7; false)
          ]}
        elif $scenario == "alias-cycle" then
          {samples: [
            sample(1; 5; 7; false),
            (sample(1; 5; 7; false) | .lod_control_fact_mask = 64),
            sample(1; 5; 7; false)
          ]}
        elif $scenario == "frame-clock-cycle" then
          {samples: [
            sample(1; 5; 7; false),
            (sample(5; 8; 7; false) |
              .lod_control_fact_mask = 33554432 |
              .lod_control_obligation_mask = 128 |
              .lod_presentation_transaction_serial = 4 |
              .lod_control_presented_frame_serial = 5),
            (sample(1; 5; 7; false) |
              .lod_presentation_transaction_serial = 5 |
              .lod_control_presented_frame_serial = 6)
          ]}
        elif $scenario == "progress-cycle" then
          {samples: [
            sample(1; 5; 7; false),
            (sample(5; 8; 7; false) |
              .lod_control_fact_mask = 33554432 |
              .lod_control_obligation_mask = 128),
            (sample(1; 5; 7; false) |
              .lod_capacity_search_phase = 3 |
              .lod_capacity_search_samples_remaining = 3 |
              .lod_capacity_search_measured_candidates = 1 |
              .lod_capacity_search_total_measured_candidates = 1 |
              .lod_capacity_search_maximum_candidates = 4 |
              .lod_capacity_search_completed_units = 7 |
              .lod_capacity_search_total_units = 20)
          ]}
        elif $scenario == "preparation-progress-cycle" then
          {samples: [
            (sample(1; 5; 7; false) |
              .lod_renderer_preparation_target_signature = 19 |
              .lod_renderer_preparation_total_units = 10 |
              .lod_renderer_preparation_completed_units = 2 |
              .lod_renderer_preparation_remaining_units = 8 |
              .lod_renderer_preparation_target_count = 1 |
              .lod_renderer_preparation_preparing_targets = 1),
            (sample(5; 8; 7; false) |
              .lod_control_fact_mask = 33554432 |
              .lod_control_obligation_mask = 128 |
              .lod_renderer_preparation_target_signature = 19 |
              .lod_renderer_preparation_total_units = 10 |
              .lod_renderer_preparation_completed_units = 4 |
              .lod_renderer_preparation_remaining_units = 6 |
              .lod_renderer_preparation_target_count = 1 |
              .lod_renderer_preparation_preparing_targets = 1),
            (sample(1; 5; 7; false) |
              .lod_renderer_preparation_target_signature = 19 |
              .lod_renderer_preparation_total_units = 10 |
              .lod_renderer_preparation_completed_units = 6 |
              .lod_renderer_preparation_remaining_units = 4 |
              .lod_renderer_preparation_target_count = 1 |
              .lod_renderer_preparation_preparing_targets = 1)
          ]}
        elif $scenario == "dense-cycle" then
          {
            samples: [sample(0; 0; 7; true)],
            control_transitions: [
              sample(1; 5; 7; false),
              (sample(5; 8; 7; false) |
                .lod_control_fact_mask = 33554432 |
                .lod_control_obligation_mask = 128),
              sample(1; 5; 7; false)
            ]
          }
        elif $scenario == "truncated" then
          {
            samples: [sample(0; 0; 7; true)],
            control_transition_trace_truncated: true,
            control_transition_trace_dropped: 1
          }
        elif $scenario == "fact-obligation-mismatch" then
          {samples: [
            (sample(1; 5; 7; false) |
              .lod_control_obligation_mask = 32)
          ]}
        elif $scenario == "fact-owner-mismatch" then
          {samples: [(sample(1; 6; 7; false))]}
        elif $scenario == "missing-fact-mask" then
          {samples: [(sample(1; 5; 7; false) |
            del(.lod_control_fact_mask))]}
        elif $scenario == "projection" then
          {samples: [(sample(1; $owner; 7; false) |
            .lod_control_fact_mask = $facts |
            .lod_control_obligation_mask = $obligations)]}
        elif $scenario == "duplicate-plan" then
          {samples: [sample(1; 5; 7; false), sample(1; 5; 8; false)]}
        elif $scenario == "bookkeeping-masked-duplicate-plan" then
          {samples: [
            sample(1; 5; 7; false),
            (sample(1; 5; 8; false) |
              .lod_cad_revision = 2 |
              .lod_resident_demand_revision = 2)
          ]}
        elif $scenario == "ownerless-nonterminal" then
          {samples: [(sample(3; 0; 0; false) |
            .lod_control_fact_mask = 0 |
            .lod_control_obligation_mask = 0 |
            .progressive_pending = false)]}
        elif $scenario == "unwitnessed-presentation" then
          {samples: [(sample(4; 6; 7; false) |
            .lod_control_fact_mask = 262144 |
            .lod_control_obligation_mask = 32 |
            .lod_control_presentation_witness_mask = 0)]}
        elif $scenario == "preparation-regression" then
          {samples: [
            (sample(1; 5; 7; false) |
              .lod_renderer_preparation_target_signature = 23 |
              .lod_renderer_preparation_total_units = 10 |
              .lod_renderer_preparation_completed_units = 6 |
              .lod_renderer_preparation_remaining_units = 4 |
              .lod_renderer_preparation_target_count = 1 |
              .lod_renderer_preparation_preparing_targets = 1),
            (sample(1; 5; 7; false) |
              .lod_renderer_preparation_target_signature = 23 |
              .lod_renderer_preparation_total_units = 10 |
              .lod_renderer_preparation_completed_units = 5 |
              .lod_renderer_preparation_remaining_units = 5 |
              .lod_renderer_preparation_target_count = 1 |
              .lod_renderer_preparation_preparing_targets = 1)
          ]}
        elif $scenario == "invalid-preparation" then
          {samples: [(sample(1; 5; 7; false) |
            .lod_renderer_preparation_target_signature = 29 |
            .lod_renderer_preparation_total_units = 10 |
            .lod_renderer_preparation_completed_units = 4 |
            .lod_renderer_preparation_remaining_units = 7 |
            .lod_renderer_preparation_target_count = 1 |
            .lod_renderer_preparation_preparing_targets = 1)]}
        elif $scenario == "invalid-capacity-rank" then
          {samples: [(sample(1; 5; 7; false) |
            .lod_capacity_search_phase = 3 |
            .lod_capacity_search_samples_remaining = 3 |
            .lod_capacity_search_maximum_candidates = 4 |
              .lod_capacity_search_completed_units = 1 |
              .lod_capacity_search_total_units = 20)]}
        elif $scenario == "orphaned-selective-scope" then
          {samples: [(sample(3; 5; 7; false) |
            .lod_control_fact_mask = 128 |
            .lod_control_obligation_mask = 16)]}
        elif $scenario == "terminal-reopen" then
          {samples: [sample(0; 0; 7; true), sample(1; 5; 7; false)]}
        elif $scenario == "input-reopen" then
          {samples: [
            sample(0; 0; 7; true),
            (sample(1; 5; 7; false) | .action = "qged_command")
          ]}
        else
          error("unknown scenario")
        end
    '
}

expect_accept()
{
    if ! make_report "$1" | "$jq_executable" -e -f "$trace_filter" \
	    >/dev/null 2>&1; then
	echo "control trace unexpectedly rejected: $1" >&2
	exit 1
    fi
}

expect_reject()
{
    if make_report "$1" | "$jq_executable" -e -f "$trace_filter" \
	    >/dev/null 2>&1; then
	echo "control trace unexpectedly accepted: $1" >&2
	exit 1
    fi
}

expect_accept valid
expect_accept input-reopen
expect_accept progress-cycle
expect_accept preparation-progress-cycle
expect_reject cycle
expect_reject alias-cycle
expect_reject frame-clock-cycle
expect_reject dense-cycle
expect_reject truncated
expect_reject fact-obligation-mismatch
expect_reject fact-owner-mismatch
expect_reject missing-fact-mask
expect_reject duplicate-plan
expect_reject bookkeeping-masked-duplicate-plan
expect_reject ownerless-nonterminal
expect_reject unwitnessed-presentation
expect_reject preparation-regression
expect_reject invalid-preparation
expect_reject invalid-capacity-rank
expect_reject orphaned-selective-scope
expect_reject terminal-reopen

# Every concrete fact and every owner-precedence boundary is an independent
# part of the C++/TLA+ refinement map.  Keep this table explicit so changing a
# numeric public diagnostic value cannot silently change the offline checker.
while read -r facts obligations owner; do
    if ! make_report projection "$facts" "$obligations" "$owner" |
	    "$jq_executable" -e -f "$trace_filter" >/dev/null 2>&1; then
	echo "control trace unexpectedly rejected projection: $facts" >&2
	exit 1
    fi
done <<'EOF'
1 1 1
2 2 2
4 4 3
8 8 4
16 8 4
32 16 5
64 16 5
160 16 5
256 16 5
512 16 5
1024 16 5
2048 16 5
4096 16 5
8192 16 5
16384 16 5
32768 16 5
65536 16 5
131072 32 6
262144 32 6
524288 32 6
1048576 32 6
2097152 32 6
4194304 32 6
8388608 32 6
16777216 64 7
33554432 128 8
67108864 256 9
134217728 16 5
268435456 32 6
131073 33 6
9 9 1
10 10 4
6 6 2
262148 36 3
262176 48 6
270336 48 5
16777248 80 5
50331648 192 7
100663296 384 8
536870911 511 6
EOF
