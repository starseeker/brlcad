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
    "$jq_executable" -n --arg scenario "$1" '
        def sample($phase; $owner; $plan; $terminal):
          {
            controller_available: true,
            lod_control_obligation_mask:
              (if $terminal then 0 else 16 end),
            lod_control_owner: $owner,
            lod_control_violation_mask: 0,
            lod_inventory_revision: 1,
            lod_availability_revision: 1,
            lod_control_view_revision: 1,
            lod_view_revision: 1,
            lod_control_policy_revision: 1,
            lod_policy_revision: 1,
            lod_capacity_revision: 1,
            lod_allocation_plan_serial: $plan,
            lod_presentation_transaction_serial: 3,
            lod_presentation_required_render_serial: 0,
            lod_control_presented_frame_serial: 4,
            lod_convergence_phase: $phase,
            lod_convergence_outcome: (if $terminal then 1 else 0 end),
            lod_convergence_terminal: $terminal,
            lod_convergence_view_ready: $terminal,
            lod_interactive: false,
            action: "wait"
          };
        if $scenario == "valid" then
          {samples: [sample(1; 5; 7; false), sample(0; 0; 7; true)]}
        elif $scenario == "cycle" then
          {samples: [
            sample(1; 5; 7; false),
            sample(5; 8; 7; false),
            sample(1; 5; 7; false)
          ]}
        elif $scenario == "duplicate-plan" then
          {samples: [sample(1; 5; 7; false), sample(1; 5; 8; false)]}
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
expect_reject cycle
expect_reject duplicate-plan
expect_reject terminal-reopen
