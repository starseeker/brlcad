#!/bin/sh

set -eu

if test "$#" -ne 2; then
    echo "Usage: $0 JQ FILTER" >&2
    exit 2
fi

jq_executable="$1"
pose_filter="$2"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/qged-pose-contract.XXXXXX")
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM

"$jq_executable" -n '
  def sample($index; $action; $name):
    {event_index: $index, action: $action,
     checkpoint: (if $name == "" then null else "images/" + $name + ".png" end),
     last_render_ms: 30, lod_interactive_target_fps: 20,
     lod_target_pixel_error: 1, lod_interactive_progressive_ceiling: -1,
     lod_scale_changing_interaction: false,
     active_progressive_cad_faces: 0, presented_cad_work_exact: true,
     presented_cad_faces: 100, presented_cad_lines: 150,
     presentation_interrupted_frames: 3, presentation_last_interrupted_ms: 50.14,
     presentation_deadline_interactive_ms: 50};
  {backend: "system_gl", samples: [
    (sample(0; "wheel"; "") | .lod_scale_changing_interaction = true),
    (sample(1; "checkpoint"; "zoom-out-motion") | .lod_scale_changing_interaction = true),
    sample(2; "checkpoint"; "zoom-return-stable"),
    sample(3; "mouse_press"; ""), sample(4; "wait"; ""),
    sample(5; "checkpoint"; "rotate-held-end"), sample(6; "mouse_release"; ""),
    sample(7; "checkpoint"; "rotate-motion")]}
' > "$test_dir/base.json"

run_case()
{
    scenario="$1"
    mutation="$2"
    expected="$3"
    mode="${4:-shaded}"
    "$jq_executable" "$mutation" "$test_dir/base.json" > "$test_dir/input.json"
    if "$jq_executable" -e --arg mode "$mode" -f "$pose_filter" \
	"$test_dir/input.json" > "$test_dir/result.json"; then
	if test "$expected" != pass; then
	    echo "ERROR: $scenario unexpectedly passed" >&2
	    exit 1
	fi
    elif test "$expected" = pass; then
	echo "ERROR: $scenario unexpectedly failed" >&2
	cat "$test_dir/result.json" >&2
	exit 1
    elif ! "$jq_executable" -e -s --arg condition "$expected" '
	.[-1] == false and any(.[0].failures[]; .condition == $condition)
	' "$test_dir/result.json" > /dev/null; then
	echo "ERROR: $scenario did not identify $expected" >&2
	cat "$test_dir/result.json" >&2
	exit 1
    fi
}

run_case responsive '.' pass
run_case historical_abort_cannot_justify_ceiling \
    '.samples[5].lod_interactive_progressive_ceiling = 0' responsive_ceiling
run_case responsive_quality_loss '.samples[5].lod_target_pixel_error = 2' responsive_pixel_error
run_case responsive_face_loss '.samples[5].presented_cad_faces = 94' responsive_population
run_case retention_boundary '.samples[5].presented_cad_faces = 95' pass
run_case empty_baseline '.samples[2].presented_cad_faces = 0' responsive_population
run_case missing_population 'del(.samples[5].presented_cad_faces)' responsive_population
run_case inexact_frame '.samples[5].presented_cad_work_exact = false' responsive_population
run_case wire '.samples[].presented_cad_faces = 0' pass wire
run_case responsive_wire_loss '.samples[].presented_cad_faces = 0 |
    .samples[5].presented_cad_lines = 0' responsive_population wire
run_case completed_miss '.samples[4].last_render_ms = 53 |
    .samples[5].lod_interactive_progressive_ceiling = 0' pass
run_case uncorrected_completed_miss '.samples[4].last_render_ms = 53' pressure_response
run_case completed_slack '.samples[4].last_render_ms = 52.5' pass
run_case fresh_abort '.samples[4:] |= map(.presentation_interrupted_frames = 4) |
    .samples[5].lod_interactive_progressive_ceiling = 0' pass
run_case uncorrected_abort '.samples[4:] |= map(.presentation_interrupted_frames = 4)' pressure_response
run_case abort_after_held '.samples[6:] |= map(.presentation_interrupted_frames = 4) |
    .samples[5].lod_interactive_progressive_ceiling = 0' responsive_ceiling
run_case slow_frame_after_held '.samples[6].last_render_ms = 60 |
    .samples[5].lod_interactive_progressive_ceiling = 0' responsive_ceiling
run_case missing_counter 'del(.samples[4].presentation_interrupted_frames)' interruption_counters
run_case decreasing_counter '.samples[4].presentation_interrupted_frames = 2' interruption_counters
run_case fractional_counter '.samples[4].presentation_interrupted_frames = 3.5' interruption_counters
run_case missing_deadline '.samples[4:] |= map(.presentation_interrupted_frames = 4) |
    del(.samples[4].presentation_deadline_interactive_ms)' interruption_deadlines
run_case premature_abort '.samples[4:] |= map(.presentation_interrupted_frames = 4) |
    .samples[4].presentation_last_interrupted_ms = 49' interruption_deadlines
run_case missing_timing 'del(.samples[4].last_render_ms)' timing
run_case invalid_fps '.samples[2].lod_interactive_target_fps = 0' timing
run_case missing_held '.samples |= map(select(.event_index != 5))' checkpoints
run_case pose_changes_scale '.samples[5].lod_scale_changing_interaction = true' pose_scale
run_case zoom_not_scale '.samples[0].lod_scale_changing_interaction = false' zoom_scale
run_case software '.backend = "osmesa" | .samples[5].lod_interactive_progressive_ceiling = 0' pass

echo "Pose completed-frame, interruption, and retained-population contract passed (28 cases)"
