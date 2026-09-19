#!/bin/sh

set -eu

if test "$#" -ne 2; then
    echo "Usage: $0 JQ FILTER" >&2
    exit 2
fi

jq_executable="$1"
lucy_filter="$2"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/qged-lucy-contract.XXXXXX")
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM

"$jq_executable" -n '
  def sample($index; $name):
    {event_index: $index, checkpoint: ("images/" + $name + ".png"),
     active_lod_cad_payloads: 1, active_progressive_cad_faces: 10000,
     presented_cad_work_exact: true, presented_cad_faces: 10000,
     presented_cad_lines: 15000, visible_structural_fallback_boxes: 0,
     active_cad_subpixel_proxy_points: 0, active_cad_subpixel_proxy_draw_points: 0,
     deep_lod_diagnostics: true, compact_lod_entries: 1,
     compact_lod_entries_with_payload: 1, lod_service_resident_assets: 1,
     lod_service_resident_bytes: 100000, lod_convergence_view_ready: true,
     lod_prominent_cad_payloads: 1, lod_prominent_cad_quality_floor_violations: 0,
     lod_convergence_performance_limited: false, lod_max_cad_normalized_error: 1};
  {samples: ([sample(0; "draw-return") |
      .active_lod_cad_payloads = 0 | .active_progressive_cad_faces = 0 |
      .presented_cad_faces = 0 | .presented_cad_lines = 0 |
      .active_cad_subpixel_proxy_points = 1 |
      .active_cad_subpixel_proxy_draw_points = 1 |
      .lod_convergence_view_ready = false] +
    (["first-coverage-ready", "first-cad-mesh-ready", "zoom-in-stable",
      "smooth-zoom-close-stable", "smooth-zoom-out-stable", "smooth-zoom-return",
      "background-cache-complete", "final-stable"] | to_entries |
      map(sample(.key + 1; .value))))}
' > "$test_dir/base.json"

run_case()
{
    scenario="$1"
    mutation="$2"
    expected="$3"
    mode="${4:-shaded}"
    "$jq_executable" "$mutation" "$test_dir/base.json" > "$test_dir/input.json"
    if "$jq_executable" -e --arg mode "$mode" -f "$lucy_filter" \
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

run_case measured '.' pass
run_case coverage_before_frame '.samples[1].presented_cad_work_exact = false' pass
run_case unmeasured '.samples[] |= (.deep_lod_diagnostics = false |
    .compact_lod_entries = null | .compact_lod_entries_with_payload = null)' pass
run_case historical_perf '.samples[] |= (.deep_lod_diagnostics = false |
    .compact_lod_entries = 0 | .compact_lod_entries_with_payload = 0)' pass
run_case measured_empty '.samples[-1].compact_lod_entries = 0' registry.entries
run_case measured_missing 'del(.samples[-1].compact_lod_entries_with_payload)' registry.payload_entries
run_case missing_availability 'del(.samples[-1].deep_lod_diagnostics)' registry.availability
run_case missing_coverage '.samples |= map(select(.event_index != 1))' coverage.checkpoint
run_case no_mesh '.samples[].presented_cad_faces = 0' first_mesh.presentation
run_case no_residency '.samples[-1].lod_service_resident_bytes = 0' final.resident_bytes
run_case early_false_ready '.samples[0].lod_convergence_view_ready = true' presentation.no_proxy_after_coverage
run_case later_proxy '.samples[2] |= (.active_cad_subpixel_proxy_points = 1 |
    .lod_convergence_view_ready = false)' presentation.no_proxy_after_coverage
run_case missing_zoom '.samples |= map(select(.event_index != 3))' quality.checkpoint.zoom-in-stable
run_case constrained_bad_zoom '.samples[3] |= (.lod_prominent_cad_quality_floor_violations = 1 |
    .lod_convergence_performance_limited = true)' quality.prominent_floor
run_case unsettled_bad_zoom '.samples[4] |= (.lod_prominent_cad_quality_floor_violations = 1 |
    .lod_convergence_view_ready = false)' quality.prominent_floor
run_case missing_floor 'del(.samples[5].lod_prominent_cad_quality_floor_violations)' quality.prominent_floor
run_case pixel_error '.samples[7].lod_max_cad_normalized_error = 2' settled.pixel_target
run_case missing_pixel_error 'del(.samples[7].lod_max_cad_normalized_error)' settled.pixel_target
run_case wire '.samples[].presented_cad_faces = 0' pass wire
run_case empty_wire '.samples[].presented_cad_lines = 0' first_mesh.presentation wire
run_case wire_floor '.samples[5].lod_prominent_cad_quality_floor_violations = 1' quality.prominent_floor wire

# A composite failure must preserve each independent diagnosis.
run_case multiple '.samples[-1].compact_lod_entries = 0 |
    .samples[3].lod_prominent_cad_quality_floor_violations = 1' registry.entries
"$jq_executable" -e -s '
    any(.[0].failures[]; .condition == "quality.prominent_floor")
' "$test_dir/result.json" > /dev/null

echo "Lucy observation, startup coverage, and settled-quality contract passed"
