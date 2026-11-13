#!/bin/sh

set -eu

if test "$#" -ne 4; then
    echo "Usage: $0 JQ IDENTIFY CONVERT CHECKER" >&2
    exit 2
fi

jq_executable="$1"
identify_executable="$2"
convert_executable="$3"
checker="$4"
test_dir=$(mktemp -d "${TMPDIR:-/tmp}/qged-hud-contract.XXXXXX")
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM

# Exact rectangles isolate the visibility contract from GL and font timing.
"$convert_executable" -size 240x200 xc:black -fill '#404040' \
    -draw 'rectangle 229,30 235,109' "$test_dir/hidden.png"
"$convert_executable" "$test_dir/hidden.png" -fill '#ffcd48' \
    -draw 'rectangle 229,30 235,42' "$test_dir/short.png"
"$convert_executable" "$test_dir/hidden.png" -fill '#ffcd48' \
    -draw 'rectangle 229,30 235,79' "$test_dir/visible.png"
"$convert_executable" "$test_dir/hidden.png" -fill '#ff5a50' \
    -draw 'rectangle 229,30 235,109' "$test_dir/error.png"

"$jq_executable" -n --arg directory "$test_dir" '
    def sample($image; $height):
        {checkpoint: ($directory + "/" + $image + ".png"),
         presented_cad_faces: 1000, lod_convergence_visible_targets: 1,
         lod_convergence_phase: 3, render_requested: false,
         cad_aggregate_work_exact: true,
         feature_presentation_revision: "10", presented_feature_revision: "10",
         lod_progress_estimate_available: true,
         lod_estimated_fraction: 0.15, lod_progress_fill_present: true,
         lod_progress_fill_minimum_y: 30,
         lod_progress_fill_maximum_y: (30 + $height),
         lod_progress_fill_line_width: 7};
    {samples: [sample("short"; 13), sample("visible"; 50)]}
' > "$test_dir/base.json"

run_case()
{
    scenario="$1"
    mutation="$2"
    expected="$3"
    "$jq_executable" "$mutation" "$test_dir/base.json" > "$test_dir/input.json"
    if sh "$checker" "$jq_executable" "$identify_executable" \
        "$convert_executable" "$test_dir/input.json" \
        > "$test_dir/result.txt" 2>&1; then
        if test "$expected" = fail; then
            echo "ERROR: $scenario unexpectedly passed" >&2
            exit 1
        fi
        if ! "$jq_executable" -e -Rs --arg expected "$expected" \
            'contains($expected)' "$test_dir/result.txt" > /dev/null; then
            echo "ERROR: $scenario did not report $expected" >&2
            cat "$test_dir/result.txt" >&2
            exit 1
        fi
    elif test "$expected" != fail; then
        echo "ERROR: $scenario unexpectedly failed" >&2
        cat "$test_dir/result.txt" >&2
        exit 1
    fi
}

run_case lagged_estimate '.' '350 pixels'
run_case zero_new_estimate '.samples[].lod_estimated_fraction = 0' '350 pixels'
run_case no_large_fill '.samples |= .[:1]' 'SKIP:'
run_case no_active_hud '.samples[].lod_convergence_phase = 0' 'SKIP:'
run_case hidden_fill '.samples[1].checkpoint |= sub("visible"; "hidden")' fail
run_case insufficient_pixels '.samples[1].checkpoint |= sub("visible"; "short")' fail
run_case no_pixel_fallback '.samples += [.samples[1]] |
    .samples[1].checkpoint |= sub("visible"; "hidden")' fail
run_case later_hidden_fill '.samples += [(.samples[1] |
    .checkpoint |= sub("visible"; "hidden"))]' fail
run_case missing_image '.samples[1].checkpoint += ".missing"' fail
run_case absent_geometry 'del(.samples[1].lod_progress_fill_minimum_y)' fail
run_case inverted_geometry '.samples[1].lod_progress_fill_maximum_y = 20' fail
run_case invalid_width '.samples[1].lod_progress_fill_line_width = 0' fail

run_case pending_unpresented_hud '.samples =
    [(.samples[1] | .checkpoint |= sub("visible"; "hidden") |
        .render_requested = true)] + .samples' '350 pixels'
run_case only_pending_hud '.samples[].render_requested = true' 'SKIP:'
run_case incomplete_frame '.samples[].cad_aggregate_work_exact = false' 'SKIP:'
run_case missing_frame_state 'del(.samples[1].render_requested)' fail
run_case missing_frame_execution 'del(.samples[1].cad_aggregate_work_exact)' fail
run_case invalid_frame_state '.samples[1].render_requested = "false"' fail
run_case interactive_frame '.samples[].lod_convergence_phase = 2' '350 pixels'
run_case wire_frame '.samples[] |=
    (.presented_cad_faces = 0 | .presented_cad_lines = 1000)' '350 pixels'
run_case no_presented_geometry '.samples[].presented_cad_faces = 0' 'SKIP:'
run_case retained_interrupted_frame '.samples[] |=
    (.checkpoint |= sub("visible"; "hidden") |
     .presented_feature_revision = "9")' 'SKIP:'
run_case provisional_frame '.samples[].presented_feature_revision = null' 'SKIP:'
run_case missing_presented_identity 'del(.samples[1].presented_feature_revision)' fail
run_case missing_current_identity 'del(.samples[1].feature_presentation_revision)' fail
run_case invalid_presented_identity '.samples[1].presented_feature_revision = false' fail
run_case invalid_current_identity '.samples[1].feature_presentation_revision = "invalid"' fail
run_case resumed_hidden_fill '.samples += [(.samples[1] |
    .checkpoint |= sub("visible"; "hidden"))] |
    .samples[1].presented_feature_revision = "9"' fail

"$jq_executable" '.samples = [(.samples[1] |
    .checkpoint |= sub("visible"; "error") |
    .lod_convergence_phase = 6 | .lod_convergence_terminal_error = true |
    .view_lod_policy = 1 | .view_lod_mesh_enabled = true |
    .lod_progress_fill_color = [255, 90, 80] |
    .lod_progress_fill_maximum_y = 110 |
    .lod_progress_label_present = true |
    .lod_progress_label_text = "View incomplete  1 geometry error" |
    .lod_progress_track_present = true |
    .lod_progress_track_minimum_y = 110 |
    .lod_progress_track_maximum_y = 110)]' "$test_dir/base.json" \
    > "$test_dir/error.json"
mv "$test_dir/error.json" "$test_dir/base.json"
run_case terminal_error '.' '560 pixels'
run_case missing_terminal_fill '.samples[].lod_progress_fill_present = false' fail
run_case unfinished_terminal_fill '.samples[].lod_progress_track_minimum_y = 30' fail
run_case missing_error_label '.samples[].lod_progress_label_present = false' fail
run_case wrong_error_color '.samples[].lod_progress_fill_color = [96, 220, 255]' fail
run_case hidden_terminal_fill '.samples[].checkpoint |= sub("error"; "hidden")' fail
run_case unpresented_terminal_fill '.samples[].presented_feature_revision = "9"' fail
run_case pending_terminal_fill '.samples[] |=
    (.render_requested = true | .lod_progress_fill_present = false)' 'SKIP:'
run_case disabled_error_hud '.samples[] |=
    (.view_lod_policy = 0 | .lod_progress_fill_present = false)' 'SKIP:'

echo "HUD captured-geometry and hidden-fill contract passed"
