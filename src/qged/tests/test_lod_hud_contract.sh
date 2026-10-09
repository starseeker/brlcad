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

"$convert_executable" -size 512x200 xc:black -fill '#101820' \
    -draw 'rectangle 62,162 503,191' "$test_dir/hidden.png"
"$convert_executable" "$test_dir/hidden.png" -fill '#ffbf47' \
    -draw 'rectangle 84,180 94,182' "$test_dir/short.png"
"$convert_executable" "$test_dir/hidden.png" -fill '#ffbf47' \
    -draw 'rectangle 84,162 85,191' "$test_dir/visible.png"
"$convert_executable" "$test_dir/hidden.png" -fill '#ff594f' \
    -draw 'rectangle 84,162 85,191' "$test_dir/error.png"

"$jq_executable" -n --arg directory "$test_dir" '
    def sample($image):
        {checkpoint: ($directory + "/" + $image + ".png"),
         presented_cad_faces: 1000, lod_convergence_visible_targets: 1,
         lod_convergence_phase: 3, render_requested: false,
         cad_aggregate_work_exact: true,
         feature_presentation_revision: "10", presented_feature_revision: "10",
         lod_progress_card_present: true,
         lod_progress_card_title: "Refining visible detail",
         lod_progress_card_detail: "6 geometry tasks running",
         lod_progress_card_color: [255, 191, 71],
         lod_progress_card_terminal: false,
         lod_progress_card_ready: false};
    {samples: [sample("short"), sample("visible")]}
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

run_case visible_card '.samples |= .[1:]' '60 pixels'
run_case no_active_hud '.samples[].lod_convergence_phase = 0' 'SKIP:'
run_case no_published_card '.samples[].lod_progress_card_present = false' 'SKIP:'
run_case hidden_card '.samples[1].checkpoint |= sub("visible"; "hidden")' fail
run_case insufficient_accent '.samples |= .[:1]' fail
run_case missing_image '.samples[1].checkpoint += ".missing"' fail
run_case pending_unpresented_card '.samples =
    [(.samples[1] | .checkpoint |= sub("visible"; "hidden") |
        .render_requested = true)] + .samples' '60 pixels'
run_case only_pending_card '.samples[].render_requested = true' 'SKIP:'
run_case incomplete_frame '.samples[].cad_aggregate_work_exact = false' 'SKIP:'
run_case missing_frame_state 'del(.samples[1].render_requested)' fail
run_case missing_frame_execution 'del(.samples[1].cad_aggregate_work_exact)' fail
run_case retained_interrupted_frame '.samples[] |=
    (.checkpoint |= sub("visible"; "hidden") |
     .presented_feature_revision = "9")' 'SKIP:'
run_case missing_presented_identity 'del(.samples[1].presented_feature_revision)' fail
run_case invalid_current_identity '.samples[1].feature_presentation_revision = "invalid"' fail

"$jq_executable" '.samples = [(.samples[1] |
    .checkpoint |= sub("visible"; "error") |
    .lod_convergence_phase = 6 | .lod_convergence_terminal_error = true |
    .view_lod_policy = 1 | .view_lod_mesh_enabled = true |
    .lod_progress_card_color = [255, 89, 79] |
    .lod_progress_card_title = "Geometry preparation failed" |
    .lod_progress_card_terminal = true)]' "$test_dir/base.json" \
    > "$test_dir/error.json"
mv "$test_dir/error.json" "$test_dir/base.json"
run_case terminal_error '.' '60 pixels'
run_case missing_terminal_card '.samples[].lod_progress_card_present = false' fail
run_case wrong_error_title '.samples[].lod_progress_card_title = "View ready"' fail
run_case wrong_error_color '.samples[].lod_progress_card_color = [97, 191, 255]' fail
run_case active_error_card '.samples[].lod_progress_card_terminal = false' fail
run_case hidden_terminal_card '.samples[].checkpoint |= sub("error"; "hidden")' fail
run_case unpresented_terminal_card '.samples[].presented_feature_revision = "9"' fail
run_case pending_terminal_card '.samples[] |=
    (.render_requested = true | .lod_progress_card_present = false)' 'SKIP:'
run_case disabled_error_card '.samples[] |=
    (.view_lod_policy = 0 | .lod_progress_card_present = false)' 'SKIP:'

echo "LoD HUD contract tests passed"
