#!/usr/bin/env bash
#
# Drive qged's actual View Select palette through Qt events.  This is an
# endpoint integration test: it intentionally does not call GED picking APIs
# directly, and runs at fractional device scale to cover logical/physical
# viewport coordinate normalization.

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_root="$(cd "$script_dir/../../.." && pwd)"
source "$script_dir/qged_test_display.sh"
build_dir="$source_root/.build"
artifact_dir=""
backend_list="system,osmesa"
run_timeout=90
scale_factor=1.5

usage()
{
    cat <<'EOF'
Usage: qged_selection_ui_matrix.sh [options]

  --build-dir DIR       Current build (default: ./.build)
  --artifact-dir DIR    New results directory (default: /tmp)
  --backends LIST       system,osmesa (default: both)
  --scale-factor VALUE  Qt device scale (default: 1.5)
  --timeout SECONDS     Per-process timeout (default: 90)
EOF
}

original_arguments=("$@")
while [[ $# -gt 0 ]]; do
    case "$1" in
	--build-dir) build_dir="$2"; shift 2 ;;
	--artifact-dir) artifact_dir="$2"; shift 2 ;;
	--backends) backend_list="$2"; shift 2 ;;
	--scale-factor) scale_factor="$2"; shift 2 ;;
	--timeout) run_timeout="$2"; shift 2 ;;
	--help|-h) usage; exit 0 ;;
	*) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

qged_test_ensure_display "$0" "${original_arguments[@]}"

build_dir="$(realpath -m "$build_dir")"
qged="$build_dir/bin/qged"
database_source="$build_dir/share/db/moss.g"
if [[ ! -f "$database_source" ]]; then
    database_source="$build_dir/regress/moss/moss.g"
fi
if [[ ! -x "$qged" || ! -f "$database_source" ]]; then
    echo "ERROR: qged and a moss.g database are required" >&2
    exit 2
fi
for required_tool in jq identify compare; do
    if ! command -v "$required_tool" >/dev/null 2>&1; then
	echo "ERROR: $required_tool is required" >&2
	exit 2
    fi
done

if [[ -z "$artifact_dir" ]]; then
    artifact_dir="/tmp/qged-selection-ui-$(date +%Y%m%d-%H%M%S)-$$"
fi
artifact_dir="$(realpath -m "$artifact_dir")"
if [[ -e "$artifact_dir" ]]; then
    echo "ERROR: artifact directory already exists: $artifact_dir" >&2
    exit 2
fi
mkdir -p "$artifact_dir"
printf 'status,backend,seconds,detail\n' > "$artifact_dir/results.csv"

canvas_target="./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas"

write_events()
{
    local file="$1"
    local images="$2"
    cat > "$file" <<EOF
{
  "schema": "brlcad.qtcad.events",
  "version": 1,
  "events": [
    {"target": ".", "action": "resize", "arguments": {"width": 1000, "height": 760, "stable_ms": 250, "timeout_ms": 5000}},
    {"target": ".", "action": "qged_command_batch", "arguments": {"commands": ["view clear", "view faceplate params 0", "in qged_select_probe.s rpp -100 100 -100 100 -100 100", "draw -m1 qged_select_probe.s", "ae 35 25", "autoview"]}},
    {"target": ".", "action": "wait_progressive_view_ready", "arguments": {"timeout_ms": 30000, "quiet_ms": 100}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$images/unselected.png"}},

    {"target": "i:org.brlcad.qged.view.select.activate", "action": "activate", "arguments": {"checked": true}},
    {"target": "i:org.brlcad.qged.view.select.add", "action": "activate", "arguments": {"checked": true}},
    {"target": "i:org.brlcad.qged.view.select.point", "action": "activate", "arguments": {"checked": true}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "$canvas_target", "action": "mouse_press", "arguments": {"x": 0.5, "y": 0.5, "button": 1, "buttons": 1, "modifiers": 0}},
    {"target": "$canvas_target", "action": "mouse_release", "arguments": {"x": 0.5, "y": 0.5, "button": 1, "buttons": 0, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 100}},
    {"target": "i:org.brlcad.qged.view.select.selection-list", "action": "assert_state", "arguments": {"count": 1}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$images/point-selected.png"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 1000}},
    {"target": "i:org.brlcad.qged.view.select.selection-list", "action": "assert_state", "arguments": {"count": 1}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$images/point-selected-stable.png"}},

    {"target": "i:org.brlcad.qged.view.select.remove", "action": "activate", "arguments": {"checked": true}},
    {"target": "$canvas_target", "action": "mouse_press", "arguments": {"x": 0.5, "y": 0.5, "button": 1, "buttons": 1, "modifiers": 0}},
    {"target": "$canvas_target", "action": "mouse_release", "arguments": {"x": 0.5, "y": 0.5, "button": 1, "buttons": 0, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 100}},
    {"target": "i:org.brlcad.qged.view.select.selection-list", "action": "assert_state", "arguments": {"count": 0}},

    {"target": "i:org.brlcad.qged.view.select.rectangle", "action": "activate", "arguments": {"checked": true}},
    {"target": "i:org.brlcad.qged.view.select.add", "action": "activate", "arguments": {"checked": true}},
    {"target": "$canvas_target", "action": "mouse_press", "arguments": {"x": 0.3, "y": 0.3, "button": 1, "buttons": 1, "modifiers": 0}},
    {"target": "$canvas_target", "action": "mouse_move", "arguments": {"x": 0.7, "y": 0.7, "button": 0, "buttons": 1, "modifiers": 0}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$images/rectangle-held.png"}},
    {"target": "$canvas_target", "action": "mouse_release", "arguments": {"x": 0.7, "y": 0.7, "button": 1, "buttons": 0, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 100}},
    {"target": "i:org.brlcad.qged.view.select.selection-list", "action": "assert_state", "arguments": {"count": 1}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$images/rectangle-selected.png"}},

    {"target": "i:org.brlcad.qged.view.select.erase-selected", "action": "activate"},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "who", "not_contains": "qged_select_probe.s"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$images/erased.png"}},
    {"target": "i:org.brlcad.qged.view.select.draw-selected", "action": "activate"},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "who", "contains": "qged_select_probe.s"}},
    {"target": ".", "action": "wait_progressive_view_ready", "arguments": {"timeout_ms": 30000, "quiet_ms": 100}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$images/redrawn.png"}}
  ]
}
EOF
}

validate_report()
{
    local report="$1"
    local images="$2"
    jq -e '
      .success == true and
      (first(.samples[] | select(.action == "wait_progressive_view_ready"))) as $stable |
	(first(.samples[] | select((.checkpoint? // "") |
	  endswith("/rectangle-selected.png")))) as $selectionEnd |
      (.samples | any(.action == "assert_state" and .selection_count == 1)) and
      (.samples | any(.action == "assert_state" and .selection_count == 0)) and
	# Selection is a semantic presentation transaction.  It may briefly own
	# one exact frame, but by the UI observation boundary it must be idle and
	# must never have upgraded that frame into renderer-capacity evidence.
	(all(.samples[];
	  if (.event_index > $stable.event_index and
	      .event_index <= $selectionEnd.event_index)
	  then ((.render_capacity_sample_requested // false) == false and
		(.lod_policy_revision == $stable.lod_policy_revision))
	  else true end)) and
	(all(.samples[] | select(.action == "assert_state");
	  (.lod_convergence_phase // -1) == 0 and
	  (.lod_convergence_view_ready // false) == true and
	  (.progressive_pending == false) and
	  (.lod_control_obligation_mask // -1) == 0 and
	  (.lod_service_pending_tasks // -1) == 0)) and
	(all(.samples[] | select(
	  .lod_convergence_semantic_presentation_frame_pending // false);
	  (.lod_progress_label_present // false) == false)) and
      (first(.samples[] | select(.action == "mouse_move" and
        .selection_rect_draw == true))) as $drag |
      ($drag.lod_gesture_active == false) and
      ($drag.lod_policy_revision == $stable.lod_policy_revision) and
      ($drag.canvas_device_pixel_ratio > 1.0) and
      ($drag.selection_rect_canvas_width > $drag.canvas_width) and
      (($drag.selection_rect_pos_x -
        ($drag.selection_rect_canvas_width * 0.3)) | fabs < 3.0) and
      ((($drag.selection_rect_pos_y + $drag.selection_rect_dim_y) -
        ($drag.selection_rect_canvas_height * 0.3)) | fabs < 3.0) and
      (($drag.selection_rect_dim_x | fabs) >
        ($drag.selection_rect_canvas_width * 0.35)) and
      (($drag.selection_rect_dim_y | fabs) >
        ($drag.selection_rect_canvas_height * 0.35)) and
      (first(.samples[] | select(.action == "mouse_release" and
        .event_index > $drag.event_index))) as $commit |
      ($commit.selection_rect_draw == false) and
      ($commit.lod_gesture_active == false) and
      ($commit.lod_policy_revision == $stable.lod_policy_revision)
    ' "$report" >/dev/null || return 1

    local name
    for name in unselected point-selected point-selected-stable rectangle-held rectangle-selected erased redrawn; do
	[[ -s "$images/$name.png" ]] || return 1
    done
    local delta
    delta="$(compare -metric AE "$images/unselected.png" \
	"$images/point-selected.png" null: 2>&1 || true)"
    delta="${delta%% *}"
    [[ "$delta" =~ ^[0-9]+$ && "$delta" -gt 0 ]] || return 1
    delta="$(compare -metric AE "$images/point-selected.png" \
	"$images/point-selected-stable.png" null: 2>&1 || true)"
    delta="${delta%% *}"
    [[ "$delta" == "0" ]] || return 1
    return 0
}

IFS=',' read -r -a backends <<< "$backend_list"
failures=0
for backend in "${backends[@]}"; do
    case "$backend" in
	system|osmesa) ;;
	*) echo "ERROR: unknown backend: $backend" >&2; exit 2 ;;
    esac
    out="$artifact_dir/$backend"
    mkdir -p "$out/images" "$out/cache"
    cp "$database_source" "$out/selection.g"
    write_events "$out/events.json" "$out/images"
    qged_args=(--test-script "$out/events.json" --test-report "$out/report.json")
    [[ "$backend" == "osmesa" ]] && qged_args=(-s "${qged_args[@]}")
    qged_args+=("$out/selection.g")
    started=$SECONDS
    if timeout --signal=TERM "$run_timeout" env \
	    QT_SCALE_FACTOR="$scale_factor" BU_DIR_CACHE="$out/cache" \
	    "$qged" "${qged_args[@]}" >"$out/stdout.log" 2>"$out/stderr.log" &&
	    validate_report "$out/report.json" "$out/images"; then
	printf 'PASS,%s,%s,\n' "$backend" "$((SECONDS - started))" \
	    >> "$artifact_dir/results.csv"
    else
	printf 'FAIL,%s,%s,selection-ui\n' "$backend" \
	    "$((SECONDS - started))" >> "$artifact_dir/results.csv"
	failures=$((failures + 1))
    fi
done

printf 'Artifacts: %s\n' "$artifact_dir"
exit "$failures"
