#!/usr/bin/env bash
#
# Exercise qged through its real Qt canvas and write reproducible screenshots,
# LoD/timing reports, logs, cache inventories, and optional perf/apitrace data.
# Cold and warm are always a pair: the cold run starts with a newly created,
# empty BU_DIR_CACHE (no format or data children), and the warm run reuses
# exactly what that run produced.

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_root="$(cd "$script_dir/../../.." && pwd)"
build_dir="$source_root/.build"
main_build_dir="/home/cyapp/brlcad/.build"
artifact_dir=""
profile="full"
case_list=""
backend_list="system,osmesa"
mode_list="shaded,wire"
swap_list="default"
run_baseline=1
perf_case=""
apitrace_case=""
run_timeout=180

usage()
{
    cat <<'EOF'
Usage: qged_gui_matrix.sh [options]

  --build-dir DIR          Current build (default: ./.build)
  --main-build-dir DIR     Production baseline build
  --artifact-dir DIR       New results directory (never cleared)
  --profile smoke|full|stress
  --cases LIST             Comma-separated case names
  --backends LIST          system,osmesa (default: both)
  --modes LIST             shaded,wire (default: both)
  --swap-intervals LIST    default,0,1,-1 (default: default)
  --no-baseline            Do not capture the production qged baseline
  --perf CASE              Record one cold System GL case with perf
  --apitrace CASE          Trace one cold System GL case with apitrace
  --timeout SECONDS        Per-process timeout (default: 180)

Cases: generic_twin, lucy, multi_lucy, multi_lucy_xpush, stanford, havoc,
       hubble

The smoke profile uses generic_twin and lucy.  Full adds havoc and Hubble.
Stress adds the shared/expanded Lucy scenes and the combined Stanford scene.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
	--build-dir) build_dir="$2"; shift 2 ;;
	--main-build-dir) main_build_dir="$2"; shift 2 ;;
	--artifact-dir) artifact_dir="$2"; shift 2 ;;
	--profile) profile="$2"; shift 2 ;;
	--cases) case_list="$2"; shift 2 ;;
	--backends) backend_list="$2"; shift 2 ;;
	--modes) mode_list="$2"; shift 2 ;;
	--swap-intervals) swap_list="$2"; shift 2 ;;
	--no-baseline) run_baseline=0; shift ;;
	--perf) perf_case="$2"; shift 2 ;;
	--apitrace) apitrace_case="$2"; shift 2 ;;
	--timeout) run_timeout="$2"; shift 2 ;;
	--help|-h) usage; exit 0 ;;
	*) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

case "$profile" in
    smoke) default_cases="generic_twin,lucy" ;;
    full) default_cases="generic_twin,lucy,havoc,hubble" ;;
    stress)
	default_cases="generic_twin,lucy,multi_lucy,multi_lucy_xpush,stanford,havoc,hubble"
	;;
    *) echo "ERROR: unknown profile '$profile'" >&2; exit 2 ;;
esac
[[ -n "$case_list" ]] || case_list="$default_cases"

build_dir="$(realpath -m "$build_dir")"
main_build_dir="$(realpath -m "$main_build_dir")"
if [[ -z "$artifact_dir" ]]; then
    artifact_dir="$build_dir/qged-gui-matrix/$(date +%Y%m%d-%H%M%S)"
fi
artifact_dir="$(realpath -m "$artifact_dir")"

qged="$build_dir/bin/qged"
gsh="$build_dir/bin/gsh"
baseline_qged="$main_build_dir/bin/qged"
if [[ ! -x "$qged" || ! -x "$gsh" ]]; then
    echo "ERROR: qged and gsh are required in $build_dir/bin" >&2
    exit 2
fi
for required_tool in jq identify convert compare; do
    if ! command -v "$required_tool" >/dev/null 2>&1; then
	echo "ERROR: $required_tool is required for GUI report validation" >&2
	exit 2
    fi
done
if [[ -e "$artifact_dir" ]]; then
    echo "ERROR: artifact directory already exists: $artifact_dir" >&2
    echo "Choose a new directory; the runner never clears prior evidence." >&2
    exit 2
fi
mkdir -p "$artifact_dir"/{cases,caches,events,inventory,baseline}

case_spec()
{
    local case_name="$1"
    local generic="$build_dir/Generic_Twin.g"
    if [[ ! -f "$generic" ]] || [[ "$(stat -c %s "$generic")" -lt 1000 ]]; then
	generic="$build_dir/share/db/faa/Generic_Twin.g"
    fi
    case "$case_name" in
	generic_twin)
	    printf '%s|%s|%s|%s|%s\n' "$generic" "all" \
		"all" "0xxx_series" "all/0xxx_series"
	    ;;
	lucy)
	    printf '%s|%s|%s|%s|%s\n' "$build_dir/lucy.g" "lucy.bot.r" \
		"lucy.bot.r" "lucy.bot" "lucy.bot.r/lucy.bot"
	    ;;
	multi_lucy)
	    printf '%s|%s|||\n' "$build_dir/stanford.g" "multi_lucy"
	    ;;
	multi_lucy_xpush)
	    printf '%s|%s|||\n' "$build_dir/stanford.g" "multi_lucy_xpush"
	    ;;
	stanford)
	    printf '%s|%s|%s|%s|%s\n' "$build_dir/stanford.g" "all" \
		"all" "Armadillo.bot.r" "all/Armadillo.bot.r"
	    ;;
	havoc)
	    printf '%s|%s|%s|%s|%s\n' "$build_dir/share/db/havoc.g" "havoc" \
		"havoc" "havoc_front" "havoc/havoc_front"
	    ;;
	hubble)
	    printf '%s|%s|%s|%s|%s\n' \
		"/home/cyapp/models/NASA/Hubble/Hubble_Space_Telescope.g" \
		"all.g" "all.g" "SA4097" "all.g/SA4097"
	    ;;
	*) return 1 ;;
    esac
}

contains_csv()
{
    local values=",$1,"
    [[ "$values" == *",$2,"* ]]
}

write_inventory()
{
    local case_name="$1"
    local db="$2"
    local object="$3"
    local hierarchy_root="$4"
    local hierarchy_child="$5"
    local hierarchy_path="$6"
    local out="$artifact_dir/inventory/$case_name.txt"
    {
	printf 'case=%s\n' "$case_name"
	printf 'database=%s\n' "$db"
	printf 'object=%s\n' "$object"
	printf 'hierarchy_root=%s\n' "$hierarchy_root"
	printf 'hierarchy_child=%s\n' "$hierarchy_child"
	printf 'hierarchy_path=%s\n' "$hierarchy_path"
	stat -c $'size=%s\nmtime=%Y\nmode=%a' "$db"
	printf 'first_mib_sha256='
	head -c 1048576 "$db" | sha256sum | cut -d' ' -f1
	printf 'tops:\n'
	"$gsh" "$db" tops 2>&1
    } > "$out"
}

canvas_target="./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas"

write_event_script()
{
    local output="$1"
    local image_dir="$2"
    local mode="$3"
    local object="$4"
    local settle_ms="$5"
    local hierarchy_root="$6"
    local hierarchy_child="$7"
    local hierarchy_path="$8"
    local draw_mode=1
    [[ "$mode" == "wire" ]] && draw_mode=0
    local hierarchy_events=""
    if [[ -n "$hierarchy_root" && -n "$hierarchy_child" &&
	    -n "$hierarchy_path" ]]; then
	hierarchy_events=$(cat <<EOF
,
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "set_expanded",
     "arguments": {"labels": ["${hierarchy_root}"], "expanded": true}},
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/tree-expanded.png"}},
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "set_current",
     "arguments": {"labels": ["${hierarchy_root}", "${hierarchy_child}"]}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/tree-selected.png"}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/selection-visible.png"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "erase ${hierarchy_path}"}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/subpath-erased.png"}},
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/tree-erased.png"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "draw -m${draw_mode} ${hierarchy_path}"}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/subpath-redraw-return.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/subpath-redraw-stable.png"}},
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/tree-redrawn.png"}},
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "clear_selection",
     "arguments": {}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/tree-cleared.png"}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/final-stable.png"}}
EOF
)
    fi

    cat > "$output" <<EOF
{
  "schema": "brlcad.qtcad.events",
  "version": 1,
  "events": [
    {"target": ".", "action": "resize",
     "arguments": {"width": 1100, "height": 800}},
    {"target": ".", "action": "wait", "arguments": {"ms": 100}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "draw -m${draw_mode} ${object}"}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/draw-return.png"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/draw-050ms.png"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "ae 90 0"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "autoview"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 200}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/ae90-0200ms.png"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 1300}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/ae90-1500ms.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/ae90-stable.png"}},

    {"target": "${canvas_target}", "action": "wheel",
     "arguments": {"x": 0.5, "y": 0.5, "pixel_x": 0, "pixel_y": 0,
                   "angle_x": 0, "angle_y": 360, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/zoom-in-motion.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/zoom-in-stable.png"}},

    {"target": "${canvas_target}", "action": "wheel",
     "arguments": {"x": 0.5, "y": 0.5, "pixel_x": 0, "pixel_y": 0,
                   "angle_x": 0, "angle_y": -720, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/zoom-out-motion.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/zoom-out-stable.png"}},

    {"target": "${canvas_target}", "action": "wheel",
     "arguments": {"x": 0.5, "y": 0.5, "pixel_x": 0, "pixel_y": 0,
                   "angle_x": 0, "angle_y": 360, "modifiers": 0}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/zoom-return-stable.png"}},

    {"target": "${canvas_target}", "action": "mouse_press",
     "arguments": {"x": 0.34, "y": 0.45, "button": 1, "buttons": 1,
                   "modifiers": 0}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.40, "y": 0.47, "button": 0, "buttons": 1,
                   "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.46, "y": 0.51, "button": 0, "buttons": 1,
                   "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.52, "y": 0.55, "button": 0, "buttons": 1,
                   "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.58, "y": 0.51, "button": 0, "buttons": 1,
                   "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/rotate-held-end.png"}},
    {"target": "${canvas_target}", "action": "mouse_release",
     "arguments": {"x": 0.58, "y": 0.51, "button": 1, "buttons": 0,
                   "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/rotate-motion.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/rotate-stable.png"}}
${hierarchy_events}
  ]
}
EOF
}

validate_report()
{
    local report="$1"
    local image_dir="$2"
    local object="$3"
    local hierarchy_path="$4"
    local validation="$5"

    if ! jq -e '
	.success == true and
	(.samples | length) > 0 and
	(any(.samples[]; (.draw_shape_count // 0) > 0)) and
	(all(.samples[]; (.failed_sources // 0) == 0)) and
	(all(.samples[]; (.cad_payloads_without_entry // 0) == 0)) and
	(all(.samples[];
	    (.superseded_fallback_presentations // 0) == 0))
	and
	(.samples[-1].progressive_pending == false) and
	(.samples[-1].lod_results_pending == false) and
	(.samples[-1].lod_submissions_pending == false) and
	(.samples[-1].lod_refinement_frame_pending == false) and
	((.samples[-1].lod_service_pending_tasks // 0) == 0) and
	((.samples[-1].lod_service_active_requests // 0) == 0) and
	((.samples[-1].lod_service_queued_results // 0) == 0) and
	((.samples[-1].lod_service_queued_cache_writes // 0) == 0) and
	((.samples[-1].active_lod_aabb_payloads // 0) == 0) and
	((.samples[-1].active_lod_obb_payloads // 0) == 0) and
	((.samples[-1].active_lod_sphere_payloads // 0) == 0) and
	((.samples[-1].compact_lod_entries // 0) ==
	 (.samples[-1].compact_lod_entries_with_payload // 0))
	' "$report" >"$validation" 2>&1; then
	return 1
    fi

    # A retained database-source node is counted before it owns drawable
    # geometry, so report counters alone cannot distinguish a useful first
    # frame from the empty gradient canvas.  Check the actual 1.5-second
    # framebuffer: the background is horizontally uniform, while boxes and
    # geometry differ from the left-edge color on enough pixels to be useful.
    # Exclude the HUD strip at the bottom.
    local first_useful="$image_dir/ae90-1500ms.png"
    local first_useful_elapsed
    first_useful_elapsed=$(jq -r '
	first(.samples[] |
	    select((.checkpoint? // "") | endswith("/ae90-1500ms.png")) |
	    (.elapsed_ms // 9223372036854775807))
	' "$report" 2>/dev/null)
    local dimensions width height crop_height background changed_pixels
    dimensions=$(identify -format '%wx%h' "$first_useful" 2>/dev/null || true)
    width="${dimensions%x*}"
    height="${dimensions#*x}"
    crop_height=0
    if [[ "$width" =~ ^[0-9]+$ && "$height" =~ ^[0-9]+$ &&
	    "$height" -gt 90 ]]; then
	crop_height=$((height - 90))
    fi
    background=$(mktemp "$artifact_dir/.qged-background.XXXXXX.png")
    changed_pixels=0
    if [[ "$crop_height" -gt 0 ]] &&
	convert "$first_useful" -crop "1x${crop_height}+0+0" +repage \
	    -scale "${width}x${crop_height}!" "$background" 2>/dev/null; then
	changed_pixels=$(compare -metric AE -fuzz 3% \
	    -crop "${width}x${crop_height}+0+0" \
	    "$first_useful" "$background" null: 2>&1 || true)
    fi
    rm -f "$background"
    if [[ ! "$first_useful_elapsed" =~ ^[0-9]+$ ||
	    "$first_useful_elapsed" -gt 5000 ||
	    ! "$changed_pixels" =~ ^[0-9]+$ ||
	    "$changed_pixels" -lt 1000 ]]; then
	printf 'no useful model pixels by first-frame checkpoint: elapsed=%s changed_pixels=%s\n' \
	    "$first_useful_elapsed" "$changed_pixels" >>"$validation"
	return 1
    fi

    if [[ -n "$hierarchy_path" ]]; then
	if ! jq -e --arg object "$object" --arg path "$hierarchy_path" '
	    (first(.samples[] |
		select(.action == "set_current") | .event_index)) as $selected |
	    (first(.samples[] |
		select(.command? == ("erase " + $path)) |
		.event_index)) as $erased |
	    (any(.samples[];
		.selection_paths? != null and
		(.selection_paths | index($path)) != null)) and
	    (all(.samples[];
		if (.event_index > $selected and .event_index < $erased)
		then (.progressive_pending == false and
		      .lod_submissions_pending == false)
		else true end)) and
	    (any(.samples[];
		.command? == ("erase " + $path) and
		(.draw_shape_count // 0) > 0 and
		(.draw_frontier_paths | index($object)) != null)) and
	    (any(.samples[];
		(.command? | type) == "string" and
		(.command | endswith(" " + $path)) and
		(.command | startswith("draw ")) and
		(.draw_shape_count // 0) > 0 and
		(.draw_frontier_paths | index($object)) != null)) and
	    (all(.samples[];
		if (.action == "set_expanded" or
		    .action == "set_current" or
		    .action == "clear_selection" or
		    ((.command? | type) == "string" and
		     ((.command == ("erase " + $path)) or
		      ((.command | startswith("draw ")) and
		       (.command | endswith(" " + $path))))))
		then ((.event_duration_us // 9223372036854775807) <= 250000)
		else true end)) and
	    ((.samples[-1].tree_notify_path_us // 9223372036854775807) <=
	     100000) and
	    ((.samples[-1].tree_notify_full_items // 0) <= 16) and
	    ((.samples[-1].selection_count // -1) == 0)
	    ' "$report" >>"$validation" 2>&1; then
	    return 1
	fi

	# Counters can agree while a stale aggregate, lost selection style, or
	# missed redraw leaves the framebuffer unchanged.  Require the real tree
	# and canvas pixels to show all four user-facing transitions.  A low
	# threshold is deliberate: Hubble contains valid selectable components
	# only a few pixels wide at the matrix's initial view.
	local tree_selection_pixels erase_pixels redraw_pixels clear_pixels
	tree_selection_pixels=$(compare -metric AE -fuzz 3% \
	    "$image_dir/tree-expanded.png" "$image_dir/tree-selected.png" \
	    null: 2>&1 || true)
	erase_pixels=$(compare -metric AE -fuzz 3% \
	    "$image_dir/selection-visible.png" "$image_dir/subpath-erased.png" \
	    null: 2>&1 || true)
	redraw_pixels=$(compare -metric AE -fuzz 3% \
	    "$image_dir/subpath-erased.png" \
	    "$image_dir/subpath-redraw-stable.png" null: 2>&1 || true)
	clear_pixels=$(compare -metric AE -fuzz 3% \
	    "$image_dir/subpath-redraw-stable.png" \
	    "$image_dir/final-stable.png" null: 2>&1 || true)
	if [[ ! "$tree_selection_pixels" =~ ^[0-9]+$ ||
		"$tree_selection_pixels" -lt 100 ||
		! "$erase_pixels" =~ ^[0-9]+$ || "$erase_pixels" -lt 12 ||
		! "$redraw_pixels" =~ ^[0-9]+$ || "$redraw_pixels" -lt 12 ||
		! "$clear_pixels" =~ ^[0-9]+$ || "$clear_pixels" -lt 12 ]]; then
	    printf 'missing hierarchy visual transition: tree=%s erase=%s redraw=%s clear=%s\n' \
		"$tree_selection_pixels" "$erase_pixels" "$redraw_pixels" \
		"$clear_pixels" >>"$validation"
	    return 1
	fi
    fi

    # A large software-rendered mesh must actually shed work while the mouse
    # is held.  Merely raising a policy number is not enough: verify the
    # submitted face cut and the resulting frame time.  Small/non-progressive
    # scenes are intentionally outside this stress assertion.
    if ! jq -e '
	if .backend != "osmesa" then true
	else
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-held-end.png")))) as $motion |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-stable.png")))) as $stable |
	    if (($stable.active_progressive_cad_faces // 0) < 100000)
	    then true
	    else
		(($motion.lod_target_pixel_error // 1) > 1) and
		(($motion.active_progressive_cad_faces // 0) <=
		 ($stable.active_progressive_cad_faces // 0)) and
		(($motion.last_render_ms // 9223372036854775807) <= 250)
	    end
	end
	' "$report" >>"$validation" 2>&1; then
	printf 'software interaction did not select and render a bounded coarse cut\n' \
	    >>"$validation"
	return 1
    fi

    if find "$image_dir" -type f -size 0 -print -quit | grep -q .; then
	printf 'zero-length checkpoint image\n' >>"$validation"
	return 1
    fi
    return 0
}

run_current()
{
    local case_name="$1"
    local db="$2"
    local object="$3"
    local backend="$4"
    local mode="$5"
    local swap="$6"
    local cache_state="$7"
    local cache_dir="$8"
    local settle_ms="$9"
    local hierarchy_root="${10}"
    local hierarchy_child="${11}"
    local hierarchy_path="${12}"
    local swap_tag="${swap//-/_}"
    local run_name="${case_name}-${backend}-${mode}-swap${swap_tag}-${cache_state}"
    local out="$artifact_dir/cases/$run_name"
    local events="$artifact_dir/events/$run_name.json"
    mkdir -p "$out/images"
    write_event_script "$events" "$out/images" "$mode" "$object" \
	"$settle_ms" "$hierarchy_root" "$hierarchy_child" "$hierarchy_path"

    local env_args=("BU_DIR_CACHE=$cache_dir")
    if [[ "$swap" != "default" ]]; then
	env_args+=("QGED_SWAP_INTERVAL=$swap")
    fi
    local qged_args=()
    [[ "$backend" == "osmesa" ]] && qged_args+=("-s")
    qged_args+=("--test-script" "$events" "--test-report" "$out/report.json" "$db")

    local command=(env "${env_args[@]}" "$qged" "${qged_args[@]}")
    if [[ "$case_name" == "$perf_case" && "$backend" == "system" &&
	    "$cache_state" == "cold" && "$mode" == "shaded" &&
	    "$swap" == "default" ]]; then
	command=(perf record -g -o "$out/perf.data" -- "${command[@]}")
    fi
    if [[ "$case_name" == "$apitrace_case" && "$backend" == "system" &&
	    "$cache_state" == "cold" && "$mode" == "shaded" &&
	    "$swap" == "default" ]]; then
	command=(apitrace trace --api=gl -o "$out/qged.trace" "${command[@]}")
    fi

    printf 'RUN %s\n' "$run_name"
    local started=$SECONDS
    if timeout --signal=TERM "$run_timeout" "${command[@]}" \
	    >"$out/stdout.log" 2>"$out/stderr.log"; then
	if validate_report "$out/report.json" "$out/images" "$object" \
		"$hierarchy_path" "$out/validation.log"; then
	    printf 'PASS,%s,%s,\n' "$run_name" "$((SECONDS - started))" \
		>> "$artifact_dir/results.csv"
	    return 0
	fi
	printf 'FAIL,%s,%s,report-validation\n' "$run_name" \
	    "$((SECONDS - started))" >> "$artifact_dir/results.csv"
	return 1
    fi
    local status=$?
    printf 'FAIL,%s,%s,status=%s\n' "$run_name" "$((SECONDS - started))" \
	"$status" >> "$artifact_dir/results.csv"
    return 1
}

capture_baseline()
{
    local case_name="$1"
    local db="$2"
    local object="$3"
    local mode="$4"
    local out="$artifact_dir/baseline/${case_name}-${mode}"
    local draw_mode=1
    [[ "$mode" == "wire" ]] && draw_mode=0
    mkdir -p "$out"

    if [[ ! -x "$baseline_qged" ]]; then
	echo "SKIP: baseline qged not found at $baseline_qged" > "$out/status.txt"
	return 0
    fi
    if [[ -z "${DISPLAY:-}" ]] || ! command -v xdotool >/dev/null ||
	    ! command -v import >/dev/null ||
	    ! command -v convert >/dev/null ||
	    ! command -v compare >/dev/null; then
	echo "SKIP: baseline capture needs DISPLAY, xdotool, import, convert, and compare" \
	    > "$out/status.txt"
	return 0
    fi

    "$baseline_qged" "$db" >"$out/stdout.log" 2>"$out/stderr.log" &
    local pid=$!
    local window=""
    for _ in $(seq 1 100); do
	window="$(xdotool search --onlyvisible --pid "$pid" 2>/dev/null |
	    tail -n 1)"
	[[ -n "$window" ]] && break
	sleep 0.1
    done
    if [[ -z "$window" ]]; then
	echo "FAIL: no baseline window" > "$out/status.txt"
	kill "$pid" 2>/dev/null || true
	wait "$pid" 2>/dev/null || true
	return 1
    fi

    eval "$(xdotool getwindowgeometry --shell "$window")"
    xdotool windowraise "$window"
    xdotool windowactivate --sync "$window"
    import -window "$window" "$out/initial.png"
    # The console is a bottom dock.  Its prompt is on the first text line,
    # not in the otherwise blank middle of the editor.  The historical main
    # qged layout restores a roughly 350-pixel console, placing that line
    # about 320 pixels above the bottom of the client.
    local console_y=$((HEIGHT - 320))
    if [[ "$console_y" -lt $((HEIGHT / 2)) ]]; then
	console_y=$((HEIGHT - 100))
    fi
    xdotool mousemove --window "$window" "$((WIDTH / 2))" "$console_y" click 1
    xdotool key End
    xdotool type --delay 1 "draw -m${draw_mode} ${object}"
    xdotool key Return
    sleep 3
    xdotool type --delay 1 "ae 90 0"
    xdotool key Return
    xdotool type --delay 1 "autoview"
    xdotool key Return
    sleep 5
    import -window "$window" "$out/ae90-stable.png"
    xdotool mousemove --window "$window" "$((WIDTH / 2))" "$((HEIGHT / 2))"
    xdotool click 4 click 4 click 4
    sleep 2
    import -window "$window" "$out/zoom-in-stable.png"

    # Window capture succeeding does not prove the commands reached qged.
    # Compare a conservative central viewport crop: drawing must change the
    # initially empty canvas, and zooming must then change the drawn canvas.
    # This caught a former false-positive path that archived three identical
    # empty screenshots and labelled the baseline PASS.
    local crop_x=$((WIDTH / 4))
    local crop_y=$((HEIGHT / 5))
    local crop_width=$((WIDTH / 2))
    local crop_height=$((HEIGHT / 2))
    local crop_geometry="${crop_width}x${crop_height}+${crop_x}+${crop_y}"
    convert "$out/initial.png" -crop "$crop_geometry" +repage \
	"$out/initial-view.png"
    convert "$out/ae90-stable.png" -crop "$crop_geometry" +repage \
	"$out/ae90-stable-view.png"
    convert "$out/zoom-in-stable.png" -crop "$crop_geometry" +repage \
	"$out/zoom-in-stable-view.png"
    if compare -metric AE "$out/initial-view.png" \
	    "$out/ae90-stable-view.png" null: >/dev/null 2>&1; then
	echo "FAIL: baseline draw did not change the viewport" > "$out/status.txt"
	kill "$pid" 2>/dev/null || true
	wait "$pid" 2>/dev/null || true
	return 1
    fi
    if compare -metric AE "$out/ae90-stable-view.png" \
	    "$out/zoom-in-stable-view.png" null: >/dev/null 2>&1; then
	echo "FAIL: baseline zoom did not change the viewport" > "$out/status.txt"
	kill "$pid" 2>/dev/null || true
	wait "$pid" 2>/dev/null || true
	return 1
    fi
    kill "$pid" 2>/dev/null || true
    wait "$pid" 2>/dev/null || true
    echo "PASS" > "$out/status.txt"
}

printf 'status,run,seconds,detail\n' > "$artifact_dir/results.csv"
{
    printf 'source_root=%s\nbuild_dir=%s\nmain_build_dir=%s\n' \
	"$source_root" "$build_dir" "$main_build_dir"
    printf 'profile=%s\ncases=%s\nbackends=%s\nmodes=%s\nswap_intervals=%s\n' \
	"$profile" "$case_list" "$backend_list" "$mode_list" "$swap_list"
    printf 'display=%s\nsession_type=%s\n' "${DISPLAY:-}" "${XDG_SESSION_TYPE:-}"
    printf 'qged_sha256='
    sha256sum "$qged" | cut -d' ' -f1
    if [[ -x "$baseline_qged" ]]; then
	printf 'baseline_qged_sha256='
	sha256sum "$baseline_qged" | cut -d' ' -f1
    fi
} > "$artifact_dir/run-info.txt"

failures=0
IFS=',' read -r -a cases <<< "$case_list"
IFS=',' read -r -a backends <<< "$backend_list"
IFS=',' read -r -a modes <<< "$mode_list"
IFS=',' read -r -a swaps <<< "$swap_list"

for case_name in "${cases[@]}"; do
    spec="$(case_spec "$case_name")" || {
	echo "ERROR: unknown case '$case_name'" >&2
	failures=$((failures + 1))
	continue
    }
    IFS='|' read -r db object hierarchy_root hierarchy_child hierarchy_path \
	<<< "$spec"
    if [[ ! -r "$db" ]]; then
	echo "ERROR: missing database for $case_name: $db" >&2
	failures=$((failures + 1))
	continue
    fi
    write_inventory "$case_name" "$db" "$object" "$hierarchy_root" \
	"$hierarchy_child" "$hierarchy_path"

    if [[ "$run_baseline" -eq 1 ]]; then
	for mode in "${modes[@]}"; do
	    capture_baseline "$case_name" "$db" "$object" "$mode" ||
		failures=$((failures + 1))
	done
    fi

    case "$profile" in
	smoke) settle_ms=15000 ;;
	full) settle_ms=30000 ;;
	stress) settle_ms=60000 ;;
    esac
    for backend in "${backends[@]}"; do
	if [[ "$backend" != "system" && "$backend" != "osmesa" ]]; then
	    echo "ERROR: unknown backend '$backend'" >&2
	    failures=$((failures + 1))
	    continue
	fi
	for mode in "${modes[@]}"; do
	    if [[ "$mode" != "shaded" && "$mode" != "wire" ]]; then
		echo "ERROR: unknown mode '$mode'" >&2
		failures=$((failures + 1))
		continue
	    fi
	    for swap in "${swaps[@]}"; do
		if [[ "$swap" != "default" && "$swap" != "0" &&
			"$swap" != "1" && "$swap" != "-1" ]]; then
		    echo "ERROR: invalid swap interval '$swap'" >&2
		    failures=$((failures + 1))
		    continue
		fi
		pair="$artifact_dir/caches/${case_name}-${backend}-${mode}-swap${swap//-/_}"
		mkdir -p "$pair"
		cache="$pair/cache"
		# BU_DIR_CACHE itself must exist or libbu correctly disables cache
		# writes.  This is still a completely cold cache: the new directory
		# contains neither format metadata nor cached payloads.
		mkdir "$cache"
		run_current "$case_name" "$db" "$object" "$backend" "$mode" \
		    "$swap" "cold" "$cache" "$settle_ms" "$hierarchy_root" \
		    "$hierarchy_child" "$hierarchy_path" ||
		    failures=$((failures + 1))
		run_current "$case_name" "$db" "$object" "$backend" "$mode" \
		    "$swap" "warm" "$cache" "$settle_ms" "$hierarchy_root" \
		    "$hierarchy_child" "$hierarchy_path" ||
		    failures=$((failures + 1))
		find "$cache" -type f -printf '%s %T@ %p\n' 2>/dev/null |
		    sort -n > "$pair/cache-files.txt"
	    done
	done
    done

    after_size="$(stat -c %s "$db")"
    before_size="$(sed -n 's/^size=//p' "$artifact_dir/inventory/$case_name.txt" |
	head -n 1)"
    if [[ "$after_size" != "$before_size" ]]; then
	echo "ERROR: input database size changed during $case_name" >&2
	failures=$((failures + 1))
    fi
done

{
    printf '# qged GUI matrix summary\n\n'
    printf -- '- Artifact directory: `%s`\n' "$artifact_dir"
    printf -- '- Failures: %s\n\n' "$failures"
    printf '| Status | Run | Seconds | Detail |\n'
    printf '|---|---|---:|---|\n'
    tail -n +2 "$artifact_dir/results.csv" |
	while IFS=',' read -r status run seconds detail; do
	    printf '| %s | %s | %s | %s |\n' "$status" "$run" "$seconds" "$detail"
	done
} > "$artifact_dir/SUMMARY.md"

echo "Artifacts: $artifact_dir"
echo "Failures: $failures"
[[ "$failures" -eq 0 ]]
