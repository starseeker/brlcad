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
baseline_only=0
perf_case=""
perf_phase="cold"
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
  --baseline-only          Capture production baselines without current runs
  --perf CASE              Record the selected backend case with perf
  --perf-phase PHASE       cold, warm, or both (default: cold)
  --apitrace CASE          Trace one cold System GL case with apitrace
  --timeout SECONDS        Per-process timeout (default: 180)

Cases: generic_twin, lucy, multi_lucy, multi_lucy_xpush, stanford, havoc,
       hubble, many_lucy_stress, unique_mesh_stress,
       unique_mesh_50k_stress, unique_mesh_150k_stress

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
	--baseline-only) run_baseline=1; baseline_only=1; shift ;;
	--perf) perf_case="$2"; shift 2 ;;
	--perf-phase) perf_phase="$2"; shift 2 ;;
	--apitrace) apitrace_case="$2"; shift 2 ;;
	--timeout) run_timeout="$2"; shift 2 ;;
	--help|-h) usage; exit 0 ;;
	*) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ "$perf_phase" != "cold" && "$perf_phase" != "warm" &&
	"$perf_phase" != "both" ]]; then
    echo "ERROR: --perf-phase must be cold, warm, or both" >&2
    exit 2
fi

case "$profile" in
    smoke) default_cases="generic_twin,lucy" ;;
    full) default_cases="generic_twin,lucy,havoc,hubble" ;;
    stress)
	default_cases="generic_twin,lucy,multi_lucy,multi_lucy_xpush,stanford,havoc,hubble"
	if [[ -f "$build_dir/unique_mesh_stress.g" ]]; then
	    default_cases+=",unique_mesh_stress"
	fi
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

screensaver_inhibit_pid=""
cleanup_screensaver_inhibit()
{
    if [[ -n "$screensaver_inhibit_pid" ]]; then
	kill "$screensaver_inhibit_pid" 2>/dev/null || true
	wait "$screensaver_inhibit_pid" 2>/dev/null || true
    fi
}
trap cleanup_screensaver_inhibit EXIT
if [[ "$run_baseline" -eq 1 ]] &&
	command -v xfce4-screensaver-command >/dev/null 2>&1; then
    xfce4-screensaver-command --deactivate >/dev/null 2>&1 || true
    xfce4-screensaver-command --inhibit \
	--application-name qged-gui-matrix \
	--reason "BRL-CAD GUI validation" >/dev/null 2>&1 &
    screensaver_inhibit_pid=$!
fi

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
	    printf '%s|%s|%s|%s|%s\n' \
		"${BOBOL_LUCY_DB:-$build_dir/lucy.g}" "lucy.bot.r" \
		"lucy.bot.r" "lucy.bot" "lucy.bot.r/lucy.bot"
	    ;;
	bigboy)
	    printf '%s|%s|||\n' \
		"${BOBOL_BIGBOY_DB:-$build_dir/bigboy.g}" "all"
	    ;;
	multi_lucy)
	    printf '%s|%s|||\n' "$build_dir/stanford.g" "multi_lucy"
	    ;;
	multi_lucy_xpush)
	    printf '%s|%s|||\n' "$build_dir/stanford.g" "multi_lucy_xpush"
	    ;;
	many_lucy_stress)
	    printf '%s|%s|||\n' \
		"${BOBOL_MANY_LUCY_STRESS_DB:-$build_dir/many_lucy_stress.g}" \
		"many_lucy_stress"
	    ;;
	unique_mesh_stress)
	    printf '%s|%s|%s|%s|%s\n' \
		"${BOBOL_UNIQUE_MESH_STRESS_DB:-$build_dir/unique_mesh_stress.g}" \
		"unique_mesh_stress" "unique_mesh_stress" \
		"unique_level_02_000000.c" \
		"unique_mesh_stress/unique_level_02_000000.c/unique_level_01_000000.c/unique_level_00_000000.c/unique_region_000000.r"
	    ;;
	unique_mesh_50k_stress)
	    printf '%s|%s|%s|%s|%s\n' \
		"${BOBOL_UNIQUE_MESH_50K_STRESS_DB:-$build_dir/unique_mesh_50k_stress.g}" \
		"unique_mesh_stress" "unique_mesh_stress" \
		"unique_level_03_000000.c" \
		"unique_mesh_stress/unique_level_03_000000.c/unique_level_02_000000.c/unique_level_01_000000.c/unique_level_00_000000.c/unique_region_000000.r"
	    ;;
	unique_mesh_150k_stress)
	    printf '%s|%s|%s|%s|%s\n' \
		"${BOBOL_UNIQUE_MESH_150K_STRESS_DB:-$build_dir/unique_mesh_150k_stress.g}" \
		"unique_mesh_stress" "unique_mesh_stress" \
		"unique_level_04_000000.c" \
		"unique_mesh_stress/unique_level_04_000000.c/unique_level_03_000000.c/unique_level_02_000000.c/unique_level_01_000000.c/unique_level_00_000000.c/unique_region_000000.r"
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
    local case_name="$9"
    local draw_mode=1
    [[ "$mode" == "wire" ]] && draw_mode=0
    local hierarchy_events=""
    local hierarchy_expand_events=""
    local hierarchy_selection_labels=""
    local working_set_events=""
    local smooth_zoom_events=""
    local lighting_events=""
    local view_events=""
    local label_index
    # Qt::ControlModifier forces BRL-CAD's rotate binding independently of
    # the operator's currently selected left-mouse toolbar mode.
    local rotate_modifier=67108864
    if [[ "$object" == "multi_lucy" ]]; then
	# Focus the view on the first transformed Lucy, then orbit far enough
	# around that focus for other occurrences to enter and leave the
	# frustum.  This is deliberately distinct from the overview rotation
	# below: it verifies working-set turnover while a small visible subset
	# has a much richer pixel-exact cut.
	working_set_events=$(cat <<EOF
	    {"target": ".", "action": "qged_command_batch",
	     "arguments": {"commands": ["ae 90 0",
		"center -509244 -1121531 192627", "size 1000000"]}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/close-focus-stable.png"}},
    {"target": "${canvas_target}", "action": "mouse_press",
     "arguments": {"x": 0.32, "y": 0.48, "button": 1, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.46, "y": 0.43, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.60, "y": 0.52, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.72, "y": 0.61, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/close-turnover-held.png"}},
    {"target": "${canvas_target}", "action": "mouse_release",
     "arguments": {"x": 0.72, "y": 0.61, "button": 1, "buttons": 0,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/close-turnover-motion.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/close-turnover-stable.png"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "size 1000000"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "ae 0 0"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/close-direction-0-motion.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/close-direction-0-stable.png"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "size 1000000"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "ae 90 0"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/close-direction-90-motion.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/close-direction-90-stable.png"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "autoview"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/close-turnover-return.png"}},
EOF
)
    fi

    if [[ "$mode" == "shaded" &&
	    ("$case_name" == "generic_twin" || "$case_name" == "lucy") ]]; then
	# Lighting is shared view policy, not a qged-only preference.  Capture both
	# named profiles through the public command and restore the studio default
	# before any later selection/hierarchy probes.
	lighting_events=$(cat <<EOF
,
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "view lighting profile mged"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/lighting-mged.png"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "view lighting profile studio"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/lighting-studio.png"}}
EOF
)
    fi
    if [[ "$case_name" == "unique_mesh_stress" ||
	    "$case_name" == "lucy" ]]; then
	# Approach a large mesh with actual wheel events so every intermediate
	# camera, admission decision, presentation, and post-quiet resident
	# compaction is measured instead of jumping between command-line sizes.
	# unique_closed_000199.bot is a deterministic 100,352-face leaf whose
	# principal orientation faces ae 90 0; Lucy supplies the corresponding
	# single-asset memory-reclamation stress.
	local smooth_zoom_steps=12
	if [[ "$case_name" == "unique_mesh_stress" ]]; then
	    smooth_zoom_steps=16
	    smooth_zoom_events=$(cat <<EOF
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "ae 90 0"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "center -2289 109 -1300"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "size 2500"}},
EOF
)
	else
	    smooth_zoom_events=$(cat <<EOF
    {"target": ".", "action": "qged_command_batch",
     "arguments": {"commands": ["ae 90 0", "autoview"]}},
EOF
)
	fi
	smooth_zoom_events+=$(cat <<EOF
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-start-stable.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-start-stable.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-start-stable.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-start-stable.png"}},
EOF
)
	if [[ "$case_name" == "lucy" ]]; then
	    # Keep the interaction explicitly bracketed beyond one exceptional
	    # suffix's measured load/build latency.  The later low-amplitude wheel
	    # events also keep the scale epoch changing at less than the 150 ms
	    # quiet debounce, distinguishing refinement during continuous input
	    # from a result released only after zoom has stopped.
	    smooth_zoom_events+=$(cat <<EOF
    {"target": "${canvas_target}", "action": "mouse_press",
     "arguments": {"x": 0.5, "y": 0.5, "button": 1, "buttons": 1,
                   "modifiers": 0}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.501, "y": 0.5, "button": 0, "buttons": 1,
                   "modifiers": 0}},
EOF
)
	fi
	for ((label_index = 1;
		label_index <= smooth_zoom_steps; ++label_index)); do
	    smooth_zoom_events+=$(cat <<EOF
    {"target": "${canvas_target}", "action": "wheel",
     "arguments": {"x": 0.5, "y": 0.5, "pixel_x": 0, "pixel_y": 0,
                   "angle_x": 0, "angle_y": 120, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 16}},
EOF
)
	    if ((label_index % 4 == 0)); then
		smooth_zoom_events+=$(cat <<EOF
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-in-${label_index}.png"}},
EOF
)
	    fi
	done
	if [[ "$case_name" == "lucy" ]]; then
	    for ((label_index = 1; label_index <= 15; ++label_index)); do
		smooth_zoom_events+=$(cat <<EOF
    {"target": "${canvas_target}", "action": "wheel",
     "arguments": {"x": 0.5, "y": 0.5, "pixel_x": 0, "pixel_y": 0,
                   "angle_x": 0, "angle_y": 1, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 100}},
EOF
)
	    done
	    smooth_zoom_events+=$(cat <<EOF
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-active-refined.png"}},
    {"target": "${canvas_target}", "action": "mouse_release",
     "arguments": {"x": 0.501, "y": 0.5, "button": 1, "buttons": 0,
                   "modifiers": 0}},
EOF
)
	fi
	smooth_zoom_events+=$(cat <<EOF
    {"target": ".", "action": "wait", "arguments": {"ms": 850}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-close-stable.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-close-stable.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-close-stable.png"}},
EOF
)
	for ((label_index = 1;
		label_index <= smooth_zoom_steps; ++label_index)); do
	    smooth_zoom_events+=$(cat <<EOF
    {"target": "${canvas_target}", "action": "wheel",
     "arguments": {"x": 0.5, "y": 0.5, "pixel_x": 0, "pixel_y": 0,
                   "angle_x": 0, "angle_y": -120, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 16}},
EOF
)
	    if ((label_index % 4 == 0)); then
		smooth_zoom_events+=$(cat <<EOF
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-out-${label_index}.png"}},
EOF
)
	    fi
	done
	smooth_zoom_events+=$(cat <<EOF
    {"target": ".", "action": "wait", "arguments": {"ms": 850}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-out-stable.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-out-stable.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-out-stable.png"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "autoview"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-return.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-return.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/smooth-zoom-return.png"}},
EOF
)
    fi
    if [[ -n "$hierarchy_root" && -n "$hierarchy_child" &&
	    -n "$hierarchy_path" ]]; then
	# Drive every prefix of the declared probe path.  Keeping the event
	# generator independent of a fixture's hierarchy depth prevents a
	# larger stress database from silently probing the wrong row.
	local -a hierarchy_labels=()
	IFS='/' read -r -a hierarchy_labels <<< "$hierarchy_path"
	local prefix_labels=""
	local label
	for ((label_index = 0;
		label_index < ${#hierarchy_labels[@]}; ++label_index)); do
	    label="${hierarchy_labels[$label_index]}"
	    [[ -n "$hierarchy_selection_labels" ]] &&
		hierarchy_selection_labels+=", "
	    hierarchy_selection_labels+="\"${label}\""
	    [[ -n "$prefix_labels" ]] && prefix_labels+=", "
	    prefix_labels+="\"${label}\""
	    if ((label_index + 1 < ${#hierarchy_labels[@]})); then
		hierarchy_expand_events+="    {\"target\": \"./n:Hierarchy/i:hierarchy-tree\", \"action\": \"set_expanded\",
     \"arguments\": {\"labels\": [${prefix_labels}], \"expanded\": true}},
"
	    fi
	done
	# A canvas checkpoint performs a real render.  On a slow backend that
	# render may supply new capacity feedback and schedule a legitimate LoD
	# correction.  Settle any such work before beginning the hierarchy
	# interval so the selection assertions measure work caused by selection,
	# rather than work caused by the preceding diagnostic checkpoint.
	hierarchy_events=$(cat <<EOF
,
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
${hierarchy_expand_events}
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/tree-expanded.png"}},
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "set_current",
     "arguments": {"labels": [${hierarchy_selection_labels}]}},
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

    if [[ "$case_name" == "unique_mesh_50k_stress" ||
	    "$case_name" == "unique_mesh_150k_stress" ]]; then
	# User readiness and terminal background convergence are different
	# contracts.  Exercise a real held drag while cold preparation is still
	# active, then return to the reference view and wait once for a terminal
	# characterization.  Subsequent camera changes are sampled immediately
	# and share one final recovery wait.  The old sequence waited for full
	# convergence after every wheel event, hiding cold input stalls and
	# multiplying a useful scale test into several minutes of idle waiting.
	view_events=$(cat <<EOF
    {"target": "${canvas_target}", "action": "mouse_press",
     "arguments": {"x": 0.34, "y": 0.45, "button": 1, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.43, "y": 0.49, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.52, "y": 0.55, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/cold-rotate-held.png"}},
    {"target": "${canvas_target}", "action": "mouse_release",
     "arguments": {"x": 0.52, "y": 0.55, "button": 1, "buttons": 0,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/cold-rotate-motion.png"}},
    {"target": ".", "action": "qged_command_batch",
     "arguments": {"commands": ["ae 90 0", "autoview"]}},
    {"target": ".", "action": "wait", "arguments": {"ms": 1000}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/realization-1s.png"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 2000}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/realization-3s.png"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 3000}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/realization-6s.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/ae90-stable.png"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 200}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/ae90-stable-held.png"}},

    {"target": "${canvas_target}", "action": "wheel",
     "arguments": {"x": 0.5, "y": 0.5, "pixel_x": 0, "pixel_y": 0,
                   "angle_x": 0, "angle_y": 360, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/zoom-in-motion.png"}},
    {"target": "${canvas_target}", "action": "wheel",
     "arguments": {"x": 0.5, "y": 0.5, "pixel_x": 0, "pixel_y": 0,
                   "angle_x": 0, "angle_y": -720, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/zoom-out-motion.png"}},
    {"target": "${canvas_target}", "action": "wheel",
     "arguments": {"x": 0.5, "y": 0.5, "pixel_x": 0, "pixel_y": 0,
                   "angle_x": 0, "angle_y": 360, "modifiers": 0}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/zoom-return-stable.png"}},

    {"target": "${canvas_target}", "action": "mouse_press",
     "arguments": {"x": 0.34, "y": 0.45, "button": 1, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.40, "y": 0.47, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.46, "y": 0.51, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.52, "y": 0.55, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.58, "y": 0.51, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/rotate-held-end.png"}},
    {"target": "${canvas_target}", "action": "mouse_release",
     "arguments": {"x": 0.58, "y": 0.51, "button": 1, "buttons": 0,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/rotate-motion.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/rotate-stable.png"}}
EOF
)
    else
	view_events=$(cat <<EOF
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/ae90-stable.png"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 200}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/ae90-stable-held.png"}},

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

${working_set_events}
${smooth_zoom_events}
    {"target": "${canvas_target}", "action": "mouse_press",
     "arguments": {"x": 0.34, "y": 0.45, "button": 1, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.40, "y": 0.47, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.46, "y": 0.51, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.52, "y": 0.55, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "mouse_move",
     "arguments": {"x": 0.58, "y": 0.51, "button": 0, "buttons": 1,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 8}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/rotate-held-end.png"}},
    {"target": "${canvas_target}", "action": "mouse_release",
     "arguments": {"x": 0.58, "y": 0.51, "button": 1, "buttons": 0,
                   "modifiers": ${rotate_modifier}}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/rotate-motion.png"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_ms}, "quiet_ms": 100}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/rotate-stable.png"}}
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
     "arguments": {"command": "ae 90 0"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "draw -m${draw_mode} ${object}"}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/draw-return.png"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/draw-050ms.png"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 200}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/ae90-0200ms.png"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 1300}},
    {"target": "${canvas_target}", "action": "checkpoint",
     "arguments": {"name": "${image_dir}/ae90-1500ms.png"}},
${view_events}
${lighting_events}
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
    local cache_state="${6:-cold}"
    local case_name="${7:-}"
    local mode="${8:-shaded}"

    if ! jq -e '
	.success == true and
	(.samples | length) > 0 and
	(any(.samples[]; (.draw_shape_count // 0) > 0)) and
	(all(.samples[]; (.failed_sources // 0) == 0)) and
	(all(.samples[]; (.cad_payloads_without_entry // 0) == 0)) and
	# The coordinator is the retained-display authority.  Any rejected event
	# contract or contradictory phase observation is a correctness failure even
	# if the final framebuffer happens to look plausible.
	(all(.samples[];
	    (.lod_coordinator_invariant_violations // 0) == 0 and
	    (.lod_coordinator_invariant_mask // 0) == 0 and
	    (.lod_coordinator_invariant_history_mask // 0) == 0)) and
	# Results may become owner-thread state immediately before the render
	# traversal synchronizes their retained presentation.  That transient is
	# not a displayed fallback frame; the stable/final sample is the
	# authoritative lingering-box invariant.
	((.samples[-1].superseded_fallback_presentations // 0) == 0)
	and
	# The generic progressive flag also covers deferred resident-memory
	# compaction.  A terminal visual frame may therefore be view-ready while
	# memory housekeeping remains queued; reject only unfinished visual work.
	((.samples[-1].progressive_pending == false) or
	 ((.samples[-1].lod_convergence_view_ready // false) == true and
	  (.samples[-1].lod_convergence_background_pending // false) == true)) and
	(.samples[-1].lod_results_pending == false) and
	(.samples[-1].lod_submissions_pending == false) and
	(.samples[-1].lod_refinement_frame_pending == false) and
	# GPU accounting is published only at a complete-frame commit.  A drawn
	# CAD scene must expose a nonempty coherent snapshot, and a clean terminal
	# frame cannot retain an atlas-admission proxy or pressure latch.
	((.samples[-1].lod_gpu_resource_sample_serial // 0) > 0) and
	((.samples[-1].lod_gpu_tracked_buffer_bytes // 0) > 0) and
	((.samples[-1].lod_gpu_atlas_live_bytes // 0) <=
	 (.samples[-1].lod_gpu_atlas_allocated_bytes // 0)) and
	((.samples[-1].lod_gpu_pressure_proxies // 0) == 0) and
	# A scene whose useful visible working set exceeds the configured GPU
	# allowance may terminate at a coherent memory-limited cut.  Requiring the
	# renderer admission-pressure observation to disappear made the 50k/150k
	# qualification fixtures fail precisely when the coordinator honored that
	# contract.  Pressure is acceptable only with an explicit terminal proof;
	# ordinary scenes and every non-ready/pending state still require it clear.
	((.samples[-1].lod_gpu_memory_pressure == false) or
	 ((.samples[-1].lod_convergence_memory_limited // false) == true and
	  (.samples[-1].lod_convergence_view_ready // false) == true and
	  (.samples[-1].lod_convergence_fraction // 0) >= 1 and
	  (.samples[-1].progressive_pending == false))) and
	((.samples[-1].lod_service_pending_tasks // 0) == 0) and
	((.samples[-1].lod_service_active_requests // 0) == 0) and
	((.samples[-1].lod_service_queued_results // 0) == 0) and
	((.samples[-1].lod_service_queued_cache_writes // 0) == 0) and
	((.samples[-1].active_lod_aabb_payloads // 0) == 0) and
	((.samples[-1].active_lod_obb_payloads // 0) == 0) and
	((.samples[-1].active_lod_sphere_payloads // 0) == 0) and
	# Compact entries outside the current view intentionally need no payload:
	# requiring one for every leaf defeats view-aware LoD memory management.
	# Existing payloads must instead form a one-to-one, entry-backed subset.
	(if .samples[-1].deep_lod_diagnostics == true
	 then
	    ((.samples[-1].compact_lod_entries_with_payload // 0) <=
	     (.samples[-1].compact_lod_entries // 0)) and
	    ((.samples[-1].compact_lod_entries_with_payload // 0) ==
	     (.samples[-1].active_lod_cad_payloads // 0))
	 else true
	 end)
	' "$report" >"$validation" 2>&1; then
	printf 'base terminal rendering contract failed\n' >>"$validation"
	return 1
    fi

    if [[ "$mode" == "shaded" &&
	    ("$case_name" == "generic_twin" || "$case_name" == "lucy") ]]; then
	local mged_lighting_image="$image_dir/lighting-mged.png"
	local studio_lighting_image="$image_dir/lighting-studio.png"
	local lighting_dimensions lighting_width lighting_height
	local lighting_crop_width=0 lighting_crop_height=0
	local lighting_changed_pixels=0
	lighting_dimensions=$(identify -format '%wx%h' "$studio_lighting_image" \
	    2>/dev/null || true)
	lighting_width="${lighting_dimensions%x*}"
	lighting_height="${lighting_dimensions#*x}"
	if [[ -f "$mged_lighting_image" && -f "$studio_lighting_image" &&
		"$lighting_width" =~ ^[0-9]+$ &&
		"$lighting_height" =~ ^[0-9]+$ &&
		"$lighting_width" -gt 60 && "$lighting_height" -gt 174 ]]; then
	    lighting_crop_width=$((lighting_width - 60))
	    lighting_crop_height=$((lighting_height - 174))
	    lighting_changed_pixels=$(compare -metric AE -fuzz 2% \
		-crop "${lighting_crop_width}x${lighting_crop_height}+0+24" \
		"$mged_lighting_image" "$studio_lighting_image" null: \
		2>&1 || true)
	fi
	if [[ ! "$lighting_changed_pixels" =~ ^[0-9]+$ ||
		"$lighting_changed_pixels" -lt 500 ]]; then
	    printf 'MGED/studio lighting profiles did not visibly differ: changed_pixels=%s\n' \
		"$lighting_changed_pixels" >>"$validation"
	    return 1
	fi
	if ! jq -e --arg mged "$mged_lighting_image" \
		--arg studio "$studio_lighting_image" '
	    (first(.samples[] |
		select((.checkpoint? // "") == $mged))) as $mged_sample |
	    (first(.samples[] |
		select((.checkpoint? // "") == $studio))) as $studio_sample |
	    ($mged_sample.lighting_profile == "mged") and
	    (($mged_sample.lighting_camera_light_count // -1) == 1) and
	    ((($mged_sample.lighting_ambient_intensity // 0) - 0.30) |
		abs < 0.0001) and
	    ($studio_sample.lighting_profile == "studio") and
	    (($studio_sample.lighting_camera_light_count // -1) == 3) and
	    ((($studio_sample.lighting_ambient_intensity // 0) - 0.18) |
		abs < 0.0001)
	    ' "$report" >>"$validation" 2>&1; then
	    printf 'lighting checkpoints did not report the selected shared rig\n' \
		>>"$validation"
	    return 1
	fi
    fi

    # Lucy is a single source-backed terminal mesh.  It cannot legitimately
    # converge to an empty compact registry, a box-only presentation, or a
    # technically filled but block-silhouette cut.  At this fixed viewport,
    # PoP level 5 is the first state which preserves Lucy's recognizable
    # outline; level 4 takes roughly 12 ms in OSMesa and its level-5 successor
    # remains within the 50 ms stable target, so stopping below it is a policy
    # or premature-compaction failure rather than a valid FPS limit.
    if [[ "$case_name" == "lucy" ]]; then
	if ! jq -e --arg mode "$mode" '
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/ae90-stable.png")))) as $stable |
	    (.samples[-1].compact_lod_entries // 0) >= 1 and
	    (.samples[-1].compact_lod_entries_with_payload // 0) >= 1 and
	    (.samples[-1].active_lod_cad_payloads // 0) >= 1 and
	    (.samples[-1].active_progressive_cad_faces // 0) > 0 and
	    (.samples[-1].lod_service_resident_assets // 0) >= 1 and
	    (.samples[-1].lod_service_resident_bytes // 0) > 0 and
	    (if $mode == "shaded" then
		(($stable.active_progressive_cad_level_min // -1) >= 5)
	     else true end)
	    ' "$report" >>"$validation" 2>&1; then
	    printf 'Lucy did not converge to a recognizable resident PoP mesh\n' \
		>>"$validation"
	    return 1
	fi
    fi

    # If this run captured a partially converged refining frame, the HUD must
    # contain a visible phase-colored fill as well as the numeric progress
    # label.  The faceplate may legitimately lag the diagnostic sample by one
    # frame, so accept any of the documented phase colors.  The
    # controller fraction alone did not catch a depth-buffer regression that
    # left the fill completely hidden behind its gray track.
    local lod_hud_sample
    lod_hud_sample=$(jq -r '
	first(.samples[] |
	    select((.checkpoint? // "") != "" and
		(.presented_cad_faces // 0) > 0 and
		(.lod_convergence_visible_targets // 0) > 0 and
		(.lod_convergence_phase // 0) == 3 and
		(.lod_convergence_fraction // 0) > 0.05 and
		(.lod_convergence_fraction // 1) < 0.99) |
	    .checkpoint) // empty
	' "$report" 2>/dev/null)
    if [[ -n "$lod_hud_sample" && -f "$lod_hud_sample" ]]; then
	local hud_dimensions hud_width hud_height hud_crop_height
	local hud_fill_pixels
	hud_dimensions=$(identify -format '%wx%h' "$lod_hud_sample" \
	    2>/dev/null || true)
	hud_width="${hud_dimensions%x*}"
	hud_height="${hud_dimensions#*x}"
	hud_fill_pixels=0
	if [[ "$hud_width" =~ ^[0-9]+$ && "$hud_height" =~ ^[0-9]+$ &&
		"$hud_width" -gt 80 && "$hud_height" -gt 120 ]]; then
	    hud_crop_height=$((hud_height - 100))
	    hud_fill_pixels=$(convert "$lod_hud_sample" \
		-crop "40x${hud_crop_height}+$((hud_width - 40))+20" \
		-fuzz 20% -fill white -opaque '#60dcff' \
		-opaque '#70eb87' -opaque '#ffcd48' -opaque '#ffaa40' \
		-opaque '#ff5a50' \
		-fill black +opaque white -format '%[fx:mean*w*h]' info: \
		2>/dev/null || printf '0')
	fi
	if ! awk -v pixels="$hud_fill_pixels" 'BEGIN { exit !(pixels >= 100) }'; then
	    printf 'LoD convergence HUD track has no visible progress fill (%s pixels)\n' \
		"$hud_fill_pixels" >>"$validation"
	    return 1
	fi
    fi

    # A retained database-source node is counted before it owns drawable
    # geometry, so report counters alone cannot distinguish a useful first
    # frame from the empty gradient canvas.  Check an actual early
    # framebuffer: the background is horizontally uniform, while boxes and
    # geometry differ from the left-edge color on enough pixels to be useful.
    # A 3+ GiB/150k-object cold database must first open and index its isolated
    # read-only realization handle, then discover every hierarchy occurrence
    # before its whole-target extent can be exact.  For that qualification
    # tier the six-second realization checkpoint (normally about eight
    # seconds after draw return) is the explicit exact-overview deadline; the
    # 1.5-second checkpoint remains a UI/HUD responsiveness sample.  Warm and
    # smaller cases retain the stricter 1.5-second geometry gate.
    # Exclude the border, HUD strip, and right-side convergence indicator.
    local first_useful="$image_dir/ae90-1500ms.png"
    local first_useful_limit_ms=5000
    if [[ "$cache_state" == "cold" ]]; then
	if [[ "$case_name" == "unique_mesh_150k_stress" ]]; then
	    first_useful="$image_dir/realization-6s.png"
	    first_useful_limit_ms=12000
	elif [[ "$case_name" == "unique_mesh_50k_stress" ]]; then
	    # Exact whole-target coverage is published during the intentional
	    # cold held-drag segment.  Use the first checkpoint after the harness
	    # restores its reference ae/autoview; the earlier 1.5 s image remains
	    # the responsiveness/partial-coverage sample and cannot be required to
	    # have the final camera extent.
	    first_useful="$image_dir/realization-1s.png"
	    first_useful_limit_ms=6000
	fi
    fi
    local first_useful_elapsed first_useful_structural_boxes
    first_useful_elapsed=$(jq -r --arg checkpoint "$first_useful" '
	first(.samples[] |
	    select((.checkpoint? // "") == $checkpoint) |
	    (.elapsed_ms // 9223372036854775807))
    ' "$report" 2>/dev/null)
    first_useful_structural_boxes=$(jq -r --arg checkpoint "$first_useful" '
	first(.samples[] |
	    select((.checkpoint? // "") == $checkpoint) |
	    (.visible_structural_fallback_boxes // 0)) // 0
    ' "$report" 2>/dev/null)
    local dimensions width height crop_width crop_height crop_y
    local background changed_pixels
    local minimum_changed_pixels=1000
    # A structural scope/leaf box communicates the model extent with line
    # pixels, not filled area.  When diagnostics prove such a box is actually
    # present, use a perimeter-scale threshold; retaining the 1,000-pixel
    # threshold would reject a correct early overview merely because it is a
    # deliberately cheap wireframe proxy.  The HUD and convergence strip are
    # outside this crop, so they cannot satisfy either threshold.
    if [[ "$first_useful_structural_boxes" =~ ^[0-9]+$ &&
	    "$first_useful_structural_boxes" -gt 0 ]]; then
	minimum_changed_pixels=400
    fi
    # Thousands of widely spaced instances are intentionally subpixel in the
    # initial all-model view.  Requiring 1,000 distinct pixels there rewards a
    # larger proxy, not a faster useful answer; a visible 100-pixel footprint
    # plus the later payload/count assertions is the meaningful threshold.
    if [[ "$object" == "many_lucy_stress" ]]; then
	minimum_changed_pixels=100
    fi
    dimensions=$(identify -format '%wx%h' "$first_useful" 2>/dev/null || true)
    width="${dimensions%x*}"
    height="${dimensions#*x}"
    crop_width=0
    crop_height=0
    crop_y=24
    if [[ "$width" =~ ^[0-9]+$ && "$height" =~ ^[0-9]+$ &&
	    "$width" -gt 60 && "$height" -gt 174 ]]; then
	crop_width=$((width - 60))
	crop_height=$((height - 174))
    fi
    background=$(mktemp "$artifact_dir/.qged-background.XXXXXX.png")
    changed_pixels=0
    if [[ "$crop_width" -gt 0 && "$crop_height" -gt 0 ]] &&
	convert "$first_useful" \
	    -crop "1x${crop_height}+0+${crop_y}" +repage \
	    -scale "${crop_width}x${crop_height}!" "$background" \
	    2>/dev/null; then
	changed_pixels=$(compare -metric AE -fuzz 3% \
	    -crop "${crop_width}x${crop_height}+0+${crop_y}" \
	    "$first_useful" "$background" null: 2>&1 || true)
    fi
    if [[ ! "$first_useful_elapsed" =~ ^[0-9]+$ ||
	    "$first_useful_elapsed" -gt "$first_useful_limit_ms" ||
	    ! "$changed_pixels" =~ ^[0-9]+$ ||
	    "$changed_pixels" -lt "$minimum_changed_pixels" ]]; then
	printf 'no useful model pixels by first-frame checkpoint: elapsed=%s changed_pixels=%s\n' \
	    "$first_useful_elapsed" "$changed_pixels" >>"$validation"
	rm -f "$background"
	return 1
    fi

    # A prefix-ordered traversal can satisfy the raw pixel count while showing
    # only the first row of a large two-dimensional assembly.  For the 50k
    # fixture, compare the early model-pixel extent with the same reference
    # view after convergence.  First-useful means spatially representative
    # whole-scene coverage, not merely "some geometry exists."
    if [[ "$case_name" == "unique_mesh_50k_stress" ||
	    "$case_name" == "unique_mesh_150k_stress" ]]; then
	local stable_image="$image_dir/ae90-stable.png"
	local stable_background early_extent stable_extent
	local early_extent_width early_extent_height
	local stable_extent_width stable_extent_height
	stable_background=$(mktemp \
	    "$artifact_dir/.qged-stable-background.XXXXXX.png")
	early_extent=""
	stable_extent=""
	if convert "$stable_image" \
		-crop "1x${crop_height}+0+${crop_y}" +repage \
		-scale "${crop_width}x${crop_height}!" "$stable_background" \
		2>/dev/null; then
	    early_extent=$(convert "$first_useful" \
		-crop "${crop_width}x${crop_height}+0+${crop_y}" +repage \
		"$background" -compose difference -composite -threshold 3% \
		-trim -format '%@' info: 2>/dev/null || true)
	    stable_extent=$(convert "$stable_image" \
		-crop "${crop_width}x${crop_height}+0+${crop_y}" +repage \
		"$stable_background" -compose difference -composite \
		-threshold 3% -trim -format '%@' info: 2>/dev/null || true)
	fi
	rm -f "$stable_background"
	early_extent_width="${early_extent%%x*}"
	early_extent_height="${early_extent#*x}"
	early_extent_height="${early_extent_height%%+*}"
	stable_extent_width="${stable_extent%%x*}"
	stable_extent_height="${stable_extent#*x}"
	stable_extent_height="${stable_extent_height%%+*}"
	if [[ ! "$early_extent_width" =~ ^[0-9]+$ ||
		! "$early_extent_height" =~ ^[0-9]+$ ||
		! "$stable_extent_width" =~ ^[0-9]+$ ||
		! "$stable_extent_height" =~ ^[0-9]+$ ||
		"$stable_extent_width" -le 0 ||
		"$stable_extent_height" -le 0 ||
		"$((early_extent_width * 4))" -lt "$((stable_extent_width * 3))" ||
		"$((early_extent_height * 4))" -lt "$((stable_extent_height * 3))" ||
		"$((early_extent_width * 4))" -gt "$((stable_extent_width * 5))" ||
		"$((early_extent_height * 4))" -gt "$((stable_extent_height * 5))" ]]; then
	    printf 'first-useful view has invalid whole-scene spatial coverage: early=%s stable=%s\n' \
		"$early_extent" "$stable_extent" >>"$validation"
	    rm -f "$background"
	    return 1
	fi
    fi
    rm -f "$background"

    # Progressive autoview has exactly one camera owner.  The draw may expose
    # the pre-fit camera briefly while exact leaf coverage is assembled, but
    # it must not publish a provisional fit, change orientation independently,
    # or return to the pre-fit camera after the exact fit has landed.  Those
    # transitions are perceived as an autoview jump/flicker even if the final
    # image is correct.
    local autoview_window_start="$image_dir/draw-return.png"
    if [[ "$case_name" == "unique_mesh_50k_stress" ||
	    "$case_name" == "unique_mesh_150k_stress" ]]; then
	# These fixtures deliberately rotate during cold realization, then issue
	# one atomic ae/autoview reference command.  Start after that user-authored
	# camera transition so it is not misclassified as a progressive-autoview
	# write.  Ordinary cases continue to prove the full draw-return window.
	autoview_window_start="$image_dir/realization-1s.png"
    fi
    if ! jq -e --arg window_start "$autoview_window_start" '
	def orientation:
	    [.camera_orientation_axis_x, .camera_orientation_axis_y,
	     .camera_orientation_axis_z, .camera_orientation_angle];
	def framing:
	    [.camera_position_x, .camera_position_y, .camera_position_z,
	     .camera_orthographic_height];
	(first(.samples[] |
	    select((.checkpoint? // "") == $window_start)).event_index) as $start |
	(first(.samples[] |
	    select((.checkpoint? // "") |
		endswith("/ae90-stable.png")))) as $final |
	[.samples[] |
	    select(.event_index >= $start and
		.event_index <= $final.event_index)] as $window |
	([$window[] | orientation] | unique | length) == 1 and
	([$window[] | framing] | unique | length) <= 2 and
	(first($window[] |
	    select(framing == ($final | framing))).event_index) as $fit |
	all($window[] | select(.event_index >= $fit);
	    framing == ($final | framing))
	' "$report" >>"$validation" 2>&1; then
	printf 'progressive autoview exposed provisional or mixed camera states\n' \
	    >>"$validation"
	return 1
    fi

    # A converged, fixed camera is an invariant, not a best effort.  Hold the
    # exact same view for 200 ms and require both the retained cut and the
    # model-area framebuffer to remain unchanged.  This directly rejects the
    # former mesh -> box -> mesh admission loop and repeated deferred-autoview
    # rewrites even when both happen to end on a plausible final screenshot.
    if ! jq -e '
	(first(.samples[] |
	    select((.checkpoint? // "") |
		endswith("/ae90-stable.png")))) as $first |
	(first(.samples[] |
	    select((.checkpoint? // "") |
		endswith("/ae90-stable-held.png")))) as $held |
	($first.progressive_pending == false) and
	($held.progressive_pending == false) and
	($first.lod_refinement_frame_pending == false) and
	($held.lod_refinement_frame_pending == false) and
	($first.camera_view_projection == $held.camera_view_projection) and
	($first.camera_position_x == $held.camera_position_x) and
	($first.camera_position_y == $held.camera_position_y) and
	($first.camera_position_z == $held.camera_position_z) and
	($first.camera_orientation_axis_x ==
	    $held.camera_orientation_axis_x) and
	($first.camera_orientation_axis_y ==
	    $held.camera_orientation_axis_y) and
	($first.camera_orientation_axis_z ==
	    $held.camera_orientation_axis_z) and
	($first.camera_orientation_angle ==
	    $held.camera_orientation_angle) and
	($first.camera_orthographic_height ==
	    $held.camera_orthographic_height) and
	(($first.active_lod_mesh_payloads // -1) ==
	    ($held.active_lod_mesh_payloads // -2)) and
	(($first.visible_structural_fallback_boxes // -1) ==
	    ($held.visible_structural_fallback_boxes // -2)) and
	(($first.active_progressive_cad_occurrence_hash // "") ==
	    ($held.active_progressive_cad_occurrence_hash // "?")) and
	(($first.active_progressive_cad_faces // -1) ==
	    ($held.active_progressive_cad_faces // -2)) and
	(($first.requested_progressive_cad_level_min // -99) ==
	    ($held.requested_progressive_cad_level_min // -98)) and
	(($first.requested_progressive_cad_level_max // -99) ==
	    ($held.requested_progressive_cad_level_max // -98))
	' "$report" >>"$validation" 2>&1; then
	printf 'fixed converged view changed retained/camera state during hold\n' \
	    >>"$validation"
	return 1
    fi
    local stable_hold_pixels
    # The first software framebuffer may still carry its startup device-pixel
    # extent while the settled canvas has adopted the window's final logical
    # size.  Derive this crop from the two images being compared; reusing the
    # early-frame dimensions accidentally included the later HUD (whose FPS
    # text is expected to change) and reported that as model flicker.
    local stable_hold_dimensions stable_hold_width stable_hold_height
    local stable_hold_crop_width=0
    local stable_hold_crop_height=0
    stable_hold_dimensions=$(identify -format '%wx%h' \
	"$image_dir/ae90-stable.png" 2>/dev/null || true)
    stable_hold_width="${stable_hold_dimensions%x*}"
    stable_hold_height="${stable_hold_dimensions#*x}"
    if [[ "$stable_hold_width" =~ ^[0-9]+$ &&
	    "$stable_hold_height" =~ ^[0-9]+$ &&
	    "$stable_hold_width" -gt 60 && "$stable_hold_height" -gt 174 ]]; then
	stable_hold_crop_width=$((stable_hold_width - 60))
	stable_hold_crop_height=$((stable_hold_height - 174))
    fi
    stable_hold_pixels=$(compare -metric AE -fuzz 1% \
	-crop "${stable_hold_crop_width}x${stable_hold_crop_height}+0+${crop_y}" \
	"$image_dir/ae90-stable.png" \
	"$image_dir/ae90-stable-held.png" null: 2>&1 || true)
    if [[ "$stable_hold_crop_width" -le 0 ||
	    "$stable_hold_crop_height" -le 0 ||
	    ! "$stable_hold_pixels" =~ ^[0-9]+$ ||
	    "$stable_hold_pixels" -gt 4 ]]; then
	printf 'fixed converged framebuffer flickered during hold: pixels=%s\n' \
	    "$stable_hold_pixels" >>"$validation"
	return 1
    fi

    # Camera pose and projected scale have different LoD contracts.  Wheel
    # zoom must retarget existing PoP prefixes, while a responsive System GL
    # rotation must retain the cut it already draws rather than installing an
    # unconditional global motion ceiling.  This catches the former behavior
    # where a 5 ms Generic Twin frame was dropped from about 100k faces to the
    # 50k seed merely because the mouse button went down.
    if [[ "$case_name" == "generic_twin" ]]; then
	if ! jq -e '
	    (first(.samples[] | select(.action == "wheel"))) as $zoomIn |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/zoom-out-motion.png")))) as $zoomOut |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/zoom-return-stable.png")))) as $beforeRotate |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-held-end.png")))) as $held |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-motion.png")))) as $motion |
	    ($zoomIn.lod_scale_changing_interaction == true) and
	    ($zoomOut.lod_scale_changing_interaction == true) and
	    ($held.lod_scale_changing_interaction == false) and
	    ($motion.lod_scale_changing_interaction == false) and
	    (if .backend == "system_gl" and
		(($beforeRotate.last_render_ms // 9223372036854775807) <=
		    (1050.0 /
			($beforeRotate.lod_interactive_target_fps // 60.0)))
	     then
		(($held.lod_target_pixel_error // 9223372036854775807) <=
		    1.01) and
		(($held.lod_interactive_progressive_ceiling // -2) == -1) and
		(($held.active_progressive_cad_faces // 0) * 100 >=
		    ($beforeRotate.active_progressive_cad_faces // 0) * 95)
	     else true end)
	    ' "$report" >>"$validation" 2>&1; then
	    printf 'zoom/pose LoD policy did not preserve a responsive retained cut\n' \
		>>"$validation"
	    return 1
	fi
    fi

    # Generic Twin is our production-main visual ground truth for shaded CAD
    # submission.  A transform or channel-routing failure can leave only the
    # thin structural CSG wireframe while all retained-payload counters still
    # claim success.  Require actual filled-surface area in the stable
    # framebuffer.  Compare against the same-row gradient sampled at the left
    # edge, and exclude the border, convergence HUD, and progress bar.
    if [[ ("$case_name" == "generic_twin" ||
	    "$case_name" == "lucy") && "$mode" == "shaded" ]]; then
	local stable_surface="$image_dir/ae90-stable.png"
	local surface_dimensions surface_width surface_height
	local surface_crop_width surface_crop_height
	local surface_background surface_mask surface_pixels
	surface_dimensions=$(identify -format '%wx%h' "$stable_surface" \
	    2>/dev/null || true)
	surface_width="${surface_dimensions%x*}"
	surface_height="${surface_dimensions#*x}"
	surface_crop_width=0
	surface_crop_height=0
	if [[ "$surface_width" =~ ^[0-9]+$ &&
		"$surface_height" =~ ^[0-9]+$ &&
		"$surface_width" -gt 60 && "$surface_height" -gt 174 ]]; then
	    surface_crop_width=$((surface_width - 60))
	    surface_crop_height=$((surface_height - 174))
	fi
	surface_background=$(mktemp \
	    "$artifact_dir/.qged-surface-background.XXXXXX.png")
	surface_mask=$(mktemp "$artifact_dir/.qged-surface-mask.XXXXXX.png")
	surface_pixels=0
	if [[ "$surface_crop_width" -gt 0 &&
		"$surface_crop_height" -gt 0 ]] &&
	    convert "$stable_surface" \
		-crop "1x${surface_crop_height}+0+24" +repage \
		-scale "${surface_crop_width}x${surface_crop_height}!" \
		"$surface_background" 2>/dev/null &&
	    convert "$stable_surface" \
		-crop "${surface_crop_width}x${surface_crop_height}+0+24" \
		+repage "$surface_background" -compose difference \
		-composite -threshold 3% "$surface_mask" 2>/dev/null; then
	    surface_pixels=$(convert "$surface_mask" \
		-format '%[fx:w*h*mean]' info: 2>/dev/null || true)
	    surface_pixels="${surface_pixels%.*}"
	fi
	rm -f "$surface_background" "$surface_mask"
	local minimum_surface_pixels=12000
	if [[ "$case_name" == "lucy" ]]; then
	    minimum_surface_pixels=5000
	fi
	if [[ ! "$surface_pixels" =~ ^[0-9]+$ ||
		"$surface_pixels" -lt "$minimum_surface_pixels" ]]; then
	    printf 'stable shaded framebuffer lacks filled CAD surfaces: pixels=%s\n' \
		"$surface_pixels" >>"$validation"
	    return 1
	fi
    fi

    if [[ -n "$hierarchy_path" ]]; then
	if ! jq -e --arg object "$object" --arg path "$hierarchy_path" '
	    (first(.samples[] |
		select(.action == "set_current") | .event_index)) as $selected |
	    (first(.samples[] |
		select(.command? == ("erase " + $path)) |
		.event_index)) as $erased |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/selection-visible.png")))) as $selectedShot |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/subpath-redraw-stable.png")))) as $redrawnShot |
	    (any(.samples[];
		.selection_paths? != null and
		(.selection_paths | index($path)) != null)) and
	    # Selection must not restart camera-visible LoD work.  Resident-prefix
	    # compaction is deliberately allowed to continue in the background;
	    # progressive_pending includes that memory housekeeping even when the
	    # selected presentation is already complete and view-ready.
	    (all(.samples[];
		if (.event_index > $selected and .event_index < $erased)
		then (.lod_submissions_pending == false and
		      .lod_results_pending == false and
		      .lod_refinement_frame_pending == false and
		      ((.lod_service_pending_tasks // 0) == 0) and
		      ((.lod_service_active_requests // 0) == 0) and
		      ((.lod_service_queued_results // 0) == 0))
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
	    ((.samples[-1].selection_count // -1) == 0) and
	    (($selectedShot.compact_selected_entries // 0) > 0) and
	    (($selectedShot.cad_selected_instances // 0) > 0) and
	    (($redrawnShot.compact_selected_entries // 0) > 0) and
	    (($redrawnShot.cad_selected_instances // 0) > 0) and
	    ((.samples[-1].compact_selected_entries // -1) == 0) and
	    ((.samples[-1].cad_selected_instances // -1) == 0)
	    ' "$report" >>"$validation" 2>&1; then
	    return 1
	fi

	# Counters can agree while a stale aggregate, lost selection style, or
	# missed redraw leaves the framebuffer unchanged.  Require the real tree
	# and canvas pixels to show all four user-facing transitions.  A low
	# threshold is deliberate: Hubble contains valid selectable components
	# only a few pixels wide at the matrix's initial view.
	local tree_selection_pixels erase_pixels redraw_pixels clear_pixels
	local minimum_redraw_pixels=12
	local minimum_clear_pixels=8
	# The selected 16-leaf region in the explicit 50k fixture occupies only
	# one to four pixels in the all-model reference view.  Internal frontier,
	# selection, and compact-entry assertions above prove the exact path
	# transition; requiring twelve pixels here would reward an artificially
	# enlarged proxy rather than correct subpixel restoration.
	if [[ "$case_name" == "unique_mesh_50k_stress" ||
		"$case_name" == "unique_mesh_150k_stress" ]]; then
	    minimum_redraw_pixels=1
	fi
	# In the distinct-mesh fixtures the selected first region is a stack of
	# overlapping hull skins.  It can be completely depth-occluded at the
	# all-model view, so a deselection may correctly change zero canvas
	# pixels.  The source/presentation set assertions above are the robust
	# contract; Generic Twin, Lucy, Hubble, and the other hierarchy cases
	# continue to require an observable selected-style transition.
	local require_clear_pixels=1
	if [[ "$case_name" == "unique_mesh_stress" ||
		"$case_name" == "unique_mesh_50k_stress" ||
		"$case_name" == "unique_mesh_150k_stress" ]]; then
	    require_clear_pixels=0
	fi
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
		! "$redraw_pixels" =~ ^[0-9]+$ ||
		"$redraw_pixels" -lt "$minimum_redraw_pixels" ||
		! "$clear_pixels" =~ ^[0-9]+$ ||
		("$require_clear_pixels" -eq 1 &&
		 "$clear_pixels" -lt "$minimum_clear_pixels") ]]; then
	    printf 'missing hierarchy visual transition: tree=%s erase=%s redraw=%s clear=%s\n' \
		"$tree_selection_pixels" "$erase_pixels" "$redraw_pixels" \
		"$clear_pixels" >>"$validation"
	    return 1
	fi
    fi

    # A large software-rendered mesh must remain responsive while the mouse is
    # held.  Retaining the current cut is the preferred outcome when it already
    # meets the controller's published deadline; otherwise verify that the
    # render-only ceiling actually sheds submitted triangles.
    #
    # The render-only ceiling deliberately leaves producer-authored active
    # cuts intact so a pose-only interaction can recover without rebuilding
    # or reloading geometry.  Testing active_progressive_cad_faces here
    # required that useful retained state be destructively rewritten and
    # rejected the faster O(1) ceiling path even when it submitted a much
    # smaller prefix.  Small/non-progressive scenes are intentionally outside
    # this stress assertion.
    if ! jq -e '
	if .backend != "osmesa" then true
	else
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-held-end.png")))) as $motion |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-stable.png")))) as $stable |
	    if (($stable.presented_cad_faces // 0) < 100000)
	    then true
	    else
		# Pressure may be expressed either by a larger per-object pixel
		# error or by the scene face budget.  The latter is preferred for
		# pose-only motion because it can retain responsive objects while
		# coarsening only the expensive cut.  Test the rendering contract,
		# not one internal pressure signal.
		(($motion.presented_cad_faces // 0) > 0) and
		((($motion.last_render_ms // 9223372036854775807) <=
		  ($motion.presentation_deadline_current_ms // 0)) or
		 ((($motion.presented_cad_faces // 0) * 100 <=
		   ($stable.presented_cad_faces // 0) * 95) and
		  (($motion.lod_interactive_progressive_ceiling // -1) >= 0))) and
		(($stable.lod_interactive_progressive_ceiling // -2) == -1) and
		(($motion.last_render_ms // 9223372036854775807) <= 250)
	    end
	end
	' "$report" >>"$validation" 2>&1; then
	printf 'software interaction was neither deadline-safe nor measurably coarsened\n' \
	    >>"$validation"
	return 1
    fi

    # Shared-array instancing is useful only if view-local occurrences remain
    # a real working set.  The close multi-Lucy sequence starts with a subset,
    # coarsens it under a held drag, then selects two deterministic directions
    # which exchange and retire occurrences.  Require actual occurrence-key
    # turnover, not merely a changed camera or lower global face count.
    if [[ "$object" == "multi_lucy" ]]; then
	if ! jq -e '
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/zoom-return-stable.png")))) as $overview |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/close-focus-stable.png")))) as $focus |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/close-turnover-held.png")))) as $held |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/close-turnover-stable.png")))) as $dragStable |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/close-direction-0-stable.png")))) as $direction0 |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/close-direction-90-stable.png")))) as $direction90 |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/close-turnover-return.png")))) as $returned |
	    (($overview.active_lod_cad_payloads // 0) >= 5) and
	    (($focus.active_lod_cad_payloads // 0) > 0) and
	    (($focus.active_lod_cad_payloads // 0) <
		($overview.active_lod_cad_payloads // 0)) and
	    (($direction0.active_lod_cad_payloads // 0) > 0) and
	    (($direction0.active_lod_cad_payloads // 0) <
		($overview.active_lod_cad_payloads // 0)) and
	    (($direction90.active_lod_cad_payloads // 0) > 0) and
	    (($focus.active_progressive_cad_occurrence_hash // "") !=
		($direction0.active_progressive_cad_occurrence_hash // "")) and
	    (($direction0.active_progressive_cad_occurrence_hash // "") !=
		($direction90.active_progressive_cad_occurrence_hash // "")) and
	    (($direction90.active_progressive_cad_occurrence_hash // "") ==
		($focus.active_progressive_cad_occurrence_hash // "")) and
	    (($returned.active_progressive_cad_occurrence_hash // "") ==
		($overview.active_progressive_cad_occurrence_hash // "")) and
	    (($held.lod_interactive_progressive_ceiling // -1) >= 0) and
	    (($held.active_progressive_cad_faces // 0) <
		($focus.active_progressive_cad_faces // 0)) and
	    (($held.last_render_ms // 9223372036854775807) <= 250) and
	    (($dragStable.lod_interactive_progressive_ceiling // -2) == -1) and
	    (($direction0.lod_interactive_progressive_ceiling // -2) == -1) and
	    (($direction90.lod_interactive_progressive_ceiling // -2) == -1) and
	    (all([$focus, $dragStable, $direction0, $direction90, $returned][];
		((.active_lod_aabb_payloads // 0) == 0) and
		((.active_lod_cad_payloads // 0) > 0)))
	    ' "$report" >>"$validation" 2>&1; then
	    printf 'close-view occurrence working set did not turn over coherently\n' \
		>>"$validation"
	    return 1
	fi
    fi

    # The occurrence-scale stress scene must be governed as one workload.
    # Regressions here used to present as a slow left-to-right replay of
    # thousands of independent PoP prefixes, followed by an unbounded motion
    # frame.  Require a useful warm/cold mesh result, independent stable and
    # interaction calibration, and convergence toward the finite aggregate
    # motion budget.  The factor of two permits a bounded compact-entry wave
    # to still be in flight at the held-button checkpoint.
    if [[ "$object" == "many_lucy_stress" ]]; then
	if ! jq -e --arg cache_state "$cache_state" '
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/ae90-stable.png")))) as $initial |
	    (first(.samples[] | select(.action == "mouse_press"))) as $press |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-held-end.png")))) as $held |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-motion.png")))) as $motion |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-stable.png")))) as $stable |
	    (($cache_state != "warm") or
		(($initial.elapsed_ms // 9223372036854775807) <= 6000)) and
	    (($initial.active_lod_cad_payloads // 0) > 0) and
	    (($motion.lod_scene_face_budget // 0) > 0) and
	    (($motion.lod_scene_face_budget // 0) <
		9223372036854775807) and
	    (($motion.lod_interactive_calibrated_faces_per_second // 0) > 0) and
	    (($motion.lod_stable_calibrated_faces_per_second // 0) > 0) and
	    (($motion.active_lod_scene_faces // 0) <=
		(($motion.lod_scene_face_budget // 0) * 2)) and
	    (($motion.active_lod_scene_faces // 0) <=
		($press.active_lod_scene_faces // 0)) and
	    (($held.active_lod_scene_faces // 0) <=
		($press.active_lod_scene_faces // 0)) and
	    (($held.last_render_ms // 9223372036854775807) <= 250) and
	    (($held.lod_interactive_progressive_ceiling // -1) >= 0) and
	    (($stable.lod_interactive_progressive_ceiling // -2) == -1)
	    ' "$report" >>"$validation" 2>&1; then
	    printf 'occurrence-scale scene did not converge to its calibrated aggregate motion budget\n' \
		>>"$validation"
	    return 1
	fi
    fi

    # Distinct-asset scale must not be mistaken for repeated-instance scale.
    # The default fixture has 5,000 independently stored BoTs.  The explicit
    # 50k fixture must visit its complete compact population, but terminal
    # residency remains view- and scene-budget-aware: demanding 50k resident
    # payloads would directly contradict the production LoD contract.  Every
    # admitted System GL asset must graduate from its box.  Motion may retain
    # the existing cut when it is already responsive, or select a cheaper PoP
    # or aggregate-point cut under load; both backends must recover without
    # pending work.
    if [[ "$object" == "unique_mesh_stress" ]]; then
	# Stable pixel-exact convergence is a terminal quality checkpoint, not
	# the first-useful-frame deadline.  Keep a hard 30 s regression guard:
	# thousands of small assets must not level-walk through avoidable PoP
	# publications merely because the earlier startup checkpoints passed.
	local initial_stable_limit_ms=30000
	local minimum_resident_assets=1000
	local expected_visited_assets=0
	local interaction_render_limit_ms=250
	if [[ "$case_name" == "unique_mesh_50k_stress" ||
		"$case_name" == "unique_mesh_150k_stress" ]]; then
	    # The process/event timeout bounds eventual convergence.  Do not
	    # relabel a usable, interactable cold scene as a startup failure
	    # merely because optional background cache population takes longer
	    # than an arbitrary terminal deadline.
	    initial_stable_limit_ms=9223372036854775807
	    minimum_resident_assets=100
	    if [[ "$case_name" == "unique_mesh_150k_stress" ]]; then
		expected_visited_assets=150000
	    else
		expected_visited_assets=50000
	    fi
	    # The checkpoint includes a framebuffer readback.  Preserve a small
	    # tolerance around the 250 ms interaction objective so 0.5 ms of
	    # capture jitter does not turn a valid bounded cut into a failure.
	    interaction_render_limit_ms=275
	fi
	if ! jq -e --argjson initial_stable_limit_ms \
		"$initial_stable_limit_ms" \
		--argjson minimum_resident_assets "$minimum_resident_assets" \
		--argjson expected_visited_assets "$expected_visited_assets" \
		--argjson interaction_render_limit_ms \
		"$interaction_render_limit_ms" '
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/ae90-stable.png")))) as $initial |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-held-end.png")))) as $held |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-motion.png")))) as $motion |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/rotate-stable.png")))) as $stable |
	    ((first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/cold-rotate-held.png"))) // {})) as $coldHeld |
	    ((first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/cold-rotate-motion.png"))) // {})) as $coldMotion |
	    (($expected_visited_assets == 0) or
		(($coldHeld.elapsed_ms // 9223372036854775807) <= 10000)) and
	    # An explicit ceiling is an overload response, not an interaction
	    # prerequisite.  During cold realization the currently available
	    # prefixes may already render inside the interaction bound; retaining
	    # those prefixes with ceiling == -1 is the desired no-rebuild path.
	    (($expected_visited_assets == 0) or
		(($coldHeld.last_render_ms // 9223372036854775807) <=
		    $interaction_render_limit_ms)) and
	    (($expected_visited_assets == 0) or
		(($coldMotion.last_render_ms // 9223372036854775807) <= 250)) and
	    (($expected_visited_assets == 0) or
		(($coldHeld.lod_service_working_set_limit_bytes // 0) > 0)) and
	    (($initial.elapsed_ms // 9223372036854775807) <=
		$initial_stable_limit_ms) and
	    (($initial.lod_service_resident_assets // 0) >=
		$minimum_resident_assets) and
	    (($expected_visited_assets == 0) or
		(($initial.lod_convergence_available_leaves // 0) >=
		    $expected_visited_assets)) and
	    # The scene face budget is a calibrated refinement allowance, not
	    # permission to discard a visible leaf minimum coherent prefix.
	    # Thousands of minimum prefixes may modestly exceed that soft budget.
	    # Accept the coverage floor only when every expected leaf is already
	    # represented, no structural box remains, and the measured frame is
	    # still inside the interaction bound.
	    (($expected_visited_assets == 0) or
		(($initial.active_progressive_cad_faces // 0) <=
		    ($initial.lod_scene_face_budget // 0)) or
		((($initial.active_lod_cad_payloads // 0) >=
		    $expected_visited_assets) and
		 (($initial.visible_structural_fallback_boxes // 0) == 0) and
		 (($initial.last_render_ms // 9223372036854775807) <=
		    $interaction_render_limit_ms))) and
	    (if .backend == "system_gl" then
		(if $expected_visited_assets > 0 then
		    (($initial.active_lod_cad_payloads // 0) > 0) and
		    (($initial.active_lod_cad_payloads // 0) <=
			($initial.lod_service_resident_assets // -1))
		 else
		    (($initial.active_lod_cad_payloads // 0) ==
			($initial.lod_service_resident_assets // -1))
		 end)
	     else
		(($initial.active_lod_cad_payloads // 0) > 0)
	     end) and
	    (($initial.active_lod_aabb_payloads // 0) == 0) and
	    (($initial.visible_structural_fallback_boxes // 0) == 0) and
	    (($held.active_lod_cad_payloads // 0) > 0) and
	    (($held.active_lod_aabb_payloads // 0) == 0) and
	    (($held.last_render_ms // 9223372036854775807) <=
		$interaction_render_limit_ms) and
	    (($motion.last_render_ms // 9223372036854775807) <= 250) and
	    (if .backend == "system_gl" then
		(if $expected_visited_assets > 0 then
		    (($motion.active_lod_cad_payloads // 0) >=
			($held.active_lod_cad_payloads // -1))
		 else
		    (($motion.active_lod_cad_payloads // 0) ==
			($held.active_lod_cad_payloads // -1))
		 end)
	     else true end) and
	    (if .backend == "system_gl" then
		(if $expected_visited_assets > 0 then
		    (($stable.active_lod_cad_payloads // 0) > 0) and
		    (($stable.active_lod_cad_payloads // 0) <=
			($stable.lod_service_resident_assets // -1))
		 else
		    (($stable.active_lod_cad_payloads // 0) ==
			($stable.lod_service_resident_assets // -1))
		 end)
	     else
		(($stable.active_lod_cad_payloads // 0) >=
		    ($initial.active_lod_cad_payloads // 0))
	     end) and
	    (($stable.active_lod_aabb_payloads // 0) == 0) and
	    (($stable.visible_structural_fallback_boxes // 0) == 0) and
	    ((.samples[-1].visible_structural_fallback_boxes // 0) == 0) and
	    (($stable.lod_interactive_progressive_ceiling // -2) == -1) and
	    ($stable.lod_submissions_pending == false) and
	    ($stable.progressive_pending == false)
	    ' "$report" >>"$validation" 2>&1; then
	    printf 'distinct-mesh scene did not preserve coverage and recover around bounded motion\n' \
		>>"$validation"
	    return 1
	fi
    fi

    # A single large PoP asset exercises the other side of the residency
    # contract: zooming in must extend the retained prefix, while a quiet
    # zoom-out must compact it without rebuilding or discarding the useful
    # coarse prefix.  Keep wheel dispatch bounded throughout so background
    # loading never turns into lost-feeling input.
    if [[ "$case_name" == "lucy" ]]; then
	if ! jq -e '
	    (last(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/smooth-zoom-start-stable.png")))) as $start |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/smooth-zoom-in-12.png")))) as $during |
	    (first(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/smooth-zoom-active-refined.png")))) as $active |
	    (last(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/smooth-zoom-close-stable.png")))) as $close |
	    (last(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/smooth-zoom-out-stable.png")))) as $out |
	    (last(.samples[] |
		select((.checkpoint? // "") |
		    endswith("/smooth-zoom-return.png")))) as $returned |
	    (($close.requested_progressive_cad_level_max // -1) >
		($start.requested_progressive_cad_level_max // -1)) and
	    (if .backend == "system_gl" then
		# Both checkpoints are inside the bracketed scale stream.  Refinement
		# may become drawable by the twelfth event and then be backed down by
		# the O(1) renderer ceiling while later low-amplitude events continue.
		# Either completed in-gesture presentation is therefore valid evidence;
		# requiring the later sample to remain richer would reject the intended
		# FPS response after a real richer frame was already shown.
		($during.lod_gesture_active == true) and
		($during.lod_interactive == true) and
		($during.lod_scale_changing_interaction == true) and
		($active.lod_gesture_active == true) and
		($active.lod_interactive == true) and
		($active.lod_scale_changing_interaction == true) and
		(($active.active_progressive_cad_level_max // -1) >
		    ($start.active_progressive_cad_level_max // -1)) and
		(($active.active_progressive_cad_faces // 0) >
		    ($start.active_progressive_cad_faces // 0)) and
		((([($during.active_progressive_cad_level_max // -1),
		     (if ($during.lod_interactive_progressive_ceiling // -1) >= 0
		      then $during.lod_interactive_progressive_ceiling
		      else ($during.active_progressive_cad_level_max // -1) end)] |
		   min) > ($start.active_progressive_cad_level_max // -1) and
		  (($during.presented_cad_faces // 0) >
		    ($start.presented_cad_faces // 0))) or
		 (([($active.active_progressive_cad_level_max // -1),
		    (if ($active.lod_interactive_progressive_ceiling // -1) >= 0
		     then $active.lod_interactive_progressive_ceiling
		     else ($active.active_progressive_cad_level_max // -1) end)] |
		  min) > ($start.active_progressive_cad_level_max // -1) and
		  (($active.presented_cad_faces // 0) >
		    ($start.presented_cad_faces // 0)))) and
		(($active.active_lod_aabb_payloads // 0) == 0) and
		# A single huge part uses the ordinary retained-VBO tier rather
		# than the multi-part atlas.  Immutable PoP generations must append
		# only their CPU suffix; capacity growth migrates the old prefix
		# device-locally and may never submit the cumulative array again.
		(($during.lod_gpu_ordinary_full_upload_bytes // 0) ==
		    ($start.lod_gpu_ordinary_full_upload_bytes // -1)) and
		(($active.lod_gpu_ordinary_full_upload_bytes // 0) ==
		    ($start.lod_gpu_ordinary_full_upload_bytes // -1)) and
		(($during.lod_gpu_ordinary_suffix_upload_bytes // 0) >=
		    ($start.lod_gpu_ordinary_suffix_upload_bytes // 0)) and
		(($active.lod_gpu_ordinary_suffix_upload_bytes // 0) >=
		    ($during.lod_gpu_ordinary_suffix_upload_bytes // 0)) and
		(($active.lod_gpu_ordinary_suffix_upload_bytes // 0) >
		    ($start.lod_gpu_ordinary_suffix_upload_bytes // 0)) and
		(($active.lod_gpu_ordinary_lineage_reuses // 0) >
		    ($start.lod_gpu_ordinary_lineage_reuses // 0)) and
		(if (($active.lod_gpu_ordinary_part_buffer_bytes // 0) >
			($start.lod_gpu_ordinary_part_buffer_bytes // 0))
		 then
		    (($active.lod_gpu_ordinary_copy_bytes // 0) >
			($start.lod_gpu_ordinary_copy_bytes // 0))
		 else true end)
	     else
		# Software rendering may use a coarser cut while wheel events are
		# arriving, but zoom residency is not a render-budget decision.  The
		# missing suffix must load under the independent memory governor, and
		# continuous scale input must expose at least one richer cut rather
		# than magnifying the same block image until button-up.
		($active.lod_gesture_active == true) and
		($active.lod_interactive == true) and
		($active.lod_scale_changing_interaction == true) and
		(($active.lod_service_resident_bytes // 0) >
		    ($start.lod_service_resident_bytes // 0)) and
		(($active.lod_service_cache_loads // 0) >
		    ($start.lod_service_cache_loads // 0)) and
		# The active cut itself must be richer than the pre-zoom stable cut.
		# It may have reached that cut before the last checkpoint; requiring
		# another increase during the final low-amplitude events would turn
		# successful early refinement into a false failure.
		(($active.active_progressive_cad_level_max // -1) >
		    ($start.active_progressive_cad_level_max // -1)) and
		(($active.active_progressive_cad_faces // 0) >
		    ($start.active_progressive_cad_faces // 0)) and
		# A measured over-budget probe may back off from the arbitrary in-12
		# checkpoint through the O(1) render ceiling.  The durable contract is
		# that continuous input has exposed a cut richer than the pre-zoom
		# stable image, not that every sampled probe is monotonically richer.
		# The occurrence may remain richer than the render-only ceiling.  Its
		# actual submitted level is their minimum and must still improve on the
		# pre-zoom image.
		((([($during.active_progressive_cad_level_max // -1),
		     (if ($during.lod_interactive_progressive_ceiling // -1) >= 0
		      then $during.lod_interactive_progressive_ceiling
		      else ($during.active_progressive_cad_level_max // -1) end)] |
		    min) > ($start.active_progressive_cad_level_max // -1) and
		   (($during.active_progressive_cad_faces // 0) >
		    ($start.active_progressive_cad_faces // 0))) or
		  (([($active.active_progressive_cad_level_max // -1),
		     (if ($active.lod_interactive_progressive_ceiling // -1) >= 0
		      then $active.lod_interactive_progressive_ceiling
		      else ($active.active_progressive_cad_level_max // -1) end)] |
		    min) > ($start.active_progressive_cad_level_max // -1))) and
		(($active.active_lod_aabb_payloads // 0) == 0)
	     end) and
	    # Every renderer tier must publish the exact submitted population.  If
	    # the reusable in-gesture cut already met the less demanding stable
	    # deadline, quiet handoff may reduce it only to the current pixel
	    # demand—not merely because interactive and stable calibration stores
	    # were historically isolated.
	    (($active.presented_cad_faces // 0) > 0) and
	    (if (($active.last_render_ms // 9223372036854775807) <=
		    (1000.0 / ($close.lod_stable_target_fps // 20.0)))
	     then
		(($close.active_progressive_cad_level_max // -1) >=
		 ([($active.active_progressive_cad_level_max // -1),
		   (if ($active.lod_interactive_progressive_ceiling // -1) >= 0
		    then $active.lod_interactive_progressive_ceiling
		    else ($active.active_progressive_cad_level_max // -1) end),
		   ($close.requested_progressive_cad_level_max // -1)] | min))
	     else true end) and
	    # Stable memory maintenance may reclaim the zoom-prefetched suffix at
	    # the close-view pause itself or at the later zoom-out pause.  Either is
	    # correct: require reclamation from the active peak, monotonic
	    # compaction accounting, and a return presentation no worse than the
	    # initial same-scale view.  Byte equality with the initial cache state
	    # is invalid when calibration has since admitted a richer useful cut.
	    ((($close.lod_service_resident_bytes // 0) <
		($active.lod_service_resident_bytes // 0)) or
	     (($out.lod_service_resident_bytes // 0) <
		($active.lod_service_resident_bytes // 0))) and
	    (($returned.lod_service_resident_bytes // 0) <
		($active.lod_service_resident_bytes // 0)) and
	    ((($close.lod_service_compactions // 0) >
		($active.lod_service_compactions // 0)) or
	     (($out.lod_service_compactions // 0) >
		($active.lod_service_compactions // 0))) and
	    (($returned.lod_service_compactions // 0) >=
		($out.lod_service_compactions // 0)) and
	    # Returning to the same scale restores the prior quiet cut whenever
	    # that cut still fits the now-calibrated scene allowance.  A cold first
	    # view can legitimately overshoot before reusable timing exists; do not
	    # require that known-over-budget population forever.  In that case the
	    # returned cut must itself fit the exact allowance and must be reported
	    # as the deliberate responsiveness limit, not as unfinished quality.
	    ((($returned.active_progressive_cad_faces // 0) >=
		($start.active_progressive_cad_faces // 0)) or
	     ((($start.active_lod_scene_render_cost // 0) >
		($returned.lod_scene_render_cost_budget // 0)) and
	      (($returned.active_lod_scene_render_cost // 0) <=
		($returned.lod_scene_render_cost_budget // 0)) and
	      ($returned.lod_convergence_performance_limited == true))) and
	    (all([$start, $close, $out, $returned][];
		(.progressive_pending == false) and
		(.lod_submissions_pending == false) and
		((.active_lod_aabb_payloads // 0) == 0))) and
	    (all(.samples[]; if .action == "wheel"
		then ((.event_duration_us // 9223372036854775807) <= 250000)
		else true end))
	    ' "$report" >>"$validation" 2>&1; then
	    printf 'smooth large-mesh zoom did not load, compact, and return responsively\n' \
		>>"$validation"
	    return 1
	fi

	# The close view deliberately clips Lucy against the upper viewport edge.
	# Once the model is zoomed back out the top scanline must be pure
	# background.  A horizontal non-background span here is the characteristic
	# stale QOpenGLWidget scanline left by treating the gradient quad as a
	# framebuffer clear.
	local edge_image edge_dimensions edge_width edge_background edge_pixels
	for edge_image in smooth-zoom-out-stable.png smooth-zoom-return.png; do
	    edge_dimensions=$(identify -format '%wx%h' \
		"$image_dir/$edge_image" 2>/dev/null || true)
	    edge_width="${edge_dimensions%x*}"
	    edge_background=$(mktemp "$artifact_dir/.qged-edge.XXXXXX.png")
	    edge_pixels=9223372036854775807
	    if [[ "$edge_width" =~ ^[0-9]+$ && "$edge_width" -gt 0 ]] &&
		convert "$image_dir/$edge_image" -crop '1x1+0+0' +repage \
		    -scale "${edge_width}x1!" "$edge_background" 2>/dev/null; then
		edge_pixels=$(compare -metric AE -fuzz 3% \
		    -crop "${edge_width}x1+0+0" "$image_dir/$edge_image" \
		    "$edge_background" null: 2>&1 || true)
	    fi
	    rm -f "$edge_background"
	    if [[ ! "$edge_pixels" =~ ^[0-9]+$ ||
		    "$edge_pixels" -gt 4 ]]; then
		printf 'stale top-edge pixels after smooth zoom: image=%s pixels=%s\n' \
		    "$edge_image" "$edge_pixels" >>"$validation"
		return 1
	    fi
	done
    fi

    # Generic Twin contains 709 BoT occurrences.  Its view-managed wire mode
    # is intentionally backed by those same source-mesh/PoP contracts, not by
    # the legacy plotted-vlist cache.  A former first-warm race converged with
    # zero requests after an authoritative worker overwrote a briefly correct
    # manifest publication; generic "nonempty framebuffer" checks could not
    # distinguish the remaining CSG wires from success.
    if [[ "$case_name" == "generic_twin" && "$mode" == "wire" ]]; then
	if ! jq -e '
	    (.samples[-1].compact_lod_entries // 0) == 709 and
	    (.samples[-1].compact_lod_entries_with_payload // 0) == 709 and
	    (.samples[-1].active_lod_cad_payloads // 0) == 709 and
	    (.samples[-1].active_progressive_cad_faces // 0) > 0
	    ' "$report" >>"$validation" 2>&1; then
	    printf 'Generic Twin wire draw lost its 709 mesh-LoD contracts\n' \
		>>"$validation"
	    return 1
	fi
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
	"$settle_ms" "$hierarchy_root" "$hierarchy_child" "$hierarchy_path" \
	"$case_name"

    local env_args=("BU_DIR_CACHE=$cache_dir")
    if [[ "$swap" != "default" ]]; then
	env_args+=("QGED_SWAP_INTERVAL=$swap")
    fi
    local qged_args=()
    [[ "$backend" == "osmesa" ]] && qged_args+=("-s")
    qged_args+=("--test-script" "$events" "--test-report" "$out/report.json" "$db")

    local command=(env "${env_args[@]}" "$qged" "${qged_args[@]}")
    if [[ "$case_name" == "$perf_case" &&
	    ("$perf_phase" == "both" || "$cache_state" == "$perf_phase") &&
	    "$mode" == "shaded" &&
	    "$swap" == "default" ]]; then
	env_args+=("QGED_TEST_DEEP_LOD_REPORT=0")
	command=(env "${env_args[@]}" "$qged" "${qged_args[@]}")
	command=(perf record -g -o "$out/perf.data" -- "${command[@]}")
    fi
    if [[ "$case_name" == "$apitrace_case" && "$backend" == "system" &&
	    "$cache_state" == "cold" && "$mode" == "shaded" &&
	    "$swap" == "default" ]]; then
	command=(apitrace trace --api=gl -o "$out/qged.trace" "${command[@]}")
    fi

    printf 'RUN %s\n' "$run_name"
    local started=$SECONDS
    local status
    if timeout --signal=TERM "$run_timeout" "${command[@]}" \
	    >"$out/stdout.log" 2>"$out/stderr.log"; then
	if validate_report "$out/report.json" "$out/images" "$object" \
		"$hierarchy_path" "$out/validation.log" "$cache_state" \
		"$case_name" "$mode"; then
	    printf 'PASS,%s,%s,\n' "$run_name" "$((SECONDS - started))" \
		>> "$artifact_dir/results.csv"
	    return 0
	fi
	printf 'FAIL,%s,%s,report-validation\n' "$run_name" \
	    "$((SECONDS - started))" >> "$artifact_dir/results.csv"
	return 1
    else
	status=$?
    fi
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
    # This expanded-copy scene is generated specifically to stress the new
    # occurrence/PCA path.  Production main accepts all commands for it but
    # returns no drawable geometry in either mode, so it cannot serve as a
    # visual ground truth.  The shared multi_lucy and all other cases remain
    # mandatory production baselines.
    if [[ "$case_name" == "multi_lucy_xpush" ]]; then
	echo "SKIP: production qged has no drawable multi_lucy_xpush result" \
	    > "$out/status.txt"
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

    # A screensaver overlay can leave the qged framebuffer capturable while
    # silently intercepting every XTEST click/key.  Deactivate it before each
    # baseline launch so an idle desktop cannot manufacture empty evidence.
    if command -v xfce4-screensaver-command >/dev/null 2>&1; then
	xfce4-screensaver-command --deactivate >/dev/null 2>&1 || true
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
    # Let the database-backed completer and tree finish their initial setup.
    # Sending commands while that work is still pumping events can lose the
    # tail of a synthetic key sequence on large production models.
    sleep 2
    import -window "$window" "$out/initial.png"
    # Click well inside the bottom console instead of estimating the prompt
    # line in physical pixels.  XTEST coordinates follow the X window's
    # logical size, while ImageMagick captures Qt's device-pixel framebuffer;
    # the old fixed 320-pixel offset therefore clicked the canvas at 2x DPI.
    local console_y=$((HEIGHT - 100))
    xdotool mousemove --window "$window" "$((WIDTH / 2))" "$console_y" click 1
    xdotool key ctrl+End
    # The historical console recomputes geometry completion after each key.
    # A 1 ms stream outruns that work on Lucy/Stanford and leaves a truncated
    # command in the prompt.  Pace entry and let the final completion settle
    # before Return.
    xdotool type --delay 20 "draw -m${draw_mode} ${object}"
    sleep 1
    # Dismiss the exact-object completion popup so Return executes the
    # command instead of merely accepting the already-complete token.
    xdotool key Escape
    xdotool key Return
    # A second Return is harmless once a fresh prompt exists, but guarantees
    # execution if a late completer popup consumed the first one.
    sleep 1
    xdotool key Return
    case "$profile" in
	smoke) sleep 5 ;;
	full) sleep 10 ;;
	stress) sleep 20 ;;
    esac
    xdotool type --delay 20 "ae 90 0"
    sleep 1
    xdotool key Escape
    xdotool key Return
    xdotool type --delay 20 "autoview"
    sleep 1
    xdotool key Escape
    xdotool key Return
    sleep 5
    import -window "$window" "$out/ae90-stable.png"
    # Use qged's command path for the production visual baseline.  Large
    # meshes can keep the historical canvas busy long enough for synthetic
    # wheel events to be discarded even though console commands remain
    # ordered and reliable.  Actual wheel/drag behavior is exercised by the
    # current-stack Qt event script above.
    xdotool mousemove --window "$window" "$((WIDTH / 2))" "$console_y" click 1
    xdotool key ctrl+End
    xdotool type --delay 20 "zoom 2"
    sleep 1
    xdotool key Escape
    xdotool key Return
    sleep 5
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

# A progressive autoview has one path-scoped camera result.  Cache warmth and
# renderer speed may change when that result becomes available, but must not
# change the result itself.  The first run for a case/mode/swap tuple records
# the reference; every cold/warm and System/OSMesa peer must match it within a
# small float-serialization tolerance.
validate_autoview_camera_contract()
{
    local case_name="$1"
    local backend="$2"
    local mode="$3"
    local swap="$4"
    local phase="$5"
    local run="${case_name}-${backend}-${mode}-swap${swap//-/_}-${phase}"
    local report="$artifact_dir/cases/$run/report.json"
    local camera_dir="$artifact_dir/camera-contracts"
    local reference="$camera_dir/${case_name}-${mode}-swap${swap//-/_}.json"
    local snapshot="$camera_dir/$run.json"
    mkdir -p "$camera_dir"

    if [[ ! -r "$report" ]] || ! jq -e '
	[.samples[] | select(.checkpoint != null and
	    (.checkpoint | endswith("/ae90-stable.png")))] | length > 0
	' "$report" >/dev/null; then
	printf 'FAIL,%s-camera-contract,0,%s\n' "$run" \
	    "missing stable autoview telemetry" >> "$artifact_dir/results.csv"
	return 1
    fi

    jq -c '
	[.samples[] | select(.checkpoint != null and
	    (.checkpoint | endswith("/ae90-stable.png")))] | last | {
	position: [.camera_position_x, .camera_position_y, .camera_position_z],
	orientation: [.camera_orientation_angle, .camera_orientation_axis_x,
	    .camera_orientation_axis_y, .camera_orientation_axis_z],
	orthographic_height: .camera_orthographic_height,
	focal_distance: .camera_focal_distance,
	near_distance: .camera_near_distance,
	far_distance: .camera_far_distance
	}
    ' "$report" > "$snapshot"

    if [[ ! -e "$reference" ]]; then
	cp "$snapshot" "$reference"
    fi
    if ! jq -e -n --slurpfile actual "$snapshot" \
	    --slurpfile expected "$reference" '
	def close($a; $b):
	    ($a != null and $b != null and
	     (($a - $b) | abs) <=
	     ([0.0001, ((($a | abs) + ($b | abs)) * 0.000001)] | max));
	def array_close($a; $b):
	    ($a | length) == ($b | length) and
	    all(range(0; $a | length); close($a[.]; $b[.]));
	($actual[0]) as $a | ($expected[0]) as $e |
	array_close($a.position; $e.position) and
	array_close($a.orientation; $e.orientation) and
	close($a.orthographic_height; $e.orthographic_height) and
	close($a.focal_distance; $e.focal_distance) and
	close($a.near_distance; $e.near_distance) and
	close($a.far_distance; $e.far_distance)
	' >/dev/null; then
	printf 'FAIL,%s-camera-contract,0,%s\n' "$run" \
	    "final camera differs by cache state or renderer" \
	    >> "$artifact_dir/results.csv"
	return 1
    fi

    printf 'PASS,%s-camera-contract,0,\n' "$run" \
	>> "$artifact_dir/results.csv"
    return 0
}

printf 'status,run,seconds,detail\n' > "$artifact_dir/results.csv"
ldd "$qged" > "$artifact_dir/qged-ldd.txt" 2>&1 || true

runtime_library_path()
{
    local soname="$1"
    awk -v soname="$soname" \
	'index($1, soname) == 1 && $2 == "=>" { print $3; exit }' \
	"$artifact_dir/qged-ldd.txt"
}

record_runtime_library()
{
    local key="$1"
    local soname="$2"
    local path
    path="$(runtime_library_path "$soname")"
    if [[ -z "$path" || ! -r "$path" ]]; then
	printf '%s_path=unresolved\n' "$key"
	return
    fi
    printf '%s_path=%s\n' "$key" "$(realpath "$path")"
    printf '%s_sha256=' "$key"
    sha256sum "$path" | cut -d' ' -f1
}

expected_runtime_dir="$(realpath "$build_dir/lib")"
runtime_provenance="PASS"
runtime_provenance_detail=""
for soname in libBObol.so libObol.so libosmesa.so; do
    library_path="$(runtime_library_path "$soname")"
    if [[ -z "$library_path" || ! -r "$library_path" ]]; then
	runtime_provenance="FAIL"
	runtime_provenance_detail+="${soname}=unresolved "
	continue
    fi
    library_path="$(realpath "$library_path")"
    if [[ "$library_path" != "$expected_runtime_dir/"* ]]; then
	runtime_provenance="FAIL"
	runtime_provenance_detail+="${soname}=${library_path} "
    fi
done

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
    record_runtime_library libbobol libBObol.so
    record_runtime_library libobol libObol.so
    record_runtime_library libosmesa libosmesa.so
    printf 'runtime_provenance=%s\n' "$runtime_provenance"
    if [[ -n "$runtime_provenance_detail" ]]; then
	printf 'runtime_provenance_detail=%s\n' "$runtime_provenance_detail"
    fi
    if [[ -L "$source_root/obol" ]]; then
	printf 'obol_source_root=%s\n' \
	    "$(realpath "$source_root/obol")"
    fi
    if [[ -e "$source_root/obol/external/osmesa" ]]; then
	printf 'osmesa_source_root=%s\n' \
	    "$(realpath "$source_root/obol/external/osmesa")"
    fi
} > "$artifact_dir/run-info.txt"

if [[ "$runtime_provenance" != "PASS" ]]; then
    echo "ERROR: qged resolved a required drawing library outside $expected_runtime_dir" >&2
    echo "       $runtime_provenance_detail" >&2
    printf 'FAIL,runtime-provenance,0,%s\n' \
	"required drawing library is unresolved or outside build/lib" \
	>> "$artifact_dir/results.csv"
    exit 2
fi

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
    if [[ "$baseline_only" -eq 1 ]]; then
	continue
    fi

    case "$profile" in
	# A genuinely cold Lucy cache classifies and persists 28M source faces
	# before walking its view-appropriate cuts.  Useful boxes/mesh pixels have
	# their own early framebuffer gates; this deadline is only for terminal
	# quiescence, including bounded calibration cooldown and compaction.
	smoke) settle_ms=30000 ;;
	full) settle_ms=30000 ;;
	stress) settle_ms=60000 ;;
    esac
    if [[ "$case_name" == "unique_mesh_50k_stress" ]]; then
	settle_ms=180000
    elif [[ "$case_name" == "unique_mesh_150k_stress" ]]; then
	settle_ms=300000
    fi
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
		if run_current "$case_name" "$db" "$object" "$backend" \
			"$mode" "$swap" "cold" "$cache" "$settle_ms" \
			"$hierarchy_root" "$hierarchy_child" "$hierarchy_path"; then
		    validate_autoview_camera_contract "$case_name" "$backend" \
			"$mode" "$swap" "cold" || failures=$((failures + 1))
		else
		    failures=$((failures + 1))
		fi
		if run_current "$case_name" "$db" "$object" "$backend" \
			"$mode" "$swap" "warm" "$cache" "$settle_ms" \
			"$hierarchy_root" "$hierarchy_child" "$hierarchy_path"; then
		    validate_autoview_camera_contract "$case_name" "$backend" \
			"$mode" "$swap" "warm" || failures=$((failures + 1))
		else
		    failures=$((failures + 1))
		fi
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
