#!/usr/bin/env bash
#
# Deterministic qged view-widget resize qualification.  This is intentionally
# separate from qged_gui_matrix.sh: resize invariants must be exercised with
# LoD both enabled and disabled without applying PoP-specific terminal checks
# to evaluated representations.

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_root="$(cd "$script_dir/../../.." && pwd)"
build_dir="$source_root/.build"
artifact_dir=""
backend_list="system,osmesa"
lod_list="auto,off"
mode_list="0,1,2,3,4,5"
case_list="moss"
run_timeout=180
capture_apng=0
minimum_selection_changed_pixels=256
initial_window_width=1100
initial_window_height=800
# A native configure acknowledgement from startup can arrive after a shorter
# apparently stable interval.  Keep the first scripted geometry authoritative
# until the event stream has crossed the delayed X11 configure window.
initial_geometry_stable_ms=1000
initial_geometry_timeout_ms=5000

usage()
{
    cat <<'EOF'
Usage: qged_resize_matrix.sh [options]

  --build-dir DIR       Current build (default: ./.build)
  --artifact-dir DIR    New results directory (never cleared)
  --backends LIST       system,osmesa (default: both)
  --lod LIST            auto,off (default: both)
  --modes LIST          comma-separated 0..5 (default: all)
  --cases LIST          moss,rook,generic_twin,lucy,nist,havoc,hubble,
                        unique_mesh_50k_stress,unique_mesh_150k_stress
                        (default: moss)
  --capture-apng        Capture every presented frame
  --timeout SECONDS     Per-process timeout (default: 180)

Each run resizes both while realization is active and after convergence, then
exercises minimize/restore, maximize/restore, fullscreen/restore, and a rapid
resize burst.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
	--build-dir) build_dir="$2"; shift 2 ;;
	--artifact-dir) artifact_dir="$2"; shift 2 ;;
	--backends) backend_list="$2"; shift 2 ;;
	--lod) lod_list="$2"; shift 2 ;;
	--modes) mode_list="$2"; shift 2 ;;
	--cases) case_list="$2"; shift 2 ;;
	--capture-apng) capture_apng=1; shift ;;
	--timeout) run_timeout="$2"; shift 2 ;;
	--help|-h) usage; exit 0 ;;
	*) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

build_dir="$(realpath -m "$build_dir")"
if [[ -z "$artifact_dir" ]]; then
    artifact_dir="$build_dir/qged-resize-matrix/$(date +%Y%m%d-%H%M%S)"
fi
artifact_dir="$(realpath -m "$artifact_dir")"
qged="$build_dir/bin/qged"
apng_encoder="$build_dir/bin/qged_apng_encode"
[[ -x "$apng_encoder" ]] || apng_encoder="$build_dir/src/qged/qged_apng_encode"

if [[ ! -x "$qged" ]]; then
    echo "ERROR: qged is required in $build_dir/bin" >&2
    exit 2
fi
for tool in jq identify; do
    command -v "$tool" >/dev/null 2>&1 || {
	echo "ERROR: $tool is required" >&2; exit 2;
    }
done
if [[ "$capture_apng" -eq 1 && ! -x "$apng_encoder" ]]; then
    echo "ERROR: APNG encoder is required: $apng_encoder" >&2
    exit 2
fi
if [[ -e "$artifact_dir" ]]; then
    echo "ERROR: artifact directory already exists: $artifact_dir" >&2
    exit 2
fi
mkdir -p "$artifact_dir"/{cases,caches,events}
printf 'status,run,seconds,reason\n' > "$artifact_dir/results.csv"
if ! ldd "$qged" > "$artifact_dir/qged-ldd.txt" 2>&1; then
    echo "ERROR: qged runtime dependency preflight failed; the build may still be linking" >&2
    sed -n '1,40p' "$artifact_dir/qged-ldd.txt" >&2
    exit 2
fi
if grep -Eq 'not found|file too short|invalid ELF' \
	"$artifact_dir/qged-ldd.txt"; then
    echo "ERROR: qged runtime dependency preflight found an incomplete build" >&2
    grep -E 'not found|file too short|invalid ELF' \
	"$artifact_dir/qged-ldd.txt" >&2
    exit 2
fi

case_spec()
{
    case "$1" in
	moss)
	    printf '%s|%s|%s\n' "$source_root/src/libged/tests/draw/moss.g" all.g '' ;;
	rook)
	    printf '%s|%s|%s\n' "$build_dir/share/db/chess/rook.g" scene.g '' ;;
	generic_twin)
	    printf '%s|%s|%s\n' "$build_dir/share/db/faa/Generic_Twin.g" all \
		'all/0xxx_series' ;;
	lucy)
	    printf '%s|%s|%s\n' "${BOBOL_LUCY_DB:-$build_dir/lucy.g}" all \
		'all/r.stl' ;;
	nist)
	    local db="${BOBOL_NIST_DB:-}"
	    if [[ -z "$db" ]]; then
		db="$(find "$build_dir/share/db" -type f -name 'NIST_MBE_PMI_*.g' -print -quit 2>/dev/null)"
	    fi
	    [[ -n "$db" ]] || return 1
	    printf '%s|%s|%s\n' "$db" Document '' ;;
	havoc)
	    printf '%s|%s|%s\n' "$build_dir/share/db/havoc.g" havoc \
		'havoc/havoc_front' ;;
	hubble)
	    printf '%s|%s|%s\n' \
		"${BOBOL_HUBBLE_DB:-/home/cyapp/models/NASA/Hubble/Hubble_Space_Telescope.g}" \
		'all.g' 'all.g/BODY' ;;
	unique_mesh_50k_stress)
	    printf '%s|%s|%s\n' \
		"${BOBOL_UNIQUE_MESH_50K_STRESS_DB:-$build_dir/unique_mesh_50k_stress.g}" \
		unique_mesh_stress '' ;;
	unique_mesh_150k_stress)
	    printf '%s|%s|%s\n' \
		"${BOBOL_UNIQUE_MESH_150K_STRESS_DB:-$build_dir/unique_mesh_150k_stress.g}" \
		unique_mesh_stress '' ;;
	*) return 1 ;;
    esac
}

write_events()
{
    local file="$1" image_dir="$2" object="$3" mode="$4" lod="$5"
    local scope_timeout_ms="$6" settle_timeout_ms="$7"
    local hierarchy_path="${8:-}"
    local lod_value=1
    local hierarchy_select_events=""
    local hierarchy_finish_events=""
    local opposite_lod_value=0
    local opposite_lod_name="disabled"
    local restored_lod_name="reenabled"
    [[ "$lod" == off ]] && lod_value=0
    if [[ "$lod" == off ]]; then
	opposite_lod_value=1
	opposite_lod_name="enabled"
	restored_lod_name="restored-off"
    fi
    if [[ -n "$hierarchy_path" ]]; then
	local -a hierarchy_labels=()
	local hierarchy_selection_labels=""
	local hierarchy_expand_events=""
	local prefix_labels=""
	local label
	local label_index
	IFS='/' read -r -a hierarchy_labels <<< "$hierarchy_path"
	for ((label_index = 0;
		label_index < ${#hierarchy_labels[@]}; ++label_index)); do
	    label="${hierarchy_labels[$label_index]}"
	    [[ -n "$hierarchy_selection_labels" ]] &&
		hierarchy_selection_labels+=", "
	    hierarchy_selection_labels+="\"${label}\""
	    [[ -n "$prefix_labels" ]] && prefix_labels+=", "
	    prefix_labels+="\"${label}\""
	    if ((label_index + 1 < ${#hierarchy_labels[@]})); then
		hierarchy_expand_events+="  {\"target\":\"./n:Hierarchy/i:hierarchy-tree\",\"action\":\"set_expanded\",\"arguments\":{\"labels\":[${prefix_labels}],\"expanded\":true}},
"
	    fi
	done
	hierarchy_select_events="${hierarchy_expand_events}  {\"target\":\"./n:Hierarchy/i:hierarchy-tree\",\"action\":\"set_current\",\"arguments\":{\"labels\":[${hierarchy_selection_labels}]}},
  {\"target\":\".\",\"action\":\"wait\",\"arguments\":{\"ms\":50}},
  {\"target\":\"./n:Hierarchy/i:hierarchy-tree\",\"action\":\"checkpoint\",\"arguments\":{\"name\":\"${image_dir}/tree_selected.png\"}},
  {\"target\":\"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas\",\"action\":\"checkpoint\",\"arguments\":{\"name\":\"${image_dir}/selected_before_resize.png\"}},"
	hierarchy_finish_events="  {\"target\":\"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas\",\"action\":\"checkpoint\",\"arguments\":{\"name\":\"${image_dir}/selected_after_storm.png\"}},
  {\"target\":\".\",\"action\":\"wait_progressive_idle\",\"arguments\":{\"timeout_ms\":${settle_timeout_ms},\"quiet_ms\":100}},
  {\"target\":\"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas\",\"action\":\"checkpoint\",\"arguments\":{\"name\":\"${image_dir}/selected_stable.png\"}},
  {\"target\":\"./n:Hierarchy/i:hierarchy-tree\",\"action\":\"clear_selection\",\"arguments\":{}},
  {\"target\":\".\",\"action\":\"wait\",\"arguments\":{\"ms\":50}},"
    fi
    cat > "$file" <<EOF
{
 "schema":"brlcad.qtcad.events",
 "version":1,
 "events":[
  {"target":".","action":"wait_canvas_ready","arguments":{"timeout_ms":10000}},
  {"target":".","action":"resize","arguments":{"width":${initial_window_width},"height":${initial_window_height},"stable_ms":${initial_geometry_stable_ms},"timeout_ms":${initial_geometry_timeout_ms}}},
  {"target":".","action":"wait","arguments":{"ms":100}},
  {"target":".","action":"qged_command_batch","arguments":{"commands":["view lod ${lod_value}","draw -m${mode} ${object}","ae 90 0","autoview"]}},
  {"target":".","action":"wait","arguments":{"ms":100}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/initial.png"}},
  {"target":".","action":"resize","arguments":{"width":760,"height":620}},
  {"target":".","action":"wait","arguments":{"ms":80}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/small.png"}},
  {"target":".","action":"resize","arguments":{"width":1420,"height":900}},
  {"target":".","action":"wait","arguments":{"ms":80}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/large.png"}},
  {"target":".","action":"window_state","arguments":{"state":"minimized"}},
  {"target":".","action":"wait","arguments":{"ms":80}},
  {"target":".","action":"window_state","arguments":{"state":"normal"}},
  {"target":".","action":"wait","arguments":{"ms":100}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/restored_early.png"}},
  {"target":".","action":"wait_progressive_scope_ready","arguments":{"timeout_ms":${scope_timeout_ms},"quiet_ms":100}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/settled_base.png"}},
${hierarchy_select_events}
  {"target":".","action":"window_state","arguments":{"state":"maximized"}},
  {"target":".","action":"wait","arguments":{"ms":100}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/maximized.png"}},
  {"target":".","action":"window_state","arguments":{"state":"normal"}},
  {"target":".","action":"wait","arguments":{"ms":100}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/restored_max.png"}},
  {"target":".","action":"window_state","arguments":{"state":"fullscreen"}},
  {"target":".","action":"wait","arguments":{"ms":100}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/fullscreen.png"}},
  {"target":".","action":"window_state","arguments":{"state":"normal"}},
  {"target":".","action":"wait","arguments":{"ms":100}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/restored_full.png"}},
  {"target":".","action":"resize","arguments":{"width":821,"height":577}},
  {"target":".","action":"resize","arguments":{"width":1399,"height":877}},
  {"target":".","action":"resize","arguments":{"width":703,"height":511}},
  {"target":".","action":"resize","arguments":{"width":1283,"height":741}},
  {"target":".","action":"resize","arguments":{"width":947,"height":693}},
  {"target":".","action":"wait","arguments":{"ms":120}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/storm.png"}},
${hierarchy_finish_events}
  {"target":".","action":"wait_progressive_idle","arguments":{"timeout_ms":${settle_timeout_ms},"quiet_ms":100}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/stable.png"}},
  {"target":".","action":"qged_command_batch","arguments":{"commands":["view lod ${opposite_lod_value}"]}},
  {"target":".","action":"wait_progressive_idle","arguments":{"timeout_ms":${settle_timeout_ms},"quiet_ms":100}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/lod-${opposite_lod_name}.png"}},
  {"target":".","action":"qged_command_batch","arguments":{"commands":["view lod ${lod_value}"]}},
  {"target":".","action":"wait_progressive_idle","arguments":{"timeout_ms":${settle_timeout_ms},"quiet_ms":100}},
  {"target":"./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas","action":"checkpoint","arguments":{"name":"${image_dir}/lod-${restored_lod_name}.png"}}
 ]
}
EOF
}

validate_run()
{
    local report="$1" image_dir="$2" log="$3" lod="$4" case_name="$5" mode="$6"
    local hierarchy_path="${7:-}"
    : > "$log"
    echo "validate: report" >>"$log"
    jq -e '.success == true and (.samples | length) >= 20' "$report" \
	>>"$log" 2>&1 || return 1
    # Every resize must be reflected in the top-level widget immediately; the
    # following wait/checkpoint must have a positive canvas whose controller
    # viewport matches its physical pixel size (within integer DPR rounding).
    echo "validate: resize-state" >>"$log"
    jq -e --arg lod "$lod" \
	--argjson initial_width "$initial_window_width" \
	--argjson initial_height "$initial_window_height" '
      def okdims:
        (.controller_available == true) and
        (.canvas_width > 0 and .canvas_height > 0) and
        (.viewport_width > 0 and .viewport_height > 0) and
        ((.viewport_width - (.canvas_width * .canvas_device_pixel_ratio)) | abs) <= 2 and
        ((.viewport_height - (.canvas_height * .canvas_device_pixel_ratio)) | abs) <= 2;
      (first(.samples[] | select(.action == "qged_command_batch" and
	any((.commands // [])[]; startswith("draw -m"))))) as $draw |
      ($draw.window_width == $initial_width and
	$draw.window_height == $initial_height) and
      all(.samples[] | select(.action == "resize");
        .window_width == .requested_width and
        .window_height == .requested_height) and
      all(.samples[] | select(.action == "window_state");
	if .requested_window_state == "minimized" then
	  .window_minimized == true
	elif .requested_window_state == "maximized" then
	  .window_maximized == true
	elif .requested_window_state == "fullscreen" then
	  .window_fullscreen == true
	else
	  .window_minimized == false and
	  .window_maximized == false and
	  .window_fullscreen == false
	end) and
      all(.samples[] | select(.action == "checkpoint");
	okdims and
	(if ((.checkpoint? // "") | endswith("/lod-disabled.png") or
	    endswith("/lod-restored-off.png")) then
	  .view_lod_policy == 0 and
	  .view_lod_mesh_enabled == false and
	  .view_lod_csg_enabled == false and
	  .lod_progress_track_present == false and
	  .lod_progress_fill_present == false and
	  .lod_progress_label_present == false
	elif ((.checkpoint? // "") | endswith("/lod-enabled.png") or
	      endswith("/lod-reenabled.png")) then
	  .view_lod_policy == 1 and
	  .view_lod_mesh_enabled == true and
	  .view_lod_csg_enabled == true
	elif $lod == "off" then
	  .view_lod_policy == 0 and
	  .view_lod_mesh_enabled == false and
	  .view_lod_csg_enabled == false and
	  .lod_progress_track_present == false and
	  .lod_progress_fill_present == false and
	  .lod_progress_label_present == false
	 else
	  .view_lod_policy == 1 and
	  .view_lod_mesh_enabled == true and
	  .view_lod_csg_enabled == true
	 end)) and
      ([.samples[] | select((.checkpoint? // "") | endswith("/settled_base.png"))][0]) as $first |
      all(.samples[] | select((.action == "resize" or .action == "checkpoint") and
	.event_index > $first.event_index);
        ((.camera_position_x - $first.camera_position_x) | abs) < 1e-4 and
        ((.camera_position_y - $first.camera_position_y) | abs) < 1e-4 and
        ((.camera_position_z - $first.camera_position_z) | abs) < 1e-4 and
	((.camera_orientation_angle - $first.camera_orientation_angle) | abs) < 1e-5 and
	((.model_view_size - $first.model_view_size) | abs) <
	  (1e-4 + (($first.model_view_size | abs) * 1e-6)) and
	.draw_scene_revision == $first.draw_scene_revision) and
      ([.samples[] | select((.checkpoint? // "") | endswith("/stable.png"))][0]) as $terminal |
      (($terminal.lod_convergence_phase != 0 or
	$terminal.lod_convergence_background_pending == true or
	$terminal.lod_convergence_performance_limited == true or
	($terminal.failed_sources // 0) > 0) as $terminal_hud_expected |
       $terminal.progressive_pending == false and
	$terminal.lod_progress_track_present == $terminal_hud_expected and
	$terminal.lod_progress_fill_present == $terminal_hud_expected and
	$terminal.lod_progress_label_present == $terminal_hud_expected and
	$terminal.visible_structural_fallback_boxes == 0)
    ' "$report" >>"$log" 2>&1 || return 1
    # Policy changes operate on the retained draw; they must not duplicate or
    # replace its semantic root, move the camera, or strand startup proxies.
    # This covers both managed -> full -> managed and full -> managed -> full
    # transitions using the same initially drawn scene.
    echo "validate: policy-toggle" >>"$log"
    jq -e --arg lod "$lod" '
      ([.samples[] | select((.checkpoint? // "") |
	endswith("/stable.png"))][0]) as $stable |
      (first(.samples[] | select(.action == "qged_command_batch" and
	any((.commands // [])[]; startswith("draw -m"))))) as $draw |
      ([.samples[] | select(.action == "qged_command_batch" and
	((.commands // []) | length) == 1 and
	((.commands[0] // "") | startswith("view lod ")))][0]) as $opposite_command |
      ([.samples[] | select(.action == "qged_command_batch" and
	((.commands // []) | length) == 1 and
	((.commands[0] // "") | startswith("view lod ")))][-1]) as $restored_command |
      ([.samples[] | select((.checkpoint? // "") |
	test("/lod-(disabled|enabled)\\.png$"))][0]) as $opposite |
      ([.samples[] | select((.checkpoint? // "") |
	test("/lod-(reenabled|restored-off)\\.png$"))][0]) as $restored |
      def same_camera($a; $b):
	(($a.camera_position_x - $b.camera_position_x) | abs) < 1e-4 and
	(($a.camera_position_y - $b.camera_position_y) | abs) < 1e-4 and
	(($a.camera_position_z - $b.camera_position_z) | abs) < 1e-4 and
	(($a.model_view_size - $b.model_view_size) | abs) <
	  (1e-4 + (($a.model_view_size | abs) * 1e-6));
      ($draw.draw_frontier_count > 0) and
      ($opposite_command.draw_frontier_count == $draw.draw_frontier_count) and
      ($restored_command.draw_frontier_count == $draw.draw_frontier_count) and
      ($opposite_command.draw_frontier_paths == $draw.draw_frontier_paths) and
      ($restored_command.draw_frontier_paths == $draw.draw_frontier_paths) and
      ($opposite_command.draw_scene_revision == $draw.draw_scene_revision) and
      ($restored_command.draw_scene_revision == $draw.draw_scene_revision) and
      ($opposite.draw_scene_revision == $stable.draw_scene_revision) and
      ($restored.draw_scene_revision == $stable.draw_scene_revision) and
      same_camera($stable; $opposite) and same_camera($stable; $restored) and
      ($opposite.visible_structural_fallback_boxes == 0) and
      ($restored.visible_structural_fallback_boxes == 0) and
      (($opposite.failed_sources // 0) == 0) and
      (($restored.failed_sources // 0) == 0) and
      (if $opposite.view_lod_policy == 0 then
	$opposite.lod_progress_track_present == false and
	$opposite.lod_progress_fill_present == false and
	$opposite.lod_progress_label_present == false
       else true end) and
      (if $restored.view_lod_policy == 0 then
	$restored.lod_progress_track_present == false and
	$restored.lod_progress_fill_present == false and
	$restored.lod_progress_label_present == false
       else true end) and
      (if $lod == "off" then
	$opposite.view_lod_policy == 1 and $restored.view_lod_policy == 0
       else
	$opposite.view_lod_policy == 0 and $restored.view_lod_policy == 1
       end)
    ' "$report" >>"$log" 2>&1 || return 1
    if [[ "$lod" == auto &&
	("$case_name" == generic_twin || "$case_name" == lucy) &&
	("$mode" == 0 || "$mode" == 1 || "$mode" == 4) ]]; then
	# A fresh cache normally leaves one of the early checkpoints in a
	# producer-active state.  Do not make that timing an invariant: fast
	# hardware or a small source can complete before the first checkpoint.  In
	# that case the terminal managed population is stronger evidence that this
	# exercised LoD rather than silently falling back to a synchronous path.
	echo "validate: generic-early-or-settled" >>"$log"
	jq -e '
	  def early_progressive:
	    any(.samples[] | select(.action == "checkpoint" and
	      ((.checkpoint? // "") | test("/(initial|small|large)\\.png$")));
	      .progressive_pending == true or
	      .source_realization_active_items > 0 or
	      .lod_submissions_pending == true or
	      .lod_results_pending == true);
	  ([.samples[] | select((.checkpoint? // "") |
	    endswith("/stable.png"))][0]) as $terminal |
	  early_progressive or
	  ($terminal.view_lod_policy == 1 and
	   $terminal.active_lod_cad_payloads > 0 and
	   $terminal.visible_structural_fallback_boxes == 0 and
	   $terminal.lod_convergence_view_ready == true)
	' "$report" >>"$log" 2>&1 || return 1
	# "Ready" means the physical view request is satisfied, unless a named
	# resource limit explains why it is not.  In particular, interaction
	# hysteresis must not silently become the quiet terminal pixel target.
	echo "validate: generic-pixel-target" >>"$log"
	jq -e '
	  ([.samples[] | select((.checkpoint? // "") |
	    endswith("/stable.png"))][0]) as $terminal |
	  ($terminal.lod_convergence_view_ready == true) and
	  (if (($terminal.lod_convergence_performance_limited == false) and
	       ($terminal.lod_convergence_memory_limited == false) and
	       ($terminal.lod_projected_demand_cad_payloads > 0)) then
	     $terminal.lod_max_cad_normalized_error <= 1.000001
	   else true end)
	' "$report" >>"$log" 2>&1 || return 1
    fi
    # A modest single-object BREP is an important backend parity oracle.  It
    # cannot legitimately finish as a recognizable-quality-floor violation;
    # that state previously exposed an adaptive-band handoff which replaced a
    # good mesh with the new hierarchy's minimum 20-face prefix after resize.
    if [[ "$lod" == auto && "$case_name" == nist && "$mode" == 1 ]]; then
	jq -e '
	  ([.samples[] | select((.checkpoint? // "") |
	    endswith("/stable.png"))][0]) as $terminal |
	  ($terminal.lod_convergence_view_ready == true) and
	  ($terminal.active_lod_cad_payloads > 0) and
	  ($terminal.active_lod_proxy_payloads == 0) and
	  ($terminal.lod_prominent_cad_quality_floor_violations == 0)
	' "$report" >>"$log" 2>&1 || return 1
    fi
    # The evaluated-points display is intentionally a randomized sample, but
    # its camera is not.  Autoview must use scene.g's authoritative 73.3892 mm
    # maximum object extent (the retained framing contract reports twice that
    # extent), never the accidental bounds of whichever samples were emitted.
    # The old sample-bounds bug produced size 96.7, an offset center, and a
    # visibly clipped rook while every generic resize invariant still passed.
    if [[ "$case_name" == rook && "$mode" == 5 ]]; then
	jq -e '
	  ([.samples[] | select((.checkpoint? // "") |
	    endswith("/stable.png"))][0]) as $terminal |
	  (($terminal.model_view_size - 146.7784) | abs) < 0.02 and
	  (($terminal.camera_position_x // 999) | abs) < 0.02 and
	  (($terminal.camera_position_z - 32.3958) | abs) < 0.02
	' "$report" >>"$log" 2>&1 || return 1
    fi
    if [[ -n "$hierarchy_path" ]]; then
	jq -e --arg path "$hierarchy_path" '
	  (. as $root |
	   first($root.samples[] | select(.action == "set_current") |
	    .event_index)) as $selected |
	  (. as $root |
	   first($root.samples[] | select(.action == "clear_selection") |
	    .event_index)) as $cleared |
	  ($selected < $cleared) and
	  (any(.samples[];
	    .event_index == $selected and
	    (.selection_paths? != null) and
	    (.selection_paths | index($path)) != null)) and
	  (all(.samples[] |
	    select(.event_index > $selected and .event_index < $cleared and
	      .action == "checkpoint");
	    ((.cad_selected_instances // 0) > 0))) and
	  (any(.samples[];
	    .event_index == $cleared and
	    ((.selection_paths // []) | length) == 0 and
	    ((.cad_selected_instances // 0) == 0)))
	' "$report" >>"$log" 2>&1 || return 1
    fi
    local transition_images="lod-disabled lod-reenabled"
    [[ "$lod" == off ]] && transition_images="lod-enabled lod-restored-off"
    for image in initial small large restored_early settled_base maximized \
	restored_max fullscreen restored_full storm stable $transition_images; do
	local file="$image_dir/$image.png"
	[[ -s "$file" ]] || { echo "missing $file" >>"$log"; return 1; }
	identify "$file" >>"$log" 2>&1 || return 1
    done
    if [[ -n "$hierarchy_path" ]]; then
	for image in tree_selected selected_before_resize selected_after_storm \
	    selected_stable; do
	    local file="$image_dir/$image.png"
	    [[ -s "$file" ]] || {
		echo "missing $file" >>"$log"; return 1;
	    }
	    identify "$file" >>"$log" 2>&1 || return 1
	done
	local selection_delta
	selection_delta="$(compare -metric AE -fuzz 2% \
	    "$image_dir/selected_stable.png" "$image_dir/stable.png" \
	    null: 2>&1 || true)"
	echo "stable selection changed pixels: $selection_delta" >>"$log"
	[[ "$selection_delta" =~ ^[0-9]+$ ]] &&
	    ((selection_delta >= minimum_selection_changed_pixels)) || return 1
    fi
}

failures=0
IFS=',' read -r -a cases <<< "$case_list"
IFS=',' read -r -a backends <<< "$backend_list"
IFS=',' read -r -a lods <<< "$lod_list"
IFS=',' read -r -a modes <<< "$mode_list"
for case_name in "${cases[@]}"; do
    spec="$(case_spec "$case_name")" || {
	echo "ERROR: unknown or unavailable case: $case_name" >&2; exit 2;
    }
    IFS='|' read -r db object hierarchy_path <<< "$spec"
    [[ -f "$db" ]] || { echo "ERROR: database not found: $db" >&2; exit 2; }
    for backend in "${backends[@]}"; do
	[[ "$backend" == system || "$backend" == osmesa ]] || exit 2
	for lod in "${lods[@]}"; do
	    [[ "$lod" == auto || "$lod" == off ]] || exit 2
	    for mode in "${modes[@]}"; do
		[[ "$mode" =~ ^[0-5]$ ]] || exit 2
		run="${case_name}-${backend}-m${mode}-lod${lod}"
		out="$artifact_dir/cases/$run"
		cache="$artifact_dir/caches/$run"
		events="$artifact_dir/events/$run.json"
		mkdir -p "$out/images" "$cache"
		scope_timeout_ms=60000
		settle_timeout_ms=60000
		if [[ "$case_name" == unique_mesh_50k_stress ]]; then
		    scope_timeout_ms=120000
		    settle_timeout_ms=180000
		elif [[ "$case_name" == unique_mesh_150k_stress ]]; then
		    scope_timeout_ms=180000
		    settle_timeout_ms=300000
		fi
		write_events "$events" "$out/images" "$object" "$mode" "$lod" \
		    "$scope_timeout_ms" "$settle_timeout_ms" "$hierarchy_path"
		args=(--test-script "$events" --test-report "$out/report.json" "$db")
		[[ "$backend" == osmesa ]] && args=(-s "${args[@]}")
		env_args=("BU_DIR_CACHE=$cache")
		if [[ "$capture_apng" -eq 1 ]]; then
		    mkdir -p "$out/frames"
		    env_args+=("QGED_TEST_FRAME_DIR=$out/frames")
		fi
		started=$SECONDS
		if timeout --signal=TERM "$run_timeout" env "${env_args[@]}" \
			"$qged" "${args[@]}" >"$out/stdout.log" \
			2>"$out/stderr.log" &&
		    validate_run "$out/report.json" "$out/images" \
			"$out/validation.log" "$lod" "$case_name" "$mode" \
			"$hierarchy_path"; then
		    if [[ "$capture_apng" -eq 1 ]]; then
			"$apng_encoder" "$out/frames/frames.tsv" \
			    "$out/presented.apng" --pad-to-max --remove-inputs \
			    >"$out/apng.stdout.log" 2>"$out/apng.stderr.log" || {
			    printf 'FAIL,%s,%s,apng\n' "$run" "$((SECONDS-started))" >>"$artifact_dir/results.csv"
			    failures=$((failures+1)); continue;
			}
		    fi
		    printf 'PASS,%s,%s,\n' "$run" "$((SECONDS-started))" >>"$artifact_dir/results.csv"
		else
		    status=$?
		    printf 'FAIL,%s,%s,status=%s\n' "$run" "$((SECONDS-started))" "$status" >>"$artifact_dir/results.csv"
		    failures=$((failures+1))
		fi
	    done
	done
    done
done

# The camera is a property of the draw target, not of its current display
# representation.  When both policies are requested, compare the terminal
# stable checkpoints so a display-derived union cannot silently make LoD-off
# frame a different model extent than managed LoD.  Position tolerance is
# scaled to the view size to accommodate float camera storage; the size itself
# is held to a much tighter relative tolerance.
if [[ ",$lod_list," == *,auto,* && ",$lod_list," == *,off,* ]]; then
    for case_name in "${cases[@]}"; do
	for backend in "${backends[@]}"; do
	    for mode in "${modes[@]}"; do
		auto_report="$artifact_dir/cases/${case_name}-${backend}-m${mode}-lodauto/report.json"
		off_report="$artifact_dir/cases/${case_name}-${backend}-m${mode}-lodoff/report.json"
		[[ -s "$auto_report" && -s "$off_report" ]] || continue
		if ! jq -e -n --slurpfile auto "$auto_report" \
			--slurpfile off "$off_report" '
                  def stable($r):
                    first($r[0].samples[] |
                      select(.action == "checkpoint" and
                        ((.checkpoint // "") | endswith("/stable.png"))));
                  def finite: type == "number" and . == . and
                    (. > -1.0e300 and . < 1.0e300);
                  stable($auto) as $a | stable($off) as $b |
                  ($a != null and $b != null) and
                  ($a.model_view_size | finite) and
                  ($b.model_view_size | finite) and
                  (($a.model_view_size - $b.model_view_size) | fabs) <=
                    (1.0e-6 * ([1.0, ($a.model_view_size | fabs),
                      ($b.model_view_size | fabs)] | max)) and
                  ([$a.camera_position_x, $a.camera_position_y,
                    $a.camera_position_z, $b.camera_position_x,
                    $b.camera_position_y, $b.camera_position_z] |
                    all(finite)) and
                  ([($a.camera_position_x - $b.camera_position_x) | fabs,
                    ($a.camera_position_y - $b.camera_position_y) | fabs,
                    ($a.camera_position_z - $b.camera_position_z) | fabs] |
                    max) <= (1.0e-6 *
                      ([1.0, ($a.model_view_size | fabs)] | max))
                ' >"$artifact_dir/cases/${case_name}-${backend}-m${mode}-policy-camera.log" \
			2>&1; then
		    run="${case_name}-${backend}-m${mode}-policy-camera"
		    printf 'FAIL,%s,0,LoD policy changed stable camera\n' "$run" \
			>>"$artifact_dir/results.csv"
		    failures=$((failures+1))
		fi
	    done
	done
    done
fi

exit "$failures"
