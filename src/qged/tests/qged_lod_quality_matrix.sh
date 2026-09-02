#!/usr/bin/env bash
#
# Compare settled qged LoD images with full-detail images produced by the same
# drawing stack.  Startup, resize, and interaction contracts belong to their
# dedicated matrices; this script isolates terminal visual quality.

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_root="$(cd "$script_dir/../../.." && pwd)"
build_dir="$source_root/.build"
artifact_dir=""
cache_root=""
shared_cache=""
control_root=""
managed_only=0
case_list="generic_twin,hubble"
backend_list="system,osmesa"
mode_list="1"
run_timeout=240
settle_timeout_ms=120000
memory_limit_mib=0
lod_memory_percent=""
lod_memory_percent_after_zoom=""
initial_view_only=0
interactive_fps=20
stable_fps=10
canvas_target="./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas"
# The progress faceplate occupies the right edge and lower caption band while
# cache maintenance continues after the current view becomes usable.  Keep the
# original checkpoint as user-visible evidence, but exclude those documented
# overlay bands from geometry scoring.  These are the same physical-pixel
# bounds used by qged_gui_matrix.sh's lighting comparison.
comparison_right_overlay_pixels=60
comparison_top_margin_pixels=24
comparison_bottom_overlay_pixels=174
silhouette_tolerance_pixels=1
# Reapplying a GED center through the float-backed camera may round by a few
# thousandths of a pixel.  Reject a visibly shifted comparison, not harmless
# model-space roundoff whose magnitude scales with the current view.
camera_alignment_tolerance_pixels=0.05
# Safe scenes are expected to reproduce the full-detail image.  For larger
# scenes, shaded SSIM remains a photometric review signal: PoP geometry can
# satisfy its screen-space bound while simplified normals change illumination.
# Topology qualification therefore uses the tolerant silhouette separately.
safe_direct_ssim_minimum=0.999
# Shaded-with-edges submits coincident filled and line primitives.  OSMesa's
# depth tie and subpixel line coverage can move a few edge samples even when
# both paths present the same pixel-target mesh.  Keep a tight, separately
# named raster tolerance for that mode; geometry, exact-presentation, proxy,
# and normalized-error gates below remain unchanged.
safe_direct_edge_ssim_minimum=0.9985
safe_direct_edge_silhouette_maximum=0.00025
pixel_target_ssim_advisory=0.97
tolerant_silhouette_maximum=0.025
normalized_error_epsilon=0.0001
# Public BObol convergence enum values used by the JSON-only verifier.
presentation_constrained_outcome=2
constraint_memory_mask=16
# Image metrics are meaningful only at one physical sampling density.  Fix the
# test process scale instead of inheriting a desktop-specific Qt setting.
quality_device_scale=1
# Named regions use the fixed 1100x800 quality canvas after its documented HUD
# crop.  The side view is the only standard Big Boy camera which cleanly
# separates running gear from the boiler and cab.
bigboy_ae90_front_wheels_crop='150x70+130+320'
bigboy_ae90_running_gear_crop='500x80+130+305'
bigboy_ae90_boiler_cab_crop='430x125+315+220'

usage()
{
    cat <<'EOF'
Usage: qged_lod_quality_matrix.sh [options]

  --build-dir DIR       Current build (default: ./.build)
  --artifact-dir DIR    New results directory (default: /tmp)
  --cache-root DIR      Reuse caches below DIR (default: results directory)
  --shared-cache DIR    Reuse one exact cache directory for every serial run
  --control-root DIR    Reuse full-detail cases from an earlier result set
  --managed-only        Run and record LoD pressure evidence without a control
  --cases LIST          generic_twin,hubble,lucy,multi_lucy,
                        multi_lucy_pair,multi_lucy_quad,multi_lucy_xpush,
                        nist,unique_mesh,unique_mesh_50k_stress,
                        unique_mesh_150k_stress,bigboy
                        (default: generic_twin,hubble)
  --backends LIST       system,osmesa (default: both)
  --modes LIST          comma-separated 0,1,2,4 (default: 1)
  --timeout SECONDS     Per-process timeout (default: 240)
  --settle-ms MSEC      Current-view convergence timeout (default: 120000)
  --memory-limit-mib N Contain a run to an N MiB user-cgroup memory ceiling
  --lod-memory-percent N
                        Set the LoD resident-mesh ceiling to N percent of the
                        available-memory snapshot
  --lod-memory-percent-after-zoom N
                        Apply the LoD resident-mesh ceiling after the first
                        close-view checkpoint, then qualify reclamation and
                        cache-backed restoration
  --interactive-fps N  Minimum interaction FPS target (default: 20)
  --stable-fps N       Minimum stable-view FPS target (default: 10)
  --initial-view-only Stop after the settled ae 90 0 checkpoint

Big Boy uses the partial big_boy.bot mesh conversion.  Its full-detail control
is that same mesh hierarchy, so the case qualifies LoD allocation quality but
not completeness relative to the original BREP model or facetization failures.
It is not run by default.

Managed-only mode is for scenes whose full-detail control would exceed the
host's safe memory or frame-time envelope.  It records allocation, constraint,
proxy, residency, and quality-debt telemetry but does not calculate or imply
matched-image quality.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
	--build-dir) build_dir="$2"; shift 2 ;;
	--artifact-dir) artifact_dir="$2"; shift 2 ;;
	--cache-root) cache_root="$2"; shift 2 ;;
	--shared-cache) shared_cache="$2"; shift 2 ;;
	--control-root) control_root="$2"; shift 2 ;;
	--managed-only) managed_only=1; shift ;;
	--cases) case_list="$2"; shift 2 ;;
	--backends) backend_list="$2"; shift 2 ;;
	--modes) mode_list="$2"; shift 2 ;;
	--timeout) run_timeout="$2"; shift 2 ;;
	--settle-ms) settle_timeout_ms="$2"; shift 2 ;;
	--memory-limit-mib) memory_limit_mib="$2"; shift 2 ;;
	--lod-memory-percent) lod_memory_percent="$2"; shift 2 ;;
	--lod-memory-percent-after-zoom)
	    lod_memory_percent_after_zoom="$2"; shift 2 ;;
	--interactive-fps) interactive_fps="$2"; shift 2 ;;
	--stable-fps) stable_fps="$2"; shift 2 ;;
	--initial-view-only) initial_view_only=1; shift ;;
	--help|-h) usage; exit 0 ;;
	*) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ ! "$run_timeout" =~ ^[1-9][0-9]*$ ]] ||
    [[ ! "$settle_timeout_ms" =~ ^[1-9][0-9]*$ ]] ||
    [[ ! "$memory_limit_mib" =~ ^[0-9]+$ ]]; then
    echo "ERROR: timeouts must be positive integers and the memory limit nonnegative" >&2
    exit 2
fi
if ! awk -v interactive="$interactive_fps" -v stable="$stable_fps" '
    BEGIN {
	exit !(interactive ~ /^[0-9]+([.][0-9]+)?$/ &&
	    stable ~ /^[0-9]+([.][0-9]+)?$/ &&
	    interactive > 0 && stable > 0)
    }
'; then
    echo "ERROR: FPS targets must be positive numbers" >&2
    exit 2
fi
for lod_percent in "$lod_memory_percent" "$lod_memory_percent_after_zoom"; do
    if [[ -n "$lod_percent" ]] && ! awk -v percent="$lod_percent" '
	BEGIN {
	    exit !(percent ~ /^[0-9]+([.][0-9]+)?$/ &&
		percent > 0 && percent <= 50)
	}
    '; then
	echo "ERROR: LoD memory percentage must be greater than 0 and at most 50" >&2
	exit 2
    fi
done
if [[ -n "$lod_memory_percent" &&
      -n "$lod_memory_percent_after_zoom" ]]; then
    echo "ERROR: startup and post-zoom LoD memory settings are mutually exclusive" >&2
    exit 2
fi
allow_memory_quality_constraint=0
if [[ -n "$lod_memory_percent" ||
      -n "$lod_memory_percent_after_zoom" ]]; then
    allow_memory_quality_constraint=1
fi
if ((memory_limit_mib > 0)) && ! command -v systemd-run >/dev/null 2>&1; then
    echo "ERROR: systemd-run is required for --memory-limit-mib" >&2
    exit 2
fi

build_dir="$(realpath -m "$build_dir")"
qged="$build_dir/bin/qged"
imgdiff="$build_dir/bin/imgdiff"
if [[ -z "$artifact_dir" ]]; then
    artifact_dir="/tmp/qged-lod-quality-$(date +%Y%m%d-%H%M%S)-$$"
fi
artifact_dir="$(realpath -m "$artifact_dir")"
if [[ -n "$cache_root" ]]; then
    cache_root="$(realpath -m "$cache_root")"
fi
if [[ -n "$shared_cache" ]]; then
    shared_cache="$(realpath -m "$shared_cache")"
    [[ -d "$shared_cache" ]] || {
	echo "ERROR: shared cache directory not found: $shared_cache" >&2
	exit 2
    }
fi
if [[ -n "$control_root" ]]; then
    control_root="$(realpath -m "$control_root")"
fi
if ((managed_only)) && [[ -n "$control_root" ]]; then
    echo "ERROR: --managed-only and --control-root are mutually exclusive" >&2
    exit 2
fi

if [[ ! -x "$qged" || ! -x "$imgdiff" ]]; then
    echo "ERROR: qged and imgdiff are required in $build_dir/bin" >&2
    exit 2
fi
if [[ ! -x /usr/bin/time ]]; then
    echo "ERROR: /usr/bin/time is required" >&2
    exit 2
fi
for tool in jq identify convert compare awk sha256sum stat; do
    command -v "$tool" >/dev/null 2>&1 || {
	echo "ERROR: $tool is required" >&2
	exit 2
    }
done

# A quality run opens its input databases read-only.  Hash each distinct
# database once for artifact provenance, then use cheap inode/size/timestamp
# metadata to prove that each qged invocation did not mutate it.  Rehashing a
# multi-gigabyte scene before and after every backend/policy combination can
# otherwise take longer than the rendering under test, especially when the
# shared test corpus is on removable storage.
declare -A database_content_hashes=()
if [[ -e "$artifact_dir" ]]; then
    echo "ERROR: artifact directory already exists: $artifact_dir" >&2
    exit 2
fi
mkdir -p "$artifact_dir"/{cases,caches,events,metrics}
printf 'status,run,seconds,reason\n' > "$artifact_dir/results.csv"
printf 'case,backend,mode,view,ssim,phash_hamming,changed_pixel_fraction,silhouette_disagreement,silhouette_disagreement_1px,control_foreground_pixels,lod_faces,lod_render_cost,lod_budget,max_normalized_error,performance_limited,memory_limited,prominent_floor_violations,scene_prominent_floor_violations,structural_boxes,terminal_proxies,point_proxy_threshold,last_render_ms,presented_work_exact,constraint_evidence_mask\n' \
    > "$artifact_dir/metrics.csv"
printf 'case,backend,mode,view,feature,crop,ssim,phash_hamming,changed_pixel_fraction,silhouette_disagreement,silhouette_disagreement_1px,control_foreground_pixels\n' \
    > "$artifact_dir/feature_metrics.csv"
printf 'case,backend,mode,view,faces,lines,render_cost,budget,max_normalized_error,max_projected_error_pixels,max_visual_footprint_pixels,performance_limited,memory_limited,prominent_candidates,prominent_payloads,prominent_floor_violations,scene_prominent_floor_violations,visual_importance_debt,scene_visual_importance_debt,structural_boxes,terminal_proxies,point_proxy_threshold,subpixel_occurrences,subpixel_points,last_render_ms,smoothed_render_ms,presented_work_exact,constraint_evidence_mask,resident_mesh_bytes,resident_mesh_limit_bytes,gpu_live_bytes,gpu_capacity_bytes,visible_targets,presented_occurrences,lod_payloads\n' \
    > "$artifact_dir/managed_metrics.csv"
printf 'case,backend,mode,view,source_path,source_name,occurrence_key,active_cut,presentation_cut,resident_cut,requested_cut,projected_diameter_pixels,projected_area_pixels,projected_perimeter_pixels,visual_footprint_pixels,target_error_pixels,projected_error_pixels,normalized_error,prominent,quality_floor_violation,faces,points\n' \
    > "$artifact_dir/managed_payload_metrics.csv"
printf 'case,backend,mode,view,source_path,source_name,occurrence_key,active_cut,presentation_cut,resident_cut,requested_cut,projected_diameter_pixels,visual_footprint_pixels,target_error_pixels,projected_error_pixels,normalized_error,faces,points\n' \
    > "$artifact_dir/managed_outlier_metrics.csv"
printf 'case,backend,mode,view,quality_class,contract_pass,photometric_target,photometric_target_met,manual_review_required,reason\n' \
    > "$artifact_dir/matched_assessment.csv"
printf 'case,backend,mode,view,quality_class,telemetry_contract_pass,manual_review_required,reason\n' \
    > "$artifact_dir/managed_assessment.csv"

active_test_pid=""
cleanup_active_test()
{
    if [[ -n "$active_test_pid" ]] && kill -0 "$active_test_pid" 2>/dev/null; then
	kill -TERM -- "-$active_test_pid" 2>/dev/null ||
	    kill -TERM "$active_test_pid" 2>/dev/null || true
	wait "$active_test_pid" 2>/dev/null || true
    fi
    active_test_pid=""
}
trap cleanup_active_test EXIT

stanford_database()
{
    if [[ -n "${BOBOL_STANFORD_DB:-}" ]]; then
	printf '%s\n' "$BOBOL_STANFORD_DB"
    elif [[ -r "$build_dir/stanford_local.g" ]]; then
	printf '%s\n' "$build_dir/stanford_local.g"
    else
	printf '%s\n' "$build_dir/stanford.g"
    fi
}

case_spec()
{
    case "$1" in
	generic_twin)
	    local db="$build_dir/Generic_Twin.g"
	    [[ -r "$db" ]] || db="$build_dir/share/db/faa/Generic_Twin.g"
	    printf '%s|%s\n' "$db" "${BOBOL_GENERIC_TWIN_OBJECT:-all}"
	    ;;
	hubble)
	    printf '%s|%s\n' \
		"${BOBOL_HUBBLE_DB:-/home/cyapp/models/NASA/Hubble/Hubble_Space_Telescope.g}" \
		"${BOBOL_HUBBLE_OBJECT:-all.g}"
	    ;;
	lucy)
	    printf '%s|%s\n' "${BOBOL_LUCY_DB:-$build_dir/lucy.g}" all
	    ;;
	multi_lucy)
	    local db
	    db="$(stanford_database)"
	    printf '%s|%s\n' \
		"$db" multi_lucy
	    ;;
	multi_lucy_pair)
	    local db
	    db="$(stanford_database)"
	    printf '%s|%s\n' "$db" "multi_lucy_01 multi_lucy_02"
	    ;;
	multi_lucy_quad)
	    local db
	    db="$(stanford_database)"
	    printf '%s|%s\n' "$db" \
		"multi_lucy_01 multi_lucy_02 multi_lucy_03 multi_lucy_04"
	    ;;
	multi_lucy_xpush)
	    local db
	    db="$(stanford_database)"
	    printf '%s|%s\n' \
		"$db" \
		multi_lucy_xpush
	    ;;
	nist)
	    local db="${BOBOL_NIST_DB:-}"
	    if [[ -z "$db" ]]; then
		local preferred_nist="$build_dir/share/db/nist/NIST_MBE_PMI_7-10.g"
		if [[ -r "$preferred_nist" ]]; then
		    db="$preferred_nist"
		else
		    db="$(find "$build_dir/share/db/nist" -type f \
			-name 'NIST_MBE_PMI_*.g' -print 2>/dev/null |
			LC_ALL=C sort | awk 'NR == 1 {candidate = $0} END {print candidate}')"
		fi
	    fi
	    [[ -n "$db" ]] || return 1
	    local object="${BOBOL_NIST_OBJECT:-}"
	    if [[ -z "$object" ]]; then
		if [[ "$(basename "$db")" == "NIST_MBE_PMI_7-10.g" ]]; then
		    object="NIST_MBE_PMI_7-10.3dm"
		else
		    object="Document"
		fi
	    fi
	    printf '%s|%s\n' "$db" "$object"
	    ;;
	unique_mesh)
	    printf '%s|%s\n' \
		"${BOBOL_UNIQUE_MESH_DB:-$build_dir/unique_mesh_stress.g}" \
		unique_mesh_stress
	    ;;
	unique_mesh_50k_stress)
	    printf '%s|%s\n' \
		"${BOBOL_UNIQUE_MESH_50K_STRESS_DB:-$build_dir/unique_mesh_50k_stress.g}" \
		unique_mesh_stress
	    ;;
	unique_mesh_150k_stress)
	    printf '%s|%s\n' \
		"${BOBOL_UNIQUE_MESH_150K_STRESS_DB:-$build_dir/unique_mesh_150k_stress.g}" \
		unique_mesh_stress
	    ;;
	bigboy)
	    printf '%s|%s\n' "${BOBOL_BIGBOY_DB:-$build_dir/bigboy.g}" \
		"${BOBOL_BIGBOY_OBJECT:-big_boy.bot}"
	    ;;
	*) return 1 ;;
    esac
}

write_events()
{
    local file="$1"
    local image_dir="$2"
    local object="$3"
    local mode="$4"
    local policy="$5"
    local ae90_center="${6:-}"
    local ae90_size="${7:-}"
    local ae35_center="${8:-}"
    local ae35_size="${9:-}"
    local lod_value=1
    [[ "$policy" == full ]] && lod_value=0
    local lod_memory_command=""
    if [[ -n "$lod_memory_percent" ]]; then
	lod_memory_command="\"view lod memory ${lod_memory_percent}\","
    fi
    local post_zoom_memory_events=""
    local post_return_idle_event=""
    local return_checkpoint_suffix=""
    local restore_zoom_events=""
    if [[ -n "$lod_memory_percent_after_zoom" ]]; then
	read -r -d '' post_zoom_memory_events <<EOF || true
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "view lod memory ${lod_memory_percent_after_zoom}"}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_timeout_ms}, "quiet_ms": 150}},
EOF
	read -r -d '' post_return_idle_event <<EOF || true
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_timeout_ms}, "quiet_ms": 150}},
EOF
	return_checkpoint_suffix=","
	read -r -d '' restore_zoom_events <<EOF || true
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "zoom 2"}},
    {"target": ".", "action": "wait_progressive_view_ready",
     "arguments": {"timeout_ms": ${settle_timeout_ms}, "quiet_ms": 150}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_timeout_ms}, "quiet_ms": 150}},
    {"target": "$canvas_target", "action": "checkpoint",
     "arguments": {"name": "$image_dir/ae35_zoom2_restore.png"}}
EOF
    fi
    local ae90_framing='"autoview"'
    local ae35_framing='"autoview"'
    if [[ -n "$ae90_center" && -n "$ae90_size" &&
	  -n "$ae35_center" && -n "$ae35_size" ]]; then
	ae90_framing="\"center $ae90_center\", \"size $ae90_size\""
	ae35_framing="\"center $ae35_center\", \"size $ae35_size\""
    fi

    cat > "$file" <<EOF
{
  "schema": "brlcad.qtcad.events",
  "version": 1,
  "events": [
    {"target": ".", "action": "resize",
     "arguments": {"width": 1100, "height": 800, "stable_ms": 250,
                   "timeout_ms": 5000}},
    {"target": ".", "action": "wait_canvas_ready",
     "arguments": {"timeout_ms": 10000}},
    {"target": ".", "action": "qged_command_batch",
     "arguments": {"commands": [
       "units mm", "view clear", "view lod ${lod_value}",
       "view lod fps ${interactive_fps} ${stable_fps}",
       ${lod_memory_command}
       "view lighting profile studio", "view faceplate center_dot 0",
       "view faceplate adc 0", "view faceplate grid 0",
       "view faceplate irect draw 0", "view faceplate navgizmo 0",
       "view faceplate model_axes 0", "view faceplate scale 0",
       "view faceplate view_axes 0", "view faceplate params 0"]}},
    {"target": "$canvas_target", "action": "checkpoint",
     "arguments": {"name": "$image_dir/empty.png"}},
    {"target": ".", "action": "qged_command_batch",
     "arguments": {"commands": ["draw -m${mode} ${object}",
                                  "ae 90 0", ${ae90_framing}]}},
    {"target": ".", "action": "wait_progressive_view_ready",
     "arguments": {"timeout_ms": ${settle_timeout_ms}, "quiet_ms": 150}},
    {"target": ".", "action": "wait_progressive_idle",
     "arguments": {"timeout_ms": ${settle_timeout_ms}, "quiet_ms": 150}},
    {"target": "$canvas_target", "action": "checkpoint",
     "arguments": {"name": "$image_dir/ae90.png"}},
    {"target": ".", "action": "qged_command_batch",
     "arguments": {"commands": ["ae 35 25", ${ae35_framing}]}},
    {"target": ".", "action": "wait_progressive_view_ready",
     "arguments": {"timeout_ms": ${settle_timeout_ms}, "quiet_ms": 150}},
    {"target": "$canvas_target", "action": "checkpoint",
     "arguments": {"name": "$image_dir/ae35.png"}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "zoom 2"}},
    {"target": ".", "action": "wait_progressive_view_ready",
     "arguments": {"timeout_ms": ${settle_timeout_ms}, "quiet_ms": 150}},
    {"target": "$canvas_target", "action": "checkpoint",
     "arguments": {"name": "$image_dir/ae35_zoom2.png"}},
${post_zoom_memory_events}
    {"target": ".", "action": "wait", "arguments": {"ms": 850}},
    {"target": ".", "action": "qged_command",
     "arguments": {"command": "zoom 0.5"}},
    {"target": ".", "action": "wait_progressive_view_ready",
     "arguments": {"timeout_ms": ${settle_timeout_ms}, "quiet_ms": 150}},
${post_return_idle_event}
    {"target": "$canvas_target", "action": "checkpoint",
     "arguments": {"name": "$image_dir/ae35_return.png"}}${return_checkpoint_suffix}
${restore_zoom_events}
  ]
}
EOF
    if ((initial_view_only)); then
	local reduced_file="${file}.initial"
	jq '.events = .events[0:8]' "$file" > "$reduced_file" || return 1
	mv "$reduced_file" "$file"
    fi
}

run_policy()
{
    local case_name="$1"
    local db="$2"
    local object="$3"
    local backend="$4"
    local mode="$5"
    local policy="$6"
    local cache_dir="$7"
    local ae90_center="${8:-}"
    local ae90_size="${9:-}"
    local ae35_center="${10:-}"
    local ae35_size="${11:-}"
    local run="${case_name}-${backend}-m${mode}-${policy}"
    local out="$artifact_dir/cases/$run"
    local events="$artifact_dir/events/$run.json"
    mkdir -p "$out/images"
    local database_path
    database_path="$(realpath -m "$db")"
    local database_identity_before
    database_identity_before="$(stat -Lc '%d:%i:%s:%Y:%Z' -- "$db")" ||
	return 1
    if [[ -z "${database_content_hashes[$database_path]+set}" ]]; then
	local database_hash
	database_hash="$(sha256sum -- "$db")" || return 1
	database_content_hashes[$database_path]="${database_hash%% *}"
    fi
    write_events "$events" "$out/images" "$object" "$mode" "$policy" \
	"$ae90_center" "$ae90_size" "$ae35_center" "$ae35_size"

    local -a qged_args=(--test-script "$events" --test-report "$out/report.json")
    [[ "$backend" == osmesa ]] && qged_args=(-s "${qged_args[@]}")
    qged_args+=("$db")
    local started=$SECONDS
    printf 'RUN %s\n' "$run"
    local -a limit_command=()
    if ((memory_limit_mib > 0)); then
	limit_command=(systemd-run --user --scope --quiet
	    -p "MemoryMax=${memory_limit_mib}M"
	    -p "MemorySwapMax=0")
    fi
    setsid timeout --signal=TERM --kill-after=5 "$run_timeout" \
	/usr/bin/time -v -o "$out/resource-usage.txt" env \
	QT_SCALE_FACTOR="$quality_device_scale" \
	BU_DIR_CACHE="$cache_dir" \
	QGED_TEST_DEEP_LOD_REPORT="${QGED_TEST_DEEP_LOD_REPORT:-0}" \
	"${limit_command[@]}" "$qged" "${qged_args[@]}" \
	>"$out/stdout.log" 2>"$out/stderr.log" &
    active_test_pid=$!
    local status=0
    wait "$active_test_pid" || status=$?
    active_test_pid=""
    local database_identity_after
    database_identity_after="$(stat -Lc '%d:%i:%s:%Y:%Z' -- "$db")" ||
	status=1
    if [[ "$database_identity_before" != "$database_identity_after" ]]; then
	# BRL-CAD may flush an unchanged writable database header at close,
	# advancing mtime/ctime without changing its bytes.  Rehash only on that
	# exceptional metadata edge; an actual content change still fails the
	# supposedly read-only quality run.
	local database_hash_after
	database_hash_after="$(sha256sum -- "$db")" || status=1
	database_hash_after="${database_hash_after%% *}"
	if [[ "${database_content_hashes[$database_path]}" != "$database_hash_after" ]]; then
	    printf 'FAIL,%s,%s,database-content-changed\n' "$run" \
		"$((SECONDS - started))" >> "$artifact_dir/results.csv"
	    return 1
	fi
    fi
    printf '%s  %s\n' "${database_content_hashes[$database_path]}" \
	"$database_path" \
	> "$out/database.sha256"
    if ((status != 0)); then
	printf 'FAIL,%s,%s,status=%s\n' "$run" "$((SECONDS - started))" \
	    "$status" >> "$artifact_dir/results.csv"
	return 1
    fi
    # Image agreement is not a substitute for a live control plane.  Refine
    # every sampled production state against the finite ownership contract so
    # a plausible framebuffer cannot hide an ownerless 99-percent stall or a
    # non-progress A/B/A control cycle.
    if ! jq -e -f "$script_dir/lod_control_trace.jq" "$out/report.json" \
	    >"$out/control-trace-validation.log" 2>&1; then
	printf 'FAIL,%s,%s,control-trace-contract\n' "$run" \
	    "$((SECONDS - started))" >> "$artifact_dir/results.csv"
	return 1
    fi
    local restore_zoom=0
    [[ -n "$lod_memory_percent_after_zoom" ]] && restore_zoom=1
    if ! jq -e --arg policy "$policy" --arg case_name "$case_name" \
	--arg mode "$mode" --argjson initial_only "$initial_view_only" \
	--argjson restore_zoom "$restore_zoom" \
	--argjson allow_memory_quality_constraint \
	"$allow_memory_quality_constraint" \
	--argjson constrained_outcome "$presentation_constrained_outcome" \
	--argjson memory_mask "$constraint_memory_mask" '
	.success == true and
	([.samples[] | select(.action == "checkpoint")] | length) ==
	    (if $initial_only == 1 then 2 else 5 + $restore_zoom end) and
	all(.samples[] | select(.action == "checkpoint");
	    .canvas_width > 0 and .canvas_height > 0 and
	    .visible_structural_fallback_boxes == 0) and
	($policy != "lod" or
	    (([.samples[] | select(.action == "checkpoint" and
		(((.checkpoint? // "") | endswith("/empty.png")) | not))] |
	      length) ==
		(if $initial_only == 1 then 1 else 4 + $restore_zoom end) and
	     all(.samples[] | select(.action == "checkpoint" and
		(((.checkpoint? // "") | endswith("/empty.png")) | not));
		.lod_convergence_terminal == true and
		.lod_convergence_view_ready == true and
		.lod_control_violation_mask == 0 and
		(if $mode == "0" then
		   ((.presented_cad_lines // 0) > 0)
		 else
		   ((.presented_cad_faces // 0) > 0)
		 end) and
		((.lod_prominent_cad_quality_floor_violations == 0 and
		  .lod_scene_prominent_quality_floor_violations == 0) or
		 ($allow_memory_quality_constraint == 1 and
		  .lod_convergence_outcome == $constrained_outcome and
		  .lod_convergence_memory_limited == true and
		  (((.lod_convergence_constraint_evidence_mask // 0) %
		    ($memory_mask * 2)) >= $memory_mask))) and
		($case_name != "generic_twin" or
		 .lod_convergence_terminal_proxy_occurrences == 0))))
	' "$out/report.json" >"$out/validation.log" 2>&1; then
	printf 'FAIL,%s,%s,report-contract\n' "$run" \
	    "$((SECONDS - started))" >> "$artifact_dir/results.csv"
	return 1
    fi
    printf 'PASS,%s,%s,\n' "$run" "$((SECONDS - started))" \
	>> "$artifact_dir/results.csv"
}

checkpoint_sample()
{
    local report="$1"
    local name="$2"
    jq -c --arg suffix "/$name.png" '
	[.samples[] | select((.checkpoint? // "") | endswith($suffix))][0]
    ' "$report"
}

validate_pair_geometry()
{
    local lod_sample
    local full_sample
    lod_sample="$(checkpoint_sample "$1" "$3")" || return 1
    full_sample="$(checkpoint_sample "$2" "$3")" || return 1
    jq -e -n --argjson lod "$lod_sample" --argjson full "$full_sample" \
	--argjson camera_pixel_tolerance \
	"$camera_alignment_tolerance_pixels" '
	def close($a; $b; $epsilon): (($a - $b) | fabs) <= $epsilon;
	def camera_target_epsilon:
	    1e-6 + (($full.model_view_size | fabs) /
		($full.canvas_height | if . > 0 then . else 1 end) *
		$camera_pixel_tolerance);
	($lod.canvas_width == $full.canvas_width) and
	($lod.canvas_height == $full.canvas_height) and
	close(($lod.canvas_device_pixel_ratio // 1);
	      ($full.canvas_device_pixel_ratio // 1); 1e-6) and
	close($lod.model_view_size; $full.model_view_size;
	      1e-4 + (($full.model_view_size | fabs) * 1e-6)) and
	close($lod.camera_target_x; $full.camera_target_x;
	      camera_target_epsilon) and
	close($lod.camera_target_y; $full.camera_target_y;
	      camera_target_epsilon) and
	close($lod.camera_target_z; $full.camera_target_z;
	      camera_target_epsilon) and
	close($lod.camera_orientation_angle; $full.camera_orientation_angle; 1e-6) and
	close($lod.camera_orientation_axis_x; $full.camera_orientation_axis_x; 1e-6) and
	close($lod.camera_orientation_axis_y; $full.camera_orientation_axis_y; 1e-6) and
	close($lod.camera_orientation_axis_z; $full.camera_orientation_axis_z; 1e-6)
    ' >/dev/null
}

validate_pair_database()
{
    local lod_hash="$1/database.sha256"
    local full_hash="$2/database.sha256"
    [[ -r "$lod_hash" && -r "$full_hash" ]] || return 1
    local lod_digest
    local full_digest
    lod_digest="$(awk '{print $1}' "$lod_hash")" || return 1
    full_digest="$(awk '{print $1}' "$full_hash")" || return 1
    [[ "$lod_digest" == "$full_digest" ]]
}

mask_foreground()
{
    convert "$1" "$2" -compose difference -composite \
	-colorspace Gray -threshold 2% "$3"
}

crop_scene_image()
{
    local input="$1"
    local output="$2"
    local width
    local height
    read -r width height <<< "$(identify -format '%w %h' "$input")"
    local crop_width=$((width - comparison_right_overlay_pixels))
    local crop_height=$((height - comparison_top_margin_pixels -
	comparison_bottom_overlay_pixels))
    if ((crop_width <= 0 || crop_height <= 0)); then
	echo "ERROR: image is too small for the qged scene comparison crop: $input" >&2
	return 1
    fi
    convert "$input" -crop \
	"${crop_width}x${crop_height}+0+${comparison_top_margin_pixels}" \
	+repage "$output"
}

mask_pixel_count()
{
    convert "$1" -format '%[fx:mean*w*h]' info: |
	awk '{printf "%.0f\n", $1}'
}

metric_value()
{
    awk -v label="$1" '
	index($0, label ": ") == 1 {print substr($0, length(label) + 3)}'
}

measure_pair_metrics()
{
    local lod_image="$1"
    local full_image="$2"
    local lod_empty="$3"
    local full_empty="$4"
    local metric_dir="$5"
    local label="$6"
    mkdir -p "$metric_dir"

    local approximate
    approximate="$("$imgdiff" -A "$lod_image" "$full_image" 2>&1)" ||
	return 1
    local ssim
    local hamming
    ssim="$(metric_value SSIM <<< "$approximate")"
    hamming="$(metric_value 'Hamming distance' <<< "$approximate")"
    [[ "$ssim" =~ ^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]] || return 1
    [[ "$hamming" =~ ^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$ ]] || return 1

    local width
    local height
    read -r width height <<< "$(identify -format '%w %h' "$full_image")"
    local changed
    changed="$(compare -metric AE -fuzz 2% "$lod_image" "$full_image" \
	null: 2>&1 || true)"
    changed="${changed%% *}"
    [[ "$changed" =~ ^[0-9]+$ ]] || return 1
    local changed_fraction
    changed_fraction="$(awk -v changed="$changed" -v pixels="$((width * height))" \
	'BEGIN {printf "%.9g", pixels ? changed / pixels : 0}')"

    local lod_mask="$metric_dir/lod-mask.png"
    local full_mask="$metric_dir/full-mask.png"
    local union_mask="$metric_dir/union-mask.png"
    mask_foreground "$lod_image" "$lod_empty" "$lod_mask" || return 1
    mask_foreground "$full_image" "$full_empty" "$full_mask" || return 1
    convert "$lod_mask" "$full_mask" -evaluate-sequence max \
	"$union_mask" || return 1
    local disagreement
    disagreement="$(compare -metric AE "$lod_mask" "$full_mask" \
	null: 2>&1 || true)"
    disagreement="${disagreement%% *}"
    [[ "$disagreement" =~ ^[0-9]+$ ]] || return 1

    local lod_dilated="$metric_dir/lod-mask-dilated-1px.png"
    local full_dilated="$metric_dir/full-mask-dilated-1px.png"
    local lod_unmatched="$metric_dir/lod-mask-unmatched-1px.png"
    local full_unmatched="$metric_dir/full-mask-unmatched-1px.png"
    local tolerant_disagreement_mask="$metric_dir/silhouette-disagreement-1px.png"
    convert "$lod_mask" -morphology Dilate \
	"Diamond:$silhouette_tolerance_pixels" "$lod_dilated" || return 1
    convert "$full_mask" -morphology Dilate \
	"Diamond:$silhouette_tolerance_pixels" "$full_dilated" || return 1
    convert "$lod_mask" "$full_dilated" -fx 'u*(1-v)' \
	"$lod_unmatched" || return 1
    convert "$full_mask" "$lod_dilated" -fx 'u*(1-v)' \
	"$full_unmatched" || return 1
    convert "$lod_unmatched" "$full_unmatched" -evaluate-sequence max \
	"$tolerant_disagreement_mask" || return 1

    local union_pixels
    local lod_pixels
    local control_pixels
    local tolerant_disagreement_pixels
    union_pixels="$(mask_pixel_count "$union_mask")" || return 1
    lod_pixels="$(mask_pixel_count "$lod_mask")" || return 1
    control_pixels="$(mask_pixel_count "$full_mask")" || return 1
    tolerant_disagreement_pixels="$(mask_pixel_count \
	"$tolerant_disagreement_mask")" || return 1
    if ((lod_pixels == 0 || control_pixels == 0)); then
	echo "ERROR: blank scene oracle for $label " \
	    "(lod=$lod_pixels full=$control_pixels)" >&2
	return 1
    fi
    local silhouette_fraction
    local tolerant_silhouette_fraction
    silhouette_fraction="$(awk -v changed="$disagreement" -v pixels="$union_pixels" \
	'BEGIN {printf "%.9g", pixels ? changed / pixels : 0}')"
    tolerant_silhouette_fraction="$(awk \
	-v changed="$tolerant_disagreement_pixels" -v pixels="$union_pixels" \
	'BEGIN {printf "%.9g", pixels ? changed / pixels : 0}')"

    printf '%s,%s,%s,%s,%s,%s' "$ssim" "$hamming" "$changed_fraction" \
	"$silhouette_fraction" "$tolerant_silhouette_fraction" \
	"$control_pixels"
}

named_feature_specs()
{
    local case_name="$1"
    local view="$2"
    if [[ "$case_name" == bigboy && "$view" == ae90 ]]; then
	printf 'front_wheels|%s\n' "$bigboy_ae90_front_wheels_crop"
	printf 'running_gear|%s\n' "$bigboy_ae90_running_gear_crop"
	printf 'boiler_cab|%s\n' "$bigboy_ae90_boiler_cab_crop"
    fi
}

compare_named_features()
{
    local case_name="$1"
    local backend="$2"
    local mode="$3"
    local view="$4"
    local lod_scene="$5"
    local full_scene="$6"
    local lod_empty="$7"
    local full_empty="$8"
    local feature
    local geometry
    while IFS='|' read -r feature geometry; do
	[[ -n "$feature" && -n "$geometry" ]] || continue
	local feature_dir="$artifact_dir/metrics/${case_name}-${backend}-m${mode}-${view}/features/$feature"
	mkdir -p "$feature_dir"
	local lod_feature="$feature_dir/lod-scene.png"
	local full_feature="$feature_dir/full-scene.png"
	local lod_feature_empty="$feature_dir/lod-empty-scene.png"
	local full_feature_empty="$feature_dir/full-empty-scene.png"
	convert "$lod_scene" -crop "$geometry" +repage "$lod_feature" || return 1
	convert "$full_scene" -crop "$geometry" +repage "$full_feature" || return 1
	convert "$lod_empty" -crop "$geometry" +repage \
	    "$lod_feature_empty" || return 1
	convert "$full_empty" -crop "$geometry" +repage \
	    "$full_feature_empty" || return 1
	local metrics
	metrics="$(measure_pair_metrics "$lod_feature" "$full_feature" \
	    "$lod_feature_empty" "$full_feature_empty" "$feature_dir" \
	    "$case_name/$backend/m$mode/$view/$feature")" || return 1
	printf '%s,%s,%s,%s,%s,%s,%s\n' "$case_name" "$backend" "$mode" \
	    "$view" "$feature" "$geometry" "$metrics" \
	    >> "$artifact_dir/feature_metrics.csv"
    done < <(named_feature_specs "$case_name" "$view")
}

compare_pair()
{
    local case_name="$1"
    local backend="$2"
    local mode="$3"
    local view="$4"
    local lod_out="$artifact_dir/cases/${case_name}-${backend}-m${mode}-lod"
    local full_out="${control_root:-$artifact_dir}/cases/${case_name}-${backend}-m${mode}-full"
    local metric_dir="$artifact_dir/metrics/${case_name}-${backend}-m${mode}-${view}"
    local lod_image="$lod_out/images/$view.png"
    local full_image="$full_out/images/$view.png"
    mkdir -p "$metric_dir"

    local lod_image_geometry
    local full_image_geometry
    lod_image_geometry="$(identify -format '%wx%h' "$lod_image")" || return 1
    full_image_geometry="$(identify -format '%wx%h' "$full_image")" || return 1
    if [[ "$lod_image_geometry" != "$full_image_geometry" ]]; then
	echo "ERROR: mismatched physical image geometry for " \
	    "$case_name/$backend/m$mode/$view " \
	    "(lod=$lod_image_geometry full=$full_image_geometry)" >&2
	return 1
    fi

    if ! validate_pair_database "$lod_out" "$full_out"; then
	echo "ERROR: mismatched or unrecorded database content for " \
	    "$case_name/$backend/m$mode/$view" >&2
	return 1
    fi
    if ! validate_pair_geometry "$lod_out/report.json" "$full_out/report.json" \
	"$view"; then
	echo "ERROR: mismatched camera or canvas for $case_name/$backend/m$mode/$view" >&2
	return 1
    fi
    local lod_comparison="$metric_dir/lod-scene.png"
    local full_comparison="$metric_dir/full-scene.png"
    local lod_empty="$metric_dir/lod-empty-scene.png"
    local full_empty="$metric_dir/full-empty-scene.png"
    crop_scene_image "$lod_image" "$lod_comparison" || return 1
    crop_scene_image "$full_image" "$full_comparison" || return 1
    crop_scene_image "$lod_out/images/empty.png" "$lod_empty" || return 1
    crop_scene_image "$full_out/images/empty.png" "$full_empty" || return 1

    local image_metrics
    image_metrics="$(measure_pair_metrics "$lod_comparison" "$full_comparison" \
	"$lod_empty" "$full_empty" "$metric_dir" \
	"$case_name/$backend/m$mode/$view")" || return 1
    compare_named_features "$case_name" "$backend" "$mode" "$view" \
	"$lod_comparison" "$full_comparison" "$lod_empty" "$full_empty" ||
	return 1

    local lod_sample
    lod_sample="$(checkpoint_sample "$lod_out/report.json" "$view")" || return 1
    local telemetry
    telemetry="$(jq -r '[
	(.active_lod_scene_faces // 0),
	(.active_lod_scene_render_cost // 0),
	(.lod_scene_render_cost_budget // 0),
	(.lod_scene_maximum_normalized_visual_error //
	 .lod_max_cad_normalized_error // 0),
	(.lod_convergence_performance_limited // false),
	(.lod_convergence_memory_limited // false),
	(.lod_prominent_cad_quality_floor_violations // 0),
	(.lod_scene_prominent_quality_floor_violations // 0),
	(.visible_structural_fallback_boxes // 0),
	(.lod_convergence_terminal_proxy_occurrences // 0),
	(.cad_point_proxy_pixel_threshold_max // 1),
	(.last_render_ms // 0),
	(.presented_cad_work_exact // false),
	(.lod_convergence_constraint_evidence_mask // 0)] | @csv' \
	<<< "$lod_sample")" ||
	return 1

    printf '%s,%s,%s,%s,%s,%s\n' "$case_name" "$backend" "$mode" \
	"$view" "$image_metrics" "$telemetry" >> "$artifact_dir/metrics.csv"
}

record_managed_metrics()
{
    local case_name="$1"
    local backend="$2"
    local mode="$3"
    local view="$4"
    local report="$artifact_dir/cases/${case_name}-${backend}-m${mode}-lod/report.json"
    local sample
    sample="$(checkpoint_sample "$report" "$view")" || return 1
    local telemetry
    telemetry="$(jq -r '[
	(.active_lod_scene_faces // 0),
	(.presented_cad_lines // 0),
	(.active_lod_scene_render_cost // 0),
	(.lod_scene_render_cost_budget // 0),
	(.lod_scene_maximum_normalized_visual_error //
	 .lod_max_cad_normalized_error // 0),
	(.lod_max_cad_projected_error_pixels // 0),
	(.lod_max_cad_visual_footprint_pixels // 0),
	(.lod_convergence_performance_limited // false),
	(.lod_convergence_memory_limited // false),
	(.lod_scene_prominent_candidates // 0),
	(.lod_prominent_cad_payloads // 0),
	(.lod_prominent_cad_quality_floor_violations // 0),
	(.lod_scene_prominent_quality_floor_violations // 0),
	(.lod_cad_visual_importance_debt // 0),
	(.lod_scene_visual_importance_debt // 0),
	(.visible_structural_fallback_boxes // 0),
	(.lod_convergence_terminal_proxy_occurrences // 0),
	(.cad_point_proxy_pixel_threshold_max // 1),
	(.lod_convergence_presented_subpixel_occurrences // 0),
	(.active_cad_subpixel_proxy_points // 0),
	(.last_render_ms // 0),
	(.smoothed_render_ms // 0),
	(.presented_cad_work_exact // false),
	(.lod_convergence_constraint_evidence_mask // 0),
	(.lod_convergence_resident_mesh_bytes // 0),
	(.lod_convergence_resident_mesh_limit_bytes // 0),
	(.lod_gpu_atlas_live_bytes // 0),
	(.lod_gpu_atlas_configured_capacity_bytes // 0),
	(.lod_convergence_visible_targets // 0),
	(.presented_cad_occurrences // 0),
	(.active_lod_cad_payloads // 0)] | @csv' <<< "$sample")" ||
	return 1
    printf '%s,%s,%s,%s,%s\n' "$case_name" "$backend" "$mode" \
	"$view" "$telemetry" >> "$artifact_dir/managed_metrics.csv"
    jq -r --arg case_name "$case_name" --arg backend "$backend" \
	--arg mode "$mode" --arg view "$view" '
	(.active_progressive_cad_payload_samples // [])[] |
	[$case_name, $backend, $mode, $view,
	 (.source_path // ""), (.source_name // ""),
	 (.occurrence_key // ""), (.active_cut // -1),
	 (.presentation_cut // .active_cut // -1), (.resident_cut // -1),
	 (.requested_cut // -1), (.projected_diameter_pixels // 0),
	 (.projected_area_pixels // 0), (.projected_perimeter_pixels // 0),
	 (.visual_footprint_pixels // 0), (.target_error_pixels // 0),
	 (.projected_error_pixels // 0), (.normalized_error // 0),
	 (.prominent // false), (.quality_floor_violation // false),
	 (.faces // 0), (.points // 0)] | @csv
	' <<< "$sample" >> "$artifact_dir/managed_payload_metrics.csv" ||
	return 1
    jq -r --arg case_name "$case_name" --arg backend "$backend" \
	--arg mode "$mode" --arg view "$view" '
	(.lod_cad_visual_importance_outliers // [])[] |
	[$case_name, $backend, $mode, $view,
	 (.source_path // ""), (.source_name // ""),
	 (.occurrence_key // ""), (.active_cut // -1),
	 (.presentation_cut // -1), (.resident_cut // -1),
	 (.requested_cut // -1), (.projected_diameter_pixels // 0),
	 (.visual_footprint_pixels // 0), (.target_error_pixels // 0),
	 (.projected_error_pixels // 0), (.normalized_error // 0),
	 (.faces // 0), (.points // 0)] | @csv
    ' <<< "$sample" >> "$artifact_dir/managed_outlier_metrics.csv" ||
	return 1
}

assess_matched_metrics()
{
    awk -F, -v OFS=, \
	-v safe_ssim="$safe_direct_ssim_minimum" \
	-v safe_edge_ssim="$safe_direct_edge_ssim_minimum" \
	-v safe_edge_silhouette="$safe_direct_edge_silhouette_maximum" \
	-v pixel_ssim="$pixel_target_ssim_advisory" \
	-v silhouette_max="$tolerant_silhouette_maximum" \
	-v error_epsilon="$normalized_error_epsilon" \
	-v allow_memory_floor="$allow_memory_quality_constraint" \
	-v memory_mask="$constraint_memory_mask" '
	NR == 1 {next}
	{
	    safe = ($1 == "generic_twin")
	    constrained = ($15 == "true" || $16 == "true" || ($24 + 0) > 0)
	    quality_class = safe ? "safe-direct" :
		(constrained ? "constrained" : "pixel-target")
	    edge_overlay = ($3 == "4")
	    target = safe ? (edge_overlay ? safe_edge_ssim : safe_ssim) :
		pixel_ssim
	    photometric = (($5 + 0) >= target)
	    silhouette_limit = safe ?
		(edge_overlay ? safe_edge_silhouette : 0) : silhouette_max
	    memory_floor = (allow_memory_floor && $16 == "true" &&
		(($24 + 0) % (memory_mask * 2)) >= memory_mask)
	    floor = (($17 + 0) == 0 && ($18 + 0) == 0) || memory_floor
	    topology = (($9 + 0) <= silhouette_limit && ($19 + 0) == 0 &&
		($20 + 0) == 0 && floor && $23 == "true")
	    error_bound = constrained || (($14 + 0) <= 1.0 + error_epsilon)
	    constraint = constrained ? (($24 + 0) > 0) : (($24 + 0) == 0)
	    contract = topology && error_bound && constraint &&
		(!safe || photometric)
	    review = constrained || !photometric
	    reason = "ok"
	    if (!topology)
		reason = "topology-or-presentation"
	    else if (!error_bound)
		reason = "pixel-error"
	    else if (!constraint)
		reason = "constraint-evidence"
	    else if (!photometric)
		reason = "photometric-advisory"
	    print $1, $2, $3, $4, quality_class,
		(contract ? "true" : "false"), target,
		(photometric ? "true" : "false"),
		(review ? "true" : "false"), reason
	    if (!contract)
		failed = 1
	}
	END {exit failed ? 1 : 0}
    ' "$artifact_dir/metrics.csv" >> "$artifact_dir/matched_assessment.csv"
}

assess_managed_metrics()
{
    awk -F, -v OFS=, -v error_epsilon="$normalized_error_epsilon" \
	-v allow_memory_floor="$allow_memory_quality_constraint" \
	-v memory_mask="$constraint_memory_mask" '
	NR == 1 {next}
	{
	    constrained = ($12 == "true" || $13 == "true" || ($28 + 0) > 0)
	    quality_class = constrained ? "constrained" : "pixel-target"
	    memory_floor = (allow_memory_floor && $13 == "true" &&
		(($28 + 0) % (memory_mask * 2)) >= memory_mask)
	    floor = (($16 + 0) == 0 && ($17 + 0) == 0) || memory_floor
	    presentation = (($20 + 0) == 0 && ($21 + 0) == 0 &&
		floor && $27 == "true")
	    error_bound = constrained || (($9 + 0) <= 1.0 + error_epsilon)
	    constraint = constrained ? (($28 + 0) > 0) : (($28 + 0) == 0)
	    contract = presentation && error_bound && constraint
	    reason = "visual-review-required"
	    if (!presentation)
		reason = "presentation-or-prominent-floor"
	    else if (!error_bound)
		reason = "pixel-error"
	    else if (!constraint)
		reason = "constraint-evidence"
	    print $1, $2, $3, $4, quality_class,
		(contract ? "true" : "false"), "true", reason
	    if (!contract)
		failed = 1
	}
	END {exit failed ? 1 : 0}
    ' "$artifact_dir/managed_metrics.csv" >> "$artifact_dir/managed_assessment.csv"
}

failures=0
quality_views=(ae90 ae35 ae35_zoom2 ae35_return)
if ((initial_view_only)); then
    quality_views=(ae90)
elif [[ -n "$lod_memory_percent_after_zoom" ]]; then
    quality_views+=(ae35_zoom2_restore)
fi
IFS=',' read -r -a cases <<< "$case_list"
IFS=',' read -r -a backends <<< "$backend_list"
IFS=',' read -r -a modes <<< "$mode_list"
for case_name in "${cases[@]}"; do
    spec="$(case_spec "$case_name")" || {
	echo "ERROR: unknown or unavailable case: $case_name" >&2
	exit 2
    }
    IFS='|' read -r db object <<< "$spec"
    [[ -r "$db" ]] || {
	echo "ERROR: database not found: $db" >&2
	exit 2
    }
    for backend in "${backends[@]}"; do
	[[ "$backend" == system || "$backend" == osmesa ]] || {
	    echo "ERROR: unknown backend: $backend" >&2
	    exit 2
	}
	if [[ -n "$shared_cache" ]]; then
	    cache_dir="$shared_cache"
	else
	    cache_dir="${cache_root:-$artifact_dir/caches}/${case_name}-${backend}"
	    mkdir -p "$cache_dir"
	fi
	for mode in "${modes[@]}"; do
	    [[ "$mode" == 0 || "$mode" == 1 || "$mode" == 2 ||
	       "$mode" == 4 ]] || {
		echo "ERROR: quality comparison supports modes 0, 1, 2, and 4" >&2
		exit 2
	    }
	    pair_valid=1
	    lod_valid=1
	    ae90_center=""
	    ae90_size=""
	    ae35_center=""
	    ae35_size=""
	    if ((managed_only)); then
		pair_valid=0
	    elif [[ -z "$control_root" ]]; then
		run_policy "$case_name" "$db" "$object" "$backend" "$mode" full \
		    "$cache_dir" || { failures=$((failures + 1)); pair_valid=0; }
	    elif [[ ! -r "$control_root/cases/${case_name}-${backend}-m${mode}-full/report.json" ]]; then
		echo "ERROR: reusable control report not found for $case_name/$backend/m$mode" >&2
		failures=$((failures + 1))
		pair_valid=0
	    fi
	    control_report="${control_root:-$artifact_dir}/cases/${case_name}-${backend}-m${mode}-full/report.json"
	    if ((pair_valid)); then
		IFS=$'\t' read -r ae90_center ae90_size <<< "$(jq -r '
		    [.samples[] |
		     select((.checkpoint? // "") | endswith("/ae90.png"))][0] |
		    if any([.camera_target_x, .camera_target_y,
			.camera_target_z, .model_view_size][]; . == null) then
		      empty
		    else
		      [([.camera_target_x, .camera_target_y, .camera_target_z] |
			map(tostring) | join(" ")), .model_view_size] | @tsv
		    end
		' "$control_report")"
		if ((!initial_view_only)); then
		    IFS=$'\t' read -r ae35_center ae35_size <<< "$(jq -r '
			[.samples[] |
			 select((.checkpoint? // "") | endswith("/ae35.png"))][0] |
			if any([.camera_target_x, .camera_target_y,
			    .camera_target_z, .model_view_size][]; . == null) then
			  empty
			else
			  [([.camera_target_x, .camera_target_y,
			      .camera_target_z] | map(tostring) | join(" ")),
			    .model_view_size] | @tsv
			end
		    ' "$control_report")"
		fi
		missing_framing=0
		if [[ -z "$ae90_center" || -z "$ae90_size" ]]; then
		    missing_framing=1
		elif ((!initial_view_only)) &&
		     [[ -z "$ae35_center" || -z "$ae35_size" ]]; then
		    missing_framing=1
		fi
		if ((missing_framing)); then
		    echo "ERROR: reusable control lacks per-view camera framing" >&2
		    failures=$((failures + 1))
		    pair_valid=0
		fi
	    fi
	    run_policy "$case_name" "$db" "$object" "$backend" "$mode" lod \
		"$cache_dir" "$ae90_center" "$ae90_size" \
		"$ae35_center" "$ae35_size" || {
		    failures=$((failures + 1)); pair_valid=0; lod_valid=0;
		}
	    if ((pair_valid)); then
		for view in "${quality_views[@]}"; do
		    if ! compare_pair "$case_name" "$backend" "$mode" "$view"; then
			printf 'FAIL,%s-%s-m%s-%s,0,comparison-contract\n' \
			    "$case_name" "$backend" "$mode" "$view" \
			    >> "$artifact_dir/results.csv"
			failures=$((failures + 1))
		    fi
		done
	    elif ((managed_only && lod_valid)); then
		for view in "${quality_views[@]}"; do
		    if ! record_managed_metrics "$case_name" "$backend" "$mode" \
			"$view"; then
			printf 'FAIL,%s-%s-m%s-%s,0,managed-metrics\n' \
			    "$case_name" "$backend" "$mode" "$view" \
			    >> "$artifact_dir/results.csv"
			failures=$((failures + 1))
		    fi
		done
	    fi
	done
    done
done

if ! assess_matched_metrics; then
    printf 'FAIL,matched-quality-assessment,0,quality-contract\n' \
	>> "$artifact_dir/results.csv"
    failures=$((failures + 1))
fi
if ! assess_managed_metrics; then
    printf 'FAIL,managed-quality-assessment,0,quality-contract\n' \
	>> "$artifact_dir/results.csv"
    failures=$((failures + 1))
fi

printf 'Artifacts: %s\n' "$artifact_dir"
exit "$failures"
