#!/usr/bin/env bash
#
# Exercise qged's complete external-raytrace framebuffer path together with
# progressive CAD, hierarchy selection, an Obol edit manipulator, faceplate
# drawing, all composition modes, resize, and fbserv client-slot reuse.

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_root="$(cd "$script_dir/../../.." && pwd)"
build_dir="$source_root/.build"
database=""
artifact_dir=""
backend_list="system,osmesa"
run_timeout=240

usage()
{
    cat <<'EOF'
Usage: qged_framebuffer_matrix.sh [options]

  --build-dir DIR       Current build (default: ./.build)
  --database FILE       Generic Twin database override
  --artifact-dir DIR    New results directory (never cleared)
  --backends LIST       system,osmesa (default: both)
  --timeout SECONDS     Per-process timeout (default: 240)
EOF
}

original_arguments=("$@")
while [[ $# -gt 0 ]]; do
    case "$1" in
	--build-dir) build_dir="$2"; shift 2 ;;
	--database) database="$2"; shift 2 ;;
	--artifact-dir) artifact_dir="$2"; shift 2 ;;
	--backends) backend_list="$2"; shift 2 ;;
	--timeout) run_timeout="$2"; shift 2 ;;
	--help|-h) usage; exit 0 ;;
	*) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "${DISPLAY:-}" ]]; then
    if ! command -v xvfb-run >/dev/null 2>&1; then
	echo "ERROR: a display or xvfb-run is required" >&2
	exit 2
    fi
    exec xvfb-run -a "$0" "${original_arguments[@]}"
fi

build_dir="$(realpath -m "$build_dir")"
qged="$build_dir/bin/qged"
if [[ ! -x "$qged" ]]; then
    echo "ERROR: qged is required: $qged" >&2
    exit 2
fi

if [[ -z "$database" ]]; then
    database="$build_dir/Generic_Twin.g"
    if [[ ! -f "$database" || "$(stat -c %s "$database")" -lt 1000 ]]; then
	database="$build_dir/share/db/faa/Generic_Twin.g"
    fi
fi
database="$(realpath -m "$database")"
if [[ ! -f "$database" ]]; then
    echo "ERROR: Generic Twin database is required: $database" >&2
    exit 2
fi

if [[ -z "$artifact_dir" ]]; then
    artifact_dir="$build_dir/qged-framebuffer-matrix/$(date +%Y%m%d-%H%M%S)"
fi
artifact_dir="$(realpath -m "$artifact_dir")"
if [[ -e "$artifact_dir" ]]; then
    echo "ERROR: artifact directory already exists: $artifact_dir" >&2
    exit 2
fi
mkdir -p "$artifact_dir"/{cases,caches,events}
printf 'status,run,seconds,detail\n' > "$artifact_dir/results.csv"

for required_tool in jq identify compare; do
    if ! command -v "$required_tool" >/dev/null 2>&1; then
	echo "ERROR: $required_tool is required" >&2
	exit 2
    fi
done

canvas_target="./i:cad-central/i:cad-quad/i:view-upper-right/i:cad-canvas"

write_events()
{
    local file="$1"
    local image_dir="$2"
    cat > "$file" <<EOF
{
  "schema": "brlcad.qtcad.events",
  "version": 1,
  "events": [
    {"target": ".", "action": "resize", "arguments": {"width": 900, "height": 700}},
    {"target": ".", "action": "qged_command_batch", "arguments": {"commands": ["ae 90 0", "view faceplate center_dot 1", "draw -m1 all"]}},
    {"target": ".", "action": "wait_progressive_scope_ready", "arguments": {"timeout_ms": 30000, "quiet_ms": 50}},
    {"target": ".", "action": "wait_progressive_idle", "arguments": {"timeout_ms": 30000, "quiet_ms": 100}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/cad.png"}},
    {"target": ".", "action": "qged_command_batch", "arguments": {"commands": ["in edit.s ell 0 0 0 150 0 0 0 100 0 0 0 75", "g edit.c edit.s", "draw -m1 edit.c"]}},
    {"target": ".", "action": "wait_progressive_idle", "arguments": {"timeout_ms": 10000, "quiet_ms": 100}},
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "set_expanded", "arguments": {"labels": ["all"], "expanded": true}},
    {"target": "./n:Hierarchy/i:hierarchy-tree", "action": "set_current", "arguments": {"labels": ["all", "0xxx_series"]}},
    {"target": ".", "action": "wait", "arguments": {"ms": 50}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/selected.png"}},
    {"target": "i:org.brlcad.qged.edit.ell.activate", "action": "activate", "arguments": {"checked": true}},
    {"target": "i:primitiveEdit.target", "action": "set_text", "arguments": {"text": "edit.c/edit.s"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "edit edit.c/edit.s status", "contains": "ell"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/edit.png"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "ert -P 1 -H 16"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 5000}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view faceplate fb", "contains": "2"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/underlay.png"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view faceplate fb 1"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/overlay.png"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view faceplate fb 3"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/interlay.png"}},
    {"target": ".", "action": "resize", "arguments": {"width": 1100, "height": 760}},
    {"target": ".", "action": "qged_command_batch", "arguments": {"commands": ["ae 35 25", "autoview"]}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "ert -P 1 -H 16"}},
    {"target": ".", "action": "wait", "arguments": {"ms": 5000}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/resized.png"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view faceplate fb 0"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/off.png"}}
  ]
}
EOF
}

pixel_delta()
{
    local output
    output="$(compare -metric AE "$1" "$2" null: 2>&1 || true)"
    output="${output%% *}"
    [[ "$output" =~ ^[0-9]+$ ]] || return 1
    printf '%s' "$output"
}

validate_case()
{
    local out="$1"
    local report="$out/report.json"
    local images="$out/images"
    local logs=("$out/stdout.log" "$out/stderr.log")
    local name

    jq -e '.success == true' "$report" >/dev/null || return 1
    for name in cad selected edit underlay overlay interlay resized off; do
	[[ -s "$images/$name.png" ]] || return 1
    done

    local initial_size resized_size
    initial_size="$(identify -format '%w %h' "$images/cad.png")"
    resized_size="$(identify -format '%w %h' "$images/resized.png")"
    local initial_width initial_height resized_width resized_height
    read -r initial_width initial_height <<< "$initial_size"
    read -r resized_width resized_height <<< "$resized_size"
    ((initial_width > 0 && initial_height > 0 &&
	resized_width > initial_width && resized_height > initial_height)) ||
	return 1

    jq -e --argjson resized_width "$resized_width" \
	--argjson resized_height "$resized_height" '
      ([.samples[] | select((.checkpoint // "") | endswith("/cad.png"))][0]) as $cad |
      ([.samples[] | select((.checkpoint // "") | endswith("/selected.png"))][0]) as $selected |
      ([.samples[] | select((.checkpoint // "") | endswith("/resized.png"))][0]) as $resized |
      ($cad.lod_convergence_visible_targets // 0) == 709 and
      ($cad.lod_convergence_view_ready // false) == true and
      ((($cad.active_lod_cad_payloads // 0) +
        ($cad.active_cad_subpixel_proxy_points // 0)) >= 709) and
      ($cad.visible_structural_fallback_boxes // 0) == 0 and
      ($selected.selection_count // 0) == 1 and
      ($selected.active_lod_cad_payloads // 0) > 0 and
      ($selected.visible_structural_fallback_boxes // 0) == 0 and
      ($resized.canvas_width // 0) == $resized_width and
      ($resized.canvas_height // 0) == $resized_height and
      ($resized.active_lod_cad_payloads // 0) > 0 and
      ($resized.visible_structural_fallback_boxes // 0) == 0
    ' "$report" >/dev/null || return 1

    local delta
    delta="$(pixel_delta "$images/cad.png" "$images/selected.png")" ||
	return 1
    ((delta >= 100)) || return 1
    delta="$(pixel_delta "$images/selected.png" "$images/edit.png")" ||
	return 1
    ((delta >= 100)) || return 1
    delta="$(pixel_delta "$images/edit.png" "$images/underlay.png")" ||
	return 1
    ((delta >= 1000)) || return 1
    delta="$(pixel_delta "$images/underlay.png" "$images/overlay.png")" ||
	return 1
    ((delta >= 1000)) || return 1
    delta="$(pixel_delta "$images/overlay.png" "$images/interlay.png")" ||
	return 1
    ((delta >= 100)) || return 1
    delta="$(pixel_delta "$images/resized.png" "$images/off.png")" ||
	return 1
    ((delta >= 1000)) || return 1

    if [[ "$(grep -h -c 'Raytrace complete' "${logs[@]}" |
	    awk '{n += $1} END {print n+0}')" -lt 2 ]] ||
	! grep -h -F "size=${initial_width}x${initial_height}" \
	    "${logs[@]}" >/dev/null ||
	! grep -h -F "size=${resized_width}x${resized_height}" \
	    "${logs[@]}" >/dev/null ||
	grep -h -Eq 'terminate called|Aborted|progressive whole-target source failed' \
	    "${logs[@]}"; then
	return 1
    fi
    return 0
}

failures=0
IFS=',' read -r -a backends <<< "$backend_list"
for backend in "${backends[@]}"; do
    if [[ "$backend" != "system" && "$backend" != "osmesa" ]]; then
	echo "ERROR: invalid backend '$backend'" >&2
	exit 2
    fi
    run="generic-twin-$backend-framebuffer"
    out="$artifact_dir/cases/$run"
    cache="$artifact_dir/caches/$run"
    events="$artifact_dir/events/$run.json"
    mkdir -p "$out/images" "$cache"
    cp "$database" "$out/model.g"
    write_events "$events" "$out/images"

    args=(--test-script "$events" --test-report "$out/report.json" "$out/model.g")
    [[ "$backend" == "osmesa" ]] && args=(-s "${args[@]}")
    echo "RUN $run"
    started=$SECONDS
    if env BU_DIR_CACHE="$cache" timeout --signal=TERM "$run_timeout" \
	    "$qged" "${args[@]}" >"$out/stdout.log" 2>"$out/stderr.log" &&
	validate_case "$out"; then
	printf 'PASS,%s,%s,\n' "$run" "$((SECONDS - started))" \
	    >> "$artifact_dir/results.csv"
    else
	printf 'FAIL,%s,%s,workflow-or-validation\n' "$run" \
	    "$((SECONDS - started))" >> "$artifact_dir/results.csv"
	failures=$((failures + 1))
    fi
done

echo "Artifacts: $artifact_dir"
echo "Failures: $failures"
exit "$failures"
