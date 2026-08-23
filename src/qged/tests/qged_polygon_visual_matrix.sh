#!/usr/bin/env bash
#
# Qualify retained polygon drawing over simultaneous shaded and wireframe CAD.
# The same event stream runs through qged's System GL and OSMesa widgets; the
# verifier checks command/store coherence, visible styling and selection
# affordances, resize/zoom behavior, cleanup, and cross-renderer placement.

set -uo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source_root="$(cd "$script_dir/../../.." && pwd)"
build_dir="$source_root/.build"
database=""
artifact_dir=""
backend_list="system,osmesa"
shaded_target="all.g"
wire_target="tor.r"
run_timeout=120

usage()
{
    cat <<'EOF'
Usage: qged_polygon_visual_matrix.sh [options]

  --build-dir DIR       Current build (default: ./.build)
  --database FILE       Small CAD database override (default: installed moss.g)
  --artifact-dir DIR    New results directory (default: a /tmp directory)
  --backends LIST       system,osmesa (default: both)
  --shaded-target PATH  CAD path drawn shaded (default: all.g)
  --wire-target PATH    CAD path drawn wireframe (default: tor.r)
  --timeout SECONDS     Per-process timeout (default: 120)
EOF
}

original_arguments=("$@")
while [[ $# -gt 0 ]]; do
    case "$1" in
	--build-dir) build_dir="$2"; shift 2 ;;
	--database) database="$2"; shift 2 ;;
	--artifact-dir) artifact_dir="$2"; shift 2 ;;
	--backends) backend_list="$2"; shift 2 ;;
	--shaded-target) shaded_target="$2"; shift 2 ;;
	--wire-target) wire_target="$2"; shift 2 ;;
	--timeout) run_timeout="$2"; shift 2 ;;
	--help|-h) usage; exit 0 ;;
	*) echo "ERROR: unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "${DISPLAY:-}" && "${QT_QPA_PLATFORM:-}" != "offscreen" ]]; then
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
    database="$build_dir/share/db/moss.g"
fi
database="$(realpath -m "$database")"
if [[ ! -f "$database" ]]; then
    echo "ERROR: CAD database is required: $database" >&2
    exit 2
fi

if [[ -z "$artifact_dir" ]]; then
    artifact_dir="/tmp/qged-polygon-visual-$(date +%Y%m%d-%H%M%S)-$$"
fi
artifact_dir="$(realpath -m "$artifact_dir")"
if [[ -e "$artifact_dir" ]]; then
    echo "ERROR: artifact directory already exists: $artifact_dir" >&2
    exit 2
fi
mkdir -p "$artifact_dir"/{cases,caches,events,masks}
printf 'status,run,seconds,detail\n' > "$artifact_dir/results.csv"

for required_tool in jq identify compare convert awk; do
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
    {"target": ".", "action": "qged_command_batch", "arguments": {"commands": ["view clear", "ae 35 25"]}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/blank.png"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "draw -m1 $shaded_target"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "draw -m0 $wire_target"}},
    {"target": ".", "action": "wait_progressive_scope_ready", "arguments": {"timeout_ms": 30000, "quiet_ms": 100}},
    {"target": ".", "action": "qged_command_batch", "arguments": {"commands": ["ae 35 25", "autoview"]}},
    {"target": ".", "action": "wait_progressive_idle", "arguments": {"timeout_ms": 30000, "quiet_ms": 100}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "who", "contains": "$shaded_target"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "who", "contains": "$wire_target"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/cad.png"}},

    {"target": ".", "action": "qged_command_batch", "arguments": {"commands": ["view polygon create visual_outer 180 140 rectangle", "view polygon update visual_outer 540 430", "view polygon create visual_hole 300 230 rectangle", "view polygon update visual_hole 420 340", "view polygon csg visual_outer - visual_hole", "view feature delete visual_hole", "view feature style set visual_outer color 255/0/255", "view polygon fill_color visual_outer 0/255/255", "view polygon fill visual_outer 1 1 8"]}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view polygon point_count visual_outer", "numeric_gt": 4}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view polygon selected visual_outer", "contains": "0"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/styled.png"}},

    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view polygon selected visual_outer 1"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view polygon selected visual_outer", "contains": "1"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view polygon select visual_outer 180 140"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/selected.png"}},

    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view polygon move visual_outer 500 120"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/moved.png"}},

    {"target": ".", "action": "resize", "arguments": {"width": 1100, "height": 800}},
    {"target": ".", "action": "wait", "arguments": {"ms": 100}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/resize-early.png"}},
    {"target": ".", "action": "wait_progressive_idle", "arguments": {"timeout_ms": 10000, "quiet_ms": 100}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/resized.png"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "zoom 2"}},
    {"target": ".", "action": "wait_progressive_idle", "arguments": {"timeout_ms": 10000, "quiet_ms": 100}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/zoomed.png"}},

    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view feature delete visual_outer"}},
    {"target": ".", "action": "qged_command_expect", "arguments": {"command": "view feature list visual_outer", "not_contains": "visual_outer"}},
    {"target": "$canvas_target", "action": "checkpoint", "arguments": {"name": "$image_dir/clean.png"}}
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

make_color_mask()
{
    local input="$1"
    local output="$2"
    local color="$3"
    convert "$input" -alpha off -fuzz 8% -fill black +opaque "$color" \
	-fill white -opaque "$color" "$output"
}

mask_count()
{
    convert "$1" -colorspace Gray -format '%[fx:mean*w*h]' info: |
	awk '{printf "%.0f", $1}'
}

validate_case()
{
    local out="$1"
    local report="$out/report.json"
    local images="$out/images"
    local name

    jq -e '.success == true' "$report" >/dev/null || return 1
    for name in blank cad styled selected moved resize-early resized zoomed clean; do
	[[ -s "$images/$name.png" ]] || return 1
    done

    local initial_width initial_height resized_width resized_height
    read -r initial_width initial_height <<< \
	"$(identify -format '%w %h' "$images/cad.png")"
    read -r resized_width resized_height <<< \
	"$(identify -format '%w %h' "$images/resized.png")"
    ((initial_width > 0 && initial_height > 0 &&
	resized_width > initial_width && resized_height > initial_height)) ||
	return 1

    jq -e --argjson width "$resized_width" --argjson height "$resized_height" '
      ([.samples[] | select((.checkpoint // "") | endswith("/cad.png"))][0]) as $cad |
      ([.samples[] | select((.checkpoint // "") | endswith("/resize-early.png"))][0]) as $resize_early |
      ([.samples[] | select((.checkpoint // "") | endswith("/resized.png"))][-1]) as $resized |
      ([.samples[] | select((.checkpoint // "") | endswith("/clean.png"))][0]) as $clean |
      ($cad.tree_top_items // 0) >= 1 and
      ($cad.lod_convergence_fraction // 0) == 1 and
      ($cad.visible_structural_fallback_boxes // 0) == 0 and
      ($cad.presented_cad_occurrences // 0) >= 2 and
      ($resize_early.presented_cad_occurrences // 0) >= 2 and
      ($resized.viewport_width // 0) == $width and
      ($resized.viewport_height // 0) == $height and
      ($clean.tree_top_items // 0) >= 1 and
      any(.samples[]; .command_output == "1")
    ' "$report" >/dev/null || return 1

    local mask_dir="$out/masks"
    mkdir -p "$mask_dir"
    convert "$images/blank.png" -gravity North -crop 100%x85%+0+0 +repage \
	"$mask_dir/blank-upper.png"
    convert "$images/cad.png" -gravity North -crop 100%x85%+0+0 +repage \
	"$mask_dir/cad-upper.png"
    local delta
    delta="$(pixel_delta "$mask_dir/blank-upper.png" "$mask_dir/cad-upper.png")" || return 1
    ((delta >= 500)) || return 1
    delta="$(pixel_delta "$images/cad.png" "$images/styled.png")" || return 1
    ((delta >= 500)) || return 1
    delta="$(pixel_delta "$images/styled.png" "$images/selected.png")" || return 1
    ((delta >= 40)) || return 1
    delta="$(pixel_delta "$images/selected.png" "$images/moved.png")" || return 1
    ((delta >= 40)) || return 1
    delta="$(pixel_delta "$images/resized.png" "$images/zoomed.png")" || return 1
    ((delta >= 200)) || return 1
    make_color_mask "$images/styled.png" "$mask_dir/magenta.png" '#ff00ff'
    make_color_mask "$images/styled.png" "$mask_dir/cyan.png" '#00ffff'
    make_color_mask "$images/clean.png" "$mask_dir/clean-magenta.png" '#ff00ff'
    make_color_mask "$images/clean.png" "$mask_dir/clean-cyan.png" '#00ffff'

    local magenta cyan clean_magenta clean_cyan
    magenta="$(mask_count "$mask_dir/magenta.png")"
    cyan="$(mask_count "$mask_dir/cyan.png")"
    clean_magenta="$(mask_count "$mask_dir/clean-magenta.png")"
    clean_cyan="$(mask_count "$mask_dir/clean-cyan.png")"

    # The source model itself has yellow CAD geometry at the polygon corner,
    # so selection appearance is verified by the selected-frame delta above.
    ((magenta >= 80 && cyan >= 40 && clean_magenta < 10 && clean_cyan < 10)) ||
	return 1

    if grep -h -Eq 'terminate called|Aborted|ASSERT|progressive whole-target source failed' \
	"$out/stdout.log" "$out/stderr.log"; then
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
    run="polygon-$backend"
    out="$artifact_dir/cases/$run"
    cache="$artifact_dir/caches/$run"
    events="$artifact_dir/events/$run.json"
    mkdir -p "$out/images" "$cache"
    write_events "$events" "$out/images"

    args=(--test-script "$events" --test-report "$out/report.json" "$database")
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

if ((failures == 0)) && [[ -d "$artifact_dir/cases/polygon-system/masks" &&
    -d "$artifact_dir/cases/polygon-osmesa/masks" ]]; then
    for color in magenta cyan; do
	system_mask="$artifact_dir/cases/polygon-system/masks/$color.png"
	osmesa_mask="$artifact_dir/cases/polygon-osmesa/masks/$color.png"
	system_count="$(mask_count "$system_mask")"
	osmesa_count="$(mask_count "$osmesa_mask")"
	larger=$((system_count > osmesa_count ? system_count : osmesa_count))
	smaller=$((system_count < osmesa_count ? system_count : osmesa_count))
	if ((smaller * 3 < larger)); then
	    echo "ERROR: $color coverage differs excessively: system=$system_count osmesa=$osmesa_count" >&2
	    failures=$((failures + 1))
	fi
    done
fi

echo "Artifacts: $artifact_dir"
echo "Failures: $failures"
exit "$failures"
