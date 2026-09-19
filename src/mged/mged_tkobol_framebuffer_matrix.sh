#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
    echo "Usage: mged_tkobol_framebuffer_matrix.sh <mged> <db> <workdir>" 1>&2
    exit 1
fi

MGED="$1"
DB="$2"
WORKDIR="$3"
DISPLAY_VALUE="${DISPLAY:-}"

if [ -z "$DISPLAY_VALUE" ]; then
    echo "MGED Tk Obol framebuffer matrix requires an X display" 1>&2
    exit 125
fi
if [ ! -x "$MGED" ] || [ ! -s "$DB" ]; then
    echo "MGED Tk Obol framebuffer matrix is missing its executable or database" 1>&2
    exit 1
fi
for tool in identify compare; do
    if ! command -v "$tool" >/dev/null 2>&1; then
	echo "MGED Tk Obol framebuffer matrix requires $tool" 1>&2
	exit 125
    fi
done

RUN_DIR=$(mktemp -d "${WORKDIR}/mged_tkobol_framebuffer.XXXXXX")
CAD="$RUN_DIR/cad.png"
UNDERLAY="$RUN_DIR/underlay.png"
OVERLAY="$RUN_DIR/overlay.png"
INTERLAY="$RUN_DIR/interlay.png"
RESIZED="$RUN_DIR/resized.png"
OFF="$RUN_DIR/off.png"
LOG="$RUN_DIR/mged_tkobol_framebuffer.log"

cleanup_success()
{
    if [ "${BRLCAD_KEEP_TEST_ARTIFACTS:-0}" = 1 ]; then
	echo "MGED Tk Obol framebuffer artifacts: $RUN_DIR"
	return
    fi
    rm -f "$CAD" "$UNDERLAY" "$OVERLAY" "$INTERLAY" "$RESIZED" \
	"$OFF" "$LOG"
    rmdir "$RUN_DIR"
}

pixel_delta()
{
    output=$(compare -metric AE "$1" "$2" null: 2>&1 || true)
    output=${output%% *}
    case "$output" in
	''|*[!0-9]*) return 1 ;;
    esac
    printf '%s' "$output"
}

fail()
{
    echo "$1" 1>&2
    echo "MGED Tk Obol framebuffer artifacts retained: $RUN_DIR" 1>&2
    cat "$LOG" 1>&2
    exit 1
}

{
printf 'attach -d %s -t 1 -n .mged_fb_matrix tkobol
wm geometry .mged_fb_matrix 480x340
update idletasks
view faceplate center_dot 1
draw -m1 all
ae 90 0
autoview
refresh
' "$DISPLAY_VALUE"
sleep 6
printf 'screengrab %s
puts OBOL_FB_INITIAL_WIDGET
winfo width .mged_fb_matrix.__obol
winfo height .mged_fb_matrix.__obol
ert -P 1 -H 16
' "$CAD"
sleep 7
printf 'puts OBOL_FB_UNDERLAY_MODE
view faceplate fb
refresh
screengrab %s
view faceplate fb 1
puts OBOL_FB_OVERLAY_MODE
view faceplate fb
refresh
screengrab %s
view faceplate fb 3
puts OBOL_FB_INTERLAY_MODE
view faceplate fb
refresh
screengrab %s
wm geometry .mged_fb_matrix 640x460
update idletasks
ae 35 25
autoview
refresh
' "$UNDERLAY" "$OVERLAY" "$INTERLAY"
sleep 1
printf 'puts OBOL_FB_RESIZED_WIDGET
winfo width .mged_fb_matrix.__obol
winfo height .mged_fb_matrix.__obol
ert -P 1 -H 16
'
sleep 15
printf 'refresh
screengrab %s
view faceplate fb 0
puts OBOL_FB_OFF_MODE
view faceplate fb
refresh
screengrab %s
release .mged_fb_matrix
quit
' "$RESIZED" "$OFF"
} | "$MGED" -c -a nu -r "$DB" > "$LOG" 2>&1

if ! grep -q '^ATTACHING tkobol (Tk Obol display endpoint)$' "$LOG" ||
    [ "$(grep -c 'Raytrace complete' "$LOG")" -lt 2 ] ||
    grep -Eq 'terminate called|Aborted|rt_gettrees.*FAILED|progressive whole-target source failed' "$LOG"; then
    fail "MGED Tk Obol framebuffer lifecycle or raytrace failed"
fi

for image in "$CAD" "$UNDERLAY" "$OVERLAY" "$INTERLAY" "$RESIZED" "$OFF"; do
    [ -s "$image" ] || fail "MGED Tk Obol framebuffer capture is missing: $image"
done

initial_size=$(identify -format '%w %h' "$CAD")
resized_size=$(identify -format '%w %h' "$RESIZED")
set -- $initial_size
initial_width=$1
initial_height=$2
set -- $resized_size
resized_width=$1
resized_height=$2
if [ "$resized_width" -le "$initial_width" ] ||
    [ "$resized_height" -le "$initial_height" ]; then
    fail "MGED Tk Obol framebuffer did not follow the native resize"
fi

initial_widget_width=$(awk '/^OBOL_FB_INITIAL_WIDGET$/ { getline; print; exit }' "$LOG")
initial_widget_height=$(awk '/^OBOL_FB_INITIAL_WIDGET$/ { getline; getline; print; exit }' "$LOG")
resized_widget_width=$(awk '/^OBOL_FB_RESIZED_WIDGET$/ { getline; print; exit }' "$LOG")
resized_widget_height=$(awk '/^OBOL_FB_RESIZED_WIDGET$/ { getline; getline; print; exit }' "$LOG")
if [ "$initial_widget_width" != "$initial_width" ] ||
    [ "$initial_widget_height" != "$initial_height" ] ||
    [ "$resized_widget_width" != "$resized_width" ] ||
    [ "$resized_widget_height" != "$resized_height" ]; then
    fail "MGED Tk Obol capture dimensions did not match the visible canvas"
fi

underlay_mode=$(awk '/^OBOL_FB_UNDERLAY_MODE$/ { getline; print; exit }' "$LOG")
overlay_mode=$(awk '/^OBOL_FB_OVERLAY_MODE$/ { getline; print; exit }' "$LOG")
interlay_mode=$(awk '/^OBOL_FB_INTERLAY_MODE$/ { getline; print; exit }' "$LOG")
off_mode=$(awk '/^OBOL_FB_OFF_MODE$/ { getline; print; exit }' "$LOG")
if [ "$underlay_mode" != 2 ] || [ "$overlay_mode" != 1 ] ||
    [ "$interlay_mode" != 3 ] || [ "$off_mode" != 0 ]; then
    fail "MGED Tk Obol framebuffer modes did not round-trip through the view state"
fi

cad_underlay=$(pixel_delta "$CAD" "$UNDERLAY") ||
    fail "MGED Tk Obol CAD/underlay comparison failed"
underlay_overlay=$(pixel_delta "$UNDERLAY" "$OVERLAY") ||
    fail "MGED Tk Obol underlay/overlay comparison failed"
overlay_interlay=$(pixel_delta "$OVERLAY" "$INTERLAY") ||
    fail "MGED Tk Obol overlay/interlay comparison failed"
resized_off=$(pixel_delta "$RESIZED" "$OFF") ||
    fail "MGED Tk Obol resized/off comparison failed"

if [ "$cad_underlay" -lt 1000 ] || [ "$underlay_overlay" -lt 1000 ] ||
    [ "$overlay_interlay" -lt 100 ] || [ "$resized_off" -lt 1000 ]; then
    fail "MGED Tk Obol framebuffer modes were not visibly distinct"
fi

echo "MGED Tk Obol framebuffer PASS size=${initial_width}x${initial_height}->${resized_width}x${resized_height} deltas=${cad_underlay},${underlay_overlay},${overlay_interlay},${resized_off}"
cleanup_success
exit 0
