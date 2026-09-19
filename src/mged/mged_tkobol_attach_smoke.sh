#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
    echo "Usage: mged_tkobol_attach_smoke.sh <mged> <db> <workdir>" 1>&2
    exit 1
fi

MGED="$1"
DB="$2"
WORKDIR="$3"
OUT="${WORKDIR}/mged_tkobol_attach_smoke.png"
RT_OUT="${WORKDIR}/mged_tkobol_attach_rt.png"
RT_MOVED_OUT="${WORKDIR}/mged_tkobol_attach_rt_moved.png"
LOG="${WORKDIR}/mged_tkobol_attach_smoke.log"
HELPER="$(dirname "$0")/mged_tkobol_attach_smoke.tcl"
DISPLAY_VALUE="${DISPLAY:-}"

if [ -z "$DISPLAY_VALUE" ]; then
    echo "MGED Tk Obol attachment smoke requires an X display" 1>&2
    exit 1
fi

rm -f "$OUT" "$RT_OUT" "$RT_MOVED_OUT" "$LOG"

{
printf 'attach -d %s -t 1 -n .mged_direct tkobol
dm host
dm size 320 240
dm size
source %s
focus .mged_direct
rset g snap 1
rset g rh 100
rset g rv 100
set perspective_mode 0
puts OBOL_PERSPECTIVE_MODE_BEFORE
set perspective_mode
$obol_input_command key press {} F3 0 0 0
after 100
puts OBOL_PERSPECTIVE_MODE_AFTER
set perspective_mode
$obol_input_command key press {} F3 0 0 0
after 100
puts OBOL_PERSPECTIVE_MODE_RESTORED
set perspective_mode
set perspective_mode 1
puts OBOL_PERSPECTIVE_CYCLE_BEFORE
set perspective
$obol_input_command key press {} F6 0 0 0
after 100
puts OBOL_PERSPECTIVE_CYCLE_AFTER
set perspective
$obol_input_command key press {} F6 0 0 0
$obol_input_command key press {} F6 0 0 0
$obol_input_command key press {} F6 0 0 0
after 100
puts OBOL_PERSPECTIVE_CYCLE_RESTORED
set perspective
set perspective_mode 0
set faceplate 1
puts OBOL_FACEPLATE_BEFORE
set faceplate
$obol_input_command key press {} F7 0 0 0
after 100
puts OBOL_FACEPLATE_AFTER
set faceplate
$obol_input_command key press {} F7 0 0 0
after 100
puts OBOL_FACEPLATE_RESTORED
set faceplate
set orig_gui 1
puts OBOL_FACEPLATE_GUI_BEFORE
set orig_gui
$obol_input_command key press {} F8 0 0 0
after 100
puts OBOL_FACEPLATE_GUI_AFTER
set orig_gui
$obol_input_command key press {} F8 0 0 0
after 100
puts OBOL_FACEPLATE_GUI_RESTORED
set orig_gui
rset ax edit_draw 0
puts OBOL_EDIT_AXES_BEFORE
rset ax edit_draw
$obol_input_command key press e e 0 0 0
after 100
puts OBOL_EDIT_AXES_AFTER
rset ax edit_draw
$obol_input_command key press e e 0 0 0
after 100
puts OBOL_EDIT_AXES_RESTORED
rset ax edit_draw
puts OBOL_MOUSE_BEHAVIOR_BINDINGS
bind .mged_direct <ButtonPress-2>
bind .mged_direct <ButtonRelease-2>
set mouse_behavior z
rset g snap 0
puts OBOL_MOUSE_BEHAVIOR_SIZE_BEFORE
view size
event generate .mged_direct.__obol <ButtonPress-2> -x 80 -y 60
event generate .mged_direct.__obol <B2-Motion> -x 240 -y 180 -time 1230
event generate .mged_direct.__obol <ButtonRelease-2> -x 240 -y 180
after 100
puts OBOL_MOUSE_BEHAVIOR_SIZE_AFTER
view size
event generate .mged_direct.__obol <Motion> -x 280 -y 200 -time 1231
after 100
puts OBOL_MOUSE_BEHAVIOR_SIZE_CANCELLED
view size
set mouse_behavior d
rset g snap 1
puts OBOL_MOUSE_ZOOM_BINDINGS
bind .mged_direct <ButtonPress-1>
bind .mged_direct <ButtonPress-3>
puts OBOL_MOUSE_ZOOM_SIZE_BEFORE
view size
event generate .mged_direct.__obol <ButtonPress-1> -x 100 -y 100
after 100
puts OBOL_MOUSE_ZOOM_SIZE_AFTER_OUT
view size
event generate .mged_direct.__obol <ButtonPress-3> -x 100 -y 100
after 100
puts OBOL_MOUSE_ZOOM_SIZE_RESTORED
view size
puts OBOL_CENTER_BEFORE
view center
puts OBOL_TRANSLATE_PRESS_BINDING
bind .mged_direct <Shift-ButtonPress-1>
puts OBOL_TRANSLATE_RELEASE_BINDING
bind .mged_direct <Shift-ButtonRelease-1>
puts OBOL_MODIFIER_RELEASE_BINDING
bind .mged_direct <KeyRelease-Shift_L>
event generate .mged_direct.__obol <Shift-ButtonPress-1> -x 100 -y 100
event generate .mged_direct.__obol <Shift-B1-Motion> -x 150 -y 100 -time 1232
after 100
puts OBOL_CENTER_AFTER
view center
event generate .mged_direct.__obol <Shift-ButtonRelease-1> -x 150 -y 100
event generate .mged_direct.__obol <Motion> -x 200 -y 100 -time 1233
after 100
puts OBOL_CENTER_CANCELLED
view center
event generate .mged_direct.__obol <Shift-ButtonPress-1> -x 150 -y 100
$obol_input_command key release {} Shift_L 0 0 0
event generate .mged_direct.__obol <B1-Motion> -x 200 -y 100 -time 1234
after 100
puts OBOL_CENTER_MODIFIER_CANCELLED
view center
event generate .mged_direct.__obol <ButtonRelease-1> -x 200 -y 100
puts OBOL_AET_BEFORE
view ae
puts OBOL_ROTATE_PRESS_BINDING
bind .mged_direct <Control-ButtonPress-1>
puts OBOL_ROTATE_RELEASE_BINDING
bind .mged_direct <Control-ButtonRelease-1>
event generate .mged_direct.__obol <Control-ButtonPress-1> -x 100 -y 100
event generate .mged_direct.__obol <Control-B1-Motion> -x 150 -y 100 -time 1235
after 100
puts OBOL_AET_AFTER
view ae
event generate .mged_direct.__obol <Control-ButtonRelease-1> -x 150 -y 100
event generate .mged_direct.__obol <Motion> -x 200 -y 100 -time 1236
after 100
puts OBOL_AET_CANCELLED
view ae
event generate .mged_direct.__obol <Control-ButtonPress-1> -x 150 -y 100
event generate .mged_direct.__obol <FocusOut>
event generate .mged_direct.__obol <B1-Motion> -x 200 -y 100 -time 1237
after 100
puts OBOL_AET_FOCUS_CANCELLED
view ae
event generate .mged_direct.__obol <FocusIn>
event generate .mged_direct.__obol <ButtonRelease-1> -x 200 -y 100
puts OBOL_SIZE_BEFORE
view size
puts OBOL_SCALE_PRESS_BINDING
bind .mged_direct <Shift-Control-ButtonPress-1>
puts OBOL_SCALE_RELEASE_BINDING
bind .mged_direct <Shift-Control-ButtonRelease-1>
event generate .mged_direct.__obol <Shift-Control-ButtonPress-1> -x 100 -y 100
event generate .mged_direct.__obol <Shift-Control-B1-Motion> -x 150 -y 100 -time 1236
after 100
puts OBOL_SIZE_AFTER
view size
event generate .mged_direct.__obol <Shift-Control-ButtonRelease-1> -x 150 -y 100
event generate .mged_direct.__obol <Motion> -x 200 -y 100 -time 1237
after 100
puts OBOL_SIZE_CANCELLED
view size
puts OBOL_CONSTRAINED_BINDINGS
bind .mged_direct <Alt-Shift-ButtonPress-1>
bind .mged_direct <Alt-Shift-ButtonRelease-1>
bind .mged_direct <Alt-Control-ButtonPress-2>
bind .mged_direct <Alt-Control-ButtonRelease-2>
bind .mged_direct <Alt-Shift-Control-ButtonPress-3>
bind .mged_direct <Alt-Shift-Control-ButtonRelease-3>
puts OBOL_CONSTRAINED_CENTER_BEFORE
view center
event generate .mged_direct.__obol <Alt-Shift-ButtonPress-1> -x 100 -y 100
event generate .mged_direct.__obol <Alt-Shift-B1-Motion> -x 140 -y 100 -time 1238
after 100
puts OBOL_CONSTRAINED_CENTER_AFTER
view center
event generate .mged_direct.__obol <Alt-Shift-ButtonRelease-1> -x 140 -y 100
event generate .mged_direct.__obol <Motion> -x 180 -y 100 -time 1239
after 100
puts OBOL_CONSTRAINED_CENTER_CANCELLED
view center
puts OBOL_CONSTRAINED_AET_BEFORE
view ae
event generate .mged_direct.__obol <Alt-Control-ButtonPress-2> -x 100 -y 100
event generate .mged_direct.__obol <Alt-Control-B2-Motion> -x 140 -y 100 -time 1240
after 100
puts OBOL_CONSTRAINED_AET_AFTER
view ae
event generate .mged_direct.__obol <Alt-Control-ButtonRelease-2> -x 140 -y 100
event generate .mged_direct.__obol <Motion> -x 180 -y 100 -time 1241
after 100
puts OBOL_CONSTRAINED_AET_CANCELLED
view ae
puts OBOL_CONSTRAINED_SIZE_BEFORE
view size
event generate .mged_direct.__obol <Alt-Shift-Control-ButtonPress-3> -x 100 -y 100
event generate .mged_direct.__obol <Alt-Shift-Control-B3-Motion> -x 140 -y 100 -time 1242
after 100
puts OBOL_CONSTRAINED_SIZE_AFTER
view size
event generate .mged_direct.__obol <Alt-Shift-Control-ButtonRelease-3> -x 140 -y 100
event generate .mged_direct.__obol <Motion> -x 180 -y 100 -time 1243
after 100
puts OBOL_CONSTRAINED_SIZE_CANCELLED
view size
rset g snap 0
adc draw 1
set transform a
puts OBOL_ADC_BINDINGS
bind .mged_direct <Shift-ButtonPress-2>
bind .mged_direct <Shift-ButtonRelease-2>
bind .mged_direct <Control-ButtonPress-2>
bind .mged_direct <Control-ButtonRelease-2>
bind .mged_direct <Shift-Control-ButtonPress-3>
bind .mged_direct <Shift-Control-ButtonRelease-3>
bind .mged_direct <Alt-Shift-ButtonPress-1>
bind .mged_direct <Alt-Shift-ButtonRelease-1>
bind .mged_direct <Alt-Control-ButtonPress-1>
bind .mged_direct <Alt-Control-ButtonRelease-1>
bind .mged_direct <Alt-Shift-Control-ButtonPress-3>
bind .mged_direct <Alt-Shift-Control-ButtonRelease-3>
adc xyz 0 0 0
puts OBOL_ADC_TRANSLATE_BEFORE
adc xyz
event generate .mged_direct.__obol <Shift-ButtonPress-2> -x 100 -y 100
event generate .mged_direct.__obol <Shift-B2-Motion> -x 160 -y 110 -time 1244
after 100
puts OBOL_ADC_TRANSLATE_AFTER
adc xyz
event generate .mged_direct.__obol <Shift-ButtonRelease-2> -x 160 -y 110
event generate .mged_direct.__obol <Motion> -x 200 -y 120 -time 1245
after 100
puts OBOL_ADC_TRANSLATE_CANCELLED
adc xyz
adc a2 0
puts OBOL_ADC_ANGLE_BEFORE
adc a2
event generate .mged_direct.__obol <Control-ButtonPress-2> -x 120 -y 100
event generate .mged_direct.__obol <Control-B2-Motion> -x 100 -y 150 -time 1246
after 100
puts OBOL_ADC_ANGLE_AFTER
adc a2
event generate .mged_direct.__obol <Control-ButtonRelease-2> -x 100 -y 150
event generate .mged_direct.__obol <Motion> -x 150 -y 150 -time 1247
after 100
puts OBOL_ADC_ANGLE_CANCELLED
adc a2
adc dst 0
puts OBOL_ADC_DISTANCE_BEFORE
adc dst
event generate .mged_direct.__obol <Shift-Control-ButtonPress-3> -x 100 -y 100
event generate .mged_direct.__obol <Shift-Control-B3-Motion> -x 160 -y 150 -time 1248
after 100
puts OBOL_ADC_DISTANCE_AFTER
adc dst
event generate .mged_direct.__obol <Shift-Control-ButtonRelease-3> -x 160 -y 150
event generate .mged_direct.__obol <Motion> -x 200 -y 180 -time 1249
after 100
puts OBOL_ADC_DISTANCE_CANCELLED
adc dst
puts OBOL_ADC_CONSTRAINED_TRANSLATE_BEFORE
adc xyz
event generate .mged_direct.__obol <Alt-Shift-ButtonPress-1> -x 100 -y 100
event generate .mged_direct.__obol <Alt-Shift-B1-Motion> -x 150 -y 140 -time 1250
after 100
puts OBOL_ADC_CONSTRAINED_TRANSLATE_AFTER
adc xyz
event generate .mged_direct.__obol <Alt-Shift-ButtonRelease-1> -x 150 -y 140
event generate .mged_direct.__obol <Motion> -x 180 -y 160 -time 1251
after 100
puts OBOL_ADC_CONSTRAINED_TRANSLATE_CANCELLED
adc xyz
puts OBOL_ADC_CONSTRAINED_ANGLE_BEFORE
adc a1
event generate .mged_direct.__obol <Alt-Control-ButtonPress-1> -x 110 -y 100
event generate .mged_direct.__obol <Alt-Control-B1-Motion> -x 150 -y 140 -time 1252
after 100
puts OBOL_ADC_CONSTRAINED_ANGLE_AFTER
adc a1
event generate .mged_direct.__obol <Alt-Control-ButtonRelease-1> -x 150 -y 140
event generate .mged_direct.__obol <Motion> -x 180 -y 160 -time 1253
after 100
puts OBOL_ADC_CONSTRAINED_ANGLE_CANCELLED
adc a1
puts OBOL_ADC_CONSTRAINED_DISTANCE_BEFORE
adc dst
event generate .mged_direct.__obol <Alt-Shift-Control-ButtonPress-3> -x 100 -y 100
event generate .mged_direct.__obol <Alt-Shift-Control-B3-Motion> -x 150 -y 150 -time 1254
after 100
puts OBOL_ADC_CONSTRAINED_DISTANCE_AFTER
adc dst
event generate .mged_direct.__obol <Alt-Shift-Control-ButtonRelease-3> -x 150 -y 150
event generate .mged_direct.__obol <Motion> -x 190 -y 170 -time 1255
after 100
puts OBOL_ADC_CONSTRAINED_DISTANCE_CANCELLED
adc dst
set transform v
adc draw 0
draw all.g
autoview
refresh
' "$DISPLAY_VALUE" "$HELPER"
# The command block above deliberately covers many application-owned input
# bindings.  Let MGED consume it and return to Tk's real event loop before the
# next producer phase; otherwise a full pipe can keep Configure/redraw events
# queued behind later commands and make a passive `dm size` value look settled
# before the native Togl child is actually resized.
sleep 6
printf 'autoview
refresh
'
sleep 1
printf 'puts OBOL_RESIZE_CAMERA_BEFORE
view center
view size
view ae
puts OBOL_RESIZE_WIDGET_BEFORE
winfo width .mged_direct.__obol
winfo height .mged_direct.__obol
wm geometry .mged_direct 520x360
update idletasks
wm geometry .mged_direct 410x300
update idletasks
wm geometry .mged_direct 600x400
update idletasks
wm geometry .mged_direct 470x330
update idletasks
'
sleep 1
printf 'puts OBOL_RESIZE_CAMERA_AFTER
view center
view size
view ae
puts OBOL_NATIVE_WIDGET_SIZE
winfo width .mged_direct.__obol
winfo height .mged_direct.__obol
puts OBOL_NATIVE_TOPLEVEL_SIZE
winfo width .mged_direct
winfo height .mged_direct
screengrab %s
dm set render.rt.samples 1
dm set render.rt.preview_scale 4
dm set render.rt.frame_budget_ms 16
dm renderer rt
' "$OUT"
sleep 4
printf 'screengrab %s
puts OBOL_RT_REVISION_0
dm get render.rt.geometry_revision
dm get render.rt.presentation_revision
ae 35 25
' "$RT_OUT"
sleep 4
printf 'screengrab %s
puts OBOL_RT_REVISION_1
dm get render.rt.geometry_revision
dm get render.rt.presentation_revision
dm renderer hw
release .mged_direct
quit
' "$RT_MOVED_OUT"
} | "$MGED" -c -a nu -r "$DB" > "$LOG" 2>&1

if ! grep -q '^ATTACHING tkobol (Tk Obol display endpoint)$' "$LOG" ||
    ! grep -qx 'tk-gl' "$LOG" || ! grep -qx '320 240' "$LOG"; then
    echo "MGED did not attach and resize the direct Tk Obol endpoint" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if grep -q 'attach(tkobol): BAD\|release: .* not found' "$LOG" ||
    [ ! -s "$OUT" ] || [ "$(wc -c < "$OUT")" -le 1000 ] ||
    [ ! -s "$RT_OUT" ] ||
    [ ! -s "$RT_MOVED_OUT" ] || [ "$(wc -c < "$RT_OUT")" -le 1000 ] ||
    [ "$(wc -c < "$RT_MOVED_OUT")" -le 1000 ]; then
    echo "MGED direct Tk Obol endpoint lifecycle failed" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

native_widget_width=$(awk '/^OBOL_NATIVE_WIDGET_SIZE$/ { getline; print; exit }' "$LOG")
native_widget_height=$(awk '/^OBOL_NATIVE_WIDGET_SIZE$/ { getline; getline; print; exit }' "$LOG")
native_toplevel_width=$(awk '/^OBOL_NATIVE_TOPLEVEL_SIZE$/ { getline; print; exit }' "$LOG")
native_toplevel_height=$(awk '/^OBOL_NATIVE_TOPLEVEL_SIZE$/ { getline; getline; print; exit }' "$LOG")
resize_widget_before_width=$(awk '/^OBOL_RESIZE_WIDGET_BEFORE$/ { getline; print; exit }' "$LOG")
resize_widget_before_height=$(awk '/^OBOL_RESIZE_WIDGET_BEFORE$/ { getline; getline; print; exit }' "$LOG")
resize_camera_before=$(awk '/^OBOL_RESIZE_CAMERA_BEFORE$/ { getline; c=$0; getline; s=$0; getline; a=$0; print c "|" s "|" a; exit }' "$LOG")
resize_camera_after=$(awk '/^OBOL_RESIZE_CAMERA_AFTER$/ { getline; c=$0; getline; s=$0; getline; a=$0; print c "|" s "|" a; exit }' "$LOG")
if [ "$resize_widget_before_width" != "320" ] ||
    [ "$resize_widget_before_height" != "240" ] ||
    [ "$native_widget_width" != "470" ] ||
    [ "$native_widget_height" != "330" ] ||
    [ "$native_toplevel_width" != "470" ] ||
    [ "$native_toplevel_height" != "330" ] ||
    [ -z "$resize_camera_before" ] ||
    [ "$resize_camera_after" != "$resize_camera_before" ]; then
    echo "MGED direct Tk Obol native window did not follow endpoint resize" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

rt_geometry_0=$(awk '/^OBOL_RT_REVISION_0$/ { getline; print; exit }' "$LOG")
rt_presentation_0=$(awk '/^OBOL_RT_REVISION_0$/ { getline; getline; print; exit }' "$LOG")
rt_geometry_1=$(awk '/^OBOL_RT_REVISION_1$/ { getline; print; exit }' "$LOG")
rt_presentation_1=$(awk '/^OBOL_RT_REVISION_1$/ { getline; getline; print; exit }' "$LOG")
if [ -z "$rt_geometry_0" ] || [ "$rt_geometry_0" -le 0 ] ||
    [ "$rt_geometry_1" != "$rt_geometry_0" ] ||
    [ -z "$rt_presentation_0" ] || [ -z "$rt_presentation_1" ] ||
    [ "$rt_presentation_1" -le "$rt_presentation_0" ] ||
    cmp -s "$RT_OUT" "$RT_MOVED_OUT"; then
    echo "MGED retained RT did not present or retain geometry across a native camera change" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

native_key_bindings=$(awk '/^OBOL_NATIVE_KEY_BINDINGS$/ { getline; print; exit }' "$LOG")
perspective_mode_before=$(awk '/^OBOL_PERSPECTIVE_MODE_BEFORE$/ { getline; print; exit }' "$LOG")
perspective_mode_after=$(awk '/^OBOL_PERSPECTIVE_MODE_AFTER$/ { getline; print; exit }' "$LOG")
perspective_mode_restored=$(awk '/^OBOL_PERSPECTIVE_MODE_RESTORED$/ { getline; print; exit }' "$LOG")
perspective_cycle_before=$(awk '/^OBOL_PERSPECTIVE_CYCLE_BEFORE$/ { getline; print; exit }' "$LOG")
perspective_cycle_after=$(awk '/^OBOL_PERSPECTIVE_CYCLE_AFTER$/ { getline; print; exit }' "$LOG")
perspective_cycle_restored=$(awk '/^OBOL_PERSPECTIVE_CYCLE_RESTORED$/ { getline; print; exit }' "$LOG")
if [ "$native_key_bindings" != "||||" ] ||
    [ "$perspective_mode_before" != "0" ] || [ "$perspective_mode_after" != "1" ] ||
    [ "$perspective_mode_restored" != "0" ] ||
    [ -z "$perspective_cycle_before" ] ||
    [ "$perspective_cycle_before" = "$perspective_cycle_after" ] ||
    [ "$perspective_cycle_restored" != "$perspective_cycle_before" ]; then
    echo "MGED perspective shortcuts did not use their application action layer" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

faceplate_before=$(awk '/^OBOL_FACEPLATE_BEFORE$/ { getline; print; exit }' "$LOG")
faceplate_after=$(awk '/^OBOL_FACEPLATE_AFTER$/ { getline; print; exit }' "$LOG")
faceplate_restored=$(awk '/^OBOL_FACEPLATE_RESTORED$/ { getline; print; exit }' "$LOG")
if [ "$faceplate_before" != "1" ] ||
    [ "$faceplate_after" != "0" ] || [ "$faceplate_restored" != "1" ]; then
    echo "MGED faceplate toggle did not use its application action layer" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

faceplate_gui_before=$(awk '/^OBOL_FACEPLATE_GUI_BEFORE$/ { getline; print; exit }' "$LOG")
faceplate_gui_after=$(awk '/^OBOL_FACEPLATE_GUI_AFTER$/ { getline; print; exit }' "$LOG")
faceplate_gui_restored=$(awk '/^OBOL_FACEPLATE_GUI_RESTORED$/ { getline; print; exit }' "$LOG")
if [ "$faceplate_gui_before" != "1" ] ||
    [ "$faceplate_gui_after" != "0" ] || [ "$faceplate_gui_restored" != "1" ]; then
    echo "MGED faceplate GUI toggle did not use its application action layer" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

edit_axes_before=$(awk '/^OBOL_EDIT_AXES_BEFORE$/ { getline; print; exit }' "$LOG")
edit_axes_after=$(awk '/^OBOL_EDIT_AXES_AFTER$/ { getline; print; exit }' "$LOG")
edit_axes_restored=$(awk '/^OBOL_EDIT_AXES_RESTORED$/ { getline; print; exit }' "$LOG")
if [ "$edit_axes_before" != "0" ] ||
    [ "$edit_axes_after" != "1" ] || [ "$edit_axes_restored" != "0" ]; then
    echo "MGED edit-axes toggle did not use its application action layer" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

mouse_behavior_bindings=$(awk '/^OBOL_MOUSE_BEHAVIOR_BINDINGS$/ { getline; first=$0; getline; print first "|" $0; exit }' "$LOG")
mouse_behavior_size_before=$(awk '/^OBOL_MOUSE_BEHAVIOR_SIZE_BEFORE$/ { getline; print; exit }' "$LOG")
mouse_behavior_size_after=$(awk '/^OBOL_MOUSE_BEHAVIOR_SIZE_AFTER$/ { getline; print; exit }' "$LOG")
mouse_behavior_size_cancelled=$(awk '/^OBOL_MOUSE_BEHAVIOR_SIZE_CANCELLED$/ { getline; print; exit }' "$LOG")
if [ "$mouse_behavior_bindings" != "break|break" ] ||
    [ -z "$mouse_behavior_size_before" ] ||
    [ "$mouse_behavior_size_before" = "$mouse_behavior_size_after" ] ||
    [ "$mouse_behavior_size_cancelled" != "$mouse_behavior_size_after" ]; then
    echo "MGED mouse-behavior lifecycle was missing, duplicated, or not cancelled" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

mouse_zoom_bindings=$(awk '/^OBOL_MOUSE_ZOOM_BINDINGS$/ { getline; first=$0; getline; print first "|" $0; exit }' "$LOG")
mouse_zoom_size_before=$(awk '/^OBOL_MOUSE_ZOOM_SIZE_BEFORE$/ { getline; print; exit }' "$LOG")
mouse_zoom_size_after_out=$(awk '/^OBOL_MOUSE_ZOOM_SIZE_AFTER_OUT$/ { getline; print; exit }' "$LOG")
mouse_zoom_size_restored=$(awk '/^OBOL_MOUSE_ZOOM_SIZE_RESTORED$/ { getline; print; exit }' "$LOG")
if [ "$mouse_zoom_bindings" != "break|break" ] ||
    [ -z "$mouse_zoom_size_before" ] ||
    [ "$mouse_zoom_size_before" = "$mouse_zoom_size_after_out" ] ||
    [ "$mouse_zoom_size_restored" != "$mouse_zoom_size_before" ]; then
    echo "MGED unmodified zoom actions were missing or duplicated" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

before_center=$(awk '/^OBOL_CENTER_BEFORE$/ { getline; print; exit }' "$LOG")
after_center=$(awk '/^OBOL_CENTER_AFTER$/ { getline; print; exit }' "$LOG")
cancelled_center=$(awk '/^OBOL_CENTER_CANCELLED$/ { getline; print; exit }' "$LOG")
modifier_cancelled_center=$(awk '/^OBOL_CENTER_MODIFIER_CANCELLED$/ { getline; print; exit }' "$LOG")
if [ -z "$before_center" ] ||
    [ "$before_center" = "$after_center" ] ||
    [ "$cancelled_center" != "$after_center" ] ||
    [ "$modifier_cancelled_center" != "$after_center" ]; then
    echo "MGED Tk Obol translation was missing, duplicated, or not cancelled" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

translate_press_binding=$(awk '/^OBOL_TRANSLATE_PRESS_BINDING$/ { getline; print; exit }' "$LOG")
translate_release_binding=$(awk '/^OBOL_TRANSLATE_RELEASE_BINDING$/ { getline; print; exit }' "$LOG")
modifier_release_binding=$(awk '/^OBOL_MODIFIER_RELEASE_BINDING$/ { getline; print; exit }' "$LOG")
if [ "$translate_press_binding" != "break" ] ||
    [ "$translate_release_binding" != "break" ] ||
    [ "$modifier_release_binding" != "break" ]; then
    echo "MGED retained a duplicate native binding for semantic translation" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

before_aet=$(awk '/^OBOL_AET_BEFORE$/ { getline; print; exit }' "$LOG")
after_aet=$(awk '/^OBOL_AET_AFTER$/ { getline; print; exit }' "$LOG")
cancelled_aet=$(awk '/^OBOL_AET_CANCELLED$/ { getline; print; exit }' "$LOG")
focus_cancelled_aet=$(awk '/^OBOL_AET_FOCUS_CANCELLED$/ { getline; print; exit }' "$LOG")
if [ "$before_aet" != "270 90 0" ] ||
    [ "$after_aet" != "180 77.5 90" ] ||
    [ "$cancelled_aet" != "$after_aet" ] ||
    [ "$focus_cancelled_aet" != "$after_aet" ]; then
    echo "MGED Tk Obol rotate gesture was missing, duplicated, or not cancelled" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

press_binding=$(awk '/^OBOL_ROTATE_PRESS_BINDING$/ { getline; print; exit }' "$LOG")
release_binding=$(awk '/^OBOL_ROTATE_RELEASE_BINDING$/ { getline; print; exit }' "$LOG")
if [ "$press_binding" != "break" ] || [ "$release_binding" != "break" ]; then
    echo "MGED retained a duplicate native binding for semantic rotate" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

before_size=$(awk '/^OBOL_SIZE_BEFORE$/ { getline; print; exit }' "$LOG")
after_size=$(awk '/^OBOL_SIZE_AFTER$/ { getline; print; exit }' "$LOG")
cancelled_size=$(awk '/^OBOL_SIZE_CANCELLED$/ { getline; print; exit }' "$LOG")
if [ -z "$before_size" ] || [ "$before_size" = "$after_size" ] ||
    [ "$cancelled_size" != "$after_size" ]; then
    echo "MGED Tk Obol scale was missing, duplicated, or not cancelled" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

scale_press_binding=$(awk '/^OBOL_SCALE_PRESS_BINDING$/ { getline; print; exit }' "$LOG")
scale_release_binding=$(awk '/^OBOL_SCALE_RELEASE_BINDING$/ { getline; print; exit }' "$LOG")
if [ "$scale_press_binding" != "break" ] ||
    [ "$scale_release_binding" != "break" ]; then
    echo "MGED retained a duplicate native binding for semantic scale" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

constrained_bindings=$(awk '/^OBOL_CONSTRAINED_BINDINGS$/ { out=""; for (i=0; i<6; i++) { getline; out=out (i ? "|" : "") $0 }; print out; exit }' "$LOG")
if [ "$constrained_bindings" != "break|break|break|break|break|break" ]; then
    echo "MGED retained duplicate native constrained-mode bindings" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

constrained_center_before=$(awk '/^OBOL_CONSTRAINED_CENTER_BEFORE$/ { getline; print; exit }' "$LOG")
constrained_center_after=$(awk '/^OBOL_CONSTRAINED_CENTER_AFTER$/ { getline; print; exit }' "$LOG")
constrained_center_cancelled=$(awk '/^OBOL_CONSTRAINED_CENTER_CANCELLED$/ { getline; print; exit }' "$LOG")
if [ -z "$constrained_center_before" ] ||
    [ "$constrained_center_before" = "$constrained_center_after" ] ||
    [ "$constrained_center_cancelled" != "$constrained_center_after" ]; then
    echo "MGED constrained translation was missing or not cancelled" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

constrained_aet_before=$(awk '/^OBOL_CONSTRAINED_AET_BEFORE$/ { getline; print; exit }' "$LOG")
constrained_aet_after=$(awk '/^OBOL_CONSTRAINED_AET_AFTER$/ { getline; print; exit }' "$LOG")
constrained_aet_cancelled=$(awk '/^OBOL_CONSTRAINED_AET_CANCELLED$/ { getline; print; exit }' "$LOG")
if [ -z "$constrained_aet_before" ] ||
    [ "$constrained_aet_before" = "$constrained_aet_after" ] ||
    [ "$constrained_aet_cancelled" != "$constrained_aet_after" ]; then
    echo "MGED constrained rotation was missing or not cancelled" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

constrained_size_before=$(awk '/^OBOL_CONSTRAINED_SIZE_BEFORE$/ { getline; print; exit }' "$LOG")
constrained_size_after=$(awk '/^OBOL_CONSTRAINED_SIZE_AFTER$/ { getline; print; exit }' "$LOG")
constrained_size_cancelled=$(awk '/^OBOL_CONSTRAINED_SIZE_CANCELLED$/ { getline; print; exit }' "$LOG")
if [ -z "$constrained_size_before" ] ||
    [ "$constrained_size_before" = "$constrained_size_after" ] ||
    [ "$constrained_size_cancelled" != "$constrained_size_after" ]; then
    echo "MGED constrained scale was missing or not cancelled" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

adc_bindings=$(awk '/^OBOL_ADC_BINDINGS$/ { out=""; for (i=0; i<12; i++) { getline; out=out (i ? "|" : "") $0 }; print out; exit }' "$LOG")
if [ "$adc_bindings" != "break|break|break|break|break|break|break|break|break|break|break|break" ]; then
    echo "MGED retained duplicate native ADC-mode bindings" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

for gesture in TRANSLATE ANGLE DISTANCE CONSTRAINED_TRANSLATE CONSTRAINED_ANGLE CONSTRAINED_DISTANCE; do
    before=$(awk -v marker="OBOL_ADC_${gesture}_BEFORE" '$0 == marker { getline; print; exit }' "$LOG")
    after=$(awk -v marker="OBOL_ADC_${gesture}_AFTER" '$0 == marker { getline; print; exit }' "$LOG")
    cancelled=$(awk -v marker="OBOL_ADC_${gesture}_CANCELLED" '$0 == marker { getline; print; exit }' "$LOG")
    if [ -z "$before" ] || [ "$before" = "$after" ] ||
        [ "$cancelled" != "$after" ]; then
        echo "MGED ADC ${gesture} gesture was missing or not cancelled" 1>&2
        cat "$LOG" 1>&2
        exit 1
    fi
done

exit 0
