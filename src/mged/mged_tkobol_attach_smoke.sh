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
LOG="${WORKDIR}/mged_tkobol_attach_smoke.log"
DISPLAY_VALUE="${DISPLAY:-}"

if [ -z "$DISPLAY_VALUE" ]; then
    echo "MGED Tk Obol attachment smoke requires an X display" 1>&2
    exit 1
fi

rm -f "$OUT" "$LOG"

printf 'attach -d %s -t 1 -n .mged_direct tkobol
dm host
dm size 320 240
dm size
rset g snap 1
rset g rh 100
rset g rv 100
puts OBOL_CENTER_BEFORE
view center
event generate .mged_direct.__obol <Shift-ButtonPress-1> -x 100 -y 100
event generate .mged_direct.__obol <Shift-B1-Motion> -x 150 -y 100
after 100
puts OBOL_CENTER_AFTER
view center
puts OBOL_AET_BEFORE
view ae
dm am r 100 100
event generate .mged_direct.__obol <Motion> -x 150 -y 100 -time 1234
after 100
puts OBOL_AET_AFTER
view ae
dm idle
draw all.g
autoview
refresh
screengrab %s
release .mged_direct
quit
' "$DISPLAY_VALUE" "$OUT" | "$MGED" -c -a nu -r "$DB" > "$LOG" 2>&1

if ! grep -q '^ATTACHING tkobol (Tk Obol display endpoint)$' "$LOG" ||
    ! grep -qx 'tk-gl' "$LOG" || ! grep -qx '320 240' "$LOG"; then
    echo "MGED did not attach and resize the direct Tk Obol endpoint" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if grep -q 'attach(tkobol): BAD\|release: .* not found' "$LOG" ||
    [ ! -s "$OUT" ]; then
    echo "MGED direct Tk Obol endpoint lifecycle failed" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

before_center=$(awk '/^OBOL_CENTER_BEFORE$/ { getline; print; exit }' "$LOG")
after_center=$(awk '/^OBOL_CENTER_AFTER$/ { getline; print; exit }' "$LOG")
if [ "$before_center" != "-0 -0 -0" ] ||
    [ "$after_center" != "-200 0 0" ]; then
    echo "MGED Tk Obol endpoint did not apply the snapped primary-pan gesture" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

before_aet=$(awk '/^OBOL_AET_BEFORE$/ { getline; print; exit }' "$LOG")
after_aet=$(awk '/^OBOL_AET_AFTER$/ { getline; print; exit }' "$LOG")
if [ "$before_aet" != "270 90 0" ] ||
    [ "$after_aet" != "180 77.5 90" ]; then
    echo "MGED Tk Obol active-mode motion was missing or applied twice" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

exit 0
