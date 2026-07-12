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

rm -f "$OUT" "$LOG"

printf 'attach -t 1 -n .mged_direct tkobol
dm host
dm size 320 240
dm size
draw all.g
autoview
refresh
screengrab %s
release .mged_direct
quit
' "$OUT" | "$MGED" -c -a nu -r "$DB" > "$LOG" 2>&1

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

exit 0
