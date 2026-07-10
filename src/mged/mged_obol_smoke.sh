#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
    echo "Usage: mged_obol_smoke.sh <mged> <db> <workdir>" 1>&2
    exit 1
fi

MGED="$1"
DB="$2"
WORKDIR="$3"
OUT="${WORKDIR}/mged_obol_smoke.png"
LOG="${WORKDIR}/mged_obol_smoke.log"

rm -f "$OUT" "$LOG"

printf 'dm type\ndm size\ndraw all.g\nautoview\nrefresh\nscreengrab %s\nquit\n' "$OUT" \
    | "$MGED" -c -a obol -r "$DB" > "$LOG" 2>&1

if ! grep -qx "obol" "$LOG"; then
    echo "MGED Obol smoke did not attach an obol display manager" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -Eq '^[1-9][0-9]* [1-9][0-9]*$' "$LOG"; then
    echo "MGED Obol smoke did not report a nonzero display size" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if [ ! -s "$OUT" ]; then
    echo "MGED Obol smoke did not create a PNG: $OUT" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

bytes=`wc -c < "$OUT"`
if [ "$bytes" -le 100 ]; then
    echo "MGED Obol smoke PNG is unexpectedly small: $bytes bytes" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

exit 0
