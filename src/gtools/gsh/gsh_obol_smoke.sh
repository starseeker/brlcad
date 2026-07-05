#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
    echo "Usage: gsh_obol_smoke.sh <gsh> <db> <workdir>" 1>&2
    exit 1
fi

GSH="$1"
DB="$2"
WORKDIR="$3"
OUT="${WORKDIR}/gsh_obol_smoke.png"
LOG="${WORKDIR}/gsh_obol_smoke.log"

rm -f "$OUT" "$LOG"

printf 'dm attach obol\ndraw all.g\nautoview\nscreengrab %s\nquit\n' "$OUT" \
    | "$GSH" --new-cmds "$DB" > "$LOG" 2>&1

if [ ! -s "$OUT" ]; then
    echo "gsh Obol smoke did not create a PNG: $OUT" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

bytes=`wc -c < "$OUT"`
if [ "$bytes" -le 100 ]; then
    echo "gsh Obol smoke PNG is unexpectedly small: $bytes bytes" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

exit 0
