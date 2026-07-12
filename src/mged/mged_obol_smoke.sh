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
PIX_OUT="${WORKDIR}/mged_obol_smoke.pix"
PS_OUT="${WORKDIR}/mged_obol_smoke.ps"
PLOT_OUT="${WORKDIR}/mged_obol_smoke.plot3"
LOG="${WORKDIR}/mged_obol_smoke.log"
PNG_PIX="`dirname "$MGED"`/png-pix"

rm -f "$OUT" "$PIX_OUT" "$PS_OUT" "$PLOT_OUT" "$LOG"

printf 'dm open --host headless --renderer sw\ndm host\ndm size 512 512\ndm size\ndm bg 0 0 32 64 0 0\ndraw all.g\nautoview\nrefresh\npostscript -c ObolSmoke -l 3 -z %s\nplot %s\nscreengrab %s\nquit\n' "$PS_OUT" "$PLOT_OUT" "$OUT" \
    | "$MGED" -c -a nu -r "$DB" > "$LOG" 2>&1

if ! grep -qx "headless" "$LOG"; then
    echo "MGED Obol smoke did not open a headless endpoint host" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -qx '512 512' "$LOG"; then
    echo "MGED Obol smoke did not preserve the requested endpoint size" 1>&2
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

if [ ! -x "$PNG_PIX" ] || ! "$PNG_PIX" "$OUT" > "$PIX_OUT"; then
    echo "MGED Obol smoke could not inspect the captured PNG" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

# PIX and endpoint captures both use bottom-left row order.  These samples
# avoid the faceplate border and prove the requested blue-bottom/red-top
# gradient was not discarded or vertically inverted by screengrab.
set -- `od -An -tu1 -N3 -j 155100 "$PIX_OUT"`
bottom_r=$1
bottom_b=$3
set -- `od -An -tu1 -N3 -j 615900 "$PIX_OUT"`
top_r=$1
top_b=$3
if [ "$bottom_b" -le "$bottom_r" ] || [ "$top_r" -le "$top_b" ]; then
    echo "MGED Obol capture lost or inverted its background gradient" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if [ ! -s "$PS_OUT" ] || ! grep -q '^%%Creator: ObolSmoke$' "$PS_OUT" ||
    ! grep -q '^3 setlinewidth$' "$PS_OUT" || ! grep -q ' lineto stroke$' "$PS_OUT"; then
    echo "MGED retained PostScript export is incomplete: $PS_OUT" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if [ ! -s "$PLOT_OUT" ]; then
    echo "MGED retained plot3 export is incomplete: $PLOT_OUT" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

exit 0
