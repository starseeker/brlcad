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

printf 'regdebug\nset\ndm open --host headless --renderer sw\ndm host\ndm size 512 512\ndm size\ndm bg 0 0 32 64 0 0\ndraw all.g\nautoview\nrefresh\nadc draw 1\nputs OBOL_ADC_ON\ndm get view.faceplate.adc.visible\nadc draw 0\nputs OBOL_ADC_OFF\ndm get view.faceplate.adc.visible\ndm set view.faceplate.center_dot.color 0/1/0\nrefresh\nputs OBOL_CENTER_GREEN\ndm get view.faceplate.center_dot.color\nrset cs center_dot 255 0 0\nputs OBOL_CENTER_RED\ndm get view.faceplate.center_dot.color\nrefresh\nputs OBOL_CENTER_AFTER_REFRESH\ndm get view.faceplate.center_dot.color\ndm set view.faceplate.model_axes.color 0/1/0\nrefresh\nputs OBOL_AXES_GREEN\ndm get view.faceplate.model_axes.color\nrset cs model_axes 255 0 0\nputs OBOL_AXES_RED\ndm get view.faceplate.model_axes.color\nrefresh\nputs OBOL_AXES_AFTER_REFRESH\ndm get view.faceplate.model_axes.color\ndm set view.faceplate.adc.line_color 0/1/0\ndm set view.faceplate.adc.tick_color 0/1/0\nrefresh\nputs OBOL_ADC_STYLE_GREEN\ndm get view.faceplate.adc.line_color\ndm get view.faceplate.adc.tick_color\nrset cs adc_line 255 0 0\nrset cs adc_tick 0 0 255\nputs OBOL_ADC_STYLE_SCHEME\ndm get view.faceplate.adc.line_color\ndm get view.faceplate.adc.tick_color\npostscript -c ObolSmoke -l 3 -z %s\nplot %s\nscreengrab %s\nquit\n' "$PS_OUT" "$PLOT_OUT" "$OUT" \
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

if grep -Eq '(^|[[:space:]])cache=' "$LOG"; then
    echo "MGED still exposes retired backend-cache policy" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -Fq 'invalid command name "regdebug"' "$LOG"; then
    echo "MGED still accepts retired display-manager debug control" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -A1 -Fx 'OBOL_ADC_ON' "$LOG" | grep -Fx 'true' >/dev/null ||
    ! grep -A1 -Fx 'OBOL_ADC_OFF' "$LOG" | grep -Fx 'false' >/dev/null; then
    echo "MGED adc command did not update the Obol visibility policy" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -A1 -Fx 'OBOL_CENTER_GREEN' "$LOG" | grep -Fx '0/1/0' >/dev/null; then
    echo "MGED Obol endpoint faceplate override did not survive refresh" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -A1 -Fx 'OBOL_CENTER_RED' "$LOG" | grep -Fx '1/0/0' >/dev/null ||
    ! grep -A1 -Fx 'OBOL_CENTER_AFTER_REFRESH' "$LOG" | grep -Fx '1/0/0' >/dev/null; then
    echo "MGED color-scheme update did not replace the endpoint faceplate override" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -A1 -Fx 'OBOL_AXES_GREEN' "$LOG" | grep -Fx '0/1/0' >/dev/null; then
    echo "MGED endpoint axes override did not survive refresh" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -A1 -Fx 'OBOL_AXES_RED' "$LOG" | grep -Fx '1/0/0' >/dev/null ||
    ! grep -A1 -Fx 'OBOL_AXES_AFTER_REFRESH' "$LOG" | grep -Fx '1/0/0' >/dev/null; then
    echo "MGED color-scheme update did not republish endpoint axes colors" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if [ "$(grep -A2 -Fx 'OBOL_ADC_STYLE_GREEN' "$LOG" | grep -Fx '0/1/0' | wc -l)" -ne 2 ]; then
    echo "MGED endpoint ADC style override did not survive refresh" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -A2 -Fx 'OBOL_ADC_STYLE_SCHEME' "$LOG" | grep -Fx '1/0/0' >/dev/null ||
    ! grep -A2 -Fx 'OBOL_ADC_STYLE_SCHEME' "$LOG" | grep -Fx '0/0/1' >/dev/null; then
    echo "MGED color-scheme update did not republish endpoint ADC style" 1>&2
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
