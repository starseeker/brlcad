#!/bin/sh
set -eu

# A non-login Git for Windows sh does not add its own core utilities to PATH.
# Absolute /usr/bin and /bin are mapped by MSYS and are harmless on Unix.
PATH="/usr/bin:/bin:${PATH}"
export PATH

if [ "$#" -ne 3 ]; then
    echo "Usage: mged_obol_smoke.sh <mged> <db> <workdir>" 1>&2
    exit 1
fi

MGED="$1"
DB="$2"
WORKDIR="$3"
OUT="${WORKDIR}/mged_obol_smoke.png"
PIX_OUT="${WORKDIR}/mged_obol_smoke.pix"
RT_OUT="${WORKDIR}/mged_obol_smoke_rt.png"
RT_PIX="${WORKDIR}/mged_obol_smoke_rt.pix"
RT_MOVED_OUT="${WORKDIR}/mged_obol_smoke_rt_moved.png"
RT_MOVED_PIX="${WORKDIR}/mged_obol_smoke_rt_moved.pix"
PS_OUT="${WORKDIR}/mged_obol_smoke.ps"
PLOT_OUT="${WORKDIR}/mged_obol_smoke.plot3"
LOG="${WORKDIR}/mged_obol_smoke.log"
MGED_DIR=${MGED%/*}
if [ "$MGED_DIR" = "$MGED" ]; then
    MGED_DIR=.
fi
PNG_PIX="${MGED_DIR}/png-pix"

rm -f "$OUT" "$PIX_OUT" "$RT_OUT" "$RT_PIX" "$RT_MOVED_OUT" \
    "$RT_MOVED_PIX" "$PS_OUT" "$PLOT_OUT" "$LOG"

HELP_OUTPUT=$("$MGED" --help 2>&1 || true)
if ! printf '%s\n' "$HELP_OUTPUT" | grep -q -- '--host'; then
    echo "MGED help does not advertise the Obol host option" 1>&2
    printf '%s\n' "$HELP_OUTPUT" 1>&2
    exit 1
fi
if printf '%s\n' "$HELP_OUTPUT" | grep -q -- '--dm-type'; then
    echo "MGED help still advertises the retired display-manager type option" 1>&2
    printf '%s\n' "$HELP_OUTPUT" 1>&2
    exit 1
fi
if printf '%s\n' "$HELP_OUTPUT" | grep -q -- '--eye-sep-dist'; then
    echo "MGED help still advertises unsupported stereo policy" 1>&2
    printf '%s\n' "$HELP_OUTPUT" 1>&2
    exit 1
fi
if "$MGED" --dm-type tkobol --help >/dev/null 2>&1; then
    echo "MGED still accepts the retired display-manager type option" 1>&2
    exit 1
fi

printf 'regdebug\nset\ndm open --host headless --renderer sw\ndm host\ndm size 512 512\ndm size\ndm bg 0 0 32 64 0 0\ndraw all.g\nautoview\nrefresh\nset perspective 45\nputs OBOL_PERSPECTIVE_ON\ndm get view.perspective\nset perspective -1\nputs OBOL_PERSPECTIVE_OFF\ndm get view.perspective\nadc draw 1\nputs OBOL_ADC_ON\ndm get view.faceplate.adc.visible\nadc draw 0\nputs OBOL_ADC_OFF\ndm get view.faceplate.adc.visible\ndm set view.faceplate.center_dot.color 0/1/0\nrefresh\nputs OBOL_CENTER_GREEN\ndm get view.faceplate.center_dot.color\nrset cs center_dot 255 0 0\nputs OBOL_CENTER_RED\ndm get view.faceplate.center_dot.color\nrefresh\nputs OBOL_CENTER_AFTER_REFRESH\ndm get view.faceplate.center_dot.color\ndm set view.faceplate.model_axes.color 0/1/0\nrefresh\nputs OBOL_AXES_GREEN\ndm get view.faceplate.model_axes.color\nrset cs model_axes 255 0 0\nputs OBOL_AXES_RED\ndm get view.faceplate.model_axes.color\nrefresh\nputs OBOL_AXES_AFTER_REFRESH\ndm get view.faceplate.model_axes.color\ndm set view.faceplate.adc.line_color 0/1/0\ndm set view.faceplate.adc.tick_color 0/1/0\nrefresh\nputs OBOL_ADC_STYLE_GREEN\ndm get view.faceplate.adc.line_color\ndm get view.faceplate.adc.tick_color\nrset cs adc_line 255 0 0\nrset cs adc_tick 0 0 255\nputs OBOL_ADC_STYLE_SCHEME\ndm get view.faceplate.adc.line_color\ndm get view.faceplate.adc.tick_color\npostscript -c ObolSmoke -l 3 -z %s\nplot %s\nscreengrab %s\ndm bg 0 0 0\ndm set renderer.headlight.color 1/0.9/0.75\ndm set renderer.headlight.intensity 0.8\ndm set render.rt.samples 1\ndm set render.rt.preview_scale 4\ndm set render.rt.frame_budget_ms 1\ndm renderer rt\nrefresh\nscreengrab %s\nputs OBOL_RT_REVISION_0\ndm get render.rt.geometry_revision\ndm get render.rt.presentation_revision\ndm set render.rt.samples 8\nae 35 25\nrefresh\nscreengrab %s\nputs OBOL_RT_REVISION_1\ndm get render.rt.geometry_revision\ndm get render.rt.presentation_revision\ndm set render.rt.samples 1\ndm renderer sw\nquit\n' "$PS_OUT" "$PLOT_OUT" "$OUT" "$RT_OUT" "$RT_MOVED_OUT" \
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

if grep -Eq '(^|[[:space:]])eye_sep_dist=' "$LOG"; then
    echo "MGED still exposes unsupported stereo state" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -Fq 'invalid command name "regdebug"' "$LOG"; then
    echo "MGED still accepts retired display-manager debug control" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -A1 -Fx 'OBOL_PERSPECTIVE_ON' "$LOG" | grep -Fx '45' >/dev/null ||
    ! grep -A1 -Fx 'OBOL_PERSPECTIVE_OFF' "$LOG" | grep -Fx '0' >/dev/null; then
    echo "MGED perspective settings did not use the endpoint policy" 1>&2
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

if [ ! -s "$RT_OUT" ] || [ ! -s "$RT_MOVED_OUT" ] ||
    ! "$PNG_PIX" "$RT_OUT" > "$RT_PIX" ||
    ! "$PNG_PIX" "$RT_MOVED_OUT" > "$RT_MOVED_PIX"; then
    echo "MGED retained RT smoke could not inspect both captures" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

rt_lit=`od -An -tu1 "$RT_PIX" | awk '{for (i=1; i<=NF; i++) if ($i != 0) n++} END {print n+0}'`
rt_moved_lit=`od -An -tu1 "$RT_MOVED_PIX" | awk '{for (i=1; i<=NF; i++) if ($i != 0) n++} END {print n+0}'`
if [ "$rt_lit" -le 100 ] || [ "$rt_moved_lit" -le 100 ] ||
    cmp -s "$RT_PIX" "$RT_MOVED_PIX"; then
    echo "MGED retained RT captures did not show camera-dependent geometry" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

rt_revisions_0=`sed -n '/^OBOL_RT_REVISION_0$/{n;p;n;p;}' "$LOG"`
rt_revisions_1=`sed -n '/^OBOL_RT_REVISION_1$/{n;p;n;p;}' "$LOG"`
rt_geometry_0=`printf '%s\n' "$rt_revisions_0" | sed -n '1p'`
rt_presentation_0=`printf '%s\n' "$rt_revisions_0" | sed -n '2p'`
rt_geometry_1=`printf '%s\n' "$rt_revisions_1" | sed -n '1p'`
rt_presentation_1=`printf '%s\n' "$rt_revisions_1" | sed -n '2p'`
for revision in "$rt_geometry_0" "$rt_presentation_0" \
    "$rt_geometry_1" "$rt_presentation_1"; do
    case "$revision" in
    ''|*[!0-9]*)
	echo "MGED retained RT smoke did not report numeric renderer revisions" 1>&2
	cat "$LOG" 1>&2
	exit 1
	;;
    esac
done
if [ "$rt_geometry_0" -le 0 ] || [ "$rt_geometry_0" -ne "$rt_geometry_1" ] ||
    [ "$rt_presentation_1" -le "$rt_presentation_0" ]; then
    echo "MGED retained RT camera motion rebuilt geometry or missed presentation revision" 1>&2
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
