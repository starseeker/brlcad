#!/bin/sh
set -eu

# A non-login Git for Windows sh does not add its own core utilities to PATH.
PATH="/usr/bin:/bin:${PATH}"
export PATH

if [ "$#" -ne 6 ]; then
    echo "Usage: gsh_obol_ert_smoke.sh <gsh> <db> <workdir> <png-pix> <pixcount> <pixdiff>" 1>&2
    exit 1
fi

GSH="$1"
DB="$2"
WORKDIR="$3"
PNG_PIX="$4"
PIXCOUNT="$5"
PIXDIFF="$6"

PARTIAL0="${WORKDIR}/gsh_obol_ert_partial_0.png"
PARTIAL1="${WORKDIR}/gsh_obol_ert_partial_1.png"
PARTIAL2="${WORKDIR}/gsh_obol_ert_partial_2.png"
FINAL="${WORKDIR}/gsh_obol_ert_final.png"
LOG="${WORKDIR}/gsh_obol_ert_smoke.log"
PARTIAL0_PIX="${WORKDIR}/gsh_obol_ert_partial_0.pix"
PARTIAL1_PIX="${WORKDIR}/gsh_obol_ert_partial_1.pix"
PARTIAL2_PIX="${WORKDIR}/gsh_obol_ert_partial_2.pix"
FINAL_PIX="${WORKDIR}/gsh_obol_ert_final.pix"
DIFF_PIX="${WORKDIR}/gsh_obol_ert_diff.pix"

rm -f "$PARTIAL0" "$PARTIAL1" "$PARTIAL2" "$FINAL" "$LOG"
rm -f "$PARTIAL0_PIX" "$PARTIAL1_PIX" "$PARTIAL2_PIX" "$FINAL_PIX" "$DIFF_PIX"

# One processor plus deterministic hypersampling keeps this render in flight
# long enough to capture a bounded early frame.  zap removes the wire overlay,
# so the intermediate image measures framebuffer progress rather than geometry.
# ERT is the behavior under test here, not deferred CAD realization.  Prepare
# this small fixture eagerly so its autoview is deterministic and the bounded
# captures measure raytrace progress rather than a race with the draw worker.
printf 'draw --eager-leaf-expansion all.g
autoview
ert -P 1 -H 8
zap
delay 0 50000
screengrab %s
delay 0 100000
screengrab %s
delay 0 150000
screengrab %s
delay 2 0
screengrab %s
ert -P 1 -H 16
delay 0 100000
quit
' "$PARTIAL0" "$PARTIAL1" "$PARTIAL2" "$FINAL" | "$GSH" --new-cmds "$DB" > "$LOG" 2>&1

if ! grep -q "rt: launching endpoint framebuffer renderer (ipc=1" "$LOG"; then
    echo "gsh Obol ert smoke did not use the Obol IPC framebuffer path" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if grep -q "rt_gettrees.*FAILED" "$LOG"; then
    echo "gsh Obol ert smoke failed to prepare raytrace geometry" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -q "SHOT:" "$LOG" || ! grep -q "Frame  *0:" "$LOG"; then
    echo "gsh Obol ert smoke did not produce raytrace frame output" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

for image in "$PARTIAL0" "$PARTIAL1" "$PARTIAL2" "$FINAL"; do
    if [ ! -s "$image" ]; then
	echo "gsh Obol ert smoke did not create a PNG: $image" 1>&2
	cat "$LOG" 1>&2
	exit 1
    fi
done

if ! "$PNG_PIX" "$PARTIAL0" > "$PARTIAL0_PIX" ||
   ! "$PNG_PIX" "$PARTIAL1" > "$PARTIAL1_PIX" ||
   ! "$PNG_PIX" "$PARTIAL2" > "$PARTIAL2_PIX" ||
   ! "$PNG_PIX" "$FINAL" > "$FINAL_PIX"; then
    echo "gsh Obol ert smoke could not decode its PNG output" 1>&2
    exit 1
fi

final_bytes=$(wc -c < "$FINAL_PIX")
if [ $((final_bytes % 3)) -ne 0 ]; then
    echo "ERT frame decoder returned incomplete RGB pixels" 1>&2
    exit 1
fi
pixel_count=$((final_bytes / 3))
# The retained renderer's nominally black background is not required to be
# bytewise 0/0/0 (the OSMesa path currently presents 0/0/1).  Use a corner
# pixel from the completed frame as the actual background and count raytrace
# progress relative to it.
set -- $(od -An -tu1 -N3 "$FINAL_PIX")
background_r=$1
background_g=$2
background_b=$3
final_background=$("$PIXCOUNT" "$FINAL_PIX" |
    awk -v r="$background_r" -v g="$background_g" -v b="$background_b" \
	'$1 == r && $2 == g && $3 == b { print $4; found = 1 } END { if (!found) print 0 }')
final_lit=$((pixel_count - final_background))
if [ "$final_lit" -le 0 ]; then
    echo "completed ERT frame contains no foreground pixels" 1>&2
    exit 1
fi

partial_pix=""
partial_lit=0
for candidate in "$PARTIAL0_PIX" "$PARTIAL1_PIX" "$PARTIAL2_PIX"; do
    candidate_bytes=$(wc -c < "$candidate")
    if [ "$candidate_bytes" -ne "$final_bytes" ]; then
	echo "partial and final ERT frames have inconsistent RGB dimensions" 1>&2
	exit 1
    fi
    candidate_background=$("$PIXCOUNT" "$candidate" |
	awk -v r="$background_r" -v g="$background_g" -v b="$background_b" \
	    '$1 == r && $2 == g && $3 == b { print $4; found = 1 } END { if (!found) print 0 }')
    candidate_lit=$((pixel_count - candidate_background))
    # Streaming scanlines replace the previously presented CAD/blank pixels;
    # foreground coverage is therefore not monotonic and may temporarily
    # exceed the completed raytrace's coverage.  A useful early frame needs
    # some non-background information and must still differ from the final
    # image.
    if [ "$candidate_lit" -gt $((final_lit / 20)) ] &&
       ! cmp -s "$candidate" "$FINAL_PIX"; then
	partial_pix="$candidate"
	partial_lit="$candidate_lit"
	break
    fi
done
if [ -z "$partial_pix" ]; then
    echo "no bounded early ERT frame was meaningfully partial" 1>&2
    exit 1
fi
"$PIXDIFF" "$partial_pix" "$FINAL_PIX" > "$DIFF_PIX"
diff_black=$("$PIXCOUNT" "$DIFF_PIX" |
    awk '$1 == 0 && $2 == 0 && $3 == 0 { print $4; found = 1 } END { if (!found) print 0 }')
changed_pixels=$((pixel_count - diff_black))
if [ "$changed_pixels" -lt $((final_lit / 10)) ]; then
    echo "partial and final ERT frames changed too little: ${changed_pixels}" 1>&2
    exit 1
fi
echo "ert_progress partial_lit=${partial_lit} final_lit=${final_lit} changed_pixels=${changed_pixels}"

exit 0
