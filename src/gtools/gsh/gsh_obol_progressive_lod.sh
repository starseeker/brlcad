#!/bin/sh
set -eu

# A non-login Git for Windows sh does not add its own core utilities to PATH.
PATH="/usr/bin:/bin:${PATH}"
export PATH

if [ "$#" -ne 6 ]; then
    echo "Usage: gsh_obol_progressive_lod.sh <gsh> <db> <workdir> <png-pix> <pixcount> <pixdiff>" 1>&2
    exit 1
fi

GSH="$1"
DB="$2"
WORKDIR="$3"
PNG_PIX="$4"
PIXCOUNT="$5"
PIXDIFF="$6"

if [ ! -d "$WORKDIR" ]; then
    echo "gsh Obol progressive LoD work directory does not exist: $WORKDIR" 1>&2
    exit 1
fi
cd "$WORKDIR"

# Keep generated paths relative after entering the work directory.  Git for
# Windows core utilities do not consistently handle drive-letter paths with
# recursive operations such as mkdir -p.
TMPDB="gsh_obol_progressive_lod.g"
CACHE="gsh_obol_progressive_lod_cache"
LOG="gsh_obol_progressive_lod.log"
FRAME0="gsh_obol_progressive_lod_0.png"
FRAME1="gsh_obol_progressive_lod_1.png"
FRAME2="gsh_obol_progressive_lod_2.png"
FRAME3="gsh_obol_progressive_lod_3.png"
PIX0="gsh_obol_progressive_lod_0.pix"
PIX1="gsh_obol_progressive_lod_1.pix"
PIX2="gsh_obol_progressive_lod_2.pix"
PIX3="gsh_obol_progressive_lod_3.pix"
DIFF01="gsh_obol_progressive_lod_01.diff.pix"
DIFF12="gsh_obol_progressive_lod_12.diff.pix"
DIFF23="gsh_obol_progressive_lod_23.diff.pix"

rm -f "$TMPDB" "$LOG" "$FRAME0" "$FRAME1" "$FRAME2" "$FRAME3"
rm -f "$PIX0" "$PIX1" "$PIX2" "$PIX3" "$DIFF01" "$DIFF12" "$DIFF23"
rm -rf "$CACHE"
mkdir -p "$CACHE"
cp "$DB" "$TMPDB"

printf 'view lod cache clear all_files
tol rel 0.0002
facetize -r all.g all.bot
view lod mesh 1
view lod scale 0.8
view lod service start 4
draw -m1 all.bot
autoview
screengrab %s
view lod service status
view lod service poll 64
delay 4 0
view lod service poll 64
screengrab %s
view lod service status
delay 4 0
view lod service poll 64
screengrab %s
view lod service status
view lod service wait 60000 64
screengrab %s
view lod service status
quit
' "$FRAME0" "$FRAME1" "$FRAME2" "$FRAME3" \
    | BU_DIR_CACHE="$CACHE" BOBOL_LOD_TASK_DELAY_MS="${BOBOL_LOD_TASK_DELAY_MS:-6000}" \
    "$GSH" --new-cmds "$TMPDB" > "$LOG" 2>&1

for img in "$FRAME0" "$FRAME1" "$FRAME2" "$FRAME3"; do
    if [ ! -s "$img" ]; then
	echo "gsh Obol progressive LoD did not create PNG: $img" 1>&2
	cat "$LOG" 1>&2
	exit 1
    fi
done

if ! grep -Eq 'last_visited_meshes: [1-9][0-9]*' "$LOG"; then
    echo "gsh Obol progressive LoD did not visit mesh LoD candidates" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -Eq 'last_submitted_tasks: [1-9][0-9]*' "$LOG"; then
    echo "gsh Obol progressive LoD did not submit LoD tasks" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -Eq 'last_applied_results: [1-9][0-9]*' "$LOG"; then
    echo "gsh Obol progressive LoD did not apply LoD results" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -Eq '^wait_timed_out=0 .*pending_submissions=0 pending_results=0 progressive_pending=0 queued=0 in_flight=0 pending=0 delayed=0$' "$LOG"; then
    echo "gsh Obol progressive LoD did not reach a terminal idle state" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

final_status_value()
{
    awk -v key="$1:" '$1 == key { value = $2; found = 1 } END {
	if (!found) exit 1
	print value
    }' "$LOG"
}

# Demand epochs may legitimately overtake an in-flight task.  Such a result is
# reported as rejected in that intermediate batch and is resubmitted for the
# current demand.  Audit the settled batch instead: it must contain neither a
# current-demand rejection nor an unmatched route, and every fixture leaf must
# own a live mesh payload.
final_rejected=$(final_status_value last_rejected_results) || {
    echo "gsh Obol progressive LoD did not report terminal rejection status" 1>&2
    cat "$LOG" 1>&2
    exit 1
}
final_unmatched=$(final_status_value last_unmatched_results) || {
    echo "gsh Obol progressive LoD did not report terminal routing status" 1>&2
    cat "$LOG" 1>&2
    exit 1
}
final_payloads=$(final_status_value active_lod_mesh_payloads) || {
    echo "gsh Obol progressive LoD did not report terminal payload status" 1>&2
    cat "$LOG" 1>&2
    exit 1
}
if [ "$final_rejected" -ne 0 ] || [ "$final_unmatched" -ne 0 ]; then
    echo "gsh Obol progressive LoD retained rejected or unmatched terminal results: rejected=${final_rejected} unmatched=${final_unmatched}" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi
if [ "$final_payloads" -ne 6 ]; then
    echo "gsh Obol progressive LoD did not activate all fixture payloads: active=${final_payloads} expected=6" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! "$PNG_PIX" "$FRAME0" > "$PIX0" ||
   ! "$PNG_PIX" "$FRAME1" > "$PIX1" ||
   ! "$PNG_PIX" "$FRAME2" > "$PIX2" ||
   ! "$PNG_PIX" "$FRAME3" > "$PIX3"; then
    echo "gsh Obol progressive LoD could not decode its PNG output" 1>&2
    exit 1
fi

bytes0=$(wc -c < "$PIX0")
for pix in "$PIX1" "$PIX2" "$PIX3"; do
    bytes=$(wc -c < "$pix")
    if [ "$bytes" -ne "$bytes0" ]; then
	echo "progressive LoD frames have inconsistent RGB dimensions" 1>&2
	exit 1
    fi
done
if [ $((bytes0 % 3)) -ne 0 ]; then
    echo "progressive LoD decoder returned incomplete RGB pixels" 1>&2
    exit 1
fi
pixel_count=$((bytes0 / 3))

lit_pixels()
{
    "$PIXCOUNT" "$1" |
	awk '$1 > 32 || $2 > 32 || $3 > 32 { lit += $4 } END { print lit + 0 }'
}

lit0=$(lit_pixels "$PIX0")
lit1=$(lit_pixels "$PIX1")
lit2=$(lit_pixels "$PIX2")
lit3=$(lit_pixels "$PIX3")
for lit in "$lit1" "$lit2" "$lit3"; do
    if [ "$lit" -lt 20 ]; then
	echo "post-poll progressive LoD frame is too dark: lit=${lit}" 1>&2
	exit 1
    fi
done

changed_pixels()
{
    "$PIXDIFF" "$1" "$2" > "$3"
    black=$("$PIXCOUNT" "$3" |
	awk '$1 == 0 && $2 == 0 && $3 == 0 { print $4; found = 1 } END { if (!found) print 0 }')
    echo $((pixel_count - black))
}

diff01=$(changed_pixels "$PIX0" "$PIX1" "$DIFF01")
diff12=$(changed_pixels "$PIX1" "$PIX2" "$DIFF12")
diff23=$(changed_pixels "$PIX2" "$PIX3" "$DIFF23")
if [ "$diff01" -le 0 ] && [ "$diff12" -le 0 ] && [ "$diff23" -le 0 ]; then
    echo "progressive LoD frames did not change: diff01=${diff01} diff12=${diff12} diff23=${diff23}" 1>&2
    exit 1
fi
echo "progressive_lod_diff diff01=${diff01} diff12=${diff12} diff23=${diff23} lit=${lit0},${lit1},${lit2},${lit3}"

exit 0
