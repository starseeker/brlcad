#!/bin/sh
set -eu

PATH="/usr/bin:/bin:${PATH}"
export PATH

if [ "$#" -ne 3 ]; then
    echo "Usage: mged_obol_edit_restore.sh <mged> <db> <workdir>" 1>&2
    exit 1
fi

MGED="$1"
DB="$2"
WORKDIR="$3"
PYTHON="${PYTHON:-}"
if [ -z "$PYTHON" ]; then
    for cand in python3.11 python3.14 python3 python; do
	if command -v "$cand" >/dev/null 2>&1; then
	    PYTHON="$cand"
	    break
	fi
    done
fi
if [ -z "$PYTHON" ]; then
    echo "MGED Obol edit/restore image test requires Python" 1>&2
    exit 125
fi

TMPDB="${WORKDIR}/mged_obol_edit_restore.g"
CACHE="${WORKDIR}/mged_obol_edit_restore_cache"
LOG="${WORKDIR}/mged_obol_edit_restore.log"
WARMUP="${WORKDIR}/mged_obol_edit_restore_warmup.png"
BEFORE="${WORKDIR}/mged_obol_edit_restore_before.png"
EDIT="${WORKDIR}/mged_obol_edit_restore_edit.png"
RESTORE0="${WORKDIR}/mged_obol_edit_restore_restore0.png"
RESTORE1="${WORKDIR}/mged_obol_edit_restore_restore1.png"
MGED_DIR=${MGED%/*}
if [ "$MGED_DIR" = "$MGED" ]; then
    MGED_DIR=.
fi
PNG_PIX="${MGED_DIR}/png-pix"

rm -f "$TMPDB" "$LOG" "$WARMUP" "$BEFORE" "$EDIT" "$RESTORE0" "$RESTORE1"
cmake -E rm -rf "$CACHE"
cmake -E make_directory "$CACHE"
cp "$DB" "$TMPDB"

# Disable both LoD channels so the fixture isolates retained presentation
# identity.  A one-occurrence solid-edit promotion must hide and restore the
# shaded instance before the very next render, without rediscovery or a blank
# intermediate frame.
printf 'dm open --host headless --renderer sw
dm size 512 512
view lod mesh 0
view lod csg 0
draw --eager-leaf-expansion -m1 all.g
autoview
screengrab %s
screengrab %s
sed all.g/tor.r/tor
puts EDIT_STATE
status state
screengrab %s
reject
puts RESTORE_STATE
status state
screengrab %s
screengrab %s
quit
' "$WARMUP" "$BEFORE" "$EDIT" "$RESTORE0" "$RESTORE1" \
    | BU_DIR_CACHE="$CACHE" "$MGED" -c -a nu "$TMPDB" > "$LOG" 2>&1

if ! awk '/^EDIT_STATE$/ {getline; if ($0 == "SOL EDIT") found=1} END {exit !found}' "$LOG" ||
   ! awk '/^RESTORE_STATE$/ {getline; if ($0 == "VIEWING") found=1} END {exit !found}' "$LOG"; then
    echo "MGED Obol edit/restore did not traverse SOL EDIT -> VIEWING" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

for img in "$BEFORE" "$EDIT" "$RESTORE0" "$RESTORE1"; do
    if [ ! -s "$img" ]; then
	echo "MGED Obol edit/restore did not create PNG: $img" 1>&2
	cat "$LOG" 1>&2
	exit 1
    fi
    "$PNG_PIX" "$img" > "${img}.pix"
done

"$PYTHON" - "$BEFORE.pix" "$EDIT.pix" "$RESTORE0.pix" "$RESTORE1.pix" <<'PY'
import sys

width = 512
height = 512
components = 3


def scene_crop(path):
    data = open(path, "rb").read()
    expected = width * height * components
    if len(data) != expected:
        raise RuntimeError("%s has %d bytes, expected %d" %
                           (path, len(data), expected))
    # Exclude MGED's faceplate text, border, FPS, and right-side progress HUD.
    # The retained test geometry and edited torus occupy this central region.
    rows = []
    for y in range(120, 390):
        first = (y * width + 120) * components
        last = (y * width + 390) * components
        rows.append(data[first:last])
    return b"".join(rows)


before, edit, restore0, restore1 = [scene_crop(path) for path in sys.argv[1:]]
edit_delta = sum(a != b for a, b in zip(before, edit))
restore0_delta = sum(a != b for a, b in zip(before, restore0))
restore1_delta = sum(a != b for a, b in zip(before, restore1))
restore_stability = sum(a != b for a, b in zip(restore0, restore1))
if edit_delta < 100:
    raise RuntimeError("solid edit did not visibly replace its retained occurrence: delta=%d" %
                       edit_delta)
if restore0_delta or restore1_delta or restore_stability:
    raise RuntimeError("solid edit restore was not exact on its first stable frame: "
                       "restore0=%d restore1=%d stability=%d" %
                       (restore0_delta, restore1_delta, restore_stability))
print("edit_restore_diff edit=%d restore0=%d restore1=%d stability=%d" %
      (edit_delta, restore0_delta, restore1_delta, restore_stability))
PY

rm -f "$TMPDB" "$WARMUP" "$BEFORE" "$EDIT" "$RESTORE0" "$RESTORE1" \
    "$BEFORE.pix" "$EDIT.pix" "$RESTORE0.pix" "$RESTORE1.pix"
cmake -E rm -rf "$CACHE"
