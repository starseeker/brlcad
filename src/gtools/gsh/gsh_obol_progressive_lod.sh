#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
    echo "Usage: gsh_obol_progressive_lod.sh <gsh> <db> <workdir>" 1>&2
    exit 1
fi

GSH="$1"
DB="$2"
WORKDIR="$3"
PYTHON="${PYTHON:-python3.11}"
TMPDB="${WORKDIR}/gsh_obol_progressive_lod.g"
CACHE="${WORKDIR}/gsh_obol_progressive_lod_cache"
LOG="${WORKDIR}/gsh_obol_progressive_lod.log"
FRAME0="${WORKDIR}/gsh_obol_progressive_lod_0.png"
FRAME1="${WORKDIR}/gsh_obol_progressive_lod_1.png"
FRAME2="${WORKDIR}/gsh_obol_progressive_lod_2.png"
FRAME3="${WORKDIR}/gsh_obol_progressive_lod_3.png"

rm -f "$TMPDB" "$LOG" "$FRAME0" "$FRAME1" "$FRAME2" "$FRAME3"
rm -rf "$CACHE"
mkdir -p "$CACHE"
cp "$DB" "$TMPDB"

printf 'dm attach obol
view lod cache clear all_files
tol rel 0.0002
facetize -r all.g all.bot
view lod mesh 1
view lod scale 0.8
view lod service start 4
draw -m1 all.bot
autoview
screengrab %s
view lod service status
delay 0 150000
view lod service poll 1
screengrab %s
view lod service status
delay 0 500000
view lod service poll 1
screengrab %s
view lod service status
view lod service wait 5000 8
screengrab %s
view lod service status
quit
' "$FRAME0" "$FRAME1" "$FRAME2" "$FRAME3" \
    | BU_DIR_CACHE="$CACHE" BRLOBOL_LOD_OBB_TASK_DELAY_MS="${BRLOBOL_LOD_OBB_TASK_DELAY_MS:-350}" BRLOBOL_LOD_TASK_DELAY_MS="${BRLOBOL_LOD_TASK_DELAY_MS:-700}" \
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

if grep -Eq 'last_rejected_results: [1-9][0-9]*' "$LOG"; then
    echo "gsh Obol progressive LoD rejected LoD results" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -Eq 'active_lod_aabb_proxies: [1-9][0-9]*' "$LOG"; then
    echo "gsh Obol progressive LoD did not activate an AABB proxy payload" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -Eq 'active_lod_obb_proxies: [1-9][0-9]*' "$LOG"; then
    echo "gsh Obol progressive LoD did not activate an OBB proxy payload" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -Eq 'active_lod_mesh_payloads: [1-9][0-9]*' "$LOG"; then
    echo "gsh Obol progressive LoD did not activate a mesh LoD payload" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

"$PYTHON" - "$FRAME0" "$FRAME1" "$FRAME2" "$FRAME3" <<'PY'
import struct
import sys
import zlib


def decode_png(path):
    data = open(path, "rb").read()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError("%s is not a PNG" % path)
    off = 8
    width = height = color_type = None
    raw = b""
    while off < len(data):
        n = struct.unpack(">I", data[off:off + 4])[0]
        typ = data[off + 4:off + 8]
        chunk = data[off + 8:off + 8 + n]
        off += 12 + n
        if typ == b"IHDR":
            width, height, bit_depth, color_type, _, _, interlace = struct.unpack(">IIBBBBB", chunk)
            if bit_depth != 8 or interlace != 0:
                raise RuntimeError("%s uses unsupported PNG encoding" % path)
        elif typ == b"IDAT":
            raw += chunk
    bpp = 3 if color_type == 2 else 4 if color_type == 6 else None
    if bpp is None:
        raise RuntimeError("%s uses unsupported PNG color type %s" % (path, color_type))
    pixels = zlib.decompress(raw)
    rows = []
    prev = bytearray(width * bpp)
    p = 0
    lit = 0
    for y in range(height):
        filt = pixels[p]
        p += 1
        row = bytearray(pixels[p:p + width * bpp])
        p += width * bpp
        for i in range(len(row)):
            left = row[i - bpp] if i >= bpp else 0
            up = prev[i]
            up_left = prev[i - bpp] if i >= bpp else 0
            if filt == 1:
                row[i] = (row[i] + left) & 255
            elif filt == 2:
                row[i] = (row[i] + up) & 255
            elif filt == 3:
                row[i] = (row[i] + ((left + up) // 2)) & 255
            elif filt == 4:
                pred = left + up - up_left
                pa = abs(pred - left)
                pb = abs(pred - up)
                pc = abs(pred - up_left)
                pr = left if pa <= pb and pa <= pc else up if pb <= pc else up_left
                row[i] = (row[i] + pr) & 255
            elif filt != 0:
                raise RuntimeError("%s uses unsupported PNG filter %s" % (path, filt))
        for x in range(width):
            r, g, b = row[x * bpp:x * bpp + 3]
            if r > 32 or g > 32 or b > 32:
                lit += 1
        rows.append(bytes(row))
        prev = row
    return width, height, bpp, rows, lit


imgs = [decode_png(p) for p in sys.argv[1:]]
for path, img in zip(sys.argv[1:], imgs):
    if img[4] < 20:
        raise RuntimeError("%s is too dark for progressive LoD validation: lit=%d" % (path, img[4]))
if any(imgs[0][0:3] != img[0:3] for img in imgs[1:]):
    raise RuntimeError("progressive LoD frames have inconsistent dimensions")
diff01 = sum(a != b for ra, rb in zip(imgs[0][3], imgs[1][3]) for a, b in zip(ra, rb))
diff12 = sum(a != b for ra, rb in zip(imgs[1][3], imgs[2][3]) for a, b in zip(ra, rb))
diff23 = sum(a != b for ra, rb in zip(imgs[2][3], imgs[3][3]) for a, b in zip(ra, rb))
diff03 = sum(a != b for ra, rb in zip(imgs[0][3], imgs[3][3]) for a, b in zip(ra, rb))
if diff03 <= 0 or (diff01 <= 0 and diff12 <= 0 and diff23 <= 0):
    raise RuntimeError("progressive LoD frames did not change: diff01=%d diff12=%d diff23=%d diff03=%d" %
                       (diff01, diff12, diff23, diff03))
print("progressive_lod_diff diff01=%d diff12=%d diff23=%d diff03=%d lit=%d,%d,%d,%d" %
      (diff01, diff12, diff23, diff03, imgs[0][4], imgs[1][4], imgs[2][4], imgs[3][4]))
PY

exit 0
