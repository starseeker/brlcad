#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
    echo "Usage: mged_obol_ert_smoke.sh <mged> <db> <workdir>" 1>&2
    exit 1
fi

MGED="$1"
DB="$2"
WORKDIR="$3"
PYTHON="${PYTHON:-python3.11}"
OUT="${WORKDIR}/mged_obol_ert_smoke.png"
LOG="${WORKDIR}/mged_obol_ert_smoke.log"

rm -f "$OUT" "$LOG"

printf 'dm type
draw all.g
autoview
ert
delay 5 0
refresh
screengrab %s
quit
' "$OUT" | "$MGED" -c -a obol -r "$DB" > "$LOG" 2>&1

if ! grep -qx "obol" "$LOG"; then
    echo "MGED Obol ert smoke did not attach an obol display manager" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if ! grep -q "Raytrace complete" "$LOG"; then
    echo "MGED Obol ert smoke did not complete the rt subprocess" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if grep -q "rt_gettrees.*FAILED" "$LOG"; then
    echo "MGED Obol ert smoke failed to prepare raytrace geometry" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

if [ ! -s "$OUT" ]; then
    echo "MGED Obol ert smoke did not create a PNG: $OUT" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

"$PYTHON" - "$OUT" <<'PY'
import struct
import sys
import zlib


def paeth(a, b, c):
    p = a + b - c
    pa = abs(p - a)
    pb = abs(p - b)
    pc = abs(p - c)
    if pa <= pb and pa <= pc:
        return a
    return b if pb <= pc else c


path = sys.argv[1]
data = open(path, "rb").read()
if data[:8] != b"\x89PNG\r\n\x1a\n":
    raise SystemExit("%s is not a PNG" % path)

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
            raise SystemExit("%s uses unsupported PNG encoding" % path)
    elif typ == b"IDAT":
        raw += chunk

bpp = 3 if color_type == 2 else 4 if color_type == 6 else None
if bpp is None:
    raise SystemExit("%s uses unsupported PNG color type %s" % (path, color_type))

pixels = zlib.decompress(raw)
prev = bytearray(width * bpp)
p = 0
unique = set()
nonblack = 0
for _ in range(height):
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
            row[i] = (row[i] + paeth(left, up, up_left)) & 255
        elif filt != 0:
            raise SystemExit("%s uses unsupported PNG filter %s" % (path, filt))
    for i in range(0, len(row), bpp):
        rgb = tuple(row[i:i + 3])
        unique.add(rgb)
        if rgb != (0, 0, 0):
            nonblack += 1
    prev = row

if len(unique) <= 1 or nonblack <= 0:
    raise SystemExit("%s has no visible framebuffer content" % path)
PY

exit 0
