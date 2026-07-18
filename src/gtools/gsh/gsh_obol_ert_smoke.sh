#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
    echo "Usage: gsh_obol_ert_smoke.sh <gsh> <db> <workdir>" 1>&2
    exit 1
fi

GSH="$1"
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
    echo "No usable Python interpreter found for PNG validation" 1>&2
    exit 1
fi

PARTIAL="${WORKDIR}/gsh_obol_ert_partial.png"
FINAL="${WORKDIR}/gsh_obol_ert_final.png"
LOG="${WORKDIR}/gsh_obol_ert_smoke.log"

rm -f "$PARTIAL" "$FINAL" "$LOG"

# One processor plus deterministic hypersampling keeps this render in flight
# long enough to capture a bounded early frame.  zap removes the wire overlay,
# so the intermediate image measures framebuffer progress rather than geometry.
printf 'draw all.g
autoview
ert -P 1 -H 8
zap
delay 1 0
screengrab %s
delay 5 0
screengrab %s
ert -P 1 -H 16
delay 0 100000
quit
' "$PARTIAL" "$FINAL" | "$GSH" --new-cmds "$DB" > "$LOG" 2>&1

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

for image in "$PARTIAL" "$FINAL"; do
    if [ ! -s "$image" ]; then
	echo "gsh Obol ert smoke did not create a PNG: $image" 1>&2
	cat "$LOG" 1>&2
	exit 1
    fi
done

"$PYTHON" - "$PARTIAL" "$FINAL" <<'PY'
import struct
import sys
import zlib


def paeth(a, b, c):
    p = a + b - c
    pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
    return a if pa <= pb and pa <= pc else b if pb <= pc else c


def decode(path):
    data = open(path, "rb").read()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        raise RuntimeError("%s is not a PNG" % path)
    off, raw = 8, b""
    width = height = color_type = None
    while off < len(data):
        n = struct.unpack(">I", data[off:off + 4])[0]
        typ = data[off + 4:off + 8]
        chunk = data[off + 8:off + 8 + n]
        off += 12 + n
        if typ == b"IHDR":
            width, height, depth, color_type, _, _, interlace = struct.unpack(">IIBBBBB", chunk)
            if depth != 8 or interlace != 0:
                raise RuntimeError("%s uses unsupported PNG encoding" % path)
        elif typ == b"IDAT":
            raw += chunk
    bpp = 3 if color_type == 2 else 4 if color_type == 6 else None
    if bpp is None:
        raise RuntimeError("%s uses unsupported PNG color type" % path)
    encoded = zlib.decompress(raw)
    prev = bytearray(width * bpp)
    pos, rgb, nonblack = 0, bytearray(), 0
    for _ in range(height):
        filt = encoded[pos]
        pos += 1
        row = bytearray(encoded[pos:pos + width * bpp])
        pos += width * bpp
        for i in range(len(row)):
            left = row[i - bpp] if i >= bpp else 0
            up = prev[i]
            upper_left = prev[i - bpp] if i >= bpp else 0
            if filt == 1:
                row[i] = (row[i] + left) & 255
            elif filt == 2:
                row[i] = (row[i] + up) & 255
            elif filt == 3:
                row[i] = (row[i] + ((left + up) // 2)) & 255
            elif filt == 4:
                row[i] = (row[i] + paeth(left, up, upper_left)) & 255
            elif filt != 0:
                raise RuntimeError("%s uses unsupported PNG filter" % path)
        for i in range(0, len(row), bpp):
            pixel = row[i:i + 3]
            rgb.extend(pixel)
            if pixel != b"\0\0\0":
                nonblack += 1
        prev = row
    return width, height, bytes(rgb), nonblack


partial = decode(sys.argv[1])
final = decode(sys.argv[2])
if partial[:2] != final[:2]:
    raise RuntimeError("partial and final ERT frames have different dimensions")
pixel_count = partial[0] * partial[1]
if not pixel_count // 20 < partial[3] < pixel_count * 19 // 20:
    raise RuntimeError("early ERT frame is not meaningfully partial: %d/%d" % (partial[3], pixel_count))
if final[3] <= partial[3]:
    raise RuntimeError("final ERT frame did not advance beyond partial frame")
changed = sum(a != b for a, b in zip(partial[2], final[2]))
if changed < pixel_count // 10:
    raise RuntimeError("partial and final ERT frames changed too little: %d" % changed)
print("ert_progress partial_lit=%d final_lit=%d changed_channels=%d" %
      (partial[3], final[3], changed))
PY

exit 0
