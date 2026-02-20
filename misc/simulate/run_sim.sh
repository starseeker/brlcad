#!/bin/bash
# run_sim.sh - Run truck-on-terrain physics simulation and render to video
#
# Usage: ./run_sim.sh [brlcad_build_dir] [output_dir]
# Example: ./run_sim.sh /home/runner/work/brlcad/brlcad-build /tmp/sim_output
#
# The script:
#  1. Creates a combined .g database (truck + terrain)
#  2. Steps through the simulation 25x 0.5s = 12.5s total
#  3. Saves a .g snapshot after each step
#  4. Renders each snapshot with rt
#  5. Assembles frames into an mp4 video

set -e

BRLCAD_BUILD=${1:-/home/runner/work/brlcad/brlcad-build}
OUTPUT_DIR=${2:-/tmp/sim_output}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export PATH="$BRLCAD_BUILD/bin:$PATH"
export LD_LIBRARY_PATH="$BRLCAD_BUILD/lib:$LD_LIBRARY_PATH"

DB_DIR="$BRLCAD_BUILD/share/db"
SIM_DB="$OUTPUT_DIR/truck_terrain_sim.g"
FRAMES_DIR="$OUTPUT_DIR/frames"

mkdir -p "$OUTPUT_DIR" "$FRAMES_DIR"

# ----------------------------------------------------------------
# Step 1: Build simulation database
# Run from DB_DIR so terra.bin is in CWD (required by DSP primitive)
# All operations (dbconcat, comb, attr, arced) are metadata-only
# and complete in milliseconds.
# ----------------------------------------------------------------
echo "=== Step 1: Building simulation database ==="
cd "$DB_DIR"
rm -f "$SIM_DB"

mged -c "$SIM_DB" < "$SCRIPT_DIR/setup_sim.tcl"
echo "  Created: $SIM_DB"

# ----------------------------------------------------------------
# Step 2: Step through simulation, saving .g snapshot each step
# 25 steps x 0.5s = 12.5 seconds total
#   t=0..10s: truck falls ~490m (free fall under gravity)
#   t=10..12.5s: truck hits and settles on terrain
# ----------------------------------------------------------------
NSTEPS=25
STEP_DUR=0.5

echo "=== Step 2: Running simulation ($NSTEPS steps x ${STEP_DUR}s = $(echo "$NSTEPS * $STEP_DUR" | bc)s) ==="
cd "$DB_DIR"

for i in $(seq 1 $NSTEPS); do
    STEPNUM=$(printf "%04d" $i)
    if [ "$i" -eq 1 ]; then
        mged -c "$SIM_DB" "simulate --steps 1 scene.c $STEP_DUR" 2>/dev/null
    else
        mged -c "$SIM_DB" "simulate --resume --steps 1 scene.c $STEP_DUR" 2>/dev/null
    fi
    cp "$SIM_DB" "$FRAMES_DIR/frame_${STEPNUM}.g"
    echo "  Step $i/$NSTEPS done"
done
echo "  Simulation frames saved to $FRAMES_DIR/"

# ----------------------------------------------------------------
# Step 3: Render each frame with rt
# Render scene.c/truck.sim path so the truck appears at its world
# position (after the physics matrix update from simulate).
# Use -a 225 -e 35 for a consistent diagonal view.
# Auto-centering produces a close-up of the truck in each frame.
# For the landing frames, the terrain surface also appears below.
# ----------------------------------------------------------------
echo "=== Step 3: Rendering $NSTEPS frames with rt ==="
cd "$DB_DIR"

IMG_W=1280
IMG_H=720

for i in $(seq 1 $NSTEPS); do
    STEPNUM=$(printf "%04d" $i)
    FRAME_G="$FRAMES_DIR/frame_${STEPNUM}.g"
    FRAME_PIX="$FRAMES_DIR/frame_${STEPNUM}.pix"
    FRAME_PNG="$FRAMES_DIR/frame_${STEPNUM}.png"

    [ -f "$FRAME_G" ] || { echo "  Frame $i: .g not found, skipping"; continue; }

    rt -w $IMG_W -n $IMG_H -a 225 -e 35 \
       -o "$FRAME_PIX" \
       "$FRAME_G" "scene.c/truck.sim" terrain.sim 2>/dev/null

    "$BRLCAD_BUILD/bin/pix-png" -w $IMG_W -n $IMG_H < "$FRAME_PIX" > "$FRAME_PNG"
    rm -f "$FRAME_PIX"
    echo "  Rendered frame $i/$NSTEPS"
done
echo "  All frames rendered."

# ----------------------------------------------------------------
# Step 4: Assemble video with ffmpeg
# 25 frames at 2 fps = 12.5s input -> 24fps output via interpolation
# ----------------------------------------------------------------
echo "=== Step 4: Assembling video ==="

VIDEO_OUT="$OUTPUT_DIR/truck_terrain_simulation.mp4"

ffmpeg -y \
    -framerate 2 \
    -pattern_type glob \
    -i "$FRAMES_DIR/frame_*.png" \
    -c:v libx264 \
    -preset slow \
    -crf 18 \
    -pix_fmt yuv420p \
    -vf "fps=24,scale=${IMG_W}:${IMG_H}:flags=lanczos" \
    "$VIDEO_OUT"

echo ""
echo "=== DONE ==="
echo "Video: $VIDEO_OUT"
ls -lh "$VIDEO_OUT"
