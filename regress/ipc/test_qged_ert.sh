#!/bin/sh
# test_qged_ert.sh — Integration test: qged ert command via IPC under Xvfb
#
# This script verifies the full Phase 3+5 code path end-to-end:
#   - qged starts under Xvfb with software-rasterizer mode (-s)
#   - A geometry object is drawn
#   - The "ert" GED command is issued via xdotool
#   - We confirm that rt was launched into the active endpoint's IPC stream
#     at the laid-out canvas size.
#
# Environment variables consumed:
#   QGED_BIN      path to the qged binary (required)
#   RT_BIN        path to the rt binary   (required)
#   TEST_DB       path to a .g file with objects (defaults to boolean-ops.g)
#   DISPLAY       verified isolated X display supplied by the CTest launcher
#
# Exit status: 0 = all checks passed, 1 = any check failed

set -e

BDIR="$(cd "$(dirname "$0")/../.." && pwd)"

QGED_BIN="${QGED_BIN:-${BDIR}/bin/qged}"
RT_BIN="${RT_BIN:-${BDIR}/bin/rt}"
TEST_DB="${TEST_DB:-${BDIR}/share/db/boolean-ops.g}"
DRAW_OBJECT="${DRAW_OBJECT:-all}"

PASS=0
FAIL=0
VERBOSE="${VERBOSE:-0}"
SKIP=77

pass() { PASS=$((PASS+1)); [ "$VERBOSE" = "1" ] && echo "PASS: $1"; }
fail() { FAIL=$((FAIL+1)); echo "FAIL: $1"; }
check() {
    label="$1"; cond="$2"
    if eval "test $cond" 2>/dev/null; then pass "$label"; else fail "$label"; fi
}

TMPDIR_TEST=$(mktemp -d /tmp/test_qged_ert.XXXXXX)
QGED_PID=""

cleanup() {
    if [ -n "$QGED_PID" ]; then kill "$QGED_PID" 2>/dev/null || true; fi
    rm -rf "$TMPDIR_TEST"
}
trap 'cleanup' EXIT

echo "=== qged ert IPC integration test ==="

# -------------------------------------------------------------------
# 1. Verify required binaries exist
# -------------------------------------------------------------------
check "qged binary exists"  "-x $QGED_BIN"
check "rt binary exists"    "-x $RT_BIN"
check "test .g file exists" "-f $TEST_DB"

if [ ! -x "$QGED_BIN" ] || [ ! -x "$RT_BIN" ] || [ ! -f "$TEST_DB" ]; then
    echo "SKIP: required binaries or test database not found"
    echo "Tests: 0/$((PASS+FAIL)) passed"
    exit "$SKIP"
fi

# -------------------------------------------------------------------
# 2. Check that xdotool is available (needed to send commands to qged)
# -------------------------------------------------------------------
if ! command -v xdotool >/dev/null 2>&1; then
    echo "SKIP: xdotool is unavailable; cannot drive QGED ERT"
    exit "$SKIP"
fi
HAVE_XDOTOOL=1
pass "xdotool available"

# -------------------------------------------------------------------
# 3. Verify the isolated X display supplied by the CTest launcher.
#    This test must not start a second X server: that hides launch failures and
#    permits the GUI to attach to a display outside the test's lifecycle.
# -------------------------------------------------------------------
display_ok=0
if [ -n "$DISPLAY" ]; then
    # Verify the display accepts connections (up to 5 s to allow startup)
    for i in $(seq 1 50); do
        if xdpyinfo -display "$DISPLAY" >/dev/null 2>&1; then
            display_ok=1; break
        fi
        sleep 0.1
    done
fi

check "X display ($DISPLAY) is usable" "$display_ok -eq 1"

if [ "$display_ok" -ne 1 ]; then
    echo "SKIP: CTest launcher did not provide a usable isolated X display"
    exit "$SKIP"
fi

# -------------------------------------------------------------------
# 4. Start qged under the X display in swrast mode
# -------------------------------------------------------------------
QGED_LOG="$TMPDIR_TEST/qged.log"
QGED_ERRLOG="$TMPDIR_TEST/qged_err.log"

export LD_LIBRARY_PATH="${BDIR}/lib:${LD_LIBRARY_PATH:-}"
# Keep saved user dock/window settings from invalidating the console coordinates.
export XDG_CONFIG_HOME="$TMPDIR_TEST/config"
mkdir -p "$XDG_CONFIG_HOME"

"$QGED_BIN" -s --exec "draw $DRAW_OBJECT;ert" "$TEST_DB" >"$QGED_LOG" 2>"$QGED_ERRLOG" &
QGED_PID=$!

# Wait up to 15 s for qged window to appear
qged_ready=0
for i in $(seq 1 150); do
    if xdotool search --pid "$QGED_PID" 2>/dev/null | grep -q .; then
        qged_ready=1; break
    fi
    sleep 0.1
done

check "qged window appeared under Xvfb" "$qged_ready -eq 1"

if [ "$qged_ready" -ne 1 ]; then
    echo "qged stdout: $(cat $QGED_LOG 2>/dev/null | head -20)"
    echo "qged stderr: $(cat $QGED_ERRLOG 2>/dev/null | head -20)"
    echo "ABORT: qged did not start"
    exit 1
fi

# Allow a little more time for qged to finish initialisation
sleep 2

check "qged process still running after init" "-d /proc/$QGED_PID"

# -------------------------------------------------------------------
# 5. Verify the startup commands launched ERT through the GUI event loop
# -------------------------------------------------------------------
if [ "$HAVE_XDOTOOL" -eq 1 ]; then
    # A small scene can complete before /proc can observe its child.  The
    # common endpoint launcher records the handoff before it starts rt, which
    # is the durable evidence that this GUI command used the endpoint stream.
    rt_appeared=0
    for i in $(seq 1 300); do
        if pgrep -x rt >/dev/null 2>&1 ||
	    grep -q "rt: launching endpoint framebuffer renderer" "$QGED_ERRLOG" 2>/dev/null; then
            rt_appeared=1; break
        fi
        sleep 0.1
    done
    check "rt process was spawned by ert" "$rt_appeared -eq 1"

    # Check that the launcher selected the IPC endpoint.  Sampling a finished
    # child process through /proc is intentionally not a test requirement.
    ipc_used=0

    if grep -q "rt: launching endpoint framebuffer renderer (ipc=1" "$QGED_ERRLOG" 2>/dev/null; then
	ipc_used=1
    fi
    check "rt launched through the endpoint IPC stream" "$ipc_used -eq 1"

    # Wait for rt to finish (up to 60 s); treat zombie as completed
    rt_done=0
    for i in $(seq 1 600); do
        RT_ALIVE=0
        for pid in $(pgrep -x rt 2>/dev/null); do
            state=$(awk '/^State/{print $2}' /proc/"$pid"/status 2>/dev/null)
            if [ "$state" != "Z" ]; then
                RT_ALIVE=1; break
            fi
        done
        if [ "$RT_ALIVE" -eq 0 ]; then
            rt_done=1; break
        fi
        sleep 0.1
    done
    check "rt process completed" "$rt_done -eq 1"

    # Reject the historical pre-layout framebuffer (typically only tens of
    # pixels wide/high) even if IPC and rendering otherwise complete.
    set -- $(sed -n 's/.*endpoint framebuffer renderer (ipc=[01] size=\([0-9][0-9]*\)x\([0-9][0-9]*\)).*/\1 \2/p' "$QGED_ERRLOG" | tail -1)
    FB_WIDTH="${1:-0}"
    FB_HEIGHT="${2:-0}"
    check "ert used laid-out canvas dimensions (${FB_WIDTH}x${FB_HEIGHT})" \
	"$FB_WIDTH -gt 200 -a $FB_HEIGHT -gt 200"

    # qged should still be alive
    check "qged process survived ert" "-d /proc/$QGED_PID"

    # Check logs for crash indicators
    if [ -f "$QGED_ERRLOG" ]; then
        if grep -q "Segmentation fault\|Aborted\|core dumped" "$QGED_ERRLOG" 2>/dev/null; then
            fail "qged did not crash during ert"
        else
            pass "qged did not crash during ert"
        fi
    fi
else
    echo "SKIP: xdotool not available; skipping interactive ert test"
fi

# -------------------------------------------------------------------
# 6. Summary
# -------------------------------------------------------------------
echo ""
echo "Tests: $PASS/$((PASS+FAIL)) passed"
if [ "$FAIL" -ne 0 ]; then
    echo "--- qged stdout ---"
    tail -40 "$QGED_LOG" 2>/dev/null || true
    echo "--- qged stderr ---"
    tail -80 "$QGED_ERRLOG" 2>/dev/null || true
fi
[ "$FAIL" -eq 0 ] && exit 0 || exit 1
