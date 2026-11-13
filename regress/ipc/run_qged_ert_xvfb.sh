#!/bin/sh
# Run the QGED ERT integration test with a CTest-owned isolated X server.
# Infrastructure failures are a CTest skip; failures after a verified display
# belongs to the QGED/ERT test and remain failures.

SKIP=77
XVFB_RUN="${1:-}"
TEST_SCRIPT="${2:-}"

if [ ! -x "$XVFB_RUN" ]; then
    echo "SKIP: xvfb-run is unavailable; cannot run the isolated QGED ERT test"
    exit "$SKIP"
fi

if [ ! -x "$TEST_SCRIPT" ]; then
    echo "FAIL: QGED ERT test script is unavailable: $TEST_SCRIPT" >&2
    exit 2
fi

TMPDIR_TEST=$(mktemp -d "${TMPDIR:-/tmp}/qged_ert_xvfb.XXXXXX") || {
    echo "SKIP: cannot create temporary directory for isolated X server"
    exit "$SKIP"
}
XERRLOG="$TMPDIR_TEST/xvfb.err"

cleanup() {
    rm -rf "$TMPDIR_TEST"
}
trap 'cleanup' EXIT

"$XVFB_RUN" -a -e "$XERRLOG" -s "-screen 0 1280x960x24" \
    /bin/sh -c '
        if ! command -v xdpyinfo >/dev/null 2>&1; then
            echo "SKIP: xdpyinfo is unavailable; cannot verify isolated X" >&2
            exit 77
        fi
        if ! xdpyinfo -display "$DISPLAY" >/dev/null 2>&1; then
            echo "SKIP: xvfb-run did not provide a usable X display" >&2
            exit 77
        fi
        exec "$@"
    ' qged-ert-xvfb "$TEST_SCRIPT"
status=$?

if [ "$status" -eq 0 ]; then
    exit 0
fi

if [ "$status" -eq "$SKIP" ]; then
    exit "$SKIP"
fi

# xvfb-run itself can fail before the test gets a DISPLAY.  Map only its known
# startup diagnostics to an infrastructure skip; all other results are product
# failures from the child test.
if [ -f "$XERRLOG" ] && grep -Eqi \
    'xvfb-run: error:|Xvfb failed to start|Cannot establish any listening sockets|Server is already active|Unable to open X display' \
    "$XERRLOG"; then
    echo "SKIP: isolated X server failed to start" >&2
    cat "$XERRLOG" >&2
    exit "$SKIP"
fi

echo "FAIL: QGED ERT test failed after isolated X startup" >&2
exit "$status"
