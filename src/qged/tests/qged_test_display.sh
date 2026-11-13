#!/usr/bin/env bash
# Shared display setup for qged integration-test scripts.

qged_test_display_available()
{
    [[ -n "${DISPLAY:-}" ]] || return 1

    # A stale DISPLAY is common after a desktop or test X server exits.  Let
    # callers use it when no probe is installed, but do not mistake a defined
    # variable for a working display when xdpyinfo can answer definitively.
    if command -v xdpyinfo >/dev/null 2>&1; then
	    xdpyinfo -display "$DISPLAY" >/dev/null 2>&1
	    return $?
    fi
    return 0
}

qged_test_ensure_display()
{
    local script="$1"
    shift

    [[ "${QT_QPA_PLATFORM:-}" == "offscreen" ]] && return 0
    qged_test_display_available && return 0

    if [[ "${QGED_TEST_XVFB_ACTIVE:-0}" == "1" ]]; then
	echo "ERROR: xvfb-run did not provide a working display" >&2
	exit 2
    fi
    if ! command -v xvfb-run >/dev/null 2>&1; then
	echo "ERROR: a working display or xvfb-run is required" >&2
	exit 2
    fi

    export QGED_TEST_XVFB_ACTIVE=1
    exec xvfb-run -a "$script" "$@"
}
