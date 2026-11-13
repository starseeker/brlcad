#!/bin/sh
#                  X V F B _ T E S T . S H
# BRL-CAD
#
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.

# A discovered xvfb-run binary is not proof that it can start an X server.
# Treat an unavailable test display as infrastructure, but leave the actual
# command's status unchanged once the preflight succeeds.
if [ "$#" -lt 2 ]; then
    echo "xvfb test wrapper requires xvfb-run and a command" >&2
    exit 2
fi

xvfb_run=$1
shift

if ! "$xvfb_run" -a /bin/sh -c 'test -n "$DISPLAY"' >/dev/null 2>&1; then
    echo "SKIP: xvfb-run cannot start a usable X display" >&2
    exit 125
fi

exec "$xvfb_run" -a "$@"
