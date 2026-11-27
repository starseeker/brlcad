#!/bin/sh
# BRL-CAD
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.
# SPDX-License-Identifier: BSD-3-Clause

mged_obol_find_python()
{
    if [ -n "${PYTHON:-}" ] &&
	"$PYTHON" -c 'import sys' >/dev/null 2>&1; then
	return 0
    fi

    PYTHON=
    for candidate in python3.14 python3.13 python3.12 python3.11 python3 python py; do
	if command -v "$candidate" >/dev/null 2>&1 &&
	    "$candidate" -c 'import sys' >/dev/null 2>&1; then
	    PYTHON="$candidate"
	    return 0
	fi
    done
    return 1
}
