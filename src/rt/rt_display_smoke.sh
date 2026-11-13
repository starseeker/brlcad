#!/bin/sh
#                  R T _ D I S P L A Y _ S M O K E . S H
# BRL-CAD

if [ "$#" -ne 3 ]; then
    echo "Usage: rt_display_smoke.sh rt model.g display-target" >&2
    exit 2
fi

# A small visible render verifies that rt's workers request image presentation
# while the caller owns the native display-session event loop.
timeout 30s "$1" -F "$3" -s 32 -P 2 "$2" all.g >/dev/null 2>&1
