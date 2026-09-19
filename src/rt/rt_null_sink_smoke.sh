#!/bin/sh
#            R T _ N U L L _ S I N K _ S M O K E . S H
# BRL-CAD

if [ "$#" -ne 2 ]; then
    echo "Usage: rt_null_sink_smoke.sh rt model.g" >&2
    exit 2
fi

# /dev/null is a diagnostic raytrace sink, not a window or image target.
# Keep every common display environment unset to catch accidental provider
# initialization in this no-presentation path.
env -u DISPLAY -u WAYLAND_DISPLAY -u QT_QPA_PLATFORM \
    timeout 30s "$1" -F /dev/null -s 32 -P 2 "$2" all.g >/dev/null 2>&1
