#!/bin/sh
#             F B S E R V _ D I S P L A Y _ S M O K E . S H
# BRL-CAD

if [ "$#" -ne 2 ]; then
    echo "Usage: fbserv_display_smoke.sh fbserv display-target" >&2
    exit 2
fi

# A live server is the success condition.  It is deliberately stopped after a
# short interval so the test does not depend on a client or a fixed TCP port.
timeout 2s "$1" -F "$2" -p 0 >/dev/null 2>&1
status=$?
test "$status" -eq 124
