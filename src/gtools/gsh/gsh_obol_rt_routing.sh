#!/bin/sh
set -eu

# A non-login Git for Windows sh does not add its own core utilities to PATH.
PATH="/usr/bin:/bin:${PATH}"
export PATH

if [ "$#" -ne 3 ]; then
    echo "Usage: gsh_obol_rt_routing.sh <gsh> <db> <workdir>" 1>&2
    exit 1
fi

GSH="$1"
DB="$2"
WORKDIR="$3"
LOG="${WORKDIR}/gsh_obol_rt_routing.log"

rm -f "$LOG"

# The first two renders must use the attached endpoint and report the
# dimensions explicitly requested by the caller.  The final render names its
# own framebuffer and must therefore bypass endpoint auto-routing.
printf 'draw all.g
autoview
ert -s 96
delay 3 0
ert -w 160 -n 80
delay 3 0
rt -s 64 -F /dev/null
delay 3 0
quit
' | "$GSH" --new-cmds "$DB" > "$LOG" 2>&1

if ! grep -Eq 'rt: launching endpoint framebuffer renderer \(ipc=[01] size=96x96\)' "$LOG"; then
    echo "endpoint square-size routing did not preserve 96x96" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi
if ! grep -Eq 'rt: launching endpoint framebuffer renderer \(ipc=[01] size=160x80\)' "$LOG"; then
    echo "endpoint width/height routing did not preserve 160x80" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

endpoint_launches=$(grep -c 'launching endpoint framebuffer renderer' "$LOG" || true)
if [ "$endpoint_launches" -ne 2 ]; then
    echo "explicit -F render was incorrectly auto-routed to the endpoint" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi
if grep -Eq 'could not open fb server|invalid embedded framebuffer dimensions|rt_gettrees.*FAILED' "$LOG"; then
    echo "framebuffer routing test reported a render failure" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

frames=$(grep -Ec 'Frame[[:space:]]+0:' "$LOG" || true)
if [ "$frames" -lt 3 ]; then
    echo "framebuffer routing test did not complete all three renders" 1>&2
    cat "$LOG" 1>&2
    exit 1
fi

exit 0
