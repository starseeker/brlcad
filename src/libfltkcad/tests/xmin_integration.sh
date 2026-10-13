#!/bin/sh

set -eu

xmin_run=$1
xminctl=$2
viewer=$3
page=$4
capture_directory=$5

"$xmin_run" -- sh -c '
  viewer=$1
  xminctl=$2
  page=$3
  before=$4/before.ppm
  linked=$4/linked.ppm
  after=$4/after.ppm
  link_x=100
  link_y=94
  scroll_steps=6
  ppm_header_bytes=15
  current_highlight_rgb="ff ae 42"

  "$viewer" "$page" needle &
  viewer_process=$!
  trap '\''kill "$viewer_process" 2>/dev/null || true'\'' EXIT HUP INT TERM

  "$xminctl" wait-window --timeout 5000 "FLTK litehtml viewer"
  "$xminctl" wait-stable --quiet 100 --timeout 5000 "FLTK litehtml viewer"
  "$xminctl" mouse-move "FLTK litehtml viewer" 300 220
  "$xminctl" capture-window "FLTK litehtml viewer" "$before"
  "$xminctl" click "FLTK litehtml viewer" "$link_x" "$link_y"
  "$xminctl" wait-stable --quiet 100 --timeout 5000 "FLTK litehtml viewer"
  "$xminctl" capture-window "FLTK litehtml viewer" "$linked"
  "$xminctl" scroll "$scroll_steps"
  "$xminctl" wait-stable --quiet 100 --timeout 5000 "FLTK litehtml viewer"
  "$xminctl" capture-window "FLTK litehtml viewer" "$after"
  "$xminctl" close "FLTK litehtml viewer"
  wait "$viewer_process"
  trap - EXIT HUP INT TERM

  test -s "$before"
  test -s "$linked"
  test -s "$after"
  if ! hexdump -v -s "$ppm_header_bytes" -e '\''3/1 "%02x " "\\n"'\'' "$before" | grep -q "$current_highlight_rgb"; then
    echo "current-match highlight color was not captured" >&2
    exit 1
  fi
  if cmp -s "$before" "$linked"; then
    echo "activating the fragment link did not change the captured window" >&2
    exit 1
  fi
  if cmp -s "$linked" "$after"; then
    echo "scrolling did not change the captured window" >&2
    exit 1
  fi
' fltkcad-xmin "$viewer" "$xminctl" "$page" "$capture_directory"
