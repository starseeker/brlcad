#!/bin/sh

set -eu

if ! command -v ldd >/dev/null 2>&1; then
  exit 77
fi

if ldd "$1" | grep -E 'libX11|libXft|libfontconfig|libfreetype'; then
  echo "Xmin build links to host X11 or font libraries" >&2
  exit 1
fi
