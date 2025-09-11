#!/bin/sh

# Test script for Terrascape DSP tessellation backend
#
# This script tests both classic and Terrascape tessellation backends
# to ensure they produce valid results.

# get common paths
LDIR=$1 . "$1/regress/dsp/dsp-common.sh"

FAILED=0

# Test data - use dsp-2-1 (2x2 DSP)
WID=2
LEN=$WID
BASE=dsp-terrascape-test

echo "=== Testing DSP Terrascape tessellation backends ==="

# Setup test files
LOG_CLASSIC=$BASE-classic.log
LOG_TERRASCAPE=$BASE-terrascape.log
TGM_CLASSIC=$BASE-classic.g
TGM_TERRASCAPE=$BASE-terrascape.g

TRASH="$TGM_CLASSIC $TGM_TERRASCAPE $LOG_CLASSIC $LOG_TERRASCAPE $BASE.pix $BASE.bw $BASE.dsp"
rm -f $TRASH

# Convert test DSP data
DSPASC=$1/regress/dsp/dsp-2-1.asc
echo "Converting DSP data from ASCII format..."
$A2P < $DSPASC > $BASE.pix
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to convert ASCII to PIX"
    FAILED=1
fi

$P2B -B1.0 $BASE.pix > $BASE.bw
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to convert PIX to BW"
    FAILED=1
fi

$CV huc nu16 $BASE.bw $BASE.dsp
if [ $? -ne 0 ]; then
    echo "ERROR: Failed to convert BW to DSP"
    FAILED=1
fi

# Test classic tessellation backend
echo "Testing classic tessellation backend..."
unset DSP_USE_TERRASCAPE
export DSP_USE_TERRASCAPE=""

$MGED -c > $LOG_CLASSIC 2>&1 <<EOF
opendb $TGM_CLASSIC y
in $BASE-classic.s dsp f $BASE.dsp $WID $LEN 0 ad 1 1
r $BASE-classic.r u $BASE-classic.s
quit
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Classic tessellation failed"
    cat $LOG_CLASSIC
    FAILED=1
else
    echo "Classic tessellation: SUCCESS"
fi

# Test Terrascape tessellation backend
echo "Testing Terrascape tessellation backend..."
export DSP_USE_TERRASCAPE=1

$MGED -c > $LOG_TERRASCAPE 2>&1 <<EOF
opendb $TGM_TERRASCAPE y
in $BASE-terrascape.s dsp f $BASE.dsp $WID $LEN 0 ad 1 1
r $BASE-terrascape.r u $BASE-terrascape.s
quit
EOF

if [ $? -ne 0 ]; then
    echo "ERROR: Terrascape tessellation failed"
    cat $LOG_TERRASCAPE
    FAILED=1
else
    echo "Terrascape tessellation: SUCCESS"
fi

# Clean up
if [ $FAILED -eq 0 ]; then
    echo "Cleaning up test files..."
    rm -f $TRASH
else
    echo "Test failed, leaving files for debugging: $TRASH"
fi

unset DSP_USE_TERRASCAPE

echo "=== DSP Terrascape tessellation test complete ==="

exit $FAILED