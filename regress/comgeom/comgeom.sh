#!/bin/sh
#                     C O M G E O M . S H
# BRL-CAD
#
# Copyright (c) 2008-2025 United States Government as represented by
# the U.S. Army Research Laboratory.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
# copyright notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials provided
# with the distribution.
#
# 3. The name of the author may not be used to endorse or promote
# products derived from this software without specific prior written
# permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
###
#
# Test conversion from .g to GIFT (v5) and back.
#
###

# Ensure /bin/sh
export PATH || (echo "This isn't sh."; sh $0 $*; kill $$)

# source common library functionality, setting ARGS, NAME_OF_THIS,
# PATH_TO_THIS, and THIS.
. "$1/regress/library.sh"

if test "x$LOGFILE" = "x" ; then
    LOGFILE=`pwd`/comgeom.log
    rm -f $LOGFILE
fi
log "=== TESTING 'comgeom-g' and 'vdeck' ==="

ASC2G="`ensearch asc2g`"
if test ! -f "$ASC2G" ; then
    log "Unable to find asc2g, aborting"
    exit 1
fi
VDECK="`ensearch vdeck`"
if test ! -f "$VDECK" ; then
    log "Unable to find vdeck, aborting"
    exit 1
fi
COMGEOM="`ensearch comgeom-g`"
if test ! -f "$COMGEOM" ; then
    log "Unable to find comgeom-g, aborting"
    exit 1
fi
GZIP="`which gzip`"
if test ! -f "$GZIP" ; then
    log "Unable to find gzip, aborting"
    exit 1
fi

FAILURES=0

# make our starting database
rm -f comgeom.m35.asc
log "... running gzip decompress"
$GZIP -d -c "$PATH_TO_THIS/../iges/m35.asc.gz" > comgeom.m35.asc
rm -f comgeom.m35.g
run $ASC2G comgeom.m35.asc comgeom.m35.g

# get our
# using vdeck interactively to convert .g to GIFT
#(following example in red.sh and mged test)
rm -f solids regions region_ids
log "... running $VDECK comgeom.m35.g"
$VDECK comgeom.m35.g >> $LOGFILE 2>&1 <<EOF
i all.g
d
q
EOF

# assemble pieces to compare with test version
rm -f comgeom.m35.cg
cat solids     >  comgeom.m35.cg
cat regions    >> comgeom.m35.cg
cat region_ids >> comgeom.m35.cg

# get test version
rm -f comgeom.m35-baseline.cg
log "... running gzip decompress"
$GZIP -d -c "$PATH_TO_THIS/m35.cg.gz" > comgeom.m35-baseline.cg

cmp comgeom.m35.cg comgeom.m35-baseline.cg
STATUS=$?

if [ X$STATUS != X0 ] ; then
    log "vdeck results differ $STATUS"
    FAILURES="`expr $FAILURES + 1`"
else
    log "vdeck test succeeded (1 of 3)"
fi

# the part 2 test checks for a known (but corrected) failures to
# convert ars solids, it also checks all solid types recognized
# by vdeck
log "...testing 'comgeom-g' command (GIFT v5)..."
rm -f comgeom.t-v5.g
run $COMGEOM "$PATH_TO_THIS/comgeom.cg" comgeom.t-v5.g
STATUS=$?

if [ X$STATUS != X0 ] ; then
    log "comgeom-g conversion errors: $STATUS"
    FAILURES="`expr $FAILURES + 1`"
else
    ERR=`grep error $LOGFILE`
    if test "x$ERR" = "x" ; then
	log "comgeom-g v5 test succeeded (2 of 3)"
    else
	log "comgeom-g v5 errors, see $LOGFILE"
	FAILURES="`expr $FAILURES + 1`"
    fi
fi

# the part 3 test checks comgeom-g against a GIFT v4 tgm
log "...testing 'comgeom-g' command (GIFT v4)..."
rm -f comgeom.t-v4.g
run $COMGEOM -v4 "$PATH_TO_THIS/comgeom.tgt-1-v4.cg" comgeom.t-v4.g
STATUS=$?

if [ X$STATUS != X0 ] ; then
    log "comgeom-g conversion errors: $STATUS"
    FAILURES="`expr $FAILURES + 1`"
else
    ERR=`grep error $LOGFILE`
    if test "x$ERR" = "x" ; then
	log "comgeom-g v4 test succeeded (3 of 3)"
    else
	log "comgeom-g v4 errors, see $LOGFILE"
	FAILURES="`expr $FAILURES + 1`"
    fi
fi

if test $FAILURES -eq 0 ; then
    log "-> vdeck/comgeom-g check succeeded"
else
    log "-> vdeck/comgeom-g check FAILED, see $LOGFILE"
    cat "$LOGFILE"
fi

exit $FAILED


# Local Variables:
# tab-width: 8
# mode: sh
# sh-indentation: 4
# sh-basic-offset: 4
# indent-tabs-mode: t
# End:
# ex: shiftwidth=4 tabstop=8
