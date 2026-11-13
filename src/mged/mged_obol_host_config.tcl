#                  M G E D _ O B O L _ H O S T _ C O N F I G . T C L
# BRL-CAD
#
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.

if {$argc != 1} {
    puts stderr "usage: mged_obol_host_config.tcl <openw.tcl>"
    exit 2
}

# openw.tcl initializes defaults and declares procedures.  Stub only its
# application/Tk setup hooks so this validates the script in a plain Tcl host.
proc check_externs args {}
proc loadtk args {}
proc bu_dir args { return /tmp }
proc color_scheme_init args {}
proc font_init args {}
proc unknown args { return "" }

source [lindex $argv 0]

if {![info exists mged_default(host_type)] ||
    $mged_default(host_type) ne "tkobol"} {
    puts stderr "MGED did not default the GUI to the tkobol host"
    exit 1
}

if {[llength [array names mged_default dm_type]] != 0} {
    puts stderr "MGED restored the retired dm_type GUI setting"
    exit 1
}

set gui_body [info body gui]
if {[string first "-host" $gui_body] < 0 ||
    [string first "-dt" $gui_body] >= 0} {
    puts stderr "MGED GUI host parsing is not Obol-only"
    exit 1
}
