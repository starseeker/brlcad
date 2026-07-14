# BRL-CAD
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.

if {$argc != 1} {
    puts stderr "Usage: new_view_cleanup.tcl database.g"
    exit 1
}

if {[llength [info commands dm_list]]} {
    puts stderr "retired dm_list command is still registered"
    exit 1
}
foreach command {dm_bestXType dm_validXType} {
    if {[llength [info commands $command]]} {
	puts stderr "retired Tcl display-manager selector $command is still registered"
	exit 1
    }
}

go_open g db [lindex $argv 0]
g new_view v1 nu
rename g {}

puts "TclCAD new_view cleanup passed"
