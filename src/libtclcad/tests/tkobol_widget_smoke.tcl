if {$argc != 2} {
    puts stderr "Usage: tkobol_widget_smoke.tcl database.g output.png"
    exit 1
}

package require Tk
go_open g db [lindex $argv 0]
g new_view v1 tkobol
g draw all.g
g autoview v1
g refresh v1
update

if {![winfo exists .v1] || ![winfo ismapped .v1] ||
        [winfo width .v1] < 100 || [winfo height .v1] < 100 ||
        [llength [winfo children .v1]] == 0} {
    puts stderr "Tk Obol view was not mapped at a usable size"
    exit 1
}

set output [lindex $argv 1]
g png v1 $output
if {![file exists $output] || [file size $output] < 1000} {
    puts stderr "Tk Obol view did not produce a non-empty image"
    exit 1
}

g delete_view v1
update
if {[winfo exists .v1]} {
    puts stderr "Tk Obol view survived delete_view"
    exit 1
}

frame .host -width 320 -height 240
pack .host -expand true -fill both
g new_view .host.v2 tkobol -t 0
pack .host.v2 -expand true -fill both
g autoview .host.v2
g refresh .host.v2
update
if {![winfo exists .host.v2] || ![winfo ismapped .host.v2] ||
        [winfo width .host.v2] < 100 || [winfo height .host.v2] < 100} {
    puts stderr "Embedded Tk Obol view was not mapped at a usable size"
    exit 1
}
g delete_view .host.v2
update
if {[winfo exists .host.v2]} {
    puts stderr "Embedded Tk Obol view survived delete_view"
    exit 1
}
destroy .host

set sw_output "[file rootname $output]_sw[file extension $output]"
g new_view vsw tkobol sw
g draw all.g
g autoview vsw
g refresh vsw
update
if {![winfo exists .vsw] || ![winfo ismapped .vsw] ||
        ![winfo exists .vsw.__obol] ||
        [winfo class .vsw.__obol] ne "Label" ||
        [winfo width .vsw] < 100 || [winfo height .vsw] < 100} {
    puts stderr "Software Tk Obol view was not mapped as a usable image widget"
    exit 1
}
g png vsw $sw_output
if {![file exists $sw_output] || [file size $sw_output] < 1000} {
    puts stderr "Software Tk Obol view did not produce a non-empty image"
    exit 1
}
g delete_view vsw
update
if {[winfo exists .vsw]} {
    puts stderr "Software Tk Obol view survived delete_view"
    exit 1
}

rename g {}
puts "TclCAD Tk Obol widget smoke passed"
exit 0
