if {$argc != 2} {
    puts stderr "Usage: tkobol_widget_smoke.tcl database.g output.png"
    exit 1
}

package require Tk
go_open g db [lindex $argv 0]

if {[g dm open --host tk-photo --renderer sw] ne "tk-photo" ||
        [g dm host] ne "tk-photo"} {
    puts stderr "GED dm open did not create a top-level Tk software endpoint"
    exit 1
}
set dm_diagnostics [g dm diagnostics]
if {[string first "host=tk-photo" $dm_diagnostics] < 0 ||
        [string first "renderer=sw" $dm_diagnostics] < 0} {
    puts stderr "GED dm diagnostics omitted the Tk endpoint state"
    exit 1
}
g dm close
if {[g dm host] ne "unbound"} {
    puts stderr "GED dm close did not release the Tk endpoint host"
    exit 1
}

g new_view v1 tkobol
g draw all.g
g autoview v1

g dm set -V v1 endpoint.title "Tk typed endpoint title"
if {[g dm get -V v1 endpoint.title] ne "Tk typed endpoint title" ||
        [wm title .v1] ne "Tk typed endpoint title"} {
    puts stderr "GED dm title property did not reach the Tk toplevel host"
    exit 1
}
g dm set -V v1 endpoint.visible 0
update
if {[winfo ismapped .v1]} {
    puts stderr "GED dm visibility property did not withdraw the Tk toplevel"
    exit 1
}
g dm set -V v1 endpoint.visible 1
update
if {![winfo ismapped .v1]} {
    puts stderr "GED dm visibility property did not restore the Tk toplevel"
    exit 1
}


g bg v1 12 34 56
if {[g bg v1] ne "12 34 56"} {
    puts stderr "TclCAD background command did not use retained view state"
    exit 1
}
foreach command {light transparency zbuffer zclip} {
    set original [g $command v1]
    set changed [expr {!$original}]
    g $command v1 $changed
    if {[g $command v1] ne "$changed"} {
        puts stderr "TclCAD $command command did not round-trip endpoint state"
        exit 1
    }
    g $command v1 $original
}
g fontsize v1 18
if {[g fontsize v1] ne "18"} {
    puts stderr "TclCAD fontsize command did not update retained HUD state"
    exit 1
}
g dm set -V v1 renderer.clip.minimum 0.25
g dm set -V v1 renderer.clip.maximum 0.75
if {[g dm get -V v1 renderer.clip.minimum] ne "0.25" ||
        [g dm get -V v1 renderer.clip.maximum] ne "0.75"} {
    puts stderr "GED dm command did not round-trip retained clip planes"
    exit 1
}
set ray [g mouse_ray v1 100 100]
if {[llength $ray] != 2 || [lindex $ray 0] eq [lindex $ray 1]} {
    puts stderr "TclCAD mouse_ray did not produce an independent model ray"
    exit 1
}

g refresh v1
update

set pick ""
foreach py {128 192 256 320 384} {
    foreach px {128 192 256 320 384} {
        set pick [g mouse_pick_detail v1 $px $py]
        if {$pick ne ""} {
            break
        }
    }
    if {$pick ne ""} {
        break
    }
}
if {$pick eq "" || [lsearch -exact $pick primitive_index] < 0} {
    puts stderr "TclCAD mouse_pick_detail did not return retained pick identity"
    exit 1
}

if {![winfo exists .v1] || ![winfo ismapped .v1] ||
        [winfo width .v1] < 100 || [winfo height .v1] < 100 ||
        [llength [winfo children .v1]] == 0} {
    puts stderr "Tk Obol view was not mapped at a usable size"
    exit 1
}
set reported_size [g view_win_size v1]
if {[lindex $reported_size 0] != [winfo width .v1] ||
        [lindex $reported_size 1] != [winfo height .v1]} {
    puts stderr "TclCAD view dimensions did not follow the Tk Obol host"
    exit 1
}

set output [lindex $argv 1]
g png v1 $output
if {![file exists $output] || [file size $output] < 1000} {
    puts stderr "Tk Obol view did not produce a non-empty image"
    exit 1
}
set pix_output "[file rootname $output].pix"
g pix v1 $pix_output
if {![file exists $pix_output] || [file size $pix_output] < 1000} {
    puts stderr "Tk Obol view did not produce a non-empty PIX image"
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
set sw_pix_output "[file rootname $sw_output].pix"
g pix vsw $sw_pix_output
if {![file exists $sw_pix_output] || [file size $sw_pix_output] < 1000} {
    puts stderr "Software Tk Obol view did not produce a non-empty PIX image"
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
