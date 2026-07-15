# BRL-CAD
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.

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
update

set dm_diagnostics [g dm diagnostics -V v1]
if {![regexp {host\.capabilities=0x([0-9a-fA-F]+)} $dm_diagnostics -> host_caps]} {
    puts stderr "Tk Obol diagnostics omitted host input capabilities"
    exit 1
}
if {[scan $host_caps %x host_cap_value] != 1 ||
        ($host_cap_value & 0x20) == 0} {
    puts stderr "Tk Obol X11 host did not advertise endpoint input capability"
    exit 1
}

# TkObol owns the shared view-key profile.  The generic Tcl bindings must not
# also apply those actions through the containing toplevel.
foreach key {3 4 f F R r l L t T b B} {
    if {[bind .v1 $key] ne ""} {
	puts stderr "TclCAD retained an overlapping Tk view binding for $key"
	exit 1
    }
}
g adc v1 draw 0
set input_tag [lindex [bindtags .v1] 0]
set input_command [lindex [bind $input_tag <KeyPress>] 0]
$input_command key press A A 0 0 0
update
if {[g dm get -V v1 view.faceplate.adc.visible] ne "true"} {
    puts stderr "Tk Obol endpoint key callback did not apply a semantic action"
    exit 1
}

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
foreach property {
    view.faceplate.params.font_size
    view.faceplate.center_dot.font_size
    view.faceplate.scale.font_size
} {
    if {[g dm get -V v1 $property] ne "18"} {
        puts stderr "TclCAD fontsize command did not update $property"
        exit 1
    }
}
foreach {kind property expected} {
    center_dot view.faceplate.center_dot.color {0 255 0}
    view_params view.faceplate.params.color {0 0 255}
    view_scale view.faceplate.scale.color {255 0 0}
} {
    g faceplate v1 $kind color {*}$expected
    set normalized [format "%g/%g/%g" \
            [expr {[lindex $expected 0] / 255.0}] \
            [expr {[lindex $expected 1] / 255.0}] \
            [expr {[lindex $expected 2] / 255.0}]]
    if {[g dm get -V v1 $property] ne $normalized} {
        puts stderr "TclCAD faceplate color command did not update $property"
        exit 1
    }
}
foreach {kind property} {
    center_dot view.faceplate.center_dot.visible
    view_params view.faceplate.params.visible
    view_scale view.faceplate.scale.visible
} {
    foreach {enabled expected} {0 false 1 true} {
        g faceplate v1 $kind draw $enabled
        if {[g dm get -V v1 $property] ne $expected} {
            puts stderr "TclCAD faceplate draw command did not set $property"
            exit 1
        }
    }
}
if {![catch {g faceplate v1 center_dot draw 2}]} {
    puts stderr "TclCAD faceplate draw command accepted an invalid endpoint Boolean"
    exit 1
}
foreach {command property} {
    model_axes view.faceplate.model_axes.visible
    view_axes view.faceplate.view_axes.visible
} {
    foreach {enabled expected} {0 false 1 true} {
        g $command v1 draw $enabled
        if {[g dm get -V v1 $property] ne $expected} {
            puts stderr "TclCAD $command draw command did not set $property"
            exit 1
        }
    }
}
if {![catch {g model_axes v1 draw 2}]} {
    puts stderr "TclCAD model_axes draw command accepted an invalid endpoint Boolean"
    exit 1
}
g model_axes v1 axes_pos 1 2 3
g model_axes v1 axes_size 4
g model_axes v1 axes_color 10 20 30
g model_axes v1 label_color 40 50 60
g model_axes v1 line_width 2
g model_axes v1 pos_only 1
g model_axes v1 triple_color 1
g model_axes v1 tick_enable 1
g model_axes v1 tick_length 3
g model_axes v1 tick_major_length 4
g model_axes v1 tick_interval 5
g model_axes v1 ticks_per_major 6
g model_axes v1 tick_threshold 7
g model_axes v1 tick_color 70 80 90
g model_axes v1 tick_major_color 100 110 120
foreach {property expected} {
    view.faceplate.model_axes.position.x 1
    view.faceplate.model_axes.position.z 3
    view.faceplate.model_axes.size 4
    view.faceplate.model_axes.line_width 2
    view.faceplate.model_axes.position_only true
    view.faceplate.model_axes.triple_color true
    view.faceplate.model_axes.ticks.visible true
    view.faceplate.model_axes.ticks.length 3
    view.faceplate.model_axes.ticks.major_length 4
    view.faceplate.model_axes.ticks.interval 5
    view.faceplate.model_axes.ticks.per_major 6
    view.faceplate.model_axes.ticks.threshold 7
} {
    if {[g dm get -V v1 $property] ne $expected} {
        puts stderr "TclCAD model axes command did not set $property"
        exit 1
    }
}
foreach {property expected} {
    view.faceplate.model_axes.color 0.039215686274509803/0.078431372549019607/0.11764705882352941
    view.faceplate.model_axes.labels.color 0.15686274509803921/0.19607843137254902/0.23529411764705882
    view.faceplate.model_axes.ticks.color 0.27450980392156865/0.31372549019607843/0.35294117647058826
    view.faceplate.model_axes.ticks.major_color 0.39215686274509803/0.43137254901960786/0.47058823529411764
} {
    if {[g dm get -V v1 $property] ne $expected} {
        puts stderr "TclCAD model axes color command did not set $property"
        exit 1
    }
}
g view_axes v1 axes_size 8
g view_axes v1 tick_enable 1
g view_axes v1 tick_color 1 2 3
if {[g dm get -V v1 view.faceplate.view_axes.size] ne "8" ||
        [g dm get -V v1 view.faceplate.view_axes.ticks.visible] ne "true" ||
        [g dm get -V v1 view.faceplate.view_axes.ticks.color] ne "0.0039215686274509803/0.0078431372549019607/0.011764705882352941"} {
    puts stderr "TclCAD view axes commands did not set endpoint policy"
    exit 1
}
foreach {enabled expected} {0 false 1 true} {
    g grid v1 draw $enabled
    if {[g dm get -V v1 view.faceplate.grid.visible] ne $expected} {
	puts stderr "TclCAD grid draw command did not set the endpoint visibility property"
	exit 1
    }
}
g grid v1 snap 1
if {[g dm get -V v1 view.faceplate.grid.snap] ne "true"} {
    puts stderr "TclCAD grid snap command did not set the endpoint policy"
    exit 1
}
g grid v1 anchor 7 8 9
if {[g dm get -V v1 view.faceplate.grid.anchor.x] ne "7" ||
        [g dm get -V v1 view.faceplate.grid.anchor.z] ne "9"} {
    puts stderr "TclCAD grid anchor command did not set endpoint coordinates"
    exit 1
}
g grid v1 rh 2
g grid v1 rv 3
g grid v1 mrh 4
g grid v1 mrv 5
if {[g dm get -V v1 view.faceplate.grid.resolution.horizontal] ne "2" ||
        [g dm get -V v1 view.faceplate.grid.resolution.vertical] ne "3" ||
        [g dm get -V v1 view.faceplate.grid.major.horizontal] ne "4" ||
        [g dm get -V v1 view.faceplate.grid.major.vertical] ne "5"} {
    puts stderr "TclCAD grid resolution commands did not set endpoint policy"
    exit 1
}
g grid v1 color 10 20 30
if {[g dm get -V v1 view.faceplate.grid.color] ne "0.039215686274509803/0.078431372549019607/0.11764705882352941"} {
    puts stderr "TclCAD grid color command did not set endpoint policy"
    exit 1
}
foreach {enabled expected} {1 true 0 false} {
    g adc v1 draw $enabled
    if {[g dm get -V v1 view.faceplate.adc.visible] ne $expected} {
	puts stderr "TclCAD adc draw command did not set the endpoint visibility property"
	exit 1
    }
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

# External renderers use the view-owned Obol/imgstream endpoint.  A uniform
# IPC framebuffer write must compose into the visible view without opening a
# second graphical framebuffer window.
g set_fb_mode v1 3
if {[g set_fb_mode v1] ne "3"} {
    puts stderr "TclCAD did not enable the Obol framebuffer overlay"
    exit 1
}
set framebuffer_ipc [g listen v1 ipc]
if {$framebuffer_ipc eq ""} {
    puts stderr "TclCAD did not provide an Obol framebuffer IPC endpoint"
    exit 1
}
if {[catch {
    exec [file join [bu_dir bin] fbclear] -F "ipc:$framebuffer_ipc" 23 45 67
} framebuffer_error]} {
    puts stderr "TclCAD Obol framebuffer IPC write failed: $framebuffer_error"
    exit 1
}
after 100
g refresh v1
update
set framebuffer_output "[file rootname [lindex $argv 1]]_framebuffer.png"
g png v1 $framebuffer_output
if {![file exists $framebuffer_output] || [file size $framebuffer_output] < 1000} {
    puts stderr "TclCAD Obol framebuffer overlay did not produce a non-empty image"
    exit 1
}
image create photo framebuffer_image
framebuffer_image read $framebuffer_output -format png
set framebuffer_pixel [framebuffer_image get 32 32]
image delete framebuffer_image
foreach actual $framebuffer_pixel expected {23 45 67} {
    if {abs($actual - $expected) > 1} {
        puts stderr "TclCAD Obol framebuffer overlay was not present in capture"
        exit 1
    }
}
g set_fb_mode v1 0
g refresh v1
update

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

# Keep a software pane alive while deleting the active hardware pane.  This
# exercises context promotion and the session framebuffer-provider handoff.
set handoff_output "[file rootname $output]_handoff[file extension $output]"
g new_view v2 tkobol sw
g refresh v2
update
if {![winfo exists .v2] || ![winfo ismapped .v2] ||
        ![winfo exists .v2.__obol] ||
        [winfo class .v2.__obol] ne "Label" ||
        [winfo width .v2] < 100 || [winfo height .v2] < 100} {
    puts stderr "Surviving software Tk Obol view was not mapped at a usable size"
    exit 1
}

g delete_view v1
update
if {[winfo exists .v1]} {
    puts stderr "Tk Obol view survived delete_view"
    exit 1
}
g autoview v2
g refresh v2
update
g png v2 $handoff_output
if {![file exists $handoff_output] || [file size $handoff_output] < 1000} {
    puts stderr "Surviving Tk Obol view did not capture after active-pane deletion"
    exit 1
}
g delete_view v2
update
if {[winfo exists .v2]} {
    puts stderr "Surviving Tk Obol view was not released after active-pane deletion"
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
