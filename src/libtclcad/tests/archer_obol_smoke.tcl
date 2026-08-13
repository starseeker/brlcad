# BRL-CAD
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.

if {$argc != 1} {
    puts stderr "Usage: archer_obol_smoke.tcl database.g"
    exit 1
}

package require Tk
package require Itcl
package require Itk
package require Iwidgets

# RaytraceWizard embeds ArcherCore as an itk::Widget.  Exercise that exact
# hosted-view configuration without bringing up the wizard's unrelated pages.
namespace eval ArcherCore {
    set parentClass itk::Widget
    set inheritFromToplevel 0
}
package require ArcherCore 1.0

set app [ArcherCore .archer_obol_smoke 1 1 1 1]
pack $app -expand true -fill both
wm geometry . 640x480
update

$app opendb [lindex $argv 0]
set ged [$app component ged]
set native_ged "${ged}_ged"
set views [$app gedCmd list_views]

if {[llength $views] != 4} {
    puts stderr "ArcherCore did not create four hosted GED views: $views"
    exit 1
}

foreach view $views {
    set diagnostics [$native_ged dm diagnostics -V $view]
    if {[string first "host=tk-" $diagnostics] < 0} {
        puts stderr "ArcherCore view $view is not hosted by Tk Obol: $diagnostics"
        exit 1
    }
    if {![winfo exists $view]} {
        puts stderr "ArcherCore view $view has no Tk widget"
        exit 1
    }
}

$app gedCmd draw all.g
$app gedCmd autoview_all
$app gedCmd refresh_all
update

# Resizing the containing application is a camera-neutral presentation
# transition.  Exercise one large resize and a burst of alternating sizes
# after drawing, then require every hosted endpoint to report its new physical
# widget extent without changing its model-space view.  This is intentionally
# performed before the lifecycle teardown: it covers the same four-pane host
# configuration used by Archer and by rtwizard's embedded ArcherCore.
set pane [$app gedCmd pane]
set active_view [$app gedCmd pane_win_name $pane]
set camera_before [list \
    [$app gedCmd pane_center $pane] \
    [$app gedCmd pane_size $pane] \
    [$app gedCmd pane_aet $pane]]
array set view_size_before {}
foreach view $views {
    set view_size_before($view) [list [winfo width $view] [winfo height $view]]
}

wm geometry . 900x700
after 100
update
foreach geometry {720x540 1040x760 760x580 920x680} {
    wm geometry . $geometry
    update idletasks
}
after 100
update

if {[winfo width .] != 920 || [winfo height .] != 680} {
    puts stderr "ArcherCore toplevel did not settle at its final resize: [winfo width .]x[winfo height .]"
    exit 1
}
set changed_views 0
foreach view $views {
    set actual [list [winfo width $view] [winfo height $view]]
    set reported [$native_ged view_win_size $view]
    # ArcherCore retains all four named endpoints even in a single-pane
    # layout.  Inactive panes are deliberately unmapped 1x1 Tk placeholders
    # with a retained nominal endpoint size; only mapped widgets have a
    # physical surface whose dimensions must agree exactly.
    if {[winfo ismapped $view] && [lindex $actual 0] >= 100 &&
	    [lindex $actual 1] >= 100} {
	if {$actual ne $reported} {
	    puts stderr "ArcherCore hosted endpoint size is stale for $view: actual=$actual reported=$reported"
	    exit 1
	}
	if {$actual ne $view_size_before($view)} {
	    incr changed_views
	}
    }
}
if {$changed_views < 1 || ![winfo ismapped $active_view] ||
	[winfo width $active_view] < 100 || [winfo height $active_view] < 100} {
    puts stderr "ArcherCore resize did not reach its mapped active pane: changed=$changed_views active=$active_view"
    exit 1
}
set camera_after [list \
    [$app gedCmd pane_center $pane] \
    [$app gedCmd pane_size $pane] \
    [$app gedCmd pane_aet $pane]]
if {$camera_after ne $camera_before} {
    puts stderr "ArcherCore resize changed the active CAD camera: before=$camera_before after=$camera_after"
    exit 1
}

$app gedCmd refresh_all
update
set resize_capture [file join [pwd] archer_obol_resize.png]
file delete -force $resize_capture
if {[catch {$app gedCmd pane_png $pane $resize_capture} capture_error] ||
        ![file exists $resize_capture] || [file size $resize_capture] < 1000} {
    puts stderr "ArcherCore resize capture failed: $capture_error"
    exit 1
}
image create photo archer_obol_resize_image
archer_obol_resize_image read $resize_capture -format png
if {[image width archer_obol_resize_image] != [winfo width $active_view] ||
        [image height archer_obol_resize_image] != [winfo height $active_view]} {
    puts stderr "ArcherCore resize capture has stale dimensions: image=[image width archer_obol_resize_image]x[image height archer_obol_resize_image] view=[winfo width $active_view]x[winfo height $active_view]"
    exit 1
}
image delete archer_obol_resize_image
file delete -force $resize_capture

if {![llength [info commands $native_ged]]} {
    puts stderr "ArcherCore did not retain its native GED command"
    exit 1
}

itcl::delete object $app
update

if {[llength [info commands $native_ged]]} {
    puts stderr "ArcherCore left its native GED command behind"
    exit 1
}
foreach view $views {
    if {[winfo exists $view] || [llength [info commands $view]]} {
        puts stderr "ArcherCore left hosted view $view behind"
        exit 1
    }
}

puts "ArcherCore Tk Obol lifecycle smoke passed"
exit 0
