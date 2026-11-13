# BRL-CAD
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.

if {$argc != 2} {
    puts stderr "Usage: archer_obol_smoke.tcl source.g working.g"
    exit 1
}

package require Tk
package require Itcl
package require Itk
package require Iwidgets

# Exercise ArcherCore as an embedded widget.  Unlike RaytraceWizard's
# tree-less instance, this production workflow keeps the hierarchy available
# so GUI selection can be checked against the shared GED selection service.
namespace eval ArcherCore {
    set parentClass itk::Widget
    set inheritFromToplevel 0
}
package require ArcherCore 1.0

proc capture_pane {app pane path} {
    file delete -force $path
    if {[catch {$app gedCmd pane_png $pane $path} capture_error] ||
	    ![file exists $path] || [file size $path] < 500} {
	puts stderr "ArcherCore capture failed for $path: $capture_error"
	exit 1
    }
}

proc image_delta {a b} {
    image create photo archer_delta_a
    image create photo archer_delta_b
    archer_delta_a read $a -format png
    archer_delta_b read $b -format png
    if {[image width archer_delta_a] != [image width archer_delta_b] ||
	    [image height archer_delta_a] != [image height archer_delta_b]} {
	image delete archer_delta_a archer_delta_b
	return -1
    }
    set delta 0
    for {set y 0} {$y < [image height archer_delta_a]} {incr y 4} {
	for {set x 0} {$x < [image width archer_delta_a]} {incr x 4} {
	    if {[archer_delta_a get $x $y] ne [archer_delta_b get $x $y]} {
		incr delta
	    }
	}
    }
    image delete archer_delta_a archer_delta_b
    return $delta
}

proc tree_child_named {tree parent name} {
    foreach node [$tree children $parent] {
	if {[$tree item $node -text] eq $name} {
	    return $node
	}
    }
    return ""
}

proc smoke_phase {name} {
    if {[info exists ::env(ARCHER_SMOKE_TRACE)]} {
	puts stderr "ARCHER_SMOKE_PHASE $name"
    }
}

set source_db [lindex $argv 0]
set working_db [lindex $argv 1]
file delete -force $working_db
file copy -force $source_db $working_db

set app [ArcherCore .archer_obol_smoke 1 1 0 1]
pack $app -expand true -fill both
wm geometry . 640x480
update

$app opendb $working_db
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
smoke_phase initial_draw

# Exercise actual hierarchy selection in both directions.  A tree click must
# publish into the authoritative GED selection set, and a command-line
# selection must move Archer's active row without any renderer polling.
set tree [$app component newtree]
$app rebuildTree
update
set all_node [tree_child_named $tree {} all.g]
if {$all_node eq ""} {
    puts stderr "ArcherCore hierarchy omitted all.g"
    exit 1
}
$tree focus $all_node
$tree item $all_node -open true
event generate $tree <<TreeviewOpen>>
update
set platform_node [tree_child_named $tree $all_node platform.r]
if {$platform_node eq ""} {
    puts stderr "ArcherCore hierarchy did not expand all.g/platform.r"
    exit 1
}
set pane [$app gedCmd pane]
set unselected_capture [file join [pwd] archer_obol_unselected.png]
set selected_capture [file join [pwd] archer_obol_selected.png]
capture_pane $app $pane $unselected_capture
$tree selection set $platform_node
event generate $tree <<TreeviewSelect>>
update
if {[string trim [$app gedCmd select list default]] ne "all.g/platform.r" ||
	[$app getSelectedTreePaths] ne "all.g/platform.r"} {
    puts stderr "ArcherCore tree selection did not reach shared GED state"
    exit 1
}
$app gedCmd refresh_all
update
capture_pane $app $pane $selected_capture
if {[image_delta $unselected_capture $selected_capture] < 1} {
    puts stderr "ArcherCore selected hierarchy path was not visibly highlighted"
    exit 1
}
$app gedCmd select clear
$app gedCmd select add all.g/cone.r
update
if {[$app getSelectedTreePaths] ne "all.g/cone.r"} {
    puts stderr "ArcherCore did not reflect command-line selection in its tree"
    exit 1
}
$app gedCmd select clear
update
if {[$app getSelectedTreePaths] ne ""} {
    puts stderr "ArcherCore retained a tree selection after GED selection clear"
    exit 1
}
smoke_phase selection

# Resizing the containing application is a camera-neutral presentation
# transition.  Exercise one large resize and a burst of alternating sizes
# after drawing, then require every hosted endpoint to report its new physical
# widget extent without changing its model-space view.  This is intentionally
# performed before the lifecycle teardown: it covers the same four-pane host
# configuration used by Archer and by rtwizard's embedded ArcherCore.
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
smoke_phase resize

# Switch the same retained scene between canonical wire and shaded modes and
# require a visible result from each, rather than accepting command success as
# evidence that an endpoint consumed the new mode.
set wire_capture [file join [pwd] archer_obol_wire.png]
set shaded_capture [file join [pwd] archer_obol_shaded.png]
$app gedCmd erase all.g
$app gedCmd draw -m0 all.g
$app gedCmd refresh_all
update
capture_pane $app $pane $wire_capture
$app gedCmd erase all.g
$app gedCmd draw -m1 all.g
$app gedCmd refresh_all
update
capture_pane $app $pane $shaded_capture
if {[image_delta $wire_capture $shaded_capture] < 10} {
    puts stderr "ArcherCore wire/shaded mode switch did not change visible output"
    exit 1
}
smoke_phase modes

# Camera rotation and zoom must preserve center while changing only the
# requested orientation and scale.
set center_before [$app gedCmd pane_center $pane]
set aet_before [$app gedCmd pane_aet $pane]
set size_before [$app gedCmd pane_size $pane]
$app gedCmd ae 90 0
$app gedCmd zoom 1.25
$app gedCmd refresh_all
update
set center_after [$app gedCmd pane_center $pane]
set aet_after [$app gedCmd pane_aet $pane]
set size_after [$app gedCmd pane_size $pane]
if {$center_after ne $center_before || $aet_after eq $aet_before ||
	$size_after == $size_before} {
    puts stderr "ArcherCore rotate/zoom violated its camera contract: center=$center_before->$center_after aet=$aet_before->$aet_after size=$size_before->$size_after"
    exit 1
}
smoke_phase camera

# Intermediate librt editing must update the retained scene, cancel must
# restore the database-backed shape, and commit must leave no active session.
set edit_path all.g/ellipse.r/ellipse.s
$app gedCmd erase all.g
$app gedCmd draw -m1 all.g/ellipse.r
$app gedCmd autoview_all
$app gedCmd refresh_all
update
set edit_before [file join [pwd] archer_obol_edit_before.png]
set edit_preview [file join [pwd] archer_obol_edit_preview.png]
capture_pane $app $pane $edit_before
set disk_before [$app gedCmd get ellipse.s]
$app gedCmd edit -i $edit_path a 24
set edit_status [$app gedCmd edit $edit_path status]
if {[string first "active" $edit_status] < 0 ||
	    [string first "dirty=1" $edit_status] < 0} {
    puts stderr "ArcherCore did not expose active intermediate edit state: $edit_status"
    exit 1
}
smoke_phase edit_preview
$app gedCmd refresh_all
update
capture_pane $app $pane $edit_preview
if {[image_delta $edit_before $edit_preview] < 10} {
    puts stderr "ArcherCore intermediate edit did not update visible geometry"
    exit 1
}
$app gedCmd edit $edit_path cancel
if {[$app gedCmd get ellipse.s] ne $disk_before} {
    puts stderr "ArcherCore edit cancel changed committed geometry"
    exit 1
}
$app gedCmd edit -i $edit_path b 18
$app gedCmd edit $edit_path commit
if {[$app gedCmd get ellipse.s] eq $disk_before ||
	    ![catch {$app gedCmd edit $edit_path status}]} {
    puts stderr "ArcherCore edit commit did not publish a clean committed state"
    exit 1
}
smoke_phase edit_commit

# Drive Archer's real RtControl endpoint and prove all framebuffer composition
# modes consume the streamed raytrace.  The model is small and the preview is
# deliberately bounded so this remains a routine GUI regression.
$app gedCmd erase $edit_path
$app gedCmd draw -m1 all.g
$app gedCmd autoview_all
$app gedCmd view faceplate center_dot 1
$app gedCmd refresh_all
update
$app rtcntrl configure -nproc 1 -hsample 0 -size 128 -other "-A 0.9"
$app rtcntrl setActivePane $pane
$app rtcntrl setFbMode 1
$app rtcntrl raytracePlus
for {set i 0} {$i < 30} {incr i} {
    after 100
    update
}
set fb_underlay [file join [pwd] archer_obol_fb_underlay.png]
set fb_interlay [file join [pwd] archer_obol_fb_interlay.png]
set fb_overlay [file join [pwd] archer_obol_fb_overlay.png]
capture_pane $app $pane $fb_underlay
$app rtcntrl setFbMode 2
$app gedCmd refresh_all
update
capture_pane $app $pane $fb_interlay
$app rtcntrl setFbMode 3
$app gedCmd refresh_all
update
capture_pane $app $pane $fb_overlay
if {[image_delta $fb_underlay $fb_overlay] < 10 ||
	    [image_delta $fb_interlay $fb_overlay] < 10} {
    puts stderr "ArcherCore framebuffer composition modes were not visibly distinct"
    exit 1
}
$app rtcntrl toggleFB
if {[$app gedCmd pane_set_fb_mode $pane] ne "0"} {
    puts stderr "ArcherCore framebuffer did not return to inactive mode"
    exit 1
}
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

foreach artifact [list $unselected_capture $selected_capture $wire_capture \
	$shaded_capture $edit_before $edit_preview $fb_underlay $fb_interlay \
	$fb_overlay] {
    file delete -force $artifact
}
file delete -force $working_db

puts "ArcherCore Tk Obol production workflow passed"
exit 0
