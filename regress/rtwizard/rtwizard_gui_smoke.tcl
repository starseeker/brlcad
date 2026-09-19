# BRL-CAD
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.
#
# Loaded through RTWIZARD_RCFILE by the production rtwizard executable.  The
# callback runs from RaytraceWizard's event loop after the historical GUI has
# created its embedded ArcherCore view.

set ::rtwizard_gui_smoke_attempts 0

proc rtwizard_gui_smoke_fail {message} {
    puts stderr "rtwizard GUI smoke failed: $message"
    exit 1
}

proc rtwizard_gui_smoke_probe {} {
    incr ::rtwizard_gui_smoke_attempts
    if {$::rtwizard_gui_smoke_attempts > 200} {
	rtwizard_gui_smoke_fail "timed out waiting for the historical GUI"
    }

    if {![info exists ::mgedObj] || ![llength [info commands $::mgedObj]] ||
	    ![info exists ::fbp] || ![llength [info commands $::fbp]]} {
	after 100 rtwizard_gui_smoke_probe
	return
    }

    if {[catch {$::mgedObj component ged} ged] || $ged eq ""} {
	after 100 rtwizard_gui_smoke_probe
	return
    }
    set native_ged "${ged}_ged"
    if {![llength [info commands $native_ged]]} {
	after 100 rtwizard_gui_smoke_probe
	return
    }
    if {[catch {$::mgedObj gedCmd list_views} views] ||
	    [llength $views] != 4} {
	after 100 rtwizard_gui_smoke_probe
	return
    }

    foreach view $views {
	if {![winfo exists $view]} {
	    rtwizard_gui_smoke_fail "embedded view $view has no Tk widget"
	}
	if {[catch {$native_ged dm diagnostics -V $view} diagnostics] ||
		[string first "host=tk-" $diagnostics] < 0} {
	    rtwizard_gui_smoke_fail \
		"embedded view $view is not hosted by Tk Obol: $diagnostics"
	}
    }

    if {![info exists ::wizardInstance] ||
	    ![llength [info commands $::wizardInstance]]} {
	rtwizard_gui_smoke_fail "historical wizard object is unavailable"
    }
    if {[catch {$::wizardInstance select fullColor} page_error]} {
	after 100 rtwizard_gui_smoke_probe
	return
    }
    update
    if {![winfo ismapped $::mgedObj] || [winfo width $::mgedObj] < 100 ||
	    [winfo height $::mgedObj] < 100} {
	rtwizard_gui_smoke_fail \
	    "rendering page did not map its embedded ArcherCore: $page_error"
    }

    if {[catch {
	$::mgedObj gedCmd draw all.g
	$::mgedObj gedCmd autoview_all
	$::mgedObj gedCmd refresh_all
	update
    } draw_error]} {
	rtwizard_gui_smoke_fail "embedded drawing failed: $draw_error"
    }

    set pane [$::mgedObj gedCmd pane]

    # Resize the actual historical wizard, not a reconstructed test widget.
    # The embedded view must consume the new physical extent while preserving
    # its CAD camera and the already-drawn scene.  A short alternating burst
    # exercises configure-event coalescing before the framebuffer endpoint is
    # opened, and the capture below proves the resized composition reaches the
    # user-visible surface.
    set pane_view [$::mgedObj gedCmd pane_win_name $pane]
    set camera_before [list \
	[$::mgedObj gedCmd pane_center $pane] \
	[$::mgedObj gedCmd pane_size $pane] \
	[$::mgedObj gedCmd pane_aet $pane]]
    set size_before [list [winfo width $pane_view] [winfo height $pane_view]]
    foreach geometry {960x720 760x560 1020x740 820x620 940x700} {
	wm geometry . $geometry
	update idletasks
    }
    after 100
    update
    set size_after [list [winfo width $pane_view] [winfo height $pane_view]]
    if {$size_after eq $size_before || [lindex $size_after 0] < 100 ||
	    [lindex $size_after 1] < 100} {
	rtwizard_gui_smoke_fail \
	    "embedded view did not consume wizard resize: before=$size_before after=$size_after"
    }
    set reported_size [$::mgedObj gedCmd pane_win_size $pane]
    if {$reported_size ne $size_after} {
	rtwizard_gui_smoke_fail \
	    "embedded endpoint retained a stale physical size: actual=$size_after reported=$reported_size"
    }
    set camera_after [list \
	[$::mgedObj gedCmd pane_center $pane] \
	[$::mgedObj gedCmd pane_size $pane] \
	[$::mgedObj gedCmd pane_aet $pane]]
    if {$camera_after ne $camera_before} {
	rtwizard_gui_smoke_fail \
	    "wizard resize changed the CAD camera: before=$camera_before after=$camera_after"
    }

    if {[catch {$::fbp getFrameBuffer 64 64} framebuffer] ||
	    ![string match "ipc:*" $framebuffer]} {
	rtwizard_gui_smoke_fail \
	    "preview did not return an Obol framebuffer IPC endpoint: $framebuffer"
    }
	if {[catch {
	exec [file join [bu_dir bin] fbclear] -F $framebuffer 23 45 67
	after 100
	if {[$::mgedObj gedCmd pane] ne $pane} {
	    error "active pane changed while opening the preview endpoint"
	}
	if {[$::mgedObj gedCmd pane_set_fb_mode $pane] ne "3"} {
	    error "active pane did not retain framebuffer overlay mode"
	}
	$::mgedObj gedCmd pane_refresh $pane
	update
    } framebuffer_error]} {
	rtwizard_gui_smoke_fail \
	    "preview framebuffer write failed: $framebuffer_error"
    }

    if {![info exists ::env(RTWIZARD_GUI_SMOKE_OUTPUT)] ||
	    $::env(RTWIZARD_GUI_SMOKE_OUTPUT) eq ""} {
	rtwizard_gui_smoke_fail "RTWIZARD_GUI_SMOKE_OUTPUT is not set"
    }
    set output $::env(RTWIZARD_GUI_SMOKE_OUTPUT)
    file delete -force $output
    if {[catch {$::mgedObj gedCmd pane_png $pane $output} capture_error] ||
	    ![file exists $output] || [file size $output] < 1000} {
	rtwizard_gui_smoke_fail "preview capture failed: $capture_error"
    }
    image create photo rtwizard_gui_smoke_image
    if {[catch {
	rtwizard_gui_smoke_image read $output -format png
	set pixel [rtwizard_gui_smoke_image get 32 32]
    } image_error]} {
	rtwizard_gui_smoke_fail "cannot inspect preview capture: $image_error"
    }
    image delete rtwizard_gui_smoke_image
    foreach actual $pixel expected {23 45 67} {
	if {abs($actual - $expected) > 1} {
	    rtwizard_gui_smoke_fail \
		"preview capture omitted framebuffer color: $pixel"
	}
    }

    itcl::delete object $::mgedObj
    update
    if {[llength [info commands $native_ged]]} {
	rtwizard_gui_smoke_fail "embedded ArcherCore left its GED command behind"
    }
    foreach view $views {
	if {[winfo exists $view] || [llength [info commands $view]]} {
	    rtwizard_gui_smoke_fail "embedded ArcherCore left view $view behind"
	}
    }

    file delete -force $output
    puts "rtwizard historical GUI Obol workflow passed"
    exit 0
}

after 100 rtwizard_gui_smoke_probe
