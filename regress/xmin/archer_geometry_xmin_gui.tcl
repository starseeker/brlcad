#    A R C H E R _ G E O M E T R Y _ X M I N _ G U I . T C L
# BRL-CAD
#
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# version 2.1 as published by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this file; see the file named COPYING for more
# information.
#
###
if {![info exists ::env(GUI_TEST_LIBRARY)]} {
    puts stderr "GUI_TEST_LIBRARY is required"
    exit 2
}
source $::env(GUI_TEST_LIBRARY)
set ::env(ARCHER_PREFS_FILE) [file join $::env(GUI_TEST_DIR) .archerrc]

proc archer_geometry_preferences {{geometry ""}} {
    set channel [open $::env(ARCHER_PREFS_FILE) w]
    if {$geometry ne ""} {
	puts $channel [list set mWindowGeometry $geometry]
    }
    close $channel
}

proc archer_geometry_saved {} {
    set channel [open $::env(ARCHER_PREFS_FILE) r]
    set contents [read $channel]
    close $channel
    if {[regexp -line {^set mWindowGeometry (.+)$} $contents -> geometry]} {
	return $geometry
    }
    return ""
}

proc archer_geometry_settle {} {
    # Allow MapNotify and the deferred geometry restore to reach the event loop.
    set settle_delay_ms 100
    after $settle_delay_ms [list set ::archer_geometry_settled 1]
    vwait ::archer_geometry_settled
    update
    ::gui::test::require_no_background_errors
}

if {[catch {
    package require Archer 1.0
    ::gui::test::capture_background_errors
    wm withdraw .

    set saved_geometry 900x650+120+90
    set startup_geometry 1000x700+20+30
    set updated_geometry 950x680+60+70
    set default_geometry 920x660+40+50

    archer_geometry_preferences $saved_geometry
    Archer .unmapped
    ::gui::test::require {[wm state .unmapped] eq "withdrawn"} \
	"startup fixture unexpectedly showed its window"
    ::itcl::delete object .unmapped
    archer_geometry_settle
    ::gui::test::require {[archer_geometry_saved] eq $saved_geometry} \
	"closing an unmapped Archer replaced its saved geometry"

    archer_geometry_preferences
    Archer .unmapped_new
    ::itcl::delete object .unmapped_new
    archer_geometry_settle
    ::gui::test::require {[archer_geometry_saved] eq ""} \
	"closing an unmapped Archer saved an unfinished layout size"

    archer_geometry_preferences $saved_geometry
    Archer .restore
    .restore configure -geometry $startup_geometry
    wm geometry .restore $startup_geometry
    event generate [.restore component hpane] <Map>
    archer_geometry_settle
    ::gui::test::require {[wm geometry .restore] eq $startup_geometry} \
	"a child Map event restored the main window's geometry prematurely"

    wm deiconify .restore
    archer_geometry_settle
    ::gui::test::require {[wm geometry .restore] eq $saved_geometry} \
	"mapping the main window did not restore its saved geometry"

    wm geometry .restore $updated_geometry
    event generate [.restore component hpane] <Map>
    archer_geometry_settle
    ::gui::test::require {[wm geometry .restore] eq $updated_geometry} \
	"a later child Map event reset the resized main window"
    wm withdraw .restore
    ::itcl::delete object .restore
    archer_geometry_settle
    ::gui::test::require {[archer_geometry_saved] eq $updated_geometry} \
	"closing a previously mapped window did not save its latest geometry"

    archer_geometry_preferences
    Archer .default
    wm geometry .default $default_geometry
    wm deiconify .default
    archer_geometry_settle
    ::gui::test::require {[wm geometry .default] eq $default_geometry} \
	"an empty geometry option reset the window to its layout size"
    wm withdraw .default
    ::itcl::delete object .default
    archer_geometry_settle
} message options]} {
    puts stderr $message
    if {[dict exists $options -errorinfo]} {
	puts stderr [dict get $options -errorinfo]
    }
    puts "FAIL: Archer window geometry regression"
    exit 1
}

puts "PASS: Archer window geometry restoration and startup cancellation"
exit 0

# Local Variables:
# tab-width: 8
# mode: Tcl
# tcl-indent-level: 4
# indent-tabs-mode: t
# End:
# ex: shiftwidth=4 tabstop=8
