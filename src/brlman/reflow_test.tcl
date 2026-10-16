#                    R E F L O W _ T E S T . T C L
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

if {[catch {package require Tk 8.6} message]} {
    puts "Tk tests skipped: $message"
    exit 77
}

proc text_y {widget marker} {
    set index [string first $marker [$widget text]]
    if {$index < 0} {
	error "Missing text marker '$marker'"
    }
    set box [$widget bbox $index]
    if {[llength $box] != 4} {
	error "Missing text box for '$marker'"
    }
    return [lindex $box 1]
}

proc text_height {widget first last} {
    return [expr {[text_y $widget $last] - [text_y $widget $first]}]
}

if {[catch {
    lassign $argv library generator
    load $library Tklitehtml

    set source {= REFLOW(1)
:doctype: manpage

== NAME
reflow - exercise manual page layout

== DESCRIPTION

}
    # Keep the paragraph long enough to wrap at every tested width.
    set paragraph_repetitions 100
    append source "PROSE_START [string repeat {Manual page content should use the available viewing width. } $paragraph_repetitions]PROSE_END\n\n"
    append source "----\nLITERAL_START\nLITERAL_END\n----\n"
    set html [exec $generator -q -b html5 -d manpage -o - << $source]

    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand 1
    set widget [$browser html]
    $widget load $html

    set window_height 600
    set window_widths {800 1200 1600 2200}
    wm maxsize . [lindex $window_widths end] $window_height
    set spans {}
    set literal_line_height {}
    # Fractional font positions may round to adjacent pixels on relayout.
    set pixel_rounding_tolerance 1
    foreach width $window_widths {
	wm geometry . ${width}x${window_height}
	update
	lappend spans [text_height $widget PROSE_START PROSE_END]
	if {[$widget xview] ne {0.0 1.0}} {
	    error "Manual page overflows horizontally at width $width"
	}
	set line_height [text_height $widget LITERAL_START LITERAL_END]
	if {$line_height <= 0 || ($literal_line_height ne "" &&
	    abs($line_height - $literal_line_height) > $pixel_rounding_tolerance)} {
	    error "Literal line breaks changed at width $width"
	}
	set literal_line_height $line_height
    }

    foreach previous [lrange $spans 0 end-1] current [lrange $spans 1 end] {
	if {$current >= $previous} {
	    error "Manual page stopped reflowing as its window widened: $spans"
	}
    }

    wm geometry . [lindex $window_widths 0]x${window_height}
    update
    if {[text_height $widget PROSE_START PROSE_END] != [lindex $spans 0]} {
	error "Manual page did not reflow back to its original width"
    }
} message options]} {
    puts stderr [dict get $options -errorinfo]
    exit 1
}

puts "Manual page reflow passed at window widths $window_widths"
exit 0

# Local Variables:
# mode: tcl
# tab-width: 8
# tcl-indent-level: 4
# indent-tabs-mode: t
# End:
# ex: shiftwidth=4 tabstop=8
