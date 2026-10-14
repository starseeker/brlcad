# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above
# copyright notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials provided
# with the distribution.
#
# 3. The name of the author may not be used to endorse or promote
# products derived from this software without specific prior written
# permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
# OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#

namespace eval ::tklitehtml {
    # Win32 defines one wheel detent as WHEEL_DELTA=120.
    variable win32WheelDelta 120
    variable browserReadChunkBytes 65536
    variable browserState
    variable browserCase
    variable browserQuery
    variable selectionState
    variable selectionScrollDelay 50
}


proc ::tklitehtml::_wheel {window axis delta} {
    variable win32WheelDelta

    if {$delta == 0} {
        return
    }

    switch -- [tk windowingsystem] {
        aqua {
            # Aqua reports small, accelerated deltas that already represent units.
            set units [expr {int(-$delta)}]
        }
        win32 {
            set units [expr {int(-$delta / $win32WheelDelta)}]
        }
        default {
            if {abs($delta) >= $win32WheelDelta} {
                set units [expr {int(-$delta / $win32WheelDelta)}]
            } else {
                set units [expr {$delta > 0 ? -1 : 1}]
            }
        }
    }

    if {$units == 0} {
        set units [expr {$delta > 0 ? -1 : 1}]
    }
    $window ${axis}view scroll $units units
}

bind TkLiteHtml <MouseWheel> {::tklitehtml::_wheel %W y %D; break}
bind TkLiteHtml <Shift-MouseWheel> {::tklitehtml::_wheel %W x %D; break}

if {[tk windowingsystem] eq "x11"} {
    bind TkLiteHtml <Button-4> {%W yview scroll -1 units; break}
    bind TkLiteHtml <Button-5> {%W yview scroll 1 units; break}
    bind TkLiteHtml <Shift-Button-4> {%W xview scroll -1 units; break}
    bind TkLiteHtml <Shift-Button-5> {%W xview scroll 1 units; break}
}

bind TkLiteHtml <Button-1> {::tklitehtml::_select_start %W %x %y character; break}
bind TkLiteHtml <Double-Button-1> {::tklitehtml::_select_start %W %x %y word; break}
bind TkLiteHtml <Triple-Button-1> {::tklitehtml::_select_start %W %x %y paragraph; break}
bind TkLiteHtml <Shift-Button-1> {::tklitehtml::_select_drag %W %x %y; break}
bind TkLiteHtml <B1-Motion> {::tklitehtml::_select_drag %W %x %y; break}
bind TkLiteHtml <ButtonRelease-1> {::tklitehtml::_select_stop %W; break}
bind TkLiteHtml <<Selection>> {::tklitehtml::_export_selection %W}
bind TkLiteHtml <<Copy>> {::tklitehtml::_copy %W; break}
bind TkLiteHtml <<SelectAll>> {%W selection set 0 end; break}
bind TkLiteHtml <Control-a> {%W selection set 0 end; break}
bind TkLiteHtml <Destroy> {::tklitehtml::_selection_destroy %W}
if {[tk windowingsystem] eq "aqua"} {
    bind TkLiteHtml <Command-a> {%W selection set 0 end; break}
}

proc ::tklitehtml::_selection_bounds {window index unit} {
    if {$unit eq "character"} {
        return [list $index $index]
    }
    set text [$window text]
    set index [expr {min($index, max(0, [string length $text] - 1))}]
    if {$unit eq "word"} {
        return [list [string wordstart $text $index] [string wordend $text $index]]
    }
    set first [expr {[string last \n $text [expr {$index - 1}]] + 1}]
    set last [string first \n $text $index]
    if {$last < 0} {
        set last [string length $text]
    } else {
        incr last
    }
    return [list $first $last]
}

proc ::tklitehtml::_select_start {window x y unit} {
    variable selectionState

    _select_stop $window
    focus $window
    set anchor [_selection_bounds $window [$window index @$x,$y] $unit]
    set selectionState($window) [dict create anchor $anchor unit $unit timer {} x $x y $y]
    $window selection set {*}$anchor
}

proc ::tklitehtml::_select_extend {window} {
    variable selectionState

    set state $selectionState($window)
    set x [expr {max(0, min([winfo width $window] - 1, [dict get $state x]))}]
    set y [expr {max(0, min([winfo height $window] - 1, [dict get $state y]))}]
    set index [$window index @$x,$y]
    lassign [dict get $state anchor] anchorFirst anchorLast
    lassign [_selection_bounds $window $index [dict get $state unit]] first last
    if {$index < $anchorFirst} {
        $window selection set $first $anchorLast
    } else {
        $window selection set $anchorFirst $last
    }
}

proc ::tklitehtml::_select_drag {window x y} {
    variable selectionState
    variable selectionScrollDelay

    if {![info exists selectionState($window)]} {
        _select_start $window $x $y character
    }
    dict set selectionState($window) x $x
    dict set selectionState($window) y $y
    _select_extend $window
    _select_stop $window
    if {$x < 0 || $y < 0 || $x >= [winfo width $window] || $y >= [winfo height $window]} {
        dict set selectionState($window) timer \
            [after $selectionScrollDelay [list ::tklitehtml::_select_scroll $window]]
    }
}

proc ::tklitehtml::_select_scroll {window} {
    variable selectionState

    dict set selectionState($window) timer {}
    set x [dict get $selectionState($window) x]
    set y [dict get $selectionState($window) y]
    foreach axis {x y} coordinate [list $x $y] size [list [winfo width $window] [winfo height $window]] {
        if {$coordinate < 0} {
            $window ${axis}view scroll -1 units
        } elseif {$coordinate >= $size} {
            $window ${axis}view scroll 1 units
        }
    }
    _select_drag $window $x $y
}

proc ::tklitehtml::_select_stop {window} {
    variable selectionState

    if {[info exists selectionState($window)]} {
        after cancel [dict get $selectionState($window) timer]
        dict set selectionState($window) timer {}
    }
}

proc ::tklitehtml::_selection_destroy {window} {
    variable selectionState

    _select_stop $window
    unset -nocomplain selectionState($window)
}

proc ::tklitehtml::_selection_chunk {window offset count} {
    return [string range [$window selection get] $offset [expr {$offset + $count - 1}]]
}

proc ::tklitehtml::_export_selection {window} {
    if {[tk windowingsystem] ne "x11"} {
        return
    }
    if {[$window selection range] eq ""} {
        if {[selection own -displayof $window] eq $window} {
            selection clear -displayof $window
        }
        return
    }
    # External clients need the target's wire format to preserve Unicode.
    foreach type {STRING UTF8_STRING} {
        selection handle -selection PRIMARY -type $type -format $type \
            $window [list ::tklitehtml::_selection_chunk $window]
    }
    selection own -selection PRIMARY -command \
        [list $window selection clear] $window
}

proc ::tklitehtml::_copy {window} {
    if {[$window selection range] ne ""} {
        clipboard clear -displayof $window
        clipboard append -displayof $window -- [$window selection get]
    }
}

bind TkLiteHtml <<DocumentChanged>> {::tklitehtml::_selection_destroy %W}
bind TkLiteHtmlFind <KeyPress> {break}

proc ::tklitehtml::_uri_encode {value {preserveSlash false}} {
    set encoded ""
    binary scan [encoding convertto utf-8 $value] c* bytes
    foreach byte $bytes {
        set byte [expr {$byte & 0xff}]
        set character [format %c $byte]
        if {($byte >= 0x30 && $byte <= 0x39) ||
            ($byte >= 0x41 && $byte <= 0x5a) ||
            ($byte >= 0x61 && $byte <= 0x7a) ||
            $character in {- . _ ~} ||
            ($preserveSlash && $character eq "/")} {
            append encoded $character
        } else {
            append encoded %[format %02X $byte]
        }
    }
    return $encoded
}

proc ::tklitehtml::_uri_decode {value} {
    set bytes ""
    set index 0
    while {$index < [string length $value]} {
        set escape [string range $value $index [expr {$index + 2}]]
        if {[regexp -- {^%([[:xdigit:]]{2})$} $escape match hex]} {
            append bytes [binary format H2 $hex]
            incr index 3
        } else {
            append bytes [string index $value $index]
            incr index
        }
    }
    if {[catch {encoding convertfrom utf-8 $bytes} decoded]} {
        return $bytes
    }
    return $decoded
}

proc ::tklitehtml::file_uri {path} {
    set normalized [string map {\\ /} [file normalize $path]]
    if {[string index $normalized 0] ne "/"} {
        set normalized /$normalized
    }
    return file://[_uri_encode $normalized true]
}

proc ::tklitehtml::local_path {location} {
    set location [lindex [split $location #] 0]
    set location [lindex [split $location ?] 0]
    if {![regexp -nocase {^file://([^/]*)(/.*)$} $location -> authority path]} {
        error "unsupported URI '$location'"
    }
    set path [_uri_decode $path]
    if {$authority ne "" && ![string equal -nocase $authority localhost]} {
        set path //$authority$path
    }
    if {$::tcl_platform(platform) eq "windows" && [regexp {^/[A-Za-z]:/} $path]} {
        set path [string range $path 1 end]
    }
    return [file normalize $path]
}

proc ::tklitehtml::browser {path args} {
    variable browserState

    if {[llength $args] % 2 != 0} {
        error "value for \"[lindex $args end]\" missing"
    }

    set linkCommand {}
    set widgetOptions {}
    foreach {option value} $args {
        if {$option eq "-linkcommand"} {
            set linkCommand $value
        } else {
            lappend widgetOptions $option $value
        }
    }

    frame $path -class TkLiteHtmlBrowser
    set hullCommand ${path}.__hull
    rename $path $hullCommand
    interp alias {} $path {} ::tklitehtml::_browser_dispatch $path
    set browserState($path) [dict create \
        hull $hullCommand \
        history {} \
        error {} \
        linkcommand $linkCommand \
        findquery {} \
        findcount 0 \
        findcurrent -1 \
        location {}]
    bind $path <Destroy> +[list ::tklitehtml::_browser_destroy $path %W]
    try {
        scrollbar $path.vscroll -orient vertical \
            -command [list $path.html yview]
        scrollbar $path.hscroll -orient horizontal \
            -command [list $path.html xview]
        label $path.error -anchor w -background #ffe8e8 \
            -foreground #800000 -justify left
        ::tklitehtml::widget $path.html {*}$widgetOptions \
            -xscrollcommand [list $path.hscroll set] \
            -yscrollcommand [list $path.vscroll set] \
            -linkcommand [list ::tklitehtml::_browser_link $path]
        _browser_controls $path

        pack $path.vscroll -side right -fill y
        pack $path.hscroll -side bottom -fill x
        pack $path.html -side left -fill both -expand true
    } on error {message options} {
        catch {destroy $path}
        return -options $options $message
    }
    return $path
}

proc ::tklitehtml::_browser_destroy {path destroyedWindow} {
    variable browserState
    variable browserCase
    variable browserQuery

    if {$path ne $destroyedWindow} {
        return
    }
    if {[info exists browserState($path)]} {
        set hullCommand [dict get $browserState($path) hull]
        trace remove variable browserQuery($path) write \
            [list ::tklitehtml::_browser_find_changed $path]
        if {[winfo exists $path.find.entry]} {
            $path.find.entry configure -textvariable {}
        }
        unset browserState($path)
        unset -nocomplain browserCase($path)
        unset -nocomplain browserQuery($path)
        catch {rename $hullCommand {}}
    }
    catch {rename $path {}}
}

proc ::tklitehtml::_browser_dispatch {path method args} {
    variable browserState

    if {![info exists browserState($path)]} {
        error "invalid browser \"$path\""
    }
    switch -- $method {
        back {
            if {[llength $args] != 0} {
                error "wrong # args: should be \"$path back\""
            }
            return [_browser_attempt $path \
                [list ::tklitehtml::_browser_back $path]]
        }
        cget {
            if {[llength $args] != 1} {
                error "wrong # args: should be \"$path cget option\""
            }
            if {[lindex $args 0] eq "-linkcommand"} {
                return [dict get $browserState($path) linkcommand]
            }
            return [$path.html cget {*}$args]
        }
        configure {
            return [_browser_configure $path {*}$args]
        }
        error {
            if {[llength $args] == 0} {
                return [dict get $browserState($path) error]
            }
            if {[llength $args] != 1} {
                error "wrong # args: should be \"$path error ?message?\""
            }
            _browser_show_error $path [lindex $args 0]
            return {}
        }
        findbar {
            if {[llength $args] > 1} {
                error "wrong # args: should be \"$path findbar ?pattern?\""
            }
            _browser_find_open $path {*}$args
            return $path.find.entry
        }
        goto {
            if {[llength $args] < 1 || [llength $args] > 2 ||
                ([llength $args] == 2 && [lindex $args 1] ne "-nosave")} {
                error "wrong # args: should be \"$path goto location ?-nosave?\""
            }
            return [_browser_attempt $path [list ::tklitehtml::_browser_goto \
                $path [lindex $args 0] [expr {[llength $args] == 1}]]]
        }
        html {
            if {[llength $args] != 0} {
                error "wrong # args: should be \"$path html\""
            }
            return $path.html
        }
        location {
            if {[llength $args] != 0} {
                error "wrong # args: should be \"$path location\""
            }
            return [dict get $browserState($path) location]
        }
        title {
            if {[llength $args] != 0} {
                error "wrong # args: should be \"$path title\""
            }
            return [$path.html title]
        }
        default {
            return [$path.html $method {*}$args]
        }
    }
}

proc ::tklitehtml::_browser_controls {path} {
    variable browserCase
    variable browserQuery

    set browserCase($path) false
    set browserQuery($path) {}
    frame $path.find
    label $path.find.label -text "Find:"
    entry $path.find.entry -width 20 -exportselection false \
        -textvariable ::tklitehtml::browserQuery($path)
    checkbutton $path.find.case -text "Match case" \
        -variable ::tklitehtml::browserCase($path) \
        -command [list ::tklitehtml::_browser_find_update $path true]
    button $path.find.previous -text "Previous" -state disabled \
        -command [list ::tklitehtml::_browser_find_step $path previous]
    button $path.find.next -text "Next" -state disabled \
        -command [list ::tklitehtml::_browser_find_step $path next]
    label $path.find.count -anchor w
    button $path.find.close -text "Close" \
        -command [list ::tklitehtml::_browser_find_hide $path]
    pack $path.find.label -side left
    pack $path.find.close $path.find.next $path.find.previous \
        $path.find.case $path.find.count -side right
    pack $path.find.entry -side left -fill x -expand true

    # Entry class bindings edit text before this tag stops application
    # shortcuts, such as the manual browser's single-key quit binding.
    bindtags $path.find.entry [linsert [bindtags $path.find.entry] 2 TkLiteHtmlFind]
    bind $path.find.entry <Return> "[list ::tklitehtml::_browser_find_step $path next]\nbreak"
    bind $path.find.entry <Shift-Return> "[list ::tklitehtml::_browser_find_step $path previous]\nbreak"
    bind $path.find.entry <Escape> "[list ::tklitehtml::_browser_find_hide $path]\nbreak"

    menu $path.menu -tearoff false
    $path.menu add command -label "Copy" -command [list ::tklitehtml::_copy $path.html]
    $path.menu add command -label "Select All" -command [list $path.html selection set 0 end]
    $path.menu add separator
    $path.menu add command -label "Find..." -accelerator "Ctrl+F" \
        -command [list ::tklitehtml::_browser_find_open $path]
    $path.menu add command -label "Find Next" -accelerator F3 \
        -command [list ::tklitehtml::_browser_find_step $path next]
    $path.menu add command -label "Find Previous" -accelerator "Shift+F3" \
        -command [list ::tklitehtml::_browser_find_step $path previous]
    $path.menu add command -label "Clear Highlights" \
        -command [list ::tklitehtml::_browser_find_clear $path]

    bind $path.html <Button-3> "[list ::tklitehtml::_browser_popup $path %X %Y]\nbreak"
    foreach sequence {<Control-f> <<Find>>} {
        bind $path.html $sequence "[list ::tklitehtml::_browser_find_open $path]\nbreak"
    }
    foreach window [list $path.html $path.find.entry] {
        bind $window <F3> "[list ::tklitehtml::_browser_find_step $path next]\nbreak"
        bind $window <Shift-F3> "[list ::tklitehtml::_browser_find_step $path previous]\nbreak"
    }
    bind $path.html <Escape> [list ::tklitehtml::_browser_find_escape $path]
    bind $path.html <<DocumentChanged>> [list ::tklitehtml::_browser_find_update $path true]
    if {[tk windowingsystem] eq "aqua"} {
        $path.menu entryconfigure "Find..." -accelerator "Command+F"
        foreach sequence {<Button-2> <Control-Button-1>} {
            bind $path.html $sequence "[list ::tklitehtml::_browser_popup $path %X %Y]\nbreak"
        }
        bind $path.html <Command-f> "[list ::tklitehtml::_browser_find_open $path]\nbreak"
    }
    trace add variable browserQuery($path) write \
        [list ::tklitehtml::_browser_find_changed $path]
}

proc ::tklitehtml::_browser_popup {path x y} {
    variable browserState

    set copyState [expr {[$path.html selection range] eq "" ? "disabled" : "normal"}]
    set findState [expr {[dict get $browserState($path) findcount] > 0 ? "normal" : "disabled"}]
    $path.menu entryconfigure Copy -state $copyState
    $path.menu entryconfigure "Find Next" -state $findState
    $path.menu entryconfigure "Find Previous" -state $findState
    tk_popup $path.menu $x $y
}

proc ::tklitehtml::_browser_find_open {path args} {
    variable browserQuery

    set selected [$path.html selection get]
    set query $browserQuery($path)
    if {[llength $args] != 0 || $selected ne ""} {
        set query [expr {[llength $args] != 0 ? [lindex $args 0] : $selected}]
    }
    pack $path.find -side top -fill x -before $path.html
    if {$query eq $browserQuery($path)} {
        _browser_find_update $path true
    } else {
        set browserQuery($path) $query
    }
    focus $path.find.entry
    $path.find.entry selection range 0 end
}

proc ::tklitehtml::_browser_find_changed {path args} {
    _browser_find_update $path
}

proc ::tklitehtml::_browser_find_update {path {force false}} {
    variable browserState
    variable browserCase
    variable browserQuery

    set query $browserQuery($path)
    if {!$force && $query eq [dict get $browserState($path) findquery]} {
        return
    }
    dict set browserState($path) findquery $query
    set count 0
    if {$query eq ""} {
        $path.html find clear
    } else {
        set options {}
        if {$browserCase($path)} {
            lappend options -case
        }
        set count [$path.html find set {*}$options $query]
    }
    dict set browserState($path) findcount $count
    dict set browserState($path) findcurrent [expr {$count > 0 ? 0 : -1}]
    _browser_find_status $path
}

proc ::tklitehtml::_browser_find_status {path} {
    variable browserState

    set count [dict get $browserState($path) findcount]
    set label {}
    if {$count > 0} {
        set label "[expr {[dict get $browserState($path) findcurrent] + 1}] of $count"
    } elseif {[dict get $browserState($path) findquery] ne ""} {
        set label "No matches"
    }
    $path.find.count configure -text $label
    set state [expr {$count > 0 ? "normal" : "disabled"}]
    $path.find.next configure -state $state
    $path.find.previous configure -state $state
}

proc ::tklitehtml::_browser_find_step {path direction} {
    variable browserState

    _browser_find_update $path
    if {[dict get $browserState($path) findquery] eq ""} {
        _browser_find_open $path
        return
    }
    set count [$path.html find $direction]
    dict set browserState($path) findcount $count
    if {$count > 0} {
        set step [expr {$direction eq "next" ? 1 : -1}]
        set current [dict get $browserState($path) findcurrent]
        dict set browserState($path) findcurrent [expr {($current + $step + $count) % $count}]
    }
    _browser_find_status $path
}

proc ::tklitehtml::_browser_find_clear {path} {
    $path.find.entry delete 0 end
    _browser_find_update $path true
}

proc ::tklitehtml::_browser_find_hide {path} {
    _browser_find_clear $path
    pack forget $path.find
    focus $path.html
}

proc ::tklitehtml::_browser_find_escape {path} {
    if {[winfo manager $path.find] ne ""} {
        _browser_find_hide $path
        return -code break
    }
}

proc ::tklitehtml::_browser_attempt {path command} {
    try {
        set result [uplevel #0 $command]
    } on error {message options} {
        _browser_show_error $path $message
        return -options $options $message
    }
    _browser_show_error $path {}
    return $result
}

proc ::tklitehtml::_browser_show_error {path message} {
    variable browserState

    if {![info exists browserState($path)]} {
        return
    }
    dict set browserState($path) error $message
    if {$message eq ""} {
        pack forget $path.error
        return
    }
    $path.error configure -text $message
    if {[winfo manager $path.error] eq ""} {
        pack $path.error -side top -fill x -before $path.html
    }
}

proc ::tklitehtml::_browser_read_document {filePath maximumBytes} {
    variable browserReadChunkBytes

    set channel [open $filePath r]
    try {
        fconfigure $channel -encoding binary -translation binary
        set bytes ""
        while {true} {
            set remaining [expr {$maximumBytes - [string length $bytes]}]
            if {$remaining == 0} {
                if {[read $channel 1] ne ""} {
                    error "HTML document '$filePath' exceeds the $maximumBytes byte limit"
                }
                break
            }
            set readBytes [expr {min($remaining, $browserReadChunkBytes)}]
            append bytes [read $channel $readBytes]
            if {[eof $channel]} {
                break
            }
        }
        return [encoding convertfrom utf-8 $bytes]
    } finally {
        close $channel
    }
}

proc ::tklitehtml::_browser_configure {path args} {
    variable browserState

    if {[llength $args] == 0} {
        return [concat [$path.html configure] \
            [list [list -linkcommand [dict get $browserState($path) linkcommand]]]]
    }
    if {[llength $args] == 1} {
        return [$path cget [lindex $args 0]]
    }
    if {[llength $args] % 2 != 0} {
        error "value for \"[lindex $args end]\" missing"
    }
    set widgetOptions {}
    foreach {option value} $args {
        if {$option eq "-linkcommand"} {
            dict set browserState($path) linkcommand $value
        } else {
            lappend widgetOptions $option $value
        }
    }
    if {[llength $widgetOptions] != 0} {
        $path.html configure {*}$widgetOptions
    }
    return {}
}

proc ::tklitehtml::_browser_split_location {path location} {
    variable browserState

    set fragment ""
    set separator [string first # $location]
    if {$separator >= 0} {
        set fragment [_uri_decode [string range $location [expr {$separator + 1}] end]]
        set location [string range $location 0 [expr {$separator - 1}]]
    }

    set current [dict get $browserState($path) location]
    if {$location eq ""} {
        if {$current eq ""} {
            error "a fragment cannot be resolved before a document is loaded"
        }
        set filePath [local_path $current]
    } elseif {[regexp {^[A-Za-z][A-Za-z0-9+.-]*:} $location]} {
        if {![string match -nocase file:* $location]} {
            return [list {} $fragment $location]
        }
        set filePath [local_path $location]
    } elseif {[file pathtype $location] eq "absolute"} {
        set filePath [file normalize [_uri_decode $location]]
    } else {
        set baseDirectory [pwd]
        if {$current ne ""} {
            set baseDirectory [file dirname [local_path $current]]
        }
        set filePath [file normalize [file join $baseDirectory [_uri_decode $location]]]
    }
    return [list $filePath $fragment {}]
}

proc ::tklitehtml::_browser_goto {path location saveHistory} {
    variable browserState

    lassign [_browser_split_location $path $location] filePath fragment externalUri
    if {$externalUri ne ""} {
        error "unsupported URI scheme in '$externalUri'"
    }

    set oldLocation [dict get $browserState($path) location]
    set sameDocument [expr {$oldLocation ne "" && [local_path $oldLocation] eq $filePath}]
    if {!$sameDocument} {
        if {![file readable $filePath] || ![file isfile $filePath]} {
            error "cannot read HTML document '$filePath'"
        }
        set maximumBytes [$path.html cget -maxdocumentbytes]
        if {[file size $filePath] > $maximumBytes} {
            error "HTML document '$filePath' exceeds the $maximumBytes byte limit"
        }
        set html [_browser_read_document $filePath $maximumBytes]
        $path.html load -baseurl [file_uri $filePath] $html
    }

    set newLocation [file_uri $filePath]
    if {$fragment ne ""} {
        append newLocation #[_uri_encode $fragment]
        $path.html fragment $fragment
    } elseif {$sameDocument} {
        $path.html yview moveto 0
    }
    if {!$sameDocument && $saveHistory && $oldLocation ne ""} {
        dict lappend browserState($path) history $oldLocation
    }
    dict set browserState($path) location $newLocation
    return $newLocation
}

proc ::tklitehtml::_browser_back {path} {
    variable browserState

    set history [dict get $browserState($path) history]
    if {[llength $history] == 0} {
        return {}
    }
    set location [lindex $history end]
    set result [_browser_goto $path $location false]
    dict set browserState($path) history [lrange $history 0 end-1]
    return $result
}

proc ::tklitehtml::_browser_link {path widget url attributes} {
    variable browserState

    if {![info exists browserState($path)]} {
        return
    }
    set command [dict get $browserState($path) linkcommand]
    if {[llength $command] != 0} {
        uplevel #0 [list {*}$command $path $url $attributes]
    } else {
        $path goto $url
    }
}
