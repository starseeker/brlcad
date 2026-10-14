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

    if {$path ne $destroyedWindow} {
        return
    }
    if {[info exists browserState($path)]} {
        set hullCommand [dict get $browserState($path) hull]
        unset browserState($path)
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
