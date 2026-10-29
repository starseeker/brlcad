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

if {![info exists ::tklitehtml_utf8_test_suite]} {
    set ::tklitehtml_utf8_test_suite 1
    source -encoding utf-8 [info script]
    return
}

package require Tcl 8.6

if {[catch {package require Tk 8.6} message]} {
    puts "Tk tests skipped: $message"
    exit 77
}

if {$argv eq {-selectionpeer}} {
    after 10000 {exit 1}
    text .paste
    pack .paste
    update
    event generate .paste <ButtonPress-2> -x 5 -y 5
    event generate .paste <ButtonRelease-2> -x 5 -y 5
    update
    set primary [.paste get 1.0 end-1c]
    .paste delete 1.0 end
    event generate .paste <<Paste>>
    update
    set clipboard [.paste get 1.0 end-1c]
    puts [list [string length $primary] \
        [zlib crc32 [encoding convertto utf-8 $primary]] \
        [string length $clipboard] \
        [zlib crc32 [encoding convertto utf-8 $clipboard]]]
    exit
}

package require tcltest 2
namespace import ::tcltest::test
::tcltest::configure -verbose {start error}
::tcltest::testConstraint x11 [expr {[tk windowingsystem] eq "x11"}]
set permissionProbe [file join [::tcltest::temporaryDirectory] \
    tklitehtml-permission-probe-[pid]]
set permissionChannel [open $permissionProbe w]
close $permissionChannel
if {$::tcl_platform(platform) eq "unix"} {
    file attributes $permissionProbe -permissions 0000
}
::tcltest::testConstraint filePermissions \
    [expr {$::tcl_platform(platform) eq "unix" && \
        ![file readable $permissionProbe]}]
catch {file attributes $permissionProbe -permissions 0600}
file delete -force $permissionProbe

set packageIndex [lsearch -exact $argv -packagepath]
if {$packageIndex < 0 || $packageIndex + 1 >= [llength $argv]} {
    error "usage: all.tcl -packagepath directory"
}
lappend auto_path [file normalize [lindex $argv [expr {$packageIndex + 1}]]]
package require tklitehtml 0.1.0

proc makeWidget {args} {
    set widget [::tklitehtml::widget .html {*}$args]
    pack $widget -fill both -expand 1
    update
    return $widget
}

proc deleteWidget {} {
    if {[winfo exists .html]} {
        destroy .html
        update
    }
}

test package-1.1 {the extension provides its package and factory command} -body {
    list [package present tklitehtml] [namespace which ::tklitehtml::widget]
} -result {0.1.0 ::tklitehtml::widget}

test widget-1.1 {widget options support configure and cget} -setup {
    set widget [makeWidget -width 320 -height 180 -background #f8f8f8]
} -body {
    $widget configure -width 360 -height 200
    list [$widget cget -width] [$widget cget -height] [$widget cget -background]
} -cleanup {
    deleteWidget
} -result {360 200 #f8f8f8}

test widget-1.2 {resource limits are configurable and reject invalid values} -setup {
    set widget [makeWidget -maxdocumentbytes 4096 \
        -maximagepixels 2048 -maxviewportpixels 1048576]
} -body {
    list \
        [$widget cget -maxdocumentbytes] \
        [$widget cget -maximagepixels] \
        [$widget cget -maxviewportpixels] \
        [catch {$widget configure -maxdocumentbytes 0}]
} -cleanup {
    deleteWidget
} -result {4096 2048 1048576 1}

test widget-1.3 {constructor errors leave no partial widget command} -body {
    set status [catch {
        ::tklitehtml::widget .invalid -maxdocumentbytes 0
    }]
    list $status \
        [winfo exists .invalid] \
        [llength [info commands .invalid]]
} -cleanup {
    catch {destroy .invalid}
} -result {1 0 0}

test document-1.1 {load exposes rendered Unicode text and character ranges} -setup {
    set widget [makeWidget]
    $widget load {<p>Hello Café world.</p><p>Hello café.</p>}
    update
} -body {
    set text [$widget text]
    set expectedFirst [string first Café $text]
    set expectedSecond [string first café $text [expr {$expectedFirst + 1}]]
    set ranges [$widget search -all -nocase café]
    list \
        [expr {[llength $ranges] == 2}] \
        [expr {[lindex $ranges 0 0] == $expectedFirst}] \
        [expr {[lindex $ranges 1 0] == $expectedSecond}] \
        [expr {[lindex $ranges 0 1] - [lindex $ranges 0 0] == 4}]
} -cleanup {
    deleteWidget
} -result {1 1 1 1}

test document-1.2 {user styles render and character boxes are queryable} -setup {
    set widget [makeWidget]
    $widget load -userstyle {body { font-size: 18px; } strong { color: #124578; }} \
        {<p>plain <strong>styled</strong> text</p>}
    update
} -body {
    set index [string first styled [$widget text]]
    lassign [$widget bbox $index] x y width height
    expr {$x >= 0 && $y >= 0 && $width > 0 && $height > 0}
} -cleanup {
    deleteWidget
} -result 1

test document-1.3 {title and fragment expose document navigation metadata} -setup {
    set widget [makeWidget -height 80]
    $widget load {<title>Navigation test</title><div style="height: 900px">top</div><h2 id="target">Target</h2>}
    update
} -body {
    $widget fragment target
    update
    list [$widget title] [expr {[lindex [$widget yview] 0] > 0.0}]
} -cleanup {
    deleteWidget
} -result {{Navigation test} 1}

test document-1.4 {text extraction preserves visual word and line boundaries} -setup {
    set widget [makeWidget]
    $widget load {<h2><span>Modeling</span> With CSG</h2><p>Second line</p>}
    update
} -body {
    set text [$widget text]
    list \
        [expr {[string first {Modeling With CSG} $text] >= 0}] \
        [expr {[string first "\n" $text] >= 0}]
} -cleanup {
    deleteWidget
} -result {1 1}

test document-1.5 {an oversized load preserves the current document} -setup {
    set widget [makeWidget -maxdocumentbytes 128]
    $widget load {<title>Original</title><p>Original text</p>}
} -body {
    set status [catch {$widget load [string repeat x 129]} message]
    list $status [string match {*exceeds the 128 byte limit} $message] \
        [$widget title] \
        [expr {[string first {Original text} [$widget text]] >= 0}]
} -cleanup {
    deleteWidget
} -result {1 1 Original 1}

test document-1.6 {exposures during page replacement and clear remain safe} -setup {
    set paintingErrors {}
    set previousBgerror [interp bgerror {}]
    interp bgerror {} {apply {{message options} {
        lappend ::paintingErrors $message
    }}}
    set widget [makeWidget]
    $widget load {<p style="font: 48px serif">Original frame</p>}
    update
} -body {
    set replacements 20
    for {set iteration 1} {$iteration <= $replacements} {incr iteration} {
        # An exposure can repaint the old frame before a queued page update.
        event generate $widget <Expose>
        $widget load "<title>Replacement $iteration</title><p style='font: 24px monospace'>Replacement frame</p>"
        update
        event generate $widget <Expose>
        $widget clear
        update
    }
    list [$widget text] [$widget title] $paintingErrors
} -cleanup {
    deleteWidget
    interp bgerror {} $previousBgerror
    unset paintingErrors previousBgerror
} -result {{} {} {}}

test search-1.1 {search option terminator permits a leading hyphen} -setup {
    set widget [makeWidget]
    $widget load {<p>-needle and needle</p>}
} -body {
    expr {[llength [$widget search -all -- -needle]] == 1}
} -cleanup {
    deleteWidget
} -result 1

test reflow-1.1 {logical text and find ranges survive wrapping and resizing} -setup {
    set originalGeometry [wm geometry .]
    wm geometry . 760x300
    set widget [makeWidget]
    $widget load {<h2>Heading</h2><p>alpha <b>Café</b> beta gamma delta epsilon zeta eta theta iota kappa lambda END</p><p>Last paragraph</p>}
    update
} -body {
    set text [$widget text]
    set expected {Heading
alpha Café beta gamma delta epsilon zeta eta theta iota kappa lambda END
Last paragraph}
    set matches [$widget search -all {lambda END}]
    $widget find set {lambda END}
    set current [$widget find current]
    set endIndex [string first END $text]
    set wideBox [$widget bbox $endIndex]
    set stable [expr {$text eq $expected && [llength $matches] == 1}]
    foreach width {140 420 760 140} {
        wm geometry . ${width}x300
        update
        if {[$widget text] ne $text ||
            [$widget search -all {lambda END}] ne $matches ||
            [$widget find current] ne $current} {
            set stable 0
        }
    }
    list $stable [expr {[lindex [$widget bbox $endIndex] 1] > [lindex $wideBox 1]}]
} -cleanup {
    deleteWidget
    wm geometry . $originalGeometry
    update
} -result {1 1}

test reflow-1.2 {structural breaks and literal whitespace are independent of layout} -setup {
    set widget [makeWidget]
    $widget load {<p>one<span>two</span>   three<br>four</p><pre>first
  second

last</pre><p style="display:none">hidden</p><p>end</p>}
    update
} -body {
    $widget text
} -cleanup {
    deleteWidget
} -result {onetwo three
four
first
  second

last
end}

test highlight-1.1 {named highlights accept, configure, and clear ranges} -setup {
    set widget [makeWidget]
    $widget load {<p>alpha beta alpha</p>}
    set ranges [$widget search -all alpha]
} -body {
    $widget highlight configure selection -background #80c0ff
    $widget highlight set selection $ranges
    update
    $widget highlight clear selection
    update
    set result ok
} -cleanup {
    deleteWidget
} -result ok

test find-1.1 {find tracks all matches and cycles the current match} -setup {
    set widget [makeWidget]
    $widget load {<p>one TWO one two</p>}
} -body {
    set count [$widget find set two]
    set first [$widget find current]
    $widget find next
    set second [$widget find current]
    $widget find previous
    list $count [expr {$first ne $second}] [expr {[$widget find current] eq $first}]
} -cleanup {
    deleteWidget
} -result {2 1 1}

test selection-1.1 {selection uses Unicode ranges and normalizes reverse drags} -setup {
    set widget [makeWidget]
    $widget load {<p>alpha Café omega</p>}
    update
} -body {
    $widget selection set 10 6
    set selected [list [$widget selection range] [$widget selection get]]
    $widget selection clear
    list $selected [$widget selection range] [$widget selection get] \
        [catch {$widget selection set -1 end}] \
        [catch {$widget index @bad,3}]
} -cleanup {
    deleteWidget
} -result {{{6 10} Café} {} {} 1 1}

test selection-1.2 {coordinate indexes follow character boxes and scrolling} -setup {
    set widget [makeWidget -height 100]
    $widget load {<div style="height:900px">top</div><p>Café target</p>}
    update
} -body {
    set index [string first target [$widget text]]
    $widget see $index
    update
    lassign [$widget bbox $index] x y width height
    list [expr {[$widget index @$x,$y] == $index}] \
        [expr {[$widget index @[expr {$x + $width}],$y] == $index + 1}]
} -cleanup {
    deleteWidget
} -result {1 1}

test selection-1.3 {selection survives reflow and find navigation} -setup {
    set originalGeometry [wm geometry .]
    wm geometry . 700x300
    set widget [makeWidget]
    $widget load {<p>one two three four five six seven eight nine ten Café target</p>}
    update
} -body {
    set first [string first Café [$widget text]]
    $widget selection set $first end
    $widget find set three
    wm geometry . 160x300
    update
    $widget find next
    update
    list [$widget selection get] [$widget selection range]
} -cleanup {
    deleteWidget
    wm geometry . $originalGeometry
    update
} -result {{Café target} {49 60}}

test selection-1.4 {Copy preserves CLIPBOARD through later selection and navigation} -setup {
    set widget [makeWidget -maxdocumentbytes 128]
    $widget load {<p>alpha Café omega</p>}
    update
} -body {
    $widget selection set 6 10
    event generate $widget <<Copy>>
    update
    $widget selection set 0 5
    update
    set copied [clipboard get -displayof $widget]
    set failed [catch {$widget load [string repeat x 129]}]
    set preserved [$widget selection get]
    $widget load {<p>replacement</p>}
    update
    list $copied $failed $preserved [$widget selection range] \
        [clipboard get -displayof $widget]
} -cleanup {
    deleteWidget
} -result {Café 1 alpha {} Café}

test selection-1.5 {X11 PRIMARY follows selection independently of CLIPBOARD} -constraints x11 -setup {
    set widget [makeWidget]
    $widget load {<p>alpha Café omega</p>}
    update
    clipboard clear -displayof $widget
    clipboard append -displayof $widget sentinel
} -body {
    $widget selection set 6 10
    update
    set primary [selection get -displayof $widget -selection PRIMARY -type UTF8_STRING]
    set other [::tklitehtml::widget .other -height 30]
    pack $other -before $widget
    $other load {<p>other</p>}
    update
    $other selection set 0 end
    update
    list $primary [clipboard get -displayof $widget] \
        [$widget selection range] \
        [selection get -displayof $other -selection PRIMARY -type UTF8_STRING]
} -cleanup {
    catch {destroy .other}
    deleteWidget
} -result {Café sentinel {} other}

proc recordSelectionLink {args} {
    set ::clickedLink $args
}

test selection-1.6 {mouse drags select link text without following the link} -setup {
    set ::clickedLink {}
    set widget [makeWidget -linkcommand recordSelectionLink]
    $widget load {<p><a href="#target">Café target</a></p><p id="target">destination</p>}
    update
} -body {
    lassign [$widget bbox 0] x y width height
    lassign [$widget bbox 3] lastX lastY lastWidth lastHeight
    incr lastX $lastWidth
    event generate $widget <Motion> -x $x -y $y
    event generate $widget <ButtonPress-1> -x $x -y $y
    event generate $widget <Motion> -state 256 -x $lastX -y $lastY
    event generate $widget <ButtonRelease-1> -x $lastX -y $lastY
    update
    list [$widget selection get] $::clickedLink
} -cleanup {
    unset -nocomplain ::clickedLink
    deleteWidget
} -result {Café {}}

test selection-1.7 {word and paragraph selection respect logical boundaries} -setup {
    set widget [makeWidget]
    $widget load {<p>alpha Café omega</p><p>second paragraph</p>}
    update
} -body {
    lassign [$widget bbox 7] x y width height
    ::tklitehtml::_select_start $widget $x $y word
    set word [$widget selection get]
    ::tklitehtml::_select_start $widget $x $y paragraph
    list $word [$widget selection get]
} -cleanup {
    deleteWidget
} -result {Café {alpha Café omega
}}

test selection-1.8 {destroying a widget cancels pending drag autoscrolling} -setup {
    set widget [makeWidget]
    $widget load {<p>drag</p>}
    update
} -body {
    ::tklitehtml::_select_start $widget 8 8 character
    ::tklitehtml::_select_drag $widget 8 -10
    set timer [dict get $::tklitehtml::selectionState($widget) timer]
    destroy $widget
    list [expr {$timer ni [after info]}] \
        [info exists ::tklitehtml::selectionState($widget)]
} -cleanup {
    deleteWidget
} -result {1 0}

test selection-1.9 {another process can paste long Unicode text after reflow} -constraints x11 -setup {
    set originalGeometry [wm geometry .]
    set widget [makeWidget]
    # Exceed a selection-transfer chunk to exercise multibyte boundaries.
    set expected [string trimright [string repeat {Café λ漢字 } 1200]]
    $widget load "<p>$expected</p>"
    update
    $widget selection set 0 end
    wm geometry . 240x200
    update
    event generate $widget <<Copy>>
    set peer {}
    set timer {}
    set ::selectionPeerResult {}
} -body {
    set peer [open [list | [info nameofexecutable] \
        [file normalize [info script]] -selectionpeer] r]
    fconfigure $peer -blocking false
    fileevent $peer readable [list apply {{channel} {
        if {[gets $channel line] >= 0} {
            set ::selectionPeerResult $line
        } elseif {[eof $channel]} {
            set ::selectionPeerResult eof
        }
    }} $peer]
    set timer [after 15000 {set ::selectionPeerResult timeout}]
    vwait ::selectionPeerResult
    set metrics [list [string length $expected] \
        [zlib crc32 [encoding convertto utf-8 $expected]]]
    if {$::selectionPeerResult ne [concat $metrics $metrics]} {
        error "paste metrics {$::selectionPeerResult}, expected {[concat $metrics $metrics]}"
    }
    set result 1
} -cleanup {
    after cancel $timer
    if {$peer ne ""} {
        fileevent $peer readable {}
        close $peer
    }
    unset -nocomplain ::selectionPeerResult
    deleteWidget
    wm geometry . $originalGeometry
    update
} -result 1

proc recordScroll {first last} {
    set ::scrollFractions [list $first $last]
}

test scrolling-1.1 {standard view commands and scrollbar callbacks stay synchronized} -setup {
    set ::scrollFractions {}
    set widget [makeWidget -height 100 -yscrollcommand recordScroll]
    set paragraphs [string repeat {<p>scrolling content</p>} 80]
    $widget load $paragraphs
    update
} -body {
    $widget yview moveto 0.5
    update
    lassign [$widget yview] first last
    list [expr {$first > 0.0}] [expr {$last < 1.0}] [expr {[llength $::scrollFractions] == 2}]
} -cleanup {
    unset -nocomplain ::scrollFractions
    deleteWidget
} -result {1 1 1}

test scrolling-1.2 {X11 wheel buttons use portable class bindings} -constraints x11 -setup {
    set widget [makeWidget -width 100 -height 100]
    $widget load "<div style='width: 2000px'>wide</div>[string repeat {<p>wheel content</p>} 80]"
    update
} -body {
    set yBefore [lindex [$widget yview] 0]
    event generate $widget <Button-5> -x 10 -y 10
    update
    set xBefore [lindex [$widget xview] 0]
    event generate $widget <Shift-Button-5> -x 10 -y 10
    update
    list \
        [expr {[lindex [$widget yview] 0] > $yBefore}] \
        [expr {[lindex [$widget xview] 0] > $xBefore}] \
        [expr {[bind TkLiteHtml <MouseWheel>] ne ""}]
} -cleanup {
    deleteWidget
} -result {1 1 1}

proc testResource {widget kind url} {
    lappend ::resourceCalls [list $kind $url]
    if {$kind eq "image"} {
        return ::testImage
    }
    if {$kind eq "stylesheet"} {
        return {p { color: #2468ac; }}
    }
    return {}
}

test resources-1.1 {resource callbacks supply images and imported stylesheets} -setup {
    set ::resourceCalls {}
    image create photo ::testImage -width 3 -height 2
    ::testImage put #33aa55 -to 0 0 3 2
    set widget [makeWidget -resourcecommand testResource]
} -body {
    $widget load -baseurl https://example.invalid/docs/page.html {
        <link rel="stylesheet" href="theme.css">
        <p>resource</p><img src="picture.png" width="9" height="6">
    }
    update
    set kinds [lmap call $::resourceCalls {lindex $call 0}]
    list [expr {"stylesheet" in $kinds}] [expr {"image" in $kinds}]
} -cleanup {
    deleteWidget
    image delete ::testImage
    unset -nocomplain ::resourceCalls
} -result {1 1}

proc failingResource {widget kind url} {
    error "fixture resource failure"
}

test resources-1.2 {a callback error preserves the current document} -setup {
    set widget [makeWidget -resourcecommand failingResource]
    $widget load {<title>Original</title><p>Original resource page</p>}
} -body {
    set status [catch {
        $widget load {<link rel="stylesheet" href="bad.css"><p>Replacement</p>}
    } message]
    list $status \
        [string match {*fixture resource failure*} $message] \
        [$widget title] \
        [expr {[string first {Original resource page} [$widget text]] >= 0}]
} -cleanup {
    deleteWidget
} -result {1 1 Original 1}

proc destroyingResource {widget kind url} {
    destroy $widget
    return {}
}

test resources-1.3 {destroying a widget from a resource callback is safe} -setup {
    set widget [makeWidget -resourcecommand destroyingResource]
} -body {
    set status [catch {
        $widget load {<link rel="stylesheet" href="destroy.css"><p>gone</p>}
    } message]
    list $status [string match {*destroyed*} $message] \
        [winfo exists .html] \
        [llength [info commands .html]]
} -cleanup {
    deleteWidget
} -result {1 1 0 0}

proc oversizedImageResource {widget kind url} {
    if {$kind eq "image"} {
        return ::oversizedImageFixture
    }
    return {}
}

test resources-1.4 {an oversized image is rejected transactionally} -setup {
    image create photo ::oversizedImageFixture -width 20 -height 20
    set widget [makeWidget -maximagepixels 100 \
        -resourcecommand oversizedImageResource]
    $widget load {<title>Original</title><p>Original image page</p>}
} -body {
    set status [catch {$widget load {<img src="large.png">}} message]
    list $status [string match {*exceeds the 100 pixel limit} $message] \
        [$widget title]
} -cleanup {
    deleteWidget
    image delete ::oversizedImageFixture
} -result {1 1 Original}

test resources-1.5 {missing local resources do not destabilize rendering} -setup {
    set widget [makeWidget]
} -body {
    set status [catch {
        $widget load -baseurl /definitely/missing/index.html {
            <link rel="stylesheet" href="missing.css">
            <p>Readable text</p><img src="missing.png">
        }
    }]
    list $status [expr {[string first {Readable text} [$widget text]] >= 0}]
} -cleanup {
    deleteWidget
} -result {0 1}

proc reentrantResource {widget kind url} {
    $widget clear
    return {}
}

test resources-1.6 {reentrant widget mutation is rejected transactionally} -setup {
    set widget [makeWidget -resourcecommand reentrantResource]
    $widget load {<title>Original</title><p>Original page</p>}
} -body {
    set status [catch {
        $widget load {<link rel="stylesheet" href="recursive.css">}
    } message]
    list $status \
        [string match {*unavailable during document loading*} $message] \
        [$widget title]
} -cleanup {
    deleteWidget
} -result {1 1 Original}

test resources-1.7 {repeated image loads release internal Tk images} -setup {
    set imagePath [file join [::tcltest::temporaryDirectory] \
        tklitehtml-image-[pid].ppm]
    image create photo ::imageFixture -width 4 -height 3
    ::imageFixture put #336699 -to 0 0 4 3
    ::imageFixture write $imagePath -format ppm
    image delete ::imageFixture
    set initialImages [lsort [image names]]
} -body {
    for {set iteration 0} {$iteration < 20} {incr iteration} {
        set widget [makeWidget]
        $widget load "<img src='[::tklitehtml::file_uri $imagePath]'>"
        update
        deleteWidget
    }
    expr {[lsort [image names]] eq $initialImages}
} -cleanup {
    deleteWidget
    file delete -force $imagePath
} -result 1

test resources-1.8 {oversized local image files preserve the document} -setup {
    set imagePath [file join [::tcltest::temporaryDirectory] \
        tklitehtml-large-image-[pid].ppm]
    image create photo ::largeImageFixture -width 20 -height 20
    ::largeImageFixture put #336699 -to 0 0 20 20
    ::largeImageFixture write $imagePath -format ppm
    image delete ::largeImageFixture
    set widget [makeWidget -maxdocumentbytes 128]
    $widget load {<title>Original</title><p>Original page</p>}
} -body {
    set baseUrl [::tklitehtml::file_uri \
        [file join [file dirname $imagePath] index.html]]
    set imageName [file tail $imagePath]
    set status [catch {
        $widget load -baseurl $baseUrl "<img src='$imageName'>"
    } message]
    list $status \
        [string match {*image resource exceeds*} $message] \
        [$widget title]
} -cleanup {
    deleteWidget
    file delete -force $imagePath
} -result {1 1 Original}

proc imageBudgetResource {widget kind url} {
    if {$kind ne "image"} {
        return {}
    }
    if {[string match {*first.png} $url]} {
        return ::firstImageFixture
    }
    return ::secondImageFixture
}

test resources-1.9 {the decoded image budget applies to the document} -setup {
    image create photo ::firstImageFixture -width 8 -height 8
    image create photo ::secondImageFixture -width 8 -height 8
    set widget [makeWidget -maximagepixels 100 \
        -resourcecommand imageBudgetResource]
    $widget load {<title>Original</title><p>Original page</p>}
} -body {
    set status [catch {
        $widget load {<img src="first.png"><img src="second.png">}
    } message]
    list $status \
        [string match {*exceeds the 100 pixel limit} $message] \
        [$widget title]
} -cleanup {
    deleteWidget
    image delete ::firstImageFixture ::secondImageFixture
} -result {1 1 Original}

proc oversizedStylesheet {widget kind url} {
    if {$kind eq "stylesheet"} {
        return [string repeat {p { color: #336699; }} 10]
    }
    return {}
}

test resources-1.10 {the byte budget includes imported stylesheets} -setup {
    set widget [makeWidget -maxdocumentbytes 128 \
        -resourcecommand oversizedStylesheet]
    $widget load {<title>Original</title><p>Original page</p>}
} -body {
    set status [catch {
        $widget load {<link rel="stylesheet" href="large.css">}
    } message]
    list $status \
        [string match {*stylesheet exceeds*} $message] \
        [$widget title]
} -cleanup {
    deleteWidget
} -result {1 1 Original}

proc recordLink {widget url attributes} {
    set ::clickedLink [list $widget $url $attributes]
}

test links-1.1 {pointer activation invokes the link callback with metadata} -setup {
    unset -nocomplain ::clickedLink
    set widget [makeWidget -linkcommand recordLink]
    $widget load -baseurl https://example.invalid/docs/page.html {
        <p><a href="next.html" title="Next page" id="next-link">Open next</a></p>
    }
    update
} -body {
    set index [string first {Open next} [$widget text]]
    lassign [$widget bbox $index] x y width height
    event generate $widget <Motion> -x [expr {$x + 1}] -y [expr {$y + 1}]
    event generate $widget <ButtonPress-1> -x [expr {$x + 1}] -y [expr {$y + 1}]
    event generate $widget <ButtonRelease-1> -x [expr {$x + 1}] -y [expr {$y + 1}]
    update
    list \
        [lindex $::clickedLink 0] \
        [lindex $::clickedLink 1] \
        [dict get [lindex $::clickedLink 2] title] \
        [dict get [lindex $::clickedLink 2] id]
} -cleanup {
    unset -nocomplain ::clickedLink
    deleteWidget
} -result {.html https://example.invalid/docs/next.html {Next page} next-link}

test links-1.2 {selection and context-menu modifier clicks do not activate links} -setup {
    set ::clickedLink {}
    set widget [makeWidget -linkcommand recordLink]
    $widget load {<p><a href="#target">link</a></p><p id="target">destination</p>}
    update
} -body {
    lassign [$widget bbox 0] x y width height
    event generate $widget <Motion> -x $x -y $y
    foreach modifier {Control Shift} {
        event generate $widget <$modifier-ButtonPress-1> -x $x -y $y
        event generate $widget <$modifier-ButtonRelease-1> -x $x -y $y
    }
    update
    set ::clickedLink
} -cleanup {
    unset -nocomplain ::clickedLink
    deleteWidget
} -result {}

test browser-1.1 {local browser follows relative links, fragments, and history} -setup {
    set testDirectory [file normalize [file join [::tcltest::temporaryDirectory] browser-test-[pid]]]
    file mkdir $testDirectory
    set firstPath [file join $testDirectory {first page.html}]
    set secondPath [file join $testDirectory second.html]
    set channel [open $firstPath w]
    puts $channel {<title>First</title><p>First page</p>}
    close $channel
    set channel [open $secondPath w]
    puts $channel {<title>Second</title><div style="height: 900px">Second page</div><h2 id="last">Last</h2>}
    close $channel
    set browser [::tklitehtml::browser .browser -width 240 -height 80]
    pack $browser -fill both -expand 1
    $browser goto $firstPath
    update
} -body {
    $browser goto second.html#last
    update
    set second [list [$browser title] [expr {[lindex [[$browser html] yview] 0] > 0.0}]]
    $browser back
    update
    list \
        $second \
        [$browser title] \
        [expr {[::tklitehtml::local_path [$browser location]] eq $firstPath}] \
        [expr {[winfo exists [$browser html]]}]
} -cleanup {
    catch {destroy .browser}
    update
    file delete -force $testDirectory
} -result {{Second 1} First 1 1}

test browser-1.2 {file URI conversion preserves Unicode and reserved path characters} -body {
    set path [file normalize [file join [::tcltest::temporaryDirectory] {manual # é.html}]]
    expr {[::tklitehtml::local_path [::tklitehtml::file_uri $path]] eq $path}
} -result 1

test browser-1.3 {constructor errors leave no partial browser commands} -body {
    set status [catch {
        ::tklitehtml::browser .broken -maxdocumentbytes 0
    }]
    list $status \
        [winfo exists .broken] \
        [llength [info commands .broken]] \
        [llength [info commands .broken.__hull]]
} -cleanup {
    catch {destroy .broken}
} -result {1 0 0 0}

test browser-1.4 {failed navigation preserves the current document} -setup {
    set testDirectory [file normalize [file join \
        [::tcltest::temporaryDirectory] browser-errors-[pid]]]
    file mkdir $testDirectory
    set firstPath [file join $testDirectory first.html]
    set channel [open $firstPath w]
    puts $channel {<title>First</title><p>First page</p>}
    close $channel
    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand 1
    $browser goto $firstPath
} -body {
    set originalLocation [$browser location]
    set statuses {}
    foreach target [list \
        [file join $testDirectory missing.html] \
        $testDirectory \
        https://example.invalid/help.html] {
        lappend statuses [catch {$browser goto $target}]
    }
    set errorVisible [expr {[$browser error] ne "" &&
        [winfo manager .browser.error] eq "pack"}]
    $browser goto $originalLocation
    list $statuses \
        [expr {[$browser location] eq $originalLocation}] \
        [$browser title] \
        $errorVisible \
        [expr {[$browser error] eq ""}] \
        [winfo manager .browser.error]
} -cleanup {
    catch {destroy .browser}
    file delete -force $testDirectory
} -result {{1 1 1} 1 First 1 1 {}}

test browser-1.5 {failed back navigation retains its history entry} -setup {
    set testDirectory [file normalize [file join \
        [::tcltest::temporaryDirectory] browser-back-[pid]]]
    file mkdir $testDirectory
    set firstPath [file join $testDirectory first.html]
    set secondPath [file join $testDirectory second.html]
    foreach {path title} [list $firstPath First $secondPath Second] {
        set channel [open $path w]
        puts $channel "<title>$title</title><p>$title page</p>"
        close $channel
    }
    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand 1
    $browser goto $firstPath
    $browser goto $secondPath
    file delete $firstPath
} -body {
    set firstFailure [catch {$browser back}]
    set secondFailure [catch {$browser back}]
    list $firstFailure $secondFailure [$browser title] \
        [expr {[::tklitehtml::local_path [$browser location]] eq $secondPath}]
} -cleanup {
    catch {destroy .browser}
    file delete -force $testDirectory
} -result {1 1 Second 1}

test browser-1.6 {file-size limits are enforced before reading} -setup {
    set testDirectory [file normalize [file join \
        [::tcltest::temporaryDirectory] browser-limit-[pid]]]
    file mkdir $testDirectory
    set firstPath [file join $testDirectory first.html]
    set largePath [file join $testDirectory large.html]
    set channel [open $firstPath w]
    puts $channel {<title>First</title><p>First page</p>}
    close $channel
    set channel [open $largePath w]
    puts -nonewline $channel [string repeat x 256]
    close $channel
    set browser [::tklitehtml::browser .browser \
        -maxdocumentbytes 128]
    pack $browser -fill both -expand 1
    $browser goto $firstPath
} -body {
    set status [catch {$browser goto $largePath} message]
    list $status \
        [string match {*exceeds the 128 byte limit} $message] \
        [$browser title]
} -cleanup {
    catch {destroy .browser}
    file delete -force $testDirectory
} -result {1 1 First}

test browser-1.7 {unreadable files preserve browser state} -constraints filePermissions -setup {
    set testDirectory [file normalize [file join \
        [::tcltest::temporaryDirectory] browser-permissions-[pid]]]
    file mkdir $testDirectory
    set firstPath [file join $testDirectory first.html]
    set unreadablePath [file join $testDirectory unreadable.html]
    foreach path [list $firstPath $unreadablePath] {
        set channel [open $path w]
        puts $channel {<title>Fixture</title><p>Fixture page</p>}
        close $channel
    }
    file attributes $unreadablePath -permissions 0000
    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand 1
    $browser goto $firstPath
} -body {
    set status [catch {$browser goto $unreadablePath}]
    list $status [$browser title]
} -cleanup {
    catch {file attributes $unreadablePath -permissions 0600}
    catch {destroy .browser}
    file delete -force $testDirectory
} -result {1 Fixture}

test browser-find-1.1 {shared find controls count, navigate, and match case} -setup {
    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand true
    $browser load {<p>alpha Alpha alpha</p>}
    update
} -body {
    set entry [$browser findbar alpha]
    update
    set first [.browser.find.count cget -text]
    .browser.find.next invoke
    set next [.browser.find.count cget -text]
    focus -force $entry
    update
    event generate $entry <Shift-F3>
    set previous [.browser.find.count cget -text]
    .browser.find.case invoke
    list $first $next $previous [.browser.find.count cget -text] \
        [.browser.find.next cget -state]
} -cleanup {
    catch {destroy .browser}
} -result {{1 of 3} {2 of 3} {1 of 3} {1 of 2} normal}

test browser-find-1.2 {find entry consumes application quit shortcuts and Escape} -setup {
    set previousQuit [bind . <q>]
    set previousEscape [bind . <Escape>]
    set ::quitKeys 0
    bind . <q> {incr ::quitKeys}
    bind . <Escape> {incr ::quitKeys}
    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand true
    $browser load {<p>query</p>}
    update
} -body {
    focus -force [$browser html]
    update
    event generate [$browser html] <Control-f>
    update
    focus -force .browser.find.entry
    update
    event generate .browser.find.entry <KeyPress-q>
    event generate .browser.find.entry <KeyRelease-q>
    set query [.browser.find.entry get]
    event generate .browser.find.entry <Escape>
    update
    set consumed $::quitKeys
    set hidden [expr {[winfo manager .browser.find] eq ""}]
    event generate [$browser html] <Escape>
    list $query $consumed $hidden $::quitKeys
} -cleanup {
    bind . <q> $previousQuit
    bind . <Escape> $previousEscape
    unset -nocomplain ::quitKeys
    catch {destroy .browser}
} -result {q 0 1 1}

test browser-find-1.3 {automatic and explicit find ownership is preserved} -setup {
    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand true
    $browser load {<p>alpha beta alpha</p>}
    update
} -body {
    focus -force [$browser html]
    update
    set originalFocus [focus]
    $browser autofind true alpha
    update
    set automatic [list [winfo manager .browser.find] \
        [.browser.find.count cget -text] [expr {[focus] eq $originalFocus}]]
    .browser.find.close invoke
    $browser autofind true beta
    update
    set suppressed [winfo manager .browser.find]
    $browser autofind false
    $browser autofind true beta
    update
    set restored [list [winfo manager .browser.find] \
        [.browser.find.entry get]]
    $browser autofind false
    $browser findbar alpha
    $browser autofind true beta
    $browser autofind false
    update
    set explicit [list [winfo manager .browser.find] \
        [.browser.find.entry get]]
    list $automatic $suppressed $restored $explicit
} -cleanup {
    catch {destroy .browser}
} -result {{pack {1 of 2} 1} {} {pack beta} {pack beta}}

test browser-find-1.4 {context menu preserves selection and exposes Copy and Find} -setup {
    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand true
    $browser load {<p>Café and Café</p>}
    update
} -body {
    $browser selection set 0 4
    update
    ::tklitehtml::_browser_prepare_menu $browser
    set selected [$browser selection get]
    set enabled [.browser.menu entrycget Copy -state]
    .browser.menu invoke Copy
    .browser.menu invoke "Find..."
    update
    set found [.browser.find.count cget -text]
    set retained [$browser selection get]
    .browser.menu invoke "Clear Highlights"
    list $selected $enabled [clipboard get -displayof $browser] \
        $found $retained [$browser find current]
} -cleanup {
    catch {grab release .browser.menu}
    catch {destroy .browser}
} -result {Café normal Café {1 of 2} Café {}}

test browser-find-1.5 {context menu clicks select all and find the selection} -constraints x11 -setup {
    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand true
    $browser load {<p>alpha alpha</p>}
    update
} -body {
    foreach label {{Select All} {Find...}} {
        event generate [$browser html] <ButtonPress-3> -x 10 -y 10
        update
        event generate .browser.menu <ButtonRelease-3> -x 2 -y 2
        update
        set y [expr {[.browser.menu yposition [.browser.menu index $label]] + 4}]
        event generate .browser.menu <Motion> -x 8 -y $y
        update
        event generate .browser.menu <ButtonPress-1> -x 8 -y $y
        event generate .browser.menu <ButtonRelease-1> -x 8 -y $y
        update
    }
    list [$browser selection get] [.browser.find.entry get] \
        [.browser.find.count cget -text]
} -cleanup {
    catch {grab release .browser.menu}
    catch {destroy .browser}
} -result {{alpha alpha} {alpha alpha} {1 of 1}}

test browser-find-1.6 {find query follows successful document changes and cleans up} -setup {
    set browser [::tklitehtml::browser .browser -maxdocumentbytes 128]
    pack $browser -fill both -expand true
    $browser load {<p>needle needle</p>}
    update
    $browser findbar needle
    update
} -body {
    set failed [catch {$browser load [string repeat x 129]}]
    update
    set preserved [.browser.find.count cget -text]
    $browser load {<p>different document</p>}
    update
    set missing [.browser.find.count cget -text]
    set disabled [.browser.find.next cget -state]
    $browser load {<p>needle</p>}
    update
    set found [.browser.find.count cget -text]
    destroy $browser
    list $failed $preserved $missing $disabled $found \
        [info exists ::tklitehtml::browserState(.browser)] \
        [info exists ::tklitehtml::browserCase(.browser)] \
        [info exists ::tklitehtml::browserQuery(.browser)]
} -cleanup {
    catch {destroy .browser}
} -result {1 {1 of 2} {No matches} disabled {1 of 1} 0 0 0}

test browser-find-1.7 {find and selection paint distinct highlights after reflow} -setup {
    set originalGeometry [wm geometry .]
    set originalImages [image names]
    wm geometry . 760x400
    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand true
    $browser load {<p>intro words before the Café target and more words after it</p>}
    update
    set photo {}
    foreach name [image names] {
        if {$name ni $originalImages && [image type $name] eq "photo"} {
            set photo $name
        }
    }
} -body {
    $browser findbar Café
    update
    set index [lindex [$browser find current] 0 0]
    lassign [$browser bbox $index] x y width height
    set findPixel [$photo get $x $y]
    $browser selection set $index [expr {$index + 4}]
    update
    set selectionPixel [$photo get $x $y]
    wm geometry . 200x400
    update
    $browser see $index
    update
    lassign [$browser bbox $index] x y width height
    list $findPixel $selectionPixel [$photo get $x $y] \
        [$browser selection get]
} -cleanup {
    catch {destroy .browser}
    wm geometry . $originalGeometry
    update
} -result {{255 191 105} {128 192 255} {128 192 255} Café}

test browser-find-1.8 {pasted and deleted queries update matches without key events} -setup {
    set browser [::tklitehtml::browser .browser]
    pack $browser -fill both -expand true
    $browser load {<p>Café Café</p>}
    update
    set entry [$browser findbar {}]
    update
    clipboard clear -displayof $browser
    clipboard append -displayof $browser Café
} -body {
    event generate $entry <<Paste>>
    update
    set found [.browser.find.count cget -text]
    $entry delete 0 end
    list $found [.browser.find.count cget -text] [$browser find current]
} -cleanup {
    catch {destroy .browser}
} -result {{1 of 2} {} {}}

test errors-1.1 {invalid arity reports an error without destabilizing the widget} -setup {
    set widget [makeWidget]
} -body {
    set results {}
    foreach command {
        {load}
        {search}
        {highlight}
        {highlight clear one two}
        {find}
        {find next extra}
        {fragment}
        {bbox -1}
        {xview moveto}
    } {
        lappend results [catch {$widget {*}$command}]
    }
    set results
} -cleanup {
    deleteWidget
} -result {1 1 1 1 1 1 1 1 1}

test lifecycle-1.1 {clear and command deletion release widget-owned resources} -setup {
    set widget [makeWidget]
    $widget load {<p>temporary document</p>}
    update
} -body {
    $widget clear
    set empty [expr {[$widget text] eq ""}]
    set noFind [expr {[$widget find current] eq ""}]
    rename $widget {}
    update
    list $empty $noFind [winfo exists .html] [llength [info commands .html]]
} -cleanup {
    deleteWidget
} -result {1 1 0 0}

test lifecycle-1.2 {window and command deletion release images and helper commands} -setup {
    set originalImages [lsort [image names]]
} -body {
    set results {}
    foreach deletion {window command} {
        set widget [makeWidget]
        $widget load {<p>Temporary frame</p>}
        update
        if {$deletion eq "window"} {
            destroy $widget
        } else {
            rename $widget {}
        }
        update
        lappend results \
            [expr {[lsort [image names]] eq $originalImages}] \
            [expr {[llength [info commands .html*]] == 0}]
    }
    set results
} -cleanup {
    deleteWidget
} -result {1 1 1 1}

proc finishTests {} {
    set failed $::tcltest::numTests(Failed)
    # cleanupTests otherwise exits directly when Tk is non-interactive,
    # preventing the failure count from becoming the process status.
    ::tcltest::testConstraint interactive true
    ::tcltest::cleanupTests
    exit [expr {$failed != 0}]
}
finishTests
