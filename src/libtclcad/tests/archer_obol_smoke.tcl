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
