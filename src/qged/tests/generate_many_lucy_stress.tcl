# Generate a balanced many-instance combination in a writable copy of lucy.g.
#
# Usage:
#   cp --reflink=auto lucy.g many_lucy.g
#   BOBOL_STRESS_INSTANCE_COUNT=5000 \
#     mged -c many_lucy.g \
#       'source generate_many_lucy_stress.tcl'
#
# A balanced tree is intentional: this stresses display occurrence planning,
# shared PoP residency, and aggregate rendering cost without conflating those
# measurements with pathological boolean-tree recursion.

set instance_count 5000
if {[info exists env(BOBOL_STRESS_INSTANCE_COUNT)]} {
    set requested $env(BOBOL_STRESS_INSTANCE_COUNT)
    if {[string is integer -strict $requested] &&
	$requested > 0 && $requested <= 100000} {
	set instance_count $requested
    }
}

proc balanced_union {nodes} {
    while {[llength $nodes] > 1} {
	set next {}
	set count [llength $nodes]
	for {set i 0} {$i < $count} {incr i 2} {
	    if {$i + 1 >= $count} {
		lappend next [lindex $nodes $i]
	    } else {
		lappend next [list u [lindex $nodes $i] \
		    [lindex $nodes [expr {$i + 1}]]]
	    }
	}
	set nodes $next
    }
    return [lindex $nodes 0]
}

set columns [expr {int(ceil(sqrt(double($instance_count))))}]
set spacing 1200000.0
set leaves {}
for {set i 0} {$i < $instance_count} {incr i} {
    set column [expr {$i % $columns}]
    set row [expr {$i / $columns}]
    set tx [expr {($column - ($columns - 1) / 2.0) * $spacing}]
    set ty [expr {($row - ($columns - 1) / 2.0) * $spacing}]
    set matrix [list \
	1 0 0 $tx \
	0 1 0 $ty \
	0 0 1 0 \
	0 0 0 1]
    lappend leaves [list l lucy.bot.r $matrix]
}

set tree [balanced_union $leaves]
catch {kill many_lucy_stress}
put many_lucy_stress comb region no tree $tree
puts "generated many_lucy_stress with $instance_count instances"
quit
