if {$argc != 1} {
    puts stderr "Usage: new_view_cleanup.tcl database.g"
    exit 1
}

go_open g db [lindex $argv 0]
g new_view v1 nu
rename g {}

puts "TclCAD new_view cleanup passed"
