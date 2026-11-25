# BRL-CAD
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.

# MGED's interactive command reader expands bracket substitutions before
# handing a line to Tcl.  Keep the small amount of Tcl introspection needed by
# the native Tk Obol smoke in a sourced script so list-valued bindtags remain a
# single Tcl argument.
if {![info exists ::brlcad::mged::display::semantic_input(.mged_direct)] ||
    [info exists ::mged_obol_semantic_input(.mged_direct)] ||
    [llength [info commands ::mged_bind_endpoint]] ||
    [bind .mged_direct <Enter>] ne
	[list ::brlcad::mged::display::activate_endpoint .mged_direct] ||
    [bind .mged_direct <ButtonRelease>] ne
	"[list ::brlcad::mged::display::release_endpoint .mged_direct]; break"} {
    error "MGED Tk Obol bindings did not use the display namespace"
}
puts OBOL_DISPLAY_NAMESPACE_OK

proc check_obol_display_released {} {
    if {[info exists ::brlcad::mged::display::semantic_input(.mged_direct)] ||
	[winfo exists .mged_direct]} {
	error "MGED retained the display widget or namespace state after release"
    }
    puts OBOL_DISPLAY_RELEASE_OK
}

foreach obol_input_tag [bindtags .mged_direct] {
    break
}
regexp {^([^ ]+)} [bind $obol_input_tag <KeyPress>] \
    obol_input_match obol_input_command

puts OBOL_NATIVE_KEY_BINDINGS
puts [string cat \
    [bind .mged_direct <F3>] "|" \
    [bind .mged_direct <F6>] "|" \
    [bind .mged_direct <F7>] "|" \
    [bind .mged_direct <F8>] "|" \
    [bind .mged_direct e]]
