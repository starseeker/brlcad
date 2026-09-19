# BRL-CAD
# Copyright (c) 2026 United States Government as represented by
# the U.S. Army Research Laboratory.

# MGED's interactive command reader expands bracket substitutions before
# handing a line to Tcl.  Keep the small amount of Tcl introspection needed by
# the native Tk Obol smoke in a sourced script so list-valued bindtags remain a
# single Tcl argument.
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
