# Create the shared pawn used by the solid-representation comparison.

title {Solid geometry representation comparison pawn}
units mm

# Deliberate overlaps keep the pawn connected while exercising curved-surface
# intersections during Boolean evaluation.
db put base.s ell V {0 0 4.5} A {17 0 0} B {0 17 0} C {0 0 4.5}
db put base_bead.s ell V {0 0 9.5} A {14 0 0} B {0 14 0} C {0 0 5}
db put body.s ell V {0 0 18} A {12 0 0} B {0 12 0} C {0 0 9}
db put stem.s ell V {0 0 29} A {7 0 0} B {0 7 0} C {0 0 13}
db put collar.s ell V {0 0 39} A {9 0 0} B {0 9 0} C {0 0 3.5}
db put neck.s ell V {0 0 43} A {6.5 0 0} B {0 6.5 0} C {0 0 5}
db put head.s ell V {0 0 50} A {9 0 0} B {0 9 0} C {0 0 9}

db put comparison.csg.r comb region yes id 1000 los 100 GIFTmater 1 rgb {235 190 90} shader {plastic} tree {u {u {u {l base.s} {l base_bead.s}} {u {l body.s} {l stem.s}}} {u {u {l collar.s} {l neck.s}} {l head.s}}}

# The revolve uses a deliberately independent outline.  Its single closed
# sketch approximates the same pawn silhouette without copying the CSG tree.
db put comparison.profile sketch V {0 0 0} A {1 0 0} B {0 0 1} VL {{0 0} {17 0} {17 2} {16.3 4} {15 6} {14 8} {15.4 7} {16 9} {15.4 11} {13.5 12.5} {11.2 14} {12 18} {11 21} {9.5 23} {8 26} {7.5 32} {6.5 36} {7.5 37} {9 39} {7.5 41} {6.5 42} {6.5 43} {7.5 44} {8.5 46} {9 50} {8.5 54} {6.7 56.5} {0 59}} SL {{line S 0 E 1} {line S 1 E 2} {line S 2 E 3} {line S 3 E 4} {line S 4 E 5} {line S 5 E 6} {line S 6 E 7} {line S 7 E 8} {line S 8 E 9} {line S 9 E 10} {line S 10 E 11} {line S 11 E 12} {line S 12 E 13} {line S 13 E 14} {line S 14 E 15} {line S 15 E 16} {line S 16 E 17} {line S 17 E 18} {line S 18 E 19} {line S 19 E 20} {line S 20 E 21} {line S 21 E 22} {line S 22 E 23} {line S 23 E 24} {line S 24 E 25} {line S 25 E 26} {line S 26 E 27} {line S 27 E 0}}
db put comparison.revolve.s revolve V {0 0 0} axis {0 0 1} R {1 0 0} ang 6.283185307179586 sk_name comparison.profile
db put comparison.revolve.r comb region yes id 1001 los 100 GIFTmater 1 rgb {235 190 90} shader {plastic} tree {l comparison.revolve.s}

q
