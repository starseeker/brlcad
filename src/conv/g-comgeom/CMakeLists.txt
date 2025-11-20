# BRL-CAD g-comgeom CMakeLists.txt

set(G_COMGEOM_SRCS
  g-comgeom.c
)

BRLCAD_ADDEXEC(g-comgeom "${G_COMGEOM_SRCS}" "librt;libwdb;libbu" NO_STRICT)

add_subdirectory(tests)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
