set(BV_TEST_LIBS "libbv;libbg;libbn;libbu;${M_LIBRARY}")

brlcad_addexec(bv_test view.c "${BV_TEST_LIBS}" TEST)
brlcad_add_test(NAME bv_basic COMMAND bv_test)

cmakefiles(
  CMakeLists.txt
  view.c
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
