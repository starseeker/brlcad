brlcad_addexec(imgstream_test stream.c "libimgstream;libicv;libbu" TEST)
brlcad_add_test(NAME imgstream_test COMMAND imgstream_test)

brlcad_addexec(imgstream_fb_compat_test fb_compat.c "libimgstream;libicv;libbu" TEST)
brlcad_add_test(NAME imgstream_fb_compat_test COMMAND imgstream_fb_compat_test)

cmakefiles(CMakeLists.txt)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
