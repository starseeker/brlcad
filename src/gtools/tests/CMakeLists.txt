brlcad_addexec(bigdb bigdb.c "libged;libwdb;librt;libbu" TEST_USESDATA)

brlcad_add_test(NAME slow-bigdb_1gb COMMAND bigdb 1)
if(TEST slow-bigdb_1gb)
  set_property(TEST slow-bigdb_1gb PROPERTY SKIP_RETURN_CODE 123)
endif(TEST slow-bigdb_1gb)

brlcad_add_test(NAME slow-bigdb_5gb COMMAND bigdb 5)
if(TEST slow-bigdb_5gb)
  set_property(TEST slow-bigdb_5gb PROPERTY SKIP_RETURN_CODE 123)
  # On Windows this test can consume most of the system commit limit.  If it
  # overlaps GUI or framebuffer process-launch tests, CreateProcess may fail
  # with STATUS_COMMITMENT_LIMIT (0xc000012d) even when both tests are valid
  # independently.
  if(WIN32)
    set_property(TEST slow-bigdb_5gb PROPERTY RUN_SERIAL TRUE)
  endif()
endif(TEST slow-bigdb_5gb)

cmakefiles(CMakeLists.txt)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
