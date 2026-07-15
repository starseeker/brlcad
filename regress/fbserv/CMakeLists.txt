# Stress test for the fbserv token-isolation + parallel-connection mechanism.
#
# Simulates multiple simultaneous rtwizard-style invocations: each worker
# launches its own fbserv on a unique port with a unique token, connects using
# libimgstream's remote framebuffer API, and verifies authenticated operations.
#
# Two ctest registrations are provided:
#
#   regress-fbserv-stress-quick   -- 8 parallel workers, run in normal CI
#   regress-fbserv-stress-heavy   -- 32 parallel workers, opt-in stress run

if(TARGET fbserv AND TARGET libimgstream AND TARGET libpkg AND TARGET libbu)

  brlcad_addexec(fbserv_stress fbserv_stress.cpp "libimgstream;libpkg;libbu;Threads::Threads;${WINSOCK_LIB}" TEST)
  if(TARGET fbserv_stress)
    set_target_properties(fbserv_stress PROPERTIES FOLDER "BRL-CAD Regression Tests/fbserv")

    # Quick variant: 8 workers.  Port allocation is PID-based (so concurrent
    # ctest instances automatically use distinct port ranges).
    brlcad_add_test(
      NAME regress-fbserv-stress-quick
      COMMAND fbserv_stress --workers 8
    )
    brlcad_regression_test(regress-fbserv-stress-quick "fbserv;fbserv_stress" STAND_ALONE TEST_DEFINED)
    distclean(${CMAKE_CURRENT_BINARY_DIR}/regress-fbserv-stress-quick.log)

    brlcad_add_test(
      NAME regress-fbserv-stress-ipc
      COMMAND fbserv_stress --workers 8 --ipc
    )
    brlcad_regression_test(regress-fbserv-stress-ipc "fbserv;fbserv_stress" STAND_ALONE TEST_DEFINED)
    distclean(${CMAKE_CURRENT_BINARY_DIR}/regress-fbserv-stress-ipc.log)

    # Heavy variant: 32 workers.  Same PID-based port isolation applies.
    brlcad_add_test(
      NAME regress-fbserv-stress-heavy
      COMMAND fbserv_stress --workers 32
    )
    brlcad_regression_test(regress-fbserv-stress-heavy "fbserv;fbserv_stress" STAND_ALONE TEST_DEFINED)
    distclean(${CMAKE_CURRENT_BINARY_DIR}/regress-fbserv-stress-heavy.log)

  endif(TARGET fbserv_stress)

endif(TARGET fbserv AND TARGET libimgstream AND TARGET libpkg AND TARGET libbu)

cmakefiles(
  CMakeLists.txt
  fbserv_stress.cpp
)

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
