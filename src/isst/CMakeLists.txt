set(QISST_DIRS)

if(BRLCAD_ENABLE_QT)
  find_package_qt(REQUIRED)
endif(BRLCAD_ENABLE_QT)

if(BRLCAD_ENABLE_QT)
  set(
    QISST_DIRS
    ${QISST_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../adrt
    ${CMAKE_CURRENT_SOURCE_DIR}/../adrt/librender
  )
endif(BRLCAD_ENABLE_QT)

# NOTE: We can't use the CMake global values for automoc et. al.
# because most of BRL-CAD's targets are not Qt targets.  We don't
# want to run the automatic moc logic for all of them.

set(
  isst_srcs
  gfile.cpp
  main.cpp
  main_window.cpp
  isstapp.cpp
  isstgl.cpp
)

set(moc_headers gfile.h main_window.h isstapp.h isstgl.h)

if(BRLCAD_ENABLE_QT)
  set(QISST_LIBS librender libgcv librt libbu)

  if(Qt6Widgets_FOUND)
    qt6_wrap_cpp(moc_srcs ${moc_headers})
    set(
      QISST_LIBS
      ${QISST_LIBS}
      Qt6::Core
      Qt6::Widgets
    )
  else(Qt6Widgets_FOUND)
    qt5_wrap_cpp(moc_srcs ${moc_headers})
    set(QISST_LIBS ${QISST_LIBS} Qt5::Core Qt5::Widgets)
  endif(Qt6Widgets_FOUND)

  brlcad_addexec(qisst "${isst_srcs};${moc_srcs}" "${QISST_LIBS}")
  target_include_directories(qisst BEFORE PRIVATE ${QISST_DIRS})
  if(UNIX AND SH_EXEC)
    find_program(QISST_XVFB_RUN_EXECUTABLE xvfb-run)
    if(QISST_XVFB_RUN_EXECUTABLE)
      brlcad_add_test(NAME qisst_qt_image_smoke
	COMMAND ${SH_EXEC} ${BRLCAD_SOURCE_DIR}/misc/CMake/scripts/xvfb_test.sh
	${QISST_XVFB_RUN_EXECUTABLE} $<TARGET_FILE:qisst>
	${CMAKE_BINARY_DIR}/share/db/axis.g all/)
      set_property(TEST qisst_qt_image_smoke APPEND PROPERTY ENVIRONMENT
	"BRLCAD_ROOT=${CMAKE_BINARY_DIR};BRLCAD_QISST_SMOKE=1")
      set_property(TEST qisst_qt_image_smoke PROPERTY SKIP_RETURN_CODE 125)
    endif(QISST_XVFB_RUN_EXECUTABLE)
  endif(UNIX AND SH_EXEC)

  distclean(
    moc_gfile.cpp_parameters
    moc_isstapp.cpp_parameters
    moc_isstgl.cpp_parameters
    moc_main_window.cpp_parameters
  )
endif(BRLCAD_ENABLE_QT)

set(isst_ignore_srcs CMakeLists.txt ${isst_srcs} ${moc_headers})

cmakefiles(${isst_ignore_srcs})

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
