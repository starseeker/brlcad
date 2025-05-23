set(QISST_DIRS)

if(BRLCAD_ENABLE_QT)
  find_package_qt(REQUIRED)
endif(BRLCAD_ENABLE_QT)

if(BRLCAD_ENABLE_QT AND BRLCAD_ENABLE_OPENGL)
  find_package(OpenGL)

  set(
    QISST_DIRS
    ${QISST_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/../adrt
    ${CMAKE_CURRENT_SOURCE_DIR}/../adrt/librender
    ${OPENGL_INCLUDE_DIR}
  )
endif(BRLCAD_ENABLE_QT AND BRLCAD_ENABLE_OPENGL)

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

if(BRLCAD_ENABLE_QT AND BRLCAD_ENABLE_OPENGL)
  set(QISST_LIBS librender libgcv librt libbu)

  if(Qt6Widgets_FOUND)
    qt6_wrap_cpp(moc_srcs ${moc_headers})
    set(
      QISST_LIBS
      ${QISST_LIBS}
      Qt6::Core
      Qt6::Widgets
      Qt6::OpenGL
      Qt6::OpenGLWidgets
    )
  else(Qt6Widgets_FOUND)
    qt5_wrap_cpp(moc_srcs ${moc_headers})
    set(QISST_LIBS ${QISST_LIBS} Qt5::Core Qt5::Widgets Qt5::OpenGL)
  endif(Qt6Widgets_FOUND)

  if(OPENGL_opengl_LIBRARY)
    set(QISST_LIBS ${QISST_LIBS} ${OPENGL_opengl_LIBRARY})
  endif(OPENGL_opengl_LIBRARY)
  if(OPENGL_gl_LIBRARY)
    set(QISST_LIBS ${QISST_LIBS} ${OPENGL_gl_LIBRARY})
  endif(OPENGL_gl_LIBRARY)

  brlcad_addexec(qisst "${isst_srcs};${moc_srcs}" "${QISST_LIBS}")
  target_include_directories(qisst BEFORE PRIVATE ${QISST_DIRS})

  distclean(
    moc_gfile.cpp_parameters
    moc_isstapp.cpp_parameters
    moc_isstgl.cpp_parameters
    moc_main_window.cpp_parameters
  )
endif(BRLCAD_ENABLE_QT AND BRLCAD_ENABLE_OPENGL)

set(isst_ignore_srcs CMakeLists.txt ${isst_srcs} ${moc_headers})

cmakefiles(${isst_ignore_srcs})

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
