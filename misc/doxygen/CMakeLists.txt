if(DOXYGEN_FOUND)
  include("${CMAKE_SOURCE_DIR}/include/BObol/api_tiers.cmake")
  set(lib_inputs "INPUT += \"${CMAKE_SOURCE_DIR}/include/common.h\"")
  set(
    lib_inputs
    "${lib_inputs}
INPUT += \"${CMAKE_BINARY_DIR}/CMakeTmp/features.dox\""
  )
  set(
    lib_inputs
    "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/misc/doxygen/vmath.dox\""
  )
  set(
    lib_inputs
    "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/include/vmath.h\""
  )
  set(
    lib_inputs
    "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/include//bv.h\""
  )
  set(
    lib_inputs
    "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/include/bio.h\""
  )
  set(
    lib_inputs
    "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/include/bnetwork.h\""
  )
  set(
    lib_inputs
    "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/include/bresource.h\""
  )
  set(
    lib_inputs
    "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/include/bsocket.h\""
  )

  set(lib_predefines)

  # Note: the ordering of libraries in this list is
  # what controls the ordering in the Modules output
  # from Doxygen
  set(
    DOX_LIBS
    common
    libbu
    libbn
    libbg
    libnmg
    libbrep
    libfft
    libicv
    libpc
    libwdb
    librt
    liboptical
    libgcv
    libanalyze
    libBObol
    libged
    libtclcad
    libpkg
  )

  foreach(libname ${DOX_LIBS})
    set(doxlibname ${libname})
    if("${libname}" STREQUAL "libBObol")
      # Preserve the installed mixed-case include directory while using the
      # actual export macro spelling.
      set(sublibname BObol)
      set(upper_name BOBOL)
    else()
      string(REPLACE "lib" "" sublibname ${libname})
      string(TOUPPER ${sublibname} upper_name)
    endif()

    set(DOXY_EXTRACT_ALL YES)
    set(DOXY_WARN_AS_ERROR NO)
    set(DOXY_CPP_PREDEFINE "__cplusplus")

    # Build local input list
    set(
      INPUTS
      "
INPUT += \"${CMAKE_SOURCE_DIR}/misc/doxygen/${libname}.dox\""
    )
    if(EXISTS "${CMAKE_SOURCE_DIR}/include/${sublibname}/")
      set(
        INPUTS
        "${INPUTS}
INPUT += \"${CMAKE_SOURCE_DIR}/include/${sublibname}/\""
      )
    endif(EXISTS "${CMAKE_SOURCE_DIR}/include/${sublibname}/")
    if(EXISTS "${CMAKE_SOURCE_DIR}/include/${sublibname}.h")
      set(
        INPUTS
        "${INPUTS}
INPUT += \"${CMAKE_SOURCE_DIR}/include/${sublibname}.h\"
    "
      )
    endif(EXISTS "${CMAKE_SOURCE_DIR}/include/${sublibname}.h")

    # libbn header doesn't fit the template
    if("${libname}" STREQUAL "libbn")
      set(
        INPUTS
        "${INPUTS}
INPUT += \"${CMAKE_SOURCE_DIR}/include/bn/vector_fpu.h\"
    "
      )
      set(
        INPUTS
        "${INPUTS}
INPUT += \"${CMAKE_SOURCE_DIR}/include/bn/dvec.h\"
    "
      )
    endif("${libname}" STREQUAL "libbn")

    # librt header doesn't fit the template
    if("${libname}" STREQUAL "librt")
      set(
        INPUTS
        "${INPUTS}
INPUT += \"${CMAKE_SOURCE_DIR}/include/raytrace.h\"
    "
      )
      set(
        INPUTS
        "${INPUTS}
INPUT += \"${CMAKE_SOURCE_DIR}/include/nmg.h\"
    "
      )
    endif("${libname}" STREQUAL "librt")

    if(DEFINED BRLCAD_DOXYGEN_ADD_SRCS)
      set(
        INPUTS
        "${INPUTS}
INPUT += \"${CMAKE_SOURCE_DIR}/src/${libname}/\""
      )
    endif(DEFINED BRLCAD_DOXYGEN_ADD_SRCS)

    set(
      PREDEFINES
      "
PREDEFINED += \"${upper_name}_EXTERN(type_and_name,args)=extern type_and_name args\"
PREDEFINED += \"${upper_name}_EXPORT=\"
  "
    )

    # Build top level input list
    set(
      lib_inputs
      "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/misc/doxygen/${libname}.dox\""
    )
    if(EXISTS "${CMAKE_SOURCE_DIR}/include/${sublibname}/")
      set(
        lib_inputs
        "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/include/${sublibname}/\""
      )
    endif(EXISTS "${CMAKE_SOURCE_DIR}/include/${sublibname}/")
    if(EXISTS "${CMAKE_SOURCE_DIR}/include/${sublibname}.h")
      set(
        lib_inputs
        "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/include/${sublibname}.h\"
    "
      )
    endif(EXISTS "${CMAKE_SOURCE_DIR}/include/${sublibname}.h")

    # librt header doesn't fit the template
    if("${libname}" STREQUAL "librt")
      set(
        lib_inputs
        "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/include/raytrace.h\"
    "
      )
      set(
        lib_inputs
        "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/include/nmg.h\"
    "
      )
    endif("${libname}" STREQUAL "librt")

    if(DEFINED BRLCAD_DOXYGEN_ADD_SRCS)
      set(
        lib_inputs
        "${lib_inputs}
INPUT += \"${CMAKE_SOURCE_DIR}/src/${libname}/\""
      )
    endif(DEFINED BRLCAD_DOXYGEN_ADD_SRCS)

    set(
      lib_predefines
      "${lib_predefines}
PREDEFINED += \"${upper_name}_EXTERN(type_and_name,args)=extern type_and_name args\"
PREDEFINED += \"${upper_name}_EXPORT=\"
  "
    )

    configure_file(Doxyfile.in "${CMAKE_BINARY_DIR}/CMakeTmp/Doxyfile_${libname}" @ONLY)
    add_custom_command(
      OUTPUT "${CMAKE_BINARY_DIR}/doc/doxygen_${libname}/html/index.html"
      COMMAND "${DOXYGEN_EXECUTABLE}" "${CMAKE_BINARY_DIR}/CMakeTmp/Doxyfile_${libname}"
      WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/doc"
    )
    add_custom_target(dox-${libname} DEPENDS "${CMAKE_BINARY_DIR}/doc/doxygen_${libname}/html/index.html")

    if("${libname}" STREQUAL "libBObol")
      # The full module remains navigable, including experimental C++ source
      # APIs.  A second configuration is the release gate: it parses only the
      # SOVERSION-backed C ABI and treats every documentation warning as an
      # error.  Keeping the inputs machine-derived prevents tier drift.
      set(_bobol_stable_inputs
        "\nINPUT += \"${CMAKE_SOURCE_DIR}/misc/doxygen/libBObol.dox\"")
      set(_bobol_stable_header_paths)
      foreach(_bobol_header IN LISTS LIBBOBOL_STABLE_C_ABI_HEADERS)
        string(APPEND _bobol_stable_inputs
          "\nINPUT += \"${CMAKE_SOURCE_DIR}/include/BObol/${_bobol_header}\"")
        list(APPEND _bobol_stable_header_paths
          "${CMAKE_SOURCE_DIR}/include/BObol/${_bobol_header}")
      endforeach()
      set(INPUTS "${_bobol_stable_inputs}")
      set(doxlibname libBObol_stable)
      set(DOXY_EXTRACT_ALL NO)
      set(DOXY_WARN_AS_ERROR YES)
      set(DOXY_CPP_PREDEFINE "BOBOL_DOXY_C_ABI")
      configure_file(Doxyfile.in
        "${CMAKE_BINARY_DIR}/CMakeTmp/Doxyfile_libBObol_stable" @ONLY)
      add_custom_command(
        OUTPUT "${CMAKE_BINARY_DIR}/doc/doxygen_libBObol_stable/html/index.html"
        COMMAND "${DOXYGEN_EXECUTABLE}"
          "${CMAKE_BINARY_DIR}/CMakeTmp/Doxyfile_libBObol_stable"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/doc"
        DEPENDS
          "${CMAKE_BINARY_DIR}/CMakeTmp/Doxyfile_libBObol_stable"
          "${CMAKE_SOURCE_DIR}/misc/doxygen/libBObol.dox"
          ${_bobol_stable_header_paths}
        VERBATIM)
      add_custom_target(dox-libBObol-stable
        DEPENDS "${CMAKE_BINARY_DIR}/doc/doxygen_libBObol_stable/html/index.html")
      add_dependencies(dox-libBObol dox-libBObol-stable)
      set(doxlibname ${libname})
    endif()
    cmakefiles(${libname}.dox)
  endforeach(libname ${DOX_LIBS})

  set(doxlibname output)
  set(INPUTS ${lib_inputs})
  set(PREDEFINES ${lib_predefines})
  set(DOXY_EXTRACT_ALL YES)
  set(DOXY_WARN_AS_ERROR NO)
  set(DOXY_CPP_PREDEFINE "__cplusplus")
  configure_file(Doxyfile.in "${CMAKE_BINARY_DIR}/CMakeTmp/Doxyfile" @ONLY)
  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/doc/doxygen_output/html/index.html"
    COMMAND "${DOXYGEN_EXECUTABLE}" "${CMAKE_BINARY_DIR}/CMakeTmp/Doxyfile"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/doc"
  )
  add_custom_target(dox DEPENDS "${CMAKE_BINARY_DIR}/doc/doxygen_output/html/index.html")
  add_custom_target(doxygen DEPENDS "${CMAKE_BINARY_DIR}/doc/doxygen_output/html/index.html")
endif(DOXYGEN_FOUND)

set(
  doxygen_ignore_files
  CMakeLists.txt
  Doxyfile.in
  common.dox
  images/README.txt
  images/rt_pattern_circ_spiral.png
  images/rt_pattern_rect_orthogrid.png
  intro.dox
  libanalyze.dox
  libBObol.dox
  libbg.dox
  libbn.dox
  libbrep.dox
  libbu.dox
  libfft.dox
  libgcv.dox
  libged.dox
  libicv.dox
  libnmg.dox
  liboptical.dox
  libpc.dox
  libpkg.dox
  librt.dox
  libtclcad.dox
  libwdb.dox
  vmath.dox
)
cmakefiles(${doxygen_ignore_files})

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
