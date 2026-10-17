function(SetTargetFolder targetname folder)
  if(TARGET ${targetname})
    set_target_properties(${targetname} PROPERTIES FOLDER "${folder}")
  endif(TARGET ${targetname})
endfunction(SetTargetFolder)

# The following are build or utility tools written by BRL-CAD developers
# specifically for use with our build and maintenance processes

# asciiquack - AsciiDoc processor used to build BRL-CAD documentation
add_subdirectory(asciiquack)

# On some platforms distclean has additional cleanup to do
distclean("${CMAKE_BINARY_DIR}/${BIN_DIR}/asciiquack.pdb")

# inactvhdrs - utility that tries to find local headers unused by other
# source code in a file hierarchy
if(BUILD_TESTING)
  add_subdirectory(inactvhdrs)
  settargetfolder(debug2c "Maintenance Utilities")
  distclean(${CMAKE_CURRENT_BINARY_DIR}/inactvhdrs/Makefile)
endif(BUILD_TESTING)

cmakefiles(
  CMakeLists.txt
  README
  asciiquack/CMakeLists.txt
  asciiquack/LICENSE
  asciiquack/README.md
  asciiquack/asciiquack.cpp
  asciiquack/attr_list.h
  asciiquack/block_scanner.h
  asciiquack/block_scanner_hand.c
  asciiquack/docbook5.hpp
  asciiquack/document.hpp
  asciiquack/html5.hpp
  asciiquack/inline_scanner.hpp
  asciiquack/manpage.hpp
  asciiquack/outbuf.hpp
  asciiquack/parser.cpp
  asciiquack/parser.hpp
  asciiquack/pdf.hpp
  asciiquack/reader.hpp
  asciiquack/substitutors.hpp
  asciiquack/vendor/cxxopts.hpp
  # brlinst is built separately; register its sources for distribution.
  brlinst/.clang-format
  brlinst/.gitignore
  brlinst/CMakeLists.txt
  brlinst/COPYING
  brlinst/README.adoc
  brlinst/THIRD_PARTY.txt
  brlinst/cmake/BRLCAD.cmake
  brlinst/cmake/Brlinst.cmake
  brlinst/cmake/Package.cmake.in
  brlinst/examples/brlcad.json.in
  brlinst/include/installer/core.hpp
  brlinst/include/installer/msi.hpp
  brlinst/include/installer/windows.hpp
  brlinst/resources/notices.hpp.in
  brlinst/resources/runtime.manifest
  brlinst/resources/runtime.rc
  brlinst/src/core.cpp
  brlinst/src/input.cpp
  brlinst/src/msi.cpp
  brlinst/src/msi_actions.cpp
  brlinst/src/msi_model.cpp
  brlinst/src/pack.cpp
  brlinst/src/package.cpp
  brlinst/src/resources.cpp
  brlinst/src/runtime.cpp
  brlinst/src/transaction.cpp
  brlinst/src/windows.cpp
  brlinst/src/wizard.cpp
  brlinst/tests/core_tests.cpp
  brlinst/tests/cpack_test.cmake
  brlinst/tests/msi_model_tests.cpp
  brlinst/tests/windows-msi-smoke.ps1
  brlinst/tests/windows-smoke.ps1
  db2adoc/db2adoc.py
  db2adoc/db2adoc.xsl
  inactvhdrs/CMakeLists.txt
  inactvhdrs/inactvhdrs.cpp
)

# Locac Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
