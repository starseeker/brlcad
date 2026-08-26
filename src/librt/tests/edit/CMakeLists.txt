brlcad_addexec(rt_edit_test_arbn arbn.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_datum datum.cpp "librt;libwdb" TEST)
brlcad_addexec(rt_edit_test_annot annot.cpp "librt;libwdb" TEST)
brlcad_addexec(rt_edit_test_revolve revolve.cpp "librt;libwdb" TEST)
brlcad_addexec(rt_edit_test_tor tor.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_ell ell.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_tgc tgc.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_epa epa.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_ehy ehy.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_eto eto.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_hyp hyp.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_rpc rpc.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_rhc rhc.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_part part.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_rec rec.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_superell superell.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_cline cline.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_extrude extrude.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_metaball metaball.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_hrt hrt.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_ars ars.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_pipe pipe.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_dsp dsp.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_vol vol.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_ebm ebm.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_sph sph.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_half half.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_arb8 arb8.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_bot bot.cpp "librt" TEST)

brlcad_addexec(rt_edit_test_sketch sketch.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_bspline bspline.cpp "librt" TEST)
brlcad_addexec(rt_edit_test_brep brep.cpp "librt;libbrep" TEST)
brlcad_addexec(rt_edit_test_nmg nmg.cpp "librt;libnmg;libwdb" TEST)
brlcad_addexec(rt_edit_test_comb comb.cpp "librt;libwdb" TEST)

brlcad_addexec(rt_edit_test_edit_desc edit_desc.cpp "librt" TEST)
brlcad_add_test(NAME rt_edit_descriptor_contract COMMAND rt_edit_test_edit_desc)
set_property(TEST rt_edit_descriptor_contract APPEND PROPERTY LABELS
  edit_contract)
brlcad_addexec(rt_edit_test_str_get_values edit_str_get_params.cpp "librt" TEST)

# TEST builds these executables, but does not register them with CTest.  Keep
# the descriptor contract and the actual primitive handlers in one CI gate:
# descriptor metadata alone cannot catch a runtime edit regression.
set(rt_edit_runtime_targets
  rt_edit_test_arbn
  rt_edit_test_datum
  rt_edit_test_annot
  rt_edit_test_revolve
  rt_edit_test_tor
  rt_edit_test_ell
  rt_edit_test_tgc
  rt_edit_test_epa
  rt_edit_test_ehy
  rt_edit_test_eto
  rt_edit_test_hyp
  rt_edit_test_rpc
  rt_edit_test_rhc
  rt_edit_test_part
  rt_edit_test_rec
  rt_edit_test_superell
  rt_edit_test_cline
  rt_edit_test_extrude
  rt_edit_test_metaball
  rt_edit_test_hrt
  rt_edit_test_ars
  rt_edit_test_pipe
  rt_edit_test_dsp
  rt_edit_test_vol
  rt_edit_test_ebm
  rt_edit_test_sph
  rt_edit_test_half
  rt_edit_test_arb8
  rt_edit_test_bot
  rt_edit_test_sketch
  rt_edit_test_bspline
  rt_edit_test_brep
  rt_edit_test_nmg
  rt_edit_test_comb
  rt_edit_test_str_get_values
)

foreach(rt_edit_runtime_target ${rt_edit_runtime_targets})
  string(REGEX REPLACE "^rt_edit_test_" "rt_edit_runtime_"
    rt_edit_runtime_name ${rt_edit_runtime_target})
  brlcad_add_test(NAME ${rt_edit_runtime_name}
    COMMAND ${rt_edit_runtime_target})
  set_property(TEST ${rt_edit_runtime_name} APPEND PROPERTY LABELS
    edit_contract;edit_runtime)
endforeach()

set(
  distcheck_files
  CMakeLists.txt
  arbn.cpp
  annot.cpp
  datum.cpp
  revolve.cpp
  tor.cpp
  ell.cpp
  tgc.cpp
  epa.cpp
  ehy.cpp
  eto.cpp
  hyp.cpp
  rpc.cpp
  rhc.cpp
  part.cpp
  rec.cpp
  superell.cpp
  cline.cpp
  extrude.cpp
  sketch.cpp
  bspline.cpp
  brep.cpp
  metaball.cpp
  hrt.cpp
  ars.cpp
  pipe.cpp
  dsp.cpp
  vol.cpp
  ebm.cpp
  sph.cpp
  half.cpp
  arb8.cpp
  bot.cpp
  nmg.cpp
  comb.cpp
  edit_desc.cpp
  edit_test_view.h
  edit_str_get_params.cpp
)

cmakefiles(${distcheck_files})

# Local Variables:
# tab-width: 8
# mode: cmake
# indent-tabs-mode: t
# End:
# ex: shiftwidth=2 tabstop=8
