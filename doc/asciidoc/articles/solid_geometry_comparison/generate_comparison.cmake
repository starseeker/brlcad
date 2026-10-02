# Regenerate the solid geometry representation comparison snapshots.
#
# The article is hand-written.  This script builds its common pawn model,
# derives each representation using BRL-CAD tools, and replaces the committed
# images.  All databases and sidecar data remain in WORKDIR.

cmake_minimum_required(VERSION 3.20)

foreach(v MGED RT GVOXEL ASC2PIX FBCLEAR PLOT3FB FBPNG
          MODEL_SCRIPT IMGDIR WORKDIR)
  if(NOT DEFINED ${v})
    message(FATAL_ERROR "generate_comparison: -D${v}=... is required")
  endif()
endforeach()

if(NOT DEFINED SIZE)
  set(SIZE 512)
endif()

set(RENDER_TIMEOUT 1800)
set(PAWN_DB "${WORKDIR}/solid_representation_pawn.g")
set(VIEW_SCRIPT "${WORKDIR}/pawn.view")
set(MGED_SCRIPT "${WORKDIR}/assemble_variants.tcl")
set(VOXELIZE_SCRIPT "${WORKDIR}/voxelize.tcl")
set(VOXEL_LIST "${WORKDIR}/voxels1.txt")
set(VOLUME_ASCII "${WORKDIR}/comparison_volume.asc")
set(VOLUME_FILE "${WORKDIR}/comparison.vol")

file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")
file(MAKE_DIRECTORY "${IMGDIR}")

function(run_process label)
  execute_process(
    COMMAND ${ARGN}
    WORKING_DIRECTORY "${WORKDIR}"
    TIMEOUT ${RENDER_TIMEOUT}
    RESULT_VARIABLE process_result
    OUTPUT_VARIABLE process_stdout
    ERROR_VARIABLE process_stderr
  )
  if(NOT process_result EQUAL 0)
    message(
      FATAL_ERROR
      "${label} failed (exit ${process_result})\n${process_stdout}${process_stderr}"
    )
  endif()
  message(STATUS "${label}")
endfunction()

function(require_brep_property object property expected_output)
  execute_process(
    COMMAND "${MGED}" --no-rc -c "${PAWN_DB}" brep "${object}" "${property}"
    WORKING_DIRECTORY "${WORKDIR}"
    TIMEOUT ${RENDER_TIMEOUT}
    RESULT_VARIABLE property_result
    OUTPUT_VARIABLE property_stdout
    ERROR_VARIABLE property_stderr
  )
  set(property_output "${property_stdout}${property_stderr}")
  if(NOT property_result EQUAL 0 OR NOT property_output MATCHES "${expected_output}")
    message(
      FATAL_ERROR
      "B-rep ${object} is not ${property}\n${property_output}"
    )
  endif()
  message(STATUS "Validated B-rep ${object}: ${property}")
endfunction()

execute_process(
  COMMAND "${MGED}" --no-rc -c "${PAWN_DB}"
  WORKING_DIRECTORY "${WORKDIR}"
  TIMEOUT ${RENDER_TIMEOUT}
  INPUT_FILE "${MODEL_SCRIPT}"
  RESULT_VARIABLE create_result
  OUTPUT_VARIABLE create_stdout
  ERROR_VARIABLE create_stderr
)
if(NOT create_result EQUAL 0 OR NOT EXISTS "${PAWN_DB}")
  message(
    FATAL_ERROR
    "Creating the comparison pawn failed (exit ${create_result})\n${create_stdout}${create_stderr}"
  )
endif()

# Preserve one Boolean hierarchy as B-rep leaves, and evaluate another into a
# single B-rep.  The default suffix is intentional: supplying an explicit
# suffix changes parsing in older MGED versions.
set(PAWN_LEAVES base base_bead body stem collar neck head)
run_process(
  "Created Boolean tree of B-rep leaves"
  "${MGED}" --no-rc -c "${PAWN_DB}" brep comparison.csg.r brep --no-evaluation
)
foreach(leaf IN LISTS PAWN_LEAVES)
  require_brep_property("${leaf}.s.brep" solid "brep is solid")
endforeach()
run_process(
  "Created evaluated B-rep"
  "${MGED}" --no-rc -c "${PAWN_DB}" brep comparison.csg.r brep comparison.brep_evaluated.s
)
require_brep_property(comparison.brep_evaluated.s valid "brep is valid")
require_brep_property(comparison.brep_evaluated.s solid "brep is solid")

# Generate both the topology-rich NMG form and the triangle-only BoT form.
run_process(
  "Created evaluated NMG"
  "${MGED}" --no-rc -c "${PAWN_DB}" facetize -n comparison.csg.r comparison.nmg.s
)
run_process(
  "Created evaluated BoT"
  "${MGED}" --no-rc -c "${PAWN_DB}" facetize --methods NMG comparison.csg.r comparison.bot.s
)

# Tessellating leaves independently preserves the editable Boolean hierarchy.
foreach(leaf IN LISTS PAWN_LEAVES)
  run_process(
    "Tessellated ${leaf}.s"
    "${MGED}" --no-rc -c "${PAWN_DB}" facetize --methods NMG "${leaf}.s" "${leaf}.bot"
  )
endforeach()

file(WRITE "${VOXELIZE_SCRIPT}" "voxelize -s {3 3 3} -d 4 -t 0.45 comparison.voxel_cells.r comparison.csg.r\nq\n")
execute_process(
  COMMAND "${MGED}" --no-rc -c "${PAWN_DB}"
  WORKING_DIRECTORY "${WORKDIR}"
  TIMEOUT ${RENDER_TIMEOUT}
  INPUT_FILE "${VOXELIZE_SCRIPT}"
  RESULT_VARIABLE voxelize_result
  OUTPUT_VARIABLE voxelize_stdout
  ERROR_VARIABLE voxelize_stderr
)
if(NOT voxelize_result EQUAL 0)
  message(
    FATAL_ERROR
    "Creating explicit voxel cells failed (exit ${voxelize_result})\n${voxelize_stdout}${voxelize_stderr}"
  )
endif()
message(STATUS "Created explicit voxel cells")

# PNTS stores point size at import time.  The generator first writes all
# sampled fields (xyzijk), then reads them back with a visible point diameter.
run_process(
  "Sampled pawn surface"
  "${MGED}" --no-rc -c "${PAWN_DB}" pnts gen -t 2 --surface --grid
  comparison.csg.r comparison.samples
)
run_process(
  "Exported pawn point samples"
  "${MGED}" --no-rc -c "${PAWN_DB}" pnts write comparison.samples comparison.xyz
)
run_process(
  "Created visible point cloud"
  "${MGED}" --no-rc -c "${PAWN_DB}" pnts read -f xyzijk --size 0.45
  comparison.xyz comparison.point_cloud.s
)

# g-voxel reports occupied cells.  Expand that sparse report into the regular,
# x-fastest byte array required by VOL.  The pawn bounds are [-17,17] in X/Y
# and [0,60] in Z; 2 mm cells therefore form a 17 x 17 x 30 grid.
file(REMOVE "${VOXEL_LIST}")
run_process(
  "Sampled VOL grid"
  "${GVOXEL}" -s "2 2 2" -d 4 -t 0.5 "${PAWN_DB}" comparison.csg.r
)
if(NOT EXISTS "${VOXEL_LIST}")
  message(FATAL_ERROR "g-voxel did not create ${VOXEL_LIST}")
endif()

file(STRINGS "${VOXEL_LIST}" voxel_lines)
list(LENGTH voxel_lines occupied_voxels)
if(occupied_voxels EQUAL 0)
  message(FATAL_ERROR "g-voxel reported no occupied cells")
endif()

foreach(line IN LISTS voxel_lines)
  string(
    REGEX MATCH
    "^\\((-?[0-9]+)\\.0+, (-?[0-9]+)\\.0+, (-?[0-9]+)\\.0+\\)[ \t]+[^ \t]+[ \t]+[0-9.]+$"
    parsed "${line}"
  )
  if(NOT parsed)
    message(FATAL_ERROR "Cannot parse g-voxel output: ${line}")
  endif()
  set(x_coord "${CMAKE_MATCH_1}")
  set(y_coord "${CMAKE_MATCH_2}")
  set(z_coord "${CMAKE_MATCH_3}")
  math(EXPR x_index "(${x_coord} + 16) / 2")
  math(EXPR y_index "(${y_coord} + 16) / 2")
  math(EXPR z_index "(${z_coord} - 1) / 2")
  if(x_index LESS 0 OR x_index GREATER 16 OR
     y_index LESS 0 OR y_index GREATER 16 OR
     z_index LESS 0 OR z_index GREATER 29)
    message(FATAL_ERROR "g-voxel coordinate outside the expected pawn grid: ${line}")
  endif()
  set("voxel_${x_index}_${y_index}_${z_index}" 255)
endforeach()

set(volume_text "")
foreach(z_index RANGE 0 29)
  foreach(y_index RANGE 0 16)
    foreach(x_index RANGE 0 16)
      if(DEFINED voxel_${x_index}_${y_index}_${z_index})
        set(sample ff)
      else()
        set(sample 00)
      endif()
      string(APPEND volume_text "${sample}\n")
    endforeach()
  endforeach()
endforeach()
file(WRITE "${VOLUME_ASCII}" "${volume_text}")

execute_process(
  COMMAND "${ASC2PIX}"
  WORKING_DIRECTORY "${WORKDIR}"
  TIMEOUT ${RENDER_TIMEOUT}
  INPUT_FILE "${VOLUME_ASCII}"
  OUTPUT_FILE "${VOLUME_FILE}"
  RESULT_VARIABLE asc2pix_result
  ERROR_VARIABLE asc2pix_error
)
if(NOT asc2pix_result EQUAL 0)
  message(FATAL_ERROR "asc2pix failed (exit ${asc2pix_result})\n${asc2pix_error}")
endif()
file(SIZE "${VOLUME_FILE}" volume_size)
if(NOT volume_size EQUAL 8670)
  message(FATAL_ERROR "VOL sidecar is ${volume_size} bytes; expected 8670")
endif()
message(STATUS "Encoded ${occupied_voxels} occupied VOL cells")

# Wrap single primitives in regions so every shaded comparison uses the same
# material.  The polygonal tree repeats the original CSG union hierarchy with
# BoT leaves.  VOL references a sidecar basename resolved beside the database.
file(WRITE "${MGED_SCRIPT}" [=[
db put comparison.brep_evaluated.r comb region yes id 1010 los 100 GIFTmater 1 rgb {235 190 90} shader {plastic} tree {l comparison.brep_evaluated.s}
db put comparison.nmg.r comb region yes id 1020 los 100 GIFTmater 1 rgb {235 190 90} shader {plastic} tree {l comparison.nmg.s}
db put comparison.bot.r comb region yes id 1030 los 100 GIFTmater 1 rgb {235 190 90} shader {plastic} tree {l comparison.bot.s}
db put comparison.polygonal_boolean.r comb region yes id 1040 los 100 GIFTmater 1 rgb {235 190 90} shader {plastic} tree {u {u {u {l base.bot} {l base_bead.bot}} {u {l body.bot} {l stem.bot}}} {u {u {l collar.bot} {l neck.bot}} {l head.bot}}}
attr set comparison.voxel_cells.r region R region_id 1050 material_id 1 los 100 shader plastic color 235/190/90
db put comparison.point_cloud.r comb region yes id 1060 los 100 GIFTmater 1 rgb {235 190 90} shader {plastic} tree {l comparison.point_cloud.s}
db put comparison.volume.s vol file {comparison.vol} name {comparison.vol} src f w 17 n 17 d 30 lo 1 hi 255 size {2 2 2} mat {1 0 0 -17 0 1 0 -17 0 0 1 0 0 0 0 1}
db put comparison.volume.r comb region yes id 1070 los 100 GIFTmater 1 rgb {235 190 90} shader {plastic} tree {l comparison.volume.s}
q
]=])

execute_process(
  COMMAND "${MGED}" --no-rc -c "${PAWN_DB}"
  WORKING_DIRECTORY "${WORKDIR}"
  TIMEOUT ${RENDER_TIMEOUT}
  INPUT_FILE "${MGED_SCRIPT}"
  RESULT_VARIABLE assemble_result
  OUTPUT_VARIABLE assemble_stdout
  ERROR_VARIABLE assemble_stderr
)
if(NOT assemble_result EQUAL 0)
  message(
    FATAL_ERROR
    "Assembling comparison regions failed (exit ${assemble_result})\n${assemble_stdout}${assemble_stderr}"
  )
endif()

file(WRITE "${VIEW_SCRIPT}" "viewsize 8.000000000000000e+01;\n")
file(APPEND "${VIEW_SCRIPT}" "orientation 2.480970000000000e-01 4.765910000000000e-01 7.480970000000000e-01 3.894350000000000e-01;\n")
file(APPEND "${VIEW_SCRIPT}" "eye_pt 4.199070000000000e+01 2.940060000000000e+01 5.329140000000000e+01;\n")
file(APPEND "${VIEW_SCRIPT}" "start 0; clean;\n")
file(APPEND "${VIEW_SCRIPT}" "end;\n")

function(copy_snapshot local_image)
  if(NOT EXISTS "${local_image}")
    message(FATAL_ERROR "Expected image was not created: ${local_image}")
  endif()
  file(SIZE "${local_image}" image_size)
  if(image_size EQUAL 0)
    message(FATAL_ERROR "Rendered image is empty: ${local_image}")
  endif()
  file(COPY "${local_image}" DESTINATION "${IMGDIR}")
  get_filename_component(image_name "${local_image}" NAME)
  file(
    CHMOD "${IMGDIR}/${image_name}"
    PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
  )
endfunction()

function(render_snapshot name object)
  set(filename "solid_rep_${name}.png")
  set(local_image "${WORKDIR}/${filename}")
  file(REMOVE "${local_image}")
  execute_process(
    COMMAND "${RT}" -B -M -W "-s${SIZE}" -P1 -o "${filename}" "${PAWN_DB}" "${object}"
    WORKING_DIRECTORY "${WORKDIR}"
    TIMEOUT ${RENDER_TIMEOUT}
    INPUT_FILE "${VIEW_SCRIPT}"
    RESULT_VARIABLE render_result
    OUTPUT_VARIABLE render_stdout
    ERROR_VARIABLE render_stderr
  )
  if(NOT render_result EQUAL 0)
    message(
      FATAL_ERROR
      "Rendering ${name} failed (exit ${render_result})\n${render_stdout}${render_stderr}"
    )
  endif()
  copy_snapshot("${local_image}")
  message(STATUS "Rendered ${filename}")
endfunction()

render_snapshot(csg comparison.csg.r)
render_snapshot(brep_boolean comparison.csg.r.brep)
render_snapshot(brep_evaluated comparison.brep_evaluated.r)
render_snapshot(polygonal_boolean comparison.polygonal_boolean.r)
render_snapshot(nmg comparison.nmg.r)
render_snapshot(bot comparison.bot.r)
render_snapshot(revolve comparison.revolve.r)
render_snapshot(voxel_cells comparison.voxel_cells.r)
render_snapshot(volume comparison.volume.r)
render_snapshot(point_cloud comparison.point_cloud.r)

# MGED exports its current wireframe display as view-plane plot3.  Rasterize
# that plot through a disk framebuffer so regeneration needs no display.
set(WIREFRAME_PLOT "${WORKDIR}/solid_rep_wireframe.plot3")
set(WIREFRAME_PNG "${WORKDIR}/solid_rep_wireframe.png")
set(WIREFRAME_FB_NAME "solid_rep_wireframe.fb")
set(WIREFRAME_SCRIPT "${WORKDIR}/wireframe.tcl")
file(REMOVE "${WIREFRAME_PLOT}" "${WORKDIR}/${WIREFRAME_FB_NAME}" "${WIREFRAME_PNG}")
file(WRITE "${WIREFRAME_SCRIPT}" [=[
draw comparison.csg.r
size 80
orientation 0.248097 0.476591 0.748097 0.389435
eye_pt 41.9907 29.4006 53.2914
plot -2d solid_rep_wireframe.plot3
q
]=])
execute_process(
  COMMAND "${MGED}" --no-rc -c "${PAWN_DB}"
  WORKING_DIRECTORY "${WORKDIR}"
  TIMEOUT ${RENDER_TIMEOUT}
  INPUT_FILE "${WIREFRAME_SCRIPT}"
  RESULT_VARIABLE wireframe_result
  OUTPUT_VARIABLE wireframe_stdout
  ERROR_VARIABLE wireframe_stderr
)
if(NOT wireframe_result EQUAL 0 OR NOT EXISTS "${WIREFRAME_PLOT}")
  message(
    FATAL_ERROR
    "MGED wireframe export failed (exit ${wireframe_result})\n${wireframe_stdout}${wireframe_stderr}"
  )
endif()
message(STATUS "Exported wireframe plot")
run_process(
  "Cleared wireframe framebuffer"
  "${FBCLEAR}" -c -F "${WIREFRAME_FB_NAME}" -w "${SIZE}" -n "${SIZE}" 24 24 24
)
execute_process(
  COMMAND "${PLOT3FB}" -F "${WIREFRAME_FB_NAME}" -o -w "${SIZE}" -n "${SIZE}"
  WORKING_DIRECTORY "${WORKDIR}"
  TIMEOUT ${RENDER_TIMEOUT}
  INPUT_FILE "${WIREFRAME_PLOT}"
  RESULT_VARIABLE plot_result
  OUTPUT_VARIABLE plot_stdout
  ERROR_VARIABLE plot_stderr
)
if(NOT plot_result EQUAL 0)
  message(FATAL_ERROR "plot3-fb failed (exit ${plot_result})\n${plot_stdout}${plot_stderr}")
endif()
run_process(
  "Rasterized wireframe snapshot"
  "${FBPNG}" -F "${WIREFRAME_FB_NAME}" -w "${SIZE}" -n "${SIZE}"
  "solid_rep_wireframe.png"
)
copy_snapshot("${WIREFRAME_PNG}")

message(STATUS "Solid representation comparison: 11 snapshots regenerated")
