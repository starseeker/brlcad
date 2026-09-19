# Data-only rule inventory for BObol_pivot_guard.cmake.
#
# Permanent invariants have no removal date.  Temporary tombstones name an
# owner and removal condition here so migration history is not mixed into the
# guard implementation.  Keep remediation prose in the guard, not in regexes.

set(BOBOL_GUARD_SOURCE_EXTENSIONS
  [[(CMakeLists\.txt|\.(c|cc|cpp|cxx|h|hh|hpp|cmake)$)]])

# Permanent invariant: the drawing stack has no source dependency on retired
# libbsg or graphical libdm APIs.
set(BOBOL_GUARD_FORBIDDEN_DEPENDENCY_PATTERNS
  [[#[ \t]*include[ \t]*[<"]bsg]]
  [[#[ \t]*include[ \t]*[<"]dm]]
  [[(^|[^A-Za-z0-9_])libbsg([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])libdm([^A-Za-z0-9_]|$)]])

# Permanent invariant: endpoint/scene ownership may not be bypassed by old
# renderer objects or immediate-mode drawing entry points.
set(BOBOL_GUARD_FORBIDDEN_OWNERSHIP_PATTERNS
  [[(^|[^A-Za-z0-9_])bsg_view([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])bsg_render_item([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])bsg_backend_scene([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])bsg_vlist([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])bsg_vlblock([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])BSG_[A-Za-z0-9_]*VLIST[A-Za-z0-9_]*([^A-Za-z0-9_]|$)]]
  [[(^|[^A-Za-z0-9_])dm_draw_]]
  [[(^|[^A-Za-z0-9_])dm_plugins([^A-Za-z0-9_]|$)]])

# Permanent API tombstones.  These names represented unsafe draft-record or
# cross-renderer ownership and must never be reintroduced under the same API.
set(BOBOL_GUARD_RETIRED_GED_SYMBOL_PATTERNS
  [[(^|[^A-Za-z0-9_])ged_draw_shape_draft_[A-Za-z0-9_]*[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_create_evaluated_path_shape_ref[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_scene_ref_create_draft_pair[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_group_ref_append_scene_ref[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_scene_ref_append_source_owner_to_group[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_scene_ref_attach_source_state_record[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_source_data_free[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_scene_ref_apply_display_settings[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_obol_database_source_publish_primitive_wireframe_for_path_internal[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_obol_database_source_publish_annotation_record_from_internal_for_path[ \t\r\n]*\(]]
  [[(^|[^A-Za-z0-9_])ged_draw_bsg_appearance_from_neutral[ \t\r\n]*\(]])

# Temporary migration tombstones.
# Owner: libged drawing maintainers.
# Removal condition: delete this category after the first release containing
# the Obol-only drawing stack, once these historical test names can no longer
# occur in supported downstream branches.
set(BOBOL_GUARD_TEMPORARY_TEST_FILES
  src/libged/tests/draw/bsg_quad_stability.cpp
  src/libged/tests/draw/bsg_render_stability.cpp
  src/libged/tests/draw/mged_bsg.cpp
  src/libged/tests/draw/mged_shaded_mode_bsg.cpp
  src/libged/tests/draw/rtwizard_bsg.cpp
  src/libged/tests/draw/tcl_overlay_bsg.cpp)

set(BOBOL_GUARD_TEMPORARY_TEST_TOKENS
  [[ged_test_bsg_quad_stability]]
  [[ged_test_bsg_render_stability]]
  [[ged_test_mged_bsg]]
  [[ged_test_mged_shaded_mode_bsg]]
  [[ged_test_rtwizard_bsg]]
  [[ged_test_tcl_overlay_bsg]]
  [[bsg_quad_stability\.cpp]]
  [[bsg_render_stability\.cpp]]
  [[mged_bsg\.cpp]]
  [[mged_shaded_mode_bsg\.cpp]]
  [[rtwizard_bsg\.cpp]]
  [[tcl_overlay_bsg\.cpp]])
