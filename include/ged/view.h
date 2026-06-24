/*                        V I E W . H
 * BRL-CAD
 *
 * Copyright (c) 2008-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @addtogroup ged_view
 *
 * Geometry EDiting Library Database View Related Functions.
 *
 */
/** @{ */
/** @file ged/view.h */

#ifndef GED_VIEW_H
#define GED_VIEW_H

#include "common.h"
#include "bn/tol.h"
#include "bg/polygon.h"
#include "rt/db_fullpath.h"
#include "rt/db_instance.h"
#include "ged/defines.h"

__BEGIN_DECLS

struct rt_view_grid_state;
struct rt_view_info;
struct rt_view_interactive_rect_state;
struct rt_view_axes_state;
struct rt_view_other_state;
struct rt_view_params_state;
struct rt_view_adc_state;
struct rt_view_knob_values;
struct rt_view_lod_policy;

typedef void (*ged_view_context_update_callback_t)(void *view_ctx, void *data);

typedef enum ged_view_clear_flags {
    GED_VIEW_CLEAR_DB = 0x01,
    GED_VIEW_CLEAR_VIEW = 0x02,
    GED_VIEW_CLEAR_LOCAL = 0x04
} ged_view_clear_flags;

GED_EXPORT extern int ged_draw_scene_available(struct ged *gedp);

struct ged_polygon_export_state {
    fastf_t scale;
    point_t origin;
    mat_t rotation;
    mat_t view2model;
    mat_t model2view;
    struct bg_polygons polygons;
    fastf_t data_vZ;
};


/** Check if a drawable exists */
#define GED_CHECK_DRAWABLE(_gedp, _flags) \
    if (!ged_draw_scene_available(_gedp)) { \
	int ged_check_drawable_quiet = (_flags) & GED_QUIET; \
	if (!ged_check_drawable_quiet) { \
	    bu_vls_trunc((_gedp)->ged_result_str, 0); \
	    bu_vls_printf((_gedp)->ged_result_str, "A drawable does not exist."); \
	} \
	return (_flags); \
    }

/** Check if a view exists */
#define GED_CHECK_VIEW(_gedp, _flags) \
    if (ged_view_active_ctx(_gedp) == GED_VIEW_NULL) { \
	int ged_check_view_quiet = (_flags) & GED_QUIET; \
	if (!ged_check_view_quiet) { \
	    bu_vls_trunc((_gedp)->ged_result_str, 0); \
	    bu_vls_printf((_gedp)->ged_result_str, "A view does not exist."); \
	} \
	return (_flags); \
    }

GED_EXPORT extern int ged_export_polygon(struct ged *gedp, const struct ged_polygon_export_state *polygon_state, size_t polygon_i, const char *sname);
GED_EXPORT extern struct bg_polygon *ged_import_polygon(struct ged *gedp, const char *sname);
GED_EXPORT extern int ged_polygons_overlap(struct ged *gedp, struct bg_polygon *polyA, struct bg_polygon *polyB);
GED_EXPORT extern void ged_polygon_fill_segments(struct ged *gedp, struct bg_polygon *poly, vect2d_t vfilldir, fastf_t vfilldelta);


/**
 * Rotate angle degrees about the specified axis
 */
GED_EXPORT extern int ged_arot_args(struct ged *gedp, int argc, const char *argv[], mat_t rmat);

/**
 * Rotate the view.
 */
GED_EXPORT extern int ged_rot_args(struct ged *gedp, int argc, const char *argv[], char *coord, mat_t rmat);

/**
 * Scale the view.
 */
GED_EXPORT extern int ged_scale_args(struct ged *gedp, int argc, const char *argv[], fastf_t *sf1, fastf_t *sf2, fastf_t *sf3);

/**
 * Opaque view-context scale helpers.
 */
GED_EXPORT extern fastf_t ged_view_context_scale_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_scale_set(void *view_ctx, fastf_t scale);
GED_EXPORT extern int ged_view_context_update(void *view_ctx);
GED_EXPORT extern int ged_view_context_is_independent(const void *view_ctx);
GED_EXPORT extern int ged_view_context_independent_scope_is_null(void *view_ctx, int create);
GED_EXPORT extern void ged_view_context_independent_scope_destroy(void *view_ctx);
GED_EXPORT extern size_t ged_view_context_clear(void *view_ctx, int flags);
GED_EXPORT extern int ged_view_context_cleared_set(void *view_ctx, int cleared);
GED_EXPORT extern int ged_view_context_refresh_drawn_count_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_refresh_request(void *view_ctx, uint32_t flags);
GED_EXPORT extern const char *ged_view_context_name_get(const void *view_ctx);
GED_EXPORT extern void *ged_view_context_create_with_set(void *view_set_ctx);
GED_EXPORT extern void *ged_view_context_create_copy_with_set(const void *src_view_ctx, void *view_set_ctx);
GED_EXPORT extern int ged_view_set_context_add(void *view_set_ctx, void *view_ctx);
GED_EXPORT extern int ged_view_context_update_callback_set(void *view_ctx, ged_view_context_update_callback_t callback, void *data);
GED_EXPORT extern void *ged_view_context_display_manager_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_display_manager_set(void *view_ctx, void *dmp);
GED_EXPORT extern int ged_view_context_width_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_height_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_dimensions_set(void *view_ctx, int width, int height);
GED_EXPORT extern fastf_t *ged_view_context_scale_storage_get(void *view_ctx);
GED_EXPORT extern int ged_view_context_framebuffer_mode_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_framebuffer_mode_set(void *view_ctx, int mode);
GED_EXPORT extern fastf_t ged_view_context_perspective_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_perspective_set(void *view_ctx, fastf_t perspective);
GED_EXPORT extern fastf_t ged_view_context_size_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_size_set(void *view_ctx, fastf_t size);
GED_EXPORT extern fastf_t ged_view_context_inverse_size_get(const void *view_ctx);
GED_EXPORT extern void ged_view_context_eye_pos_get(point_t eye_pos, const void *view_ctx);
GED_EXPORT extern int ged_view_context_eye_pos_set(void *view_ctx, const point_t eye_pos);
GED_EXPORT extern void ged_view_context_keypoint_get(point_t keypoint, const void *view_ctx);
GED_EXPORT extern int ged_view_context_keypoint_set(void *view_ctx, const point_t keypoint);
GED_EXPORT extern char ged_view_context_rotate_about_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_rotate_about_set(void *view_ctx, char rotate_about);
GED_EXPORT extern char ged_view_context_coord_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_coord_set(void *view_ctx, char coord);
GED_EXPORT extern int ged_view_context_zclip_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_zclip_set(void *view_ctx, int zclip);
GED_EXPORT extern int ged_view_context_autoview(void *view_ctx, fastf_t scale, int all_view_objs);
GED_EXPORT extern int ged_view_context_autoview_bounds(void *view_ctx, fastf_t scale, const point_t min, const point_t max);
GED_EXPORT extern int ged_view_context_screen_to_view(fastf_t *fx, fastf_t *fy, void *view_ctx, fastf_t x, fastf_t y);
GED_EXPORT extern int ged_view_context_screen_point(point_t point, void *view_ctx, fastf_t x, fastf_t y);
GED_EXPORT extern int ged_view_context_mouse_state_set(void *view_ctx, int x, int y);
GED_EXPORT extern void ged_view_context_model2view_get(mat_t model2view, const void *view_ctx);
GED_EXPORT extern void ged_view_context_view2model_get(mat_t view2model, const void *view_ctx);
GED_EXPORT extern int ged_view_context_view2model_set(void *view_ctx, const mat_t view2model);
GED_EXPORT extern void ged_view_context_pmodel2view_get(mat_t pmodel2view, const void *view_ctx);
GED_EXPORT extern void ged_view_context_pmat_get(mat_t pmat, const void *view_ctx);
GED_EXPORT extern int ged_view_context_pmat_set(void *view_ctx, const mat_t pmat);
GED_EXPORT extern void ged_view_context_center_get(mat_t center, const void *view_ctx);
GED_EXPORT extern int ged_view_context_center_vec_set(void *view_ctx, const point_t center);
GED_EXPORT extern void ged_view_context_rotation_get(mat_t rotation, const void *view_ctx);
GED_EXPORT extern int ged_view_context_rotation_set(void *view_ctx, const mat_t rotation);
GED_EXPORT extern void ged_view_context_aet_get(vect_t aet, const void *view_ctx);
GED_EXPORT extern int ged_view_context_aet_set(void *view_ctx, const vect_t aet);
GED_EXPORT extern int ged_view_context_plane_get(plane_t *plane, const void *view_ctx);
GED_EXPORT extern void ged_view_context_info_get(struct rt_view_info *info, const void *view_ctx);
GED_EXPORT extern int ged_view_context_orientation_quat_get(quat_t orientation, const void *view_ctx);
GED_EXPORT extern int ged_view_context_adc_state_get(struct rt_view_adc_state *adc, const void *view_ctx);
GED_EXPORT extern int ged_view_context_adc_state_set(void *view_ctx, const struct rt_view_adc_state *adc);
GED_EXPORT extern int ged_view_context_knob_values_get(struct rt_view_knob_values *values, const void *view_ctx);
GED_EXPORT extern int ged_view_context_knobs_reset(void *view_ctx, int category);
GED_EXPORT extern int ged_view_context_knobs_calibrate(void *view_ctx);
GED_EXPORT extern int ged_view_context_knobs_cmd_process(vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, void *view_ctx, const char *cmd, fastf_t factor, char origin, int model_flag, int incr_flag);
GED_EXPORT extern int ged_view_context_knobs_translate(void *view_ctx, const vect_t tvec, int model_flag);
GED_EXPORT extern int ged_view_context_knobs_rotate(void *view_ctx, const vect_t rvec, char origin, char coords, const matp_t obj_rot, const pointp_t pvt_pt);
GED_EXPORT extern int ged_view_context_knobs_update_rate_flags(void *view_ctx);
GED_EXPORT extern int ged_view_context_grid_state_get(struct rt_view_grid_state *grid, const void *view_ctx);
GED_EXPORT extern int ged_view_context_grid_state_set(void *view_ctx, const struct rt_view_grid_state *grid);
GED_EXPORT extern int ged_view_context_model_axes_state_get(struct rt_view_axes_state *axes, const void *view_ctx);
GED_EXPORT extern int ged_view_context_model_axes_state_set(void *view_ctx, const struct rt_view_axes_state *axes);
GED_EXPORT extern int ged_view_context_view_axes_state_get(struct rt_view_axes_state *axes, const void *view_ctx);
GED_EXPORT extern int ged_view_context_view_axes_state_set(void *view_ctx, const struct rt_view_axes_state *axes);
GED_EXPORT extern int ged_view_context_center_dot_state_get(struct rt_view_other_state *center_dot, const void *view_ctx);
GED_EXPORT extern int ged_view_context_center_dot_state_set(void *view_ctx, const struct rt_view_other_state *center_dot);
GED_EXPORT extern int ged_view_context_scale_overlay_state_get(struct rt_view_other_state *scale_state, const void *view_ctx);
GED_EXPORT extern int ged_view_context_scale_overlay_state_set(void *view_ctx, const struct rt_view_other_state *scale_state);
GED_EXPORT extern int ged_view_context_params_state_get(struct rt_view_params_state *params, const void *view_ctx);
GED_EXPORT extern int ged_view_context_params_state_set(void *view_ctx, const struct rt_view_params_state *params);
GED_EXPORT extern int ged_view_context_snap_grid_2d(void *view_ctx, fastf_t *vx, fastf_t *vy);
GED_EXPORT extern int ged_view_context_interactive_rect_state_get(struct rt_view_interactive_rect_state *rect, const void *view_ctx);
GED_EXPORT extern int ged_view_context_interactive_rect_state_set(void *view_ctx, const struct rt_view_interactive_rect_state *rect);
GED_EXPORT extern int ged_view_context_snap_lines_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_snap_lines_set(void *view_ctx, int enabled);
GED_EXPORT extern double ged_view_context_snap_tolerance_factor_get(const void *view_ctx);
GED_EXPORT extern int ged_view_context_snap_tolerance_factor_set(void *view_ctx, double factor);
GED_EXPORT extern int ged_view_context_unit_conversion_set(void *view_ctx, fastf_t local2base, fastf_t base2local);
GED_EXPORT extern int ged_view_context_lod_policy_get(struct rt_view_lod_policy *policy, const void *view_ctx);
GED_EXPORT extern int ged_view_context_lod_policy_apply(void *view_ctx, const struct rt_view_lod_policy *policy);

/**
 * Translate the view.
 */
GED_EXPORT extern int ged_tra_args(struct ged *gedp, int argc, const char *argv[], char *coord, vect_t tvec);



/* Geometry EDiting Library Object Selection Functions.
 *
 *
 * TODO - both the drawn state and the selection state need a rethink, particularly
 * since they need to sync when an interface is visually indicating one or both
 * sets of information.
 *
 * There are three things we may want from either:
 *
 * 1.  "minimal" list of active entities - the shallowest paths that fully
 * express the set of active objects
 *
 * 2.  "input" list - what the user has actually specified to create the
 * current state
 *
 * 3.  "solids" list - the full paths to the solids that will be the sources of
 * the actual geometry shown in the scene.  A complication here is drawing
 * modes that evaluate the geometry, which will not have child solids but will
 * themselves be holding a temporary set of generated view data.
 *
 * 4.  "active" list - the set of all paths between specified/minimal and solids
 * that are implicit active
 *
 * Both #1 and #3 can be generated from #2, although how quickly this can be
 * done is a concern as hierarchies get large.  #3 is necessary for graphical
 * interrogation of scenes to build selection sets, since it is those solid
 * objects that the user will be interacting with.  (In the case of a selection
 * set built solely from graphical selection #2 and #3 will be the same.
 * However, since the drawn scene will more likely have been specified with one
 * or a small number of higher level objects, the #2 draw list needs to be used
 * to generate a #3 list to support building up the selection list.)
 *
 * When a "who" command is run the user is probably looking for #1, but that is
 * not as clear in the select case - in a tree view, selecting a comb and all
 * of the children of that comb have different implications for editing.  The
 * former, unless the comb is a top level object, implies editing the instance of
 * the comb in its own parent.  The latter implies editing the comb definition,
 * altering the position of its children.  Unlike the former, the latter will
 * impact ALL uses of comb, not just a single instance.  That would imply select
 * should default to #2, but that's most likely not what is desired if the user
 * has selected large numbers of individual solids from #3 - they MAY want that,
 * but they may also be seeking to specify groups of regions or higher level objects
 * to manipulate.  That implies the need to expose some "collapse" operations to
 * the user so they may identify wholly selected parent objects and "promote"
 * to them.
 *
 * We also need an efficient way for both drawing and selection codes to look up
 * and manipulate the state of their corresponding information in the other container.
 * When drawing, the interface may need to illuminate selected objects (which can
 * be selected even if not drawn, and so will need to "appear" selected.)  This requires
 * being able to interrogate the #3 select list from the drawing code to identify if
 * a given full path is selected.  Likewise, when a selection occurs, drawn objects
 * will need to have their illuminate state updated - the #2 select list will need
 * to generate a #3 select list and then update (illuminate or de-illuminate, as the
 * case may be) corresponding solids on the #3 drawn list.  All this needs to also
 * happen quickly, to handle large selections and de-selections on very large geometry
 * hierarchies.  The lists may also be invalidated by geometry editing operations,
 * and so will need to be quickly validated against the .g state or updated in response
 * to .g changes.
 *
 * A complication on the #2 lists is what happens if an erase or deselect
 * command removes part of a previously specified object.  In that case a new
 * #2 list must be generated which captures as much of the original semantic
 * intent of the user specification as possible - "expanding" the relevant
 * entries from the #2 list to their solids, removing the solids from the user
 * specified removal parents, and then collapsing the remaining solids back to
 * a minimal set.
 *
 * Another #2 alteration case is when a higher level specification "subsumes" previous
 * #2 paths into it.  The existing paths must be recognized as subpaths of the
 * specified path, and removed in favor of it.
 *
 * At the moment, we have explicit logic in the drawing code for handling the above
 * cases, using db_full_path top matching.  Syncing between the selection and drawing
 * codes is imperfect - draw does not currently check for selected status when
 * creating scene objects.
 *
 * Design thoughts - what if for both the ged drawn and selected lists we
 * stash unordered sets of path hashes to indicate activity, and then add/remove
 * those hashes based on specified paths and the related solid tree walks?  To
 * generate the #1 set from the #2 we get the parent of each #2, check if all the
 * parent's children are in #4, and if so add it to #4 and process its parent and
 * so on.  If we find a parent without all of its children in #4, remove it from #4
 * if present and report it as a #1.  To go the other way, do a full path tree walk
 * calculating hashes as we go and populating #4.  Any solid we reach that way is
 * part of #3.  If the view containers maintain unordered maps of hashes to scene
 * objects, the selection code can directly identify any active objects.  The
 * highlighting update pass would be one iteration to clear all flags, and then
 * a second over the #3 hashes from the selection to set illum flags on the
 * currently selected objects.  As long as the same hashes are used for both,
 * no further processing would be necessary to ID selection changes.
 *
 * Draw objects also still need to be synced with the .g state in response to
 * update notifications from the .g I/O callbacks, which report directory
 * pointers... we need a way to flag solids based on that info as well, which
 * probably means we still need to keep the db_full_path associated with the
 * scene object.  May want to re-examine the current two-tier storage system
 * and just make all solid objs top level entries.  If we're willing to re-calculate
 * the "drawn" #1 list on the fly, and not worry about #2 for the drawn objects
 * case, that may simplify some things.
 */
 
/**
 * Return ged selections for specified object. Created if it doesn't
 * exist.
 */
GED_EXPORT struct rt_object_selections *ged_get_object_selections(struct ged *gedp,
                                                                  const char *object_name);

/**
 * Return ged selections of specified kind for specified object.
 * Created if it doesn't exist.
 */
GED_EXPORT struct rt_selection_set *ged_get_selection_set(struct ged *gedp,
                                                          const char *object_name,
                                                          const char *selection_name);



typedef void (*ged_drawable_notify_func_t)(int);

GED_EXPORT void
ged_dl_notify_func_set(struct ged *gedp, ged_drawable_notify_func_t f);
GED_EXPORT ged_drawable_notify_func_t
ged_dl_notify_func_get(struct ged *gedp);




/* NMG specific visualizations that (currently) need libged routines.
 *
 * This will almost certainly move elsewhere - its presence here should be
 * considered temporary and not relied on from an API design perspective.
 */
GED_EXPORT extern void nmg_plot_eu(struct ged *gedp, struct edgeuse *es_eu, const struct bn_tol *tol);


__END_DECLS

#endif /* GED_VIEW_H */

/** @} */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
