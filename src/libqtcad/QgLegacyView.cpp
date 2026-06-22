/*                    Q G L E G A C Y V I E W . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
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
/** @file QgLegacyView.cpp
 *
 * Neutral qtcad helpers for the staged legacy view handle.
 */

#include "common.h"

#include "qtcad/QgLegacyView.h"
#include "qtcad/QgLegacyViewBsg.h"

extern "C" {
#include "rt/view_legacy_bsg.h"
}

int
qg_legacy_view_dimensions_set(qg_legacy_view *view, int width, int height)
{
    return rt_view_dimensions_set_bsg(qg_legacy_view_to_bsg(view), width,
	    height);
}

int
qg_legacy_view_width_get(const qg_legacy_view *view)
{
    return rt_view_width_from_bsg(qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_height_get(const qg_legacy_view *view)
{
    return rt_view_height_from_bsg(qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_unit_conversion_set(qg_legacy_view *view,
	double local2base,
	double base2local)
{
    return rt_view_unit_conversion_set_bsg(qg_legacy_view_to_bsg(view),
	    local2base, base2local);
}

double
qg_legacy_view_local2base_get(const qg_legacy_view *view)
{
    return rt_view_local2base_from_bsg(qg_legacy_view_to_bsg(view));
}

double
qg_legacy_view_base2local_get(const qg_legacy_view *view)
{
    return rt_view_base2local_from_bsg(qg_legacy_view_to_bsg(view));
}

void
qg_legacy_view_info_get(const qg_legacy_view *view, struct rt_view_info *info)
{
    rt_view_info_from_bsg(info, qg_legacy_view_to_bsg(view));
}

const char *
qg_legacy_view_name_get(const qg_legacy_view *view)
{
    return rt_view_name_from_bsg(qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_aet_get(const qg_legacy_view *view, vect_t aet)
{
    return rt_view_aet_from_bsg(aet, qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_center_get(const qg_legacy_view *view, mat_t center)
{
    return rt_view_center_from_bsg(center, qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_unique_object_name(struct bu_vls *name, const char *seed,
	qg_legacy_view *view)
{
    return rt_view_unique_object_name_bsg(name, seed,
	    qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_framebuffer_mode_get(const qg_legacy_view *view)
{
    return rt_view_framebuffer_mode_from_bsg(qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_framebuffer_mode_set(qg_legacy_view *view, int mode)
{
    return rt_view_framebuffer_mode_set_bsg(qg_legacy_view_to_bsg(view), mode);
}

int
qg_legacy_view_lod_policy_get(const qg_legacy_view *view,
	struct rt_view_lod_policy *policy)
{
    return rt_view_lod_policy_from_bsg(policy, qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_lod_policy_apply(qg_legacy_view *view,
	const struct rt_view_lod_policy *policy)
{
    return rt_view_lod_policy_apply_bsg(qg_legacy_view_to_bsg(view), policy);
}

int
qg_legacy_view_adc_state_get(const qg_legacy_view *view,
	struct rt_view_adc_state *state)
{
    return rt_view_adc_state_from_bsg(state, qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_adc_state_set(qg_legacy_view *view,
	const struct rt_view_adc_state *state)
{
    return rt_view_adc_state_set_bsg(qg_legacy_view_to_bsg(view), state);
}

int
qg_legacy_view_center_dot_state_get(const qg_legacy_view *view,
	struct rt_view_other_state *state)
{
    return rt_view_center_dot_state_from_bsg(state,
	    qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_center_dot_state_set(qg_legacy_view *view,
	const struct rt_view_other_state *state)
{
    return rt_view_center_dot_state_set_bsg(qg_legacy_view_to_bsg(view),
	    state);
}

int
qg_legacy_view_grid_state_get(const qg_legacy_view *view,
	struct rt_view_grid_state *state)
{
    return rt_view_grid_state_from_bsg(state, qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_grid_state_set(qg_legacy_view *view,
	const struct rt_view_grid_state *state)
{
    return rt_view_grid_state_set_bsg(qg_legacy_view_to_bsg(view), state);
}

int
qg_legacy_view_model_axes_state_get(const qg_legacy_view *view,
	struct rt_view_axes_state *state)
{
    return rt_view_model_axes_state_from_bsg(state,
	    qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_model_axes_state_set(qg_legacy_view *view,
	const struct rt_view_axes_state *state)
{
    return rt_view_model_axes_state_set_bsg(qg_legacy_view_to_bsg(view),
	    state);
}

int
qg_legacy_view_scale_overlay_state_get(const qg_legacy_view *view,
	struct rt_view_other_state *state)
{
    return rt_view_scale_overlay_state_from_bsg(state,
	    qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_scale_overlay_state_set(qg_legacy_view *view,
	const struct rt_view_other_state *state)
{
    return rt_view_scale_overlay_state_set_bsg(qg_legacy_view_to_bsg(view),
	    state);
}

int
qg_legacy_view_view_axes_state_get(const qg_legacy_view *view,
	struct rt_view_axes_state *state)
{
    return rt_view_view_axes_state_from_bsg(state,
	    qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_view_axes_state_set(qg_legacy_view *view,
	const struct rt_view_axes_state *state)
{
    return rt_view_view_axes_state_set_bsg(qg_legacy_view_to_bsg(view),
	    state);
}

int
qg_legacy_view_params_state_get(const qg_legacy_view *view,
	struct rt_view_params_state *state)
{
    return rt_view_params_state_from_bsg(state, qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_params_state_set(qg_legacy_view *view,
	const struct rt_view_params_state *state)
{
    return rt_view_params_state_set_bsg(qg_legacy_view_to_bsg(view), state);
}

int
qg_legacy_view_snap_source_view_only_set(qg_legacy_view *view)
{
    return rt_view_snap_source_flags_set_bsg(qg_legacy_view_to_bsg(view),
	    RT_VIEW_SNAP_VIEW_BSG);
}

int
qg_legacy_view_snap_lines_get(const qg_legacy_view *view)
{
    return rt_view_snap_lines_from_bsg(qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_snap_lines_set(qg_legacy_view *view, int enabled)
{
    return rt_view_snap_lines_set_bsg(qg_legacy_view_to_bsg(view), enabled);
}

int
qg_legacy_view_snap_exclude_clear(qg_legacy_view *view)
{
    return rt_view_snap_exclude_feature_clear_bsg(qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_center_vec_set(qg_legacy_view *view, const point_t center)
{
    return rt_view_center_vec_set_bsg(qg_legacy_view_to_bsg(view), center);
}

int
qg_legacy_view_aet_set(qg_legacy_view *view, const vect_t aet)
{
    return rt_view_aet_set_bsg(qg_legacy_view_to_bsg(view), aet);
}

int
qg_legacy_view_update(qg_legacy_view *view)
{
    return rt_view_update_bsg(qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_screen_point_get(qg_legacy_view *view, point_t point,
	fastf_t sx, fastf_t sy)
{
    return rt_view_screen_point_from_bsg(point, qg_legacy_view_to_bsg(view),
	    sx, sy);
}

int
qg_legacy_view_polygon_ref_is_null(rt_view_polygon_ref ref)
{
    return rt_view_polygon_ref_is_null_bsg(ref);
}

int
qg_legacy_view_polygon_record_get(rt_view_polygon_ref ref,
	struct rt_view_polygon_record *record)
{
    return rt_view_polygon_record_get_bsg(ref, record);
}

rt_view_polygon_ref
qg_legacy_view_polygon_create(qg_legacy_view *view, int type,
	point_t *first_point)
{
    return rt_view_polygon_create_bsg(qg_legacy_view_to_bsg(view), type,
	    first_point);
}

rt_view_polygon_ref
qg_legacy_view_polygon_select(qg_legacy_view *view, point_t *current_point)
{
    return rt_view_polygon_select_bsg(qg_legacy_view_to_bsg(view),
	    current_point);
}

rt_view_polygon_ref
qg_legacy_view_polygon_find(qg_legacy_view *view, const char *name)
{
    return rt_view_polygon_find_bsg(qg_legacy_view_to_bsg(view), name);
}

rt_view_polygon_ref
qg_legacy_view_polygon_dup(qg_legacy_view *view, const char *name,
	const char *new_name)
{
    return rt_view_polygon_dup_bsg(qg_legacy_view_to_bsg(view), name, new_name);
}

void
qg_legacy_view_polygon_visit_records(qg_legacy_view *view,
	qg_legacy_view_polygon_record_callback_t callback,
	void *data)
{
    rt_view_polygon_visit_records_bsg(qg_legacy_view_to_bsg(view), callback,
	    data);
}

size_t
qg_legacy_view_polygon_snap_count(qg_legacy_view *view,
	rt_view_polygon_ref exclude)
{
    return rt_view_polygon_snap_count_bsg(qg_legacy_view_to_bsg(view),
	    exclude);
}

int
qg_legacy_view_polygon_clear_point_selection(qg_legacy_view *view)
{
    return rt_view_polygon_clear_point_selection_bsg(
	    qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_polygon_update(rt_view_polygon_ref ref, qg_legacy_view *view,
	int utype)
{
    return rt_view_polygon_update_bsg(ref, qg_legacy_view_to_bsg(view), utype);
}

int
qg_legacy_view_polygon_update_screen_point(rt_view_polygon_ref ref,
	qg_legacy_view *view,
	int x,
	int y,
	int utype)
{
    return rt_view_polygon_update_screen_pt_bsg(ref,
	    qg_legacy_view_to_bsg(view), x, y, utype);
}

int
qg_legacy_view_polygon_move(rt_view_polygon_ref ref, point_t *current_point,
	point_t *previous_point)
{
    return rt_view_polygon_move_bsg(ref, current_point, previous_point);
}

int
qg_legacy_view_polygon_set_name(rt_view_polygon_ref ref, const char *name)
{
    return rt_view_polygon_set_name_bsg(ref, name);
}

int
qg_legacy_view_polygon_set_view(rt_view_polygon_ref ref, qg_legacy_view *view)
{
    return rt_view_polygon_set_view_bsg(ref, qg_legacy_view_to_bsg(view));
}

int
qg_legacy_view_polygon_set_visual(rt_view_polygon_ref ref,
	const struct bu_color *edge_color,
	const struct bu_color *fill_color,
	fastf_t fill_slope_x,
	fastf_t fill_slope_y,
	fastf_t fill_density,
	fastf_t vZ,
	int fill_flag)
{
    return rt_view_polygon_set_visual_bsg(ref, edge_color, fill_color,
	    fill_slope_x, fill_slope_y, fill_density, vZ, fill_flag);
}

int
qg_legacy_view_polygon_set_open(rt_view_polygon_ref ref, int open)
{
    return rt_view_polygon_set_open_bsg(ref, open);
}

int
qg_legacy_view_polygon_close(rt_view_polygon_ref ref)
{
    return rt_view_polygon_close_bsg(ref);
}

int
qg_legacy_view_polygon_clear_selected_point(rt_view_polygon_ref ref)
{
    return rt_view_polygon_clear_selected_point_bsg(ref);
}

int
qg_legacy_view_polygon_remove(rt_view_polygon_ref ref)
{
    return rt_view_polygon_remove_bsg(ref);
}

void *
qg_legacy_view_polygon_user_data(rt_view_polygon_ref ref)
{
    return rt_view_polygon_user_data_bsg(ref);
}

int
qg_legacy_view_polygon_user_data_set(rt_view_polygon_ref ref, void *user_data)
{
    return rt_view_polygon_user_data_set_bsg(ref, user_data);
}

int
qg_legacy_view_polygon_csg(rt_view_polygon_ref target,
	rt_view_polygon_ref stencil,
	bg_clip_t op)
{
    return rt_view_polygon_csg_bsg(target, stencil, op);
}

rt_view_polygon_ref
qg_legacy_view_polygon_import_sketch(const char *name, struct db_i *dbip,
	struct directory *dp,
	qg_legacy_view *view)
{
    return rt_view_polygon_import_sketch_bsg(name, dbip, dp,
	    qg_legacy_view_to_bsg(view));
}

struct directory *
qg_legacy_view_polygon_export_sketch(struct db_i *dbip, const char *name,
	rt_view_polygon_ref ref)
{
    return rt_view_polygon_export_sketch_bsg(dbip, name, ref);
}

int
qg_legacy_view_polygon_snap_exclude_set(qg_legacy_view *view,
	rt_view_polygon_ref ref)
{
    return rt_view_polygon_snap_exclude_set_bsg(qg_legacy_view_to_bsg(view),
	    ref);
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
