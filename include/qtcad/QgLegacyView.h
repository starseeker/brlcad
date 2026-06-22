/*                    Q G L E G A C Y V I E W . H
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
/** @file QgLegacyView.h
 *
 * Opaque transitional qtcad view handle.
 */

#ifndef QGLEGACYVIEW_H
#define QGLEGACYVIEW_H

#include "bg/polygon.h"
#include "qtcad/defines.h"
#include "rt/view.h"
#include "vmath.h"

struct bu_color;
struct bu_vls;
struct db_i;
struct directory;
struct qg_legacy_view;

typedef int (*qg_legacy_view_polygon_record_callback_t)(
	rt_view_polygon_ref ref,
	const struct rt_view_polygon_record *record,
	void *data);

QTCAD_EXPORT extern int qg_legacy_view_dimensions_set(qg_legacy_view *view,
	int width,
	int height);

QTCAD_EXPORT extern int qg_legacy_view_width_get(const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_height_get(const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_unit_conversion_set(qg_legacy_view *view,
	double local2base,
	double base2local);

QTCAD_EXPORT extern double qg_legacy_view_local2base_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern double qg_legacy_view_base2local_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern void qg_legacy_view_info_get(const qg_legacy_view *view,
	struct rt_view_info *info);

QTCAD_EXPORT extern const char *qg_legacy_view_name_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_aet_get(const qg_legacy_view *view,
	vect_t aet);

QTCAD_EXPORT extern int qg_legacy_view_center_get(const qg_legacy_view *view,
	mat_t center);

QTCAD_EXPORT extern int qg_legacy_view_unique_object_name(struct bu_vls *name,
	const char *seed,
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_framebuffer_mode_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_framebuffer_mode_set(
	qg_legacy_view *view,
	int mode);

QTCAD_EXPORT extern int qg_legacy_view_lod_policy_get(
	const qg_legacy_view *view,
	struct rt_view_lod_policy *policy);

QTCAD_EXPORT extern int qg_legacy_view_lod_policy_apply(qg_legacy_view *view,
	const struct rt_view_lod_policy *policy);

QTCAD_EXPORT extern int qg_legacy_view_adc_state_get(
	const qg_legacy_view *view,
	struct rt_view_adc_state *state);

QTCAD_EXPORT extern int qg_legacy_view_adc_state_set(qg_legacy_view *view,
	const struct rt_view_adc_state *state);

QTCAD_EXPORT extern int qg_legacy_view_center_dot_state_get(
	const qg_legacy_view *view,
	struct rt_view_other_state *state);

QTCAD_EXPORT extern int qg_legacy_view_center_dot_state_set(
	qg_legacy_view *view,
	const struct rt_view_other_state *state);

QTCAD_EXPORT extern int qg_legacy_view_grid_state_get(
	const qg_legacy_view *view,
	struct rt_view_grid_state *state);

QTCAD_EXPORT extern int qg_legacy_view_grid_state_set(qg_legacy_view *view,
	const struct rt_view_grid_state *state);

QTCAD_EXPORT extern int qg_legacy_view_model_axes_state_get(
	const qg_legacy_view *view,
	struct rt_view_axes_state *state);

QTCAD_EXPORT extern int qg_legacy_view_model_axes_state_set(
	qg_legacy_view *view,
	const struct rt_view_axes_state *state);

QTCAD_EXPORT extern int qg_legacy_view_scale_overlay_state_get(
	const qg_legacy_view *view,
	struct rt_view_other_state *state);

QTCAD_EXPORT extern int qg_legacy_view_scale_overlay_state_set(
	qg_legacy_view *view,
	const struct rt_view_other_state *state);

QTCAD_EXPORT extern int qg_legacy_view_view_axes_state_get(
	const qg_legacy_view *view,
	struct rt_view_axes_state *state);

QTCAD_EXPORT extern int qg_legacy_view_view_axes_state_set(
	qg_legacy_view *view,
	const struct rt_view_axes_state *state);

QTCAD_EXPORT extern int qg_legacy_view_params_state_get(
	const qg_legacy_view *view,
	struct rt_view_params_state *state);

QTCAD_EXPORT extern int qg_legacy_view_params_state_set(qg_legacy_view *view,
	const struct rt_view_params_state *state);

QTCAD_EXPORT extern int qg_legacy_view_snap_source_view_only_set(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_snap_lines_get(
	const qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_snap_lines_set(qg_legacy_view *view,
	int enabled);

QTCAD_EXPORT extern int qg_legacy_view_snap_exclude_clear(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_center_vec_set(qg_legacy_view *view,
	const point_t center);

QTCAD_EXPORT extern int qg_legacy_view_aet_set(qg_legacy_view *view,
	const vect_t aet);

QTCAD_EXPORT extern int qg_legacy_view_update(qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_screen_to_view(
	qg_legacy_view *view,
	fastf_t *vx,
	fastf_t *vy,
	fastf_t sx,
	fastf_t sy);

QTCAD_EXPORT extern int qg_legacy_view_screen_point_get(
	qg_legacy_view *view,
	point_t point,
	fastf_t sx,
	fastf_t sy);

QTCAD_EXPORT extern int qg_legacy_view_view2model_get(
	const qg_legacy_view *view,
	mat_t view2model);

QTCAD_EXPORT extern int qg_legacy_view_model2view_get(
	const qg_legacy_view *view,
	mat_t model2view);

QTCAD_EXPORT extern int qg_legacy_view_polygon_ref_is_null(
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_polygon_record_get(
	rt_view_polygon_ref ref,
	struct rt_view_polygon_record *record);

QTCAD_EXPORT extern rt_view_polygon_ref qg_legacy_view_polygon_create(
	qg_legacy_view *view,
	int type,
	point_t *first_point);

QTCAD_EXPORT extern rt_view_polygon_ref qg_legacy_view_polygon_select(
	qg_legacy_view *view,
	point_t *current_point);

QTCAD_EXPORT extern rt_view_polygon_ref qg_legacy_view_polygon_find(
	qg_legacy_view *view,
	const char *name);

QTCAD_EXPORT extern rt_view_polygon_ref qg_legacy_view_polygon_dup(
	qg_legacy_view *view,
	const char *name,
	const char *new_name);

QTCAD_EXPORT extern void qg_legacy_view_polygon_visit_records(
	qg_legacy_view *view,
	qg_legacy_view_polygon_record_callback_t callback,
	void *data);

QTCAD_EXPORT extern size_t qg_legacy_view_polygon_snap_count(
	qg_legacy_view *view,
	rt_view_polygon_ref exclude);

QTCAD_EXPORT extern int qg_legacy_view_polygon_clear_point_selection(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_polygon_update(
	rt_view_polygon_ref ref,
	qg_legacy_view *view,
	int utype);

QTCAD_EXPORT extern int qg_legacy_view_polygon_update_screen_point(
	rt_view_polygon_ref ref,
	qg_legacy_view *view,
	int x,
	int y,
	int utype);

QTCAD_EXPORT extern int qg_legacy_view_polygon_move(
	rt_view_polygon_ref ref,
	point_t *current_point,
	point_t *previous_point);

QTCAD_EXPORT extern int qg_legacy_view_polygon_set_name(
	rt_view_polygon_ref ref,
	const char *name);

QTCAD_EXPORT extern int qg_legacy_view_polygon_set_view(
	rt_view_polygon_ref ref,
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_polygon_set_visual(
	rt_view_polygon_ref ref,
	const struct bu_color *edge_color,
	const struct bu_color *fill_color,
	fastf_t fill_slope_x,
	fastf_t fill_slope_y,
	fastf_t fill_density,
	fastf_t vZ,
	int fill_flag);

QTCAD_EXPORT extern int qg_legacy_view_polygon_set_open(
	rt_view_polygon_ref ref,
	int open);

QTCAD_EXPORT extern int qg_legacy_view_polygon_close(
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_polygon_clear_selected_point(
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_polygon_remove(
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern void *qg_legacy_view_polygon_user_data(
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_polygon_user_data_set(
	rt_view_polygon_ref ref,
	void *user_data);

QTCAD_EXPORT extern int qg_legacy_view_polygon_csg(
	rt_view_polygon_ref target,
	rt_view_polygon_ref stencil,
	bg_clip_t op);

QTCAD_EXPORT extern rt_view_polygon_ref qg_legacy_view_polygon_import_sketch(
	const char *name,
	struct db_i *dbip,
	struct directory *dp,
	qg_legacy_view *view);

QTCAD_EXPORT extern struct directory *qg_legacy_view_polygon_export_sketch(
	struct db_i *dbip,
	const char *name,
	rt_view_polygon_ref ref);

QTCAD_EXPORT extern int qg_legacy_view_polygon_snap_exclude_set(
	qg_legacy_view *view,
	rt_view_polygon_ref ref);

#endif /* QGLEGACYVIEW_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
