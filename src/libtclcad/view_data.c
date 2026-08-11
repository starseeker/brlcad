/*                    V I E W _ D A T A . C
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file libtclcad/view_data.c
 *
 * Private TclCAD view-data initialization helpers.
 */

#include "common.h"

#include <string.h>

#include "bg/polygon.h"
#include "ged/view.h"
#include "./tclcad_private.h"

#define TCLCAD_DEFAULT_FONT_SIZE 20

static void
tclcad_polygon_state_init(tclcad_polygon_state *state)
{
    if (!state)
	return;

    memset(state, 0, sizeof(*state));
    state->gdps_clip_type = BG_POLYGON_BOOLEAN_UNION;
    MAT_ZERO(state->gdps_rotation);
    MAT_ZERO(state->gdps_view2model);
    MAT_ZERO(state->gdps_model2view);
}

static void
tclcad_view_state_init(tclcad_view_state *state)
{
    if (!state)
	return;

    memset(state, 0, sizeof(*state));
    tclcad_polygon_state_init(&state->gv_data_polygons);
    tclcad_polygon_state_init(&state->gv_sdata_polygons);
    state->gv_prim_labels.gos_font_size = TCLCAD_DEFAULT_FONT_SIZE;
}

void
tclcad_view_data_init(struct tclcad_view_data *tvd, struct ged *gedp)
{
    if (!tvd)
	return;

    bu_vls_init(&tvd->gdv_pathname);
    bu_vls_init(&tvd->gdv_edit_motion_delta_callback);
    tvd->gdv_edit_motion_delta_callback_cnt = 0;
    bu_vls_init(&tvd->gdv_callback);
    tvd->gdv_callback_cnt = 0;
    tvd->gedp = gedp;

    tclcad_view_state_init(&tvd->tcl_data);
}

static void
tclcad_arrow_state_clear(struct tclcad_data_arrow_state *state)
{
    if (state && state->gdas_points)
	bu_free(state->gdas_points, "TclCAD arrow points");
}

static void
tclcad_axes_state_clear(struct tclcad_data_axes_state *state)
{
    if (state && state->points)
	bu_free(state->points, "TclCAD axes points");
}

static void
tclcad_label_state_clear(tclcad_label_state *state)
{
    if (!state)
	return;
    if (state->gdls_labels) {
	for (int i = 0; i < state->gdls_num_labels; i++) {
	    if (state->gdls_labels[i])
		bu_free(state->gdls_labels[i], "TclCAD label");
	}
	bu_free(state->gdls_labels, "TclCAD labels");
    }
    if (state->gdls_points)
	bu_free(state->gdls_points, "TclCAD label points");
}

static void
tclcad_line_state_clear(struct tclcad_data_line_state *state)
{
    if (state && state->gdls_points)
	bu_free(state->gdls_points, "TclCAD line points");
}

static void
tclcad_polygon_state_clear(tclcad_polygon_state *state)
{
    if (!state)
	return;
    bg_polygons_clear(&state->gdps_polygons);
    if (state->gdps_styles)
	bu_free(state->gdps_styles, "TclCAD polygon styles");
}

void
tclcad_view_data_clear(struct tclcad_view_data *tvd)
{
    if (!tvd)
	return;
    tclcad_view_state *state = &tvd->tcl_data;
    tclcad_arrow_state_clear(&state->gv_data_arrows);
    tclcad_axes_state_clear(&state->gv_data_axes);
    tclcad_label_state_clear(&state->gv_data_labels);
    tclcad_line_state_clear(&state->gv_data_lines);
    tclcad_polygon_state_clear(&state->gv_data_polygons);
    tclcad_arrow_state_clear(&state->gv_sdata_arrows);
    tclcad_axes_state_clear(&state->gv_sdata_axes);
    tclcad_label_state_clear(&state->gv_sdata_labels);
    tclcad_line_state_clear(&state->gv_sdata_lines);
    tclcad_polygon_state_clear(&state->gv_sdata_polygons);
    memset(state, 0, sizeof(*state));
}

struct bu_vls *
tclcad_view_pathname_vls(const struct ged_view_context *view_ctx)
{
    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);
    return tvd && BU_VLS_IS_INITIALIZED(&tvd->gdv_pathname) ?
	&tvd->gdv_pathname : NULL;
}

struct tclcad_view_data *
tclcad_view_data_from_view_ctx(const struct ged_view_context *view_ctx)
{
    void *tcl_data = ged_view_context_tclcad_data_get(view_ctx);
    return tcl_data ? (struct tclcad_view_data *)((char *)tcl_data -
	offsetof(struct tclcad_view_data, tcl_data)) : NULL;
}

tclcad_view_state *
tclcad_view_tcl_data_from_view_ctx(struct ged_view_context *view_ctx)
{
    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);
    return tvd ? &tvd->tcl_data : NULL;
}

tclcad_polygon_state *
tclcad_view_polygon_state_from_view_ctx(struct ged_view_context *view_ctx, int staged)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!tcl_data)
	return NULL;

    return staged ? &tcl_data->gv_sdata_polygons : &tcl_data->gv_data_polygons;
}

int
tclcad_view_polygon_mode_from_view_ctx(struct ged_view_context *view_ctx)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    return tcl_data ? tcl_data->gv_polygon_mode : TCLCAD_IDLE_MODE;
}

int
tclcad_view_polygon_mode_set(struct ged_view_context *view_ctx, int mode)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!tcl_data)
	return 0;

    tcl_data->gv_polygon_mode = mode;
    return 1;
}

fastf_t
tclcad_view_data_vZ_from_view_ctx(struct ged_view_context *view_ctx)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    return tcl_data ? tcl_data->gv_data_vZ : 0.0;
}

int
tclcad_view_data_vZ_set(struct ged_view_context *view_ctx, fastf_t vZ)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!tcl_data)
	return 0;

    tcl_data->gv_data_vZ = vZ;
    return 1;
}

int
tclcad_view_hide_from_view_ctx(struct ged_view_context *view_ctx)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    return tcl_data ? tcl_data->gv_hide : 0;
}

int
tclcad_view_polygon_cflag_from_view_ctx(struct ged_view_context *view_ctx, int staged)
{
    tclcad_polygon_state *polygon_state =
	tclcad_view_polygon_state_from_view_ctx(view_ctx, staged);
    return polygon_state ? polygon_state->gdps_cflag : 0;
}

int
tclcad_view_polygon_cflag_clear(struct ged_view_context *view_ctx, int staged)
{
    tclcad_polygon_state *polygon_state =
	tclcad_view_polygon_state_from_view_ctx(view_ctx, staged);
    if (!polygon_state)
	return 0;

    polygon_state->gdps_cflag = 0;
    return 1;
}

tclcad_label_state *
tclcad_view_label_state_from_view_ctx(struct ged_view_context *view_ctx, int staged)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!tcl_data)
	return NULL;

    return staged ? &tcl_data->gv_sdata_labels : &tcl_data->gv_data_labels;
}

int
tclcad_view_prim_labels_state_from_view_ctx(struct bv_other_state *state, struct ged_view_context *view_ctx)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!state || !tcl_data)
	return 0;

    state->gos_draw = tcl_data->gv_prim_labels.gos_draw;
    for (int i = 0; i < 3; i++) {
	state->gos_line_color[i] = tcl_data->gv_prim_labels.gos_line_color[i];
	state->gos_text_color[i] = tcl_data->gv_prim_labels.gos_text_color[i];
    }
    state->gos_font_size = tcl_data->gv_prim_labels.gos_font_size;
    return 1;
}

int
tclcad_view_prim_labels_state_set(struct ged_view_context *view_ctx, const struct bv_other_state *state)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!state || !tcl_data)
	return 0;

    tcl_data->gv_prim_labels.gos_draw = state->gos_draw;
    for (int i = 0; i < 3; i++) {
	tcl_data->gv_prim_labels.gos_line_color[i] = state->gos_line_color[i];
	tcl_data->gv_prim_labels.gos_text_color[i] = state->gos_text_color[i];
    }
    tcl_data->gv_prim_labels.gos_font_size = state->gos_font_size;
    return 1;
}

int
tclcad_view_data_bind_view_ctx(struct ged_view_context *view_ctx, struct tclcad_view_data *tvd)
{
    return ged_view_context_tclcad_data_set(view_ctx,
	tvd ? (void *)&tvd->tcl_data : NULL);
}

void
tclcad_view_data_unbind_view_ctx(struct ged_view_context *view_ctx)
{
    (void)ged_view_context_tclcad_data_set(view_ctx, NULL);
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
