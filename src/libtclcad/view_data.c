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

#include <stddef.h>
#include <string.h>

#include "bsg/tcl_data.h"
#include "ged/view.h"
#include "./tclcad_private.h"

#define TCLCAD_LAYOUT_ASSERT(_expr, _name) \
    typedef char tclcad_layout_assert_ ## _name[(_expr) ? 1 : -1]
#define TCLCAD_LAYOUT_SIZE(_tcl_type, _bsg_type, _name) \
    TCLCAD_LAYOUT_ASSERT(sizeof(_tcl_type) == sizeof(_bsg_type), _name)
#define TCLCAD_LAYOUT_OFFSET(_tcl_type, _bsg_type, _member, _name) \
    TCLCAD_LAYOUT_ASSERT(offsetof(_tcl_type, _member) == offsetof(_bsg_type, _member), _name)

TCLCAD_LAYOUT_SIZE(struct tclcad_data_axes_state, struct bsg_data_axes_state, axes_size);
TCLCAD_LAYOUT_SIZE(struct tclcad_data_arrow_state, struct bsg_data_arrow_state, arrow_size);
TCLCAD_LAYOUT_SIZE(tclcad_label_state, struct bsg_data_label_state, label_size);
TCLCAD_LAYOUT_SIZE(struct tclcad_data_line_state, struct bsg_data_line_state, line_size);
TCLCAD_LAYOUT_SIZE(tclcad_polygon_state, bsg_data_polygon_state, polygon_size);
TCLCAD_LAYOUT_SIZE(struct bv_other_state, struct bsg_other_state, other_size);
TCLCAD_LAYOUT_SIZE(tclcad_view_state, struct bsg_data_tclcad, view_state_size);

TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_polygon_mode, view_polygon_mode);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_hide, view_hide);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_data_vZ, view_data_vz);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_data_arrows, view_data_arrows);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_data_axes, view_data_axes);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_data_labels, view_data_labels);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_data_lines, view_data_lines);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_data_polygons, view_data_polygons);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_sdata_arrows, view_sdata_arrows);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_sdata_axes, view_sdata_axes);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_sdata_labels, view_sdata_labels);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_sdata_lines, view_sdata_lines);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_sdata_polygons, view_sdata_polygons);
TCLCAD_LAYOUT_OFFSET(tclcad_view_state, struct bsg_data_tclcad, gv_prim_labels, view_prim_labels);

#define TCLCAD_DEFAULT_FONT_SIZE 20

static void
tclcad_polygon_state_init(tclcad_polygon_state *state)
{
    if (!state)
	return;

    memset(state, 0, sizeof(*state));
    state->gdps_clip_type = bg_Union;
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

    bu_vls_init(&tvd->gdv_edit_motion_delta_callback);
    tvd->gdv_edit_motion_delta_callback_cnt = 0;
    bu_vls_init(&tvd->gdv_callback);
    tvd->gdv_callback_cnt = 0;
    tvd->gedp = gedp;

    tclcad_view_state_init(&tvd->tcl_data);
}

struct tclcad_view_data *
tclcad_view_data_from_view_ctx(void *view_ctx)
{
    void *tcl_data = ged_view_context_tclcad_data_get(view_ctx);
    return tcl_data ? (struct tclcad_view_data *)((char *)tcl_data -
	offsetof(struct tclcad_view_data, tcl_data)) : NULL;
}

tclcad_view_state *
tclcad_view_tcl_data_from_view_ctx(void *view_ctx)
{
    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);
    return tvd ? &tvd->tcl_data : NULL;
}

tclcad_polygon_state *
tclcad_view_polygon_state_from_view_ctx(void *view_ctx, int staged)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!tcl_data)
	return NULL;

    return staged ? &tcl_data->gv_sdata_polygons : &tcl_data->gv_data_polygons;
}

int
tclcad_view_polygon_mode_from_view_ctx(void *view_ctx)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    return tcl_data ? tcl_data->gv_polygon_mode : TCLCAD_IDLE_MODE;
}

int
tclcad_view_polygon_mode_set(void *view_ctx, int mode)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!tcl_data)
	return 0;

    tcl_data->gv_polygon_mode = mode;
    return 1;
}

fastf_t
tclcad_view_data_vZ_from_view_ctx(void *view_ctx)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    return tcl_data ? tcl_data->gv_data_vZ : 0.0;
}

int
tclcad_view_data_vZ_set(void *view_ctx, fastf_t vZ)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!tcl_data)
	return 0;

    tcl_data->gv_data_vZ = vZ;
    return 1;
}

int
tclcad_view_hide_from_view_ctx(void *view_ctx)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    return tcl_data ? tcl_data->gv_hide : 0;
}

int
tclcad_view_polygon_cflag_from_view_ctx(void *view_ctx, int staged)
{
    tclcad_polygon_state *polygon_state =
	tclcad_view_polygon_state_from_view_ctx(view_ctx, staged);
    return polygon_state ? polygon_state->gdps_cflag : 0;
}

int
tclcad_view_polygon_cflag_clear(void *view_ctx, int staged)
{
    tclcad_polygon_state *polygon_state =
	tclcad_view_polygon_state_from_view_ctx(view_ctx, staged);
    if (!polygon_state)
	return 0;

    polygon_state->gdps_cflag = 0;
    return 1;
}

tclcad_label_state *
tclcad_view_label_state_from_view_ctx(void *view_ctx, int staged)
{
    tclcad_view_state *tcl_data = tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!tcl_data)
	return NULL;

    return staged ? &tcl_data->gv_sdata_labels : &tcl_data->gv_data_labels;
}

int
tclcad_view_prim_labels_state_from_view_ctx(struct bv_other_state *state, void *view_ctx)
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
tclcad_view_prim_labels_state_set(void *view_ctx, const struct bv_other_state *state)
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
tclcad_view_data_bind_view_ctx(void *view_ctx, struct tclcad_view_data *tvd)
{
    return ged_view_context_tclcad_data_set(view_ctx,
	tvd ? (void *)&tvd->tcl_data : NULL);
}

void
tclcad_view_data_unbind_view_ctx(void *view_ctx)
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
