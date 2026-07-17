/*                     C O M M A N D S . C
 * BRL-CAD
 *
 * Copyright (c) 2000-2026 United States Government as represented by
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
/** @addtogroup libtclcad */
/** @{ */
/** @file libtclcad/commands.c
 *
 * A quasi-object-oriented database interface.
 *
 * A GED object contains the attributes and methods for controlling a
 * BRL-CAD geometry edit object.
 *
 */
/** @} */

#include "common.h"

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <assert.h>

#include "png.h"

#include "tcl.h"

#include "bio.h"

#include "bn.h"
#include "bv.h"
#include "bu/cmd.h"
#include "bu/path.h"
#include "bu/process.h"
#include "bu/units.h"
#include "vmath.h"
#include "rt/db4.h"
#include "rt/geom.h"
#include "rt/view.h"
#include "wdb.h"
#include "raytrace.h"
#include "ged.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "ged/event_txn.h"
#include "tclcad.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BHostFactory.h"

// tclcad.h pulls in OpenNURBS in C++ compilation mode, which defines None,
// which will conflict with Tk.h's Xlib None if we include tk.h before tclcad.h
#include "tcl.h"
#ifdef HAVE_TK
#  include "tk.h"
#endif

#include "imgstream/fbserv.h"
#include "bg/lseg.h"

#include "icv/io.h"
#include "icv/ops.h"
#include "icv/crop.h"

#ifdef HAVE_GL_GL_H
#  include <GL/gl.h>
#endif

/* For the moment call internal libged functions - a cleaner
 * solution will be needed eventually */
#include "ged/draw.h"
#include "../libged/ged_private.h"

/* Private headers */
#include "brlcad_version.h"
#include "./tclcad_private.h"
#include "./view/view.h"
#include "./draw_view_move_helpers.h"

static bobol_display_endpoint_t *
tclcad_commands_endpoint(const void *view_ctx)
{
    return ged_view_context_display_endpoint_get(view_ctx);
}

/* Only replace Tcl bindings when the selected native host can actually
 * normalize and deliver events.  In particular, the WGL host intentionally
 * does not claim INPUT until its native adapter exists. */
static int
tclcad_commands_endpoint_input_enabled(const void *view_ctx)
{
    bobol_display_endpoint_t *endpoint =
	tclcad_commands_endpoint(view_ctx);
    return endpoint && (bobol_display_endpoint_host_capabilities(endpoint) &
	BOBOL_HOST_CAP_INPUT);
}

static int
tclcad_commands_endpoint_dimension_get(const void *view_ctx,
	const char *name)
{
    bobol_display_endpoint_t *endpoint = tclcad_commands_endpoint(view_ctx);
    struct bobol_endpoint_property_value property =
	BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    if (!endpoint || bobol_display_endpoint_property_get(endpoint, name,
	    &property) != BOBOL_ENDPOINT_PROPERTY_OK ||
	property.type != BOBOL_ENDPOINT_PROPERTY_UINT || !property.uint_value)
	return 0;
    return (int)property.uint_value;
}

static void
tclcad_commands_dimensions_sync(void *view_ctx)
{
    int width = tclcad_commands_endpoint_dimension_get(view_ctx,
	"endpoint.width");
    int height = tclcad_commands_endpoint_dimension_get(view_ctx,
	"endpoint.height");
    if (width > 0 && height > 0)
	bv_context_dimensions_set((struct bv_context *)view_ctx, width, height);
}

static int
tclcad_commands_width(const void *view_ctx)
{
    int width = tclcad_commands_endpoint_dimension_get(view_ctx,
	"endpoint.width");
    return width > 0 ? width :
	bv_context_width_get((const struct bv_context *)view_ctx);
}

static int
tclcad_commands_height(const void *view_ctx)
{
    int height = tclcad_commands_endpoint_dimension_get(view_ctx,
	"endpoint.height");
    return height > 0 ? height :
	bv_context_height_get((const struct bv_context *)view_ctx);
}

static int
tclcad_commands_endpoint_bool_get(const void *view_ctx, const char *name,
	int *value)
{
    bobol_display_endpoint_t *endpoint =
	tclcad_commands_endpoint(view_ctx);
    struct bobol_endpoint_property_value property =
	BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    if (!endpoint || !value ||
	bobol_display_endpoint_property_get(endpoint, name, &property) !=
	    BOBOL_ENDPOINT_PROPERTY_OK ||
	property.type != BOBOL_ENDPOINT_PROPERTY_BOOL)
	return 0;
    *value = property.bool_value ? 1 : 0;
    return 1;
}

static int
tclcad_commands_endpoint_bool_set(void *view_ctx, const char *name, int value)
{
    bobol_display_endpoint_t *endpoint =
	tclcad_commands_endpoint(view_ctx);
    struct bobol_endpoint_property_value property =
	BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    property.type = BOBOL_ENDPOINT_PROPERTY_BOOL;
    property.bool_value = value ? 1 : 0;
    return endpoint &&
	bobol_display_endpoint_property_set(endpoint, name, &property) ==
	    BOBOL_ENDPOINT_PROPERTY_OK;
}

static int
tclcad_obol_input_action(void *user_data, BObolInputAction action,
	const BObolInputEvent *event)
{
    void *view_ctx = user_data;
    if (!view_ctx || !event || event->type != BOBOL_INPUT_KEY_PRESS)
	return 0;

    int handled = 0;
    switch (action) {
	case BOBOL_ACTION_TOGGLE_ADC:
	case BOBOL_ACTION_TOGGLE_MODEL_AXES:
	case BOBOL_ACTION_TOGGLE_VIEW_AXES:
	    handled = bobol_display_endpoint_input_faceplate_toggle_apply(
		tclcad_commands_endpoint(view_ctx), view_ctx, action, NULL);
	    break;
	case BOBOL_ACTION_VIEW_2:
	case BOBOL_ACTION_VIEW_3:
	case BOBOL_ACTION_VIEW_4:
	case BOBOL_ACTION_VIEW_5:
	case BOBOL_ACTION_VIEW_6:
	case BOBOL_ACTION_VIEW_7:
	case BOBOL_ACTION_VIEW_FRONT:
	case BOBOL_ACTION_VIEW_TOP:
	case BOBOL_ACTION_VIEW_BOTTOM:
	case BOBOL_ACTION_VIEW_LEFT:
	case BOBOL_ACTION_VIEW_REAR:
	case BOBOL_ACTION_VIEW_RIGHT:
	    handled = bobol_input_view_orientation_apply(view_ctx, action);
	    break;
	default:
	    return 0;
	}
    if (handled)
	to_refresh_view(view_ctx);
    return handled;
}

static const char *
tclcad_commands_view_name(const void *view_ctx)
{
    const char *name = bv_context_name_get(
	    (const struct bv_context *)view_ctx);
    return name ? name : "";
}

static void
tclcad_commands_sync_dimensions(void *target_ctx, void *source_ctx)
{
    tclcad_commands_dimensions_sync(source_ctx);
    bv_context_dimensions_set((struct bv_context *)target_ctx,
	bv_context_width_get((const struct bv_context *)source_ctx),
	bv_context_height_get((const struct bv_context *)source_ctx));
}

static struct bv *
tclcad_commands_bv(void *view_ctx)
{
    return bv_context_view((struct bv_context *)view_ctx);
}

static const struct bv *
tclcad_commands_bv_const(const void *view_ctx)
{
    return bv_context_view_const((const struct bv_context *)view_ctx);
}

static unsigned long long
tclcad_commands_prepare_snap(void *view_ctx)
{
    struct bv *view = tclcad_commands_bv(view_ctx);
    int flags;

    if (!view)
	return 0ULL;

    flags = bv_snap_source_flags_get(view);
    (void)bv_snap_source_flags_set(view, flags | BV_SNAP_TCL);
    return bv_snap_kind_mask_get(view);
}

static int
tclcad_commands_snap_point_2d(void *view_ctx, fastf_t *vx, fastf_t *vy,
	unsigned long long kinds)
{
    return (kinds & BV_SNAP_KIND_GRID) ?
	bv_snap_grid_2d(tclcad_commands_bv_const(view_ctx), vx, vy) : 0;
}

static int to_base2local(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_bg(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_configure(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_constrain_rmode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_constrain_tmode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_copy(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_data_move(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_data_move_func(struct ged *gedp,
	void *gdvp,
	int argc,
	const char *argv[],
	const char *usage);
static int to_data_move_object_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_data_move_object_mode_func(struct ged *gedp,
	void *gdvp,
	int argc,
	const char *argv[],
	const char *usage);
static int to_data_move_point_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_data_move_point_mode_func(struct ged *gedp,
	void *gdvp,
	int argc,
	const char *argv[],
	const char *usage);
static int to_data_pick(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int
to_data_pick_func(struct ged *gedp,
	void *gdvp,
	int argc,
	const char *argv[],
	const char *usage);
static int to_data_vZ(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_dplot(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_bot_edge_split(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_bot_face_split(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_fontsize(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_fit_png_image(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_init_view_bindings(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_delete_view(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_handle_expose(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_hide_view(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_idle_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_light(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_list_views(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_local2base(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_lod(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_make(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_mirror(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_edit_motion_delta_callback(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_more_args_callback(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_move_arb_edge_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_move_arb_face_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_bot_move_pnt(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_bot_move_pnts(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_bot_move_pnt_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_bot_move_pnts_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_metaball_move_pnt_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_pipe_move_pnt_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_move_pnt_common(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_new_view(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_orotate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_oscale_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_otranslate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_paint_rect_area(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_pix(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_png(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_rect_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_rotate_arb_face_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_rotate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_rt_end_callback(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_rt_gettrees(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_protate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_pscale_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_ptranslate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_data_scale_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_scale_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_screen2model(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_screen2view(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_set_coord(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_snap_view(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_translate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_transparency(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_view_callback(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_view_win_size(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_view2screen(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_vmake(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_vslew(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_zbuffer(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);
static int to_zclip(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *usage,
	int maxargs);

/* Utility Functions */

static void to_rt_end_callback_internal(int aborted);

static void to_output_handler(struct ged *gedp, char *line);



static struct to_cmdtab ged_cmds[] = {
    {"3ptarb",	(char *)0, TO_UNLIMITED, to_more_args_func, ged_exec_3ptarb},
    {"adc",	"args", 7, to_view_func, ged_exec_adc},
    {"adjust",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_adjust},
    {"ae2dir",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_ae2dir},
    {"aet",	"[[-i] az el [tw]]", 6, to_view_func_plus, ged_exec_aet},
    {"analyze",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_analyze},
    {"annotate", (char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_annotate},
    {"pipe_append_pnt",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_pipe_append_pnt},
    {"arb",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_arb},
    {"arced",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_arced},
    {"arot",	"x y z angle", 6, to_view_func_plus, ged_exec_arot},
    {"art",	"art test", TO_UNLIMITED, to_view_func, ged_exec_art},
    {"attr",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_attr},
    {"bb",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bb},
    {"bev",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bev},
    {"bo",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bo},
    {"bot",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot},
    {"bot_condense",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_condense},
    {"bot_decimate",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_decimate},
    {"bot_dump",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_dump},
    {"bot_exterior",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_exterior},
    {"bot_face_fuse",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_face_fuse},
    {"bot_face_sort",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_face_sort},
    {"bot_flip",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_flip},
    {"bot_fuse",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_fuse},
    {"bot_merge",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_merge},
    {"bot_smooth",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_smooth},
    {"bot_split",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_split},
    {"bot_sync",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_sync},
    {"bot_vertex_fuse",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_bot_vertex_fuse},
    {"brep",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_brep},
    {"c",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_c},
    {"cat",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_cat},
    {"center",	"[x y z]", 5, to_view_func_plus, ged_exec_center},
    {"check",	(char *)0, TO_UNLIMITED, to_view_func, ged_exec_check},
    {"clear",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_clear},
    {"clone",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_clone},
    {"coil",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_coil},
    {"color",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_color},
    {"comb",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_comb},
    {"comb_color",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_comb_color},
    {"combmem",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_combmem},
    {"constraint", (char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_constraint},
    {"copyeval",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_copyeval},
    {"copymat",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_copymat},
    {"cpi",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_cpi},
    {"d",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_d},
    {"dbconcat",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_dbconcat},
    {"dbfind",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_dbfind},
    {"dbip",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_dbip}, // TODO - this needs to go away
    {"dbot_dump",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_dbot_dump},
    {"debug", 	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_debug},
    {"debugbu", 	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_debugbu},
    {"debugdir",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_debugdir},
    {"debuglib",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_debuglib},
    {"debugnmg",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_debugnmg},
    {"decompose",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_decompose},
    {"delay",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_delay},
    {"dm",	"[options] subcommand [args]", TO_UNLIMITED, to_pass_through_func, ged_exec_dm},
    {"dplot",	"dplot_logfile", 1, to_dplot, ged_exec_dplot},
    {"metaball_delete_pnt",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_metaball_delete_pnt},
    {"pipe_delete_pnt",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_pipe_delete_pnt},
    {"dir2ae",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_dir2ae},
    {"draw",	(char *)0, TO_UNLIMITED, to_autoview_func, ged_exec_draw},
    {"dump",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_dump},
    {"dup",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_dup},
    {"E",	(char *)0, TO_UNLIMITED, to_autoview_func, ged_exec_E},
    {"e",	(char *)0, TO_UNLIMITED, to_autoview_func, ged_exec_e},
    {"eac",	(char *)0, TO_UNLIMITED, to_autoview_func, ged_exec_eac},
    {"echo",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_echo},
    {"edarb",	(char *)0, TO_UNLIMITED, to_more_args_func, ged_exec_edarb},
    {"edcodes",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_edcodes},
    {"edcolor",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_edcolor},
    {"edcomb",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_edcomb},
    {"edit",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_edit},
    {"edmater",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_edmater},
    {"env",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_env},
    {"erase",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_erase},
    {"ert",	"[args]", TO_UNLIMITED, to_view_func, ged_exec_ert},
    {"ev",	(char *)0, TO_UNLIMITED, to_autoview_func, ged_exec_ev},
    {"expand",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_expand},
    {"eye",	"[x y z]", 5, to_view_func_plus, ged_exec_eye},
    {"eye_pos",	"[x y z]", 5, to_view_func_plus, ged_exec_eye_pos},
    {"eye_pt",	"[x y z]", 5, to_view_func_plus, ged_exec_eye_pt},
    {"exists",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_exists},
    {"facetize",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_facetize},
    {"voxelize",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_voxelize},
    {"fb2pix",  	"[-h -i -c] [-s squaresize] [-w width] [-n height] [file.pix]", TO_UNLIMITED, to_view_func, ged_exec_fb2pix},
    {"fbclear",  	"[r g b]", TO_UNLIMITED, to_view_func, ged_exec_fbclear},
    {"find_arb_edge",	"arb vx vy ptol", 5, to_view_func, ged_exec_find_arb_edge},
    {"find_bot_edge",	"bot vx vy", 5, to_view_func, ged_exec_find_bot_edge},
    {"find_bot_pnt",	"bot vx vy", 5, to_view_func, ged_exec_find_bot_pnt},
    {"find_pipe_pnt",	"pipe x y z", 6, to_view_func, ged_exec_find_pipe_pnt},
    {"form",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_form},
    {"fracture",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_fracture},
    {"g",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_g},
    {"garbage_collect",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_garbage_collect},
    {"gdiff",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_gdiff},
    {"get",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_get},
    {"get_autoview",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_get_autoview},
    {"get_bot_edges",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_get_bot_edges},
    {"get_comb",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_get_comb},
    {"get_eyemodel",	"vname", 2, to_view_func, ged_exec_get_eyemodel},
    {"get_type",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_get_type},
    {"glob",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_glob},
    {"gqa",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_gqa},
    {"graph",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_graph},
    {"grid",	"args", 6, to_view_func, ged_exec_grid},
    {"grid2model_lu",	"x y", 4, to_view_func_less, ged_exec_grid2model_lu},
    {"grid2view_lu",	"x y", 4, to_view_func_less, ged_exec_grid2view_lu},
    {"heal",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_heal},
    {"hide",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_hide},
    {"how",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_how},
    {"human",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_human},
    {"i",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_i},
    {"idents",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_idents},
    {"illum",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_illum},
    {"importFg4Section",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_importFg4Section},
    {"in",	(char *)0, TO_UNLIMITED, to_more_args_func, ged_exec_in},
    {"inside",	(char *)0, TO_UNLIMITED, to_more_args_func, ged_exec_inside},
    {"isize",	"vname", 2, to_view_func, ged_exec_isize},
    {"item",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_item},
    {"joint",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_joint},
    {"joint2",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_joint2},
    {"keep",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_keep},
    {"keypoint",	"[x y z]", 5, to_view_func, ged_exec_keypoint},
    {"kill",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_kill},
    {"killall",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_killall},
    {"killrefs",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_killrefs},
    {"killtree",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_killtree},
    {"l",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_l},
    {"lc",      "[-d|-s|-r] [-z] [-0|-1|-2|-3|-4|-5] [-f {FileName}] {GroupName}", TO_UNLIMITED, to_pass_through_func, ged_exec_lc},
    {"listeval",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_listeval},
    {"loadview",	"filename", 3, to_view_func, ged_exec_loadview},
    {"lod",	(char *)0, TO_UNLIMITED, to_lod, ged_exec_lod},
    {"log",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_log},
    {"lookat",	"x y z", 5, to_view_func_plus, ged_exec_lookat},
    {"ls",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_ls},
    {"lt",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_lt},
    {"m2v_point",	"x y z", 5, to_view_func, ged_exec_m2v_point},
    {"make_name",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_make_name},
    {"make_pnts",	(char *)0, TO_UNLIMITED, to_more_args_func, ged_exec_make_pnts},
    {"mat4x3pnt",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_mat4x3pnt},
    {"match",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_match},
    {"mater",	(char *)0, TO_UNLIMITED, to_more_args_func, ged_exec_mater},
    {"material",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_material},
    {"model2grid_lu",	"x y z", 5, to_view_func_less, ged_exec_model2grid_lu},
    {"model2view",	"vname", 2, to_view_func, ged_exec_model2view},
    {"model2view_lu",	"x y z", 5, to_view_func_less, ged_exec_model2view_lu},
    {"move_arb_edge",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_move_arb_edge},
    {"move_arb_face",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_move_arb_face},
    {"metaball_move_pnt",	(char *)0, TO_UNLIMITED, to_move_pnt_common, ged_exec_metaball_move_pnt},
    {"pipe_move_pnt",	(char *)0, TO_UNLIMITED, to_move_pnt_common, ged_exec_pipe_move_pnt},
    {"mouse_add_metaball_pnt",	"obj mx my", TO_UNLIMITED, to_mouse_append_pnt_common, ged_exec_mouse_add_metaball_pnt},
    {"mouse_append_pipe_pnt",	"obj mx my", TO_UNLIMITED, to_mouse_append_pnt_common, ged_exec_mouse_append_pipe_pnt},
    {"mouse_move_metaball_pnt",	"obj i mx my", TO_UNLIMITED, to_mouse_move_pnt_common, ged_exec_mouse_move_metaball_pnt},
    {"mouse_move_pipe_pnt",	"obj i mx my", TO_UNLIMITED, to_mouse_move_pnt_common, ged_exec_mouse_move_pipe_pnt},
    {"mouse_prepend_pipe_pnt",	"obj mx my", TO_UNLIMITED, to_mouse_append_pnt_common, ged_exec_mouse_prepend_pipe_pnt},
    {"mv",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_mv},
    {"mvall",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_mvall},
    {"nirt",	"[args]", TO_UNLIMITED, to_view_func, ged_exec_nirt},
    {"nmg_collapse",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_nmg_collapse},
    {"nmg_fix_normals",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_nmg_fix_normals},
    {"nmg_simplify",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_nmg_simplify},
    {"npush",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_npush},
    {"ocenter",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_ocenter},
    {"open",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_open},
    {"orient",	"quat", 6, to_view_func_plus, ged_exec_orient},
    {"orientation",	"quat", 6, to_view_func_plus, ged_exec_orientation},
    {"orotate",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_orotate},
    {"oscale",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_oscale},
    {"otranslate",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_otranslate},
    {"overlay",	(char *)0, TO_UNLIMITED, to_autoview_func, ged_exec_overlay},
    {"pathlist",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_pathlist},
    {"paths",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_paths},
    {"perspective",	"[angle]", 3, to_view_func_plus, ged_exec_perspective},
    {"pix2fb",  	"[options] [file.pix]", TO_UNLIMITED, to_view_func, ged_exec_pix2fb},
    {"plot",	"[options] file.pl", 16, to_view_func, ged_exec_plot},
    {"pmat",	"[mat]", 3, to_view_func, ged_exec_pmat},
    {"pmodel2view",	"vname", 2, to_view_func, ged_exec_pmodel2view},
    {"png2fb",  "[options] [file.png]", TO_UNLIMITED, to_view_func, ged_exec_png2fb},
    {"pngwf",	"[options] file.png", 16, to_view_func, ged_exec_pngwf},
    {"prcolor",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_prcolor},
    {"prefix",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_prefix},
    {"pipe_prepend_pnt",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_pipe_prepend_pnt},
    {"preview",	"[options] script", TO_UNLIMITED, to_view_context_func, ged_exec_preview},
    {"protate",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_protate},
    {"postscript", "[options] file.ps", 16, to_view_func, ged_exec_postscript},
    {"pscale",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_pscale},
    {"pset",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_pset},
    {"ptranslate",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_ptranslate},
    {"push",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_push},
    {"put",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_put},
    {"put_comb",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_put_comb},
    {"putmat",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_putmat},
    {"qray",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_qray},
    {"quat",	"a b c d", 6, to_view_func_plus, ged_exec_quat},
    {"qvrot",	"x y z angle", 6, to_view_func_plus, ged_exec_qvrot},
    {"r",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_r},
    {"rcodes",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_rcodes},
    {"rect",	"args", 6, to_view_func, ged_exec_rect},
    {"red",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_red},
    {"regdef",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_regdef},
    {"regions",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_regions},
    {"solid_report",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_solid_report},
    {"rfarb",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_rfarb},
    {"rm",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_rm},
    {"rmap",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_rmap},
    {"rmat",	"[mat]", 3, to_view_func, ged_exec_rmat},
    {"rmater",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_rmater},
    {"rot",	"[-m|-v] x y z", 6, to_view_func_plus, ged_exec_rot},
    {"rot_about",	"[e|k|m|v]", 3, to_view_func, ged_exec_rot_about},
    {"rot_point",	"x y z", 5, to_view_func, ged_exec_rot_point},
    {"rotate_arb_face",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_rotate_arb_face},
    {"rrt",	"[args]", TO_UNLIMITED, to_view_func, ged_exec_rrt},
    {"rselect",		(char *)0, TO_UNLIMITED, to_view_func, ged_exec_rselect},
    {"rt",	"[args]", TO_UNLIMITED, to_view_func, ged_exec_rt},
    {"rtabort",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_rtabort},
    {"rtarea",	"[args]", TO_UNLIMITED, to_view_func, ged_exec_rtarea},
    {"rtcheck",	"[args]", TO_UNLIMITED, to_view_func, ged_exec_rtcheck},
    {"rtedge",	"[args]", TO_UNLIMITED, to_view_func, ged_exec_rtedge},
    {"rtweight", "[args]", TO_UNLIMITED, to_view_func, ged_exec_rtweight},
    {"rtwizard", "[args]", TO_UNLIMITED, to_view_func, ged_exec_rtwizard},
    {"savekey",	"filename", 3, to_view_func, ged_exec_savekey},
    {"saveview", (char *)0, TO_UNLIMITED, to_view_func, ged_exec_saveview},
    {"sca",	"sf", 3, to_view_func_plus, ged_exec_sca},
    {"screengrab",	"imagename.ext", TO_UNLIMITED, to_view_context_func, ged_exec_screengrab},
    {"search",		(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_search},
    {"select",		(char *)0, TO_UNLIMITED, to_view_func, ged_exec_select},
    {"set_output_script",	"[script]", TO_UNLIMITED, to_pass_through_func, ged_exec_set_output_script},
    {"set_transparency",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_set_transparency},
    {"set_uplotOutputMode",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_set_uplotOutputMode},
    {"setview",	"x y z", 5, to_view_func_plus, ged_exec_setview},
    {"shaded_mode",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_shaded_mode},
    {"shader",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_shader},
    {"shells",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_shells},
    {"showmats",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_showmats},
    {"size",	"[size]", 3, to_view_func_plus, ged_exec_size},
    {"slew",	"x y [z]", 5, to_view_func_plus, ged_exec_slew},
    {"solids",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_solids},
    {"solids_on_ray",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_solids_on_ray},
    {"summary",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_summary},
    {"sv",	"x y [z]", 5, to_view_func_plus, ged_exec_sv},
    {"sync",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_sync},
    {"t",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_t},
    {"tire",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_tire},
    {"title",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_title},
    {"tol",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_tol},
    {"tops",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_tops},
    {"tra",	"[-m|-v] x y z", 6, to_view_func_plus, ged_exec_tra},
    {"track",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_track},
    {"tree",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_tree},
    {"unhide",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_unhide},
    {"units",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_units},
    {"v2m_point",	"x y z", 5, to_view_func, ged_exec_v2m_point},
    {"vdraw",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_vdraw},
    {"version",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_version},
    {"view",	"quat|ypr|aet|center|eye|size [args]", 7, to_view_func_plus, ged_exec_view},
    {"view2grid_lu",	"x y z", 5, to_view_func_less, ged_exec_view2grid_lu},
    {"view2model",	"", 2, to_view_func_less, ged_exec_view2model},
    {"view2model_lu",	"x y z", 5, to_view_func_less, ged_exec_view2model_lu},
    {"view2model_vec",	"x y z", 5, to_view_func_less, ged_exec_view2model_vec},
    {"viewdir",	"[-i]", 3, to_view_func_less, ged_exec_viewdir},
    {"vnirt",	"[args]", TO_UNLIMITED, to_view_func, ged_exec_vnirt},
    {"wcodes",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_wcodes},
    {"whatid",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_whatid},
    {"which_shader",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_which_shader},
    {"whichair",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_whichair},
    {"whichid",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_whichid},
    {"who",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_who},
    {"wmater",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_wmater},
    {"x",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_x},
    {"xpush",	(char *)0, TO_UNLIMITED, to_pass_through_func, ged_exec_xpush},
    {"ypr",	"yaw pitch roll", 5, to_view_func_plus, ged_exec_ypr},
    {"zap",	(char *)0, TO_UNLIMITED, to_pass_through_and_refresh_func, ged_exec_zap},
    {"zoom",	"sf", 3, to_view_func_plus, ged_exec_zoom},
    {(char *)0,	(char *)0, 0, TO_WRAPPER_FUNC_PTR_NULL, GED_FUNC_PTR_NULL}
};


struct to_cmdtab to_cmds[] = {
    {"autoview",	"vname", TO_UNLIMITED, to_autoview, GED_FUNC_PTR_NULL},
    {"base2local",	(char *)0, TO_UNLIMITED, to_base2local, GED_FUNC_PTR_NULL},
    {"bg",	"[r g b]", TO_UNLIMITED, to_bg, GED_FUNC_PTR_NULL},
    {"blast",	(char *)0, TO_UNLIMITED, to_blast, GED_FUNC_PTR_NULL},
    {"bot_edge_split",	"bot face", TO_UNLIMITED, to_bot_edge_split, GED_FUNC_PTR_NULL},
    {"bot_face_split",	"bot face", TO_UNLIMITED, to_bot_face_split, GED_FUNC_PTR_NULL},
    {"bot_move_pnt",	(char *)0, TO_UNLIMITED, to_bot_move_pnt, GED_FUNC_PTR_NULL},
    {"bot_move_pnt_mode",	"obj i mx my", TO_UNLIMITED, to_bot_move_pnt_mode, GED_FUNC_PTR_NULL},
    {"bot_move_pnts",	(char *)0, TO_UNLIMITED, to_bot_move_pnts, GED_FUNC_PTR_NULL},
    {"bot_move_pnts_mode",	"mx my obj i1 [i2 ... iN]", TO_UNLIMITED, to_bot_move_pnts_mode, GED_FUNC_PTR_NULL},
    {"configure",	"vname", TO_UNLIMITED, to_configure, GED_FUNC_PTR_NULL},
    {"constrain_rmode",	"x|y|z x y", TO_UNLIMITED, to_constrain_rmode, GED_FUNC_PTR_NULL},
    {"constrain_tmode",	"x|y|z x y", TO_UNLIMITED, to_constrain_tmode, GED_FUNC_PTR_NULL},
    {"cp",	"[-f] from to", TO_UNLIMITED, to_copy, GED_FUNC_PTR_NULL},
    {"data_arrows",	"???", TO_UNLIMITED, to_data_arrows, GED_FUNC_PTR_NULL},
    {"data_axes",	"???", TO_UNLIMITED, to_data_axes, GED_FUNC_PTR_NULL},
    {"data_labels",	"???", TO_UNLIMITED, to_data_labels, GED_FUNC_PTR_NULL},
    {"data_lines",	"???", TO_UNLIMITED, to_data_lines, GED_FUNC_PTR_NULL},
    {"data_move",	"???", TO_UNLIMITED, to_data_move, GED_FUNC_PTR_NULL},
    {"data_move_object_mode",	"x y", TO_UNLIMITED, to_data_move_object_mode, GED_FUNC_PTR_NULL},
    {"data_move_point_mode",	"x y", TO_UNLIMITED, to_data_move_point_mode, GED_FUNC_PTR_NULL},
    {"data_pick",	"???", TO_UNLIMITED, to_data_pick, GED_FUNC_PTR_NULL},
    {"data_polygons",	"???", TO_UNLIMITED, to_data_polygons, GED_FUNC_PTR_NULL},
    {"data_scale_mode",	"x y", TO_UNLIMITED, to_data_scale_mode, GED_FUNC_PTR_NULL},
    {"data_vZ",	"[z]", TO_UNLIMITED, to_data_vZ, GED_FUNC_PTR_NULL},
    {"delete_view",	"vname", TO_UNLIMITED, to_delete_view, GED_FUNC_PTR_NULL},
    {"edit_motion_delta_callback",	"vname [args]", TO_UNLIMITED, to_edit_motion_delta_callback, GED_FUNC_PTR_NULL},
    {"faceplate",	"center_dot|prim_labels|view_params|view_scale color|draw [val(s)]", TO_UNLIMITED, to_faceplate, GED_FUNC_PTR_NULL},
    {"fit_png_image",	"image_file_in req_width req_height scale image_file_out", 6, to_fit_png_image, GED_FUNC_PTR_NULL},
    {"fontsize",	"[fontsize]", 3, to_fontsize, GED_FUNC_PTR_NULL},
    {"get_prev_mouse",	"vname", TO_UNLIMITED, to_get_prev_mouse, GED_FUNC_PTR_NULL},
    {"handle_expose",	"vname count", TO_UNLIMITED, to_handle_expose, GED_FUNC_PTR_NULL},
    {"hide_view",	"vname [0|1]", 3, to_hide_view, GED_FUNC_PTR_NULL},
    {"idle_mode",	"vname", TO_UNLIMITED, to_idle_mode, GED_FUNC_PTR_NULL},
    {"init_view_bindings",	"vname", TO_UNLIMITED, to_init_view_bindings, GED_FUNC_PTR_NULL},
    {"light",	"[0|1]", TO_UNLIMITED, to_light, GED_FUNC_PTR_NULL},
    {"list_views",	(char *)0, TO_UNLIMITED, to_list_views, GED_FUNC_PTR_NULL},
    {"listen",	"[port]", TO_UNLIMITED, to_listen, GED_FUNC_PTR_NULL},
    {"local2base",	(char *)0, TO_UNLIMITED, to_local2base, GED_FUNC_PTR_NULL},
    {"make",	(char *)0, TO_UNLIMITED, to_make, GED_FUNC_PTR_NULL},
    {"metaball_move_pnt_mode",	"obj pt_i mx my", TO_UNLIMITED, to_metaball_move_pnt_mode, GED_FUNC_PTR_NULL},
    {"mirror",	(char *)0, TO_UNLIMITED, to_mirror, GED_FUNC_PTR_NULL},
    {"model_axes",	"???", TO_UNLIMITED, to_model_axes, GED_FUNC_PTR_NULL},
    {"more_args_callback",	"set/get the \"more args\" callback", TO_UNLIMITED, to_more_args_callback, GED_FUNC_PTR_NULL},
    {"mouse_brep_selection_append", "obj mx my", 5, to_mouse_brep_selection_append, GED_FUNC_PTR_NULL},
    {"mouse_brep_selection_translate", "obj mx my", 5, to_mouse_brep_selection_translate, GED_FUNC_PTR_NULL},
    {"mouse_constrain_rot",	"coord mx my", TO_UNLIMITED, to_mouse_constrain_rot, GED_FUNC_PTR_NULL},
    {"mouse_constrain_trans",	"coord mx my", TO_UNLIMITED, to_mouse_constrain_trans, GED_FUNC_PTR_NULL},
    {"mouse_data_scale",	"mx my", TO_UNLIMITED, to_mouse_data_scale, GED_FUNC_PTR_NULL},
    {"mouse_find_arb_edge",	"obj mx my ptol", TO_UNLIMITED, to_mouse_find_arb_edge, GED_FUNC_PTR_NULL},
    {"mouse_find_bot_edge",	"obj mx my", TO_UNLIMITED, to_mouse_find_bot_edge, GED_FUNC_PTR_NULL},
    {"mouse_find_bot_pnt",	"obj mx my", TO_UNLIMITED, to_mouse_find_bot_pnt, GED_FUNC_PTR_NULL},
    {"mouse_find_metaball_pnt",	"obj mx my", TO_UNLIMITED, to_mouse_find_metaball_pnt, GED_FUNC_PTR_NULL},
    {"mouse_find_pipe_pnt",	"obj mx my", TO_UNLIMITED, to_mouse_find_pipe_pnt, GED_FUNC_PTR_NULL},
    {"mouse_joint_select", "obj mx my", 5, to_mouse_joint_select, GED_FUNC_PTR_NULL},
    {"mouse_joint_selection_translate", "obj mx my", 5, to_mouse_joint_selection_translate, GED_FUNC_PTR_NULL},
    {"mouse_move_arb_edge",	"obj edge mx my", TO_UNLIMITED, to_mouse_move_arb_edge, GED_FUNC_PTR_NULL},
    {"mouse_move_arb_face",	"obj face mx my", TO_UNLIMITED, to_mouse_move_arb_face, GED_FUNC_PTR_NULL},
    {"mouse_move_bot_pnt",	"[-r] obj i mx my", TO_UNLIMITED, to_mouse_move_bot_pnt, GED_FUNC_PTR_NULL},
    {"mouse_move_bot_pnts",	"mx my obj i1 [i2 ... iN]", TO_UNLIMITED, to_mouse_move_bot_pnts, GED_FUNC_PTR_NULL},
    {"mouse_orotate",	"obj mx my", TO_UNLIMITED, to_mouse_orotate, GED_FUNC_PTR_NULL},
    {"mouse_oscale",	"obj mx my", TO_UNLIMITED, to_mouse_oscale, GED_FUNC_PTR_NULL},
    {"mouse_otranslate",	"obj mx my", TO_UNLIMITED, to_mouse_otranslate, GED_FUNC_PTR_NULL},
    {"mouse_poly_circ",	"mx my", TO_UNLIMITED, to_mouse_poly_circ, GED_FUNC_PTR_NULL},
    {"mouse_poly_cont",	"mx my", TO_UNLIMITED, to_mouse_poly_cont, GED_FUNC_PTR_NULL},
    {"mouse_poly_ell",	"mx my", TO_UNLIMITED, to_mouse_poly_ell, GED_FUNC_PTR_NULL},
    {"mouse_poly_rect",	"mx my", TO_UNLIMITED, to_mouse_poly_rect, GED_FUNC_PTR_NULL},
    {"mouse_protate",	"obj attribute mx my", TO_UNLIMITED, to_mouse_protate, GED_FUNC_PTR_NULL},
    {"mouse_pscale",	"obj attribute mx my", TO_UNLIMITED, to_mouse_pscale, GED_FUNC_PTR_NULL},
    {"mouse_ptranslate",	"obj attribute mx my", TO_UNLIMITED, to_mouse_ptranslate, GED_FUNC_PTR_NULL},
    {"mouse_ray",	"mx my", TO_UNLIMITED, to_mouse_ray, GED_FUNC_PTR_NULL},
    {"mouse_pick_detail", "mx my [object]", TO_UNLIMITED, to_mouse_pick_detail, GED_FUNC_PTR_NULL},
    {"mouse_rect",	"mx my", TO_UNLIMITED, to_mouse_rect, GED_FUNC_PTR_NULL},
    {"mouse_rot",	"mx my", TO_UNLIMITED, to_mouse_rot, GED_FUNC_PTR_NULL},
    {"mouse_rotate_arb_face",	"obj face v mx my", TO_UNLIMITED, to_mouse_rotate_arb_face, GED_FUNC_PTR_NULL},
    {"mouse_scale",	"mx my", TO_UNLIMITED, to_mouse_scale, GED_FUNC_PTR_NULL},
    {"mouse_trans",	"mx my", TO_UNLIMITED, to_mouse_trans, GED_FUNC_PTR_NULL},
    {"move_arb_edge_mode",	"obj edge x y", TO_UNLIMITED, to_move_arb_edge_mode, GED_FUNC_PTR_NULL},
    {"move_arb_face_mode",	"obj face x y", TO_UNLIMITED, to_move_arb_face_mode, GED_FUNC_PTR_NULL},
    {"new_view",	"vname type [args]", TO_UNLIMITED, to_new_view, GED_FUNC_PTR_NULL},
    {"orotate_mode",	"obj x y", TO_UNLIMITED, to_orotate_mode, GED_FUNC_PTR_NULL},
    {"oscale_mode",	"obj x y", TO_UNLIMITED, to_oscale_mode, GED_FUNC_PTR_NULL},
    {"otranslate_mode",	"obj x y", TO_UNLIMITED, to_otranslate_mode, GED_FUNC_PTR_NULL},
    {"paint_rect_area",	"vname", TO_UNLIMITED, to_paint_rect_area, GED_FUNC_PTR_NULL},
    {"pipe_pnt_mode",	"obj seg_i mx my", TO_UNLIMITED, to_pipe_move_pnt_mode, GED_FUNC_PTR_NULL},
    {"pix",	"file", TO_UNLIMITED, to_pix, GED_FUNC_PTR_NULL},
    {"png",	"file", TO_UNLIMITED, to_png, GED_FUNC_PTR_NULL},
    {"poly_circ_mode",	"x y", TO_UNLIMITED, to_poly_circ_mode, GED_FUNC_PTR_NULL},
    {"poly_cont_build",	"x y", TO_UNLIMITED, to_poly_cont_build, GED_FUNC_PTR_NULL},
    {"poly_cont_build_end",	"y", TO_UNLIMITED, to_poly_cont_build_end, GED_FUNC_PTR_NULL},
    {"poly_ell_mode",	"x y", TO_UNLIMITED, to_poly_ell_mode, GED_FUNC_PTR_NULL},
    {"poly_rect_mode",	"x y [s]", TO_UNLIMITED, to_poly_rect_mode, GED_FUNC_PTR_NULL},
    {"prim_label",	"[prim_1 prim_2 ... prim_N]", TO_UNLIMITED, to_prim_label, GED_FUNC_PTR_NULL},
    {"protate_mode",	"obj attribute x y", TO_UNLIMITED, to_protate_mode, GED_FUNC_PTR_NULL},
    {"pscale_mode",	"obj attribute x y", TO_UNLIMITED, to_pscale_mode, GED_FUNC_PTR_NULL},
    {"ptranslate_mode",	"obj attribute x y", TO_UNLIMITED, to_ptranslate_mode, GED_FUNC_PTR_NULL},
    {"rect_mode",	"x y", TO_UNLIMITED, to_rect_mode, GED_FUNC_PTR_NULL},
    {"redraw",	"obj", 2, to_redraw, GED_FUNC_PTR_NULL},
    {"refresh",	"vname", TO_UNLIMITED, to_refresh, GED_FUNC_PTR_NULL},
    {"refresh_all",	(char *)0, TO_UNLIMITED, to_refresh_all, GED_FUNC_PTR_NULL},
    {"refresh_on",	"[0|1]", TO_UNLIMITED, to_refresh_on, GED_FUNC_PTR_NULL},
    {"rotate_arb_face_mode",	"obj face v x y", TO_UNLIMITED, to_rotate_arb_face_mode, GED_FUNC_PTR_NULL},
    {"rotate_mode",	"x y", TO_UNLIMITED, to_rotate_mode, GED_FUNC_PTR_NULL},
    {"rt_end_callback",	"[args]", TO_UNLIMITED, to_rt_end_callback, GED_FUNC_PTR_NULL},
    {"rt_gettrees",	"[-i] [-u] pname object", TO_UNLIMITED, to_rt_gettrees, GED_FUNC_PTR_NULL},
    {"scale_mode",	"x y", TO_UNLIMITED, to_scale_mode, GED_FUNC_PTR_NULL},
    {"screen2model",	"x y", TO_UNLIMITED, to_screen2model, GED_FUNC_PTR_NULL},
    {"screen2view",	"x y", TO_UNLIMITED, to_screen2view, GED_FUNC_PTR_NULL},
    {"sdata_arrows",	"???", TO_UNLIMITED, to_data_arrows, GED_FUNC_PTR_NULL},
    {"sdata_axes",	"???", TO_UNLIMITED, to_data_axes, GED_FUNC_PTR_NULL},
    {"sdata_labels",	"???", TO_UNLIMITED, to_data_labels, GED_FUNC_PTR_NULL},
    {"sdata_lines",	"???", TO_UNLIMITED, to_data_lines, GED_FUNC_PTR_NULL},
    {"sdata_polygons",	"???", TO_UNLIMITED, to_data_polygons, GED_FUNC_PTR_NULL},
    {"set_coord",	"[m|v]", TO_UNLIMITED, to_set_coord, GED_FUNC_PTR_NULL},
    {"set_fb_mode",	"[mode]", TO_UNLIMITED, to_set_fb_mode, GED_FUNC_PTR_NULL},
    {"snap_view",	"vx vy", 4, to_snap_view, GED_FUNC_PTR_NULL},
    {"translate_mode",	"x y", TO_UNLIMITED, to_translate_mode, GED_FUNC_PTR_NULL},
    {"transparency",	"[val]", TO_UNLIMITED, to_transparency, GED_FUNC_PTR_NULL},
    {"view2screen",	"", 2, to_view2screen, GED_FUNC_PTR_NULL},
    {"view_axes",	"vname [args]", TO_UNLIMITED, to_view_axes, GED_FUNC_PTR_NULL},
    {"view_callback",	"vname [args]", TO_UNLIMITED, to_view_callback, GED_FUNC_PTR_NULL},
    {"view_win_size",	"[s] | [x y]", 4, to_view_win_size, GED_FUNC_PTR_NULL},
    {"vmake",	"pname ptype", TO_UNLIMITED, to_vmake, GED_FUNC_PTR_NULL},
    {"vslew",	"x y", TO_UNLIMITED, to_vslew, GED_FUNC_PTR_NULL},
    {"zbuffer",	"[0|1]", TO_UNLIMITED, to_zbuffer, GED_FUNC_PTR_NULL},
    {"zclip",	"[0|1]", TO_UNLIMITED, to_zclip, GED_FUNC_PTR_NULL},
    {(char *)0,	(char *)0, 0, TO_WRAPPER_FUNC_PTR_NULL, GED_FUNC_PTR_NULL}
};


/**
 * @brief create the Tcl command for to_open
 *
 */
TCLCAD_EXPORT int
Ged_Init(Tcl_Interp *interp)
{

    if (library_initialized(0))
	return TCL_OK;

    {
	const char *version_str = brlcad_version();
	tclcad_eval_noresult(interp, "set brlcad_version", 1, &version_str);
    }

    BU_LIST_INIT(&HeadTclcadObj.l);
    (void)Tcl_CreateCommand(interp, (const char *)"go_open", to_open_tcl,
	    (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    (void)library_initialized(1);

    return TCL_OK;
}


/**
 * @brief
 * Generic interface for database commands.
 *
 * @par Usage:
 * procname cmd ?args?
 *
 * @return result of ged command.
 */
static int
to_cmd(ClientData clientData,
	Tcl_Interp *interp,
	int argc,
	char **argv)
{
    struct to_cmdtab *ctp;
    struct tclcad_obj *top = (struct tclcad_obj *)clientData;
    Tcl_DString ds;
    int ret = BRLCAD_ERROR;

    Tcl_DStringInit(&ds);

    if (argc < 2) {
	Tcl_DStringAppend(&ds, "subcommand not specified; must be one of: ", -1);
	for (ctp = ged_cmds; ctp->to_name != (char *)NULL; ctp++) {
	    Tcl_DStringAppend(&ds, " ", -1);
	    Tcl_DStringAppend(&ds, ctp->to_name, -1);
	}
	for (ctp = to_cmds; ctp->to_name != (char *)NULL; ctp++) {
	    Tcl_DStringAppend(&ds, " ", -1);
	    Tcl_DStringAppend(&ds, ctp->to_name, -1);
	}
	Tcl_DStringAppend(&ds, "\n", -1);
	Tcl_DStringResult(interp, &ds);

	return TCL_ERROR;
    }

    current_top = top;

    for (ctp = to_cmds; ctp->to_name != (char *)0; ctp++) {
	if (BU_STR_EQUAL(ctp->to_name, argv[1])) {
	    struct ged *gedp = top->to_gedp;
	    ret = (*ctp->to_wrapper_func)(gedp, argc-1, (const char **)argv+1, ctp->to_func, ctp->to_usage, ctp->to_maxargs);
	    break;
	}
    }
    if (ctp->to_name == (char *)0) {
	for (ctp = ged_cmds; ctp->to_name != (char *)0; ctp++) {
	    if (BU_STR_EQUAL(ctp->to_name, argv[1])) {
		struct ged *gedp = top->to_gedp;
		ret = (*ctp->to_wrapper_func)(gedp, argc-1, (const char **)argv+1, ctp->to_func, ctp->to_usage, ctp->to_maxargs);
		break;
	    }
	}
    }

    /* Command not found. */
    if (ctp->to_name == (char *)0) {
	Tcl_DStringAppend(&ds, "unknown subcommand: ", -1);
	Tcl_DStringAppend(&ds, argv[1], -1);
	Tcl_DStringAppend(&ds, "; must be one of: ", -1);

	for (ctp = ged_cmds; ctp->to_name != (char *)NULL; ctp++) {
	    Tcl_DStringAppend(&ds, " ", -1);
	    Tcl_DStringAppend(&ds, ctp->to_name, -1);
	}

	for (ctp = to_cmds; ctp->to_name != (char *)NULL; ctp++) {
	    Tcl_DStringAppend(&ds, " ", -1);
	    Tcl_DStringAppend(&ds, ctp->to_name, -1);
	}

	Tcl_DStringAppend(&ds, "\n", -1);
	Tcl_DStringResult(interp, &ds);

	return TCL_ERROR;
    }

    Tcl_DStringAppend(&ds, bu_vls_addr(top->to_gedp->ged_result_str), -1);
    Tcl_DStringResult(interp, &ds);

    if (ret & BRLCAD_ERROR)
	return TCL_ERROR;

    return TCL_OK;
}


static void
free_path_edit_params(struct bu_hash_tbl *t)
{
    struct bu_hash_entry *entry = bu_hash_next(t, NULL);
    while (entry) {
	struct path_edit_params *pp = (struct path_edit_params *)bu_hash_value(entry, NULL);
	BU_PUT(pp, struct path_edit_params);
	entry = bu_hash_next(t, entry);
    }
}


static void
tclcad_view_host_destroy(struct tclcad_obj *top, void *view_ctx)
{
    if (!top || !top->to_gedp || !view_ctx)
	return;

    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(view_ctx);

    if (tvd) {
	(void)to_close_fbs(view_ctx);
    }

    struct bu_vls *view_command = tclcad_view_pathname_vls(view_ctx);
    if (view_command && bu_vls_strlen(view_command))
	Tcl_DeleteCommand(top->to_interp, bu_vls_cstr(view_command));
    bobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(view_ctx);
    if (endpoint)
	(void)bobol_display_endpoint_input_action_handler_clear_if(endpoint,
	    tclcad_obol_input_action, view_ctx);
    (void)ged_view_context_display_endpoint_set(view_ctx, NULL, 0);

    if (tvd) {
	bu_vls_free(&tvd->gdv_pathname);
	bu_vls_free(&tvd->gdv_edit_motion_delta_callback);
	bu_vls_free(&tvd->gdv_callback);
	tclcad_view_data_unbind_view_ctx(view_ctx);
	BU_PUT(tvd, struct tclcad_view_data);
    }
}


/**
 * @brief
 * Called by Tcl when the object is destroyed.
 */
void
to_deleteProc(ClientData clientData)
{
    struct tclcad_obj *top = (struct tclcad_obj *)clientData;
    BU_LIST_DEQUEUE(&top->l);

    if (current_top == top)
	current_top = TCLCAD_OBJ_NULL;

    if (top->to_gedp) {

	// Clean up the libtclcad view data.
	void *gdvp = NULL;
	struct bu_ptbl *views = ged_view_set_views_ctx(top->to_gedp);
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    gdvp = BU_PTBL_GET(views, i);
	    tclcad_view_host_destroy(top, gdvp);
	}

	// Clean up the other libtclcad data
	if (top->to_gedp->u_data) {
	    struct tclcad_ged_data *tgd = (struct tclcad_ged_data *)top->to_gedp->u_data;
	    bu_vls_free(&tgd->go_rt_end_callback);
	    bu_vls_free(&tgd->go_more_args_callback);
	    free_path_edit_params(tgd->go_edited_paths);
	    bu_hash_destroy(tgd->go_edited_paths);
	    if (tgd->go_prim_label_list) {
		for (int i = 0; i < tgd->go_prim_label_list_size; ++i)
		    bu_vls_free(&tgd->go_prim_label_list[i]);
		bu_free(tgd->go_prim_label_list, "prim_label");
	    }
	    BU_PUT(tgd, struct tclcad_ged_data);
	    top->to_gedp->u_data = NULL;
	}
	if (top->to_gedp->ged_io_data) {
	    struct tclcad_io_data *t_iod = (struct tclcad_io_data *)top->to_gedp->ged_io_data;
	    tclcad_destroy_io_data(t_iod);
	    top->to_gedp->ged_io_data = NULL;
	}

	// Got the libtclcad cleanup done, have libged do its up.
	ged_close(top->to_gedp);
    }

    bu_free((void *)top, "struct ged_obj");
}


/**
 * @brief
 * Create a command named "oname" in "interp" using "gedp" as its state.
 *
 */
int
to_create_cmd(Tcl_Interp *interp,
	struct tclcad_obj *top,	/* pointer to object */
	const char *oname)	/* object name */
{
    if (top == TCLCAD_OBJ_NULL) {
	Tcl_AppendResult(interp, "to_create_cmd ", oname, " failed", NULL);
	return TCL_ERROR;
    }

    /* Instantiate the newprocname, with clientData of top */
    /* Beware, returns a "token", not TCL_OK. */
    (void)Tcl_CreateCommand(interp, oname, (Tcl_CmdProc *)to_cmd,
	    (ClientData)top, to_deleteProc);

    /* Return new function name as result */
    Tcl_AppendResult(interp, oname, (char *)NULL);

    return TCL_OK;
}

/**
 * @brief
 * A TCL interface to wdb_fopen() and wdb_dbopen().
 *
 * @par Implicit return -
 * Creates a new TCL proc which responds to get/put/etc. arguments
 * when invoked.  clientData of that proc will be ged_obj pointer for
 * this instance of the database.  Easily allows keeping track of
 * multiple databases.
 *
 * @return wdb pointer, for more traditional C-style interfacing.
 *
 * @par Example -
 * set top [go_open .inmem inmem $dbip]
 *@n	.inmem get box.s
 *@n	.inmem close
 *
 *@n go_open db file "bob.g"
 *@n db get white.r
 *@n db close
 */
int
to_open_tcl(ClientData UNUSED(clientData),
	Tcl_Interp *interp,
	int argc,
	const char **argv)
{
    struct tclcad_obj *top = NULL;
    struct ged *gedp = NULL;
    const char *dbname = NULL;

    if (argc == 1) {
	/* get list of database objects */
	for (BU_LIST_FOR(top, tclcad_obj, &HeadTclcadObj.l))
	    Tcl_AppendResult(interp, bu_vls_addr(&top->to_gedp->go_name), " ", (char *)NULL);

	return TCL_OK;
    }

    if (argc < 3 || 4 < argc) {
	Tcl_AppendResult(interp, "\
		Usage: go_open\n\
		go_open newprocname file filename\n\
		go_open newprocname disk $dbip\n\
		go_open newprocname disk_append $dbip\n\
		go_open newprocname inmem $dbip\n\
		go_open newprocname inmem_append $dbip\n\
		go_open newprocname db filename\n\
		go_open newprocname filename\n",
		NULL);
	return TCL_ERROR;
    }

    /* Delete previous proc (if any) to release all that memory, first */
    (void)Tcl_DeleteCommand(interp, argv[1]);

    if (argc == 3 || BU_STR_EQUAL(argv[2], "db")) {
	if (argc == 3) {
	    dbname = argv[2];
	    gedp = ged_open("filename", dbname, 0);
	} else {
	    dbname = argv[3];
	    gedp = ged_open("db", dbname, 0);
	}
    } else {
	dbname = argv[3];
	gedp = ged_open(argv[2], dbname, 0);
    }

    if (gedp == GED_NULL) {
	Tcl_AppendResult(interp, "Unable to open geometry database: ", dbname, (char *)NULL);
	return TCL_ERROR;
    }
    const char *disable_events = Tcl_GetVar(interp, "tclcad_disable_events", TCL_GLOBAL_ONLY);
    if (disable_events && atoi(disable_events))
	ged_event_txn_disable(gedp);
    gedp->ged_interp = (void *)interp;

    disable_events = Tcl_GetVar(interp, "tclcad_disable_events",
	    TCL_GLOBAL_ONLY);
    if (disable_events && disable_events[0] &&
	    !BU_STR_EQUAL(disable_events, "0"))
	(void)ged_event_bulk_begin(gedp);

    /* Set the Tcl specific I/O handlers for asynchronous subprocess I/O */
    struct tclcad_io_data *t_iod = tclcad_create_io_data();
    t_iod->io_mode  = TCL_READABLE;
    t_iod->interp = interp;
    gedp->ged_io_data = (void *)t_iod;
    gedp->ged_create_io_handler = &tclcad_create_io_handler;
    gedp->ged_delete_io_handler = &tclcad_delete_io_handler;

    /* initialize tclcad_obj */
    BU_ALLOC(top, struct tclcad_obj);
    top->to_interp = interp;

    BU_ASSERT(gedp != NULL);
    top->to_gedp = gedp;

    top->to_gedp->ged_output_handler = to_output_handler;
    top->to_gedp->ged_refresh_handler = to_refresh_handler;

    ged_dl_notify_func_set(top->to_gedp, to_rt_end_callback_internal);

    // Initialize libtclcad GED data container
    struct tclcad_ged_data *tgd;
    BU_GET(tgd, struct tclcad_ged_data);
    bu_vls_init(&tgd->go_rt_end_callback);
    tgd->go_rt_end_callback_cnt = 0;
    bu_vls_init(&tgd->go_more_args_callback);
    tgd->go_more_args_callback_cnt = 0;
    tgd->go_edited_paths = bu_hash_create(0);
    tgd->go_prim_label_list = NULL;
    tgd->go_prim_label_list_size = 0;
    tgd->gedp = top->to_gedp;
    gedp->u_data = (void *)tgd;

    bu_vls_strcpy(&top->to_gedp->go_name, argv[1]);

    /* append to list of tclcad_obj */
    BU_LIST_APPEND(&HeadTclcadObj.l, &top->l);

    return to_create_cmd(interp, top, argv[1]);
}


/*************************** Local Command Functions ***************************/

static int
to_base2local(struct ged *gedp,
	int UNUSED(argc),
	const char *UNUSED(argv[]),
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    bu_vls_printf(gedp->ged_result_str, "%lf", current_top->to_gedp->dbip->dbi_base2local);

    return BRLCAD_OK;
}


static int
to_bg(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int r, g, b;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 2 && argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    bobol_display_endpoint_t *endpoint = tclcad_commands_endpoint(gdvp);
    if (!endpoint)
	return BRLCAD_ERROR;

    /* The endpoint controller is the rendering authority. */
    if (argc == 2) {
	struct bobol_endpoint_property_value property =
	    BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	if (bobol_display_endpoint_property_get(endpoint,
		"controller.background.bottom", &property) !=
	    BOBOL_ENDPOINT_PROPERTY_OK)
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str, "%d %d %d",
		(int)lrint(property.color3[0] * 255.0),
		(int)lrint(property.color3[1] * 255.0),
		(int)lrint(property.color3[2] * 255.0));
	return BRLCAD_OK;
    }

    /* set background color */
    if (bu_sscanf(argv[2], "%d", &r) != 1 ||
	    bu_sscanf(argv[3], "%d", &g) != 1 ||
	    bu_sscanf(argv[4], "%d", &b) != 1)
	goto bad_color;

    /* validate color */
    if (r < 0 || 255 < r ||
	    g < 0 || 255 < g ||
	    b < 0 || 255 < b)
	goto bad_color;

    struct bobol_endpoint_property_value property =
	BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    property.type = BOBOL_ENDPOINT_PROPERTY_COLOR3;
    VSET(property.color3, r / 255.0, g / 255.0, b / 255.0);
    if (bobol_display_endpoint_property_set(endpoint,
	    "controller.background.bottom", &property) !=
		BOBOL_ENDPOINT_PROPERTY_OK ||
	bobol_display_endpoint_property_set(endpoint,
	    "controller.background.top", &property) !=
		BOBOL_ENDPOINT_PROPERTY_OK)
	return BRLCAD_ERROR;

    to_refresh_view(gdvp);

    return BRLCAD_OK;

bad_color:
    bu_vls_printf(gedp->ged_result_str, "%s: %s %s %s", argv[0], argv[2], argv[3], argv[4]);
    return BRLCAD_ERROR;
}


static int
to_configure(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int status;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* configure the display manager window */
    tclcad_commands_dimensions_sync(gdvp);
    status = TCL_OK;

    {
	char cdimX[32];
	char cdimY[32];
	const char *av[5];
	void *active_view_ctx = ged_view_active_ctx(gedp);

	snprintf(cdimX, 32, "%d", tclcad_commands_width(gdvp));
	snprintf(cdimY, 32, "%d", tclcad_commands_height(gdvp));

	av[0] = "rect";
	av[1] = "cdim";
	av[2] = cdimX;
	av[3] = cdimY;
	av[4] = NULL;

	ged_view_active_ctx_set(gedp, gdvp);
	(void)ged_exec_rect(gedp, 4, (const char **)av);
	ged_view_active_ctx_set(gedp, active_view_ctx);

	/* A configure event updates only this view's geometry.  It must not make
	 * a newly created auxiliary pane the session framebuffer target. */
	if (active_view_ctx == gdvp)
	    (void)to_open_fbs(gdvp, (Tcl_Interp *)gedp->ged_interp);
    }

    if (status == TCL_OK) {
	to_refresh_view(gdvp);
	return BRLCAD_OK;
    }

    return BRLCAD_ERROR;
}


static int
to_constrain_rmode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if ((argv[2][0] != 'x' &&
		argv[2][0] != 'y' &&
		argv[2][0] != 'z') || argv[2][1] != '\0') {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_OK;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	    bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_CONSTRAINED_ROTATE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_constrain_rot %s %s %%x %%y}; break",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_constrain_tmode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if ((argv[2][0] != 'x' &&
		argv[2][0] != 'y' &&
		argv[2][0] != 'z') || argv[2][1] != '\0') {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_OK;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	    bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_CONSTRAINED_TRANSLATE_MODE);

    if (tclcad_view_pathname_vls(gdvp)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_constrain_trans %s %s %%x %%y}; break",
		bu_vls_addr(tclcad_view_pathname_vls(gdvp)),
		bu_vls_addr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2]);
	Tcl_Eval(current_top->to_interp, bu_vls_addr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_copy(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct ged *from_gedp = GED_NULL;
    struct ged *to_gedp = GED_NULL;
    int ret;
    const char *cp;
    struct tclcad_obj *top;
    struct bu_vls db_vls = BU_VLS_INIT_ZERO;
    struct bu_vls from_vls = BU_VLS_INIT_ZERO;
    struct bu_vls to_vls = BU_VLS_INIT_ZERO;
    int fflag;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 3 || 4 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    if (argc == 4) {
	if (argv[1][0] != '-' || argv[1][1] != 'f' ||  argv[1][2] != '\0') {
	    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	    return BRLCAD_ERROR;
	}

	fflag = 1;

	/* Advance past the -f option */
	--argc;
	++argv;
    } else
	fflag = 0;

    cp = strchr(argv[1], ':');
    if (cp) {
	bu_vls_strncpy(&db_vls, argv[1], cp-argv[1]);
	bu_vls_strcpy(&from_vls, cp+1);

	for (BU_LIST_FOR(top, tclcad_obj, &HeadTclcadObj.l)) {
	    if (BU_STR_EQUAL(bu_vls_addr(&top->to_gedp->go_name), bu_vls_addr(&db_vls))) {
		from_gedp = top->to_gedp;
		break;
	    }
	}

	bu_vls_free(&db_vls);

	if (from_gedp == GED_NULL) {
	    bu_vls_free(&from_vls);
	    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	    return BRLCAD_ERROR;
	}
    } else {
	bu_vls_strcpy(&from_vls, argv[1]);
	from_gedp = gedp;
    }

    cp = strchr(argv[2], ':');
    if (cp) {
	bu_vls_trunc(&db_vls, 0);
	bu_vls_strncpy(&db_vls, argv[2], cp-argv[2]);
	bu_vls_strcpy(&to_vls, cp+1);

	for (BU_LIST_FOR(top, tclcad_obj, &HeadTclcadObj.l)) {
	    if (BU_STR_EQUAL(bu_vls_addr(&top->to_gedp->go_name), bu_vls_addr(&db_vls))) {
		to_gedp = top->to_gedp;
		break;
	    }
	}

	bu_vls_free(&db_vls);

	if (to_gedp == GED_NULL) {
	    bu_vls_free(&from_vls);
	    bu_vls_free(&to_vls);
	    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	    return BRLCAD_ERROR;
	}
    } else {
	bu_vls_strcpy(&to_vls, argv[2]);
	to_gedp = gedp;
    }

    if (from_gedp == to_gedp) {
	ret = ged_dbcopy(from_gedp, to_gedp,
		bu_vls_addr(&from_vls),
		bu_vls_addr(&to_vls),
		fflag);

	if (ret != BRLCAD_OK && from_gedp != gedp)
	    bu_vls_strcpy(gedp->ged_result_str, bu_vls_addr(from_gedp->ged_result_str));
    } else {
	ret = ged_dbcopy(from_gedp, to_gedp,
		bu_vls_addr(&from_vls),
		bu_vls_addr(&to_vls),
		fflag);

	if (ret != BRLCAD_OK) {
	    if (bu_vls_strlen(from_gedp->ged_result_str)) {
		if (from_gedp != gedp)
		    bu_vls_strcpy(gedp->ged_result_str, bu_vls_addr(from_gedp->ged_result_str));
	    } else if (to_gedp != gedp && bu_vls_strlen(to_gedp->ged_result_str))
		bu_vls_strcpy(gedp->ged_result_str, bu_vls_addr(to_gedp->ged_result_str));
	}
    }

    bu_vls_free(&from_vls);
    bu_vls_free(&to_vls);

    return ret;
}


int
go_data_move(Tcl_Interp *UNUSED(interp),
	struct ged *gedp,
	void *draw_view_ctx,
	int argc,
	const char *argv[],
	const char *usage)
{
    void *gdvp = draw_view_ctx;
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 4 || 5 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    to_refresh_suppress_all_begin(current_top);

    int ret = to_data_move_func(gedp, gdvp, argc, argv, usage);
    to_refresh_suppress_all_end(current_top);
    return ret;
}


/*
 * Usage: data_move vname dtype dindex mx my
 */
static int
to_data_move(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 5 || 6 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* shift the command name to argv[1] before calling to_data_arrows_func */
    argv[1] = argv[0];
    return to_data_move_func(gedp, gdvp, argc-1, argv+1, usage);
}


static int
to_data_move_func(struct ged *gedp,
	void *gdvp,
	int argc,
	const char *argv[],
	const char *usage)
{
    int mx, my;
    int width, height;
    int dindex;
    fastf_t cx, cy;
    fastf_t vx, vy;
    fastf_t sf;
    point_t mpoint, vpoint;
    mat_t model2view, view2model;

    if (bu_sscanf(argv[2], "%d", &dindex) != 1 || dindex < 0)
	goto bad;

    if (argc == 4) {
	if (bu_sscanf(argv[3], "%d %d", &mx, &my) != 2)
	    goto bad;
    } else {
	if (bu_sscanf(argv[3], "%d", &mx) != 1)
	    goto bad;

	if (bu_sscanf(argv[4], "%d", &my) != 1)
	    goto bad;
    }

    width = tclcad_commands_width(gdvp);
    cx = 0.5 * (fastf_t)width;
    height = tclcad_commands_height(gdvp);
    cy = 0.5 * (fastf_t)height;
    sf = 2.0 / width;
    vx = (mx - cx) * sf;
    vy = (cy - my) * sf;
    bv_model2view_get(model2view, tclcad_commands_bv_const(gdvp));
    bv_view2model_get(view2model, tclcad_commands_bv_const(gdvp));

    if (BU_STR_EQUAL(argv[1], "data_polygons")) {
	size_t i, j, k;
	tclcad_polygon_state *gdpsp = tclcad_view_polygon_state_from_view_ctx(gdvp, 0);

	if (!gdpsp)
	    return BRLCAD_ERROR;

	if (bu_sscanf(argv[2], "%zu %zu %zu", &i, &j, &k) != 3)
	    goto bad;

	/* Silently ignore */
	if (i >= gdpsp->gdps_polygons.num_polygons ||
		j >= gdpsp->gdps_polygons.polygon[i].num_contours ||
		k >= gdpsp->gdps_polygons.polygon[i].contour[j].num_points)
	    return BRLCAD_OK;

	/* This section is for moving more than a single point on a contour */
	if (tclcad_view_polygon_mode_from_view_ctx(gdvp) == TCLCAD_DATA_MOVE_OBJECT_MODE) {
	    point_t old_mpoint, new_mpoint;
	    vect_t diff;

	    VMOVE(old_mpoint, gdpsp->gdps_polygons.polygon[i].contour[j].point[k]);

	    MAT4X3PNT(vpoint, model2view, gdpsp->gdps_polygons.polygon[i].contour[j].point[k]);
	    vpoint[X] = vx;
	    vpoint[Y] = vy;
	    MAT4X3PNT(new_mpoint, view2model, vpoint);
	    VSUB2(diff, new_mpoint, old_mpoint);

	    /* Move all polygons and all their respective contours. */
	    if (gdpsp->gdps_moveAll) {
		size_t p, c;
		for (p = 0; p < gdpsp->gdps_polygons.num_polygons; ++p) {
		    for (c = 0; c < gdpsp->gdps_polygons.polygon[p].num_contours; ++c) {
			for (k = 0; k < gdpsp->gdps_polygons.polygon[p].contour[c].num_points; ++k) {
			    VADD2(gdpsp->gdps_polygons.polygon[p].contour[c].point[k],
				    gdpsp->gdps_polygons.polygon[p].contour[c].point[k],
				    diff);
			}
		    }
		}
	    } else {
		/* Move only the contour. */
		for (k = 0; k < gdpsp->gdps_polygons.polygon[i].contour[j].num_points; ++k) {
		    VADD2(gdpsp->gdps_polygons.polygon[i].contour[j].point[k],
			    gdpsp->gdps_polygons.polygon[i].contour[j].point[k],
			    diff);
		}
	    }
	} else {
	    /* This section is for moving a single point on a contour */
	    MAT4X3PNT(vpoint, model2view, gdpsp->gdps_polygons.polygon[i].contour[j].point[k]);
	    vpoint[X] = vx;
	    vpoint[Y] = vy;
	    MAT4X3PNT(gdpsp->gdps_polygons.polygon[i].contour[j].point[k], view2model, vpoint);
	}

	to_refresh_view(gdvp);
	return BRLCAD_OK;
    }

	if (BU_STR_EQUAL(argv[1], "data_arrows")) {
	    /* Operate on typed draw-view data-arrow state instead of gv_tcl. */
	    const char *feature_name = "_tcl_data_arrows";

	    point_t *_pts = NULL;
	    int _npts = _tclcad_draw_view_data_arrows_points_copy(gdvp, feature_name, &_pts);

	    if (!_npts || dindex >= _npts) { if (_pts) bu_free(_pts, "TclCAD draw-view points"); return BRLCAD_OK; }

	if (tclcad_view_polygon_mode_from_view_ctx(gdvp) == TCLCAD_DATA_MOVE_OBJECT_MODE) {
	    int dindexA = dindex;
	    int dindexB = (dindex % 2) ? dindex - 1 : dindex + 1;

		    if (dindexB < 0 || dindexB >= _npts) { if (_pts) bu_free(_pts, "TclCAD draw-view points"); return BRLCAD_OK; }

	    point_t old_mpoint, new_mpoint;
	    vect_t diff;
	    VMOVE(old_mpoint, _pts[dindexA]);
	    MAT4X3PNT(vpoint, model2view, _pts[dindexA]);
	    vpoint[X] = vx; vpoint[Y] = vy;
	    MAT4X3PNT(new_mpoint, view2model, vpoint);
	    VSUB2(diff, new_mpoint, old_mpoint);
	    VMOVE(_pts[dindexA], new_mpoint);
	    VADD2(_pts[dindexB], _pts[dindexB], diff);
	} else {
	    MAT4X3PNT(vpoint, model2view, _pts[dindex]);
	    vpoint[X] = vx; vpoint[Y] = vy;
	    MAT4X3PNT(mpoint, view2model, vpoint);
	    VMOVE(_pts[dindex], mpoint);
	    }

	    int _color[3]; int _lw, _tl, _tw, _vis;
	    _tclcad_draw_view_data_arrows_style_read(gdvp, feature_name, _color, &_lw, &_tl, &_tw, &_vis);
	    _tclcad_draw_view_data_arrows_replace(gdvp, feature_name, _pts, _npts, _color, _lw, _tl, _tw, _vis);
	    bu_free(_pts, "TclCAD draw-view points");
	    to_refresh_view(gdvp);
	    return BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "sdata_arrows")) {
	    /* Operate on typed draw-view data-arrow state instead of gv_tcl. */
	    const char *feature_name = "_tcl_sdata_arrows";

	    point_t *_pts = NULL;
	    int _npts = _tclcad_draw_view_data_arrows_points_copy(gdvp, feature_name, &_pts);

	    if (!_npts || dindex >= _npts) { if (_pts) bu_free(_pts, "TclCAD draw-view points"); return BRLCAD_OK; }

	if (tclcad_view_polygon_mode_from_view_ctx(gdvp) == TCLCAD_DATA_MOVE_OBJECT_MODE) {
	    int dindexA = dindex;
	    int dindexB = (dindex % 2) ? dindex - 1 : dindex + 1;

		    if (dindexB < 0 || dindexB >= _npts) { if (_pts) bu_free(_pts, "TclCAD draw-view points"); return BRLCAD_OK; }

	    point_t old_mpoint, new_mpoint;
	    vect_t diff;
	    VMOVE(old_mpoint, _pts[dindexA]);
	    MAT4X3PNT(vpoint, model2view, _pts[dindexA]);
	    vpoint[X] = vx; vpoint[Y] = vy;
	    MAT4X3PNT(new_mpoint, view2model, vpoint);
	    VSUB2(diff, new_mpoint, old_mpoint);
	    VMOVE(_pts[dindexA], new_mpoint);
	    VADD2(_pts[dindexB], _pts[dindexB], diff);
	} else {
	    MAT4X3PNT(vpoint, model2view, _pts[dindex]);
	    vpoint[X] = vx; vpoint[Y] = vy;
	    MAT4X3PNT(mpoint, view2model, vpoint);
	    VMOVE(_pts[dindex], mpoint);
	    }

	    int _color[3]; int _lw, _tl, _tw, _vis;
	    _tclcad_draw_view_data_arrows_style_read(gdvp, feature_name, _color, &_lw, &_tl, &_tw, &_vis);
	    _tclcad_draw_view_data_arrows_replace(gdvp, feature_name, _pts, _npts, _color, _lw, _tl, _tw, _vis);
	    bu_free(_pts, "TclCAD draw-view points");
	    to_refresh_view(gdvp);
	    return BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "data_axes")) {
	    /* Extract typed draw-view center points, move one, rebuild. */
	    const char *feature_name = "_tcl_data_axes";

	    point_t *_cpts = NULL;
	    int _ncpts = _tclcad_draw_view_data_axes_centers_copy(gdvp, feature_name, &_cpts);

	    if (!_ncpts || dindex >= _ncpts) { if (_cpts) bu_free(_cpts, "TclCAD draw-view axes points"); return BRLCAD_OK; }

	MAT4X3PNT(vpoint, model2view, _cpts[dindex]);
	vpoint[X] = vx; vpoint[Y] = vy;
	MAT4X3PNT(mpoint, view2model, vpoint);
	VMOVE(_cpts[dindex], mpoint);

	    fastf_t _half = 1.0;
	    ged_draw_view_context_data_axes_half_size_get(gdvp, feature_name, &_half);

	    int _color[3]; int _lw, _vis;
	    _tclcad_draw_view_data_axes_style_read(gdvp, feature_name, _color, &_lw, &_vis);
	    _tclcad_draw_view_data_axes_replace(gdvp, feature_name, _cpts, _ncpts, _half, _color, _lw, _vis);
	    bu_free(_cpts, "TclCAD draw-view axes points");
	    to_refresh_view(gdvp);
	    return BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "sdata_axes")) {
	    /* Extract typed draw-view center points, move one, rebuild. */
	    const char *feature_name = "_tcl_sdata_axes";

	    point_t *_cpts = NULL;
	    int _ncpts = _tclcad_draw_view_data_axes_centers_copy(gdvp, feature_name, &_cpts);

	    if (!_ncpts || dindex >= _ncpts) { if (_cpts) bu_free(_cpts, "TclCAD draw-view axes points"); return BRLCAD_OK; }

	MAT4X3PNT(vpoint, model2view, _cpts[dindex]);
	vpoint[X] = vx; vpoint[Y] = vy;
	MAT4X3PNT(mpoint, view2model, vpoint);
	VMOVE(_cpts[dindex], mpoint);

	    fastf_t _half = 1.0;
	    ged_draw_view_context_data_axes_half_size_get(gdvp, feature_name, &_half);

	    int _color[3]; int _lw, _vis;
	    _tclcad_draw_view_data_axes_style_read(gdvp, feature_name, _color, &_lw, &_vis);
	    _tclcad_draw_view_data_axes_replace(gdvp, feature_name, _cpts, _ncpts, _half, _color, _lw, _vis);
	    bu_free(_cpts, "TclCAD draw-view axes points");
	    to_refresh_view(gdvp);
	    return BRLCAD_OK;
    }


	if (BU_STR_EQUAL(argv[1], "data_labels")) {
	    /* Modify typed draw-view label payloads through the data-label facade. */
	    const char *label_name = "_tcl_data_labels";
	    if ((size_t)dindex >= _tclcad_draw_view_data_labels_count(gdvp, label_name)) return BRLCAD_OK;

	    point_t _label_pt;
	    if (!_tclcad_draw_view_data_label_copy(gdvp, label_name, (size_t)dindex, NULL, _label_pt, NULL))
		return BRLCAD_OK;
	    MAT4X3PNT(vpoint, model2view, _label_pt);
	    vpoint[X] = vx; vpoint[Y] = vy;
	    MAT4X3PNT(mpoint, view2model, vpoint);
	    (void)_tclcad_draw_view_data_label_point_set(gdvp, label_name, (size_t)dindex, mpoint);

	to_refresh_view(gdvp);
	return BRLCAD_OK;
    }

	if (BU_STR_EQUAL(argv[1], "sdata_labels")) {
	    /* Modify typed draw-view label payloads through the data-label facade. */
	    const char *label_name = "_tcl_sdata_labels";
	    if ((size_t)dindex >= _tclcad_draw_view_data_labels_count(gdvp, label_name)) return BRLCAD_OK;

	    point_t _label_pt;
	    if (!_tclcad_draw_view_data_label_copy(gdvp, label_name, (size_t)dindex, NULL, _label_pt, NULL))
		return BRLCAD_OK;
	    MAT4X3PNT(vpoint, model2view, _label_pt);
	    vpoint[X] = vx; vpoint[Y] = vy;
	    MAT4X3PNT(mpoint, view2model, vpoint);
	    (void)_tclcad_draw_view_data_label_point_set(gdvp, label_name, (size_t)dindex, mpoint);

	to_refresh_view(gdvp);
	return BRLCAD_OK;
    }

	if (BU_STR_EQUAL(argv[1], "data_lines")) {
	    /* Operate on typed draw-view data-line points instead of gv_tcl. */
	    const char *feature_name = "_tcl_data_lines";

	    point_t *_pts = NULL;
	    int _npts = _tclcad_draw_view_data_lines_points_copy(gdvp, feature_name, &_pts);

	    if (!_npts || dindex >= _npts) { if (_pts) bu_free(_pts, "TclCAD draw-view points"); return BRLCAD_OK; }

	if (tclcad_view_polygon_mode_from_view_ctx(gdvp) == TCLCAD_DATA_MOVE_OBJECT_MODE) {
	    int dindexA = dindex;
	    int dindexB = (dindex % 2) ? dindex - 1 : dindex + 1;

		    if (dindexB < 0 || dindexB >= _npts) { if (_pts) bu_free(_pts, "TclCAD draw-view points"); return BRLCAD_OK; }

	    point_t old_mpoint, new_mpoint;
	    vect_t diff;
	    VMOVE(old_mpoint, _pts[dindexA]);
	    MAT4X3PNT(vpoint, model2view, _pts[dindexA]);
	    vpoint[X] = vx; vpoint[Y] = vy;
	    MAT4X3PNT(new_mpoint, view2model, vpoint);
	    VSUB2(diff, new_mpoint, old_mpoint);
	    VMOVE(_pts[dindexA], new_mpoint);
	    VADD2(_pts[dindexB], _pts[dindexB], diff);
	} else {
	    MAT4X3PNT(vpoint, model2view, _pts[dindex]);
	    vpoint[X] = vx; vpoint[Y] = vy;
	    MAT4X3PNT(mpoint, view2model, vpoint);
	    VMOVE(_pts[dindex], mpoint);
	    }

	    int _color[3]; int _lw, _vis;
	    _tclcad_draw_view_data_lines_style_read(gdvp, feature_name, _color, &_lw, &_vis);
	    _tclcad_draw_view_data_lines_replace(gdvp, feature_name, _pts, _npts, _color, _lw, _vis);
	    bu_free(_pts, "TclCAD draw-view points");
	    to_refresh_view(gdvp);
	    return BRLCAD_OK;
	}

	if (BU_STR_EQUAL(argv[1], "sdata_lines")) {
	    /* Operate on typed draw-view data-line points instead of gv_tcl. */
	    const char *feature_name = "_tcl_sdata_lines";

	    point_t *_pts = NULL;
	    int _npts = _tclcad_draw_view_data_lines_points_copy(gdvp, feature_name, &_pts);

	    if (!_npts || dindex >= _npts) { if (_pts) bu_free(_pts, "TclCAD draw-view points"); return BRLCAD_OK; }

	if (tclcad_view_polygon_mode_from_view_ctx(gdvp) == TCLCAD_DATA_MOVE_OBJECT_MODE) {
	    int dindexA = dindex;
	    int dindexB = (dindex % 2) ? dindex - 1 : dindex + 1;

		    if (dindexB < 0 || dindexB >= _npts) { if (_pts) bu_free(_pts, "TclCAD draw-view points"); return BRLCAD_OK; }

	    point_t old_mpoint, new_mpoint;
	    vect_t diff;
	    VMOVE(old_mpoint, _pts[dindexA]);
	    MAT4X3PNT(vpoint, model2view, _pts[dindexA]);
	    vpoint[X] = vx; vpoint[Y] = vy;
	    MAT4X3PNT(new_mpoint, view2model, vpoint);
	    VSUB2(diff, new_mpoint, old_mpoint);
	    VMOVE(_pts[dindexA], new_mpoint);
	    VADD2(_pts[dindexB], _pts[dindexB], diff);
	} else {
	    MAT4X3PNT(vpoint, model2view, _pts[dindex]);
	    vpoint[X] = vx; vpoint[Y] = vy;
	    MAT4X3PNT(mpoint, view2model, vpoint);
	    VMOVE(_pts[dindex], mpoint);
	    }

	    int _color[3]; int _lw, _vis;
	    _tclcad_draw_view_data_lines_style_read(gdvp, feature_name, _color, &_lw, &_vis);
	    _tclcad_draw_view_data_lines_replace(gdvp, feature_name, _pts, _npts, _color, _lw, _vis);
	    bu_free(_pts, "TclCAD draw-view points");
	    to_refresh_view(gdvp);
	    return BRLCAD_OK;
    }

bad:
    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
    return BRLCAD_ERROR;
}


int
go_data_move_object_mode(Tcl_Interp *UNUSED(interp),
	struct ged *gedp,
	void *draw_view_ctx,
	int argc,
	const char *argv[],
	const char *usage)
{
    void *gdvp = draw_view_ctx;
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    to_refresh_suppress_all_begin(current_top);

    int ret = to_data_move_object_mode_func(gedp, gdvp, argc, argv, usage);
    to_refresh_suppress_all_end(current_top);
    return ret;
}


static int
to_data_move_object_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* shift the command name to argv[1] before calling to_data_move_point_mode_func */
    argv[1] = argv[0];
    return to_data_move_object_mode_func(gedp, gdvp, argc-1, argv+1, usage);
}


static int
to_data_move_object_mode_func(struct ged *gedp,
	void *gdvp,
	int UNUSED(argc),
	const char *argv[],
	const char *usage)
{
    int x, y;

    ged_view_active_ctx_set(gedp, gdvp);

    if (bu_sscanf(argv[1], "%d", &x) != 1 ||
	    bu_sscanf(argv[2], "%d", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    /* At the moment, only the TclCAD polygon mode is being used. */
    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_DATA_MOVE_OBJECT_MODE);

    return BRLCAD_OK;
}


int
go_data_move_point_mode(Tcl_Interp *UNUSED(interp),
	struct ged *gedp,
	void *draw_view_ctx,
	int argc,
	const char *argv[],
	const char *usage)
{
    void *gdvp = draw_view_ctx;
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    to_refresh_suppress_all_begin(current_top);

    int ret = to_data_move_point_mode_func(gedp, gdvp, argc, argv, usage);
    to_refresh_suppress_all_end(current_top);
    return ret;
}


static int
to_data_move_point_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* shift the command name to argv[1] before calling to_data_move_point_mode_func */
    argv[1] = argv[0];
    return to_data_move_point_mode_func(gedp, gdvp, argc-1, argv+1, usage);
}


static int
to_data_move_point_mode_func(struct ged *gedp,
	void *gdvp,
	int UNUSED(argc),
	const char *argv[],
	const char *usage)
{
    int x, y;

    ged_view_active_ctx_set(gedp, gdvp);

    if (bu_sscanf(argv[1], "%d", &x) != 1 ||
	    bu_sscanf(argv[2], "%d", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    /* At the moment, only the TclCAD polygon mode is being used. */
    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_DATA_MOVE_POINT_MODE);

    return BRLCAD_OK;
}


int
go_data_pick(struct ged *gedp,
	void *draw_view_ctx,
	int argc,
	const char *argv[],
	const char *usage)
{
    void *gdvp = draw_view_ctx;
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 2 || 3 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    to_refresh_suppress_all_begin(current_top);

    int ret = to_data_pick_func(gedp, gdvp, argc, argv, usage);
    to_refresh_suppress_all_end(current_top);
    return ret;
}


static int
to_data_pick(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 3 || 4 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* shift the command name to argv[1] before calling to_data_pick_func */
    argv[1] = argv[0];
    return to_data_pick_func(gedp, gdvp, argc-1, argv+1, usage);
}


static int
to_data_pick_func(struct ged *gedp,
	void *gdvp,
	int argc,
	const char *argv[],
	const char *usage)
{
    int mx, my, width, height;
    fastf_t cx, cy;
    fastf_t vx, vy;
    fastf_t sf;
    point_t dpoint, vpoint;
    mat_t model2view;
    int i;
    fastf_t top_z = -MAX_FASTF;
    point_t top_point = VINIT_ZERO;
    size_t top_i = 0;
    size_t top_j = 0;
    size_t top_k = 0;
    int found_top = 0;
    const char *top_data_str = NULL;
    const char *top_data_label = NULL;
    struct bu_vls top_label_vls = BU_VLS_INIT_ZERO;
    static fastf_t tol = 0.015;
    static const char *data_polygons_str = "data_polygons";
    static const char *data_labels_str = "data_labels";
    static const char *sdata_labels_str = "sdata_labels";
    static const char *sdata_lines_str = "sdata_lines";
    static const char *data_arrows_str = "data_arrows";
    static const char *sdata_arrows_str = "sdata_arrows";
    static const char *data_axes_str = "data_axes";
    static const char *sdata_axes_str = "sdata_axes";

    if (argc == 2) {
	if (bu_sscanf(argv[1], "%d %d", &mx, &my) != 2)
	    goto bad;
    } else {
	if (bu_sscanf(argv[1], "%d", &mx) != 1)
	    goto bad;

	if (bu_sscanf(argv[2], "%d", &my) != 1)
	    goto bad;
    }

    width = tclcad_commands_width(gdvp);
    cx = 0.5 * (fastf_t)width;
    height = tclcad_commands_height(gdvp);
    cy = 0.5 * (fastf_t)height;
    sf = 2.0 / width;
    vx = (mx - cx) * sf;
    vy = (cy - my) * sf;
    bv_model2view_get(model2view, tclcad_commands_bv_const(gdvp));

    /* check for polygon points */
    tclcad_polygon_state *data_polygons = tclcad_view_polygon_state_from_view_ctx(gdvp, 0);
    if (data_polygons && data_polygons->gdps_draw &&
	    data_polygons->gdps_polygons.num_polygons) {
	size_t si, sj, sk;

	tclcad_polygon_state *gdpsp = data_polygons;

	for (si = 0; si < gdpsp->gdps_polygons.num_polygons; ++si)
	    for (sj = 0; sj < gdpsp->gdps_polygons.polygon[si].num_contours; ++sj)
		for (sk = 0; sk < gdpsp->gdps_polygons.polygon[si].contour[sj].num_points; ++sk) {
		    fastf_t minX, maxX;
		    fastf_t minY, maxY;

		    MAT4X3PNT(vpoint, model2view, gdpsp->gdps_polygons.polygon[si].contour[sj].point[sk]);
		    minX = vpoint[X] - tol;
		    maxX = vpoint[X] + tol;
		    minY = vpoint[Y] - tol;
		    maxY = vpoint[Y] + tol;

		    if (minX < vx && vx < maxX &&
			    minY < vy && vy < maxY) {
			if (!found_top || top_z < vpoint[Z]) {
			    top_z = vpoint[Z];
			    top_data_str = data_polygons_str;
			    top_i = si;
			    top_j = sj;
			    top_k = sk;
			    VMOVE(top_point, gdpsp->gdps_polygons.polygon[si].contour[sj].point[sk]);
			    found_top = 1;
			}
		    }
		}
    }

    if (found_top) {
	bu_vls_printf(gedp->ged_result_str, "%s {%zu %zu %zu} {%lf %lf %lf}",
		top_data_str, top_i, top_j, top_k, V3ARGS(top_point));
	bu_vls_free(&top_label_vls);
	return BRLCAD_OK;
    }

    /* check for label points */
    {
	const char *label_name = "_tcl_data_labels";
	size_t _child_cnt = _tclcad_draw_view_data_labels_count(gdvp, label_name);
	if (_child_cnt > 0) {
	    for (size_t _k = 0; _k < _child_cnt; _k++) {
		fastf_t minX, maxX;
		fastf_t minY, maxY;

		struct bu_vls label = BU_VLS_INIT_ZERO;
		if (!_tclcad_draw_view_data_label_copy(gdvp, label_name, _k, &label, dpoint, NULL)) {
		    bu_vls_free(&label);
		    continue;
		}

		MAT4X3PNT(vpoint, model2view, dpoint);

		minX = vpoint[X];
		maxX = vpoint[X] + (2 * tol);
		minY = vpoint[Y];
		maxY = vpoint[Y] + (2 * tol);

		if (minX < vx && vx < maxX &&
			minY < vy && vy < maxY) {
		    if (!found_top || top_z < vpoint[Z]) {
			top_z = vpoint[Z];
			top_data_str = data_labels_str;
			top_i = _k;
			bu_vls_sprintf(&top_label_vls, "%s", bu_vls_cstr(&label));
			top_data_label = bu_vls_cstr(&top_label_vls);
			VMOVE(top_point, dpoint);
			found_top = 1;
		    }
		}
		bu_vls_free(&label);
	    }
	}
    }

    /* check for selected label points */
    {
	const char *label_name = "_tcl_sdata_labels";
	size_t _child_cnt = _tclcad_draw_view_data_labels_count(gdvp, label_name);
	if (_child_cnt > 0) {
	    for (size_t _k = 0; _k < _child_cnt; _k++) {
		fastf_t minX, maxX;
		fastf_t minY, maxY;

		struct bu_vls label = BU_VLS_INIT_ZERO;
		if (!_tclcad_draw_view_data_label_copy(gdvp, label_name, _k, &label, dpoint, NULL)) {
		    bu_vls_free(&label);
		    continue;
		}

		MAT4X3PNT(vpoint, model2view, dpoint);

		minX = vpoint[X];
		maxX = vpoint[X] + (2 * tol);
		minY = vpoint[Y];
		maxY = vpoint[Y] + (2 * tol);

		if (minX < vx && vx < maxX &&
			minY < vy && vy < maxY) {
		    if (!found_top || top_z < vpoint[Z]) {
			top_z = vpoint[Z];
			top_data_str = sdata_labels_str;
			top_i = _k;
			bu_vls_sprintf(&top_label_vls, "%s", bu_vls_cstr(&label));
			top_data_label = bu_vls_cstr(&top_label_vls);
			VMOVE(top_point, dpoint);
			found_top = 1;
		    }
		}
		bu_vls_free(&label);
	    }
	}
    }

    if (found_top) {
	bu_vls_printf(gedp->ged_result_str, "%s %zu {{%s} {%lf %lf %lf}}",
		top_data_str, top_i, top_data_label, V3ARGS(top_point));
	bu_vls_free(&top_label_vls);
	return BRLCAD_OK;
    }

    /* check for line points */
    {
	const char *feature_name = "_tcl_data_lines";
	point_t *_lpts = NULL;
	int _lnpts = _tclcad_draw_view_data_lines_points_copy(gdvp, feature_name, &_lpts);
	for (i = 0; i < _lnpts; ++i) {
	    fastf_t minX, maxX;
	    fastf_t minY, maxY;
	    VMOVE(dpoint, _lpts[i]);
	    MAT4X3PNT(vpoint, model2view, dpoint);
	    minX = vpoint[X] - tol; maxX = vpoint[X] + tol;
	    minY = vpoint[Y] - tol; maxY = vpoint[Y] + tol;
	    if (minX < vx && vx < maxX && minY < vy && vy < maxY) {
		bu_vls_printf(gedp->ged_result_str, "data_lines %d {%lf %lf %lf}", i, V3ARGS(dpoint));
		bu_free(_lpts, "TclCAD draw-view points");
		bu_vls_free(&top_label_vls);
		return BRLCAD_OK;
	    }
	}
	if (_lpts)
	    bu_free(_lpts, "TclCAD draw-view points");
    }

    /* check for selected line points */
    {
	const char *feature_name = "_tcl_sdata_lines";
	point_t *_lpts = NULL;
	int _lnpts = _tclcad_draw_view_data_lines_points_copy(gdvp, feature_name, &_lpts);
	for (i = 0; i < _lnpts; ++i) {
	    fastf_t minX, maxX;
	    fastf_t minY, maxY;
	    VMOVE(dpoint, _lpts[i]);
	    MAT4X3PNT(vpoint, model2view, dpoint);
	    minX = vpoint[X] - tol; maxX = vpoint[X] + tol;
	    minY = vpoint[Y] - tol; maxY = vpoint[Y] + tol;
	    if (minX < vx && vx < maxX && minY < vy && vy < maxY) {
		if (!found_top || top_z < vpoint[Z]) {
		    top_z = vpoint[Z];
		    top_data_str = sdata_lines_str;
		    top_i = (size_t)i;
		    VMOVE(top_point, dpoint);
		    found_top = 1;
		}
	    }
	}
	if (_lpts)
	    bu_free(_lpts, "TclCAD draw-view points");
    }

    if (found_top) {
	bu_vls_printf(gedp->ged_result_str, "%s %zu {%lf %lf %lf}",
		top_data_str, top_i, V3ARGS(top_point));
	bu_vls_free(&top_label_vls);
	return BRLCAD_OK;
    }

    /* check for arrow points */
    {
	const char *feature_name = "_tcl_data_arrows";
	point_t *_apts = NULL;
	int _anpts = _tclcad_draw_view_data_arrows_points_copy(gdvp, feature_name, &_apts);
	for (i = 0; i < _anpts; ++i) {
	    fastf_t minX, maxX, minY, maxY;
	    VMOVE(dpoint, _apts[i]);
	    MAT4X3PNT(vpoint, model2view, dpoint);
	    minX = vpoint[X] - tol; maxX = vpoint[X] + tol;
	    minY = vpoint[Y] - tol; maxY = vpoint[Y] + tol;
	    if (minX < vx && vx < maxX && minY < vy && vy < maxY) {
		if (!found_top || top_z < vpoint[Z]) {
		    top_z = vpoint[Z];
		    top_data_str = data_arrows_str;
		    top_i = (size_t)i;
		    VMOVE(top_point, dpoint);
		    found_top = 1;
		}
	    }
	}
	if (_apts)
	    bu_free(_apts, "TclCAD draw-view points");
    }

    /* check for selected arrow points */
    {
	const char *feature_name = "_tcl_sdata_arrows";
	point_t *_apts = NULL;
	int _anpts = _tclcad_draw_view_data_arrows_points_copy(gdvp, feature_name, &_apts);
	for (i = 0; i < _anpts; ++i) {
	    fastf_t minX, maxX, minY, maxY;
	    VMOVE(dpoint, _apts[i]);
	    MAT4X3PNT(vpoint, model2view, dpoint);
	    minX = vpoint[X] - tol; maxX = vpoint[X] + tol;
	    minY = vpoint[Y] - tol; maxY = vpoint[Y] + tol;
	    if (minX < vx && vx < maxX && minY < vy && vy < maxY) {
		if (!found_top || top_z < vpoint[Z]) {
		    top_z = vpoint[Z];
		    top_data_str = sdata_arrows_str;
		    top_i = (size_t)i;
		    VMOVE(top_point, dpoint);
		    found_top = 1;
		}
	    }
	}
	if (_apts)
	    bu_free(_apts, "TclCAD draw-view points");
    }

    if (found_top) {
	bu_vls_printf(gedp->ged_result_str, "%s %zu {%lf %lf %lf}",
		top_data_str, top_i, V3ARGS(top_point));
	bu_vls_free(&top_label_vls);
	return BRLCAD_OK;
    }

    /* check for axes points */
    {
	const char *feature_name = "_tcl_data_axes";
	point_t *_cpts = NULL;
	int _ncpts = _tclcad_draw_view_data_axes_centers_copy(gdvp, feature_name, &_cpts);
	for (i = 0; i < _ncpts; ++i) {
	    fastf_t minX, maxX, minY, maxY;
	    VMOVE(dpoint, _cpts[i]);
	    MAT4X3PNT(vpoint, model2view, dpoint);
	    minX = vpoint[X] - tol; maxX = vpoint[X] + tol;
	    minY = vpoint[Y] - tol; maxY = vpoint[Y] + tol;
	    if (minX < vx && vx < maxX && minY < vy && vy < maxY) {
		if (!found_top || top_z < vpoint[Z]) {
		    top_z = vpoint[Z];
		    top_i = (size_t)i;
		    top_data_str = data_axes_str;
		    VMOVE(top_point, dpoint);
		    found_top = 1;
		}
	    }
	}
	if (_cpts)
	    bu_free(_cpts, "TclCAD draw-view axes points");
    }

    /* check for selected axes points */
    {
	const char *feature_name = "_tcl_sdata_axes";
	point_t *_cpts = NULL;
	int _ncpts = _tclcad_draw_view_data_axes_centers_copy(gdvp, feature_name, &_cpts);
	for (i = 0; i < _ncpts; ++i) {
	    fastf_t minX, maxX, minY, maxY;
	    VMOVE(dpoint, _cpts[i]);
	    MAT4X3PNT(vpoint, model2view, dpoint);
	    minX = vpoint[X] - tol; maxX = vpoint[X] + tol;
	    minY = vpoint[Y] - tol; maxY = vpoint[Y] + tol;
	    if (minX < vx && vx < maxX && minY < vy && vy < maxY) {
		if (!found_top || top_z < vpoint[Z]) {
		    top_z = vpoint[Z];
		    top_i = (size_t)i;
		    top_data_str = sdata_axes_str;
		    VMOVE(top_point, dpoint);
		    found_top = 1;
		}
	    }
	}
	if (_cpts)
	    bu_free(_cpts, "TclCAD draw-view axes points");
    }

    if (found_top)
	bu_vls_printf(gedp->ged_result_str, "%s %zu {%lf %lf %lf}",
		top_data_str, top_i, V3ARGS(top_point));

    bu_vls_free(&top_label_vls);
    return BRLCAD_OK;

bad:
    bu_vls_free(&top_label_vls);
    bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
    return BRLCAD_ERROR;
}


static int
to_data_vZ(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    /* must be double for scanf */
    double vZ;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 2 && argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* Get the data vZ */
    if (argc == 2) {
	bu_vls_printf(gedp->ged_result_str, "%lf", tclcad_view_data_vZ_from_view_ctx(gdvp));
	return BRLCAD_OK;
    }

    /* Set the data vZ */
    if (bu_sscanf(argv[2], "%lf", &vZ) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    (void)tclcad_view_data_vZ_set(gdvp, vZ);

    return BRLCAD_OK;
}


static void
to_init_default_bindings(void *gdvp)
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;
    const int endpoint_keyboard_input =
	tclcad_commands_endpoint_input_enabled(gdvp);

    if (tclcad_view_pathname_vls(gdvp)) {
	struct bu_vls *pathvls = tclcad_view_pathname_vls(gdvp);
	if (pathvls) {
	    bu_vls_printf(&bindings, "bind %s <Configure> {%s configure %s; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Enter> {focus %s; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(pathvls));
	    bu_vls_printf(&bindings, "bind %s <Expose> {%s handle_expose %s %%c; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "catch {wm protocol %s WM_DELETE_WINDOW {%s delete_view %s; break}}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));

	    /* Mouse Bindings */
	    bu_vls_printf(&bindings, "bind %s <2> {%s vslew %s %%x %%y; focus %s; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp),
		    bu_vls_addr(pathvls));
	    bu_vls_printf(&bindings, "bind %s <1> {%s zoom %s 0.5; focus %s; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp),
		    bu_vls_addr(pathvls));
	    bu_vls_printf(&bindings, "bind %s <3> {%s zoom %s 2.0; focus %s;  break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp),
		    bu_vls_addr(pathvls));
	    bu_vls_printf(&bindings, "bind %s <4> {%s zoom %s 1.1; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <5> {%s zoom %s 0.9; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <MouseWheel> {if {%%D < 0} {%s zoom %s 0.9} else {%s zoom %s 1.1}; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));

	    /* Idle Mode */
	    bu_vls_printf(&bindings, "bind %s <ButtonRelease> {%s idle_mode %s}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <KeyRelease-Control_L> {%s idle_mode %s}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <KeyRelease-Control_R> {%s idle_mode %s}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <KeyRelease-Shift_L> {%s idle_mode %s}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <KeyRelease-Shift_R> {%s idle_mode %s}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <KeyRelease-Alt_L> {%s idle_mode %s; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <KeyRelease-Alt_R> {%s idle_mode %s; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));

	    /* Rotate Mode */
	    bu_vls_printf(&bindings, "bind %s <Control-ButtonRelease-1> {%s idle_mode %s}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Control-ButtonPress-1> {%s rotate_mode %s %%x %%y}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Control-ButtonPress-2> {%s rotate_mode %s %%x %%y}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Control-ButtonPress-3> {%s rotate_mode %s %%x %%y}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));

	    /* Translate Mode */
	    bu_vls_printf(&bindings, "bind %s <Shift-ButtonRelease-1> {%s idle_mode %s}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Shift-ButtonPress-1> {%s translate_mode %s %%x %%y}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Shift-ButtonPress-2> {%s translate_mode %s %%x %%y}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Shift-ButtonPress-3> {%s translate_mode %s %%x %%y}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));

	    /* Scale Mode */
	    bu_vls_printf(&bindings, "bind %s <Control-Shift-ButtonRelease-1> {%s idle_mode %s}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Control-Shift-ButtonPress-1> {%s scale_mode %s %%x %%y}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Control-Shift-ButtonPress-2> {%s scale_mode %s %%x %%y}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Control-Shift-ButtonPress-3> {%s scale_mode %s %%x %%y}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));

	    /* Constrained Rotate Mode */
	    bu_vls_printf(&bindings, "bind %s <Control-Lock-ButtonRelease-1> {%s idle_mode %s}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Control-Lock-ButtonPress-1> {%s constrain_rmode %s x %%x %%y; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Control-Lock-ButtonPress-2> {%s constrain_rmode %s y %%x %%y; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Control-Lock-ButtonPress-3> {%s constrain_rmode %s z %%x %%y; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));

	    /* Constrained Translate Mode */
	    bu_vls_printf(&bindings, "bind %s <Shift-Lock-ButtonRelease-1> {%s idle_mode %s; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Shift-Lock-ButtonPress-1> {%s constrain_tmode %s x %%x %%y; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Shift-Lock-ButtonPress-2> {%s constrain_tmode %s y %%x %%y; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Shift-Lock-ButtonPress-3> {%s constrain_tmode %s z %%x %%y; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));

	    /* The endpoint owns its keyboard-only profile.  Keep pointer and
	     * application-mode bindings in TclCAD until they have equivalent
	     * semantic actions, but do not execute shared view shortcuts twice. */
	    if (!endpoint_keyboard_input) {
	    /* Key Bindings */
	    bu_vls_printf(&bindings, "bind %s 3 {%s aet %s 35 25; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s 4 {%s aet %s 45 45; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s f {%s aet %s 0 0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s F {%s aet %s 0 0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s R {%s aet %s 180 0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s r {%s aet %s 270 0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s l {%s aet %s 90 0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s L {%s aet %s 90 0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s t {%s aet %s 270 90; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s T {%s aet %s 270 90; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s b {%s aet %s 270 -90; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s B {%s aet %s 270 -90; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    }
	    bu_vls_printf(&bindings, "bind %s + {%s zoom %s 2.0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s = {%s zoom %s 2.0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s _ {%s zoom %s 0.5; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s - {%s zoom %s 0.5; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Key-Left> {%s rot %s -v 0 1 0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Key-Right> {%s rot %s -v 0 -1 0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Key-Up> {%s rot %s -v 1 0 0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));
	    bu_vls_printf(&bindings, "bind %s <Key-Down> {%s rot %s -v -1 0 0; break}; ",
		    bu_vls_addr(pathvls),
		    bu_vls_addr(&current_top->to_gedp->go_name),
		    tclcad_commands_view_name(gdvp));

	    Tcl_Eval(current_top->to_interp, bu_vls_addr(&bindings));
	}
    }
    bu_vls_free(&bindings);
}


static int
to_dplot(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    int ac;
    int ret;
    char *av[256];
    struct bu_vls callback_cmd = BU_VLS_INIT_ZERO;
    struct bu_vls temp = BU_VLS_INIT_ZERO;
    struct bu_vls result_copy = BU_VLS_INIT_ZERO;
    void *gdvp;
    struct tclcad_ged_data *tgd = (struct tclcad_ged_data *)current_top->to_gedp->u_data;
    const char *who_av[3] = {"who", "b", NULL};
    int first = 1;
    int aflag = 0;

    /* copy all args */
    ac = argc;
    for (int i = 0; i < ac; ++i)
	av[i] = bu_strdup((char *)argv[i]);
    av[ac] = (char *)0;

    /* check for displayed objects */
    ret = ged_exec_who(gedp, 2, (const char **)who_av);
    if (ret == BRLCAD_OK && strlen(bu_vls_addr(gedp->ged_result_str)) == 0)
	aflag = 1;
    bu_vls_trunc(gedp->ged_result_str, 0);

    while ((*func)(gedp, ac, (const char **)av) & GED_MORE) {
	int ac_more;
	const char **avmp;
	const char **av_more = NULL;

	/* save result string */
	bu_vls_substr(&result_copy, gedp->ged_result_str, 0, bu_vls_strlen(gedp->ged_result_str));
	bu_vls_trunc(gedp->ged_result_str, 0);

	ret = ged_exec_who(gedp, 1, (const char **)who_av);
	if (ret == BRLCAD_OK && strlen(bu_vls_addr(gedp->ged_result_str)) == 0)
	    aflag = 1;

	bu_vls_trunc(gedp->ged_result_str, 0);

	struct bu_ptbl *views = ged_view_set_views_ctx(current_top->to_gedp);
	for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	    gdvp = BU_PTBL_GET(views, i);
	    if (to_is_viewable(gdvp)) {
		tclcad_commands_sync_dimensions(ged_view_active_ctx(gedp), gdvp);
	    }
	}

	if (first && aflag) {
	    first = 0;
	    to_autoview_all_views(current_top);
	} else {
	    to_refresh_all_views(current_top);
	}
	/* restore result string */
	bu_vls_substr(gedp->ged_result_str, &result_copy, 0, bu_vls_strlen(&result_copy));
	bu_vls_free(&result_copy);

	if (0 < bu_vls_strlen(&tgd->go_more_args_callback)) {
	    bu_vls_trunc(&callback_cmd, 0);
	    bu_vls_printf(&callback_cmd, "%s [string range {%s} 0 end]",
		    bu_vls_addr(&tgd->go_more_args_callback),
		    bu_vls_addr(gedp->ged_result_str));

	    if (Tcl_Eval(current_top->to_interp, bu_vls_addr(&callback_cmd)) != TCL_OK) {
		bu_vls_trunc(gedp->ged_result_str, 0);
		bu_vls_printf(gedp->ged_result_str, "%s", Tcl_GetStringResult(current_top->to_interp));
		Tcl_ResetResult(current_top->to_interp);
		return BRLCAD_ERROR;
	    }

	    bu_vls_trunc(&temp, 0);
	    bu_vls_printf(&temp, "%s", Tcl_GetStringResult(current_top->to_interp));
	    Tcl_ResetResult(current_top->to_interp);
	} else {
	    bu_log("\r%s", bu_vls_addr(gedp->ged_result_str));
	    bu_vls_trunc(&temp, 0);
	    if (bu_vls_gets(&temp, stdin) < 0) {
		break;
	    }
	}

	if (Tcl_SplitList(current_top->to_interp, bu_vls_addr(&temp), &ac_more, &av_more) != TCL_OK) {
	    continue;
	}

	if (ac_more < 1) {
	    /* space has still been allocated */
	    Tcl_Free((char *)av_more);

	    continue;
	}

	/* skip first element if empty */
	avmp = av_more;
	if (*avmp[0] == '\0') {
	    --ac_more;
	    ++avmp;
	}

	/* ignore last element if empty */
	if (*avmp[ac_more-1] == '\0')
	    --ac_more;

	/* copy additional args */
	for (int i = 0; i < ac_more; ++i)
	    av[ac++] = bu_strdup(avmp[i]);
	av[ac+1] = (char *)0;

	Tcl_Free((char *)av_more);
    }

    struct bu_ptbl *views = ged_view_set_views_ctx(current_top->to_gedp);
    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	gdvp = BU_PTBL_GET(views, i);
	if (to_is_viewable(gdvp)) {
	    tclcad_commands_sync_dimensions(ged_view_active_ctx(gedp), gdvp);
	}
    }
    to_refresh_all_views(current_top);

    bu_vls_free(&callback_cmd);
    bu_vls_free(&temp);

    bu_vls_printf(gedp->ged_result_str, "BUILT_BY_MORE_ARGS");
    for (int i = 0; i < ac; ++i) {
	bu_vls_printf(gedp->ged_result_str, "%s ", av[i]);
	bu_free((void *)av[i], "to_more_args_func");
    }
    return BRLCAD_OK;
}

static int
to_fontsize(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int fontsize;
    struct bv *view;
    struct bv_params_state params = BV_PARAMS_STATE_INIT;
    struct bv_other_state center_dot = BV_OTHER_STATE_INIT;
    struct bv_other_state scale_overlay = BV_OTHER_STATE_INIT;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 2 || 3 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }
    view = tclcad_commands_bv(gdvp);
    if (!view || !bv_params_state_get(&params, view) ||
	!bv_center_dot_state_get(&center_dot, view) ||
	!bv_scale_overlay_state_get(&scale_overlay, view))
	return BRLCAD_ERROR;

    /* get the font size */
    if (argc == 2) {
	bu_vls_printf(gedp->ged_result_str, "%d", params.font_size);
	return BRLCAD_OK;
    }

    /* set background color */
    if (bu_sscanf(argv[2], "%d", &fontsize) != 1)
	goto bad_fontsize;

    if (fontsize == 0)
	fontsize = 20;
    if (fontsize < 5 || fontsize > 96)
	goto bad_fontsize;

    bobol_display_endpoint_t *endpoint = tclcad_commands_endpoint(gdvp);
    if (endpoint) {
	const char *font_properties[] = {
	    "view.faceplate.params.font_size",
	    "view.faceplate.center_dot.font_size",
	    "view.faceplate.scale.font_size"
	};
	struct bobol_endpoint_property_value property =
	    BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	property.type = BOBOL_ENDPOINT_PROPERTY_UINT;
	property.uint_value = (uint64_t)fontsize;
	for (size_t i = 0; i < sizeof(font_properties) / sizeof(font_properties[0]); i++) {
	    if (bobol_display_endpoint_property_set(endpoint,
		font_properties[i], &property) != BOBOL_ENDPOINT_PROPERTY_OK)
		return BRLCAD_ERROR;
	}
    } else {
	/* An unattached view has no renderer policy surface yet. */
	params.font_size = fontsize;
	center_dot.gos_font_size = fontsize;
	scale_overlay.gos_font_size = fontsize;
	if (!bv_params_state_set(view, &params) ||
	    !bv_center_dot_state_set(view, &center_dot) ||
	    !bv_scale_overlay_state_set(view, &scale_overlay))
	    return BRLCAD_ERROR;
    }
    to_refresh_view(gdvp);
    return BRLCAD_OK;

bad_fontsize:
    bu_vls_printf(gedp->ged_result_str, "%s: %s", argv[0], argv[2]);
    return BRLCAD_ERROR;
}


static int
to_fit_png_image(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    icv_image_t *img;
    size_t o_w_requested, o_n_requested;
    fastf_t sf;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6 ||
	    bu_sscanf(argv[2], "%zu", &o_w_requested) != 1 ||
	    bu_sscanf(argv[3], "%zu", &o_n_requested) != 1 ||
	    bu_sscanf(argv[4], "%lf", &sf) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    img = icv_read(argv[1], BU_MIME_IMAGE_PNG, 0, 0);
    if (img == NULL) {
	bu_vls_printf(gedp->ged_result_str, "%s: %s does not exist or is not a png file", argv[0], argv[1]);
	return BRLCAD_ERROR;
    }

    int ret = icv_fit(img, gedp->ged_result_str, o_w_requested, o_n_requested, sf);
    if (ret == BRLCAD_ERROR) {
	icv_destroy(img);
	return ret;
    }

    /* icv_write should return < 0 for errors but doesn't */
    if (icv_write(img, argv[5], BU_MIME_IMAGE_PNG) == 0) {
	icv_destroy(img);
	return BRLCAD_ERROR;
    }

    icv_destroy(img);
    return BRLCAD_OK;
}


static int
to_init_view_bindings(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    to_init_default_bindings(gdvp);

    return BRLCAD_OK;
}


static int
to_delete_view(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (!current_top || current_top->to_gedp != gedp) {
	bu_vls_printf(gedp->ged_result_str, "View host is unavailable - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    const int was_active = (ged_view_active_ctx(gedp) == gdvp);
    tclcad_view_host_destroy(current_top, gdvp);
    (void)ged_view_set_context_remove(ged_view_set_ctx(gedp), gdvp);

    if (was_active) {
	ged_view_active_ctx_set(gedp, NULL);
	struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
	if (views && BU_PTBL_LEN(views))
	    ged_view_active_ctx_set(gedp, BU_PTBL_GET(views, 0));
    }

    /* The fbserv bridge is session-owned.  Keep its single capture provider
     * on the active surviving endpoint after any pane is removed. */
    void *active_view_ctx = ged_view_active_ctx(gedp);
    if (active_view_ctx)
	(void)to_open_fbs(active_view_ctx, current_top->to_interp);

    ged_view_context_free(gdvp);
    return BRLCAD_OK;
}


static int
to_handle_expose(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int count;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3 ||
	    bu_sscanf(argv[2], "%d", &count) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s, argv[2] - %s", argv[0], usage, argv[2]);
	return BRLCAD_ERROR;
    }

    /* There are more expose events to come so ignore this one */
    if (count)
	return BRLCAD_OK;

    return to_handle_refresh(gedp, argv[1]);
}

// TODO - does this do anything?  It sscanfs the value into hide_view,
// but then doesn't do anything with it...
static int
to_hide_view(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int hide_view;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc > 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* return the hide view setting */
    if (argc == 2) {
	bu_vls_printf(gedp->ged_result_str, "%d", tclcad_view_hide_from_view_ctx(gdvp));
	return BRLCAD_OK;
    }

    if (bu_sscanf(argv[2], "%d", &hide_view) != 1) {
	bu_vls_printf(gedp->ged_result_str, "%s: bad value - %s, should be 0 or 1", argv[1], argv[2]);
	return BRLCAD_ERROR;
    }

    return BRLCAD_OK;
}


struct redraw_edited_path_data {
    struct ged *gedp;
    void *gdvp;
    int *need_refresh;
};


static void
redraw_edited_paths(struct bu_hash_tbl *t, void *udata)
{
    const char *av[5] = {0};
    uint8_t *key;
    char *draw_path;
    struct redraw_edited_path_data *data;
    int ret, draw_mode = 0;
    struct bu_vls path_dmode = BU_VLS_INIT_ZERO;
    struct tclcad_path_edit_params *params;
    struct bu_hash_entry *entry = bu_hash_next(t, NULL);

    while (entry) {

	bu_hash_key(entry, &key, NULL);
	draw_path = (char *)key;

	data = (struct redraw_edited_path_data *)udata;

	params = (struct tclcad_path_edit_params *)bu_hash_value(entry, NULL);
	if (params->edit_mode == TCLCAD_OTRANSLATE_MODE) {
	    struct bu_vls tcl_cmd = BU_VLS_INIT_ZERO;
	    struct bu_vls tran_x_vls = BU_VLS_INIT_ZERO;
	    struct bu_vls tran_y_vls = BU_VLS_INIT_ZERO;
	    struct bu_vls tran_z_vls = BU_VLS_INIT_ZERO;
	    mat_t dvec;

	    MAT_DELTAS_GET(dvec, params->edit_mat);
	    VSCALE(dvec, dvec, data->gedp->dbip->dbi_base2local);

	    bu_vls_printf(&tran_x_vls, "%lf", dvec[X]);
	    bu_vls_printf(&tran_y_vls, "%lf", dvec[Y]);
	    bu_vls_printf(&tran_z_vls, "%lf", dvec[Z]);
	    MAT_IDN(params->edit_mat);

	    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(data->gdvp);
	    bu_vls_printf(&tcl_cmd, "%s otranslate %s %s %s",
		    bu_vls_addr(&tvd->gdv_edit_motion_delta_callback),
		    bu_vls_addr(&tran_x_vls), bu_vls_addr(&tran_y_vls),
		    bu_vls_addr(&tran_z_vls));
	    tvd->gdv_edit_motion_delta_callback_cnt++;
	    if (tvd->gdv_edit_motion_delta_callback_cnt > 1) {
		bu_log("Warning - recursive gdv_edit_motion_delta_callback call\n");
	    }

	    Tcl_Eval(current_top->to_interp, bu_vls_addr(&tcl_cmd));
	    tvd->gdv_edit_motion_delta_callback_cnt++;
	    bu_vls_free(&tcl_cmd);
	    bu_vls_free(&tran_x_vls);
	    bu_vls_free(&tran_y_vls);
	    bu_vls_free(&tran_z_vls);
	}

	av[0] = "how";
	av[1] = draw_path;
	av[2] = NULL;
	ret = ged_exec_how(data->gedp, 2, av);
	if (ret == BRLCAD_OK) {
	    bu_sscanf(bu_vls_cstr(data->gedp->ged_result_str), "%d", &draw_mode);
	}
	if (draw_mode == 5) {
	    bu_vls_printf(&path_dmode, "-h");
	} else {
	    bu_vls_printf(&path_dmode, "-m%d", draw_mode);
	}

	av[0] = "erase";
	ret = ged_exec_erase(data->gedp, 2, av);

	if (ret == BRLCAD_OK) {
	    av[0] = "draw";
	    av[1] = "-R";
	    av[2] = bu_vls_cstr(&path_dmode);
	    av[3] = draw_path;
	    av[4] = NULL;
	    ged_exec_draw(data->gedp, 4, av);
	}

	*data->need_refresh = 1;

	entry = bu_hash_next(t, entry);
    }
}


static int
to_idle_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int mode, need_refresh = 0;
    struct redraw_edited_path_data data;
    struct tclcad_ged_data *tgd = (struct tclcad_ged_data *)current_top->to_gedp->u_data;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    mode = tclcad_view_polygon_mode_from_view_ctx(gdvp);

    ged_draw_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
    if (ged_draw_view_context_lod_policy_get(&lod_policy, gdvp) &&
	    lod_policy.csg_enabled &&
	    lod_policy.zoom_refresh &&
	    mode == TCLCAD_SCALE_MODE)
    {
	const char *av[] = {"redraw", NULL};

	ged_exec_redraw(gedp, 1, (const char **)av);

	need_refresh = 1;
    }

    if (mode != TCLCAD_POLY_CONTOUR_MODE ||
	    tclcad_view_polygon_cflag_from_view_ctx(gdvp, 0) == 0)
    {
	struct bu_vls bindings = BU_VLS_INIT_ZERO;

	struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
	if (pathname && bu_vls_strlen(pathname)) {
	    bu_vls_printf(&bindings, "bind %s <Motion> {}", bu_vls_cstr(pathname));
	    Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
	}
	bu_vls_free(&bindings);
    }

    struct bv_grid_state grid;
    (void)bv_grid_state_get(&grid, tclcad_commands_bv_const(gdvp));
    if (grid.snap &&
	    (mode == TCLCAD_TRANSLATE_MODE ||
	     mode == TCLCAD_CONSTRAINED_TRANSLATE_MODE))
    {
	const char *av[3];

	tclcad_commands_sync_dimensions(gdvp, gdvp);

	ged_view_active_ctx_set(gedp, gdvp);
	av[0] = "grid";
	av[1] = "vsnap";
	av[2] = NULL;
	ged_exec_grid(gedp, 2, (const char **)av);

	struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
	if (0 < bu_vls_strlen(&tvd->gdv_callback)) {
	    tvd->gdv_callback_cnt++;
	    if (tvd->gdv_callback_cnt > 1) {
		bu_log("Warning - recursive gvd_callback call\n");
	    }
	    tclcad_eval_noresult(current_top->to_interp, bu_vls_addr(&tvd->gdv_callback), 0, NULL);
	    tvd->gdv_callback_cnt--;
	}

	need_refresh = 1;
    }

    /* redraw any edited paths, then clear them from our table */
    Tcl_Eval(current_top->to_interp, "SetWaitCursor $::ArcherCore::application");
    data.gedp = gedp;
    data.gdvp = gdvp;
    data.need_refresh = &need_refresh;

    redraw_edited_paths(tgd->go_edited_paths, &data);

    free_path_edit_params(tgd->go_edited_paths);
    bu_hash_destroy(tgd->go_edited_paths);
    tgd->go_edited_paths = bu_hash_create(0);
    Tcl_Eval(current_top->to_interp, "SetNormalCursor $::ArcherCore::application");

    if (need_refresh) {
	to_refresh_all_views(current_top);
    }

    tclcad_view_polygon_mode_set(gdvp, TCLCAD_IDLE_MODE);
    tclcad_view_polygon_cflag_clear(gdvp, 1);

    return BRLCAD_OK;
}


static int
to_light(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int light;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (3 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* get light flag */
    if (argc == 2) {
	if (!tclcad_commands_endpoint_bool_get(gdvp,
		"renderer.lighting", &light))
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str, "%d", light);
	return BRLCAD_OK;
    }

    /* set light flag */
    if (bu_sscanf(argv[2], "%d", &light) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    if (light < 0)
	light = 0;

    if (!tclcad_commands_endpoint_bool_set(gdvp, "renderer.lighting",
	    light))
	return BRLCAD_ERROR;
    to_refresh_view(gdvp);

    return BRLCAD_OK;
}


static int
to_list_views(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    void *gdvp;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s", argv[0]);
	return BRLCAD_ERROR;
    }

    struct bu_ptbl *views = ged_view_set_views_ctx(current_top->to_gedp);
    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	gdvp = BU_PTBL_GET(views, i);
	bu_vls_printf(gedp->ged_result_str, "%s ", tclcad_commands_view_name(gdvp));
    }

    return BRLCAD_OK;
}


static int
to_local2base(struct ged *gedp,
	int UNUSED(argc),
	const char *UNUSED(argv[]),
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    bu_vls_printf(gedp->ged_result_str, "%lf", current_top->to_gedp->dbip->dbi_local2base);

    return BRLCAD_OK;
}


static int
to_lod(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    void *gdvp;

    struct bu_ptbl *views = ged_view_set_views_ctx(current_top->to_gedp);
    for (size_t i = 0; i < BU_PTBL_LEN(views); i++) {
	gdvp = BU_PTBL_GET(views, i);
	ged_view_active_ctx_set(gedp, gdvp);
	(*func)(gedp, argc, (const char **)argv);
    }

    return BRLCAD_OK;
}


static int
to_make(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    int ret;
    const char *av[3];

    ret = ged_exec(gedp, argc, argv);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[argc-2];
	av[2] = (char *)0;
	to_autoview_func(gedp, 2, (const char **)av, ged_exec, (char *)0, TO_UNLIMITED);
    }

    return ret;
}


static int
to_mirror(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    int ret;
    const char *av[3];

    ret = ged_exec(gedp, argc, argv);

    if (ret == BRLCAD_OK) {
	av[0] = "draw";
	av[1] = (char *)argv[argc-1];
	av[2] = (char *)0;
	to_autoview_func(gedp, 2, (const char **)av, ged_exec, (char *)0, TO_UNLIMITED);
    }

    return ret;
}


static int
to_edit_motion_delta_callback(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int i;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* get the callback string */
    if (argc == 2) {
	struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
	bu_vls_printf(gedp->ged_result_str, "%s", bu_vls_addr(&tvd->gdv_edit_motion_delta_callback));

	return BRLCAD_OK;
    }

    /* set the callback string */
    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
    bu_vls_trunc(&tvd->gdv_edit_motion_delta_callback, 0);
    for (i = 2; i < argc; ++i)
	bu_vls_printf(&tvd->gdv_edit_motion_delta_callback, "%s ", argv[i]);

    return BRLCAD_OK;
}


static int
to_more_args_callback(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    int i;
    struct tclcad_ged_data *tgd = (struct tclcad_ged_data *)current_top->to_gedp->u_data;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* get the callback string */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "%s", bu_vls_addr(&tgd->go_more_args_callback));

	return BRLCAD_OK;
    }

    /* set the callback string */
    bu_vls_trunc(&tgd->go_more_args_callback, 0);
    for (i = 1; i < argc; ++i)
	bu_vls_printf(&tgd->go_more_args_callback, "%s ", argv[i]);

    return BRLCAD_OK;
}


static int
to_view_cmd(ClientData UNUSED(clientData),
	Tcl_Interp *UNUSED(interp),
	int UNUSED(argc),
	char **UNUSED(argv))
{
    return TCL_OK;
}


static int
to_move_arb_edge_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	    bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_MOVE_ARB_EDGE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_move_arb_edge %s %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2],
		argv[3]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_move_arb_face_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	    bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_MOVE_ARB_FACE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_move_arb_face %s %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2],
		argv[3]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_bot_move_pnt(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    int ret;

    if ((ret = ged_exec(gedp, argc, argv)) == BRLCAD_OK) {
	const char *av[3];
	int i;

	if (argc == 4)
	    i = 1;
	else
	    i = 2;

	av[0] = "draw";
	av[1] = (char *)argv[i];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);

	return BRLCAD_OK;
    }

    return ret;
}


static int
to_bot_move_pnts(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    int ret;

    if ((ret = ged_exec(gedp, argc, argv)) == BRLCAD_OK) {
	const char *av[3];

	av[0] = "draw";
	av[1] = (char *)argv[1];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);

	return BRLCAD_OK;
    }

    return ret;
}


static int
to_bot_move_pnt_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	    bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_MOVE_BOT_POINT_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_bot_move_pnt -r %s %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2],
		argv[3]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_bot_move_pnts_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int i;
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &x) != 1 ||
	    bu_sscanf(argv[3], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_MOVE_BOT_POINTS_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_bot_move_pnts %s %%x %%y %s ",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[4]);
    }
    for (i = 5; i < argc; ++i)
	bu_vls_printf(&bindings, "%s ", argv[i]);
    bu_vls_printf(&bindings, "}");

    Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_metaball_move_pnt_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	    bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_MOVE_METABALL_POINT_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_metaball_move_pnt %s %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2],
		argv[3]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_pipe_move_pnt_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	    bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_MOVE_PIPE_POINT_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_pipe_move_pnt %s %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2],
		argv[3]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_move_pnt_common(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr func,
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    int ret;

    if ((ret = (*func)(gedp, argc, argv)) == BRLCAD_OK) {
	const char *av[3];
	int i;

	if (argc == 4)
	    i = 1;
	else
	    i = 2;

	av[0] = "draw";
	av[1] = (char *)argv[i];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);

	return BRLCAD_OK;
    }

    return ret;
}

static int
tclcad_view_host_open(void *view_ctx, struct tclcad_view_data *tvd,
	Tcl_Interp *interp, int argc, const char *argv[], struct bu_vls *result)
{
    if (!view_ctx || !tvd || !interp || argc < 3 || !argv)
	return 0;

    const char *name = argv[1];
    if (name[0] == '.')
	bu_vls_strcpy(&tvd->gdv_pathname, name);
    else
	bu_vls_printf(&tvd->gdv_pathname, ".%s", name);

    if (BU_STR_EQUAL(argv[2], "nu"))
	return 1;
    if (!BU_STR_EQUAL(argv[2], "tkobol"))
	return 0;

#ifndef HAVE_TKOBOL_HOST
    bu_vls_printf(result, "Tk Obol host support is unavailable\n");
    return 0;
#else
    int width = bv_width_get(tclcad_commands_bv(view_ctx));
    int height = bv_height_get(tclcad_commands_bv(view_ctx));
    int toplevel = 1;
    int software = 0;
    for (int i = 3; i < argc; i++) {
	if (BU_STR_EQUAL(argv[i], "sw")) {
	    software = 1;
	    continue;
	}
	if (BU_STR_EQUAL(argv[i], "hw")) {
	    software = 0;
	    continue;
	}
	if (i + 1 >= argc)
	    continue;
	if (BU_STR_EQUAL(argv[i], "-W"))
	    width = atoi(argv[++i]);
	else if (BU_STR_EQUAL(argv[i], "-N"))
	    height = atoi(argv[++i]);
	else if (BU_STR_EQUAL(argv[i], "-S") ||
		BU_STR_EQUAL(argv[i], "-s"))
	    width = height = atoi(argv[++i]);
	else if (BU_STR_EQUAL(argv[i], "-t"))
	    toplevel = atoi(argv[++i]) ? 1 : 0;
    }
    if (width <= 0)
	width = 512;
    if (height <= 0)
	height = 512;
    (void)bv_dimensions_set(tclcad_commands_bv(view_ctx), width, height);

    if (!tclcad_obol_host_factories_register()) {
	bu_vls_printf(result, "Tk Obol host factory registration failed\n");
	return 0;
    }
    bobol_display_endpoint_t *endpoint =
	bobol_display_endpoint_create(NULL, 0);
    if (!endpoint) {
	bu_vls_printf(result, "Obol display endpoint creation failed\n");
	return 0;
    }
    if (!bobol_display_endpoint_render_engine_set(endpoint,
	    software ? BOBOL_RENDER_ENGINE_SW : BOBOL_RENDER_ENGINE_HW) ||
	!ged_view_context_display_endpoint_set(view_ctx, endpoint, 1)) {
	bobol_display_endpoint_destroy(endpoint);
	bu_vls_printf(result, "Obol view endpoint attachment failed\n");
	return 0;
    }

    struct bobol_host_desc desc = {0};
    desc.struct_size = sizeof(desc);
    desc.mode = toplevel ? BOBOL_HOST_MODE_TOPLEVEL :
	BOBOL_HOST_MODE_EMBEDDED;
    desc.width = (unsigned int)width;
    desc.height = (unsigned int)height;
    desc.device_pixel_ratio = 1.0;
    desc.visible = 1;
    desc.required_capabilities = software ?
	BOBOL_HOST_CAP_PIXEL_PRESENT : BOBOL_HOST_CAP_SYSTEM_GL;
    desc.title = bu_vls_cstr(&tvd->gdv_pathname);
    desc.native_id_hint = bu_vls_cstr(&tvd->gdv_pathname);
    desc.application_context = interp;
    if (!bobol_display_endpoint_host_open(endpoint,
	    software ? "tk-photo" : "tk-gl", &desc)) {
	(void)ged_view_context_display_endpoint_set(view_ctx, NULL, 0);
	bu_vls_printf(result, "%s host open failed\n",
	    software ? "TkPhoto" : "Tk OpenGL");
	return 0;
    }
    if (bobol_display_endpoint_host_capabilities(endpoint) &
	BOBOL_HOST_CAP_INPUT) {
	if (!bobol_display_endpoint_input_profile_set(endpoint,
		bobol_input_keyboard_view_profile()) ||
	    !bobol_display_endpoint_input_action_handler_set(endpoint,
		tclcad_obol_input_action, view_ctx)) {
	    (void)ged_view_context_display_endpoint_set(view_ctx, NULL, 0);
	    bu_vls_printf(result, "Tk Obol input endpoint setup failed\n");
	    return 0;
	}
    }
    return 1;
#endif
}


static int
to_new_view(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    void *view_set_ctx = NULL;
    void *active_view_ctx = NULL;
    void *new_view_ctx = NULL;
    int reuse_active_view = 0;
    static const int name_index = 1;
    const char *type = NULL;
    struct bu_vls event_vls = BU_VLS_INIT_ZERO;

    GED_CHECK_DATABASE_OPEN(gedp, BRLCAD_ERROR);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    if (BU_STR_EQUAL(argv[2], "nu") || BU_STR_EQUAL(argv[2], "tkobol"))
	type = argv[2];

    if (!type) {
	bu_vls_printf(gedp->ged_result_str,
		"Unsupported view presentation type '%s'; expected tkobol or nu.\n",
		argv[2]);
	return BRLCAD_ERROR;
    }

    view_set_ctx = ged_view_set_ctx(current_top->to_gedp);
    active_view_ctx = ged_view_active_ctx(gedp);
    if (active_view_ctx && !tclcad_view_data_from_view_ctx(active_view_ctx)) {
	new_view_ctx = active_view_ctx;
	reuse_active_view = 1;
    } else {
	new_view_ctx = ged_view_context_create_with_set(view_set_ctx);
	if (!new_view_ctx) {
	    bu_vls_printf(gedp->ged_result_str, "Failed to create %s\n", argv[1]);
	    return BRLCAD_ERROR;
	}
	ged_view_set_context_add(view_set_ctx, new_view_ctx);
	ged_view_context_owned_add(gedp, new_view_ctx);
    }

    struct bu_ptbl *callbacks;
    BU_GET(callbacks, struct bu_ptbl);
    bu_ptbl_init(callbacks, 8, "bv callbacks");

    struct tclcad_view_data *tvd;
    BU_GET(tvd, struct tclcad_view_data);
    tclcad_view_data_init(tvd, current_top->to_gedp);

    if (!bv_name_set(tclcad_commands_bv(new_view_ctx), argv[name_index]) ||
	    !tclcad_view_data_bind_view_ctx(new_view_ctx, tvd) ||
	    !ged_view_context_callbacks_set(new_view_ctx, callbacks)) {
	tclcad_view_data_unbind_view_ctx(new_view_ctx);
	bu_ptbl_free(callbacks);
	BU_PUT(callbacks, struct bu_ptbl);
	BU_PUT(tvd, struct tclcad_view_data);
	if (!reuse_active_view)
	    ged_view_context_free(new_view_ctx);

	bu_vls_printf(gedp->ged_result_str, "Failed to initialize %s\n", argv[1]);
	return BRLCAD_ERROR;
    }
    callbacks = NULL;

    if (!tclcad_view_host_open(new_view_ctx, tvd, current_top->to_interp,
	    argc, argv, gedp->ged_result_str)) {
	tclcad_view_data_unbind_view_ctx(new_view_ctx);
	bu_vls_free(&tvd->gdv_pathname);
	bu_vls_free(&tvd->gdv_edit_motion_delta_callback);
	bu_vls_free(&tvd->gdv_callback);
	BU_PUT(tvd, struct tclcad_view_data);
	if (!reuse_active_view)
	    ged_view_context_free(new_view_ctx);
	bu_vls_printf(gedp->ged_result_str, "Failed to create %s\n", argv[1]);
	return BRLCAD_ERROR;
    }

    ged_draw_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
    if (ged_draw_view_context_lod_policy_get(&lod_policy, new_view_ctx)) {
	lod_policy.point_scale = 1.0;
	lod_policy.curve_scale = 1.0;
	ged_draw_view_context_lod_policy_apply(new_view_ctx, &lod_policy);
    }

    /* Set default bindings */
    to_init_default_bindings(new_view_ctx);

    struct bu_vls *pathname = tclcad_view_pathname_vls(new_view_ctx);
    const char *view_name = bv_name_get(tclcad_commands_bv_const(new_view_ctx));
    if (!view_name)
	view_name = argv[name_index];
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&event_vls, "event generate %s <Configure>; %s autoview %s",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		view_name);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&event_vls));
    }
    bu_vls_free(&event_vls);

    if (pathname && bu_vls_strlen(pathname)) {
	(void)Tcl_CreateCommand(current_top->to_interp,
		bu_vls_cstr(pathname),
		(Tcl_CmdProc *)to_view_cmd,
		(ClientData)new_view_ctx,
		NULL);
    }

    // If we don't already have a default GED view, set the one we just
    // created as the new default
    if (!ged_view_active_ctx(gedp))
	ged_view_active_ctx_set(gedp, new_view_ctx);

    /* A new auxiliary pane must not retarget the session stream. */
    if (ged_view_active_ctx(gedp) == new_view_ctx)
	(void)to_open_fbs(new_view_ctx, current_top->to_interp);

    bu_vls_printf(gedp->ged_result_str, "%s", view_name);
    return BRLCAD_OK;
}


static int
to_orotate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	    bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_OROTATE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_orotate %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_oscale_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	    bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_OSCALE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_oscale %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


/* to_model_delta_mode */
static int
to_otranslate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 5) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[3], "%lf", &x) != 1 ||
	    bu_sscanf(argv[4], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_OTRANSLATE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_otranslate %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_paint_rect_area(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }
    if (argc < 2 || 7 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    struct bv_interactive_rect_state rect = BV_INTERACTIVE_RECT_STATE_INIT;
    (void)bv_interactive_rect_state_get(&rect, tclcad_commands_bv_const(gdvp));
    (void)rect;
    (void)ged_draw_obol_framebuffer_present(gedp);
    to_refresh_view(gdvp);

    return BRLCAD_OK;
}


static int
tclcad_commands_capture_rgb(struct ged *gedp, void *view_ctx,
	unsigned char **pixels, unsigned int *width, unsigned int *height)
{
    if (!gedp || !view_ctx || !pixels || !width || !height)
	return 0;
    *pixels = NULL;
    *width = 0;
    *height = 0;

    bobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(view_ctx);
    if (!endpoint)
	return 0;

    (void)ged_draw_obol_framebuffer_present(gedp);
    if (!bobol_display_endpoint_view_sync(endpoint, view_ctx))
	return 0;

    unsigned char *source = NULL;
    size_t source_size = 0;
    unsigned int components = 0;
    if (!bobol_display_endpoint_capture(endpoint, &source, &source_size,
	    width, height, &components) || !source || !*width || !*height ||
	    (components != 3 && components != 4) ||
	    source_size < (size_t)(*width) * (*height) * components) {
	if (source)
	    bu_free(source, "TclCAD endpoint capture");
	return 0;
    }

    const size_t pixel_count = (size_t)(*width) * (*height);
    unsigned char *rgb = (unsigned char *)bu_malloc(pixel_count * 3,
	"TclCAD RGB endpoint capture");
    for (size_t i = 0; i < pixel_count; i++) {
	rgb[i * 3] = source[i * components];
	rgb[i * 3 + 1] = source[i * components + 1];
	rgb[i * 3 + 2] = source[i * components + 2];
    }
    bu_free(source, "TclCAD endpoint capture");
    *pixels = rgb;
    return 1;
}


static int
to_pix(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    FILE *fp = NULL;
    unsigned char *pixels = NULL;
    unsigned int width = 0;
    unsigned int height = 0;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if ((fp = fopen(argv[2], "wb")) == NULL) {
	bu_vls_printf(gedp->ged_result_str,
		"%s: cannot open \"%s\" for writing.",
		argv[0], argv[2]);
	return BRLCAD_ERROR;
    }

    if (!tclcad_commands_capture_rgb(gedp, gdvp, &pixels, &width, &height)) {
	bu_vls_printf(gedp->ged_result_str,
		"%s: view endpoint did not return an RGB image\n", argv[0]);
	fclose(fp);
	return BRLCAD_ERROR;
    }

    const size_t image_size = (size_t)width * height * 3;
    if (fwrite(pixels, image_size, 1, fp) != 1) {
	bu_vls_printf(gedp->ged_result_str,
		"%s: failed writing RGB image to %s\n", argv[0], argv[2]);
	bu_free(pixels, "pixels");
	fclose(fp);
	return BRLCAD_ERROR;
    }

    bu_free(pixels, "pixels");
    fclose(fp);

    return BRLCAD_OK;
}


static int
to_png(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    png_structp png_p;
    png_infop info_p;
    FILE *fp = NULL;
    unsigned char **rows = NULL;
    unsigned char *pixels = NULL;
    static int bits_per_channel = 8;  /* bits per color channel */
    unsigned int width = 0;
    unsigned int height = 0;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if ((fp = fopen(argv[2], "wb")) == NULL) {
	bu_vls_printf(gedp->ged_result_str,
		"%s: cannot open \"%s\" for writing.",
		argv[0], argv[2]);
	return BRLCAD_ERROR;
    }

    png_p = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_p) {
	bu_vls_printf(gedp->ged_result_str, "%s: could not create PNG write structure.", argv[0]);
	fclose(fp);
	return BRLCAD_ERROR;
    }

    info_p = png_create_info_struct(png_p);
    if (!info_p) {
	bu_vls_printf(gedp->ged_result_str, "%s: could not create PNG info structure.", argv[0]);
	png_destroy_write_struct(&png_p, NULL);
	fclose(fp);
	return BRLCAD_ERROR;
    }

    if (!tclcad_commands_capture_rgb(gedp, gdvp, &pixels, &width, &height)) {
	bu_vls_printf(gedp->ged_result_str,
		"%s: view endpoint did not return an RGB image\n", argv[0]);
	png_destroy_write_struct(&png_p, &info_p);
	fclose(fp);
	return BRLCAD_ERROR;
    }

    rows = (unsigned char **)bu_calloc(height, sizeof(unsigned char *), "rows");

    const size_t bytes_per_line = (size_t)width * 3;
    for (unsigned int i = 0; i < height; ++i)
	rows[i] = pixels + ((size_t)(height-i-1) * bytes_per_line);

    png_init_io(png_p, fp);
    png_set_filter(png_p, 0, PNG_FILTER_NONE);
    png_set_compression_level(png_p, 9);
    png_set_IHDR(png_p, info_p, width, height, bits_per_channel,
	    PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	    PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_gAMA(png_p, info_p, 0.5);
    png_write_info(png_p, info_p);
    png_write_image(png_p, rows);
    png_write_end(png_p, NULL);

    bu_free(rows, "rows");
    bu_free(pixels, "pixels");
    png_destroy_write_struct(&png_p, &info_p);
    fclose(fp);

    return BRLCAD_OK;
}


static int
to_rect_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int ac;
    int x, y;
    const char *av[5];
    struct bu_vls bindings = BU_VLS_INIT_ZERO;
    struct bu_vls x_vls = BU_VLS_INIT_ZERO;
    struct bu_vls y_vls = BU_VLS_INIT_ZERO;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    ged_view_active_ctx_set(gedp, gdvp);

    if (bu_sscanf(argv[2], "%d", &x) != 1 ||
	    bu_sscanf(argv[3], "%d", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x,
	    tclcad_commands_height(gdvp) - y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_RECTANGLE_MODE);

    ac = 4;
    av[0] = "rect";
    av[1] = "dim";
    av[2] = "0";
    av[3] = "0";
    av[4] = (char *)0;
    (void)ged_exec_rect(gedp, ac, (const char **)av);

    fastf_t prev_x = 0.0;
    fastf_t prev_y = 0.0;
    (void)bv_previous_mouse_get(&prev_x, &prev_y,
	    tclcad_commands_bv_const(gdvp));
    bu_vls_printf(&x_vls, "%d", (int)prev_x);
    bu_vls_printf(&y_vls, "%d", (int)prev_y);
    av[1] = "pos";
    av[2] = bu_vls_addr(&x_vls);
    av[3] = bu_vls_addr(&y_vls);
    (void)ged_exec_rect(gedp, ac, (const char **)av);
    bu_vls_free(&x_vls);
    bu_vls_free(&y_vls);

    ac = 3;
    av[0] = "rect";
    av[1] = "draw";
    av[2] = "1";
    av[3] = (char *)0;
    (void)ged_exec_rect(gedp, ac, (const char **)av);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_rect %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp));
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    to_refresh_view(gdvp);

    return BRLCAD_OK;
}


static int
to_rotate_arb_face_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 7) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[5], "%lf", &x) != 1 ||
	    bu_sscanf(argv[6], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_ROTATE_ARB_FACE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_rotate_arb_face %s %s %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2],
		argv[3],
		argv[4]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_rotate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &x) != 1 ||
	    bu_sscanf(argv[3], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_ROTATE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_rot %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp));
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


/**
 * Called when the named proc created by rt_gettrees() is destroyed.
 */
static void
to_deleteProc_rt(ClientData clientData)
{
    struct application *ap = (struct application *)clientData;
    struct rt_i *rtip;

    RT_AP_CHECK(ap);
    rtip = ap->a_rt_i;
    RT_CK_RTI(rtip);

    rt_i_destroy(rtip);
    ap->a_rt_i = (struct rt_i *)NULL;

    bu_free((void *)ap, "struct application");
}


static int
to_rt_end_callback(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    int i;
    struct tclcad_ged_data *tgd = (struct tclcad_ged_data *)current_top->to_gedp->u_data;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* get the callback string */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "%s", bu_vls_addr(&tgd->go_rt_end_callback));

	return BRLCAD_OK;
    }

    /* set the callback string */
    bu_vls_trunc(&tgd->go_rt_end_callback, 0);
    for (i = 1; i < argc; ++i)
	bu_vls_printf(&tgd->go_rt_end_callback, "%s ", argv[i]);

    return BRLCAD_OK;
}


/**
 * @brief
 * Given an instance of a database and the name of some treetops,
 * create a named "ray-tracing" object (proc) which will respond to
 * subsequent operations.  Returns new proc name as result.
 *
 * @par Example:
 *	.inmem rt_gettrees .rt all.g light.r
 */
int
to_rt_gettrees(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct application *ap;
    char *newprocname;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc < 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    newprocname = (char *)argv[1];

    /* Delete previous proc (if any) to release all that memory, first */
    (void)Tcl_DeleteCommand(current_top->to_interp, newprocname);

    /* Skip past newprocname when calling to_rt_gettrees_application */
    if ((ap=to_rt_gettrees_application(gedp, argc-2, argv+2)) == RT_APPLICATION_NULL) {
	return BRLCAD_ERROR;
    }

    /* Instantiate the proc, with clientData of wdb */
    /* Beware, returns a "token", not TCL_OK. */
    (void)Tcl_CreateCommand(current_top->to_interp, newprocname, tclcad_rt,
	    (ClientData)ap, to_deleteProc_rt);

    /* Return new function name as result */
    bu_vls_printf(gedp->ged_result_str, "%s", newprocname);

    return BRLCAD_OK;
}


static int
to_protate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	    bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_PROTATE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_protate %s %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2],
		argv[3]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_pscale_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	    bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_PSCALE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_pscale %s %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2],
		argv[3]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_ptranslate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 6) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[4], "%lf", &x) != 1 ||
	    bu_sscanf(argv[5], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_PTRANSLATE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_ptranslate %s %s %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp),
		argv[2],
		argv[3]);
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_data_scale_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &x) != 1 ||
	    bu_sscanf(argv[3], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_DATA_SCALE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_data_scale %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp));
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_scale_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &x) != 1 ||
	    bu_sscanf(argv[3], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_SCALE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_scale %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp));
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_screen2model(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    point_t view;
    point_t model;
    mat_t view2model;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &x) != 1 ||
	    bu_sscanf(argv[3], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    tclcad_commands_sync_dimensions(gdvp, gdvp);
    bv_screen_to_view(&x, &y, tclcad_commands_bv_const(gdvp), x, y);
    VSET(view, x, y, 0.0);
    bv_view2model_get(view2model, tclcad_commands_bv_const(gdvp));
    MAT4X3PNT(model, view2model, view);

    bu_vls_printf(gedp->ged_result_str, "%lf %lf %lf", V3ARGS(model));

    return BRLCAD_OK;
}


static int
to_screen2view(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    point_t view;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &x) != 1 ||
	    bu_sscanf(argv[3], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    tclcad_commands_sync_dimensions(gdvp, gdvp);
    bv_screen_to_view(&x, &y, tclcad_commands_bv_const(gdvp), x, y);
    VSET(view, x, y, 0.0);

    bu_vls_printf(gedp->ged_result_str, "%lf %lf %lf", V3ARGS(view));

    return BRLCAD_OK;
}


static int
to_set_coord(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (3 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* Get coord */
    if (argc == 2) {
	bu_vls_printf(gedp->ged_result_str, "%c",
		bv_coord_get(tclcad_commands_bv_const(gdvp)));
	return BRLCAD_OK;
    }

    /* Set coord */
    if ((argv[2][0] != 'm' && argv[2][0] != 'v') || argv[2][1] != '\0') {
	bu_vls_printf(gedp->ged_result_str, "set_coord: bad value - %s\n", argv[2]);
	return BRLCAD_ERROR;
    }

    bv_coord_set(tclcad_commands_bv(gdvp), argv[2][0]);

    return BRLCAD_OK;
}


static int
to_snap_view(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    fastf_t fvx, fvy;

    /* must be double for scanf */
    double vx, vy;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &vx) != 1 ||
	    bu_sscanf(argv[3], "%lf", &vy) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }
    /* convert double to fastf_t */
    fvx = vx;
    fvy = vy;

    tclcad_commands_sync_dimensions(gdvp, gdvp);
    bv_unit_conversion_set(tclcad_commands_bv(gdvp),
	    bv_local2base_get(tclcad_commands_bv_const(gdvp)),
	    gedp->dbip->dbi_base2local);

    ged_view_active_ctx_set(gedp, gdvp);
    struct bv_grid_state grid;
    (void)bv_grid_state_get(&grid, tclcad_commands_bv_const(gdvp));
    if (!bv_snap_lines_get(tclcad_commands_bv_const(gdvp)) && !grid.snap) {
	bu_vls_printf(gedp->ged_result_str, "%lf %lf", fvx, fvy);
	return BRLCAD_OK;
    }

    {
	unsigned long long snap_kinds = tclcad_commands_prepare_snap(gdvp);
	if (snap_kinds)
	    tclcad_commands_snap_point_2d(gdvp, &fvx, &fvy, snap_kinds);
    }

    bu_vls_printf(gedp->ged_result_str, "%lf %lf", fvx, fvy);

    return BRLCAD_OK;
}


static int
to_bot_edge_split(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    int ret;

    if ((ret = ged_exec(gedp, argc, argv)) == BRLCAD_OK) {
	const char *av[3];
	struct bu_vls save_result;

	bu_vls_init(&save_result);

	av[0] = "draw";
	av[1] = (char *)argv[1];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);

	bu_vls_trunc(gedp->ged_result_str, 0);
	bu_vls_printf(gedp->ged_result_str, "%s", bu_vls_addr(&save_result));
	bu_vls_free(&save_result);

	return BRLCAD_OK;
    }

    return ret;
}


static int
to_bot_face_split(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *UNUSED(usage),
	int UNUSED(maxargs))
{
    int ret;

    if ((ret = ged_exec(gedp, argc, argv)) == BRLCAD_OK) {
	const char *av[3];
	struct bu_vls save_result;

	bu_vls_init(&save_result);

	av[0] = "draw";
	av[1] = (char *)argv[1];
	av[2] = (char *)0;
	to_edit_redraw(gedp, 2, (const char **)av);

	bu_vls_trunc(gedp->ged_result_str, 0);
	bu_vls_printf(gedp->ged_result_str, "%s", bu_vls_addr(&save_result));
	bu_vls_free(&save_result);

	return BRLCAD_OK;
    }

    return ret;
}


static int
to_translate_mode(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    struct bu_vls bindings = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double x, y;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &x) != 1 ||
	    bu_sscanf(argv[3], "%lf", &y) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    bv_previous_mouse_set(tclcad_commands_bv(gdvp), x, y);
    tclcad_view_polygon_mode_set(gdvp, TCLCAD_TRANSLATE_MODE);

    struct bu_vls *pathname = tclcad_view_pathname_vls(gdvp);
    if (pathname && bu_vls_strlen(pathname)) {
	bu_vls_printf(&bindings, "bind %s <Motion> {%s mouse_trans %s %%x %%y}",
		bu_vls_cstr(pathname),
		bu_vls_cstr(&current_top->to_gedp->go_name),
		tclcad_commands_view_name(gdvp));
	Tcl_Eval(current_top->to_interp, bu_vls_cstr(&bindings));
    }
    bu_vls_free(&bindings);

    return BRLCAD_OK;
}


static int
to_transparency(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int transparency;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 2 && argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* get transparency flag */
    if (argc == 2) {
	if (!tclcad_commands_endpoint_bool_get(gdvp,
		"renderer.transparency", &transparency))
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str, "%d", transparency);
	return BRLCAD_OK;
    }

    /* set transparency flag */
    if (argc == 3) {
	if (bu_sscanf(argv[2], "%d", &transparency) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "%s: invalid transparency value - %s", argv[0], argv[2]);
	    return BRLCAD_ERROR;
	}

	if (!tclcad_commands_endpoint_bool_set(gdvp,
		"renderer.transparency", transparency))
	    return BRLCAD_ERROR;
	to_refresh_view(gdvp);
	return BRLCAD_OK;
    }

    return BRLCAD_OK;
}


static int
to_view_callback(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int i;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* get the callback string */
    if (argc == 2) {
	struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
	bu_vls_printf(gedp->ged_result_str, "%s", bu_vls_addr(&tvd->gdv_callback));

	return BRLCAD_OK;
    }

    /* set the callback string */
    struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
    bu_vls_trunc(&tvd->gdv_callback, 0);
    for (i = 2; i < argc; ++i)
	bu_vls_printf(&tvd->gdv_callback, "%s ", argv[i]);

    return BRLCAD_OK;
}


static int
to_view_win_size(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int width, height;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc < 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc > 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (argc == 2) {
	bu_vls_printf(gedp->ged_result_str, "%d %d",
	    tclcad_commands_width(gdvp), tclcad_commands_height(gdvp));
	return BRLCAD_OK;
    }

    if (argc == 3) {
	if (bu_sscanf(argv[2], "%d", &width) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "%s: bad size %s", argv[0], argv[2]);
	    return BRLCAD_ERROR;
	}

	height = width;
    } else {
	if (bu_sscanf(argv[2], "%d", &width) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "%s: bad width %s", argv[0], argv[2]);
	    return BRLCAD_ERROR;
	}

	if (bu_sscanf(argv[3], "%d", &height) != 1) {
	    bu_vls_printf(gedp->ged_result_str, "%s: bad height %s", argv[0], argv[3]);
	    return BRLCAD_ERROR;
	}
    }

    if (width <= 0 || height <= 0 ||
	!bobol_display_endpoint_resize(tclcad_commands_endpoint(gdvp),
	    (unsigned int)width, (unsigned int)height, 1.0))
	return BRLCAD_ERROR;
    bv_context_dimensions_set((struct bv_context *)gdvp, width, height);

    return BRLCAD_OK;
}


static int
to_view2screen(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int width, height;
    fastf_t x, y;
    fastf_t aspect;

    /* must be double for scanf */
    double view[ELEMENTS_PER_POINT];

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf %lf", &view[X], &view[Y]) != 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    width = tclcad_commands_width(gdvp);
    height = tclcad_commands_height(gdvp);
    aspect = (fastf_t)width/(fastf_t)height;
    x = (view[X] + 1.0) * 0.5 * (fastf_t)width;
    y = (view[Y] * aspect - 1.0) * -0.5 * (fastf_t)height;

    bu_vls_printf(gedp->ged_result_str, "%d %d", (int)x, (int)y);

    return BRLCAD_OK;
}


static int
to_vmake(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    {
	int ret;
	const char *av[8];
	char center[512];
	char scale[128];
	mat_t view_center;
	fastf_t view_scale;

	bv_center_mat_get(view_center, tclcad_commands_bv_const(gdvp));
	view_scale = bv_scale_get(tclcad_commands_bv_const(gdvp));
	sprintf(center, "%f %f %f",
		-view_center[MDX],
		-view_center[MDY],
		-view_center[MDZ]);
	sprintf(scale, "%f", view_scale * 2.0);

	av[0] = (char *)argv[0];
	av[1] = "-o";
	av[2] = center;
	av[3] = "-s";
	av[4] = scale;
	av[5] = (char *)argv[2];
	av[6] = (char *)argv[3];
	av[7] = (char *)0;

	ret = ged_exec(gedp, 7, (const char **)av);

	if (ret == BRLCAD_OK) {
	    av[0] = "draw";
	    av[1] = (char *)argv[2];
	    av[2] = (char *)0;
	    to_autoview_func(gedp, 2, (const char **)av, ged_exec, (char *)0, TO_UNLIMITED);
	}

	return ret;
    }
}


static int
to_vslew(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int ret, width, height;
    int ac;
    const char *av[3];
    fastf_t xpos2, ypos2;
    fastf_t sf;
    struct bu_vls slew_vec = BU_VLS_INIT_ZERO;

    /* must be double for scanf */
    double xpos1, ypos1;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 4) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    if (bu_sscanf(argv[2], "%lf", &xpos1) != 1 ||
	    bu_sscanf(argv[3], "%lf", &ypos1) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    width = tclcad_commands_width(gdvp);
    xpos2 = 0.5 * (fastf_t)width;
    height = tclcad_commands_height(gdvp);
    ypos2 = 0.5 * (fastf_t)height;
    sf = 2.0 / width;

    bu_vls_printf(&slew_vec, "%lf %lf", (xpos1 - xpos2) * sf, (ypos2 - ypos1) * sf);

    ged_view_active_ctx_set(gedp, gdvp);
    ac = 2;
    av[0] = (char *)argv[0];
    av[1] = bu_vls_addr(&slew_vec);
    av[2] = (char *)0;

    ret = ged_exec(gedp, ac, (const char **)av);
    bu_vls_free(&slew_vec);

    if (ret == BRLCAD_OK) {
	struct bv_grid_state grid;
	(void)bv_grid_state_get(&grid, tclcad_commands_bv_const(gdvp));
	if (grid.snap) {

	    tclcad_commands_sync_dimensions(gdvp, gdvp);

	    ged_view_active_ctx_set(gedp, gdvp);
	    av[0] = "grid";
	    av[1] = "vsnap";
	    av[2] = NULL;
	    ged_exec_grid(gedp, 2, (const char **)av);
	}

	/* Legacy line-snap centering is retired; Obol scene snapping must own any
	 * future line-snap recenter behavior. */

	struct tclcad_view_data *tvd = tclcad_view_data_from_view_ctx(gdvp);
	if (0 < bu_vls_strlen(&tvd->gdv_callback)) {
	    tvd->gdv_callback_cnt++;
	    if (tvd->gdv_callback_cnt > 1) {
		bu_log("Warning - recursive gvd_callback call\n");
	    }
	    Tcl_Eval(current_top->to_interp, bu_vls_addr(&tvd->gdv_callback));
	    tvd->gdv_callback_cnt--;
	}

	to_refresh_view(gdvp);
    }

    return ret;
}


static int
to_zbuffer(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int zbuffer;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (3 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* get zbuffer flag */
    if (argc == 2) {
	if (!tclcad_commands_endpoint_bool_get(gdvp,
		"renderer.depth_test", &zbuffer))
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str, "%d", zbuffer);
	return BRLCAD_OK;
    }

    /* set zbuffer flag */
    if (bu_sscanf(argv[2], "%d", &zbuffer) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    if (zbuffer < 0)
	zbuffer = 0;
    else if (1 < zbuffer)
	zbuffer = 1;

    if (!tclcad_commands_endpoint_bool_set(gdvp, "renderer.depth_test",
	    zbuffer))
	return BRLCAD_ERROR;
    to_refresh_view(gdvp);

    return BRLCAD_OK;
}


static int
to_zclip(struct ged *gedp,
	int argc,
	const char *argv[],
	ged_func_ptr UNUSED(func),
	const char *usage,
	int UNUSED(maxargs))
{
    int zclip;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (3 < argc) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    void *gdvp = ged_view_find_ctx(gedp, argv[1]);
    if (!gdvp) {
	bu_vls_printf(gedp->ged_result_str, "View not found - %s", argv[1]);
	return BRLCAD_ERROR;
    }

    /* get zclip flag */
    if (argc == 2) {
	bu_vls_printf(gedp->ged_result_str, "%d",
		bv_zclip_get(tclcad_commands_bv_const(gdvp)));
	return BRLCAD_OK;
    }

    /* set zclip flag */
    if (bu_sscanf(argv[2], "%d", &zclip) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    if (zclip < 0)
	zclip = 0;
    else if (1 < zclip)
	zclip = 1;

    if (!tclcad_commands_endpoint_bool_set(gdvp, "view.zclip", zclip))
	return BRLCAD_ERROR;
    to_refresh_view(gdvp);

    return BRLCAD_OK;
}


/*************************** Local Utility Functions ***************************/

static void
to_rt_end_callback_internal(int aborted)
{
    struct tclcad_ged_data *tgd = (struct tclcad_ged_data *)current_top->to_gedp->u_data;
    if (0 < bu_vls_strlen(&tgd->go_rt_end_callback)) {
	struct bu_vls callback_cmd = BU_VLS_INIT_ZERO;
	tgd->go_rt_end_callback_cnt++;
	if (tgd->go_rt_end_callback_cnt > 1) {
	    bu_log("Warning - recursive go_rt_end_callback call\n");
	}
	bu_vls_printf(&callback_cmd, "%s %d",
		bu_vls_addr(&tgd->go_rt_end_callback), aborted);
	Tcl_Eval(current_top->to_interp, bu_vls_addr(&callback_cmd));
	tgd->go_rt_end_callback_cnt--;
    }
}


static void
to_output_handler(struct ged *gedp, char *line)
{
    const char *script;

    if (gedp->ged_output_script != (char *)0)
	script = gedp->ged_output_script;
    else
	script = "puts";

    tclcad_eval_noresult(current_top->to_interp, script, 1, (const char **)&line);
}


int
go_run_tclscript(Tcl_Interp *interp,
	const char *tclscript,
	struct bu_vls *result_str)
{
    int ret;

    /* initialize result */
    bu_vls_trunc(result_str, 0);

    ret = Tcl_Eval(interp, tclscript);
    bu_vls_printf(result_str, "%s", Tcl_GetStringResult(interp));
    Tcl_ResetResult(interp);

    return ret;
}


struct application *
to_rt_gettrees_application(struct ged *gedp,
	int argc,
	const char *argv[])
{
    struct rt_i *rtip;
    struct application *ap;
    static struct resource resp = RT_RESOURCE_INIT_ZERO;

    if (argc < 1) {
	return RT_APPLICATION_NULL;
    }

    rtip = rt_i_create(gedp->dbip);

    while (0 < argc && argv[0][0] == '-') {
	if (BU_STR_EQUAL(argv[0], "-i")) {
	    rtip->rti_dont_instance = 1;
	    argc--;
	    argv++;
	    continue;
	}
	if (BU_STR_EQUAL(argv[0], "-u")) {
	    rtip->useair = 1;
	    argc--;
	    argv++;
	    continue;
	}
	break;
    }

    if (rt_gettrees(rtip, argc, (const char **)&argv[0], 1) < 0) {
	bu_vls_printf(gedp->ged_result_str, "rt_gettrees() returned error");
	rt_i_destroy(rtip);
	return RT_APPLICATION_NULL;
    }

    /* Establish defaults for this rt_i */
    rtip->rti_hasty_prep = 1;	/* Tcl isn't going to fire many rays */

    /*
     * In case of multiple instances of the library, make sure that
     * each instance has a separate resource structure, because the
     * bit vector lengths depend on # of solids.  And the "overwrite"
     * sequence in Tcl is to create the new proc before running the
     * Tcl_CmdDeleteProc on the old one, which in this case would
     * trash rt_uniresource.  Once on the rti_resources list,
     * rt_clean() will clean 'em up.
     */
    rt_init_resource(&resp, 0, rtip);
    BU_ASSERT(BU_PTBL_GET(&rtip->rti_resources, 0) != NULL);

    BU_ALLOC(ap, struct application);
    RT_APPLICATION_INIT(ap);
    ap->a_magic = RT_AP_MAGIC;
    ap->a_resource = &resp;
    ap->a_rt_i = rtip;
    ap->a_purpose = "Conquest!";

    rt_ck(rtip);

    return ap;
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
