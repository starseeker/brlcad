/*		M G E D _ D I S P L A Y . H
 * BRL-CAD
 *
 * Copyright (c) 1985-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file mged/mged_display.h
 *
 * Internal MGED display, view, input, and GUI-state record.
 *
 */

#ifndef MGED_MGED_DISPLAY_H
#define MGED_MGED_DISPLAY_H

#include "common.h"

#include <string.h>

#include "bv.h"
#include "imgstream/fbserv.h"

#include "pkg.h" /* struct pkg_conn */
#include "ged.h"
#include "ged/view.h"
#include "rt/view.h"

#include "menu.h"
#include "mged.h"

struct scroll_item {
    char *scroll_string;
    void (*scroll_func)(struct scroll_item *, double);
    int scroll_val;
    char *scroll_cmd;
};

#ifndef COMMA
#  define COMMA ','
#endif

#define MGED_DISPLAY_VAR "mged_display"

/* +-2048 to +-1 */
#define GED2PM1(x) (((fastf_t)(x))*BV_INV_VIEW)

#define LAST_SOLID(_sp)       DB_FULL_PATH_CUR_DIR( &(_sp)->s_fullpath )

#define AMM_IDLE 0
#define AMM_ROT 1
#define AMM_TRAN 2
#define AMM_SCALE 3
#define AMM_ADC_ANG1 4
#define AMM_ADC_ANG2 5
#define AMM_ADC_TRAN 6
#define AMM_ADC_DIST 7
#define AMM_CON_ROT_X 8
#define AMM_CON_ROT_Y 9
#define AMM_CON_ROT_Z 10
#define AMM_CON_TRAN_X 11
#define AMM_CON_TRAN_Y 12
#define AMM_CON_TRAN_Z 13
#define AMM_CON_SCALE_X 14
#define AMM_CON_SCALE_Y 15
#define AMM_CON_SCALE_Z 16
#define AMM_CON_XADC 17
#define AMM_CON_YADC 18
#define AMM_CON_ANG1 19
#define AMM_CON_ANG2 20
#define AMM_CON_DIST 21

struct view_ring {
    struct bu_list	l;
    mat_t			vr_rot_mat;
    mat_t			vr_tvc_mat;
    fastf_t		vr_scale;
    int			vr_id;
};


#ifdef MAX_CLIENTS
#	undef MAX_CLIENTS
#	define MAX_CLIENTS 32
#else
#	define MAX_CLIENTS 32
#endif

/* mged command variables for affecting the user environment */
struct _mged_variables {
    int		mv_rc;
    int		mv_autosize;
    int		mv_rateknobs;
    int		mv_sliders;
    int		mv_faceplate;
    int		mv_orig_gui;
    int		mv_linewidth;
    char	mv_linestyle;
    int		mv_hot_key;
    int		mv_context;
    int		mv_use_air;
    int		mv_listen;			/* nonzero to listen on port */
    int		mv_port;			/* port to listen on */
    int		mv_fb;				/* toggle image on/off */
    int		mv_fb_all;			/* 0 - use part of image as defined by the rectangle
						   1 - use the entire image */
    int		mv_fb_overlay;			/* 0 - underlay    1 - interlay    2 - overlay */
    char	mv_mouse_behavior;
    char	mv_coords;
    char	mv_rotate_about;
    char	mv_transform;
    double	mv_perspective;			/* used to directly set the perspective angle */
    int		mv_perspective_mode;		/* used to toggle perspective viewing on/off */
    int		mv_toggle_perspective;		/* used to toggle through values in perspective_table[] */
    double	mv_nmg_eu_dist;
    double	mv_eye_sep_dist;		/* >0 implies stereo.  units = "room" mm */
    char	mv_union_lexeme[1024];
    char	mv_intersection_lexeme[1024];
    char	mv_difference_lexeme[1024];
};


struct _axes_state {
    int		ax_rc;
    int		ax_model_draw;			/* model axes */
    int		ax_model_size;
    int		ax_model_linewidth;
    fastf_t	ax_model_pos[3];
    int		ax_view_draw;			/* view axes */
    int		ax_view_size;
    int		ax_view_linewidth;
    int		ax_view_pos[2];
    int		ax_edit_draw;			/* edit axes */
    int		ax_edit_size1;
    int		ax_edit_size2;
    int		ax_edit_linewidth1;
    int		ax_edit_linewidth2;
};


struct _rubber_band {
    int		rb_rc;
    int		rb_active;	/* 1 - actively drawing a rubber band */
    int		rb_draw;	/* draw rubber band rectangle */
    int		rb_linewidth;
    char	rb_linestyle;
    int		rb_pos[2];	/* Position in image coordinates */
    int		rb_dim[2];	/* Rectangle dimension in image coordinates */
    fastf_t	rb_x;		/* Corner of rectangle in normalized     */
    fastf_t	rb_y;		/* ------ view coordinates (i.e. +-1.0). */
    fastf_t	rb_width;	/* Width and height of rectangle in      */
    fastf_t	rb_height;	/* ------ normalized view coordinates.   */
};


struct _view_state {
    int		vs_rc;

    struct bv_context *vs_gvp;
    mat_t	vs_model2objview;
    mat_t	vs_objview2model;
    mat_t	vs_ModelDelta;		/* changes to Viewrot this frame */

    struct view_ring	vs_headView;
    struct view_ring	*vs_current_view;
    struct view_ring	*vs_last_view;

    /* Rate stuff */
    struct bv_knobs k;

    /* Virtual trackball stuff */
    point_t	vs_orig_pos;
};


struct _color_scheme {
    int	cs_rc;
    int	cs_mode;
    int	cs_bg[3];
    int	cs_bg_a[3];
    int	cs_bg_ia[3];
    int	cs_adc_line[3];
    int	cs_adc_line_a[3];
    int	cs_adc_line_ia[3];
    int	cs_adc_tick[3];
    int	cs_adc_tick_a[3];
    int	cs_adc_tick_ia[3];
    int	cs_geo_def[3];
    int	cs_geo_def_a[3];
    int	cs_geo_def_ia[3];
    int	cs_geo_hl[3];
    int	cs_geo_hl_a[3];
    int	cs_geo_hl_ia[3];
    int	cs_geo_label[3];
    int	cs_geo_label_a[3];
    int	cs_geo_label_ia[3];
    int	cs_model_axes[3];
    int	cs_model_axes_a[3];
    int	cs_model_axes_ia[3];
    int	cs_model_axes_label[3];
    int	cs_model_axes_label_a[3];
    int	cs_model_axes_label_ia[3];
    int	cs_view_axes[3];
    int	cs_view_axes_a[3];
    int	cs_view_axes_ia[3];
    int	cs_view_axes_label[3];
    int	cs_view_axes_label_a[3];
    int	cs_view_axes_label_ia[3];
    int	cs_edit_axes1[3];
    int	cs_edit_axes1_a[3];
    int	cs_edit_axes1_ia[3];
    int	cs_edit_axes_label1[3];
    int	cs_edit_axes_label1_a[3];
    int	cs_edit_axes_label1_ia[3];
    int	cs_edit_axes2[3];
    int	cs_edit_axes2_a[3];
    int	cs_edit_axes2_ia[3];
    int	cs_edit_axes_label2[3];
    int	cs_edit_axes_label2_a[3];
    int	cs_edit_axes_label2_ia[3];
    int	cs_rubber_band[3];
    int	cs_rubber_band_a[3];
    int	cs_rubber_band_ia[3];
    int	cs_grid[3];
    int	cs_grid_a[3];
    int	cs_grid_ia[3];
    int	cs_menu_line[3];
    int	cs_menu_line_a[3];
    int	cs_menu_line_ia[3];
    int	cs_slider_line[3];
    int	cs_slider_line_a[3];
    int	cs_slider_line_ia[3];
    int	cs_other_line[3];
    int	cs_other_line_a[3];
    int	cs_other_line_ia[3];
    int	cs_status_text1[3];
    int	cs_status_text1_a[3];
    int	cs_status_text1_ia[3];
    int	cs_status_text2[3];
    int	cs_status_text2_a[3];
    int	cs_status_text2_ia[3];
    int	cs_slider_text1[3];
    int	cs_slider_text1_a[3];
    int	cs_slider_text1_ia[3];
    int	cs_slider_text2[3];
    int	cs_slider_text2_a[3];
    int	cs_slider_text2_ia[3];
    int	cs_menu_text1[3];
    int	cs_menu_text1_a[3];
    int	cs_menu_text1_ia[3];
    int	cs_menu_text2[3];
    int	cs_menu_text2_a[3];
    int	cs_menu_text2_ia[3];
    int	cs_menu_title[3];
    int	cs_menu_title_a[3];
    int	cs_menu_title_ia[3];
    int	cs_menu_arrow[3];
    int	cs_menu_arrow_a[3];
    int	cs_menu_arrow_ia[3];
    int	cs_state_text1[3];
    int	cs_state_text1_a[3];
    int	cs_state_text1_ia[3];
    int	cs_state_text2[3];
    int	cs_state_text2_a[3];
    int	cs_state_text2_ia[3];
    int	cs_edit_info[3];
    int	cs_edit_info_a[3];
    int	cs_edit_info_ia[3];
    int	cs_center_dot[3];
    int	cs_center_dot_a[3];
    int	cs_center_dot_ia[3];
};


struct _menu_state {
    int	ms_rc;
    int	ms_flag;
    int	ms_top;
    int	ms_cur_menu;
    int	ms_cur_item;
    struct rt_edit_menu_item *ms_menus[NMENU];    /* base of menu items array */
};


struct mged_display {
    struct bu_vls	display_pathname;
    unsigned long	display_native_id;
    int			display_hosted;

    /* Direct endpoint input is per display, not a process-global Tk hook. */
    struct mged_state	*display_input_state;
    int			display_input_motion_pending;
    unsigned long	display_input_motion_timestamp;
    int			display_input_motion_x;
    int			display_input_motion_y;
    int			display_input_snap_pan_active;

    int			repaint_pending;		/* true if this display needs another paint */
    unsigned int	repaint_reasons;		/* audit-only mged_repaint_reason flags */
    int			display_mapped;
    int			display_owner;			/* true if owner of the view info */
    int			display_alternate_mode;			/* alternate mouse mode */
    int			display_perspective_angle;
    int			*display_zclip_ptr;
    struct cmd_list	*display_tie;

    int			display_adc_auto;
    int			display_grid_auto_size;
    int			display_mouse_dx;
    int			display_mouse_dy;
    int			display_pointer_x;
    int			display_pointer_y;
    int			display_knobs[8];
    point_t		display_work_point;

    /* Tcl variable names for display info */
    struct bu_vls	display_fps_name;
    struct bu_vls	display_aet_name;
    struct bu_vls	display_ang_name;
    struct bu_vls	display_center_name;
    struct bu_vls	display_size_name;
    struct bu_vls	display_adc_name;

    /* Slider stuff */
    int			display_scroll_top;
    int			display_scroll_active;
    int			display_scroll_y;
    struct scroll_item	*display_scroll_array[6];

    /* Shareable Resources */
    struct _view_state	*display_view_state;
    struct _menu_state	*display_menu_state;
    struct _rubber_band	*display_rubber_band;
    struct _mged_variables *display_variables;
    struct _color_scheme	*display_color_scheme;
    int			display_color_scheme_dirty;
    struct _axes_state	*display_axes_state;
    int			display_axes_state_dirty;
    int			display_adc_style_dirty;

    /* Hooks */
    int			(*display_command_hook)(int, const char **, void *);
    void			(*display_viewpoint_hook)(void);
    int			(*display_event_handler)(void);
};

enum mged_repaint_reason {
    MGED_REPAINT_NATIVE_EVENT = 1u << 0,
    MGED_REPAINT_VIEW_RECORD = 1u << 1,
    MGED_REPAINT_FRAMEBUFFER = 1u << 2,
    MGED_REPAINT_DEVICE_SETTING = 1u << 3,
    MGED_REPAINT_INTERACTION = 1u << 4
};

static inline void
mged_display_repaint_request(struct mged_display *mdmp, unsigned int reason)
{
    if (!mdmp)
	return;
    mdmp->repaint_pending = 1;
    mdmp->repaint_reasons |= reason;
}

static inline int
mged_display_repaint_pending(const struct mged_display *mdmp)
{
    return mdmp ? mdmp->repaint_pending : 0;
}

static inline const char *
mged_display_pathname(const struct mged_display *mdmp)
{
    return (mdmp && BU_VLS_IS_INITIALIZED(&mdmp->display_pathname) &&
	bu_vls_strlen(&mdmp->display_pathname)) ?
	bu_vls_cstr(&mdmp->display_pathname) : NULL;
}

static inline struct bu_vls *
mged_display_pathname_vls(struct mged_display *mdmp)
{
    return mdmp && BU_VLS_IS_INITIALIZED(&mdmp->display_pathname) ?
	&mdmp->display_pathname : NULL;
}

static inline void
mged_display_repaint_consume(struct mged_display *mdmp)
{
    if (!mdmp)
	return;
    mdmp->repaint_pending = 0;
    mdmp->repaint_reasons = 0;
}

/* Keep GED's active view context synchronized with the active display. */
__BEGIN_DECLS
extern void mged_current_display_set(struct mged_state *s, struct mged_display *nl);
extern void mged_display_adc_state_set(struct mged_display *dm, const struct bv_adc_state *adc);
extern int mged_display_adc_visibility_set(struct mged_display *dm, int enabled);
extern void mged_color_scheme_changed(struct mged_state *s,
	struct _color_scheme *scheme);
extern void mged_obol_faceplate_color_scheme_sync(struct mged_state *s,
	struct mged_display *display);
extern void mged_obol_faceplate_sync(struct mged_state *s,
	struct mged_display *display);
extern int mged_obol_input_motion_consumed(struct mged_display *dm,
	unsigned long timestamp, int x, int y);
extern int mged_obol_framebuffer_ensure(struct mged_state *s);
__END_DECLS

static inline struct bv *
mged_view_state_view(struct _view_state *view_state)
{
    return view_state ? bv_context_view(view_state->vs_gvp) : NULL;
}

static inline const struct bv *
mged_view_state_view_const(const struct _view_state *view_state)
{
    return view_state ? bv_context_view_const(view_state->vs_gvp) : NULL;
}

static inline struct bv *
mged_view_context_view(void *view_ctx)
{
    return bv_context_view((struct bv_context *)view_ctx);
}

static inline const struct bv *
mged_view_context_view_const(const void *view_ctx)
{
    return bv_context_view_const((const struct bv_context *)view_ctx);
}

static inline int
mged_display_adc_state_get(struct mged_display *dm, struct bv_adc_state *adc)
{
    struct bv_adc_state bv_adc;

    if (!dm || !dm->display_view_state || !dm->display_view_state->vs_gvp)
	return 0;
    if (!bv_adc_state_get(&bv_adc,
	    mged_view_state_view_const(dm->display_view_state)))
	return 0;
    memcpy(adc, &bv_adc, sizeof(*adc));
    return 1;
}

static inline int
mged_display_grid_state_get(struct mged_display *dm, struct bv_grid_state *grid)
{
    struct bv_grid_state bv_grid;

    if (!dm || !dm->display_view_state || !dm->display_view_state->vs_gvp)
	return 0;
    if (!bv_grid_state_get(&bv_grid,
	    mged_view_state_view_const(dm->display_view_state)))
	return 0;
    memcpy(grid, &bv_grid, sizeof(*grid));
    return 1;
}

static inline void
mged_display_grid_state_set(struct mged_display *dm, const struct bv_grid_state *grid)
{
    struct bv_grid_state bv_grid;

    if (!dm || !dm->display_view_state || !dm->display_view_state->vs_gvp)
	return;
    memcpy(&bv_grid, grid, sizeof(bv_grid));
    bv_grid_state_set(mged_view_state_view(dm->display_view_state), &bv_grid);
}

static inline int
mged_display_view_settings_shared(struct mged_display *a, struct mged_display *b)
{
    if (!a || !a->display_view_state || !a->display_view_state->vs_gvp ||
	    !b || !b->display_view_state || !b->display_view_state->vs_gvp)
	return 0;
    return bv_context_settings_shared(a->display_view_state->vs_gvp,
	    b->display_view_state->vs_gvp);
}

#define MGED_DISPLAY_NULL ((struct mged_display *)NULL)
#define mapped s->mged_curr_display->display_mapped
#define am_mode s->mged_curr_display->display_alternate_mode
#define perspective_angle s->mged_curr_display->display_perspective_angle
#define zclip_ptr s->mged_curr_display->display_zclip_ptr

#define view_state s->mged_curr_display->display_view_state
#define menu_state s->mged_curr_display->display_menu_state
#define rubber_band s->mged_curr_display->display_rubber_band
#define mged_variables s->mged_curr_display->display_variables
#define color_scheme s->mged_curr_display->display_color_scheme
#define axes_state s->mged_curr_display->display_axes_state

#define cmd_hook s->mged_curr_display->display_command_hook
#define viewpoint_hook s->mged_curr_display->display_viewpoint_hook
#define eventHandler s->mged_curr_display->display_event_handler

#define adc_auto s->mged_curr_display->display_adc_auto
#define grid_auto_size s->mged_curr_display->display_grid_auto_size

/* Shortcuts intentionally differ from the display-record field names. */
#define mouse_dx s->mged_curr_display->display_mouse_dx
#define mouse_dy s->mged_curr_display->display_mouse_dy
#define pointer_x s->mged_curr_display->display_pointer_x
#define pointer_y s->mged_curr_display->display_pointer_y
#define knob_values s->mged_curr_display->display_knobs
#define work_point s->mged_curr_display->display_work_point

#define scroll_top s->mged_curr_display->display_scroll_top
#define scroll_active s->mged_curr_display->display_scroll_active
#define scroll_y s->mged_curr_display->display_scroll_y
#define scroll_array s->mged_curr_display->display_scroll_array

#define VIEWSIZE	(bv_size_get(mged_view_state_view_const(view_state)))	/* Width of viewing cube */
#define VIEWFACTOR	(1/bv_scale_get(mged_view_state_view_const(view_state)))

#define RATE_ROT_FACTOR 6.0
#define ABS_ROT_FACTOR 180.0
#define ADC_ANGLE_FACTOR 45.0

/*
 * Definitions for dealing with the buttons and lights.
 * BV are for viewing, and BE are for editing functions.
 */
#define LIGHT_OFF	0
#define LIGHT_ON	1
#define LIGHT_RESET	2		/* all lights out */

/* Function button/light codes.  Note that code 0 is reserved */
#define BV_TOP		15+16
#define BV_BOTTOM	14+16
#define BV_RIGHT	13+16
#define BV_LEFT		12+16
#define BV_FRONT	11+16
#define BV_REAR		10+16
#define BV_VRESTORE	9+16
#define BV_VSAVE	8+16
#define BE_O_ILLUMINATE	7+16
#define BE_O_SCALE	6+16
#define BE_O_X		5+16
#define BE_O_Y		4+16
#define BE_O_XY		3+16
#define BE_O_ROTATE	2+16
#define BE_ACCEPT	1+16
#define BE_REJECT	0+16

#define BE_S_EDIT	14
#define BE_S_ROTATE	13
#define BE_S_TRANS	12
#define BE_S_SCALE	11
#define BE_MENU		10
#define BV_ADCURSOR	9
#define BV_RESET	8
#define BE_S_ILLUMINATE	7
#define BE_O_XSCALE	6
#define BE_O_YSCALE	5
#define BE_O_ZSCALE	4
#define BV_ZOOM_IN	3
#define BV_ZOOM_OUT	2
#define BV_45_45	1
#define BV_35_25	0+32

#define BV_RATE_TOGGLE	1+32
#define BV_EDIT_TOGGLE  2+32
#define BV_EYEROT_TOGGLE 3+32
#define BE_S_CONTEXT    4+32

#define BV_MAXFUNC	64	/* largest code used */

#define GET_MGED_DISPLAY(p, id) { \
    \
    (p) = MGED_DISPLAY_NULL; \
    for (size_t display_index = 0; display_index < BU_PTBL_LEN(&active_display_set); display_index++) { \
	struct mged_display *tp = (struct mged_display *)BU_PTBL_GET(&active_display_set, display_index); \
	if ((id) == tp->display_native_id) { \
	    (p) = tp; \
	    break; \
	} \
    } \
    \
}

extern double frametime;		/* defined in mged.c */
extern struct bu_ptbl active_display_set;	/* defined in attach.c */
extern struct mged_display *mged_initial_display;

/* defined in doevent.c */
#ifdef HAVE_X11_TYPES
extern int doEvent(ClientData, XEvent *);
#else
extern int doEvent(ClientData, void *);
#endif

/* defined in attach.c */
extern void mged_display_var_init(struct mged_state *s,
	struct mged_display *target_display);

/* defined in display-command.c */
extern int mged_display_command_common(struct mged_state *s, int argc, const char *argv[]);
extern int mged_display_motion(struct mged_state *s, int x, int y);

/* external sp_hook functions */
extern void cs_set_bg(const struct bu_structparse *, const char *, void *, const char *, void *); /* defined in color_scheme.c */

/* defined in setup.c */
extern void mged_rtCmdNotify(int);

int mged_display_command(int argc, const char *argv[], void *data);


#endif /* MGED_MGED_DISPLAY_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
