/*                T C L C A D _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2012-2026 United States Government as represented by
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
/**
 *
 * Private header for libtclcad.
 *
 */

#ifndef LIBTCLCAD_TCLCAD_PRIVATE_H
#define LIBTCLCAD_TCLCAD_PRIVATE_H

#include "common.h"

#include <string.h>
#include <tcl.h>

#include "bg/polygon_types.h"
#include "bu/vls.h"
#include "dm/fbserv.h"
#include "ged/defines.h"
#include "rt/view.h"
#include "tclcad/defines.h"

__BEGIN_DECLS

struct tclcad_data_axes_state {
    int       draw;
    int       color[3];
    int       line_width;
    fastf_t   size;
    int       num_points;
    point_t   *points;
};

struct tclcad_data_arrow_state {
    int       gdas_draw;
    int       gdas_color[3];
    int       gdas_line_width;
    int       gdas_tip_length;
    int       gdas_tip_width;
    int       gdas_num_points;
    point_t   *gdas_points;
};

typedef struct tclcad_data_label_state {
    int         gdls_draw;
    int         gdls_color[3];
    int         gdls_num_labels;
    int         gdls_size;
    char        **gdls_labels;
    point_t     *gdls_points;
} tclcad_label_state;

struct tclcad_data_line_state {
    int       gdls_draw;
    int       gdls_color[3];
    int       gdls_line_width;
    int       gdls_num_points;
    point_t   *gdls_points;
};

typedef struct tclcad_polygon_state {
    int                 gdps_draw;
    int                 gdps_moveAll;
    int                 gdps_color[3];
    int                 gdps_line_width;
    int                 gdps_line_style;
    int                 gdps_cflag;
    size_t              gdps_target_polygon_i;
    size_t              gdps_curr_polygon_i;
    size_t              gdps_curr_point_i;
    point_t             gdps_prev_point;
    bg_clip_t           gdps_clip_type;
    fastf_t             gdps_scale;
    point_t             gdps_origin;
    mat_t               gdps_rotation;
    mat_t               gdps_view2model;
    mat_t               gdps_model2view;
    struct bg_polygons  gdps_polygons;
    fastf_t             gdps_data_vZ;
} tclcad_polygon_state;

typedef struct tclcad_view_state {
    int				  gv_polygon_mode;
    int				  gv_hide;
    fastf_t			  gv_data_vZ;
    struct tclcad_data_arrow_state	  gv_data_arrows;
    struct tclcad_data_axes_state	  gv_data_axes;
    tclcad_label_state			  gv_data_labels;
    struct tclcad_data_line_state	  gv_data_lines;
    tclcad_polygon_state		  gv_data_polygons;
    struct tclcad_data_arrow_state	  gv_sdata_arrows;
    struct tclcad_data_axes_state	  gv_sdata_axes;
    tclcad_label_state			  gv_sdata_labels;
    struct tclcad_data_line_state	  gv_sdata_lines;
    tclcad_polygon_state		  gv_sdata_polygons;
    struct bv_other_state		  gv_prim_labels;
} tclcad_view_state;


#define TO_UNLIMITED -1

typedef int (*to_wrapper_func_ptr)(struct ged *, int, const char *[], ged_func_ptr, const char *, int);
#define TO_WRAPPER_FUNC_PTR_NULL (to_wrapper_func_ptr)0

struct to_cmdtab {
    const char *to_name;
    const char *to_usage;
    int to_maxargs;
    to_wrapper_func_ptr to_wrapper_func;
    ged_func_ptr to_func;
};

// For the test program check_tclcad_cmds
TCLCAD_EXPORT extern struct to_cmdtab to_cmds[];

extern struct tclcad_obj HeadTclcadObj;
extern struct tclcad_obj *current_top;

// Data specific to an individual view rather than the geometry database
// instance.
struct tclcad_view_data {
    struct ged *gedp;
    struct bu_vls gdv_edit_motion_delta_callback;
    int gdv_edit_motion_delta_callback_cnt;
    struct bu_vls gdv_callback;
    int gdv_callback_cnt;
    struct fbserv_obj gdv_fbs;
    /* Tcl-specific overlay state owned by TclCAD view data.  The view's
     * Tcl pointer is bound to this record while the view is live and cleared
     * when the view is freed. */
    tclcad_view_state tcl_data;
};

/**
 * function returns truthfully whether the library has been
 * initialized.  calling this routine with setit true considers the
 * library henceforth initialized.  there is presently no way to unset
 * or reset initialization.
 */
extern int library_initialized(int setit);


/**
 * Evaluates a TCL command, escaping the list of arguments.
 */
extern int tclcad_eval(Tcl_Interp *interp, const char *command, size_t num_args, const char * const *args);


/**
 * Evaluates a TCL command, escaping the list of arguments and preserving the TCL result object.
 */
extern int tclcad_eval_noresult(Tcl_Interp *interp, const char *command, size_t num_args, const char * const *args);

extern void tclcad_view_data_init(struct tclcad_view_data *tvd, struct ged *gedp);
extern struct tclcad_view_data *tclcad_view_data_from_view_ctx(void *view_ctx);
extern tclcad_view_state *tclcad_view_tcl_data_from_view_ctx(void *view_ctx);
extern tclcad_polygon_state *tclcad_view_polygon_state_from_view_ctx(void *view_ctx, int staged);
extern int tclcad_view_polygon_mode_from_view_ctx(void *view_ctx);
extern int tclcad_view_polygon_mode_set(void *view_ctx, int mode);
extern fastf_t tclcad_view_data_vZ_from_view_ctx(void *view_ctx);
extern int tclcad_view_data_vZ_set(void *view_ctx, fastf_t vZ);
extern int tclcad_view_hide_from_view_ctx(void *view_ctx);
extern int tclcad_view_polygon_cflag_from_view_ctx(void *view_ctx, int staged);
extern int tclcad_view_polygon_cflag_clear(void *view_ctx, int staged);
extern tclcad_label_state *tclcad_view_label_state_from_view_ctx(void *view_ctx, int staged);
extern int tclcad_view_prim_labels_state_from_view_ctx(struct bv_other_state *state, void *view_ctx);
extern int tclcad_view_prim_labels_state_set(void *view_ctx, const struct bv_other_state *state);
extern int tclcad_view_data_bind_view_ctx(void *view_ctx, struct tclcad_view_data *tvd);
extern void tclcad_view_data_unbind_view_ctx(void *view_ctx);


/* Tcl initialization routines */
TCLCAD_EXPORT extern int Bu_Init(Tcl_Interp *interp);
TCLCAD_EXPORT extern int Bn_Init(Tcl_Interp *interp);
TCLCAD_EXPORT extern int Dm_Init(Tcl_Interp *interp);
TCLCAD_EXPORT extern int Fbo_Init(Tcl_Interp *interp);
TCLCAD_EXPORT extern int Ged_Init(Tcl_Interp *interp);

/* Fb functions */
extern int to_close_fbs(void *view_ctx);
extern void to_fbs_callback(void *);
extern int to_open_fbs(void *view_ctx, Tcl_Interp *interp);
extern int to_set_fb_mode(struct ged *gedp,
			  int argc,
			  const char *argv[],
			  ged_func_ptr func,
			  const char *usage,
			  int maxargs);
extern int to_listen(struct ged *gedp,
		     int argc,
		     const char *argv[],
		     ged_func_ptr func,
		     const char *usage,
		     int maxargs);

/* Dm functions */
extern int
dm_list_tcl(ClientData clientData,
	    Tcl_Interp *interp,
	    int argc,
	    const char **argv);

/* Tclcad mouse routines */
extern int to_get_prev_mouse(struct ged *gedp,
                            int argc,
                            const char *argv[],
                            ged_func_ptr func,
                            const char *usage,
                            int maxargs);
extern int to_mouse_append_pnt_common(struct ged *gedp,
                                     int argc,
                                     const char *argv[],
                                     ged_func_ptr func,
                                     const char *usage,
                                     int maxargs);
extern int to_mouse_brep_selection_append(struct ged *gedp,
                                         int argc,
                                         const char *argv[],
                                         ged_func_ptr func,
                                         const char *usage,
                                         int maxargs);
extern int to_mouse_brep_selection_translate(struct ged *gedp,
                                            int argc,
                                            const char *argv[],
                                            ged_func_ptr func,
                                            const char *usage,
                                            int maxargs);
extern int to_mouse_constrain_rot(struct ged *gedp,
                                 int argc,
                                 const char *argv[],
                                 ged_func_ptr func,
                                 const char *usage,
                                 int maxargs);
extern int to_mouse_constrain_trans(struct ged *gedp,
                                   int argc,
                                   const char *argv[],
                                   ged_func_ptr func,
                                   const char *usage,
                                   int maxargs);
extern int to_mouse_find_arb_edge(struct ged *gedp,
                                 int argc,
                                 const char *argv[],
                                 ged_func_ptr func,
                                 const char *usage,
                                 int maxargs);
extern int to_mouse_find_bot_edge(struct ged *gedp,
                                 int argc,
                                 const char *argv[],
                                 ged_func_ptr func,
                                 const char *usage,
                                 int maxargs);
extern int to_mouse_find_bot_pnt(struct ged *gedp,
                                int argc,
                                const char *argv[],
                                ged_func_ptr func,
                                const char *usage,
                                int maxargs);
extern int to_mouse_find_metaball_pnt(struct ged *gedp,
                                     int argc,
                                     const char *argv[],
                                     ged_func_ptr func,
                                     const char *usage,
                                     int maxargs);
extern int to_mouse_find_pipe_pnt(struct ged *gedp,
                                 int argc,
                                 const char *argv[],
                                 ged_func_ptr func,
                                 const char *usage,
                                 int maxargs);
extern int to_mouse_joint_select(struct ged *gedp,
                                int argc,
                                const char *argv[],
                                ged_func_ptr func,
                                const char *usage,
                                int maxargs);
extern int to_mouse_joint_selection_translate(struct ged *gedp,
                                             int argc,
                                             const char *argv[],
                                             ged_func_ptr func,
                                             const char *usage,
                                             int maxargs);
extern int to_mouse_move_arb_edge(struct ged *gedp,
                                 int argc,
                                 const char *argv[],
                                 ged_func_ptr func,
                                 const char *usage,
                                 int maxargs);
extern int to_mouse_move_arb_face(struct ged *gedp,
                                 int argc,
                                 const char *argv[],
                                 ged_func_ptr func,
                                 const char *usage,
                                 int maxargs);
extern int to_mouse_move_bot_pnt(struct ged *gedp,
                                int argc,
                                const char *argv[],
                                ged_func_ptr func,
                                const char *usage,
                                int maxargs);
extern int to_mouse_move_bot_pnts(struct ged *gedp,
                                 int argc,
                                 const char *argv[],
                                 ged_func_ptr func,
                                 const char *usage,
                                 int maxargs);
extern int to_mouse_move_pnt_common(struct ged *gedp,
                                   int argc,
                                   const char *argv[],
                                   ged_func_ptr func,
                                   const char *usage,
                                   int maxargs);
extern int to_mouse_orotate(struct ged *gedp,
                           int argc,
                           const char *argv[],
                           ged_func_ptr func,
                           const char *usage,
                           int maxargs);
extern int to_mouse_oscale(struct ged *gedp,
                          int argc,
                          const char *argv[],
                          ged_func_ptr func,
                          const char *usage,
                          int maxargs);
extern int to_mouse_otranslate(struct ged *gedp,
                              int argc,
                              const char *argv[],
                              ged_func_ptr func,
                              const char *usage,
                              int maxargs);
extern int to_mouse_poly_circ(struct ged *gedp,
                             int argc,
                             const char *argv[],
                             ged_func_ptr func,
                             const char *usage,
                             int maxargs);
extern int to_mouse_poly_circ_func(Tcl_Interp *interp,
                                  struct ged *gedp,
                                  void *view_ctx,
                                  int argc,
                                  const char *argv[],
                                  const char *usage);
extern int to_mouse_poly_cont(struct ged *gedp,
                             int argc,
                             const char *argv[],
                             ged_func_ptr func,
                             const char *usage,
                             int maxargs);
extern int to_mouse_poly_cont_func(Tcl_Interp *interp,
                                  struct ged *gedp,
                                  void *view_ctx,
                                  int argc,
                                  const char *argv[],
                                  const char *usage);
extern int to_mouse_poly_ell(struct ged *gedp,
                            int argc,
                            const char *argv[],
                            ged_func_ptr func,
                            const char *usage,
                            int maxargs);
extern int to_mouse_poly_ell_func(Tcl_Interp *interp,
                                 struct ged *gedp,
                                 void *view_ctx,
                                 int argc,
                                 const char *argv[],
                                 const char *usage);
extern int to_mouse_poly_rect(struct ged *gedp,
                             int argc,
                             const char *argv[],
                             ged_func_ptr func,
                             const char *usage,
                             int maxargs);
extern int to_mouse_poly_rect_func(Tcl_Interp *interp,
                                  struct ged *gedp,
                                  void *view_ctx,
                                  int argc,
                                  const char *argv[],
                                  const char *usage);
extern int to_mouse_ray(struct ged *gedp,
                       int argc,
                       const char *argv[],
                       ged_func_ptr func,
                       const char *usage,
                       int maxargs);
extern int to_mouse_rect(struct ged *gedp,
                        int argc,
                        const char *argv[],
                        ged_func_ptr func,
                        const char *usage,
                        int maxargs);
extern int to_mouse_rot(struct ged *gedp,
                       int argc,
                       const char *argv[],
                       ged_func_ptr func,
                       const char *usage,
                       int maxargs);
extern int to_mouse_rotate_arb_face(struct ged *gedp,
                                   int argc,
                                   const char *argv[],
                                   ged_func_ptr func,
                                   const char *usage,
                                   int maxargs);
extern int to_mouse_data_scale(struct ged *gedp,
                              int argc,
                              const char *argv[],
                              ged_func_ptr func,
                              const char *usage,
                              int maxargs);
extern int to_mouse_scale(struct ged *gedp,
                         int argc,
                         const char *argv[],
                         ged_func_ptr func,
                         const char *usage,
                         int maxargs);
extern int to_mouse_protate(struct ged *gedp,
                           int argc,
                           const char *argv[],
                           ged_func_ptr func,
                           const char *usage,
                           int maxargs);
extern int to_mouse_pscale(struct ged *gedp,
                          int argc,
                          const char *argv[],
                          ged_func_ptr func,
                          const char *usage,
                          int maxargs);
extern int to_mouse_ptranslate(struct ged *gedp,
                              int argc,
                              const char *argv[],
                              ged_func_ptr func,
                              const char *usage,
                              int maxargs);
extern int to_mouse_trans(struct ged *gedp,
                         int argc,
                         const char *argv[],
                         ged_func_ptr func,
                         const char *usage,
                         int maxargs);

/* Tclcad polygon routines */
extern int to_data_polygons_func(Tcl_Interp *interp,
                                 struct ged *gedp,
                                 void *view_ctx,
                                 int argc,
                                 const char *argv[]);
extern int to_data_polygons(struct ged *gedp,
			    int argc,
			    const char *argv[],
			    ged_func_ptr func,
			    const char *usage,
			    int maxargs);

extern int to_poly_circ_mode(struct ged *gedp,
			     int argc,
			     const char *argv[],
			     ged_func_ptr func,
			     const char *usage,
			     int maxargs);
extern int to_poly_circ_mode_func(Tcl_Interp *interp,
				  struct ged *gedp,
				  void *view_ctx,
				  int argc,
				  const char *argv[],
				  const char *usage);
extern int to_poly_cont_build(struct ged *gedp,
			      int argc,
			      const char *argv[],
			      ged_func_ptr func,
			      const char *usage,
			      int maxargs);
extern int to_poly_cont_build_end(struct ged *gedp,
				  int argc,
				  const char *argv[],
				  ged_func_ptr func,
				  const char *usage,
				  int maxargs);
extern int to_poly_cont_build_end_func(void *view_ctx,
				       int argc,
				       const char *argv[]);
extern int to_poly_ell_mode(struct ged *gedp,
			    int argc,
			    const char *argv[],
			    ged_func_ptr func,
			    const char *usage,
			    int maxargs);
extern int to_poly_ell_mode_func(Tcl_Interp *interp,
				 struct ged *gedp,
				 void *view_ctx,
				 int argc,
				 const char *argv[],
				 const char *usage);
extern int to_poly_rect_mode(struct ged *gedp,
			     int argc,
			     const char *argv[],
			     ged_func_ptr func,
			     const char *usage,
			     int maxargs);
extern int to_poly_rect_mode_func(Tcl_Interp *interp,
				  struct ged *gedp,
				  void *view_ctx,
				  int argc,
				  const char *argv[],
				  const char *usage);


/* Tclcad obj wrapper routines */
extern int to_pass_through_func(struct ged *gedp,
				int argc,
				const char *argv[],
				ged_func_ptr func,
				const char *usage,
				int maxargs);
extern int to_pass_through_and_refresh_func(struct ged *gedp,
					    int argc,
					    const char *argv[],
					    ged_func_ptr func,
					    const char *usage,
					    int maxargs);
extern int to_view_func(struct ged *gedp,
			int argc,
			const char *argv[],
			ged_func_ptr func,
			const char *usage,
			int maxargs);
extern int to_view_func_common(struct ged *gedp,
			       int argc,
			       const char *argv[],
			       ged_func_ptr func,
			       const char *usage,
			       int maxargs,
			       int cflag,
			       int rflag);
extern int to_view_func_less(struct ged *gedp,
			     int argc,
			     const char *argv[],
			     ged_func_ptr func,
			     const char *usage,
			     int maxargs);
extern int to_view_func_plus(struct ged *gedp,
			     int argc,
			     const char *argv[],
			     ged_func_ptr func,
			     const char *usage,
			     int maxargs);
extern int to_dm_func(struct ged *gedp,
		      int argc,
		      const char *argv[],
		      ged_func_ptr func,
		      const char *usage,
		      int maxargs);

extern int to_more_args_func(struct ged *gedp,
                             int argc,
                             const char *argv[],
                             ged_func_ptr func,
                             const char *usage,
                             int maxargs);

extern int
tclcad_is_listening(struct fbserv_obj *fbsp);
extern int
tclcad_listen_on_port(struct fbserv_obj *fbsp, int available_port);
extern void
tclcad_open_server_handler(struct fbserv_obj *fbsp);
extern void
tclcad_close_server_handler(struct fbserv_obj *fbsp);
extern void
tclcad_open_client_handler(struct fbserv_obj *fbsp, int i, void *data);
extern void
tclcad_close_client_handler(struct fbserv_obj *fbsp, int sub);
TCLCAD_EXPORT extern void
tclcad_fbserv_set_transport(struct fbserv_obj *fbsp);
TCLCAD_EXPORT extern int
tclcad_listen_ipc(struct fbserv_obj *fbsp, Tcl_Interp *interp);

__END_DECLS

#endif /* LIBTCLCAD_TCLCAD_PRIVATE_H */

/** @} */
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
