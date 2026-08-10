/*                           I S S T  . C
 * BRL-CAD
 *
 * Copyright (c) 2005-2026 United States Government as represented by
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
/** @file isst.c
 *
 *  Tcl/Tk version of ISST.
 *
 */

#include "common.h"

#include "bnetwork.h"
#include "bio.h"

#include <limits.h>
#include <stdint.h>

#include "tcl.h"

#include "bu/app.h"
#include "bu/parallel.h"
#include "bu/datetime.h"

#include "rt/tie.h"
#include "adrt.h"
#include "adrt_struct.h"
#include "camera.h"
#include "raytrace.h"
#include "tclcad.h"

// Tclcad pulls in OpenNURBS in C++ compilation mode, which defines None, which
// will conflict with Tk.h's Xlib None if we include tk.h before tclcad.h
#include "tk.h"

#ifdef HAVE_STRING_H
#include <string.h>
#endif

struct isst_s {
    struct tie_s *tie;
    struct render_camera_s camera;
    struct camera_tile_s tile;
    struct adrt_mesh_s *meshes;
    tienet_buffer_t buffer_image;
    int w, h, gs;
    Tk_PhotoHandle photo;
    Tcl_Interp *photo_interp;
    unsigned char *display_pixels;
    size_t display_pixels_size;
    int photo_width;
    int photo_height;
    vect_t camera_pos_init;
    vect_t camera_focus_init;
    int64_t t1;
    int64_t t2;
    int dirty;
};

static struct isst_s *isst;

static int resize_isst(struct isst_s *);

static int
isst_error(Tcl_Interp *interp, const char *message)
{
    Tcl_SetObjResult(interp, Tcl_NewStringObj(message, -1));
    return TCL_ERROR;
}

static int
isst_display_buffer_resize(struct isst_s *isstp, size_t width, size_t height)
{
    if (!isstp || !width || !height || width > SIZE_MAX / height ||
	width * height > SIZE_MAX / 3u)
	return 0;

    const size_t bytes = width * height * 3u;
    if (bytes > isstp->display_pixels_size) {
	isstp->display_pixels = (unsigned char *)bu_realloc(
	    isstp->display_pixels, bytes, "isst Tk photo pixels");
	isstp->display_pixels_size = bytes;
    }
    return isstp->display_pixels ? 1 : 0;
}

static int
isst_photo_present(struct isst_s *isstp)
{
    if (!isstp || !isstp->photo || !isstp->photo_interp ||
	isstp->camera.w <= 0 || isstp->camera.h <= 0 ||
	isstp->w <= 0 || isstp->h <= 0 || !isstp->buffer_image.data)
	return 0;

    const size_t source_width = (size_t)isstp->camera.w;
    const size_t source_height = (size_t)isstp->camera.h;
    const size_t target_width = (size_t)isstp->w;
    const size_t target_height = (size_t)isstp->h;
    const unsigned char *source = isstp->buffer_image.data +
	sizeof(camera_tile_t);
    unsigned char *pixels = (unsigned char *)source;

    /* ISST has always rendered at a selectable lower resolution and let its
     * presentation layer scale the image to the Tk window.  Keep that quality
     * control without requiring an OpenGL texture or display manager. */
    if (source_width != target_width || source_height != target_height) {
	if (!isst_display_buffer_resize(isstp, target_width, target_height))
	    return 0;
	for (size_t y = 0; y < target_height; y++) {
	    const size_t source_y = y * source_height / target_height;
	    for (size_t x = 0; x < target_width; x++) {
		const size_t source_x = x * source_width / target_width;
		const unsigned char *src = source +
		    (source_y * source_width + source_x) * 3u;
		unsigned char *dst = isstp->display_pixels +
		    (y * target_width + x) * 3u;
		dst[0] = src[0];
		dst[1] = src[1];
		dst[2] = src[2];
	    }
	}
	pixels = isstp->display_pixels;
    }

    Tk_PhotoImageBlock block;
    block.pixelPtr = pixels;
    block.width = isstp->w;
    block.height = isstp->h;
    block.pitch = isstp->w * 3;
    block.pixelSize = 3;
    block.offset[0] = 0;
    block.offset[1] = 1;
    block.offset[2] = 2;
    block.offset[3] = 0;

#if TK_MAJOR_VERSION < 9
    if ((isstp->photo_width != block.width ||
	 isstp->photo_height != block.height) &&
	Tk_PhotoSetSize(isstp->photo_interp, isstp->photo,
	    block.width, block.height) != TCL_OK)
	return 0;
    if (Tk_PhotoPutBlock(isstp->photo_interp, isstp->photo, &block,
	0, 0, block.width, block.height, TK_PHOTO_COMPOSITE_SET) != TCL_OK)
	return 0;
#else
    if (isstp->photo_width != block.width ||
	isstp->photo_height != block.height)
	Tk_PhotoSetSize(isstp->photo, block.width, block.height);
    Tk_PhotoPutBlock(isstp->photo, &block, 0, 0, block.width, block.height,
	TK_PHOTO_COMPOSITE_SET);
#endif
    isstp->photo_width = block.width;
    isstp->photo_height = block.height;
    return 1;
}

/* new window size or exposure */
static int
reshape(ClientData UNUSED(clientData), Tcl_Interp *interp, int objc,
	Tcl_Obj *const *objv)
{
    int width, height;
    if (objc != 3) {
	Tcl_WrongNumArgs(interp, 1, objv, "width height");
	return TCL_ERROR;
    }

    if (Tcl_GetIntFromObj(interp, objv[1], &width) != TCL_OK ||
	Tcl_GetIntFromObj(interp, objv[2], &height) != TCL_OK)
	return TCL_ERROR;
    if (width <= 0 || height <= 0)
	return isst_error(interp, "ISST window dimensions must be positive");

    isst->w = width;
    isst->h = height;

    if (!resize_isst(isst))
	return isst_error(interp, "ISST image dimensions are too large");

    return TCL_OK;
}

static int
resize_isst(struct isst_s *isstp)
{
    size_t render_width;
    size_t render_height;

    if (!isstp || isstp->w <= 0 || isstp->h <= 0)
	return 0;

    /* Tk's photo pitch is an int and the ADRT tile dimensions are uint16_t. */
    if (isstp->w > INT_MAX / 3)
	return 0;

    switch (isstp->gs) {
	case 0:
	    render_width = (size_t)isstp->w;
	    render_height = (size_t)isstp->h;
	    break;
	default:
	    render_width = (size_t)isstp->gs;
	    if (render_width > SIZE_MAX / (size_t)isstp->h)
		return 0;
	    render_height = render_width * (size_t)isstp->h /
		(size_t)isstp->w;
	    break;
    }


    if (!render_width || !render_height || render_width > UINT16_MAX ||
	render_height > UINT16_MAX || render_width > SIZE_MAX / render_height ||
	render_width * render_height >
	(SIZE_MAX - sizeof(camera_tile_t)) / 3u)
	return 0;

    isstp->camera.w = isstp->tile.size_x = (uint16_t)render_width;
    isstp->camera.h = isstp->tile.size_y = (uint16_t)render_height;

    isstp->tile.format = RENDER_CAMERA_BIT_DEPTH_24;
    const size_t render_bytes = (size_t)isstp->camera.w *
	(size_t)isstp->camera.h * 3u + sizeof(camera_tile_t);
    if (render_bytes > UINT32_MAX ||
	!isst_display_buffer_resize(isstp, (size_t)isstp->w,
	    (size_t)isstp->h))
	return 0;
    TIENET_BUFFER_SIZE(isstp->buffer_image, (uint32_t)render_bytes);
    isstp->dirty = 1;
    return 1;
}


static int
isst_load_g(ClientData UNUSED(clientData), Tcl_Interp *interp, int objc,
	    Tcl_Obj *const *objv)
{
    char **argv;
    int argc;
    double az, el;
    struct bu_vls tclstr = BU_VLS_INIT_ZERO;
    vect_t vec;

    if (objc < 4) {
	Tcl_WrongNumArgs(interp, 1, objv, "load_g pathname object");
	return TCL_ERROR;
    }

    argv = (char **)malloc(sizeof(char *) * (strlen(Tcl_GetString(objv[3])) + 1));	/* allocate way too much. */
    argc = (int)bu_argv_from_string(argv, strlen(Tcl_GetString(objv[3])), Tcl_GetString(objv[3]));

    load_g(isst->tie, Tcl_GetString(objv[2]), argc, (const char **)argv, &(isst->meshes));
    free(argv);

    VSETALL(isst->camera.pos, isst->tie->radius);
    VMOVE(isst->camera.focus, isst->tie->mid);
    VMOVE(isst->camera_pos_init, isst->camera.pos);
    VMOVE(isst->camera_focus_init, isst->camera.focus);

    /* Set the initial az and el values in Tcl/Tk */
    VSUB2(vec, isst->camera.pos, isst->camera.focus);
    VUNITIZE(vec);
    bn_ae_vec(&az, &el, vec);
    az = az * -DEG2RAD;
    el = el * -DEG2RAD;
    bu_vls_sprintf(&tclstr, "%f", az);
    Tcl_SetVar(interp, "az", bu_vls_addr(&tclstr), 0);
    bu_vls_sprintf(&tclstr, "%f", el);
    Tcl_SetVar(interp, "el", bu_vls_addr(&tclstr), 0);
    bu_vls_free(&tclstr);

    render_phong_init(&isst->camera.render, NULL);

    isst->w = 800;
    isst->h = 600;

    if (!resize_isst(isst))
	return isst_error(interp, "ISST image dimensions are too large");

    isst->t1 = bu_gettime();
    isst->t2 = bu_gettime();

    return TCL_OK;
}

static int
list_geometry(ClientData UNUSED(clientData), Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
{
    static struct db_i *dbip;
    struct directory *dp;
    struct bu_vls tclstr = BU_VLS_INIT_ZERO;

    if (objc < 3) {
	Tcl_WrongNumArgs(interp, 1, objv, "file varname");
	return TCL_ERROR;
    }

    if ((dbip = db_open(Tcl_GetString(objv[1]), DB_OPEN_READONLY)) == DBI_NULL) {
	bu_log("Unable to open geometry database file (%s)\n", Tcl_GetString(objv[1]));
	return TCL_ERROR;
    }
    db_dirbuild(dbip);
    FOR_ALL_DIRECTORY_START(dp, dbip)
	if (dp->d_flags & RT_DIR_HIDDEN) continue;
	bu_vls_sprintf(&tclstr, "set %s [concat $%s [list %s]]", Tcl_GetString(objv[2]), Tcl_GetString(objv[2]), dp->d_namep);
	Tcl_Eval(interp, bu_vls_addr(&tclstr));
    FOR_ALL_DIRECTORY_END;
    db_close(dbip);
    bu_vls_free(&tclstr);
    return TCL_OK;
}


static int
paint_window(ClientData UNUSED(clientData), Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
{
    double dt = 1.0;

    if (objc != 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "pathName");
	return TCL_ERROR;
    }

    isst->t2 = bu_gettime();

    dt = isst->t2 - isst->t1;

    if (dt > 1e6*0.08 && isst->dirty) {
	isst->buffer_image.ind = 0;

	render_camera_prep(&isst->camera);
	render_camera_render(&isst->camera, isst->tie, &isst->tile, &isst->buffer_image);

	isst->t1 = bu_gettime();

	if (!isst_photo_present(isst))
	    return isst_error(interp, "ISST could not present the rendered image");

	isst->dirty = 0;
    }
    return TCL_OK;
}

static int
set_resolution(ClientData UNUSED(clientData), Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
{
    int resolution = 0;

    if (objc != 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "resolution");
	return TCL_ERROR;
    }

    if (Tcl_GetIntFromObj(interp, objv[1], &resolution) != TCL_OK) {
	return TCL_ERROR;
    }

    CLAMP(resolution, 1, 20);
    if (resolution == 20)
	isst->gs = 0;
    else
	isst->gs = lrint(floor(isst->w * .05 * resolution));

    if (!resize_isst(isst))
	return isst_error(interp, "ISST image dimensions are too large");

    return TCL_OK;
}


static int
isst_init(ClientData UNUSED(clientData), Tcl_Interp *interp, int UNUSED(objc), Tcl_Obj *const *UNUSED(objv))
{
    BU_ALLOC(isst, struct isst_s);
    isst->w = 800;
    isst->h = 600;

    BU_ALLOC(isst->tie, struct tie_s);
    TIENET_BUFFER_INIT(isst->buffer_image);
    render_camera_init(&isst->camera, bu_avail_cpus());

    isst->camera.type = RENDER_CAMERA_PERSPECTIVE;
    isst->camera.fov = 25;

    Tcl_LinkVar(interp, "pos_x", (char *)&isst->camera.pos[0], TCL_LINK_DOUBLE);
    Tcl_LinkVar(interp, "pos_y", (char *)&isst->camera.pos[1], TCL_LINK_DOUBLE);
    Tcl_LinkVar(interp, "pos_z", (char *)&isst->camera.pos[2], TCL_LINK_DOUBLE);
    Tcl_LinkVar(interp, "focus_x", (char *)&isst->camera.focus[0], TCL_LINK_DOUBLE);
    Tcl_LinkVar(interp, "focus_y", (char *)&isst->camera.focus[1], TCL_LINK_DOUBLE);
    Tcl_LinkVar(interp, "focus_z", (char *)&isst->camera.focus[2], TCL_LINK_DOUBLE);

    return TCL_OK;
}

static int
isst_bind_photo(ClientData UNUSED(clientData), Tcl_Interp *interp, int objc,
	Tcl_Obj *const *objv)
{
    if (objc != 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "imageName");
	return TCL_ERROR;
    }
    if (!isst)
	return isst_error(interp, "ISST is not initialized");

    isst->photo = Tk_FindPhoto(interp, Tcl_GetString(objv[1]));
    if (!isst->photo)
	return isst_error(interp, "ISST requires a Tk photo image");
    isst->photo_interp = interp;
    isst->photo_width = 0;
    isst->photo_height = 0;
    return TCL_OK;
}

static int
isst_zap(ClientData UNUSED(clientData), Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
{
    if (objc != 2) {
	Tcl_WrongNumArgs(interp, 1, objv, "pathName");
	return TCL_ERROR;
    }

    if (isst) {
	if (isst->buffer_image.data)
	    bu_free(isst->buffer_image.data, "isst render pixels");
	if (isst->display_pixels)
	    bu_free(isst->display_pixels, "isst Tk photo pixels");
	bu_free(isst, "isst free");
    }
    isst = NULL;

    return TCL_OK;
}

static int
render_mode(ClientData UNUSED(clientData), Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
{
    char *buf = NULL;
    char *mode;

    if (objc < 3) {
	Tcl_WrongNumArgs(interp, 1, objv, "pathName mode [arguments]");
	return TCL_ERROR;
    }

    mode = Tcl_GetString(objv[2]);
    if (objc == 4)
	buf = Tcl_GetString(objv[3]);

    /* pack the 'rest' into buf - probably should use a vls for this*/
    if ( strlen(mode) == 3 && bu_strncmp("cut", mode, 3) == 0 ) {
	struct adrt_mesh_s *mesh;

	/* clear all the hit list */
	for (BU_LIST_FOR(mesh, adrt_mesh_s, &isst->meshes->l))
	    mesh->flags &= ~ADRT_MESH_HIT;
    }

    if (render_shader_init(&isst->camera.render, mode, buf) != 0)
	return TCL_ERROR;

    isst->dirty = 1;

    return TCL_OK;
}


static int
zero_view(ClientData UNUSED(clientData), Tcl_Interp *UNUSED(interp), int UNUSED(objc), Tcl_Obj *const *UNUSED(objv))
{
    vect_t vec;
    double mag_vec;

    mag_vec = DIST_PNT_PNT(isst->camera.pos, isst->camera.focus);

    VSUB2(vec, isst->camera_focus_init, isst->camera.pos);
    VUNITIZE(vec);
    VSCALE(vec, vec, mag_vec);
    VADD2(isst->camera.focus, isst->camera.pos, vec);

    isst->dirty = 1;
    return TCL_OK;
}


static int
move_walk(ClientData UNUSED(clientData), Tcl_Interp *interp, int UNUSED(objc), Tcl_Obj *const *objv)
{
    vect_t vec;
    int flag;

    if (Tcl_GetIntFromObj(interp, objv[2], &flag) != TCL_OK)
	return TCL_ERROR;

    if (flag >= 0) {
	VSUB2(vec, isst->camera.focus, isst->camera.pos);
	VSCALE(vec, vec, 0.1 * isst->tie->radius);
	VADD2(isst->camera.pos, isst->camera.pos, vec);
	VADD2(isst->camera.focus, isst->camera.focus, vec);
    } else {
	VSUB2(vec, isst->camera.pos, isst->camera.focus);
	VSCALE(vec, vec, 0.1 * isst->tie->radius);
	VADD2(isst->camera.pos, isst->camera.pos, vec);
	VADD2(isst->camera.focus, isst->camera.focus, vec);
    }
    isst->dirty = 1;
    return TCL_OK;
}

static int
move_strafe(ClientData UNUSED(clientData), Tcl_Interp *interp, int UNUSED(objc), Tcl_Obj *const *objv)
{
    vect_t vec, dir, up;

    int flag;

    if (Tcl_GetIntFromObj(interp, objv[2], &flag) != TCL_OK)
	return TCL_ERROR;

    VSET(up, 0, 0, 1);

    if (flag >= 0) {
	VSUB2(dir, isst->camera.focus, isst->camera.pos);
	VCROSS(vec, dir, up);
	VSCALE(vec, vec, 0.1 * isst->tie->radius);
	VADD2(isst->camera.pos, isst->camera.pos, vec);
	VADD2(isst->camera.focus, isst->camera.pos, dir);
    } else {
	VSUB2(dir, isst->camera.focus, isst->camera.pos);
	VCROSS(vec, dir, up);
	VSCALE(vec, vec, -0.1 * isst->tie->radius);
	VADD2(isst->camera.pos, isst->camera.pos, vec);
	VADD2(isst->camera.focus, isst->camera.pos, dir);
    }

    isst->dirty = 1;
    return TCL_OK;
}

static int
move_float(ClientData UNUSED(clientData), Tcl_Interp *UNUSED(interp), int UNUSED(objc), Tcl_Obj *const *UNUSED(objv))
{
    isst->camera.pos[2] += 0.05;
    isst->camera.focus[2] += 0.05;
    isst->dirty = 1;
    return TCL_OK;
}


static int
aetolookat(ClientData UNUSED(clientData), Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
{
    vect_t vecdfoc;
    double x, y;
    double az, el;
    double mag_vec;

    if (objc < 4) {
	Tcl_WrongNumArgs(interp, 1, objv, "pathName az el");
	return TCL_ERROR;
    }

    if (Tcl_GetDoubleFromObj(interp, objv[2], &x) != TCL_OK)
	return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[3], &y) != TCL_OK)
	return TCL_ERROR;

    mag_vec = DIST_PNT_PNT(isst->camera.pos, isst->camera.focus);

    VSUB2(vecdfoc, isst->camera.pos, isst->camera.focus);
    VUNITIZE(vecdfoc);
    bn_ae_vec(&az, &el, vecdfoc);
    az = az * -DEG2RAD + x;
    el = el * -DEG2RAD + y;
    bn_vec_ae(vecdfoc, az, el);
    VUNITIZE(vecdfoc);
    VSCALE(vecdfoc, vecdfoc, mag_vec);
    VADD2(isst->camera.focus, isst->camera.pos, vecdfoc);

    isst->dirty = 1;
    return TCL_OK;
}

static int
aerotate(ClientData UNUSED(clientData), Tcl_Interp *interp, int objc, Tcl_Obj *const *objv)
{
    vect_t vec, vecdpos, vecdfoc;
    double x, y;
    double az, el;
    double mag_pos, mag_focus;
    struct bu_vls tclstr = BU_VLS_INIT_ZERO;

    if (objc < 4) {
	Tcl_WrongNumArgs(interp, 1, objv, "pathName x y");
	return TCL_ERROR;
    }

    if (Tcl_GetDoubleFromObj(interp, objv[2], &x) != TCL_OK)
	return TCL_ERROR;
    if (Tcl_GetDoubleFromObj(interp, objv[3], &y) != TCL_OK)
	return TCL_ERROR;

    mag_pos = DIST_PNT_PNT(isst->camera.pos, isst->camera_focus_init);

    mag_focus = DIST_PNT_PNT(isst->camera.focus, isst->camera_focus_init);

    VSUB2(vecdpos, isst->camera_focus_init, isst->camera.pos);
    VUNITIZE(vecdpos);
    bn_ae_vec(&az, &el, vecdpos);
    az = az * -DEG2RAD - x;
    el = el * -DEG2RAD + y;

    /* clamp to sane values */
    while (az > M_2PI) az -= M_2PI;
    while (az < 0) az += M_2PI;
    if (el>M_PI_2) el=M_PI_2 - 0.001;
    if (el<-M_PI_2) el=-M_PI_2 + 0.001;

    bn_vec_ae(vecdfoc, az, el);
    VSCALE(vecdpos, vecdpos, mag_pos);
    VADD2(isst->camera.pos, isst->camera_focus_init, vecdpos);
    if (mag_focus > 0) {
	VSUB2(vecdfoc, isst->camera_focus_init, isst->camera.focus);
	VUNITIZE(vecdfoc);
	bn_ae_vec(&az, &el, vecdfoc);
	az = az * -DEG2RAD - x;
	el = el * -DEG2RAD + y;

	/* clamp to sane values */
	while (az > M_2PI) az -= M_2PI;
	while (az < 0) az += M_2PI;
	if (el>M_PI_2) el=M_PI_2 - 0.001;
	if (el<-M_PI_2) el=-M_PI_2 + 0.001;

	bn_vec_ae(vecdfoc, az, el);
	VSCALE(vecdfoc, vecdfoc, mag_focus);
	VADD2(isst->camera.focus, isst->camera_focus_init, vecdfoc);
    }
    /* Update the tcl copies of the az/el vars */
    VSUB2(vec, isst->camera.focus, isst->camera.pos);
    VUNITIZE(vec);
    bn_ae_vec(&az, &el, vec);
    bu_vls_sprintf(&tclstr, "%f", az);
    Tcl_SetVar(interp, "az", bu_vls_addr(&tclstr), 0);
    bu_vls_sprintf(&tclstr, "%f", el);
    Tcl_SetVar(interp, "el", bu_vls_addr(&tclstr), 0);
    bu_vls_free(&tclstr);
    isst->dirty = 1;
    return TCL_OK;
}

/* this function needs to be Isst_Init() for fbsd and mac, but may need to be
 * Issttcltk_Init on other platforms (I'm looking at you, windows). Needs more
 * investigation.
 */
int
Isst_Init(Tcl_Interp *interp)
{
    if (Tcl_PkgProvide(interp, "isst", "0.1") != TCL_OK) {
	return TCL_ERROR;
    }

    /*
     * Initialize Tcl
     */
    if (Tcl_InitStubs(interp, "8.1", 0) == NULL) {
	return TCL_ERROR;
    }

    /*
     * Specify the C callback functions for widget creation, display,
     * and reshape.
     */
    Tcl_CreateObjCommand(interp, "isst_init", isst_init, NULL, NULL);
    Tcl_CreateObjCommand(interp, "isst_zap", isst_zap, NULL, NULL);
    Tcl_CreateObjCommand(interp, "refresh_isst", paint_window, NULL, NULL);
    Tcl_CreateObjCommand(interp, "reshape", reshape, NULL, NULL);
    Tcl_CreateObjCommand(interp, "isst_bind_photo", isst_bind_photo, NULL, NULL);
    Tcl_CreateObjCommand(interp, "load_g", isst_load_g, NULL, NULL);
    Tcl_CreateObjCommand(interp, "list_g", list_geometry, NULL, NULL);
    Tcl_CreateObjCommand(interp, "aetolookat", aetolookat, NULL, NULL);
    Tcl_CreateObjCommand(interp, "aerotate", aerotate, NULL, NULL);
    Tcl_CreateObjCommand(interp, "walk", move_walk, NULL, NULL);
    Tcl_CreateObjCommand(interp, "strafe", move_strafe, NULL, NULL);
    Tcl_CreateObjCommand(interp, "float", move_float, NULL, NULL);
    Tcl_CreateObjCommand(interp, "reset", zero_view, NULL, NULL);
    Tcl_CreateObjCommand(interp, "set_resolution", set_resolution, NULL, NULL);
    Tcl_CreateObjCommand(interp, "render_mode", render_mode, NULL, NULL);
    return TCL_OK;
}

#ifdef HAVE_WINDOWS_H
int APIENTRY
WinMain(HINSTANCE hInstance,
	HINSTANCE hPrevInstance,
	LPSTR lpszCmdLine,
	int nCmdShow)
{
    char **argv;
    int argc;
#else
int
main(int argc, const char **argv)
{
#endif
    Tcl_DString temp;
const char *fullname;
    int status;
    const char *isst_tcl = NULL;
    struct bu_vls tlog = BU_VLS_INIT_ZERO;
    Tcl_Interp *interp = Tcl_CreateInterp();

#ifdef HAVE_WINDOWS_H
    /* Get our args from the c-runtime. Ignore lpszCmdLine. */
    argc = __argc;
    argv = __argv;
#endif

    /* initialize progname for run-tim resource finding */
    bu_setprogname(argv[0]);

#ifdef HAVE_WINDOWS_H
   Tk_InitConsoleChannels(interp);
#endif

    status = tclcad_init(interp, 1, &tlog);
    if (status == TCL_ERROR) {
	bu_log("Isst tclcad init failure:\n%s\n", bu_vls_addr(&tlog));
	bu_vls_free(&tlog);
	bu_exit(1, "tcl init error");
    }

    status = Isst_Init(interp);
    if (status == TCL_ERROR) {
	bu_vls_free(&tlog);
	bu_exit(1, "Isst Tcl/Tk init error");
    }
    bu_vls_free(&tlog);

    /* Skip first arg */
    argv++; argc--;
    tclcad_set_argv(interp, argc, argv);

    Tcl_DStringInit(&temp);
    isst_tcl = bu_dir(NULL, 0, BU_DIR_DATA, "tclscripts", "isst", "isst.tcl", NULL);
    fullname = Tcl_TranslateFileName(interp, isst_tcl, &temp);
    status = Tcl_EvalFile(interp, fullname);
    Tcl_DStringFree(&temp);
    Tcl_DeleteInterp(interp);
    return status;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
