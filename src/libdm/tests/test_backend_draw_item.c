/*        T E S T _ B A C K E N D _ D R A W _ I T E M . C
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
/** @file libdm/tests/test_backend_draw_item.c
 *
 * Regression coverage for the modern libdm backend draw contract.
 */

#include "common.h"

#include <stdio.h>
#include <string.h>

#include "bu/app.h"
#include "bu/log.h"
#include "bu/vls.h"
#include "dm.h"
#include "rt/view_legacy_bsg.h"
#include "../include/private.h"

static int g_fail = 0;
static int g_draw_item_count = 0;
static int g_resource_free_count = 0;

static int file_has_token(const char *path, const char *token);

#define DMCHECK(_cond, _msg) \
    do { \
	if (!(_cond)) { \
	    bu_log("FAIL [%s:%d] %s\n", __FILE__, __LINE__, (_msg)); \
	    g_fail++; \
	} \
    } while (0)

struct draw_capture {
    int line2d_count;
    int line3d_count;
    int string2d_count;
    point_t first_3d_a;
    point_t first_3d_b;
    char first_string[64];
    fastf_t first_string_x;
    fastf_t first_string_y;
    int first_string_size;
    int first_string_use_aspect;
    int string2d_rot_count;
    fastf_t first_string_rotation;
};

static struct draw_capture g_capture;

static int
capture_draw_line_2d(struct dm *UNUSED(dmp), fastf_t UNUSED(x1),
	fastf_t UNUSED(y1), fastf_t UNUSED(x2), fastf_t UNUSED(y2))
{
    g_capture.line2d_count++;
    return 0;
}

static int
capture_draw_line_3d(struct dm *UNUSED(dmp), point_t pt1, point_t pt2)
{
    if (!g_capture.line3d_count) {
	VMOVE(g_capture.first_3d_a, pt1);
	VMOVE(g_capture.first_3d_b, pt2);
    }
    g_capture.line3d_count++;
    return 0;
}

static int
capture_draw_string_2d(struct dm *UNUSED(dmp), const char *str,
	fastf_t x, fastf_t y, int size, int use_aspect)
{
    if (!g_capture.string2d_count) {
	snprintf(g_capture.first_string, sizeof(g_capture.first_string), "%s",
		str ? str : "");
	g_capture.first_string_x = x;
	g_capture.first_string_y = y;
	g_capture.first_string_size = size;
	g_capture.first_string_use_aspect = use_aspect;
    }
    g_capture.string2d_count++;
    return 0;
}

static int
capture_draw_string_2d_rot(struct dm *UNUSED(dmp), const char *str,
	fastf_t x, fastf_t y, int size, int use_aspect, fastf_t angle)
{
    if (!g_capture.string2d_count) {
	snprintf(g_capture.first_string, sizeof(g_capture.first_string), "%s",
		str ? str : "");
	g_capture.first_string_x = x;
	g_capture.first_string_y = y;
	g_capture.first_string_size = size;
	g_capture.first_string_use_aspect = use_aspect;
	g_capture.first_string_rotation = angle;
    }
    g_capture.string2d_count++;
    g_capture.string2d_rot_count++;
    return 0;
}

static struct dm *
open_null_dm(void)
{
    const char *av0 = "attach";
    return dm_open(NULL, NULL, "nu", 1, &av0);
}

static int
count_draw_item(struct dm *UNUSED(dmp), const void *render_item_ctx)
{
    if (!render_item_ctx)
	return -1;
    g_draw_item_count++;
    return 0;
}

static void
test_requires_draw_item(void)
{
    bu_log("=== backend draw item: no node-only fallback ===\n");

    struct dm *dmp = open_null_dm();
    DMCHECK(dmp != NULL, "opened null display manager");
    if (!dmp)
	return;

    void *v = rt_view_context_create_bsg();
    DMCHECK(v != NULL, "created retained test view context");
    if (!v) {
	dm_close(dmp);
	return;
    }

    void *item = rt_view_line_render_item_create_bsg(v, 42, "node-free-item",
	    84, 1);
    DMCHECK(item != NULL, "created opaque retained render item");
    if (!item) {
	rt_view_context_free_bsg(v);
	dm_close(dmp);
	return;
    }

    dm_set_backend_ops(dmp, NULL);
    int ret = dm_backend_draw_item(dmp, item);
    DMCHECK(ret != 0, "draw_item without backend op must fail");
    DMCHECK(g_draw_item_count == 0, "missing draw_item must not draw");

    static const struct dm_backend_ops count_ops = {
	0,
	NULL,
	count_draw_item,
	NULL,
	NULL,
	NULL
    };
    dm_set_backend_ops(dmp, &count_ops);
    ret = dm_backend_draw_item(dmp, item);
    DMCHECK(ret == 0, "registered draw_item op is called");
    DMCHECK(g_draw_item_count == 1, "draw_item count incremented once");

    dm_set_backend_ops(dmp, NULL);
    rt_view_render_item_free_bsg(item);
    rt_view_context_free_bsg(v);
    dm_close(dmp);
}

static void
count_resource_free(struct dm *UNUSED(dmp), struct dm_backend_resource *r)
{
    g_resource_free_count++;
    r->handle = NULL;
}

static void
test_backend_resource_cache_contract(void)
{
    bu_log("=== backend resource cache: render/backend record keys ===\n");

    struct dm *dmp = open_null_dm();
    DMCHECK(dmp != NULL, "opened null display manager");
    if (!dmp)
	return;

    g_resource_free_count = 0;

    struct dm_backend_resource *r1 = dm_backend_resource_get(dmp,
	    0x1001, 0x501, 0x9001, 1, count_resource_free);
    DMCHECK(r1 != NULL, "created first backend resource");
    r1->handle = (void *)0x1;

    struct dm_backend_resource *r2 = dm_backend_resource_get(dmp,
	    0x1001, 0x501, 0x9001, 0, count_resource_free);
    DMCHECK(r2 == r1, "same cache/source/backend identities reuse resource");

    struct dm_backend_resource *r3 = dm_backend_resource_get(dmp,
	    0x1001, 0x501, 0x9002, 1, count_resource_free);
    DMCHECK(r3 != NULL && r3 != r1, "backend capability/context identity separates resources");
    r3->handle = (void *)0x1;

    struct dm_backend_resource *r4 = dm_backend_resource_get(dmp,
	    0x1001, 0x502, 0x9001, 1, count_resource_free);
    DMCHECK(r4 != NULL && r4 != r1, "source identity separates resources");
    r4->handle = (void *)0x1;

    dm_backend_resource_invalidate(dmp, 0x501);
    DMCHECK(r1->stale == 1, "invalidate marks source resources stale");
    DMCHECK(r3->stale == 1, "invalidate marks same-source backend variants stale");
    DMCHECK(r4->stale == 0, "invalidate does not stale unrelated source");

    struct dm_backend_resource *stale = dm_backend_resource_get(dmp,
	    0x1001, 0x501, 0x9001, 0, count_resource_free);
    DMCHECK(stale == NULL, "stale resource is not returned without create");

    struct dm_backend_resource *r5 = dm_backend_resource_get(dmp,
	    0x1001, 0x501, 0x9001, 1, count_resource_free);
    DMCHECK(r5 != NULL, "create after stale invalidation creates replacement");
    DMCHECK(r5->stale == 0, "replacement resource is current");
    DMCHECK(g_resource_free_count >= 1, "stale resource was released on replacement");
    r5->handle = (void *)0x1;

    dm_backend_resource_release_source(dmp, 0x501);
    DMCHECK(g_resource_free_count >= 3, "release_source frees all source variants");
    dm_backend_resource_release_source(dmp, 0x502);

    dm_close(dmp);
}

static void
test_annotation_curve_dm_output(void)
{
    bu_log("=== backend draw item: annotation curve DM output ===\n");

    struct dm *dmp = open_null_dm();
    DMCHECK(dmp != NULL, "opened null display manager");
    if (!dmp)
	return;

    int (*save_line2d)(struct dm *, fastf_t, fastf_t, fastf_t, fastf_t) =
	dmp->i->dm_drawLine2D;
    int (*save_line3d)(struct dm *, point_t, point_t) =
	dmp->i->dm_drawLine3D;
    int (*save_string2d)(struct dm *, const char *, fastf_t, fastf_t, int, int) =
	dmp->i->dm_drawString2D;
    dmp->i->dm_drawLine2D = capture_draw_line_2d;
    dmp->i->dm_drawLine3D = capture_draw_line_3d;
    dmp->i->dm_drawString2D = capture_draw_string_2d;

    dm_set_width(dmp, 512);
    dm_set_height(dmp, 512);
    dmp->i->dm_aspect = 1.0;

    void *v = rt_view_context_create_bsg();
    DMCHECK(v != NULL, "created retained annotation view context");
    if (!v) {
	dmp->i->dm_drawLine2D = save_line2d;
	dmp->i->dm_drawLine3D = save_line3d;
	dmp->i->dm_drawString2D = save_string2d;
	dm_close(dmp);
	return;
    }

    DMCHECK(rt_view_context_display_manager_set_bsg(v, dmp),
	    "attached display manager to retained annotation view context");
    DMCHECK(rt_view_context_dimensions_set_bsg(v, 512, 512),
	    "set retained annotation view dimensions");
    DMCHECK(rt_view_context_model_matrices_identity_bsg(v),
	    "set retained annotation view matrices");

    point_t pts[2] = {VINIT_ZERO, VINIT_ZERO};
    VSET(pts[0], 0.0, 0.0, 0.0);
    VSET(pts[1], 1.0, 0.0, 0.0);
    DMCHECK(rt_view_annotation_curves_add_bsg(v, "dm_annotation_curves"),
	    "configured model-space annotation curves");

    memset(&g_capture, 0, sizeof(g_capture));
    dm_draw_objs(v);

    DMCHECK(g_capture.line3d_count >= 60,
	    "annotation CARC/NURB/BEZIER emitted 3D line output");
    DMCHECK(VNEAR_EQUAL(g_capture.first_3d_a, pts[0], SMALL_FASTF) &&
	    VNEAR_EQUAL(g_capture.first_3d_b, pts[1], SMALL_FASTF),
	    "annotation line segment preserved model-space endpoints");

    rt_view_context_free_bsg(v);
    dmp->i->dm_drawLine2D = save_line2d;
    dmp->i->dm_drawLine3D = save_line3d;
    dmp->i->dm_drawString2D = save_string2d;
    dm_close(dmp);
}

static void
test_annotation_display_text_position(void)
{
    bu_log("=== backend draw item: annotation display text position ===\n");

    struct dm *dmp = open_null_dm();
    DMCHECK(dmp != NULL, "opened null display manager");
    if (!dmp)
	return;

    int (*save_line2d)(struct dm *, fastf_t, fastf_t, fastf_t, fastf_t) =
	dmp->i->dm_drawLine2D;
    int (*save_line3d)(struct dm *, point_t, point_t) =
	dmp->i->dm_drawLine3D;
    int (*save_string2d)(struct dm *, const char *, fastf_t, fastf_t, int, int) =
	dmp->i->dm_drawString2D;
    int (*save_string2d_rot)(struct dm *, const char *, fastf_t, fastf_t, int, int, fastf_t) =
	dmp->i->dm_drawString2DRot;
    dmp->i->dm_drawLine2D = capture_draw_line_2d;
    dmp->i->dm_drawLine3D = capture_draw_line_3d;
    dmp->i->dm_drawString2D = capture_draw_string_2d;
    dmp->i->dm_drawString2DRot = capture_draw_string_2d_rot;

    dm_set_width(dmp, 200);
    dm_set_height(dmp, 200);
    dmp->i->dm_aspect = 1.0;

    void *v = rt_view_context_create_bsg();
    DMCHECK(v != NULL, "created retained text annotation view context");
    if (!v) {
	dmp->i->dm_drawLine2D = save_line2d;
	dmp->i->dm_drawLine3D = save_line3d;
	dmp->i->dm_drawString2D = save_string2d;
	dmp->i->dm_drawString2DRot = save_string2d_rot;
	dm_close(dmp);
	return;
    }

    DMCHECK(rt_view_context_display_manager_set_bsg(v, dmp),
	    "attached display manager to retained text annotation view context");
    DMCHECK(rt_view_context_dimensions_set_bsg(v, 200, 200),
	    "set retained text annotation view dimensions");
    DMCHECK(rt_view_context_model_matrices_identity_bsg(v),
	    "set retained text annotation view matrices");
    DMCHECK(rt_view_annotation_display_text_add_bsg(v, "dm_annotation_text",
	    "ABCD", 20.0, 45.0),
	    "configured display-space annotation text");

    memset(&g_capture, 0, sizeof(g_capture));
    dm_draw_objs(v);

    fastf_t expected_width = 2.0 * 20.0 * 0.6 * 4.0 / 200.0;
    fastf_t expected_height = 2.0 * 20.0 / 200.0;
    DMCHECK(g_capture.string2d_count == 1,
	    "annotation text emitted one native string draw");
    DMCHECK(strcmp(g_capture.first_string, "ABCD") == 0,
	    "annotation text preserves string");
    DMCHECK(NEAR_EQUAL(g_capture.first_string_x, -expected_width, 1.0e-6) &&
	    NEAR_EQUAL(g_capture.first_string_y, -expected_height, 1.0e-6),
	    "annotation top-right relative position offsets text from anchor");
    DMCHECK(g_capture.first_string_size == 20 &&
	    g_capture.first_string_use_aspect == 1,
	    "annotation text forwards size and aspect intent");
    DMCHECK(g_capture.string2d_rot_count == 1 &&
	    NEAR_EQUAL(g_capture.first_string_rotation, 45.0, SMALL_FASTF),
	    "annotation text forwards rotation intent");
    DMCHECK(g_capture.line2d_count == 0 && g_capture.line3d_count == 0,
	    "annotation text draw does not emit line fallback geometry");

    rt_view_context_free_bsg(v);
    dmp->i->dm_drawLine2D = save_line2d;
    dmp->i->dm_drawLine3D = save_line3d;
    dmp->i->dm_drawString2D = save_string2d;
    dmp->i->dm_drawString2DRot = save_string2d_rot;
    dm_close(dmp);
}

static void
test_generic_rotated_text_stroke_fallback(void)
{
    bu_log("=== generic DM rotated text stroke fallback ===\n");

    struct dm *dmp = open_null_dm();
    DMCHECK(dmp != NULL, "opened null display manager");
    if (!dmp)
	return;

    int (*save_line2d)(struct dm *, fastf_t, fastf_t, fastf_t, fastf_t) =
	dmp->i->dm_drawLine2D;
    int (*save_string2d)(struct dm *, const char *, fastf_t, fastf_t, int, int) =
	dmp->i->dm_drawString2D;
    int (*save_string2d_rot)(struct dm *, const char *, fastf_t, fastf_t, int, int, fastf_t) =
	dmp->i->dm_drawString2DRot;
    dmp->i->dm_drawLine2D = capture_draw_line_2d;
    dmp->i->dm_drawString2D = capture_draw_string_2d;
    dmp->i->dm_drawString2DRot = NULL;

    dm_set_width(dmp, 320);
    dm_set_height(dmp, 240);
    dmp->i->dm_aspect = 320.0 / 240.0;
    dm_set_fontsize(dmp, 12);

    memset(&g_capture, 0, sizeof(g_capture));
    int ret = dm_draw_string_2d_rot(dmp, "A", 0.0, 0.0, 20, 1, 30.0);
    DMCHECK(ret == BRLCAD_OK, "stroke fallback reports success");
    DMCHECK(g_capture.line2d_count > 0,
	    "rotated fallback emits vector line strokes");
    DMCHECK(g_capture.string2d_count == 0,
	    "rotated fallback does not emit unrotated native text");

    memset(&g_capture, 0, sizeof(g_capture));
    ret = dm_draw_string_2d_rot(dmp, "A", 0.0, 0.0, 20, 1, 0.0);
    DMCHECK(ret == BRLCAD_OK, "zero-degree fallback reports success");
    DMCHECK(g_capture.string2d_count == 1 && g_capture.line2d_count == 0,
	    "zero-degree text keeps native backend string path");

    dmp->i->dm_drawLine2D = save_line2d;
    dmp->i->dm_drawString2D = save_string2d;
    dmp->i->dm_drawString2DRot = save_string2d_rot;
    dm_close(dmp);
}

static void
test_postscript_rotated_text_output(const char *source_root)
{
    if (!source_root)
	return;

    bu_log("=== postscript backend: rotated native text output ===\n");

    struct bu_vls path = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&path, "%s/src/libdm/postscript/dm-ps.c", source_root);

    DMCHECK(file_has_token(bu_vls_cstr(&path), "ps_drawString2DRot") == 1,
	    "PostScript backend defines native rotated string draw");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "translate %g rotate") == 1,
	    "PostScript backend writes native rotate transform");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "ps_drawString2DRot,") == 1,
	    "PostScript backend registers rotated string draw hook");

    bu_vls_free(&path);
}

static void
test_plot_rotated_text_output(const char *source_root)
{
    if (!source_root)
	return;

    bu_log("=== plot backend: rotated vector text output ===\n");

    struct bu_vls path = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&path, "%s/src/libdm/plot/dm-plot.c", source_root);

    DMCHECK(file_has_token(bu_vls_cstr(&path), "plot_drawString2DRot") == 1,
	    "Plot backend defines native rotated string draw");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "pd_symbol") == 1,
	    "Plot backend writes rotated vector text");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "plot_drawString2DRot,") == 1,
	    "Plot backend registers rotated string draw hook");

    bu_vls_free(&path);
}

static int
file_has_token(const char *path, const char *token)
{
    FILE *fp = fopen(path, "rb");
    if (!fp)
	return -1;

    int found = 0;
    char buf[4096];
    while (bu_fgets(buf, sizeof(buf), fp)) {
	if (strstr(buf, token)) {
	    found = 1;
	    break;
	}
    }
    fclose(fp);
    return found;
}

static void
test_private_dm_vlist_hygiene(const char *source_root)
{
    if (!source_root)
	return;

    bu_log("=== private dm direct-vlist hygiene ===\n");

    const char *files[] = {
	"src/libdm/include/private.h",
	"src/libdm/include/calltable.h",
	"src/libdm/dm-generic.c",
	"src/libdm/dm-gl.h",
	"src/libdm/dm-gl.c",
	"src/libdm/X/dm-X.c",
	"src/libdm/plot/dm-plot.c",
	"src/libdm/postscript/dm-ps.c",
	"src/libdm/null/dm-Null.h",
	"src/libdm/null/dm-Null.c",
	"src/libdm/txt/dm-txt.c"
    };
    const char *forbidden[] = {
	"bsg/vlist.h",
	"bsg_vlist",
	"BSG_VLIST"
    };

    struct bu_vls path = BU_VLS_INIT_ZERO;
    for (size_t i = 0; i < sizeof(files) / sizeof(files[0]); i++) {
	bu_vls_sprintf(&path, "%s/%s", source_root, files[i]);
	DMCHECK(file_has_token(bu_vls_cstr(&path), "bg_vlist") == 1,
		"private dm direct-vlist source uses neutral bg_vlist spelling");
	for (size_t j = 0; j < sizeof(forbidden) / sizeof(forbidden[0]); j++) {
	    int ret = file_has_token(bu_vls_cstr(&path), forbidden[j]);
	    if (ret < 0) {
		DMCHECK(0, "could not open private dm direct-vlist source");
		break;
	    }
	    if (ret) {
		struct bu_vls msg = BU_VLS_INIT_ZERO;
		bu_vls_sprintf(&msg, "private dm direct-vlist source contains legacy token: %s", forbidden[j]);
		DMCHECK(0, bu_vls_cstr(&msg));
		bu_vls_free(&msg);
	    }
	}
    }

    bu_vls_sprintf(&path, "%s/src/libdm/plot/dm-plot.c", source_root);
    DMCHECK(file_has_token(bu_vls_cstr(&path), "rt_vlist_to_uplot") == 1,
	    "plot backend routes floating vlist export through the RT vlist bridge");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "bsg_vlist_to_uplot") == 0,
	    "plot backend does not call the BSG vlist uplot wrapper");

    bu_vls_free(&path);
}

static void
test_retained_gl_source_hygiene(const char *source_root)
{
    if (!source_root)
	return;

    bu_log("=== retained GL source hygiene ===\n");

    struct bu_vls path = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&path, "%s/src/libdm/dm-gl_lod.cpp", source_root);

    const char *needles[] = {
	"dm_draw_device_vlist(",
	"dm_draw_device_vlist_hidden_line(",
	"dlist",
	"display list",
	"fallback"
    };
    for (size_t i = 0; i < sizeof(needles) / sizeof(needles[0]); i++) {
	int ret = file_has_token(bu_vls_cstr(&path), needles[i]);
	if (ret < 0) {
	    DMCHECK(0, "could not open retained GL source file");
	    break;
	}
	if (ret) {
	    struct bu_vls msg = BU_VLS_INIT_ZERO;
	    bu_vls_sprintf(&msg, "retained GL source contains forbidden token: %s", needles[i]);
	    DMCHECK(0, bu_vls_cstr(&msg));
	    bu_vls_free(&msg);
	}
    }

    bu_vls_free(&path);
}

static void
test_public_dm_header_hygiene(const char *source_root)
{
    if (!source_root)
	return;

    bu_log("=== public dm.h vlist hygiene ===\n");

    struct bu_vls path = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&path, "%s/include/dm.h", source_root);

    const char *needles[] = {
	"bsg/vlist.h",
	"bsg_vlist",
	"drawVList",
	"dm_draw("
    };
    for (size_t i = 0; i < sizeof(needles) / sizeof(needles[0]); i++) {
	int ret = file_has_token(bu_vls_cstr(&path), needles[i]);
	if (ret < 0) {
	    DMCHECK(0, "could not open public dm.h");
	    break;
	}
	if (ret) {
	    struct bu_vls msg = BU_VLS_INIT_ZERO;
	    bu_vls_sprintf(&msg, "public dm.h contains legacy vlist token: %s", needles[i]);
	    DMCHECK(0, bu_vls_cstr(&msg));
	    bu_vls_free(&msg);
	}
    }

    bu_vls_sprintf(&path, "%s/include/dm/vlist.h", source_root);
    DMCHECK(file_has_token(bu_vls_cstr(&path), "dm_draw(") == 1,
	    "legacy dm_draw prototype is isolated in dm/vlist.h");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "bg/vlist.h") == 1,
	    "legacy dm_draw prototype uses the neutral BG vlist spelling");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "bsg/vlist.h") == 0,
	    "legacy dm_draw prototype does not include the BSG vlist compatibility header");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "bsg_vlist") == 0,
	    "legacy dm_draw prototype does not expose the BSG vlist compatibility typedef");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "drawVList") == 0,
	    "direct drawVList hooks remain backend-local and out of dm/vlist.h");

    bu_vls_sprintf(&path, "%s/include/dm/view.h", source_root);
    DMCHECK(file_has_token(bu_vls_cstr(&path), "dm_draw_faceplate(void *") == 1,
	    "legacy faceplate draw prototype uses an opaque view context");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "dm_draw_objs(void *") == 1,
	    "legacy draw prototype uses an opaque view context");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "struct bsg_view") == 0,
	    "legacy draw prototypes do not expose BSG view spelling");
    DMCHECK(file_has_token(bu_vls_cstr(&path), "bsg/") == 0,
	    "legacy draw prototypes do not include BSG headers");

    bu_vls_free(&path);
}

int
main(int argc, const char **argv)
{
    bu_setprogname(argv[0]);

    test_requires_draw_item();
    test_backend_resource_cache_contract();
    test_annotation_curve_dm_output();
    test_annotation_display_text_position();
    test_generic_rotated_text_stroke_fallback();
    test_postscript_rotated_text_output(argc > 1 ? argv[1] : NULL);
    test_plot_rotated_text_output(argc > 1 ? argv[1] : NULL);
    test_private_dm_vlist_hygiene(argc > 1 ? argv[1] : NULL);
    test_retained_gl_source_hygiene(argc > 1 ? argv[1] : NULL);
    test_public_dm_header_hygiene(argc > 1 ? argv[1] : NULL);

    if (g_fail) {
	bu_log("%d backend draw item test failure(s)\n", g_fail);
	return 1;
    }

    bu_log("backend draw item tests passed\n");
    return 0;
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
