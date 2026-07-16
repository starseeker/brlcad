/*                            D M . C
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
/** @file libged/dm.c
 *
 * The dm and screengrab commands.
 *
 */

#include "common.h"

#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "brlobol/display_endpoint.h"
#include "brlobol/host_factory.h"
#include "bu/cmd.h"
#include "bu/log.h"
#include "bu/opt.h"
#include "bu/vls.h"
#include "bv.h"
#include "ged/draw_obol.h"

#include "../ged_private.h"

#define HELPFLAG "--print-help"
#define PURPOSEFLAG "--print-purpose"

struct _ged_dm_info {
    struct ged *gedp;
    long verbosity;
    const struct bu_cmdtab *cmds;
    struct bu_opt_desc *gopts;
};

static const char *
_dm_view_name(const void *view_ctx)
{
    const struct bv *view = bv_context_view_const((const struct bv_context *)view_ctx);
    const char *name = bv_name_get(view);
    return name ? name : "";
}

static void
_dm_cmd_during_clbk(struct _ged_dm_info *gd, int argc, const char **argv)
{
    if (!gd || !gd->gedp || argc < 1 || !argv)
        return;

    bu_clbk_t clbk = NULL;
    void *u2 = NULL;
    if (ged_clbk_get(&clbk, &u2, gd->gedp, "dm", BU_CLBK_DURING) != BRLCAD_OK || !clbk)
        return;

    const char *dbg = getenv("GED_DISPLAY_DURING_DEBUG");
    if (dbg && BU_STR_EQUAL(dbg, "1")) {
        bu_log("ged dm during callback: ");
        for (int i = 0; i < argc; i++) {
            bu_log("%s%s", argv[i], (i + 1 < argc) ? " " : "\n");
        }
    }

    int cbret = ged_clbk_exec(gd->gedp->ged_result_str, gd->gedp, GED_CMD_RECURSION_LIMIT, clbk, argc, argv, (void *)gd->gedp, u2);
    if (cbret != BRLCAD_OK) {
        bu_vls_printf(gd->gedp->ged_result_str, "\nwarning: dm during callback returned %d", cbret);
    }
}

static int
_dm_cmd_msgs(void *bs, int argc, const char **argv, const char *us, const char *ps)
{
    struct _ged_dm_info *gd = (struct _ged_dm_info *)bs;
    if (argc == 2 && BU_STR_EQUAL(argv[1], HELPFLAG)) {
	bu_vls_printf(gd->gedp->ged_result_str, "%s\n%s\n", us, ps);
	return 1;
    }
    if (argc == 2 && BU_STR_EQUAL(argv[1], PURPOSEFLAG)) {
	bu_vls_printf(gd->gedp->ged_result_str, "%s\n", ps);
	return 1;
    }
    return 0;
}

static const char *
_dm_render_engine_name(enum brlobol_render_engine engine)
{
    switch (engine) {
	case BRLOBOL_RENDER_ENGINE_HW: return "hw";
	case BRLOBOL_RENDER_ENGINE_SW: return "sw";
	case BRLOBOL_RENDER_ENGINE_RT: return "rt";
	case BRLOBOL_RENDER_ENGINE_NONE: return "none";
	case BRLOBOL_RENDER_ENGINE_DIAGNOSTIC: return "diagnostic";
	default: return "auto";
    }
}

static int
_dm_render_engine_parse(const char *name, enum brlobol_render_engine *engine)
{
    if (!name || !engine)
	return 0;
    if (BU_STR_EQUAL(name, "auto"))
	*engine = BRLOBOL_RENDER_ENGINE_AUTO;
    else if (BU_STR_EQUAL(name, "hw"))
	*engine = BRLOBOL_RENDER_ENGINE_HW;
    else if (BU_STR_EQUAL(name, "sw"))
	*engine = BRLOBOL_RENDER_ENGINE_SW;
    else if (BU_STR_EQUAL(name, "rt"))
	*engine = BRLOBOL_RENDER_ENGINE_RT;
    else if (BU_STR_EQUAL(name, "none"))
	*engine = BRLOBOL_RENDER_ENGINE_NONE;
    else if (BU_STR_EQUAL(name, "diagnostic"))
	*engine = BRLOBOL_RENDER_ENGINE_DIAGNOSTIC;
    else
	return 0;
    return 1;
}

static void *
_dm_endpoint_view(struct _ged_dm_info *gd, const struct bu_vls *view_name)
{
    if (!gd || !gd->gedp)
	return NULL;
    if (view_name && bu_vls_strlen(view_name))
	return ged_view_find_ctx(gd->gedp, bu_vls_cstr(view_name));
    return ged_view_active_ctx(gd->gedp);
}

static brlobol_display_endpoint_t *
_dm_endpoint(struct _ged_dm_info *gd, const struct bu_vls *view_name,
	void **view_ctx)
{
    void *view = _dm_endpoint_view(gd, view_name);
    if (view_ctx)
	*view_ctx = view;
    if (!view) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"target view does not exist\n");
	return NULL;
    }
    brlobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(view);
    if (!endpoint) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"view '%s' has no Obol display endpoint\n", _dm_view_name(view));
	return NULL;
    }
    return endpoint;
}

static void
_dm_property_print(struct bu_vls *out,
	const struct brlobol_endpoint_property_value *value)
{
    switch (value->type) {
	case BRLOBOL_ENDPOINT_PROPERTY_BOOL:
	    bu_vls_printf(out, "%s", value->bool_value ? "true" : "false");
	    break;
	case BRLOBOL_ENDPOINT_PROPERTY_INT:
	    bu_vls_printf(out, "%" PRId64 "", value->int_value);
	    break;
	case BRLOBOL_ENDPOINT_PROPERTY_UINT:
	    bu_vls_printf(out, "%" PRIu64 "", value->uint_value);
	    break;
	case BRLOBOL_ENDPOINT_PROPERTY_DOUBLE:
	    bu_vls_printf(out, "%.17g", value->double_value);
	    break;
	case BRLOBOL_ENDPOINT_PROPERTY_COLOR3:
	    bu_vls_printf(out, "%.17g/%.17g/%.17g", value->color3[0],
		    value->color3[1], value->color3[2]);
	    break;
	case BRLOBOL_ENDPOINT_PROPERTY_STRING:
	case BRLOBOL_ENDPOINT_PROPERTY_ENUM:
	    bu_vls_printf(out, "%s", value->string_value ? value->string_value : "");
	    break;
	default:
	    break;
    }
}

static int
_dm_property_value_parse(struct brlobol_endpoint_property_value *value,
	const char *str)
{
    char *end = NULL;
    if (!value || !str)
	return 0;
    errno = 0;
    switch (value->type) {
	case BRLOBOL_ENDPOINT_PROPERTY_BOOL:
	    if (BU_STR_EQUAL(str, "!"))
		value->bool_value = value->bool_value ? 0 : 1;
	    else if (BU_STR_EQUAL(str, "true") || BU_STR_EQUAL(str, "1") ||
		    BU_STR_EQUAL(str, "on"))
		value->bool_value = 1;
	    else if (BU_STR_EQUAL(str, "false") || BU_STR_EQUAL(str, "0") ||
		    BU_STR_EQUAL(str, "off"))
		value->bool_value = 0;
	    else
		return 0;
	    return 1;
	case BRLOBOL_ENDPOINT_PROPERTY_INT:
	    value->int_value = strtoimax(str, &end, 0);
	    return !errno && end && !*end;
	case BRLOBOL_ENDPOINT_PROPERTY_UINT:
	    if (*str == '-')
		return 0;
	    value->uint_value = strtoumax(str, &end, 0);
	    return !errno && end && !*end;
	case BRLOBOL_ENDPOINT_PROPERTY_DOUBLE:
	    value->double_value = strtod(str, &end);
	    return !errno && end && !*end && isfinite(value->double_value);
	case BRLOBOL_ENDPOINT_PROPERTY_COLOR3: {
	    double r = 0.0, g = 0.0, b = 0.0;
	    char trailing = '\0';
	    if (sscanf(str, "%lf/%lf/%lf%c", &r, &g, &b, &trailing) != 3 &&
		    sscanf(str, "%lf,%lf,%lf%c", &r, &g, &b, &trailing) != 3)
		return 0;
	    if (!isfinite(r) || !isfinite(g) || !isfinite(b))
		return 0;
	    value->color3[0] = r;
	    value->color3[1] = g;
	    value->color3[2] = b;
	    return 1;
	}
	case BRLOBOL_ENDPOINT_PROPERTY_STRING:
	case BRLOBOL_ENDPOINT_PROPERTY_ENUM:
	    value->string_value = str;
	    return 1;
	default:
	    return 0;
    }
}

static const char *
_dm_property_error(int ret)
{
    switch (ret) {
	case BRLOBOL_ENDPOINT_PROPERTY_UNKNOWN: return "unknown property";
	case BRLOBOL_ENDPOINT_PROPERTY_INVALID: return "invalid value";
	case BRLOBOL_ENDPOINT_PROPERTY_READ_ONLY: return "read-only property";
	case BRLOBOL_ENDPOINT_PROPERTY_UNSUPPORTED: return "unsupported by this endpoint";
	default: return "property operation failed";
    }
}

static const char *
_dm_property_type_name(enum brlobol_endpoint_property_type type)
{
    switch (type) {
	case BRLOBOL_ENDPOINT_PROPERTY_BOOL: return "bool";
	case BRLOBOL_ENDPOINT_PROPERTY_INT: return "int";
	case BRLOBOL_ENDPOINT_PROPERTY_UINT: return "uint";
	case BRLOBOL_ENDPOINT_PROPERTY_DOUBLE: return "double";
	case BRLOBOL_ENDPOINT_PROPERTY_STRING: return "string";
	case BRLOBOL_ENDPOINT_PROPERTY_COLOR3: return "color3";
	case BRLOBOL_ENDPOINT_PROPERTY_ENUM: return "enum";
	default: return "unknown";
    }
}

int
_dm_cmd_bg(void *ds, int argc, const char **argv)
{
    const char *usage_string =
	"dm [options] bg [-V view] [r/g/b [r/g/b]]";
    const char *purpose_string = "get or set the endpoint background gradient";
    if (_dm_cmd_msgs(ds, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    struct _ged_dm_info *gd = (struct _ged_dm_info *)ds;
    struct ged *gedp = gd->gedp;
    argc--; argv++;
    struct bu_vls view_name = BU_VLS_INIT_ZERO;
    struct bu_opt_desc d[2];
    BU_OPT(d[0], "V", "view", "name", &bu_opt_vls, &view_name,
	    "view endpoint to inspect or modify");
    BU_OPT_NULL(d[1]);
    int ac = bu_opt_parse(NULL, argc, argv, d);
    brlobol_display_endpoint_t *endpoint =
	_dm_endpoint(gd, &view_name, NULL);
    bu_vls_free(&view_name);
    if (!endpoint)
	return BRLCAD_ERROR;

    if (!ac) {
	struct brlobol_endpoint_property_value bottom =
	    BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	struct brlobol_endpoint_property_value top =
	    BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	if (brlobol_display_endpoint_property_get(endpoint,
		"controller.background.bottom", &bottom) !=
		BRLOBOL_ENDPOINT_PROPERTY_OK ||
	    brlobol_display_endpoint_property_get(endpoint,
		"controller.background.top", &top) !=
		BRLOBOL_ENDPOINT_PROPERTY_OK)
	    return BRLCAD_ERROR;
	unsigned char b[3] = {
	    (unsigned char)lrint(bottom.color3[0] * 255.0),
	    (unsigned char)lrint(bottom.color3[1] * 255.0),
	    (unsigned char)lrint(bottom.color3[2] * 255.0)
	};
	unsigned char t[3] = {
	    (unsigned char)lrint(top.color3[0] * 255.0),
	    (unsigned char)lrint(top.color3[1] * 255.0),
	    (unsigned char)lrint(top.color3[2] * 255.0)
	};
	if (memcmp(b, t, sizeof(b)) == 0)
	    bu_vls_printf(gedp->ged_result_str, "%u/%u/%u\n",
		    b[0], b[1], b[2]);
	else
	    bu_vls_printf(gedp->ged_result_str,
		    "%u/%u/%u->%u/%u/%u\n",
		    b[0], b[1], b[2], t[0], t[1], t[2]);
	return BRLCAD_OK;
    }

    struct bu_color bottom_color = BU_COLOR_INIT_ZERO;
    struct bu_color top_color = BU_COLOR_INIT_ZERO;
    int ac_used = bu_opt_color(NULL, ac, argv, &bottom_color);
    if (ac_used == -1) {
	bu_vls_printf(gedp->ged_result_str, "invalid color specification\n");
	return BRLCAD_ERROR;
    }
    int remaining = ac - ac_used;
    top_color = bottom_color;
    if (remaining) {
	int top_used = bu_opt_color(NULL, remaining, argv + ac_used,
		&top_color);
	if (top_used == -1 || top_used != remaining) {
	    bu_vls_printf(gedp->ged_result_str, "invalid color specification\n");
	    return BRLCAD_ERROR;
	}
    }

    fastf_t bottom_rgb[3] = {0.0, 0.0, 0.0};
    fastf_t top_rgb[3] = {0.0, 0.0, 0.0};
    if (!bu_color_to_rgb_floats(&bottom_color, bottom_rgb) ||
	!bu_color_to_rgb_floats(&top_color, top_rgb))
	return BRLCAD_ERROR;
    struct brlobol_endpoint_property_value value =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    value.type = BRLOBOL_ENDPOINT_PROPERTY_COLOR3;
    for (int i = 0; i < 3; i++)
	value.color3[i] = bottom_rgb[i];
    if (brlobol_display_endpoint_property_set(endpoint,
	    "controller.background.bottom", &value) !=
	    BRLOBOL_ENDPOINT_PROPERTY_OK)
	return BRLCAD_ERROR;
    for (int i = 0; i < 3; i++)
	value.color3[i] = top_rgb[i];
    if (brlobol_display_endpoint_property_set(endpoint,
	    "controller.background.top", &value) !=
	    BRLOBOL_ENDPOINT_PROPERTY_OK)
	return BRLCAD_ERROR;

    const char *cbav[2] = {"dm", "bg"};
    _dm_cmd_during_clbk(gd, 2, cbav);
    (void)brlobol_display_endpoint_request_frame(endpoint,
	    "dm background changed");

    return BRLCAD_OK;
}

int
_dm_cmd_list(void *ds, int argc, const char **argv)
{
    const char *usage_string = "dm [options] list [endpoints|hosts|renderers]";
    const char *purpose_string = "list Obol endpoints, hosts, or renderers.";
    if (_dm_cmd_msgs(ds, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--; argv++;

    struct _ged_dm_info *gd = (struct _ged_dm_info *)ds;
    struct ged *gedp = gd->gedp;

    if (argc > 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s", usage_string);
	return BRLCAD_ERROR;
    }
    if (argc == 1 && BU_STR_EQUAL(argv[0], "renderers")) {
	bu_vls_printf(gedp->ged_result_str,
		"auto\nhw\nsw\nrt\nnone\ndiagnostic\n");
	return BRLCAD_OK;
    }
    if (argc == 1 && BU_STR_EQUAL(argv[0], "hosts")) {
	const size_t count = brlobol_host_factory_registry_count();
	for (size_t i = 0; i < count; i++) {
	    char name[256] = {0};
	    if (!brlobol_host_factory_registry_name(i, name, sizeof(name)))
		continue;
	    if (gd->verbosity)
		bu_vls_printf(gedp->ged_result_str, "%s capabilities=0x%016"
			PRIx64 "\n", name,
			brlobol_host_factory_registry_capabilities(i));
	    else
		bu_vls_printf(gedp->ged_result_str, "%s\n", name);
	}
	return BRLCAD_OK;
    }
    if (argc == 0 || (argc == 1 && BU_STR_EQUAL(argv[0], "endpoints"))) {
	struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
	const size_t count = views ? BU_PTBL_LEN(views) : 0;
	for (size_t i = 0; i < count; i++) {
	    void *view = (void *)BU_PTBL_GET(views, i);
	    brlobol_display_endpoint_t *endpoint =
		ged_view_context_display_endpoint_get(view);
	    if (!endpoint)
		continue;
	    bu_vls_printf(gedp->ged_result_str, "%s%s\n",
		    view == ged_view_active_ctx(gedp) ? "*" : "",
		    _dm_view_name(view));
	}
	return BRLCAD_OK;
    }
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str,
		"unknown list target '%s'; expected endpoints, hosts, or renderers\n",
		argv[0]);
	return BRLCAD_ERROR;
    }
    return BRLCAD_ERROR;
}

int
_dm_cmd_status(void *ds, int argc, const char **argv)
{
    const char *usage_string = "dm [options] status [-V view]";
    const char *purpose_string = "report the selected Obol display endpoint.";
    if (_dm_cmd_msgs(ds, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;
    struct _ged_dm_info *gd = (struct _ged_dm_info *)ds;
    struct bu_vls view_name = BU_VLS_INIT_ZERO;
    struct bu_opt_desc d[2];
    BU_OPT(d[0], "V", "view", "name", &bu_opt_vls, &view_name,
	    "view endpoint to report");
    BU_OPT_NULL(d[1]);
    int ac = bu_opt_parse(NULL, argc, argv, d);
    if (ac) {
	bu_vls_printf(gd->gedp->ged_result_str, "Usage: %s", usage_string);
	bu_vls_free(&view_name);
	return BRLCAD_ERROR;
    }

    void *view = NULL;
    brlobol_display_endpoint_t *endpoint =
	_dm_endpoint(gd, &view_name, &view);
    bu_vls_free(&view_name);
    if (!endpoint)
	return BRLCAD_ERROR;

    struct brlobol_endpoint_property_value width =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    struct brlobol_endpoint_property_value height =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    struct brlobol_endpoint_property_value dpr =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    (void)brlobol_display_endpoint_property_get(endpoint, "endpoint.width",
	    &width);
    (void)brlobol_display_endpoint_property_get(endpoint, "endpoint.height",
	    &height);
    (void)brlobol_display_endpoint_property_get(endpoint,
	    "endpoint.device_pixel_ratio", &dpr);
    const char *host = brlobol_display_endpoint_host_factory_name(endpoint);
    bu_vls_printf(gd->gedp->ged_result_str,
	    "view=%s host=%s renderer=%s resolved=%s width=%" PRIu64
	    " height=%" PRIu64 " device_pixel_ratio=%.17g\n",
	    _dm_view_name(view), host ? host : "unbound",
	    _dm_render_engine_name(
		brlobol_display_endpoint_render_engine_get(endpoint)),
	    _dm_render_engine_name(
		brlobol_display_endpoint_render_engine_resolved_get(endpoint)),
	    width.uint_value, height.uint_value, dpr.double_value);
    return BRLCAD_OK;
}

int
_dm_cmd_renderer(void *ds, int argc, const char **argv)
{
    const char *usage_string =
	"dm [options] renderer [-V view] [auto|hw|sw|rt|none|diagnostic]";
    const char *purpose_string = "get or select the Obol endpoint renderer.";
    if (_dm_cmd_msgs(ds, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;
    struct _ged_dm_info *gd = (struct _ged_dm_info *)ds;
    struct bu_vls view_name = BU_VLS_INIT_ZERO;
    struct bu_opt_desc d[2];
    BU_OPT(d[0], "V", "view", "name", &bu_opt_vls, &view_name,
	    "view endpoint to inspect or modify");
    BU_OPT_NULL(d[1]);
    int ac = bu_opt_parse(NULL, argc, argv, d);
    if (ac > 1) {
	bu_vls_printf(gd->gedp->ged_result_str, "Usage: %s", usage_string);
	bu_vls_free(&view_name);
	return BRLCAD_ERROR;
    }

    brlobol_display_endpoint_t *endpoint =
	_dm_endpoint(gd, &view_name, NULL);
    bu_vls_free(&view_name);
    if (!endpoint)
	return BRLCAD_ERROR;
    if (!ac) {
	bu_vls_printf(gd->gedp->ged_result_str, "%s",
		_dm_render_engine_name(
		    brlobol_display_endpoint_render_engine_get(endpoint)));
	return BRLCAD_OK;
    }

    enum brlobol_render_engine engine = BRLOBOL_RENDER_ENGINE_AUTO;
    if (!_dm_render_engine_parse(argv[0], &engine)) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"unknown renderer '%s'; expected auto, hw, sw, rt, none, or diagnostic\n",
		argv[0]);
	return BRLCAD_ERROR;
    }
    if (!brlobol_display_endpoint_render_engine_set(endpoint, engine)) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"renderer '%s' is unsupported by this endpoint host\n", argv[0]);
	return BRLCAD_ERROR;
    }
    const char *cbav[3] = {"dm", "renderer", argv[0]};
    _dm_cmd_during_clbk(gd, 3, cbav);
    return BRLCAD_OK;
}

int
_dm_cmd_open(void *ds, int argc, const char **argv)
{
    const char *usage_string =
	"dm [options] open [-V view] --host name "
	"[--renderer auto|hw|sw|rt|none|diagnostic]";
    const char *purpose_string =
	"open a registered host for an Obol display endpoint.";
    if (_dm_cmd_msgs(ds, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;
    struct _ged_dm_info *gd = (struct _ged_dm_info *)ds;
    struct bu_vls view_name = BU_VLS_INIT_ZERO;
    struct bu_vls host_name = BU_VLS_INIT_ZERO;
    struct bu_vls renderer_name = BU_VLS_INIT_ZERO;
    struct bu_opt_desc d[4];
    BU_OPT(d[0], "V", "view", "name", &bu_opt_vls, &view_name,
	    "view endpoint to open");
    BU_OPT(d[1], "H", "host", "name", &bu_opt_vls, &host_name,
	    "registered host factory");
    BU_OPT(d[2], "R", "renderer", "name", &bu_opt_vls, &renderer_name,
	    "render engine");
    BU_OPT_NULL(d[3]);
    int ac = bu_opt_parse(NULL, argc, argv, d);
    if (ac || !bu_vls_strlen(&host_name)) {
	bu_vls_printf(gd->gedp->ged_result_str, "Usage: %s", usage_string);
	bu_vls_free(&view_name);
	bu_vls_free(&host_name);
	bu_vls_free(&renderer_name);
	return BRLCAD_ERROR;
    }

    void *view = _dm_endpoint_view(gd, &view_name);
    if (!view) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"target view does not exist\n");
	bu_vls_free(&view_name);
	bu_vls_free(&host_name);
	bu_vls_free(&renderer_name);
	return BRLCAD_ERROR;
    }
    brlobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(view);
    if (!endpoint &&
	ged_draw_obol_render_endpoint_ensure_for_view(gd->gedp, view, 1))
	endpoint = ged_view_context_display_endpoint_get(view);
    if (!endpoint) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"failed to create Obol display endpoint\n");
	bu_vls_free(&view_name);
	bu_vls_free(&host_name);
	bu_vls_free(&renderer_name);
	return BRLCAD_ERROR;
    }

    const enum brlobol_render_engine old_engine =
	brlobol_display_endpoint_render_engine_get(endpoint);
    enum brlobol_render_engine requested_engine = old_engine;
    if (bu_vls_strlen(&renderer_name) &&
	!_dm_render_engine_parse(bu_vls_cstr(&renderer_name),
	    &requested_engine)) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"unknown renderer '%s'; expected auto, hw, sw, rt, none, or diagnostic\n",
		bu_vls_cstr(&renderer_name));
	bu_vls_free(&view_name);
	bu_vls_free(&host_name);
	bu_vls_free(&renderer_name);
	return BRLCAD_ERROR;
    }
    struct brlobol_endpoint_property_value dpr_value =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    (void)brlobol_display_endpoint_property_get(endpoint,
	    "endpoint.device_pixel_ratio", &dpr_value);
    struct bv *bview = bv_context_view((struct bv_context *)view);
    unsigned int host_width = bview && bv_width_get(bview) > 0 ?
	(unsigned int)bv_width_get(bview) : 512u;
    unsigned int host_height = bview && bv_height_get(bview) > 0 ?
	(unsigned int)bv_height_get(bview) : 512u;
    if (bview && (bv_width_get(bview) <= 0 || bv_height_get(bview) <= 0))
	(void)bv_dimensions_set(bview, (int)host_width, (int)host_height);
    struct brlobol_host_desc desc = {0};
    desc.struct_size = sizeof(desc);
    desc.mode = BU_STR_EQUAL(bu_vls_cstr(&host_name), "headless") ?
	BRLOBOL_HOST_MODE_HEADLESS : BRLOBOL_HOST_MODE_TOPLEVEL;
    desc.width = host_width;
    desc.height = host_height;
    desc.device_pixel_ratio =
	dpr_value.type == BRLOBOL_ENDPOINT_PROPERTY_DOUBLE ?
	dpr_value.double_value : 1.0;
    desc.visible = desc.mode == BRLOBOL_HOST_MODE_TOPLEVEL ? 1 : 0;
    desc.title = _dm_view_name(view);
    /* Tk is the only current host needing a toolkit context.  The old
     * DM-type map had no writers, so its lookup could never contribute. */
    if (BU_STR_EQUIV(bu_vls_cstr(&host_name), "tk-gl") ||
	BU_STR_EQUIV(bu_vls_cstr(&host_name), "tk-photo"))
	desc.application_context = gd->gedp->ged_interp;
    if (requested_engine == BRLOBOL_RENDER_ENGINE_HW)
	desc.required_capabilities |= BRLOBOL_HOST_CAP_SYSTEM_GL;
    else if (requested_engine == BRLOBOL_RENDER_ENGINE_SW)
	desc.required_capabilities |= BRLOBOL_HOST_CAP_PIXEL_PRESENT;
    else if (requested_engine == BRLOBOL_RENDER_ENGINE_RT)
	desc.required_capabilities |= BRLOBOL_HOST_CAP_PROGRESSIVE_PRESENT;

    /* Host creation is staged by display_endpoint_host_open.  AUTO avoids
     * rejecting the old host before its replacement is ready; the descriptor
     * still constrains selection for an explicit requested engine. */
    if (!brlobol_display_endpoint_render_engine_set(endpoint,
	    BRLOBOL_RENDER_ENGINE_AUTO) ||
	!brlobol_display_endpoint_host_open(endpoint,
	    bu_vls_cstr(&host_name), &desc)) {
	(void)brlobol_display_endpoint_render_engine_set(endpoint, old_engine);
	bu_vls_printf(gd->gedp->ged_result_str,
		"host '%s' is unavailable or incompatible with renderer '%s'\n",
		bu_vls_cstr(&host_name),
		_dm_render_engine_name(requested_engine));
	bu_vls_free(&view_name);
	bu_vls_free(&host_name);
	bu_vls_free(&renderer_name);
	return BRLCAD_ERROR;
    }
    if (!brlobol_display_endpoint_render_engine_set(endpoint,
	    requested_engine)) {
	brlobol_display_endpoint_host_detach(endpoint);
	(void)brlobol_display_endpoint_render_engine_set(endpoint, old_engine);
	bu_vls_printf(gd->gedp->ged_result_str,
		"host '%s' does not support renderer '%s'\n",
		bu_vls_cstr(&host_name),
		_dm_render_engine_name(requested_engine));
	bu_vls_free(&view_name);
	bu_vls_free(&host_name);
	bu_vls_free(&renderer_name);
	return BRLCAD_ERROR;
    }
    const char *cbav[5] = {"dm", "open", "--host",
	bu_vls_cstr(&host_name), _dm_render_engine_name(requested_engine)};
    _dm_cmd_during_clbk(gd, 5, cbav);
    bu_vls_printf(gd->gedp->ged_result_str, "%s",
	    bu_vls_cstr(&host_name));
    bu_vls_free(&view_name);
    bu_vls_free(&host_name);
    bu_vls_free(&renderer_name);
    return BRLCAD_OK;
}

int
_dm_cmd_close(void *ds, int argc, const char **argv)
{
    const char *usage_string = "dm [options] close [-V view]";
    const char *purpose_string =
	"close the selected endpoint host while retaining its scene.";
    if (_dm_cmd_msgs(ds, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;
    struct _ged_dm_info *gd = (struct _ged_dm_info *)ds;
    struct bu_vls view_name = BU_VLS_INIT_ZERO;
    struct bu_opt_desc d[2];
    BU_OPT(d[0], "V", "view", "name", &bu_opt_vls, &view_name,
	    "view endpoint to close");
    BU_OPT_NULL(d[1]);
    int ac = bu_opt_parse(NULL, argc, argv, d);
    if (ac) {
	bu_vls_printf(gd->gedp->ged_result_str, "Usage: %s", usage_string);
	bu_vls_free(&view_name);
	return BRLCAD_ERROR;
    }
    brlobol_display_endpoint_t *endpoint =
	_dm_endpoint(gd, &view_name, NULL);
    bu_vls_free(&view_name);
    if (!endpoint)
	return BRLCAD_ERROR;
    brlobol_display_endpoint_host_detach(endpoint);
    const char *cbav[2] = {"dm", "close"};
    _dm_cmd_during_clbk(gd, 2, cbav);
    return BRLCAD_OK;
}

int
_dm_cmd_host(void *ds, int argc, const char **argv)
{
    const char *usage_string = "dm [options] host [-V view]";
    const char *purpose_string = "report the selected endpoint host.";
    if (_dm_cmd_msgs(ds, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;
    struct _ged_dm_info *gd = (struct _ged_dm_info *)ds;
    struct bu_vls view_name = BU_VLS_INIT_ZERO;
    struct bu_opt_desc d[2];
    BU_OPT(d[0], "V", "view", "name", &bu_opt_vls, &view_name,
	    "view endpoint to report");
    BU_OPT_NULL(d[1]);
    int ac = bu_opt_parse(NULL, argc, argv, d);
    if (ac) {
	bu_vls_printf(gd->gedp->ged_result_str, "Usage: %s", usage_string);
	bu_vls_free(&view_name);
	return BRLCAD_ERROR;
    }
    brlobol_display_endpoint_t *endpoint =
	_dm_endpoint(gd, &view_name, NULL);
    bu_vls_free(&view_name);
    if (!endpoint)
	return BRLCAD_ERROR;
    const char *host = brlobol_display_endpoint_host_factory_name(endpoint);
    bu_vls_printf(gd->gedp->ged_result_str, "%s", host ? host : "unbound");
    return BRLCAD_OK;
}

int
_dm_cmd_diagnostics(void *ds, int argc, const char **argv)
{
    const char *usage_string = "dm [options] diagnostics [-V view]";
    const char *purpose_string =
	"report endpoint identity, capabilities, and typed properties.";
    if (_dm_cmd_msgs(ds, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;
    struct _ged_dm_info *gd = (struct _ged_dm_info *)ds;
    struct bu_vls view_name = BU_VLS_INIT_ZERO;
    struct bu_opt_desc d[2];
    BU_OPT(d[0], "V", "view", "name", &bu_opt_vls, &view_name,
	    "view endpoint to diagnose");
    BU_OPT_NULL(d[1]);
    int ac = bu_opt_parse(NULL, argc, argv, d);
    if (ac) {
	bu_vls_printf(gd->gedp->ged_result_str, "Usage: %s", usage_string);
	bu_vls_free(&view_name);
	return BRLCAD_ERROR;
    }
    void *view = NULL;
    brlobol_display_endpoint_t *endpoint =
	_dm_endpoint(gd, &view_name, &view);
    bu_vls_free(&view_name);
    if (!endpoint)
	return BRLCAD_ERROR;

    if (brlobol_display_endpoint_render_engine_get(endpoint) ==
	BRLOBOL_RENDER_ENGINE_DIAGNOSTIC)
	(void)brlobol_display_endpoint_diagnostic_refresh(endpoint);

    const char *host = brlobol_display_endpoint_host_factory_name(endpoint);
    bu_vls_printf(gd->gedp->ged_result_str,
	    "view=%s\nhost=%s\nhost.capabilities=0x%016" PRIx64
	    "\nrenderer=%s\nrenderer.capabilities=0x%016" PRIx64
	    "\nrenderer.resolved=%s\ncontroller=%s\n",
	    _dm_view_name(view), host ? host : "unbound",
	    brlobol_display_endpoint_host_capabilities(endpoint),
	    _dm_render_engine_name(
		brlobol_display_endpoint_render_engine_get(endpoint)),
	    brlobol_display_endpoint_render_engine_capabilities(endpoint),
	    _dm_render_engine_name(
		brlobol_display_endpoint_render_engine_resolved_get(endpoint)),
	    brlobol_display_endpoint_controller(endpoint) ? "available" :
	    "missing");
    const size_t count = brlobol_display_endpoint_property_count();
    for (size_t i = 0; i < count; i++) {
	struct brlobol_endpoint_property_desc desc = {0};
	desc.struct_size = sizeof(desc);
	if (brlobol_display_endpoint_property_descriptor(i, &desc) !=
		BRLOBOL_ENDPOINT_PROPERTY_OK)
	    continue;
	bu_vls_printf(gd->gedp->ged_result_str,
		"property.%s type=%s access=%s%s required=0x%016" PRIx64,
		desc.name, _dm_property_type_name(desc.type),
		(desc.access & BRLOBOL_ENDPOINT_PROPERTY_READ) ? "r" : "",
		(desc.access & BRLOBOL_ENDPOINT_PROPERTY_WRITE) ? "w" : "",
		desc.required_host_capabilities);
	if (desc.access & BRLOBOL_ENDPOINT_PROPERTY_READ) {
	    struct brlobol_endpoint_property_value value =
		BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	    if (brlobol_display_endpoint_property_get(endpoint, desc.name,
		    &value) == BRLOBOL_ENDPOINT_PROPERTY_OK) {
		bu_vls_printf(gd->gedp->ged_result_str, " value=");
		_dm_property_print(gd->gedp->ged_result_str, &value);
	    }
	}
	if (desc.allowed_values)
	    bu_vls_printf(gd->gedp->ged_result_str, " allowed=%s",
		    desc.allowed_values);
	bu_vls_printf(gd->gedp->ged_result_str, "\n");
    }
    return BRLCAD_OK;
}

int
_dm_cmd_get(void *ds, int argc, const char **argv)
{
    const char *usage_string = "dm [options] get [-V view] [property ...]";
    const char *purpose_string = "report typed Obol endpoint properties.";
    if (_dm_cmd_msgs(ds, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--; argv++;

    struct bu_vls view_name = BU_VLS_INIT_ZERO;

    struct bu_opt_desc d[2];
    BU_OPT(d[0], "V", "view", "name", &bu_opt_vls, &view_name, "view endpoint to report");
    BU_OPT_NULL(d[1]);

    int ac = bu_opt_parse(NULL, argc, argv, d);

    struct _ged_dm_info *gd = (struct _ged_dm_info *)ds;
    void *view = _dm_endpoint_view(gd, &view_name);
    brlobol_display_endpoint_t *endpoint = view ?
	ged_view_context_display_endpoint_get(view) : NULL;
    if (endpoint) {
	if (!ac) {
	    const size_t count = brlobol_display_endpoint_property_count();
	    for (size_t i = 0; i < count; i++) {
		struct brlobol_endpoint_property_desc desc = {0};
		desc.struct_size = sizeof(desc);
		if (brlobol_display_endpoint_property_descriptor(i, &desc) !=
			BRLOBOL_ENDPOINT_PROPERTY_OK ||
			!(desc.access & BRLOBOL_ENDPOINT_PROPERTY_READ))
		    continue;
		struct brlobol_endpoint_property_value value =
		    BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
		if (brlobol_display_endpoint_property_get(endpoint, desc.name,
			    &value) != BRLOBOL_ENDPOINT_PROPERTY_OK)
		    continue;
		bu_vls_printf(gd->gedp->ged_result_str, "%s=", desc.name);
		_dm_property_print(gd->gedp->ged_result_str, &value);
		bu_vls_printf(gd->gedp->ged_result_str, "\n");
	    }
	    bu_vls_free(&view_name);
	    return BRLCAD_OK;
	}
	for (int i = 0; i < ac; i++) {
	    const char *name = argv[i];
	    struct brlobol_endpoint_property_value value =
		BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	    int ret = brlobol_display_endpoint_property_get(endpoint, name,
		    &value);
	    if (ret != BRLOBOL_ENDPOINT_PROPERTY_OK) {
		bu_vls_printf(gd->gedp->ged_result_str, "%s: %s\n", name,
			_dm_property_error(ret));
		bu_vls_free(&view_name);
		return BRLCAD_ERROR;
	    }
	    if (gd->verbosity || ac > 1)
		bu_vls_printf(gd->gedp->ged_result_str, "%s=", name);
	    if (value.type == BRLOBOL_ENDPOINT_PROPERTY_BOOL &&
		    !strchr(argv[i], '.'))
		bu_vls_printf(gd->gedp->ged_result_str, "%d", value.bool_value);
	    else
		_dm_property_print(gd->gedp->ged_result_str, &value);
	    if (i + 1 < ac)
		bu_vls_printf(gd->gedp->ged_result_str, "\n");
	}
	bu_vls_free(&view_name);
	return BRLCAD_OK;
    }
    (void)_dm_endpoint(gd, &view_name, NULL);
    bu_vls_free(&view_name);
    return BRLCAD_ERROR;
}

int
_dm_cmd_set(void *ds, int argc, const char **argv)
{
    const char *usage_string = "dm [options] set [-V view] property [value]";
    const char *purpose_string = "inspect or assign a typed Obol endpoint property.";
    if (_dm_cmd_msgs(ds, argc, argv, usage_string, purpose_string)) {
	return BRLCAD_OK;
    }

    argc--; argv++;

    struct bu_vls view_name = BU_VLS_INIT_ZERO;

    struct bu_opt_desc d[2];
    BU_OPT(d[0], "V", "view", "name", &bu_opt_vls, &view_name, "view endpoint to modify");
    BU_OPT_NULL(d[1]);

    int ac = bu_opt_parse(NULL, argc, argv, d);

    struct _ged_dm_info *gd = (struct _ged_dm_info *)ds;
    {
	const char *name = ac > 0 ? argv[0] : NULL;
	brlobol_display_endpoint_t *endpoint =
	    _dm_endpoint(gd, &view_name, NULL);
	if (!endpoint) {
	    bu_vls_free(&view_name);
	    return BRLCAD_ERROR;
	}
	if (ac < 1 || ac > 2) {
	    bu_vls_printf(gd->gedp->ged_result_str, "Usage: %s", usage_string);
	    bu_vls_free(&view_name);
	    return BRLCAD_ERROR;
	}
	struct brlobol_endpoint_property_value value =
	    BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	int ret = brlobol_display_endpoint_property_get(endpoint, name, &value);
	if (ret != BRLOBOL_ENDPOINT_PROPERTY_OK) {
	    bu_vls_printf(gd->gedp->ged_result_str, "%s: %s\n", name,
		    _dm_property_error(ret));
	    bu_vls_free(&view_name);
	    return BRLCAD_ERROR;
	}
	if (ac == 1) {
	    if (value.type == BRLOBOL_ENDPOINT_PROPERTY_BOOL &&
		    !strchr(argv[0], '.'))
		bu_vls_printf(gd->gedp->ged_result_str, "%d", value.bool_value);
	    else
		_dm_property_print(gd->gedp->ged_result_str, &value);
	    bu_vls_free(&view_name);
	    return BRLCAD_OK;
	}
	if (!_dm_property_value_parse(&value, argv[1])) {
	    bu_vls_printf(gd->gedp->ged_result_str,
		    "%s: invalid value '%s'\n", name, argv[1]);
	    bu_vls_free(&view_name);
	    return BRLCAD_ERROR;
	}
	ret = brlobol_display_endpoint_property_set(endpoint, name, &value);
	if (ret != BRLOBOL_ENDPOINT_PROPERTY_OK) {
	    bu_vls_printf(gd->gedp->ged_result_str, "%s: %s\n", name,
		    _dm_property_error(ret));
	    bu_vls_free(&view_name);
	    return BRLCAD_ERROR;
	}
	const char *cbav[4] = {"dm", "set", argv[0], argv[1]};
	_dm_cmd_during_clbk(gd, 4, cbav);
	bu_vls_free(&view_name);
	return BRLCAD_OK;
    }
}

const struct bu_cmdtab _dm_cmds[] = {
    { "bg",              _dm_cmd_bg},
    { "close",           _dm_cmd_close},
    { "diagnostics",     _dm_cmd_diagnostics},
    { "get",             _dm_cmd_get},
    { "host",            _dm_cmd_host},
    { "list",            _dm_cmd_list},
    { "open",            _dm_cmd_open},
    { "renderer",        _dm_cmd_renderer},
    { "set",             _dm_cmd_set},
    { "status",          _dm_cmd_status},
    { (char *)NULL,      NULL}
};

int
ged_dm_core(struct ged *gedp, int argc, const char *argv[])
{
    int help = 0;
    struct _ged_dm_info gd;
    gd.gedp = gedp;
    gd.cmds = _dm_cmds;
    gd.verbosity = 0;

    // Sanity
    if (UNLIKELY(!gedp || !argc || !argv)) {
	return BRLCAD_ERROR;
    }

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    // We know we're the dm command - start processing args
    argc--; argv++;

    // See if we have any high level options set
    struct bu_opt_desc d[3];
    BU_OPT(d[0], "h", "help",    "",  NULL,               &help,         "Print help");
    BU_OPT(d[1], "v", "verbose", "",  &bu_opt_incr_long,  &gd.verbosity, "Verbose output");
    BU_OPT_NULL(d[2]);

    gd.gopts = d;

    int ac = bu_opt_parse(NULL, argc, argv, d);

    if (!ac || help) {
	_ged_subcmd_help(gedp, (struct bu_opt_desc *)d, (const struct bu_cmdtab *)_dm_cmds, "dm", "[options] subcommand [args]", &gd, 0, NULL);
	return BRLCAD_OK;
    }

    int ret;
    if (bu_cmd(_dm_cmds, ac, argv, 0, (void *)&gd, &ret) == BRLCAD_OK) {
	return ret;
    } else {
	bu_vls_printf(gedp->ged_result_str, "subcommand %s not defined", argv[0]);
    }

    return BRLCAD_ERROR;
}


#include "../include/plugin.h"

extern int ged_ert_core(struct ged *gedp, int argc, const char *argv[]);
extern int ged_screen_grab_core(struct ged *gedp, int argc, const char *argv[]);

#define GED_DM_COMMANDS(X, XID) \
    X(ert, ged_ert_core, GED_CMD_DEFAULT) \
    X(dm, ged_dm_core, GED_CMD_DEFAULT) \
    X(screen_grab, ged_screen_grab_core, GED_CMD_DEFAULT) \
    X(screengrab, ged_screen_grab_core, GED_CMD_DEFAULT) \

GED_DECLARE_COMMAND_SET(GED_DM_COMMANDS)
GED_DECLARE_PLUGIN_MANIFEST("libged_dm", 1, GED_DM_COMMANDS)

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
