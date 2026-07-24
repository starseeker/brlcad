/*                    E X T E R N A L _ R T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libged/external_rt.cpp
 *
 * Launch an external image renderer into an endpoint-owned imgstream
 * framebuffer.  This is common process plumbing, not an interactive renderer
 * engine: the authoritative scene remains Obol-owned and the child remains a
 * normal rt-family executable.
 */

#include "common.h"

#include <climits>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "bu/app.h"
#include "bu/env.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "bv.h"
#include "imgstream/fbserv.h"
#include "pkg.h"

#include "./ged_private.h"


static bool
positive_dimension(const char *value, int *dimension)
{
    if (!value || !value[0] || !dimension)
	return false;

    char *end = NULL;
    long parsed = strtol(value, &end, 10);
    if (!end || end == value || end[0] || parsed <= 0 || parsed > INT_MAX)
	return false;

    *dimension = (int)parsed;
    return true;
}


static void
requested_dimensions(int argc, const char *argv[], int *width, int *height)
{
    for (int i = 1; i < argc; i++) {
	if (BU_STR_EQUAL(argv[i], "--"))
	    break;

	const char *value = NULL;
	char dimension = '\0';
	if (BU_STR_EQUAL(argv[i], "-s") || BU_STR_EQUAL(argv[i], "--size")) {
	    dimension = 's';
	    if (++i < argc)
		value = argv[i];
	} else if (BU_STR_EQUAL(argv[i], "-w") ||
		BU_STR_EQUAL(argv[i], "--width")) {
	    dimension = 'w';
	    if (++i < argc)
		value = argv[i];
	} else if (BU_STR_EQUAL(argv[i], "-n") ||
		BU_STR_EQUAL(argv[i], "--height")) {
	    dimension = 'n';
	    if (++i < argc)
		value = argv[i];
	} else if (argv[i][0] == '-' && argv[i][1] && argv[i][2] &&
		(argv[i][1] == 's' || argv[i][1] == 'w' || argv[i][1] == 'n')) {
	    dimension = argv[i][1];
	    value = argv[i] + 2;
	} else if (bu_strncmp(argv[i], "--size=", 7) == 0) {
	    dimension = 's';
	    value = argv[i] + 7;
	} else if (bu_strncmp(argv[i], "--width=", 8) == 0) {
	    dimension = 'w';
	    value = argv[i] + 8;
	} else if (bu_strncmp(argv[i], "--height=", 9) == 0) {
	    dimension = 'n';
	    value = argv[i] + 9;
	}

	int parsed = 0;
	if (!positive_dimension(value, &parsed))
	    continue;
	if (dimension == 's') {
	    *width = parsed;
	    *height = parsed;
	} else if (dimension == 'w') {
	    *width = parsed;
	} else if (dimension == 'n') {
	    *height = parsed;
	}
    }
}


extern "C" GED_EXPORT int
_ged_external_rt_to_endpoint(struct ged *gedp, int argc, const char *argv[],
	const char *program, const char *callback_command)
{
    if (!program || !program[0]) {
	if (gedp)
	    bu_vls_printf(gedp->ged_result_str,
		    "no external renderer program specified\n");
	return BRLCAD_ERROR;
    }

    GED_CHECK_DATABASE_OPEN(gedp, BRLCAD_ERROR);
    GED_CHECK_DRAWABLE(gedp, BRLCAD_ERROR);
    GED_CHECK_VIEW(gedp, BRLCAD_ERROR);
    GED_CHECK_ARGC_GT_0(gedp, argc, BRLCAD_ERROR);

    bu_vls_trunc(gedp->ged_result_str, 0);

    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx || !ged_view_context_display_endpoint_get(view_ctx)) {
	bu_vls_printf(gedp->ged_result_str,
		"active view has no Obol display endpoint\n");
	return BRLCAD_ERROR;
    }
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    if (!view) {
	bu_vls_printf(gedp->ged_result_str, "active view is unavailable\n");
	return BRLCAD_ERROR;
    }

    struct fbserv_obj *fbs = gedp->ged_fbs;
    if (!fbs) {
	bu_vls_printf(gedp->ged_result_str,
		"no framebuffer server configured\n");
	return BRLCAD_ERROR;
    }

    /* The GED owns one framebuffer bridge, while a GUI may own several
     * independently sized views.  Rebind unconditionally so this render uses
     * the active endpoint and its current toolkit-reported dimensions rather
     * than a construction-time size or the previously active pane. */
    struct fbserv_fb_info fbinfo;
    int have_backend =
	ged_draw_obol_framebuffer_backend_ensure_for_view(gedp, view_ctx) ==
	    BRLCAD_OK && fbs_framebuffer_info(fbs, &fbinfo) == 0;
    if (!have_backend) {
	bu_vls_printf(gedp->ged_result_str,
		"view has no endpoint-backed framebuffer stream\n");
	return BRLCAD_ERROR;
    }

    if (!ged_who_argc(gedp)) {
	bu_vls_printf(gedp->ged_result_str, "no objects displayed\n");
	return BRLCAD_ERROR;
    }

    struct bobol_endpoint_property_value framebuffer_mode =
	BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    if (ged_view_context_display_property_get(view_ctx,
		"composition.framebuffer.mode", &framebuffer_mode) !=
	BOBOL_ENDPOINT_PROPERTY_OK) {
	bu_vls_printf(gedp->ged_result_str,
		"active view has no Obol framebuffer composition policy\n");
	return BRLCAD_ERROR;
    }
    if (BU_STR_EQUAL(framebuffer_mode.string_value, "off")) {
	framebuffer_mode.type = BOBOL_ENDPOINT_PROPERTY_ENUM;
	framebuffer_mode.string_value = "underlay";
	if (ged_view_context_display_property_set(view_ctx,
		    "composition.framebuffer.mode", &framebuffer_mode) !=
	    BOBOL_ENDPOINT_PROPERTY_OK) {
	    bu_vls_printf(gedp->ged_result_str,
		    "unable to enable the Obol framebuffer underlay\n");
	    return BRLCAD_ERROR;
	}
    }

    bool using_ipc = false;
    if (fbs_can_open_ipc(fbs) && fbs_open_ipc(fbs) == BRLCAD_OK) {
	using_ipc = true;
    } else if (!fbs_can_open_network(fbs) || fbs_open(fbs, 0) != BRLCAD_OK) {
	bu_vls_printf(gedp->ged_result_str, "could not open fb server\n");
	return BRLCAD_ERROR;
    }

    int width = fbinfo.width;
    int height = fbinfo.height;
    if (width <= 0 || height <= 0) {
	bu_vls_printf(gedp->ged_result_str,
		"invalid embedded framebuffer dimensions\n");
	return BRLCAD_ERROR;
    }
    requested_dimensions(argc, argv, &width, &height);

    char executable[MAXPATHLEN] = {0};
    bu_dir(executable, MAXPATHLEN, BU_DIR_BIN, program, BU_DIR_EXT, NULL);

    std::vector<std::string> args;
    args.emplace_back(executable);
    args.emplace_back("-F");
    args.emplace_back(using_ipc ? "0" :
	    std::to_string(fbs_listener_port(fbs)));
    args.emplace_back("-M");
    args.emplace_back("-w");
    args.emplace_back(std::to_string(width));
    args.emplace_back("-n");
    args.emplace_back(std::to_string(height));
    args.emplace_back("-V");
    args.emplace_back(std::to_string((double)width / (double)height));
    const fastf_t perspective = bv_perspective_get(view);
    if (perspective > 0) {
	args.emplace_back("-p");
	args.emplace_back(std::to_string(perspective));
    }

    int units_supplied = 0;
    int object_arg = 1;
    for (; object_arg < argc; object_arg++) {
	if (argv[object_arg][0] == '-' && argv[object_arg][1] == '-' &&
	    argv[object_arg][2] == '\0') {
	    object_arg++;
	    break;
	}
	if (BU_STR_EQUAL(argv[object_arg], "-u") ||
	    BU_STR_EQUAL(argv[object_arg], "--units") ||
	    (argv[object_arg][0] == '-' && argv[object_arg][1] == 'u' &&
	     strlen(argv[object_arg]) > 2) ||
	    bu_strncmp(argv[object_arg], "--units=", 8) == 0)
	    units_supplied = 1;
	args.emplace_back(argv[object_arg]);
    }
    if (!units_supplied) {
	args.emplace_back("-u");
	args.emplace_back("model");
    }
    args.emplace_back(gedp->dbip->dbi_filename);

    std::vector<const char *> child_argv;
    child_argv.reserve(args.size() + 1);
    for (const std::string &arg : args)
	child_argv.push_back(arg.c_str());
    child_argv.push_back(NULL);

    bu_clbk_t linger_callback = NULL;
    void *linger_data = NULL;
    int child_pid = -1;
    if (callback_command)
	ged_clbk_get(&linger_callback, &linger_data, gedp, callback_command,
		BU_CLBK_LINGER);

    const char *saved_pkg_addr = getenv(PKG_ADDR_ENVVAR);
    const bool had_pkg_addr = saved_pkg_addr != NULL;
    const std::string saved_pkg_addr_value = had_pkg_addr ?
	std::string(saved_pkg_addr) : std::string();
    if (using_ipc) {
	const char *addr_env = fbs_ipc_child_addr_env(fbs);
	if (!addr_env) {
	    bu_vls_printf(gedp->ged_result_str,
		    "unable to get the endpoint IPC address\n");
	    return BRLCAD_ERROR;
	} else {
	    const char *eq = strchr(addr_env, '=');
	    if (!eq || !eq[1]) {
		bu_vls_printf(gedp->ged_result_str,
			"invalid endpoint IPC address\n");
		return BRLCAD_ERROR;
	    }
	    (void)bu_setenv(PKG_ADDR_ENVVAR, eq + 1, 1);
	}
    }

    bu_log("%s: launching endpoint framebuffer renderer (ipc=%d size=%dx%d)\n",
	program, using_ipc ? 1 : 0, width, height);
    const int ret = _ged_run_rt(gedp, (int)args.size(), child_argv.data(),
	argc - object_arg, &(argv[object_arg]), 0,
	callback_command ? &child_pid : NULL, linger_callback, linger_data);

    if (using_ipc)
	(void)bu_setenv(PKG_ADDR_ENVVAR, had_pkg_addr ?
		saved_pkg_addr_value.c_str() : "", 1);

    if (callback_command) {
	bu_clbk_t during_callback = NULL;
	void *during_data = NULL;
	ged_clbk_get(&during_callback, &during_data, gedp, callback_command,
		BU_CLBK_DURING);
	if (during_callback)
	    (*during_callback)(argc, argv, (void *)&child_pid, during_data);
    }

    return ret;
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
