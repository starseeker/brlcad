/*                     G E D _ D R A W _ P E R F . C P P
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
/** @file ged_draw_perf.cpp
 *
 * Minimal libged/libBObol draw-stack timing probe.  This executable is
 * intentionally not a user shell: it opens a database, creates a GED-owned
 * Obol endpoint, runs one draw command, optionally autoviews/renders, reports
 * stage timings and a stable image hash, and exits.
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "bu/app.h"

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <bu.h>
#include <bu/hash.h>
#include <icv.h>

#include "bv.h"
#include <ged.h>
#include "BObol.h"
#include "BObol/BPerformance.h"
#include "BObol/BVListShape.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BPerformance.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "ged/view.h"

#include <Inventor/nodes/SoGroup.h>

struct options {
    bool autoview = false;
    bool render = false;
    bool profile = false;
    bool redraw = false;
    int width = 512;
    int height = 512;
    int repeat = 1;
    int progressive_frames = 1;
    int render_frames = 1;
    int style_updates = 0;
    bool clear_cache = false;
    std::string cache_dir;
    std::string screengrab;
    std::string db_path;
    std::vector<std::string> draw_args;
    int lod_mesh = -1;
    int lod_csg = -1;
    int lod_bot_threshold = -1;
    BObolViewController::SoftwareWireMode software_wire =
	BObolViewController::SOFTWARE_WIRE_AUTO;
};

struct timings {
    double open_ms = 0.0;
    double endpoint_ms = 0.0;
    double draw_ms = 0.0;
    double redraw_ms = 0.0;
    double autoview_ms = 0.0;
    double render_ms = 0.0;
    double render_first_ms = 0.0;
    double render_steady_ms = 0.0;
    double style_update_ms = 0.0;
    double style_render_ms = 0.0;
    double screengrab_ms = 0.0;
};

static double
elapsed_ms(std::chrono::steady_clock::time_point start)
{
    auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

static void
collect_database_sources(SoNode *node,
	std::vector<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;
    SoBRLDatabaseSource *source = dynamic_cast<SoBRLDatabaseSource *>(node);
    if (source) {
	sources.push_back(source);
	return;
    }
    SoGroup *group = dynamic_cast<SoGroup *>(node);
    if (!group)
	return;
    for (int i = 0; i < group->getNumChildren(); i++)
	collect_database_sources(group->getChild(i), sources);
}

static void
usage(FILE *fp)
{
    std::fprintf(fp,
	"Usage: ged_draw_perf [options] <model.g> <draw-arg> [draw-arg ...]\n"
	"\n"
	"Options:\n"
	"  --autoview           run autoview after draw\n"
	"  --render             render the GED-owned Obol endpoint to RGB\n"
	"  --profile            enable Obol and evaluated-wire stage counters\n"
	"  --redraw             run the front-end redraw-all transaction after draw\n"
	"  --screengrab <file>  write a screengrab after draw/autoview\n"
	"  --size <WxH>         render size (default: 512x512)\n"
	"  --repeat <N>         repeat full open/endpoint/draw/close cycle\n"
	"  --progressive-frames <N>\n"
	"                       render up to N progressive frames (default: 1)\n"
	"  --render-frames <N>  render at least N frames for steady-state timing\n"
	"  --style-updates <N>  mutate aggregate appearance and render N times\n"
	"  --cache-dir <dir>    set BU_DIR_CACHE for this run\n"
	"  --clear-cache        clear --cache-dir before the first iteration\n"
	"  --lod-mesh <0|1>     run 'view lod mesh <0|1>' before draw\n"
	"  --lod-csg <0|1>      run 'view lod csg <0|1>' before draw\n"
	"  --lod-bot-threshold <N>\n"
	"                       run 'view lod bot_threshold <N>' before draw\n"
	"  --software-wire <auto|quality|fast>\n"
	"                       choose the OSMesa wireframe path (default: auto)\n"
	"  --help               print this help\n"
	"\n"
	"Examples:\n"
	"  ged_draw_perf build/share/db/havoc.g havoc\n"
	"  ged_draw_perf --autoview --render build/share/db/havoc.g havoc\n"
	"  ged_draw_perf -- build/share/db/faa/Generic_Twin.g all\n");
}

static bool
parse_positive_int(const char *str, int *out)
{
    if (!str || !out)
	return false;
    char *endp = NULL;
    long val = std::strtol(str, &endp, 10);
    if (!endp || *endp != '\0' || val <= 0 || val > 1000000)
	return false;
    *out = (int)val;
    return true;
}

static bool
parse_nonnegative_int(const char *str, int *out)
{
    if (!str || !out)
	return false;
    char *endp = NULL;
    long val = std::strtol(str, &endp, 10);
    if (!endp || *endp != '\0' || val < 0 || val > 1000000)
	return false;
    *out = (int)val;
    return true;
}

static bool
parse_bool_int(const char *str, int *out)
{
    int val = -1;
    if (!parse_nonnegative_int(str, &val) || (val != 1 && val != 0))
	return false;
    *out = val;
    return true;
}

static bool
parse_size(const char *str, int *width, int *height)
{
    if (!str || !width || !height)
	return false;
    const char *x = std::strchr(str, 'x');
    if (!x)
	x = std::strchr(str, 'X');
    if (!x)
	return false;
    std::string ws(str, (size_t)(x - str));
    std::string hs(x + 1);
    int w = 0;
    int h = 0;
    if (!parse_positive_int(ws.c_str(), &w) ||
	    !parse_positive_int(hs.c_str(), &h))
	return false;
    *width = w;
    *height = h;
    return true;
}

static bool
parse_software_wire_mode(const char *str,
	BObolViewController::SoftwareWireMode *mode)
{
    if (!str || !mode)
	return false;
    if (BU_STR_EQUAL(str, "auto"))
	*mode = BObolViewController::SOFTWARE_WIRE_AUTO;
    else if (BU_STR_EQUAL(str, "quality"))
	*mode = BObolViewController::SOFTWARE_WIRE_QUALITY;
    else if (BU_STR_EQUAL(str, "fast"))
	*mode = BObolViewController::SOFTWARE_WIRE_FAST;
    else
	return false;
    return true;
}

static bool
parse_args(int argc, const char **argv, struct options *opts)
{
    if (!opts)
	return false;

    int i = 1;
    for (; i < argc; i++) {
	const char *arg = argv[i];
	if (!arg)
	    return false;
	if (BU_STR_EQUAL(arg, "--")) {
	    i++;
	    break;
	}
	if (arg[0] != '-')
	    break;
	if (BU_STR_EQUAL(arg, "--help") || BU_STR_EQUAL(arg, "-h")) {
	    usage(stdout);
	    std::exit(0);
	}
	if (BU_STR_EQUAL(arg, "--autoview")) {
	    opts->autoview = true;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--render")) {
	    opts->render = true;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--profile")) {
	    opts->profile = true;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--redraw")) {
	    opts->redraw = true;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--clear-cache")) {
	    opts->clear_cache = true;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--screengrab")) {
	    if (++i >= argc)
		return false;
	    opts->screengrab = argv[i];
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--cache-dir")) {
	    if (++i >= argc)
		return false;
	    opts->cache_dir = argv[i];
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--size")) {
	    if (++i >= argc ||
		    !parse_size(argv[i], &opts->width, &opts->height))
		return false;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--repeat")) {
	    if (++i >= argc || !parse_positive_int(argv[i], &opts->repeat))
		return false;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--progressive-frames")) {
	    if (++i >= argc ||
		    !parse_positive_int(argv[i], &opts->progressive_frames))
		return false;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--render-frames")) {
	    if (++i >= argc || !parse_positive_int(argv[i], &opts->render_frames))
		return false;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--style-updates")) {
	    if (++i >= argc || !parse_positive_int(argv[i], &opts->style_updates))
		return false;
	    opts->render = true;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--lod-mesh")) {
	    if (++i >= argc || !parse_bool_int(argv[i], &opts->lod_mesh))
		return false;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--lod-csg")) {
	    if (++i >= argc || !parse_bool_int(argv[i], &opts->lod_csg))
		return false;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--lod-bot-threshold")) {
	    if (++i >= argc ||
		    !parse_nonnegative_int(argv[i], &opts->lod_bot_threshold))
		return false;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--software-wire")) {
	    if (++i >= argc ||
		    !parse_software_wire_mode(argv[i], &opts->software_wire))
		return false;
	    continue;
	}
	return false;
    }

    if (i >= argc)
	return false;
    opts->db_path = argv[i++];
    for (; i < argc; i++)
	opts->draw_args.emplace_back(argv[i]);

    return !opts->db_path.empty() && !opts->draw_args.empty();
}

static bool
cache_dir_dangerous(const std::string &cache_dir)
{
    if (cache_dir.empty() || cache_dir == "/" || cache_dir == "." ||
	    cache_dir == "..")
	return true;
    return false;
}

static int
prepare_cache(const struct options &opts)
{
    if (opts.cache_dir.empty())
	return BRLCAD_OK;
    if (opts.clear_cache) {
	if (cache_dir_dangerous(opts.cache_dir)) {
	    std::fprintf(stderr,
		"ged_draw_perf: refusing to clear unsafe cache path '%s'\n",
		opts.cache_dir.c_str());
	    return BRLCAD_ERROR;
	}
	bu_dirclear(opts.cache_dir.c_str());
    }
    bu_mkdir(opts.cache_dir.c_str());
    bu_setenv("BU_DIR_CACHE", opts.cache_dir.c_str(), 1);
    return BRLCAD_OK;
}

static const char *
ged_result(struct ged *gedp)
{
    if (!gedp || !gedp->ged_result_str)
	return "";
    return bu_vls_cstr(gedp->ged_result_str);
}

static int
run_view_lod_setting(struct ged *gedp, const char *name, const char *value)
{
    const char *av[4] = {"view", "lod", name, value};
    return ged_exec_view(gedp, 4, av);
}

static int
apply_lod_options(struct ged *gedp, const struct options &opts)
{
    char bot_threshold[64] = {0};
    if (opts.lod_mesh >= 0) {
	const char *val = opts.lod_mesh ? "1" : "0";
	if (run_view_lod_setting(gedp, "mesh", val) != BRLCAD_OK)
	    return BRLCAD_ERROR;
    }
    if (opts.lod_csg >= 0) {
	const char *val = opts.lod_csg ? "1" : "0";
	if (run_view_lod_setting(gedp, "csg", val) != BRLCAD_OK)
	    return BRLCAD_ERROR;
    }
    if (opts.lod_bot_threshold >= 0) {
	std::snprintf(bot_threshold, sizeof(bot_threshold), "%d",
	    opts.lod_bot_threshold);
	if (run_view_lod_setting(gedp, "bot_threshold", bot_threshold) !=
		BRLCAD_OK)
	    return BRLCAD_ERROR;
    }
    return BRLCAD_OK;
}

static SoDB::ContextManager *
performance_context_manager(void)
{
    static SoDB::ContextManager *manager = SoDB::createOSMesaContextManager();
    return manager;
}

static BObolViewController *
initialize_endpoint(struct ged *gedp, const struct options &opts)
{
    if (!gedp || !gedp->dbip)
	return NULL;
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    struct bv *view = view_ctx ?
	bv_context_view((struct bv_context *)view_ctx) : NULL;
    if (!view_ctx || !view)
	return NULL;
    bv_dimensions_set(view, opts.width, opts.height);
    bv_unit_conversion_set(view,
	gedp->dbip->dbi_local2base,
	gedp->dbip->dbi_base2local);
    bobol_display_endpoint_t *endpoint =
	ged_view_context_obol_endpoint_get(view_ctx);
    if (!endpoint) {
	endpoint = bobol_display_endpoint_create(NULL, 0);
	if (!endpoint)
	    return NULL;
	if (!ged_view_context_obol_endpoint_set(view_ctx, endpoint, 1)) {
	    bobol_display_endpoint_destroy(endpoint);
	    return NULL;
	}
    }
    BObolViewController *controller = static_cast<BObolViewController *>(
	bobol_display_endpoint_controller(endpoint));
    if (!controller)
	return NULL;
    controller->setRenderContextManager(performance_context_manager());
    controller->setViewportSize((unsigned int)opts.width,
	(unsigned int)opts.height);
    controller->setSoftwareWireMode(opts.software_wire);
    if (!controller->syncCameraFromViewContext(view_ctx))
	return NULL;
    return controller;
}

static int
run_draw(struct ged *gedp, const struct options &opts)
{
    std::vector<const char *> av;
    av.reserve(opts.draw_args.size() + 1);
    av.push_back("draw");
    for (const std::string &arg : opts.draw_args)
	av.push_back(arg.c_str());
    return ged_exec_draw(gedp, (int)av.size(), av.data());
}

static int
run_autoview(struct ged *gedp)
{
    const char *av[2] = {"autoview", NULL};
    return ged_exec_autoview(gedp, 1, av);
}

static int
run_render(struct ged *gedp, BObolViewController *controller,
	   unsigned char **image, BObolProgressiveStatus *status)
{
    struct ged_view_context *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx || !controller || !image)
	return BRLCAD_ERROR;
    if (!controller->syncCameraFromViewContext(view_ctx))
	return BRLCAD_ERROR;
    (void)ged_view_faceplate_sync(gedp, view_ctx);
    return controller->renderToImage(image, 1, 0, NULL,
	performance_context_manager(), status);
}

static int
write_screengrab(const unsigned char *image, int width, int height,
		 const char *filename)
{
    if (!image || width <= 0 || height <= 0 || !filename || !filename[0])
	return BRLCAD_ERROR;
    icv_image_t *out = icv_create((size_t)width, (size_t)height,
	ICV_COLOR_SPACE_RGB);
    if (!out)
	return BRLCAD_ERROR;
    int ret = BRLCAD_OK;
    const size_t stride = (size_t)width * 3;
    for (int y = 0; y < height; y++) {
	if (icv_writeline(out, (size_t)y,
		(void *)(image + (size_t)y * stride), ICV_DATA_UCHAR) != 0) {
	    ret = BRLCAD_ERROR;
	    break;
	}
    }
    if (ret == BRLCAD_OK &&
	icv_write(out, filename, BU_MIME_IMAGE_PNG) != 0)
	ret = BRLCAD_ERROR;
    icv_destroy(out);
    return ret;
}

static int
run_once(const struct options &opts, int iter)
{
    struct timings t;
    struct BObolPerformanceCounters perf;
    bobol_performance_counters_init(&perf);
    int endpoint_ret = BRLCAD_OK;
    int lod_ret = BRLCAD_OK;
    int draw_ret = BRLCAD_OK;
    int redraw_ret = BRLCAD_OK;
    int autoview_ret = BRLCAD_OK;
    int render_ret = BRLCAD_OK;
    int style_ret = BRLCAD_OK;
    int screengrab_ret = BRLCAD_OK;
    int source_count_ret = 0;
    int profile_started = 0;
    size_t source_count = 0;
    unsigned char *image = NULL;
    uint64_t image_hash = 0;
    size_t image_nonzero = 0;
    int progressive_frames = 0;
    int progressive_more = 0;
    size_t progressive_expanded = 0;
    size_t progressive_remaining = 0;
    size_t profile_sources_realized = 0;
    size_t profile_sources_stale = 0;
    size_t profile_sources_view_dependent = 0;
    size_t profile_sources_auxiliary = 0;
    size_t profile_sources_compact = 0;
    size_t profile_sources_compiled = 0;
    size_t profile_primary_shapes = 0;
    size_t profile_primary_points = 0;
    size_t profile_auxiliary_shapes = 0;
    size_t profile_auxiliary_points = 0;
    size_t profile_point_shapes = 0;
    size_t profile_point_commands = 0;
    size_t profile_styled_shapes = 0;
    size_t profile_emphasis_shapes = 0;
    BObolViewController *controller = NULL;

    auto start = std::chrono::steady_clock::now();
    struct ged *gedp = ged_open("db", opts.db_path.c_str(), 1);
    t.open_ms = elapsed_ms(start);
    if (!gedp) {
	std::fprintf(stderr, "ged_draw_perf: failed to open %s\n",
	    opts.db_path.c_str());
	return 1;
    }

    start = std::chrono::steady_clock::now();
    controller = initialize_endpoint(gedp, opts);
    endpoint_ret = controller ? BRLCAD_OK : BRLCAD_ERROR;
    t.endpoint_ms = elapsed_ms(start);

    if (endpoint_ret == BRLCAD_OK) {
	lod_ret = apply_lod_options(gedp, opts);
    }

    if (endpoint_ret == BRLCAD_OK && lod_ret == BRLCAD_OK) {
	if (opts.profile) {
	    bobol_performance_counters_reset();
	    bobol_performance_counters_set_enabled(1);
	    profile_started = 1;
	}
	start = std::chrono::steady_clock::now();
	draw_ret = run_draw(gedp, opts);
	t.draw_ms = elapsed_ms(start);
    }

    if (draw_ret == BRLCAD_OK) {
	if (opts.redraw) {
	    struct ged_draw_transaction txn =
		ged_draw_transaction_make(GED_DRAW_TXN_REDRAW, NULL);
	    struct ged_draw_transaction_result result;
	    ged_draw_transaction_result_init(&result);
	    start = std::chrono::steady_clock::now();
	    const int txn_ret = ged_draw_apply_transaction(gedp, &txn, &result);
	    t.redraw_ms = elapsed_ms(start);
	    ged_draw_transaction_result_free(&result);
	    redraw_ret = txn_ret < 0 ? BRLCAD_ERROR : BRLCAD_OK;
	}

	const int semantic_shape_count = ged_draw_shape_count(gedp);
	if (semantic_shape_count >= 0) {
	    source_count = (size_t)semantic_shape_count;
	    source_count_ret = BRLCAD_OK;
	} else {
	    source_count_ret = BRLCAD_ERROR;
	}
    }

    if (draw_ret == BRLCAD_OK && opts.autoview) {
	start = std::chrono::steady_clock::now();
	autoview_ret = run_autoview(gedp);
	t.autoview_ms = elapsed_ms(start);
    }

    if (draw_ret == BRLCAD_OK && autoview_ret == BRLCAD_OK &&
	(opts.render || !opts.screengrab.empty())) {
	start = std::chrono::steady_clock::now();
	const int frameLimit = std::max(opts.progressive_frames,
	    opts.render_frames);
	for (int frame = 0; frame < frameLimit; frame++) {
	    if (image) {
		bu_free(image, "ged_draw_perf intermediate image");
		image = NULL;
	    }
	    BObolProgressiveStatus status;
	    const auto frameStart = std::chrono::steady_clock::now();
	    render_ret = run_render(gedp, controller, &image, &status);
	    const double frameMs = elapsed_ms(frameStart);
	    if (render_ret != BRLCAD_OK)
		break;
	    if (!progressive_frames)
		t.render_first_ms = frameMs;
	    else
		t.render_steady_ms += frameMs;
	    progressive_frames++;
	    progressive_more = status.hasMore ? 1 : 0;
	    progressive_expanded += status.expanded;
	    progressive_remaining = status.remaining;
	    if (progressive_frames >= opts.render_frames && !status.hasMore)
		break;
	}
	t.render_ms = elapsed_ms(start);
	if (progressive_frames > 1)
	    t.render_steady_ms /= progressive_frames - 1;
    }

    if (opts.style_updates > 0 && draw_ret == BRLCAD_OK &&
	autoview_ret == BRLCAD_OK && render_ret == BRLCAD_OK) {
	std::vector<SoBRLDatabaseSource *> sources;
	if (controller)
	    collect_database_sources(controller->getRenderSceneRoot(), sources);
	if (sources.empty() || progressive_more) {
	    style_ret = BRLCAD_ERROR;
	} else {
	    if (profile_started)
		bobol_performance_counters_reset();
	    for (int update = 0; update < opts.style_updates; update++) {
		const bool alternate = (update & 1) != 0;
		const auto updateStart = std::chrono::steady_clock::now();
		for (SoBRLDatabaseSource *source : sources) {
		    if (!source)
			continue;
		    BObolDatabaseSourceDisplayPatch patch;
		    patch.colorOverrideValid = TRUE;
		    patch.colorOverride = TRUE;
		    patch.colorValid = TRUE;
		    patch.color = alternate ? SbColor(0.12f, 0.72f, 0.24f) :
			SbColor(0.88f, 0.16f, 0.10f);
		    patch.lineWidthValid = TRUE;
		    patch.lineWidth = alternate ? 2 : 1;
		    patch.lineStyleValid = TRUE;
		    patch.lineStyle = alternate ? 1 : 0;
		    (void)source->applyDisplayPatch(patch);
		    (void)source->setCompactInstanceDisplayStateForPath("", TRUE,
			0, FALSE, 1, alternate ? TRUE : FALSE,
			1, alternate ? FALSE : TRUE);
		}
		t.style_update_ms += elapsed_ms(updateStart);

		if (image) {
		    bu_free(image, "ged_draw_perf style image");
		    image = NULL;
		}
		BObolProgressiveStatus styleStatus;
		const auto renderStart = std::chrono::steady_clock::now();
		style_ret = run_render(gedp, controller, &image, &styleStatus);
		t.style_render_ms += elapsed_ms(renderStart);
		if (style_ret != BRLCAD_OK || styleStatus.hasMore) {
		    style_ret = BRLCAD_ERROR;
		    break;
		}
	    }
	    t.style_update_ms /= opts.style_updates;
	    t.style_render_ms /= opts.style_updates;
	}
    }

    if (draw_ret == BRLCAD_OK && autoview_ret == BRLCAD_OK &&
	    render_ret == BRLCAD_OK && style_ret == BRLCAD_OK &&
	    !opts.screengrab.empty()) {
	start = std::chrono::steady_clock::now();
	screengrab_ret = write_screengrab(image, opts.width, opts.height,
	    opts.screengrab.c_str());
	t.screengrab_ms = elapsed_ms(start);
    }

    if (image) {
	const size_t image_bytes = (size_t)opts.width * (size_t)opts.height * 3;
	image_hash = bu_data_hash(image, image_bytes);
	for (size_t i = 0; i < image_bytes; i++)
	    image_nonzero += image[i] != 0;
    }

    if (profile_started) {
	bobol_performance_counters_set_enabled(0);
	bobol_performance_counters_get(&perf);
        if (controller) {
	    std::vector<SoBRLDatabaseSource *> sources;
	    collect_database_sources(controller->getRenderSceneRoot(), sources);
	    for (SoBRLDatabaseSource *source : sources) {
		if (!source)
		    continue;
		profile_sources_realized += source->realizationStatus.getValue() ==
		    SoBRLDatabaseSource::REALIZED;
		profile_sources_stale += source->needsRealization() ? 1 : 0;
		profile_sources_view_dependent +=
		    source->realizationViewDependent.getValue() ? 1 : 0;
		profile_sources_auxiliary += source->auxiliarySource.getValue() ? 1 : 0;
		profile_sources_compact += source->hasCompactInstanceIndex() ? 1 : 0;
		profile_sources_compiled += source->hasCompiledAssembly() ? 1 : 0;
		for (int j = 0; j < source->getRealizedShapeSummaryCount(); j++) {
		    SoBRLVListShape *shape = source->getRealizedShape(j);
		    BObolRealizedShapeSummary summary;
		    if (!source->getRealizedShapeSummary(j, summary))
			continue;
		    const bool auxiliary = BU_STR_EQUAL(
			summary.recordRole.getString(), "auxiliary");
		    if (auxiliary) {
			profile_auxiliary_shapes++;
			profile_auxiliary_points += summary.pointCount;
		    } else {
			profile_primary_shapes++;
			profile_primary_points += summary.pointCount;
		    }
		    if (shape) {
			bool hasPoints = false;
			for (int k = 0; k < shape->command.getNum(); k++) {
			    if (shape->command[k] == SoBRLVListShape::POINT) {
				hasPoints = true;
				profile_point_commands++;
			    }
			}
			profile_point_shapes += hasPoints ? 1 : 0;
			profile_styled_shapes += shape->lineStyle.getValue() != 0 ? 1 : 0;
			profile_emphasis_shapes +=
			    (shape->hiddenLine.getValue() ||
			     shape->editEmphasis.getValue() ||
			     shape->selectedPrimitive.getNum() > 0 ||
			     shape->highlightedPrimitive.getNum() > 0) ? 1 : 0;
		    }
		}
	    }
	}
    }

    const char *result = ged_result(gedp);
    if (endpoint_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: endpoint setup failed: %s\n", result);
    if (lod_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: LoD setup failed: %s\n", result);
    if (draw_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: draw failed: %s\n", result);
    if (redraw_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: redraw failed: %s\n", result);
    if (autoview_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: autoview failed: %s\n", result);
    if (render_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: render failed: %s\n", result);
    if (style_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: style churn failed: %s\n", result);
    if (screengrab_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: screengrab failed: %s\n", result);

    std::printf("iter=%d open_ms=%.3f endpoint_ms=%.3f draw_ms=%.3f",
	iter, t.open_ms, t.endpoint_ms, t.draw_ms);
    if (opts.redraw)
	std::printf(" redraw_ms=%.3f", t.redraw_ms);
    if (opts.autoview)
	std::printf(" autoview_ms=%.3f", t.autoview_ms);
    if (opts.render)
	std::printf(" render_ms=%.3f render_frame_ms=%.3f render_first_ms=%.3f render_steady_ms=%.3f render_frames=%d progressive_more=%d progressive_expanded=%zu progressive_remaining=%zu",
	    t.render_ms, progressive_frames ? t.render_ms / progressive_frames : 0.0,
	    t.render_first_ms, t.render_steady_ms, progressive_frames, progressive_more,
	    progressive_expanded, progressive_remaining);
    if (opts.style_updates > 0)
	std::printf(" style_updates=%d style_update_ms=%.3f style_render_ms=%.3f",
	    opts.style_updates, t.style_update_ms, t.style_render_ms);
    if (!opts.screengrab.empty())
	std::printf(" screengrab_ms=%.3f", t.screengrab_ms);
    if (image)
	std::printf(" image_hash=%016" PRIx64 " image_nonzero=%zu",
	    image_hash, image_nonzero);
    if (source_count_ret)
	std::printf(" obol_sources=%zu", source_count);
    if (profile_started) {
	std::printf(" perf_realize_ms=%.3f perf_seed_ms=%.3f perf_walk_ms=%.3f",
	    perf.realize_total_us / 1000.0, perf.realize_seed_us / 1000.0,
	    perf.realize_walk_us / 1000.0);
	std::printf(" perf_sources=%llu/%llu/%llu",
	    (unsigned long long)perf.sources_visited,
	    (unsigned long long)perf.sources_realized,
	    (unsigned long long)perf.sources_failed);
	std::printf(" perf_wire_ms=%.3f perf_wire_calls=%llu",
	    perf.wire_realize_us / 1000.0,
	    (unsigned long long)perf.wire_realize_calls);
	std::printf(" perf_mesh_ms=%.3f perf_mesh_calls=%llu",
	    perf.mesh_realize_us / 1000.0,
	    (unsigned long long)perf.mesh_realize_calls);
	std::printf(" perf_direct_leaf=%llu/%llu/%llu",
	    (unsigned long long)perf.direct_leaf_realized,
	    (unsigned long long)perf.direct_leaf_failed,
	    (unsigned long long)perf.direct_leaf_fallback);
	std::printf(" perf_wire_cache=%llu/%llu",
	    (unsigned long long)perf.wire_cache_hits,
	    (unsigned long long)perf.wire_cache_misses);
	std::printf(" perf_mesh_cache=%llu/%llu",
	    (unsigned long long)perf.mesh_cache_hits,
	    (unsigned long long)perf.mesh_cache_misses);
	std::printf(" perf_plot_ms=%.3f perf_plot_calls=%llu",
	    perf.plot_us / 1000.0,
	    (unsigned long long)perf.plot_calls);
	std::printf(" perf_vlist_ms=%.3f perf_vlist_calls=%llu perf_vlist_points=%llu",
	    perf.vlist_convert_us / 1000.0,
	    (unsigned long long)perf.vlist_convert_calls,
	    (unsigned long long)perf.vlist_points);
	std::printf(" perf_nodes=%llu perf_node_ms=%.3f",
	    (unsigned long long)perf.realized_instance_nodes,
	    perf.realized_instance_node_us / 1000.0);
	std::printf(" perf_cad_compact=%llu/%llu/%llu perf_cad_compact_ms=%.3f",
	    (unsigned long long)perf.cad_compact_attempts,
	    (unsigned long long)perf.cad_compact_sources,
	    (unsigned long long)perf.cad_compact_instances,
	    perf.cad_compact_us / 1000.0);
	std::printf(" perf_replace=%llu/%.3f perf_state=%llu/%.3f perf_move=%llu/%.3f perf_index=%llu/%.3f",
	    (unsigned long long)perf.source_replace_calls,
	    perf.source_replace_us / 1000.0,
	    (unsigned long long)perf.source_state_calls,
	    perf.source_state_us / 1000.0,
	    (unsigned long long)perf.source_move_calls,
	    perf.source_move_us / 1000.0,
	    (unsigned long long)perf.source_index_rebuild_calls,
	    perf.source_index_rebuild_us / 1000.0);
	std::printf(" perf_source_state=%zu/%zu/%zu/%zu perf_source_residency=%zu/%zu perf_shape_points=%zu/%zu/%zu/%zu",
	    profile_sources_realized, profile_sources_stale,
	    profile_sources_view_dependent, profile_sources_auxiliary,
	    profile_sources_compact, profile_sources_compiled,
	    profile_primary_shapes, profile_primary_points,
	    profile_auxiliary_shapes, profile_auxiliary_points);
	std::printf(" perf_shape_features=%zu/%zu/%zu/%zu",
	    profile_point_shapes, profile_point_commands,
	    profile_styled_shapes, profile_emphasis_shapes);
    }
    std::printf(" status=%s\n",
	(endpoint_ret == BRLCAD_OK && draw_ret == BRLCAD_OK &&
	 redraw_ret == BRLCAD_OK &&
	 lod_ret == BRLCAD_OK &&
	 autoview_ret == BRLCAD_OK && render_ret == BRLCAD_OK &&
	 style_ret == BRLCAD_OK &&
	 screengrab_ret == BRLCAD_OK) ? "ok" : "error");
    if (controller && controller->getLastDiagnostics().getLength() > 0)
	std::fprintf(stderr, "ged_draw_perf: Obol realization diagnostics: %s\n",
	    controller->getLastDiagnostics().getString());

    int ret = (endpoint_ret == BRLCAD_OK && lod_ret == BRLCAD_OK &&
	draw_ret == BRLCAD_OK &&
	redraw_ret == BRLCAD_OK &&
	autoview_ret == BRLCAD_OK && render_ret == BRLCAD_OK &&
	style_ret == BRLCAD_OK &&
	screengrab_ret == BRLCAD_OK) ? 0 : 1;

    if (image)
	bu_free(image, "ged_draw_perf image");
    ged_close(gedp);
    return ret;
}

int
main(int argc, const char **argv)
{
    bu_setprogname(argv[0]);
    struct options opts;
    if (!parse_args(argc, argv, &opts)) {
	usage(stderr);
	return 1;
    }
    if (prepare_cache(opts) != BRLCAD_OK)
	return 1;
    if (opts.profile)
	bu_setenv("RT_EVAL_WIREFRAME_PROFILE", "1", 1);

    SoDB::ContextManager *manager = performance_context_manager();
    if (!manager) {
	std::fprintf(stderr, "ged_draw_perf: OSMesa context manager unavailable\n");
	return 1;
    }
    bobol_init(NULL);

    std::printf("ged_draw_perf db=%s endpoint=obol draw_args=",
	opts.db_path.c_str());
    for (size_t i = 0; i < opts.draw_args.size(); i++)
	std::printf("%s%s", i ? " " : "", opts.draw_args[i].c_str());
    std::printf(" repeat=%d size=%dx%d autoview=%d render=%d profile=%d redraw=%d style_updates=%d screengrab=%s",
	opts.repeat, opts.width, opts.height, opts.autoview ? 1 : 0,
	opts.render ? 1 : 0, opts.profile ? 1 : 0,
	opts.redraw ? 1 : 0, opts.style_updates,
	opts.screengrab.empty() ? "none" : opts.screengrab.c_str());
    std::printf(" cache_dir=%s clear_cache=%d lod_mesh=%d lod_csg=%d lod_bot_threshold=%d\n",
	opts.cache_dir.empty() ? "ambient" : opts.cache_dir.c_str(),
	opts.clear_cache ? 1 : 0, opts.lod_mesh, opts.lod_csg,
	opts.lod_bot_threshold);

    int ret = 0;
    for (int i = 1; i <= opts.repeat; i++) {
	if (run_once(opts, i) != 0)
	    ret = 1;
    }

    return ret;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
