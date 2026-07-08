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
 * Minimal libged/libbrlobol draw-stack timing probe.  This executable is
 * intentionally not a user shell: it opens a database, attaches an Obol DM by
 * default, runs one draw command, optionally autoviews/renders, reports stage
 * timings, and exits.
 */

#include "common.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <bu.h>

#define DM_WITH_RT
#include "dm.h"

#include "bv.h"
#include <ged.h>
#include "brlobol/performance.h"
#include "brlobol/scene_controller.h"
#include "ged/draw_obol.h"

struct options {
    bool attach_dm = true;
    bool autoview = false;
    bool render = false;
    bool profile = false;
    int width = 512;
    int height = 512;
    int repeat = 1;
    bool clear_cache = false;
    std::string dm_backend = "obol";
    std::string dm_name = "OBOL";
    std::string cache_dir;
    std::string screengrab;
    std::string db_path;
    std::vector<std::string> draw_args;
    int lod_mesh = -1;
    int lod_csg = -1;
    int lod_bot_threshold = -1;
    int cad_compact = -1;
};

struct timings {
    double open_ms = 0.0;
    double attach_ms = 0.0;
    double draw_ms = 0.0;
    double autoview_ms = 0.0;
    double render_ms = 0.0;
    double screengrab_ms = 0.0;
};

static double
elapsed_ms(std::chrono::steady_clock::time_point start)
{
    auto end = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

static void
usage(FILE *fp)
{
    std::fprintf(fp,
	"Usage: ged_draw_perf [options] <model.g> <draw-arg> [draw-arg ...]\n"
	"\n"
	"Options:\n"
	"  --no-dm              run draw without attaching a display manager\n"
	"  --dm <backend>       display backend to attach (default: obol)\n"
	"  --dm-name <name>     attached display name (default: OBOL)\n"
	"  --autoview           run autoview after draw\n"
	"  --render             call dm_draw_end after draw/autoview\n"
	"  --profile            enable libbrlobol internal timing counters\n"
	"  --screengrab <file>  write a screengrab after draw/autoview\n"
	"  --size <WxH>         display size for DM-backed runs (default: 512x512)\n"
	"  --repeat <N>         repeat full open/attach/draw/close cycle\n"
	"  --cache-dir <dir>    set BU_DIR_CACHE for this run\n"
	"  --clear-cache        clear --cache-dir before the first iteration\n"
	"  --lod-mesh <0|1>     run 'view lod mesh <0|1>' before draw\n"
	"  --lod-csg <0|1>      run 'view lod csg <0|1>' before draw\n"
	"  --lod-bot-threshold <N>\n"
	"                       run 'view lod bot_threshold <N>' before draw\n"
	"  --cad-compact <0|1>  compact eligible realized sources into CAD assemblies\n"
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
	if (BU_STR_EQUAL(arg, "--no-dm")) {
	    opts->attach_dm = false;
	    continue;
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
	if (BU_STR_EQUAL(arg, "--clear-cache")) {
	    opts->clear_cache = true;
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--dm")) {
	    if (++i >= argc)
		return false;
	    opts->dm_backend = argv[i];
	    continue;
	}
	if (BU_STR_EQUAL(arg, "--dm-name")) {
	    if (++i >= argc)
		return false;
	    opts->dm_name = argv[i];
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
	if (BU_STR_EQUAL(arg, "--cad-compact")) {
	    if (++i >= argc || !parse_bool_int(argv[i], &opts->cad_compact))
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

static int
apply_cad_options(struct ged *gedp, const struct options &opts)
{
    if (opts.cad_compact < 0)
	return BRLCAD_OK;

    SoBRLSceneController *scene =
	ged_draw_obol_scene_controller_ensure(gedp, 0);
    if (!scene)
	return BRLCAD_ERROR;

    scene->setCompactCadRealizationEnabled(opts.cad_compact ? TRUE : FALSE);
    return BRLCAD_OK;
}

static int
attach_display(struct ged *gedp, const struct options &opts)
{
    const char *av[5] = {
	"dm",
	"attach",
	opts.dm_backend.c_str(),
	opts.dm_name.c_str(),
	NULL
    };

    if (ged_exec_dm(gedp, 4, av) != BRLCAD_OK)
	return BRLCAD_ERROR;

    void *view_ctx = ged_view_active_ctx(gedp);
    struct dm *dmp = view_ctx ?
	(struct dm *)ged_view_context_display_manager_get(view_ctx) : NULL;
    if (!dmp)
	return BRLCAD_ERROR;

    dm_set_width(dmp, opts.width);
    dm_set_height(dmp, opts.height);
    dm_configure_win(dmp, 0);
    dm_set_zbuffer(dmp, 1);

    fastf_t windowbounds[6] = {-1, 1, -1, 1, -100, 100};
    dm_set_win_bounds(dmp, windowbounds);
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    dm_set_vp(dmp, bv_scale_storage_get(view));
    ged_view_context_display_manager_set(view_ctx, dmp);
    bv_dimensions_set(view, dm_get_width(dmp),
	dm_get_height(dmp));
    bv_unit_conversion_set(view,
	gedp->dbip->dbi_local2base,
	gedp->dbip->dbi_base2local);

    return BRLCAD_OK;
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
run_render(struct ged *gedp)
{
    void *view_ctx = ged_view_active_ctx(gedp);
    struct dm *dmp = view_ctx ?
	(struct dm *)ged_view_context_display_manager_get(view_ctx) : NULL;
    if (!dmp)
	return BRLCAD_ERROR;
    if (dm_draw_begin(dmp) != BRLCAD_OK)
	return BRLCAD_ERROR;
    return dm_draw_end(dmp);
}

static int
run_screengrab(struct ged *gedp, const char *filename)
{
    const char *av[3] = {"screengrab", filename, NULL};
    return ged_exec_screengrab(gedp, 2, av);
}

static int
run_once(const struct options &opts, int iter)
{
    struct timings t;
    struct BRLObolPerformanceCounters perf;
    brlobol_performance_counters_init(&perf);
    int attach_ret = BRLCAD_OK;
    int lod_ret = BRLCAD_OK;
    int cad_ret = BRLCAD_OK;
    int draw_ret = BRLCAD_OK;
    int autoview_ret = BRLCAD_OK;
    int render_ret = BRLCAD_OK;
    int screengrab_ret = BRLCAD_OK;
    int source_count_ret = 0;
    int profile_started = 0;
    size_t source_count = 0;

    auto start = std::chrono::steady_clock::now();
    struct ged *gedp = ged_open("db", opts.db_path.c_str(), 1);
    t.open_ms = elapsed_ms(start);
    if (!gedp) {
	std::fprintf(stderr, "ged_draw_perf: failed to open %s\n",
	    opts.db_path.c_str());
	return 1;
    }

    if (opts.attach_dm) {
	start = std::chrono::steady_clock::now();
	attach_ret = attach_display(gedp, opts);
	t.attach_ms = elapsed_ms(start);
    }

    if (attach_ret == BRLCAD_OK) {
	lod_ret = apply_lod_options(gedp, opts);
    }

    if (attach_ret == BRLCAD_OK && lod_ret == BRLCAD_OK) {
	cad_ret = apply_cad_options(gedp, opts);
    }

    if (attach_ret == BRLCAD_OK && lod_ret == BRLCAD_OK &&
	    cad_ret == BRLCAD_OK) {
	if (opts.profile) {
	    brlobol_performance_counters_reset();
	    brlobol_performance_counters_set_enabled(1);
	    profile_started = 1;
	}
	start = std::chrono::steady_clock::now();
	draw_ret = run_draw(gedp, opts);
	t.draw_ms = elapsed_ms(start);
    }

    if (draw_ret == BRLCAD_OK) {
	source_count_ret = ged_draw_obol_database_source_count(gedp, 0,
	    &source_count);
    }

    if (draw_ret == BRLCAD_OK && opts.autoview) {
	start = std::chrono::steady_clock::now();
	autoview_ret = run_autoview(gedp);
	t.autoview_ms = elapsed_ms(start);
    }

    if (draw_ret == BRLCAD_OK && autoview_ret == BRLCAD_OK && opts.render) {
	start = std::chrono::steady_clock::now();
	render_ret = run_render(gedp);
	t.render_ms = elapsed_ms(start);
    }

    if (draw_ret == BRLCAD_OK && autoview_ret == BRLCAD_OK &&
	    opts.attach_dm && !opts.screengrab.empty()) {
	start = std::chrono::steady_clock::now();
	screengrab_ret = run_screengrab(gedp, opts.screengrab.c_str());
	t.screengrab_ms = elapsed_ms(start);
    }

    if (profile_started) {
	brlobol_performance_counters_set_enabled(0);
	brlobol_performance_counters_get(&perf);
    }

    const char *result = ged_result(gedp);
    if (attach_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: dm attach failed: %s\n", result);
    if (lod_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: LoD setup failed: %s\n", result);
    if (cad_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: CAD setup failed: %s\n", result);
    if (draw_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: draw failed: %s\n", result);
    if (autoview_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: autoview failed: %s\n", result);
    if (render_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: render failed: %s\n", result);
    if (screengrab_ret != BRLCAD_OK)
	std::fprintf(stderr, "ged_draw_perf: screengrab failed: %s\n", result);

    std::printf("iter=%d open_ms=%.3f attach_ms=%.3f draw_ms=%.3f",
	iter, t.open_ms, t.attach_ms, t.draw_ms);
    if (opts.autoview)
	std::printf(" autoview_ms=%.3f", t.autoview_ms);
    if (opts.render)
	std::printf(" render_ms=%.3f", t.render_ms);
    if (!opts.screengrab.empty())
	std::printf(" screengrab_ms=%.3f", t.screengrab_ms);
    if (source_count_ret)
	std::printf(" obol_sources=%zu", source_count);
    if (perf.realize_calls || perf.source_replace_calls) {
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
    }
    std::printf(" status=%s\n",
	(attach_ret == BRLCAD_OK && draw_ret == BRLCAD_OK &&
	 lod_ret == BRLCAD_OK && cad_ret == BRLCAD_OK &&
	 autoview_ret == BRLCAD_OK && render_ret == BRLCAD_OK &&
	 screengrab_ret == BRLCAD_OK) ? "ok" : "error");

    int ret = (attach_ret == BRLCAD_OK && lod_ret == BRLCAD_OK &&
	cad_ret == BRLCAD_OK && draw_ret == BRLCAD_OK &&
	autoview_ret == BRLCAD_OK && render_ret == BRLCAD_OK &&
	screengrab_ret == BRLCAD_OK) ? 0 : 1;

    ged_close(gedp);
    return ret;
}

int
main(int argc, const char **argv)
{
    struct options opts;
    if (!parse_args(argc, argv, &opts)) {
	usage(stderr);
	return 1;
    }
    if (prepare_cache(opts) != BRLCAD_OK)
	return 1;

    std::printf("ged_draw_perf db=%s dm=%s draw_args=",
	opts.db_path.c_str(), opts.attach_dm ? opts.dm_backend.c_str() : "none");
    for (size_t i = 0; i < opts.draw_args.size(); i++)
	std::printf("%s%s", i ? " " : "", opts.draw_args[i].c_str());
    std::printf(" repeat=%d size=%dx%d autoview=%d render=%d profile=%d screengrab=%s",
	opts.repeat, opts.width, opts.height, opts.autoview ? 1 : 0,
	opts.render ? 1 : 0, opts.profile ? 1 : 0,
	opts.screengrab.empty() ? "none" : opts.screengrab.c_str());
    std::printf(" cache_dir=%s clear_cache=%d lod_mesh=%d lod_csg=%d lod_bot_threshold=%d\n",
	opts.cache_dir.empty() ? "ambient" : opts.cache_dir.c_str(),
	opts.clear_cache ? 1 : 0, opts.lod_mesh, opts.lod_csg,
	opts.lod_bot_threshold);
    std::printf("cad_compact=%d\n", opts.cad_compact);

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
