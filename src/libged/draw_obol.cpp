/*                  D R A W _ O B O L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file draw_obol.cpp
 *
 * libged bridge from neutral GED draw transactions to Obol scene state.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BADC.h"
#include "BObol/BDrawCache.h"
#include "BObol/BGrid.h"
#include "BObol/BInit.h"
#include "BObol/BImageSource.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BLodRealization.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSceneController.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BSnapAction.h"
#include "BObol/BViewportImage.h"
#include "BObol/BVListShape.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "BObol/BViewQuery.h"
#include "BObol/BViewStore.h"
#include "bg/line_layer.h"
#include "bg/plane.h"
#include "bg/polygon.h"
#include "bv.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/path.h"
#include "bu/str.h"
#include "bu/time.h"
#include "bu/units.h"
#include "bu/vls.h"
#include "ged.h"
#include "ged/db_index.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "icv.h"
#include "rt/db5.h"
#include "rt/db_fullpath.h"
#include "rt/search.h"
#include "rt/view.h"
#include "vmath.h"

#include "./ged_bobol_private.hpp"
#include "./draw_obol_bridge_private.hpp"
#include "./ged_draw_view_private.h"
#include "./ged_private.h"

#include <algorithm>
#include <Inventor/SbVec2f.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Obol/cad/SoCADAssembly.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <atomic>
#include <memory>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static const char ged_obol_group_intent_prefix[] = "ged-draw-group:";
static const char ged_obol_database_source_mode_key_marker[] =
    ":ged-draw-mode:";
static std::atomic<bool> ged_obol_test_force_face_set_failure(false);

extern "C" void
ged_draw_test_force_primitive_face_set_failure(int enable)
{
    ged_obol_test_force_face_set_failure.store(enable != 0,
	std::memory_order_relaxed);
}

extern "C" int
ged_draw_test_primitive_face_set_failure_enabled(void)
{
    return ged_obol_test_force_face_set_failure.load(
	std::memory_order_relaxed) ? 1 : 0;
}

static void ged_obol_transaction_observer(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    void *client_data);
static SbColor ged_obol_color_from_rgb(const unsigned char rgb[3]);
static float ged_obol_transparency_from_appearance_opacity(fastf_t opacity);
static fastf_t ged_obol_appearance_opacity_from_transparency(float transparency);
static int ged_obol_observer_ensure(struct ged *gedp,
				    struct ged_drawable *gdp);
static int ged_obol_bind_view_render_root(
    struct ged_view_context *view_ctx,
    BObolSceneController *shared_scene,
    BObolViewController *view_controller);
static int ged_obol_node_contains(SoNode *root, SoNode *target);
static size_t ged_obol_submit_region_bounds_async(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const std::unordered_set<std::string> &names,
    int draw_mode);
static int ged_obol_progressive_advance_provider(
    BObolViewController *controller,
    void *user_data,
    const BObolProgressiveOptions *options,
    BObolProgressiveStatus *status);
static int ged_obol_apply_draw_paths_to_scene(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char **paths,
    int path_count,
    const struct ged_draw_appearance_settings *settings,
    struct ged_draw_transaction_result *result,
    BObolSceneController *scene,
    uint32_t source_revision,
    int preserve_display_state);

/* Stage timing instrumentation.  Set BOBOL_DRAW_TIMING=1 to log where time is
 * spent across the progressive/coarse-first pipeline (draw command, structural
 * snapshot + region bounds, deferred realization job, per-pump deepening,
 * provider pump).  Off by default (single getenv, cached). */
static int
ged_obol_timing_enabled(void)
{
    static int enabled = -1;
    if (enabled < 0)
	enabled = getenv("BOBOL_DRAW_TIMING") ? 1 : 0;
    return enabled;
}

static void
ged_obol_timing_log(const char *label, int64_t start_us, long count)
{
    if (!ged_obol_timing_enabled())
	return;
    const double ms = (double)(bu_gettime() - start_us) / 1000.0;
    if (count >= 0)
	bu_log("[obol-timing] %-34s %8.1f ms  (n=%ld)\n", label, ms, count);
    else
	bu_log("[obol-timing] %-34s %8.1f ms\n", label, ms);
}

struct ged_obol_scoped_timer {
    const char *label;
    int64_t start;
    double minimumMilliseconds;
    explicit ged_obol_scoped_timer(const char *l,
	double minimumMs = 0.0) :
	label(l), start(bu_gettime()), minimumMilliseconds(minimumMs) {}
    ~ged_obol_scoped_timer()
    {
	if (minimumMilliseconds > 0.0 &&
	    (double)(bu_gettime() - start) / 1000.0 <
		minimumMilliseconds)
	    return;
	ged_obol_timing_log(label, start, -1);
    }
};

static struct bv *
ged_obol_bv(struct ged_view_context *view_ctx)
{
    return bv_context_view((struct bv_context *)view_ctx);
}

static const struct bv *
ged_obol_bv_const(const struct ged_view_context *view_ctx)
{
    return bv_context_view_const((const struct bv_context *)view_ctx);
}

struct ged_obol_deferred_realization_job;

struct ged_obol_progressive_provider_data {
    ged_obol_progressive_provider_data(void) :
	gedp(NULL),
	view_ctx(NULL),
	pending_autoview(0),
	expected_view_revision(0),
	autoview_factor(BV_AUTOVIEW_SCALE_DEFAULT),
	deferred_refine_stage(0)
    {
	deferred_appearance = ged_draw_appearance_settings
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
    }

    ~ged_obol_progressive_provider_data(void);

    struct ged *gedp;
    struct ged_view_context *view_ctx;
    int pending_autoview;
    uint64_t expected_view_revision;
    fastf_t autoview_factor;
    int deferred_refine_stage;
    struct ged_draw_appearance_settings deferred_appearance;
    std::vector<std::string> deferred_paths;
    std::shared_ptr<ged_obol_deferred_realization_job> deferred_job;
    /* Draw modes may be added while a prior root is still being realized.
     * Keep those independent jobs alive rather than leaving their root proxy
     * permanently visible. */
    std::vector<std::shared_ptr<ged_obol_deferred_realization_job>>
	pending_jobs;
    std::vector<std::shared_ptr<ged_obol_deferred_realization_job>>
	retired_jobs;
};

struct ged_obol_deferred_realization_item {
    ged_obol_deferred_realization_item(void) :
	source(NULL), database(NULL), primaryScene(FALSE), allowWireFallback(FALSE)
    {
    }

    ~ged_obol_deferred_realization_item(void)
    {
	if (source)
	    source->unref();
	if (database)
	    db_close(database);
	if (snapshotPath.getLength() > 0)
	    (void)bu_file_delete(snapshotPath.getString());
    }

    SoBRLDatabaseSource *source;
    struct db_i *database;
    SbString snapshotPath;
    std::string instanceKey;
    SbBool primaryScene;
    SbBool allowWireFallback;
    /* Worker->pump stream of completed per-leaf occurrences.  Shared with the
     * worker lambda so it outlives every push (the job join blocks on the
     * worker before items are destroyed). */
    std::shared_ptr<BObolCompactOccurrenceStream> stream;
};

struct ged_obol_deferred_realization_job {
    enum State {
	PENDING = 0,
	RUNNING = 1,
	COMPLETE = 2,
	FAILED = 3,
	CANCELLED = 4
    };

    ged_obol_deferred_realization_job(void) :
	state(PENDING), cancelRequested(false)
    {
    }

    ~ged_obol_deferred_realization_job(void)
    {
	cancelRequested.store(true, std::memory_order_release);
	if (worker.joinable())
	    worker.join();
    }

    void cancel(void)
    {
	cancelRequested.store(true, std::memory_order_release);
	for (const std::unique_ptr<ged_obol_deferred_realization_item> &item :
	     items) {
	    if (item && item->stream)
		item->stream->requestCancel();
	}
    }

    bool start(void)
    {
	if (items.empty() || state.load(std::memory_order_acquire) != PENDING)
	    return false;
	try {
	    worker = std::thread([this]() {
		state.store(RUNNING, std::memory_order_release);
		bool success = true;
		for (const std::unique_ptr<ged_obol_deferred_realization_item> &item :
		     items) {
		    if (cancelRequested.load(std::memory_order_acquire)) {
			success = false;
			break;
		    }
		    SoBRLDatabaseSource *source = item ? item->source : NULL;
		    if (!source) {
			success = false;
			break;
		    }
		    const int representation = source->representationMode.getValue();
		    const bool mesh =
			(source->realizationRoleFlags.getValue() &
			 SoBRLDatabaseSource::REALIZATION_ROLE_MESH) ||
			source->drawMode.getValue() ==
			SoBRLDatabaseSource::SHADED ||
			representation == SoBRLDatabaseSource::REPRESENTATION_SHADED ||
			representation == SoBRLDatabaseSource::REPRESENTATION_SHADED_BOTS ||
			representation == SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE;
		    const int64_t realize_start = bu_gettime();
		    BObolCompactOccurrenceStream *stream =
			item->stream ? item->stream.get() : NULL;
		    SbBool realized = mesh ? source->realizeDatabaseMesh(stream) :
			source->realizeDatabaseWireframe(stream);
		    if (!realized && mesh && item->allowWireFallback)
			realized = source->realizeDatabaseWireframe(stream);
		    ged_obol_timing_log(mesh ?
			"job: realizeDatabaseMesh (worker)" :
			"job: realizeDatabaseWireframe (worker)",
			realize_start, -1);
		    if (!realized &&
			!cancelRequested.load(std::memory_order_acquire)) {
			bu_log("Obol deferred realization failed for '%s': %s\n",
			    item->instanceKey.c_str(),
			    source->realizationDiagnostic.getValue().getString());
			success = false;
			break;
		    }
		    if (!realized) {
			success = false;
			break;
		    }
		}
		if (cancelRequested.load(std::memory_order_acquire))
		    state.store(CANCELLED, std::memory_order_release);
		else
		    state.store(success ? COMPLETE : FAILED,
			std::memory_order_release);
	    });
	} catch (...) {
	    state.store(FAILED, std::memory_order_release);
	    return false;
	}
	return true;
    }

    void join(void)
    {
	if (worker.joinable())
	    worker.join();
    }

    std::atomic<int> state;
    std::atomic<bool> cancelRequested;
    std::vector<std::unique_ptr<ged_obol_deferred_realization_item>> items;
    std::thread worker;
};

ged_obol_progressive_provider_data::~ged_obol_progressive_provider_data(void)
{
    if (deferred_job)
	deferred_job->cancel();
    for (const std::shared_ptr<ged_obol_deferred_realization_job> &job :
	 pending_jobs) {
	if (job)
	    job->cancel();
    }
    for (const std::shared_ptr<ged_obol_deferred_realization_job> &job :
	 retired_jobs) {
	if (job)
	    job->cancel();
    }
    deferred_job.reset();
    pending_jobs.clear();
    retired_jobs.clear();
}

class ged_obol_scene_mutation_batch_scope
{
public:
    ged_obol_scene_mutation_batch_scope(BObolSceneController *scene,
					size_t expected_database_sources = 0,
					size_t expected_groups = 0) :
	m_scene(scene),
	m_active(scene != NULL)
    {
	if (m_active)
	    m_scene->beginSceneMutationBatch(expected_database_sources,
					     expected_groups);
    }

    ~ged_obol_scene_mutation_batch_scope(void)
    {
	if (m_active)
	    m_scene->endSceneMutationBatch();
    }

    ged_obol_scene_mutation_batch_scope(
	const ged_obol_scene_mutation_batch_scope &) = delete;
    ged_obol_scene_mutation_batch_scope &operator=(
	const ged_obol_scene_mutation_batch_scope &) = delete;

private:
    BObolSceneController *m_scene;
    bool m_active;
};

struct ged_obol_source_state {
    ged_obol_source_state(void) :
	valid(false),
	viewMatched(false),
	groupValid(false),
	groupDrawMode(BOBOL_LOD_DRAW_WIRE),
	groupVisible(true),
	groupOverlay(false),
	groupTransparency(0.0f),
	sourceRevision(0),
	inputsRevision(0),
	drawMode(SoBRLDatabaseSource::WIREFRAME),
	visible(true),
	selected(false),
	highlighted(false),
	lineStyle(0),
	lineWidth(0),
	transparency(0.0f),
	colorOverride(false),
	color(1.0f, 1.0f, 1.0f),
	materialColorValid(false),
	materialColor(1.0f, 1.0f, 1.0f),
	materialRevision(0),
	drawMatrixValid(false),
	drawMatrix(SbMatrix::identity())
    {
    }

    bool valid;
    bool viewMatched;
    bool groupValid;
    std::string groupPath;
    int groupDrawMode;
    bool groupVisible;
    bool groupOverlay;
    float groupTransparency;
    uint32_t sourceRevision;
    uint32_t inputsRevision;
    int drawMode;
    bool visible;
    bool selected;
    bool highlighted;
    int lineStyle;
    int lineWidth;
    float transparency;
    bool colorOverride;
    SbColor color;
    bool materialColorValid;
    SbColor materialColor;
    uint32_t materialRevision;
    bool drawMatrixValid;
    SbMatrix drawMatrix;
};

static void
ged_obol_source_state_apply_appearance(
    ged_obol_source_state &state,
    const struct ged_draw_appearance_settings *settings)
{
    if (!settings)
	return;

    state.lineWidth = settings->s_line_width;
    state.transparency =
	ged_obol_transparency_from_appearance_opacity(settings->transparency);
    state.colorOverride = settings->color_override ? true : false;
    state.color = ged_obol_color_from_rgb(settings->color);
    if (state.colorOverride) {
	state.materialColorValid = true;
	state.materialColor = state.color;
    }
    if (state.groupValid)
	state.groupTransparency = state.transparency;
}

static SbMatrix
ged_obol_sbmatrix_from_mat(const mat_t mat);

/* The sole per-GED Obol owner.  Per-view controllers and their LoD/progressive
 * services belong to display endpoints and are always derived from the GED
 * view registry; they are never mirrored here. */
struct ged_obol_state {
    ged_obol_state(void) :
	observer_token(0),
	shared_controller(NULL),
	full_sync(0)
    {
    }

    ~ged_obol_state(void)
    {
	delete shared_controller;
    }

    ged_draw_observer_token observer_token;
    BObolViewController *shared_controller;
    int full_sync;
};

static struct ged_drawable *
ged_obol_gdp(struct ged *gedp)
{
    if (!gedp || !gedp->i)
	return NULL;
    return gedp->i->ged_gdp;
}

static ged_obol_state *
ged_obol_state_get(struct ged_drawable *gdp, int create)
{
    if (!gdp)
	return NULL;
    if (!gdp->gd_obol_state && create)
	gdp->gd_obol_state = new (std::nothrow) ged_obol_state;
    return static_cast<ged_obol_state *>(gdp->gd_obol_state);
}

static void
ged_obol_append_unique_path(std::vector<std::string> &paths, const char *path)
{
    if (!path || !path[0])
	return;
    std::string spath(path);
    if (std::find(paths.begin(), paths.end(), spath) == paths.end())
	paths.push_back(spath);
}

static void
ged_obol_append_unique_paths_from_words(std::vector<std::string> &paths,
					const char *words)
{
    if (!words || !words[0])
	return;

    std::string names(words);
    size_t pos = 0;
    while (pos < names.size()) {
	pos = names.find_first_not_of(" \t\r\n", pos);
	if (pos == std::string::npos)
	    break;
	size_t end = names.find_first_of(" \t\r\n", pos);
	std::string path = names.substr(pos,
					(end == std::string::npos) ? std::string::npos : end - pos);
	ged_obol_append_unique_path(paths, path.c_str());
	if (end == std::string::npos)
	    break;
	pos = end + 1;
    }
}

static std::vector<std::string>
ged_obol_transaction_paths(const struct ged_draw_transaction *txn,
			   const struct ged_draw_transaction_result *result)
{
    std::vector<std::string> paths;
    if (result && BU_VLS_IS_INITIALIZED(&result->names) &&
	bu_vls_strlen(&result->names)) {
	ged_obol_append_unique_paths_from_words(paths,
						bu_vls_cstr(&result->names));
    }

    if (txn && txn->path)
	ged_obol_append_unique_path(paths, txn->path);

    int path_count = (txn && txn->paths && txn->path_count > 0) ?
		     txn->path_count : 0;
    for (int i = 0; i < path_count; i++)
	ged_obol_append_unique_path(paths, txn->paths[i]);

    return paths;
}

static uint32_t
ged_obol_fold_revision(uint64_t revision)
{
    if (!revision)
	return 0;

    uint32_t folded = static_cast<uint32_t>(revision ^ (revision >> 32));
    return folded ? folded : 1;
}

static uint32_t
ged_obol_transaction_source_revision(
    const struct ged_draw_transaction_result *result)
{
    return ged_obol_fold_revision(result ? result->scene_revision_after : 0);
}

static int
ged_obol_database_draw_mode_from_ged(int mode)
{
    if (mode == GED_DRAW_MODE_SHADED ||
	mode == GED_DRAW_MODE_SHADED_BOTS)
	return SoBRLDatabaseSource::SHADED;
    return SoBRLDatabaseSource::WIREFRAME;
}

static int
ged_obol_database_draw_mode_to_ged(int mode)
{
    if (mode == SoBRLDatabaseSource::SHADED)
	return GED_DRAW_MODE_SHADED;
    return GED_DRAW_MODE_WIRE;
}

static int
ged_obol_database_representation_mode_from_ged(int mode)
{
    switch (mode) {
	case GED_DRAW_MODE_SHADED_BOTS:
	    return SoBRLDatabaseSource::REPRESENTATION_SHADED_BOTS;
	case GED_DRAW_MODE_SHADED:
	    return SoBRLDatabaseSource::REPRESENTATION_SHADED;
	case GED_DRAW_MODE_EVAL_WIRE:
	    return SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE;
	case GED_DRAW_MODE_HIDDEN_LINE:
	    return SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE;
	case GED_DRAW_MODE_EVAL_POINTS:
	    return SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS;
	case GED_DRAW_MODE_WIRE:
	default:
	    return SoBRLDatabaseSource::REPRESENTATION_WIRE;
    }
}

static int
ged_obol_draw_mode_uses_database_source_provider(int mode)
{
    switch (mode) {
	case GED_DRAW_MODE_WIRE:
	case GED_DRAW_MODE_SHADED_BOTS:
	case GED_DRAW_MODE_SHADED:
	case GED_DRAW_MODE_HIDDEN_LINE:
	case GED_DRAW_MODE_EVAL_WIRE:
	case GED_DRAW_MODE_EVAL_POINTS:
	    return 1;
	default:
	    return 0;
    }
}

static int
ged_obol_material_policy_from_ged(int policy)
{
    if (policy == GED_DRAW_OBOL_DATABASE_SOURCE_MATERIAL_DATABASE)
	return SoBRLDatabaseSource::MATERIAL_DATABASE;
    return SoBRLDatabaseSource::MATERIAL_INHERIT;
}

static int
ged_obol_material_policy_to_ged(int policy)
{
    if (policy == SoBRLDatabaseSource::MATERIAL_DATABASE)
	return GED_DRAW_OBOL_DATABASE_SOURCE_MATERIAL_DATABASE;
    return GED_DRAW_OBOL_DATABASE_SOURCE_MATERIAL_INHERIT;
}

static int
ged_obol_lod_draw_mode_from_ged(int mode)
{
    if (mode == GED_DRAW_MODE_SHADED_BOTS)
	return BOBOL_LOD_DRAW_SHADED_BOTS;
    if (mode == GED_DRAW_MODE_SHADED)
	return BOBOL_LOD_DRAW_SHADED;
    if (mode == GED_DRAW_MODE_HIDDEN_LINE)
	return BOBOL_LOD_DRAW_HIDDEN_LINE;
    if (mode == GED_DRAW_MODE_EVAL_POINTS)
	return BOBOL_LOD_DRAW_POINTS;
    return BOBOL_LOD_DRAW_WIRE;
}

static int
ged_obol_lod_draw_mode_to_ged(int mode)
{
    if (mode == BOBOL_LOD_DRAW_SHADED_BOTS)
	return GED_DRAW_MODE_SHADED_BOTS;
    if (mode == BOBOL_LOD_DRAW_SHADED)
	return GED_DRAW_MODE_SHADED;
    if (mode == BOBOL_LOD_DRAW_HIDDEN_LINE)
	return GED_DRAW_MODE_HIDDEN_LINE;
    if (mode == BOBOL_LOD_DRAW_POINTS)
	return GED_DRAW_MODE_EVAL_POINTS;
    return GED_DRAW_MODE_WIRE;
}

static int
ged_obol_lod_draw_mode_from_database_source(const SoBRLDatabaseSource *source)
{
    if (!source)
	return BOBOL_LOD_DRAW_WIRE;

    switch (source->representationMode.getValue()) {
	case SoBRLDatabaseSource::REPRESENTATION_SHADED_BOTS:
	    return BOBOL_LOD_DRAW_SHADED_BOTS;
	case SoBRLDatabaseSource::REPRESENTATION_SHADED:
	    return BOBOL_LOD_DRAW_SHADED;
	case SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE:
	    return BOBOL_LOD_DRAW_HIDDEN_LINE;
	case SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS:
	    return BOBOL_LOD_DRAW_POINTS;
	case SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE:
	case SoBRLDatabaseSource::REPRESENTATION_WIRE:
	    return BOBOL_LOD_DRAW_WIRE;
	default:
	    break;
    }

    return source->drawMode.getValue() == SoBRLDatabaseSource::SHADED ?
	   BOBOL_LOD_DRAW_SHADED : BOBOL_LOD_DRAW_WIRE;
}

static int
ged_obol_transaction_ged_draw_mode(struct ged *gedp,
				   const struct ged_draw_transaction *txn)
{
    int mode = -1;
    if (txn && txn->appearance) {
	const struct ged_draw_appearance_settings *appearance =
	    (const struct ged_draw_appearance_settings *)txn->appearance;
	mode = appearance->draw_mode;
    } else if (txn && txn->kind == GED_DRAW_TXN_DRAW && txn->mode >= 0) {
	mode = txn->mode;
    }
    if (mode < 0 && gedp)
	mode = ged_draw_default_mode(gedp);
    return mode < 0 ? GED_DRAW_MODE_WIRE : mode;
}

static int
ged_obol_transaction_is_completed_database_source_draw(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result)
{
    return (txn && result && result->status > 0 &&
	    txn->kind == GED_DRAW_TXN_DRAW &&
	    ged_obol_draw_mode_uses_database_source_provider(
		ged_obol_transaction_ged_draw_mode(gedp, txn))) ? 1 : 0;
}

static int
ged_obol_transaction_defer_leaf_expansion(
    const struct ged_draw_transaction *txn)
{
    if (!txn || txn->kind != GED_DRAW_TXN_DRAW || !txn->appearance)
	return 0;

    const struct ged_draw_appearance_settings *appearance =
	(const struct ged_draw_appearance_settings *)txn->appearance;
    return appearance->defer_leaf_expansion ? 1 : 0;
}

static int
ged_obol_drawn_path_mode(struct ged *gedp, struct ged_view_context *view_ctx,
			 const char *path)
{
    if (ged_draw_path_state(gedp, view_ctx, path,
			    GED_DRAW_MODE_HIDDEN_LINE) > 0)
	return GED_DRAW_MODE_HIDDEN_LINE;
    if (ged_draw_path_state(gedp, view_ctx, path,
			    GED_DRAW_MODE_SHADED_BOTS) > 0 ||
	ged_draw_path_state(gedp, view_ctx, path,
			    GED_DRAW_MODE_SHADED) > 0)
	return GED_DRAW_MODE_SHADED;
    return GED_DRAW_MODE_WIRE;
}

static const char *
ged_obol_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static std::string
ged_obol_publish_leaf_name_from_path(const char *path)
{
    const char *name = ged_obol_skip_leading_slash(path);
    if (!name || !name[0])
	return std::string();

    const char *slash = strrchr(name, '/');
    std::string leaf((slash && slash[1]) ? slash + 1 : name);
    const size_t instance_specifier = leaf.find('@');
    if (instance_specifier != std::string::npos)
	leaf.erase(instance_specifier);
    return leaf;
}

static int
ged_obol_path_can_publish_database_source(struct db_i *dbip, const char *path)
{
    if (!dbip || !path || !path[0])
	return 0;

    const std::string leaf_name = ged_obol_publish_leaf_name_from_path(path);
    if (leaf_name.empty())
	return 0;

    struct directory *dp = db_lookup(dbip, leaf_name.c_str(), LOOKUP_QUIET);
    if (!dp || dp->d_addr == RT_DIR_PHONY_ADDR)
	return 0;

    if (dp->d_flags & RT_DIR_COMB)
	return 1;
    return (dp->d_flags & RT_DIR_SOLID) ? 1 : 0;
}

static int
ged_obol_path_equal(const char *a, const char *b)
{
    if (!a || !b)
	return 0;
    if (BU_STR_EQUAL(a, b))
	return 1;
    return BU_STR_EQUAL(ged_obol_skip_leading_slash(a),
			ged_obol_skip_leading_slash(b));
}

static int
ged_obol_view_scope_is_independent(struct ged_view_context *view_ctx)
{
    return view_ctx && ged_view_context_is_independent(view_ctx);
}

std::string
ged_obol_view_scope_name(struct ged_view_context *view_ctx)
{
    if (!view_ctx)
	return "shared";

    const char *name = bv_name_get(ged_obol_bv_const(view_ctx));
    if (name && name[0])
	return std::string(name);

    char fallback[64] = {0};
    snprintf(fallback, sizeof(fallback), "%p", (void *)view_ctx);
    return std::string(fallback);
}

static std::string
ged_obol_database_source_instance_key(struct ged_view_context *view_ctx, const char *path)
{
    if (!path || !path[0])
	return std::string();

    if (!ged_obol_view_scope_is_independent(view_ctx)) {
	std::string key("/");
	key += ged_obol_skip_leading_slash(path);
	return key;
    }

    std::string key("ged-view:");
    key += ged_obol_view_scope_name(view_ctx);
    key += ":";
    key += ged_obol_skip_leading_slash(path);
    return key;
}

static std::string
ged_obol_database_source_instance_key_for_mode(
    struct ged_view_context *view_ctx,
    const char *path,
    int draw_mode)
{
    std::string key = ged_obol_database_source_instance_key(view_ctx, path);
    if (key.empty() || draw_mode < 0 ||
	draw_mode == GED_DRAW_MODE_WIRE ||
	ged_obol_view_scope_is_independent(view_ctx))
	return key;

    char mode_buf[64] = {0};
    snprintf(mode_buf, sizeof(mode_buf), "%s%d",
	     ged_obol_database_source_mode_key_marker, draw_mode);
    key += mode_buf;
    return key;
}

static int
ged_obol_database_source_mode_from_instance_key(const char *instance_key)
{
    if (!instance_key || !instance_key[0])
	return -1;

    const std::string key(instance_key);
    const size_t marker_pos =
	key.rfind(ged_obol_database_source_mode_key_marker);
    if (marker_pos == std::string::npos)
	return -1;

    const char *mode_str = instance_key + marker_pos +
			   strlen(ged_obol_database_source_mode_key_marker);
    if (!mode_str[0])
	return -1;

    char *endptr = NULL;
    long mode = strtol(mode_str, &endptr, 10);
    if (endptr == mode_str || (endptr && *endptr != '\0') ||
	mode < 0 || mode > INT_MAX)
	return -1;

    return static_cast<int>(mode);
}

static std::string
ged_obol_database_source_base_instance_key(const char *instance_key)
{
    if (!instance_key)
	return std::string();

    const std::string key(instance_key);
    const size_t marker_pos =
	key.rfind(ged_obol_database_source_mode_key_marker);
    if (marker_pos == std::string::npos)
	return key;

    return key.substr(0, marker_pos);
}

static int
ged_obol_database_source_summary_ged_mode(
    const BObolDatabaseSourceSummary &summary)
{
    if (summary.representationMode >= 0)
	return summary.representationMode;

    const int keyed_mode =
	ged_obol_database_source_mode_from_instance_key(
	    summary.instanceKey.getString());
    if (keyed_mode >= 0)
	return keyed_mode;

    return ged_obol_database_draw_mode_to_ged(summary.drawMode);
}

static int
ged_obol_database_source_summary_matches_mode(
    const BObolDatabaseSourceSummary &summary,
    int draw_mode)
{
    return draw_mode < 0 ||
	   ged_obol_database_source_summary_ged_mode(summary) == draw_mode;
}

static std::string
ged_obol_database_source_instance_prefix(struct ged_view_context *view_ctx)
{
    if (!ged_obol_view_scope_is_independent(view_ctx))
	return std::string();

    std::string prefix("ged-view:");
    prefix += ged_obol_view_scope_name(view_ctx);
    prefix += ":";
    return prefix;
}

static int
ged_obol_database_source_instance_in_scope(
    const BObolDatabaseSourceSummary &summary,
    struct ged_view_context *view_ctx)
{
    if (!ged_obol_view_scope_is_independent(view_ctx)) {
	if (summary.instanceKey.getLength() == 0)
	    return 1;
	const std::string base_instance_key =
	    ged_obol_database_source_base_instance_key(
		summary.instanceKey.getString());
	if (base_instance_key.compare(0, 9, "ged-view:") == 0)
	    return 0;
	/* Shared occurrence keys carry hierarchy and Boolean identity and are
	 * intentionally not path strings.  Only an explicit view prefix makes
	 * a source view-local. */
	return 1;
    }

    const std::string prefix =
	ged_obol_database_source_instance_prefix(view_ctx);
    const char *key = summary.instanceKey.getString();
    return key && bu_strncmp(key, prefix.c_str(), prefix.size()) == 0;
}

static std::string
ged_obol_database_source_owner_group_path_from_summary(
    const BObolDatabaseSourceSummary &summary)
{
    const char *parent_path = summary.parentGroupPath.getString();
    const char *parent_norm = ged_obol_skip_leading_slash(parent_path);
    if (parent_norm && parent_norm[0])
	return std::string(parent_norm);

    return std::string(ged_obol_skip_leading_slash(
			   summary.path.getString()));
}

static int
ged_obol_shape_record_matches_path(const struct ged_draw_shape_record *record,
				   const char *path)
{
    if (!record || !path || !path[0])
	return 0;
    if (ged_obol_path_equal(record->display_name, path) ||
	ged_obol_path_equal(record->leaf_name, path))
	return 1;

    if (!record->fullpath || record->fullpath->fp_len <= 0)
	return 0;

    char *record_path = db_path_to_string(record->fullpath);
    if (!record_path)
	return 0;
    int matched = ged_obol_path_equal(record_path, path);
    bu_free(record_path, "GED Obol draw sync record path");
    return matched;
}

static void
ged_obol_source_state_add_group(ged_obol_source_state &state,
				struct ged *gedp,
				ged_draw_group_ref group)
{
    if (!gedp || ged_draw_group_ref_is_null(group))
	return;

    struct ged_draw_group_record group_record;
    if (!ged_draw_group_record_get(gedp, group, &group_record) ||
	!group_record.path || !group_record.path[0])
	return;

    state.groupValid = true;
    state.groupPath = group_record.path;
    state.groupDrawMode = ged_obol_lod_draw_mode_from_ged(
			      group_record.draw_mode);
    state.groupVisible = group_record.visible ? true : false;
    state.groupOverlay = group_record.is_overlay ? true : false;
    state.groupTransparency =
	static_cast<float>(group_record.transparency);
    if (group_record.appearance.color_override) {
	state.colorOverride = true;
	state.color = ged_obol_color_from_rgb(group_record.appearance.color);
	state.materialColorValid = true;
	state.materialColor = state.color;
    }
}

static void
ged_obol_source_state_from_record(ged_obol_source_state &state,
				  struct ged *gedp,
				  const struct ged_draw_shape_record *record,
				  int view_matched)
{
    if (!record)
	return;

    state.valid = true;
    state.viewMatched = view_matched ? true : false;
    state.sourceRevision = ged_obol_fold_revision(record->source_revision);
    state.inputsRevision = ged_obol_fold_revision(record->inputs_revision);
    state.drawMode = ged_obol_database_draw_mode_from_ged(record->draw_mode);
    state.visible = record->visible ? true : false;
    state.selected = record->selected ? true : false;
    state.highlighted = record->highlighted ? true : false;
    state.lineWidth = record->line_width;
    state.transparency = static_cast<float>(record->transparency);
    ged_obol_source_state_add_group(state, gedp, record->group);

    if ((!state.materialColorValid || state.materialRevision == 0) &&
	!state.colorOverride && gedp && gedp->dbip &&
	record->fullpath && record->fullpath->fp_len > 0) {
	SbColor db_material_color;
	if (bobol_database_source_fullpath_material_color(
		gedp->dbip, record->fullpath, db_material_color)) {
	    state.materialColorValid = true;
	    state.materialColor = db_material_color;
	}
    }
}

static int
ged_obol_source_state_from_database_source(ged_obol_source_state &state,
	struct ged *gedp,
	const char *path,
	int draw_mode)
{
    if (!gedp || !path || !path[0])
	return 0;

    struct ged_draw_database_source_summary source_summary;
    memset(&source_summary, 0, sizeof(source_summary));
    struct ged_draw_scene_display_summary display_summary;
    memset(&display_summary, 0, sizeof(display_summary));
    const int have_source =
	ged_draw_obol_database_source_summary_for_path(gedp, path,
	    &source_summary) && source_summary.valid;
    const int have_display =
	ged_draw_obol_database_source_display_summary_for_path(gedp, path,
	    &display_summary) && display_summary.valid;
    if (!have_source && !have_display)
	return 0;

    state.valid = true;
    state.viewMatched = true;
    state.drawMode = ged_obol_database_draw_mode_from_ged(draw_mode);
    if (have_source) {
	state.sourceRevision =
	    ged_obol_fold_revision(source_summary.source_revision);
	state.inputsRevision =
	    ged_obol_fold_revision(source_summary.inputs_revision);
    }
    if (have_display) {
	state.visible = display_summary.visible ? true : false;
	state.selected = display_summary.selected ? true : false;
	state.highlighted = display_summary.highlighted ? true : false;
	state.lineStyle = display_summary.line_style;
	state.lineWidth = display_summary.line_width;
	state.transparency =
	    static_cast<float>(display_summary.transparency);
	state.materialColorValid =
	    display_summary.material_valid ? true : false;
	if (state.materialColorValid) {
	    state.materialColor = SbColor(
				      static_cast<float>(display_summary.material_color[0]) / 255.0f,
				      static_cast<float>(display_summary.material_color[1]) / 255.0f,
				      static_cast<float>(display_summary.material_color[2]) / 255.0f);
	}
    }

    struct ged_draw_shape_material_summary material_summary;
    if (ged_draw_obol_database_source_material_summary_for_path(gedp, path,
	    &material_summary) && material_summary.valid)
	state.materialRevision =
	    ged_obol_fold_revision(material_summary.material_revision);

    struct ged_draw_obol_draw_state_summary draw_state;
    memset(&draw_state, 0, sizeof(draw_state));
    if (ged_draw_obol_database_source_draw_state_for_path(gedp, path,
	    &draw_state) && draw_state.valid && draw_state.draw_mat_valid) {
	state.drawMatrixValid = true;
	state.drawMatrix = ged_obol_sbmatrix_from_mat(draw_state.draw_mat);
    }

    if ((!state.materialColorValid || state.materialRevision == 0) &&
	!state.colorOverride) {
	SbColor db_material_color;
	if (gedp && bobol_database_source_path_material_color(gedp->dbip, path,
		db_material_color)) {
	    state.materialColorValid = true;
	    state.materialColor = db_material_color;
	}
    }

    struct bu_vls owner_group_path = BU_VLS_INIT_ZERO;
    if (ged_draw_obol_database_source_owner_group_path_for_path(gedp, path,
	    &owner_group_path) && bu_vls_strlen(&owner_group_path) > 0) {
	const char *group_path = bu_vls_cstr(&owner_group_path);
	state.groupValid = true;
	state.groupPath = group_path;
	state.groupDrawMode = ged_obol_lod_draw_mode_from_ged(draw_mode);
	state.groupVisible = true;
	state.groupOverlay = false;
	state.groupTransparency = state.transparency;
	struct ged_draw_group_record_summary group_summary;
	memset(&group_summary, 0, sizeof(group_summary));
	if (ged_draw_obol_group_record_summary_for_path(gedp, group_path,
		&group_summary)) {
	    state.groupDrawMode =
		ged_obol_lod_draw_mode_from_ged(group_summary.draw_mode);
	    state.groupVisible = group_summary.visible ? true : false;
	    state.groupOverlay = group_summary.is_overlay ? true : false;
	    state.groupTransparency =
		static_cast<float>(group_summary.transparency);
	}
    }
    bu_vls_free(&owner_group_path);

    return 1;
}

static void
ged_obol_source_state_resolve_database_material(ged_obol_source_state &state,
	struct db_i *dbip,
	const char *path)
{
    if (state.materialColorValid || state.colorOverride || !dbip || !path ||
	!path[0])
	return;

    SbColor db_material_color;
    if (bobol_database_source_path_material_color(dbip, path,
	    db_material_color)) {
	state.materialColorValid = true;
	state.materialColor = db_material_color;
    }
}

struct ged_obol_find_source_state_context {
    struct ged *gedp;
    struct ged_view_context *viewCtx;
    const char *path;
    int mode;
    ged_obol_source_state state;
};

static int
ged_obol_find_source_state_cb(const struct ged_draw_shape_record *record,
			      void *userdata)
{
    ged_obol_find_source_state_context *ctx =
	static_cast<ged_obol_find_source_state_context *>(userdata);
    if (!ctx || !record || !ged_obol_shape_record_matches_path(record,
	    ctx->path))
	return 1;

    struct ged_draw_group_record group_record;
    if (!ged_draw_group_record_get(ctx->gedp, record->group,
				   &group_record))
	return 1;

    const int view_matched =
	ged_draw_group_record_in_view(&group_record, ctx->viewCtx) ? 1 : 0;
    if (!ctx->state.valid || view_matched || !ctx->state.viewMatched)
	ged_obol_source_state_from_record(ctx->state, ctx->gedp, record,
					  view_matched);

    return view_matched ? 0 : 1;
}

static ged_obol_source_state
ged_obol_find_source_state(struct ged *gedp,
			   struct ged_view_context *view_ctx,
			   const char *path,
			   int draw_mode)
{
    if (!gedp || !path || !path[0])
	return ged_obol_source_state();

    ged_obol_source_state direct_state;
    if (ged_obol_source_state_from_database_source(direct_state, gedp, path,
	    draw_mode))
	return direct_state;

    if (draw_mode >= 0 &&
	ged_draw_path_state(gedp, view_ctx, path, draw_mode) <= 0)
	return ged_obol_source_state();

    ged_obol_find_source_state_context ctx;
    ctx.gedp = gedp;
    ctx.viewCtx = view_ctx;
    ctx.path = path;
    ctx.mode = draw_mode;
    ged_draw_foreach_shape_record(gedp, ged_obol_find_source_state_cb, &ctx);
    return ctx.state;
}

static int
ged_obol_intent_is_ged_draw_group(const SbString &intent)
{
    const char *value = intent.getString();
    if (!value)
	return 0;
    return bu_strncmp(value, ged_obol_group_intent_prefix,
		   sizeof(ged_obol_group_intent_prefix) - 1) == 0;
}

static std::string
ged_obol_group_intent_path(const char *group_path)
{
    std::string intent(ged_obol_group_intent_prefix);
    if (group_path)
	intent += ged_obol_skip_leading_slash(group_path);
    return intent;
}

static std::string
ged_obol_group_path_from_record_path(const char *path)
{
    if (!path)
	return std::string();

    const size_t prefix_len = sizeof(ged_obol_group_intent_prefix) - 1;
    if (bu_strncmp(path, ged_obol_group_intent_prefix, prefix_len) == 0)
	path += prefix_len;
    return std::string(ged_obol_skip_leading_slash(path));
}


static std::string
ged_obol_top_group_path_from_record_path(const char *path)
{
    std::string group_path = ged_obol_group_path_from_record_path(path);
    const size_t slash = group_path.find('/');
    if (slash == std::string::npos)
	return group_path;
    return group_path.substr(0, slash);
}

static const char *
ged_obol_group_record_path(const SoBRLSceneGroup *scene_group)
{
    if (!scene_group)
	return NULL;

    if (!scene_group->drawIntentValid.getValue())
	return NULL;

    const char *path = scene_group->drawIntentPath.getValue().getString();
    if (!path)
	return NULL;

    const size_t prefix_len = sizeof(ged_obol_group_intent_prefix) - 1;
    if (bu_strncmp(path, ged_obol_group_intent_prefix, prefix_len) == 0)
	path += prefix_len;
    return ged_obol_skip_leading_slash(path);
}

static bool
ged_obol_group_parent_leaf(const std::string &path,
			   std::string &parent,
			   std::string &leaf)
{
    if (path.empty())
	return false;

    const size_t slash = path.find_last_of('/');
    if (slash == std::string::npos) {
	parent.clear();
	leaf = path;
    } else {
	parent = path.substr(0, slash);
	leaf = path.substr(slash + 1);
    }
    return !leaf.empty();
}

static int
ged_obol_sync_group_state(BObolSceneController *scene,
			  const ged_obol_source_state &state,
			  const char *source_instance_key)
{
    if (!scene || !state.groupValid || state.groupPath.empty())
	return 0;

    int changed = 0;
    SoGroup *existing_group = scene->findGroup(state.groupPath.c_str());
    SoGroup *group = scene->ensureGroup(state.groupPath.c_str());
    if (!group)
	return 0;
    if (!existing_group)
	changed = 1;

    const std::string intent_path =
	ged_obol_group_intent_path(state.groupPath.c_str());
    if (scene->setGroupDrawIntent(state.groupPath.c_str(),
				  intent_path.c_str(),
				  state.groupDrawMode,
				  BOBOL_LOD_DRAW_WIRE,
				  state.groupOverlay ? TRUE : FALSE,
				  state.sourceRevision) > 0)
	changed = 1;

    const SbBool group_color_override =
	state.colorOverride ? TRUE : FALSE;
    const SbColor group_color = state.colorOverride ?
				state.color : SbColor(1.0f, 1.0f, 1.0f);
    const SbBool group_material_valid =
	(state.colorOverride && state.materialColorValid) ? TRUE : FALSE;
    const SbColor group_material = group_material_valid ?
				   state.materialColor : SbColor(1.0f, 1.0f, 1.0f);
    if (scene->setGroupDisplayState(state.groupPath.c_str(),
				    state.groupVisible ? TRUE : FALSE,
				    FALSE,
				    FALSE,
				    state.lineStyle,
				    state.lineWidth,
				    state.groupTransparency,
				    group_color_override,
				    group_color,
				    group_material_valid,
				    group_material,
				    state.materialRevision) > 0)
	changed = 1;

    if (source_instance_key && source_instance_key[0]) {
	int source_already_grouped = 0;
	BObolDatabaseSourceSummary summary;
	if (scene->getDatabaseSourceSummaryForInstance(source_instance_key,
		summary) && summary.valid) {
	    const std::string owner_group_path =
		ged_obol_database_source_owner_group_path_from_summary(summary);
	    source_already_grouped =
		ged_obol_path_equal(owner_group_path.c_str(),
				    state.groupPath.c_str());
	}
	if (!source_already_grouped &&
	    scene->moveDatabaseSourceInstanceToGroup(source_instance_key,
		    state.groupPath.c_str()) > 0)
	    changed = 1;
    }
    return changed;
}

static int
ged_obol_prune_empty_groups(BObolSceneController *scene)
{
    if (!scene)
	return 0;

    int changed = 0;
    int pass_removed = 1;
    while (pass_removed) {
	pass_removed = 0;
	std::vector<std::string> prune_paths;
	const int summary_count = scene->getSceneTreeSummaryCount();
	for (int i = summary_count - 1; i >= 0; i--) {
	    BObolSceneTreeSummary tree_summary;
	    if (!scene->getSceneTreeSummary(i, tree_summary) ||
		!tree_summary.valid ||
		tree_summary.nodeKind !=
		BObolSceneTreeSummary::NODE_GROUP ||
		tree_summary.childCount != 0 ||
		tree_summary.path.getLength() == 0 ||
		BU_STR_EQUAL(tree_summary.path.getString(), "/"))
		continue;

	    SoGroup *group = scene->findGroup(tree_summary.path.getString());
	    if (group && group->isOfType(SoBRLSceneGroup::getClassTypeId()) &&
		static_cast<SoBRLSceneGroup *>(group)->
		overlayIntent.getValue())
		continue;

	    prune_paths.push_back(tree_summary.path.getString());
	}
	for (const std::string &path : prune_paths) {
	    if (scene->removeGroup(path.c_str()) > 0) {
		pass_removed = 1;
		changed = 1;
	    }
	}
    }
    return changed;
}

static int
ged_obol_path_has_prefix(const char *path, const char *prefix)
{
    if (!path || !path[0] || !prefix || !prefix[0])
	return 0;

    path = ged_obol_skip_leading_slash(path);
    prefix = ged_obol_skip_leading_slash(prefix);
    const size_t prefix_len = strlen(prefix);
    if (prefix_len == 0)
	return 0;
    if (bu_strncmp(path, prefix, prefix_len) != 0)
	return 0;
    return path[prefix_len] == '\0' || path[prefix_len] == '/';
}

static int
ged_obol_path_has_component_name(const char *path,
				 const char *name,
				 size_t first_idx)
{
    if (!path || !name)
	return 0;

    path = ged_obol_skip_leading_slash(path);
    name = ged_obol_skip_leading_slash(name);
    const size_t name_len = strlen(name);
    if (!name_len)
	return 0;

    size_t idx = 0;
    const char *cursor = path;
    while (*cursor) {
	while (*cursor == '/')
	    cursor++;
	if (!*cursor)
	    break;
	const char *slash = strchr(cursor, '/');
	size_t component_len = slash ?
			       static_cast<size_t>(slash - cursor) : strlen(cursor);
	if (idx >= first_idx && component_len == name_len &&
	    bu_strncmp(cursor, name, component_len) == 0)
	    return 1;
	if (!slash)
	    break;
	cursor = slash + 1;
	idx++;
    }
    return 0;
}

static std::string
ged_obol_normalized_path_string(const char *path)
{
    if (!path || !path[0])
	return std::string();

    std::string normalized(ged_obol_skip_leading_slash(path));
    while (!normalized.empty() && normalized.back() == '/')
	normalized.pop_back();
    return normalized;
}

static std::string
ged_obol_strip_comb_instance_suffix(const std::string &component)
{
    const size_t at_pos = component.rfind('@');
    if (at_pos == std::string::npos || at_pos == 0 ||
	    at_pos + 1 >= component.size())
	return component;

    for (size_t i = at_pos + 1; i < component.size(); i++) {
	if (component[i] < '0' || component[i] > '9')
	    return component;
    }
    return component.substr(0, at_pos);
}

static std::string
ged_obol_semantic_path_string(const char *path)
{
    const char *cursor = ged_obol_skip_leading_slash(path);
    if (!cursor || !cursor[0])
	return std::string();

    std::string semantic;
    while (*cursor) {
	while (*cursor == '/')
	    cursor++;
	if (!*cursor)
	    break;
	const char *slash = strchr(cursor, '/');
	const size_t component_len = slash ?
				     static_cast<size_t>(slash - cursor) : strlen(cursor);
	std::string component;
	component.reserve(component_len);
	for (size_t i = 0; i < component_len; i++)
	    component.push_back(cursor[i]);
	if (!semantic.empty())
	    semantic += "/";
	semantic += ged_obol_strip_comb_instance_suffix(component);
	if (!slash)
	    break;
	cursor = slash + 1;
    }
    return semantic;
}

static int
ged_obol_path_has_semantic_prefix(const char *path, const char *prefix)
{
    const std::string semantic_path = ged_obol_semantic_path_string(path);
    const std::string semantic_prefix =
	ged_obol_semantic_path_string(prefix);
    if (semantic_path.empty() || semantic_prefix.empty())
	return 0;

    const size_t prefix_len = semantic_prefix.size();
    if (semantic_path.compare(0, prefix_len, semantic_prefix) != 0)
	return 0;
    return semantic_path[prefix_len] == '\0' ||
	   semantic_path[prefix_len] == '/';
}

static int
ged_obol_path_has_semantic_component_name(const char *path,
					  const char *name,
					  size_t first_idx)
{
    if (!path || !name)
	return 0;

    const std::string semantic_name = ged_obol_semantic_path_string(name);
    if (semantic_name.empty() || semantic_name.find('/') != std::string::npos)
	return 0;

    const char *cursor = ged_obol_skip_leading_slash(path);
    size_t idx = 0;
    while (cursor && *cursor) {
	while (*cursor == '/')
	    cursor++;
	if (!*cursor)
	    break;
	const char *slash = strchr(cursor, '/');
	const size_t component_len = slash ?
				     static_cast<size_t>(slash - cursor) : strlen(cursor);
	std::string component(cursor, component_len);
	if (idx >= first_idx &&
		ged_obol_strip_comb_instance_suffix(component) ==
		semantic_name)
	    return 1;
	if (!slash)
	    break;
	cursor = slash + 1;
	idx++;
    }
    return 0;
}

static std::vector<std::string>
ged_obol_unshadowed_source_paths(const std::vector<std::string> &paths)
{
    std::vector<std::string> filtered;
    if (paths.empty())
	return filtered;

    std::vector<std::pair<std::string, size_t>> ordered;
    ordered.reserve(paths.size());
    for (size_t i = 0; i < paths.size(); i++) {
	std::string normalized =
	    ged_obol_normalized_path_string(paths[i].c_str());
	if (!normalized.empty())
	    ordered.push_back(std::make_pair(normalized, i));
    }
    std::sort(ordered.begin(), ordered.end(),
	      [](const std::pair<std::string, size_t> &a,
    const std::pair<std::string, size_t> &b) {
	if (a.first != b.first)
	    return a.first < b.first;
	return a.second < b.second;
    });

    filtered.reserve(paths.size());
    std::unordered_set<std::string> emitted;
    emitted.reserve(paths.size());
    for (const std::pair<std::string, size_t> &entry : ordered) {
	const std::string descendant_prefix = entry.first + "/";
	auto descendant = std::lower_bound(ordered.begin(), ordered.end(),
					   descendant_prefix,
					   [](const std::pair<std::string, size_t> &candidate,
	const std::string &prefix) {
	    return candidate.first < prefix;
	});
	if (descendant != ordered.end() &&
	    descendant->first.compare(0, descendant_prefix.size(),
				      descendant_prefix) == 0)
	    continue;
	if (emitted.insert(entry.first).second)
	    filtered.push_back(paths[entry.second]);
    }
    return filtered;
}

static void
ged_obol_remove_shadowed_source_paths(std::vector<std::string> &paths)
{
    if (paths.size() < 2)
	return;

    paths = ged_obol_unshadowed_source_paths(paths);
}

struct ged_obol_drawn_source_path_mode {
    std::string path;
    int mode;
};

static void
ged_obol_append_unique_path_mode(
    std::vector<ged_obol_drawn_source_path_mode> &path_modes,
    const char *path,
    int mode)
{
    if (!path || !path[0])
	return;

    for (const ged_obol_drawn_source_path_mode &entry : path_modes) {
	if (entry.mode == mode &&
	    ged_obol_path_equal(entry.path.c_str(), path))
	    return;
    }

    ged_obol_drawn_source_path_mode entry;
    entry.path = path;
    entry.mode = mode;
    path_modes.push_back(entry);
}

static void
ged_obol_remove_shadowed_source_path_modes(
    std::vector<ged_obol_drawn_source_path_mode> &path_modes)
{
    if (path_modes.size() < 2)
	return;

    std::vector<std::string> paths;
    for (const ged_obol_drawn_source_path_mode &entry : path_modes)
	ged_obol_append_unique_path(paths, entry.path.c_str());

    std::vector<ged_obol_drawn_source_path_mode> filtered;
    std::vector<std::string> unshadowed =
	ged_obol_unshadowed_source_paths(paths);
    std::unordered_set<std::string> keep;
    keep.reserve(unshadowed.size());
    for (const std::string &path : unshadowed)
	keep.insert(ged_obol_normalized_path_string(path.c_str()));

    for (const ged_obol_drawn_source_path_mode &entry : path_modes) {
	if (keep.find(ged_obol_normalized_path_string(entry.path.c_str())) !=
	    keep.end())
	    ged_obol_append_unique_path_mode(filtered, entry.path.c_str(),
					     entry.mode);
    }
    path_modes.swap(filtered);
}

static int
ged_obol_record_matches_any_target(const struct ged_draw_shape_record *record,
				   const char *record_path,
				   const std::vector<std::string> *targets)
{
    if (!targets || targets->empty())
	return 1;
    if (!record_path || !record_path[0])
	return 0;

    for (const std::string &target : *targets) {
	if (ged_obol_path_has_prefix(record_path, target.c_str()) ||
	    ged_obol_shape_record_matches_path(record, target.c_str()) ||
	    ged_obol_path_has_component_name(record_path, target.c_str(), 0))
	    return 1;
    }

    return 0;
}

static int
ged_obol_group_summary_in_view(
    const struct ged_draw_group_record_summary *group_summary,
    struct ged_view_context *view_ctx)
{
    if (!group_summary)
	return 0;
    if (!view_ctx)
	return !group_summary->in_view_scope;
    if (ged_view_context_is_independent(view_ctx))
	return group_summary->in_view_scope &&
	       group_summary->view_ctx == view_ctx;
    if (!group_summary->in_view_scope)
	return 1;
    return group_summary->view_ctx == view_ctx;
}

struct ged_obol_drawn_source_path_ctx {
    struct ged *gedp;
    struct ged_view_context *viewCtx;
    int mode;
    const std::vector<std::string> *targets;
    std::vector<std::string> paths;
    std::vector<ged_obol_drawn_source_path_mode> pathModes;
    std::unordered_map<std::string, int> groupVisible;
    std::unordered_set<std::string> pathKeys;
    std::unordered_set<std::string> pathModeKeys;
};

static void
ged_obol_drawn_source_append_path(
    ged_obol_drawn_source_path_ctx *ctx,
    const char *path)
{
    if (!ctx || !path || !path[0])
	return;

    const std::string normalized = ged_obol_normalized_path_string(path);
    if (normalized.empty() || !ctx->pathKeys.insert(normalized).second)
	return;
    ctx->paths.push_back(path);
}

static void
ged_obol_drawn_source_append_path_mode(
    ged_obol_drawn_source_path_ctx *ctx,
    const char *path,
    int mode)
{
    if (!ctx || !path || !path[0])
	return;

    const std::string normalized = ged_obol_normalized_path_string(path);
    if (normalized.empty())
	return;

    std::string key = normalized;
    key += '\037';
    key += std::to_string(mode);
    if (!ctx->pathModeKeys.insert(key).second)
	return;

    ged_obol_drawn_source_path_mode entry;
    entry.path = path;
    entry.mode = mode;
    ctx->pathModes.push_back(entry);
}

static int
ged_obol_drawn_source_group_visible(
    ged_obol_drawn_source_path_ctx *ctx,
    const BObolDatabaseSourceSummary &summary)
{
    if (!ctx)
	return 0;

    const std::string owner_group_path =
	ged_obol_database_source_owner_group_path_from_summary(summary);
    if (owner_group_path.empty())
	return 1;

    auto cached = ctx->groupVisible.find(owner_group_path);
    if (cached != ctx->groupVisible.end())
	return cached->second;

    int visible = 1;
    struct ged_draw_group_record_summary group_summary;
    memset(&group_summary, 0, sizeof(group_summary));
    if (ged_draw_obol_group_record_summary_for_path(ctx->gedp,
	    owner_group_path.c_str(), &group_summary)) {
	visible = (!group_summary.is_overlay &&
		   group_summary.visible &&
		   ged_obol_group_summary_in_view(&group_summary,
		       ctx->viewCtx)) ? 1 : 0;
    }

    ctx->groupVisible[owner_group_path] = visible;
    return visible;
}

static void
ged_obol_drawn_source_summary_append(
    ged_obol_drawn_source_path_ctx *ctx,
    const BObolDatabaseSourceSummary &summary)
{
    if (!ctx || !summary.valid || !summary.visible ||
	!ged_obol_database_source_instance_in_scope(summary, ctx->viewCtx) ||
	!ged_obol_database_source_summary_matches_mode(summary, ctx->mode) ||
	!ged_obol_drawn_source_group_visible(ctx, summary))
	return;

    const char *source_path = summary.path.getString();
    if (!source_path || !source_path[0] ||
	!ged_obol_record_matches_any_target(NULL, source_path, ctx->targets))
	return;

    const int source_mode = ged_obol_database_source_summary_ged_mode(summary);
    ged_obol_drawn_source_append_path_mode(ctx, source_path, source_mode);
    ged_obol_drawn_source_append_path(ctx, source_path);
}

static void
ged_obol_collect_drawn_source_paths(ged_obol_drawn_source_path_ctx *ctx)
{
    if (!ctx || !ctx->gedp)
	return;

    BObolSceneController *scene = ged_draw_obol_scene_controller(ctx->gedp);
    if (!scene)
	return;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;
	ged_obol_drawn_source_summary_append(ctx, summary);
    }
}

static std::vector<ged_obol_drawn_source_path_mode>
ged_obol_drawn_source_path_modes(struct ged *gedp,
				 struct ged_view_context *view_ctx,
				 int mode,
				 const std::vector<std::string> *targets)
{
    std::vector<ged_obol_drawn_source_path_mode> path_modes;
    if (!gedp)
	return path_modes;

    ged_obol_drawn_source_path_ctx ctx;
    ctx.gedp = gedp;
    ctx.viewCtx = view_ctx;
    ctx.mode = mode;
    ctx.targets = targets;
    ged_obol_collect_drawn_source_paths(&ctx);
    path_modes.swap(ctx.pathModes);
    ged_obol_remove_shadowed_source_path_modes(path_modes);
    return path_modes;
}

static std::vector<std::string>
ged_obol_drawn_source_paths(struct ged *gedp,
			    struct ged_view_context *view_ctx,
			    int mode,
			    const std::vector<std::string> *targets)
{
    std::vector<std::string> paths;
    if (!gedp)
	return paths;

    ged_obol_drawn_source_path_ctx ctx;
    ctx.gedp = gedp;
    ctx.viewCtx = view_ctx;
    ctx.mode = mode;
    ctx.targets = targets;
    ged_obol_collect_drawn_source_paths(&ctx);
    paths.swap(ctx.paths);
    ged_obol_remove_shadowed_source_paths(paths);
    return paths;
}

static int
ged_obol_direct_draw_path_modes(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const std::vector<ged_obol_drawn_source_path_mode> &path_modes,
    uint32_t source_revision,
    BObolSceneController *scene,
    int preserve_display_state = 0)
{
    if (!gedp || !view_ctx || path_modes.empty() || !scene)
	return 0;

    std::unordered_map<int, std::vector<std::string>> paths_by_mode;
    for (const ged_obol_drawn_source_path_mode &entry : path_modes) {
	const int mode = entry.mode >= 0 ? entry.mode : GED_DRAW_MODE_WIRE;
	ged_obol_append_unique_path(paths_by_mode[mode], entry.path.c_str());
    }

    int changed = 0;
    for (auto &mode_paths : paths_by_mode) {
	if (mode_paths.second.empty())
	    continue;
	struct ged_draw_appearance_settings settings =
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
	settings.draw_mode = mode_paths.first;
	settings.mixed_modes = 1;

	std::vector<const char *> path_av;
	path_av.reserve(mode_paths.second.size());
	for (const std::string &path : mode_paths.second)
	    path_av.push_back(path.c_str());

	const int ret = ged_obol_apply_draw_paths_to_scene(gedp, view_ctx,
		path_av.data(), static_cast<int>(path_av.size()),
		&settings, NULL, scene, source_revision,
		preserve_display_state);
	if (ret < 0)
	    return -1;
	if (ret > 0)
	    changed = 1;
    }

    return changed;
}

static std::vector<std::string>
ged_obol_shadowed_target_source_paths(
    const std::vector<std::string> &targets,
    const std::vector<std::string> &source_paths)
{
    std::vector<std::string> shadowed;
    for (const std::string &target : targets) {
	for (const std::string &source_path : source_paths) {
	    if (!ged_obol_path_equal(source_path.c_str(), target.c_str()) &&
		ged_obol_path_has_prefix(source_path.c_str(),
					 target.c_str())) {
		ged_obol_append_unique_path(shadowed, target.c_str());
		break;
	    }
	}
    }
    return shadowed;
}

static int
ged_obol_database_source_mark_published_current(BObolSceneController *scene,
	SoBRLDatabaseSource *source);

struct ged_obol_view_lod_policy_state {
    ged_obol_view_lod_policy_state(void) :
	valid(false),
	viewDependent(false),
	csgLodEnabled(false),
	meshLodEnabled(false),
	viewScale(0.0f),
	lodScale(1.0f),
	viewWidth(0),
	viewHeight(0),
	botThreshold(0),
	curveScale(0.0f),
	pointScale(0.0f),
	meshEnabled(false)
    {
    }

    bool valid;
    bool viewDependent;
    bool csgLodEnabled;
    bool meshLodEnabled;
    float viewScale;
    float lodScale;
    int viewWidth;
    int viewHeight;
    uint32_t botThreshold;
    float curveScale;
    float pointScale;
    bool meshEnabled;
};


struct ged_obol_preserved_source_display_state {
    std::string instanceKey;
    SbBool visible;
    SbBool selected;
    SbBool highlighted;
    int lineStyle;
    int lineWidth;
    float transparency;
    SbBool colorOverride;
    SbColor color;
    SbBool materialColorValid;
    SbColor materialColor;
    uint32_t materialRevision;
};


static int
ged_obol_source_path_matches_any_target(
    const char *source_path,
    const std::vector<std::string> &targets)
{
    if (!source_path || !source_path[0])
	return 0;

    for (const std::string &target : targets) {
	if (ged_obol_path_equal(source_path, target.c_str()) ||
	    ged_obol_path_has_prefix(source_path, target.c_str()) ||
	    ged_obol_path_has_semantic_prefix(source_path, target.c_str()))
	    return 1;
    }
    return 0;
}


static std::vector<ged_obol_preserved_source_display_state>
ged_obol_preserve_source_display_states(
    BObolSceneController *scene,
    struct ged_view_context *view_ctx,
    const std::vector<std::string> &targets,
    int draw_mode)
{
    std::vector<ged_obol_preserved_source_display_state> states;
    if (!scene || targets.empty())
	return states;

    const int source_count = scene->getDatabaseSourceCount();
    states.reserve(static_cast<size_t>(source_count));
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary,
		    draw_mode) ||
	    !ged_obol_source_path_matches_any_target(
		summary.path.getString(), targets))
	    continue;

	ged_obol_preserved_source_display_state state;
	state.instanceKey = summary.instanceKey.getString();
	state.visible = summary.visible;
	state.selected = summary.selected;
	state.highlighted = summary.highlighted;
	state.lineStyle = summary.lineStyle;
	state.lineWidth = summary.lineWidth;
	state.transparency = summary.transparency;
	state.colorOverride = summary.colorOverride;
	state.color = summary.color;
	state.materialColorValid = summary.materialColorValid;
	state.materialColor = summary.materialColor;
	state.materialRevision = summary.materialRevision;
	states.push_back(state);
    }
    return states;
}


static int
ged_obol_restore_source_display_states(
    BObolSceneController *scene,
    const std::vector<ged_obol_preserved_source_display_state> &states)
{
    if (!scene || states.empty())
	return 0;

    int changed = 0;
    for (const ged_obol_preserved_source_display_state &state : states) {
	SoBRLDatabaseSource *source =
	    scene->findDatabaseSourceInstance(state.instanceKey.c_str());
	if (!source)
	    continue;

	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    continue;

	const int source_changed = scene->setDatabaseSourceInstanceState(
	    state.instanceKey.c_str(),
	    TRUE,
	    summary.sourceRevision,
	    summary.inputsRevision,
	    state.visible,
	    state.selected,
	    state.highlighted,
	    state.lineStyle,
	    state.lineWidth,
	    state.transparency,
	    state.colorOverride,
	    state.color,
	    state.materialColorValid,
	    state.materialColor,
	    state.materialRevision);
	if (source_changed > 0)
	    changed = 1;
    }
    return changed;
}


struct ged_obol_publish_placement_state {
    ged_obol_publish_placement_state(void) :
	valid(false),
	drawMatrixValid(false),
	drawMatrix(SbMatrix::identity()),
	drawCenterValid(false),
	drawCenter(0.0f, 0.0f, 0.0f),
	drawSizeValid(false),
	drawSize(0.0f)
    {
    }

    bool valid;
    bool drawMatrixValid;
    SbMatrix drawMatrix;
    bool drawCenterValid;
    SbVec3f drawCenter;
    bool drawSizeValid;
    float drawSize;
};

static ged_obol_view_lod_policy_state
ged_obol_view_lod_policy_state_for_source(struct ged *gedp, struct ged_view_context *view_ctx)
{
    ged_obol_view_lod_policy_state state;
    if (!gedp)
	return state;

    struct ged_view_context *policy_view = view_ctx ? view_ctx :
	ged_view_active_ctx(gedp);
    if (!policy_view)
	return state;

    ged_draw_source_lod_policy policy;
    if (!ged_draw_source_lod_policy_get(&policy, policy_view))
	return state;

    state.valid = true;
    state.meshEnabled = policy.mesh_enabled ? true : false;
    state.csgLodEnabled = policy.csg_enabled ? true : false;
    state.meshLodEnabled = policy.mesh_enabled ? true : false;
    if (policy.mesh_enabled)
	state.botThreshold = policy.bot_threshold > 0 ?
			     (uint32_t)std::min(policy.bot_threshold,
						static_cast<size_t>(UINT32_MAX)) : 1;
    else if (policy.csg_enabled && policy.bot_threshold > 0)
	state.botThreshold =
	    (uint32_t)std::min(policy.bot_threshold,
			       static_cast<size_t>(UINT32_MAX));
    state.viewDependent =
	(policy.csg_enabled || policy.mesh_enabled) ? true : false;
    state.viewScale =
	static_cast<float>(bv_scale_get(ged_obol_bv_const(policy_view)));
    state.lodScale = static_cast<float>(policy.scale);
    state.viewWidth = bv_width_get(ged_obol_bv_const(policy_view));
    state.viewHeight = bv_height_get(ged_obol_bv_const(policy_view));
    state.curveScale = static_cast<float>(policy.curve_scale);
    state.pointScale = static_cast<float>(policy.point_scale);
    return state;
}

static int
ged_obol_apply_view_lod_policy(struct ged *gedp,
			       struct ged_view_context *view_ctx,
			       BObolSceneController *scene,
			       const char *source_instance_key)
{
    if (!gedp || !scene || !source_instance_key || !source_instance_key[0])
	return 0;

    ged_obol_view_lod_policy_state policy_state =
	ged_obol_view_lod_policy_state_for_source(gedp, view_ctx);
    if (!policy_state.valid)
	return 0;

    return scene->setDatabaseSourceInstanceRealizationViewPolicy(
	       source_instance_key,
	       policy_state.viewDependent ? TRUE : FALSE,
	       policy_state.csgLodEnabled ? TRUE : FALSE,
	       policy_state.meshLodEnabled ? TRUE : FALSE,
	       policy_state.viewScale,
	       policy_state.lodScale,
	       policy_state.viewWidth,
	       policy_state.viewHeight,
	       policy_state.botThreshold,
	       policy_state.curveScale,
	       policy_state.pointScale);
}

static int
ged_obol_replace_path(struct ged *gedp,
		      struct ged_view_context *view_ctx,
		      struct db_i *dbip,
		      const char *path,
		      int draw_mode,
		      uint32_t source_revision,
		      BObolSceneController *scene,
		      int use_retained_state = 1,
		      int preserve_existing_revision = 0,
		      const struct ged_draw_appearance_settings *appearance_settings = NULL,
		      const ged_obol_publish_placement_state *placement_state = NULL,
		      const char *target_group_path = NULL,
		      const char *source_instance_key_override = NULL,
		      const struct BObolDrawMetadataRecord *source_metadata = NULL,
		      const struct ged_bobol_publication_context *publication = NULL)
{
    if (!dbip || !path || !path[0] || !scene)
	return 0;
    if (!ged_obol_path_can_publish_database_source(dbip, path))
	return 0;

    ged_obol_source_state source_state = use_retained_state ?
					 ged_obol_find_source_state(gedp, view_ctx, path, draw_mode) :
					 ged_obol_source_state();
    const bool retained_state_valid = source_state.valid;
    const uint32_t retained_source_revision = source_state.sourceRevision;
    const uint32_t retained_inputs_revision = source_state.inputsRevision;
    if (source_state.valid) {
	if (source_state.sourceRevision != 0)
	    source_revision = source_state.sourceRevision;
    }

    const std::string instance_key =
	(source_instance_key_override && source_instance_key_override[0]) ?
	std::string(source_instance_key_override) :
	ged_obol_database_source_instance_key_for_mode(view_ctx, path,
	    draw_mode);
    const std::string representation_key = instance_key;
    const int next_draw_mode = ged_obol_database_draw_mode_from_ged(
				   draw_mode);
    if (draw_mode >= 0 &&
	(!source_instance_key_override || !source_instance_key_override[0])) {
	const std::string base_instance_key =
	    ged_obol_database_source_instance_key(view_ctx, path);
	if (!base_instance_key.empty() &&
	    base_instance_key != instance_key) {
	    SoBRLDatabaseSource *base_source =
		scene->findDatabaseSourceInstance(base_instance_key.c_str());
	    BObolDatabaseSourceSummary base_summary;
	    if (base_source && base_source->getSummary(base_summary) &&
		base_summary.valid &&
		base_summary.representationMode < 0) {
		if (!scene->findDatabaseSourceInstance(instance_key.c_str()) &&
		    base_summary.drawMode == next_draw_mode) {
		    uint32_t promoted_revision =
			base_summary.sourceRevision ?
			base_summary.sourceRevision : source_revision;
		    if (scene->renameDatabaseSourceInstance(
			    base_instance_key.c_str(), instance_key.c_str(),
			    path, promoted_revision) >= 0)
			(void)scene->setDatabaseSourceInstanceRepresentation(
			    instance_key.c_str(), representation_key.c_str(),
			    draw_mode);
		} else {
		    (void)scene->removeDatabaseSourceInstance(
			base_instance_key.c_str());
		}
	    }
	}
    }
    SoBRLDatabaseSource *existing_source =
	scene->findDatabaseSourceInstance(instance_key.c_str());
    bool preserve_external_current = false;
    bool replace_deferred_proxy = false;
    BObolDatabaseSourceSummary existing_summary;
    if (existing_source && existing_source->getSummary(existing_summary) &&
	existing_summary.valid) {
	const bool same_database_source =
	    preserve_existing_revision &&
	    !existing_summary.stale &&
	    existing_source->getDatabase() == dbip &&
	    ged_obol_path_equal(existing_summary.path.getString(), path) &&
	    BU_STR_EQUAL(existing_summary.instanceKey.getString(),
			 instance_key.c_str()) &&
	    BU_STR_EQUAL(existing_summary.representationKey.getString(),
			 representation_key.c_str()) &&
	    existing_summary.representationMode == draw_mode &&
	    existing_summary.drawMode == next_draw_mode;
	if (same_database_source)
	    source_revision = existing_summary.sourceRevision;
	if ((existing_summary.realizationRoleFlags &
	     SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL) &&
	    existing_summary.realizationStatus ==
	    SoBRLDatabaseSource::REALIZED &&
	    !existing_summary.stale) {
	    /* External primaries may be proxy line sets for shaded draws; preserve
	     * them by source lifetime, not by final draw-mode geometry type. */
	    const bool has_external_primary =
		existing_summary.realizedShapeCount > 0 ||
		existing_summary.realizedMeshCount > 0;
	    for (int i = 0; has_external_primary &&
		 i < existing_source->getRealizedShapeSummaryCount(); i++) {
		BObolRealizedShapeSummary shape;
		if (existing_source->getRealizedShapeSummary(i, shape) &&
		    shape.valid &&
		    BU_STR_EQUAL(shape.sourceType.getString(), "proxy") &&
		    BU_STR_EQUAL(shape.geometryKind.getString(), "aabb")) {
		    replace_deferred_proxy = true;
		    break;
		}
	    }
	    if (has_external_primary && !replace_deferred_proxy)
		preserve_external_current = true;
	}
	if (!source_state.valid) {
	    source_state.valid = true;
	    source_state.sourceRevision = existing_summary.sourceRevision;
	    source_state.inputsRevision = existing_summary.inputsRevision;
	    source_state.visible = existing_summary.visible ? true : false;
	    source_state.selected = existing_summary.selected ? true : false;
	    source_state.highlighted =
		existing_summary.highlighted ? true : false;
	    source_state.lineStyle = existing_summary.lineStyle;
	    source_state.lineWidth = existing_summary.lineWidth;
	    source_state.transparency =
		static_cast<float>(existing_summary.transparency);
	    source_state.colorOverride =
		existing_summary.colorOverride ? true : false;
	    source_state.color = existing_summary.color;
	    source_state.materialColorValid =
		existing_summary.materialColorValid ? true : false;
	    source_state.materialColor = existing_summary.materialColor;
	    source_state.materialRevision = existing_summary.materialRevision;
	    source_state.drawMatrixValid =
		existing_summary.drawMatrixValid ? true : false;
	    source_state.drawMatrix = existing_summary.drawMatrix;
	}
	if (retained_state_valid) {
	    source_state.sourceRevision = retained_source_revision;
	    source_state.inputsRevision = retained_inputs_revision;
	} else if (source_revision != 0) {
	    source_state.sourceRevision = source_revision;
	}
    }

    if (replace_deferred_proxy) {
	(void)existing_source->clearRealizedGeometry(TRUE);
	existing_source->markStale(SoBRLDatabaseSource::STALE_SOURCE);
    }
    const bool had_group_state = source_state.groupValid;
    const struct ged_draw_appearance_settings *publication_appearance =
	(publication && publication->active) ? publication->appearance : NULL;
    if (draw_mode >= 0) {
	source_state.valid = true;
	if (!publication_appearance) {
	    source_state.groupValid = true;
	}
	if (source_state.groupValid && source_state.groupPath.empty())
	    source_state.groupPath = ged_obol_skip_leading_slash(path);
	source_state.groupDrawMode =
	    ged_obol_lod_draw_mode_from_ged(draw_mode);
	if (source_state.groupValid && !had_group_state) {
	    source_state.groupVisible = true;
	    source_state.groupTransparency = source_state.transparency;
	}
	source_state.groupOverlay = false;
    }
    if (publication_appearance) {
	source_state.lineWidth =
	    publication_appearance->s_line_width;
	source_state.transparency =
	    ged_obol_transparency_from_appearance_opacity(
		publication_appearance->transparency);
	source_state.colorOverride =
	    publication_appearance->color_override ? true : false;
	source_state.color = ged_obol_color_from_rgb(
				 publication_appearance->color);
	if (source_state.colorOverride) {
	    source_state.materialColorValid = true;
	    source_state.materialColor = source_state.color;
	}
	if (source_state.groupValid)
	    source_state.groupTransparency = source_state.transparency;
    }
    ged_obol_source_state_apply_appearance(source_state, appearance_settings);
    if (!source_state.materialColorValid && !source_state.colorOverride &&
	source_metadata && source_metadata->directoryFound &&
	source_metadata->hasColor) {
	source_state.materialColorValid = true;
	source_state.materialColor = SbColor(
	    static_cast<float>(source_metadata->color[0]) / 255.0f,
	    static_cast<float>(source_metadata->color[1]) / 255.0f,
	    static_cast<float>(source_metadata->color[2]) / 255.0f);
    }
    if (publication && publication->active &&
	!publication->group_path.empty() &&
	(!source_state.groupValid || !retained_state_valid)) {
	source_state.groupValid = true;
	source_state.groupPath = publication->group_path;
	source_state.groupDrawMode =
	    ged_obol_lod_draw_mode_from_ged(draw_mode);
	source_state.groupVisible = true;
	source_state.groupOverlay = false;
	source_state.groupTransparency = source_state.transparency;
    }
    if (target_group_path && target_group_path[0]) {
	const std::string explicit_group_path =
	    ged_obol_group_path_from_record_path(target_group_path);
	if (!explicit_group_path.empty()) {
	    source_state.groupValid = true;
	    source_state.groupPath = explicit_group_path;
	    source_state.groupDrawMode =
		ged_obol_lod_draw_mode_from_ged(draw_mode);
	    source_state.groupVisible = true;
	    source_state.groupOverlay = false;
	    source_state.groupTransparency = source_state.transparency;
	}
    }

	/* Tree metadata only supersedes the database material when it provides a
	 * color.  Root-backed primitive draws have directory metadata but no tree
	 * material, and must still resolve their effective database color. */
    if (!source_metadata || !source_metadata->hasColor)
	ged_obol_source_state_resolve_database_material(source_state, dbip, path);

    ged_obol_view_lod_policy_state policy_state =
	ged_obol_view_lod_policy_state_for_source(gedp, view_ctx);

    BObolDatabaseSourcePublishState publish_state;
    publish_state.sourceInstanceKey = instance_key.c_str();
    publish_state.sourcePath = path;
    publish_state.sourceRepresentationKey = representation_key.c_str();
    publish_state.targetGroupPath =
	(source_state.groupValid && !source_state.groupPath.empty()) ?
	source_state.groupPath.c_str() : NULL;
    publish_state.database = dbip;
    publish_state.drawMode = next_draw_mode;
    publish_state.representationMode = draw_mode;
    publish_state.sourceRevisionValid = source_revision != 0 ? TRUE : FALSE;
    publish_state.sourceRevision = source_revision;
    if (source_state.sourceRevision != 0) {
	publish_state.sourceRevisionValid = TRUE;
	publish_state.sourceRevision = source_state.sourceRevision;
    }
    publish_state.inputsRevision = source_state.inputsRevision;
    publish_state.visible = source_state.visible ? TRUE : FALSE;
    publish_state.selected = source_state.selected ? TRUE : FALSE;
    publish_state.highlighted = source_state.highlighted ? TRUE : FALSE;
    publish_state.lineStyle = source_state.lineStyle;
    publish_state.lineWidth = source_state.lineWidth;
    publish_state.transparency = source_state.transparency;
    publish_state.colorOverride = source_state.colorOverride ? TRUE : FALSE;
    publish_state.color = source_state.color;
    publish_state.materialColorValid =
	source_state.materialColorValid ? TRUE : FALSE;
    publish_state.materialColor = source_state.materialColor;
    publish_state.materialRevision = source_state.materialRevision;
    if (!preserve_external_current) {
	publish_state.roleFlagsValid = TRUE;
	publish_state.roleFlags = SoBRLDatabaseSource::REALIZATION_ROLE_NONE;
	/*
	 * A wire presentation is still a mesh realization when mesh LoD is
	 * active.  The compact mesh path preserves plotted wire geometry for
	 * non-BoTs, while BoTs use the same producer-authored PoP prefix as
	 * shaded drawing.  Leaving this as a legacy wire realization bypasses
	 * PoP submission entirely.
	 */
	if (draw_mode == GED_DRAW_MODE_EVAL_POINTS ||
	    (draw_mode == GED_DRAW_MODE_WIRE &&
	     policy_state.valid && policy_state.viewDependent &&
	     policy_state.meshLodEnabled))
	    publish_state.roleFlags |=
		SoBRLDatabaseSource::REALIZATION_ROLE_MESH;
    }
    if (policy_state.valid && policy_state.viewDependent) {
	publish_state.viewPolicyValid = TRUE;
	publish_state.viewDependent =
	    policy_state.viewDependent ? TRUE : FALSE;
	publish_state.csgLodEnabled =
	    policy_state.csgLodEnabled ? TRUE : FALSE;
	publish_state.meshLodEnabled =
	    policy_state.meshLodEnabled ? TRUE : FALSE;
	publish_state.viewScale = policy_state.viewScale;
	publish_state.lodScale = policy_state.lodScale;
	publish_state.viewWidth = policy_state.viewWidth;
	publish_state.viewHeight = policy_state.viewHeight;
	publish_state.botThreshold = policy_state.botThreshold;
	publish_state.curveScale = policy_state.curveScale;
	publish_state.pointScale = policy_state.pointScale;
    }
    if (placement_state && placement_state->valid) {
	publish_state.placementValid = TRUE;
	publish_state.drawMatrixValid =
	    placement_state->drawMatrixValid ? TRUE : FALSE;
	publish_state.drawMatrix = placement_state->drawMatrix;
	publish_state.drawCenterValid =
	    placement_state->drawCenterValid ? TRUE : FALSE;
	publish_state.drawCenter = placement_state->drawCenter;
	publish_state.drawSizeValid =
	    placement_state->drawSizeValid ? TRUE : FALSE;
	publish_state.drawSize = placement_state->drawSize;
    } else if (source_state.drawMatrixValid) {
	publish_state.placementValid = TRUE;
	publish_state.drawMatrixValid = TRUE;
	publish_state.drawMatrix = source_state.drawMatrix;
    }

    int changed = scene->publishDatabaseSourceInstance(publish_state);
    if (changed < 0)
	return changed;

    /* Database-source draws preserve material resolved for individual tree
     * occurrences.  This also keeps independently synchronized views from
     * applying an aggregate root color to every compact instance. */
    const int material_policy_changed =
	scene->setDatabaseSourceInstanceMaterialPolicy(instance_key.c_str(),
	    SoBRLDatabaseSource::MATERIAL_DATABASE);
    if (material_policy_changed < 0)
	return material_policy_changed;
    if (material_policy_changed > 0)
	changed = 1;
    if (ged_obol_sync_group_state(scene, source_state,
				  instance_key.c_str()))
	changed = 1;
    if (preserve_external_current) {
	SoBRLDatabaseSource *source =
	    scene->findDatabaseSourceInstance(instance_key.c_str());
	BObolDatabaseSourceSummary summary;
	if (source && source->getSummary(summary) && summary.valid &&
	    (summary.realizedShapeCount > 0 || summary.realizedMeshCount > 0) &&
	    ged_obol_database_source_mark_published_current(scene, source))
	    changed = 1;
    }
    return changed;
}

static int
ged_obol_replace_path_and_realize(struct ged *gedp,
				  struct ged_view_context *view_ctx,
				  struct db_i *dbip,
				  const char *path,
				  int draw_mode,
				  uint32_t source_revision,
				  BObolSceneController *scene,
				  int use_retained_state = 1,
				  int preserve_existing_revision = 0,
				  const struct ged_draw_appearance_settings *appearance_settings = NULL)
{
    const int changed = ged_obol_replace_path(gedp, view_ctx, dbip, path,
			draw_mode, source_revision, scene,
			use_retained_state, preserve_existing_revision,
			appearance_settings);
    if (changed < 0)
	return changed;

    if ((draw_mode == GED_DRAW_MODE_EVAL_WIRE ||
	 draw_mode == GED_DRAW_MODE_EVAL_POINTS) &&
	ged_draw_obol_database_source_realize_for_path(gedp, path))
	return 1;

    return changed;
}

static int
ged_obol_replace_paths(struct db_i *dbip,
		       const std::vector<std::string> &paths,
		       int draw_mode,
		       uint32_t source_revision,
		       struct ged *gedp,
		       struct ged_view_context *view_ctx,
		       BObolSceneController *scene,
		       int use_retained_state = 1,
		       int preserve_existing_revision = 0,
		       const struct ged_draw_appearance_settings *appearance_settings = NULL)
{
    if (!dbip || paths.empty() || !scene)
	return 0;

    ged_obol_scene_mutation_batch_scope batch(scene, paths.size(),
	    paths.size());
    int changed = 0;
    for (const std::string &path : paths) {
	if (ged_obol_replace_path_and_realize(gedp, view_ctx, dbip,
					      path.c_str(),
					      draw_mode, source_revision, scene,
					      use_retained_state,
					      preserve_existing_revision,
					      appearance_settings) > 0)
	    changed = 1;
    }

    if (changed)
	scene->realizePending();
    return changed;
}

static int
ged_obol_replace_path_modes(struct db_i *dbip,
			    const std::vector<ged_obol_drawn_source_path_mode> &path_modes,
			    uint32_t source_revision,
			    struct ged *gedp,
			    struct ged_view_context *view_ctx,
			    BObolSceneController *scene,
			    int use_retained_state = 1,
			    int preserve_existing_revision = 0)
{
    if (!dbip || path_modes.empty() || !scene)
	return 0;

    ged_obol_scene_mutation_batch_scope batch(scene, path_modes.size(),
	    path_modes.size());
    int changed = 0;
    for (const ged_obol_drawn_source_path_mode &entry : path_modes) {
	if (ged_obol_replace_path_and_realize(gedp, view_ctx, dbip,
					      entry.path.c_str(), entry.mode, source_revision,
					      scene, use_retained_state,
					      preserve_existing_revision) > 0)
	    changed = 1;
    }

    if (changed)
	scene->realizePending();
    return changed;
}

static void
ged_obol_append_database_source_instance_key(
    std::vector<std::string> &instance_keys,
    const BObolDatabaseSourceSummary &summary);

static int
ged_obol_append_exact_database_source_instance_keys(
    std::vector<std::string> &instance_keys,
    BObolSceneController *scene,
    struct ged_view_context *view_ctx,
    const char *path,
    int draw_mode)
{
    if (!scene || !path || !path[0])
	return 0;

    const int count = scene->getDatabaseSourceInstanceCountForPath(path);
    int appended = 0;
    for (int i = 0; i < count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceInstanceSummaryForPath(path, i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, draw_mode))
	    continue;
	const size_t previous_size = instance_keys.size();
	ged_obol_append_database_source_instance_key(instance_keys, summary);
	if (instance_keys.size() > previous_size)
	    appended++;
    }
    return appended;
}

static void
ged_obol_collect_database_sources_matching(
    std::vector<std::string> &paths,
    BObolSceneController *scene,
    struct ged_view_context *view_ctx,
    const char *target_path,
    size_t component_first_idx,
    int allow_path_prefix,
    int draw_mode)
{
    if (!scene || !target_path || !target_path[0])
	return;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (!source)
	    continue;

	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary,
		    draw_mode))
	    continue;

	const char *source_path = source->path.getValue().getString();
	if (ged_obol_path_equal(source_path, target_path) ||
	    ged_obol_semantic_path_string(source_path) ==
	    ged_obol_semantic_path_string(target_path) ||
	    (allow_path_prefix &&
	     (ged_obol_path_has_prefix(source_path, target_path) ||
	      ged_obol_path_has_semantic_prefix(source_path, target_path))) ||
	    ged_obol_path_has_component_name(source_path, target_path,
					     component_first_idx) ||
	    ged_obol_path_has_semantic_component_name(source_path, target_path,
		component_first_idx))
	    ged_obol_append_unique_path(paths, source_path);
    }
}

static void
ged_obol_collect_database_source_instance_keys_matching(
    std::vector<std::string> &instance_keys,
    BObolSceneController *scene,
    struct ged_view_context *view_ctx,
    const char *target_path,
    size_t component_first_idx,
    int allow_path_prefix,
    int draw_mode)
{
    if (!scene || !target_path || !target_path[0])
	return;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary,
		    draw_mode))
	    continue;

	const char *source_path = summary.path.getString();
	if (ged_obol_path_equal(source_path, target_path) ||
	    ged_obol_semantic_path_string(source_path) ==
	    ged_obol_semantic_path_string(target_path) ||
	    (allow_path_prefix &&
	     (ged_obol_path_has_prefix(source_path, target_path) ||
	      ged_obol_path_has_semantic_prefix(source_path, target_path))) ||
	    ged_obol_path_has_component_name(source_path, target_path,
					     component_first_idx) ||
	    ged_obol_path_has_semantic_component_name(source_path, target_path,
		component_first_idx))
	    ged_obol_append_database_source_instance_key(instance_keys,
		    summary);
    }
}

static std::vector<std::string>
ged_obol_matching_database_source_paths(
    BObolSceneController *scene,
    struct ged_view_context *view_ctx,
    const std::vector<std::string> &targets,
    size_t component_first_idx,
    int allow_path_prefix,
    int draw_mode)
{
    std::vector<std::string> paths;
    if (!scene)
	return paths;

    for (const std::string &target : targets)
	ged_obol_collect_database_sources_matching(paths, scene, view_ctx,
		target.c_str(), component_first_idx, allow_path_prefix,
		draw_mode);
    return paths;
}

static std::vector<std::string>
ged_obol_primary_matching_database_source_paths(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const std::vector<std::string> &targets,
    int draw_mode)
{
    std::vector<std::string> source_paths;
    if (!gedp || targets.empty())
	return source_paths;

    BObolSceneController *primary_scene =
	ged_draw_obol_scene_controller(gedp);
    if (!primary_scene)
	return source_paths;

    for (const std::string &target : targets) {
	ged_obol_collect_database_sources_matching(source_paths,
		primary_scene, view_ctx, target.c_str(), 0, 1, draw_mode);
    }
    ged_obol_remove_shadowed_source_paths(source_paths);
    return source_paths;
}

static std::vector<std::string>
ged_obol_matching_database_source_instance_keys(
    BObolSceneController *scene,
    struct ged_view_context *view_ctx,
    const std::vector<std::string> &targets,
    size_t component_first_idx,
    int allow_path_prefix,
    int draw_mode)
{
    std::vector<std::string> instance_keys;
    if (!scene)
	return instance_keys;

    for (const std::string &target : targets)
	ged_obol_collect_database_source_instance_keys_matching(
	    instance_keys, scene, view_ctx, target.c_str(),
	    component_first_idx, allow_path_prefix, draw_mode);
    return instance_keys;
}

static int
ged_obol_replace_matching_database_sources(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const std::vector<std::string> &targets,
    size_t component_first_idx,
    int allow_path_prefix,
    uint32_t source_revision,
    BObolSceneController *scene,
    int draw_mode)
{
    if (!gedp || !gedp->dbip || !scene || targets.empty())
	return 0;

    std::vector<std::string> source_paths =
	ged_obol_matching_database_source_paths(scene, view_ctx, targets,
	    component_first_idx, allow_path_prefix, draw_mode);
    ged_obol_remove_shadowed_source_paths(source_paths);
    if (source_paths.empty())
	return 0;

    std::vector<ged_obol_drawn_source_path_mode> source_path_modes =
	ged_obol_drawn_source_path_modes(gedp, view_ctx, draw_mode,
					 &source_paths);
    if (source_path_modes.empty()) {
	for (const std::string &source_path : source_paths) {
	    const int selected_mode = draw_mode >= 0 ? draw_mode :
				      ged_obol_drawn_path_mode(gedp, view_ctx,
					  source_path.c_str());
	    ged_obol_append_unique_path_mode(source_path_modes,
		    source_path.c_str(), selected_mode);
	}
    }

    const int changed = ged_obol_direct_draw_path_modes(gedp, view_ctx,
			source_path_modes, source_revision, scene, 1);
    (void)source_revision;
    return changed >= 0 ? 1 : 0;
}

static int
ged_obol_compact_source_matches_target(
    SoBRLDatabaseSource *source, const char *target)
{
    if (!source || !source->hasCompactInstanceIndex() ||
	!target || !target[0])
	return 0;
    int count = source->getCompactInstanceCountForPath(target, TRUE);
    if (count > 0)
	return count;

    const char *targetName = strrchr(target, '/');
    targetName = targetName && targetName[1] ? targetName + 1 : target;
    for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	BObolCompactInstanceHandle handle;
	BObolCompactInstanceSummary entry;
	if (source->getCompactInstanceHandle(i, handle) &&
	    source->getCompactInstanceSummary(handle, entry) &&
	    BU_STR_EQUAL(entry.sourceName.getString(), targetName))
	    count++;
    }
    return count;
}

static int
ged_obol_mark_matching_database_sources_stale(
    struct ged_view_context *view_ctx,
    const std::vector<std::string> &targets,
    size_t component_first_idx,
    int allow_path_prefix,
    uint32_t stale_reason,
    BObolSceneController *scene,
    int draw_mode)
{
    if (!scene || targets.empty())
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_matching_database_source_instance_keys(scene, view_ctx,
	    targets, component_first_idx, allow_path_prefix,
	    draw_mode);

    /* Nested compact occurrences are not independent source nodes.  Include
     * their owning source for geometry invalidation only; generic source
     * matching is also used by visibility and highlight transactions, where
     * promoting a nested path to the root would be incorrect. */
    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	BObolDatabaseSourceSummary summary;
	if (!source || !source->hasCompactInstanceIndex() ||
	    !source->getSummary(summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, draw_mode))
	    continue;
	for (const std::string &target : targets) {
	    if (ged_obol_compact_source_matches_target(source,
		    target.c_str()) <= 0)
		continue;
	    ged_obol_append_database_source_instance_key(instance_keys, summary);
	    break;
	}
    }
    if (instance_keys.empty())
	return 0;

    for (const std::string &instance_key : instance_keys) {
	bool refreshed = false;
	bool refreshFailed = false;
	if (stale_reason == SoBRLDatabaseSource::STALE_SOURCE) {
	    for (const std::string &target : targets) {
		const int ret = scene->refreshDatabaseSourceInstanceObject(
		    instance_key.c_str(), target.c_str(), 0);
		if (ret < 0) {
		    refreshFailed = true;
		    break;
		}
		if (ret > 0)
		    refreshed = true;
	    }
	}
	if (!refreshed || refreshFailed)
	    (void)scene->markDatabaseSourceInstanceStale(instance_key.c_str(),
		stale_reason);
    }
    return 1;
}

static int
ged_obol_refresh_matching_compact_parts(
    BObolSceneController *scene, struct ged_view_context *view_ctx,
    const std::vector<std::string> &targets, int draw_mode,
    uint32_t source_revision)
{
    if (!scene || targets.empty())
	return 0;

    int refreshed = 0;
    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	BObolDatabaseSourceSummary summary;
	if (!source || !source->isCompactOccurrenceRegistry() ||
	    !source->getSummary(summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, draw_mode))
	    continue;

	int sourceRefreshed = 0;
	const char *key = summary.instanceKey.getLength() > 0 ?
	    summary.instanceKey.getString() : summary.path.getString();
	for (const std::string &target : targets) {
	    if (ged_obol_compact_source_matches_target(source,
		    target.c_str()) <= 0)
		continue;
	    const int ret = scene->refreshDatabaseSourceInstanceObject(key,
		target.c_str(), source_revision);
	    if (ret > 0)
		sourceRefreshed += ret;
	}
	if (sourceRefreshed > 0)
	    refreshed += sourceRefreshed;
    }
    return refreshed;
}

static int
ged_obol_set_database_source_visible(BObolSceneController *scene,
				     const char *source_instance_key,
				     int visible)
{
    if (!scene || !source_instance_key || !source_instance_key[0])
	return 0;

    SoBRLDatabaseSource *source =
	scene->findDatabaseSourceInstance(source_instance_key);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    (void)scene->setDatabaseSourceInstanceState(source_instance_key,
	    TRUE,
	    summary.sourceRevision,
	    summary.inputsRevision,
	    visible ? TRUE : FALSE,
	    summary.selected,
	    summary.highlighted,
	    summary.lineStyle,
	    summary.lineWidth,
	    summary.transparency,
	    summary.colorOverride,
	    summary.color,
	    summary.materialColorValid,
	    summary.materialColor,
	    summary.materialRevision);
    return 1;
}

static int
ged_obol_set_group_visible(BObolSceneController *scene,
			   const char *group_path,
			   int visible)
{
    if (!scene || !group_path || !group_path[0])
	return 0;

    SoGroup *group = scene->findGroup(group_path);
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    (void)scene->setGroupDisplayState(group_path,
				      visible ? TRUE : FALSE,
				      scene_group->selected.getValue(),
				      scene_group->highlighted.getValue(),
				      scene_group->lineStyle.getValue(),
				      scene_group->lineWidth.getValue(),
				      scene_group->transparency.getValue(),
				      scene_group->colorOverride.getValue(),
				      scene_group->color.getValue(),
				      scene_group->materialColorValid.getValue(),
				      scene_group->materialColor.getValue(),
				      scene_group->materialRevision.getValue());
    return 1;
}

static int
ged_obol_set_database_source_highlighted(BObolSceneController *scene,
	const char *source_instance_key,
	int highlighted)
{
    if (!scene || !source_instance_key || !source_instance_key[0])
	return 0;

    SoBRLDatabaseSource *source =
	scene->findDatabaseSourceInstance(source_instance_key);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    const SbBool next_highlighted = highlighted ? TRUE : FALSE;
    int changed = scene->setDatabaseSourceInstanceState(source_instance_key,
	    TRUE,
	    summary.sourceRevision,
	    summary.inputsRevision,
	    summary.visible,
	    summary.selected,
	    next_highlighted,
	    summary.lineStyle,
	    summary.lineWidth,
	    summary.transparency,
	    summary.colorOverride,
	    summary.color,
	    summary.materialColorValid,
	    summary.materialColor,
	    summary.materialRevision);

    /* A compact occurrence registry deliberately preserves per-occurrence
     * selection styles when its source state changes.  This operation is the
     * explicit global highlight command, so it must update every occurrence. */
    if (source->hasCompactInstanceIndex()) {
	const int compact_changed = source->setCompactInstanceDisplayStateForPath(
		"", TRUE, 0, FALSE, 0, FALSE, 1, next_highlighted);
	if (compact_changed > 0)
	    changed = 1;
    }

    return changed > 0 ? 1 : 0;
}

static int
ged_obol_set_group_highlighted(BObolSceneController *scene,
			       const char *group_path,
			       int highlighted)
{
    if (!scene || !group_path || !group_path[0])
	return 0;

    SoGroup *group = scene->findGroup(group_path);
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    (void)scene->setGroupDisplayState(group_path,
				      scene_group->visible.getValue(),
				      scene_group->selected.getValue(),
				      highlighted ? TRUE : FALSE,
				      scene_group->lineStyle.getValue(),
				      scene_group->lineWidth.getValue(),
				      scene_group->transparency.getValue(),
				      scene_group->colorOverride.getValue(),
				      scene_group->color.getValue(),
				      scene_group->materialColorValid.getValue(),
				      scene_group->materialColor.getValue(),
				      scene_group->materialRevision.getValue());
    return 1;
}

static void
ged_obol_append_database_source_instance_key(
    std::vector<std::string> &instance_keys,
    const BObolDatabaseSourceSummary &summary);

extern "C" int
ged_draw_obol_highlight_state_set(struct ged *gedp, int highlighted)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid)
	    continue;
	ged_obol_append_database_source_instance_key(instance_keys, summary);
    }

    int changed = 0;
    for (const std::string &instance_key : instance_keys) {
	if (ged_obol_set_database_source_highlighted(scene,
		instance_key.c_str(), highlighted))
	    changed = 1;
    }

    const int tree_count = scene->getSceneTreeSummaryCount();
    for (int i = 0; i < tree_count; i++) {
	BObolSceneTreeSummary tree_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !tree_summary.valid ||
	    tree_summary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	    tree_summary.path.getLength() == 0 ||
	    BU_STR_EQUAL(tree_summary.path.getString(), "/"))
	    continue;
	if (ged_obol_set_group_highlighted(scene,
					   tree_summary.path.getString(), highlighted))
	    changed = 1;
    }

    if (changed)
	scene->realizePending();
    return changed;
}

static int
ged_obol_apply_highlight_transaction(
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    BObolSceneController *scene)
{
    if (!txn || !scene)
	return 0;

    std::vector<std::string> targets =
	ged_obol_transaction_paths(txn, result);
    if (targets.empty())
	return 0;

    const int highlighted = ZERO(txn->value) ? 0 : 1;
    int handled = 0;
    for (const std::string &target : targets) {
	for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	    SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	    if (source && source->setCompactInstanceDisplayStateForPath(
		    target.c_str(), TRUE, 0, FALSE, 0, FALSE, 1,
		    highlighted ? TRUE : FALSE) > 0)
		handled = 1;
	}
	std::vector<std::string> instance_keys;
	ged_obol_collect_database_source_instance_keys_matching(
	    instance_keys, scene, txn->view, target.c_str(), 0, 1,
	    txn->mode);
	for (const std::string &instance_key : instance_keys) {
	    if (ged_obol_set_database_source_highlighted(scene,
		    instance_key.c_str(), highlighted))
		handled = 1;
	}

	const int summary_count = scene->getSceneTreeSummaryCount();
	for (int i = 0; i < summary_count; i++) {
	    BObolSceneTreeSummary tree_summary;
	    if (!scene->getSceneTreeSummary(i, tree_summary) ||
		!tree_summary.valid ||
		tree_summary.nodeKind !=
		BObolSceneTreeSummary::NODE_GROUP ||
		tree_summary.path.getLength() == 0 ||
		BU_STR_EQUAL(tree_summary.path.getString(), "/"))
		continue;
	    const char *group_path = tree_summary.path.getString();
	    if (!ged_obol_path_has_prefix(group_path, target.c_str()) &&
		!ged_obol_path_has_component_name(group_path,
						  target.c_str(), 0))
		continue;
	    if (ged_obol_set_group_highlighted(scene, group_path,
					       highlighted))
		handled = 1;
	}
    }

    return handled;
}

static int
ged_obol_apply_visibility_transaction(
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    BObolSceneController *scene)
{
    if (!txn || !scene)
	return 0;

    std::vector<std::string> targets =
	ged_obol_transaction_paths(txn, result);
    if (targets.empty())
	return 0;

    const int visible = ZERO(txn->value) ? 0 : 1;
    int handled = 0;
    for (const std::string &target : targets) {
	for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	    SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	    if (source && source->setCompactInstanceDisplayStateForPath(
		    target.c_str(), TRUE, 1, visible ? TRUE : FALSE,
		    0, FALSE, 0, FALSE) > 0)
		handled = 1;
	}
	std::vector<std::string> instance_keys;
	ged_obol_collect_database_source_instance_keys_matching(
	    instance_keys, scene, txn->view, target.c_str(), 0, 1,
	    txn->mode);
	for (const std::string &instance_key : instance_keys) {
	    if (ged_obol_set_database_source_visible(scene,
		    instance_key.c_str(), visible))
		handled = 1;
	}

	const int summary_count = scene->getSceneTreeSummaryCount();
	for (int i = 0; i < summary_count; i++) {
	    BObolSceneTreeSummary tree_summary;
	    if (!scene->getSceneTreeSummary(i, tree_summary) ||
		!tree_summary.valid ||
		tree_summary.nodeKind !=
		BObolSceneTreeSummary::NODE_GROUP ||
		tree_summary.path.getLength() == 0 ||
		BU_STR_EQUAL(tree_summary.path.getString(), "/"))
		continue;
	    const char *group_path = tree_summary.path.getString();
	    if (!ged_obol_path_has_prefix(group_path, target.c_str()) &&
		!ged_obol_path_has_component_name(group_path,
						  target.c_str(), 0))
		continue;
	    if (ged_obol_set_group_visible(scene, group_path, visible))
		handled = 1;
	}
    }

    return handled;
}

static int
ged_obol_remove_paths(const std::vector<std::string> &paths,
		      struct ged_view_context *view_ctx,
		      BObolSceneController *scene,
		      int draw_mode = -1)
{
    if (paths.empty() || !scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary,
		    draw_mode))
	    continue;

	if (ged_obol_source_path_matches_any_target(summary.path.getString(),
		paths)) {
	    ged_obol_append_database_source_instance_key(instance_keys,
		    summary);
	}
    }

    if (instance_keys.empty()) {
	for (const std::string &path : paths) {
	    const std::string instance_key =
		ged_obol_database_source_instance_key_for_mode(view_ctx,
		    path.c_str(), draw_mode);
	    if (!instance_key.empty())
		ged_obol_append_unique_path(instance_keys,
					    instance_key.c_str());
	}
    }

    int changed = 0;
    for (const std::string &instance_key : instance_keys) {
	int removed = scene->removeDatabaseSourceInstance(instance_key.c_str());
	if (removed > 0)
	    changed = 1;
    }

    if (changed && ged_obol_prune_empty_groups(scene))
	changed = 1;

    return changed;
}

static void
ged_obol_append_database_source_instance_key(
    std::vector<std::string> &instance_keys,
    const BObolDatabaseSourceSummary &summary)
{
    if (!summary.valid)
	return;

    const char *instance_key = summary.instanceKey.getString();
    if (instance_key && instance_key[0]) {
	ged_obol_append_unique_path(instance_keys, instance_key);
	return;
    }

    const char *path = summary.path.getString();
    if (path && path[0])
	ged_obol_append_unique_path(instance_keys, path);
}

static int
ged_obol_append_existing_database_source_instance_key(
    std::vector<std::string> &instance_keys,
    BObolSceneController *scene,
    const char *source_instance_key,
    int draw_mode)
{
    if (!scene || !source_instance_key || !source_instance_key[0])
	return 0;

    SoBRLDatabaseSource *source =
	scene->findDatabaseSourceInstance(source_instance_key);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (source->getSummary(summary) && summary.valid) {
	if (!ged_obol_database_source_summary_matches_mode(summary,
		draw_mode))
	    return 0;
	ged_obol_append_database_source_instance_key(instance_keys, summary);
	return 1;
    }

    if (draw_mode >= 0)
	return 0;
    ged_obol_append_unique_path(instance_keys, source_instance_key);
    return 1;
}

static int
ged_obol_append_database_source_instance_key_for_source(
    std::vector<std::string> &instance_keys,
    SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (source->getSummary(summary) && summary.valid) {
	ged_obol_append_database_source_instance_key(instance_keys, summary);
	return 1;
    }

    const char *source_path = source->path.getValue().getString();
    if (!source_path || !source_path[0])
	return 0;

    ged_obol_append_unique_path(instance_keys, source_path);
    return 1;
}

static std::vector<std::string>
ged_obol_database_source_instance_keys_for_path(
    struct ged *gedp,
    const char *path,
    int draw_mode,
    int allow_path_prefix,
    const struct ged_bobol_publication_context *publication = NULL)
{
    std::vector<std::string> instance_keys;
    if (!gedp || !path || !path[0])
	return instance_keys;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return instance_keys;

    (void)ged_obol_append_existing_database_source_instance_key(
	instance_keys, scene, path, draw_mode);


    const bool explicit_publication = publication && publication->active;
    if (explicit_publication) {
	struct ged_view_context *publication_view = publication->view_ctx;
	const int publication_mode =
	    publication->draw_mode >= 0 ? publication->draw_mode : draw_mode;
	if (publication_mode >= 0) {
	    const std::string mode_key =
		ged_obol_database_source_instance_key_for_mode(
		    publication_view, path,
		    publication_mode);
	    (void)ged_obol_append_existing_database_source_instance_key(
		instance_keys, scene, mode_key.c_str(), publication_mode);
	}

	const std::string base_key =
	    ged_obol_database_source_instance_key(
		publication_view, path);
	(void)ged_obol_append_existing_database_source_instance_key(
	    instance_keys, scene, base_key.c_str(), publication_mode);
	if (!allow_path_prefix) {
	    const int exact_count =
		ged_obol_append_exact_database_source_instance_keys(instance_keys,
		    scene, publication_view, path,
		    publication_mode);
	    if (exact_count > 0)
		return instance_keys;
	}

	ged_obol_collect_database_source_instance_keys_matching(
	    instance_keys, scene, publication_view,
	    path, 0, allow_path_prefix, publication_mode);
	return instance_keys;
    }

    if (!allow_path_prefix) {
	const int exact_count =
	    ged_obol_append_exact_database_source_instance_keys(instance_keys,
		scene, NULL, path, draw_mode);
	if (exact_count > 0)
	    return instance_keys;
    }

    ged_obol_collect_database_source_instance_keys_matching(
	    instance_keys, scene, NULL, path, 0, allow_path_prefix,
	    draw_mode);

    SoBRLDatabaseSource *path_indexed_source =
	scene->findDatabaseSource(path);
    (void)ged_obol_append_database_source_instance_key_for_source(
	instance_keys, path_indexed_source);
    return instance_keys;
}

static std::vector<std::string>
ged_obol_database_source_instance_keys_for_record_or_path(
    struct ged *gedp,
    const char *path,
    const struct ged_draw_obol_database_source_record *record,
    int draw_mode,
    int allow_path_prefix)
{
    std::vector<std::string> instance_keys;
    if (!gedp)
	return instance_keys;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return instance_keys;

    if (record && record->instance_key && record->instance_key[0]) {
	(void)ged_obol_append_existing_database_source_instance_key(
	    instance_keys, scene, record->instance_key, draw_mode);
	if (!instance_keys.empty())
	    return instance_keys;
    }

    const char *record_path =
	record && record->database_path && record->database_path[0] ?
	record->database_path : NULL;
    const char *target_path = (path && path[0]) ? path : record_path;
    if (!target_path)
	return instance_keys;

    return ged_obol_database_source_instance_keys_for_path(gedp,
	    target_path, draw_mode, allow_path_prefix);
}

static SoBRLDatabaseSource *
ged_obol_database_source_for_instance_key(
    BObolSceneController *scene,
    const std::string &source_instance_key)
{
    if (!scene || source_instance_key.empty())
	return NULL;
    return scene->findDatabaseSourceInstance(source_instance_key.c_str());
}

static std::string
ged_obol_renamed_database_source_instance_key(
    const std::string &source_instance_key,
    const char *old_path,
    const char *new_path)
{
    if (source_instance_key.empty() || !old_path || !old_path[0] ||
	!new_path || !new_path[0])
	return source_instance_key;

    if (source_instance_key.compare(0, 14, "brlcad-direct:") == 0)
	return source_instance_key;

    std::string mode_suffix;
    std::string key_base(source_instance_key);
    const size_t marker_pos =
	key_base.rfind(ged_obol_database_source_mode_key_marker);
    if (marker_pos != std::string::npos) {
	mode_suffix = key_base.substr(marker_pos);
	key_base.erase(marker_pos);
    }

    const std::string old_norm = ged_obol_normalized_path_string(old_path);
    const std::string new_norm = ged_obol_normalized_path_string(new_path);
    if (old_norm.empty() || new_norm.empty())
	return source_instance_key;

    if (key_base == old_norm)
	return new_norm + mode_suffix;
    if (key_base == "/" + old_norm)
	return std::string("/") + new_norm + mode_suffix;

    if (key_base.size() > old_norm.size() &&
	key_base.compare(key_base.size() - old_norm.size(), old_norm.size(),
	    old_norm) == 0) {
	const size_t prefix_len = key_base.size() - old_norm.size();
	if (prefix_len > 0 &&
	    (key_base[prefix_len - 1] == ':' ||
	     key_base[prefix_len - 1] == '/'))
	    return key_base.substr(0, prefix_len) + new_norm + mode_suffix;
    }

    return source_instance_key;
}

static int
ged_obol_remove_instance_keys(const std::vector<std::string> &instance_keys,
			      BObolSceneController *scene)
{
    if (instance_keys.empty() || !scene)
	return 0;

    int changed = 0;
    for (const std::string &instance_key : instance_keys) {
	int removed = scene->removeDatabaseSourceInstance(
			  instance_key.c_str());
	if (removed > 0)
	    changed = 1;
    }

    if (changed && ged_obol_prune_empty_groups(scene))
	changed = 1;

    return changed;
}

static int
ged_obol_clear_database_sources_in_scope(BObolSceneController *scene,
	struct ged_view_context *view_ctx)
{
    if (!scene)
	return 0;

    if (!ged_obol_view_scope_is_independent(view_ctx)) {
	int changed = (scene->clearDatabaseSources() > 0) ? 1 : 0;
	if (changed && ged_obol_prune_empty_groups(scene))
	    changed = 1;
	return changed;
    }

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx))
	    continue;

	ged_obol_append_database_source_instance_key(instance_keys, summary);
    }

    return ged_obol_remove_instance_keys(instance_keys, scene);
}

BObolSceneController *
ged_draw_obol_scene_controller(struct ged *gedp)
{
    ged_obol_state *state = ged_obol_state_get(ged_obol_gdp(gedp), 0);
    return (state && state->shared_controller) ?
	state->shared_controller->getSceneController() : NULL;
}

BObolViewController *
ged_bobol_view_controller(struct ged_view_context *view_ctx)
{
    bobol_display_endpoint_t *endpoint = view_ctx ?
	ged_view_context_display_endpoint_get(view_ctx) : NULL;
    return endpoint ? static_cast<BObolViewController *>(
	bobol_display_endpoint_controller(endpoint)) : NULL;
}

BObolSceneController *
ged_bobol_scene(struct ged *gedp)
{
    return ged_draw_obol_scene_controller(gedp);
}

BObolViewController *
ged_bobol_shared_view_controller(struct ged *gedp)
{
    return ged_draw_obol_controller(gedp);
}

BObolViewController *
ged_draw_obol_controller(struct ged *gedp)
{
    ged_obol_state *state = ged_obol_state_get(ged_obol_gdp(gedp), 0);
    return state ? state->shared_controller : NULL;
}

static void
ged_obol_advance_lod_policy_revision(BObolViewController *controller)
{
    if (!controller)
	return;

    const uint64_t current_revision = controller->getLodPolicyRevision();
    controller->setLodPolicyRevision(
	current_revision == UINT64_MAX ? 1 : current_revision + 1);
}

extern "C" int
ged_draw_obol_view_lod_policy_changed(struct ged *gedp, struct ged_view_context *view_ctx)
{
    BObolViewController *view_controller =
	ged_bobol_view_controller(view_ctx);
    if (!view_controller)
	return 0;

    view_controller->clearViewLodState();
    ged_obol_advance_lod_policy_revision(view_controller);

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	scene = view_controller->getSceneController();

    int changed = 1;
    if (scene) {
	const ged_obol_view_lod_policy_state policy_state =
	    ged_obol_view_lod_policy_state_for_source(gedp, view_ctx);
	const int independent_view =
	    ged_obol_view_scope_is_independent(view_ctx);
	const int source_count = scene->getDatabaseSourceCount();
	for (int i = 0; i < source_count; i++) {
	    BObolDatabaseSourceSummary summary;
	    if (!scene->getDatabaseSourceSummary(i, summary) ||
		!summary.valid || summary.instanceKey.getLength() == 0)
		continue;
	    if (independent_view &&
		!ged_obol_database_source_instance_in_scope(summary, view_ctx))
		continue;
	    if (ged_obol_apply_view_lod_policy(gedp, view_ctx, scene,
					       summary.instanceKey.getString()) > 0)
		changed = 1;
	    /* A compact registry owns stable occurrence state.  LoD changes only
	     * its view-local backing payloads; marking it stale would rebuild the
	     * scene merely to turn LoD off. */
	    if (policy_state.valid && !policy_state.meshEnabled &&
		!scene->getDatabaseSource(i)->isCompactOccurrenceRegistry() &&
		scene->markDatabaseSourceInstanceStale(
		    summary.instanceKey.getString(),
		    SoBRLDatabaseSource::STALE_VIEW) > 0)
		changed = 1;
	}
	if (changed)
	    (void)scene->realizePending();
    }

    view_controller->requestRender("view-lod-policy");
    return changed;
}

static size_t
ged_obol_lod_default_worker_count(size_t worker_count)
{
    if (worker_count != 0)
	return worker_count;

    size_t cpus = bu_avail_cpus();
    /* A display service shares the machine with the UI, database snapshot
     * worker, renderer, and exact operations.  Consuming every logical CPU on
     * a large workstation made first-frame latency worse through contention. */
    return std::max(static_cast<size_t>(1),
	std::min(static_cast<size_t>(8), cpus));
}

static BObolLodService *
ged_obol_lod_service_ensure(BObolViewController *controller)
{
    if (!controller)
	return NULL;

    BObolLodService *service = controller->getLodService();
    if (service && service->isRunning()) {
	controller->setLodAutoSubmit(TRUE);
	return service;
    }

    service = controller->ensureManagedLodService(
	ged_obol_lod_default_worker_count(0));
    if (!service)
	return NULL;

    /* The view policy, not the call site that first needed a worker, owns
     * automatic display LoD.  Keeping this service "passive" made qged depend
     * on a host-specific startup command and left MGED/gsh on a different
     * path.  AUTO is the default; an OFF policy force-realizes and therefore
     * produces no display-LoD requests even though the worker pool exists. */
    controller->setLodAutoSubmit(TRUE);
    return service;
}

static void
ged_obol_lod_status_fill(BObolViewController *controller,
			 ged_draw_obol_lod_service_status_t *status)
{
    if (!status)
	return;

    memset(status, 0, sizeof(*status));
    if (!controller)
	return;

    BObolLodService *service = controller->getLodService();
    status->attached = 1;
    status->auto_submit = controller->isLodAutoSubmitEnabled() ? 1 : 0;
    status->worker_count = controller->getManagedLodWorkerCount();
    status->last_visited_mesh_count =
	controller->getLastLodVisitedMeshCount();
    status->last_submitted_task_count =
	controller->getLastLodSubmittedTaskCount();
    status->last_updated_cut_count =
	controller->getLastLodUpdatedCutCount();
    status->last_skipped_mesh_count =
	controller->getLastLodSkippedMeshCount();
    status->last_result_count =
	controller->getLastLodResultCount();
    status->last_matched_result_count =
	controller->getLastLodMatchedResultCount();
    status->last_applied_result_count =
	controller->getLastLodAppliedResultCount();
    status->last_rejected_result_count =
	controller->getLastLodRejectedResultCount();
    status->last_unmatched_result_count =
	controller->getLastLodUnmatchedResultCount();
    status->active_mesh_payloads =
	controller->getActiveLodMeshPayloadCount();
    status->active_aabb_proxy_payloads =
	controller->getActiveLodProxyPayloadCount(
	    BOBOL_LOD_PROXY_AABB);
    status->active_obb_proxy_payloads =
	controller->getActiveLodProxyPayloadCount(
	    BOBOL_LOD_PROXY_OBB);
    bu_strlcpy(status->last_diagnostics,
	       controller->getLastLodDiagnostics().getString(),
	       sizeof(status->last_diagnostics));
    if (!service)
	return;

    status->running = service->isRunning() ? 1 : 0;
    status->in_flight = service->inFlightCount();
    status->pending_tasks = service->pendingTaskCountForDiagnostics();
    status->queued_results = service->queuedResultCountForDiagnostics();
    status->queued_cache_writes =
	service->queuedCacheWriteCountForDiagnostics();
    status->delayed_tasks = service->delayedTaskCountForDiagnostics();
}

extern "C" int
ged_draw_obol_lod_service_start(struct ged *gedp,
				struct ged_view_context *view_ctx,
				size_t worker_count)
{
    (void)gedp;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	return 0;

    worker_count = ged_obol_lod_default_worker_count(worker_count);
    if (!controller->ensureManagedLodService(worker_count))
	return 0;

    controller->setLodAutoSubmit(TRUE);
    controller->requestRender("lod-service-start");
    return 1;
}

extern "C" int
ged_draw_obol_lod_service_stop(struct ged *gedp,
			       struct ged_view_context *view_ctx)
{
    (void)gedp;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	return 0;

    controller->stopManagedLodService();
    controller->requestRender("lod-service-stop");
    return 1;
}

extern "C" int
ged_draw_obol_lod_service_poll(struct ged *gedp,
			       struct ged_view_context *view_ctx,
			       size_t max_results,
			       ged_draw_obol_lod_service_status_t *status)
{
    (void)gedp;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	return 0;

    if (controller->hasPendingLodResults())
	(void)controller->processPendingLodResults(max_results);
    if (controller->isLodAutoSubmitEnabled())
	(void)controller->submitLodRequestsIfNeeded();
    ged_obol_lod_status_fill(controller, status);
    return 1;
}

extern "C" int
ged_draw_obol_lod_service_status(struct ged *gedp,
				 struct ged_view_context *view_ctx,
				 ged_draw_obol_lod_service_status_t *status)
{
    (void)gedp;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	return 0;

    ged_obol_lod_status_fill(controller, status);
    return 1;
}

static int
ged_obol_lod_prewarm_eligible(const struct directory *dp)
{
    if (!dp)
	return 0;
    if (dp->d_major_type != DB5_MAJORTYPE_BRLCAD)
	return 0;
    if (dp->d_addr == RT_DIR_PHONY_ADDR)
	return 0;
    if (dp->d_flags & RT_DIR_COMB)
	return 0;
    return dp->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT ? 1 : 0;
}

static size_t
ged_obol_lod_prewarm_submit_name(BObolLodService *service,
				 struct db_i *dbip,
				 const char *database_id,
				 uint64_t generation,
				 const char *name)
{
    if (!service || !service->isRunning() || !dbip || !name || !name[0])
	return 0;

    struct directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (!ged_obol_lod_prewarm_eligible(dp))
	return 0;

    BObolMeshLodProvider *provider = new BObolMeshLodProvider;
    provider->service = service;
    provider->dbip = dbip;
    provider->useView = FALSE;
    provider->refreshMissing = TRUE;
    provider->shrinkAfterCopy = FALSE;

    BObolLodTask task;
    task.generation = generation;
    task.request.databaseId = database_id ? database_id : "";
    task.request.objectPath = dp->d_namep;
    task.request.objectName = dp->d_namep;
    task.request.providerId = "bobol_mesh_lod_cache";
    task.request.providerVersion = "bobol-cache-v1";
    task.request.qualityTier = BOBOL_LOD_QUALITY_FAST_DISPLAY;
    task.realize = bobol_mesh_lod_cache_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = bobol_mesh_lod_provider_free;
    task.publishResult = FALSE;

    uint64_t taskId = service->submitIfNotActive(task);
    if (taskId == 0) {
	bobol_mesh_lod_provider_free(provider);
	return 0;
    }

    return 1;
}

extern "C" size_t
ged_draw_obol_lod_service_prewarm(struct ged *gedp,
				  struct ged_view_context *view_ctx,
				  int argc,
				  const char * const *argv,
				  ged_draw_obol_lod_service_status_t *status)
{
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	return 0;

    BObolLodService *service = controller->getLodService();
    if (!service || !service->isRunning()) {
	service = ged_obol_lod_service_ensure(controller);
	if (!service || !service->isRunning())
	    return 0;
    }

    struct db_i *dbip = gedp ? gedp->dbip : NULL;
    if (!dbip)
	return 0;

    const char *database_id = dbip->dbi_filename ? dbip->dbi_filename : "";
    const uint64_t previous_generation = service->currentGeneration();
    if (previous_generation != 0)
	service->cancelGeneration(previous_generation);
    uint64_t generation = service->beginGeneration();
    size_t submitted = 0;
    int submit_all = argc <= 0 || !argv;
    if (!submit_all) {
	for (int i = 0; i < argc; i++) {
	    if (argv[i] && (BU_STR_EQUAL(argv[i], "all") ||
			    BU_STR_EQUAL(argv[i], "*"))) {
		submit_all = 1;
		break;
	    }
	}
    }

    if (submit_all) {
	struct directory **paths = NULL;
	size_t path_cnt = db_ls(dbip, DB_LS_HIDDEN, NULL, &paths);
	for (size_t i = 0; i < path_cnt; i++) {
	    struct directory *dp = paths[i];
	    if (!ged_obol_lod_prewarm_eligible(dp))
		continue;
	    submitted += ged_obol_lod_prewarm_submit_name(
			     service, dbip, database_id, generation,
			     dp->d_namep);
	}
	if (paths)
	    bu_free(paths, "free obol lod prewarm db_ls output");
    } else {
	for (int i = 0; i < argc; i++)
	    submitted += ged_obol_lod_prewarm_submit_name(
			     service, dbip, database_id, generation, argv[i]);
    }

    ged_obol_lod_status_fill(controller, status);
    return submitted;
}

extern "C" int
ged_draw_obol_scene_controller_owned_internal(struct ged *gedp)
{
    ged_obol_state *state = ged_obol_state_get(ged_obol_gdp(gedp), 0);
    return (state && state->shared_controller) ? 1 : 0;
}

int
ged_draw_obol_scene_controller_owned(struct ged *gedp)
{
    return ged_draw_obol_scene_controller_owned_internal(gedp);
}

extern "C" int
ged_draw_obol_scene_controller_attached(struct ged *gedp)
{
    return ged_draw_obol_scene_controller(gedp) ? 1 : 0;
}

extern "C" int
ged_draw_obol_scene_controller_ensure_owned(struct ged *gedp,
	int sync_current_scene)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp)
	return 0;
    if (!ged_obol_observer_ensure(gedp, gdp))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (scene) {
	ged_obol_state *state = ged_obol_state_get(gdp, 0);
	if (sync_current_scene && state && !state->full_sync) {
	    (void)ged_draw_obol_scene_sync_full_scene(gedp, NULL, 0, scene);
	    state->full_sync = 1;
	}
	return 1;
    }

    bobol_init(NULL);

    SoBRLSceneGroup *root = new SoBRLSceneGroup;
    BObolViewController *owned_controller = new BObolViewController(root);
    BObolSceneController *owned_scene = owned_controller->getSceneController();
    ged_obol_state *state = ged_obol_state_get(gdp, 1);
    if (!state) {
	delete owned_controller;
	return 0;
    }
    state->shared_controller = owned_controller;
    state->full_sync = 0;

    if (sync_current_scene) {
	(void)ged_draw_obol_scene_sync_full_scene(gedp, NULL, 0,
		owned_scene);
	state->full_sync = 1;
    }

    return 1;
}

extern "C" int
ged_draw_obol_scene_controller_full_synced(struct ged *gedp)
{
    ged_obol_state *state = ged_obol_state_get(ged_obol_gdp(gedp), 0);
    return (state && state->shared_controller && state->full_sync) ? 1 : 0;
}

static unsigned char
ged_obol_rgb_channel(float value)
{
    if (value <= 0.0f)
	return 0;
    if (value >= 1.0f)
	return 255;
    return static_cast<unsigned char>(value * 255.0f + 0.5f);
}

static void
ged_obol_rgb_from_color(const SbColor &color, unsigned char rgb[3])
{
    rgb[0] = ged_obol_rgb_channel(color.getValue()[0]);
    rgb[1] = ged_obol_rgb_channel(color.getValue()[1]);
    rgb[2] = ged_obol_rgb_channel(color.getValue()[2]);
}

static void
ged_obol_mat_from_sbmatrix(const SbMatrix &matrix, mat_t mat)
{
    const SbMat &m = matrix.getValue();
    mat[0] = m[0][0];
    mat[4] = m[0][1];
    mat[8] = m[0][2];
    mat[12] = m[0][3];
    mat[1] = m[1][0];
    mat[5] = m[1][1];
    mat[9] = m[1][2];
    mat[13] = m[1][3];
    mat[2] = m[2][0];
    mat[6] = m[2][1];
    mat[10] = m[2][2];
    mat[14] = m[2][3];
    mat[3] = m[3][0];
    mat[7] = m[3][1];
    mat[11] = m[3][2];
    mat[15] = m[3][3];
}

static SbMatrix
ged_obol_sbmatrix_from_mat(const mat_t mat)
{
    return SbMatrix(
	       static_cast<float>(mat[0]), static_cast<float>(mat[4]),
	       static_cast<float>(mat[8]), static_cast<float>(mat[12]),
	       static_cast<float>(mat[1]), static_cast<float>(mat[5]),
	       static_cast<float>(mat[9]), static_cast<float>(mat[13]),
	       static_cast<float>(mat[2]), static_cast<float>(mat[6]),
	       static_cast<float>(mat[10]), static_cast<float>(mat[14]),
	       static_cast<float>(mat[3]), static_cast<float>(mat[7]),
	       static_cast<float>(mat[11]), static_cast<float>(mat[15]));
}

static SbColor
ged_obol_color_from_rgb(const unsigned char rgb[3])
{
    return SbColor(
	       static_cast<float>(rgb[0]) / 255.0f,
	       static_cast<float>(rgb[1]) / 255.0f,
	       static_cast<float>(rgb[2]) / 255.0f);
}

/* ged_draw_appearance_settings predates the retained scene and calls its
 * opacity field "transparency".  Keep that command-facing contract at this
 * boundary; Obol display state uses conventional transparency. */
static float
ged_obol_transparency_from_appearance_opacity(fastf_t opacity)
{
    if (opacity <= 0.0)
	return 1.0f;
    if (opacity >= 1.0)
	return 0.0f;
    return 1.0f - static_cast<float>(opacity);
}

static fastf_t
ged_obol_appearance_opacity_from_transparency(float transparency)
{
    if (transparency <= 0.0f)
	return 1.0;
    if (transparency >= 1.0f)
	return 0.0;
    return static_cast<fastf_t>(1.0f - transparency);
}

static fastf_t
ged_obol_reported_transparency(float transparency)
{
    const double scaled = (double)transparency * 1000000.0;
    return (fastf_t)(floor(scaled + 0.5) / 1000000.0);
}

static BObolViewController *
ged_obol_view_controller_for_context(struct ged_view_context *view_ctx)
{
    return ged_bobol_view_controller(view_ctx);
}

static BObolViewController *
ged_obol_view_controller_ensure_for_context(struct ged_view_context *view_ctx,
	int sync_current_scene)
{
    struct ged *gedp = view_ctx ?
			   static_cast<struct ged *>(ged_view_context_user_data_get(view_ctx)) :
		       NULL;
    if (!gedp)
	return NULL;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (controller)
	return controller;
    if (!ged_draw_obol_render_endpoint_ensure_for_view(gedp, view_ctx,
	    sync_current_scene))
	return NULL;
    return ged_bobol_view_controller(view_ctx);
}

static BObolViewController *
ged_obol_shared_view_controller_ensure_for_context(struct ged_view_context *view_ctx,
	int sync_current_scene)
{
    struct ged *gedp = view_ctx ?
			   static_cast<struct ged *>(ged_view_context_user_data_get(view_ctx)) :
		       NULL;
    if (!gedp)
	return NULL;
    BObolViewController *controller = ged_draw_obol_controller(gedp);
    if (controller)
	return controller;
    if (!ged_draw_obol_scene_controller_ensure_owned(gedp, sync_current_scene))
	return NULL;
    return ged_draw_obol_controller(gedp);
}

static BObolViewController *
ged_obol_shared_view_controller_for_context(struct ged_view_context *view_ctx)
{
    struct ged *gedp = view_ctx ?
			   static_cast<struct ged *>(ged_view_context_user_data_get(view_ctx)) :
		       NULL;
    return gedp ? ged_draw_obol_controller(gedp) : NULL;
}

BObolViewController *
ged_obol_view_controller_for_scope(struct ged_view_context *view_ctx,
				   int local,
				   int sync_current_scene)
{
    return local ?
	   ged_obol_view_controller_ensure_for_context(view_ctx,
		   sync_current_scene) :
	   ged_obol_shared_view_controller_ensure_for_context(view_ctx,
		   sync_current_scene);
}

struct ged_obol_pick_candidate {
    ged_obol_pick_candidate(void) :
	distance(FLT_MAX),
	source_id(0),
	primitive_kind(0),
	primitive_index(-1),
	material_id(0),
	nearest_face_vertex_index(-1),
	model_point(0.0f, 0.0f, 0.0f)
    {
	face_vertex_index[0] = -1;
	face_vertex_index[1] = -1;
	face_vertex_index[2] = -1;
    }

    std::string path;
    float distance;
    uint32_t source_id;
    int primitive_kind;
    int primitive_index;
    int material_id;
    int face_vertex_index[3];
    int nearest_face_vertex_index;
    SbVec3f model_point;
};

static std::string
ged_obol_normalized_pick_path(const std::string &path)
{
    size_t start = 0;
    while (start < path.size() && path[start] == '/')
	start++;
    return path.substr(start);
}

static std::string
ged_obol_pick_candidate_key(const ged_obol_pick_candidate &candidate)
{
    char buffer[128] = {0};
    snprintf(buffer, sizeof(buffer), ":%u:%d:%d:%d",
	    candidate.source_id, candidate.primitive_kind,
	    candidate.primitive_index, candidate.material_id);
    return ged_obol_normalized_pick_path(candidate.path) + buffer;
}

static BObolFeatureOwner
ged_obol_pick_view_owner(struct ged_view_context *view_ctx)
{
    BObolFeatureOwner owner;
    owner.ownerToken = view_ctx;
    owner.ownerRole = "view";

    const struct bv *view = ged_obol_bv_const(view_ctx);
    const char *name = bv_name_get(view);
    if (name && name[0]) {
	owner.ownerId = name;
	return owner;
    }

    char fallback[64] = {0};
    snprintf(fallback, sizeof(fallback), "%p", (void *)view_ctx);
    owner.ownerId = fallback;
    return owner;
}

static std::string
ged_obol_pick_resolved_path(BObolViewController *controller,
			    const BObolFeatureOwner *owner,
			    const SoBRLPickDetail *detail)
{
    if (!detail)
	return std::string();

    const char *path = detail->getPath().getString();
    const char *source_name = detail->getSourceName().getString();
    const std::string picked_name =
	(path && path[0]) ? std::string(path) :
	(source_name ? std::string(source_name) : std::string());
    if (!controller || picked_name.empty() || detail->getPrimitiveIndex() < 0)
	return picked_name;

    BObolFeaturePrimitivePick pick;
    if ((owner && controller->features().resolvePrimitivePick(
	    picked_name.c_str(), detail->getPrimitiveIndex(), pick,
	    BOBOL_FEATURE_SCOPE_LOCAL, owner)) ||
	    controller->features().resolvePrimitivePick(picked_name.c_str(),
		detail->getPrimitiveIndex(), pick,
		BOBOL_FEATURE_SCOPE_SHARED, NULL))
	return std::string(pick.featureName.getString());

    return picked_name;
}

static ged_obol_pick_candidate
ged_obol_pick_candidate_from_detail(BObolViewController *controller,
				    const BObolFeatureOwner *owner,
				    const SoBRLPickDetail *detail,
				    const SbVec3f &point,
				    float distance)
{
    ged_obol_pick_candidate candidate;
    candidate.path = ged_obol_pick_resolved_path(controller, owner, detail);
    candidate.distance = distance;
    candidate.source_id = detail ? detail->getSourceId() : 0;
    candidate.primitive_kind = detail ?
	static_cast<int>(detail->getPrimitiveKind()) : 0;
    candidate.primitive_index = detail ? detail->getPrimitiveIndex() : -1;
    candidate.material_id = detail ? detail->getMaterialId() : 0;
    if (detail) {
	candidate.face_vertex_index[0] = detail->getFaceVertexIndexA();
	candidate.face_vertex_index[1] = detail->getFaceVertexIndexB();
	candidate.face_vertex_index[2] = detail->getFaceVertexIndexC();
	candidate.nearest_face_vertex_index =
	    detail->getNearestFaceVertexIndex();
	candidate.model_point = detail->getModelPoint();
    } else {
	candidate.model_point = point;
    }
    return candidate;
}

static ged_obol_pick_candidate
ged_obol_pick_candidate_from_record(BObolViewController *controller,
				    const BObolFeatureOwner *owner,
				    const BObolViewPickRecord &record)
{
    return ged_obol_pick_candidate_from_detail(controller, owner,
	&record.detail, record.point, record.distance);
}

static void
ged_obol_pick_insert(std::vector<ged_obol_pick_candidate> &candidates,
		     const ged_obol_pick_candidate &candidate,
		     bool pick_all)
{
    if (candidate.path.empty())
	return;

    if (pick_all) {
	candidates.push_back(candidate);
	return;
    }

    if (candidates.empty() || candidate.distance < candidates[0].distance)
	candidates.assign(1, candidate);
}

static bool
ged_obol_pick_candidate_nearer(const ged_obol_pick_candidate &a,
			       const ged_obol_pick_candidate &b)
{
    return a.distance < b.distance;
}

static void
ged_obol_pick_sort(std::vector<ged_obol_pick_candidate> &candidates)
{
    if (candidates.size() > 1)
	std::stable_sort(candidates.begin(), candidates.end(),
		ged_obol_pick_candidate_nearer);
}

static struct ged_pick_result *
ged_obol_pick_result_from_candidates(
	const std::vector<ged_obol_pick_candidate> &candidates)
{
    struct ged_pick_result *result = ged_pick_result_create();
    if (!result)
	return NULL;

    for (size_t i = 0; i < candidates.size(); i++) {
	struct ged_pick_detail detail = GED_PICK_DETAIL_INIT;
	detail.source_id = candidates[i].source_id;
	detail.primitive_kind = candidates[i].primitive_kind;
	detail.primitive_index = candidates[i].primitive_index;
	detail.material_id = candidates[i].material_id;
	detail.face_vertex_index[0] = candidates[i].face_vertex_index[0];
	detail.face_vertex_index[1] = candidates[i].face_vertex_index[1];
	detail.face_vertex_index[2] = candidates[i].face_vertex_index[2];
	detail.nearest_face_vertex_index =
	    candidates[i].nearest_face_vertex_index;
	detail.model_point[X] = candidates[i].model_point[0];
	detail.model_point[Y] = candidates[i].model_point[1];
	detail.model_point[Z] = candidates[i].model_point[2];
	detail.model_point_valid = 1;
	(void)ged_pick_result_append_detail(result,
		candidates[i].path.c_str(), candidates[i].distance, &detail);
    }

    return result;
}

static int
ged_obol_pick_sync_view_controller(BObolViewController *controller,
				   struct ged_view_context *view_ctx)
{
    if (!controller || !view_ctx)
	return 0;

    const struct bv *view = ged_obol_bv_const(view_ctx);
    const int width = bv_width_get(view);
    const int height = bv_height_get(view);
    if (width > 0 && height > 0)
	controller->setViewportSize(static_cast<unsigned int>(width),
		static_cast<unsigned int>(height));

    controller->syncCameraFromViewContext(view_ctx, TRUE);
    return 1;
}

static int
ged_obol_pick_point_candidates(
	struct ged_view_context *view_ctx,
	int x,
	int y,
	float radius_pixels,
	bool pick_all,
	std::vector<ged_obol_pick_candidate> &candidates)
{
    candidates.clear();
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;
    ged_obol_pick_sync_view_controller(controller, view_ctx);

    std::vector<BObolViewPickRecord> records;
    (void)bobol_view_pick_point(controller, x, y, radius_pixels,
	pick_all, records, NULL);
    BObolFeatureOwner owner = ged_obol_pick_view_owner(view_ctx);
    candidates.reserve(records.size());
    for (const BObolViewPickRecord &record : records)
	ged_obol_pick_insert(candidates,
	    ged_obol_pick_candidate_from_record(controller, &owner, record),
	    pick_all);

    if (pick_all)
	ged_obol_pick_sort(candidates);

    return static_cast<int>(candidates.size());
}

extern "C" struct ged_pick_result *
ged_draw_obol_view_context_pick_point(struct ged_view_context *view_ctx,
				      int x,
				      int y,
				      int first_only)
{
    std::vector<ged_obol_pick_candidate> candidates;
    (void)ged_obol_pick_point_candidates(view_ctx, x, y, 1.0f,
	    first_only ? false : true, candidates);
    return ged_obol_pick_result_from_candidates(candidates);
}

extern "C" struct ged_pick_result *
ged_draw_obol_view_context_pick_nearest(struct ged_view_context *view_ctx, int x, int y)
{
    return ged_draw_obol_view_context_pick_point(view_ctx, x, y, 1);
}

extern "C" struct ged_pick_result *
ged_draw_obol_view_context_pick_rect(struct ged_view_context *view_ctx,
				     int x0,
				     int y0,
				     int x1,
				     int y1)
{
    std::vector<ged_obol_pick_candidate> candidates;
    std::unordered_set<std::string> seen;
    int min_x = std::min(x0, x1);
    int max_x = std::max(x0, x1);
    int min_y = std::min(y0, y1);
    int max_y = std::max(y0, y1);

    const struct bv *view = ged_obol_bv_const(view_ctx);
    const int width_px = bv_width_get(view);
    const int height_px = bv_height_get(view);
    if (width_px > 0) {
	min_x = std::max(0, std::min(min_x, width_px - 1));
	max_x = std::max(0, std::min(max_x, width_px - 1));
    }
    if (height_px > 0) {
	min_y = std::max(0, std::min(min_y, height_px - 1));
	max_y = std::max(0, std::min(max_y, height_px - 1));
    }

    int width = std::max(1, max_x - min_x);
    int height = std::max(1, max_y - min_y);
    int x_steps = std::max(1, std::min(6, width / 16));
    int y_steps = std::max(1, std::min(6, height / 16));

    for (int yi = 0; yi <= y_steps; yi++) {
	int y = min_y + (height * yi) / y_steps;
	for (int xi = 0; xi <= x_steps; xi++) {
	    int x = min_x + (width * xi) / x_steps;
	    std::vector<ged_obol_pick_candidate> sampled;
	    ged_obol_pick_point_candidates(view_ctx, x, y, 1.0f, true,
		    sampled);
	    for (size_t i = 0; i < sampled.size(); i++) {
		const std::string key = ged_obol_pick_candidate_key(sampled[i]);
		if (!seen.insert(key).second)
		    continue;
		candidates.push_back(sampled[i]);
	    }
	}
    }

    ged_obol_pick_sort(candidates);
    return ged_obol_pick_result_from_candidates(candidates);
}

static uint32_t
ged_obol_snap_kind(enum ged_selection_snap_kind kind)
{
    switch (kind) {
	case GED_SELECTION_SNAP_GRID:
	    return static_cast<uint32_t>(SoBRLSnapAction::GRID);
	case GED_SELECTION_SNAP_ENDPOINT:
	    return static_cast<uint32_t>(SoBRLSnapAction::ENDPOINT);
	default:
	    return 0;
    }
}

extern "C" int
ged_draw_obol_view_context_snap_first_candidate(
	struct ged_view_context *view_ctx,
	const point_t sample,
	enum ged_selection_snap_kind kind,
	point_t candidate)
{
    if (candidate)
	VSETALL(candidate, 0.0);
    if (!view_ctx || !sample || !candidate)
	return 0;

    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;
    if (!ged_obol_pick_sync_view_controller(controller, view_ctx))
	return 0;

    const uint32_t enabled_kinds = ged_obol_snap_kind(kind);
    if (!enabled_kinds)
	return 0;

    BObolViewSnapRecord snap;
    if (!bobol_view_snap_point(controller,
	SbVec3f(static_cast<float>(sample[X]), static_cast<float>(sample[Y]),
	    static_cast<float>(sample[Z])), FLT_MAX, enabled_kinds,
	SoBRLSnapAction::FULL_DETAIL, false, snap))
	return 0;

    const SbVec3f &point = snap.point;
    VSET(candidate, point[0], point[1], point[2]);
    return 1;
}

BObolFeatureStyle
ged_obol_feature_style_from_ged(
    const struct ged_view_feature_style *style)
{
    BObolFeatureStyle out;
    if (!style)
	return out;

    if (style->visible != -1) {
	out.hasVisible = TRUE;
	out.visible = style->visible ? TRUE : FALSE;
    }
    if (style->selectable != -1) {
	out.hasSelectable = TRUE;
	out.selectable = style->selectable ? TRUE : FALSE;
    }
    if (style->color_valid) {
	out.hasColor = TRUE;
	out.color = ged_obol_color_from_rgb(style->color);
    }
    if (style->line_width >= 0) {
	out.hasLineWidth = TRUE;
	out.lineWidth = style->line_width;
    }
    if (style->line_style >= 0) {
	out.hasLineStyle = TRUE;
	out.lineStyle = style->line_style;
    }
    if (style->arrow != -1) {
	out.hasArrow = TRUE;
	out.arrow = style->arrow ? TRUE : FALSE;
    }
    if (style->arrow_tip_length >= 0.0 || style->arrow_tip_width >= 0.0) {
	out.hasArrowTip = TRUE;
	out.arrowTipLength = style->arrow_tip_length >= 0.0 ?
			     static_cast<float>(style->arrow_tip_length) : 0.0f;
	out.arrowTipWidth = style->arrow_tip_width >= 0.0 ?
			    static_cast<float>(style->arrow_tip_width) : 0.0f;
    }
    return out;
}

int32_t
ged_obol_line_command_from_ged(int command, size_t index)
{
    if (command == GED_DRAW_VIEW_LINE_DRAW ||
	command == BG_GEOMETRY_LINE_DRAW)
	return static_cast<int32_t>(BObolLineCommand::Draw);
    if (command == GED_DRAW_VIEW_LINE_POINT_DRAW ||
	command == BG_GEOMETRY_POINT_DRAW)
	return static_cast<int32_t>(BObolLineCommand::Point);
    if (command == GED_DRAW_VIEW_LINE_MOVE ||
	command == BG_GEOMETRY_LINE_MOVE)
	return static_cast<int32_t>(BObolLineCommand::Move);
    return static_cast<int32_t>(index ? BObolLineCommand::Draw :
				BObolLineCommand::Move);
}

static int
ged_obol_line_command_to_ged(int32_t command)
{
    if (command == static_cast<int32_t>(BObolLineCommand::Draw))
	return GED_DRAW_VIEW_LINE_DRAW;
    if (command == static_cast<int32_t>(BObolLineCommand::Point))
	return GED_DRAW_VIEW_LINE_POINT_DRAW;
    return GED_DRAW_VIEW_LINE_MOVE;
}

std::vector<SbVec3f>
ged_obol_points_from_ged(const point_t *points, size_t point_count)
{
    std::vector<SbVec3f> out;
    if (!points || !point_count)
	return out;

    out.reserve(point_count);
    for (size_t i = 0; i < point_count; i++)
	out.push_back(SbVec3f(
			  static_cast<float>(points[i][X]),
			  static_cast<float>(points[i][Y]),
			  static_cast<float>(points[i][Z])));
    return out;
}

std::vector<int32_t>
ged_obol_commands_from_ged(const int *cmds, size_t point_count)
{
    std::vector<int32_t> out;
    out.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	const int command = cmds ? cmds[i] : -1;
	out.push_back(ged_obol_line_command_from_ged(command, i));
    }
    return out;
}

std::vector<int32_t>
ged_obol_indices_from_ged(const int *indices, size_t index_count)
{
    std::vector<int32_t> out;
    if (!indices || !index_count)
	return out;

    out.reserve(index_count);
    for (size_t i = 0; i < index_count; i++)
	out.push_back(static_cast<int32_t>(indices[i]));
    return out;
}

std::vector<SbVec3f>
ged_obol_vectors_from_ged(const vect_t *vectors, size_t vector_count)
{
    std::vector<SbVec3f> out;
    if (!vectors || !vector_count)
	return out;

    out.reserve(vector_count);
    for (size_t i = 0; i < vector_count; i++)
	out.push_back(SbVec3f(
			  static_cast<float>(vectors[i][X]),
			  static_cast<float>(vectors[i][Y]),
			  static_cast<float>(vectors[i][Z])));
    return out;
}

static BObolFeatureHandle
ged_obol_feature_handle(BObolViewController *controller,
			struct ged_view_context *view_ctx,
			const char *name)
{
    if (!controller || !name)
	return BObolFeatureHandle();

    BObolFeatureOwner owner;
    owner.ownerToken = view_ctx;
    owner.ownerId = ged_obol_view_scope_name(view_ctx).c_str();
    owner.ownerRole = "view";

    BObolFeatureHandle local =
	controller->features().findOwned(name, BOBOL_FEATURE_SCOPE_LOCAL,
					 &owner);
    if (local.isValid())
	return local;

    BObolFeatureHandle shared =
	controller->features().find(name, BOBOL_FEATURE_SCOPE_SHARED);
    if (shared.isValid())
	return shared;

    return controller->features().find(name);
}

struct ged_obol_feature_lookup {
    BObolViewController *controller;
    BObolFeatureHandle handle;
};

static ged_obol_feature_lookup
ged_obol_feature_lookup_for_context(struct ged_view_context *view_ctx, const char *name)
{
    ged_obol_feature_lookup out;
    out.controller = NULL;
    out.handle = BObolFeatureHandle();
    if (!name)
	return out;

    BObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (local_controller) {
	BObolFeatureHandle handle =
	    ged_obol_feature_handle(local_controller, view_ctx, name);
	if (handle.isValid()) {
	    out.controller = local_controller;
	    out.handle = handle;
	    return out;
	}
    }

    BObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if (shared_controller && shared_controller != local_controller) {
	BObolFeatureHandle handle =
	    shared_controller->features().find(name,
					       BOBOL_FEATURE_SCOPE_SHARED);
	if (handle.isValid()) {
	    out.controller = shared_controller;
	    out.handle = handle;
	}
    }

    return out;
}

static BObolFeatureScope
ged_obol_feature_scope(int local)
{
    return local ? BObolFeatureScope::Local : BObolFeatureScope::Shared;
}

static BObolFeatureOwner
ged_obol_feature_owner(struct ged_view_context *view_ctx, int local)
{
    BObolFeatureOwner owner;
    if (!local)
	return owner;

    owner.ownerToken = view_ctx;
    owner.ownerId = ged_obol_view_scope_name(view_ctx).c_str();
    owner.ownerRole = "view";
    return owner;
}

BObolOverlayInfo
ged_obol_model_overlay_info(struct ged_view_context *view_ctx,
			    BObolOverlayClass overlay_class,
			    BObolOverlayLifecycle lifecycle,
			    BObolOverlayOrder order,
			    const char *source_path)
{
    BObolOverlayInfo info;
    info.isOverlay = TRUE;
    info.ownerToken = view_ctx;
    info.role = BObolOverlayRole::Model;
    info.overlayClass = overlay_class;
    info.lifecycle = lifecycle;
    info.order = order;
    info.sortOrder = 0;
    info.sourcePath = source_path ? source_path : "";
    return info;
}

int
ged_obol_feature_mark_overlay(BObolViewController *controller,
			      BObolFeatureHandle handle,
			      const BObolOverlayInfo &overlay)
{
    if (!controller || !handle.isValid())
	return 0;
    return controller->features().setOverlayInfo(handle, overlay) ? 1 : 0;
}


static BObolFeatureHandle
ged_obol_publish_line_set(BObolViewController *controller,
			  struct ged_view_context *view_ctx,
			  const char *name,
			  int local,
			  const std::vector<SbVec3f> &points,
			  const std::vector<int32_t> &commands,
			  const BObolFeatureStyle *style)
{
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishLineSet(name,
	    ged_obol_feature_scope(local), points, commands, style,
	    local ? &owner : NULL);
}

static int
ged_obol_remove_feature(BObolViewController *controller,
			struct ged_view_context *view_ctx,
			const char *name,
			int local_mode)
{
    if (!controller || !name)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle;
    if (local_mode > 0) {
	handle = controller->features().findOwned(name,
		 BOBOL_FEATURE_SCOPE_LOCAL, &owner);
    } else if (local_mode == 0) {
	handle = controller->features().find(name,
					     BOBOL_FEATURE_SCOPE_SHARED);
    } else {
	handle = ged_obol_feature_handle(controller, view_ctx, name);
    }

    return handle.isValid() ? (controller->features().remove(handle) ? 1 : 0) : 0;
}

static BObolFeatureHandle
ged_obol_publish_indexed_face_set(BObolViewController *controller,
				  struct ged_view_context *view_ctx,
				  const char *name,
				  int local,
				  const std::vector<SbVec3f> &points,
				  const std::vector<SbVec3f> &normals,
				  const std::vector<int32_t> &indices,
				  const BObolFeatureStyle *style)
{
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishIndexedFaceSet(name,
	    ged_obol_feature_scope(local), points, normals, indices, style,
	    local ? &owner : NULL);
}

static BObolFeatureHandle
ged_obol_publish_labels(BObolViewController *controller,
			struct ged_view_context *view_ctx,
			const char *name,
			int local,
			const std::vector<BObolLabel> &labels,
			const BObolFeatureStyle *style = NULL)
{
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishLabels(name,
	    ged_obol_feature_scope(local), labels, style,
	    local ? &owner : NULL);
}

static BObolFeatureHandle
ged_obol_publish_arrow(BObolViewController *controller,
		       struct ged_view_context *view_ctx,
		       const char *name,
		       int local,
		       const std::vector<SbVec3f> &points,
		       const BObolFeatureStyle *style)
{
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishArrow(name,
	    ged_obol_feature_scope(local), points, style,
	    local ? &owner : NULL);
}

static BObolFeatureHandle
ged_obol_publish_axes(BObolViewController *controller,
		      struct ged_view_context *view_ctx,
		      const char *name,
		      int local,
		      const std::vector<SbVec3f> &centers,
		      float half_axes_size,
		      const BObolFeatureStyle *style)
{
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishAxes(name,
	    ged_obol_feature_scope(local), centers, half_axes_size, style,
	    local ? &owner : NULL);
}

static void
ged_obol_point_from_sb(point_t dst, const SbVec3f &src)
{
    VSET(dst, src[0], src[1], src[2]);
}

static BObolLabel
ged_obol_label_from_ged(const struct ged_annotation_label &label)
{
    BObolLabel out;
    out.text = label.text ? label.text : "";
    out.point = SbVec3f(
		    static_cast<float>(label.point[X]),
		    static_cast<float>(label.point[Y]),
		    static_cast<float>(label.point[Z]));
    if (label.color_valid) {
	out.hasColor = TRUE;
	out.color = ged_obol_color_from_rgb(label.color);
    }
    if (label.line_flag) {
	out.hasLeader = TRUE;
	out.target = SbVec3f(
			 static_cast<float>(label.target[X]),
			 static_cast<float>(label.target[Y]),
			 static_cast<float>(label.target[Z]));
    }
    out.anchor = label.anchor;
    out.arrow = label.arrow ? TRUE : FALSE;
    if (label.font_size > 0.0)
	out.fontSize = static_cast<float>(label.font_size);
    return out;
}

BObolLabel
ged_obol_label_from_hud(const struct ged_diagnostic_hud_label &label)
{
    BObolLabel out;
    out.text = label.text ? label.text : "";
    out.point = SbVec3f(
		    static_cast<float>(label.position[0]),
		    static_cast<float>(label.position[1]),
		    0.0f);
    unsigned char rgb[3] = {
	static_cast<unsigned char>(label.color[0] < 0 ? 0 :
	(label.color[0] > 255 ? 255 : label.color[0])),
	static_cast<unsigned char>(label.color[1] < 0 ? 0 :
	(label.color[1] > 255 ? 255 : label.color[1])),
	static_cast<unsigned char>(label.color[2] < 0 ? 0 :
	(label.color[2] > 255 ? 255 : label.color[2]))
    };
    out.hasColor = TRUE;
    out.color = ged_obol_color_from_rgb(rgb);
    out.fontSize = label.font_size > 0.0 ?
		   static_cast<float>(label.font_size) : 12.0f;
    out.sourceId = label.source_id;
    return out;
}

static std::vector<BObolLabel>
ged_obol_labels_from_ged(const struct ged_annotation_label *labels,
			 size_t label_count)
{
    std::vector<BObolLabel> out;
    if (!labels || !label_count)
	return out;

    out.reserve(label_count);
    for (size_t i = 0; i < label_count; i++)
	out.push_back(ged_obol_label_from_ged(labels[i]));
    return out;
}

static int
ged_obol_rgb_is_zero(const int rgb[3])
{
    return !rgb || (rgb[0] == 0 && rgb[1] == 0 && rgb[2] == 0);
}

static unsigned char
ged_obol_clamp_color_int(int v)
{
    return static_cast<unsigned char>(v < 0 ? 0 : (v > 255 ? 255 : v));
}

static SbColor
ged_obol_color_from_int_rgb(const int rgb[3],
			    int fallback_r,
			    int fallback_g,
			    int fallback_b)
{
    int r = fallback_r;
    int g = fallback_g;
    int b = fallback_b;
    if (!ged_obol_rgb_is_zero(rgb)) {
	r = rgb[0];
	g = rgb[1];
	b = rgb[2];
    }
    const unsigned char crgb[3] = {
	ged_obol_clamp_color_int(r),
	ged_obol_clamp_color_int(g),
	ged_obol_clamp_color_int(b)
    };
    return ged_obol_color_from_rgb(crgb);
}

static BObolFeatureStyle
ged_obol_faceplate_style(const int rgb[3],
			 int fallback_r,
			 int fallback_g,
			 int fallback_b,
			 int line_width)
{
    BObolFeatureStyle style;
    style.hasVisible = TRUE;
    style.visible = TRUE;
    style.hasSelectable = TRUE;
    style.selectable = FALSE;
    style.hasColor = TRUE;
    style.color = ged_obol_color_from_int_rgb(rgb, fallback_r, fallback_g,
		  fallback_b);
    style.hasLineWidth = TRUE;
    style.lineWidth = line_width > 0 ? line_width : 1;
    return style;
}

static BObolOverlayInfo
ged_obol_faceplate_overlay_info(struct ged_view_context *view_ctx,
				BObolOverlayOrder order =
				    BObolOverlayOrder::Screen)
{
    BObolOverlayInfo info;
    info.isOverlay = TRUE;
    info.ownerToken = view_ctx;
    info.role = BObolOverlayRole::Screen;
    info.overlayClass = BObolOverlayClass::Faceplate;
    info.lifecycle = BObolOverlayLifecycle::PerView;
    info.order = order;
    info.sortOrder = 0;
    info.sourcePath = "_faceplate";
    return info;
}

static int
ged_obol_view_to_model_point(SbVec3f &out,
			     struct ged_view_context *view_ctx,
			     fastf_t x,
			     fastf_t y,
			     fastf_t z = 0.0)
{
    mat_t view2model;
    if (!bv_view2model_get(view2model, ged_obol_bv_const(view_ctx)))
	return 0;

    point_t vpt;
    point_t mpt;
    VSET(vpt, x, y, z);
    MAT4X3PNT(mpt, view2model, vpt);
    out = SbVec3f(static_cast<float>(mpt[X]),
		  static_cast<float>(mpt[Y]),
		  static_cast<float>(mpt[Z]));
    return 1;
}

/* Map a faceplate overlay coordinate (GED2PM1 space, -1..1 in both axes) to the
 * pixel coordinates an SoHUDKit expects (origin bottom-left, 0..width/height),
 * matching how HUD text labels are placed.  This keeps the yellow faceplate
 * lines screen-locked -- like legacy main's fixed glOrtho(-1,1,-1,1) overlay --
 * instead of routing them through the view->model transform (which made them
 * move and skew with the camera/aspect). */
static void
ged_obol_faceplate_to_pixel(SbVec3f &out,
			    struct ged_view_context *view_ctx,
			    fastf_t x,
			    fastf_t y)
{
    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    out = SbVec3f(static_cast<float>((x + 1.0) * 0.5 * width),
		  static_cast<float>((y + 1.0) * 0.5 * height),
		  0.0f);
}

static void
ged_obol_faceplate_append_line(std::vector<SbVec3f> &points,
			       std::vector<int32_t> &commands,
			       struct ged_view_context *view_ctx,
			       fastf_t x1,
			       fastf_t y1,
			       fastf_t x2,
			       fastf_t y2)
{
    SbVec3f a;
    SbVec3f b;
    if (!ged_obol_view_to_model_point(a, view_ctx, x1, y1) ||
	!ged_obol_view_to_model_point(b, view_ctx, x2, y2))
	return;

    points.push_back(a);
    commands.push_back(static_cast<int32_t>(BObolLineCommand::Move));
    points.push_back(b);
    commands.push_back(static_cast<int32_t>(BObolLineCommand::Draw));
}

static BObolFeatureHandle
ged_obol_faceplate_publish_lines(BObolViewController *controller,
				 struct ged_view_context *view_ctx,
				 const char *name,
				 const std::vector<SbVec3f> &points,
				 const std::vector<int32_t> &commands,
				 const BObolFeatureStyle &style)
{
    if (!controller || !name || points.empty() || commands.empty())
	return BObolFeatureHandle();

    BObolFeatureHandle handle = ged_obol_publish_line_set(controller,
				  view_ctx, name, 1, points, commands, &style);
    (void)ged_obol_feature_mark_overlay(controller, handle,
					ged_obol_faceplate_overlay_info(view_ctx));
    return handle;
}

static BObolFeatureHandle
ged_obol_faceplate_publish_hud_labels(BObolViewController *controller,
				      struct ged_view_context *view_ctx,
				      const char *name,
				      const std::vector<BObolLabel> &labels,
				      const BObolFeatureStyle &style)
{
    if (!controller || !name || labels.empty())
	return BObolFeatureHandle();

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle =
	controller->features().publishHudLabels(name,
	    BObolFeatureScope::Local, labels, &style, &owner);
    (void)ged_obol_feature_mark_overlay(controller, handle,
					ged_obol_faceplate_overlay_info(view_ctx));
    return handle;
}

static void
ged_obol_faceplate_remove(BObolViewController *controller,
			  struct ged_view_context *view_ctx,
			  const char *name)
{
    (void)ged_obol_remove_feature(controller, view_ctx, name, 1);
}

static BObolLabel
ged_obol_faceplate_label(struct ged_view_context *view_ctx,
			 const char *text,
			 fastf_t x,
			 fastf_t y,
			 const int rgb[3],
			 int fallback_r,
			 int fallback_g,
			 int fallback_b,
			 int font_size,
			 int anchor = 0)
{
    BObolLabel label;
    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    fastf_t px = x;
    fastf_t py = y;
    if (width > 0 && height > 0) {
	px = (x + 1.0) * 0.5 * (fastf_t)width;
	py = (y + 1.0) * 0.5 * (fastf_t)height;
    }
    label.text = text ? text : "";
    label.point = SbVec3f(static_cast<float>(px),
			  static_cast<float>(py),
			  0.0f);
    label.hasColor = TRUE;
    label.color = ged_obol_color_from_int_rgb(rgb, fallback_r, fallback_g,
		  fallback_b);
    label.fontSize = static_cast<float>(font_size > 0 ? font_size : 20);
    label.anchor = anchor;
    return label;
}

static void
ged_obol_faceplate_sync_center_dot(BObolViewController *controller,
				   struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/center_dot";
    struct bv_other_state state = BV_OTHER_STATE_INIT;
    if (!bv_center_dot_state_get(&state, ged_obol_bv_const(view_ctx)) ||
	!state.gos_draw) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(4);
    commands.reserve(4);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   -0.01, 0.0, 0.01, 0.0);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   0.0, -0.01, 0.0, 0.01);
    BObolFeatureStyle style = ged_obol_faceplate_style(
				    state.gos_line_color, 255, 255, 0, 1);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx, name,
					   points, commands, style);
}

static void
ged_obol_faceplate_sync_interactive_rect(BObolViewController *controller,
					 struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/interactive_rect";
    struct bv_interactive_rect_state state = BV_INTERACTIVE_RECT_STATE_INIT;
    if (!bv_interactive_rect_state_get(&state, ged_obol_bv_const(view_ctx)) ||
	!state.draw || (ZERO(state.width) && ZERO(state.height))) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    const fastf_t aspect = width > 0 && height > 0 ?
	(fastf_t)width / (fastf_t)height : 1.0;
    const fastf_t x0 = state.x;
    const fastf_t x1 = state.x + state.width;
    const fastf_t y0 = state.y * aspect;
    const fastf_t y1 = (state.y + state.height) * aspect;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(8);
    commands.reserve(8);
    ged_obol_faceplate_append_line(points, commands, view_ctx, x0, y0, x0, y1);
    ged_obol_faceplate_append_line(points, commands, view_ctx, x0, y1, x1, y1);
    ged_obol_faceplate_append_line(points, commands, view_ctx, x1, y1, x1, y0);
    ged_obol_faceplate_append_line(points, commands, view_ctx, x1, y0, x0, y0);

    BObolFeatureStyle style = ged_obol_faceplate_style(state.color,
	255, 255, 255, state.line_width > 0 ? state.line_width : 1);
    style.lineStyle = state.line_style;
    (void)ged_obol_faceplate_publish_lines(controller, view_ctx, name,
					   points, commands, style);
}

static void
ged_obol_faceplate_sync_grid(BObolViewController *controller,
			     struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/grid";
    struct bv_grid_state grid = BV_GRID_STATE_INIT;
    if (!bv_grid_state_get(&grid, ged_obol_bv_const(view_ctx)) || !grid.draw) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    SoBRLGrid *node = new SoBRLGrid;
    node->ref();
    node->overlayId = name;
    if (!bobol_grid_configure_from_view_context(node, &grid, view_ctx) ||
	node->getTotalSegmentCount() <= 0) {
	node->unref();
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    BObolFeatureStyle style;
    style.hasVisible = TRUE;
    style.visible = TRUE;
    style.hasSelectable = TRUE;
    style.selectable = FALSE;
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle =
	controller->features().publishCustomNode(name,
	    BObolFeatureScope::Local, node, &style, &owner);
    node->unref();
    (void)ged_obol_feature_mark_overlay(controller, handle,
					ged_obol_faceplate_overlay_info(view_ctx));
}

static void
ged_obol_faceplate_sync_adc(BObolViewController *controller,
			    struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/adc";
    struct bv_adc_state state = BV_ADC_STATE_INIT;
    if (!bv_adc_state_get(&state, ged_obol_bv_const(view_ctx)) ||
	!state.draw) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    SoBRLADC *node = new SoBRLADC;
    node->ref();
    node->overlayId = name;
    node->center = SbVec3f(static_cast<float>(state.pos_model[X]),
			  static_cast<float>(state.pos_model[Y]),
			  static_cast<float>(state.pos_model[Z]));
    node->angleDegrees = static_cast<float>(state.a1);
    node->distance = static_cast<float>(state.dst > SMALL_FASTF ?
					state.dst : 1.0);
    node->lineColor = SbColor(state.line_color[0] / 255.0f,
	state.line_color[1] / 255.0f, state.line_color[2] / 255.0f);
    node->tickColor = SbColor(state.tick_color[0] / 255.0f,
	state.tick_color[1] / 255.0f, state.tick_color[2] / 255.0f);
    node->lineWidth = state.line_width > 0 ? state.line_width : 1;
    node->visible = TRUE;
    node->rebuildGeometry();

    BObolFeatureStyle style = ged_obol_faceplate_style(
	state.line_color, 255, 255, 255, state.line_width);
    style.hasSelectable = TRUE;
    style.selectable = FALSE;
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle =
	controller->features().publishCustomNode(name,
	    BObolFeatureScope::Local, node, &style, &owner);
    node->unref();
    (void)ged_obol_feature_mark_overlay(controller, handle,
	ged_obol_faceplate_overlay_info(view_ctx));
}

static std::vector<std::string>
ged_obol_faceplate_params_parts(struct ged_view_context *view_ctx,
				const struct bv_params_state *params)
{
    std::vector<std::string> parts;
    if (!view_ctx || !params)
	return parts;

    point_t center = VINIT_ZERO;
    mat_t center_mat;
    if (bv_center_mat_get(center_mat, ged_obol_bv_const(view_ctx)))
	MAT_DELTAS_GET_NEG(center, center_mat);
    VSCALE(center, center, bv_base2local_get(ged_obol_bv_const(view_ctx)));

    const char *ustr = bu_units_string(bv_local2base_get(ged_obol_bv_const(view_ctx)));
    if (!ustr)
	ustr = "";

    vect_t aet = VINIT_ZERO;
    (void)bv_aet_get(aet, ged_obol_bv_const(view_ctx));

    struct bu_vls text = BU_VLS_INIT_ZERO;
    if (params->draw_size) {
	bu_vls_sprintf(&text, "size[%s]: %.2f", ustr,
		      bv_size_get(ged_obol_bv_const(view_ctx)) *
		      bv_base2local_get(ged_obol_bv_const(view_ctx)));
	parts.push_back(bu_vls_cstr(&text));
    }
    if (params->draw_center) {
	bu_vls_sprintf(&text, "center[%s]: (%.2f, %.2f, %.2f)",
		      ustr, V3ARGS(center));
	parts.push_back(bu_vls_cstr(&text));
    }
    if (params->draw_az) {
	bu_vls_sprintf(&text, "az:%.2f", aet[0]);
	parts.push_back(bu_vls_cstr(&text));
    }
    if (params->draw_el) {
	bu_vls_sprintf(&text, "el:%.2f", aet[1]);
	parts.push_back(bu_vls_cstr(&text));
    }
    if (params->draw_tw) {
	bu_vls_sprintf(&text, "tw:%.2f", aet[2]);
	parts.push_back(bu_vls_cstr(&text));
    }

    const uint64_t frametime = bv_frametime_get(ged_obol_bv_const(view_ctx));
    if (params->draw_fps && frametime > 0) {
	bu_vls_sprintf(&text, "FPS:%.1f",
		1000000000.0 / (fastf_t)frametime);
	parts.push_back(bu_vls_cstr(&text));
    }
    bu_vls_free(&text);
    return parts;
}

static void
ged_obol_faceplate_sync_params(BObolViewController *controller,
			       struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/params";
    struct bv_params_state params = BV_PARAMS_STATE_INIT;
    if (!bv_params_state_get(&params, ged_obol_bv_const(view_ctx)) ||
	!params.draw) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    if (!params.draw_size && !params.draw_center && !params.draw_az &&
	!params.draw_el && !params.draw_tw && !params.draw_fps) {
	params.draw_size = 1;
	params.draw_center = 1;
	params.draw_az = 1;
	params.draw_el = 1;
	params.draw_tw = 1;
    }

    const int font_size = params.font_size > 0 ? params.font_size : 20;
    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    const int available_width = width > 16 ? width - 16 : width;
    const size_t max_chars = available_width > 0 ?
	std::max((size_t)1, (size_t)((fastf_t)available_width /
		(0.58 * (fastf_t)font_size))) : (size_t)80;

    const std::vector<std::string> parts =
	ged_obol_faceplate_params_parts(view_ctx, &params);
    std::vector<std::string> lines;
    std::string line;
    for (size_t i = 0; i < parts.size(); i++) {
	const size_t joined_len = line.size() + (line.empty() ? 0 : 1) +
	    parts[i].size();
	if (!line.empty() && joined_len > max_chars) {
	    lines.push_back(line);
	    line.clear();
	}
	if (!line.empty())
	    line.append(" ");
	line.append(parts[i]);
    }
    if (!line.empty())
	lines.push_back(line);

    std::vector<BObolLabel> labels;
    const fastf_t line_step = height > 0 ?
	(2.0 * (fastf_t)(font_size + 4) / (fastf_t)height) : 0.10;
    for (size_t i = 0; i < lines.size(); i++) {
	int line_font_size = font_size;
	if (available_width > 0 && !lines[i].empty()) {
	    const int fit_size = (int)((fastf_t)available_width /
		(0.58 * (fastf_t)lines[i].size()));
	    line_font_size = std::max(6, std::min(font_size, fit_size));
	}
	labels.push_back(ged_obol_faceplate_label(view_ctx, lines[i].c_str(),
		-0.98, -0.90 + (fastf_t)i * line_step,
		params.color, 255, 255, 0, line_font_size));
    }

    if (labels.empty()) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    BObolFeatureStyle style = ged_obol_faceplate_style(
				    params.color, 255, 255, 0, 1);
    (void)ged_obol_faceplate_publish_hud_labels(controller, view_ctx, name,
	    labels, style);
}

static void
ged_obol_faceplate_sync_scale(BObolViewController *controller,
			      struct ged_view_context *view_ctx)
{
    static const char line_name[] = "_faceplate/scale";
    static const char label_name[] = "_faceplate/scale_labels";
    struct bv_other_state state = BV_OTHER_STATE_INIT;
    if (!bv_scale_overlay_state_get(&state, ged_obol_bv_const(view_ctx)) ||
	!state.gos_draw) {
	ged_obol_faceplate_remove(controller, view_ctx, line_name);
	ged_obol_faceplate_remove(controller, view_ctx, label_name);
	return;
    }

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(6);
    commands.reserve(6);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   -0.5, -0.8, 0.5, -0.8);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   -0.5, -0.79, -0.5, -0.81);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   0.5, -0.79, 0.5, -0.81);

    BObolFeatureStyle line_style = ged_obol_faceplate_style(
					 state.gos_line_color, 255, 255, 0, 1);
    (void)ged_obol_faceplate_publish_lines(controller, view_ctx, line_name,
					   points, commands, line_style);

    struct bu_vls scale = BU_VLS_INIT_ZERO;
    const fastf_t base2local = bv_base2local_get(ged_obol_bv_const(view_ctx));
    const char *unit = !ZERO(base2local) ? bu_units_string(1.0 / base2local) :
		       NULL;
    if (!unit)
	unit = "";
    bu_vls_printf(&scale, "%g%s",
		  bv_size_get(ged_obol_bv_const(view_ctx)) * 0.5 *
		  base2local,
		  unit);
    const int soffset = (int)(strlen(bu_vls_cstr(&scale)) * 0.5);
    std::vector<BObolLabel> labels;
    labels.push_back(ged_obol_faceplate_label(view_ctx, "0", -0.505, -0.78,
		     state.gos_line_color, 255, 255, 0,
		     state.gos_font_size > 0 ? state.gos_font_size : 20));
    labels.push_back(ged_obol_faceplate_label(view_ctx, bu_vls_cstr(&scale),
		     0.5 - (soffset * 0.015), -0.78,
		     state.gos_line_color, 255, 255, 0,
		     state.gos_font_size > 0 ? state.gos_font_size : 20));
    bu_vls_free(&scale);

    BObolFeatureStyle text_style = ged_obol_faceplate_style(
					 state.gos_line_color, 255, 255, 0, 1);
    (void)ged_obol_faceplate_publish_hud_labels(controller, view_ctx,
	    label_name, labels, text_style);
}

static void
ged_obol_faceplate_append_axis(std::vector<SbVec3f> &points,
			       std::vector<int32_t> &commands,
			       std::vector<BObolLabel> &labels,
			       struct ged_view_context *view_ctx,
			       const mat_t rmat,
			       const struct bv_axes_state *axes,
			       int axis,
			       fastf_t aspect,
			       const char *label_text)
{
    point_t v2 = VINIT_ZERO;
    v2[axis] = axes->axes_size > 0.0 ? axes->axes_size * 0.5 : 0.1;

    point_t rv2;
    point_t rv1;
    MAT4X3PNT(rv2, rmat, v2);
    if (axes->pos_only) {
	VSET(rv1, 0.0, 0.0, 0.0);
    } else {
	VSCALE(rv1, rv2, -1.0);
    }

    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   rv1[X] + axes->axes_pos[X],
				   (rv1[Y] + axes->axes_pos[Y]) * aspect,
				   rv2[X] + axes->axes_pos[X],
				   (rv2[Y] + axes->axes_pos[Y]) * aspect);

    if (axes->label_flag) {
	point_t lv;
	point_t lrv;
	VSET(lv, v2[X] + 0.0078125, v2[Y] + 0.0078125,
	     v2[Z] + 0.0078125);
	MAT4X3PNT(lrv, rmat, lv);
	labels.push_back(ged_obol_faceplate_label(view_ctx, label_text,
			 lrv[X] + axes->axes_pos[X],
			 lrv[Y] + axes->axes_pos[Y],
			 axes->label_color, 255, 255, 0, 20));
    }
}

static void
ged_obol_faceplate_axis_triple_color(int axis, int rgb[3])
{
    if (!rgb)
	return;

    switch (axis) {
	case X:
	    VSET(rgb, 255, 0, 0);
	    break;
	case Y:
	    VSET(rgb, 0, 255, 0);
	    break;
	case Z:
	    VSET(rgb, 0, 0, 255);
	    break;
	default:
	    VSET(rgb, 255, 255, 255);
	    break;
    }
}

static const char *
ged_obol_faceplate_axis_suffix(int axis)
{
    switch (axis) {
	case X:
	    return "/x";
	case Y:
	    return "/y";
	case Z:
	    return "/z";
	default:
	    return "/axis";
    }
}

static void
ged_obol_faceplate_remove_axis_variants(BObolViewController *controller,
					struct ged_view_context *view_ctx,
					const std::string &line_name)
{
    ged_obol_faceplate_remove(controller, view_ctx, line_name.c_str());
    for (int axis = X; axis <= Z; axis++) {
	std::string axis_name = line_name + ged_obol_faceplate_axis_suffix(axis);
	ged_obol_faceplate_remove(controller, view_ctx, axis_name.c_str());
    }
}

static void
ged_obol_faceplate_append_tick_segment(std::vector<SbVec3f> &points,
				       std::vector<int32_t> &commands,
				       struct ged_view_context *view_ctx,
				       const fastf_t axes_pos[3],
				       const point_t t1,
				       const point_t t2,
				       fastf_t aspect)
{
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   t1[X] + axes_pos[X], (t1[Y] + axes_pos[Y]) * aspect,
				   t2[X] + axes_pos[X], (t2[Y] + axes_pos[Y]) * aspect);
}

static void
ged_obol_faceplate_append_axis_ticks(std::vector<SbVec3f> &tick_points,
				     std::vector<int32_t> &tick_commands,
				     std::vector<SbVec3f> &major_points,
				     std::vector<int32_t> &major_commands,
				     struct ged_view_context *view_ctx,
				     const mat_t rmat,
				     const struct bv_axes_state *axes,
				     fastf_t aspect)
{
    if (!view_ctx || !axes || !axes->tick_enabled ||
	axes->tick_interval <= SMALL_FASTF)
	return;

    const fastf_t view_size = bv_size_get(ged_obol_bv_const(view_ctx));
    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const fastf_t half_axes_size =
	axes->axes_size > 0.0 ? axes->axes_size * 0.5 : 0.1;
    if (view_size <= SMALL_FASTF || width <= 0 ||
	half_axes_size <= SMALL_FASTF)
	return;

    int num_ticks = static_cast<int>(
			view_size / axes->tick_interval * 0.5 * half_axes_size);
    if (num_ticks <= 0)
	return;

    int do_major_only = 0;
    if (axes->tick_threshold > 0 &&
	width <= num_ticks / half_axes_size * axes->tick_threshold * 2) {
	const int ticks_per_major = axes->ticks_per_major > 0 ?
				    axes->ticks_per_major : 1;
	const int num_major_ticks = num_ticks / ticks_per_major;
	if (width <= num_major_ticks / half_axes_size *
	    axes->tick_threshold * 2)
	    return;
	do_major_only = 1;
    }

    const fastf_t interval = axes->tick_interval / view_size * 2.0;
    const fastf_t tlen = axes->tick_length / (fastf_t)width * 2.0;
    const fastf_t maj_tlen =
	axes->tick_major_length / (fastf_t)width * 2.0;

    vect_t xend1 = VINIT_ZERO;
    vect_t xend2 = VINIT_ZERO;
    vect_t yend1 = VINIT_ZERO;
    vect_t yend2 = VINIT_ZERO;
    vect_t zend1 = VINIT_ZERO;
    vect_t zend2 = VINIT_ZERO;
    vect_t maj_xend1 = VINIT_ZERO;
    vect_t maj_xend2 = VINIT_ZERO;
    vect_t maj_yend1 = VINIT_ZERO;
    vect_t maj_yend2 = VINIT_ZERO;
    vect_t maj_zend1 = VINIT_ZERO;
    vect_t maj_zend2 = VINIT_ZERO;
    vect_t rxdir = VINIT_ZERO;
    vect_t neg_rxdir = VINIT_ZERO;
    vect_t rydir = VINIT_ZERO;
    vect_t neg_rydir = VINIT_ZERO;
    vect_t rzdir = VINIT_ZERO;
    vect_t neg_rzdir = VINIT_ZERO;
    vect_t dir = VINIT_ZERO;

    if (!do_major_only) {
	VSET(dir, tlen, 0.0, 0.0);
	MAT4X3PNT(xend1, rmat, dir);
	VSCALE(xend2, xend1, -1.0);
	VSET(dir, 0.0, tlen, 0.0);
	MAT4X3PNT(yend1, rmat, dir);
	VSCALE(yend2, yend1, -1.0);
	VSET(dir, 0.0, 0.0, tlen);
	MAT4X3PNT(zend1, rmat, dir);
	VSCALE(zend2, zend1, -1.0);
    }

    VSET(dir, maj_tlen, 0.0, 0.0);
    MAT4X3PNT(maj_xend1, rmat, dir);
    VSCALE(maj_xend2, maj_xend1, -1.0);
    VSET(dir, 0.0, maj_tlen, 0.0);
    MAT4X3PNT(maj_yend1, rmat, dir);
    VSCALE(maj_yend2, maj_yend1, -1.0);
    VSET(dir, 0.0, 0.0, maj_tlen);
    MAT4X3PNT(maj_zend1, rmat, dir);
    VSCALE(maj_zend2, maj_zend1, -1.0);

    VSET(dir, interval, 0.0, 0.0);
    MAT4X3PNT(rxdir, rmat, dir);
    VSCALE(neg_rxdir, rxdir, -1.0);
    VSET(dir, 0.0, interval, 0.0);
    MAT4X3PNT(rydir, rmat, dir);
    VSCALE(neg_rydir, rydir, -1.0);
    VSET(dir, 0.0, 0.0, interval);
    MAT4X3PNT(rzdir, rmat, dir);
    VSCALE(neg_rzdir, rzdir, -1.0);

    auto append_tick_pair = [&](const vect_t e1, const vect_t e2,
    const vect_t along, int major) {
	point_t t1;
	point_t t2;
	VADD2(t1, e1, along);
	VADD2(t2, e2, along);
	if (major)
	    ged_obol_faceplate_append_tick_segment(major_points,
						   major_commands, view_ctx, axes->axes_pos, t1, t2,
						   aspect);
	else
	    ged_obol_faceplate_append_tick_segment(tick_points,
						   tick_commands, view_ctx, axes->axes_pos, t1, t2,
						   aspect);
    };

    for (int i = 1; i <= num_ticks; i++) {
	const int major = axes->ticks_per_major > 0 ?
			  (i % axes->ticks_per_major == 0) : 0;
	if (!major && do_major_only)
	    continue;

	const vect_t *x1 = major ? &maj_xend1 : &xend1;
	const vect_t *x2 = major ? &maj_xend2 : &xend2;
	const vect_t *y1 = major ? &maj_yend1 : &yend1;
	const vect_t *y2 = major ? &maj_yend2 : &yend2;
	const vect_t *z1 = major ? &maj_zend1 : &zend1;
	const vect_t *z2 = major ? &maj_zend2 : &zend2;

	vect_t tvec;
	VSCALE(tvec, rxdir, i);
	append_tick_pair(*y1, *y2, tvec, major);
	append_tick_pair(*z1, *z2, tvec, major);
	if (!axes->pos_only) {
	    VSCALE(tvec, neg_rxdir, i);
	    append_tick_pair(*y1, *y2, tvec, major);
	    append_tick_pair(*z1, *z2, tvec, major);
	}

	VSCALE(tvec, rydir, i);
	append_tick_pair(*x1, *x2, tvec, major);
	append_tick_pair(*z1, *z2, tvec, major);
	if (!axes->pos_only) {
	    VSCALE(tvec, neg_rydir, i);
	    append_tick_pair(*x1, *x2, tvec, major);
	    append_tick_pair(*z1, *z2, tvec, major);
	}

	VSCALE(tvec, rzdir, i);
	append_tick_pair(*x1, *x2, tvec, major);
	append_tick_pair(*y1, *y2, tvec, major);
	if (!axes->pos_only) {
	    VSCALE(tvec, neg_rzdir, i);
	    append_tick_pair(*x1, *x2, tvec, major);
	    append_tick_pair(*y1, *y2, tvec, major);
	}
    }
}

static void
ged_obol_faceplate_sync_axes_one(BObolViewController *controller,
				 struct ged_view_context *view_ctx,
				 const char *prefix,
				 struct bv_axes_state axes,
				 int position_mode,
				 const mat_t supplied_rotation = NULL)
{
    const int model_axes = position_mode == 2;
    std::string line_name = std::string(prefix) + "/lines";
    std::string label_name = std::string(prefix) + "/labels";
    std::string tick_name = std::string(prefix) + "/ticks";
    std::string major_tick_name = std::string(prefix) + "/major_ticks";
    if (!axes.draw) {
	ged_obol_faceplate_remove_axis_variants(controller, view_ctx, line_name);
	ged_obol_faceplate_remove(controller, view_ctx, label_name.c_str());
	ged_obol_faceplate_remove(controller, view_ctx, tick_name.c_str());
	ged_obol_faceplate_remove(controller, view_ctx,
				  major_tick_name.c_str());
	return;
    }

    if (axes.axes_size <= 0.0)
	axes.axes_size = model_axes ? 2.0 : 0.2;
    if (ged_obol_rgb_is_zero(axes.axes_color))
	VSET(axes.axes_color, 255, 255, 255);
    if (ged_obol_rgb_is_zero(axes.label_color))
	VSET(axes.label_color, 255, 255, 0);
    if (position_mode == 1 && VNEAR_ZERO(axes.axes_pos, SMALL_FASTF)) {
	VSET(axes.axes_pos, 0.80, -0.80, 0.0);
	axes.pos_only = 1;
	axes.triple_color = 1;
	axes.label_flag = 1;
    }
    if (model_axes)
	axes.label_flag = 1;

    mat_t rmat;
    if (supplied_rotation) {
	MAT_COPY(rmat, supplied_rotation);
    } else if (!bv_rotation_get(rmat, ged_obol_bv_const(view_ctx))) {
	ged_obol_faceplate_remove_axis_variants(controller, view_ctx, line_name);
	ged_obol_faceplate_remove(controller, view_ctx, label_name.c_str());
	ged_obol_faceplate_remove(controller, view_ctx, tick_name.c_str());
	ged_obol_faceplate_remove(controller, view_ctx,
				  major_tick_name.c_str());
	return;
    }
    if (model_axes) {
	point_t map;
	mat_t model2view;
	if (bv_model2view_get(model2view, ged_obol_bv_const(view_ctx))) {
	    VSCALE(map, axes.axes_pos,
		   bv_local2base_get(ged_obol_bv_const(view_ctx)));
	    MAT4X3PNT(axes.axes_pos, model2view, map);
	}
    }

    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    const fastf_t aspect = (width > 0 && height > 0) ?
			   (fastf_t)width / (fastf_t)height : 1.0;
    if (!model_axes)
	axes.axes_pos[Y] *= (aspect > SMALL_FASTF) ? 1.0 / aspect : 1.0;

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    std::vector<BObolLabel> labels;
    points.reserve(6);
    commands.reserve(6);
    labels.reserve(3);

    if (axes.triple_color) {
	ged_obol_faceplate_remove(controller, view_ctx, line_name.c_str());
	for (int axis = X; axis <= Z; axis++) {
	    struct bv_axes_state axis_axes = axes;
	    int axis_color[3] = {0, 0, 0};
	    ged_obol_faceplate_axis_triple_color(axis, axis_color);
	    VMOVE(axis_axes.axes_color, axis_color);
	    VMOVE(axis_axes.label_color, axis_color);

	    std::vector<SbVec3f> axis_points;
	    std::vector<int32_t> axis_commands;
	    std::vector<BObolLabel> axis_labels;
	    axis_points.reserve(2);
	    axis_commands.reserve(2);
	    ged_obol_faceplate_append_axis(axis_points, axis_commands,
					   axis_labels, view_ctx, rmat, &axis_axes, axis, aspect,
					   axis == X ? "X" : axis == Y ? "Y" : "Z");

	    BObolFeatureStyle axis_style = ged_obol_faceplate_style(
						 axis_axes.axes_color, 255, 255, 255,
						 axis_axes.line_width);
	    std::string axis_name =
		line_name + ged_obol_faceplate_axis_suffix(axis);
	    (void)ged_obol_faceplate_publish_lines(controller, view_ctx,
						   axis_name.c_str(), axis_points, axis_commands,
						   axis_style);
	    labels.insert(labels.end(), axis_labels.begin(),
			  axis_labels.end());
	}
    } else {
	for (int axis = X; axis <= Z; axis++) {
	    std::string axis_name =
		line_name + ged_obol_faceplate_axis_suffix(axis);
	    ged_obol_faceplate_remove(controller, view_ctx,
				      axis_name.c_str());
	}
	ged_obol_faceplate_append_axis(points, commands, labels, view_ctx,
				       rmat, &axes, X, aspect, "X");
	ged_obol_faceplate_append_axis(points, commands, labels, view_ctx,
				       rmat, &axes, Y, aspect, "Y");
	ged_obol_faceplate_append_axis(points, commands, labels, view_ctx,
				       rmat, &axes, Z, aspect, "Z");

	BObolFeatureStyle line_style = ged_obol_faceplate_style(
					     axes.axes_color, 255, 255, 255, axes.line_width);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx,
					       line_name.c_str(), points, commands, line_style);
    }

    std::vector<SbVec3f> tick_points;
    std::vector<int32_t> tick_commands;
    std::vector<SbVec3f> major_tick_points;
    std::vector<int32_t> major_tick_commands;
    ged_obol_faceplate_append_axis_ticks(tick_points, tick_commands,
					 major_tick_points, major_tick_commands, view_ctx, rmat, &axes,
					 aspect);
    if (!tick_points.empty()) {
	BObolFeatureStyle tick_style = ged_obol_faceplate_style(
					     axes.tick_color, 255, 255, 0, axes.line_width);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx,
					       tick_name.c_str(), tick_points, tick_commands, tick_style);
    } else {
	ged_obol_faceplate_remove(controller, view_ctx, tick_name.c_str());
    }
    if (!major_tick_points.empty()) {
	BObolFeatureStyle major_tick_style = ged_obol_faceplate_style(
		axes.tick_major_color, 255, 0, 0, axes.line_width);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx,
					       major_tick_name.c_str(), major_tick_points,
					       major_tick_commands, major_tick_style);
    } else {
	ged_obol_faceplate_remove(controller, view_ctx,
				  major_tick_name.c_str());
    }

    if (!labels.empty()) {
	BObolFeatureStyle label_style = ged_obol_faceplate_style(
					      axes.label_color, 255, 255, 0, 1);
	(void)ged_obol_faceplate_publish_hud_labels(controller, view_ctx,
		label_name.c_str(), labels, label_style);
    } else {
	ged_obol_faceplate_remove(controller, view_ctx, label_name.c_str());
    }
}

static void
ged_obol_faceplate_sync_axes(BObolViewController *controller,
			     struct ged_view_context *view_ctx)
{
    struct bv_axes_state axes = BV_AXES_STATE_INIT;
    if (bv_view_axes_state_get(&axes, ged_obol_bv_const(view_ctx)))
	ged_obol_faceplate_sync_axes_one(controller, view_ctx,
					 "_faceplate/view_axes", axes, 1);
    if (bv_model_axes_state_get(&axes, ged_obol_bv_const(view_ctx)))
	ged_obol_faceplate_sync_axes_one(controller, view_ctx,
					 "_faceplate/model_axes", axes, 2);
}

extern "C" int
ged_draw_obol_view_context_hud_axes_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct bv_axes_state *axes,
    const mat_t rotation)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !view_ctx || !name || !name[0])
	return 0;

    struct bv_axes_state state = BV_AXES_STATE_INIT;
    if (axes)
	state = *axes;
    ged_obol_faceplate_sync_axes_one(controller, view_ctx, name, state, 0,
	rotation);
    return 1;
}

extern "C" int
ged_draw_obol_view_context_hud_lines_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !view_ctx || !name || !name[0])
	return 0;
    if (!points || !point_count) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return 1;
    }

    std::vector<SbVec3f> screen_points;
    screen_points.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	SbVec3f point;
	ged_obol_faceplate_to_pixel(point, view_ctx, points[i][X], points[i][Y]);
	screen_points.push_back(point);
    }
    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    /* Screen-lock: pixel-space geometry rendered through an SoHUDKit, matching
     * the HUD text labels (store_hud_wrap_if_needed). */
    obol_style.hud = TRUE;
    BObolFeatureHandle handle = ged_obol_faceplate_publish_lines(controller,
	view_ctx, name, screen_points,
	ged_obol_commands_from_ged(cmds, point_count), obol_style);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_hud_labels_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct ged_annotation_label *labels,
    size_t label_count,
    const struct ged_view_feature_style *style)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !view_ctx || !name || !name[0])
	return 0;
    if (!labels || !label_count) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return 1;
    }

    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    std::vector<BObolLabel> hud_labels;
    hud_labels.reserve(label_count);
    for (size_t i = 0; i < label_count; i++) {
	if (!labels[i].text || !labels[i].text[0])
	    continue;
	BObolLabel label;
	label.text = labels[i].text;
	label.point = SbVec3f(
	    static_cast<float>((labels[i].point[X] + 1.0) * 0.5 * width),
	    static_cast<float>((labels[i].point[Y] + 1.0) * 0.5 * height),
	    0.0f);
	label.anchor = labels[i].anchor;
	label.fontSize = static_cast<float>(labels[i].font_size > 0.0 ?
	    labels[i].font_size : 20.0);
	if (labels[i].color_valid) {
	    label.hasColor = TRUE;
	    label.color = SbColor(labels[i].color[0] / 255.0f,
		labels[i].color[1] / 255.0f,
		labels[i].color[2] / 255.0f);
	}
	hud_labels.push_back(label);
    }
    if (hud_labels.empty()) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return 1;
    }

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle = ged_obol_faceplate_publish_hud_labels(
	controller, view_ctx, name, hud_labels, obol_style);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_hud_line_layers_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct ged_annotation_line_layer *layers,
    size_t layer_count,
    const struct ged_view_feature_style *style)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !view_ctx || !name || !name[0])
	return 0;

    std::vector<BObolLineLayer> obol_layers;
    obol_layers.reserve(layer_count);
    for (size_t i = 0; layers && i < layer_count; i++) {
	if (!layers[i].points || !layers[i].point_count)
	    continue;
	BObolLineLayer layer;
	layer.name = layers[i].name ? layers[i].name : name;
	layer.style = ged_obol_feature_style_from_ged(&layers[i].style);
	layer.points.reserve(layers[i].point_count);
	layer.commands.reserve(layers[i].point_count);
	for (size_t j = 0; j < layers[i].point_count; j++) {
	    SbVec3f point;
	    ged_obol_faceplate_to_pixel(point, view_ctx,
		    layers[i].points[j][X], layers[i].points[j][Y]);
	    layer.points.push_back(point);
	    layer.commands.push_back(ged_obol_line_command_from_ged(
		layers[i].commands ? layers[i].commands[j] : -1, j));
	}
	obol_layers.push_back(layer);
    }
    if (obol_layers.empty()) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return 1;
    }

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    /* Screen-lock the faceplate lines: pixel-space geometry rendered through an
     * SoHUDKit, matching the HUD text labels (store_hud_wrap_if_needed). */
    obol_style.hud = TRUE;
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle = controller->features().publishLineLayers(
	name, BObolFeatureScope::Local, obol_layers, &obol_style, &owner);
    (void)ged_obol_feature_mark_overlay(controller, handle,
	ged_obol_faceplate_overlay_info(view_ctx));
    return handle.isValid() ? 1 : 0;
}

static void
ged_obol_faceplate_sync_framebuffer(BObolViewController *controller,
			    struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/framebuffer";
    /* Framebuffer images are owned by BObolFramebufferStream.  Remove the
     * retired display-manager snapshot feature if an older host created it. */
    ged_obol_faceplate_remove(controller, view_ctx, name);
}

static int
ged_obol_view_context_faceplate_sync(struct ged *gedp, struct ged_view_context *view_ctx,
	BObolViewController *controller)
{
    if (!gedp || !view_ctx || !controller)
	return BRLCAD_OK;

    ged_obol_faceplate_sync_center_dot(controller, view_ctx);
	ged_obol_faceplate_sync_interactive_rect(controller, view_ctx);
    ged_obol_faceplate_sync_grid(controller, view_ctx);
    ged_obol_faceplate_sync_adc(controller, view_ctx);
    /* Controller render timing is useful performance telemetry, but it is not
     * an observed frame rate.  Hosts publish their completed-presentation
     * cadence separately, including offscreen/headless hosts. */
    const uint64_t presentation_interval =
	controller->getDisplayedPresentationIntervalNanoseconds();
    if (presentation_interval)
	(void)bv_frametime_set(ged_obol_bv(view_ctx), presentation_interval);
    ged_obol_faceplate_sync_params(controller, view_ctx);
    ged_obol_faceplate_sync_scale(controller, view_ctx);
    ged_obol_faceplate_sync_axes(controller, view_ctx);
    ged_obol_faceplate_sync_framebuffer(controller, view_ctx);

	const int framebuffer_mode =
	bv_framebuffer_mode_get(ged_obol_bv_const(view_ctx));
    const int framebuffer_visible = framebuffer_mode ==
	BV_FRAMEBUFFER_MODE_OVERLAY || framebuffer_mode ==
	BV_FRAMEBUFFER_MODE_UNDERLAY || framebuffer_mode ==
	BV_FRAMEBUFFER_MODE_INTERLAY;
    (void)ged_obol_fbserv_composition_set(gedp, framebuffer_mode);
    if (framebuffer_visible)
	(void)ged_draw_obol_framebuffer_present(gedp);

    return BRLCAD_OK;
}

extern "C" int
ged_draw_obol_faceplate_sync(struct ged *gedp, struct ged_view_context *view_ctx)
{
    return ged_obol_view_context_faceplate_sync(gedp, view_ctx,
	ged_obol_view_controller_for_context(view_ctx));
}

/* First shader token == "light"? */
static bool
ged_obol_shader_is_light(const char *shader)
{
    if (!shader)
	return false;
    while (*shader == ' ' || *shader == '\t')
	shader++;
    if (bu_strncmp(shader, "light", 5) != 0)
	return false;
    const char c = shader[5];
    return (c == '\0' || c == ' ' || c == '\t' || c == ';');
}

/* Collect in-scene lights from the database's "light"-shader regions.  Done in
 * the GED layer (not the Obol geometry realize walk) so it is independent of the
 * LoD/geometry cache, which otherwise skips the walk on a warm cache.  Positions
 * are model/world-space (region bbox centers); parameters mirror sh_light.c. */
static void
ged_obol_collect_scene_lights(struct ged *gedp,
	std::vector<BObolSceneLightRealization> &out)
{
    out.clear();
    if (!gedp || !gedp->dbip)
	return;
    struct db_i *dbip = gedp->dbip;

    struct directory *dp = RT_DIR_NULL;
    FOR_ALL_DIRECTORY_START(dp, dbip) {
	    if ((dp->d_flags & RT_DIR_HIDDEN) || !(dp->d_flags & RT_DIR_COMB))
		continue;

	    struct rt_db_internal intern;
	    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
		continue;
	    struct rt_comb_internal *comb =
		(struct rt_comb_internal *)intern.idb_ptr;
	    if (!comb || comb->magic != RT_COMB_MAGIC ||
		!ged_obol_shader_is_light(bu_vls_cstr(&comb->shader))) {
		rt_db_free_internal(&intern);
		continue;
	    }

	    /* World-space light position: region geometry bbox center. */
	    struct bu_vls msgs = BU_VLS_INIT_ZERO;
	    point_t rpp_min = VINIT_ZERO;
	    point_t rpp_max = VINIT_ZERO;
	    const char *av[2] = {dp->d_namep, NULL};
	    if (rt_obj_bounds(&msgs, dbip, 1, av, 0, rpp_min, rpp_max) < 0) {
		bu_vls_free(&msgs);
		rt_db_free_internal(&intern);
		continue;
	    }
	    bu_vls_free(&msgs);
	    point_t center;
	    VADD2SCALE(center, rpp_min, rpp_max, 0.5);

	    /* Parameters (defaults mirror sh_light.c): intensity 1, angle 180
	     * (omni), not infinite. */
	    double intensity = 1.0;
	    double angle = 180.0;
	    double fraction = 1.0;
	    double target[3] = {0.0, 0.0, 0.0};
	    int infinite = 0;
	    int exaim = 0;
	    const char *shader = bu_vls_cstr(&comb->shader);
	    struct bu_vls argcopy = BU_VLS_INIT_ZERO;
	    bu_vls_strcpy(&argcopy, shader + 5); /* skip "light" */
	    for (char *s = bu_vls_addr(&argcopy); *s; s++) {
		if (*s == '=' || *s == ';' || *s == ',')
		    *s = ' ';
	    }
	    char *lav[64];
	    size_t lac = bu_argv_from_string(lav, 63, bu_vls_addr(&argcopy));
	    for (size_t k = 0; k < lac; k++) {
		const char *key = lav[k];
		if ((BU_STR_EQUAL(key, "bright") || BU_STR_EQUAL(key, "b") ||
			BU_STR_EQUAL(key, "inten")) && k + 1 < lac)
		    intensity = atof(lav[++k]);
		else if ((BU_STR_EQUAL(key, "angle") || BU_STR_EQUAL(key, "a")) &&
			k + 1 < lac)
		    angle = atof(lav[++k]);
		else if ((BU_STR_EQUAL(key, "fract") || BU_STR_EQUAL(key, "f")) &&
			k + 1 < lac)
		    fraction = atof(lav[++k]);
		else if ((BU_STR_EQUAL(key, "infinite") || BU_STR_EQUAL(key, "i")) &&
			k + 1 < lac)
		    infinite = atoi(lav[++k]);
		else if ((BU_STR_EQUAL(key, "target") || BU_STR_EQUAL(key, "t") ||
			BU_STR_EQUAL(key, "aim") || BU_STR_EQUAL(key, "d") ||
			BU_STR_EQUAL(key, "dir")) && k + 3 < lac) {
		    target[0] = atof(lav[++k]);
		    target[1] = atof(lav[++k]);
		    target[2] = atof(lav[++k]);
		    exaim = 1;
		}
	    }
	    bu_vls_free(&argcopy);

	    vect_t aim = {0.0, 0.0, -1.0};
	    if (exaim)
		VSUB2(aim, target, center);
	    VUNITIZE(aim);

	    BObolSceneLightRealization L;
	    L.name = dp->d_namep ? dp->d_namep : "";
	    if (comb->rgb_valid) {
		L.color = SbColor(comb->rgb[0] / 255.0f, comb->rgb[1] / 255.0f,
		    comb->rgb[2] / 255.0f);
	    } else {
		L.color = SbColor(1.0f, 1.0f, 1.0f);
	    }
	    double norm = intensity * (fraction >= 0.0 ? fraction : 1.0);
	    if (norm < 0.0)
		norm = 0.0;
	    if (norm > 1.0)
		norm = 1.0;
	    L.intensity = (float)norm;
	    L.direction = SbVec3f((float)aim[0], (float)aim[1], (float)aim[2]);
	    if (infinite) {
		L.kind = BOBOL_SCENE_LIGHT_DIRECTIONAL;
	    } else if (angle < 180.0) {
		L.kind = BOBOL_SCENE_LIGHT_SPOT;
		L.position = SbVec3f((float)center[0], (float)center[1],
		    (float)center[2]);
		L.coneAngleDeg = (float)angle;
	    } else {
		L.kind = BOBOL_SCENE_LIGHT_POINT;
		L.position = SbVec3f((float)center[0], (float)center[1],
		    (float)center[2]);
	    }
	    out.push_back(L);
	    rt_db_free_internal(&intern);
    } FOR_ALL_DIRECTORY_END;
}

extern "C" int
ged_draw_obol_lighting_sync(struct ged *gedp,
	struct ged_view_context *view_ctx)
{
    if (!view_ctx)
	return BRLCAD_ERROR;
    const struct bv *view = ged_obol_bv_const(view_ctx);
    if (!view)
	return BRLCAD_ERROR;
    struct bv_lighting_state lighting;
    if (!bv_lighting_state_get(&lighting, view))
	return BRLCAD_ERROR;

    /* bv_lighting_state is the persistent source of truth; push it to the live
     * Obol controller when one is attached (headless views just keep the bv
     * state).  Direction tracking must be applied before the enable so a
     * freshly-enabled headlight picks up the correct direction. */
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return BRLCAD_OK;

    controller->setHeadlightOffset(SbVec3f(
	(float)lighting.headlight_offset[0],
	(float)lighting.headlight_offset[1],
	(float)lighting.headlight_offset[2]));
    controller->setHeadlightCameraTracked(
	lighting.headlight_tracks_camera ? TRUE : FALSE);
    controller->setHeadlightEnabled(
	lighting.headlight_enabled ? TRUE : FALSE);
    /* Supply the in-scene lights from the database (cache-immune) and apply the
     * enable state. */
    std::vector<BObolSceneLightRealization> scene_lights;
    ged_obol_collect_scene_lights(gedp, scene_lights);
    controller->setSceneLights(scene_lights);
    controller->setSceneLightsEnabled(
	lighting.scene_lights_enabled ? TRUE : FALSE);
    return BRLCAD_OK;
}

extern "C" int
ged_draw_obol_shading_sync(struct ged *gedp,
	struct ged_view_context *view_ctx)
{
    (void)gedp;
    if (!view_ctx)
	return BRLCAD_ERROR;
    const struct bv *view = ged_obol_bv_const(view_ctx);
    if (!view)
	return BRLCAD_ERROR;
    struct bv_shading_state shading;
    if (!bv_shading_state_get(&shading, view))
	return BRLCAD_ERROR;

    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return BRLCAD_OK;

    BObolViewLodState::NormalStyle style =
	BObolViewLodState::NORMAL_AUTHORED;
    if (shading.normal_style == BV_NORMAL_FLAT)
	style = BObolViewLodState::NORMAL_FLAT;
    else if (shading.normal_style == BV_NORMAL_SMOOTH)
	style = BObolViewLodState::NORMAL_SMOOTH;
    controller->setNormalStyle(style,
	static_cast<float>(shading.normal_crease_angle));
    return BRLCAD_OK;
}

extern "C" int
ged_draw_obol_view_context_feature_store_active(struct ged_view_context *view_ctx)
{
    return ged_obol_view_controller_for_context(view_ctx) ? 1 : 0;
}

extern "C" size_t
ged_draw_obol_view_context_clear(struct ged_view_context *view_ctx, int flags)
{
    if (!(flags & GED_VIEW_CLEAR_VIEW))
	return 0;

    size_t removed = 0;
    if (flags & GED_VIEW_CLEAR_LOCAL) {
	BObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	if (!controller)
	    return 0;
	BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
	removed += controller->features().removeScope(
		       BOBOL_FEATURE_SCOPE_LOCAL, &owner);
	removed += controller->polygons().removeScope(
		       BOBOL_FEATURE_SCOPE_LOCAL);
	controller->selection().clear(&owner, BOBOL_SELECTION_ALL);
    } else {
	BObolViewController *controller =
	    ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
	if (!controller)
	    return 0;
	removed += controller->features().removeScope(
		       BOBOL_FEATURE_SCOPE_SHARED, NULL);
	removed += controller->polygons().removeScope(
		       BOBOL_FEATURE_SCOPE_SHARED);
	controller->selection().clear(NULL, BOBOL_SELECTION_ALL);
    }
    return removed;
}

static enum ged_draw_obol_view_feature_kind
ged_obol_view_feature_kind(BObolFeatureKind kind) {
    switch (kind)
    {
	case BObolFeatureKind::Lines:
		    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINES;
	case BObolFeatureKind::IndexedLines:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_LINES;
	case BObolFeatureKind::Points:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_POINTS;
	case BObolFeatureKind::Labels:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_LABELS;
	case BObolFeatureKind::Arrow:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_ARROW;
	case BObolFeatureKind::Axes:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_AXES;
	case BObolFeatureKind::LineLayer:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINE_LAYER;
	case BObolFeatureKind::EditPreview:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_EDIT_PREVIEW;
	case BObolFeatureKind::IndexedFaceSet:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_FACE_SET;
	case BObolFeatureKind::PolygonOverlay:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_POLYGON_OVERLAY;
	case BObolFeatureKind::HudLabel:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_HUD_LABEL;
	case BObolFeatureKind::CustomNode:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_CUSTOM_NODE;
	case BObolFeatureKind::Unknown:
	default:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_UNKNOWN;
    }
}

struct ged_obol_feature_records_visit {
    struct ged_view_context *view_ctx;
    unsigned int query_flags;
    const char *glob;
    ged_draw_obol_view_feature_record_cb cb;
    void *userdata;
    int count;
};

static BObolFeatureStyle
ged_obol_line_layer_effective_style(const BObolFeatureStyle &parent,
				    const BObolFeatureStyle &layer)
{
    BObolFeatureStyle out = parent;
    if (layer.hasVisible) {
	out.hasVisible = TRUE;
	out.visible = layer.visible;
    }
    if (layer.hasSelectable) {
	out.hasSelectable = TRUE;
	out.selectable = layer.selectable;
    }
    if (layer.hasColor) {
	out.hasColor = TRUE;
	out.color = layer.color;
    }
    if (layer.hasLineWidth) {
	out.hasLineWidth = TRUE;
	out.lineWidth = layer.lineWidth;
    }
    if (layer.hasLineStyle) {
	out.hasLineStyle = TRUE;
	out.lineStyle = layer.lineStyle;
    }
    if (layer.hasArrow) {
	out.hasArrow = TRUE;
	out.arrow = layer.arrow;
    }
    if (layer.hasArrowTip) {
	out.hasArrowTip = TRUE;
	out.arrowTipLength = layer.arrowTipLength;
	out.arrowTipWidth = layer.arrowTipWidth;
    }
    return out;
}

static int
ged_obol_feature_record_visible(const BObolFeatureStyle &parent,
				const BObolFeatureStyle *child)
{
    if (parent.hasVisible && !parent.visible)
	return 0;
    if (child && child->hasVisible && !child->visible)
	return 0;
    return 1;
}

static int
ged_obol_feature_record_glob_matches(const char *glob, const SbString &name)
{
    if (!glob || !glob[0])
	return 1;

    const char *str = name.getString();
    return (str && bu_path_match(glob, str, 0) == 0) ? 1 : 0;
}

static int
ged_obol_feature_record_emit(struct ged_obol_feature_records_visit *ctx,
			     const BObolFeatureRecord &record,
			     const SbString &name,
			     enum ged_draw_obol_view_feature_kind kind,
			     const BObolFeatureStyle &style,
			     int visible,
			     size_t point_count,
			     size_t command_count,
			     size_t index_count,
			     size_t normal_count,
			     size_t label_count,
			     size_t axes_center_count,
			     size_t child_count,
			     const char *line_layer_parent_name,
			     size_t line_layer_index)
{
    if (!ctx || !ctx->cb)
	return 0;

    if ((ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY) &&
	!visible)
	return 1;

    if (!ged_obol_feature_record_glob_matches(ctx->glob, name))
	return 1;

    struct ged_draw_obol_view_feature_record out;
    memset(&out, 0, sizeof(out));
    out.name = name.getString();
    out.kind = kind;
    out.local = record.scope == BObolFeatureScope::Local ? 1 : 0;
    out.visible = visible;
    out.realized = record.realized ? 1 : 0;
    out.color[0] = 255;
    out.color[1] = 255;
    out.color[2] = 255;
    if (style.hasColor)
	ged_obol_rgb_from_color(style.color, out.color);
    out.line_style = style.hasLineStyle ? style.lineStyle : 0;
    out.line_width = style.hasLineWidth ? style.lineWidth : 1;
    out.point_count = point_count;
    out.command_count = command_count;
    out.index_count = index_count;
    out.normal_count = normal_count;
    out.label_count = label_count;
    out.axes_center_count = axes_center_count;
    out.child_count = child_count;
    out.line_layer_parent_name = line_layer_parent_name;
    out.line_layer_index = line_layer_index;

    ctx->count++;
    return ctx->cb(&out, ctx->userdata);
}

static int
ged_obol_feature_record_visit_cb(const BObolFeatureRecord &record,
				 void *userData)
{
    struct ged_obol_feature_records_visit *ctx =
	(struct ged_obol_feature_records_visit *)userData;
    if (!ctx || !ctx->cb)
	return 0;

    const int wants_db =
	(ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS) ? 1 : 0;
    const int wants_view =
	(ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS) ? 1 : 0;
    if (wants_db && !wants_view)
	return 1;

    if (record.kind == BObolFeatureKind::LineLayer &&
	!record.layers.empty()) {
	for (size_t i = 0; i < record.layers.size(); i++) {
	    const BObolLineLayer &layer = record.layers[i];
	    BObolFeatureStyle style =
		ged_obol_line_layer_effective_style(record.style,
						    layer.style);
	    int visible = ged_obol_feature_record_visible(record.style,
			  &layer.style);
	    size_t child_count = layer.points.empty() ? 0 : 1;
	    if (!ged_obol_feature_record_emit(ctx, record,
					      layer.name.getLength() ? layer.name : record.name,
					      GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINES, style, visible,
					      layer.points.size(), layer.commands.size(), 0, 0, 0, 0,
					      child_count, record.name.getString(), i))
		return 0;
	}
	return 1;
    }

    size_t child_count = 0;
    if (!record.layers.empty())
	child_count = record.layers.size();
    else if (!record.labels.empty())
	child_count = record.labels.size();
    else if (!record.axesCenters.empty())
	child_count = record.axesCenters.size();
    else if (!record.points.empty())
	child_count = 1;
    else if (record.kind == BObolFeatureKind::CustomNode &&
	     record.realized)
	child_count = 1;

    return ged_obol_feature_record_emit(ctx, record, record.name,
					ged_obol_view_feature_kind(record.kind), record.style,
					ged_obol_feature_record_visible(record.style, NULL),
					record.points.size(), record.commands.size(),
					record.indices.size(), record.normals.size(), record.labels.size(),
					record.axesCenters.size(), child_count, NULL, 0);
}

extern "C" int
ged_draw_obol_view_context_feature_records_foreach(
    struct ged_view_context *view_ctx,
    unsigned int query_flags,
    const char *glob,
    ged_draw_obol_view_feature_record_cb cb,
    void *userdata)
{
    BObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    BObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if ((!local_controller && !shared_controller) || !cb)
	return 0;

    const int wants_db =
	(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS) ? 1 : 0;
    const int wants_view =
	(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS) ? 1 : 0;
    if (wants_db && !wants_view)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    struct ged_obol_feature_records_visit ctx;
    ctx.view_ctx = view_ctx;
    ctx.query_flags = query_flags;
    ctx.glob = glob;
    ctx.cb = cb;
    ctx.userdata = userdata;
    ctx.count = 0;
    if (local_controller) {
	local_controller->features().visitRecords(
		ged_obol_feature_record_visit_cb, &ctx,
		BOBOL_FEATURE_SCOPE_LOCAL, &owner);
    }
    if (!(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_LOCAL_ONLY) &&
	shared_controller) {
	shared_controller->features().visitRecords(
		ged_obol_feature_record_visit_cb, &ctx,
		BOBOL_FEATURE_SCOPE_SHARED, NULL);
    }
    return ctx.count;
}

extern "C" int
ged_draw_obol_view_context_indexed_face_set_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count,
    const struct ged_view_feature_style *style)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name || !points || !point_count || !indices ||
	!index_count)
	return 0;

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle =
	ged_obol_publish_indexed_face_set(controller, view_ctx, name, local,
					  ged_obol_points_from_ged(points, point_count),
					  ged_obol_vectors_from_ged(normals, normal_count),
					  ged_obol_indices_from_ged(indices, index_count),
					  style ? &obol_style : NULL);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_lines_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;
    if (!points || !point_count)
	return ged_obol_remove_feature(controller, view_ctx, name,
				       local ? 1 : 0) ? 1 : 1;

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle =
	ged_obol_publish_line_set(controller, view_ctx, name, local,
				  ged_obol_points_from_ged(points, point_count),
				  ged_obol_commands_from_ged(cmds, point_count),
				  style ? &obol_style : NULL);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_tcl_polygons_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_view_feature_style *style)
{
    if (!points || !cmds || !point_count) {
	BObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	return (controller && name) ?
	       (ged_obol_remove_feature(controller, view_ctx, name, 1) ? 1 : 1) :
	       0;
    }
    int ret = ged_draw_obol_view_context_lines_replace(view_ctx, name, 1,
	      points, cmds, point_count, style);
    if (!ret)
	return 0;

    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BObolOverlayClass::PolygonEdit,
						 BObolOverlayLifecycle::PerTool,
						 BObolOverlayOrder::PostTransparent,
						 name));
}

static BObolViewController *
ged_obol_feature_controller_from_ged_ref(ged_view_feature_ref ref)
{
    int local = 0;
    struct ged_view_context *view_ctx = ged_view_context_from_reference_owner(ref.owner, &local);
    if (!view_ctx || !ref.id || !ref.generation)
	return NULL;

    BObolViewController *controller = local ?
	ged_obol_view_controller_for_context(view_ctx) :
	ged_obol_shared_view_controller_for_context(view_ctx);
    if (!controller ||
	controller->features().referenceGeneration() != ref.generation)
	return NULL;

    BObolFeatureRecord record;
    if (!controller->features().record(BObolFeatureHandle(ref.id, 1),
	    record))
	return NULL;
    return controller;
}

static BObolFeatureHandle
ged_obol_feature_handle_from_ged_ref(ged_view_feature_ref ref)
{
    return (ref.id && ref.generation) ? BObolFeatureHandle(ref.id, 1) :
	   BObolFeatureHandle();
}

static ged_view_feature_ref
ged_obol_ged_feature_ref(struct ged_view_context *view_ctx,
			 BObolViewController *controller,
			 BObolFeatureHandle handle,
			 int local)
{
    ged_view_feature_ref ref = GED_VIEW_FEATURE_REF_NULL_INIT;
    if (!controller || !handle.isValid())
	return ref;
    ref.owner = ged_view_context_reference_owner(view_ctx, local);
    ref.id = handle.id;
    ref.generation = controller->features().referenceGeneration();
    if (!ref.owner)
	return GED_VIEW_FEATURE_REF_NULL;
    return ref;
}

static BObolOverlayInfo
ged_obol_edit_overlay_info(struct ged_view_context *view_ctx,
			      const void *owner,
			      const char *source_path,
			      int sort_order)
{
    BObolOverlayInfo overlay = ged_obol_model_overlay_info(view_ctx,
				 BObolOverlayClass::EditHandle,
				 BObolOverlayLifecycle::PerTool,
				 BObolOverlayOrder::PostTransparent,
				 source_path);
    overlay.ownerToken = owner ? owner : view_ctx;
    overlay.sortOrder = sort_order;
    return overlay;
}

static std::vector<BObolLabel>
ged_obol_labels_from_ged_feature(
    const struct ged_view_feature_label *labels,
    size_t label_count)
{
    std::vector<BObolLabel> out;
    if (!labels || !label_count)
	return out;

    out.reserve(label_count);
    for (size_t i = 0; i < label_count; i++) {
	BObolLabel label;
	label.text = labels[i].text ? labels[i].text : "";
	label.point = SbVec3f(
			  static_cast<float>(labels[i].point[X]),
			  static_cast<float>(labels[i].point[Y]),
			  static_cast<float>(labels[i].point[Z]));
	if (labels[i].color_valid) {
	    label.hasColor = TRUE;
	    label.color = ged_obol_color_from_rgb(labels[i].color);
	}
	if (labels[i].font_size > 0.0)
	    label.fontSize = static_cast<float>(labels[i].font_size);
	out.push_back(label);
    }
    return out;
}

static unsigned char ged_obol_rgb_byte_from_int(int value);

extern "C" int
ged_draw_obol_view_feature_ref_is_null(ged_view_feature_ref ref)
{
    return (!ref.owner || !ref.id || !ref.generation) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_remove_ref(ged_view_feature_ref ref)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(ref);
    return controller ?
	   (controller->features().remove(
		ged_obol_feature_handle_from_ged_ref(ref)) ? 1 : 0) : 0;
}

extern "C" int
ged_draw_obol_view_context_edit_preview_publish_event(
    struct ged_view_context *UNUSED(view_ctx),
    ged_view_feature_ref feature,
    enum ged_view_edit_preview_event event,
    const char *source_path)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(feature);
    if (!controller)
	return 0;

    BObolFeatureHandle handle = ged_obol_feature_handle_from_ged_ref(feature);
    BObolFeatureRecord record;
    if (!controller->features().record(handle, record) ||
	record.kind != BObolFeatureKind::EditPreview)
	return 0;

    if (source_path && source_path[0] && record.overlay.isOverlay) {
	BObolOverlayInfo overlay = record.overlay;
	if (overlay.sourcePath != source_path) {
	    overlay.sourcePath = source_path;
	}
	(void)controller->features().setOverlayInfo(handle, overlay);
    }

    switch (event) {
	case GED_VIEW_EDIT_PREVIEW_COMMIT:
	case GED_VIEW_EDIT_PREVIEW_CANCEL:
	case GED_VIEW_EDIT_PREVIEW_DISCARD:
	    return controller->features().remove(handle) ? 1 : 0;
	case GED_VIEW_EDIT_PREVIEW_BEGIN:
	case GED_VIEW_EDIT_PREVIEW_UPDATE:
	case GED_VIEW_EDIT_PREVIEW_REPLACE_SOURCE:
	default:
	    break;
    }

    {
	return controller->features().touch(handle) ? 1 : 0;
    }
}

extern "C" ged_view_feature_ref
ged_draw_obol_view_context_feature_overlay_ensure(
    struct ged_view_context *view_ctx,
    const char *name,
    const void *owner,
    const char *source_path)
{
    BObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller || !name || !name[0])
	return GED_VIEW_FEATURE_REF_NULL;

    BObolFeatureOwner feature_owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle = controller->features().findOwned(name,
				  BOBOL_FEATURE_SCOPE_LOCAL, &feature_owner);

    BObolFeatureRecord record;
    const int needs_publish = !handle.isValid() ||
			      !controller->features().record(handle, record) ||
			      record.kind != BObolFeatureKind::EditPreview;
    if (needs_publish) {
	std::vector<SbVec3f> points;
	std::vector<int32_t> commands;
	handle = controller->features().publishEditPreview(name,
		 source_path && source_path[0] ? source_path : name,
		 name,
		 "edit-handle",
		 points,
		 commands,
		 0,
		 0,
		 &feature_owner);
    }

    if (!handle.isValid())
	return GED_VIEW_FEATURE_REF_NULL;

    (void)controller->features().setOverlayInfo(handle,
	    ged_obol_edit_overlay_info(view_ctx, owner, source_path, 0));
    return ged_obol_ged_feature_ref(view_ctx, controller, handle, 1);
}

extern "C" ged_view_feature_ref
ged_draw_obol_view_context_feature_label_ensure(
    struct ged_view_context *view_ctx,
    const char *name,
    const void *owner)
{
    BObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller || !name || !name[0])
	return GED_VIEW_FEATURE_REF_NULL;

    BObolFeatureOwner feature_owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle = controller->features().findOwned(name,
				  BOBOL_FEATURE_SCOPE_LOCAL, &feature_owner);

    BObolFeatureRecord record;
    if (!handle.isValid() ||
	!controller->features().record(handle, record) ||
	record.kind != BObolFeatureKind::Labels) {
	std::vector<BObolLabel> labels;
	handle = controller->features().publishLabels(name,
		 BObolFeatureScope::Local, labels, NULL, &feature_owner);
    }

    if (!handle.isValid())
	return GED_VIEW_FEATURE_REF_NULL;

    (void)controller->features().setOverlayInfo(handle,
	    ged_obol_edit_overlay_info(view_ctx, owner, name, 1));
    return ged_obol_ged_feature_ref(view_ctx, controller, handle, 1);
}

extern "C" int
ged_draw_obol_view_feature_set_context(ged_view_feature_ref ref,
				       struct ged_view_context *view_ctx)
{
    int local = 0;
    void *owner_ctx = ged_view_context_from_reference_owner(ref.owner,
	&local);
    return (view_ctx && owner_ctx == view_ctx &&
	ged_obol_feature_controller_from_ged_ref(ref)) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_set_visible(ged_view_feature_ref ref,
				       int visible)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(ref);
    return controller ?
	   (controller->features().setVisible(
		ged_obol_feature_handle_from_ged_ref(ref),
		visible ? TRUE : FALSE) ? 1 : 0) : 0;
}

extern "C" int
ged_draw_obol_view_feature_set_color(ged_view_feature_ref ref,
				     int r,
				     int g,
				     int b)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(ref);
    if (!controller)
	return 0;

    const unsigned char rgb[3] = {
	ged_obol_rgb_byte_from_int(r),
	ged_obol_rgb_byte_from_int(g),
	ged_obol_rgb_byte_from_int(b)
    };
    return controller->features().setColor(
	       ged_obol_feature_handle_from_ged_ref(ref),
	       ged_obol_color_from_rgb(rgb)) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_touch(ged_view_feature_ref ref)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(ref);
    if (!controller)
	return 0;

    return controller->features().touch(
	ged_obol_feature_handle_from_ged_ref(ref)) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_labels_replace(
    ged_view_feature_ref ref,
    const struct ged_view_feature_label *labels,
    size_t label_count)
{
    if (label_count && !labels)
	return 0;

    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(ref);
    if (!controller)
	return 0;

    BObolFeatureRecord record;
    BObolFeatureHandle handle = ged_obol_feature_handle_from_ged_ref(ref);
    if (!controller->features().record(handle, record))
	return 0;

    std::vector<BObolLabel> obol_labels =
	ged_obol_labels_from_ged_feature(labels, label_count);
    if (record.kind == BObolFeatureKind::Labels ||
	record.kind == BObolFeatureKind::HudLabel)
	return controller->features().replaceLabels(handle,
		obol_labels) ? 1 : 0;

    BObolFeatureOwner owner = record.owner;
    BObolFeatureHandle labels_handle =
	controller->features().publishLabels(record.name,
	    record.scope, obol_labels, &record.style,
	    record.scope == BObolFeatureScope::Local ? &owner : NULL);
    return labels_handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_points_replace(
    ged_view_feature_ref ref,
    enum ged_view_feature_family family,
    const point_t *points,
    const int *cmds,
    size_t point_count)
{
    if (point_count && !points)
	return 0;

    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(ref);
    if (!controller)
	return 0;

    BObolFeatureRecord record;
    BObolFeatureHandle handle = ged_obol_feature_handle_from_ged_ref(ref);
    if (!controller->features().record(handle, record))
	return 0;

    std::vector<SbVec3f> obol_points =
	ged_obol_points_from_ged(points, point_count);
    std::vector<int32_t> obol_commands =
	ged_obol_commands_from_ged(cmds, point_count);
    if (family == GED_VIEW_FEATURE_TRANSIENT_PREVIEW ||
	record.kind == BObolFeatureKind::EditPreview) {
	if (record.kind == BObolFeatureKind::EditPreview)
	    return controller->features().replaceEditPreviewGeometry(
		       handle,
		       record.name,
		       obol_points,
		       obol_commands) ? 1 : 0;

	BObolFeatureOwner owner = record.owner;
	BObolFeatureHandle preview_handle =
	    controller->features().publishEditPreview(record.name,
		record.name,
		record.name,
		"edit-handle",
		obol_points,
		obol_commands,
		0,
		0,
		record.scope == BObolFeatureScope::Local ? &owner : NULL);
	return preview_handle.isValid() ? 1 : 0;
    }

    BObolFeatureOwner owner = record.owner;
    BObolFeatureHandle line_handle =
	controller->features().publishLineSet(record.name,
	    record.scope, obol_points, obol_commands, &record.style,
	    record.scope == BObolFeatureScope::Local ? &owner : NULL);
    return line_handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_edit_preview_replace(
    ged_view_feature_ref ref,
    const char *source_path,
    const char *edit_intent_id,
    const char *edit_intent_role,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    uint32_t source_revision,
    uint32_t inputs_revision)
{
    if (point_count && !points)
	return 0;

    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(ref);
    if (!controller)
	return 0;

    BObolFeatureRecord record;
    BObolFeatureHandle handle = ged_obol_feature_handle_from_ged_ref(ref);
    if (!controller->features().record(handle, record))
	return 0;

    std::vector<SbVec3f> obol_points =
	ged_obol_points_from_ged(points, point_count);
    std::vector<int32_t> obol_commands =
	ged_obol_commands_from_ged(cmds, point_count);

    BObolFeatureOwner owner = record.owner;
    const SbString preview_name = record.name;
    const SbString preview_path =
	(source_path && source_path[0]) ? SbString(source_path) :
	record.identity.getLength() > 0 ? record.identity : record.name;
    const SbString intent_id =
	(edit_intent_id && edit_intent_id[0]) ? SbString(edit_intent_id) :
	record.editIntentId.getLength() > 0 ? record.editIntentId : record.name;
    const SbString intent_role =
	(edit_intent_role && edit_intent_role[0]) ? SbString(edit_intent_role) :
	record.editIntentRole.getLength() > 0 ? record.editIntentRole :
	SbString("preview");
    const uint32_t next_source_revision =
	source_revision ? source_revision :
	record.sourceRevision ? record.sourceRevision + 1 : 0;
    const uint32_t next_inputs_revision =
	inputs_revision ? inputs_revision :
	record.inputsRevision ? record.inputsRevision + 1 : 0;

    BObolFeatureHandle preview_handle =
	controller->features().publishEditPreview(preview_name,
	    preview_path,
	    intent_id,
	    intent_role,
	    obol_points,
	    obol_commands,
	    next_source_revision,
	    next_inputs_revision,
	    record.scope == BObolFeatureScope::Local ? &owner : NULL);
    return preview_handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_clear_geometry(ged_view_feature_ref ref)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(ref);
    return controller ?
	   (controller->features().clearGeometry(
		ged_obol_feature_handle_from_ged_ref(ref)) ? 1 : 0) : 0;
}

static unsigned char
ged_obol_rgb_byte_from_int(int value)
{
    if (value < 0)
	return 0;
    if (value > 255)
	return 255;
    return static_cast<unsigned char>(value);
}

static BObolPolygonType
ged_obol_polygon_type_from_ged(int type)
{
    switch (type) {
	case GED_POLYGON_CIRCLE:
	    return BObolPolygonType::Circle;
	case GED_POLYGON_ELLIPSE:
	    return BObolPolygonType::Ellipse;
	case GED_POLYGON_RECTANGLE:
	    return BObolPolygonType::Rectangle;
	case GED_POLYGON_SQUARE:
	    return BObolPolygonType::Square;
	default:
	    return BObolPolygonType::General;
    }
}

static BObolPolygonUpdate
ged_obol_polygon_update_from_ged(int op)
{
    switch (op) {
	case GED_POLYGON_UPDATE_PROPS_ONLY:
	    return BObolPolygonUpdate::PropsOnly;
	case GED_POLYGON_UPDATE_PT_SELECT:
	    return BObolPolygonUpdate::PointSelect;
	case GED_POLYGON_UPDATE_PT_SELECT_CLEAR:
	    return BObolPolygonUpdate::PointSelectClear;
	case GED_POLYGON_UPDATE_PT_MOVE:
	    return BObolPolygonUpdate::PointMove;
	case GED_POLYGON_UPDATE_PT_APPEND:
	    return BObolPolygonUpdate::PointAppend;
	default:
	    return BObolPolygonUpdate::Default;
    }
}

static BObolFeatureScope
ged_obol_polygon_scope(int local)
{
    return local ? BObolFeatureScope::Local : BObolFeatureScope::Shared;
}

static BObolViewController *
ged_obol_polygon_controller_from_ged_ref(ged_polygon_ref ref)
{
    int local = 0;
    struct ged_view_context *view_ctx = ged_view_context_from_reference_owner(ref.owner, &local);
    if (!view_ctx || !ref.id || !ref.generation)
	return NULL;

    BObolViewController *controller = local ?
	ged_obol_view_controller_for_context(view_ctx) :
	ged_obol_shared_view_controller_for_context(view_ctx);
    if (!controller ||
	controller->polygons().referenceGeneration() != ref.generation)
	return NULL;

    BObolPolygonRecord record;
    if (!controller->polygons().record(BObolPolygonHandle(ref.id, 1),
	    record))
	return NULL;
    return controller;
}

static BObolPolygonHandle
ged_obol_polygon_handle_from_ged_ref(ged_polygon_ref ref)
{
    return (ref.id && ref.generation) ? BObolPolygonHandle(ref.id, 1) :
	   BObolPolygonHandle();
}

static ged_polygon_ref
ged_obol_ged_polygon_ref(struct ged_view_context *view_ctx,
			 BObolViewController *controller,
			 BObolPolygonHandle handle,
			 int local)
{
    ged_polygon_ref ged_ref = GED_POLYGON_REF_NULL_INIT;
    if (!controller || !handle.isValid())
	return ged_ref;
    ged_ref.owner = ged_view_context_reference_owner(view_ctx, local);
    ged_ref.id = handle.id;
    ged_ref.generation = controller->polygons().referenceGeneration();
    if (!ged_ref.owner)
	return GED_POLYGON_REF_NULL;
    return ged_ref;
}

static ged_polygon_ref
ged_obol_ged_polygon_ref_for_handle(struct ged_view_context *view_ctx,
				      BObolViewController *controller,
				      BObolPolygonHandle handle)
{
    BObolPolygonRecord record;
    if (!controller || !controller->polygons().record(handle, record))
	return GED_POLYGON_REF_NULL;
    return ged_obol_ged_polygon_ref(view_ctx, controller, handle,
	record.scope == BObolFeatureScope::Local ? 1 : 0);
}

static void
ged_obol_polygon_point_from_ged(SbVec3f &dst, const point_t src)
{
    dst = SbVec3f(static_cast<float>(src[X]),
		  static_cast<float>(src[Y]),
		  static_cast<float>(src[Z]));
}

static void
ged_obol_polygon_project_point(point_t dst, struct ged_view_context *view_ctx, const point_t src,
			       plane_t *view_plane)
{
    VMOVE(dst, src);
    if (!view_plane)
	return;
    HSET(*view_plane, 0.0, 0.0, 1.0, src[Z]);

    if (bv_plane_get(view_plane, ged_obol_bv_const(view_ctx)) != 0)
	return;

    fastf_t fx = 0.0;
    fastf_t fy = 0.0;
    bg_plane_closest_pt(&fx, &fy, view_plane, (point_t *)src);

    point_t projected = VINIT_ZERO;
    bg_plane_pt_at(&projected, view_plane, fx, fy);
    VMOVE(dst, projected);
}

static SbColor
ged_obol_polygon_color_from_bu(const struct bu_color *color)
{
    fastf_t rgb[3] = {0.0, 0.0, 0.0};
    if (color)
	bu_color_to_rgb_floats(color, rgb);
    return SbColor(static_cast<float>(rgb[0]),
		   static_cast<float>(rgb[1]),
		   static_cast<float>(rgb[2]));
}

static void
ged_obol_polygon_color_to_bu(struct bu_color *dst, const SbColor &src)
{
    if (!dst)
	return;
    fastf_t rgb[3] = {
	static_cast<fastf_t>(src[0]),
	static_cast<fastf_t>(src[1]),
	static_cast<fastf_t>(src[2])
    };
    bu_color_from_rgb_floats(dst, rgb);
}

static int
ged_obol_polygon_record_to_ged(
    BObolViewController *controller,
    ged_polygon_ref ref,
    const BObolPolygonRecord &src,
    struct ged_polygon_record *dst)
{
    if (!controller || !dst)
	return 0;

    const char *name = controller->polygons().name(src.handle);
    if (!name)
	return 0;

    memset(dst, 0, sizeof(*dst));
    dst->ref = ref;
    dst->name = name;
    dst->type = static_cast<int>(src.type);
    dst->fill_flag = src.fill ? 1 : 0;
    V2SET(dst->fill_dir, src.fillSlope[0], src.fillSlope[1]);
    dst->fill_delta = src.fillSpacing;
    ged_obol_polygon_color_to_bu(&dst->fill_color, src.fillColor);
    ged_obol_rgb_from_color(src.edgeColor, dst->edge_color);
    dst->curr_contour_i = src.currentContour;
    dst->curr_point_i = src.currentPoint;
    dst->first_contour_open = src.firstContourOpen ? 1 : 0;
    dst->contour_count = src.contourCount;
    dst->point_count = src.pointCount;
    VSET(dst->origin_point, src.originPoint[0], src.originPoint[1],
	 src.originPoint[2]);
    HMOVE(dst->vp, src.viewPlane);
    dst->vZ = src.viewZ;
    dst->user_data = src.userData;
    return 1;
}

static BObolPolygonHandle
ged_obol_polygon_create_named(
    BObolViewController *controller,
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    int type,
    const point_t input_point)
{
    if (!controller || !name || !name[0] || !input_point)
	return BObolPolygonHandle();

    point_t origin = VINIT_ZERO;
    plane_t view_plane = HINIT_ZERO;
    ged_obol_polygon_project_point(origin, view_ctx, input_point,
				   &view_plane);

    SbVec3f obol_origin;
    ged_obol_polygon_point_from_ged(obol_origin, origin);
    return controller->polygons().create(name,
					 ged_obol_polygon_scope(local),
					 ged_obol_polygon_type_from_ged(type),
					 obol_origin,
					 view_plane,
					 0.0f);
}

extern "C" ged_polygon_ref
ged_draw_obol_view_context_polygon_find(struct ged_view_context *view_ctx, const char *name)
{
    if (!name)
	return GED_POLYGON_REF_NULL;

    BObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (local_controller) {
	BObolPolygonHandle local =
	    local_controller->polygons().find(name,
					      BOBOL_FEATURE_SCOPE_LOCAL);
	if (local.isValid())
	    return ged_obol_ged_polygon_ref(view_ctx, local_controller, local,
		    1);
    }

    BObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if (!shared_controller)
	return GED_POLYGON_REF_NULL;
    return ged_obol_ged_polygon_ref(view_ctx, shared_controller,
				    shared_controller->polygons().find(name,
					    BOBOL_FEATURE_SCOPE_SHARED), 0);
}

extern "C" ged_polygon_ref
ged_draw_obol_view_context_polygon_find_scoped(
    struct ged_view_context *view_ctx,
    const char *name,
    int local_only)
{
    if (!name)
	return GED_POLYGON_REF_NULL;

    if (local_only) {
	BObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	if (!controller)
	    return GED_POLYGON_REF_NULL;
	return ged_obol_ged_polygon_ref(view_ctx, controller,
					controller->polygons().find(name, BOBOL_FEATURE_SCOPE_LOCAL), 1);
    }

    return ged_draw_obol_view_context_polygon_find(view_ctx, name);
}

extern "C" ged_polygon_ref
ged_draw_obol_view_context_polygon_create(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    int type,
    const point_t screen_point)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    return ged_obol_ged_polygon_ref(view_ctx, controller,
				    ged_obol_polygon_create_named(controller, view_ctx, name, local,
					    type, screen_point), local);
}

extern "C" ged_polygon_ref
ged_draw_obol_view_context_polygon_import_sketch(
    const char *name,
    struct db_i *dbip,
    struct directory *dp,
    struct ged_view_context *view_ctx,
    int local)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return GED_POLYGON_REF_NULL;
    return ged_obol_ged_polygon_ref(view_ctx, controller,
				    controller->polygons().importSketch(name,
					    ged_obol_polygon_scope(local), dbip, dp), local);
}

extern "C" int
ged_draw_obol_view_context_polygon_set_current(
    ged_polygon_ref ref,
    long contour_i,
    long point_i)
{
    ged_polygon_ref ged_ref = ref;
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ged_ref);
    if (!controller)
	return 0;
    return controller->polygons().setCurrent(
	       ged_obol_polygon_handle_from_ged_ref(ged_ref),
	       contour_i, point_i) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_polygon_set_contour_open(
    ged_polygon_ref ref,
    long contour_i,
    int open)
{
    ged_polygon_ref ged_ref = ref;
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ged_ref);
    if (!controller)
	return 0;
    return controller->polygons().setContourOpen(
	       ged_obol_polygon_handle_from_ged_ref(ged_ref),
	       contour_i, open ? TRUE : FALSE) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_polygon_area(
    ged_polygon_ref ref,
    struct ged_view_context *view_ctx,
    fastf_t *area)
{
    if (area)
	*area = 0.0;
    if (!area || !view_ctx)
	return 0;

    ged_polygon_ref ged_ref = ref;
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ged_ref);
    if (!controller)
	return 0;

    *area = controller->polygons().area(
	ged_obol_polygon_handle_from_ged_ref(ged_ref),
		bv_scale_get(ged_obol_bv_const(view_ctx)));
    return 1;
}

extern "C" int
ged_draw_obol_view_context_polygon_overlap(
    ged_polygon_ref ref,
    struct ged_view_context *view_ctx,
    const char *other_name,
    const struct bn_tol *tol,
    int *overlap)
{
    if (overlap)
	*overlap = 0;
    if (!view_ctx || !other_name || !tol || !overlap)
	return 0;

    ged_polygon_ref ged_ref = ref;
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ged_ref);
    if (!controller)
	return 0;

    BObolPolygonHandle other =
	controller->polygons().find(other_name);
    if (!other.isValid())
	return 0;

    *overlap = controller->polygons().overlaps(
		   ged_obol_polygon_handle_from_ged_ref(ged_ref),
		   other,
		   *tol,
		   bv_scale_get(ged_obol_bv_const(view_ctx))) ? 1 : 0;
    return 1;
}

static size_t
ged_draw_obol_polygon_snap_count(
    struct ged_view_context *view_ctx,
    ged_polygon_ref exclude,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;
    BObolPolygonHandle exclude_handle;
    if (ged_obol_polygon_controller_from_ged_ref(exclude) == controller)
	exclude_handle = ged_obol_polygon_handle_from_ged_ref(exclude);
    return controller->polygons().snapCount(exclude_handle);
}

static int
ged_draw_obol_polygon_clear_point_selection(
    struct ged_view_context *view_ctx,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    return controller ?
	   (controller->polygons().clearAllPointSelections() ? 1 : 0) : 0;
}

static int
ged_draw_obol_polygon_update(
    ged_polygon_ref ref,
    struct ged_view_context *UNUSED(view_ctx),
    int utype,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    if (!controller)
	return 0;
    return controller->polygons().update(
	       ged_obol_polygon_handle_from_ged_ref(ref),
	       ged_obol_polygon_update_from_ged(utype)) ? 1 : 0;
}

static int
ged_draw_obol_polygon_update_screen_pt(
    ged_polygon_ref ref,
    struct ged_view_context *view_ctx,
    int x,
    int y,
    int utype,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    if (!controller || !view_ctx)
	return 0;

    point_t model_point = VINIT_ZERO;
    if (!bv_screen_to_model(model_point, ged_obol_bv_const(view_ctx),
				       (fastf_t)x, (fastf_t)y))
	return 0;

    return controller->polygons().updateModelPoint(
	       ged_obol_polygon_handle_from_ged_ref(ref),
	       SbVec3f(static_cast<float>(model_point[X]),
		       static_cast<float>(model_point[Y]),
		       static_cast<float>(model_point[Z])),
	       ged_obol_polygon_update_from_ged(utype)) ? 1 : 0;
}

static int
ged_draw_obol_polygon_move(
    ged_polygon_ref ref,
    point_t *current_point,
    point_t *previous_point,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    if (!controller || !current_point || !previous_point)
	return 0;
    return controller->polygons().move(
	       ged_obol_polygon_handle_from_ged_ref(ref),
	       SbVec3f(static_cast<float>((*current_point)[X]),
		       static_cast<float>((*current_point)[Y]),
		       static_cast<float>((*current_point)[Z])),
	       SbVec3f(static_cast<float>((*previous_point)[X]),
		       static_cast<float>((*previous_point)[Y]),
		       static_cast<float>((*previous_point)[Z]))) ? 1 : 0;
}

static int
ged_draw_obol_polygon_set_name(
    ged_polygon_ref ref,
    const char *name,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    if (!controller || !name)
	return 0;
    return controller->polygons().rename(
	       ged_obol_polygon_handle_from_ged_ref(ref), name) ? 1 : 0;
}

static int
ged_draw_obol_polygon_set_visual(
    ged_polygon_ref ref,
    const struct bu_color *edge_color,
    const struct bu_color *fill_color,
    fastf_t fill_slope_x,
    fastf_t fill_slope_y,
    fastf_t fill_density,
    fastf_t vZ,
    int fill_flag,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    if (!controller)
	return 0;

    BObolPolygonVisual visual;
    (void)controller->polygons().visual(
	ged_obol_polygon_handle_from_ged_ref(ref), visual);
    if (edge_color)
	visual.edgeColor = ged_obol_polygon_color_from_bu(edge_color);
    if (fill_color)
	visual.fillColor = ged_obol_polygon_color_from_bu(fill_color);
    visual.fillSlope = SbVec2f(static_cast<float>(fill_slope_x),
			       static_cast<float>(fill_slope_y));
    visual.fillSpacing = static_cast<float>(fill_density);
    visual.viewZ = static_cast<float>(vZ);
    if (fill_flag)
	visual.fillFlags |= BOBOL_POLYGON_FILL_HATCH;
    else
	visual.fillFlags &= ~BOBOL_POLYGON_FILL_HATCH;
    visual.fill = (visual.fillFlags & BOBOL_POLYGON_FILL_HATCH) ?
		  TRUE : FALSE;
    return controller->polygons().setVisual(
	       ged_obol_polygon_handle_from_ged_ref(ref), visual) ? 1 : 0;
}

static int
ged_draw_obol_polygon_set_open(
    ged_polygon_ref ref,
    int open,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    if (!controller)
	return 0;
    return controller->polygons().setAllContoursOpen(
	       ged_obol_polygon_handle_from_ged_ref(ref),
	       open ? TRUE : FALSE) ? 1 : 0;
}

static int
ged_draw_obol_polygon_close(
    ged_polygon_ref ref,
    void *data)
{
    return ged_draw_obol_polygon_set_open(ref, 0, data);
}

static int
ged_draw_obol_polygon_clear_selected_point(
    ged_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    if (!controller)
	return 0;
    return controller->polygons().clearSelectedPoint(
	       ged_obol_polygon_handle_from_ged_ref(ref)) ? 1 : 0;
}

static int
ged_draw_obol_polygon_remove(
    ged_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    if (!controller)
	return 0;
    return controller->polygons().remove(
	       ged_obol_polygon_handle_from_ged_ref(ref)) ? 1 : 0;
}

static void *
ged_draw_obol_polygon_user_data(
    ged_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    return controller ?
	   controller->polygons().userData(
	       ged_obol_polygon_handle_from_ged_ref(ref)) : NULL;
}

static int
ged_draw_obol_polygon_user_data_set(
    ged_polygon_ref ref,
    void *user_data,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    if (!controller)
	return 0;
    return controller->polygons().setUserData(
	       ged_obol_polygon_handle_from_ged_ref(ref), user_data) ? 1 : 0;
}

static int
ged_draw_obol_polygon_csg(
    ged_polygon_ref target,
    ged_polygon_ref stencil,
    bg_clip_t op,
    void *UNUSED(data))
{
    BObolViewController *target_controller =
	ged_obol_polygon_controller_from_ged_ref(target);
    BObolViewController *stencil_controller =
	ged_obol_polygon_controller_from_ged_ref(stencil);
    if (!target_controller || target_controller != stencil_controller)
	return 0;
    return target_controller->polygons().csg(
	       ged_obol_polygon_handle_from_ged_ref(target),
	       ged_obol_polygon_handle_from_ged_ref(stencil),
	       op) ? 1 : 0;
}

static struct directory *
ged_draw_obol_polygon_export_sketch(
    struct db_i *dbip,
    const char *name,
    ged_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ref);
    if (!controller || !dbip || !name)
	return NULL;
    return controller->polygons().exportSketch(
	       ged_obol_polygon_handle_from_ged_ref(ref), dbip, name) ?
	   db_lookup(dbip, name, LOOKUP_QUIET) : NULL;
}

static int
ged_draw_obol_polygon_snap_exclude_set(
    struct ged_view_context *view_ctx,
    ged_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller ||
	ged_obol_polygon_controller_from_ged_ref(ref) != controller)
	return 0;
    return controller->polygons().setSnapExclude(
	       ged_obol_polygon_handle_from_ged_ref(ref)) ? 1 : 0;
}

extern "C" ged_polygon_ref
ged_draw_obol_view_context_polygon_select(
    struct ged_view_context *view_ctx,
    const point_t model_point)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !model_point)
	return GED_POLYGON_REF_NULL;
    return ged_obol_ged_polygon_ref_for_handle(view_ctx, controller,
	    controller->polygons().selectAtModelPoint(SbVec3f(
		    static_cast<float>(model_point[X]),
		    static_cast<float>(model_point[Y]),
		    static_cast<float>(model_point[Z]))));
}

extern "C" ged_polygon_ref
ged_draw_obol_view_context_polygon_dup(
    struct ged_view_context *view_ctx,
    const char *name,
    const char *new_name)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !new_name)
	return GED_POLYGON_REF_NULL;
    BObolPolygonHandle src = controller->polygons().find(name);
    return ged_obol_ged_polygon_ref_for_handle(view_ctx, controller,
	    controller->polygons().duplicate(src, new_name));
}

extern "C" int
ged_draw_obol_view_polygon_record_get(
    ged_polygon_ref ref,
    struct ged_polygon_record *record)
{
    ged_polygon_ref ged_ref = ref;
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(ged_ref);
    if (!controller || !record)
	return 0;

    BObolPolygonRecord obol_record;
    if (!controller->polygons().record(
	    ged_obol_polygon_handle_from_ged_ref(ged_ref), obol_record))
	return 0;
    return ged_obol_polygon_record_to_ged(controller, ref,
	    obol_record, record);
}

struct ged_obol_polygon_ged_visit_context {
    struct ged_view_context *view_ctx;
    BObolViewController *controller;
    ged_polygon_record_cb callback;
    void *data;
    BObolFeatureScope scope;
};

static int
ged_draw_obol_polygon_visit_ged_cb(const BObolPolygonRecord &obol_record,
				   void *data)
{
    ged_obol_polygon_ged_visit_context *ctx =
	static_cast<ged_obol_polygon_ged_visit_context *>(data);
    if (!ctx || !ctx->controller || !ctx->callback)
	return 0;
    if (obol_record.scope != ctx->scope)
	return 1;

    ged_polygon_ref ref =
	ged_obol_ged_polygon_ref(ctx->view_ctx, ctx->controller,
		obol_record.handle,
		obol_record.scope == BObolFeatureScope::Local ? 1 : 0);
    struct ged_polygon_record record;
    if (!ged_obol_polygon_record_to_ged(ctx->controller, ref, obol_record,
	    &record))
	return 1;
    return ctx->callback(ref, &record, ctx->data);
}

extern "C" void
ged_draw_obol_view_context_polygon_visit_records(
    struct ged_view_context *view_ctx,
    ged_polygon_record_cb callback,
    void *data)
{
    BObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    BObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if ((!local_controller && !shared_controller) || !callback)
	return;

    ged_obol_polygon_ged_visit_context ctx;
    ctx.view_ctx = view_ctx;
    ctx.callback = callback;
    ctx.data = data;
    if (local_controller) {
	ctx.controller = local_controller;
	ctx.scope = BObolFeatureScope::Local;
	local_controller->polygons().visitRecords(
		ged_draw_obol_polygon_visit_ged_cb, &ctx);
    }
    if (shared_controller) {
	ctx.controller = shared_controller;
	ctx.scope = BObolFeatureScope::Shared;
	shared_controller->polygons().visitRecords(
		ged_draw_obol_polygon_visit_ged_cb, &ctx);
    }
}

extern "C" size_t
ged_draw_obol_view_context_polygon_snap_count(
    struct ged_view_context *view_ctx,
    ged_polygon_ref exclude)
{
    return ged_draw_obol_polygon_snap_count(view_ctx,
	    exclude, NULL);
}

extern "C" int
ged_draw_obol_view_context_polygon_clear_point_selection(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_polygon_clear_point_selection(view_ctx, NULL);
}

extern "C" int
ged_draw_obol_view_context_polygon_snap_exclude_set(
    struct ged_view_context *view_ctx,
    ged_polygon_ref ref)
{
    return ged_draw_obol_polygon_snap_exclude_set(view_ctx,
	    ref, NULL);
}

extern "C" struct directory *
ged_draw_obol_view_polygon_export_sketch(
    struct db_i *dbip,
    const char *name,
    ged_polygon_ref ref)
{
    return ged_draw_obol_polygon_export_sketch(dbip, name,
	    ref, NULL);
}

extern "C" int
ged_draw_obol_view_context_polygon_update(
    ged_polygon_ref ref,
    struct ged_view_context *view_ctx,
    int op)
{
    return ged_draw_obol_polygon_update(
	    ref, view_ctx, op, NULL);
}

extern "C" int
ged_draw_obol_view_context_polygon_update_screen_pt(
    ged_polygon_ref ref,
    struct ged_view_context *view_ctx,
    int x,
    int y,
    int op)
{
    return ged_draw_obol_polygon_update_screen_pt(
	    ref, view_ctx, x, y, op, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_move(
    ged_polygon_ref ref,
    point_t *current_point,
    point_t *previous_point)
{
    return ged_draw_obol_polygon_move(
	    ref, current_point,
	    previous_point, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_set_name(
    ged_polygon_ref ref,
    const char *name)
{
    return ged_draw_obol_polygon_set_name(
	    ref, name, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_set_visual(
    ged_polygon_ref ref,
    const struct bu_color *edge_color,
    const struct bu_color *fill_color,
    fastf_t fill_slope_x,
    fastf_t fill_slope_y,
    fastf_t fill_density,
    fastf_t vZ,
    int fill_flag)
{
    return ged_draw_obol_polygon_set_visual(
	    ref, edge_color, fill_color,
	    fill_slope_x, fill_slope_y, fill_density, vZ, fill_flag, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_set_all_contours_open(
    ged_polygon_ref ref,
    int open)
{
    return ged_draw_obol_polygon_set_open(
	    ref, open, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_close(ged_polygon_ref ref)
{
    return ged_draw_obol_polygon_close(
	    ref, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_clear_selected_point(
    ged_polygon_ref ref)
{
    return ged_draw_obol_polygon_clear_selected_point(
	    ref, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_remove(ged_polygon_ref ref)
{
    return ged_draw_obol_polygon_remove(
	    ref, NULL);
}

extern "C" void *
ged_draw_obol_view_polygon_user_data(ged_polygon_ref ref)
{
    return ged_draw_obol_polygon_user_data(
	    ref, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_user_data_set(
    ged_polygon_ref ref,
    void *user_data)
{
    return ged_draw_obol_polygon_user_data_set(
	    ref, user_data, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_csg(
    ged_polygon_ref target,
    ged_polygon_ref stencil,
    bg_clip_t op)
{
    return ged_draw_obol_polygon_csg(
	    target,
	    stencil, op, NULL);
}

static int
ged_obol_selection_kind_from_ged(int kind)
{
    switch (kind) {
	case -1:
	    return BOBOL_SELECTION_ALL;
	case 0:
	    return BOBOL_SELECTION_SELECTED_PATH;
	case 1:
	    return BOBOL_SELECTION_HIGHLIGHTED_REF;
	default:
	    return INT_MIN;
    }
}

static int
ged_draw_obol_selection_available_impl(struct ged_view_context *view_ctx)
{
    return ged_obol_view_controller_ensure_for_context(view_ctx, 1) ? 1 : 0;
}

static size_t
ged_draw_obol_selection_count_impl(struct ged_view_context *view_ctx)
{
    BObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller)
	return 0;
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    return controller->selection().count(&owner, BOBOL_SELECTION_ALL);
}

static int
ged_draw_obol_selection_clear_impl(struct ged_view_context *view_ctx)
{
    BObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller)
	return 0;
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    controller->selection().clear(&owner, BOBOL_SELECTION_ALL);
    return 1;
}

extern "C" int
ged_draw_obol_view_context_selection_available(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_selection_available_impl(view_ctx);
}

extern "C" size_t
ged_draw_obol_view_context_selection_count(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_selection_count_impl(view_ctx);
}

struct ged_obol_selection_foreach_context {
    struct ged_view_context *view_ctx;
    ged_view_selection_visit_cb cb;
    void *data;
    int ok;
};

static int
ged_obol_selection_foreach_callback(const SbString &path, void *data)
{
    struct ged_obol_selection_foreach_context *ctx =
	    static_cast<struct ged_obol_selection_foreach_context *>(data);
    const char *path_str = path.getString();
    if (!ctx || !ctx->cb)
	return 0;
    if (!path_str || !path_str[0])
	return 1;
    if (!ctx->cb(ctx->view_ctx, path_str, ctx->data)) {
	ctx->ok = 0;
	return 0;
    }
    return 1;
}

extern "C" int
ged_draw_obol_view_context_selection_path_foreach(
    struct ged_view_context *view_ctx,
    ged_view_selection_visit_cb cb,
    void *data)
{
    if (!cb)
	return 0;
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    struct ged_obol_selection_foreach_context ctx;
    ctx.view_ctx = view_ctx;
    ctx.cb = cb;
    ctx.data = data;
    ctx.ok = 1;
    controller->selection().visitPaths(ged_obol_selection_foreach_callback,
				       &ctx, &owner, BOBOL_SELECTION_ALL);
    return ctx.ok;
}

extern "C" int
ged_draw_obol_view_context_selection_clear(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_selection_clear_impl(view_ctx);
}

extern "C" int
ged_draw_obol_view_context_selection_contains_path(
    struct ged_view_context *view_ctx,
    int kind,
    const char *path)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !path || !path[0])
	return 0;

    int obol_kind = ged_obol_selection_kind_from_ged(kind);
    if (obol_kind == INT_MIN)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    return controller->selection().containsPath(
	       ged_obol_skip_leading_slash(path), obol_kind, &owner) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_selection_add_path(
    struct ged_view_context *view_ctx,
    int kind,
    const char *path)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !path || !path[0])
	return 0;

    int obol_kind = ged_obol_selection_kind_from_ged(kind);
    if (obol_kind == INT_MIN || obol_kind == BOBOL_SELECTION_ALL)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    return controller->selection().addPath(ged_obol_skip_leading_slash(path),
					   obol_kind, &owner) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_selection_set_path(
    struct ged_view_context *view_ctx,
    int kind,
    const char *path)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;

    int obol_kind = ged_obol_selection_kind_from_ged(kind);
    if (obol_kind == INT_MIN)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    if (obol_kind == BOBOL_SELECTION_ALL) {
	controller->selection().clear(&owner, BOBOL_SELECTION_ALL);
	return path && path[0] ? 0 : 1;
    }

    return controller->selection().setPath(
	       ged_obol_skip_leading_slash(path), obol_kind, &owner) ? 1 : 0;
}


extern "C" int
ged_draw_obol_view_context_line_layer_builder_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const struct bg_line_layer_builder *builder)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;
    if (!builder || !bg_line_layer_builder_point_count(builder))
	return ged_obol_remove_feature(controller, view_ctx, name,
				       local ? 1 : 0) ? 1 : 1;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    BObolFeatureHandle handle =
	controller->features().publishLineLayerBuilder(name,
	    local ? BObolFeatureScope::Local :
	    BObolFeatureScope::Shared, builder, NULL,
	    local ? &owner : NULL);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_diagnostic_line_layer_builder_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct bg_line_layer_builder *builder)
{
    int ret = ged_draw_obol_view_context_line_layer_builder_replace(
		  view_ctx, name, 0, builder);
    if (!ret || !builder || !bg_line_layer_builder_point_count(builder))
	return ret;

    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, 0, 0);
    BObolFeatureHandle handle = controller ?
				  controller->features().find(name, BOBOL_FEATURE_SCOPE_SHARED) :
				  BObolFeatureHandle();
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BObolOverlayClass::Diagnostic,
						 BObolOverlayLifecycle::PerCommand,
						 BObolOverlayOrder::PostTransparent,
						 name));
}

extern "C" int
ged_draw_obol_view_context_line_layers_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const struct ged_annotation_line_layer *layers,
    size_t layer_count,
    const struct ged_view_feature_style *style)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;

    size_t live_layers = 0;
    for (size_t i = 0; layers && i < layer_count; i++) {
	if (layers[i].points && layers[i].point_count)
	    live_layers++;
    }
    if (!layers || !live_layers)
	return ged_obol_remove_feature(controller, view_ctx, name,
				       local ? 1 : 0) ? 1 : 1;

    std::vector<BObolLineLayer> obol_layers;
    obol_layers.reserve(live_layers);
    for (size_t i = 0; i < layer_count; i++) {
	if (!layers[i].points || !layers[i].point_count)
	    continue;

	BObolLineLayer layer;
	layer.name = layers[i].name ? layers[i].name : name;
	layer.style = ged_obol_feature_style_from_ged(&layers[i].style);
	layer.points.reserve(layers[i].point_count);
	layer.commands.reserve(layers[i].point_count);
	for (size_t j = 0; j < layers[i].point_count; j++) {
	    layer.points.push_back(SbVec3f(
				       static_cast<float>(layers[i].points[j][0]),
				       static_cast<float>(layers[i].points[j][1]),
				       static_cast<float>(layers[i].points[j][2])));
	    const int command = layers[i].commands ?
				layers[i].commands[j] : -1;
	    layer.commands.push_back(ged_obol_line_command_from_ged(command,
				     j));
	}
	obol_layers.push_back(layer);
    }

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    BObolFeatureHandle handle =
	controller->features().publishLineLayers(name,
	    local ? BObolFeatureScope::Local :
	    BObolFeatureScope::Shared, obol_layers,
	    style ? &obol_style : NULL, local ? &owner : NULL);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_lines_create_model_annotation(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const point_t point)
{
    if (!point)
	return 0;

    point_t points[1];
    VMOVE(points[0], point);
    int cmds[1] = {GED_DRAW_VIEW_LINE_MOVE};
    struct ged_view_feature_style style =
	    GED_VIEW_FEATURE_STYLE_INIT;
    style.color_valid = 1;
    style.color[0] = 255;
    style.color[1] = 0;
    style.color[2] = 0;
    int ret = ged_draw_obol_view_context_lines_replace(view_ctx, name, local,
	      points, cmds, 1, &style);
    if (!ret)
	return 0;

    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 0);
    BObolFeatureHandle handle = local ?
				  ged_obol_feature_handle(controller, view_ctx, name) :
				  (controller ? controller->features().find(name,
				      BOBOL_FEATURE_SCOPE_SHARED) : BObolFeatureHandle());
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BObolOverlayClass::UserAnnotation,
						 BObolOverlayLifecycle::Persistent,
						 BObolOverlayOrder::Model,
						 name));
}

extern "C" int
ged_draw_obol_view_context_lines_append_point(
    struct ged_view_context *view_ctx,
    const char *name,
    const point_t point)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    if (!lookup.handle.isValid() || !lookup.controller || !point)
	return 0;
    return lookup.controller->features().appendLinePoint(lookup.handle, SbVec3f(
		static_cast<float>(point[X]),
		static_cast<float>(point[Y]),
		static_cast<float>(point[Z]))) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_label_create(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const char *text,
    const point_t point,
    const point_t target,
    int has_target)
{
    if (!text || !point || (has_target && !target))
	return 0;

    struct ged_annotation_label label = GED_ANNOTATION_LABEL_INIT;
    label.text = text;
    VMOVE(label.point, point);
    if (has_target) {
	label.line_flag = 1;
	VMOVE(label.target, target);
    }

    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;

    BObolFeatureStyle style;
    style.hasColor = TRUE;
    style.color = SbColor(1.0f, 1.0f, 0.0f);

    BObolFeatureHandle handle =
	ged_obol_publish_labels(controller, view_ctx, name, local,
				ged_obol_labels_from_ged(&label, 1), &style);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_labels_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const struct ged_annotation_label *labels,
    size_t label_count)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;
    if (!labels || !label_count)
	return ged_obol_remove_feature(controller, view_ctx, name,
				       local ? 1 : 0) ? 1 : 1;

    BObolFeatureHandle handle =
	ged_obol_publish_labels(controller, view_ctx, name, local,
				ged_obol_labels_from_ged(labels, label_count));
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_tcl_labels_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    int draw,
    const struct ged_annotation_label *labels,
    size_t label_count)
{
    if (!draw || !labels || !label_count) {
	BObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	return (controller && name) ?
	       (ged_obol_remove_feature(controller, view_ctx, name, 1) ? 1 : 1) :
	       0;
    }
    return ged_draw_obol_view_context_labels_replace(view_ctx, name, 1,
	    labels, label_count);
}

extern "C" size_t
ged_draw_obol_view_context_label_count(struct ged_view_context *view_ctx, const char *name)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<BObolLabel> labels;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().labels(lookup.handle, labels))
	return 0;
    return labels.size();
}

extern "C" int
ged_draw_obol_view_context_label_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t index,
    struct bu_vls *text,
    point_t point,
    unsigned char rgb[3])
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<BObolLabel> labels;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().labels(lookup.handle, labels) ||
	index >= labels.size())
	return 0;

    const BObolLabel &label = labels[index];
    if (text) {
	bu_vls_trunc(text, 0);
	bu_vls_strcat(text, label.text.getString());
    }
    if (point)
	ged_obol_point_from_sb(point, label.point);
    if (rgb) {
	BObolFeatureStyle style;
	const SbColor color = label.hasColor ? label.color :
			      (lookup.controller->features().style(lookup.handle, style) &&
			       style.hasColor ?
			       style.color : SbColor(1.0f, 1.0f, 1.0f));
	ged_obol_rgb_from_color(color, rgb);
    }
    return 1;
}

extern "C" int
ged_draw_obol_view_context_label_point_set(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t index,
    const point_t point)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<BObolLabel> labels;
    if (!lookup.handle.isValid() || !lookup.controller || !point ||
	!lookup.controller->features().labels(lookup.handle, labels) ||
	index >= labels.size())
	return 0;
    labels[index].point = SbVec3f(
			      static_cast<float>(point[X]),
			      static_cast<float>(point[Y]),
			      static_cast<float>(point[Z]));
    return lookup.controller->features().replaceLabels(lookup.handle,
	    labels) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_line_style_get(
    struct ged_view_context *view_ctx,
    const char *name,
    struct ged_draw_view_line_style *style)
{
    if (!style)
	return 0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    BObolFeatureStyle obol_style;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().style(lookup.handle, obol_style))
	return 0;
    unsigned char rgb[3] = {255, 255, 255};
    if (obol_style.hasColor)
	ged_obol_rgb_from_color(obol_style.color, rgb);
    style->color[0] = rgb[0];
    style->color[1] = rgb[1];
    style->color[2] = rgb[2];
    style->line_width = obol_style.hasLineWidth ? obol_style.lineWidth : 1;
    return 1;
}

extern "C" int
ged_draw_obol_view_context_line_color_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int r,
    int g,
    int b)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    if (!lookup.handle.isValid() || !lookup.controller)
	return 0;
    unsigned char rgb[3] = {
	static_cast<unsigned char>(std::max(0, std::min(255, r))),
	static_cast<unsigned char>(std::max(0, std::min(255, g))),
	static_cast<unsigned char>(std::max(0, std::min(255, b)))
    };
    return lookup.controller->features().setColor(lookup.handle,
	    ged_obol_color_from_rgb(rgb)) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_line_width_set(
    struct ged_view_context *view_ctx,
    const char *name,
    int line_width)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    if (!lookup.handle.isValid() || !lookup.controller)
	return 0;
    return lookup.controller->features().setLineWidth(lookup.handle,
	    line_width) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_feature_points_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t **points,
    size_t *point_count)
{
    if (points)
	*points = NULL;
    if (point_count)
	*point_count = 0;
    if (!points || !point_count)
	return 0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<SbVec3f> obol_points;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().points(lookup.handle,
		obol_points))
	return 0;

    if (obol_points.empty())
	return 1;
    *points = (point_t *)bu_calloc(obol_points.size(), sizeof(point_t),
				   "GED Obol feature points copy");
    for (size_t i = 0; i < obol_points.size(); i++)
	ged_obol_point_from_sb((*points)[i], obol_points[i]);
    *point_count = obol_points.size();
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_indices_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    int **indices,
    size_t *index_count)
{
    if (indices)
	*indices = NULL;
    if (index_count)
	*index_count = 0;
    if (!indices || !index_count)
	return 0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<int32_t> obol_indices;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().indices(lookup.handle,
		obol_indices))
	return 0;

    if (obol_indices.empty())
	return 1;
    *indices = (int *)bu_calloc(obol_indices.size(), sizeof(int),
				"GED Obol feature indices copy");
    for (size_t i = 0; i < obol_indices.size(); i++)
	(*indices)[i] = static_cast<int>(obol_indices[i]);
    *index_count = obol_indices.size();
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_line_command_at(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t index,
    int *out)
{
    if (out)
	*out = 0;
    if (!out)
	return 0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    int32_t command = 0;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().lineCommandAt(lookup.handle,
		index, command))
	return 0;
    *out = ged_obol_line_command_to_ged(command);
    return 1;
}

static int
ged_obol_feature_layer_lookup(struct ged_view_context *view_ctx,
			      const char *name,
			      size_t layer_index,
			      BObolLineLayer &layer)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    BObolFeatureRecord record;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().record(lookup.handle,
		record) || record.kind != BObolFeatureKind::LineLayer ||
	layer_index >= record.layers.size())
	return 0;

    layer = record.layers[layer_index];
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_layer_points_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t layer_index,
    point_t **points,
    size_t *point_count)
{
    if (points)
	*points = NULL;
    if (point_count)
	*point_count = 0;
    if (!points || !point_count)
	return 0;

    BObolLineLayer layer;
    if (!ged_obol_feature_layer_lookup(view_ctx, name, layer_index, layer))
	return 0;

    if (layer.points.empty())
	return 1;
    *points = (point_t *)bu_calloc(layer.points.size(), sizeof(point_t),
				   "GED Obol feature layer points copy");
    for (size_t i = 0; i < layer.points.size(); i++)
	ged_obol_point_from_sb((*points)[i], layer.points[i]);
    *point_count = layer.points.size();
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_layer_line_command_at(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t layer_index,
    size_t point_index,
    int *out)
{
    if (out)
	*out = 0;
    if (!out)
	return 0;

    BObolLineLayer layer;
    if (!ged_obol_feature_layer_lookup(view_ctx, name, layer_index, layer) ||
	point_index >= layer.commands.size())
	return 0;

    *out = ged_obol_line_command_to_ged(layer.commands[point_index]);
    return 1;
}

extern "C" int
ged_draw_obol_view_context_tcl_lines_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_draw_view_line_style *style)
{
    if (point_count % 2)
	return 0;

    if (!points || point_count < 2) {
	BObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	return (controller && name) ?
	       (ged_obol_remove_feature(controller, view_ctx, name, 1) ? 1 : 1) :
	       0;
    }

    std::vector<int> cmds(point_count, GED_DRAW_VIEW_LINE_DRAW);
    for (size_t i = 0; i + 1 < point_count; i += 2) {
	cmds[i] = GED_DRAW_VIEW_LINE_MOVE;
	cmds[i + 1] = GED_DRAW_VIEW_LINE_DRAW;
    }

    struct ged_view_feature_style feature_style =
	    GED_VIEW_FEATURE_STYLE_INIT;
    const struct ged_view_feature_style *feature_stylep = NULL;
    if (style) {
	feature_style.visible = 1;
	feature_style.color_valid = 1;
	feature_style.color[0] = static_cast<unsigned char>(
				     std::max(0, std::min(255, style->color[0])));
	feature_style.color[1] = static_cast<unsigned char>(
				     std::max(0, std::min(255, style->color[1])));
	feature_style.color[2] = static_cast<unsigned char>(
				     std::max(0, std::min(255, style->color[2])));
	feature_style.line_width = style->line_width;
	feature_stylep = &feature_style;
    }

    int ret = ged_draw_obol_view_context_lines_replace(view_ctx, name, 1,
	      points, cmds.data(), point_count, feature_stylep);
    if (!ret)
	return 0;

    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BObolOverlayClass::TclOverlay,
						 BObolOverlayLifecycle::PerCommand,
						 BObolOverlayOrder::PostTransparent,
						 name));
}

extern "C" int
ged_draw_obol_view_context_arrow_tip_get(
    struct ged_view_context *view_ctx,
    const char *name,
    fastf_t *tip_length,
    fastf_t *tip_width)
{
    if (tip_length)
	*tip_length = 0.0;
    if (tip_width)
	*tip_width = 0.0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    float length = 0.0f;
    float width = 0.0f;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().arrowTip(lookup.handle,
		length, width))
	return 0;
    if (tip_length)
	*tip_length = length;
    if (tip_width)
	*tip_width = width;
    return 1;
}

extern "C" int
ged_draw_obol_view_context_arrow_tip_set(
    struct ged_view_context *view_ctx,
    const char *name,
    fastf_t tip_length,
    fastf_t tip_width)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    if (!lookup.handle.isValid() || !lookup.controller)
	return 0;
    return lookup.controller->features().setArrowTip(lookup.handle,
	    static_cast<float>(tip_length),
	    static_cast<float>(tip_width)) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_tcl_arrows_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t *points,
    size_t point_count,
    const struct ged_view_feature_style *style)
{
    if (point_count % 2)
	return 0;
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name)
	return 0;
    if (!points || !point_count)
	return ged_obol_remove_feature(controller, view_ctx, name, 1) ? 1 : 1;

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle =
	ged_obol_publish_arrow(controller, view_ctx, name, 1,
			       ged_obol_points_from_ged((const point_t *)points, point_count),
			       style ? &obol_style : NULL);
    if (!handle.isValid())
	return 0;
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BObolOverlayClass::TclOverlay,
						 BObolOverlayLifecycle::PerCommand,
						 BObolOverlayOrder::PostTransparent,
						 name));
}

extern "C" int
ged_draw_obol_view_context_feature_axes_centers_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t **centers,
    size_t *center_count)
{
    if (centers)
	*centers = NULL;
    if (center_count)
	*center_count = 0;
    if (!centers || !center_count)
	return 0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<SbVec3f> obol_centers;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().axesCenters(lookup.handle,
	    obol_centers))
	return 0;

    if (obol_centers.empty())
	return 1;
    *centers = (point_t *)bu_calloc(obol_centers.size(), sizeof(point_t),
				    "GED Obol feature axes centers copy");
    for (size_t i = 0; i < obol_centers.size(); i++)
	ged_obol_point_from_sb((*centers)[i], obol_centers[i]);
    *center_count = obol_centers.size();
    return 1;
}

extern "C" int
ged_draw_obol_view_context_tcl_axes_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t *centers,
    size_t center_count,
    fastf_t half_axes_size,
    const struct ged_view_feature_style *style)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !centers || !center_count)
	return 0;

    BObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BObolFeatureHandle handle =
	ged_obol_publish_axes(controller, view_ctx, name, 1,
			       ged_obol_points_from_ged((const point_t *)centers, center_count),
			      static_cast<float>(half_axes_size),
			      style ? &obol_style : NULL);
    if (!handle.isValid())
	return 0;
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BObolOverlayClass::TclOverlay,
						 BObolOverlayLifecycle::PerCommand,
						 BObolOverlayOrder::PostTransparent,
						 name));
}

static BObolFeatureStyle
ged_obol_axes_style_from_state(const struct ged_annotation_axes *state)
{
    BObolFeatureStyle style;
    if (!state)
	return style;
    style.hasColor = TRUE;
    unsigned char rgb[3] = {
	static_cast<unsigned char>(std::max(0, std::min(255,
	state->color[0]))),
	static_cast<unsigned char>(std::max(0, std::min(255,
	state->color[1]))),
	static_cast<unsigned char>(std::max(0, std::min(255,
	state->color[2])))
    };
    style.color = ged_obol_color_from_rgb(rgb);
    style.hasLineWidth = TRUE;
    style.lineWidth = state->line_width;
    return style;
}

extern "C" int
ged_draw_obol_view_context_axes_create(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    const struct ged_annotation_axes *state)
{
    if (!state)
	return 0;

    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;

    BObolFeatureStyle style = ged_obol_axes_style_from_state(state);
    point_t centers[1];
    VMOVE(centers[0], state->position);
    BObolFeatureHandle handle =
	ged_obol_publish_axes(controller, view_ctx, name, local,
			      ged_obol_points_from_ged(centers, 1),
			      static_cast<float>(state->size),
			      &style);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_axes_state_get(
    struct ged_view_context *view_ctx,
    const char *name,
    struct ged_annotation_axes *state)
{
    if (!state)
	return 0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<SbVec3f> centers;
    float half_axes_size = 0.0f;
    if (!lookup.controller || !lookup.handle.isValid() ||
	!lookup.controller->features().axesCenters(lookup.handle,
	    centers, &half_axes_size) || centers.empty())
	return 0;

    memset(state, 0, sizeof(*state));
    ged_obol_point_from_sb(state->position, centers[0]);
    state->size = half_axes_size;
    BObolFeatureStyle style;
    if (lookup.controller->features().style(lookup.handle, style)) {
	unsigned char rgb[3] = {255, 255, 255};
	if (style.hasColor)
	    ged_obol_rgb_from_color(style.color, rgb);
	state->color[0] = rgb[0];
	state->color[1] = rgb[1];
	state->color[2] = rgb[2];
	state->line_width = style.hasLineWidth ? style.lineWidth : 1;
    }
    return 1;
}

extern "C" int
ged_draw_obol_view_context_axes_state_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct ged_annotation_axes *state)
{
    if (!state)
	return 0;

    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    if (!handle.isValid())
	return ged_draw_obol_view_context_axes_create(view_ctx, name, 1, state);

    BObolFeatureStyle style = ged_obol_axes_style_from_state(state);
    point_t centers[1];
    VMOVE(centers[0], state->position);
    BObolFeatureHandle updated =
	ged_obol_publish_axes(controller, view_ctx, name, 1,
			      ged_obol_points_from_ged(centers, 1),
			      static_cast<float>(state->size),
			      &style);
    return updated.isValid() ? 1 : 0;
}

static std::string
ged_obol_overlay_group_path_for_name(const char *name)
{
    std::string path("_overlays");
    if (name && name[0]) {
	path += "/";
	path += name;
    }
    return path;
}

static std::string
ged_obol_overlay_shape_path_for_name(const char *name)
{
    std::string path = ged_obol_overlay_group_path_for_name(name);
    path += "/geometry";
    return path;
}

static void
ged_obol_overlay_group_state_set(BObolSceneController *scene,
				 const char *path,
				 const unsigned char *rgb,
				 fastf_t transparency,
				 int draw_mode)
{
    if (!scene || !path || !path[0])
	return;

    const SbColor color = rgb ? ged_obol_color_from_rgb(rgb) :
			  SbColor(1.0f, 1.0f, 1.0f);
    const int obol_draw_mode = draw_mode ?
			       ged_obol_lod_draw_mode_from_ged(draw_mode) :
			       BOBOL_LOD_DRAW_DIAGNOSTIC;
    (void)scene->setGroupDrawIntent(path, path, obol_draw_mode,
				    BOBOL_LOD_DRAW_WIRE, TRUE, 0);
    (void)scene->setGroupDisplayState(path, TRUE, FALSE, FALSE, 0, 0,
				      static_cast<float>(transparency), rgb ? TRUE : FALSE, color,
				      rgb ? TRUE : FALSE, color, 0);
}

template <typename ShapeT>
static void
ged_obol_overlay_shape_common_set(ShapeT *shape,
				  const char *name,
				  const std::string &shape_path,
				  const unsigned char basecolor[3],
				  fastf_t transparency,
				  int draw_mode,
				  const char *source_type,
				  const char *geometry_kind)
{
    if (!shape)
	return;

    const char *display_name = name && name[0] ? name : "overlay";
    shape->sourcePath = shape_path.c_str();
    shape->sourceName = display_name;
    shape->sourceType = source_type;
    shape->displayName = display_name;
    shape->geometryName = display_name;
    shape->sourceIdentity = shape_path.c_str();
    shape->cacheIdentity = shape_path.c_str();
    shape->ownerSourcePath = "_overlays";
    shape->databaseIntent = FALSE;
    shape->overlayIntent = TRUE;
    shape->hudIntent = FALSE;
    shape->localSource = TRUE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = draw_mode ?
		      ged_obol_lod_draw_mode_from_ged(draw_mode) :
		      BOBOL_LOD_DRAW_DIAGNOSTIC;
    shape->recordRole = "overlay";
    shape->geometryKind = geometry_kind;
    shape->visible = TRUE;
    shape->selectable = TRUE;
    shape->transparency = static_cast<float>(transparency);
    shape->colorOverride = TRUE;
    shape->color = ged_obol_color_from_rgb(basecolor);
    shape->materialColorValid = TRUE;
    shape->materialColor = ged_obol_color_from_rgb(basecolor);
}

static std::vector<SbVec3f>
ged_obol_overlay_points(const struct ged_draw_overlay_geometry *geometry)
{
    std::vector<SbVec3f> points;
    if (!geometry || !geometry->points || !geometry->point_count)
	return points;

    points.reserve(geometry->point_count);
    for (size_t i = 0; i < geometry->point_count; i++)
	points.push_back(SbVec3f(
			     static_cast<float>(geometry->points[i][X]),
			     static_cast<float>(geometry->points[i][Y]),
			     static_cast<float>(geometry->points[i][Z])));
    return points;
}

static SoBRLVListShape *
ged_obol_overlay_vlist_shape_create(const char *name,
				    const std::string &shape_path,
				    const struct ged_draw_overlay_geometry *geometry,
				    const unsigned char basecolor[3],
				    fastf_t transparency,
				    int draw_mode)
{
    if (!geometry)
	return NULL;

    std::vector<SbVec3f> points = ged_obol_overlay_points(geometry);
    if (points.empty())
	return NULL;

    std::vector<int32_t> commands;
    const char *source_type = "line-set";
    const char *geometry_kind = "line";
    if (geometry->kind == GED_DRAW_OVERLAY_GEOMETRY_POINT_SET) {
	commands.assign(points.size(),
			static_cast<int32_t>(BObolLineCommand::Point));
	source_type = "point-set";
	geometry_kind = "point";
    } else {
	if (geometry->kind != GED_DRAW_OVERLAY_GEOMETRY_LINE_SET ||
	    !geometry->commands ||
	    geometry->command_count != geometry->point_count)
	    return NULL;
	commands.reserve(geometry->command_count);
	for (size_t i = 0; i < geometry->command_count; i++)
	    commands.push_back(ged_obol_line_command_from_ged(
				   geometry->commands[i], i));
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    ged_obol_overlay_shape_common_set(shape, name, shape_path, basecolor,
				      transparency, draw_mode, source_type, geometry_kind);
    shape->setLineSet(points.data(), commands.data(),
		      static_cast<int>(points.size()));
    return shape;
}

static void
ged_obol_overlay_append_face_triangles(const std::vector<int32_t> &face,
				       size_t point_count,
				       std::vector<int32_t> &triangles)
{
    if (face.size() < 3)
	return;

    std::vector<int32_t> clean_face = face;
    if (clean_face.size() > 3 && clean_face.front() == clean_face.back())
	clean_face.pop_back();
    if (clean_face.size() < 3)
	return;

    for (size_t i = 0; i < clean_face.size(); i++) {
	if (clean_face[i] < 0 ||
	    static_cast<size_t>(clean_face[i]) >= point_count)
	    return;
    }

    for (size_t i = 1; i + 1 < clean_face.size(); i++) {
	triangles.push_back(clean_face[0]);
	triangles.push_back(clean_face[i]);
	triangles.push_back(clean_face[i + 1]);
    }
}

static std::vector<int32_t>
ged_obol_overlay_triangles(const struct ged_draw_overlay_geometry *geometry)
{
    std::vector<int32_t> triangles;
    if (!geometry || !geometry->indices || !geometry->index_count ||
	!geometry->point_count)
	return triangles;

    std::vector<int32_t> face;
    for (size_t i = 0; i < geometry->index_count; i++) {
	const int32_t idx = static_cast<int32_t>(geometry->indices[i]);
	if (idx < 0) {
	    ged_obol_overlay_append_face_triangles(face,
						   geometry->point_count, triangles);
	    face.clear();
	} else {
	    face.push_back(idx);
	}
    }
    ged_obol_overlay_append_face_triangles(face, geometry->point_count,
					   triangles);
    return triangles;
}

static SoBRLMeshShape *
ged_obol_overlay_mesh_shape_create(const char *name,
				   const std::string &shape_path,
				   const struct ged_draw_overlay_geometry *geometry,
				   const unsigned char basecolor[3],
				   fastf_t transparency,
				   int draw_mode)
{
    if (!geometry || geometry->kind != GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET)
	return NULL;

    std::vector<SbVec3f> points = ged_obol_overlay_points(geometry);
    std::vector<int32_t> triangles = ged_obol_overlay_triangles(geometry);
    if (points.empty() || triangles.empty())
	return NULL;

    SoBRLMeshShape *shape = new SoBRLMeshShape;
    ged_obol_overlay_shape_common_set(shape, name, shape_path, basecolor,
				      transparency, draw_mode, "indexed-face-set", "surface");
    shape->setIndexedTriangles(points.data(), static_cast<int>(points.size()),
			       triangles.data(), static_cast<int>(triangles.size()));
    return shape;
}

extern "C" int
ged_draw_obol_overlay_erase_name_context(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *name)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !name[0])
	return 0;

    BObolSceneController *scene = controller->getSceneController();
    if (!scene)
	return 0;

    const std::string overlay_path = ged_obol_overlay_group_path_for_name(name);
    const int removed = scene->removeGroup(overlay_path.c_str());
    if (scene->getGroupChildCount("_overlays") == 0)
	(void)scene->removeGroup("_overlays");
    if (removed > 0) {
	scene->realizePending();
	controller->requestRender("overlay-remove");
    }
    (void)gedp;
    return removed > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_overlay_geometry_insert_context(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *name,
    const struct ged_draw_overlay_geometry *geometry,
    const unsigned char basecolor[3],
    fastf_t transparency,
    int draw_mode,
    char **shape_path_out)
{
    if (shape_path_out)
	*shape_path_out = NULL;
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !name[0] || !geometry || !basecolor)
	return 0;

    BObolSceneController *scene = controller->getSceneController();
    if (!scene)
	return 0;

    const std::string overlay_path = ged_obol_overlay_group_path_for_name(name);
    const std::string shape_path = ged_obol_overlay_shape_path_for_name(name);
    (void)scene->removeGroup(overlay_path.c_str());
    if (!scene->ensureGroup("_overlays") ||
	!scene->ensureGroup(overlay_path.c_str()))
	return 0;

    ged_obol_overlay_group_state_set(scene, "_overlays", NULL, 0.0, 0);
    ged_obol_overlay_group_state_set(scene, overlay_path.c_str(), basecolor,
				     transparency, draw_mode);

    SoNode *shape = NULL;
    if (geometry->kind == GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET)
	shape = ged_obol_overlay_mesh_shape_create(name, shape_path, geometry,
		basecolor, transparency, draw_mode);
    else
	shape = ged_obol_overlay_vlist_shape_create(name, shape_path, geometry,
		basecolor, transparency, draw_mode);
    if (!shape) {
	(void)scene->removeGroup(overlay_path.c_str());
	if (scene->getGroupChildCount("_overlays") == 0)
	    (void)scene->removeGroup("_overlays");
	return 0;
    }

    if (scene->appendChildToGroup(overlay_path.c_str(), shape) < 0) {
	shape->ref();
	shape->unref();
	(void)scene->removeGroup(overlay_path.c_str());
	if (scene->getGroupChildCount("_overlays") == 0)
	    (void)scene->removeGroup("_overlays");
	return 0;
    }

    if (shape_path_out)
	*shape_path_out = bu_strdup(shape_path.c_str());
    controller->requestRender("overlay-insert");
    (void)gedp;
    return 1;
}

static SbColor
ged_obol_summary_material_color(
    const BObolDatabaseSourceSummary &summary)
{
    if (summary.materialColorValid)
	return summary.materialColor;
    if (summary.databaseMaterialColorValid)
	return summary.databaseMaterialColor;
    if (summary.colorOverride)
	return summary.color;
    return SbColor(1.0f, 1.0f, 1.0f);
}

static SoBRLDatabaseSource *
ged_obol_owned_database_source_for_path(struct ged *gedp, const char *path)
{
    if (!gedp || !path || !path[0])
	return NULL;
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return NULL;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (!instance_keys.empty())
	return scene->findDatabaseSourceInstance(instance_keys.front().c_str());

    return scene->findDatabaseSource(path);
}

bool
ged_bobol_publication_begin(
    struct ged_bobol_publication_context *publication,
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    int draw_mode)
{
    if (!publication || publication->active || !gedp)
	return false;

    BObolSceneController *scene = ged_bobol_scene(gedp);
    if (!scene)
	return false;

    publication->gedp = gedp;
    publication->view_ctx = view_ctx;
    publication->scene = scene;
    publication->draw_mode = draw_mode;
    publication->group_path.clear();
    publication->appearance = NULL;
    publication->active = true;
    scene->beginSceneMutationBatch();
    return true;
}

void
ged_bobol_publication_appearance_set(
    struct ged_bobol_publication_context *publication,
    const struct ged_draw_appearance_settings *settings)
{
    if (publication && publication->active)
	publication->appearance = settings;
}

void
ged_bobol_publication_group_set(
    struct ged_bobol_publication_context *publication,
    const char *group_path)
{
    if (!publication || !publication->active)
	return;
    publication->group_path = group_path ?
	ged_obol_group_path_from_record_path(group_path) : std::string();
}

void
ged_bobol_publication_end(
    struct ged_bobol_publication_context *publication)
{
    if (!publication || !publication->active)
	return;

    BObolSceneController *scene = publication->scene;
    publication->active = false;
    publication->gedp = NULL;
    publication->view_ctx = NULL;
    publication->scene = NULL;
    publication->draw_mode = -1;
    publication->group_path.clear();
    publication->appearance = NULL;
    if (scene)
	scene->endSceneMutationBatch();
}

static int
ged_obol_database_source_instance_key_for_source(
    SoBRLDatabaseSource *source,
    std::string &source_instance_key)
{
    source_instance_key.clear();
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    const char *instance_key = summary.instanceKey.getString();
    if (instance_key && instance_key[0]) {
	source_instance_key = instance_key;
	return 1;
    }

    const char *source_path = summary.path.getString();
    if (source_path && source_path[0]) {
	source_instance_key = source_path;
	return 1;
    }

    return 0;
}

static int
ged_obol_database_source_scene_instance_for_path(
    struct ged *gedp,
    const char *path,
    BObolSceneController **scene_out,
    std::string &source_instance_key,
    const struct ged_bobol_publication_context *publication = NULL)
{
    source_instance_key.clear();
    if (scene_out)
	*scene_out = NULL;

    BObolSceneController *scene =
	(publication && publication->active) ? publication->scene :
	ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int publication_mode =
	(publication && publication->active) ? publication->draw_mode : -1;
    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path,
		publication_mode, 0, publication);
    SoBRLDatabaseSource *source = !instance_keys.empty() ?
	scene->findDatabaseSourceInstance(instance_keys.front().c_str()) :
	NULL;
    if (!source)
	return 0;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key) && source_instance_key.empty())
	return 0;

    if (scene_out)
	*scene_out = scene;
    return 1;
}

static int
ged_obol_database_source_mark_published_current(BObolSceneController *scene,
	SoBRLDatabaseSource *source)
{
    if (!scene || !source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    const char *instance_key = summary.instanceKey.getString();
    if (!instance_key || !instance_key[0])
	instance_key = source->instanceKey.getValue().getString();
    if (!instance_key || !instance_key[0])
	instance_key = summary.path.getString();
    if (!instance_key || !instance_key[0])
	return 0;

    if (scene->setDatabaseSourceInstanceRealizationRoleFlags(instance_key,
	    summary.realizationRoleFlags |
	    SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL) < 0)
	return 0;

    const int changed = scene->setDatabaseSourceInstanceRealizationState(
			    instance_key,
			    SoBRLDatabaseSource::REALIZED,
			    summary.sourceRevision,
			    summary.inputsRevision,
			    SoBRLDatabaseSource::STALE_NONE,
			    NULL);
    return changed >= 0 ? 1 : 0;
}

static bool
ged_obol_database_source_controller_summary_for_path(
    BObolSceneController *scene,
    const char *path,
    BObolDatabaseSourceSummary &summary)
{
    summary = BObolDatabaseSourceSummary();
    if (!scene || !path || !path[0])
	return false;

    return scene->getDatabaseSourceSummaryForPath(path, summary) &&
	   summary.valid;
}

static bool
ged_obol_database_source_controller_summary_for_source(
    BObolSceneController *scene,
    SoBRLDatabaseSource *source,
    BObolDatabaseSourceSummary &summary)
{
    summary = BObolDatabaseSourceSummary();
    if (!source)
	return false;

    BObolDatabaseSourceSummary source_summary;
    if (!source->getSummary(source_summary) || !source_summary.valid)
	return false;

    const char *source_instance_key =
	source_summary.instanceKey.getString();
    if (scene && source_instance_key && source_instance_key[0]) {
	if (scene->getDatabaseSourceSummaryForInstance(source_instance_key,
		summary) && summary.valid)
	    return true;
    }

    summary = source_summary;
    return true;
}

static ged_draw_stale_reason
ged_obol_database_source_stale_reason(
    const BObolDatabaseSourceSummary &summary)
{
    if (!summary.stale)
	return GED_DRAW_STALE_NONE;

    if (summary.realizationStatus == SoBRLDatabaseSource::FAILED)
	return GED_DRAW_STALE_UPDATE_FAILED;

    if (summary.staleReason &
	(SoBRLDatabaseSource::STALE_SOURCE |
	 SoBRLDatabaseSource::STALE_DATABASE))
	return GED_DRAW_STALE_SOURCE_CHANGED;

    if (summary.staleReason &
	(SoBRLDatabaseSource::STALE_INPUTS |
	 SoBRLDatabaseSource::STALE_VIEW))
	return GED_DRAW_STALE_VIEW_INPUT_CHANGED;

    if (summary.staleReason &
	(SoBRLDatabaseSource::STALE_DRAW |
	 SoBRLDatabaseSource::STALE_TESSELLATION))
	return GED_DRAW_STALE_SETTINGS_CHANGED;

    return GED_DRAW_STALE_SOURCE_CHANGED;
}

static uint32_t
ged_obol_stale_reason_from_ged(ged_draw_stale_reason reason)
{
    switch (reason) {
	case GED_DRAW_STALE_VIEW_INPUT_CHANGED:
	    return SoBRLDatabaseSource::STALE_INPUTS;
	case GED_DRAW_STALE_SETTINGS_CHANGED:
	    return SoBRLDatabaseSource::STALE_DRAW;
	case GED_DRAW_STALE_FORCED:
	    return SoBRLDatabaseSource::STALE_DRAW |
		   SoBRLDatabaseSource::STALE_TESSELLATION;
	case GED_DRAW_STALE_UPDATE_FAILED:
	    return SoBRLDatabaseSource::STALE_SOURCE;
	case GED_DRAW_STALE_SOURCE_CHANGED:
	case GED_DRAW_STALE_NONE:
	default:
	    return SoBRLDatabaseSource::STALE_SOURCE;
    }
}

static void
ged_obol_bounds_to_vmath(const SbBox3f &bounds, vect_t *min, vect_t *max);

extern "C" int
ged_draw_obol_scene_database_bounds(
    struct ged *gedp,
    vect_t *min,
    vect_t *max,
    int *empty_out)
{
    if (!gedp || !min || !max || !empty_out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    VSETALL(*min, INFINITY);
    VSETALL(*max, -INFINITY);
    *empty_out = 1;

    SbBox3f bounds;
    if (scene->getDatabaseSourceBounds(bounds, FALSE) && !bounds.isEmpty()) {
	ged_obol_bounds_to_vmath(bounds, min, max);
	*empty_out = 0;
    }

    return 1;
}

static int
ged_obol_database_record_all_halfspaces(
    struct ged *gedp,
    const struct ged_db_index_record &record,
    std::unordered_set<ged_db_index_id> &visiting)
{
    if (!gedp || !record.valid || !record.dp)
	return 0;
    if (record.dp->d_flags & RT_DIR_SOLID)
	return record.dp->d_minor_type == DB5_MINORTYPE_BRLCAD_HALF;
    if (!(record.dp->d_flags & RT_DIR_COMB))
	return 0;

    const ged_db_index_id object_id = record.object_id ?
	record.object_id : record.id;
    if (!object_id || !visiting.insert(object_id).second)
	return 0;

    const size_t child_count = ged_db_index_child_count(gedp, object_id);
    int all_halfspaces = child_count ? 1 : 0;
    for (size_t i = 0; all_halfspaces && i < child_count; i++) {
	struct ged_db_index_child child;
	memset(&child, 0, sizeof(child));
	if (!ged_db_index_child_at(gedp, object_id, i, &child) ||
	    !ged_obol_database_record_all_halfspaces(gedp, child.record,
		visiting))
	    all_halfspaces = 0;
    }
    visiting.erase(object_id);
    return all_halfspaces;
}

static int
ged_obol_database_source_member_autoview_bounds(
    struct ged *gedp,
    const char *path,
    vect_t *min,
    vect_t *max)
{
    if (!gedp || !gedp->dbip || !path || !path[0] || !min || !max)
	return 0;

    const char *normalized = ged_obol_skip_leading_slash(path);
    const size_t path_length =
	ged_db_index_path_resolve(gedp, normalized, NULL, 0);
    if (!path_length)
	return 0;

    std::vector<ged_db_index_id> ids(path_length);
    if (ged_db_index_path_resolve(gedp, normalized, ids.data(), ids.size()) !=
	path_length)
	return 0;

    const size_t child_count = ged_db_index_child_count(gedp, ids.back());
    if (!child_count)
	return 0;

    VSETALL(*min, INFINITY);
    VSETALL(*max, -INFINITY);
    int have_bounds = 0;
    std::unordered_set<ged_db_index_id> halfspace_visiting;
    for (size_t i = 0; i < child_count; i++) {
	struct ged_db_index_child child;
	memset(&child, 0, sizeof(child));
	if (!ged_db_index_child_at(gedp, ids.back(), i, &child) ||
	    !child.record.valid || child.bool_op == DB_OP_SUBTRACT)
	    continue;
	if (ged_obol_database_record_all_halfspaces(gedp, child.record,
		halfspace_visiting))
	    continue;

	const char *child_name = child.record.name;
	if ((!child_name || !child_name[0]) && child.record.dp)
	    child_name = child.record.dp->d_namep;
	if (!child_name || !child_name[0])
	    continue;

	struct bu_vls child_path = BU_VLS_INIT_ZERO;
	bu_vls_printf(&child_path, "%s/%s", normalized, child_name);
	const char *bounds_path = bu_vls_cstr(&child_path);
	point_t child_min;
	point_t child_max;
	struct bu_vls msgs = BU_VLS_INIT_ZERO;
	const int bounds_ret = rt_obj_bounds(&msgs, gedp->dbip, 1,
	    &bounds_path, 0, child_min, child_max);
	bu_vls_free(&msgs);
	bu_vls_free(&child_path);
	if (bounds_ret != BRLCAD_OK ||
	    !isfinite(child_min[X]) || !isfinite(child_min[Y]) ||
	    !isfinite(child_min[Z]) || !isfinite(child_max[X]) ||
	    !isfinite(child_max[Y]) || !isfinite(child_max[Z]))
	    continue;

	point_t center;
	VADD2SCALE(center, child_min, child_max, 0.5);
	fastf_t extent = child_max[X] - child_min[X];
	V_MAX(extent, child_max[Y] - child_min[Y]);
	V_MAX(extent, child_max[Z] - child_min[Z]);
	if (extent < SQRT_SMALL_FASTF)
	    extent = SQRT_SMALL_FASTF;

	/* Frame each finite, positive member around its own center before
	 * combining the results.  This keeps an infinite halfspace or tiny
	 * helper object from shifting the center of a much larger visible model,
	 * while retaining enough margin for rotation. */
	point_t padded_min;
	point_t padded_max;
	VSET(padded_min, center[X] - extent, center[Y] - extent,
	    center[Z] - extent);
	VSET(padded_max, center[X] + extent, center[Y] + extent,
	    center[Z] + extent);
	VMIN(*min, padded_min);
	VMAX(*max, padded_max);
	have_bounds = 1;
    }

    return have_bounds;
}

extern "C" int
ged_draw_obol_scene_database_autoview_bounds(
    struct ged *gedp,
    vect_t *min,
    vect_t *max,
    int *empty_out,
    int allow_member_bounds)
{
    if (!gedp || !min || !max || !empty_out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    VSETALL(*min, INFINITY);
    VSETALL(*max, -INFINITY);
    *empty_out = 1;

    /* Use database members when a displayed source is a combination.  The
     * retained renderer may flatten that combination into one aggregate
     * shape, whose raw AABB can be dominated by non-finite or tiny helper
     * geometry rather than the model the user expects to frame. */
    const int source_count = scene->getDatabaseSourceCount();
    int have_source_bounds = 0;
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary source;
	if (!scene->getDatabaseSourceSummary(i, source) || !source.valid ||
	    !source.visible)
	    continue;

	/* Prefer the source's own already-computed AABB (from the coarse proxy
	 * occurrence union or realized geometry): it is essentially free.  Only
	 * fall back to the expensive per-child rt_obj_bounds member walk when the
	 * cheap bounds are unavailable or non-finite (e.g. a halfspace-dominated
	 * aggregate the member walk deliberately skips).  Framing the view must
	 * never prep the whole assembly synchronously -- that was ~11s of the
	 * first-draw stall on large models. */
	const SbBool cheap_bounds_usable = source.sourceBoundsValid &&
	    !source.sourceBounds.isEmpty() &&
	    isfinite(source.sourceBounds.getMin()[0]) &&
	    isfinite(source.sourceBounds.getMin()[1]) &&
	    isfinite(source.sourceBounds.getMin()[2]) &&
	    isfinite(source.sourceBounds.getMax()[0]) &&
	    isfinite(source.sourceBounds.getMax()[1]) &&
	    isfinite(source.sourceBounds.getMax()[2]);

	if (!cheap_bounds_usable) {
	    /* The only remaining bound is the expensive per-child rt_obj_bounds
	     * walk.  Skip it on the coarse-first / progressive path (caller passes
	     * allow_member_bounds=0): framing there is refined progressively as
	     * cheap proxy bounds stream in, and prepping the whole assembly just to
	     * frame the first frame was ~11s of stall.  Explicit autoview commands
	     * pass 1 for precise framing. */
	    if (!allow_member_bounds)
		continue;
	    vect_t member_min;
	    vect_t member_max;
	    if (source.path.getLength() > 0 &&
		ged_obol_database_source_member_autoview_bounds(gedp,
		    source.path.getString(), &member_min, &member_max)) {
		VMIN(*min, member_min);
		VMAX(*max, member_max);
		have_source_bounds = 1;
	    }
	    continue;
	}

	const SbVec3f source_min = source.sourceBounds.getMin();
	const SbVec3f source_max = source.sourceBounds.getMax();
	const fastf_t center_x = (source_min[X] + source_max[X]) * 0.5;
	const fastf_t center_y = (source_min[Y] + source_max[Y]) * 0.5;
	const fastf_t center_z = (source_min[Z] + source_max[Z]) * 0.5;
	fastf_t extent = source_max[X] - source_min[X];
	extent = std::max(extent,
		static_cast<fastf_t>(source_max[Y] - source_min[Y]));
	extent = std::max(extent,
		static_cast<fastf_t>(source_max[Z] - source_min[Z]));
	if (extent < SQRT_SMALL_FASTF)
	    extent = SQRT_SMALL_FASTF;

	vect_t source_bounds_min;
	vect_t source_bounds_max;
	VSET(source_bounds_min, center_x - extent, center_y - extent,
	     center_z - extent);
	VSET(source_bounds_max, center_x + extent, center_y + extent,
	     center_z + extent);
	VMIN(*min, source_bounds_min);
	VMAX(*max, source_bounds_max);
	have_source_bounds = 1;
    }
    if (have_source_bounds) {
	*empty_out = 0;
	return 1;
    }

    SbBox3f bounds;
    if (scene->getDatabaseSourceBounds(bounds, TRUE) && !bounds.isEmpty()) {
	ged_obol_bounds_to_vmath(bounds, min, max);
	*empty_out = 0;
    }

    return 1;
}

static void
ged_obol_bounds_to_vmath(const SbBox3f &bounds, vect_t *min, vect_t *max)
{
    VSETALL(*min, INFINITY);
    VSETALL(*max, -INFINITY);
    if (bounds.isEmpty())
	return;

    const SbVec3f &bmin = bounds.getMin();
    const SbVec3f &bmax = bounds.getMax();
    VMOVE(*min, bmin);
    VMOVE(*max, bmax);
}

static int
ged_draw_obol_scene_subtree_bounds_for_path(
    struct ged *gedp,
    const char *path,
    vect_t *min,
    vect_t *max,
    int include_overlays,
    int *empty_out)
{
    if (empty_out)
	*empty_out = 1;
    if (min)
	VSETALL(*min, INFINITY);
    if (max)
	VSETALL(*max, -INFINITY);
    if (!gedp || !path || !path[0] || !min || !max || !empty_out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SbBox3f bounds;
    if (!scene->getSceneSubtreeBounds(path, include_overlays ? TRUE : FALSE,
				      bounds))
	return 0;

    ged_obol_bounds_to_vmath(bounds, min, max);
    *empty_out = bounds.isEmpty() ? 1 : 0;
    return 1;
}

extern "C" int
ged_draw_obol_database_source_bounds_for_path(
    struct ged *gedp,
    const char *path,
    vect_t *min,
    vect_t *max,
    int include_overlays,
    int *empty_out)
{
    return ged_draw_obol_scene_subtree_bounds_for_path(gedp, path, min, max,
	    include_overlays, empty_out);
}

extern "C" int
ged_draw_obol_group_subtree_bounds_for_path(
    struct ged *gedp,
    const char *path,
    vect_t *min,
    vect_t *max,
    int include_overlays,
    int *empty_out)
{
    if (!path || !path[0])
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    if (group_path.empty() && !BU_STR_EQUAL(path, "/"))
	return 0;
    return ged_draw_obol_scene_subtree_bounds_for_path(gedp,
	    group_path.empty() ? "/" : group_path.c_str(), min, max,
	    include_overlays, empty_out);
}

extern "C" int
ged_draw_obol_scene_root_child_count(struct ged *gedp, size_t *out)
{
    if (out)
	*out = 0;
    if (!out || !gedp)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int child_count = scene->getGroupChildCount("/");
    if (child_count < 0)
	return 0;

    *out = static_cast<size_t>(child_count);
    return 1;
}

static std::string
ged_obol_context_leaf_name_from_path(const char *path)
{
    if (!path || !path[0])
	return std::string();
    if (BU_STR_EQUAL(path, "/"))
	return std::string("/");

    const char *normalized = ged_obol_skip_leading_slash(path);
    const char *slash = strrchr(normalized, '/');
    return std::string((slash && slash[1]) ? slash + 1 : normalized);
}

extern "C" void
ged_draw_obol_scene_context_info_free(
    struct ged_draw_obol_scene_context_info *info)
{
    if (!info)
	return;
    if (info->path)
	bu_free(info->path, "GED Obol scene context path");
    if (info->instance_key)
	bu_free(info->instance_key, "GED Obol scene context instance key");
    if (info->name)
	bu_free(info->name, "GED Obol scene context name");
    memset(info, 0, sizeof(*info));
}

static int
ged_obol_scene_context_info_from_summary(
    const BObolSceneTreeSummary &summary,
    struct ged_draw_obol_scene_context_info *out)
{
    if (!out || !summary.valid || summary.path.getLength() == 0)
	return 0;

    ged_draw_obol_scene_context_info_free(out);
    const char *path = summary.path.getString();
    const std::string name = summary.displayName.getLength() > 0 ?
			     std::string(summary.displayName.getString()) :
			     ged_obol_context_leaf_name_from_path(path);

    out->path = bu_strdup(path);
    out->name = bu_strdup(name.c_str());
    out->node_kind = summary.nodeKind;
    out->is_group = summary.isGroup ? 1 : 0;
    out->is_database_source =
	(summary.isDatabaseSource ||
	 summary.nodeKind == BObolSceneTreeSummary::NODE_DATABASE_SOURCE) ?
	1 : 0;
    out->is_shape =
	(!out->is_database_source &&
	 (summary.isShape ||
	  summary.nodeKind == BObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	  summary.nodeKind == BObolSceneTreeSummary::NODE_MESH_SHAPE)) ?
	1 : 0;
    out->has_parent = summary.hasParent ? 1 : 0;
    out->draw_tree_depth = summary.drawTreeDepth;
    out->child_count = summary.childCount > 0 ?
		       static_cast<size_t>(summary.childCount) : 0;
    return 1;
}

extern "C" int
ged_draw_obol_scene_context_info_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_obol_scene_context_info *out)
{
    if (!out)
	return 0;
    memset(out, 0, sizeof(*out));
    if (!gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    BObolSceneTreeSummary summary;
    if (!scene->getSceneTreeSummaryForPath(path, summary))
	return 0;
    return ged_obol_scene_context_info_from_summary(summary, out);
}

extern "C" int
ged_draw_obol_scene_child_context_info_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    struct ged_draw_obol_scene_context_info *out)
{
    if (!out)
	return 0;
    memset(out, 0, sizeof(*out));
    if (!gedp || !path || !path[0] ||
	index > static_cast<size_t>(INT_MAX) ||
	!ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    BObolSceneTreeSummary summary;
    if (!scene->getSceneChildTreeSummary(path, static_cast<int>(index),
					 summary))
	return 0;
    return ged_obol_scene_context_info_from_summary(summary, out);
}

extern "C" int
ged_draw_obol_database_source_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_database_source_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    SoBRLDatabaseSource *source =
	scene ? scene->findDatabaseSourceInstance(source_instance_key.c_str()) :
	NULL;
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    out->valid = 1;
    out->is_database_source = 1;
    out->has_state = 1;
    out->stale = summary.stale ? 1 : 0;
    out->database_path = source->path.getValue().getString();
    out->source_revision = summary.sourceRevision;
    out->inputs_revision = summary.inputsRevision;
    out->realized_source_revision = summary.realizedSourceRevision;
    out->realized_inputs_revision = summary.realizedInputsRevision;

    const char *identity = summary.realizationIdentity.getString();
    if (identity && identity[0])
	out->realization_identity = bu_data_hash(identity,
				    strlen(identity) * sizeof(char));

    ged_draw_stale_reason reason =
	ged_obol_database_source_stale_reason(summary);
    out->stale_reason = ged_draw_source_stale_reason_name(reason);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_mark_stale_for_path(
    struct ged *gedp,
    const char *path,
    ged_draw_stale_reason reason)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int changed = scene->markDatabaseSourceInstanceStale(
				source_instance_key.c_str(),
				ged_obol_stale_reason_from_ged(reason));
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

static int
ged_obol_database_source_record_from_summary(
    struct ged_draw_obol_database_source_record *out,
    const BObolDatabaseSourceSummary &summary)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!summary.valid)
	return 0;

    out->valid = 1;
    out->database_path = summary.path.getString();
    out->source_path = summary.path.getString();
    out->instance_key = summary.instanceKey.getString();
    const char *owner_group_path = summary.parentGroupPath.getString();
    out->owner_group_path = (owner_group_path && owner_group_path[0]) ?
			    owner_group_path : summary.path.getString();
    out->visible = summary.visible ? 1 : 0;
    out->selected = summary.selected ? 1 : 0;
    out->highlighted = summary.highlighted ? 1 : 0;
    out->transparency = ged_obol_reported_transparency(summary.transparency);
    out->draw_mode = summary.representationMode >= 0 ?
		     summary.representationMode :
		     ged_obol_database_draw_mode_to_ged(summary.drawMode);
    out->material_policy =
	ged_obol_material_policy_to_ged(summary.materialPolicy);
    out->source_revision = summary.sourceRevision;
    out->inputs_revision = summary.inputsRevision;
    out->realized_source_revision = summary.realizedSourceRevision;
    out->realized_inputs_revision = summary.realizedInputsRevision;
    out->stale_reason = ged_obol_database_source_stale_reason(summary);

    const char *identity = summary.realizationIdentity.getString();
    if (identity && identity[0])
	out->realization_identity = bu_data_hash(identity,
				    strlen(identity) * sizeof(char));

    out->realization_status =
	summary.realizationStatus == SoBRLDatabaseSource::REALIZED ?
	GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_CURRENT :
	GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_STALE;
    out->realization_role_flags = summary.realizationRoleFlags;
    out->realization_view_dependent =
	summary.realizationViewDependent ? 1 : 0;
    out->realization_csg_lod_enabled =
	summary.realizationCsgLodEnabled ? 1 : 0;
    out->realization_mesh_lod_enabled =
	summary.realizationMeshLodEnabled ? 1 : 0;
    out->realization_view_scale = (fastf_t)summary.realizationViewScale;
    out->realization_lod_scale = (fastf_t)summary.realizationLodScale;
    out->realization_bot_threshold =
	(uint64_t)summary.realizationBotThreshold;
    out->realization_curve_scale = (fastf_t)summary.realizationCurveScale;
    out->realization_point_scale = (fastf_t)summary.realizationPointScale;
    return 1;
}

extern "C" int
ged_draw_obol_database_source_record_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_obol_database_source_record *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    if (!ged_obol_database_source_record_from_summary(out, summary))
	return 0;

    /* The record borrows strings.  Point it at source-owned fields rather
     * than the temporary summary's SbString storage. */
    out->database_path = source->path.getValue().getString();
    out->source_path = out->database_path;
    out->instance_key = source->instanceKey.getValue().getString();
    out->owner_group_path = out->database_path;
    return 1;
}

extern "C" int
ged_draw_obol_database_source_runtime_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_obol_database_source_runtime *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    out->valid = 1;
    out->database_path = source->path.getValue().getString();
    out->dbip = source->getDatabase();
    out->tessellation_abs_tol =
	(fastf_t)source->tessellationAbsTol.getValue();
    out->tessellation_rel_tol =
	(fastf_t)source->tessellationRelTol.getValue();
    out->tessellation_norm_tol =
	(fastf_t)source->tessellationNormTol.getValue();
    out->lod_bot_threshold =
	(uint64_t)source->lodBotThreshold.getValue();
    BObolDatabaseSourceSummary summary;
    if (source->getSummary(summary) && summary.valid) {
	out->draw_size_valid = summary.drawSizeValid ? 1 : 0;
	out->draw_size = (fastf_t)summary.drawSize;
    }
    out->mesh_lod = source->getMeshLod();
    SbVec3f lod_bmin;
    SbVec3f lod_bmax;
    if (source->getMeshLodBounds(lod_bmin, lod_bmax)) {
	out->mesh_lod_bounds_valid = 1;
	VSET(out->mesh_lod_bmin, lod_bmin[0], lod_bmin[1], lod_bmin[2]);
	VSET(out->mesh_lod_bmax, lod_bmax[0], lod_bmax[1], lod_bmax[2]);
    }
    return out->dbip ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_set_mesh_lod_for_path(
    struct ged *gedp,
    const char *path,
    struct BObolMeshLod *lod)
{
    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    SoBRLDatabaseSource *source = scene ?
	ged_obol_database_source_for_instance_key(scene, instance_keys.front()) :
	NULL;
    if (!source)
	return 0;

    source->setMeshLod(lod);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_set_mesh_lod_bounds_for_path(
    struct ged *gedp,
    const char *path,
    const point_t bmin,
    const point_t bmax)
{
    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;
	const int changed = source->setMeshLodBounds(
				SbVec3f((float)bmin[X], (float)bmin[Y],
				    (float)bmin[Z]),
				SbVec3f((float)bmax[X], (float)bmax[Y],
				    (float)bmax[Z]));
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_apply_record_for_path(
    struct ged *gedp,
    const char *path,
    const struct ged_draw_obol_database_source_record *record)
{
    if (!record)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_record_or_path(gedp,
	    path, record, record->draw_mode, 0);
    if (instance_keys.empty())
	return 0;

    const std::string source_instance_key = instance_keys.front();
    SoBRLDatabaseSource *source =
	ged_obol_database_source_for_instance_key(scene, source_instance_key);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    const int next_draw_mode =
	ged_obol_database_draw_mode_from_ged(record->draw_mode);
    int draw_mode_changed =
	scene->setDatabaseSourceInstanceDrawMode(source_instance_key.c_str(),
	    next_draw_mode);
    if (draw_mode_changed < 0)
	return 0;

    const char *representation_key = summary.representationKey.getString();
    if (!representation_key || !representation_key[0])
	representation_key = source_instance_key.c_str();
    int representation_changed =
	scene->setDatabaseSourceInstanceRepresentation(
	    source_instance_key.c_str(),
	    representation_key,
	    ged_obol_database_representation_mode_from_ged(
		record->draw_mode));
    if (representation_changed < 0)
	return 0;

    int state_changed = scene->setDatabaseSourceInstanceState(
			    source_instance_key.c_str(),
			    TRUE,
			    ged_obol_fold_revision(record->source_revision),
			    ged_obol_fold_revision(record->inputs_revision),
			    summary.visible,
			    summary.selected,
			    summary.highlighted,
			    summary.lineStyle,
			    summary.lineWidth,
			    summary.transparency,
			    summary.colorOverride,
			    summary.color,
			    summary.materialColorValid,
			    summary.materialColor,
			    summary.materialRevision);
    if (state_changed < 0)
	return 0;

    int material_policy_changed =
	scene->setDatabaseSourceInstanceMaterialPolicy(
	    source_instance_key.c_str(),
	    ged_obol_material_policy_from_ged(record->material_policy));
    if (material_policy_changed < 0)
	return 0;

    int realization_status = SoBRLDatabaseSource::UNREALIZED;
    if (record->realization_status ==
	GED_DRAW_OBOL_DATABASE_SOURCE_REALIZATION_CURRENT) {
	realization_status = SoBRLDatabaseSource::REALIZED;
    } else if (record->stale_reason == GED_DRAW_STALE_UPDATE_FAILED) {
	realization_status = SoBRLDatabaseSource::FAILED;
    }

    uint32_t stale_reason = SoBRLDatabaseSource::STALE_SOURCE;
    if (realization_status == SoBRLDatabaseSource::REALIZED)
	stale_reason = SoBRLDatabaseSource::STALE_NONE;
    else if (record->stale_reason != GED_DRAW_STALE_NONE)
	stale_reason = ged_obol_stale_reason_from_ged(record->stale_reason);

    int realization_changed =
	scene->setDatabaseSourceInstanceRealizationState(
	    source_instance_key.c_str(),
	    realization_status,
	    ged_obol_fold_revision(record->realized_source_revision),
	    ged_obol_fold_revision(record->realized_inputs_revision),
	    stale_reason,
	    realization_status == SoBRLDatabaseSource::FAILED ?
	    "GED source realization failed" : NULL);
    if (realization_changed < 0)
	return 0;

    int role_changed =
	scene->setDatabaseSourceInstanceRealizationRoleFlags(
	    source_instance_key.c_str(),
	    record->realization_role_flags);
    if (role_changed < 0)
	return 0;

    const uint32_t clamped_bot_threshold =
	record->realization_bot_threshold > UINT32_MAX ? UINT32_MAX :
	(uint32_t)record->realization_bot_threshold;
    int policy_changed =
	scene->setDatabaseSourceInstanceRealizationViewPolicy(
	    source_instance_key.c_str(),
	    record->realization_view_dependent ? TRUE : FALSE,
	    record->realization_csg_lod_enabled ? TRUE : FALSE,
	    record->realization_mesh_lod_enabled ? TRUE : FALSE,
	    (float)record->realization_view_scale,
	    (float)record->realization_lod_scale,
	    0,
	    0,
	    clamped_bot_threshold,
	    (float)record->realization_curve_scale,
	    (float)record->realization_point_scale);
    if (policy_changed < 0)
	return 0;

    return 1;
}

extern "C" int
ged_draw_obol_database_source_set_realization_for_path(
    struct ged *gedp,
    const char *path,
    int current,
    int failed,
    uint64_t realized_source_revision,
    uint64_t realized_inputs_revision,
    ged_draw_stale_reason reason)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int status = SoBRLDatabaseSource::UNREALIZED;
    if (failed)
	status = SoBRLDatabaseSource::FAILED;
    else if (current)
	status = SoBRLDatabaseSource::REALIZED;

    uint32_t stale_reason = SoBRLDatabaseSource::STALE_SOURCE;
    if (current && !failed)
	stale_reason = SoBRLDatabaseSource::STALE_NONE;
    else if (reason != GED_DRAW_STALE_NONE)
	stale_reason = ged_obol_stale_reason_from_ged(reason);

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int changed = scene->setDatabaseSourceInstanceRealizationState(
				source_instance_key.c_str(),
				status,
				ged_obol_fold_revision(realized_source_revision),
				ged_obol_fold_revision(realized_inputs_revision),
				stale_reason,
				failed ? "GED source realization failed" : NULL);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_realization_policy_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_obol_realization_policy_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    out->valid = 1;
    out->role_flags = summary.realizationRoleFlags;
    out->view_dependent = summary.realizationViewDependent ? 1 : 0;
    out->csg_lod_enabled = summary.realizationCsgLodEnabled ? 1 : 0;
    out->mesh_lod_enabled = summary.realizationMeshLodEnabled ? 1 : 0;
    out->view_scale = (fastf_t)summary.realizationViewScale;
    out->lod_scale = (fastf_t)summary.realizationLodScale;
    out->bot_threshold = (uint64_t)summary.realizationBotThreshold;
    out->curve_scale = (fastf_t)summary.realizationCurveScale;
    out->point_scale = (fastf_t)summary.realizationPointScale;
    return 1;
}

extern "C" int
ged_draw_obol_database_source_set_realization_roles_for_path(
    struct ged *gedp,
    const char *path,
    int role_flags)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int changed =
	    scene->setDatabaseSourceInstanceRealizationRoleFlags(
		source_instance_key.c_str(), role_flags);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_set_realization_view_policy_for_path(
    struct ged *gedp,
    const char *path,
    int view_dependent,
    int csg_lod_enabled,
    int mesh_lod_enabled,
    fastf_t view_scale,
    fastf_t lod_scale,
    int view_width,
    int view_height,
    uint64_t bot_threshold,
    fastf_t curve_scale,
    fastf_t point_scale)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    const uint32_t clamped_bot_threshold =
	bot_threshold > UINT32_MAX ? UINT32_MAX : (uint32_t)bot_threshold;
    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int changed =
	    scene->setDatabaseSourceInstanceRealizationViewPolicy(
		source_instance_key.c_str(),
		view_dependent ? TRUE : FALSE,
		csg_lod_enabled ? TRUE : FALSE,
		mesh_lod_enabled ? TRUE : FALSE,
		(float)view_scale,
		(float)lod_scale,
		view_width,
		view_height,
		clamped_bot_threshold,
		(float)curve_scale,
		(float)point_scale);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_realize_for_path(struct ged *gedp,
	const char *path)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	source = scene->findDatabaseSource(path);
    if (!source)
	return 0;

    if (source->needsRealization())
	(void)scene->realizePending();

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    return summary.realizationStatus == SoBRLDatabaseSource::REALIZED &&
	   !source->needsRealization() ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_sources_redraw(struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *path,
	int draw_mode)
{
    if (!gedp)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    std::vector<std::string> instance_keys;
    if (path && path[0]) {
	(void)ged_obol_append_exact_database_source_instance_keys(
	    instance_keys, scene, view_ctx, path, draw_mode);
	if (instance_keys.empty())
	    ged_obol_collect_database_source_instance_keys_matching(
		instance_keys, scene, view_ctx, path, 0, 1, draw_mode);
    } else {
	const int source_count = scene->getDatabaseSourceCount();
	if (!view_ctx && draw_mode < 0) {
	    int needs_realization = 0;
	    for (int i = 0; i < source_count; i++) {
		SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
		if (source && source->needsRealization()) {
		    needs_realization = 1;
		    break;
		}
	    }
	    if (needs_realization)
		(void)scene->realizePending();
	    BObolViewController *controller = ged_draw_obol_controller(gedp);
	    if (controller && source_count > 0)
		controller->requestRender(needs_realization ?
		    "ged-redraw-realized" : "ged-redraw-retained");
	    return source_count > 0 ? 1 : 0;
	}
	for (int i = 0; i < source_count; i++) {
	    BObolDatabaseSourceSummary summary;
	    if (!scene->getDatabaseSourceSummary(i, summary) ||
		!summary.valid ||
		!ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
		!ged_obol_database_source_summary_matches_mode(summary,
		    draw_mode))
		continue;
	    ged_obol_append_database_source_instance_key(instance_keys,
		summary);
	}
    }

    if (instance_keys.empty())
	return 0;

    int needs_realization = 0;
    for (const std::string &instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    scene->findDatabaseSourceInstance(instance_key.c_str());
	if (source && source->needsRealization()) {
	    needs_realization = 1;
	    break;
	}
    }
    if (needs_realization)
	(void)scene->realizePending();

    BObolViewController *controller = ged_draw_obol_controller(gedp);
    if (controller)
	controller->requestRender(needs_realization ?
	    "ged-redraw-realized" : "ged-redraw-retained");
    return 1;
}

extern "C" int
ged_draw_obol_database_source_realize_pending(struct ged *gedp)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    scene->realizePending();
    return 1;
}

static int
ged_bobol_database_source_ensure_for_path_impl(
    struct ged_bobol_publication_context *publication,
    struct ged *gedp,
    const char *path,
    struct db_i *dbip,
    int draw_mode,
    uint64_t source_revision)
{
    if (!gedp || !path || !path[0] || !dbip ||
	!ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const uint32_t folded_revision =
	source_revision ? ged_obol_fold_revision(source_revision) : 0;
    const bool explicit_publication = publication && publication->active;
    const int publication_mode = explicit_publication ?
	publication->draw_mode : -1;
    struct ged_view_context *publication_view = explicit_publication ?
	publication->view_ctx : NULL;
    if (explicit_publication && publication_mode >= 0) {
	/* Batched draw publication already has a canonical full path from the
	 * database walk.  Do not scan retained source records for every leaf:
	 * large initial draws otherwise become O(N^2) as the Obol tree grows.
	 * Exact existing instances are still preserved by ged_obol_replace_path's
	 * direct instance lookup.
	 */
	const int changed = ged_obol_replace_path(gedp,
			    publication_view, dbip, path,
			    publication_mode, folded_revision,
			    scene, 0, 1, NULL, NULL, NULL, NULL, NULL,
			    publication);
	return changed >= 0 ? 1 : 0;
    }

    const std::string instance_key =
	explicit_publication ?
	ged_obol_database_source_instance_key(
	    publication_view, path) :
	std::string(path);
    const int changed = (instance_key == path) ?
			scene->replaceDatabaseSource(path, dbip,
			    ged_obol_database_draw_mode_from_ged(draw_mode),
			    folded_revision) :
			scene->replaceDatabaseSourceInstance(instance_key.c_str(), path,
			    dbip, ged_obol_database_draw_mode_from_ged(draw_mode),
			    folded_revision);
    (void)ged_obol_apply_view_lod_policy(gedp, NULL, scene,
					 instance_key.c_str());

    BObolDatabaseSourceSummary summary;
    SoBRLDatabaseSource *source =
	scene->findDatabaseSourceInstance(instance_key.c_str());
    if ((source && source->getSummary(summary) && summary.valid) ||
	ged_obol_database_source_controller_summary_for_path(scene, path,
		summary)) {
	const char *parent_norm =
	    ged_obol_skip_leading_slash(
		summary.parentGroupPath.getString());
	if (!parent_norm || !parent_norm[0]) {
	    const std::string owner_group_path =
		ged_obol_database_source_owner_group_path_from_summary(summary);
	    if (!owner_group_path.empty()) {
		SoGroup *group = scene->ensureGroup(owner_group_path.c_str());
		if (group) {
		    const std::string intent_path =
			ged_obol_group_intent_path(owner_group_path.c_str());
		    (void)scene->setGroupDrawIntent(owner_group_path.c_str(),
						    intent_path.c_str(),
						    ged_obol_lod_draw_mode_from_ged(draw_mode),
						    BOBOL_LOD_DRAW_WIRE,
						    FALSE,
						    folded_revision);
		    (void)scene->moveDatabaseSourceInstanceToGroup(
			instance_key.c_str(),
			owner_group_path.c_str());
		}
	    }
	}
    }

    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_ensure_for_path(
    struct ged *gedp,
    const char *path,
    struct db_i *dbip,
    int draw_mode,
    uint64_t source_revision)
{
    return ged_bobol_database_source_ensure_for_path_impl(NULL, gedp, path,
	    dbip, draw_mode, source_revision);
}

static int
ged_bobol_database_source_ensure_for_path_with_placement_impl(
    struct ged_bobol_publication_context *publication,
    struct ged *gedp,
    const char *path,
    struct db_i *dbip,
    int draw_mode,
    uint64_t source_revision,
    int draw_mat_valid,
    const mat_t draw_mat,
    int draw_center_valid,
    const point_t draw_center,
    int draw_size_valid,
    fastf_t draw_size)
{
    if ((draw_mat_valid && !draw_mat) ||
	(draw_center_valid && !draw_center))
	return 0;

    ged_obol_publish_placement_state placement;
    if (draw_mat_valid || draw_center_valid || draw_size_valid) {
	placement.valid = true;
	if (draw_mat_valid) {
	    placement.drawMatrixValid = true;
	    placement.drawMatrix = ged_obol_sbmatrix_from_mat(draw_mat);
	}
	if (draw_center_valid) {
	    placement.drawCenterValid = true;
	    placement.drawCenter =
		SbVec3f(static_cast<float>(draw_center[0]),
			static_cast<float>(draw_center[1]),
			static_cast<float>(draw_center[2]));
	}
	if (draw_size_valid) {
	    placement.drawSizeValid = true;
	    placement.drawSize = static_cast<float>(draw_size);
	}
    }

    if (!placement.valid)
	return ged_bobol_database_source_ensure_for_path_impl(publication,
	    gedp, path, dbip, draw_mode, source_revision);

    if (!gedp || !path || !path[0] || !dbip ||
	!ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const uint32_t folded_revision =
	source_revision ? ged_obol_fold_revision(source_revision) : 0;
    const bool explicit_publication = publication && publication->active;
    const int publication_mode = explicit_publication ?
	publication->draw_mode : -1;
    struct ged_view_context *publication_view = explicit_publication ?
	publication->view_ctx : NULL;
    if (explicit_publication && publication_mode >= 0) {
	const int changed = ged_obol_replace_path(gedp,
			    publication_view, dbip, path,
			    publication_mode, folded_revision,
			    scene, 0, 1, NULL, &placement, NULL, NULL, NULL,
			    publication);
	return changed >= 0 ? 1 : 0;
    }

    if (!ged_bobol_database_source_ensure_for_path_impl(publication, gedp,
	    path, dbip, draw_mode, source_revision))
	return 0;

    return ged_draw_obol_database_source_set_placement_for_path(
	       gedp, path,
	       draw_mat_valid, draw_mat,
	       draw_center_valid, draw_center,
	       draw_size_valid, draw_size);
}

extern "C" int
ged_draw_obol_database_source_ensure_for_path_with_placement(
    struct ged *gedp,
    const char *path,
    struct db_i *dbip,
    int draw_mode,
    uint64_t source_revision,
    int draw_mat_valid,
    const mat_t draw_mat,
    int draw_center_valid,
    const point_t draw_center,
    int draw_size_valid,
    fastf_t draw_size)
{
    return ged_bobol_database_source_ensure_for_path_with_placement_impl(
	NULL, gedp, path, dbip, draw_mode, source_revision,
	draw_mat_valid, draw_mat, draw_center_valid, draw_center,
	draw_size_valid, draw_size);
}

int
ged_bobol_database_source_ensure_for_path(
    struct ged_bobol_publication_context *publication,
    const char *path,
    struct db_i *dbip,
    int draw_mode,
    uint64_t source_revision)
{
    if (!publication || !publication->active || !publication->gedp)
	return 0;
    return ged_bobol_database_source_ensure_for_path_impl(publication,
	publication->gedp, path, dbip, draw_mode, source_revision);
}

int
ged_bobol_database_source_ensure_for_path_with_placement(
    struct ged_bobol_publication_context *publication,
    const char *path,
    struct db_i *dbip,
    int draw_mode,
    uint64_t source_revision,
    int draw_mat_valid,
    const mat_t draw_mat,
    int draw_center_valid,
    const point_t draw_center,
    int draw_size_valid,
    fastf_t draw_size)
{
    if (!publication || !publication->active || !publication->gedp)
	return 0;
    return ged_bobol_database_source_ensure_for_path_with_placement_impl(
	publication, publication->gedp, path, dbip, draw_mode,
	source_revision, draw_mat_valid, draw_mat, draw_center_valid,
	draw_center, draw_size_valid, draw_size);
}

extern "C" int
ged_draw_obol_database_source_rename_for_path(
    struct ged *gedp,
    const char *path,
    const char *new_path,
    uint64_t source_revision)
{
    if (!gedp || !path || !path[0] || !new_path || !new_path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    BObolDatabaseSourceSummary source_summary;
    std::string old_owner_group_path =
	ged_obol_top_group_path_from_record_path(path);
    SoBRLDatabaseSource *first_source =
	ged_obol_database_source_for_instance_key(scene, instance_keys.front());
    if (first_source &&
	first_source->getSummary(source_summary) && source_summary.valid) {
	const std::string candidate_owner =
	    ged_obol_database_source_owner_group_path_from_summary(
		source_summary);
	if (old_owner_group_path.empty() &&
	    ged_obol_path_equal(candidate_owner.c_str(), path))
	    old_owner_group_path = candidate_owner;
    }

    int owner_draw_mode = BOBOL_LOD_DRAW_WIRE;
    int owner_fallback_draw_mode = BOBOL_LOD_DRAW_WIRE;
    SbBool owner_overlay = FALSE;
    SbBool owner_visible = TRUE;
    SbBool owner_selected = FALSE;
    SbBool owner_highlighted = FALSE;
    int owner_line_style = 0;
    int owner_line_width = 0;
    float owner_transparency = 0.0f;
    SbBool owner_color_override = FALSE;
    SbColor owner_color(1.0f, 1.0f, 1.0f);
    SbBool owner_material_color_valid = FALSE;
    SbColor owner_material_color(1.0f, 1.0f, 1.0f);
    uint32_t owner_material_revision = 0;
    uint32_t owner_revalidation_revision = 0;
    if (!old_owner_group_path.empty()) {
	SoGroup *old_owner_group =
	    scene->findGroup(old_owner_group_path.c_str());
	if (old_owner_group &&
	    old_owner_group->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	    SoBRLSceneGroup *scene_group =
		static_cast<SoBRLSceneGroup *>(old_owner_group);
	    owner_draw_mode = scene_group->drawMode.getValue();
	    owner_fallback_draw_mode =
		scene_group->fallbackDrawMode.getValue();
	    if (owner_fallback_draw_mode == BOBOL_LOD_DRAW_UNKNOWN)
		owner_fallback_draw_mode = BOBOL_LOD_DRAW_WIRE;
	    owner_overlay = scene_group->overlayIntent.getValue();
	    owner_visible = scene_group->visible.getValue();
	    owner_selected = scene_group->selected.getValue();
	    owner_highlighted = scene_group->highlighted.getValue();
	    owner_line_style = scene_group->lineStyle.getValue();
	    owner_line_width = scene_group->lineWidth.getValue();
	    owner_transparency = scene_group->transparency.getValue();
	    owner_color_override = scene_group->colorOverride.getValue();
	    owner_color = scene_group->color.getValue();
	    owner_material_color_valid =
		scene_group->materialColorValid.getValue();
	    owner_material_color = scene_group->materialColor.getValue();
	    owner_material_revision =
		scene_group->materialRevision.getValue();
	    owner_revalidation_revision =
		scene_group->revalidationRevision.getValue();
	}
    }

    const std::string new_owner_group_path =
	ged_obol_top_group_path_from_record_path(new_path);
    if (!old_owner_group_path.empty() && !new_owner_group_path.empty()) {
	std::string old_parent, old_leaf, new_parent, new_leaf;
	if (ged_obol_group_parent_leaf(old_owner_group_path, old_parent,
				       old_leaf) &&
	    ged_obol_group_parent_leaf(new_owner_group_path, new_parent,
				       new_leaf) &&
	    old_parent == new_parent &&
	    scene->renameGroup(old_owner_group_path.c_str(),
			       new_leaf.c_str()) > 0) {
	    old_owner_group_path = new_owner_group_path;
	}
    }

    const uint32_t folded_revision =
	source_revision ? ged_obol_fold_revision(source_revision) : 0;
    int changed = 0;
    std::vector<std::string> renamed_instance_keys;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	BObolDatabaseSourceSummary per_source_summary;
	const int have_per_source_summary =
	    source && source->getSummary(per_source_summary) &&
	    per_source_summary.valid;
	const SbBool per_source_visible =
	    (!have_per_source_summary || per_source_summary.visible) ?
	    TRUE : FALSE;
	const SbBool per_source_selected =
	    (have_per_source_summary && per_source_summary.selected) ?
	    TRUE : FALSE;
	const SbBool per_source_highlighted =
	    (have_per_source_summary && per_source_summary.highlighted) ?
	    TRUE : FALSE;
	const int per_source_line_style = have_per_source_summary ?
					  per_source_summary.lineStyle : 0;
	const int per_source_line_width = have_per_source_summary ?
					  per_source_summary.lineWidth : 0;
	const float per_source_transparency = have_per_source_summary ?
					      per_source_summary.transparency :
					      0.0f;
	const SbBool per_source_color_override =
	    (have_per_source_summary && per_source_summary.colorOverride) ?
	    TRUE : FALSE;
	const SbColor per_source_color = have_per_source_summary ?
					 per_source_summary.color :
					 SbColor(1.0f, 1.0f, 1.0f);
	const SbBool per_source_material_color_valid =
	    (have_per_source_summary && per_source_summary.materialColorValid) ?
	    TRUE : FALSE;
	const SbColor per_source_material_color = have_per_source_summary ?
						  per_source_summary.materialColor :
						  SbColor(1.0f, 1.0f, 1.0f);
	const uint32_t per_source_material_revision =
	    have_per_source_summary ? per_source_summary.materialRevision : 0;

	const std::string new_source_instance_key =
	    ged_obol_renamed_database_source_instance_key(source_instance_key,
		path, new_path);
	const int renamed = scene->renameDatabaseSourceInstance(
				source_instance_key.c_str(),
				new_source_instance_key.c_str(),
				new_path,
				folded_revision);
	if (renamed < 0)
	    return 0;
	if (renamed > 0) {
	    changed = 1;
	    ged_obol_append_unique_path(renamed_instance_keys,
		new_source_instance_key.c_str());
	    BObolDatabaseSourceSummary renamed_summary;
	    SoBRLDatabaseSource *renamed_source =
		ged_obol_database_source_for_instance_key(scene,
		    new_source_instance_key);
	    if (renamed_source &&
		renamed_source->getSummary(renamed_summary) &&
		renamed_summary.valid) {
		(void)scene->setDatabaseSourceInstanceState(
		    new_source_instance_key.c_str(),
		    FALSE,
		    renamed_summary.sourceRevision,
		    renamed_summary.inputsRevision,
		    per_source_visible,
		    per_source_selected,
		    per_source_highlighted,
		    per_source_line_style,
		    per_source_line_width,
		    per_source_transparency,
		    per_source_color_override,
		    per_source_color,
		    per_source_material_color_valid,
		    per_source_material_color,
		    per_source_material_revision);
	    }
	}
    }
    if (changed > 0 && !old_owner_group_path.empty()) {
	if (!new_owner_group_path.empty()) {
	    SoGroup *new_owner_group =
		scene->ensureGroup(new_owner_group_path.c_str());
	    if (new_owner_group &&
		new_owner_group->isOfType(
		    SoBRLSceneGroup::getClassTypeId())) {
		const std::string intent_path =
		    ged_obol_group_intent_path(new_owner_group_path.c_str());
		(void)scene->setGroupDrawIntent(
		    new_owner_group_path.c_str(),
		    intent_path.c_str(),
		    owner_draw_mode,
		    owner_fallback_draw_mode,
		    owner_overlay,
		    owner_revalidation_revision);
		(void)scene->setGroupDisplayState(
		    new_owner_group_path.c_str(),
		    owner_visible,
		    owner_selected,
		    owner_highlighted,
		    owner_line_style,
		    owner_line_width,
		    owner_transparency,
		    owner_color_override,
		    owner_color,
		    owner_material_color_valid,
		    owner_material_color,
		    owner_material_revision);
		for (const std::string &renamed_instance_key :
		     renamed_instance_keys) {
		    (void)scene->moveDatabaseSourceInstanceToGroup(
			renamed_instance_key.c_str(),
			new_owner_group_path.c_str());
		}
		(void)ged_obol_prune_empty_groups(scene);
	    }
	}
    }
    if (changed > 0)
	scene->realizePending();
    return changed > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_move_to_group_for_path(
    struct ged *gedp,
    const char *source_path,
    const char *group_path)
{
    if (!gedp || !source_path || !source_path[0] ||
	!group_path || !group_path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string target_group =
	ged_obol_group_path_from_record_path(group_path);
    if (target_group.empty())
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp,
	    source_path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int moved = scene->moveDatabaseSourceInstanceToGroup(
			      source_instance_key.c_str(),
			      target_group.c_str());
	if (moved < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

static bool
ged_obol_group_path_is_overlay(BObolSceneController *scene,
			       const char *group_path)
{
    if (!scene || !group_path || !group_path[0] ||
	BU_STR_EQUAL(group_path, "/"))
	return false;

    std::string current = ged_obol_skip_leading_slash(group_path);
    while (!current.empty()) {
	SoGroup *group = scene->findGroup(current.c_str());
	if (group && group->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	    SoBRLSceneGroup *scene_group =
		static_cast<SoBRLSceneGroup *>(group);
	    if (scene_group->overlayIntent.getValue())
		return true;
	}

	const size_t slash = current.rfind('/');
	if (slash == std::string::npos)
	    break;
	current = current.substr(0, slash);
    }

    return false;
}

extern "C" int
ged_draw_obol_database_source_count(
    struct ged *gedp,
    int skip_overlay_groups,
    size_t *out)
{
    if (out)
	*out = 0;
    if (!out || !gedp)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int source_count = scene->getDatabaseSourceCount();
    if (source_count <= 0)
	return 1;

    if (!skip_overlay_groups) {
	*out = (size_t)source_count;
	return 1;
    }

    size_t count = 0;
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;
	if (ged_obol_group_path_is_overlay(scene,
					   summary.parentGroupPath.getString()))
	    continue;
	count++;
    }

    *out = count;
    return 1;
}

extern "C" int
ged_draw_obol_database_source_paths_foreach(
    struct ged *gedp,
    int skip_overlay_groups,
    ged_draw_obol_database_source_path_cb cb,
    void *userdata)
{
    if (!gedp || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    if (source_count <= 0)
	return 1;

    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	if (skip_overlay_groups &&
	    ged_obol_group_path_is_overlay(scene,
					   summary.parentGroupPath.getString()))
	    continue;

	const char *source_path = summary.path.getString();
	if (source_path && source_path[0] && !(*cb)(gedp, source_path,
		userdata))
	    return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_database_source_records_foreach(
    struct ged *gedp,
    int skip_overlay_groups,
    ged_draw_obol_database_source_record_cb cb,
    void *userdata)
{
    if (!gedp || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    if (source_count <= 0)
	return 1;

    std::vector<BObolDatabaseSourceSummary> summaries;
    summaries.reserve(static_cast<size_t>(source_count));
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	if (skip_overlay_groups &&
	    ged_obol_group_path_is_overlay(scene,
					   summary.parentGroupPath.getString()))
	    continue;

	summaries.push_back(summary);
    }

    std::stable_sort(summaries.begin(), summaries.end(),
	[scene](const BObolDatabaseSourceSummary &a,
		const BObolDatabaseSourceSummary &b) {
	const auto active_in_parent = [scene](
	    const BObolDatabaseSourceSummary &summary) -> int {
	    const char *group_path = summary.parentGroupPath.getString();
	    SoGroup *group = (group_path && group_path[0]) ?
			     scene->findGroup(group_path) : NULL;
	    if (!group ||
		!group->isOfType(SoBRLSceneGroup::getClassTypeId()))
		return 0;
	    SoBRLSceneGroup *scene_group =
		static_cast<SoBRLSceneGroup *>(group);
	    const int group_mode =
		ged_obol_lod_draw_mode_to_ged(
		    scene_group->drawMode.getValue());
	    return group_mode ==
		   ged_obol_database_source_summary_ged_mode(summary) ? 1 : 0;
	};
	const int a_active = active_in_parent(a);
	const int b_active = active_in_parent(b);
	if (a_active != b_active)
	    return a_active > b_active;
	return false;
    });

    for (const BObolDatabaseSourceSummary &summary : summaries) {
	struct ged_draw_obol_database_source_record record;
	if (!ged_obol_database_source_record_from_summary(&record, summary))
	    continue;
	if (record.database_path && record.database_path[0] &&
	    !(*cb)(gedp, &record, userdata))
	    return 0;

	SoBRLDatabaseSource *source = scene->findDatabaseSourceInstance(
	    summary.instanceKey.getString());
	if (!source || !source->hasCompactInstanceIndex())
	    continue;
	for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary occurrence;
	    if (!source->getCompactInstanceHandle(i, handle) ||
		!source->getCompactInstanceSummary(handle, occurrence) ||
		occurrence.path.getLength() == 0)
		continue;
	    struct ged_draw_obol_database_source_record occurrence_record =
		record;
	    occurrence_record.database_path = occurrence.path.getString();
	    occurrence_record.visible = record.visible && occurrence.visible;
	    occurrence_record.selected = occurrence.selected ? 1 : 0;
	    occurrence_record.highlighted = occurrence.highlighted ? 1 : 0;
	    if (!(*cb)(gedp, &occurrence_record, userdata))
		return 0;
	}
    }

    return 1;
}

extern "C" int
ged_draw_obol_visible_database_source_records_foreach_fast(
    struct ged *gedp,
    ged_draw_obol_database_source_record_cb cb,
    void *userdata)
{
    if (!gedp || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (!source || !source->visible.getValue())
	    continue;

	const char *path = source->path.getValue().getString();
	const char *instance_key = source->instanceKey.getValue().getString();
	if (!path || !path[0])
	    continue;
	if (!instance_key || !instance_key[0])
	    instance_key = path;

	struct ged_draw_obol_database_source_record record;
	memset(&record, 0, sizeof(record));
	record.valid = 1;
	record.database_path = path;
	record.source_path = path;
	record.instance_key = instance_key;
	record.visible = 1;
	const int representation_mode = source->representationMode.getValue();
	record.draw_mode = representation_mode >= 0 ? representation_mode :
	    ged_obol_database_draw_mode_to_ged(source->drawMode.getValue());
	if (!(*cb)(gedp, &record, userdata))
	    return 0;
    }

    return 1;
}

static SbBool
ged_obol_shape_node_bool_field(SoNode *node,
			       const char *field_name,
			       SbBool fallback)
{
    if (!node || !field_name)
	return fallback;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	if (BU_STR_EQUAL(field_name, "databaseIntent"))
	    return shape->databaseIntent.getValue();
	if (BU_STR_EQUAL(field_name, "overlayIntent"))
	    return shape->overlayIntent.getValue();
	if (BU_STR_EQUAL(field_name, "localSource"))
	    return shape->localSource.getValue();
	if (BU_STR_EQUAL(field_name, "sharedSource"))
	    return shape->sharedSource.getValue();
	if (BU_STR_EQUAL(field_name, "nonDatabaseSource"))
	    return shape->nonDatabaseSource.getValue();
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	if (BU_STR_EQUAL(field_name, "databaseIntent"))
	    return shape->databaseIntent.getValue();
	if (BU_STR_EQUAL(field_name, "overlayIntent"))
	    return shape->overlayIntent.getValue();
	if (BU_STR_EQUAL(field_name, "localSource"))
	    return shape->localSource.getValue();
	if (BU_STR_EQUAL(field_name, "sharedSource"))
	    return shape->sharedSource.getValue();
	if (BU_STR_EQUAL(field_name, "nonDatabaseSource"))
	    return shape->nonDatabaseSource.getValue();
    }

    return fallback;
}

static const char *
ged_obol_shape_node_record_role(SoNode *node)
{
    if (!node)
	return "";

    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<SoBRLVListShape *>(node)->
	       recordRole.getValue().getString();

    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<SoBRLMeshShape *>(node)->
	       recordRole.getValue().getString();

    return "";
}

static int
ged_obol_shape_node_is_database_realization(SoNode *node)
{
    if (!node)
	return 0;

    const char *role = ged_obol_shape_node_record_role(node);
    return ged_obol_shape_node_bool_field(node, "databaseIntent", FALSE) &&
	   !ged_obol_shape_node_bool_field(node, "localSource", FALSE) &&
	   !ged_obol_shape_node_bool_field(node, "sharedSource", FALSE) &&
	   !ged_obol_shape_node_bool_field(node, "nonDatabaseSource", FALSE) &&
	   role && BU_STR_EQUAL(role, "database");
}

static int
ged_obol_shape_node_is_auxiliary_record(SoNode *node)
{
    const char *role = ged_obol_shape_node_record_role(node);
    return role && BU_STR_EQUAL(role, "auxiliary");
}

static int
ged_draw_obol_shape_paths_foreach_node(
    struct ged *gedp,
    SoNode *node,
    int skip_overlay_groups,
    ged_draw_obol_database_source_path_cb cb,
    void *userdata)
{
    if (!node)
	return 1;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return 1;

    if (skip_overlay_groups &&
	node->isOfType(SoBRLSceneGroup::getClassTypeId()) &&
	static_cast<SoBRLSceneGroup *>(node)->overlayIntent.getValue())
	return 1;

    if (node->isOfType(SoBRLVListShape::getClassTypeId()) ||
	node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	if (ged_obol_shape_node_is_database_realization(node) ||
	    ged_obol_shape_node_is_auxiliary_record(node))
	    return 1;
	if (skip_overlay_groups &&
	    ged_obol_shape_node_bool_field(node, "overlayIntent", FALSE))
	    return 1;

	const char *shape_path = NULL;
	if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	    shape_path = static_cast<SoBRLVListShape *>(node)->
			 sourcePath.getValue().getString();
	else
	    shape_path = static_cast<SoBRLMeshShape *>(node)->
			 sourcePath.getValue().getString();

	if (shape_path && shape_path[0] && !(*cb)(gedp, shape_path, userdata))
	    return 0;
	return 1;
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return 1;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (!ged_draw_obol_shape_paths_foreach_node(gedp, group->getChild(i),
	    skip_overlay_groups, cb, userdata))
	    return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_shape_paths_foreach(
    struct ged *gedp,
    int skip_overlay_groups,
    ged_draw_obol_database_source_path_cb cb,
    void *userdata)
{
    if (!gedp || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    return ged_draw_obol_shape_paths_foreach_node(gedp, scene->getSceneRoot(),
	    skip_overlay_groups, cb, userdata);
}

extern "C" int
ged_draw_obol_group_database_source_paths_foreach(
    struct ged *gedp,
    const char *group_path,
    ged_draw_obol_database_source_path_cb cb,
    void *userdata)
{
    if (!gedp || !group_path || !group_path[0] || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const std::string target_group =
	ged_obol_group_path_from_record_path(group_path);
    if (target_group.empty())
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	const std::string owner_group_path =
	    ged_obol_database_source_owner_group_path_from_summary(summary);
	if (!ged_obol_path_equal(owner_group_path.c_str(),
				 target_group.c_str()) &&
	    !ged_obol_path_has_prefix(owner_group_path.c_str(),
				      target_group.c_str()))
	    continue;

	const char *source_path = summary.path.getString();
	if (source_path && source_path[0] && !(*cb)(gedp, source_path,
		userdata))
	    return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_group_database_source_records_foreach(
    struct ged *gedp,
    const char *group_path,
    ged_draw_obol_database_source_record_cb cb,
    void *userdata)
{
    if (!gedp || !group_path || !group_path[0] || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const std::string target_group =
	ged_obol_group_path_from_record_path(group_path);
    if (target_group.empty())
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	const std::string owner_group_path =
	    ged_obol_database_source_owner_group_path_from_summary(summary);
	if (!ged_obol_path_equal(owner_group_path.c_str(),
				 target_group.c_str()) &&
	    !ged_obol_path_has_prefix(owner_group_path.c_str(),
				      target_group.c_str()))
	    continue;

	struct ged_draw_obol_database_source_record record;
	if (!ged_obol_database_source_record_from_summary(&record, summary))
	    continue;
	if (record.database_path && record.database_path[0] &&
	    !(*cb)(gedp, &record, userdata))
	    return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_database_source_owner_group_path_for_path(
    struct ged *gedp,
    const char *path,
    struct bu_vls *out)
{
    if (out)
	bu_vls_trunc(out, 0);
    if (!gedp || !path || !path[0] || !out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    BObolDatabaseSourceSummary summary;
    if (!ged_obol_database_source_controller_summary_for_path(scene, path,
	    summary))
	return 0;

    const std::string owner_group_path =
	ged_obol_database_source_owner_group_path_from_summary(summary);
    if (owner_group_path.empty())
	return 0;

    bu_vls_strcpy(out, owner_group_path.c_str());
    return 1;
}

extern "C" int
ged_draw_obol_database_source_remove_for_path(
    struct ged *gedp,
    const char *path)
{
    if (!gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> paths;
    ged_obol_append_unique_path(paths, path);
    return ged_obol_remove_paths(paths, NULL, scene);
}

extern "C" int
ged_draw_obol_database_sources_remove_for_path_prefix(
    struct ged *gedp,
    const char *path_prefix)
{
    if (!gedp || !path_prefix || !path_prefix[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> paths;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (!source)
	    continue;
	const char *source_path = source->path.getValue().getString();
	if (ged_obol_path_has_prefix(source_path, path_prefix))
	    ged_obol_append_unique_path(paths, source_path);
    }
    return ged_obol_remove_paths(paths, NULL, scene);
}

extern "C" int
ged_draw_obol_database_sources_remove_for_path_prefix_in_scope(
    struct ged *gedp,
    const char *path_prefix,
    struct ged_view_context *view_ctx)
{
    return ged_draw_obol_database_sources_remove_for_path_prefix_in_scope_mode(
	       gedp, path_prefix, view_ctx, -1);
}

extern "C" int
ged_draw_obol_database_sources_remove_for_path_prefix_in_scope_mode(
    struct ged *gedp,
    const char *path_prefix,
    struct ged_view_context *view_ctx,
    int mode)
{
    if (!gedp || !path_prefix || !path_prefix[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	const char *source_path = summary.path.getString();
	if (ged_obol_path_has_prefix(source_path, path_prefix))
	    ged_obol_append_database_source_instance_key(instance_keys,
		    summary);
    }
    return ged_obol_remove_instance_keys(instance_keys, scene);
}

extern "C" int
ged_draw_obol_active_database_sources_remove_for_path_prefix(
    struct ged *gedp,
    const char *path_prefix)
{
    return ged_draw_obol_active_database_sources_remove_for_path_prefix_mode(
	       gedp, path_prefix, -1);
}

extern "C" int
ged_draw_obol_active_database_sources_remove_for_path_prefix_mode(
    struct ged *gedp,
    const char *path_prefix,
    int mode)
{
    if (!gedp || !path_prefix || !path_prefix[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	const char *source_path = summary.path.getString();
	if (ged_obol_path_has_prefix(source_path, path_prefix))
	    ged_obol_append_database_source_instance_key(instance_keys,
		    summary);
    }
    return ged_obol_remove_instance_keys(instance_keys, scene);
}

extern "C" int
ged_draw_obol_database_sources_remove_for_component_name(
    struct ged *gedp,
    const char *name,
    int nonroot_only,
    int mode)
{
    if (!gedp || !name || !name[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> targets;
    ged_obol_append_unique_path(targets, name);

    std::vector<std::string> paths;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;
	if (!ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	for (const std::string &target : targets) {
	    if (ged_obol_path_has_component_name(summary.path.getString(),
						 target.c_str(), nonroot_only ? 1 : 0)) {
		ged_obol_append_unique_path(paths, summary.path.getString());
		break;
	    }
	}
    }

    return ged_obol_remove_paths(paths, NULL, scene, mode);
}

extern "C" int
ged_draw_obol_database_sources_clear(struct ged *gedp)
{
    if (!gedp)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    int changed = (scene->clearDatabaseSources() > 0) ? 1 : 0;
    if (ged_obol_prune_empty_groups(scene))
	changed = 1;
    if (changed)
	scene->realizePending();
    return changed;
}

extern "C" int
ged_draw_obol_database_sources_clear_in_scope(struct ged *gedp,
	struct ged_view_context *view_ctx)
{
    if (!gedp)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    int changed = ged_obol_clear_database_sources_in_scope(scene, view_ctx);
    if (changed)
	scene->realizePending();
    return changed;
}

static int
ged_obol_scene_clear_controller(BObolSceneController *scene)
{
    if (!scene)
	return 0;

    int changed = (scene->clearDatabaseSources() > 0) ? 1 : 0;

    std::vector<std::string> group_paths;
    const int summary_count = scene->getSceneTreeSummaryCount();
    for (int i = 0; i < summary_count; i++) {
	BObolSceneTreeSummary tree_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !tree_summary.valid ||
	    tree_summary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	    tree_summary.path.getLength() == 0 ||
	    BU_STR_EQUAL(tree_summary.path.getString(), "/"))
	    continue;
	group_paths.push_back(tree_summary.path.getString());
    }

    std::sort(group_paths.begin(), group_paths.end(),
    [](const std::string &a, const std::string &b) {
	return a.size() > b.size();
    });
    for (const std::string &path : group_paths) {
	if (scene->removeGroup(path.c_str()) > 0)
	    changed = 1;
    }

    if (changed)
	scene->realizePending();
    return changed;
}

extern "C" int
ged_draw_obol_scene_clear(struct ged *gedp)
{
    if (!gedp)
	return 0;

    return ged_obol_scene_clear_controller(
	       ged_draw_obol_scene_controller(gedp));
}

extern "C" int
ged_draw_obol_active_scene_clear(struct ged *gedp)
{
    if (!gedp)
	return 0;

    return ged_obol_scene_clear_controller(
	       ged_draw_obol_scene_controller(gedp));
}

extern "C" int
ged_draw_obol_groups_remove_for_component_name(
    struct ged *gedp,
    const char *name)
{
    if (!gedp || !name || !name[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> paths;
    const int summary_count = scene->getSceneTreeSummaryCount();
    for (int i = summary_count - 1; i >= 0; i--) {
	BObolSceneTreeSummary tree_summary;
	BObolSceneDisplaySummary display_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !scene->getSceneDisplaySummary(i, display_summary) ||
	    !tree_summary.valid ||
	    !display_summary.valid ||
	    tree_summary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	    tree_summary.path.getLength() == 0 ||
	    BU_STR_EQUAL(tree_summary.path.getString(), "/") ||
	    !display_summary.hasDrawIntent ||
	    !ged_obol_intent_is_ged_draw_group(
		display_summary.intentPath) ||
	    !ged_obol_path_has_component_name(
		tree_summary.path.getString(), name, 0))
	    continue;

	SoGroup *group = scene->findGroup(tree_summary.path.getString());
	if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    continue;
	SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
	if (scene_group->overlayIntent.getValue())
	    continue;

	ged_obol_append_unique_path(paths, tree_summary.path.getString());
    }

    int changed = 0;
    for (const std::string &path : paths) {
	if (scene->removeGroup(path.c_str()) > 0)
	    changed = 1;
    }
    if (changed)
	scene->realizePending();
    return changed;
}

extern "C" int
ged_draw_obol_group_remove_for_path(
    struct ged *gedp,
    const char *path)
{
    if (!gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    if (group_path.empty())
	return 0;

    const int removed = scene->removeGroup(group_path.c_str());
    if (removed > 0)
	scene->realizePending();
    return removed > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_clear_for_path(
    struct ged *gedp,
    const char *path)
{
    if (!gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    if (group_path.empty())
	return 0;

    const int cleared = scene->clearGroup(group_path.c_str());
    if (cleared > 0)
	scene->realizePending();
    return cleared > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_rename_for_path(
    struct ged *gedp,
    const char *path,
    const char *new_path)
{
    if (!gedp || !path || !path[0] || !new_path || !new_path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    const std::string target_group_path =
	ged_obol_group_path_from_record_path(new_path);
    if (group_path.empty() || target_group_path.empty())
	return 0;

    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    if (group_path == target_group_path)
	return 1;

    std::string old_parent, old_leaf, new_parent, new_leaf;
    if (!ged_obol_group_parent_leaf(group_path, old_parent, old_leaf) ||
	!ged_obol_group_parent_leaf(target_group_path, new_parent,
				    new_leaf) ||
	old_parent != new_parent)
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    int group_draw_mode = scene_group->drawMode.getValue();
    int fallback_draw_mode = scene_group->fallbackDrawMode.getValue();
    if (fallback_draw_mode == BOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BOBOL_LOD_DRAW_WIRE;
    SbBool overlay = scene_group->overlayIntent.getValue();
    uint32_t revalidation_revision =
	scene_group->revalidationRevision.getValue();

    if (scene->renameGroup(group_path.c_str(), new_leaf.c_str()) <= 0)
	return 0;

    group = scene->findGroup(target_group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    const std::string intent_path =
	ged_obol_group_intent_path(target_group_path.c_str());
    int changed = scene->setGroupDrawIntent(target_group_path.c_str(),
					    intent_path.c_str(), group_draw_mode, fallback_draw_mode, overlay,
					    revalidation_revision);
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_erase_subpath_for_path(
    struct ged *gedp,
    const char *parent_path,
    const char *subpath)
{
    if (!gedp || !parent_path || !parent_path[0] ||
	!subpath || !subpath[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path =
	ged_obol_group_path_from_record_path(parent_path);
    const std::string relative_path =
	ged_obol_group_path_from_record_path(subpath);
    if (group_path.empty() || relative_path.empty())
	return 0;

    const int erased = scene->eraseGroupSubpath(group_path.c_str(),
		       relative_path.c_str());
    if (erased > 0)
	scene->realizePending();
    return erased > 0 ? 1 : 0;
}


static int ged_obol_collect_view_database_sources(
    struct ged_view_context *view_ctx,
    BObolViewController *controller,
    void *userdata);

static std::set<SoBRLDatabaseSource *>
ged_obol_attached_database_sources(struct ged *gedp)
{
    std::set<SoBRLDatabaseSource *> sources;
    BObolSceneController *owned = ged_draw_obol_scene_controller(gedp);
    if (owned) {
	for (int i = 0; i < owned->getDatabaseSourceCount(); i++)
	    sources.insert(owned->getDatabaseSource(i));
    }
    ged_bobol_view_controllers_foreach(gedp,
	ged_obol_collect_view_database_sources, &sources);
    return sources;
}


extern "C" int
ged_draw_obol_source_visibility_frontier_set(
    struct ged *gedp,
    const char *root_path,
    struct ged_view_context *view_ctx,
    int mode,
    const char *const *paths,
    size_t path_count)
{
    if (!gedp || !root_path || !root_path[0])
	return 0;

    std::vector<SbString> frontier;
    frontier.reserve(path_count);
    for (size_t i = 0; paths && i < path_count; i++) {
	if (paths[i] && paths[i][0])
	    frontier.push_back(SbString(paths[i]));
    }

    int changed = 0;
    const std::set<SoBRLDatabaseSource *> sources =
	ged_obol_attached_database_sources(gedp);
    for (SoBRLDatabaseSource *source : sources) {
	BObolDatabaseSourceSummary summary;
	if (!source || !source->getSummary(summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	const char *source_path = summary.path.getString();
	if (!ged_obol_path_equal(source_path, root_path) &&
	    ged_obol_semantic_path_string(source_path) !=
		ged_obol_semantic_path_string(root_path))
	    continue;
	if (source->setCompactInstanceVisibilityFrontier(frontier) > 0)
	    changed++;
    }
    return changed;
}


extern "C" int
ged_draw_obol_source_visibility_overrides_set(
    struct ged *gedp,
    const char *root_path,
    struct ged_view_context *view_ctx,
    int mode,
    const char *const *paths,
    const int *visible,
    size_t rule_count)
{
    if (!gedp || !root_path || !root_path[0] ||
	(rule_count && (!paths || !visible)))
	return 0;

    std::vector<SbString> rule_paths;
    std::vector<SbBool> rule_states;
    rule_paths.reserve(rule_count);
    rule_states.reserve(rule_count);
    for (size_t i = 0; i < rule_count; i++) {
	if (!paths[i] || !paths[i][0])
	    continue;
	rule_paths.push_back(SbString(paths[i]));
	rule_states.push_back(visible[i] ? TRUE : FALSE);
    }

    int changed = 0;
    const std::set<SoBRLDatabaseSource *> sources =
	ged_obol_attached_database_sources(gedp);
    for (SoBRLDatabaseSource *source : sources) {
	BObolDatabaseSourceSummary summary;
	if (!source || !source->getSummary(summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	const char *source_path = summary.path.getString();
	if (!ged_obol_path_equal(source_path, root_path) &&
	    ged_obol_semantic_path_string(source_path) !=
		ged_obol_semantic_path_string(root_path))
	    continue;
	if (source->setCompactInstanceVisibilityOverrides(
		rule_paths, rule_states) > 0)
	    changed++;
    }
    return changed;
}


extern "C" int
ged_draw_obol_source_visibility_frontier_clear(
    struct ged *gedp,
    const char *root_path,
    struct ged_view_context *view_ctx,
    int mode)
{
    if (!gedp || !root_path || !root_path[0])
	return 0;

    int changed = 0;
    const std::set<SoBRLDatabaseSource *> sources =
	ged_obol_attached_database_sources(gedp);
    for (SoBRLDatabaseSource *source : sources) {
	BObolDatabaseSourceSummary summary;
	if (!source || !source->getSummary(summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary, mode))
	    continue;
	const char *source_path = summary.path.getString();
	if (!ged_obol_path_equal(source_path, root_path) &&
	    ged_obol_semantic_path_string(source_path) !=
		ged_obol_semantic_path_string(root_path))
	    continue;
	if (source->clearCompactInstanceVisibilityFrontier() > 0)
	    changed++;
    }
    return changed;
}

static int
ged_bobol_database_source_update_display_for_path_impl(
    struct ged_bobol_publication_context *publication,
    struct ged *gedp,
    const char *path,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision)
{
    if ((color_valid && !color) ||
	(material_color_valid && !material_color))
	return 0;

    BObolSceneController *scene =
	(publication && publication->active) ? publication->scene :
	ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0,
	    publication);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;

	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    return 0;

	SbBool nextVisible = visible_valid ?
			     (visible ? TRUE : FALSE) : summary.visible;
	SbBool nextSelected = selected_valid ?
			      (selected ? TRUE : FALSE) : summary.selected;
	SbBool nextHighlighted = highlighted_valid ?
				 (highlighted ? TRUE : FALSE) :
				 summary.highlighted;
	int nextDrawMode = draw_mode_valid ?
			   ged_obol_database_draw_mode_from_ged(draw_mode) :
			   summary.drawMode;
	int nextLineStyle = line_style_valid ?
			    line_style : summary.lineStyle;
	int nextLineWidth = line_width_valid ?
			    line_width : summary.lineWidth;
	float nextTransparency = transparency_valid ?
				 static_cast<float>(transparency) :
				 summary.transparency;

	SbBool nextColorOverride = summary.colorOverride;
	SbColor nextColor = summary.color;
	if (color_valid) {
	    nextColorOverride = TRUE;
	    nextColor = ged_obol_color_from_rgb(color);
	}

	SbBool nextMaterialColorValid = summary.materialColorValid;
	SbColor nextMaterialColor = summary.materialColor;
	uint32_t nextMaterialRevision = summary.materialRevision;
	if (material_color_valid) {
	    nextMaterialColorValid = TRUE;
	    nextMaterialColor = ged_obol_color_from_rgb(material_color);
	    if (material_revision_valid) {
		nextMaterialRevision =
		    ged_obol_fold_revision(material_revision);
	    } else {
		nextMaterialRevision++;
		if (!nextMaterialRevision)
		    nextMaterialRevision = 1;
	    }
	}

	int draw_mode_changed = scene->setDatabaseSourceInstanceDrawMode(
				    source_instance_key.c_str(), nextDrawMode);
	if (draw_mode_changed < 0)
	    return 0;

	if (draw_mode_valid) {
	    const char *representation_key =
		summary.representationKey.getString();
	    if (!representation_key || !representation_key[0])
		representation_key = source_instance_key.c_str();
	    int representation_changed =
		scene->setDatabaseSourceInstanceRepresentation(
		    source_instance_key.c_str(),
		    representation_key,
		    ged_obol_database_representation_mode_from_ged(
			draw_mode));
	    if (representation_changed < 0)
		return 0;
	}

	int changed = scene->setDatabaseSourceInstanceState(
			  source_instance_key.c_str(),
			  FALSE,
			  summary.sourceRevision,
			  summary.inputsRevision,
			  nextVisible,
			  nextSelected,
			  nextHighlighted,
			  nextLineStyle,
			  nextLineWidth,
			  nextTransparency,
			  nextColorOverride,
			  nextColor,
			  nextMaterialColorValid,
			  nextMaterialColor,
			  nextMaterialRevision);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_update_display_for_path(
    struct ged *gedp,
    const char *path,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision)
{
    return ged_bobol_database_source_update_display_for_path_impl(NULL, gedp,
	path, visible_valid, visible, selected_valid, selected,
	highlighted_valid, highlighted, draw_mode_valid, draw_mode,
	line_style_valid, line_style, line_width_valid, line_width,
	transparency_valid, transparency, color_valid, color,
	material_color_valid, material_color, material_revision_valid,
	material_revision);
}

int
ged_bobol_database_source_update_display_for_path(
    struct ged_bobol_publication_context *publication,
    const char *path,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision)
{
    if (!publication || !publication->active || !publication->gedp)
	return 0;
    return ged_bobol_database_source_update_display_for_path_impl(publication,
	publication->gedp, path, visible_valid, visible, selected_valid,
	selected, highlighted_valid, highlighted, draw_mode_valid, draw_mode,
	line_style_valid, line_style, line_width_valid, line_width,
	transparency_valid, transparency, color_valid, color,
	material_color_valid, material_color, material_revision_valid,
	material_revision);
}

extern "C" int
ged_draw_obol_database_source_set_selected_for_instance_key(
    struct ged *gedp,
    const char *instance_key,
    int selected)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !instance_key || !instance_key[0])
	return 0;

    SoBRLDatabaseSource *source =
	scene->findDatabaseSourceInstance(instance_key);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    const int changed = scene->setDatabaseSourceInstanceState(
			    instance_key,
			    FALSE,
			    summary.sourceRevision,
			    summary.inputsRevision,
			    summary.visible,
			    selected ? TRUE : FALSE,
			    summary.highlighted,
			    summary.lineStyle,
			    summary.lineWidth,
			    summary.transparency,
			    summary.colorOverride,
			    summary.color,
			    summary.materialColorValid,
			    summary.materialColor,
			    summary.materialRevision);
    return changed > 0 ? 1 : 0;
}

static void
ged_obol_collect_database_sources(SoNode *node,
	std::set<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	sources.insert(static_cast<SoBRLDatabaseSource *>(node));
	return;
    }
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	ged_obol_collect_database_sources(group->getChild(i), sources);
}

static int
ged_obol_collect_view_database_sources(
    struct ged_view_context *UNUSED(view_ctx),
    BObolViewController *controller, void *userdata)
{
    std::set<SoBRLDatabaseSource *> *sources =
	static_cast<std::set<SoBRLDatabaseSource *> *>(userdata);
    if (controller && sources)
	ged_obol_collect_database_sources(controller->getRenderSceneRoot(),
	    *sources);
    return 1;
}

static bool
ged_obol_path_in_selected_set(const char *candidate,
	const char *const *paths, size_t path_count)
{
    if (!candidate || !candidate[0])
	return false;
    for (size_t i = 0; i < path_count; i++) {
	if (paths[i] && paths[i][0] &&
	    (ged_obol_path_has_prefix(candidate, paths[i]) ||
	     ged_obol_path_has_semantic_prefix(candidate, paths[i])))
	    return true;
    }
    return false;
}

extern "C" int
ged_draw_obol_database_sources_sync_selected_paths(
    struct ged *gedp,
    const char *const *paths,
    size_t path_count)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp)
	return 0;

    std::set<SoBRLDatabaseSource *> sources;
    BObolSceneController *owned = ged_draw_obol_scene_controller(gedp);
    if (owned) {
	for (int i = 0; i < owned->getDatabaseSourceCount(); i++)
	    sources.insert(owned->getDatabaseSource(i));
    }

    ged_bobol_view_controllers_foreach(gedp,
	ged_obol_collect_view_database_sources, &sources);
    int applied = 0;
    std::vector<SbString> compact_selected_paths;
    compact_selected_paths.reserve(path_count);
    for (size_t i = 0; i < path_count; i++) {
	if (paths[i] && paths[i][0])
	    compact_selected_paths.push_back(SbString(paths[i]));
    }
    for (SoBRLDatabaseSource *source : sources) {
	BObolDatabaseSourceSummary sourceSummary;
	if (!source || !source->getSummary(sourceSummary) ||
	    !sourceSummary.valid)
	    continue;
	/* Retain the selection rule before deferred compact geometry arrives.
	 * Index installation and streamed additions will then apply it directly,
	 * rather than performing a delayed full-scene catch-up during the next
	 * unrelated draw or erase operation. */
	applied += source->syncCompactInstanceSelectedPaths(
	    compact_selected_paths);

	const SbBool sourceSelected = ged_obol_path_in_selected_set(
	    sourceSummary.path.getString(), paths, path_count) ? TRUE : FALSE;
	if (sourceSummary.selected == sourceSelected)
	    continue;
	const int changed = source->setDisplayState(FALSE,
	    sourceSummary.sourceRevision, sourceSummary.inputsRevision,
	    sourceSummary.visible, sourceSelected, sourceSummary.highlighted,
	    sourceSummary.lineStyle, sourceSummary.lineWidth,
	    sourceSummary.transparency, sourceSummary.colorOverride,
	    sourceSummary.color, sourceSummary.materialColorValid,
	    sourceSummary.materialColor, sourceSummary.materialRevision);
	if (changed > 0)
	    applied++;
    }
    return applied;
}

template <typename ShapeT>
static int
ged_obol_shape_update_display_typed(
    BObolSceneController *scene,
    const char *path,
    ShapeT *shape,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision)
{
    if (!scene || !path || !path[0] || !shape)
	return 0;

    if (draw_mode_valid) {
	const int draw_changed = scene->setShapeDrawState(path,
				 ged_obol_lod_draw_mode_from_ged(draw_mode),
				 shape->databaseIntent.getValue(),
				 shape->overlayIntent.getValue(),
				 shape->hudIntent.getValue());
	if (draw_changed < 0)
	    return 0;
    }

    const SbBool next_visible = visible_valid ?
				(visible ? TRUE : FALSE) : shape->visible.getValue();
    const SbBool next_selected = selected_valid ?
				 (selected ? TRUE : FALSE) : shape->selected.getValue();
    const SbBool next_highlighted = highlighted_valid ?
				    (highlighted ? TRUE : FALSE) : shape->highlighted.getValue();
    const int next_line_style = line_style_valid ?
				line_style : shape->lineStyle.getValue();
    const int next_line_width = line_width_valid ?
				line_width : shape->lineWidth.getValue();
    const float next_transparency = transparency_valid ?
				    static_cast<float>(transparency) : shape->transparency.getValue();

    SbBool next_color_override = shape->colorOverride.getValue();
    SbColor next_color = shape->color.getValue();
    if (color_valid) {
	next_color_override = TRUE;
	next_color = ged_obol_color_from_rgb(color);
    }

    SbBool next_material_color_valid = shape->materialColorValid.getValue();
    SbColor next_material_color = shape->materialColor.getValue();
    uint32_t next_material_revision = shape->materialRevision.getValue();
    if (material_color_valid) {
	next_material_color_valid = TRUE;
	next_material_color = ged_obol_color_from_rgb(material_color);
	if (material_revision_valid) {
	    next_material_revision =
		ged_obol_fold_revision(material_revision);
	} else {
	    next_material_revision++;
	    if (!next_material_revision)
		next_material_revision = 1;
	}
    }

    const int changed = scene->setShapeDisplayState(path,
			next_visible,
			next_selected,
			next_highlighted,
			next_line_style,
			next_line_width,
			next_transparency,
			next_color_override,
			next_color,
			next_material_color_valid,
			next_material_color,
			next_material_revision);
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_shape_update_display_for_path(
    struct ged *gedp,
    const char *path,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int draw_mode,
    int line_style_valid,
    int line_style,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_valid,
    const unsigned char color[3],
    int material_color_valid,
    const unsigned char material_color[3],
    int material_revision_valid,
    uint64_t material_revision)
{
    if ((color_valid && !color) ||
	(material_color_valid && !material_color))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node)
	return 0;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	return ged_obol_shape_update_display_typed(scene, path,
		static_cast<SoBRLVListShape *>(node),
		visible_valid, visible,
		selected_valid, selected,
		highlighted_valid, highlighted,
		draw_mode_valid, draw_mode,
		line_style_valid, line_style,
		line_width_valid, line_width,
		transparency_valid, transparency,
		color_valid, color,
		material_color_valid, material_color,
		material_revision_valid, material_revision);
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	return ged_obol_shape_update_display_typed(scene, path,
		static_cast<SoBRLMeshShape *>(node),
		visible_valid, visible,
		selected_valid, selected,
		highlighted_valid, highlighted,
		draw_mode_valid, draw_mode,
		line_style_valid, line_style,
		line_width_valid, line_width,
		transparency_valid, transparency,
		color_valid, color,
		material_color_valid, material_color,
		material_revision_valid, material_revision);
    }

    return 0;
}

extern "C" int
ged_draw_obol_database_source_update_appearance_for_path(
    struct ged *gedp,
    const char *path,
    int line_width_valid,
    int line_width,
    int transparency_valid,
    fastf_t transparency,
    int color_override_valid,
    int color_override,
    int color_valid,
    const unsigned char color[3])
{
    if (color_valid && !color)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    BObolDatabaseSourceDisplayPatch patch;
    patch.lineWidthValid = line_width_valid ? TRUE : FALSE;
    patch.lineWidth = line_width;
    patch.transparencyValid = transparency_valid ? TRUE : FALSE;
    patch.transparency = static_cast<float>(transparency);
    patch.colorOverrideValid = color_override_valid ? TRUE : FALSE;
    patch.colorOverride = color_override ? TRUE : FALSE;
    patch.colorValid = color_valid ? TRUE : FALSE;
    if (color_valid)
	patch.color = ged_obol_color_from_rgb(color);

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	int changed = scene->setDatabaseSourceInstanceDisplayPatch(
			  source_instance_key.c_str(), patch);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_group_update_display_for_path(
    struct ged *gedp,
    const char *path,
    int visible_valid,
    int visible)
{
    if (!gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    int changed = scene->setGroupDisplayState(group_path.c_str(),
		  visible_valid ? (visible ? TRUE : FALSE) :
		  scene_group->visible.getValue(),
		  scene_group->selected.getValue(),
		  scene_group->highlighted.getValue(),
		  scene_group->lineStyle.getValue(),
		  scene_group->lineWidth.getValue(),
		  scene_group->transparency.getValue(),
		  scene_group->colorOverride.getValue(),
		  scene_group->color.getValue(),
		  scene_group->materialColorValid.getValue(),
		  scene_group->materialColor.getValue(),
		  scene_group->materialRevision.getValue());
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_ensure_for_path(
    struct ged *gedp,
    const char *path,
    const char *intent_path,
    int mode,
    int overlay)
{
    if (!gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    if (group_path.empty())
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoGroup *group = scene->ensureGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    const int next_draw_mode = mode >= 0 ?
			       ged_obol_lod_draw_mode_from_ged(mode) :
			       scene_group->drawMode.getValue();
    int fallback_draw_mode = scene_group->fallbackDrawMode.getValue();
    if (fallback_draw_mode == BOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BOBOL_LOD_DRAW_WIRE;
    const SbBool next_overlay = overlay >= 0 ?
				(overlay ? TRUE : FALSE) :
				scene_group->overlayIntent.getValue();

    const std::string target_path =
	(intent_path && intent_path[0]) ?
	ged_obol_group_path_from_record_path(intent_path) : group_path;
    const std::string draw_intent_path =
	ged_obol_group_intent_path(
	    target_path.empty() ? group_path.c_str() :
	    target_path.c_str());

    const int changed = scene->setGroupDrawIntent(group_path.c_str(),
			draw_intent_path.c_str(),
			next_draw_mode,
			fallback_draw_mode,
			next_overlay,
			scene_group->revalidationRevision.getValue());
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_record_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_group_record_summary *out)
{
    if (!out || !gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    const char *record_path = ged_obol_group_record_path(scene_group);
    if (record_path && record_path[0])
	out->path = record_path;
    out->draw_mode = ged_obol_lod_draw_mode_to_ged(
			 scene_group->drawMode.getValue());
    out->transparency = scene_group->transparency.getValue();
    out->visible = scene_group->visible.getValue() ? 1 : 0;
    out->is_overlay = scene_group->overlayIntent.getValue() ? 1 : 0;
    return 1;
}


static int
ged_draw_obol_group_paths_foreach_node(
    struct ged *gedp,
    BObolSceneController *scene,
    SoNode *node,
    int skip_overlay_groups,
    ged_draw_obol_group_path_cb cb,
    void *userdata)
{
    if (!node)
	return 1;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return 1;

    if (node->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	SoBRLSceneGroup *group = static_cast<SoBRLSceneGroup *>(node);
	const char *group_path = group->groupPath.getValue().getString();
	if ((!skip_overlay_groups ||
	     !ged_obol_group_path_is_overlay(scene, group_path)) &&
	    group_path && group_path[0] &&
	    !BU_STR_EQUAL(group_path, "/")) {
	    const char *record_path = ged_obol_group_record_path(group);
	    if (record_path && record_path[0] &&
		!(*cb)(gedp, record_path, userdata))
		return 0;
	}
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return 1;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (!ged_draw_obol_group_paths_foreach_node(gedp, scene,
		group->getChild(i), skip_overlay_groups, cb, userdata))
	    return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_group_paths_foreach(
    struct ged *gedp,
    int skip_overlay_groups,
    ged_draw_obol_group_path_cb cb,
    void *userdata)
{
    if (!gedp || !cb)
	return -1;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    return ged_draw_obol_group_paths_foreach_node(gedp, scene,
	scene->getSceneRoot(), skip_overlay_groups, cb, userdata);
}


extern "C" int
ged_draw_obol_group_display_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_scene_display_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    int source_like_local_group = 0;
    for (int i = 0; i < scene_group->getNumChildren(); i++) {
	SoNode *child = scene_group->getChild(i);
	if (child &&
	    (child->isOfType(SoBRLVListShape::getClassTypeId()) ||
	     child->isOfType(SoBRLMeshShape::getClassTypeId())) &&
	    ged_obol_shape_node_bool_field(child, "databaseIntent", FALSE) &&
	    ged_obol_shape_node_bool_field(child, "localSource", FALSE)) {
	    source_like_local_group = 1;
	    break;
	}
    }

    out->valid = 1;
    out->is_database_source = source_like_local_group;
    out->has_draw_intent =
	scene_group->drawIntentValid.getValue() ? 1 : 0;
    out->intent_path = scene_group->drawIntentPath.getValue().getString();
    out->intent_draw_mode = ged_obol_lod_draw_mode_to_ged(
				scene_group->drawMode.getValue());
    out->visible = scene_group->visible.getValue() ? 1 : 0;
    out->selected = scene_group->selected.getValue() ? 1 : 0;
    out->highlighted = scene_group->highlighted.getValue() ? 1 : 0;
    out->line_style = scene_group->lineStyle.getValue();
    out->line_width = scene_group->lineWidth.getValue();
    out->transparency = scene_group->transparency.getValue();
    out->draw_mode = ged_obol_lod_draw_mode_to_ged(
			 scene_group->drawMode.getValue());
    out->material_valid =
	(scene_group->materialColorValid.getValue() ||
	 scene_group->colorOverride.getValue()) ? 1 : 0;
    if (scene_group->materialColorValid.getValue())
	ged_obol_rgb_from_color(scene_group->materialColor.getValue(),
				out->material_color);
    else if (scene_group->colorOverride.getValue())
	ged_obol_rgb_from_color(scene_group->color.getValue(),
				out->material_color);
    return 1;
}

extern "C" int
ged_draw_obol_group_shape_count_for_path(
    struct ged *gedp,
    const char *path,
    int *out)
{
    if (out)
	*out = 0;
    if (!out || !gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    int count = scene->getGroupDatabaseSourceCount(group_path.c_str());
    if (count < 0)
	count = 0;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	const std::string owner_group_path =
	    ged_obol_database_source_owner_group_path_from_summary(summary);
	if (ged_obol_path_equal(owner_group_path.c_str(),
				group_path.c_str()))
	    continue;
	if (ged_obol_path_has_prefix(owner_group_path.c_str(),
				     group_path.c_str()))
	    count++;
    }

    *out = count;
    return 1;
}


extern "C" int
ged_draw_obol_group_descendant_group_count_for_path(
    struct ged *gedp,
    const char *path,
    size_t *out)
{
    if (out)
	*out = 0;
    if (!out || !gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    const int count = scene->getGroupDescendantGroupCount(
			  group_path.c_str());
    if (count < 0)
	return 0;

    *out = static_cast<size_t>(count);
    return 1;
}


extern "C" int
ged_draw_obol_group_child_count_for_path(
    struct ged *gedp,
    const char *path,
    size_t *out)
{
    if (out)
	*out = 0;
    if (!out || !gedp || !path || !path[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    const int count = scene->getGroupChildCount(group_path.c_str());
    if (count < 0)
	return 0;

    *out = static_cast<size_t>(count);
    return 1;
}


extern "C" int
ged_draw_obol_group_update_appearance_for_path(
    struct ged *gedp,
    const char *path,
    const struct ged_draw_appearance_settings *settings)
{
    if (!gedp || !path || !path[0] || !settings)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    const char *retained_intent =
	scene_group->drawIntentPath.getValue().getString();
    const std::string next_intent_path =
	(retained_intent && retained_intent[0]) ?
	std::string(retained_intent) :
	ged_obol_group_intent_path(group_path.c_str());
    int fallback_draw_mode = scene_group->fallbackDrawMode.getValue();
    if (fallback_draw_mode == BOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BOBOL_LOD_DRAW_WIRE;
    int intent_changed = scene->setGroupDrawIntent(group_path.c_str(),
			 next_intent_path.c_str(),
			 ged_obol_lod_draw_mode_from_ged(settings->draw_mode),
			 fallback_draw_mode,
			 scene_group->overlayIntent.getValue(),
			 scene_group->revalidationRevision.getValue());
    if (intent_changed < 0)
	return 0;

    const SbColor next_color = ged_obol_color_from_rgb(settings->color);
    int changed = scene->setGroupDisplayState(group_path.c_str(),
		  scene_group->visible.getValue(),
		  scene_group->selected.getValue(),
		  scene_group->highlighted.getValue(),
	  scene_group->lineStyle.getValue(),
	  settings->s_line_width,
	  ged_obol_transparency_from_appearance_opacity(settings->transparency),
		  settings->color_override ? TRUE : FALSE,
		  next_color,
		  scene_group->materialColorValid.getValue(),
		  scene_group->materialColor.getValue(),
		  scene_group->materialRevision.getValue());
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_appearance_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_appearance_settings *settings)
{
    if (!gedp || !path || !path[0] || !settings)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    struct ged_draw_appearance_settings next =
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
    next.draw_mode = ged_obol_lod_draw_mode_to_ged(
			 scene_group->drawMode.getValue());
    next.transparency = ged_obol_appearance_opacity_from_transparency(
	scene_group->transparency.getValue());
    next.color_override = scene_group->colorOverride.getValue() ? 1 : 0;
    ged_obol_rgb_from_color(scene_group->color.getValue(), next.color);
    next.s_line_width = scene_group->lineWidth.getValue();
    *settings = next;
    return 1;
}

extern "C" int
ged_draw_obol_group_update_draw_intent_for_path(
    struct ged *gedp,
    const char *path,
    const char *intent_path,
    int mode_valid,
    int mode,
    int overlay_valid,
    int overlay)
{
    if (!gedp)
	return 0;

    std::string group_path = ged_obol_group_path_from_record_path(path);
    std::string target_group_path =
	(intent_path && intent_path[0]) ?
	ged_obol_group_path_from_record_path(intent_path) : group_path;
    if (group_path.empty())
	group_path = target_group_path;
    if (target_group_path.empty())
	target_group_path = group_path;
    if (group_path.empty())
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoGroup *group = scene->findGroup(group_path.c_str());
    if (!group && target_group_path != group_path) {
	group = scene->findGroup(target_group_path.c_str());
	if (group)
	    group_path = target_group_path;
    }
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;

    if (target_group_path != group_path) {
	std::string old_parent, old_leaf, new_parent, new_leaf;
	if (ged_obol_group_parent_leaf(group_path, old_parent, old_leaf) &&
	    ged_obol_group_parent_leaf(target_group_path, new_parent,
				       new_leaf) &&
	    old_parent == new_parent &&
	    scene->renameGroup(group_path.c_str(), new_leaf.c_str()) > 0) {
	    group_path = target_group_path;
	    group = scene->findGroup(group_path.c_str());
	}
	if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    return 0;
    }

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    int next_draw_mode = mode_valid ?
			 ged_obol_lod_draw_mode_from_ged(mode) :
			 scene_group->drawMode.getValue();
    int fallback_draw_mode = scene_group->fallbackDrawMode.getValue();
    if (fallback_draw_mode == BOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BOBOL_LOD_DRAW_WIRE;
    SbBool next_overlay = overlay_valid ?
			  (overlay ? TRUE : FALSE) : scene_group->overlayIntent.getValue();
    uint32_t revalidation_revision =
	scene_group->revalidationRevision.getValue();

    std::string next_intent_path;
    if (intent_path && intent_path[0]) {
	next_intent_path =
	    ged_obol_group_intent_path(target_group_path.c_str());
    } else {
	const char *retained_intent =
	    scene_group->drawIntentPath.getValue().getString();
	if (retained_intent && retained_intent[0])
	    next_intent_path = retained_intent;
	else
	    next_intent_path = ged_obol_group_intent_path(group_path.c_str());
    }

    int changed = scene->setGroupDrawIntent(group_path.c_str(),
					    next_intent_path.c_str(),
					    next_draw_mode,
					    fallback_draw_mode,
					    next_overlay,
					    revalidation_revision);
    return changed >= 0 ? 1 : 0;
}

static int
ged_obol_database_source_exact_draw_mode_to_ged(
    struct ged *gedp,
    const BObolDatabaseSourceSummary &summary,
    SoBRLDatabaseSource *source)
{
    if (summary.representationMode >= 0)
	return summary.representationMode;

    const int source_ged_mode =
	ged_obol_database_draw_mode_to_ged(summary.drawMode);
    int exact_ged_mode = source_ged_mode;
    if (source_ged_mode == GED_DRAW_MODE_SHADED && source &&
	source->hasRealizedMeshGeometry())
	exact_ged_mode = GED_DRAW_MODE_SHADED_BOTS;

    const std::string owner_group_path =
	ged_obol_database_source_owner_group_path_from_summary(summary);
    if (owner_group_path.empty())
	return exact_ged_mode;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    SoGroup *group = scene ? scene->findGroup(owner_group_path.c_str()) :
		     NULL;
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return exact_ged_mode;

    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    const int group_ged_mode =
	ged_obol_lod_draw_mode_to_ged(scene_group->drawMode.getValue());
    if (group_ged_mode != GED_DRAW_MODE_WIRE ||
	source_ged_mode == GED_DRAW_MODE_WIRE)
	return group_ged_mode;
    return exact_ged_mode;
}

extern "C" int
ged_draw_obol_database_source_display_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_scene_display_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!ged_obol_database_source_controller_summary_for_source(scene,
	    source, summary))
	return 0;
    if (!summary.valid)
	return 0;

    int exact_draw_mode =
	ged_obol_database_source_exact_draw_mode_to_ged(gedp, summary,
	    source);

    out->valid = 1;
    out->is_database_source = 1;
    out->has_draw_intent = 1;
    out->intent_path = source->path.getValue().getString();
    out->intent_draw_mode = exact_draw_mode;
    out->visible = summary.visible ? 1 : 0;
    out->selected = summary.selected ? 1 : 0;
    out->highlighted = summary.highlighted ? 1 : 0;
    out->line_style = summary.lineStyle;
    out->line_width = summary.lineWidth;
    out->transparency = ged_obol_reported_transparency(summary.transparency);
    out->draw_mode = exact_draw_mode;
    out->material_valid = (summary.materialColorValid ||
			   summary.databaseMaterialColorValid ||
			   summary.colorOverride) ? 1 : 0;
    if (out->material_valid)
	ged_obol_rgb_from_color(ged_obol_summary_material_color(summary),
				out->material_color);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_draw_state_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_obol_draw_state_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    SoBRLDatabaseSource *source = NULL;
    BObolDatabaseSourceSummary source_summary;
    int have_source = 0;
    int best_score = -1;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *candidate =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	BObolDatabaseSourceSummary candidate_summary;
	if (!ged_obol_database_source_controller_summary_for_source(scene,
		candidate, candidate_summary) || !candidate_summary.valid)
	    continue;

	int score = 0;
	if (candidate_summary.drawMatrixValid)
	    score += 4;
	if (candidate_summary.lineStyle)
	    score += 2;
	if (candidate_summary.lineWidth)
	    score += 1;
	if (!have_source || score > best_score) {
	    source = candidate;
	    source_summary = candidate_summary;
	    have_source = 1;
	    best_score = score;
	}
    }
    if (!have_source || !source)
	return 0;

    out->valid = 1;
    out->draw_mode_valid = 1;
    out->draw_mode = ged_obol_database_source_exact_draw_mode_to_ged(gedp,
		     source_summary, source);
    out->line_style = source_summary.lineStyle;
    if (source_summary.drawMatrixValid) {
	out->draw_mat_valid = 1;
	ged_obol_mat_from_sbmatrix(source_summary.drawMatrix, out->draw_mat);
    }

    if (!out->draw_mat_valid) {
	const int count = source->getRealizedDisplaySummaryCount();
	for (int i = 0; i < count; i++) {
	    BObolSceneDisplaySummary display_summary;
	    if (!source->getRealizedDisplaySummary(i, display_summary) ||
		!display_summary.valid ||
		!display_summary.drawMatrixValid)
		continue;
	    if (display_summary.nodeKind !=
		BObolSceneTreeSummary::NODE_VLIST_SHAPE &&
		display_summary.nodeKind !=
		BObolSceneTreeSummary::NODE_MESH_SHAPE)
		continue;
	    out->draw_mat_valid = 1;
	    ged_obol_mat_from_sbmatrix(display_summary.drawMatrix,
				       out->draw_mat);
	    break;
	}
    }

    const auto matrix_delta = [](const mat_t mat) -> fastf_t {
	fastf_t delta = 0.0;
	for (int i = 0; i < 16; i++) {
	    const fastf_t expected =
		(i == 0 || i == 5 || i == 10 || i == 15) ? 1.0 : 0.0;
	    delta += fabs(mat[i] - expected);
	}
	return delta;
    };

    fastf_t best_matrix_delta =
	out->draw_mat_valid ? matrix_delta(out->draw_mat) : -1.0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *candidate =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	BObolDatabaseSourceSummary candidate_summary;
	if (!ged_obol_database_source_controller_summary_for_source(scene,
		candidate, candidate_summary) || !candidate_summary.valid)
	    continue;
	if (candidate_summary.lineStyle > out->line_style)
	    out->line_style = candidate_summary.lineStyle;
	if (candidate_summary.drawMatrixValid) {
	    mat_t candidate_mat;
	    ged_obol_mat_from_sbmatrix(candidate_summary.drawMatrix,
				       candidate_mat);
	    const fastf_t candidate_delta = matrix_delta(candidate_mat);
	    if (!out->draw_mat_valid ||
		candidate_delta > best_matrix_delta + VUNITIZE_TOL) {
		out->draw_mat_valid = 1;
		MAT_COPY(out->draw_mat, candidate_mat);
		best_matrix_delta = candidate_delta;
	    }
	}
    }

    return 1;
}

static const char *ged_obol_leaf_name_from_path(const char *path);

static SoBRLVListShape *
ged_obol_owned_vlist_shape_for_source(SoBRLDatabaseSource *source,
				      const char *fallback_path,
				      int create)
{
    if (!source)
	return NULL;

    SoBRLVListShape *fallback = NULL;
    const int count = source->getRealizedShapeCount();
    for (int i = 0; i < count; i++) {
	SoBRLVListShape *shape = source->getRealizedShape(i);
	if (!shape)
	    continue;
	if (!fallback)
	    fallback = shape;
	const SoBRLVListShape *geom = shape->getGeometrySource();
	if (geom->point.getNum() > 0 || geom->command.getNum() > 0)
	    return shape;
    }

    if (fallback || !create)
	return fallback;

    const char *source_path = source->path.getValue().getString();
    if (!source_path || !source_path[0])
	source_path = fallback_path ? fallback_path : "";
    const char *source_name = ged_obol_leaf_name_from_path(source_path);

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = source_path;
    shape->sourceName = source_name;
    shape->sourceType = "line-set";
    shape->sourceId = source->realizedRevision.getValue();
    shape->displayName = source_name;
    shape->geometryName = source_name;
    shape->sourceIdentity = source_path;
    shape->cacheIdentity = source_path;
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = FALSE;
    shape->drawMode = ged_obol_lod_draw_mode_from_database_source(source);
    shape->recordRole = "database";
    shape->geometryKind = "line-set";
    shape->visible = source->visible.getValue();
    shape->highlighted = source->highlighted.getValue();
    shape->lineStyle = source->lineStyle.getValue();
    shape->lineWidth = source->lineWidth.getValue();
    shape->transparency = source->transparency.getValue();
    shape->hiddenLine = shape->drawMode.getValue() ==
			BOBOL_LOD_DRAW_HIDDEN_LINE ? TRUE : FALSE;
    shape->colorOverride = source->colorOverride.getValue();
    shape->color = source->color.getValue();
    shape->materialColorValid = source->materialColorValid.getValue();
    shape->materialColor = source->materialColor.getValue();
    shape->materialRevision = source->materialRevision.getValue();
    shape->drawMatrixValid = source->drawMatrixValid.getValue();
    shape->drawMatrix = source->drawMatrix.getValue();
    shape->drawCenterValid = source->drawCenterValid.getValue();
    shape->drawCenter = source->drawCenter.getValue();
    shape->drawSizeValid = source->drawSizeValid.getValue();
    shape->drawSize = source->drawSize.getValue();
    source->addChild(shape);
    return shape;
}

static SoBRLVListShape *
ged_obol_owned_vlist_shape_for_path(struct ged *gedp, const char *path)
{
    return ged_obol_owned_vlist_shape_for_source(
	       ged_obol_owned_database_source_for_path(gedp, path), path, 0);
}


static int
ged_obol_vlist_shape_is_annotation(SoBRLVListShape *shape)
{
    if (!shape)
	return 0;

    const char *kind = shape->geometryKind.getValue().getString();
    if (kind && BU_STR_EQUAL(kind, "annotation"))
	return 1;

    const char *source_type = shape->sourceType.getValue().getString();
    return (source_type && BU_STR_EQUAL(source_type, "annotation")) ? 1 : 0;
}


static int
ged_obol_vlist_shape_has_annotation_record(SoBRLVListShape *shape)
{
    if (!shape)
	return 0;
    const SoBRLVListShape *geom = shape->getGeometrySource();
    return (geom->annotationPoint.getNum() > 0 ||
	    geom->annotationSegmentKind.getNum() > 0) ? 1 : 0;
}


static SoBRLVListShape *
ged_obol_owned_annotation_vlist_shape_for_source(
    SoBRLDatabaseSource *source,
    const char *fallback_path)
{
    if (!source)
	return NULL;

    SoBRLVListShape *annotation_shape = NULL;
    const int count = source->getRealizedShapeCount();
    for (int i = 0; i < count; i++) {
	SoBRLVListShape *shape = source->getRealizedShape(i);
	if (!ged_obol_vlist_shape_is_annotation(shape))
	    continue;
	if (ged_obol_vlist_shape_has_annotation_record(shape))
	    return shape;
	if (!annotation_shape)
	    annotation_shape = shape;
    }

    return annotation_shape ? annotation_shape :
	   ged_obol_owned_vlist_shape_for_source(source, fallback_path, 0);
}


static SoBRLVListShape *
ged_obol_owned_annotation_vlist_shape_for_path(struct ged *gedp,
	const char *path)
{
    return ged_obol_owned_annotation_vlist_shape_for_source(
	       ged_obol_owned_database_source_for_path(gedp, path), path);
}

static const char *
ged_obol_leaf_name_from_path(const char *path)
{
    if (!path || !path[0])
	return "";
    const char *leaf = strrchr(path, '/');
    return (leaf && leaf[1]) ? leaf + 1 : path;
}

static SoBRLMeshShape *
ged_obol_owned_mesh_shape_for_source(SoBRLDatabaseSource *source,
				     const char *fallback_path,
				     int create)
{
    if (!source)
	return NULL;

    SoBRLMeshShape *shape = source->getRealizedMesh();
    if (shape || !create)
	return shape;

    const char *source_path = source->path.getValue().getString();
    if (!source_path || !source_path[0])
	source_path = fallback_path ? fallback_path : "";
    const char *source_name = ged_obol_leaf_name_from_path(source_path);

    shape = new SoBRLMeshShape;
    shape->sourcePath = source_path;
    shape->sourceName = source_name;
    shape->sourceType = "indexed-face-set";
    shape->sourceId = source->realizedRevision.getValue();
    shape->displayName = source_name;
    shape->geometryName = source_name;
    shape->sourceIdentity = source_path;
    shape->cacheIdentity = source_path;
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = FALSE;
    shape->drawMode = ged_obol_lod_draw_mode_from_database_source(source);
    shape->recordRole = "database";
    shape->geometryKind = "surface";
    shape->visible = source->visible.getValue();
    shape->highlighted = source->highlighted.getValue();
    shape->lineStyle = source->lineStyle.getValue();
    shape->lineWidth = source->lineWidth.getValue();
    shape->transparency = source->transparency.getValue();
    shape->hiddenLine = shape->drawMode.getValue() ==
			BOBOL_LOD_DRAW_HIDDEN_LINE ? TRUE : FALSE;
    shape->colorOverride = source->colorOverride.getValue();
    shape->color = source->color.getValue();
    shape->materialColorValid = source->materialColorValid.getValue();
    shape->materialColor = source->materialColor.getValue();
    shape->materialRevision = source->materialRevision.getValue();
    shape->drawMatrixValid = source->drawMatrixValid.getValue();
    shape->drawMatrix = source->drawMatrix.getValue();
    shape->drawCenterValid = source->drawCenterValid.getValue();
    shape->drawCenter = source->drawCenter.getValue();
    shape->drawSizeValid = source->drawSizeValid.getValue();
    shape->drawSize = source->drawSize.getValue();
    source->addChild(shape);
    return shape;
}

static SoBRLMeshShape *
ged_obol_owned_mesh_shape_for_path(struct ged *gedp, const char *path,
				   int create)
{
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    return ged_obol_owned_mesh_shape_for_source(source, path, create);
}

static uint64_t
ged_obol_hash_sb_string(const SbString &value)
{
    const char *str = value.getString();
    if (!str || !str[0])
	return 0;
    return bu_data_hash(str, strlen(str) * sizeof(char));
}

static uint64_t
ged_obol_hash_cstr(const char *str)
{
    if (!str || !str[0])
	return 0;
    return bu_data_hash(str, strlen(str) * sizeof(char));
}

static int
ged_obol_vlist_command_to_ged(int command)
{
    if (command == SoBRLVListShape::MOVE)
	return GED_DRAW_VIEW_LINE_MOVE;
    if (command == SoBRLVListShape::DRAW)
	return GED_DRAW_VIEW_LINE_DRAW;
    if (command == SoBRLVListShape::POINT)
	return GED_DRAW_VIEW_LINE_POINT_DRAW;
    return command;
}

static int32_t
ged_obol_vlist_command_from_ged(int command, size_t index)
{
    if (command == GED_DRAW_VIEW_LINE_MOVE)
	return SoBRLVListShape::MOVE;
    if (command == GED_DRAW_VIEW_LINE_DRAW)
	return SoBRLVListShape::DRAW;
    if (command == GED_DRAW_VIEW_LINE_POINT_DRAW)
	return SoBRLVListShape::POINT;
    if (command < 0 && (index % 2) == 0)
	return SoBRLVListShape::MOVE;
    if (command < 0)
	return SoBRLVListShape::DRAW;
    return -1;
}

static void
ged_obol_vlist_shape_set_precise_points(SoBRLVListShape *shape,
					const point_t *points,
					size_t point_count)
{
    if (!shape || !points || point_count > static_cast<size_t>(INT_MAX))
	return;

    std::vector<double> precise_points;
    precise_points.reserve(point_count * 3);
    for (size_t i = 0; i < point_count; i++) {
	precise_points.push_back(points[i][0]);
	precise_points.push_back(points[i][1]);
	precise_points.push_back(points[i][2]);
    }

    shape->setPrecisePoints(precise_points.empty() ? NULL :
			    precise_points.data(), static_cast<int>(point_count));
}

template <typename ShapeT>
static void
ged_obol_local_shape_apply_common_state(
    ShapeT *shape,
    const char *shape_path,
    const char *display_name,
    const char *source_type,
    const char *geometry_kind,
    const struct ged_draw_scene_display_summary *display_state)
{
    if (!shape || !shape_path)
	return;

    const char *leaf = ged_obol_leaf_name_from_path(shape_path);
    const char *name = (display_name && display_name[0]) ? display_name : leaf;
    const int local_draw_mode = (display_state && display_state->valid) ?
				display_state->draw_mode : GED_DRAW_MODE_WIRE;

    shape->sourcePath = shape_path;
    shape->sourceName = leaf;
    shape->sourceType = source_type ? source_type : "local";
    shape->sourceId = static_cast<uint32_t>(ged_obol_hash_cstr(shape_path));
    shape->displayName = name ? name : leaf;
    shape->geometryName = name ? name : leaf;
    shape->sourceIdentity = shape_path;
    shape->cacheIdentity = shape_path;
    const char *intent_path = (display_state && display_state->valid &&
			       display_state->intent_path && display_state->intent_path[0]) ?
			      display_state->intent_path : "";
    shape->ownerSourcePath = intent_path;
    shape->databaseIntent = intent_path[0] ? TRUE : FALSE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = TRUE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = ged_obol_lod_draw_mode_from_ged(local_draw_mode);
    shape->recordRole = "local";
    shape->geometryKind = geometry_kind ? geometry_kind : "local";
    shape->visible = (!display_state || !display_state->valid ||
		      display_state->visible) ? TRUE : FALSE;
    shape->highlighted = (display_state && display_state->valid &&
			  display_state->highlighted) ? TRUE : FALSE;
    shape->lineStyle = (display_state && display_state->valid) ?
		       display_state->line_style : 0;
    shape->lineWidth = (display_state && display_state->valid) ?
		       display_state->line_width : 0;
    shape->transparency = (display_state && display_state->valid) ?
			  static_cast<float>(display_state->transparency) : 0.0f;
    shape->colorOverride = FALSE;
    shape->color = SbColor(1.0f, 1.0f, 1.0f);
    shape->materialColorValid = (display_state && display_state->valid &&
				 display_state->material_valid) ? TRUE : FALSE;
    if (display_state && display_state->valid &&
	display_state->material_valid) {
	const SbColor material =
	    ged_obol_color_from_rgb(display_state->material_color);
	shape->materialColor = material;
	shape->color = material;
	shape->colorOverride = TRUE;
    }
    shape->materialRevision = 0;
}

static SoBRLVListShape *
ged_obol_local_vlist_shape_for_path(
    BObolSceneController *scene,
    const char *group_path,
    const char *shape_path)
{
    if (!scene || !group_path || !group_path[0] ||
	!shape_path || !shape_path[0])
	return NULL;

    if (!scene->ensureGroup(group_path))
	return NULL;

    SoNode *node = scene->findShape(shape_path);
    if (node && !node->isOfType(SoBRLVListShape::getClassTypeId())) {
	(void)scene->removeShape(shape_path);
	node = NULL;
    }
    if (node) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	(void)scene->moveShapeToGroup(shape_path, group_path);
	return shape;
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->sourcePath = shape_path;
    const int appended = scene->appendChildToGroup(group_path, shape);
    if (appended < 0)
	return NULL;
    return shape;
}

static SoBRLMeshShape *
ged_obol_local_mesh_shape_for_path(
    BObolSceneController *scene,
    const char *group_path,
    const char *shape_path)
{
    if (!scene || !group_path || !group_path[0] ||
	!shape_path || !shape_path[0])
	return NULL;

    if (!scene->ensureGroup(group_path))
	return NULL;

    SoNode *node = scene->findShape(shape_path);
    if (node && !node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	(void)scene->removeShape(shape_path);
	node = NULL;
    }
    if (node) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	(void)scene->moveShapeToGroup(shape_path, group_path);
	return shape;
    }

    SoBRLMeshShape *shape = new SoBRLMeshShape;
    shape->sourcePath = shape_path;
    const int appended = scene->appendChildToGroup(group_path, shape);
    if (appended < 0)
	return NULL;
    return shape;
}

extern "C" int
ged_draw_obol_local_shape_publish_line_set_for_path(
    struct ged *gedp,
    const char *group_path,
    const char *shape_path,
    const char *display_name,
    const point_t *points,
    const int *commands,
    size_t point_count,
    const struct ged_draw_scene_display_summary *display_state)
{
    if (!group_path || !group_path[0] || !shape_path || !shape_path[0] ||
	point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_local_vlist_shape_for_path(scene, group_path, shape_path);
    if (!shape)
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	const int command = commands ? commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	obol_commands.push_back(obol_command);
    }

    ged_obol_local_shape_apply_common_state(shape, shape_path, display_name,
					    "local-line-set", "line", display_state);
    shape->setLineSet(obol_points.empty() ? NULL : obol_points.data(),
		      obol_commands.empty() ? NULL : obol_commands.data(),
		      static_cast<int>(point_count));
    ged_obol_vlist_shape_set_precise_points(shape, points, point_count);
    return 1;
}

extern "C" int
ged_draw_obol_local_shape_set_record_role_for_path(
    struct ged *gedp,
    const char *shape_path,
    const char *record_role)
{
    if (!shape_path || !shape_path[0] || !record_role)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoNode *node = scene->findShape(shape_path);
    if (!node)
	return 0;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	shape->recordRole = record_role;
	return 1;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	shape->recordRole = record_role;
	return 1;
    }

    return 0;
}

extern "C" int
ged_draw_obol_shape_display_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_scene_display_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node)
	return 0;
    const char *record_role = ged_obol_shape_node_record_role(node);
    if (record_role && (BU_STR_EQUAL(record_role, "view-feature") ||
			BU_STR_EQUAL(record_role, "view-polygon")))
	return 0;

    out->valid = 1;
    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	const char *intent_path = shape->ownerSourcePath.getValue().getString();
	out->is_database_source = shape->databaseIntent.getValue() ? 1 : 0;
	out->has_draw_intent = 1;
	out->intent_path = (intent_path && intent_path[0]) ? intent_path : path;
	out->intent_draw_mode =
	    ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
	out->visible = shape->visible.getValue() ? 1 : 0;
	out->selected = shape->selected.getValue() ? 1 : 0;
	out->highlighted = shape->highlighted.getValue() ? 1 : 0;
	out->line_style = shape->lineStyle.getValue();
	out->line_width = shape->lineWidth.getValue();
	out->transparency = shape->transparency.getValue();
	out->draw_mode =
	    ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
	const SbColor material = shape->materialColorValid.getValue() ?
				 shape->materialColor.getValue() : shape->color.getValue();
	ged_obol_rgb_from_color(material, out->material_color);
	out->material_valid = 1;
	return 1;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	const char *intent_path = shape->ownerSourcePath.getValue().getString();
	out->is_database_source = shape->databaseIntent.getValue() ? 1 : 0;
	out->has_draw_intent = 1;
	out->intent_path = (intent_path && intent_path[0]) ? intent_path : path;
	out->intent_draw_mode =
	    ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
	out->visible = shape->visible.getValue() ? 1 : 0;
	out->selected = shape->selected.getValue() ? 1 : 0;
	out->highlighted = shape->highlighted.getValue() ? 1 : 0;
	out->line_style = shape->lineStyle.getValue();
	out->line_width = shape->lineWidth.getValue();
	out->transparency = shape->transparency.getValue();
	out->draw_mode =
	    ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
	const SbColor material = shape->materialColorValid.getValue() ?
				 shape->materialColor.getValue() : shape->color.getValue();
	ged_obol_rgb_from_color(material, out->material_color);
	out->material_valid = 1;
	return 1;
    }

    memset(out, 0, sizeof(*out));
    return 0;
}

extern "C" int
ged_draw_obol_shape_geometry_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_shape_geometry_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node)
	return 0;

    out->valid = 1;
    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	const SoBRLVListShape *geom = shape->getGeometrySource();
	const char *kind = shape->geometryKind.getValue().getString();
	out->geometry_name = (kind && BU_STR_EQUAL(kind, "annotation")) ?
			     "annotation" : "line-set";
	out->point_count = static_cast<size_t>(geom->point.getNum());
	out->index_count = 0;
	return 1;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	const SoBRLMeshShape *geom = shape->getGeometrySource();
	out->geometry_name = "indexed-face-set";
	out->point_count = static_cast<size_t>(geom->point.getNum());
	out->index_count =
	    static_cast<size_t>(shape->getTriangleCount()) * 4;
	return 1;
    }

    memset(out, 0, sizeof(*out));
    return 0;
}

extern "C" int
ged_draw_obol_shape_surface_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_view_surface_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return 0;

    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
    const size_t triangle_count =
	static_cast<size_t>(shape->getTriangleCount());
    const SoBRLMeshShape *geom = shape->getGeometrySource();

    out->valid = 1;
    out->point_count = static_cast<size_t>(geom->point.getNum());
    out->normal_count = triangle_count * 3;
    out->index_count = triangle_count * 4;
    out->face_count = triangle_count;
    out->normals_per_index = 1;
    out->material_valid = 1;
    out->material_draw_mode =
	ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
    out->material_transparency = shape->transparency.getValue();
    out->material_highlighted = shape->highlighted.getValue() ? 1 : 0;
    const SbColor material = shape->materialColorValid.getValue() ?
			     shape->materialColor.getValue() : shape->color.getValue();
    ged_obol_rgb_from_color(material, out->material_color);
    out->cache_identity = ged_obol_hash_cstr(path);
    const char *source_path = shape->sourcePath.getValue().getString();
    if (!source_path || !source_path[0])
	source_path = shape->ownerSourcePath.getValue().getString();
    out->source_identity = ged_obol_hash_cstr(source_path);
    if (!out->source_identity)
	out->source_identity = out->cache_identity;
    return 1;
}

extern "C" int
ged_draw_obol_shape_surface_point_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    point_t out)
{
    if (!out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return 0;

    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    if (index >= static_cast<size_t>(geom->point.getNum()))
	return 0;

    const SbVec3f &point = geom->point[static_cast<int>(index)];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}

extern "C" int
ged_draw_obol_shape_surface_index_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    int *out)
{
    if (!out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return 0;

    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    const size_t triangle_count =
	static_cast<size_t>(shape->getTriangleCount());
    if (index >= triangle_count * 4)
	return 0;

    const size_t face_offset = index % 4;
    if (face_offset == 3) {
	*out = -1;
	return 1;
    }

    const size_t coord_index = (index / 4) * 3 + face_offset;
    if (coord_index >= static_cast<size_t>(geom->coordIndex.getNum()))
	return 0;
    *out = geom->coordIndex[static_cast<int>(coord_index)];
    return 1;
}

extern "C" int
ged_draw_obol_shape_line_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_view_line_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	return 0;

    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
    const SoBRLVListShape *geom = shape->getGeometrySource();
    out->valid = 1;
    out->point_count = static_cast<size_t>(geom->point.getNum());
    return 1;
}

extern "C" int
ged_draw_obol_shape_line_point_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    point_t out)
{
    if (!out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	return 0;

    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
    const SoBRLVListShape *geom = shape->getGeometrySource();
    if (index >= static_cast<size_t>(geom->point.getNum()))
	return 0;

    double precise[3];
    if (shape->getPrecisePoint(static_cast<int>(index), precise)) {
	out[0] = precise[0];
	out[1] = precise[1];
	out[2] = precise[2];
	return 1;
    }

    const SbVec3f &point = geom->point[static_cast<int>(index)];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}

extern "C" int
ged_draw_obol_shape_line_command_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    int *out)
{
    if (!out)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !path || !path[0])
	return 0;

    SoNode *node = scene->findShape(path);
    if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	return 0;

    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
    const SoBRLVListShape *geom = shape->getGeometrySource();
    if (index >= static_cast<size_t>(geom->command.getNum()))
	return 0;

    *out = ged_obol_vlist_command_to_ged(
	       geom->command[static_cast<int>(index)]);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_last_point_for_path(
    struct ged *gedp,
    const char *path,
    point_t out)
{
    if (!out)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    const SoBRLVListShape *geom = shape ? shape->getGeometrySource() : NULL;
    if (!geom || geom->point.getNum() <= 0)
	return 0;

    const SbVec3f &point = geom->point[geom->point.getNum() - 1];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}

extern "C" int
ged_draw_obol_database_source_line_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_view_line_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_source(source, path);

    out->valid = 1;
    if (!shape)
	return 1;

    const SoBRLVListShape *geom = shape->getGeometrySource();
    out->point_count = static_cast<size_t>(geom->point.getNum());
    out->cache_identity =
	ged_obol_hash_sb_string(shape->cacheIdentity.getValue());
    out->source_identity =
	ged_obol_hash_sb_string(shape->sourceIdentity.getValue());
    return 1;
}

extern "C" int
ged_draw_obol_database_source_line_point_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    point_t out)
{
    if (!out)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    const SoBRLVListShape *geom = shape ? shape->getGeometrySource() : NULL;
    if (!geom || index >= static_cast<size_t>(geom->point.getNum()))
	return 0;

    double precise[3];
    if (shape->getPrecisePoint(static_cast<int>(index), precise)) {
	out[0] = precise[0];
	out[1] = precise[1];
	out[2] = precise[2];
	return 1;
    }

    const SbVec3f &point = geom->point[static_cast<int>(index)];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}

extern "C" int
ged_draw_obol_database_source_line_command_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    int *out)
{
    if (!out)
	return 0;

    SoBRLVListShape *shape = ged_obol_owned_vlist_shape_for_path(gedp, path);
    const SoBRLVListShape *geom = shape ? shape->getGeometrySource() : NULL;
    if (!geom || index >= static_cast<size_t>(geom->command.getNum()))
	return 0;

    *out = ged_obol_vlist_command_to_ged(
	       geom->command[static_cast<int>(index)]);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_line_data_copy_for_path(
    struct ged *gedp,
    const char *path,
    point_t **points,
    int **commands,
    size_t *point_count)
{
    if (points)
	*points = NULL;
    if (commands)
	*commands = NULL;
    if (point_count)
	*point_count = 0;
    if (!points || !commands || !point_count)
	return 0;

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source || !source->hasCompactInstanceIndex())
	return 0;

    std::vector<SbVec3f> compactPoints;
    std::vector<int32_t> compactCommands;
    if (!source->copyCompactWireGeometry(compactPoints, compactCommands) ||
	compactPoints.size() != compactCommands.size())
	return 0;

    point_t *copiedPoints = (point_t *)bu_calloc(compactPoints.size(),
	sizeof(point_t), "GED Obol compact export points");
    int *copiedCommands = (int *)bu_calloc(compactCommands.size(),
	sizeof(int), "GED Obol compact export commands");
    for (size_t i = 0; i < compactPoints.size(); i++) {
	copiedPoints[i][X] = compactPoints[i][0];
	copiedPoints[i][Y] = compactPoints[i][1];
	copiedPoints[i][Z] = compactPoints[i][2];
	copiedCommands[i] = compactCommands[i] == 0 ?
	    GED_DRAW_VIEW_LINE_MOVE : GED_DRAW_VIEW_LINE_DRAW;
    }

    *points = copiedPoints;
    *commands = copiedCommands;
    *point_count = compactPoints.size();
    return 1;
}

extern "C" int
ged_draw_obol_database_source_surface_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_view_surface_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    SoBRLMeshShape *shape = ged_obol_owned_mesh_shape_for_source(source,
			    path, 0);
    if (!shape) {
	out->valid = 1;
	return 1;
    }

    int material_draw_mode =
	ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
    BObolDatabaseSourceSummary source_summary;
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (ged_obol_database_source_controller_summary_for_source(scene,
	    source, source_summary) && source_summary.valid)
	material_draw_mode =
	    ged_obol_database_source_exact_draw_mode_to_ged(gedp,
		source_summary, source);

    const size_t triangle_count =
	static_cast<size_t>(shape->getTriangleCount());
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    out->valid = 1;
    out->point_count = static_cast<size_t>(geom->point.getNum());
    out->normal_count = triangle_count * 3;
    out->index_count = triangle_count * 4;
    out->face_count = triangle_count;
    out->normals_per_index = 1;
    out->material_valid = 1;
    out->material_draw_mode = material_draw_mode;
    out->material_transparency = shape->transparency.getValue();
    out->material_highlighted = shape->highlighted.getValue() ? 1 : 0;
    const SbColor material = shape->materialColorValid.getValue() ?
			     shape->materialColor.getValue() : shape->color.getValue();
    ged_obol_rgb_from_color(material, out->material_color);
    out->cache_identity =
	ged_obol_hash_sb_string(shape->cacheIdentity.getValue());
    out->source_identity =
	ged_obol_hash_sb_string(shape->sourceIdentity.getValue());
    return 1;
}

extern "C" int
ged_draw_obol_database_source_surface_point_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    point_t out)
{
    if (!out)
	return 0;

    SoBRLMeshShape *shape = ged_obol_owned_mesh_shape_for_path(gedp, path, 0);
    const SoBRLMeshShape *geom = shape ? shape->getGeometrySource() : NULL;
    if (!geom || index >= static_cast<size_t>(geom->point.getNum()))
	return 0;

    const SbVec3f &point = geom->point[static_cast<int>(index)];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}

extern "C" int
ged_draw_obol_database_source_surface_index_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    int *out)
{
    if (!out)
	return 0;

    SoBRLMeshShape *shape = ged_obol_owned_mesh_shape_for_path(gedp, path, 0);
    if (!shape)
	return 0;

    const SoBRLMeshShape *geom = shape->getGeometrySource();
    const size_t triangle_count =
	static_cast<size_t>(shape->getTriangleCount());
    if (index >= triangle_count * 4)
	return 0;

    const size_t face_offset = index % 4;
    if (face_offset == 3) {
	*out = -1;
	return 1;
    }

    const size_t coord_index = (index / 4) * 3 + face_offset;
    if (coord_index >= static_cast<size_t>(geom->coordIndex.getNum()))
	return 0;
    *out = geom->coordIndex[static_cast<int>(coord_index)];
    return 1;
}

extern "C" int
ged_draw_obol_database_source_translate_vlist_for_path(
    struct ged *gedp,
    const char *path,
    const vect_t xlate)
{
    if (!xlate)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    if (!shape)
	return 0;

    return shape->translatePoints(SbVec3f(
				      static_cast<float>(xlate[0]),
				      static_cast<float>(xlate[1]),
				      static_cast<float>(xlate[2]))) ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_clear_vlist_for_path(
    struct ged *gedp,
    const char *path)
{
    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    const int cleared =
	scene->clearDatabaseSourceInstanceExternalPrimaryGeometry(
	    source_instance_key.c_str());
    return cleared >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_publish_primitive_wireframe_for_path(
    struct ged *gedp,
    const char *path,
    struct rt_db_internal *ip,
    const struct bg_tess_tol *ttol,
    const struct bn_tol *tol)
{
    if (!gedp || !path || !path[0] || !ip)
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    return scene->publishDatabaseSourceInstancePrimitiveWireframe(
	       source_instance_key.c_str(), ip, ttol, tol);
}

extern "C" int
ged_draw_obol_database_source_publish_line_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    const int *commands,
    size_t point_count)
{
    if (point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    std::vector<double> precise_points;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    precise_points.reserve(point_count * 3);
    for (size_t i = 0; i < point_count; i++) {
	const int command = commands ? commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	obol_commands.push_back(obol_command);
	precise_points.push_back(points[i][0]);
	precise_points.push_back(points[i][1]);
	precise_points.push_back(points[i][2]);
    }

    BObolExternalLineSet line_set;
    line_set.points = obol_points.empty() ? NULL : obol_points.data();
    line_set.commands = obol_commands.empty() ? NULL : obol_commands.data();
    line_set.precisePoints = precise_points.empty() ? NULL :
			     precise_points.data();
    line_set.count = static_cast<int>(point_count);
    const int published =
	scene->publishDatabaseSourceInstanceExternalLineSet(
	    source_instance_key.c_str(), line_set);
    return published > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_publish_annotation_line_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    const int *commands,
    size_t point_count)
{
    if (point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    std::vector<double> precise_points;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    precise_points.reserve(point_count * 3);
    for (size_t i = 0; i < point_count; i++) {
	const int command = commands ? commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	obol_commands.push_back(obol_command);
	precise_points.push_back(points[i][0]);
	precise_points.push_back(points[i][1]);
	precise_points.push_back(points[i][2]);
    }

    BObolExternalLineSet line_set;
    line_set.points = obol_points.empty() ? NULL : obol_points.data();
    line_set.commands = obol_commands.empty() ? NULL : obol_commands.data();
    line_set.precisePoints = precise_points.empty() ? NULL :
			     precise_points.data();
    line_set.count = static_cast<int>(point_count);
    line_set.sourceType = "annotation";
    line_set.geometryKind = "annotation";
    const int published =
	scene->publishDatabaseSourceInstanceExternalLineSet(
	    source_instance_key.c_str(), line_set);
    return published > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_publish_annotation_record_for_path(
    struct ged *gedp,
    const char *path,
    const point_t base_point,
    const point_t *annotation_points,
    size_t annotation_point_count,
    const struct ged_draw_obol_annotation_segment *segments,
    size_t segment_count,
    const point_t *line_points,
    const int *line_commands,
    size_t line_point_count)
{
    if ((annotation_point_count && !annotation_points) ||
	(segment_count && !segments) ||
	segment_count > static_cast<size_t>(INT_MAX) ||
	annotation_point_count > static_cast<size_t>(INT_MAX) ||
	line_point_count > static_cast<size_t>(INT_MAX) ||
	(line_point_count && !line_points))
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    std::vector<SbVec3f> obol_line_points;
    std::vector<int32_t> obol_line_commands;
    std::vector<double> precise_line_points;
    obol_line_points.reserve(line_point_count);
    obol_line_commands.reserve(line_point_count);
    precise_line_points.reserve(line_point_count * 3);
    for (size_t i = 0; i < line_point_count; i++) {
	const int command = line_commands ? line_commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_line_points.push_back(SbVec3f(
				       static_cast<float>(line_points[i][0]),
				       static_cast<float>(line_points[i][1]),
				       static_cast<float>(line_points[i][2])));
	obol_line_commands.push_back(obol_command);
	precise_line_points.push_back(line_points[i][0]);
	precise_line_points.push_back(line_points[i][1]);
	precise_line_points.push_back(line_points[i][2]);
    }

    std::vector<SbVec3f> obol_annotation_points;
    std::vector<double> precise_annotation_points;
    obol_annotation_points.reserve(annotation_point_count);
    precise_annotation_points.reserve(annotation_point_count * 3);
    if (annotation_point_count && annotation_points) {
	for (size_t i = 0; i < annotation_point_count; i++) {
	    obol_annotation_points.push_back(SbVec3f(
						 static_cast<float>(annotation_points[i][0]),
						 static_cast<float>(annotation_points[i][1]),
						 static_cast<float>(annotation_points[i][2])));
	    precise_annotation_points.push_back(annotation_points[i][0]);
	    precise_annotation_points.push_back(annotation_points[i][1]);
	    precise_annotation_points.push_back(annotation_points[i][2]);
	}
    }

    std::vector<BObolExternalAnnotationSegment> obol_segments;
    obol_segments.reserve(segment_count);
    for (size_t i = 0; i < segment_count; i++) {
	const struct ged_draw_obol_annotation_segment *seg = &segments[i];
	BObolExternalAnnotationSegment obol_seg;
	if (seg->kind == GED_DRAW_OBOL_ANNOTATION_SEGMENT_LINE)
	    obol_seg.kind = BObolExternalAnnotationSegment::SEGMENT_LINE;
	else if (seg->kind == GED_DRAW_OBOL_ANNOTATION_SEGMENT_TEXT)
	    obol_seg.kind = BObolExternalAnnotationSegment::SEGMENT_TEXT;
	else
	    obol_seg.kind = BObolExternalAnnotationSegment::SEGMENT_NONE;
	obol_seg.lineStart = seg->line_start;
	obol_seg.lineEnd = seg->line_end;
	obol_seg.textRefPoint = seg->text_ref_point;
	obol_seg.text = seg->text;
	obol_segments.push_back(obol_seg);
    }

    BObolExternalAnnotation annotation;
    annotation.basePoint = base_point ? SbVec3f(
			       static_cast<float>(base_point[0]),
			       static_cast<float>(base_point[1]),
			       static_cast<float>(base_point[2])) :
			   SbVec3f(0.0f, 0.0f, 0.0f);
    annotation.linePoints = obol_line_points.empty() ? NULL :
			    obol_line_points.data();
    annotation.lineCommands = obol_line_commands.empty() ? NULL :
			      obol_line_commands.data();
    annotation.preciseLinePoints = precise_line_points.empty() ? NULL :
				   precise_line_points.data();
    annotation.linePointCount = static_cast<int>(line_point_count);
    annotation.annotationPoints = obol_annotation_points.empty() ? NULL :
				  obol_annotation_points.data();
    annotation.preciseAnnotationPoints =
	precise_annotation_points.empty() ? NULL :
	precise_annotation_points.data();
    annotation.annotationPointCount =
	static_cast<int>(annotation_point_count);
    annotation.segments = obol_segments.empty() ? NULL : obol_segments.data();
    annotation.segmentCount = static_cast<int>(segment_count);
    const int published =
	scene->publishDatabaseSourceInstanceExternalAnnotation(
	    source_instance_key.c_str(), annotation);
    return published > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_annotation_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_obol_annotation_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    if (!shape)
	return 0;

    const char *kind = shape->geometryKind.getValue().getString();
    if (!kind || !BU_STR_EQUAL(kind, "annotation"))
	return 0;

    const SoBRLVListShape *geom = shape->getGeometrySource();
    out->valid = 1;
    out->point_count = static_cast<size_t>(geom->annotationPoint.getNum());
    out->segment_count =
	static_cast<size_t>(geom->annotationSegmentKind.getNum());
    out->cache_identity =
	ged_obol_hash_sb_string(shape->cacheIdentity.getValue());
    out->source_identity =
	ged_obol_hash_sb_string(shape->sourceIdentity.getValue());

    for (int i = 0; i < geom->annotationSegmentKind.getNum(); i++) {
	const int kind_value = geom->annotationSegmentKind[i];
	if (!out->line_segment_valid &&
	    kind_value == GED_DRAW_OBOL_ANNOTATION_SEGMENT_LINE) {
	    out->line_segment_valid = 1;
	    out->line_start = (i < geom->annotationSegmentStart.getNum()) ?
			      geom->annotationSegmentStart[i] : 0;
	    out->line_end = (i < geom->annotationSegmentEnd.getNum()) ?
			    geom->annotationSegmentEnd[i] : 0;
	}
	if (!out->text_segment_valid &&
	    kind_value == GED_DRAW_OBOL_ANNOTATION_SEGMENT_TEXT) {
	    out->text_segment_valid = 1;
	    out->text_ref_point = (i < geom->annotationTextRefPoint.getNum()) ?
				  geom->annotationTextRefPoint[i] : 0;
	    out->text = (i < geom->annotationText.getNum()) ?
			geom->annotationText[i].getString() : "";
	}
    }

    return 1;
}

extern "C" int
ged_draw_obol_database_source_annotation_point_at_for_path(
    struct ged *gedp,
    const char *path,
    size_t index,
    point_t out)
{
    if (!out)
	return 0;

    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    const SoBRLVListShape *geom = shape ? shape->getGeometrySource() : NULL;
    if (!geom || index >= static_cast<size_t>(geom->annotationPoint.getNum()))
	return 0;

    double precise[3];
    if (shape->getPreciseAnnotationPoint(static_cast<int>(index), precise)) {
	out[0] = precise[0];
	out[1] = precise[1];
	out[2] = precise[2];
	return 1;
    }

    const SbVec3f &point = geom->annotationPoint[static_cast<int>(index)];
    out[0] = point[0];
    out[1] = point[1];
    out[2] = point[2];
    return 1;
}


extern "C" int
ged_draw_obol_database_source_publish_auxiliary_line_set_for_path(
    struct ged *gedp,
    const char *path,
    const char *name,
    const point_t *points,
    const int *commands,
    size_t point_count,
    const struct ged_draw_scene_display_summary *display_state)
{
    if (!name || !name[0] || point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	const int command = commands ? commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	obol_commands.push_back(obol_command);
    }

    BObolAuxiliaryLineSetDisplayState obol_display;
    const BObolAuxiliaryLineSetDisplayState *obol_display_ptr = NULL;
    if (display_state && display_state->valid) {
	obol_display.valid = TRUE;
	obol_display.drawMode =
	    ged_obol_lod_draw_mode_from_ged(display_state->draw_mode);
	obol_display.visible = display_state->visible ? TRUE : FALSE;
	obol_display.highlighted = display_state->highlighted ? TRUE : FALSE;
	obol_display.lineStyle = display_state->line_style;
	obol_display.lineWidth = display_state->line_width;
	obol_display.transparency =
	    static_cast<float>(display_state->transparency);
	obol_display.materialColorValid =
	    display_state->material_valid ? TRUE : FALSE;
	if (display_state->material_valid)
	    obol_display.materialColor =
		ged_obol_color_from_rgb(display_state->material_color);
	obol_display.materialRevision = 0;
	obol_display_ptr = &obol_display;
    }

    int changed_any = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int changed =
	    scene->publishDatabaseSourceInstanceAuxiliaryLineSet(
		source_instance_key.c_str(),
		name,
		obol_points.empty() ? NULL : obol_points.data(),
		obol_commands.empty() ? NULL : obol_commands.data(),
		static_cast<int>(point_count),
		obol_display_ptr);
	if (changed < 0)
	    return 0;
	if (changed > 0)
	    changed_any = 1;
    }
    return changed_any;
}

extern "C" int
ged_draw_obol_database_source_publish_auxiliary_source_line_set_for_path(
    struct ged *gedp,
    const char *path,
    const char *source_path,
    const char *display_name,
    const point_t *points,
    const int *commands,
    size_t point_count,
    const struct ged_draw_scene_display_summary *display_state)
{
    if (!source_path || !source_path[0] ||
	point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> owner_instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (owner_instance_keys.empty())
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	const int command = commands ? commands[i] : -1;
	const int32_t obol_command =
	    ged_obol_vlist_command_from_ged(command, i);
	if (obol_command < 0)
	    return 0;

	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	obol_commands.push_back(obol_command);
    }

    BObolAuxiliaryLineSetDisplayState obol_display;
    const BObolAuxiliaryLineSetDisplayState *obol_display_ptr = NULL;
    if (display_state && display_state->valid) {
	obol_display.valid = TRUE;
	obol_display.drawMode =
	    ged_obol_lod_draw_mode_from_ged(display_state->draw_mode);
	obol_display.visible = display_state->visible ? TRUE : FALSE;
	obol_display.highlighted = display_state->highlighted ? TRUE : FALSE;
	obol_display.lineStyle = display_state->line_style;
	obol_display.lineWidth = display_state->line_width;
	obol_display.transparency =
	    static_cast<float>(display_state->transparency);
	obol_display.materialColorValid =
	    display_state->material_valid ? TRUE : FALSE;
	if (display_state->material_valid)
	    obol_display.materialColor =
		ged_obol_color_from_rgb(display_state->material_color);
	obol_display.materialRevision = 0;
	obol_display_ptr = &obol_display;
    }

    int changed_any = 0;
    for (const std::string &owner_instance_key : owner_instance_keys) {
	const int changed =
	    scene->publishDatabaseSourceInstanceAuxiliarySourceLineSet(
		owner_instance_key.c_str(),
		source_path,
		display_name ? display_name : source_path,
		obol_points.empty() ? NULL : obol_points.data(),
		obol_commands.empty() ? NULL : obol_commands.data(),
		static_cast<int>(point_count),
		obol_display_ptr);
	if (changed < 0)
	    return 0;
	if (changed > 0)
	    changed_any = 1;
    }
    if (!changed_any)
	return 0;

    SoBRLDatabaseSource *source = scene->findDatabaseSource(source_path);
    SoBRLVListShape *shape = source ? source->getRealizedShape(0) : NULL;
    if (shape)
	ged_obol_vlist_shape_set_precise_points(shape, points, point_count);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_clear_auxiliary_shapes_for_path(
    struct ged *gedp,
    const char *path)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int cleared = scene->clearDatabaseSourceInstanceAuxiliaryShapes(
				source_instance_key.c_str());
	if (cleared < 0)
	    return 0;
	if (cleared > 0)
	    applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_publish_point_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count)
{
    if (point_count > static_cast<size_t>(INT_MAX))
	return 0;
    if (point_count && !points)
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    std::vector<SbVec3f> obol_points;
    std::vector<double> precise_points;
    obol_points.reserve(point_count);
    precise_points.reserve(point_count * 3);
    for (size_t i = 0; i < point_count; i++) {
	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
	precise_points.push_back(points[i][0]);
	precise_points.push_back(points[i][1]);
	precise_points.push_back(points[i][2]);
    }

    BObolExternalPointSet point_set;
    point_set.points = obol_points.empty() ? NULL : obol_points.data();
    point_set.precisePoints = precise_points.empty() ? NULL :
			      precise_points.data();
    point_set.count = static_cast<int>(point_count);
    const int published =
	scene->publishDatabaseSourceInstanceExternalPointSet(
	    source_instance_key.c_str(), point_set);
    return published > 0 ? 1 : 0;
}

static int
ged_obol_indexed_face_finish(std::vector<int32_t> &face,
			     size_t point_count,
			     std::vector<int32_t> &triangles,
			     size_t *face_count,
			     unsigned int *face_stamp,
			     std::vector<unsigned int> &seen)
{
    if (face.size() < 3)
	return 0;

    for (size_t i = 1; i + 1 < face.size(); i++) {
	triangles.push_back(face[0]);
	triangles.push_back(face[i]);
	triangles.push_back(face[i + 1]);
    }

    face.clear();
    if (face_count)
	(*face_count)++;
    if (face_stamp && seen.size() == point_count) {
	if (*face_stamp == UINT_MAX) {
	    std::fill(seen.begin(), seen.end(), 0);
	    *face_stamp = 1;
	} else {
	    (*face_stamp)++;
	}
    }
    return 1;
}

static int
ged_obol_indexed_faces_to_triangles(const int *indices,
				    size_t index_count,
				    size_t point_count,
				    std::vector<int32_t> &triangles,
				    size_t *face_count_out,
				    size_t *vertex_index_count_out)
{
    if (!indices || !index_count || !point_count ||
	point_count > static_cast<size_t>(INT_MAX) ||
	index_count > static_cast<size_t>(INT_MAX))
	return 0;

    size_t face_count = 0;
    size_t vertex_index_count = 0;
    unsigned int face_stamp = 1;
    std::vector<unsigned int> seen(point_count, 0);
    std::vector<int32_t> face;

    for (size_t i = 0; i < index_count; i++) {
	const int idx = indices[i];
	if (idx < 0) {
	    if (idx != -1 || !ged_obol_indexed_face_finish(face,
		    point_count, triangles, &face_count, &face_stamp, seen))
		return 0;
	    continue;
	}

	if (static_cast<size_t>(idx) >= point_count)
	    return 0;
	if (seen[static_cast<size_t>(idx)] == face_stamp)
	    return 0;
	seen[static_cast<size_t>(idx)] = face_stamp;
	vertex_index_count++;
	face.push_back(static_cast<int32_t>(idx));
    }

    if (!face.empty() && !ged_obol_indexed_face_finish(face,
	    point_count, triangles, &face_count, &face_stamp, seen))
	return 0;
    if (!face_count || triangles.empty())
	return 0;

    if (face_count_out)
	*face_count_out = face_count;
    if (vertex_index_count_out)
	*vertex_index_count_out = vertex_index_count;
    return 1;
}

static int
ged_obol_indexed_face_triangle_normals(const vect_t *normals,
	size_t normal_count, const int *indices, size_t index_count,
	size_t point_count, size_t face_count, size_t vertex_index_count,
	const std::vector<int32_t> &triangles,
	std::vector<SbVec3f> &triangle_normals)
{
    triangle_normals.clear();
    if (!normal_count)
	return 1;
    if (!normals || !indices || !index_count)
	return 0;

    enum normal_binding {
	NORMAL_PER_CORNER,
	NORMAL_PER_RAW_INDEX,
	NORMAL_PER_POINT,
	NORMAL_PER_FACE
    };
    normal_binding binding;
    if (normal_count == vertex_index_count)
	binding = NORMAL_PER_CORNER;
    else if (normal_count == index_count)
	binding = NORMAL_PER_RAW_INDEX;
    else if (normal_count == point_count)
	binding = NORMAL_PER_POINT;
    else if (normal_count == face_count)
	binding = NORMAL_PER_FACE;
    else
	return 0;

    std::vector<int> face_indices;
    std::vector<SbVec3f> face_normals;
    size_t corner_index = 0;
    size_t face_index = 0;
    auto append_face = [&]() -> int {
	if (face_indices.size() < 3 || face_normals.size() != face_indices.size())
	    return 0;
	for (size_t i = 1; i + 1 < face_indices.size(); i++) {
	    triangle_normals.push_back(face_normals[0]);
	    triangle_normals.push_back(face_normals[i]);
	    triangle_normals.push_back(face_normals[i + 1]);
	}
	face_indices.clear();
	face_normals.clear();
	face_index++;
	return 1;
    };

    for (size_t i = 0; i < index_count; i++) {
	const int point_index = indices[i];
	if (point_index < 0) {
	    if (point_index != -1 || !append_face())
		return 0;
	    continue;
	}
	if (static_cast<size_t>(point_index) >= point_count)
	    return 0;
	size_t normal_index = 0;
	switch (binding) {
	    case NORMAL_PER_CORNER:
		normal_index = corner_index;
		break;
	    case NORMAL_PER_RAW_INDEX:
		normal_index = i;
		break;
	    case NORMAL_PER_POINT:
		normal_index = static_cast<size_t>(point_index);
		break;
	    case NORMAL_PER_FACE:
		normal_index = face_index;
		break;
	}
	if (normal_index >= normal_count)
	    return 0;
	face_indices.push_back(point_index);
	face_normals.push_back(SbVec3f(
	    static_cast<float>(normals[normal_index][X]),
	    static_cast<float>(normals[normal_index][Y]),
	    static_cast<float>(normals[normal_index][Z])));
	corner_index++;
    }
    if (!face_indices.empty() && !append_face())
	return 0;
    return face_index == face_count &&
	triangle_normals.size() == triangles.size();
}

extern "C" int
ged_draw_obol_local_shape_publish_indexed_face_set_for_path(
    struct ged *gedp,
    const char *group_path,
    const char *shape_path,
    const char *display_name,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count,
    const struct ged_draw_scene_display_summary *display_state)
{
    if (!group_path || !group_path[0] || !shape_path || !shape_path[0])
	return 0;
    if (!point_count || !index_count)
	return 0;
    if (!points || !indices || (normal_count && !normals) ||
	point_count > static_cast<size_t>(INT_MAX))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<int32_t> triangles;
    size_t face_count = 0;
    size_t vertex_index_count = 0;
    if (!ged_obol_indexed_faces_to_triangles(indices, index_count,
	    point_count, triangles, &face_count, &vertex_index_count))
	return 0;
    if (normal_count && normal_count != vertex_index_count &&
	normal_count != index_count &&
	normal_count != point_count && normal_count != face_count)
	return 0;
    if (triangles.size() > static_cast<size_t>(INT_MAX))
	return 0;

    SoBRLMeshShape *shape =
	ged_obol_local_mesh_shape_for_path(scene, group_path, shape_path);
    if (!shape)
	return 0;

    std::vector<SbVec3f> obol_points;
    obol_points.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
    }

    std::vector<SbVec3f> triangle_normals;
    if (!ged_obol_indexed_face_triangle_normals(normals, normal_count,
	indices, index_count, point_count, face_count, vertex_index_count,
	triangles, triangle_normals))
	return 0;

    ged_obol_local_shape_apply_common_state(shape, shape_path, display_name,
					    "local-indexed-face-set", "surface", display_state);
    shape->setIndexedTriangles(obol_points.data(),
			       static_cast<int>(obol_points.size()),
			       triangles.data(),
			       static_cast<int>(triangles.size()),
			       triangle_normals.empty() ? NULL : triangle_normals.data(),
			       static_cast<int>(triangle_normals.size()));
    return 1;
}

extern "C" int
ged_draw_obol_database_source_clear_mesh_for_path(
    struct ged *gedp,
    const char *path)
{
    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    const int cleared =
	scene->clearDatabaseSourceInstanceExternalPrimaryGeometry(
	    source_instance_key.c_str());
    return cleared >= 0 ? 1 : 0;
}

static int
ged_obol_database_source_publish_indexed_face_set_for_path(
    struct ged_bobol_publication_context *publication,
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count,
    int lod_backed)
{
    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key, publication))
	return 0;

    if (!point_count || !index_count) {
	BObolExternalTriangleMesh empty_mesh;
	const int published =
	    scene->publishDatabaseSourceInstanceExternalTriangleMesh(
		source_instance_key.c_str(), empty_mesh);
	return published > 0 ? 1 : 0;
    }

    if (!points || !indices || (normal_count && !normals) ||
	point_count > static_cast<size_t>(INT_MAX))
	return 0;

    std::vector<int32_t> triangles;
    size_t face_count = 0;
    size_t vertex_index_count = 0;
    if (!ged_obol_indexed_faces_to_triangles(indices, index_count,
	    point_count, triangles, &face_count, &vertex_index_count))
	return 0;
    if (normal_count && normal_count != vertex_index_count &&
	normal_count != index_count &&
	normal_count != point_count && normal_count != face_count)
	return 0;
    if (triangles.size() > static_cast<size_t>(INT_MAX))
	return 0;

    std::vector<SbVec3f> obol_points;
    obol_points.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	obol_points.push_back(SbVec3f(
				  static_cast<float>(points[i][0]),
				  static_cast<float>(points[i][1]),
				  static_cast<float>(points[i][2])));
    }

    BObolExternalTriangleMesh triangle_mesh;
    triangle_mesh.points = obol_points.empty() ? NULL : obol_points.data();
    triangle_mesh.pointCount = static_cast<int>(obol_points.size());
    triangle_mesh.indices = triangles.empty() ? NULL : triangles.data();
    triangle_mesh.indexCount = static_cast<int>(triangles.size());
    std::vector<SbVec3f> triangle_normals;
    if (!ged_obol_indexed_face_triangle_normals(normals, normal_count,
	indices, index_count, point_count, face_count, vertex_index_count,
	triangles, triangle_normals))
	return 0;
    triangle_mesh.normals = triangle_normals.empty() ? NULL :
	triangle_normals.data();
    triangle_mesh.normalCount = static_cast<int>(triangle_normals.size());
    triangle_mesh.lodBacked = lod_backed ? TRUE : FALSE;
    const int published =
	scene->publishDatabaseSourceInstanceExternalTriangleMesh(
	    source_instance_key.c_str(), triangle_mesh);
    return published > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_publish_indexed_face_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count)
{
	return ged_obol_database_source_publish_indexed_face_set_for_path(NULL,
	    gedp, path, points, point_count, normals, normal_count, indices,
	    index_count, 0);
}

int
ged_bobol_database_source_publish_indexed_face_set_for_path(
    struct ged_bobol_publication_context *publication,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count)
{
    if (!publication || !publication->active || !publication->gedp)
	return 0;
    return ged_obol_database_source_publish_indexed_face_set_for_path(
	publication, publication->gedp, path, points, point_count, normals,
	normal_count, indices, index_count, 0);
}

extern "C" int
ged_draw_obol_database_source_publish_lod_indexed_face_set_for_path(
    struct ged *gedp,
    const char *path,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count)
{
	return ged_obol_database_source_publish_indexed_face_set_for_path(NULL,
	    gedp, path, points, point_count, normals, normal_count, indices,
	    index_count, 1);
}

extern "C" int
ged_draw_obol_database_source_set_vlist_center_for_path(
    struct ged *gedp,
    const char *path,
    const point_t center)
{
    if (!center)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;

	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    return 0;

	const int changed = scene->setDatabaseSourceInstancePlacementState(
				source_instance_key.c_str(),
				summary.drawMatrixValid,
				summary.drawMatrix,
				TRUE,
				SbVec3f(
				    static_cast<float>(center[0]),
				    static_cast<float>(center[1]),
				    static_cast<float>(center[2])),
				summary.drawSizeValid,
				summary.drawSize);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_update_vlist_bounds_for_path(
    struct ged *gedp,
    const char *path)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;

	SoBRLVListShape *shape =
	    ged_obol_owned_vlist_shape_for_source(source, path, 0);
	if (!shape)
	    return 0;

	if (!shape->updateDrawBoundsFromPoints())
	    return 0;

	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid) {
	    applied = 1;
	    continue;
	}

	SbBool nextDrawMatrixValid = summary.drawMatrixValid ? TRUE : FALSE;
	SbMatrix nextDrawMatrix = summary.drawMatrix;
	if (!nextDrawMatrixValid && shape->drawCenterValid.getValue() &&
	    shape->drawSizeValid.getValue()) {
	    nextDrawMatrixValid = TRUE;
	    nextDrawMatrix = SbMatrix::identity();
	}

	(void)scene->setDatabaseSourceInstancePlacementState(
	    source_instance_key.c_str(),
	    nextDrawMatrixValid,
	    nextDrawMatrix,
	    shape->drawCenterValid.getValue(),
	    shape->drawCenter.getValue(),
	    shape->drawSizeValid.getValue(),
	    shape->drawSize.getValue());
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_set_placement_for_path(
    struct ged *gedp,
    const char *path,
    int draw_mat_valid,
    const mat_t draw_mat,
    int draw_center_valid,
    const point_t draw_center,
    int draw_size_valid,
    fastf_t draw_size)
{
    if ((draw_mat_valid && !draw_mat) ||
	(draw_center_valid && !draw_center))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;

	BObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    return 0;

	SbBool nextDrawMatrixValid = draw_mat_valid ?
				     TRUE : summary.drawMatrixValid;
	SbMatrix nextDrawMatrix = draw_mat_valid ?
				  ged_obol_sbmatrix_from_mat(draw_mat) :
				  summary.drawMatrix;
	SbBool nextDrawCenterValid = draw_center_valid ?
				     TRUE : summary.drawCenterValid;
	SbVec3f nextDrawCenter = draw_center_valid ?
				 SbVec3f(static_cast<float>(draw_center[0]),
					 static_cast<float>(draw_center[1]),
					 static_cast<float>(draw_center[2])) :
				 summary.drawCenter;
	SbBool nextDrawSizeValid = draw_size_valid ?
				   TRUE : summary.drawSizeValid;
	float nextDrawSize = draw_size_valid ?
			     static_cast<float>(draw_size) : summary.drawSize;

	const int changed = scene->setDatabaseSourceInstancePlacementState(
				source_instance_key.c_str(),
				nextDrawMatrixValid, nextDrawMatrix,
				nextDrawCenterValid,
				nextDrawCenter, nextDrawSizeValid,
				nextDrawSize);
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

static void
ged_obol_source_expansion_status_clear(
    struct ged_draw_obol_source_expansion_status *status)
{
    if (status)
	memset(status, 0, sizeof(*status));
}

static void
ged_obol_source_prewarm_status_clear(
    struct ged_draw_obol_source_prewarm_status *status)
{
    if (status)
	memset(status, 0, sizeof(*status));
}

static void
ged_obol_source_expansion_status_accumulate(
    struct ged_draw_obol_source_expansion_status *dst,
    const struct ged_draw_obol_source_expansion_status *src)
{
    if (!dst || !src)
	return;
    dst->child_count += src->child_count;
    dst->considered += src->considered;
    dst->expanded += src->expanded;
    dst->existing += src->existing;
    dst->skipped_non_union += src->skipped_non_union;
    dst->skipped_duplicate_instance += src->skipped_duplicate_instance;
    dst->expanded_non_union += src->expanded_non_union;
    dst->expanded_duplicate_instance +=
	src->expanded_duplicate_instance;
    dst->skipped_invalid += src->skipped_invalid;
    dst->remaining += src->remaining;
    dst->proxy_published += src->proxy_published;
    dst->metadata_applied += src->metadata_applied;
    dst->comb_sources += src->comb_sources;
    dst->leaf_sources += src->leaf_sources;
}

static void
ged_obol_source_prewarm_status_accumulate(
    struct ged_draw_obol_source_prewarm_status *dst,
    const struct ged_draw_obol_source_prewarm_status *src)
{
    if (!dst || !src)
	return;
    dst->child_count += src->child_count;
    dst->considered += src->considered;
    dst->submitted += src->submitted;
    dst->already_cached += src->already_cached;
    dst->skipped_non_union += src->skipped_non_union;
    dst->skipped_duplicate_instance += src->skipped_duplicate_instance;
    dst->shared_request += src->shared_request;
    dst->non_union_children += src->non_union_children;
    dst->duplicate_instances += src->duplicate_instances;
    dst->skipped_invalid += src->skipped_invalid;
    dst->remaining += src->remaining;
    dst->comb_sources += src->comb_sources;
    dst->leaf_sources += src->leaf_sources;
}

static const char *ged_obol_child_object_name(
    const struct ged_db_index_child *child);

struct ged_obol_source_frontier_entry {
    std::string path;
    std::string instance_key;
    size_t depth;
    size_t order;
    int priority;
};

static std::string ged_obol_child_instance_key(
    struct ged_view_context *view_ctx,
    const char *parent_instance_key,
    const struct ged_db_index_child *child,
    int draw_mode);

static std::vector<ged_obol_source_frontier_entry>
ged_obol_current_source_frontier(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *root_path,
    int draw_mode)
{
    std::vector<ged_obol_source_frontier_entry> frontier;
    if (!root_path || !root_path[0])
	return frontier;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return frontier;

    std::vector<std::string> instance_keys;
    ged_obol_collect_database_source_instance_keys_matching(instance_keys,
	scene, view_ctx, root_path, 0, 1, draw_mode);
    for (const std::string &instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    scene->findDatabaseSourceInstance(instance_key.c_str());
	BObolDatabaseSourceSummary summary;
	if (!source || !source->getSummary(summary) || !summary.valid ||
	    summary.path.getLength() == 0)
	    continue;

	size_t path_len = ged_db_index_path_resolve(gedp,
	    summary.path.getString(), NULL, 0);
	if (!path_len)
	    continue;
	std::vector<ged_db_index_id> path_ids(path_len);
	if (ged_db_index_path_resolve(gedp, summary.path.getString(),
		path_ids.data(), path_ids.size()) != path_len)
	    continue;
	const ged_db_index_id parent_id = path_ids.back();
	const size_t child_count = ged_db_index_child_count(gedp, parent_id);
	if (!child_count)
	    continue;

	int missing_priority = 0;
	for (size_t row = 0; row < child_count; row++) {
	    struct ged_db_index_child child;
	    memset(&child, 0, sizeof(child));
	    if (!ged_db_index_child_at(gedp, parent_id, row, &child))
		continue;
	    const std::string child_key = ged_obol_child_instance_key(view_ctx,
		instance_key.c_str(), &child, draw_mode);
	    SoBRLDatabaseSource *child_source = child_key.empty() ? NULL :
		scene->findDatabaseSourceInstance(child_key.c_str());
	    if (!child_key.empty() && !child_source)
		missing_priority = std::max(missing_priority,
		    child.record.is_comb ? 1 : 2);
	    if (child_source && !child.record.is_comb) {
		BObolDatabaseSourceSummary child_summary;
		if (!child_source->getSummary(child_summary) ||
		    !child_summary.valid ||
		    (child_summary.realizedShapeCount == 0 &&
		     child_summary.realizedMeshCount == 0)) {
		    missing_priority = std::max(missing_priority, 3);
		}
	    }
	}
	if (!missing_priority)
	    continue;

	ged_obol_source_frontier_entry entry;
	entry.path = summary.path.getString();
	entry.instance_key = instance_key;
	entry.depth = 0;
	entry.order = frontier.size();
	entry.priority = missing_priority;
	int in_component = 0;
	for (const char c : entry.path) {
	    if (c == '/') {
		in_component = 0;
	    } else if (!in_component) {
		entry.depth++;
		in_component = 1;
	    }
	}
	frontier.push_back(entry);
    }

    std::stable_sort(frontier.begin(), frontier.end(),
	[](const ged_obol_source_frontier_entry &a,
	   const ged_obol_source_frontier_entry &b) {
	    if (a.priority != b.priority)
		return a.priority > b.priority;
	    if (a.depth != b.depth)
		return a.depth > b.depth;
	    return a.order > b.order;
	});
    return frontier;
}

static std::string
ged_obol_child_instance_key(struct ged_view_context *view_ctx,
	const char *parent_instance_key,
	const struct ged_db_index_child *child,
	int draw_mode)
{
    if (!child || !child->record.valid)
	return std::string();

    std::string key = parent_instance_key && parent_instance_key[0] ?
	ged_obol_database_source_base_instance_key(parent_instance_key) :
	std::string();
    if (key.empty())
	return std::string();

    const char *occurrence_name = child->record.name;
    if (!occurrence_name || !occurrence_name[0])
	occurrence_name = ged_obol_child_object_name(child);
    if (!occurrence_name || !occurrence_name[0])
	return std::string();

    key += "/";
    key += occurrence_name;
    key += "#bool=";
    key += std::to_string(child->bool_op);
    if (draw_mode >= 0 && draw_mode != GED_DRAW_MODE_WIRE &&
	!ged_obol_view_scope_is_independent(view_ctx)) {
	key += ged_obol_database_source_mode_key_marker;
	key += std::to_string(draw_mode);
    }
    return key;
}

static void
ged_obol_aabb_proxy_line_set(const point_t bmin,
			     const point_t bmax,
			     point_t points[24],
			     int commands[24])
{
    point_t corners[8];
    static const int edges[12][2] = {
	{0, 1}, {1, 2}, {2, 3}, {3, 0},
	{4, 5}, {5, 6}, {6, 7}, {7, 4},
	{0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    VSET(corners[0], bmin[X], bmin[Y], bmin[Z]);
    VSET(corners[1], bmax[X], bmin[Y], bmin[Z]);
    VSET(corners[2], bmax[X], bmax[Y], bmin[Z]);
    VSET(corners[3], bmin[X], bmax[Y], bmin[Z]);
    VSET(corners[4], bmin[X], bmin[Y], bmax[Z]);
    VSET(corners[5], bmax[X], bmin[Y], bmax[Z]);
    VSET(corners[6], bmax[X], bmax[Y], bmax[Z]);
    VSET(corners[7], bmin[X], bmax[Y], bmax[Z]);

    for (size_t i = 0; i < 12; i++) {
	VMOVE(points[i * 2], corners[edges[i][0]]);
	VMOVE(points[i * 2 + 1], corners[edges[i][1]]);
	commands[i * 2] = GED_DRAW_VIEW_LINE_MOVE;
	commands[i * 2 + 1] = GED_DRAW_VIEW_LINE_DRAW;
    }
}

static int
ged_obol_apply_cached_path_metadata(struct ged *gedp,
				    const char *path,
				    struct ged_draw_obol_source_expansion_status *status)
{
    if (!gedp || !gedp->dbip || !path || !path[0])
	return 0;

    struct BObolDrawMetadataRecord metadata;
    bobol_draw_metadata_record_init(&metadata);
    const char *metadata_path = ged_obol_skip_leading_slash(path);
    if (!metadata_path || !metadata_path[0])
	return 0;

    if (bobol_draw_path_metadata_cache_get(gedp->dbip, metadata_path,
	    &metadata) != BRLCAD_OK)
	return 0;
    if (!metadata.directoryFound)
	return 0;

    if (!ged_draw_obol_database_source_apply_draw_metadata_for_path(gedp,
	    path, &metadata))
	return 0;

    if (status)
	status->metadata_applied++;
    return 1;
}

static int
ged_obol_apply_index_record_metadata(
    struct ged *gedp,
    const char *path,
    const struct ged_db_index_record *record,
    struct ged_draw_obol_source_expansion_status *status)
{
    if (!record || !record->valid || !record->dp)
	return ged_obol_apply_cached_path_metadata(gedp, path, status);

    struct BObolDrawMetadataRecord metadata;
    bobol_draw_metadata_record_init(&metadata);
    metadata.directoryFound = 1;
    metadata.isPhony = (record->dp->d_addr == RT_DIR_PHONY_ADDR) ? 1 : 0;
    metadata.flags = record->dp->d_flags;
    metadata.majorType = record->dp->d_major_type;
    metadata.minorType = record->dp->d_minor_type;
    metadata.isSolid = (record->dp->d_flags & RT_DIR_SOLID) ? 1 : 0;
    metadata.isComb = (record->dp->d_flags & RT_DIR_COMB) ? 1 : 0;
    metadata.isRegion = (record->dp->d_flags & RT_DIR_REGION) ? 1 : 0;
    metadata.isHidden = (record->dp->d_flags & RT_DIR_HIDDEN) ? 1 : 0;

    if (ged_draw_obol_database_source_apply_draw_metadata_for_path(gedp,
	    path, &metadata)) {
	if (status)
	    status->metadata_applied++;
	return 1;
    }

    return ged_obol_apply_cached_path_metadata(gedp, path, status);
}

static int
ged_obol_publish_aabb_bounds_for_instance(
    BObolSceneController *scene,
    const char *source_instance_key,
    const point_t bmin,
    const point_t bmax)
{
    if (!scene || !source_instance_key || !source_instance_key[0] ||
	!bmin || !bmax)
	return 0;

    point_t points[24];
    int commands[24];
    ged_obol_aabb_proxy_line_set(bmin, bmax, points, commands);

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    std::vector<double> precise_points;
    obol_points.reserve(24);
    obol_commands.reserve(24);
    precise_points.reserve(72);
    for (size_t i = 0; i < 24; i++) {
	obol_points.push_back(SbVec3f(
	    static_cast<float>(points[i][X]),
	    static_cast<float>(points[i][Y]),
	    static_cast<float>(points[i][Z])));
	obol_commands.push_back(ged_obol_vlist_command_from_ged(commands[i], i));
	precise_points.push_back(points[i][X]);
	precise_points.push_back(points[i][Y]);
	precise_points.push_back(points[i][Z]);
    }

    BObolExternalLineSet line_set;
    line_set.points = obol_points.data();
    line_set.commands = obol_commands.data();
    line_set.precisePoints = precise_points.data();
    line_set.count = 24;
    line_set.sourceType = "proxy";
    line_set.geometryKind = "aabb";
    return scene->publishDatabaseSourceInstanceExternalLineSet(
	       source_instance_key, line_set);
}


static int
ged_obol_publish_aabb_proxy_for_path(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *path,
    const char *source_instance_key,
    const char *cache_name,
    int draw_mode,
    int refresh_missing,
    struct ged_draw_obol_source_expansion_status *status)
{
    if (!gedp || !gedp->dbip || !path || !path[0] ||
	!cache_name || !cache_name[0])
	return 0;
    (void)view_ctx;
    (void)draw_mode;

    struct BObolDrawProxyRecord record;
    bobol_draw_proxy_record_init(&record);
    if (bobol_draw_proxy_cache_get(gedp->dbip, cache_name,
				     BOBOL_DRAW_CACHE_PROXY_AABB, &record) != BRLCAD_OK) {
	if (!refresh_missing)
	    return 0;
	if (bobol_draw_proxy_cache_refresh(gedp->dbip, cache_name,
					     BOBOL_DRAW_CACHE_PROXY_AABB, NULL) != BRLCAD_OK ||
	    bobol_draw_proxy_cache_get(gedp->dbip, cache_name,
					 BOBOL_DRAW_CACHE_PROXY_AABB, &record) != BRLCAD_OK)
	    return 0;
    }

    if (record.pointCount != 2)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    int published = ged_obol_publish_aabb_bounds_for_instance(scene,
	source_instance_key, record.points[0], record.points[1]);
    if (published && status)
	status->proxy_published++;
    return published;
}

struct ged_obol_cheap_local_bounds {
    point_t bmin;
    point_t bmax;
};

struct ged_obol_cheap_bounds_context {
    struct ged *gedp;
    struct bg_tess_tol ttol;
    struct bn_tol tol;
    /* Bounds are object-local and may be reused by every occurrence in the
     * initial structural proxy snapshot. */
    std::unordered_map<ged_db_index_id, ged_obol_cheap_local_bounds>
	objectBounds;
    std::unordered_set<ged_db_index_id> activeObjects;
    /* When cacheOnly is set, never compute a fresh AABB with rt_bound_instance
     * on the calling (UI) thread.  Treat an object with no cached proxy AABB as
     * "not yet available", record its name in missingNames, and return failure
     * so the caller can compute it on the LoD worker pool instead.  Keeps the
     * coarse outline's bounds work off the UI thread on a cold first draw. */
    int cacheOnly;
    std::unordered_set<std::string> missingNames;
};

static int
ged_obol_transform_bounds(point_t outMin, point_t outMax,
    const point_t localMin, const point_t localMax, const mat_t matrix)
{
    VSETALL(outMin, INFINITY);
    VSETALL(outMax, -INFINITY);
    for (int corner = 0; corner < 8; corner++) {
	point_t local = {
	    (corner & 1) ? localMax[X] : localMin[X],
	    (corner & 2) ? localMax[Y] : localMin[Y],
	    (corner & 4) ? localMax[Z] : localMin[Z]
	};
	point_t transformed;
	MAT4X3PNT(transformed, matrix, local);
	if (!isfinite(transformed[X]) || !isfinite(transformed[Y]) ||
	    !isfinite(transformed[Z]))
	    return 0;
	VMINMAX(outMin, outMax, transformed);
    }
    return 1;
}

static int
ged_obol_cheap_local_bounds_object(ged_obol_cheap_bounds_context &ctx,
    ged_db_index_id objectId, point_t boundsMin, point_t boundsMax)
{
    struct ged_db_index_record record;
    memset(&record, 0, sizeof(record));
    if (!ged_db_index_record_get(ctx.gedp, objectId, &record) ||
	!record.valid || !record.dp)
	return 0;

    const ged_db_index_id canonicalId = record.object_id ?
	record.object_id : record.id;

    auto cached = ctx.objectBounds.find(canonicalId);
    if (cached != ctx.objectBounds.end()) {
	VMOVE(boundsMin, cached->second.bmin);
	VMOVE(boundsMax, cached->second.bmax);
	return 1;
    }

    ged_obol_cheap_local_bounds local;
    const char *objectName = record.dp->d_namep;
    struct BObolDrawProxyRecord cachedProxy;
    bobol_draw_proxy_record_init(&cachedProxy);
    if (objectName && objectName[0] &&
	bobol_draw_proxy_cache_get(ctx.gedp->dbip, objectName,
	    BOBOL_DRAW_CACHE_PROXY_AABB, &cachedProxy) == BRLCAD_OK &&
	cachedProxy.pointCount == 2) {
	VMOVE(local.bmin, cachedProxy.points[0]);
	VMOVE(local.bmax, cachedProxy.points[1]);
	ctx.objectBounds.emplace(canonicalId, local);
	VMOVE(boundsMin, local.bmin);
	VMOVE(boundsMax, local.bmax);
	return 1;
    }

    /* Cache-only: do not compute a fresh AABB on this (UI) thread.  Record the
     * name so the caller can schedule the computation on the LoD worker pool,
     * and report the object as not-yet-available. */
    if (ctx.cacheOnly) {
	if (objectName && objectName[0])
	    ctx.missingNames.insert(objectName);
	return 0;
    }

    if (!record.is_comb) {
	mat_t identity;
	MAT_IDN(identity);
	if (rt_bound_instance(&local.bmin, &local.bmax, record.dp,
	    ctx.gedp->dbip, &ctx.ttol, &ctx.tol, &identity) < 0)
	    return 0;
	ctx.objectBounds.emplace(canonicalId, local);
	VMOVE(boundsMin, local.bmin);
	VMOVE(boundsMax, local.bmax);
	return 1;
    }

    if (!ctx.activeObjects.insert(canonicalId).second)
	return 0;
    VSETALL(boundsMin, INFINITY);
    VSETALL(boundsMax, -INFINITY);
    int haveBounds = 0;
    const size_t childCount = ged_db_index_child_count(ctx.gedp,
	canonicalId);
    for (size_t row = 0; row < childCount; row++) {
	struct ged_db_index_child child;
	memset(&child, 0, sizeof(child));
	if (!ged_db_index_child_at(ctx.gedp, canonicalId, row, &child) ||
	    !child.record.valid)
	    continue;
	mat_t childMatrix;
	point_t childMin;
	point_t childMax;
	const ged_db_index_id childId = child.record.object_id ?
	    child.record.object_id : child.record.id;
	if (!ged_obol_cheap_local_bounds_object(ctx, childId, childMin,
		childMax))
	    continue;
	if (child.matrix_valid)
	    MAT_COPY(childMatrix, child.matrix);
	else
	    MAT_IDN(childMatrix);
	point_t transformedMin;
	point_t transformedMax;
	if (!ged_obol_transform_bounds(transformedMin, transformedMax,
		childMin, childMax, childMatrix))
	    continue;
	VMINMAX(boundsMin, boundsMax, transformedMin);
	VMINMAX(boundsMin, boundsMax, transformedMax);
	haveBounds = 1;
    }
    ctx.activeObjects.erase(canonicalId);
    if (haveBounds) {
	VMOVE(local.bmin, boundsMin);
	VMOVE(local.bmax, boundsMax);
	ctx.objectBounds.emplace(canonicalId, local);
    }
    return haveBounds;
}

/* Initial proxy publication must remain bounded: it is the synchronous part
 * of an interactive draw.  These limits provide a useful, colored structural
 * outline without turning a large assembly walk into a full realization. */
static const size_t ged_obol_structural_proxy_max_nodes = 128;
static const size_t ged_obol_structural_proxy_max_proxies = 64;
static const size_t ged_obol_structural_proxy_max_depth = 2;
static const size_t ged_obol_structural_proxy_max_metadata = 64;

struct ged_obol_structural_proxy_node {
    std::string path;
    std::string objectName;
    std::string instanceKey;
    std::string parentInstanceKey;
    struct ged_db_index_record record;
    int boolOp;
    size_t row;
    /* This is relative to the root source.  The root source owns its own
     * placement matrix, so compact occurrences must not bake it in. */
    SbMatrix localMatrix;
    point_t boundsMin;
    point_t boundsMax;
    int publishBounds;
    int metadataValid;
    struct BObolDrawMetadataRecord metadata;
};

struct ged_obol_structural_proxy_context {
    struct ged *gedp;
    struct ged_view_context *viewCtx;
    int drawMode;
    uint32_t sourceRevision;
    ged_obol_cheap_bounds_context boundsContext;
    std::vector<ged_obol_structural_proxy_node> nodes;
    std::unordered_set<ged_db_index_id> persistedBounds;
    size_t proxyCount;
    size_t metadataCount;
    /* Traversal caps.  Defaulted to the bounded synchronous limits; the
     * progressive provider raises them over successive frames to deepen the
     * compact-occurrence coarse outline without exploding scene-graph nodes. */
    size_t maxNodes;
    size_t maxProxies;
    size_t maxDepth;
    size_t maxMetadata;
};

static void
ged_obol_structural_proxy_context_init_caps(
	struct ged_obol_structural_proxy_context &ctx)
{
    ctx.maxNodes = ged_obol_structural_proxy_max_nodes;
    ctx.maxProxies = ged_obol_structural_proxy_max_proxies;
    ctx.maxDepth = ged_obol_structural_proxy_max_depth;
    ctx.maxMetadata = ged_obol_structural_proxy_max_metadata;
}

static int
ged_obol_structural_proxy_metadata(struct ged_obol_structural_proxy_context &ctx,
	const char *path, struct BObolDrawMetadataRecord *metadata)
{
    if (!metadata || !ctx.gedp || !ctx.gedp->dbip || !path || !path[0] ||
	ctx.metadataCount >= ctx.maxMetadata)
	return 0;

    ctx.metadataCount++;
    bobol_draw_metadata_record_init(metadata);
    const char *cache_path = ged_obol_skip_leading_slash(path);
    if (!cache_path || !cache_path[0])
	return 0;
    if (bobol_draw_path_metadata_cache_get(ctx.gedp->dbip, cache_path,
	metadata) != BRLCAD_OK &&
	bobol_draw_path_metadata_cache_refresh(ctx.gedp->dbip, cache_path,
	NULL) != BRLCAD_OK)
	return 0;
    if (!metadata->directoryFound)
	return bobol_draw_path_metadata_cache_get(ctx.gedp->dbip, cache_path,
	metadata) == BRLCAD_OK && metadata->directoryFound;
    return 1;
}

static int
ged_obol_collect_structural_proxy_children(
	struct ged_obol_structural_proxy_context &ctx,
	ged_db_index_id parentId,
	const std::string &parentPath,
	const std::string &parentInstanceKey,
	const SbMatrix &parentLocalMatrix,
	size_t depth)
{
    if (!ctx.gedp || !ctx.gedp->dbip || parentPath.empty() ||
	parentInstanceKey.empty())
	return 0;

    const size_t childCount = ged_db_index_child_count(ctx.gedp, parentId);
    int collected = 0;
    for (size_t row = 0; row < childCount; row++) {
	if (ctx.nodes.size() >= ctx.maxNodes) {
	    break;
	}

	struct ged_db_index_child child;
	memset(&child, 0, sizeof(child));
	if (!ged_db_index_child_at(ctx.gedp, parentId, row, &child) ||
	    !child.record.valid)
	    continue;
	const char *childName = ged_obol_child_object_name(&child);
	if (!childName || !childName[0])
	    continue;

	std::string childPath(parentPath);
	if (!childPath.empty())
	    childPath += "/";
	childPath += childName;
	const std::string childInstanceKey = ged_obol_child_instance_key(
	    ctx.viewCtx, parentInstanceKey.c_str(), &child, ctx.drawMode);
	if (childInstanceKey.empty())
	    continue;

	const ged_db_index_id childId = child.record.object_id ?
	    child.record.object_id : child.record.id;

	ged_obol_structural_proxy_node node;
	node.path = childPath;
	node.objectName = childName;
	node.instanceKey = childInstanceKey;
	node.parentInstanceKey = parentInstanceKey;
	node.record = child.record;
	node.boolOp = child.bool_op;
	node.row = child.row;
	node.localMatrix = child.matrix_valid ?
	    ged_obol_sbmatrix_from_mat(child.matrix) : SbMatrix::identity();
	node.localMatrix.multRight(parentLocalMatrix);
	VSETALL(node.boundsMin, 0.0);
	VSETALL(node.boundsMax, 0.0);
	node.publishBounds = 0;
	node.metadataValid = ged_obol_structural_proxy_metadata(ctx,
	    childPath.c_str(), &node.metadata);

	/* Box each primitive leaf: the coarse outline is per-leaf so a leaf's
	 * box is replaced in place by its own realized geometry as it streams
	 * in (mergeCompactOccurrences matches by leaf path).  Combs are recursed
	 * through -- regions and material boundaries no longer stop the walk --
	 * except at the depth cap, where a comb is boxed as one node until a
	 * deeper deepening pass reaches its leaves. */
	const int semanticBoundary = !child.record.is_comb ||
	    depth >= ctx.maxDepth;
	if (semanticBoundary) {
	    if (ctx.proxyCount >= ctx.maxProxies) {
		break;
	    }
	    point_t boundsMin;
	    point_t boundsMax;
	    if (!ged_obol_cheap_local_bounds_object(ctx.boundsContext, childId,
		    boundsMin, boundsMax))
		continue;
	    VMOVE(node.boundsMin, boundsMin);
	    VMOVE(node.boundsMax, boundsMax);
	    /* The recursive walk is memoized in memory.  Only the bounded set of
	     * proxy roots published by this snapshot needs persistent warm-start
	     * storage; persisting every descendant made first draw dominated by
	     * thousands of individual LMDB writes. */
	    if (ctx.persistedBounds.insert(childId).second) {
		point_t cacheBounds[2];
		VMOVE(cacheBounds[0], boundsMin);
		VMOVE(cacheBounds[1], boundsMax);
		(void)bobol_draw_proxy_cache_store(ctx.gedp->dbip, childName,
		    BOBOL_DRAW_CACHE_PROXY_AABB, cacheBounds, 2, NULL);
	    }
	    node.publishBounds = 1;
	    ctx.proxyCount++;
	}

	ctx.nodes.push_back(node);
	collected = 1;
	if (!semanticBoundary &&
	    !ged_obol_collect_structural_proxy_children(ctx, childId,
		childPath, childInstanceKey, node.localMatrix, depth + 1))
	    continue;
    }
    return collected;
}

static std::shared_ptr<const Obol::PartGeometry>
ged_obol_aabb_proxy_geometry(const point_t bounds_min, const point_t bounds_max)
{
    if (!bounds_min || !bounds_max)
	return std::shared_ptr<const Obol::PartGeometry>();

    std::shared_ptr<Obol::PartGeometry> geometry(
	new Obol::PartGeometry);
    Obol::WireRep wire;
    const SbVec3f corners[8] = {
	SbVec3f(static_cast<float>(bounds_min[X]),
		static_cast<float>(bounds_min[Y]),
		static_cast<float>(bounds_min[Z])),
	SbVec3f(static_cast<float>(bounds_max[X]),
		static_cast<float>(bounds_min[Y]),
		static_cast<float>(bounds_min[Z])),
	SbVec3f(static_cast<float>(bounds_min[X]),
		static_cast<float>(bounds_max[Y]),
		static_cast<float>(bounds_min[Z])),
	SbVec3f(static_cast<float>(bounds_max[X]),
		static_cast<float>(bounds_max[Y]),
		static_cast<float>(bounds_min[Z])),
	SbVec3f(static_cast<float>(bounds_min[X]),
		static_cast<float>(bounds_min[Y]),
		static_cast<float>(bounds_max[Z])),
	SbVec3f(static_cast<float>(bounds_max[X]),
		static_cast<float>(bounds_min[Y]),
		static_cast<float>(bounds_max[Z])),
	SbVec3f(static_cast<float>(bounds_min[X]),
		static_cast<float>(bounds_max[Y]),
		static_cast<float>(bounds_max[Z])),
	SbVec3f(static_cast<float>(bounds_max[X]),
		static_cast<float>(bounds_max[Y]),
		static_cast<float>(bounds_max[Z]))
    };
    static const unsigned int edges[12][2] = {
	{0, 1}, {1, 3}, {3, 2}, {2, 0},
	{4, 5}, {5, 7}, {7, 6}, {6, 4},
	{0, 4}, {1, 5}, {2, 6}, {3, 7}
    };
    wire.segmentPoints.reserve(24);
    wire.segmentIds.reserve(12);
    for (unsigned int edge = 0; edge < 12; edge++) {
	wire.segmentPoints.push_back(corners[edges[edge][0]]);
	wire.segmentPoints.push_back(corners[edges[edge][1]]);
	wire.segmentIds.push_back(edge + 1);
    }
    wire.bounds = SbBox3f(corners[0], corners[7]);
    geometry->wire = std::move(wire);
    /* Structural bounds are conservative LoD proxies, not authored wire.
     * SoCADAssembly may render a depth-tested point when every AABB corner
     * projects into one pixel, while retaining the box for bounds and picks. */
    geometry->subpixelProxyEligible = true;
    return geometry;
}

/* A bounded structural proxy may represent an assembly rather than a region.
 * In that case db_full_path_color quite correctly reports the assembly's own
 * (usually white) state, while the visible geometry below it is colored by a
 * region-id table.  Pick a representative descendant region so the temporary
 * proxy reflects the database appearance it stands in for. */
static int
ged_obol_structural_proxy_descendant_material_color(
    const ged_obol_structural_proxy_context &ctx,
    ged_db_index_id parentId,
    const std::string &parentPath,
    SbColor &color,
    size_t &visited)
{
    if (!ctx.gedp || !ctx.gedp->dbip || parentPath.empty() ||
	visited >= ged_obol_structural_proxy_max_nodes)
	return 0;

    const size_t childCount = ged_db_index_child_count(ctx.gedp, parentId);
    for (size_t row = 0; row < childCount; row++) {
	if (visited++ >= ged_obol_structural_proxy_max_nodes)
	    return 0;
	struct ged_db_index_child child;
	memset(&child, 0, sizeof(child));
	if (!ged_db_index_child_at(ctx.gedp, parentId, row, &child) ||
	    !child.record.valid)
	    continue;
	const char *childName = ged_obol_child_object_name(&child);
	if (!childName || !childName[0])
	    continue;

	std::string childPath(parentPath);
	childPath += "/";
	childPath += childName;
	if (child.record.dp &&
	    (child.record.dp->d_flags & RT_DIR_REGION) &&
	    bobol_database_source_path_material_color(ctx.gedp->dbip,
		childPath.c_str(), color))
	    return 1;

	if (child.record.is_comb) {
	    const ged_db_index_id childId = child.record.object_id ?
		child.record.object_id : child.record.id;
	    if (ged_obol_structural_proxy_descendant_material_color(ctx,
		    childId, childPath, color, visited))
		return 1;
	}
    }
    return 0;
}

static int
ged_obol_structural_proxy_material_color(
    const ged_obol_structural_proxy_context &ctx,
    const ged_obol_structural_proxy_node &node,
    SbColor &color)
{
    if (!ctx.gedp || !ctx.gedp->dbip || node.path.empty())
	return 0;

    const size_t pathLength = ged_db_index_path_resolve(ctx.gedp,
	node.path.c_str(), NULL, 0);
    if (pathLength) {
	std::vector<ged_db_index_id> ids(pathLength);
	if (ged_db_index_path_resolve(ctx.gedp, node.path.c_str(), ids.data(),
		ids.size()) == pathLength) {
	    size_t visited = 0;
	    if (ged_obol_structural_proxy_descendant_material_color(ctx,
		    ids.back(), node.path, color, visited))
		return 1;
	}
    }

    return bobol_database_source_path_material_color(ctx.gedp->dbip,
	node.path.c_str(), color) ? 1 : 0;
}

static BObolCompactOccurrence
ged_obol_structural_proxy_occurrence(
    const ged_obol_structural_proxy_context &ctx,
    const ged_obol_structural_proxy_node &node)
{
    BObolCompactOccurrence occurrence;
    occurrence.geometry = ged_obol_aabb_proxy_geometry(node.boundsMin,
	node.boundsMax);
    occurrence.localTransform = node.localMatrix;
    occurrence.lodBacked = TRUE;
    occurrence.occurrenceIndex = static_cast<uint32_t>(node.row);
    occurrence.booleanOperation = node.boolOp == DB_OP_SUBTRACT ?
	SoBRLDatabaseSource::BOOLEAN_SUBTRACT :
	(node.boolOp == DB_OP_INTERSECT ?
	 SoBRLDatabaseSource::BOOLEAN_INTERSECT :
	 SoBRLDatabaseSource::BOOLEAN_UNION);

    BObolRealizedShapeSummary &summary = occurrence.summary;
    summary.valid = occurrence.geometry ? TRUE : FALSE;
    summary.shapeKind = BObolRealizedShapeSummary::SHAPE_VLIST;
    summary.path = node.path.c_str();
    summary.sourceName = node.objectName.c_str();
    summary.sourceType = "proxy";
    summary.sourceId = ctx.sourceRevision;
    summary.sourceIdentity = node.instanceKey.c_str();
    summary.databaseIntent = TRUE;
    summary.localSource = TRUE;
    summary.drawMode = ged_obol_lod_draw_mode_from_ged(ctx.drawMode);
    summary.recordRole = "lod-aabb";
    summary.geometryKind = "aabb";
    summary.visible = TRUE;
    summary.selectable = TRUE;
    summary.lodAvailable = TRUE;
    summary.lodActiveLevel = BOBOL_LOD_QUALITY_PROXY;
    summary.lodBoundsMin = SbVec3f(static_cast<float>(node.boundsMin[X]),
	static_cast<float>(node.boundsMin[Y]), static_cast<float>(node.boundsMin[Z]));
    summary.lodBoundsMax = SbVec3f(static_cast<float>(node.boundsMax[X]),
	static_cast<float>(node.boundsMax[Y]), static_cast<float>(node.boundsMax[Z]));
    summary.pointCount = 24;
    summary.segmentCount = 12;
    summary.commandCount = 24;
    summary.boundsValid = TRUE;
    summary.bounds = SbBox3f(summary.lodBoundsMin, summary.lodBoundsMax);
    if (node.metadataValid) {
	if (node.metadata.hasRegionId)
	    summary.regionId = node.metadata.regionId;
	if (node.metadata.hasAircode)
	    summary.airCode = node.metadata.aircode;
	if (node.metadata.hasMaterialId)
	    summary.materialId = node.metadata.materialId;
	if (node.metadata.hasLos)
	    summary.los = node.metadata.los;
	if (node.metadata.hasColor) {
	    summary.materialColorValid = TRUE;
	    summary.materialColor = ged_obol_color_from_rgb(node.metadata.color);
	}
	if (node.metadata.hasShader)
	    summary.materialShader = node.metadata.shader;
    }
    /* Cached path metadata describes authored object attributes, but it does
     * not by itself resolve an active region-id color table.  Structural
     * proxies are visible draw geometry and must use the same effective
     * full-path database color as their eventual detailed occurrences. */
    SbColor databaseColor;
    if (ged_obol_structural_proxy_material_color(ctx, node, databaseColor)) {
	summary.materialColorValid = TRUE;
	summary.materialColor = databaseColor;
    }
    return occurrence;
}

static BObolCompactOccurrence
ged_obol_structural_proxy_manifest_occurrence(
    const ged_obol_structural_proxy_context &ctx,
    const struct BObolDrawManifestOccurrence &cached)
{
    ged_obol_structural_proxy_node node;
    node.path = cached.path ? cached.path : "";
    node.objectName = cached.sourceName ? cached.sourceName : "";
    /* The compact index derives the stable occurrence key from the root
     * source, path, and registry order.  This path identity is retained for
     * diagnostics and future request construction. */
    node.instanceKey = node.path;
    node.boolOp = cached.booleanOperation;
    node.row = cached.occurrenceIndex;
    node.localMatrix = ged_obol_sbmatrix_from_mat(cached.localMatrix);
    VMOVE(node.boundsMin, cached.boundsMin);
    VMOVE(node.boundsMax, cached.boundsMax);
    node.publishBounds = 1;
    node.metadataValid = cached.metadataValid;
    if (node.metadataValid)
	node.metadata = cached.metadata;
    BObolCompactOccurrence occurrence =
	ged_obol_structural_proxy_occurrence(ctx, node);
    if (cached.sourceMeshRequestValid && cached.meshAssetPath &&
	cached.meshAssetPath[0] && cached.meshAssetName &&
	cached.meshAssetName[0]) {
	occurrence.sourceMeshRequestValid = TRUE;
	BObolSourceMeshRequest &request = occurrence.sourceMeshRequest;
	request.path = node.path.c_str();
	request.sourceName = node.objectName.c_str();
	request.sourceType = "bot";
	request.meshAssetPath = cached.meshAssetPath;
	request.meshAssetName = cached.meshAssetName;
	request.meshAssetBounds = SbBox3f(
	    SbVec3f(static_cast<float>(cached.meshAssetBoundsMin[X]),
		static_cast<float>(cached.meshAssetBoundsMin[Y]),
		static_cast<float>(cached.meshAssetBoundsMin[Z])),
	    SbVec3f(static_cast<float>(cached.meshAssetBoundsMax[X]),
		static_cast<float>(cached.meshAssetBoundsMax[Y]),
		static_cast<float>(cached.meshAssetBoundsMax[Z])));
	request.faceCount = cached.sourceFaceCount;
	request.pointCount = cached.sourcePointCount;
	request.bounds = occurrence.summary.bounds;
	request.databaseIntent = 1;
	request.localSource = 1;
	request.drawMode = occurrence.summary.drawMode;
	request.recordRole = "lod-source";
	request.geometryKind = "surface";
	request.lodAvailable = 1;
	const SbVec3f bmin = request.bounds.getMin();
	const SbVec3f bmax = request.bounds.getMax();
	request.lodBoundsMin = bmin;
	request.lodBoundsMax = bmax;
    }
    return occurrence;
}

static int
ged_obol_store_structural_proxy_manifest(
    struct ged *gedp,
    const std::string &rootPath,
    const ged_obol_structural_proxy_context &ctx)
{
    if (!gedp || !gedp->dbip || rootPath.empty())
	return 0;

    size_t count = 0;
    for (const ged_obol_structural_proxy_node &node : ctx.nodes) {
	if (node.publishBounds)
	    count++;
    }
    if (!count)
	return 0;

    struct BObolDrawManifest manifest;
    bobol_draw_manifest_init(&manifest);
    manifest.occurrenceCount = count;
    manifest.occurrences = static_cast<BObolDrawManifestOccurrence *>(
	bu_calloc(count, sizeof(*manifest.occurrences),
	    "Obol structural proxy manifest"));
    if (!manifest.occurrences)
	return 0;

    size_t index = 0;
    int valid = 1;
    for (const ged_obol_structural_proxy_node &node : ctx.nodes) {
	if (!node.publishBounds)
	    continue;
	struct BObolDrawManifestOccurrence &cached =
	    manifest.occurrences[index++];
	cached.path = bu_strdup(node.path.c_str());
	cached.sourceName = bu_strdup(node.objectName.c_str());
	if (!cached.path || !cached.sourceName) {
	    valid = 0;
	    break;
	}
	ged_obol_mat_from_sbmatrix(node.localMatrix, cached.localMatrix);
	VMOVE(cached.boundsMin, node.boundsMin);
	VMOVE(cached.boundsMax, node.boundsMax);
	cached.booleanOperation = node.boolOp;
	cached.occurrenceIndex = static_cast<uint32_t>(node.row);
	cached.metadataValid = node.metadataValid;
	if (cached.metadataValid)
	    cached.metadata = node.metadata;
    }

    const int stored = valid && index == count &&
	bobol_draw_manifest_cache_store(gedp->dbip, rootPath.c_str(),
	    &manifest) == BRLCAD_OK;
    bobol_draw_manifest_free(&manifest);
    return stored ? 1 : 0;
}

/* Persist the authoritative per-leaf occurrence bounds as one batched warm
 * start record.  Region/root boxes are derivable unions; leaf-local bounds and
 * occurrence transforms are the reusable facts needed by progressive drawing.
 * Keeping this as one manifest write avoids the thousands of tiny cache
 * transactions that motivated the old region-only shortcut. */
static int
ged_obol_store_leaf_proxy_manifest(
    struct ged *gedp,
    const SoBRLDatabaseSource *source)
{
    if (!gedp || !gedp->dbip || !source)
	return 0;
    const char *sourcePath = source->path.getValue().getString();
    const std::string rootPath = ged_obol_skip_leading_slash(sourcePath);
    if (rootPath.empty())
	return 0;

    std::vector<BObolCompactOccurrence> occurrences;
    const int occurrenceCount = source->getCompactInstanceCount();
    occurrences.reserve(occurrenceCount > 0 ?
	static_cast<size_t>(occurrenceCount) : 0);
    for (int i = 0; i < occurrenceCount; i++) {
	BObolCompactOccurrence occurrence;
	if (!source->getCompactOccurrence(i, occurrence) ||
	    !occurrence.summary.valid || !occurrence.summary.boundsValid ||
	    occurrence.summary.bounds.isEmpty() ||
	    occurrence.summary.path.getLength() == 0)
	    continue;
	occurrences.push_back(std::move(occurrence));
    }
    if (occurrences.empty())
	return 0;

    struct BObolDrawManifest manifest;
    bobol_draw_manifest_init(&manifest);
    manifest.occurrenceCount = occurrences.size();
    manifest.occurrences = static_cast<BObolDrawManifestOccurrence *>(
	bu_calloc(manifest.occurrenceCount, sizeof(*manifest.occurrences),
	    "Obol leaf proxy manifest"));
    if (!manifest.occurrences)
	return 0;

    int valid = 1;
    for (size_t i = 0; i < occurrences.size(); i++) {
	const BObolCompactOccurrence &occurrence = occurrences[i];
	struct BObolDrawManifestOccurrence &cached =
	    manifest.occurrences[i];
	cached.path = bu_strdup(occurrence.summary.path.getString());
	cached.sourceName =
	    bu_strdup(occurrence.summary.sourceName.getString());
	if (!cached.path || !cached.sourceName) {
	    valid = 0;
	    break;
	}
	ged_obol_mat_from_sbmatrix(occurrence.localTransform,
	    cached.localMatrix);
	const SbVec3f bmin = occurrence.summary.bounds.getMin();
	const SbVec3f bmax = occurrence.summary.bounds.getMax();
	VSET(cached.boundsMin, bmin[0], bmin[1], bmin[2]);
	VSET(cached.boundsMax, bmax[0], bmax[1], bmax[2]);
	cached.booleanOperation = occurrence.booleanOperation;
	cached.occurrenceIndex = occurrence.occurrenceIndex;
	if (occurrence.sourceMeshRequestValid) {
	    const BObolSourceMeshRequest &request =
		occurrence.sourceMeshRequest;
	    const char *assetPath = request.meshAssetPath.getLength() > 0 ?
		request.meshAssetPath.getString() :
		occurrence.summary.path.getString();
	    const char *assetName = request.meshAssetName.getLength() > 0 ?
		request.meshAssetName.getString() :
		occurrence.summary.sourceName.getString();
	    const SbBox3f assetBounds = !request.meshAssetBounds.isEmpty() ?
		request.meshAssetBounds :
		(!request.bounds.isEmpty() ? request.bounds :
		    occurrence.summary.bounds);
	    if (!assetPath || !assetPath[0] || !assetName ||
		!assetName[0] || assetBounds.isEmpty()) {
		valid = 0;
		break;
	    }
	    cached.sourceMeshRequestValid = 1;
	    cached.meshAssetPath = bu_strdup(assetPath);
	    cached.meshAssetName = bu_strdup(assetName);
	    if (!cached.meshAssetPath || !cached.meshAssetName) {
		valid = 0;
		break;
	    }
	    const SbVec3f assetMin = assetBounds.getMin();
	    const SbVec3f assetMax = assetBounds.getMax();
	    VSET(cached.meshAssetBoundsMin,
		assetMin[0], assetMin[1], assetMin[2]);
	    VSET(cached.meshAssetBoundsMax,
		assetMax[0], assetMax[1], assetMax[2]);
	    cached.sourceFaceCount = request.faceCount;
	    cached.sourcePointCount = request.pointCount;
	}
	/* Geometry realization already owns the full material semantics.  The
	 * leaf-box manifest is intentionally presentation-only; fabricating a
	 * partial directory record here could override correct inherited state. */
	cached.metadataValid = 0;
    }

    const int stored = valid &&
	bobol_draw_manifest_cache_store(gedp->dbip, rootPath.c_str(),
	    &manifest) == BRLCAD_OK;
    bobol_draw_manifest_free(&manifest);
    return stored ? 1 : 0;
}

static int
ged_obol_publish_deferred_structural_proxy_snapshot(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	const char *path,
	const char *root_instance_key,
	int draw_mode,
	size_t capProxies = 0,
	size_t capDepth = 0,
	size_t capNodes = 0,
	int cachedManifestOnly = 0)
{
    if (!gedp || !gedp->dbip || !path || !path[0] ||
	!root_instance_key || !root_instance_key[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    ged_obol_structural_proxy_context ctx;
    ctx.gedp = gedp;
    ctx.viewCtx = view_ctx;
    ctx.drawMode = draw_mode;
    ctx.sourceRevision = ged_obol_fold_revision(ged_draw_scene_revision(gedp));
    const struct bg_tess_tol defaultTtol = BG_TESS_TOL_INIT_TOL;
    ctx.boundsContext.gedp = gedp;
    ctx.boundsContext.ttol = defaultTtol;
    BN_TOL_INIT(&ctx.boundsContext.tol);
    ctx.boundsContext.cacheOnly = 1;
    ctx.proxyCount = 0;
    ctx.metadataCount = 0;
    /* Default to the bounded synchronous caps; the progressive provider passes
     * larger caps to deepen the flat region-level coarse outline over frames.
     * All occurrences live on the one root compact node, so raising these grows
     * the outline without adding scene-graph nodes. */
    ged_obol_structural_proxy_context_init_caps(ctx);
    if (capProxies)
	ctx.maxProxies = capProxies;
    if (capDepth)
	ctx.maxDepth = capDepth;
    if (capNodes)
	ctx.maxNodes = capNodes;
    else if (capProxies && capNodes == 0)
	ctx.maxNodes = capProxies > ctx.maxNodes ? capProxies * 4 : ctx.maxNodes;
    if (capProxies && ctx.maxMetadata < capProxies)
	ctx.maxMetadata = capProxies;

    const std::string rootPath = ged_obol_skip_leading_slash(path);
    SoBRLDatabaseSource *root = scene->findDatabaseSourceInstance(
	root_instance_key);
    if (!root)
	return 0;

    struct BObolDrawManifest manifest;
    bobol_draw_manifest_init(&manifest);
    if (bobol_draw_manifest_cache_get(gedp->dbip, rootPath.c_str(),
	&manifest) == BRLCAD_OK) {
	std::vector<BObolCompactOccurrence> occurrences;
	occurrences.reserve(manifest.occurrenceCount);
	for (size_t i = 0; i < manifest.occurrenceCount; i++) {
	    BObolCompactOccurrence occurrence =
		ged_obol_structural_proxy_manifest_occurrence(ctx,
		    manifest.occurrences[i]);
	    if (!occurrence.geometry) {
		occurrences.clear();
		break;
	    }
	    occurrences.push_back(std::move(occurrence));
	}
	bobol_draw_manifest_free(&manifest);
	if (!occurrences.empty()) {
	    ged_obol_scene_mutation_batch_scope batch(scene, 1,
		occurrences.size());
	    if (root->setCompactOccurrenceRegistry(occurrences) > 0)
		return ged_obol_database_source_mark_published_current(scene, root);
	}
    }

    /* The draw-command path may read one already-batched leaf manifest, but it
     * must never recurse through the hierarchy or calculate bounds on the UI
     * thread.  A cold miss proceeds directly to the detached leaf worker,
     * which streams leaf-local AABBs/geometry as each occurrence is visited. */
    if (cachedManifestOnly)
	return 0;

    if (!ged_db_index_available(gedp))
	return 0;
    const size_t pathLength = ged_db_index_path_resolve(gedp, path, NULL, 0);
    if (!pathLength)
	return 0;
    std::vector<ged_db_index_id> ids(pathLength);
    if (ged_db_index_path_resolve(gedp, path, ids.data(), ids.size()) !=
	pathLength)
	return 0;

    const int64_t collect_start = bu_gettime();
    const int collected_ok = ged_obol_collect_structural_proxy_children(ctx,
	ids.back(), rootPath, root_instance_key, SbMatrix::identity(), 1);
    ged_obol_timing_log("structural collect+bounds (cap)",
	collect_start, (long)ctx.proxyCount);
    /* Cache-only collection above skipped regions whose AABB was not cached;
     * compute those on the LoD worker pool.  A later provider pump re-runs this
     * snapshot and publishes them as the cache warms -- keeping rt_bound_instance
     * off the UI thread on a cold first draw. */
    const size_t scheduled_bounds = ged_obol_submit_region_bounds_async(
	gedp, view_ctx, ctx.boundsContext.missingNames, draw_mode);
    ged_obol_timing_log("structural: async bounds scheduled",
	collect_start, (long)ctx.boundsContext.missingNames.size());
    if (!collected_ok || !ctx.proxyCount) {
	/* Nothing publishable yet, but async work may now be in flight; tell the
	 * caller so it does not fall back to a synchronous whole-model bounds. */
	return scheduled_bounds > 0 ? -1 : 0;
    }

    std::vector<BObolCompactOccurrence> occurrences;
    occurrences.reserve(ctx.proxyCount);
    for (const ged_obol_structural_proxy_node &node : ctx.nodes) {
	if (!node.publishBounds)
	    continue;
	BObolCompactOccurrence occurrence =
	    ged_obol_structural_proxy_occurrence(ctx, node);
	if (occurrence.geometry)
	    occurrences.push_back(std::move(occurrence));
    }

    if (occurrences.empty())
	return scheduled_bounds > 0 ? -1 : 0;
    ged_obol_scene_mutation_batch_scope batch(scene, 1, occurrences.size());
    if (root->setCompactOccurrenceRegistry(occurrences) <= 0) {
	return scheduled_bounds > 0 ? -1 : 0;
	}
    /* Persist the manifest for instant warm-start only when the walk captured
     * the whole region set.  A truncated (capped) pass -- e.g. the small
     * synchronous first snapshot before the provider deepens it -- or a pass
     * with region bounds still being computed asynchronously must not freeze a
     * partial outline into the cache. */
    const int truncated =
	(ctx.proxyCount >= ctx.maxProxies) ||
	(ctx.nodes.size() >= ctx.maxNodes) ||
	!ctx.boundsContext.missingNames.empty();
    if (!truncated)
	(void)ged_obol_store_structural_proxy_manifest(gedp, rootPath, ctx);
    (void)ged_obol_database_source_mark_published_current(scene, root);
    /* Report async-pending so the caller keeps refining and does not fall back
     * to a synchronous whole-model bounds computation. */
    return scheduled_bounds > 0 ? -1 : 1;
}

static int
ged_obol_publish_deferred_root_proxy(struct ged *gedp,
				     struct ged_view_context *view_ctx,
				     const char *path,
				     const char *source_instance_key,
				     int draw_mode)
{
    if (!gedp || !gedp->dbip || !path || !path[0] ||
	!source_instance_key || !source_instance_key[0])
	return 0;

    const int structural_snapshot =
	ged_obol_publish_deferred_structural_proxy_snapshot(gedp, view_ctx,
	    path, source_instance_key, draw_mode, 0, 0, 0, 1);
    if (structural_snapshot != 0)
	return 1;

    /* A cached root AABB is a useful last-resort cold-start presentation.  A
     * cache miss is intentionally not refreshed here: recursive bounds work
     * belongs on a worker, and the detached realization will shortly stream
     * the first leaf presentation. */
    const char *cache_name = strrchr(path, '/');
    cache_name = cache_name ? cache_name + 1 : path;
    struct ged_draw_obol_source_expansion_status proxy_status;
    ged_obol_source_expansion_status_clear(&proxy_status);
    if (cache_name[0] &&
	ged_obol_publish_aabb_proxy_for_path(gedp, view_ctx, path,
	    source_instance_key, cache_name, draw_mode, 0, &proxy_status))
	return 1;
    return 0;
}

static const char *
ged_obol_child_object_name(const struct ged_db_index_child *child)
{
    if (!child || !child->record.valid)
	return NULL;
    if (child->record.dp && child->record.dp->d_namep &&
	child->record.dp->d_namep[0])
	return child->record.dp->d_namep;
    if (child->record.name && child->record.name[0] &&
	!strchr(child->record.name, '@'))
	return child->record.name;
    return NULL;
}

static size_t
ged_obol_submit_child_aabb_prewarm(
    BObolLodService *service,
    struct db_i *dbip,
    const char *database_id,
    uint64_t generation,
    const char *child_path,
    const char *child_name,
    int draw_mode,
    uint32_t source_revision)
{
    if (!service || !service->isRunning() || !dbip || !child_path ||
	!child_path[0] || !child_name || !child_name[0])
	return 0;

    BObolRtProxyProvider *provider = new (std::nothrow)
    BObolRtProxyProvider;
    if (!provider)
	return 0;

    provider->dbip = dbip;
    provider->proxyKind = BOBOL_LOD_PROXY_AABB;
    provider->useRequestBounds = FALSE;

    BObolLodTask task;
    task.generation = generation;
    task.request.databaseId = database_id ? database_id : "";
    task.request.sourceRevision = source_revision;
    task.request.objectPath = child_path;
    task.request.objectName = child_name;
    task.request.drawMode = ged_obol_lod_draw_mode_from_ged(draw_mode);
    task.request.providerId = "bobol_draw_aabb_cache";
    task.request.providerVersion = "bobol-cache-v1";
    task.request.qualityTier = BOBOL_LOD_QUALITY_PROXY;
    task.realize = bobol_rt_proxy_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = bobol_rt_proxy_provider_free;
    task.publishResult = FALSE;

    uint64_t task_id = service->submitIfNotActive(task);
    if (task_id == 0) {
	bobol_rt_proxy_provider_free(provider);
	return 0;
    }

    return 1;
}

static uint64_t
ged_obol_lod_service_incremental_generation(BObolLodService *service)
{
    if (!service)
	return 0;
    const uint64_t generation = service->currentGeneration();
    if (generation != 0 && !service->isGenerationCancelled(generation))
	return generation;
    return service->beginGeneration();
}

/* Schedule region-level AABB computation on the LoD worker pool for objects that
 * had no cached proxy bound during a cache-only structural pass.  Each task
 * computes rt_obj_bounds for one object and writes it to the draw proxy cache
 * (publishResult=FALSE), so a later cache-only pass on the UI thread can publish
 * that region's coarse box without ever calling rt_bound_instance itself. */
static size_t
ged_obol_submit_region_bounds_async(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const std::unordered_set<std::string> &names,
    int draw_mode)
{
    if (!gedp || !gedp->dbip || names.empty())
	return 0;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	return 0;
    BObolLodService *service = ged_obol_lod_service_ensure(controller);
    if (!service || !service->isRunning())
	return 0;
    const char *database_id = gedp->dbip->dbi_filename ?
	gedp->dbip->dbi_filename : "";
    const uint64_t generation =
	ged_obol_lod_service_incremental_generation(service);
    const uint32_t source_revision =
	ged_obol_fold_revision(ged_draw_scene_revision(gedp));
    size_t submitted = 0;
    static const size_t max_bounds_submissions_per_pass = 128;
    for (const std::string &name : names) {
	if (name.empty())
	    continue;
	submitted += ged_obol_submit_child_aabb_prewarm(service, gedp->dbip,
	    database_id, generation, name.c_str(), name.c_str(),
	    draw_mode, source_revision);
	if (submitted >= max_bounds_submissions_per_pass)
	    break;
    }
    return submitted;
}

extern "C" size_t
ged_draw_obol_database_source_prewarm_child_aabb_proxies(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *path,
    int draw_mode,
    size_t max_children,
    struct ged_draw_obol_source_prewarm_status *status)
{
    ged_obol_source_prewarm_status_clear(status);
    if (!gedp || !gedp->dbip || !path || !path[0])
	return 0;

    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	return 0;

    BObolLodService *service = controller->getLodService();
    if (!service || !service->isRunning()) {
	service = ged_obol_lod_service_ensure(controller);
	if (!service || !service->isRunning())
	    return 0;
    }

    size_t path_len = ged_db_index_path_resolve(gedp, path, NULL, 0);
    if (!path_len)
	return 0;

    std::vector<ged_db_index_id> path_ids(path_len);
    if (ged_db_index_path_resolve(gedp, path, path_ids.data(),
				  path_ids.size()) != path_len)
	return 0;

    const ged_db_index_id parent_id = path_ids.back();
    size_t child_count = ged_db_index_child_count(gedp, parent_id);
    if (status)
	status->child_count = child_count;
    if (!child_count)
	return 0;

    if (!max_children)
	max_children = child_count;

    const char *database_id = gedp->dbip->dbi_filename ?
			      gedp->dbip->dbi_filename : "";
    const uint64_t generation =
	ged_obol_lod_service_incremental_generation(service);
    const uint32_t source_revision =
	ged_obol_fold_revision(ged_draw_scene_revision(gedp));
    std::string parent_path = ged_obol_skip_leading_slash(path);
    size_t submitted = 0;
    std::unordered_set<std::string> requested_objects;

    for (size_t row = 0; row < child_count; row++) {
	struct ged_db_index_child child;
	memset(&child, 0, sizeof(child));
	if (!ged_db_index_child_at(gedp, parent_id, row, &child)) {
	    if (status)
		status->skipped_invalid++;
	    continue;
	}
	if (status)
	    status->considered++;

	if (status && child.bool_op != DB_OP_UNION)
	    status->non_union_children++;
	if (status && child.record.is_instance)
	    status->duplicate_instances++;

	const char *child_name = ged_obol_child_object_name(&child);
	if (!child_name || !child_name[0]) {
	    if (status)
		status->skipped_invalid++;
	    continue;
	}

	if (submitted >= max_children) {
	    if (status)
		status->remaining += child_count - row;
	    break;
	}

	if (child.record.is_comb) {
	    if (status)
		status->comb_sources++;
	    continue;
	}
	if (!requested_objects.insert(child_name).second) {
	    if (status)
		status->shared_request++;
	    continue;
	}

	struct BObolDrawCacheStatus cache_status;
	bobol_draw_cache_status_init(&cache_status);
	if (bobol_draw_proxy_cache_status(gedp->dbip, child_name,
					    BOBOL_DRAW_CACHE_PROXY_AABB, &cache_status) == BRLCAD_OK &&
	    cache_status.hasCachedPayload) {
	    if (status)
		status->already_cached++;
	    continue;
	}

	std::string child_path = parent_path;
	if (!child_path.empty())
	    child_path += "/";
	child_path += child_name;

	const size_t task_submitted = ged_obol_submit_child_aabb_prewarm(
					  service, gedp->dbip,
					  database_id, generation,
					  child_path.c_str(), child_name,
					  draw_mode, source_revision);
	if (!task_submitted) {
	    if (status)
		status->skipped_invalid++;
	    continue;
	}

	submitted += task_submitted;
	if (status) {
	    status->submitted += task_submitted;
	    if (child.record.is_comb)
		status->comb_sources += task_submitted;
	    else
		status->leaf_sources += task_submitted;
	}
    }

    return submitted;
}

extern "C" size_t
ged_draw_obol_database_source_prewarm_visible_child_aabb_proxies(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *root_path,
    int draw_mode,
    size_t max_sources,
    size_t max_children_per_source,
    struct ged_draw_obol_source_prewarm_status *status)
{
    ged_obol_source_prewarm_status_clear(status);
    if (!gedp || !root_path || !root_path[0])
	return 0;

    std::vector<ged_obol_source_frontier_entry> frontier =
	ged_obol_current_source_frontier(gedp, view_ctx, root_path, draw_mode);
    if (frontier.empty())
	return 0;

    size_t submitted = 0;
    size_t source_count = 0;
    for (const ged_obol_source_frontier_entry &entry : frontier) {
	if (max_sources && source_count >= max_sources)
	    break;
	source_count++;

	struct ged_draw_obol_source_prewarm_status path_status;
	ged_obol_source_prewarm_status_clear(&path_status);
	submitted += ged_draw_obol_database_source_prewarm_child_aabb_proxies(
			 gedp, view_ctx, entry.path.c_str(), draw_mode,
			 max_children_per_source, &path_status);
	ged_obol_source_prewarm_status_accumulate(status, &path_status);
    }

    return submitted;
}

static int
ged_obol_database_source_expand_children_impl(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *path,
    int draw_mode,
    size_t max_children,
    int refresh_missing_proxy,
    int require_cached_leaf_proxy,
    const char *parent_instance_key,
    struct ged_draw_obol_source_expansion_status *status)
{
    ged_obol_source_expansion_status_clear(status);
    if (!gedp || !gedp->dbip || !path || !path[0] ||
	!ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    size_t path_len = ged_db_index_path_resolve(gedp, path, NULL, 0);
    if (!path_len)
	return 0;

    std::vector<ged_db_index_id> path_ids(path_len);
    if (ged_db_index_path_resolve(gedp, path, path_ids.data(),
				  path_ids.size()) != path_len)
	return 0;

    const ged_db_index_id parent_id = path_ids.back();
    size_t child_count = ged_db_index_child_count(gedp, parent_id);
    if (status)
	status->child_count = child_count;
    if (!child_count)
	return 0;

    if (!max_children)
	max_children = child_count;

    const uint32_t source_revision =
	ged_obol_fold_revision(ged_draw_scene_revision(gedp));
    std::string parent_path = ged_obol_skip_leading_slash(path);
    const std::string effective_parent_instance_key =
	(parent_instance_key && parent_instance_key[0]) ?
	std::string(parent_instance_key) :
	ged_obol_database_source_instance_key_for_mode(view_ctx, path, draw_mode);
    int changed = 0;
    size_t expanded = 0;

    {
	ged_obol_scene_mutation_batch_scope batch(scene, max_children,
		max_children);
	for (size_t row = 0; row < child_count; row++) {
	    struct ged_db_index_child child;
	    memset(&child, 0, sizeof(child));
	    if (!ged_db_index_child_at(gedp, parent_id, row, &child)) {
		if (status)
		    status->skipped_invalid++;
		continue;
	    }
	    if (status)
		status->considered++;

	    const char *child_name = ged_obol_child_object_name(&child);
	    if (!child_name || !child_name[0]) {
		if (status)
		    status->skipped_invalid++;
		continue;
	    }

	    if (!child.record.is_comb && require_cached_leaf_proxy) {
		struct BObolDrawCacheStatus cache_status;
		bobol_draw_cache_status_init(&cache_status);
		if (bobol_draw_proxy_cache_status(gedp->dbip, child_name,
						    BOBOL_DRAW_CACHE_PROXY_AABB, &cache_status) != BRLCAD_OK ||
		    !cache_status.hasCachedPayload) {
		    if (status)
			status->remaining++;
		    continue;
		}
	    }

	    if (expanded >= max_children) {
		if (status)
		    status->remaining += child_count - row;
		break;
	    }

	    std::string child_path = parent_path;
	    if (!child_path.empty())
		child_path += "/";
	    child_path += child_name;

	    const std::string child_instance_key =
		ged_obol_child_instance_key(view_ctx,
		    effective_parent_instance_key.c_str(), &child, draw_mode);
	    if (child_instance_key.empty()) {
		if (status)
		    status->skipped_invalid++;
		continue;
	    }
	    SoBRLDatabaseSource *existing_child =
		scene->findDatabaseSourceInstance(child_instance_key.c_str());
	    if (existing_child) {
		if (status)
		    status->existing++;
		if (!child.record.is_comb) {
		    BObolDatabaseSourceSummary existing_summary;
		    const int has_geometry =
			existing_child->getSummary(existing_summary) &&
			existing_summary.valid &&
			(existing_summary.realizedShapeCount > 0 ||
			 existing_summary.realizedMeshCount > 0);
		    if (!has_geometry) {
			const int published =
			    ged_obol_publish_aabb_proxy_for_path(gedp, view_ctx,
				child_path.c_str(), child_instance_key.c_str(),
				child_name, draw_mode, refresh_missing_proxy,
				status);
			if (published)
			    changed = 1;
			else if (status)
			    status->remaining++;
		    }
		}
		continue;
	    }

	    const int replace_changed = ged_obol_replace_path(gedp, view_ctx,
					gedp->dbip, child_path.c_str(), draw_mode, source_revision,
					scene, 0, 0, NULL, NULL, NULL,
					child_instance_key.c_str());
	    if (replace_changed < 0) {
		if (status)
		    status->skipped_invalid++;
		continue;
	    }
	    if (replace_changed > 0)
		changed = 1;

	    SoBRLDatabaseSource *parent_source = scene->findDatabaseSourceInstance(
		effective_parent_instance_key.c_str());
	    BObolDatabaseSourceSummary parent_summary;
	    SbMatrix parent_matrix = SbMatrix::identity();
	    if (parent_source && parent_source->getSummary(parent_summary) &&
		parent_summary.valid && parent_summary.drawMatrixValid)
		parent_matrix = parent_summary.drawMatrix;
	    SbMatrix child_matrix = child.matrix_valid ?
		ged_obol_sbmatrix_from_mat(child.matrix) : SbMatrix::identity();
	    child_matrix.multRight(parent_matrix);
	    (void)scene->setDatabaseSourceInstancePlacementState(
		child_instance_key.c_str(), TRUE, child_matrix,
		FALSE, SbVec3f(0.0f, 0.0f, 0.0f), FALSE, 0.0f);
	    const int obol_bool_op = child.bool_op == DB_OP_SUBTRACT ?
		SoBRLDatabaseSource::BOOLEAN_SUBTRACT :
		(child.bool_op == DB_OP_INTERSECT ?
		 SoBRLDatabaseSource::BOOLEAN_INTERSECT :
		 SoBRLDatabaseSource::BOOLEAN_UNION);
	    (void)scene->setDatabaseSourceInstanceHierarchyState(
		child_instance_key.c_str(),
		effective_parent_instance_key.c_str(),
		static_cast<uint32_t>(child.row), obol_bool_op);
	    (void)ged_obol_apply_index_record_metadata(gedp, child_path.c_str(),
		    &child.record, status);

	    if (child.bool_op != DB_OP_UNION || child.record.is_instance) {
		SoBRLDatabaseSource *child_source =
		    scene->findDatabaseSourceInstance(child_instance_key.c_str());
		BObolDatabaseSourceSummary child_summary;
		if (child_source && child_source->getSummary(child_summary) &&
		    child_summary.valid) {
		    const int line_style = child.bool_op == DB_OP_SUBTRACT ? 1 :
			(child.bool_op == DB_OP_INTERSECT ? 2 :
			 child_summary.lineStyle);
		    (void)scene->setDatabaseSourceInstanceState(
			child_instance_key.c_str(), FALSE,
			child_summary.sourceRevision,
			child_summary.inputsRevision, child_summary.visible,
			child_summary.selected, child_summary.highlighted,
			line_style, child_summary.lineWidth,
			child_summary.transparency,
			child_summary.colorOverride, child_summary.color,
			child_summary.materialColorValid,
			child_summary.materialColor,
			child_summary.materialRevision);
		}
	    }

	    if (!child.record.is_comb) {
		const int published = ged_obol_publish_aabb_proxy_for_path(
		    gedp, view_ctx, child_path.c_str(),
		    child_instance_key.c_str(), child_name, draw_mode,
		    refresh_missing_proxy, status);
		if (published)
		    changed = 1;
		else if (status)
		    status->remaining++;
	    }

	    if (status) {
		status->expanded++;
		if (child.bool_op != DB_OP_UNION)
		    status->expanded_non_union++;
		if (child.record.is_instance)
		    status->expanded_duplicate_instance++;
		if (child.record.is_comb)
		    status->comb_sources++;
		else
		    status->leaf_sources++;
	    }
	    expanded++;
	}
    }

    if (changed) {
	BObolViewController *controller = ged_bobol_view_controller(view_ctx);
	if (controller) {
	    controller->clearViewLodState();
	    controller->requestRender("database-source-expand");
	}
    }

    return expanded > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_expand_children(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *path,
    int draw_mode,
    size_t max_children,
    struct ged_draw_obol_source_expansion_status *status)
{
    return ged_obol_database_source_expand_children_impl(gedp, view_ctx,
	    path, draw_mode, max_children, 1, 0, NULL, status);
}

static int
ged_obol_database_source_expand_visible_children_impl(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *root_path,
    int draw_mode,
    size_t max_sources,
    size_t max_children_per_source,
    int refresh_missing_proxy,
    int require_cached_leaf_proxy,
    struct ged_draw_obol_source_expansion_status *status)
{
    ged_obol_source_expansion_status_clear(status);
    if (!gedp || !root_path || !root_path[0])
	return 0;

    int changed = 0;
    size_t source_count = 0;
    while (!max_sources || source_count < max_sources) {
	std::vector<ged_obol_source_frontier_entry> frontier =
	    ged_obol_current_source_frontier(gedp, view_ctx, root_path,
		draw_mode);
	if (frontier.empty())
	    break;

	int pass_changed = 0;
	for (const ged_obol_source_frontier_entry &entry : frontier) {
	    if (max_sources && source_count >= max_sources)
		break;
	    source_count++;

	    struct ged_draw_obol_source_expansion_status path_status;
	    ged_obol_source_expansion_status_clear(&path_status);
	    const int path_changed =
		ged_obol_database_source_expand_children_impl(gedp, view_ctx,
		    entry.path.c_str(), draw_mode, max_children_per_source,
		    refresh_missing_proxy, require_cached_leaf_proxy,
		    entry.instance_key.c_str(), &path_status);
	    if (path_changed) {
		changed = 1;
		pass_changed = 1;
	    }
	    ged_obol_source_expansion_status_accumulate(status,
		&path_status);
	}

	if (!pass_changed)
	    break;
    }

    return changed;
}

extern "C" int
ged_draw_obol_database_source_expand_visible_children(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *root_path,
    int draw_mode,
    size_t max_sources,
    size_t max_children_per_source,
    struct ged_draw_obol_source_expansion_status *status)
{
    return ged_obol_database_source_expand_visible_children_impl(gedp,
	    view_ctx, root_path, draw_mode, max_sources,
	    max_children_per_source, 0, 1, status);
}

static int
ged_obol_progressive_autoview_apply(
    ged_obol_progressive_provider_data *data)
{
    if (!data || !data->gedp || !data->view_ctx ||
	!data->pending_autoview)
	return 0;

    struct bv *view = ged_obol_bv(data->view_ctx);
    if (!view) {
	data->pending_autoview = 0;
	return 0;
    }
    if (bv_frame_revision_get(view) != data->expected_view_revision) {
	data->pending_autoview = 0;
	return 0;
    }

    vect_t bmin;
    vect_t bmax;
    int empty = 1;
    if (!ged_draw_obol_scene_database_autoview_bounds(data->gedp, &bmin,
	&bmax, &empty, 0) || empty)
	return 0;

    if (!bv_autoview_bounds(view, data->autoview_factor, bmin, bmax))
	return 0;
    data->expected_view_revision = bv_frame_revision_get(view);
    bv_refresh_request(view, GED_VIEW_REFRESH_DRAW);
    return 1;
}

static int
ged_obol_progressive_autoview_arm(
	ged_obol_progressive_provider_data *data,
	fastf_t factor)
{
    if (!data || !data->view_ctx || data->deferred_refine_stage != 1 ||
	!data->deferred_job)
	return 0;

    const struct bv *view = ged_obol_bv_const(data->view_ctx);
    if (!view)
	return 0;

    data->pending_autoview = 1;
    data->expected_view_revision = bv_frame_revision_get(view);
    data->autoview_factor = factor;
    return 1;
}

extern "C" int
ged_draw_obol_progressive_autoview_follow(
	struct ged *gedp,
	struct ged_view_context *view_ctx,
	fastf_t factor)
{
    if (!gedp || !view_ctx)
	return 0;

    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	return 0;

    ged_obol_progressive_provider_data *data =
	static_cast<ged_obol_progressive_provider_data *>(
	    controller->findProgressiveProviderData(
		ged_obol_progressive_advance_provider));
    if (!ged_obol_progressive_autoview_arm(data, factor))
	return 0;
    controller->markProgressiveWorkPending();
    return 1;
}

static void
ged_obol_retire_deferred_job(ged_obol_progressive_provider_data *data)
{
    if (!data || !data->deferred_job)
	return;
    data->deferred_job->cancel();
    data->retired_jobs.push_back(data->deferred_job);
    data->deferred_job.reset();
}

static void
ged_obol_retire_all_deferred_jobs(ged_obol_progressive_provider_data *data)
{
    if (!data)
	return;

    ged_obol_retire_deferred_job(data);
    for (const std::shared_ptr<ged_obol_deferred_realization_job> &job :
	 data->pending_jobs) {
	if (!job)
	    continue;
	job->cancel();
	data->retired_jobs.push_back(job);
    }
    data->pending_jobs.clear();
}

static void
ged_obol_cleanup_retired_jobs(ged_obol_progressive_provider_data *data)
{
    if (!data)
	return;
    data->retired_jobs.erase(std::remove_if(data->retired_jobs.begin(),
	data->retired_jobs.end(),
	[](const std::shared_ptr<ged_obol_deferred_realization_job> &job) {
	    if (!job)
		return true;
	    const int state = job->state.load(std::memory_order_acquire);
	    return state == ged_obol_deferred_realization_job::COMPLETE ||
		state == ged_obol_deferred_realization_job::FAILED ||
		state == ged_obol_deferred_realization_job::CANCELLED;
	}), data->retired_jobs.end());
}

static SoBRLDatabaseSource *
ged_obol_find_deferred_source(BObolSceneController *scene,
    const std::string &instanceKey, const std::string &path,
    int representationMode)
{
    if (!scene)
	return NULL;
    SoBRLDatabaseSource *source = scene->findDatabaseSourceInstance(
	instanceKey.c_str());
    if (source && source->representationMode.getValue() == representationMode)
	return source;
    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid ||
	    summary.representationMode != representationMode ||
	    !ged_obol_path_equal(summary.path.getString(), path.c_str()))
	    continue;
	return scene->getDatabaseSource(i);
    }
    return NULL;
}

static int
ged_obol_start_deferred_realization(
    ged_obol_progressive_provider_data *data,
    BObolViewController *controller, int drawMode)
{
    if (!data || !controller || data->deferred_paths.empty())
	return 0;
    BObolSceneController *scene = controller->getSceneController();
    BObolSceneController *primaryScene = ged_draw_obol_scene_controller(
	data->gedp);
    if (!scene && !primaryScene) {
	bu_log("Obol deferred realization has no attached scene controller\n");
	return 0;
    }
    if (!scene)
	scene = primaryScene;

    if (data->deferred_job) {
	/* A second deferred mode is additive.  Its realization must not cancel
	 * an earlier mode, otherwise mixed shaded/wire draws get stuck at the
	 * earlier root proxy. */
	data->pending_jobs.push_back(data->deferred_job);
	data->deferred_job.reset();
    }
    std::shared_ptr<ged_obol_deferred_realization_job> job =
	std::make_shared<ged_obol_deferred_realization_job>();
    const int representationMode =
	ged_obol_database_representation_mode_from_ged(drawMode);
    for (const std::string &path : data->deferred_paths) {
	const std::string key = ged_obol_database_source_instance_key_for_mode(
	    data->view_ctx, path.c_str(), drawMode);
	SoBRLDatabaseSource *live = NULL;
	SbBool usePrimaryScene = FALSE;
	if (primaryScene) {
	    live = ged_obol_find_deferred_source(primaryScene, key, path,
		representationMode);
	    if (live)
		usePrimaryScene = TRUE;
	}
	if (!live) {
	    live = ged_obol_find_deferred_source(scene, key, path,
		representationMode);
	}
	if (!live) {
	    bu_log("Obol deferred realization source '%s' was not found\n",
		key.c_str());
	    return 0;
	}
	std::unique_ptr<ged_obol_deferred_realization_item> item(
	    new ged_obol_deferred_realization_item);
	const int64_t snapshot_start = bu_gettime();
	item->source = live->createDetachedRealizationSource(&item->database,
	    &item->snapshotPath);
	ged_obol_timing_log("deferred: DB snapshot (detach)", snapshot_start, -1);
	if (!item->source || !item->database) {
	    bu_log("Obol deferred realization could not snapshot source '%s'\n",
		key.c_str());
	    return 0;
	}
	item->instanceKey = live->instanceKey.getValue().getString();
	item->primaryScene = usePrimaryScene;
	item->allowWireFallback =
	    data->deferred_appearance.strict_fallback ? FALSE : TRUE;
	item->stream = std::make_shared<BObolCompactOccurrenceStream>();
	job->items.push_back(std::move(item));
    }
    if (!job->start()) {
	bu_log("Obol deferred realization worker could not be started\n");
	return 0;
    }
    data->deferred_job = job;
    return 1;
}

static int
ged_obol_remove_database_source_descendants(
    BObolSceneController *scene,
    const std::string &root_instance_key)
{
    if (!scene || root_instance_key.empty())
	return 0;

    /* Proxy-frontier sources are independent scene nodes, but their hierarchy
     * key records which deferred root owns them.  Remove deepest nodes first
     * so a completed compact root cannot leave proxy geometry behind. */
    std::unordered_set<std::string> known_keys;
    known_keys.insert(root_instance_key);
    std::vector<std::string> descendants;
    if (ged_obol_timing_enabled()) {
	bu_log("[obol-timing] retire proxy descendants: root=%s sources=%d\n",
	    root_instance_key.c_str(), scene->getDatabaseSourceCount());
	for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	    BObolDatabaseSourceSummary summary;
	    if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
		continue;
	    bu_log("[obol-timing] scene source: key=%s parent=%s path=%s "
		"meshes=%d shapes=%d\n",
		summary.instanceKey.getString(),
		summary.parentInstanceKey.getString(),
		summary.path.getString(), summary.realizedMeshCount,
		summary.realizedShapeCount);
	}
    }
    int found = 1;
    while (found) {
	found = 0;
	const int source_count = scene->getDatabaseSourceCount();
	for (int i = 0; i < source_count; i++) {
	    BObolDatabaseSourceSummary summary;
	    if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
		continue;
	    const char *instance_key = summary.instanceKey.getString();
	    const char *parent_key = summary.parentInstanceKey.getString();
	    if (!instance_key || !instance_key[0] || !parent_key || !parent_key[0] ||
		known_keys.find(parent_key) == known_keys.end())
		continue;
	    if (known_keys.insert(instance_key).second) {
		descendants.push_back(instance_key);
		found = 1;
	    }
	}
    }

    if (descendants.empty())
	return 0;
    std::reverse(descendants.begin(), descendants.end());
    return ged_obol_remove_instance_keys(descendants, scene);
}

static int
ged_obol_publish_deferred_realization(
    ged_obol_progressive_provider_data *data,
    BObolViewController *controller,
    const std::shared_ptr<ged_obol_deferred_realization_job> &job)
{
    if (!data || !data->gedp || !controller || !job)
	return 0;
    job->join();
    if (job->state.load(std::memory_order_acquire) !=
	ged_obol_deferred_realization_job::COMPLETE)
	return 0;

    BObolSceneController *scene = controller->getSceneController();
    BObolSceneController *primaryScene = ged_draw_obol_scene_controller(
	data->gedp);
    if (!scene && !primaryScene)
	return 0;
    if (!scene)
	scene = primaryScene;
    std::vector<SoBRLDatabaseSource *> liveSources;
    liveSources.reserve(job->items.size());
    for (const std::unique_ptr<ged_obol_deferred_realization_item> &item :
	 job->items) {
	BObolSceneController *itemScene = item && item->primaryScene ?
	    primaryScene : scene;
	SoBRLDatabaseSource *live = item && itemScene ?
	    itemScene->findDatabaseSourceInstance(item->instanceKey.c_str()) :
	    NULL;
	SoBRLDatabaseSource *detached = item ? item->source : NULL;
	if (!live || !detached)
	    return 0;
	if (live->sourceRevision.getValue() !=
		detached->sourceRevision.getValue() ||
	    live->inputsRevision.getValue() !=
		detached->inputsRevision.getValue() ||
	    live->viewRevision.getValue() != detached->viewRevision.getValue() ||
	    live->drawMode.getValue() != detached->drawMode.getValue() ||
	    live->representationMode.getValue() !=
		detached->representationMode.getValue()) {
	    if (ged_obol_timing_enabled()) {
		bu_log("[obol-timing] deferred adoption rejected for %s: "
		    "source=%u/%u inputs=%u/%u view=%u/%u mode=%d/%d "
		    "representation=%d/%d\n",
		    item->instanceKey.c_str(),
		    live->sourceRevision.getValue(),
		    detached->sourceRevision.getValue(),
		    live->inputsRevision.getValue(),
		    detached->inputsRevision.getValue(),
		    live->viewRevision.getValue(),
		    detached->viewRevision.getValue(),
		    live->drawMode.getValue(), detached->drawMode.getValue(),
		    live->representationMode.getValue(),
		    detached->representationMode.getValue());
	    }
	    return 0;
	}
	liveSources.push_back(live);
    }

    int adopted = 0;
    for (size_t i = 0; i < liveSources.size(); i++) {
	const int adopted_count = liveSources[i]->adoptDetachedCompactRealization(
	    job->items[i]->source);
	if (ged_obol_timing_enabled())
	    bu_log("[obol-timing] deferred adoption for %s: n=%d\n",
		job->items[i]->instanceKey.c_str(), adopted_count);
	if (adopted_count > 0) {
	    adopted++;
	    const int64_t manifest_start = bu_gettime();
	    const int manifest_stored =
		ged_obol_store_leaf_proxy_manifest(data->gedp,
		    liveSources[i]);
	    ged_obol_timing_log("deferred: store leaf manifest",
		manifest_start, manifest_stored ? adopted_count : 0);
	    /* A per-view worker can realize a synchronized scene while the
	     * structural frontier is owned by the primary scene.  Retire matching
	     * proxies from both; either side may be the adopted live source. */
	    if (primaryScene)
		(void)ged_obol_remove_database_source_descendants(primaryScene,
		    job->items[i]->instanceKey);
	    if (scene && scene != primaryScene)
		(void)ged_obol_remove_database_source_descendants(scene,
		    job->items[i]->instanceKey);
	}
    }
    return adopted;
}

/* A streamed merge is an atomic scene mutation, so the controller's time
 * budget can only be checked between batches.  Keep each batch bounded even
 * when a caller explicitly removes the total item cap: workers can fill a
 * queue while the GUI is busy, and draining that entire backlog in one paint
 * turns "progressive" work into a multi-second input stall. */
static const size_t GED_OBOL_STREAM_BATCH_QUANTUM = 64;

/* Drain completed per-leaf occurrences from every running realization job and
 * merge them onto the live compact root as they arrive, so geometry streams in
 * incrementally instead of appearing all at once when the whole walk finishes.
 * Returns 1 if any geometry was merged; sets *has_more when a stream still has
 * queued occurrences beyond this tick's cap. */
static int
ged_obol_drain_streamed_realizations(
    ged_obol_progressive_provider_data *data,
    BObolViewController *controller,
    size_t max_items,
    uint64_t max_microseconds,
    int *has_more)
{
    if (!data || !controller)
	return 0;
    BObolSceneController *scene = controller->getSceneController();
    BObolSceneController *primaryScene = ged_draw_obol_scene_controller(
	data->gedp);
    if (!scene && !primaryScene)
	return 0;
    if (!scene)
	scene = primaryScene;

    std::vector<std::shared_ptr<ged_obol_deferred_realization_job>> jobs;
    if (data->deferred_job)
	jobs.push_back(data->deferred_job);
    for (const std::shared_ptr<ged_obol_deferred_realization_job> &job :
	 data->pending_jobs)
	jobs.push_back(job);

    const int64_t started = bu_gettime();
    size_t drained_count = 0;
    int merged = 0;
    for (const std::shared_ptr<ged_obol_deferred_realization_job> &job : jobs) {
	if (!job)
	    continue;
	for (const std::unique_ptr<ged_obol_deferred_realization_item> &item :
	     job->items) {
	    if (!item || !item->stream)
		continue;
	    BObolSceneController *itemScene = item->primaryScene ?
		primaryScene : scene;
	    SoBRLDatabaseSource *live = itemScene ?
		itemScene->findDatabaseSourceInstance(item->instanceKey.c_str()) :
		NULL;
	    SoBRLDatabaseSource *detached = item->source;
	    if (!live || !detached)
		continue;
	    /* A source/draw-mode change after this job launched makes its stream
	     * stale; the job will be retired by the pump's revision drain.  Do not
	     * merge stale geometry onto the current live source. */
	    if (live->sourceRevision.getValue() !=
		    detached->sourceRevision.getValue() ||
		live->inputsRevision.getValue() !=
		    detached->inputsRevision.getValue() ||
		live->drawMode.getValue() != detached->drawMode.getValue() ||
		live->representationMode.getValue() !=
		    detached->representationMode.getValue())
		continue;

	    size_t batch_cap = GED_OBOL_STREAM_BATCH_QUANTUM;
	    if (max_items) {
		if (drained_count >= max_items) {
		    if (has_more)
			*has_more = 1;
		    return merged;
		}
		batch_cap = std::min(batch_cap, max_items - drained_count);
	    }
	    std::vector<BObolCompactOccurrence> batch;
	    (void)item->stream->drain(batch, batch_cap);
	    if (batch.empty())
		continue;
	    drained_count += batch.size();
	    ged_obol_scene_mutation_batch_scope mutation(itemScene, 1,
		batch.size());
	    const int mergedCount = live->mergeCompactOccurrences(batch);
	    if (mergedCount > 0) {
		merged = 1;
		ged_obol_timing_log("stream: merged occurrences",
		    bu_gettime(), (long)mergedCount);
	    }
	    /* More may have been queued while we drained; ask for another tick. */
	    if (has_more && item->stream->size() > 0)
		*has_more = 1;
	    const int64_t elapsed = bu_gettime() - started;
	    if ((max_items && drained_count >= max_items) ||
		(max_microseconds && elapsed >= 0 &&
		 static_cast<uint64_t>(elapsed) >= max_microseconds)) {
		if (has_more)
		    *has_more = 1;
		return merged;
	    }
	}
    }
    return merged;
}

static int
ged_obol_progressive_advance_provider(
    BObolViewController *controller,
    void *user_data,
    const BObolProgressiveOptions *options,
    BObolProgressiveStatus *status)
{
    /* Idle pumps are expected at frame cadence.  Logging thousands of 0.0 ms
     * samples obscures the stalls timing mode exists to expose and can itself
     * perturb interactive traces. */
    ged_obol_scoped_timer _pump_timer("provider: pump (total)", 1.0);
    BObolProgressiveStatus local_status;
    if (status)
	local_status.providerCount = status->providerCount;

    ged_obol_progressive_provider_data *data =
	static_cast<ged_obol_progressive_provider_data *>(user_data);
    if (!controller || !data || !data->gedp)
	return -1;

    struct ged_view_context *view_ctx = data->view_ctx;
    (void)options;

    if (data->pending_autoview) {
	const struct bv *view = ged_obol_bv_const(view_ctx);
	if (!view ||
	    bv_frame_revision_get(view) != data->expected_view_revision)
	    data->pending_autoview = 0;
    }

    /* There is one production progression: compact per-leaf boxes followed by
     * streamed view-appropriate geometry. */
    int has_pending_job = 0;
    ged_obol_cleanup_retired_jobs(data);

    /* Drain before retiring completed jobs.  A worker can finish between GUI
     * pumps with its final (or entire) occurrence set still queued.  Erasing
     * that job first discarded the stream whenever detached adoption was
     * rejected (for example after a harmless view-revision change), leaving
     * the live source permanently at its boxes. */
    int stream_more = 0;
    int refined =
	ged_obol_drain_streamed_realizations(data, controller,
	    options ? options->maxProviderItems : 0,
	    options ? options->maxProviderMicroseconds : 0, &stream_more) ?
	1 : 0;
    for (std::vector<std::shared_ptr<ged_obol_deferred_realization_job>>::iterator
	    it = data->pending_jobs.begin(); it != data->pending_jobs.end();) {
	    const std::shared_ptr<ged_obol_deferred_realization_job> job = *it;
	    const int job_state = job ?
		job->state.load(std::memory_order_acquire) :
		ged_obol_deferred_realization_job::FAILED;
	    if (job_state == ged_obol_deferred_realization_job::PENDING ||
		job_state == ged_obol_deferred_realization_job::RUNNING) {
		has_pending_job = 1;
		++it;
		continue;
	    }
	    if (job_state == ged_obol_deferred_realization_job::COMPLETE)
		refined += ged_obol_publish_deferred_realization(data,
			controller, job);
	    it = data->pending_jobs.erase(it);
	}
	if (data->deferred_refine_stage == 1) {
	    const int jobState = data->deferred_job ?
		data->deferred_job->state.load(std::memory_order_acquire) :
		ged_obol_deferred_realization_job::FAILED;
	    if (jobState == ged_obol_deferred_realization_job::PENDING ||
		jobState == ged_obol_deferred_realization_job::RUNNING) {
		has_pending_job = 1;
	    } else {
		if (jobState == ged_obol_deferred_realization_job::COMPLETE)
		    refined += ged_obol_publish_deferred_realization(data,
			controller, data->deferred_job);
		data->deferred_refine_stage = 3;
		data->deferred_paths.clear();
		data->deferred_job.reset();
	    }
	}

	/* A concurrent worker may have queued more after the pre-retirement
	 * drain.  Keep the provider alive for another pump in that case. */
	if (stream_more)
	    has_pending_job = 1;

    local_status.changed = refined > 0 ? 1 : 0;
    local_status.hasMore = has_pending_job;
    local_status.inFlight = has_pending_job ? 1 : 0;
    if (has_pending_job && ged_obol_timing_enabled()) {
	const int deferred_state = data->deferred_job ?
	    data->deferred_job->state.load(std::memory_order_acquire) : -1;
	bu_log("[obol-timing] provider pending: stage=%d deferred=%d "
	       "pending_jobs=%zu stream_more=%d\n",
	       data->deferred_refine_stage, deferred_state,
	       data->pending_jobs.size(), stream_more);
    }

    /* The deferred realization job streams leaf-local boxes/geometry onto the
     * root occurrence registry and atomically adopts the completed index. */
	if (local_status.changed && data->pending_autoview &&
	    ged_obol_progressive_autoview_apply(data))
	    controller->syncCameraFromViewContext(view_ctx, TRUE);
	if (!local_status.hasMore)
	    data->pending_autoview = 0;
	if (local_status.changed || local_status.hasMore)
	    controller->requestRender(local_status.changed ?
		"ged-deferred-full-detail" : "ged-deferred-root-building");
	if (status)
	    *status = local_status;
	return (local_status.changed || local_status.hasMore) ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_set_bounds_for_path(
    struct ged *gedp,
    const char *path,
    const point_t bmin,
    const point_t bmax)
{
    if (!bmin || !bmax)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int changed = scene->setDatabaseSourceInstanceBoundsState(
				source_instance_key.c_str(), TRUE,
				SbVec3f(static_cast<float>(bmin[0]),
					static_cast<float>(bmin[1]),
					static_cast<float>(bmin[2])),
				SbVec3f(static_cast<float>(bmax[0]),
					static_cast<float>(bmax[1]),
					static_cast<float>(bmax[2])));
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

extern "C" int
ged_draw_obol_database_source_set_display_name_for_path(
    struct ged *gedp,
    const char *path,
    const char *display_name)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int applied = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int changed = scene->setDatabaseSourceInstanceDisplayName(
				source_instance_key.c_str(),
				display_name ? display_name : "");
	if (changed < 0)
	    return 0;
	applied = 1;
    }
    return applied;
}

static const char *
ged_obol_realized_geometry_name(const BObolRealizedShapeSummary &summary)
{
    if (summary.shapeKind == BObolRealizedShapeSummary::SHAPE_MESH)
	return "indexed-face-set";
    if (summary.shapeKind == BObolRealizedShapeSummary::SHAPE_VLIST) {
	const char *kind = summary.geometryKind.getString();
	if (kind && BU_STR_EQUAL(kind, "annotation"))
	    return "annotation";
	if (kind && (BU_STR_EQUAL(kind, "point") ||
		     BU_STR_EQUAL(kind, "point-set")))
	    return "point-set";
	return "line-set";
    }
    return NULL;
}

static int
ged_obol_database_source_geometry_summary_for_source(
    SoBRLDatabaseSource *source,
    const char *fallback_path,
    struct ged_draw_shape_geometry_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!source)
	return 0;

    SoBRLVListShape *annotation_shape =
	ged_obol_owned_annotation_vlist_shape_for_source(source, fallback_path);
    const SoBRLVListShape *annotation_geom = annotation_shape ?
	annotation_shape->getGeometrySource() : NULL;
    if (ged_obol_vlist_shape_is_annotation(annotation_shape) &&
	(annotation_geom->point.getNum() > 0 ||
	 ged_obol_vlist_shape_has_annotation_record(annotation_shape))) {
	out->valid = 1;
	out->geometry_name = "annotation";
	out->point_count =
	    static_cast<size_t>(annotation_geom->point.getNum());
	out->index_count = 0;
	return 1;
    }

    if (source->hasCompactInstanceIndex()) {
	for (int i = 0; i < source->getCompactInstanceCount(); i++) {
	    BObolCompactOccurrence occurrence;
	    if (!source->getCompactOccurrence(i, occurrence) ||
		!occurrence.geometry)
		continue;
	    const char *geometryName =
		ged_obol_realized_geometry_name(occurrence.summary);
	    if (!geometryName)
		continue;
	    size_t pointCount = occurrence.summary.pointCount > 0 ?
		static_cast<size_t>(occurrence.summary.pointCount) : 0;
	    if (BU_STR_EQUAL(geometryName, "line-set") &&
		occurrence.geometry->wire) {
		pointCount = occurrence.geometry->wire->segmentPoints.size();
		for (const Obol::WirePolyline &polyline :
		     occurrence.geometry->wire->polylines)
		    pointCount += polyline.points.size();
	    }
	    if (!pointCount)
		continue;
	    out->valid = 1;
	    out->geometry_name = geometryName;
	    out->point_count = pointCount;
	    out->index_count = occurrence.summary.indexCount > 0 ?
		static_cast<size_t>(occurrence.summary.indexCount) : 0;
	    return 1;
	}
    }

    const int count = source->getRealizedShapeSummaryCount();
    struct ged_draw_shape_geometry_summary empty_summary;
    memset(&empty_summary, 0, sizeof(empty_summary));
    int have_empty_summary = 0;
    for (int i = 0; i < count; i++) {
	BObolRealizedShapeSummary summary;
	if (!source->getRealizedShapeSummary(i, summary) || !summary.valid)
	    continue;

	const char *geometry_name = ged_obol_realized_geometry_name(summary);
	if (!geometry_name)
	    continue;

	struct ged_draw_shape_geometry_summary current_summary;
	memset(&current_summary, 0, sizeof(current_summary));
	current_summary.valid = 1;
	current_summary.geometry_name = geometry_name;
	current_summary.point_count = (summary.pointCount > 0) ?
				      static_cast<size_t>(summary.pointCount) : 0;
	current_summary.index_count = (summary.indexCount > 0) ?
				      static_cast<size_t>(summary.indexCount) : 0;
	if (current_summary.point_count || current_summary.index_count) {
	    *out = current_summary;
	    return 1;
	}
	if (!have_empty_summary) {
	    empty_summary = current_summary;
	    have_empty_summary = 1;
	}
    }

    if (have_empty_summary) {
	*out = empty_summary;
	return 1;
    }

    out->valid = 1;
    out->geometry_name = "empty";
    return 1;
}

extern "C" int
ged_draw_obol_database_source_geometry_summary_for_path_mode(
    struct ged *gedp,
    const char *path,
    int draw_mode_valid,
    int draw_mode,
    struct ged_draw_shape_geometry_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int mode = draw_mode_valid ? draw_mode : -1;
    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, mode, 0);
    if (instance_keys.empty())
	return 0;

    struct ged_draw_shape_geometry_summary fallback_summary;
    memset(&fallback_summary, 0, sizeof(fallback_summary));
    int have_fallback = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	struct ged_draw_shape_geometry_summary current_summary;
	memset(&current_summary, 0, sizeof(current_summary));
	if (!ged_obol_database_source_geometry_summary_for_source(source,
		path, &current_summary) || !current_summary.valid)
	    continue;
	if ((current_summary.point_count || current_summary.index_count) &&
	    current_summary.geometry_name) {
	    *out = current_summary;
	    return 1;
	}
	if (!have_fallback) {
	    fallback_summary = current_summary;
	    have_fallback = 1;
	}
    }

    if (have_fallback) {
	*out = fallback_summary;
	return 1;
    }
    return 0;
}

extern "C" int
ged_draw_obol_database_source_geometry_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_shape_geometry_summary *out)
{
    return ged_draw_obol_database_source_geometry_summary_for_path_mode(gedp,
	    path, 0, GED_DRAW_MODE_WIRE, out);
}

extern "C" int
ged_draw_obol_database_source_material_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_shape_material_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    out->valid = 1;
    out->material_revision = summary.materialRevision;
    ged_obol_rgb_from_color(ged_obol_summary_material_color(summary),
			    out->material_color);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_refresh_material_color_for_path(
    struct ged *gedp,
    const char *path,
    struct db_i *dbip,
    uint64_t material_revision)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int changed_any = 0;
    for (const std::string &source_instance_key : instance_keys) {
	const int changed =
	    scene->refreshDatabaseSourceInstanceMaterialColorFromDatabase(
		source_instance_key.c_str(),
		ged_obol_fold_revision(material_revision),
		dbip);
	if (changed < 0)
	    return 0;
	if (changed > 0)
	    changed_any = 1;
    }
    return changed_any;
}

extern "C" int
ged_draw_obol_database_sources_refresh_material_colors(
    struct ged *gedp,
    struct db_i *dbip,
    uint64_t material_revision)
{
    if (!gedp)
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed =
	scene->refreshDatabaseSourceMaterialColorsFromDatabase(
	    ged_obol_fold_revision(material_revision), dbip);
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_evaluated_region_for_path(
    struct ged *gedp,
    const char *path,
    int *out)
{
    if (out)
	*out = 0;
    if (!out)
	return 0;

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    const int count = source->getRealizedShapeSummaryCount();
    for (int i = 0; i < count; i++) {
	BObolRealizedShapeSummary summary;
	if (!source->getRealizedShapeSummary(i, summary) || !summary.valid)
	    continue;
	*out = summary.regionId ? 1 : 0;
	return 1;
    }

    return 0;
}

extern "C" int
ged_draw_obol_database_source_set_evaluated_region_for_path(
    struct ged *gedp,
    const char *path,
    int evaluated_region)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    const int next_region = evaluated_region ? 1 : 0;
    int updated = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;
	if (source->hasCompactInstanceIndex()) {
	    updated |= source->setCompactInstanceRegionIdForPath(path, TRUE,
		next_region);
	    continue;
	}

	const int vlist_count = source->getRealizedShapeCount();
	for (int i = 0; i < vlist_count; i++) {
	    SoBRLVListShape *shape = source->getRealizedShape(i);
	    if (!shape)
		continue;
	    shape->regionId = next_region;
	    updated = 1;
	}

	SoBRLMeshShape *mesh = source->getRealizedMesh();
	if (mesh) {
	    mesh->regionId = next_region;
	    updated = 1;
	}
    }

    return updated;
}

extern "C" int
ged_draw_obol_database_source_set_region_metadata_for_path(
    struct ged *gedp,
    const char *path,
    int region_id,
    int aircode,
    int los,
    int material_id)
{
    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_database_source_instance_keys_for_path(gedp, path, -1, 0);
    if (instance_keys.empty())
	return 0;

    int updated = 0;
    for (const std::string &source_instance_key : instance_keys) {
	SoBRLDatabaseSource *source =
	    ged_obol_database_source_for_instance_key(scene,
		source_instance_key);
	if (!source)
	    return 0;
	if (source->hasCompactInstanceIndex()) {
	    updated |= source->setCompactInstanceRegionMetadataForPath(path,
		TRUE, region_id, aircode, material_id, los);
	    continue;
	}

	const int vlist_count = source->getRealizedShapeCount();
	for (int i = 0; i < vlist_count; i++) {
	    SoBRLVListShape *shape = source->getRealizedShape(i);
	    if (!shape)
		continue;
	    shape->regionId = region_id;
	    shape->airCode = aircode;
	    shape->los = los;
	    shape->materialId = material_id;
	    updated = 1;
	}

	SoBRLMeshShape *mesh = source->getRealizedMesh();
	if (mesh) {
	    mesh->regionId = region_id;
	    mesh->airCode = aircode;
	    mesh->los = los;
	    mesh->materialId = material_id;
	    updated = 1;
	}
    }

    return updated;
}

template <typename ShapeT>
static void
ged_obol_apply_draw_metadata_to_shape(
    ShapeT *shape,
    const SoBRLDatabaseSource *source,
    const struct BObolDrawMetadataRecord *record)
{
    if (!shape || !record)
	return;

    shape->regionId = record->hasRegionId ? record->regionId : 0;
    shape->airCode = record->hasAircode ? record->aircode : 0;
    shape->los = record->hasLos ? record->los : 0;
    shape->materialId = record->hasMaterialId ? record->materialId : 0;

    shape->materialColorValid = record->hasColor ? TRUE : FALSE;
    if (record->hasColor) {
	shape->materialColor = SbColor(
				   static_cast<float>(record->color[0]) / 255.0f,
				   static_cast<float>(record->color[1]) / 255.0f,
				   static_cast<float>(record->color[2]) / 255.0f);
    } else {
	shape->materialColor = SbColor(1.0f, 1.0f, 1.0f);
    }

    if (source && source->materialColorValid.getValue() &&
	(source->materialPolicy.getValue() !=
	 SoBRLDatabaseSource::MATERIAL_DATABASE ||
	 !shape->materialColorValid.getValue())) {
	shape->materialColorValid = TRUE;
	shape->materialColor = source->materialColor.getValue();
	shape->materialRevision = source->materialRevision.getValue();
    }

    shape->materialShader = record->hasShader ? record->shader : "";
}

static int
ged_obol_apply_draw_metadata_to_source(
    SoBRLDatabaseSource *source,
    const struct BObolDrawMetadataRecord *record)
{
    if (!source || !record || !record->directoryFound)
	return 0;

    SbColor metadata_color(1.0f, 1.0f, 1.0f);
    if (record->hasColor) {
	metadata_color = SbColor(
			     static_cast<float>(record->color[0]) / 255.0f,
			     static_cast<float>(record->color[1]) / 255.0f,
			     static_cast<float>(record->color[2]) / 255.0f);
    }
    (void)source->setDatabaseMetadataState(
	TRUE,
	record->hasRegionId ? record->regionId : 0,
	record->hasAircode ? record->aircode : 0,
	record->hasMaterialId ? record->materialId : 0,
	record->hasLos ? record->los : 0,
	record->hasColor ? TRUE : FALSE,
	metadata_color,
	record->hasShader ? SbString(record->shader) : SbString(""));
    return 1;
}

extern "C" int
ged_draw_obol_database_source_apply_draw_metadata_for_path(
    struct ged *gedp,
    const char *path,
    const struct BObolDrawMetadataRecord *record)
{
    if (!record || !record->directoryFound)
	return 0;

    BObolSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    SoBRLDatabaseSource *source = scene ?
				  scene->findDatabaseSourceInstance(source_instance_key.c_str()) : NULL;
    if (!source)
	return 0;

    if (source->hasCompactInstanceIndex()) {
	SbColor metadata_color(1.0f, 1.0f, 1.0f);
	if (record->hasColor) {
	    metadata_color = SbColor(
		static_cast<float>(record->color[0]) / 255.0f,
		static_cast<float>(record->color[1]) / 255.0f,
		static_cast<float>(record->color[2]) / 255.0f);
	}
	return source->setCompactInstanceMetadataForPath(path, TRUE,
	    record->hasRegionId ? record->regionId : 0,
	    record->hasAircode ? record->aircode : 0,
	    record->hasMaterialId ? record->materialId : 0,
	    record->hasLos ? record->los : 0,
	    record->hasColor ? TRUE : FALSE, metadata_color,
	    record->hasShader ? SbString(record->shader) : SbString(""));
    }

    int updated = ged_obol_apply_draw_metadata_to_source(source, record);

    const int vlist_count = source->getRealizedShapeCount();
    for (int i = 0; i < vlist_count; i++) {
	SoBRLVListShape *shape = source->getRealizedShape(i);
	if (!shape)
	    continue;
	ged_obol_apply_draw_metadata_to_shape(shape, source, record);
	updated = 1;
    }

    const int mesh_count = source->getRealizedMeshCount();
    for (int i = 0; i < mesh_count; i++) {
	SoBRLMeshShape *mesh = source->getRealizedMesh(i);
	if (!mesh)
	    continue;
	ged_obol_apply_draw_metadata_to_shape(mesh, source, record);
	updated = 1;
    }

    return updated;
}

static int
ged_obol_draw_mode_uses_database_source(int ged_mode)
{
    switch (ged_mode) {
	case GED_DRAW_MODE_WIRE:
	case GED_DRAW_MODE_SHADED_BOTS:
	case GED_DRAW_MODE_SHADED:
	case GED_DRAW_MODE_HIDDEN_LINE:
	case GED_DRAW_MODE_EVAL_WIRE:
	case GED_DRAW_MODE_EVAL_POINTS:
	    return 1;
	default:
	    return 0;
    }
}

static int
ged_obol_draw_mode_uses_root_source(int ged_mode)
{
    return ged_mode == GED_DRAW_MODE_WIRE ||
	   ged_mode == GED_DRAW_MODE_SHADED_BOTS ||
	   ged_mode == GED_DRAW_MODE_SHADED ||
	   ged_mode == GED_DRAW_MODE_HIDDEN_LINE ||
	   ged_mode == GED_DRAW_MODE_EVAL_WIRE ||
	   ged_mode == GED_DRAW_MODE_EVAL_POINTS;
}

static int
ged_obol_draw_mode_publishes_direct_face_set(int ged_mode)
{
    return ged_mode == GED_DRAW_MODE_SHADED_BOTS ||
	   ged_mode == GED_DRAW_MODE_SHADED ||
	    ged_mode == GED_DRAW_MODE_HIDDEN_LINE;
}

static int
ged_obol_existing_path_is_nondrawable(struct db_i *dbip, const char *path)
{
    if (!dbip || !path || !path[0])
	return 0;

    struct db_full_path dfp;
    db_full_path_init(&dfp);
    int valid = (db_string_to_path(&dfp, dbip, path) == 0);
    if (!valid) {
	db_free_full_path(&dfp);
	return 0;
    }

    struct directory *dp = DB_FULL_PATH_CUR_DIR(&dfp);
    int nondrawable = (!dp ||
	    !(dp->d_flags & (RT_DIR_SOLID | RT_DIR_COMB))) ? 1 : 0;
    db_free_full_path(&dfp);
    return nondrawable;
}

static union tree *
ged_obol_direct_draw_nop_tree(void)
{
    union tree *curtree;

    BU_GET(curtree, union tree);
    RT_TREE_INIT(curtree);
    curtree->tr_op = OP_NOP;
    return curtree;
}

static std::string
ged_obol_direct_leaf_instance_key(const struct db_tree_state *tsp,
				  const struct db_full_path *pathp)
{
    std::string key;
    const uint32_t key_version = 5;
    const int sofar = tsp ?
		      (tsp->ts_sofar &
		       (TS_SOFAR_MINUS | TS_SOFAR_INTER | TS_SOFAR_REGION)) : 0;
    mat_t mat;

    if (!pathp)
	return key;

    RT_CK_FULL_PATH(pathp);
    if (tsp)
	MAT_COPY(mat, tsp->ts_mat);
    else
	MAT_IDN(mat);

    auto append_bytes = [&key](const void *ptr, size_t len) {
	key.append(static_cast<const char *>(ptr), len);
    };

    append_bytes(&key_version, sizeof(key_version));
    append_bytes(&pathp->fp_len, sizeof(pathp->fp_len));
    append_bytes(&sofar, sizeof(sofar));
    append_bytes(mat, sizeof(mat));
    for (size_t i = 0; i < pathp->fp_len; i++) {
	const struct directory *dp = pathp->fp_names ?
				     pathp->fp_names[i] : NULL;
	const char *name = (dp && dp->d_namep) ? dp->d_namep : "";
	const int cinst = pathp->fp_cinst ?
			  DB_FULL_PATH_GET_COMB_INST(pathp, i) : 0;
	const size_t name_len = strlen(name) + 1;
	append_bytes(&name_len, sizeof(name_len));
	append_bytes(name, name_len);
	append_bytes(&cinst, sizeof(cinst));
    }

    return key;
}

static std::string
ged_obol_direct_leaf_source_instance_key(struct ged_view_context *view_ctx,
					int ged_mode,
					const struct db_tree_state *tsp,
					const struct db_full_path *pathp)
{
    const std::string identity =
	ged_obol_direct_leaf_instance_key(tsp, pathp);
    if (identity.empty())
	return std::string();

    bu_h128_t hval = bu_data_hash128(identity.data(), identity.size());
    char hash_buf[128] = {0};
    snprintf(hash_buf, sizeof(hash_buf),
	     "%016" PRIx64 "%016" PRIx64,
	     hval.w[1], hval.w[0]);

    std::string key;
    if (ged_obol_view_scope_is_independent(view_ctx)) {
	key = "ged-view:";
	key += ged_obol_view_scope_name(view_ctx);
	key += ":";
    }
    key += "brlcad-direct:";
    key += hash_buf;
    if (ged_mode >= 0) {
	char mode_buf[64] = {0};
	snprintf(mode_buf, sizeof(mode_buf), "%s%d",
		 ged_obol_database_source_mode_key_marker, ged_mode);
	key += mode_buf;
    }
    return key;
}

struct ged_obol_direct_draw_ctx {
    struct ged *gedp;
    struct ged_view_context *view_ctx;
    BObolSceneController *scene;
    const struct ged_draw_appearance_settings *settings;
    uint32_t source_revision;
    std::string target_group_path;
    std::unordered_set<std::string> seen_instances;
    int realized;
    int failed;
};

static int
ged_obol_direct_apply_leaf_state(struct ged_obol_direct_draw_ctx *ctx,
				 const char *path,
				 const char *source_instance_key,
				 const struct BObolDrawMetadataRecord *metadata,
				 int line_style)
{
    if (!ctx || !ctx->scene || !path || !path[0] ||
	!source_instance_key || !source_instance_key[0])
	return 0;

    SoBRLDatabaseSource *source =
	ctx->scene->findDatabaseSourceInstance(source_instance_key);
    if (!source)
	return 0;

    if (metadata && metadata->directoryFound)
	(void)ged_obol_apply_draw_metadata_to_source(source, metadata);

    BObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    if (summary.lineStyle == line_style)
	return 1;

    const int changed = ctx->scene->setDatabaseSourceInstanceState(
			    source_instance_key,
			    FALSE,
			    summary.sourceRevision,
			    summary.inputsRevision,
			    summary.visible,
			    summary.selected,
			    summary.highlighted,
			    line_style,
			    summary.lineWidth,
			    summary.transparency,
			    summary.colorOverride,
			    summary.color,
			    summary.materialColorValid,
			    summary.materialColor,
			    summary.materialRevision);
    return changed >= 0 ? 1 : 0;
}

static union tree *
ged_obol_direct_draw_leaf_cb(struct db_tree_state *tsp,
			     const struct db_full_path *pathp,
			     struct directory *dp,
			     void *client_data)
{
    struct ged_obol_direct_draw_ctx *ctx =
	static_cast<struct ged_obol_direct_draw_ctx *>(client_data);
    if (!ctx || !ctx->gedp || !ctx->scene || !ctx->settings || !pathp ||
	!dp || ctx->failed)
	return TREE_NULL;

    const std::string instance_key =
	ged_obol_direct_leaf_source_instance_key(ctx->view_ctx,
	    ctx->settings->draw_mode, tsp, pathp);
    if (instance_key.empty()) {
	ctx->failed = 1;
	return TREE_NULL;
    }
    if (!ctx->seen_instances.insert(instance_key).second)
	return ged_obol_direct_draw_nop_tree();

    char *path = db_path_to_string(pathp);
    if (!path || !path[0]) {
	if (path)
	    bu_free(path, "Obol direct draw leaf source path");
	ctx->failed = 1;
	return TREE_NULL;
    }

    struct BObolDrawMetadataRecord metadata;
    bobol_draw_metadata_record_from_tree_state(&metadata, tsp, dp);

    ged_obol_publish_placement_state placement;
    if (tsp) {
	placement.valid = true;
	placement.drawMatrixValid = true;
	placement.drawMatrix = ged_obol_sbmatrix_from_mat(tsp->ts_mat);
    }

    const int changed = ged_obol_replace_path(ctx->gedp, ctx->view_ctx,
			    ctx->gedp->dbip, path, ctx->settings->draw_mode,
			    ctx->source_revision, ctx->scene, 0, 1,
			    ctx->settings, placement.valid ? &placement : NULL,
			    ctx->target_group_path.c_str(),
			    instance_key.c_str(), &metadata);
    if (changed < 0) {
	bu_free(path, "Obol direct draw leaf source path");
	ctx->failed = 1;
	return TREE_NULL;
    }

    if (ged_draw_test_primitive_face_set_failure_enabled() &&
	ged_obol_draw_mode_publishes_direct_face_set(ctx->settings->draw_mode)) {
	(void)ged_draw_obol_database_source_publish_indexed_face_set_for_path(
	    ctx->gedp, path, NULL, 0, NULL, 0, NULL, 0);
    }

    const int subtract_style =
	(ctx->settings->draw_solid_lines_only || !tsp ||
	 !(tsp->ts_sofar & TS_SOFAR_MINUS)) ? 0 : 1;
    if (!ged_obol_direct_apply_leaf_state(ctx, path,
	    instance_key.c_str(),
	    metadata.directoryFound ? &metadata : NULL, subtract_style)) {
	bu_free(path, "Obol direct draw leaf source path");
	ctx->failed = 1;
	return TREE_NULL;
    }

    bu_free(path, "Obol direct draw leaf source path");
    ctx->realized++;
    return ged_obol_direct_draw_nop_tree();
}

static int
ged_obol_direct_draw_root_source(struct ged *gedp,
				 struct ged_view_context *view_ctx,
				 BObolSceneController *scene,
				 const char *path,
				 const struct ged_draw_appearance_settings *settings,
				 uint32_t source_revision)
{
    if (!gedp || !gedp->dbip || !scene || !path || !path[0] || !settings)
	return 0;

    (void)ged_draw_obol_group_ensure_for_path(gedp, path, path,
	    settings->draw_mode, 0);
    const std::string group_path = ged_obol_group_path_from_record_path(path);
    struct BObolDrawMetadataRecord root_metadata;
    const struct BObolDrawMetadataRecord *metadata = NULL;
    if (settings->draw_mode == GED_DRAW_MODE_WIRE ||
	settings->draw_mode == GED_DRAW_MODE_SHADED_BOTS ||
	settings->draw_mode == GED_DRAW_MODE_SHADED ||
	settings->draw_mode == GED_DRAW_MODE_HIDDEN_LINE) {
	bobol_draw_metadata_record_init(&root_metadata);
	root_metadata.directoryFound = 1;
	metadata = &root_metadata;
    }
    const int changed = ged_obol_replace_path(gedp, view_ctx, gedp->dbip,
			    path, settings->draw_mode, source_revision, scene,
			    0, 1, settings, NULL, group_path.c_str(), NULL,
			    metadata);
    if (changed < 0)
	return 0;

    return 1;
}

static int ged_obol_apply_cached_path_metadata_to_scene(
    struct ged *gedp,
    BObolSceneController *scene,
    struct ged_view_context *view_ctx,
    const char *path,
    int draw_mode);

static int
ged_obol_apply_draw_paths_to_scene(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char **paths,
    int path_count,
    const struct ged_draw_appearance_settings *settings,
    struct ged_draw_transaction_result *result,
    BObolSceneController *scene,
    uint32_t source_revision,
    int preserve_display_state)
{
    if (!gedp || !gedp->dbip || !view_ctx || !paths || path_count <= 0 ||
	!settings || !scene || !ged_obol_draw_mode_uses_database_source(
	    settings->draw_mode))
	return -1;
    ged_obol_scoped_timer _draw_timer("draw: apply_draw_paths (total)");

    std::vector<std::string> draw_paths;
    draw_paths.reserve(static_cast<size_t>(path_count));
    for (int i = 0; i < path_count; i++) {
	if (paths[i] && paths[i][0])
	    ged_obol_append_unique_path(draw_paths, paths[i]);
    }
    if (draw_paths.empty())
	return 0;

    if (!source_revision)
	source_revision = ged_obol_fold_revision(ged_draw_scene_revision(gedp));
    int changed = 0;
    int realized_roots = 0;
    int realized_sources = 0;
    ged_obol_scene_mutation_batch_scope batch(scene, draw_paths.size(),
	    draw_paths.size());
    std::vector<ged_obol_preserved_source_display_state>
	preserved_display_states;
    if (preserve_display_state) {
	preserved_display_states = ged_obol_preserve_source_display_states(
	    scene, view_ctx, draw_paths, settings->draw_mode);
    }

    const int remove_mode = settings->mixed_modes ? settings->draw_mode : -1;
    const bool retained_root_mode = settings->defer_leaf_expansion ||
	ged_obol_draw_mode_uses_root_source(settings->draw_mode);
    if (!retained_root_mode)
	(void)ged_obol_remove_paths(draw_paths, view_ctx, scene, remove_mode);

    auto remove_obsolete_root_sources = [&]() {
	std::set<std::string> expected;
	std::map<std::string, SoBRLDatabaseSource *> broad_sources;
	for (const std::string &path : draw_paths)
	{
	    const std::string key =
		ged_obol_database_source_instance_key_for_mode(
		    view_ctx, path.c_str(), settings->draw_mode);
	    expected.insert(key);
	    SoBRLDatabaseSource *broad = scene->findDatabaseSourceInstance(
		key.c_str());
	    if (broad)
		broad_sources[path] = broad;
	}
	std::vector<std::string> candidates =
	    ged_obol_matching_database_source_instance_keys(scene, view_ctx,
		draw_paths, 0, 1, remove_mode);
	std::vector<std::string> obsolete;
	for (const std::string &key : candidates) {
	    if (expected.find(key) != expected.end())
		continue;

	    /* A broader draw replaces a narrower source, but it must not regress
	     * the visible data to a root proxy while progressive realization
	     * restarts.  Share every resident immutable occurrence with the broad
	     * owner before retiring the narrow presentation.  The broad source's
	     * normal realization then fills only missing/richer data. */
	    SoBRLDatabaseSource *narrow =
		scene->findDatabaseSourceInstance(key.c_str());
	    BObolDatabaseSourceSummary narrow_summary;
	    if (narrow && narrow->getSummary(narrow_summary) &&
		narrow_summary.valid &&
		ged_obol_database_source_summary_matches_mode(narrow_summary,
		    settings->draw_mode)) {
		const char *narrow_path = narrow_summary.path.getString();
		for (const auto &broad_entry : broad_sources) {
		    if (!narrow_path || !narrow_path[0] ||
			ged_obol_path_equal(narrow_path,
			    broad_entry.first.c_str()))
			continue;
		    if (!ged_obol_path_has_prefix(narrow_path,
			    broad_entry.first.c_str()) &&
			!ged_obol_path_has_semantic_prefix(narrow_path,
			    broad_entry.first.c_str()))
			continue;
		    (void)broad_entry.second->adoptCompactOccurrencesFrom(narrow);
		    break;
		}
	    }
	    obsolete.push_back(key);
	}
	return ged_obol_remove_instance_keys(obsolete, scene);
    };

    if (settings->defer_leaf_expansion) {
	for (const std::string &path : draw_paths) {
	    const int64_t root_start = bu_gettime();
	    if (ged_obol_direct_draw_root_source(gedp, view_ctx, scene,
		    path.c_str(), settings, source_revision)) {
		changed = 1;
		realized_roots++;
		realized_sources++;
	    }
	    ged_obol_timing_log("draw: direct_draw_root_source", root_start, -1);
	    const std::string instance_key =
		ged_obol_database_source_instance_key_for_mode(view_ctx,
		    path.c_str(), settings->draw_mode);
	    const int64_t proxy_start = bu_gettime();
	    if (!instance_key.empty() &&
		ged_obol_publish_deferred_root_proxy(gedp, view_ctx,
		    path.c_str(), instance_key.c_str(), settings->draw_mode))
		changed = 1;
	    ged_obol_timing_log("draw: publish_root_proxy (coarse)",
		proxy_start, -1);
	    if (ged_obol_apply_cached_path_metadata_to_scene(gedp, scene,
		    view_ctx, path.c_str(), settings->draw_mode))
		changed = 1;
	}
	if (changed)
	    changed |= ged_obol_restore_source_display_states(scene,
		preserved_display_states);
	/* A deferred root is still the sole representation for its draw path.
	 * Normal draw replaces every existing representation; only --add-mode
	 * retains another mode alongside this root. */
	if (remove_obsolete_root_sources())
	    changed = 1;
	/* Deferred roots own a structural proxy frontier.  Its descendants are
	 * not obsolete roots; full-detail adoption retires them atomically. */
	if (result) {
	    result->affected_groups = realized_roots;
	    result->affected_shapes = realized_sources;
	}
	return 1;
    }

    if (ged_obol_draw_mode_uses_root_source(settings->draw_mode)) {
	for (const std::string &path : draw_paths) {
	    if (ged_obol_direct_draw_root_source(gedp, view_ctx, scene,
		    path.c_str(), settings, source_revision)) {
		changed = 1;
		realized_roots++;
		realized_sources++;
	    }
	}
	if (changed)
	    changed |= ged_obol_restore_source_display_states(scene,
		preserved_display_states);
	if (remove_obsolete_root_sources())
	    changed = 1;
	if (changed) {
	    const SbBool realized = scene->realizePending();
	    if (!realized && settings->strict_fallback)
		return -1;
	    if (!realized && !settings->strict_fallback &&
		ged_obol_draw_mode_publishes_direct_face_set(
		    settings->draw_mode)) {
		int fallbackRealized = 0;
		for (const std::string &path : draw_paths) {
		    const std::string instance_key =
			ged_obol_database_source_instance_key_for_mode(
			    view_ctx, path.c_str(), settings->draw_mode);
		    SoBRLDatabaseSource *source =
			scene->findDatabaseSourceInstance(instance_key.c_str());
		    if (source && source->realizationStatus.getValue() ==
			    SoBRLDatabaseSource::FAILED &&
			source->realizeDatabaseWireframe())
			fallbackRealized++;
		}
		if (fallbackRealized > 0)
		    (void)scene->realizePending();
	    }
	}
	for (const std::string &path : draw_paths) {
	    const std::string instance_key =
		ged_obol_database_source_instance_key_for_mode(view_ctx,
		    path.c_str(), settings->draw_mode);
	    SoBRLDatabaseSource *source =
		scene->findDatabaseSourceInstance(instance_key.c_str());
	    if (source)
		(void)source->setCompactSubtractLineStyle(
		    settings->draw_solid_lines_only ? 0 : 1);
	}
	if (result) {
	    result->affected_groups = realized_roots;
	    result->affected_shapes = realized_sources;
	}
	return 1;
    }

    int old_use_comb_instance_ids = db_comb_instance_ids_set(gedp->dbip, 1);
    if (old_use_comb_instance_ids < 0)
	return 0;
    for (const std::string &path : draw_paths) {
	(void)ged_draw_obol_group_ensure_for_path(gedp, path.c_str(),
		path.c_str(), settings->draw_mode, 0);

	struct db_tree_state init_state;
	db_init_db_tree_state(&init_state, gedp->dbip);
	init_state.ts_stop_at_regions = 0;

	struct ged_obol_direct_draw_ctx ctx;
	ctx.gedp = gedp;
	ctx.view_ctx = view_ctx;
	ctx.scene = scene;
	ctx.settings = settings;
	ctx.source_revision = source_revision;
	ctx.target_group_path = ged_obol_group_path_from_record_path(
	    path.c_str());
	ctx.realized = 0;
	ctx.failed = 0;

	const char *av[1] = { ged_draw_dbpath_skip_lead_slash(path.c_str()) };
	int walk_ret = db_walk_tree_leaf_instances(gedp->dbip, 1, av, 1,
			   &init_state, NULL, NULL,
			   ged_obol_direct_draw_leaf_cb,
			   static_cast<void *>(&ctx));
	db_free_db_tree_state(&init_state);

	if (walk_ret < 0 || ctx.failed) {
	    if (walk_ret < 0 && !ctx.failed &&
		    ged_obol_existing_path_is_nondrawable(gedp->dbip,
			path.c_str()))
		continue;
	    (void)db_comb_instance_ids_set(gedp->dbip,
		old_use_comb_instance_ids);
	    return -1;
	}
	if (ctx.realized > 0) {
	    changed = 1;
	    realized_roots++;
	    realized_sources += ctx.realized;
	} else if (ged_obol_direct_draw_root_source(gedp, view_ctx, scene,
		path.c_str(), settings, source_revision)) {
	    changed = 1;
	    realized_roots++;
	    realized_sources++;
	}
    }
    (void)db_comb_instance_ids_set(gedp->dbip,
	old_use_comb_instance_ids);

    if (changed)
	changed |= ged_obol_restore_source_display_states(scene,
	    preserved_display_states);
    if (changed)
	scene->realizePending();

    if (result) {
	result->affected_groups = realized_roots;
	result->affected_shapes = realized_sources;
    }

    return 1;
}

extern "C" int
ged_draw_obol_apply_draw_paths(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char **paths,
    int path_count,
    const struct ged_draw_appearance_settings *settings,
    struct ged_draw_transaction_result *result)
{
    BObolSceneController *scene = gedp ? ged_draw_obol_scene_controller(gedp) : NULL;
    return ged_obol_apply_draw_paths_to_scene(gedp, view_ctx, paths, path_count,
	    settings, result, scene, 0, 0);
}

static int
ged_obol_scene_sync_full_scene_impl(struct ged *gedp,
				    struct ged_view_context *view_ctx,
				    uint32_t source_revision,
				    BObolSceneController *controller,
				    int preserve_existing_revision)
{
    BObolSceneController *scene = controller ?
				  controller : ged_draw_obol_scene_controller(gedp);
    if (!gedp || !scene)
	return 0;

    if (!source_revision)
	source_revision = ged_obol_fold_revision(ged_draw_scene_revision(gedp));

    struct db_i *dbip = gedp->dbip;
    std::vector<ged_obol_drawn_source_path_mode> path_modes =
	ged_obol_drawn_source_path_modes(gedp, view_ctx, -1, NULL);

    if (!dbip || path_modes.empty()) {
	int changed = ged_obol_clear_database_sources_in_scope(scene, view_ctx);
	if (changed)
	    scene->realizePending();
	return changed ? 1 : 0;
    }

    std::set<std::string> expected_keys;
    for (const ged_obol_drawn_source_path_mode &entry : path_modes) {
	expected_keys.insert(ged_obol_database_source_instance_key_for_mode(
	    view_ctx, entry.path.c_str(), entry.mode));
    }

    int changed = 0;
    if (ged_obol_replace_path_modes(dbip, path_modes, source_revision,
				    gedp, view_ctx, scene, 1,
				    preserve_existing_revision))
	changed = 1;

    std::vector<std::string> obsolete_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx))
	    continue;
	const std::string key = summary.instanceKey.getString();
	if (!key.empty() && expected_keys.find(key) == expected_keys.end())
	    obsolete_keys.push_back(key);
    }
    if (ged_obol_remove_instance_keys(obsolete_keys, scene))
	changed = 1;

    if (changed)
	scene->realizePending();
    return changed ? 1 : 0;
}

int
ged_draw_obol_scene_sync_full_scene(struct ged *gedp,
				    struct ged_view_context *view_ctx,
				    uint32_t source_revision,
				    BObolSceneController *controller)
{
    return ged_obol_scene_sync_full_scene_impl(gedp, view_ctx,
	    source_revision, controller, 0);
}

static int
ged_obol_apply_source_update_transaction(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    uint32_t source_revision,
    BObolSceneController *scene)
{
    if (!gedp || !txn || !scene)
	return 0;

    std::vector<std::string> targets;
    ged_obol_append_unique_path(targets, txn->path);
    for (int i = 0; txn->paths && i < txn->path_count; i++)
	ged_obol_append_unique_path(targets, txn->paths[i]);
    if (targets.empty())
	targets = ged_obol_transaction_paths(txn, result);
    if (targets.empty())
	return 0;

    if (!source_revision)
	source_revision = ged_obol_fold_revision(ged_draw_scene_revision(gedp));

    if (txn->removed) {
	std::vector<std::string> source_paths =
	    ged_obol_matching_database_source_paths(scene, view_ctx, targets,
		0, 1, txn->mode);
	if (source_paths.empty())
	    return (result && result->status > 0) ? 1 : 0;
	(void)ged_obol_remove_paths(source_paths, view_ctx, scene,
				    txn->mode);
	return 1;
    }

    if (txn->redraw) {
	int handled = ged_obol_replace_matching_database_sources(gedp,
		      view_ctx, targets, 0, 1, source_revision, scene,
		      txn->mode);
	if (handled)
	    return 1;
    }

    const ged_draw_stale_reason stale_reason = txn->stale_reason ?
	txn->stale_reason : GED_DRAW_STALE_SOURCE_CHANGED;
    return ged_obol_mark_matching_database_sources_stale(view_ctx, targets,
	    0, 1,
	    ged_obol_stale_reason_from_ged(stale_reason), scene,
	    txn->mode);
}

static int
ged_obol_apply_source_references_removed_transaction(
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    BObolSceneController *scene)
{
    if (!txn || !scene)
	return 0;

    std::vector<std::string> targets =
	ged_obol_transaction_paths(txn, result);
    if (targets.empty())
	return 0;

    std::vector<std::string> source_paths =
	ged_obol_matching_database_source_paths(scene, txn->view, targets, 1,
	    0, txn->mode);
    if (source_paths.empty())
	return (result && result->status > 0) ? 1 : 0;
    (void)ged_obol_remove_paths(source_paths, txn->view, scene, txn->mode);
    return 1;
}

static int
ged_obol_apply_stale_source_transaction(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    uint32_t source_revision,
    BObolSceneController *scene)
{
    if (!txn || !scene)
	return 0;

    std::vector<std::string> targets =
	ged_obol_transaction_paths(txn, result);
    if (targets.empty())
	return (result && result->status >= 0) ? 1 : 0;

    const ged_draw_stale_reason stale_reason = txn->stale_reason ?
	txn->stale_reason : GED_DRAW_STALE_SOURCE_CHANGED;
    if (ged_obol_stale_reason_from_ged(stale_reason) ==
	SoBRLDatabaseSource::STALE_SOURCE &&
	ged_obol_refresh_matching_compact_parts(scene, view_ctx, targets,
	    txn->mode, source_revision) > 0)
	return 1;

    if (view_ctx) {
	if (!source_revision)
	    source_revision =
		ged_obol_fold_revision(ged_draw_scene_revision(gedp));
	int refreshed = ged_obol_replace_matching_database_sources(gedp,
			view_ctx, targets, 0, 1, source_revision, scene,
			txn->mode);
	if (refreshed)
	    return 1;
    }

    int handled = ged_obol_mark_matching_database_sources_stale(view_ctx,
		  targets, 0, 1,
		  ged_obol_stale_reason_from_ged(stale_reason), scene,
		  txn->mode);
    if (!handled && result && result->status >= 0)
	return 1;
    return handled;
}

static int
ged_obol_remove_groups_by_path_prefix(
    const std::vector<std::string> &targets,
    BObolSceneController *scene)
{
    if (!scene || targets.empty())
	return 0;

    std::vector<std::string> group_paths;
    const int summary_count = scene->getSceneTreeSummaryCount();
    for (int i = summary_count - 1; i >= 0; i--) {
	BObolSceneTreeSummary tree_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !tree_summary.valid ||
	    tree_summary.nodeKind != BObolSceneTreeSummary::NODE_GROUP ||
	    tree_summary.path.getLength() == 0 ||
	    BU_STR_EQUAL(tree_summary.path.getString(), "/"))
	    continue;

	const char *group_path = tree_summary.path.getString();
	int matches = 0;
	for (const std::string &target : targets) {
	    if (ged_obol_path_has_prefix(group_path, target.c_str())) {
		matches = 1;
		break;
	    }
	}
	if (!matches)
	    continue;

	SoGroup *group = scene->findGroup(group_path);
	if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    continue;
	SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
	if (scene_group->overlayIntent.getValue())
	    continue;

	ged_obol_append_unique_path(group_paths, group_path);
    }

    int changed = 0;
    for (const std::string &group_path : group_paths) {
	if (scene->removeGroup(group_path.c_str()) > 0)
	    changed = 1;
    }
    return changed;
}

extern "C" int
ged_draw_obol_groups_remove_for_path_prefix(
    struct ged *gedp,
    const char *path_prefix)
{
    if (!gedp || !path_prefix || !path_prefix[0])
	return 0;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> targets;
    ged_obol_append_unique_path(targets, path_prefix);
    const int changed = ged_obol_remove_groups_by_path_prefix(targets,
			scene);
    if (changed)
	scene->realizePending();
    return changed;
}

static int
ged_obol_apply_erase_prefix_transaction(
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    BObolSceneController *scene)
{
    if (!txn || !scene)
	return 0;

    std::vector<std::string> targets =
	ged_obol_transaction_paths(txn, result);
    if (targets.empty())
	return (result && result->status >= 0) ? 1 : 0;

    std::vector<std::string> source_instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, txn->view))
	    continue;
	const char *source_path = summary.path.getString();
	for (const std::string &target : targets) {
	    if (ged_obol_path_has_prefix(source_path, target.c_str())) {
		ged_obol_append_database_source_instance_key(
		    source_instance_keys, summary);
		break;
	    }
	}
    }

    int handled = !source_instance_keys.empty() ? 1 : 0;
    if (!source_instance_keys.empty())
	(void)ged_obol_remove_instance_keys(source_instance_keys, scene);
    if (!ged_obol_view_scope_is_independent(txn->view) &&
	ged_obol_remove_groups_by_path_prefix(targets, scene))
	handled = 1;
    if (handled)
	scene->realizePending();
    if (!handled && result && result->status >= 0)
	return 1;
    return handled;
}

static int
ged_obol_apply_redraw_transaction(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    uint32_t UNUSED(source_revision),
    BObolSceneController *scene)
{
    if (!gedp || !txn || !scene || !gedp->dbip)
	return 0;
    (void)view_ctx;
    const std::vector<std::string> targets =
	ged_obol_transaction_paths(txn, result);
    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (!source)
	    continue;
	if (targets.empty()) {
	    (void)source->setCompactInstanceDisplayStateForPath("", TRUE,
		1, TRUE, 0, FALSE, 0, FALSE);
	    continue;
	}
	for (const std::string &target : targets)
	    (void)source->setCompactInstanceDisplayStateForPath(target.c_str(),
		TRUE, 1, TRUE, 0, FALSE, 0, FALSE);
    }
    BObolViewController *controller = ged_draw_obol_controller(gedp);
    if (controller)
	controller->requestRender("ged-redraw-transaction");
    return 1;
}

static int
ged_obol_apply_cached_path_metadata_to_scene(
    struct ged *gedp,
    BObolSceneController *scene,
    struct ged_view_context *view_ctx,
    const char *path,
    int draw_mode)
{
    if (!gedp || !gedp->dbip || !scene || !path || !path[0])
	return 0;

    struct BObolDrawMetadataRecord metadata;
    bobol_draw_metadata_record_init(&metadata);
    const char *metadata_path = ged_obol_skip_leading_slash(path);
    if (!metadata_path || !metadata_path[0])
	return 0;

    if (bobol_draw_path_metadata_cache_get(gedp->dbip, metadata_path,
	    &metadata) != BRLCAD_OK &&
	bobol_draw_path_metadata_cache_refresh(gedp->dbip,
		metadata_path, NULL) == BRLCAD_OK) {
	(void)bobol_draw_path_metadata_cache_get(gedp->dbip,
		metadata_path, &metadata);
    }
    if (!metadata.directoryFound)
	return 0;

    const std::string mode_key =
	ged_obol_database_source_instance_key_for_mode(view_ctx, path,
	    draw_mode);
    SoBRLDatabaseSource *source =
	scene->findDatabaseSourceInstance(mode_key.c_str());
    if (!source) {
	const std::string base_key =
	    ged_obol_database_source_instance_key(view_ctx, path);
	source = scene->findDatabaseSourceInstance(base_key.c_str());
    }

    return ged_obol_apply_draw_metadata_to_source(source, &metadata);
}

static int
ged_obol_replace_deferred_paths(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    struct db_i *dbip,
    const std::vector<std::string> &paths,
    int draw_mode,
    uint32_t source_revision,
    BObolSceneController *scene,
    int preserve_existing_revision = 0,
    const struct ged_draw_appearance_settings *appearance_settings = NULL)
{
    if (!gedp || !dbip || paths.empty() || !scene)
	return 0;

    ged_obol_scene_mutation_batch_scope batch(scene, paths.size(),
	    paths.size());
    int changed = 0;
    for (const std::string &path : paths) {
	if (ged_obol_replace_path(gedp, view_ctx, dbip, path.c_str(),
				  draw_mode, source_revision, scene, 0,
				  preserve_existing_revision,
				  appearance_settings) > 0)
	    changed = 1;
	const std::string instance_key =
	    ged_obol_database_source_instance_key_for_mode(view_ctx,
		path.c_str(), draw_mode);
	if (!instance_key.empty() &&
	    ged_obol_publish_deferred_root_proxy(gedp, view_ctx,
		path.c_str(), instance_key.c_str(), draw_mode))
	    changed = 1;
	if (ged_obol_apply_cached_path_metadata_to_scene(gedp, scene,
		view_ctx, path.c_str(), draw_mode))
	    changed = 1;
    }

    return changed;
}

static int
ged_obol_summary_matches_transaction_targets(
    const BObolDatabaseSourceSummary &summary,
    struct ged_view_context *view_ctx,
    const std::vector<std::string> &targets,
    int draw_mode)
{
    if (!summary.valid ||
	!ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	!ged_obol_database_source_summary_matches_mode(summary, draw_mode))
	return 0;

    const char *source_path = summary.path.getString();
    if (!source_path || !source_path[0])
	return 0;

    for (const std::string &target : targets) {
	if (ged_obol_path_has_prefix(source_path, target.c_str()) ||
		ged_obol_path_has_semantic_prefix(source_path,
		    target.c_str()) ||
		ged_obol_path_has_component_name(source_path,
		    target.c_str(), 0) ||
		ged_obol_path_has_semantic_component_name(source_path,
		    target.c_str(), 0))
	    return 1;
    }
    return 0;
}

static int
ged_obol_publish_database_source_summary_to_scene(
    struct ged *gedp,
    BObolSceneController *scene,
    const BObolDatabaseSourceSummary &summary)
{
    if (!gedp || !gedp->dbip || !scene || !summary.valid)
	return 0;

    const char *instance_key = summary.instanceKey.getString();
    const char *source_path = summary.path.getString();
    if (!instance_key || !instance_key[0] || !source_path || !source_path[0])
	return 0;

    const char *representation_key = summary.representationKey.getString();
    if (!representation_key || !representation_key[0])
	representation_key = instance_key;

    const char *target_group_path = summary.parentGroupPath.getString();
    BObolDatabaseSourcePublishState state;
    state.sourceInstanceKey = instance_key;
    state.sourcePath = source_path;
    state.sourceRepresentationKey = representation_key;
    state.targetGroupPath =
	(target_group_path && target_group_path[0]) ? target_group_path : NULL;
    state.database = gedp->dbip;
    state.drawMode = summary.drawMode;
    state.representationMode = summary.representationMode;
    state.sourceRevisionValid = TRUE;
    state.sourceRevision = summary.sourceRevision;
    state.inputsRevision = summary.inputsRevision;
    state.visible = summary.visible;
    state.selected = summary.selected;
    state.highlighted = summary.highlighted;
    state.lineStyle = summary.lineStyle;
    state.lineWidth = summary.lineWidth;
    state.transparency = summary.transparency;
    state.colorOverride = summary.colorOverride;
    state.color = summary.color;
    state.materialColorValid = summary.materialColorValid;
    state.materialColor = summary.materialColor;
    state.materialRevision = summary.materialRevision;
    state.roleFlagsValid = TRUE;
    state.roleFlags = summary.realizationRoleFlags;
    state.viewPolicyValid = TRUE;
    state.viewDependent = summary.realizationViewDependent;
    state.csgLodEnabled = summary.realizationCsgLodEnabled;
    state.meshLodEnabled = summary.realizationMeshLodEnabled;
    state.viewScale = summary.realizationViewScale;
    state.lodScale = summary.realizationLodScale;
    state.viewWidth = summary.realizationViewWidth;
    state.viewHeight = summary.realizationViewHeight;
    state.botThreshold = summary.realizationBotThreshold;
    state.curveScale = summary.realizationCurveScale;
    state.pointScale = summary.realizationPointScale;
    state.placementValid = summary.drawMatrixValid ||
			   summary.drawCenterValid ||
			   summary.drawSizeValid;
    state.drawMatrixValid = summary.drawMatrixValid;
    state.drawMatrix = summary.drawMatrix;
    state.drawCenterValid = summary.drawCenterValid;
    state.drawCenter = summary.drawCenter;
    state.drawSizeValid = summary.drawSizeValid;
    state.drawSize = summary.drawSize;

    int changed = scene->publishDatabaseSourceInstance(state);
    if (changed < 0)
	return 0;

    int material_changed = scene->setDatabaseSourceInstanceMaterialPolicy(
			       instance_key, summary.materialPolicy);
    if (material_changed < 0)
	return 0;
    return (changed > 0 || material_changed > 0) ? 1 : 0;
}

static int
ged_obol_copy_matching_primary_sources_to_scene(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const std::vector<std::string> &targets,
    int draw_mode,
    int remove_draw_mode,
    BObolSceneController *scene)
{
    if (!gedp || targets.empty() || !scene)
	return 0;

    BObolSceneController *primary_scene =
	ged_draw_obol_scene_controller(gedp);
    if (!primary_scene || primary_scene == scene)
	return 0;

    std::vector<BObolDatabaseSourceSummary> summaries;
    const int source_count = primary_scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!primary_scene->getDatabaseSourceSummary(i, summary) ||
	    !ged_obol_summary_matches_transaction_targets(summary, view_ctx,
		targets, draw_mode))
	    continue;
	const char *instance_key = summary.instanceKey.getString();
	if (!instance_key || !instance_key[0])
	    continue;
	summaries.push_back(summary);
    }
    if (summaries.empty())
	return 0;

    std::set<std::string> retained_keys;
    for (const BObolDatabaseSourceSummary &summary : summaries)
	retained_keys.insert(summary.instanceKey.getString());
    std::vector<std::string> obsolete_keys;
    const int destination_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < destination_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid ||
	    !ged_obol_summary_matches_transaction_targets(summary, view_ctx,
		targets, remove_draw_mode))
	    continue;
	const std::string key = summary.instanceKey.getString();
	if (!key.empty() && retained_keys.find(key) == retained_keys.end())
	    obsolete_keys.push_back(key);
    }

    ged_obol_scene_mutation_batch_scope batch(scene, summaries.size(),
	    summaries.size() + obsolete_keys.size());
    int changed = 0;
    for (const BObolDatabaseSourceSummary &summary : summaries) {
	if (ged_obol_publish_database_source_summary_to_scene(gedp, scene,
		summary))
	    changed = 1;
    }
    for (const std::string &key : obsolete_keys) {
	if (scene->removeDatabaseSourceInstance(key.c_str()) > 0)
	    changed = 1;
    }

    if (changed)
	scene->realizePending();
    return 1;
}

int
ged_draw_obol_scene_sync_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    BObolSceneController *controller)
{
    if (!gedp || !txn || (result && result->status < 0))
	return 0;

    /* A retained-frontier operation has already updated compact instance
     * visibility on the owning source.  Replaying it as an ordinary
     * draw/erase would remove and republish that source, discard its active
     * view cuts, and restart deferred realization. */
    if (result && result->presentation_only)
	return 1;

    BObolSceneController *scene = controller ?
				  controller : ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    struct ged_view_context *view_ctx = txn->view;
    uint32_t source_revision = ged_obol_transaction_source_revision(result);
    int changed = 0;

    switch (txn->kind) {
	case GED_DRAW_TXN_DRAW: {
	    std::vector<std::string> paths =
		ged_obol_transaction_paths(txn, result);
	    if (result && result->affected_shapes <= 0 &&
		result->affected_groups <= 0 && paths.empty())
		break;
	    const int draw_mode =
		ged_obol_transaction_ged_draw_mode(gedp, txn);
	    const struct ged_draw_appearance_settings *appearance =
		    static_cast<const struct ged_draw_appearance_settings *>(
			txn->appearance);
	    const int remove_draw_mode =
		(appearance && !appearance->mixed_modes) ? -1 : draw_mode;
	    if (ged_obol_transaction_defer_leaf_expansion(txn)) {
		if (!paths.empty() &&
		    ged_obol_remove_paths(paths, view_ctx, scene,
			remove_draw_mode))
		    changed = 1;
		if (!paths.empty())
		    changed |= ged_obol_replace_deferred_paths(gedp, view_ctx,
			      gedp->dbip, paths, draw_mode, source_revision,
			      scene, 1, appearance);
		break;
	    }
	    if (scene != ged_draw_obol_scene_controller(gedp) &&
		    !paths.empty()) {
		changed = ged_obol_copy_matching_primary_sources_to_scene(
			      gedp, view_ctx, paths, draw_mode,
			      remove_draw_mode, scene);
		if (changed)
		    break;
	    }
	    if (paths.empty()) {
		changed = ged_obol_scene_sync_full_scene_impl(gedp,
			  view_ctx, source_revision, scene, 1);
	    } else {
		std::vector<std::string> source_paths =
		    ged_obol_primary_matching_database_source_paths(gedp,
			view_ctx, paths, draw_mode);
		if (source_paths.empty())
		    source_paths =
			ged_obol_drawn_source_paths(gedp, view_ctx, -1, &paths);
		if (!source_paths.empty()) {
		    std::vector<std::string> shadowed_targets =
			ged_obol_shadowed_target_source_paths(paths,
			    source_paths);
		    if (ged_obol_remove_paths(shadowed_targets, view_ctx,
					      scene))
			changed = 1;
		    if (ged_obol_remove_paths(paths, view_ctx, scene,
			    remove_draw_mode))
			changed = 1;
		    if (ged_obol_replace_paths(gedp->dbip, source_paths,
					       draw_mode, source_revision,
					       gedp, view_ctx, scene, 1, 1,
					       appearance))
			changed = 1;
		} else {
		    if (ged_obol_remove_paths(paths, view_ctx, scene,
			    remove_draw_mode))
			changed = 1;
		    changed |= ged_obol_replace_paths(gedp->dbip, paths,
						     draw_mode, source_revision,
						     gedp, view_ctx, scene, 1, 1,
						     appearance);
		}
	    }
	    break;
	}
	case GED_DRAW_TXN_ERASE: {
	    std::vector<std::string> paths =
		ged_obol_transaction_paths(txn, result);
	    if (paths.empty()) {
		changed = ged_draw_obol_scene_sync_full_scene(gedp,
			  view_ctx, source_revision, scene);
	} else {
	    std::vector<std::string> matching_instance_keys;
	    int compact_changed = 0;
	    const int source_count = scene->getDatabaseSourceCount();
	    for (int i = 0; i < source_count; i++) {
		SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
		/* A retained source frontier is the authority for a nested
		 * structural erase.  Do not also write the transient mask into
		 * authored occurrence visibility: that made a later draw/collapse
		 * unable to reveal the occurrence again. */
		if (source &&
		    !source->hasCompactInstanceVisibilityFrontier()) {
		    for (const std::string &target : paths) {
			if (source->setCompactInstanceDisplayStateForPath(
				target.c_str(), TRUE, 1, FALSE,
				0, FALSE, 0, FALSE) > 0)
			    compact_changed = 1;
		    }
		}
		BObolDatabaseSourceSummary summary;
		    if (!scene->getDatabaseSourceSummary(i, summary) ||
			!summary.valid ||
			!ged_obol_database_source_instance_in_scope(summary,
				view_ctx) ||
			!ged_obol_database_source_summary_matches_mode(
			    summary, txn->mode))
			continue;

		    const char *source_path = summary.path.getString();
		    const std::string owner_group_path =
			ged_obol_database_source_owner_group_path_from_summary(
			    summary);
		    for (const std::string &target : paths) {
			if (ged_obol_path_equal(source_path, target.c_str()) ||
			    ged_obol_path_equal(owner_group_path.c_str(),
						target.c_str())) {
			    ged_obol_append_database_source_instance_key(
				matching_instance_keys, summary);
			    break;
			}
		    }
		}
	    changed = matching_instance_keys.empty() ?
		      ged_obol_remove_paths(paths, view_ctx, scene, txn->mode) :
		      ged_obol_remove_instance_keys(matching_instance_keys,
						scene);
	    if (!matching_instance_keys.empty())
		scene->clearRealizationRepository();
	    if (compact_changed)
		changed = 1;
	}
	    break;
	}
	case GED_DRAW_TXN_CLEAR:
	case GED_DRAW_TXN_TEARDOWN:
	    changed = ged_obol_scene_clear_controller(scene);
	    scene->clearRealizationRepository();
	    break;
	case GED_DRAW_TXN_CLEAR_SCOPE:
	    changed = ged_obol_clear_database_sources_in_scope(scene,
		      view_ctx);
	    if (changed)
		scene->clearRealizationRepository();
	    if (changed)
		scene->realizePending();
	    break;
	case GED_DRAW_TXN_VISIBILITY:
	    changed = ged_obol_apply_visibility_transaction(txn, result,
		      scene);
	    break;
	case GED_DRAW_TXN_HIGHLIGHT:
	    changed = ged_obol_apply_highlight_transaction(txn, result,
		      scene);
	    break;
	case GED_DRAW_TXN_MATERIAL_CHANGED:
	    /* Core libged has already advanced the material revision.  The
	     * following REFRESH_MATERIAL_COLORS transaction applies colors. */
	    changed = (result && result->status >= 0) ? 1 : 0;
	    break;
	case GED_DRAW_TXN_REFRESH_MATERIAL_COLORS: {
	    const int refreshed =
		scene->refreshDatabaseSourceMaterialColorsFromDatabase(
		    ged_obol_fold_revision(ged_draw_material_revision(gedp)),
		    gedp->dbip);
	    changed = refreshed >= 0 ? 1 : 0;
	    break;
	}
	case GED_DRAW_TXN_STALE_SOURCE:
	    changed = ged_obol_apply_stale_source_transaction(gedp,
		      view_ctx, txn, result, source_revision, scene);
	    break;
	case GED_DRAW_TXN_ERASE_PREFIX:
	    changed = ged_obol_apply_erase_prefix_transaction(txn, result,
		      scene);
	    break;
	case GED_DRAW_TXN_REDRAW:
	    changed = ged_obol_apply_redraw_transaction(gedp, view_ctx,
		      txn, result, source_revision, scene);
	    break;
	case GED_DRAW_TXN_SOURCE_UPDATED:
	    changed = ged_obol_apply_source_update_transaction(gedp,
		      view_ctx, txn, result, source_revision, scene);
	    if (changed > 0)
		break;
	    changed = ged_draw_obol_scene_sync_full_scene(gedp, view_ctx,
		      source_revision, scene);
	    break;
	case GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
	    changed = ged_obol_apply_source_references_removed_transaction(
			  txn, result, scene);
	    if (changed > 0)
		break;
	    changed = ged_draw_obol_scene_sync_full_scene(gedp, view_ctx,
		      source_revision, scene);
	    break;
	case GED_DRAW_TXN_SOURCE_RENAMED:
	    if (txn->path && txn->new_path) {
		scene->renameRealizationObject(txn->path, txn->new_path);
		/* Retarget aggregate occurrence semantics explicitly before
		 * rebuilding.  This preserves exact-path handle/state matching
		 * without the ambiguous legacy fallback of matching any leaf with
		 * the same source name or list position. */
		for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
		    SoBRLDatabaseSource *source =
			scene->getDatabaseSource(i);
		    if (source)
			(void)source->retargetCompactOccurrencePaths(
			    txn->path, txn->new_path);
		}
		const int renamed = ged_draw_obol_database_source_rename_for_path(gedp,
			  txn->path, txn->new_path, source_revision);
		std::vector<std::string> renamed_targets;
		renamed_targets.push_back(txn->path);
		const int aggregates =
		    ged_obol_mark_matching_database_sources_stale(view_ctx,
			renamed_targets, 0, 1,
			SoBRLDatabaseSource::STALE_SOURCE, scene, txn->mode);
		if (aggregates > 0)
		    scene->realizePending();
		changed = renamed > 0 || aggregates > 0 ? 1 : 0;
		if (changed > 0)
		    break;
		if (scene->findDatabaseSource(txn->new_path)) {
		    changed = 1;
		    break;
		}
	    }
	    changed = ged_draw_obol_scene_sync_full_scene(gedp, view_ctx,
		      source_revision, scene);
	    break;
	case GED_DRAW_TXN_DEFAULT_DRAW_MODE:
	case GED_DRAW_TXN_TRANSPARENCY:
	case GED_DRAW_TXN_NONE:
	default:
	    break;
    }

    return changed ? 1 : 0;
}

static int
ged_obol_transaction_invalidates_view_lod(
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    int full_sync)
{
    if (full_sync)
	return 1;
    if (!txn)
	return 0;
    if (result && result->presentation_only)
	return 0;

    switch (txn->kind) {
	case GED_DRAW_TXN_DRAW:
	case GED_DRAW_TXN_ERASE:
	case GED_DRAW_TXN_CLEAR:
	case GED_DRAW_TXN_TEARDOWN:
	case GED_DRAW_TXN_CLEAR_SCOPE:
	case GED_DRAW_TXN_STALE_SOURCE:
	case GED_DRAW_TXN_ERASE_PREFIX:
	case GED_DRAW_TXN_SOURCE_UPDATED:
	case GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
	case GED_DRAW_TXN_SOURCE_RENAMED:
	    return 1;
	default:
	    return 0;
    }
}

extern "C" int
ged_draw_obol_scene_sync_attached_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result)
{
    if (!gedp || !txn)
	return 0;

    const int completed_direct_draw =
	ged_obol_transaction_is_completed_database_source_draw(gedp, txn,
	    result);
    BObolSceneController *primary_scene = ged_draw_obol_scene_controller(gedp);
    int changed = (!completed_direct_draw && primary_scene) ?
	ged_draw_obol_scene_sync_transaction(gedp, txn, result,
	    primary_scene) : 0;

    struct ged_obol_endpoint_sync_context {
	struct ged *gedp;
	const struct ged_draw_transaction *txn;
	const struct ged_draw_transaction_result *result;
	BObolSceneController *shared_scene;
	int changed;
    } sync_ctx = {gedp, txn, result, primary_scene, changed};

    const auto sync_endpoint = [](struct ged_view_context *view_ctx,
	BObolViewController *controller, void *userdata) -> int {
	struct ged_obol_endpoint_sync_context *ctx =
	    static_cast<struct ged_obol_endpoint_sync_context *>(userdata);
	if (!ctx || !view_ctx || !controller)
	    return 1;

	const int independent = ged_obol_view_scope_is_independent(view_ctx);
	SoNode *shared_root = ctx->shared_scene ?
	    ctx->shared_scene->getSceneRoot() : NULL;
	const int shared_bound = ged_obol_node_contains(
	    controller->getRenderSceneRoot(), shared_root);
	if (ctx->shared_scene && shared_bound == independent)
	    (void)ged_obol_bind_view_render_root(view_ctx, ctx->shared_scene,
		controller);

	if (!independent) {
	    if (ged_obol_transaction_invalidates_view_lod(ctx->txn,
		    ctx->result, 0))
		controller->clearViewLodState();
	    return 1;
	}
	if (ctx->txn->kind == GED_DRAW_TXN_CLEAR ||
	    ctx->txn->kind == GED_DRAW_TXN_TEARDOWN)
	    return 1;
	if (ged_obol_view_scope_is_independent(ctx->txn->view) &&
	    ctx->txn->view != view_ctx)
	    return 1;

	int full_sync = 0;
	if (ctx->txn->view != view_ctx) {
	    switch (ctx->txn->kind) {
		case GED_DRAW_TXN_SOURCE_UPDATED:
		case GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
		case GED_DRAW_TXN_SOURCE_RENAMED:
		case GED_DRAW_TXN_STALE_SOURCE:
		case GED_DRAW_TXN_MATERIAL_CHANGED:
		case GED_DRAW_TXN_REFRESH_MATERIAL_COLORS:
		case GED_DRAW_TXN_HIGHLIGHT:
		    full_sync = 1;
		    break;
		case GED_DRAW_TXN_REDRAW:
		    if (ctx->txn->view)
			return 1;
		    full_sync = 1;
		    break;
		default:
		    return 1;
	    }
	}

	struct ged_draw_transaction local_txn = *ctx->txn;
	local_txn.view = view_ctx;
	BObolSceneController *scene = controller->getSceneController();
	const int endpoint_changed = full_sync ?
	    ged_draw_obol_scene_sync_full_scene(ctx->gedp, view_ctx,
		ged_obol_transaction_source_revision(ctx->result), scene) :
	    ged_draw_obol_scene_sync_transaction(ctx->gedp, &local_txn,
		ctx->result, scene);
	if (ged_obol_transaction_invalidates_view_lod(&local_txn,
		ctx->result, full_sync))
	    controller->clearViewLodState();
	if (endpoint_changed)
	    ctx->changed = 1;
	return 1;
    };

    ged_bobol_view_controllers_foreach(gedp, sync_endpoint, &sync_ctx);
    return sync_ctx.changed ? 1 : 0;
}

static void
ged_obol_progressive_autoview_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result)
{
    if (!gedp || !txn)
	return;

    const int successful = !result || result->status >= 0;
    const int deferred = successful &&
	!(result && result->presentation_only) &&
	txn->kind == GED_DRAW_TXN_DRAW &&
	ged_obol_transaction_defer_leaf_expansion(txn);
    const int arm_autoview = deferred && txn->autoview;
    const int invalidate =
	ged_obol_transaction_invalidates_view_lod(txn, result, 0);

    struct ged_obol_progressive_transaction_context {
	struct ged *gedp;
	const struct ged_draw_transaction *txn;
	const struct ged_draw_transaction_result *result;
	int deferred;
	int arm_autoview;
	int invalidate;
    } progressive_ctx = {gedp, txn, result, deferred, arm_autoview,
	invalidate};

    const auto update_endpoint = [](struct ged_view_context *view_ctx,
	BObolViewController *controller, void *userdata) -> int {
	struct ged_obol_progressive_transaction_context *ctx =
	    static_cast<struct ged_obol_progressive_transaction_context *>(userdata);
	if (!ctx || !view_ctx || !controller ||
	    (ctx->txn->view && ctx->txn->view != view_ctx))
	    return 1;
	ged_obol_progressive_provider_data *data =
	    static_cast<ged_obol_progressive_provider_data *>(
		controller->findProgressiveProviderData(
		    ged_obol_progressive_advance_provider));
	if (!data)
	    return 1;

	if (ctx->invalidate)
	    controller->clearViewLodState();
	if (!ctx->deferred) {
	    if (ctx->invalidate) {
		ged_obol_retire_all_deferred_jobs(data);
		ged_obol_cleanup_retired_jobs(data);
		data->pending_autoview = 0;
		data->deferred_refine_stage = 0;
		data->deferred_paths.clear();
	    }
	    return 1;
	}

	data->deferred_paths = ged_obol_transaction_paths(ctx->txn, ctx->result);
	if (ctx->txn->appearance) {
	    data->deferred_appearance =
		*static_cast<const struct ged_draw_appearance_settings *>(
		    ctx->txn->appearance);
	} else {
	    data->deferred_appearance = ged_draw_appearance_settings
		GED_DRAW_APPEARANCE_SETTINGS_INIT;
	    data->deferred_appearance.draw_mode =
		ged_obol_transaction_ged_draw_mode(ctx->gedp, ctx->txn);
	}
	data->deferred_appearance.defer_leaf_expansion = 0;
	data->deferred_refine_stage =
	    data->deferred_paths.empty() ? 0 : 1;

	if (data->deferred_refine_stage == 1 &&
	    !ged_obol_start_deferred_realization(data, controller,
		data->deferred_appearance.draw_mode)) {
	    data->deferred_refine_stage = 3;
	    data->deferred_paths.clear();
	}
	if (ctx->arm_autoview)
	    (void)ged_obol_progressive_autoview_arm(data,
		BV_AUTOVIEW_SCALE_DEFAULT);
	controller->markProgressiveWorkPending();
	return 1;
    };
    ged_bobol_view_controllers_foreach(gedp, update_endpoint,
	&progressive_ctx);
}

static void
ged_obol_transaction_observer(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    void *client_data)
{
    (void)client_data;
    (void)ged_draw_obol_scene_sync_attached_transaction(gedp, txn, result);
    ged_obol_progressive_autoview_transaction(gedp, txn, result);
}

static int
ged_obol_observer_ensure(struct ged *gedp, struct ged_drawable *gdp)
{
    if (!gedp || !gdp)
	return 0;
    ged_obol_state *state = ged_obol_state_get(gdp, 1);
    if (!state)
	return 0;
    if (state->observer_token)
	return 1;

    state->observer_token = ged_draw_observer_add(gedp,
	ged_obol_transaction_observer, NULL);
    return state->observer_token ? 1 : 0;
}

static void
ged_obol_configure_compact_realization(struct ged *gedp,
	struct ged_view_context *view_ctx, BObolSceneController *scene_controller)
{
    if (!scene_controller)
	return;

    (void)gedp;
    (void)view_ctx;
}

static int
ged_obol_register_progressive_provider(struct ged *gedp,
			       struct ged_view_context *view_ctx,
			       BObolViewController *controller)
{
    if (!gedp || !view_ctx || !controller)
	return 0;
    ged_obol_progressive_provider_data *existing =
	static_cast<ged_obol_progressive_provider_data *>(
	    controller->findProgressiveProviderData(
		ged_obol_progressive_advance_provider));
    if (existing)
	return (existing->gedp == gedp && existing->view_ctx == view_ctx) ?
	    1 : 0;

    ged_obol_progressive_provider_data *data =
	new (std::nothrow) ged_obol_progressive_provider_data;
    if (!data)
	return 0;

    data->gedp = gedp;
    data->view_ctx = view_ctx;
    uint64_t token = controller->registerProgressiveProvider(
	ged_obol_progressive_advance_provider, data,
	[](void *userdata) {
	    delete static_cast<ged_obol_progressive_provider_data *>(userdata);
	});
    if (!token) {
	delete data;
	return 0;
    }
    return 1;
}

extern "C" int
ged_draw_obol_progressive_available(struct ged *gedp, struct ged_view_context *view_ctx)
{
    (void)gedp;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    return (controller && controller->findProgressiveProviderData(
	ged_obol_progressive_advance_provider)) ? 1 : 0;
}

/* Progressive (coarse-first / deferred-leaf) drawing is driven by the part of
 * the view policy used by the requested representation.  Shaded and
 * hidden-line modes consume mesh LoD.  Wire may consume CSG LoD for primitive
 * plotting and mesh LoD for PoP-backed BoTs, so either policy can make its
 * production path progressive. */
extern "C" int
ged_draw_obol_view_lod_enabled(struct ged *gedp,
			      struct ged_view_context *view_ctx,
			      int draw_mode)
{
    (void)gedp;
    if (!ged_bobol_view_controller(view_ctx))
	return 0;

    ged_draw_source_lod_policy policy;
    if (!ged_draw_source_lod_policy_get(&policy, view_ctx) ||
	policy.policy == BV_LOD_OFF)
	return 0;

    switch (draw_mode) {
	case GED_DRAW_MODE_SHADED_BOTS:
	case GED_DRAW_MODE_SHADED:
	case GED_DRAW_MODE_HIDDEN_LINE:
	    return policy.mesh_enabled ? 1 : 0;
	case GED_DRAW_MODE_WIRE:
	    return (policy.csg_enabled || policy.mesh_enabled) ? 1 : 0;
	case GED_DRAW_MODE_EVAL_WIRE:
	case GED_DRAW_MODE_EVAL_POINTS:
	default:
	    return 0;
    }
}

static int
ged_obol_bind_view_render_root(struct ged_view_context *view_ctx,
			       BObolSceneController *shared_scene,
			       BObolViewController *view_controller)
{
    if (!shared_scene || !view_controller)
	return 0;

    SoNode *shared_root = shared_scene->getSceneRoot();
    if (!shared_root)
	return 0;

    view_controller->getSceneController()->shareRealizationRepository(
	shared_scene);

    SoNode *local_root = view_controller->getSceneRoot();
    if (!local_root || local_root == shared_root) {
	SoBRLSceneGroup *view_group = new SoBRLSceneGroup;
	std::string group_path("_view/");
	group_path += ged_obol_view_scope_name(view_ctx);
	view_group->groupPath = group_path.c_str();
	view_controller->setSceneRoot(view_group);
	local_root = view_group;
    }

    SoSeparator *render_root = new SoSeparator;
    if (!ged_obol_view_scope_is_independent(view_ctx))
	render_root->addChild(shared_root);
    /* The interlay belongs after model geometry but before the view-local
     * faceplate and interactive feature root. */
    SoGroup *interlay = view_controller->getFramebufferInterlayRoot();
    if (!interlay)
	return 0;
    render_root->addChild(interlay);
    if (local_root && local_root != shared_root)
	render_root->addChild(local_root);
    view_controller->setRenderSceneRoot(render_root);
    return 1;
}

static int
ged_obol_node_contains(SoNode *root, SoNode *target)
{
    if (!root || !target)
	return 0;
    if (root == target)
	return 1;
    if (!root->isOfType(SoGroup::getClassTypeId()))
	return 0;

    SoGroup *group = static_cast<SoGroup *>(root);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (ged_obol_node_contains(group->getChild(i), target))
	    return 1;
    }
    return 0;
}

static int
ged_draw_obol_attach_view_common(struct ged *gedp,
				 struct ged_view_context *view_ctx,
				 BObolSceneController *scene_controller,
				 BObolViewController *view_controller,
				 int sync_current_scene)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp || !scene_controller || !view_controller)
	return 0;

    ged_obol_configure_compact_realization(gedp, view_ctx,
	    scene_controller);

    if (!ged_obol_observer_ensure(gedp, gdp))
	return 0;

    if (!ged_draw_obol_scene_controller_ensure_owned(gedp,
	    sync_current_scene))
	return 0;
    BObolSceneController *shared_scene =
	ged_draw_obol_scene_controller(gedp);
    if (!shared_scene)
	return 0;

    (void)ged_view_context_host_attach(gedp, view_ctx);
    if (!ged_view_context_obol_attachment_bind(view_ctx,
	    view_controller->getViewAttachment()))
	return 0;

    /* The shared controller may itself be a single-view endpoint.  Its scene
     * root already is the shared composition, so replacing it with a local
     * render composition would make the shared scene self-referential. */
    if (view_controller != ged_draw_obol_controller(gedp)) {
	ged_obol_progressive_provider_data *existing =
	    static_cast<ged_obol_progressive_provider_data *>(
		view_controller->findProgressiveProviderData(
		    ged_obol_progressive_advance_provider));
	SoNode *shared_root = shared_scene->getSceneRoot();
	const int shared_visible =
	    !ged_obol_view_scope_is_independent(view_ctx);
	const int shared_bound = ged_obol_node_contains(
	    view_controller->getRenderSceneRoot(), shared_root);
	if (!existing || shared_visible != shared_bound) {
	    if (!ged_obol_bind_view_render_root(view_ctx, shared_scene,
		    view_controller)) {
		(void)ged_view_context_obol_attachment_unbind(view_ctx,
		    view_controller->getViewAttachment());
		return 0;
	    }
	}

	if (sync_current_scene &&
		ged_obol_view_scope_is_independent(view_ctx))
	    (void)ged_draw_obol_scene_sync_full_scene(gedp, view_ctx,
		ged_obol_fold_revision(ged_draw_scene_revision(gedp)),
		view_controller->getSceneController());
    }

    if (!ged_obol_register_progressive_provider(gedp, view_ctx,
	    view_controller)) {
	(void)ged_view_context_obol_attachment_unbind(view_ctx,
	    view_controller->getViewAttachment());
	return 0;
    }

    /* Every graphical endpoint uses the same automatic display-LoD service.
     * Starting it only when cold structural bounds happened to be missing made
     * cold cache misses work while warm caches (and small models with directly
     * available bounds) never traversed their valid source-mesh requests.
     * Failure is non-fatal: the standing compact boxes remain a useful
     * fallback and an explicit service start may be retried later. */
    (void)ged_obol_lod_service_ensure(view_controller);

    return 1;
}

extern "C" int
ged_draw_obol_controller_attach_opaque_for_view(struct ged *gedp,
	struct ged_view_context *view_ctx,
	void *controller,
	int sync_current_scene)
{
    BObolViewController *view_controller =
	static_cast<BObolViewController *>(controller);
    if (!view_controller)
	return 0;
    return ged_draw_obol_attach_view_common(gedp, view_ctx,
	view_controller->getSceneController(), view_controller,
	sync_current_scene);
}

extern "C" int
ged_draw_obol_render_endpoint_ensure_for_view(struct ged *gedp,
	struct ged_view_context *view_ctx,
	int sync_current_scene)
{
    (void)sync_current_scene;
    if (!gedp || !view_ctx)
	return 0;
    if (ged_view_context_display_endpoint_get(view_ctx))
	return 1;

    bobol_init(NULL);
    BObolViewController *controller =
	new (std::nothrow) BObolViewController(new SoBRLSceneGroup);
    if (!controller)
	return 0;

    bobol_display_endpoint_t *endpoint =
	bobol_display_endpoint_create(controller,
	    BOBOL_ENDPOINT_OWN_CONTROLLER);
    if (!endpoint) {
	delete controller;
	return 0;
    }
    if (!ged_view_context_display_endpoint_set(view_ctx, endpoint, 1)) {
	bobol_display_endpoint_destroy(endpoint);
	return 0;
    }
    return 1;
}

void
ged_draw_obol_scene_controller_detach(struct ged *gedp)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp)
	return;

    ged_obol_state *state = ged_obol_state_get(gdp, 0);
    if (state && state->observer_token) {
	(void)ged_draw_observer_remove(gedp, state->observer_token);
	state->observer_token = 0;
    }

    if (state) {
	delete state;
	gdp->gd_obol_state = NULL;
    }
}

extern "C" void
ged_draw_obol_controller_detach_for_view(struct ged *gedp,
					 struct ged_view_context *view_ctx)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp || !view_ctx)
	return;

    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (!controller)
	return;
    const uint64_t token = controller->findProgressiveProviderToken(
	ged_obol_progressive_advance_provider);
    if (token)
	controller->unregisterProgressiveProvider(token);
    controller->stopManagedLodService();
    (void)ged_view_context_obol_attachment_unbind(view_ctx,
	controller->getViewAttachment());
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
