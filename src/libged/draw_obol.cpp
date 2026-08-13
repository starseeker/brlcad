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

#include "ged/display_obol_private.h"

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
#include "BObol/BSourceRealization.h"
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
#include "bu/datetime.h"
#include "bu/units.h"
#include "bu/vls.h"
#include "ged.h"
#include "ged/db_index.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "ged/selection.h"
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
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static const char ged_obol_group_intent_prefix[] = "ged-draw-group:";
static const char ged_obol_database_source_mode_key_marker[] =
    ":ged-draw-mode:";
static std::atomic<bool> ged_obol_test_force_face_set_failure(false);

/* A serialized database object is only the floor for cold preparation
 * memory: import owns decoded arrays, and PoP/PCA construction may also own
 * topology and publication buffers.  Give every database-reading task an
 * explicit conservative reservation so an unknown (zero) source count never
 * bypasses the service's byte governor.  Version-4 directory lengths count
 * 128-byte records rather than bytes. */
static size_t
ged_obol_cold_object_working_set_bytes(const struct db_i *dbip,
	const struct directory *dp, size_t expansion, size_t fixedBytes)
{
    if (!dp)
	return fixedBytes ? fixedBytes : 1;

    size_t encoded = dp->d_len;
    if (dbip && db_version(dbip) < 5) {
	if (encoded > SIZE_MAX / sizeof(union record))
	    encoded = SIZE_MAX;
	else
	    encoded *= sizeof(union record);
    }
    if (encoded == SIZE_MAX || (expansion && encoded > SIZE_MAX / expansion))
	return SIZE_MAX;
    size_t estimate = expansion ? encoded * expansion : encoded;
    if (estimate > SIZE_MAX - fixedBytes)
	return SIZE_MAX;
    estimate += fixedBytes;
    return estimate ? estimate : 1;
}

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
static int ged_obol_stream_cached_leaf_manifest(
    SoBRLDatabaseSource *source,
    struct db_i *database,
    int draw_mode,
    uint32_t source_revision,
    BObolCompactOccurrenceStream *stream,
    void *user_data);
static int ged_obol_store_leaf_proxy_manifest(
    struct db_i *database,
    const SoBRLDatabaseSource *source,
    BObolCompactOccurrenceStream *stream,
    void *user_data);
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
	pending_autoview_bounds_complete(0),
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
    int pending_autoview_bounds_complete;
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
	source(NULL), snapshotSourceDatabase(NULL),
	ownsSource(FALSE),
	sourceDrawMode(SoBRLDatabaseSource::WIREFRAME),
	primaryScene(FALSE), allowWireFallback(FALSE),
	streamBatchQuantum(256),
	streamMergeMicrosecondsPerItem(0.0)
    {
    }

    ~ged_obol_deferred_realization_item(void)
    {
	if (ownsSource && source)
	    source->unref();
	if (snapshotSourceDatabase)
	    db_close(snapshotSourceDatabase);
    }

    SoBRLDatabaseSource *source;
    /* Ref-counted live database handle used only to construct the immutable
     * closure snapshot on the worker. */
    struct db_i *snapshotSourceDatabase;
    SbBool ownsSource;
    std::string instanceKey;
    int drawMode = 0;
    /* drawMode is the GED presentation mode passed to manifest callbacks.
     * SoBRLDatabaseSource deliberately has only wire/shaded draw modes, with
     * shaded-bots versus shaded-all carried by representationMode.  Preserve
     * that distinct source field as the asynchronous result-admission stamp. */
    int sourceDrawMode;
    uint32_t sourceRevision = 0;
    uint32_t inputsRevision = 0;
    uint32_t viewRevision = 0;
    int representationMode = 0;
    SbBool primaryScene;
    SbBool allowWireFallback;
    /* GUI-thread adaptive publication quantum. */
    size_t streamBatchQuantum;
    /* Conservative observed merge cost.  Batch size alone is not a useful
     * predictor once retained occurrence maps grow from hundreds to tens of
     * thousands of records. */
    double streamMergeMicrosecondsPerItem;
    /* Worker->pump stream of completed per-leaf occurrences.  Shared with the
     * realization service so it outlives every push. */
    std::shared_ptr<BObolCompactOccurrenceStream> stream;
};

struct ged_obol_deferred_realization_job {
    enum State {
	PENDING = BOBOL_SOURCE_REALIZATION_PENDING,
	RUNNING = BOBOL_SOURCE_REALIZATION_RUNNING,
	COMPLETE = BOBOL_SOURCE_REALIZATION_COMPLETE,
	FAILED = BOBOL_SOURCE_REALIZATION_FAILED,
	CANCELLED = BOBOL_SOURCE_REALIZATION_CANCELLED
    };

    ged_obol_deferred_realization_job(void) :
	stagingLeaseReleaseAfter(0)
    {
    }

    ~ged_obol_deferred_realization_job(void)
    {
	cancel();
    }

    void cancel(void)
    {
	if (realization)
	    realization->cancel();
	for (const std::unique_ptr<ged_obol_deferred_realization_item> &item :
	     items) {
	    if (item && item->stream)
		item->stream->requestCancel();
	}
    }

    bool start(void)
    {
	if (items.empty() || realization)
	    return false;
	std::vector<BObolSourceRealizationRequest> requests(items.size());
	for (size_t i = 0; i < items.size(); i++) {
	    ged_obol_deferred_realization_item *item = items[i].get();
	    if (!item || !item->source || !item->snapshotSourceDatabase)
		return false;
	    BObolSourceRealizationRequest &request = requests[i];
	    request.source = item->source;
	    request.snapshotSourceDatabase = item->snapshotSourceDatabase;
	    request.stream = item->stream;
	    request.clientToken = i;
	    request.estimatedWorkingSetBytes = 256ULL * 1024ULL * 1024ULL;
	    request.drawMode = item->drawMode;
	    request.allowWireFallback = item->allowWireFallback;
	    request.probeWarmManifest =
		ged_obol_stream_cached_leaf_manifest;
	    request.storeManifest = ged_obol_store_leaf_proxy_manifest;
	}
	realization = BObolSourceRealizationCoordinator::global().submit(
	    requests);
	for (size_t i = 0; i < items.size(); i++) {
	    if (!requests[i].source) {
		items[i]->ownsSource = FALSE;
		items[i]->snapshotSourceDatabase = NULL;
	    }
	}
	return realization ? true : false;
    }

    int stateValue(void) const
    {
	return realization ? realization->state() : PENDING;
    }

    /* Completed streams retain their bounded cold-source leases briefly after
     * adoption so the following view-LoD traversal can capture them.  Once a
     * task owns a lease its provider shared_ptr, not this grace period, keeps
     * the arrays alive. */
    int64_t stagingLeaseReleaseAfter;
    std::vector<std::unique_ptr<ged_obol_deferred_realization_item>> items;
    std::shared_ptr<BObolSourceRealizationJob> realization;
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
	shared_controller(NULL),
	full_sync(0)
    {
    }

    ~ged_obol_state(void)
    {
	delete shared_controller;
    }

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

static std::vector<std::string>
ged_obol_transaction_paths(const struct ged_draw_transaction *txn,
			   const struct ged_draw_transaction_result *result)
{
    std::vector<std::string> paths;
    if (result)
	for (size_t i = 0; i < result->path_count; i++)
	    ged_obol_append_unique_path(paths, result->paths[i]);

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

/* Complete a physical-source match with semantic roots that do not yet have
 * a source.  The transaction reducer records those roots before this backend
 * adapter runs, so an existing matching source must not cause a sibling
 * request in the same transaction to be dropped.  Keep the shallowest owner:
 * drawing a broader path subsumes previously drawn descendants. */
static void
ged_obol_complete_source_roots(std::vector<std::string> &source_paths,
	const std::vector<std::string> &targets)
{
    for (const std::string &target : targets) {
	bool covered = false;
	for (const std::string &source_path : source_paths) {
	    if (ged_obol_path_has_prefix(target.c_str(),
		    source_path.c_str()) ||
		ged_obol_path_has_semantic_prefix(target.c_str(),
		    source_path.c_str())) {
		covered = true;
		break;
	    }
	}
	if (!covered)
	    ged_obol_append_unique_path(source_paths, target.c_str());
    }

    std::vector<std::string> roots;
    roots.reserve(source_paths.size());
    for (size_t i = 0; i < source_paths.size(); i++) {
	bool shadowed = false;
	for (size_t j = 0; j < source_paths.size(); j++) {
	    if (i == j)
		continue;
	    if (!ged_obol_path_equal(source_paths[i].c_str(),
		    source_paths[j].c_str()) &&
		(ged_obol_path_has_prefix(source_paths[i].c_str(),
		     source_paths[j].c_str()) ||
		 ged_obol_path_has_semantic_prefix(source_paths[i].c_str(),
		     source_paths[j].c_str()))) {
		shadowed = true;
		break;
	    }
	}
	if (!shadowed)
	    ged_obol_append_unique_path(roots, source_paths[i].c_str());
    }
    source_paths.swap(roots);
}

struct ged_obol_drawn_source_path_mode {
    std::string path;
    int mode;
    bool appearanceValid = false;
    struct ged_draw_appearance_settings appearance =
	GED_DRAW_APPEARANCE_SETTINGS_INIT;
};

static void
ged_obol_append_unique_path_mode(
    std::vector<ged_obol_drawn_source_path_mode> &path_modes,
    const char *path,
    int mode,
    const struct ged_draw_appearance_settings *appearance = NULL)
{
    if (!path || !path[0])
	return;

    for (ged_obol_drawn_source_path_mode &entry : path_modes) {
	if (entry.mode == mode &&
	    ged_obol_path_equal(entry.path.c_str(), path)) {
	    if (appearance) {
		entry.appearance = *appearance;
		entry.appearanceValid = true;
	    }
	    return;
	}
    }

    ged_obol_drawn_source_path_mode entry;
    entry.path = path;
    entry.mode = mode;
    if (appearance) {
	entry.appearance = *appearance;
	entry.appearanceValid = true;
    }
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
				     entry.mode,
				     entry.appearanceValid ?
				     &entry.appearance : NULL);
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
    int mode,
    const struct ged_draw_appearance_settings *appearance = NULL)
{
    if (!ctx || !path || !path[0])
	return;

    const std::string normalized = ged_obol_normalized_path_string(path);
    if (normalized.empty())
	return;

    std::string key = normalized;
    key += '\037';
    key += std::to_string(mode);
    if (ctx->pathModeKeys.find(key) != ctx->pathModeKeys.end()) {
	if (appearance) {
	    for (ged_obol_drawn_source_path_mode &entry : ctx->pathModes) {
		if (entry.mode == mode &&
		    ged_obol_path_equal(entry.path.c_str(), path)) {
		    entry.appearance = *appearance;
		    entry.appearanceValid = true;
		    break;
		}
	    }
	}
	return;
    }
    ctx->pathModeKeys.insert(key);

    ged_obol_drawn_source_path_mode entry;
    entry.path = path;
    entry.mode = mode;
    if (appearance) {
	entry.appearance = *appearance;
	entry.appearanceValid = true;
    }
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

static int
ged_obol_collect_frontier_root(
    const struct ged_draw_frontier_root_record *record,
    void *userdata)
{
    ged_obol_drawn_source_path_ctx *ctx =
	static_cast<ged_obol_drawn_source_path_ctx *>(userdata);
    if (!ctx || !record || !record->path || !record->path[0])
	return 1;

    bool matched = !ctx->targets || ctx->targets->empty();
    if (!matched) {
	for (const std::string &target : *ctx->targets) {
	    if (ged_obol_path_has_prefix(record->path, target.c_str()) ||
		ged_obol_path_has_prefix(target.c_str(), record->path) ||
		ged_obol_path_has_semantic_prefix(record->path,
		    target.c_str()) ||
		ged_obol_path_has_semantic_prefix(target.c_str(),
		    record->path)) {
		matched = true;
		break;
	    }
	}
    }
    if (!matched)
	return 1;

    ged_obol_drawn_source_append_path_mode(ctx, record->path, record->mode,
	record->appearance);
    ged_obol_drawn_source_append_path(ctx, record->path);
    return 1;
}

static void
ged_obol_collect_drawn_frontier_roots(ged_obol_drawn_source_path_ctx *ctx)
{
    if (!ctx || !ctx->gedp)
	return;
    (void)ged_draw_frontier_foreach_root(ctx->gedp, ctx->viewCtx,
	ctx->mode, ged_obol_collect_frontier_root, ctx);
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
    ged_obol_collect_drawn_frontier_roots(&ctx);
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
    ged_obol_collect_drawn_frontier_roots(&ctx);
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

    ged_view_lod_policy policy;
    if (!ged_view_lod_policy_get(&policy, policy_view))
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
    publish_state.materialPolicyValid = TRUE;
    publish_state.materialPolicy = SoBRLDatabaseSource::MATERIAL_DATABASE;
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
				      preserve_existing_revision,
				      entry.appearanceValid ?
				      &entry.appearance : NULL) > 0)
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

static int
ged_obol_set_group_transparency(BObolSceneController *scene,
				const char *group_path,
				float transparency)
{
    if (!scene || !group_path || !group_path[0])
	return 0;
    SoGroup *group = scene->findGroup(group_path);
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return 0;
    SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
    (void)scene->setGroupDisplayState(group_path,
	scene_group->visible.getValue(), scene_group->selected.getValue(),
	scene_group->highlighted.getValue(),
	scene_group->lineStyle.getValue(), scene_group->lineWidth.getValue(),
	transparency, scene_group->colorOverride.getValue(),
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

static int
ged_obol_scene_highlight_state_set(BObolSceneController *scene,
	int highlighted)
{
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (source && !highlighted)
	    (void)source->clearCompactInstanceHighlightOverrides();
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
ged_obol_transaction_path_matches(const char *candidate, const char *target,
				  enum ged_scene_path_match match)
{
    if (!candidate || !target)
	return 0;
    if (ged_obol_path_equal(candidate, target) ||
	ged_obol_semantic_path_string(candidate) ==
	ged_obol_semantic_path_string(target))
	return 1;
    if (match == GED_SCENE_PATH_MATCH_SUBTREE)
	return ged_obol_path_has_prefix(candidate, target) ||
	    ged_obol_path_has_semantic_prefix(candidate, target);
    if (match == GED_SCENE_PATH_MATCH_OBJECT)
	return ged_obol_path_has_component_name(candidate, target, 0) ||
	    ged_obol_path_has_semantic_component_name(candidate, target, 0);
    return 0;
}

static BObolCompactPathMatch
ged_obol_compact_path_match(enum ged_scene_path_match match)
{
    switch (match) {
	case GED_SCENE_PATH_MATCH_SUBTREE:
	    return BOBOL_COMPACT_PATH_SUBTREE;
	case GED_SCENE_PATH_MATCH_OBJECT:
	    return BOBOL_COMPACT_PATH_OBJECT;
	case GED_SCENE_PATH_MATCH_EXACT:
	default:
	    return BOBOL_COMPACT_PATH_EXACT;
    }
}


static int
ged_obol_database_source_in_transaction_scope(
    SoBRLDatabaseSource *source,
    const struct ged_draw_transaction *txn)
{
    if (!source || !txn)
	return 0;
    BObolDatabaseSourceSummary summary;
    return source->getSummary(summary) && summary.valid &&
	ged_obol_database_source_instance_in_scope(summary, txn->view) &&
	ged_obol_database_source_summary_matches_mode(summary, txn->mode);
}

static void
ged_obol_collect_transaction_instance_keys(
    std::vector<std::string> &instance_keys,
    BObolSceneController *scene,
    const struct ged_draw_transaction *txn,
    const char *target)
{
    if (!scene || !txn || !target)
	return;
    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, txn->view) ||
	    !ged_obol_database_source_summary_matches_mode(summary, txn->mode) ||
	    !ged_obol_transaction_path_matches(summary.path.getString(), target,
		txn->match))
	    continue;
	ged_obol_append_database_source_instance_key(instance_keys, summary);
    }
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
    const BObolCompactPathMatch compact_match =
	ged_obol_compact_path_match(txn->match);
    int handled = 0;
    for (const std::string &target : targets) {
	for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	    SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	    if (ged_obol_database_source_in_transaction_scope(source, txn) &&
		source->setCompactInstanceHighlightOverrideForPathMatch(
		    target.c_str(), compact_match,
		    highlighted ? TRUE : FALSE) > 0)
		handled = 1;
	}
	std::vector<std::string> instance_keys;
	ged_obol_collect_transaction_instance_keys(instance_keys, scene, txn,
	    target.c_str());
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
	    if (!ged_obol_transaction_path_matches(group_path, target.c_str(),
		txn->match))
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
    const BObolCompactPathMatch compact_match =
	ged_obol_compact_path_match(txn->match);
    int handled = 0;
    for (const std::string &target : targets) {
	for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	    SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	    if (ged_obol_database_source_in_transaction_scope(source, txn) &&
		source->setCompactInstanceVisibilityOverrideForPathMatch(
		    target.c_str(), compact_match,
		    visible ? TRUE : FALSE) > 0)
		handled = 1;
	}
	std::vector<std::string> instance_keys;
	ged_obol_collect_transaction_instance_keys(instance_keys, scene, txn,
	    target.c_str());
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
	    if (!ged_obol_transaction_path_matches(group_path, target.c_str(),
		txn->match))
		continue;
	    if (ged_obol_set_group_visible(scene, group_path, visible))
		handled = 1;
	}
    }

    return handled;
}

static int
ged_obol_apply_transparency_transaction(
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

    const float transparency = static_cast<float>(
	std::max(0.0, std::min(1.0, txn->value)));
    const BObolCompactPathMatch compact_match =
	ged_obol_compact_path_match(txn->match);
    int handled = 0;
    for (const std::string &target : targets) {
	for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	    SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	    if (ged_obol_database_source_in_transaction_scope(source, txn) &&
		source->setCompactInstanceTransparencyOverrideForPathMatch(
		    target.c_str(), compact_match, transparency) > 0)
		handled = 1;
	}

	std::vector<std::string> instance_keys;
	ged_obol_collect_transaction_instance_keys(instance_keys, scene, txn,
	    target.c_str());
	BObolDatabaseSourceDisplayPatch patch;
	patch.transparencyValid = TRUE;
	patch.transparency = transparency;
	for (const std::string &instance_key : instance_keys) {
	    if (scene->setDatabaseSourceInstanceDisplayPatch(
		    instance_key.c_str(), patch) >= 0)
		handled = 1;
	}

	const int summary_count = scene->getSceneTreeSummaryCount();
	for (int i = 0; i < summary_count; i++) {
	    BObolSceneTreeSummary tree_summary;
	    if (!scene->getSceneTreeSummary(i, tree_summary) ||
		!tree_summary.valid || tree_summary.nodeKind !=
		BObolSceneTreeSummary::NODE_GROUP ||
		tree_summary.path.getLength() == 0 ||
		BU_STR_EQUAL(tree_summary.path.getString(), "/"))
		continue;
	    const char *group_path = tree_summary.path.getString();
	    if (ged_obol_transaction_path_matches(group_path,
		    target.c_str(), txn->match) &&
		ged_obol_set_group_transparency(scene, group_path, transparency))
		handled = 1;
	}
    }
    return handled;
}


struct ged_obol_presentation_snapshot_context {
    struct ged_view_context *view_ctx;
    BObolSceneController *scene;
};


static int
ged_obol_frontier_presentation_apply_cb(
    const struct ged_draw_frontier_presentation_record *record,
    void *userdata)
{
    ged_obol_presentation_snapshot_context *ctx =
	static_cast<ged_obol_presentation_snapshot_context *>(userdata);
    if (!record || !record->path || !record->path[0] || !ctx || !ctx->scene)
	return 1;

    if (record->view && ctx->view_ctx && record->view != ctx->view_ctx)
	return 1;

    struct ged_draw_transaction txn = ged_draw_transaction_make(
	record->kind, record->path);
    txn.view = record->view;
    txn.mode = record->mode;
    txn.match = record->match;
    txn.value = record->value;
    switch (record->kind) {
	case GED_DRAW_TXN_VISIBILITY:
	    (void)ged_obol_apply_visibility_transaction(&txn, nullptr,
		ctx->scene);
	    break;
	case GED_DRAW_TXN_HIGHLIGHT:
	    (void)ged_obol_apply_highlight_transaction(&txn, nullptr,
		ctx->scene);
	    break;
	case GED_DRAW_TXN_TRANSPARENCY:
	    (void)ged_obol_apply_transparency_transaction(&txn, nullptr,
		ctx->scene);
	    break;
	default:
	    break;
    }
    return 1;
}


static void
ged_obol_frontier_presentation_snapshot_apply(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    BObolSceneController *scene)
{
    if (!gedp || !scene)
	return;
    ged_obol_presentation_snapshot_context ctx = {view_ctx, scene};
    (void)ged_draw_frontier_presentation_snapshot_foreach(gedp,
	ged_obol_frontier_presentation_apply_cb, &ctx);
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
	ged_view_context_obol_endpoint_get(view_ctx) : NULL;
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
     * a large workstation made first-frame latency worse through contention.
     * The byte governor controls memory, but it cannot recover UI CPU time. */
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
    if (!provider->setDatabase(dbip)) {
	delete provider;
	return 0;
    }
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
    /* Cache construction owns source arrays plus PoP topology/prefix buffers.
     * This task does not yet have decoded face/point counts, so reserve from
     * its serialized size instead of letting a zero estimate mean unlimited. */
    task.estimatedWorkingSetBytes =
	ged_obol_cold_object_working_set_bytes(dbip, dp, 16,
	    8 * 1024 * 1024);
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
			   ged_view_context_owner(view_ctx) :
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
			   ged_view_context_owner(view_ctx) :
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

extern "C" int
ged_database_path_member_autoview_bounds(
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

	/* Combine the actual finite member boxes first.  Applying framing to each
	 * member independently magnifies the margin by the number and aspect
	 * ratios of the parts. */
	VMIN(*min, child_min);
	VMAX(*max, child_max);
	have_bounds = 1;
    }

    if (have_bounds) {
	/* Match the retained-source framing contract: form one rotation-stable
	 * cube around the completed union.  This keeps a later aet change from
	 * clipping a long assembly without accumulating per-part padding. */
	point_t center;
	VADD2SCALE(center, *min, *max, 0.5);
	fastf_t extent = (*max)[X] - (*min)[X];
	V_MAX(extent, (*max)[Y] - (*min)[Y]);
	V_MAX(extent, (*max)[Z] - (*min)[Z]);
	if (extent < SQRT_SMALL_FASTF)
	    extent = SQRT_SMALL_FASTF;
	VSET(*min, center[X] - extent, center[Y] - extent,
	    center[Z] - extent);
	VSET(*max, center[X] + extent, center[Y] + extent,
	    center[Z] + extent);
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
		ged_database_path_member_autoview_bounds(gedp,
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
ged_draw_group_erase_subpath_for_path(
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
ged_obol_frontier_visibility_apply_cb(
    const struct ged_draw_frontier_visibility_record *record,
    void *userdata)
{
    struct ged *gedp = static_cast<struct ged *>(userdata);
    if (!gedp || !record || !record->root_path || !record->root_path[0])
	return 1;
    if (record->clear) {
	(void)ged_draw_obol_source_visibility_frontier_clear(gedp,
	    record->root_path, record->view, record->mode);
	return 1;
    }
    (void)ged_draw_obol_source_visibility_overrides_set(gedp,
	record->root_path, record->view, record->mode, record->paths,
	record->visible, record->rule_count);
    return 1;
}


static void
ged_obol_frontier_visibility_changes_apply(struct ged *gedp)
{
    if (!gedp || !ged_draw_obol_scene_controller_attached(gedp))
	return;
    (void)ged_draw_frontier_visibility_changes_foreach(gedp,
	ged_obol_frontier_visibility_apply_cb, gedp);
}


static void
ged_obol_frontier_visibility_snapshot_apply(struct ged *gedp)
{
    if (!gedp || !ged_draw_obol_scene_controller_attached(gedp))
	return;
    (void)ged_draw_frontier_visibility_snapshot_foreach(gedp,
	ged_obol_frontier_visibility_apply_cb, gedp);
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

	/* Compact occurrences own their individual selection presentation.  Once
	 * they exist, also selecting their aggregate source styles the retired
	 * overview/root proxy a second time and can make a white box reappear over
	 * an otherwise correct selected model. */
	const SbBool sourceSelected = source->getCompactInstanceCount() > 0 ?
	    FALSE : (ged_obol_path_in_selected_set(
		sourceSummary.path.getString(), paths, path_count) ? TRUE : FALSE);
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

struct ged_obol_selection_scene_collect_ctx {
    std::set<BObolSceneController *> *scenes;
};

static int
ged_obol_selection_scene_collect(struct ged_view_context *UNUSED(view_ctx),
    BObolViewController *controller, void *userdata)
{
    ged_obol_selection_scene_collect_ctx *ctx =
	static_cast<ged_obol_selection_scene_collect_ctx *>(userdata);
    if (ctx && ctx->scenes && controller && controller->getSceneController())
	ctx->scenes->insert(controller->getSceneController());
    return 1;
}

static void
ged_obol_selection_sources_for_path(BObolSceneController *scene,
    const char *path, std::set<SoBRLDatabaseSource *> &sources)
{
    if (!scene || !path || !path[0])
	return;

    std::string candidate(ged_obol_skip_leading_slash(path));
    while (!candidate.empty()) {
	SoBRLDatabaseSource *source =
	    scene->findDatabaseSource(candidate.c_str());
	if (source)
	    sources.insert(source);

	const int count =
	    scene->getDatabaseSourceInstanceCountForPath(candidate.c_str());
	for (int i = 0; i < count; i++) {
	    BObolDatabaseSourceSummary summary;
	    if (!scene->getDatabaseSourceInstanceSummaryForPath(
		    candidate.c_str(), i, summary) || !summary.valid ||
		summary.instanceKey.getLength() <= 0)
		continue;
	    source = scene->findDatabaseSourceInstance(
		summary.instanceKey.getString());
	    if (source)
		sources.insert(source);
	}

	const size_t slash = candidate.rfind('/');
	if (slash == std::string::npos)
	    break;
	candidate.erase(slash);
    }
}

extern "C" int
ged_draw_obol_database_sources_apply_selection_delta(
    struct ged *gedp,
    const char *const *added_paths,
    size_t added_count,
    const char *const *removed_paths,
    size_t removed_count,
    const char *const *selected_paths,
    size_t selected_count)
{
    if (!gedp || !gedp->i || !gedp->i->ged_gdp)
	return 0;

    std::set<BObolSceneController *> scenes;
    BObolSceneController *owned = ged_draw_obol_scene_controller(gedp);
    if (owned)
	scenes.insert(owned);
    ged_obol_selection_scene_collect_ctx scene_ctx = {&scenes};
    ged_bobol_view_controllers_foreach(gedp,
	ged_obol_selection_scene_collect, &scene_ctx);

    std::vector<SbString> added;
    std::vector<SbString> removed;
    added.reserve(added_count);
    removed.reserve(removed_count);
    for (size_t i = 0; added_paths && i < added_count; i++) {
	if (added_paths[i] && added_paths[i][0])
	    added.push_back(SbString(added_paths[i]));
    }
    for (size_t i = 0; removed_paths && i < removed_count; i++) {
	if (removed_paths[i] && removed_paths[i][0])
	    removed.push_back(SbString(removed_paths[i]));
    }

    int applied = 0;
    for (BObolSceneController *scene : scenes) {
	std::set<SoBRLDatabaseSource *> sources;
	for (size_t i = 0; added_paths && i < added_count; i++)
	    ged_obol_selection_sources_for_path(scene, added_paths[i], sources);
	for (size_t i = 0; removed_paths && i < removed_count; i++)
	    ged_obol_selection_sources_for_path(scene, removed_paths[i], sources);

	for (SoBRLDatabaseSource *source : sources) {
	    BObolDatabaseSourceSummary summary;
	    if (!source || !source->getSummary(summary) || !summary.valid)
		continue;
	    applied += source->applyCompactInstanceSelectionDelta(added,
		removed);
	    const SbBool selected = source->getCompactInstanceCount() > 0 ?
		FALSE : (ged_obol_path_in_selected_set(summary.path.getString(),
		    selected_paths, selected_count) ? TRUE : FALSE);
	    if (summary.selected != selected) {
		const int changed = source->setDisplayState(FALSE,
		    summary.sourceRevision, summary.inputsRevision,
		    summary.visible, selected, summary.highlighted,
		    summary.lineStyle, summary.lineWidth,
		    summary.transparency, summary.colorOverride,
		    summary.color, summary.materialColorValid,
		    summary.materialColor, summary.materialRevision);
		if (changed > 0)
		    applied++;
	    }
	}
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
    const int published = scene->publishDatabaseSourceInstanceExternalLineSet(
	       source_instance_key, line_set);
    /* The top-level draw-proxy cache record is the complete path AABB, not an
     * append-only leaf union.  Publish the same exact fact through the source
     * bounds contract so an autoview transaction can frame this useful
     * overview immediately instead of leaving an off-center box invisible
     * until detached leaf discovery finishes.  Do this even when the line set
     * was already present: an older retained source may have identical proxy
     * geometry but no explicit exact-bounds witness. */
    (void)scene->setDatabaseSourceInstanceBoundsState(
	source_instance_key, TRUE,
	SbVec3f(static_cast<float>(bmin[X]),
	    static_cast<float>(bmin[Y]),
	    static_cast<float>(bmin[Z])),
	SbVec3f(static_cast<float>(bmax[X]),
	    static_cast<float>(bmax[Y]),
	    static_cast<float>(bmax[Z])));
    return published;
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
ged_obol_aabb_proxy_geometry(const point_t bounds_min, const point_t bounds_max,
    SbMatrix *geometry_transform)
{
    if (!bounds_min || !bounds_max)
	return std::shared_ptr<const Obol::PartGeometry>();

    if (geometry_transform)
	*geometry_transform = SbMatrix::identity();

    /*
     * Structural AABBs are instances of one twelve-edge primitive.  Retain
     * one unit cube and keep its geometry-to-object transform distinct from
     * the hierarchy placement.  The distinction is required when a box is
     * upgraded in place to native object-local mesh data.
     */
    const double sx = bounds_max[X] - bounds_min[X];
    const double sy = bounds_max[Y] - bounds_min[Y];
    const double sz = bounds_max[Z] - bounds_min[Z];
    if (geometry_transform && isfinite(sx) && isfinite(sy) && isfinite(sz) &&
	sx > SMALL_FASTF && sy > SMALL_FASTF && sz > SMALL_FASTF) {
	static const std::shared_ptr<const Obol::PartGeometry> unit_geometry = []() {
	    std::shared_ptr<Obol::PartGeometry> geometry(
		new Obol::PartGeometry);
	    Obol::WireRep wire;
	    const SbVec3f corners[8] = {
		SbVec3f(0.0f, 0.0f, 0.0f), SbVec3f(1.0f, 0.0f, 0.0f),
		SbVec3f(0.0f, 1.0f, 0.0f), SbVec3f(1.0f, 1.0f, 0.0f),
		SbVec3f(0.0f, 0.0f, 1.0f), SbVec3f(1.0f, 0.0f, 1.0f),
		SbVec3f(0.0f, 1.0f, 1.0f), SbVec3f(1.0f, 1.0f, 1.0f)
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
	    geometry->subpixelProxyEligible = true;
	    geometry->structuralProxy = true;
	    return std::shared_ptr<const Obol::PartGeometry>(geometry);
	}();
	geometry_transform->setScale(
	    SbVec3f(static_cast<float>(sx), static_cast<float>(sy),
		static_cast<float>(sz)));
	SbMatrix translation;
	translation.setTranslate(
	    SbVec3f(static_cast<float>(bounds_min[X]),
		static_cast<float>(bounds_min[Y]),
		static_cast<float>(bounds_min[Z])));
	geometry_transform->multRight(translation);
	return unit_geometry;
    }

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
    geometry->structuralProxy = true;
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
	node.boundsMax, &occurrence.geometryTransform);
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
    summary.lodActiveCut = BOBOL_LOD_QUALITY_PROXY;
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
	request.sourceType = cached.meshAssetKind ==
	    BOBOL_DRAW_CACHE_MESH_ASSET_BREP ? "brep" : "bot";
	request.meshAssetPath = cached.meshAssetPath;
	request.meshAssetName = cached.meshAssetName;
	request.meshAssetContentHash = cached.meshAssetContentHash;
	request.meshAssetTessellationAbsTol =
	    cached.meshAssetTessellationAbsTol;
	request.meshAssetTessellationRelTol =
	    cached.meshAssetTessellationRelTol;
	request.meshAssetTessellationNormTol =
	    cached.meshAssetTessellationNormTol;
	request.meshAssetBounds = SbBox3f(
	    SbVec3f(static_cast<float>(cached.meshAssetBoundsMin[X]),
		static_cast<float>(cached.meshAssetBoundsMin[Y]),
		static_cast<float>(cached.meshAssetBoundsMin[Z])),
	    SbVec3f(static_cast<float>(cached.meshAssetBoundsMax[X]),
		static_cast<float>(cached.meshAssetBoundsMax[Y]),
		static_cast<float>(cached.meshAssetBoundsMax[Z])));
	request.meshAssetTransform =
	    ged_obol_sbmatrix_from_mat(cached.meshAssetMatrix);
	request.faceCount = cached.sourceFaceCount;
	request.pointCount = cached.sourcePointCount;
	request.bounds = occurrence.summary.bounds;
	request.databaseIntent = 1;
	request.localSource = 1;
	request.drawMode = occurrence.summary.drawMode;
	request.recordRole = "lod-source";
	request.geometryKind = "surface";
	/* The manifest summary contains the authoritative effective material and
	 * region state for this occurrence.  Keep the source-backed refinement
	 * request semantically identical to its standing proxy: otherwise the
	 * first PoP/BREP replacement loses full-path color inheritance and falls
	 * back to the aggregate source material. */
	request.sourceId = occurrence.summary.sourceId;
	request.sourceIdentity = occurrence.summary.sourceIdentity;
	request.regionId = occurrence.summary.regionId;
	request.airCode = occurrence.summary.airCode;
	request.materialId = occurrence.summary.materialId;
	request.los = occurrence.summary.los;
	request.materialColorValid =
	    occurrence.summary.materialColorValid ? 1 : 0;
	request.materialColor = occurrence.summary.materialColor;
	request.materialShader = occurrence.summary.materialShader;
	request.colorOverride = occurrence.summary.colorOverride ? 1 : 0;
	request.color = occurrence.summary.color;
	request.transparency = occurrence.summary.transparency;
	request.lodAvailable = 1;
	const SbVec3f bmin = request.bounds.getMin();
	const SbVec3f bmax = request.bounds.getMax();
	request.lodBoundsMin = bmin;
	request.lodBoundsMax = bmax;
    }
    return occurrence;
}

/* A final leaf registry is representation data, not merely hierarchy data.
 * Keep it out of the structural-manifest namespace and distinguish every
 * source policy which can change the immutable wire/shaded asset contract.
 * The draw-cache payload still carries and validates this complete identity,
 * so a hash collision remains a cache miss rather than a correctness risk. */
static std::string
ged_obol_leaf_manifest_cache_identity(const SoBRLDatabaseSource *source,
    const std::string &rootPath)
{
    if (!source || rootPath.empty())
	return std::string();
    char variant[320] = {0};
    snprintf(variant, sizeof(variant),
	"|leaf-v2|draw=%d|representation=%d|bot=%u|tess=%.17g,%.17g,%.17g",
	source->drawMode.getValue(), source->representationMode.getValue(),
	source->lodBotThreshold.getValue(),
	static_cast<double>(source->tessellationAbsTol.getValue()),
	static_cast<double>(source->tessellationRelTol.getValue()),
	static_cast<double>(source->tessellationNormTol.getValue()));
    return rootPath + variant;
}

struct ged_obol_warm_manifest_stream_context {
    ged_obol_structural_proxy_context *proxyContext = NULL;
    BObolCompactOccurrenceStream *stream = NULL;
    std::string rootPath;
    std::string cacheName;
    bool exactOverviewQueued = false;
    bool complete = false;
    size_t pushed = 0;
};

static int
ged_obol_warm_manifest_begin(const BObolDrawManifest *manifest,
	void *userData)
{
    ged_obol_warm_manifest_stream_context *context =
	static_cast<ged_obol_warm_manifest_stream_context *>(userData);
    if (!context || !context->proxyContext || !context->stream ||
	!manifest || context->stream->isCancelled())
	return 0;
    context->stream->setExpectedCount(manifest->occurrenceCount);
    if (!manifest->coverageBoundsValid)
	return 1;

    const SbBox3f exactBounds(
	SbVec3f(static_cast<float>(manifest->coverageBoundsMin[X]),
	    static_cast<float>(manifest->coverageBoundsMin[Y]),
	    static_cast<float>(manifest->coverageBoundsMin[Z])),
	SbVec3f(static_cast<float>(manifest->coverageBoundsMax[X]),
	    static_cast<float>(manifest->coverageBoundsMax[Y]),
	    static_cast<float>(manifest->coverageBoundsMax[Z])));
    if (exactBounds.isEmpty())
	return 1;

    context->stream->setCoverageBounds(exactBounds);
    ged_obol_structural_proxy_node node;
    /* The overview is presentation data for the real draw root, not a
     * synthetic database descendant.  Keep recordRole/geometryKind as its
     * internal discriminator; a fake path is rejected by path visibility and
     * material resolution in some frames. */
    node.path = context->rootPath;
    node.objectName = context->cacheName;
    node.instanceKey = node.path;
    node.boolOp = DB_OP_UNION;
    node.row = 0;
    node.localMatrix = SbMatrix::identity();
    VMOVE(node.boundsMin, manifest->coverageBoundsMin);
    VMOVE(node.boundsMax, manifest->coverageBoundsMax);
    node.publishBounds = 1;
    node.metadataValid = 0;
    BObolCompactOccurrence overview =
	ged_obol_structural_proxy_occurrence(*context->proxyContext, node);
    overview.summary.recordRole = "lod-overview";
    overview.summary.geometryKind = "overview-aabb";
    overview.summary.selectable = FALSE;
    if (!overview.geometry)
	return 1;
    context->stream->pushPriority(overview);
    /* The persisted extent is exact and the priority lane is drained before
     * ordinary leaves.  This establishes correct autoview and screen-space
     * demand while the cache mapping is still streaming its first leaves. */
    context->stream->setCoverageBoundsComplete(true);
    context->exactOverviewQueued = true;
    context->complete = manifest->occurrenceCount > 0;
    return 1;
}

static int
ged_obol_warm_manifest_occurrence(
	const BObolDrawManifestOccurrence *cached, size_t occurrenceIndex,
	void *userData)
{
    (void)occurrenceIndex;
    ged_obol_warm_manifest_stream_context *context =
	static_cast<ged_obol_warm_manifest_stream_context *>(userData);
    if (!context || !context->proxyContext || !context->stream || !cached ||
	context->stream->isCancelled())
	return 0;
    BObolCompactOccurrence occurrence =
	ged_obol_structural_proxy_manifest_occurrence(
	    *context->proxyContext, *cached);
    if (!occurrence.geometry) {
	context->complete = false;
	return 1;
    }
    if (!occurrence.sourceMeshRequestValid)
	context->complete = false;
    context->stream->push(std::move(occurrence));
    context->pushed++;
    return 1;
}

/*
 * Warm leaf coverage is a producer-side concern.  Deserialize and enqueue it
 * on the detached worker rather than constructing tens of thousands of CAD
 * records in the command callback.  The live pump retains its bounded merge
 * deadline, so a ready 50k manifest cannot monopolize input handling.
 *
 * This is a presentation seed.  The following database walk still produces
 * authoritative material and hierarchy semantics, but it can skip the
 * redundant full-BoT coverage import when every manifest occurrence carries a
 * source-mesh request.
 */
static int
ged_obol_stream_cached_leaf_manifest(
    SoBRLDatabaseSource *source,
    struct db_i *database,
    int draw_mode,
    uint32_t source_revision,
    BObolCompactOccurrenceStream *stream,
    void *user_data)
{
    (void)user_data;
    if (!source || !database || !stream || stream->isCancelled())
	return 0;
    const std::string rootPath = ged_obol_skip_leading_slash(
	source->path.getValue().getString());
    if (rootPath.empty())
	return 0;

    const int64_t started = bu_gettime();
    ged_obol_structural_proxy_context ctx{};
    ctx.gedp = NULL;
    ctx.viewCtx = NULL;
    ctx.drawMode = draw_mode;
    ctx.sourceRevision = source_revision;

    const std::string manifestIdentity =
	ged_obol_leaf_manifest_cache_identity(source, rootPath);
    if (manifestIdentity.empty())
	return 0;
    const char *cacheName = strrchr(rootPath.c_str(), '/');
    cacheName = cacheName ? cacheName + 1 : rootPath.c_str();
    ged_obol_warm_manifest_stream_context streamContext;
    streamContext.proxyContext = &ctx;
    streamContext.stream = stream;
    streamContext.rootPath = rootPath;
    streamContext.cacheName = cacheName;
    const int streamStatus = bobol_draw_manifest_cache_stream(database,
	manifestIdentity.c_str(), ged_obol_warm_manifest_begin,
	ged_obol_warm_manifest_occurrence, &streamContext);
    const bool complete = streamStatus == BRLCAD_OK &&
	streamContext.complete && streamContext.exactOverviewQueued &&
	streamContext.pushed;
    if (complete)
	stream->setWarmCoverageComplete(true);
    ged_obol_timing_log(complete ?
	"job: stream warm leaf manifest" :
	"job: stream partial leaf manifest",
	started, static_cast<long>(streamContext.pushed));
    return complete ? 2 : (streamContext.pushed ? 1 : 0);
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
    struct db_i *database,
    const SoBRLDatabaseSource *source,
    BObolCompactOccurrenceStream *stream,
    void *user_data)
{
    (void)user_data;
    if (!database || !source)
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
	    occurrence.summary.path.getLength() == 0 ||
	    BU_STR_EQUAL(occurrence.summary.recordRole.getString(),
		"lod-overview"))
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

    SbBox3f exactCoverageBounds;
    if (stream && stream->getCoverageBounds(exactCoverageBounds)) {
	manifest.coverageBoundsValid = 1;
	const SbVec3f exactMin = exactCoverageBounds.getMin();
	const SbVec3f exactMax = exactCoverageBounds.getMax();
	VSET(manifest.coverageBoundsMin, exactMin[0], exactMin[1],
	    exactMin[2]);
	VSET(manifest.coverageBoundsMax, exactMax[0], exactMax[1],
	    exactMax[2]);
    }

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
	cached.booleanOperation =
	    occurrence.booleanOperation ==
		SoBRLDatabaseSource::BOOLEAN_SUBTRACT ? DB_OP_SUBTRACT :
	    (occurrence.booleanOperation ==
		SoBRLDatabaseSource::BOOLEAN_INTERSECT ? DB_OP_INTERSECT :
		DB_OP_UNION);
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
	    cached.meshAssetKind = BU_STR_EQUAL(
		request.sourceType.getString(), "brep") ?
		BOBOL_DRAW_CACHE_MESH_ASSET_BREP :
		BOBOL_DRAW_CACHE_MESH_ASSET_BOT;
	    cached.meshAssetContentHash = request.meshAssetContentHash;
	    cached.meshAssetTessellationAbsTol =
		request.meshAssetTessellationAbsTol;
	    cached.meshAssetTessellationRelTol =
		request.meshAssetTessellationRelTol;
	    cached.meshAssetTessellationNormTol =
		request.meshAssetTessellationNormTol;
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
	    ged_obol_mat_from_sbmatrix(request.meshAssetTransform,
		cached.meshAssetMatrix);
	}
	/*
	 * A complete leaf manifest is now also the authoritative warm-start
	 * registry.  Preserve the effective path material values already resolved
	 * by the cold database walk so the worker need not repeat a full geometry
	 * import merely to recover colors and region semantics.
	 */
	bobol_draw_metadata_record_init(&cached.metadata);
	cached.metadataValid = 1;
	cached.metadata.hasRegionId = 1;
	cached.metadata.regionId = occurrence.summary.regionId;
	cached.metadata.hasAircode = 1;
	cached.metadata.aircode = occurrence.summary.airCode;
	cached.metadata.hasLos = 1;
	cached.metadata.los = occurrence.summary.los;
	cached.metadata.hasMaterialId = 1;
	cached.metadata.materialId = occurrence.summary.materialId;
	if (occurrence.summary.materialColorValid) {
	    const SbColor color = occurrence.summary.materialColor;
	    cached.metadata.hasColor = 1;
	    for (size_t channel = 0; channel < 3; channel++) {
		const float component = std::max(0.0f,
		    std::min(1.0f, color[static_cast<int>(channel)]));
		cached.metadata.color[channel] =
		    static_cast<unsigned char>(component * 255.0f + 0.5f);
	    }
	}
	const char *shader =
	    occurrence.summary.materialShader.getString();
	if (shader && shader[0]) {
	    cached.metadata.hasShader = 1;
	    bu_strlcpy(cached.metadata.shader, shader,
		sizeof(cached.metadata.shader));
	}
    }

    const std::string manifestIdentity =
	ged_obol_leaf_manifest_cache_identity(source, rootPath);
    const int stored = valid && !manifestIdentity.empty() &&
	bobol_draw_manifest_cache_store(database, manifestIdentity.c_str(),
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

    /* A warm draw reads only the fixed-size manifest description here.  Its
     * exact coverage extent is enough to publish an immediate whole-target
     * overview; the detached worker retains ownership of occurrence decoding
     * and streams the leaf records through the bounded provider. */
    if (cachedManifestOnly) {
	if (root->getCompactInstanceCount() > 0)
	    return 0;
	const std::string manifestIdentity =
	    ged_obol_leaf_manifest_cache_identity(root, rootPath);
	struct BObolDrawManifest description;
	bobol_draw_manifest_init(&description);
	if (manifestIdentity.empty() ||
	    bobol_draw_manifest_cache_describe(gedp->dbip,
		manifestIdentity.c_str(), &description) != BRLCAD_OK ||
	    !description.coverageBoundsValid ||
	    !description.occurrenceCount)
	    return 0;

	const char *cacheName = strrchr(rootPath.c_str(), '/');
	cacheName = cacheName ? cacheName + 1 : rootPath.c_str();
	ged_obol_structural_proxy_node node;
	node.path = rootPath;
	node.objectName = cacheName;
	node.instanceKey = node.path;
	node.boolOp = DB_OP_UNION;
	node.row = 0;
	node.localMatrix = SbMatrix::identity();
	VMOVE(node.boundsMin, description.coverageBoundsMin);
	VMOVE(node.boundsMax, description.coverageBoundsMax);
	node.publishBounds = 1;
	node.metadataValid = 0;
	BObolCompactOccurrence overview =
	    ged_obol_structural_proxy_occurrence(ctx, node);
	overview.summary.recordRole = "lod-overview";
	overview.summary.geometryKind = "overview-aabb";
	overview.summary.selectable = FALSE;
	if (!overview.geometry)
	    return 0;

	std::vector<BObolCompactOccurrence> occurrences;
	occurrences.push_back(std::move(overview));
	ged_obol_scene_mutation_batch_scope batch(scene, 1, 1);
	if (root->setCompactOccurrenceRegistry(occurrences) <= 0)
	    return 0;
	(void)root->setSourceBoundsState(TRUE,
	    SbVec3f(static_cast<float>(description.coverageBoundsMin[X]),
		static_cast<float>(description.coverageBoundsMin[Y]),
		static_cast<float>(description.coverageBoundsMin[Z])),
	    SbVec3f(static_cast<float>(description.coverageBoundsMax[X]),
		static_cast<float>(description.coverageBoundsMax[Y]),
		static_cast<float>(description.coverageBoundsMax[Z])));
	return ged_obol_database_source_mark_published_current(scene, root) ?
	    1 : 0;
    }

    struct BObolDrawManifest manifest;
    bobol_draw_manifest_init(&manifest);
    const int64_t manifest_load_start = bu_gettime();
    const int manifest_status = bobol_draw_manifest_cache_get(gedp->dbip,
	rootPath.c_str(), &manifest);
    ged_obol_timing_log(manifest_status == BRLCAD_OK ?
	"structural: leaf manifest hit" : "structural: leaf manifest miss",
	manifest_load_start,
	manifest_status == BRLCAD_OK ?
	    static_cast<long>(manifest.occurrenceCount) : 0);
    if (manifest_status == BRLCAD_OK) {
	const int64_t manifest_publish_start = bu_gettime();
	std::vector<BObolCompactOccurrence> occurrences;
	occurrences.reserve(manifest.occurrenceCount);
	bool allOccurrencesMeshBacked = manifest.occurrenceCount > 0;
	for (size_t i = 0; i < manifest.occurrenceCount; i++) {
	    BObolCompactOccurrence occurrence =
		ged_obol_structural_proxy_manifest_occurrence(ctx,
		    manifest.occurrences[i]);
	    if (!occurrence.geometry) {
		occurrences.clear();
		break;
	    }
	    if (!occurrence.sourceMeshRequestValid)
		allOccurrencesMeshBacked = false;
	    occurrences.push_back(std::move(occurrence));
	}
	bobol_draw_manifest_free(&manifest);
	if (!occurrences.empty()) {
	    ged_obol_scene_mutation_batch_scope batch(scene, 1,
		occurrences.size());
	    if (root->setCompactOccurrenceRegistry(occurrences) > 0) {
		ged_obol_timing_log("structural: publish leaf manifest",
		    manifest_publish_start,
		    static_cast<long>(occurrences.size()));
		if (!ged_obol_database_source_mark_published_current(scene, root))
		    return 0;
		/*
		 * A complete leaf-bound manifest is only terminal for this
		 * realization stage when every occurrence carries an independent
		 * mesh-source request.  Otherwise its non-mesh AABBs still need
		 * the detached native-geometry pass.  Treating structural
		 * completeness as drawing-data completeness made mixed models
		 * (notably Generic_Twin) work cold but remain boxes forever warm.
		 */
		return allOccurrencesMeshBacked ? 2 : 1;
	    }
	}
	ged_obol_timing_log("structural: reject leaf manifest",
	    manifest_publish_start, static_cast<long>(occurrences.size()));
    }

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

    /*
     * The aggregate is deliberately a separate small cache record.  Loading
     * and constructing a 50k-occurrence leaf manifest on the UI thread took
     * seconds even though its whole-target extent was already known.  Publish
     * that O(1) overview first; its ordinary progressive job will stream the
     * leaf manifest/authoritative data without blocking command input.
     */
    const char *cache_name = strrchr(path, '/');
    cache_name = cache_name ? cache_name + 1 : path;
    struct ged_draw_obol_source_expansion_status proxy_status;
    ged_obol_source_expansion_status_clear(&proxy_status);
    if (cache_name[0] &&
	ged_obol_publish_aabb_proxy_for_path(gedp, view_ctx, path,
	    source_instance_key, cache_name, draw_mode, 0, &proxy_status))
	return 1;

    const int structural_snapshot =
	ged_obol_publish_deferred_structural_proxy_snapshot(gedp, view_ctx,
	    path, source_instance_key, draw_mode, 0, 0, 0, 1);
    if (structural_snapshot != 0)
	return structural_snapshot;

    /* A cache miss is intentionally not refreshed here: recursive bounds work
     * belongs on a worker, and detached realization shortly publishes an
     * expanding aggregate followed by per-leaf presentation. */
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

    if (!provider->setDatabase(dbip)) {
	delete provider;
	return 0;
    }
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
    struct directory *dp = db_lookup(dbip, child_name, LOOKUP_QUIET);
    if (dp && (dp->d_flags & RT_DIR_COMB)) {
	/* A recursive combination bound has no useful reservation in the tiny
	 * combination record itself.  Serialize it against other expensive work;
	 * rt_obj_bounds releases descendant internals as it walks them. */
	const size_t limit = service->getWorkingSetLimit();
	task.estimatedWorkingSetBytes =
	    limit && limit != SIZE_MAX ? limit : 1024ULL * 1024ULL * 1024ULL;
    } else {
	/* Bounding a BoT still deserializes its complete vertex/face arrays. */
	task.estimatedWorkingSetBytes =
	    ged_obol_cold_object_working_set_bytes(dbip, dp, 4,
		2 * 1024 * 1024);
    }
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
ged_obol_deferred_job_coverage_bounds_complete(
    const std::shared_ptr<ged_obol_deferred_realization_job> &job)
{
    if (!job || job->items.empty())
	return 0;
    for (const std::unique_ptr<ged_obol_deferred_realization_item> &item :
	 job->items) {
	if (!item || !item->stream ||
	    !item->stream->hasCoverageBoundsComplete())
	    return 0;
    }
    return 1;
}

static int
ged_obol_progressive_autoview_apply(
    ged_obol_progressive_provider_data *data)
{
    if (!data || !data->gedp || !data->view_ctx ||
	!data->pending_autoview ||
	!data->pending_autoview_bounds_complete)
	return 0;

    struct bv *view = ged_obol_bv(data->view_ctx);
    if (!view) {
	data->pending_autoview = 0;
	data->pending_autoview_bounds_complete = 0;
	return 0;
    }
    if (bv_frame_revision_get(view) != data->expected_view_revision) {
	data->pending_autoview = 0;
	data->pending_autoview_bounds_complete = 0;
	return 0;
    }

    vect_t bmin;
    vect_t bmax;
    int empty = 1;
    if (!ged_draw_obol_scene_database_autoview_bounds(data->gedp, &bmin,
	&bmax, &empty, 0) || empty) {
	return 0;
	}

    if (!bv_autoview_bounds(view, data->autoview_factor, bmin, bmax)) {
	return 0;
	}
    data->expected_view_revision = bv_frame_revision_get(view);
    bv_refresh_request(view, GED_VIEW_REFRESH_DRAW);
    /*
     * Exact coverage bounds are immutable for this source revision.  One
     * successful application fulfills the deferred request; leaving it armed
     * made every later leaf/mesh publication rewrite an identical camera and
     * visibly flash some GL backends.
     */
    data->pending_autoview = 0;
    data->pending_autoview_bounds_complete = 0;
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
    data->pending_autoview_bounds_complete = 0;
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
    const int64_t now = bu_gettime();
    data->retired_jobs.erase(std::remove_if(data->retired_jobs.begin(),
	data->retired_jobs.end(),
	[now](const std::shared_ptr<ged_obol_deferred_realization_job> &job) {
	    if (!job)
		return true;
	    const int state = job->stateValue();
	    if (state == ged_obol_deferred_realization_job::COMPLETE &&
		job->stagingLeaseReleaseAfter > now)
		return false;
	    return state == ged_obol_deferred_realization_job::COMPLETE ||
		   state == ged_obol_deferred_realization_job::FAILED ||
		   state == ged_obol_deferred_realization_job::CANCELLED;
	}), data->retired_jobs.end());
}

static void
ged_obol_hold_completed_job(
    ged_obol_progressive_provider_data *data,
    const std::shared_ptr<ged_obol_deferred_realization_job> &job)
{
    if (!data || !job)
	return;
    job->stagingLeaseReleaseAfter = bu_gettime() + 500000;
    data->retired_jobs.push_back(job);
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

/* A deferred source may be attached after the draw transaction's ordinary
 * selection-delta synchronization has completed.  Initialize that new
 * renderer owner from the semantic selection service before any streamed
 * occurrence arrives.  The source retains this path rule while its compact
 * index is empty and applies it incrementally to later batches, avoiding both
 * a missed selected occurrence and an O(N) catch-up scan. */
static void
ged_obol_initialize_deferred_source_selection(struct ged *gedp,
	SoBRLDatabaseSource *source)
{
    if (!gedp || !source || !ged_selection_state_available(gedp))
	return;

    struct bu_vls listed = BU_VLS_INIT_ZERO;
    (void)ged_selection_list_paths(gedp, NULL, &listed);
    const char *text = bu_vls_cstr(&listed);
    std::vector<SbString> paths;
    if (text && text[0]) {
	const char *begin = text;
	for (const char *cursor = text;; ++cursor) {
	    if (*cursor != '\n' && *cursor != '\0')
		continue;
	    if (cursor > begin)
		paths.push_back(SbString(std::string(begin,
		    static_cast<size_t>(cursor - begin)).c_str()));
	    if (*cursor == '\0')
		break;
	    begin = cursor + 1;
	}
    }
    (void)source->syncCompactInstanceSelectedPaths(paths);
    bu_vls_free(&listed);
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
	ged_obol_initialize_deferred_source_selection(data->gedp, live);
	std::unique_ptr<ged_obol_deferred_realization_item> item(
	    new ged_obol_deferred_realization_item);
	struct db_i *liveDatabase = live->getDatabase();
	item->source = live->createDetachedRealizationTemplate();
	item->ownsSource = item->source ? TRUE : FALSE;
	item->snapshotSourceDatabase = liveDatabase ?
	    db_clone_dbi(liveDatabase, NULL) : NULL;
	if (!item->source || !item->snapshotSourceDatabase) {
	    bu_log("Obol deferred realization could not capture source '%s'\n",
		key.c_str());
	    return 0;
	}
	/* These fields are the immutable result-admission stamp.  The detached
	 * Coin source becomes worker-exclusive after start(); the GUI streaming
	 * pump must not read it again until the job's release-published terminal
	 * state makes final adoption safe. */
	item->sourceRevision = item->source->sourceRevision.getValue();
	item->inputsRevision = item->source->inputsRevision.getValue();
	item->viewRevision = item->source->viewRevision.getValue();
	item->representationMode = item->source->representationMode.getValue();
	item->instanceKey = live->instanceKey.getValue().getString();
	item->drawMode = drawMode;
	item->sourceDrawMode = item->source->drawMode.getValue();
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

static int ged_obol_apply_stream_coverage_bounds(
    SoBRLDatabaseSource *source,
    BObolCompactOccurrenceStream *stream);

static int
ged_obol_publish_deferred_realization(
    ged_obol_progressive_provider_data *data,
    BObolViewController *controller,
    const std::shared_ptr<ged_obol_deferred_realization_job> &job)
{
    if (!data || !data->gedp || !controller || !job)
	return 0;
    if (job->stateValue() !=
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
	if (live->sourceRevision.getValue() != item->sourceRevision ||
	    live->inputsRevision.getValue() != item->inputsRevision ||
	    live->drawMode.getValue() != item->sourceDrawMode ||
	    live->representationMode.getValue() !=
		item->representationMode) {
	    if (ged_obol_timing_enabled()) {
		bu_log("[obol-timing] deferred adoption rejected for %s: "
		    "source=%u/%u inputs=%u/%u view(non-binding)=%u/%u "
		    "mode=%d/%d "
		    "representation=%d/%d\n",
		    item->instanceKey.c_str(),
		    live->sourceRevision.getValue(),
		    item->sourceRevision,
		    live->inputsRevision.getValue(),
		    item->inputsRevision,
		    live->viewRevision.getValue(),
		    item->viewRevision,
		    live->drawMode.getValue(), item->sourceDrawMode,
		    live->representationMode.getValue(),
		    item->representationMode);
	    }
	    return 0;
	}
	liveSources.push_back(live);
    }

    int adopted = 0;
    for (size_t i = 0; i < liveSources.size(); i++) {
	const int adopted_count = liveSources[i]->adoptDetachedCompactRealization(
	    job->items[i]->source, TRUE);
	/* Detached adoption rebuilds the registry-derived union.  Preserve the
	 * stream's representation-independent path extent for future explicit
	 * autoviews as well as the pending progressive fit. */
	(void)ged_obol_apply_stream_coverage_bounds(liveSources[i],
	    job->items[i]->stream.get());
	if (ged_obol_timing_enabled())
	    bu_log("[obol-timing] deferred adoption for %s: n=%d\n",
		job->items[i]->instanceKey.c_str(), adopted_count);
	if (adopted_count > 0) {
	    adopted++;
	    BObolSourceRealizationItemResult realization_result;
	    const SbBool have_realization_result = job->realization &&
		job->realization->itemResult(i, realization_result);
	    if (ged_obol_timing_enabled()) {
		const SbBool warm_manifest = have_realization_result ?
		    realization_result.warmManifest : FALSE;
		const SbBool manifest_stored = have_realization_result ?
		    realization_result.manifestStored : FALSE;
		bu_log("[obol-timing] deferred manifest: %s n=%d\n",
		    warm_manifest ? "warm-reused" :
		    (manifest_stored ?
			"worker-stored" : "worker-store-failed"),
		    adopted_count);
	    }
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
/* One provider is normally one large retained root.  A fixed 64-item quantum
 * made the 16 ms host timer, rather than actual work, the throughput limit;
 * a fixed 512-item quantum occasionally exceeded an input frame when later
 * merges had become more expensive.  Start each stream at 256 and adapt
 * between 64 and this ceiling from its measured atomic merge time. */
static const size_t GED_OBOL_STREAM_BATCH_QUANTUM_MAX = 2048;

/*
 * A compact registry's derived union is useful while no stronger information
 * exists, but it is not the autoview contract.  Structural and mesh records
 * may use asset-local bounds plus placement transforms, and they arrive in a
 * renderer-dependent order.  Recomputing sourceBounds from each partial
 * registry therefore made a warm, fast renderer occasionally frame a
 * different extent than a cold or slower renderer.
 *
 * The realization stream publishes one immutable path-scoped coverage bound.
 * Once available, keep that value authoritative across every incremental
 * merge and final detached adoption.  setSourceBoundsState stores source-local
 * bounds; getEffectiveSourceBounds applies the source placement exactly once.
 */
static int
ged_obol_apply_stream_coverage_bounds(
    SoBRLDatabaseSource *source,
    BObolCompactOccurrenceStream *stream)
{
    if (!source || !stream || !stream->hasCoverageBoundsComplete())
	return 0;

    SbBox3f bounds;
    if (!stream->getCoverageBounds(bounds) || bounds.isEmpty())
	return 0;

    (void)source->setSourceBoundsState(TRUE, bounds.getMin(), bounds.getMax());
    (void)source->setSourceBoundsExactState(TRUE);
    return 1;
}

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
	    if (!live || !item->source)
		continue;
	    /* A source/draw-mode change after this job launched makes its stream
	     * stale; the job will be retired by the pump's revision drain.  Do not
	     * merge stale geometry onto the current live source. */
	    if (live->sourceRevision.getValue() != item->sourceRevision ||
		live->inputsRevision.getValue() != item->inputsRevision ||
		live->drawMode.getValue() != item->sourceDrawMode ||
		live->representationMode.getValue() !=
		    item->representationMode)
		continue;
	    /* Exact coverage is a producer-certified fact and the priority lane
	     * coalesces to one current overview.  Publish that fact before draining
	     * so this same provider tick can merge the final overview and fulfill
	     * the one-shot autoview; waiting for full job adoption unnecessarily
	     * left a 150k cold scene at its pre-fit camera for tens of seconds. */
	    (void)ged_obol_apply_stream_coverage_bounds(live,
		item->stream.get());
	    const size_t expectedCount = item->stream->getExpectedCount();
	    if (expectedCount)
		live->reserveCompactOccurrenceCapacity(expectedCount);

	    /* A provider budget is a per-pump contract, not a one-batch limit.
	     * Consume successive adaptive batches from this root while there is
	     * cheap work and budget left.  Returning after the first 256 records
	     * coupled realization throughput to paint cadence: an 80 ms OSMesa
	     * frame advanced a large root far more slowly than a 2 ms GL frame even
	     * though both had the same 4 ms owner-thread merge allowance. */
	    for (;;) {
		size_t batch_cap = std::min(item->streamBatchQuantum,
		    GED_OBOL_STREAM_BATCH_QUANTUM_MAX);
		if (max_items) {
		    if (drained_count >= max_items) {
			if (has_more)
			    *has_more = 1;
			return merged;
		    }
		    batch_cap = std::min(batch_cap,
			max_items - drained_count);
		}
		if (max_microseconds &&
		    item->streamMergeMicrosecondsPerItem > 0.0) {
		    const int64_t elapsed = bu_gettime() - started;
		    const uint64_t consumed = elapsed > 0 ?
			static_cast<uint64_t>(elapsed) : 0;
		    const uint64_t remaining =
			consumed < max_microseconds ?
			max_microseconds - consumed : 0;
		    /* Leave headroom for map growth, notifications, and the
		     * controller epilogue.  One occurrence is the irreducible
		     * forward-progress unit if the estimate is sub-quantum. */
		    const double estimatedItems =
			0.70 * static_cast<double>(remaining) /
			item->streamMergeMicrosecondsPerItem;
		    const size_t timedCap = std::max<size_t>(
			1, static_cast<size_t>(estimatedItems));
		    batch_cap = std::min(batch_cap, timedCap);
		}
		std::vector<BObolCompactOccurrence> batch;
		(void)item->stream->drain(batch, batch_cap);
		if (batch.empty())
		    break;
		drained_count += batch.size();
		const int64_t mergeStarted = bu_gettime();
		ged_obol_scene_mutation_batch_scope mutation(itemScene, 1,
		    batch.size());
		const int mergedCount = live->mergeCompactOccurrences(batch);
		/* mergeCompactOccurrences derives a partial registry union.
		 * Restore the immutable whole-target extent before any camera
		 * observer can inspect this atomic scene mutation. */
		(void)ged_obol_apply_stream_coverage_bounds(live,
		    item->stream.get());
		if (mergedCount > 0) {
		    merged = 1;
		    ged_obol_timing_log("stream: merged occurrences",
			mergeStarted, (long)mergedCount);
		}
		const int64_t mergeElapsed = bu_gettime() - mergeStarted;
		if (mergeElapsed >= 0) {
		    const double observed =
			static_cast<double>(mergeElapsed) /
			static_cast<double>(batch.size());
		    if (item->streamMergeMicrosecondsPerItem <= 0.0 ||
			observed > item->streamMergeMicrosecondsPerItem) {
			/* React to worsening costs in one sample. */
			item->streamMergeMicrosecondsPerItem = observed;
		    } else {
			/* Permit gradual recovery when later batches are cheaper. */
			item->streamMergeMicrosecondsPerItem =
			    0.875 * item->streamMergeMicrosecondsPerItem +
			    0.125 * observed;
		    }
		    const double targetMicroseconds =
			max_microseconds ?
			0.70 * static_cast<double>(max_microseconds) : 4000.0;
		    size_t targetItems = static_cast<size_t>(
			targetMicroseconds /
			std::max(0.001,
			    item->streamMergeMicrosecondsPerItem));
		    targetItems = std::max<size_t>(16, std::min(
			targetItems, GED_OBOL_STREAM_BATCH_QUANTUM_MAX));
		    /* Cost can fall sharply on warm records, but bounded growth keeps
		     * one cheap batch from scheduling a nonlinear next-frame jump. */
		    item->streamBatchQuantum = std::min(
			targetItems, std::max<size_t>(
			    16, item->streamBatchQuantum * 2));
		}
		const bool streamHasMore = item->stream->size() > 0;
		if (has_more && streamHasMore)
		    *has_more = 1;
		const int64_t elapsed = bu_gettime() - started;
		if ((max_items && drained_count >= max_items) ||
		    (max_microseconds && elapsed >= 0 &&
		     static_cast<uint64_t>(elapsed) >= max_microseconds)) {
		    if (has_more)
			*has_more = 1;
		    return merged;
		}
		if (!streamHasMore)
		    break;
	    }
	}
    }
    return merged;
}

static bool
ged_obol_deferred_streams_pending(
    const std::shared_ptr<ged_obol_deferred_realization_job> &job)
{
    if (!job)
	return false;
    for (const std::unique_ptr<ged_obol_deferred_realization_item> &item :
	 job->items) {
	if (item && item->stream && item->stream->size() > 0)
	    return true;
    }
    return false;
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
	    bv_frame_revision_get(view) != data->expected_view_revision) {
	    data->pending_autoview = 0;
	    data->pending_autoview_bounds_complete = 0;
	}
    }
    if (data->pending_autoview &&
	ged_obol_deferred_job_coverage_bounds_complete(
	    data->deferred_job)) {
	data->pending_autoview_bounds_complete = 1;
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
    /* Producer completion certifies an immutable exact extent.  The drain
     * above publishes that bound and consumes the stream's single coalesced
     * overview before autoview is applied at the end of this provider tick.
     * Recheck immediately before a completed job is retired below. */
    if (data->pending_autoview &&
	ged_obol_deferred_job_coverage_bounds_complete(data->deferred_job))
	data->pending_autoview_bounds_complete = 1;
    for (std::vector<std::shared_ptr<ged_obol_deferred_realization_job>>::iterator
	    it = data->pending_jobs.begin(); it != data->pending_jobs.end();) {
	    const std::shared_ptr<ged_obol_deferred_realization_job> job = *it;
	    const int job_state = job ?
		job->stateValue() :
		ged_obol_deferred_realization_job::FAILED;
	    if (job_state == ged_obol_deferred_realization_job::PENDING ||
		job_state == ged_obol_deferred_realization_job::RUNNING) {
		has_pending_job = 1;
		++it;
		continue;
	    }
	    /* Completion means the producer is done, not that every completed
	     * occurrence has reached the live scene.  Keep draining at the
	     * provider's item/time budget before adopting the detached result.
	     * Otherwise a warm 50k-leaf manifest can finish between two GUI
	     * ticks and turn the next tick into one unbounded scene mutation. */
	    if (job_state == ged_obol_deferred_realization_job::COMPLETE &&
		ged_obol_deferred_streams_pending(job)) {
		has_pending_job = 1;
		++it;
		continue;
	    }
	    if (job_state == ged_obol_deferred_realization_job::COMPLETE) {
		refined += ged_obol_publish_deferred_realization(data,
			controller, job);
		ged_obol_hold_completed_job(data, job);
	    }
	    it = data->pending_jobs.erase(it);
	}
	if (data->deferred_refine_stage == 1) {
	    const int jobState = data->deferred_job ?
		data->deferred_job->stateValue() :
		ged_obol_deferred_realization_job::FAILED;
	    if (jobState == ged_obol_deferred_realization_job::PENDING ||
		jobState == ged_obol_deferred_realization_job::RUNNING) {
		has_pending_job = 1;
	    } else if (jobState ==
		    ged_obol_deferred_realization_job::COMPLETE &&
		ged_obol_deferred_streams_pending(data->deferred_job)) {
		has_pending_job = 1;
	    } else {
		if (jobState == ged_obol_deferred_realization_job::COMPLETE) {
		    refined += ged_obol_publish_deferred_realization(data,
			controller, data->deferred_job);
		    ged_obol_hold_completed_job(data, data->deferred_job);
		}
		data->deferred_refine_stage = 3;
		data->deferred_paths.clear();
		data->deferred_job.reset();
	    }
	}

    /* Keep pumping only through the short lease handoff window.  This both
     * gives the render-triggered LoD traversal a chance to capture the source
     * and guarantees completed offscreen imports are released without waiting
     * for unrelated future progressive work. */
    if (!data->retired_jobs.empty())
	has_pending_job = 1;

    /* A concurrent worker may have queued more after the pre-retirement
     * drain.  Keep the provider alive for another pump in that case. */
    if (stream_more)
	has_pending_job = 1;

    local_status.changed = refined > 0 ? 1 : 0;
    local_status.hasMore = has_pending_job;
    local_status.inFlight = has_pending_job ? 1 : 0;
    if (has_pending_job && getenv("BOBOL_DRAW_TIMING_VERBOSE")) {
	const int deferred_state = data->deferred_job ?
	    data->deferred_job->stateValue() : -1;
	bu_log("[obol-timing] provider pending: stage=%d deferred=%d "
	       "pending_jobs=%zu stream_more=%d\n",
	       data->deferred_refine_stage, deferred_state,
	       data->pending_jobs.size(), stream_more);
    }

    /* The deferred realization job streams leaf-local boxes/geometry onto the
     * root occurrence registry and atomically adopts the completed index. */
	if (data->pending_autoview &&
	    data->pending_autoview_bounds_complete &&
	    ged_obol_progressive_autoview_apply(data)) {
	    /*
	     * Exact bounds may have been merged by the preceding provider tick,
	     * leaving no new occurrence in this one.  Camera fulfillment is
	     * itself a presentation change and must not wait for an unrelated
	     * later mesh publication.
	     */
	    local_status.changed = 1;
	    controller->syncCameraFromViewContext(view_ctx, TRUE);
	}
	if (!local_status.hasMore) {
	    data->pending_autoview = 0;
	    data->pending_autoview_bounds_complete = 0;
	}
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
	size_t complete_progressive_roots = 0;
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
	    if (!instance_key.empty()) {
		const int proxy_status =
		    ged_obol_publish_deferred_root_proxy(gedp, view_ctx,
			path.c_str(), instance_key.c_str(),
			settings->draw_mode);
		if (proxy_status != 0)
		    changed = 1;
		if (proxy_status == 2)
		    complete_progressive_roots++;
	    }
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
	    result->progressive_data_complete =
		complete_progressive_roots == draw_paths.size() ? 1 : 0;
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

    ged_obol_frontier_visibility_snapshot_apply(gedp);
    ged_obol_frontier_presentation_snapshot_apply(gedp, view_ctx, scene);
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
	/* Database mutation events and explicit stale-source notifications
	 * describe the same retained-source problem.  Try the compact
	 * occurrence diff first so a transform or one changed leaf updates only
	 * that part instead of replacing an aggregate source. */
	if (ged_obol_refresh_matching_compact_parts(scene, view_ctx, targets,
		txn->mode, source_revision) > 0)
	    return 1;
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
	ged_obol_complete_source_roots(source_paths, paths);
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
	    /* Removing a presentation owner is not geometry invalidation.  The
	     * shared repository tracks every source owner and evicts an object's
	     * payload when its last source is released.  Clearing it here made an
	     * erase in one view discard geometry still owned by other views. */
	    if (compact_changed)
		changed = 1;
	}
	    break;
	}
	case GED_DRAW_TXN_CLEAR:
	case GED_DRAW_TXN_TEARDOWN:
	    changed = ged_obol_scene_clear_controller(scene);
	    break;
	case GED_DRAW_TXN_CLEAR_SCOPE:
	    changed = ged_obol_clear_database_sources_in_scope(scene,
		      view_ctx);
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
	case GED_DRAW_TXN_HIGHLIGHTS_CLEAR:
	    changed = ged_obol_scene_highlight_state_set(scene, 0);
	    break;
	case GED_DRAW_TXN_HIGHLIGHT_OCCURRENCE:
	    if (ged_draw_shape_ref_is_null(txn->shape_ref)) {
		changed = ged_obol_scene_highlight_state_set(scene, 0);
		break;
	    }
	    if (!ZERO(txn->value))
		(void)ged_obol_scene_highlight_state_set(scene, 0);
	    changed = ged_draw_shape_ref_set_highlighted(gedp,
		txn->shape_ref, !ZERO(txn->value));
	    break;
	case GED_DRAW_TXN_TRANSPARENCY:
	    changed = ged_obol_apply_transparency_transaction(txn, result,
		      scene);
	    break;
	case GED_DRAW_TXN_MATERIAL_CHANGED: {
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

    BObolSceneController *primary_scene = ged_draw_obol_scene_controller(gedp);
    int changed = primary_scene ?
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
		case GED_DRAW_TXN_HIGHLIGHT:
		case GED_DRAW_TXN_HIGHLIGHTS_CLEAR:
		case GED_DRAW_TXN_HIGHLIGHT_OCCURRENCE:
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
ged_obol_progressive_autoview_apply_exact_proxy_bounds(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    BObolViewController *controller,
    ged_obol_progressive_provider_data *data)
{
    if (!gedp || !view_ctx || !controller || !data ||
	!data->pending_autoview)
	return;

    BObolSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return;

    /* A cached root proxy is useful as an immediate frame only when every
     * displayed source has an exact source bound.  Never frame a partial leaf
     * union: doing so and reframing after discovery is the center/scale jump
     * this transaction contract is designed to prevent. */
    int visible_sources = 0;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid ||
	    !summary.visible)
	    continue;
	visible_sources++;
	if (!summary.sourceBoundsValid || summary.sourceBounds.isEmpty() ||
	    !isfinite(summary.sourceBounds.getMin()[0]) ||
	    !isfinite(summary.sourceBounds.getMin()[1]) ||
	    !isfinite(summary.sourceBounds.getMin()[2]) ||
	    !isfinite(summary.sourceBounds.getMax()[0]) ||
	    !isfinite(summary.sourceBounds.getMax()[1]) ||
	    !isfinite(summary.sourceBounds.getMax()[2]))
	    return;
    }
    if (!visible_sources)
	return;

    data->pending_autoview_bounds_complete = 1;
    if (!ged_obol_progressive_autoview_apply(data))
	return;
    controller->syncCameraFromViewContext(view_ctx, TRUE);
    controller->requestRender("ged-exact-root-proxy-autoview");
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
	!(result && result->progressive_data_complete) &&
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
		data->pending_autoview_bounds_complete = 0;
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
	if (ctx->arm_autoview &&
	    ged_obol_progressive_autoview_arm(data,
		BV_AUTOVIEW_SCALE_DEFAULT))
	    ged_obol_progressive_autoview_apply_exact_proxy_bounds(ctx->gedp,
		view_ctx, controller, data);
	controller->markProgressiveWorkPending();
	return 1;
    };
    ged_bobol_view_controllers_foreach(gedp, update_endpoint,
	&progressive_ctx);
}

extern "C" int
ged_draw_obol_backend_apply_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result)
{
    /* A compact nested draw/erase is already represented by the frontier
     * change stream.  Applying the ordinary transaction as well briefly
     * mutates occurrence state before the authoritative frontier arrives,
     * which can expose a one-frame flicker during resize or refinement. */
    const int frontier_only = result && result->presentation_only &&
	(txn->kind == GED_DRAW_TXN_DRAW || txn->kind == GED_DRAW_TXN_ERASE ||
	 txn->kind == GED_DRAW_TXN_ERASE_PREFIX);
    const int changed = frontier_only ? 0 :
	ged_draw_obol_scene_sync_attached_transaction(gedp, txn, result);
    ged_obol_frontier_visibility_changes_apply(gedp);
    if (txn->kind == GED_DRAW_TXN_CLEAR ||
	txn->kind == GED_DRAW_TXN_TEARDOWN)
	ged_draw_registry_free(gedp);
    ged_obol_progressive_autoview_transaction(gedp, txn, result);
    return changed || frontier_only;
}

static int
ged_obol_observer_ensure(struct ged *gedp, struct ged_drawable *gdp)
{
    if (!gedp || !gdp)
	return 0;
    ged_obol_state *state = ged_obol_state_get(gdp, 1);
    return state ? 1 : 0;
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

    ged_view_lod_policy policy;
    if (!ged_view_lod_policy_get(&policy, view_ctx) ||
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

    /* Keep the renderer-neutral view record synchronized with the retained
     * root it now renders.  This is an identity attachment only: all
     * non-independent views continue to share the same source nodes and
     * geometry repository. */
    (void)ged_draw_source_root_attach_view_contexts(gedp, view_ctx, NULL);

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
    if (ged_view_context_obol_endpoint_get(view_ctx))
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
    if (!ged_view_context_obol_endpoint_set(view_ctx, endpoint, 1)) {
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
