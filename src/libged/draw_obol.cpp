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

#include "brlobol/database_source.h"
#include "brlobol/draw_cache.h"
#include "brlobol/grid.h"
#include "brlobol/init.h"
#include "brlobol/image_source.h"
#include "brlobol/lod_realization.h"
#include "brlobol/lod_service.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/scene_controller.h"
#include "brlobol/scene_group.h"
#include "brlobol/viewport_image.h"
#include "brlobol/vlist_shape.h"
#include "brlobol/view_controller.h"
#include "brlobol/view_store.h"
#include "bg/line_layer.h"
#include "bg/plane.h"
#include "bg/polygon.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/path.h"
#include "bu/str.h"
#include "bu/units.h"
#include "bu/vls.h"
#include "dm.h"
#include "ged.h"
#include "ged/db_index.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "icv.h"
#include "rt/db_fullpath.h"
#include "rt/search.h"
#include "rt/view.h"
#include "vmath.h"

#include "./ged_private.h"

#include <algorithm>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoSeparator.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static const char ged_obol_group_intent_prefix[] = "ged-draw-group:";
static const char ged_obol_database_source_mode_key_marker[] =
    ":ged-draw-mode:";
static thread_local int ged_obol_source_summary_force_adapter = 0;
static thread_local int ged_obol_database_source_publication_depth = 0;
static thread_local void *ged_obol_database_source_publication_view = NULL;
static thread_local int ged_obol_database_source_publication_mode = -1;
static thread_local SoBRLSceneController *ged_obol_database_source_publication_scene = NULL;
static thread_local std::string ged_obol_database_source_publication_group_path;
static thread_local int ged_obol_database_source_publication_appearance = 0;
static thread_local int ged_obol_database_source_publication_line_width = 0;
static thread_local float ged_obol_database_source_publication_transparency = 0.0f;
static thread_local int ged_obol_database_source_publication_color_override = 0;
static thread_local unsigned char ged_obol_database_source_publication_color[3] = {255, 255, 255};

static void ged_obol_transaction_observer(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    void *client_data);
static SbColor ged_obol_color_from_rgb(const unsigned char rgb[3]);
static int ged_obol_observer_ensure(struct ged *gedp,
				    struct ged_drawable *gdp);
static int ged_obol_bind_view_render_root(
    void *view_ctx,
    SoBRLSceneController *shared_scene,
    BRLObolViewController *view_controller);
static int ged_obol_progressive_advance_provider(
    BRLObolViewController *controller,
    void *user_data,
    const BRLObolProgressiveOptions *options,
    BRLObolProgressiveStatus *status);

struct ged_obol_progressive_provider_data {
    struct ged *gedp;
    void *view_ctx;
};

class ged_obol_source_summary_adapter_scope
{
public:
    ged_obol_source_summary_adapter_scope(void)
    {
	ged_obol_source_summary_force_adapter++;
    }

    ~ged_obol_source_summary_adapter_scope(void)
    {
	ged_obol_source_summary_force_adapter--;
    }
};

class ged_obol_scene_mutation_batch_scope
{
public:
    ged_obol_scene_mutation_batch_scope(SoBRLSceneController *scene,
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
    SoBRLSceneController *m_scene;
    bool m_active;
};

struct ged_obol_source_state {
    ged_obol_source_state(void) :
	valid(false),
	viewMatched(false),
	groupValid(false),
	groupDrawMode(BRLOBOL_LOD_DRAW_WIRE),
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
    state.transparency = static_cast<float>(settings->transparency);
    state.colorOverride = settings->color_override ? true : false;
    state.color = ged_obol_color_from_rgb(settings->color);
    if (state.colorOverride) {
	state.materialColorValid = true;
	state.materialColor = state.color;
    }
    if (state.groupValid)
	state.groupTransparency = state.transparency;
}

struct ged_obol_preserved_source {
    std::string path;
    int drawMode;
};

static SbMatrix
ged_obol_sbmatrix_from_mat(const mat_t mat);

struct ged_obol_attached_controller {
    ged_obol_attached_controller(void) :
	view_ctx(NULL),
	scene_controller(NULL),
	view_controller(NULL),
	lod_service(NULL),
	lod_worker_count(0),
	progressive_provider_data(NULL),
	progressive_provider_token(0),
	owned_scene_controller(0),
	owned_view_controller(0),
	owned_lod_service(0),
	full_sync(0),
	use_attached_view_scope(0),
	render_endpoint_only(0),
	render_shared_root_visible(-1)
    {
    }

    void *view_ctx;
    SoBRLSceneController *scene_controller;
    BRLObolViewController *view_controller;
    BRLObolLodService *lod_service;
    size_t lod_worker_count;
    void *progressive_provider_data;
    uint64_t progressive_provider_token;
    int owned_scene_controller;
    int owned_view_controller;
    int owned_lod_service;
    int full_sync;
    int use_attached_view_scope;
    int render_endpoint_only;
    int render_shared_root_visible;
};

static struct ged_drawable *
ged_obol_gdp(struct ged *gedp)
{
    if (!gedp || !gedp->i)
	return NULL;
    return gedp->i->ged_gdp;
}

static std::vector<ged_obol_attached_controller> *
ged_obol_attached_controllers(struct ged_drawable *gdp, int create)
{
    if (!gdp)
	return NULL;
    if (!gdp->gd_obol_attached_controllers && create)
	gdp->gd_obol_attached_controllers =
	    new std::vector<ged_obol_attached_controller>();
    return static_cast<std::vector<ged_obol_attached_controller> *>(
	       gdp->gd_obol_attached_controllers);
}

static void
ged_obol_unregister_progressive_provider(
    const ged_obol_attached_controller &entry)
{
    if (entry.view_controller && entry.progressive_provider_token)
	entry.view_controller->unregisterProgressiveProvider(
	    entry.progressive_provider_token);

    ged_obol_progressive_provider_data *data =
	static_cast<ged_obol_progressive_provider_data *>(
	    entry.progressive_provider_data);
    delete data;
}

static void
ged_obol_stop_attached_controller_lod(const ged_obol_attached_controller &entry)
{
    if (entry.view_controller && entry.view_controller->getLodService() ==
	entry.lod_service) {
	entry.view_controller->setLodAutoSubmit(FALSE);
	entry.view_controller->setLodService(NULL);
    }

    if (entry.owned_lod_service && entry.lod_service) {
	entry.lod_service->stop();
    }
}

static void
ged_obol_delete_owned_attached_controller(
    const ged_obol_attached_controller &entry)
{
    ged_obol_unregister_progressive_provider(entry);
    ged_obol_stop_attached_controller_lod(entry);

    if (entry.owned_lod_service && entry.lod_service) {
	delete entry.lod_service;
    }

    if (entry.owned_view_controller && entry.view_controller) {
	delete entry.view_controller;
	return;
    }

    if (entry.owned_scene_controller && entry.scene_controller)
	delete entry.scene_controller;
}

static void
ged_obol_primary_clear(struct ged_drawable *gdp)
{
    if (!gdp)
	return;
    gdp->gd_obol_scene_controller = NULL;
    gdp->gd_obol_controller = NULL;
    gdp->gd_obol_scene_controller_owned = 0;
    gdp->gd_obol_controller_owned = 0;
    gdp->gd_obol_scene_controller_full_sync = 0;
}

static void
ged_obol_primary_set(struct ged_drawable *gdp,
		     const ged_obol_attached_controller *entry)
{
    if (!gdp || !entry) {
	ged_obol_primary_clear(gdp);
	return;
    }

    gdp->gd_obol_scene_controller = entry->scene_controller;
    gdp->gd_obol_controller = entry->view_controller;
    gdp->gd_obol_scene_controller_owned = entry->owned_scene_controller ? 1 : 0;
    gdp->gd_obol_controller_owned = entry->owned_view_controller ? 1 : 0;
    gdp->gd_obol_scene_controller_full_sync = entry->full_sync ? 1 : 0;
}

static void
ged_obol_primary_refresh_from_registry(struct ged_drawable *gdp)
{
    std::vector<ged_obol_attached_controller> *entries =
	ged_obol_attached_controllers(gdp, 0);
    if (!entries || entries->empty()) {
	ged_obol_primary_clear(gdp);
	return;
    }

    for (const ged_obol_attached_controller &entry : *entries) {
	if (!entry.render_endpoint_only) {
	    ged_obol_primary_set(gdp, &entry);
	    return;
	}
    }

    ged_obol_primary_set(gdp, &entries->front());
}

static void
ged_obol_attached_controllers_free(struct ged_drawable *gdp)
{
    std::vector<ged_obol_attached_controller> *entries =
	ged_obol_attached_controllers(gdp, 0);
    if (!entries)
	return;

    for (const ged_obol_attached_controller &entry : *entries)
	ged_obol_stop_attached_controller_lod(entry);
    for (const ged_obol_attached_controller &entry : *entries)
	ged_obol_delete_owned_attached_controller(entry);
    delete entries;
    gdp->gd_obol_attached_controllers = NULL;
    ged_obol_primary_clear(gdp);
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
	if (!paths.empty())
	    return paths;
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
	return BRLOBOL_LOD_DRAW_SHADED_BOTS;
    if (mode == GED_DRAW_MODE_SHADED)
	return BRLOBOL_LOD_DRAW_SHADED;
    if (mode == GED_DRAW_MODE_HIDDEN_LINE)
	return BRLOBOL_LOD_DRAW_HIDDEN_LINE;
    if (mode == GED_DRAW_MODE_EVAL_POINTS)
	return BRLOBOL_LOD_DRAW_POINTS;
    return BRLOBOL_LOD_DRAW_WIRE;
}

static int
ged_obol_lod_draw_mode_to_ged(int mode)
{
    if (mode == BRLOBOL_LOD_DRAW_SHADED_BOTS)
	return GED_DRAW_MODE_SHADED_BOTS;
    if (mode == BRLOBOL_LOD_DRAW_SHADED)
	return GED_DRAW_MODE_SHADED;
    if (mode == BRLOBOL_LOD_DRAW_HIDDEN_LINE)
	return GED_DRAW_MODE_HIDDEN_LINE;
    if (mode == BRLOBOL_LOD_DRAW_POINTS)
	return GED_DRAW_MODE_EVAL_POINTS;
    return GED_DRAW_MODE_WIRE;
}

static int
ged_obol_lod_draw_mode_from_database_source(const SoBRLDatabaseSource *source)
{
    if (!source)
	return BRLOBOL_LOD_DRAW_WIRE;

    switch (source->representationMode.getValue()) {
	case SoBRLDatabaseSource::REPRESENTATION_SHADED_BOTS:
	    return BRLOBOL_LOD_DRAW_SHADED_BOTS;
	case SoBRLDatabaseSource::REPRESENTATION_SHADED:
	    return BRLOBOL_LOD_DRAW_SHADED;
	case SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE:
	    return BRLOBOL_LOD_DRAW_HIDDEN_LINE;
	case SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS:
	    return BRLOBOL_LOD_DRAW_POINTS;
	case SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE:
	case SoBRLDatabaseSource::REPRESENTATION_WIRE:
	    return BRLOBOL_LOD_DRAW_WIRE;
	default:
	    break;
    }

    return source->drawMode.getValue() == SoBRLDatabaseSource::SHADED ?
	   BRLOBOL_LOD_DRAW_SHADED : BRLOBOL_LOD_DRAW_WIRE;
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
ged_obol_drawn_path_mode(struct ged *gedp, void *view_ctx,
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
ged_obol_view_scope_is_independent(void *view_ctx)
{
    return view_ctx && rt_view_context_is_independent(view_ctx);
}

static std::string
ged_obol_view_scope_name(void *view_ctx)
{
    if (!view_ctx)
	return "shared";

    const char *name = ged_view_context_name_get(view_ctx);
    if (name && name[0])
	return std::string(name);

    char fallback[64] = {0};
    snprintf(fallback, sizeof(fallback), "%p", view_ctx);
    return std::string(fallback);
}

static std::string
ged_obol_database_source_instance_key(void *view_ctx, const char *path)
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
    void *view_ctx,
    const char *path,
    int ged_draw_mode)
{
    std::string key = ged_obol_database_source_instance_key(view_ctx, path);
    if (key.empty() || ged_draw_mode < 0 ||
	ged_draw_mode == GED_DRAW_MODE_WIRE ||
	ged_obol_view_scope_is_independent(view_ctx))
	return key;

    char mode_buf[64] = {0};
    snprintf(mode_buf, sizeof(mode_buf), "%s%d",
	     ged_obol_database_source_mode_key_marker, ged_draw_mode);
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
    const BRLObolDatabaseSourceSummary &summary)
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
    const BRLObolDatabaseSourceSummary &summary,
    int ged_draw_mode)
{
    return ged_draw_mode < 0 ||
	   ged_obol_database_source_summary_ged_mode(summary) == ged_draw_mode;
}

static std::string
ged_obol_database_source_instance_prefix(void *view_ctx)
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
    const BRLObolDatabaseSourceSummary &summary,
    void *view_ctx)
{
    if (!ged_obol_view_scope_is_independent(view_ctx)) {
	if (summary.instanceKey.getLength() == 0)
	    return 1;
	const std::string base_instance_key =
	    ged_obol_database_source_base_instance_key(
		summary.instanceKey.getString());
	return ged_obol_path_equal(base_instance_key.c_str(),
				   summary.path.getString());
    }

    const std::string prefix =
	ged_obol_database_source_instance_prefix(view_ctx);
    const char *key = summary.instanceKey.getString();
    return key && strncmp(key, prefix.c_str(), prefix.size()) == 0;
}

static std::string
ged_obol_database_source_owner_group_path_from_summary(
    const BRLObolDatabaseSourceSummary &summary)
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
    state.highlighted = record->highlighted ? true : false;
    state.lineWidth = record->line_width;
    state.transparency = static_cast<float>(record->transparency);
    ged_obol_source_state_add_group(state, gedp, record->group);

    {
	ged_obol_source_summary_adapter_scope adapter_summary_scope;
	struct ged_draw_scene_display_summary display_summary;
	if (ged_draw_shape_ref_display_summary(gedp, record->ref,
					       &display_summary) &&
	    display_summary.valid) {
	    state.visible = display_summary.visible ? true : false;
	    state.selected = display_summary.selected ? true : false;
	    state.highlighted = display_summary.highlighted ? true : false;
	    state.lineStyle = display_summary.line_style;
	    state.lineWidth = display_summary.line_width;
	    state.transparency = static_cast<float>(display_summary.transparency);
	    state.materialColorValid =
		display_summary.material_valid ? true : false;
	    if (state.materialColorValid) {
		state.materialColor = SbColor(
					  static_cast<float>(display_summary.material_color[0]) /
					  255.0f,
					  static_cast<float>(display_summary.material_color[1]) /
					  255.0f,
					  static_cast<float>(display_summary.material_color[2]) /
					  255.0f);
	    }
	}

	struct ged_draw_shape_material_summary material_summary;
	if (ged_draw_shape_ref_material_summary(gedp, record->ref,
						&material_summary) && material_summary.valid)
	    state.materialRevision =
		ged_obol_fold_revision(material_summary.material_revision);

    }

    if ((!state.materialColorValid || state.materialRevision == 0) &&
	!state.colorOverride && gedp && gedp->dbip &&
	record->fullpath && record->fullpath->fp_len > 0) {
	SbColor db_material_color;
	if (brlobol_database_source_fullpath_material_color(
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
	int ged_draw_mode)
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
    state.drawMode = ged_obol_database_draw_mode_from_ged(ged_draw_mode);
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
	if (gedp && brlobol_database_source_path_material_color(gedp->dbip, path,
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
	state.groupDrawMode = ged_obol_lod_draw_mode_from_ged(ged_draw_mode);
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
    if (brlobol_database_source_path_material_color(dbip, path,
	    db_material_color)) {
	state.materialColorValid = true;
	state.materialColor = db_material_color;
    }
}

struct ged_obol_find_source_state_context {
    struct ged *gedp;
    void *viewCtx;
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
			   void *view_ctx,
			   const char *path,
			   int ged_draw_mode)
{
    if (!gedp || !path || !path[0])
	return ged_obol_source_state();

    ged_obol_source_state direct_state;
    if (ged_obol_source_state_from_database_source(direct_state, gedp, path,
	    ged_draw_mode))
	return direct_state;

    if (ged_draw_mode >= 0 &&
	ged_draw_path_state(gedp, view_ctx, path, ged_draw_mode) <= 0)
	return ged_obol_source_state();

    ged_obol_find_source_state_context ctx;
    ctx.gedp = gedp;
    ctx.viewCtx = view_ctx;
    ctx.path = path;
    ctx.mode = ged_draw_mode;
    ged_draw_foreach_shape_record(gedp, ged_obol_find_source_state_cb, &ctx);
    return ctx.state;
}

static int
ged_obol_intent_is_ged_draw_group(const SbString &intent)
{
    const char *value = intent.getString();
    if (!value)
	return 0;
    return strncmp(value, ged_obol_group_intent_prefix,
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
    if (strncmp(path, ged_obol_group_intent_prefix, prefix_len) == 0)
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
    if (strncmp(path, ged_obol_group_intent_prefix, prefix_len) == 0)
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
ged_obol_sync_group_state(SoBRLSceneController *scene,
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
				  BRLOBOL_LOD_DRAW_WIRE,
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
	BRLObolDatabaseSourceSummary summary;
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
ged_obol_prune_empty_groups(SoBRLSceneController *scene)
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
	    BRLObolSceneTreeSummary tree_summary;
	    if (!scene->getSceneTreeSummary(i, tree_summary) ||
		!tree_summary.valid ||
		tree_summary.nodeKind !=
		BRLObolSceneTreeSummary::NODE_GROUP ||
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
    if (strncmp(path, prefix, prefix_len) != 0)
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
	    strncmp(cursor, name, component_len) == 0)
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
    void *view_ctx)
{
    if (!group_summary)
	return 0;
    if (!view_ctx)
	return !group_summary->in_view_scope;
    if (rt_view_context_is_independent(view_ctx))
	return group_summary->in_view_scope &&
	       group_summary->view_ctx == view_ctx;
    if (!group_summary->in_view_scope)
	return 1;
    return group_summary->view_ctx == view_ctx;
}

struct ged_obol_drawn_source_path_ctx {
    struct ged *gedp;
    void *viewCtx;
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
    const BRLObolDatabaseSourceSummary &summary)
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
    const BRLObolDatabaseSourceSummary &summary)
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(ctx->gedp);
    if (!scene)
	return;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;
	ged_obol_drawn_source_summary_append(ctx, summary);
    }
}

static std::vector<ged_obol_drawn_source_path_mode>
ged_obol_drawn_source_path_modes(struct ged *gedp,
				 void *view_ctx,
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
			    void *view_ctx,
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
ged_obol_database_source_mark_published_current(SoBRLSceneController *scene,
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
ged_obol_view_lod_policy_state_for_source(struct ged *gedp, void *view_ctx)
{
    ged_obol_view_lod_policy_state state;
    if (!gedp)
	return state;

    void *policy_view = view_ctx ? view_ctx : ged_view_active_ctx(gedp);
    if (!policy_view)
	return state;

    ged_draw_view_lod_policy policy;
    if (!ged_draw_view_context_lod_policy_get(&policy, policy_view))
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
	static_cast<float>(ged_view_context_scale_get(policy_view));
    state.lodScale = static_cast<float>(policy.scale);
    state.viewWidth = rt_view_context_width_get(policy_view);
    state.viewHeight = rt_view_context_height_get(policy_view);
    state.curveScale = static_cast<float>(policy.curve_scale);
    state.pointScale = static_cast<float>(policy.point_scale);
    return state;
}

static int
ged_obol_apply_view_lod_policy(struct ged *gedp,
			       void *view_ctx,
			       SoBRLSceneController *scene,
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
		      void *view_ctx,
		      struct db_i *dbip,
		      const char *path,
		      int ged_draw_mode,
		      uint32_t source_revision,
		      SoBRLSceneController *scene,
		      int use_retained_state = 1,
		      int preserve_existing_revision = 0,
		      const struct ged_draw_appearance_settings *appearance_settings = NULL,
		      const ged_obol_publish_placement_state *placement_state = NULL)
{
    if (!dbip || !path || !path[0] || !scene)
	return 0;
    if (!ged_obol_path_can_publish_database_source(dbip, path))
	return 0;

    ged_obol_source_state source_state = use_retained_state ?
					 ged_obol_find_source_state(gedp, view_ctx, path, ged_draw_mode) :
					 ged_obol_source_state();
    const bool retained_state_valid = source_state.valid;
    const uint32_t retained_source_revision = source_state.sourceRevision;
    const uint32_t retained_inputs_revision = source_state.inputsRevision;
    if (source_state.valid) {
	if (source_state.sourceRevision != 0)
	    source_revision = source_state.sourceRevision;
    }

    const std::string instance_key =
	ged_obol_database_source_instance_key_for_mode(view_ctx, path,
	    ged_draw_mode);
    const std::string representation_key = instance_key;
    const int next_draw_mode = ged_obol_database_draw_mode_from_ged(
				   ged_draw_mode);
    if (ged_draw_mode >= 0) {
	const std::string base_instance_key =
	    ged_obol_database_source_instance_key(view_ctx, path);
	if (!base_instance_key.empty() &&
	    base_instance_key != instance_key) {
	    SoBRLDatabaseSource *base_source =
		scene->findDatabaseSourceInstance(base_instance_key.c_str());
	    BRLObolDatabaseSourceSummary base_summary;
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
			    ged_draw_mode);
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
    BRLObolDatabaseSourceSummary existing_summary;
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
	    existing_summary.representationMode == ged_draw_mode &&
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
	    if (has_external_primary)
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
    const bool had_group_state = source_state.groupValid;
    if (ged_draw_mode >= 0) {
	source_state.valid = true;
	if (!ged_obol_database_source_publication_appearance) {
	    source_state.groupValid = true;
	}
	if (source_state.groupValid && source_state.groupPath.empty())
	    source_state.groupPath = ged_obol_skip_leading_slash(path);
	source_state.groupDrawMode =
	    ged_obol_lod_draw_mode_from_ged(ged_draw_mode);
	if (source_state.groupValid && !had_group_state) {
	    source_state.groupVisible = true;
	    source_state.groupTransparency = source_state.transparency;
	}
	source_state.groupOverlay = false;
    }
    if (ged_obol_database_source_publication_appearance) {
	source_state.lineWidth =
	    ged_obol_database_source_publication_line_width;
	source_state.transparency =
	    ged_obol_database_source_publication_transparency;
	source_state.colorOverride =
	    ged_obol_database_source_publication_color_override ? true : false;
	source_state.color = ged_obol_color_from_rgb(
				 ged_obol_database_source_publication_color);
	if (source_state.colorOverride) {
	    source_state.materialColorValid = true;
	    source_state.materialColor = source_state.color;
	}
	if (source_state.groupValid)
	    source_state.groupTransparency = source_state.transparency;
    }
    ged_obol_source_state_apply_appearance(source_state, appearance_settings);
    if (ged_obol_database_source_publication_depth > 0 &&
	!ged_obol_database_source_publication_group_path.empty() &&
	(!source_state.groupValid || !retained_state_valid)) {
	source_state.groupValid = true;
	source_state.groupPath = ged_obol_database_source_publication_group_path;
	source_state.groupDrawMode =
	    ged_obol_lod_draw_mode_from_ged(ged_draw_mode);
	source_state.groupVisible = true;
	source_state.groupOverlay = false;
	source_state.groupTransparency = source_state.transparency;
    }

    ged_obol_source_state_resolve_database_material(source_state, dbip, path);

    ged_obol_view_lod_policy_state policy_state =
	ged_obol_view_lod_policy_state_for_source(gedp, view_ctx);

    BRLObolDatabaseSourcePublishState publish_state;
    publish_state.sourceInstanceKey = instance_key.c_str();
    publish_state.sourcePath = path;
    publish_state.sourceRepresentationKey = representation_key.c_str();
    publish_state.targetGroupPath =
	(source_state.groupValid && !source_state.groupPath.empty()) ?
	source_state.groupPath.c_str() : NULL;
    publish_state.database = dbip;
    publish_state.drawMode = next_draw_mode;
    publish_state.representationMode = ged_draw_mode;
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
	if (ged_draw_mode == GED_DRAW_MODE_EVAL_POINTS)
	    publish_state.roleFlags |=
		SoBRLDatabaseSource::REALIZATION_ROLE_MESH;
    }
    if (policy_state.valid) {
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
	BRLObolDatabaseSourceSummary summary;
	if (source && source->getSummary(summary) && summary.valid &&
	    (summary.realizedShapeCount > 0 || summary.realizedMeshCount > 0) &&
	    ged_obol_database_source_mark_published_current(scene, source))
	    changed = 1;
    }
    return changed;
}

static int
ged_obol_replace_path_and_realize(struct ged *gedp,
				  void *view_ctx,
				  struct db_i *dbip,
				  const char *path,
				  int ged_draw_mode,
				  uint32_t source_revision,
				  SoBRLSceneController *scene,
				  int use_retained_state = 1,
				  int preserve_existing_revision = 0,
				  const struct ged_draw_appearance_settings *appearance_settings = NULL)
{
    const int changed = ged_obol_replace_path(gedp, view_ctx, dbip, path,
			ged_draw_mode, source_revision, scene,
			use_retained_state, preserve_existing_revision,
			appearance_settings);
    if (changed < 0)
	return changed;

    if ((ged_draw_mode == GED_DRAW_MODE_EVAL_WIRE ||
	 ged_draw_mode == GED_DRAW_MODE_EVAL_POINTS) &&
	ged_draw_obol_database_source_realize_for_path(gedp, path))
	return 1;

    return changed;
}

static int
ged_obol_replace_paths(struct db_i *dbip,
		       const std::vector<std::string> &paths,
		       int ged_draw_mode,
		       uint32_t source_revision,
		       struct ged *gedp,
		       void *view_ctx,
		       SoBRLSceneController *scene,
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
					      ged_draw_mode, source_revision, scene,
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
			    void *view_ctx,
			    SoBRLSceneController *scene,
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
    const BRLObolDatabaseSourceSummary &summary);

static void
ged_obol_collect_database_sources_matching(
    std::vector<std::string> &paths,
    SoBRLSceneController *scene,
    void *view_ctx,
    const char *target_path,
    size_t component_first_idx,
    int allow_path_prefix,
    int ged_draw_mode)
{
    if (!scene || !target_path || !target_path[0])
	return;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (!source)
	    continue;

	BRLObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary,
		    ged_draw_mode))
	    continue;

	const char *source_path = source->path.getValue().getString();
	if ((allow_path_prefix &&
	     ged_obol_path_has_prefix(source_path, target_path)) ||
	    ged_obol_path_has_component_name(source_path, target_path,
					     component_first_idx))
	    ged_obol_append_unique_path(paths, source_path);
    }
}

static void
ged_obol_collect_database_source_instance_keys_matching(
    std::vector<std::string> &instance_keys,
    SoBRLSceneController *scene,
    void *view_ctx,
    const char *target_path,
    size_t component_first_idx,
    int allow_path_prefix,
    int ged_draw_mode)
{
    if (!scene || !target_path || !target_path[0])
	return;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary,
		    ged_draw_mode))
	    continue;

	const char *source_path = summary.path.getString();
	if ((allow_path_prefix &&
	     ged_obol_path_has_prefix(source_path, target_path)) ||
	    ged_obol_path_has_component_name(source_path, target_path,
					     component_first_idx))
	    ged_obol_append_database_source_instance_key(instance_keys,
		    summary);
    }
}

static std::vector<std::string>
ged_obol_matching_database_source_paths(
    SoBRLSceneController *scene,
    void *view_ctx,
    const std::vector<std::string> &targets,
    size_t component_first_idx,
    int allow_path_prefix,
    int ged_draw_mode)
{
    std::vector<std::string> paths;
    if (!scene)
	return paths;

    for (const std::string &target : targets)
	ged_obol_collect_database_sources_matching(paths, scene, view_ctx,
		target.c_str(), component_first_idx, allow_path_prefix,
		ged_draw_mode);
    return paths;
}

static std::vector<std::string>
ged_obol_primary_matching_database_source_paths(
    struct ged *gedp,
    void *view_ctx,
    const std::vector<std::string> &targets,
    int ged_draw_mode)
{
    std::vector<std::string> source_paths;
    if (!gedp || targets.empty())
	return source_paths;

    SoBRLSceneController *primary_scene =
	ged_draw_obol_scene_controller(gedp);
    if (!primary_scene)
	return source_paths;

    for (const std::string &target : targets) {
	ged_obol_collect_database_sources_matching(source_paths,
		primary_scene, view_ctx, target.c_str(), 0, 1, ged_draw_mode);
    }
    ged_obol_remove_shadowed_source_paths(source_paths);
    return source_paths;
}

static std::vector<std::string>
ged_obol_matching_database_source_instance_keys(
    SoBRLSceneController *scene,
    void *view_ctx,
    const std::vector<std::string> &targets,
    size_t component_first_idx,
    int allow_path_prefix,
    int ged_draw_mode)
{
    std::vector<std::string> instance_keys;
    if (!scene)
	return instance_keys;

    for (const std::string &target : targets)
	ged_obol_collect_database_source_instance_keys_matching(
	    instance_keys, scene, view_ctx, target.c_str(),
	    component_first_idx, allow_path_prefix, ged_draw_mode);
    return instance_keys;
}

static int
ged_obol_replace_matching_database_sources(
    struct ged *gedp,
    void *view_ctx,
    const std::vector<std::string> &targets,
    size_t component_first_idx,
    int allow_path_prefix,
    uint32_t source_revision,
    SoBRLSceneController *scene,
    int ged_draw_mode)
{
    if (!gedp || !gedp->dbip || !scene || targets.empty())
	return 0;

    std::vector<std::string> source_paths =
	ged_obol_matching_database_source_paths(scene, view_ctx, targets,
	    component_first_idx, allow_path_prefix, ged_draw_mode);
    ged_obol_remove_shadowed_source_paths(source_paths);
    if (source_paths.empty())
	return 0;

    ged_obol_scene_mutation_batch_scope batch(scene, source_paths.size(),
	    source_paths.size());
    int changed = 0;
    std::vector<ged_obol_drawn_source_path_mode> source_path_modes =
	ged_obol_drawn_source_path_modes(gedp, view_ctx, ged_draw_mode,
					 &source_paths);
    if (!source_path_modes.empty()) {
	if (ged_obol_replace_path_modes(gedp->dbip, source_path_modes,
					source_revision, gedp, view_ctx, scene, 0))
	    changed = 1;
    } else {
	for (const std::string &source_path : source_paths) {
	    const int selected_mode = ged_draw_mode >= 0 ? ged_draw_mode :
				      ged_obol_drawn_path_mode(gedp, view_ctx,
					  source_path.c_str());
	    if (ged_obol_replace_path_and_realize(gedp, view_ctx, gedp->dbip,
						  source_path.c_str(), selected_mode,
						  source_revision, scene, 0) > 0)
		changed = 1;
	}
    }

    if (changed)
	scene->realizePending();
    return 1;
}

static int
ged_obol_mark_matching_database_sources_stale(
    void *view_ctx,
    const std::vector<std::string> &targets,
    size_t component_first_idx,
    int allow_path_prefix,
    uint32_t stale_reason,
    SoBRLSceneController *scene,
    int ged_draw_mode)
{
    if (!scene || targets.empty())
	return 0;

    std::vector<std::string> instance_keys =
	ged_obol_matching_database_source_instance_keys(scene, view_ctx,
	    targets, component_first_idx, allow_path_prefix,
	    ged_draw_mode);
    if (instance_keys.empty())
	return 0;

    for (const std::string &instance_key : instance_keys)
	(void)scene->markDatabaseSourceInstanceStale(instance_key.c_str(),
		stale_reason);
    return 1;
}

static int
ged_obol_set_database_source_visible(SoBRLSceneController *scene,
				     const char *source_instance_key,
				     int visible)
{
    if (!scene || !source_instance_key || !source_instance_key[0])
	return 0;

    SoBRLDatabaseSource *source =
	scene->findDatabaseSourceInstance(source_instance_key);
    if (!source)
	return 0;

    BRLObolDatabaseSourceSummary summary;
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
ged_obol_set_group_visible(SoBRLSceneController *scene,
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
ged_obol_set_database_source_highlighted(SoBRLSceneController *scene,
	const char *source_instance_key,
	int highlighted)
{
    if (!scene || !source_instance_key || !source_instance_key[0])
	return 0;

    SoBRLDatabaseSource *source =
	scene->findDatabaseSourceInstance(source_instance_key);
    if (!source)
	return 0;

    BRLObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    (void)scene->setDatabaseSourceInstanceState(source_instance_key,
	    TRUE,
	    summary.sourceRevision,
	    summary.inputsRevision,
	    summary.visible,
	    summary.selected,
	    highlighted ? TRUE : FALSE,
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
ged_obol_set_group_highlighted(SoBRLSceneController *scene,
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
    const BRLObolDatabaseSourceSummary &summary);

extern "C" int
ged_draw_obol_highlight_state_set(struct ged *gedp, int highlighted)
{
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
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
	BRLObolSceneTreeSummary tree_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !tree_summary.valid ||
	    tree_summary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
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
    SoBRLSceneController *scene)
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
	    BRLObolSceneTreeSummary tree_summary;
	    if (!scene->getSceneTreeSummary(i, tree_summary) ||
		!tree_summary.valid ||
		tree_summary.nodeKind !=
		BRLObolSceneTreeSummary::NODE_GROUP ||
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
    SoBRLSceneController *scene)
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
	    BRLObolSceneTreeSummary tree_summary;
	    if (!scene->getSceneTreeSummary(i, tree_summary) ||
		!tree_summary.valid ||
		tree_summary.nodeKind !=
		BRLObolSceneTreeSummary::NODE_GROUP ||
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
		      void *view_ctx,
		      SoBRLSceneController *scene,
		      int ged_draw_mode = -1)
{
    if (paths.empty() || !scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx) ||
	    !ged_obol_database_source_summary_matches_mode(summary,
		    ged_draw_mode))
	    continue;

	const char *source_path = summary.path.getString();
	for (const std::string &path : paths) {
	    if (ged_obol_path_equal(source_path, path.c_str())) {
		ged_obol_append_database_source_instance_key(instance_keys,
			summary);
		break;
	    }
	}
    }

    if (instance_keys.empty()) {
	for (const std::string &path : paths) {
	    const std::string instance_key =
		ged_obol_database_source_instance_key_for_mode(view_ctx,
		    path.c_str(), ged_draw_mode);
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
    const BRLObolDatabaseSourceSummary &summary)
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
ged_obol_remove_instance_keys(const std::vector<std::string> &instance_keys,
			      SoBRLSceneController *scene)
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
ged_obol_clear_database_sources_in_scope(SoBRLSceneController *scene,
	void *view_ctx)
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
	BRLObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    !summary.valid ||
	    !ged_obol_database_source_instance_in_scope(summary, view_ctx))
	    continue;

	ged_obol_append_database_source_instance_key(instance_keys, summary);
    }

    return ged_obol_remove_instance_keys(instance_keys, scene);
}

SoBRLSceneController *
ged_draw_obol_scene_controller(struct ged *gedp)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp)
	return NULL;
    return static_cast<SoBRLSceneController *>(
	       gdp->gd_obol_scene_controller);
}

BRLObolViewController *
ged_draw_obol_controller(struct ged *gedp)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp)
	return NULL;
    return static_cast<BRLObolViewController *>(gdp->gd_obol_controller);
}

static ged_obol_attached_controller *
ged_obol_attached_controller_for_view(struct ged *gedp, void *view_ctx)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    std::vector<ged_obol_attached_controller> *entries =
	ged_obol_attached_controllers(gdp, 0);
    if (!entries || !view_ctx)
	return NULL;

    for (ged_obol_attached_controller &entry : *entries) {
	if (entry.view_ctx == view_ctx && entry.view_controller)
	    return &entry;
    }

    return NULL;
}

static void
ged_obol_advance_lod_policy_revision(BRLObolViewController *controller)
{
    if (!controller)
	return;

    const uint64_t current_revision = controller->getLodPolicyRevision();
    controller->setLodPolicyRevision(
	current_revision == UINT64_MAX ? 1 : current_revision + 1);
}

extern "C" int
ged_draw_obol_view_lod_policy_changed(struct ged *gedp, void *view_ctx)
{
    ged_obol_attached_controller *entry =
	ged_obol_attached_controller_for_view(gedp, view_ctx);
    BRLObolViewController *view_controller = entry ?
	entry->view_controller :
	ged_draw_obol_controller(gedp);
    if (!view_controller)
	return 0;

    view_controller->clearViewLodState();
    ged_obol_advance_lod_policy_revision(view_controller);

    SoBRLSceneController *scene = entry ? entry->scene_controller : NULL;
    if (!scene)
	scene = view_controller->getSceneController();
    if (!scene)
	scene = ged_draw_obol_scene_controller(gedp);

    int changed = 1;
    if (scene) {
	const int independent_view =
	    ged_obol_view_scope_is_independent(view_ctx);
	const ged_obol_view_lod_policy_state policy_state =
	    ged_obol_view_lod_policy_state_for_source(gedp, view_ctx);
	const int source_count = scene->getDatabaseSourceCount();
	for (int i = 0; i < source_count; i++) {
	    BRLObolDatabaseSourceSummary summary;
	    if (!scene->getDatabaseSourceSummary(i, summary) ||
		!summary.valid || summary.instanceKey.getLength() == 0)
		continue;
	    if (independent_view &&
		!ged_obol_database_source_instance_in_scope(summary, view_ctx))
		continue;
	    if (ged_obol_apply_view_lod_policy(gedp, view_ctx, scene,
					       summary.instanceKey.getString()) > 0)
		changed = 1;
	    if (policy_state.valid && !policy_state.meshEnabled &&
		scene->markDatabaseSourceInstanceStale(
		    summary.instanceKey.getString(),
		    SoBRLDatabaseSource::STALE_VIEW) > 0)
		changed = 1;
	}
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
    return cpus > 0 ? cpus : 1;
}

static void
ged_obol_lod_status_fill(const ged_obol_attached_controller *entry,
			 struct ged_draw_obol_lod_service_status *status)
{
    if (!status)
	return;

    memset(status, 0, sizeof(*status));
    if (!entry || !entry->view_controller)
	return;

    BRLObolLodService *service = entry->lod_service;
    status->attached = 1;
    status->auto_submit =
	entry->view_controller->isLodAutoSubmitEnabled() ? 1 : 0;
    status->worker_count = entry->lod_worker_count;
    status->last_visited_mesh_count =
	entry->view_controller->getLastLodVisitedMeshCount();
    status->last_submitted_task_count =
	entry->view_controller->getLastLodSubmittedTaskCount();
    status->last_skipped_mesh_count =
	entry->view_controller->getLastLodSkippedMeshCount();
    status->last_result_count =
	entry->view_controller->getLastLodResultCount();
    status->last_matched_result_count =
	entry->view_controller->getLastLodMatchedResultCount();
    status->last_applied_result_count =
	entry->view_controller->getLastLodAppliedResultCount();
    status->last_rejected_result_count =
	entry->view_controller->getLastLodRejectedResultCount();
    status->last_unmatched_result_count =
	entry->view_controller->getLastLodUnmatchedResultCount();
    status->active_mesh_payloads =
	entry->view_controller->getActiveLodMeshPayloadCount();
    status->active_aabb_proxy_payloads =
	entry->view_controller->getActiveLodProxyPayloadCount(
	    BRLOBOL_LOD_PROXY_AABB);
    status->active_obb_proxy_payloads =
	entry->view_controller->getActiveLodProxyPayloadCount(
	    BRLOBOL_LOD_PROXY_OBB);
    bu_strlcpy(status->last_diagnostics,
	       entry->view_controller->getLastLodDiagnostics().getString(),
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
				void *view_ctx,
				size_t worker_count)
{
    ged_obol_attached_controller *entry =
	ged_obol_attached_controller_for_view(gedp, view_ctx);
    if (!entry || !entry->view_controller)
	return 0;

    worker_count = ged_obol_lod_default_worker_count(worker_count);
    if (!entry->lod_service) {
	entry->lod_service = new BRLObolLodService;
	entry->owned_lod_service = 1;
    }

    if (!entry->lod_service)
	return 0;

    if (entry->lod_service->isRunning() &&
	entry->lod_worker_count != worker_count)
	entry->lod_service->stop();

    if (!entry->lod_service->isRunning() &&
	!entry->lod_service->start(worker_count, TRUE))
	return 0;

    entry->lod_worker_count = worker_count;
    entry->view_controller->setLodService(entry->lod_service);
    entry->view_controller->setLodAutoSubmit(TRUE);
    entry->view_controller->requestRender("lod-service-start");
    return 1;
}

extern "C" int
ged_draw_obol_lod_service_stop(struct ged *gedp,
			       void *view_ctx)
{
    ged_obol_attached_controller *entry =
	ged_obol_attached_controller_for_view(gedp, view_ctx);
    if (!entry || !entry->view_controller)
	return 0;

    entry->view_controller->setLodAutoSubmit(FALSE);
    if (entry->view_controller->getLodService() == entry->lod_service)
	entry->view_controller->setLodService(NULL);
    if (entry->lod_service)
	entry->lod_service->stop();
    entry->lod_worker_count = 0;
    entry->view_controller->requestRender("lod-service-stop");
    return 1;
}

extern "C" int
ged_draw_obol_lod_service_poll(struct ged *gedp,
			       void *view_ctx,
			       size_t max_results,
			       struct ged_draw_obol_lod_service_status *status)
{
    ged_obol_attached_controller *entry =
	ged_obol_attached_controller_for_view(gedp, view_ctx);
    if (!entry || !entry->view_controller)
	return 0;

    if (entry->view_controller->hasPendingLodResults())
	(void)entry->view_controller->processPendingLodResults(max_results);
    if (entry->view_controller->isLodAutoSubmitEnabled())
	(void)entry->view_controller->submitLodRequestsIfNeeded();
    ged_obol_lod_status_fill(entry, status);
    return 1;
}

extern "C" int
ged_draw_obol_lod_service_status(struct ged *gedp,
				 void *view_ctx,
				 struct ged_draw_obol_lod_service_status *status)
{
    ged_obol_attached_controller *entry =
	ged_obol_attached_controller_for_view(gedp, view_ctx);
    if (!entry || !entry->view_controller)
	return 0;

    ged_obol_lod_status_fill(entry, status);
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
ged_obol_lod_prewarm_submit_name(BRLObolLodService *service,
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

    BRLObolMeshLodProvider *provider = new BRLObolMeshLodProvider;
    provider->dbip = dbip;
    provider->useView = FALSE;
    provider->refreshMissing = TRUE;
    provider->shrinkAfterCopy = TRUE;

    BRLObolLodTask task;
    task.generation = generation;
    task.request.databaseId = database_id ? database_id : "";
    task.request.objectPath = dp->d_namep;
    task.request.objectName = dp->d_namep;
    task.request.providerId = "brlobol_mesh_lod_cache";
    task.request.providerVersion = "brlobol-cache-v1";
    task.request.qualityTier = BRLOBOL_LOD_QUALITY_FAST_DISPLAY;
    task.realize = brlobol_mesh_lod_cache_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = brlobol_mesh_lod_provider_free;
    task.publishResult = FALSE;

    uint64_t taskId = service->submitIfNotActive(task);
    if (taskId == 0) {
	brlobol_mesh_lod_provider_free(provider);
	return 0;
    }

    return 1;
}

extern "C" size_t
ged_draw_obol_lod_service_prewarm(struct ged *gedp,
				  void *view_ctx,
				  int argc,
				  const char * const *argv,
				  struct ged_draw_obol_lod_service_status *status)
{
    ged_obol_attached_controller *entry =
	ged_obol_attached_controller_for_view(gedp, view_ctx);
    if (!entry || !entry->view_controller)
	return 0;

    if (!entry->lod_service || !entry->lod_service->isRunning()) {
	if (!ged_draw_obol_lod_service_start(gedp, view_ctx, 0))
	    return 0;
	entry = ged_obol_attached_controller_for_view(gedp, view_ctx);
	if (!entry || !entry->lod_service || !entry->lod_service->isRunning())
	    return 0;
    }

    struct db_i *dbip = gedp ? gedp->dbip : NULL;
    if (!dbip)
	return 0;

    const char *database_id = dbip->dbi_filename ? dbip->dbi_filename : "";
    uint64_t generation = entry->lod_service->beginGeneration();
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
			     entry->lod_service, dbip, database_id, generation,
			     dp->d_namep);
	}
	if (paths)
	    bu_free(paths, "free obol lod prewarm db_ls output");
    } else {
	for (int i = 0; i < argc; i++)
	    submitted += ged_obol_lod_prewarm_submit_name(
			     entry->lod_service, dbip, database_id, generation, argv[i]);
    }

    ged_obol_lod_status_fill(entry, status);
    return submitted;
}

int
ged_draw_obol_scene_controller_owned(struct ged *gedp)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    return (gdp && gdp->gd_obol_scene_controller &&
	    (gdp->gd_obol_scene_controller_owned ||
	     gdp->gd_obol_controller_owned)) ? 1 : 0;
}

extern "C" int
ged_draw_obol_scene_controller_attached(struct ged *gedp)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    std::vector<ged_obol_attached_controller> *entries =
	ged_obol_attached_controllers(gdp, 0);
    return (ged_draw_obol_scene_controller(gedp) ||
	    (entries && !entries->empty())) ? 1 : 0;
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (scene && ged_draw_obol_scene_controller_owned(gedp)) {
	if (sync_current_scene && gdp &&
	    !gdp->gd_obol_scene_controller_full_sync) {
	    (void)ged_draw_obol_scene_sync_full_scene(gedp, NULL, 0, scene);
	    gdp->gd_obol_scene_controller_full_sync = 1;
	}
	return 1;
    }

    brlobol_init(NULL);

    SoBRLSceneGroup *root = new SoBRLSceneGroup;
    BRLObolViewController *owned_controller = new BRLObolViewController(root);
    SoBRLSceneController *owned_scene = owned_controller->getSceneController();
    std::vector<ged_obol_attached_controller> *entries =
	ged_obol_attached_controllers(gdp, 1);
    if (!entries) {
	delete owned_controller;
	return 0;
    }

    ged_obol_attached_controller entry;
    entry.scene_controller = owned_scene;
    entry.view_controller = owned_controller;
    entry.owned_scene_controller = 0;
    entry.owned_view_controller = 1;
    entry.full_sync = 0;
    entry.use_attached_view_scope = 0;
    entry.render_endpoint_only = 0;
    entries->insert(entries->begin(), entry);
    ged_obol_primary_set(gdp, &entries->front());

    if (sync_current_scene) {
	(void)ged_draw_obol_scene_sync_full_scene(gedp, NULL, 0,
		owned_scene);
	entries->front().full_sync = 1;
	ged_obol_primary_set(gdp, &entries->front());
    }

    return 1;
}

extern "C" int
ged_draw_obol_scene_controller_full_synced(struct ged *gedp)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (gdp && gdp->gd_obol_scene_controller &&
	gdp->gd_obol_scene_controller_full_sync &&
	ged_draw_obol_scene_controller_owned(gedp))
	return 1;
    return 0;
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

static fastf_t
ged_obol_reported_transparency(float transparency)
{
    const double scaled = (double)transparency * 1000000.0;
    return (fastf_t)(floor(scaled + 0.5) / 1000000.0);
}

static BRLObolViewController *
ged_obol_view_controller_for_context(void *view_ctx)
{
    struct ged *gedp = view_ctx ?
			   static_cast<struct ged *>(ged_view_context_user_data_get(view_ctx)) :
		       NULL;
    if (!gedp)
	return NULL;
    ged_obol_attached_controller *entry =
	ged_obol_attached_controller_for_view(gedp, view_ctx);
    if (entry && entry->view_controller)
	return entry->view_controller;
    return ged_draw_obol_controller(gedp);
}

static BRLObolViewController *
ged_obol_view_controller_ensure_for_context(void *view_ctx,
	int sync_current_scene)
{
    struct ged *gedp = view_ctx ?
			   static_cast<struct ged *>(ged_view_context_user_data_get(view_ctx)) :
		       NULL;
    if (!gedp)
	return NULL;
    ged_obol_attached_controller *entry =
	ged_obol_attached_controller_for_view(gedp, view_ctx);
    if (entry && entry->view_controller)
	return entry->view_controller;
    BRLObolViewController *controller = ged_draw_obol_controller(gedp);
    if (controller)
	return controller;
    if (!ged_draw_obol_scene_controller_ensure_owned(gedp, sync_current_scene))
	return NULL;
    return ged_draw_obol_controller(gedp);
}

static BRLObolViewController *
ged_obol_shared_view_controller_ensure_for_context(void *view_ctx,
	int sync_current_scene)
{
    struct ged *gedp = view_ctx ?
			   static_cast<struct ged *>(ged_view_context_user_data_get(view_ctx)) :
		       NULL;
    if (!gedp)
	return NULL;
    BRLObolViewController *controller = ged_draw_obol_controller(gedp);
    if (controller)
	return controller;
    if (!ged_draw_obol_scene_controller_ensure_owned(gedp, sync_current_scene))
	return NULL;
    return ged_draw_obol_controller(gedp);
}

static BRLObolViewController *
ged_obol_view_controller_for_scope(void *view_ctx,
				   int local,
				   int sync_current_scene)
{
    return local ?
	   ged_obol_view_controller_ensure_for_context(view_ctx,
		   sync_current_scene) :
	   ged_obol_shared_view_controller_ensure_for_context(view_ctx,
		   sync_current_scene);
}

static BRLObolFeatureStyle
ged_obol_feature_style_from_ged(
    const struct ged_draw_view_feature_style *style)
{
    BRLObolFeatureStyle out;
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

static void
ged_obol_feature_style_to_ged(
    struct ged_draw_view_feature_style *dst,
    const BRLObolFeatureStyle &src)
{
    if (!dst)
	return;

    struct ged_draw_view_feature_style init =
	    GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    *dst = init;
    if (src.hasVisible)
	dst->visible = src.visible ? 1 : 0;
    if (src.hasSelectable)
	dst->selectable = src.selectable ? 1 : 0;
    if (src.hasColor) {
	dst->color_valid = 1;
	ged_obol_rgb_from_color(src.color, dst->color);
    }
    if (src.hasLineWidth)
	dst->line_width = src.lineWidth;
    if (src.hasLineStyle)
	dst->line_style = src.lineStyle;
    if (src.hasArrow)
	dst->arrow = src.arrow ? 1 : 0;
    if (src.hasArrowTip) {
	dst->arrow_tip_length = src.arrowTipLength;
	dst->arrow_tip_width = src.arrowTipWidth;
    }
}

static int32_t
ged_obol_line_command_from_ged(int command, size_t index)
{
    if (command == GED_DRAW_VIEW_LINE_DRAW ||
	command == BG_GEOMETRY_LINE_DRAW)
	return static_cast<int32_t>(BRLObolLineCommand::Draw);
    if (command == GED_DRAW_VIEW_LINE_POINT_DRAW ||
	command == BG_GEOMETRY_POINT_DRAW)
	return static_cast<int32_t>(BRLObolLineCommand::Point);
    if (command == GED_DRAW_VIEW_LINE_MOVE ||
	command == BG_GEOMETRY_LINE_MOVE)
	return static_cast<int32_t>(BRLObolLineCommand::Move);
    return static_cast<int32_t>(index ? BRLObolLineCommand::Draw :
				BRLObolLineCommand::Move);
}

static int
ged_obol_line_command_to_ged(int32_t command)
{
    if (command == static_cast<int32_t>(BRLObolLineCommand::Draw))
	return GED_DRAW_VIEW_LINE_DRAW;
    if (command == static_cast<int32_t>(BRLObolLineCommand::Point))
	return GED_DRAW_VIEW_LINE_POINT_DRAW;
    return GED_DRAW_VIEW_LINE_MOVE;
}

static std::vector<SbVec3f>
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

static std::vector<int32_t>
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

static std::vector<int32_t>
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

static std::vector<SbVec3f>
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

static BRLObolFeatureHandle
ged_obol_feature_handle(BRLObolViewController *controller,
			void *view_ctx,
			const char *name)
{
    if (!controller || !name)
	return BRLObolFeatureHandle();

    BRLObolFeatureOwner owner;
    owner.ownerToken = view_ctx;
    owner.ownerId = ged_obol_view_scope_name(view_ctx).c_str();
    owner.ownerRole = "view";

    BRLObolFeatureHandle local =
	controller->features().findOwned(name, BRLOBOL_FEATURE_SCOPE_LOCAL,
					 &owner);
    if (local.isValid())
	return local;

    BRLObolFeatureHandle shared =
	controller->features().find(name, BRLOBOL_FEATURE_SCOPE_SHARED);
    if (shared.isValid())
	return shared;

    return controller->features().find(name);
}

struct ged_obol_feature_lookup {
    BRLObolViewController *controller;
    BRLObolFeatureHandle handle;
};

static ged_obol_feature_lookup
ged_obol_feature_lookup_for_context(void *view_ctx, const char *name)
{
    ged_obol_feature_lookup out;
    out.controller = NULL;
    out.handle = BRLObolFeatureHandle();
    if (!name)
	return out;

    BRLObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (local_controller) {
	BRLObolFeatureHandle handle =
	    ged_obol_feature_handle(local_controller, view_ctx, name);
	if (handle.isValid()) {
	    out.controller = local_controller;
	    out.handle = handle;
	    return out;
	}
    }

    BRLObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if (shared_controller && shared_controller != local_controller) {
	BRLObolFeatureHandle handle =
	    shared_controller->features().find(name,
					       BRLOBOL_FEATURE_SCOPE_SHARED);
	if (handle.isValid()) {
	    out.controller = shared_controller;
	    out.handle = handle;
	}
    }

    return out;
}

static BRLObolFeatureScope
ged_obol_feature_scope(int local)
{
    return local ? BRLObolFeatureScope::Local : BRLObolFeatureScope::Shared;
}

static BRLObolFeatureOwner
ged_obol_feature_owner(void *view_ctx, int local)
{
    BRLObolFeatureOwner owner;
    if (!local)
	return owner;

    owner.ownerToken = view_ctx;
    owner.ownerId = ged_obol_view_scope_name(view_ctx).c_str();
    owner.ownerRole = "view";
    return owner;
}

static BRLObolOverlayInfo
ged_obol_model_overlay_info(void *view_ctx,
			    BRLObolOverlayClass overlay_class,
			    BRLObolOverlayLifecycle lifecycle,
			    BRLObolOverlayOrder order,
			    const char *source_path)
{
    BRLObolOverlayInfo info;
    info.isOverlay = TRUE;
    info.ownerToken = view_ctx;
    info.role = BRLObolOverlayRole::Model;
    info.overlayClass = overlay_class;
    info.lifecycle = lifecycle;
    info.order = order;
    info.sortOrder = 0;
    info.sourcePath = source_path ? source_path : "";
    return info;
}

static int
ged_obol_feature_mark_overlay(BRLObolViewController *controller,
			      BRLObolFeatureHandle handle,
			      const BRLObolOverlayInfo &overlay)
{
    if (!controller || !handle.isValid())
	return 0;
    return controller->features().setOverlayInfo(handle, overlay) ? 1 : 0;
}

#define GED_DRAW_COMMAND_SCENE_MAGIC 0x67646373u

struct ged_draw_command_scene {
    uint32_t magic;
    void *view_ctx;
    BRLObolViewController *controller;
    BRLObolFeatureOwner owner;
    BRLObolFeatureScope scope;
    ged_draw_command_result_cb result_cb;
    void *result_cb_data;
    int changed;
    int failed;
};

static int
ged_obol_command_scene_valid(const struct ged_draw_command_scene *scene)
{
    return scene && scene->magic == GED_DRAW_COMMAND_SCENE_MAGIC &&
	   scene->controller;
}

static void
ged_obol_command_scene_notify(
    struct ged_draw_command_scene *scene,
    int status,
    const char *name,
    const char *command,
    const char *diagnostic,
    BRLObolFeatureHandle handle = BRLObolFeatureHandle())
{
    if (!ged_obol_command_scene_valid(scene) || !scene->result_cb)
	return;

    struct ged_draw_command_result result =
	    GED_DRAW_COMMAND_RESULT_INIT;
    result.status = status;
    result.feature_name = name;
    result.command = command;
    result.diagnostic = diagnostic;
    if (handle.isValid()) {
	result.feature_id = handle.id;
	result.feature_revision = handle.revision;
    }
    scene->result_cb(&result, scene->result_cb_data);
}

static int
ged_obol_command_scene_writable(
    struct ged_draw_command_scene *scene,
    const char *name,
    const char *command)
{
    if (!ged_obol_command_scene_valid(scene))
	return 0;

    if (!scene->controller->features().commandOwnerGenerationCurrent(
	    scene->owner)) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name, command,
				      "stale command-scene generation");
	return 0;
    }

    return 1;
}

static BRLObolFeatureOwner
ged_obol_command_scene_owner(
    void *view_ctx,
    const struct ged_draw_command_scene_desc *desc)
{
    BRLObolFeatureOwner owner;
    const char *owner_id = (desc && desc->owner_id && desc->owner_id[0]) ?
			   desc->owner_id : "ged-command";
    const char *owner_role =
	(desc && desc->owner_role && desc->owner_role[0]) ?
	desc->owner_role : "command-scene";
    const char *run_id = (desc && desc->run_id && desc->run_id[0]) ?
			 desc->run_id : NULL;

    std::string id(owner_id);
    if (run_id) {
	id += "#";
	id += run_id;
    }
    owner.ownerToken = NULL;
    owner.ownerId = id.c_str();
    owner.ownerRole = owner_role;
    owner.generation = desc ? desc->generation : 0;
    if (owner.ownerId.getLength() == 0)
	owner.ownerId = ged_obol_view_scope_name(view_ctx).c_str();
    return owner;
}

static BRLObolOverlayInfo
ged_obol_command_scene_overlay_info(
    const struct ged_draw_command_scene *scene,
    const char *source_path)
{
    BRLObolOverlayInfo info =
	ged_obol_model_overlay_info(scene ? scene->view_ctx : NULL,
				    BRLObolOverlayClass::CommandResult,
				    BRLObolOverlayLifecycle::PerCommand,
				    BRLObolOverlayOrder::PostTransparent,
				    source_path);
    return info;
}

static int
ged_obol_command_scene_remove_feature(
    struct ged_draw_command_scene *scene,
    const char *name,
    const char *command)
{
    if (!ged_obol_command_scene_writable(scene, name, command) || !name)
	return 0;

    const unsigned int scope_mask =
	scene->scope == BRLObolFeatureScope::Local ?
	BRLOBOL_FEATURE_SCOPE_LOCAL : BRLOBOL_FEATURE_SCOPE_SHARED;
    BRLObolFeatureHandle handle =
	scene->controller->features().findOwned(name, scope_mask,
	    &scene->owner);
    const int removed = scene->controller->features().removeOwned(name,
			scope_mask, &scene->owner) ? 1 : 0;
    if (removed)
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_REMOVED, name,
				      command ? command : "remove", NULL, handle);
    return removed;
}

static BRLObolFeatureHandle
ged_obol_publish_line_set(BRLObolViewController *controller,
			  void *view_ctx,
			  const char *name,
			  int local,
			  const std::vector<SbVec3f> &points,
			  const std::vector<int32_t> &commands,
			  const BRLObolFeatureStyle *style)
{
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishLineSet(name,
	    ged_obol_feature_scope(local), points, commands, style,
	    local ? &owner : NULL);
}

static int
ged_obol_remove_feature(BRLObolViewController *controller,
			void *view_ctx,
			const char *name,
			int local_mode)
{
    if (!controller || !name)
	return 0;

    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BRLObolFeatureHandle handle;
    if (local_mode > 0) {
	handle = controller->features().findOwned(name,
		 BRLOBOL_FEATURE_SCOPE_LOCAL, &owner);
    } else if (local_mode == 0) {
	handle = controller->features().find(name,
					     BRLOBOL_FEATURE_SCOPE_SHARED);
    } else {
	handle = ged_obol_feature_handle(controller, view_ctx, name);
    }

    return handle.isValid() ? (controller->features().remove(handle) ? 1 : 0) : 0;
}

static BRLObolFeatureHandle
ged_obol_publish_indexed_face_set(BRLObolViewController *controller,
				  void *view_ctx,
				  const char *name,
				  int local,
				  const std::vector<SbVec3f> &points,
				  const std::vector<SbVec3f> &normals,
				  const std::vector<int32_t> &indices,
				  const BRLObolFeatureStyle *style)
{
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishIndexedFaceSet(name,
	    ged_obol_feature_scope(local), points, normals, indices, style,
	    local ? &owner : NULL);
}

static BRLObolFeatureHandle
ged_obol_publish_labels(BRLObolViewController *controller,
			void *view_ctx,
			const char *name,
			int local,
			const std::vector<BRLObolLabel> &labels,
			const BRLObolFeatureStyle *style = NULL)
{
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishLabels(name,
	    ged_obol_feature_scope(local), labels, style,
	    local ? &owner : NULL);
}

static BRLObolFeatureHandle
ged_obol_publish_arrow(BRLObolViewController *controller,
		       void *view_ctx,
		       const char *name,
		       int local,
		       const std::vector<SbVec3f> &points,
		       const BRLObolFeatureStyle *style)
{
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishArrow(name,
	    ged_obol_feature_scope(local), points, style,
	    local ? &owner : NULL);
}

static BRLObolFeatureHandle
ged_obol_publish_axes(BRLObolViewController *controller,
		      void *view_ctx,
		      const char *name,
		      int local,
		      const std::vector<SbVec3f> &centers,
		      float half_axes_size,
		      const BRLObolFeatureStyle *style)
{
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishAxes(name,
	    ged_obol_feature_scope(local), centers, half_axes_size, style,
	    local ? &owner : NULL);
}

static void
ged_obol_point_from_sb(point_t dst, const SbVec3f &src)
{
    VSET(dst, src[0], src[1], src[2]);
}

static BRLObolLabel
ged_obol_label_from_ged(const struct ged_draw_view_label_data &label)
{
    BRLObolLabel out;
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

static BRLObolLabel
ged_obol_label_from_hud(const struct ged_diagnostic_hud_label &label)
{
    BRLObolLabel out;
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

static std::vector<BRLObolLabel>
ged_obol_labels_from_ged(const struct ged_draw_view_label_data *labels,
			 size_t label_count)
{
    std::vector<BRLObolLabel> out;
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

static BRLObolFeatureStyle
ged_obol_faceplate_style(const int rgb[3],
			 int fallback_r,
			 int fallback_g,
			 int fallback_b,
			 int line_width)
{
    BRLObolFeatureStyle style;
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

static BRLObolOverlayInfo
ged_obol_faceplate_overlay_info(void *view_ctx,
				BRLObolOverlayOrder order =
				    BRLObolOverlayOrder::Screen)
{
    BRLObolOverlayInfo info;
    info.isOverlay = TRUE;
    info.ownerToken = view_ctx;
    info.role = BRLObolOverlayRole::Screen;
    info.overlayClass = BRLObolOverlayClass::Faceplate;
    info.lifecycle = BRLObolOverlayLifecycle::PerView;
    info.order = order;
    info.sortOrder = 0;
    info.sourcePath = "_faceplate";
    return info;
}

static int
ged_obol_view_to_model_point(SbVec3f &out,
			     void *view_ctx,
			     fastf_t x,
			     fastf_t y,
			     fastf_t z = 0.0)
{
    mat_t view2model;
    if (!rt_view_context_view2model_get(view2model, view_ctx))
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

static void
ged_obol_faceplate_append_line(std::vector<SbVec3f> &points,
			       std::vector<int32_t> &commands,
			       void *view_ctx,
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
    commands.push_back(static_cast<int32_t>(BRLObolLineCommand::Move));
    points.push_back(b);
    commands.push_back(static_cast<int32_t>(BRLObolLineCommand::Draw));
}

static BRLObolFeatureHandle
ged_obol_faceplate_publish_lines(BRLObolViewController *controller,
				 void *view_ctx,
				 const char *name,
				 const std::vector<SbVec3f> &points,
				 const std::vector<int32_t> &commands,
				 const BRLObolFeatureStyle &style)
{
    if (!controller || !name || points.empty() || commands.empty())
	return BRLObolFeatureHandle();

    BRLObolFeatureHandle handle = ged_obol_publish_line_set(controller,
				  view_ctx, name, 1, points, commands, &style);
    (void)ged_obol_feature_mark_overlay(controller, handle,
					ged_obol_faceplate_overlay_info(view_ctx));
    return handle;
}

static BRLObolFeatureHandle
ged_obol_faceplate_publish_hud_labels(BRLObolViewController *controller,
				      void *view_ctx,
				      const char *name,
				      const std::vector<BRLObolLabel> &labels,
				      const BRLObolFeatureStyle &style)
{
    if (!controller || !name || labels.empty())
	return BRLObolFeatureHandle();

    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BRLObolFeatureHandle handle =
	controller->features().publishHudLabels(name,
	    BRLObolFeatureScope::Local, labels, &style, &owner);
    (void)ged_obol_feature_mark_overlay(controller, handle,
					ged_obol_faceplate_overlay_info(view_ctx));
    return handle;
}

static void
ged_obol_faceplate_remove(BRLObolViewController *controller,
			  void *view_ctx,
			  const char *name)
{
    (void)ged_obol_remove_feature(controller, view_ctx, name, 1);
}

static BRLObolLabel
ged_obol_faceplate_label(void *view_ctx,
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
    BRLObolLabel label;
    const int width = rt_view_context_width_get(view_ctx);
    const int height = rt_view_context_height_get(view_ctx);
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
ged_obol_faceplate_sync_center_dot(BRLObolViewController *controller,
				   void *view_ctx)
{
    static const char name[] = "_faceplate/center_dot";
    struct rt_view_other_state state = RT_VIEW_OTHER_STATE_INIT;
    if (!rt_view_context_center_dot_state_get(&state, view_ctx) ||
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
    BRLObolFeatureStyle style = ged_obol_faceplate_style(
				    state.gos_line_color, 255, 255, 0, 1);
    (void)ged_obol_faceplate_publish_lines(controller, view_ctx, name,
					   points, commands, style);
}

static void
ged_obol_faceplate_sync_grid(BRLObolViewController *controller,
			     void *view_ctx)
{
    static const char name[] = "_faceplate/grid";
    struct rt_view_grid_state grid = RT_VIEW_GRID_STATE_INIT;
    if (!rt_view_context_grid_state_get(&grid, view_ctx) || !grid.draw) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    SoBRLGrid *node = new SoBRLGrid;
    node->ref();
    node->overlayId = name;
    if (!brlobol_grid_configure_from_view_context(node, &grid, view_ctx) ||
	node->getTotalSegmentCount() <= 0) {
	node->unref();
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    BRLObolFeatureStyle style;
    style.hasVisible = TRUE;
    style.visible = TRUE;
    style.hasSelectable = TRUE;
    style.selectable = FALSE;
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BRLObolFeatureHandle handle =
	controller->features().publishCustomNode(name,
	    BRLObolFeatureScope::Local, node, &style, &owner);
    node->unref();
    (void)ged_obol_feature_mark_overlay(controller, handle,
					ged_obol_faceplate_overlay_info(view_ctx));
}

static void
ged_obol_faceplate_params_string(void *view_ctx,
				 const struct rt_view_params_state *params,
				 struct bu_vls *vls)
{
    if (!view_ctx || !params || !vls)
	return;

    point_t center = VINIT_ZERO;
    mat_t center_mat;
    if (rt_view_context_center_get(center_mat, view_ctx))
	MAT_DELTAS_GET_NEG(center, center_mat);
    VSCALE(center, center, rt_view_context_base2local_get(view_ctx));

    const char *ustr = bu_units_string(rt_view_context_local2base_get(view_ctx));
    if (!ustr)
	ustr = "";

    vect_t aet = VINIT_ZERO;
    (void)rt_view_context_aet_get(aet, view_ctx);

    if (params->draw_size) {
	if (bu_vls_strlen(vls) > 0)
	    bu_vls_putc(vls, ' ');
	bu_vls_printf(vls, "size[%s]: %.2f", ustr,
		      rt_view_context_size_get(view_ctx) *
		      rt_view_context_base2local_get(view_ctx));
    }
    if (params->draw_center) {
	if (bu_vls_strlen(vls) > 0)
	    bu_vls_putc(vls, ' ');
	bu_vls_printf(vls, "center[%s]: (%.2f, %.2f, %.2f)",
		      ustr, V3ARGS(center));
    }
    if (params->draw_az) {
	if (bu_vls_strlen(vls) > 0)
	    bu_vls_putc(vls, ' ');
	bu_vls_printf(vls, "az:%.2f", aet[0]);
    }
    if (params->draw_el) {
	if (bu_vls_strlen(vls) > 0)
	    bu_vls_putc(vls, ' ');
	bu_vls_printf(vls, "el:%.2f", aet[1]);
    }
    if (params->draw_tw) {
	if (bu_vls_strlen(vls) > 0)
	    bu_vls_putc(vls, ' ');
	bu_vls_printf(vls, "tw:%.2f", aet[2]);
    }

    const uint64_t frametime = rt_view_context_frametime_get(view_ctx);
    if (params->draw_fps && frametime > 0) {
	if (bu_vls_strlen(vls) > 0)
	    bu_vls_putc(vls, ' ');
	bu_vls_printf(vls, "FPS:%.2f", 1.0 / (fastf_t)frametime);
    }
}

static void
ged_obol_faceplate_sync_params(BRLObolViewController *controller,
			       void *view_ctx)
{
    static const char name[] = "_faceplate/params";
    struct rt_view_params_state params = RT_VIEW_PARAMS_STATE_INIT;
    if (!rt_view_context_params_state_get(&params, view_ctx) ||
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

    struct bu_vls text = BU_VLS_INIT_ZERO;
    ged_obol_faceplate_params_string(view_ctx, &params, &text);
    std::vector<BRLObolLabel> labels;
    labels.push_back(ged_obol_faceplate_label(view_ctx, bu_vls_cstr(&text),
		     -0.98, -0.965, params.color, 255, 255, 0,
		     params.font_size > 0 ? params.font_size : 20));
    bu_vls_free(&text);

    BRLObolFeatureStyle style = ged_obol_faceplate_style(
				    params.color, 255, 255, 0, 1);
    (void)ged_obol_faceplate_publish_hud_labels(controller, view_ctx, name,
	    labels, style);
}

static void
ged_obol_faceplate_sync_scale(BRLObolViewController *controller,
			      void *view_ctx)
{
    static const char line_name[] = "_faceplate/scale";
    static const char label_name[] = "_faceplate/scale_labels";
    struct rt_view_other_state state = RT_VIEW_OTHER_STATE_INIT;
    if (!rt_view_context_scale_overlay_state_get(&state, view_ctx) ||
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

    BRLObolFeatureStyle line_style = ged_obol_faceplate_style(
					 state.gos_line_color, 255, 255, 0, 1);
    (void)ged_obol_faceplate_publish_lines(controller, view_ctx, line_name,
					   points, commands, line_style);

    struct bu_vls scale = BU_VLS_INIT_ZERO;
    const fastf_t base2local = rt_view_context_base2local_get(view_ctx);
    const char *unit = !ZERO(base2local) ? bu_units_string(1.0 / base2local) :
		       NULL;
    if (!unit)
	unit = "";
    bu_vls_printf(&scale, "%g%s",
		  rt_view_context_size_get(view_ctx) * 0.5 *
		  base2local,
		  unit);
    const int soffset = (int)(strlen(bu_vls_cstr(&scale)) * 0.5);
    std::vector<BRLObolLabel> labels;
    labels.push_back(ged_obol_faceplate_label(view_ctx, "0", -0.505, -0.78,
		     state.gos_line_color, 255, 255, 0,
		     state.gos_font_size > 0 ? state.gos_font_size : 20));
    labels.push_back(ged_obol_faceplate_label(view_ctx, bu_vls_cstr(&scale),
		     0.5 - (soffset * 0.015), -0.78,
		     state.gos_line_color, 255, 255, 0,
		     state.gos_font_size > 0 ? state.gos_font_size : 20));
    bu_vls_free(&scale);

    BRLObolFeatureStyle text_style = ged_obol_faceplate_style(
					 state.gos_line_color, 255, 255, 0, 1);
    (void)ged_obol_faceplate_publish_hud_labels(controller, view_ctx,
	    label_name, labels, text_style);
}

static void
ged_obol_faceplate_append_axis(std::vector<SbVec3f> &points,
			       std::vector<int32_t> &commands,
			       std::vector<BRLObolLabel> &labels,
			       void *view_ctx,
			       const mat_t rmat,
			       const struct rt_view_axes_state *axes,
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
ged_obol_faceplate_remove_axis_variants(BRLObolViewController *controller,
					void *view_ctx,
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
				       void *view_ctx,
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
				     void *view_ctx,
				     const mat_t rmat,
				     const struct rt_view_axes_state *axes,
				     fastf_t aspect)
{
    if (!view_ctx || !axes || !axes->tick_enabled ||
	axes->tick_interval <= SMALL_FASTF)
	return;

    const fastf_t view_size = rt_view_context_size_get(view_ctx);
    const int width = rt_view_context_width_get(view_ctx);
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
ged_obol_faceplate_sync_axes_one(BRLObolViewController *controller,
				 void *view_ctx,
				 const char *prefix,
				 struct rt_view_axes_state axes,
				 int model_axes)
{
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
    if (!model_axes && VNEAR_ZERO(axes.axes_pos, SMALL_FASTF)) {
	VSET(axes.axes_pos, 0.80, -0.80, 0.0);
	axes.pos_only = 1;
	axes.triple_color = 1;
	axes.label_flag = 1;
    }
    if (model_axes)
	axes.label_flag = 1;

    mat_t rmat;
    if (!rt_view_context_rotation_get(rmat, view_ctx)) {
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
	if (rt_view_context_model2view_get(model2view, view_ctx)) {
	    VSCALE(map, axes.axes_pos,
		   rt_view_context_local2base_get(view_ctx));
	    MAT4X3PNT(axes.axes_pos, model2view, map);
	}
    }

    const int width = rt_view_context_width_get(view_ctx);
    const int height = rt_view_context_height_get(view_ctx);
    const fastf_t aspect = (width > 0 && height > 0) ?
			   (fastf_t)width / (fastf_t)height : 1.0;
    if (!model_axes)
	axes.axes_pos[Y] *= (aspect > SMALL_FASTF) ? 1.0 / aspect : 1.0;

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    std::vector<BRLObolLabel> labels;
    points.reserve(6);
    commands.reserve(6);
    labels.reserve(3);

    if (axes.triple_color) {
	ged_obol_faceplate_remove(controller, view_ctx, line_name.c_str());
	for (int axis = X; axis <= Z; axis++) {
	    struct rt_view_axes_state axis_axes = axes;
	    int axis_color[3] = {0, 0, 0};
	    ged_obol_faceplate_axis_triple_color(axis, axis_color);
	    VMOVE(axis_axes.axes_color, axis_color);
	    VMOVE(axis_axes.label_color, axis_color);

	    std::vector<SbVec3f> axis_points;
	    std::vector<int32_t> axis_commands;
	    std::vector<BRLObolLabel> axis_labels;
	    axis_points.reserve(2);
	    axis_commands.reserve(2);
	    ged_obol_faceplate_append_axis(axis_points, axis_commands,
					   axis_labels, view_ctx, rmat, &axis_axes, axis, aspect,
					   axis == X ? "X" : axis == Y ? "Y" : "Z");

	    BRLObolFeatureStyle axis_style = ged_obol_faceplate_style(
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

	BRLObolFeatureStyle line_style = ged_obol_faceplate_style(
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
	BRLObolFeatureStyle tick_style = ged_obol_faceplate_style(
					     axes.tick_color, 255, 255, 0, axes.line_width);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx,
					       tick_name.c_str(), tick_points, tick_commands, tick_style);
    } else {
	ged_obol_faceplate_remove(controller, view_ctx, tick_name.c_str());
    }
    if (!major_tick_points.empty()) {
	BRLObolFeatureStyle major_tick_style = ged_obol_faceplate_style(
		axes.tick_major_color, 255, 0, 0, axes.line_width);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx,
					       major_tick_name.c_str(), major_tick_points,
					       major_tick_commands, major_tick_style);
    } else {
	ged_obol_faceplate_remove(controller, view_ctx,
				  major_tick_name.c_str());
    }

    if (!labels.empty()) {
	BRLObolFeatureStyle label_style = ged_obol_faceplate_style(
					      axes.label_color, 255, 255, 0, 1);
	(void)ged_obol_faceplate_publish_hud_labels(controller, view_ctx,
		label_name.c_str(), labels, label_style);
    } else {
	ged_obol_faceplate_remove(controller, view_ctx, label_name.c_str());
    }
}

static void
ged_obol_faceplate_sync_axes(BRLObolViewController *controller,
			     void *view_ctx)
{
    struct rt_view_axes_state axes = RT_VIEW_AXES_STATE_INIT;
    if (rt_view_context_view_axes_state_get(&axes, view_ctx))
	ged_obol_faceplate_sync_axes_one(controller, view_ctx,
					 "_faceplate/view_axes", axes, 0);
    if (rt_view_context_model_axes_state_get(&axes, view_ctx))
	ged_obol_faceplate_sync_axes_one(controller, view_ctx,
					 "_faceplate/model_axes", axes, 1);
}

static void
ged_obol_faceplate_sync_framebuffer(BRLObolViewController *controller,
				    void *view_ctx)
{
    static const char name[] = "_faceplate/framebuffer";
    if (rt_view_context_framebuffer_mode_get(view_ctx) == 0) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    struct dm *dmp =
	    static_cast<struct dm *>(rt_view_context_display_manager_get(view_ctx));
    struct fb *fbp = dmp ? dm_get_fb(dmp) : NULL;
    if (!fbp) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    const int width = fb_getwidth(fbp);
    const int height = fb_getheight(fbp);
    if (width <= 0 || height <= 0) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    icv_image_t *image = fb_write_icv(fbp, 0, 0, width, height);
    if (!image) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    SoBRLImageSource *source = new SoBRLImageSource;
    source->ref();
    const int source_ok = source->setImage(image);
    icv_destroy(image);
    if (source_ok != 0) {
	source->unref();
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    SoBRLViewportImage *viewport = new SoBRLViewportImage;
    viewport->ref();
    viewport->overlayId = name;
    viewport->imageSource = source;
    source->unref();
    viewport->visible = TRUE;
    viewport->layer = SoBRLViewportImage::OVERLAY;
    viewport->zOrder = 0;
    viewport->anchor = SoBRLViewportImage::LOWER_LEFT;
    viewport->units = SoBRLViewportImage::PIXELS;
    viewport->position = SbVec2f(0.0f, 0.0f);
    viewport->size = SbVec2f(
			 static_cast<float>(rt_view_context_width_get(view_ctx) > 0 ?
					    rt_view_context_width_get(view_ctx) : width),
			 static_cast<float>(rt_view_context_height_get(view_ctx) > 0 ?
					    rt_view_context_height_get(view_ctx) : height));
    viewport->fit = SoBRLViewportImage::STRETCH;
    viewport->preserveAspect = FALSE;
    viewport->opacity = 1.0f;
    if (viewport->rebuildGeometry() != 0) {
	viewport->unref();
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    BRLObolFeatureStyle style;
    style.hasVisible = TRUE;
    style.visible = TRUE;
    style.hasSelectable = TRUE;
    style.selectable = FALSE;
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BRLObolFeatureHandle handle =
	controller->features().publishCustomNode(name,
	    BRLObolFeatureScope::Local, viewport, &style, &owner);
    viewport->unref();
    (void)ged_obol_feature_mark_overlay(controller, handle,
					ged_obol_faceplate_overlay_info(view_ctx,
						BRLObolOverlayOrder::PostTransparent));
}

extern "C" int
ged_draw_obol_view_context_faceplate_sync(struct ged *gedp, void *view_ctx)
{
    if (!gedp || !view_ctx)
	return BRLCAD_OK;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return BRLCAD_OK;

    ged_obol_faceplate_sync_center_dot(controller, view_ctx);
    ged_obol_faceplate_sync_grid(controller, view_ctx);
    ged_obol_faceplate_sync_params(controller, view_ctx);
    ged_obol_faceplate_sync_scale(controller, view_ctx);
    ged_obol_faceplate_sync_axes(controller, view_ctx);
    ged_obol_faceplate_sync_framebuffer(controller, view_ctx);

    if (rt_view_context_framebuffer_mode_get(view_ctx) != 0)
	(void)ged_draw_obol_framebuffer_present(gedp);

    return BRLCAD_OK;
}

extern "C" int
ged_draw_obol_view_context_feature_store_active(void *view_ctx)
{
    return ged_obol_view_controller_for_context(view_ctx) ? 1 : 0;
}

extern "C" size_t
ged_draw_obol_view_context_clear(void *view_ctx, int flags)
{
    if (!(flags & GED_VIEW_CLEAR_VIEW))
	return 0;

    size_t removed = 0;
    if (flags & GED_VIEW_CLEAR_LOCAL) {
	BRLObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	if (!controller)
	    return 0;
	BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
	removed += controller->features().removeScope(
		       BRLOBOL_FEATURE_SCOPE_LOCAL, &owner);
	removed += controller->polygons().removeScope(
		       BRLOBOL_FEATURE_SCOPE_LOCAL);
	controller->selection().clear(&owner, BRLOBOL_SELECTION_ALL);
    } else {
	BRLObolViewController *controller =
	    ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
	if (!controller)
	    return 0;
	removed += controller->features().removeScope(
		       BRLOBOL_FEATURE_SCOPE_SHARED, NULL);
	removed += controller->polygons().removeScope(
		       BRLOBOL_FEATURE_SCOPE_SHARED);
	controller->selection().clear(NULL, BRLOBOL_SELECTION_ALL);
    }
    return removed;
}

extern "C" int
ged_draw_obol_view_context_feature_exists(void *view_ctx, const char *name)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    if (!lookup.controller || !name)
	return 0;
    return lookup.handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_feature_remove(void *view_ctx, const char *name)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    if (!lookup.controller || !lookup.handle.isValid())
	return 0;
    return lookup.controller->features().remove(lookup.handle) ? 1 : 0;
}

extern "C" size_t
ged_draw_obol_view_context_features_remove_prefix(
    void *view_ctx,
    const char *prefix)
{
    if (!prefix || !prefix[0])
	return 0;

    size_t removed = 0;
    BRLObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (local_controller) {
	BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
	removed += local_controller->features().removePrefix(prefix,
		   BRLOBOL_FEATURE_SCOPE_LOCAL, &owner);
    }

    BRLObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if (shared_controller)
	removed += shared_controller->features().removePrefix(prefix,
		   BRLOBOL_FEATURE_SCOPE_SHARED, NULL);
    return removed;
}

extern "C" int
ged_draw_obol_view_context_feature_summary(
    void *view_ctx,
    const char *name,
    struct ged_draw_view_feature_summary *summary)
{
    if (!summary)
	return 0;

    memset(summary, 0, sizeof(*summary));
    if (!name)
	return 0;

    BRLObolFeatureSummary obol_summary;
    BRLObolViewController *summary_controller = NULL;

    BRLObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (local_controller) {
	BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
	if (!local_controller->features().summaryOwned(name, obol_summary,
		BRLOBOL_FEATURE_SCOPE_LOCAL, &owner))
	    return 0;
	if (obol_summary.exists)
	    summary_controller = local_controller;
    }

    if (!obol_summary.exists) {
	BRLObolViewController *shared_controller =
	    ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
	if (!shared_controller ||
	    !shared_controller->features().summary(name, obol_summary,
		    BRLOBOL_FEATURE_SCOPE_SHARED))
	    return 0;
	if (obol_summary.exists)
	    summary_controller = shared_controller;
    }

    if (!obol_summary.exists)
	return 1;

    summary->exists = 1;
    summary->visible = obol_summary.visible ? 1 : 0;
    summary->child_count = obol_summary.childCount;
    summary->geometry_command_count = obol_summary.commandCount;
    summary->metadata_count = obol_summary.metadataCount;
    summary->primitive_metadata_count = obol_summary.primitiveMetadataCount;
    summary->selected_primitive_count =
	obol_summary.selectedPrimitiveCount;
    summary->highlighted_primitive_count =
	obol_summary.highlightedPrimitiveCount;
    summary->is_label =
	(obol_summary.kind == BRLObolFeatureKind::Labels ||
	 obol_summary.kind == BRLObolFeatureKind::HudLabel) ? 1 : 0;
    summary->is_transient_preview =
	(obol_summary.kind == BRLObolFeatureKind::EditPreview) ? 1 : 0;
    summary->is_command_result =
	(obol_summary.overlay.overlayClass ==
	 BRLObolOverlayClass::CommandResult) ? 1 : 0;
    summary->is_overlay = (obol_summary.overlay.isOverlay ||
			   (!summary->is_label &&
			    !summary->is_transient_preview)) ? 1 : 0;

    BRLObolFeatureHandle handle = summary_controller ?
				  ged_obol_feature_handle(summary_controller, view_ctx, name) :
				  BRLObolFeatureHandle();
    BRLObolFeatureStyle style;
    if (summary_controller && handle.isValid() &&
	summary_controller->features().style(handle, style) &&
	style.hasColor)
	ged_obol_rgb_from_color(style.color, summary->color);

    return 1;
}

static std::vector<int32_t>
ged_obol_primitives_from_ged(const int *primitives, size_t primitive_count)
{
    std::vector<int32_t> out;
    if (!primitives || !primitive_count)
	return out;

    out.reserve(primitive_count);
    for (size_t i = 0; i < primitive_count; i++)
	if (primitives[i] >= 0)
	    out.push_back(static_cast<int32_t>(primitives[i]));
    return out;
}

extern "C" int
ged_draw_obol_view_context_feature_selected_primitives_replace(
    void *view_ctx,
    const char *name,
    const int *primitives,
    size_t primitive_count)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    if (!handle.isValid())
	return 0;

    std::vector<int32_t> obol_primitives =
	ged_obol_primitives_from_ged(primitives, primitive_count);
    return controller->features().replaceSelectedPrimitives(handle,
	    obol_primitives) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_feature_highlighted_primitives_replace(
    void *view_ctx,
    const char *name,
    const int *primitives,
    size_t primitive_count)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    if (!handle.isValid())
	return 0;

    std::vector<int32_t> obol_primitives =
	ged_obol_primitives_from_ged(primitives, primitive_count);
    return controller->features().replaceHighlightedPrimitives(handle,
	    obol_primitives) ? 1 : 0;
}

extern "C" size_t
ged_draw_obol_view_context_feature_selected_primitive_count(
    void *view_ctx,
    const char *name)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    std::vector<int32_t> primitives;
    if (!handle.isValid() ||
	!controller->features().selectedPrimitives(handle, primitives))
	return 0;
    return primitives.size();
}

extern "C" size_t
ged_draw_obol_view_context_feature_highlighted_primitive_count(
    void *view_ctx,
    const char *name)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    std::vector<int32_t> primitives;
    if (!handle.isValid() ||
	!controller->features().highlightedPrimitives(handle, primitives))
	return 0;
    return primitives.size();
}

extern "C" int
ged_draw_obol_view_context_feature_selected_primitive_at(
    void *view_ctx,
    const char *name,
    size_t index,
    int *primitive)
{
    if (primitive)
	*primitive = -1;
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    std::vector<int32_t> primitives;
    if (!primitive || !handle.isValid() ||
	!controller->features().selectedPrimitives(handle, primitives) ||
	index >= primitives.size())
	return 0;
    *primitive = static_cast<int>(primitives[index]);
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_highlighted_primitive_at(
    void *view_ctx,
    const char *name,
    size_t index,
    int *primitive)
{
    if (primitive)
	*primitive = -1;
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    std::vector<int32_t> primitives;
    if (!primitive || !handle.isValid() ||
	!controller->features().highlightedPrimitives(handle,
		primitives) ||
	index >= primitives.size())
	return 0;
    *primitive = static_cast<int>(primitives[index]);
    return 1;
}

extern "C" size_t
ged_draw_obol_view_context_feature_metadata_count(void *view_ctx,
	const char *name)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    std::vector<BRLObolFeatureMetadata> metadata;
    if (!handle.isValid() ||
	!controller->features().metadata(handle, metadata))
	return 0;
    return metadata.size();
}

extern "C" int
ged_draw_obol_view_context_feature_metadata_copy(
    void *view_ctx,
    const char *name,
    size_t index,
    struct bu_vls *key,
    struct bu_vls *value)
{
    if (key)
	bu_vls_trunc(key, 0);
    if (value)
	bu_vls_trunc(value, 0);

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    std::vector<BRLObolFeatureMetadata> metadata;
    if (!handle.isValid() ||
	!controller->features().metadata(handle, metadata) ||
	index >= metadata.size())
	return 0;

    if (key)
	bu_vls_strcat(key, metadata[index].key.getString());
    if (value)
	bu_vls_strcat(value, metadata[index].value.getString());
    return 1;
}

extern "C" size_t
ged_draw_obol_view_context_feature_primitive_metadata_count(
    void *view_ctx,
    const char *name,
    int primitive)
{
    if (primitive < 0)
	return 0;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    std::vector<BRLObolFeatureMetadata> metadata;
    if (!handle.isValid() ||
	!controller->features().primitiveMetadata(handle,
		static_cast<int32_t>(primitive), metadata))
	return 0;
    return metadata.size();
}

extern "C" int
ged_draw_obol_view_context_feature_primitive_metadata_copy(
    void *view_ctx,
    const char *name,
    int primitive,
    size_t index,
    struct bu_vls *key,
    struct bu_vls *value)
{
    if (key)
	bu_vls_trunc(key, 0);
    if (value)
	bu_vls_trunc(value, 0);
    if (primitive < 0)
	return 0;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    std::vector<BRLObolFeatureMetadata> metadata;
    if (!handle.isValid() ||
	!controller->features().primitiveMetadata(handle,
		static_cast<int32_t>(primitive), metadata) ||
	index >= metadata.size())
	return 0;

    if (key)
	bu_vls_strcat(key, metadata[index].key.getString());
    if (value)
	bu_vls_strcat(value, metadata[index].value.getString());
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_pick_primitive_resolve(
    void *view_ctx,
    const char *picked_feature_name,
    int picked_primitive,
    int select,
    int highlight,
    struct bu_vls *feature_name,
    int *feature_primitive)
{
    if (feature_name)
	bu_vls_trunc(feature_name, 0);
    if (feature_primitive)
	*feature_primitive = -1;
    if (!picked_feature_name || picked_primitive < 0 ||
	!feature_primitive)
	return 0;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;

    BRLObolFeaturePrimitivePick pick;
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    if (!controller->features().resolvePrimitivePick(picked_feature_name,
	    static_cast<int32_t>(picked_primitive), pick,
	    BRLOBOL_FEATURE_SCOPE_LOCAL, &owner) &&
	!controller->features().resolvePrimitivePick(picked_feature_name,
		static_cast<int32_t>(picked_primitive), pick,
		BRLOBOL_FEATURE_SCOPE_SHARED, NULL))
	return 0;

    std::vector<int32_t> primitives;
    primitives.push_back(pick.primitiveIndex);
    if (select && !controller->features().replaceSelectedPrimitives(
	    pick.handle, primitives))
	return 0;
    if (highlight && !controller->features().replaceHighlightedPrimitives(
	    pick.handle, primitives))
	return 0;

    if (feature_name)
	bu_vls_strcat(feature_name, pick.featureName.getString());
    *feature_primitive = static_cast<int>(pick.primitiveIndex);
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_visible(void *view_ctx, const char *name)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    BRLObolFeatureStyle style;
    if (!handle.isValid() || !controller->features().style(handle, style))
	return 0;
    return (!style.hasVisible || style.visible) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_feature_visible_set(
    void *view_ctx,
    const char *name,
    int visible)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    if (!handle.isValid())
	return 0;
    return controller->features().setVisible(handle,
	    visible ? TRUE : FALSE) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_feature_style_get(
    void *view_ctx,
    const char *name,
    struct ged_draw_view_feature_style *style)
{
    if (!style)
	return 0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    BRLObolFeatureStyle obol_style;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().style(lookup.handle,
		obol_style))
	return 0;
    ged_obol_feature_style_to_ged(style, obol_style);
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_style_apply(
    void *view_ctx,
    const char *name,
    const struct ged_draw_view_feature_style *style,
    int recursive)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    if (!lookup.controller || !lookup.handle.isValid() || !style)
	return 0;

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    return lookup.controller->features().applyStyle(lookup.handle, obol_style,
	    recursive ? TRUE : FALSE) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_feature_realize(void *view_ctx,
	const char *name,
	int recursive)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    if (!lookup.handle.isValid() || !lookup.controller)
	return 0;
    return lookup.controller->features().realize(lookup.handle,
	    recursive ? TRUE : FALSE) ? 1 : 0;
}

static enum ged_draw_obol_view_feature_kind
ged_obol_view_feature_kind(BRLObolFeatureKind kind) {
    switch (kind)
    {
	case BRLObolFeatureKind::Lines:
		    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINES;
	case BRLObolFeatureKind::IndexedLines:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_LINES;
	case BRLObolFeatureKind::Points:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_POINTS;
	case BRLObolFeatureKind::Labels:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_LABELS;
	case BRLObolFeatureKind::Arrow:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_ARROW;
	case BRLObolFeatureKind::Axes:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_AXES;
	case BRLObolFeatureKind::LineLayer:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINE_LAYER;
	case BRLObolFeatureKind::EditPreview:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_EDIT_PREVIEW;
	case BRLObolFeatureKind::IndexedFaceSet:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_FACE_SET;
	case BRLObolFeatureKind::PolygonOverlay:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_POLYGON_OVERLAY;
	case BRLObolFeatureKind::HudLabel:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_HUD_LABEL;
	case BRLObolFeatureKind::CustomNode:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_CUSTOM_NODE;
	case BRLObolFeatureKind::Unknown:
	default:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_UNKNOWN;
    }
}

static void
ged_obol_depth_consider(void *view_ctx,
			const SbVec3f &point,
			int mode,
			int &have_depth,
			fastf_t &depth)
{
    mat_t model2view;
    ged_view_context_model2view_get(model2view, view_ctx);

    point_t model_pt;
    point_t view_pt;
    ged_obol_point_from_sb(model_pt, point);
    MAT4X3PNT(view_pt, model2view, model_pt);

    if (!have_depth) {
	depth = view_pt[Z];
	have_depth = 1;
	return;
    }

    if (mode) {
	if (view_pt[Z] > depth)
	    depth = view_pt[Z];
    } else if (view_pt[Z] < depth) {
	depth = view_pt[Z];
    }
}

static int
ged_obol_feature_record_depth(void *view_ctx,
			      const BRLObolFeatureRecord &record,
			      int mode,
			      fastf_t *depth)
{
    if (depth)
	*depth = 0.0;
    if (!view_ctx || !depth)
	return 0;

    int have_depth = 0;
    fastf_t calc_depth = mode ? -DBL_MAX : DBL_MAX;
    for (size_t i = 0; i < record.points.size(); i++)
	ged_obol_depth_consider(view_ctx, record.points[i], mode,
				have_depth, calc_depth);
    for (size_t i = 0; i < record.labels.size(); i++) {
	ged_obol_depth_consider(view_ctx, record.labels[i].point, mode,
				have_depth, calc_depth);
	if (record.labels[i].hasLeader)
	    ged_obol_depth_consider(view_ctx, record.labels[i].target, mode,
				    have_depth, calc_depth);
    }
    for (size_t i = 0; i < record.axesCenters.size(); i++)
	ged_obol_depth_consider(view_ctx, record.axesCenters[i], mode,
				have_depth, calc_depth);

    if (!have_depth)
	return 0;

    *depth = calc_depth;
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_depth(
    void *view_ctx,
    const char *name,
    int mode,
    fastf_t *depth)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    BRLObolFeatureRecord record;
    if (!handle.isValid() || !controller->features().record(handle, record))
	return 0;
    return ged_obol_feature_record_depth(view_ctx, record, mode, depth);
}

struct ged_obol_feature_depth_visit {
    void *view_ctx;
    int mode;
    ged_draw_view_feature_depth_cb cb;
    void *data;
    int count;
};

static int
ged_obol_feature_depth_visit_cb(const BRLObolFeatureRecord &record,
				void *userData)
{
    struct ged_obol_feature_depth_visit *ctx =
	(struct ged_obol_feature_depth_visit *)userData;
    if (!ctx || !ctx->cb)
	return 0;

    fastf_t depth = 0.0;
    if (!ged_obol_feature_record_depth(ctx->view_ctx, record, ctx->mode,
				       &depth))
	return 1;

    ctx->count++;
    return ctx->cb(depth, ctx->data);
}

extern "C" int
ged_draw_obol_view_context_feature_depth_foreach(
    void *view_ctx,
    int mode,
    ged_draw_view_feature_depth_cb cb,
    void *data)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !cb)
	return 0;

    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    struct ged_obol_feature_depth_visit ctx;
    ctx.view_ctx = view_ctx;
    ctx.mode = mode;
    ctx.cb = cb;
    ctx.data = data;
    ctx.count = 0;
    controller->features().visitRecords(ged_obol_feature_depth_visit_cb,
					&ctx, BRLOBOL_FEATURE_SCOPE_ALL, &owner);
    return ctx.count;
}

struct ged_obol_feature_records_visit {
    void *view_ctx;
    unsigned int query_flags;
    const char *glob;
    ged_draw_obol_view_feature_record_cb cb;
    void *userdata;
    int count;
};

static BRLObolFeatureStyle
ged_obol_line_layer_effective_style(const BRLObolFeatureStyle &parent,
				    const BRLObolFeatureStyle &layer)
{
    BRLObolFeatureStyle out = parent;
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
ged_obol_feature_record_visible(const BRLObolFeatureStyle &parent,
				const BRLObolFeatureStyle *child)
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
			     const BRLObolFeatureRecord &record,
			     const SbString &name,
			     enum ged_draw_obol_view_feature_kind kind,
			     const BRLObolFeatureStyle &style,
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
    out.local = record.scope == BRLObolFeatureScope::Local ? 1 : 0;
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
ged_obol_feature_record_visit_cb(const BRLObolFeatureRecord &record,
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

    if (record.kind == BRLObolFeatureKind::LineLayer &&
	!record.layers.empty()) {
	for (size_t i = 0; i < record.layers.size(); i++) {
	    const BRLObolLineLayer &layer = record.layers[i];
	    BRLObolFeatureStyle style =
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
    else if (record.kind == BRLObolFeatureKind::CustomNode &&
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
    void *view_ctx,
    unsigned int query_flags,
    const char *glob,
    ged_draw_obol_view_feature_record_cb cb,
    void *userdata)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !cb)
	return 0;

    const int wants_db =
	(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS) ? 1 : 0;
    const int wants_view =
	(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS) ? 1 : 0;
    if (wants_db && !wants_view)
	return 0;

    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    unsigned int scope_mask = BRLOBOL_FEATURE_SCOPE_ALL;
    if (query_flags & GED_DRAW_VIEW_EXPORT_QUERY_LOCAL_ONLY)
	scope_mask = BRLOBOL_FEATURE_SCOPE_LOCAL;

    struct ged_obol_feature_records_visit ctx;
    ctx.view_ctx = view_ctx;
    ctx.query_flags = query_flags;
    ctx.glob = glob;
    ctx.cb = cb;
    ctx.userdata = userdata;
    ctx.count = 0;
    controller->features().visitRecords(ged_obol_feature_record_visit_cb,
					&ctx, scope_mask, &owner);
    return ctx.count;
}

extern "C" int
ged_draw_obol_view_context_indexed_face_set_replace(
    void *view_ctx,
    const char *name,
    int local,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count,
    const struct ged_draw_view_feature_style *style)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name || !points || !point_count || !indices ||
	!index_count)
	return 0;

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureHandle handle =
	ged_obol_publish_indexed_face_set(controller, view_ctx, name, local,
					  ged_obol_points_from_ged(points, point_count),
					  ged_obol_vectors_from_ged(normals, normal_count),
					  ged_obol_indices_from_ged(indices, index_count),
					  style ? &obol_style : NULL);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_lines_replace(
    void *view_ctx,
    const char *name,
    int local,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_draw_view_feature_style *style)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;
    if (!points || !point_count)
	return ged_obol_remove_feature(controller, view_ctx, name,
				       local ? 1 : 0) ? 1 : 1;

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureHandle handle =
	ged_obol_publish_line_set(controller, view_ctx, name, local,
				  ged_obol_points_from_ged(points, point_count),
				  ged_obol_commands_from_ged(cmds, point_count),
				  style ? &obol_style : NULL);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_tcl_polygons_replace(
    void *view_ctx,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_draw_view_feature_style *style)
{
    if (!points || !cmds || !point_count) {
	BRLObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	return (controller && name) ?
	       (ged_obol_remove_feature(controller, view_ctx, name, 1) ? 1 : 1) :
	       0;
    }
    int ret = ged_draw_obol_view_context_lines_replace(view_ctx, name, 1,
	      points, cmds, point_count, style);
    if (!ret)
	return 0;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BRLObolOverlayClass::PolygonEdit,
						 BRLObolOverlayLifecycle::PerTool,
						 BRLObolOverlayOrder::PostTransparent,
						 name));
}

static const uint32_t GED_OBOL_FEATURE_TOKEN_MAGIC = 0x474f4654U;

struct ged_obol_feature_ref_token {
    uint32_t magic;
    BRLObolViewController *controller;
    uint64_t id;
    std::string name;
};

static std::vector<ged_obol_feature_ref_token *> ged_obol_feature_ref_tokens;

static ged_obol_feature_ref_token *
ged_obol_feature_token_find(BRLObolViewController *controller, uint64_t id)
{
    if (!controller || !id)
	return NULL;

    for (size_t i = 0; i < ged_obol_feature_ref_tokens.size(); i++) {
	ged_obol_feature_ref_token *token = ged_obol_feature_ref_tokens[i];
	if (token && token->controller == controller && token->id == id)
	    return token;
    }
    return NULL;
}

static ged_obol_feature_ref_token *
ged_obol_feature_token_get(BRLObolViewController *controller,
			   BRLObolFeatureHandle handle,
			   const char *name = NULL)
{
    if (!controller || !handle.isValid())
	return NULL;

    ged_obol_feature_ref_token *token =
	ged_obol_feature_token_find(controller, handle.id);
    if (!token) {
	token = new ged_obol_feature_ref_token;
	token->magic = GED_OBOL_FEATURE_TOKEN_MAGIC;
	token->controller = controller;
	token->id = handle.id;
	ged_obol_feature_ref_tokens.push_back(token);
    }
    if (name)
	token->name = name;
    return token;
}

static ged_obol_feature_ref_token *
ged_obol_feature_token_from_rt_ref(rt_view_feature_ref ref)
{
    if (!ref.token)
	return NULL;

    for (size_t i = 0; i < ged_obol_feature_ref_tokens.size(); i++) {
	ged_obol_feature_ref_token *token = ged_obol_feature_ref_tokens[i];
	if (token && reinterpret_cast<uintptr_t>(token) == ref.token &&
	    token->magic == GED_OBOL_FEATURE_TOKEN_MAGIC)
	    return token;
    }
    return NULL;
}

static BRLObolFeatureHandle
ged_obol_feature_handle_from_rt_ref(rt_view_feature_ref ref)
{
    ged_obol_feature_ref_token *token =
	ged_obol_feature_token_from_rt_ref(ref);
    return token ? BRLObolFeatureHandle(token->id, ref.revision) :
	   BRLObolFeatureHandle();
}

static rt_view_feature_ref
ged_obol_rt_feature_ref(BRLObolViewController *controller,
			BRLObolFeatureHandle handle,
			const char *name = NULL)
{
    ged_obol_feature_ref_token *token =
	ged_obol_feature_token_get(controller, handle, name);
    rt_view_feature_ref ref = RT_VIEW_FEATURE_REF_NULL_INIT;
    if (!token)
	return ref;
    ref.token = reinterpret_cast<uintptr_t>(token);
    ref.revision = handle.revision;
    return ref;
}

static BRLObolOverlayInfo
ged_obol_rt_edit_overlay_info(void *view_ctx,
			      const void *owner,
			      const char *source_path,
			      int sort_order)
{
    BRLObolOverlayInfo overlay = ged_obol_model_overlay_info(view_ctx,
				 BRLObolOverlayClass::EditHandle,
				 BRLObolOverlayLifecycle::PerTool,
				 BRLObolOverlayOrder::PostTransparent,
				 source_path);
    overlay.ownerToken = owner ? owner : view_ctx;
    overlay.sortOrder = sort_order;
    return overlay;
}

static std::vector<BRLObolLabel>
ged_obol_labels_from_rt(const struct rt_view_feature_label *labels,
			size_t label_count)
{
    std::vector<BRLObolLabel> out;
    if (!labels || !label_count)
	return out;

    out.reserve(label_count);
    for (size_t i = 0; i < label_count; i++) {
	BRLObolLabel label;
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

static BRLObolEditPreviewCallbacks
ged_obol_edit_preview_callbacks_from_rt(
    const struct rt_view_edit_preview_callbacks *callbacks,
    void *preview_ctx)
{
    BRLObolEditPreviewCallbacks out;
    out.previewContext = preview_ctx;
    if (callbacks) {
	out.revisionCallback = callbacks->revision_cb;
	out.updateCallback = callbacks->update_cb;
	out.pickCallback = callbacks->pick_cb;
    }
    return out;
}

static int
ged_draw_rt_obol_feature_owns_ref(rt_view_feature_ref ref,
				  void *UNUSED(data))
{
    return ged_obol_feature_token_from_rt_ref(ref) ? 1 : 0;
}

static int
ged_draw_rt_obol_edit_preview_publish_event(
    void *UNUSED(view_ctx),
    rt_view_feature_ref feature,
    enum rt_view_edit_preview_event UNUSED(event),
    const char *UNUSED(source_path),
    void *UNUSED(data))
{
    return ged_obol_feature_token_from_rt_ref(feature) ? 1 : 0;
}

static rt_view_feature_ref
ged_draw_rt_obol_feature_overlay_ensure(
    void *view_ctx,
    const char *name,
    const void *owner,
    void *preview_ctx,
    const struct rt_view_edit_preview_callbacks *callbacks,
    const char *source_path,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller || !name || !name[0])
	return RT_VIEW_FEATURE_REF_NULL;

    BRLObolFeatureOwner feature_owner = ged_obol_feature_owner(view_ctx, 1);
    BRLObolFeatureHandle handle = controller->features().findOwned(name,
				  BRLOBOL_FEATURE_SCOPE_LOCAL, &feature_owner);

    BRLObolFeatureRecord record;
    const int needs_publish = !handle.isValid() ||
			      !controller->features().record(handle, record) ||
			      record.kind != BRLObolFeatureKind::EditPreview;
    if (needs_publish) {
	std::vector<SbVec3f> points;
	std::vector<int32_t> commands;
	BRLObolEditPreviewCallbacks obol_callbacks =
	    ged_obol_edit_preview_callbacks_from_rt(callbacks, preview_ctx);
	handle = controller->features().publishEditPreview(name,
		 source_path && source_path[0] ? source_path : name,
		 name,
		 "edit-handle",
		 points,
		 commands,
		 0,
		 0,
		 callbacks ? &obol_callbacks : NULL,
		 &feature_owner);
    }

    if (!handle.isValid())
	return RT_VIEW_FEATURE_REF_NULL;

    (void)controller->features().setOverlayInfo(handle,
	    ged_obol_rt_edit_overlay_info(view_ctx, owner, source_path, 0));
    return ged_obol_rt_feature_ref(controller, handle, name);
}

static rt_view_feature_ref
ged_draw_rt_obol_feature_label_ensure(
    void *view_ctx,
    const char *name,
    const void *owner,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller || !name || !name[0])
	return RT_VIEW_FEATURE_REF_NULL;

    BRLObolFeatureOwner feature_owner = ged_obol_feature_owner(view_ctx, 1);
    BRLObolFeatureHandle handle = controller->features().findOwned(name,
				  BRLOBOL_FEATURE_SCOPE_LOCAL, &feature_owner);

    BRLObolFeatureRecord record;
    if (!handle.isValid() ||
	!controller->features().record(handle, record) ||
	record.kind != BRLObolFeatureKind::Labels) {
	std::vector<BRLObolLabel> labels;
	handle = controller->features().publishLabels(name,
		 BRLObolFeatureScope::Local, labels, NULL, &feature_owner);
    }

    if (!handle.isValid())
	return RT_VIEW_FEATURE_REF_NULL;

    (void)controller->features().setOverlayInfo(handle,
	    ged_obol_rt_edit_overlay_info(view_ctx, owner, name, 1));
    return ged_obol_rt_feature_ref(controller, handle, name);
}

static int
ged_draw_rt_obol_feature_remove(void *view_ctx,
				const char *name,
				void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    return ged_obol_remove_feature(controller, view_ctx, name, -1);
}

static int
ged_draw_rt_obol_feature_set_context(rt_view_feature_ref ref,
				     void *UNUSED(ctx),
				     void *UNUSED(data))
{
    return ged_obol_feature_token_from_rt_ref(ref) ? 1 : 0;
}

static int
ged_draw_rt_obol_feature_set_visible(rt_view_feature_ref ref,
				     int visible,
				     void *UNUSED(data))
{
    ged_obol_feature_ref_token *token =
	ged_obol_feature_token_from_rt_ref(ref);
    return (token && token->controller) ?
	   (token->controller->features().setVisible(
		ged_obol_feature_handle_from_rt_ref(ref),
		visible ? TRUE : FALSE) ? 1 : 0) : 0;
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

static int
ged_draw_rt_obol_feature_set_color(rt_view_feature_ref ref,
				   int r,
				   int g,
				   int b,
				   void *UNUSED(data))
{
    ged_obol_feature_ref_token *token =
	ged_obol_feature_token_from_rt_ref(ref);
    if (!token || !token->controller)
	return 0;

    const unsigned char rgb[3] = {
	ged_obol_rgb_byte_from_int(r),
	ged_obol_rgb_byte_from_int(g),
	ged_obol_rgb_byte_from_int(b)
    };
    return token->controller->features().setColor(
	       ged_obol_feature_handle_from_rt_ref(ref),
	       ged_obol_color_from_rgb(rgb)) ? 1 : 0;
}

static int
ged_draw_rt_obol_feature_touch(rt_view_feature_ref ref,
			       void *UNUSED(data))
{
    ged_obol_feature_ref_token *token =
	ged_obol_feature_token_from_rt_ref(ref);
    if (!token || !token->controller)
	return 0;

    const int ret = token->controller->features().updateEditPreview(
			ged_obol_feature_handle_from_rt_ref(ref));
    return ret >= 0 ? ret : 1;
}

static int
ged_draw_rt_obol_feature_labels_replace(
    rt_view_feature_ref ref,
    const struct rt_view_feature_label *labels,
    size_t label_count,
    void *UNUSED(data))
{
    if (label_count && !labels)
	return 0;

    ged_obol_feature_ref_token *token =
	ged_obol_feature_token_from_rt_ref(ref);
    if (!token || !token->controller)
	return 0;

    BRLObolFeatureRecord record;
    BRLObolFeatureHandle handle = ged_obol_feature_handle_from_rt_ref(ref);
    if (!token->controller->features().record(handle, record))
	return 0;

    std::vector<BRLObolLabel> obol_labels =
	ged_obol_labels_from_rt(labels, label_count);
    if (record.kind == BRLObolFeatureKind::Labels ||
	record.kind == BRLObolFeatureKind::HudLabel)
	return token->controller->features().replaceLabels(handle,
		obol_labels) ? 1 : 0;

    BRLObolFeatureOwner owner = record.owner;
    BRLObolFeatureHandle labels_handle =
	token->controller->features().publishLabels(record.name,
	    record.scope, obol_labels, &record.style,
	    record.scope == BRLObolFeatureScope::Local ? &owner : NULL);
    return labels_handle.isValid() ? 1 : 0;
}

static int
ged_draw_rt_obol_feature_points_replace(
    rt_view_feature_ref ref,
    enum rt_view_feature_family family,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    void *UNUSED(data))
{
    if (point_count && !points)
	return 0;

    ged_obol_feature_ref_token *token =
	ged_obol_feature_token_from_rt_ref(ref);
    if (!token || !token->controller)
	return 0;

    BRLObolFeatureRecord record;
    BRLObolFeatureHandle handle = ged_obol_feature_handle_from_rt_ref(ref);
    if (!token->controller->features().record(handle, record))
	return 0;

    std::vector<SbVec3f> obol_points =
	ged_obol_points_from_ged(points, point_count);
    std::vector<int32_t> obol_commands =
	ged_obol_commands_from_ged(cmds, point_count);
    if (family == RT_VIEW_FEATURE_TRANSIENT_PREVIEW ||
	record.kind == BRLObolFeatureKind::EditPreview) {
	if (record.kind == BRLObolFeatureKind::EditPreview)
	    return token->controller->features().replaceEditPreviewGeometry(
		       handle,
		       record.name,
		       obol_points,
		       obol_commands) ? 1 : 0;

	BRLObolFeatureOwner owner = record.owner;
	BRLObolFeatureHandle preview_handle =
	    token->controller->features().publishEditPreview(record.name,
		record.name,
		record.name,
		"edit-handle",
		obol_points,
		obol_commands,
		0,
		0,
		NULL,
		record.scope == BRLObolFeatureScope::Local ? &owner : NULL);
	return preview_handle.isValid() ? 1 : 0;
    }

    BRLObolFeatureOwner owner = record.owner;
    BRLObolFeatureHandle line_handle =
	token->controller->features().publishLineSet(record.name,
	    record.scope, obol_points, obol_commands, &record.style,
	    record.scope == BRLObolFeatureScope::Local ? &owner : NULL);
    return line_handle.isValid() ? 1 : 0;
}

static int
ged_draw_rt_obol_feature_clear_geometry(rt_view_feature_ref ref,
					void *UNUSED(data))
{
    ged_obol_feature_ref_token *token =
	ged_obol_feature_token_from_rt_ref(ref);
    return (token && token->controller) ?
	   (token->controller->features().clearGeometry(
		ged_obol_feature_handle_from_rt_ref(ref)) ? 1 : 0) : 0;
}

extern "C" int
ged_draw_view_context_obol_feature_adapter_attach(
    struct ged *gedp,
    void *view_ctx)
{
    if (!gedp || !view_ctx)
	return 0;

    struct rt_view_context_feature_adapter adapter;
    memset(&adapter, 0, sizeof(adapter));
    adapter.owns_ref = ged_draw_rt_obol_feature_owns_ref;
    adapter.edit_preview_publish_event =
	ged_draw_rt_obol_edit_preview_publish_event;
    adapter.overlay_ensure = ged_draw_rt_obol_feature_overlay_ensure;
    adapter.label_ensure = ged_draw_rt_obol_feature_label_ensure;
    adapter.remove = ged_draw_rt_obol_feature_remove;
    adapter.set_context = ged_draw_rt_obol_feature_set_context;
    adapter.set_visible = ged_draw_rt_obol_feature_set_visible;
    adapter.set_color = ged_draw_rt_obol_feature_set_color;
    adapter.touch = ged_draw_rt_obol_feature_touch;
    adapter.labels_replace = ged_draw_rt_obol_feature_labels_replace;
    adapter.points_replace = ged_draw_rt_obol_feature_points_replace;
    adapter.clear_geometry = ged_draw_rt_obol_feature_clear_geometry;
    adapter.data = gedp;
    return rt_view_context_feature_adapter_set(view_ctx, &adapter);
}

static const uint32_t GED_OBOL_POLYGON_TOKEN_MAGIC = 0x474f504cU;

struct ged_obol_polygon_ref_token {
    uint32_t magic;
    BRLObolViewController *controller;
    uint64_t id;
    std::string name;
};

static std::vector<ged_obol_polygon_ref_token *> ged_obol_polygon_ref_tokens;
static uint64_t ged_obol_polygon_auto_id = 1;

static BRLObolPolygonType
ged_obol_polygon_type_from_ged(int type)
{
    switch (type) {
	case GED_DRAW_VIEW_POLYGON_CIRCLE:
	    return BRLObolPolygonType::Circle;
	case GED_DRAW_VIEW_POLYGON_ELLIPSE:
	    return BRLObolPolygonType::Ellipse;
	case GED_DRAW_VIEW_POLYGON_RECTANGLE:
	    return BRLObolPolygonType::Rectangle;
	case GED_DRAW_VIEW_POLYGON_SQUARE:
	    return BRLObolPolygonType::Square;
	default:
	    return BRLObolPolygonType::General;
    }
}

static BRLObolPolygonUpdate
ged_obol_polygon_update_from_ged(int op)
{
    switch (op) {
	case GED_DRAW_VIEW_POLYGON_UPDATE_PROPS_ONLY:
	    return BRLObolPolygonUpdate::PropsOnly;
	case GED_DRAW_VIEW_POLYGON_UPDATE_PT_SELECT:
	    return BRLObolPolygonUpdate::PointSelect;
	case GED_DRAW_VIEW_POLYGON_UPDATE_PT_SELECT_CLEAR:
	    return BRLObolPolygonUpdate::PointSelectClear;
	case GED_DRAW_VIEW_POLYGON_UPDATE_PT_MOVE:
	    return BRLObolPolygonUpdate::PointMove;
	case GED_DRAW_VIEW_POLYGON_UPDATE_PT_APPEND:
	    return BRLObolPolygonUpdate::PointAppend;
	default:
	    return BRLObolPolygonUpdate::Default;
    }
}

static BRLObolFeatureScope
ged_obol_polygon_scope(int local)
{
    return local ? BRLObolFeatureScope::Local : BRLObolFeatureScope::Shared;
}

static ged_obol_polygon_ref_token *
ged_obol_polygon_token_find(BRLObolViewController *controller,
			    uint64_t id)
{
    if (!controller || !id)
	return NULL;
    for (size_t i = 0; i < ged_obol_polygon_ref_tokens.size(); i++) {
	ged_obol_polygon_ref_token *token = ged_obol_polygon_ref_tokens[i];
	if (token && token->controller == controller && token->id == id)
	    return token;
    }
    return NULL;
}

static ged_obol_polygon_ref_token *
ged_obol_polygon_token_get(BRLObolViewController *controller,
			   BRLObolPolygonHandle handle)
{
    if (!controller || !handle.isValid())
	return NULL;

    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_find(controller, handle.id);
    if (token)
	return token;

    token = new ged_obol_polygon_ref_token;
    token->magic = GED_OBOL_POLYGON_TOKEN_MAGIC;
    token->controller = controller;
    token->id = handle.id;
    ged_obol_polygon_ref_tokens.push_back(token);
    return token;
}

static ged_obol_polygon_ref_token *
ged_obol_polygon_token_from_rt_ref(rt_view_polygon_ref ref)
{
    if (!ref.token)
	return NULL;

    for (size_t i = 0; i < ged_obol_polygon_ref_tokens.size(); i++) {
	ged_obol_polygon_ref_token *token = ged_obol_polygon_ref_tokens[i];
	if (token && reinterpret_cast<uintptr_t>(token) == ref.token &&
	    token->magic == GED_OBOL_POLYGON_TOKEN_MAGIC)
	    return token;
    }

    return NULL;
}

static BRLObolPolygonHandle
ged_obol_polygon_handle_from_rt_ref(rt_view_polygon_ref ref)
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    return token ? BRLObolPolygonHandle(token->id, ref.revision) :
	   BRLObolPolygonHandle();
}

static rt_view_polygon_ref
ged_obol_rt_polygon_ref(BRLObolViewController *controller,
			BRLObolPolygonHandle handle)
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_get(controller, handle);
    rt_view_polygon_ref ref = RT_VIEW_POLYGON_REF_NULL_INIT;
    if (!token)
	return ref;
    ref.token = reinterpret_cast<uintptr_t>(token);
    ref.revision = handle.revision;
    return ref;
}

static ged_draw_view_polygon_ref
ged_obol_ged_polygon_ref(BRLObolViewController *controller,
			 BRLObolPolygonHandle handle)
{
    rt_view_polygon_ref rt_ref = ged_obol_rt_polygon_ref(controller, handle);
    ged_draw_view_polygon_ref ged_ref = GED_DRAW_VIEW_POLYGON_REF_NULL_INIT;
    ged_ref.token = rt_ref.token;
    ged_ref.revision = rt_ref.revision;
    return ged_ref;
}

static rt_view_polygon_ref
ged_obol_rt_polygon_ref_from_ged(ged_draw_view_polygon_ref ref)
{
    rt_view_polygon_ref rt_ref = RT_VIEW_POLYGON_REF_NULL_INIT;
    rt_ref.token = ref.token;
    rt_ref.revision = ref.revision;
    return rt_ref;
}

static void
ged_obol_polygon_point_from_ged(SbVec3f &dst, const point_t src)
{
    dst = SbVec3f(static_cast<float>(src[X]),
		  static_cast<float>(src[Y]),
		  static_cast<float>(src[Z]));
}

static void
ged_obol_polygon_project_point(point_t dst, void *view_ctx, const point_t src,
			       plane_t *view_plane)
{
    VMOVE(dst, src);
    if (!view_plane)
	return;
    HSET(*view_plane, 0.0, 0.0, 1.0, src[Z]);

    if (rt_view_context_plane_get(view_plane, view_ctx) != 0)
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
ged_obol_polygon_record_to_rt(
    BRLObolViewController *controller,
    rt_view_polygon_ref ref,
    const BRLObolPolygonRecord &src,
    struct rt_view_polygon_record *dst)
{
    if (!controller || !dst)
	return 0;

    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_get(controller, src.handle);
    if (!token)
	return 0;

    token->name = src.name.getString() ? src.name.getString() : "";
    memset(dst, 0, sizeof(*dst));
    dst->ref = ref;
    dst->ref.revision = src.handle.revision;
    dst->name = token->name.c_str();
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

static BRLObolPolygonHandle
ged_obol_polygon_create_named(
    BRLObolViewController *controller,
    void *view_ctx,
    const char *name,
    int local,
    int type,
    const point_t input_point)
{
    if (!controller || !name || !name[0] || !input_point)
	return BRLObolPolygonHandle();

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

static std::string
ged_obol_polygon_auto_name(BRLObolViewController *controller)
{
    if (!controller)
	return std::string();

    for (;;) {
	char buf[64];
	snprintf(buf, sizeof(buf), "__rt_polygon_%llu",
		 (unsigned long long)ged_obol_polygon_auto_id++);
	if (!controller->polygons().find(buf).isValid())
	    return std::string(buf);
    }
}

extern "C" ged_draw_view_polygon_ref
ged_draw_obol_view_context_polygon_find(void *view_ctx, const char *name)
{
    if (!name)
	return GED_DRAW_VIEW_POLYGON_REF_NULL;

    BRLObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (local_controller) {
	BRLObolPolygonHandle local =
	    local_controller->polygons().find(name,
					      BRLOBOL_FEATURE_SCOPE_LOCAL);
	if (local.isValid())
	    return ged_obol_ged_polygon_ref(local_controller, local);
    }

    BRLObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if (!shared_controller)
	return GED_DRAW_VIEW_POLYGON_REF_NULL;
    return ged_obol_ged_polygon_ref(shared_controller,
				    shared_controller->polygons().find(name,
					    BRLOBOL_FEATURE_SCOPE_SHARED));
}

extern "C" ged_draw_view_polygon_ref
ged_draw_obol_view_context_polygon_find_scoped(
    void *view_ctx,
    const char *name,
    int local_only)
{
    if (!name)
	return GED_DRAW_VIEW_POLYGON_REF_NULL;

    if (local_only) {
	BRLObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	if (!controller)
	    return GED_DRAW_VIEW_POLYGON_REF_NULL;
	return ged_obol_ged_polygon_ref(controller,
					controller->polygons().find(name, BRLOBOL_FEATURE_SCOPE_LOCAL));
    }

    return ged_draw_obol_view_context_polygon_find(view_ctx, name);
}

extern "C" ged_draw_view_polygon_ref
ged_draw_obol_view_context_polygon_create(
    void *view_ctx,
    const char *name,
    int local,
    int type,
    const point_t screen_point)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    return ged_obol_ged_polygon_ref(controller,
				    ged_obol_polygon_create_named(controller, view_ctx, name, local,
					    type, screen_point));
}

extern "C" ged_draw_view_polygon_ref
ged_draw_obol_view_context_polygon_import_sketch(
    const char *name,
    struct db_i *dbip,
    struct directory *dp,
    void *view_ctx,
    int local)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return GED_DRAW_VIEW_POLYGON_REF_NULL;
    return ged_obol_ged_polygon_ref(controller,
				    controller->polygons().importSketch(name,
					    ged_obol_polygon_scope(local), dbip, dp));
}

extern "C" int
ged_draw_obol_view_context_polygon_set_current(
    ged_draw_view_polygon_ref ref,
    long contour_i,
    long point_i)
{
    rt_view_polygon_ref rt_ref = ged_obol_rt_polygon_ref_from_ged(ref);
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(rt_ref);
    if (!token || !token->controller)
	return 0;
    return token->controller->polygons().setCurrent(
	       ged_obol_polygon_handle_from_rt_ref(rt_ref),
	       contour_i, point_i) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_polygon_set_contour_open(
    ged_draw_view_polygon_ref ref,
    long contour_i,
    int open)
{
    rt_view_polygon_ref rt_ref = ged_obol_rt_polygon_ref_from_ged(ref);
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(rt_ref);
    if (!token || !token->controller)
	return 0;
    return token->controller->polygons().setContourOpen(
	       ged_obol_polygon_handle_from_rt_ref(rt_ref),
	       contour_i, open ? TRUE : FALSE) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_polygon_area(
    ged_draw_view_polygon_ref ref,
    void *view_ctx,
    fastf_t *area)
{
    if (area)
	*area = 0.0;
    if (!area || !view_ctx)
	return 0;

    rt_view_polygon_ref rt_ref = ged_obol_rt_polygon_ref_from_ged(ref);
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(rt_ref);
    if (!token || !token->controller)
	return 0;

    *area = token->controller->polygons().area(
		ged_obol_polygon_handle_from_rt_ref(rt_ref),
		ged_view_context_scale_get(view_ctx));
    return 1;
}

extern "C" int
ged_draw_obol_view_context_polygon_overlap(
    ged_draw_view_polygon_ref ref,
    void *view_ctx,
    const char *other_name,
    const struct bn_tol *tol,
    int *overlap)
{
    if (overlap)
	*overlap = 0;
    if (!view_ctx || !other_name || !tol || !overlap)
	return 0;

    rt_view_polygon_ref rt_ref = ged_obol_rt_polygon_ref_from_ged(ref);
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(rt_ref);
    if (!token || !token->controller)
	return 0;

    BRLObolPolygonHandle other =
	token->controller->polygons().find(other_name);
    if (!other.isValid())
	return 0;

    *overlap = token->controller->polygons().overlaps(
		   ged_obol_polygon_handle_from_rt_ref(rt_ref),
		   other,
		   *tol,
		   ged_view_context_scale_get(view_ctx)) ? 1 : 0;
    return 1;
}

static int
ged_draw_rt_obol_polygon_owns_ref(rt_view_polygon_ref ref, void *UNUSED(data))
{
    return ged_obol_polygon_token_from_rt_ref(ref) ? 1 : 0;
}

static int
ged_draw_rt_obol_polygon_record_get(
    rt_view_polygon_ref ref,
    struct rt_view_polygon_record *record,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller || !record)
	return 0;

    BRLObolPolygonRecord obol_record;
    if (!token->controller->polygons().record(
	    ged_obol_polygon_handle_from_rt_ref(ref), obol_record))
	return 0;
    return ged_obol_polygon_record_to_rt(token->controller, ref,
					 obol_record, record);
}

static rt_view_polygon_ref
ged_draw_rt_obol_polygon_create(
    void *view_ctx,
    int type,
    point_t *fp,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!fp)
	return RT_VIEW_POLYGON_REF_NULL;
    std::string name = ged_obol_polygon_auto_name(controller);
    return ged_obol_rt_polygon_ref(controller,
				   ged_obol_polygon_create_named(controller, view_ctx, name.c_str(),
					   0, type, *fp));
}

static rt_view_polygon_ref
ged_draw_rt_obol_polygon_select(
    void *view_ctx,
    point_t *cp,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !cp)
	return RT_VIEW_POLYGON_REF_NULL;
    return ged_obol_rt_polygon_ref(controller,
				   controller->polygons().selectAtModelPoint(SbVec3f(
					   static_cast<float>((*cp)[X]),
					   static_cast<float>((*cp)[Y]),
					   static_cast<float>((*cp)[Z]))));
}

static rt_view_polygon_ref
ged_draw_rt_obol_polygon_find(
    void *view_ctx,
    const char *name,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name)
	return RT_VIEW_POLYGON_REF_NULL;
    return ged_obol_rt_polygon_ref(controller,
				   controller->polygons().find(name));
}

static rt_view_polygon_ref
ged_draw_rt_obol_polygon_dup(
    void *view_ctx,
    const char *name,
    const char *new_name,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !new_name)
	return RT_VIEW_POLYGON_REF_NULL;
    BRLObolPolygonHandle src = controller->polygons().find(name);
    return ged_obol_rt_polygon_ref(controller,
				   controller->polygons().duplicate(src, new_name));
}

struct ged_obol_polygon_visit_context {
    BRLObolViewController *controller;
    rt_view_polygon_record_callback_t callback;
    void *data;
};

static int
ged_draw_rt_obol_polygon_visit_cb(const BRLObolPolygonRecord &obol_record,
				  void *data)
{
    ged_obol_polygon_visit_context *ctx =
	static_cast<ged_obol_polygon_visit_context *>(data);
    if (!ctx || !ctx->controller || !ctx->callback)
	return 0;

    rt_view_polygon_ref ref =
	ged_obol_rt_polygon_ref(ctx->controller, obol_record.handle);
    struct rt_view_polygon_record record;
    if (!ged_obol_polygon_record_to_rt(ctx->controller, ref, obol_record,
				       &record))
	return 1;
    return ctx->callback(ref, &record, ctx->data);
}

static void
ged_draw_rt_obol_polygon_visit_records(
    void *view_ctx,
    rt_view_polygon_record_callback_t callback,
    void *callback_data,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !callback)
	return;

    ged_obol_polygon_visit_context ctx;
    ctx.controller = controller;
    ctx.callback = callback;
    ctx.data = callback_data;
    controller->polygons().visitRecords(
	ged_draw_rt_obol_polygon_visit_cb, &ctx);
}

static size_t
ged_draw_rt_obol_polygon_snap_count(
    void *view_ctx,
    rt_view_polygon_ref exclude,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;
    BRLObolPolygonHandle exclude_handle;
    if (ged_obol_polygon_token_from_rt_ref(exclude))
	exclude_handle = ged_obol_polygon_handle_from_rt_ref(exclude);
    return controller->polygons().snapCount(exclude_handle);
}

static int
ged_draw_rt_obol_polygon_clear_point_selection(
    void *view_ctx,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    return controller ?
	   (controller->polygons().clearAllPointSelections() ? 1 : 0) : 0;
}

static int
ged_draw_rt_obol_polygon_update(
    rt_view_polygon_ref ref,
    void *UNUSED(view_ctx),
    int utype,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller)
	return 0;
    return token->controller->polygons().update(
	       ged_obol_polygon_handle_from_rt_ref(ref),
	       ged_obol_polygon_update_from_ged(utype)) ? 1 : 0;
}

static int
ged_draw_rt_obol_polygon_update_screen_pt(
    rt_view_polygon_ref ref,
    void *view_ctx,
    int x,
    int y,
    int utype,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller || !view_ctx)
	return 0;

    point_t model_point = VINIT_ZERO;
    if (!ged_view_context_screen_point(model_point, view_ctx,
				       (fastf_t)x, (fastf_t)y))
	return 0;

    return token->controller->polygons().updateModelPoint(
	       ged_obol_polygon_handle_from_rt_ref(ref),
	       SbVec3f(static_cast<float>(model_point[X]),
		       static_cast<float>(model_point[Y]),
		       static_cast<float>(model_point[Z])),
	       ged_obol_polygon_update_from_ged(utype)) ? 1 : 0;
}

static int
ged_draw_rt_obol_polygon_move(
    rt_view_polygon_ref ref,
    point_t *current_point,
    point_t *previous_point,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller || !current_point || !previous_point)
	return 0;
    return token->controller->polygons().move(
	       ged_obol_polygon_handle_from_rt_ref(ref),
	       SbVec3f(static_cast<float>((*current_point)[X]),
		       static_cast<float>((*current_point)[Y]),
		       static_cast<float>((*current_point)[Z])),
	       SbVec3f(static_cast<float>((*previous_point)[X]),
		       static_cast<float>((*previous_point)[Y]),
		       static_cast<float>((*previous_point)[Z]))) ? 1 : 0;
}

static int
ged_draw_rt_obol_polygon_set_name(
    rt_view_polygon_ref ref,
    const char *name,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller || !name)
	return 0;
    return token->controller->polygons().rename(
	       ged_obol_polygon_handle_from_rt_ref(ref), name) ? 1 : 0;
}

static int
ged_draw_rt_obol_polygon_set_context(
    rt_view_polygon_ref ref,
    void *UNUSED(ctx),
    void *UNUSED(data))
{
    return ged_obol_polygon_token_from_rt_ref(ref) ? 1 : 0;
}

static int
ged_draw_rt_obol_polygon_set_visual(
    rt_view_polygon_ref ref,
    const struct bu_color *edge_color,
    const struct bu_color *fill_color,
    fastf_t fill_slope_x,
    fastf_t fill_slope_y,
    fastf_t fill_density,
    fastf_t vZ,
    int fill_flag,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller)
	return 0;

    BRLObolPolygonVisual visual;
    (void)token->controller->polygons().visual(
	ged_obol_polygon_handle_from_rt_ref(ref), visual);
    if (edge_color)
	visual.edgeColor = ged_obol_polygon_color_from_bu(edge_color);
    if (fill_color)
	visual.fillColor = ged_obol_polygon_color_from_bu(fill_color);
    visual.fillSlope = SbVec2f(static_cast<float>(fill_slope_x),
			       static_cast<float>(fill_slope_y));
    visual.fillSpacing = static_cast<float>(fill_density);
    visual.viewZ = static_cast<float>(vZ);
    if (fill_flag)
	visual.fillFlags |= BRLOBOL_POLYGON_FILL_HATCH;
    else
	visual.fillFlags &= ~BRLOBOL_POLYGON_FILL_HATCH;
    visual.fill = (visual.fillFlags & BRLOBOL_POLYGON_FILL_HATCH) ?
		  TRUE : FALSE;
    return token->controller->polygons().setVisual(
	       ged_obol_polygon_handle_from_rt_ref(ref), visual) ? 1 : 0;
}

static int
ged_draw_rt_obol_polygon_set_open(
    rt_view_polygon_ref ref,
    int open,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller)
	return 0;
    return token->controller->polygons().setAllContoursOpen(
	       ged_obol_polygon_handle_from_rt_ref(ref),
	       open ? TRUE : FALSE) ? 1 : 0;
}

static int
ged_draw_rt_obol_polygon_close(
    rt_view_polygon_ref ref,
    void *data)
{
    return ged_draw_rt_obol_polygon_set_open(ref, 0, data);
}

static int
ged_draw_rt_obol_polygon_clear_selected_point(
    rt_view_polygon_ref ref,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller)
	return 0;
    return token->controller->polygons().clearSelectedPoint(
	       ged_obol_polygon_handle_from_rt_ref(ref)) ? 1 : 0;
}

static int
ged_draw_rt_obol_polygon_remove(
    rt_view_polygon_ref ref,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller)
	return 0;
    return token->controller->polygons().remove(
	       ged_obol_polygon_handle_from_rt_ref(ref)) ? 1 : 0;
}

static void *
ged_draw_rt_obol_polygon_user_data(
    rt_view_polygon_ref ref,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    return (token && token->controller) ?
	   token->controller->polygons().userData(
	       ged_obol_polygon_handle_from_rt_ref(ref)) : NULL;
}

static int
ged_draw_rt_obol_polygon_user_data_set(
    rt_view_polygon_ref ref,
    void *user_data,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller)
	return 0;
    return token->controller->polygons().setUserData(
	       ged_obol_polygon_handle_from_rt_ref(ref), user_data) ? 1 : 0;
}

static int
ged_draw_rt_obol_polygon_csg(
    rt_view_polygon_ref target,
    rt_view_polygon_ref stencil,
    bg_clip_t op,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *target_token =
	ged_obol_polygon_token_from_rt_ref(target);
    ged_obol_polygon_ref_token *stencil_token =
	ged_obol_polygon_token_from_rt_ref(stencil);
    if (!target_token || !stencil_token || !target_token->controller ||
	target_token->controller != stencil_token->controller)
	return 0;
    return target_token->controller->polygons().csg(
	       ged_obol_polygon_handle_from_rt_ref(target),
	       ged_obol_polygon_handle_from_rt_ref(stencil),
	       op) ? 1 : 0;
}

static rt_view_polygon_ref
ged_draw_rt_obol_polygon_import_sketch_context(
    const char *name,
    struct db_i *dbip,
    struct directory *dp,
    void *view_ctx,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name)
	return RT_VIEW_POLYGON_REF_NULL;
    return ged_obol_rt_polygon_ref(controller,
				   controller->polygons().importSketch(name,
					   BRLObolFeatureScope::Shared, dbip, dp));
}

static struct directory *
ged_draw_rt_obol_polygon_export_sketch(
    struct db_i *dbip,
    const char *name,
    rt_view_polygon_ref ref,
    void *UNUSED(data))
{
    ged_obol_polygon_ref_token *token =
	ged_obol_polygon_token_from_rt_ref(ref);
    if (!token || !token->controller || !dbip || !name)
	return NULL;
    return token->controller->polygons().exportSketch(
	       ged_obol_polygon_handle_from_rt_ref(ref), dbip, name) ?
	   db_lookup(dbip, name, LOOKUP_QUIET) : NULL;
}

static int
ged_draw_rt_obol_polygon_snap_exclude_set(
    void *view_ctx,
    rt_view_polygon_ref ref,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !ged_obol_polygon_token_from_rt_ref(ref))
	return 0;
    return controller->polygons().setSnapExclude(
	       ged_obol_polygon_handle_from_rt_ref(ref)) ? 1 : 0;
}

extern "C" int
ged_draw_view_context_obol_polygon_adapter_attach(
    struct ged *gedp,
    void *view_ctx)
{
    if (!gedp || !view_ctx)
	return 0;

    struct rt_view_context_polygon_adapter adapter;
    memset(&adapter, 0, sizeof(adapter));
    adapter.owns_ref = ged_draw_rt_obol_polygon_owns_ref;
    adapter.record_get = ged_draw_rt_obol_polygon_record_get;
    adapter.create = ged_draw_rt_obol_polygon_create;
    adapter.select = ged_draw_rt_obol_polygon_select;
    adapter.find = ged_draw_rt_obol_polygon_find;
    adapter.dup = ged_draw_rt_obol_polygon_dup;
    adapter.visit_records = ged_draw_rt_obol_polygon_visit_records;
    adapter.snap_count = ged_draw_rt_obol_polygon_snap_count;
    adapter.clear_point_selection =
	ged_draw_rt_obol_polygon_clear_point_selection;
    adapter.update = ged_draw_rt_obol_polygon_update;
    adapter.update_screen_pt = ged_draw_rt_obol_polygon_update_screen_pt;
    adapter.move = ged_draw_rt_obol_polygon_move;
    adapter.set_name = ged_draw_rt_obol_polygon_set_name;
    adapter.set_context = ged_draw_rt_obol_polygon_set_context;
    adapter.set_visual = ged_draw_rt_obol_polygon_set_visual;
    adapter.set_open = ged_draw_rt_obol_polygon_set_open;
    adapter.close = ged_draw_rt_obol_polygon_close;
    adapter.clear_selected_point =
	ged_draw_rt_obol_polygon_clear_selected_point;
    adapter.remove = ged_draw_rt_obol_polygon_remove;
    adapter.user_data = ged_draw_rt_obol_polygon_user_data;
    adapter.user_data_set = ged_draw_rt_obol_polygon_user_data_set;
    adapter.csg = ged_draw_rt_obol_polygon_csg;
    adapter.import_sketch_context =
	ged_draw_rt_obol_polygon_import_sketch_context;
    adapter.export_sketch = ged_draw_rt_obol_polygon_export_sketch;
    adapter.snap_exclude_set = ged_draw_rt_obol_polygon_snap_exclude_set;
    adapter.data = gedp;
    return rt_view_context_polygon_adapter_set(view_ctx, &adapter);
}

struct ged_obol_selection_path_callback_context {
    rt_view_selection_path_callback_t callback;
    void *data;
};

static int
ged_obol_selection_kind_from_ged(int kind)
{
    switch (kind) {
	case -1:
	    return BRLOBOL_SELECTION_ALL;
	case 0:
	    return BRLOBOL_SELECTION_SELECTED_PATH;
	case 1:
	    return BRLOBOL_SELECTION_HIGHLIGHTED_REF;
	default:
	    return INT_MIN;
    }
}

static void
ged_obol_selection_path_callback(const SbString &path, void *data)
{
    struct ged_obol_selection_path_callback_context *ctx =
	    static_cast<struct ged_obol_selection_path_callback_context *>(data);
    const char *path_str = path.getString();
    if (ctx && ctx->callback && path_str && path_str[0])
	ctx->callback(path_str, ctx->data);
}

static int
ged_draw_rt_obol_selection_available(void *view_ctx, void *UNUSED(data))
{
    return ged_obol_view_controller_ensure_for_context(view_ctx, 1) ? 1 : 0;
}

static size_t
ged_draw_rt_obol_selection_count(void *view_ctx, void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller)
	return 0;
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    return controller->selection().count(&owner, BRLOBOL_SELECTION_ALL);
}

static int
ged_draw_rt_obol_selection_set_pick_result_context(
    void *view_ctx,
    const void *result_ctx,
    rt_view_selection_path_callback_t callback,
    void *callback_data,
    void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller)
	return 0;

    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    std::vector<BRLObolSelectionRecord> records;

    size_t result_count = result_ctx ?
			  rt_view_pick_result_context_count(result_ctx) : 0;
    for (size_t i = 0; i < result_count; i++) {
	struct bu_vls path = BU_VLS_INIT_ZERO;
	if (!rt_view_pick_result_context_path(result_ctx, i, &path)) {
	    bu_vls_free(&path);
	    continue;
	}
	const char *normalized =
	    ged_obol_skip_leading_slash(bu_vls_cstr(&path));
	if (normalized && normalized[0]) {
	    BRLObolSelectionRecord record;
	    record.path = normalized;
	    record.kind = BRLOBOL_SELECTION_SELECTED_PATH;
	    record.hitDistance =
		rt_view_pick_result_context_hit_dist(result_ctx, i);
	    record.owner = owner;
	    records.push_back(record);
	}
	bu_vls_free(&path);
    }

    struct ged_obol_selection_path_callback_context cb_ctx;
    cb_ctx.callback = callback;
    cb_ctx.data = callback_data;
    return controller->selection().applyPickResults(records,
	    callback ? ged_obol_selection_path_callback : NULL,
	    &cb_ctx, &owner) ? 1 : 0;
}

static int
ged_draw_rt_obol_selection_clear(void *view_ctx, void *UNUSED(data))
{
    BRLObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller)
	return 0;
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    controller->selection().clear(&owner, BRLOBOL_SELECTION_ALL);
    return 1;
}

extern "C" int
ged_draw_view_context_obol_selection_adapter_attach(
    struct ged *gedp,
    void *view_ctx)
{
    if (!gedp || !view_ctx)
	return 0;

    struct rt_view_context_selection_adapter adapter;
    memset(&adapter, 0, sizeof(adapter));
    adapter.available = ged_draw_rt_obol_selection_available;
    adapter.count = ged_draw_rt_obol_selection_count;
    adapter.set_pick_result_context =
	ged_draw_rt_obol_selection_set_pick_result_context;
    adapter.clear = ged_draw_rt_obol_selection_clear;
    adapter.data = gedp;
    return rt_view_context_selection_adapter_set(view_ctx, &adapter);
}

extern "C" int
ged_draw_obol_view_context_selection_available(void *view_ctx)
{
    return ged_draw_rt_obol_selection_available(view_ctx, NULL);
}

extern "C" size_t
ged_draw_obol_view_context_selection_count(void *view_ctx)
{
    return ged_draw_rt_obol_selection_count(view_ctx, NULL);
}

struct ged_obol_selection_foreach_context {
    void *view_ctx;
    ged_draw_view_context_selection_path_cb cb;
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
    void *view_ctx,
    ged_draw_view_context_selection_path_cb cb,
    void *data)
{
    if (!cb)
	return 0;
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;

    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    struct ged_obol_selection_foreach_context ctx;
    ctx.view_ctx = view_ctx;
    ctx.cb = cb;
    ctx.data = data;
    ctx.ok = 1;
    controller->selection().visitPaths(ged_obol_selection_foreach_callback,
				       &ctx, &owner, BRLOBOL_SELECTION_ALL);
    return ctx.ok;
}

extern "C" int
ged_draw_obol_view_context_selection_clear(void *view_ctx)
{
    return ged_draw_rt_obol_selection_clear(view_ctx, NULL);
}

extern "C" int
ged_draw_obol_view_context_selection_contains_path(
    void *view_ctx,
    int kind,
    const char *path)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !path || !path[0])
	return 0;

    int obol_kind = ged_obol_selection_kind_from_ged(kind);
    if (obol_kind == INT_MIN)
	return 0;

    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    return controller->selection().containsPath(
	       ged_obol_skip_leading_slash(path), obol_kind, &owner) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_selection_add_path(
    void *view_ctx,
    int kind,
    const char *path)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !path || !path[0])
	return 0;

    int obol_kind = ged_obol_selection_kind_from_ged(kind);
    if (obol_kind == INT_MIN || obol_kind == BRLOBOL_SELECTION_ALL)
	return 0;

    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    return controller->selection().addPath(ged_obol_skip_leading_slash(path),
					   obol_kind, &owner) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_selection_set_path(
    void *view_ctx,
    int kind,
    const char *path)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;

    int obol_kind = ged_obol_selection_kind_from_ged(kind);
    if (obol_kind == INT_MIN)
	return 0;

    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    if (obol_kind == BRLOBOL_SELECTION_ALL) {
	controller->selection().clear(&owner, BRLOBOL_SELECTION_ALL);
	return path && path[0] ? 0 : 1;
    }

    return controller->selection().setPath(
	       ged_obol_skip_leading_slash(path), obol_kind, &owner) ? 1 : 0;
}

extern "C" struct ged_draw_command_scene *
ged_draw_command_scene_begin(
    void *view_ctx,
    const struct ged_draw_command_scene_desc *desc)
{
    const int local = desc && desc->local;
    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller)
	return NULL;

    struct ged_draw_command_scene *scene = new ged_draw_command_scene;
    scene->magic = GED_DRAW_COMMAND_SCENE_MAGIC;
    scene->view_ctx = view_ctx;
    scene->controller = controller;
    scene->owner = ged_obol_command_scene_owner(view_ctx, desc);
    scene->controller->features().markCommandOwnerGeneration(scene->owner);
    scene->scope = local ?
		   BRLObolFeatureScope::Local : BRLObolFeatureScope::Shared;
    scene->result_cb = desc ? desc->result_cb : NULL;
    scene->result_cb_data = desc ? desc->result_cb_data : NULL;
    scene->changed = 0;
    scene->failed = 0;
    ged_obol_command_scene_notify(scene, GED_DRAW_COMMAND_RESULT_ACCEPTED,
				  NULL, "begin", NULL);
    return scene;
}

extern "C" size_t
ged_draw_command_scene_features_remove_prefix(
    struct ged_draw_command_scene *scene,
    const char *prefix)
{
    if (!prefix) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, NULL, "removePrefix",
					  "missing feature prefix");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, prefix, "removePrefix"))
	return 0;

    size_t removed = 0;
    if (scene->scope == BRLObolFeatureScope::Local) {
	removed = scene->controller->features().removePrefix(prefix,
		  BRLOBOL_FEATURE_SCOPE_LOCAL, &scene->owner);
    } else {
	removed = scene->controller->features().removePrefix(prefix,
		  BRLOBOL_FEATURE_SCOPE_SHARED, &scene->owner);
    }
    scene->changed += static_cast<int>(removed);
    if (removed) {
	char diagnostic[64];
	snprintf(diagnostic, sizeof(diagnostic), "%zu feature(s) removed",
		 removed);
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_REMOVED, prefix, "removePrefix",
				      diagnostic);
    }
    return removed;
}

extern "C" int
ged_draw_command_scene_line_layer_builder_replace(
    struct ged_draw_command_scene *scene,
    const char *name,
    const struct bg_line_layer_builder *builder,
    const struct ged_draw_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, NULL,
					  "lineLayerBuilderReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name,
					 "lineLayerBuilderReplace"))
	return 0;

    if (!builder || !bg_line_layer_builder_point_count(builder)) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"lineLayerBuilderReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureHandle handle =
	scene->controller->features().publishLineLayerBuilder(name,
	    scene->scope, builder, style ? &obol_style : NULL,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name,
				      "lineLayerBuilderReplace", "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_DRAW_COMMAND_RESULT_UPDATED, name,
				  "lineLayerBuilderReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_draw_command_scene_line_layers_replace(
    struct ged_draw_command_scene *scene,
    const char *name,
    const struct ged_draw_view_line_layer_data *layers,
    size_t layer_count,
    const struct ged_draw_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, NULL,
					  "lineLayersReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name,
					 "lineLayersReplace"))
	return 0;

    size_t live_layers = 0;
    for (size_t i = 0; layers && i < layer_count; i++) {
	if (layers[i].points && layers[i].point_count)
	    live_layers++;
    }
    if (!layers || !live_layers) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"lineLayersReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    std::vector<BRLObolLineLayer> obol_layers;
    obol_layers.reserve(live_layers);
    for (size_t i = 0; i < layer_count; i++) {
	if (!layers[i].points || !layers[i].point_count)
	    continue;

	BRLObolLineLayer layer;
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

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureHandle handle =
	scene->controller->features().publishLineLayers(name,
	    scene->scope, obol_layers, style ? &obol_style : NULL,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name,
				      "lineLayersReplace", "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_DRAW_COMMAND_RESULT_UPDATED, name, "lineLayersReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_draw_command_scene_line_set_replace(
    struct ged_draw_command_scene *scene,
    const char *name,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    const struct ged_draw_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, NULL,
					  "lineSetReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name, "lineSetReplace"))
	return 0;

    if (point_count && !points) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name, "lineSetReplace",
				      "missing line payload");
	return 0;
    }

    if (!points || !point_count) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"lineSetReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureHandle handle =
	scene->controller->features().publishLineSet(name, scene->scope,
	    ged_obol_points_from_ged(points, point_count),
	    ged_obol_commands_from_ged(cmds, point_count),
	    style ? &obol_style : NULL, &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name, "lineSetReplace",
				      "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_DRAW_COMMAND_RESULT_UPDATED, name, "lineSetReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_draw_command_scene_point_set_replace(
    struct ged_draw_command_scene *scene,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_draw_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, NULL,
					  "pointSetReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name, "pointSetReplace"))
	return 0;

    if (point_count && !points) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name, "pointSetReplace",
				      "missing point payload");
	return 0;
    }

    if (!points || !point_count) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"pointSetReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureHandle handle =
	scene->controller->features().publishPointSet(name, scene->scope,
	    ged_obol_points_from_ged(points, point_count),
	    style ? &obol_style : NULL, &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name, "pointSetReplace",
				      "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_DRAW_COMMAND_RESULT_UPDATED, name, "pointSetReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_draw_command_scene_indexed_face_set_replace(
    struct ged_draw_command_scene *scene,
    const char *name,
    const point_t *points,
    size_t point_count,
    const vect_t *normals,
    size_t normal_count,
    const int *indices,
    size_t index_count,
    const struct ged_draw_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, NULL,
					  "indexedFaceSetReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name,
					 "indexedFaceSetReplace"))
	return 0;

    if ((point_count && !points) || (index_count && !indices)) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name,
				      "indexedFaceSetReplace", "missing indexed-face payload");
	return 0;
    }

    if (!points || !point_count || !indices || !index_count) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"indexedFaceSetReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureHandle handle =
	scene->controller->features().publishIndexedFaceSet(name,
	    scene->scope,
	    ged_obol_points_from_ged(points, point_count),
	    ged_obol_vectors_from_ged(normals, normal_count),
	    ged_obol_indices_from_ged(indices, index_count),
	    style ? &obol_style : NULL, &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name,
				      "indexedFaceSetReplace", "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_DRAW_COMMAND_RESULT_UPDATED, name, "indexedFaceSetReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_draw_command_scene_hud_label_replace(
    struct ged_draw_command_scene *scene,
    const char *name,
    const struct ged_diagnostic_hud_label *label)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, NULL, "hudLabelReplace",
					  "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name, "hudLabelReplace"))
	return 0;

    if (!label || !label->text || !label->text[0]) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"hudLabelReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    BRLObolFeatureStyle obol_style;
    obol_style.hasVisible = TRUE;
    obol_style.visible = TRUE;
    obol_style.hasSelectable = TRUE;
    obol_style.selectable = TRUE;

    std::vector<BRLObolLabel> labels;
    labels.push_back(ged_obol_label_from_hud(*label));
    BRLObolFeatureHandle handle =
	scene->controller->features().publishHudLabels(name,
	    scene->scope, labels, &obol_style, &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name, "hudLabelReplace",
				      "feature publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_DRAW_COMMAND_RESULT_UPDATED, name, "hudLabelReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_draw_command_scene_custom_node_replace(
    struct ged_draw_command_scene *scene,
    const char *name,
    ged_draw_command_scene_custom_node_cb node_cb,
    void *node_cb_data,
    const struct ged_draw_view_feature_style *style)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, NULL,
					  "customNodeReplace", "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name,
					 "customNodeReplace"))
	return 0;

    if (!node_cb) {
	(void)ged_obol_command_scene_remove_feature(scene, name,
		"customNodeReplace");
	if (scene->failed)
	    return 0;
	scene->changed++;
	return 1;
    }

    struct ged_draw_command_scene_custom_node_request request =
	    GED_DRAW_COMMAND_SCENE_CUSTOM_NODE_REQUEST_INIT;
    request.feature_name = name;
    request.owner_id = scene->owner.ownerId.getString();
    request.owner_role = scene->owner.ownerRole.getString();
    request.generation = scene->owner.generation;
    request.local = scene->scope == BRLObolFeatureScope::Local ? 1 : 0;

    SoNode *node = static_cast<SoNode *>(node_cb(&request,
					 node_cb_data));
    if (!node) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name,
				      "customNodeReplace", "custom Coin node provider returned null");
	return 0;
    }

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureHandle handle =
	scene->controller->features().publishCustomNode(name,
	    scene->scope, node, style ? &obol_style : NULL,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name,
				      "customNodeReplace", "custom Coin node publish failed");
	return 0;
    }

    (void)ged_obol_feature_mark_overlay(scene->controller, handle,
					ged_obol_command_scene_overlay_info(scene, name));
    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_DRAW_COMMAND_RESULT_UPDATED, name,
				  "customNodeReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_draw_command_scene_feature_metadata_replace(
    struct ged_draw_command_scene *scene,
    const char *name,
    const struct ged_draw_command_scene_metadata *metadata,
    size_t metadata_count)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, NULL, "metadataReplace",
					  "missing feature name");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name, "metadataReplace"))
	return 0;

    BRLObolFeatureHandle handle =
	scene->controller->features().findOwned(name,
	    scene->scope == BRLObolFeatureScope::Local ?
	    BRLOBOL_FEATURE_SCOPE_LOCAL : BRLOBOL_FEATURE_SCOPE_SHARED,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name, "metadataReplace",
				      "owned feature not found");
	return 0;
    }

    std::vector<BRLObolFeatureMetadata> obol_metadata;
    obol_metadata.reserve(metadata_count);
    for (size_t i = 0; metadata && i < metadata_count; i++) {
	if (!metadata[i].key || !metadata[i].key[0])
	    continue;
	BRLObolFeatureMetadata item;
	item.key = metadata[i].key;
	item.value = metadata[i].value ? metadata[i].value : "";
	obol_metadata.push_back(item);
    }

    if (!scene->controller->features().replaceMetadata(handle,
	    obol_metadata)) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name, "metadataReplace",
				      "metadata replace failed", handle);
	return 0;
    }

    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_DRAW_COMMAND_RESULT_UPDATED, name, "metadataReplace",
				  NULL, handle);
    return 1;
}

extern "C" int
ged_draw_command_scene_feature_primitive_metadata_replace(
    struct ged_draw_command_scene *scene,
    const char *name,
    int primitive,
    const struct ged_draw_command_scene_metadata *metadata,
    size_t metadata_count)
{
    if (!name) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, NULL,
					  "primitiveMetadataReplace", "missing feature name");
	}
	return 0;
    }
    if (primitive < 0) {
	if (ged_obol_command_scene_valid(scene)) {
	    scene->failed = 1;
	    ged_obol_command_scene_notify(scene,
					  GED_DRAW_COMMAND_RESULT_FAILED, name,
					  "primitiveMetadataReplace",
					  "invalid primitive index");
	}
	return 0;
    }

    if (!ged_obol_command_scene_writable(scene, name,
					 "primitiveMetadataReplace"))
	return 0;

    BRLObolFeatureHandle handle =
	scene->controller->features().findOwned(name,
	    scene->scope == BRLObolFeatureScope::Local ?
	    BRLOBOL_FEATURE_SCOPE_LOCAL : BRLOBOL_FEATURE_SCOPE_SHARED,
	    &scene->owner);
    if (!handle.isValid()) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name,
				      "primitiveMetadataReplace", "owned feature not found");
	return 0;
    }

    std::vector<BRLObolFeatureMetadata> obol_metadata;
    obol_metadata.reserve(metadata_count);
    for (size_t i = 0; metadata && i < metadata_count; i++) {
	if (!metadata[i].key || !metadata[i].key[0])
	    continue;
	BRLObolFeatureMetadata item;
	item.key = metadata[i].key;
	item.value = metadata[i].value ? metadata[i].value : "";
	obol_metadata.push_back(item);
    }

    if (!scene->controller->features().replacePrimitiveMetadata(handle,
	    static_cast<int32_t>(primitive), obol_metadata)) {
	scene->failed = 1;
	ged_obol_command_scene_notify(scene,
				      GED_DRAW_COMMAND_RESULT_FAILED, name,
				      "primitiveMetadataReplace",
				      "primitive metadata replace failed", handle);
	return 0;
    }

    scene->changed++;
    ged_obol_command_scene_notify(scene,
				  GED_DRAW_COMMAND_RESULT_UPDATED, name,
				  "primitiveMetadataReplace", NULL, handle);
    return 1;
}

extern "C" int
ged_draw_command_scene_commit(struct ged_draw_command_scene *scene)
{
    if (!ged_obol_command_scene_valid(scene))
	return 0;

    const int ret = (scene->failed ||
		     !scene->controller->features().commandOwnerGenerationCurrent(
			 scene->owner)) ? 0 : 1;
    ged_obol_command_scene_notify(scene,
				  ret ? GED_DRAW_COMMAND_RESULT_ACCEPTED :
				  GED_DRAW_COMMAND_RESULT_FAILED,
				  NULL, "commit",
				  ret ? NULL : "command-scene commit rejected");
    scene->magic = 0;
    delete scene;
    return ret;
}

extern "C" void
ged_draw_command_scene_abort(struct ged_draw_command_scene *scene)
{
    if (!ged_obol_command_scene_valid(scene))
	return;

    scene->magic = 0;
    delete scene;
}

extern "C" int
ged_draw_obol_view_context_line_layer_builder_replace(
    void *view_ctx,
    const char *name,
    int local,
    const struct bg_line_layer_builder *builder)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;
    if (!builder || !bg_line_layer_builder_point_count(builder))
	return ged_obol_remove_feature(controller, view_ctx, name,
				       local ? 1 : 0) ? 1 : 1;

    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    BRLObolFeatureHandle handle =
	controller->features().publishLineLayerBuilder(name,
	    local ? BRLObolFeatureScope::Local :
	    BRLObolFeatureScope::Shared, builder, NULL,
	    local ? &owner : NULL);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_diagnostic_line_layer_builder_replace(
    void *view_ctx,
    const char *name,
    const struct bg_line_layer_builder *builder)
{
    int ret = ged_draw_obol_view_context_line_layer_builder_replace(
		  view_ctx, name, 0, builder);
    if (!ret || !builder || !bg_line_layer_builder_point_count(builder))
	return ret;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, 0, 0);
    BRLObolFeatureHandle handle = controller ?
				  controller->features().find(name, BRLOBOL_FEATURE_SCOPE_SHARED) :
				  BRLObolFeatureHandle();
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BRLObolOverlayClass::Diagnostic,
						 BRLObolOverlayLifecycle::PerCommand,
						 BRLObolOverlayOrder::PostTransparent,
						 name));
}

extern "C" int
ged_draw_obol_view_context_line_layers_replace(
    void *view_ctx,
    const char *name,
    int local,
    const struct ged_draw_view_line_layer_data *layers,
    size_t layer_count,
    const struct ged_draw_view_feature_style *style)
{
    BRLObolViewController *controller =
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

    std::vector<BRLObolLineLayer> obol_layers;
    obol_layers.reserve(live_layers);
    for (size_t i = 0; i < layer_count; i++) {
	if (!layers[i].points || !layers[i].point_count)
	    continue;

	BRLObolLineLayer layer;
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

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    BRLObolFeatureHandle handle =
	controller->features().publishLineLayers(name,
	    local ? BRLObolFeatureScope::Local :
	    BRLObolFeatureScope::Shared, obol_layers,
	    style ? &obol_style : NULL, local ? &owner : NULL);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_lines_create_model_annotation(
    void *view_ctx,
    const char *name,
    int local,
    const point_t point)
{
    if (!point)
	return 0;

    point_t points[1];
    VMOVE(points[0], point);
    int cmds[1] = {GED_DRAW_VIEW_LINE_MOVE};
    struct ged_draw_view_feature_style style =
	    GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    style.color_valid = 1;
    style.color[0] = 255;
    style.color[1] = 0;
    style.color[2] = 0;
    int ret = ged_draw_obol_view_context_lines_replace(view_ctx, name, local,
	      points, cmds, 1, &style);
    if (!ret)
	return 0;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 0);
    BRLObolFeatureHandle handle = local ?
				  ged_obol_feature_handle(controller, view_ctx, name) :
				  (controller ? controller->features().find(name,
				      BRLOBOL_FEATURE_SCOPE_SHARED) : BRLObolFeatureHandle());
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BRLObolOverlayClass::UserAnnotation,
						 BRLObolOverlayLifecycle::Persistent,
						 BRLObolOverlayOrder::Model,
						 name));
}

extern "C" int
ged_draw_obol_view_context_lines_append_point(
    void *view_ctx,
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
    void *view_ctx,
    const char *name,
    int local,
    const char *text,
    const point_t point,
    const point_t target,
    int has_target)
{
    if (!text || !point || (has_target && !target))
	return 0;

    struct ged_draw_view_label_data label = GED_DRAW_VIEW_LABEL_DATA_INIT;
    label.text = text;
    VMOVE(label.point, point);
    if (has_target) {
	label.line_flag = 1;
	VMOVE(label.target, target);
    }

    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;

    BRLObolFeatureStyle style;
    style.hasColor = TRUE;
    style.color = SbColor(1.0f, 1.0f, 0.0f);

    BRLObolFeatureHandle handle =
	ged_obol_publish_labels(controller, view_ctx, name, local,
				ged_obol_labels_from_ged(&label, 1), &style);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_labels_replace(
    void *view_ctx,
    const char *name,
    int local,
    const struct ged_draw_view_label_data *labels,
    size_t label_count)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;
    if (!labels || !label_count)
	return ged_obol_remove_feature(controller, view_ctx, name,
				       local ? 1 : 0) ? 1 : 1;

    BRLObolFeatureHandle handle =
	ged_obol_publish_labels(controller, view_ctx, name, local,
				ged_obol_labels_from_ged(labels, label_count));
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_tcl_labels_replace(
    void *view_ctx,
    const char *name,
    int draw,
    const struct ged_draw_view_label_data *labels,
    size_t label_count)
{
    if (!draw || !labels || !label_count) {
	BRLObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	return (controller && name) ?
	       (ged_obol_remove_feature(controller, view_ctx, name, 1) ? 1 : 1) :
	       0;
    }
    return ged_draw_obol_view_context_labels_replace(view_ctx, name, 1,
	    labels, label_count);
}

extern "C" size_t
ged_draw_obol_view_context_label_count(void *view_ctx, const char *name)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<BRLObolLabel> labels;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().labels(lookup.handle, labels))
	return 0;
    return labels.size();
}

extern "C" int
ged_draw_obol_view_context_label_copy(
    void *view_ctx,
    const char *name,
    size_t index,
    struct bu_vls *text,
    point_t point,
    unsigned char rgb[3])
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<BRLObolLabel> labels;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().labels(lookup.handle, labels) ||
	index >= labels.size())
	return 0;

    const BRLObolLabel &label = labels[index];
    if (text) {
	bu_vls_trunc(text, 0);
	bu_vls_strcat(text, label.text.getString());
    }
    if (point)
	ged_obol_point_from_sb(point, label.point);
    if (rgb) {
	BRLObolFeatureStyle style;
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
    void *view_ctx,
    const char *name,
    size_t index,
    const point_t point)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<BRLObolLabel> labels;
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
    void *view_ctx,
    const char *name,
    struct ged_draw_view_line_style *style)
{
    if (!style)
	return 0;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    BRLObolFeatureStyle obol_style;
    if (!handle.isValid() || !controller->features().style(handle,
	    obol_style))
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
    void *view_ctx,
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
    void *view_ctx,
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
    void *view_ctx,
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
    void *view_ctx,
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
    void *view_ctx,
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
ged_obol_feature_layer_lookup(void *view_ctx,
			      const char *name,
			      size_t layer_index,
			      BRLObolLineLayer &layer)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    BRLObolFeatureRecord record;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().record(lookup.handle,
		record) || record.kind != BRLObolFeatureKind::LineLayer ||
	layer_index >= record.layers.size())
	return 0;

    layer = record.layers[layer_index];
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_layer_points_copy(
    void *view_ctx,
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

    BRLObolLineLayer layer;
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
    void *view_ctx,
    const char *name,
    size_t layer_index,
    size_t point_index,
    int *out)
{
    if (out)
	*out = 0;
    if (!out)
	return 0;

    BRLObolLineLayer layer;
    if (!ged_obol_feature_layer_lookup(view_ctx, name, layer_index, layer) ||
	point_index >= layer.commands.size())
	return 0;

    *out = ged_obol_line_command_to_ged(layer.commands[point_index]);
    return 1;
}

extern "C" int
ged_draw_obol_view_context_tcl_lines_replace(
    void *view_ctx,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_draw_view_line_style *style)
{
    if (point_count % 2)
	return 0;

    if (!points || point_count < 2) {
	BRLObolViewController *controller =
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

    struct ged_draw_view_feature_style feature_style =
	    GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    const struct ged_draw_view_feature_style *feature_stylep = NULL;
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

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BRLObolOverlayClass::TclOverlay,
						 BRLObolOverlayLifecycle::PerCommand,
						 BRLObolOverlayOrder::PostTransparent,
						 name));
}

extern "C" int
ged_draw_obol_view_context_arrow_tip_get(
    void *view_ctx,
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
    void *view_ctx,
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
    void *view_ctx,
    const char *name,
    const point_t *points,
    size_t point_count,
    const struct ged_draw_view_feature_style *style)
{
    if (point_count % 2)
	return 0;
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name)
	return 0;
    if (!points || !point_count)
	return ged_obol_remove_feature(controller, view_ctx, name, 1) ? 1 : 1;

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureHandle handle =
	ged_obol_publish_arrow(controller, view_ctx, name, 1,
			       ged_obol_points_from_ged(points, point_count),
			       style ? &obol_style : NULL);
    if (!handle.isValid())
	return 0;
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BRLObolOverlayClass::TclOverlay,
						 BRLObolOverlayLifecycle::PerCommand,
						 BRLObolOverlayOrder::PostTransparent,
						 name));
}

extern "C" int
ged_draw_obol_view_context_feature_axes_centers_copy(
    void *view_ctx,
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

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    std::vector<SbVec3f> obol_centers;
    if (!handle.isValid() || !controller->features().axesCenters(handle,
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
    void *view_ctx,
    const char *name,
    const point_t *centers,
    size_t center_count,
    fastf_t half_axes_size,
    const struct ged_draw_view_feature_style *style)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !centers || !center_count)
	return 0;

    BRLObolFeatureStyle obol_style =
	ged_obol_feature_style_from_ged(style);
    BRLObolFeatureHandle handle =
	ged_obol_publish_axes(controller, view_ctx, name, 1,
			      ged_obol_points_from_ged(centers, center_count),
			      static_cast<float>(half_axes_size),
			      style ? &obol_style : NULL);
    if (!handle.isValid())
	return 0;
    return ged_obol_feature_mark_overlay(controller, handle,
					 ged_obol_model_overlay_info(view_ctx,
						 BRLObolOverlayClass::TclOverlay,
						 BRLObolOverlayLifecycle::PerCommand,
						 BRLObolOverlayOrder::PostTransparent,
						 name));
}

static BRLObolFeatureStyle
ged_obol_axes_style_from_state(const struct ged_draw_view_axes_state *state)
{
    BRLObolFeatureStyle style;
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
    void *view_ctx,
    const char *name,
    int local,
    const struct ged_draw_view_axes_state *state)
{
    if (!state)
	return 0;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return 0;

    BRLObolFeatureStyle style = ged_obol_axes_style_from_state(state);
    point_t centers[1];
    VMOVE(centers[0], state->position);
    BRLObolFeatureHandle handle =
	ged_obol_publish_axes(controller, view_ctx, name, local,
			      ged_obol_points_from_ged(centers, 1),
			      static_cast<float>(state->size),
			      &style);
    return handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_axes_state_get(
    void *view_ctx,
    const char *name,
    struct ged_draw_view_axes_state *state)
{
    if (!state)
	return 0;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    std::vector<SbVec3f> centers;
    float half_axes_size = 0.0f;
    if (!handle.isValid() || !controller->features().axesCenters(handle,
	    centers, &half_axes_size) || centers.empty())
	return 0;

    memset(state, 0, sizeof(*state));
    ged_obol_point_from_sb(state->position, centers[0]);
    state->size = half_axes_size;
    BRLObolFeatureStyle style;
    if (controller->features().style(handle, style)) {
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
    void *view_ctx,
    const char *name,
    const struct ged_draw_view_axes_state *state)
{
    if (!state)
	return 0;

    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    BRLObolFeatureHandle handle = ged_obol_feature_handle(controller,
				  view_ctx, name);
    if (!handle.isValid())
	return ged_draw_obol_view_context_axes_create(view_ctx, name, 1, state);

    BRLObolFeatureStyle style = ged_obol_axes_style_from_state(state);
    point_t centers[1];
    VMOVE(centers[0], state->position);
    BRLObolFeatureHandle updated =
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
ged_obol_overlay_group_state_set(SoBRLSceneController *scene,
				 const char *path,
				 const unsigned char *rgb,
				 fastf_t transparency,
				 int ged_draw_mode)
{
    if (!scene || !path || !path[0])
	return;

    const SbColor color = rgb ? ged_obol_color_from_rgb(rgb) :
			  SbColor(1.0f, 1.0f, 1.0f);
    const int obol_draw_mode = ged_draw_mode ?
			       ged_obol_lod_draw_mode_from_ged(ged_draw_mode) :
			       BRLOBOL_LOD_DRAW_DIAGNOSTIC;
    (void)scene->setGroupDrawIntent(path, path, obol_draw_mode,
				    BRLOBOL_LOD_DRAW_WIRE, TRUE, 0);
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
				  int ged_draw_mode,
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
    shape->drawMode = ged_draw_mode ?
		      ged_obol_lod_draw_mode_from_ged(ged_draw_mode) :
		      BRLOBOL_LOD_DRAW_DIAGNOSTIC;
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
				    int ged_draw_mode)
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
			static_cast<int32_t>(BRLObolLineCommand::Point));
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
				      transparency, ged_draw_mode, source_type, geometry_kind);
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
				   int ged_draw_mode)
{
    if (!geometry || geometry->kind != GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET)
	return NULL;

    std::vector<SbVec3f> points = ged_obol_overlay_points(geometry);
    std::vector<int32_t> triangles = ged_obol_overlay_triangles(geometry);
    if (points.empty() || triangles.empty())
	return NULL;

    SoBRLMeshShape *shape = new SoBRLMeshShape;
    ged_obol_overlay_shape_common_set(shape, name, shape_path, basecolor,
				      transparency, ged_draw_mode, "indexed-face-set", "surface");
    shape->setIndexedTriangles(points.data(), static_cast<int>(points.size()),
			       triangles.data(), static_cast<int>(triangles.size()));
    return shape;
}

extern "C" int
ged_draw_obol_overlay_erase_name_context(
    struct ged *gedp,
    void *view_ctx,
    const char *name)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !name[0])
	return 0;

    SoBRLSceneController *scene = controller->getSceneController();
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

static SoBRLSceneGroup *
ged_draw_obol_overlay_group_for_name(struct ged *gedp,
				     const char *name,
				     std::string *group_path)
{
    if (!gedp || !name || !name[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return NULL;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return NULL;

    std::string path = ged_obol_overlay_group_path_for_name(name);
    SoGroup *group = scene->findGroup(path.c_str());
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return NULL;

    if (group_path)
	*group_path = path;
    return static_cast<SoBRLSceneGroup *>(group);
}

static int
ged_draw_obol_overlay_group_style_apply(struct ged *gedp,
					const std::string &group_path,
					SoBRLSceneGroup *group,
					const struct ged_draw_view_feature_style *style,
					struct bu_vls *result,
					const char *name)
{
    if (!gedp || group_path.empty() || !group || !style)
	return 0;

    if (style->arrow >= 0 || style->arrow_tip_length >= 0.0 ||
	style->arrow_tip_width >= 0.0) {
	if (result)
	    bu_vls_printf(result,
			  "View feature %s does not support arrow settings\n", name);
	return 0;
    }

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SbBool visible = group->visible.getValue();
    int line_style = group->lineStyle.getValue();
    int line_width = group->lineWidth.getValue();
    SbBool color_override = group->colorOverride.getValue();
    SbColor color = group->color.getValue();
    if (style->visible >= 0)
	visible = style->visible ? TRUE : FALSE;
    if (style->line_style >= 0)
	line_style = style->line_style;
    if (style->line_width >= 0)
	line_width = style->line_width;
    if (style->color_valid) {
	color_override = TRUE;
	color = ged_obol_color_from_rgb(style->color);
    }

    int changed = scene->setGroupDisplayState(group_path.c_str(),
		  visible,
		  group->selected.getValue(),
		  group->highlighted.getValue(),
		  line_style,
		  line_width,
		  group->transparency.getValue(),
		  color_override,
		  color,
		  group->materialColorValid.getValue(),
		  group->materialColor.getValue(),
		  group->materialRevision.getValue());
    return changed >= 0 ? 1 : 0;
}

static int
ged_draw_obol_database_group_display_summary(
    struct ged *gedp,
    void *view_ctx,
    const char *name,
    struct ged_draw_scene_display_summary *summary)
{
    if (!gedp || !name || !name[0] ||
	!ged_draw_obol_group_display_summary_for_path(gedp, name, summary))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string target_group =
	ged_obol_group_path_from_record_path(name);
    if (target_group.empty())
	return 0;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary source_summary;
	if (!scene->getDatabaseSourceSummary(i, source_summary) ||
	    !source_summary.valid ||
	    !ged_obol_database_source_instance_in_scope(source_summary, view_ctx))
	    continue;

	const std::string owner_group_path =
	    ged_obol_database_source_owner_group_path_from_summary(
		source_summary);
	if (ged_obol_path_equal(owner_group_path.c_str(),
				target_group.c_str()) ||
	    ged_obol_path_has_prefix(owner_group_path.c_str(),
				     target_group.c_str()))
	    return 1;
    }

    return 0;
}


static int
ged_draw_obol_database_group_remove(struct ged *gedp, void *view_ctx,
				    const char *name)
{
    if (!gedp || !name || !name[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string target_group =
	ged_obol_group_path_from_record_path(name);
    if (target_group.empty())
	return 0;

    struct ged_draw_scene_display_summary summary;
    if (!ged_draw_obol_database_group_display_summary(gedp, view_ctx, name,
	    &summary))
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary source_summary;
	if (!scene->getDatabaseSourceSummary(i, source_summary) ||
	    !source_summary.valid ||
	    !ged_obol_database_source_instance_in_scope(source_summary, view_ctx))
	    continue;

	const std::string owner_group_path =
	    ged_obol_database_source_owner_group_path_from_summary(
		source_summary);
	if (ged_obol_path_equal(owner_group_path.c_str(),
				target_group.c_str()) ||
	    ged_obol_path_has_prefix(owner_group_path.c_str(),
				     target_group.c_str()))
	    ged_obol_append_database_source_instance_key(instance_keys,
		    source_summary);
    }

    int changed = ged_obol_remove_instance_keys(instance_keys, scene);
    if (scene->removeGroup(target_group.c_str()) > 0)
	changed = 1;
    if (changed)
	scene->realizePending();
    return changed ? 1 : 0;
}


static int
ged_draw_obol_database_group_style_get(
    struct ged *gedp,
    void *view_ctx,
    const char *name,
    struct ged_draw_view_feature_style *style)
{
    if (!style)
	return 0;

    struct ged_draw_scene_display_summary summary;
    if (!ged_draw_obol_database_group_display_summary(gedp, view_ctx, name,
	    &summary))
	return 0;

    struct ged_draw_appearance_settings appearance =
	    GED_DRAW_APPEARANCE_SETTINGS_INIT;
    (void)ged_draw_obol_group_appearance_for_path(gedp, name, &appearance);

    style->visible = summary.visible ? 1 : 0;
    style->color_valid = appearance.color_override ? 1 : summary.material_valid;
    if (appearance.color_override)
	VMOVE(style->color, appearance.color);
    else if (summary.material_valid)
	VMOVE(style->color, summary.material_color);
    else
	VMOVE(style->color, appearance.color);
    style->line_width = summary.line_width;
    style->line_style = summary.line_style;
    return 1;
}


static int
ged_draw_obol_database_group_style_apply(
    struct ged *gedp,
    void *view_ctx,
    const char *name,
    const struct ged_draw_view_feature_style *style,
    struct bu_vls *result)
{
    if (!gedp || !name || !name[0] || !style)
	return 0;

    struct ged_draw_scene_display_summary summary;
    if (!ged_draw_obol_database_group_display_summary(gedp, view_ctx, name,
	    &summary))
	return 0;

    if (style->arrow >= 0 || style->arrow_tip_length >= 0.0 ||
	style->arrow_tip_width >= 0.0) {
	if (result)
	    bu_vls_printf(result,
			  "View feature %s does not support arrow settings\n", name);
	return 0;
    }

    int ok = 1;
    if (style->visible >= 0)
	ok = ged_draw_obol_group_update_display_for_path(gedp, name, 1,
	     style->visible ? 1 : 0) && ok;

    if (style->color_valid || style->line_width >= 0) {
	struct ged_draw_appearance_settings appearance =
		GED_DRAW_APPEARANCE_SETTINGS_INIT;
	if (!ged_draw_obol_group_appearance_for_path(gedp, name, &appearance))
	    return 0;
	if (style->color_valid) {
	    appearance.color_override = 1;
	    VMOVE(appearance.color, style->color);
	}
	if (style->line_width >= 0)
	    appearance.s_line_width = style->line_width;
	ok = ged_draw_obol_group_update_appearance_for_path(gedp, name,
	     &appearance) && ok;
    }

    return ok ? 1 : 0;
}


extern "C" int
ged_draw_view_context_managed_feature_remove(
    struct ged *gedp,
    void *view_ctx,
    const char *name,
    ged_draw_shape_ref shape_ref,
    struct bu_vls *result)
{
    if (!gedp || !view_ctx || !name || !name[0])
	return 0;

    if (ged_draw_view_context_feature_exists(view_ctx, name)) {
	if (!ged_draw_view_context_feature_remove(view_ctx, name)) {
	    if (result)
		bu_vls_printf(result, "No view feature named %s\n", name);
	    return 0;
	}
	return 1;
    }

    if (ged_draw_shape_ref_is_null(shape_ref)) {
	if (ged_draw_obol_overlay_erase_name_context(gedp, view_ctx, name) ||
	    ged_draw_obol_database_source_remove_for_path(gedp, name) ||
	    ged_draw_obol_database_group_remove(gedp, view_ctx, name))
	    return 1;
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }

    struct ged_draw_shape_record rec;
    if (!ged_draw_shape_record_get(gedp, shape_ref, &rec)) {
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }
    if (rec.fullpath) {
	if (result)
	    bu_vls_printf(result,
			  "View feature %s is associated with a database object - use 'erase' cmd to clear\n",
			  name);
	return 0;
    }

    if (rec.display_name && rec.display_name[0] &&
	ged_draw_obol_database_source_remove_for_path(gedp, rec.display_name))
	return 1;

    if (result)
	bu_vls_printf(result, "No Obol-managed view feature named %s\n", name);
    return 0;
}

extern "C" int
ged_draw_view_context_managed_feature_visible_get(
    struct ged *gedp,
    void *view_ctx,
    const char *name,
    ged_draw_shape_ref shape_ref,
    int *visible,
    struct bu_vls *result)
{
    if (!gedp || !view_ctx || !name || !name[0] || !visible)
	return 0;

    if (ged_draw_view_context_feature_exists(view_ctx, name)) {
	struct ged_draw_view_feature_style style =
		GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	if (!ged_draw_view_context_feature_style_get(view_ctx, name, &style)) {
	    if (result)
		bu_vls_printf(result, "No view feature named %s\n", name);
	    return 0;
	}
	*visible = style.visible ? 1 : 0;
	return 1;
    }

    if (ged_draw_shape_ref_is_null(shape_ref)) {
	std::string group_path;
	SoBRLSceneGroup *group =
	    ged_draw_obol_overlay_group_for_name(gedp, name, &group_path);
	if (group) {
	    *visible = group->visible.getValue() ? 1 : 0;
	    return 1;
	}
	struct ged_draw_scene_display_summary summary;
	if (ged_draw_obol_database_group_display_summary(gedp, view_ctx, name,
		&summary)) {
	    *visible = summary.visible ? 1 : 0;
	    return 1;
	}
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }

    struct ged_draw_shape_record rec;
    if (!ged_draw_shape_record_get(gedp, shape_ref, &rec)) {
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }
    *visible = rec.visible ? 1 : 0;
    return 1;
}

extern "C" int
ged_draw_view_context_managed_feature_visible_set(
    struct ged *gedp,
    void *view_ctx,
    const char *name,
    ged_draw_shape_ref shape_ref,
    int visible,
    struct bu_vls *result)
{
    if (!gedp || !view_ctx || !name || !name[0])
	return 0;

    if (ged_draw_view_context_feature_exists(view_ctx, name)) {
	struct ged_draw_view_feature_style style =
		GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	if (!ged_draw_view_context_feature_style_get(view_ctx, name, &style)) {
	    if (result)
		bu_vls_printf(result, "No view feature named %s\n", name);
	    return 0;
	}
	style.visible = visible ? 1 : 0;
	return ged_draw_view_context_feature_style_apply(view_ctx, name,
		&style, 0);
    }

    if (ged_draw_shape_ref_is_null(shape_ref)) {
	std::string group_path;
	SoBRLSceneGroup *group =
	    ged_draw_obol_overlay_group_for_name(gedp, name, &group_path);
	if (group) {
	    struct ged_draw_view_feature_style style =
		    GED_DRAW_VIEW_FEATURE_STYLE_INIT;
	    style.visible = visible ? 1 : 0;
	    return ged_draw_obol_overlay_group_style_apply(gedp, group_path,
		    group, &style, result, name);
	}
	struct ged_draw_scene_display_summary summary;
	if (ged_draw_obol_database_group_display_summary(gedp, view_ctx, name,
		&summary))
	    return ged_draw_obol_group_update_display_for_path(gedp, name, 1,
		    visible ? 1 : 0);
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }

    return ged_draw_shape_ref_set_visible(gedp, shape_ref, visible ? 1 : 0);
}

struct ged_draw_managed_feature_color_prefix_state {
    struct ged *gedp;
    const struct db_full_path *prefix;
    const unsigned char *rgb;
};

static int
ged_draw_managed_feature_color_prefix_cb(const struct ged_draw_shape_record *rec,
	void *ud)
{
    struct ged_draw_managed_feature_color_prefix_state *ctx =
	(struct ged_draw_managed_feature_color_prefix_state *)ud;
    if (!ctx || !rec || !rec->fullpath || !ctx->prefix ||
	!db_full_path_match_top(ctx->prefix, rec->fullpath))
	return 1;
    if (ctx->rgb)
	(void)ged_draw_shape_ref_set_color(ctx->gedp, rec->ref, ctx->rgb);
    return 1;
}

extern "C" int
ged_draw_view_context_managed_feature_style_get(
    struct ged *gedp,
    void *view_ctx,
    const char *name,
    ged_draw_shape_ref shape_ref,
    struct ged_draw_view_feature_style *style,
    struct bu_vls *result)
{
    if (!gedp || !view_ctx || !name || !name[0] || !style)
	return 0;

    struct ged_draw_view_feature_style init =
	    GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    *style = init;

    if (ged_draw_view_context_feature_exists(view_ctx, name)) {
	if (!ged_draw_view_context_feature_style_get(view_ctx, name, style)) {
	    if (result)
		bu_vls_printf(result, "No view feature named %s\n", name);
	    return 0;
	}
	if (style->arrow < 0)
	    style->arrow = 0;
	return 1;
    }

    if (ged_draw_shape_ref_is_null(shape_ref)) {
	std::string group_path;
	SoBRLSceneGroup *group =
	    ged_draw_obol_overlay_group_for_name(gedp, name, &group_path);
	if (group) {
	    style->visible = group->visible.getValue() ? 1 : 0;
	    style->color_valid = group->colorOverride.getValue() ? 1 : 0;
	    ged_obol_rgb_from_color(group->color.getValue(), style->color);
	    style->line_width = group->lineWidth.getValue();
	    style->line_style = group->lineStyle.getValue();
	    return 1;
	}
	if (ged_draw_obol_database_group_style_get(gedp, view_ctx, name, style))
	    return 1;
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }

    struct ged_draw_shape_record rec;
    if (!ged_draw_shape_record_get(gedp, shape_ref, &rec)) {
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }

    style->visible = rec.visible ? 1 : 0;
    style->color_valid = 1;
    (void)ged_draw_shape_ref_get_color(gedp, shape_ref, style->color);
    style->line_width = rec.line_width;
    return 1;
}

extern "C" int
ged_draw_view_context_managed_feature_style_apply(
    struct ged *gedp,
    void *view_ctx,
    const char *name,
    ged_draw_shape_ref shape_ref,
    const struct ged_draw_view_feature_style *style,
    int recursive,
    struct bu_vls *result)
{
    if (!gedp || !view_ctx || !name || !name[0] || !style)
	return 0;

    if (ged_draw_view_context_feature_exists(view_ctx, name))
	return ged_draw_view_context_feature_style_apply(view_ctx, name, style,
		recursive);

    if (ged_draw_shape_ref_is_null(shape_ref)) {
	std::string group_path;
	SoBRLSceneGroup *group =
	    ged_draw_obol_overlay_group_for_name(gedp, name, &group_path);
	if (group)
	    return ged_draw_obol_overlay_group_style_apply(gedp, group_path,
		    group, style, result, name);
	if (ged_draw_obol_database_group_style_apply(gedp, view_ctx, name, style,
		result))
	    return 1;
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }

    struct ged_draw_shape_record rec;
    if (!ged_draw_shape_record_get(gedp, shape_ref, &rec)) {
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }

    if (style->arrow >= 0 || style->arrow_tip_length >= 0.0 ||
	style->arrow_tip_width >= 0.0) {
	if (result)
	    bu_vls_printf(result,
			  "View feature %s does not support arrow settings\n", name);
	return 0;
    }
    if (style->line_width >= 0 || style->line_style >= 0) {
	if (result)
	    bu_vls_printf(result,
			  "View feature %s does not support line style settings\n", name);
	return 0;
    }

    int ok = 1;
    if (style->visible >= 0)
	ok = ged_draw_shape_ref_set_visible(gedp, shape_ref,
					    style->visible ? 1 : 0) && ok;
    if (style->color_valid) {
	ok = ged_draw_shape_ref_set_color(gedp, shape_ref, style->color) && ok;
	if (recursive && rec.fullpath) {
	    struct ged_draw_managed_feature_color_prefix_state ctx;
	    ctx.gedp = gedp;
	    ctx.prefix = rec.fullpath;
	    ctx.rgb = style->color;
	    ged_draw_foreach_shape_record(gedp,
					  ged_draw_managed_feature_color_prefix_cb, &ctx);
	}
    }

    return ok ? 1 : 0;
}

extern "C" int
ged_draw_view_context_managed_feature_realize(
    struct ged *gedp,
    void *view_ctx,
    const char *name,
    ged_draw_shape_ref shape_ref,
    struct bu_vls *result)
{
    if (!gedp || !view_ctx || !name || !name[0])
	return 0;

    if (ged_draw_view_context_feature_exists(view_ctx, name)) {
	(void)ged_draw_view_context_feature_realize(view_ctx, name, 1);
	return 1;
    }

    if (ged_draw_shape_ref_is_null(shape_ref)) {
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }

    struct ged_draw_shape_record rec;
    memset(&rec, 0, sizeof(rec));
    if (!ged_draw_shape_record_get(gedp, shape_ref, &rec)) {
	if (result)
	    bu_vls_printf(result, "No view feature named %s\n", name);
	return 0;
    }

    std::string path;
    if (rec.fullpath) {
	char *fullpath = db_path_to_string(rec.fullpath);
	if (fullpath) {
	    path = fullpath;
	    bu_free(fullpath, "view feature realize fullpath");
	}
    }
    if (path.empty() && rec.display_name && rec.display_name[0])
	path = rec.display_name;

    if (path.empty() ||
	!ged_draw_obol_database_source_realize_for_path(gedp,
		path.c_str())) {
	if (result)
	    bu_vls_printf(result, "No Obol database source named %s\n",
			  name);
	return 0;
    }

    return 1;
}

extern "C" int
ged_draw_obol_overlay_geometry_insert_context(
    struct ged *gedp,
    void *view_ctx,
    const char *name,
    const struct ged_draw_overlay_geometry *geometry,
    const unsigned char basecolor[3],
    fastf_t transparency,
    int ged_draw_mode,
    char **shape_path_out)
{
    if (shape_path_out)
	*shape_path_out = NULL;
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !name[0] || !geometry || !basecolor)
	return 0;

    SoBRLSceneController *scene = controller->getSceneController();
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
				     transparency, ged_draw_mode);

    SoNode *shape = NULL;
    if (geometry->kind == GED_DRAW_OVERLAY_GEOMETRY_INDEXED_FACE_SET)
	shape = ged_obol_overlay_mesh_shape_create(name, shape_path, geometry,
		basecolor, transparency, ged_draw_mode);
    else
	shape = ged_obol_overlay_vlist_shape_create(name, shape_path, geometry,
		basecolor, transparency, ged_draw_mode);
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
    const BRLObolDatabaseSourceSummary &summary)
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
    if (!gedp || !path || !path[0] ||
	ged_obol_source_summary_force_adapter > 0 ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return NULL;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return NULL;

    if (ged_obol_database_source_publication_depth > 0) {
	if (ged_obol_database_source_publication_mode >= 0) {
	    const std::string mode_key =
		ged_obol_database_source_instance_key_for_mode(
		    ged_obol_database_source_publication_view, path,
		    ged_obol_database_source_publication_mode);
	    SoBRLDatabaseSource *source =
		scene->findDatabaseSourceInstance(mode_key.c_str());
	    if (source)
		return source;
	}

	const std::string base_key =
	    ged_obol_database_source_instance_key(
		ged_obol_database_source_publication_view, path);
	SoBRLDatabaseSource *source =
	    scene->findDatabaseSourceInstance(base_key.c_str());
	return source;
    }

    return scene->findDatabaseSource(path);
}

extern "C" int
ged_draw_obol_database_source_publication_begin(
    struct ged *gedp,
    void *view_ctx,
    int publication_draw_mode)
{
    if (!gedp || !ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    if (ged_obol_database_source_publication_depth == 0) {
	ged_obol_database_source_publication_scene = scene;
	scene->beginSceneMutationBatch();
    }
    ged_obol_database_source_publication_depth++;
    ged_obol_database_source_publication_view = view_ctx;
    ged_obol_database_source_publication_mode = publication_draw_mode;
    return 1;
}

extern "C" void
ged_draw_obol_database_source_publication_appearance_set(
    struct ged *gedp,
    const struct ged_draw_appearance_settings *settings)
{
    if (!gedp || !settings || ged_obol_database_source_publication_depth <= 0)
	return;

    ged_obol_database_source_publication_appearance = 1;
    ged_obol_database_source_publication_line_width = settings->s_line_width;
    ged_obol_database_source_publication_transparency =
	static_cast<float>(settings->transparency);
    ged_obol_database_source_publication_color_override =
	settings->color_override ? 1 : 0;
    ged_obol_database_source_publication_color[0] = settings->color[0];
    ged_obol_database_source_publication_color[1] = settings->color[1];
    ged_obol_database_source_publication_color[2] = settings->color[2];
}

extern "C" int
ged_draw_obol_database_source_publication_appearance_active(struct ged *gedp)
{
    return (gedp && ged_obol_database_source_publication_depth > 0 &&
	    ged_obol_database_source_publication_appearance) ? 1 : 0;
}

extern "C" void
ged_draw_obol_database_source_publication_group_set(
    struct ged *gedp,
    const char *group_path)
{
    if (!gedp || ged_obol_database_source_publication_depth <= 0)
	return;

    ged_obol_database_source_publication_group_path =
	group_path ? ged_obol_group_path_from_record_path(group_path) :
	std::string();
}

extern "C" void
ged_draw_obol_database_source_publication_end(struct ged *gedp)
{
    if (ged_obol_database_source_publication_depth > 0)
	ged_obol_database_source_publication_depth--;
    if (ged_obol_database_source_publication_depth == 0) {
	if (ged_obol_database_source_publication_scene) {
	    ged_obol_database_source_publication_scene->
	    endSceneMutationBatch();
	    ged_obol_database_source_publication_scene = NULL;
	} else if (gedp) {
	    SoBRLSceneController *scene =
		ged_draw_obol_scene_controller(gedp);
	    if (scene)
		scene->endSceneMutationBatch();
	}
	ged_obol_database_source_publication_view = NULL;
	ged_obol_database_source_publication_mode = -1;
	ged_obol_database_source_publication_group_path.clear();
	ged_obol_database_source_publication_appearance = 0;
    }
}

static int
ged_obol_database_source_instance_key_for_source(
    SoBRLDatabaseSource *source,
    std::string &source_instance_key)
{
    source_instance_key.clear();
    if (!source)
	return 0;

    BRLObolDatabaseSourceSummary summary;
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
    SoBRLSceneController **scene_out,
    std::string &source_instance_key)
{
    source_instance_key.clear();
    if (scene_out)
	*scene_out = NULL;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source && ged_obol_database_source_publication_depth > 0) {
	if (ged_obol_database_source_publication_mode >= 0) {
	    source_instance_key =
		ged_obol_database_source_instance_key_for_mode(
		    ged_obol_database_source_publication_view, path,
		    ged_obol_database_source_publication_mode);
	    source = scene->findDatabaseSourceInstance(source_instance_key.c_str());
	}
	if (!source) {
	    source_instance_key = ged_obol_database_source_instance_key(
				      ged_obol_database_source_publication_view, path);
	    source = scene->findDatabaseSourceInstance(source_instance_key.c_str());
	}
    }
    if (!source) {
	source_instance_key = path;
	source = scene->findDatabaseSourceInstance(source_instance_key.c_str());
	if (!source)
	    source = scene->findDatabaseSource(path);
    }
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
ged_obol_database_source_mark_published_current(SoBRLSceneController *scene,
	SoBRLDatabaseSource *source)
{
    if (!scene || !source)
	return 0;

    BRLObolDatabaseSourceSummary summary;
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
    SoBRLSceneController *scene,
    const char *path,
    BRLObolDatabaseSourceSummary &summary)
{
    summary = BRLObolDatabaseSourceSummary();
    if (!scene || !path || !path[0])
	return false;

    return scene->getDatabaseSourceSummaryForPath(path, summary) &&
	   summary.valid;
}

static bool
ged_obol_database_source_controller_summary_for_source(
    SoBRLSceneController *scene,
    SoBRLDatabaseSource *source,
    BRLObolDatabaseSourceSummary &summary)
{
    summary = BRLObolDatabaseSourceSummary();
    if (!source)
	return false;

    BRLObolDatabaseSourceSummary source_summary;
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
    const BRLObolDatabaseSourceSummary &summary)
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
    if (!gedp || !min || !max || !empty_out ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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

extern "C" int
ged_draw_obol_scene_database_autoview_bounds(
    struct ged *gedp,
    vect_t *min,
    vect_t *max,
    int *empty_out)
{
    if (!gedp || !min || !max || !empty_out ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    VSETALL(*min, INFINITY);
    VSETALL(*max, -INFINITY);
    *empty_out = 1;

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
    if (!gedp || !path || !path[0] || !min || !max || !empty_out ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (!out || !gedp || !ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (info->name)
	bu_free(info->name, "GED Obol scene context name");
    memset(info, 0, sizeof(*info));
}

static int
ged_obol_scene_context_info_from_summary(
    const BRLObolSceneTreeSummary &summary,
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
	 summary.nodeKind == BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE) ?
	1 : 0;
    out->is_shape =
	(!out->is_database_source &&
	 (summary.isShape ||
	  summary.nodeKind == BRLObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	  summary.nodeKind == BRLObolSceneTreeSummary::NODE_MESH_SHAPE)) ?
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    BRLObolSceneTreeSummary summary;
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    BRLObolSceneTreeSummary summary;
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
    SoBRLSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    SoBRLDatabaseSource *source =
	scene ? scene->findDatabaseSourceInstance(source_instance_key.c_str()) :
	NULL;
    if (!source)
	return 0;

    BRLObolDatabaseSourceSummary summary;
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
    out->stale_reason = ged_draw_stale_reason_name(reason);
    return 1;
}

extern "C" int
ged_draw_obol_database_source_mark_stale_for_path(
    struct ged *gedp,
    const char *path,
    ged_draw_stale_reason reason)
{
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed = scene->markDatabaseSourceInstanceStale(
			    source_instance_key.c_str(),
			    ged_obol_stale_reason_from_ged(reason));
    return changed >= 0 ? 1 : 0;
}

static int
ged_obol_database_source_record_from_summary(
    struct ged_draw_obol_database_source_record *out,
    const BRLObolDatabaseSourceSummary &summary)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    if (!summary.valid)
	return 0;

    out->valid = 1;
    out->database_path = summary.path.getString();
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

    BRLObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    return ged_obol_database_source_record_from_summary(out, summary);
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
    BRLObolDatabaseSourceSummary summary;
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
    struct BRLObolMeshLod *lod)
{
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    return source->setMeshLodBounds(
	       SbVec3f((float)bmin[X], (float)bmin[Y], (float)bmin[Z]),
	       SbVec3f((float)bmax[X], (float)bmax[Y], (float)bmax[Z]));
}

extern "C" int
ged_draw_obol_database_source_apply_record_for_path(
    struct ged *gedp,
    const char *path,
    const struct ged_draw_obol_database_source_record *record)
{
    if (!record)
	return 0;

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    BRLObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
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

    const int changed = scene->setDatabaseSourceInstanceRealizationState(
			    source_instance_key.c_str(),
			    status,
			    ged_obol_fold_revision(realized_source_revision),
			    ged_obol_fold_revision(realized_inputs_revision),
			    stale_reason,
			    failed ? "GED source realization failed" : NULL);
    return changed >= 0 ? 1 : 0;
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

    BRLObolDatabaseSourceSummary summary;
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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed =
	scene->setDatabaseSourceInstanceRealizationRoleFlags(
	    source_instance_key.c_str(), role_flags);
    return changed >= 0 ? 1 : 0;
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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const uint32_t clamped_bot_threshold =
	bot_threshold > UINT32_MAX ? UINT32_MAX : (uint32_t)bot_threshold;
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
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_realize_for_path(struct ged *gedp,
	const char *path)
{
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	source = scene->findDatabaseSource(path);
    if (!source)
	return 0;

    (void)scene->realizePending();

    BRLObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    return summary.realizationStatus == SoBRLDatabaseSource::REALIZED &&
	   !source->needsRealization() ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_realize_pending(struct ged *gedp)
{
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    scene->realizePending();
    return 1;
}

extern "C" int
ged_draw_obol_database_source_ensure_for_path(
    struct ged *gedp,
    const char *path,
    struct db_i *dbip,
    int ged_draw_mode,
    uint64_t source_revision)
{
    if (!gedp || !path || !path[0] || !dbip ||
	!ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const uint32_t folded_revision =
	source_revision ? ged_obol_fold_revision(source_revision) : 0;
    if (ged_obol_database_source_publication_depth > 0 &&
	ged_obol_database_source_publication_mode >= 0) {
	/* Batched draw publication already has a canonical full path from the
	 * database walk.  Do not scan retained source records for every leaf:
	 * large initial draws otherwise become O(N^2) as the Obol tree grows.
	 * Exact existing instances are still preserved by ged_obol_replace_path's
	 * direct instance lookup.
	 */
	const int changed = ged_obol_replace_path(gedp,
			    ged_obol_database_source_publication_view, dbip, path,
			    ged_obol_database_source_publication_mode, folded_revision,
			    scene, 0, 1);
	return changed >= 0 ? 1 : 0;
    }

    const std::string instance_key =
	ged_obol_database_source_publication_depth > 0 ?
	ged_obol_database_source_instance_key(
	    ged_obol_database_source_publication_view, path) :
	std::string(path);
    const int changed = (instance_key == path) ?
			scene->replaceDatabaseSource(path, dbip,
			    ged_obol_database_draw_mode_from_ged(ged_draw_mode),
			    folded_revision) :
			scene->replaceDatabaseSourceInstance(instance_key.c_str(), path,
			    dbip, ged_obol_database_draw_mode_from_ged(ged_draw_mode),
			    folded_revision);
    (void)ged_obol_apply_view_lod_policy(gedp, NULL, scene,
					 instance_key.c_str());

    BRLObolDatabaseSourceSummary summary;
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
						    ged_obol_lod_draw_mode_from_ged(ged_draw_mode),
						    BRLOBOL_LOD_DRAW_WIRE,
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
ged_draw_obol_database_source_ensure_for_path_with_placement(
    struct ged *gedp,
    const char *path,
    struct db_i *dbip,
    int ged_draw_mode,
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
	return ged_draw_obol_database_source_ensure_for_path(gedp, path, dbip,
		ged_draw_mode, source_revision);

    if (!gedp || !path || !path[0] || !dbip ||
	!ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const uint32_t folded_revision =
	source_revision ? ged_obol_fold_revision(source_revision) : 0;
    if (ged_obol_database_source_publication_depth > 0 &&
	ged_obol_database_source_publication_mode >= 0) {
	const int changed = ged_obol_replace_path(gedp,
			    ged_obol_database_source_publication_view, dbip, path,
			    ged_obol_database_source_publication_mode, folded_revision,
			    scene, 0, 1, NULL, &placement);
	return changed >= 0 ? 1 : 0;
    }

    if (!ged_draw_obol_database_source_ensure_for_path(gedp, path, dbip,
	    ged_draw_mode, source_revision))
	return 0;

    return ged_draw_obol_database_source_set_placement_for_path(
	       gedp, path,
	       draw_mat_valid, draw_mat,
	       draw_center_valid, draw_center,
	       draw_size_valid, draw_size);
}

extern "C" int
ged_draw_obol_database_source_rename_for_path(
    struct ged *gedp,
    const char *path,
    const char *new_path,
    uint64_t source_revision)
{
    if (!gedp || !path || !path[0] || !new_path || !new_path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !scene->findDatabaseSource(path))
	return 0;

    BRLObolDatabaseSourceSummary source_summary;
    std::string old_owner_group_path =
	ged_obol_top_group_path_from_record_path(path);
    int have_source_summary = 0;
    if (ged_obol_database_source_controller_summary_for_path(scene, path,
	    source_summary)) {
	have_source_summary = 1;
	const std::string candidate_owner =
	    ged_obol_database_source_owner_group_path_from_summary(
		source_summary);
	if (old_owner_group_path.empty() &&
	    ged_obol_path_equal(candidate_owner.c_str(), path))
	    old_owner_group_path = candidate_owner;
    }

    const SbBool source_visible =
	(!have_source_summary || source_summary.visible) ? TRUE : FALSE;
    const SbBool source_selected =
	(have_source_summary && source_summary.selected) ? TRUE : FALSE;
    const SbBool source_highlighted =
	(have_source_summary && source_summary.highlighted) ? TRUE : FALSE;
    const int source_line_style = have_source_summary ?
				  source_summary.lineStyle : 0;
    const int source_line_width = have_source_summary ?
				  source_summary.lineWidth : 0;
    const float source_transparency = have_source_summary ?
				      source_summary.transparency : 0.0f;
    const SbBool source_color_override =
	(have_source_summary && source_summary.colorOverride) ? TRUE : FALSE;
    const SbColor source_color = have_source_summary ?
				 source_summary.color : SbColor(1.0f, 1.0f, 1.0f);
    const SbBool source_material_color_valid =
	(have_source_summary && source_summary.materialColorValid) ?
	TRUE : FALSE;
    const SbColor source_material_color = have_source_summary ?
					  source_summary.materialColor : SbColor(1.0f, 1.0f, 1.0f);
    const uint32_t source_material_revision =
	have_source_summary ? source_summary.materialRevision : 0;

    int owner_draw_mode = BRLOBOL_LOD_DRAW_WIRE;
    int owner_fallback_draw_mode = BRLOBOL_LOD_DRAW_WIRE;
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
	    if (owner_fallback_draw_mode == BRLOBOL_LOD_DRAW_UNKNOWN)
		owner_fallback_draw_mode = BRLOBOL_LOD_DRAW_WIRE;
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
    const int changed = scene->renameDatabaseSource(path, new_path,
			folded_revision);
    if (changed > 0) {
	BRLObolDatabaseSourceSummary renamed_summary;
	if (ged_obol_database_source_controller_summary_for_path(scene,
		new_path, renamed_summary)) {
	    (void)scene->setDatabaseSourceState(new_path,
						FALSE,
						renamed_summary.sourceRevision,
						renamed_summary.inputsRevision,
						source_visible,
						source_selected,
						source_highlighted,
						source_line_style,
						source_line_width,
						source_transparency,
						source_color_override,
						source_color,
						source_material_color_valid,
						source_material_color,
						source_material_revision);
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
		(void)scene->moveDatabaseSourceToGroup(new_path,
						       new_owner_group_path.c_str());
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
	!group_path || !group_path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::string source_instance_key(source_path);
    if (ged_obol_database_source_publication_depth > 0 &&
	ged_obol_database_source_publication_mode >= 0) {
	const std::string mode_key =
	    ged_obol_database_source_instance_key_for_mode(
		ged_obol_database_source_publication_view, source_path,
		ged_obol_database_source_publication_mode);
	if (scene->findDatabaseSourceInstance(mode_key.c_str()))
	    source_instance_key = mode_key;
    }
    if (!scene->findDatabaseSourceInstance(source_instance_key.c_str()) &&
	!scene->findDatabaseSource(source_path))
	return 0;

    const std::string target_group =
	ged_obol_group_path_from_record_path(group_path);
    if (target_group.empty())
	return 0;

    const int moved = source_instance_key == source_path ?
		      scene->moveDatabaseSourceToGroup(source_path, target_group.c_str()) :
		      scene->moveDatabaseSourceInstanceToGroup(source_instance_key.c_str(),
			  target_group.c_str());
    return moved >= 0 ? 1 : 0;
}

static bool
ged_obol_group_path_is_overlay(SoBRLSceneController *scene,
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
    if (!out || !gedp || !ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
	BRLObolDatabaseSourceSummary summary;
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
    if (!gedp || !cb || !ged_draw_obol_scene_controller_owned(gedp))
	return -1;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    if (source_count <= 0)
	return 1;

    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
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
    if (!gedp || !cb || !ged_draw_obol_scene_controller_owned(gedp))
	return -1;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    if (source_count <= 0)
	return 1;

    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;

	if (skip_overlay_groups &&
	    ged_obol_group_path_is_overlay(scene,
					   summary.parentGroupPath.getString()))
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
    if (!gedp || !cb || !ged_draw_obol_scene_controller_owned(gedp))
	return -1;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (!gedp || !group_path || !group_path[0] || !cb ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return -1;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const std::string target_group =
	ged_obol_group_path_from_record_path(group_path);
    if (target_group.empty())
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
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
    if (!gedp || !group_path || !group_path[0] || !cb ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return -1;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const std::string target_group =
	ged_obol_group_path_from_record_path(group_path);
    if (target_group.empty())
	return -1;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
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
    if (!gedp || !path || !path[0] || !out ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    BRLObolDatabaseSourceSummary summary;
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
    if (!gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (!gedp || !path_prefix || !path_prefix[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    void *view_ctx)
{
    return ged_draw_obol_database_sources_remove_for_path_prefix_in_scope_mode(
	       gedp, path_prefix, view_ctx, -1);
}

extern "C" int
ged_draw_obol_database_sources_remove_for_path_prefix_in_scope_mode(
    struct ged *gedp,
    const char *path_prefix,
    void *view_ctx,
    int mode)
{
    if (!gedp || !path_prefix || !path_prefix[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> instance_keys;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
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
    if (!gedp || !name || !name[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> targets;
    ged_obol_append_unique_path(targets, name);

    std::vector<std::string> paths;
    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
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
    if (!gedp || !ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
	void *view_ctx)
{
    if (!gedp || !ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    int changed = ged_obol_clear_database_sources_in_scope(scene, view_ctx);
    if (changed)
	scene->realizePending();
    return changed;
}

static int
ged_obol_scene_clear_controller(SoBRLSceneController *scene)
{
    if (!scene)
	return 0;

    int changed = (scene->clearDatabaseSources() > 0) ? 1 : 0;

    std::vector<std::string> group_paths;
    const int summary_count = scene->getSceneTreeSummaryCount();
    for (int i = 0; i < summary_count; i++) {
	BRLObolSceneTreeSummary tree_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !tree_summary.valid ||
	    tree_summary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
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
    if (!gedp || !ged_draw_obol_scene_controller_owned(gedp))
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
    if (!gedp || !name || !name[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    std::vector<std::string> paths;
    const int summary_count = scene->getSceneTreeSummaryCount();
    for (int i = summary_count - 1; i >= 0; i--) {
	BRLObolSceneTreeSummary tree_summary;
	BRLObolSceneDisplaySummary display_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !scene->getSceneDisplaySummary(i, display_summary) ||
	    !tree_summary.valid ||
	    !display_summary.valid ||
	    tree_summary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
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
    if (!gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (!gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (!gedp || !path || !path[0] || !new_path || !new_path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (fallback_draw_mode == BRLOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BRLOBOL_LOD_DRAW_WIRE;
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
	!subpath || !subpath[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    int ged_draw_mode,
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

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    BRLObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SbBool nextVisible = visible_valid ?
			 (visible ? TRUE : FALSE) : summary.visible;
    SbBool nextSelected = selected_valid ?
			  (selected ? TRUE : FALSE) : summary.selected;
    SbBool nextHighlighted = highlighted_valid ?
			     (highlighted ? TRUE : FALSE) : summary.highlighted;
    int nextDrawMode = draw_mode_valid ?
		       ged_obol_database_draw_mode_from_ged(ged_draw_mode) :
		       summary.drawMode;
    int nextLineStyle = line_style_valid ? line_style : summary.lineStyle;
    int nextLineWidth = line_width_valid ? line_width : summary.lineWidth;
    float nextTransparency = transparency_valid ?
			     static_cast<float>(transparency) : summary.transparency;

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
	    nextMaterialRevision = ged_obol_fold_revision(material_revision);
	} else {
	    nextMaterialRevision++;
	    if (!nextMaterialRevision)
		nextMaterialRevision = 1;
	}
    }

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    int draw_mode_changed = scene->setDatabaseSourceInstanceDrawMode(
				source_instance_key.c_str(), nextDrawMode);
    if (draw_mode_changed < 0)
	return 0;

    if (draw_mode_valid) {
	const char *representation_key = summary.representationKey.getString();
	if (!representation_key || !representation_key[0])
	    representation_key = source_instance_key.c_str();
	int representation_changed =
	    scene->setDatabaseSourceInstanceRepresentation(
		source_instance_key.c_str(),
		representation_key,
		ged_obol_database_representation_mode_from_ged(
		    ged_draw_mode));
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
    return draw_mode_changed >= 0 ? 1 : 0;
}

template <typename ShapeT>
static int
ged_obol_shape_update_display_typed(
    SoBRLSceneController *scene,
    const char *path,
    ShapeT *shape,
    int visible_valid,
    int visible,
    int selected_valid,
    int selected,
    int highlighted_valid,
    int highlighted,
    int draw_mode_valid,
    int ged_draw_mode,
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
				 ged_obol_lod_draw_mode_from_ged(ged_draw_mode),
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
    int ged_draw_mode,
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
		draw_mode_valid, ged_draw_mode,
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
		draw_mode_valid, ged_draw_mode,
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    BRLObolDatabaseSourceDisplayPatch patch;
    patch.lineWidthValid = line_width_valid ? TRUE : FALSE;
    patch.lineWidth = line_width;
    patch.transparencyValid = transparency_valid ? TRUE : FALSE;
    patch.transparency = static_cast<float>(transparency);
    patch.colorOverrideValid = color_override_valid ? TRUE : FALSE;
    patch.colorOverride = color_override ? TRUE : FALSE;
    patch.colorValid = color_valid ? TRUE : FALSE;
    if (color_valid)
	patch.color = ged_obol_color_from_rgb(color);

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    int changed = scene->setDatabaseSourceInstanceDisplayPatch(
		      source_instance_key.c_str(), patch);
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_group_update_display_for_path(
    struct ged *gedp,
    const char *path,
    int visible_valid,
    int visible)
{
    if (!gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (fallback_draw_mode == BRLOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BRLOBOL_LOD_DRAW_WIRE;
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
    if (!out || !gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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


extern "C" int
ged_draw_obol_group_paths_foreach(
    struct ged *gedp,
    int skip_overlay_groups,
    ged_draw_obol_group_path_cb cb,
    void *userdata)
{
    if (!gedp || !cb || !ged_draw_obol_scene_controller_owned(gedp))
	return -1;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return -1;

    const int summary_count = scene->getSceneTreeSummaryCount();
    for (int i = 0; i < summary_count; i++) {
	BRLObolSceneTreeSummary tree_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !tree_summary.valid ||
	    tree_summary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
	    tree_summary.path.getLength() == 0 ||
	    BU_STR_EQUAL(tree_summary.path.getString(), "/"))
	    continue;

	if (skip_overlay_groups &&
	    ged_obol_group_path_is_overlay(scene,
					   tree_summary.path.getString()))
	    continue;

	SoGroup *group = scene->findGroup(tree_summary.path.getString());
	if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    continue;

	const char *record_path = ged_obol_group_record_path(
				      static_cast<SoBRLSceneGroup *>(group));
	if (!record_path || !record_path[0])
	    continue;

	if (!(*cb)(gedp, record_path, userdata))
	    return 0;
    }

    return 1;
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
    if (!gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (!out || !gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const std::string group_path = ged_obol_group_path_from_record_path(path);
    int count = scene->getGroupDatabaseSourceCount(group_path.c_str());
    if (count < 0)
	count = 0;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
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
    if (!out || !gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (!out || !gedp || !path || !path[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (!gedp || !path || !path[0] || !settings ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (fallback_draw_mode == BRLOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BRLOBOL_LOD_DRAW_WIRE;
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
		  static_cast<float>(settings->transparency),
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
    if (!gedp || !path || !path[0] || !settings ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    next.transparency = scene_group->transparency.getValue();
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
    if (!gedp || !ged_draw_obol_scene_controller_owned(gedp))
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    if (fallback_draw_mode == BRLOBOL_LOD_DRAW_UNKNOWN)
	fallback_draw_mode = BRLOBOL_LOD_DRAW_WIRE;
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
    const BRLObolDatabaseSourceSummary &summary,
    SoBRLDatabaseSource *source)
{
    if (summary.representationMode >= 0)
	return summary.representationMode;

    const int source_ged_mode =
	ged_obol_database_draw_mode_to_ged(summary.drawMode);
    int exact_ged_mode = source_ged_mode;
    if (source_ged_mode == GED_DRAW_MODE_SHADED && source &&
	source->getRealizedMeshCount() > 0)
	exact_ged_mode = GED_DRAW_MODE_SHADED_BOTS;

    const std::string owner_group_path =
	ged_obol_database_source_owner_group_path_from_summary(summary);
    if (owner_group_path.empty())
	return exact_ged_mode;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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

    BRLObolDatabaseSourceSummary summary;
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    BRLObolDatabaseSourceSummary source_summary;
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!ged_obol_database_source_controller_summary_for_source(scene,
	    source, source_summary))
	return 0;
    if (!source_summary.valid)
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
	    BRLObolSceneDisplaySummary display_summary;
	    if (!source->getRealizedDisplaySummary(i, display_summary) ||
		!display_summary.valid ||
		!display_summary.drawMatrixValid)
		continue;
	    if (display_summary.nodeKind !=
		BRLObolSceneTreeSummary::NODE_VLIST_SHAPE &&
		display_summary.nodeKind !=
		BRLObolSceneTreeSummary::NODE_MESH_SHAPE)
		continue;
	    out->draw_mat_valid = 1;
	    ged_obol_mat_from_sbmatrix(display_summary.drawMatrix,
				       out->draw_mat);
	    break;
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
			BRLOBOL_LOD_DRAW_HIDDEN_LINE ? TRUE : FALSE;
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
			BRLOBOL_LOD_DRAW_HIDDEN_LINE ? TRUE : FALSE;
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
    SoBRLSceneController *scene,
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
    SoBRLSceneController *scene,
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !ged_draw_obol_scene_controller_owned(gedp))
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !ged_draw_obol_scene_controller_owned(gedp))
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
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    BRLObolDatabaseSourceSummary source_summary;
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    SoBRLSceneController *scene = NULL;
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

    SoBRLSceneController *scene = NULL;
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

    SoBRLSceneController *scene = NULL;
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

    BRLObolExternalLineSet line_set;
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

    SoBRLSceneController *scene = NULL;
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

    BRLObolExternalLineSet line_set;
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

    SoBRLSceneController *scene = NULL;
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

    std::vector<BRLObolExternalAnnotationSegment> obol_segments;
    obol_segments.reserve(segment_count);
    for (size_t i = 0; i < segment_count; i++) {
	const struct ged_draw_obol_annotation_segment *seg = &segments[i];
	BRLObolExternalAnnotationSegment obol_seg;
	if (seg->kind == GED_DRAW_OBOL_ANNOTATION_SEGMENT_LINE)
	    obol_seg.kind = BRLObolExternalAnnotationSegment::SEGMENT_LINE;
	else if (seg->kind == GED_DRAW_OBOL_ANNOTATION_SEGMENT_TEXT)
	    obol_seg.kind = BRLObolExternalAnnotationSegment::SEGMENT_TEXT;
	else
	    obol_seg.kind = BRLObolExternalAnnotationSegment::SEGMENT_NONE;
	obol_seg.lineStart = seg->line_start;
	obol_seg.lineEnd = seg->line_end;
	obol_seg.textRefPoint = seg->text_ref_point;
	obol_seg.text = seg->text;
	obol_segments.push_back(obol_seg);
    }

    BRLObolExternalAnnotation annotation;
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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
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

    BRLObolAuxiliaryLineSetDisplayState obol_display;
    const BRLObolAuxiliaryLineSetDisplayState *obol_display_ptr = NULL;
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

    const int changed =
	scene->publishDatabaseSourceInstanceAuxiliaryLineSet(
	    source_instance_key.c_str(),
	    name,
	    obol_points.empty() ? NULL : obol_points.data(),
	    obol_commands.empty() ? NULL : obol_commands.data(),
	    static_cast<int>(point_count),
	    obol_display_ptr);
    return changed > 0 ? 1 : 0;
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
    SoBRLDatabaseSource *owner_source =
	ged_obol_owned_database_source_for_path(gedp, path);
    std::string owner_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(owner_source,
	    owner_instance_key))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
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

    BRLObolAuxiliaryLineSetDisplayState obol_display;
    const BRLObolAuxiliaryLineSetDisplayState *obol_display_ptr = NULL;
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

    const int changed =
	scene->publishDatabaseSourceInstanceAuxiliarySourceLineSet(
	    owner_instance_key.c_str(),
	    source_path,
	    display_name ? display_name : source_path,
	    obol_points.empty() ? NULL : obol_points.data(),
	    obol_commands.empty() ? NULL : obol_commands.data(),
	    static_cast<int>(point_count),
	    obol_display_ptr);
    if (changed <= 0)
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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int cleared = scene->clearDatabaseSourceInstanceAuxiliaryShapes(
			    source_instance_key.c_str());
    return cleared > 0 ? 1 : 0;
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

    SoBRLSceneController *scene = NULL;
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

    BRLObolExternalPointSet point_set;
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene || !ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    std::vector<int32_t> triangles;
    size_t face_count = 0;
    size_t vertex_index_count = 0;
    if (!ged_obol_indexed_faces_to_triangles(indices, index_count,
	    point_count, triangles, &face_count, &vertex_index_count))
	return 0;
    if (normal_count && normal_count != vertex_index_count &&
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

    ged_obol_local_shape_apply_common_state(shape, shape_path, display_name,
					    "local-indexed-face-set", "surface", display_state);
    shape->setIndexedTriangles(obol_points.data(),
			       static_cast<int>(obol_points.size()),
			       triangles.data(),
			       static_cast<int>(triangles.size()));
    return 1;
}

extern "C" int
ged_draw_obol_database_source_clear_mesh_for_path(
    struct ged *gedp,
    const char *path)
{
    SoBRLSceneController *scene = NULL;
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
    SoBRLSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    if (!point_count || !index_count) {
	BRLObolExternalTriangleMesh empty_mesh;
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

    BRLObolExternalTriangleMesh triangle_mesh;
    triangle_mesh.points = obol_points.empty() ? NULL : obol_points.data();
    triangle_mesh.pointCount = static_cast<int>(obol_points.size());
    triangle_mesh.indices = triangles.empty() ? NULL : triangles.data();
    triangle_mesh.indexCount = static_cast<int>(triangles.size());
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
    return ged_obol_database_source_publish_indexed_face_set_for_path(gedp,
	    path, points, point_count, normals, normal_count, indices,
	    index_count, 0);
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
    return ged_obol_database_source_publish_indexed_face_set_for_path(gedp,
	    path, points, point_count, normals, normal_count, indices,
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

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    BRLObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
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
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_update_vlist_bounds_for_path(
    struct ged *gedp,
    const char *path)
{
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    SoBRLVListShape *shape = ged_obol_owned_vlist_shape_for_path(gedp, path);
    if (!shape)
	return 0;

    if (!shape->updateDrawBoundsFromPoints())
	return 0;

    BRLObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 1;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 1;

    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 1;

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
    return 1;
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

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    BRLObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SbBool nextDrawMatrixValid = draw_mat_valid ?
				 TRUE : summary.drawMatrixValid;
    SbMatrix nextDrawMatrix = draw_mat_valid ?
			      ged_obol_sbmatrix_from_mat(draw_mat) : summary.drawMatrix;
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed = scene->setDatabaseSourceInstancePlacementState(
			    source_instance_key.c_str(),
			    nextDrawMatrixValid, nextDrawMatrix, nextDrawCenterValid,
			    nextDrawCenter, nextDrawSizeValid, nextDrawSize);
    return changed >= 0 ? 1 : 0;
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
    dst->skipped_invalid += src->skipped_invalid;
    dst->remaining += src->remaining;
    dst->comb_sources += src->comb_sources;
    dst->leaf_sources += src->leaf_sources;
}

static std::vector<std::string>
ged_obol_current_source_frontier_paths(
    struct ged *gedp,
    void *view_ctx,
    const char *root_path,
    int ged_draw_mode)
{
    std::vector<std::string> targets;
    if (!root_path || !root_path[0])
	return targets;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return targets;

    targets.push_back(std::string(root_path));
    return ged_obol_matching_database_source_paths(scene, view_ctx, targets,
	    0, 1, ged_draw_mode);
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
ged_obol_apply_path_placement(struct ged *gedp,
			      void *view_ctx,
			      const char *path,
			      int ged_draw_mode)
{
    if (!gedp || !gedp->dbip || !path || !path[0])
	return 0;

    struct db_tree_state state;
    db_init_db_tree_state(&state, gedp->dbip);
    struct db_full_path full_path;
    db_full_path_init(&full_path);

    int ret = db_follow_path_for_state(&state, &full_path, path,
				       LOOKUP_QUIET);
    int applied = 0;
    if (ret == 0) {
	int scoped = ged_draw_obol_database_source_publication_begin(gedp,
		     view_ctx, ged_draw_mode);
	applied = ged_draw_obol_database_source_set_placement_for_path(gedp,
		  path, 1, state.ts_mat, 0, NULL, 0, 0);
	if (scoped)
	    ged_draw_obol_database_source_publication_end(gedp);
    }

    db_free_full_path(&full_path);
    db_free_db_tree_state(&state);
    return applied;
}

static int
ged_obol_apply_cached_path_metadata(struct ged *gedp,
				    const char *path,
				    struct ged_draw_obol_source_expansion_status *status)
{
    if (!gedp || !gedp->dbip || !path || !path[0])
	return 0;

    struct BRLObolDrawMetadataRecord metadata;
    brlobol_draw_metadata_record_init(&metadata);
    const char *metadata_path = ged_obol_skip_leading_slash(path);
    if (!metadata_path || !metadata_path[0])
	return 0;

    if (brlobol_draw_path_metadata_cache_get(gedp->dbip, metadata_path,
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

    struct BRLObolDrawMetadataRecord metadata;
    brlobol_draw_metadata_record_init(&metadata);
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
ged_obol_publish_aabb_proxy_for_path(
    struct ged *gedp,
    void *view_ctx,
    const char *path,
    const char *cache_name,
    int ged_draw_mode,
    int refresh_missing,
    struct ged_draw_obol_source_expansion_status *status)
{
    if (!gedp || !gedp->dbip || !path || !path[0] ||
	!cache_name || !cache_name[0])
	return 0;

    struct BRLObolDrawProxyRecord record;
    brlobol_draw_proxy_record_init(&record);
    if (brlobol_draw_proxy_cache_get(gedp->dbip, cache_name,
				     BRLOBOL_DRAW_CACHE_PROXY_AABB, &record) != BRLCAD_OK) {
	if (!refresh_missing)
	    return 0;
	if (brlobol_draw_proxy_cache_refresh(gedp->dbip, cache_name,
					     BRLOBOL_DRAW_CACHE_PROXY_AABB, NULL) != BRLCAD_OK ||
	    brlobol_draw_proxy_cache_get(gedp->dbip, cache_name,
					 BRLOBOL_DRAW_CACHE_PROXY_AABB, &record) != BRLCAD_OK)
	    return 0;
    }

    if (record.pointCount != 2)
	return 0;

    point_t points[24];
    int commands[24];
    ged_obol_aabb_proxy_line_set(record.points[0], record.points[1],
				 points, commands);

    int scoped = ged_draw_obol_database_source_publication_begin(gedp,
		 view_ctx, ged_draw_mode);
    int published = ged_draw_obol_database_source_publish_line_set_for_path(
			gedp, path, (const point_t *)points, commands, 24);
    if (published && status)
	status->proxy_published++;

    if (scoped)
	ged_draw_obol_database_source_publication_end(gedp);
    return published;
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
    BRLObolLodService *service,
    struct db_i *dbip,
    const char *database_id,
    uint64_t generation,
    const char *child_path,
    const char *child_name,
    int ged_draw_mode,
    uint32_t source_revision)
{
    if (!service || !service->isRunning() || !dbip || !child_path ||
	!child_path[0] || !child_name || !child_name[0])
	return 0;

    BRLObolRtProxyProvider *provider = new (std::nothrow)
    BRLObolRtProxyProvider;
    if (!provider)
	return 0;

    provider->dbip = dbip;
    provider->proxyKind = BRLOBOL_LOD_PROXY_AABB;
    provider->useRequestBounds = FALSE;

    BRLObolLodTask task;
    task.generation = generation;
    task.request.databaseId = database_id ? database_id : "";
    task.request.sourceRevision = source_revision;
    task.request.objectPath = child_path;
    task.request.objectName = child_name;
    task.request.drawMode = ged_obol_lod_draw_mode_from_ged(ged_draw_mode);
    task.request.providerId = "brlobol_draw_aabb_cache";
    task.request.providerVersion = "brlobol-cache-v1";
    task.request.qualityTier = BRLOBOL_LOD_QUALITY_PROXY;
    task.realize = brlobol_rt_proxy_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = brlobol_rt_proxy_provider_free;
    task.publishResult = FALSE;

    uint64_t task_id = service->submitIfNotActive(task);
    if (task_id == 0) {
	brlobol_rt_proxy_provider_free(provider);
	return 0;
    }

    return 1;
}

extern "C" size_t
ged_draw_obol_database_source_prewarm_child_aabb_proxies(
    struct ged *gedp,
    void *view_ctx,
    const char *path,
    int ged_draw_mode,
    size_t max_children,
    struct ged_draw_obol_source_prewarm_status *status)
{
    ged_obol_source_prewarm_status_clear(status);
    if (!gedp || !gedp->dbip || !path || !path[0])
	return 0;

    ged_obol_attached_controller *entry =
	ged_obol_attached_controller_for_view(gedp, view_ctx);
    if (!entry || !entry->view_controller)
	return 0;

    if (!entry->lod_service || !entry->lod_service->isRunning()) {
	if (!ged_draw_obol_lod_service_start(gedp, view_ctx, 0))
	    return 0;
	entry = ged_obol_attached_controller_for_view(gedp, view_ctx);
	if (!entry || !entry->lod_service || !entry->lod_service->isRunning())
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
    const uint64_t generation = entry->lod_service->beginGeneration();
    const uint32_t source_revision =
	ged_obol_fold_revision(ged_draw_scene_revision(gedp));
    std::string parent_path = ged_obol_skip_leading_slash(path);
    size_t submitted = 0;

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

	if (child.bool_op != DB_OP_UNION) {
	    if (status)
		status->skipped_non_union++;
	    continue;
	}
	if (child.record.is_instance) {
	    if (status)
		status->skipped_duplicate_instance++;
	    continue;
	}

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

	struct BRLObolDrawCacheStatus cache_status;
	brlobol_draw_cache_status_init(&cache_status);
	if (brlobol_draw_proxy_cache_status(gedp->dbip, child_name,
					    BRLOBOL_DRAW_CACHE_PROXY_AABB, &cache_status) == BRLCAD_OK &&
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
					  entry->lod_service, gedp->dbip,
					  database_id, generation,
					  child_path.c_str(), child_name,
					  ged_draw_mode, source_revision);
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
    void *view_ctx,
    const char *root_path,
    int ged_draw_mode,
    size_t max_sources,
    size_t max_children_per_source,
    struct ged_draw_obol_source_prewarm_status *status)
{
    ged_obol_source_prewarm_status_clear(status);
    if (!gedp || !root_path || !root_path[0])
	return 0;

    std::vector<std::string> frontier =
	ged_obol_current_source_frontier_paths(gedp, view_ctx, root_path,
	    ged_draw_mode);
    if (frontier.empty())
	return 0;

    size_t submitted = 0;
    size_t source_count = 0;
    for (const std::string &path : frontier) {
	if (max_sources && source_count >= max_sources)
	    break;
	source_count++;

	struct ged_draw_obol_source_prewarm_status path_status;
	ged_obol_source_prewarm_status_clear(&path_status);
	submitted += ged_draw_obol_database_source_prewarm_child_aabb_proxies(
			 gedp, view_ctx, path.c_str(), ged_draw_mode,
			 max_children_per_source, &path_status);
	ged_obol_source_prewarm_status_accumulate(status, &path_status);
    }

    return submitted;
}

static int
ged_obol_database_source_expand_children_impl(
    struct ged *gedp,
    void *view_ctx,
    const char *path,
    int ged_draw_mode,
    size_t max_children,
    int refresh_missing_proxy,
    int require_cached_leaf_proxy,
    struct ged_draw_obol_source_expansion_status *status)
{
    ged_obol_source_expansion_status_clear(status);
    if (!gedp || !gedp->dbip || !path || !path[0] ||
	!ged_draw_obol_scene_controller_attached(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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

	    if (child.bool_op != DB_OP_UNION) {
		if (status)
		    status->skipped_non_union++;
		continue;
	    }
	    if (child.record.is_instance) {
		if (status)
		    status->skipped_duplicate_instance++;
		continue;
	    }

	    const char *child_name = ged_obol_child_object_name(&child);
	    if (!child_name || !child_name[0]) {
		if (status)
		    status->skipped_invalid++;
		continue;
	    }

	    if (!child.record.is_comb && require_cached_leaf_proxy) {
		struct BRLObolDrawCacheStatus cache_status;
		brlobol_draw_cache_status_init(&cache_status);
		if (brlobol_draw_proxy_cache_status(gedp->dbip, child_name,
						    BRLOBOL_DRAW_CACHE_PROXY_AABB, &cache_status) != BRLCAD_OK ||
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
		ged_obol_database_source_instance_key_for_mode(
		    view_ctx, child_path.c_str(), ged_draw_mode);
	    if (scene->findDatabaseSourceInstance(child_instance_key.c_str())) {
		if (status)
		    status->existing++;
		continue;
	    }

	    const int replace_changed = ged_obol_replace_path(gedp, view_ctx,
					gedp->dbip, child_path.c_str(), ged_draw_mode, source_revision,
					scene, 0);
	    if (replace_changed < 0) {
		if (status)
		    status->skipped_invalid++;
		continue;
	    }
	    if (replace_changed > 0)
		changed = 1;

	    (void)ged_obol_apply_path_placement(gedp, view_ctx,
						child_path.c_str(), ged_draw_mode);
	    (void)ged_obol_apply_index_record_metadata(gedp, child_path.c_str(),
		    &child.record, status);

	    if (!child.record.is_comb) {
		(void)ged_obol_publish_aabb_proxy_for_path(gedp, view_ctx,
			child_path.c_str(), child_name, ged_draw_mode,
			refresh_missing_proxy, status);
	    }

	    if (status) {
		status->expanded++;
		if (child.record.is_comb)
		    status->comb_sources++;
		else
		    status->leaf_sources++;
	    }
	    expanded++;
	}
    }

    if (changed) {
	ged_obol_attached_controller *entry =
	    ged_obol_attached_controller_for_view(gedp, view_ctx);
	if (entry && entry->view_controller) {
	    entry->view_controller->clearViewLodState();
	    entry->view_controller->requestRender("database-source-expand");
	}
    }

    return expanded > 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_expand_children(
    struct ged *gedp,
    void *view_ctx,
    const char *path,
    int ged_draw_mode,
    size_t max_children,
    struct ged_draw_obol_source_expansion_status *status)
{
    return ged_obol_database_source_expand_children_impl(gedp, view_ctx,
	    path, ged_draw_mode, max_children, 1, 0, status);
}

static int
ged_obol_database_source_expand_visible_children_impl(
    struct ged *gedp,
    void *view_ctx,
    const char *root_path,
    int ged_draw_mode,
    size_t max_sources,
    size_t max_children_per_source,
    int refresh_missing_proxy,
    int require_cached_leaf_proxy,
    struct ged_draw_obol_source_expansion_status *status)
{
    ged_obol_source_expansion_status_clear(status);
    if (!gedp || !root_path || !root_path[0])
	return 0;

    std::vector<std::string> frontier =
	ged_obol_current_source_frontier_paths(gedp, view_ctx, root_path,
	    ged_draw_mode);
    if (frontier.empty())
	return 0;

    int changed = 0;
    size_t source_count = 0;
    for (const std::string &path : frontier) {
	if (max_sources && source_count >= max_sources)
	    break;
	source_count++;

	struct ged_draw_obol_source_expansion_status path_status;
	ged_obol_source_expansion_status_clear(&path_status);
	const int path_changed =
	    ged_obol_database_source_expand_children_impl(gedp, view_ctx,
		path.c_str(), ged_draw_mode, max_children_per_source,
		refresh_missing_proxy, require_cached_leaf_proxy, &path_status);
	if (path_changed)
	    changed = 1;
	ged_obol_source_expansion_status_accumulate(status, &path_status);
    }

    return changed;
}

extern "C" int
ged_draw_obol_database_source_expand_visible_children(
    struct ged *gedp,
    void *view_ctx,
    const char *root_path,
    int ged_draw_mode,
    size_t max_sources,
    size_t max_children_per_source,
    struct ged_draw_obol_source_expansion_status *status)
{
    return ged_obol_database_source_expand_visible_children_impl(gedp,
	    view_ctx, root_path, ged_draw_mode, max_sources,
	    max_children_per_source, 1, 0, status);
}

static size_t
ged_obol_remaining_budget(size_t budget, size_t used)
{
    if (!budget)
	return 0;
    return used >= budget ? 0 : budget - used;
}

static int
ged_obol_progressive_advance_provider(
    BRLObolViewController *controller,
    void *user_data,
    const BRLObolProgressiveOptions *options,
    BRLObolProgressiveStatus *status)
{
    BRLObolProgressiveStatus local_status;
    if (status)
	local_status.providerCount = status->providerCount;

    ged_obol_progressive_provider_data *data =
	static_cast<ged_obol_progressive_provider_data *>(user_data);
    if (!controller || !data || !data->gedp)
	return -1;

    struct ged *gedp = data->gedp;
    void *view_ctx = data->view_ctx;
    if (!options)
	options = &controller->getDefaultProgressiveOptions();

    const uint32_t flags = options->flags;
    const int refresh_missing_proxy =
	(flags & BRLOBOL_PROGRESSIVE_REFRESH_MISSING_PROXIES) ? 1 : 0;
    const int require_cached_leaf_proxy =
	(!refresh_missing_proxy &&
	 (flags & BRLOBOL_PROGRESSIVE_REQUIRE_CACHED_PROXIES)) ? 1 : 0;

    std::vector<ged_obol_drawn_source_path_mode> roots =
	ged_obol_drawn_source_path_modes(gedp, view_ctx, -1, NULL);
    if (roots.empty())
	return 0;

    size_t used_sources = 0;
    size_t used_submissions = 0;
    int changed = 0;
    for (const ged_obol_drawn_source_path_mode &root : roots) {
	if (options->maxSources && used_sources >= options->maxSources) {
	    local_status.hasMore = 1;
	    break;
	}
	if (options->maxSubmissions &&
	    used_submissions >= options->maxSubmissions) {
	    local_status.hasMore = 1;
	    break;
	}

	size_t max_sources = ged_obol_remaining_budget(options->maxSources,
			     used_sources);
	size_t max_children = options->maxChildrenPerSource;
	if (options->maxSubmissions) {
	    size_t remaining_submissions =
		ged_obol_remaining_budget(options->maxSubmissions,
					  used_submissions);
	    if (max_children == 0 || max_children > remaining_submissions)
		max_children = remaining_submissions;
	}

	struct ged_draw_obol_source_prewarm_status prewarm_status;
	ged_obol_source_prewarm_status_clear(&prewarm_status);
	size_t submitted = 0;
	if (flags & BRLOBOL_PROGRESSIVE_VISIBLE_FRONTIER) {
	    submitted =
		ged_draw_obol_database_source_prewarm_visible_child_aabb_proxies(
		    gedp, view_ctx, root.path.c_str(), root.mode, max_sources,
		    max_children, &prewarm_status);
	} else {
	    submitted =
		ged_draw_obol_database_source_prewarm_child_aabb_proxies(
		    gedp, view_ctx, root.path.c_str(), root.mode, max_children,
		    &prewarm_status);
	}
	used_submissions += submitted;
	local_status.submitted += prewarm_status.submitted;
	local_status.alreadyCached += prewarm_status.already_cached;
	local_status.remaining += prewarm_status.remaining;
	if (submitted > 0)
	    local_status.hasMore = 1;

	struct ged_draw_obol_source_expansion_status expansion_status;
	ged_obol_source_expansion_status_clear(&expansion_status);
	int root_changed = 0;
	if (flags & BRLOBOL_PROGRESSIVE_VISIBLE_FRONTIER) {
	    root_changed =
		ged_obol_database_source_expand_visible_children_impl(gedp,
		    view_ctx, root.path.c_str(), root.mode, max_sources,
		    max_children, refresh_missing_proxy,
		    require_cached_leaf_proxy, &expansion_status);
	} else {
	    root_changed =
		ged_obol_database_source_expand_children_impl(gedp, view_ctx,
		    root.path.c_str(), root.mode, max_children,
		    refresh_missing_proxy, require_cached_leaf_proxy,
		    &expansion_status);
	}
	if (root_changed) {
	    changed = 1;
	    local_status.changed = 1;
	}

	local_status.expanded += expansion_status.expanded;
	local_status.existing += expansion_status.existing;
	local_status.remaining += expansion_status.remaining;
	local_status.proxyPublished += expansion_status.proxy_published;
	local_status.metadataApplied += expansion_status.metadata_applied;
	used_sources += max_sources ? max_sources : 1;

	if (expansion_status.expanded > 0)
	    local_status.hasMore = 1;
	if (expansion_status.remaining > 0)
	    local_status.hasMore = 1;
    }

    struct ged_draw_obol_lod_service_status service_status;
    memset(&service_status, 0, sizeof(service_status));
    if (ged_draw_obol_lod_service_status(gedp, view_ctx, &service_status)) {
	local_status.pendingTasks = service_status.pending_tasks;
	local_status.inFlight = service_status.in_flight;
	local_status.queuedResults = service_status.queued_results;
	local_status.queuedCacheWrites = service_status.queued_cache_writes;
	if (service_status.pending_tasks || service_status.in_flight ||
	    service_status.queued_results ||
	    service_status.queued_cache_writes)
	    local_status.hasMore = 1;
    }

    if (local_status.hasMore && view_ctx)
	ged_view_context_refresh_request(view_ctx, GED_VIEW_REFRESH_DRAW);
    if (local_status.changed || local_status.hasMore)
	controller->requestRender(local_status.changed ?
				  "ged-progressive-update" : "ged-progressive-pending");

    if (status)
	*status = local_status;
    return (changed || local_status.hasMore) ? 1 : 0;
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

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed = scene->setDatabaseSourceInstanceBoundsState(
			    source_instance_key.c_str(), TRUE,
			    SbVec3f(static_cast<float>(bmin[0]),
				    static_cast<float>(bmin[1]),
				    static_cast<float>(bmin[2])),
			    SbVec3f(static_cast<float>(bmax[0]),
				    static_cast<float>(bmax[1]),
				    static_cast<float>(bmax[2])));
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_set_display_name_for_path(
    struct ged *gedp,
    const char *path,
    const char *display_name)
{
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    std::string source_instance_key;
    if (!ged_obol_database_source_instance_key_for_source(source,
	    source_instance_key))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed = scene->setDatabaseSourceInstanceDisplayName(
			    source_instance_key.c_str(),
			    display_name ? display_name : "");
    return changed >= 0 ? 1 : 0;
}

static const char *
ged_obol_realized_geometry_name(const BRLObolRealizedShapeSummary &summary)
{
    if (summary.shapeKind == BRLObolRealizedShapeSummary::SHAPE_MESH)
	return "indexed-face-set";
    if (summary.shapeKind == BRLObolRealizedShapeSummary::SHAPE_VLIST) {
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

extern "C" int
ged_draw_obol_database_source_geometry_summary_for_path(
    struct ged *gedp,
    const char *path,
    struct ged_draw_shape_geometry_summary *out)
{
    if (!out)
	return 0;

    memset(out, 0, sizeof(*out));
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    SoBRLVListShape *annotation_shape =
	ged_obol_owned_annotation_vlist_shape_for_source(source, path);
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

    const int count = source->getRealizedShapeSummaryCount();
    struct ged_draw_shape_geometry_summary empty_summary;
    memset(&empty_summary, 0, sizeof(empty_summary));
    int have_empty_summary = 0;
    for (int i = 0; i < count; i++) {
	BRLObolRealizedShapeSummary summary;
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

    BRLObolDatabaseSourceSummary summary;
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
    SoBRLSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    const int changed =
	scene->refreshDatabaseSourceInstanceMaterialColorFromDatabase(
	    source_instance_key.c_str(),
	    ged_obol_fold_revision(material_revision),
	    dbip);
    return changed > 0 ? 1 : 0;
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
	BRLObolRealizedShapeSummary summary;
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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    const int next_region = evaluated_region ? 1 : 0;
    int updated = 0;

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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    int updated = 0;

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

    return updated;
}

template <typename ShapeT>
static void
ged_obol_apply_draw_metadata_to_shape(
    ShapeT *shape,
    const SoBRLDatabaseSource *source,
    const struct BRLObolDrawMetadataRecord *record)
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
    const struct BRLObolDrawMetadataRecord *record)
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
    const struct BRLObolDrawMetadataRecord *record)
{
    if (!record || !record->directoryFound)
	return 0;

    SoBRLSceneController *scene = NULL;
    std::string source_instance_key;
    if (!ged_obol_database_source_scene_instance_for_path(gedp, path, &scene,
	    source_instance_key))
	return 0;

    SoBRLDatabaseSource *source = scene ?
				  scene->findDatabaseSourceInstance(source_instance_key.c_str()) : NULL;
    if (!source)
	return 0;

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
ged_obol_scene_sync_full_scene_impl(struct ged *gedp,
				    void *view_ctx,
				    uint32_t source_revision,
				    SoBRLSceneController *controller,
				    int preserve_existing_revision)
{
    SoBRLSceneController *scene = controller ?
				  controller : ged_draw_obol_scene_controller(gedp);
    if (!gedp || !scene)
	return 0;

    if (!source_revision)
	source_revision = ged_obol_fold_revision(ged_draw_scene_revision(gedp));

    struct db_i *dbip = gedp->dbip;
    std::vector<ged_obol_drawn_source_path_mode> path_modes =
	ged_obol_drawn_source_path_modes(gedp, view_ctx, -1, NULL);

    int changed = ged_obol_clear_database_sources_in_scope(scene, view_ctx);
    if (!dbip || path_modes.empty()) {
	if (changed)
	    scene->realizePending();
	return changed ? 1 : 0;
    }

    if (ged_obol_replace_path_modes(dbip, path_modes, source_revision,
				    gedp, view_ctx, scene, 1,
				    preserve_existing_revision))
	changed = 1;

    if (changed)
	scene->realizePending();
    return changed ? 1 : 0;
}

int
ged_draw_obol_scene_sync_full_scene(struct ged *gedp,
				    void *view_ctx,
				    uint32_t source_revision,
				    SoBRLSceneController *controller)
{
    return ged_obol_scene_sync_full_scene_impl(gedp, view_ctx,
	    source_revision, controller, 0);
}

static int
ged_obol_refresh_scene_material_colors(struct ged *gedp,
				       SoBRLSceneController *scene)
{
    if (!gedp || !gedp->dbip || !scene)
	return 0;

    const int source_count = scene->getDatabaseSourceCount();
    int changed = 0;
    for (int i = 0; i < source_count; i++) {
	SoBRLDatabaseSource *source = scene->getDatabaseSource(i);
	if (!source)
	    continue;

	BRLObolDatabaseSourceSummary summary;
	if (!source->getSummary(summary) || !summary.valid)
	    continue;

	const char *path = summary.path.getString();
	const char *metadata_path = ged_obol_skip_leading_slash(path);
	if (!metadata_path || !metadata_path[0])
	    continue;

	struct BRLObolDrawMetadataRecord metadata;
	brlobol_draw_metadata_record_init(&metadata);
	if (brlobol_draw_path_metadata_cache_refresh(gedp->dbip,
		metadata_path, NULL) != BRLCAD_OK)
	    continue;
	if (brlobol_draw_path_metadata_cache_get(gedp->dbip,
		metadata_path, &metadata) != BRLCAD_OK ||
	    !metadata.directoryFound)
	    continue;

	if (ged_obol_apply_draw_metadata_to_source(source, &metadata))
	    changed = 1;
    }

    return changed ? 1 : 0;
}

static int
ged_obol_apply_source_update_transaction(
    struct ged *gedp,
    void *view_ctx,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    uint32_t source_revision,
    SoBRLSceneController *scene)
{
    if (!gedp || !txn || !scene)
	return 0;

    std::vector<std::string> targets =
	ged_obol_transaction_paths(txn, result);
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
    SoBRLSceneController *scene)
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
    void *view_ctx,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    uint32_t source_revision,
    SoBRLSceneController *scene)
{
    if (!txn || !scene)
	return 0;

    std::vector<std::string> targets =
	ged_obol_transaction_paths(txn, result);
    if (targets.empty())
	return (result && result->status >= 0) ? 1 : 0;

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

    const ged_draw_stale_reason stale_reason = txn->stale_reason ?
	txn->stale_reason : GED_DRAW_STALE_SOURCE_CHANGED;
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
    SoBRLSceneController *scene)
{
    if (!scene || targets.empty())
	return 0;

    std::vector<std::string> group_paths;
    const int summary_count = scene->getSceneTreeSummaryCount();
    for (int i = summary_count - 1; i >= 0; i--) {
	BRLObolSceneTreeSummary tree_summary;
	if (!scene->getSceneTreeSummary(i, tree_summary) ||
	    !tree_summary.valid ||
	    tree_summary.nodeKind != BRLObolSceneTreeSummary::NODE_GROUP ||
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
    if (!gedp || !path_prefix || !path_prefix[0] ||
	!ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
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
    SoBRLSceneController *scene)
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
	BRLObolDatabaseSourceSummary summary;
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
    void *view_ctx,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    uint32_t source_revision,
    SoBRLSceneController *scene)
{
    if (!gedp || !txn || !scene || !gedp->dbip)
	return 0;

    std::vector<std::string> targets =
	ged_obol_transaction_paths(txn, result);
    if (targets.empty())
	targets = ged_obol_drawn_source_paths(gedp, view_ctx, -1, NULL);

    if (targets.empty())
	return (result && result->status >= 0) ? 1 : 0;

    if (!source_revision)
	source_revision = ged_obol_fold_revision(ged_draw_scene_revision(gedp));

    int handled = 0;
    ged_obol_scene_mutation_batch_scope batch(scene, targets.size(),
	    targets.size());
    for (const std::string &target : targets) {
	std::vector<std::string> matching_sources;
	ged_obol_collect_database_sources_matching(matching_sources, scene,
		view_ctx, target.c_str(), 0, 1, txn->mode);
	if (!matching_sources.empty()) {
	    if (ged_obol_replace_matching_database_sources(gedp, view_ctx,
		    matching_sources, 0, 1, source_revision, scene,
		    txn->mode))
		handled = 1;
	    continue;
	}
	std::vector<std::string> one_target;
	ged_obol_append_unique_path(one_target, target.c_str());
	std::vector<ged_obol_drawn_source_path_mode> target_source_modes =
	    ged_obol_drawn_source_path_modes(gedp, view_ctx, txn->mode,
					     &one_target);
	if (!target_source_modes.empty()) {
	    if (ged_obol_replace_path_modes(gedp->dbip, target_source_modes,
					    source_revision, gedp, view_ctx, scene, 0))
		handled = 1;
	    continue;
	}
	if (ged_obol_replace_path_and_realize(gedp, view_ctx, gedp->dbip,
					      target.c_str(),
					      ged_obol_drawn_path_mode(gedp, view_ctx, target.c_str()),
					      source_revision, scene, 0) > 0)
	    handled = 1;
    }

    if (handled)
	scene->realizePending();
    if (!handled && result && result->status >= 0)
	return 1;
    return handled;
}

static int
ged_obol_apply_cached_path_metadata_to_scene(
    struct ged *gedp,
    SoBRLSceneController *scene,
    void *view_ctx,
    const char *path,
    int ged_draw_mode)
{
    if (!gedp || !gedp->dbip || !scene || !path || !path[0])
	return 0;

    struct BRLObolDrawMetadataRecord metadata;
    brlobol_draw_metadata_record_init(&metadata);
    const char *metadata_path = ged_obol_skip_leading_slash(path);
    if (!metadata_path || !metadata_path[0])
	return 0;

    if (brlobol_draw_path_metadata_cache_get(gedp->dbip, metadata_path,
	    &metadata) != BRLCAD_OK &&
	brlobol_draw_path_metadata_cache_refresh(gedp->dbip,
		metadata_path, NULL) == BRLCAD_OK) {
	(void)brlobol_draw_path_metadata_cache_get(gedp->dbip,
		metadata_path, &metadata);
    }
    if (!metadata.directoryFound)
	return 0;

    const std::string mode_key =
	ged_obol_database_source_instance_key_for_mode(view_ctx, path,
	    ged_draw_mode);
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
    void *view_ctx,
    struct db_i *dbip,
    const std::vector<std::string> &paths,
    int ged_draw_mode,
    uint32_t source_revision,
    SoBRLSceneController *scene,
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
				  ged_draw_mode, source_revision, scene, 0,
				  preserve_existing_revision,
				  appearance_settings) > 0)
	    changed = 1;
	if (ged_obol_apply_cached_path_metadata_to_scene(gedp, scene,
		view_ctx, path.c_str(), ged_draw_mode))
	    changed = 1;
    }

    return changed;
}

int
ged_draw_obol_scene_sync_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    SoBRLSceneController *controller)
{
    if (!gedp || !txn || (result && result->status < 0))
	return 0;

    SoBRLSceneController *scene = controller ?
				  controller : ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    void *view_ctx = txn->view;
    uint32_t source_revision = ged_obol_transaction_source_revision(result);
    int changed = 0;

    switch (txn->kind) {
	case GED_DRAW_TXN_DRAW: {
	    std::vector<std::string> paths =
		ged_obol_transaction_paths(txn, result);
	    if (result && result->affected_shapes <= 0 &&
		result->affected_groups <= 0 && paths.empty())
		break;
	    const int ged_draw_mode =
		ged_obol_transaction_ged_draw_mode(gedp, txn);
	    const struct ged_draw_appearance_settings *appearance =
		    static_cast<const struct ged_draw_appearance_settings *>(
			txn->appearance);
	    if (ged_obol_transaction_defer_leaf_expansion(txn)) {
		if (!paths.empty())
		    changed = ged_obol_replace_deferred_paths(gedp, view_ctx,
			      gedp->dbip, paths, ged_draw_mode, source_revision,
			      scene, 1, appearance);
		break;
	    }
	    if (paths.empty()) {
		changed = ged_obol_scene_sync_full_scene_impl(gedp,
			  view_ctx, source_revision, scene, 1);
	    } else {
		std::vector<std::string> source_paths =
		    ged_obol_primary_matching_database_source_paths(gedp,
			view_ctx, paths, ged_draw_mode);
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
		    if (ged_obol_replace_paths(gedp->dbip, source_paths,
					       ged_draw_mode, source_revision,
					       gedp, view_ctx, scene, 1, 1,
					       appearance))
			changed = 1;
		} else {
		    changed = ged_obol_replace_paths(gedp->dbip, paths,
						     ged_draw_mode, source_revision,
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
		const int source_count = scene->getDatabaseSourceCount();
		for (int i = 0; i < source_count; i++) {
		    BRLObolDatabaseSourceSummary summary;
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
	case GED_DRAW_TXN_MATERIAL_CHANGED:
	case GED_DRAW_TXN_REFRESH_MATERIAL_COLORS:
	    changed = ged_obol_refresh_scene_material_colors(gedp, scene);
	    if (!changed && result && result->status >= 0)
		changed = 1;
	    break;
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
		changed = ged_draw_obol_database_source_rename_for_path(gedp,
			  txn->path, txn->new_path, source_revision);
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
    int full_sync)
{
    if (full_sync)
	return 1;
    if (!txn)
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

    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    std::vector<ged_obol_attached_controller> *entries =
	ged_obol_attached_controllers(gdp, 0);
    const int completed_direct_draw =
	ged_obol_transaction_is_completed_database_source_draw(gedp, txn,
	    result);
    if (!entries || entries->empty()) {
	if (!ged_draw_obol_scene_controller(gedp))
	    return 0;
	if (completed_direct_draw)
	    return 0;
	return ged_draw_obol_scene_sync_transaction(gedp, txn, result, NULL);
    }

    SoBRLSceneController *primary_scene =
	completed_direct_draw ? ged_draw_obol_scene_controller(gedp) : NULL;
    int changed = 0;
    for (ged_obol_attached_controller &entry : *entries) {
	if (!entry.scene_controller)
	    continue;
	if (completed_direct_draw && entry.scene_controller == primary_scene)
	    continue;

	if (entry.render_endpoint_only && entry.view_controller) {
	    const int shared_root_visible =
		ged_obol_view_scope_is_independent(entry.view_ctx) ? 0 : 1;
	    if (entry.render_shared_root_visible != shared_root_visible) {
		SoBRLSceneController *shared_scene =
		    ged_draw_obol_scene_controller(gedp);
		if (shared_scene &&
		    ged_obol_bind_view_render_root(entry.view_ctx,
						   shared_scene, entry.view_controller))
		    entry.render_shared_root_visible = shared_root_visible;
	    }
	}

	if (entry.render_endpoint_only) {
	    if (!ged_obol_view_scope_is_independent(entry.view_ctx))
		continue;
	    if (txn->kind == GED_DRAW_TXN_CLEAR ||
		txn->kind == GED_DRAW_TXN_TEARDOWN)
		continue;
	}

	const struct ged_draw_transaction *sync_txn = txn;
	struct ged_draw_transaction local_txn;
	int full_sync = 0;

	if (entry.use_attached_view_scope && entry.view_ctx) {
	    const int txn_independent =
		ged_obol_view_scope_is_independent(txn->view);
	    const int entry_independent =
		ged_obol_view_scope_is_independent(entry.view_ctx);

	    if (txn_independent && txn->view != entry.view_ctx)
		continue;

	    if (entry_independent && txn->view != entry.view_ctx) {
		switch (txn->kind) {
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
			/* A redraw for one concrete shared view is a render
			 * refresh, not an instruction to rebuild independent
			 * attached scenes.  Database edits still use the
			 * source/material/highlight cases above to refresh all
			 * scopes that contain the affected sources. */
			if (txn->view)
			    continue;
			full_sync = 1;
			break;
		    case GED_DRAW_TXN_CLEAR:
		    case GED_DRAW_TXN_TEARDOWN:
			break;
		    default:
			continue;
		}
	    }

	    local_txn = *txn;
	    local_txn.view = entry.view_ctx;
	    sync_txn = &local_txn;
	}

	int entry_changed = full_sync ?
			    ged_draw_obol_scene_sync_full_scene(gedp, entry.view_ctx,
				ged_obol_transaction_source_revision(result),
				entry.scene_controller) :
			    ged_draw_obol_scene_sync_transaction(gedp, sync_txn, result,
				entry.scene_controller);
	if (entry.view_controller &&
	    ged_obol_transaction_invalidates_view_lod(sync_txn,
		    full_sync))
	    entry.view_controller->clearViewLodState();
	if (entry_changed)
	    changed = 1;
    }

    return changed ? 1 : 0;
}

int
ged_draw_obol_sync_full_scene(struct ged *gedp,
			      void *view_ctx,
			      uint32_t source_revision,
			      BRLObolViewController *controller)
{
    SoBRLSceneController *scene = controller ?
				  controller->getSceneController() : ged_draw_obol_scene_controller(gedp);
    return ged_draw_obol_scene_sync_full_scene(gedp, view_ctx,
	    source_revision, scene);
}

int
ged_draw_obol_sync_transaction(
    struct ged *gedp,
    const struct ged_draw_transaction *txn,
    const struct ged_draw_transaction_result *result,
    BRLObolViewController *controller)
{
    SoBRLSceneController *scene = controller ?
				  controller->getSceneController() : ged_draw_obol_scene_controller(gedp);
    return ged_draw_obol_scene_sync_transaction(gedp, txn, result, scene);
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
}

static void
ged_obol_collect_preserved_sources(
    SoBRLSceneController *scene,
    std::vector<ged_obol_preserved_source> &sources)
{
    if (!scene)
	return;

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) ||
	    summary.path.getLength() == 0)
	    continue;

	ged_obol_preserved_source source;
	source.path = summary.path.getString();
	source.drawMode = ged_obol_lod_draw_mode_to_ged(summary.drawMode);
	sources.push_back(source);
    }

    if (sources.size() < 2)
	return;

    std::vector<std::string> source_paths;
    for (const ged_obol_preserved_source &source : sources)
	ged_obol_append_unique_path(source_paths, source.path.c_str());
    ged_obol_remove_shadowed_source_paths(source_paths);

    std::vector<ged_obol_preserved_source> filtered;
    for (const ged_obol_preserved_source &source : sources) {
	for (const std::string &path : source_paths) {
	    if (ged_obol_path_equal(source.path.c_str(), path.c_str())) {
		filtered.push_back(source);
		break;
	    }
	}
    }
    sources.swap(filtered);
}

static void
ged_obol_preserved_sources_clear(struct ged_drawable *gdp)
{
    if (!gdp || !gdp->gd_obol_preserved_sources)
	return;
    std::vector<ged_obol_preserved_source> *sources =
	static_cast<std::vector<ged_obol_preserved_source> *>(
	    gdp->gd_obol_preserved_sources);
    delete sources;
    gdp->gd_obol_preserved_sources = NULL;
}

extern "C" void
ged_draw_obol_preserved_sources_free(struct ged *gedp)
{
    ged_obol_preserved_sources_clear(ged_obol_gdp(gedp));
}

static void
ged_obol_preserved_sources_store(
    struct ged_drawable *gdp,
    const std::vector<ged_obol_preserved_source> &sources)
{
    if (!gdp)
	return;

    ged_obol_preserved_sources_clear(gdp);
    if (sources.empty())
	return;

    gdp->gd_obol_preserved_sources =
	new std::vector<ged_obol_preserved_source>(sources);
}

static void
ged_obol_preserved_sources_take(
    struct ged_drawable *gdp,
    std::vector<ged_obol_preserved_source> &sources)
{
    sources.clear();
    if (!gdp || !gdp->gd_obol_preserved_sources)
	return;

    std::vector<ged_obol_preserved_source> *stored =
	static_cast<std::vector<ged_obol_preserved_source> *>(
	    gdp->gd_obol_preserved_sources);
    sources = *stored;
    delete stored;
    gdp->gd_obol_preserved_sources = NULL;
}

static int
ged_obol_replay_preserved_sources(
    struct ged *gedp,
    const std::vector<ged_obol_preserved_source> &sources,
    SoBRLSceneController *scene)
{
    if (!gedp || !gedp->dbip || !scene || sources.empty())
	return 0;

    uint32_t source_revision =
	ged_obol_fold_revision(ged_draw_scene_revision(gedp));
    int changed = 0;
    ged_obol_scene_mutation_batch_scope batch(scene, sources.size(),
	    sources.size());
    for (const ged_obol_preserved_source &source : sources) {
	if (ged_obol_replace_path_and_realize(gedp, NULL, gedp->dbip,
					      source.path.c_str(), source.drawMode, source_revision,
					      scene) > 0)
	    changed = 1;
    }

    if (changed)
	scene->realizePending();
    return changed;
}

static int
ged_obol_replay_preserved_sources_not_current(
    struct ged *gedp,
    const std::vector<ged_obol_preserved_source> &sources,
    SoBRLSceneController *scene)
{
    if (!gedp || !scene || sources.empty())
	return 0;

    std::vector<std::string> current_paths =
	ged_obol_drawn_source_paths(gedp, NULL, -1, NULL);
    std::vector<ged_obol_preserved_source> replay_sources;
    for (const ged_obol_preserved_source &source : sources) {
	int already_current = 0;
	for (const std::string &path : current_paths) {
	    if (ged_obol_path_equal(source.path.c_str(), path.c_str())) {
		already_current = 1;
		break;
	    }
	}
	if (!already_current)
	    replay_sources.push_back(source);
    }

    return ged_obol_replay_preserved_sources(gedp, replay_sources, scene);
}

static int
ged_obol_observer_ensure(struct ged *gedp, struct ged_drawable *gdp)
{
    if (!gedp || !gdp)
	return 0;
    if (gdp->gd_obol_observer_token)
	return 1;

    gdp->gd_obol_observer_token = ged_draw_observer_add(gedp,
				  ged_obol_transaction_observer, NULL);
    return gdp->gd_obol_observer_token ? 1 : 0;
}

static int
ged_obol_register_progressive_provider(struct ged *gedp,
				       void *view_ctx,
				       ged_obol_attached_controller &entry)
{
    if (!gedp || !entry.view_controller)
	return 0;

    ged_obol_progressive_provider_data *data =
	new (std::nothrow) ged_obol_progressive_provider_data;
    if (!data)
	return 0;

    data->gedp = gedp;
    data->view_ctx = view_ctx;
    uint64_t token = entry.view_controller->registerProgressiveProvider(
			 ged_obol_progressive_advance_provider, data);
    if (!token) {
	delete data;
	return 0;
    }

    entry.progressive_provider_data = data;
    entry.progressive_provider_token = token;
    return 1;
}

static int
ged_draw_obol_attach_common(struct ged *gedp,
			    SoBRLSceneController *scene_controller,
			    BRLObolViewController *view_controller,
			    int sync_current_scene,
			    int owned_scene_controller,
			    int owned_view_controller)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp || !scene_controller)
	return 0;

    ged_draw_obol_scene_controller_detach(gedp);

    if (!ged_obol_observer_ensure(gedp, gdp))
	return 0;

    std::vector<ged_obol_attached_controller> *entries =
	ged_obol_attached_controllers(gdp, 1);
    if (!entries)
	return 0;

    ged_obol_attached_controller entry;
    entry.scene_controller = scene_controller;
    entry.view_controller = view_controller;
    entry.owned_scene_controller = owned_scene_controller ? 1 : 0;
    entry.owned_view_controller = owned_view_controller ? 1 : 0;
    entry.full_sync = 0;
    entry.use_attached_view_scope = 0;
    entry.render_endpoint_only = 0;
    entries->push_back(entry);
    ged_obol_primary_set(gdp, &entries->back());

    std::vector<ged_obol_preserved_source> preserved_sources;
    if (sync_current_scene)
	ged_obol_preserved_sources_take(gdp, preserved_sources);
    else
	ged_obol_preserved_sources_clear(gdp);

    if (sync_current_scene) {
	(void)ged_draw_obol_scene_sync_full_scene(gedp, NULL, 0,
		scene_controller);
	(void)ged_obol_replay_preserved_sources_not_current(gedp,
		preserved_sources, scene_controller);
    }
    if (sync_current_scene) {
	entries->back().full_sync = 1;
	ged_obol_primary_set(gdp, &entries->back());
    }

    return 1;
}

static int
ged_obol_bind_view_render_root(void *view_ctx,
			       SoBRLSceneController *shared_scene,
			       BRLObolViewController *view_controller)
{
    if (!shared_scene || !view_controller)
	return 0;

    SoNode *shared_root = shared_scene->getSceneRoot();
    if (!shared_root)
	return 0;

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
    if (local_root && local_root != shared_root)
	render_root->addChild(local_root);
    view_controller->setRenderSceneRoot(render_root);
    return 1;
}

static int
ged_draw_obol_attach_view_common(struct ged *gedp,
				 void *view_ctx,
				 SoBRLSceneController *scene_controller,
				 BRLObolViewController *view_controller,
				 int sync_current_scene)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp || !scene_controller || !view_controller)
	return 0;

    if (!ged_obol_observer_ensure(gedp, gdp))
	return 0;

    SoBRLSceneController *shared_scene =
	ged_draw_obol_scene_controller_ensure(gedp, sync_current_scene);
    if (!shared_scene ||
	!ged_obol_bind_view_render_root(view_ctx, shared_scene,
					view_controller))
	return 0;

    std::vector<ged_obol_attached_controller> *entries =
	ged_obol_attached_controllers(gdp, 1);
    if (!entries)
	return 0;

    for (size_t i = 0; i < entries->size(); i++) {
	ged_obol_attached_controller &entry = (*entries)[i];
	if (entry.view_ctx == view_ctx ||
	    entry.view_controller == view_controller ||
	    entry.scene_controller == scene_controller) {
	    ged_obol_delete_owned_attached_controller(entry);
	    entries->erase(entries->begin() +
			   static_cast<std::vector<ged_obol_attached_controller>::difference_type>(i));
	    break;
	}
    }

    ged_obol_attached_controller entry;
    entry.view_ctx = view_ctx;
    entry.scene_controller = view_controller->getSceneController();
    entry.view_controller = view_controller;
    entry.owned_scene_controller = 0;
    entry.owned_view_controller = 0;
    entry.full_sync = 0;
    entry.use_attached_view_scope = 1;
    entry.render_endpoint_only = 1;
    entry.render_shared_root_visible =
	ged_obol_view_scope_is_independent(view_ctx) ? 0 : 1;
    if (!ged_obol_register_progressive_provider(gedp, view_ctx, entry)) {
	ged_obol_delete_owned_attached_controller(entry);
	return 0;
    }
    entries->push_back(entry);

    ged_obol_primary_refresh_from_registry(gdp);

    return 1;
}

int
ged_draw_obol_scene_controller_attach(struct ged *gedp,
				      SoBRLSceneController *controller,
				      int sync_current_scene)
{
    return ged_draw_obol_attach_common(gedp, controller, NULL,
				       sync_current_scene, 0, 0);
}

int
ged_draw_obol_controller_attach(struct ged *gedp,
				BRLObolViewController *controller,
				int sync_current_scene)
{
    if (!controller)
	return 0;
    return ged_draw_obol_attach_common(gedp, controller->getSceneController(),
				       controller, sync_current_scene, 0, 0);
}

extern "C" int
ged_draw_obol_controller_attach_opaque(struct ged *gedp,
				       void *controller,
				       int sync_current_scene)
{
    return ged_draw_obol_controller_attach(gedp,
					   static_cast<BRLObolViewController *>(controller),
					   sync_current_scene);
}

int
ged_draw_obol_controller_attach_for_view(struct ged *gedp,
	void *view_ctx,
	BRLObolViewController *controller,
	int sync_current_scene)
{
    if (!controller)
	return 0;
    return ged_draw_obol_attach_view_common(gedp, view_ctx,
					    controller->getSceneController(), controller,
					    sync_current_scene);
}

extern "C" int
ged_draw_obol_controller_attach_opaque_for_view(struct ged *gedp,
	void *view_ctx,
	void *controller,
	int sync_current_scene)
{
    return ged_draw_obol_controller_attach_for_view(gedp, view_ctx,
	    static_cast<BRLObolViewController *>(controller),
	    sync_current_scene);
}

SoBRLSceneController *
ged_draw_obol_scene_controller_ensure(struct ged *gedp,
				      int sync_current_scene)
{
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (scene) {
	struct ged_drawable *gdp = ged_obol_gdp(gedp);
	if (sync_current_scene && gdp) {
	    std::vector<ged_obol_preserved_source> preserved_sources;
	    ged_obol_preserved_sources_take(gdp, preserved_sources);
	    if (!gdp->gd_obol_scene_controller_full_sync) {
		(void)ged_draw_obol_scene_sync_full_scene(gedp, NULL, 0,
			scene);
		gdp->gd_obol_scene_controller_full_sync = 1;
	    }
	    (void)ged_obol_replay_preserved_sources_not_current(gedp,
		    preserved_sources, scene);
	}
	return scene;
    }

    if (!ged_obol_gdp(gedp))
	return NULL;

    brlobol_init(NULL);

    SoSeparator *root = new SoSeparator;
    BRLObolViewController *owned_controller = new BRLObolViewController(root);
    SoBRLSceneController *owned_scene = owned_controller->getSceneController();
    if (!ged_draw_obol_attach_common(gedp, owned_scene, owned_controller,
				     sync_current_scene, 0, 1)) {
	delete owned_controller;
	return NULL;
    }

    return owned_scene;
}

void
ged_draw_obol_scene_controller_detach(struct ged *gedp)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    if (!gdp)
	return;

    SoBRLSceneController *preserve_scene =
	static_cast<SoBRLSceneController *>(gdp->gd_obol_scene_controller);
    if (preserve_scene) {
	std::vector<ged_obol_preserved_source> preserved_sources;
	ged_obol_collect_preserved_sources(preserve_scene, preserved_sources);
	ged_obol_preserved_sources_store(gdp, preserved_sources);
    }

    if (gdp->gd_obol_observer_token) {
	(void)ged_draw_observer_remove(gedp, gdp->gd_obol_observer_token);
	gdp->gd_obol_observer_token = 0;
    }

    ged_obol_attached_controllers_free(gdp);
}

void
ged_draw_obol_controller_detach(struct ged *gedp)
{
    ged_draw_obol_scene_controller_detach(gedp);
}

extern "C" void
ged_draw_obol_controller_detach_opaque(struct ged *gedp,
				       void *controller)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    BRLObolViewController *view_controller =
	static_cast<BRLObolViewController *>(controller);
    if (!gdp || !view_controller)
	return;

    std::vector<ged_obol_attached_controller> *entries =
	ged_obol_attached_controllers(gdp, 0);
    if (!entries)
	return;

    for (size_t i = 0; i < entries->size(); i++) {
	ged_obol_attached_controller &entry = (*entries)[i];
	if (entry.view_controller != view_controller)
	    continue;
	ged_obol_delete_owned_attached_controller(entry);
	entries->erase(entries->begin() +
		       static_cast<std::vector<ged_obol_attached_controller>::difference_type>(i));
	break;
    }

    if (entries->empty()) {
	delete entries;
	gdp->gd_obol_attached_controllers = NULL;
	ged_obol_primary_clear(gdp);
	if (gdp->gd_obol_observer_token) {
	    (void)ged_draw_observer_remove(gedp,
					   gdp->gd_obol_observer_token);
	    gdp->gd_obol_observer_token = 0;
	}
	return;
    }

    ged_obol_primary_refresh_from_registry(gdp);
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
