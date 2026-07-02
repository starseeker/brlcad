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
#include "brlobol/init.h"
#include "brlobol/lod_realization.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/scene_controller.h"
#include "brlobol/scene_group.h"
#include "brlobol/vlist_shape.h"
#include "brlobol/view_controller.h"
#include "brlobol/view_store.h"
#include "bg/line_layer.h"
#include "bu/hash.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "ged.h"
#include "ged/draw.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "rt/view.h"
#include "vmath.h"

#include "./bsg_ged_draw_private.h"
#include "./ged_private.h"

#include <algorithm>
#include <Inventor/nodes/SoSeparator.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <string>
#include <vector>

static const char ged_obol_group_intent_prefix[] = "ged-draw-group:";
static thread_local int ged_obol_source_summary_force_adapter = 0;

class ged_obol_source_summary_adapter_scope {
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
	highlighted(false),
	lineStyle(0),
	lineWidth(0),
	transparency(0.0f),
	materialColorValid(false),
	materialColor(1.0f, 1.0f, 1.0f),
	materialRevision(0)
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
    bool highlighted;
    int lineStyle;
    int lineWidth;
    float transparency;
    bool materialColorValid;
    SbColor materialColor;
    uint32_t materialRevision;
};

struct ged_obol_preserved_source {
    std::string path;
    int drawMode;
};

static struct ged_drawable *
ged_obol_gdp(struct ged *gedp)
{
    if (!gedp || !gedp->i)
	return NULL;
    return gedp->i->ged_gdp;
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

    if (!ged_obol_view_scope_is_independent(view_ctx))
	return std::string(path);

    std::string key("ged-view:");
    key += ged_obol_view_scope_name(view_ctx);
    key += ":";
    key += ged_obol_skip_leading_slash(path);
    return key;
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
    if (!ged_obol_view_scope_is_independent(view_ctx))
	return summary.instanceKey.getLength() == 0 ||
	    ged_obol_path_equal(summary.instanceKey.getString(),
		summary.path.getString());

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

static int
ged_obol_group_matches_view(struct ged *gedp,
	ged_draw_group_ref group,
	void *view_ctx)
{
    if (!gedp || ged_draw_group_ref_is_null(group))
	return 0;

    struct ged_draw_group_record group_record;
    if (!ged_draw_group_record_get(gedp, group, &group_record))
	return 0;
    return ged_draw_group_record_in_view(&group_record, view_ctx) ? 1 : 0;
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

    struct ged_draw_scene_display_summary display_summary;
    if (ged_draw_shape_ref_display_summary(gedp, record->ref,
	    &display_summary) &&
	    display_summary.valid) {
	state.visible = display_summary.visible ? true : false;
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

struct ged_obol_find_source_state_context {
    struct ged *gedp;
    void *viewCtx;
    const char *path;
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

    const int view_matched = ged_obol_group_matches_view(ctx->gedp,
	    record->group, ctx->viewCtx);
    if (!ctx->state.valid || view_matched || !ctx->state.viewMatched)
	ged_obol_source_state_from_record(ctx->state, ctx->gedp, record,
		view_matched);

    return view_matched ? 0 : 1;
}

static ged_obol_source_state
ged_obol_find_source_state(struct ged *gedp,
	void *view_ctx,
	const char *path)
{
    if (!gedp || !path || !path[0])
	return ged_obol_source_state();

    ged_obol_find_source_state_context ctx;
    ctx.gedp = gedp;
    ctx.viewCtx = view_ctx;
    ctx.path = path;
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

    if (scene->setGroupDisplayState(state.groupPath.c_str(),
	    state.groupVisible ? TRUE : FALSE,
	    FALSE,
	    FALSE,
	    0,
	    state.lineWidth,
	    state.groupTransparency,
	    FALSE,
	    SbColor(1.0f, 1.0f, 1.0f),
	    FALSE,
	    SbColor(1.0f, 1.0f, 1.0f),
	    0) > 0)
	changed = 1;

    if (source_instance_key && source_instance_key[0] &&
	    scene->moveDatabaseSourceInstanceToGroup(source_instance_key,
		state.groupPath.c_str()) > 0)
	changed = 1;
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

static int
ged_obol_path_has_descendant_in_list(const std::vector<std::string> &paths,
	const char *path)
{
    if (!path || !path[0])
	return 0;

    for (const std::string &candidate : paths) {
	if (ged_obol_path_equal(candidate.c_str(), path))
	    continue;
	if (ged_obol_path_has_prefix(candidate.c_str(), path))
	    return 1;
    }

    return 0;
}

static void
ged_obol_remove_shadowed_source_paths(std::vector<std::string> &paths)
{
    if (paths.size() < 2)
	return;

    std::vector<std::string> filtered;
    for (const std::string &path : paths) {
	if (!ged_obol_path_has_descendant_in_list(paths, path.c_str()))
	    ged_obol_append_unique_path(filtered, path.c_str());
    }
    paths.swap(filtered);
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

struct ged_obol_drawn_source_path_ctx {
    struct ged *gedp;
    void *viewCtx;
    int mode;
    const std::vector<std::string> *targets;
    std::vector<std::string> paths;
};

static int
ged_obol_drawn_source_path_cb(const struct ged_draw_shape_record *record,
	void *userdata)
{
    ged_obol_drawn_source_path_ctx *ctx =
	static_cast<ged_obol_drawn_source_path_ctx *>(userdata);
    if (!ctx || !record || !record->fullpath || !record->visible)
	return 1;
    if (ctx->mode >= 0 && record->draw_mode != ctx->mode)
	return 1;

    struct ged_draw_group_record group_record;
    if (!ged_draw_group_record_get(ctx->gedp, record->group, &group_record) ||
	    group_record.is_overlay ||
	    !group_record.visible ||
	    !ged_draw_group_record_in_view(&group_record, ctx->viewCtx))
	return 1;

    char *record_path = db_path_to_string(record->fullpath);
    if (!record_path)
	return 1;

    if (ged_obol_record_matches_any_target(record, record_path, ctx->targets))
	ged_obol_append_unique_path(ctx->paths, record_path);

    bu_free(record_path, "GED Obol drawn source path");
    return 1;
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
    ged_draw_foreach_shape_record(gedp, ged_obol_drawn_source_path_cb, &ctx);
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
ged_obol_apply_source_state(SoBRLSceneController *scene,
	const char *source_instance_key,
	const ged_obol_source_state &state)
{
    if (!scene || !source_instance_key || !source_instance_key[0] ||
	    !state.valid)
	return 0;

    return scene->setDatabaseSourceInstanceState(source_instance_key,
	    state.sourceRevision != 0 ? TRUE : FALSE,
	    state.sourceRevision,
	    state.inputsRevision,
	    state.visible ? TRUE : FALSE,
	    state.highlighted ? TRUE : FALSE,
	    state.lineStyle,
	    state.lineWidth,
	    state.transparency,
	    FALSE,
	    SbColor(1.0f, 1.0f, 1.0f),
	    state.materialColorValid ? TRUE : FALSE,
	    state.materialColor,
	    state.materialRevision);
}

static int
ged_obol_replace_path(struct ged *gedp,
	void *view_ctx,
	struct db_i *dbip,
	const char *path,
	int ged_draw_mode,
	uint32_t source_revision,
	SoBRLSceneController *scene)
{
    if (!dbip || !path || !path[0] || !scene)
	return 0;

    ged_obol_source_state source_state =
	ged_obol_find_source_state(gedp, view_ctx, path);
    const bool retained_state_valid = source_state.valid;
    const uint32_t retained_source_revision = source_state.sourceRevision;
    const uint32_t retained_inputs_revision = source_state.inputsRevision;
    if (source_state.valid) {
	if (source_state.sourceRevision != 0)
	    source_revision = source_state.sourceRevision;
    }

    const std::string instance_key =
	ged_obol_database_source_instance_key(view_ctx, path);
    SoBRLDatabaseSource *existing_source =
	scene->findDatabaseSourceInstance(instance_key.c_str());
    BRLObolDatabaseSourceSummary existing_summary;
    if (existing_source && existing_source->getSummary(existing_summary) &&
	    existing_summary.valid) {
	source_state.valid = true;
	source_state.sourceRevision = existing_summary.sourceRevision;
	source_state.inputsRevision = existing_summary.inputsRevision;
	source_state.visible = existing_summary.visible ? true : false;
	source_state.highlighted =
	    existing_summary.highlighted ? true : false;
	source_state.lineStyle = existing_summary.lineStyle;
	source_state.lineWidth = existing_summary.lineWidth;
	source_state.transparency =
	    static_cast<float>(existing_summary.transparency);
	source_state.materialColorValid =
	    existing_summary.materialColorValid ? true : false;
	source_state.materialColor = existing_summary.materialColor;
	source_state.materialRevision = existing_summary.materialRevision;
	if (retained_state_valid) {
	    source_state.sourceRevision = retained_source_revision;
	    source_state.inputsRevision = retained_inputs_revision;
	} else if (source_revision != 0) {
	    source_state.sourceRevision = source_revision;
	}
    }
    if (!source_state.groupValid && ged_draw_mode >= 0) {
	source_state.valid = true;
	source_state.groupValid = true;
	source_state.groupPath = ged_obol_skip_leading_slash(path);
	source_state.groupDrawMode =
	    ged_obol_lod_draw_mode_from_ged(ged_draw_mode);
	source_state.groupVisible = true;
	source_state.groupOverlay = false;
	source_state.groupTransparency = source_state.transparency;
    }
    int changed = scene->replaceDatabaseSourceInstance(instance_key.c_str(),
	    path, dbip,
	    ged_obol_database_draw_mode_from_ged(ged_draw_mode),
	    source_revision);
    if (ged_obol_apply_source_state(scene, instance_key.c_str(),
	    source_state))
	changed = 1;
    if (ged_obol_sync_group_state(scene, source_state,
	    instance_key.c_str()))
	changed = 1;
    return changed;
}

static int
ged_obol_replace_paths(struct db_i *dbip,
	const std::vector<std::string> &paths,
	int ged_draw_mode,
	uint32_t source_revision,
	struct ged *gedp,
	void *view_ctx,
	SoBRLSceneController *scene)
{
    if (!dbip || paths.empty() || !scene)
	return 0;

    int changed = 0;
    for (const std::string &path : paths) {
	if (ged_obol_replace_path(gedp, view_ctx, dbip, path.c_str(),
		ged_draw_mode, source_revision, scene) > 0)
	    changed = 1;
    }

    if (changed)
	scene->realizePending();
    return changed;
}

static void
ged_obol_collect_database_sources_matching(
	std::vector<std::string> &paths,
	SoBRLSceneController *scene,
	void *view_ctx,
	const char *target_path,
	size_t component_first_idx,
	int allow_path_prefix)
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
		!ged_obol_database_source_instance_in_scope(summary, view_ctx))
	    continue;

	const char *source_path = source->path.getValue().getString();
	if ((allow_path_prefix &&
		ged_obol_path_has_prefix(source_path, target_path)) ||
		ged_obol_path_has_component_name(source_path, target_path,
		    component_first_idx))
	    ged_obol_append_unique_path(paths, source_path);
    }
}

static std::vector<std::string>
ged_obol_matching_database_source_paths(
	SoBRLSceneController *scene,
	void *view_ctx,
	const std::vector<std::string> &targets,
	size_t component_first_idx,
	int allow_path_prefix)
{
    std::vector<std::string> paths;
    if (!scene)
	return paths;

    for (const std::string &target : targets)
	ged_obol_collect_database_sources_matching(paths, scene, view_ctx,
		target.c_str(), component_first_idx, allow_path_prefix);
    return paths;
}

static int
ged_obol_replace_matching_database_sources(
	struct ged *gedp,
	void *view_ctx,
	const std::vector<std::string> &targets,
	size_t component_first_idx,
	int allow_path_prefix,
	uint32_t source_revision,
	SoBRLSceneController *scene)
{
    if (!gedp || !gedp->dbip || !scene || targets.empty())
	return 0;

    std::vector<std::string> source_paths =
	ged_obol_matching_database_source_paths(scene, view_ctx, targets,
		component_first_idx, allow_path_prefix);
    ged_obol_remove_shadowed_source_paths(source_paths);
    if (source_paths.empty())
	return 0;

    int changed = 0;
    for (const std::string &source_path : source_paths) {
	if (ged_obol_replace_path(gedp, view_ctx, gedp->dbip,
		source_path.c_str(),
		ged_obol_drawn_path_mode(gedp, view_ctx, source_path.c_str()),
		source_revision, scene) > 0)
	    changed = 1;
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
	SoBRLSceneController *scene)
{
    if (!scene || targets.empty())
	return 0;

    std::vector<std::string> source_paths =
	ged_obol_matching_database_source_paths(scene, view_ctx, targets,
		component_first_idx, allow_path_prefix);
    if (source_paths.empty())
	return 0;

    for (const std::string &source_path : source_paths) {
	const std::string instance_key =
	    ged_obol_database_source_instance_key(view_ctx,
		    source_path.c_str());
	(void)scene->markDatabaseSourceInstanceStale(instance_key.c_str(),
		stale_reason);
    }
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
	return (result && result->status >= 0) ? 1 : 0;

    const int visible = ZERO(txn->value) ? 0 : 1;
    int handled = 0;
    for (const std::string &target : targets) {
	std::vector<std::string> source_paths;
	ged_obol_collect_database_sources_matching(source_paths, scene,
		txn->view,
		target.c_str(), 0, 1);
	for (const std::string &source_path : source_paths) {
	    const std::string instance_key =
		ged_obol_database_source_instance_key(txn->view,
			source_path.c_str());
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

    if (!handled && result && result->status >= 0)
	return 1;
    return handled;
}

static int
ged_obol_remove_paths(const std::vector<std::string> &paths,
	void *view_ctx,
	SoBRLSceneController *scene)
{
    if (paths.empty() || !scene)
	return 0;

    int changed = 0;
    for (const std::string &path : paths) {
	const std::string instance_key =
	    ged_obol_database_source_instance_key(view_ctx, path.c_str());
	int removed = scene->removeDatabaseSourceInstance(instance_key.c_str());
	if (removed > 0)
	    changed = 1;
    }

    if (ged_obol_prune_empty_groups(scene))
	changed = 1;

    if (changed)
	scene->realizePending();
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
	if (scene->removeDatabaseSourceInstance(instance_key.c_str()) > 0)
	    changed = 1;
    }

    if (ged_obol_prune_empty_groups(scene))
	changed = 1;

    if (changed)
	scene->realizePending();
    return changed;
}

static int
ged_obol_clear_database_sources_in_scope(SoBRLSceneController *scene,
	void *view_ctx)
{
    if (!scene)
	return 0;

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

int
ged_draw_obol_scene_controller_owned(struct ged *gedp)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    return (gdp && gdp->gd_obol_scene_controller &&
	    (gdp->gd_obol_scene_controller_owned ||
	     gdp->gd_obol_controller_owned)) ? 1 : 0;
}

extern "C" int
ged_draw_obol_scene_controller_ensure_owned(struct ged *gedp,
					    int sync_current_scene)
{
    SoBRLSceneController *scene =
	ged_draw_obol_scene_controller_ensure(gedp, sync_current_scene);
    return scene && ged_draw_obol_scene_controller_owned(gedp) ? 1 : 0;
}

extern "C" int
ged_draw_obol_scene_controller_full_synced(struct ged *gedp)
{
    struct ged_drawable *gdp = ged_obol_gdp(gedp);
    return (gdp && gdp->gd_obol_scene_controller &&
	    gdp->gd_obol_scene_controller_full_sync) ? 1 : 0;
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
    mat[0] = m[0][0];  mat[4] = m[0][1];  mat[8] = m[0][2];   mat[12] = m[0][3];
    mat[1] = m[1][0];  mat[5] = m[1][1];  mat[9] = m[1][2];   mat[13] = m[1][3];
    mat[2] = m[2][0];  mat[6] = m[2][1];  mat[10] = m[2][2];  mat[14] = m[2][3];
    mat[3] = m[3][0];  mat[7] = m[3][1];  mat[11] = m[3][2];  mat[15] = m[3][3];
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
    if (!gedp || !ged_draw_obol_scene_controller_owned(gedp))
	return NULL;
    return ged_draw_obol_controller(gedp);
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

extern "C" int
ged_draw_obol_view_context_feature_store_active(void *view_ctx)
{
    return ged_obol_view_controller_for_context(view_ctx) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_feature_exists(void *view_ctx, const char *name)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name)
	return 0;
    return controller->features().exists(name) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_feature_remove(void *view_ctx, const char *name)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name)
	return 0;
    return controller->features().remove(name) ? 1 : 0;
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
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name)
	return 0;

    BRLObolFeatureSummary obol_summary;
    if (!controller->features().summary(name, obol_summary))
	return 0;
    if (!obol_summary.exists)
	return 1;

    summary->exists = 1;
    summary->visible = obol_summary.visible ? 1 : 0;
    summary->child_count = obol_summary.childCount;
    summary->geometry_command_count = obol_summary.commandCount;
    summary->is_label =
	(obol_summary.kind == BRLObolFeatureKind::Labels ||
	 obol_summary.kind == BRLObolFeatureKind::HudLabel) ? 1 : 0;
    summary->is_transient_preview =
	(obol_summary.kind == BRLObolFeatureKind::EditPreview) ? 1 : 0;
    summary->is_overlay = (!summary->is_label &&
	    !summary->is_transient_preview) ? 1 : 0;

    BRLObolFeatureHandle handle = controller->features().find(name);
    BRLObolFeatureStyle style;
    if (handle.isValid() && controller->features().style(handle, style) &&
	    style.hasColor)
	ged_obol_rgb_from_color(style.color, summary->color);

    return 1;
}

extern "C" int
ged_draw_obol_view_context_line_layer_builder_replace(
	void *view_ctx,
	const char *name,
	int local,
	const struct bg_line_layer_builder *builder)
{
    BRLObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name)
	return 0;
    if (!builder || !bg_line_layer_builder_point_count(builder))
	return controller->features().remove(name) ? 1 : 1;

    BRLObolFeatureHandle handle =
	controller->features().publishLineLayerBuilder(name,
		local ? BRLObolFeatureScope::Local :
		BRLObolFeatureScope::Shared, builder);
    return handle.isValid() ? 1 : 0;
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
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name)
	return 0;

    size_t live_layers = 0;
    for (size_t i = 0; layers && i < layer_count; i++) {
	if (layers[i].points && layers[i].point_count)
	    live_layers++;
    }
    if (!layers || !live_layers)
	return controller->features().remove(name) ? 1 : 1;

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
	controller->features().publishLineLayers(name,
		local ? BRLObolFeatureScope::Local :
		BRLObolFeatureScope::Shared, obol_layers,
	style ? &obol_style : NULL);
    return handle.isValid() ? 1 : 0;
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
    if (removed > 0)
	controller->requestRender("overlay-remove");
    (void)gedp;
    return removed >= 0 ? 1 : 0;
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

    return scene->findDatabaseSource(path);
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

    const int source_count = scene->getDatabaseSourceCount();
    for (int i = 0; i < source_count; i++) {
	BRLObolDatabaseSourceSummary candidate;
	if (!scene->getDatabaseSourceSummary(i, candidate) ||
		!candidate.valid)
	    continue;
	if (ged_obol_path_equal(candidate.path.getString(), path)) {
	    summary = candidate;
	    return true;
	}
    }

    return false;
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

    const int count = scene->getSceneBoundsSummaryCount();
    for (int i = 0; i < count; i++) {
	BRLObolSceneBoundsSummary summary;
	if (!scene->getSceneBoundsSummary(i, summary) ||
		!summary.valid ||
		summary.nodeKind != BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE ||
		!summary.boundsValid ||
		summary.bounds.isEmpty())
	    continue;

	const SbVec3f &bmin = summary.bounds.getMin();
	const SbVec3f &bmax = summary.bounds.getMax();
	VMINMAX(*min, *max, bmin);
	VMINMAX(*min, *max, bmax);
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
	    !ged_draw_obol_scene_controller_owned(gedp))
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
	    !ged_draw_obol_scene_controller_owned(gedp))
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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
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
    if (!ged_obol_owned_database_source_for_path(gedp, path))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed = scene->markDatabaseSourceStale(path,
	    ged_obol_stale_reason_from_ged(reason));
    return changed >= 0 ? 1 : 0;
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

    out->valid = 1;
    out->database_path = source->path.getValue().getString();
    out->draw_mode = ged_obol_database_draw_mode_to_ged(summary.drawMode);
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
    out->realization_view_scale = (fastf_t)summary.realizationViewScale;
    out->realization_bot_threshold =
	(uint64_t)summary.realizationBotThreshold;
    out->realization_curve_scale = (fastf_t)summary.realizationCurveScale;
    out->realization_point_scale = (fastf_t)summary.realizationPointScale;
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
    BRLObolDatabaseSourceSummary summary;
    if (source->getSummary(summary) && summary.valid) {
	out->draw_size_valid = summary.drawSizeValid ? 1 : 0;
	out->draw_size = (fastf_t)summary.drawSize;
    }
    out->rt_mesh_lod = source->getMeshLod();
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
	struct rt_mesh_lod *lod)
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

    const int next_draw_mode =
	ged_obol_database_draw_mode_from_ged(record->draw_mode);
    int draw_mode_changed =
	scene->setDatabaseSourceDrawMode(path, next_draw_mode);
    if (draw_mode_changed < 0)
	return 0;

    int state_changed = scene->setDatabaseSourceState(path,
	    TRUE,
	    ged_obol_fold_revision(record->source_revision),
	    ged_obol_fold_revision(record->inputs_revision),
	    summary.visible,
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

    int material_policy_changed = scene->setDatabaseSourceMaterialPolicy(
	    path, ged_obol_material_policy_from_ged(record->material_policy));
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

    int realization_changed = scene->setDatabaseSourceRealizationState(path,
	    realization_status,
	    ged_obol_fold_revision(record->realized_source_revision),
	    ged_obol_fold_revision(record->realized_inputs_revision),
	    stale_reason,
	    realization_status == SoBRLDatabaseSource::FAILED ?
	    "GED source realization failed" : NULL);
    if (realization_changed < 0)
	return 0;

    int role_changed = scene->setDatabaseSourceRealizationRoleFlags(path,
	    record->realization_role_flags);
    if (role_changed < 0)
	return 0;

    const uint32_t clamped_bot_threshold =
	record->realization_bot_threshold > UINT32_MAX ? UINT32_MAX :
	(uint32_t)record->realization_bot_threshold;
    int policy_changed = scene->setDatabaseSourceRealizationViewPolicy(path,
	    record->realization_view_dependent ? TRUE : FALSE,
	    (float)record->realization_view_scale,
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
    if (!ged_obol_owned_database_source_for_path(gedp, path))
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

    const int changed = scene->setDatabaseSourceRealizationState(path,
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
    out->view_scale = (fastf_t)summary.realizationViewScale;
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
    if (!ged_obol_owned_database_source_for_path(gedp, path))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed = scene->setDatabaseSourceRealizationRoleFlags(path,
	    role_flags);
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_set_realization_view_policy_for_path(
	struct ged *gedp,
	const char *path,
	int view_dependent,
	fastf_t view_scale,
	uint64_t bot_threshold,
	fastf_t curve_scale,
	fastf_t point_scale)
{
    if (!ged_obol_owned_database_source_for_path(gedp, path))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const uint32_t clamped_bot_threshold =
	bot_threshold > UINT32_MAX ? UINT32_MAX : (uint32_t)bot_threshold;
    const int changed = scene->setDatabaseSourceRealizationViewPolicy(path,
	    view_dependent ? TRUE : FALSE,
	    (float)view_scale,
	    clamped_bot_threshold,
	    (float)curve_scale,
	    (float)point_scale);
    return changed >= 0 ? 1 : 0;
}

extern "C" int
ged_draw_obol_database_source_realize_for_path(struct ged *gedp,
					       const char *path)
{
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    (void)scene->realizePending();

    BRLObolDatabaseSourceSummary summary;
    if (!source->getSummary(summary) || !summary.valid)
	return 0;

    return summary.realizationStatus == SoBRLDatabaseSource::REALIZED &&
	!source->needsRealization() ? 1 : 0;
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
	    !ged_draw_obol_scene_controller_owned(gedp))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const uint32_t folded_revision =
	source_revision ? ged_obol_fold_revision(source_revision) : 0;
    const int changed = scene->replaceDatabaseSource(path, dbip,
	    ged_obol_database_draw_mode_from_ged(ged_draw_mode),
	    folded_revision);

    BRLObolDatabaseSourceSummary summary;
    if (ged_obol_database_source_controller_summary_for_path(scene, path,
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
		    (void)scene->moveDatabaseSourceToGroup(path,
			    owner_group_path.c_str());
		}
	    }
	}
    }

    return changed >= 0 ? 1 : 0;
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
    std::string old_owner_group_path;
    int have_source_summary = 0;
    if (ged_obol_database_source_controller_summary_for_path(scene, path,
	    source_summary)) {
	have_source_summary = 1;
	const std::string candidate_owner =
	    ged_obol_database_source_owner_group_path_from_summary(
		    source_summary);
	if (ged_obol_path_equal(candidate_owner.c_str(), path))
	    old_owner_group_path = candidate_owner;
    }

    const SbBool source_visible =
	(!have_source_summary || source_summary.visible) ? TRUE : FALSE;
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
	ged_obol_group_path_from_record_path(new_path);
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
    if (!scene || !scene->findDatabaseSource(source_path))
	return 0;

    const std::string target_group =
	ged_obol_group_path_from_record_path(group_path);
    if (target_group.empty())
	return 0;

    const int moved = scene->moveDatabaseSourceToGroup(source_path,
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

    const int summary_count = scene->getSceneTreeSummaryCount();
    for (int i = 0; i < summary_count; i++) {
	BRLObolSceneTreeSummary summary;
	if (!scene->getSceneTreeSummary(i, summary) ||
		!summary.valid ||
		(summary.nodeKind != BRLObolSceneTreeSummary::NODE_VLIST_SHAPE &&
		 summary.nodeKind != BRLObolSceneTreeSummary::NODE_MESH_SHAPE) ||
		summary.path.getLength() == 0)
	    continue;

	const char *shape_path = summary.path.getString();
	SoNode *node = scene->findShape(shape_path);
	if (!node)
	    continue;
	if (ged_obol_shape_node_is_database_realization(node))
	    continue;
	if (ged_obol_shape_node_is_auxiliary_record(node))
	    continue;
	if (skip_overlay_groups &&
		(ged_obol_shape_node_bool_field(node, "overlayIntent", FALSE) ||
		 ged_obol_group_path_is_overlay(scene,
		     summary.path.getString())))
	    continue;

	if (!(*cb)(gedp, shape_path, userdata))
	    return 0;
    }

    return 1;
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
		!ged_obol_database_source_instance_in_scope(summary, view_ctx))
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
    if (!gedp || !path_prefix || !path_prefix[0])
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
	if (mode >= 0 &&
		ged_obol_database_draw_mode_to_ged(summary.drawMode) != mode)
	    continue;
	for (const std::string &target : targets) {
	    if (ged_obol_path_has_component_name(summary.path.getString(),
		    target.c_str(), nonroot_only ? 1 : 0)) {
		ged_obol_append_unique_path(paths, summary.path.getString());
		break;
	    }
	}
    }

    return ged_obol_remove_paths(paths, NULL, scene);
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

    SbBool nextVisible = visible_valid ?
	(visible ? TRUE : FALSE) : summary.visible;
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

    int draw_mode_changed = scene->setDatabaseSourceDrawMode(path,
	    nextDrawMode);
    if (draw_mode_changed < 0)
	return 0;

    int changed = scene->setDatabaseSourceState(path,
	    FALSE,
	    summary.sourceRevision,
	    summary.inputsRevision,
	    nextVisible,
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
    const SbBool next_selected = shape->selected.getValue();
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

    int changed = scene->setDatabaseSourceDisplayPatch(path, patch);
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
	    !ged_draw_obol_scene_controller_owned(gedp))
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
    const int count = scene->getGroupDatabaseSourceCount(group_path.c_str());
    if (count < 0)
	return 0;

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
    if (!scene ||
	    !ged_obol_database_source_controller_summary_for_path(scene, path,
		summary)) {
	if (!source->getSummary(summary) || !summary.valid)
	    return 0;
    }
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
    out->highlighted = summary.highlighted ? 1 : 0;
    out->line_style = summary.lineStyle;
    out->line_width = summary.lineWidth;
    out->transparency = ged_obol_reported_transparency(summary.transparency);
    out->draw_mode = exact_draw_mode;
    out->material_valid = 1;
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
    if (!source->getSummary(source_summary) || !source_summary.valid)
	return 0;

    out->valid = 1;
    out->draw_mode_valid = 1;
    out->draw_mode =
	ged_obol_database_draw_mode_to_ged(source_summary.drawMode);
    out->line_style = source_summary.lineStyle;

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
	ged_obol_mat_from_sbmatrix(display_summary.drawMatrix, out->draw_mat);
	break;
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
	if (shape->point.getNum() > 0 || shape->command.getNum() > 0)
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
    shape->drawMode = source->drawMode.getValue() == SoBRLDatabaseSource::SHADED ?
	BRLOBOL_LOD_DRAW_SHADED : BRLOBOL_LOD_DRAW_WIRE;
    shape->recordRole = "database";
    shape->geometryKind = "line-set";
    shape->visible = source->visible.getValue();
    shape->highlighted = source->highlighted.getValue();
    shape->lineStyle = source->lineStyle.getValue();
    shape->lineWidth = source->lineWidth.getValue();
    shape->transparency = source->transparency.getValue();
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
    return (shape && (shape->annotationPoint.getNum() > 0 ||
	    shape->annotationSegmentKind.getNum() > 0)) ? 1 : 0;
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
    shape->drawMode = source->drawMode.getValue() == SoBRLDatabaseSource::SHADED ?
	BRLOBOL_LOD_DRAW_SHADED : BRLOBOL_LOD_DRAW_WIRE;
    shape->recordRole = "database";
    shape->geometryKind = "surface";
    shape->visible = source->visible.getValue();
    shape->highlighted = source->highlighted.getValue();
    shape->lineStyle = source->lineStyle.getValue();
    shape->lineWidth = source->lineWidth.getValue();
    shape->transparency = source->transparency.getValue();
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
	const char *kind = shape->geometryKind.getValue().getString();
	out->geometry_name = (kind && BU_STR_EQUAL(kind, "annotation")) ?
	    "annotation" : "line-set";
	out->point_count = static_cast<size_t>(shape->point.getNum());
	out->index_count = 0;
	return 1;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	out->geometry_name = "indexed-face-set";
	out->point_count = static_cast<size_t>(shape->point.getNum());
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

    out->valid = 1;
    out->point_count = static_cast<size_t>(shape->point.getNum());
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
    if (index >= static_cast<size_t>(shape->point.getNum()))
	return 0;

    const SbVec3f &point = shape->point[static_cast<int>(index)];
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
    if (coord_index >= static_cast<size_t>(shape->coordIndex.getNum()))
	return 0;
    *out = shape->coordIndex[static_cast<int>(coord_index)];
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
    out->valid = 1;
    out->point_count = static_cast<size_t>(shape->point.getNum());
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
    if (index >= static_cast<size_t>(shape->point.getNum()))
	return 0;

    double precise[3];
    if (shape->getPrecisePoint(static_cast<int>(index), precise)) {
	out[0] = precise[0];
	out[1] = precise[1];
	out[2] = precise[2];
	return 1;
    }

    const SbVec3f &point = shape->point[static_cast<int>(index)];
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
    if (index >= static_cast<size_t>(shape->command.getNum()))
	return 0;

    *out = ged_obol_vlist_command_to_ged(
	    shape->command[static_cast<int>(index)]);
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
    if (!shape || shape->point.getNum() <= 0)
	return 0;

    const SbVec3f &point = shape->point[shape->point.getNum() - 1];
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
    SoBRLVListShape *shape =
	ged_obol_owned_annotation_vlist_shape_for_path(gedp, path);
    if (!shape)
	return 0;

    out->valid = 1;
    out->point_count = static_cast<size_t>(shape->point.getNum());
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
    if (!shape || index >= static_cast<size_t>(shape->point.getNum()))
	return 0;

    double precise[3];
    if (shape->getPrecisePoint(static_cast<int>(index), precise)) {
	out[0] = precise[0];
	out[1] = precise[1];
	out[2] = precise[2];
	return 1;
    }

    const SbVec3f &point = shape->point[static_cast<int>(index)];
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
    if (!shape || index >= static_cast<size_t>(shape->command.getNum()))
	return 0;

    *out = ged_obol_vlist_command_to_ged(
	    shape->command[static_cast<int>(index)]);
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
    if (!shape)
	return 0;

    int material_draw_mode =
	ged_obol_lod_draw_mode_to_ged(shape->drawMode.getValue());
    BRLObolDatabaseSourceSummary source_summary;
    if (source->getSummary(source_summary) && source_summary.valid)
	material_draw_mode =
	    ged_obol_database_source_exact_draw_mode_to_ged(gedp,
		    source_summary, source);

    const size_t triangle_count =
	static_cast<size_t>(shape->getTriangleCount());
    out->valid = 1;
    out->point_count = static_cast<size_t>(shape->point.getNum());
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
    if (!shape || index >= static_cast<size_t>(shape->point.getNum()))
	return 0;

    const SbVec3f &point = shape->point[static_cast<int>(index)];
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
    if (coord_index >= static_cast<size_t>(shape->coordIndex.getNum()))
	return 0;
    *out = shape->coordIndex[static_cast<int>(coord_index)];
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
    SoBRLVListShape *shape = ged_obol_owned_vlist_shape_for_path(gedp, path);
    if (!shape)
	return 0;

    shape->setLineSet(NULL, NULL, 0);
    return 1;
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

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    SoBRLVListShape *shape =
	ged_obol_owned_vlist_shape_for_source(source, path, 1);
    if (!shape)
	return 0;

    if (point_count == 0) {
	shape->setLineSet(NULL, NULL, 0);
	return 1;
    }

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

    shape->sourceType = "line-set";
    shape->geometryKind = "line";
    shape->setLineSet(obol_points.data(), obol_commands.data(),
	    static_cast<int>(point_count));
    ged_obol_vlist_shape_set_precise_points(shape, points, point_count);
    SoBRLMeshShape *mesh_shape = source->getRealizedMesh();
    if (mesh_shape)
	mesh_shape->setIndexedTriangles(NULL, 0, NULL, 0);
    return 1;
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

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    SoBRLVListShape *shape =
	ged_obol_owned_vlist_shape_for_source(source, path, 1);
    if (!shape)
	return 0;

    if (point_count == 0) {
	shape->setLineSet(NULL, NULL, 0);
	return 1;
    }

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

    shape->sourceType = "annotation";
    shape->geometryKind = "annotation";
    shape->setLineSet(obol_points.data(), obol_commands.data(),
	    static_cast<int>(point_count));
    ged_obol_vlist_shape_set_precise_points(shape, points, point_count);
    SoBRLMeshShape *mesh_shape = source->getRealizedMesh();
    if (mesh_shape)
	mesh_shape->setIndexedTriangles(NULL, 0, NULL, 0);
    return 1;
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
	    annotation_point_count > static_cast<size_t>(INT_MAX))
	return 0;

    if (!ged_draw_obol_database_source_publish_annotation_line_set_for_path(
	    gedp, path, line_points, line_commands, line_point_count))
	return 0;

    SoBRLVListShape *shape = ged_obol_owned_vlist_shape_for_path(gedp, path);
    if (!shape)
	return 0;

    shape->sourceType = "annotation";
    shape->geometryKind = "annotation";
    shape->annotationBasePoint = base_point ? SbVec3f(
	    static_cast<float>(base_point[0]),
	    static_cast<float>(base_point[1]),
	    static_cast<float>(base_point[2])) :
	SbVec3f(0.0f, 0.0f, 0.0f);

    if (annotation_point_count && annotation_points) {
	std::vector<double> precise_points(annotation_point_count * 3);
	for (size_t i = 0; i < annotation_point_count; i++) {
	    const size_t offset = i * 3;
	    precise_points[offset + 0] = annotation_points[i][0];
	    precise_points[offset + 1] = annotation_points[i][1];
	    precise_points[offset + 2] = annotation_points[i][2];
	}
	shape->setPreciseAnnotationPoints(precise_points.data(),
		static_cast<int>(annotation_point_count));
    } else {
	shape->setPreciseAnnotationPoints(NULL, 0);
    }

    shape->annotationSegmentKind.setNum(static_cast<int>(segment_count));
    shape->annotationSegmentStart.setNum(static_cast<int>(segment_count));
    shape->annotationSegmentEnd.setNum(static_cast<int>(segment_count));
    shape->annotationTextRefPoint.setNum(static_cast<int>(segment_count));
    shape->annotationText.setNum(static_cast<int>(segment_count));
    shape->annotationSegmentTextValid.setNum(static_cast<int>(segment_count));
    for (size_t i = 0; i < segment_count; i++) {
	const struct ged_draw_obol_annotation_segment *seg = &segments[i];
	const int idx = static_cast<int>(i);
	shape->annotationSegmentKind.set1Value(idx, seg->kind);
	shape->annotationSegmentStart.set1Value(idx, seg->line_start);
	shape->annotationSegmentEnd.set1Value(idx, seg->line_end);
	shape->annotationTextRefPoint.set1Value(idx, seg->text_ref_point);
	shape->annotationText.set1Value(idx,
		(seg->text && seg->text[0]) ? seg->text : "");
	shape->annotationSegmentTextValid.set1Value(idx,
		(seg->kind == GED_DRAW_OBOL_ANNOTATION_SEGMENT_TEXT &&
		 seg->text && seg->text[0]) ? TRUE : FALSE);
    }

    return 1;
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

    out->valid = 1;
    out->point_count = static_cast<size_t>(shape->annotationPoint.getNum());
    out->segment_count =
	static_cast<size_t>(shape->annotationSegmentKind.getNum());
    out->cache_identity =
	ged_obol_hash_sb_string(shape->cacheIdentity.getValue());
    out->source_identity =
	ged_obol_hash_sb_string(shape->sourceIdentity.getValue());

    for (int i = 0; i < shape->annotationSegmentKind.getNum(); i++) {
	const int kind_value = shape->annotationSegmentKind[i];
	if (!out->line_segment_valid &&
		kind_value == GED_DRAW_OBOL_ANNOTATION_SEGMENT_LINE) {
	    out->line_segment_valid = 1;
	    out->line_start = (i < shape->annotationSegmentStart.getNum()) ?
		shape->annotationSegmentStart[i] : 0;
	    out->line_end = (i < shape->annotationSegmentEnd.getNum()) ?
		shape->annotationSegmentEnd[i] : 0;
	}
	if (!out->text_segment_valid &&
		kind_value == GED_DRAW_OBOL_ANNOTATION_SEGMENT_TEXT) {
	    out->text_segment_valid = 1;
	    out->text_ref_point = (i < shape->annotationTextRefPoint.getNum()) ?
		shape->annotationTextRefPoint[i] : 0;
	    out->text = (i < shape->annotationText.getNum()) ?
		shape->annotationText[i].getString() : "";
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
    if (!shape || index >= static_cast<size_t>(shape->annotationPoint.getNum()))
	return 0;

    double precise[3];
    if (shape->getPreciseAnnotationPoint(static_cast<int>(index), precise)) {
	out[0] = precise[0];
	out[1] = precise[1];
	out[2] = precise[2];
	return 1;
    }

    const SbVec3f &point = shape->annotationPoint[static_cast<int>(index)];
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
    if (!ged_obol_owned_database_source_for_path(gedp, path))
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

    const int changed = scene->publishDatabaseSourceAuxiliaryLineSet(path,
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
    if (!ged_obol_owned_database_source_for_path(gedp, path))
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

    const int changed = scene->publishDatabaseSourceAuxiliarySourceLineSet(
	    path,
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
    if (!ged_obol_owned_database_source_for_path(gedp, path))
	return 0;

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int cleared = scene->clearDatabaseSourceAuxiliaryShapes(path);
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

    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    SoBRLVListShape *shape =
	ged_obol_owned_vlist_shape_for_source(source, path, 1);
    if (!shape)
	return 0;

    if (point_count == 0) {
	shape->setLineSet(NULL, NULL, 0);
	return 1;
    }

    std::vector<SbVec3f> obol_points;
    std::vector<int32_t> obol_commands;
    obol_points.reserve(point_count);
    obol_commands.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	obol_points.push_back(SbVec3f(
		static_cast<float>(points[i][0]),
		static_cast<float>(points[i][1]),
		static_cast<float>(points[i][2])));
	obol_commands.push_back(SoBRLVListShape::POINT);
    }

    shape->sourceType = "point-set";
    shape->geometryKind = "point";
    shape->setLineSet(obol_points.data(), obol_commands.data(),
	    static_cast<int>(point_count));
    ged_obol_vlist_shape_set_precise_points(shape, points, point_count);
    SoBRLMeshShape *mesh_shape = source->getRealizedMesh();
    if (mesh_shape)
	mesh_shape->setIndexedTriangles(NULL, 0, NULL, 0);
    return 1;
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
    SoBRLMeshShape *shape = ged_obol_owned_mesh_shape_for_path(gedp,
	    path, 0);
    if (!shape)
	return 0;

    shape->setIndexedTriangles(NULL, 0, NULL, 0);
    return 1;
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
    SoBRLDatabaseSource *source =
	ged_obol_owned_database_source_for_path(gedp, path);
    if (!source)
	return 0;

    if (!point_count || !index_count) {
	SoBRLMeshShape *shape =
	    ged_obol_owned_mesh_shape_for_source(source, path, 0);
	if (!shape)
	    return 0;
	shape->setIndexedTriangles(NULL, 0, NULL, 0);
	return 1;
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

    SoBRLMeshShape *shape =
	ged_obol_owned_mesh_shape_for_source(source, path, 1);
    if (!shape)
	return 0;

    SoBRLVListShape *vlist_shape = source->getRealizedShape();
    if (vlist_shape)
	vlist_shape->setLineSet(NULL, NULL, 0);
    shape->setIndexedTriangles(obol_points.data(),
	    static_cast<int>(obol_points.size()),
	    triangles.data(),
	    static_cast<int>(triangles.size()));
    return 1;
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed = scene->setDatabaseSourcePlacementState(path,
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

    SbBool nextDrawMatrixValid = summary.drawMatrixValid ? TRUE : FALSE;
    SbMatrix nextDrawMatrix = summary.drawMatrix;
    if (!nextDrawMatrixValid && shape->drawCenterValid.getValue() &&
	    shape->drawSizeValid.getValue()) {
	nextDrawMatrixValid = TRUE;
	nextDrawMatrix = SbMatrix::identity();
    }

    (void)scene->setDatabaseSourcePlacementState(path,
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

    const int changed = scene->setDatabaseSourcePlacementState(path,
	    nextDrawMatrixValid, nextDrawMatrix, nextDrawCenterValid,
	    nextDrawCenter, nextDrawSizeValid, nextDrawSize);
    return changed >= 0 ? 1 : 0;
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

    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed = scene->setDatabaseSourceBoundsState(path, TRUE,
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
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (!scene)
	return 0;

    const int changed = scene->setDatabaseSourceDisplayName(path,
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
    if (ged_obol_vlist_shape_is_annotation(annotation_shape) &&
	    (annotation_shape->point.getNum() > 0 ||
	     ged_obol_vlist_shape_has_annotation_record(annotation_shape))) {
	out->valid = 1;
	out->geometry_name = "annotation";
	out->point_count =
	    static_cast<size_t>(annotation_shape->point.getNum());
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

    return 0;
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

int
ged_draw_obol_scene_sync_full_scene(struct ged *gedp,
	void *view_ctx,
	uint32_t source_revision,
	SoBRLSceneController *controller)
{
    SoBRLSceneController *scene = controller ?
	controller : ged_draw_obol_scene_controller(gedp);
    if (!gedp || !scene)
	return 0;

    if (!source_revision)
	source_revision = ged_obol_fold_revision(ged_draw_scene_revision(gedp));

    struct db_i *dbip = gedp->dbip;
    std::vector<std::string> paths =
	ged_obol_drawn_source_paths(gedp, view_ctx, -1, NULL);

    int changed = ged_obol_clear_database_sources_in_scope(scene, view_ctx);
    if (!dbip || paths.empty()) {
	if (changed)
	    scene->realizePending();
	return changed ? 1 : 0;
    }

    for (const std::string &path : paths) {
	int ged_draw_mode = ged_obol_drawn_path_mode(gedp, view_ctx,
		path.c_str());
	if (ged_obol_replace_path(gedp, view_ctx, dbip, path.c_str(),
		ged_draw_mode, source_revision, scene) > 0)
	    changed = 1;
    }

    if (changed)
	scene->realizePending();
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
		    0, 1);
	if (source_paths.empty())
	    return (result && result->status > 0) ? 1 : 0;
	(void)ged_obol_remove_paths(source_paths, view_ctx, scene);
	return 1;
    }

    if (txn->redraw) {
	int handled = ged_obol_replace_matching_database_sources(gedp,
		view_ctx, targets, 0, 1, source_revision, scene);
	if (handled)
	    return 1;
    }

    const ged_draw_stale_reason stale_reason = txn->stale_reason ?
	txn->stale_reason : GED_DRAW_STALE_SOURCE_CHANGED;
    return ged_obol_mark_matching_database_sources_stale(view_ctx, targets,
	    0, 1,
	    ged_obol_stale_reason_from_ged(stale_reason), scene);
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
		0);
    if (source_paths.empty())
	return (result && result->status > 0) ? 1 : 0;
    (void)ged_obol_remove_paths(source_paths, txn->view, scene);
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
		view_ctx, targets, 0, 1, source_revision, scene);
	if (refreshed)
	    return 1;
    }

    const ged_draw_stale_reason stale_reason = txn->stale_reason ?
	txn->stale_reason : GED_DRAW_STALE_SOURCE_CHANGED;
    int handled = ged_obol_mark_matching_database_sources_stale(view_ctx,
	    targets, 0, 1,
	    ged_obol_stale_reason_from_ged(stale_reason), scene);
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
    for (const std::string &target : targets) {
	std::vector<std::string> matching_sources;
	ged_obol_collect_database_sources_matching(matching_sources, scene,
		view_ctx, target.c_str(), 0, 1);
	if (!matching_sources.empty()) {
	    if (ged_obol_replace_matching_database_sources(gedp, view_ctx,
		    matching_sources, 0, 1, source_revision, scene))
		handled = 1;
	    continue;
	}
	std::vector<std::string> one_target;
	ged_obol_append_unique_path(one_target, target.c_str());
	std::vector<std::string> target_sources =
	    ged_obol_drawn_source_paths(gedp, view_ctx, -1, &one_target);
	if (!target_sources.empty()) {
	    if (ged_obol_replace_paths(gedp->dbip, target_sources,
		    ged_obol_drawn_path_mode(gedp, view_ctx, target.c_str()),
		    source_revision, gedp, view_ctx, scene))
		handled = 1;
	    continue;
	}
	if (ged_obol_replace_path(gedp, view_ctx, gedp->dbip,
		target.c_str(),
		ged_obol_drawn_path_mode(gedp, view_ctx, target.c_str()),
		source_revision, scene) > 0)
	    handled = 1;
    }

    if (handled)
	scene->realizePending();
    if (!handled && result && result->status >= 0)
	return 1;
    return handled;
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
	case GED_DRAW_TXN_DRAW:
	{
	    if (result && result->affected_shapes <= 0 &&
		    result->affected_groups <= 0)
		break;
	    std::vector<std::string> paths =
		ged_obol_transaction_paths(txn, result);
	    if (paths.empty()) {
		changed = ged_draw_obol_scene_sync_full_scene(gedp,
			view_ctx, source_revision, scene);
	    } else {
		std::vector<std::string> source_paths =
		    ged_obol_drawn_source_paths(gedp, view_ctx, -1, &paths);
		if (!source_paths.empty()) {
		    std::vector<std::string> shadowed_targets =
			ged_obol_shadowed_target_source_paths(paths,
				source_paths);
		    if (ged_obol_remove_paths(shadowed_targets, view_ctx,
			    scene))
			changed = 1;
		    if (ged_obol_replace_paths(gedp->dbip, source_paths,
			    ged_obol_transaction_ged_draw_mode(gedp, txn),
			    source_revision, gedp, view_ctx, scene))
			changed = 1;
		} else {
		    changed = ged_obol_replace_paths(gedp->dbip, paths,
			    ged_obol_transaction_ged_draw_mode(gedp, txn),
			    source_revision, gedp, view_ctx, scene);
		}
	    }
	    break;
	}
	case GED_DRAW_TXN_ERASE:
	{
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
				view_ctx))
			continue;

		    const char *source_path = summary.path.getString();
		    for (const std::string &target : paths) {
			if (ged_obol_path_has_prefix(source_path,
				target.c_str()) ||
				ged_obol_path_has_component_name(source_path,
				    target.c_str(), 0) ||
				(!ged_obol_path_equal(source_path,
				    target.c_str()) &&
				    ged_obol_path_has_prefix(target.c_str(),
					source_path))) {
			    ged_obol_append_database_source_instance_key(
				    matching_instance_keys, summary);
			    break;
			}
		    }
		}
		changed = matching_instance_keys.empty() ?
		    ged_obol_remove_paths(paths, view_ctx, scene) :
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
	case GED_DRAW_TXN_MATERIAL_CHANGED:
	case GED_DRAW_TXN_REFRESH_MATERIAL_COLORS:
	    changed = (result && result->status >= 0) ? 1 : 0;
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
	case GED_DRAW_TXN_HIGHLIGHT:
	case GED_DRAW_TXN_NONE:
	default:
	    break;
    }

    return changed ? 1 : 0;
}

extern "C" int
ged_draw_obol_scene_sync_attached_transaction(
	struct ged *gedp,
	const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result)
{
    if (!ged_draw_obol_scene_controller(gedp))
	return 0;
    return ged_draw_obol_scene_sync_transaction(gedp, txn, result, NULL);
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
    (void)ged_draw_obol_scene_sync_transaction(gedp, txn, result, NULL);
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
    for (const ged_obol_preserved_source &source : sources) {
	if (ged_obol_replace_path(gedp, NULL, gedp->dbip,
		source.path.c_str(), source.drawMode, source_revision,
		scene) > 0)
	    changed = 1;
    }

    if (changed)
	scene->realizePending();
    return changed;
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

    gdp->gd_obol_scene_controller = scene_controller;
    gdp->gd_obol_controller = view_controller;
    gdp->gd_obol_scene_controller_owned = owned_scene_controller ? 1 : 0;
    gdp->gd_obol_controller_owned = owned_view_controller ? 1 : 0;
    gdp->gd_obol_scene_controller_full_sync = 0;
    gdp->gd_obol_observer_token = ged_draw_observer_add(gedp,
	    ged_obol_transaction_observer, NULL);
    if (!gdp->gd_obol_observer_token) {
	gdp->gd_obol_scene_controller = NULL;
	gdp->gd_obol_controller = NULL;
	gdp->gd_obol_scene_controller_owned = 0;
	gdp->gd_obol_controller_owned = 0;
	gdp->gd_obol_scene_controller_full_sync = 0;
	return 0;
    }

    std::vector<ged_obol_preserved_source> preserved_sources;
    if (sync_current_scene)
	ged_obol_preserved_sources_take(gdp, preserved_sources);
    else
	ged_obol_preserved_sources_clear(gdp);

    int replayed_preserved = 0;
    if (sync_current_scene && !preserved_sources.empty())
	replayed_preserved = ged_obol_replay_preserved_sources(gedp,
		preserved_sources, scene_controller);

    if (sync_current_scene && !replayed_preserved)
	(void)ged_draw_obol_scene_sync_full_scene(gedp, NULL, 0,
		scene_controller);
    if (sync_current_scene)
	gdp->gd_obol_scene_controller_full_sync = 1;

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

SoBRLSceneController *
ged_draw_obol_scene_controller_ensure(struct ged *gedp,
	int sync_current_scene)
{
    SoBRLSceneController *scene = ged_draw_obol_scene_controller(gedp);
    if (scene) {
	struct ged_drawable *gdp = ged_obol_gdp(gedp);
	if (sync_current_scene && gdp &&
		!gdp->gd_obol_scene_controller_full_sync) {
	    (void)ged_draw_obol_scene_sync_full_scene(gedp, NULL, 0, scene);
	    gdp->gd_obol_scene_controller_full_sync = 1;
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

    if (gdp->gd_obol_observer_token) {
	(void)ged_draw_observer_remove(gedp, gdp->gd_obol_observer_token);
	gdp->gd_obol_observer_token = 0;
    }

    BRLObolViewController *owned_controller =
	gdp->gd_obol_controller_owned ?
	static_cast<BRLObolViewController *>(gdp->gd_obol_controller) :
	NULL;
    SoBRLSceneController *owned_scene =
	(!owned_controller && gdp->gd_obol_scene_controller_owned) ?
	static_cast<SoBRLSceneController *>(gdp->gd_obol_scene_controller) :
	NULL;
    SoBRLSceneController *preserve_scene =
	(gdp->gd_obol_scene_controller_owned || gdp->gd_obol_controller_owned) ?
	static_cast<SoBRLSceneController *>(gdp->gd_obol_scene_controller) :
	NULL;
    if (preserve_scene) {
	std::vector<ged_obol_preserved_source> preserved_sources;
	ged_obol_collect_preserved_sources(preserve_scene, preserved_sources);
	ged_obol_preserved_sources_store(gdp, preserved_sources);
    }
    gdp->gd_obol_scene_controller = NULL;
    gdp->gd_obol_controller = NULL;
    gdp->gd_obol_scene_controller_owned = 0;
    gdp->gd_obol_controller_owned = 0;
    gdp->gd_obol_scene_controller_full_sync = 0;
    if (owned_controller)
	delete owned_controller;
    else if (owned_scene)
	delete owned_scene;
}

void
ged_draw_obol_controller_detach(struct ged *gedp)
{
    ged_draw_obol_scene_controller_detach(gedp);
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
