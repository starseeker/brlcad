/*              Q G O B O L D R A W S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolDrawSync.cpp */

#include "common.h"

#include "QgLegacyViewContext.h"
#include "QgObolDatabaseSyncPrivate.h"
#include "QgObolDrawSyncPrivate.h"
#include "brlobol/database_source.h"
#include "brlobol/lod_realization.h"
#include "brlobol/view_controller.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "ged/draw.h"
#include "qtcad/QgLegacyView.h"
#include "qtcad/QgView.h"
#include "rt/db_fullpath.h"

#include <algorithm>
#include <stdint.h>
#include <string.h>
#include <string>
#include <vector>

static const char qg_obol_ged_group_intent_prefix[] = "ged-draw-group:";

struct qg_obol_ged_source_state {
    qg_obol_ged_source_state(void) :
	valid(false),
	viewMatched(false),
	groupValid(false),
	groupDrawMode(BRLOBOL_LOD_DRAW_WIRE),
	groupVisible(true),
	groupOverlay(false),
	groupTransparency(0.0f),
	sourceRevision(0),
	inputsRevision(0),
	drawMode(QG_OBOL_DATABASE_WIREFRAME),
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

static void
qg_obol_append_unique_path(std::vector<std::string> &paths, const char *path)
{
    if (!path || !path[0])
	return;
    std::string spath(path);
    if (std::find(paths.begin(), paths.end(), spath) == paths.end())
	paths.push_back(spath);
}

static void
qg_obol_append_unique_paths_from_words(std::vector<std::string> &paths,
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
	qg_obol_append_unique_path(paths, path.c_str());
	if (end == std::string::npos)
	    break;
	pos = end + 1;
    }
}

static std::vector<std::string>
qg_obol_ged_transaction_paths(const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result)
{
    std::vector<std::string> paths;
    if (result && BU_VLS_IS_INITIALIZED(&result->names) &&
	    bu_vls_strlen(&result->names)) {
	qg_obol_append_unique_paths_from_words(paths,
		bu_vls_cstr(&result->names));
	if (!paths.empty())
	    return paths;
    }

    if (txn && txn->path)
	qg_obol_append_unique_path(paths, txn->path);

    int path_count = (txn && txn->paths && txn->path_count > 0) ?
	txn->path_count : 0;
    for (int i = 0; i < path_count; i++)
	qg_obol_append_unique_path(paths, txn->paths[i]);

    return paths;
}

static uint32_t
qg_obol_fold_revision(uint64_t revision)
{
    if (!revision)
	return 0;

    uint32_t folded = static_cast<uint32_t>(revision ^ (revision >> 32));
    return folded ? folded : 1;
}

static uint32_t
qg_obol_ged_source_revision(const struct ged_draw_transaction_result *result)
{
    return qg_obol_fold_revision(result ? result->scene_revision_after : 0);
}

static int
qg_obol_draw_mode_from_ged(int mode)
{
    if (mode == GED_DRAW_MODE_SHADED ||
	    mode == GED_DRAW_MODE_SHADED_BOTS)
	return QG_OBOL_DATABASE_SHADED;
    return QG_OBOL_DATABASE_WIREFRAME;
}

static int
qg_obol_lod_draw_mode_from_ged(int mode)
{
    if (mode == GED_DRAW_MODE_SHADED ||
	    mode == GED_DRAW_MODE_SHADED_BOTS)
	return BRLOBOL_LOD_DRAW_SHADED;
    if (mode == GED_DRAW_MODE_EVAL_POINTS)
	return BRLOBOL_LOD_DRAW_POINTS;
    return BRLOBOL_LOD_DRAW_WIRE;
}

static int
qg_obol_database_source_draw_mode(int drawMode)
{
    return drawMode == QG_OBOL_DATABASE_SHADED ?
	SoBRLDatabaseSource::SHADED : SoBRLDatabaseSource::WIREFRAME;
}

static int
qg_obol_ged_transaction_draw_mode(struct ged *gedp,
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
    return qg_obol_draw_mode_from_ged(mode);
}

static int
qg_obol_drawn_path_mode(struct ged *gedp, void *view_ctx,
	const char *path)
{
    if (ged_draw_path_state(gedp, view_ctx, path,
	    GED_DRAW_MODE_SHADED_BOTS) > 0 ||
	    ged_draw_path_state(gedp, view_ctx, path,
		GED_DRAW_MODE_SHADED) > 0)
	return QG_OBOL_DATABASE_SHADED;
    return QG_OBOL_DATABASE_WIREFRAME;
}

static const char *
qg_obol_draw_sync_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static int
qg_obol_draw_sync_path_equal(const char *a, const char *b)
{
    if (!a || !b)
	return 0;
    if (BU_STR_EQUAL(a, b))
	return 1;
    return BU_STR_EQUAL(qg_obol_draw_sync_skip_leading_slash(a),
	    qg_obol_draw_sync_skip_leading_slash(b));
}

static int
qg_obol_shape_record_matches_path(const struct ged_draw_shape_record *record,
	const char *path)
{
    if (!record || !path || !path[0])
	return 0;
    if (qg_obol_draw_sync_path_equal(record->display_name, path) ||
	    qg_obol_draw_sync_path_equal(record->leaf_name, path))
	return 1;

    if (!record->fullpath || record->fullpath->fp_len <= 0)
	return 0;

    char *recordPath = db_path_to_string(record->fullpath);
    if (!recordPath)
	return 0;
    int matched = qg_obol_draw_sync_path_equal(recordPath, path);
    bu_free(recordPath, "qtcad Obol draw sync record path");
    return matched;
}

static int
qg_obol_group_matches_view(struct ged *gedp,
	ged_draw_group_ref group,
	void *viewCtx)
{
    if (!gedp || !viewCtx || ged_draw_group_ref_is_null(group))
	return 0;

    struct ged_draw_group_record groupRecord;
    if (!ged_draw_group_record_get(gedp, group, &groupRecord))
	return 0;
    return ged_draw_group_record_in_view(&groupRecord, viewCtx) ? 1 : 0;
}

static void
qg_obol_source_state_add_group(qg_obol_ged_source_state &state,
	struct ged *gedp,
	ged_draw_group_ref group)
{
    if (!gedp || ged_draw_group_ref_is_null(group))
	return;

    struct ged_draw_group_record groupRecord;
    if (!ged_draw_group_record_get(gedp, group, &groupRecord) ||
	    !groupRecord.path || !groupRecord.path[0])
	return;

    state.groupValid = true;
    state.groupPath = groupRecord.path;
    state.groupDrawMode = qg_obol_lod_draw_mode_from_ged(
	    groupRecord.draw_mode);
    state.groupVisible = groupRecord.visible ? true : false;
    state.groupOverlay = groupRecord.is_overlay ? true : false;
    state.groupTransparency =
	static_cast<float>(groupRecord.transparency);
}

static void
qg_obol_source_state_from_record(qg_obol_ged_source_state &state,
	struct ged *gedp,
	const struct ged_draw_shape_record *record,
	int viewMatched)
{
    if (!record)
	return;

    state.valid = true;
    state.viewMatched = viewMatched ? true : false;
    state.sourceRevision = qg_obol_fold_revision(record->source_revision);
    state.inputsRevision = qg_obol_fold_revision(record->inputs_revision);
    state.drawMode = qg_obol_draw_mode_from_ged(record->draw_mode);
    state.visible = record->visible ? true : false;
    state.highlighted = record->highlighted ? true : false;
    state.lineWidth = record->line_width;
    state.transparency = static_cast<float>(record->transparency);
    qg_obol_source_state_add_group(state, gedp, record->group);

    void *shapeCtx = ged_draw_shape_ref_context(gedp, record->ref);
    struct ged_draw_scene_display_summary displaySummary;
    if (ged_draw_scene_context_display_summary(shapeCtx, &displaySummary) &&
	    displaySummary.valid) {
	state.lineStyle = displaySummary.line_style;
	state.lineWidth = displaySummary.line_width;
	state.transparency = static_cast<float>(displaySummary.transparency);
	state.materialColorValid = displaySummary.material_valid ? true : false;
	if (state.materialColorValid) {
	    state.materialColor = SbColor(
		    static_cast<float>(displaySummary.material_color[0]) /
			255.0f,
		    static_cast<float>(displaySummary.material_color[1]) /
			255.0f,
		    static_cast<float>(displaySummary.material_color[2]) /
			255.0f);
	}
    }

    struct ged_draw_shape_material_summary materialSummary;
    if (ged_draw_shape_ref_material_summary(gedp, record->ref,
	    &materialSummary) && materialSummary.valid)
	state.materialRevision =
	    qg_obol_fold_revision(materialSummary.material_revision);
}

struct qg_obol_find_source_state_context {
    struct ged *gedp;
    void *viewCtx;
    const char *path;
    qg_obol_ged_source_state state;
};

static int
qg_obol_find_source_state_cb(const struct ged_draw_shape_record *record,
	void *userdata)
{
    qg_obol_find_source_state_context *ctx =
	static_cast<qg_obol_find_source_state_context *>(userdata);
    if (!ctx || !record || !qg_obol_shape_record_matches_path(record,
	    ctx->path))
	return 1;

    const int viewMatched = qg_obol_group_matches_view(ctx->gedp,
	    record->group, ctx->viewCtx);
    if (!ctx->state.valid || viewMatched || !ctx->state.viewMatched)
	qg_obol_source_state_from_record(ctx->state, ctx->gedp, record,
		viewMatched);

    return viewMatched ? 0 : 1;
}

static qg_obol_ged_source_state
qg_obol_find_ged_source_state(struct ged *gedp,
	void *viewCtx,
	const char *path)
{
    if (!gedp || !path || !path[0])
	return qg_obol_ged_source_state();

    qg_obol_find_source_state_context ctx;
    ctx.gedp = gedp;
    ctx.viewCtx = viewCtx;
    ctx.path = path;
    ged_draw_foreach_shape_record(gedp, qg_obol_find_source_state_cb, &ctx);
    return ctx.state;
}

static int
qg_obol_intent_is_ged_draw_group(const SbString &intent)
{
    const char *value = intent.getString();
    if (!value)
	return 0;
    return strncmp(value, qg_obol_ged_group_intent_prefix,
	    sizeof(qg_obol_ged_group_intent_prefix) - 1) == 0;
}

static std::string
qg_obol_group_intent_path(const char *groupPath)
{
    std::string intent(qg_obol_ged_group_intent_prefix);
    if (groupPath)
	intent += qg_obol_draw_sync_skip_leading_slash(groupPath);
    return intent;
}

static int
qg_obol_sync_ged_group_state(BRLObolViewController *obol,
	const qg_obol_ged_source_state &state,
	const char *sourcePath)
{
    if (!obol || !state.groupValid || state.groupPath.empty())
	return 0;

    int changed = 0;
    SoGroup *existingGroup = obol->findGroup(state.groupPath.c_str());
    SoGroup *group = obol->ensureGroup(state.groupPath.c_str());
    if (!group)
	return 0;
    if (!existingGroup)
	changed = 1;

    const std::string intentPath =
	qg_obol_group_intent_path(state.groupPath.c_str());
    if (obol->setGroupDrawIntent(state.groupPath.c_str(),
	    intentPath.c_str(),
	    state.groupDrawMode,
	    BRLOBOL_LOD_DRAW_WIRE,
	    state.groupOverlay ? TRUE : FALSE,
	    state.sourceRevision) > 0)
	changed = 1;

    if (obol->setGroupDisplayState(state.groupPath.c_str(),
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

    if (sourcePath && sourcePath[0] &&
	    obol->moveDatabaseSourceToGroup(sourcePath,
		state.groupPath.c_str()) > 0)
	changed = 1;
    return changed;
}

static int
qg_obol_prune_empty_ged_groups(BRLObolViewController *obol)
{
    if (!obol || !obol->getSceneController())
	return 0;

    int changed = 0;
    int passRemoved = 1;
    while (passRemoved) {
	passRemoved = 0;
	std::vector<std::string> prunePaths;
	SoBRLSceneController *scene = obol->getSceneController();
	const int summaryCount = scene->getSceneTreeSummaryCount();
	for (int i = summaryCount - 1; i >= 0; i--) {
	    BRLObolSceneTreeSummary treeSummary;
	    BRLObolSceneDisplaySummary displaySummary;
	    if (!scene->getSceneTreeSummary(i, treeSummary) ||
		    !scene->getSceneDisplaySummary(i, displaySummary) ||
		    !treeSummary.valid ||
		    !displaySummary.valid ||
		    treeSummary.nodeKind !=
			BRLObolSceneTreeSummary::NODE_GROUP ||
		    treeSummary.childCount != 0 ||
		    treeSummary.path.getLength() == 0 ||
		    BU_STR_EQUAL(treeSummary.path.getString(), "/") ||
		    !displaySummary.hasDrawIntent ||
		    !qg_obol_intent_is_ged_draw_group(
			displaySummary.intentPath))
		continue;
	    prunePaths.push_back(treeSummary.path.getString());
	}
	for (const std::string &path : prunePaths) {
	    if (obol->removeGroup(path.c_str()) > 0) {
		passRemoved = 1;
		changed = 1;
	    }
	}
    }
    return changed;
}

static int
qg_obol_apply_ged_source_state(BRLObolViewController *obol,
	const char *sourcePath,
	const qg_obol_ged_source_state &state)
{
    if (!obol || !sourcePath || !sourcePath[0] || !state.valid)
	return 0;

    return obol->setDatabaseSourceState(sourcePath,
	    state.sourceRevision != 0 ? TRUE : FALSE,
	    state.sourceRevision,
	    state.inputsRevision,
	    state.visible ? TRUE : FALSE,
	    state.highlighted ? TRUE : FALSE,
	    state.lineStyle,
	    state.lineWidth,
	    state.transparency,
	    state.materialColorValid ? TRUE : FALSE,
	    state.materialColor,
	    state.materialRevision);
}

static int
qg_obol_replace_path(struct ged *gedp,
	void *viewCtx,
	struct db_i *dbip,
	const char *path,
	int drawMode,
	uint32_t sourceRevision,
	BRLObolViewController *obol)
{
    if (!dbip || !path || !path[0] || !obol)
	return 0;

    qg_obol_ged_source_state sourceState =
	qg_obol_find_ged_source_state(gedp, viewCtx, path);
    if (sourceState.valid) {
	drawMode = sourceState.drawMode;
	if (sourceState.sourceRevision != 0)
	    sourceRevision = sourceState.sourceRevision;
    }

    int changed = obol->replaceDatabaseSource(path, dbip,
	    qg_obol_database_source_draw_mode(drawMode), sourceRevision);
    if (qg_obol_apply_ged_source_state(obol, path, sourceState))
	changed = 1;
    if (qg_obol_sync_ged_group_state(obol, sourceState, path))
	changed = 1;
    return changed;
}

static int
qg_obol_replace_paths(struct db_i *dbip,
	const std::vector<std::string> &paths,
	int drawMode,
	uint32_t sourceRevision,
	struct ged *gedp,
	void *viewCtx,
	QgView *display)
{
    if (!dbip || paths.empty() || !display)
	return 0;

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    int changed = 0;
    for (const std::string &path : paths) {
	if (qg_obol_replace_path(gedp, viewCtx, dbip, path.c_str(),
		drawMode, sourceRevision, obol) > 0)
	    changed = 1;
    }

    if (changed) {
	obol->realizePending();
	display->need_update(QG_VIEW_REFRESH);
    }
    return changed;
}

static int
qg_obol_remove_paths(const std::vector<std::string> &paths,
	QgView *display)
{
    if (paths.empty() || !display)
	return 0;

    std::vector<const char *> cpaths;
    cpaths.reserve(paths.size());
    for (const std::string &path : paths)
	cpaths.push_back(path.c_str());

    int changed = qg_obol_remove_database_sources(cpaths.data(),
	static_cast<int>(cpaths.size()), display);
    if (qg_obol_prune_empty_ged_groups(display->obolViewController())) {
	changed = 1;
	display->need_update(QG_VIEW_REFRESH);
    }
    return changed;
}

static int
qg_obol_sync_full_scene(struct ged *gedp,
	void *view_ctx,
	uint32_t sourceRevision,
	QgView *display)
{
    if (!display)
	return 0;

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    int changed = qg_obol_clear_database_sources(display);
    if (qg_obol_prune_empty_ged_groups(obol))
	changed = 1;
    struct db_i *dbip = qg_legacy_view_ged_database(gedp);
    if (!dbip)
	return changed;

    struct bu_vls listed = BU_VLS_INIT_ZERO;
    (void)ged_draw_list_paths(gedp, view_ctx, -1, 0, &listed);

    std::vector<std::string> paths;
    qg_obol_append_unique_paths_from_words(paths, bu_vls_cstr(&listed));
    bu_vls_free(&listed);

    for (const std::string &path : paths) {
	int drawMode = qg_obol_drawn_path_mode(gedp, view_ctx, path.c_str());
	if (qg_obol_replace_path(gedp, view_ctx, dbip, path.c_str(),
		drawMode, sourceRevision, obol) > 0)
	    changed = 1;
    }

    if (changed) {
	obol->realizePending();
	display->need_update(QG_VIEW_REFRESH);
    }
    return changed;
}

int
qg_obol_display_accepts_ged_draw_transaction_view(
	const struct ged_draw_transaction *txn,
	QgView *display)
{
    if (!display)
	return 0;
    if (!txn || !txn->view)
	return 1;
    return qg_legacy_view_to_context(display->view()) == txn->view;
}

int
qg_obol_ged_draw_transaction_has_view(
	const struct ged_draw_transaction *txn)
{
    return (txn && txn->view) ? 1 : 0;
}

int
qg_obol_sync_ged_draw_transaction(struct ged *gedp,
	const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result,
	QgView *display)
{
    if (!txn || (result && result->status < 0) || !display)
	return 0;

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    void *view_ctx = txn->view;
    if (!view_ctx)
	view_ctx = qg_legacy_view_to_context(display->view());
    uint32_t sourceRevision = qg_obol_ged_source_revision(result);
    int changed = 0;

    switch (txn->kind) {
	case GED_DRAW_TXN_DRAW:
	{
	    std::vector<std::string> paths =
		qg_obol_ged_transaction_paths(txn, result);
	    if (paths.empty()) {
		changed = qg_obol_sync_full_scene(gedp, view_ctx,
			sourceRevision, display);
	    } else {
		changed = qg_obol_replace_paths(
			qg_legacy_view_ged_database(gedp), paths,
			qg_obol_ged_transaction_draw_mode(gedp, txn),
			sourceRevision, gedp, view_ctx, display);
	    }
	    break;
	}
	case GED_DRAW_TXN_ERASE:
	{
	    std::vector<std::string> paths =
		qg_obol_ged_transaction_paths(txn, result);
	    if (paths.empty())
		changed = qg_obol_sync_full_scene(gedp, view_ctx,
			sourceRevision, display);
	    else
		changed = qg_obol_remove_paths(paths, display);
	    break;
	}
	case GED_DRAW_TXN_CLEAR:
	case GED_DRAW_TXN_CLEAR_SCOPE:
	case GED_DRAW_TXN_TEARDOWN:
	    changed = qg_obol_clear_database_sources(display);
	    if (qg_obol_prune_empty_ged_groups(obol)) {
		changed = 1;
		display->need_update(QG_VIEW_REFRESH);
	    }
	    break;
	case GED_DRAW_TXN_VISIBILITY:
	case GED_DRAW_TXN_MATERIAL_CHANGED:
	case GED_DRAW_TXN_REFRESH_MATERIAL_COLORS:
	case GED_DRAW_TXN_STALE_SOURCE:
	case GED_DRAW_TXN_SOURCE_UPDATED:
	case GED_DRAW_TXN_SOURCE_RENAMED:
	case GED_DRAW_TXN_SOURCE_REFERENCES_REMOVED:
	case GED_DRAW_TXN_ERASE_PREFIX:
	case GED_DRAW_TXN_REDRAW:
	    changed = qg_obol_sync_full_scene(gedp, view_ctx,
		    sourceRevision, display);
	    break;
	case GED_DRAW_TXN_DEFAULT_DRAW_MODE:
	case GED_DRAW_TXN_TRANSPARENCY:
	case GED_DRAW_TXN_HIGHLIGHT:
	case GED_DRAW_TXN_NONE:
	default:
	    break;
    }

    return changed;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
