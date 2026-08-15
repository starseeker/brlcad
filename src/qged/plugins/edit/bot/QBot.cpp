/*                       Q B O T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QBot.cpp
 *
 * BoT editing adapter.  GED owns the sole mutable primitive; this class keeps
 * only a revision-tagged retained presentation cache.  Interactive vertex,
 * edge, and face moves patch the affected retained vertices in place and
 * never copy or rebuild the full mesh on a pointer-motion event.
 */

#include "common.h"

#include <QComboBox>
#include <QEvent>
#include <QGroupBox>
#include <QLabel>
#include <QMouseEvent>
#include <QSignalBlocker>
#include <QVBoxLayout>

#include "BObol/BMeshShape.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "bn/mat.h"
#include "bu/malloc.h"
#include "bv.h"
#include "ged.h"
#include "ged/edit.h"
#include "ged/plugin/obol.h"
#include "ged/selection.h"
#include "ged/view_feature.h"
#include "ged/view_feature_batch.h"
#include "rt/db_fullpath.h"
#include "rt/edit.h"
#include "rt/geom.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgPrimitiveEdit.h"
#include "qtcad/QgSignalFlags.h"
#include "../qged_edit_preview_util.h"
#include "QBot.h"

#include <algorithm>
#include <cmath>
#include <vector>


namespace {

static const char *qbot_surface_name = "_bot_edit_surface";
static const char *qbot_handle_name = "_bot_edit_handle";

static bool
same_session(ged_edit_session_ref a, ged_edit_session_ref b)
{
    return a.owner == b.owner && a.id == b.id &&
	a.generation == b.generation;
}

static fastf_t
point_segment_distance_sq(const point_t point, const SbVec3f &a,
	const SbVec3f &b)
{
    const SbVec3f p(static_cast<float>(point[X]),
	static_cast<float>(point[Y]), static_cast<float>(point[Z]));
    const SbVec3f ab = b - a;
    const float lengthSq = ab.sqrLength();
    float t = lengthSq > static_cast<float>(SMALL_FASTF) ?
	(p - a).dot(ab) / lengthSq : 0.0f;
    t = std::max(0.0f, std::min(1.0f, t));
    return static_cast<fastf_t>((p - (a + ab * t)).sqrLength());
}

}


struct QBotPresentationState {
    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    uint64_t revision = UINT64_MAX;
    std::vector<int32_t> faces;
    std::vector<int> selected_vertices;
    int selected_face = -1;
    int pick_commands[3] = {0, 0, 0};
    int move_commands[3] = {0, 0, 0};
    mat_t path_matrix = MAT_INIT_IDN;
    mat_t inverse_path_matrix = MAT_INIT_IDN;
    bool matrix_valid = false;
    SoBRLMeshShape *geometry = nullptr;
    std::vector<struct ged_view_context *> view_contexts;
    bool dragging = false;
    struct ged_view_context *drag_view_context = nullptr;
    point_t drag_start = VINIT_ZERO;
    point_t drag_anchor = VINIT_ZERO;
};


QBot::QBot()
    : QWidget(), state(new QBotPresentationState)
{
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);

    QGroupBox *modeBox = new QGroupBox(tr("Viewport selection"), this);
    QVBoxLayout *modeLayout = new QVBoxLayout(modeBox);
    edit_mode = new QComboBox(modeBox);
    edit_mode->setObjectName(QStringLiteral("bot.editMode"));
    edit_mode->setProperty("qgTestId", QStringLiteral("bot.editMode"));
    edit_mode->addItem(tr("Vertex"));
    edit_mode->addItem(tr("Face"));
    edit_mode->addItem(tr("Edge"));
    modeLayout->addWidget(edit_mode);
    layout->addWidget(modeBox);

    editor = new QgPrimitiveEdit(this);
    editor->setObjectName(QStringLiteral("bot.sharedPrimitiveEditor"));
    layout->addWidget(editor);

    connect(editor, &QgPrimitiveEdit::targetChanged,
	this, &QBot::target_changed);
    connect(editor, &QgPrimitiveEdit::sessionEvent,
	this, &QBot::update_preview);
}


QBot::~QBot()
{
    clear_preview();
    delete state;
    state = nullptr;
}


void
QBot::setContext(QgPluginContext *ctx)
{
    m_ctx = ctx;
    editor->setGed(getGed());
}


struct ged *
QBot::getGed() const
{
    return m_ctx ? m_ctx->getGed() : nullptr;
}


int
QBot::command_id(const char *name) const
{
    if (!name || !editor || ged_edit_session_ref_is_null(editor->session()))
	return 0;
    const struct rt_edit_prim_desc *descriptor = nullptr;
    if (ged_edit_session_descriptor_get(getGed(), editor->session(),
	    &descriptor) != GED_EDIT_OK || !descriptor)
	return 0;
    for (int i = 0; i < descriptor->ncmd; i++) {
	if (descriptor->cmds[i].name &&
	    BU_STR_EQUAL(descriptor->cmds[i].name, name))
	    return descriptor->cmds[i].cmd_id;
    }
    return 0;
}


void
QBot::clear_view_presentations()
{
    if (!state)
	return;

    for (struct ged_view_context *viewContext : state->view_contexts) {
	if (!viewContext)
	    continue;
	(void)ged_view_feature_remove(viewContext, qbot_surface_name);
	(void)ged_view_feature_remove(viewContext, qbot_handle_name);
    }
    state->view_contexts.clear();
}


bool
QBot::sync_view_presentations(bool replaceAll)
{
    if (!state || !state->geometry)
	return false;

    const std::vector<struct ged_view_context *> contexts =
	qged_edit_ged_view_contexts(m_ctx);
    for (auto it = state->view_contexts.begin();
	it != state->view_contexts.end();) {
	struct ged_view_context *viewContext = *it;
	if (replaceAll || std::find(contexts.begin(), contexts.end(),
		viewContext) == contexts.end() ||
	    !ged_plugin_obol_view_controller(viewContext)) {
	    if (viewContext) {
		(void)ged_view_feature_remove(viewContext, qbot_surface_name);
		(void)ged_view_feature_remove(viewContext, qbot_handle_name);
	    }
	    it = state->view_contexts.erase(it);
	    continue;
	}
	++it;
    }

    size_t eligible = 0;
    for (struct ged_view_context *viewContext : contexts) {
	BObolViewController *controller =
	    ged_plugin_obol_view_controller(viewContext);
	if (!controller)
	    continue;
	eligible++;
	if (std::find(state->view_contexts.begin(),
		state->view_contexts.end(), viewContext) !=
		state->view_contexts.end())
	    continue;

	SoBRLMeshShape *node = new SoBRLMeshShape;
	node->ref();
	node->sourcePath = qbot_surface_name;
	node->sourceName = qbot_surface_name;
	node->sourceType = "indexed-face-set";
	node->displayName = qbot_surface_name;
	node->geometryName = qbot_surface_name;
	node->sourceIdentity = qbot_surface_name;
	node->cacheIdentity = qbot_surface_name;
	node->databaseIntent = FALSE;
	node->overlayIntent = TRUE;
	node->localSource = TRUE;
	node->sharedSource = FALSE;
	node->nonDatabaseSource = TRUE;
	node->drawMode = BOBOL_LOD_DRAW_SHADED;
	node->recordRole = "edit-surface";
	node->geometryKind = "surface";
	node->colorOverride = TRUE;
	node->color = SbColor(72.0f / 255.0f, 126.0f / 255.0f,
	    168.0f / 255.0f);
	node->visible = TRUE;
	node->selectable = TRUE;
	node->setSharedGeometry(state->geometry);

	BObolFeatureStyle style;
	style.hasVisible = TRUE;
	style.visible = TRUE;
	style.hasSelectable = TRUE;
	style.selectable = TRUE;
	style.hasColor = TRUE;
	style.color = SbColor(72.0f / 255.0f, 126.0f / 255.0f,
	    168.0f / 255.0f);
	BObolFeatureOwner owner;
	owner.ownerToken = this;
	owner.ownerId = "qged::bot-edit";
	owner.ownerRole = "edit-surface";
	const BObolFeatureHandle feature =
	    controller->features().publishCustomNode(qbot_surface_name,
		BObolFeatureScope::Local, node, &style, &owner);
	node->unref();
	if (!feature.isValid())
	    continue;

	BObolOverlayInfo overlay;
	overlay.isOverlay = TRUE;
	overlay.ownerToken = this;
	overlay.role = BObolOverlayRole::XRay;
	overlay.overlayClass = BObolOverlayClass::EditHandle;
	overlay.lifecycle = BObolOverlayLifecycle::PerTool;
	overlay.order = BObolOverlayOrder::XRay;
	overlay.sortOrder = 90;
	overlay.sourcePath = preview_path.toUtf8().constData();
	(void)controller->features().setOverlayInfo(feature, overlay);
	state->view_contexts.push_back(viewContext);
    }
    return eligible > 0 && state->view_contexts.size() == eligible;
}


void
QBot::clear_preview()
{
    if (!state)
	return;

    clear_view_presentations();
    if (state->geometry) {
	state->geometry->unref();
	state->geometry = nullptr;
    }
    state->session = GED_EDIT_SESSION_REF_NULL;
    state->revision = UINT64_MAX;
    state->faces.clear();
    state->selected_vertices.clear();
    state->selected_face = -1;
    state->matrix_valid = false;
    state->dragging = false;
    state->drag_view_context = nullptr;
    preview_path.clear();
}


void
QBot::target_changed(const QString &path)
{
    if (path != preview_path)
	clear_preview();
}


bool
QBot::rebuild_preview(uint64_t revision)
{
    struct ged *gedp = getGed();
    const ged_edit_session_ref session = editor->session();
    const QString path = editor->targetPath();
    if (!state || !gedp || !gedp->dbip || path.isEmpty() ||
	ged_edit_session_ref_is_null(session))
	return false;

    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    if (ged_edit_session_internal_copy(gedp, session, &intern) != GED_EDIT_OK ||
	intern.idb_type != ID_BOT || !intern.idb_ptr) {
	rt_db_free_internal(&intern);
	return false;
    }
    struct rt_bot_internal *bot =
	static_cast<struct rt_bot_internal *>(intern.idb_ptr);
    RT_BOT_CK_MAGIC(bot);

    struct db_full_path fullPath;
    db_full_path_init(&fullPath);
    mat_t pathMatrix;
    MAT_IDN(pathMatrix);
    const QByteArray pathBytes = path.toUtf8();
    const bool pathOk = db_string_to_path(&fullPath, gedp->dbip,
	pathBytes.constData()) == 0 && db_path_to_mat(gedp->dbip, &fullPath,
	pathMatrix, static_cast<int>(fullPath.fp_len) - 1);
    db_free_full_path(&fullPath);
    if (!pathOk) {
	rt_db_free_internal(&intern);
	return false;
    }
    std::vector<SbVec3f> displayPoints(bot->num_vertices);
    for (size_t i = 0; i < bot->num_vertices; i++) {
	point_t displayPoint;
	MAT4X3PNT(displayPoint, pathMatrix, &bot->vertices[i * 3]);
	displayPoints[i] = SbVec3f(
	    static_cast<float>(displayPoint[X]),
	    static_cast<float>(displayPoint[Y]),
	    static_cast<float>(displayPoint[Z]));
    }
    state->faces.resize(bot->num_faces * 3);
    for (size_t i = 0; i < state->faces.size(); i++)
	state->faces[i] = static_cast<int32_t>(bot->faces[i]);

    clear_view_presentations();
    if (state->geometry)
	state->geometry->unref();
    state->geometry = new SoBRLMeshShape;
    state->geometry->ref();
    if (!displayPoints.empty() && !state->faces.empty()) {
	state->geometry->setIndexedTriangles(displayPoints.data(),
	    static_cast<int>(displayPoints.size()), state->faces.data(),
	    static_cast<int>(state->faces.size()));
    }

    state->session = session;
    state->revision = revision;
    MAT_COPY(state->path_matrix, pathMatrix);
    bn_mat_inv(state->inverse_path_matrix, pathMatrix);
    state->matrix_valid = true;
    state->selected_vertices.clear();
    state->selected_face = -1;
    state->pick_commands[0] = command_id("ECMD_BOT_PICKV");
    state->pick_commands[1] = command_id("ECMD_BOT_PICKT");
    state->pick_commands[2] = command_id("ECMD_BOT_PICKE");
    state->move_commands[0] = command_id("ECMD_BOT_MOVEV");
    state->move_commands[1] = command_id("ECMD_BOT_MOVET");
    state->move_commands[2] = command_id("ECMD_BOT_MOVEE");
    preview_path = path;

    if (!sync_view_presentations(true)) {
	rt_db_free_internal(&intern);
	return false;
    }
    struct ged_edit_session_info info = {};
    if (ged_edit_session_info_get(gedp, session, &info) == GED_EDIT_OK &&
	info.command_id)
	(void)sync_session_selection(info.command_id);
    else
	publish_selection();

    rt_db_free_internal(&intern);
    emit view_updated(QG_VIEW_REFRESH);
    return true;
}


bool
QBot::sync_session_selection(int commandId)
{
    if (!state || ged_edit_session_ref_is_null(editor->session()))
	return false;

    int mode = -1;
    for (int i = 0; i < 3; i++) {
	if (commandId == state->pick_commands[i] ||
	    commandId == state->move_commands[i]) {
	    mode = i;
	    break;
	}
    }
    if (mode < 0 || !state->pick_commands[mode])
	return false;

    struct rt_edit_cmd_values values;
    if (ged_edit_session_command_values_get(getGed(), editor->session(),
	    state->pick_commands[mode], &values) != GED_EDIT_OK)
	return false;
    const size_t required = mode == 0 ? 1u : mode == 2 ? 2u : 3u;
    if (values.value_count < required)
	return false;

    std::vector<int> selected;
    for (size_t i = 0; i < required; i++) {
	if (!values.value_valid[i])
	    return false;
	const int vertex = static_cast<int>(values.values[i]);
	if (vertex < 0 || !state->geometry ||
	    vertex >= state->geometry->point.getNum())
	    return false;
	selected.push_back(vertex);
    }
    state->selected_vertices = selected;
    state->selected_face = -1;
    const size_t faceCount = state->faces.size() / 3;
    for (size_t fi = 0; fi < faceCount; fi++) {
	const int32_t *face = &state->faces[fi * 3];
	bool match = false;
	if (mode == 0) {
	    match = face[0] == selected[0] || face[1] == selected[0] ||
		face[2] == selected[0];
	} else if (mode == 1) {
	    match = face[0] == selected[0] && face[1] == selected[1] &&
		face[2] == selected[2];
	} else {
	    for (int edge = 0; edge < 3 && !match; edge++) {
		const int a = face[edge];
		const int b = face[(edge + 1) % 3];
		match = (a == selected[0] && b == selected[1]) ||
		    (a == selected[1] && b == selected[0]);
	    }
	}
	if (match) {
	    state->selected_face = static_cast<int>(fi);
	    break;
	}
    }

    {
	QSignalBlocker blocker(edit_mode);
	edit_mode->setCurrentIndex(mode);
    }
    publish_selection();
    return true;
}


bool
QBot::patch_session_move(int commandId)
{
    if (!state || !state->geometry || state->view_contexts.empty() ||
	!state->matrix_valid)
	return false;
    int mode = -1;
    for (int i = 0; i < 3; i++) {
	if (commandId == state->move_commands[i]) {
	    mode = i;
	    break;
	}
    }
    if (mode < 0)
	return false;
    const size_t required = mode == 0 ? 1u : mode == 2 ? 2u : 3u;
    if (state->selected_vertices.size() != required &&
	!sync_session_selection(commandId))
	return false;

    struct rt_edit_cmd_values values;
    if (ged_edit_session_command_values_get(getGed(), editor->session(),
	    commandId, &values) != GED_EDIT_OK || values.value_count < 3 ||
	!values.value_valid[0] || !values.value_valid[1] ||
	!values.value_valid[2])
	return false;

    const fastf_t local2base = getGed()->dbip->dbi_local2base;
    point_t sourceTarget = {
	values.values[0] * local2base,
	values.values[1] * local2base,
	values.values[2] * local2base
    };
    point_t displayTarget;
    MAT4X3PNT(displayTarget, state->path_matrix, sourceTarget);

    point_t currentAnchor = VINIT_ZERO;
    for (const int vertex : state->selected_vertices) {
	const SbVec3f &point = state->geometry->point[vertex];
	currentAnchor[X] += point[0];
	currentAnchor[Y] += point[1];
	currentAnchor[Z] += point[2];
    }
    VSCALE(currentAnchor, currentAnchor,
	1.0 / static_cast<fastf_t>(state->selected_vertices.size()));
    vect_t delta;
    VSUB2(delta, displayTarget, currentAnchor);

    int patchIndices[3] = {-1, -1, -1};
    SbVec3f updatedPoints[3];
    for (size_t i = 0; i < state->selected_vertices.size(); i++) {
	const int vertex = state->selected_vertices[i];
	updatedPoints[i] = state->geometry->point[vertex] +
	    SbVec3f(static_cast<float>(delta[X]),
	    static_cast<float>(delta[Y]), static_cast<float>(delta[Z]));
	patchIndices[i] = vertex;
    }
    const SbBool notify = state->geometry->enableNotify(FALSE);
    for (size_t i = 0; i < state->selected_vertices.size(); i++) {
	state->geometry->point.set1Value(patchIndices[i], updatedPoints[i]);
    }
    state->geometry->normal.setNum(0);
    state->geometry->sourceId = state->geometry->sourceId.getValue() + 1;
    state->geometry->enableNotify(notify);
    if (notify)
	state->geometry->touch();
    for (struct ged_view_context *viewContext : state->view_contexts) {
	BObolViewController *controller =
	    ged_plugin_obol_view_controller(viewContext);
	if (!controller)
	    continue;
	const BObolFeatureHandle feature = controller->features().find(
	    qbot_surface_name, BOBOL_FEATURE_SCOPE_LOCAL);
	SoNode *rawNode = feature.isValid() ?
	    controller->features().node(feature) : nullptr;
	if (rawNode && rawNode->isOfType(SoBRLMeshShape::getClassTypeId())) {
	    SoBRLMeshShape *node = static_cast<SoBRLMeshShape *>(rawNode);
	    node->clearRenderLists();
	    node->touch();
	}
	controller->requestPresentationRender("bot-edit-points");
    }
    publish_selection();
    return true;
}


void
QBot::publish_selection()
{
    if (!state || state->view_contexts.empty())
	return;

    point_t points[3] = {VINIT_ZERO, VINIT_ZERO, VINIT_ZERO};
    const size_t count = std::min(state->selected_vertices.size(),
	static_cast<size_t>(3));
    for (size_t i = 0; i < count; i++) {
	const int vertex = state->selected_vertices[i];
	if (vertex < 0 || !state->geometry ||
	    vertex >= state->geometry->point.getNum())
	    return;
	const SbVec3f &point = state->geometry->point[vertex];
	VSET(points[i], point[0], point[1], point[2]);
    }

    struct ged_view_feature_style style = ged_view_feature_style_default();
    style.visible = 1;
    style.selectable = 1;
    style.color_valid = 1;
    style.color[0] = 255;
    style.color[1] = 196;
    style.color[2] = 32;
    struct ged_view_feature_batch_desc desc =
	ged_view_feature_batch_desc_default();
    desc.owner_id = "qged-bot-edit";
    desc.owner_role = "edit-handle";
    desc.overlay_class = GED_VIEW_FEATURE_OVERLAY_CLASS_EDIT_HANDLE;
    desc.lifecycle = GED_VIEW_FEATURE_LIFECYCLE_PER_TOOL;
    desc.local = 1;
    for (struct ged_view_context *viewContext : state->view_contexts) {
	if (state->selected_face >= 0)
	    (void)ged_view_feature_set_selection(viewContext,
		qbot_surface_name, &state->selected_face, 1);
	else
	    (void)ged_view_feature_set_selection(viewContext,
		qbot_surface_name, nullptr, 0);

	if (state->selected_vertices.empty()) {
	    (void)ged_view_feature_remove(viewContext, qbot_handle_name);
	    continue;
	}
	struct ged_view_feature_batch *batch =
	    ged_view_feature_batch_begin(viewContext, &desc);
	if (batch && ged_view_feature_batch_point_set_replace(batch,
		qbot_handle_name, points, count, &style))
	    (void)ged_view_feature_batch_commit(batch);
	else if (batch)
	    ged_view_feature_batch_abort(batch);
    }
}


void
QBot::update_preview(int kindValue, qulonglong revisionValue)
{
    const enum ged_edit_session_event_kind kind =
	static_cast<enum ged_edit_session_event_kind>(kindValue);
    if (kind == GED_EDIT_SESSION_COMMIT || kind == GED_EDIT_SESSION_CANCEL ||
	kind == GED_EDIT_SESSION_INVALIDATE) {
	clear_preview();
	emit view_updated(kind == GED_EDIT_SESSION_COMMIT ?
	    QG_VIEW_DB : QG_VIEW_REFRESH);
	return;
    }

    struct ged *gedp = getGed();
    const ged_edit_session_ref session = editor->session();
    if (!state || !gedp || ged_edit_session_ref_is_null(session))
	return;

    struct ged_edit_session_info info = {};
    if (ged_edit_session_info_get(gedp, session, &info) != GED_EDIT_OK)
	return;
    const uint64_t revision = static_cast<uint64_t>(revisionValue);
    const bool replacement = !same_session(state->session, session);
    if (!replacement && state->revision == info.revision &&
	kind != GED_EDIT_SESSION_REVERT) {
	if (sync_view_presentations(false))
	    publish_selection();
	return;
    }

    if (replacement || kind == GED_EDIT_SESSION_BEGIN ||
	kind == GED_EDIT_SESSION_REVERT || !state->geometry ||
	state->geometry->point.getNum() == 0) {
	(void)rebuild_preview(info.revision);
	return;
    }

    bool handled = false;
    for (int i = 0; i < 3; i++) {
	if (info.command_id == state->pick_commands[i]) {
	    handled = sync_session_selection(info.command_id);
	    break;
	}
	if (info.command_id == state->move_commands[i]) {
	    handled = patch_session_move(info.command_id);
	    break;
	}
    }
    if (!handled) {
	const struct rt_edit_prim_desc *descriptor = nullptr;
	const struct rt_edit_cmd_desc *command = nullptr;
	if (ged_edit_session_descriptor_get(gedp, session, &descriptor) ==
		GED_EDIT_OK && descriptor) {
	    for (int i = 0; i < descriptor->ncmd; i++) {
		if (descriptor->cmds[i].cmd_id == info.command_id) {
		    command = &descriptor->cmds[i];
		    break;
		}
	    }
	}
	/* Property changes do not alter vertex/topology presentation.  Unknown
	 * and topology commands require one full refresh, never one per motion. */
	if (!command || !command->category ||
	    !BU_STR_EQUAL(command->category, "properties")) {
	    (void)rebuild_preview(info.revision);
	    return;
	}
    }

    state->revision = revision ? revision : info.revision;
    emit view_updated(QG_VIEW_REFRESH);
}


void
QBot::refresh_preview()
{
    editor->refreshFromSession();
    update_preview(static_cast<int>(GED_EDIT_SESSION_UPDATE), 0);
}


void
QBot::sync_selection()
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip)
	return;

    QString selectedPath;
    if (ged_selection_count(gedp, nullptr) == 1) {
	struct bu_vls paths = BU_VLS_INIT_ZERO;
	(void)ged_selection_list_paths(gedp, nullptr, &paths);
	selectedPath = QString::fromUtf8(bu_vls_cstr(&paths)).trimmed();
	bu_vls_free(&paths);

	struct db_full_path fullPath;
	db_full_path_init(&fullPath);
	const QByteArray pathBytes = selectedPath.toUtf8();
	const bool isBot = !selectedPath.isEmpty() &&
	    db_string_to_path(&fullPath, gedp->dbip,
		pathBytes.constData()) == 0 &&
	    DB_FULL_PATH_CUR_DIR(&fullPath) &&
	    DB_FULL_PATH_CUR_DIR(&fullPath)->d_minor_type ==
		DB5_MINORTYPE_BRLCAD_BOT;
	db_free_full_path(&fullPath);
	if (!isBot)
	    selectedPath.clear();
    }

    if (selectedPath.isEmpty()) {
	if (!selection_path.isEmpty() &&
	    editor->targetPath() == selection_path)
	    editor->setTargetPath(QString());
	selection_path.clear();
	return;
    }
    selection_path = selectedPath;
    editor->setTargetPath(selectedPath);
}


void
QBot::reset_for_database()
{
    selection_path.clear();
    clear_preview();
    editor->setGed(nullptr);
    editor->setTargetPath(QString());
    editor->setGed(getGed());
    emit view_updated(QG_VIEW_REFRESH);
}


bool
QBot::eventFilter(QObject *, QEvent *event)
{
    if (!event || !state || !m_ctx || !state->geometry ||
	state->geometry->point.getNum() == 0)
	return false;
    const QEvent::Type type = event->type();
    if (type != QEvent::MouseMove && type != QEvent::MouseButtonPress &&
	type != QEvent::MouseButtonRelease)
	return false;
    QMouseEvent *mouse = static_cast<QMouseEvent *>(event);
    struct ged_view_context *viewContext = qged_edit_ged_view_context(m_ctx);
    if (!viewContext || std::find(state->view_contexts.begin(),
	    state->view_contexts.end(), viewContext) ==
	    state->view_contexts.end())
	return false;

    if (type == QEvent::MouseMove && state->dragging) {
	if (state->drag_view_context != viewContext)
	    return false;
	point_t current;
	const struct bv *view = bv_context_view_const(
	    ged_view_context_bv_const(viewContext));
	if (!view || !bv_screen_to_model(current, view,
		mouse->position().x(), mouse->position().y()))
	    return true;
	vect_t delta;
	VSUB2(delta, current, state->drag_start);
	point_t displayTarget;
	VADD2(displayTarget, state->drag_anchor, delta);
	point_t sourceTarget;
	MAT4X3PNT(sourceTarget, state->inverse_path_matrix, displayTarget);
	const fastf_t base2local = getGed()->dbip->dbi_base2local;
	fastf_t values[3] = {
	    sourceTarget[X] * base2local,
	    sourceTarget[Y] * base2local,
	    sourceTarget[Z] * base2local
	};
	const int mode = edit_mode->currentIndex();
	if (mode < 0 || mode >= 3 || !state->move_commands[mode])
	    return true;
	struct ged_edit_command_input input = {};
	input.command_id = state->move_commands[mode];
	input.values = values;
	input.value_count = 3;
	input.view = viewContext;
	const enum ged_edit_status result = ged_edit_session_apply(getGed(),
	    editor->session(), &input);
	if (result != GED_EDIT_OK) {
	    state->dragging = false;
	    state->drag_view_context = nullptr;
	    if (result == GED_EDIT_REJECTED || result == GED_EDIT_ERROR)
		(void)ged_edit_session_revert(getGed(), editor->session());
	    else
		editor->refreshFromSession();
	}
	return true;
    }

    if (type == QEvent::MouseButtonRelease && state->dragging &&
	mouse->button() == Qt::LeftButton) {
	state->dragging = false;
	state->drag_view_context = nullptr;
	emit view_updated(QG_VIEW_REFRESH);
	return true;
    }

    if (type != QEvent::MouseButtonPress ||
	mouse->button() != Qt::LeftButton)
	return false;

    /* The editable surface deliberately overlaps the promoted database
     * occurrence.  Ask for all ordered hits and select our retained feature;
     * using only the globally nearest hit made mesh editing dependent on
     * coincident-node traversal order. */
    struct ged_pick_result *result = ged_pick_point(viewContext,
	static_cast<int>(mouse->position().x()),
	static_cast<int>(mouse->position().y()), 0);
    struct ged_pick_detail detail = ged_pick_detail_default();
    struct bu_vls path = BU_VLS_INIT_ZERO;
    bool picked = false;
    const size_t resultCount = ged_pick_result_count(result);
    for (size_t i = 0; i < resultCount; i++) {
	bu_vls_trunc(&path, 0);
	struct ged_pick_detail candidate = ged_pick_detail_default();
	if (!ged_pick_result_path(result, i, &path) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&path), qbot_surface_name) ||
	    !ged_pick_result_detail(result, i, &candidate) ||
	    candidate.primitive_kind != 3 ||
	    candidate.primitive_index < 0)
	    continue;
	detail = candidate;
	picked = true;
	break;
    }
    if (!picked) {
	bu_vls_free(&path);
	ged_pick_result_free(result);
	return false;
    }

    const int faceIndex = detail.primitive_index;
    if (static_cast<size_t>(faceIndex) >= state->faces.size() / 3) {
	bu_vls_free(&path);
	ged_pick_result_free(result);
	return false;
    }
    const int32_t *face =
	&state->faces[static_cast<size_t>(faceIndex) * 3];
    const int mode = edit_mode->currentIndex();
    fastf_t values[3] = {0.0, 0.0, 0.0};
    size_t valueCount = 0;
    if (mode == 0) {
	int vertex = detail.nearest_face_vertex_index;
	if (vertex < 0)
	    vertex = face[0];
	values[0] = vertex;
	valueCount = 1;
    } else if (mode == 1) {
	values[0] = face[0];
	values[1] = face[1];
	values[2] = face[2];
	valueCount = 3;
    } else if (mode == 2) {
	int edge = 0;
	if (detail.model_point_valid) {
	    fastf_t best = INFINITY;
	    for (int i = 0; i < 3; i++) {
		const fastf_t distance = point_segment_distance_sq(
		    detail.model_point,
		    state->geometry->point[face[i]],
		    state->geometry->point[face[(i + 1) % 3]]);
		if (distance < best) {
		    best = distance;
		    edge = i;
		}
	    }
	}
	values[0] = face[edge];
	values[1] = face[(edge + 1) % 3];
	valueCount = 2;
    } else {
	bu_vls_free(&path);
	ged_pick_result_free(result);
	return false;
    }

    struct ged_edit_command_input input = {};
    input.command_id = state->pick_commands[mode];
    input.values = values;
    input.value_count = valueCount;
    input.view = viewContext;
    if (!input.command_id || ged_edit_session_apply(getGed(),
	editor->session(), &input) != GED_EDIT_OK ||
	ged_edit_session_checkpoint(getGed(), editor->session()) != GED_EDIT_OK) {
	bu_vls_free(&path);
	ged_pick_result_free(result);
	return true;
    }

    const struct bv *view = bv_context_view_const(
	ged_view_context_bv_const(viewContext));
    if (!view || !bv_screen_to_model(state->drag_start, view,
	    mouse->position().x(), mouse->position().y()))
	VSETALL(state->drag_start, 0.0);
    VSETALL(state->drag_anchor, 0.0);
    for (const int vertex : state->selected_vertices) {
	const SbVec3f &point = state->geometry->point[vertex];
	state->drag_anchor[X] += point[0];
	state->drag_anchor[Y] += point[1];
	state->drag_anchor[Z] += point[2];
    }
    VSCALE(state->drag_anchor, state->drag_anchor,
	1.0 / static_cast<fastf_t>(state->selected_vertices.size()));
    state->dragging = true;
    state->drag_view_context = viewContext;

    bu_vls_free(&path);
    ged_pick_result_free(result);
    emit view_updated(QG_VIEW_REFRESH);
    return true;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
