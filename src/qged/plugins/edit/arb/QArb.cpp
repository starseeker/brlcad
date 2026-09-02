/*                         Q A R B . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QArb.cpp
 *
 * ARB editing adapter.  Librt supplies subtype-specific legal topology, GED
 * owns the sole mutable session, and one retained Obol node per view presents
 * all indexed handles.  No widget or scene node owns a primitive copy.
 */

#include "common.h"

#include <QComboBox>
#include <QEvent>
#include <QGroupBox>
#include <QSignalBlocker>
#include <QVBoxLayout>

#include "BObol/BEditManipulator.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "bn/mat.h"
#include "bv.h"
#include "ged.h"
#include "ged/edit.h"
#include "ged/plugin/obol.h"
#include "ged/selection.h"
#include "rt/db_fullpath.h"
#include "rt/edit.h"
#include "rt/geom.h"
#include "rt/primitives/arb8.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgPrimitiveEdit.h"
#include "qtcad/QgSignalFlags.h"
#include "../qged_edit_preview_util.h"
#include "QArb.h"

#include <Inventor/nodes/SoCamera.h>

#include <algorithm>
#include <vector>


namespace {

enum QArbManipulatorAction {
    QARB_MANIPULATOR_PRESS = 0x41524201u,
    QARB_MANIPULATOR_RELEASE = 0x41524202u,
    QARB_MANIPULATOR_MOTION = 0x41524203u
};

static const char *qarb_label_name = "_arb_edit_labels";
static const char *qarb_manipulator_name = "_arb_edit_manipulator";

}


struct QArbManipulatorState {
    QArb *owner = nullptr;
    struct ged_view_context *view_context = nullptr;
    bobol_display_endpoint_t *endpoint = nullptr;
    BObolViewController *controller = nullptr;
    BObolFeatureHandle feature;
    SoBRLIndexedEditManipulator *node = nullptr;
    int active = -1;
    int active_domain = -1;
    int active_edit_index = -1;
    point_t drag_start = VINIT_ZERO;
    point_t drag_anchor = VINIT_ZERO;
};


struct QArbPresentationState {
    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    uint64_t revision = UINT64_MAX;
    struct rt_arb_edit_topology topology = {};
    std::vector<SbVec3f> points;
    std::vector<int32_t> edges;
    std::vector<int> edge_topology_indices;
    std::vector<int32_t> faces;
    std::vector<int32_t> face_counts;
    std::vector<int> face_topology_indices;
    mat_t path_matrix = MAT_INIT_IDN;
    mat_t inverse_path_matrix = MAT_INIT_IDN;
    bool matrix_valid = false;
    int select_commands[3] = {0, 0, 0};
    int move_commands[3] = {0, 0, 0};
    int selection_domain = -1;
    int selected_feature = -1;
    std::vector<QArbManipulatorState *> manipulators;
};


static void
qarb_manipulator_removed(const BObolCommandResult &result, void *userData)
{
    if (result.status != BObolCommandResultStatus::Removed)
	return;
    QArbManipulatorState *state =
	static_cast<QArbManipulatorState *>(userData);
    if (!state)
	return;
    if (state->endpoint)
	(void)bobol_display_endpoint_input_action_layer_clear_if(
	    state->endpoint, state);
    state->view_context = nullptr;
    state->endpoint = nullptr;
    state->controller = nullptr;
    state->feature = BObolFeatureHandle();
    state->node = nullptr;
    state->active = -1;
    state->active_domain = -1;
    state->active_edit_index = -1;
}


static void
qarb_manipulator_clear(QArbManipulatorState *state)
{
    if (!state)
	return;
    if (state->endpoint)
	(void)bobol_display_endpoint_input_action_layer_clear_if(
	    state->endpoint, state);
    BObolViewController *controller = state->controller;
    const BObolFeatureHandle feature = state->feature;
    state->active = -1;
    state->active_domain = -1;
    state->active_edit_index = -1;
    if (controller && feature.isValid())
	(void)controller->features().remove(feature);
    state->view_context = nullptr;
    state->endpoint = nullptr;
    state->controller = nullptr;
    state->feature = BObolFeatureHandle();
    state->node = nullptr;
}


QArb::QArb() : QWidget(), state(new QArbPresentationState)
{
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);

    QGroupBox *modeBox = new QGroupBox(tr("Viewport handles"), this);
    QVBoxLayout *modeLayout = new QVBoxLayout(modeBox);
    edit_mode = new QComboBox(modeBox);
    edit_mode->setObjectName(QStringLiteral("arb.editMode"));
    edit_mode->setProperty("qgTestId", QStringLiteral("arb.editMode"));
    edit_mode->addItem(tr("Vertex"));
    edit_mode->addItem(tr("Edge"));
    edit_mode->addItem(tr("Face"));
    modeLayout->addWidget(edit_mode);
    layout->addWidget(modeBox);

    editor = new QgPrimitiveEdit(this);
    editor->setObjectName(QStringLiteral("arb.sharedPrimitiveEditor"));
    layout->addWidget(editor);

    connect(edit_mode, QOverload<int>::of(&QComboBox::currentIndexChanged),
	this, &QArb::mode_changed);
    connect(editor, &QgPrimitiveEdit::targetChanged,
	this, &QArb::target_changed);
    connect(editor, &QgPrimitiveEdit::sessionEvent,
	this, &QArb::update_preview);
}


QArb::~QArb()
{
    clear_preview();
    delete state;
    state = nullptr;
}


void
QArb::setContext(QgPluginContext *ctx)
{
    m_ctx = ctx;
    editor->setGed(getGed());
}


struct ged *
QArb::getGed() const
{
    return m_ctx ? m_ctx->getGed() : nullptr;
}


int
QArb::command_id(const char *name) const
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
QArb::clear_labels()
{
    (void)qged_edit_feature_remove_all(m_ctx, qarb_label_name);
}


void
QArb::clear_manipulators()
{
    if (!state)
	return;
    for (QArbManipulatorState *manipulator : state->manipulators) {
	qarb_manipulator_clear(manipulator);
	delete manipulator;
    }
    state->manipulators.clear();
}


void
QArb::clear_preview()
{
    if (!state)
	return;
    clear_labels();
    clear_manipulators();
    state->session = GED_EDIT_SESSION_REF_NULL;
    state->revision = UINT64_MAX;
    state->points.clear();
    state->edges.clear();
    state->edge_topology_indices.clear();
    state->faces.clear();
    state->face_counts.clear();
    state->face_topology_indices.clear();
    state->matrix_valid = false;
    state->selection_domain = -1;
    state->selected_feature = -1;
    preview_path.clear();
}


void
QArb::target_changed(const QString &path)
{
    if (path != preview_path)
	clear_preview();
}


void
QArb::mode_changed(int mode)
{
    if (!state || mode < 0 || mode > 2)
	return;
    state->selection_domain = mode;
    state->selected_feature = -1;
    const SoBRLIndexedEditManipulator::Domain domain =
	static_cast<SoBRLIndexedEditManipulator::Domain>(mode + 1);
    for (QArbManipulatorState *manipulator : state->manipulators) {
	if (!manipulator || !manipulator->node)
	    continue;
	manipulator->node->setSelectionDomain(domain);
	manipulator->node->setSelectedIndex(-1);
	if (manipulator->controller)
	    manipulator->controller->requestPresentationRender("arb-edit-mode");
    }
    emit view_updated(QG_VIEW_REFRESH);
}


void
QArb::update_manipulators(uint64_t revision)
{
    if (!state)
	return;
    const std::vector<struct ged_view_context *> contexts =
	qged_edit_ged_view_contexts(m_ctx);
    for (auto it = state->manipulators.begin();
	it != state->manipulators.end();) {
	QArbManipulatorState *manipulator = *it;
	bobol_display_endpoint_t *endpoint = manipulator &&
	    manipulator->view_context ?
	    ged_plugin_obol_endpoint_get(manipulator->view_context) : nullptr;
	BObolViewController *controller = manipulator &&
	    manipulator->view_context ?
	    ged_plugin_obol_view_controller(manipulator->view_context) : nullptr;
	if (!manipulator || std::find(contexts.begin(), contexts.end(),
		manipulator->view_context) == contexts.end() ||
	    endpoint != manipulator->endpoint ||
	    controller != manipulator->controller) {
	    qarb_manipulator_clear(manipulator);
	    delete manipulator;
	    it = state->manipulators.erase(it);
	    continue;
	}
	++it;
    }

    for (struct ged_view_context *viewContext : contexts) {
	bobol_display_endpoint_t *endpoint =
	    ged_plugin_obol_endpoint_get(viewContext);
	BObolViewController *controller =
	    ged_plugin_obol_view_controller(viewContext);
	if (!endpoint || !controller)
	    continue;
	QArbManipulatorState *manipulator = nullptr;
	for (QArbManipulatorState *candidate : state->manipulators) {
	    if (candidate && candidate->view_context == viewContext) {
		manipulator = candidate;
		break;
	    }
	}
	if (!manipulator) {
	    manipulator = new QArbManipulatorState;
	    manipulator->owner = this;
	    manipulator->view_context = viewContext;
	    manipulator->endpoint = endpoint;
	    manipulator->controller = controller;

	    SoBRLIndexedEditManipulator *node =
		new SoBRLIndexedEditManipulator;
	    node->ref();
	    node->manipulatorId = "arb.topology";
	    BObolFeatureStyle style;
	    style.hasVisible = TRUE;
	    style.visible = TRUE;
	    style.hasSelectable = TRUE;
	    style.selectable = FALSE;
	    BObolFeatureOwner owner;
	    owner.ownerToken = manipulator;
	    owner.ownerId = "qged::arb-edit";
	    owner.ownerRole = "qged.primitive-manipulator";
	    owner.resultCallback = qarb_manipulator_removed;
	    owner.callbackUserData = manipulator;
	    manipulator->feature = controller->features().publishCustomNode(
		qarb_manipulator_name, BObolFeatureScope::Local, node,
		&style, &owner);
	    manipulator->node = manipulator->feature.isValid() ? node : nullptr;
	    node->unref();
	    if (!manipulator->node) {
		qarb_manipulator_clear(manipulator);
		delete manipulator;
		continue;
	    }

	    BObolOverlayInfo overlay;
	    overlay.isOverlay = TRUE;
	    overlay.ownerToken = manipulator;
	    overlay.role = BObolOverlayRole::XRay;
	    overlay.overlayClass = BObolOverlayClass::EditHandle;
	    overlay.lifecycle = BObolOverlayLifecycle::PerTool;
	    overlay.order = BObolOverlayOrder::XRay;
	    overlay.sortOrder = 100;
	    overlay.sourcePath = preview_path.toUtf8().constData();
	    (void)controller->features().setOverlayInfo(
		manipulator->feature, overlay);

	    static const unsigned int modifiers = BOBOL_INPUT_MOD_SHIFT |
		BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
		BOBOL_INPUT_MOD_META;
	    static const BObolInputBinding bindings[] = {
		{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0, 0,
		 modifiers, 1200, QARB_MANIPULATOR_PRESS},
		{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0, 0,
		 modifiers, 1200, QARB_MANIPULATOR_RELEASE},
		{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY,
		 BOBOL_INPUT_ANY, 0, modifiers, 1200,
		 QARB_MANIPULATOR_MOTION}
	    };
	    static const BObolInputActionLayer layer = {
		"qged-arb-edit-manipulator", bindings,
		sizeof(bindings) / sizeof(bindings[0]),
		QArb::manipulator_input
	    };
	    if (!bobol_display_endpoint_input_action_layer_set(endpoint,
		    &layer, manipulator, manipulator)) {
		qarb_manipulator_clear(manipulator);
		delete manipulator;
		continue;
	    }
	    state->manipulators.push_back(manipulator);
	}

	manipulator->node->sessionRevision = static_cast<uint32_t>(revision);
	manipulator->node->setTopology(state->points.data(),
	    static_cast<int>(state->points.size()), state->edges.data(),
	    static_cast<int>(state->edge_topology_indices.size()),
	    state->faces.data(), state->face_counts.data(),
	    static_cast<int>(state->face_topology_indices.size()));
	const int mode = edit_mode->currentIndex();
	manipulator->node->setSelectionDomain(
	    static_cast<SoBRLIndexedEditManipulator::Domain>(mode + 1));
	manipulator->node->setSelectedIndex(
	    state->selection_domain == mode ? state->selected_feature : -1);
	if (manipulator->active >= 0)
	    manipulator->node->setActiveIndex(manipulator->active);
	(void)bobol_display_endpoint_request_presentation_frame(endpoint,
	    "arb-edit-manipulator-update");
    }
}


void
QArb::sync_session_selection(int commandId)
{
    if (!state || ged_edit_session_ref_is_null(editor->session()))
	return;
    int domain = -1;
    for (int i = 0; i < 3; i++) {
	if (commandId == state->select_commands[i] ||
	    commandId == state->move_commands[i]) {
	    domain = i;
	    break;
	}
    }
    if (domain < 0)
	return;
    struct rt_edit_cmd_values values;
    if (ged_edit_session_command_values_get(getGed(), editor->session(),
	    commandId, &values) != GED_EDIT_OK || values.value_count < 1 ||
	!values.value_valid[0])
	return;
    const int editIndex = static_cast<int>(values.values[0]);
    int feature = -1;
    if (domain == 0) {
	for (int i = 0; i < state->topology.vertex_count; i++) {
	    if (state->topology.vertices[i].edit_index == editIndex) {
		feature = i;
		break;
	    }
	}
    } else if (domain == 1) {
	for (size_t i = 0; i < state->edge_topology_indices.size(); i++) {
	    const int topologyIndex = state->edge_topology_indices[i];
	    if (state->topology.edges[topologyIndex].edit_index == editIndex) {
		feature = static_cast<int>(i);
		break;
	    }
	}
    } else {
	for (size_t i = 0; i < state->face_topology_indices.size(); i++) {
	    const int topologyIndex = state->face_topology_indices[i];
	    if (state->topology.faces[topologyIndex].edit_index == editIndex) {
		feature = static_cast<int>(i);
		break;
	    }
	}
    }
    if (feature < 0)
	return;
    state->selection_domain = domain;
    state->selected_feature = feature;
    {
	QSignalBlocker blocker(edit_mode);
	edit_mode->setCurrentIndex(domain);
    }
    for (QArbManipulatorState *manipulator : state->manipulators) {
	if (!manipulator || !manipulator->node)
	    continue;
	manipulator->node->setSelectionDomain(
	    static_cast<SoBRLIndexedEditManipulator::Domain>(domain + 1));
	manipulator->node->setSelectedIndex(feature);
	if (manipulator->controller)
	    manipulator->controller->requestPresentationRender(
		"arb-edit-selection-sync");
    }
}


bool
QArb::rebuild_preview(uint64_t revision)
{
    struct ged *gedp = getGed();
    const ged_edit_session_ref session = editor->session();
    const QString path = editor->targetPath();
    if (!state || !gedp || !gedp->dbip || path.isEmpty() ||
	ged_edit_session_ref_is_null(session))
	return false;

    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    if (ged_edit_session_internal_copy(gedp, session, &intern) != GED_EDIT_OK ||
	intern.idb_type != ID_ARB8 || !intern.idb_ptr) {
	rt_db_free_internal(&intern);
	return false;
    }
    struct rt_arb_internal *arb =
	static_cast<struct rt_arb_internal *>(intern.idb_ptr);
    RT_ARB_CK_MAGIC(arb);
    struct bn_tol tolerance = BN_TOL_INIT_TOL;
    if (rt_arb_edit_topology_get(&state->topology, &intern, &tolerance) !=
	BRLCAD_OK) {
	rt_db_free_internal(&intern);
	return false;
    }

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

    state->points.resize(static_cast<size_t>(state->topology.vertex_count));
    for (int i = 0; i < state->topology.vertex_count; i++) {
	point_t displayPoint;
	MAT4X3PNT(displayPoint, pathMatrix,
	    arb->pt[state->topology.vertices[i].edit_index]);
	state->points[static_cast<size_t>(i)] = SbVec3f(
	    static_cast<float>(displayPoint[X]),
	    static_cast<float>(displayPoint[Y]),
	    static_cast<float>(displayPoint[Z]));
    }
    state->edges.clear();
    state->edge_topology_indices.clear();
    for (int i = 0; i < state->topology.edge_count; i++) {
	const struct rt_arb_edit_edge *edge = &state->topology.edges[i];
	if (edge->edit_index < 0)
	    continue;
	state->edges.push_back(edge->vertices[0]);
	state->edges.push_back(edge->vertices[1]);
	state->edge_topology_indices.push_back(i);
    }
    state->faces.clear();
    state->face_counts.clear();
    state->face_topology_indices.clear();
    for (int i = 0; i < state->topology.face_count; i++) {
	const struct rt_arb_edit_face *face = &state->topology.faces[i];
	if (!face->movable)
	    continue;
	state->face_counts.push_back(face->vertex_count);
	for (int j = 0; j < face->vertex_count; j++)
	    state->faces.push_back(face->vertices[j]);
	state->face_topology_indices.push_back(i);
    }
    MAT_COPY(state->path_matrix, pathMatrix);
    bn_mat_inv(state->inverse_path_matrix, pathMatrix);
    state->matrix_valid = true;
    state->session = session;
    state->revision = revision;
    state->select_commands[0] = command_id("ECMD_ARB_SELECT_VERTEX");
    state->select_commands[1] = command_id("ECMD_ARB_SELECT_EDGE");
    state->select_commands[2] = command_id("ECMD_ARB_SELECT_FACE");
    state->move_commands[0] = command_id("ECMD_ARB_MOVE_VERTEX");
    state->move_commands[1] = command_id("ECMD_ARB_MOVE_EDGE");
    state->move_commands[2] = command_id("ECMD_ARB_MOVE_FACE");

    struct rt_arb_internal display = *arb;
    for (int i = 0; i < 8; i++)
	MAT4X3PNT(display.pt[i], pathMatrix, arb->pt[i]);
    struct rt_db_internal displayIntern = intern;
    displayIntern.idb_ptr = &display;
    struct rt_point_labels pointLabels[9];
    int labelCount = 0;
    if (displayIntern.idb_meth && displayIntern.idb_meth->ft_labels)
	labelCount = displayIntern.idb_meth->ft_labels(pointLabels, 8,
	    bn_mat_identity, &displayIntern, &tolerance);
    const unsigned char labelColor[3] = {255, 255, 0};
    (void)qged_edit_feature_labels_replace_all(m_ctx, qarb_label_name,
	pointLabels, labelCount, labelColor);

    preview_path = path;
    update_manipulators(revision);
    struct ged_edit_session_info info = {};
    if (ged_edit_session_info_get(gedp, session, &info) == GED_EDIT_OK &&
	info.command_id)
	sync_session_selection(info.command_id);
    rt_db_free_internal(&intern);
    emit view_updated(QG_VIEW_REFRESH);
    return true;
}


void
QArb::update_preview(int kindValue, qulonglong revisionValue)
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
    (void)rebuild_preview(static_cast<uint64_t>(revisionValue));
}


void
QArb::refresh_preview()
{
    editor->refreshFromSession();
    (void)rebuild_preview(0);
}


int
QArb::manipulator_input(void *userData, BObolInputAction action,
	const BObolInputEvent *event)
{
    QArbManipulatorState *state =
	static_cast<QArbManipulatorState *>(userData);
    return state && state->owner ?
	state->owner->handle_manipulator_input(state, action, event) :
	BOBOL_INPUT_RESULT_UNHANDLED;
}


int
QArb::handle_manipulator_input(QArbManipulatorState *manipulator,
	BObolInputAction action, const BObolInputEvent *event)
{
    if (!state || !manipulator || !manipulator->node ||
	!manipulator->controller || !event || !state->matrix_valid)
	return BOBOL_INPUT_RESULT_UNHANDLED;
    int width = 0;
    int height = 0;
    const struct bv_context *context =
	reinterpret_cast<const struct bv_context *>(manipulator->view_context);
    if (context) {
	width = bv_context_width_get(context);
	height = bv_context_height_get(context);
    }
    if (width <= 0 || height <= 0) {
	const SbVec2s viewport = manipulator->controller->getViewportRegion().
	    getViewportSizePixels();
	width = static_cast<int>(viewport[0]);
	height = static_cast<int>(viewport[1]);
    }
    SoCamera *camera = manipulator->controller->getCamera();
    if (!camera || width <= 0 || height <= 0)
	return BOBOL_INPUT_RESULT_UNHANDLED;

    const int domainIndex = edit_mode->currentIndex();
    const SoBRLIndexedEditManipulator::Domain domain =
	static_cast<SoBRLIndexedEditManipulator::Domain>(domainIndex + 1);
    if (action == QARB_MANIPULATOR_PRESS) {
	const int feature = manipulator->node->hitTest(domain, event->x,
	    event->y, width, height, camera);
	if (domainIndex < 0 || domainIndex > 2 || feature < 0 ||
	    !state->select_commands[domainIndex] ||
	    !state->move_commands[domainIndex])
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	int editIndex = -1;
	if (domainIndex == 0) {
	    editIndex = state->topology.vertices[feature].edit_index;
	} else if (domainIndex == 1) {
	    if (static_cast<size_t>(feature) >=
		state->edge_topology_indices.size())
		return BOBOL_INPUT_RESULT_UNHANDLED;
	    editIndex = state->topology.edges[
		state->edge_topology_indices[feature]].edit_index;
	} else {
	    if (static_cast<size_t>(feature) >=
		state->face_topology_indices.size())
		return BOBOL_INPUT_RESULT_UNHANDLED;
	    editIndex = state->topology.faces[
		state->face_topology_indices[feature]].edit_index;
	}
	SbVec3f anchor;
	if (editIndex < 0 ||
	    !manipulator->node->featurePosition(domain, feature, anchor) ||
	    ged_edit_session_checkpoint(getGed(), editor->session()) !=
		GED_EDIT_OK)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	fastf_t selectionValue[1] = {static_cast<fastf_t>(editIndex)};
	struct ged_edit_command_input selection = {};
	selection.command_id = state->select_commands[domainIndex];
	selection.values = selectionValue;
	selection.value_count = 1;
	selection.view = manipulator->view_context;
	if (ged_edit_session_apply(getGed(), editor->session(), &selection) !=
	    GED_EDIT_OK)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	const struct bv *view = bv_context_view_const(
	    ged_view_context_bv_const(manipulator->view_context));
	if (!view || !bv_screen_to_model(manipulator->drag_start, view,
		event->x, event->y))
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	VSET(manipulator->drag_anchor, anchor[0], anchor[1], anchor[2]);
	manipulator->active = feature;
	manipulator->active_domain = domainIndex;
	manipulator->active_edit_index = editIndex;
	manipulator->node->setActiveIndex(feature);
	manipulator->node->setHoverIndex(feature);
	manipulator->controller->requestPresentationRender(
	    "arb-edit-manipulator-press");
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (action == QARB_MANIPULATOR_MOTION) {
	if (manipulator->active < 0) {
	    const int hover = manipulator->node->hitTest(domain, event->x,
		event->y, width, height, camera);
	    const bool changed = manipulator->node->hoverIndex.getValue() !=
		hover;
	    manipulator->node->setHoverIndex(hover);
	    if (changed)
		manipulator->controller->requestPresentationRender(
		    "arb-edit-manipulator-hover");
	    return hover < 0 ? BOBOL_INPUT_RESULT_UNHANDLED :
		BOBOL_INPUT_RESULT_HANDLED;
	}
	const int activeDomain = manipulator->active_domain;
	if (activeDomain < 0 || activeDomain > 2)
	    return BOBOL_INPUT_RESULT_HANDLED;
	const struct bv *view = bv_context_view_const(
	    ged_view_context_bv_const(manipulator->view_context));
	point_t current;
	if (!view || !bv_screen_to_model(current, view, event->x, event->y))
	    return BOBOL_INPUT_RESULT_HANDLED;
	vect_t delta;
	VSUB2(delta, current, manipulator->drag_start);
	point_t displayTarget;
	VADD2(displayTarget, manipulator->drag_anchor, delta);
	point_t sourceTarget;
	MAT4X3PNT(sourceTarget, state->inverse_path_matrix, displayTarget);
	const fastf_t base2local = getGed()->dbip->dbi_base2local;
	fastf_t values[4] = {
	    static_cast<fastf_t>(manipulator->active_edit_index),
	    sourceTarget[X] * base2local,
	    sourceTarget[Y] * base2local,
	    sourceTarget[Z] * base2local
	};
	struct ged_edit_command_input input = {};
	input.command_id = state->move_commands[activeDomain];
	input.values = values;
	input.value_count = 4;
	input.view = manipulator->view_context;
	const enum ged_edit_status result = ged_edit_session_apply(getGed(),
	    editor->session(), &input);
	if (result != GED_EDIT_OK) {
	    manipulator->active = -1;
	    manipulator->active_domain = -1;
	    manipulator->active_edit_index = -1;
	    manipulator->node->setActiveIndex(-1);
	    manipulator->node->setHoverIndex(-1);
	    if (result == GED_EDIT_REJECTED || result == GED_EDIT_ERROR)
		(void)ged_edit_session_revert(getGed(), editor->session());
	    else
		editor->refreshFromSession();
	    manipulator->controller->requestLodCapacityRender(
		"arb-edit-manipulator-rejected");
	}
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (action == QARB_MANIPULATOR_RELEASE) {
	if (manipulator->active < 0)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	manipulator->active = -1;
	manipulator->active_domain = -1;
	manipulator->active_edit_index = -1;
	manipulator->node->setActiveIndex(-1);
	manipulator->node->setHoverIndex(manipulator->node->hitTest(domain,
	    event->x, event->y, width, height, camera));
	manipulator->controller->requestPresentationRender(
	    "arb-edit-manipulator-release");
	return BOBOL_INPUT_RESULT_HANDLED;
    }
    return BOBOL_INPUT_RESULT_UNHANDLED;
}


void
QArb::sync_selection()
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
	const bool isArb = !selectedPath.isEmpty() &&
	    db_string_to_path(&fullPath, gedp->dbip,
		pathBytes.constData()) == 0 &&
	    DB_FULL_PATH_CUR_DIR(&fullPath) &&
	    DB_FULL_PATH_CUR_DIR(&fullPath)->d_minor_type ==
		DB5_MINORTYPE_BRLCAD_ARB8;
	db_free_full_path(&fullPath);
	if (!isArb)
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
QArb::reset_for_database()
{
    selection_path.clear();
    clear_preview();
    editor->setGed(nullptr);
    editor->setTargetPath(QString());
    editor->setGed(getGed());
    emit view_updated(QG_VIEW_REFRESH);
}


bool
QArb::eventFilter(QObject *, QEvent *)
{
    return false;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
