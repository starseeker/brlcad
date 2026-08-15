/*                     Q S K E T C H . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QSketch.cpp
 *
 * Sketch editing adapter.  The GED session owns all mutable state, librt
 * supplies exact vertices and grouped curve samples, and one retained Obol
 * node per view presents the topology.  Qt and Coin hold no private sketch.
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
#include "rt/functab.h"
#include "rt/geom.h"
#include "rt/primitives/sketch.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgPrimitiveEdit.h"
#include "qtcad/QgSketchEdit.h"
#include "qtcad/QgSignalFlags.h"
#include "../qged_edit_preview_util.h"
#include "QSketch.h"

#include <Inventor/nodes/SoCamera.h>

#include <algorithm>
#include <cmath>
#include <vector>


namespace {

enum QSketchManipulatorAction {
    QSKETCH_MANIPULATOR_PRESS = 0x534b5401u,
    QSKETCH_MANIPULATOR_RELEASE = 0x534b5402u,
    QSKETCH_MANIPULATOR_MOTION = 0x534b5403u
};

static const char *qsketch_manipulator_name = "_sketch_edit_manipulator";

}


struct QSketchManipulatorState {
    QSketch *owner = nullptr;
    struct ged_view_context *view_context = nullptr;
    bobol_display_endpoint_t *endpoint = nullptr;
    BObolViewController *controller = nullptr;
    BObolFeatureHandle feature;
    SoBRLIndexedEditManipulator *node = nullptr;
    int active = -1;
    int active_domain = -1;
    point_t drag_start = VINIT_ZERO;
    point_t drag_anchor = VINIT_ZERO;
};


struct QSketchPresentationState {
    ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
    uint64_t revision = UINT64_MAX;
    std::vector<SbVec3f> points;
    std::vector<int32_t> edges;
    std::vector<int32_t> edge_features;
    int vertex_count = 0;
    int segment_count = 0;
    int select_commands[2] = {0, 0};
    int move_commands[2] = {0, 0};
    int selection_domain = -1;
    int selected_feature = -1;
    std::vector<QSketchManipulatorState *> manipulators;
};


static void
qsketch_manipulator_removed(const BObolCommandResult &result, void *userData)
{
    if (result.status != BObolCommandResultStatus::Removed)
	return;
    QSketchManipulatorState *state =
	static_cast<QSketchManipulatorState *>(userData);
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
}


static void
qsketch_manipulator_clear(QSketchManipulatorState *state)
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
    if (controller && feature.isValid())
	(void)controller->features().remove(feature);
    state->view_context = nullptr;
    state->endpoint = nullptr;
    state->controller = nullptr;
    state->feature = BObolFeatureHandle();
    state->node = nullptr;
}


QSketch::QSketch() : QWidget(), interaction(new QgSketchEdit),
    state(new QSketchPresentationState)
{
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);

    QGroupBox *modeBox = new QGroupBox(tr("Viewport handles"), this);
    QVBoxLayout *modeLayout = new QVBoxLayout(modeBox);
    edit_mode = new QComboBox(modeBox);
    edit_mode->setObjectName(QStringLiteral("sketch.editMode"));
    edit_mode->setProperty("qgTestId", QStringLiteral("sketch.editMode"));
    edit_mode->addItem(tr("Vertex"));
    edit_mode->addItem(tr("Segment"));
    modeLayout->addWidget(edit_mode);
    layout->addWidget(modeBox);

    editor = new QgPrimitiveEdit(this);
    editor->setObjectName(QStringLiteral("sketch.sharedPrimitiveEditor"));
    layout->addWidget(editor);

    connect(edit_mode, QOverload<int>::of(&QComboBox::currentIndexChanged),
	this, &QSketch::mode_changed);
    connect(editor, &QgPrimitiveEdit::targetChanged,
	this, &QSketch::target_changed);
    connect(editor, &QgPrimitiveEdit::sessionEvent,
	this, &QSketch::update_preview);
}


QSketch::~QSketch()
{
    clear_preview();
    delete interaction;
    interaction = nullptr;
    delete state;
    state = nullptr;
}


void
QSketch::setContext(QgPluginContext *ctx)
{
    m_ctx = ctx;
    editor->setGed(getGed());
}


struct ged *
QSketch::getGed() const
{
    return m_ctx ? m_ctx->getGed() : nullptr;
}


void
QSketch::clear_manipulators()
{
    if (!state)
	return;
    for (QSketchManipulatorState *manipulator : state->manipulators) {
	qsketch_manipulator_clear(manipulator);
	delete manipulator;
    }
    state->manipulators.clear();
}


void
QSketch::clear_preview()
{
    if (!state)
	return;
    clear_manipulators();
    state->session = GED_EDIT_SESSION_REF_NULL;
    state->revision = UINT64_MAX;
    state->points.clear();
    state->edges.clear();
    state->edge_features.clear();
    state->vertex_count = 0;
    state->segment_count = 0;
    state->selection_domain = -1;
    state->selected_feature = -1;
    if (interaction)
	interaction->clear();
    preview_path.clear();
}


void
QSketch::target_changed(const QString &path)
{
    if (path != preview_path)
	clear_preview();
}


void
QSketch::mode_changed(int mode)
{
    if (!state || mode < 0 || mode > 1)
	return;
    state->selection_domain = mode;
    state->selected_feature = -1;
    const SoBRLIndexedEditManipulator::Domain domain = mode == 0 ?
	SoBRLIndexedEditManipulator::DOMAIN_VERTEX :
	SoBRLIndexedEditManipulator::DOMAIN_EDGE;
    for (QSketchManipulatorState *manipulator : state->manipulators) {
	if (!manipulator || !manipulator->node)
	    continue;
	manipulator->node->setSelectionDomain(domain);
	manipulator->node->setSelectedIndex(-1);
	if (manipulator->controller)
	    manipulator->controller->requestRender("sketch-edit-mode");
    }
    emit view_updated(QG_VIEW_REFRESH);
}


void
QSketch::update_manipulators(uint64_t revision)
{
    if (!state)
	return;
    const std::vector<struct ged_view_context *> contexts =
	qged_edit_ged_view_contexts(m_ctx);
    for (auto it = state->manipulators.begin();
	it != state->manipulators.end();) {
	QSketchManipulatorState *manipulator = *it;
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
	    qsketch_manipulator_clear(manipulator);
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
	QSketchManipulatorState *manipulator = nullptr;
	for (QSketchManipulatorState *candidate : state->manipulators) {
	    if (candidate && candidate->view_context == viewContext) {
		manipulator = candidate;
		break;
	    }
	}
	if (!manipulator) {
	    manipulator = new QSketchManipulatorState;
	    manipulator->owner = this;
	    manipulator->view_context = viewContext;
	    manipulator->endpoint = endpoint;
	    manipulator->controller = controller;

	    SoBRLIndexedEditManipulator *node =
		new SoBRLIndexedEditManipulator;
	    node->ref();
	    node->manipulatorId = "sketch.topology";
	    BObolFeatureStyle style;
	    style.hasVisible = TRUE;
	    style.visible = TRUE;
	    style.hasSelectable = TRUE;
	    style.selectable = FALSE;
	    BObolFeatureOwner owner;
	    owner.ownerToken = manipulator;
	    owner.ownerId = "qged::sketch-edit";
	    owner.ownerRole = "qged.primitive-manipulator";
	    owner.resultCallback = qsketch_manipulator_removed;
	    owner.callbackUserData = manipulator;
	    manipulator->feature = controller->features().publishCustomNode(
		qsketch_manipulator_name, BObolFeatureScope::Local, node,
		&style, &owner);
	    manipulator->node = manipulator->feature.isValid() ? node : nullptr;
	    node->unref();
	    if (!manipulator->node) {
		qsketch_manipulator_clear(manipulator);
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
		 modifiers, 1200, QSKETCH_MANIPULATOR_PRESS},
		{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0, 0,
		 modifiers, 1200, QSKETCH_MANIPULATOR_RELEASE},
		{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY,
		 BOBOL_INPUT_ANY, 0, modifiers, 1200,
		 QSKETCH_MANIPULATOR_MOTION}
	    };
	    static const BObolInputActionLayer layer = {
		"qged-sketch-edit-manipulator", bindings,
		sizeof(bindings) / sizeof(bindings[0]),
		QSketch::manipulator_input
	    };
	    if (!bobol_display_endpoint_input_action_layer_set(endpoint,
		    &layer, manipulator, manipulator)) {
		qsketch_manipulator_clear(manipulator);
		delete manipulator;
		continue;
	    }
	    state->manipulators.push_back(manipulator);
	}

	manipulator->node->sessionRevision = static_cast<uint32_t>(revision);
	manipulator->node->setTopology(state->points.data(),
	    static_cast<int>(state->points.size()), state->edges.data(),
	    static_cast<int>(state->edge_features.size()), nullptr, nullptr, 0,
	    state->edge_features.data(), state->vertex_count);
	const int mode = edit_mode->currentIndex();
	manipulator->node->setSelectionDomain(mode == 0 ?
	    SoBRLIndexedEditManipulator::DOMAIN_VERTEX :
	    SoBRLIndexedEditManipulator::DOMAIN_EDGE);
	manipulator->node->setSelectedIndex(
	    state->selection_domain == mode ? state->selected_feature : -1);
	if (manipulator->active >= 0)
	    manipulator->node->setActiveIndex(manipulator->active);
	controller->requestRender("sketch-edit-manipulator-update");
	(void)bobol_display_endpoint_request_frame(endpoint,
	    "sketch-edit-manipulator-update");
    }
}


void
QSketch::sync_session_selection(int commandId)
{
    if (!state || ged_edit_session_ref_is_null(editor->session()))
	return;
    int domain = -1;
    for (int i = 0; i < 2; i++) {
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
	    state->select_commands[domain], &values) != GED_EDIT_OK ||
	values.value_count < 1 || !values.value_valid[0])
	return;
    const int feature = static_cast<int>(values.values[0]);
    const int limit = domain == 0 ? state->vertex_count :
	state->segment_count;
    if (feature < 0 || feature >= limit)
	return;
    state->selection_domain = domain;
    state->selected_feature = feature;
    {
	QSignalBlocker blocker(edit_mode);
	edit_mode->setCurrentIndex(domain);
    }
    for (QSketchManipulatorState *manipulator : state->manipulators) {
	if (!manipulator || !manipulator->node)
	    continue;
	manipulator->node->setSelectionDomain(domain == 0 ?
	    SoBRLIndexedEditManipulator::DOMAIN_VERTEX :
	    SoBRLIndexedEditManipulator::DOMAIN_EDGE);
	manipulator->node->setSelectedIndex(feature);
	if (manipulator->controller)
	    manipulator->controller->requestRender(
		"sketch-edit-selection-sync");
    }
}


bool
QSketch::rebuild_preview(uint64_t revision)
{
    struct ged *gedp = getGed();
    const ged_edit_session_ref session = editor->session();
    const QString path = editor->targetPath();
    if (!state || !gedp || !gedp->dbip || path.isEmpty() ||
	ged_edit_session_ref_is_null(session))
	return false;

    if (!interaction || !interaction->refresh(gedp, session, path))
	return false;

    state->points.clear();
    state->edges.clear();
    state->edge_features.clear();
    state->points.reserve(interaction->vertexCount() +
	interaction->linePointCount());
    for (size_t i = 0; i < interaction->vertexCount(); i++) {
	point_t displayPoint;
	if (!interaction->vertexDisplayPoint(i, displayPoint))
	    return false;
	state->points.emplace_back(static_cast<float>(displayPoint[X]),
	    static_cast<float>(displayPoint[Y]),
	    static_cast<float>(displayPoint[Z]));
    }
    for (size_t i = 0; i < interaction->linePointCount(); i++) {
	point_t displayPoint;
	if (!interaction->lineDisplayPoint(i, displayPoint, nullptr, nullptr))
	    return false;
	state->points.emplace_back(static_cast<float>(displayPoint[X]),
	    static_cast<float>(displayPoint[Y]),
	    static_cast<float>(displayPoint[Z]));
    }
    int previous = -1;
    const int lineOffset = static_cast<int>(interaction->vertexCount());
    for (size_t i = 0; i < interaction->linePointCount(); i++) {
	int command = 0;
	int segment = -1;
	point_t ignored;
	(void)interaction->lineDisplayPoint(i, ignored, &command, &segment);
	const int pointIndex = lineOffset + static_cast<int>(i);
	if (command == RT_PRIMITIVE_LINE_MOVE) {
	    previous = pointIndex;
	    continue;
	}
	if (command == RT_PRIMITIVE_LINE_DRAW &&
	    previous >= 0) {
	    state->edges.push_back(previous);
	    state->edges.push_back(pointIndex);
	    state->edge_features.push_back(segment);
	    previous = pointIndex;
	} else {
	    previous = -1;
	}
    }
    state->vertex_count = static_cast<int>(interaction->vertexCount());
    state->segment_count = static_cast<int>(interaction->segmentCount());
    state->session = session;
    state->revision = revision;
    state->select_commands[0] = interaction->selectCommand(
	QgSketchEdit::Vertex);
    state->select_commands[1] = interaction->selectCommand(
	QgSketchEdit::Segment);
    state->move_commands[0] = interaction->moveCommand(QgSketchEdit::Vertex);
    state->move_commands[1] = interaction->moveCommand(QgSketchEdit::Segment);

    preview_path = path;
    update_manipulators(revision);
    struct ged_edit_session_info info = {};
    if (ged_edit_session_info_get(gedp, session, &info) == GED_EDIT_OK &&
	info.command_id)
	sync_session_selection(info.command_id);
    emit view_updated(QG_VIEW_REFRESH);
    return true;
}


void
QSketch::update_preview(int kindValue, qulonglong revisionValue)
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
QSketch::refresh_preview()
{
    editor->refreshFromSession();
    (void)rebuild_preview(0);
}


int
QSketch::manipulator_input(void *userData, BObolInputAction action,
	const BObolInputEvent *event)
{
    QSketchManipulatorState *state =
	static_cast<QSketchManipulatorState *>(userData);
    return state && state->owner ?
	state->owner->handle_manipulator_input(state, action, event) :
	BOBOL_INPUT_RESULT_UNHANDLED;
}


int
QSketch::handle_manipulator_input(QSketchManipulatorState *manipulator,
	BObolInputAction action, const BObolInputEvent *event)
{
    if (!state || !manipulator || !manipulator->node ||
	!manipulator->controller || !event || !interaction ||
	!interaction->isValid())
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
    const SoBRLIndexedEditManipulator::Domain domain = domainIndex == 0 ?
	SoBRLIndexedEditManipulator::DOMAIN_VERTEX :
	SoBRLIndexedEditManipulator::DOMAIN_EDGE;
    if (action == QSKETCH_MANIPULATOR_PRESS) {
	const int feature = manipulator->node->hitTest(domain, event->x,
	    event->y, width, height, camera);
	if (domainIndex < 0 || domainIndex > 1 || feature < 0 ||
	    !state->select_commands[domainIndex] ||
	    !state->move_commands[domainIndex])
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	SbVec3f anchor;
	if (!manipulator->node->featurePosition(domain, feature, anchor) ||
	    ged_edit_session_checkpoint(getGed(), editor->session()) !=
		GED_EDIT_OK)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	if (interaction->selectFeature(domainIndex == 0 ? QgSketchEdit::Vertex :
	    QgSketchEdit::Segment, feature, manipulator->view_context) != GED_EDIT_OK)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	const struct bv *view = bv_context_view_const(
	    ged_view_context_bv_const(manipulator->view_context));
	if (!view || !bv_screen_to_model(manipulator->drag_start, view,
		event->x, event->y))
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	VSET(manipulator->drag_anchor, anchor[0], anchor[1], anchor[2]);
	manipulator->active = feature;
	manipulator->active_domain = domainIndex;
	manipulator->node->setActiveIndex(feature);
	manipulator->node->setHoverIndex(feature);
	manipulator->controller->requestRender("sketch-edit-press");
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (action == QSKETCH_MANIPULATOR_MOTION) {
	if (manipulator->active < 0) {
	    const int hover = manipulator->node->hitTest(domain, event->x,
		event->y, width, height, camera);
	    const bool changed = manipulator->node->hoverIndex.getValue() != hover;
	    manipulator->node->setHoverIndex(hover);
	    if (changed)
		manipulator->controller->requestRender("sketch-edit-hover");
	    return hover < 0 ? BOBOL_INPUT_RESULT_UNHANDLED :
		BOBOL_INPUT_RESULT_HANDLED;
	}
	const int activeDomain = manipulator->active_domain;
	if (activeDomain < 0 || activeDomain > 1)
	    return BOBOL_INPUT_RESULT_HANDLED;
	const struct bv *view = bv_context_view_const(
	    ged_view_context_bv_const(manipulator->view_context));
	point_t current;
	if (!view || !bv_screen_to_model(current, view, event->x, event->y))
	    return BOBOL_INPUT_RESULT_HANDLED;

	point_t secondDisplay;
	vect_t displayDelta;
	if (activeDomain == 0) {
	    VSUB2(displayDelta, current, manipulator->drag_start);
	    VADD2(secondDisplay, manipulator->drag_anchor, displayDelta);
	} else {
	    VSUB2(displayDelta, current, manipulator->drag_start);
	}
	const enum ged_edit_status result = activeDomain == 0 ?
	    interaction->moveVertexToDisplayPoint(secondDisplay,
		manipulator->view_context) :
	    interaction->moveSegmentByDisplayDelta(displayDelta,
		manipulator->view_context);
	if (result == GED_EDIT_OK && activeDomain == 1) {
	    VMOVE(manipulator->drag_start, current);
	} else if (result != GED_EDIT_OK) {
	    manipulator->active = -1;
	    manipulator->active_domain = -1;
	    manipulator->node->setActiveIndex(-1);
	    manipulator->node->setHoverIndex(-1);
	    if (result == GED_EDIT_REJECTED || result == GED_EDIT_ERROR)
		(void)ged_edit_session_revert(getGed(), editor->session());
	    else
		editor->refreshFromSession();
	    manipulator->controller->requestRender("sketch-edit-rejected");
	}
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (action == QSKETCH_MANIPULATOR_RELEASE) {
	if (manipulator->active < 0)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	manipulator->active = -1;
	manipulator->active_domain = -1;
	manipulator->node->setActiveIndex(-1);
	manipulator->node->setHoverIndex(manipulator->node->hitTest(domain,
	    event->x, event->y, width, height, camera));
	manipulator->controller->requestRender("sketch-edit-release");
	return BOBOL_INPUT_RESULT_HANDLED;
    }
    return BOBOL_INPUT_RESULT_UNHANDLED;
}


void
QSketch::sync_selection()
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
	const bool isSketch = !selectedPath.isEmpty() &&
	    db_string_to_path(&fullPath, gedp->dbip,
		pathBytes.constData()) == 0 &&
	    DB_FULL_PATH_CUR_DIR(&fullPath) &&
	    DB_FULL_PATH_CUR_DIR(&fullPath)->d_minor_type ==
		DB5_MINORTYPE_BRLCAD_SKETCH;
	db_free_full_path(&fullPath);
	if (!isSketch)
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
QSketch::reset_for_database()
{
    selection_path.clear();
    clear_preview();
    editor->setGed(nullptr);
    editor->setTargetPath(QString());
    editor->setGed(getGed());
    emit view_updated(QG_VIEW_REFRESH);
}


bool
QSketch::eventFilter(QObject *, QEvent *)
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
