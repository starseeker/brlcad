/*                         Q E L L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2022-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file QEll.cpp
 *
 * Ellipsoid adapter for the shared descriptor-generated primitive editor.
 * Geometry state and transactions live in GED; this class only presents the
 * current session snapshot as retained edit feedback.
 */

#include "common.h"

#include <QEvent>
#include <QVBoxLayout>

#include "BObol/BEditManipulator.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"
#include "bv.h"
#include "ged.h"
#include "ged/edit.h"
#include "ged/plugin/obol.h"
#include "ged/selection.h"
#include "rt/db_fullpath.h"
#include "rt/edit.h"
#include "rt/geom.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgPrimitiveEdit.h"
#include "qtcad/QgSignalFlags.h"
#include "../qged_edit_preview_util.h"
#include "QEll.h"

#include <Inventor/nodes/SoCamera.h>

#include <algorithm>
#include <cmath>


namespace {

enum QEllManipulatorAction {
    QELL_MANIPULATOR_PRESS = 0x454c4c01u,
    QELL_MANIPULATOR_RELEASE = 0x454c4c02u,
    QELL_MANIPULATOR_MOTION = 0x454c4c03u
};

static const char *qell_manipulator_name = "_ell_edit_manipulator";

}


struct QEllManipulatorState {
    QEll *owner = nullptr;
    struct ged_view_context *view_context = nullptr;
    bobol_display_endpoint_t *endpoint = nullptr;
    BObolViewController *controller = nullptr;
    BObolFeatureHandle feature;
    SoBRLEditManipulator *node = nullptr;
    SoBRLEditManipulator::Handle active =
	SoBRLEditManipulator::HANDLE_NONE;
    fastf_t local_lengths[3] = {0.0, 0.0, 0.0};
    int command_ids[3] = {0, 0, 0};
    float drag_base_length = 0.0f;
    float drag_f0 = 0.0f;
    float drag_fx = 0.0f;
    float drag_fy = 0.0f;
};


static void
qell_manipulator_removed(const BObolCommandResult &result, void *user_data)
{
    if (result.status != BObolCommandResultStatus::Removed)
	return;
    QEllManipulatorState *state =
	static_cast<QEllManipulatorState *>(user_data);
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
    state->active = SoBRLEditManipulator::HANDLE_NONE;
}


static void
qell_manipulator_clear(QEllManipulatorState *state)
{
    if (!state)
	return;
    if (state->endpoint)
	(void)bobol_display_endpoint_input_action_layer_clear_if(
	    state->endpoint, state);
    BObolViewController *controller = state->controller;
    const BObolFeatureHandle feature = state->feature;
    state->active = SoBRLEditManipulator::HANDLE_NONE;
    if (controller && feature.isValid())
	(void)controller->features().remove(feature);
    state->view_context = nullptr;
    state->endpoint = nullptr;
    state->controller = nullptr;
    state->feature = BObolFeatureHandle();
    state->node = nullptr;
}


static int
qell_manipulator_command(const struct rt_edit_prim_desc *descriptor,
	const char *parameter_name)
{
    if (!descriptor || !parameter_name)
	return 0;
    for (int i = 0; i < descriptor->ncmd; i++) {
	const struct rt_edit_cmd_desc *command = &descriptor->cmds[i];
	if (command->nparam != 1 || !command->params ||
	    command->params[0].type != RT_EDIT_PARAM_SCALAR ||
	    !command->params[0].name)
	    continue;
	if (BU_STR_EQUAL(command->params[0].name, parameter_name))
	    return command->cmd_id;
    }
    return 0;
}


QEll::QEll()
    : QWidget()
{
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);
    editor = new QgPrimitiveEdit(this);
    editor->setObjectName(QStringLiteral("ell.sharedPrimitiveEditor"));
    layout->addWidget(editor);

    connect(editor, &QgPrimitiveEdit::targetChanged,
	this, &QEll::target_changed);
    connect(editor, &QgPrimitiveEdit::sessionEvent,
	this, &QEll::update_preview);
}


QEll::~QEll()
{
    clear_preview();
}


void
QEll::setContext(QgPluginContext *ctx)
{
    m_ctx = ctx;
    editor->setGed(getGed());
}


struct ged *
QEll::getGed() const
{
    return m_ctx ? m_ctx->getGed() : nullptr;
}


void
QEll::clear_labels()
{
    (void)qged_edit_feature_remove_all(m_ctx, "_ell_edit_labels");
}


void
QEll::clear_preview()
{
    clear_labels();
    clear_manipulator();
    preview_path.clear();
}


void
QEll::clear_manipulator()
{
    for (QEllManipulatorState *state : manipulator_states) {
	qell_manipulator_clear(state);
	delete state;
    }
    manipulator_states.clear();
}


void
QEll::update_manipulator(const point_t center, const vect_t axis_a,
	const vect_t axis_b, const vect_t axis_c,
	const fastf_t local_lengths[3], uint64_t revision)
{
    const std::vector<struct ged_view_context *> contexts =
	qged_edit_ged_view_contexts(m_ctx);
    for (auto it = manipulator_states.begin();
	it != manipulator_states.end();) {
	QEllManipulatorState *state = *it;
	bobol_display_endpoint_t *endpoint = state && state->view_context ?
	    ged_plugin_obol_endpoint_get(state->view_context) : nullptr;
	BObolViewController *controller = state && state->view_context ?
	    ged_plugin_obol_view_controller(state->view_context) : nullptr;
	if (!state || std::find(contexts.begin(), contexts.end(),
		state->view_context) == contexts.end() ||
	    endpoint != state->endpoint || controller != state->controller) {
	    qell_manipulator_clear(state);
	    delete state;
	    it = manipulator_states.erase(it);
	    continue;
	}
	++it;
    }

    const struct rt_edit_prim_desc *descriptor = nullptr;
    int commandIds[3] = {0, 0, 0};
    if (ged_edit_session_descriptor_get(getGed(), editor->session(),
	    &descriptor) == GED_EDIT_OK) {
	commandIds[0] = qell_manipulator_command(descriptor, "a");
	commandIds[1] = qell_manipulator_command(descriptor, "b");
	commandIds[2] = qell_manipulator_command(descriptor, "c");
    }

    for (struct ged_view_context *viewContext : contexts) {
	bobol_display_endpoint_t *endpoint =
	    ged_plugin_obol_endpoint_get(viewContext);
	BObolViewController *controller =
	    ged_plugin_obol_view_controller(viewContext);
	if (!endpoint || !controller)
	    continue;
	QEllManipulatorState *state = nullptr;
	for (QEllManipulatorState *candidate : manipulator_states) {
	    if (candidate && candidate->view_context == viewContext) {
		state = candidate;
		break;
	    }
	}
	if (!state) {
	    state = new QEllManipulatorState;
	    state->owner = this;
	    state->view_context = viewContext;
	    state->endpoint = endpoint;
	    state->controller = controller;

	    SoBRLEditManipulator *node = new SoBRLEditManipulator;
	    node->ref();
	    node->manipulatorId = "ell.axes";
	    BObolFeatureStyle style;
	    style.hasVisible = TRUE;
	    style.visible = TRUE;
	    style.hasSelectable = TRUE;
	    style.selectable = FALSE;
	    BObolFeatureOwner owner;
	    owner.ownerToken = state;
	    owner.ownerId = "qged::ell-edit";
	    owner.ownerRole = "qged.primitive-manipulator";
	    owner.resultCallback = qell_manipulator_removed;
	    owner.callbackUserData = state;
	    state->feature = controller->features().publishCustomNode(
		qell_manipulator_name, BObolFeatureScope::Local, node,
		&style, &owner);
	    state->node = state->feature.isValid() ? node : nullptr;
	    node->unref();
	    if (!state->node) {
		qell_manipulator_clear(state);
		delete state;
		continue;
	    }

	    BObolOverlayInfo overlay;
	    overlay.isOverlay = TRUE;
	    overlay.ownerToken = state;
	    overlay.role = BObolOverlayRole::XRay;
	    overlay.overlayClass = BObolOverlayClass::EditHandle;
	    overlay.lifecycle = BObolOverlayLifecycle::PerTool;
	    overlay.order = BObolOverlayOrder::XRay;
	    overlay.sortOrder = 100;
	    overlay.sourcePath = preview_path.toUtf8().constData();
	    (void)controller->features().setOverlayInfo(state->feature, overlay);

	    static const unsigned int modifiers = BOBOL_INPUT_MOD_SHIFT |
		BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
		BOBOL_INPUT_MOD_META;
	    static const BObolInputBinding bindings[] = {
		{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0, 0,
		 modifiers, 1200, QELL_MANIPULATOR_PRESS},
		{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0, 0,
		 modifiers, 1200, QELL_MANIPULATOR_RELEASE},
		{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY,
		 BOBOL_INPUT_ANY, 0, modifiers, 1200,
		 QELL_MANIPULATOR_MOTION}
	    };
	    static const BObolInputActionLayer layer = {
		"qged-ell-edit-manipulator", bindings,
		sizeof(bindings) / sizeof(bindings[0]),
		QEll::manipulator_input
	    };
	    if (!bobol_display_endpoint_input_action_layer_set(endpoint,
		    &layer, state, state)) {
		qell_manipulator_clear(state);
		delete state;
		continue;
	    }
	    manipulator_states.push_back(state);
	}

	state->local_lengths[0] = local_lengths[0];
	state->local_lengths[1] = local_lengths[1];
	state->local_lengths[2] = local_lengths[2];
	state->command_ids[0] = commandIds[0];
	state->command_ids[1] = commandIds[1];
	state->command_ids[2] = commandIds[2];
	state->node->sessionRevision = static_cast<uint32_t>(revision);
	state->node->setEllipsoidAxes(
	SbVec3f(static_cast<float>(center[0]), static_cast<float>(center[1]),
	    static_cast<float>(center[2])),
	SbVec3f(static_cast<float>(axis_a[0]), static_cast<float>(axis_a[1]),
	    static_cast<float>(axis_a[2])),
	SbVec3f(static_cast<float>(axis_b[0]), static_cast<float>(axis_b[1]),
	    static_cast<float>(axis_b[2])),
	SbVec3f(static_cast<float>(axis_c[0]), static_cast<float>(axis_c[1]),
	    static_cast<float>(axis_c[2])));
	(void)bobol_display_endpoint_request_presentation_frame(endpoint,
	    "ell-edit-manipulator-update");
    }
}


int
QEll::manipulator_input(void *user_data, BObolInputAction action,
	const BObolInputEvent *event)
{
    QEllManipulatorState *state =
	static_cast<QEllManipulatorState *>(user_data);
    return state && state->owner ?
	state->owner->handle_manipulator_input(state, action, event) :
	BOBOL_INPUT_RESULT_UNHANDLED;
}


int
QEll::handle_manipulator_input(QEllManipulatorState *state,
	BObolInputAction action,
	const BObolInputEvent *event)
{
    if (!state || !state->node || !state->controller || !event)
	return BOBOL_INPUT_RESULT_UNHANDLED;
    int width = 0;
    int height = 0;
    const struct bv_context *context =
	reinterpret_cast<const struct bv_context *>(state->view_context);
    if (context) {
	width = bv_context_width_get(context);
	height = bv_context_height_get(context);
    }
    if (width <= 0 || height <= 0) {
	const SbVec2s viewport = state->controller->getViewportRegion().
	    getViewportSizePixels();
	width = static_cast<int>(viewport[0]);
	height = static_cast<int>(viewport[1]);
    }
    SoCamera *camera = state->controller->getCamera();
    if (!camera || width <= 0 || height <= 0)
	return BOBOL_INPUT_RESULT_UNHANDLED;

    if (action == QELL_MANIPULATOR_PRESS) {
	const SoBRLEditManipulator::Handle handle = state->node->hitTest(
	    event->x, event->y, width, height, camera);
	const int index = static_cast<int>(handle);
	if (handle == SoBRLEditManipulator::HANDLE_NONE || index < 0 ||
	    index >= 3 || !state->command_ids[index])
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	float f0 = 0.0f;
	float f10 = 0.0f;
	float f01 = 0.0f;
	if (!state->node->projectedScale(handle, 0, 0, width, height,
		camera, f0) ||
	    !state->node->projectedScale(handle, 1, 0, width, height,
		camera, f10) ||
	    !state->node->projectedScale(handle, 0, 1, width, height,
		camera, f01))
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	if (ged_edit_session_checkpoint(getGed(), editor->session()) !=
	    GED_EDIT_OK)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	state->active = handle;
	state->drag_base_length =
	    static_cast<float>(state->local_lengths[index]);
	state->drag_f0 = f0;
	state->drag_fx = f10 - f0;
	state->drag_fy = f01 - f0;
	state->node->setActiveHandle(handle);
	state->node->setHoverHandle(handle);
	state->controller->requestPresentationRender(
	    "ell-edit-manipulator-press");
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (action == QELL_MANIPULATOR_MOTION) {
	if (state->active == SoBRLEditManipulator::HANDLE_NONE) {
	    const SoBRLEditManipulator::Handle hover = state->node->hitTest(
		event->x, event->y, width, height, camera);
	    const bool changed = state->node->hoverHandle.getValue() !=
		static_cast<int>(hover);
	    state->node->setHoverHandle(hover);
	    if (changed)
		state->controller->requestPresentationRender(
		    "ell-edit-manipulator-hover");
	    return hover == SoBRLEditManipulator::HANDLE_NONE ?
		BOBOL_INPUT_RESULT_UNHANDLED : BOBOL_INPUT_RESULT_HANDLED;
	}

	const int index = static_cast<int>(state->active);
	const float factor = std::max(1.0e-9f,
	    state->drag_f0 + state->drag_fx * static_cast<float>(event->x) +
	    state->drag_fy * static_cast<float>(event->y));
	fastf_t value[1] = {
	    static_cast<fastf_t>(state->drag_base_length * factor)
	};
	struct ged_edit_command_input input = {};
	input.command_id = state->command_ids[index];
	input.values = value;
	input.value_count = 1;
	input.view = state->view_context;
	const enum ged_edit_status result = ged_edit_session_apply(getGed(),
	    editor->session(), &input);
	if (result != GED_EDIT_OK) {
	    /* Press established a checkpoint.  A primitive handler rejection may
	     * have touched private edit bookkeeping even when geometry validation
	     * fails, so terminate the drag and restore the authoritative session. */
	    state->active = SoBRLEditManipulator::HANDLE_NONE;
	    state->node->setActiveHandle(SoBRLEditManipulator::HANDLE_NONE);
	    state->node->setHoverHandle(SoBRLEditManipulator::HANDLE_NONE);
	    if (result == GED_EDIT_REJECTED || result == GED_EDIT_ERROR)
		(void)ged_edit_session_revert(getGed(), editor->session());
	    else
		editor->refreshFromSession();
	    state->controller->requestLodCapacityRender(
		"ell-edit-manipulator-rejected");
	}
	return BOBOL_INPUT_RESULT_HANDLED;
    }

    if (action == QELL_MANIPULATOR_RELEASE) {
	if (state->active == SoBRLEditManipulator::HANDLE_NONE)
	    return BOBOL_INPUT_RESULT_UNHANDLED;
	state->active = SoBRLEditManipulator::HANDLE_NONE;
	state->node->setActiveHandle(SoBRLEditManipulator::HANDLE_NONE);
	state->node->setHoverHandle(state->node->hitTest(event->x, event->y,
	    width, height, camera));
	state->controller->requestPresentationRender(
	    "ell-edit-manipulator-release");
	return BOBOL_INPUT_RESULT_HANDLED;
    }
    return BOBOL_INPUT_RESULT_UNHANDLED;
}


void
QEll::target_changed(const QString &path)
{
    if (path != preview_path)
	clear_preview();
}


void
QEll::refresh_preview()
{
    editor->refreshFromSession();
    update_preview(static_cast<int>(GED_EDIT_SESSION_UPDATE), 0);
}


void
QEll::update_preview(int kindValue, qulonglong UNUSED(revision))
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
    const QString path = editor->targetPath();
    if (!gedp || !gedp->dbip || path.isEmpty() ||
	ged_edit_session_ref_is_null(session))
	return;

    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    if (ged_edit_session_internal_copy(gedp, session, &intern) != GED_EDIT_OK ||
	intern.idb_type != ID_ELL || !intern.idb_ptr) {
	rt_db_free_internal(&intern);
	return;
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
	return;
    }

    struct rt_ell_internal *source =
	static_cast<struct rt_ell_internal *>(intern.idb_ptr);
    RT_ELL_CK_MAGIC(source);
    struct rt_ell_internal display = *source;
    MAT4X3PNT(display.v, pathMatrix, source->v);
    MAT4X3VEC(display.a, pathMatrix, source->a);
    MAT4X3VEC(display.b, pathMatrix, source->b);
    MAT4X3VEC(display.c, pathMatrix, source->c);

    struct bn_tol tolerance = BN_TOL_INIT_TOL;

    struct rt_db_internal displayIntern = intern;
    displayIntern.idb_ptr = &display;
    struct rt_point_labels pointLabels[9];
    int labelCount = 0;
    if (displayIntern.idb_meth && displayIntern.idb_meth->ft_labels)
	labelCount = displayIntern.idb_meth->ft_labels(pointLabels, 8,
	    bn_mat_identity, &displayIntern, &tolerance);
    const unsigned char labelColor[3] = {255, 255, 0};
    (void)qged_edit_feature_labels_replace_all(m_ctx, "_ell_edit_labels",
	pointLabels, labelCount, labelColor);

    preview_path = path;
    const fastf_t localLengths[3] = {
	MAGNITUDE(source->a), MAGNITUDE(source->b), MAGNITUDE(source->c)
    };
    struct ged_edit_session_info sessionInfo = {};
    const uint64_t sessionRevision = ged_edit_session_info_get(gedp,
	session, &sessionInfo) == GED_EDIT_OK ? sessionInfo.revision : 0;
    update_manipulator(display.v, display.a, display.b, display.c,
	localLengths, sessionRevision);
    rt_db_free_internal(&intern);
    emit view_updated(QG_VIEW_REFRESH);
}


void
QEll::sync_selection()
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
	const bool isEll = !selectedPath.isEmpty() &&
	    db_string_to_path(&fullPath, gedp->dbip,
		pathBytes.constData()) == 0 &&
	    DB_FULL_PATH_CUR_DIR(&fullPath) &&
	    DB_FULL_PATH_CUR_DIR(&fullPath)->d_minor_type ==
		DB5_MINORTYPE_BRLCAD_ELL;
	db_free_full_path(&fullPath);
	if (!isEll)
	    selectedPath.clear();
    }

    if (selectedPath.isEmpty()) {
	/* Preserve a manually entered target across unrelated selection changes. */
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
QEll::reset_for_database()
{
    selection_path.clear();
    clear_preview();
    editor->setGed(nullptr);
    editor->setTargetPath(QString());
    editor->setGed(getGed());
    emit view_updated(QG_VIEW_REFRESH);
}


bool
QEll::eventFilter(QObject *, QEvent *)
{
    return false;
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
