/*                   Q G C A N V A S I N P U T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2021-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgCanvasInput.cpp
 *
 * Qt event normalization and endpoint-local BRL-CAD view actions.
 */

#include "common.h"

#include <QtGlobal>

#include <cmath>

extern "C" {
#include "bn/str.h"
}

#include "qtcad/defines.h"
#include "QgCanvasInput.h"
#include "BObol/BDisplayEndpoint.h"
#include "bv.h"

struct QgCanvasInput::Impl {
    BObolInputContext context;
    bobol_display_endpoint_t *endpoint = NULL;
    struct bv_context *dispatch_view = NULL;
    double press_x = 0.0;
    double press_y = 0.0;
    int mouse_mode = BV_ADJUST_SCALE;
};

static unsigned int
qgcanvasinput_modifiers(Qt::KeyboardModifiers modifiers)
{
    unsigned int result = BOBOL_INPUT_MOD_NONE;
    if (modifiers.testFlag(Qt::ShiftModifier))
	result |= BOBOL_INPUT_MOD_SHIFT;
    if (modifiers.testFlag(Qt::ControlModifier))
	result |= BOBOL_INPUT_MOD_CONTROL;
    if (modifiers.testFlag(Qt::AltModifier))
	result |= BOBOL_INPUT_MOD_ALT;
    if (modifiers.testFlag(Qt::MetaModifier))
	result |= BOBOL_INPUT_MOD_META;
    return result;
}

static int
qgcanvasinput_key(int key)
{
    const int first_function_key = static_cast<int>(Qt::Key_F1);
    const int last_function_key = static_cast<int>(Qt::Key_F12);
    if (key >= first_function_key && key <= last_function_key)
	return BOBOL_INPUT_KEY_F1 + (key - first_function_key);
    if (key == static_cast<int>(Qt::Key_Shift))
	return BOBOL_INPUT_KEY_SHIFT;
    if (key == static_cast<int>(Qt::Key_Control))
	return BOBOL_INPUT_KEY_CONTROL;
    if (key == static_cast<int>(Qt::Key_Alt))
	return BOBOL_INPUT_KEY_ALT;
    if (key == static_cast<int>(Qt::Key_Meta) ||
	key == static_cast<int>(Qt::Key_Super_L) ||
	key == static_cast<int>(Qt::Key_Super_R))
	return BOBOL_INPUT_KEY_META;
    return key;
}

static int
qgcanvasinput_button(Qt::MouseButton button)
{
    switch (button) {
	case Qt::LeftButton:
	    return 0;
	case Qt::MiddleButton:
	    return 1;
	case Qt::RightButton:
	    return 2;
	default:
	    return BOBOL_INPUT_ANY;
    }
}

static unsigned int
qgcanvasinput_buttons(Qt::MouseButtons buttons)
{
    unsigned int result = 0;
    if (buttons.testFlag(Qt::LeftButton))
	result |= 1u << 0;
    if (buttons.testFlag(Qt::MiddleButton))
	result |= 1u << 1;
    if (buttons.testFlag(Qt::RightButton))
	result |= 1u << 2;
    return result;
}

static void
qgcanvasinput_position(const QMouseEvent *event, int &x, int &y)
{
    if (!event) {
	x = 0;
	y = 0;
	return;
    }
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    x = event->x();
    y = event->y();
#else
    x = static_cast<int>(event->position().x());
    y = static_cast<int>(event->position().y());
#endif
}

QgCanvasInput::QgCanvasInput() :
    m(new Impl)
{
    m->context.setProfile(&BObolInputContext::defaultViewProfile());
    m->context.setActionHandler(QgCanvasInput::actionDispatch, this);
}

QgCanvasInput::~QgCanvasInput()
{
	this->setEndpoint(NULL);
    delete m;
}

void
QgCanvasInput::setEndpoint(bobol_display_endpoint_t *endpoint)
{
    if (!m || m->endpoint == endpoint)
	return;
    if (m->endpoint)
	(void)bobol_display_endpoint_input_action_handler_clear_if(
	    m->endpoint, QgCanvasInput::actionDispatch, this);
    m->endpoint = endpoint;
    if (!m->endpoint)
	return;
    (void)bobol_display_endpoint_input_profile_set(m->endpoint,
	bobol_input_default_view_profile());
    (void)bobol_display_endpoint_input_action_handler_set(m->endpoint,
	QgCanvasInput::actionDispatch, this);
}

int
QgCanvasInput::actionDispatch(void *userData, BObolInputAction action,
	const BObolInputEvent *event)
{
    QgCanvasInput *input = static_cast<QgCanvasInput *>(userData);
    return input ? input->applyAction(action, event) : -1;
}

int
QgCanvasInput::applyAction(BObolInputAction action,
	const BObolInputEvent *event)
{
    if (!m || !m->dispatch_view || !event)
	return 0;
    struct bv_context *view_ctx = m->dispatch_view;
    struct bv *view = bv_context_view(view_ctx);
    const struct bv *cview = bv_context_view_const(view_ctx);

    switch (action) {
	case BOBOL_ACTION_TOGGLE_ADC:
	case BOBOL_ACTION_TOGGLE_MODEL_AXES:
	case BOBOL_ACTION_TOGGLE_VIEW_AXES:
	    return bobol_display_endpoint_input_faceplate_toggle_apply(
		m->endpoint, view_ctx, action, NULL);
	case BOBOL_ACTION_VIEW_2:
	case BOBOL_ACTION_VIEW_3:
	case BOBOL_ACTION_VIEW_4:
	case BOBOL_ACTION_VIEW_5:
	case BOBOL_ACTION_VIEW_6:
	case BOBOL_ACTION_VIEW_7:
	case BOBOL_ACTION_VIEW_FRONT:
	case BOBOL_ACTION_VIEW_TOP:
	case BOBOL_ACTION_VIEW_BOTTOM:
	case BOBOL_ACTION_VIEW_LEFT:
	case BOBOL_ACTION_VIEW_REAR:
	case BOBOL_ACTION_VIEW_RIGHT:
	    return bobol_input_view_orientation_apply(view_ctx, action);
	case BOBOL_ACTION_VIEW_PRIMARY_RELEASE:
	case BOBOL_ACTION_VIEW_SECONDARY_RELEASE:
	case BOBOL_ACTION_VIEW_CENTER_RELEASE: {
	    if (event->buttons ||
	std::fabs(static_cast<double>(event->x) - m->press_x) > 10.0 ||
	std::fabs(static_cast<double>(event->y) - m->press_y) > 10.0)
		return 0;
	    int dx = 1;
	    int dy = 1;
	    unsigned long long viewFlags = BV_ADJUST_IDLE;
	    if (action == BOBOL_ACTION_VIEW_PRIMARY_RELEASE) {
		if (m->mouse_mode == BV_ADJUST_CENTER) {
		    viewFlags = BV_ADJUST_CENTER;
		    dx = event->x;
		    dy = event->y;
		} else {
		    viewFlags = BV_ADJUST_SCALE;
		    dx = 10;
		    dy = 5;
		}
	    } else if (action == BOBOL_ACTION_VIEW_SECONDARY_RELEASE) {
		if (m->mouse_mode == BV_ADJUST_CENTER)
		    return 0;
		viewFlags = BV_ADJUST_SCALE;
		dx = 1;
		dy = 2;
	    } else {
		viewFlags = BV_ADJUST_CENTER;
		dx = event->x;
		dy = event->y;
	    }
	    point_t keypt = VINIT_ZERO;
	    return bv_adjust(view, dx, dy, keypt, 0,
		viewFlags);
	}
	case BOBOL_ACTION_VIEW_ROTATE:
	case BOBOL_ACTION_VIEW_PAN:
	case BOBOL_ACTION_VIEW_ZOOM:
	case BOBOL_ACTION_VIEW_ADJUST: {
	    if (event->type == BOBOL_INPUT_WHEEL) {
		const int dx = 100 + event->wheelDelta;
		point_t origin = VINIT_ZERO;
		return bv_adjust(view, dx, 100, origin, 0,
		    BV_ADJUST_SCALE);
	    }
	    unsigned long long viewFlags = BV_ADJUST_SCALE;
	    if (action == BOBOL_ACTION_VIEW_ROTATE)
		viewFlags = BV_ADJUST_ROT;
	    else if (action == BOBOL_ACTION_VIEW_PAN)
		viewFlags = BV_ADJUST_TRANS;
	    else if (action == BOBOL_ACTION_VIEW_ADJUST) {
		viewFlags = m->mouse_mode;
		if (viewFlags == BV_ADJUST_CENTER)
		    viewFlags = BV_ADJUST_SCALE;
	    }
	    point_t center;
	    mat_t viewCenter;
	    bv_center_mat_get(viewCenter, cview);
	    MAT_DELTAS_GET_NEG(center, viewCenter);
	    return bv_mouse_delta_adjust(view, event->dx, event->dy, center,
		viewFlags);
	}
	default:
	    return 0;
    }
}

int
QgCanvasInput::keyPressEvent(struct bv_context *view_ctx, int UNUSED(x_prev),
	int UNUSED(y_prev), QKeyEvent *event)
{
    QTCAD_EVENT("keyPress", 1);
	if (!view_ctx || !event)
	return 0;
    BObolInputEvent input;
    input.type = BOBOL_INPUT_KEY_PRESS;
    input.key = qgcanvasinput_key(event->key());
    input.modifiers = qgcanvasinput_modifiers(event->modifiers());
    m->dispatch_view = view_ctx;
    const int ret = m->endpoint ? bobol_display_endpoint_input_dispatch(
	m->endpoint, &input) : m->context.dispatch(&input);
    m->dispatch_view = NULL;
    return ret;
}

int
QgCanvasInput::mousePressEvent(struct bv_context *view_ctx, int UNUSED(x_prev),
	int UNUSED(y_prev), QMouseEvent *event)
{
    QTCAD_EVENT("mousePress", 1);
	if (!view_ctx || !event)
	return 0;
    BObolInputEvent input;
    input.type = BOBOL_INPUT_POINTER_PRESS;
    qgcanvasinput_position(event, input.x, input.y);
    input.button = qgcanvasinput_button(event->button());
    input.buttons = qgcanvasinput_buttons(event->buttons());
    input.modifiers = qgcanvasinput_modifiers(event->modifiers());
    m->dispatch_view = view_ctx;
    const int ret = m->endpoint ? bobol_display_endpoint_input_dispatch(
	m->endpoint, &input) : m->context.dispatch(&input);
    m->dispatch_view = NULL;
    return ret;
}

int
QgCanvasInput::mouseReleaseEvent(struct bv_context *view_ctx, double x_press,
	double y_press, int UNUSED(x_prev), int UNUSED(y_prev), QMouseEvent *event,
	int mode)
{
    QTCAD_EVENT("mouseRelease", 1);
	if (!view_ctx || !event)
	return 0;
    BObolInputEvent input;
    input.type = BOBOL_INPUT_POINTER_RELEASE;
    qgcanvasinput_position(event, input.x, input.y);
    input.button = qgcanvasinput_button(event->button());
    input.buttons = qgcanvasinput_buttons(event->buttons());
    input.modifiers = qgcanvasinput_modifiers(event->modifiers());
    m->dispatch_view = view_ctx;
    m->press_x = x_press;
    m->press_y = y_press;
    m->mouse_mode = mode;
    const int ret = m->endpoint ? bobol_display_endpoint_input_dispatch(
	m->endpoint, &input) : m->context.dispatch(&input);
    m->dispatch_view = NULL;
    return ret;
}

int
QgCanvasInput::mouseMoveEvent(struct bv_context *view_ctx, int x_prev, int y_prev,
	QMouseEvent *event, int mode)
{
    QTCAD_EVENT("mouseMove", 2);
    if (!view_ctx || !event)
	return 0;

    /* An application-layer gesture must see its first drag event.  View
     * navigation can still use a zero delta until the canvas has established
     * a baseline. */
    const bool first_motion = x_prev == -INT_MAX || y_prev == -INT_MAX;
    BObolInputEvent input;
    input.type = BOBOL_INPUT_POINTER_MOTION;
    qgcanvasinput_position(event, input.x, input.y);
    input.dx = first_motion ? 0 : input.x - x_prev;
    input.dy = first_motion ? 0 : input.y - y_prev;
    input.button = event->buttons().testFlag(Qt::LeftButton) ? 0 :
	BOBOL_INPUT_ANY;
    input.buttons = qgcanvasinput_buttons(event->buttons());
    input.modifiers = qgcanvasinput_modifiers(event->modifiers());
    m->dispatch_view = view_ctx;
    m->mouse_mode = mode;
    const int ret = m->endpoint ? bobol_display_endpoint_input_dispatch(
	m->endpoint, &input) : m->context.dispatch(&input);
    m->dispatch_view = NULL;
    return ret;
}

int
QgCanvasInput::wheelEvent(struct bv_context *view_ctx, QWheelEvent *event)
{
    QTCAD_EVENT("mouseWheel", 1);
	if (!view_ctx || !event)
	return 0;
    BObolInputEvent input;
    input.type = BOBOL_INPUT_WHEEL;
    input.wheelDelta = -event->angleDelta().y() / 8;
    input.modifiers = qgcanvasinput_modifiers(event->modifiers());
    m->dispatch_view = view_ctx;
    const int ret = m->endpoint ? bobol_display_endpoint_input_dispatch(
	m->endpoint, &input) : m->context.dispatch(&input);
    m->dispatch_view = NULL;
    return ret;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
