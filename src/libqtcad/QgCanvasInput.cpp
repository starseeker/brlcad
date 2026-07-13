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

#include <chrono>
#include <cmath>
#include <unordered_map>

extern "C" {
#include "bn/str.h"
}

#include "qtcad/defines.h"
#include "QgCanvasInput.h"
#include "QgLegacyViewContext.h"
#include "brlobol/display_endpoint.h"
#include "bv.h"

struct QgCanvasInput::Impl {
    std::unordered_map<qg_legacy_view *, long long> drag_update_ts;
    BRLObolInputContext context;
    brlobol_display_endpoint_t *endpoint = NULL;
    qg_legacy_view *dispatch_view = NULL;
    double press_x = 0.0;
    double press_y = 0.0;
    int mouse_mode = BV_ADJUST_SCALE;
};

static int
qgcanvasinput_set_aet(qg_legacy_view *v, const char *aet_string)
{
    if (!v || !aet_string)
	return 0;
    vect_t aet_vec;
    bn_decode_vect(aet_vec, aet_string);
    (void)bv_aet_set(qg_legacy_view_bv(v), aet_vec);
    (void)bv_context_update(qg_legacy_view_context(v), BV_CONTEXT_CHANGED_VIEW);
    return 1;
}

static unsigned int
qgcanvasinput_modifiers(Qt::KeyboardModifiers modifiers)
{
    unsigned int result = BRLOBOL_INPUT_MOD_NONE;
    if (modifiers.testFlag(Qt::ShiftModifier))
	result |= BRLOBOL_INPUT_MOD_SHIFT;
    if (modifiers.testFlag(Qt::ControlModifier))
	result |= BRLOBOL_INPUT_MOD_CONTROL;
    if (modifiers.testFlag(Qt::AltModifier))
	result |= BRLOBOL_INPUT_MOD_ALT;
    if (modifiers.testFlag(Qt::MetaModifier))
	result |= BRLOBOL_INPUT_MOD_META;
    return result;
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
	    return BRLOBOL_INPUT_ANY;
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
    m->context.setProfile(&BRLObolInputContext::defaultViewProfile());
    m->context.setActionHandler(QgCanvasInput::actionDispatch, this);
}

QgCanvasInput::~QgCanvasInput()
{
	this->setEndpoint(NULL);
    delete m;
}

void
QgCanvasInput::setEndpoint(brlobol_display_endpoint_t *endpoint)
{
    if (!m || m->endpoint == endpoint)
	return;
    if (m->endpoint)
	(void)brlobol_display_endpoint_input_action_handler_clear_if(
	    m->endpoint, QgCanvasInput::actionDispatch, this);
    m->endpoint = endpoint;
    if (!m->endpoint)
	return;
    (void)brlobol_display_endpoint_input_profile_set(m->endpoint,
	brlobol_input_default_view_profile());
    (void)brlobol_display_endpoint_input_action_handler_set(m->endpoint,
	QgCanvasInput::actionDispatch, this);
}

int
QgCanvasInput::actionDispatch(void *userData, BRLObolInputAction action,
	const BRLObolInputEvent *event)
{
    QgCanvasInput *input = static_cast<QgCanvasInput *>(userData);
    return input ? input->applyAction(action, event) : -1;
}

int
QgCanvasInput::applyAction(BRLObolInputAction action,
	const BRLObolInputEvent *event)
{
    if (!m || !m->dispatch_view || !event)
	return 0;
    qg_legacy_view *view = m->dispatch_view;

    switch (action) {
	case BRLOBOL_ACTION_TOGGLE_ADC: {
	    struct bv_adc_state adc;
	    if (!bv_adc_state_get(&adc, qg_legacy_view_bv_const(view)))
		return 0;
	    adc.draw = !adc.draw;
	    bv_adc_state_set(qg_legacy_view_bv(view), &adc);
	    return 1;
	}
	case BRLOBOL_ACTION_TOGGLE_MODEL_AXES: {
	    struct bv_axes_state axes;
	    if (!bv_model_axes_state_get(&axes,
		qg_legacy_view_bv_const(view)))
		return 0;
	    axes.draw = !axes.draw;
	    bv_model_axes_state_set(qg_legacy_view_bv(view), &axes);
	    return 1;
	}
	case BRLOBOL_ACTION_TOGGLE_VIEW_AXES: {
	    struct bv_axes_state axes;
	    if (!bv_view_axes_state_get(&axes,
		qg_legacy_view_bv_const(view)))
		return 0;
	    axes.draw = !axes.draw;
	    bv_view_axes_state_set(qg_legacy_view_bv(view), &axes);
	    return 1;
	}
	case BRLOBOL_ACTION_VIEW_2:
	    return qgcanvasinput_set_aet(view, "35 -25 0");
	case BRLOBOL_ACTION_VIEW_3:
	    return qgcanvasinput_set_aet(view, "35 25 0");
	case BRLOBOL_ACTION_VIEW_4:
	    return qgcanvasinput_set_aet(view, "45 45 0");
	case BRLOBOL_ACTION_VIEW_5:
	    return qgcanvasinput_set_aet(view, "145 25 0");
	case BRLOBOL_ACTION_VIEW_6:
	    return qgcanvasinput_set_aet(view, "215 25 0");
	case BRLOBOL_ACTION_VIEW_7:
	    return qgcanvasinput_set_aet(view, "325 25 0");
	case BRLOBOL_ACTION_VIEW_FRONT:
	    return qgcanvasinput_set_aet(view, "0 0 0");
	case BRLOBOL_ACTION_VIEW_TOP:
	    return qgcanvasinput_set_aet(view, "270 90 0");
	case BRLOBOL_ACTION_VIEW_BOTTOM:
	    return qgcanvasinput_set_aet(view, "270 -90 0");
	case BRLOBOL_ACTION_VIEW_LEFT:
	    return qgcanvasinput_set_aet(view, "90 0 0");
	case BRLOBOL_ACTION_VIEW_REAR:
	    return qgcanvasinput_set_aet(view, "180 0 0");
	case BRLOBOL_ACTION_VIEW_RIGHT:
	    return qgcanvasinput_set_aet(view, "270 0 0");
	case BRLOBOL_ACTION_VIEW_PRIMARY_RELEASE:
	case BRLOBOL_ACTION_VIEW_SECONDARY_RELEASE:
	case BRLOBOL_ACTION_VIEW_CENTER_RELEASE: {
	    if (event->buttons ||
	std::fabs(static_cast<double>(event->x) - m->press_x) > 10.0 ||
	std::fabs(static_cast<double>(event->y) - m->press_y) > 10.0)
		return 0;
	    int dx = 1;
	    int dy = 1;
	    unsigned long long viewFlags = BV_ADJUST_IDLE;
	    if (action == BRLOBOL_ACTION_VIEW_PRIMARY_RELEASE) {
		if (m->mouse_mode == BV_ADJUST_CENTER) {
		    viewFlags = BV_ADJUST_CENTER;
		    dx = event->x;
		    dy = event->y;
		} else {
		    viewFlags = BV_ADJUST_SCALE;
		    dx = 10;
		    dy = 5;
		}
	    } else if (action == BRLOBOL_ACTION_VIEW_SECONDARY_RELEASE) {
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
	    return bv_adjust(qg_legacy_view_bv(view), dx, dy, keypt, 0,
		viewFlags);
	}
	case BRLOBOL_ACTION_VIEW_ROTATE:
	case BRLOBOL_ACTION_VIEW_PAN:
	case BRLOBOL_ACTION_VIEW_ZOOM:
	case BRLOBOL_ACTION_VIEW_ADJUST: {
	    if (event->type == BRLOBOL_INPUT_WHEEL) {
		const int dx = 100 + event->wheelDelta;
		point_t origin = VINIT_ZERO;
		return bv_adjust(qg_legacy_view_bv(view), dx, 100, origin, 0,
		    BV_ADJUST_SCALE);
	    }
	    unsigned long long viewFlags = BV_ADJUST_SCALE;
	    if (action == BRLOBOL_ACTION_VIEW_ROTATE)
		viewFlags = BV_ADJUST_ROT;
	    else if (action == BRLOBOL_ACTION_VIEW_PAN)
		viewFlags = BV_ADJUST_TRANS;
	    else if (action == BRLOBOL_ACTION_VIEW_ADJUST) {
		viewFlags = m->mouse_mode;
		if (viewFlags == BV_ADJUST_CENTER)
		    viewFlags = BV_ADJUST_SCALE;
	    }
	    int dx = event->dx;
	    int dy = event->dy;
	    if (viewFlags == BV_ADJUST_SCALE) {
		const int mdelta = std::abs(dx) > std::abs(dy) ? dx : -dy;
		const int height = std::max(1, bv_height_get(
		    qg_legacy_view_bv_const(view)));
		const int factor = static_cast<int>(200.0 * std::abs(mdelta) /
		    static_cast<double>(height));
		dy = mdelta > 0 ? 101 + factor : 99 - factor;
		dx = 100;
	    }
	    point_t center;
	    mat_t viewCenter;
	    bv_center_mat_get(viewCenter, qg_legacy_view_bv_const(view));
	    MAT_DELTAS_GET_NEG(center, viewCenter);
	    return bv_adjust(qg_legacy_view_bv(view), dx, dy, center, 0,
		viewFlags);
	}
	default:
	    return 0;
    }
}

int
QgCanvasInput::keyPressEvent(qg_legacy_view *view, int UNUSED(x_prev),
	int UNUSED(y_prev), QKeyEvent *event)
{
    QTCAD_EVENT("keyPress", 1);
    if (!view || !event)
	return 0;
    BRLObolInputEvent input;
    input.type = BRLOBOL_INPUT_KEY_PRESS;
    input.key = event->key();
    input.modifiers = qgcanvasinput_modifiers(event->modifiers());
    m->dispatch_view = view;
    const int ret = m->endpoint ? brlobol_display_endpoint_input_dispatch(
	m->endpoint, &input) : m->context.dispatch(&input);
    m->dispatch_view = NULL;
    return ret;
}

int
QgCanvasInput::mousePressEvent(qg_legacy_view *view, int UNUSED(x_prev),
	int UNUSED(y_prev), QMouseEvent *event)
{
    QTCAD_EVENT("mousePress", 1);
    if (!view || !event)
	return 0;
    BRLObolInputEvent input;
    input.type = BRLOBOL_INPUT_POINTER_PRESS;
    qgcanvasinput_position(event, input.x, input.y);
    input.button = qgcanvasinput_button(event->button());
    input.buttons = qgcanvasinput_buttons(event->buttons());
    input.modifiers = qgcanvasinput_modifiers(event->modifiers());
    m->dispatch_view = view;
    const int ret = m->endpoint ? brlobol_display_endpoint_input_dispatch(
	m->endpoint, &input) : m->context.dispatch(&input);
    m->dispatch_view = NULL;
    return ret;
}

int
QgCanvasInput::mouseReleaseEvent(qg_legacy_view *view, double x_press,
	double y_press, int UNUSED(x_prev), int UNUSED(y_prev), QMouseEvent *event,
	int mode)
{
    QTCAD_EVENT("mouseRelease", 1);
    if (!view || !event)
	return 0;
    m->drag_update_ts.erase(view);

    BRLObolInputEvent input;
    input.type = BRLOBOL_INPUT_POINTER_RELEASE;
    qgcanvasinput_position(event, input.x, input.y);
    input.button = qgcanvasinput_button(event->button());
    input.buttons = qgcanvasinput_buttons(event->buttons());
    input.modifiers = qgcanvasinput_modifiers(event->modifiers());
    m->dispatch_view = view;
    m->press_x = x_press;
    m->press_y = y_press;
    m->mouse_mode = mode;
    const int ret = m->endpoint ? brlobol_display_endpoint_input_dispatch(
	m->endpoint, &input) : m->context.dispatch(&input);
    m->dispatch_view = NULL;
    return ret;
}

int
QgCanvasInput::mouseMoveEvent(qg_legacy_view *view, int x_prev, int y_prev,
	QMouseEvent *event, int mode)
{
    QTCAD_EVENT("mouseMove", 2);
    if (!view || !event || x_prev == -INT_MAX ||
	!event->buttons().testFlag(Qt::LeftButton))
	return 0;

    const long long now = std::chrono::duration_cast<std::chrono::milliseconds>(
	std::chrono::steady_clock::now().time_since_epoch()).count();
    const auto previous = m->drag_update_ts.find(view);
    if (previous != m->drag_update_ts.end() && now - previous->second <
	s_drag_update_interval_ms)
	return -1;
    m->drag_update_ts[view] = now;

    BRLObolInputEvent input;
    input.type = BRLOBOL_INPUT_POINTER_MOTION;
    qgcanvasinput_position(event, input.x, input.y);
    input.dx = input.x - x_prev;
    input.dy = input.y - y_prev;
    input.button = 0;
    input.buttons = qgcanvasinput_buttons(event->buttons());
    input.modifiers = qgcanvasinput_modifiers(event->modifiers());
    m->dispatch_view = view;
    m->mouse_mode = mode;
    const int ret = m->endpoint ? brlobol_display_endpoint_input_dispatch(
	m->endpoint, &input) : m->context.dispatch(&input);
    m->dispatch_view = NULL;
    return ret;
}

int
QgCanvasInput::wheelEvent(qg_legacy_view *view, QWheelEvent *event)
{
    QTCAD_EVENT("mouseWheel", 1);
    if (!view || !event)
	return 0;
    BRLObolInputEvent input;
    input.type = BRLOBOL_INPUT_WHEEL;
    input.wheelDelta = -event->angleDelta().y() / 8;
    input.modifiers = qgcanvasinput_modifiers(event->modifiers());
    m->dispatch_view = view;
    const int ret = m->endpoint ? brlobol_display_endpoint_input_dispatch(
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
