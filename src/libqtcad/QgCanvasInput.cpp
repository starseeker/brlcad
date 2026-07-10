/*                   Q G C A N V A S I N P U T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2021-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file QgCanvasInput.cpp
 *
 * CAD-specific mouse/keyboard bindings, held as per-canvas-instance
 * state.
 */

#include "common.h"

#include <QtGlobal>
#include <chrono>
#include <unordered_map>

extern "C" {
#include "bn/str.h"
}

#include "qtcad/defines.h"
#include "QgCanvasInput.h"
#include "QgLegacyViewContext.h"
#include "bv.h"

struct QgCanvasInput::Impl {
	std::unordered_map<qg_legacy_view *, long long> drag_update_ts;
};

QgCanvasInput::QgCanvasInput() :
	m(new Impl)
{
}

QgCanvasInput::~QgCanvasInput()
{
	delete m;
}

static int
qgcanvasinput_set_aet(qg_legacy_view *v, const char *aet_string)
{
	vect_t aet_vec;
	bn_decode_vect(aet_vec, aet_string);
	(void)bv_aet_set(qg_legacy_view_bv(v), aet_vec);
	(void)bv_context_update(qg_legacy_view_context(v), BV_CONTEXT_CHANGED_VIEW);
	return 1;
}

// TODO - look into QShortcut, see if it might be a better way
// to manage this
int
QgCanvasInput::keyPressEvent(qg_legacy_view *v, int UNUSED(x_prev),
                             int UNUSED(y_prev), QKeyEvent *k)
{
	QTCAD_EVENT("keyPress", 1);
	if (!v)
		return 0;
#if 0
	QString kstr = QKeySequence(k->key()).toString();
	bu_log("%s\n", kstr.toStdString().c_str());
#endif
	switch (k->key()) {
		case 'A': {
			struct bv_adc_state adc;
			if (!bv_adc_state_get(&adc, qg_legacy_view_bv_const(v)))
				return 0;
			adc.draw = !adc.draw;
			bv_adc_state_set(qg_legacy_view_bv(v), &adc);
			return 1;
		}
		case 'M': {
			struct bv_axes_state axes;
			if (!bv_model_axes_state_get(&axes, qg_legacy_view_bv_const(v)))
				return 0;
			axes.draw = !axes.draw;
			bv_model_axes_state_set(qg_legacy_view_bv(v), &axes);
			return 1;
		}
		case 'V': {
			struct bv_axes_state axes;
			if (!bv_view_axes_state_get(&axes, qg_legacy_view_bv_const(v)))
				return 0;
			axes.draw = !axes.draw;
			bv_view_axes_state_set(qg_legacy_view_bv(v), &axes);
			return 1;
		}
		case '2': {
				return qgcanvasinput_set_aet(v, "35 -25 0");
			}
		case '3': {
				return qgcanvasinput_set_aet(v, "35 25 0");
			}
		case '4': {
				return qgcanvasinput_set_aet(v, "45 45 0");
			}
		case '5': {
				return qgcanvasinput_set_aet(v, "145 25 0");
			}
		case '6': {
				return qgcanvasinput_set_aet(v, "215 25 0");
			}
		case '7': {
				return qgcanvasinput_set_aet(v, "325 25 0");
			}
		case 'F': {
				return qgcanvasinput_set_aet(v, "0 0 0");
			}
		case 'T': {
				return qgcanvasinput_set_aet(v, "270 90 0");
			}
		case 'B': {
				return qgcanvasinput_set_aet(v, "270 -90 0");
			}
		case 'L': {
				return qgcanvasinput_set_aet(v, "90 0 0");
			}
		case 'R': {
				if (k->modifiers().testFlag(Qt::ShiftModifier) == true) {
					return qgcanvasinput_set_aet(v, "180 0 0");
				}
				return qgcanvasinput_set_aet(v, "270 0 0");
			}
		default:
			break;
	}
	return 0;
}

int
QgCanvasInput::mousePressEvent(qg_legacy_view *v, int UNUSED(x_prev),
                               int UNUSED(y_prev), QMouseEvent *e)
{
	QTCAD_EVENT("mousePress", 1);

	if (!v)
		return 0;

	// If we're intending the mouse motion to do the work,
	// then the press has to be a no-op.  If we're going
	// to do configurable key bindings, this will take some
	// thought - if we want unmodded left button to be a
	// rotation, and Ctrl+Left to do something else, these
	// checks are all going to have to be exact-flag-combo-only
	// actions.
	if (e->modifiers()) {
		return 0;
	}

	if (e->buttons().testFlag(Qt::LeftButton)) {
		//bu_log("Press Left\n");
	}
	if (e->buttons().testFlag(Qt::RightButton)) {
		//bu_log("Press Right\n");
	}
	if (e->buttons().testFlag(Qt::MiddleButton)) {
		//bu_log("Press Middle\n");
	}

	return 0;
}

int
QgCanvasInput::mouseReleaseEvent(qg_legacy_view *v, double x_press,
                                 double y_press, int UNUSED(x_prev),
                                 int UNUSED(y_prev), QMouseEvent *e, int mode)
{
	QTCAD_EVENT("mouseRelease", 1);

	if (!v)
		return 0;

	m->drag_update_ts.erase(v);

	// If we're intending the mouse motion to do the work,
	// then the release has to be a no-op.  If we're going
	// to do configurable key bindings, this will take some
	// thought - if we want unmodded left button to be a
	// rotation, and Ctrl+Left to do something else, these
	// checks are all going to have to be exact-flag-combo-only
	// actions.
	if (e->modifiers()) {
		return 0;
	}
	if (e->buttons().testFlag(Qt::LeftButton) || e->buttons().testFlag(Qt::RightButton)) {
		return 0;
	}

	double cx, cy;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	cx = (double)e->x();
	cy = (double)e->y();
#else
	cx = e->position().x();
	cy = e->position().y();
#endif
	if ((fabs(cx - x_press) > 10) || (fabs(cy - y_press) > 10))
		return 0;

	int dx = 1;
	int dy = 1;
	unsigned long long view_flags = BV_ADJUST_IDLE;

	if (e->button() == Qt::LeftButton) {
		//bu_log("Release Left\n");
		if (mode != BV_ADJUST_CENTER) {
			view_flags = BV_ADJUST_SCALE;
			dx = 10;
			dy = 5;
		}
		else {
			view_flags = BV_ADJUST_CENTER;
			dx = (int)cx;
			dy = (int)cy;
		}
	}
	if (e->button() == Qt::RightButton) {
		//bu_log("Release Right\n");
		if (mode == BV_ADJUST_CENTER)
			return 0;
		view_flags = BV_ADJUST_SCALE;
		dx = 1;
		dy = 2;
	}

	if (e->button() == Qt::MiddleButton) {
		//bu_log("Release Center\n");
		view_flags = BV_ADJUST_CENTER;
		dx = (int)cx;
		dy = (int)cy;
	}

	point_t keypt = VINIT_ZERO;
	return bv_adjust(qg_legacy_view_bv(v), dx, dy, keypt, 0, view_flags);
}

int
QgCanvasInput::mouseMoveEvent(qg_legacy_view *v, int x_prev, int y_prev,
                              QMouseEvent *e, int mode)
{
	QTCAD_EVENT("mouseMove", 2);

	if (!v)
		return 0;

	unsigned long long view_flags = BV_ADJUST_IDLE;

	if (x_prev == -INT_MAX) {
		//x_prev = e->x();
		//y_prev = e->y();
		return 0;
	}

	view_flags = mode;
	if (mode == BV_ADJUST_CENTER)
		view_flags = BV_ADJUST_SCALE;

	if (e->buttons().testFlag(Qt::LeftButton)) {
		//bu_log("Left\n");

		if (e->modifiers().testFlag(Qt::ControlModifier)) {
			//bu_log("Ctrl+Left\n");
			view_flags = BV_ADJUST_ROT;
		}

		if (e->modifiers().testFlag(Qt::ShiftModifier)) {
			//bu_log("Shift+Left\n");
			view_flags = BV_ADJUST_TRANS;
		}

		if (e->modifiers().testFlag(Qt::ShiftModifier) && e->modifiers().testFlag(Qt::ControlModifier)) {
			//bu_log("Ctrl+Shift+Left\n");
			view_flags = BV_ADJUST_SCALE;
		}
	}

	if (e->buttons().testFlag(Qt::MiddleButton)) {
		//bu_log("Middle\n");
	}

	if (e->buttons().testFlag(Qt::RightButton)) {
		//bu_log("Right\n");
	}

	if (!e->buttons().testFlag(Qt::LeftButton))
		return 0;

	long long now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
	                           std::chrono::steady_clock::now().time_since_epoch()).count();
	auto ts_it = m->drag_update_ts.find(v);
	if (ts_it != m->drag_update_ts.end() && (now_ms - ts_it->second) < s_drag_update_interval_ms)
		return -1;
	m->drag_update_ts[v] = now_ms;

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
	int dx = e->x() - x_prev;
	int dy = e->y() - y_prev;
#else
	int dx = e->position().x() - x_prev;
	int dy = e->position().y() - y_prev;
#endif

	if (view_flags == BV_ADJUST_SCALE) {
		// Build in some sensitivity to how much the mouse moved when doing
		// a motion based scale
		int mdelta = (abs(dx) > abs(dy)) ? dx : -dy;
		int f = (int)(2*100*(double)abs(mdelta) /
			(double)bv_height_get(qg_legacy_view_bv_const(v)));

		if (mdelta > 0) {
			dy = 101 + f;
			dx = 100;
		}
		else {
			dy = 99 - f;
			dx = 100;
		}
	}


	// TODO - the key point and the mode/flags are all hardcoded
	// right now, but eventually for shift grips they will need to
	// respond to the various mod keys.  The intent is to set flags
	// based on which mod keys are set to allow bv_adjust to
	// do the correct math.
	point_t center;
	mat_t view_center;
	bv_center_mat_get(view_center, qg_legacy_view_bv_const(v));
	MAT_DELTAS_GET_NEG(center, view_center);

	return bv_adjust(qg_legacy_view_bv(v), dx, dy, center, 0, view_flags);
}

int
QgCanvasInput::wheelEvent(qg_legacy_view *v, QWheelEvent *e)
{
	QTCAD_EVENT("mouseWheel", 1);

	if (!v)
		return 0;

	QPoint delta = e->angleDelta();
	int mdelta = -1 * delta.y() / 8;

	int dx = 100 + mdelta;
	int dy = 100;

	point_t origin = VINIT_ZERO;
	return bv_adjust(qg_legacy_view_bv(v), dx, dy, origin, 0, BV_ADJUST_SCALE);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
