/*                   Q G C A N V A S I N P U T . H
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
/** @file QgCanvasInput.h
 *
 * Per-canvas CAD-specific mouse/keyboard binding handler.
 *
 * This class replaces the deleted legacy free-function input API and its
 * global static drag-tracking maps.  Moving the state inside a class
 * instance eliminates cross-canvas interference when multiple canvases are
 * simultaneously active (Phase 3 of the libqtcad refactor plan).
 */

#ifndef QGCANVASINPUT_H
#define QGCANVASINPUT_H

#include "common.h"

#include <QKeyEvent>
#include <QMouseEvent>
#include <QWheelEvent>

#include "brlobol/input.h"
#include "qtcad/QgLegacyView.h"

/**
 * Per-canvas input-binding handler.
 *
 * Each QgGL / QgSW canvas embeds one instance of this class (via
 * QgCanvasState::input).  The drag-tracking maps that were previously
 * global statics in the deleted legacy binding layer are now per-instance,
 * so concurrent canvases cannot interfere with each other.
 */
class QgCanvasInput {
public:
	QgCanvasInput();
	~QgCanvasInput();
	void setEndpoint(struct brlobol_display_endpoint *endpoint);

	int keyPressEvent(qg_legacy_view *v, int x_prev, int y_prev,
	                  QKeyEvent *k);

	int mousePressEvent(qg_legacy_view *v, int x_prev, int y_prev,
	                    QMouseEvent *e);

	int mouseReleaseEvent(qg_legacy_view *v, double x_press, double y_press,
	                      int x_prev, int y_prev, QMouseEvent *e, int mode);

	int mouseMoveEvent(qg_legacy_view *v, int x_prev, int y_prev,
	                   QMouseEvent *e, int mode);

	int wheelEvent(qg_legacy_view *v, QWheelEvent *e);

private:
	struct Impl;
	static int actionDispatch(void *userData, BRLObolInputAction action,
		const BRLObolInputEvent *event);
	int applyAction(BRLObolInputAction action,
		const BRLObolInputEvent *event);

	Impl *m = nullptr;
	static const long long s_drag_update_interval_ms = 16;
};

#endif /* QGCANVASINPUT_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
