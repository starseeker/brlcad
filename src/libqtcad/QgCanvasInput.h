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
 * Each canvas owns its drag state, so simultaneous canvases cannot interfere
 * with each other.
 */

#ifndef QGCANVASINPUT_H
#define QGCANVASINPUT_H

#include "common.h"

#include <QKeyEvent>
#include <QMouseEvent>
#include <QWheelEvent>

#include "brlobol/input.h"
#include "bv.h"

/**
 * Per-canvas input-binding handler.
 *
 * Each QgGL / QgSW canvas embeds one instance through QgCanvasState::input.
 */
class QgCanvasInput {
public:
	QgCanvasInput();
	~QgCanvasInput();
	void setEndpoint(struct brlobol_display_endpoint *endpoint);

	int keyPressEvent(struct bv_context *view_ctx, int x_prev, int y_prev,
	                  QKeyEvent *k);

	int mousePressEvent(struct bv_context *view_ctx, int x_prev, int y_prev,
	                    QMouseEvent *e);

	int mouseReleaseEvent(struct bv_context *view_ctx, double x_press, double y_press,
	                      int x_prev, int y_prev, QMouseEvent *e, int mode);

	int mouseMoveEvent(struct bv_context *view_ctx, int x_prev, int y_prev,
	                   QMouseEvent *e, int mode);

	int wheelEvent(struct bv_context *view_ctx, QWheelEvent *e);

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
