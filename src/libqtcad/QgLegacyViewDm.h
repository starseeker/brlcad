/*                  Q G L E G A C Y V I E W D M . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
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
/** @file QgLegacyViewDm.h
 *
 * Private libqtcad declarations for transitional framebuffer compatibility
 * services.  Ordinary qtcad canvases render via Obol; these helpers remain
 * only while qged/fbserv and libdm's standalone Qt framebuffer windows are
 * moved to the Obol window-host path.
 */

#ifndef QGLEGACYVIEWDM_H
#define QGLEGACYVIEWDM_H

#include "qtcad/QgLegacyView.h"
#include "vmath.h"

class QgGL;
class QgSW;
class QWidget;

struct qg_legacy_dm;
struct qg_legacy_fb;

class QTCAD_EXPORT QgCanvasBridgeFactory {
public:
	static QgGL *create_qtgl(QWidget *parent, qg_legacy_fb *fbp);
	static QgSW *create_swrast(QWidget *parent, qg_legacy_fb *fbp);
};

QTCAD_EXPORT extern qg_legacy_dm *qg_legacy_view_dm_open_qtgl(void *context);

QTCAD_EXPORT extern qg_legacy_dm *qg_legacy_view_dm_open_swrast(
	qg_legacy_view *view, void *context);

QTCAD_EXPORT extern int qg_legacy_view_dm_close(qg_legacy_dm *dmp);

QTCAD_EXPORT extern qg_legacy_fb *qg_legacy_view_framebuffer_handle_from_raw(
	void *ifp);

QTCAD_EXPORT extern QgGL *qg_legacy_view_framebuffer_qtgl_canvas_get(
	qg_legacy_fb *ifp);

QTCAD_EXPORT extern QgSW *qg_legacy_view_framebuffer_swrast_canvas_get(
	qg_legacy_fb *ifp);

QTCAD_EXPORT extern int qg_legacy_view_framebuffer_release(qg_legacy_fb *ifp,
	int initialized);

QTCAD_EXPORT extern int qg_legacy_view_framebuffer_standalone_get(
	qg_legacy_fb *ifp);

QTCAD_EXPORT extern const char *qg_legacy_view_dm_init_messages(void);

#endif /* QGLEGACYVIEWDM_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
