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
 * Private libqtcad declarations for the retained display-manager fallback
 * bridge.  These helpers intentionally stay out of the installed qtcad
 * legacy-view API while qged/qtcad still carries the temporary DM path.
 */

#ifndef QGLEGACYVIEWDM_H
#define QGLEGACYVIEWDM_H

#include "qtcad/QgLegacyView.h"

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

QTCAD_EXPORT extern fastf_t *qg_legacy_view_scale_storage_get(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_display_manager_set(qg_legacy_view *view,
	qg_legacy_dm *dmp);

QTCAD_EXPORT extern qg_legacy_dm *qg_legacy_view_display_manager_get(
	qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_dm_bind(qg_legacy_view *view,
	qg_legacy_dm *dmp);

QTCAD_EXPORT extern int qg_legacy_view_dm_sync_dimensions(qg_legacy_view *view,
	qg_legacy_dm *dmp);

QTCAD_EXPORT extern unsigned long long qg_legacy_view_dm_hash(
	qg_legacy_dm *dmp);

QTCAD_EXPORT extern int qg_legacy_view_dm_native_repaint_pending_get(
	qg_legacy_dm *dmp);

QTCAD_EXPORT extern void qg_legacy_view_dm_native_repaint_pending_set(
	qg_legacy_dm *dmp, int pending);

QTCAD_EXPORT extern int qg_legacy_view_dm_configure_window(qg_legacy_dm *dmp,
	int force);

QTCAD_EXPORT extern int qg_legacy_view_dm_width_get(qg_legacy_dm *dmp);

QTCAD_EXPORT extern int qg_legacy_view_dm_height_get(qg_legacy_dm *dmp);

QTCAD_EXPORT extern void qg_legacy_view_dm_dimensions_set(qg_legacy_dm *dmp,
	int width, int height);

QTCAD_EXPORT extern int qg_legacy_view_dm_framebuffer_setup_existing(
	qg_legacy_fb *ifp, qg_legacy_dm *dmp);

QTCAD_EXPORT extern qg_legacy_dm *qg_legacy_view_dm_open_qtgl(void *context);

QTCAD_EXPORT extern qg_legacy_dm *qg_legacy_view_dm_open_swrast(
	qg_legacy_view *view, void *context);

QTCAD_EXPORT extern int qg_legacy_view_dm_close(qg_legacy_dm *dmp);

QTCAD_EXPORT extern int qg_legacy_view_dm_setup_qtgl(qg_legacy_dm *dmp,
	fastf_t zmin, fastf_t zmax);

QTCAD_EXPORT extern int qg_legacy_view_dm_setup_swrast(qg_legacy_dm *dmp,
	fastf_t zmin, fastf_t zmax);

QTCAD_EXPORT extern int qg_legacy_view_dm_background_get(qg_legacy_dm *dmp,
	unsigned char bg1[3],
	unsigned char bg2[3]);

QTCAD_EXPORT extern int qg_legacy_view_dm_background_set(qg_legacy_dm *dmp,
	const unsigned char bg1[3],
	const unsigned char bg2[3]);

QTCAD_EXPORT extern int qg_legacy_view_dm_background_restore(qg_legacy_dm *dmp);

QTCAD_EXPORT extern int qg_legacy_view_dm_draw_begin(qg_legacy_dm *dmp);

QTCAD_EXPORT extern int qg_legacy_view_dm_draw_end(qg_legacy_dm *dmp);

QTCAD_EXPORT extern int qg_legacy_view_dm_display_image_get(qg_legacy_dm *dmp,
	unsigned char **image,
	int copy,
	int release);

QTCAD_EXPORT extern void qg_legacy_view_dm_display_image_release(
	unsigned char *image);

QTCAD_EXPORT extern int qg_legacy_view_dm_load_current_model2view(
	qg_legacy_dm *dmp, const qg_legacy_view *view, int which_eye_or_mode);

QTCAD_EXPORT extern qg_legacy_fb *qg_legacy_view_framebuffer_raw_create(
	const char *type);

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

QTCAD_EXPORT extern int qg_legacy_view_framebuffer_configure(qg_legacy_fb *ifp,
	int width,
	int height);

QTCAD_EXPORT extern void *qg_legacy_view_dm_open_context(qg_legacy_view *view);

QTCAD_EXPORT extern int qg_legacy_view_dm_draw(qg_legacy_view *view);

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
