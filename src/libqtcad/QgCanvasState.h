/*                   Q G C A N V A S S T A T E . H
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
/** @file QgCanvasState.h
 *
 * Private (libqtcad-internal) pimpl struct that consolidates the state
 * shared between QgGL and QgSW, together with inline helper operations
 * that eliminate textual duplication in those two canvas implementation
 * files (Phase 3, §2.4 of the refactor plan).
 *
 * This header is NOT part of the installed libqtcad public API.
 */

#ifndef QGCANVASSTATE_H
#define QGCANVASSTATE_H

#include "common.h"

#include <climits>
#include <cmath>
#include <cstring>
#include <QImage>
#include <QSize>
#include <QWidget>

#include "brlobol/init.h"
#include "brlobol/adc.h"
#include "brlobol/axes.h"
#include "brlobol/grid.h"
#include "brlobol/view_controller.h"
#include "QgObolContextManager.h"
#include "QgLegacyViewDm.h"
#include "QgLegacyViewContext.h"
#include "qtcad/QgLegacyView.h"
#include "rt/view.h"

#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/SoViewport.h>
#include <Inventor/SbColor.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoSeparator.h>

#include "QgCanvasInput.h"

/**
 * Plain-data struct that consolidates the private state shared between
 * QgGL and QgSW.  It is held as a pimpl-style raw pointer in both canvas
 * class headers so that implementation details are not part of the
 * installed public interface.
 *
 * Ownership summary
 * ─────────────────
 * local_v  – Opaque handle for the view created by
 *            qg_legacy_view_local_create in the canvas constructor; released
 *            by qg_legacy_view_local_free in the canvas destructor.  The
 *            canvas owns it unconditionally.
 *
 * v        – Normally points to local_v (the widget-owned view).  When
 *            QgCanvasBase::set_view() is called with a non-null external
 *            view, v points to that caller-owned qg_legacy_view instead.
 *            Passing nullptr to set_view() reverts v back to local_v.  The
 *            canvas never frees v — it only frees local_v.
 *
 * dmp      – Opened lazily in the first paint event through the opaque
 *            qg_legacy_view_dm_open_* helpers; closed by
 *            qg_legacy_view_dm_close() in the canvas destructor.  The canvas
 *            owns it.
 *
 * ifp      – When the canvas is constructed without a caller-supplied fb*,
 *            it allocates a raw framebuffer through
 *            qg_legacy_view_framebuffer_raw_create() and owns it.  The
 *            destructor releases it through qg_legacy_view_framebuffer_release(),
 *            which preserves the retained initialized-vs-raw lifetime split.
 *            When the caller supplies a standalone fb* at construction time,
 *            the canvas does not close it on destruction.
 */
struct QgCanvasState {
	/* ---- view / dm / fb plumbing ---- */
	qg_legacy_view *v = nullptr;           /* active view: normally == local_v,
	                                       set_view() can redirect to an
	                                       external caller-owned legacy view     */
	qg_legacy_dm   *dmp = nullptr;        /* retained display manager (owned) */
	qg_legacy_fb   *ifp = nullptr;        /* framebuffer (see ownership note) */
	qg_legacy_view *local_v = nullptr;     /* widget-owned view (canvas owns)   */
	BRLObolViewController *obol = nullptr; /* Obol-canonical view controller */

	/* ---- hash tracking for incremental updates ---- */
	unsigned long long prev_dhash = 0;
	unsigned long long prev_vhash = 0;

	/* ---- input-binding flags ---- */
	bool use_default_keybindings   = true;
	bool use_default_mousebindings = true;
	int  lmouse_mode = -1;  /* set to RT_VIEW_ADJUST_SCALE in canvas constructor */

	/* ---- widget-level tracking ---- */
	int    current = 1;     /* 1 = this view is active */
	bool   m_init = false;  /* DM has been opened and configured */
	int    x_prev = -INT_MAX;
	int    y_prev = -INT_MAX;
	double x_press_pos = -INT_MAX;
	double y_press_pos = -INT_MAX;
	bool   obol_paint_initialized = false;
	bool   fb_update_queued = false;

	/* ---- per-canvas input handler ---- */
	QgCanvasInput input;
};

/* ------------------------------------------------------------------ */
/* Shared inline helpers (static so each TU gets its own copy)        */
/* ------------------------------------------------------------------ */

/** Compute the physical render size of a widget (accounting for DPR). */
static inline QSize
qgcanvas_render_size(const QWidget *w)
{
	if (!w)
		return QSize();
	qreal dpr = w->devicePixelRatioF();
	return QSize(qMax(1, static_cast<int>(std::ceil(w->width()  * dpr))),
	             qMax(1, static_cast<int>(std::ceil(w->height() * dpr))));
}

/** Keep the Obol view controller viewport aligned with the Qt canvas. */
static inline void
qgcanvas_sync_obol_viewport(QgCanvasState &s, const QWidget *w)
{
	if (!s.obol)
		return;
	QSize rsize = qgcanvas_render_size(w);
	s.obol->setViewportSize(static_cast<unsigned int>(rsize.width()),
	    static_cast<unsigned int>(rsize.height()));
}

/** Application-wide qtcad Obol context manager. */
static inline SoDB::ContextManager *
qgcanvas_obol_context_manager(void)
{
	static QgObolContextManager manager;
	return &manager;
}

static inline void
qgcanvas_request_obol_render_if_idle(QgCanvasState &s, const char *reason)
{
	if (s.obol && !s.obol->isRenderRequested())
		s.obol->requestRender(reason);
}

/** Mirror the current RT view state into the Obol direct camera. */
static inline void
qgcanvas_sync_obol_camera(QgCanvasState &s)
{
	if (!s.obol || !s.v)
		return;

	const void *view_ctx = qg_legacy_view_to_context(s.v);
	(void)s.obol->syncCameraFromRtViewContext(view_ctx);
}

/** Render the Obol scene through SoOffscreenRenderer into a QImage. */
static inline void
qgcanvas_get_obol_viewport_image(QgCanvasState &s, const QWidget *w, QImage &img,
	bool consumeRenderRequest = false)
{
	img = QImage();
	if (!s.obol || !s.obol->getViewport())
		return;

	qgcanvas_sync_obol_viewport(s, w);

	const SbViewportRegion &region = s.obol->getViewportRegion();
	SbVec2s size = region.getViewportSizePixels();
	if (size[0] <= 0 || size[1] <= 0)
		return;

	SoOffscreenRenderer renderer(qgcanvas_obol_context_manager(), region);
	renderer.setComponents(SoOffscreenRenderer::RGB);
	renderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
	if (!s.obol->getViewport()->render(&renderer))
		return;

	unsigned char *buffer = renderer.getBuffer();
	if (!buffer)
		return;

	QImage raw(buffer, size[0], size[1], size[0] * 3, QImage::Format_RGB888);
	img = QImage(size[0], size[1], QImage::Format_RGB888);
	for (int y = 0; y < size[1]; y++)
		std::memcpy(img.scanLine(size[1] - 1 - y), raw.constScanLine(y),
		    static_cast<size_t>(size[0]) * 3);
	if (w)
		img.setDevicePixelRatio(w->devicePixelRatioF());
	if (consumeRenderRequest && s.obol->isRenderRequested())
		(void)s.obol->consumeRenderRequest(NULL);
}

/** Create the Obol view state every qtcad canvas exposes. */
static inline void
qgcanvas_init_obol(QgCanvasState &s, const QWidget *w)
{
	brlobol_init(qgcanvas_obol_context_manager());
	s.obol = new BRLObolViewController();

	SoSeparator *root = new SoSeparator;
	SoOrthographicCamera *camera = new SoOrthographicCamera;

	s.obol->setSceneRoot(root);
	s.obol->setCamera(camera);
	qgcanvas_sync_obol_viewport(s, w);
	qgcanvas_sync_obol_camera(s);
}

/** Destroy the Obol view state every qtcad canvas exposes. */
static inline void
qgcanvas_destroy_obol(QgCanvasState &s)
{
	delete s.obol;
	s.obol = nullptr;
}

/** Return true when the Obol scene root contains drawable/migrated content. */
static inline bool
qgcanvas_obol_scene_has_content(QgCanvasState &s)
{
	if (!s.obol)
		return false;
	SoNode *root = s.obol->getSceneRoot();
	if (!root || !root->isOfType(SoGroup::getClassTypeId()))
		return false;
	return static_cast<SoGroup *>(root)->getNumChildren() > 0;
}

static inline int
qgcanvas_find_obol_axes_child(SoGroup *group, const char *overlayId)
{
	if (!group || !overlayId)
		return -1;
	for (int i = 0; i < group->getNumChildren(); i++) {
		SoNode *node = group->getChild(i);
		if (!node || !node->isOfType(SoBRLAxes::getClassTypeId()))
			continue;
		SoBRLAxes *axes = static_cast<SoBRLAxes *>(node);
		if (std::strcmp(axes->overlayId.getValue().getString(), overlayId) == 0)
			return i;
	}
	return -1;
}

static inline int
qgcanvas_find_obol_grid_child(SoGroup *group, const char *overlayId)
{
	if (!group || !overlayId)
		return -1;
	for (int i = 0; i < group->getNumChildren(); i++) {
		SoNode *node = group->getChild(i);
		if (!node || !node->isOfType(SoBRLGrid::getClassTypeId()))
			continue;
		SoBRLGrid *grid = static_cast<SoBRLGrid *>(node);
		if (std::strcmp(grid->overlayId.getValue().getString(), overlayId) == 0)
			return i;
	}
	return -1;
}

static inline int
qgcanvas_find_obol_adc_child(SoGroup *group, const char *overlayId)
{
	if (!group || !overlayId)
		return -1;
	for (int i = 0; i < group->getNumChildren(); i++) {
		SoNode *node = group->getChild(i);
		if (!node || !node->isOfType(SoBRLADC::getClassTypeId()))
			continue;
		SoBRLADC *adc = static_cast<SoBRLADC *>(node);
		if (std::strcmp(adc->overlayId.getValue().getString(), overlayId) == 0)
			return i;
	}
	return -1;
}

static inline void
qgcanvas_sync_obol_axes(QgCanvasState &s,
	SoGroup *group,
	const char *overlayId,
	const struct rt_view_axes_state &state)
{
	const int childIndex = qgcanvas_find_obol_axes_child(group, overlayId);
	if (!state.draw) {
		if (childIndex >= 0) {
			group->removeChild(childIndex);
			qgcanvas_request_obol_render_if_idle(s, "faceplate");
		}
		return;
	}

	SoBRLAxes *axes = childIndex >= 0 ?
		static_cast<SoBRLAxes *>(group->getChild(childIndex)) :
		new SoBRLAxes;
	axes->overlayId = overlayId;
	axes->origin = SbVec3f(static_cast<float>(state.axes_pos[X]),
		static_cast<float>(state.axes_pos[Y]),
		static_cast<float>(state.axes_pos[Z]));
	axes->size = static_cast<float>(state.axes_size > SMALL_FASTF ?
		state.axes_size : 1.0);
	axes->visible = TRUE;
	axes->rebuildGeometry();
	if (childIndex < 0)
		group->addChild(axes);
	qgcanvas_request_obol_render_if_idle(s, "faceplate");
}

static inline void
qgcanvas_sync_obol_grid(QgCanvasState &s,
	SoGroup *group,
	const struct rt_view_grid_state &state)
{
	const char *overlayId = "faceplate::grid";
	const int childIndex = qgcanvas_find_obol_grid_child(group, overlayId);
	if (!state.draw) {
		if (childIndex >= 0) {
			group->removeChild(childIndex);
			qgcanvas_request_obol_render_if_idle(s, "faceplate");
		}
		return;
	}

	SoBRLGrid *grid = childIndex >= 0 ?
		static_cast<SoBRLGrid *>(group->getChild(childIndex)) :
		new SoBRLGrid;
	grid->overlayId = overlayId;
	grid->center = SbVec3f(static_cast<float>(state.anchor[X]),
		static_cast<float>(state.anchor[Y]),
		static_cast<float>(state.anchor[Z]));
	double spacing = state.res_h > SMALL_FASTF ? state.res_h : state.res_v;
	if (spacing <= SMALL_FASTF)
		spacing = 1.0;
	grid->spacing = static_cast<float>(spacing);
	int divisions = state.res_major_h > state.res_major_v ?
		state.res_major_h : state.res_major_v;
	grid->divisions = divisions > 0 ? divisions : 5;
	grid->visible = TRUE;
	grid->rebuildGeometry();
	if (childIndex < 0)
		group->addChild(grid);
	qgcanvas_request_obol_render_if_idle(s, "faceplate");
}

static inline void
qgcanvas_sync_obol_adc(QgCanvasState &s,
	SoGroup *group,
	const struct rt_view_adc_state &state)
{
	const char *overlayId = "faceplate::adc";
	const int childIndex = qgcanvas_find_obol_adc_child(group, overlayId);
	if (!state.draw) {
		if (childIndex >= 0) {
			group->removeChild(childIndex);
			qgcanvas_request_obol_render_if_idle(s, "faceplate");
		}
		return;
	}

	SoBRLADC *adc = childIndex >= 0 ?
		static_cast<SoBRLADC *>(group->getChild(childIndex)) :
		new SoBRLADC;
	adc->overlayId = overlayId;
	adc->center = SbVec3f(static_cast<float>(state.pos_model[X]),
		static_cast<float>(state.pos_model[Y]),
		static_cast<float>(state.pos_model[Z]));
	adc->angleDegrees = static_cast<float>(state.a1);
	adc->distance = static_cast<float>(state.dst > SMALL_FASTF ?
		state.dst : 1.0);
	adc->visible = TRUE;
	adc->rebuildGeometry();
	if (childIndex < 0)
		group->addChild(adc);
	qgcanvas_request_obol_render_if_idle(s, "faceplate");
}

static inline void
qgcanvas_sync_obol_faceplate(QgCanvasState &s)
{
	if (!s.obol || !s.v)
		return;

	SoNode *root = s.obol->getSceneRoot();
	if (!root || !root->isOfType(SoGroup::getClassTypeId()))
		return;
	SoGroup *group = static_cast<SoGroup *>(root);

	struct rt_view_grid_state grid = {};
	struct rt_view_axes_state modelAxes = {};
	struct rt_view_axes_state viewAxes = {};
	struct rt_view_adc_state adc = {};
	void *view_ctx = qg_legacy_view_to_context(s.v);
	(void)rt_view_context_grid_state_get(&grid, view_ctx);
	(void)rt_view_context_model_axes_state_get(&modelAxes, view_ctx);
	(void)rt_view_context_view_axes_state_get(&viewAxes, view_ctx);
	(void)rt_view_context_adc_state_get(&adc, view_ctx);

	qgcanvas_sync_obol_grid(s, group, grid);
	qgcanvas_sync_obol_axes(s, group, "faceplate::model_axes", modelAxes);
	qgcanvas_sync_obol_axes(s, group, "faceplate::view_axes", viewAxes);
	qgcanvas_sync_obol_adc(s, group, adc);
}

/** Render queued Obol work from a caller-owned current GL context. */
static inline SbBool
qgcanvas_render_obol_pending(QgCanvasState &s,
	SbBool clearWindow = FALSE,
	SbBool clearZBuffer = FALSE)
{
	if (!s.obol)
		return FALSE;
	return s.obol->renderPending(clearWindow, clearZBuffer, NULL);
}

/** Store current DM and view hashes in @p s. */
static inline void
qgcanvas_stash_hashes(QgCanvasState &s)
{
	s.prev_dhash = qg_legacy_view_dm_hash(s.dmp);
	s.prev_vhash = rt_view_context_hash(qg_legacy_view_to_context(s.v));
}

/** Request a semantic view refresh and wake the canvas backend. */
static inline void
qgcanvas_request_update(QgCanvasState &s, uint32_t flags)
{
	uint32_t requested = flags ? flags : RT_VIEW_REFRESH_ALL;
	if (requested & RT_VIEW_REFRESH_VIEW)
		qgcanvas_sync_obol_camera(s);
	qgcanvas_sync_obol_faceplate(s);
	if (s.v)
		rt_view_context_refresh_request(qg_legacy_view_to_context(s.v),
			requested);
	if (s.dmp)
		qg_legacy_view_dm_native_repaint_pending_set(s.dmp, 1);
}

/**
 * Compare current DM/view hashes against the stored values and update the
 * view refresh record and DM wakeup flag when differences are found.
 *
 * Returns true if any value changed.  The caller is responsible for calling
 * need_update() and emitting changed() — signal emission requires a QObject
 * context that QgCanvasState does not have.
 */
static inline bool
qgcanvas_diff_hashes_check(QgCanvasState &s)
{
	bool ret = false;
	unsigned long long c_dhash = qg_legacy_view_dm_hash(s.dmp);
	unsigned long long c_vhash =
		rt_view_context_hash(qg_legacy_view_to_context(s.v));

	if (qg_legacy_view_dm_native_repaint_pending_get(s.dmp)) {
		if (s.v)
			rt_view_context_refresh_request(qg_legacy_view_to_context(s.v),
			    RT_VIEW_REFRESH_FORCE);
		ret = true;
	}

	if (s.prev_dhash != c_dhash) {
		qgcanvas_request_update(s, RT_VIEW_REFRESH_FRAMEBUFFER | RT_VIEW_REFRESH_FORCE);
		ret = true;
	}
	if (s.prev_vhash != c_vhash) {
		qgcanvas_request_update(s, RT_VIEW_REFRESH_VIEW | RT_VIEW_REFRESH_DRAW);
		ret = true;
	}
	return ret;
}

/** Set the azimuth/elevation/twist of the view stored in @p s. */
static inline void
qgcanvas_aet(QgCanvasState &s, double a, double e, double t)
{
	if (!s.v)
		return;
	fastf_t aet_v[3];
	double  aetd[3] = {a, e, t};
	VMOVE(aet_v, aetd);
	void *view_ctx = qg_legacy_view_to_context(s.v);
	rt_view_context_aet_set(view_ctx, aet_v);
	rt_view_context_update(view_ctx);
	qgcanvas_sync_obol_camera(s);
}

/** Bind an external legacy view (or nullptr to revert to the widget-local view). */
static inline void
qgcanvas_set_view(QgCanvasState &s, qg_legacy_view *nv)
{
	if (!nv) {
		/* Revert to the widget-owned local view. */
		s.v = s.local_v;
		qgcanvas_sync_obol_camera(s);
		return;
	}
	s.v = nv;
	if (!s.dmp) {
		qgcanvas_sync_obol_camera(s);
		return;
	}
	qg_legacy_view_dm_bind(s.v, s.dmp);
	qg_legacy_view_dm_configure_window(s.dmp, 0);
	qg_legacy_view_dm_sync_dimensions(s.v, s.dmp);
	qgcanvas_sync_obol_camera(s);
}

#endif /* QGCANVASSTATE_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
