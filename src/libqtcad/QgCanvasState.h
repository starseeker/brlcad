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
 * files.
 *
 * This header is NOT part of the installed libqtcad public API.
 */

#ifndef QGCANVASSTATE_H
#define QGCANVASSTATE_H

#include "common.h"

#include "bu/str.h"
#include <climits>
#include <chrono>
#include <cmath>
#include <cstring>
#include <QImage>
#include <QOpenGLContext>
#include <QOpenGLExtraFunctions>
#include <QOpenGLFramebufferObject>
#include <QOpenGLWidget>
#include <QSize>
#include <QTimer>
#include <QWidget>

#include "BObol/BInit.h"
#include "BObol/BADC.h"
#include "BObol/BAxes.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BGrid.h"
#include "BObol/BHUDLabelOverlay.h"
#include "BObol/BLineLayerOverlay.h"
#include "BObol/BViewController.h"
#include "BObol/BViewLod.h"
#include "bv.h"
#include "ged/display.h"
#include "ged/view.h"
#include "QgObolContextManager.h"

#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>
#include <Inventor/actions/SoGLRenderAction.h>
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
 * v        – GED/bv view context created and owned by this canvas.  A canvas
 *            never swaps in externally-owned view state; the owning QgView
 *            exposes this context to GED and endpoint-facing clients.
 *
 */
struct QgCanvasState {
    /* ---- view plumbing ---- */
    struct bv_context *v = nullptr;         /* widget-owned view context */
    BObolViewController *obol = nullptr; /* Obol-canonical view controller */
    bool owns_obol = false;

    /* ---- hash tracking for incremental updates ---- */
    unsigned long long prev_dhash = 0;
    unsigned long long prev_vhash = 0;

    /* ---- input-binding flags ---- */
    bool use_default_keybindings   = true;
    bool use_default_mousebindings = true;
    int  lmouse_mode = -1;  /* set to BV_ADJUST_SCALE in canvas constructor */

    /* ---- widget-level tracking ---- */
    int    current = 1;     /* 1 = this view is active */
    int    x_prev = -INT_MAX;
    int    y_prev = -INT_MAX;
    double x_press_pos = -INT_MAX;
    double y_press_pos = -INT_MAX;
    bool   obol_paint_initialized = false;
    bool   fb_update_queued = false;
    bool   progressive_update_queued = false;
    bool   lod_progress_last_pending = false;
    bool   software_backend = false;
    QWidget *frame_request_widget = nullptr;
    SoOffscreenRenderer *offscreen_renderer = nullptr;
    QOpenGLFramebufferObject *presentation_fbo = nullptr;
    bool presentation_fbo_has_completed_frame = false;
    QImage last_completed_software_frame;
    std::chrono::steady_clock::time_point fps_last_publish;
    std::chrono::steady_clock::time_point lod_progress_last_publish;

    /* ---- per-canvas input handler ---- */
    QgCanvasInput input;
};

/* ------------------------------------------------------------------ */
/* Shared inline helpers (static so each TU gets its own copy)        */
/* ------------------------------------------------------------------ */

static inline struct bv_context *
qgcanvas_view_context_create(const char *name)
{
    struct bv_context *view_ctx = ged_view_context_bv(
	ged_view_context_create());
    if (!view_ctx)
	return NULL;

    struct bv_background_state background = BV_BACKGROUND_STATE_INIT;
    VSET(background.bottom, 110, 110, 110);
    VSET(background.top, 0, 0, 50);
    (void)bv_background_state_set(bv_context_view(view_ctx), &background);
    if (name)
	(void)bv_context_name_set(view_ctx, name);
    return view_ctx;
}

static inline void
qgcanvas_view_context_destroy(struct bv_context *view_ctx)
{
    if (view_ctx)
	ged_view_context_free(ged_view_context_from_bv(view_ctx));
}

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
qgcanvas_obol_context_manager(bool software = false)
{
    static QgObolContextManager hardwareManager(false);
    static QgObolContextManager softwareManager(true);
    return software ? static_cast<SoDB::ContextManager *>(&softwareManager) :
	static_cast<SoDB::ContextManager *>(&hardwareManager);
}

static inline void
qgcanvas_bind_obol_render_context(QgCanvasState &s)
{
    if (!s.obol)
	return;
    s.obol->setRenderContextManager(
	qgcanvas_obol_context_manager(s.software_backend));
}

static inline void
qgcanvas_request_obol_render_if_idle(QgCanvasState &s, const char *reason)
{
    if (s.obol && !s.obol->isRenderRequested())
	s.obol->requestPresentationRender(reason);
}

/* LoD completion may be reported by a worker thread after Qt's last paint.
 * Marshal the controller's frame request back to the canvas event loop so a
 * completed payload cannot remain hidden behind its startup proxy until the
 * next unrelated mouse or console event. */
static inline void
qgcanvas_obol_frame_requested(void *user_data, const char *UNUSED(reason))
{
    QgCanvasState *s = static_cast<QgCanvasState *>(user_data);
    QWidget *w = s ? s->frame_request_widget : nullptr;
    if (!w)
	return;
    const bool softwareBackend = s->software_backend;
    /* A software canvas has no swap/expose wake source of its own, and Qt may
     * coalesce update() indefinitely while qged's command/test layer owns a
     * nested event loop.  Queue the call onto the widget thread in both
     * cases, but make the OSMesa endpoint synchronous once it gets there.
     * Its Coin traversal still owns the controller's hard abort deadline. */

    QMetaObject::invokeMethod(w, [s, w, softwareBackend]() {
	/* update() is only a hint and may be folded into the paint whose
	 * completion raised this request.  A refinement barrier is a stronger
	 * contract: the controller cannot admit its successor cut until a new
	 * frame has actually traversed the published population.  Present that
	 * frame synchronously once this queued callback reaches the widget thread,
	 * on either backend.  Ordinary System GL refreshes remain coalescible. */
	if (softwareBackend ||
	    (s->obol && s->obol->hasPendingLodRefinementFrame()))
	    w->repaint();
	else
	    w->update();
    }, Qt::QueuedConnection);
}

static inline void
qgcanvas_bind_obol_frame_requests(QgCanvasState &s, QWidget *w)
{
    if (s.obol && w) {
	s.frame_request_widget = w;
	s.obol->setFrameRequestCallback(qgcanvas_obol_frame_requested, &s);
    }
}

static inline void
qgcanvas_unbind_obol_frame_requests(QgCanvasState &s, QWidget *w)
{
    if (s.obol && w) {
	s.obol->clearFrameRequestCallback(&s);
	s.frame_request_widget = nullptr;
    }
}

static inline void
qgcanvas_advance_obol_progressive(QgCanvasState &s)
{
    if (s.obol)
	(void)s.obol->advanceProgressiveWork(NULL, NULL);
}

static inline void
qgcanvas_queue_obol_progressive_update(QgCanvasState &s, QWidget *w)
{
    if (!s.obol || !w ||
	(!s.obol->hasProgressiveWorkPending() &&
	 !s.obol->isRenderRequested()) ||
	s.progressive_update_queued)
	return;

    s.progressive_update_queued = true;
    QTimer::singleShot(16, w, [&s, w]() {
	s.progressive_update_queued = false;
	if (!s.obol ||
	    (!s.obol->hasProgressiveWorkPending() &&
	     !s.obol->isRenderRequested()))
	    return;

	/* Poll background providers and LoD service state without forcing an
	 * expensive duplicate paint.  advanceProgressiveWork() requests a
	 * render only when presentation data or a retained PoP cut actually
	 * changes; otherwise keep this lightweight timer pump alive. */
	if (s.obol->hasProgressiveWorkPending())
	    (void)s.obol->advanceProgressiveWork(NULL, NULL);
	const bool needsPresentation = s.obol->isRenderRequested() ||
	    (s.lod_progress_last_pending &&
	     !s.obol->hasProgressiveWorkPending());
	if (needsPresentation) {
	    /* QWidget::update() may remain coalesced indefinitely while the GUI
	     * test/command layer is running a nested event loop and no unrelated
	     * expose event occurs.  A direct-GL canvas has Qt's swap machinery as
	     * an additional wake source; the software canvas does not.  Its
	     * bounded 16 ms endpoint pump therefore presents synchronously from
	     * this timer callback.  The callback is outside paintEvent, and each
	     * traversal retains Obol's hard frame deadline, so this cannot recurse
	     * or turn an expensive frame into an uninterruptible GUI loop. */
	    if (s.software_backend ||
		s.obol->hasPendingLodRefinementFrame()) {
		w->repaint();
	    } else
		w->update();
	}
	/* Do not depend on Qt delivering that paint to keep either the provider
	 * pump or an explicit frame request alive.  update() may be coalesced with
	 * the paint whose completion raised a calibration request.  That request
	 * does not necessarily represent background progressive work, so a timer
	 * gated only by hasProgressiveWorkPending() can strand it forever.  The
	 * next timer remains lightweight and retires itself as soon as a paint
	 * consumes the request and no provider work remains. */
	qgcanvas_queue_obol_progressive_update(s, w);
    });
}

/** Mirror the current RT view state into the Obol direct camera. */
static inline void
qgcanvas_sync_obol_camera(QgCanvasState &s)
{
    if (!s.obol || !s.v)
	return;

    (void)s.obol->syncCameraFromViewContext(s.v);
}

/* Seed a new controller from passive view state once.  Thereafter the
 * endpoint/controller owns background policy. */
static inline void
qgcanvas_initialize_obol_background(QgCanvasState &s)
{
	if (!s.obol || !s.v)
	return;
    struct bv_background_state background = BV_BACKGROUND_STATE_INIT;
    if (!bv_background_state_get(&background, bv_context_view_const(s.v)))
	return;
    s.obol->setBackgroundColors(
	SbColor(background.bottom[0] / 255.0f,
		background.bottom[1] / 255.0f,
		background.bottom[2] / 255.0f),
	SbColor(background.top[0] / 255.0f,
		background.top[1] / 255.0f,
		background.top[2] / 255.0f));
}

/** Render the Obol scene through SoOffscreenRenderer into a QImage. */
static inline void
qgcanvas_get_obol_viewport_image(QgCanvasState &s, const QWidget *w, QImage &img,
				 bool consumeRenderRequest = false,
				 bool borrowRendererBuffer = false,
				 bool recordPresentationTiming = false,
				 bool *completedPresentation = nullptr)
{
    img = QImage();
    if (completedPresentation)
	*completedPresentation = false;
    if (!s.obol || !s.obol->getViewport())
	return;

    qgcanvas_sync_obol_viewport(s, w);
    /* Only the endpoint's actual paint may drive the progressive state
     * machine.  A diagnostic/export traversal must be observational: making
     * a checkpoint advance admission can reopen a budget probe after the
     * scripted idle barrier and leave the supposedly final report pending. */
    if (recordPresentationTiming)
	qgcanvas_advance_obol_progressive(s);
    /* Endpoint-owned image producers, including the retained librt engine,
     * publish completed worker frames through this host-thread hook.  The Qt
     * canvases render directly instead of using
     * BObolViewController::renderToImage(), so they must perform the same
     * presentation synchronization explicitly before traversing the scene. */
    s.obol->synchronizePresentation();
    /* The software canvas traverses directly rather than calling
     * BObolViewController::renderPending().  Consume the request being
     * attempted before traversal, just as renderPending() does.  Leaving the
     * old boolean set until after completeRenderTiming() prevents a follow-up
     * calibration request from producing an empty-to-pending wake edge: the
     * serial preserves that request, but Qt never learns it is runnable.
     * Consuming first also makes requests published during traversal distinct
     * and self-waking. */
    SbBool lodCapacityRelevant = TRUE;
    if (recordPresentationTiming && consumeRenderRequest)
	(void)s.obol->consumeRenderRequest(NULL, &lodCapacityRelevant);

    const SbViewportRegion &region = s.obol->getViewportRegion();
    SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return;

    if (!s.offscreen_renderer) {
	s.offscreen_renderer = new SoOffscreenRenderer(
	    qgcanvas_obol_context_manager(s.software_backend), region);
    } else {
	s.offscreen_renderer->setViewportRegion(region);
    }
    SoOffscreenRenderer &renderer = *s.offscreen_renderer;
    renderer.setComponents(SoOffscreenRenderer::RGB_TRANSPARENCY);
    SoGLRenderAction *action = renderer.getGLRenderAction();
    if (action) {
	action->setSmoothing(s.obol->isAntialiasingEnabled());
	action->setNumPasses(1);
    }
    renderer.setBackgroundColor(s.obol->getBackgroundBottomColor());
    if (s.obol->getBackgroundBottomColor() !=
	s.obol->getBackgroundTopColor())
	renderer.setBackgroundGradient(s.obol->getBackgroundBottomColor(),
	    s.obol->getBackgroundTopColor());
    struct DeadlineContext {
	BObolViewController *controller = nullptr;
	uint64_t deadline = 0;
	SoGLRenderAction::SoGLRenderAbortCB *previous = nullptr;
	void *previousData = nullptr;
    } deadlineContext;
    const auto deadlineCallback = [](void *userData) {
	DeadlineContext *context = static_cast<DeadlineContext *>(userData);
	if (!context)
	    return SoGLRenderAction::CONTINUE;
	if (context->previous) {
	    const SoGLRenderAction::AbortCode prior =
		(*context->previous)(context->previousData);
	    if (prior != SoGLRenderAction::CONTINUE)
		return prior;
	}
	return context->controller && context->deadline &&
	    context->controller->beginRenderTiming() >= context->deadline ?
	    SoGLRenderAction::ABORT : SoGLRenderAction::CONTINUE;
    };
    const uint64_t started = s.obol->beginRenderTiming();
    const uint64_t deadlineDuration =
	recordPresentationTiming && lodCapacityRelevant ?
	s.obol->getCurrentPresentationFrameDeadline() : 0;
    if (action && deadlineDuration) {
	deadlineContext.controller = s.obol;
	deadlineContext.deadline =
	    started > UINT64_MAX - deadlineDuration ? UINT64_MAX :
	    started + deadlineDuration;
	action->getAbortCallback(
	    deadlineContext.previous, deadlineContext.previousData);
	action->setAbortCallback(deadlineCallback, &deadlineContext);
    }
    BObolViewLodState *presentationState =
	s.obol ? s.obol->getViewLodState() : nullptr;
    const uint64_t cadExecutionBefore = presentationState ?
	presentationState->cadPresentationExecutionSerial() : 0;
    const uint64_t cadPreparationBefore = presentationState ?
	presentationState->cadPresentationPreparationSerial() : 0;
    if (presentationState)
	presentationState->beginCadPresentationFrame();
    const SbBool rendered = renderer.render(s.obol->getRenderRoot());
    const uint64_t completed = s.obol->beginRenderTiming();
    const uint64_t cadExecutionAfter = presentationState ?
	presentationState->cadPresentationExecutionSerial() : 0;
    const uint64_t cadPreparationAfter = presentationState ?
	presentationState->cadPresentationPreparationSerial() : 0;
    if (presentationState)
	presentationState->refreshCadPresentationFrameStatus();
    const SbBool cadFrameIncomplete = presentationState &&
	!presentationState->lastCadPresentationFrameExact();
    const SbBool interrupted =
	(action && action->hasTerminated()) || cadFrameIncomplete;
    if (action && deadlineDuration)
	action->setAbortCallback(
	    deadlineContext.previous, deadlineContext.previousData);
    /* Image export/checkpoint readback is a second traversal, not a frame the
     * viewport needed to present.  Feeding it into the scene LoD capacity
     * estimator makes screenshot frequency and PNG test checkpoints alter
     * the terminal PoP cut.  QgSW's paint path opts in below; diagnostic
     * image producers deliberately do not. */
    if (recordPresentationTiming && !interrupted)
	s.obol->completeRenderTiming(started, lodCapacityRelevant);
    if (interrupted) {
	if (recordPresentationTiming) {
	    s.obol->notePresentationRenderInterrupted(
		completed > started ? completed - started : 1,
		cadExecutionAfter != cadExecutionBefore ? TRUE : FALSE,
		cadPreparationAfter != cadPreparationBefore ? TRUE : FALSE,
		lodCapacityRelevant);
	}
	if (!s.last_completed_software_frame.isNull()) {
	    img = s.last_completed_software_frame;
	} else if (rendered) {
	    /* Before the first exact frame exists, a deadline-bounded structural
	     * prefix is still the fastest useful cold-start presentation.  Show
	     * that provisional image, but do not promote it to the retained
	     * completed frame: later interruptions must never regress a mesh view
	     * back to this box-only prefix. */
	    unsigned char *buffer = renderer.getBuffer();
	    if (buffer) {
		QImage raw(buffer, size[0], size[1], size[0] * 4,
		    QImage::Format_RGBX8888);
		if (borrowRendererBuffer) {
		    img = raw.copy();
		} else {
		    img = QImage(size[0], size[1], QImage::Format_RGBX8888);
		    for (int y = 0; y < size[1]; y++)
			std::memcpy(img.scanLine(size[1] - 1 - y),
			    raw.constScanLine(y),
			    static_cast<size_t>(size[0]) * 4);
		}
	    }
	}
	if (w && !img.isNull())
	    img.setDevicePixelRatio(w->devicePixelRatioF());
	return;
    }
    if (!rendered)
	return;

    unsigned char *buffer = renderer.getBuffer();
    if (!buffer)
	return;

    QImage raw(buffer, size[0], size[1], size[0] * 4,
	QImage::Format_RGBX8888);
    if (borrowRendererBuffer) {
	/* Keep one immutable completed frame.  A deadline-aborted OSMesa
	 * traversal has already overwritten its renderer-owned buffer, so a
	 * borrowed view cannot preserve the last coherent presentation. */
	s.last_completed_software_frame = raw.copy();
	img = s.last_completed_software_frame;
    } else {
	img = QImage(size[0], size[1], QImage::Format_RGBX8888);
	for (int y = 0; y < size[1]; y++)
	    std::memcpy(img.scanLine(size[1] - 1 - y), raw.constScanLine(y),
		static_cast<size_t>(size[0]) * 4);
    }
    if (w)
	img.setDevicePixelRatio(w->devicePixelRatioF());
    if (completedPresentation)
	*completedPresentation = true;
}

/** Create the Obol view state every qtcad canvas exposes. */
static inline void
qgcanvas_init_obol(QgCanvasState &s, QWidget *w,
	bool software_backend, BObolViewController *controller = nullptr,
	bool create_controller = true)
{
    s.software_backend = software_backend;
    /* Coin has one process-global fallback manager.  Keep it stable and use
    * explicit per-renderer managers below: QgGL binds its direct-rendering
    * manager to its render action, while QgSW supplies its private OSMesa manager to
     * every SoOffscreenRenderer.  Replacing the global manager as canvases
     * are constructed can otherwise cross-contaminate those backends. */
    bobol_init(NULL);
    s.obol = controller;
    s.owns_obol = false;
    if (!s.obol && create_controller) {
	s.obol = new BObolViewController();
	s.owns_obol = true;
    }
    if (!s.obol)
	return;
    qgcanvas_bind_obol_render_context(s);
	qgcanvas_bind_obol_frame_requests(s, w);

    SoSeparator *root = new SoSeparator;
    SoOrthographicCamera *camera = new SoOrthographicCamera;

    s.obol->setSceneRoot(root);
    s.obol->setCamera(camera);
    qgcanvas_sync_obol_viewport(s, w);
    qgcanvas_sync_obol_camera(s);
}

/** Destroy the Obol view state every qtcad canvas exposes. */
static inline void
qgcanvas_destroy_obol(QgCanvasState &s, QWidget *w)
{
    delete s.offscreen_renderer;
    s.offscreen_renderer = nullptr;
    delete s.presentation_fbo;
    s.presentation_fbo = nullptr;
    s.last_completed_software_frame = QImage();
    if (s.obol && s.obol->getRenderContextManager() ==
	    qgcanvas_obol_context_manager(s.software_backend))
	s.obol->setRenderContextManager(NULL);
    qgcanvas_unbind_obol_frame_requests(s, w);
    if (s.owns_obol)
	delete s.obol;
    s.obol = nullptr;
    s.owns_obol = false;
}

/** Replace the canvas-owned controller with a borrowed endpoint controller. */
static inline void
qgcanvas_bind_obol_controller(QgCanvasState &s, QWidget *w,
	BObolViewController *controller)
{
    if (s.obol == controller) {
	qgcanvas_bind_obol_render_context(s);
	qgcanvas_bind_obol_frame_requests(s, w);
	return;
    }

    if (s.obol) {
	qgcanvas_unbind_obol_frame_requests(s, w);
	s.obol->setRenderContextManager(NULL);
    }
    delete s.offscreen_renderer;
    s.offscreen_renderer = nullptr;
    if (s.owns_obol)
	delete s.obol;
    s.obol = controller;
    s.owns_obol = false;
    s.obol_paint_initialized = false;

    if (!s.obol)
	return;
    qgcanvas_bind_obol_render_context(s);
	qgcanvas_bind_obol_frame_requests(s, w);
    if (!s.obol->getSceneRoot())
	s.obol->setSceneRoot(new SoSeparator);
    if (!s.obol->getCamera())
	s.obol->setCamera(new SoOrthographicCamera);
    qgcanvas_sync_obol_viewport(s, w);
    qgcanvas_sync_obol_camera(s);
    /* The endpoint's GED view seeds renderer policy at attachment.  Host
     * canvases have their own passive view state, so applying it here would
     * overwrite a retained endpoint property during host replacement. */
    s.obol->requestRender("qt-controller-bind");
}

static inline bool
qgcanvas_obol_node_has_drawable_content(SoNode *node)
{
    if (!node)
	return false;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()) ||
	node->isOfType(SoBRLLineLayerOverlay::getClassTypeId()) ||
	node->isOfType(SoBRLHUDLabelOverlay::getClassTypeId()) ||
	node->isOfType(SoBRLGrid::getClassTypeId()) ||
	node->isOfType(SoBRLAxes::getClassTypeId()) ||
	node->isOfType(SoBRLADC::getClassTypeId()))
	return true;

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return true;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (qgcanvas_obol_node_has_drawable_content(group->getChild(i)))
	    return true;
    }

    return false;
}

/** Return true when the Obol scene contains drawable/migrated content. */
static inline bool
qgcanvas_obol_scene_has_content(QgCanvasState &s)
{
    if (!s.obol)
	return false;
    if (qgcanvas_obol_node_has_drawable_content(s.obol->getRenderSceneRoot()))
	return true;
    return qgcanvas_obol_node_has_drawable_content(s.obol->getSceneRoot());
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
	if (bu_strcmp(axes->overlayId.getValue().getString(), overlayId) == 0)
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
	if (bu_strcmp(grid->overlayId.getValue().getString(), overlayId) == 0)
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
	if (bu_strcmp(adc->overlayId.getValue().getString(), overlayId) == 0)
	    return i;
    }
    return -1;
}

static inline void
qgcanvas_sync_obol_axes(QgCanvasState &s,
			SoGroup *group,
			const char *overlayId,
			const struct bv_axes_state &state)
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
			const struct bv_grid_state &state,
			struct ged_view_context *view_ctx)
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
    (void)bobol_grid_configure_from_view_context(grid, &state, view_ctx);
    if (childIndex < 0)
	group->addChild(grid);
    qgcanvas_request_obol_render_if_idle(s, "faceplate");
}

static inline void
qgcanvas_sync_obol_adc(QgCanvasState &s,
		       SoGroup *group,
		       const struct bv_adc_state &state)
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
    adc->lineColor = SbColor(state.line_color[0] / 255.0f,
	state.line_color[1] / 255.0f, state.line_color[2] / 255.0f);
    adc->tickColor = SbColor(state.tick_color[0] / 255.0f,
	state.tick_color[1] / 255.0f, state.tick_color[2] / 255.0f);
    adc->lineWidth = state.line_width > 0 ? state.line_width : 1;
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

    struct ged_view_context *view_ctx = ged_view_context_from_bv(s.v);
    struct ged *gedp = ged_view_context_owner(view_ctx);
    if (gedp) {
	(void)ged_view_faceplate_sync(gedp, view_ctx);
	return;
    }

    SoNode *root = s.obol->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return;
    SoGroup *group = static_cast<SoGroup *>(root);

    struct bv_grid_state grid = {};
    struct bv_axes_state modelAxes = {};
    struct bv_axes_state viewAxes = {};
    struct bv_adc_state adc = {};
    const struct bv *view = bv_context_view_const(s.v);
    (void)bv_grid_state_get(&grid, view);
    (void)bv_model_axes_state_get(&modelAxes, view);
    (void)bv_view_axes_state_get(&viewAxes, view);
    (void)bv_adc_state_get(&adc, view);

    qgcanvas_sync_obol_grid(s, group, grid, view_ctx);
    qgcanvas_sync_obol_axes(s, group, "faceplate::model_axes", modelAxes);
    qgcanvas_sync_obol_axes(s, group, "faceplate::view_axes", viewAxes);
    qgcanvas_sync_obol_adc(s, group, adc);
}

/** Record actual Qt presentation cadence and refresh the label sparingly. */
static inline void
qgcanvas_frame_complete(QgCanvasState &s, QWidget *w)
{
    if (!s.v || !s.obol || !w)
	return;

    struct bv *view = bv_context_view(s.v);
    s.obol->noteFramePresented();
    const uint64_t presentation_interval =
	s.obol->getDisplayedPresentationIntervalNanoseconds();
    if (presentation_interval)
	(void)bv_frametime_set(view, presentation_interval);

    const auto now = std::chrono::steady_clock::now();
    const bool lod_pending =
	s.obol->hasProgressiveWorkPending() ||
	s.obol->isLodInteractionActive();
    const bool lod_state_changed =
	lod_pending != s.lod_progress_last_pending;
    const bool lod_first =
	s.lod_progress_last_publish.time_since_epoch().count() == 0;
    const bool lod_publish = lod_state_changed ||
	(lod_pending && (lod_first ||
	    std::chrono::duration_cast<std::chrono::milliseconds>(
		now - s.lod_progress_last_publish).count() >= 100));
    if (lod_publish) {
	s.lod_progress_last_pending = lod_pending;
	s.lod_progress_last_publish = now;
	qgcanvas_sync_obol_faceplate(s);
    }

    /*
     * completeRenderTiming() may install an unchanged calibration replay
     * after the progressive pump has already transitioned to idle.  Queue it
     * before the optional FPS/HUD reporting exits below; frame scheduling is
     * a renderer contract and must not depend on whether the user enabled an
     * FPS label.
     */
    if (s.obol->isRenderRequested() &&
	(lod_publish || s.obol->hasPendingLodRefinementFrame()))
	w->update();

    struct bv_params_state params = BV_PARAMS_STATE_INIT;
    if (!bv_params_state_get(&params, view) || !params.draw ||
	!params.draw_fps)
	return;

    if (!presentation_interval)
	return;

    const bool first =
	s.fps_last_publish.time_since_epoch().count() == 0;
    const bool publish = first ||
	std::chrono::duration_cast<std::chrono::milliseconds>(
	    now - s.fps_last_publish).count() >= 250;
    if (publish) {
	s.fps_last_publish = now;
	qgcanvas_sync_obol_faceplate(s);
    }

    /* FPS faceplate publication may itself install one final render. */
    if (s.obol->isRenderRequested())
	w->update();

}

/** Render queued Obol work from a caller-owned current GL context. */
static inline SbBool
qgcanvas_render_obol_pending(QgCanvasState &s, QOpenGLWidget *widget,
			     SbBool clearWindow = FALSE,
			     SbBool clearZBuffer = FALSE)
{
    if (!s.obol || !widget)
	return FALSE;
    if (!s.software_backend && !QOpenGLContext::currentContext())
	return FALSE;
    qgcanvas_bind_obol_render_context(s);
    const QSize renderSize = qgcanvas_render_size(widget);
    if (renderSize.isEmpty())
	return FALSE;
    if (!s.presentation_fbo || s.presentation_fbo->size() != renderSize) {
	delete s.presentation_fbo;
	s.presentation_fbo_has_completed_frame = false;
	QOpenGLFramebufferObjectFormat format;
	format.setAttachment(QOpenGLFramebufferObject::CombinedDepthStencil);
	format.setSamples(0);
	s.presentation_fbo = new QOpenGLFramebufferObject(renderSize, format);
    }
    if (!s.presentation_fbo || !s.presentation_fbo->isValid()) {
	delete s.presentation_fbo;
	s.presentation_fbo = nullptr;
	s.presentation_fbo_has_completed_frame = false;
	return s.obol->renderPending(clearWindow, clearZBuffer, NULL);
    }

    s.presentation_fbo->bind();
    const SbBool rendered =
	s.obol->renderPending(clearWindow, clearZBuffer, NULL);
    QOpenGLContext *context = QOpenGLContext::currentContext();
    QOpenGLExtraFunctions *gl = context ? context->extraFunctions() : nullptr;
    if (gl) {
	gl->glBindFramebuffer(GL_READ_FRAMEBUFFER,
	    s.presentation_fbo->handle());
	gl->glBindFramebuffer(GL_DRAW_FRAMEBUFFER,
	    widget->defaultFramebufferObject());
	if (rendered || !s.presentation_fbo_has_completed_frame)
	    gl->glBlitFramebuffer(0, 0, renderSize.width(), renderSize.height(),
		0, 0, renderSize.width(), renderSize.height(),
		GL_COLOR_BUFFER_BIT, GL_NEAREST);
	gl->glBindFramebuffer(GL_FRAMEBUFFER,
	    widget->defaultFramebufferObject());
    } else {
	s.presentation_fbo->release();
    }
    if (rendered)
	s.presentation_fbo_has_completed_frame = true;
    return rendered;
}

/** Store current view hashes in @p s. */
static inline void
qgcanvas_stash_hashes(QgCanvasState &s)
{
    s.prev_dhash = 0;
    s.prev_vhash = bv_hash(bv_context_view_const(s.v));
}

/** Request a semantic view refresh and wake the canvas backend. */
static inline void
qgcanvas_request_update(QgCanvasState &s, uint32_t flags)
{
    uint32_t requested = flags ? flags : BV_REFRESH_ALL;
    if (requested & BV_REFRESH_VIEW)
	qgcanvas_sync_obol_camera(s);
    qgcanvas_sync_obol_faceplate(s);
    if (s.v)
	bv_refresh_request(bv_context_view(s.v), requested);
}

/**
 * Compare current view hashes against the stored values and update the
 * view refresh record when differences are found.
 *
 * Returns true if any value changed.  The caller is responsible for calling
 * need_update() and emitting changed() — signal emission requires a QObject
 * context that QgCanvasState does not have.
 */
static inline bool
qgcanvas_diff_hashes_check(QgCanvasState &s)
{
    bool ret = false;
    const struct bv *view = bv_context_view_const(s.v);
    unsigned long long c_vhash = bv_hash(view);

    /* Commands may change renderer policy without changing the camera hash
     * (lighting, normal presentation, overlays, and similar state).  Honor
     * libbv's explicit refresh request as well as transform differences so
     * qged's post-command DisplayDiff path actually schedules a canvas paint. */
    if (s.prev_vhash != c_vhash || bv_refresh_dirty_get(view)) {
	qgcanvas_request_update(s, BV_REFRESH_VIEW | BV_REFRESH_DRAW);
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
    bv_aet_set(bv_context_view(s.v), aet_v);
    bv_context_update(s.v, BV_CONTEXT_CHANGED_VIEW);
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
