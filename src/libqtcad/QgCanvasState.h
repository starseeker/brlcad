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
#include <atomic>
#include <climits>
#include <chrono>
#include <cmath>
#include <cstring>
#include <exception>
#include <optional>
#include <utility>
#include <vector>
#include <QImage>
#include <QColor>
#include <QFontMetrics>
#include <QOpenGLContext>
#include <QOpenGLExtraFunctions>
#include <QOpenGLFramebufferObject>
#include <QOpenGLWidget>
#include <QPainter>
#include <QSize>
#include <QString>
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
#include "BObol/BViewStore.h"
#include "bv.h"
#include "ged/display.h"
#include "ged/view.h"
#include "QgObolContextManager.h"

#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/SbColor.h>
#include <Inventor/misc/SoChildList.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/tools/SbModernUtils.h>

#include "QgCanvasInput.h"

static inline QImage
qgcanvas_flip_vertical(const QImage &image)
{
#if QT_VERSION >= QT_VERSION_CHECK(6, 9, 0)
    return image.flipped(Qt::Vertical);
#else
    return image.mirrored(false, true);
#endif
}

struct QgLodProgressOverlayState {
    bool visible = false;
    bool determinate = false;
    bool etaVisible = false;
    bool refinementCycleBased = false;
    bool terminalReady = false;
    bool resourceLimited = false;
    int percent = 0;
    uint64_t estimatedRemainingMilliseconds = 0;
    uint64_t remainingRefinementCycles = 0;
    unsigned int animationStep = 0;
    int displayClass = BOBOL_LOD_PROGRESS_DISPLAY_IDLE;
    QString title;
    QString detail;

    bool operator==(const QgLodProgressOverlayState &other) const
    {
	return this->visible == other.visible &&
	    this->determinate == other.determinate &&
	    this->etaVisible == other.etaVisible &&
	    this->refinementCycleBased == other.refinementCycleBased &&
	    this->terminalReady == other.terminalReady &&
	    this->resourceLimited == other.resourceLimited &&
	    this->percent == other.percent &&
	    this->estimatedRemainingMilliseconds ==
		other.estimatedRemainingMilliseconds &&
	    this->remainingRefinementCycles ==
		other.remainingRefinementCycles &&
	    this->animationStep == other.animationStep &&
	    this->displayClass == other.displayClass &&
	    this->title == other.title && this->detail == other.detail;
    }

    bool operator!=(const QgLodProgressOverlayState &other) const
    {
	return !(*this == other);
    }
};

/* Presentation-only continuity for the native progress card.  Controller
 * convergence remains authoritative; this state smooths qualified estimates
 * and volatile diagnostics without carrying a percentage across a loss of
 * forecast confidence. */
struct QgLodProgressPresentationState {
    bool active = false;
    bool determinateLatched = false;
    bool etaVisible = false;
    bool finalizing = false;
    uint64_t episodeRevision = 0;
    uint64_t lastElapsedMilliseconds = 0;
    uint64_t etaConfidenceStartMilliseconds = 0;
    uint64_t smoothedCompletionMilliseconds = 0;
    uint64_t stableDetailUpdateMilliseconds = 0;
    unsigned int consistentEtaSamples = 0;
    int percentFloor = 0;
    QString stableDetail;

    void reset(void)
    {
	*this = QgLodProgressPresentationState();
    }
};

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
    uint64_t prev_frame_revision = 0;
    unsigned long long faceplate_view_hash = 0;
    uint64_t faceplate_frame_revision = 0;
    bool faceplate_sync_initialized = false;

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
    QTimer progressive_update_timer;
    bool   progressive_update_connected = false;
    int    progressive_update_delay_msec = 16;
    /* Set only by the built-in camera drag action.  Selection/edit gestures
     * use the same Qt buttons but must not coarsen or restart scene LoD. */
    bool   lod_pointer_interaction_active = false;
    bool   lod_progress_idle_tail_pending = false;
    std::atomic<bool> frame_request_dispatch_queued {false};
    BObolLodProgressDisplayStatus lod_progress_last_state;
    QgLodProgressOverlayState lod_progress_overlay;
    QgLodProgressPresentationState lod_progress_presentation;
    bool lod_progress_overlay_dirty = false;
    std::chrono::steady_clock::time_point lod_progress_overlay_last_request;
    bool   software_backend = false;
    QWidget *frame_request_widget = nullptr;
    SoOffscreenRenderer *offscreen_renderer = nullptr;
    /* Direct GL renders into staging and promotes it only after an exact
     * traversal.  A deadline-aborted traversal may overwrite staging, but it
     * can never corrupt the completed framebuffer presented to Qt. */
    QOpenGLFramebufferObject *presentation_fbo = nullptr;
    QOpenGLFramebufferObject *presentation_staging_fbo = nullptr;
    bool presentation_fbo_has_completed_frame = false;
    /* Immutable renderer-native (bottom-up) pixels from the last completed
     * software presentation.  QgSW's paint hot path consumes this orientation
     * directly; observational Qt image callers must receive a flipped copy. */
    QImage last_completed_software_frame;
    /* The visible image may be a provisional first frame or a failure clear.
     * Keep its identity separate from the exact frame used for abort fallback.
     * Completed presentations share their QImage storage with that fallback. */
    QImage last_presented_software_frame;
    /* Completion feedback can publish newer HUD records before the next
     * paint. Keep the traversed revision with the completed pixels; a clear
     * or provisional first image has no presented feature certificate. */
    uint64_t completed_feature_revision = 0;
    std::optional<uint64_t> presented_feature_revision;
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

/** True when the software endpoint can repaint Qt from its immutable last
 * completed presentation without entering Coin/OSMesa again. */
static inline bool
qgcanvas_has_completed_software_frame(const QgCanvasState &s,
	const QWidget *w)
{
    if (!w || s.last_completed_software_frame.isNull())
	return false;
    const QSize renderSize = qgcanvas_render_size(w);
    /* Pixel dimensions alone are insufficient.  A retained image produced
     * for DPR 1.5 and later painted with Qt's default DPR 1 contract is
     * interpreted as a logical-size image: the model and faceplate both grow
     * and the viewport appears to jump.  Treat the scale metadata as part of
     * the immutable presentation's endpoint identity. */
    return s.last_completed_software_frame.size() == renderSize &&
	qFuzzyCompare(s.last_completed_software_frame.devicePixelRatio(),
	    w->devicePixelRatioF());
}

/** True when the direct-GL endpoint has a completed presentation for the
 * current physical viewport. */
static inline bool
qgcanvas_has_completed_gl_frame(const QgCanvasState &s, const QWidget *w)
{
    if (!w || !s.presentation_fbo_has_completed_frame ||
	!s.presentation_fbo)
	return false;
    return s.presentation_fbo->size() == qgcanvas_render_size(w);
}

static inline void qgcanvas_queue_obol_progressive_update(
    QgCanvasState &s, QWidget *w);
static inline bool qgcanvas_sync_obol_lod_progress(
    QgCanvasState &s, bool allowPeriodic = true);
static inline void qgcanvas_request_lod_overlay_repaint(
    QgCanvasState &s, QWidget *w, bool force = false);

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
    /* Controller transitions may originate on worker and owner threads, and
     * one bounded pump can publish several internal obligation edges.  The
     * callback is only a wake hint for the standing level below: coalesce it
     * before posting to Qt so a 50k result wave cannot flood the GUI queue
     * with thousands of stale callbacks. */
    if (s->frame_request_dispatch_queued.exchange(true,
	std::memory_order_acq_rel))
	return;
    const bool softwareBackend = s->software_backend;
    /* A progressive wake is not necessarily a presentation request.  It may
     * mean only that another bounded provider/planning slice can run.  Queue
     * the callback onto the widget thread in both cases; the controller's
     * explicit render latch below remains the sole authority for repainting.
     * This distinction is essential for software rendering: one OSMesa
     * traversal per 2k-entry planning window turned a 50k current-view delta
     * into hundreds of unchanged 40--80 ms frames. */

	QMetaObject::invokeMethod(w, [s, w, softwareBackend]() {
	s->frame_request_dispatch_queued.store(false,
	    std::memory_order_release);
	/* Frame callbacks are the producer-to-host wake edge.  Service one
	 * bounded provider slice directly on the widget thread and explicitly
	 * arm the continuing timer; relying on repaint() to do both is not
	 * sufficient when Qt suppresses a paint for an obscured/initializing
	 * widget or while a parent temporarily disables updates. */
	BObolHostWorkSnapshot work = s->obol ?
	    s->obol->getHostWorkSnapshot() : BObolHostWorkSnapshot();
	if (s->obol && work.pumpPending())
	    (void)s->obol->advanceProgressiveWork(NULL, NULL);
	/* The wake itself may be the edge from a terminal retained frame into
	 * background compaction or a newly opened refinement obligation.  Publish
	 * that semantic transition before deciding whether geometry needs a paint;
	 * otherwise the retained faceplate can continue to promise "View ready"
	 * until the following timer slice even though the progress bar/controller
	 * has already re-entered active work. */
	const bool lodProgressPublish =
	    qgcanvas_sync_obol_lod_progress(*s, false);
	if (lodProgressPublish)
	    qgcanvas_request_obol_render_if_idle(*s, "lod-progress-wake");
	/* A partial immutable result may install a refinement barrier while its
	 * adaptive publication deadline is still accumulating a batch.  The
	 * barrier deliberately has no render request in that interval.  Once the
	 * controller replaces the timer witness with an explicit request, present
	 * synchronously for OSMesa and for a strong refinement barrier; ordinary
	 * System GL refreshes remain coalescible. */
	work = s->obol ? s->obol->getHostWorkSnapshot() :
	    BObolHostWorkSnapshot();
	if (s->obol && work.renderPending()) {
	    if (softwareBackend ||
		s->obol->hasPendingLodRefinementFrame())
		w->repaint();
	    else
		w->update();
	} else
	    qgcanvas_request_lod_overlay_repaint(*s, w);
	qgcanvas_queue_obol_progressive_update(*s, w);
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
    if (!s.obol)
	return;
    const BObolHostWorkSnapshot work = s.obol->getHostWorkSnapshot();
    /* A presentation-only refresh is already complete work: it makes a
     * retained style/overlay mutation visible.  Advancing an otherwise idle
     * coordinator here made selection drags reopen LoD policy even though
     * neither the camera nor the managed population changed. */
    if (work.pumpPending() || work.capacitySampleRequested())
	(void)s.obol->advanceProgressiveWork(NULL, NULL);
}

static inline void
qgcanvas_queue_obol_progressive_update(QgCanvasState &s, QWidget *w)
{
    const BObolHostWorkSnapshot initialWork = s.obol ?
	s.obol->getHostWorkSnapshot() : BObolHostWorkSnapshot();
    if (!s.obol || !w ||
	(initialWork.flags == BOBOL_HOST_WORK_NONE &&
	 !s.lod_progress_idle_tail_pending) ||
	s.progressive_update_timer.isActive())
	return;

    /* A controller pump is itself bounded and returns to Qt between slices.
     * Keep the ordinary 16 ms cadence for worker polling, debounce, and frame
     * requests, but let measured owner-thread planning resume after a 1 ms
     * event-loop yield.  Large compact-source scans otherwise paid 16 ms of
     * idle time after every 8 ms slice and could miss a 60-second liveness
     * gate despite having no worker, renderer, or cache work outstanding. */
    const int delay = s.progressive_update_delay_msec <= 1 ? 1 : 16;
    if (!s.progressive_update_connected) {
	s.progressive_update_timer.setSingleShot(true);
	QObject::connect(&s.progressive_update_timer, &QTimer::timeout, w,
	    [&s, w]() {
	BObolHostWorkSnapshot work = s.obol ?
	    s.obol->getHostWorkSnapshot() : BObolHostWorkSnapshot();
	if (!s.obol)
	    return;
	if (work.flags == BOBOL_HOST_WORK_NONE) {
	    s.progressive_update_delay_msec = 16;
	    /* One delayed no-work observation closes the host/HUD race.  The last
	     * pump may clear its work flag before the coordinator publishes IDLE;
	     * this tail never advances geometry and never reschedules itself. */
	    s.lod_progress_idle_tail_pending = false;
	    const bool lodProgressPublish =
		qgcanvas_sync_obol_lod_progress(s, false);
	    if (lodProgressPublish) {
		qgcanvas_request_obol_render_if_idle(s, "lod-progress-idle");
		if (s.software_backend)
		    w->repaint();
		else
		    w->update();
		/* Synchronizing the terminal HUD is allowed to create a fresh
		 * presentation request (normally "view-feature-store").  That is
		 * new level-triggered host work, so retain one timer witness until
		 * a paint actually consumes it.  System GL update() can be coalesced
		 * while a nested command/test event loop is active; retiring the
		 * idle tail here used to strand both that request and a refinement
		 * barrier waiting for the following completed frame. */
		qgcanvas_queue_obol_progressive_update(s, w);
	    } else
		qgcanvas_request_lod_overlay_repaint(s, w, true);
	    return;
	}

	/* Poll background providers and LoD service state without forcing an
	 * expensive duplicate paint.  advanceProgressiveWork() requests a
	 * render only when presentation data or a retained PoP cut actually
	 * changes; otherwise keep this lightweight timer pump alive. */
	if (work.pumpPending()) {
	    const std::chrono::steady_clock::time_point pumpStarted =
		std::chrono::steady_clock::now();
	    (void)s.obol->advanceProgressiveWork(NULL, NULL);
	    const std::chrono::microseconds pumpElapsed =
		std::chrono::duration_cast<std::chrono::microseconds>(
		    std::chrono::steady_clock::now() - pumpStarted);
	    /* An expensive bounded slice is positive evidence that immediately
	     * runnable owner-thread work remains.  Cheap no-progress polling keeps
	     * the normal cadence so worker waits and time-based debounce cannot
	     * busy-spin the GUI thread. */
	    s.progressive_update_delay_msec =
		pumpElapsed.count() >= 1000 ? 1 : 16;
	} else {
	    s.progressive_update_delay_msec = 16;
	}
	work = s.obol->getHostWorkSnapshot();
	/* A bounded pump can publish the terminal convergence transition without
	 * changing geometry.  Synchronize that user-facing state here, not only
	 * after a rendered frame: otherwise the work latch is already empty and
	 * no future paint exists to remove the last LoD progress HUD. */
	const bool lodProgressPublish =
	    qgcanvas_sync_obol_lod_progress(s, false);
	if (lodProgressPublish)
	    qgcanvas_request_obol_render_if_idle(s, "lod-progress-state");
	work = s.obol->getHostWorkSnapshot();
	const bool needsPresentation = work.renderPending() ||
	    lodProgressPublish;
	if (needsPresentation) {
	    /* QWidget::update() may remain coalesced indefinitely while the GUI
	     * test/command layer is running a nested event loop and no unrelated
	     * expose event occurs.  A direct-GL canvas has Qt's swap machinery as
	     * an additional wake source; the software canvas does not.  Its
	     * bounded endpoint pump therefore presents synchronously from
	     * this timer callback.  The callback is outside paintEvent, and each
	     * traversal retains Obol's hard frame deadline, so this cannot recurse
	     * or turn an expensive frame into an uninterruptible GUI loop. */
	    /* A render-only terminal transition has no remaining pump wake.  In a
	     * nested command/test loop, update() may remain coalesced after the
	     * producer has gone idle, leaving the final HUD/presentation request
	     * latched indefinitely.  Present that narrow state synchronously; an
	     * active progressive pass remains queued and bounded as before. */
	    const bool terminalPresentation = work.renderPending() &&
		!work.pumpPending();
	    if (s.software_backend ||
		s.obol->hasPendingLodRefinementFrame() || terminalPresentation) {
		w->repaint();
	    } else
		w->update();
	} else
	    qgcanvas_request_lod_overlay_repaint(s, w);
	/* A synchronous software paint may consume the last render request and
	 * publish IDLE before returning here.  Arm the no-work tail from the
	 * post-presentation snapshot, not the stale pre-paint work record. */
	work = s.obol->getHostWorkSnapshot();
	if (work.flags == BOBOL_HOST_WORK_NONE &&
	    s.lod_progress_last_state.visible)
	    s.lod_progress_idle_tail_pending = true;
	/* Do not depend on Qt delivering that paint to keep either the provider
	 * pump or an explicit frame request alive.  update() may be coalesced with
	 * the paint whose completion raised a calibration request.  That request
	 * does not necessarily represent background progressive work, so a timer
	 * gated only by hasProgressiveWorkPending() can strand it forever.  The
	 * next timer remains lightweight and retires itself as soon as a paint
	 * consumes the request and no provider work remains. */
	qgcanvas_queue_obol_progressive_update(s, w);
    });
	s.progressive_update_connected = true;
    }
    s.progressive_update_timer.start(delay);
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
				 bool *completedPresentation = nullptr,
				 std::optional<uint64_t> *imageFeatureRevision = nullptr)
{
    img = QImage();
    if (completedPresentation)
	*completedPresentation = false;
    if (imageFeatureRevision)
	imageFeatureRevision->reset();
    if (!s.obol || !s.obol->getViewport())
	return;

    qgcanvas_sync_obol_viewport(s, w);
    /* Only the endpoint's actual paint may drive the progressive state
     * machine.  A diagnostic/export traversal must be observational: making
     * a checkpoint advance admission can reopen a budget probe after the
     * scripted idle barrier and leave the supposedly final report pending. */
    if (recordPresentationTiming)
	qgcanvas_advance_obol_progressive(s);
    /* The progressive pump may have crossed into or out of convergence since
     * the last endpoint paint.  Synchronize that transition before traversing
     * the retained faceplate so the first frame in the new state does not
     * present the previous progress indicator and merely queue a corrective
     * second frame. */
    if (recordPresentationTiming)
	(void)qgcanvas_sync_obol_lod_progress(s,
	    s.obol->isRenderRequested() != FALSE);
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
    SbBool lodPlanningRelevant = TRUE;
    if (recordPresentationTiming && consumeRenderRequest)
	(void)s.obol->consumeRenderRequest(NULL, &lodCapacityRelevant,
	    &lodPlanningRelevant);
    if (!recordPresentationTiming || !consumeRenderRequest) {
	lodCapacityRelevant = FALSE;
	lodPlanningRelevant = FALSE;
    } else if (!s.obol->isLodPresentationCapacityRelevant()) {
	lodCapacityRelevant = FALSE;
    }

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
    if (recordPresentationTiming && presentationState)
	presentationState->beginCadPresentationFrame();
    const uint64_t traversedFeatureRevision = s.obol->features().presentationRevision();
    const SbBool rendered = renderer.render(s.obol->getRenderRoot());
    const uint64_t completed = s.obol->beginRenderTiming();
    const uint64_t cadExecutionAfter = presentationState ?
	presentationState->cadPresentationExecutionSerial() : 0;
    if (recordPresentationTiming && presentationState)
	presentationState->refreshCadPresentationFrameStatus();
    const BObolCadPreparationProgress cadPreparation = presentationState ?
	presentationState->cadPresentationPreparationProgress() :
	BOBOL_CAD_PREPARATION_NONE;
    const SbBool cadFrameIncomplete = presentationState &&
	!presentationState->lastCadPresentationFrameExact();
    /* The action callback is sampled between Coin nodes.  A callback node or
     * one large native draw can return after its last sample.  Such an exact
     * late frame is coherent and publishable; its actual duration lets the
     * ordinary capacity reducer classify the result. */
    const SbBool traversalInterrupted = !rendered ||
	(action && action->hasTerminated()) || cadFrameIncomplete;
    if (action && deadlineDuration)
	action->setAbortCallback(
	    deadlineContext.previous, deadlineContext.previousData);
    const BObolPresentationTimingContext timingContext(
	lodCapacityRelevant ?
	    BObolLodCapacityRelevance::RELEVANT :
	    BObolLodCapacityRelevance::EXCLUDED,
	lodPlanningRelevant ?
	    BObolLodPlanningRelevance::RELEVANT :
	    BObolLodPlanningRelevance::EXCLUDED,
	cadExecutionAfter != cadExecutionBefore ?
	    BObolCadPresentationExecution::EXECUTED :
	    BObolCadPresentationExecution::NOT_EXECUTED,
	cadPreparation,
	cadFrameIncomplete ?
	    BObolCadPresentationCompleteness::INCOMPLETE :
	    BObolCadPresentationCompleteness::EXACT);
    /* Image export/checkpoint readback is a second traversal, not a frame the
     * viewport needed to present.  Feeding it into the scene LoD capacity
     * estimator makes screenshot frequency and PNG test checkpoints alter
     * the terminal PoP cut.  QgSW's paint path opts in below; diagnostic
     * image producers deliberately do not. */
    const SbBool coherentPresentation = recordPresentationTiming ?
	s.obol->finishPresentationRenderTiming(started,
	    completed > started ? completed - started : 1,
	    traversalInterrupted, timingContext) :
	!traversalInterrupted;
    const SbBool interrupted = !coherentPresentation;
    if (interrupted) {
	/* A retained frame is useful only for the same physical viewport.  In
	 * particular, never present pre-resize pixels while an interrupted
	 * traversal is preparing the first frame at the new dimensions. */
	const bool retainedFrameMatches =
	    qgcanvas_has_completed_software_frame(s, w) &&
	    s.last_completed_software_frame.width() == size[0] &&
	    s.last_completed_software_frame.height() == size[1];
	if (retainedFrameMatches) {
	    if (imageFeatureRevision)
		*imageFeatureRevision = s.completed_feature_revision;
	    /* The paint path requests renderer-native bottom-up pixels and applies
	     * its inverse QPainter transform.  Export/readback callers request the
	     * normal top-down Qt image contract.  Keeping this distinction here
	     * avoids adding a full-frame copy to every OSMesa presentation. */
	    img = borrowRendererBuffer ?
		s.last_completed_software_frame :
		qgcanvas_flip_vertical(s.last_completed_software_frame);
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
	s.completed_feature_revision = traversedFeatureRevision;
	/* Store the endpoint's logical-pixel contract on the retained image,
	 * rather than only on the temporary return value below.  QImage is
	 * implicitly shared; setting DPR on img after copying it from this member
	 * detaches the metadata and leaves idle QgSW repaints at DPR 1. */
	if (w)
	    s.last_completed_software_frame.setDevicePixelRatio(
		w->devicePixelRatioF());
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
    if (imageFeatureRevision)
	*imageFeatureRevision = traversedFeatureRevision;
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
    s.faceplate_sync_initialized = false;
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
    s.progressive_update_timer.stop();
    delete s.offscreen_renderer;
    s.offscreen_renderer = nullptr;
    delete s.presentation_fbo;
    s.presentation_fbo = nullptr;
    delete s.presentation_staging_fbo;
    s.presentation_staging_fbo = nullptr;
    s.last_completed_software_frame = QImage();
    s.last_presented_software_frame = QImage();
    s.lod_progress_overlay = QgLodProgressOverlayState();
    s.lod_progress_presentation.reset();
    s.lod_progress_overlay_dirty = false;
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
    s.faceplate_sync_initialized = false;
    s.lod_progress_last_state = BObolLodProgressDisplayStatus();
    s.lod_progress_last_publish =
	std::chrono::steady_clock::time_point();
    s.lod_progress_overlay = QgLodProgressOverlayState();
    s.lod_progress_presentation.reset();
    s.lod_progress_overlay_dirty = true;
    s.lod_progress_overlay_last_request =
	std::chrono::steady_clock::time_point();
    /* A completed image belongs to its controller/scene identity.  Never
     * preserve pixels from the previous endpoint through a borrowed-controller
     * replacement. */
    s.presentation_fbo_has_completed_frame = false;
    s.last_completed_software_frame = QImage();
    s.last_presented_software_frame = QImage();
    s.presented_feature_revision.reset();

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
    s.obol->requestLodCapacityRender("qt-controller-bind");
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

template <typename NodeType>
static inline int
qgcanvas_find_obol_faceplate_child(SoGroup *group, const char *overlayId)
{
    if (!group || !overlayId)
	return -1;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node || !node->isOfType(NodeType::getClassTypeId()))
	    continue;
	NodeType *faceplateNode = static_cast<NodeType *>(node);
	if (bu_strcmp(faceplateNode->overlayId.getValue().getString(),
		overlayId) == 0)
	    return i;
    }
    return -1;
}

static inline void
qgcanvas_publish_obol_faceplate_child(QgCanvasState &s, SoGroup *group,
	int childIndex, SoNode *candidate)
{
    if (!group || (childIndex < 0 && !candidate))
	return;

    const int childCount = group->getNumChildren();
    const size_t appendedChildCount = childIndex < 0 && candidate ? 1u : 0u;
    std::vector<SoNode *> children;
    children.reserve(static_cast<size_t>(childCount) + appendedChildCount);
    for (int i = 0; i < childCount; i++) {
	if (i == childIndex) {
	    if (candidate)
		children.push_back(candidate);
	    continue;
	}
	children.push_back(group->getChild(i));
    }
    if (childIndex < 0)
	children.push_back(candidate);

    /* The no-GED canvas fallback has no feature-store record to stage. Build
     * its complete node off graph, then use the same prepared child-list
     * primitive as the retained publishers so observers see one root cut. */
    auto replacement = group->getChildren()->prepareReplacement(children);
    replacement->commit();

    std::exception_ptr failure;
    try {
	replacement->notify();
    } catch (...) {
	failure = std::current_exception();
    }
    try {
	qgcanvas_request_obol_render_if_idle(s, "faceplate");
    } catch (...) {
	if (!failure)
	    failure = std::current_exception();
    }
    if (failure)
	std::rethrow_exception(failure);
}

static inline void
qgcanvas_sync_obol_axes(QgCanvasState &s,
			SoGroup *group,
			const char *overlayId,
			const struct bv_axes_state &state)
{
    const int childIndex =
	qgcanvas_find_obol_faceplate_child<SoBRLAxes>(group, overlayId);
    if (!state.draw) {
	qgcanvas_publish_obol_faceplate_child(s, group, childIndex, NULL);
	return;
    }

    SbModernUtils::SoNodeRef axesOwner(new SoBRLAxes);
    SoBRLAxes *axes = static_cast<SoBRLAxes *>(axesOwner.get());
    axes->overlayId = overlayId;
    if (!bobol_axes_configure_from_view(axes, &state))
	return;
    qgcanvas_publish_obol_faceplate_child(s, group, childIndex, axes);
}

static inline void
qgcanvas_sync_obol_grid(QgCanvasState &s,
			SoGroup *group,
			const struct bv_grid_state &state,
			struct ged_view_context *view_ctx)
{
    const char *overlayId = "faceplate::grid";
    const int childIndex =
	qgcanvas_find_obol_faceplate_child<SoBRLGrid>(group, overlayId);
    if (!state.draw && !state.snap) {
	qgcanvas_publish_obol_faceplate_child(s, group, childIndex, NULL);
	return;
    }

    SbModernUtils::SoNodeRef gridOwner(new SoBRLGrid);
    SoBRLGrid *grid = static_cast<SoBRLGrid *>(gridOwner.get());
    grid->overlayId = overlayId;
    if (!bobol_grid_configure_from_view_context(grid, &state, view_ctx))
	return;
    qgcanvas_publish_obol_faceplate_child(s, group, childIndex, grid);
}

static inline void
qgcanvas_sync_obol_adc(QgCanvasState &s,
		       SoGroup *group,
		       const struct bv_adc_state &state)
{
    const char *overlayId = "faceplate::adc";
    const int childIndex =
	qgcanvas_find_obol_faceplate_child<SoBRLADC>(group, overlayId);
    if (!state.draw) {
	qgcanvas_publish_obol_faceplate_child(s, group, childIndex, NULL);
	return;
    }

    SbModernUtils::SoNodeRef adcOwner(new SoBRLADC);
    SoBRLADC *adc = static_cast<SoBRLADC *>(adcOwner.get());
    adc->overlayId = overlayId;
    if (!bobol_adc_configure_from_view(adc, &state))
	return;
    qgcanvas_publish_obol_faceplate_child(s, group, childIndex, adc);
}

static inline void
qgcanvas_sync_obol_faceplate(QgCanvasState &s)
{
    if (!s.obol || !s.v)
	return;

    const struct bv *view = bv_context_view_const(s.v);
    const unsigned long long viewHash = bv_hash(view);
    const uint64_t frameRevision = bv_frame_revision_get(view);
    if (s.faceplate_sync_initialized &&
	s.faceplate_view_hash == viewHash &&
	s.faceplate_frame_revision == frameRevision)
	return;

    struct ged_view_context *view_ctx = ged_view_context_from_bv(s.v);
    struct ged *gedp = ged_view_context_owner(view_ctx);
    if (gedp) {
	(void)ged_view_faceplate_sync(gedp, view_ctx);
	s.faceplate_view_hash = viewHash;
	s.faceplate_frame_revision = frameRevision;
	s.faceplate_sync_initialized = true;
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
    (void)bv_grid_state_get(&grid, view);
    (void)bv_model_axes_state_get(&modelAxes, view);
    (void)bv_view_axes_state_get(&viewAxes, view);
    (void)bv_adc_state_get(&adc, view);

    qgcanvas_sync_obol_grid(s, group, grid, view_ctx);
    qgcanvas_sync_obol_axes(s, group, "faceplate::model_axes", modelAxes);
    qgcanvas_sync_obol_axes(s, group, "faceplate::view_axes", viewAxes);
    qgcanvas_sync_obol_adc(s, group, adc);
    s.faceplate_view_hash = viewHash;
    s.faceplate_frame_revision = frameRevision;
    s.faceplate_sync_initialized = true;
}

static inline QString
qgcanvas_lod_overlay_count(uint64_t value)
{
    if (value >= 1000000)
	return QStringLiteral("%1M").arg(
	    static_cast<double>(value) / 1000000.0, 0, 'f', 1);
    if (value >= 1000)
	return QStringLiteral("%1k").arg(
	    static_cast<double>(value) / 1000.0, 0, 'f', 1);
    return QString::number(static_cast<qulonglong>(value));
}

static inline QString
qgcanvas_lod_overlay_duration(uint64_t microseconds)
{
    if (microseconds >= 1000000)
	return QStringLiteral("%1 s").arg(
	    static_cast<double>(microseconds) / 1000000.0, 0, 'f', 1);
    if (microseconds >= 1000)
	return QStringLiteral("%1 ms").arg(
	    static_cast<double>(microseconds) / 1000.0, 0, 'f',
	    microseconds < 10000 ? 1 : 0);
    return QStringLiteral("<1 ms");
}

static inline QString
qgcanvas_lod_overlay_bytes(uint64_t bytes)
{
    static constexpr double kibibyte = 1024.0;
    static constexpr double mebibyte = 1024.0 * kibibyte;
    static constexpr double gibibyte = 1024.0 * mebibyte;
    if (bytes >= static_cast<uint64_t>(gibibyte))
	return QStringLiteral("%1 GiB").arg(
	    static_cast<double>(bytes) / gibibyte, 0, 'f', 1);
    if (bytes >= static_cast<uint64_t>(mebibyte))
	return QStringLiteral("%1 MiB").arg(
	    static_cast<double>(bytes) / mebibyte, 0, 'f', 1);
    if (bytes >= static_cast<uint64_t>(kibibyte))
	return QStringLiteral("%1 KiB").arg(
	    static_cast<double>(bytes) / kibibyte, 0, 'f', 1);
    return QStringLiteral("%1 B").arg(static_cast<qulonglong>(bytes));
}

static inline QString
qgcanvas_lod_producer_stage_short_title(int stage)
{
    switch (stage) {
	case BOBOL_LOD_PRODUCER_STAGE_ASSET_SERIALIZATION:
	    return QStringLiteral("waiting for shared asset");
	case BOBOL_LOD_PRODUCER_STAGE_CACHE_LOOKUP:
	    return QStringLiteral("checking cache");
	case BOBOL_LOD_PRODUCER_STAGE_SOURCE_PREPARATION:
	    return QStringLiteral("preparing source");
	case BOBOL_LOD_PRODUCER_STAGE_COVERAGE_PREVIEW:
	    return QStringLiteral("sampling coverage");
	case BOBOL_LOD_PRODUCER_STAGE_SOURCE_HASHING:
	    return QStringLiteral("hashing");
	case BOBOL_LOD_PRODUCER_STAGE_BOUNDS_ANALYSIS:
	    return QStringLiteral("analyzing bounds");
	case BOBOL_LOD_PRODUCER_STAGE_FACE_CLASSIFICATION:
	    return QStringLiteral("classifying");
	case BOBOL_LOD_PRODUCER_STAGE_PREFIX_MATERIALIZATION:
	    return QStringLiteral("building meshes");
	case BOBOL_LOD_PRODUCER_STAGE_SPATIAL_CONSTRUCTION:
	    return QStringLiteral("building pages");
	case BOBOL_LOD_PRODUCER_STAGE_CACHE_PERSISTENCE:
	    return QStringLiteral("saving cache");
	default:
	    return QString();
    }
}

static inline QgLodProgressOverlayState
qgcanvas_lod_progress_overlay_state(
    const BObolLodConvergenceStatus &status)
{
    QgLodProgressOverlayState result;
    const BObolLodProgressDisplayStatus display =
	status.progressDisplayStatus();
    result.visible = display.visible != FALSE;
    if (!result.visible)
	return result;

    result.terminalReady = display.terminalReady != FALSE;
    result.resourceLimited = status.memoryLimited != FALSE ||
	status.gpuMemoryPressure != FALSE;
    result.displayClass = display.publicationClass;
    const float fraction = status.progressEstimateAvailable ?
	status.estimatedFraction : status.fraction;
    result.percent = static_cast<int>(std::floor(
	std::max(0.0f, std::min(1.0f, fraction)) * 100.0f + 0.5f));
    result.determinate = status.progressEstimateAvailable != FALSE;
    result.refinementCycleBased =
	status.progressEstimateRefinementCycleBased != FALSE;
    result.estimatedRemainingMilliseconds =
	status.estimatedRemainingMilliseconds;
    result.remainingRefinementCycles =
	status.estimatedRemainingRefinementCycles;
    static constexpr uint64_t animationIntervalMilliseconds = 100;
    static constexpr unsigned int animationStepCount = 20;
    result.animationStep = static_cast<unsigned int>(
	(status.episode.elapsedMilliseconds / animationIntervalMilliseconds) %
	animationStepCount);
    if (status.terminalError || status.failedSourceCount > 0) {
	result.title = QStringLiteral("Geometry preparation failed");
    } else if (status.terminal &&
	(status.memoryLimited || status.gpuMemoryPressure)) {
	result.title = QStringLiteral("Detail limited by memory");
    } else if (status.terminal && status.performanceLimited) {
	result.title = QStringLiteral("Detail limited by frame budget");
    } else if (status.episode.stableViewReached && !status.terminal) {
	result.title = QStringLiteral("Finalizing view");
    } else if (!status.episode.firstMeshReached && !status.terminal) {
	result.title = QStringLiteral("Preparing geometry");
	result.displayClass = BOBOL_LOD_PROGRESS_DISPLAY_PREPARING;
    } else if (display.publicationClass ==
	BOBOL_LOD_PROGRESS_DISPLAY_INTERACTIVE) {
	result.title = QStringLiteral("Adjusting view detail");
	result.displayClass = BOBOL_LOD_PROGRESS_DISPLAY_INTERACTIVE;
    } else if (status.sourcePreparationPending ||
	status.activeProducerCount > 0) {
	result.title = QStringLiteral("Preparing geometry");
	result.displayClass = BOBOL_LOD_PROGRESS_DISPLAY_PREPARING;
    } else {
	result.title = QStringLiteral("Refining visible detail");
	result.displayClass = BOBOL_LOD_PROGRESS_DISPLAY_SETTLING;
    }

    const size_t temporaryProxyCount =
	status.proxyReasons.temporaryStructuralOccurrenceCount();
    const size_t budgetProxyCount =
	status.proxyReasons.budgetLimitedStructuralOccurrenceCount();
    const size_t failedProxyCount =
	status.proxyReasons.terminalFailureOccurrenceCount;
    const auto appendDetail = [&](const QString &text) {
	if (text.isEmpty())
	    return;
	if (!result.detail.isEmpty())
	    result.detail += QStringLiteral(" | ");
	result.detail += text;
    };
    /* Put latency first so it survives elision on a narrow canvas.  Queue
     * age means no worker has acquired that task; producer and stage time
     * mean actual source work is in progress. */
    if (status.oldestPendingTaskAgeMicroseconds > 0)
	appendDetail(QStringLiteral("oldest queued %1").arg(
	    qgcanvas_lod_overlay_duration(
		status.oldestPendingTaskAgeMicroseconds)));
    if (status.runnableQueuedTasks > 0)
	appendDetail(QStringLiteral("%1 ready for worker").arg(
	    qgcanvas_lod_overlay_count(status.runnableQueuedTasks)));
    if (status.dependencyBlockedTasks > 0)
	appendDetail(QStringLiteral("%1 waiting on prerequisites").arg(
	    qgcanvas_lod_overlay_count(status.dependencyBlockedTasks)));
    if (status.transientMemoryBlockedTasks > 0)
	appendDetail(QStringLiteral("%1 waiting on transient memory").arg(
	    qgcanvas_lod_overlay_count(status.transientMemoryBlockedTasks)));
    if (status.cpuAdmissionWaitingTasks > 0)
	appendDetail(QStringLiteral("%1 waiting for CPU").arg(
	    qgcanvas_lod_overlay_count(status.cpuAdmissionWaitingTasks)));
    if (status.taskSubmissionCapacityBlocked)
	appendDetail(QStringLiteral("producer task slots full"));
    if (status.resultSubmissionCapacityBlocked)
	appendDetail(QStringLiteral("result queue slots full"));
    if (status.sourcePreparationPending &&
	status.sourcePreparationTotalUnits >
	    status.sourcePreparationCompletedUnits) {
	const uint64_t remaining = status.sourcePreparationTotalUnits -
	    status.sourcePreparationCompletedUnits;
	appendDetail(QStringLiteral("%1 source preparation%2 remaining")
	    .arg(qgcanvas_lod_overlay_count(remaining))
	    .arg(remaining == 1 ? QString() : QStringLiteral("s")));
    }
    if (status.activeProducerCount > 0 &&
	status.producerStageElapsedMicroseconds > 0)
	appendDetail(QStringLiteral("current stage %1").arg(
	    qgcanvas_lod_overlay_duration(
		status.producerStageElapsedMicroseconds)));
    if (status.activeProducerCount > 0 &&
	status.maximumProducerElapsedMicroseconds >
	    status.producerStageElapsedMicroseconds &&
	status.maximumProducerElapsedMicroseconds -
	    status.producerStageElapsedMicroseconds > 100000)
	appendDetail(QStringLiteral("active %1 total").arg(
	    qgcanvas_lod_overlay_duration(
		status.maximumProducerElapsedMicroseconds)));
    if (status.maximumProducerQueueWaitMicroseconds >= 100000)
	appendDetail(QStringLiteral("queued %1 before start").arg(
	    qgcanvas_lod_overlay_duration(
		status.maximumProducerQueueWaitMicroseconds)));
    if (status.activeProducerSourceFaceCount > 0 ||
	status.activeProducerSourcePointCount > 0 ||
	status.activeProducerSourceByteCount > 0) {
	QString sourceSize;
	if (status.activeProducerSourceFaceCount > 0)
	    sourceSize = QStringLiteral("%1 faces").arg(
		qgcanvas_lod_overlay_count(
		    status.activeProducerSourceFaceCount));
	else if (status.activeProducerSourcePointCount > 0)
	    sourceSize = QStringLiteral("%1 points").arg(
		qgcanvas_lod_overlay_count(
		    status.activeProducerSourcePointCount));
	if (status.activeProducerSourceByteCount > 0) {
	    if (!sourceSize.isEmpty())
		sourceSize += QStringLiteral(", ");
	    sourceSize += qgcanvas_lod_overlay_bytes(
		status.activeProducerSourceByteCount);
	}
	appendDetail(QStringLiteral("source %1").arg(sourceSize));
    }
    if (status.producerStage != BOBOL_LOD_PRODUCER_STAGE_NONE) {
	QString stageDistribution;
	const int producerStageOrder[] = {
	    BOBOL_LOD_PRODUCER_STAGE_ASSET_SERIALIZATION,
	    BOBOL_LOD_PRODUCER_STAGE_CACHE_LOOKUP,
	    BOBOL_LOD_PRODUCER_STAGE_SOURCE_PREPARATION,
	    BOBOL_LOD_PRODUCER_STAGE_BOUNDS_ANALYSIS,
	    BOBOL_LOD_PRODUCER_STAGE_COVERAGE_PREVIEW,
	    BOBOL_LOD_PRODUCER_STAGE_SOURCE_HASHING,
	    BOBOL_LOD_PRODUCER_STAGE_FACE_CLASSIFICATION,
	    BOBOL_LOD_PRODUCER_STAGE_PREFIX_MATERIALIZATION,
	    BOBOL_LOD_PRODUCER_STAGE_SPATIAL_CONSTRUCTION,
	    BOBOL_LOD_PRODUCER_STAGE_CACHE_PERSISTENCE
	};
	for (int stage : producerStageOrder) {
	    const size_t count = status.producerStageTaskCounts[stage];
	    if (!count)
		continue;
	    if (!stageDistribution.isEmpty())
		stageDistribution += QStringLiteral(", ");
	    stageDistribution += QStringLiteral("%1 %2")
		.arg(static_cast<qulonglong>(count))
		.arg(qgcanvas_lod_producer_stage_short_title(stage));
	}
	if (!stageDistribution.isEmpty())
	    appendDetail(stageDistribution);
	if (status.producerStageTotalUnits > 0) {
	    const uint64_t completed = std::min(
		status.producerStageCompletedUnits,
		status.producerStageTotalUnits);
	    const int stagePercent = static_cast<int>(std::floor(
		100.0 * static_cast<double>(completed) /
		static_cast<double>(status.producerStageTotalUnits) + 0.5));
	    appendDetail(QStringLiteral("%1 %2% (%3/%4)")
		.arg(qgcanvas_lod_producer_stage_short_title(
		    status.producerStage))
		.arg(stagePercent)
		.arg(qgcanvas_lod_overlay_count(
		    completed))
		.arg(qgcanvas_lod_overlay_count(
		    status.producerStageTotalUnits)));
	}
    }
    const size_t unresolvedDetailCount = status.activePayloadCount >
	status.satisfiedPayloadCount ?
	status.activePayloadCount - status.satisfiedPayloadCount : 0;
    if (!status.terminal && unresolvedDetailCount > 0 &&
	status.renderCostBudget > 0 && status.retainedRenderCost > 0) {
	const long double percent = 100.0L * static_cast<long double>(
	    status.retainedRenderCost) /
	    static_cast<long double>(status.renderCostBudget);
	if (percent >= 95.0L) {
	    const unsigned int rounded = static_cast<unsigned int>(
		std::min<long double>(999.0L, std::floor(percent + 0.5L)));
	    appendDetail(QStringLiteral("render budget %1% allocated").arg(
		rounded));
	}
    }
    if (status.temporaryCoverageOccurrenceCount > 0)
	appendDetail(QStringLiteral("%1 temporary coverage %2").arg(
	    qgcanvas_lod_overlay_count(
		status.temporaryCoverageOccurrenceCount),
	    status.temporaryCoverageOccurrenceCount == 1 ?
		QStringLiteral("preview") : QStringLiteral("previews")));
    if (status.activeSourceFaces > 0)
	appendDetail(QStringLiteral("%1 source %2").arg(
	    qgcanvas_lod_overlay_count(status.activeSourceFaces),
	    status.activeSourceFaces == 1 ? QStringLiteral("triangle") :
		QStringLiteral("triangles")));
    else if (status.sourceMeshOccurrenceCount > 0)
	appendDetail(QStringLiteral("%1 source %2").arg(
	    qgcanvas_lod_overlay_count(status.sourceMeshOccurrenceCount),
	    status.sourceMeshOccurrenceCount == 1 ? QStringLiteral("mesh") :
		QStringLiteral("meshes")));
    if (temporaryProxyCount > 0)
	appendDetail(QStringLiteral("%1 %2").arg(
	    qgcanvas_lod_overlay_count(temporaryProxyCount),
	    temporaryProxyCount == 1 ? QStringLiteral("temporary box") :
		QStringLiteral("temporary boxes")));
    if (budgetProxyCount > 0)
	appendDetail(QStringLiteral("%1 %2").arg(
	    qgcanvas_lod_overlay_count(budgetProxyCount),
	    budgetProxyCount == 1 ? QStringLiteral("budget-limited box") :
		QStringLiteral("budget-limited boxes")));
    if (failedProxyCount > 0)
	appendDetail(QStringLiteral("%1 failed %2").arg(
	    qgcanvas_lod_overlay_count(failedProxyCount),
	    failedProxyCount == 1 ? QStringLiteral("box") :
		QStringLiteral("boxes")));
    const size_t subpixelProxyCount =
	status.proxyReasons.intentionalSubpixelOccurrenceCount;
    if (subpixelProxyCount > 0)
	appendDetail(QStringLiteral("%1 subpixel %2").arg(
	    qgcanvas_lod_overlay_count(subpixelProxyCount),
	    subpixelProxyCount == 1 ? QStringLiteral("point") :
		QStringLiteral("points")));
    if (!status.episode.firstMeshReached && !status.terminal)
	appendDetail(QStringLiteral("first mesh pending"));
    if (status.failedSourceCount > 0)
	appendDetail(QStringLiteral("%1 failed source%2")
	    .arg(static_cast<qulonglong>(status.failedSourceCount))
	    .arg(status.failedSourceCount == 1 ? QString() :
		QStringLiteral("s")));
    if (status.terminal &&
	(status.memoryLimited || status.gpuMemoryPressure ||
	 status.performanceLimited))
	appendDetail(QStringLiteral("best available under current budget"));
    if (status.episode.elapsedMilliseconds >= 1000)
	appendDetail(QStringLiteral("%1 s elapsed").arg(
	    static_cast<double>(status.episode.elapsedMilliseconds) / 1000.0,
	    0, 'f', 1));
    return result;
}

static inline QString
qgcanvas_lod_eta_text(uint64_t milliseconds)
{
    static constexpr uint64_t millisecondsPerSecond = 1000;
    static constexpr uint64_t secondsPerMinute = 60;
    const uint64_t seconds = milliseconds / millisecondsPerSecond +
	(milliseconds % millisecondsPerSecond != 0 ? 1 : 0);
    if (seconds <= 5)
	return QStringLiteral("under 5 s remaining");
    if (seconds < secondsPerMinute) {
	static constexpr uint64_t secondsPerBucket = 10;
	const uint64_t rounded = ((seconds + secondsPerBucket / 2) /
	    secondsPerBucket) * secondsPerBucket;
	return QStringLiteral("about %1 s remaining").arg(
	    static_cast<qulonglong>(std::max<uint64_t>(secondsPerBucket,
		rounded)));
    }
    const uint64_t minutes = (seconds + secondsPerMinute / 2) /
	secondsPerMinute;
    return QStringLiteral("about %1 min remaining").arg(
	static_cast<qulonglong>(std::max<uint64_t>(1, minutes)));
}

/* An ETA is intentionally unavailable while the visible frontier grows or a
 * completed frame fails to reduce it.  Preserve the useful exact facts in
 * that state instead of replacing them with an uninformative "estimating"
 * message.  These are current work counts, not a prediction, so they remain
 * truthful even when budget rebalancing makes the count non-monotonic. */
static inline QString
qgcanvas_lod_unestimated_work_text(
    const BObolLodConvergenceStatus &status)
{
    if (status.sourcePreparationPending &&
	status.sourcePreparationTotalUnits >
	    status.sourcePreparationCompletedUnits) {
	const uint64_t remaining = status.sourcePreparationTotalUnits -
	    status.sourcePreparationCompletedUnits;
	return QStringLiteral("%1 source preparation%2 remaining")
	    .arg(qgcanvas_lod_overlay_count(remaining))
	    .arg(remaining == 1 ? QString() : QStringLiteral("s"));
    }

    if (status.rendererPreparationRemainingUnits > 0)
	return QStringLiteral("%1 renderer work units remaining").arg(
	    qgcanvas_lod_overlay_count(
		status.rendererPreparationRemainingUnits));

    const size_t unresolved = status.activePayloadCount >
	status.satisfiedPayloadCount ?
	status.activePayloadCount - status.satisfiedPayloadCount : 0;
    if (unresolved > 0)
	return QStringLiteral("%1 visible item%2 still refining")
	    .arg(qgcanvas_lod_overlay_count(unresolved))
	    .arg(unresolved == 1 ? QString() : QStringLiteral("s"));

    if (status.inFlight > 0)
	return QStringLiteral("%1 geometry task%2 running")
	    .arg(qgcanvas_lod_overlay_count(status.inFlight))
	    .arg(status.inFlight == 1 ? QString() : QStringLiteral("s"));
    if (status.queuedResults > 0)
	return QStringLiteral("%1 completed result%2 awaiting publication")
	    .arg(qgcanvas_lod_overlay_count(status.queuedResults))
	    .arg(status.queuedResults == 1 ? QString() : QStringLiteral("s"));
    if (status.pendingTasks > 0)
	return QStringLiteral("%1 geometry task%2 queued")
	    .arg(qgcanvas_lod_overlay_count(status.pendingTasks))
	    .arg(status.pendingTasks == 1 ? QString() : QStringLiteral("s"));
    if (status.queuedCacheWrites > 0)
	return QStringLiteral("%1 cache write%2 pending")
	    .arg(qgcanvas_lod_overlay_count(status.queuedCacheWrites))
	    .arg(status.queuedCacheWrites == 1 ? QString() :
		QStringLiteral("s"));
    if (status.budgetCalibrationPending ||
	status.pointProxyCalibrationPending ||
	status.stablePointProxyCalibrationPending)
	return QStringLiteral("Measuring the render budget");
    if (status.residentGrowthReallocationPending)
	return QStringLiteral("Rebalancing newly available detail");
    if (status.stablePresentationHandoffPending)
	return QStringLiteral("Reconciling the final presentation");
    if (status.publicationFramePending)
	return QStringLiteral("Publishing refined geometry");
    if (status.refinementFramePending)
	return QStringLiteral("Waiting for a refinement frame");
    return QString();
}

static inline QgLodProgressOverlayState
qgcanvas_stabilize_lod_progress_overlay(
    QgLodProgressPresentationState &presentation,
    const BObolLodConvergenceStatus &status,
    QgLodProgressOverlayState overlay)
{
    const uint64_t elapsed = status.episode.elapsedMilliseconds;
    const bool newEpisode = !presentation.active ||
	presentation.episodeRevision != status.episodeRevision ||
	elapsed < presentation.lastElapsedMilliseconds;
    if (!overlay.visible) {
	presentation.reset();
	return overlay;
    }
    if (newEpisode) {
	presentation.reset();
	presentation.active = true;
	presentation.episodeRevision = status.episodeRevision;
    }
    presentation.lastElapsedMilliseconds = elapsed;
    const bool wasDeterminate = presentation.determinateLatched;

    /* The estimator is exact about whether all work ranks are measurable, but
     * that bit can legitimately drop for a few samples while control moves
     * between otherwise continuous ranks.  Once an episode has a measured
     * fraction, retain its monotonic floor and reserve the final five percent
     * for exact presentation/certificate reconciliation. */
    if (overlay.terminalReady) {
	overlay.determinate = true;
	overlay.percent = 100;
	presentation.determinateLatched = true;
	presentation.percentFloor = 100;
    } else if (status.episode.stableViewReached) {
	/* The first visually stable frame is an exact milestone even when the
	 * remaining certificate/presentation work has no measurable rank.  Move
	 * into the reserved final five percent instead of leaving a low or
	 * indeterminate bar beside the "Finalizing view" title. */
	static constexpr int minimumFinalizingPercent = 95;
	static constexpr int maximumFinalizingPercent = 99;
	const int candidate = overlay.determinate ? overlay.percent :
	    minimumFinalizingPercent;
	presentation.percentFloor = std::min(maximumFinalizingPercent,
	    std::max(presentation.percentFloor,
		std::max(minimumFinalizingPercent, candidate)));
	overlay.determinate = true;
	overlay.percent = presentation.percentFloor;
	presentation.determinateLatched = true;
    } else if (overlay.determinate) {
	static constexpr int maximumUnstablePercent = 95;
	overlay.percent = std::min(maximumUnstablePercent, overlay.percent);
	presentation.percentFloor = std::max(
	    presentation.percentFloor, overlay.percent);
	overlay.percent = presentation.percentFloor;
	presentation.determinateLatched = true;
    } else {
	/* An unqualified cycle forecast means the frontier is still changing or
	 * completed frames have stopped improving.  Showing the predecessor's
	 * percentage in that state looks authoritative while conveying no useful
	 * progress, so return to honest indeterminate activity. */
	presentation.determinateLatched = false;
	presentation.percentFloor = 0;
    }
    if (overlay.determinate)
	overlay.animationStep = 0;

    /* Estimate a completion timestamp rather than smoothing the raw remaining
     * duration.  A coherent estimate then counts down naturally.  Large target
     * changes reset confidence instead of displaying an implausible jump from
     * "under a second" to many seconds. */
    static constexpr uint64_t etaMinimumConfidenceMilliseconds = 1500;
    static constexpr uint64_t etaMinimumToleranceMilliseconds = 2000;
    const bool etaSampleUsable = status.progressEstimateAvailable &&
	status.estimatedRemainingMilliseconds > 0 &&
	!status.terminal;
    if (etaSampleUsable) {
	const uint64_t remaining = status.estimatedRemainingMilliseconds;
	const uint64_t candidateCompletion = remaining > UINT64_MAX - elapsed ?
	    UINT64_MAX : elapsed + remaining;
	if (!presentation.smoothedCompletionMilliseconds) {
	    presentation.smoothedCompletionMilliseconds = candidateCompletion;
	    presentation.etaConfidenceStartMilliseconds = elapsed;
	    presentation.consistentEtaSamples = 1;
	    presentation.etaVisible = false;
	} else {
	    const uint64_t priorCompletion =
		presentation.smoothedCompletionMilliseconds;
	    const uint64_t difference = candidateCompletion > priorCompletion ?
		candidateCompletion - priorCompletion :
		priorCompletion - candidateCompletion;
	    const uint64_t priorRemaining = priorCompletion > elapsed ?
		priorCompletion - elapsed : 0;
	    const uint64_t tolerance = std::max(
		etaMinimumToleranceMilliseconds, priorRemaining / 3);
	    if (difference > tolerance) {
		presentation.smoothedCompletionMilliseconds =
		    candidateCompletion;
		presentation.etaConfidenceStartMilliseconds = elapsed;
		presentation.consistentEtaSamples = 1;
		presentation.etaVisible = false;
	    } else {
		if (candidateCompletion >= priorCompletion)
		    presentation.smoothedCompletionMilliseconds =
			priorCompletion +
			(candidateCompletion - priorCompletion) / 4;
		else
		    presentation.smoothedCompletionMilliseconds =
			priorCompletion -
			(priorCompletion - candidateCompletion) / 4;
		if (presentation.consistentEtaSamples < UINT_MAX)
		    presentation.consistentEtaSamples++;
		presentation.etaVisible =
		    presentation.consistentEtaSamples >= 3 &&
		    elapsed >= presentation.etaConfidenceStartMilliseconds &&
		    elapsed - presentation.etaConfidenceStartMilliseconds >=
			etaMinimumConfidenceMilliseconds;
	    }
	}
    } else {
	presentation.etaVisible = false;
	presentation.smoothedCompletionMilliseconds = 0;
	presentation.consistentEtaSamples = 0;
    }

    overlay.etaVisible = presentation.etaVisible &&
	presentation.smoothedCompletionMilliseconds > elapsed;
    overlay.estimatedRemainingMilliseconds = overlay.etaVisible ?
	presentation.smoothedCompletionMilliseconds - elapsed : 0;

    const QString unestimatedWork =
	qgcanvas_lod_unestimated_work_text(status);
    const QString pendingEstimate = status.progressEstimateAvailable ?
	QStringLiteral("estimating remaining time") : unestimatedWork;
    QString progressDetail;
    if (overlay.terminalReady) {
	progressDetail = QStringLiteral("Stable");
    } else if (overlay.determinate) {
	progressDetail = QStringLiteral("%1%").arg(overlay.percent);
	if (status.episode.stableViewReached || overlay.percent >= 95) {
	    progressDetail += QStringLiteral(" | finalizing");
	    if (overlay.etaVisible)
		progressDetail += overlay.refinementCycleBased &&
		    overlay.remainingRefinementCycles > 0 ?
		    QStringLiteral(" | ~%1 refinement cycle%2, %3")
			.arg(static_cast<qulonglong>(
			    overlay.remainingRefinementCycles))
			.arg(overlay.remainingRefinementCycles == 1 ?
			    QString() : QStringLiteral("s"))
			.arg(qgcanvas_lod_eta_text(
			    overlay.estimatedRemainingMilliseconds)) :
		    QStringLiteral(" | ") + qgcanvas_lod_eta_text(
			overlay.estimatedRemainingMilliseconds);
	    else
		progressDetail += QStringLiteral(" | ") + pendingEstimate;
	} else if (overlay.etaVisible)
	    progressDetail += overlay.refinementCycleBased &&
		overlay.remainingRefinementCycles > 0 ?
		QStringLiteral(" | ~%1 refinement cycle%2, %3")
		    .arg(static_cast<qulonglong>(
			overlay.remainingRefinementCycles))
		    .arg(overlay.remainingRefinementCycles == 1 ?
			QString() : QStringLiteral("s"))
		    .arg(qgcanvas_lod_eta_text(
			overlay.estimatedRemainingMilliseconds)) :
		QStringLiteral(" | ") + qgcanvas_lod_eta_text(
		    overlay.estimatedRemainingMilliseconds);
	else
	    progressDetail += QStringLiteral(" | ") + pendingEstimate;
    } else {
	progressDetail = unestimatedWork;
    }
    if (!overlay.detail.isEmpty()) {
	if (!progressDetail.isEmpty())
	    progressDetail += QStringLiteral(" | ");
	progressDetail += overlay.detail;
    }

    /* Keep volatile queue/stage diagnostics useful without changing the text
     * on every 100 ms animation sample.  Terminal and episode transitions are
     * published immediately. */
    static constexpr uint64_t detailRefreshMilliseconds = 750;
    const bool finalizingChanged = presentation.finalizing !=
	(status.episode.stableViewReached != FALSE);
    const bool refreshDetail = newEpisode ||
	presentation.stableDetail.isEmpty() || overlay.terminalReady ||
	status.terminalError || finalizingChanged ||
	(etaSampleUsable && presentation.consistentEtaSamples == 1) ||
	wasDeterminate != presentation.determinateLatched ||
	elapsed < presentation.stableDetailUpdateMilliseconds ||
	elapsed - presentation.stableDetailUpdateMilliseconds >=
	    detailRefreshMilliseconds;
    if (refreshDetail) {
	presentation.stableDetail = progressDetail;
	presentation.stableDetailUpdateMilliseconds = elapsed;
    }
    presentation.finalizing = status.episode.stableViewReached != FALSE;
    overlay.detail = presentation.stableDetail;
    return overlay;
}

static inline bool
qgcanvas_native_lod_progress_selected(const QgCanvasState &s)
{
    if (!s.v)
	return true;
    struct ged_view_context *view_ctx = ged_view_context_from_bv(s.v);
    if (!ged_view_context_owner(view_ctx))
	return true;
    enum ged_view_lod_progress_presentation_mode mode =
	GED_VIEW_LOD_PROGRESS_PRESENTATION_RETAINED;
    return ged_view_lod_progress_presentation_mode_get(&mode, view_ctx) &&
	mode == GED_VIEW_LOD_PROGRESS_PRESENTATION_NATIVE_HOST;
}

static inline void
qgcanvas_request_lod_overlay_repaint(QgCanvasState &s, QWidget *w,
	bool force)
{
    if (!w || !s.lod_progress_overlay_dirty)
	return;
    const std::chrono::steady_clock::time_point now =
	std::chrono::steady_clock::now();
    static constexpr std::chrono::milliseconds minimumInterval(100);
    const bool urgent = force || !s.lod_progress_overlay.visible;
    if (!urgent &&
	s.lod_progress_overlay_last_request.time_since_epoch().count() != 0 &&
	now - s.lod_progress_overlay_last_request < minimumInterval)
	return;
    s.lod_progress_overlay_last_request = now;
    w->update();
}

static inline QColor
qgcanvas_lod_overlay_color(const QgLodProgressOverlayState &overlay)
{
    /* A frame-budget-limited view is complete and usable even though its
     * terminal notice remains visible.  Give that state the same green
     * completion cue as the retained HUD so it cannot look like frozen
     * orange refinement.  Resource limits remain warnings. */
    if (overlay.terminalReady && !overlay.resourceLimited &&
	overlay.displayClass == BOBOL_LOD_PROGRESS_DISPLAY_IDLE)
	return QColor(112, 235, 135);

    switch (overlay.displayClass) {
	case BOBOL_LOD_PROGRESS_DISPLAY_PREPARING:
	case BOBOL_LOD_PROGRESS_DISPLAY_DISCOVERING:
	    return QColor(96, 190, 255);
	case BOBOL_LOD_PROGRESS_DISPLAY_ERROR:
	case BOBOL_LOD_PROGRESS_DISPLAY_TERMINAL_ERROR:
	    return QColor(255, 90, 80);
	case BOBOL_LOD_PROGRESS_DISPLAY_BACKGROUND:
	    return QColor(112, 235, 135);
	default:
	    return QColor(255, 190, 72);
    }
}

static inline void
qgcanvas_paint_lod_progress_overlay(QgCanvasState &s, QWidget *w,
	QPainter &painter)
{
    s.lod_progress_overlay_dirty = false;
    const QgLodProgressOverlayState &overlay = s.lod_progress_overlay;
    if (!w || !overlay.visible || !qgcanvas_native_lod_progress_selected(s))
	return;

    static constexpr int margin = 12;
    static constexpr int horizontalPadding = 10;
    static constexpr int verticalPadding = 7;
    static constexpr int textGap = 3;
    static constexpr int progressHeight = 4;
    static constexpr int cornerRadius = 5;
    static constexpr int preferredCardWidth = 320;
    QFont titleFont = w->font();
    titleFont.setBold(true);
    QFont detailFont = w->font();
    detailFont.setPointSizeF(std::max(7.0,
	detailFont.pointSizeF() > 0.0 ? detailFont.pointSizeF() - 1.0 : 8.0));
    const QFontMetrics titleMetrics(titleFont);
    const QFontMetrics detailMetrics(detailFont);
    const int maximumWidth = std::max(1, w->width() - 2 * margin);
    const int boxWidth = std::min(maximumWidth,
	preferredCardWidth);
    const int boxHeight = 2 * verticalPadding + titleMetrics.height() +
	textGap + detailMetrics.height() + textGap + progressHeight;
    const QRectF box(margin, margin, boxWidth, boxHeight);
    const QColor accent = qgcanvas_lod_overlay_color(overlay);

    painter.save();
    painter.resetTransform();
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(18, 22, 28, 220));
    painter.drawRoundedRect(box, cornerRadius, cornerRadius);

    const int availableTextWidth = std::max(1,
	boxWidth - 2 * horizontalPadding);
    painter.setFont(titleFont);
    painter.setPen(QColor(245, 247, 250));
    painter.drawText(margin + horizontalPadding,
	margin + verticalPadding + titleMetrics.ascent(),
	titleMetrics.elidedText(overlay.title, Qt::ElideRight,
	    availableTextWidth));
    int progressTop = margin + verticalPadding + titleMetrics.height();
    painter.setFont(detailFont);
    painter.setPen(QColor(205, 211, 219));
    const int detailTop = progressTop + textGap;
    painter.drawText(margin + horizontalPadding,
	detailTop + detailMetrics.ascent(),
	detailMetrics.elidedText(overlay.detail, Qt::ElideRight,
	    availableTextWidth));
    progressTop = detailTop + detailMetrics.height();
    progressTop += textGap;
    const QRectF track(margin + horizontalPadding, progressTop,
	availableTextWidth, progressHeight);
    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(72, 78, 86));
    painter.drawRoundedRect(track, progressHeight / 2.0,
	progressHeight / 2.0);
    QRectF fill = track;
    if (overlay.determinate) {
	fill.setWidth(track.width() *
	    static_cast<qreal>(overlay.percent) / 100.0);
    } else {
	static constexpr qreal segmentFraction = 0.25;
	static constexpr unsigned int halfSweep = 10;
	const unsigned int step = overlay.animationStep <= halfSweep ?
	    overlay.animationStep : 2 * halfSweep - overlay.animationStep;
	const qreal offset = static_cast<qreal>(step) /
	    static_cast<qreal>(halfSweep) *
	    track.width() * (1.0 - segmentFraction);
	fill.setLeft(track.left() + offset);
	fill.setWidth(track.width() * segmentFraction);
    }
    if (fill.width() > 0.0) {
	painter.setBrush(accent);
	painter.drawRoundedRect(fill, progressHeight / 2.0,
	    progressHeight / 2.0);
    }
    painter.restore();
}

/**
 * Synchronize LoD progress on state transitions and at a bounded cadence
 * while work is active.  A retained owner publishes scene records; a native
 * Qt owner only updates the lightweight widget overlay.  This helper is
 * deliberately usable both immediately before a presentation and after its
 * timing feedback.  The post-render call must permit state transitions only:
 * a periodic retained HUD mutation there would request its own next frame.
 */
static inline bool
qgcanvas_sync_obol_lod_progress(QgCanvasState &s, bool allowPeriodic)
{
    if (!s.v || !s.obol)
	return false;

    const auto now = std::chrono::steady_clock::now();
    const bool lod_pending =
	s.obol->hasProgressiveWorkPending() ||
	s.obol->isLodInteractionActive();
    BObolLodConvergenceStatus lod_status;
    s.obol->getLodConvergenceStatus(lod_status);
    const BObolLodProgressDisplayStatus display =
	lod_status.progressDisplayStatus();
    const QgLodProgressOverlayState overlay =
	qgcanvas_stabilize_lod_progress_overlay(
	    s.lod_progress_presentation, lod_status,
	    qgcanvas_lod_progress_overlay_state(lod_status));
    if (overlay != s.lod_progress_overlay) {
	s.lod_progress_overlay = overlay;
	s.lod_progress_overlay_dirty = true;
    }
    const bool lod_first =
	s.lod_progress_last_publish.time_since_epoch().count() == 0;
    /* The host-work latch may clear one coordinator transition before the
     * convergence state machine publishes IDLE.  The controller-owned display
     * classification makes that final HUD removal a visible state change even
     * when no geometry work remains. */
    const bool lod_state_changed = lod_first ||
	display != s.lod_progress_last_state;
    const bool lod_publish = lod_state_changed ||
	(allowPeriodic && lod_pending && (lod_first ||
	    std::chrono::duration_cast<std::chrono::milliseconds>(
		now - s.lod_progress_last_publish).count() >= 100));
    if (!lod_publish)
	return false;

    s.lod_progress_last_state = display;
    s.lod_progress_last_publish = now;
    struct ged_view_context *view_ctx = ged_view_context_from_bv(s.v);
    struct ged *gedp = ged_view_context_owner(view_ctx);
    enum ged_view_lod_progress_presentation_mode presentation_mode =
	GED_VIEW_LOD_PROGRESS_PRESENTATION_NATIVE_HOST;
    const bool retained = gedp &&
	(!ged_view_lod_progress_presentation_mode_get(&presentation_mode,
	    view_ctx) || presentation_mode ==
	    GED_VIEW_LOD_PROGRESS_PRESENTATION_RETAINED);
    if (retained)
	(void)ged_view_lod_progress_sync(gedp, view_ctx);
    return retained;
}

/* Publish the interaction transition synchronously with the input event.
 * Waiting for the next paint leaves the retained terminal label visible
 * beside a controller which has already re-entered interactive/refining
 * state.  The paint itself remains queued and coalescible; only the cheap HUD
 * record mutation is immediate. */
static inline void
qgcanvas_set_obol_pointer_interaction(
    QgCanvasState &s, QWidget *w, bool active)
{
    if (!s.obol || !w || s.lod_pointer_interaction_active == active)
	return;

    if (active)
	s.obol->beginLodInteraction();
    else
	s.obol->endLodInteraction();
    s.lod_pointer_interaction_active = active;
    if (qgcanvas_sync_obol_lod_progress(s, false)) {
	qgcanvas_request_obol_render_if_idle(s, "lod-interaction-hud");
	w->update();
    } else
	qgcanvas_request_lod_overlay_repaint(s, w, true);
    qgcanvas_queue_obol_progressive_update(s, w);
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

    /* completeRenderTiming() may change settling -> ready, which needs one
     * final HUD frame.  It must not start the next periodic progress frame;
     * periodic samples are folded into independently requested renders by the
     * pre-render call above. */
    const bool lod_publish = qgcanvas_sync_obol_lod_progress(s, false);
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
    else
	qgcanvas_request_lod_overlay_repaint(s, w, true);

    /* Presenting a frame is itself a control transition.  The completed
     * timing sample may open capacity-search, handoff, or demand-rescan work
     * after the previous progressive timer has gone idle.  Always re-evaluate
     * the level-triggered host-work predicate here so that work cannot depend
     * on an unrelated mouse event or paint request to make progress. */
    qgcanvas_queue_obol_progressive_update(s, w);

    /* Never make the FPS label a periodic scene-render timer.  The pre-render
	 * faceplate synchronization folds its newest value into every independently
	 * requested geometry/HUD frame, and the state-transition publication above
	 * supplies the final settled label.  Republishing from frame completion
	 * requested another whole-scene traversal; on OSMesa that feedback loop
	 * consumed the owner thread while immutable mesh results waited to be
	 * admitted. */

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
    if (!s.presentation_fbo || !s.presentation_staging_fbo ||
	s.presentation_fbo->size() != renderSize ||
	s.presentation_staging_fbo->size() != renderSize) {
	delete s.presentation_fbo;
	delete s.presentation_staging_fbo;
	s.presentation_fbo_has_completed_frame = false;
	QOpenGLFramebufferObjectFormat format;
	format.setAttachment(QOpenGLFramebufferObject::CombinedDepthStencil);
	format.setSamples(0);
	s.presentation_fbo = new QOpenGLFramebufferObject(renderSize, format);
	s.presentation_staging_fbo =
	    new QOpenGLFramebufferObject(renderSize, format);
    }
    if (!s.presentation_fbo || !s.presentation_fbo->isValid() ||
	!s.presentation_staging_fbo ||
	!s.presentation_staging_fbo->isValid()) {
	delete s.presentation_fbo;
	s.presentation_fbo = nullptr;
	delete s.presentation_staging_fbo;
	s.presentation_staging_fbo = nullptr;
	s.presentation_fbo_has_completed_frame = false;
	return s.obol->renderPending(clearWindow, clearZBuffer, NULL);
    }

    /* Coin is allowed to mutate only staging.  renderPending() returns false
     * for a deadline abort or an incomplete resumable CAD traversal; in that
     * case retain and blit the last completed framebuffer instead. */
    QOpenGLFramebufferObject *renderingFbo =
	s.presentation_staging_fbo;
    renderingFbo->bind();
    uint64_t featureRevision = 0;
    const SbBool rendered =
	s.obol->renderPending(clearWindow, clearZBuffer, NULL, &featureRevision);
    if (rendered) {
	std::swap(s.presentation_fbo, s.presentation_staging_fbo);
	s.presentation_fbo_has_completed_frame = true;
	s.completed_feature_revision = featureRevision;
    }
    QOpenGLContext *context = QOpenGLContext::currentContext();
    QOpenGLExtraFunctions *gl = context ? context->extraFunctions() : nullptr;
    if (gl) {
	/* On the first interrupted draw there is no prior frame.  Preserve the
	 * useful prefix Coin completed before the deadline, as the former single
	 * FBO path did.  Once an exact frame exists, staging is never observable. */
	QOpenGLFramebufferObject *source =
	    s.presentation_fbo_has_completed_frame ? s.presentation_fbo :
	    s.presentation_staging_fbo;
	gl->glBindFramebuffer(GL_READ_FRAMEBUFFER,
	    source->handle());
	gl->glBindFramebuffer(GL_DRAW_FRAMEBUFFER,
	    widget->defaultFramebufferObject());
	gl->glBlitFramebuffer(0, 0, renderSize.width(), renderSize.height(),
	    0, 0, renderSize.width(), renderSize.height(),
	    GL_COLOR_BUFFER_BIT, GL_NEAREST);
	gl->glBindFramebuffer(GL_FRAMEBUFFER,
	    widget->defaultFramebufferObject());
	if (s.presentation_fbo_has_completed_frame)
	    s.presented_feature_revision = s.completed_feature_revision;
	else
	    s.presented_feature_revision.reset();
    } else {
	renderingFbo->release();
    }
    return rendered;
}

/** Store current view hashes in @p s. */
static inline void
qgcanvas_stash_hashes(QgCanvasState &s)
{
    s.prev_dhash = 0;
    const struct bv *view = bv_context_view_const(s.v);
    s.prev_vhash = bv_hash(view);
    s.prev_frame_revision = bv_frame_revision_get(view);
}

/** Request a semantic view refresh and wake the canvas backend. */
static inline void
qgcanvas_request_update(QgCanvasState &s, uint32_t flags)
{
    uint32_t requested = flags ? flags : BV_REFRESH_ALL;
    if (requested & BV_REFRESH_VIEW)
	qgcanvas_sync_obol_camera(s);
    qgcanvas_sync_obol_faceplate(s);

    /* Qt update() only asks the widget to paint.  When a completed Obol
     * framebuffer is retained, painting without a controller request is an
     * intentional zero-copy replay of that frame.  Semantic scene, overlay,
     * and edit changes therefore must latch one presentation render here.
     * Camera-only refreshes remain excluded: sync_obol_camera requests a
     * capacity-relevant frame itself when (and only when) the camera changed.
     * This preserves passive retained-frame repaint without allowing GED
     * erase/redraw or selection changes to leave obsolete pixels onscreen. */
    if (s.obol && (requested &
	    (BV_REFRESH_DRAW | BV_REFRESH_OVERLAY | BV_REFRESH_EDIT)))
	s.obol->requestPresentationRender("qt-semantic-refresh");

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
    const uint64_t frameRevision = bv_frame_revision_get(view);

    /* Commands may change retained presentation state without changing the
     * camera hash (lighting, faceplate overlays, and similar policy).  The
     * frame revision is the durable record of that semantic change.  The
     * dirty latch remains useful for explicit redraws with unchanged state,
     * but a fast renderer is allowed to consume it before this comparison. */
    if (s.prev_vhash != c_vhash ||
	s.prev_frame_revision != frameRevision ||
	bv_refresh_dirty_get(view)) {
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
