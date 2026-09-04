/*                          Q G S W . C P P
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
/** @file QgSW.cpp
 *
 * Qt widget for visualizing Obol/Coin output using offscreen readback.
 */

#define USE_MGL_NAMESPACE 1

#include "common.h"

#include "bv.h"

#include <QImage>
#include <QMetaMethod>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QResizeEvent>
#include <QtGlobal>
#include <QWheelEvent>

#include "QgCanvasState.h"   /* pimpl definition + shared helpers */
#include "qtcad/QgSW.h"

QgSW::QgSW(QWidget *parent, BObolViewController *controller,
    bool create_controller)
    : QWidget(parent)
{
    d = new QgCanvasState();
    qgcanvas_init_obol(*d, this, true, controller, create_controller);
    d->lmouse_mode = BV_ADJUST_SCALE;

    // Provide a view specific to this widget - set gedp->ged_gvp to v
    // if this is the current view
    d->v = qgcanvas_view_context_create("swrast");
    qgcanvas_sync_obol_camera(*d);
    qgcanvas_initialize_obol_background(*d);

    /* Obol always supplies a complete RGB frame.  Mark the canvas opaque so
     * Qt does not repaint the parent behind every software readback blit. */
    setAttribute(Qt::WA_OpaquePaintEvent, true);
    setAutoFillBackground(false);

    // This is an important Qt setting for interactivity - it allowing key
    // bindings to propagate to this widget and trigger actions such as
    // resolution scaling, rotation, etc.
    setFocusPolicy(Qt::WheelFocus);
    setMouseTracking(true);
}

QgSW::~QgSW()
{
	d->input.setEndpoint(NULL);
    qgcanvas_destroy_obol(*d, this);
    qgcanvas_view_context_destroy(d->v);
    d->v = nullptr;
    delete d;
    d = nullptr;
}

struct bv_context *
QgSW::viewContext() const
{
    return d ? d->v : nullptr;
}

BObolViewController *
QgSW::obolViewController() const
{
    return d->obol;
}

void
QgSW::setObolViewController(BObolViewController *controller)
{
    qgcanvas_bind_obol_controller(*d, this, controller);
}

void
QgSW::setObolInputEndpoint(struct bobol_display_endpoint *endpoint)
{
	d->input.setEndpoint(endpoint);
}

int
QgSW::currentView() const
{
    return d->current;
}

void
QgSW::set_current(int active)
{
    d->current = active;
}

bool
QgSW::isPresentationInitialized() const
{
    return d->obol_paint_initialized;
}

void QgSW::request_update(uint32_t refresh_flags)
{
    uint32_t requested = refresh_flags ? refresh_flags : BV_REFRESH_ALL;
    qgcanvas_request_update(*d, requested | BV_REFRESH_FRAMEBUFFER | BV_REFRESH_FORCE);
    if (d->fb_update_queued)
return;
    d->fb_update_queued = true;
    QMetaObject::invokeMethod(this, "queued_update", Qt::QueuedConnection);
}

void QgSW::need_update()
{
    QTCAD_SLOT("QgSW::need_update", 1);
    /* diff_hashes() reaches this path after a GED command.  A framebuffer
     * repaint alone preserves the previous Obol camera, so an A/E command can
     * appear to take effect only when later progressive work happens to sync
     * it.  Match the GL canvas contract: synchronize semantic view state
     * before presenting the retained framebuffer. */
    request_update(BV_REFRESH_VIEW | BV_REFRESH_FRAMEBUFFER |
	BV_REFRESH_FORCE);
}

void QgSW::present_frame()
{
    QTCAD_SLOT("QgSW::present_frame", 1);
    /* A software frame can be requested solely by a background realization;
     * unlike a direct GL expose it may be the first host callback following a
     * GED command.  Reconcile the authoritative view at that presentation
     * boundary so a warm realization cannot show the controller's startup
     * camera until a later mouse event. */
    request_update(BV_REFRESH_VIEW | BV_REFRESH_FRAMEBUFFER |
	BV_REFRESH_FORCE);
}

void QgSW::queued_update()
{
    QTCAD_SLOT("QgSW::queued_update", 1);
    d->fb_update_queued = false;
    update();
}

void QgSW::paintEvent(QPaintEvent *e)
{
    QTCAD_SLOT("QgSW::paintEvent", 1);
    if (!d->v)
return;

    QSize rsize = qgcanvas_render_size(this);
    bv_context_dimensions_set(d->v, rsize.width(), rsize.height());
    qgcanvas_sync_obol_viewport(*d, this);
    qgcanvas_sync_obol_camera(*d);

    /* QWidget may repaint because it was exposed, uncovered, or composited;
     * none of those events changed the CAD presentation.  Re-entering Coin
     * and OSMesa for every such paint used to manufacture a continuous
     * render/calibration loop and made a pose-only LoD restore look like a
     * coarse-to-fine replay.  The controller's render latch is the authority
     * for semantic changes.  Only cold start (or a size mismatch) needs an
     * otherwise-idle presentation request. */
    const bool retainedFrameMatches =
	qgcanvas_has_completed_software_frame(*d, this);
    if (!retainedFrameMatches)
	qgcanvas_request_obol_render_if_idle(*d, "qtsw-initial-paint");

    QImage image;
    bool completedPresentation = false;
    if (retainedFrameMatches && d->obol &&
	!d->obol->isRenderRequested()) {
	/* Renderer-native pixels are bottom-up, matching the inverse QPainter
	 * transform below.  QImage is implicitly shared, so this is a zero-copy
	 * presentation of the immutable completed frame. */
	image = d->last_completed_software_frame;
    } else {
	qgcanvas_get_obol_viewport_image(
	    *d, this, image, true, true, true, &completedPresentation);
    }
    if (image.isNull()) {
	/* Preserve the opaque-widget contract if an offscreen render fails. */
	{
	    QPainter painter(this);
	    painter.fillRect(e->rect(), Qt::black);
	}
	QImage black(qgcanvas_render_size(this), QImage::Format_RGBA8888);
	black.fill(Qt::black);
	black.setDevicePixelRatio(devicePixelRatioF());
	if (isSignalConnected(QMetaMethod::fromSignal(
		&QgSW::frame_presented)))
	    emit frame_presented(black);
	/*
	 * A cold progressive source may not have published even its scope
	 * bounds by the first traversal.  That is a valid transient blank
	 * presentation, not an idle pipeline.  Keep the provider pump armed or
	 * this early return loses the only wakeup that can replace the blank
	 * frame with the newly published boxes and meshes.
	 */
	qgcanvas_queue_obol_progressive_update(*d, this);
	if (!d->obol_paint_initialized) {
	    d->obol_paint_initialized = true;
	    emit init_done();
	}
	return;
    }
    {
	QPainter painter(this);
	painter.translate(0, height());
	painter.scale(1, -1);
	painter.drawImage(QPoint(0, 0), image);
    }
    if (completedPresentation) {
	(void)bv_refresh_consume(bv_context_view(d->v));
	bv_refresh_complete(bv_context_view(d->v));
	qgcanvas_frame_complete(*d, this);
    }
    qgcanvas_queue_obol_progressive_update(*d, this);
    if (!d->obol_paint_initialized) {
	d->obol_paint_initialized = true;
	emit init_done();
    }
    if (isSignalConnected(QMetaMethod::fromSignal(
	    &QgSW::frame_presented))) {
	/* paintEvent presents the renderer's bottom-up buffer through the
	 * inverse QPainter transform above.  Publish the exact displayed
	 * orientation to diagnostic observers. */
	QImage presented = image.flipped(Qt::Vertical);
	presented.setDevicePixelRatio(devicePixelRatioF());
	emit frame_presented(presented);
    }
}

void QgSW::resizeEvent(QResizeEvent *e)
{
    QWidget::resizeEvent(e);
    qgcanvas_sync_obol_viewport(*d, this);
    if (d->v) {
	QSize rsize = qgcanvas_render_size(this);
	bv_context_dimensions_set(d->v, rsize.width(), rsize.height());
	qgcanvas_request_update(*d, BV_REFRESH_VIEW | BV_REFRESH_FRAMEBUFFER);
	emit changed();
    }
}

void QgSW::keyPressEvent(QKeyEvent *k)
{

    if (!d->v || !d->current || !d->use_default_keybindings) {
QWidget::keyPressEvent(k);
return;
    }

    // Let bv know what the current view width and height are, in
    // case the dx/dy mouse translations need that information
    QSize rsize = qgcanvas_render_size(this);
    bv_context_dimensions_set(d->v, rsize.width(), rsize.height());

    if (d->input.keyPressEvent(d->v, d->x_prev,
	    d->y_prev, k)) {
qgcanvas_request_update(*d, BV_REFRESH_VIEW);
update();
emit changed();
    }

    QWidget::keyPressEvent(k);
}

void QgSW::mousePressEvent(QMouseEvent *e)
{

    if (!d->v || !d->current || !d->use_default_mousebindings) {
QWidget::mousePressEvent(e);
return;
    }

    // Let bv know what the current view width and height are, in
    // case the dx/dy mouse translations need that information
    QSize rsize = qgcanvas_render_size(this);
    bv_context_dimensions_set(d->v, rsize.width(), rsize.height());
    d->input.setDevicePixelRatio(devicePixelRatioF());

    if (d->input.mousePressEvent(d->v, d->x_prev,
	    d->y_prev, e)) {
qgcanvas_request_update(*d, BV_REFRESH_VIEW);
update();
emit changed();
    }

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    d->x_press_pos = (double)e->x();
    d->y_press_pos = (double)e->y();
#else
    d->x_press_pos = e->position().x();
    d->y_press_pos = e->position().y();
#endif
    //bu_log("X,Y: %g, %g\n", d->x_press_pos, d->y_press_pos);

    QWidget::mousePressEvent(e);
}

void QgSW::mouseReleaseEvent(QMouseEvent *e)
{

    if (d->obol && d->lod_pointer_interaction_active &&
	e->button() == Qt::LeftButton) {
	qgcanvas_set_obol_pointer_interaction(*d, this, false);
    }
    if (!d->v) {
QWidget::mouseReleaseEvent(e);
return;
    }

    // To avoid an abrupt jump in scene motion the next time movement is
    // started with the mouse, after we release we return to the default state.
    d->x_prev = -INT_MAX;
    d->y_prev = -INT_MAX;
    d->input.setDevicePixelRatio(devicePixelRatioF());

    if (d->input.mouseReleaseEvent(d->v,
	    d->x_press_pos, d->y_press_pos, d->x_prev, d->y_prev, e,
	    d->lmouse_mode)) {
qgcanvas_request_update(*d, BV_REFRESH_VIEW);
update();
emit changed();
    }

    QWidget::mouseReleaseEvent(e);
}


void QgSW::mouseMoveEvent(QMouseEvent *e)
{
    if (!d->v || !d->current || !d->use_default_mousebindings) {
QWidget::mouseMoveEvent(e);
return;
    }

    // Let bv know what the current view width and height are, in
    // case the dx/dy mouse translations need that information
    QSize rsize = qgcanvas_render_size(this);
    bv_context_dimensions_set(d->v, rsize.width(), rsize.height());
    d->input.setDevicePixelRatio(devicePixelRatioF());

    int mret = d->input.mouseMoveEvent(d->v,
	    d->x_prev, d->y_prev, e, d->lmouse_mode);
    if (d->obol && e->buttons().testFlag(Qt::LeftButton) &&
	d->input.lastDispatchWasViewMotion()) {
	qgcanvas_set_obol_pointer_interaction(*d, this, true);
    }
    if (mret > 0) {
qgcanvas_request_update(*d, BV_REFRESH_VIEW);
update();
emit changed();
    }

    // Current positions are the new previous positions
    if (mret != -1) {
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
d->x_prev = e->x();
d->y_prev = e->y();
#else
d->x_prev = e->position().x();
d->y_prev = e->position().y();
#endif
    }

    QWidget::mouseMoveEvent(e);
}

void QgSW::wheelEvent(QWheelEvent *e)
{

    if (!d->v || !d->current || !d->use_default_mousebindings) {
QWidget::wheelEvent(e);
return;
    }

    // Let bv know what the current view width and height are, in
    // case the dx/dy mouse translations need that information
    QSize rsize = qgcanvas_render_size(this);
    bv_context_dimensions_set(d->v, rsize.width(), rsize.height());

    if (d->input.wheelEvent(d->v, e)) {
qgcanvas_request_update(*d, BV_REFRESH_VIEW);
update();
emit changed();
    }

    QWidget::wheelEvent(e);
}

void QgSW::stash_hashes()
{
    qgcanvas_stash_hashes(*d);
}

bool QgSW::diff_hashes()
{
    bool ret = qgcanvas_diff_hashes_check(*d);
    if (ret) {
need_update();
emit changed();
    }
    return ret;
}

void QgSW::save_image()
{
    QImage image;
    get_viewport_image(image);
    if (!image.isNull())
	image.save("file.png");
}

/* Render the current view to a file without relying on Qt paint events.
 * This is safe to call in headless / off-screen mode. */
void QgSW::render_to_file(const QString &filename)
{
    QImage img;
    get_viewport_image(img);
    if (!img.isNull())
img.convertToFormat(QImage::Format_RGB32).save(filename);
}

/* Render the current Obol view and return the pixel data. */
void QgSW::get_viewport_image(QImage &img)
{
    img = QImage();  /* null sentinel */
    if (!d->v) return;

    qgcanvas_get_obol_viewport_image(*d, this, img);
}

void QgSW::get_obol_viewport_image(QImage &img)
{
    qgcanvas_get_obol_viewport_image(*d, this, img);
}

void QgSW::aet(double a, double e, double t)
{
    qgcanvas_aet(*d, a, e, t);
}

void
QgSW::enableDefaultKeyBindings()
{
    d->use_default_keybindings = true;
}

void
QgSW::disableDefaultKeyBindings()
{
    d->use_default_keybindings = false;
}

void
QgSW::enableDefaultMouseBindings()
{
    d->use_default_mousebindings = true;
}

void
QgSW::disableDefaultMouseBindings()
{
    d->use_default_mousebindings = false;
}

int
QgSW::lmouseMoveDefault() const
{
    return d->lmouse_mode;
}

void
QgSW::set_lmouse_move_default(int mm)
{
    QTCAD_SLOT("QgSW::set_lmouse_move_default", 1);
    d->lmouse_mode = mm;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
