/*                      Q G G L . C P P
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
/** @file QgGL.cpp
 *
 * Use a QOpenGLWidget to display Obol/Coin drawing content.
 *
 */

#include "common.h"

#include "bv.h"

#include <QImage>
#include <QKeyEvent>
#include <QMetaMethod>
#include <QMouseEvent>
#include <QOpenGLWidget>
#include <QResizeEvent>
#include <QWheelEvent>
#include <QtGlobal>

#include "QgCanvasState.h"   /* pimpl definition + shared helpers */
#include "qtcad/QgGL.h"

QgGL::QgGL(QWidget *parent, BObolViewController *controller,
    bool create_controller)
    : QOpenGLWidget(parent)
{
    d = new QgCanvasState();
    qgcanvas_init_obol(*d, this, false, controller, create_controller);
    d->lmouse_mode = BV_ADJUST_SCALE;

    // Provide a view specific to this widget - set gedp->ged_gvp to v
    // if this is the current view
    d->v = qgcanvas_view_context_create("qtgl");
    qgcanvas_sync_obol_camera(*d);
    qgcanvas_initialize_obol_background(*d);

    /* A Qt update may reach paintGL before asynchronous Obol work has
     * published a replacement frame.  Preserve the last complete image in
     * that case.  NoPartialUpdate invalidates the QOpenGLWidget FBO before
     * every paint; a renderPending() no-op would then hand Qt undefined color
     * storage, observed on NVIDIA as black/RGB-dotted zoom frames. */
    setUpdateBehavior(QOpenGLWidget::PartialUpdate);

    // This is an important Qt setting for interactivity - it allowing key
    // bindings to propagate to this widget and trigger actions such as
    // resolution scaling, rotation, etc.
    setFocusPolicy(Qt::WheelFocus);
    setMouseTracking(true);
}

QgGL::~QgGL()
{
	d->input.setEndpoint(NULL);
    qgcanvas_destroy_obol(*d, this);
    qgcanvas_view_context_destroy(d->v);
    d->v = nullptr;
    delete d;
    d = nullptr;
}

struct bv_context *
QgGL::viewContext() const
{
    return d ? d->v : nullptr;
}

BObolViewController *
QgGL::obolViewController() const
{
    return d->obol;
}

void
QgGL::setObolViewController(BObolViewController *controller)
{
    qgcanvas_bind_obol_controller(*d, this, controller);
}

void
QgGL::setObolInputEndpoint(struct bobol_display_endpoint *endpoint)
{
	d->input.setEndpoint(endpoint);
}

int
QgGL::currentView() const
{
    return d->current;
}

void
QgGL::set_current(int active)
{
    d->current = active;
}

void QgGL::paintGL()
{
    int w = width();
    int h = height();
    // Zero size == nothing to do
    if (!w || !h)
return;

    qgcanvas_sync_obol_viewport(*d, this);
    qgcanvas_sync_obol_camera(*d);
    initializeOpenGLFunctions();
    /* An expose/compositor paint is not a semantic scene update.  Retain and
     * blit the completed presentation FBO unless the controller has an
     * explicit render request; manufacture a request only when this viewport
     * has no usable completed frame yet. */
    if (!qgcanvas_has_completed_gl_frame(*d, this))
	qgcanvas_request_obol_render_if_idle(*d, "qtgl-initial-paint");
    /* A provider/idle transition may have occurred since the preceding
     * frame.  Apply its retained HUD delta before rendering this frame so a
     * completed view never presents one final stale progress indicator. */
    (void)qgcanvas_sync_obol_lod_progress(*d);

    const SbBool rendered =
	qgcanvas_render_obol_pending(*d, this, TRUE, TRUE);
    if (rendered && d->v) {
	(void)bv_refresh_consume(bv_context_view(d->v));
	bv_refresh_complete(bv_context_view(d->v));
    }
    if (rendered)
	qgcanvas_frame_complete(*d, this);
    qgcanvas_queue_obol_progressive_update(*d, this);
    /*
     * A test recorder must observe the framebuffer already produced by this
     * paint.  QOpenGLWidget::grabFramebuffer() from frameSwapped can request
     * another paint and thereby change the behavior being measured.  Do the
     * readback here, while Qt's FBO is current, and pay nothing unless a
     * diagnostic observer is connected.
     */
    if (isSignalConnected(QMetaMethod::fromSignal(
	    &QgGL::frame_presented))) {
	const QSize renderSize = qgcanvas_render_size(this);
	QImage image(renderSize, QImage::Format_RGBA8888);
	if (!image.isNull()) {
	    GLint oldPackAlignment = 4;
	    glGetIntegerv(GL_PACK_ALIGNMENT, &oldPackAlignment);
	    glPixelStorei(GL_PACK_ALIGNMENT, 4);
	    glReadPixels(0, 0, renderSize.width(), renderSize.height(),
		GL_RGBA, GL_UNSIGNED_BYTE, image.bits());
	    glPixelStorei(GL_PACK_ALIGNMENT, oldPackAlignment);
	    image = image.flipped(Qt::Vertical);
	    image.setDevicePixelRatio(devicePixelRatioF());
	    emit frame_presented(image);
	}
    }
    if (!d->obol_paint_initialized) {
	d->obol_paint_initialized = true;
	emit init_done();
    }
}

void QgGL::resizeGL(int, int)
{
    qgcanvas_sync_obol_viewport(*d, this);
    if (d->v)
	qgcanvas_request_update(*d, BV_REFRESH_VIEW);
    emit changed();
}

void QgGL::resizeEvent(QResizeEvent *e)
{
    QOpenGLWidget::resizeEvent(e);
    qgcanvas_sync_obol_viewport(*d, this);
    if (!d->v)
return;
    QSize rsize = qgcanvas_render_size(this);
    bv_context_dimensions_set(d->v, rsize.width(), rsize.height());
    qgcanvas_request_update(*d, BV_REFRESH_VIEW);
    emit changed();
}

void QgGL::request_update(uint32_t refresh_flags)
{
    uint32_t requested = refresh_flags ? refresh_flags : BV_REFRESH_ALL;
    qgcanvas_request_update(*d, requested);
    if (d->fb_update_queued)
return;
    d->fb_update_queued = true;
    QMetaObject::invokeMethod(this, "queued_update", Qt::QueuedConnection);
}

void QgGL::need_update()
{
    QTCAD_SLOT("QgGL::need_update", 1);
    request_update(BV_REFRESH_VIEW);
}

void QgGL::present_frame()
{
    QTCAD_SLOT("QgGL::present_frame", 1);
    request_update(BV_REFRESH_FRAMEBUFFER | BV_REFRESH_FORCE);
}

void QgGL::queued_update()
{
    d->fb_update_queued = false;
    update();
}

void QgGL::keyPressEvent(QKeyEvent *k)
{

    if (!d->v || !d->current || !d->use_default_keybindings) {
QOpenGLWidget::keyPressEvent(k);
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

    QOpenGLWidget::keyPressEvent(k);
}

void QgGL::mousePressEvent(QMouseEvent *e)
{

    if (!d->v || !d->current || !d->use_default_mousebindings) {
QOpenGLWidget::mousePressEvent(e);
return;
    }

    // Let bv know what the current view width and height are, in
    // case the dx/dy mouse translations need that information
    QSize rsize = qgcanvas_render_size(this);
    bv_context_dimensions_set(d->v, rsize.width(), rsize.height());

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

    QOpenGLWidget::mousePressEvent(e);
}

void QgGL::mouseReleaseEvent(QMouseEvent *e)
{
    if (d->obol && e->button() == Qt::LeftButton)
	d->obol->endLodInteraction();

    if (!d->v) {
QOpenGLWidget::mouseReleaseEvent(e);
return;
    }

    // To avoid an abrupt jump in scene motion the next time movement is
    // started with the mouse, after we release we return to the default state.
    d->x_prev = -INT_MAX;
    d->y_prev = -INT_MAX;

    if (d->input.mouseReleaseEvent(d->v,
	    d->x_press_pos, d->y_press_pos, d->x_prev, d->y_prev, e,
	    d->lmouse_mode)) {
qgcanvas_request_update(*d, BV_REFRESH_VIEW);
update();
emit changed();
    }

    QOpenGLWidget::mouseReleaseEvent(e);
}

void QgGL::mouseMoveEvent(QMouseEvent *e)
{
    if (!d->v || !d->current || !d->use_default_mousebindings) {
QOpenGLWidget::mouseMoveEvent(e);
return;
    }

    // Let bv know what the current view width and height are, in
    // case the dx/dy mouse translations need that information
    QSize rsize = qgcanvas_render_size(this);
    bv_context_dimensions_set(d->v, rsize.width(), rsize.height());

    if (d->obol && e->buttons().testFlag(Qt::LeftButton))
	d->obol->beginLodInteraction();
    int mret = d->input.mouseMoveEvent(d->v,
	    d->x_prev, d->y_prev, e, d->lmouse_mode);
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

    QOpenGLWidget::mouseMoveEvent(e);
}

void QgGL::wheelEvent(QWheelEvent *e)
{

    if (!d->v || !d->current || !d->use_default_mousebindings) {
QOpenGLWidget::wheelEvent(e);
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

    QOpenGLWidget::wheelEvent(e);
}

void QgGL::stash_hashes()
{
    qgcanvas_stash_hashes(*d);
}

bool QgGL::diff_hashes()
{
    bool ret = qgcanvas_diff_hashes_check(*d);
    if (ret) {
need_update();
emit changed();
    }
    return ret;
}

void QgGL::save_image()
{
    QImage image;
    get_viewport_image(image);
    if (!image.isNull())
	image.save("file.png");
}

void QgGL::render_to_file(const QString &filename)
{
    QImage img;
    get_viewport_image(img);
    if (!img.isNull())
img.convertToFormat(QImage::Format_RGB32).save(filename);
}

void QgGL::get_viewport_image(QImage &img)
{
    img = QImage();
    if (!d->v)
	return;

    qgcanvas_get_obol_viewport_image(*d, this, img);
}

void QgGL::get_obol_viewport_image(QImage &img)
{
    qgcanvas_get_obol_viewport_image(*d, this, img);
}

void QgGL::aet(double a, double e, double t)
{
    qgcanvas_aet(*d, a, e, t);
}

void
QgGL::enableDefaultKeyBindings()
{
    d->use_default_keybindings = true;
}

void
QgGL::disableDefaultKeyBindings()
{
    d->use_default_keybindings = false;
}

void
QgGL::enableDefaultMouseBindings()
{
    d->use_default_mousebindings = true;
}

void
QgGL::disableDefaultMouseBindings()
{
    d->use_default_mousebindings = false;
}

int
QgGL::lmouseMoveDefault() const
{
    return d->lmouse_mode;
}

void
QgGL::set_lmouse_move_default(int mm)
{
    QTCAD_SLOT("QgGL::set_lmouse_move_default", 1);
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
