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

#include <QImage>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QOpenGLWidget>
#include <QResizeEvent>
#include <QWheelEvent>
#include <QtGlobal>

#include "QgCanvasState.h"   /* pimpl definition + shared helpers */
#include "QgLegacyViewContext.h"
#include "qtcad/QgGL.h"

QgGL::QgGL(QWidget *parent)
    : QOpenGLWidget(parent)
{
    d = new QgCanvasState();
    qgcanvas_init_obol(*d, this);
    d->lmouse_mode = BV_ADJUST_SCALE;

    // Provide a view specific to this widget - set gedp->ged_gvp to v
    // if this is the current view
    d->local_v = qg_legacy_view_local_create("qtgl");
    d->v = d->local_v;
    qgcanvas_sync_obol_camera(*d);

    // This is an important Qt setting for interactivity - it allowing key
    // bindings to propagate to this widget and trigger actions such as
    // resolution scaling, rotation, etc.
    setFocusPolicy(Qt::WheelFocus);
}

QgGL::~QgGL()
{
    qgcanvas_destroy_obol(*d);
    qg_legacy_view_local_free(d->local_v);
    d->local_v = nullptr;
    d->v = nullptr;
    delete d;
    d = nullptr;
}

qg_legacy_view *
QgGL::view() const
{
    return d->v;
}

bool
QgGL::legacyBackendInitialized() const
{
    return false;
}

BRLObolViewController *
QgGL::obolViewController() const
{
    return d->obol;
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

void
QgGL::set_view(qg_legacy_view *nv)
{
    qgcanvas_set_view(*d, nv);
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
    qgcanvas_request_obol_render_if_idle(*d, "qtgl-paint");

    if (qgcanvas_render_obol_pending(*d, TRUE, TRUE) && d->v) {
	(void)bv_refresh_consume(qg_legacy_view_bv(d->v));
	bv_refresh_complete(qg_legacy_view_bv(d->v));
    }
    qgcanvas_frame_complete(*d, this);
    qgcanvas_queue_obol_progressive_update(*d, this);
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
    qg_legacy_view_dimensions_set(d->v, rsize.width(), rsize.height());
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
    qg_legacy_view_dimensions_set(d->v, rsize.width(), rsize.height());

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
    qg_legacy_view_dimensions_set(d->v, rsize.width(), rsize.height());

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
    qg_legacy_view_dimensions_set(d->v, rsize.width(), rsize.height());

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
    qg_legacy_view_dimensions_set(d->v, rsize.width(), rsize.height());

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

    qgcanvas_get_obol_viewport_image(*d, this, img, true);
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
