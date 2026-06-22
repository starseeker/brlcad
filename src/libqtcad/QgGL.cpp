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
 * Use a QOpenGLWidget to display libdm drawing content.
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
#include "qtcad/QgGL.h"

// FROM MGED
#define XMIN            (-2048)
#define XMAX            (2047)
#define YMIN            (-2048)
#define YMAX            (2047)

// from BSG_VIEW_MIN and BSG_VIEW_MAX
#define QTGL_ZMIN -2048
#define QTGL_ZMAX 2047

static thread_local qg_legacy_fb *qggl_bridge_framebuffer = nullptr;

struct QgGLBridgeFramebufferScope {
    explicit QgGLBridgeFramebufferScope(qg_legacy_fb *fbp)
	: previous(qggl_bridge_framebuffer)
    {
	qggl_bridge_framebuffer = fbp;
    }

    ~QgGLBridgeFramebufferScope()
    {
	qggl_bridge_framebuffer = previous;
    }

    qg_legacy_fb *previous = nullptr;
};

QgGL::QgGL(QWidget *parent)
    : QOpenGLWidget(parent)
{
    d = new QgCanvasState();
    qgcanvas_init_obol(*d, this);
    d->ifp = qggl_bridge_framebuffer;
    d->lmouse_mode = QG_LEGACY_VIEW_ADJUST_SCALE;

    // Provide a view specific to this widget - set gedp->ged_gvp to v
    // if this is the current view
    d->local_v = qg_legacy_view_local_create("qtgl");
    d->v = d->local_v;
    qgcanvas_sync_obol_camera(*d);

    // We can't initialize dmp successfully until more of the OpenGL
    // initialization is complete
    d->dmp = nullptr;

    // If we weren't supplied with a framebuffer, allocate one.
    // We don't open it until we have the dmp.
    if (!d->ifp)
d->ifp = qg_legacy_view_framebuffer_raw_create("qtgl");

    // This is an important Qt setting for interactivity - it allowing key
    // bindings to propagate to this widget and trigger actions such as
    // resolution scaling, rotation, etc.
    setFocusPolicy(Qt::WheelFocus);
}

QgGL *
QgCanvasBridgeFactory::create_qtgl(QWidget *parent, qg_legacy_fb *fbp)
{
    QgGLBridgeFramebufferScope scope(fbp);
    return new QgGL(parent);
}

QgGL::~QgGL()
{
    qg_legacy_view_dm_close(d->dmp);
    qg_legacy_view_framebuffer_release(d->ifp, d->m_init);
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
    return d->dmp != nullptr;
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

    if (qgcanvas_obol_scene_has_content(*d)) {
	qgcanvas_sync_obol_viewport(*d, this);
	qgcanvas_sync_obol_camera(*d);
	initializeOpenGLFunctions();
	if (qgcanvas_render_obol_pending(*d, TRUE, TRUE)) {
	    if (d->v) {
		(void)qg_legacy_view_refresh_consume(d->v);
		qg_legacy_view_refresh_complete(d->v);
	    }
	    if (d->dmp)
		qg_legacy_view_dm_native_repaint_pending_set(d->dmp, 0);
	    if (!d->obol_paint_initialized) {
		d->obol_paint_initialized = true;
		emit init_done();
	    }
	    return;
	}
    }

    if (!d->m_init) {

if (!d->dmp) {

    // This is needed so we can work with Qt's OpenGL widget
    // using standard OpenGL functions.
    initializeOpenGLFunctions();

    d->dmp = qg_legacy_view_dm_open_qtgl((void *)this);
    if (!d->dmp)
return;

    // If we have a framebuffer, now we can open it
    (void)qg_legacy_view_dm_framebuffer_setup_existing(d->ifp, d->dmp);
}

// QTGL_ZMIN and QTGL_ZMAX are historical - need better
// documentation on why those specific values are used.
(void)qg_legacy_view_dm_setup_qtgl(d->dmp, QTGL_ZMIN, QTGL_ZMAX);

if (d->v) {
    qg_legacy_view_dm_bind(d->v, d->dmp);
    qg_legacy_view_dm_sync_dimensions(d->v, d->dmp);
}

// Ready to go
d->m_init = true;
emit init_done();
    }

    if (!d->m_init || !d->dmp || !d->v)
return;

    QSize rsize = qgcanvas_render_size(this);
    if (qg_legacy_view_dm_width_get(d->dmp) != rsize.width() ||
	    qg_legacy_view_dm_height_get(d->dmp) != rsize.height()) {
qg_legacy_view_dm_dimensions_set(d->dmp, rsize.width(), rsize.height());
qg_legacy_view_dm_configure_window(d->dmp, 0);
if (d->ifp)
    qg_legacy_view_framebuffer_configure(d->ifp, rsize.width(),
	    rsize.height());
    }
    qg_legacy_view_dm_sync_dimensions(d->v, d->dmp);

    // Re-draw the background to clear any previous drawing
    qg_legacy_view_dm_background_restore(d->dmp);

    // Go ahead and set the flag, but (unlike the rendering thread
    // implementation) we need to do the draw routine every time in paintGL, or
    // we end up with unrendered frames.
    (void)qg_legacy_view_refresh_consume(d->v);
    qg_legacy_view_dm_native_repaint_pending_set(d->dmp, 0);
    qg_legacy_view_dm_draw(d->v);
    qg_legacy_view_dm_draw_end(d->dmp);
    qg_legacy_view_refresh_complete(d->v);
}

void QgGL::resizeGL(int, int)
{
    qgcanvas_sync_obol_viewport(*d, this);
    if (!d->dmp || !d->v)
return;
    qg_legacy_view_dm_configure_window(d->dmp, 0);
    qg_legacy_view_dm_sync_dimensions(d->v, d->dmp);
    if (d->ifp) {
qg_legacy_view_framebuffer_configure(d->ifp,
	qg_legacy_view_width_get(d->v),
	qg_legacy_view_height_get(d->v));
    }
    if (d->dmp)
qgcanvas_request_update(*d, QG_LEGACY_VIEW_REFRESH_VIEW);
    emit changed();
}

void QgGL::resizeEvent(QResizeEvent *e)
{
    QOpenGLWidget::resizeEvent(e);
    qgcanvas_sync_obol_viewport(*d, this);
    if (!d->dmp || !d->v)
return;
    QSize rsize = qgcanvas_render_size(this);
    qg_legacy_view_dm_dimensions_set(d->dmp, rsize.width(), rsize.height());
    qg_legacy_view_dimensions_set(d->v, rsize.width(), rsize.height());
    qg_legacy_view_dm_configure_window(d->dmp, 0);
    if (d->ifp) {
qg_legacy_view_framebuffer_configure(d->ifp,
	qg_legacy_view_width_get(d->v),
	qg_legacy_view_height_get(d->v));
    }
    if (d->dmp)
qgcanvas_request_update(*d, QG_LEGACY_VIEW_REFRESH_VIEW);
    emit changed();
}

void QgGL::request_update(uint32_t refresh_flags)
{
    uint32_t requested = refresh_flags ? refresh_flags : QG_LEGACY_VIEW_REFRESH_ALL;
    qgcanvas_request_update(*d, requested);
    if (d->fb_update_queued)
return;
    d->fb_update_queued = true;
    QMetaObject::invokeMethod(this, "queued_update", Qt::QueuedConnection);
}

void QgGL::need_update()
{
    QTCAD_SLOT("QgGL::need_update", 1);
    request_update(QG_LEGACY_VIEW_REFRESH_VIEW);
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
qgcanvas_request_update(*d, QG_LEGACY_VIEW_REFRESH_VIEW);
update();
emit changed();
    }

    QOpenGLWidget::keyPressEvent(k);
}

void QgGL::mousePressEvent(QMouseEvent *e)
{

    if (qg_legacy_view_framebuffer_standalone_get(d->ifp) &&
	    e->button() == Qt::RightButton) {
if (window())
    window()->close();
return;
    }

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
qgcanvas_request_update(*d, QG_LEGACY_VIEW_REFRESH_VIEW);
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
qgcanvas_request_update(*d, QG_LEGACY_VIEW_REFRESH_VIEW);
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
qgcanvas_request_update(*d, QG_LEGACY_VIEW_REFRESH_VIEW);
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
qgcanvas_request_update(*d, QG_LEGACY_VIEW_REFRESH_VIEW);
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

    if (qgcanvas_obol_scene_has_content(*d)) {
	qgcanvas_get_obol_viewport_image(*d, this, img, true);
	if (!img.isNull())
	    return;
    }

    if (!d->m_init || !d->dmp)
	return;
    img = grabFramebuffer();
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
