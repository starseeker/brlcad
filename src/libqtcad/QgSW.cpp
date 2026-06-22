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
 * Qt widget for visualizing libosmesa OpenGL software rasterizer output.
 */

#define USE_MGL_NAMESPACE 1

#include "common.h"

#include <QImage>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QPaintEvent>
#include <QResizeEvent>
#include <QtGlobal>
#include <QWheelEvent>

#include "QgCanvasState.h"   /* pimpl definition + shared helpers */
#include "qtcad/QgSW.h"

// Using the full BSG_VIEW_MIN/BSG_VIEW_MAX was causing drawing artifacts with moss I
// in shaded mode (I think I was seeing the "Z-fighting" problem:
// https://www.sjbaker.org/steve/omniv/love_your_z_buffer.html )
//
// Setting to (-1,1) clips geometry too quickly as we start to zoom in.
// -100,100 seems to work, but may need a better long term solution to
// this... maybe basing it on the currently visible object bounds?
#define QTSW_ZMIN -100
#define QTSW_ZMAX 100
/* Background grey level used when capturing the viewport as an image.
 * Dark but not black, so the yellow wireframe is clearly visible. */
#define QTSW_SCREENSHOT_BG_GREY 40

static thread_local qg_legacy_fb *qgsw_bridge_framebuffer = nullptr;

struct QgSWBridgeFramebufferScope {
    explicit QgSWBridgeFramebufferScope(qg_legacy_fb *fbp)
	: previous(qgsw_bridge_framebuffer)
    {
	qgsw_bridge_framebuffer = fbp;
    }

    ~QgSWBridgeFramebufferScope()
    {
	qgsw_bridge_framebuffer = previous;
    }

    qg_legacy_fb *previous = nullptr;
};

QgSW::QgSW(QWidget *parent)
    : QWidget(parent)
{
    d = new QgCanvasState();
    qgcanvas_init_obol(*d, this);
    d->ifp = qgsw_bridge_framebuffer;
    d->lmouse_mode = QG_LEGACY_VIEW_ADJUST_SCALE;

    // Provide a view specific to this widget - set gedp->ged_gvp to v
    // if this is the current view
    d->local_v = qg_legacy_view_local_create("swrast");
    d->v = d->local_v;
    qgcanvas_sync_obol_camera(*d);

    // Don't dm_open until we have the view.
    d->dmp = nullptr;

    // If we weren't supplied with a framebuffer, allocate one.
    // We don't open it until we have the dmp.
    if (!d->ifp)
d->ifp = qg_legacy_view_framebuffer_raw_create("swrast");

    // This is an important Qt setting for interactivity - it allowing key
    // bindings to propagate to this widget and trigger actions such as
    // resolution scaling, rotation, etc.
    setFocusPolicy(Qt::WheelFocus);
}

QgSW *
QgCanvasBridgeFactory::create_swrast(QWidget *parent, qg_legacy_fb *fbp)
{
    QgSWBridgeFramebufferScope scope(fbp);
    return new QgSW(parent);
}

QgSW::~QgSW()
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
QgSW::view() const
{
    return d->v;
}

bool
QgSW::legacyBackendInitialized() const
{
    return d->dmp != nullptr;
}

BRLObolViewController *
QgSW::obolViewController() const
{
    return d->obol;
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

void
QgSW::set_view(qg_legacy_view *nv)
{
    qgcanvas_set_view(*d, nv);
}

void QgSW::request_update(uint32_t refresh_flags)
{
    uint32_t requested = refresh_flags ? refresh_flags : QG_LEGACY_VIEW_REFRESH_ALL;
    qgcanvas_request_update(*d, requested | QG_LEGACY_VIEW_REFRESH_FRAMEBUFFER | QG_LEGACY_VIEW_REFRESH_FORCE);
    if (d->fb_update_queued)
return;
    d->fb_update_queued = true;
    QMetaObject::invokeMethod(this, "queued_update", Qt::QueuedConnection);
}

void QgSW::need_update()
{
    QTCAD_SLOT("QgSW::need_update", 1);
    request_update(QG_LEGACY_VIEW_REFRESH_FRAMEBUFFER | QG_LEGACY_VIEW_REFRESH_FORCE);
}

void QgSW::queued_update()
{
    d->fb_update_queued = false;
    update();
}

void QgSW::paintEvent(QPaintEvent *e)
{
    // Without a view, SWrast can't work
    if (!d->v)
return;

    if (qgcanvas_obol_scene_has_content(*d)) {
	QImage obolImage;
	qgcanvas_get_obol_viewport_image(*d, this, obolImage, true);
	if (!obolImage.isNull()) {
	    QPainter painter(this);
	    painter.drawImage(QPoint(0, 0), obolImage);
	    QWidget::paintEvent(e);
	    return;
	}
    }

    if (!d->m_init) {

if (!d->dmp) {
    // swrast will need to know the window size
    QSize rsize = qgcanvas_render_size(this);
    qg_legacy_view_dimensions_set(d->v, rsize.width(), rsize.height());

    d->dmp = qg_legacy_view_dm_open_swrast(d->v, this);
    if (!d->dmp)
return;

    // If we have a framebuffer, now we can open it
    (void)qg_legacy_view_dm_framebuffer_setup_existing(d->ifp, d->dmp);
}

(void)qg_legacy_view_dm_setup_swrast(d->dmp, QTSW_ZMIN, QTSW_ZMAX);

qg_legacy_view_dm_bind(d->v, d->dmp);
qg_legacy_view_dm_sync_dimensions(d->v, d->dmp);

// Ready to go
d->m_init = true;

emit init_done();
    }

    if (!d->m_init || !d->dmp)
return;

    // Go ahead and set the flag, but (unlike the rendering thread
    // implementation) we need to do the draw routine every time in paintGL, or
    // we end up with unrendered frames.
    qg_legacy_view_dm_native_repaint_pending_set(d->dmp, 0);

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

    qg_legacy_view_dm_background_restore(d->dmp);

    (void)qg_legacy_view_refresh_consume(d->v);
    qg_legacy_view_dm_draw_begin(d->dmp);
    qg_legacy_view_dm_draw(d->v);
    qg_legacy_view_dm_draw_end(d->dmp);
    qg_legacy_view_refresh_complete(d->v);

    // Set up a QImage with the rendered output..
    unsigned char *dm_image;
    if (qg_legacy_view_dm_display_image_get(d->dmp, &dm_image, 0, 1)) {
return;
    }
    QImage image(dm_image, qg_legacy_view_dm_width_get(d->dmp),
	    qg_legacy_view_dm_height_get(d->dmp), QImage::Format_RGBX8888);
    image.setDevicePixelRatio(devicePixelRatioF());
    QPainter painter(this);
    painter.translate(0, height());
    painter.scale(1.0, -1.0);
    painter.drawImage(QPoint(0, 0), image);
    if (qgcanvas_obol_scene_has_content(*d)) {
	QImage obolImage;
	qgcanvas_get_obol_viewport_image(*d, this, obolImage, true);
	if (!obolImage.isNull()) {
	    painter.resetTransform();
	    painter.drawImage(QPoint(0, 0), obolImage);
	}
    }
    QWidget::paintEvent(e);
}

void QgSW::resizeEvent(QResizeEvent *e)
{
    QWidget::resizeEvent(e);
    qgcanvas_sync_obol_viewport(*d, this);
    if (d->dmp && d->v) {
	QSize rsize = qgcanvas_render_size(this);
	qg_legacy_view_dm_dimensions_set(d->dmp, rsize.width(), rsize.height());
	qg_legacy_view_dimensions_set(d->v, rsize.width(), rsize.height());
	qg_legacy_view_dm_configure_window(d->dmp, 0);
	if (d->ifp) {
	    qg_legacy_view_framebuffer_configure(d->ifp,
		    qg_legacy_view_width_get(d->v),
		    qg_legacy_view_height_get(d->v));
	}
	qgcanvas_request_update(*d, QG_LEGACY_VIEW_REFRESH_VIEW | QG_LEGACY_VIEW_REFRESH_FRAMEBUFFER);
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
    qg_legacy_view_dimensions_set(d->v, rsize.width(), rsize.height());

    if (d->input.keyPressEvent(d->v, d->x_prev,
	    d->y_prev, k)) {
qgcanvas_request_update(*d, QG_LEGACY_VIEW_REFRESH_VIEW);
update();
emit changed();
    }

    QWidget::keyPressEvent(k);
}

void QgSW::mousePressEvent(QMouseEvent *e)
{

    if (qg_legacy_view_framebuffer_standalone_get(d->ifp) &&
	    e->button() == Qt::RightButton) {
if (window())
    window()->close();
return;
    }

    if (!d->v || !d->current || !d->use_default_mousebindings) {
QWidget::mousePressEvent(e);
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

    QWidget::mousePressEvent(e);
}

void QgSW::mouseReleaseEvent(QMouseEvent *e)
{
    if (!d->v) {
QWidget::mouseReleaseEvent(e);
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
    qg_legacy_view_dimensions_set(d->v, rsize.width(), rsize.height());

    if (d->input.wheelEvent(d->v, e)) {
qgcanvas_request_update(*d, QG_LEGACY_VIEW_REFRESH_VIEW);
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
    // Set up a QImage with the rendered output..
    unsigned char *dm_image;
    if (qg_legacy_view_dm_display_image_get(d->dmp, &dm_image, 0, 1)) {
return;
    }
    QImage image(dm_image, qg_legacy_view_dm_width_get(d->dmp),
	    qg_legacy_view_dm_height_get(d->dmp), QImage::Format_RGBX8888);
#if QT_VERSION >= QT_VERSION_CHECK(6, 9, 0)
    image.flipped(Qt::Vertical).save("file.png");
#else
    image.mirrored(false, true).save("file.png");
#endif
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

/* Render the current view and return the pixel data.  Obol content is captured
 * through Obol readback; legacy DM capture remains the fallback. */
void QgSW::get_viewport_image(QImage &img)
{
    img = QImage();  /* null sentinel */
    if (!d->v) return;

    if (qgcanvas_obol_scene_has_content(*d)) {
	qgcanvas_get_obol_viewport_image(*d, this, img, true);
	if (!img.isNull())
	    return;
    }

    /* Ensure DM is initialised (reuse render_to_file init logic) */
    if (!d->m_init) {
if (!d->dmp) {
    int rw = (width()  > 50) ? width()  : 800;
    int rh = (height() > 50) ? height() : 600;
    qg_legacy_view_dimensions_set(d->v, rw, rh);
    d->dmp = qg_legacy_view_dm_open_swrast(d->v, this);
    if (!d->dmp) return;
}
(void)qg_legacy_view_dm_setup_swrast(d->dmp, QTSW_ZMIN, QTSW_ZMAX);
qg_legacy_view_dm_bind(d->v, d->dmp);
qg_legacy_view_dm_sync_dimensions(d->v, d->dmp);
d->m_init = true;
    }
    if (!d->dmp) return;

    /* Render */
    unsigned char bg1[3] = {0, 0, 0};
    unsigned char bg2[3] = {0, 0, 0};
    (void)qg_legacy_view_dm_background_get(d->dmp, bg1, bg2);
    /* Use a dark-grey background for better visibility in screenshots;
     * fall through to the stored background if it is already non-black. */
    if (bg1[0] == 0 && bg1[1] == 0 && bg1[2] == 0 &&
    bg2[0] == 0 && bg2[1] == 0 && bg2[2] == 0) {
/* Default black: override with a neutral dark background */
bg1[0] = bg1[1] = bg1[2] = QTSW_SCREENSHOT_BG_GREY;
bg2[0] = bg2[1] = bg2[2] = QTSW_SCREENSHOT_BG_GREY;
    }
    qg_legacy_view_dm_background_set(d->dmp, bg1, bg2);
    qg_legacy_view_dm_load_current_model2view(d->dmp, d->v, 0);
    (void)qg_legacy_view_refresh_consume(d->v);
    qg_legacy_view_dm_draw_begin(d->dmp);
    qg_legacy_view_dm_draw(d->v);
    qg_legacy_view_dm_draw_end(d->dmp);
    qg_legacy_view_refresh_complete(d->v);

    unsigned char *vp_image = nullptr;
    if (qg_legacy_view_dm_display_image_get(d->dmp, &vp_image, 1, 1) ||
	    !vp_image)
	return;
    /* Copy pixel data into a QImage (QImage doesn't own vp_image) */
    img = QImage(vp_image, qg_legacy_view_dm_width_get(d->dmp),
	    qg_legacy_view_dm_height_get(d->dmp),
 QImage::Format_RGBA8888).copy();
    qg_legacy_view_dm_display_image_release(vp_image);
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
