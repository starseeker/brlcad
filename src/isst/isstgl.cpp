/*                      I S S T G L . C P P
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

#include "common.h"

#include <QKeyEvent>
#include <QPainter>
#include <QResizeEvent>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

#include "bu/parallel.h"
#include "isstgl.h"


TIERenderer::TIERenderer()
{
    camera.type = RENDER_CAMERA_PERSPECTIVE;
    camera.fov = 25;
    render_camera_init(&camera, bu_avail_cpus());
    render_phong_init(&camera.render, NULL);

    TIENET_BUFFER_INIT(buffer_image);
    tile.orig_x = 0;
    tile.orig_y = 0;
    tile.format = RENDER_CAMERA_BIT_DEPTH_24;
}

TIERenderer::~TIERenderer()
{
    TIENET_BUFFER_FREE(buffer_image);
}

bool
TIERenderer::resize()
{
    size_t render_width;
    size_t render_height;

    if (viewport_width <= 0 || viewport_height <= 0)
	return false;

    if (resolution_factor == 0) {
	render_width = (size_t)viewport_width;
	render_height = (size_t)viewport_height;
    } else {
	render_width = (size_t)resolution_factor;
	if (render_width > SIZE_MAX / (size_t)viewport_height)
	    return false;
	render_height = render_width * (size_t)viewport_height /
	    (size_t)viewport_width;
    }

    if (!render_width || !render_height || render_width > UINT16_MAX ||
	render_height > UINT16_MAX || render_width > SIZE_MAX / render_height ||
	render_width * render_height >
	(std::numeric_limits<uint32_t>::max() - sizeof(camera_tile_t)) / 3u)
	return false;

    camera.w = tile.size_x = (uint16_t)render_width;
    camera.h = tile.size_y = (uint16_t)render_height;
    const size_t render_bytes = render_width * render_height * 3u +
	sizeof(camera_tile_t);
    TIENET_BUFFER_SIZE(buffer_image, (uint32_t)render_bytes);
    return true;
}

void
TIERenderer::res_incr()
{
    setResolution(resolution + 1);
}

void
TIERenderer::res_decr()
{
    setResolution(resolution - 1);
}

void
TIERenderer::setResolution(int in_resolution)
{
    resolution = std::max(1, std::min(in_resolution, 20));
    resolution_factor = (resolution == 20) ? 0 :
	(int)std::lround((double)viewport_width * .05 * resolution);
    changed = true;
    render();
}

void
TIERenderer::setSize(int width, int height)
{
    if (viewport_width == width && viewport_height == height)
	return;

    viewport_width = width;
    viewport_height = height;
    if (resolution != 20)
	resolution_factor = (int)std::lround((double)viewport_width * .05 * resolution);
    changed = true;
    render();
}

void
TIERenderer::setTie(struct tie_s *in_tie)
{
    tie = in_tie;
    if (!tie)
	return;

    VSETALL(camera.pos, tie->radius);
    VMOVE(camera.focus, tie->mid);
    VSETALL(camera_pos_init, tie->radius);
    VMOVE(camera_focus_init, tie->mid);
    changed = true;
    render();
}

void
TIERenderer::clearTie()
{
    tie = NULL;
    changed = false;
}

void
TIERenderer::render()
{
    if (m_exiting.load() || !changed || !tie)
	return;
    if (!resize())
	return;

    changed = false;
    buffer_image.ind = 0;
    render_camera_prep(&camera);
    render_camera_render(&camera, tie, &tile, &buffer_image);
    if (m_exiting.load())
	return;

    const QImage frame(buffer_image.data + sizeof(camera_tile_t), camera.w,
	camera.h, camera.w * 3, QImage::Format_RGB888);
    emit imageReady(frame.copy());
}

isstView::isstView(QWidget *parent)
    : QWidget(parent)
{
    m_thread = new QThread;
    m_renderer = new TIERenderer();
    m_renderer->moveToThread(m_thread);
    connect(m_thread, &QThread::finished, m_renderer, &QObject::deleteLater);

    connect(this, &isstView::resolutionIncreased, m_renderer,
	&TIERenderer::res_incr);
    connect(this, &isstView::resolutionDecreased, m_renderer,
	&TIERenderer::res_decr);
    connect(this, &isstView::resolutionRequested, m_renderer,
	&TIERenderer::setResolution);
    connect(this, &isstView::sceneChanged, m_renderer, &TIERenderer::setTie);
    connect(this, &isstView::sizeChanged, m_renderer, &TIERenderer::setSize);
    connect(m_renderer, &TIERenderer::imageReady, this, &isstView::setImage);

    m_thread->start();
    setFocusPolicy(Qt::WheelFocus);
}

isstView::~isstView()
{
    m_renderer->prepareExit();
    m_thread->quit();
    m_thread->wait();
    delete m_thread;
}

void
isstView::set_tie(struct tie_s *in_tie)
{
    emit sceneChanged(in_tie);
}

void
isstView::clear_tie()
{
    QMetaObject::invokeMethod(m_renderer, "clearTie", Qt::BlockingQueuedConnection);
}

void
isstView::setPreviewResolution(int resolution)
{
    emit resolutionRequested(resolution);
}

void
isstView::setImage(const QImage &image)
{
    m_image = image;
    update();
    emit imagePresented();
}

void
isstView::paintEvent(QPaintEvent *event)
{
    QWidget::paintEvent(event);
    QPainter painter(this);
    painter.fillRect(rect(), Qt::black);
    if (!m_image.isNull()) {
	painter.setRenderHint(QPainter::SmoothPixmapTransform, false);
	painter.drawImage(rect(), m_image);
    }
}

void
isstView::resizeEvent(QResizeEvent *event)
{
    QWidget::resizeEvent(event);
    emit sizeChanged(event->size().width(), event->size().height());
}

void
isstView::keyPressEvent(QKeyEvent *key)
{
    switch (key->key()) {
	case Qt::Key_Equal:
	    emit resolutionIncreased();
	    return;
	case Qt::Key_Minus:
	    emit resolutionDecreased();
	    return;
	default:
	    break;
    }
    QWidget::keyPressEvent(key);
}

void
isstView::save_image()
{
    if (!m_image.isNull())
	m_image.save("file.png");
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
