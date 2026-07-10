/*              Q G O B O L W I N D O W H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "qtcad/QgObolWindowHost.h"

#include "brlobol/view_controller.h"
#include "qtcad/QgCanvasBase.h"
#include "qtcad/QgSW.h"

#include <QByteArray>
#include <QImage>
#include <QSize>
#include <QString>
#include <QWidget>

#include <cmath>

struct QgObolWindowHostPrivate {
    QgObolWindowHostPrivate(void) :
	canvas(NULL),
	ownsCanvas(FALSE),
	pollRate(0),
	renderCount(0),
	lastRenderReason("")
    {
    }

    QgCanvasBase *canvas;
    SbBool ownsCanvas;
    long pollRate;
    QImage lastFrame;
    int renderCount;
    SbString lastRenderReason;
};

static QSize
qg_obol_window_host_render_size(const QWidget *widget)
{
    if (!widget)
	return QSize(1, 1);

    qreal dpr = widget->devicePixelRatioF();
    int width = static_cast<int>(std::ceil(widget->width() * dpr));
    int height = static_cast<int>(std::ceil(widget->height() * dpr));
    return QSize(width > 0 ? width : 1, height > 0 ? height : 1);
}

static BRLObolWindowDesc
qg_obol_window_host_desc(const BRLObolWindowDesc *input, const QWidget *widget)
{
    BRLObolWindowDesc desc;
    if (input) {
	desc = *input;
    } else {
	desc.mode = BRLOBOL_WINDOW_EMBEDDED;
	desc.backend = BRLOBOL_WINDOW_BACKEND_QT;
	desc.width = 0;
	desc.height = 0;
	desc.title = "BRL-CAD Obol Qt";
	desc.display = "";
	desc.nativeIdHint = "";
	desc.visible = FALSE;
    }

    QSize renderSize = qg_obol_window_host_render_size(widget);
    if (desc.width == 0)
	desc.width = static_cast<unsigned int>(renderSize.width());
    if (desc.height == 0)
	desc.height = static_cast<unsigned int>(renderSize.height());

    desc.backend = BRLOBOL_WINDOW_BACKEND_QT;
    if (!input)
	desc.mode = BRLOBOL_WINDOW_EMBEDDED;
    if (widget) {
	desc.visible = widget->isVisible() ? TRUE : FALSE;
	if (desc.title.getLength() == 0) {
	    QByteArray title = widget->windowTitle().toUtf8();
	    desc.title = title.isEmpty() ? "BRL-CAD Obol Qt" : title.constData();
	}
    }
    return desc;
}

static int
qg_obol_window_host_can_create_canvas(const BRLObolWindowDesc *desc)
{
    if (!desc)
	return 0;
    if (desc->mode != BRLOBOL_WINDOW_TOPLEVEL)
	return 0;
    return desc->backend == BRLOBOL_WINDOW_BACKEND_QT ||
	desc->backend == BRLOBOL_WINDOW_BACKEND_AUTO;
}

static QgCanvasBase *
qg_obol_window_host_create_canvas(const BRLObolWindowDesc *desc)
{
    if (!qg_obol_window_host_can_create_canvas(desc))
	return NULL;

    QgSW *canvas = new QgSW(NULL);
    QWidget *widget = canvas->canvasWidget();
    if (desc->width || desc->height) {
	int width = desc->width ? (int)desc->width : widget->width();
	int height = desc->height ? (int)desc->height : widget->height();
	widget->resize(width > 0 ? width : 1, height > 0 ? height : 1);
    }
    if (desc->title.getLength() > 0)
	widget->setWindowTitle(QString::fromUtf8(desc->title.getString()));
    if (desc->visible)
	widget->show();

    return canvas;
}

static void
qg_obol_window_host_destroy_owned_canvas(QgObolWindowHostPrivate *qp)
{
    if (!qp || !qp->ownsCanvas || !qp->canvas)
	return;

    QgCanvasBase *canvas = qp->canvas;
    qp->canvas = NULL;
    qp->ownsCanvas = FALSE;

    QWidget *widget = canvas->canvasWidget();
    if (widget)
	widget->close();
    delete canvas;
}

QgObolWindowHost::QgObolWindowHost(QgCanvasBase *canvas) :
    BRLObolWindowHost(),
    qp(new QgObolWindowHostPrivate)
{
    this->setCanvas(canvas);
}

QgObolWindowHost::~QgObolWindowHost(void)
{
    this->close();
    delete this->qp;
    this->qp = NULL;
}

void
QgObolWindowHost::setCanvas(QgCanvasBase *canvas)
{
    if (this->qp->canvas == canvas)
	return;

    if (this->qp->ownsCanvas) {
	qg_obol_window_host_destroy_owned_canvas(this->qp);
    }

    this->qp->canvas = canvas;
    this->qp->ownsCanvas = FALSE;
    this->attachController(canvas ? canvas->obolViewController() : NULL, FALSE);
}

QgCanvasBase *
QgObolWindowHost::canvas(void) const
{
    return this->qp->canvas;
}

int
QgObolWindowHost::open(const BRLObolWindowDesc *desc)
{
    if (!this->qp->canvas) {
	QgCanvasBase *ownedCanvas = qg_obol_window_host_create_canvas(desc);
	if (ownedCanvas) {
	    this->qp->canvas = ownedCanvas;
	    this->qp->ownsCanvas = TRUE;
	}
    }

    if (!this->qp->canvas || !this->qp->canvas->obolViewController())
	return -1;

    SbBool ownsCanvas = this->qp->ownsCanvas;
    this->qp->ownsCanvas = FALSE;
    this->attachController(this->qp->canvas->obolViewController(), FALSE);
    this->qp->ownsCanvas = ownsCanvas;
    BRLObolWindowDesc actual =
	qg_obol_window_host_desc(desc, this->qp->canvas->canvasWidget());
    int ret = BRLObolWindowHost::open(&actual);
    if (ret != 0 && this->qp->ownsCanvas) {
	qg_obol_window_host_destroy_owned_canvas(this->qp);
	this->attachController(NULL, FALSE);
    }
    return ret;
}

void
QgObolWindowHost::close(void)
{
    BRLObolWindowHost::close();
    this->qp->lastFrame = QImage();
    this->qp->lastRenderReason = "";
    if (this->qp->ownsCanvas) {
	qg_obol_window_host_destroy_owned_canvas(this->qp);
	this->attachController(NULL, FALSE);
    }
}

int
QgObolWindowHost::poll(const BRLObolInputProfile *UNUSED(profile))
{
    if (!this->qp->canvas || !this->isOpen())
	return -1;

    BRLObolViewController *controller = this->getController();
    if (!controller)
	return -1;
    if (!controller->isRenderRequested())
	return 0;

    SbString reason = controller->getRenderReason();
    QImage image;
    this->qp->canvas->get_obol_viewport_image(image);
    if (image.isNull())
	return -1;

    controller->consumeRenderRequest(NULL);
    this->qp->lastRenderReason = reason;
    this->qp->lastFrame = image.copy();
    this->qp->renderCount++;

    QWidget *widget = this->qp->canvas->canvasWidget();
    if (widget)
	widget->update();
    return 1;
}

long
QgObolWindowHost::pollRate(void) const
{
    return this->qp->pollRate;
}

void
QgObolWindowHost::setPollRate(long rate)
{
    this->qp->pollRate = rate > 0 ? rate : 0;
}

const QImage &
QgObolWindowHost::lastFrameImage(void) const
{
    return this->qp->lastFrame;
}

int
QgObolWindowHost::lastFrameWidth(void) const
{
    return this->qp->lastFrame.width();
}

int
QgObolWindowHost::lastFrameHeight(void) const
{
    return this->qp->lastFrame.height();
}

int
QgObolWindowHost::renderCount(void) const
{
    return this->qp->renderCount;
}

const SbString &
QgObolWindowHost::lastRenderReason(void) const
{
    return this->qp->lastRenderReason;
}
