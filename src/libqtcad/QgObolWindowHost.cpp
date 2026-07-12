/*              Q G O B O L W I N D O W H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "qtcad/QgObolWindowHost.h"

#include "brlobol/view_controller.h"
#include "bu/malloc.h"
#include "qtcad/QgCanvasBase.h"
#ifdef BRLCAD_OPENGL
#  include "qtcad/QgGL.h"
#endif
#include "qtcad/QgSW.h"

#include <QByteArray>
#include <QCoreApplication>
#include <QImage>
#include <QMetaObject>
#include <QSize>
#include <QString>
#include <QThread>
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

QgObolWindowHost::QgObolWindowHost(QgCanvasBase *canvas,
	bool takeCanvasOwnership) :
    BRLObolWindowHost(),
    qp(new QgObolWindowHostPrivate)
{
    this->setCanvas(canvas);
    this->qp->ownsCanvas = canvas && takeCanvasOwnership ? TRUE : FALSE;
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

void
QgObolWindowHost::bindController(BRLObolViewController *controller)
{
    if (this->qp->canvas)
	this->qp->canvas->setObolViewController(controller);
    this->attachController(controller, FALSE);
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

    BRLObolViewController *controller = this->getController();
    if (controller && this->qp->canvas->obolViewController() != controller)
	this->qp->canvas->setObolViewController(controller);
    if (!controller)
	this->attachController(this->qp->canvas->obolViewController(), FALSE);

    QWidget *widget = this->qp->canvas->canvasWidget();
    /* Embedded canvases are sized and shown by their owning application.
     * Re-entering QWidget::show/resize from endpoint or framebuffer setup can
     * synchronously recurse through view and fbserv callbacks. */
    if (widget && desc && this->qp->ownsCanvas) {
	widget->resize(desc->width ? (int)desc->width : 1,
	    desc->height ? (int)desc->height : 1);
	if (desc->title.getLength())
	    widget->setWindowTitle(QString::fromUtf8(desc->title.getString()));
	if (desc->visible)
	    widget->show();
    }
    BRLObolWindowDesc actual =
	qg_obol_window_host_desc(desc, widget);
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

struct QgObolFactoryKind {
    bool software;
};

static int
qg_factory_probe(const struct brlobol_host_desc *desc,
	void *UNUSED(user_data))
{
    QCoreApplication *app = QCoreApplication::instance();
    if (!app || QThread::currentThread() != app->thread())
	return 0;
    if (!desc)
	return 0;
    if (desc->mode == BRLOBOL_HOST_MODE_TOPLEVEL)
	return 1;
    return desc->mode == BRLOBOL_HOST_MODE_EMBEDDED &&
	desc->application_context;
}

static void *
qg_factory_create(const struct brlobol_host_desc *desc, void *data)
{
    QgObolFactoryKind *kind = static_cast<QgObolFactoryKind *>(data);
    if (!kind || !desc)
	return NULL;

    QgCanvasBase *canvas = NULL;
    bool owns_canvas = true;
    if (desc->mode == BRLOBOL_HOST_MODE_EMBEDDED) {
	canvas = static_cast<QgCanvasBase *>(desc->application_context);
	owns_canvas = false;
	if (!canvas || !canvas->canvasWidget())
	    return NULL;
	if (kind->software && !dynamic_cast<QgSW *>(canvas))
	    return NULL;
#ifdef BRLCAD_OPENGL
	if (!kind->software && !dynamic_cast<QgGL *>(canvas))
	    return NULL;
#else
	if (!kind->software)
	    return NULL;
#endif
    }

    if (canvas)
	return new QgObolWindowHost(canvas, owns_canvas);
    if (kind->software) {
	canvas = new QgSW(NULL);
    } else {
#ifdef BRLCAD_OPENGL
	canvas = new QgGL(NULL);
#else
	return NULL;
#endif
    }
    return new QgObolWindowHost(canvas, owns_canvas);
}

static void
qg_factory_destroy(void *instance, void *UNUSED(user_data))
{
    delete static_cast<QgObolWindowHost *>(instance);
}

static int
qg_factory_bind(void *instance, void *controller, void *UNUSED(user_data))
{
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    if (!host)
	return 0;
    host->bindController(static_cast<BRLObolViewController *>(controller));
    return 1;
}

static BRLObolWindowDesc
qg_factory_desc(const struct brlobol_host_desc *desc)
{
    BRLObolWindowDesc actual;
    actual.mode = desc && desc->mode == BRLOBOL_HOST_MODE_EMBEDDED ?
	BRLOBOL_WINDOW_EMBEDDED : BRLOBOL_WINDOW_TOPLEVEL;
    actual.backend = BRLOBOL_WINDOW_BACKEND_QT;
    actual.width = desc && desc->width ? desc->width : 1;
    actual.height = desc && desc->height ? desc->height : 1;
    actual.title = desc && desc->title ? desc->title : "BRL-CAD Obol Qt";
    actual.display = desc && desc->display ? desc->display : "";
    actual.nativeIdHint = desc && desc->native_id_hint ?
	desc->native_id_hint : "";
    actual.visible = desc && desc->visible ? TRUE : FALSE;
    return actual;
}

static int
qg_factory_open(void *instance, const struct brlobol_host_desc *desc,
	void *UNUSED(user_data))
{
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    BRLObolWindowDesc actual = qg_factory_desc(desc);
    return host && host->open(&actual) == 0 ? 1 : 0;
}

static void
qg_factory_close(void *instance, void *UNUSED(user_data))
{
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    if (host)
	host->close();
}

static int
qg_factory_request_frame(void *instance, const char *UNUSED(reason),
	void *UNUSED(user_data))
{
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    if (!host || !host->canvas() || !host->canvas()->asQObject())
	return 0;
    return QMetaObject::invokeMethod(host->canvas()->asQObject(), "need_update",
	Qt::QueuedConnection) ? 1 : 0;
}

static int
qg_factory_resize(void *instance, unsigned int width, unsigned int height,
	double device_pixel_ratio, void *UNUSED(user_data))
{
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    if (!host || !host->canvas() || !host->canvas()->canvasWidget() ||
	!width || !height || device_pixel_ratio <= 0.0)
	return 0;
    QWidget *widget = host->canvas()->canvasWidget();
    widget->resize((int)std::ceil(width / device_pixel_ratio),
	(int)std::ceil(height / device_pixel_ratio));
    return 1;
}

static int
qg_factory_capture(void *instance, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components,
	void *UNUSED(user_data))
{
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    if (!host || !host->canvas() || !pixels || !size || !width || !height ||
	!components)
	return 0;

    QImage image;
    host->canvas()->get_obol_viewport_image(image);
    if (image.isNull())
	return 0;
    image = image.convertToFormat(QImage::Format_RGBA8888);
    *width = (unsigned int)image.width();
    *height = (unsigned int)image.height();
    *components = 4;
    *size = (size_t)(*width) * (size_t)(*height) * 4;
    *pixels = static_cast<unsigned char *>(bu_malloc(*size,
	"Qt endpoint capture"));
    const size_t row_size = (size_t)(*width) * 4;
    for (unsigned int y = 0; y < *height; y++)
	memcpy(*pixels + (size_t)y * row_size,
	    image.constScanLine((int)(*height - y - 1)),
	    row_size);
    return 1;
}

static brlobol_host_factory_token_t *
qg_factory_register(const char *name, int priority, uint64_t capabilities,
	QgObolFactoryKind *kind)
{
    struct brlobol_host_factory factory;
    memset(&factory, 0, sizeof(factory));
    factory.abi_version = BRLOBOL_HOST_FACTORY_ABI_VERSION;
    factory.struct_size = sizeof(factory);
    factory.name = name;
    factory.priority = priority;
    factory.capabilities = capabilities;
    factory.user_data = kind;
    factory.probe = qg_factory_probe;
    factory.create = qg_factory_create;
    factory.destroy = qg_factory_destroy;
    factory.bind_controller = qg_factory_bind;
    factory.open = qg_factory_open;
    factory.close = qg_factory_close;
    factory.request_frame = qg_factory_request_frame;
    factory.resize = qg_factory_resize;
    factory.capture = qg_factory_capture;
    return brlobol_host_factory_register(&factory);
}

int
qtcad_obol_host_factories_register(void)
{
    if (!QCoreApplication::instance())
	return 0;

    static QgObolFactoryKind sw_kind = {true};
    static brlobol_host_factory_token_t *sw_token = qg_factory_register(
	"qt-sw", 50, BRLOBOL_HOST_CAP_TOPLEVEL | BRLOBOL_HOST_CAP_EMBEDDED |
	BRLOBOL_HOST_CAP_PIXEL_PRESENT |
	BRLOBOL_HOST_CAP_PROGRESSIVE_PRESENT |
	BRLOBOL_HOST_CAP_INPUT | BRLOBOL_HOST_CAP_READBACK |
	BRLOBOL_HOST_CAP_FRAMEBUFFER_PRESENT |
	BRLOBOL_HOST_CAP_MULTI_VIEW | BRLOBOL_HOST_CAP_THREAD_AFFINE,
	&sw_kind);

#ifdef BRLCAD_OPENGL
    static QgObolFactoryKind gl_kind = {false};
    static brlobol_host_factory_token_t *gl_token = qg_factory_register(
	"qt-gl", 60, BRLOBOL_HOST_CAP_TOPLEVEL | BRLOBOL_HOST_CAP_EMBEDDED |
	BRLOBOL_HOST_CAP_SYSTEM_GL |
	BRLOBOL_HOST_CAP_PROGRESSIVE_PRESENT |
	BRLOBOL_HOST_CAP_INPUT | BRLOBOL_HOST_CAP_READBACK |
	BRLOBOL_HOST_CAP_FRAMEBUFFER_PRESENT |
	BRLOBOL_HOST_CAP_MULTI_VIEW | BRLOBOL_HOST_CAP_THREAD_AFFINE,
	&gl_kind);
    return sw_token && gl_token ? 1 : 0;
#else
    return sw_token ? 1 : 0;
#endif
}
