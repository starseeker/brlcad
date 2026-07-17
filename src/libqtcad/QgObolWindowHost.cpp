/*              Q G O B O L W I N D O W H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "qtcad/QgObolWindowHost.h"

#include "BObol/BDisplayEndpoint.h"
#include "BObol/BViewController.h"
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
#ifdef BRLCAD_OPENGL
#  include <QOpenGLContext>
#  include <QSurfaceFormat>
#endif
#include <QSize>
#include <QString>
#include <QThread>
#include <QWidget>

#include <cmath>

struct QgObolWindowHostPrivate {
    QgObolWindowHostPrivate(void) :
	canvas(NULL),
	ownsCanvas(FALSE),
	preserveCanvasOnBind(FALSE),
	pollRate(0),
	renderCount(0),
	lastRenderReason("")
    {
    }

    QgCanvasBase *canvas;
    SbBool ownsCanvas;
    SbBool preserveCanvasOnBind;
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

static QSize
qg_obol_window_host_logical_size(const QWidget *widget, unsigned int width,
	unsigned int height)
{
    const qreal dpr = widget && widget->devicePixelRatioF() > 0.0 ?
	widget->devicePixelRatioF() : 1.0;
    return QSize(qMax(1, static_cast<int>(std::ceil(width / dpr))),
	qMax(1, static_cast<int>(std::ceil(height / dpr))));
}

static BObolWindowDesc
qg_obol_window_host_desc(const BObolWindowDesc *input, const QWidget *widget)
{
    BObolWindowDesc desc;
    if (input) {
	desc = *input;
    } else {
	desc.mode = BOBOL_WINDOW_EMBEDDED;
	desc.backend = BOBOL_WINDOW_BACKEND_QT;
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

    desc.backend = BOBOL_WINDOW_BACKEND_QT;
    if (!input)
	desc.mode = BOBOL_WINDOW_EMBEDDED;
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
qg_obol_window_host_can_create_canvas(const BObolWindowDesc *desc)
{
    if (!desc)
	return 0;
    if (desc->mode != BOBOL_WINDOW_TOPLEVEL)
	return 0;
    return desc->backend == BOBOL_WINDOW_BACKEND_QT ||
	desc->backend == BOBOL_WINDOW_BACKEND_AUTO;
}

static QgCanvasBase *
qg_obol_window_host_create_canvas(const BObolWindowDesc *desc)
{
    if (!qg_obol_window_host_can_create_canvas(desc))
	return NULL;

    QgSW *canvas = new QgSW(NULL);
    QWidget *widget = canvas->canvasWidget();
    if (desc->width || desc->height) {
	const QSize renderSize = qg_obol_window_host_render_size(widget);
	const QSize logicalSize = qg_obol_window_host_logical_size(widget,
	    desc->width ? desc->width : (unsigned int)renderSize.width(),
	    desc->height ? desc->height : (unsigned int)renderSize.height());
	widget->resize(logicalSize);
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
    BObolWindowHost(),
    qp(new QgObolWindowHostPrivate)
{
    this->setCanvas(canvas);
    this->qp->ownsCanvas = canvas && takeCanvasOwnership ? TRUE : FALSE;
}

QgObolWindowHost::~QgObolWindowHost(void)
{
    this->close();
    if (this->qp->canvas)
	this->qp->canvas->setObolInputEndpoint(NULL);
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
QgObolWindowHost::bindController(BObolViewController *controller)
{
    if (this->qp->canvas)
	this->qp->canvas->setObolViewController(controller);

    /* BObolWindowHost::attachController closes the current host before it
     * swaps controller ownership.  A factory bind is not a terminal close:
     * retain the already-created Qt canvas so a qt-gl factory cannot reopen
     * as a software QgSW canvas. */
    this->qp->preserveCanvasOnBind = TRUE;
    this->attachController(controller, FALSE);
    this->qp->preserveCanvasOnBind = FALSE;
}

int
QgObolWindowHost::open(const BObolWindowDesc *desc)
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

    BObolViewController *controller = this->getController();
    if (controller && this->qp->canvas->obolViewController() != controller)
	this->qp->canvas->setObolViewController(controller);
    if (!controller)
	this->attachController(this->qp->canvas->obolViewController(), FALSE);
    controller = this->getController();
    if (!controller || !controller->getRenderContextManager())
	return -1;

    QWidget *widget = this->qp->canvas->canvasWidget();
    /* Embedded canvases are sized and shown by their owning application.
     * Re-entering QWidget::show/resize from endpoint or framebuffer setup can
     * synchronously recurse through view and fbserv callbacks. */
    if (widget && desc && this->qp->ownsCanvas) {
	const QSize logicalSize = qg_obol_window_host_logical_size(widget,
	    desc->width ? desc->width : 1,
	    desc->height ? desc->height : 1);
	widget->resize(logicalSize);
	if (desc->title.getLength())
	    widget->setWindowTitle(QString::fromUtf8(desc->title.getString()));
	if (desc->visible)
	    widget->show();
    }
    BObolWindowDesc actual =
	qg_obol_window_host_desc(desc, widget);
    int ret = BObolWindowHost::open(&actual);
    if (ret != 0 && this->qp->ownsCanvas) {
	qg_obol_window_host_destroy_owned_canvas(this->qp);
	this->attachController(NULL, FALSE);
    }
    return ret;
}

void
QgObolWindowHost::close(void)
{
    BObolWindowHost::close();
    this->qp->lastFrame = QImage();
    this->qp->lastRenderReason = "";
    if (this->qp->ownsCanvas && !this->qp->preserveCanvasOnBind) {
	qg_obol_window_host_destroy_owned_canvas(this->qp);
	this->attachController(NULL, FALSE);
    }
}

int
QgObolWindowHost::poll(const BObolInputProfile *UNUSED(profile))
{
    if (!this->qp->canvas || !this->isOpen())
	return -1;

    BObolViewController *controller = this->getController();
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
qg_factory_probe(const struct bobol_host_desc *desc,
	void *UNUSED(user_data))
{
    QCoreApplication *app = QCoreApplication::instance();
    if (!app || QThread::currentThread() != app->thread())
	return 0;
    if (!desc)
	return 0;
    if (desc->mode == BOBOL_HOST_MODE_TOPLEVEL)
	return 1;
    return desc->mode == BOBOL_HOST_MODE_EMBEDDED &&
	desc->application_context;
}

static void *
qg_factory_create(const struct bobol_host_desc *desc, void *data)
{
    QgObolFactoryKind *kind = static_cast<QgObolFactoryKind *>(data);
    if (!kind || !desc)
	return NULL;

    QgCanvasBase *canvas = NULL;
    bool owns_canvas = true;
    if (desc->mode == BOBOL_HOST_MODE_EMBEDDED) {
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

    if (canvas) {
#ifdef BRLCAD_OPENGL
	if (!kind->software && desc->vsync != BOBOL_HOST_VSYNC_AUTO) {
	    QgGL *gl_canvas = dynamic_cast<QgGL *>(canvas);
	    if (!gl_canvas)
		return NULL;
	    QOpenGLContext *context = gl_canvas->context();
	    const bool enabled = desc->vsync != BOBOL_HOST_VSYNC_OFF;
	    if (context &&
		((context->format().swapInterval() != 0) != enabled))
		return NULL;
	    if (!context) {
		QSurfaceFormat format = gl_canvas->format();
		format.setSwapInterval(enabled ? 1 : 0);
		gl_canvas->setFormat(format);
	    }
	}
#endif
	canvas->setObolInputEndpoint(
	    static_cast<bobol_display_endpoint_t *>(desc->input_dispatch_data));
	return new QgObolWindowHost(canvas, owns_canvas);
    }
    if (kind->software) {
	canvas = new QgSW(NULL, nullptr, false);
    } else {
#ifdef BRLCAD_OPENGL
	QgGL *gl_canvas = new QgGL(NULL, nullptr, false);
	if (desc->vsync != BOBOL_HOST_VSYNC_AUTO) {
	    QSurfaceFormat format = gl_canvas->format();
	    format.setSwapInterval(desc->vsync == BOBOL_HOST_VSYNC_OFF ?
		0 : 1);
	    gl_canvas->setFormat(format);
	}
	canvas = gl_canvas;
#else
	return NULL;
#endif
    }
	canvas->setObolInputEndpoint(
	static_cast<bobol_display_endpoint_t *>(desc->input_dispatch_data));
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
    host->bindController(static_cast<BObolViewController *>(controller));
    BObolViewController *bound = host->getController();
    return !controller || (bound == controller &&
	bound->getRenderContextManager() != NULL);
}

static BObolWindowDesc
qg_factory_desc(const struct bobol_host_desc *desc)
{
    BObolWindowDesc actual;
    actual.mode = desc && desc->mode == BOBOL_HOST_MODE_EMBEDDED ?
	BOBOL_WINDOW_EMBEDDED : BOBOL_WINDOW_TOPLEVEL;
    actual.backend = BOBOL_WINDOW_BACKEND_QT;
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
qg_factory_open(void *instance, const struct bobol_host_desc *desc,
	void *UNUSED(user_data))
{
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    BObolWindowDesc actual = qg_factory_desc(desc);
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
qg_factory_dimensions(void *instance, unsigned int *width,
	unsigned int *height, double *device_pixel_ratio,
	void *UNUSED(user_data))
{
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    QWidget *widget = host && host->canvas() ?
	host->canvas()->canvasWidget() : NULL;
    if (!widget || !width || !height || !device_pixel_ratio)
	return 0;

    const QSize render_size = qg_obol_window_host_render_size(widget);
    *width = static_cast<unsigned int>(render_size.width());
    *height = static_cast<unsigned int>(render_size.height());
    const qreal dpr = widget->devicePixelRatioF();
    *device_pixel_ratio = dpr > 0.0 ? static_cast<double>(dpr) : 1.0;
    return 1;
}

static int
qg_factory_set_title(void *instance, const char *title,
	void *UNUSED(user_data))
{
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    QWidget *widget = host && host->canvas() ?
	host->canvas()->canvasWidget() : NULL;
    if (!widget || !title)
	return 0;
    widget->setWindowTitle(QString::fromUtf8(title));
    return 1;
}

static int
qg_factory_set_visible(void *instance, int visible,
	void *UNUSED(user_data))
{
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    QWidget *widget = host && host->canvas() ?
	host->canvas()->canvasWidget() : NULL;
    if (!widget)
	return 0;
    widget->setVisible(visible != 0);
    return 1;
}

static int
qg_factory_set_vsync(void *instance, int enabled,
	void *UNUSED(user_data))
{
#ifdef BRLCAD_OPENGL
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(instance);
    QgGL *canvas = host && host->canvas() ?
	dynamic_cast<QgGL *>(host->canvas()) : NULL;
    if (!canvas)
	return 0;
    QOpenGLContext *context = canvas->context();
    if (context)
	return (context->format().swapInterval() != 0) == (enabled != 0) ?
	    1 : 0;
    QSurfaceFormat format = canvas->format();
    format.setSwapInterval(enabled ? 1 : 0);
    canvas->setFormat(format);
    return 1;
#else
    (void)instance;
    (void)enabled;
    return 0;
#endif
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

static void *
qg_factory_framebuffer_window_host(void *instance, void *UNUSED(user_data))
{
    return static_cast<QgObolWindowHost *>(instance);
}

static bobol_host_factory_token_t *
qg_factory_register(const char *name, int priority, uint64_t capabilities,
	QgObolFactoryKind *kind)
{
    struct bobol_host_factory factory;
    memset(&factory, 0, sizeof(factory));
    factory.abi_version = BOBOL_HOST_FACTORY_ABI_VERSION;
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
    factory.dimensions = qg_factory_dimensions;
    factory.capture = qg_factory_capture;
    factory.set_title = qg_factory_set_title;
    factory.set_visible = qg_factory_set_visible;
    factory.set_vsync = qg_factory_set_vsync;
    factory.framebuffer_window_host = qg_factory_framebuffer_window_host;
    return bobol_host_factory_register(&factory);
}

int
qtcad_obol_host_factories_register(void)
{
    if (!QCoreApplication::instance())
	return 0;

    static QgObolFactoryKind sw_kind = {true};
    static bobol_host_factory_token_t *sw_token = qg_factory_register(
	"qt-sw", 50, BOBOL_HOST_CAP_TOPLEVEL | BOBOL_HOST_CAP_EMBEDDED |
	BOBOL_HOST_CAP_PIXEL_PRESENT |
	BOBOL_HOST_CAP_PROGRESSIVE_PRESENT |
	BOBOL_HOST_CAP_INPUT | BOBOL_HOST_CAP_READBACK |
	BOBOL_HOST_CAP_FRAMEBUFFER_PRESENT |
	BOBOL_HOST_CAP_MULTI_VIEW | BOBOL_HOST_CAP_THREAD_AFFINE,
	&sw_kind);

#ifdef BRLCAD_OPENGL
    static QgObolFactoryKind gl_kind = {false};
    static bobol_host_factory_token_t *gl_token = qg_factory_register(
	"qt-gl", 60, BOBOL_HOST_CAP_TOPLEVEL | BOBOL_HOST_CAP_EMBEDDED |
	BOBOL_HOST_CAP_SYSTEM_GL |
	BOBOL_HOST_CAP_PIXEL_PRESENT |
	BOBOL_HOST_CAP_PRESENT_VSYNC |
	BOBOL_HOST_CAP_PROGRESSIVE_PRESENT |
	BOBOL_HOST_CAP_INPUT | BOBOL_HOST_CAP_READBACK |
	BOBOL_HOST_CAP_FRAMEBUFFER_PRESENT |
	BOBOL_HOST_CAP_MULTI_VIEW | BOBOL_HOST_CAP_THREAD_AFFINE,
	&gl_kind);
    return sw_token && gl_token ? 1 : 0;
#else
    return sw_token ? 1 : 0;
#endif
}
