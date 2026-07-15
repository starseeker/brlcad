/*        T E S T _ Q T C A D _ O B O L _ W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bv.h"

#include "qtcad/QgObolWindowHost.h"
#include "qtcad/QgQuadView.h"
#include "qtcad/QgSession.h"
#ifdef BRLCAD_OPENGL
#  include "qtcad/QgGL.h"
#endif
#include "qtcad/QgSW.h"
#include "qtcad/QgView.h"

#include "brlobol.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "bu/file.h"
#include "ged.h"
#include "ged/view.h"
#include "wdb.h"

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wshadow"
#endif
#include "imgstream/fb_compat.h"
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/SoViewport.h>

#include <QApplication>
#include <QElapsedTimer>
#include <QImage>
#ifdef BRLCAD_OPENGL
#  include <QOpenGLContext>
#  include <QOpenGLWidget>
#endif
#include <QThread>
#include <QWidget>

#include <string.h>
#include <thread>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)

class ThreadAffinityCanvas : public QgSW {
public:
    explicit ThreadAffinityCanvas(QWidget *parent = nullptr) : QgSW(parent) {}

    void need_update() override
    {
	need_update_calls++;
	update_on_application_thread &=
	    QThread::currentThread() == QCoreApplication::instance()->thread();
	QgSW::need_update();
    }

    int need_update_calls = 0;
    bool update_on_application_thread = true;
};

static void
add_visible_obol_content(BRLObolViewController *controller)
{
    if (!controller || !controller->getSceneRoot() ||
	    !controller->getSceneRoot()->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *root = static_cast<SoGroup *>(controller->getSceneRoot());
    root->addChild(new SoDirectionalLight);
    SoMaterial *material = new SoMaterial;
    material->emissiveColor.setValue(0.9f, 0.3f, 0.1f);
    material->diffuseColor.setValue(0.9f, 0.3f, 0.1f);
    root->addChild(material);
    SoCube *cube = new SoCube;
    cube->width = 80.0f;
    cube->height = 80.0f;
    cube->depth = 80.0f;
    root->addChild(cube);
    controller->getViewport()->viewAll();
}

static int
lit_pixel_count(const QImage &image)
{
    QImage rgba = image.convertToFormat(QImage::Format_RGBA8888);
    int count = 0;
    for (int y = 0; y < rgba.height(); y++) {
	const unsigned char *line = rgba.constScanLine(y);
	for (int x = 0; x < rgba.width(); x++) {
	    const unsigned char *p = line + x * 4;
	    if (p[0] > 16 || p[1] > 16 || p[2] > 16)
		count++;
	}
    }
    return count;
}

static BRLObolWindowDesc
qt_toplevel_desc(unsigned int width, unsigned int height, SbBool visible)
{
    BRLObolWindowDesc desc;
    desc.mode = BRLOBOL_WINDOW_TOPLEVEL;
    desc.backend = BRLOBOL_WINDOW_BACKEND_QT;
    desc.width = width;
    desc.height = height;
    desc.title = "BRL-CAD Obol Qt Test";
    desc.display = "";
    desc.nativeIdHint = "";
    desc.visible = visible;
    return desc;
}

static int
test_qtcad_window_host_contract(void)
{
    QgSW canvas;
    canvas.resize(96, 72);
    add_visible_obol_content(canvas.obolViewController());

    QgObolWindowHost host(&canvas);
    CHECK(host.open() == 0, "qtcad window host opens against canvas");
    CHECK(host.canvas() == &canvas, "qtcad window host records canvas");
    CHECK(host.getController() == canvas.obolViewController(),
	"qtcad window host attaches canvas Obol controller");
    CHECK(host.getDesc().backend == BRLOBOL_WINDOW_BACKEND_QT,
	"qtcad window host forces Qt backend");
    CHECK(host.getDesc().mode == BRLOBOL_WINDOW_EMBEDDED,
	"qtcad window host defaults to embedded mode");
    CHECK(host.getDesc().width >= 96 && host.getDesc().height >= 72,
	"qtcad window host derives size from canvas");

    host.setPollRate(30);
    CHECK(host.pollRate() == 30, "qtcad window host stores poll rate");
    host.setPollRate(-1);
    CHECK(host.pollRate() == 0, "qtcad window host clamps negative poll rate");

    canvas.obolViewController()->requestRender("qtcad-window-host-render");
    CHECK(host.poll(NULL) == 1, "qtcad window host drains pending render");
    CHECK(!canvas.obolViewController()->isRenderRequested(),
	"qtcad window host consumes render request after readback");
    CHECK(host.renderCount() == 1, "qtcad window host records render count");
    CHECK(strcmp(host.lastRenderReason().getString(), "qtcad-window-host-render") == 0,
	"qtcad window host records render reason");
    CHECK(!host.lastFrameImage().isNull() && lit_pixel_count(host.lastFrameImage()) > 0,
	"qtcad window host records visible readback pixels");
    CHECK(host.poll(NULL) == 0, "qtcad window host poll is idle without pending render");
    return 0;
}

static int
test_qtcad_owned_window_host(void)
{
    BRLObolWindowDesc desc = qt_toplevel_desc(88, 66, TRUE);
    QgObolWindowHost host;
    CHECK(host.open(&desc) == 0, "qtcad window host creates owned Qt canvas");
    CHECK(host.canvas() != NULL, "owned qtcad window host records created canvas");
    CHECK(host.getController() == host.canvas()->obolViewController(),
	"owned qtcad window host attaches created canvas controller");
    CHECK(host.getDesc().mode == BRLOBOL_WINDOW_TOPLEVEL &&
	    host.getDesc().backend == BRLOBOL_WINDOW_BACKEND_QT,
	"owned qtcad window host preserves Qt toplevel descriptor");
    CHECK(host.getDesc().width >= 88 && host.getDesc().height >= 66,
	"owned qtcad window host applies requested size");

    QWidget *widget = host.canvas()->canvasWidget();
    CHECK(widget != NULL && widget->isWindow(),
	"owned qtcad window host creates toplevel widget");
    QApplication::processEvents();
    CHECK(widget->isVisible(), "owned qtcad window host shows requested visible widget");

    add_visible_obol_content(host.getController());
    host.getController()->requestRender("owned-qtcad-window-host-render");
    CHECK(host.poll(NULL) == 1, "owned qtcad window host drains pending render");
    CHECK(!host.lastFrameImage().isNull() && lit_pixel_count(host.lastFrameImage()) > 0,
	"owned qtcad window host records visible readback pixels");

    host.close();
    CHECK(host.canvas() == NULL, "owned qtcad window host releases owned canvas on close");
    return 0;
}

static int
test_qtcad_display_host_bridge(void)
{
    QgSW canvas;
    canvas.resize(64, 48);

    QgObolWindowHost host(&canvas);
    imgstream_fb_t *fb = brlobol_window_host_open_display_framebuffer(&host, "/dev/qtgl", 8, 6);
    CHECK(fb != NULL, "qtgl framebuffer opens through qtcad window host");
    CHECK(host.isOpen(), "qtgl framebuffer opens qtcad host");
    CHECK(host.getFramebufferCount() == 1,
	"qtcad host records framebuffer attachment");
    CHECK(host.getFramebufferViewportImage(fb) != NULL,
	"qtcad host creates framebuffer viewport image");
    CHECK(host.getDesc().backend == BRLOBOL_WINDOW_BACKEND_QT &&
	    host.getDesc().mode == BRLOBOL_WINDOW_TOPLEVEL,
	"qtgl framebuffer preserves display-capable Qt descriptor");

    unsigned char rgb[3] = {255, 255, 255};
    CHECK(imgstream_fb_writerect(fb, 0, 0, 1, 1, rgb) == 1,
	"qtgl framebuffer stream accepts writes");
    CHECK(imgstream_fb_flush(fb) == 0,
	"qtgl framebuffer flush syncs Obol viewport image");
    CHECK(imgstream_fb_poll(fb) == 1,
	"qtgl framebuffer poll drains qtcad Obol readback");
    CHECK(!host.lastFrameImage().isNull(),
	"qtcad display host records readback image");
    CHECK(strcmp(host.lastRenderReason().getString(), "fb-flush") == 0,
	"qtcad display host records framebuffer render reason");

    imgstream_fb_close(fb);
    CHECK(host.getFramebufferCount() == 0,
	"qtcad display host closes framebuffer attachment");
    return 0;
}

static int
test_qtcad_owned_display_host_bridge(void)
{
    QgObolWindowHost host;
    imgstream_fb_t *fb = brlobol_window_host_open_display_framebuffer(&host, "/dev/qtgl", 10, 7);
    CHECK(fb != NULL, "qtgl framebuffer creates owned qtcad host");
    CHECK(host.isOpen(), "owned qtgl framebuffer opens qtcad host");
    CHECK(host.canvas() != NULL, "owned qtgl framebuffer creates qtcad canvas");
    CHECK(host.canvas()->canvasWidget()->isWindow(),
	"owned qtgl framebuffer canvas is a toplevel widget");
    CHECK(host.getDesc().backend == BRLOBOL_WINDOW_BACKEND_QT &&
	    host.getDesc().mode == BRLOBOL_WINDOW_TOPLEVEL &&
	    host.getDesc().visible == TRUE,
	"owned qtgl framebuffer preserves visible Qt descriptor");
    CHECK(host.getFramebufferCount() == 1,
	"owned qtcad host records framebuffer attachment");

    unsigned char rgb[3] = {64, 192, 255};
    CHECK(imgstream_fb_writerect(fb, 0, 0, 1, 1, rgb) == 1,
	"owned qtgl framebuffer stream accepts writes");
    CHECK(imgstream_fb_flush(fb) == 0,
	"owned qtgl framebuffer flush syncs Obol viewport image");
    CHECK(imgstream_fb_poll(fb) == 1,
	"owned qtgl framebuffer poll drains qtcad Obol readback");
    CHECK(!host.lastFrameImage().isNull(),
	"owned qtcad display host records readback image");

    imgstream_fb_close(fb);
    CHECK(host.getFramebufferCount() == 0,
	"owned qtcad display host closes framebuffer attachment");
    host.close();
    CHECK(host.canvas() == NULL, "owned display host releases qtcad canvas on close");
    return 0;
}

static int
test_qtcad_factory_endpoint(void)
{
    CHECK(qtcad_obol_host_factories_register(),
	"qtcad registers Obol host factories");

    struct brlobol_host_desc desc;
    memset(&desc, 0, sizeof(desc));
    desc.struct_size = sizeof(desc);
    desc.mode = BRLOBOL_HOST_MODE_TOPLEVEL;
    desc.width = 80;
    desc.height = 60;
    desc.device_pixel_ratio = 1.0;
    desc.visible = 0;
    desc.required_capabilities = BRLOBOL_HOST_CAP_PIXEL_PRESENT |
	BRLOBOL_HOST_CAP_READBACK;

    brlobol_display_endpoint_t *endpoint =
	brlobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL, "qtcad factory test creates endpoint");
    CHECK(brlobol_display_endpoint_host_open(endpoint, "qt-sw", &desc),
	"endpoint opens Qt software factory");
    CHECK(strcmp(brlobol_display_endpoint_host_factory_name(endpoint),
	"qt-sw") == 0,
	"endpoint reports Qt software factory identity");

    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(
	brlobol_display_endpoint_host(endpoint));
    BRLObolViewController *controller =
	static_cast<BRLObolViewController *>(
	    brlobol_display_endpoint_controller(endpoint));
    CHECK(host && host->canvas() &&
	host->canvas()->obolViewController() == controller &&
	host->getController() == controller,
	"Qt factory canvas borrows the endpoint controller");

    struct brlobol_endpoint_property_value property =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    property.type = BRLOBOL_ENDPOINT_PROPERTY_STRING;
    property.string_value = "Qt typed endpoint title";
    CHECK(brlobol_display_endpoint_property_set(endpoint, "endpoint.title",
	&property) == BRLOBOL_ENDPOINT_PROPERTY_OK &&
	host->canvas()->canvasWidget()->windowTitle() ==
	    QStringLiteral("Qt typed endpoint title"),
	"Qt toplevel factory applies the typed title property");
    property = BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    property.type = BRLOBOL_ENDPOINT_PROPERTY_BOOL;
    property.bool_value = 1;
    CHECK(brlobol_display_endpoint_property_set(endpoint, "endpoint.visible",
	&property) == BRLOBOL_ENDPOINT_PROPERTY_OK &&
	host->canvas()->canvasWidget()->isVisible(),
	"Qt toplevel factory applies the typed visibility property");

    add_visible_obol_content(controller);
    CHECK(brlobol_display_endpoint_request_frame(endpoint, "qt-factory-test"),
	"Qt factory queues endpoint frame request");
    QApplication::processEvents();

    unsigned char *pixels = NULL;
    size_t size = 0;
    unsigned int width = 0;
    unsigned int height = 0;
    unsigned int components = 0;
    CHECK(brlobol_display_endpoint_capture(endpoint, &pixels, &size, &width,
	&height, &components),
	"Qt software factory captures through endpoint API");
    CHECK(pixels && width == 80 && height == 60 && components == 4 &&
	size == 80u * 60u * 4u,
	"Qt software factory capture reports packed endpoint dimensions");
    bu_free(pixels, "Qt software factory capture");

    property.bool_value = 0;
    CHECK(brlobol_display_endpoint_property_set(endpoint, "endpoint.visible",
	&property) == BRLOBOL_ENDPOINT_PROPERTY_OK &&
	!host->canvas()->canvasWidget()->isVisible(),
	"Qt typed visibility property hides the toplevel host");

    brlobol_display_endpoint_destroy(endpoint);
    return 0;
}

static int
test_qtcad_embedded_factory_endpoint(void)
{
    QWidget parent;
    parent.resize(120, 90);
    QgSW canvas(&parent);
    canvas.resize(72, 54);

    struct brlobol_host_desc desc;
    memset(&desc, 0, sizeof(desc));
    desc.struct_size = sizeof(desc);
    desc.mode = BRLOBOL_HOST_MODE_EMBEDDED;
    desc.width = 80;
    desc.height = 60;
    desc.device_pixel_ratio = 1.0;
    desc.visible = 1;
    desc.required_capabilities = BRLOBOL_HOST_CAP_EMBEDDED |
	BRLOBOL_HOST_CAP_PIXEL_PRESENT | BRLOBOL_HOST_CAP_READBACK;
    desc.application_context = static_cast<QgCanvasBase *>(&canvas);

    brlobol_display_endpoint_t *endpoint =
	brlobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL, "embedded Qt factory test creates endpoint");
    CHECK(brlobol_display_endpoint_host_open(endpoint, "qt-sw", &desc),
	"endpoint opens Qt software factory on existing canvas");

    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(
	brlobol_display_endpoint_host(endpoint));
    BRLObolViewController *controller =
	static_cast<BRLObolViewController *>(
	    brlobol_display_endpoint_controller(endpoint));
    CHECK(host && host->canvas() == &canvas &&
	canvas.obolViewController() == controller,
	"embedded Qt factory borrows canvas and endpoint controller");
    CHECK(!canvas.isWindow() && canvas.parentWidget() == &parent,
	"embedded Qt factory preserves canvas parent and window mode");
    CHECK(!canvas.isVisible() && canvas.size() == QSize(72, 54),
	"embedded Qt host does not show or resize its borrowed canvas");
    CHECK(host->getDesc().mode == BRLOBOL_WINDOW_EMBEDDED,
	"embedded Qt host preserves endpoint mode");

    add_visible_obol_content(controller);
    unsigned char *pixels = NULL;
    size_t size = 0;
    unsigned int width = 0;
    unsigned int height = 0;
    unsigned int components = 0;
    CHECK(brlobol_display_endpoint_capture(endpoint, &pixels, &size, &width,
	&height, &components),
	"embedded Qt software host captures through endpoint");
    CHECK(pixels && width == 72 && height == 54 && components == 4,
	"embedded Qt capture follows the application-owned canvas dimensions");
    bu_free(pixels, "embedded Qt factory capture");

    brlobol_display_endpoint_destroy(endpoint);
    CHECK(canvas.obolViewController() == NULL,
	"embedded Qt factory unbinds canvas before controller destruction");
    CHECK(canvas.parentWidget() == &parent,
	"embedded Qt endpoint teardown does not destroy borrowed canvas");

    desc.application_context = NULL;
    endpoint = brlobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL,
	"missing embedded Qt canvas test creates endpoint");
    CHECK(!brlobol_display_endpoint_host_open(endpoint, "qt-sw", &desc),
	"embedded Qt factory rejects a missing borrowed canvas");
    brlobol_display_endpoint_destroy(endpoint);
    return 0;
}

static int
test_qtcad_factory_frame_thread_affinity(void)
{
    QWidget parent;
    ThreadAffinityCanvas canvas(&parent);
    canvas.resize(72, 54);

    struct brlobol_host_desc desc;
    memset(&desc, 0, sizeof(desc));
    desc.struct_size = sizeof(desc);
    desc.mode = BRLOBOL_HOST_MODE_EMBEDDED;
    desc.width = 72;
    desc.height = 54;
    desc.device_pixel_ratio = 1.0;
    desc.required_capabilities = BRLOBOL_HOST_CAP_EMBEDDED |
	BRLOBOL_HOST_CAP_PIXEL_PRESENT;
    desc.application_context = static_cast<QgCanvasBase *>(&canvas);

    brlobol_display_endpoint_t *endpoint =
	brlobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL, "Qt affinity test creates endpoint");
    CHECK(brlobol_display_endpoint_host_open(endpoint, "qt-sw", &desc),
	"Qt affinity test opens embedded software factory");

    QApplication::processEvents();
    const int calls_before = canvas.need_update_calls;
    canvas.update_on_application_thread = true;
    int request_ok = 0;
    std::thread worker([&]() {
	request_ok = brlobol_display_endpoint_request_frame(endpoint,
	    "worker-thread-frame");
    });
    worker.join();

    QElapsedTimer timer;
    timer.start();
    while (canvas.need_update_calls == calls_before && timer.elapsed() < 2000)
	QApplication::processEvents();
    CHECK(request_ok && canvas.need_update_calls > calls_before &&
	canvas.update_on_application_thread,
	"Qt factory queues worker frame requests on the application thread");

    brlobol_display_endpoint_destroy(endpoint);
    return 0;
}

#ifdef BRLCAD_OPENGL
static int
test_qtcad_system_gl_factory_endpoint(void)
{
    const QString platform_name = QGuiApplication::platformName();
    if (platform_name == QStringLiteral("offscreen") ||
	platform_name == QStringLiteral("minimal"))
	return 0;

    CHECK(qtcad_obol_host_factories_register(),
	"Qt system-GL test registers Obol host factories");

    struct brlobol_host_desc desc = {};
    desc.struct_size = sizeof(desc);
    desc.mode = BRLOBOL_HOST_MODE_TOPLEVEL;
    desc.width = 96;
    desc.height = 72;
    desc.device_pixel_ratio = 1.0;
    desc.visible = 1;
    desc.required_capabilities = BRLOBOL_HOST_CAP_SYSTEM_GL |
	BRLOBOL_HOST_CAP_PIXEL_PRESENT | BRLOBOL_HOST_CAP_READBACK;

    brlobol_display_endpoint_t *endpoint =
	brlobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint && brlobol_display_endpoint_render_engine_set(endpoint,
	BRLOBOL_RENDER_ENGINE_HW) &&
	brlobol_display_endpoint_host_open(endpoint, "qt-gl", &desc),
	"Qt system-GL endpoint opens an explicit hardware host");

    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(
	brlobol_display_endpoint_host(endpoint));
    CHECK(host, "Qt system-GL endpoint exposes its Qt host");
    CHECK(host->canvas(), "Qt system-GL endpoint creates a canvas");
    QOpenGLWidget *canvas = qobject_cast<QOpenGLWidget *>(
	host->canvas()->asQObject());
    CHECK(canvas, "Qt system-GL endpoint owns a QOpenGLWidget");

    QElapsedTimer timer;
    timer.start();
	while ((!canvas->isVisible() || !canvas->isValid()) && timer.elapsed() < 2000) {
	QApplication::processEvents();
	QThread::msleep(1);
    }
    CHECK(canvas->isVisible() && canvas->isValid() && canvas->context(),
	"Qt system-GL endpoint creates a visible widget with a valid context");

    canvas->makeCurrent();
    CHECK(QOpenGLContext::currentContext() == canvas->context(),
	"Qt system-GL endpoint makes its widget context current");
    canvas->doneCurrent();

    BRLObolViewController *controller = static_cast<BRLObolViewController *>(
	brlobol_display_endpoint_controller(endpoint));
    add_visible_obol_content(controller);
    CHECK(brlobol_display_endpoint_request_frame(endpoint, "qt-system-gl"),
	"Qt system-GL endpoint queues an Obol frame");

    timer.restart();
    while (controller->isRenderRequested() && timer.elapsed() < 2000) {
	QApplication::processEvents();
	QThread::msleep(1);
    }
    CHECK(!controller->isRenderRequested(),
	"Qt system-GL widget consumes the queued Obol frame");

    QImage presented = canvas->grabFramebuffer();
    CHECK(!presented.isNull() && lit_pixel_count(presented) > 0,
	"Qt system-GL widget presents visible Obol pixels");

    unsigned char *pixels = NULL;
    size_t size = 0;
    unsigned int width = 0;
    unsigned int height = 0;
    unsigned int components = 0;
    CHECK(brlobol_display_endpoint_capture(endpoint, &pixels, &size, &width,
	&height, &components) && pixels && width == 96 && height == 72 &&
	components == 4 && size == 96u * 72u * 4u,
	"Qt system-GL endpoint captures through its active host");
    bu_free(pixels, "Qt system-GL endpoint capture");

    brlobol_display_endpoint_destroy(endpoint);
    return 0;
}

static int
test_qtcad_gl_vsync_policy(void)
{
    const QString platform_name = QGuiApplication::platformName();
    if (platform_name == QStringLiteral("offscreen") ||
	platform_name == QStringLiteral("minimal"))
	return 0;

    struct brlobol_host_desc desc = {};
    desc.struct_size = sizeof(desc);
    desc.mode = BRLOBOL_HOST_MODE_TOPLEVEL;
    desc.width = 32;
    desc.height = 24;
    desc.device_pixel_ratio = 1.0;
    desc.vsync = BRLOBOL_HOST_VSYNC_OFF;
    desc.required_capabilities = BRLOBOL_HOST_CAP_SYSTEM_GL |
	BRLOBOL_HOST_CAP_PRESENT_VSYNC;

    brlobol_display_endpoint_t *endpoint =
	brlobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint && brlobol_display_endpoint_render_engine_set(endpoint,
	BRLOBOL_RENDER_ENGINE_HW) &&
	brlobol_display_endpoint_host_open(endpoint, "qt-gl", &desc),
	"Qt GL factory accepts an explicit construction-time vsync policy");
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(
	brlobol_display_endpoint_host(endpoint));
    QgGL *canvas = host && host->canvas() ?
	dynamic_cast<QgGL *>(host->canvas()) : NULL;
    CHECK(canvas && canvas->format().swapInterval() == 0,
	"Qt GL canvas applies disabled vsync before context creation");

    struct brlobol_endpoint_property_value property =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    property.type = BRLOBOL_ENDPOINT_PROPERTY_BOOL;
    property.bool_value = 1;
    CHECK(brlobol_display_endpoint_property_set(endpoint, "endpoint.vsync",
	&property) == BRLOBOL_ENDPOINT_PROPERTY_OK,
	"Qt GL host applies typed vsync while its context is unrealized");
    property = BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    CHECK(brlobol_display_endpoint_property_get(endpoint, "endpoint.vsync",
	&property) == BRLOBOL_ENDPOINT_PROPERTY_OK && property.bool_value &&
	canvas->format().swapInterval() != 0,
	"Qt GL host reports the applied typed vsync policy");
    brlobol_display_endpoint_destroy(endpoint);
    return 0;
}
#endif

static int
test_qtcad_quad_view_endpoint_association(void)
{
    QgSession session;
    QgQuadView quad(NULL, &session, QgViewType::SW);
    QgView *view = quad.curr_view();
    CHECK(view && view->displayEndpoint(),
	"qtcad quad view creates an endpoint-hosted viewport");
    void *view_ctx = view->viewContext();
    CHECK(ged_view_context_display_endpoint_get(view_ctx) ==
	view->displayEndpoint(),
	"qtcad quad view associates its endpoint with the GED view record");
    CHECK(session.activeViewContext() == view_ctx,
	"qtcad session records the active endpoint view context");
    CHECK(brlobol_display_endpoint_controller(view->displayEndpoint()) ==
	view->obolViewController(),
	"qtcad GED view and visible canvas share one endpoint controller");

    quad.changeToQuadFrame();
    const QgQuadrantId quadrants[] = {
	QgQuadrantId::UpperRight,
	QgQuadrantId::UpperLeft,
	QgQuadrantId::LowerLeft,
	QgQuadrantId::LowerRight
    };
    for (QgQuadrantId quadrant : quadrants) {
	QgView *quadrant_view = quad.get(quadrant);
	CHECK(quadrant_view && quadrant_view->displayEndpoint(),
	    "qtcad lazy quad creation gives every pane an endpoint");
	CHECK(ged_view_context_display_endpoint_get(
	    quadrant_view->viewContext()) == quadrant_view->displayEndpoint(),
	    "qtcad lazy quad pane has one GED-associated endpoint");
	CHECK(brlobol_display_endpoint_controller(
	    quadrant_view->displayEndpoint()) ==
	    quadrant_view->obolViewController(),
	    "qtcad lazy quad pane shares its endpoint controller with its canvas");
    }

    QgView *upper_right = quad.get(QgQuadrantId::UpperRight);
    QgView *upper_left = quad.get(QgQuadrantId::UpperLeft);
    CHECK(upper_right && upper_left && upper_right != upper_left &&
	upper_right->viewContext() != upper_left->viewContext() &&
	upper_right->displayEndpoint() != upper_left->displayEndpoint(),
	"qtcad quad panes retain independent context and endpoint identities");
    quad.select(QgQuadrantId::UpperLeft);
    session.setActiveViewContext(upper_left->viewContext());
    quad.changeToSingleFrame();
    CHECK(quad.curr_view() == upper_right &&
	session.activeViewContext() == upper_right->viewContext(),
	"qtcad quad-to-single switch selects a surviving GED view before teardown");
    CHECK(quad.get(QgQuadrantId::UpperLeft) == upper_right,
	"qtcad quad-to-single switch deletes lazy secondary panes");

    quad.changeToQuadFrame();
    for (QgQuadrantId quadrant : quadrants) {
	QgView *quadrant_view = quad.get(quadrant);
	CHECK(quadrant_view && quadrant_view->displayEndpoint() &&
	    ged_view_context_display_endpoint_get(quadrant_view->viewContext()) ==
	    quadrant_view->displayEndpoint(),
	    "qtcad lazy quad recreation restores one endpoint per pane");
    }
    return 0;
}

static int
test_qtcad_dm_open_command(void)
{
    const char *dbpath = "qtcad_dm_open.g";
    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    CHECK(wdbp != NULL, "Qt dm open test creates a database");
    mk_id(wdbp, "Qt endpoint dm open test");
    wdb_close(wdbp);

    struct ged *gedp = ged_open("db", dbpath, 1);
    CHECK(gedp != NULL, "Qt dm open test opens GED");
    void *view = ged_view_active_ctx(gedp);
    bv_dimensions_set(bv_context_view((struct bv_context *)view), 96, 72);
    CHECK(!ged_view_context_display_endpoint_get(view),
	"Qt dm open test starts without an endpoint");

    const char *open_av[7] = {
	"dm", "open", "--host", "qt-sw", "--renderer", "sw", NULL
    };
    int ret = ged_exec_dm(gedp, 6, open_av);
    brlobol_display_endpoint_t *endpoint =
	ged_view_context_display_endpoint_get(view);
    CHECK(ret == BRLCAD_OK && endpoint &&
	brlobol_display_endpoint_host_factory_name(endpoint) &&
	strcmp(brlobol_display_endpoint_host_factory_name(endpoint),
	    "qt-sw") == 0,
	"dm open creates a Qt software endpoint through the registered factory");
    QgObolWindowHost *host = static_cast<QgObolWindowHost *>(
	brlobol_display_endpoint_host(endpoint));
    CHECK(host && host->canvas() && host->canvas()->canvasWidget()->isWindow(),
	"dm open creates a top-level Qt canvas");
    void *controller = brlobol_display_endpoint_controller(endpoint);

    struct brlobol_endpoint_property_value background =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    background.type = BRLOBOL_ENDPOINT_PROPERTY_COLOR3;
    VSET(background.color3, 0.10, 0.20, 0.30);
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	"controller.background.bottom", &background) ==
	BRLOBOL_ENDPOINT_PROPERTY_OK,
	"Qt dm open test updates the endpoint-owned lower background");
    VSET(background.color3, 0.40, 0.50, 0.60);
    CHECK(brlobol_display_endpoint_property_set(endpoint,
	"controller.background.top", &background) ==
	BRLOBOL_ENDPOINT_PROPERTY_OK,
	"Qt dm open test updates the endpoint-owned upper background");

    const char *close_av[3] = {"dm", "close", NULL};
    CHECK(ged_exec_dm(gedp, 2, close_av) == BRLCAD_OK &&
	ged_view_context_display_endpoint_get(view) == endpoint &&
	brlobol_display_endpoint_controller(endpoint) == controller &&
	!brlobol_display_endpoint_host_factory_name(endpoint),
	"dm close releases the Qt host while retaining the endpoint");
    CHECK(ged_exec_dm(gedp, 6, open_av) == BRLCAD_OK &&
	ged_view_context_display_endpoint_get(view) == endpoint &&
	brlobol_display_endpoint_controller(endpoint) == controller &&
	brlobol_display_endpoint_host_factory_name(endpoint) &&
	strcmp(brlobol_display_endpoint_host_factory_name(endpoint),
	"qt-sw") == 0,
	"dm open reattaches the existing Qt endpoint to a fresh host");
    host = static_cast<QgObolWindowHost *>(
	brlobol_display_endpoint_host(endpoint));
    CHECK(host && host->canvas() && host->canvas()->canvasWidget()->isWindow(),
	"reopened Qt endpoint owns a new top-level canvas");
    BRLObolViewController *endpoint_controller =
	static_cast<BRLObolViewController *>(controller);
    CHECK(endpoint_controller->getBackgroundBottomColor() ==
	SbColor(0.10f, 0.20f, 0.30f) &&
	endpoint_controller->getBackgroundTopColor() ==
	SbColor(0.40f, 0.50f, 0.60f),
	"replacement Qt host preserves endpoint-owned background policy");
    CHECK(ged_exec_dm(gedp, 2, close_av) == BRLCAD_OK &&
	!brlobol_display_endpoint_host_factory_name(endpoint),
	"reopened Qt endpoint closes its replacement host");
    ged_close(gedp);
    bu_file_delete(dbpath);
    return 0;
}

int
main(int argc, char **argv)
{
    QApplication app(argc, argv);
    (void)app;

#ifdef BRLCAD_OPENGL
    if (argc == 2 && strcmp(argv[1], "--system-gl") == 0)
	return test_qtcad_system_gl_factory_endpoint();
#endif

    if (test_qtcad_window_host_contract())
	return 1;
    if (test_qtcad_owned_window_host())
	return 1;
    if (test_qtcad_display_host_bridge())
	return 1;
    if (test_qtcad_owned_display_host_bridge())
	return 1;
    if (test_qtcad_factory_endpoint())
	return 1;
    if (test_qtcad_embedded_factory_endpoint())
	return 1;
    if (test_qtcad_factory_frame_thread_affinity())
	return 1;
#ifdef BRLCAD_OPENGL
    if (test_qtcad_system_gl_factory_endpoint())
	return 1;
    if (test_qtcad_gl_vsync_policy())
	return 1;
#endif
    if (test_qtcad_quad_view_endpoint_association())
	return 1;
    if (test_qtcad_dm_open_command())
	return 1;
    return 0;
}
