/*        T E S T _ Q T C A D _ O B O L _ W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "qtcad/QgObolWindowHost.h"
#include "qtcad/QgSW.h"

#include "brlobol.h"
#include "bu/log.h"
#include "dm.h"

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
#include <QImage>
#include <QWidget>

#include <string.h>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)

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
    brlobol_window_host_unregister_display_host();
    CHECK(brlobol_window_host_register_display_host(&host) == 0,
	"qtcad window host registers as display host");

    imgstream_fb_t *fb = imgstream_fb_open("/dev/qtgl", 8, 6);
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
    brlobol_window_host_unregister_display_host();
    return 0;
}

static int
test_qtcad_legacy_fb_open_display_host_bridge(void)
{
    if (dm_valid_type("qtgl", NULL) != 1) {
	bu_log("SKIP: qtgl display manager plugin is unavailable\n");
	return 0;
    }

    QgSW canvas;
    canvas.resize(64, 48);

    QgObolWindowHost host(&canvas);
    brlobol_window_host_unregister_display_host();
    CHECK(brlobol_window_host_register_display_host(&host) == 0,
	"qtcad window host registers for legacy fb_open bridge");

    struct fb *fbp = fb_open("/dev/qtgl", 9, 7);
    CHECK(fbp != NULL, "legacy fb_open uses registered Obol display host");
    CHECK(fb_get_dm(fbp) == NULL,
	"legacy fb_open host bridge does not create nested libdm backend");
    CHECK(host.getFramebufferCount() == 1,
	"legacy fb_open attaches one Obol framebuffer");

    unsigned char rgb[3] = {11, 22, 33};
    CHECK(fb_writerect(fbp, 1, 1, 1, 1, rgb) == 1,
	"legacy fb_open bridge accepts RGB rectangle writes");
    unsigned char readback[3] = {0, 0, 0};
    CHECK(fb_read(fbp, 1, 1, readback, 1) == 1 &&
	    readback[0] == 11 && readback[1] == 22 && readback[2] == 33,
	"legacy fb_open bridge reads from imgstream storage");
    CHECK(fb_view(fbp, 4, 3, 2, 2) == 0,
	"legacy fb_open bridge forwards view state");
    CHECK(fb_cursor(fbp, 1, 5, 4) == 0,
	"legacy fb_open bridge forwards cursor state");
    CHECK(fb_flush(fbp) == 0,
	"legacy fb_open bridge flushes Obol viewport image");
    CHECK(fb_poll(fbp) == 1,
	"legacy fb_open bridge drains qtcad Obol readback");
    CHECK(!host.lastFrameImage().isNull(),
	"legacy fb_open bridge records qtcad readback image");

    CHECK(fb_close(fbp) == 0,
	"legacy fb_open bridge closes cleanly");
    CHECK(host.getFramebufferCount() == 0,
	"legacy fb_open bridge detaches Obol framebuffer on close");
    brlobol_window_host_unregister_display_host();
    return 0;
}

static int
test_qtcad_owned_display_host_bridge(void)
{
    QgObolWindowHost host;
    brlobol_window_host_unregister_display_host();
    CHECK(brlobol_window_host_register_display_host(&host) == 0,
	"owned qtcad window host registers as display host");

    imgstream_fb_t *fb = imgstream_fb_open("/dev/qtgl", 10, 7);
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
    brlobol_window_host_unregister_display_host();
    return 0;
}

int
main(int argc, char **argv)
{
    QApplication app(argc, argv);
    (void)app;

    if (test_qtcad_window_host_contract())
	return 1;
    if (test_qtcad_owned_window_host())
	return 1;
    if (test_qtcad_display_host_bridge())
	return 1;
    if (test_qtcad_legacy_fb_open_display_host_bridge())
	return 1;
    if (test_qtcad_owned_display_host_bridge())
	return 1;
    return 0;
}
