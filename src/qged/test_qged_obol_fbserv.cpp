/*              T E S T _ Q G E D _ O B O L _ F B S E R V . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file test_qged_obol_fbserv.cpp
 *
 * Exercise qged's Obol/imgstream framebuffer bridge through the same fbserv
 * backend operation table used by ert framebuffer packets.
 */

#include "common.h"

#include "fbserv.h"

#include "brlobol/display_endpoint.h"
#include "brlobol/image_source.h"
#include "brlobol/view_controller.h"
#include "brlobol/viewport_image.h"
#include "bu/app.h"
#include "bu/log.h"
#include "imgstream/fbserv.h"
#include "ged.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "qtcad/QgCanvasBase.h"
#include "qtcad/QgObolWindowHost.h"
#include "qtcad/QgView.h"

#include <Inventor/nodes/SoGroup.h>

#include <QApplication>
#include <QImage>
#include <QWidget>

#include <cmath>
#include <string.h>
#include <vector>

#define FAIL(_msg) \
    do { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

#define CHECK(_expr, _msg) \
    do { \
	if (!(_expr)) \
	    FAIL(_msg); \
    } while (0)

static int
lit_pixel_count(const QImage &image)
{
    QImage rgba = image.convertToFormat(QImage::Format_RGBA8888);
    int count = 0;
    for (int y = 0; y < rgba.height(); y++) {
	const unsigned char *line = rgba.constScanLine(y);
	for (int x = 0; x < rgba.width(); x++) {
	    const unsigned char *p = line + x * 4;
	    if (p[0] > 32 || p[1] > 32 || p[2] > 32)
		count++;
	}
    }
    return count;
}

static int
orange_pixel_count(const QImage &image)
{
    QImage rgba = image.convertToFormat(QImage::Format_RGBA8888);
    int count = 0;
    for (int y = 0; y < rgba.height(); y++) {
	const unsigned char *line = rgba.constScanLine(y);
	for (int x = 0; x < rgba.width(); x++) {
	    const unsigned char *p = line + x * 4;
	    if (p[0] > 200 && p[1] > 32 && p[1] < 100 && p[2] < 64)
		count++;
	}
    }
    return count;
}

static void
find_framebuffer_nodes(SoNode *node, SoBRLImageSource **source,
		       SoBRLViewportImage **viewport)
{
    if (!node || ((source == NULL || *source) &&
		  (viewport == NULL || *viewport))) {
	return;
    }

    if (source && !*source &&
	node->isOfType(SoBRLImageSource::getClassTypeId()))
	*source = static_cast<SoBRLImageSource *>(node);
    if (viewport && !*viewport &&
	node->isOfType(SoBRLViewportImage::getClassTypeId()))
	*viewport = static_cast<SoBRLViewportImage *>(node);

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	find_framebuffer_nodes(group->getChild(i), source, viewport);
}

static int
test_qged_obol_fbserv_backend(void)
{
    QgView view(NULL, QgViewType::SW);
    view.resize(160, 120);
    view.show();
    QApplication::processEvents();

    BRLObolViewController *controller = view.obolViewController();
    CHECK(controller != NULL, "qged test view must expose an Obol controller");

    struct ged *gedp = ged_open("inmem", "0x0", 0);
    CHECK(gedp != NULL, "qged fbserv test must create an in-memory GED");

#define GED_CHECK(_expr, _msg) \
    do { \
	if (!(_expr)) { \
	    bu_log("FAIL: %s\n", _msg); \
	    qged_fbserv_release_ged_handlers(gedp); \
	    ged_close(gedp); \
	    return 1; \
	} \
    } while (0)

    ged_view_active_ctx_set(gedp, view.viewContext());
    GED_CHECK(ged_view_context_host_attach(gedp, view.viewContext()),
	"qged framebuffer GED view must attach its endpoint context");

    qged_fbserv_configure_ged_handlers(gedp, &view);
    void *view_ctx = ged_view_active_ctx(gedp);
    GED_CHECK(ged_draw_obol_controller_opaque_for_view(view_ctx) == controller,
	"qged framebuffer GED view must share the visible endpoint controller");
    GED_CHECK(view.displayEndpoint() != NULL &&
	ged_view_context_display_endpoint_get(view_ctx) == view.displayEndpoint(),
	"qged framebuffer uses the visible view display endpoint");
    GED_CHECK(brlobol_display_endpoint_host(view.displayEndpoint()) != NULL,
	"qged visible endpoint retains its Qt host");
    struct fbserv_obj *fbs = gedp->ged_fbs;
    GED_CHECK(fbs != NULL, "GED must own an fbserv object");
    GED_CHECK(fbs_framebuffer_backend_installed(fbs),
	      "qged must install an Obol fbserv backend");
    GED_CHECK(fbs_can_open_ipc(fbs),
	      "qged must install an IPC fbserv client handler for ert");

    struct fbserv_fb_info info;
    GED_CHECK(fbs_framebuffer_info(fbs, &info) == 0,
	      "qged Obol fbserv backend must report framebuffer dimensions");
    GED_CHECK(info.width > 0 && info.height > 0,
	      "qged Obol framebuffer dimensions must be positive");

    view.resize(300, 220);
    QApplication::processEvents();
    qged_fbserv_configure_ged_handlers(gedp, &view);
    GED_CHECK(fbs_framebuffer_info(fbs, &info) == 0,
	      "qged must report resized framebuffer dimensions");
    QWidget *canvasWidget = view.canvasBase()->canvasWidget();
    int expectedWidth = (int)std::ceil(canvasWidget->width() *
	canvasWidget->devicePixelRatioF());
    int expectedHeight = (int)std::ceil(canvasWidget->height() *
	canvasWidget->devicePixelRatioF());
    GED_CHECK(info.width == expectedWidth && info.height == expectedHeight,
	      "qged framebuffer must track the current physical canvas size");
    SbVec2s viewportSize = controller->getViewportRegion().getViewportSizePixels();
    GED_CHECK(viewportSize[0] == expectedWidth &&
	      viewportSize[1] == expectedHeight,
	      "qged Obol controller must track the physical canvas size");

    struct brlobol_endpoint_property_value composition =
	BRLOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    composition.type = BRLOBOL_ENDPOINT_PROPERTY_ENUM;
    composition.string_value = "overlay";
    GED_CHECK(brlobol_display_endpoint_property_set(view.displayEndpoint(),
	      "composition.framebuffer.mode", &composition) ==
	      BRLOBOL_ENDPOINT_PROPERTY_OK,
	      "qged framebuffer endpoint must enable overlay composition");

    SoBRLImageSource *source = NULL;
    SoBRLViewportImage *viewport = NULL;
    /* Framebuffers are retained presentation layers, not model-scene nodes. */
    find_framebuffer_nodes(controller->getRenderRoot(), NULL, &viewport);
    if (viewport)
	source = viewport->getImageSource();
    GED_CHECK(source != NULL && viewport != NULL,
	      "qged Obol fbserv backend must create framebuffer scene nodes");
    GED_CHECK(viewport->getImageSource() == source,
	      "qged framebuffer viewport must reference its image source");
    GED_CHECK(NEAR_EQUAL(viewport->size.getValue()[0],
		(float)expectedWidth, SMALL_FASTF) &&
	      NEAR_EQUAL(viewport->size.getValue()[1],
		(float)expectedHeight, SMALL_FASTF),
	      "qged framebuffer viewport must fill the resized canvas");
    GED_CHECK(source->hasPendingStreamUpdate() == FALSE,
	      "fresh qged framebuffer source must start without pending pixels");

    std::vector<unsigned char> pixels((size_t)info.width *
	(size_t)info.height * 3);
    for (size_t i = 0; i < pixels.size(); i += 3) {
	pixels[i] = 255;
	pixels[i + 1] = 64;
	pixels[i + 2] = 32;
    }
    int firstBandHeight = info.height / 2;
    GED_CHECK(fbs_framebuffer_writerect(fbs, 0, 0, info.width, firstBandHeight,
		pixels.data()) == info.width * firstBandHeight,
	      "qged Obol fbserv backend must accept a partial scan band");
    GED_CHECK(source->hasPendingStreamUpdate() == TRUE,
	      "qged framebuffer writes must mark source data pending only");
    GED_CHECK(source->dirtyRevision.getValue() == 0,
	      "qged framebuffer writes must not mutate Coin dirty fields");

    GED_CHECK(fbs_framebuffer_view(fbs, 3, 4, 2, 2) == 0,
	      "qged Obol fbserv backend must record framebuffer view state");
    GED_CHECK(fbs_framebuffer_cursor(fbs, 1, 7, 8) == 0,
	      "qged Obol fbserv backend must record cursor state");
    GED_CHECK(fbs_framebuffer_flush(fbs) == 0,
	      "qged Obol fbserv backend must present pending stream state");

    GED_CHECK(source->hasPendingStreamUpdate() == FALSE,
	      "qged framebuffer flush must consume pending stream updates");
    GED_CHECK(source->dirtyRevision.getValue() > 0,
	      "qged framebuffer flush must publish dirty generation to Coin");
    GED_CHECK(viewport->realizedDirtyRevision.getValue() ==
	      source->dirtyRevision.getValue(),
	      "qged framebuffer viewport must sync to source dirty generation");
    GED_CHECK(NEAR_EQUAL(viewport->sourceCenter.getValue()[0], 3.0f, SMALL_FASTF) &&
	      NEAR_EQUAL(viewport->sourceCenter.getValue()[1], 4.0f, SMALL_FASTF) &&
	      NEAR_EQUAL(viewport->sourceZoom.getValue(), 2.0f, SMALL_FASTF),
	      "qged framebuffer flush must publish view state to viewport image");
    GED_CHECK(viewport->cursorVisible.getValue() == TRUE &&
	      NEAR_EQUAL(viewport->cursorImagePosition.getValue()[0], 7.0f, SMALL_FASTF) &&
	      NEAR_EQUAL(viewport->cursorImagePosition.getValue()[1], 8.0f, SMALL_FASTF),
	      "qged framebuffer flush must publish cursor state to viewport image");
    GED_CHECK(controller->isRenderRequested(),
	      "qged framebuffer flush must request an Obol render");
    GED_CHECK(fbs_framebuffer_poll(fbs) == 1,
	      "qged framebuffer poll must drain the requested Obol render");
    GED_CHECK(!controller->isRenderRequested(),
	      "qged framebuffer poll must consume the render request");

    QImage intermediate;
    view.get_obol_viewport_image(intermediate);
    int intermediateOrange = orange_pixel_count(intermediate);
    GED_CHECK(intermediateOrange > intermediate.width() * intermediate.height() / 4 &&
	intermediateOrange < intermediate.width() * intermediate.height() * 3 / 4,
	      "qged must visibly present a meaningful partial framebuffer update");

    size_t secondBandOffset = (size_t)info.width * (size_t)firstBandHeight * 3;
    int secondBandHeight = info.height - firstBandHeight;
    GED_CHECK(fbs_framebuffer_writerect(fbs, 0, firstBandHeight,
		info.width, secondBandHeight, pixels.data() + secondBandOffset) ==
	      info.width * secondBandHeight,
	      "qged Obol fbserv backend must accept the final scan band");
    GED_CHECK(fbs_framebuffer_flush(fbs) == 0 &&
	      fbs_framebuffer_poll(fbs) == 1,
	      "qged must present the completed framebuffer update");

    QImage image;
    view.get_obol_viewport_image(image);
    GED_CHECK(!image.isNull() && lit_pixel_count(image) >
	      image.width() * image.height() * 8 / 10,
	      "qged Obol framebuffer content must fill the view");

    /* A GED session has one framebuffer stream, but QGED may move it to a
     * different endpoint when the active quad pane changes.  Libged must move
     * the retained image node with the stream rather than leave it attached to
     * the old controller. */
    QgView secondView(NULL, QgViewType::SW);
    secondView.resize(120, 90);
    secondView.show();
    QApplication::processEvents();
    BRLObolViewController *secondController = secondView.obolViewController();
    GED_CHECK(secondController != NULL,
	"qged secondary view must expose an Obol controller");
    ged_view_active_ctx_set(gedp, secondView.viewContext());
    GED_CHECK(ged_view_context_host_attach(gedp, secondView.viewContext()),
	"qged secondary framebuffer view must attach its endpoint context");
    qged_fbserv_configure_ged_handlers(gedp, &secondView);
    GED_CHECK(ged_draw_obol_controller_opaque_for_view(
	secondView.viewContext()) == secondController,
	"qged framebuffer stream must adopt the secondary endpoint controller");
    GED_CHECK(fbs_framebuffer_info(fbs, &info) == 0,
	"qged framebuffer stream must remain valid after an active-view switch");
    QWidget *secondCanvasWidget = secondView.canvasBase()->canvasWidget();
    int secondExpectedWidth = (int)std::ceil(secondCanvasWidget->width() *
	secondCanvasWidget->devicePixelRatioF());
    int secondExpectedHeight = (int)std::ceil(secondCanvasWidget->height() *
	secondCanvasWidget->devicePixelRatioF());
    GED_CHECK(info.width == secondExpectedWidth &&
	info.height == secondExpectedHeight,
	"qged framebuffer stream must adopt the secondary physical canvas size");

    composition.string_value = "overlay";
    GED_CHECK(brlobol_display_endpoint_property_set(
	secondView.displayEndpoint(), "composition.framebuffer.mode",
	&composition) == BRLOBOL_ENDPOINT_PROPERTY_OK,
	"qged framebuffer switch must enable composition on the new endpoint");
    QgObolWindowHost *secondHost = static_cast<QgObolWindowHost *>(
	brlobol_display_endpoint_host(secondView.displayEndpoint()));
    GED_CHECK(secondHost && secondHost->getFramebufferCount() == 1,
	"qged framebuffer switch must bind the stream to the new endpoint host");
    GED_CHECK(secondHost->getController() == secondController &&
	secondController->getFramebufferOverlayRoot()->getNumChildren() == 1,
	"qged framebuffer switch must retain its image in the new controller overlay");

    SoBRLImageSource *firstSource = NULL;
    SoBRLViewportImage *firstViewport = NULL;
    find_framebuffer_nodes(controller->getRenderRoot(), &firstSource,
	&firstViewport);
    GED_CHECK(firstSource == NULL && firstViewport == NULL,
	"qged framebuffer switch must detach retained nodes from the old endpoint");
    SoBRLImageSource *secondSource = NULL;
    SoBRLViewportImage *secondViewport = NULL;
    find_framebuffer_nodes(secondController->getRenderRoot(), NULL,
	&secondViewport);
    if (secondViewport)
	secondSource = secondViewport->getImageSource();
    GED_CHECK(secondSource != NULL && secondViewport != NULL,
	"qged framebuffer switch must attach retained nodes to the new endpoint");

    std::vector<unsigned char> secondPixels((size_t)info.width *
	(size_t)info.height * 3);
    for (size_t i = 0; i < secondPixels.size(); i += 3) {
	secondPixels[i] = 255;
	secondPixels[i + 1] = 64;
	secondPixels[i + 2] = 32;
    }
    GED_CHECK(fbs_framebuffer_writerect(fbs, 0, 0, info.width, info.height,
		secondPixels.data()) == info.width * info.height &&
	fbs_framebuffer_flush(fbs) == 0 && fbs_framebuffer_poll(fbs) == 1,
	"qged framebuffer switch must publish new pixels through the new endpoint");
    QImage switchedImage;
    secondView.get_obol_viewport_image(switchedImage);
    GED_CHECK(!switchedImage.isNull() &&
	orange_pixel_count(switchedImage) >
	switchedImage.width() * switchedImage.height() * 8 / 10,
	"qged framebuffer switch must visibly present the new endpoint frame");

    qged_fbserv_release_ged_handlers(gedp);
    ged_close(gedp);
#undef GED_CHECK
    return 0;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    QApplication app(argc, argv);

    return test_qged_obol_fbserv_backend();
}
