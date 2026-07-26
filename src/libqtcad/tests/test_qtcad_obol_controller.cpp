/*             T E S T _ Q T C A D _ O B O L _ C O N T R O L L E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"

#include "bu/str.h"

#include "bv.h"

#include "BObol/BDisplayEndpoint.h"
#include "BObol/BViewController.h"
#ifdef BRLCAD_OPENGL
#include "qtcad/QgGL.h"
#endif
#include "qtcad/QgSW.h"
#include "qtcad/QgView.h"

#include <Inventor/SoDB.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoMaterial.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoSeparator.h>

#include <QApplication>
#include <QImage>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPainter>
#include <QWheelEvent>

#include <math.h>
#include <stdio.h>
#include <string.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

#ifdef BRLCAD_OPENGL
class TestQgGL : public QgGL {
public:
    explicit TestQgGL(QWidget *parent = NULL) : QgGL(parent) {}
    void runKeyPressForTest(QKeyEvent *event) { this->keyPressEvent(event); }
    void runMouseMoveForTest(QMouseEvent *event) { this->mouseMoveEvent(event); }
    void runPaintGLForTest(void) { this->paintGL(); }
    void runWheelForTest(QWheelEvent *event) { this->wheelEvent(event); }
};
#endif

class TestQgSW : public QgSW {
public:
    explicit TestQgSW(QWidget *parent = NULL) : QgSW(parent) {}
    void runKeyPressForTest(QKeyEvent *event) { this->keyPressEvent(event); }
    void runMouseMoveForTest(QMouseEvent *event) { this->mouseMoveEvent(event); }
    void runWheelForTest(QWheelEvent *event) { this->wheelEvent(event); }
};

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

static QMouseEvent
mouse_move_event(int x, int y, Qt::MouseButtons buttons)
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    return QMouseEvent(QEvent::MouseMove, QPoint(x, y), Qt::NoButton,
	    buttons, Qt::NoModifier);
#else
    return QMouseEvent(QEvent::MouseMove, QPointF(x, y), QPointF(x, y),
	    Qt::NoButton, buttons, Qt::NoModifier);
#endif
}

static QWheelEvent
wheel_event(int x, int y, int angleDeltaY)
{
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    return QWheelEvent(QPointF(x, y), QPointF(x, y), QPoint(),
	    QPoint(0, angleDeltaY), angleDeltaY, Qt::Vertical, Qt::NoButton,
	    Qt::NoModifier);
#else
    return QWheelEvent(QPointF(x, y), QPointF(x, y), QPoint(),
	    QPoint(0, angleDeltaY), Qt::NoButton, Qt::NoModifier,
	    Qt::NoScrollPhase, false);
#endif
}

static int
camera_moved(SoCamera *camera, const SbVec3f &before)
{
    if (!camera)
	return 0;
    SbVec3f after = camera->position.getValue();
    return fabsf(before[0] - after[0]) >= 0.001f ||
	fabsf(before[1] - after[1]) >= 0.001f ||
	fabsf(before[2] - after[2]) >= 0.001f;
}

static int
camera_state_changed(SoCamera *camera, const SbVec3f &beforePosition,
	float beforeFocalDistance)
{
    if (camera_moved(camera, beforePosition))
	return 1;
    return fabsf(camera->focalDistance.getValue() - beforeFocalDistance) >= 0.001f;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    QApplication app(argc, argv);

    QgView view(NULL, QgViewType::SW);
    view.resize(160, 120);
    if (!SoDB::getContextManager())
	FAIL("QgView should install an Obol context manager");

    BObolViewController *controller = view.obolViewController();
    if (!controller)
	FAIL("QgView should expose an Obol view controller");
    if (!view.displayEndpoint() ||
	bobol_display_endpoint_controller(view.displayEndpoint()) != controller)
	FAIL("QgView should expose its endpoint-owned Obol controller");
    if (!bobol_display_endpoint_host(view.displayEndpoint()) ||
	bu_strcmp(bobol_display_endpoint_host_factory_name(view.displayEndpoint()),
	    "qt-sw") != 0)
	FAIL("QgView should host its software canvas through the Qt endpoint factory");
    QWidget *swWidget = view.canvasBase() ? view.canvasBase()->canvasWidget() : NULL;
    if (!swWidget || !swWidget->testAttribute(Qt::WA_OpaquePaintEvent) ||
	swWidget->autoFillBackground())
	FAIL("QgSW should own opaque Obol presentation without parent background painting");

    if (controller != view.obolViewController())
	FAIL("QgView should return a stable Obol view controller");
    if (!controller->getSceneRoot())
	FAIL("QgView Obol controller should have a default scene root");
    if (!controller->getRenderRoot())
	FAIL("QgView Obol controller should expose an Obol render root");
    if (!controller->getCamera())
	FAIL("QgView Obol controller should have a default camera");
    SoCamera *camera = controller->getCamera();
    if (!camera->isOfType(SoOrthographicCamera::getClassTypeId()))
	FAIL("QgView Obol controller should default to an orthographic camera");
    if (!controller->getViewport() ||
	    !controller->getRenderSceneRoot() ||
	    controller->getViewport()->getSceneGraph() != controller->getRenderSceneRoot() ||
	    controller->getViewport()->getCamera() != camera)
	FAIL("QgView Obol controller should wire render scene and camera through SoViewport");
    if (!controller->getRenderManager() ||
	    controller->getRenderManager()->getSceneGraph() != controller->getRenderRoot() ||
	    controller->getRenderManager()->getCamera() != camera)
	FAIL("QgView Obol controller should wire scene and camera through SoRenderManager");
    controller->getRenderManager()->setRenderMode(SoRenderManager::WIREFRAME);
    if (controller->getRenderManager()->getRenderMode() != SoRenderManager::WIREFRAME)
	FAIL("QgView Obol controller should expose native SoRenderManager features");
    if (controller->getViewportRegion().getWindowSize()[0] <= 0 ||
	    controller->getViewportRegion().getWindowSize()[1] <= 0)
	FAIL("Obol view controller should have a valid viewport");
    if (controller->getViewportRegion().getViewportOriginPixels() !=
	    SbVec2s(0, 0) ||
	controller->getViewportRegion().getViewportSizePixels() !=
	    controller->getViewportRegion().getWindowSize())
	FAIL("Obol view controller should render the full window viewport");
    struct bv_background_state background = BV_BACKGROUND_STATE_INIT;
    if (!bv_background_state_get(&background,
	    bv_context_view_const(static_cast<const struct bv_context *>(view.viewContext()))) ||
	background.bottom[0] != 110 || background.bottom[1] != 110 ||
	background.bottom[2] != 110 || background.top[0] != 0 ||
	background.top[1] != 0 || background.top[2] != 50)
	FAIL("qtcad local views should preserve the historical qged gradient");
    view.need_update(QG_VIEW_DRAWN);
    if (controller->getBackgroundBottomColor() !=
	    SbColor(110.0f / 255.0f, 110.0f / 255.0f, 110.0f / 255.0f) ||
	controller->getBackgroundTopColor() !=
	    SbColor(0.0f, 0.0f, 50.0f / 255.0f))
	FAIL("qtcad should synchronize its default gradient to Obol");

    struct bobol_endpoint_property_value background_property =
	BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
    background_property.type = BOBOL_ENDPOINT_PROPERTY_COLOR3;
    VSET(background_property.color3, 16.0 / 255.0, 32.0 / 255.0,
	48.0 / 255.0);
    if (bobol_display_endpoint_property_set(view.displayEndpoint(),
	"controller.background.bottom", &background_property) !=
	BOBOL_ENDPOINT_PROPERTY_OK)
	FAIL("qtcad background bottom should use the endpoint property");
    VSET(background_property.color3, 64.0 / 255.0, 80.0 / 255.0,
	96.0 / 255.0);
    if (bobol_display_endpoint_property_set(view.displayEndpoint(),
	"controller.background.top", &background_property) !=
	BOBOL_ENDPOINT_PROPERTY_OK)
	FAIL("qtcad background top should use the endpoint property");
    view.need_update(QG_VIEW_DRAWN);
    if (controller->getBackgroundBottomColor() !=
	    SbColor(16.0f / 255.0f, 32.0f / 255.0f, 48.0f / 255.0f) ||
	controller->getBackgroundTopColor() !=
	    SbColor(64.0f / 255.0f, 80.0f / 255.0f, 96.0f / 255.0f))
	FAIL("qtcad should preserve endpoint-owned gradient colors");

    /* A stale passive view value must not overwrite an endpoint policy during
     * a normal canvas refresh. */
    background = BV_BACKGROUND_STATE_INIT;
    VSET(background.bottom, 1, 2, 3);
    VSET(background.top, 4, 5, 6);
    (void)bv_background_state_set(bv_context_view(
	static_cast<struct bv_context *>(view.viewContext())), &background);
    view.need_update(QG_VIEW_DRAWN);
    if (controller->getBackgroundBottomColor() !=
	    SbColor(16.0f / 255.0f, 32.0f / 255.0f, 48.0f / 255.0f) ||
	controller->getBackgroundTopColor() !=
	    SbColor(64.0f / 255.0f, 80.0f / 255.0f, 96.0f / 255.0f))
	FAIL("qtcad refresh should not overwrite endpoint background policy");

    /* Renderer policy commands explicitly dirty the view without changing
     * its camera hash.  qged relies on diff_hashes() after each command to
     * schedule the actual Qt paint. */
    struct bv *localView = bv_context_view(
	static_cast<struct bv_context *>(view.viewContext()));
    (void)bv_refresh_complete(localView);
    view.stash_hashes();
    (void)bv_refresh_request(localView, BV_REFRESH_VIEW);
    if (!view.diff_hashes())
	FAIL("qtcad view diffs should honor explicit non-camera refresh requests");

    controller->clearRenderRequest();
    /* isRenderRequested() intentionally also reports independently queued
     * progressive work.  Test the explicit render-request latch directly so
     * an asynchronous LoD completion cannot make this assertion flaky. */
    if (controller->consumeRenderRequest(NULL))
	FAIL("Obol view controller should clear the explicit render request");
    controller->requestRender("qtcad-test");
    if (!controller->isRenderRequested() ||
	    bu_strcmp(controller->getRenderReason().getString(), "qtcad-test") != 0)
	FAIL("Obol view controller should retain render requests");

    SbVec3f beforeAET = camera->position.getValue();
    controller->clearRenderRequest();
    view.aet(90.0, 0.0, 0.0);
    SbVec3f afterAET = camera->position.getValue();
    if (fabsf(beforeAET[0] - afterAET[0]) < 0.001f &&
	    fabsf(beforeAET[1] - afterAET[1]) < 0.001f &&
	    fabsf(beforeAET[2] - afterAET[2]) < 0.001f)
	FAIL("QgView AET updates should synchronize the Obol camera");
    if (!controller->isRenderRequested() ||
	    bu_strcmp(controller->getRenderReason().getString(), "rt-view-camera") != 0)
	FAIL("Obol camera synchronization should request a render");
    SbString renderReason;
    if (!controller->consumeRenderRequest(&renderReason) ||
	    bu_strcmp(renderReason.getString(), "rt-view-camera") != 0)
	FAIL("Obol render requests should be consumable by render managers");
    if (controller->consumeRenderRequest(NULL))
	FAIL("Consumed Obol render requests should stay clear");
    if (controller->renderPending(FALSE, FALSE, NULL))
	FAIL("Obol render drain should be idle when no render is pending");

    const SbVec2s stableSize =
	controller->getViewportRegion().getWindowSize();
    controller->clearRenderRequest();
    controller->setViewportSize(static_cast<unsigned int>(stableSize[0]),
	static_cast<unsigned int>(stableSize[1]));
    if (controller->isRenderRequested())
	FAIL("unchanged viewport size should not request another render");

    void *viewCtx = view.viewContext();
    (void)controller->syncCameraFromViewContext(viewCtx);
    controller->clearRenderRequest();
    (void)controller->syncCameraFromViewContext(viewCtx);
    if (controller->isRenderRequested())
	FAIL("unchanged view camera should not request another render");

    if (!controller->getSceneRoot()->isOfType(SoSeparator::getClassTypeId()))
	FAIL("QgView default Obol scene root should be a separator");
    SoSeparator *sceneRoot = static_cast<SoSeparator *>(controller->getSceneRoot());
    sceneRoot->addChild(new SoDirectionalLight);
    SoMaterial *material = new SoMaterial;
    material->emissiveColor.setValue(1.0f, 0.2f, 0.1f);
    material->diffuseColor.setValue(1.0f, 0.2f, 0.1f);
    sceneRoot->addChild(material);
    SoCube *cube = new SoCube;
    cube->width = 250.0f;
    cube->height = 250.0f;
    cube->depth = 250.0f;
    sceneRoot->addChild(cube);
    controller->getViewport()->viewAll();
    if (camera->nearDistance.getValue() <= 0.0f ||
	    camera->farDistance.getValue() <= camera->nearDistance.getValue())
	FAIL("Obol orthographic viewAll should produce a valid clipping range");

    QImage obolImage;
    view.get_obol_viewport_image(obolImage);
    if (obolImage.isNull())
	FAIL("QgView should capture the Obol scene through SoOffscreenRenderer");
    int litPixels = lit_pixel_count(obolImage);
    if (obolImage.width() <= 0 || obolImage.height() <= 0 ||
	    litPixels < 10) {
	SbVec3f cpos = camera->position.getValue();
	fprintf(stderr, "Obol capture diagnostics: width=%d height=%d lit=%d\n",
		obolImage.width(), obolImage.height(), litPixels);
	fprintf(stderr, "Obol camera diagnostics: pos=(%g,%g,%g) focal=%g near=%g far=%g\n",
		cpos[0], cpos[1], cpos[2],
		camera->focalDistance.getValue(),
		camera->nearDistance.getValue(),
		camera->farDistance.getValue());
	FAIL("QgView Obol capture should produce non-empty pixels");
    }
    if (obolImage.pixelColor(2, 2) ==
	obolImage.pixelColor(2, obolImage.height() - 3))
	FAIL("QgView Obol capture should preserve configured gradient endpoints");

    controller->requestRender("sw-visible-readback");
    QImage visibleImage;
    view.get_viewport_image(visibleImage);
    if (visibleImage.isNull() || lit_pixel_count(visibleImage) < 10)
	FAIL("QgView visible SW capture should use populated Obol scenes");
    if (controller->isRenderRequested())
	FAIL("QgView visible SW capture should consume Obol render requests");

    struct bv *fpsView = bv_context_view(static_cast<struct bv_context *>(view.viewContext()));
    struct bv_params_state fpsParams = BV_PARAMS_STATE_INIT;
    (void)bv_params_state_get(&fpsParams, fpsView);
    fpsParams.draw = 1;
    fpsParams.draw_fps = 1;
    (void)bv_params_state_set(fpsView, &fpsParams);
    controller->requestRender("sw-visible-paint");
    view.show();
    QImage paintTarget(view.size(), QImage::Format_RGBA8888);
    paintTarget.fill(0);
    QPainter painter(&paintTarget);
    view.render(&painter);
    painter.end();
    if (controller->isRenderRequested())
	FAIL("QgView visible SW paint should consume Obol render requests");
    if (lit_pixel_count(paintTarget) < 10)
	FAIL("QgView visible SW paint should draw populated Obol scenes");
    if (controller->getLastRenderTimeNanoseconds() == 0 ||
	controller->getSmoothedRenderTimeNanoseconds() == 0)
	FAIL("Obol controller should record rendered frame telemetry");

    TestQgSW swCanvas(NULL);
    swCanvas.resize(160, 120);
    BObolViewController *swController = swCanvas.obolViewController();
    if (!swController ||
	    !swController->getSceneRoot()->isOfType(SoSeparator::getClassTypeId()))
	FAIL("QgSW input test should expose an Obol scene before first paint");
    SoSeparator *swRoot = static_cast<SoSeparator *>(swController->getSceneRoot());
    swRoot->addChild(new SoDirectionalLight);
    SoMaterial *swMaterial = new SoMaterial;
    swMaterial->emissiveColor.setValue(0.2f, 0.7f, 1.0f);
    swMaterial->diffuseColor.setValue(0.2f, 0.7f, 1.0f);
    swRoot->addChild(swMaterial);
    SoCube *swCube = new SoCube;
    swCube->width = 180.0f;
    swCube->height = 180.0f;
    swCube->depth = 180.0f;
    swRoot->addChild(swCube);
    swController->getViewport()->viewAll();
    SoCamera *swCamera = swController->getCamera();

    /* Exercise the canvas through an endpoint rather than its local fallback
     * context.  These bindings own view-local faceplate state, so they must
     * not quietly fall back to an application-global input path. */
    bobol_display_endpoint_t *swEndpoint =
	bobol_display_endpoint_create(swController, 0);
    if (!swEndpoint)
	FAIL("QgSW input test should create an endpoint");
    swCanvas.setObolInputEndpoint(swEndpoint);
    struct bv *swView = bv_context_view(swCanvas.viewContext());
    struct bv_adc_state swAdc = BV_ADC_STATE_INIT;
    struct bv_axes_state swModelAxes = BV_AXES_STATE_INIT;
    struct bv_axes_state swViewAxes = BV_AXES_STATE_INIT;
    (void)bv_adc_state_get(&swAdc, swView);
    (void)bv_model_axes_state_get(&swModelAxes, swView);
    (void)bv_view_axes_state_get(&swViewAxes, swView);
    const int beforeAdc = swAdc.draw;
    const int beforeModelAxes = swModelAxes.draw;
    const int beforeViewAxes = swViewAxes.draw;
    QKeyEvent adcKey(QEvent::KeyPress, Qt::Key_A, Qt::NoModifier, "A");
    QKeyEvent modelAxesKey(QEvent::KeyPress, Qt::Key_M, Qt::NoModifier, "M");
    QKeyEvent viewAxesKey(QEvent::KeyPress, Qt::Key_V, Qt::NoModifier, "V");
    swCanvas.runKeyPressForTest(&adcKey);
    swCanvas.runKeyPressForTest(&modelAxesKey);
    swCanvas.runKeyPressForTest(&viewAxesKey);
    (void)bv_adc_state_get(&swAdc, swView);
    (void)bv_model_axes_state_get(&swModelAxes, swView);
    (void)bv_view_axes_state_get(&swViewAxes, swView);
    if (swAdc.draw == beforeAdc ||
	swModelAxes.draw == beforeModelAxes ||
	swViewAxes.draw == beforeViewAxes)
	FAIL("QgSW endpoint input should toggle local faceplate state");

    SbVec3f beforeKey = swCamera->position.getValue();
    swController->clearRenderRequest();
    QKeyEvent keyEvent(QEvent::KeyPress, Qt::Key_T, Qt::NoModifier, "T");
    swCanvas.runKeyPressForTest(&keyEvent);
    if (!camera_moved(swCamera, beforeKey))
	FAIL("QgSW default key navigation should update the Obol camera");
    if (!swController->isRenderRequested())
	FAIL("QgSW default key navigation should request an Obol render");

    SbVec3f beforeWheel = swCamera->position.getValue();
    float beforeWheelFocal = swCamera->focalDistance.getValue();
    swController->clearRenderRequest();
    QWheelEvent swWheel = wheel_event(80, 60, 120);
    swCanvas.runWheelForTest(&swWheel);
    if (!camera_state_changed(swCamera, beforeWheel, beforeWheelFocal))
	FAIL("QgSW wheel navigation should update the Obol camera");
    if (!swController->isRenderRequested())
	FAIL("QgSW wheel navigation should request an Obol render");

    SbVec3f beforeDrag = swCamera->position.getValue();
    float beforeDragFocal = swCamera->focalDistance.getValue();
    swController->clearRenderRequest();
    QMouseEvent swMoveStart = mouse_move_event(80, 60, Qt::LeftButton);
    QMouseEvent swMoveDrag = mouse_move_event(115, 90, Qt::LeftButton);
    swCanvas.runMouseMoveForTest(&swMoveStart);
    swCanvas.runMouseMoveForTest(&swMoveDrag);
    if (!camera_state_changed(swCamera, beforeDrag, beforeDragFocal))
	FAIL("QgSW drag navigation should update the Obol camera");
    if (!swController->isRenderRequested())
	FAIL("QgSW drag navigation should request an Obol render");

    swCanvas.setObolInputEndpoint(NULL);
    bobol_display_endpoint_destroy(swEndpoint);

    QgView glView(NULL, QgViewType::GL);
    glView.resize(128, 96);
    if (glView.view_type() == QgViewType::GL) {
	QgGL *glPresentationCanvas = glView.canvasBase() ?
	    dynamic_cast<QgGL *>(glView.canvasBase()) : NULL;
	if (!glPresentationCanvas || glPresentationCanvas->updateBehavior() !=
		QOpenGLWidget::PartialUpdate)
	    FAIL("QgGL should preserve the last completed Obol frame while asynchronous replacement work is pending");
	BObolViewController *glController = glView.obolViewController();
	if (!glController ||
		!glController->getSceneRoot()->isOfType(SoSeparator::getClassTypeId()))
	    FAIL("QgGL should expose an Obol scene before first paint");
	SoGLRenderAction *glAction =
	    glController->getRenderManager()->getGLRenderAction();
	SoGLRenderAction *swAction =
	    controller->getRenderManager()->getGLRenderAction();
	if (!glAction || !glAction->getContextManager() || !swAction ||
		!swAction->getContextManager() ||
		glAction->getContextManager() == swAction->getContextManager())
	    FAIL("QgGL and QgSW should retain separate Obol context managers");

	SoSeparator *glRoot = static_cast<SoSeparator *>(glController->getSceneRoot());
	glRoot->addChild(new SoDirectionalLight);
	SoMaterial *glMaterial = new SoMaterial;
	glMaterial->emissiveColor.setValue(0.1f, 0.9f, 0.3f);
	glMaterial->diffuseColor.setValue(0.1f, 0.9f, 0.3f);
	glRoot->addChild(glMaterial);
	SoCube *glCube = new SoCube;
	glCube->width = 180.0f;
	glCube->height = 180.0f;
	glCube->depth = 180.0f;
	glRoot->addChild(glCube);
	glController->getViewport()->viewAll();

	glController->requestRender("gl-visible-readback");
	QImage glVisibleImage;
	glView.get_viewport_image(glVisibleImage);
	if (glVisibleImage.isNull() || lit_pixel_count(glVisibleImage) < 10)
	    FAIL("QgGL visible capture should use Obol readback before first paint");
	if (glController->isRenderRequested())
	    FAIL("QgGL visible capture should consume Obol render requests");

#ifdef BRLCAD_OPENGL
	TestQgGL glCanvas(NULL);
	glCanvas.resize(128, 96);
	glCanvas.show();
	app.processEvents();
	if (glCanvas.isValid()) {
	    BObolViewController *paintController = glCanvas.obolViewController();
	    if (!paintController ||
		    !paintController->getSceneRoot()->isOfType(SoSeparator::getClassTypeId()))
		FAIL("QgGL paint test should expose an Obol scene");
	    SoGLRenderAction *paintAction =
		paintController->getRenderManager()->getGLRenderAction();
	    if (!paintAction || !paintAction->getContextManager() ||
		paintAction->getContextManager() == swAction->getContextManager())
		FAIL("QgGL direct rendering should retain its Qt context manager");
	    SoSeparator *paintRoot = static_cast<SoSeparator *>(paintController->getSceneRoot());
	    paintRoot->addChild(new SoDirectionalLight);
	    SoMaterial *paintMaterial = new SoMaterial;
	    paintMaterial->emissiveColor.setValue(0.9f, 0.8f, 0.1f);
	    paintMaterial->diffuseColor.setValue(0.9f, 0.8f, 0.1f);
	    paintRoot->addChild(paintMaterial);
	    SoCube *paintCube = new SoCube;
	    paintCube->width = 180.0f;
	    paintCube->height = 180.0f;
	    paintCube->depth = 180.0f;
	    paintRoot->addChild(paintCube);
	    paintController->getViewport()->viewAll();
	    SoCamera *paintCamera = paintController->getCamera();
	    SbVec3f beforeGLKey = paintCamera->position.getValue();
	    paintController->clearRenderRequest();
	    QKeyEvent glKeyEvent(QEvent::KeyPress, Qt::Key_T, Qt::NoModifier, "T");
	    glCanvas.runKeyPressForTest(&glKeyEvent);
	    if (!camera_moved(paintCamera, beforeGLKey))
		FAIL("QgGL default key navigation should update the Obol camera");
	    if (!paintController->isRenderRequested())
		FAIL("QgGL default key navigation should request an Obol render");

	    SbVec3f beforeGLWheel = paintCamera->position.getValue();
	    float beforeGLWheelFocal = paintCamera->focalDistance.getValue();
	    paintController->clearRenderRequest();
	    QWheelEvent glWheel = wheel_event(64, 48, 120);
	    glCanvas.runWheelForTest(&glWheel);
	    if (!camera_state_changed(paintCamera, beforeGLWheel, beforeGLWheelFocal))
		FAIL("QgGL wheel navigation should update the Obol camera");
	    if (!paintController->isRenderRequested())
		FAIL("QgGL wheel navigation should request an Obol render");

	    SbVec3f beforeGLDrag = paintCamera->position.getValue();
	    float beforeGLDragFocal = paintCamera->focalDistance.getValue();
	    paintController->clearRenderRequest();
	    QMouseEvent glMoveStart = mouse_move_event(64, 48, Qt::LeftButton);
	    QMouseEvent glMoveDrag = mouse_move_event(100, 80, Qt::LeftButton);
	    glCanvas.runMouseMoveForTest(&glMoveStart);
	    glCanvas.runMouseMoveForTest(&glMoveDrag);
	    if (!camera_state_changed(paintCamera, beforeGLDrag, beforeGLDragFocal))
		FAIL("QgGL drag navigation should update the Obol camera");
	    if (!paintController->isRenderRequested())
		FAIL("QgGL drag navigation should request an Obol render");
	    paintController->requestRender("gl-visible-paint");

	    glCanvas.makeCurrent();
	    glCanvas.runPaintGLForTest();
	    glCanvas.doneCurrent();
	    if (paintController->isRenderRequested())
		FAIL("QgGL visible paint should consume Obol render requests");
	}
#endif
    }

    return 0;
}
