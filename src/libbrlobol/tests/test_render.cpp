/*                     T E S T _ R E N D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol.h"

#include <Inventor/SbColor.h>
#include <Inventor/SbString.h>
#include <Inventor/SbViewVolume.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/SoPrimitiveVertex.h>
#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/SoSceneCollector.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/annex/HUD/nodes/SoHUDButton.h>
#include <Inventor/annex/HUD/nodes/SoHUDLabel.h>
#include <Inventor/gl.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoDrawStyle.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoShape.h>
#include <Inventor/nodes/SoSeparator.h>
#include <obol/cad/SoCADAssembly.h>
#include <OSMesa/osmesa.h>

#include <algorithm>
#include <math.h>
#include <memory>
#include <stdio.h>
#include <string.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

class BRLOBOLRenderContext
{
public:
    BRLOBOLRenderContext(unsigned int w, unsigned int h) :
	context(NULL),
	buffer(NULL),
	width(w),
	height(h),
	previousContext(NULL),
	previousBuffer(NULL),
	previousWidth(0),
	previousHeight(0),
	previousBytesPerRow(0),
	previousFormat(0)
    {
	context = OSMesaCreateContextExt(OSMESA_RGBA, 24, 0, 0, NULL);
	if (context)
	    buffer.reset(new unsigned char[width * height * 4]);
    }

    ~BRLOBOLRenderContext(void)
    {
	if (context)
	    OSMesaDestroyContext(context);
    }

    bool isValid(void) const
    {
	return context != NULL && buffer.get() != NULL;
    }

    SbBool makeCurrent(void)
    {
	if (!this->isValid())
	    return FALSE;

	previousContext = OSMesaGetCurrentContext();
	previousBuffer = NULL;
	previousWidth = 0;
	previousHeight = 0;
	previousBytesPerRow = 0;
	previousFormat = 0;
	if (previousContext) {
	    GLint fmt = 0;
	    OSMesaGetColorBuffer(previousContext, &previousWidth, &previousHeight,
				 &fmt, &previousBuffer);
	    previousFormat = (GLenum)fmt;
	}

	return OSMesaMakeCurrent(context, buffer.get(), GL_UNSIGNED_BYTE,
				 static_cast<GLsizei>(width), static_cast<GLsizei>(height)) ? TRUE : FALSE;
    }

    void restorePrevious(void)
    {
	if (previousContext && previousBuffer) {
	    OSMesaMakeCurrent(previousContext, previousBuffer, GL_UNSIGNED_BYTE,
			      previousWidth, previousHeight);
	} else {
	    OSMesaMakeCurrent(NULL, NULL, 0, 0, 0);
	}
    }

    OSMesaContext context;
    std::unique_ptr<unsigned char[]> buffer;
    unsigned int width;
    unsigned int height;
    OSMesaContext previousContext;
    void *previousBuffer;
    GLsizei previousWidth;
    GLsizei previousHeight;
    GLsizei previousBytesPerRow;
    GLenum previousFormat;
};

class BRLOBOLRenderContextManager : public SoDB::ContextManager
{
public:
    virtual void *createOffscreenContext(unsigned int width, unsigned int height)
    {
	BRLOBOLRenderContext *ctx = new BRLOBOLRenderContext(width, height);
	if (ctx->isValid())
	    return ctx;
	delete ctx;
	return NULL;
    }

    virtual SbBool isOSMesaContext(void *UNUSED(context))
    {
	return TRUE;
    }

    virtual SbBool makeContextCurrent(void *context)
    {
	return context ? static_cast<BRLOBOLRenderContext *>(context)->makeCurrent() : FALSE;
    }

    virtual void restorePreviousContext(void *context)
    {
	if (context)
	    static_cast<BRLOBOLRenderContext *>(context)->restorePrevious();
    }

    virtual void destroyContext(void *context)
    {
	BRLOBOLRenderContext *ctx = static_cast<BRLOBOLRenderContext *>(context);
	if (ctx && ctx->context && OSMesaGetCurrentContext() == ctx->context)
	    OSMesaMakeCurrent(NULL, NULL, 0, 0, 0);
	delete ctx;
    }

    virtual void *getProcAddress(const char *funcName)
    {
	return reinterpret_cast<void *>(OSMesaGetProcAddress(funcName));
    }
};

struct BRLOBOLSoftwareLineState {
    unsigned char *pixels;
    unsigned int width;
    unsigned int height;
    unsigned int components;
    SbViewVolume viewVolume;
    int lineCount;
};

static void
set_pixel(BRLOBOLSoftwareLineState *state, int x, int y)
{
    if (!state || !state->pixels || x < 0 || y < 0 ||
	x >= static_cast<int>(state->width) ||
	y >= static_cast<int>(state->height))
	return;

    unsigned char *p = state->pixels + ((y * state->width + x) * state->components);
    if (state->components == 1) {
	p[0] = 255;
    } else if (state->components == 2) {
	p[0] = 255;
	p[1] = 255;
    } else {
	p[0] = 255;
	p[1] = 64;
	p[2] = 32;
	if (state->components > 3)
	    p[3] = 255;
    }
}

static void
draw_line(BRLOBOLSoftwareLineState *state, const SbVec3f &a, const SbVec3f &b)
{
    SbVec3f sa;
    SbVec3f sb;
    state->viewVolume.projectToScreen(a, sa);
    state->viewVolume.projectToScreen(b, sb);

    int x0 = static_cast<int>(floorf(sa[0] * static_cast<float>(state->width - 1) + 0.5f));
    int y0 = static_cast<int>(floorf(sa[1] * static_cast<float>(state->height - 1) + 0.5f));
    int x1 = static_cast<int>(floorf(sb[0] * static_cast<float>(state->width - 1) + 0.5f));
    int y1 = static_cast<int>(floorf(sb[1] * static_cast<float>(state->height - 1) + 0.5f));

    int dx = abs(x1 - x0);
    int sx = x0 < x1 ? 1 : -1;
    int dy = -abs(y1 - y0);
    int sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;

    for (;;) {
	set_pixel(state, x0, y0);
	set_pixel(state, x0 + 1, y0);
	set_pixel(state, x0 - 1, y0);
	set_pixel(state, x0, y0 + 1);
	set_pixel(state, x0, y0 - 1);
	if (x0 == x1 && y0 == y1)
	    break;
	int e2 = 2 * err;
	if (e2 >= dy) {
	    err += dy;
	    x0 += sx;
	}
	if (e2 <= dx) {
	    err += dx;
	    y0 += sy;
	}
    }
}

static void
software_line_cb(void *userdata, SoCallbackAction *action,
		 const SoPrimitiveVertex *v1, const SoPrimitiveVertex *v2)
{
    BRLOBOLSoftwareLineState *state = static_cast<BRLOBOLSoftwareLineState *>(userdata);
    if (!state || !v1 || !v2)
	return;

    SbVec3f a;
    SbVec3f b;
    const SbMatrix &modelMatrix = action->getModelMatrix();
    modelMatrix.multVecMatrix(v1->getPoint(), a);
    modelMatrix.multVecMatrix(v2->getPoint(), b);
    draw_line(state, a, b);
    state->lineCount++;
}

class BRLOBOLSoftwareLineContextManager : public SoDB::ContextManager
{
public:
    BRLOBOLSoftwareLineContextManager(void) :
	lastLineCount(0),
	lastOverlayCount(0)
    {
    }

    virtual void *createOffscreenContext(unsigned int UNUSED(width), unsigned int UNUSED(height))
    {
	return NULL;
    }

    virtual SbBool makeContextCurrent(void *UNUSED(context))
    {
	return FALSE;
    }

    virtual void restorePreviousContext(void *UNUSED(context))
    {
    }

    virtual void destroyContext(void *UNUSED(context))
    {
    }

    virtual SbBool renderScene(SoNode *scene, unsigned int width, unsigned int height,
			       unsigned char *pixels, unsigned int nrcomponents,
			       const float background_rgb[3])
    {
	(void)background_rgb;
	if (!scene || !pixels || width == 0 || height == 0 ||
	    nrcomponents < 1 || nrcomponents > 4)
	    return FALSE;

	SoSearchAction search;
	search.setType(SoCamera::getClassTypeId());
	search.setInterest(SoSearchAction::FIRST);
	search.apply(scene);
	SoPath *cameraPath = search.getPath();
	if (!cameraPath || cameraPath->getLength() < 1)
	    return FALSE;

	SoNode *tail = cameraPath->getTail();
	if (!tail || !tail->isOfType(SoCamera::getClassTypeId()))
	    return FALSE;
	SoCamera *camera = static_cast<SoCamera *>(tail);

	BRLOBOLSoftwareLineState state;
	state.pixels = pixels;
	state.width = width;
	state.height = height;
	state.components = nrcomponents;
	state.viewVolume = camera->getViewVolume(static_cast<float>(width) / static_cast<float>(height));
	state.lineCount = 0;

	SoCallbackAction callback;
	callback.addLineSegmentCallback(SoShape::getClassTypeId(), software_line_cb, &state);
	callback.apply(scene);

	SoSceneCollector collector;
	collector.collectOverlaysOnly(scene, SbViewportRegion(width, height));
	this->lastLineCount = state.lineCount;
	this->lastOverlayCount = static_cast<int>(collector.getOverlays().size());
	collector.compositeOverlays(pixels, width, height, nrcomponents);

	return (state.lineCount > 0 || !collector.getOverlays().empty()) ? TRUE : FALSE;
    }

    int lastLineCount;
    int lastOverlayCount;
};

class BRLOBOLDiagnosticContextManager : public SoDB::ContextManager
{
public:
    BRLOBOLDiagnosticContextManager(void) :
	diagnosticCount(0),
	lastComponent(""),
	lastMessage("")
    {
    }

    virtual void *createOffscreenContext(unsigned int UNUSED(width), unsigned int UNUSED(height))
    {
	return NULL;
    }

    virtual SbBool makeContextCurrent(void *UNUSED(context))
    {
	return FALSE;
    }

    virtual void restorePreviousContext(void *UNUSED(context))
    {
    }

    virtual void destroyContext(void *UNUSED(context))
    {
    }

    virtual SbBool renderScene(SoNode *UNUSED(scene), unsigned int UNUSED(width),
			       unsigned int UNUSED(height), unsigned char *UNUSED(pixels),
			       unsigned int UNUSED(nrcomponents), const float background_rgb[3])
    {
	(void)background_rgb;
	return FALSE;
    }

    virtual void reportDiagnostic(const char *component, const char *message)
    {
	this->diagnosticCount++;
	this->lastComponent = component ? component : "";
	this->lastMessage = message ? message : "";
    }

    int diagnosticCount;
    SbString lastComponent;
    SbString lastMessage;
};

static int
count_lit_pixels(const unsigned char *buffer, int width, int height)
{
    int count = 0;
    for (int i = 0; i < width * height; i++) {
	const unsigned char *p = buffer + i * 3;
	if (p[0] > 32 || p[1] > 32 || p[2] > 32)
	    count++;
    }
    return count;
}

static int
count_green_pixels(const unsigned char *buffer, int width, int height)
{
    int count = 0;
    for (int i = 0; i < width * height; i++) {
	const unsigned char *p = buffer + i * 3;
	if (p[1] > 96 && p[0] < 80 && p[2] < 80)
	    count++;
    }
    return count;
}

static int
count_blue_pixels(const unsigned char *buffer, int width, int height)
{
    int count = 0;
    for (int i = 0; i < width * height; i++) {
	const unsigned char *p = buffer + i * 3;
	if (p[2] > 96 && p[0] < 80 && p[1] < 80)
	    count++;
    }
    return count;
}

static int
count_yellow_pixels(const unsigned char *buffer, int width, int height)
{
    int count = 0;
    for (int i = 0; i < width * height; i++) {
	const unsigned char *p = buffer + i * 3;
	if (p[0] > 96 && p[1] > 96 && p[2] < 80)
	    count++;
    }
    return count;
}

int
main(int UNUSED(argc), const char **UNUSED(argv))
{
    const int width = 160;
    const int height = 120;

    BRLOBOLRenderContextManager contextManager;
    brlobol_init(&contextManager);
    if (SoDB::getContextManager() != &contextManager)
	FAIL("brlobol_init should install the application-provided context manager");

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->position = SbVec3f(0.0f, 0.0f, 10.0f);
    camera->height = 6.0f;
    camera->nearDistance = 1.0f;
    camera->farDistance = 20.0f;
    root->addChild(camera);

    SoBaseColor *color = new SoBaseColor;
    color->rgb = SbColor(1.0f, 0.2f, 0.1f);
    root->addChild(color);

    SoDrawStyle *style = new SoDrawStyle;
    style->lineWidth = 4.0f;
    root->addChild(style);

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->path = "/prototype/render-square";
    source->sourceRevision = 1;
    root->addChild(source);

    SoBRLHUDLabelOverlay *brlLabel = new SoBRLHUDLabelOverlay;
    brlLabel->labelId = "render::hud-label";
    brlLabel->sourceId = 7;
    brlLabel->text = "BRL-CAD Obol label";
    brlLabel->position = SbVec2f(8.0f, 86.0f);
    brlLabel->color = SbColor(0.0f, 1.0f, 0.0f);
    brlLabel->fontSize = 12.0f;
    if (!brlLabel->rebuildGeometry() || !brlLabel->getHUDLabel())
	FAIL("BRL-CAD HUD label overlay should produce direct Obol HUD label nodes");
    root->addChild(brlLabel);

    SoHUDKit *hud = new SoHUDKit;
    root->addChild(hud);

    SoHUDLabel *label = new SoHUDLabel;
    label->position = SbVec2f(8.0f, 104.0f);
    label->string.set1Value(0, "BRL-CAD Obol HUD");
    label->color = SbColor(0.0f, 1.0f, 0.0f);
    label->fontSize = 14.0f;
    hud->addWidget(label);

    SoHUDLabel *status = new SoHUDLabel;
    status->position = SbVec2f(8.0f, 12.0f);
    status->string.set1Value(0, "FB: ready  scale 1:1");
    status->color = SbColor(0.0f, 0.75f, 1.0f);
    status->fontSize = 12.0f;
    hud->addWidget(status);

    SoHUDButton *rubberBand = new SoHUDButton;
    rubberBand->position = SbVec2f(86.0f, 18.0f);
    rubberBand->size = SbVec2f(62.0f, 28.0f);
    rubberBand->string = "select";
    rubberBand->color = SbColor(1.0f, 1.0f, 0.0f);
    rubberBand->borderColor = SbColor(1.0f, 1.0f, 0.0f);
    rubberBand->fontSize = 10.0f;
    hud->addWidget(rubberBand);

    SoBRLRealizeAction realize;
    realize.apply(root);
    if (source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED)
	FAIL("render source should realize before rendering");
    SoBRLVListShape *shape = source->getRealizedShape();
    if (!shape)
	FAIL("render source should expose realized line geometry");
    shape->colorOverride = TRUE;
    shape->color = SbColor(0.0f, 0.0f, 1.0f);

    SbViewportRegion viewport(width, height);
    SoOffscreenRenderer renderer(&contextManager, viewport);
    renderer.setComponents(SoOffscreenRenderer::RGB);
    renderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (!renderer.render(root))
	FAIL("SoOffscreenRenderer should render BRL-CAD Obol geometry through OSMesa");

    const unsigned char *buffer = renderer.getBuffer();
    if (!buffer)
	FAIL("SoOffscreenRenderer should expose a render buffer");

    int litPixels = count_lit_pixels(buffer, width, height);
    if (litPixels < 100)
	FAIL("rendered BRL-CAD Obol geometry should produce visible non-background pixels");
    int bluePixels = count_blue_pixels(buffer, width, height);
    if (bluePixels < 20)
	FAIL("OSMesa GL render should honor BRL-CAD shape color override");
    int greenPixels = count_green_pixels(buffer, width, height);
    if (greenPixels < 20)
	FAIL("OSMesa GL render should render Obol HUD labels");

    {
	SoOffscreenRenderer secondRenderer(&contextManager, viewport);
	secondRenderer.setComponents(SoOffscreenRenderer::RGB);
	secondRenderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
	if (!secondRenderer.render(root) || !secondRenderer.getBuffer())
	    FAIL("retained CAD scene should render in a second GL context");
    }

    SoSeparator *pointRoot = new SoSeparator;
    pointRoot->ref();
    SoOrthographicCamera *pointCamera = new SoOrthographicCamera;
    pointCamera->position = SbVec3f(0.0f, 0.0f, 10.0f);
    pointCamera->height = 6.0f;
    pointCamera->nearDistance = 1.0f;
    pointCamera->farDistance = 20.0f;
    pointRoot->addChild(pointCamera);
    SoCADAssembly *pointAssembly = new SoCADAssembly;
    obol::PartGeometry pointGeometry;
    obol::PointRep points;
    points.positions = {SbVec3f(-1.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f)};
    points.pointIds = {7u, 11u};
    points.colorValid = {1u, 1u};
    points.colors = {SbColor(1.0f, 0.0f, 0.0f),
	SbColor(0.0f, 1.0f, 0.0f)};
    points.bounds.setBounds(SbVec3f(-1.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f));
    pointGeometry.points = points;
    const obol::PartId pointPart = obol::CadIdBuilder::hash128("render-points");
    pointAssembly->upsertPart(pointPart, pointGeometry);
    obol::InstanceRecord pointInstance;
    pointInstance.part = pointPart;
    pointInstance.parent = obol::CadIdBuilder::Root();
    pointInstance.childName = "render-points";
    pointInstance.localToRoot.makeIdentity();
    pointInstance.style.lineWidth = 9.0f;
    pointAssembly->upsertInstanceAuto(pointInstance);
    pointRoot->addChild(pointAssembly);

    SoOffscreenRenderer pointRenderer(&contextManager, viewport);
    pointRenderer.setComponents(SoOffscreenRenderer::RGB);
    pointRenderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (!pointRenderer.render(pointRoot) || !pointRenderer.getBuffer())
	FAIL("retained Obol point geometry should render through OSMesa");
    if (count_lit_pixels(pointRenderer.getBuffer(), width, height) < 20 ||
	count_green_pixels(pointRenderer.getBuffer(), width, height) < 5)
	FAIL("retained Obol point rendering should preserve point size and color");
    pointRoot->unref();

    BRLOBOLSoftwareLineContextManager softwareManager;
    SoOffscreenRenderer softwareRenderer(&softwareManager, viewport);
    softwareRenderer.setComponents(SoOffscreenRenderer::RGB);
    softwareRenderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (!softwareRenderer.render(root))
	FAIL("software ContextManager should render BRL-CAD Obol line geometry");

    buffer = softwareRenderer.getBuffer();
    if (!buffer)
	FAIL("software ContextManager should expose a render buffer");

    litPixels = count_lit_pixels(buffer, width, height);
    if (litPixels < 100)
	FAIL("software ContextManager should produce visible non-background pixels");

    greenPixels = count_green_pixels(buffer, width, height);
    if (greenPixels < 20)
	FAIL("software ContextManager should composite Obol HUD label overlays");
    int yellowPixels = count_yellow_pixels(buffer, width, height);
    if (yellowPixels < 20)
	FAIL("software ContextManager should composite direct Obol HUD rubber-band/button overlays");
    if (softwareManager.lastOverlayCount < 3)
	FAIL("software ContextManager should collect multiple direct Obol HUD overlays");

    BRLOBOLDiagnosticContextManager diagnosticManager;
    SoOffscreenRenderer diagnosticRenderer(&diagnosticManager, viewport);
    diagnosticRenderer.setComponents(SoOffscreenRenderer::RGB);
    diagnosticRenderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (diagnosticRenderer.render(root))
	FAIL("diagnostic ContextManager should decline renderScene and fail without GL context");
    if (diagnosticManager.diagnosticCount <= 0 ||
	strstr(diagnosticManager.lastComponent.getString(), "SoOffscreenRenderer") == NULL ||
	strstr(diagnosticManager.lastMessage.getString(), "renderScene") == NULL)
	FAIL("ContextManager diagnostic hook should report backend renderScene fallback");

    root->unref();
    return 0;
}
