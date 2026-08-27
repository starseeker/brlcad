/*                     T E S T _ R E N D E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include <OSMesa/osmesa.h>

#include "bu/app.h"

#include "BObol.h"
#include "BObol/BHUDLabelOverlay.h"
#include "BObol/BRealizeAction.h"
#include "BObol/BVListShape.h"

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
#include <Obol/cad/SoCADAssembly.h>
#include <Obol/cad/SoCADViewState.h>
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

static Obol::CadGeometryValidation
admitAndUpsertPart(SoCADAssembly *assembly, Obol::PartId part,
	Obol::PartGeometryBuilder geometry)
{
    const Obol::CadGeometryAdmission admission =
	Obol::cadAdmitPartGeometry(std::move(geometry));
    if (!admission)
	return admission.validation;
    return assembly->upsertParts({{part, admission.geometry}});
}

class BOBOLRenderContext
{
public:
    BOBOLRenderContext(unsigned int w, unsigned int h) :
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

    ~BOBOLRenderContext(void)
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

class BOBOLRenderContextManager : public SoDB::ContextManager
{
public:
    BOBOLRenderContextManager(void) :
	manager(SoDB::createOSMesaContextManager()),
	softwareFramebufferQueries(0)
    {
    }

    ~BOBOLRenderContextManager(void)
    {
	delete manager;
    }

    virtual void *createOffscreenContext(unsigned int width, unsigned int height)
    {
	return manager ? manager->createOffscreenContext(width, height) : NULL;
    }

    virtual SbBool isOSMesaContext(void *context)
    {
	return manager ? manager->isOSMesaContext(context) : FALSE;
    }

    virtual void maxOffscreenDimensions(unsigned int &width,
	unsigned int &height) const
    {
	if (manager)
	    manager->maxOffscreenDimensions(width, height);
    }

    virtual SbBool makeContextCurrent(void *context)
    {
	return manager ? manager->makeContextCurrent(context) : FALSE;
    }

    virtual void restorePreviousContext(void *context)
    {
	if (manager)
	    manager->restorePreviousContext(context);
    }

    virtual void destroyContext(void *context)
    {
	if (manager)
	    manager->destroyContext(context);
    }

    virtual SbBool getCurrentSoftwareFramebuffer(unsigned char *&pixels,
	unsigned int &width, unsigned int &height, unsigned int &components)
    {
	this->softwareFramebufferQueries++;
	return manager ? manager->getCurrentSoftwareFramebuffer(
	    pixels, width, height, components) : FALSE;
    }

    virtual void *getProcAddress(const char *funcName)
    {
	return manager ? manager->getProcAddress(funcName) : NULL;
    }

    SoDB::ContextManager *manager;
    int softwareFramebufferQueries;
};

struct BOBOLSoftwareLineState {
    unsigned char *pixels;
    unsigned int width;
    unsigned int height;
    unsigned int components;
    SbViewVolume viewVolume;
    int lineCount;
};

static void
set_pixel(BOBOLSoftwareLineState *state, int x, int y)
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
draw_line(BOBOLSoftwareLineState *state, const SbVec3f &a, const SbVec3f &b)
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
    BOBOLSoftwareLineState *state = static_cast<BOBOLSoftwareLineState *>(userdata);
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

class BOBOLSoftwareLineContextManager : public SoDB::ContextManager
{
public:
    BOBOLSoftwareLineContextManager(void) :
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

	BOBOLSoftwareLineState state;
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

class BOBOLDiagnosticContextManager : public SoDB::ContextManager
{
public:
    BOBOLDiagnosticContextManager(void) :
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

static int
count_red_pixels_in_box(const unsigned char *buffer, int width, int height,
	int xmin, int ymin, int xmax, int ymax)
{
    int count = 0;
    xmin = (std::max)(0, xmin);
    ymin = (std::max)(0, ymin);
    xmax = (std::min)(width - 1, xmax);
    ymax = (std::min)(height - 1, ymax);
    for (int y = ymin; y <= ymax; y++) {
	for (int x = xmin; x <= xmax; x++) {
	    const unsigned char *p = buffer + (y * width + x) * 3;
	    if (p[0] > 96 && p[1] < 80 && p[2] < 80)
		count++;
	}
    }
    return count;
}

static int
test_independent_line_segments(SoDB::ContextManager *contextManager)
{
    const int width = 160;
    const int height = 120;
    SoSeparator *root = new SoSeparator;
    root->ref();

    SoOrthographicCamera *camera = new SoOrthographicCamera;
    camera->position = SbVec3f(0.0f, 0.0f, 10.0f);
    camera->height = 6.0f;
    camera->nearDistance = 1.0f;
    camera->farDistance = 20.0f;
    root->addChild(camera);

    const SbVec3f points[] = {
	SbVec3f(-2.0f, -2.0f, 0.0f), SbVec3f(-2.0f,  2.0f, 0.0f),
	SbVec3f(-2.0f,  2.0f, 0.0f), SbVec3f( 2.0f,  2.0f, 0.0f),
	SbVec3f( 2.0f,  2.0f, 0.0f), SbVec3f( 2.0f, -2.0f, 0.0f),
	SbVec3f( 2.0f, -2.0f, 0.0f), SbVec3f(-2.0f, -2.0f, 0.0f)
    };
    const int32_t commands[] = {
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW,
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW,
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW,
	SoBRLVListShape::MOVE, SoBRLVListShape::DRAW
    };
    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->setLineSet(points, commands, 8);
    shape->sourceType = "view-polygon-edge";
    shape->colorOverride = TRUE;
    shape->color = SbColor(1.0f, 0.0f, 0.0f);
    shape->lineWidth = 1;
    root->addChild(shape);

    SbViewportRegion viewport(width, height);
    SoOffscreenRenderer renderer(contextManager, viewport);
    renderer.setComponents(SoOffscreenRenderer::RGB);
    renderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (!renderer.render(root) || !renderer.getBuffer()) {
	root->unref();
	FAIL("independent polygon edge pairs should render through OSMesa");
    }

    const unsigned char *buffer = renderer.getBuffer();
    /* With height 6 and a 4:3 viewport, x +/-2 maps to 40/120 and
	* y +/-2 maps to 20/100.  Test the middle of every edge separately
	* so a three-sided outline cannot satisfy a total-pixel assertion. */
    if (count_red_pixels_in_box(buffer, width, height, 37, 30, 43, 90) < 20 ||
	count_red_pixels_in_box(buffer, width, height, 50, 97, 110, 103) < 20 ||
	count_red_pixels_in_box(buffer, width, height, 117, 30, 123, 90) < 20 ||
	count_red_pixels_in_box(buffer, width, height, 50, 17, 110, 23) < 20) {
	root->unref();
	FAIL("independent polygon rendering should preserve every outline edge");
    }

    root->unref();
    return 0;
}

int
main(int UNUSED(argc), const char **UNUSED(argv))
{
    bu_setprogname("test_bobol_render");
    const int width = 160;
    const int height = 120;

    BOBOLRenderContextManager contextManager;
    bobol_init(&contextManager);
    if (SoDB::getContextManager() != &contextManager)
	FAIL("bobol_init should install the application-provided context manager");
    {
	/* Keep this renderer's context lifecycle independent of the policy
	 * counters exercised below. */
	BOBOLRenderContextManager segmentContextManager;
	if (test_independent_line_segments(&segmentContextManager) != 0)
	    return 1;
    }

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

    SoCADViewState *cadViewState = new SoCADViewState;
    cadViewState->softwareWireMode = SoCADViewState::SOFTWARE_WIRE_QUALITY;
    root->addChild(cadViewState);

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
	FAIL("QUALITY software-wire policy should render BRL-CAD Obol geometry through OSMesa GL");
    const int qualityFramebufferQueries =
	contextManager.softwareFramebufferQueries;

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

    cadViewState->softwareWireMode = SoCADViewState::SOFTWARE_WIRE_FAST;
    contextManager.softwareFramebufferQueries = 0;
    SoOffscreenRenderer fastRenderer(&contextManager, viewport);
    fastRenderer.setComponents(SoOffscreenRenderer::RGB);
    fastRenderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (!fastRenderer.render(root) || !fastRenderer.getBuffer())
	FAIL("FAST software-wire policy should render through the OSMesa framebuffer");
    if (contextManager.softwareFramebufferQueries <=
	qualityFramebufferQueries)
	FAIL("FAST software-wire policy should add a direct writable-framebuffer raster pass");
    if (count_blue_pixels(fastRenderer.getBuffer(), width, height) < 20)
	FAIL("FAST software-wire policy should preserve visible database wire color");
    if (count_green_pixels(fastRenderer.getBuffer(), width, height) < 20)
	FAIL("FAST software-wire policy should preserve retained HUD composition");
    cadViewState->softwareWireMode = SoCADViewState::SOFTWARE_WIRE_QUALITY;

    /* A direct software-framebuffer executor is still part of the retained
     * CAD renderer contract.  In particular it must publish exact work and
     * advance the completed resource-frame token; otherwise the view-LoD
     * controller observes a successful image as zero work and resubmits it
     * forever. */
    SoSeparator *directRoot = new SoSeparator;
    directRoot->ref();
    SoOrthographicCamera *directCamera = new SoOrthographicCamera;
    directCamera->position = SbVec3f(0.0f, 0.0f, 10.0f);
    directCamera->height = 4.0f;
    directCamera->nearDistance = 1.0f;
    directCamera->farDistance = 20.0f;
    directRoot->addChild(directCamera);
    SoCADViewState *directViewState = new SoCADViewState;
    directViewState->softwareWireMode = SoCADViewState::SOFTWARE_WIRE_FAST;
    directRoot->addChild(directViewState);
    SoCADAssembly *directAssembly = new SoCADAssembly;
    Obol::WireRep directWire;
    directWire.segmentPoints = {
	SbVec3f(-1.0f, -1.0f, 0.0f), SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(1.0f, -1.0f, 0.0f), SbVec3f(1.0f, 1.0f, 0.0f),
	SbVec3f(1.0f, 1.0f, 0.0f), SbVec3f(-1.0f, 1.0f, 0.0f),
	SbVec3f(-1.0f, 1.0f, 0.0f), SbVec3f(-1.0f, -1.0f, 0.0f)
    };
    directWire.bounds.setBounds(
	SbVec3f(-1.0f, -1.0f, 0.0f), SbVec3f(1.0f, 1.0f, 0.0f));
    Obol::PartGeometryBuilder directGeometry;
    directGeometry.wire = directWire;
    const Obol::PartId directPart =
	Obol::CadIdBuilder::partId("render-direct-wire");
    if (!admitAndUpsertPart(directAssembly, directPart, directGeometry))
	FAIL("direct wire geometry publication failed");
    Obol::InstanceRecord directInstance;
    directInstance.part = directPart;
    directInstance.parent = Obol::CadIdBuilder::rootInstance();
    directInstance.childName = "render-direct-wire";
    directInstance.localToRoot.makeIdentity();
    if (!directAssembly->upsertInstanceAuto(directInstance))
	FAIL("direct wire instance publication failed");
    directRoot->addChild(directAssembly);

    SoOffscreenRenderer directRenderer(&contextManager, viewport);
    directRenderer.setComponents(SoOffscreenRenderer::RGB);
    directRenderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (!directRenderer.render(directRoot) ||
	!directAssembly->lastRenderUsedDirectSoftwareWire())
	FAIL("FAST retained wire geometry should use the direct OSMesa executor");
    const Obol::CadRenderedWork directWork =
	directAssembly->lastRenderedWork();
    const Obol::CadGpuResourceSnapshot directResources =
	directAssembly->gpuResourceSnapshot();
    if (!directWork.exact || directWork.lineCount != 4 ||
	directWork.positionCount != 8 || directWork.occurrenceCount != 1 ||
	!directResources.frameSerial)
	FAIL("direct OSMesa executor should publish exact completed CAD work");
    const uint64_t directFrameSerial = directResources.frameSerial;
    if (!directRenderer.render(directRoot) ||
	directAssembly->gpuResourceSnapshot().frameSerial <= directFrameSerial)
	FAIL("direct OSMesa executor should advance the completed frame token");
    directRoot->unref();

    /* A cut's declared positionCount is its complete index domain.  Keeping
     * richer positions resident must not make a malformed coarse prefix look
     * valid: a GL implementation is otherwise free to fetch beyond the
     * uploaded VBO prefix and produce transient vertices far out in space. */
    SoSeparator *invalidPrefixRoot = new SoSeparator;
    invalidPrefixRoot->ref();
    SoOrthographicCamera *invalidPrefixCamera = new SoOrthographicCamera;
    invalidPrefixCamera->position = SbVec3f(0.0f, 0.0f, 10.0f);
    invalidPrefixCamera->height = 4.0f;
    invalidPrefixCamera->nearDistance = 1.0f;
    invalidPrefixCamera->farDistance = 20.0f;
    invalidPrefixRoot->addChild(invalidPrefixCamera);
    SoCADViewState *invalidPrefixViewState = new SoCADViewState;
    invalidPrefixViewState->softwareWireMode =
	SoCADViewState::SOFTWARE_WIRE_FAST;
    invalidPrefixRoot->addChild(invalidPrefixViewState);
    SoCADAssembly *invalidPrefixAssembly = new SoCADAssembly;
    std::shared_ptr<Obol::TriMesh> invalidPrefixMesh(
	new Obol::TriMesh);
    invalidPrefixMesh->positions = {
	SbVec3f(-1.0f, -1.0f, 0.0f),
	SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(0.0f, 1.0f, 0.0f),
	SbVec3f(2000.0f, 2000.0f, 0.0f)
    };
    invalidPrefixMesh->indices = {0u, 1u, 3u};
    invalidPrefixMesh->bounds.setBounds(
	SbVec3f(-1.0f, -1.0f, 0.0f),
	SbVec3f(2000.0f, 2000.0f, 0.0f));
    invalidPrefixMesh->progressiveCuts.resize(1);
    invalidPrefixMesh->progressiveCuts[0].indexCount = 3u;
    invalidPrefixMesh->progressiveCuts[0].positionCount = 3u;
    invalidPrefixMesh->progressiveMinimumCut = 0u;
    invalidPrefixMesh->progressiveResidentCut = 0u;
    invalidPrefixMesh->progressiveQuantizationMinimum =
	invalidPrefixMesh->bounds.getMin();
    invalidPrefixMesh->progressiveQuantizationMaximum =
	invalidPrefixMesh->bounds.getMax();
    Obol::PartGeometryBuilder invalidPrefixGeometry;
    invalidPrefixGeometry.shaded = std::move(*invalidPrefixMesh);
    const Obol::PartId invalidPrefixPart =
	Obol::CadIdBuilder::partId("render-invalid-progressive-prefix");
    const Obol::CadGeometryValidation invalidPrefixResult =
	admitAndUpsertPart(invalidPrefixAssembly,
	    invalidPrefixPart, invalidPrefixGeometry);
    invalidPrefixRoot->addChild(invalidPrefixAssembly);
    if (invalidPrefixResult.error !=
	    Obol::CadGeometryError::InvalidProgressiveCut ||
	invalidPrefixAssembly->partCount() != 0)
	FAIL("CAD admission accepted an index beyond the active position prefix");
    invalidPrefixRoot->unref();

    /* The OSMesa quality path needs one compact position/edge-index pair,
     * not an additional triangle VBO copy which no software wire executor
     * submits.  This invariant is material for multi-gigabyte wire scenes. */
    SoSeparator *compactWireRoot = new SoSeparator;
    compactWireRoot->ref();
    SoOrthographicCamera *compactWireCamera = new SoOrthographicCamera;
    compactWireCamera->position = SbVec3f(0.0f, 0.0f, 10.0f);
    compactWireCamera->height = 4.0f;
    compactWireCamera->nearDistance = 1.0f;
    compactWireCamera->farDistance = 20.0f;
    compactWireRoot->addChild(compactWireCamera);
    SoCADViewState *compactWireViewState = new SoCADViewState;
    compactWireViewState->softwareWireMode =
	SoCADViewState::SOFTWARE_WIRE_QUALITY;
    compactWireRoot->addChild(compactWireViewState);
    SoCADAssembly *compactWireAssembly = new SoCADAssembly;
    Obol::TriMesh compactWireMesh;
    compactWireMesh.positions = {
	SbVec3f(-1.0f, -1.0f, 0.0f), SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(1.0f, 1.0f, 0.0f), SbVec3f(-1.0f, 1.0f, 0.0f)
    };
    compactWireMesh.indices = {0u, 1u, 2u, 0u, 2u, 3u};
    compactWireMesh.bounds.setBounds(
	SbVec3f(-1.0f, -1.0f, 0.0f), SbVec3f(1.0f, 1.0f, 0.0f));
    Obol::PartGeometryBuilder compactWireSource;
    compactWireSource.shaded = compactWireMesh;
    const Obol::CadGeometryAdmission compactWireSourceAdmission =
	Obol::cadAdmitPartGeometry(std::move(compactWireSource));
    if (!compactWireSourceAdmission)
	FAIL("compact wire source admission failed");
    Obol::WireRep compactWire;
    compactWire.triangleEdgeGeometry =
	compactWireSourceAdmission.geometry.shared();
    compactWire.triangleEdgeSegmentCount = compactWireMesh.indices.size();
    compactWire.bounds = compactWireMesh.bounds;
    Obol::PartGeometryBuilder compactWireGeometry;
    compactWireGeometry.wire = std::move(compactWire);
    const Obol::PartId compactWirePart =
	Obol::CadIdBuilder::partId("render-compact-derived-wire");
    if (!admitAndUpsertPart(compactWireAssembly,
	    compactWirePart, compactWireGeometry))
	FAIL("compact wire geometry publication failed");
    Obol::InstanceRecord compactWireInstance;
    compactWireInstance.part = compactWirePart;
    compactWireInstance.parent = Obol::CadIdBuilder::rootInstance();
    compactWireInstance.childName = "render-compact-derived-wire";
    compactWireInstance.localToRoot.makeIdentity();
    if (!compactWireAssembly->upsertInstanceAuto(compactWireInstance))
	FAIL("compact wire instance publication failed");
    compactWireRoot->addChild(compactWireAssembly);

    SoOffscreenRenderer compactWireRenderer(&contextManager, viewport);
    compactWireRenderer.setComponents(SoOffscreenRenderer::RGB);
    compactWireRenderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (!compactWireRenderer.render(compactWireRoot) ||
	compactWireAssembly->lastRenderUsedDirectSoftwareWire() ||
	compactWireAssembly->lastRenderedWork().lineCount != 6u)
	FAIL("compact derived-wire OSMesa quality render failed");
    const Obol::CadGpuResourceSnapshot compactWireResources =
	compactWireAssembly->gpuResourceSnapshot();
    const size_t expectedCompactWireBytes =
	4u * 3u * sizeof(float) + 12u * sizeof(uint32_t);
    if (compactWireResources.ordinaryPartBufferBytes !=
	expectedCompactWireBytes)
	FAIL("OSMesa derived wire retained an unused triangle VBO copy");
    compactWireRoot->unref();

    /* Wire and shaded streams are independently authored channels.  A
     * stable shaded lineage must not accidentally certify a changed wire
     * prefix merely because both live in one PartGeometry generation. */
    SoSeparator *lineageRoot = new SoSeparator;
    lineageRoot->ref();
    SoOrthographicCamera *lineageCamera = new SoOrthographicCamera;
    lineageCamera->position = SbVec3f(0.0f, 0.0f, 10.0f);
    lineageCamera->height = 4.0f;
    lineageCamera->nearDistance = 1.0f;
    lineageCamera->farDistance = 20.0f;
    lineageRoot->addChild(lineageCamera);
    SoCADViewState *lineageViewState = new SoCADViewState;
    lineageViewState->softwareWireMode =
	SoCADViewState::SOFTWARE_WIRE_QUALITY;
    lineageViewState->drawMode = SoCADViewState::SHADED_WITH_EDGES;
    lineageRoot->addChild(lineageViewState);
    SoCADAssembly *lineageAssembly = new SoCADAssembly;
    Obol::PartGeometryBuilder lineageGeometry;
    Obol::TriMesh lineageShaded;
    lineageShaded.positions = {
	SbVec3f(-1.0f, -1.0f, 0.0f), SbVec3f(1.0f, -1.0f, 0.0f),
	SbVec3f(0.0f, 1.0f, 0.0f)
    };
    lineageShaded.indices = {0u, 1u, 2u};
    lineageShaded.bounds.setBounds(
	SbVec3f(-1.0f, -1.0f, 0.0f), SbVec3f(1.0f, 1.0f, 0.0f));
    lineageShaded.progressiveCuts.resize(1);
    lineageShaded.progressiveCuts[0].indexCount = 3u;
    lineageShaded.progressiveCuts[0].positionCount = 3u;
    lineageShaded.progressiveMinimumCut = 0u;
    lineageShaded.progressiveResidentCut = 0u;
    lineageShaded.progressiveLineage = UINT64_C(0x534841444544);
    lineageGeometry.shaded = lineageShaded;
    Obol::WireRep lineageWire;
    lineageWire.segmentPoints = {
	SbVec3f(-1.0f, -0.75f, 0.05f), SbVec3f(1.0f, -0.75f, 0.05f)
    };
    lineageWire.bounds.setBounds(
	SbVec3f(-1.0f, -1.0f, 0.05f), SbVec3f(1.0f, 1.0f, 0.05f));
    lineageWire.progressiveCuts.resize(1);
    lineageWire.progressiveCuts[0].segmentCount = 1u;
    lineageWire.progressiveMinimumCut = 0u;
    lineageWire.progressiveResidentCut = 0u;
    lineageWire.progressiveLineage = UINT64_C(0x574952453031);
    lineageGeometry.wire = lineageWire;
    const Obol::PartId lineagePart =
	Obol::CadIdBuilder::partId("render-independent-channel-lineage");
    if (!admitAndUpsertPart(lineageAssembly, lineagePart, lineageGeometry))
	FAIL("independent channel geometry publication failed");
    Obol::InstanceRecord lineageInstance;
    lineageInstance.part = lineagePart;
    lineageInstance.parent = Obol::CadIdBuilder::rootInstance();
    lineageInstance.childName = "render-independent-channel-lineage";
    lineageInstance.localToRoot.makeIdentity();
    lineageInstance.lodCut = 0u;
    if (!lineageAssembly->upsertInstanceAuto(lineageInstance))
	FAIL("independent channel instance publication failed");
    lineageRoot->addChild(lineageAssembly);

    SoOffscreenRenderer lineageRenderer(&contextManager, viewport);
    lineageRenderer.setComponents(SoOffscreenRenderer::RGB);
    lineageRenderer.setBackgroundColor(SbColor(0.0f, 0.0f, 0.0f));
    if (!lineageRenderer.render(lineageRoot))
	FAIL("independent channel lineage setup render failed");
    const Obol::CadGpuResourceSnapshot coarseLineageResources =
	lineageAssembly->gpuResourceSnapshot();
    lineageWire.segmentPoints = {
	SbVec3f(0.0f, -1.0f, 0.05f), SbVec3f(0.0f, 1.0f, 0.05f),
	SbVec3f(-1.0f, 0.75f, 0.05f), SbVec3f(1.0f, 0.75f, 0.05f)
    };
    lineageWire.progressiveCuts[0].segmentCount = 2u;
    lineageWire.progressiveLineage = UINT64_C(0x574952453032);
    lineageGeometry.wire = lineageWire;
    if (!admitAndUpsertPart(lineageAssembly, lineagePart, lineageGeometry))
	FAIL("independent channel replacement publication failed");
    if (!lineageRenderer.render(lineageRoot) ||
	lineageAssembly->lastRenderedWork().lineCount != 2u)
	FAIL("changed wire lineage did not render its replacement stream");
    const Obol::CadGpuResourceSnapshot richLineageResources =
	lineageAssembly->gpuResourceSnapshot();
    const uint64_t richWireBytes = 4u * 3u * sizeof(float);
    if (richLineageResources.ordinaryPartFullUploadBytes <
	    coarseLineageResources.ordinaryPartFullUploadBytes + richWireBytes ||
	richLineageResources.ordinaryPartLineageReplacementCount <=
	    coarseLineageResources.ordinaryPartLineageReplacementCount)
	FAIL("stable shaded lineage incorrectly certified a changed wire prefix");
    lineageRoot->unref();

    SoSeparator *pointRoot = new SoSeparator;
    pointRoot->ref();
    SoOrthographicCamera *pointCamera = new SoOrthographicCamera;
    pointCamera->position = SbVec3f(0.0f, 0.0f, 10.0f);
    pointCamera->height = 6.0f;
    pointCamera->nearDistance = 1.0f;
    pointCamera->farDistance = 20.0f;
    pointRoot->addChild(pointCamera);
    SoCADAssembly *pointAssembly = new SoCADAssembly;
    Obol::PartGeometryBuilder pointGeometry;
    Obol::PointRep points;
    points.positions = {SbVec3f(-1.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f)};
    points.pointIds = {7u, 11u};
    points.colorValid = {1u, 1u};
    points.colors = {SbColor(1.0f, 0.0f, 0.0f),
	SbColor(0.0f, 1.0f, 0.0f)};
    points.bounds.setBounds(SbVec3f(-1.0f, 0.0f, 0.0f),
	SbVec3f(1.0f, 0.0f, 0.0f));
    pointGeometry.points = points;
    const Obol::PartId pointPart = Obol::CadIdBuilder::partId("render-points");
    if (!admitAndUpsertPart(pointAssembly, pointPart, pointGeometry))
	FAIL("point geometry publication failed");
    Obol::InstanceRecord pointInstance;
    pointInstance.part = pointPart;
    pointInstance.parent = Obol::CadIdBuilder::rootInstance();
    pointInstance.childName = "render-points";
    pointInstance.localToRoot.makeIdentity();
    pointInstance.style.lineWidth = 9.0f;
    if (!pointAssembly->upsertInstanceAuto(pointInstance))
	FAIL("point instance publication failed");
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

    BOBOLSoftwareLineContextManager softwareManager;
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

    BOBOLDiagnosticContextManager diagnosticManager;
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
