/*      H E A D L E S S _ W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/headless_window_host.h"
#include "brlobol/view_controller.h"

#include <Inventor/SoDB.h>
#include <Inventor/SoOffscreenRenderer.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoOrthographicCamera.h>

#include <vector>

struct BRLObolHeadlessWindowHostPrivate {
    BRLObolHeadlessWindowHostPrivate(void) :
	pollRate(0),
	backgroundColor(0.0f, 0.0f, 0.0f),
	outputComponents(3),
	frameWidth(0),
	frameHeight(0),
	frameComponents(0),
	renderCount(0),
	lastRenderReason("")
    {
    }

    long pollRate;
    SbColor backgroundColor;
    int outputComponents;
    std::vector<unsigned char> frame;
    unsigned int frameWidth;
    unsigned int frameHeight;
    int frameComponents;
    int renderCount;
    SbString lastRenderReason;
};

static BRLObolWindowDesc
headless_desc(const BRLObolWindowDesc *input)
{
    BRLObolWindowDesc desc;
    if (input) {
	desc = *input;
    } else {
	desc.mode = BRLOBOL_WINDOW_HEADLESS;
	desc.backend = BRLOBOL_WINDOW_BACKEND_OFFSCREEN;
	desc.width = 1;
	desc.height = 1;
	desc.title = "BRL-CAD Obol Headless";
	desc.display = "";
	desc.nativeIdHint = "";
	desc.visible = FALSE;
    }

    if (desc.width == 0)
	desc.width = 1;
    if (desc.height == 0)
	desc.height = 1;
    desc.mode = BRLOBOL_WINDOW_HEADLESS;
    desc.backend = BRLOBOL_WINDOW_BACKEND_OFFSCREEN;
    desc.visible = FALSE;
    return desc;
}

static SoCamera *
find_camera(SoNode *root)
{
    if (!root)
	return NULL;

    SoSearchAction search;
    search.setType(SoCamera::getClassTypeId());
    search.setInterest(SoSearchAction::FIRST);
    search.apply(root);

    SoPath *path = search.getPath();
    if (!path || path->getLength() < 1)
	return NULL;

    SoNode *tail = path->getTail();
    if (!tail || !tail->isOfType(SoCamera::getClassTypeId()))
	return NULL;
    return static_cast<SoCamera *>(tail);
}

static int
ensure_headless_camera(BRLObolViewController *controller)
{
    if (!controller)
	return -1;
    if (controller->getCamera())
	return 0;

    SoNode *root = controller->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return -1;

    SoCamera *camera = find_camera(root);
    if (!camera) {
	SoOrthographicCamera *newCamera = new SoOrthographicCamera;
	newCamera->position.setValue(0.0f, 0.0f, 10.0f);
	newCamera->height = 2.0f;
	newCamera->nearDistance = 1.0f;
	newCamera->farDistance = 20.0f;
	newCamera->focalDistance = 10.0f;
	static_cast<SoGroup *>(root)->addChild(newCamera);
	camera = newCamera;
    }

    controller->setCamera(camera);
    return 0;
}

static SoOffscreenRenderer::Components
renderer_components(int components)
{
    switch (components) {
	case 1:
	    return SoOffscreenRenderer::LUMINANCE;
	case 2:
	    return SoOffscreenRenderer::LUMINANCE_TRANSPARENCY;
	case 4:
	    return SoOffscreenRenderer::RGB_TRANSPARENCY;
	case 3:
	default:
	    return SoOffscreenRenderer::RGB;
    }
}

static int
valid_components(int components)
{
    return components >= 1 && components <= 4 ? components : 3;
}

BRLObolHeadlessWindowHost::BRLObolHeadlessWindowHost(void) :
    hp(new BRLObolHeadlessWindowHostPrivate)
{
}

BRLObolHeadlessWindowHost::~BRLObolHeadlessWindowHost(void)
{
    delete this->hp;
    this->hp = NULL;
}

int
BRLObolHeadlessWindowHost::open(const BRLObolWindowDesc *desc)
{
    BRLObolWindowDesc actual = headless_desc(desc);
    if (BRLObolWindowHost::open(&actual) != 0)
	return -1;
    return ensure_headless_camera(this->getController());
}

int
BRLObolHeadlessWindowHost::poll(const BRLObolInputProfile *UNUSED(profile))
{
    return this->renderPending();
}

long
BRLObolHeadlessWindowHost::pollRate(void) const
{
    return this->hp->pollRate;
}

void
BRLObolHeadlessWindowHost::setPollRate(long rate)
{
    this->hp->pollRate = rate > 0 ? rate : 0;
}

int
BRLObolHeadlessWindowHost::renderPending(void)
{
    BRLObolViewController *controller = this->getController();
    if (!this->isOpen() || !controller)
	return -1;
    if (!controller->isRenderRequested())
	return 0;
    if (ensure_headless_camera(controller) != 0)
	return -1;

    SoNode *root = controller->getRenderRoot();
    if (!root)
	root = controller->getSceneRoot();
    if (!root)
	return -1;

    const SbViewportRegion &region = controller->getViewportRegion();
    SbVec2s windowSize = region.getWindowSize();
    if (windowSize[0] <= 0 || windowSize[1] <= 0)
	return -1;

    SoDB::ContextManager *manager = SoDB::getContextManager();
    if (!manager)
	return -1;

    SoOffscreenRenderer renderer(manager, region);
    renderer.setComponents(renderer_components(this->hp->outputComponents));
    renderer.setBackgroundColor(this->hp->backgroundColor);
    if (!renderer.render(root))
	return -1;

    unsigned char *pixels = renderer.getBuffer();
    if (!pixels)
	return -1;

    const unsigned int width = (unsigned int)windowSize[0];
    const unsigned int height = (unsigned int)windowSize[1];
    const int components = valid_components(this->hp->outputComponents);
    const size_t rowPixels = (size_t)width * (size_t)components;
    if (width == 0 || height == 0 || rowPixels / (size_t)components != (size_t)width)
	return -1;
    const size_t byteCount = rowPixels * (size_t)height;
    if (height != 0 && byteCount / (size_t)height != rowPixels)
	return -1;

    this->hp->frame.assign(pixels, pixels + byteCount);
    this->hp->frameWidth = width;
    this->hp->frameHeight = height;
    this->hp->frameComponents = components;

    SbString reason;
    controller->consumeRenderRequest(&reason);
    this->hp->lastRenderReason = reason;
    this->hp->renderCount++;
    return 1;
}

void
BRLObolHeadlessWindowHost::setBackgroundColor(const SbColor &color)
{
    this->hp->backgroundColor = color;
}

const SbColor &
BRLObolHeadlessWindowHost::getBackgroundColor(void) const
{
    return this->hp->backgroundColor;
}

void
BRLObolHeadlessWindowHost::setOutputComponents(int components)
{
    this->hp->outputComponents = valid_components(components);
}

int
BRLObolHeadlessWindowHost::getOutputComponents(void) const
{
    return this->hp->outputComponents;
}

const unsigned char *
BRLObolHeadlessWindowHost::getLastFrameBuffer(void) const
{
    return this->hp->frame.empty() ? NULL : &this->hp->frame[0];
}

size_t
BRLObolHeadlessWindowHost::getLastFrameBufferSize(void) const
{
    return this->hp->frame.size();
}

unsigned int
BRLObolHeadlessWindowHost::getLastFrameWidth(void) const
{
    return this->hp->frameWidth;
}

unsigned int
BRLObolHeadlessWindowHost::getLastFrameHeight(void) const
{
    return this->hp->frameHeight;
}

int
BRLObolHeadlessWindowHost::getLastFrameComponents(void) const
{
    return this->hp->frameComponents;
}

int
BRLObolHeadlessWindowHost::getRenderCount(void) const
{
    return this->hp->renderCount;
}

const SbString &
BRLObolHeadlessWindowHost::getLastRenderReason(void) const
{
    return this->hp->lastRenderReason;
}
