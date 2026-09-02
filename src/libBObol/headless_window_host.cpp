/*      H E A D L E S S _ W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BHeadlessWindowHost.h"
#include "BObol/BInit.h"
#include "BObol/BViewController.h"

#include "bu/malloc.h"

#include <Inventor/SoDB.h>
#include <Inventor/actions/SoSearchAction.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoOrthographicCamera.h>

#include <cstring>
#include <new>
#include <vector>

struct BObolHeadlessWindowHostPrivate {
    BObolHeadlessWindowHostPrivate(void) :
	pollRate(0),
	contextManager(bobol_headless_context_manager()),
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
    SoDB::ContextManager *contextManager;
    SbColor backgroundColor;
    int outputComponents;
    std::vector<unsigned char> frame;
    unsigned int frameWidth;
    unsigned int frameHeight;
    int frameComponents;
    int renderCount;
    SbString lastRenderReason;
};

static BObolWindowDesc
headless_desc(const BObolWindowDesc *input)
{
    BObolWindowDesc desc;
    if (input) {
	desc = *input;
    } else {
	desc.mode = BOBOL_WINDOW_HEADLESS;
	desc.backend = BOBOL_WINDOW_BACKEND_OFFSCREEN;
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
    desc.mode = BOBOL_WINDOW_HEADLESS;
    desc.backend = BOBOL_WINDOW_BACKEND_OFFSCREEN;
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
ensure_headless_camera(BObolViewController *controller)
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

static int
valid_components(int components)
{
    return components >= 1 && components <= 4 ? components : 3;
}

BObolHeadlessWindowHost::BObolHeadlessWindowHost(void) :
    hp(new BObolHeadlessWindowHostPrivate)
{
}

BObolHeadlessWindowHost::~BObolHeadlessWindowHost(void)
{
    this->close();
    delete this->hp;
    this->hp = NULL;
}

int
BObolHeadlessWindowHost::open(const BObolWindowDesc *desc)
{
    BObolViewController *controller = this->getController();
    if (!controller || !this->hp->contextManager)
	return -1;
    controller->setRenderContextManager(this->hp->contextManager);
    if (controller->getRenderContextManager() != this->hp->contextManager)
	return -1;

    BObolWindowDesc actual = headless_desc(desc);
    if (BObolWindowHost::open(&actual) != 0)
	return -1;
    return ensure_headless_camera(controller);
}

void
BObolHeadlessWindowHost::close(void)
{
    BObolViewController *controller = this->getController();
    if (controller && controller->getRenderContextManager() ==
	    this->hp->contextManager)
	controller->setRenderContextManager(NULL);
    BObolWindowHost::close();
}

int
BObolHeadlessWindowHost::poll(const BObolInputProfile *UNUSED(profile))
{
    return this->renderPending();
}

long
BObolHeadlessWindowHost::pollRate(void) const
{
    return this->hp->pollRate;
}

void
BObolHeadlessWindowHost::setPollRate(long rate)
{
    this->hp->pollRate = rate > 0 ? rate : 0;
}

void
BObolHeadlessWindowHost::setContextManager(SoDB::ContextManager *manager)
{
    SoDB::ContextManager *oldManager = this->hp->contextManager;
    this->hp->contextManager = manager;
    BObolViewController *controller = this->getController();
    if (controller && (this->isOpen() ||
	    controller->getRenderContextManager() == oldManager))
	controller->setRenderContextManager(manager);
}

SoDB::ContextManager *
BObolHeadlessWindowHost::getContextManager(void) const
{
    return this->hp->contextManager;
}

int
BObolHeadlessWindowHost::renderPending(void)
{
    BObolViewController *controller = this->getController();
    if (!this->isOpen() || !controller)
	return -1;
    if (!controller->isRenderRequested())
	return 0;
    if (ensure_headless_camera(controller) != 0)
	return -1;

    const SbString requestedReason = controller->getRenderReason();

    SoDB::ContextManager *manager = this->hp->contextManager;
    if (!manager)
	return -1;
    if (controller->getRenderContextManager() != manager)
	return -1;
    const SbViewportRegion &region = controller->getViewportRegion();
    const SbVec2s windowSize = region.getWindowSize();
    if (windowSize[0] <= 0 || windowSize[1] <= 0)
	return -1;

    unsigned char *pixels = NULL;
    BObolProgressiveStatus progressiveStatus;
    const int components = valid_components(this->hp->outputComponents);
    const SbBool alpha = components == 2 || components == 4 ? TRUE : FALSE;
    if (controller->renderToImage(&pixels, FALSE, alpha,
	    NULL, manager, &progressiveStatus) != BRLCAD_OK || !pixels) {
	if (pixels)
	    bu_free(pixels, "headless endpoint interrupted frame");
	return -1;
    }

    const unsigned int width = (unsigned int)windowSize[0];
    const unsigned int height = (unsigned int)windowSize[1];
    const size_t rowPixels = (size_t)width * (size_t)components;
    if (width == 0 || height == 0 || rowPixels / (size_t)components != (size_t)width)
	return -1;
    const size_t byteCount = rowPixels * (size_t)height;
    if (height != 0 && byteCount / (size_t)height != rowPixels)
	return -1;

    const int sourceComponents = alpha ? 4 : 3;
    if (components == sourceComponents) {
	this->hp->frame.assign(pixels, pixels + byteCount);
    } else {
	this->hp->frame.resize(byteCount);
	const size_t pixelCount = (size_t)width * (size_t)height;
	for (size_t i = 0; i < pixelCount; i++) {
	    const unsigned char *source = pixels + i * sourceComponents;
	    /* Match the conventional integer Rec. 601 luminance conversion.
	     * renderToImage owns the canonical RGB(A) capture; monochrome output
	     * is a host transport format, not a second scene traversal. */
	    const unsigned char luminance = (unsigned char)(
		(77u * source[0] + 150u * source[1] + 29u * source[2] + 128u) >> 8);
	    unsigned char *target = &this->hp->frame[i * components];
	    target[0] = luminance;
	    if (components == 2)
		target[1] = source[3];
	}
    }
    bu_free(pixels, "headless endpoint frame");
    this->hp->frameWidth = width;
    this->hp->frameHeight = height;
    this->hp->frameComponents = components;

    /* renderToImage owns serial-checked request retirement.  Consuming here
     * would also clear a newer result/frame request published while the image
     * was rendering, leaving that work invisible until unrelated input. */
    controller->noteFramePresented();
    this->hp->lastRenderReason = requestedReason;
    if (controller->hasProgressiveWorkPending())
	controller->requestLodCapacityRender("progressive-pending");
    this->hp->renderCount++;
    return 1;
}

void
BObolHeadlessWindowHost::setBackgroundColor(const SbColor &color)
{
    this->hp->backgroundColor = color;
    if (this->getController())
	this->getController()->setBackgroundColors(color, color);
}

const SbColor &
BObolHeadlessWindowHost::getBackgroundColor(void) const
{
    return this->hp->backgroundColor;
}

void
BObolHeadlessWindowHost::setOutputComponents(int components)
{
    this->hp->outputComponents = valid_components(components);
}

int
BObolHeadlessWindowHost::getOutputComponents(void) const
{
    return this->hp->outputComponents;
}

const unsigned char *
BObolHeadlessWindowHost::getLastFrameBuffer(void) const
{
    return this->hp->frame.empty() ? NULL : &this->hp->frame[0];
}

size_t
BObolHeadlessWindowHost::getLastFrameBufferSize(void) const
{
    return this->hp->frame.size();
}

unsigned int
BObolHeadlessWindowHost::getLastFrameWidth(void) const
{
    return this->hp->frameWidth;
}

unsigned int
BObolHeadlessWindowHost::getLastFrameHeight(void) const
{
    return this->hp->frameHeight;
}

int
BObolHeadlessWindowHost::getLastFrameComponents(void) const
{
    return this->hp->frameComponents;
}

int
BObolHeadlessWindowHost::getRenderCount(void) const
{
    return this->hp->renderCount;
}

const SbString &
BObolHeadlessWindowHost::getLastRenderReason(void) const
{
    return this->hp->lastRenderReason;
}

static BObolWindowDesc
headless_factory_window_desc(const struct bobol_host_desc *desc)
{
    BObolWindowDesc actual;
    actual.mode = BOBOL_WINDOW_HEADLESS;
    actual.backend = BOBOL_WINDOW_BACKEND_OFFSCREEN;
    actual.width = desc && desc->width ? desc->width : 1;
    actual.height = desc && desc->height ? desc->height : 1;
    actual.title = desc && desc->title ? desc->title : "BRL-CAD Obol Headless";
    actual.display = desc && desc->display ? desc->display : "";
    actual.nativeIdHint = desc && desc->native_id_hint ?
	desc->native_id_hint : "";
    actual.visible = FALSE;
    return actual;
}

static int
headless_factory_probe(const struct bobol_host_desc *desc,
	void *UNUSED(user_data))
{
    return !desc || desc->mode == BOBOL_HOST_MODE_HEADLESS ||
	desc->mode == BOBOL_HOST_MODE_DIAGNOSTIC;
}

static void *
headless_factory_create(const struct bobol_host_desc *UNUSED(desc),
	void *UNUSED(user_data))
{
    BObolHeadlessWindowHost *host =
	new (std::nothrow) BObolHeadlessWindowHost;
    if (host)
	host->setContextManager(bobol_headless_context_manager());
    return host;
}

static void
headless_factory_destroy(void *instance, void *UNUSED(user_data))
{
    delete static_cast<BObolHeadlessWindowHost *>(instance);
}

static int
headless_factory_bind(void *instance, void *controller,
	void *UNUSED(user_data))
{
    BObolHeadlessWindowHost *host =
	static_cast<BObolHeadlessWindowHost *>(instance);
    if (!host)
	return 0;
    BObolViewController *newController =
	static_cast<BObolViewController *>(controller);
    BObolViewController *oldController = host->getController();
    if (oldController && oldController != newController &&
	oldController->getRenderContextManager() == host->getContextManager())
	oldController->setRenderContextManager(NULL);
    if (oldController != newController)
	host->attachController(newController, FALSE);
    if (!newController)
	return 1;
    newController->setRenderContextManager(host->getContextManager());
    return host->getContextManager() &&
	newController->getRenderContextManager() == host->getContextManager();
}

static int
headless_factory_open(void *instance, const struct bobol_host_desc *desc,
	void *UNUSED(user_data))
{
    BObolHeadlessWindowHost *host =
	static_cast<BObolHeadlessWindowHost *>(instance);
    BObolWindowDesc actual = headless_factory_window_desc(desc);
    return host && host->open(&actual) == 0 ? 1 : 0;
}

static void
headless_factory_close(void *instance, void *UNUSED(user_data))
{
    BObolHeadlessWindowHost *host =
	static_cast<BObolHeadlessWindowHost *>(instance);
    if (host)
	host->close();
}

static int
headless_factory_request_frame(void *instance, const char *reason,
	void *UNUSED(user_data))
{
    BObolHeadlessWindowHost *host =
	static_cast<BObolHeadlessWindowHost *>(instance);
    if (!host || !host->getController())
	return 0;
    host->getController()->requestLodCapacityRender(reason);
    return 1;
}

static int
headless_factory_resize(void *instance, unsigned int width,
	unsigned int height, double device_pixel_ratio, void *UNUSED(user_data))
{
    BObolHeadlessWindowHost *host =
	static_cast<BObolHeadlessWindowHost *>(instance);
    if (!host || !host->getController() || !width || !height ||
	device_pixel_ratio <= 0.0)
	return 0;
    host->getController()->setViewportSize(width, height);
    return 1;
}

static int
headless_factory_capture(void *instance, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components,
	void *UNUSED(user_data))
{
    BObolHeadlessWindowHost *host =
	static_cast<BObolHeadlessWindowHost *>(instance);
    if (!host || !pixels || !size || !width || !height || !components)
	return 0;

    if (host->getController())
	host->getController()->requestLodCapacityRender("capture");
    if (host->renderPending() < 0 || !host->getLastFrameBuffer())
	return 0;

    *size = host->getLastFrameBufferSize();
    *width = host->getLastFrameWidth();
    *height = host->getLastFrameHeight();
    *components = (unsigned int)host->getLastFrameComponents();
    *pixels = static_cast<unsigned char *>(bu_malloc(*size,
	"headless endpoint capture"));
    memcpy(*pixels, host->getLastFrameBuffer(), *size);
    return 1;
}

bobol_host_factory_token_t *
bobol_headless_host_factory_register(void)
{
    static bobol_host_factory_token_t *token = []() {
	struct bobol_host_factory factory;
	memset(&factory, 0, sizeof(factory));
	factory.abi_version = BOBOL_HOST_FACTORY_ABI_VERSION;
	factory.struct_size = sizeof(factory);
	factory.name = "headless";
	factory.priority = 0;
	factory.capabilities = BOBOL_HOST_CAP_PIXEL_PRESENT |
	    BOBOL_HOST_CAP_PROGRESSIVE_PRESENT |
	    BOBOL_HOST_CAP_READBACK |
	    BOBOL_HOST_CAP_FRAMEBUFFER_PRESENT |
	    BOBOL_HOST_CAP_MULTI_VIEW;
	factory.probe = headless_factory_probe;
	factory.create = headless_factory_create;
	factory.destroy = headless_factory_destroy;
	factory.bind_controller = headless_factory_bind;
	factory.open = headless_factory_open;
	factory.close = headless_factory_close;
	factory.request_frame = headless_factory_request_frame;
	factory.resize = headless_factory_resize;
	factory.capture = headless_factory_capture;
	return bobol_host_factory_register(&factory);
    }();
    return token;
}
