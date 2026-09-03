/*              D I S P L A Y _ E N D P O I N T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/malloc.h"
#include "bu/str.h"

#include "BObol/BDisplayEndpoint.h"
#include "BObol/BImageSource.h"
#include "BObol/BInit.h"
#include "BObol/BRtRender.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BViewController.h"
#include "BObol/BViewportImage.h"
#include "BObol/BWindowHost.h"
#include "identity_counter_private.h"

#include "imgstream/stream.h"

#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/nodes/SoGroup.h>

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <mutex>
#include <new>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

struct EndpointRtState;

struct bobol_display_endpoint {
    bobol_display_endpoint(void) :
	controller(NULL),
	owns_controller(false),
	host(NULL),
	owns_host(false),
	factory(NULL),
	factory_instance(NULL),
	engine(BOBOL_RENDER_ENGINE_AUTO),
	host_mode(BOBOL_HOST_MODE_HEADLESS),
	width(1),
	height(1),
	device_pixel_ratio(1.0),
	visible(false),
	vsync(true),
	framebuffer_capture_callback(NULL),
	framebuffer_capture_user_data(NULL),
	factory_frame_callback_bound(false),
	property_get_callback(NULL),
	property_set_callback(NULL),
	property_user_data(NULL),
	diagnostic_revision(0),
	diagnostic_visited_sources(0),
	diagnostic_realized_sources(0),
	diagnostic_failed_sources(0),
	diagnostic_progressive_pending(false),
	rt(NULL)
    {
	input.setProfile(bobol_input_default_view_profile());
    }

    BObolViewController *controller;
    bool owns_controller;
    BObolWindowHost *host;
    bool owns_host;
    bobol_host_factory_token_t *factory;
    void *factory_instance;
    enum bobol_render_engine engine;
    enum bobol_host_mode host_mode;
    unsigned int width;
    unsigned int height;
    double device_pixel_ratio;
    std::string title;
    bool visible;
    bool vsync;
    bobol_endpoint_framebuffer_capture_callback framebuffer_capture_callback;
    void *framebuffer_capture_user_data;
    bool factory_frame_callback_bound;
    BObolInputContext input;
    bv_display_property_get_callback property_get_callback;
    bv_display_property_set_callback property_set_callback;
    void *property_user_data;
    mutable std::string supported_engines;
    std::string diagnostic_summary;
    uint64_t diagnostic_revision;
    uint64_t diagnostic_visited_sources;
    uint64_t diagnostic_realized_sources;
    uint64_t diagnostic_failed_sources;
    bool diagnostic_progressive_pending;
    EndpointRtState *rt;

    void graphicalRenderingSet(bool enabled)
    {
	if (this->controller)
	    this->controller->setEndpointGraphicalRenderingEnabled(
		enabled ? TRUE : FALSE);
    }

    void rendererPerformanceChanged(void)
    {
	if (this->controller)
	    this->controller->invalidateRendererPerformanceHistory();
    }

    void retireWork(void)
    {
	if (this->controller)
	    this->controller->retireDisplayEndpointWork();
    }

    void resumeWork(void)
    {
	if (this->controller)
	    this->controller->resumeDisplayEndpointWork();
    }
};

/* librt owns no Coin or GUI state.  The worker stores complete RGB frames
 * here; the controller's presentation hook applies them on the host thread. */
struct EndpointRtState {
    EndpointRtState(void) :
	cancelled(false),
	generation(0),
	readyWidth(0),
	readyHeight(0),
	ready(false),
	presentedWidth(0),
	presentedHeight(0),
	workers(0),
	samples(1),
	previewScale(4),
	frameBudgetMilliseconds(33),
	interactiveQuality(true),
	presentationLayer(SoBRLViewportImage::INTERLAY),
	stream(NULL),
	source(NULL),
	viewport(NULL),
	presentationRoot(NULL),
	presentationAttached(false)
    {
    }

    BObolRtRenderer renderer;
    std::thread worker;
    std::mutex readyMutex;
    std::atomic_bool cancelled;
    std::atomic<uint64_t> generation;
    std::vector<unsigned char> readyPixels;
    BObolRtRenderPlanes readyPlanes;
    unsigned int readyWidth;
    unsigned int readyHeight;
    bool ready;
    BObolRtRenderPlanes presentedPlanes;
    unsigned int presentedWidth;
    unsigned int presentedHeight;
    unsigned int workers;
    unsigned int samples;
    unsigned int previewScale;
    unsigned int frameBudgetMilliseconds;
    bool interactiveQuality;
    int presentationLayer;
    imgstream_t *stream;
    SoBRLImageSource *source;
    SoBRLViewportImage *viewport;
    SoGroup *presentationRoot;
    bool presentationAttached;
};

static bool valid_engine(enum bobol_render_engine engine);
static int endpoint_rt_start(bobol_display_endpoint_t *endpoint);
static int endpoint_diagnostic_refresh(bobol_display_endpoint_t *endpoint);

static void
endpoint_frame_requested(void *user_data, const char *reason)
{
    bobol_display_endpoint_t *endpoint =
	static_cast<bobol_display_endpoint_t *>(user_data);
    if (!endpoint || !endpoint->factory || !endpoint->factory_instance)
	return;
    (void)bobol_host_factory_instance_request_frame(endpoint->factory,
	endpoint->factory_instance, reason);
}

static int
endpoint_input_dispatch(void *user_data, const BObolInputEvent *event)
{
    return bobol_display_endpoint_input_dispatch(
	static_cast<bobol_display_endpoint_t *>(user_data), event);
}

static void
endpoint_dimensions_refresh(bobol_display_endpoint_t *endpoint)
{
    if (!endpoint || !endpoint->factory || !endpoint->factory_instance)
	return;
    unsigned int width = 0;
    unsigned int height = 0;
    double device_pixel_ratio = 0.0;
    if (!bobol_host_factory_instance_dimensions(endpoint->factory,
	endpoint->factory_instance, &width, &height, &device_pixel_ratio))
	return;
    if (width)
	endpoint->width = width;
    if (height)
	endpoint->height = height;
    if (device_pixel_ratio > 0.0)
	endpoint->device_pixel_ratio = device_pixel_ratio;
}

static void
endpoint_rt_stop(bobol_display_endpoint_t *endpoint)
{
    if (!endpoint || !endpoint->rt)
	return;
    EndpointRtState *state = endpoint->rt;
    state->cancelled.store(true, std::memory_order_release);
    if (state->worker.joinable())
	state->worker.join();
}

static void
endpoint_rt_remove_presentation(bobol_display_endpoint_t *endpoint)
{
    if (!endpoint || !endpoint->rt)
	return;
    EndpointRtState *state = endpoint->rt;
    if (state->presentationAttached && state->presentationRoot &&
	state->viewport) {
	for (int i = 0; i < state->presentationRoot->getNumChildren(); ++i) {
	    if (state->presentationRoot->getChild(i) == state->viewport) {
		state->presentationRoot->removeChild(i);
		break;
	    }
	}
    }
    state->presentationAttached = false;
    state->presentationRoot = NULL;
    if (state->viewport) {
	state->viewport->unref();
	state->viewport = NULL;
    }
    if (state->source) {
	state->source->unref();
	state->source = NULL;
    }
    if (state->stream) {
	imgstream_destroy(state->stream);
	state->stream = NULL;
    }
}

static SoGroup *
endpoint_rt_presentation_root(BObolViewController *controller, int layer)
{
    if (!controller)
	return NULL;
    switch (layer) {
	case SoBRLViewportImage::UNDERLAY:
	return controller->getFramebufferUnderlayRoot();
	case SoBRLViewportImage::OVERLAY:
	return controller->getFramebufferOverlayRoot();
	case SoBRLViewportImage::INTERLAY:
	return controller->getFramebufferInterlayRoot();
	default:
	return NULL;
    }
}

static const char *
endpoint_rt_presentation_layer_name(int layer)
{
    switch (layer) {
	case SoBRLViewportImage::UNDERLAY:
	return "underlay";
	case SoBRLViewportImage::OVERLAY:
	return "overlay";
	case SoBRLViewportImage::INTERLAY:
	return "interlay";
	default:
	return NULL;
    }
}

static int
endpoint_rt_presentation_layer_from_name(const char *name, int *layer)
{
    if (!name || !layer)
	return 0;
    if (bu_strcmp(name, "underlay") == 0)
	*layer = SoBRLViewportImage::UNDERLAY;
    else if (bu_strcmp(name, "interlay") == 0)
	*layer = SoBRLViewportImage::INTERLAY;
    else if (bu_strcmp(name, "overlay") == 0)
	*layer = SoBRLViewportImage::OVERLAY;
    else
	return 0;
    return 1;
}

static int
endpoint_rt_presentation_layer_apply(bobol_display_endpoint_t *endpoint)
{
    if (!endpoint || !endpoint->controller || !endpoint->rt)
	return 0;
    EndpointRtState *state = endpoint->rt;
    if (!state->viewport)
	return 1;
    SoGroup *root = endpoint_rt_presentation_root(endpoint->controller,
	state->presentationLayer);
    if (!root)
	return 0;

    if (state->presentationAttached && state->presentationRoot != root &&
	state->presentationRoot) {
	for (int i = 0; i < state->presentationRoot->getNumChildren(); ++i) {
	    if (state->presentationRoot->getChild(i) == state->viewport) {
		state->presentationRoot->removeChild(i);
		break;
	    }
	}
	state->presentationAttached = false;
    }

    state->viewport->layer = state->presentationLayer;
    if (!state->presentationAttached) {
	if (state->viewport->rebuildGeometry() != 0)
	    return 0;
	/* The retained RT image is renderer output.  Host framebuffer images in
	 * the same named layer composite after it, independent of attachment or
	 * reparent timing.  Cross-layer root ordering still takes precedence. */
	root->insertChild(state->viewport, 0);
	state->presentationRoot = root;
	state->presentationAttached = true;
    }
    return 1;
}

static int
endpoint_rt_ensure_presentation(bobol_display_endpoint_t *endpoint,
	unsigned int width, unsigned int height)
{
    if (!endpoint || !endpoint->controller || !endpoint->rt || !width ||
	!height)
	return 0;
    EndpointRtState *state = endpoint->rt;
    SoGroup *root = endpoint_rt_presentation_root(endpoint->controller,
	state->presentationLayer);
    if (!root)
	return 0;

    if (!state->stream) {
	state->stream = imgstream_create(width, height, IMGSTREAM_PIXEL_RGB8);
	if (!state->stream)
	    return 0;
	state->source = new (std::nothrow) SoBRLImageSource;
	if (state->source)
	    state->source->ref();
	state->viewport = new (std::nothrow) SoBRLViewportImage;
	if (state->viewport)
	    state->viewport->ref();
	if (!state->source || !state->viewport) {
	    endpoint_rt_remove_presentation(endpoint);
	    return 0;
	}
	state->source->imageId = "renderer::rt";
	state->source->sourceUri = "librt:retained";
	if (state->source->setStream(state->stream) != 0) {
	    endpoint_rt_remove_presentation(endpoint);
	    return 0;
	}
	state->viewport->overlayId = "renderer::rt";
	state->viewport->imageSource.setValue(state->source);
	state->viewport->layer = state->presentationLayer;
	state->viewport->anchor = SoBRLViewportImage::LOWER_LEFT;
	state->viewport->fit = SoBRLViewportImage::STRETCH;
	state->viewport->preserveAspect = FALSE;
	state->viewport->position.setValue(0.0f, 0.0f);
	state->viewport->cursorVisible = FALSE;
    }

    if (imgstream_width(state->stream) != width ||
	imgstream_height(state->stream) != height) {
	if (imgstream_resize(state->stream, width, height) != 0)
	    return 0;
    }
    state->viewport->size.setValue(static_cast<float>(width),
	static_cast<float>(height));
    state->viewport->sourceCenter.setValue(static_cast<float>(width) * 0.5f,
	static_cast<float>(height) * 0.5f);
    state->viewport->sourceZoom = 1.0f;

    return endpoint_rt_presentation_layer_apply(endpoint);
}

static bool
endpoint_rt_scene_request(const char *reason)
{
    if (!reason || !reason[0])
	return false;
    return bu_strcmp(reason, "rt-view-camera") == 0 ||
	bu_strcmp(reason, "camera") == 0 ||
	bu_strcmp(reason, "viewport") == 0 ||
	bu_strcmp(reason, "viewport-size") == 0 ||
	bu_strcmp(reason, "background") == 0 ||
	bu_strcmp(reason, "clip-bounds") == 0 ||
	bu_strcmp(reason, "lighting") == 0 ||
	bu_strcmp(reason, "transparency") == 0 ||
	bu_strcmp(reason, "scene-root") == 0 ||
	bu_strcmp(reason, "render-scene-root") == 0 ||
	bu_strcmp(reason, "view-attachment") == 0 ||
	bu_strcmp(reason, "database-source") == 0 ||
	bu_strcmp(reason, "external property") == 0;
}

static void
endpoint_rt_clear_presentation(bobol_display_endpoint_t *endpoint)
{
    if (!endpoint || !endpoint->controller || !endpoint->rt ||
	!endpoint->rt->stream || !endpoint->rt->source ||
	!endpoint->rt->viewport)
	return;
    EndpointRtState *state = endpoint->rt;
    const unsigned int width = static_cast<unsigned int>(
	imgstream_width(state->stream));
    const unsigned int height = static_cast<unsigned int>(
	imgstream_height(state->stream));
    if (!width || !height)
	return;
    const SbColor bottom = endpoint->controller->getBackgroundBottomColor();
    const SbColor top = endpoint->controller->getBackgroundTopColor();
    std::vector<unsigned char> pixels(static_cast<size_t>(width) * height * 3u,
	0u);
    for (unsigned int y = 0; y < height; ++y) {
	const float t = height > 1u ? static_cast<float>(y) /
	    static_cast<float>(height - 1u) : 0.0f;
	const SbColor background = bottom * (1.0f - t) + top * t;
	for (unsigned int x = 0; x < width; ++x) {
	    unsigned char *pixel = &pixels[(static_cast<size_t>(y) * width + x) * 3u];
	    pixel[0] = static_cast<unsigned char>(std::lround(
		std::max(0.0f, std::min(1.0f, background[0])) * 255.0f));
	    pixel[1] = static_cast<unsigned char>(std::lround(
		std::max(0.0f, std::min(1.0f, background[1])) * 255.0f));
	    pixel[2] = static_cast<unsigned char>(std::lround(
		std::max(0.0f, std::min(1.0f, background[2])) * 255.0f));
	}
    }
    if (imgstream_write(state->stream, pixels.data(),
	static_cast<size_t>(width) * 3u) == 0) {
	(void)state->source->refreshFromStream();
	(void)state->viewport->syncFromSource();
    }
}

static void
endpoint_rt_presentation_sync(void *userData)
{
    bobol_display_endpoint_t *endpoint =
	static_cast<bobol_display_endpoint_t *>(userData);
    if (!endpoint || endpoint->engine != BOBOL_RENDER_ENGINE_RT ||
	!endpoint->rt)
	return;
    if (endpoint_rt_scene_request(endpoint->controller->getRenderReason().getString())) {
	(void)endpoint_rt_start(endpoint);
	return;
    }
    EndpointRtState *state = endpoint->rt;
    std::vector<unsigned char> pixels;
    BObolRtRenderPlanes planes;
    unsigned int width = 0;
    unsigned int height = 0;
    {
	std::lock_guard<std::mutex> lock(state->readyMutex);
	if (!state->ready)
	    return;
	pixels.swap(state->readyPixels);
	planes = std::move(state->readyPlanes);
	width = state->readyWidth;
	height = state->readyHeight;
	state->ready = false;
    }
    if (!endpoint_rt_ensure_presentation(endpoint, width, height) ||
	!state->stream || pixels.size() !=
	static_cast<size_t>(width) * static_cast<size_t>(height) * 3u)
	return;
    if (imgstream_write(state->stream, pixels.data(),
	static_cast<size_t>(width) * 3u) != 0 ||
	state->source->refreshFromStream() != 0 ||
	state->viewport->syncFromSource() != 0)
	return;
    {
	std::lock_guard<std::mutex> lock(state->readyMutex);
	state->presentedPlanes = std::move(planes);
	state->presentedWidth = width;
	state->presentedHeight = height;
    }
}

static void
endpoint_rt_publish(bobol_display_endpoint_t *endpoint,
	EndpointRtState *state, uint64_t generation,
	std::vector<unsigned char> &&pixels, BObolRtRenderPlanes &&planes,
	unsigned int width, unsigned int height)
{
    if (!endpoint || !state || state->cancelled.load(std::memory_order_acquire) ||
	state->generation.load(std::memory_order_acquire) != generation)
	return;
    {
	std::lock_guard<std::mutex> lock(state->readyMutex);
	if (state->cancelled.load(std::memory_order_acquire) ||
	    state->generation.load(std::memory_order_acquire) != generation)
	    return;
	state->readyPixels = std::move(pixels);
	state->readyPlanes = std::move(planes);
	state->readyWidth = width;
	state->readyHeight = height;
	state->ready = true;
    }
    endpoint->controller->requestLodCapacityRender("rt-frame");
    endpoint_frame_requested(endpoint, "rt-frame");
}

static std::vector<unsigned char>
endpoint_rt_upscale(const std::vector<unsigned char> &pixels,
	unsigned int sourceWidth, unsigned int sourceHeight, unsigned int width,
	unsigned int height)
{
    std::vector<unsigned char> enlarged(static_cast<size_t>(width) * height * 3u,
	0u);
    if (pixels.size() != static_cast<size_t>(sourceWidth) * sourceHeight * 3u ||
	!sourceWidth || !sourceHeight)
	return enlarged;
    for (unsigned int y = 0; y < height; ++y) {
	const unsigned int sourceY = std::min(sourceHeight - 1u,
	    static_cast<unsigned int>((static_cast<uint64_t>(y) * sourceHeight) /
	    height));
	for (unsigned int x = 0; x < width; ++x) {
	    const unsigned int sourceX = std::min(sourceWidth - 1u,
		static_cast<unsigned int>((static_cast<uint64_t>(x) * sourceWidth) /
		width));
	    const size_t sourceOffset =
		(static_cast<size_t>(sourceY) * sourceWidth + sourceX) * 3u;
	    const size_t offset = (static_cast<size_t>(y) * width + x) * 3u;
	    enlarged[offset] = pixels[sourceOffset];
	    enlarged[offset + 1u] = pixels[sourceOffset + 1u];
	    enlarged[offset + 2u] = pixels[sourceOffset + 2u];
	}
    }
    return enlarged;
}

static void
endpoint_rt_worker(bobol_display_endpoint_t *endpoint, EndpointRtState *state,
	uint64_t generation, BObolRtRenderSettings settings,
	unsigned int previewScale, unsigned int frameBudgetMilliseconds,
	bool interactiveQuality)
{
    if (!endpoint || !state)
	return;
    if (interactiveQuality && previewScale > 1u) {
	BObolRtRenderSettings preview = settings;
	preview.width = std::max(1u, settings.width / previewScale);
	preview.height = std::max(1u, settings.height / previewScale);
	preview.samples = 1u;
	std::vector<unsigned char> previewPixels;
	if (state->renderer.render(preview, previewPixels, NULL,
	    &state->cancelled)) {
	    endpoint_rt_publish(endpoint, state, generation,
		endpoint_rt_upscale(previewPixels, preview.width, preview.height,
		settings.width, settings.height), BObolRtRenderPlanes(),
		settings.width, settings.height);
	}
    }
    if (state->cancelled.load(std::memory_order_acquire))
	return;

    /* Render complete rows into one retained image.  This keeps every
     * publication self-contained while allowing a slow final image to make
     * visible progress at a caller-selected cadence.  The row batch adapts to
     * its measured cost; a single row remains the lower bound when a scene is
     * expensive enough that it alone exceeds the requested budget. */
    std::vector<unsigned char> pixels;
	BObolRtRenderPlanes planes;
    unsigned int firstRow = 0;
    unsigned int rowBatch = std::max(1u, std::min(settings.workers, 32u));
    const unsigned int budget = std::max(1u, frameBudgetMilliseconds);
    std::chrono::steady_clock::time_point lastPublish =
	std::chrono::steady_clock::now();
    while (firstRow < settings.height &&
	!state->cancelled.load(std::memory_order_acquire)) {
	const unsigned int rows = std::min(rowBatch, settings.height - firstRow);
	const std::chrono::steady_clock::time_point batchStart =
	    std::chrono::steady_clock::now();
	if (!state->renderer.renderRowsWithPlanes(settings, pixels, planes,
		firstRow, rows, NULL, &state->cancelled))
	    return;
	firstRow += rows;
	const std::chrono::steady_clock::time_point now =
	    std::chrono::steady_clock::now();
	const unsigned int batchMilliseconds = static_cast<unsigned int>(
	    std::chrono::duration_cast<std::chrono::milliseconds>(
		now - batchStart).count());
	if (batchMilliseconds > budget && rowBatch > 1u) {
	    rowBatch = std::max(1u, rowBatch / 2u);
	} else if (batchMilliseconds < budget / 2u) {
	    rowBatch = std::min(128u, rowBatch * 2u);
	}

	const unsigned int elapsedMilliseconds = static_cast<unsigned int>(
	    std::chrono::duration_cast<std::chrono::milliseconds>(
		now - lastPublish).count());
	if (firstRow < settings.height && elapsedMilliseconds >= budget) {
	    std::vector<unsigned char> partial(pixels);
	    BObolRtRenderPlanes partialPlanes(planes);
	    endpoint_rt_publish(endpoint, state, generation, std::move(partial),
		std::move(partialPlanes), settings.width, settings.height);
	    lastPublish = now;
	}
    }
    if (!state->cancelled.load(std::memory_order_acquire))
	endpoint_rt_publish(endpoint, state, generation, std::move(pixels),
	    std::move(planes), settings.width, settings.height);
}

static int
endpoint_rt_start(bobol_display_endpoint_t *endpoint)
{
    if (!endpoint || !endpoint->controller ||
	endpoint->engine != BOBOL_RENDER_ENGINE_RT)
	return 0;
    if (!endpoint->rt)
	endpoint->rt = new (std::nothrow) EndpointRtState;
    if (!endpoint->rt)
	return 0;
    EndpointRtState *state = endpoint->rt;
    endpoint_rt_stop(endpoint);
	{
	std::lock_guard<std::mutex> lock(state->readyMutex);
	state->readyPixels.clear();
	state->readyPlanes.clear();
	state->readyWidth = 0;
	state->readyHeight = 0;
	state->ready = false;
	state->presentedPlanes.clear();
	state->presentedWidth = 0;
	state->presentedHeight = 0;
	}
    endpoint_dimensions_refresh(endpoint);
    if (!endpoint_rt_ensure_presentation(endpoint, endpoint->width,
	endpoint->height)) {
	return 0;
	}
	endpoint_rt_clear_presentation(endpoint);
	endpoint->controller->requestLodCapacityRender("rt-restart");
	endpoint_frame_requested(endpoint, "rt-restart");
	/* Renderer policy is valid before GED has supplied a camera.  Keep the
	 * opaque presentation surface installed and start work on the next view
	 * synchronization rather than reporting a supported engine as unavailable. */
    if (!state->renderer.synchronize(endpoint->controller))
	return 1;
    state->cancelled.store(false, std::memory_order_release);
    const uint64_t generation =
	bobol_atomic_identity_advance(state->generation,
	    std::memory_order_acq_rel, std::memory_order_acquire);
    BObolRtRenderSettings settings;
    settings.width = endpoint->width;
    settings.height = endpoint->height;
    settings.workers = state->workers ? state->workers :
	std::max(1u, std::thread::hardware_concurrency());
    settings.samples = state->samples;
    settings.backgroundBottom = endpoint->controller->getBackgroundBottomColor();
    settings.backgroundTop = endpoint->controller->getBackgroundTopColor();
    state->worker = std::thread(endpoint_rt_worker, endpoint, state,
	generation, settings, state->previewScale, state->frameBudgetMilliseconds,
	state->interactiveQuality);
    return 1;
}

static void
endpoint_rt_destroy(bobol_display_endpoint_t *endpoint)
{
    if (!endpoint || !endpoint->rt)
	return;
    endpoint_rt_stop(endpoint);
    if (endpoint->controller)
	endpoint->controller->clearPresentationSyncCallback(endpoint);
    endpoint_rt_remove_presentation(endpoint);
    delete endpoint->rt;
    endpoint->rt = NULL;
}

#define FACEPLATE_AXES_STYLE_PROPERTIES(_axis) \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".position.x", \
	BV_DISPLAY_PROPERTY_DOUBLE, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, -1.0e12, 1.0e12, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".position.y", \
	BV_DISPLAY_PROPERTY_DOUBLE, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, -1.0e12, 1.0e12, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".position.z", \
	BV_DISPLAY_PROPERTY_DOUBLE, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, -1.0e12, 1.0e12, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".size", \
	BV_DISPLAY_PROPERTY_DOUBLE, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1.0e12, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".line_width", \
	BV_DISPLAY_PROPERTY_UINT, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 1.0, 1048576.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".color", \
	BV_DISPLAY_PROPERTY_COLOR3, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".position_only", \
	BV_DISPLAY_PROPERTY_BOOL, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".labels.visible", \
	BV_DISPLAY_PROPERTY_BOOL, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".labels.color", \
	BV_DISPLAY_PROPERTY_COLOR3, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".triple_color", \
	BV_DISPLAY_PROPERTY_BOOL, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".ticks.visible", \
	BV_DISPLAY_PROPERTY_BOOL, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".ticks.length", \
	BV_DISPLAY_PROPERTY_UINT, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 1.0, 1048576.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".ticks.major_length", \
	BV_DISPLAY_PROPERTY_UINT, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 1.0, 1048576.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".ticks.interval", \
	BV_DISPLAY_PROPERTY_DOUBLE, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1.0e12, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".ticks.per_major", \
	BV_DISPLAY_PROPERTY_UINT, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1048576.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".ticks.threshold", \
	BV_DISPLAY_PROPERTY_UINT, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 2147483647.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".ticks.color", \
	BV_DISPLAY_PROPERTY_COLOR3, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}, \
    {sizeof(bv_display_property_desc), "view.faceplate." _axis ".ticks.major_color", \
	BV_DISPLAY_PROPERTY_COLOR3, \
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE, \
	0, 0.0, 1.0, NULL}

static const bv_display_property_desc endpoint_properties[] = {
    {sizeof(bv_display_property_desc), "endpoint.renderer",
	BV_DISPLAY_PROPERTY_ENUM,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 0.0, "auto,hw,sw,rt,none,diagnostic"},
    {sizeof(bv_display_property_desc), "endpoint.renderer.resolved",
	BV_DISPLAY_PROPERTY_ENUM, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 0.0, "auto,hw,sw,rt,none,diagnostic"},
    {sizeof(bv_display_property_desc), "endpoint.renderer.supported",
	BV_DISPLAY_PROPERTY_STRING, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 0.0, NULL},
    {sizeof(bv_display_property_desc), "endpoint.width",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 1.0, 32767.0, NULL},
    {sizeof(bv_display_property_desc), "endpoint.height",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 1.0, 32767.0, NULL},
    {sizeof(bv_display_property_desc), "endpoint.device_pixel_ratio",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.01, 64.0, NULL},
    {sizeof(bv_display_property_desc), "endpoint.host",
	BV_DISPLAY_PROPERTY_STRING, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 0.0, NULL},
    {sizeof(bv_display_property_desc), "endpoint.title",
	BV_DISPLAY_PROPERTY_STRING,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	BOBOL_HOST_CAP_TOPLEVEL, 0.0, 0.0, NULL},
    {sizeof(bv_display_property_desc), "endpoint.visible",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	BOBOL_HOST_CAP_TOPLEVEL, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "endpoint.vsync",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	BOBOL_HOST_CAP_PRESENT_VSYNC, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "controller.background.bottom",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "controller.background.top",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "controller.software_wire",
	BV_DISPLAY_PROPERTY_ENUM,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 0.0, "auto,quality,fast"},
    {sizeof(bv_display_property_desc), "renderer.depth_test",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.lighting",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.lighting.profile",
	BV_DISPLAY_PROPERTY_ENUM,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 0.0, "studio,mged"},
    {sizeof(bv_display_property_desc), "renderer.headlight",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.headlight.tracking",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.scene_lights",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.headlight.color",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.headlight.intensity",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.headlight.direction",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, -1.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.transparency",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.antialiasing",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.diagnostic.revision",
	BV_DISPLAY_PROPERTY_UINT, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 0.0, NULL},
    {sizeof(bv_display_property_desc),
	"renderer.diagnostic.visited_sources",
	BV_DISPLAY_PROPERTY_UINT, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 0.0, NULL},
    {sizeof(bv_display_property_desc),
	"renderer.diagnostic.realized_sources",
	BV_DISPLAY_PROPERTY_UINT, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 0.0, NULL},
    {sizeof(bv_display_property_desc),
	"renderer.diagnostic.failed_sources",
	BV_DISPLAY_PROPERTY_UINT, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 0.0, NULL},
    {sizeof(bv_display_property_desc),
	"renderer.diagnostic.progressive_pending",
	BV_DISPLAY_PROPERTY_BOOL, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "renderer.diagnostic.summary",
	BV_DISPLAY_PROPERTY_STRING, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 0.0, NULL},
    {sizeof(bv_display_property_desc), "render.rt.workers",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 1.0, 64.0, NULL},
    {sizeof(bv_display_property_desc), "render.rt.samples",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 1.0, 64.0, NULL},
    {sizeof(bv_display_property_desc), "render.rt.preview_scale",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 1.0, 16.0, NULL},
    {sizeof(bv_display_property_desc), "render.rt.frame_budget_ms",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 1.0, 1000.0, NULL},
    {sizeof(bv_display_property_desc), "render.rt.quality",
	BV_DISPLAY_PROPERTY_ENUM,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 0.0, "interactive,final"},
    {sizeof(bv_display_property_desc), "render.rt.geometry_revision",
	BV_DISPLAY_PROPERTY_UINT, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 0.0, NULL},
    {sizeof(bv_display_property_desc), "render.rt.presentation_revision",
	BV_DISPLAY_PROPERTY_UINT, BV_DISPLAY_PROPERTY_READ,
	0, 0.0, 0.0, NULL},
    {sizeof(bv_display_property_desc), "composition.rt.layer",
	BV_DISPLAY_PROPERTY_ENUM,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 0.0, "underlay,interlay,overlay"},
    {sizeof(bv_display_property_desc), "renderer.clip.minimum",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, -1.0e12, 1.0e12, NULL},
    {sizeof(bv_display_property_desc), "renderer.clip.maximum",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, -1.0e12, 1.0e12, NULL},
    {sizeof(bv_display_property_desc), "renderer.depth_cue",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.perspective",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 179.999, NULL},
    {sizeof(bv_display_property_desc), "view.zclip",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.navigation.min_delta",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, -10000.0, 0.0, NULL},
    {sizeof(bv_display_property_desc), "view.navigation.max_delta",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 10000.0, NULL},
    {sizeof(bv_display_property_desc), "view.navigation.rotate_scale",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.001, 10.0, NULL},
    {sizeof(bv_display_property_desc), "view.navigation.scale_scale",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.001, 100.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.adc.visible",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.adc.line_color",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.adc.tick_color",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.adc.line_width",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 1.0, 1048576.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.center_dot.visible",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.visible",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.adaptive",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.snap",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.anchor.x",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, -1.0e12, 1.0e12, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.anchor.y",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, -1.0e12, 1.0e12, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.anchor.z",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, -1.0e12, 1.0e12, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.resolution.horizontal",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0e12, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.resolution.vertical",
	BV_DISPLAY_PROPERTY_DOUBLE,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0e12, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.major.horizontal",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 2147483647.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.major.vertical",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 2147483647.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.grid.color",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc),
	"view.interactive.rectangle.visible",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc),
	"view.interactive.rectangle.line_width",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1048576.0, NULL},
    {sizeof(bv_display_property_desc),
	"view.interactive.rectangle.line_style",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc),
	"view.interactive.rectangle.color",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.model_axes.visible",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    FACEPLATE_AXES_STYLE_PROPERTIES("model_axes"),
    {sizeof(bv_display_property_desc), "view.faceplate.scale.visible",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.view_axes.visible",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    FACEPLATE_AXES_STYLE_PROPERTIES("view_axes"),
    {sizeof(bv_display_property_desc), "view.faceplate.params.visible",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.params.size",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.params.center",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.params.azimuth",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.params.elevation",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.params.twist",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.params.fps",
	BV_DISPLAY_PROPERTY_BOOL,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.params.font_size",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 5.0, 96.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.center_dot.font_size",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 5.0, 96.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.scale.font_size",
	BV_DISPLAY_PROPERTY_UINT,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 5.0, 96.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.params.color",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.center_dot.color",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "view.faceplate.scale.color",
	BV_DISPLAY_PROPERTY_COLOR3,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 1.0, NULL},
    {sizeof(bv_display_property_desc), "composition.framebuffer.mode",
	BV_DISPLAY_PROPERTY_ENUM,
	BV_DISPLAY_PROPERTY_READ | BV_DISPLAY_PROPERTY_WRITE,
	0, 0.0, 0.0, "off,overlay,underlay,interlay"}
};

#undef FACEPLATE_AXES_STYLE_PROPERTIES

static const bv_display_property_desc *
endpoint_property(const char *name)
{
    if (!name || !name[0])
	return NULL;
    for (const bv_display_property_desc &property : endpoint_properties) {
	if (bu_strcmp(property.name, name) == 0)
	    return &property;
    }
    return NULL;
}

static int
endpoint_property_host_supported(const bobol_display_endpoint_t *endpoint,
	const bv_display_property_desc *property)
{
    if (!endpoint || !property)
	return 0;
    if (property->required_host_capabilities &&
	(bobol_display_endpoint_host_capabilities(endpoint) &
	 property->required_host_capabilities) !=
	 property->required_host_capabilities)
	return 0;
    if (bu_strncmp(property->name, "renderer.diagnostic.", 20) == 0 &&
	endpoint->engine != BOBOL_RENDER_ENGINE_DIAGNOSTIC)
	return 0;
    if ((bu_strcmp(property->name, "endpoint.title") == 0 ||
	 bu_strcmp(property->name, "endpoint.visible") == 0) &&
	endpoint->host_mode != BOBOL_HOST_MODE_TOPLEVEL)
	return 0;
    return 1;
}

static const char *
render_engine_name(enum bobol_render_engine engine)
{
    static const char *names[] = {
	"auto", "hw", "sw", "rt", "none", "diagnostic"
    };
    return valid_engine(engine) ? names[static_cast<int>(engine)] : "auto";
}

static int
render_engine_from_name(const char *name, enum bobol_render_engine *engine)
{
    if (!name || !engine)
	return 0;
    for (int i = BOBOL_RENDER_ENGINE_AUTO;
	    i <= BOBOL_RENDER_ENGINE_DIAGNOSTIC; i++) {
	if (bu_strcmp(name, render_engine_name(
		static_cast<enum bobol_render_engine>(i))) == 0) {
	    *engine = static_cast<enum bobol_render_engine>(i);
	    return 1;
	}
    }
    return 0;
}

static bool
valid_engine(enum bobol_render_engine engine)
{
    return engine >= BOBOL_RENDER_ENGINE_AUTO &&
	engine <= BOBOL_RENDER_ENGINE_DIAGNOSTIC;
}

static bool
engine_host_compatible(enum bobol_render_engine engine,
	uint64_t capabilities)
{
    switch (engine) {
	case BOBOL_RENDER_ENGINE_HW:
	    return (capabilities & BOBOL_HOST_CAP_SYSTEM_GL) != 0;
	case BOBOL_RENDER_ENGINE_SW:
	    return (capabilities & BOBOL_HOST_CAP_PIXEL_PRESENT) != 0;
	case BOBOL_RENDER_ENGINE_RT:
	    return (capabilities & BOBOL_HOST_CAP_PROGRESSIVE_PRESENT) != 0;
	default:
	    return true;
    }
}

static uint64_t
render_engine_capability(enum bobol_render_engine engine)
{
    switch (engine) {
	case BOBOL_RENDER_ENGINE_AUTO:
	    return BOBOL_RENDER_ENGINE_CAP_AUTO;
	case BOBOL_RENDER_ENGINE_HW:
	    return BOBOL_RENDER_ENGINE_CAP_HW;
	case BOBOL_RENDER_ENGINE_SW:
	    return BOBOL_RENDER_ENGINE_CAP_SW;
	case BOBOL_RENDER_ENGINE_RT:
	    return BOBOL_RENDER_ENGINE_CAP_RT;
	case BOBOL_RENDER_ENGINE_NONE:
	    return BOBOL_RENDER_ENGINE_CAP_NONE;
	case BOBOL_RENDER_ENGINE_DIAGNOSTIC:
	    return BOBOL_RENDER_ENGINE_CAP_DIAGNOSTIC;
	default:
	    return 0;
    }
}

extern "C" int
bobol_display_endpoint_render_engine_supported(
	const bobol_display_endpoint_t *endpoint,
	enum bobol_render_engine engine)
{
    if (!endpoint || !valid_engine(engine))
	return 0;
    if (endpoint->factory)
	return engine_host_compatible(engine,
	    bobol_host_factory_capabilities(endpoint->factory)) ? 1 : 0;
    if (endpoint->host && (engine == BOBOL_RENDER_ENGINE_HW ||
	engine == BOBOL_RENDER_ENGINE_SW))
	return 0;
    return 1;
}

extern "C" uint64_t
bobol_display_endpoint_render_engine_capabilities(
	const bobol_display_endpoint_t *endpoint)
{
    if (!endpoint)
	return 0;
    uint64_t capabilities = 0;
    for (int i = BOBOL_RENDER_ENGINE_AUTO;
	i <= BOBOL_RENDER_ENGINE_DIAGNOSTIC; ++i) {
	const enum bobol_render_engine engine =
	    static_cast<enum bobol_render_engine>(i);
	if (bobol_display_endpoint_render_engine_supported(endpoint, engine))
	    capabilities |= render_engine_capability(engine);
    }
    return capabilities;
}

static const char *
endpoint_supported_engines(const bobol_display_endpoint_t *endpoint)
{
    if (!endpoint)
	return "";
    endpoint->supported_engines.clear();
    for (int i = BOBOL_RENDER_ENGINE_AUTO;
	i <= BOBOL_RENDER_ENGINE_DIAGNOSTIC; ++i) {
	const enum bobol_render_engine engine =
	    static_cast<enum bobol_render_engine>(i);
	if (!bobol_display_endpoint_render_engine_supported(endpoint, engine))
	    continue;
	if (!endpoint->supported_engines.empty())
	    endpoint->supported_engines += ",";
	endpoint->supported_engines += render_engine_name(engine);
    }
    return endpoint->supported_engines.c_str();
}

static int
endpoint_diagnostic_refresh(bobol_display_endpoint_t *endpoint)
{
    if (!endpoint || !endpoint->controller ||
	endpoint->engine != BOBOL_RENDER_ENGINE_DIAGNOSTIC)
	return 0;

    BObolViewController *controller = endpoint->controller;
    (void)controller->realizePending();
    BObolProgressiveStatus progressive;
    (void)controller->advanceProgressiveWork(NULL, &progressive);

    SoNode *root = controller->getSceneRoot();
    if (!root)
	root = controller->getRenderSceneRoot();
    std::ostringstream summary;
    if (!root) {
	summary << "scene=missing bounds=unavailable";
    } else {
	SoGetBoundingBoxAction bounds(controller->getViewportRegion());
	bounds.apply(root);
	const SbBox3f box = bounds.getBoundingBox();
	summary << "scene=available";
	if (box.isEmpty()) {
	    summary << " bounds=empty";
	} else {
	    const SbVec3f minimum = box.getMin();
	    const SbVec3f maximum = box.getMax();
	    summary << " bounds=" << minimum[0] << "," << minimum[1] << ","
		<< minimum[2] << ":" << maximum[0] << "," << maximum[1]
		<< "," << maximum[2];
	}
    }
    endpoint->diagnostic_visited_sources =
	controller->getLastVisitedSourceCount();
    endpoint->diagnostic_realized_sources =
	controller->getLastRealizedSourceCount();
    endpoint->diagnostic_failed_sources =
	controller->getLastFailedSourceCount();
    endpoint->diagnostic_progressive_pending = progressive.hasMore != 0 ||
	controller->hasProgressiveWorkPending();
    summary << " visited_sources=" << endpoint->diagnostic_visited_sources
	<< " realized_sources=" << endpoint->diagnostic_realized_sources
	<< " failed_sources=" << endpoint->diagnostic_failed_sources
	<< " progressive=" <<
	(endpoint->diagnostic_progressive_pending ? "pending" : "idle");
    const SbString &detail = controller->getLastDiagnostics();
    if (detail.getLength() > 0)
	summary << " detail=" << detail.getString();
    endpoint->diagnostic_summary = summary.str();
    bobol_identity_advance(endpoint->diagnostic_revision);
    controller->clearRenderRequest();
    return 1;
}

extern "C" int
bobol_display_endpoint_diagnostic_refresh(
	bobol_display_endpoint_t *endpoint)
{
    return endpoint_diagnostic_refresh(endpoint);
}

extern "C" bobol_display_endpoint_t *
bobol_display_endpoint_create(void *controller, unsigned int flags)
{
    bobol_init(NULL);

    bobol_display_endpoint_t *endpoint =
	new (std::nothrow) bobol_display_endpoint_t;
    if (!endpoint)
	return NULL;

    endpoint->controller = static_cast<BObolViewController *>(controller);
    endpoint->owns_controller =
	(flags & BOBOL_ENDPOINT_OWN_CONTROLLER) != 0;
    if (!endpoint->controller) {
	endpoint->controller = new (std::nothrow)
	    BObolViewController(new SoBRLSceneGroup);
	endpoint->owns_controller = true;
    }

    if (!endpoint->controller) {
	delete endpoint;
	return NULL;
    }

    const SbVec2s viewport =
	endpoint->controller->getViewportRegion().getWindowSize();
    endpoint->width = viewport[0] > 0 ?
	static_cast<unsigned int>(viewport[0]) : 1;
    endpoint->height = viewport[1] > 0 ?
	static_cast<unsigned int>(viewport[1]) : 1;

    return endpoint;
}

extern "C" void
bobol_display_endpoint_host_detach(bobol_display_endpoint_t *endpoint)
{
    if (!endpoint)
	return;

    endpoint_rt_stop(endpoint);

	/* The factory callback is endpoint-owned.  Once it is cleared, a second
	 * detach must not touch a controller that GED may already have released. */
    if (endpoint->factory_frame_callback_bound && endpoint->controller) {
	endpoint->controller->clearFrameRequestCallback(endpoint);
    }
    endpoint->factory_frame_callback_bound = false;

    if (endpoint->factory) {
	bobol_host_factory_instance_destroy(endpoint->factory,
	    endpoint->factory_instance);
	/* Do not let a defective factory leave its provider installed after the
	 * instance (and possibly the provider itself) has been destroyed. */
	if (endpoint->controller)
	    endpoint->controller->setRenderContextManager(NULL);
	bobol_host_factory_release(endpoint->factory);
	endpoint->factory = NULL;
	endpoint->factory_instance = NULL;
    }

    endpoint->host_mode = BOBOL_HOST_MODE_HEADLESS;
    endpoint->title.clear();
    endpoint->visible = false;
    endpoint->vsync = true;

    if (endpoint->host) {
	BObolWindowHost *host = endpoint->host;
	const bool owns_host = endpoint->owns_host;
	endpoint->host = NULL;
	endpoint->owns_host = false;

	host->attachController(NULL, FALSE);
	if (owns_host)
	    delete host;
    }

    endpoint->retireWork();
}

extern "C" void
bobol_display_endpoint_destroy(bobol_display_endpoint_t *endpoint)
{
    if (!endpoint)
	return;

    endpoint_rt_destroy(endpoint);
    bobol_display_endpoint_host_detach(endpoint);
    /* A borrowed controller must not retain an endpoint-only NONE or
     * DIAGNOSTIC render suppression after its endpoint is gone. */
    endpoint->graphicalRenderingSet(true);
    if (endpoint->owns_controller)
	delete endpoint->controller;
    endpoint->controller = NULL;
    delete endpoint;
}

extern "C" void *
bobol_display_endpoint_controller(const bobol_display_endpoint_t *endpoint)
{
    return endpoint ? endpoint->controller : NULL;
}

extern "C" int
bobol_display_endpoint_view_sync(bobol_display_endpoint_t *endpoint,
	const void *view_ctx)
{
    if (!endpoint || !endpoint->controller || !view_ctx)
	return 0;
    const SbString pending_reason = endpoint->controller->getRenderReason();
    SbBool camera_changed = FALSE;
    if (!endpoint->controller->syncCameraFromViewContext(view_ctx, TRUE,
	&camera_changed))
	return 0;
    if (endpoint->engine == BOBOL_RENDER_ENGINE_DIAGNOSTIC)
	return endpoint_diagnostic_refresh(endpoint);
    if (endpoint->engine != BOBOL_RENDER_ENGINE_RT)
	return 1;

    /* MGED and TclCAD synchronize passive view state before every native
     * refresh.  An unchanged sync must not cancel and restart an in-flight RT
     * frame; only a newly changed camera or pre-existing scene request needs
     * a new generation. */
    return (!camera_changed && !endpoint_rt_scene_request(
	pending_reason.getString())) || endpoint_rt_start(endpoint) ? 1 : 0;
}

extern "C" int
bobol_display_endpoint_host_bind(bobol_display_endpoint_t *endpoint,
	void *host_ptr, unsigned int flags)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(host_ptr);
    if (!endpoint || !endpoint->controller || !host) {
	return 0;
    }
    if (endpoint->engine == BOBOL_RENDER_ENGINE_HW ||
	endpoint->engine == BOBOL_RENDER_ENGINE_SW) {
	return 0;
    }

    if (endpoint->host == host) {
	endpoint->owns_host = endpoint->owns_host ||
	    (flags & BOBOL_ENDPOINT_OWN_HOST) != 0;
	if (host->getController() != endpoint->controller)
	    host->attachController(endpoint->controller, FALSE);
	return 1;
    }

    bobol_display_endpoint_host_detach(endpoint);
    host->attachController(endpoint->controller, FALSE);
	endpoint->host = host;
	endpoint->owns_host = (flags & BOBOL_ENDPOINT_OWN_HOST) != 0;
	return endpoint->engine != BOBOL_RENDER_ENGINE_RT ||
	    endpoint_rt_start(endpoint);
}

extern "C" void *
bobol_display_endpoint_host(const bobol_display_endpoint_t *endpoint)
{
    if (!endpoint)
	return NULL;
    return endpoint->factory_instance ? endpoint->factory_instance :
	endpoint->host;
}

extern "C" void *
bobol_display_endpoint_framebuffer_window_host(
	const bobol_display_endpoint_t *endpoint)
{
    if (!endpoint)
	return NULL;
    if (endpoint->factory && endpoint->factory_instance)
	return bobol_host_factory_instance_framebuffer_window_host(
	    endpoint->factory, endpoint->factory_instance);
    return endpoint->host;
}

extern "C" int
bobol_display_endpoint_host_open(bobol_display_endpoint_t *endpoint,
	const char *factory_name, const struct bobol_host_desc *desc)
{
    if (!endpoint || !endpoint->controller)
	return 0;

    struct bobol_host_desc required_desc = {};
    if (desc) {
	const size_t copy_size = desc->struct_size < sizeof(required_desc) ?
	    desc->struct_size : sizeof(required_desc);
	std::memcpy(&required_desc, desc, copy_size);
	required_desc.struct_size = sizeof(required_desc);
    } else {
	required_desc.struct_size = sizeof(required_desc);
	required_desc.mode = BOBOL_HOST_MODE_HEADLESS;
	required_desc.width = endpoint->width;
	required_desc.height = endpoint->height;
	required_desc.device_pixel_ratio = endpoint->device_pixel_ratio;
	required_desc.visible = 0;
    }
    if (endpoint->engine == BOBOL_RENDER_ENGINE_HW)
	required_desc.required_capabilities |= BOBOL_HOST_CAP_SYSTEM_GL;
    else if (endpoint->engine == BOBOL_RENDER_ENGINE_SW)
	required_desc.required_capabilities |= BOBOL_HOST_CAP_PIXEL_PRESENT;
	else if (endpoint->engine == BOBOL_RENDER_ENGINE_RT)
	required_desc.required_capabilities |= BOBOL_HOST_CAP_PROGRESSIVE_PRESENT;
    if (required_desc.vsync != BOBOL_HOST_VSYNC_AUTO)
	required_desc.required_capabilities |= BOBOL_HOST_CAP_PRESENT_VSYNC;
    required_desc.input_dispatch = endpoint_input_dispatch;
    required_desc.input_dispatch_data = endpoint;

    bobol_host_factory_token_t *factory =
	bobol_host_factory_acquire(factory_name, &required_desc);
    if (!factory)
	return 0;
    if (!engine_host_compatible(endpoint->engine,
	    bobol_host_factory_capabilities(factory))) {
	bobol_host_factory_release(factory);
	return 0;
    }

    const uint64_t factoryCapabilities =
	bobol_host_factory_capabilities(factory);
    const bool providerRequired =
	endpoint->engine == BOBOL_RENDER_ENGINE_HW ||
	endpoint->engine == BOBOL_RENDER_ENGINE_SW ||
	endpoint->engine == BOBOL_RENDER_ENGINE_RT ||
	(endpoint->engine == BOBOL_RENDER_ENGINE_AUTO &&
	 (factoryCapabilities & (BOBOL_HOST_CAP_SYSTEM_GL |
	 BOBOL_HOST_CAP_PIXEL_PRESENT |
	 BOBOL_HOST_CAP_PROGRESSIVE_PRESENT)) != 0);

    SoDB::ContextManager *previousProvider =
	endpoint->controller->getRenderContextManager();
    void *instance = NULL;
	if (!bobol_host_factory_instance_create(factory, &required_desc,
	endpoint->controller, &instance)) {
	/* Provisional binding is allowed to change controller policy while the
	 * candidate opens.  A rejected candidate must not disturb the live host. */
	endpoint->controller->setRenderContextManager(previousProvider);
	bobol_host_factory_release(factory);
	return 0;
    }
    if (providerRequired &&
	!endpoint->controller->getRenderContextManager()) {
	bobol_host_factory_instance_destroy(factory, instance);
	endpoint->controller->setRenderContextManager(previousProvider);
	bobol_host_factory_release(factory);
	return 0;
    }

    /* Keep the old host authoritative until its transactional detach.  The
     * candidate is reasserted immediately after that detach below. */
    endpoint->controller->setRenderContextManager(previousProvider);
    bobol_display_endpoint_host_detach(endpoint);
    /* The old host's teardown clears its provider from the shared controller.
     * Reassert the new host after that transactional replacement step. */
    if (!bobol_host_factory_instance_bind_controller(factory, instance,
	    endpoint->controller) ||
	(providerRequired &&
	 !endpoint->controller->getRenderContextManager())) {
	bobol_host_factory_instance_destroy(factory, instance);
	endpoint->controller->setRenderContextManager(NULL);
	bobol_host_factory_release(factory);
	return 0;
    }
    endpoint->factory = factory;
    endpoint->factory_instance = instance;
    endpoint->host_mode = required_desc.mode;
    endpoint->width = required_desc.width ? required_desc.width :
	endpoint->width;
    endpoint->height = required_desc.height ? required_desc.height :
	endpoint->height;
    endpoint->device_pixel_ratio = required_desc.device_pixel_ratio > 0.0 ?
	required_desc.device_pixel_ratio : endpoint->device_pixel_ratio;
    endpoint->title = required_desc.title ? required_desc.title : "";
    endpoint->visible = required_desc.visible != 0;
    endpoint->vsync = required_desc.vsync != BOBOL_HOST_VSYNC_OFF;
	endpoint->controller->setFrameRequestCallback(endpoint_frame_requested,
	endpoint);
	endpoint->factory_frame_callback_bound = true;
	endpoint->resumeWork();
	if (endpoint->engine == BOBOL_RENDER_ENGINE_RT) {
	    if (!endpoint_rt_start(endpoint))
		return 0;
	} else if (endpoint->engine == BOBOL_RENDER_ENGINE_DIAGNOSTIC) {
	    (void)endpoint_diagnostic_refresh(endpoint);
	} else if (endpoint->engine == BOBOL_RENDER_ENGINE_NONE) {
	    endpoint->controller->clearRenderRequest();
	} else if (endpoint->controller->isRenderRequested() ||
	endpoint->controller->hasProgressiveWorkPending())
	endpoint_frame_requested(endpoint, "endpoint-host-open");
    return 1;
}

extern "C" const char *
bobol_display_endpoint_host_factory_name(
	const bobol_display_endpoint_t *endpoint)
{
    return endpoint && endpoint->factory ?
	bobol_host_factory_name(endpoint->factory) : NULL;
}

extern "C" uint64_t
bobol_display_endpoint_host_capabilities(
	const bobol_display_endpoint_t *endpoint)
{
    return endpoint && endpoint->factory ?
	bobol_host_factory_capabilities(endpoint->factory) : 0;
}

static int
endpoint_request_frame(bobol_display_endpoint_t *endpoint,
	const char *reason, bool capacity_relevant)
{
    if (!endpoint || !endpoint->controller)
	return 0;
    if (endpoint->engine == BOBOL_RENDER_ENGINE_NONE) {
	endpoint->controller->clearRenderRequest();
	return 1;
    }
    if (endpoint->engine == BOBOL_RENDER_ENGINE_DIAGNOSTIC)
	return endpoint_diagnostic_refresh(endpoint);
    const bool render_already_pending =
	endpoint->controller->isRenderRequested() != FALSE;
    const SbString pending_reason = endpoint->controller->getRenderReason();

    if (capacity_relevant)
	endpoint->controller->requestLodCapacityRender(reason);
    else
	endpoint->controller->requestPresentationRender(reason);

    if (capacity_relevant && endpoint->engine == BOBOL_RENDER_ENGINE_RT &&
	(endpoint_rt_scene_request(pending_reason.getString()) ||
	 endpoint_rt_scene_request(reason)) && !endpoint_rt_start(endpoint))
	return 0;
    if (!endpoint->factory)
	return 1;
    /* A new request wakes the bound host through the controller callback.
     * Explicitly replay only an already-standing level; doing both for a new
     * request queues duplicate paints and can expose a transient LoD state to
     * one of them. */
    if (!render_already_pending)
	return 1;
    return bobol_host_factory_instance_request_frame(endpoint->factory,
	endpoint->factory_instance, reason);
}

extern "C" int
bobol_display_endpoint_request_frame(bobol_display_endpoint_t *endpoint,
	const char *reason)
{
    return endpoint_request_frame(endpoint, reason, true);
}

extern "C" int
bobol_display_endpoint_request_presentation_frame(
	bobol_display_endpoint_t *endpoint, const char *reason)
{
    return endpoint_request_frame(endpoint, reason, false);
}

extern "C" int
bobol_display_endpoint_resize(bobol_display_endpoint_t *endpoint,
	unsigned int width, unsigned int height, double device_pixel_ratio)
{
    if (!endpoint || !endpoint->controller || !width || !height ||
	device_pixel_ratio <= 0.0)
	return 0;
    const bool devicePixelRatioChanged =
	std::fabs(endpoint->device_pixel_ratio - device_pixel_ratio) >
	    std::numeric_limits<double>::epsilon();
    if (endpoint->factory && !bobol_host_factory_instance_resize(
	endpoint->factory, endpoint->factory_instance, width, height,
	device_pixel_ratio))
	return 0;
    endpoint->controller->setViewportSize(width, height);
    endpoint->width = width;
    endpoint->height = height;
    endpoint->device_pixel_ratio = device_pixel_ratio;
	if (devicePixelRatioChanged)
	    endpoint->rendererPerformanceChanged();
	if (endpoint->engine == BOBOL_RENDER_ENGINE_NONE) {
	    endpoint->controller->clearRenderRequest();
	    return 1;
	}
	if (endpoint->engine == BOBOL_RENDER_ENGINE_DIAGNOSTIC)
	    return endpoint_diagnostic_refresh(endpoint);
	return endpoint->engine != BOBOL_RENDER_ENGINE_RT ||
	    endpoint_rt_start(endpoint);
}

extern "C" int
bobol_display_endpoint_input_profile_set(
	bobol_display_endpoint_t *endpoint,
	const BObolInputProfile *profile)
{
    if (!endpoint)
	return 0;
    endpoint->input.setProfile(profile ? profile :
	bobol_input_default_view_profile());
    return 1;
}

extern "C" int
bobol_display_endpoint_input_action_handler_set(
	bobol_display_endpoint_t *endpoint, BObolInputActionHandler handler,
	void *user_data)
{
    if (!endpoint)
	return 0;
    endpoint->input.setActionHandler(handler, user_data);
    return 1;
}

extern "C" int
bobol_display_endpoint_input_action_handler_clear_if(
	bobol_display_endpoint_t *endpoint, BObolInputActionHandler handler,
	void *user_data)
{
    return endpoint ? endpoint->input.clearActionHandlerIf(handler,
	user_data) : 0;
}

extern "C" int
bobol_display_endpoint_input_action_layer_set(
	bobol_display_endpoint_t *endpoint,
	const BObolInputActionLayer *layer, void *owner, void *user_data)
{
    return endpoint ? endpoint->input.setActionLayer(layer, owner, user_data) :
	0;
}

extern "C" int
bobol_display_endpoint_input_action_layer_clear_if(
	bobol_display_endpoint_t *endpoint, void *owner)
{
    return endpoint ? endpoint->input.clearActionLayerIf(owner) : 0;
}

extern "C" int
bobol_display_endpoint_input_dispatch(bobol_display_endpoint_t *endpoint,
	const BObolInputEvent *event)
{
    return endpoint ? endpoint->input.dispatch(event) : -1;
}

extern "C" int
bobol_display_endpoint_capture(bobol_display_endpoint_t *endpoint,
	unsigned char **pixels, size_t *size, unsigned int *width,
	unsigned int *height, unsigned int *components)
{
    return bobol_display_endpoint_capture_plane(endpoint,
	BOBOL_CAPTURE_COMPOSITE, pixels, size, width, height, components);
}

extern "C" int
bobol_display_endpoint_capture_plane(bobol_display_endpoint_t *endpoint,
	enum bobol_capture_plane plane, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components)
{
    if (pixels)
	*pixels = NULL;
    if (size)
	*size = 0;
    if (width)
	*width = 0;
    if (height)
	*height = 0;
    if (components)
	*components = 0;
    if (!endpoint || !endpoint->controller || !pixels || !size || !width ||
	!height || !components)
	return 0;

    if (plane == BOBOL_CAPTURE_FRAMEBUFFER) {
	return endpoint->framebuffer_capture_callback ?
	    endpoint->framebuffer_capture_callback(
		endpoint->framebuffer_capture_user_data, pixels, size, width,
		height, components) > 0 ? 1 : 0 : 0;
    }
    if (plane != BOBOL_CAPTURE_COMPOSITE)
	return 0;

    /* NONE retains state without rendering.  DIAGNOSTIC traverses and
     * reports; neither is a graphical backend or a source of composite
     * pixels.  The independent framebuffer plane remains capturable above. */
    if (endpoint->engine == BOBOL_RENDER_ENGINE_NONE ||
	endpoint->engine == BOBOL_RENDER_ENGINE_DIAGNOSTIC)
	return 0;

    if (endpoint->engine == BOBOL_RENDER_ENGINE_RT && endpoint->rt) {
	if (endpoint->rt->worker.joinable())
	    endpoint->rt->worker.join();
	endpoint->controller->synchronizePresentation();
    }

    if (endpoint->factory)
	return bobol_host_factory_instance_capture(endpoint->factory,
	    endpoint->factory_instance, pixels, size, width, height,
	    components) > 0 ? 1 : 0;

    if (endpoint->controller->renderToImage(pixels, 0, 0) != BRLCAD_OK ||
	!(*pixels))
	return 0;
    const SbVec2s viewport =
	endpoint->controller->getViewportRegion().getWindowSize();
    *width = viewport[0] > 0 ? (unsigned int)viewport[0] : 1;
    *height = viewport[1] > 0 ? (unsigned int)viewport[1] : 1;
    *components = 3;
    *size = (size_t)(*width) * (size_t)(*height) * (*components);
    return 1;
}

extern "C" int
bobol_display_endpoint_rt_plane_capture(bobol_display_endpoint_t *endpoint,
	enum bobol_rt_output_plane plane, void **samples, size_t *size,
	unsigned int *width, unsigned int *height)
{
    if (samples)
	*samples = NULL;
    if (size)
	*size = 0;
    if (width)
	*width = 0;
    if (height)
	*height = 0;
    if (!endpoint || endpoint->engine != BOBOL_RENDER_ENGINE_RT ||
	!endpoint->rt || !samples || !size || !width || !height)
	return 0;

    EndpointRtState *state = endpoint->rt;
    if (state->worker.joinable())
	state->worker.join();
    endpoint->controller->synchronizePresentation();

    std::lock_guard<std::mutex> lock(state->readyMutex);
    const void *source = NULL;
    size_t bytes = 0;
    if (plane == BOBOL_RT_OUTPUT_DEPTH) {
	if (state->presentedPlanes.depth.empty())
	    return 0;
	source = state->presentedPlanes.depth.data();
	bytes = state->presentedPlanes.depth.size() * sizeof(float);
    } else if (plane == BOBOL_RT_OUTPUT_SOURCE_ID) {
	if (state->presentedPlanes.sourceIdentity.empty())
	    return 0;
	source = state->presentedPlanes.sourceIdentity.data();
	bytes = state->presentedPlanes.sourceIdentity.size() * sizeof(uint32_t);
    } else {
	return 0;
    }
    if (!state->presentedWidth || !state->presentedHeight ||
	bytes / (plane == BOBOL_RT_OUTPUT_DEPTH ? sizeof(float) :
	    sizeof(uint32_t)) != static_cast<size_t>(state->presentedWidth) *
	static_cast<size_t>(state->presentedHeight))
	return 0;

    void *copy = bu_malloc(bytes, "retained rt output plane");
    std::memcpy(copy, source, bytes);
    *samples = copy;
    *size = bytes;
    *width = state->presentedWidth;
    *height = state->presentedHeight;
    return 1;
}

extern "C" int
bobol_display_endpoint_rt_source_identity_get(
	const bobol_display_endpoint_t *endpoint, uint32_t identifier,
	struct bobol_rt_source_identity *out)
{
    if (!endpoint || endpoint->engine != BOBOL_RENDER_ENGINE_RT ||
	!endpoint->rt || !identifier || !out ||
	out->struct_size < sizeof(*out) || out->instance_key || out->path)
	return 0;
    BObolRtSourceIdentity identity;
    if (!endpoint->rt->renderer.getSourceIdentity(identifier, identity))
	return 0;
    out->database = identity.database;
    out->instance_key = bu_strdup(identity.instanceKey.getString());
    out->path = bu_strdup(identity.path.getString());
    out->source_revision = identity.sourceRevision;
    return 1;
}

extern "C" void
bobol_display_endpoint_rt_source_identity_clear(
	struct bobol_rt_source_identity *identity)
{
    if (!identity)
	return;
    if (identity->instance_key)
	bu_free(identity->instance_key, "retained rt instance key");
    if (identity->path)
	bu_free(identity->path, "retained rt path");
    identity->database = NULL;
    identity->instance_key = NULL;
    identity->path = NULL;
    identity->source_revision = 0;
}

extern "C" int
bobol_display_endpoint_framebuffer_capture_provider_set(
	bobol_display_endpoint_t *endpoint,
	bobol_endpoint_framebuffer_capture_callback callback,
	void *user_data)
{
    if (!endpoint || (callback && !user_data))
	return 0;
    if (!callback && endpoint->framebuffer_capture_user_data != user_data)
	return 0;
    endpoint->framebuffer_capture_callback = callback;
    endpoint->framebuffer_capture_user_data = callback ? user_data : NULL;
    return 1;
}

extern "C" int
bobol_display_endpoint_render_engine_set(
	bobol_display_endpoint_t *endpoint,
	enum bobol_render_engine engine)
{
    if (!endpoint || !valid_engine(engine)) {
	return 0;
    }
    if (!bobol_display_endpoint_render_engine_supported(endpoint, engine))
	return 0;

    const enum bobol_render_engine previous = endpoint->engine;
    if (previous == engine) {
	if (engine == BOBOL_RENDER_ENGINE_DIAGNOSTIC)
	    return endpoint_diagnostic_refresh(endpoint);
	return 1;
    }
    if (previous == BOBOL_RENDER_ENGINE_RT)
	endpoint_rt_destroy(endpoint);

    endpoint->engine = engine;
    endpoint->rendererPerformanceChanged();
    endpoint->graphicalRenderingSet(
	engine != BOBOL_RENDER_ENGINE_NONE &&
	engine != BOBOL_RENDER_ENGINE_DIAGNOSTIC);
    if (engine == BOBOL_RENDER_ENGINE_RT) {
	if (!endpoint->rt)
	    endpoint->rt = new (std::nothrow) EndpointRtState;
	if (endpoint->rt) {
	    endpoint->controller->setPresentationSyncCallback(
		endpoint_rt_presentation_sync, endpoint);
	}
	if (!endpoint->rt || !endpoint_rt_start(endpoint)) {
	    endpoint_rt_destroy(endpoint);
	    endpoint->engine = previous;
	    endpoint->rendererPerformanceChanged();
	    endpoint->graphicalRenderingSet(
		previous != BOBOL_RENDER_ENGINE_NONE &&
		previous != BOBOL_RENDER_ENGINE_DIAGNOSTIC);
	    if (previous == BOBOL_RENDER_ENGINE_DIAGNOSTIC)
		(void)endpoint_diagnostic_refresh(endpoint);
	    return 0;
	}
    } else if (engine == BOBOL_RENDER_ENGINE_DIAGNOSTIC) {
	return endpoint_diagnostic_refresh(endpoint);
    } else if (engine == BOBOL_RENDER_ENGINE_NONE) {
	endpoint->controller->clearRenderRequest();
    } else {
	endpoint->controller->requestLodCapacityRender("render-engine");
	endpoint_frame_requested(endpoint, "render-engine");
    }
    return 1;
}

extern "C" enum bobol_render_engine
bobol_display_endpoint_render_engine_get(
	const bobol_display_endpoint_t *endpoint)
{
    return endpoint ? endpoint->engine : BOBOL_RENDER_ENGINE_AUTO;
}

extern "C" enum bobol_render_engine
bobol_display_endpoint_render_engine_resolved_get(
	const bobol_display_endpoint_t *endpoint)
{
    if (!endpoint || endpoint->engine != BOBOL_RENDER_ENGINE_AUTO)
	return endpoint ? endpoint->engine : BOBOL_RENDER_ENGINE_AUTO;
    if (!endpoint->factory)
	return BOBOL_RENDER_ENGINE_AUTO;

    const uint64_t capabilities =
	bobol_host_factory_capabilities(endpoint->factory);
    if (capabilities & BOBOL_HOST_CAP_SYSTEM_GL)
	return BOBOL_RENDER_ENGINE_HW;
    if (capabilities & BOBOL_HOST_CAP_PIXEL_PRESENT)
	return BOBOL_RENDER_ENGINE_SW;
    return BOBOL_RENDER_ENGINE_AUTO;
}

extern "C" size_t
bobol_display_endpoint_property_count(void)
{
    return sizeof(endpoint_properties) / sizeof(endpoint_properties[0]);
}

extern "C" int
bobol_display_endpoint_property_descriptor(size_t index,
	struct bv_display_property_desc *out)
{
    if (!out || out->struct_size < sizeof(*out) ||
	index >= bobol_display_endpoint_property_count())
	return BV_DISPLAY_PROPERTY_INVALID;
    *out = endpoint_properties[index];
    return BV_DISPLAY_PROPERTY_OK;
}

static void
property_value_prepare(bv_display_property_value *value,
	enum bv_display_property_type type)
{
    const uint32_t struct_size = value->struct_size;
    std::memset(value, 0, sizeof(*value));
    value->struct_size = struct_size;
    value->type = type;
}

extern "C" int
bobol_display_endpoint_property_get(
	const bobol_display_endpoint_t *endpoint, const char *name,
	struct bv_display_property_value *out)
{
    if (!endpoint || !endpoint->controller || !out ||
	out->struct_size < sizeof(*out))
	return BV_DISPLAY_PROPERTY_INVALID;
    const bv_display_property_desc *property = endpoint_property(name);
    if (!property)
	return BV_DISPLAY_PROPERTY_UNKNOWN;
    if (!endpoint_property_host_supported(endpoint, property))
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    property_value_prepare(out, property->type);

    if (bu_strcmp(name, "endpoint.width") == 0 ||
	bu_strcmp(name, "endpoint.height") == 0 ||
	bu_strcmp(name, "endpoint.device_pixel_ratio") == 0)
	endpoint_dimensions_refresh(
	    const_cast<bobol_display_endpoint_t *>(endpoint));

    if (bu_strcmp(name, "endpoint.renderer") == 0) {
	out->string_value = render_engine_name(endpoint->engine);
    } else if (bu_strcmp(name, "endpoint.renderer.resolved") == 0) {
	out->string_value = render_engine_name(
	    bobol_display_endpoint_render_engine_resolved_get(endpoint));
    } else if (bu_strcmp(name, "endpoint.renderer.supported") == 0) {
	out->string_value = endpoint_supported_engines(endpoint);
    } else if (bu_strcmp(name, "endpoint.width") == 0) {
	out->uint_value = endpoint->width;
    } else if (bu_strcmp(name, "endpoint.height") == 0) {
	out->uint_value = endpoint->height;
    } else if (bu_strcmp(name, "endpoint.device_pixel_ratio") == 0) {
	out->double_value = endpoint->device_pixel_ratio;
    } else if (bu_strcmp(name, "endpoint.host") == 0) {
	const char *factory_name =
	    bobol_display_endpoint_host_factory_name(endpoint);
	out->string_value = factory_name ? factory_name :
	    (endpoint->host ? "bound" : "");
    } else if (bu_strcmp(name, "endpoint.title") == 0) {
	out->string_value = endpoint->title.c_str();
    } else if (bu_strcmp(name, "endpoint.visible") == 0) {
	out->bool_value = endpoint->visible ? 1 : 0;
    } else if (bu_strcmp(name, "endpoint.vsync") == 0) {
	out->bool_value = endpoint->vsync ? 1 : 0;
    } else if (bu_strcmp(name, "controller.background.bottom") == 0 ||
	bu_strcmp(name, "controller.background.top") == 0) {
	const SbColor color = bu_strcmp(name,
	    "controller.background.bottom") == 0 ?
	    endpoint->controller->getBackgroundBottomColor() :
	    endpoint->controller->getBackgroundTopColor();
	out->color3[0] = color[0];
	out->color3[1] = color[1];
	out->color3[2] = color[2];
    } else if (bu_strcmp(name, "controller.software_wire") == 0) {
	switch (endpoint->controller->getSoftwareWireMode()) {
	    case BObolViewController::SOFTWARE_WIRE_QUALITY:
		out->string_value = "quality";
		break;
	    case BObolViewController::SOFTWARE_WIRE_FAST:
		out->string_value = "fast";
		break;
	    default:
		out->string_value = "auto";
		break;
	}
    } else if (bu_strcmp(name, "renderer.depth_test") == 0) {
	out->bool_value = endpoint->controller->isDepthTestEnabled() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.lighting") == 0) {
	out->bool_value = endpoint->controller->isLightingEnabled() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.lighting.profile") == 0) {
	out->string_value = endpoint->controller->getLightingProfile() ==
	    BObolViewController::LIGHTING_MGED ? "mged" : "studio";
    } else if (bu_strcmp(name, "renderer.headlight") == 0) {
	out->bool_value = endpoint->controller->isHeadlightEnabled() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.headlight.tracking") == 0) {
	out->bool_value =
	    endpoint->controller->isHeadlightCameraTracked() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.scene_lights") == 0) {
	out->bool_value = endpoint->controller->isSceneLightsEnabled() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.headlight.color") == 0) {
	const SbColor color = endpoint->controller->getHeadlightColor();
	out->color3[0] = color[0];
	out->color3[1] = color[1];
	out->color3[2] = color[2];
    } else if (bu_strcmp(name, "renderer.headlight.intensity") == 0) {
	out->double_value = endpoint->controller->getHeadlightIntensity();
    } else if (bu_strcmp(name, "renderer.headlight.direction") == 0) {
	const SbVec3f dir = endpoint->controller->getHeadlightOffset();
	out->color3[0] = dir[0];
	out->color3[1] = dir[1];
	out->color3[2] = dir[2];
    } else if (bu_strcmp(name, "renderer.transparency") == 0) {
	out->bool_value = endpoint->controller->isTransparencyEnabled() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.antialiasing") == 0) {
	out->bool_value = endpoint->controller->isAntialiasingEnabled() ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.diagnostic.revision") == 0) {
	out->uint_value = endpoint->diagnostic_revision;
    } else if (bu_strcmp(name,
	"renderer.diagnostic.visited_sources") == 0) {
	out->uint_value = endpoint->diagnostic_visited_sources;
    } else if (bu_strcmp(name,
	"renderer.diagnostic.realized_sources") == 0) {
	out->uint_value = endpoint->diagnostic_realized_sources;
    } else if (bu_strcmp(name,
	"renderer.diagnostic.failed_sources") == 0) {
	out->uint_value = endpoint->diagnostic_failed_sources;
    } else if (bu_strcmp(name,
	"renderer.diagnostic.progressive_pending") == 0) {
	out->bool_value = endpoint->diagnostic_progressive_pending ? 1 : 0;
    } else if (bu_strcmp(name, "renderer.diagnostic.summary") == 0) {
	out->string_value = endpoint->diagnostic_summary.c_str();
    } else if (bu_strcmp(name, "render.rt.workers") == 0) {
	out->uint_value = endpoint->rt && endpoint->rt->workers ?
	    endpoint->rt->workers : std::max(1u, std::thread::hardware_concurrency());
    } else if (bu_strcmp(name, "render.rt.samples") == 0) {
	out->uint_value = endpoint->rt ? endpoint->rt->samples : 1;
    } else if (bu_strcmp(name, "render.rt.preview_scale") == 0) {
	out->uint_value = endpoint->rt ? endpoint->rt->previewScale : 4;
	    } else if (bu_strcmp(name, "render.rt.frame_budget_ms") == 0) {
	out->uint_value = endpoint->rt ? endpoint->rt->frameBudgetMilliseconds : 33;
    } else if (bu_strcmp(name, "render.rt.quality") == 0) {
	out->string_value = endpoint->rt && !endpoint->rt->interactiveQuality ?
	    "final" : "interactive";
    } else if (bu_strcmp(name, "render.rt.geometry_revision") == 0) {
	out->uint_value = endpoint->rt ?
	    endpoint->rt->renderer.getGeometryRevision() : 0u;
    } else if (bu_strcmp(name, "render.rt.presentation_revision") == 0) {
	out->uint_value = endpoint->rt ?
	    endpoint->rt->renderer.getPresentationRevision() : 0u;
    } else if (bu_strcmp(name, "composition.rt.layer") == 0) {
	out->string_value = endpoint->rt ? endpoint_rt_presentation_layer_name(
	    endpoint->rt->presentationLayer) : "interlay";
    } else if (bu_strcmp(name, "renderer.clip.minimum") == 0 ||
	bu_strcmp(name, "renderer.clip.maximum") == 0) {
	double minimum = 0.0;
	double maximum = 0.0;
	endpoint->controller->getClipBounds(minimum, maximum);
	out->double_value = bu_strcmp(name, "renderer.clip.minimum") == 0 ?
	    minimum : maximum;
    } else if (bu_strcmp(name, "renderer.depth_cue") == 0) {
	out->bool_value = endpoint->controller->isDepthCueEnabled() ? 1 : 0;
    } else if (endpoint->property_get_callback) {
	return endpoint->property_get_callback(endpoint->property_user_data,
	    name, out);
    } else {
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    }
    return BV_DISPLAY_PROPERTY_OK;
}

static int
valid_color3(const double color[3])
{
    return color && std::isfinite(color[0]) && std::isfinite(color[1]) &&
	std::isfinite(color[2]) && color[0] >= 0.0 && color[0] <= 1.0 &&
	color[1] >= 0.0 && color[1] <= 1.0 && color[2] >= 0.0 &&
	color[2] <= 1.0;
}

extern "C" int
bobol_display_endpoint_property_set(
	bobol_display_endpoint_t *endpoint, const char *name,
	const struct bv_display_property_value *value)
{
    if (!endpoint || !endpoint->controller || !value ||
	value->struct_size < sizeof(*value))
	return BV_DISPLAY_PROPERTY_INVALID;
    const bv_display_property_desc *property = endpoint_property(name);
    if (!property)
	return BV_DISPLAY_PROPERTY_UNKNOWN;
    if (!(property->access & BV_DISPLAY_PROPERTY_WRITE))
	return BV_DISPLAY_PROPERTY_READ_ONLY;
    if (value->type != property->type)
	return BV_DISPLAY_PROPERTY_INVALID;
    if (property->type == BV_DISPLAY_PROPERTY_BOOL &&
	value->bool_value != 0 && value->bool_value != 1)
	return BV_DISPLAY_PROPERTY_INVALID;
    if (property->type == BV_DISPLAY_PROPERTY_DOUBLE &&
	(!std::isfinite(value->double_value) ||
	 (property->minimum < property->maximum &&
	  (value->double_value < property->minimum ||
	   value->double_value > property->maximum))))
	return BV_DISPLAY_PROPERTY_INVALID;
    if (property->type == BV_DISPLAY_PROPERTY_UINT &&
	property->minimum < property->maximum &&
	(value->uint_value < static_cast<uint64_t>(property->minimum) ||
	 value->uint_value > static_cast<uint64_t>(property->maximum)))
	return BV_DISPLAY_PROPERTY_INVALID;
    if (property->type == BV_DISPLAY_PROPERTY_COLOR3 &&
	!valid_color3(value->color3))
	return BV_DISPLAY_PROPERTY_INVALID;
    if (!endpoint_property_host_supported(endpoint, property))
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;

    if (bu_strcmp(name, "endpoint.renderer") == 0) {
	enum bobol_render_engine engine = BOBOL_RENDER_ENGINE_AUTO;
	if (!render_engine_from_name(value->string_value, &engine))
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (!bobol_display_endpoint_render_engine_set(endpoint, engine))
	    return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    } else if (bu_strcmp(name, "endpoint.width") == 0) {
	if (value->uint_value < 1 || value->uint_value > 32767 ||
	    !bobol_display_endpoint_resize(endpoint,
		static_cast<unsigned int>(value->uint_value), endpoint->height,
		endpoint->device_pixel_ratio))
	    return BV_DISPLAY_PROPERTY_INVALID;
    } else if (bu_strcmp(name, "endpoint.height") == 0) {
	if (value->uint_value < 1 || value->uint_value > 32767 ||
	    !bobol_display_endpoint_resize(endpoint, endpoint->width,
		static_cast<unsigned int>(value->uint_value),
		endpoint->device_pixel_ratio))
	    return BV_DISPLAY_PROPERTY_INVALID;
    } else if (bu_strcmp(name, "endpoint.device_pixel_ratio") == 0) {
	if (!std::isfinite(value->double_value) ||
	    value->double_value < 0.01 || value->double_value > 64.0 ||
	    !bobol_display_endpoint_resize(endpoint, endpoint->width,
		endpoint->height, value->double_value))
	    return BV_DISPLAY_PROPERTY_INVALID;
    } else if (bu_strcmp(name, "endpoint.title") == 0) {
	if (!value->string_value)
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (!bobol_host_factory_instance_set_title(endpoint->factory,
	    endpoint->factory_instance, value->string_value))
	    return BV_DISPLAY_PROPERTY_UNSUPPORTED;
	endpoint->title = value->string_value;
    } else if (bu_strcmp(name, "endpoint.visible") == 0) {
	if (!bobol_host_factory_instance_set_visible(endpoint->factory,
	    endpoint->factory_instance, value->bool_value))
	    return BV_DISPLAY_PROPERTY_UNSUPPORTED;
	endpoint->visible = value->bool_value != 0;
    } else if (bu_strcmp(name, "endpoint.vsync") == 0) {
	const bool changed = endpoint->vsync != (value->bool_value != 0);
	if (!bobol_host_factory_instance_set_vsync(endpoint->factory,
	    endpoint->factory_instance, value->bool_value))
	    return BV_DISPLAY_PROPERTY_UNSUPPORTED;
	endpoint->vsync = value->bool_value != 0;
	if (changed)
	    endpoint->rendererPerformanceChanged();
    } else if (bu_strcmp(name, "controller.background.bottom") == 0 ||
	bu_strcmp(name, "controller.background.top") == 0) {
	if (!valid_color3(value->color3))
	    return BV_DISPLAY_PROPERTY_INVALID;
	const SbColor updated(static_cast<float>(value->color3[0]),
	    static_cast<float>(value->color3[1]),
	    static_cast<float>(value->color3[2]));
	const SbColor bottom = bu_strcmp(name,
	    "controller.background.bottom") == 0 ? updated :
	    endpoint->controller->getBackgroundBottomColor();
	const SbColor top = bu_strcmp(name,
	    "controller.background.top") == 0 ? updated :
	    endpoint->controller->getBackgroundTopColor();
	endpoint->controller->setBackgroundColors(bottom, top);
    } else if (bu_strcmp(name, "controller.software_wire") == 0) {
	BObolViewController::SoftwareWireMode mode =
	    BObolViewController::SOFTWARE_WIRE_AUTO;
	if (!value->string_value)
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (bu_strcmp(value->string_value, "quality") == 0)
	    mode = BObolViewController::SOFTWARE_WIRE_QUALITY;
	else if (bu_strcmp(value->string_value, "fast") == 0)
	    mode = BObolViewController::SOFTWARE_WIRE_FAST;
	else if (bu_strcmp(value->string_value, "auto") != 0)
	    return BV_DISPLAY_PROPERTY_INVALID;
	endpoint->controller->setSoftwareWireMode(mode);
    } else if (bu_strcmp(name, "renderer.depth_test") == 0) {
	endpoint->controller->setDepthTestEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "renderer.lighting") == 0) {
	endpoint->controller->setLightingEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "renderer.lighting.profile") == 0) {
	if (!value->string_value)
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (bu_strcmp(value->string_value, "studio") == 0)
	    endpoint->controller->setLightingProfile(
		BObolViewController::LIGHTING_STUDIO);
	else if (bu_strcmp(value->string_value, "mged") == 0)
	    endpoint->controller->setLightingProfile(
		BObolViewController::LIGHTING_MGED);
	else
	    return BV_DISPLAY_PROPERTY_INVALID;
    } else if (bu_strcmp(name, "renderer.headlight") == 0) {
	endpoint->controller->setHeadlightEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "renderer.headlight.tracking") == 0) {
	endpoint->controller->setHeadlightCameraTracked(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "renderer.scene_lights") == 0) {
	endpoint->controller->setSceneLightsEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "renderer.headlight.color") == 0) {
	endpoint->controller->setHeadlightColor(SbColor(
	    static_cast<float>(value->color3[0]),
	    static_cast<float>(value->color3[1]),
	    static_cast<float>(value->color3[2])));
    } else if (bu_strcmp(name, "renderer.headlight.intensity") == 0) {
	endpoint->controller->setHeadlightIntensity(
	    static_cast<float>(value->double_value));
    } else if (bu_strcmp(name, "renderer.headlight.direction") == 0) {
	endpoint->controller->setHeadlightOffset(SbVec3f(
	    static_cast<float>(value->color3[0]),
	    static_cast<float>(value->color3[1]),
	    static_cast<float>(value->color3[2])));
    } else if (bu_strcmp(name, "renderer.transparency") == 0) {
	endpoint->controller->setTransparencyEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "renderer.antialiasing") == 0) {
	endpoint->controller->setAntialiasingEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (bu_strcmp(name, "render.rt.workers") == 0 ||
	bu_strcmp(name, "render.rt.samples") == 0 ||
	bu_strcmp(name, "render.rt.preview_scale") == 0 ||
	bu_strcmp(name, "render.rt.frame_budget_ms") == 0 ||
	bu_strcmp(name, "render.rt.quality") == 0 ||
	bu_strcmp(name, "composition.rt.layer") == 0) {
	if (!endpoint->rt)
	    endpoint->rt = new (std::nothrow) EndpointRtState;
	if (!endpoint->rt)
	    return BV_DISPLAY_PROPERTY_UNSUPPORTED;
	if (bu_strcmp(name, "render.rt.workers") == 0)
	    endpoint->rt->workers = static_cast<unsigned int>(value->uint_value);
	else if (bu_strcmp(name, "render.rt.samples") == 0)
	    endpoint->rt->samples = static_cast<unsigned int>(value->uint_value);
	else if (bu_strcmp(name, "render.rt.preview_scale") == 0)
	    endpoint->rt->previewScale = static_cast<unsigned int>(value->uint_value);
	else if (bu_strcmp(name, "render.rt.frame_budget_ms") == 0)
	    endpoint->rt->frameBudgetMilliseconds =
		static_cast<unsigned int>(value->uint_value);
	else if (bu_strcmp(name, "render.rt.quality") == 0) {
	    if (!value->string_value ||
		(bu_strcmp(value->string_value, "interactive") != 0 &&
		 bu_strcmp(value->string_value, "final") != 0))
		return BV_DISPLAY_PROPERTY_INVALID;
	    endpoint->rt->interactiveQuality =
		bu_strcmp(value->string_value, "interactive") == 0;
	} else {
	    int layer = SoBRLViewportImage::INTERLAY;
	    if (!endpoint_rt_presentation_layer_from_name(value->string_value,
		&layer))
		return BV_DISPLAY_PROPERTY_INVALID;
	    endpoint->rt->presentationLayer = layer;
	    if (!endpoint_rt_presentation_layer_apply(endpoint))
		return BV_DISPLAY_PROPERTY_UNSUPPORTED;
	}
    } else if (bu_strcmp(name, "renderer.clip.minimum") == 0 ||
	bu_strcmp(name, "renderer.clip.maximum") == 0) {
	if (!std::isfinite(value->double_value) ||
	    value->double_value < -1.0e12 || value->double_value > 1.0e12)
	    return BV_DISPLAY_PROPERTY_INVALID;
	double minimum = 0.0;
	double maximum = 0.0;
	endpoint->controller->getClipBounds(minimum, maximum);
	if (bu_strcmp(name, "renderer.clip.minimum") == 0)
	    minimum = value->double_value;
	else
	    maximum = value->double_value;
	if (!endpoint->controller->setClipBounds(minimum, maximum))
	    return BV_DISPLAY_PROPERTY_INVALID;
    } else if (bu_strcmp(name, "renderer.depth_cue") == 0) {
	endpoint->controller->setDepthCueEnabled(
	    value->bool_value ? TRUE : FALSE);
    } else if (endpoint->property_set_callback) {
	const int ret = endpoint->property_set_callback(
	    endpoint->property_user_data, name, value);
	if (ret == BV_DISPLAY_PROPERTY_OK) {
	    if (bu_strncmp(name, "renderer.", 9) == 0)
		endpoint->rendererPerformanceChanged();
	    endpoint->controller->requestLodCapacityRender("external property");
	}
	return ret;
    } else {
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    }
    if (endpoint->engine == BOBOL_RENDER_ENGINE_RT &&
	!endpoint_rt_start(endpoint))
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    return BV_DISPLAY_PROPERTY_OK;
}

extern "C" int
bobol_display_endpoint_property_provider_set(
	bobol_display_endpoint_t *endpoint,
	bv_display_property_get_callback get_callback,
	bv_display_property_set_callback set_callback,
	void *user_data)
{
    if (!endpoint || ((get_callback || set_callback) && !user_data))
	return 0;
    endpoint->property_get_callback = get_callback;
    endpoint->property_set_callback = set_callback;
    endpoint->property_user_data = user_data;
    return 1;
}

BObolDisplayEndpoint::BObolDisplayEndpoint(
	BObolViewController *controller, bool takeOwnership) :
    endpoint(bobol_display_endpoint_create(controller,
	takeOwnership ? BOBOL_ENDPOINT_OWN_CONTROLLER : 0))
{
}

BObolDisplayEndpoint::~BObolDisplayEndpoint(void)
{
    bobol_display_endpoint_destroy(this->endpoint);
    this->endpoint = NULL;
}

bool
BObolDisplayEndpoint::isValid(void) const
{
    return this->endpoint != NULL;
}

BObolViewController *
BObolDisplayEndpoint::controller(void) const
{
    return static_cast<BObolViewController *>(
	bobol_display_endpoint_controller(this->endpoint));
}

int
BObolDisplayEndpoint::bindHost(BObolWindowHost *host, bool takeOwnership)
{
    return bobol_display_endpoint_host_bind(this->endpoint, host,
	takeOwnership ? BOBOL_ENDPOINT_OWN_HOST : 0);
}

void
BObolDisplayEndpoint::detachHost(void)
{
    bobol_display_endpoint_host_detach(this->endpoint);
}

BObolWindowHost *
BObolDisplayEndpoint::host(void) const
{
    return static_cast<BObolWindowHost *>(
	bobol_display_endpoint_host(this->endpoint));
}

int
BObolDisplayEndpoint::setRenderEngine(enum bobol_render_engine engine)
{
    return bobol_display_endpoint_render_engine_set(this->endpoint, engine);
}

int
BObolDisplayEndpoint::supportsRenderEngine(
	enum bobol_render_engine engine) const
{
    return bobol_display_endpoint_render_engine_supported(this->endpoint,
	engine);
}

uint64_t
BObolDisplayEndpoint::renderEngineCapabilities(void) const
{
    return bobol_display_endpoint_render_engine_capabilities(this->endpoint);
}

int
BObolDisplayEndpoint::refreshDiagnostics(void)
{
    return bobol_display_endpoint_diagnostic_refresh(this->endpoint);
}

enum bobol_render_engine
BObolDisplayEndpoint::renderEngine(void) const
{
    return bobol_display_endpoint_render_engine_get(this->endpoint);
}

enum bobol_render_engine
BObolDisplayEndpoint::resolvedRenderEngine(void) const
{
    return bobol_display_endpoint_render_engine_resolved_get(this->endpoint);
}

int
BObolDisplayEndpoint::propertyGet(const char *name,
	struct bv_display_property_value *value) const
{
    return bobol_display_endpoint_property_get(this->endpoint, name, value);
}

int
BObolDisplayEndpoint::propertySet(const char *name,
	const struct bv_display_property_value *value)
{
    return bobol_display_endpoint_property_set(this->endpoint, name, value);
}

bobol_display_endpoint_t *
BObolDisplayEndpoint::release(void)
{
    bobol_display_endpoint_t *released = this->endpoint;
    this->endpoint = NULL;
    return released;
}

bobol_display_endpoint_t *
BObolDisplayEndpoint::get(void) const
{
    return this->endpoint;
}
