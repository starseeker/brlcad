/*      T E S T _ H E A D L E S S _ W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"

#include "bu/str.h"

#include "BObol.h"

#include "bu/log.h"
#include "bu/malloc.h"

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wshadow"
#endif
#include "imgstream/fb_compat.h"
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

#include <Inventor/SoDB.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoGroup.h>

#include <string.h>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)

class HeadlessTestContextManager : public SoDB::ContextManager
{
public:
    HeadlessTestContextManager(void) :
	renderCount(0),
	lastWidth(0),
	lastHeight(0),
	lastComponents(0),
	lastScene(NULL),
	lastBackground(0.0f, 0.0f, 0.0f)
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
	if (!scene || !pixels || width == 0 || height == 0 ||
	    nrcomponents < 1 || nrcomponents > 4)
	    return FALSE;

	this->renderCount++;
	this->lastWidth = width;
	this->lastHeight = height;
	this->lastComponents = nrcomponents;
	this->lastScene = scene;
	this->lastBackground.setValue(background_rgb);

	const size_t byteCount = (size_t)width * (size_t)height * (size_t)nrcomponents;
	for (size_t i = 0; i < byteCount; i += nrcomponents) {
	    pixels[i] = 31;
	    if (nrcomponents > 1)
		pixels[i + 1] = 63;
	    if (nrcomponents > 2)
		pixels[i + 2] = 127;
	    if (nrcomponents > 3)
		pixels[i + 3] = 255;
	}
	return TRUE;
    }

    int renderCount;
    unsigned int lastWidth;
    unsigned int lastHeight;
    unsigned int lastComponents;
    SoNode *lastScene;
    SbColor lastBackground;
};

struct HeadlessProgressiveState {
    int advances;
};

static int
headless_progressive_advance(BObolViewController *UNUSED(controller),
	void *userData,
	const BObolProgressiveOptions *UNUSED(options),
	BObolProgressiveStatus *UNUSED(status))
{
    HeadlessProgressiveState *state =
	static_cast<HeadlessProgressiveState *>(userData);
    state->advances++;
    return 1;
}

static BObolWindowDesc
visible_desc(unsigned int width, unsigned int height)
{
    BObolWindowDesc desc;
    desc.mode = BOBOL_WINDOW_TOPLEVEL;
    desc.backend = BOBOL_WINDOW_BACKEND_QT;
    desc.width = width;
    desc.height = height;
    desc.title = "headless-window-host-test";
    desc.display = ":99";
    desc.nativeIdHint = "native-test";
    desc.visible = TRUE;
    return desc;
}

static int
test_headless_contract(HeadlessTestContextManager *manager)
{
    BObolHeadlessWindowHost host;
    BObolWindowDesc desc = visible_desc(7, 5);

    host.setContextManager(NULL);
    CHECK(host.open(&desc) != 0,
	  "headless host rejects an ambiguous rendering provider");
    host.setContextManager(manager);
    CHECK(host.open(&desc) == 0, "headless host opens");
    CHECK(host.isOpen(), "headless host records open state");
    CHECK(host.getDesc().mode == BOBOL_WINDOW_HEADLESS,
	  "headless host forces headless window mode");
    CHECK(host.getDesc().backend == BOBOL_WINDOW_BACKEND_OFFSCREEN,
	  "headless host forces offscreen backend");
    CHECK(host.getDesc().visible == FALSE,
	  "headless host disables visible windows");
    CHECK(host.getController() != NULL, "headless host has a controller");
    CHECK(host.getController()->getCamera() != NULL,
	  "headless host seeds a default camera");
    CHECK(host.getController()->getViewportRegion().getWindowSize()[0] == 7 &&
	  host.getController()->getViewportRegion().getWindowSize()[1] == 5,
	  "headless host applies requested viewport size");

    host.setPollRate(42);
    CHECK(host.pollRate() == 42, "headless host stores positive poll rate");
    host.setPollRate(-1);
    CHECK(host.pollRate() == 0, "headless host clamps negative poll rate");
    host.setOutputComponents(4);
    CHECK(host.getOutputComponents() == 4, "headless host accepts RGBA output");
    host.setOutputComponents(99);
    CHECK(host.getOutputComponents() == 3, "headless host rejects invalid components");
    host.setBackgroundColor(SbColor(0.25f, 0.5f, 0.75f));
    CHECK(host.getController()->getBackgroundBottomColor() ==
	  SbColor(0.25f, 0.5f, 0.75f) &&
	  host.getController()->getBackgroundTopColor() ==
	  SbColor(0.25f, 0.5f, 0.75f),
	  "headless host background setter updates controller policy");
    return 0;
}

static int
test_headless_render_pending(HeadlessTestContextManager *manager)
{
    BObolHeadlessWindowHost host;
    host.setContextManager(manager);
    BObolWindowDesc desc = visible_desc(6, 4);
    CHECK(host.open(&desc) == 0, "headless host opens for render");

    HeadlessProgressiveState progressiveState = {0};
    const uint64_t progressiveToken =
	host.getController()->registerProgressiveProvider(
	    headless_progressive_advance, &progressiveState);
    CHECK(progressiveToken != 0, "headless host registers progressive provider");

    SoGroup *root = static_cast<SoGroup *>(host.getController()->getSceneRoot());
    root->addChild(new SoCube);
    host.getController()->setBackgroundColors(
	SbColor(0.125f, 0.25f, 0.5f), SbColor(0.75f, 0.5f, 0.25f));

    CHECK(host.renderPending() == 1, "open-requested render is drained");
    CHECK(!host.getController()->isRenderRequested(),
	  "successful headless render consumes request");
    CHECK(host.getRenderCount() == 1 && manager->renderCount == 1,
	  "headless host records render count");
    CHECK(progressiveState.advances == 1,
	  "headless host advances progressive work before rendering");
    CHECK(manager->lastBackground == SbColor(0.125f, 0.25f, 0.5f),
	  "headless render uses controller background policy");
    CHECK(bu_strcmp(host.getLastRenderReason().getString(), "background") == 0,
	  "headless host records the latest consumed render reason");
    CHECK(host.getLastFrameWidth() == 6 && host.getLastFrameHeight() == 4 &&
	  host.getLastFrameComponents() == 3,
	  "headless host records frame dimensions");
    CHECK(host.getLastFrameBufferSize() == 6 * 4 * 3,
	  "headless host records frame byte count");
    const unsigned char *pixels = host.getLastFrameBuffer();
    CHECK(pixels != NULL && pixels[0] == 31 && pixels[1] == 63 && pixels[2] == 127,
	  "headless host preserves rendered pixels");

    CHECK(host.poll(NULL) == 0, "headless poll is idle without pending render");
    host.getController()->requestRender("poll-render");
    CHECK(host.poll(NULL) == 1, "headless poll drains pending render");
    CHECK(bu_strcmp(host.getLastRenderReason().getString(), "poll-render") == 0,
	  "headless poll records render reason");
    CHECK(host.getController()->getSmoothedPresentationIntervalNanoseconds() > 0,
	  "headless host records completed-presentation cadence separately from render cost");
    CHECK(host.getController()->getDisplayedPresentationIntervalNanoseconds() > 0,
	  "headless host records a separately smoothed human-facing presentation cadence");
    host.getController()->unregisterProgressiveProvider(progressiveToken);
    return 0;
}

static int
test_imgstream_headless_poll(HeadlessTestContextManager *manager)
{
    BObolHeadlessWindowHost host;
    host.setContextManager(manager);
    host.setOutputComponents(4);

    imgstream_fb_t *fb = bobol_window_host_open_display_framebuffer(&host, "/dev/swrast", 3, 2);
    CHECK(fb != NULL, "swrast framebuffer opens through headless host");
    CHECK(host.isOpen(), "swrast framebuffer opens headless host");
    CHECK(host.getDesc().mode == BOBOL_WINDOW_HEADLESS &&
	  host.getDesc().backend == BOBOL_WINDOW_BACKEND_OFFSCREEN,
	  "swrast framebuffer uses headless offscreen descriptor");

    unsigned char rgb[3] = {255, 255, 255};
    CHECK(imgstream_fb_writerect(fb, 0, 0, 1, 1, rgb) == 1,
	  "swrast framebuffer write succeeds");
    CHECK(imgstream_fb_flush(fb) == 0,
	  "swrast framebuffer flush syncs Obol image nodes");
    CHECK(imgstream_fb_poll(fb) == 1,
	  "swrast framebuffer poll drains headless render");
    CHECK(host.getLastFrameWidth() == 3 && host.getLastFrameHeight() == 2 &&
	  host.getLastFrameComponents() == 4,
	  "swrast poll records RGBA headless frame");
    CHECK(bu_strcmp(host.getLastRenderReason().getString(), "fb-flush") == 0,
	  "swrast poll records framebuffer render reason");

    imgstream_fb_close(fb);
    return 0;
}

static int
test_controller_provider_switch(void)
{
    HeadlessTestContextManager first;
    HeadlessTestContextManager second;
    HeadlessTestContextManager oneShot;
    BObolHeadlessWindowHost host;
    BObolWindowDesc desc = visible_desc(4, 3);
    host.setContextManager(&first);
    CHECK(host.open(&desc) == 0, "provider-switch host opens");
    SoGroup *root = static_cast<SoGroup *>(host.getController()->getSceneRoot());
    root->addChild(new SoCube);

    unsigned char *image = NULL;
    CHECK(host.getController()->renderToImage(&image, 0, 0) == BRLCAD_OK &&
	image && first.renderCount == 1 && second.renderCount == 0,
	"implicit controller capture uses its explicitly bound provider");
    bu_free(image, "first explicit-provider image");

    host.setContextManager(&second);
    image = NULL;
    CHECK(host.getController()->renderToImage(&image, 0, 0) == BRLCAD_OK &&
	image && first.renderCount == 1 && second.renderCount == 1,
	"provider switch invalidates the cached offscreen renderer");
    bu_free(image, "second explicit-provider image");

    image = NULL;
    CHECK(host.getController()->renderToImage(&image, 0, 0, NULL,
	&oneShot) == BRLCAD_OK && image && oneShot.renderCount == 1 &&
	second.renderCount == 1,
	"per-call provider override is explicit and isolated");
    bu_free(image, "one-shot provider image");
    image = NULL;
    CHECK(host.getController()->renderToImage(&image, 0, 0) == BRLCAD_OK &&
	image && second.renderCount == 2 && oneShot.renderCount == 1,
	"one-shot override does not replace the controller provider");
    bu_free(image, "post-override provider image");

    host.close();
    image = NULL;
    CHECK(host.getController()->getRenderContextManager() == NULL &&
	host.getController()->renderToImage(&image, 0, 0) == BRLCAD_ERROR &&
	!image,
	"controller capture fails after its host releases the provider");
    return 0;
}

static int
test_headless_factory_endpoint(void)
{
    CHECK(bobol_headless_host_factory_register() != NULL,
	  "built-in headless factory is registered");

    struct bobol_host_desc desc;
    memset(&desc, 0, sizeof(desc));
    desc.struct_size = sizeof(desc);
    desc.mode = BOBOL_HOST_MODE_HEADLESS;
    desc.width = 5;
    desc.height = 3;
    desc.device_pixel_ratio = 1.0;
    desc.required_capabilities = BOBOL_HOST_CAP_READBACK;

    bobol_display_endpoint_t *endpoint =
	bobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL, "headless factory test creates endpoint");
    CHECK(bobol_display_endpoint_render_engine_set(endpoint,
	  BOBOL_RENDER_ENGINE_HW),
	  "endpoint accepts hardware policy before host selection");
    CHECK(!bobol_display_endpoint_host_open(endpoint, "headless", &desc),
	  "hardware renderer rejects a pixel-presentation-only host");
    CHECK(bobol_display_endpoint_render_engine_set(endpoint,
	  BOBOL_RENDER_ENGINE_SW),
	  "endpoint accepts software policy before compatible host selection");
    CHECK(bobol_display_endpoint_host_open(endpoint, "headless", &desc),
	  "endpoint opens the built-in headless factory");
    CHECK(bu_strcmp(bobol_display_endpoint_host_factory_name(endpoint),
	  "headless") == 0,
	  "endpoint reports built-in headless factory identity");
    CHECK(!bobol_display_endpoint_render_engine_set(endpoint,
	  BOBOL_RENDER_ENGINE_HW) &&
	  bobol_display_endpoint_render_engine_get(endpoint) ==
	  BOBOL_RENDER_ENGINE_SW,
	  "endpoint rejects an incompatible renderer change after host open");

    BObolViewController *controller =
	static_cast<BObolViewController *>(
	    bobol_display_endpoint_controller(endpoint));
    SoGroup *root = static_cast<SoGroup *>(controller->getSceneRoot());
    root->addChild(new SoCube);
    CHECK(bobol_display_endpoint_request_frame(endpoint, "factory-test"),
	  "headless factory accepts endpoint frame request");

    unsigned char *pixels = NULL;
    size_t size = 0;
    unsigned int width = 0;
    unsigned int height = 0;
    unsigned int components = 0;
    CHECK(bobol_display_endpoint_capture(endpoint, &pixels, &size, &width,
	  &height, &components),
	  "headless factory captures through endpoint API");
    CHECK(pixels && width == 5 && height == 3 && components == 3 &&
	  size == 5u * 3u * 3u,
	  "headless factory capture reports endpoint dimensions");
    bu_free(pixels, "headless factory endpoint capture");
    bobol_display_endpoint_destroy(endpoint);
    return 0;
}

int
main(int ac, char **av)
{
    bu_setprogname(av[0]);
    (void)ac;
    (void)av;

    HeadlessTestContextManager manager;
    bobol_init(&manager);
    CHECK(SoDB::getContextManager() == &manager,
	  "bobol_init installs test context manager");

    if (test_headless_contract(&manager))
	return 1;
    if (test_headless_render_pending(&manager))
	return 1;
    if (test_imgstream_headless_poll(&manager))
	return 1;
    if (test_controller_provider_switch())
	return 1;
    if (test_headless_factory_endpoint())
	return 1;

    return 0;
}
