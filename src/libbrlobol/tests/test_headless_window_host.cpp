/*      T E S T _ H E A D L E S S _ W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol.h"

#include "bu/log.h"

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
	lastScene(NULL)
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

	this->renderCount++;
	this->lastWidth = width;
	this->lastHeight = height;
	this->lastComponents = nrcomponents;
	this->lastScene = scene;

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
};

static BRLObolWindowDesc
visible_desc(unsigned int width, unsigned int height)
{
    BRLObolWindowDesc desc;
    desc.mode = BRLOBOL_WINDOW_TOPLEVEL;
    desc.backend = BRLOBOL_WINDOW_BACKEND_QT;
    desc.width = width;
    desc.height = height;
    desc.title = "headless-window-host-test";
    desc.display = ":99";
    desc.nativeIdHint = "native-test";
    desc.visible = TRUE;
    return desc;
}

static int
test_headless_contract(void)
{
    BRLObolHeadlessWindowHost host;
    BRLObolWindowDesc desc = visible_desc(7, 5);

    CHECK(host.open(&desc) == 0, "headless host opens");
    CHECK(host.isOpen(), "headless host records open state");
    CHECK(host.getDesc().mode == BRLOBOL_WINDOW_HEADLESS,
	  "headless host forces headless window mode");
    CHECK(host.getDesc().backend == BRLOBOL_WINDOW_BACKEND_OFFSCREEN,
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
    return 0;
}

static int
test_headless_render_pending(HeadlessTestContextManager *manager)
{
    BRLObolHeadlessWindowHost host;
    BRLObolWindowDesc desc = visible_desc(6, 4);
    CHECK(host.open(&desc) == 0, "headless host opens for render");

    SoGroup *root = static_cast<SoGroup *>(host.getController()->getSceneRoot());
    root->addChild(new SoCube);

    CHECK(host.renderPending() == 1, "open-requested render is drained");
    CHECK(!host.getController()->isRenderRequested(),
	  "successful headless render consumes request");
    CHECK(host.getRenderCount() == 1 && manager->renderCount == 1,
	  "headless host records render count");
    CHECK(strcmp(host.getLastRenderReason().getString(), "camera") == 0,
	  "headless host records consumed render reason");
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
    CHECK(strcmp(host.getLastRenderReason().getString(), "poll-render") == 0,
	  "headless poll records render reason");
    return 0;
}

static int
test_imgstream_headless_poll(HeadlessTestContextManager *UNUSED(manager))
{
    BRLObolHeadlessWindowHost host;
    host.setOutputComponents(4);

    brlobol_window_host_unregister_display_host();
    CHECK(brlobol_window_host_register_display_host(&host) == 0,
	  "headless host registers as display host");

    imgstream_fb_t *fb = imgstream_fb_open("/dev/swrast", 3, 2);
    CHECK(fb != NULL, "swrast framebuffer opens through headless host");
    CHECK(host.isOpen(), "swrast framebuffer opens headless host");
    CHECK(host.getDesc().mode == BRLOBOL_WINDOW_HEADLESS &&
	  host.getDesc().backend == BRLOBOL_WINDOW_BACKEND_OFFSCREEN,
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
    CHECK(strcmp(host.getLastRenderReason().getString(), "fb-flush") == 0,
	  "swrast poll records framebuffer render reason");

    imgstream_fb_close(fb);
    brlobol_window_host_unregister_display_host();
    return 0;
}

int
main(int ac, char **av)
{
    (void)ac;
    (void)av;

    HeadlessTestContextManager manager;
    brlobol_init(&manager);
    CHECK(SoDB::getContextManager() == &manager,
	  "brlobol_init installs test context manager");

    if (test_headless_contract())
	return 1;
    if (test_headless_render_pending(&manager))
	return 1;
    if (test_imgstream_headless_poll(&manager))
	return 1;

    return 0;
}
