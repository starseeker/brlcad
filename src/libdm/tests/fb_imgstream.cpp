/*                 F B _ I M G S T R E A M . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libdm/tests/fb_imgstream.cpp */

#include "common.h"

#include "brlobol.h"
#include "bu/app.h"
#include "bu/log.h"
#include "dm.h"

#include <Inventor/SoDB.h>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)

class HeadlessTestContextManager : public SoDB::ContextManager
{
public:
    HeadlessTestContextManager(void) : renderCount(0) {}

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
	const size_t byteCount = (size_t)width * (size_t)height *
	    (size_t)nrcomponents;
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
};

static int
test_legacy_swrast_fb_open_headless_host(void)
{
    if (dm_valid_type("swrast", NULL) != 1) {
	bu_log("SKIP: swrast display manager plugin is unavailable\n");
	return 0;
    }

    BRLObolHeadlessWindowHost host;
    host.setOutputComponents(4);

    brlobol_window_host_unregister_display_host();
    CHECK(brlobol_window_host_register_display_host(&host) == 0,
	"headless host registers for legacy swrast fb_open");

    struct fb *fbp = fb_open("/dev/swrast", 4, 3);
    CHECK(fbp != NULL, "legacy swrast fb_open uses registered Obol display host");
    CHECK(fb_get_dm(fbp) == NULL,
	"legacy swrast fb_open host bridge does not create nested libdm backend");
    CHECK(host.isOpen(), "legacy swrast fb_open opens headless host");
    CHECK(host.getFramebufferCount() == 1,
	"legacy swrast fb_open attaches one Obol framebuffer");

    unsigned char rgb[3] = {128, 64, 32};
    CHECK(fb_writerect(fbp, 0, 0, 1, 1, rgb) == 1,
	"legacy swrast fb_open bridge accepts RGB writes");
    unsigned char readback[3] = {0, 0, 0};
    CHECK(fb_read(fbp, 0, 0, readback, 1) == 1 &&
	    readback[0] == 128 && readback[1] == 64 && readback[2] == 32,
	"legacy swrast fb_open bridge reads from imgstream storage");
    CHECK(fb_flush(fbp) == 0,
	"legacy swrast fb_open bridge flushes Obol image nodes");
    CHECK(fb_poll(fbp) == 1,
	"legacy swrast fb_open bridge drains headless render");
    CHECK(host.getLastFrameWidth() == 4 && host.getLastFrameHeight() == 3 &&
	    host.getLastFrameComponents() == 4,
	"legacy swrast fb_open poll records RGBA headless frame");

    CHECK(fb_close(fbp) == 0,
	"legacy swrast fb_open bridge closes cleanly");
    CHECK(host.getFramebufferCount() == 0,
	"legacy swrast fb_open bridge detaches Obol framebuffer on close");
    brlobol_window_host_unregister_display_host();
    return 0;
}

int
main(int ac, char **av)
{
    (void)ac;
    bu_setprogname(av[0]);

    HeadlessTestContextManager manager;
    brlobol_init(&manager);
    CHECK(SoDB::getContextManager() == &manager,
	"brlobol_init installs test context manager");

    if (test_legacy_swrast_fb_open_headless_host())
	return 1;

    return 0;
}
