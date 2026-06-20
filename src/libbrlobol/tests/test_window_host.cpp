/*              T E S T _ W I N D O W _ H O S T . C P P
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

#include <Inventor/SbViewportRegion.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoTexture2.h>

#include <math.h>
#include <string.h>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)

static int
float_equal(float a, float b)
{
    return fabs((double)a - (double)b) <= 1.0e-6;
}

static int
test_window_host_contract(void)
{
    BRLObolWindowHost host;
    BRLObolWindowDesc desc;

    desc.mode = BRLOBOL_WINDOW_HEADLESS;
    desc.backend = BRLOBOL_WINDOW_BACKEND_OFFSCREEN;
    desc.width = 6;
    desc.height = 4;
    desc.title = "window-host-test";
    desc.display = "";
    desc.nativeIdHint = "";
    desc.visible = FALSE;

    CHECK(host.open(&desc) == 0, "window host opens neutral controller");
    CHECK(host.isOpen(), "window host records open state");
    CHECK(host.getController() != NULL, "window host owns a view controller");
    CHECK(host.getController()->getSceneRoot() != NULL, "window host creates scene root");
    CHECK(host.getController()->getSceneRoot()->isOfType(SoGroup::getClassTypeId()),
	"window host scene root is a group");
    CHECK(host.getController()->getViewportRegion().getWindowSize()[0] == 6 &&
	    host.getController()->getViewportRegion().getWindowSize()[1] == 4,
	"window host applies requested viewport size");

    BRLObolInputBinding binding;
    binding.eventType = BRLOBOL_INPUT_KEY;
    binding.key = 'f';
    binding.button = 0;
    binding.modifiers = 0;
    binding.action = BRLOBOL_ACTION_VIEW_FRONT;

    BRLObolInputProfile profile;
    profile.name = "test";
    profile.bindings = &binding;
    profile.bindingCount = 1;

    BRLObolInputEvent event;
    memset(&event, 0, sizeof(event));
    event.type = BRLOBOL_INPUT_KEY;
    event.key = 'f';
    CHECK(host.handleInputEvent(&event, &profile) == 1,
	"window host applies semantic input profile action");
    CHECK(host.getController()->isRenderRequested(),
	"semantic input action requests render");
    host.getController()->clearRenderRequest();

    host.close();
    CHECK(!host.isOpen(), "window host closes");
    return 0;
}

static int
texture_matches(SoTexture2 *texture, int width, int height, int channels)
{
    if (!texture)
	return 0;
    int texWidth = 0;
    int texHeight = 0;
    int texChannels = 0;
    const unsigned char *pixels = texture->getImageData(texWidth, texHeight, texChannels);
    return pixels && texWidth == width && texHeight == height && texChannels == channels;
}

static int
test_imgstream_display_host_bridge(void)
{
    BRLObolWindowHost host;

    brlobol_window_host_unregister_display_host();
    CHECK(brlobol_window_host_register_display_host(NULL) == -1,
	"null Obol display host registration is rejected");
    CHECK(brlobol_window_host_register_display_host(&host) == 0,
	"Obol display host registered with libimgstream");

    imgstream_fb_t *fb = imgstream_fb_open("/dev/qtgl", 5, 4);
    CHECK(fb != NULL, "display framebuffer opens through Obol host");
    CHECK(host.isOpen(), "display framebuffer opens host");
    CHECK(host.getFramebufferCount() == 1, "host records one framebuffer");

    SoBRLImageSource *source = host.getFramebufferImageSource(fb);
    SoBRLViewportImage *viewport = host.getFramebufferViewportImage(fb);
    CHECK(source != NULL, "display host creates image source");
    CHECK(viewport != NULL, "display host creates viewport image");
    CHECK(source->getStream() == imgstream_fb_stream(fb),
	"display host image source borrows framebuffer stream");
    CHECK(viewport->getImageSource() == source,
	"display host viewport image references source");
    CHECK(viewport->getTextureNode() != NULL &&
	    texture_matches(viewport->getTextureNode(), 5, 4, 3),
	"display host realizes framebuffer texture");
    CHECK(host.getController()->getSceneRoot()->isOfType(SoGroup::getClassTypeId()),
	"display host controller root is a group");
    CHECK(static_cast<SoGroup *>(host.getController()->getSceneRoot())->getNumChildren() >= 2,
	"display host attaches image nodes to controller root");

    unsigned char red[3] = {255, 0, 0};
    CHECK(imgstream_fb_writerect(fb, 2, 1, 1, 1, red) == 1,
	"framebuffer write updates stream");
    CHECK(imgstream_fb_flush(fb) == 0, "display host flush syncs stream");
    CHECK(viewport->realizedDirtyRevision.getValue() == source->dirtyRevision.getValue(),
	"display host flush refreshes viewport dirty revision");

    CHECK(imgstream_fb_view(fb, 3, 2, 4, 4) == 0,
	"display host accepts framebuffer view transform");
    CHECK(float_equal(viewport->sourceCenter.getValue()[0], 3.0f) &&
	    float_equal(viewport->sourceCenter.getValue()[1], 2.0f) &&
	    float_equal(viewport->sourceZoom.getValue(), 4.0f),
	"display host maps framebuffer view to viewport image");

    CHECK(imgstream_fb_cursor(fb, 1, 4, 3) == 0,
	"display host accepts cursor state");
    CHECK(viewport->cursorVisible.getValue() == TRUE &&
	    float_equal(viewport->cursorImagePosition.getValue()[0], 4.0f) &&
	    float_equal(viewport->cursorImagePosition.getValue()[1], 3.0f),
	"display host maps cursor state to viewport image");

    CHECK(imgstream_fb_viewport(fb, 1, 2, 11, 8) == 0,
	"display host accepts viewport state");
    CHECK(host.getController()->getViewportRegion().getWindowSize()[0] == 10 &&
	    host.getController()->getViewportRegion().getWindowSize()[1] == 6,
	"display host maps viewport state to controller size");

    unsigned char cursorBits[1] = {0xff};
    CHECK(imgstream_fb_setcursor(fb, cursorBits, 8, 1, 0, 0) == 0,
	"display host accepts custom cursor shape");
    CHECK(viewport->cursorShape.getValue() == SoBRLViewportImage::CURSOR_CUSTOM,
	"display host records custom cursor policy");

    CHECK(imgstream_fb_poll(fb) == 0, "display host poll succeeds");
    CHECK(imgstream_fb_poll_rate(fb) == 0, "display host poll rate defaults to zero");

    imgstream_fb_close(fb);
    CHECK(host.getFramebufferCount() == 0, "display host closes framebuffer attachment");

    brlobol_window_host_unregister_display_host();
    CHECK(imgstream_fb_open("/dev/qtgl", 2, 2) == NULL,
	"display framebuffer fails explicitly without registered host");
    return 0;
}

int
main(int ac, char **av)
{
    (void)ac;
    (void)av;

    brlobol_init(NULL);

    if (test_window_host_contract())
	return 1;
    if (test_imgstream_display_host_bridge())
	return 1;

    return 0;
}
