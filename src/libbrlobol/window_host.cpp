/*                W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/image_source.h"
#include "brlobol/view_controller.h"
#include "brlobol/viewport_image.h"
#include "brlobol/window_host.h"

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wshadow"
#endif
#include "imgstream/fb_compat.h"
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoSeparator.h>

#include <algorithm>
#include <string.h>
#include <vector>

struct BRLObolFramebufferAttachment {
    imgstream_fb_t *fb;
    SoBRLImageSource *source;
    SoBRLViewportImage *viewport;
};

struct BRLObolWindowHostPrivate {
    BRLObolWindowHostPrivate(void) :
	controller(new BRLObolViewController),
	ownsController(TRUE),
	open(FALSE),
	pollRate(0)
    {
	desc.mode = BRLOBOL_WINDOW_HEADLESS;
	desc.backend = BRLOBOL_WINDOW_BACKEND_OFFSCREEN;
	desc.width = 1;
	desc.height = 1;
	desc.title = "BRL-CAD Obol";
	desc.display = "";
	desc.nativeIdHint = "";
	desc.visible = FALSE;
    }

    BRLObolViewController *controller;
    SbBool ownsController;
    SbBool open;
    BRLObolWindowDesc desc;
    long pollRate;
    std::vector<BRLObolFramebufferAttachment> framebuffers;
};

static BRLObolWindowDesc
default_desc(void)
{
    BRLObolWindowDesc desc;
    desc.mode = BRLOBOL_WINDOW_HEADLESS;
    desc.backend = BRLOBOL_WINDOW_BACKEND_OFFSCREEN;
    desc.width = 1;
    desc.height = 1;
    desc.title = "BRL-CAD Obol";
    desc.display = "";
    desc.nativeIdHint = "";
    desc.visible = FALSE;
    return desc;
}

static void
sanitize_desc(BRLObolWindowDesc *desc)
{
    if (!desc)
	return;
    if (desc->width == 0)
	desc->width = 1;
    if (desc->height == 0)
	desc->height = 1;
}

static BRLObolWindowBackend
backend_from_fb_display(enum imgstream_fb_display_kind display)
{
    switch (display) {
	case IMGSTREAM_FB_DISPLAY_QTGL:
	    return BRLOBOL_WINDOW_BACKEND_QT;
	case IMGSTREAM_FB_DISPLAY_SWRAST:
	    return BRLOBOL_WINDOW_BACKEND_OFFSCREEN;
	case IMGSTREAM_FB_DISPLAY_X:
	    return BRLOBOL_WINDOW_BACKEND_TK;
	case IMGSTREAM_FB_DISPLAY_OGL:
	case IMGSTREAM_FB_DISPLAY_WGL:
	    return BRLOBOL_WINDOW_BACKEND_OPENGL;
	case IMGSTREAM_FB_DISPLAY_NONE:
	default:
	    return BRLOBOL_WINDOW_BACKEND_AUTO;
    }
}

static SoGroup *
controller_root_group(BRLObolViewController *controller)
{
    if (!controller)
	return NULL;

    SoNode *root = controller->getSceneRoot();
    if (!root) {
	SoSeparator *newRoot = new SoSeparator;
	newRoot->ref();
	controller->setSceneRoot(newRoot);
	newRoot->unref();
	root = controller->getSceneRoot();
    }

    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return NULL;
    return static_cast<SoGroup *>(root);
}

static void
remove_child_if_present(SoGroup *group, SoNode *node)
{
    if (!group || !node)
	return;
    for (int i = 0; i < group->getNumChildren(); i++) {
	if (group->getChild(i) == node) {
	    group->removeChild(i);
	    return;
	}
    }
}

static float
positive_zoom(int xzoom, int yzoom)
{
    int x = xzoom > 0 ? xzoom : 1;
    int y = yzoom > 0 ? yzoom : 1;
    return (float)((x < y) ? x : y);
}

static BRLObolFramebufferAttachment *
find_attachment(std::vector<BRLObolFramebufferAttachment> &attachments,
	imgstream_fb_t *fb)
{
    for (size_t i = 0; i < attachments.size(); i++) {
	if (attachments[i].fb == fb)
	    return &attachments[i];
    }
    return NULL;
}

static const BRLObolFramebufferAttachment *
find_attachment_const(const std::vector<BRLObolFramebufferAttachment> &attachments,
	imgstream_fb_t *fb)
{
    for (size_t i = 0; i < attachments.size(); i++) {
	if (attachments[i].fb == fb)
	    return &attachments[i];
    }
    return NULL;
}

BRLObolWindowHost::BRLObolWindowHost(void) :
    p(new BRLObolWindowHostPrivate)
{
}

BRLObolWindowHost::~BRLObolWindowHost(void)
{
    this->close();
    if (this->p->ownsController)
	delete this->p->controller;
    this->p->controller = NULL;
    delete this->p;
}

int
BRLObolWindowHost::open(const BRLObolWindowDesc *desc)
{
    this->p->desc = desc ? *desc : default_desc();
    sanitize_desc(&this->p->desc);

    if (!this->p->controller)
	return -1;

    if (!controller_root_group(this->p->controller))
	return -1;

    this->p->controller->setViewportSize(this->p->desc.width, this->p->desc.height);
    this->p->open = TRUE;
    return 0;
}

void
BRLObolWindowHost::close(void)
{
    while (!this->p->framebuffers.empty())
	this->closeFramebuffer(this->p->framebuffers.back().fb);
    this->p->open = FALSE;
}

SbBool
BRLObolWindowHost::isOpen(void) const
{
    return this->p->open;
}

const BRLObolWindowDesc &
BRLObolWindowHost::getDesc(void) const
{
    return this->p->desc;
}

void
BRLObolWindowHost::attachController(BRLObolViewController *controller,
	SbBool takeOwnership)
{
    this->close();
    if (this->p->ownsController)
	delete this->p->controller;
    this->p->controller = controller;
    this->p->ownsController = takeOwnership;
}

BRLObolViewController *
BRLObolWindowHost::getController(void) const
{
    return this->p->controller;
}

int
BRLObolWindowHost::handleInputEvent(const BRLObolInputEvent *event,
	const BRLObolInputProfile *profile)
{
    if (!event)
	return -1;
    if (!profile || !profile->bindings || profile->bindingCount == 0)
	return 0;

    for (size_t i = 0; i < profile->bindingCount; i++) {
	const BRLObolInputBinding &binding = profile->bindings[i];
	if (binding.eventType != event->type)
	    continue;
	if (binding.key && binding.key != event->key)
	    continue;
	if (binding.button && binding.button != event->button)
	    continue;
	if (binding.modifiers != event->modifiers)
	    continue;
	return this->applyInputAction(binding.action, event);
    }
    return 0;
}

int
BRLObolWindowHost::applyInputAction(BRLObolInputAction action,
	const BRLObolInputEvent *UNUSED(event))
{
    if (action == BRLOBOL_ACTION_NONE)
	return 0;
    if (!this->p->controller)
	return -1;

    this->p->controller->requestRender("input-action");
    return 1;
}

int
BRLObolWindowHost::poll(const BRLObolInputProfile *UNUSED(profile))
{
    return 0;
}

long
BRLObolWindowHost::pollRate(void) const
{
    return this->p->pollRate;
}

int
BRLObolWindowHost::openFramebuffer(imgstream_fb_t *fb,
	const struct imgstream_fb_spec_info *info)
{
    if (!fb)
	return -1;
    if (find_attachment(this->p->framebuffers, fb))
	return 0;

    BRLObolWindowDesc desc = this->p->open ? this->p->desc : default_desc();
    desc.width = (unsigned int)std::max<size_t>(imgstream_fb_width(fb), 1);
    desc.height = (unsigned int)std::max<size_t>(imgstream_fb_height(fb), 1);
    desc.mode = (info && info->display == IMGSTREAM_FB_DISPLAY_SWRAST) ?
	BRLOBOL_WINDOW_HEADLESS : BRLOBOL_WINDOW_TOPLEVEL;
    desc.backend = info ? backend_from_fb_display(info->display) :
	BRLOBOL_WINDOW_BACKEND_AUTO;
    desc.visible = desc.mode == BRLOBOL_WINDOW_TOPLEVEL;
    if (info && info->host && info->host_len)
	desc.display = SbString(info->host, 0, (int)info->host_len - 1);
    desc.title = imgstream_fb_name(fb) ? imgstream_fb_name(fb) : "BRL-CAD framebuffer";

    if (this->open(&desc) != 0)
	return -1;

    imgstream_t *stream = imgstream_fb_stream(fb);
    if (!stream)
	return -1;

    SoGroup *group = controller_root_group(this->p->controller);
    if (!group)
	return -1;

    SoBRLImageSource *source = new SoBRLImageSource;
    source->ref();
    source->imageId = imgstream_fb_name(fb) ? imgstream_fb_name(fb) : "framebuffer";
    source->sourceUri = source->imageId.getValue();
    if (source->setStream(stream) != 0) {
	source->unref();
	return -1;
    }

    SoBRLViewportImage *viewport = new SoBRLViewportImage;
    viewport->ref();
    viewport->overlayId = source->imageId.getValue();
    viewport->imageSource.setValue(source);
    viewport->layer = SoBRLViewportImage::HUD;
    viewport->anchor = SoBRLViewportImage::LOWER_LEFT;
    viewport->fit = SoBRLViewportImage::STRETCH;
    viewport->preserveAspect = FALSE;
    viewport->position.setValue(0.0f, 0.0f);
    viewport->size.setValue((float)desc.width, (float)desc.height);
    viewport->sourceCenter.setValue((float)desc.width * 0.5f,
	    (float)desc.height * 0.5f);
    viewport->sourceZoom = 1.0f;
    viewport->cursorVisible = FALSE;
    if (viewport->rebuildGeometry() != 0) {
	viewport->unref();
	source->unref();
	return -1;
    }

    group->addChild(source);
    group->addChild(viewport);

    BRLObolFramebufferAttachment attachment;
    attachment.fb = fb;
    attachment.source = source;
    attachment.viewport = viewport;
    this->p->framebuffers.push_back(attachment);
    this->p->controller->requestRender("fb-open");
    return 0;
}

void
BRLObolWindowHost::closeFramebuffer(imgstream_fb_t *fb)
{
    for (size_t i = 0; i < this->p->framebuffers.size(); i++) {
	if (this->p->framebuffers[i].fb != fb)
	    continue;

	SoGroup *group = controller_root_group(this->p->controller);
	remove_child_if_present(group, this->p->framebuffers[i].viewport);
	remove_child_if_present(group, this->p->framebuffers[i].source);
	this->p->framebuffers[i].viewport->unref();
	this->p->framebuffers[i].source->unref();
	this->p->framebuffers.erase(this->p->framebuffers.begin() + (ptrdiff_t)i);
	if (this->p->controller)
	    this->p->controller->requestRender("fb-close");
	return;
    }
}

int
BRLObolWindowHost::flushFramebuffer(imgstream_fb_t *fb)
{
    BRLObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment)
	return -1;

    if (attachment->source->refreshFromStream() != 0 ||
	    attachment->viewport->syncFromSource() != 0)
	return -1;
    if (this->p->controller)
	this->p->controller->requestRender("fb-flush");
    return 0;
}

int
BRLObolWindowHost::resetFramebuffer(imgstream_fb_t *fb)
{
    BRLObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment)
	return -1;

    attachment->viewport->sourceCenter.setValue(
	    (float)imgstream_fb_width(fb) * 0.5f,
	    (float)imgstream_fb_height(fb) * 0.5f);
    attachment->viewport->sourceZoom = 1.0f;
    attachment->viewport->cursorVisible = FALSE;
    if (attachment->viewport->syncFromSource() != 0)
	return -1;
    if (this->p->controller)
	this->p->controller->requestRender("fb-reset");
    return 0;
}

int
BRLObolWindowHost::setFramebufferViewport(imgstream_fb_t *fb,
	int left, int top, int right, int bottom)
{
    BRLObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment)
	return -1;

    int width = right - left;
    int height = bottom - top;
    if (width <= 0)
	width = (int)imgstream_fb_width(fb);
    if (height <= 0)
	height = (int)imgstream_fb_height(fb);
    if (width <= 0)
	width = 1;
    if (height <= 0)
	height = 1;

    if (this->p->controller)
	this->p->controller->setViewportSize((unsigned int)width,
		(unsigned int)height);
    attachment->viewport->position.setValue((float)left, (float)top);
    attachment->viewport->size.setValue((float)width, (float)height);
    if (attachment->viewport->syncFromSource() != 0)
	return -1;
    if (this->p->controller)
	this->p->controller->requestRender("fb-viewport");
    return 0;
}

int
BRLObolWindowHost::setFramebufferView(imgstream_fb_t *fb,
	const struct imgstream_fb_view *view)
{
    BRLObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment || !view)
	return -1;

    attachment->viewport->sourceCenter.setValue((float)view->xcenter,
	    (float)view->ycenter);
    attachment->viewport->sourceZoom = positive_zoom(view->xzoom, view->yzoom);
    if (this->p->controller)
	this->p->controller->requestRender("fb-view");
    return 0;
}

int
BRLObolWindowHost::setFramebufferCursor(imgstream_fb_t *fb,
	const struct imgstream_fb_cursor *cursor)
{
    BRLObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment || !cursor)
	return -1;

    attachment->viewport->cursorVisible = cursor->mode ? TRUE : FALSE;
    attachment->viewport->cursorImagePosition.setValue((float)cursor->x,
	    (float)cursor->y);
    attachment->viewport->cursorShape = cursor->mode ?
	SoBRLViewportImage::CURSOR_DEFAULT : SoBRLViewportImage::CURSOR_NONE;
    if (this->p->controller)
	this->p->controller->requestRender("fb-cursor");
    return 0;
}

int
BRLObolWindowHost::setFramebufferScreenCursor(imgstream_fb_t *fb,
	int mode, int x, int y)
{
    BRLObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment)
	return -1;

    attachment->viewport->cursorVisible = mode ? TRUE : FALSE;
    attachment->viewport->cursorImagePosition.setValue((float)x, (float)y);
    attachment->viewport->cursorShape = mode ?
	SoBRLViewportImage::CURSOR_DEFAULT : SoBRLViewportImage::CURSOR_NONE;
    if (this->p->controller)
	this->p->controller->requestRender("fb-screen-cursor");
    return 0;
}

int
BRLObolWindowHost::setFramebufferCursorShape(imgstream_fb_t *fb,
	const unsigned char *bits, int xbits, int ybits,
	int UNUSED(xorig), int UNUSED(yorig))
{
    BRLObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment)
	return -1;

    if (bits && xbits > 0 && ybits > 0) {
	attachment->viewport->cursorVisible = TRUE;
	attachment->viewport->cursorShape = SoBRLViewportImage::CURSOR_CUSTOM;
    } else {
	attachment->viewport->cursorShape = attachment->viewport->cursorVisible.getValue() ?
	    SoBRLViewportImage::CURSOR_DEFAULT : SoBRLViewportImage::CURSOR_NONE;
    }
    if (this->p->controller)
	this->p->controller->requestRender("fb-cursor-shape");
    return 0;
}

int
BRLObolWindowHost::getFramebufferCount(void) const
{
    return (int)this->p->framebuffers.size();
}

SoBRLImageSource *
BRLObolWindowHost::getFramebufferImageSource(imgstream_fb_t *fb) const
{
    const BRLObolFramebufferAttachment *attachment =
	find_attachment_const(this->p->framebuffers, fb);
    return attachment ? attachment->source : NULL;
}

SoBRLViewportImage *
BRLObolWindowHost::getFramebufferViewportImage(imgstream_fb_t *fb) const
{
    const BRLObolFramebufferAttachment *attachment =
	find_attachment_const(this->p->framebuffers, fb);
    return attachment ? attachment->viewport : NULL;
}

static int
brlobol_fb_host_open(imgstream_fb_t *fb,
	const struct imgstream_fb_spec_info *info,
	void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    return host ? host->openFramebuffer(fb, info) : -1;
}

static void
brlobol_fb_host_close(imgstream_fb_t *fb, void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    if (host)
	host->closeFramebuffer(fb);
}

static int
brlobol_fb_host_flush(imgstream_fb_t *fb, void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    return host ? host->flushFramebuffer(fb) : -1;
}

static int
brlobol_fb_host_reset(imgstream_fb_t *fb, void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    return host ? host->resetFramebuffer(fb) : -1;
}

static int
brlobol_fb_host_viewport(imgstream_fb_t *fb,
	int left, int top, int right, int bottom,
	void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    return host ? host->setFramebufferViewport(fb, left, top, right, bottom) : -1;
}

static int
brlobol_fb_host_view(imgstream_fb_t *fb,
	const struct imgstream_fb_view *view,
	void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    return host ? host->setFramebufferView(fb, view) : -1;
}

static int
brlobol_fb_host_cursor(imgstream_fb_t *fb,
	const struct imgstream_fb_cursor *cursor,
	void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    return host ? host->setFramebufferCursor(fb, cursor) : -1;
}

static int
brlobol_fb_host_scursor(imgstream_fb_t *fb,
	int mode, int x, int y,
	void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    return host ? host->setFramebufferScreenCursor(fb, mode, x, y) : -1;
}

static int
brlobol_fb_host_setcursor(imgstream_fb_t *fb,
	const unsigned char *bits, int xbits, int ybits,
	int xorig, int yorig,
	void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    return host ? host->setFramebufferCursorShape(fb, bits, xbits, ybits,
	    xorig, yorig) : -1;
}

static int
brlobol_fb_host_poll(imgstream_fb_t *UNUSED(fb), void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    return host ? host->poll(NULL) : -1;
}

static long
brlobol_fb_host_poll_rate(const imgstream_fb_t *UNUSED(fb), void *data)
{
    BRLObolWindowHost *host = static_cast<BRLObolWindowHost *>(data);
    return host ? host->pollRate() : 0;
}

int
brlobol_window_host_register_display_host(BRLObolWindowHost *host)
{
    if (!host)
	return -1;

    struct imgstream_fb_display_host displayHost;
    memset(&displayHost, 0, sizeof(displayHost));
    displayHost.open = brlobol_fb_host_open;
    displayHost.close = brlobol_fb_host_close;
    displayHost.flush = brlobol_fb_host_flush;
    displayHost.reset = brlobol_fb_host_reset;
    displayHost.viewport = brlobol_fb_host_viewport;
    displayHost.view = brlobol_fb_host_view;
    displayHost.cursor = brlobol_fb_host_cursor;
    displayHost.scursor = brlobol_fb_host_scursor;
    displayHost.setcursor = brlobol_fb_host_setcursor;
    displayHost.poll = brlobol_fb_host_poll;
    displayHost.poll_rate = brlobol_fb_host_poll_rate;
    return imgstream_fb_display_host_set(&displayHost, host);
}

void
brlobol_window_host_unregister_display_host(void)
{
    imgstream_fb_display_host_clear();
}
