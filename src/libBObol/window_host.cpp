/*                W I N D O W _ H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol/BImageSource.h"
#include "BObol/BFramebuffer.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BViewController.h"
#include "BObol/BViewportImage.h"
#include "BObol/BWindowHost.h"

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
#include <cmath>
#include <string.h>
#include <vector>

struct BObolFramebufferAttachment {
    imgstream_fb_t *fb;
    SoBRLImageSource *source;
    SoBRLViewportImage *viewport;
    BObolFramebufferComposition composition;
};

struct BObolWindowHostPrivate {
    BObolWindowHostPrivate(void) :
	controller(new BObolViewController),
	ownsController(TRUE),
	open(FALSE),
	pollRate(0)
    {
	desc.mode = BOBOL_WINDOW_HEADLESS;
	desc.backend = BOBOL_WINDOW_BACKEND_OFFSCREEN;
	desc.width = 1;
	desc.height = 1;
	desc.title = "BRL-CAD Obol";
	desc.display = "";
	desc.nativeIdHint = "";
	desc.visible = FALSE;
    }

    BObolViewController *controller;
    SbBool ownsController;
    SbBool open;
    BObolWindowDesc desc;
    long pollRate;
    BObolInputContext input;
    std::vector<BObolFramebufferAttachment> framebuffers;
    std::vector<BObolFramebufferStream *> framebufferStreams;
    std::vector<imgstream_fb_t *> displayFramebuffers;
};

static BObolWindowDesc
default_desc(void)
{
    BObolWindowDesc desc;
    desc.mode = BOBOL_WINDOW_HEADLESS;
    desc.backend = BOBOL_WINDOW_BACKEND_OFFSCREEN;
    desc.width = 1;
    desc.height = 1;
    desc.title = "BRL-CAD Obol";
    desc.display = "";
    desc.nativeIdHint = "";
    desc.visible = FALSE;
    return desc;
}

static void
sanitize_desc(BObolWindowDesc *desc)
{
    if (!desc)
	return;
    if (desc->width == 0)
	desc->width = 1;
    if (desc->height == 0)
	desc->height = 1;
}

static BObolWindowBackend
backend_from_fb_display(enum imgstream_fb_display_kind display)
{
    switch (display) {
	case IMGSTREAM_FB_DISPLAY_QTGL:
	    return BOBOL_WINDOW_BACKEND_QT;
	case IMGSTREAM_FB_DISPLAY_SWRAST:
	    return BOBOL_WINDOW_BACKEND_OFFSCREEN;
	case IMGSTREAM_FB_DISPLAY_X:
	    return BOBOL_WINDOW_BACKEND_TK;
	case IMGSTREAM_FB_DISPLAY_OGL:
	case IMGSTREAM_FB_DISPLAY_WGL:
	    return BOBOL_WINDOW_BACKEND_OPENGL;
	case IMGSTREAM_FB_DISPLAY_NONE:
	default:
	    return BOBOL_WINDOW_BACKEND_AUTO;
    }
}

static SoGroup *
controller_root_group(BObolViewController *controller)
{
    if (!controller)
	return NULL;

    SoNode *root = controller->getSceneRoot();
    if (!root) {
	SoBRLSceneGroup *newRoot = new SoBRLSceneGroup;
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

static void
remove_framebuffer_viewport(BObolViewController *controller,
	SoBRLViewportImage *viewport)
{
    if (!controller || !viewport)
	return;
    remove_child_if_present(controller->getFramebufferUnderlayRoot(), viewport);
    remove_child_if_present(controller->getFramebufferInterlayRoot(), viewport);
    remove_child_if_present(controller->getFramebufferOverlayRoot(), viewport);
}

static float
positive_zoom(int xzoom, int yzoom)
{
    int x = xzoom > 0 ? xzoom : 1;
    int y = yzoom > 0 ? yzoom : 1;
    return (float)((x < y) ? x : y);
}

static bool
same_float(float a, float b)
{
    return std::fabs(a - b) <= 1.0e-6f;
}

static BObolFramebufferAttachment *
find_attachment(std::vector<BObolFramebufferAttachment> &attachments,
		imgstream_fb_t *fb)
{
    for (size_t i = 0; i < attachments.size(); i++) {
	if (attachments[i].fb == fb)
	    return &attachments[i];
    }
    return NULL;
}

static const BObolFramebufferAttachment *
find_attachment_const(const std::vector<BObolFramebufferAttachment> &attachments,
		      imgstream_fb_t *fb)
{
    for (size_t i = 0; i < attachments.size(); i++) {
	if (attachments[i].fb == fb)
	    return &attachments[i];
    }
    return NULL;
}

BObolWindowHost::BObolWindowHost(void) :
    p(new BObolWindowHostPrivate)
{
}

BObolWindowHost::~BObolWindowHost(void)
{
    this->detachDisplayFramebuffers();
    this->detachFramebufferStreams();
    this->close();
    if (this->p->ownsController)
	delete this->p->controller;
    this->p->controller = NULL;
    delete this->p;
}

void
BObolWindowHost::registerFramebufferStream(BObolFramebufferStream *stream)
{
    if (!stream || std::find(this->p->framebufferStreams.begin(),
	this->p->framebufferStreams.end(), stream) !=
	this->p->framebufferStreams.end())
	return;
    this->p->framebufferStreams.push_back(stream);
}

void
BObolWindowHost::unregisterFramebufferStream(BObolFramebufferStream *stream)
{
    this->p->framebufferStreams.erase(std::remove_if(
	this->p->framebufferStreams.begin(), this->p->framebufferStreams.end(),
	[stream](BObolFramebufferStream *candidate) {
	    return candidate == stream;
	}),
	this->p->framebufferStreams.end());
}

void
BObolWindowHost::registerDisplayFramebuffer(imgstream_fb_t *fb)
{
    if (!fb || std::find(this->p->displayFramebuffers.begin(),
	this->p->displayFramebuffers.end(), fb) != this->p->displayFramebuffers.end())
	return;
    this->p->displayFramebuffers.push_back(fb);
}

void
BObolWindowHost::unregisterDisplayFramebuffer(imgstream_fb_t *fb)
{
    this->p->displayFramebuffers.erase(std::remove_if(
	this->p->displayFramebuffers.begin(), this->p->displayFramebuffers.end(),
	[fb](imgstream_fb_t *candidate) {
	    return candidate == fb;
	}),
	this->p->displayFramebuffers.end());
}

void
BObolWindowHost::detachFramebufferStreams(void)
{
    for (size_t i = 0; i < this->p->framebufferStreams.size(); i++)
	this->p->framebufferStreams[i]->detachHost(this);
    this->p->framebufferStreams.clear();
}

void
BObolWindowHost::detachDisplayFramebuffers(void)
{
    while (!this->p->displayFramebuffers.empty()) {
	imgstream_fb_t *fb = this->p->displayFramebuffers.back();
	this->p->displayFramebuffers.pop_back();
	(void)imgstream_fb_detach_display_host(fb);
    }
}

int
BObolWindowHost::open(const BObolWindowDesc *desc)
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
BObolWindowHost::close(void)
{
    while (!this->p->framebuffers.empty())
	this->closeFramebuffer(this->p->framebuffers.back().fb);
    this->p->open = FALSE;
}

SbBool
BObolWindowHost::isOpen(void) const
{
    return this->p->open;
}

const BObolWindowDesc &
BObolWindowHost::getDesc(void) const
{
    return this->p->desc;
}

void
BObolWindowHost::attachController(BObolViewController *controller,
				    SbBool takeOwnership)
{
    this->close();
    if (this->p->ownsController)
	delete this->p->controller;
    this->p->controller = controller;
    this->p->ownsController = takeOwnership;
}

BObolViewController *
BObolWindowHost::getController(void) const
{
    return this->p->controller;
}

static int
window_host_input_action(void *data, BObolInputAction action,
	const BObolInputEvent *event)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->applyInputAction(action, event) : -1;
}

int
BObolWindowHost::handleInputEvent(const BObolInputEvent *event,
				    const BObolInputProfile *profile)
{
    if (!event)
	return -1;
    if (!profile || !profile->bindings || profile->bindingCount == 0)
	return 0;

    this->p->input.setProfile(profile);
    this->p->input.setActionHandler(window_host_input_action, this);
    return this->p->input.dispatch(event);
}

int
BObolWindowHost::applyInputAction(BObolInputAction action,
				    const BObolInputEvent *UNUSED(event))
{
    if (action == BOBOL_ACTION_NONE)
	return 0;
    if (!this->p->controller)
	return -1;

    this->p->controller->requestRender("input-action");
    return 1;
}

int
BObolWindowHost::poll(const BObolInputProfile *UNUSED(profile))
{
    return 0;
}

long
BObolWindowHost::pollRate(void) const
{
    return this->p->pollRate;
}

int
BObolWindowHost::openFramebuffer(imgstream_fb_t *fb,
				   const imgstream_fb_spec_info_t *info)
{
    if (!fb)
	return -1;
    if (find_attachment(this->p->framebuffers, fb))
	return 0;

    BObolWindowDesc desc = this->p->open ? this->p->desc : default_desc();
    desc.width = (unsigned int)std::max<size_t>(imgstream_fb_width(fb), 1);
    desc.height = (unsigned int)std::max<size_t>(imgstream_fb_height(fb), 1);
    desc.mode = (info && info->display == IMGSTREAM_FB_DISPLAY_SWRAST) ?
		BOBOL_WINDOW_HEADLESS : BOBOL_WINDOW_TOPLEVEL;
    desc.backend = info ? backend_from_fb_display(info->display) :
		   BOBOL_WINDOW_BACKEND_AUTO;
    desc.visible = desc.mode == BOBOL_WINDOW_TOPLEVEL;
    if (info && info->host && info->host_len)
	desc.display = SbString(info->host, 0, (int)info->host_len - 1);
    desc.title = imgstream_fb_name(fb) ? imgstream_fb_name(fb) : "BRL-CAD framebuffer";

    if (this->open(&desc) != 0)
	return -1;

    imgstream_t *stream = imgstream_fb_stream(fb);
    if (!stream)
	return -1;

    if (!controller_root_group(this->p->controller) ||
	!this->p->controller->getFramebufferOverlayRoot())
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
    viewport->layer = SoBRLViewportImage::OVERLAY;
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

    BObolFramebufferAttachment attachment;
    attachment.fb = fb;
    attachment.source = source;
    attachment.viewport = viewport;
    attachment.composition = BOBOL_FRAMEBUFFER_COMPOSITION_OVERLAY;
    this->p->framebuffers.push_back(attachment);
	this->p->controller->getFramebufferOverlayRoot()->addChild(viewport);
    this->p->controller->requestRender("fb-open");
    return 0;
}

void
BObolWindowHost::closeFramebuffer(imgstream_fb_t *fb)
{
    for (size_t i = 0; i < this->p->framebuffers.size(); i++) {
	if (this->p->framebuffers[i].fb != fb)
	    continue;

	remove_framebuffer_viewport(this->p->controller,
		this->p->framebuffers[i].viewport);
	this->p->framebuffers[i].viewport->unref();
	this->p->framebuffers[i].source->unref();
	this->p->framebuffers.erase(this->p->framebuffers.begin() + (ptrdiff_t)i);
	if (this->p->controller)
	    this->p->controller->requestRender("fb-close");
	return;
    }
}

int
BObolWindowHost::setFramebufferComposition(imgstream_fb_t *fb,
	BObolFramebufferComposition composition)
{
    if (composition < BOBOL_FRAMEBUFFER_COMPOSITION_OFF ||
	composition > BOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY)
	return -1;

    BObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    BObolViewController *controller = this->p->controller;
    if (!attachment || !controller)
	return -1;
    /*
     * The framebuffer bridge may synchronize its policy on every canvas
     * presentation.  Reapplying an unchanged mode used to remove, rebuild,
     * and reinsert the viewport and leave another "fb-composition" render
     * request behind each OSMesa frame.  Besides wasting owner-thread work,
     * that makes a fully settled progressive view observably never idle.
     */
    if (attachment->composition == composition)
	return 0;

    SoBRLViewportImage *viewport = attachment->viewport;
    remove_framebuffer_viewport(controller, viewport);
    if (composition == BOBOL_FRAMEBUFFER_COMPOSITION_OFF) {
	viewport->visible = FALSE;
	if (viewport->rebuildGeometry() != 0)
	    return -1;
	attachment->composition = composition;
	controller->requestRender("fb-composition");
	return 0;
    }

    SoGroup *layer = NULL;
    int viewport_layer = SoBRLViewportImage::OVERLAY;
    switch (composition) {
	case BOBOL_FRAMEBUFFER_COMPOSITION_UNDERLAY:
	    layer = controller->getFramebufferUnderlayRoot();
	    viewport_layer = SoBRLViewportImage::UNDERLAY;
	    break;
	case BOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY:
	    layer = controller->getFramebufferInterlayRoot();
	    viewport_layer = SoBRLViewportImage::INTERLAY;
	    break;
	case BOBOL_FRAMEBUFFER_COMPOSITION_OVERLAY:
	    layer = controller->getFramebufferOverlayRoot();
	    break;
	default:
	    return -1;
    }
    if (!layer)
	return -1;
    viewport->visible = TRUE;
    viewport->layer = viewport_layer;
    if (viewport->rebuildGeometry() != 0)
	return -1;
    layer->addChild(viewport);
    attachment->composition = composition;
    controller->requestRender("fb-composition");
    return 0;
}

int
BObolWindowHost::flushFramebuffer(imgstream_fb_t *fb)
{
    BObolFramebufferAttachment *attachment =
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
BObolWindowHost::resetFramebuffer(imgstream_fb_t *fb)
{
    BObolFramebufferAttachment *attachment =
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
BObolWindowHost::setFramebufferViewport(imgstream_fb_t *fb,
	int left, int top, int right, int bottom)
{
    BObolFramebufferAttachment *attachment =
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
BObolWindowHost::setFramebufferView(imgstream_fb_t *fb,
				      const imgstream_fb_view_t *view)
{
    BObolFramebufferAttachment *attachment =
	find_attachment(this->p->framebuffers, fb);
    if (!attachment || !view)
	return -1;

    float zoom = positive_zoom(view->xzoom, view->yzoom);
    SbVec2f oldCenter = attachment->viewport->sourceCenter.getValue();
    if (same_float(oldCenter[0], (float)view->xcenter) &&
	    same_float(oldCenter[1], (float)view->ycenter) &&
	    same_float(attachment->viewport->sourceZoom.getValue(), zoom))
	return 0;

    attachment->viewport->sourceCenter.setValue((float)view->xcenter,
				      (float)view->ycenter);
    attachment->viewport->sourceZoom = zoom;
    if (attachment->viewport->rebuildGeometry() != 0)
	return -1;
    if (this->p->controller)
	this->p->controller->requestRender("fb-view");
    return 0;
}

int
BObolWindowHost::setFramebufferCursor(imgstream_fb_t *fb,
					const imgstream_fb_cursor_t *cursor)
{
    BObolFramebufferAttachment *attachment =
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
BObolWindowHost::setFramebufferScreenCursor(imgstream_fb_t *fb,
	int mode, int x, int y)
{
    BObolFramebufferAttachment *attachment =
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
BObolWindowHost::setFramebufferCursorShape(imgstream_fb_t *fb,
	const unsigned char *bits, int xbits, int ybits,
	int UNUSED(xorig), int UNUSED(yorig))
{
    BObolFramebufferAttachment *attachment =
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
BObolWindowHost::getFramebufferCount(void) const
{
    return (int)this->p->framebuffers.size();
}

SoBRLImageSource *
BObolWindowHost::getFramebufferImageSource(imgstream_fb_t *fb) const
{
    const BObolFramebufferAttachment *attachment =
	find_attachment_const(this->p->framebuffers, fb);
    return attachment ? attachment->source : NULL;
}

SoBRLViewportImage *
BObolWindowHost::getFramebufferViewportImage(imgstream_fb_t *fb) const
{
    const BObolFramebufferAttachment *attachment =
	find_attachment_const(this->p->framebuffers, fb);
    return attachment ? attachment->viewport : NULL;
}

struct BObolFramebufferBridge {
    static int
    open(imgstream_fb_t *fb, const imgstream_fb_spec_info_t *info, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    if (!host || host->openFramebuffer(fb, info) != 0)
	return -1;
    host->registerDisplayFramebuffer(fb);
    return 0;
}

    static void
    close(imgstream_fb_t *fb, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    if (host) {
	host->unregisterDisplayFramebuffer(fb);
	host->closeFramebuffer(fb);
	    }
}

    static int
    flush(imgstream_fb_t *fb, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->flushFramebuffer(fb) : -1;
}

    static int
    reset(imgstream_fb_t *fb, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->resetFramebuffer(fb) : -1;
}

    static int
    viewport(imgstream_fb_t *fb, int left, int top, int right, int bottom,
	     void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->setFramebufferViewport(fb, left, top, right, bottom) : -1;
}

    static int
    view(imgstream_fb_t *fb, const imgstream_fb_view_t *view, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->setFramebufferView(fb, view) : -1;
}

    static int
    cursor(imgstream_fb_t *fb, const imgstream_fb_cursor_t *cursor, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->setFramebufferCursor(fb, cursor) : -1;
}

    static int
    scursor(imgstream_fb_t *fb, int mode, int x, int y, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->setFramebufferScreenCursor(fb, mode, x, y) : -1;
}

    static int
    setcursor(imgstream_fb_t *fb, const unsigned char *bits, int xbits, int ybits,
	      int xorig, int yorig, void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->setFramebufferCursorShape(fb, bits, xbits, ybits,
	    xorig, yorig) : -1;
}

    static int
    poll(imgstream_fb_t *UNUSED(fb), void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->poll(NULL) : -1;
}

    static long
    pollRate(const imgstream_fb_t *UNUSED(fb), void *data)
{
    BObolWindowHost *host = static_cast<BObolWindowHost *>(data);
    return host ? host->pollRate() : 0;
}
};

imgstream_fb_t *
bobol_window_host_open_display_framebuffer(BObolWindowHost *host,
	const char *spec, size_t width, size_t height)
{
    if (!host)
	return NULL;

    struct imgstream_fb_display_host displayHost;
    memset(&displayHost, 0, sizeof(displayHost));
    displayHost.open = BObolFramebufferBridge::open;
    displayHost.close = BObolFramebufferBridge::close;
    displayHost.flush = BObolFramebufferBridge::flush;
    displayHost.reset = BObolFramebufferBridge::reset;
    displayHost.viewport = BObolFramebufferBridge::viewport;
    displayHost.view = BObolFramebufferBridge::view;
    displayHost.cursor = BObolFramebufferBridge::cursor;
    displayHost.scursor = BObolFramebufferBridge::scursor;
    displayHost.setcursor = BObolFramebufferBridge::setcursor;
    displayHost.poll = BObolFramebufferBridge::poll;
    displayHost.poll_rate = BObolFramebufferBridge::pollRate;
    return imgstream_fb_open_display(spec, width, height, &displayHost, host);
}
