/*                     O B O L _ F B S E R V . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libged/obol_fbserv.cpp
 *
 * Shared libged bridge from the fbserv protocol backend table to an
 * Obol/image-stream framebuffer.  This is intentionally a libged helper rather
 * than a renderer behavior: fbserv remains a protocol transport for now, while Obol
 * owns presentation of the resulting image stream.
 */

#include "common.h"

#include <atomic>
#include <climits>
#include <cerrno>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include <thread>

#ifndef _WIN32
#  include <sys/select.h>
#  include <sys/time.h>
#  include <unistd.h>
#endif

#include "brlobol/framebuffer.h"
#include "brlobol/display_endpoint.h"
#include "brlobol/view_controller.h"
#include "brlobol/viewport_image.h"
#include "brlobol/window_host.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "bv.h"
#include "imgstream/fbserv.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "imgstream/fb_compat.h"

#include "./ged_private.h"

class BRLObolViewController;

namespace {

class GedObolFbservBridge;

static int ged_obol_capture_framebuffer(void *ctx, unsigned char **pixels,
	size_t *size, unsigned int *width, unsigned int *height,
	unsigned int *components);

static int ged_obol_fb_info(void *ctx, struct fbserv_fb_info *info);
static int ged_obol_fb_clear(void *ctx, const unsigned char rgb[3]);
static ssize_t ged_obol_fb_read(void *ctx, int x, int y,
	unsigned char *rgb, size_t count);
static ssize_t ged_obol_fb_write(void *ctx, int x, int y,
	const unsigned char *rgb, size_t count);
static int ged_obol_fb_readrect(void *ctx, int xmin, int ymin, int width,
	int height, unsigned char *rgb);
static int ged_obol_fb_writerect(void *ctx, int xmin, int ymin, int width,
	int height, const unsigned char *rgb);
static int ged_obol_fb_bwreadrect(void *ctx, int xmin, int ymin, int width,
	int height, unsigned char *bw);
static int ged_obol_fb_bwwriterect(void *ctx, int xmin, int ymin, int width,
	int height, const unsigned char *bw);
static int ged_obol_fb_cursor(void *ctx, int mode, int x, int y);
static int ged_obol_fb_getcursor(void *ctx, int *mode, int *x, int *y);
static int ged_obol_fb_setcursor(void *ctx, const unsigned char *bits,
	int xbits, int ybits, int xorig, int yorig);
static int ged_obol_fb_scursor(void *ctx, int mode, int x, int y);
static int ged_obol_fb_window(void *ctx, int xcenter, int ycenter);
static int ged_obol_fb_zoom(void *ctx, int xzoom, int yzoom);
static int ged_obol_fb_view(void *ctx, int xcenter, int ycenter, int xzoom,
	int yzoom);
static int ged_obol_fb_getview(void *ctx, int *xcenter, int *ycenter,
	int *xzoom, int *yzoom);
static int ged_obol_fb_rmap(void *ctx, struct fbserv_colormap *cmap);
static int ged_obol_fb_wmap(void *ctx, const struct fbserv_colormap *cmap);
static int ged_obol_fb_flush(void *ctx);
static int ged_obol_fb_poll(void *ctx);
static int ged_obol_fb_help(void *ctx);

static const struct fbserv_fb_ops ged_obol_fb_ops = {
    ged_obol_fb_info,
    ged_obol_fb_clear,
    ged_obol_fb_read,
    ged_obol_fb_write,
    ged_obol_fb_readrect,
    ged_obol_fb_writerect,
    ged_obol_fb_bwreadrect,
    ged_obol_fb_bwwriterect,
    ged_obol_fb_cursor,
    ged_obol_fb_getcursor,
    ged_obol_fb_setcursor,
    ged_obol_fb_scursor,
    ged_obol_fb_window,
    ged_obol_fb_zoom,
    ged_obol_fb_view,
    ged_obol_fb_getview,
    ged_obol_fb_rmap,
    ged_obol_fb_wmap,
    ged_obol_fb_flush,
    ged_obol_fb_poll,
    ged_obol_fb_help
};

struct GedObolFbservWatcher {
    std::atomic<bool> stop{false};
    std::thread thread;
    int slot = -1;
};

class GedObolFbservBridge {
public:
    explicit GedObolFbservBridge(struct fbserv_obj *fbs_obj) :
	fbs(fbs_obj),
	capture_endpoint(NULL),
	framebuffer(&default_host)
    {
    }

    ~GedObolFbservBridge()
    {
	stopWatchers();
	clearCaptureProvider();
	framebuffer.close();
    }

    int configure(struct ged *new_gedp, void *new_view_ctx,
	    BRLObolViewController *controller,
	    BRLObolWindowHost *window_host,
	    int requested_width,
	    int requested_height,
	    int requested_present_on_flush)
    {
	std::lock_guard<std::mutex> guard(lock);
	if (!new_gedp || !new_view_ctx || !controller)
	    return BRLCAD_ERROR;

	/* Keep the old capture provider bound until the target host has accepted
	 * the live stream.  This makes a rejected replacement transactional across
	 * both retained image nodes and capture ownership. */
	const bool view_changed = view_ctx != new_view_ctx;
	const struct bv *view = bv_context_view_const(
	    (const struct bv_context *)new_view_ctx);
	int composition = view ? bv_framebuffer_mode_get(view) : 0;
	if (composition < BRLOBOL_FRAMEBUFFER_COMPOSITION_OFF ||
	    composition > BRLOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY)
	    composition = BRLOBOL_FRAMEBUFFER_COMPOSITION_OFF;
	BRLObolWindowHost *active_host = window_host ? window_host :
	    &default_host;
	if (active_host->getController() != controller) {
	    /* Changing a host controller closes its retained framebuffer nodes.
	     * Close the stream first so it can reopen and attach to the new
	     * controller even when the generic host pointer itself is unchanged. */
	    framebuffer.close();
	    active_host->attachController(controller, FALSE);
	}
	if (framebuffer.setHost(active_host,
		static_cast<BRLObolFramebufferComposition>(composition)) != 0)
	    return BRLCAD_ERROR;

	int width = requested_width;
	int height = requested_height;
	if (view) {
	    if (width <= 0)
		width = bv_width_get(view);
	    if (height <= 0)
		height = bv_height_get(view);
	}
	if (width <= 0)
	    width = 512;
	if (height <= 0)
	    height = 512;

	if (framebuffer.configure(width, height) != 0)
	    return BRLCAD_ERROR;
	if (framebuffer.ensure() != 0)
	    return BRLCAD_ERROR;
	if (setCompositionLocked(composition) != BRLCAD_OK)
	    return BRLCAD_ERROR;

	if (view_changed)
	    clearCaptureProviderLocked();
	gedp = new_gedp;
	view_ctx = new_view_ctx;
	present_on_flush = requested_present_on_flush ? 1 : 0;
	bindCaptureProviderLocked();
	return BRLCAD_OK;
    }

    int info(struct fbserv_fb_info *fbinfo)
    {
	std::lock_guard<std::mutex> guard(lock);
	if (!fbinfo || framebuffer.ensure() != 0)
	    return -1;

	BRLObolFramebufferInfo info;
	if (framebuffer.info(&info) != 0)
	    return -1;
	fbinfo->max_width = info.max_width;
	fbinfo->max_height = info.max_height;
	fbinfo->width = info.width;
	fbinfo->height = info.height;
	return 0;
    }

    int clear(const unsigned char rgb[3])
    {
	std::lock_guard<std::mutex> guard(lock);
	int ret = framebuffer.clear(rgb);
	if (ret == 0)
	    notifyUpdatedLocked();
	return ret;
    }

    ssize_t read(int x, int y, unsigned char *rgb, size_t count)
    {
	std::lock_guard<std::mutex> guard(lock);
	return framebuffer.read(x, y, rgb, count);
    }

    ssize_t write(int x, int y, const unsigned char *rgb, size_t count)
    {
	std::lock_guard<std::mutex> guard(lock);
	ssize_t ret = framebuffer.write(x, y, rgb, count);
	if (ret >= 0)
	    notifyUpdatedLocked();
	return ret;
    }

    int readrect(int xmin, int ymin, int width, int height,
	    unsigned char *rgb)
    {
	std::lock_guard<std::mutex> guard(lock);
	return framebuffer.readrect(xmin, ymin, width, height, rgb);
    }

    int writerect(int xmin, int ymin, int width, int height,
	    const unsigned char *rgb)
    {
	std::lock_guard<std::mutex> guard(lock);
	int ret = framebuffer.writerect(xmin, ymin, width, height, rgb);
	if (ret >= 0)
	    notifyUpdatedLocked();
	return ret;
    }

    int bwreadrect(int xmin, int ymin, int width, int height,
	    unsigned char *bw)
    {
	std::lock_guard<std::mutex> guard(lock);
	return framebuffer.bwreadrect(xmin, ymin, width, height, bw);
    }

    int bwwriterect(int xmin, int ymin, int width, int height,
	    const unsigned char *bw)
    {
	std::lock_guard<std::mutex> guard(lock);
	int ret = framebuffer.bwwriterect(xmin, ymin, width, height, bw);
	if (ret >= 0)
	    notifyUpdatedLocked();
	return ret;
    }

    int cursor(int mode, int x, int y)
    {
	std::lock_guard<std::mutex> guard(lock);
	int ret = framebuffer.cursor(mode, x, y);
	if (ret == 0)
	    notifyUpdatedLocked();
	return ret;
    }

    int getcursor(int *mode, int *x, int *y)
    {
	std::lock_guard<std::mutex> guard(lock);
	return framebuffer.getcursor(mode, x, y);
    }

    int setcursor(const unsigned char *bits, int xbits, int ybits,
	    int xorig, int yorig)
    {
	std::lock_guard<std::mutex> guard(lock);
	int ret = framebuffer.setcursor(bits, xbits, ybits, xorig, yorig);
	if (ret == 0)
	    notifyUpdatedLocked();
	return ret;
    }

    int scursor(int mode, int x, int y)
    {
	std::lock_guard<std::mutex> guard(lock);
	int ret = framebuffer.scursor(mode, x, y);
	if (ret == 0)
	    notifyUpdatedLocked();
	return ret;
    }

    int window(int xcenter, int ycenter)
    {
	std::lock_guard<std::mutex> guard(lock);
	int ret = framebuffer.window(xcenter, ycenter);
	if (ret == 0)
	    notifyUpdatedLocked();
	return ret;
    }

    int zoom(int xzoom, int yzoom)
    {
	std::lock_guard<std::mutex> guard(lock);
	int ret = framebuffer.zoom(xzoom, yzoom);
	if (ret == 0)
	    notifyUpdatedLocked();
	return ret;
    }

    int view(int xcenter, int ycenter, int xzoom, int yzoom)
    {
	std::lock_guard<std::mutex> guard(lock);
	int ret = framebuffer.view(xcenter, ycenter, xzoom, yzoom);
	if (ret == 0)
	    notifyUpdatedLocked();
	return ret;
    }

    int getview(int *xcenter, int *ycenter, int *xzoom, int *yzoom)
    {
	std::lock_guard<std::mutex> guard(lock);
	return framebuffer.getview(xcenter, ycenter, xzoom, yzoom);
    }

    int rmap(struct fbserv_colormap *out)
    {
	std::lock_guard<std::mutex> guard(lock);
	if (!out)
	    return -1;

	struct imgstream_fb_colormap imap;
	if (framebuffer.rmap(&imap) != 0)
	    return -1;
	for (int i = 0; i < 256; i++) {
	    out->red[i] = imap.red[i];
	    out->green[i] = imap.green[i];
	    out->blue[i] = imap.blue[i];
	}
	return 0;
    }

    int wmap(const struct fbserv_colormap *in)
    {
	std::lock_guard<std::mutex> guard(lock);
	if (!in)
	    return framebuffer.wmap(NULL);

	struct imgstream_fb_colormap imap;
	for (int i = 0; i < 256; i++) {
	    imap.red[i] = in->red[i];
	    imap.green[i] = in->green[i];
	    imap.blue[i] = in->blue[i];
	}
	return framebuffer.wmap(&imap);
    }

    int flush()
    {
	std::lock_guard<std::mutex> guard(lock);
	int ret = present_on_flush ?
	    framebuffer.flush() :
	    (framebuffer.ensure() == 0 ?
		imgstream_fb_flush(framebuffer.framebuffer()) : -1);
	if (ret == 0)
	    notifyUpdatedLocked();
	return ret;
    }

    int present()
    {
	std::lock_guard<std::mutex> guard(lock);
	bindCaptureProviderLocked();
	return framebuffer.present();
    }

    int apply(ged_draw_obol_framebuffer_operation_t operation,
	    void *userdata, int publish)
    {
	if (!operation)
	    return BRLCAD_ERROR;
	std::lock_guard<std::mutex> guard(lock);
	if (framebuffer.ensure() != 0)
	    return BRLCAD_ERROR;
	int ret = operation(framebuffer.framebuffer(), userdata);
	if (ret != BRLCAD_OK || !publish)
	    return ret;
	bindCaptureProviderLocked();
	/* Publication marks the image stream complete and schedules endpoint
	 * presentation.  Refresh host-side retained image data now, but tolerate a
	 * host that can only finish presentation from its normal update thread. */
	(void)framebuffer.present();
	if (imgstream_fb_flush(framebuffer.framebuffer()) != 0)
	    return BRLCAD_ERROR;
	notifyUpdatedLocked();
	return BRLCAD_OK;
    }

    int captureFramebuffer(unsigned char **pixels, size_t *size,
	    unsigned int *width, unsigned int *height,
	    unsigned int *components)
    {
	if (!pixels || !size || !width || !height || !components)
	    return 0;
	*pixels = NULL;
	*size = 0;
	*width = 0;
	*height = 0;
	*components = 0;

	std::lock_guard<std::mutex> guard(lock);
	BRLObolFramebufferInfo info;
	if (framebuffer.info(&info) != 0 || info.width <= 0 ||
	    info.height <= 0 || info.width > INT_MAX / info.height)
	    return 0;
	const size_t byte_count = (size_t)info.width *
	    (size_t)info.height * 3;
	unsigned char *data = static_cast<unsigned char *>(bu_malloc(
	    byte_count, "Obol framebuffer-plane capture"));
	if (framebuffer.readrect(0, 0, info.width, info.height, data) !=
	    info.width * info.height) {
	    bu_free(data, "Obol framebuffer-plane capture");
	    return 0;
	}
	*pixels = data;
	*size = byte_count;
	*width = (unsigned int)info.width;
	*height = (unsigned int)info.height;
	*components = 3;
	return 1;
    }

    int setComposition(int composition)
    {
	std::lock_guard<std::mutex> guard(lock);
	return setCompositionLocked(composition);
    }

    int setCompositionLocked(int composition)
    {
	if (composition < BRLOBOL_FRAMEBUFFER_COMPOSITION_OFF ||
	    composition > BRLOBOL_FRAMEBUFFER_COMPOSITION_INTERLAY)
	    return BRLCAD_ERROR;
	return framebuffer.setComposition(
	    static_cast<BRLObolFramebufferComposition>(composition)) == 0 ?
	    BRLCAD_OK : BRLCAD_ERROR;
    }

    int poll()
    {
	std::lock_guard<std::mutex> guard(lock);
	return framebuffer.poll();
    }

    int help()
    {
	return 0;
    }

    void detachEndpoint(brlobol_display_endpoint_t *endpoint)
    {
	std::lock_guard<std::mutex> guard(lock);
	if (!endpoint || capture_endpoint != endpoint)
	    return;

	/* The stream image belongs to the host behind this endpoint.  Close it
	 * before the endpoint/controller changes, then let a later ensure bind a
	 * fresh retained image to the replacement host. */
	clearCaptureProviderLocked();
	framebuffer.close();
    }

    void openClient(int slot)
    {
#ifdef _WIN32
	(void)slot;
#else
	if (fbs_client_fd(fbs, slot) <= 0)
	    return;

	std::shared_ptr<GedObolFbservWatcher> watcher =
	    std::make_shared<GedObolFbservWatcher>();
	watcher->slot = slot;

	{
	    std::lock_guard<std::mutex> guard(watcher_lock);
	    watchers[slot] = watcher;
	}

	watcher->thread = std::thread([this, watcher]() {
	    runWatcher(watcher);
	});
#endif
    }

    void closeClient(int slot)
    {
	std::shared_ptr<GedObolFbservWatcher> watcher;
	{
	    std::lock_guard<std::mutex> guard(watcher_lock);
	    std::map<int, std::shared_ptr<GedObolFbservWatcher>>::iterator it =
		watchers.find(slot);
	    if (it == watchers.end())
		return;
	    watcher = it->second;
	    watchers.erase(it);
	}

	watcher->stop = true;
	if (watcher->thread.joinable()) {
	    if (watcher->thread.get_id() == std::this_thread::get_id())
		watcher->thread.detach();
	    else
		watcher->thread.join();
	}
    }

    void stopWatchers()
    {
	std::map<int, std::shared_ptr<GedObolFbservWatcher>> local;
	{
	    std::lock_guard<std::mutex> guard(watcher_lock);
	    local.swap(watchers);
	}

	for (std::map<int, std::shared_ptr<GedObolFbservWatcher>>::iterator it =
		local.begin(); it != local.end(); ++it) {
	    it->second->stop = true;
	}
	for (std::map<int, std::shared_ptr<GedObolFbservWatcher>>::iterator it =
		local.begin(); it != local.end(); ++it) {
	    if (it->second->thread.joinable()) {
		if (it->second->thread.get_id() == std::this_thread::get_id())
		    it->second->thread.detach();
		else
		    it->second->thread.join();
	    }
	}
    }

private:
    void bindCaptureProviderLocked()
    {
	brlobol_display_endpoint_t *endpoint = view_ctx ?
	    ged_view_context_display_endpoint_get(view_ctx) : NULL;
	if (capture_endpoint == endpoint)
	    return;

	clearCaptureProviderLocked();
	if (endpoint &&
	    brlobol_display_endpoint_framebuffer_capture_provider_set(
		endpoint, ged_obol_capture_framebuffer, this))
	    capture_endpoint = endpoint;
    }

    void clearCaptureProviderLocked()
    {
	if (capture_endpoint)
	    (void)brlobol_display_endpoint_framebuffer_capture_provider_set(
		capture_endpoint, NULL, this);
	capture_endpoint = NULL;
    }

    void clearCaptureProvider()
    {
	std::lock_guard<std::mutex> guard(lock);
	clearCaptureProviderLocked();
    }

    void notifyUpdatedLocked()
    {
	if (view_ctx)
	    (void)bv_refresh_request(bv_context_view((struct bv_context *)view_ctx),
		    GED_VIEW_REFRESH_FRAMEBUFFER);
    }

#ifndef _WIN32
    void runWatcher(std::shared_ptr<GedObolFbservWatcher> watcher)
    {
	while (!watcher->stop) {
	    int fd = fbs_client_fd(fbs, watcher->slot);
	    if (fd <= 0)
		break;

	    fd_set read_set;
	    FD_ZERO(&read_set);
	    FD_SET(fd, &read_set);
	    struct timeval timeout;
	    timeout.tv_sec = 0;
	    timeout.tv_usec = 100000;
	    int ret = select(fd + 1, &read_set, NULL, NULL, &timeout);
	    if (watcher->stop)
		break;
	    if (ret < 0) {
		if (errno == EINTR)
		    continue;
		break;
	    }
	    if (ret == 0)
		continue;
	    if (FD_ISSET(fd, &read_set)) {
		void *client_data =
		    fbs_client_handler_data(fbs, watcher->slot);
		if (!client_data)
		    break;
		fbs_existing_client_handler(client_data, 0);
		notifyUpdated();
	    }
	}
    }
#endif

    void notifyUpdated()
    {
	std::lock_guard<std::mutex> guard(lock);
	notifyUpdatedLocked();
    }

    struct ged *gedp = NULL;
    void *view_ctx = NULL;
    struct fbserv_obj *fbs = NULL;
    brlobol_display_endpoint_t *capture_endpoint;
    int present_on_flush = 0;
    BRLObolWindowHost default_host;
    BRLObolFramebufferStream framebuffer;
    std::mutex lock;
    std::mutex watcher_lock;
    std::map<int, std::shared_ptr<GedObolFbservWatcher>> watchers;
};

static GedObolFbservBridge *
bridge_from_ctx(void *ctx)
{
    return static_cast<GedObolFbservBridge *>(ctx);
}

static GedObolFbservBridge *
bridge_from_fbs(struct fbserv_obj *fbsp)
{
    if (!fbsp || fbsp->fbs_fb_ops != &ged_obol_fb_ops)
	return NULL;
    return bridge_from_ctx(fbsp->fbs_fb_ctx);
}

static int
ged_obol_is_listening(struct fbserv_obj *fbsp)
{
    return fbs_listener_fd(fbsp) >= 0;
}

static int
ged_obol_listen_on_port(struct fbserv_obj *, int)
{
    return 0;
}

static void
ged_obol_server_noop(struct fbserv_obj *)
{
}

static void
ged_obol_open_client(struct fbserv_obj *fbsp, int slot, void *)
{
    GedObolFbservBridge *bridge = bridge_from_fbs(fbsp);
    if (bridge)
	bridge->openClient(slot);
}

static void
ged_obol_close_client(struct fbserv_obj *fbsp, int slot)
{
    GedObolFbservBridge *bridge = bridge_from_fbs(fbsp);
    if (bridge)
	bridge->closeClient(slot);
}

static const struct fbserv_transport_ops ged_obol_transport_ops = {
    ged_obol_is_listening,
    ged_obol_listen_on_port,
    ged_obol_server_noop,
    ged_obol_server_noop,
    ged_obol_open_client,
    ged_obol_close_client,
    ged_obol_open_client,
    ged_obol_close_client
};

static int ged_obol_fb_info(void *ctx, struct fbserv_fb_info *info) { return bridge_from_ctx(ctx)->info(info); }
static int ged_obol_fb_clear(void *ctx, const unsigned char rgb[3]) { return bridge_from_ctx(ctx)->clear(rgb); }
static ssize_t ged_obol_fb_read(void *ctx, int x, int y, unsigned char *rgb, size_t count) { return bridge_from_ctx(ctx)->read(x, y, rgb, count); }
static ssize_t ged_obol_fb_write(void *ctx, int x, int y, const unsigned char *rgb, size_t count) { return bridge_from_ctx(ctx)->write(x, y, rgb, count); }
static int ged_obol_fb_readrect(void *ctx, int xmin, int ymin, int width, int height, unsigned char *rgb) { return bridge_from_ctx(ctx)->readrect(xmin, ymin, width, height, rgb); }
static int ged_obol_fb_writerect(void *ctx, int xmin, int ymin, int width, int height, const unsigned char *rgb) { return bridge_from_ctx(ctx)->writerect(xmin, ymin, width, height, rgb); }
static int ged_obol_fb_bwreadrect(void *ctx, int xmin, int ymin, int width, int height, unsigned char *bw) { return bridge_from_ctx(ctx)->bwreadrect(xmin, ymin, width, height, bw); }
static int ged_obol_fb_bwwriterect(void *ctx, int xmin, int ymin, int width, int height, const unsigned char *bw) { return bridge_from_ctx(ctx)->bwwriterect(xmin, ymin, width, height, bw); }
static int ged_obol_fb_cursor(void *ctx, int mode, int x, int y) { return bridge_from_ctx(ctx)->cursor(mode, x, y); }
static int ged_obol_fb_getcursor(void *ctx, int *mode, int *x, int *y) { return bridge_from_ctx(ctx)->getcursor(mode, x, y); }
static int ged_obol_fb_setcursor(void *ctx, const unsigned char *bits, int xbits, int ybits, int xorig, int yorig) { return bridge_from_ctx(ctx)->setcursor(bits, xbits, ybits, xorig, yorig); }
static int ged_obol_fb_scursor(void *ctx, int mode, int x, int y) { return bridge_from_ctx(ctx)->scursor(mode, x, y); }
static int ged_obol_fb_window(void *ctx, int xcenter, int ycenter) { return bridge_from_ctx(ctx)->window(xcenter, ycenter); }
static int ged_obol_fb_zoom(void *ctx, int xzoom, int yzoom) { return bridge_from_ctx(ctx)->zoom(xzoom, yzoom); }
static int ged_obol_fb_view(void *ctx, int xcenter, int ycenter, int xzoom, int yzoom) { return bridge_from_ctx(ctx)->view(xcenter, ycenter, xzoom, yzoom); }
static int ged_obol_fb_getview(void *ctx, int *xcenter, int *ycenter, int *xzoom, int *yzoom) { return bridge_from_ctx(ctx)->getview(xcenter, ycenter, xzoom, yzoom); }
static int ged_obol_fb_rmap(void *ctx, struct fbserv_colormap *cmap) { return bridge_from_ctx(ctx)->rmap(cmap); }
static int ged_obol_fb_wmap(void *ctx, const struct fbserv_colormap *cmap) { return bridge_from_ctx(ctx)->wmap(cmap); }
static int ged_obol_fb_flush(void *ctx) { return bridge_from_ctx(ctx)->flush(); }
static int ged_obol_fb_poll(void *ctx) { return bridge_from_ctx(ctx)->poll(); }
static int ged_obol_fb_help(void *ctx) { return bridge_from_ctx(ctx)->help(); }

static int
ged_obol_capture_framebuffer(void *ctx, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components)
{
    GedObolFbservBridge *bridge = bridge_from_ctx(ctx);
    return bridge ? bridge->captureFramebuffer(pixels, size, width, height,
	components) : 0;
}

} // namespace

static int
ged_obol_fbserv_configure_for_view(struct ged *gedp,
	void *view_ctx,
	BRLObolWindowHost *window_host,
	int width,
	int height,
	int present_on_flush,
	int install_default_handlers)
{
    if (!gedp || !gedp->ged_fbs || !view_ctx)
	return BRLCAD_ERROR;

    /* Most toolkit factories keep their instance opaque.  A factory that owns
     * a BRLObolWindowHost may explicitly expose it for retained framebuffer
     * attachment; otherwise the bridge uses its controller-owned host. */
    if (!window_host) {
	brlobol_display_endpoint_t *endpoint =
	    ged_view_context_display_endpoint_get(view_ctx);
	if (endpoint) {
	    window_host = static_cast<BRLObolWindowHost *>(
		brlobol_display_endpoint_framebuffer_window_host(endpoint));
	    /* A live endpoint must publish flushes to its retained image source.
	     * Headless/direct callers retain the requested deferred behavior. */
	    present_on_flush = 1;
	}
    }

    void *controller =
	ged_draw_obol_controller_ensure_opaque_for_view(view_ctx, 1);
    if (!controller)
	return BRLCAD_ERROR;

    struct fbserv_obj *fbs = gedp->ged_fbs;
    GedObolFbservBridge *bridge = bridge_from_fbs(fbs);
    if (!bridge) {
	if (fbs->fbs_fb_ops || fbs->fbs_fb_ctx)
	    return BRLCAD_ERROR;

	bridge = new GedObolFbservBridge(fbs);
	if (fbs_set_backend(fbs, &ged_obol_fb_ops, bridge) != 0) {
	    delete bridge;
	    return BRLCAD_ERROR;
	}
    }


    /* A direct host install may supply its transport immediately afterward
     * (Qt/Tk), but the bridge must still be independently closable during
     * that interval.  Preserve any application transport already installed. */
    if (install_default_handlers || !fbs_can_close(fbs))
	fbs_set_transport(fbs, &ged_obol_transport_ops);

    return bridge->configure(gedp, view_ctx,
	    static_cast<BRLObolViewController *>(controller),
	    window_host, width, height, present_on_flush);
}

extern "C" int
ged_obol_fbserv_ensure_for_view(struct ged *gedp, void *view_ctx)
{
    return ged_draw_obol_framebuffer_backend_ensure_for_view(gedp, view_ctx);
}

extern "C" GED_EXPORT int
ged_draw_obol_framebuffer_backend_ensure_for_view(struct ged *gedp,
	void *view_ctx)
{
    return ged_obol_fbserv_configure_for_view(gedp, view_ctx, NULL, 0, 0,
	0, 1);
}

extern "C" GED_EXPORT int
ged_draw_obol_framebuffer_apply_for_view(struct ged *gedp,
	void *view_ctx,
	ged_draw_obol_framebuffer_operation_t operation,
	void *userdata,
	int publish)
{
    if (!gedp || !operation)
	return BRLCAD_ERROR;
    if (!view_ctx)
	view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx ||
	ged_draw_obol_framebuffer_backend_ensure_for_view(gedp, view_ctx) !=
	    BRLCAD_OK)
	return BRLCAD_ERROR;

    GedObolFbservBridge *bridge = bridge_from_fbs(gedp->ged_fbs);
    return bridge ? bridge->apply(operation, userdata, publish) :
	BRLCAD_ERROR;
}

extern "C" GED_EXPORT int
ged_draw_obol_framebuffer_backend_install_for_view(struct ged *gedp,
	void *view_ctx,
	void *window_host,
	int width,
	int height,
	int present_on_flush)
{
    return ged_obol_fbserv_configure_for_view(gedp, view_ctx,
	    static_cast<BRLObolWindowHost *>(window_host), width, height,
	    present_on_flush, 0);
}

extern "C" void
ged_obol_fbserv_release(struct ged *gedp)
{
    if (!gedp || !gedp->ged_fbs)
	return;

    struct fbserv_obj *fbs = gedp->ged_fbs;
    GedObolFbservBridge *bridge = bridge_from_fbs(fbs);
    if (!bridge)
	return;

    fbs_close(fbs);
    fbs_clear_backend(fbs);
    delete bridge;

    fbs_clear_transport(fbs);
}

extern "C" int
ged_obol_fbserv_composition_set(struct ged *gedp, int mode)
{
    if (!gedp || !gedp->ged_fbs)
	return BRLCAD_ERROR;
    GedObolFbservBridge *bridge = bridge_from_fbs(gedp->ged_fbs);
    /* Without a live Obol framebuffer bridge this remains passive desired
     * view state.  Once a bridge exists its host must accept the update. */
    return bridge ? bridge->setComposition(mode) : BRLCAD_OK;
}

extern "C" GED_EXPORT int
ged_draw_obol_view_display_image(struct ged *gedp,
				 void *view_ctx,
				 unsigned char **image,
				 int flip,
				 int alpha)
{
    if (!image)
	return -1;
    *image = NULL;

    if (!gedp)
	return -1;
    if (!view_ctx)
	view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx)
	return -1;

    BRLObolViewController *controller =
	static_cast<BRLObolViewController *>(
	    ged_draw_obol_controller_opaque_for_view(view_ctx));
    if (!controller)
	return 0;

    int width = 0;
    int height = 0;
    const struct bv *view =
	bv_context_view_const((const struct bv_context *)view_ctx);
    if (view) {
	int vw = bv_width_get(view);
	int vh = bv_height_get(view);
	if (width <= 0 && vw > 0)
	    width = vw;
	if (height <= 0 && vh > 0)
	    height = vh;
    }
    if (width <= 0)
	width = 512;
    if (height <= 0)
	height = 512;

    controller->setViewportSize((unsigned int)width,
	(unsigned int)height);
    if (!controller->syncCameraFromViewContext(view_ctx))
	return -1;

    BRLObolProgressiveStatus progressiveStatus;
    int ret = controller->renderToImage(image, flip, alpha, NULL,
	NULL, &progressiveStatus);
    if (progressiveStatus.hasMore && view)
	(void)bv_refresh_request(bv_context_view((struct bv_context *)view_ctx),
		GED_VIEW_REFRESH_VIEW);
    return ret == BRLCAD_OK ? 1 : -1;
}

extern "C" int
ged_obol_fbserv_present(struct ged *gedp)
{
    if (!gedp || !gedp->ged_fbs)
	return BRLCAD_ERROR;

    GedObolFbservBridge *bridge = bridge_from_fbs(gedp->ged_fbs);
    if (!bridge)
	return BRLCAD_ERROR;

    return bridge->present() == 0 ? BRLCAD_OK : BRLCAD_ERROR;
}

extern "C" int
ged_draw_obol_framebuffer_present(struct ged *gedp)
{
    return ged_obol_fbserv_present(gedp);
}

extern "C" void
ged_draw_obol_framebuffer_release(struct ged *gedp)
{
    ged_obol_fbserv_release(gedp);
}

extern "C" GED_EXPORT void
ged_draw_obol_framebuffer_endpoint_detach(struct ged *gedp,
	brlobol_display_endpoint_t *endpoint)
{
    if (!gedp || !gedp->ged_fbs || !endpoint)
	return;

    GedObolFbservBridge *bridge = bridge_from_fbs(gedp->ged_fbs);
    if (bridge)
	bridge->detachEndpoint(endpoint);
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
