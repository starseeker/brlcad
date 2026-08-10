/*              D I S P L A Y _ S E S S I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BDisplaySession.h"
#include "BObol/BViewController.h"
#include "BObol/BWindowHost.h"
#include "bu/str.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <mutex>
#include <new>
#include <thread>
#include <vector>

struct bobol_display_session {
    bobol_display_endpoint_t *endpoint = NULL;
    const bobol_display_provider_t *provider = NULL;
    void *provider_instance = NULL;
    BObolWindowHost attachment_host;
    imgstream_fb_t *framebuffer = NULL;
    std::atomic<bool> present_pending = false;
};

static std::mutex provider_lock;
static std::vector<const bobol_display_provider_t *> providers;

static int
session_present(bobol_display_session_t *session, imgstream_fb_t *fb,
	bool force)
{
    if (!session || !fb)
	return -1;
    if (!force && !session->present_pending.exchange(false,
		std::memory_order_acq_rel))
	return 0;
    session->present_pending.store(false, std::memory_order_release);
    if (session->attachment_host.flushFramebuffer(fb) != 0)
	return -1;
    return bobol_display_endpoint_request_frame(session->endpoint,
	"framebuffer-present") ? 0 : -1;
}

static int
session_request_present(bobol_display_session_t *session)
{
    if (!session)
	return -1;
    bool expected = false;
    if (!session->present_pending.compare_exchange_strong(expected, true,
	std::memory_order_acq_rel))
	return 0;
    if (bobol_display_endpoint_request_frame(session->endpoint,
	"framebuffer-dirty"))
	return 0;
    session->present_pending.store(false, std::memory_order_release);
    return -1;
}

static int
session_open(imgstream_fb_t *fb, const imgstream_fb_spec_info_t *info,
	void *data)
{
    bobol_display_session_t *session =
	static_cast<bobol_display_session_t *>(data);
    BObolViewController *controller = session ?
	static_cast<BObolViewController *>(
	    bobol_display_endpoint_controller(session->endpoint)) : NULL;
    if (!session || !controller || !fb || !info)
	return -1;
    session->attachment_host.attachController(controller, FALSE);
    return session->attachment_host.openFramebuffer(fb, info);
}

static void
session_close(imgstream_fb_t *fb, void *data)
{
    bobol_display_session_t *session =
	static_cast<bobol_display_session_t *>(data);
    if (session && fb)
	session->attachment_host.closeFramebuffer(fb);
}

static int
session_flush(imgstream_fb_t *fb, void *data)
{
    if (!fb)
	return -1;
    /* Writers may run on tracing or transport threads.  They only request a
     * frame here; the caller-owned poll path refreshes Coin image nodes. */
    return session_request_present(static_cast<bobol_display_session_t *>(data));
}

static int
session_reset(imgstream_fb_t *fb, void *data)
{
    bobol_display_session_t *session =
	static_cast<bobol_display_session_t *>(data);
    if (!session || session->attachment_host.resetFramebuffer(fb) != 0)
	return -1;
    return bobol_display_endpoint_request_frame(session->endpoint,
	"framebuffer-reset") ? 0 : -1;
}

static int
session_viewport(imgstream_fb_t *fb, int left, int top, int right, int bottom,
	void *data)
{
    bobol_display_session_t *session =
	static_cast<bobol_display_session_t *>(data);
    return session && session->attachment_host.setFramebufferViewport(fb,
	left, top, right, bottom) == 0 ? 0 : -1;
}

static int
session_view(imgstream_fb_t *fb, const imgstream_fb_view_t *view, void *data)
{
    bobol_display_session_t *session =
	static_cast<bobol_display_session_t *>(data);
    return session && session->attachment_host.setFramebufferView(fb, view) == 0 ?
	0 : -1;
}

static int
session_cursor(imgstream_fb_t *fb, const imgstream_fb_cursor_t *cursor,
	void *data)
{
    bobol_display_session_t *session =
	static_cast<bobol_display_session_t *>(data);
    return session && session->attachment_host.setFramebufferCursor(fb, cursor) == 0 ?
	0 : -1;
}

static int
session_scursor(imgstream_fb_t *fb, int mode, int x, int y, void *data)
{
    bobol_display_session_t *session =
	static_cast<bobol_display_session_t *>(data);
    return session && session->attachment_host.setFramebufferScreenCursor(fb,
	mode, x, y) == 0 ? 0 : -1;
}

static int
session_setcursor(imgstream_fb_t *fb, const unsigned char *bits, int xbits,
	int ybits, int xorig, int yorig, void *data)
{
    bobol_display_session_t *session =
	static_cast<bobol_display_session_t *>(data);
    return session && session->attachment_host.setFramebufferCursorShape(fb,
	bits, xbits, ybits, xorig, yorig) == 0 ? 0 : -1;
}

static int
session_poll(imgstream_fb_t *fb, void *data)
{
    bobol_display_session_t *session =
	static_cast<bobol_display_session_t *>(data);
    if (!session)
	return -1;
    /* Polling supplies progressive output for writers that do not explicitly
     * flush every scan line, while keeping toolkit work on the caller thread. */
    if (session_present(session, fb, false) != 0)
	return -1;
    return session->provider && session->provider->poll ?
	session->provider->poll(session->provider_instance,
	    session->provider->user_data) : 0;
}

static long
session_poll_rate(const imgstream_fb_t *UNUSED(fb), void *data)
{
    const bobol_display_session_t *session =
	static_cast<const bobol_display_session_t *>(data);
    if (!session || !session->provider || !session->provider->poll_rate)
	return 16000;
    long rate = session->provider->poll_rate(session->provider_instance,
	session->provider->user_data);
    return rate > 0 ? rate : 16000;
}

extern "C" int
bobol_display_provider_register(const bobol_display_provider_t *provider)
{
    if (!provider || provider->abi_version != BOBOL_DISPLAY_PROVIDER_ABI_VERSION ||
	provider->struct_size < sizeof(bobol_display_provider_t) ||
	!provider->name || !provider->name[0] || !provider->open)
	return 0;
    std::lock_guard<std::mutex> guard(provider_lock);
    for (const bobol_display_provider_t *candidate : providers) {
	if (bu_strcmp(candidate->name, provider->name) == 0)
	    return candidate == provider;
    }
    providers.push_back(provider);
    std::stable_sort(providers.begin(), providers.end(),
	[](const bobol_display_provider_t *a,
	   const bobol_display_provider_t *b) {
	    return a->priority > b->priority;
	});
    return 1;
}

extern "C" bobol_display_session_t *
bobol_display_session_open(const char *spec, size_t width, size_t height,
	const char *title)
{
    imgstream_fb_spec_info_t info = {};
    if (!spec || imgstream_fb_spec_info(spec, &info) != 0 ||
	info.kind != IMGSTREAM_FB_SPEC_DISPLAY)
	return NULL;
    /* imgstream accepts zero dimensions as its historical 512-pixel default.
     * A native top-level host needs concrete dimensions before the framebuffer
     * is created, so apply that compatibility default at this one boundary. */
    width = width ? width : 512u;
    height = height ? height : 512u;

    std::vector<const bobol_display_provider_t *> candidates;
    {
	std::lock_guard<std::mutex> guard(provider_lock);
	candidates = providers;
    }
    for (const bobol_display_provider_t *provider : candidates) {
	bobol_display_endpoint_t *endpoint =
	    bobol_display_endpoint_create(NULL, BOBOL_ENDPOINT_OWN_CONTROLLER);
	if (!endpoint)
	    continue;
	void *provider_instance = NULL;
	if (!provider->open(endpoint, &info, width, height, title,
		&provider_instance, provider->user_data)) {
	    bobol_display_endpoint_destroy(endpoint);
	    continue;
	}

	bobol_display_session_t *session = new (std::nothrow)
	    bobol_display_session_t;
	if (!session) {
	    bobol_display_endpoint_destroy(endpoint);
	    if (provider->close)
		provider->close(provider_instance, provider->user_data);
	    continue;
	}
	session->endpoint = endpoint;
	session->provider = provider;
	session->provider_instance = provider_instance;
	struct imgstream_fb_display_host host = {};
	host.open = session_open;
	host.close = session_close;
	host.flush = session_flush;
	host.reset = session_reset;
	host.viewport = session_viewport;
	host.view = session_view;
	host.cursor = session_cursor;
	host.scursor = session_scursor;
	host.setcursor = session_setcursor;
	host.poll = session_poll;
	host.poll_rate = session_poll_rate;
	session->framebuffer = imgstream_fb_open_display(spec, width, height,
	    &host, session);
	if (session->framebuffer)
	    return session;
	bobol_display_endpoint_destroy(endpoint);
	if (provider->close)
	    provider->close(provider_instance, provider->user_data);
	delete session;
    }
    return NULL;
}

extern "C" void
bobol_display_session_close(bobol_display_session_t *session)
{
    if (!session)
	return;
    if (session->framebuffer)
	imgstream_fb_close(session->framebuffer);
    session->framebuffer = NULL;
    bobol_display_endpoint_destroy(session->endpoint);
    if (session->provider && session->provider->close)
	session->provider->close(session->provider_instance,
	    session->provider->user_data);
    delete session;
}

extern "C" imgstream_fb_t *
bobol_display_session_framebuffer(const bobol_display_session_t *session)
{
    return session ? session->framebuffer : NULL;
}

extern "C" int
bobol_display_session_present(bobol_display_session_t *session)
{
    return session && session->framebuffer ?
	session_present(session, session->framebuffer, true) : -1;
}

extern "C" int
bobol_display_session_poll(bobol_display_session_t *session)
{
    return session && session->framebuffer ? session_poll(session->framebuffer,
	session) : -1;
}

extern "C" long
bobol_display_session_poll_rate(const bobol_display_session_t *session)
{
    return session ? session_poll_rate(session->framebuffer,
	const_cast<bobol_display_session_t *>(session)) : 0;
}

extern "C" int
bobol_display_session_run(bobol_display_session_t *session,
	bobol_display_session_task_t task, void *task_data)
{
    if (!session || !session->framebuffer || !task)
	return -1;

    std::atomic<bool> complete(false);
    int task_result = -1;
    std::thread worker([&]() {
	task_result = task(task_data);
	complete.store(true, std::memory_order_release);
    });

    while (!complete.load(std::memory_order_acquire)) {
	(void)bobol_display_session_poll(session);
	const long rate = std::max(1000L,
	    bobol_display_session_poll_rate(session));
	std::this_thread::sleep_for(std::chrono::microseconds(rate));
    }
    worker.join();
    (void)bobol_display_session_poll(session);
    return task_result;
}
