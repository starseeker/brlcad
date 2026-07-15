/*                    F B _ C O M P A T . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file libimgstream/tests/fb_compat.c */

#include "common.h"

#include "bu/str.h"

#include <stdio.h>
#include <string.h>

#include "bu/app.h"
#include "bu/file.h"
#include "bu/log.h"
#include "imgstream.h"


#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	return 1; \
    } \
} while (0)


struct display_host_test_state {
    int open_count;
    int close_count;
    int flush_count;
    int reset_count;
    int viewport_count;
    int view_count;
    int cursor_count;
    int scursor_count;
    int setcursor_count;
    int poll_count;
    int poll_rate_count;
    int open_return;
    int poll_return;
    long poll_rate_return;
    enum imgstream_fb_display_kind open_display;
    size_t open_width;
    size_t open_height;
    imgstream_fb_view_t last_view;
    imgstream_fb_cursor_t last_cursor;
    int viewport[4];
    int scursor[3];
    int setcursor_args[4];
};


static int
display_host_open(imgstream_fb_t *fb, const imgstream_fb_spec_info_t *info, void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (!state || !fb || !info)
	return -1;

    state->open_count++;
    state->open_display = info->display;
    state->open_width = imgstream_fb_width(fb);
    state->open_height = imgstream_fb_height(fb);
    return state->open_return;
}


static void
display_host_close(imgstream_fb_t *UNUSED(fb), void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (state)
	state->close_count++;
}


static int
display_host_flush(imgstream_fb_t *UNUSED(fb), void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (!state)
	return -1;
    state->flush_count++;
    return 0;
}


static int
display_host_reset(imgstream_fb_t *UNUSED(fb), void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (!state)
	return -1;
    state->reset_count++;
    return 0;
}


static int
display_host_viewport(imgstream_fb_t *UNUSED(fb), int left, int top, int right, int bottom, void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (!state)
	return -1;
    state->viewport_count++;
    state->viewport[0] = left;
    state->viewport[1] = top;
    state->viewport[2] = right;
    state->viewport[3] = bottom;
    return 0;
}


static int
display_host_view(imgstream_fb_t *UNUSED(fb), const imgstream_fb_view_t *view, void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (!state || !view)
	return -1;
    state->view_count++;
    state->last_view = *view;
    return 0;
}


static int
display_host_cursor(imgstream_fb_t *UNUSED(fb), const imgstream_fb_cursor_t *cursor, void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (!state || !cursor)
	return -1;
    state->cursor_count++;
    state->last_cursor = *cursor;
    return 0;
}


static int
display_host_scursor(imgstream_fb_t *UNUSED(fb), int mode, int x, int y, void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (!state)
	return -1;
    state->scursor_count++;
    state->scursor[0] = mode;
    state->scursor[1] = x;
    state->scursor[2] = y;
    return 0;
}


static int
display_host_setcursor(imgstream_fb_t *UNUSED(fb), const unsigned char *UNUSED(bits), int xbits, int ybits, int xorig, int yorig, void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (!state)
	return -1;
    state->setcursor_count++;
    state->setcursor_args[0] = xbits;
    state->setcursor_args[1] = ybits;
    state->setcursor_args[2] = xorig;
    state->setcursor_args[3] = yorig;
    return 0;
}


static int
display_host_poll(imgstream_fb_t *UNUSED(fb), void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (!state)
	return -1;
    state->poll_count++;
    return state->poll_return;
}


static long
display_host_poll_rate(const imgstream_fb_t *UNUSED(fb), void *data)
{
    struct display_host_test_state *state = (struct display_host_test_state *)data;
    if (!state)
	return 0;
    state->poll_rate_count++;
    return state->poll_rate_return;
}


static int
test_spec_classification(void)
{
    imgstream_fb_spec_info_t info;

    CHECK(imgstream_fb_spec_kind(NULL) == IMGSTREAM_FB_SPEC_MEMORY, "null spec maps to memory");
    CHECK(imgstream_fb_spec_kind("") == IMGSTREAM_FB_SPEC_MEMORY, "empty spec maps to memory");
    CHECK(imgstream_fb_spec_kind("/dev/mem") == IMGSTREAM_FB_SPEC_MEMORY, "memory spec maps to memory");
    CHECK(imgstream_fb_spec_supported("/dev/mem") == 1, "memory spec is currently supported");
    CHECK(imgstream_fb_spec_info(NULL, &info) == 0 && info.kind == IMGSTREAM_FB_SPEC_MEMORY,
	    "null spec info maps to memory");

    CHECK(imgstream_fb_spec_kind("/dev/disk") == IMGSTREAM_FB_SPEC_FILE, "disk spec maps to file");
    CHECK(imgstream_fb_spec_supported("/dev/disk") == 0, "bare disk spec needs an explicit file path");
    CHECK(imgstream_fb_spec_kind("out.pix") == IMGSTREAM_FB_SPEC_FILE, "plain path maps to file");
    CHECK(imgstream_fb_spec_supported("out.pix") == 1, "plain path file spec is supported");
    CHECK(imgstream_fb_spec_kind("/dev/disk:out.pix") == IMGSTREAM_FB_SPEC_FILE, "explicit disk path maps to file");
    CHECK(imgstream_fb_spec_supported("/dev/disk:out.pix") == 1, "explicit disk path file spec is supported");
    CHECK(imgstream_fb_spec_info("/dev/disk:out.pix", &info) == 0 &&
	    info.kind == IMGSTREAM_FB_SPEC_FILE &&
	    info.target && info.target_len == 7 &&
	    memcmp(info.target, "out.pix", info.target_len) == 0,
	    "explicit disk path spec info records file target");
    CHECK(imgstream_fb_spec_kind("/dev/stack") == IMGSTREAM_FB_SPEC_FANOUT, "stack spec maps to fanout");
    CHECK(imgstream_fb_spec_info("/dev/stack /dev/mem; /dev/null", &info) == 0 &&
	    info.kind == IMGSTREAM_FB_SPEC_FANOUT &&
	    info.target && info.target_len == strlen("/dev/mem; /dev/null") &&
	    memcmp(info.target, "/dev/mem; /dev/null", info.target_len) == 0,
	    "stack spec info records child target list");
    CHECK(imgstream_fb_spec_kind("/dev/remote") == IMGSTREAM_FB_SPEC_REMOTE, "remote device spec maps to remote");
    CHECK(imgstream_fb_spec_kind("localhost:5555") == IMGSTREAM_FB_SPEC_REMOTE, "host:port spec maps to remote");
    CHECK(imgstream_fb_spec_kind("ipc:/tmp/fb") == IMGSTREAM_FB_SPEC_REMOTE, "ipc spec maps to remote");
    CHECK(imgstream_fb_spec_kind("0") == IMGSTREAM_FB_SPEC_REMOTE, "numeric port spec maps to remote");
    CHECK(imgstream_fb_spec_supported("0") == 1, "numeric remote spec is stream supported");
    CHECK(imgstream_fb_spec_info("0", &info) == 0 &&
	    info.kind == IMGSTREAM_FB_SPEC_REMOTE &&
	    info.remote == IMGSTREAM_FB_REMOTE_TCP &&
	    info.host && info.host_len == 9 &&
	    memcmp(info.host, "localhost", info.host_len) == 0 &&
	    info.port == 5559,
	    "numeric remote spec info maps through legacy port offset");
    CHECK(imgstream_fb_spec_info("localhost:5555", &info) == 0 &&
	    info.kind == IMGSTREAM_FB_SPEC_REMOTE &&
	    info.remote == IMGSTREAM_FB_REMOTE_TCP &&
	    info.host && info.host_len == 9 &&
	    memcmp(info.host, "localhost", info.host_len) == 0 &&
	    info.port == 5555 &&
	    info.device_len == 0,
	    "host:port spec info records TCP target");
    CHECK(imgstream_fb_spec_info("host:0:/dev/ogl", &info) == 0 &&
	    info.kind == IMGSTREAM_FB_SPEC_REMOTE &&
	    info.remote == IMGSTREAM_FB_REMOTE_TCP &&
	    info.host && info.host_len == 4 &&
	    memcmp(info.host, "host", info.host_len) == 0 &&
	    info.port == 5559 &&
	    info.device && info.device_len == 8 &&
	    memcmp(info.device, "/dev/ogl", info.device_len) == 0,
	    "host:port:device spec info records remote device");
    CHECK(imgstream_fb_spec_info("unix:/dev/X", &info) == 0 &&
	    info.kind == IMGSTREAM_FB_SPEC_REMOTE &&
	    info.remote == IMGSTREAM_FB_REMOTE_TCP &&
	    info.host && info.host_len == 9 &&
	    memcmp(info.host, "localhost", info.host_len) == 0 &&
	    info.port == 5558 &&
	    info.device && info.device_len == 6 &&
	    memcmp(info.device, "/dev/X", info.device_len) == 0,
	    "unix host alias maps to localhost");
    CHECK(imgstream_fb_spec_info("ipc:/tmp/fb", &info) == 0 &&
	    info.kind == IMGSTREAM_FB_SPEC_REMOTE &&
	    info.remote == IMGSTREAM_FB_REMOTE_IPC &&
	    info.target && info.target_len == 7 &&
	    memcmp(info.target, "/tmp/fb", info.target_len) == 0,
	    "ipc remote spec info records IPC target");
    CHECK(imgstream_fb_spec_info("tcp:localhost:5555", &info) == 0 &&
	    info.kind == IMGSTREAM_FB_SPEC_REMOTE &&
	    info.remote == IMGSTREAM_FB_REMOTE_TCP &&
	    info.port == 5555,
	    "tcp-prefixed remote spec info records TCP target");
    CHECK(imgstream_fb_spec_kind("/dev/qtgl") == IMGSTREAM_FB_SPEC_DISPLAY, "display spec maps to display");
    CHECK(imgstream_fb_spec_info("/dev/qtgl", &info) == 0 &&
	    info.kind == IMGSTREAM_FB_SPEC_DISPLAY &&
	    info.display == IMGSTREAM_FB_DISPLAY_QTGL &&
	    info.target && info.target_len == 4 &&
	    memcmp(info.target, "qtgl", info.target_len) == 0,
	    "display spec info records display kind");
    CHECK(imgstream_fb_spec_kind("/dev/null") == IMGSTREAM_FB_SPEC_DIAGNOSTIC, "null spec maps to diagnostic");
    CHECK(imgstream_fb_spec_info("/dev/null", &info) == 0 &&
	    info.kind == IMGSTREAM_FB_SPEC_DIAGNOSTIC &&
	    info.diagnostic == IMGSTREAM_FB_DIAGNOSTIC_NULL,
	    "diagnostic spec info records diagnostic kind");
    CHECK(imgstream_fb_spec_kind("/dev/debug4") == IMGSTREAM_FB_SPEC_DIAGNOSTIC, "debug spec with flags maps to diagnostic");
    CHECK(imgstream_fb_spec_kind("/dev/debugger") == IMGSTREAM_FB_SPEC_UNSUPPORTED, "debug prefix without numeric flags is unsupported");
    CHECK(imgstream_fb_spec_kind("/dev/txt") == IMGSTREAM_FB_SPEC_DIAGNOSTIC, "txt spec maps to diagnostic");
    CHECK(imgstream_fb_spec_supported("/dev/null") == 1, "null diagnostic spec is supported");
    CHECK(imgstream_fb_spec_supported("/dev/debug4") == 1, "debug diagnostic spec is supported");
    CHECK(imgstream_fb_spec_kind("/dev/not-real") == IMGSTREAM_FB_SPEC_UNSUPPORTED, "unknown device spec is unsupported");
    CHECK(imgstream_fb_spec_supported("/dev/qtgl") == 0, "display spec is not core-stream supported");
    CHECK(imgstream_fb_open("/dev/qtgl", 4, 4) == NULL, "display spec open fails explicitly for now");

    CHECK(bu_strcmp(imgstream_fb_spec_kind_name(IMGSTREAM_FB_SPEC_MEMORY), "memory") == 0, "memory kind name");
    CHECK(bu_strcmp(imgstream_fb_spec_kind_name(IMGSTREAM_FB_SPEC_DISPLAY), "display") == 0, "display kind name");

    return 0;
}


static int
test_display_host_compat_stream(void)
{
    struct display_host_test_state state;
    struct display_host_test_state state2;
    struct imgstream_fb_display_host host;
    memset(&state, 0, sizeof(state));
    memset(&state2, 0, sizeof(state2));
    memset(&host, 0, sizeof(host));

    CHECK(imgstream_fb_open_display("/dev/qtgl", 4, 4, NULL, &state) == NULL,
	    "null display host rejected");
    CHECK(imgstream_fb_open_display("/dev/qtgl", 4, 4, &host, &state) == NULL,
	    "display host without open callback rejected");

    host.open = display_host_open;
    host.close = display_host_close;
    host.flush = display_host_flush;
    host.reset = display_host_reset;
    host.viewport = display_host_viewport;
    host.view = display_host_view;
    host.cursor = display_host_cursor;
    host.scursor = display_host_scursor;
    host.setcursor = display_host_setcursor;
    host.poll = display_host_poll;
    host.poll_rate = display_host_poll_rate;
    state.poll_return = 12;
    state.poll_rate_return = 30;

    CHECK(imgstream_fb_spec_supported("/dev/qtgl") == 0,
	    "display spec requires an explicit host");

    imgstream_fb_t *fb = imgstream_fb_open_display("/dev/qtgl", 5, 4, &host, &state);
    CHECK(fb != NULL, "display compatibility stream opens through host");
    CHECK(state.open_count == 1 && state.open_display == IMGSTREAM_FB_DISPLAY_QTGL,
	    "display host open callback saw qtgl display spec");
    CHECK(state.open_width == 5 && state.open_height == 4, "display host open callback saw requested size");
    CHECK(imgstream_fb_stream(fb) != NULL, "display compatibility stream owns image stream");

    unsigned char clear_rgb[3] = {21, 22, 23};
    unsigned char pixel[3] = {0, 0, 0};
    CHECK(imgstream_fb_clear(fb, clear_rgb) == 0, "display compatibility clear updates stream");
    CHECK(imgstream_fb_read(fb, 0, 0, pixel, 1) == 1, "display compatibility read accepted");
    CHECK(pixel[0] == 21 && pixel[1] == 22 && pixel[2] == 23, "display compatibility stream stores pixels");

    CHECK(imgstream_fb_view(fb, 3, 2, 4, 5) == 0, "display view update accepted");
    CHECK(state.view_count == 1 && state.last_view.xcenter == 3 &&
	    state.last_view.ycenter == 2 && state.last_view.xzoom == 4 &&
	    state.last_view.yzoom == 5, "display host view callback received view state");
    CHECK(imgstream_fb_cursor(fb, 1, 4, 3) == 0, "display cursor update accepted");
    CHECK(state.cursor_count == 1 && state.last_cursor.mode == 1 &&
	    state.last_cursor.x == 4 && state.last_cursor.y == 3,
	    "display host cursor callback received cursor state");
    CHECK(imgstream_fb_scursor(fb, 2, 7, 8) == 0, "display screen cursor update accepted");
    CHECK(state.scursor_count == 1 && state.scursor[0] == 2 &&
	    state.scursor[1] == 7 && state.scursor[2] == 8,
	    "display host screen cursor callback received arguments");
    CHECK(imgstream_fb_setcursor(fb, NULL, 16, 15, 8, 7) == 0, "display setcursor accepted");
    CHECK(state.setcursor_count == 1 && state.setcursor_args[0] == 16 &&
	    state.setcursor_args[1] == 15 && state.setcursor_args[2] == 8 &&
	    state.setcursor_args[3] == 7,
	    "display host setcursor callback received arguments");
    CHECK(imgstream_fb_viewport(fb, -1, -2, 10, 11) == 0, "display viewport accepted");
    CHECK(state.viewport_count == 1 && state.viewport[0] == -1 &&
	    state.viewport[1] == -2 && state.viewport[2] == 10 &&
	    state.viewport[3] == 11,
	    "display host viewport callback received arguments");
    CHECK(imgstream_fb_reset(fb) == 0 && state.reset_count == 1,
	    "display reset forwarded to host");
    CHECK(imgstream_fb_flush(fb) == 0 && state.flush_count == 1,
	    "display flush forwarded to host");
    CHECK(imgstream_fb_poll(fb) == 12 && state.poll_count == 1,
	    "display poll forwarded to host");
    CHECK(imgstream_fb_poll_rate(fb) == 30 && state.poll_rate_count == 1,
	    "display poll rate forwarded to host");

    state2.poll_return = 24;
    state2.poll_rate_return = 60;
    imgstream_fb_t *fb2 = imgstream_fb_open_display("/dev/swrast", 7, 6, &host, &state2);
    CHECK(fb2 != NULL && state2.open_count == 1 && state.open_count == 1,
	    "two display framebuffers open with isolated host state");
    CHECK(imgstream_fb_poll(fb2) == 24 && state2.poll_count == 1 && state.poll_count == 1,
	    "second display framebuffer dispatches only to its own host data");
    imgstream_fb_close(fb);
    CHECK(state.close_count == 1 && state2.close_count == 0,
	    "closing one display framebuffer leaves the other host active");
    CHECK(imgstream_fb_poll_rate(fb2) == 60 && state2.poll_rate_count == 1,
	    "remaining display framebuffer stays usable after peer teardown");
    imgstream_fb_close(fb2);
    CHECK(state2.close_count == 1, "second display framebuffer closes its own host");

    memset(&state, 0, sizeof(state));
    state.open_return = -1;
    CHECK(imgstream_fb_open_display("/dev/qtgl", 2, 2, &host, &state) == NULL,
	    "display host open failure rejects fb open");
    CHECK(state.open_count == 1 && state.close_count == 0, "failed display open does not report close");

    return 0;
}


static int
test_display_host_detach(void)
{
    struct display_host_test_state state;
    struct imgstream_fb_display_host host;
    memset(&state, 0, sizeof(state));
    memset(&host, 0, sizeof(host));
    host.open = display_host_open;
    host.close = display_host_close;
    host.flush = display_host_flush;
    host.view = display_host_view;

    imgstream_fb_t *fb = imgstream_fb_open_display("/dev/qtgl", 3, 2,
	&host, &state);
    CHECK(fb != NULL && state.open_count == 1,
	"display host opens before detachment");
    CHECK(imgstream_fb_detach_display_host(fb) == 0,
	"display host detaches without closing the pixel stream");
    CHECK(imgstream_fb_flush(fb) == 0 && state.flush_count == 0,
	"detached display flush does not access the former host");
    CHECK(imgstream_fb_view(fb, 1, 1, 2, 2) == 0 && state.view_count == 0,
	"detached display view does not access the former host");
    imgstream_fb_close(fb);
    CHECK(state.close_count == 1,
	"display host closes exactly once during detachment");
    CHECK(imgstream_fb_detach_display_host(NULL) == -1,
	"null display host detachment is rejected");
    return 0;
}


static int
test_memory_compat_stream(void)
{
    imgstream_fb_t *fb = imgstream_fb_open("/dev/mem", 4, 3);
    CHECK(fb != NULL, "opened memory compatibility stream");
    CHECK(imgstream_fb_width(fb) == 4 && imgstream_fb_height(fb) == 3, "compat stream dimensions recorded");
    CHECK(imgstream_fb_stream(fb) != NULL, "compat stream exposes mutable stream");
    CHECK(imgstream_fb_cstream(fb) != NULL, "compat stream exposes const stream");
    CHECK(bu_strcmp(imgstream_fb_name(fb), "/dev/mem") == 0, "compat stream name recorded");

    int xcenter = 0;
    int ycenter = 0;
    int xzoom = 0;
    int yzoom = 0;
    CHECK(imgstream_fb_getview(fb, &xcenter, &ycenter, &xzoom, &yzoom) == 0, "default view query accepted");
    CHECK(xcenter == 2 && ycenter == 1 && xzoom == 1 && yzoom == 1, "default view initialized from size");
    CHECK(imgstream_fb_view(fb, 3, 2, 4, 5) == 0, "view update accepted");
    CHECK(imgstream_fb_getview(fb, &xcenter, &ycenter, &xzoom, &yzoom) == 0, "updated view query accepted");
    CHECK(xcenter == 3 && ycenter == 2 && xzoom == 4 && yzoom == 5, "updated view preserved");
    CHECK(imgstream_fb_window(fb, 1, 2) == 0, "legacy window update accepted");
    CHECK(imgstream_fb_getview(fb, &xcenter, &ycenter, &xzoom, &yzoom) == 0, "window view query accepted");
    CHECK(xcenter == 1 && ycenter == 2 && xzoom == 4 && yzoom == 5, "legacy window preserves zoom");
    CHECK(imgstream_fb_zoom(fb, 6, 7) == 0, "legacy zoom update accepted");
    CHECK(imgstream_fb_getview(fb, &xcenter, &ycenter, &xzoom, &yzoom) == 0, "zoom view query accepted");
    CHECK(xcenter == 1 && ycenter == 2 && xzoom == 6 && yzoom == 7, "legacy zoom preserves center");

    int mode = 0;
    int cx = 0;
    int cy = 0;
    CHECK(imgstream_fb_cursor(fb, 1, 2, 1) == 0, "cursor update accepted");
    CHECK(imgstream_fb_getcursor(fb, &mode, &cx, &cy) == 0, "cursor query accepted");
    CHECK(mode == 1 && cx == 2 && cy == 1, "cursor state preserved");
    CHECK(imgstream_fb_scursor(fb, 0, 3, 2) == 0, "legacy screen cursor spelling accepted");
    CHECK(imgstream_fb_getcursor(fb, &mode, &cx, &cy) == 0, "screen cursor cursor query accepted");
    CHECK(mode == 1 && cx == 2 && cy == 1, "legacy screen cursor does not change cursor state");
    CHECK(imgstream_fb_setcursor(fb, NULL, 16, 16, 8, 8) == 0, "legacy cursor-shape spelling accepted");
    CHECK(imgstream_fb_reset(fb) == 0, "legacy reset spelling accepted");
    CHECK(imgstream_fb_viewport(fb, -1, -2, 20, 30) == 0, "legacy viewport spelling accepted");
    CHECK(imgstream_fb_poll(fb) == 0, "non-display poll is a no-op success");
    CHECK(imgstream_fb_poll_rate(fb) == 0, "non-display poll rate is zero");

    unsigned char clear_rgb[3] = {7, 8, 9};
    CHECK(imgstream_fb_clear(fb, clear_rgb) == 0, "compat clear accepted");

    unsigned char pixel[3] = {0, 0, 0};
    CHECK(imgstream_fb_read(fb, 0, 0, pixel, 1) == 1, "single pixel read accepted");
    CHECK(pixel[0] == 7 && pixel[1] == 8 && pixel[2] == 9, "clear value read back");

    unsigned char linear[6 * 3];
    for (size_t i = 0; i < sizeof(linear); i++)
	linear[i] = (unsigned char)(30 + i);
    CHECK(imgstream_fb_write(fb, 2, 0, linear, 6) == 6, "linear write crossing rows accepted");

    unsigned char linear_out[6 * 3];
    memset(linear_out, 0, sizeof(linear_out));
    CHECK(imgstream_fb_read(fb, 2, 0, linear_out, 6) == 6, "linear read crossing rows accepted");
    CHECK(memcmp(linear, linear_out, sizeof(linear)) == 0, "linear read returns written bytes");

    unsigned char rect[2 * 2 * 3];
    for (size_t i = 0; i < sizeof(rect); i++)
	rect[i] = (unsigned char)(100 + i);
    CHECK(imgstream_fb_writerect(fb, 1, 1, 2, 2, rect) == 4, "rect write accepted");

    unsigned char rect_out[2 * 2 * 3];
    memset(rect_out, 0, sizeof(rect_out));
    CHECK(imgstream_fb_readrect(fb, 1, 1, 2, 2, rect_out) == 4, "rect read accepted");
    CHECK(memcmp(rect, rect_out, sizeof(rect)) == 0, "rect read returns written bytes");

    unsigned char bw_out[2 * 2];
    memset(bw_out, 0, sizeof(bw_out));
    CHECK(imgstream_fb_bwreadrect(fb, 1, 1, 2, 2, bw_out) == 4, "bw rect read accepted");
    CHECK(bw_out[0] == rect[1] && bw_out[1] == rect[4] &&
	    bw_out[2] == rect[7] && bw_out[3] == rect[10], "bw rect read returns green channel");

    unsigned char bw_in[2 * 2] = {11, 22, 33, 44};
    CHECK(imgstream_fb_bwwriterect(fb, 1, 1, 2, 2, bw_in) == 4, "bw rect write accepted");
    memset(rect_out, 0, sizeof(rect_out));
    CHECK(imgstream_fb_readrect(fb, 1, 1, 2, 2, rect_out) == 4, "bw-written rect RGB read accepted");
    for (size_t i = 0; i < sizeof(bw_in); i++) {
	CHECK(rect_out[i * 3 + 0] == bw_in[i] &&
		rect_out[i * 3 + 1] == bw_in[i] &&
		rect_out[i * 3 + 2] == bw_in[i], "bw rect write expands to RGB");
    }

    CHECK(imgstream_fb_read(fb, -1, 0, pixel, 1) == -1, "negative read coordinate rejected");
    CHECK(imgstream_fb_writerect(fb, 3, 2, 2, 1, rect) == -1, "out-of-bounds rect rejected");
    CHECK(imgstream_fb_bwreadrect(fb, 3, 2, 2, 1, bw_out) == -1, "out-of-bounds bw rect read rejected");
    CHECK(imgstream_fb_bwwriterect(fb, 3, 2, 2, 1, bw_in) == -1, "out-of-bounds bw rect write rejected");

    struct imgstream_info info;
    CHECK(imgstream_get_info(imgstream_fb_stream(fb), &info) == 0, "underlying stream info query accepted");
    CHECK(info.format == IMGSTREAM_PIXEL_RGB8 && info.channels == 3, "compat stream is RGB");
    CHECK(info.generation > 0, "compat writes advanced stream generation");

    imgstream_fb_close(fb);
    return 0;
}


static int
read_exact_file(const char *path, unsigned char *bytes, size_t byte_count)
{
    FILE *fp = fopen(path, "rb");
    if (!fp)
	return 0;

    size_t got = fread(bytes, 1, byte_count, fp);
    int extra = fgetc(fp);
    int ok = (got == byte_count && extra == EOF && !ferror(fp));
    fclose(fp);
    return ok;
}


static int
test_colormap_compat_state(void)
{
    struct imgstream_fb_colormap cmap;
    struct imgstream_fb_colormap readback;

    imgstream_fb_colormap_linear(NULL);
    CHECK(imgstream_fb_colormap_is_linear(NULL) == 0, "null colormap is not linear");
    imgstream_fb_colormap_linear(&cmap);
    CHECK(imgstream_fb_colormap_is_linear(&cmap) == 1, "generated linear colormap is linear");
    cmap.red[3] = 65535;
    CHECK(imgstream_fb_colormap_is_linear(&cmap) == 0, "edited colormap is not linear");

    imgstream_fb_t *fb = imgstream_fb_open("/dev/mem", 2, 2);
    CHECK(fb != NULL, "opened memory stream for colormap test");
    CHECK(imgstream_fb_rmap(fb, &readback) == 0, "default memory colormap read accepted");
    CHECK(imgstream_fb_colormap_is_linear(&readback) == 1, "default memory colormap is linear");

    CHECK(imgstream_fb_wmap(fb, &cmap) == 0, "memory colormap write accepted");
    memset(&readback, 0, sizeof(readback));
    CHECK(imgstream_fb_rmap(fb, &readback) == 0, "memory colormap reread accepted");
    CHECK(readback.red[3] == 65535 && readback.green[3] == (3 << 8) &&
	    readback.blue[3] == (3 << 8), "memory colormap state preserved");

    CHECK(imgstream_fb_wmap(fb, NULL) == 0, "memory null colormap resets to linear");
    memset(&readback, 0, sizeof(readback));
    CHECK(imgstream_fb_rmap(fb, &readback) == 0, "memory reset colormap read accepted");
    CHECK(imgstream_fb_colormap_is_linear(&readback) == 1, "memory reset colormap is linear");
    CHECK(imgstream_fb_rmap(NULL, &readback) == -1, "null fb colormap read rejected");
    CHECK(imgstream_fb_wmap(NULL, &readback) == -1, "null fb colormap write rejected");
    imgstream_fb_close(fb);

    fb = imgstream_fb_open("/dev/null", 2, 2);
    CHECK(fb != NULL, "opened diagnostic stream for colormap test");
    CHECK(imgstream_fb_wmap(fb, &cmap) == 0, "diagnostic colormap write accepted");
    memset(&readback, 0, sizeof(readback));
    CHECK(imgstream_fb_rmap(fb, &readback) == 0, "diagnostic colormap read accepted");
    CHECK(readback.red[3] == 65535, "diagnostic colormap state preserved");
    imgstream_fb_close(fb);

    char path_one[2048] = {0};
    char path_two[2048] = {0};
    FILE *tmp = bu_temp_file(path_one, sizeof(path_one));
    CHECK(tmp != NULL && path_one[0] != '\0', "first temporary colormap fanout file path created");
    fclose(tmp);
    (void)bu_file_delete(path_one);
    tmp = bu_temp_file(path_two, sizeof(path_two));
    CHECK(tmp != NULL && path_two[0] != '\0', "second temporary colormap fanout file path created");
    fclose(tmp);
    (void)bu_file_delete(path_two);

    char spec[5000] = {0};
    int len = snprintf(spec, sizeof(spec), "/dev/stack %s; /dev/disk:%s", path_one, path_two);
    CHECK(len > 0 && (size_t)len < sizeof(spec), "colormap fanout stack spec formatted");
    fb = imgstream_fb_open(spec, 2, 2);
    CHECK(fb != NULL, "opened fanout stream for colormap test");
    CHECK(imgstream_fb_wmap(fb, &cmap) == 0, "fanout colormap write accepted");
    memset(&readback, 0, sizeof(readback));
    CHECK(imgstream_fb_rmap(fb, &readback) == 0, "fanout colormap read accepted");
    CHECK(readback.red[3] == 65535, "fanout colormap read uses first child state");
    CHECK(imgstream_fb_wmap(fb, NULL) == 0, "fanout null colormap reset accepted");
    memset(&readback, 0, sizeof(readback));
    CHECK(imgstream_fb_rmap(fb, &readback) == 0, "fanout reset colormap read accepted");
    CHECK(imgstream_fb_colormap_is_linear(&readback) == 1, "fanout reset colormap is linear");
    imgstream_fb_close(fb);
    (void)bu_file_delete(path_one);
    (void)bu_file_delete(path_two);

    return 0;
}


static int
test_file_compat_stream(void)
{
    char path[2048] = {0};
    FILE *tmp = bu_temp_file(path, sizeof(path));
    CHECK(tmp != NULL && path[0] != '\0', "temporary file path created");
    fclose(tmp);
    (void)bu_file_delete(path);

    unsigned char frame[3 * 2 * 3];
    for (size_t i = 0; i < sizeof(frame); i++)
	frame[i] = (unsigned char)(10 + i);

    imgstream_fb_t *fb = imgstream_fb_open(path, 3, 2);
    CHECK(fb != NULL, "opened file compatibility stream");
    CHECK(imgstream_fb_width(fb) == 3 && imgstream_fb_height(fb) == 2, "file stream dimensions recorded");
    CHECK(bu_strcmp(imgstream_fb_name(fb), path) == 0, "file stream name recorded");
    CHECK(imgstream_fb_write(fb, 0, 0, frame, 6) == 6, "file stream linear write accepted");
    CHECK(imgstream_fb_flush(fb) == 0, "file stream flush accepted");

    unsigned char file_bytes[sizeof(frame)];
    memset(file_bytes, 0, sizeof(file_bytes));
    CHECK(read_exact_file(path, file_bytes, sizeof(file_bytes)), "file stream flushed exact raw RGB payload");
    CHECK(memcmp(file_bytes, frame, sizeof(frame)) == 0, "flushed file payload matches stream bytes");
    imgstream_fb_close(fb);

    fb = imgstream_fb_open(path, 3, 2);
    CHECK(fb != NULL, "reopened file compatibility stream");
    unsigned char reopened[sizeof(frame)];
    memset(reopened, 0, sizeof(reopened));
    CHECK(imgstream_fb_read(fb, 0, 0, reopened, 6) == 6, "reopened file stream read accepted");
    CHECK(memcmp(reopened, frame, sizeof(frame)) == 0, "reopened file stream loaded raw RGB payload");
    imgstream_fb_close(fb);

    char spec[2300] = {0};
    int len = snprintf(spec, sizeof(spec), "/dev/disk:%s", path);
    CHECK(len > 0 && (size_t)len < sizeof(spec), "explicit disk path spec formatted");
    unsigned char short_frame[2 * 1 * 3] = {201, 202, 203, 204, 205, 206};
    fb = imgstream_fb_open(spec, 2, 1);
    CHECK(fb != NULL, "opened explicit disk path compatibility stream");
    CHECK(imgstream_fb_writerect(fb, 0, 0, 2, 1, short_frame) == 2, "explicit disk path rect write accepted");
    imgstream_fb_close(fb);

    memset(file_bytes, 0, sizeof(file_bytes));
    CHECK(read_exact_file(path, file_bytes, sizeof(short_frame)), "explicit disk path close flushed exact raw RGB payload");
    CHECK(memcmp(file_bytes, short_frame, sizeof(short_frame)) == 0, "explicit disk path payload matches stream bytes");

    (void)bu_file_delete(path);
    return 0;
}


static int
test_fanout_compat_stream(void)
{
    char path_one[2048] = {0};
    char path_two[2048] = {0};
    FILE *tmp = bu_temp_file(path_one, sizeof(path_one));
    CHECK(tmp != NULL && path_one[0] != '\0', "first temporary fanout file path created");
    fclose(tmp);
    (void)bu_file_delete(path_one);

    tmp = bu_temp_file(path_two, sizeof(path_two));
    CHECK(tmp != NULL && path_two[0] != '\0', "second temporary fanout file path created");
    fclose(tmp);
    (void)bu_file_delete(path_two);

    char spec[5000] = {0};
    int len = snprintf(spec, sizeof(spec), "/dev/stack %s; /dev/disk:%s", path_one, path_two);
    CHECK(len > 0 && (size_t)len < sizeof(spec), "fanout stack spec formatted");
    CHECK(imgstream_fb_spec_supported("/dev/stack") == 0, "bare stack spec needs child sinks");
    CHECK(imgstream_fb_spec_supported(spec) == 1, "stack spec with supported children is supported");
    CHECK(imgstream_fb_spec_supported("/dev/stack /dev/qtgl") == 0, "stack spec with display child is not stream supported");

    unsigned char frame[3 * 2 * 3];
    for (size_t i = 0; i < sizeof(frame); i++)
	frame[i] = (unsigned char)(50 + i);

    imgstream_fb_t *fb = imgstream_fb_open(spec, 3, 2);
    CHECK(fb != NULL, "opened fanout compatibility stream");
    CHECK(imgstream_fb_width(fb) == 3 && imgstream_fb_height(fb) == 2, "fanout stream dimensions recorded");
    CHECK(imgstream_fb_stream(fb) != NULL, "fanout exposes first child stream");
    CHECK(bu_strcmp(imgstream_fb_name(fb), spec) == 0, "fanout stream name recorded");
    CHECK(imgstream_fb_writerect(fb, 0, 0, 3, 2, frame) == 6, "fanout rect write accepted");

    unsigned char readback[sizeof(frame)];
    memset(readback, 0, sizeof(readback));
    CHECK(imgstream_fb_read(fb, 0, 0, readback, 6) == 6, "fanout read uses first child");
    CHECK(memcmp(readback, frame, sizeof(frame)) == 0, "fanout readback matches written bytes");

    int xcenter = 0;
    int ycenter = 0;
    int xzoom = 0;
    int yzoom = 0;
    CHECK(imgstream_fb_view(fb, 2, 1, 3, 4) == 0, "fanout view update accepted");
    CHECK(imgstream_fb_getview(fb, &xcenter, &ycenter, &xzoom, &yzoom) == 0, "fanout view query accepted");
    CHECK(xcenter == 2 && ycenter == 1 && xzoom == 3 && yzoom == 4, "fanout view state preserved");
    CHECK(imgstream_fb_window(fb, 1, 1) == 0, "fanout window update accepted");
    CHECK(imgstream_fb_zoom(fb, 5, 6) == 0, "fanout zoom update accepted");
    CHECK(imgstream_fb_getview(fb, &xcenter, &ycenter, &xzoom, &yzoom) == 0, "fanout window/zoom view query accepted");
    CHECK(xcenter == 1 && ycenter == 1 && xzoom == 5 && yzoom == 6, "fanout window/zoom state preserved");

    int mode = 0;
    int cx = 0;
    int cy = 0;
    CHECK(imgstream_fb_cursor(fb, 1, 1, 1) == 0, "fanout cursor update accepted");
    CHECK(imgstream_fb_getcursor(fb, &mode, &cx, &cy) == 0, "fanout cursor query accepted");
    CHECK(mode == 1 && cx == 1 && cy == 1, "fanout cursor state preserved");
    CHECK(imgstream_fb_scursor(fb, 0, 2, 2) == 0, "fanout screen cursor spelling accepted");
    CHECK(imgstream_fb_setcursor(fb, NULL, 8, 8, 4, 4) == 0, "fanout cursor-shape spelling accepted");
    CHECK(imgstream_fb_reset(fb) == 0, "fanout reset accepted");
    CHECK(imgstream_fb_viewport(fb, 0, 0, 2, 1) == 0, "fanout viewport accepted");
    CHECK(imgstream_fb_poll(fb) == 0, "fanout poll accepted");
    CHECK(imgstream_fb_poll_rate(fb) == 0, "fanout poll rate is zero");

    unsigned char fanout_bw[2] = {12, 34};
    CHECK(imgstream_fb_bwwriterect(fb, 1, 1, 2, 1, fanout_bw) == 2, "fanout bw rect write accepted");
    unsigned char fanout_bw_out[2] = {0, 0};
    CHECK(imgstream_fb_bwreadrect(fb, 1, 1, 2, 1, fanout_bw_out) == 2, "fanout bw rect read accepted");
    CHECK(fanout_bw_out[0] == fanout_bw[0] && fanout_bw_out[1] == fanout_bw[1], "fanout bw rect read returns written values");
    for (size_t c = 0; c < 3; c++) {
	frame[(4 * 3) + c] = fanout_bw[0];
	frame[(5 * 3) + c] = fanout_bw[1];
    }
    CHECK(imgstream_fb_flush(fb) == 0, "fanout flush accepted");

    unsigned char file_bytes[sizeof(frame)];
    memset(file_bytes, 0, sizeof(file_bytes));
    CHECK(read_exact_file(path_one, file_bytes, sizeof(file_bytes)), "first fanout sink flushed exact raw RGB payload");
    CHECK(memcmp(file_bytes, frame, sizeof(frame)) == 0, "first fanout sink payload matches stream bytes");
    memset(file_bytes, 0, sizeof(file_bytes));
    CHECK(read_exact_file(path_two, file_bytes, sizeof(file_bytes)), "second fanout sink flushed exact raw RGB payload");
    CHECK(memcmp(file_bytes, frame, sizeof(frame)) == 0, "second fanout sink payload matches stream bytes");

    unsigned char clear_rgb[3] = {9, 8, 7};
    unsigned char clear_frame[sizeof(frame)];
    for (size_t i = 0; i < sizeof(clear_frame); i += 3) {
	clear_frame[i] = clear_rgb[0];
	clear_frame[i + 1] = clear_rgb[1];
	clear_frame[i + 2] = clear_rgb[2];
    }
    CHECK(imgstream_fb_clear(fb, clear_rgb) == 0, "fanout clear accepted");
    CHECK(imgstream_fb_flush(fb) == 0, "fanout clear flush accepted");
    memset(file_bytes, 0, sizeof(file_bytes));
    CHECK(read_exact_file(path_one, file_bytes, sizeof(file_bytes)), "first fanout sink flushed clear payload");
    CHECK(memcmp(file_bytes, clear_frame, sizeof(clear_frame)) == 0, "first fanout clear payload matches");
    memset(file_bytes, 0, sizeof(file_bytes));
    CHECK(read_exact_file(path_two, file_bytes, sizeof(file_bytes)), "second fanout sink flushed clear payload");
    CHECK(memcmp(file_bytes, clear_frame, sizeof(clear_frame)) == 0, "second fanout clear payload matches");

    imgstream_fb_close(fb);
    (void)bu_file_delete(path_one);
    (void)bu_file_delete(path_two);
    return 0;
}


static int
test_diagnostic_compat_stream(void)
{
    imgstream_fb_t *fb = imgstream_fb_open("/dev/null", 4, 3);
    CHECK(fb != NULL, "opened null diagnostic compatibility stream");
    CHECK(imgstream_fb_width(fb) == 4 && imgstream_fb_height(fb) == 3, "null diagnostic dimensions recorded");
    CHECK(imgstream_fb_stream(fb) == NULL, "null diagnostic retains no image stream");
    CHECK(imgstream_fb_cstream(fb) == NULL, "null diagnostic has no const image stream");
    CHECK(bu_strcmp(imgstream_fb_name(fb), "/dev/null") == 0, "null diagnostic name recorded");

    unsigned char pixel[3] = {200, 201, 202};
    CHECK(imgstream_fb_write(fb, -10, -20, pixel, 1) == 1, "null diagnostic write discards out-of-range pixels");
    CHECK(imgstream_fb_read(fb, -10, -20, NULL, 2) == 2, "null diagnostic read accepts null output buffer");

    unsigned char readback[2 * 3];
    memset(readback, 0xff, sizeof(readback));
    CHECK(imgstream_fb_read(fb, -10, -20, readback, 2) == 2, "null diagnostic read accepted");
    for (size_t i = 0; i < sizeof(readback); i++)
	CHECK(readback[i] == 0, "null diagnostic read returns zero pixels");

    unsigned char rect[2 * 1 * 3] = {1, 2, 3, 4, 5, 6};
    CHECK(imgstream_fb_writerect(fb, -1, -1, 2, 1, rect) == 2, "null diagnostic rect write accepted");
    memset(readback, 0xff, sizeof(readback));
    CHECK(imgstream_fb_readrect(fb, -1, -1, 2, 1, readback) == 2, "null diagnostic rect read accepted");
    for (size_t i = 0; i < sizeof(readback); i++)
	CHECK(readback[i] == 0, "null diagnostic rect read returns zero pixels");
    unsigned char bw[2] = {12, 34};
    CHECK(imgstream_fb_bwwriterect(fb, -1, -1, 2, 1, bw) == 2, "null diagnostic bw rect write accepted");
    memset(bw, 0xff, sizeof(bw));
    CHECK(imgstream_fb_bwreadrect(fb, -1, -1, 2, 1, bw) == 2, "null diagnostic bw rect read accepted");
    CHECK(bw[0] == 0 && bw[1] == 0, "null diagnostic bw rect read returns zero pixels");
    CHECK(imgstream_fb_readrect(fb, 0, 0, 0, 1, readback) == 0, "null diagnostic zero-width rect read is empty");
    CHECK(imgstream_fb_bwreadrect(fb, 0, 0, 0, 1, bw) == 0, "null diagnostic zero-width bw rect read is empty");

    int xcenter = 0;
    int ycenter = 0;
    int xzoom = 0;
    int yzoom = 0;
    CHECK(imgstream_fb_view(fb, 5, 6, 2, 3) == 0, "null diagnostic view update accepted");
    CHECK(imgstream_fb_getview(fb, &xcenter, &ycenter, &xzoom, &yzoom) == 0, "null diagnostic view query accepted");
    CHECK(xcenter == 5 && ycenter == 6 && xzoom == 2 && yzoom == 3, "null diagnostic view state preserved");
    CHECK(imgstream_fb_window(fb, 7, 8) == 0, "null diagnostic window update accepted");
    CHECK(imgstream_fb_zoom(fb, 4, 5) == 0, "null diagnostic zoom update accepted");
    CHECK(imgstream_fb_getview(fb, &xcenter, &ycenter, &xzoom, &yzoom) == 0, "null diagnostic window/zoom view query accepted");
    CHECK(xcenter == 7 && ycenter == 8 && xzoom == 4 && yzoom == 5, "null diagnostic window/zoom state preserved");

    int mode = 0;
    int cx = 0;
    int cy = 0;
    CHECK(imgstream_fb_cursor(fb, 1, 8, 9) == 0, "null diagnostic cursor update accepted");
    CHECK(imgstream_fb_getcursor(fb, &mode, &cx, &cy) == 0, "null diagnostic cursor query accepted");
    CHECK(mode == 1 && cx == 8 && cy == 9, "null diagnostic cursor state preserved");
    CHECK(imgstream_fb_scursor(fb, 0, 0, 0) == 0, "null diagnostic screen cursor spelling accepted");
    CHECK(imgstream_fb_setcursor(fb, NULL, 1, 1, 0, 0) == 0, "null diagnostic cursor-shape spelling accepted");
    CHECK(imgstream_fb_reset(fb) == 0, "null diagnostic reset accepted");
    CHECK(imgstream_fb_viewport(fb, -20, -30, 20, 30) == 0, "null diagnostic viewport accepted");
    CHECK(imgstream_fb_poll_rate(fb) == 0, "null diagnostic poll rate is zero");

    CHECK(imgstream_fb_flush(fb) == 0, "null diagnostic flush accepted");
    imgstream_fb_close(fb);

    fb = imgstream_fb_open("/dev/debug4", 2, 1);
    CHECK(fb != NULL, "opened debug diagnostic compatibility stream");
    CHECK(imgstream_fb_write(fb, 0, 0, pixel, 1) == 1, "debug diagnostic write accepted");
    CHECK(imgstream_fb_read(fb, 0, 0, readback, 1) == 1, "debug diagnostic read accepted");
    CHECK(readback[0] == 0 && readback[1] == 0 && readback[2] == 0, "debug diagnostic read returns zero pixel");
    imgstream_fb_close(fb);

    fb = imgstream_fb_open("/dev/txt", 2, 1);
    CHECK(fb != NULL, "opened txt diagnostic compatibility stream");
    CHECK(imgstream_fb_clear(fb, pixel) == 0, "txt diagnostic clear accepted");
    imgstream_fb_close(fb);

    return 0;
}


static int
test_image_io(void)
{
    static const unsigned char pixels[36] = {
	1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
	21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
	41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52
    };
    char input_path[MAXPATHLEN] = {0};
    FILE *input = bu_temp_file(input_path, sizeof(input_path));
    CHECK(input != NULL, "allocated PIX import file");
    CHECK(fwrite(pixels, 1, sizeof(pixels), input) == sizeof(pixels),
	"wrote PIX import fixture");
    CHECK(fflush(input) == 0 && fseek(input, 0, SEEK_SET) == 0,
	"rewound PIX import fixture");

    imgstream_fb_t *fb = imgstream_fb_open(NULL, 4, 3);
    CHECK(fb != NULL, "opened image-I/O framebuffer");
    struct imgstream_fb_import_options options =
	IMGSTREAM_FB_IMPORT_OPTIONS_INIT;
    CHECK(imgstream_fb_import_pix_fd(fb, fileno(input), input_path, 4, 3,
	0, &options) == 0, "imported PIX fixture");
    unsigned char readback[36] = {0};
    CHECK(imgstream_fb_readrect(fb, 0, 0, 4, 3, readback) == 12 &&
	memcmp(readback, pixels, sizeof(pixels)) == 0,
	"PIX import preserved bottom-up rows");

    char output_path[MAXPATHLEN] = {0};
    FILE *output = bu_temp_file(output_path, sizeof(output_path));
    CHECK(output != NULL, "allocated PIX export file");
    CHECK(imgstream_fb_export_pix_fp(fb, output, 4, 3, 0, 0) == 0 &&
	fflush(output) == 0 && fseek(output, 0, SEEK_SET) == 0,
	"exported PIX fixture");
    memset(readback, 0, sizeof(readback));
    CHECK(fread(readback, 1, sizeof(readback), output) == sizeof(readback) &&
	memcmp(readback, pixels, sizeof(pixels)) == 0,
	"PIX export round trip was exact");

    CHECK(fseek(input, 0, SEEK_SET) == 0, "rewound inverse PIX fixture");
    options.inverse = 1;
    CHECK(imgstream_fb_import_pix_fd(fb, fileno(input), input_path, 4, 3,
	0, &options) == 0, "imported inverse PIX fixture");
    memset(readback, 0, sizeof(readback));
    CHECK(imgstream_fb_readrect(fb, 0, 0, 4, 3, readback) == 12 &&
	memcmp(readback, pixels + 24, 12) == 0 &&
	memcmp(readback + 24, pixels, 12) == 0,
	"inverse PIX import flipped row placement");

    icv_image_t *image = icv_image_create(2, 1, ICV_COLOR_SPACE_RGB);
    CHECK(image != NULL, "created ICV import fixture");
    image->data[0] = 1.0;
    image->data[1] = 0.0;
    image->data[2] = 0.0;
    image->data[3] = 0.0;
    image->data[4] = 1.0;
    image->data[5] = 0.0;
    options = (struct imgstream_fb_import_options)
	IMGSTREAM_FB_IMPORT_OPTIONS_INIT;
    options.screen_xoff = 1;
    options.screen_yoff = 1;
    CHECK(imgstream_fb_import_icv(fb, image, &options) == 0,
	"imported ICV fixture");
    unsigned char image_row[6] = {0};
    CHECK(imgstream_fb_read(fb, 1, 1, image_row, 2) == 2 &&
	image_row[0] == 255 && image_row[1] == 0 && image_row[2] == 0 &&
	image_row[3] == 0 && image_row[4] == 255 && image_row[5] == 0,
	"ICV import preserved RGB values and placement");

    icv_destroy(image);
    imgstream_fb_close(fb);
    fclose(input);
    fclose(output);
    bu_file_delete(input_path);
    bu_file_delete(output_path);
    return 0;
}


static int
test_defaults(void)
{
    imgstream_fb_t *fb = imgstream_fb_open(NULL, 0, 0);
    CHECK(fb != NULL, "default memory compatibility stream opened");
    CHECK(imgstream_fb_width(fb) == 512 && imgstream_fb_height(fb) == 512, "default dimensions match legacy memory framebuffer");
    CHECK(imgstream_fb_reset(NULL) == 0, "null reset is a no-op success");
    CHECK(imgstream_fb_viewport(NULL, 0, 0, 0, 0) == 0, "null viewport is a no-op success");
    CHECK(imgstream_fb_poll_rate(NULL) == 0, "null poll rate is zero");
    imgstream_fb_close(fb);
    return 0;
}


int
main(int ac, char **av)
{
    bu_setprogname(av[0]);

    if (ac != 1)
	return 1;

    if (test_spec_classification())
	return 1;
    if (test_display_host_compat_stream())
	return 1;
    if (test_display_host_detach())
	return 1;
    if (test_memory_compat_stream())
	return 1;
    if (test_colormap_compat_state())
	return 1;
    if (test_file_compat_stream())
	return 1;
    if (test_fanout_compat_stream())
	return 1;
    if (test_diagnostic_compat_stream())
	return 1;
    if (test_image_io())
	return 1;
    if (test_defaults())
	return 1;

    return 0;
}
