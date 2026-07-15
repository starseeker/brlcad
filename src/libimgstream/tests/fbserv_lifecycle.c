/*              F B S E R V _ L I F E C Y C L E . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libimgstream/tests/fbserv_lifecycle.c */

#include "common.h"

#include "bu/app.h"

#include <stdlib.h>
#include <string.h>

#include "bu/log.h"
#include "imgstream/fb_compat.h"
#include "imgstream/fbserv.h"
#include "pkg.h"

struct test_backend {
    int writes;
    int x;
    int y;
    int width;
    int height;
    unsigned char pixels[12];
};

static int close_count = 0;
static int update_count = 0;

static void comm_error(const char *msg) { if (msg) bu_log("%s", msg); }
static int not_listening(struct fbserv_obj *UNUSED(fbsp)) { return 0; }
static int no_listen(struct fbserv_obj *UNUSED(fbsp), int UNUSED(port)) { return 0; }
static void no_server_handler(struct fbserv_obj *UNUSED(fbsp)) { }
static void no_client_handler(struct fbserv_obj *UNUSED(fbsp),
	int UNUSED(sub), void *UNUSED(data)) { }
static void close_client_handler(struct fbserv_obj *UNUSED(fbsp),
	int UNUSED(sub)) { close_count++; }
static void update_callback(void *UNUSED(data)) { update_count++; }

static int
backend_info(void *UNUSED(ctx), struct fbserv_fb_info *info)
{
    if (!info)
	return -1;
    info->max_width = info->width = 2;
    info->max_height = info->height = 2;
    return 0;
}

static int
backend_writerect(void *ctx, int x, int y, int width, int height,
	const unsigned char *pixels)
{
    struct test_backend *backend = (struct test_backend *)ctx;
    size_t pixel_bytes = (size_t)width * (size_t)height * 3;
    if (!backend || !pixels || pixel_bytes > sizeof(backend->pixels))
	return -1;
    backend->writes++;
    backend->x = x;
    backend->y = y;
    backend->width = width;
    backend->height = height;
    memcpy(backend->pixels, pixels, pixel_bytes);
    return width * height;
}

static size_t
append_packet(unsigned char *wire, size_t offset, unsigned short type,
	const unsigned char *payload, size_t payload_size)
{
    (void)pkg_pshort((char *)&wire[offset], PKG_MAGIC);
    (void)pkg_pshort((char *)&wire[offset + 2], type);
    (void)pkg_plong((char *)&wire[offset + 4], payload_size);
    if (payload_size)
	memcpy(&wire[offset + sizeof(struct pkg_header)], payload, payload_size);
    return offset + sizeof(struct pkg_header) + payload_size;
}

int
main(void)
{
    bu_setprogname("fbserv_lifecycle");
    static const struct fbserv_fb_ops backend_ops = {
	.info = backend_info,
	.writerect = backend_writerect
    };
    static const struct fbserv_transport_ops transport_ops = {
	not_listening,
	no_listen,
	no_server_handler,
	no_server_handler,
	no_client_handler,
	close_client_handler,
	no_client_handler,
	close_client_handler
    };
    struct test_backend backend = {0};
    struct fbserv_obj fbs;
    struct pkg_conn *server_end = PKC_NULL;
    struct pkg_conn *client_end = PKC_NULL;
    unsigned char rect[4 * NET_LONG_LEN + 12] = {0};
    unsigned char wire[256] = {0};
    unsigned char *remaining = NULL;
    size_t remaining_size = 0;
    size_t wire_size = 0;
    int slot;

    if (fbs_init(&fbs) != BRLCAD_OK ||
	fbs_set_backend(&fbs, &backend_ops, &backend) != BRLCAD_OK)
	return 1;
    fbs_set_transport(&fbs, &transport_ops);
    fbs.fbs_require_auth = 1;
    memset(fbs.fbs_auth_token, 'a', FBSERV_AUTH_TOKEN_LEN);
    fbs.fbs_auth_token[FBSERV_AUTH_TOKEN_LEN] = '\0';
    fbs.fbs_callback = update_callback;

    if (pkg_pair(&server_end, &client_end, fbs_pkg_switch(), comm_error) != 0)
	return 1;
    slot = fbs_new_client(&fbs, server_end, NULL);
    if (slot < 0)
	return 1;

    wire_size = append_packet(wire, wire_size, FBSERV_MSG_FBAUTH,
	(const unsigned char *)fbs.fbs_auth_token, FBSERV_AUTH_TOKEN_LEN);
    (void)pkg_plong((char *)&rect[0 * NET_LONG_LEN], 3);
    (void)pkg_plong((char *)&rect[1 * NET_LONG_LEN], 4);
    (void)pkg_plong((char *)&rect[2 * NET_LONG_LEN], 2);
    (void)pkg_plong((char *)&rect[3 * NET_LONG_LEN], 2);
    for (size_t i = 4 * NET_LONG_LEN; i < sizeof(rect); i++)
	rect[i] = (unsigned char)i;
    wire_size = append_packet(wire, wire_size,
	FBSERV_MSG_FBWRITERECT + FBSERV_MSG_NORETURN, rect, sizeof(rect));

    /* Deliver one byte at a time, preserving the unconsumed fragment exactly
     * as a toolkit host must when a packet spans event-loop notifications. */
    for (size_t i = 0; i < wire_size; i++) {
	struct fbserv_process_result result;
	size_t aggregate_size = remaining_size + 1;
	unsigned char *aggregate = (unsigned char *)malloc(aggregate_size);
	if (!aggregate)
	    return 1;
	if (remaining_size)
	    memcpy(aggregate, remaining, remaining_size);
	aggregate[remaining_size] = wire[i];
	free(remaining);
	remaining = NULL;
	remaining_size = 0;
	if (fbs_process_client_bytes(&fbs, slot, aggregate, aggregate_size,
		&remaining, &remaining_size, &result) != BRLCAD_OK) {
	    free(aggregate);
	    return 1;
	}
	free(aggregate);
    }
    free(remaining);

    if (remaining_size != 0 || backend.writes != 1 || backend.x != 3 ||
	backend.y != 4 || backend.width != 2 || backend.height != 2 ||
	backend.pixels[0] != 4 * NET_LONG_LEN ||
	backend.pixels[11] != 4 * NET_LONG_LEN + 11 || update_count < 2) {
	bu_log("fragmented authenticated framebuffer dispatch failed\n");
	return 1;
    }

    pkg_close(client_end);
    client_end = PKC_NULL;
    struct fbserv_process_result drain_result;
    if (fbs_drain_client(&fbs, slot, 4096, 1, 10000,
	    &drain_result) != BRLCAD_OK || !drain_result.disconnected ||
	fbs_client_active(&fbs, slot) || close_count != 1) {
	bu_log("EOF did not deterministically close the framebuffer client\n");
	return 1;
    }

    /* Server shutdown is the host-side cancellation path. */
    if (pkg_pair(&server_end, &client_end, fbs_pkg_switch(), comm_error) != 0)
	return 1;
    slot = fbs_new_client(&fbs, server_end, NULL);
    if (slot < 0 || fbs_close(&fbs) != BRLCAD_OK ||
	fbs_client_active(&fbs, slot) || close_count != 2) {
	bu_log("server cancellation did not close the active client\n");
	return 1;
    }
    pkg_close(client_end);

    fbs_clear_backend(&fbs);
    fbs_clear_transport(&fbs);

    struct fbserv_obj native_fbs;
    imgstream_fb_t *native_fb = imgstream_fb_open("/dev/mem", 2, 2);
    unsigned char native_pixels[12] = {
	1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
    };
    unsigned char native_readback[12] = {0};
    struct fbserv_fb_info native_info = {0};
    if (!native_fb || fbs_init(&native_fbs) != BRLCAD_OK ||
	imgstream_fbserv_set_framebuffer(&native_fbs, native_fb) != BRLCAD_OK ||
	!fbs_framebuffer_backend_installed(&native_fbs) ||
	fbserv_backend_info(&native_fbs, &native_info) != 0 ||
	native_info.width != 2 || native_info.height != 2 ||
	fbserv_backend_writerect(&native_fbs, 0, 0, 2, 2,
	    native_pixels) != 4 ||
	fbserv_backend_readrect(&native_fbs, 0, 0, 2, 2,
	    native_readback) != 4 ||
	memcmp(native_pixels, native_readback, sizeof(native_pixels)) != 0) {
	bu_log("imgstream-native framebuffer backend installation failed\n");
	return 1;
    }
    if (imgstream_fbserv_set_framebuffer(&native_fbs, NULL) != BRLCAD_OK ||
	fbs_framebuffer_backend_installed(&native_fbs)) {
	bu_log("imgstream-native framebuffer backend detach failed\n");
	return 1;
    }
    imgstream_fb_close(native_fb);
    return 0;
}
