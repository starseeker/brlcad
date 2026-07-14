/*                   F B S E R V _ A U T H . C
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
/** @file libimgstream/tests/fbserv_auth.c */

#include "common.h"

#include "bu/app.h"

#include <stdio.h>
#include <string.h>

#include "bu/defines.h"
#include "imgstream/fbserv.h"


static int dummy_is_listening(struct fbserv_obj *UNUSED(fbsp)) { return 0; }
static int dummy_listen_on_port(struct fbserv_obj *UNUSED(fbsp), int UNUSED(port)) { return 1; }
static void dummy_server_handler(struct fbserv_obj *UNUSED(fbsp)) {}
static void dummy_client_handler(struct fbserv_obj *UNUSED(fbsp), int UNUSED(sub), void *UNUSED(data)) {}
static void dummy_close_client_handler(struct fbserv_obj *UNUSED(fbsp), int UNUSED(sub)) {}


static int
dummy_info(void *UNUSED(ctx), struct fbserv_fb_info *info)
{
    if (!info)
	return -1;
    info->max_width = 10;
    info->max_height = 20;
    info->width = 5;
    info->height = 6;
    return 0;
}


static int
dummy_clear(void *ctx, const unsigned char rgb[3])
{
    int *marker = (int *)ctx;
    if (!marker || !rgb)
	return -1;
    *marker += rgb[0] + rgb[1] + rgb[2];
    return 7;
}


static ssize_t
dummy_write(void *ctx,
	int x,
	int y,
	const unsigned char *UNUSED(rgb),
	size_t count)
{
    int *marker = (int *)ctx;
    if (!marker)
	return -1;
    *marker += x + y + (int)count;
    return (ssize_t)count;
}


static int
dummy_flush(void *ctx)
{
    int *marker = (int *)ctx;
    if (!marker)
	return -1;
    *marker += 100;
    return 9;
}


static int
is_hex_token(const char *token)
{
    if (!token || strlen(token) != FBSERV_AUTH_TOKEN_LEN)
	return 0;

    for (size_t i = 0; i < FBSERV_AUTH_TOKEN_LEN; i++) {
	char c = token[i];
	if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f')))
	    return 0;
    }

    return 1;
}


int
main(void)
{
    bu_setprogname("fbserv_auth");
    char token[FBSERV_AUTH_TOKEN_LEN + 1] = {0};
    char mismatch[FBSERV_AUTH_TOKEN_LEN + 1] = {0};
    struct fbserv_obj fbs;
    int legacy_marker = 0;
    int backend_marker = 0;
    struct fbserv_fb_info info;
    unsigned char rgb[3] = {1, 2, 3};
    unsigned char pixels[6] = {0};
    static const struct fbserv_fb_ops ops = {
	dummy_info, dummy_clear, NULL, dummy_write, NULL, NULL, NULL, NULL,
	NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	dummy_flush, NULL, NULL
    };
    static const struct fbserv_transport_ops transport_ops = {
	dummy_is_listening,
	dummy_listen_on_port,
	dummy_server_handler,
	dummy_server_handler,
	dummy_client_handler,
	dummy_close_client_handler,
	dummy_client_handler,
	dummy_close_client_handler
    };

    fbserv_generate_token(token);
    if (!is_hex_token(token)) {
	fprintf(stderr, "generated token is not a %d-char lowercase hex string\n",
	    FBSERV_AUTH_TOKEN_LEN);
	return 1;
    }

    if (!fbserv_verify_token(token, token)) {
	fprintf(stderr, "token did not verify against itself\n");
	return 1;
    }

    memcpy(mismatch, token, sizeof(mismatch));
    mismatch[0] = (mismatch[0] == '0') ? '1' : '0';
    if (fbserv_verify_token(mismatch, token)) {
	fprintf(stderr, "mismatched token verified\n");
	return 1;
    }

    if (fbserv_verify_token(NULL, token) || fbserv_verify_token(token, NULL) ||
	fbserv_verify_token("short", token)) {
	fprintf(stderr, "invalid token inputs verified\n");
	return 1;
    }

    if (!FBSERV_AUTH_TOKEN_ENVVAR[0]) {
	fprintf(stderr, "auth token environment variable name is empty\n");
	return 1;
    }

    if (fbserv_obj_init(&fbs) != BRLCAD_OK) {
	fprintf(stderr, "fbserv object init failed\n");
	return 1;
    }
    if (fbserv_listener_fd(&fbs) != -1 || fbserv_listener_port(&fbs) != -1 ||
	fbserv_listener_owner(fbserv_listener_handler_data(&fbs)) != &fbs) {
	fprintf(stderr, "fbserv object init did not set listener sentinels/owner\n");
	return 1;
    }

    if (fbserv_set_backend(&fbs, &ops, &backend_marker) != BRLCAD_OK ||
	!fbserv_framebuffer_backend_installed(&fbs)) {
	fprintf(stderr, "fbserv backend registration failed\n");
	return 1;
    }
    if (!fbserv_backend_op_installed(&fbs, FBSERV_BACKEND_OP_INFO) ||
	!fbserv_backend_op_installed(&fbs, FBSERV_BACKEND_OP_CLEAR) ||
	!fbserv_backend_op_installed(&fbs, FBSERV_BACKEND_OP_WRITE) ||
	fbserv_backend_op_installed(&fbs, FBSERV_BACKEND_OP_READ)) {
	fprintf(stderr, "fbserv backend operation capability check failed\n");
	return 1;
    }
    if (fbserv_backend_info(&fbs, &info) != 0 ||
	info.max_width != 10 || info.max_height != 20 ||
	info.width != 5 || info.height != 6) {
	fprintf(stderr, "fbserv backend info dispatch failed\n");
	return 1;
    }
    if (fbserv_backend_clear(&fbs, rgb) != 7 || backend_marker != 6) {
	fprintf(stderr, "fbserv backend clear dispatch failed\n");
	return 1;
    }
    if (fbserv_backend_write(&fbs, 3, 4, pixels, 2) != 2 ||
	backend_marker != 15) {
	fprintf(stderr, "fbserv backend write dispatch failed\n");
	return 1;
    }
    if (fbserv_backend_read(&fbs, 0, 0, pixels, 1) !=
	FBSERV_FRAMEBUFFER_UNHANDLED) {
	fprintf(stderr, "fbserv backend unsupported dispatch failed\n");
	return 1;
    }
    if (fbserv_backend_flush(&fbs) != 9 || backend_marker != 115) {
	fprintf(stderr, "fbserv backend flush dispatch failed\n");
	return 1;
    }
    fbserv_clear_backend(&fbs);
    if (fbserv_framebuffer_backend_installed(&fbs) ||
	fbserv_backend_op_installed(&fbs, FBSERV_BACKEND_OP_INFO)) {
	fprintf(stderr, "fbserv backend clear failed\n");
	return 1;
    }

    fbserv_set_transport(&fbs, &transport_ops);
    if (!fbserv_can_open_ipc(&fbs) || !fbserv_can_open_network(&fbs) ||
	!fbserv_can_close(&fbs)) {
	fprintf(stderr, "fbserv transport capability checks failed\n");
	return 1;
    }
    fbserv_clear_transport(&fbs);
    if (fbserv_can_open_ipc(&fbs) || fbserv_can_open_network(&fbs) ||
	fbserv_can_close(&fbs)) {
	fprintf(stderr, "fbserv transport clear failed\n");
	return 1;
    }

    fbserv_set_listener_fd(&fbs, 123);
    fbserv_set_listener_channel(&fbs, &legacy_marker);
    if (fbserv_listener_fd(&fbs) != 123 ||
	fbserv_listener_channel(&fbs) != &legacy_marker ||
	fbserv_listener_data_fd(fbserv_listener_handler_data(&fbs)) != 123) {
	fprintf(stderr, "fbserv listener accessors failed\n");
	return 1;
    }

    fbserv_set_client_channel(&fbs, 0, &legacy_marker);
    fbserv_set_client_handler(&fbs, 0, &backend_marker);
    fbs.fbs_clients[0].fbsc_fd = 321;
    if (!fbserv_client_active(&fbs, 0) ||
	fbserv_client_fd(&fbs, 0) != 321 ||
	fbserv_client_channel(&fbs, 0) != &legacy_marker ||
	fbserv_client_handler(&fbs, 0) != &backend_marker ||
	fbserv_client_handler_data(&fbs, 0) != &fbs.fbs_clients[0]) {
	fprintf(stderr, "fbserv client accessors failed\n");
	return 1;
    }

    fbserv_set_legacy_framebuffer(&fbs, (struct fb *)(void *)&legacy_marker);
    if (fbserv_legacy_framebuffer(&fbs) != (struct fb *)(void *)&legacy_marker) {
	fprintf(stderr, "fbserv legacy framebuffer accessor failed\n");
	return 1;
    }

    return 0;
}


/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
