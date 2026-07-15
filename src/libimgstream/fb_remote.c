/*                 F B _ R E M O T E . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libimgstream/fb_remote.c */

#include "common.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bu/log.h"
#include "bu/malloc.h"
#include "imgstream/fbserv.h"
#include "pkg.h"

#include "fb_remote.h"

#define FBSERV_TLS_IMPL
#define FBSERV_TLS_CLIENT
#include "../fbserv/tls_wrap.h"

#define REMOTE_CMAP_BYTES (256 * 3 * 2)

struct imgstream_fb_remote {
    struct pkg_conn *connection;
};

static void
remote_log(const char *message)
{
    if (message)
	bu_log("%s", message);
}

static void
remote_error(struct pkg_conn *UNUSED(connection), char *message)
{
    remote_log(message);
    free(message);
}

static struct pkg_switch remote_switch[] = {
    {FBSERV_MSG_ERROR, remote_error, "fbserv error", NULL},
    {0, NULL, NULL, NULL}
};

static char *
span_copy(const char *value, size_t length)
{
    char *copy = (char *)bu_malloc(length + 1, "remote framebuffer span");
    if (length)
	memcpy(copy, value, length);
    copy[length] = '\0';
    return copy;
}

static int
remote_command(struct imgstream_fb_remote *remote, int type,
	const void *payload, size_t payload_size)
{
    char response[FBSERV_NET_LONG_LEN] = {0};
    if (!remote || !remote->connection || payload_size > INT_MAX ||
	pkg_send(type, (const char *)payload, payload_size,
	    remote->connection) != (int)payload_size ||
	pkg_waitfor(FBSERV_MSG_RETURN, response, sizeof(response),
	    remote->connection) < FBSERV_NET_LONG_LEN)
	return -1;
    return (int)(int32_t)pkg_glong(response);
}

static int
remote_pack_rect(char payload[4 * FBSERV_NET_LONG_LEN], int x, int y,
	int width, int height, size_t bytes_per_pixel, size_t *byte_count)
{
    if (!byte_count || width < 0 || height < 0 ||
	(size_t)width > SIZE_MAX / (size_t)(height ? height : 1) ||
	(size_t)width * (size_t)height > SIZE_MAX / bytes_per_pixel)
	return 0;
    *byte_count = (size_t)width * (size_t)height * bytes_per_pixel;
    (void)pkg_plong(payload + 0 * FBSERV_NET_LONG_LEN, (uint32_t)x);
    (void)pkg_plong(payload + 1 * FBSERV_NET_LONG_LEN, (uint32_t)y);
    (void)pkg_plong(payload + 2 * FBSERV_NET_LONG_LEN, (uint32_t)width);
    (void)pkg_plong(payload + 3 * FBSERV_NET_LONG_LEN, (uint32_t)height);
    return 1;
}

int
imgstream_fb_remote_open(const imgstream_fb_spec_info_t *info,
	size_t requested_width, size_t requested_height,
	const struct imgstream_fb_remote_options *options,
	struct imgstream_fb_remote **remote_out, size_t *width, size_t *height)
{
    struct pkg_conn *connection = PKC_NULL;
    char *host = NULL;
    char *target = NULL;
    char *device = NULL;
    char service[32] = {0};
    char response[5 * FBSERV_NET_LONG_LEN] = {0};
    char *request = NULL;
    const char *token = options ? options->auth_token : NULL;
    struct imgstream_fb_remote *remote = NULL;
    size_t request_size = 0;
    int ret = 0;

    if (!info || info->kind != IMGSTREAM_FB_SPEC_REMOTE || !remote_out ||
	!width || !height || requested_width > UINT32_MAX ||
	requested_height > UINT32_MAX)
	return 0;
    *remote_out = NULL;
    *width = 0;
    *height = 0;

    if (info->remote == IMGSTREAM_FB_REMOTE_IPC) {
	if (!info->target || !info->target_len)
	    return 0;
	target = span_copy(info->target, info->target_len);
	connection = pkg_connect_addr(target, remote_switch, remote_log);
    } else if (info->remote == IMGSTREAM_FB_REMOTE_TCP) {
	const char *inherited = getenv(PKG_ADDR_ENVVAR);
	if (inherited && inherited[0])
	    connection = pkg_connect_addr(inherited, remote_switch, remote_log);
	if (connection == PKC_NULL || connection == PKC_ERROR) {
	    if (!info->host || !info->host_len || info->port < 0)
		goto cleanup;
	    host = span_copy(info->host, info->host_len);
	    snprintf(service, sizeof(service), "%d", info->port);
	    connection = pkg_open(host, service, NULL, NULL, NULL,
		remote_switch, remote_log);
	}
    } else {
	goto cleanup;
    }
    if (connection == PKC_NULL || connection == PKC_ERROR)
	goto cleanup;

#ifdef HAVE_OPENSSL_SSL_H
    if (options && options->use_tls) {
	SSL_CTX *tls_context = fbserv_tls_client_ctx();
	if (!tls_context ||
	    fbserv_tls_connect(tls_context, connection) != FBSERV_TLS_OK) {
	    if (tls_context)
		SSL_CTX_free(tls_context);
	    goto cleanup;
	}
	SSL_CTX_free(tls_context);
    }
#else
    if (options && options->use_tls)
	goto cleanup;
#endif

    if (token && token[0] &&
	pkg_send(FBSERV_MSG_FBAUTH, token, strlen(token), connection) !=
	    (int)strlen(token))
	goto cleanup;

    if (info->device && info->device_len)
	device = span_copy(info->device, info->device_len);
    request_size = 2 * FBSERV_NET_LONG_LEN + (device ? strlen(device) : 0);
    request = (char *)bu_calloc(request_size + 1, 1,
	"remote framebuffer open request");
    (void)pkg_plong(request + 0 * FBSERV_NET_LONG_LEN,
	(uint32_t)requested_width);
    (void)pkg_plong(request + 1 * FBSERV_NET_LONG_LEN,
	(uint32_t)requested_height);
    if (device)
	memcpy(request + 2 * FBSERV_NET_LONG_LEN, device, strlen(device));
    if (pkg_send(FBSERV_MSG_FBOPEN, request, request_size, connection) !=
	    (int)request_size ||
	pkg_waitfor(FBSERV_MSG_RETURN, response, sizeof(response), connection) <
	    5 * FBSERV_NET_LONG_LEN ||
	(int32_t)pkg_glong(response) != 0)
	goto cleanup;

    *width = (size_t)pkg_glong(response + 3 * FBSERV_NET_LONG_LEN);
    *height = (size_t)pkg_glong(response + 4 * FBSERV_NET_LONG_LEN);
    if (!*width || !*height)
	goto cleanup;

    BU_ALLOC(remote, struct imgstream_fb_remote);
    remote->connection = connection;
    *remote_out = remote;
    connection = PKC_NULL;
    ret = 1;

cleanup:
    if (connection != PKC_NULL && connection != PKC_ERROR)
	pkg_close(connection);
    if (request)
	bu_free(request, "remote framebuffer open request");
    if (device)
	bu_free(device, "remote framebuffer span");
    if (target)
	bu_free(target, "remote framebuffer span");
    if (host)
	bu_free(host, "remote framebuffer span");
    return ret;
}

void
imgstream_fb_remote_close(struct imgstream_fb_remote *remote)
{
    if (!remote)
	return;
    if (remote->connection) {
	char response[FBSERV_NET_LONG_LEN] = {0};
	(void)pkg_send(FBSERV_MSG_FBCLOSE, NULL, 0, remote->connection);
	(void)pkg_waitfor(FBSERV_MSG_RETURN, response, sizeof(response),
	    remote->connection);
	pkg_close(remote->connection);
    }
    BU_FREE(remote, struct imgstream_fb_remote);
}

int
imgstream_fb_remote_clear(struct imgstream_fb_remote *remote,
	const unsigned char *rgb)
{
    const unsigned char black[3] = {0, 0, 0};
    return remote_command(remote, FBSERV_MSG_FBCLEAR, rgb ? rgb : black, 3);
}

ssize_t
imgstream_fb_remote_read(struct imgstream_fb_remote *remote, int x, int y,
	unsigned char *rgb, size_t count)
{
    char payload[3 * FBSERV_NET_LONG_LEN] = {0};
    if (!remote || !rgb || count > UINT32_MAX || count > INT_MAX / 3)
	return -1;
    if (!count)
	return 0;
    (void)pkg_plong(payload + 0 * FBSERV_NET_LONG_LEN, (uint32_t)x);
    (void)pkg_plong(payload + 1 * FBSERV_NET_LONG_LEN, (uint32_t)y);
    (void)pkg_plong(payload + 2 * FBSERV_NET_LONG_LEN, (uint32_t)count);
    if (pkg_send(FBSERV_MSG_FBREAD, payload, sizeof(payload),
	    remote->connection) != (int)sizeof(payload))
	return -1;
    int bytes = pkg_waitfor(FBSERV_MSG_RETURN, (char *)rgb, count * 3,
	remote->connection);
    return bytes > 0 ? bytes / 3 : -1;
}

ssize_t
imgstream_fb_remote_write(struct imgstream_fb_remote *remote, int x, int y,
	const unsigned char *rgb, size_t count)
{
    char payload[3 * FBSERV_NET_LONG_LEN] = {0};
    if (!remote || (!rgb && count) || count > UINT32_MAX ||
	count > (size_t)INT_MAX / 3)
	return -1;
    if (!count)
	return 0;
    (void)pkg_plong(payload + 0 * FBSERV_NET_LONG_LEN, (uint32_t)x);
    (void)pkg_plong(payload + 1 * FBSERV_NET_LONG_LEN, (uint32_t)y);
    (void)pkg_plong(payload + 2 * FBSERV_NET_LONG_LEN, (uint32_t)count);
    int sent = pkg_2send(FBSERV_MSG_FBWRITE + FBSERV_MSG_NORETURN,
	payload, sizeof(payload), (const char *)rgb, count * 3,
	remote->connection);
    return sent == (int)(sizeof(payload) + count * 3) ? (ssize_t)count : -1;
}

int
imgstream_fb_remote_readrect(struct imgstream_fb_remote *remote, int x, int y,
	int width, int height, unsigned char *rgb)
{
    char payload[4 * FBSERV_NET_LONG_LEN] = {0};
    size_t byte_count = 0;
    if (!remote || !rgb || !remote_pack_rect(payload, x, y, width, height,
	3, &byte_count) || byte_count > INT_MAX)
	return -1;
    if (!byte_count)
	return 0;
    if (pkg_send(FBSERV_MSG_FBREADRECT, payload, sizeof(payload),
	    remote->connection) != (int)sizeof(payload))
	return -1;
    int bytes = pkg_waitfor(FBSERV_MSG_RETURN, (char *)rgb, byte_count,
	remote->connection);
    return bytes > 0 ? bytes / 3 : -1;
}

int
imgstream_fb_remote_writerect(struct imgstream_fb_remote *remote, int x,
	int y, int width, int height, const unsigned char *rgb)
{
    char payload[4 * FBSERV_NET_LONG_LEN] = {0};
    size_t byte_count = 0;
    if (!remote || !rgb || !remote_pack_rect(payload, x, y, width, height,
	3, &byte_count) || byte_count > INT_MAX)
	return -1;
    if (!byte_count)
	return 0;
    int sent = pkg_2send(FBSERV_MSG_FBWRITERECT + FBSERV_MSG_NORETURN,
	payload, sizeof(payload), (const char *)rgb, byte_count,
	remote->connection);
    return sent == (int)(sizeof(payload) + byte_count) ? width * height : -1;
}

int
imgstream_fb_remote_bwreadrect(struct imgstream_fb_remote *remote, int x,
	int y, int width, int height, unsigned char *bw)
{
    char payload[4 * FBSERV_NET_LONG_LEN] = {0};
    size_t byte_count = 0;
    if (!remote || !bw || !remote_pack_rect(payload, x, y, width, height,
	1, &byte_count) || byte_count > INT_MAX)
	return -1;
    if (!byte_count)
	return 0;
    if (pkg_send(FBSERV_MSG_FBBWREADRECT, payload, sizeof(payload),
	    remote->connection) != (int)sizeof(payload))
	return -1;
    int bytes = pkg_waitfor(FBSERV_MSG_RETURN, (char *)bw, byte_count,
	remote->connection);
    return bytes > 0 ? bytes : -1;
}

int
imgstream_fb_remote_bwwriterect(struct imgstream_fb_remote *remote, int x,
	int y, int width, int height, const unsigned char *bw)
{
    char payload[4 * FBSERV_NET_LONG_LEN] = {0};
    size_t byte_count = 0;
    if (!remote || !bw || !remote_pack_rect(payload, x, y, width, height,
	1, &byte_count) || byte_count > INT_MAX)
	return -1;
    if (!byte_count)
	return 0;
    int sent = pkg_2send(FBSERV_MSG_FBBWWRITERECT + FBSERV_MSG_NORETURN,
	payload, sizeof(payload), (const char *)bw, byte_count,
	remote->connection);
    return sent == (int)(sizeof(payload) + byte_count) ? (int)byte_count : -1;
}

int
imgstream_fb_remote_view(struct imgstream_fb_remote *remote, int xcenter,
	int ycenter, int xzoom, int yzoom)
{
    char payload[4 * FBSERV_NET_LONG_LEN] = {0};
    (void)pkg_plong(payload + 0 * FBSERV_NET_LONG_LEN, (uint32_t)xcenter);
    (void)pkg_plong(payload + 1 * FBSERV_NET_LONG_LEN, (uint32_t)ycenter);
    (void)pkg_plong(payload + 2 * FBSERV_NET_LONG_LEN, (uint32_t)xzoom);
    (void)pkg_plong(payload + 3 * FBSERV_NET_LONG_LEN, (uint32_t)yzoom);
    return remote_command(remote, FBSERV_MSG_FBVIEW, payload, sizeof(payload));
}

int
imgstream_fb_remote_getview(struct imgstream_fb_remote *remote, int *xcenter,
	int *ycenter, int *xzoom, int *yzoom)
{
    char response[5 * FBSERV_NET_LONG_LEN] = {0};
    if (!remote || !xcenter || !ycenter || !xzoom || !yzoom ||
	pkg_send(FBSERV_MSG_FBGETVIEW, NULL, 0, remote->connection) < 0 ||
	pkg_waitfor(FBSERV_MSG_RETURN, response, sizeof(response),
	    remote->connection) < 5 * FBSERV_NET_LONG_LEN ||
	(int32_t)pkg_glong(response) != 0)
	return -1;
    *xcenter = (int)(int32_t)pkg_glong(response + 1 * FBSERV_NET_LONG_LEN);
    *ycenter = (int)(int32_t)pkg_glong(response + 2 * FBSERV_NET_LONG_LEN);
    *xzoom = (int)(int32_t)pkg_glong(response + 3 * FBSERV_NET_LONG_LEN);
    *yzoom = (int)(int32_t)pkg_glong(response + 4 * FBSERV_NET_LONG_LEN);
    return 0;
}

int
imgstream_fb_remote_cursor(struct imgstream_fb_remote *remote, int mode,
	int x, int y)
{
    char payload[3 * FBSERV_NET_LONG_LEN] = {0};
    (void)pkg_plong(payload + 0 * FBSERV_NET_LONG_LEN, (uint32_t)mode);
    (void)pkg_plong(payload + 1 * FBSERV_NET_LONG_LEN, (uint32_t)x);
    (void)pkg_plong(payload + 2 * FBSERV_NET_LONG_LEN, (uint32_t)y);
    return remote_command(remote, FBSERV_MSG_FBCURSOR, payload, sizeof(payload));
}

int
imgstream_fb_remote_getcursor(struct imgstream_fb_remote *remote, int *mode,
	int *x, int *y)
{
    char response[4 * FBSERV_NET_LONG_LEN] = {0};
    if (!remote || !mode || !x || !y ||
	pkg_send(FBSERV_MSG_FBGETCURSOR, NULL, 0, remote->connection) < 0 ||
	pkg_waitfor(FBSERV_MSG_RETURN, response, sizeof(response),
	    remote->connection) < 4 * FBSERV_NET_LONG_LEN ||
	(int32_t)pkg_glong(response) != 0)
	return -1;
    *mode = (int)(int32_t)pkg_glong(response + 1 * FBSERV_NET_LONG_LEN);
    *x = (int)(int32_t)pkg_glong(response + 2 * FBSERV_NET_LONG_LEN);
    *y = (int)(int32_t)pkg_glong(response + 3 * FBSERV_NET_LONG_LEN);
    return 0;
}

int
imgstream_fb_remote_scursor(struct imgstream_fb_remote *remote, int mode,
	int x, int y)
{
    char payload[3 * FBSERV_NET_LONG_LEN] = {0};
    (void)pkg_plong(payload + 0 * FBSERV_NET_LONG_LEN, (uint32_t)mode);
    (void)pkg_plong(payload + 1 * FBSERV_NET_LONG_LEN, (uint32_t)x);
    (void)pkg_plong(payload + 2 * FBSERV_NET_LONG_LEN, (uint32_t)y);
    return remote_command(remote, FBSERV_MSG_FBSCURSOR, payload,
	sizeof(payload));
}

int
imgstream_fb_remote_setcursor(struct imgstream_fb_remote *remote,
	const unsigned char *bits, int xbits, int ybits, int xorig, int yorig)
{
    char payload[4 * FBSERV_NET_LONG_LEN] = {0};
    if (!remote || xbits < 0 || ybits < 0 ||
	(size_t)xbits > SIZE_MAX / (size_t)(ybits ? ybits : 1))
	return -1;
    size_t bytes = ((size_t)xbits * (size_t)ybits + 7) / 8;
    if (bytes > INT_MAX || (bytes && !bits))
	return -1;
    (void)pkg_plong(payload + 0 * FBSERV_NET_LONG_LEN, (uint32_t)xbits);
    (void)pkg_plong(payload + 1 * FBSERV_NET_LONG_LEN, (uint32_t)ybits);
    (void)pkg_plong(payload + 2 * FBSERV_NET_LONG_LEN, (uint32_t)xorig);
    (void)pkg_plong(payload + 3 * FBSERV_NET_LONG_LEN, (uint32_t)yorig);
    int sent = pkg_2send(FBSERV_MSG_FBSETCURSOR + FBSERV_MSG_NORETURN,
	payload, sizeof(payload), (const char *)bits, bytes,
	remote->connection);
    return sent == (int)(sizeof(payload) + bytes) ? 0 : -1;
}

int
imgstream_fb_remote_rmap(struct imgstream_fb_remote *remote,
	struct imgstream_fb_colormap *cmap)
{
    unsigned char packed[REMOTE_CMAP_BYTES] = {0};
    char response[FBSERV_NET_LONG_LEN] = {0};
    if (!remote || !cmap ||
	pkg_send(FBSERV_MSG_FBRMAP, NULL, 0, remote->connection) < 0 ||
	pkg_waitfor(FBSERV_MSG_DATA, (char *)packed, sizeof(packed),
	    remote->connection) < REMOTE_CMAP_BYTES ||
	pkg_waitfor(FBSERV_MSG_RETURN, response, sizeof(response),
	    remote->connection) < FBSERV_NET_LONG_LEN)
	return -1;
    for (int i = 0; i < 256; i++) {
	cmap->red[i] = pkg_gshort((char *)packed + 2 * (0 + i));
	cmap->green[i] = pkg_gshort((char *)packed + 2 * (256 + i));
	cmap->blue[i] = pkg_gshort((char *)packed + 2 * (512 + i));
    }
    return (int)(int32_t)pkg_glong(response);
}

int
imgstream_fb_remote_wmap(struct imgstream_fb_remote *remote,
	const struct imgstream_fb_colormap *cmap)
{
    unsigned char packed[REMOTE_CMAP_BYTES] = {0};
    if (!remote)
	return -1;
    size_t size = 0;
    if (cmap) {
	for (int i = 0; i < 256; i++) {
	    (void)pkg_pshort((char *)packed + 2 * (0 + i), cmap->red[i]);
	    (void)pkg_pshort((char *)packed + 2 * (256 + i), cmap->green[i]);
	    (void)pkg_pshort((char *)packed + 2 * (512 + i), cmap->blue[i]);
	}
	size = sizeof(packed);
    }
    return remote_command(remote, FBSERV_MSG_FBWMAP, packed, size);
}

int
imgstream_fb_remote_flush(struct imgstream_fb_remote *remote)
{
    return remote_command(remote, FBSERV_MSG_FBFLUSH, NULL, 0);
}

int
imgstream_fb_remote_poll(struct imgstream_fb_remote *remote)
{
    if (!remote || !remote->connection)
	return -1;
    return pkg_send(FBSERV_MSG_FBPOLL, NULL, 0, remote->connection) < 0 ?
	-1 : 0;
}
