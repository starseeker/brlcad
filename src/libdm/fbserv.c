/*                        F B S E R V . C
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
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
/** @addtogroup libdm */
/** @{ */
/**
 *
 * A framebuffer server object contains the attributes and
 * methods for implementing an fbserv. This code was developed
 * in large part by modifying the stand-alone version of fbserv.
 */
/** @} */

#include "common.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "bio.h"
#include "bnetwork.h"

#include "bu/str.h"
#include "raytrace.h"
#include "dm.h"
#include "./include/private.h"

/* Enable TLS server-side functions */
#define FBSERV_TLS_IMPL
#include "../fbserv/tls_wrap.h"

static struct fbserv_obj *
_fbs_conn_obj(struct pkg_conn *pcp)
{
    struct fbserv_client *fbscp;
    if (!pcp || !pcp->pkc_server_data)
	return NULL;
    fbscp = (struct fbserv_client *)pcp->pkc_server_data;
    return fbscp->fbsc_fbsp;
}

static int
fbs_has_backend(const struct fbserv_obj *fbsp)
{
    return fbserv_framebuffer_backend_installed(fbsp);
}

static int
fbs_has_framebuffer(const struct fbserv_obj *fbsp)
{
    return fbs_has_backend(fbsp) || (fbsp && fbsp->fbs_fbp);
}

static void
fbs_colormap_from_dm(struct fbserv_colormap *dst, const ColorMap *src)
{
    int i;
    if (!dst || !src)
	return;
    for (i = 0; i < 256; i++) {
	dst->red[i] = src->cm_red[i];
	dst->green[i] = src->cm_green[i];
	dst->blue[i] = src->cm_blue[i];
    }
}

static void
fbs_colormap_to_dm(ColorMap *dst, const struct fbserv_colormap *src)
{
    int i;
    if (!dst || !src)
	return;
    for (i = 0; i < 256; i++) {
	dst->cm_red[i] = src->red[i];
	dst->cm_green[i] = src->green[i];
	dst->cm_blue[i] = src->blue[i];
    }
}

static int
fbs_fb_info(struct fbserv_obj *fbsp, struct fbserv_fb_info *info)
{
    if (!fbsp || !info)
	return -1;
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_INFO))
	return fbserv_backend_info(fbsp, info);
    if (!fbsp->fbs_fbp || !fbsp->fbs_fbp->i)
	return -1;
    info->max_width = fbsp->fbs_fbp->i->if_max_width;
    info->max_height = fbsp->fbs_fbp->i->if_max_height;
    info->width = fbsp->fbs_fbp->i->if_width;
    info->height = fbsp->fbs_fbp->i->if_height;
    return 0;
}

static int
fbs_fb_clear(struct fbserv_obj *fbsp, const unsigned char rgb[3])
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_CLEAR))
	return fbserv_backend_clear(fbsp, rgb);
    return fbsp && fbsp->fbs_fbp ? fb_clear(fbsp->fbs_fbp, (unsigned char *)rgb) : -1;
}

static ssize_t
fbs_fb_read(struct fbserv_obj *fbsp, int x, int y, unsigned char *rgb, size_t count)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_READ))
	return fbserv_backend_read(fbsp, x, y, rgb, count);
    return fbsp && fbsp->fbs_fbp ? fb_read(fbsp->fbs_fbp, x, y, rgb, count) : -1;
}

static ssize_t
fbs_fb_write(struct fbserv_obj *fbsp, int x, int y, const unsigned char *rgb, size_t count)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_WRITE))
	return fbserv_backend_write(fbsp, x, y, rgb, count);
    return fbsp && fbsp->fbs_fbp ? fb_write(fbsp->fbs_fbp, x, y, rgb, count) : -1;
}

static int
fbs_fb_readrect(struct fbserv_obj *fbsp, int xmin, int ymin, int width, int height, unsigned char *rgb)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_READRECT))
	return fbserv_backend_readrect(fbsp, xmin, ymin, width, height, rgb);
    return fbsp && fbsp->fbs_fbp ? fb_readrect(fbsp->fbs_fbp, xmin, ymin, width, height, rgb) : -1;
}

static int
fbs_fb_writerect(struct fbserv_obj *fbsp, int xmin, int ymin, int width, int height, const unsigned char *rgb)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_WRITERECT))
	return fbserv_backend_writerect(fbsp, xmin, ymin, width, height, rgb);
    return fbsp && fbsp->fbs_fbp ? fb_writerect(fbsp->fbs_fbp, xmin, ymin, width, height, rgb) : -1;
}

static int
fbs_fb_bwreadrect(struct fbserv_obj *fbsp, int xmin, int ymin, int width, int height, unsigned char *bw)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_BWREADRECT))
	return fbserv_backend_bwreadrect(fbsp, xmin, ymin, width, height, bw);
    return fbsp && fbsp->fbs_fbp ? fb_bwreadrect(fbsp->fbs_fbp, xmin, ymin, width, height, bw) : -1;
}

static int
fbs_fb_bwwriterect(struct fbserv_obj *fbsp, int xmin, int ymin, int width, int height, const unsigned char *bw)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_BWWRITERECT))
	return fbserv_backend_bwwriterect(fbsp, xmin, ymin, width, height, bw);
    return fbsp && fbsp->fbs_fbp ? fb_bwwriterect(fbsp->fbs_fbp, xmin, ymin, width, height, bw) : -1;
}

static int
fbs_fb_cursor(struct fbserv_obj *fbsp, int mode, int x, int y)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_CURSOR))
	return fbserv_backend_cursor(fbsp, mode, x, y);
    return fbsp && fbsp->fbs_fbp ? fb_cursor(fbsp->fbs_fbp, mode, x, y) : -1;
}

static int
fbs_fb_getcursor(struct fbserv_obj *fbsp, int *mode, int *x, int *y)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_GETCURSOR))
	return fbserv_backend_getcursor(fbsp, mode, x, y);
    return fbsp && fbsp->fbs_fbp ? fb_getcursor(fbsp->fbs_fbp, mode, x, y) : -1;
}

static int
fbs_fb_setcursor(struct fbserv_obj *fbsp, const unsigned char *bits, int xbits, int ybits, int xorig, int yorig)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_SETCURSOR))
	return fbserv_backend_setcursor(fbsp, bits, xbits, ybits, xorig, yorig);
    return fbsp && fbsp->fbs_fbp ? fb_setcursor(fbsp->fbs_fbp, bits, xbits, ybits, xorig, yorig) : -1;
}

static int
fbs_fb_scursor(struct fbserv_obj *fbsp, int mode, int x, int y)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_SCURSOR))
	return fbserv_backend_scursor(fbsp, mode, x, y);
    return fbsp && fbsp->fbs_fbp ? fb_scursor(fbsp->fbs_fbp, mode, x, y) : -1;
}

static int
fbs_fb_window(struct fbserv_obj *fbsp, int xcenter, int ycenter)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_WINDOW))
	return fbserv_backend_window(fbsp, xcenter, ycenter);
    return fbsp && fbsp->fbs_fbp ? fb_window(fbsp->fbs_fbp, xcenter, ycenter) : -1;
}

static int
fbs_fb_zoom(struct fbserv_obj *fbsp, int xzoom, int yzoom)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_ZOOM))
	return fbserv_backend_zoom(fbsp, xzoom, yzoom);
    return fbsp && fbsp->fbs_fbp ? fb_zoom(fbsp->fbs_fbp, xzoom, yzoom) : -1;
}

static int
fbs_fb_view(struct fbserv_obj *fbsp, int xcenter, int ycenter, int xzoom, int yzoom)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_VIEW))
	return fbserv_backend_view(fbsp, xcenter, ycenter, xzoom, yzoom);
    return fbsp && fbsp->fbs_fbp ? fb_view(fbsp->fbs_fbp, xcenter, ycenter, xzoom, yzoom) : -1;
}

static int
fbs_fb_getview(struct fbserv_obj *fbsp, int *xcenter, int *ycenter, int *xzoom, int *yzoom)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_GETVIEW))
	return fbserv_backend_getview(fbsp, xcenter, ycenter, xzoom, yzoom);
    return fbsp && fbsp->fbs_fbp ? fb_getview(fbsp->fbs_fbp, xcenter, ycenter, xzoom, yzoom) : -1;
}

static int
fbs_fb_rmap(struct fbserv_obj *fbsp, struct fbserv_colormap *cmap)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_RMAP))
	return fbserv_backend_rmap(fbsp, cmap);
    if (fbsp && fbsp->fbs_fbp) {
	int ret;
	ColorMap dm_map;
	ret = fb_rmap(fbsp->fbs_fbp, &dm_map);
	fbs_colormap_from_dm(cmap, &dm_map);
	return ret;
    }
    return -1;
}

static int
fbs_fb_wmap(struct fbserv_obj *fbsp, const struct fbserv_colormap *cmap)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_WMAP))
	return fbserv_backend_wmap(fbsp, cmap);
    if (fbsp && fbsp->fbs_fbp) {
	if (!cmap)
	    return fb_wmap(fbsp->fbs_fbp, COLORMAP_NULL);
	ColorMap dm_map;
	fbs_colormap_to_dm(&dm_map, cmap);
	return fb_wmap(fbsp->fbs_fbp, &dm_map);
    }
    return -1;
}

static int
fbs_fb_flush(struct fbserv_obj *fbsp)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_FLUSH))
	return fbserv_backend_flush(fbsp);
    return fbsp && fbsp->fbs_fbp ? fb_flush(fbsp->fbs_fbp) : -1;
}

static int
fbs_fb_poll(struct fbserv_obj *fbsp)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_POLL))
	return fbserv_backend_poll(fbsp);
    return fbsp && fbsp->fbs_fbp ? fb_poll(fbsp->fbs_fbp) : -1;
}

static int
fbs_fb_help(struct fbserv_obj *fbsp)
{
    if (fbserv_backend_op_installed(fbsp, FBSERV_BACKEND_OP_HELP))
	return fbserv_backend_help(fbsp);
    return fbsp && fbsp->fbs_fbp ? fb_help(fbsp->fbs_fbp) : 0;
}

static void
drop_client(struct fbserv_obj *fbsp, int sub)
{
    if (fbsp->fbs_clients[sub].fbsc_pkg != PKC_NULL) {
	pkg_close(fbsp->fbs_clients[sub].fbsc_pkg);
	fbsp->fbs_clients[sub].fbsc_pkg = PKC_NULL;
    }

    if (fbsp->fbs_clients[sub].fbsc_fd != 0) {
	/* Phase A1 (ert reliability): release the framebuffer-dimension
	 * lock now that this client is no longer streaming. */
	if (fbsp->fbs_fbp && fbsp->fbs_fbp->i &&
	    fbsp->fbs_fbp->i->if_active_clients > 0)
	    fbsp->fbs_fbp->i->if_active_clients--;
	/* Phase A2 (ert reliability): when the lock has just been
	 * released (count == 0), finalize any window resize that was
	 * deferred while a client was streaming.  swrast/qtgl
	 * configureWindow only updated vp_width/vp_height while the
	 * lock was held, leaving if_width/if_height (and if_mem) at
	 * their old values.  Re-driving fb_configure_window with the
	 * dm's current dimensions takes the unlocked branch and
	 * properly reallocs if_mem, resets if_xcenter/if_ycenter,
	 * and runs fb_clipper.  Must run BEFORE the close handler
	 * (which may deleteLater() the toolkit channel object) and
	 * AFTER the decrement so the unlocked branch is taken. */
	if (fbsp->fbs_fbp && fbsp->fbs_fbp->i &&
	    fbsp->fbs_fbp->i->if_active_clients == 0 &&
	    fbsp->fbs_fbp->i->dmp) {
	    struct dm *dmp = fbsp->fbs_fbp->i->dmp;
	    int dw = dm_get_width(dmp);
	    int dh = dm_get_height(dmp);
	    bu_log("drop_client: sub=%d dm=(%d,%d) if=(%d,%d)%s\n",
		   sub, dw, dh,
		   fbsp->fbs_fbp->i->if_width, fbsp->fbs_fbp->i->if_height,
		   (dw > 0 && dh > 0 && (dw != fbsp->fbs_fbp->i->if_width || dh != fbsp->fbs_fbp->i->if_height))
		   ? " => calling fb_configure_window" : " => no resize needed");
	    if (dw > 0 && dh > 0 &&
		(dw != fbsp->fbs_fbp->i->if_width ||
		 dh != fbsp->fbs_fbp->i->if_height)) {
		fb_configure_window(fbsp->fbs_fbp, dw, dh);
	    }
	}
	/* Use the IPC-specific close handler if the client was opened via IPC
	 * and the caller registered one; otherwise fall back to the generic
	 * TCP close handler. */
	if (fbsp->fbs_clients[sub].fbsc_is_ipc && fbsp->fbs_close_ipc_client_handler)
	    (*fbsp->fbs_close_ipc_client_handler)(fbsp, sub);
	else
	    (*fbsp->fbs_close_client_handler)(fbsp, sub);
	fbsp->fbs_clients[sub].fbsc_fd = 0;
    }
    fbsp->fbs_clients[sub].fbsc_auth_ok = 0;
    fbsp->fbs_clients[sub].fbsc_pending_drop = 0;
    fbsp->fbs_clients[sub].fbsc_is_ipc = 0;
}


/**
 * Guard for data-op handlers: check auth + non-NULL framebuffer.
 * Returns 0 on success.  On failure sends a -1 reply, schedules a
 * deferred drop when auth fails, frees buf, and returns -1.
 */
static int
fbs_data_guard(struct pkg_conn *pcp, char *buf)
{
    char erbuf[NET_LONG_LEN+1];
    struct fbserv_client *fbscp;
    struct fbserv_obj *fbsp;

    if (!pcp || !pcp->pkc_server_data) {
	if (buf) (void)free(buf);
	return -1;
    }
    fbscp = (struct fbserv_client *)pcp->pkc_server_data;
    fbsp  = fbscp->fbsc_fbsp;

    if (fbsp && fbsp->fbs_require_auth && !fbscp->fbsc_auth_ok) {
	bu_log("fbserv: unauthenticated data request (type %d) rejected\n",
	       pcp->pkc_type);
	(void)pkg_plong(erbuf, -1);
	pkg_send(MSG_RETURN, erbuf, NET_LONG_LEN, pcp);
	fbscp->fbsc_pending_drop = 1;
	if (buf) (void)free(buf);
	return -1;
    }

    if (!fbs_has_framebuffer(fbsp)) {
	bu_log("fbserv: data request (type %d) with null framebuffer\n",
	       pcp->pkc_type);
	(void)pkg_plong(erbuf, -1);
	pkg_send(MSG_RETURN, erbuf, NET_LONG_LEN, pcp);
	if (buf) (void)free(buf);
	return -1;
    }

    return 0;
}


/*
 * This is where we go for message types we don't understand.
 */
static void
fbs_rfbunknown(struct pkg_conn *pcp, char *buf)
{
    bu_log("fbserv: unable to handle message type %d\n", pcp->pkc_type);
    if (buf) {
	(void)free(buf);
    }
}


/******** Here's where the hooks lead *********/

/**
 * MSG_FBAUTH — session token authentication.
 *
 * The client sends a FBSERV_AUTH_TOKEN_LEN-byte hex token string.
 * If it matches the server's session token the connection is marked
 * authenticated and all subsequent requests are allowed.  If the token
 * is wrong the connection is closed immediately.
 *
 * Old clients that do not send MSG_FBAUTH are still accepted unless
 * the server is running in strict mode (fbs_require_auth != 0).
 */
static void
fbs_rfbauth(struct pkg_conn *pcp, char *buf)
{
    struct fbserv_client *fbscp;
    struct fbserv_obj *fbsp;
    char provided[FBSERV_AUTH_TOKEN_LEN + 1] = {0};

    if (!pcp || !pcp->pkc_server_data) {
	if (buf) (void)free(buf);
	return;
    }

    fbscp = (struct fbserv_client *)pcp->pkc_server_data;
    fbsp = fbscp->fbsc_fbsp;

    if (!fbsp || fbsp->fbs_auth_token[0] == '\0') {
	/* No token configured — mark auth as satisfied */
	fbscp->fbsc_auth_ok = 1;
	if (buf) (void)free(buf);
	return;
    }

    if (buf && pcp->pkc_len >= FBSERV_AUTH_TOKEN_LEN) {
	bu_strlcpy(provided, buf, sizeof(provided));
    }

    if (fbserv_verify_token(provided, fbsp->fbs_auth_token)) {
	fbscp->fbsc_auth_ok = 1;
    } else {
	bu_log("fbserv: MSG_FBAUTH token mismatch from client — dropping\n");
	/* Use deferred drop: pkg_process still holds a reference to pcp.
	 * Setting fbsc_pending_drop causes fbs_existing_client_handler to
	 * call drop_client() after pkg_process() returns. */
	fbscp->fbsc_pending_drop = 1;
    }

    if (buf) (void)free(buf);
}


static void
fbs_rfbopen(struct pkg_conn *pcp, char *buf)
{
    struct fbserv_client *fbscp;
    struct fbserv_obj *fbsp;
    struct fbserv_fb_info info;
    char rbuf[5*NET_LONG_LEN+1] = {0};
    int want;

    /* Auth check: if the server requires authentication and this
     * connection has not sent a valid MSG_FBAUTH, reject with an
     * error return code and close the connection. */
    if (pcp && pcp->pkc_server_data) {
	fbscp = (struct fbserv_client *)pcp->pkc_server_data;
	fbsp  = fbscp->fbsc_fbsp;
	if (fbsp && fbsp->fbs_require_auth && !fbscp->fbsc_auth_ok) {
	    bu_log("fbserv: unauthenticated MSG_FBOPEN rejected (strict mode)\n");
	    (void)pkg_plong(&rbuf[0*NET_LONG_LEN], -1);	/* failure */
	    (void)pkg_plong(&rbuf[1*NET_LONG_LEN], 0);
	    (void)pkg_plong(&rbuf[2*NET_LONG_LEN], 0);
	    (void)pkg_plong(&rbuf[3*NET_LONG_LEN], 0);
	    (void)pkg_plong(&rbuf[4*NET_LONG_LEN], 0);
	    pkg_send(MSG_RETURN, rbuf, 5*NET_LONG_LEN, pcp);
	    /* Deferred drop: pkg_process still holds pcp; close after return */
	    fbscp->fbsc_pending_drop = 1;
	    if (buf) (void)free(buf);
	    return;
	}
    }
    fbsp = _fbs_conn_obj(pcp);
    if (fbs_fb_info(fbsp, &info) != 0) {
	(void)pkg_plong(&rbuf[0*NET_LONG_LEN], -1);
	(void)pkg_plong(&rbuf[1*NET_LONG_LEN], 0);
	(void)pkg_plong(&rbuf[2*NET_LONG_LEN], 0);
	(void)pkg_plong(&rbuf[3*NET_LONG_LEN], 0);
	(void)pkg_plong(&rbuf[4*NET_LONG_LEN], 0);
	pkg_send(MSG_RETURN, rbuf, 5*NET_LONG_LEN, pcp);
	if (buf) (void)free(buf);
	return;
    }

    /* Don't really open a new framebuffer --- use existing one */
    (void)pkg_plong(&rbuf[0*NET_LONG_LEN], 0);	/* ret */
    (void)pkg_plong(&rbuf[1*NET_LONG_LEN], info.max_width);
    (void)pkg_plong(&rbuf[2*NET_LONG_LEN], info.max_height);
    (void)pkg_plong(&rbuf[3*NET_LONG_LEN], info.width);
    (void)pkg_plong(&rbuf[4*NET_LONG_LEN], info.height);

    want = 5*NET_LONG_LEN;
    if (pkg_send(MSG_RETURN, rbuf, want, pcp) != want)
	bu_log("pkg_send fb_open reply\n");

    if (buf) {
	(void)free(buf);
    }
}


void
fbs_rfbclose(struct pkg_conn *pcp, char *buf)
{
    struct fbserv_obj *fbsp;
    char rbuf[NET_LONG_LEN+1] = {0};

    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    /*
     * We are playing FB server so we don't really close the frame
     * buffer.  We should flush output however.
     */
    (void)fbs_fb_flush(fbsp);
    (void)pkg_plong(&rbuf[0], 0);		/* return success */

    /* Don't check for errors, linger mode or other events may have
     * already closed down all the file descriptors.  If communication
     * has broken, other end will know we are gone.
     */
    (void)pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);

    if (buf) {
	(void)free(buf);
    }
}


void
fbs_rfbfree(struct pkg_conn *pcp, char *buf)
{
    char rbuf[NET_LONG_LEN+1] = {0};

    if (fbs_data_guard(pcp, buf) < 0) return;
    /* Don't really free framebuffer */
    if (pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp) != NET_LONG_LEN)
	bu_log("pkg_send fb_free reply\n");

    if (buf) {
	(void)free(buf);
    }
}


void
fbs_rfbclear(struct pkg_conn *pcp, char *buf)
{
    struct fbserv_obj *fbsp;
    RGBpixel bg;
    char rbuf[NET_LONG_LEN+1] = {0};

    if (!buf) {
	bu_log("fbs_rfbclear: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    bg[RED] = buf[0];
    bg[GRN] = buf[1];
    bg[BLU] = buf[2];

    (void)pkg_plong(rbuf, fbs_fb_clear(fbsp, bg));
    pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);

    (void)free(buf);
}


void
fbs_rfbread(struct pkg_conn *pcp, char *buf)
{
    int x, y;
    size_t num;
    int ret;
    static unsigned char *scanbuf = NULL;
    static size_t buflen = 0;
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbread: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    x = pkg_glong(&buf[0*NET_LONG_LEN]);
    y = pkg_glong(&buf[1*NET_LONG_LEN]);
    num = (size_t)pkg_glong(&buf[2*NET_LONG_LEN]);

    if (num*sizeof(RGBpixel) > buflen) {
	if (scanbuf != NULL)
	    free((char *)scanbuf);
	buflen = num*sizeof(RGBpixel);
	if (buflen < 1024*sizeof(RGBpixel))
	    buflen = 1024*sizeof(RGBpixel);
	if ((scanbuf = (unsigned char *)malloc(buflen)) == NULL) {
	    fb_log("fb_read: malloc failed!");
	    (void)free(buf);
	    buflen = 0;
	    return;
	}
    }

    ret = fbs_fb_read(fbsp, x, y, scanbuf, num);
    if (ret < 0) ret = 0;		/* map error indications */
    /* sending a 0-length package indicates error */
    pkg_send(MSG_RETURN, (char *)scanbuf, ret*sizeof(RGBpixel), pcp);
    (void)free(buf);
}


void
fbs_rfbwrite(struct pkg_conn *pcp, char *buf)
{
    long x, y, num;
    char rbuf[NET_LONG_LEN+1] = {0};
    int ret;
    int type;
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbwrite: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    x = pkg_glong(&buf[0*NET_LONG_LEN]);
    y = pkg_glong(&buf[1*NET_LONG_LEN]);
    num = pkg_glong(&buf[2*NET_LONG_LEN]);
    type = pcp->pkc_type;
    ret = fbs_fb_write(fbsp, x, y, (unsigned char *)&buf[3*NET_LONG_LEN], (size_t)num);

    if (type < MSG_NORETURN) {
	(void)pkg_plong(&rbuf[0*NET_LONG_LEN], ret);
	pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    }
    (void)free(buf);
}


void
fbs_rfbreadrect(struct pkg_conn *pcp, char *buf)
{
    int xmin, ymin;
    int width, height;
    size_t num;
    int ret;
    static unsigned char *scanbuf = NULL;
    static size_t buflen = 0;
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbreadrect: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    xmin = pkg_glong(&buf[0*NET_LONG_LEN]);
    ymin = pkg_glong(&buf[1*NET_LONG_LEN]);
    width = pkg_glong(&buf[2*NET_LONG_LEN]);
    height = pkg_glong(&buf[3*NET_LONG_LEN]);
    num = width * height;

    if (num*sizeof(RGBpixel) > buflen) {
	if (scanbuf != NULL)
	    free((char *)scanbuf);
	buflen = num*sizeof(RGBpixel);
	if (buflen < 1024*sizeof(RGBpixel))
	    buflen = 1024*sizeof(RGBpixel);
	if ((scanbuf = (unsigned char *)malloc(buflen)) == NULL) {
	    fb_log("fb_read: malloc failed!");
	    (void)free(buf);
	    buflen = 0;
	    return;
	}
    }

    ret = fbs_fb_readrect(fbsp, xmin, ymin, width, height, scanbuf);
    if (ret < 0) ret = 0;		/* map error indications */
    /* sending a 0-length package indicates error */
    pkg_send(MSG_RETURN, (char *)scanbuf, ret*sizeof(RGBpixel), pcp);
    (void)free(buf);
}


void
fbs_rfbwriterect(struct pkg_conn *pcp, char *buf)
{
    int x, y;
    int width, height;
    char rbuf[NET_LONG_LEN+1] = {0};
    int ret;
    int type;
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbwriterect: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    x = pkg_glong(&buf[0*NET_LONG_LEN]);
    y = pkg_glong(&buf[1*NET_LONG_LEN]);
    width = pkg_glong(&buf[2*NET_LONG_LEN]);
    height = pkg_glong(&buf[3*NET_LONG_LEN]);

    type = pcp->pkc_type;
    ret = fbs_fb_writerect(fbsp, x, y, width, height,
	    (unsigned char *)&buf[4*NET_LONG_LEN]);

    if (type < MSG_NORETURN) {
	(void)pkg_plong(&rbuf[0*NET_LONG_LEN], ret);
	pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    }
    (void)free(buf);
}


void
fbs_rfbbwreadrect(struct pkg_conn *pcp, char *buf)
{
    int xmin, ymin;
    int width, height;
    int num;
    int ret;
    static unsigned char *scanbuf = NULL;
    static int buflen = 0;
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbbwreadrect: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    xmin = pkg_glong(&buf[0*NET_LONG_LEN]);
    ymin = pkg_glong(&buf[1*NET_LONG_LEN]);
    width = pkg_glong(&buf[2*NET_LONG_LEN]);
    height = pkg_glong(&buf[3*NET_LONG_LEN]);
    num = width * height;

    if (num > buflen) {
	if (scanbuf != NULL)
	    free((char *)scanbuf);
	buflen = num;
	if (buflen < 1024)
	    buflen = 1024;
	if ((scanbuf = (unsigned char *)malloc(buflen)) == NULL) {
	    fb_log("fbs_rfbbwreadrect: malloc failed!");
	    (void)free(buf);
	    buflen = 0;
	    return;
	}
    }

    ret = fbs_fb_bwreadrect(fbsp, xmin, ymin, width, height, scanbuf);
    if (ret < 0) ret = 0;		/* map error indications */
    /* sending a 0-length package indicates error */
    pkg_send(MSG_RETURN, (char *)scanbuf, ret, pcp);
    (void)free(buf);
}


void
fbs_rfbbwwriterect(struct pkg_conn *pcp, char *buf)
{
    int x, y;
    int width, height;
    char rbuf[NET_LONG_LEN+1] = {0};
    int ret;
    int type;
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbbwwriterect: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    x = pkg_glong(&buf[0*NET_LONG_LEN]);
    y = pkg_glong(&buf[1*NET_LONG_LEN]);
    width = pkg_glong(&buf[2*NET_LONG_LEN]);
    height = pkg_glong(&buf[3*NET_LONG_LEN]);

    type = pcp->pkc_type;
    ret = fbs_fb_bwwriterect(fbsp, x, y, width, height,
	    (unsigned char *)&buf[4*NET_LONG_LEN]);

    if (type < MSG_NORETURN) {
	(void)pkg_plong(&rbuf[0*NET_LONG_LEN], ret);
	pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    }
    (void)free(buf);
}


void
fbs_rfbcursor(struct pkg_conn *pcp, char *buf)
{
    int mode, x, y;
    char rbuf[NET_LONG_LEN+1] = {0};
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbcursor: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    mode = pkg_glong(&buf[0*NET_LONG_LEN]);
    x = pkg_glong(&buf[1*NET_LONG_LEN]);
    y = pkg_glong(&buf[2*NET_LONG_LEN]);

    (void)pkg_plong(&rbuf[0], fbs_fb_cursor(fbsp, mode, x, y));
    pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    (void)free(buf);
}


void
fbs_rfbgetcursor(struct pkg_conn *pcp, char *buf)
{
    int ret;
    int mode, x, y;
    char rbuf[4*NET_LONG_LEN+1] = {0};
    struct fbserv_obj *fbsp;

    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);
    ret = fbs_fb_getcursor(fbsp, &mode, &x, &y);
    (void)pkg_plong(&rbuf[0*NET_LONG_LEN], ret);
    (void)pkg_plong(&rbuf[1*NET_LONG_LEN], mode);
    (void)pkg_plong(&rbuf[2*NET_LONG_LEN], x);
    (void)pkg_plong(&rbuf[3*NET_LONG_LEN], y);
    pkg_send(MSG_RETURN, rbuf, 4*NET_LONG_LEN, pcp);

    if (buf) {
	(void)free(buf);
    }
}


void
fbs_rfbsetcursor(struct pkg_conn *pcp, char *buf)
{
    char rbuf[NET_LONG_LEN+1] = {0};
    int ret;
    int xbits, ybits;
    int xorig, yorig;
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfsetcursor: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    xbits = pkg_glong(&buf[0*NET_LONG_LEN]);
    ybits = pkg_glong(&buf[1*NET_LONG_LEN]);
    xorig = pkg_glong(&buf[2*NET_LONG_LEN]);
    yorig = pkg_glong(&buf[3*NET_LONG_LEN]);

    ret = fbs_fb_setcursor(fbsp, (unsigned char *)&buf[4*NET_LONG_LEN],
	    xbits, ybits, xorig, yorig);

    if (pcp->pkc_type < MSG_NORETURN) {
	(void)pkg_plong(&rbuf[0*NET_LONG_LEN], ret);
	pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    }
    (void)free(buf);
}


/*OLD*/
void
fbs_rfbscursor(struct pkg_conn *pcp, char *buf)
{
    int mode, x, y;
    char rbuf[NET_LONG_LEN+1] = {0};
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbscursor: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    mode = pkg_glong(&buf[0*NET_LONG_LEN]);
    x = pkg_glong(&buf[1*NET_LONG_LEN]);
    y = pkg_glong(&buf[2*NET_LONG_LEN]);

    (void)pkg_plong(&rbuf[0], fbs_fb_scursor(fbsp, mode, x, y));
    pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    (void)free(buf);
}


/*OLD*/
void
fbs_rfbwindow(struct pkg_conn *pcp, char *buf)
{
    int x, y;
    char rbuf[NET_LONG_LEN+1] = {0};
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbwindow: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    x = pkg_glong(&buf[0*NET_LONG_LEN]);
    y = pkg_glong(&buf[1*NET_LONG_LEN]);

    (void)pkg_plong(&rbuf[0], fbs_fb_window(fbsp, x, y));
    pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);

    (void)free(buf);
}


/*OLD*/
void
fbs_rfbzoom(struct pkg_conn *pcp, char *buf)
{
    int x, y;
    char rbuf[NET_LONG_LEN+1] = {0};
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbzoom: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    x = pkg_glong(&buf[0*NET_LONG_LEN]);
    y = pkg_glong(&buf[1*NET_LONG_LEN]);

    (void)pkg_plong(&rbuf[0], fbs_fb_zoom(fbsp, x, y));
    pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    (void)free(buf);
}


void
fbs_rfbview(struct pkg_conn *pcp, char *buf)
{
    int ret;
    int xcenter, ycenter, xzoom, yzoom;
    char rbuf[NET_LONG_LEN+1] = {0};
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbv: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    xcenter = pkg_glong(&buf[0*NET_LONG_LEN]);
    ycenter = pkg_glong(&buf[1*NET_LONG_LEN]);
    xzoom = pkg_glong(&buf[2*NET_LONG_LEN]);
    yzoom = pkg_glong(&buf[3*NET_LONG_LEN]);

    ret = fbs_fb_view(fbsp, xcenter, ycenter, xzoom, yzoom);
    (void)pkg_plong(&rbuf[0], ret);
    pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    (void)free(buf);
}


void
fbs_rfbgetview(struct pkg_conn *pcp, char *buf)
{
    int ret;
    int xcenter, ycenter, xzoom, yzoom;
    char rbuf[5*NET_LONG_LEN+1] = {0};
    struct fbserv_obj *fbsp;

    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);
    ret = fbs_fb_getview(fbsp, &xcenter, &ycenter, &xzoom, &yzoom);
    (void)pkg_plong(&rbuf[0*NET_LONG_LEN], ret);
    (void)pkg_plong(&rbuf[1*NET_LONG_LEN], xcenter);
    (void)pkg_plong(&rbuf[2*NET_LONG_LEN], ycenter);
    (void)pkg_plong(&rbuf[3*NET_LONG_LEN], xzoom);
    (void)pkg_plong(&rbuf[4*NET_LONG_LEN], yzoom);
    pkg_send(MSG_RETURN, rbuf, 5*NET_LONG_LEN, pcp);

    if (buf) {
	(void)free(buf);
    }
}


void
fbs_rfbrmap(struct pkg_conn *pcp, char *buf)
{
    register int i;
    char rbuf[NET_LONG_LEN+1] = {0};
    struct fbserv_colormap map;
    unsigned char cm[256*2*3];
    struct fbserv_obj *fbsp;

    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);
    (void)pkg_plong(&rbuf[0*NET_LONG_LEN], fbs_fb_rmap(fbsp, &map));
    for (i = 0; i < 256; i++) {
	(void)pkg_pshort((char *)(cm+2*(0+i)), map.red[i]);
	(void)pkg_pshort((char *)(cm+2*(256+i)), map.green[i]);
	(void)pkg_pshort((char *)(cm+2*(512+i)), map.blue[i]);
    }
    pkg_send(MSG_DATA, (char *)cm, sizeof(cm), pcp);
    pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);

    if (buf) {
	(void)free(buf);
    }
}


/*
 * Accept a color map sent by the client, and write it to the
 * framebuffer.  Network format is to send each entry as a network
 * (IBM) order 2-byte short, 256 red shorts, followed by 256 green and
 * 256 blue, for a total of 3*256*2 bytes.
 */
void
fbs_rfbwmap(struct pkg_conn *pcp, char *buf)
{
    int i;
    char rbuf[NET_LONG_LEN+1] = {0};
    long ret;
    struct fbserv_colormap map;
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbwmap: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    if (pcp->pkc_len == 0) {
	ret = fbs_fb_wmap(fbsp, NULL);
    } else {
	for (i = 0; i < 256; i++) {
	    map.red[i] = pkg_gshort(buf+2*(0+i));
	    map.green[i] = pkg_gshort(buf+2*(256+i));
	    map.blue[i] = pkg_gshort(buf+2*(512+i));
	}
	ret = fbs_fb_wmap(fbsp, &map);
    }
    (void)pkg_plong(&rbuf[0], ret);
    pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    (void)free(buf);
}


void
fbs_rfbflush(struct pkg_conn *pcp, char *buf)
{
    int ret;
    char rbuf[NET_LONG_LEN+1] = {0};
    struct fbserv_obj *fbsp;

    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);
    ret = fbs_fb_flush(fbsp);

    if (pcp->pkc_type < MSG_NORETURN) {
	(void)pkg_plong(rbuf, ret);
	pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    }

    if (buf) {
	(void)free(buf);
    }
}


void
fbs_rfbpoll(struct pkg_conn *pcp, char *buf)
{
    struct fbserv_obj *fbsp;

    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    (void)fbs_fb_poll(fbsp);
    if (buf) {
	(void)free(buf);
    }
}


/*
 * At one time at least we couldn't send a zero length PKG message
 * back and forth, so we receive a dummy long here.
 */
void
fbs_rfbhelp(struct pkg_conn *pcp, char *buf)
{
    long ret;
    char rbuf[NET_LONG_LEN+1] = {0};
    struct fbserv_obj *fbsp;

    if (!buf) {
	bu_log("fbs_rfbhelp: null buffer\n");
	return;
    }
    if (fbs_data_guard(pcp, buf) < 0) return;
    fbsp = _fbs_conn_obj(pcp);

    (void)pkg_glong(&buf[0*NET_LONG_LEN]);

    ret = fbs_fb_help(fbsp);
    (void)pkg_plong(&rbuf[0], ret);
    pkg_send(MSG_RETURN, rbuf, NET_LONG_LEN, pcp);
    (void)free(buf);
}


/**
 * Initialise fbsp->fbs_auth_token for session authentication.
 *
 * If the FBSERV_TOKEN environment variable is already set to a valid
 * 64-hex-char token, that value is used directly.  This lets the
 * hosting application (MGED, qged, etc.) pre-supply a known token and
 * pass the same value to child processes (e.g. via setenv() before
 * fork/exec of rt).  Token authentication works regardless of whether
 * TLS is enabled — it provides session isolation even on plain TCP.
 *
 * If FBSERV_TOKEN is not set or is the wrong length, a fresh random
 * token is generated.
 *
 * Should be called before fbs_open() so the token is ready for the
 * first connecting client.  Returns a pointer to fbsp->fbs_auth_token.
 */
const char *
fbs_generate_token(struct fbserv_obj *fbsp)
{
    return fbserv_obj_generate_token(fbsp);
}


int
fbs_init(struct fbserv_obj *fbsp)
{
    return fbserv_obj_init(fbsp);
}


int
fbs_set_backend(struct fbserv_obj *fbsp,
	const struct fbserv_fb_ops *ops,
	void *ctx)
{
    return fbserv_set_backend(fbsp, ops, ctx);
}


void
fbs_clear_backend(struct fbserv_obj *fbsp)
{
    fbserv_clear_backend(fbsp);
}


void
fbs_set_transport(struct fbserv_obj *fbsp, const struct fbserv_transport_ops *ops)
{
    fbserv_set_transport(fbsp, ops);
}


void
fbs_clear_transport(struct fbserv_obj *fbsp)
{
    fbserv_clear_transport(fbsp);
}


int
fbs_can_open_ipc(const struct fbserv_obj *fbsp)
{
    return fbserv_can_open_ipc(fbsp);
}


int
fbs_can_open_network(const struct fbserv_obj *fbsp)
{
    return fbserv_can_open_network(fbsp);
}


int
fbs_can_close(const struct fbserv_obj *fbsp)
{
    return fbserv_can_close(fbsp);
}


int
fbs_listener_port(const struct fbserv_obj *fbsp)
{
    return fbserv_listener_port(fbsp);
}


int
fbs_listener_fd(const struct fbserv_obj *fbsp)
{
    return fbserv_listener_fd(fbsp);
}


void
fbs_set_listener_fd(struct fbserv_obj *fbsp, int fd)
{
    fbserv_set_listener_fd(fbsp, fd);
}


void *
fbs_listener_channel(const struct fbserv_obj *fbsp)
{
    return fbserv_listener_channel(fbsp);
}


void
fbs_set_listener_channel(struct fbserv_obj *fbsp, void *chan)
{
    fbserv_set_listener_channel(fbsp, chan);
}


void *
fbs_listener_handler_data(struct fbserv_obj *fbsp)
{
    return fbserv_listener_handler_data(fbsp);
}


struct fbserv_obj *
fbs_listener_owner(void *listener_data)
{
    return fbserv_listener_owner(listener_data);
}


int
fbs_listener_data_fd(const void *listener_data)
{
    return fbserv_listener_data_fd(listener_data);
}


struct pkg_listener *
fbs_listener_data_pkg_listener(void *listener_data)
{
    return fbserv_listener_data_pkg_listener(listener_data);
}


void
fbs_set_listener_pkg_listener(struct fbserv_obj *fbsp,
	struct pkg_listener *listener)
{
    fbserv_set_listener_pkg_listener(fbsp, listener);
}


int
fbs_client_active(const struct fbserv_obj *fbsp, int sub)
{
    return fbserv_client_active(fbsp, sub);
}


int
fbs_client_fd(const struct fbserv_obj *fbsp, int sub)
{
    return fbserv_client_fd(fbsp, sub);
}


struct pkg_conn *
fbs_client_pkg(const struct fbserv_obj *fbsp, int sub)
{
    struct pkg_conn *pcp = fbserv_client_pkg(fbsp, sub);
    return pcp ? pcp : PKC_NULL;
}


struct pkg_conn *
fbs_client_data_pkg(void *client_data)
{
    struct pkg_conn *pcp = fbserv_client_data_pkg(client_data);
    return pcp ? pcp : PKC_NULL;
}


int
fbs_client_data_fd(const void *client_data)
{
    return fbserv_client_data_fd(client_data);
}


void *
fbs_client_channel(const struct fbserv_obj *fbsp, int sub)
{
    return fbserv_client_channel(fbsp, sub);
}


void
fbs_set_client_channel(struct fbserv_obj *fbsp, int sub, void *chan)
{
    fbserv_set_client_channel(fbsp, sub, chan);
}


void *
fbs_client_handler(const struct fbserv_obj *fbsp, int sub)
{
    return fbserv_client_handler(fbsp, sub);
}


void
fbs_set_client_handler(struct fbserv_obj *fbsp, int sub, void *handler)
{
    fbserv_set_client_handler(fbsp, sub, handler);
}


void *
fbs_client_handler_data(struct fbserv_obj *fbsp, int sub)
{
    return fbserv_client_handler_data(fbsp, sub);
}


void
fbs_set_client_data_channel(void *client_data, void *chan)
{
    fbserv_set_client_data_channel(client_data, chan);
}


void
fbs_set_legacy_framebuffer(struct fbserv_obj *fbsp, struct fb *fbp)
{
    fbserv_set_legacy_framebuffer(fbsp, fbp);
}


struct fb *
fbs_legacy_framebuffer(struct fbserv_obj *fbsp)
{
    return fbserv_legacy_framebuffer(fbsp);
}


int
fbs_framebuffer_backend_installed(const struct fbserv_obj *fbsp)
{
    return fbserv_framebuffer_backend_installed(fbsp);
}


int
fbs_framebuffer_info(struct fbserv_obj *fbsp, struct fbserv_fb_info *info)
{
    return fbs_fb_info(fbsp, info);
}


int
fbs_framebuffer_writerect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	const unsigned char *rgb)
{
    return fbs_fb_writerect(fbsp, xmin, ymin, width, height, rgb);
}


int
fbs_framebuffer_view(struct fbserv_obj *fbsp,
	int xcenter,
	int ycenter,
	int xzoom,
	int yzoom)
{
    return fbs_fb_view(fbsp, xcenter, ycenter, xzoom, yzoom);
}


int
fbs_framebuffer_cursor(struct fbserv_obj *fbsp, int mode, int x, int y)
{
    return fbs_fb_cursor(fbsp, mode, x, y);
}


int
fbs_framebuffer_flush(struct fbserv_obj *fbsp)
{
    return fbs_fb_flush(fbsp);
}


int
fbs_framebuffer_poll(struct fbserv_obj *fbsp)
{
    return fbs_fb_poll(fbsp);
}


int
fbs_open(struct fbserv_obj *fbsp, int port)
{
    int i;
    int available_port = port;

    /* Already listening; nothing more to do. */
    if ((*fbsp->fbs_is_listening)(fbsp)) {
	return BRLCAD_OK;
    }

    if (available_port < 0) {
	available_port = 5559;
    } else if (available_port < 1024) {
	available_port += 5559;
    }

    /* Try a reasonable number of times to hang a listen */
    int have_listen = 0;
    for (i = 0; i < MAX_PORT_TRIES; ++i) {
	/*
	 * Hang an unending listen for PKG connections
	 */
	if (!(*fbsp->fbs_listen_on_port)(fbsp, available_port)) {
	    ++available_port;
	} else {
	    have_listen = 1;
	    break;
	}
    }

    if (!have_listen) {
	if (fbsp->msgs)
	    bu_vls_printf(fbsp->msgs, "fbs_open: failed to hang a listen on ports %d - %d\n", port, available_port);
	fbsp->fbs_listener.fbsl_port = -1;
	return BRLCAD_ERROR;
    }

    fbsp->fbs_listener.fbsl_port = available_port;

    (*fbsp->fbs_open_server_handler)(fbsp);

    return BRLCAD_OK;
}


int
fbs_close(struct fbserv_obj *fbsp)
{
    int i;

    /* first drop all clients */
    for (i = 0; i < MAX_CLIENTS; ++i)
	drop_client(fbsp, i);

    (*fbsp->fbs_close_server_handler)(fbsp);

    /* Close the TCP listener if one was created by fbs_listen_on_port(). */
    if (fbsp->fbs_listener.fbsl_listener) {
	pkg_listener_close(fbsp->fbs_listener.fbsl_listener);
	fbsp->fbs_listener.fbsl_listener = NULL;
    } else if (0 <= fbsp->fbs_listener.fbsl_fd) {
	close(fbsp->fbs_listener.fbsl_fd);
    }
    fbsp->fbs_listener.fbsl_fd = -1;
    fbsp->fbs_listener.fbsl_port = -1;

    /* Close the IPC child-end channel if one was created by fbs_open_ipc(). */
    if (fbsp->fbs_listener.fbsl_ipc_child) {
	pkg_close(fbsp->fbs_listener.fbsl_ipc_child);
	fbsp->fbs_listener.fbsl_ipc_child = NULL;
    }

    return BRLCAD_OK;
}

struct pkg_switch *
fbs_pkg_switch(void)
{
    static struct pkg_switch pswitch[FBSERV_PKG_SWITCH_COUNT];
    static const struct fbserv_pkg_handlers handlers = {
	fbs_rfbunknown,
	fbs_rfbauth,
	fbs_rfbopen,
	fbs_rfbclose,
	fbs_rfbfree,
	fbs_rfbclear,
	fbs_rfbread,
	fbs_rfbwrite,
	fbs_rfbreadrect,
	fbs_rfbwriterect,
	fbs_rfbbwreadrect,
	fbs_rfbbwwriterect,
	fbs_rfbcursor,
	fbs_rfbgetcursor,
	fbs_rfbsetcursor,
	fbs_rfbscursor,
	fbs_rfbwindow,
	fbs_rfbzoom,
	fbs_rfbview,
	fbs_rfbgetview,
	fbs_rfbrmap,
	fbs_rfbwmap,
	fbs_rfbflush,
	fbs_rfbpoll,
	fbs_rfbhelp
    };
    static int initialized = 0;

    if (!initialized) {
	if (fbserv_pkg_switch_init(pswitch, FBSERV_PKG_SWITCH_COUNT,
		&handlers) != 0)
	    return NULL;
	initialized = 1;
    }

    return pswitch;
}

void
fbs_setup_socket(int fd)
{
    int on     = 1;
    int retval = 0;

#if defined(SO_KEEPALIVE)
    /* FIXME: better to show an error message but need thread considerations for strerror */
    if ((retval = setsockopt(fd, SOL_SOCKET, SO_KEEPALIVE, (char *)&on, sizeof(on))) < 0) {
	bu_log("setsockopt (SO_KEEPALIVE) error return: %d", retval);
    }
#endif
#if defined(SO_RCVBUF)
    /* try to set our buffers up larger */
    {
	int m = -1;
	int n = -1;
	int val;
	int size;

	for (size = 256; size > 16; size /= 2) {
	    val = size * 1024;
	    m = setsockopt(fd, SOL_SOCKET, SO_RCVBUF,
			   (char *)&val, sizeof(val));
	    val = size * 1024;
	    n = setsockopt(fd, SOL_SOCKET, SO_SNDBUF,
			   (char *)&val, sizeof(val));
	    if (m >= 0 && n >= 0) break;
	}

	if (m < 0 || n < 0)
	    bu_log("setup_socket: setsockopt()");
    }
#endif
}

/*
 * Process arrivals from existing clients.
 */
void
fbs_existing_client_handler(void *clientData, int UNUSED(mask))
{
    register int i;
    struct fbserv_client *fbscp = (struct fbserv_client *)clientData;
    struct fbserv_obj *fbsp = fbscp->fbsc_fbsp;
    int fd = fbscp->fbsc_fd;

    for (i = MAX_CLIENTS - 1; i >= 0; i--) {
	if (fbsp->fbs_clients[i].fbsc_fd == 0)
	    continue;

	fbsp->fbs_clients[i].fbsc_pkg->pkc_server_data = (void *)&fbsp->fbs_clients[i];

	if ((pkg_process(fbsp->fbs_clients[i].fbsc_pkg)) < 0)
	    bu_log("pkg_process error encountered (1)\n");

	/* Act on deferred drop requested by a handler (e.g. token mismatch).
	 * We must not call pkg_close() from inside pkg_process's loop. */
	if (fbsp->fbs_clients[i].fbsc_pending_drop) {
	    drop_client(fbsp, i);
	    continue;
	}

	if (fbsp->fbs_clients[i].fbsc_fd != fd)
	    continue;

	if (pkg_suckin(fbsp->fbs_clients[i].fbsc_pkg) <= 0) {
	    /* Probably EOF */
	    drop_client(fbsp, i);

	    continue;
	}

	if ((pkg_process(fbsp->fbs_clients[i].fbsc_pkg)) < 0)
	    bu_log("pkg_process error encountered (2)\n");

	/* Deferred drop from second-pass handler */
	if (fbsp->fbs_clients[i].fbsc_pending_drop) {
	    drop_client(fbsp, i);
	    continue;
	}
    }

    if (fbsp->fbs_callback != (void (*)(void *))FBS_CALLBACK_NULL) {
	/* need to cast func pointer explicitly to get the function call */
	void (*cfp)(void *);
	cfp = (void (*)(void *))fbsp->fbs_callback;
	cfp(fbsp->fbs_clientData);
    }
}


int
fbs_new_client(struct fbserv_obj *fbsp, struct pkg_conn *pcp, void *data)
{
    if (pcp == PKC_ERROR)
	return -1;

    for (int i = MAX_CLIENTS - 1; i >= 0; i--) {
	/* this slot is being used */
	if (fbsp->fbs_clients[i].fbsc_fd != 0)
	    continue;

	/* Found an available slot */
	fbsp->fbs_clients[i].fbsc_fd = pkg_get_read_fd(pcp);
	fbsp->fbs_clients[i].fbsc_pkg = pcp;
	fbsp->fbs_clients[i].fbsc_fbsp = fbsp;
	fbsp->fbs_clients[i].fbsc_auth_ok = 0;
	fbsp->fbs_clients[i].fbsc_pending_drop = 0;
	fbs_setup_socket(pkg_get_read_fd(pcp));
	/* Phase A1 (ert reliability): mark the framebuffer as having an
	 * active streaming client so its dimensions won't be mutated under
	 * the client's feet by GUI events (e.g. window resize). */
	if (fbsp->fbs_fbp && fbsp->fbs_fbp->i)
	    fbsp->fbs_fbp->i->if_active_clients++;

	/* Point pkc_server_data at the fbserv_client so handlers can
	 * reach back to the fbserv_obj (needed for auth checks). */
	pcp->pkc_server_data = (void *)&fbsp->fbs_clients[i];

#ifdef HAVE_OPENSSL_SSL_H
	/* Optional TLS: if the server has a TLS context, perform the
	 * server-side handshake before the first PKG message is read. */
	if (fbsp->fbs_tls_ctx) {
	    if (fbserv_tls_accept((SSL_CTX *)fbsp->fbs_tls_ctx, pcp) != FBSERV_TLS_OK) {
		bu_log("fbs_new_client: TLS handshake failed — dropping client\n");
		pkg_close(pcp);
		fbsp->fbs_clients[i].fbsc_pkg = PKC_NULL;
		fbsp->fbs_clients[i].fbsc_fd = 0;
		return -1;
	    }
	}
#endif

	(*fbsp->fbs_open_client_handler)(fbsp, i, data);

	return i;
    }

    bu_log("fbs_new_client: too many clients\n");
    pkg_close(pcp);

    return -1;
}


static void
_fbs_ipc_comm_error(const char *msg)
{
    bu_log("%s", msg);
}


int
fbs_open_ipc(struct fbserv_obj *fbsp)
{
    struct pkg_conn *pe = NULL, *ce = NULL;
    struct pkg_conn *pc;
    int i;

    if (pkg_pair(&pe, &ce, fbs_pkg_switch(), _fbs_ipc_comm_error) != 0) {
	if (fbsp->msgs)
	    bu_vls_printf(fbsp->msgs, "fbs_open_ipc: pkg_pair failed\n");
	return BRLCAD_ERROR;
    }

    /* Move child-end fd(s) above bu_process_create()'s close(3..19) sweep
     * so the fd survives into the spawned subprocess after exec().          */
    if (pkg_move_high_fd(ce, 64) != 0)
	bu_log("fbs_open_ipc: pkg_move_high_fd failed; fd may be swept\n");

    bu_log("fbs_open_ipc: pe rfd=%d wfd=%d  ce rfd=%d wfd=%d  ce_addr='%s'\n",
	   pkg_get_read_fd(pe), pkg_get_write_fd(pe),
	   pkg_get_read_fd(ce), pkg_get_write_fd(ce),
	   pkg_child_addr_env(ce) ? pkg_child_addr_env(ce) : "(null)");

    pc = pe;

    /* Find an empty client slot and register the pre-connected pkg_conn.
     * The transport may be a unidirectional pipe pair internally, so the
     * fd to monitor for readability must come from pkg_get_read_fd(),
     * not from struct internals.  fbsc_fd is used by callers
     * (Tcl_CreateFileHandler, select-based loops, etc.) as the fd to
     * monitor for readability, so it must always hold a valid (>= 0) fd. */
    int effective_fd = pkg_get_read_fd(pc);

    for (i = MAX_CLIENTS - 1; i >= 0; i--) {
	if (fbsp->fbs_clients[i].fbsc_fd != 0)
	    continue;

	fbsp->fbs_clients[i].fbsc_fd      = effective_fd;
	fbsp->fbs_clients[i].fbsc_pkg     = pc;
	fbsp->fbs_clients[i].fbsc_fbsp    = fbsp;
	fbsp->fbs_clients[i].fbsc_auth_ok = 1; /* IPC client is implicitly trusted */
	fbsp->fbs_clients[i].fbsc_pending_drop = 0;
	fbsp->fbs_clients[i].fbsc_is_ipc  = 1;
	pc->pkc_server_data = (void *)&fbsp->fbs_clients[i];
	/* Phase A1 (ert reliability): mark the fb dimension lock active. */
	if (fbsp->fbs_fbp && fbsp->fbs_fbp->i)
	    fbsp->fbs_fbp->i->if_active_clients++;

	/* Call the IPC-specific open handler if one is registered, otherwise
	 * fall back to the generic TCP client handler.  Callers that use the
	 * generic handler must ensure it handles NULL data gracefully.        */
	if (fbsp->fbs_open_ipc_client_handler)
	    (*fbsp->fbs_open_ipc_client_handler)(fbsp, i, NULL);
	else if (fbsp->fbs_open_client_handler)
	    (*fbsp->fbs_open_client_handler)(fbsp, i, NULL);

	break;
    }

    if (i < 0) {
	bu_log("fbs_open_ipc: too many clients\n");
	pkg_close(ce);
	pkg_close(pc);
	return BRLCAD_ERROR;
    }

    /* Store the child end so fbs_ipc_child_addr_env() can retrieve it, and
     * so fbs_close() can close it when the session ends.                    */
    if (fbsp->fbs_listener.fbsl_ipc_child) {
	bu_log("fbs_open_ipc: closing old ipc_child rfd=%d wfd=%d\n",
	       pkg_get_read_fd(fbsp->fbs_listener.fbsl_ipc_child),
	       pkg_get_write_fd(fbsp->fbs_listener.fbsl_ipc_child));
	pkg_close(fbsp->fbs_listener.fbsl_ipc_child);
    }
    fbsp->fbs_listener.fbsl_ipc_child = ce;
    bu_log("fbs_open_ipc: new ipc_child stored (rfd=%d wfd=%d) slot=%d effective_fd=%d\n",
	   pkg_get_read_fd(ce), pkg_get_write_fd(ce), i, effective_fd);

#ifndef _WIN32
    /* Phase C1 (ert reliability): make the parent's read+write fds
     * non-blocking so that pkg_suckin() driven from a Qt event-loop
     * notifier can never block the GUI thread when the kernel buffer
     * is empty (e.g. on a level-triggered notifier re-fire after the
     * fd has already been drained).  pkg_suckin() returns 0 with
     * pc->pkc_would_block == 1 in that case, which the IPC handler
     * uses to short-circuit cleanly.  The child end (ce) intentionally
     * remains blocking — rt's writes against it must back-pressure on
     * a slow consumer rather than spinning. */
    {
	int rfd = pkg_get_read_fd(pc);
	int wfd = pkg_get_write_fd(pc);
	int flags;
	if (rfd >= 0 && (flags = fcntl(rfd, F_GETFL, 0)) != -1)
	    (void)fcntl(rfd, F_SETFL, flags | O_NONBLOCK);
	if (wfd >= 0 && wfd != rfd && (flags = fcntl(wfd, F_GETFL, 0)) != -1)
	    (void)fcntl(wfd, F_SETFL, flags | O_NONBLOCK);
    }
#endif

    return BRLCAD_OK;
}


/**
 * Public wrapper around drop_client() so toolkit-specific client
 * handlers (e.g. qged's QFBIPCSocket) can request a clean teardown
 * after detecting EOF / error on the IPC fd without depending on the
 * select-based fbs_existing_client_handler path.
 */
void
fbs_drop_client(struct fbserv_obj *fbsp, int sub)
{
    if (!fbsp || sub < 0 || sub >= MAX_CLIENTS)
	return;
    drop_client(fbsp, sub);
}


int
fbs_drop_ipc_clients(struct fbserv_obj *fbsp)
{
    int dropped = 0;

    if (!fbsp)
	return 0;

    for (int i = 0; i < MAX_CLIENTS; ++i) {
	if (fbsp->fbs_clients[i].fbsc_fd == 0 ||
	    !fbsp->fbs_clients[i].fbsc_is_ipc)
	    continue;
	bu_log("fbs_drop_ipc_clients: dropping IPC client slot %d fd=%d\n",
	       i, fbsp->fbs_clients[i].fbsc_fd);
	drop_client(fbsp, i);
	++dropped;
    }

    return dropped;
}


const char *
fbs_ipc_child_addr_env(struct fbserv_obj *fbsp)
{
    if (!fbsp || !fbsp->fbs_listener.fbsl_ipc_child)
	return NULL;
    return pkg_child_addr_env(fbsp->fbs_listener.fbsl_ipc_child);
}


/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
