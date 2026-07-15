/*                    F B S E R V _ S T A T E . C
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
/** @file libimgstream/fbserv_state.c */

#include "common.h"

#include <stdlib.h>
#include <string.h>

#include "bu/defines.h"
#include "bu/str.h"
#include "imgstream/fbserv.h"


static int
fbserv_has_backend(const struct fbserv_obj *fbsp)
{
    return fbsp && fbsp->fbs_fb_ops && fbsp->fbs_fb_ctx;
}


const char *
fbserv_obj_generate_token(struct fbserv_obj *fbsp)
{
    const char *env_token;

    if (!fbsp)
	return NULL;

    env_token = getenv(FBSERV_AUTH_TOKEN_ENVVAR);
    if (env_token && strlen(env_token) == FBSERV_AUTH_TOKEN_LEN) {
	bu_strlcpy(fbsp->fbs_auth_token, env_token,
	    sizeof(fbsp->fbs_auth_token));
    } else {
	fbserv_generate_token(fbsp->fbs_auth_token);
    }

    return fbsp->fbs_auth_token;
}


int
fbserv_obj_init(struct fbserv_obj *fbsp)
{
    if (!fbsp)
	return BRLCAD_ERROR;

    memset(fbsp, 0, sizeof(*fbsp));
    fbsp->fbs_listener.fbsl_fd = -1;
    fbsp->fbs_listener.fbsl_port = -1;
    fbsp->fbs_listener.fbsl_fbsp = fbsp;
    return BRLCAD_OK;
}


int
fbserv_set_backend(struct fbserv_obj *fbsp,
	const struct fbserv_fb_ops *ops,
	void *ctx)
{
    if (!fbsp || !ops || !ctx || !ops->info)
	return BRLCAD_ERROR;

    fbsp->fbs_fb_ops = ops;
    fbsp->fbs_fb_ctx = ctx;
    return BRLCAD_OK;
}


void
fbserv_clear_backend(struct fbserv_obj *fbsp)
{
    if (!fbsp)
	return;

    fbsp->fbs_fb_ops = NULL;
    fbsp->fbs_fb_ctx = NULL;
}


void
fbserv_set_transport(struct fbserv_obj *fbsp, const struct fbserv_transport_ops *ops)
{
    if (!fbsp)
	return;

    if (!ops) {
	fbserv_clear_transport(fbsp);
	return;
    }

    fbsp->fbs_is_listening = ops->is_listening;
    fbsp->fbs_listen_on_port = ops->listen_on_port;
    fbsp->fbs_open_server_handler = ops->open_server_handler;
    fbsp->fbs_close_server_handler = ops->close_server_handler;
    fbsp->fbs_open_client_handler = ops->open_client_handler;
    fbsp->fbs_close_client_handler = ops->close_client_handler;
    fbsp->fbs_open_ipc_client_handler = ops->open_ipc_client_handler;
    fbsp->fbs_close_ipc_client_handler = ops->close_ipc_client_handler;
}


void
fbserv_clear_transport(struct fbserv_obj *fbsp)
{
    if (!fbsp)
	return;

    fbsp->fbs_is_listening = NULL;
    fbsp->fbs_listen_on_port = NULL;
    fbsp->fbs_open_server_handler = NULL;
    fbsp->fbs_close_server_handler = NULL;
    fbsp->fbs_open_client_handler = NULL;
    fbsp->fbs_close_client_handler = NULL;
    fbsp->fbs_open_ipc_client_handler = NULL;
    fbsp->fbs_close_ipc_client_handler = NULL;
}


int
fbserv_can_open_ipc(const struct fbserv_obj *fbsp)
{
    return fbsp && fbsp->fbs_open_ipc_client_handler &&
	(fbsp->fbs_close_ipc_client_handler || fbsp->fbs_close_client_handler);
}


int
fbserv_can_open_network(const struct fbserv_obj *fbsp)
{
    return fbsp && fbsp->fbs_is_listening &&
	fbsp->fbs_listen_on_port && fbsp->fbs_open_server_handler;
}


int
fbserv_can_close(const struct fbserv_obj *fbsp)
{
    return fbsp && fbsp->fbs_close_server_handler;
}


int
fbserv_listener_port(const struct fbserv_obj *fbsp)
{
    return fbsp ? fbsp->fbs_listener.fbsl_port : -1;
}


int
fbserv_listener_fd(const struct fbserv_obj *fbsp)
{
    return fbsp ? fbsp->fbs_listener.fbsl_fd : -1;
}


void
fbserv_set_listener_fd(struct fbserv_obj *fbsp, int fd)
{
    if (!fbsp)
	return;
    fbsp->fbs_listener.fbsl_fd = fd;
}


void *
fbserv_listener_channel(const struct fbserv_obj *fbsp)
{
    return fbsp ? fbsp->fbs_listener.fbsl_chan : NULL;
}


void
fbserv_set_listener_channel(struct fbserv_obj *fbsp, void *chan)
{
    if (!fbsp)
	return;
    fbsp->fbs_listener.fbsl_chan = chan;
}


void *
fbserv_listener_handler_data(struct fbserv_obj *fbsp)
{
    return fbsp ? (void *)&fbsp->fbs_listener : NULL;
}


struct fbserv_obj *
fbserv_listener_owner(void *listener_data)
{
    struct fbserv_listener *fbslp = (struct fbserv_listener *)listener_data;
    return fbslp ? fbslp->fbsl_fbsp : NULL;
}


int
fbserv_listener_data_fd(const void *listener_data)
{
    const struct fbserv_listener *fbslp =
	(const struct fbserv_listener *)listener_data;
    return fbslp ? fbslp->fbsl_fd : -1;
}


struct pkg_listener *
fbserv_listener_data_pkg_listener(void *listener_data)
{
    struct fbserv_listener *fbslp = (struct fbserv_listener *)listener_data;
    return fbslp ? fbslp->fbsl_listener : NULL;
}


void
fbserv_set_listener_pkg_listener(struct fbserv_obj *fbsp,
	struct pkg_listener *listener)
{
    if (!fbsp)
	return;
    fbsp->fbs_listener.fbsl_listener = listener;
}


int
fbserv_client_active(const struct fbserv_obj *fbsp, int sub)
{
    return fbsp && sub >= 0 && sub < FBSERV_MAX_CLIENTS &&
	fbsp->fbs_clients[sub].fbsc_fd != 0;
}


int
fbserv_client_fd(const struct fbserv_obj *fbsp, int sub)
{
    if (!fbsp || sub < 0 || sub >= FBSERV_MAX_CLIENTS)
	return -1;
    return fbsp->fbs_clients[sub].fbsc_fd;
}


struct pkg_conn *
fbserv_client_pkg(const struct fbserv_obj *fbsp, int sub)
{
    if (!fbsp || sub < 0 || sub >= FBSERV_MAX_CLIENTS)
	return NULL;
    return fbsp->fbs_clients[sub].fbsc_pkg;
}


struct pkg_conn *
fbserv_client_data_pkg(void *client_data)
{
    struct fbserv_client *fbscp = (struct fbserv_client *)client_data;
    return fbscp ? fbscp->fbsc_pkg : NULL;
}


int
fbserv_client_data_fd(const void *client_data)
{
    const struct fbserv_client *fbscp =
	(const struct fbserv_client *)client_data;
    return fbscp ? fbscp->fbsc_fd : -1;
}


void *
fbserv_client_channel(const struct fbserv_obj *fbsp, int sub)
{
    if (!fbsp || sub < 0 || sub >= FBSERV_MAX_CLIENTS)
	return NULL;
    return fbsp->fbs_clients[sub].fbsc_chan;
}


void
fbserv_set_client_channel(struct fbserv_obj *fbsp, int sub, void *chan)
{
    if (!fbsp || sub < 0 || sub >= FBSERV_MAX_CLIENTS)
	return;
    fbsp->fbs_clients[sub].fbsc_chan = chan;
}


void *
fbserv_client_handler(const struct fbserv_obj *fbsp, int sub)
{
    if (!fbsp || sub < 0 || sub >= FBSERV_MAX_CLIENTS)
	return NULL;
    return fbsp->fbs_clients[sub].fbsc_handler;
}


void
fbserv_set_client_handler(struct fbserv_obj *fbsp, int sub, void *handler)
{
    if (!fbsp || sub < 0 || sub >= FBSERV_MAX_CLIENTS)
	return;
    fbsp->fbs_clients[sub].fbsc_handler = handler;
}


void *
fbserv_client_handler_data(struct fbserv_obj *fbsp, int sub)
{
    if (!fbsp || sub < 0 || sub >= FBSERV_MAX_CLIENTS)
	return NULL;
    if (fbsp->fbs_clients[sub].fbsc_fd == 0)
	return NULL;
    return (void *)&fbsp->fbs_clients[sub];
}


void
fbserv_set_client_data_channel(void *client_data, void *chan)
{
    struct fbserv_client *fbscp = (struct fbserv_client *)client_data;
    if (!fbscp)
	return;
    fbscp->fbsc_chan = chan;
}


int
fbserv_framebuffer_backend_installed(const struct fbserv_obj *fbsp)
{
    return fbserv_has_backend(fbsp);
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
