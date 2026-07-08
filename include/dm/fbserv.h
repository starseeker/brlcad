/*                        F B S E R V . H
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
/** @file fbserv.h
 *
 * @brief
 * This header exposes the transitional libdm-hosted fbserv packet transport
 * entry points used for TCP based (and, via libpkg, local-IPC based)
 * communication between a framebuffer and a remote process.  The shared
 * protocol constants, auth helpers, backend operation table, and server object
 * layout live in imgstream/fbserv.h.
 *
 * Asynchronous interprocess communication and event monitoring is (as of 2021)
 * still very much platform and toolkit specific.  Hence, these data structures
 * contain some void pointers which are used by individual applications to
 * connect their own specific methods (for example, Tcl_Channel) to handle this
 * problem.  Improving this to be more generic and less dependent on specific
 * toolkits and/or platform mechanisms would be a laudable goal, if practical.
 *
 * pkg IPC path (fbs_open_ipc):
 *   Instead of binding a TCP listen socket, creates a pkg_pair() and
 *   immediately wraps the parent end as a pre-connected pkg_conn client.
 *   The child-end address is retrieved via fbs_ipc_child_addr_env() and
 *   passed to the spawned rt subprocess via the PKG_ADDR environment
 *   variable (set with bu_setenv() before the fork, cleared after).
 *
 */

#ifndef DM_FBSERV_H
#define DM_FBSERV_H

#include "common.h"
#include <stddef.h>
#include "imgstream/fbserv.h"
#include "pkg.h"
#include "dm/defines.h"

__BEGIN_DECLS

#define MAX_PORT_TRIES 100

/**
 * Initialise a newly allocated fbserv object to a closed, empty state.
 *
 * This clears all fields and sets listener fd/port sentinels.  It is intended
 * for freshly allocated storage, not live servers with active clients.
 */
DM_EXPORT extern int fbs_init(struct fbserv_obj *fbsp);
DM_EXPORT extern int fbs_open(struct fbserv_obj *fbsp, int port);
DM_EXPORT extern int fbs_close(struct fbserv_obj *fbsp);
DM_EXPORT extern int fbs_set_backend(struct fbserv_obj *fbsp,
	const struct fbserv_fb_ops *ops,
	void *ctx);
DM_EXPORT extern void fbs_clear_backend(struct fbserv_obj *fbsp);

/**
 * Install or clear toolkit/application transport callbacks.
 *
 * These helpers are the preferred access path for migrated callers that only
 * need to configure fbserv transport behavior.  Toolkit hosts that own their
 * event-loop integration should use the narrow listener/client helpers below
 * until the packet layer moves out of libdm.
 */
DM_EXPORT extern void fbs_set_transport(struct fbserv_obj *fbsp,
	const struct fbserv_transport_ops *ops);
DM_EXPORT extern void fbs_clear_transport(struct fbserv_obj *fbsp);

/**
 * Narrow accessors for migrated callers.  They avoid depending on the
 * transitional fbserv_obj listener/client/callback layout.
 */
DM_EXPORT extern int fbs_can_open_ipc(const struct fbserv_obj *fbsp);
DM_EXPORT extern int fbs_can_open_network(const struct fbserv_obj *fbsp);
DM_EXPORT extern int fbs_can_close(const struct fbserv_obj *fbsp);
DM_EXPORT extern int fbs_listener_port(const struct fbserv_obj *fbsp);
DM_EXPORT extern int fbs_listener_fd(const struct fbserv_obj *fbsp);
DM_EXPORT extern void fbs_set_listener_fd(struct fbserv_obj *fbsp, int fd);
DM_EXPORT extern void *fbs_listener_channel(const struct fbserv_obj *fbsp);
DM_EXPORT extern void fbs_set_listener_channel(struct fbserv_obj *fbsp,
	void *chan);
DM_EXPORT extern void *fbs_listener_handler_data(struct fbserv_obj *fbsp);
DM_EXPORT extern struct fbserv_obj *fbs_listener_owner(void *listener_data);
DM_EXPORT extern int fbs_listener_data_fd(const void *listener_data);
DM_EXPORT extern struct pkg_listener *fbs_listener_data_pkg_listener(
	void *listener_data);
DM_EXPORT extern void fbs_set_listener_pkg_listener(struct fbserv_obj *fbsp,
	struct pkg_listener *listener);
DM_EXPORT extern int fbs_client_active(const struct fbserv_obj *fbsp,
	int sub);
DM_EXPORT extern int fbs_client_fd(const struct fbserv_obj *fbsp, int sub);
DM_EXPORT extern struct pkg_conn *fbs_client_pkg(const struct fbserv_obj *fbsp,
	int sub);
DM_EXPORT extern struct pkg_conn *fbs_client_data_pkg(void *client_data);
DM_EXPORT extern int fbs_client_data_fd(const void *client_data);
DM_EXPORT extern void *fbs_client_channel(const struct fbserv_obj *fbsp,
	int sub);
DM_EXPORT extern void fbs_set_client_channel(struct fbserv_obj *fbsp,
	int sub,
	void *chan);
DM_EXPORT extern void *fbs_client_handler(const struct fbserv_obj *fbsp,
	int sub);
DM_EXPORT extern void fbs_set_client_handler(struct fbserv_obj *fbsp,
	int sub,
	void *handler);
DM_EXPORT extern void *fbs_client_handler_data(struct fbserv_obj *fbsp,
	int sub);
DM_EXPORT extern void fbs_set_client_data_channel(void *client_data,
	void *chan);
DM_EXPORT extern void fbs_set_legacy_framebuffer(struct fbserv_obj *fbsp,
	struct fb *fbp);
DM_EXPORT extern struct fb *fbs_legacy_framebuffer(struct fbserv_obj *fbsp);
DM_EXPORT extern int fbs_framebuffer_backend_installed(
	const struct fbserv_obj *fbsp);
DM_EXPORT extern int fbs_framebuffer_info(struct fbserv_obj *fbsp,
	struct fbserv_fb_info *info);
DM_EXPORT extern int fbs_framebuffer_writerect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	const unsigned char *rgb);
DM_EXPORT extern int fbs_framebuffer_view(struct fbserv_obj *fbsp,
	int xcenter,
	int ycenter,
	int xzoom,
	int yzoom);
DM_EXPORT extern int fbs_framebuffer_cursor(struct fbserv_obj *fbsp,
	int mode,
	int x,
	int y);
DM_EXPORT extern int fbs_framebuffer_flush(struct fbserv_obj *fbsp);
DM_EXPORT extern int fbs_framebuffer_poll(struct fbserv_obj *fbsp);
DM_EXPORT extern struct pkg_switch *fbs_pkg_switch(void);
DM_EXPORT extern void fbs_setup_socket(int fd);
DM_EXPORT extern int fbs_new_client(struct fbserv_obj *fbsp, struct pkg_conn *pcp, void *data);
DM_EXPORT extern void fbs_existing_client_handler(void *clientData, int mask);

/**
 * @brief Tear down an fbserv client slot.
 *
 * Closes the libpkg connection, calls the registered TCP or IPC close
 * handler, and resets the slot fields.  Safe to call from a
 * Qt/Tcl/etc. handler that has detected EOF or an error on its
 * fd-monitor channel and does not (or cannot) rely on the
 * select-based fbs_existing_client_handler() path to do the drop.
 *
 * @p sub is the client index returned by fbs_new_client() / stored in
 * fbsp->fbs_clients[*].  Out-of-range or unused indices are no-ops.
 */
DM_EXPORT extern void fbs_drop_client(struct fbserv_obj *fbsp, int sub);
DM_EXPORT extern int fbs_drop_ipc_clients(struct fbserv_obj *fbsp);

/**
 * @brief Open an IPC-based framebuffer server (no TCP listen socket).
 *
 * Creates a pkg_pair(), wraps the parent end as a pre-connected pkg_conn
 * client (bypassing the TCP accept loop entirely), and registers it via
 * fbs_open_ipc_client_handler (or fbs_open_client_handler if the former is
 * NULL).  The child end's address is stored in fbsp->fbs_listener.fbsl_ipc_child
 * and can be retrieved with fbs_ipc_child_addr_env().
 *
 * Callers should:
 *   1. Call fbs_open_ipc() to start the server.
 *   2. Call fbs_ipc_child_addr_env() to get "PKG_ADDR=<addr>".
 *   3. Set that variable in the parent env (bu_setenv) before spawning rt.
 *   4. Clear the variable after bu_process_create() returns.
 *   5. Pass "-F 0" (or any port spec) to rt so it opens a remote framebuffer;
 *      if_remote.c will detect PKG_ADDR and use the IPC channel instead.
 *
 * The child end is closed (pkg_close) when fbs_close() is called.
 *
 * @return BRLCAD_OK on success, BRLCAD_ERROR if pkg_pair fails.
 */
DM_EXPORT extern int fbs_open_ipc(struct fbserv_obj *fbsp);

/**
 * @brief Return the "PKG_ADDR=<addr>" env string for the spawned child.
 *
 * Valid after a successful fbs_open_ipc() call and until fbs_close() is
 * called.  Returns NULL if no IPC channel is active.
 *
 * The string is owned by the internal channel struct and must not be freed
 * by the caller.
 */
DM_EXPORT extern const char *fbs_ipc_child_addr_env(struct fbserv_obj *fbsp);

/**
 * Initialise @p fbsp->fbs_auth_token for session authentication.
 *
 * If the FBSERV_TOKEN environment variable is already set to a valid
 * 64-hex-char token, that value is used directly so that the hosting
 * application can pre-supply a known token and pass the same value to
 * child processes (e.g. set FBSERV_TOKEN before execing rt/pix-fb).
 * Token authentication works regardless of whether TLS is enabled.
 *
 * If FBSERV_TOKEN is not set or is the wrong length, a fresh random
 * 256-bit token is generated.  Falls back to /dev/urandom or a
 * time+PID PRNG when OpenSSL is not available.
 *
 * Call this before fbs_open() so the token is ready for the first
 * connecting client.  Returns a pointer to fbsp->fbs_auth_token.
 */
DM_EXPORT extern const char *fbs_generate_token(struct fbserv_obj *fbsp);


__END_DECLS

#endif /* DM_FBSERV_H */
/** @} */
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
