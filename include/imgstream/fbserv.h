/*                        F B S E R V . H
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
/** @addtogroup libimgstream */
/** @{ */
/** @file imgstream/fbserv.h
 *
 * Framebuffer-server protocol contracts for image-stream producers.
 *
 * These types describe the framebuffer operation table and transitional
 * server object consumed by the fbserv packet transport.  They intentionally
 * live with libimgstream so Obol/libged/qtcad hosts can publish image-stream
 * backed framebuffer behavior without depending on a scene-rendering
 * subsystem.
 */
#ifndef IMGSTREAM_FBSERV_H
#define IMGSTREAM_FBSERV_H

#include "common.h"

#include <stddef.h>
#include <stdint.h>
#if defined(HAVE_SYS_TYPES_H)
#  include <sys/types.h>
#endif

__BEGIN_DECLS

#ifndef IMGSTREAM_EXPORT
#  if defined(IMGSTREAM_DLL_EXPORTS) && defined(IMGSTREAM_DLL_IMPORTS)
#    error "Only IMGSTREAM_DLL_EXPORTS or IMGSTREAM_DLL_IMPORTS can be defined, not both."
#  elif defined(IMGSTREAM_DLL_EXPORTS)
#    define IMGSTREAM_EXPORT COMPILER_DLLEXPORT
#  elif defined(IMGSTREAM_DLL_IMPORTS)
#    define IMGSTREAM_EXPORT COMPILER_DLLIMPORT
#  else
#    define IMGSTREAM_EXPORT
#  endif
#endif

/*
 * fbserv wire-protocol message IDs.
 *
 * The historical MSG_* spellings remain as compatibility aliases, but the
 * authoritative constants live here with the framebuffer/image-stream
 * protocol contract rather than in the display-manager umbrella header.
 */
#define FBSERV_MSG_FBOPEN        1
#define FBSERV_MSG_FBCLOSE       2
#define FBSERV_MSG_FBCLEAR       3
#define FBSERV_MSG_FBREAD        4
#define FBSERV_MSG_FBWRITE       5
#define FBSERV_MSG_FBCURSOR      6  /**< @brief fb_cursor(void) */
#define FBSERV_MSG_FBWINDOW      7  /**< @brief OLD */
#define FBSERV_MSG_FBZOOM        8  /**< @brief OLD */
#define FBSERV_MSG_FBSCURSOR     9  /**< @brief OLD */
#define FBSERV_MSG_FBVIEW        10 /**< @brief NEW */
#define FBSERV_MSG_FBGETVIEW     11 /**< @brief NEW */
#define FBSERV_MSG_FBRMAP        12
#define FBSERV_MSG_FBWMAP        13
#define FBSERV_MSG_FBHELP        14
#define FBSERV_MSG_FBREADRECT    15
#define FBSERV_MSG_FBWRITERECT   16
#define FBSERV_MSG_FBFLUSH       17
#define FBSERV_MSG_FBFREE        18
#define FBSERV_MSG_FBGETCURSOR   19 /**< @brief NEW */
#define FBSERV_MSG_DATA          20
#define FBSERV_MSG_RETURN        21
#define FBSERV_MSG_CLOSE         22
#define FBSERV_MSG_ERROR         23
#define FBSERV_MSG_FBPOLL        30 /**< @brief NEW */
#define FBSERV_MSG_FBSETCURSOR   31 /**< @brief NEW in Release 4.4 */
#define FBSERV_MSG_FBBWREADRECT  32 /**< @brief NEW in Release 4.6 */
#define FBSERV_MSG_FBBWWRITERECT 33 /**< @brief NEW in Release 4.6 */
#define FBSERV_MSG_FBAUTH        34 /**< @brief Session token authentication */
#define FBSERV_MSG_NORETURN      100

#ifndef MSG_FBOPEN
#  define MSG_FBOPEN FBSERV_MSG_FBOPEN
#endif
#ifndef MSG_FBCLOSE
#  define MSG_FBCLOSE FBSERV_MSG_FBCLOSE
#endif
#ifndef MSG_FBCLEAR
#  define MSG_FBCLEAR FBSERV_MSG_FBCLEAR
#endif
#ifndef MSG_FBREAD
#  define MSG_FBREAD FBSERV_MSG_FBREAD
#endif
#ifndef MSG_FBWRITE
#  define MSG_FBWRITE FBSERV_MSG_FBWRITE
#endif
#ifndef MSG_FBCURSOR
#  define MSG_FBCURSOR FBSERV_MSG_FBCURSOR
#endif
#ifndef MSG_FBWINDOW
#  define MSG_FBWINDOW FBSERV_MSG_FBWINDOW
#endif
#ifndef MSG_FBZOOM
#  define MSG_FBZOOM FBSERV_MSG_FBZOOM
#endif
#ifndef MSG_FBSCURSOR
#  define MSG_FBSCURSOR FBSERV_MSG_FBSCURSOR
#endif
#ifndef MSG_FBVIEW
#  define MSG_FBVIEW FBSERV_MSG_FBVIEW
#endif
#ifndef MSG_FBGETVIEW
#  define MSG_FBGETVIEW FBSERV_MSG_FBGETVIEW
#endif
#ifndef MSG_FBRMAP
#  define MSG_FBRMAP FBSERV_MSG_FBRMAP
#endif
#ifndef MSG_FBWMAP
#  define MSG_FBWMAP FBSERV_MSG_FBWMAP
#endif
#ifndef MSG_FBHELP
#  define MSG_FBHELP FBSERV_MSG_FBHELP
#endif
#ifndef MSG_FBREADRECT
#  define MSG_FBREADRECT FBSERV_MSG_FBREADRECT
#endif
#ifndef MSG_FBWRITERECT
#  define MSG_FBWRITERECT FBSERV_MSG_FBWRITERECT
#endif
#ifndef MSG_FBFLUSH
#  define MSG_FBFLUSH FBSERV_MSG_FBFLUSH
#endif
#ifndef MSG_FBFREE
#  define MSG_FBFREE FBSERV_MSG_FBFREE
#endif
#ifndef MSG_FBGETCURSOR
#  define MSG_FBGETCURSOR FBSERV_MSG_FBGETCURSOR
#endif
#ifndef MSG_DATA
#  define MSG_DATA FBSERV_MSG_DATA
#endif
#ifndef MSG_RETURN
#  define MSG_RETURN FBSERV_MSG_RETURN
#endif
#ifndef MSG_CLOSE
#  define MSG_CLOSE FBSERV_MSG_CLOSE
#endif
#ifndef MSG_ERROR
#  define MSG_ERROR FBSERV_MSG_ERROR
#endif
#ifndef MSG_FBPOLL
#  define MSG_FBPOLL FBSERV_MSG_FBPOLL
#endif
#ifndef MSG_FBSETCURSOR
#  define MSG_FBSETCURSOR FBSERV_MSG_FBSETCURSOR
#endif
#ifndef MSG_FBBWREADRECT
#  define MSG_FBBWREADRECT FBSERV_MSG_FBBWREADRECT
#endif
#ifndef MSG_FBBWWRITERECT
#  define MSG_FBBWWRITERECT FBSERV_MSG_FBBWWRITERECT
#endif
#ifndef MSG_FBAUTH
#  define MSG_FBAUTH FBSERV_MSG_FBAUTH
#endif
#ifndef MSG_NORETURN
#  define MSG_NORETURN FBSERV_MSG_NORETURN
#endif

/*
 * fbserv session-token authentication.
 *
 * Authentication is protocol-level behavior shared by standalone fbserv,
 * embedded fbserv hosts, and remote framebuffer clients.
 */
#define FBSERV_AUTH_TOKEN_LEN 64
#define FBSERV_AUTH_TOKEN_ENVVAR "FBSERV_TOKEN"

IMGSTREAM_EXPORT void fbserv_generate_token(char *buf);
IMGSTREAM_EXPORT int fbserv_verify_token(const char *provided, const char *expected);

/*
 * Protocol/server state constants.
 *
 * The FBSERV_* names are authoritative.  The historical generic aliases remain
 * for protocol-source compatibility.
 */
#define FBSERV_NET_LONG_LEN 4
#define FBSERV_MAX_CLIENTS 32

#ifndef NET_LONG_LEN
#  define NET_LONG_LEN FBSERV_NET_LONG_LEN
#endif
#ifndef MAX_CLIENTS
#  define MAX_CLIENTS FBSERV_MAX_CLIENTS
#endif
#ifndef FBS_CALLBACK_NULL
#  define FBS_CALLBACK_NULL (void (*)(void))NULL
#endif
#ifndef FBSERV_OBJ_NULL
#  define FBSERV_OBJ_NULL (struct fbserv_obj *)NULL
#endif

#define FBSERV_PKG_SWITCH_COUNT 32
#define FBSERV_FRAMEBUFFER_UNHANDLED (-2)

struct bu_vls;
struct imgstream_fb;
typedef struct imgstream_fb imgstream_fb_t;
struct pkg_conn;
struct pkg_listener;
struct pkg_switch;
struct fbserv_obj;

struct fbserv_fb_info {
    int max_width;
    int max_height;
    int width;
    int height;
};

struct fbserv_colormap {
    uint16_t red[256];
    uint16_t green[256];
    uint16_t blue[256];
};

enum fbserv_backend_op {
    FBSERV_BACKEND_OP_INFO,
    FBSERV_BACKEND_OP_CLEAR,
    FBSERV_BACKEND_OP_READ,
    FBSERV_BACKEND_OP_WRITE,
    FBSERV_BACKEND_OP_READRECT,
    FBSERV_BACKEND_OP_WRITERECT,
    FBSERV_BACKEND_OP_BWREADRECT,
    FBSERV_BACKEND_OP_BWWRITERECT,
    FBSERV_BACKEND_OP_CURSOR,
    FBSERV_BACKEND_OP_GETCURSOR,
    FBSERV_BACKEND_OP_SETCURSOR,
    FBSERV_BACKEND_OP_SCURSOR,
    FBSERV_BACKEND_OP_WINDOW,
    FBSERV_BACKEND_OP_ZOOM,
    FBSERV_BACKEND_OP_VIEW,
    FBSERV_BACKEND_OP_GETVIEW,
    FBSERV_BACKEND_OP_RMAP,
    FBSERV_BACKEND_OP_WMAP,
    FBSERV_BACKEND_OP_FLUSH,
    FBSERV_BACKEND_OP_POLL,
    FBSERV_BACKEND_OP_HELP
};

struct fbserv_fb_ops {
    int (*info)(void *ctx, struct fbserv_fb_info *info);
    int (*clear)(void *ctx, const unsigned char rgb[3]);
    ssize_t (*read)(void *ctx, int x, int y, unsigned char *rgb, size_t count);
    ssize_t (*write)(void *ctx, int x, int y, const unsigned char *rgb, size_t count);
    int (*readrect)(void *ctx, int xmin, int ymin, int width, int height, unsigned char *rgb);
    int (*writerect)(void *ctx, int xmin, int ymin, int width, int height, const unsigned char *rgb);
    int (*bwreadrect)(void *ctx, int xmin, int ymin, int width, int height, unsigned char *bw);
    int (*bwwriterect)(void *ctx, int xmin, int ymin, int width, int height, const unsigned char *bw);
    int (*cursor)(void *ctx, int mode, int x, int y);
    int (*getcursor)(void *ctx, int *mode, int *x, int *y);
    int (*setcursor)(void *ctx, const unsigned char *bits, int xbits, int ybits, int xorig, int yorig);
    int (*scursor)(void *ctx, int mode, int x, int y);
    int (*window)(void *ctx, int xcenter, int ycenter);
    int (*zoom)(void *ctx, int xzoom, int yzoom);
    int (*view)(void *ctx, int xcenter, int ycenter, int xzoom, int yzoom);
    int (*getview)(void *ctx, int *xcenter, int *ycenter, int *xzoom, int *yzoom);
    int (*rmap)(void *ctx, struct fbserv_colormap *cmap);
    int (*wmap)(void *ctx, const struct fbserv_colormap *cmap);
    int (*flush)(void *ctx);
    int (*poll)(void *ctx);
    int (*help)(void *ctx);
};

struct fbserv_transport_ops {
    int (*is_listening)(struct fbserv_obj *);
    int (*listen_on_port)(struct fbserv_obj *, int);
    void (*open_server_handler)(struct fbserv_obj *);
    void (*close_server_handler)(struct fbserv_obj *);
    void (*open_client_handler)(struct fbserv_obj *, int, void *);
    void (*close_client_handler)(struct fbserv_obj *, int);
    void (*open_ipc_client_handler)(struct fbserv_obj *, int, void *);
    void (*close_ipc_client_handler)(struct fbserv_obj *, int);
};

struct fbserv_listener {
    int fbsl_fd;                        /**< @brief socket fd to listen for connections (copy of listener fd) */
    void *fbsl_chan;                    /**< @brief platform/toolkit specific channel */
    int fbsl_port;                      /**< @brief port number to listen on */
    int fbsl_listen;                    /**< @brief !0 means listen for connections */
    struct fbserv_obj *fbsl_fbsp;       /**< @brief points to its fbserv object */
    struct pkg_conn *fbsl_ipc_child;    /**< @brief IPC child-end channel (NULL when using TCP) */
    struct pkg_listener *fbsl_listener; /**< @brief TCP listener (NULL when using IPC or Tcl channel) */
};

struct fbserv_client {
    int fbsc_fd;                        /**< @brief socket to send data down */
    void *fbsc_chan;                    /**< @brief platform/toolkit specific channel */
    void *fbsc_handler;                 /**< @brief platform/toolkit specific handler */
    struct pkg_conn *fbsc_pkg;
    struct fbserv_obj *fbsc_fbsp;       /**< @brief points to its fbserv object */
    int fbsc_auth_ok;                   /**< @brief !0 = client has sent a valid MSG_FBAUTH */
    int fbsc_pending_drop;              /**< @brief !0 = drop this client after pkg_process() returns */
    int fbsc_is_ipc;                    /**< @brief !0 = client is connected via IPC (not TCP) */
};

struct fbserv_obj {
    const struct fbserv_fb_ops *fbs_fb_ops;        /**< @brief framebuffer implementation */
    void *fbs_fb_ctx;                              /**< @brief user data for fbs_fb_ops */
    void *fbs_interp;                              /**< @brief interpreter */
    struct fbserv_listener fbs_listener;           /**< @brief data for listening */
    struct fbserv_client fbs_clients[FBSERV_MAX_CLIENTS]; /**< @brief connected clients */

    int (*fbs_is_listening)(struct fbserv_obj *);          /**< @brief return 1 if listening, else 0 */
    int (*fbs_listen_on_port)(struct fbserv_obj *, int);  /**< @brief return 1 on success, 0 on failure */
    void (*fbs_open_server_handler)(struct fbserv_obj *);   /**< @brief platform/toolkit method to open listener handler */
    void (*fbs_close_server_handler)(struct fbserv_obj *);   /**< @brief platform/toolkit method to close handler listener */
    void (*fbs_open_client_handler)(struct fbserv_obj *, int, void *);   /**< @brief platform/toolkit specific client handler setup */
    void (*fbs_close_client_handler)(struct fbserv_obj *, int);   /**< @brief platform/toolkit method to close handler for client at index client_id */
    void (*fbs_open_ipc_client_handler)(struct fbserv_obj *, int, void *);
    void (*fbs_close_ipc_client_handler)(struct fbserv_obj *, int);

    void (*fbs_callback)(void *);                  /**< @brief callback function */
    void *fbs_clientData;
    struct bu_vls *msgs;
    int fbs_mode;                                  /**< @brief 0-off, 1-underlay, 2-interlay, 3-overlay */

    char fbs_auth_token[FBSERV_AUTH_TOKEN_LEN + 1]; /**< @brief session token; empty = no auth required */
    int fbs_require_auth;     /**< @brief !0 = reject clients that don't send MSG_FBAUTH */
    void *fbs_tls_ctx;        /**< @brief opaque SSL_CTX* for TLS; NULL = no TLS */
};

struct fbserv_pkg_handlers {
    void (*unknown)(struct pkg_conn *, char *);
    void (*auth)(struct pkg_conn *, char *);
    void (*open)(struct pkg_conn *, char *);
    void (*close)(struct pkg_conn *, char *);
    void (*free)(struct pkg_conn *, char *);
    void (*clear)(struct pkg_conn *, char *);
    void (*read)(struct pkg_conn *, char *);
    void (*write)(struct pkg_conn *, char *);
    void (*readrect)(struct pkg_conn *, char *);
    void (*writerect)(struct pkg_conn *, char *);
    void (*bwreadrect)(struct pkg_conn *, char *);
    void (*bwwriterect)(struct pkg_conn *, char *);
    void (*cursor)(struct pkg_conn *, char *);
    void (*getcursor)(struct pkg_conn *, char *);
    void (*setcursor)(struct pkg_conn *, char *);
    void (*scursor)(struct pkg_conn *, char *);
    void (*window)(struct pkg_conn *, char *);
    void (*zoom)(struct pkg_conn *, char *);
    void (*view)(struct pkg_conn *, char *);
    void (*getview)(struct pkg_conn *, char *);
    void (*rmap)(struct pkg_conn *, char *);
    void (*wmap)(struct pkg_conn *, char *);
    void (*flush)(struct pkg_conn *, char *);
    void (*poll)(struct pkg_conn *, char *);
    void (*help)(struct pkg_conn *, char *);
};

struct fbserv_process_result {
    size_t bytes_read;
    int messages_processed;
    int would_block;
    int disconnected;
    int error;
};

typedef void (*fbserv_drop_client_fn)(struct fbserv_obj *fbsp, int sub);

IMGSTREAM_EXPORT int fbserv_pkg_switch_init(
	struct pkg_switch *pswitch,
	size_t switch_count,
	const struct fbserv_pkg_handlers *handlers);
IMGSTREAM_EXPORT struct fbserv_obj *fbserv_conn_obj(struct pkg_conn *pcp);
IMGSTREAM_EXPORT int fbserv_framebuffer_available(
	const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_data_guard(struct pkg_conn *pcp, char *buf);
IMGSTREAM_EXPORT void fbserv_handle_unknown(struct pkg_conn *pcp, char *buf);
IMGSTREAM_EXPORT void fbserv_handle_auth(struct pkg_conn *pcp, char *buf);
IMGSTREAM_EXPORT void fbserv_process_existing_client(
	void *client_data,
	int mask,
	fbserv_drop_client_fn drop_client);
IMGSTREAM_EXPORT const char *fbserv_obj_generate_token(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_obj_init(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_set_backend(struct fbserv_obj *fbsp,
	const struct fbserv_fb_ops *ops,
	void *ctx);
IMGSTREAM_EXPORT void fbserv_clear_backend(struct fbserv_obj *fbsp);
/** Install or detach an imgstream framebuffer as the protocol backend. */
IMGSTREAM_EXPORT int imgstream_fbserv_set_framebuffer(
	struct fbserv_obj *fbsp, imgstream_fb_t *fb);
IMGSTREAM_EXPORT void fbserv_set_transport(struct fbserv_obj *fbsp,
	const struct fbserv_transport_ops *ops);
IMGSTREAM_EXPORT void fbserv_clear_transport(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_can_open_ipc(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_can_open_network(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_can_close(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_listener_port(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_listener_fd(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT void fbserv_set_listener_fd(struct fbserv_obj *fbsp, int fd);
IMGSTREAM_EXPORT void *fbserv_listener_channel(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT void fbserv_set_listener_channel(struct fbserv_obj *fbsp,
	void *chan);
IMGSTREAM_EXPORT void *fbserv_listener_handler_data(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT struct fbserv_obj *fbserv_listener_owner(void *listener_data);
IMGSTREAM_EXPORT int fbserv_listener_data_fd(const void *listener_data);
IMGSTREAM_EXPORT struct pkg_listener *fbserv_listener_data_pkg_listener(
	void *listener_data);
IMGSTREAM_EXPORT void fbserv_set_listener_pkg_listener(struct fbserv_obj *fbsp,
	struct pkg_listener *listener);
IMGSTREAM_EXPORT int fbserv_client_active(const struct fbserv_obj *fbsp,
	int sub);
IMGSTREAM_EXPORT int fbserv_client_fd(const struct fbserv_obj *fbsp, int sub);
IMGSTREAM_EXPORT struct pkg_conn *fbserv_client_pkg(const struct fbserv_obj *fbsp,
	int sub);
IMGSTREAM_EXPORT struct pkg_conn *fbserv_client_data_pkg(void *client_data);
IMGSTREAM_EXPORT int fbserv_client_data_fd(const void *client_data);
IMGSTREAM_EXPORT void *fbserv_client_channel(const struct fbserv_obj *fbsp,
	int sub);
IMGSTREAM_EXPORT void fbserv_set_client_channel(struct fbserv_obj *fbsp,
	int sub,
	void *chan);
IMGSTREAM_EXPORT void *fbserv_client_handler(const struct fbserv_obj *fbsp,
	int sub);
IMGSTREAM_EXPORT void fbserv_set_client_handler(struct fbserv_obj *fbsp,
	int sub,
	void *handler);
IMGSTREAM_EXPORT void *fbserv_client_handler_data(struct fbserv_obj *fbsp,
	int sub);
IMGSTREAM_EXPORT void fbserv_set_client_data_channel(void *client_data,
	void *chan);
IMGSTREAM_EXPORT int fbserv_framebuffer_backend_installed(
	const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_backend_op_installed(
	const struct fbserv_obj *fbsp,
	enum fbserv_backend_op op);
IMGSTREAM_EXPORT int fbserv_backend_info(struct fbserv_obj *fbsp,
	struct fbserv_fb_info *info);
IMGSTREAM_EXPORT int fbserv_backend_clear(struct fbserv_obj *fbsp,
	const unsigned char rgb[3]);
IMGSTREAM_EXPORT ssize_t fbserv_backend_read(struct fbserv_obj *fbsp,
	int x,
	int y,
	unsigned char *rgb,
	size_t count);
IMGSTREAM_EXPORT ssize_t fbserv_backend_write(struct fbserv_obj *fbsp,
	int x,
	int y,
	const unsigned char *rgb,
	size_t count);
IMGSTREAM_EXPORT int fbserv_backend_readrect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	unsigned char *rgb);
IMGSTREAM_EXPORT int fbserv_backend_writerect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	const unsigned char *rgb);
IMGSTREAM_EXPORT int fbserv_backend_bwreadrect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	unsigned char *bw);
IMGSTREAM_EXPORT int fbserv_backend_bwwriterect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	const unsigned char *bw);
IMGSTREAM_EXPORT int fbserv_backend_cursor(struct fbserv_obj *fbsp,
	int mode,
	int x,
	int y);
IMGSTREAM_EXPORT int fbserv_backend_getcursor(struct fbserv_obj *fbsp,
	int *mode,
	int *x,
	int *y);
IMGSTREAM_EXPORT int fbserv_backend_setcursor(struct fbserv_obj *fbsp,
	const unsigned char *bits,
	int xbits,
	int ybits,
	int xorig,
	int yorig);
IMGSTREAM_EXPORT int fbserv_backend_scursor(struct fbserv_obj *fbsp,
	int mode,
	int x,
	int y);
IMGSTREAM_EXPORT int fbserv_backend_window(struct fbserv_obj *fbsp,
	int xcenter,
	int ycenter);
IMGSTREAM_EXPORT int fbserv_backend_zoom(struct fbserv_obj *fbsp,
	int xzoom,
	int yzoom);
IMGSTREAM_EXPORT int fbserv_backend_view(struct fbserv_obj *fbsp,
	int xcenter,
	int ycenter,
	int xzoom,
	int yzoom);
IMGSTREAM_EXPORT int fbserv_backend_getview(struct fbserv_obj *fbsp,
	int *xcenter,
	int *ycenter,
	int *xzoom,
	int *yzoom);
IMGSTREAM_EXPORT int fbserv_backend_rmap(struct fbserv_obj *fbsp,
	struct fbserv_colormap *cmap);
IMGSTREAM_EXPORT int fbserv_backend_wmap(struct fbserv_obj *fbsp,
	const struct fbserv_colormap *cmap);
IMGSTREAM_EXPORT int fbserv_backend_flush(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_backend_poll(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbserv_backend_help(struct fbserv_obj *fbsp);

/* Complete framebuffer server lifecycle and packet transport API. */
IMGSTREAM_EXPORT int fbs_init(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbs_open(struct fbserv_obj *fbsp, int port);
IMGSTREAM_EXPORT int fbs_close(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT void *fbs_tls_server_context_create(const char *certfile,
	const char *keyfile);
IMGSTREAM_EXPORT void fbs_tls_server_context_destroy(void *ctx);
IMGSTREAM_EXPORT int fbs_set_backend(struct fbserv_obj *fbsp,
	const struct fbserv_fb_ops *ops,
	void *ctx);
IMGSTREAM_EXPORT void fbs_clear_backend(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT void fbs_set_transport(struct fbserv_obj *fbsp,
	const struct fbserv_transport_ops *ops);
IMGSTREAM_EXPORT void fbs_clear_transport(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbs_can_open_ipc(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbs_can_open_network(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbs_can_close(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbs_listener_port(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbs_listener_fd(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT void fbs_set_listener_fd(struct fbserv_obj *fbsp, int fd);
IMGSTREAM_EXPORT void *fbs_listener_channel(const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT void fbs_set_listener_channel(struct fbserv_obj *fbsp,
	void *chan);
IMGSTREAM_EXPORT void *fbs_listener_handler_data(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT struct fbserv_obj *fbs_listener_owner(void *listener_data);
IMGSTREAM_EXPORT int fbs_listener_data_fd(const void *listener_data);
IMGSTREAM_EXPORT struct pkg_listener *fbs_listener_data_pkg_listener(
	void *listener_data);
IMGSTREAM_EXPORT void fbs_set_listener_pkg_listener(struct fbserv_obj *fbsp,
	struct pkg_listener *listener);
IMGSTREAM_EXPORT int fbs_client_active(const struct fbserv_obj *fbsp, int sub);
IMGSTREAM_EXPORT int fbs_client_fd(const struct fbserv_obj *fbsp, int sub);
IMGSTREAM_EXPORT struct pkg_conn *fbs_client_pkg(const struct fbserv_obj *fbsp,
	int sub);
IMGSTREAM_EXPORT struct pkg_conn *fbs_client_data_pkg(void *client_data);
IMGSTREAM_EXPORT int fbs_client_data_fd(const void *client_data);
IMGSTREAM_EXPORT void *fbs_client_channel(const struct fbserv_obj *fbsp,
	int sub);
IMGSTREAM_EXPORT void fbs_set_client_channel(struct fbserv_obj *fbsp,
	int sub,
	void *chan);
IMGSTREAM_EXPORT void *fbs_client_handler(const struct fbserv_obj *fbsp,
	int sub);
IMGSTREAM_EXPORT void fbs_set_client_handler(struct fbserv_obj *fbsp,
	int sub,
	void *handler);
IMGSTREAM_EXPORT void *fbs_client_handler_data(struct fbserv_obj *fbsp,
	int sub);
IMGSTREAM_EXPORT void fbs_set_client_data_channel(void *client_data,
	void *chan);
IMGSTREAM_EXPORT int fbs_framebuffer_backend_installed(
	const struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbs_framebuffer_info(struct fbserv_obj *fbsp,
	struct fbserv_fb_info *info);
IMGSTREAM_EXPORT int fbs_framebuffer_writerect(struct fbserv_obj *fbsp,
	int xmin,
	int ymin,
	int width,
	int height,
	const unsigned char *rgb);
IMGSTREAM_EXPORT int fbs_framebuffer_view(struct fbserv_obj *fbsp,
	int xcenter,
	int ycenter,
	int xzoom,
	int yzoom);
IMGSTREAM_EXPORT int fbs_framebuffer_cursor(struct fbserv_obj *fbsp,
	int mode,
	int x,
	int y);
IMGSTREAM_EXPORT int fbs_framebuffer_flush(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbs_framebuffer_poll(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT struct pkg_switch *fbs_pkg_switch(void);
IMGSTREAM_EXPORT void fbs_setup_socket(int fd);
IMGSTREAM_EXPORT int fbs_new_client(struct fbserv_obj *fbsp,
	struct pkg_conn *pcp,
	void *data);
IMGSTREAM_EXPORT int fbs_process_client_bytes(struct fbserv_obj *fbsp,
	int sub,
	const unsigned char *data,
	size_t data_size,
	unsigned char **remaining,
	size_t *remaining_size,
	struct fbserv_process_result *result);
IMGSTREAM_EXPORT int fbs_drain_client(struct fbserv_obj *fbsp,
	int sub,
	size_t byte_budget,
	int iteration_limit,
	int64_t time_budget_usec,
	struct fbserv_process_result *result);
IMGSTREAM_EXPORT void fbs_existing_client_handler(void *client_data, int mask);
IMGSTREAM_EXPORT void fbs_drop_client(struct fbserv_obj *fbsp, int sub);
IMGSTREAM_EXPORT int fbs_drop_ipc_clients(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT int fbs_open_ipc(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT const char *fbs_ipc_child_addr_env(struct fbserv_obj *fbsp);
IMGSTREAM_EXPORT const char *fbs_generate_token(struct fbserv_obj *fbsp);

__END_DECLS

#endif /* IMGSTREAM_FBSERV_H */
/** @} */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
