/*                        F B S E R V . C
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file fbserv/fbserv.c
 *
 * Standalone transport host for libimgstream's framebuffer server.
 */

#include "common.h"

#include <errno.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_SYS_SELECT_H
#  include <sys/select.h>
#endif
#ifdef HAVE_SYS_SOCKET_H
#  include <sys/socket.h>
#endif

#include "bio.h"
#include "bsocket.h"
#include "bu/app.h"
#include "bu/exit.h"
#include "bu/getopt.h"
#include "bu/log.h"
#include "brlobol/display_session.h"
#if defined(HAVE_QTCAD_OBOL_DISPLAY_PROVIDER)
#  include "qtcad/display_provider.h"
#elif defined(HAVE_TCLCAD_OBOL_DISPLAY_PROVIDER)
#  include "tclcad/setup.h"
#endif
#include "imgstream/fb_compat.h"
#include "imgstream/fbserv.h"
#include "pkg.h"

#if defined(HAVE_QTCAD_OBOL_DISPLAY_PROVIDER) || \
    defined(HAVE_TCLCAD_OBOL_DISPLAY_PROVIDER)
#  define HAVE_BRLCAD_OBOL_DISPLAY_PROVIDER 1
#endif

#ifdef HAVE_BRLCAD_OBOL_DISPLAY_PROVIDER
static int
fbserv_obol_display_provider_register(void)
{
#ifdef HAVE_QTCAD_OBOL_DISPLAY_PROVIDER
    return qtcad_obol_display_provider_register();
#elif defined(HAVE_TCLCAD_OBOL_DISPLAY_PROVIDER)
    return tclcad_obol_display_provider_register();
#else
    return 0;
#endif
}
#endif

static struct fbserv_obj server;
static imgstream_fb_t *framebuffer = NULL;
#ifdef HAVE_BRLCAD_OBOL_DISPLAY_PROVIDER
static brlobol_display_session_t *display_session = NULL;
#endif
static const char *framebuffer_name = NULL;
static const char *ipc_address = NULL;
static int framebuffer_width = 0;
static int framebuffer_height = 0;
static int port = -1;
static int require_auth = 0;
static int use_tls = 0;
static int verbose = 0;

static const char usage[] =
    "Usage: fbserv [-v] [-T] [-A] [-F framebuffer] [-s size]\n"
    "              [-w width] [-n height] (-p port | -I ipc_address)\n"
    "       fbserv port [framebuffer]\n"
    "\n"
    "The standalone server owns one image-stream framebuffer for its lifetime.\n"
    "If -F is omitted, a memory framebuffer is used.\n"
    "  -T  enable TLS (requires an OpenSSL build)\n"
    "  -A  require FBSERV_TOKEN authentication\n"
    "  -I  listen on a libpkg IPC address\n";

static void
communications_error(const char *msg)
{
    if (msg)
	fputs(msg, stderr);
}

static int
transport_is_listening(struct fbserv_obj *fbsp)
{
    return fbsp && fbs_listener_fd(fbsp) >= 0;
}

static int
transport_listen_on_port(struct fbserv_obj *fbsp, int listen_port)
{
    char service[32] = {0};
    struct pkg_listener *listener;

    if (!fbsp)
	return 0;
    snprintf(service, sizeof(service), "%d", listen_port);
    listener = pkg_listen(service, NULL, 16, communications_error);
    if (!listener)
	return 0;
    fbs_set_listener_pkg_listener(fbsp, listener);
    fbs_set_listener_fd(fbsp, pkg_get_listener_fd(listener));
    return 1;
}

static void transport_open_server(struct fbserv_obj *UNUSED(fbsp)) { }
static void transport_close_server(struct fbserv_obj *UNUSED(fbsp)) { }
static void transport_open_client(struct fbserv_obj *UNUSED(fbsp),
	int UNUSED(sub), void *UNUSED(data)) { }
static void transport_close_client(struct fbserv_obj *UNUSED(fbsp),
	int UNUSED(sub)) { }

static const struct fbserv_transport_ops transport_ops = {
    transport_is_listening,
    transport_listen_on_port,
    transport_open_server,
    transport_close_server,
    transport_open_client,
    transport_close_client,
    transport_open_client,
    transport_close_client
};

static int
parse_args(int argc, char **argv)
{
    int c;
    int port_set = 0;

    while ((c = bu_getopt(argc, argv, "vTAI:F:s:w:n:S:W:N:p:h?")) != -1) {
	switch (c) {
	    case 'v': verbose = 1; break;
	    case 'T': use_tls = 1; break;
	    case 'A': require_auth = 1; break;
	    case 'I': ipc_address = bu_optarg; break;
	    case 'F': framebuffer_name = bu_optarg; break;
	    case 's':
	    case 'S': framebuffer_width = framebuffer_height = atoi(bu_optarg); break;
	    case 'w':
	    case 'W': framebuffer_width = atoi(bu_optarg); break;
	    case 'n':
	    case 'N': framebuffer_height = atoi(bu_optarg); break;
	    case 'p': port = atoi(bu_optarg); port_set = 1; break;
	    default: return 0;
	}
    }

    if (bu_optind < argc && !port_set && !ipc_address) {
	port = atoi(argv[bu_optind++]);
	port_set = 1;
    }
    if (bu_optind < argc && !framebuffer_name)
	framebuffer_name = argv[bu_optind++];

    return bu_optind == argc && (port_set || ipc_address);
}

static int
active_clients(const struct fbserv_obj *fbsp)
{
    int count = 0;
    for (int i = 0; i < FBSERV_MAX_CLIENTS; i++)
	if (fbs_client_active(fbsp, i))
	    count++;
    return count;
}

static int
service_client(int sub)
{
    struct fbserv_process_result result;
    int ret = fbs_drain_client(&server, sub, 1024 * 1024, 1, 10000, &result);
    if (ret != BRLCAD_OK && verbose)
	bu_log("fbserv: client %d protocol or transport error\n", sub);
    return result.disconnected;
}

static int
service_loop(struct pkg_listener *listener)
{
    for (;;) {
	fd_set read_fds;
	struct timeval timeout = {60, 0};
	int max_fd = -1;
	int listener_fd = listener ? pkg_get_listener_fd(listener) : -1;

	FD_ZERO(&read_fds);
	if (listener_fd >= 0) {
	    FD_SET(listener_fd, &read_fds);
	    max_fd = listener_fd;
	}

	for (int i = 0; i < FBSERV_MAX_CLIENTS; i++) {
	    int fd = fbs_client_fd(&server, i);
	    if (fbs_client_active(&server, i) && fd >= 0) {
		FD_SET(fd, &read_fds);
		if (fd > max_fd)
		    max_fd = fd;
	    }
	}

	if (framebuffer) {
	    if (imgstream_fb_poll_rate(framebuffer) > 0) {
		long usec = imgstream_fb_poll_rate(framebuffer);
		timeout.tv_sec = usec / 1000000;
		timeout.tv_usec = usec % 1000000;
	    }
	}

#ifdef _WIN32
	/* Named-pipe listeners do not expose a select-compatible descriptor. */
	if (listener && listener_fd < 0) {
	    struct pkg_conn *pcp = pkg_accept(listener, fbs_pkg_switch(),
		communications_error, 0);
	    int sub = fbs_new_client(&server, pcp, NULL);
	    while (sub >= 0 && fbs_client_active(&server, sub))
		(void)service_client(sub);
	    continue;
	}
#endif

	if (max_fd < 0)
	    return BRLCAD_ERROR;

	int selected = select(max_fd + 1, &read_fds, NULL, NULL, &timeout);
	if (selected < 0) {
	    if (errno == EINTR)
		continue;
	    return BRLCAD_ERROR;
	}
	if (selected == 0) {
	    if (framebuffer && imgstream_fb_poll(framebuffer) != 0)
		return BRLCAD_OK;
	    continue;
	}

	if (listener_fd >= 0 && FD_ISSET(listener_fd, &read_fds)) {
	    struct pkg_conn *pcp = pkg_accept(listener, fbs_pkg_switch(),
		communications_error, 0);
	    if (pcp != PKC_ERROR && pcp != PKC_NULL)
		(void)fbs_new_client(&server, pcp, NULL);
	}

	for (int i = 0; i < FBSERV_MAX_CLIENTS; i++) {
	    int fd = fbs_client_fd(&server, i);
	    if (fbs_client_active(&server, i) && fd >= 0 &&
		FD_ISSET(fd, &read_fds))
		(void)service_client(i);
	}

	if (!listener && active_clients(&server) == 0)
	    return BRLCAD_OK;
    }
}

int
main(int argc, char **argv)
{
    struct pkg_listener *listener = NULL;
    int ret = BRLCAD_ERROR;

    bu_setprogname(argv[0]);
#ifdef SIGPIPE
    (void)signal(SIGPIPE, SIG_IGN);
#endif

    if (!parse_args(argc, argv)) {
	fputs(usage, stderr);
	return 1;
    }

    if (fbs_init(&server) != BRLCAD_OK)
	return 1;
    fbs_set_transport(&server, &transport_ops);
    server.fbs_require_auth = require_auth;
    fprintf(stderr, "fbserv: Session token: %s\n", fbs_generate_token(&server));

    if (use_tls || getenv("FBSERV_TLS")) {
	server.fbs_tls_ctx = fbs_tls_server_context_create(getenv("FBSERV_TLS_CERT"),
		getenv("FBSERV_TLS_KEY"));
	if (!server.fbs_tls_ctx)
	    bu_exit(1, "fbserv: TLS is unavailable or failed to initialize\n");
    }

    const size_t fb_width = framebuffer_width > 0 ?
	(size_t)framebuffer_width : 0;
    const size_t fb_height = framebuffer_height > 0 ?
	(size_t)framebuffer_height : 0;
    if (imgstream_fb_spec_kind(framebuffer_name) == IMGSTREAM_FB_SPEC_DISPLAY) {
#ifdef HAVE_BRLCAD_OBOL_DISPLAY_PROVIDER
	if (!fbserv_obol_display_provider_register())
	    bu_exit(1, "fbserv: selected UI toolkit has no usable Obol display provider\n");
	display_session = brlobol_display_session_open(framebuffer_name,
	    fb_width, fb_height, "BRL-CAD fbserv");
	framebuffer = brlobol_display_session_framebuffer(display_session);
#else
	bu_exit(1, "fbserv: display targets require a Qt or Tcl/Tk Obol provider; use an imgstream target\n");
#endif
    } else {
	framebuffer = imgstream_fb_open(framebuffer_name, fb_width, fb_height);
    }
    if (!framebuffer)
	bu_exit(1, "fbserv: unsupported or invalid image-stream framebuffer '%s'\n",
	    framebuffer_name ? framebuffer_name : "/dev/mem");
    if (imgstream_fbserv_set_framebuffer(&server, framebuffer) != BRLCAD_OK)
	bu_exit(1, "fbserv: unable to install framebuffer backend\n");

    if (ipc_address) {
	if (pkg_addr_is_ipc_listener(ipc_address)) {
	    listener = pkg_listen(ipc_address, NULL, 16, communications_error);
	    if (!listener)
		bu_exit(1, "fbserv: unable to listen on %s\n", ipc_address);
	    fbs_set_listener_pkg_listener(&server, listener);
	    fbs_set_listener_fd(&server, pkg_get_listener_fd(listener));
	} else {
	    struct pkg_conn *pcp = pkg_connect_addr(ipc_address,
		fbs_pkg_switch(), communications_error);
	    if (pcp == PKC_ERROR || pcp == PKC_NULL)
		bu_exit(1, "fbserv: unable to connect to %s\n", ipc_address);
	    (void)fbs_new_client(&server, pcp, NULL);
	}
    } else {
	if (fbs_open(&server, port) != BRLCAD_OK)
	    bu_exit(1, "fbserv: unable to open a network listener\n");
	listener = fbs_listener_data_pkg_listener(fbs_listener_handler_data(&server));
	fprintf(stderr, "fbserv: Listening on port %d\n", fbs_listener_port(&server));
    }

    ret = service_loop(listener);
    (void)fbs_close(&server);
    (void)imgstream_fbserv_set_framebuffer(&server, NULL);
#ifdef HAVE_BRLCAD_OBOL_DISPLAY_PROVIDER
    if (display_session) {
	brlobol_display_session_close(display_session);
	display_session = NULL;
	framebuffer = NULL;
    } else
#endif
    imgstream_fb_close(framebuffer);

    if (server.fbs_tls_ctx)
	fbs_tls_server_context_destroy(server.fbs_tls_ctx);
    return ret == BRLCAD_OK ? 0 : 1;
}
