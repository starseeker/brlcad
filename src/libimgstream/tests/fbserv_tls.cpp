/*                    F B S E R V _ T L S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libimgstream/tests/fbserv_tls.cpp */

#include "common.h"

#include "bu/app.h"

#include <atomic>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <thread>

#include "imgstream/fb_compat.h"
#include "imgstream/fbserv.h"
#include "pkg.h"


static void
communications_error(const char *message)
{
    if (message)
	std::fprintf(stderr, "%s", message);
}


static int not_listening(struct fbserv_obj *UNUSED(fbsp)) { return 0; }
static int no_listen(struct fbserv_obj *UNUSED(fbsp), int UNUSED(port)) { return 0; }
static void no_server_handler(struct fbserv_obj *UNUSED(fbsp)) { }
static void no_client_handler(struct fbserv_obj *UNUSED(fbsp),
	int UNUSED(sub), void *UNUSED(data)) { }
static void no_close_client_handler(struct fbserv_obj *UNUSED(fbsp),
	int UNUSED(sub)) { }


int
main(int UNUSED(argc), char **argv)
{
    bu_setprogname(argv[0]);
#ifdef SIGPIPE
    (void)std::signal(SIGPIPE, SIG_IGN);
#endif

    static const struct fbserv_transport_ops transport_ops = {
	not_listening,
	no_listen,
	no_server_handler,
	no_server_handler,
	no_client_handler,
	no_close_client_handler,
	no_client_handler,
	no_close_client_handler
    };
    static const char token[] =
	"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
    struct fbserv_obj server;
    imgstream_fb_t *server_fb = NULL;
    pkg_listener_t *listener = NULL;
    void *tls_context = NULL;
    char fingerprint[FBSERV_AUTH_TOKEN_LEN + 1] = {0};
    char wrong_fingerprint[FBSERV_AUTH_TOKEN_LEN + 1] = {0};
    char remote_spec[80] = {0};
    std::atomic<int> server_failures(0);
    int failures = 0;

    if (std::strlen(token) != FBSERV_AUTH_TOKEN_LEN ||
	fbs_init(&server) != BRLCAD_OK) {
	std::fprintf(stderr, "unable to initialize TLS test server\n");
	return 1;
    }
    fbs_set_transport(&server, &transport_ops);
    server.fbs_require_auth = 1;
    std::memcpy(server.fbs_auth_token, token, sizeof(token));

    server_fb = imgstream_fb_open("/dev/mem", 2, 2);
    if (!server_fb ||
	imgstream_fbserv_set_framebuffer(&server, server_fb) != BRLCAD_OK) {
	std::fprintf(stderr, "unable to install TLS test framebuffer\n");
	if (server_fb)
	    imgstream_fb_close(server_fb);
	return 1;
    }

    tls_context = fbs_tls_server_context_create(NULL, NULL);
    if (!tls_context ||
	fbs_tls_server_sha256(tls_context, fingerprint) != BRLCAD_OK) {
	std::fprintf(stderr, "unable to create TLS test certificate\n");
	(void)imgstream_fbserv_set_framebuffer(&server, NULL);
	imgstream_fb_close(server_fb);
	if (tls_context)
	    fbs_tls_server_context_destroy(tls_context);
	return 1;
    }
    server.fbs_tls_ctx = tls_context;

    listener = pkg_listen("0", "127.0.0.1", 2, communications_error);
    int port = listener ? pkg_get_listener_port(listener) : -1;
    if (!listener || port <= 0) {
	std::fprintf(stderr, "unable to create TLS test listener\n");
	(void)imgstream_fbserv_set_framebuffer(&server, NULL);
	imgstream_fb_close(server_fb);
	fbs_tls_server_context_destroy(tls_context);
	return 1;
    }
    std::snprintf(remote_spec, sizeof(remote_spec), "127.0.0.1:%d", port);

    std::thread server_thread([&]() {
	for (int attempt = 0; attempt < 2; attempt++) {
	    struct pkg_conn *connection = pkg_accept(listener, fbs_pkg_switch(),
		communications_error, 0);
	    if (connection == PKC_NULL || connection == PKC_ERROR) {
		server_failures++;
		continue;
	    }
	    int slot = fbs_new_client(&server, connection, NULL);
	    if (slot < 0) {
		server_failures++;
		continue;
	    }
	    while (fbs_client_active(&server, slot)) {
		struct fbserv_process_result result;
		if (fbs_drain_client(&server, slot, 1024 * 1024, 1, 100000,
			&result) != BRLCAD_OK && !result.disconnected) {
		    server_failures++;
		    fbs_drop_client(&server, slot);
		}
	    }
	}
    });

    std::memcpy(wrong_fingerprint, fingerprint, sizeof(fingerprint));
    wrong_fingerprint[0] = wrong_fingerprint[0] == '0' ? '1' : '0';
    struct imgstream_fb_remote_options options =
	IMGSTREAM_FB_REMOTE_OPTIONS_INIT;
    options.auth_token = token;
    options.use_tls = 1;
    options.tls_server_sha256 = wrong_fingerprint;

    imgstream_fb_t *remote = imgstream_fb_open_remote(remote_spec, 2, 2,
	&options);
    if (remote) {
	std::fprintf(stderr, "TLS client accepted the wrong certificate pin\n");
	failures++;
	imgstream_fb_close(remote);
    }

    options.tls_server_sha256 = fingerprint;
    remote = imgstream_fb_open_remote(remote_spec, 2, 2, &options);
    if (!remote) {
	std::fprintf(stderr, "TLS client rejected the correct certificate pin\n");
	failures++;
    } else {
	const unsigned char pixels[12] = {
	    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
	};
	unsigned char readback[12] = {0};
	if (imgstream_fb_writerect(remote, 0, 0, 2, 2, pixels) != 4 ||
	    imgstream_fb_readrect(remote, 0, 0, 2, 2, readback) != 4 ||
	    std::memcmp(pixels, readback, sizeof(pixels)) != 0) {
	    std::fprintf(stderr, "authenticated TLS framebuffer exchange failed\n");
	    failures++;
	}
	imgstream_fb_close(remote);
    }

    server_thread.join();
    if (server_failures.load() != 0) {
	std::fprintf(stderr, "TLS server reported %d lifecycle failure(s)\n",
	    server_failures.load());
	failures++;
    }

    pkg_listener_close(listener);
    (void)fbs_close(&server);
    (void)imgstream_fbserv_set_framebuffer(&server, NULL);
    imgstream_fb_close(server_fb);
    server.fbs_tls_ctx = NULL;
    fbs_tls_server_context_destroy(tls_context);
    return failures ? 1 : 0;
}

