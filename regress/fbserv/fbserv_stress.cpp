/*              F B S E R V _ S T R E S S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2025 United States Government as represented by
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
/** @file regress/fbserv/fbserv_stress.cpp
 *
 * Stress test for the fbserv token-isolation mechanism.
 *
 * Simulates the rtwizard / fbserv connection scenario: multiple concurrent
 * workers each launch an isolated fbserv instance, connect through
 * libimgstream's authenticated remote-framebuffer client, and verify pixel
 * round trips.  Both TCP and local IPC transports use the same public client
 * contract exercised by applications.
 *
 * The test is designed to flush out races such as:
 *   - Client connecting before fbserv is ready to accept (ECONNREFUSED)
 *   - Token cross-talk when two instances share a port
 *   - Port collisions between concurrent workers
 *
 * Each worker retries the TCP connect with exponential back-off for up
 * to MAX_CONNECT_WAIT_SEC seconds, matching the behaviour expected of
 * a robust rtwizard front-end.
 *
 * Usage:
 *   fbserv_stress [--workers N] [--base-port P] [--ipc]
 *
 * Exits 0 on full success, 1 if any worker fails.
 */

#include "common.h"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

/* system headers */
#include <fcntl.h>
#include <time.h>

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#include "bu/app.h"
#include "bu/log.h"
#include "bu/process.h"
#include "bu/snooze.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "imgstream/fb_compat.h"
#include "pkg.h"

/* --------------------------------------------------------------------------
 * Tunables
 * -------------------------------------------------------------------------- */

/** Maximum seconds a worker will keep retrying the TCP connect.
 *  Overridable via --timeout on the command line. */
static int g_connect_timeout_sec = 15;

/** Initial retry interval in milliseconds. */
static const int RETRY_INTERVAL_MS = 50;

/** Maximum retry interval cap in milliseconds. */
static const int RETRY_MAX_MS = 500;

/** Length (in hex chars) of the authentication token. */
static const int TOKEN_LEN = 64;

/** Framebuffer dimensions used by every isolated worker. */
static const size_t FRAMEBUFFER_SIZE = 512;

/* --------------------------------------------------------------------------
 * Global state
 * -------------------------------------------------------------------------- */

static std::atomic<int> g_port_counter;
static std::mutex g_print_mutex;
/** Serialize the environment update and process creation so each child
 * inherits its worker's token.  The lock is released as soon as the child
 * is created; server startup and all subsequent worker activity overlap. */
static std::mutex g_token_env_mutex;

/* --------------------------------------------------------------------------
 * Logging helpers
 * -------------------------------------------------------------------------- */

static void
worker_log(int wid, const char *fmt, ...)
{
    std::lock_guard<std::mutex> lk(g_print_mutex);
    va_list ap;
    char buf[512];
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    fprintf(stderr, "[worker %d] %s\n", wid, buf);
}

/* --------------------------------------------------------------------------
 * Token generation
 * -------------------------------------------------------------------------- */

/** Fill buf (TOKEN_LEN+1 bytes) with a freshly-generated lower-case hex token.
 *  Uses /dev/urandom where available; falls back to seeded rand(). */
static void
generate_token(char *buf, int wid)
{
    unsigned char raw[32];
    bool ok = false;

#if !defined(_WIN32)
    {
        int fd = open("/dev/urandom", O_RDONLY);
        if (fd >= 0) {
            ssize_t n = read(fd, raw, sizeof(raw));
            if (n == (ssize_t)sizeof(raw))
                ok = true;
            close(fd);
        }
    }
#endif

    if (!ok) {
        /* Last-resort: seed from time + PID + worker id */
        unsigned int seed = (unsigned int)(time(NULL)) ^ ((unsigned int)bu_pid() << 8) ^ (unsigned int)wid;
        srand(seed);
        for (size_t i = 0; i < sizeof(raw); i++)
            raw[i] = (unsigned char)(rand() & 0xff);
    }

    for (size_t i = 0; i < sizeof(raw); i++)
        snprintf(buf + i * 2, 3, "%02x", (unsigned int)raw[i]);
    buf[TOKEN_LEN] = '\0';
}

/* --------------------------------------------------------------------------
 * Port probe: check if a TCP port is already in use
 * -------------------------------------------------------------------------- */

/** Return true if something is already listening on localhost:port. */
static bool
port_in_use(int port)
{
    char ps[16];
    snprintf(ps, sizeof(ps), "%d", port);
    struct pkg_conn *pc = pkg_open("localhost", ps, 0, 0, 0, NULL, NULL);
    if (pc != PKC_ERROR) {
        pkg_close(pc);
        return true;
    }
    return false;
}

/* --------------------------------------------------------------------------
 * Worker
 * -------------------------------------------------------------------------- */

struct WorkerResult {
    int worker_id = 0;
    bool passed = false;
    std::string message;
};

static void
stop_worker(struct bu_process **fbserv_proc)
{
    if (!fbserv_proc || !*fbserv_proc)
        return;

    int fbserv_pid = bu_process_pid(*fbserv_proc);
    if (fbserv_pid > 0 && bu_process_alive(*fbserv_proc))
        (void)bu_pid_terminate(fbserv_pid);
    (void)bu_process_wait_n(fbserv_proc, 0);
}

static void
run_worker(int wid, int base_port, bool use_ipc, const char *fbserv_path,
	WorkerResult &result)
{
    result.worker_id = wid;

    char endpoint[MAXPATHLEN] = {0};
    char remote_spec[MAXPATHLEN + 8] = {0};
    char port_str[16] = {0};
    if (use_ipc) {
	if (pkg_ipc_addr(endpoint, sizeof(endpoint), "fbserv-stress") != 0) {
	    result.message = "failed to allocate an IPC endpoint";
	    return;
	}
	snprintf(remote_spec, sizeof(remote_spec), "ipc:%s", endpoint);
	worker_log(wid, "using IPC endpoint %s", endpoint);
    } else {
	int port = base_port + g_port_counter.fetch_add(1);
	for (int tries = 0; tries < 200 && port_in_use(port); tries++)
	    port = base_port + g_port_counter.fetch_add(1);
	snprintf(port_str, sizeof(port_str), "%d", port);
	snprintf(remote_spec, sizeof(remote_spec), "tcp:localhost:%d", port);
	worker_log(wid, "using port %d", port);
    }

    /* ------------------------------------------------------------------
     * 2. Generate a unique session token for this worker.
     * ------------------------------------------------------------------ */
    char token[TOKEN_LEN + 1];
    generate_token(token, wid);
    worker_log(wid, "token %.16s...", token);

    /* ------------------------------------------------------------------
     * 3. Build fbserv command line and set the token in the environment.
     *    We pass -A (strict auth) so that a connection without the token
     *    is rejected — this verifies isolation between workers.
     *
     *    fbserv -A -w 512 -n 512 -F /dev/mem (-p port | -I address)
     *
     *    NOTE: /dev/mem is the in-memory framebuffer (null device);
     *    the -F form launches a single-frame-buffer server which stays
     *    alive until killed, which is what rtwizard uses.
     * ------------------------------------------------------------------ */
    const char *argv_fbserv[12] = {
	fbserv_path, "-A", "-w", "512", "-n", "512", "-F", "/dev/mem",
	use_ipc ? "-I" : "-p", use_ipc ? endpoint : port_str, NULL, NULL
    };

    /* Pre-supply the token via the environment variable.  Environment state
     * is process-wide, so protect the set/create/clear sequence while still
     * allowing the created servers to start simultaneously. */
    struct bu_process *fbserv_proc = NULL;
    {
        std::lock_guard<std::mutex> lk(g_token_env_mutex);
#if defined(HAVE_UNISTD_H) && !defined(_WIN32)
        setenv("FBSERV_TOKEN", token, 1 /* overwrite */);
#else
        /* Windows fallback: _putenv_s makes its own copy. */
        _putenv_s("FBSERV_TOKEN", token);
#endif
        bu_process_create(&fbserv_proc, argv_fbserv, BU_PROCESS_DEFAULT);

        /* Immediately clear so subsequent processes do not inherit this token. */
#if defined(HAVE_UNISTD_H) && !defined(_WIN32)
        unsetenv("FBSERV_TOKEN");
#else
        _putenv_s("FBSERV_TOKEN", "");
#endif
    }

    int fbserv_pid = bu_process_pid(fbserv_proc);
    if (fbserv_pid == -1) {
        result.message = "failed to launch fbserv";
        worker_log(wid, "ERROR: %s", result.message.c_str());
        stop_worker(&fbserv_proc);
        return;
    }
    worker_log(wid, "fbserv launched (pid %d)", fbserv_pid);

    /* fbserv needs a moment to bind its endpoint.  Retry the same public
     * libimgstream open that applications use, with bounded backoff. */
    struct imgstream_fb_remote_options remote_options =
	IMGSTREAM_FB_REMOTE_OPTIONS_INIT;
    remote_options.auth_token = token;
    imgstream_fb_t *remote = NULL;
    int interval_ms = RETRY_INTERVAL_MS;
    int64_t deadline = bu_gettime() + BU_SEC2USEC(g_connect_timeout_sec);
    int attempt = 0;

    while (!remote && bu_gettime() < deadline) {
	attempt++;
	remote = imgstream_fb_open_remote(remote_spec, FRAMEBUFFER_SIZE,
	    FRAMEBUFFER_SIZE, &remote_options);
	if (!remote) {
	    worker_log(wid, "connect attempt %d failed, retrying in %dms...",
		attempt, interval_ms);
            bu_snooze((int64_t)interval_ms * (int64_t)1000);
            interval_ms = (interval_ms * 2 < RETRY_MAX_MS) ? interval_ms * 2 : RETRY_MAX_MS;
        }
    }

    if (!remote) {
	result.message = "timed out connecting to fbserv";
	worker_log(wid, "ERROR: %s (tried %d times)", result.message.c_str(), attempt);
	stop_worker(&fbserv_proc);
	return;
    }
    worker_log(wid, "connected after %d attempt(s)", attempt);

    const unsigned char pixels[12] = {
	(unsigned char)wid, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
    };
    unsigned char readback[sizeof(pixels)] = {0};
    if (imgstream_fb_writerect(remote, 0, 0, 2, 2, pixels) != 4 ||
	imgstream_fb_readrect(remote, 0, 0, 2, 2, readback) != 4 ||
	memcmp(pixels, readback, sizeof(pixels)) != 0) {
	result.message = "authenticated framebuffer pixel round trip failed";
	worker_log(wid, "ERROR: %s", result.message.c_str());
	imgstream_fb_close(remote);
	stop_worker(&fbserv_proc);
	return;
    }
    worker_log(wid, "authenticated pixel round trip completed");

    imgstream_fb_close(remote);
    stop_worker(&fbserv_proc);

    result.passed = true;
    result.message = "OK";
    worker_log(wid, "PASSED");
}

/* --------------------------------------------------------------------------
 * main
 * -------------------------------------------------------------------------- */

static void
print_usage(const char *prog)
{
	fprintf(stderr,
	    "Usage: %s [options]\n"
	    "  --workers N      number of parallel workers (default: 8)\n"
	    "  --base-port P    first TCP port to try (PID-derived by default)\n"
	    "  --ipc            use isolated local IPC endpoints instead of TCP\n"
	    "  --timeout T      seconds to wait for fbserv to become ready (default: 15)\n"
        "  --help           show this help\n",
        prog);
}

int
main(int argc, const char *argv[])
{
    bu_setprogname(argv[0]);

    int num_workers = 8;
    /* Default base port: derive a per-invocation offset from the PID so that
     * two concurrent or rapidly sequential test runs don't collide on the same
     * port range.  Stays in the 10000-19999 range (well above reserved ports
     * and away from the typical ephemeral-port range). */
    int pid_offset = (bu_pid() % 1000) * 10;
    int base_port  = 10000 + pid_offset;
    bool use_ipc = false;

    for (int i = 1; i < argc; i++) {
        if (BU_STR_EQUAL(argv[i], "--workers") && i + 1 < argc) {
            num_workers = atoi(argv[++i]);
        } else if (BU_STR_EQUAL(argv[i], "--base-port") && i + 1 < argc) {
            base_port = atoi(argv[++i]);
	} else if (BU_STR_EQUAL(argv[i], "--timeout") && i + 1 < argc) {
	    g_connect_timeout_sec = atoi(argv[++i]);
	} else if (BU_STR_EQUAL(argv[i], "--ipc")) {
	    use_ipc = true;
        } else if (BU_STR_EQUAL(argv[i], "--help") || BU_STR_EQUAL(argv[i], "-h")) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    if (num_workers < 1 || num_workers > 256) {
        fprintf(stderr, "num_workers must be between 1 and 256\n");
        return 1;
    }

    g_port_counter.store(0);

    fprintf(stderr, "fbserv_stress: launching %d parallel workers (%s)\n",
	    num_workers, use_ipc ? "local IPC" : "loopback TCP");

    std::vector<WorkerResult> results(num_workers);
    std::vector<std::thread> threads;
    threads.reserve(num_workers);

    /* bu_dir consults process-global application path state.  Resolve the
     * invariant executable location before starting concurrent workers. */
    char fbserv_path[MAXPATHLEN];
    bu_dir(fbserv_path, MAXPATHLEN, BU_DIR_BIN, "fbserv", BU_DIR_EXT, NULL);

    for (int w = 0; w < num_workers; w++)
	threads.emplace_back(run_worker, w, base_port, use_ipc, fbserv_path,
	    std::ref(results[w]));

    for (auto &t : threads)
        t.join();

    int passed = 0, failed = 0;
    for (const auto &r : results) {
        if (r.passed)
            passed++;
        else {
            failed++;
            fprintf(stderr, "FAILED worker %d: %s\n", r.worker_id, r.message.c_str());
        }
    }

    fprintf(stderr, "fbserv_stress: %d/%d workers passed\n", passed, num_workers);

    return (failed == 0) ? 0 : 1;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
