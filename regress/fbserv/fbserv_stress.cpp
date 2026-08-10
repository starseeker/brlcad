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
 * Simulates the rtwizard / fbserv connection scenario: multiple
 * concurrent workers each launch their own fbserv instance on a
 * unique port with a unique token, connect through libimgstream's remote
 * framebuffer client, and verify authenticated state and pixel operations.
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
 *   fbserv_stress [--workers N] [--base-port P]
 *
 * Exits 0 on full success, 1 if any worker fails.
 */

#include "common.h"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <sstream>
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

#include "bnetwork.h"

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
static bool g_use_ipc = false;

/** Initial retry interval in milliseconds. */
static const int RETRY_INTERVAL_MS = 50;

/** Maximum retry interval cap in milliseconds. */
static const int RETRY_MAX_MS = 500;

/** Length (in hex chars) of the authentication token. */
static const int TOKEN_LEN = 64;

/* --------------------------------------------------------------------------
 * Global state
 * -------------------------------------------------------------------------- */

static std::atomic<int> g_port_counter;
static std::mutex g_print_mutex;
static char g_fbserv_path[MAXPATHLEN] = {0};
/** Serialise the putenv(FBSERV_TOKEN) + fork + putenv-clear sequence so that
 *  concurrent workers do not race on the process-wide environment.  Each
 *  invocation of fbserv reads FBSERV_TOKEN at start-up; by holding this
 *  mutex until the child has been forked we guarantee it sees the token
 *  that *this* worker intends. */
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
run_worker(int wid, int base_port, WorkerResult &result)
{
    result.worker_id = wid;

    int port = base_port + g_port_counter.fetch_add(1);
    char ipc_addr[MAXPATHLEN] = {0};
    if (g_use_ipc) {
	char hint[64] = {0};
	snprintf(hint, sizeof(hint), "fbserv-stress-%d-%d", bu_pid(), wid);
	if (pkg_ipc_addr(ipc_addr, sizeof(ipc_addr), hint) != 0) {
	    result.message = "failed to allocate an IPC listener address";
	    worker_log(wid, "ERROR: %s", result.message.c_str());
	    return;
	}
	worker_log(wid, "using IPC address %s", ipc_addr);
    } else {
	/* Scan forward until we find a free port (max 200 attempts). */
	for (int tries = 0; tries < 200 && port_in_use(port);
		tries++, port = base_port + g_port_counter.fetch_add(1))
	    ;
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
     *    fbserv -A -w 512 -n 512 -F /dev/mem -p <port>
     *
     *    NOTE: /dev/mem is the in-memory framebuffer (null device);
     *    the -F form launches a single-frame-buffer server which stays
     *    alive until killed, which is what rtwizard uses.
     * ------------------------------------------------------------------ */
    char port_str[16];
    snprintf(port_str, sizeof(port_str), "%d", port);

    const char *argv_fbserv_tcp[] = {
	g_fbserv_path,
        "-A",
        "-w", "512",
        "-n", "512",
        "-F", "/dev/mem",
        "-p", port_str,
        NULL
    };
    const char *argv_fbserv_ipc[] = {
	g_fbserv_path,
	"-A",
	"-w", "512",
	"-n", "512",
	"-F", "/dev/mem",
	"-I", ipc_addr,
	NULL
    };
    const char **argv_fbserv = g_use_ipc ? argv_fbserv_ipc : argv_fbserv_tcp;

    /* Pre-supply the token via the environment variable.
     * Lock the mutex so that concurrent workers don't stomp on each
     * other's FBSERV_TOKEN while the child process is being forked.
     *
     * setenv/unsetenv copy the string internally (unlike putenv), so
     * there is no dangling-pointer hazard. */
    struct bu_process *fbserv_proc = NULL;
    {
        std::lock_guard<std::mutex> lk(g_token_env_mutex);
#if defined(HAVE_UNISTD_H) && !defined(_WIN32)
        setenv("FBSERV_TOKEN", token, 1 /* overwrite */);
#else
        /* Windows fallback: _putenv_s makes its own copy. */
        _putenv_s("FBSERV_TOKEN", token);
#endif
        struct bu_process *fbserv_proc_tmp = NULL;
        bu_process_create(&fbserv_proc_tmp, argv_fbserv, BU_PROCESS_DEFAULT);
        fbserv_proc = fbserv_proc_tmp;

        /* Immediately clear so that other workers (or any subsequent code
         * in this worker) don't accidentally inherit this token. */
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
        return;
    }
    worker_log(wid, "fbserv launched (pid %d)", fbserv_pid);

    /* fbserv needs a moment to bind its socket.  Retrying the public open
     * also exercises cleanup of connections which fail before FBOPEN. */
    imgstream_fb_t *fb = NULL;
    int interval_ms = RETRY_INTERVAL_MS;
    int64_t deadline = bu_gettime() + BU_SEC2USEC(g_connect_timeout_sec);
    int attempt = 0;
    std::string remote_spec = g_use_ipc ? std::string("ipc:") + ipc_addr :
	std::string("localhost:") + port_str;
    struct imgstream_fb_remote_options options =
	IMGSTREAM_FB_REMOTE_OPTIONS_INIT;
    options.auth_token = token;

    while (!fb && bu_gettime() < deadline) {
        attempt++;
        fb = imgstream_fb_open_remote(remote_spec.c_str(), 512, 512,
		&options);
        if (!fb) {
	    if (!bu_process_alive(fbserv_proc)) {
		char child_error[1024] = {0};
		int error_length = bu_process_read_n(fbserv_proc,
		    BU_PROCESS_STDERR, (int)sizeof(child_error) - 1,
		    child_error);
		if (error_length > 0)
		    child_error[error_length] = '\0';
		int child_status = bu_process_wait_n(&fbserv_proc, 1);
		std::ostringstream ss;
		ss << "fbserv exited before accepting a connection (status "
		   << child_status << ")";
		if (error_length > 0)
		    ss << ": " << child_error;
		result.message = ss.str();
		break;
	    }
            worker_log(wid, "connect attempt %d failed, retrying in %dms...",
                       attempt, interval_ms);
            bu_snooze((int64_t)interval_ms * (int64_t)1000);
            interval_ms = (interval_ms * 2 < RETRY_MAX_MS) ? interval_ms * 2 : RETRY_MAX_MS;
        }
    }

    if (!fb) {
	if (result.message.empty())
	    result.message = "timed out connecting to fbserv";
        worker_log(wid, "ERROR: %s (tried %d times)", result.message.c_str(), attempt);
	if (fbserv_proc) {
	    bu_pid_terminate(fbserv_pid);
	    (void)bu_process_wait_n(&fbserv_proc, 1000);
	}
        return;
    }
    worker_log(wid, "connected after %d attempt(s)", attempt);

    bool operations_ok = imgstream_fb_width(fb) == 512 &&
	imgstream_fb_height(fb) == 512;
    const unsigned char clear_rgb[3] = {
	(unsigned char)(17 + wid), (unsigned char)(41 + wid),
	(unsigned char)(89 + wid)
    };
    unsigned char written[4 * 3 * 3];
    unsigned char readback[sizeof(written)];
    for (size_t i = 0; i < sizeof(written); i++)
	written[i] = (unsigned char)((i * 13 + wid * 7) & 0xff);
    memset(readback, 0, sizeof(readback));

    operations_ok = operations_ok && imgstream_fb_clear(fb, clear_rgb) == 0;
    operations_ok = operations_ok &&
	imgstream_fb_writerect(fb, 11, 13, 4, 3, written) == 12;
    operations_ok = operations_ok && imgstream_fb_flush(fb) == 0;
    operations_ok = operations_ok &&
	imgstream_fb_readrect(fb, 11, 13, 4, 3, readback) == 12;
    operations_ok = operations_ok &&
	memcmp(written, readback, sizeof(written)) == 0;

    int xcenter = 0, ycenter = 0, xzoom = 0, yzoom = 0;
    operations_ok = operations_ok && imgstream_fb_view(fb, 101, 103, 2, 3) == 0;
    operations_ok = operations_ok &&
	imgstream_fb_getview(fb, &xcenter, &ycenter, &xzoom, &yzoom) == 0 &&
	xcenter == 101 && ycenter == 103 && xzoom == 2 && yzoom == 3;

    int mode = 0, cursor_x = 0, cursor_y = 0;
    operations_ok = operations_ok && imgstream_fb_cursor(fb, 1, 23, 29) == 0;
    operations_ok = operations_ok &&
	imgstream_fb_getcursor(fb, &mode, &cursor_x, &cursor_y) == 0 &&
	mode == 1 && cursor_x == 23 && cursor_y == 29;

    struct imgstream_fb_colormap written_map;
    struct imgstream_fb_colormap read_map;
    imgstream_fb_colormap_linear(&written_map);
    memset(&read_map, 0, sizeof(read_map));
    written_map.red[17] = (uint16_t)(1000 + wid);
    written_map.green[17] = (uint16_t)(2000 + wid);
    written_map.blue[17] = (uint16_t)(3000 + wid);
    operations_ok = operations_ok && imgstream_fb_wmap(fb, &written_map) == 0;
    operations_ok = operations_ok && imgstream_fb_rmap(fb, &read_map) == 0;
    operations_ok = operations_ok &&
	memcmp(&written_map, &read_map, sizeof(written_map)) == 0;

    imgstream_fb_close(fb);
    bu_pid_terminate(fbserv_pid);
    (void)bu_process_wait_n(&fbserv_proc, 1000);

    if (!operations_ok) {
	result.message = "authenticated remote framebuffer operation failed";
	worker_log(wid, "ERROR: %s", result.message.c_str());
	return;
    }

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
        "  --base-port P    first port to try (default: 5600)\n"
        "  --ipc            use reusable local IPC listeners instead of TCP\n"
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

    for (int i = 1; i < argc; i++) {
        if (BU_STR_EQUAL(argv[i], "--workers") && i + 1 < argc) {
            num_workers = atoi(argv[++i]);
        } else if (BU_STR_EQUAL(argv[i], "--base-port") && i + 1 < argc) {
            base_port = atoi(argv[++i]);
        } else if (BU_STR_EQUAL(argv[i], "--timeout") && i + 1 < argc) {
            g_connect_timeout_sec = atoi(argv[++i]);
        } else if (BU_STR_EQUAL(argv[i], "--ipc")) {
	    g_use_ipc = true;
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
    bu_dir(g_fbserv_path, sizeof(g_fbserv_path), BU_DIR_BIN, "fbserv",
	BU_DIR_EXT, NULL);
    if (!g_fbserv_path[0]) {
	fprintf(stderr, "fbserv_stress: unable to locate fbserv\n");
	return 1;
    }

    fprintf(stderr, "fbserv_stress: launching %d parallel workers (%s)\n",
	    num_workers, g_use_ipc ? "IPC" : "TCP");

    std::vector<WorkerResult> results(num_workers);
    std::vector<std::thread> threads;
    threads.reserve(num_workers);

    for (int w = 0; w < num_workers; w++) {
        threads.emplace_back(run_worker, w, base_port, std::ref(results[w]));
    }

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
