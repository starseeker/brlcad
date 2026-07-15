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
/** @addtogroup libtclcad */
/** @{ */
/**
 *
 *  These are the Tcl specific callbacks used for I/O between client
 *  and server.  So far, experiments to simply use file handlers for
 *  Windows as well as other platforms haven't been successful, which
 *  is why we've had to accept the increased complexity of this code.
 */
/** @} */

#include "common.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <tcl.h>
#include "bio.h"
#include "bnetwork.h"
#include "imgstream/fbserv.h"
#include "pkg.h"
#include "tclcad.h"

/* We need to use different Tcl I/O mechanisms on different
 * platforms.  Make a define for that. */
#if defined(_WIN32) && !defined(__CYGWIN__)
#define USE_TCL_CHAN
#endif

/*
 * Communication error.  An error occurred on the PKG link.
 */
static void
comm_error(const char *str)
{
    bu_log("%s", str);
}

#ifdef USE_TCL_CHAN
/* Per SVN commit r29234, this function is needed because Windows uses
 * Tcl_OpenTcpServer rather than pkg_permserver to set up communication
 * on Windows (Tcl_OpenTcpServer's output is saved on fbsl_chan rather
 * than fbsl_fd.)  This function set up a pkg_conn structure and
 * assigns the fbsl_fd number for that scenario. */
static struct pkg_conn *
fbs_makeconn(int fd, const struct pkg_switch *switchp)
{
    struct pkg_conn *pc = pkg_adopt_socket(fd, switchp, comm_error);
    if (pc == PKC_ERROR) {
	comm_error("fbs_makeconn: pkg_adopt_socket failure\n");
    }
    return pc;
}
#endif

/*
 * Accept any new client connections.
 */
static void
#ifdef USE_TCL_CHAN
new_client_handler(ClientData clientData, Tcl_Channel chan, char *UNUSED(host), int UNUSED(port))
#else
new_client_handler(ClientData clientData, int UNUSED(port))
#endif
{
    struct fbserv_obj *fbsp = fbs_listener_owner(clientData);
    struct pkg_switch *pswitch = fbs_pkg_switch();
    void *cdata = NULL;
    struct pkg_conn *pcp = NULL;

#ifdef USE_TCL_CHAN
    uintptr_t pfd = (uintptr_t)fbs_listener_data_fd(clientData);
    if (Tcl_GetChannelHandle(chan, TCL_READABLE, (ClientData *)&pfd) != TCL_OK)
	return;
    cdata = (void *)chan;
    // TODO - we need the uintptr_t above, or we introduce stack corruption on
    // Windows.  Less certain about why we're casting back to int and using
    // it for fbs_makeconn...  presumably it's not really a valid pkc_fd, so
    // does it make sense to assign it there?
    pcp = fbs_makeconn((int)pfd, pswitch);
#else
    pcp = pkg_accept(fbs_listener_data_pkg_listener(clientData), pswitch, comm_error, 0);
#endif

    fbs_new_client(fbsp, pcp, cdata);
}

/* Check if we're already listening. */
C_DECL int
tclcad_is_listening(struct fbserv_obj *fbsp)
{
#ifdef USE_TCL_CHAN
    if (fbs_listener_channel(fbsp) != NULL) {
#else
    if (fbs_listener_fd(fbsp) >= 0) {
#endif
	return 1;
    }
    return 0;
}

C_DECL int
tclcad_listen_on_port(struct fbserv_obj *fbsp, int available_port)
{
    char hostname[32] = {0};
    /* XXX hardwired for now */
    sprintf(hostname, "localhost");

#ifdef USE_TCL_CHAN
    fbs_set_listener_channel(fbsp, Tcl_OpenTcpServer((Tcl_Interp *)fbsp->fbs_interp, available_port, hostname, new_client_handler, (ClientData)fbs_listener_handler_data(fbsp)));
    if (fbs_listener_channel(fbsp) == NULL) {
	/* This clobbers the result string which probably has junk
	 * related to the failed open.
	 */
	Tcl_DString ds;
	Tcl_DStringInit(&ds);
	Tcl_DStringResult((Tcl_Interp *)fbsp->fbs_interp, &ds);
	return 0;
    } else {
	return 1;
    }
#else
    char portname[32] = {0};
    sprintf(portname, "%d", available_port);
    struct pkg_listener *listener = pkg_listen(portname, NULL, 0, comm_error);
    fbs_set_listener_pkg_listener(fbsp, listener);
    if (listener) {
	fbs_set_listener_fd(fbsp, pkg_get_listener_fd(listener));
	return 1;
    }
#endif
    return 0;
}

C_DECL void
tclcad_open_server_handler(struct fbserv_obj *fbsp)
{
#ifdef USE_TCL_CHAN
    uintptr_t pfd = (uintptr_t)fbs_listener_fd(fbsp);
    Tcl_GetChannelHandle((Tcl_Channel)fbs_listener_channel(fbsp), TCL_READABLE, (ClientData *)&pfd);
    fbs_set_listener_fd(fbsp, (int)pfd);
#else
    Tcl_CreateFileHandler(fbs_listener_fd(fbsp), TCL_READABLE, (Tcl_FileProc *)new_client_handler, (ClientData)fbs_listener_handler_data(fbsp));
#endif
}

C_DECL void
tclcad_close_server_handler(struct fbserv_obj *fbsp)
{
#ifdef USE_TCL_CHAN
    Tcl_Channel chan = (Tcl_Channel)fbs_listener_channel(fbsp);
    if (chan != NULL) {
	Tcl_ChannelProc *callback = (Tcl_ChannelProc *)new_client_handler;
	Tcl_DeleteChannelHandler(chan, callback, (ClientData)(uintptr_t)fbs_listener_fd(fbsp));
	Tcl_Close((Tcl_Interp *)fbsp->fbs_interp, chan);
	fbs_set_listener_channel(fbsp, NULL);
    }
#else
    Tcl_DeleteFileHandler(fbs_listener_fd(fbsp));
#endif
}

C_DECL void
#ifdef USE_TCL_CHAN
tclcad_open_client_handler(struct fbserv_obj *fbsp, int i, void *data)
#else
tclcad_open_client_handler(struct fbserv_obj *fbsp, int i, void *UNUSED(data))
#endif
{
#ifdef USE_TCL_CHAN
    fbs_set_client_channel(fbsp, i, data);
    fbs_set_client_handler(fbsp, i, (void *)fbs_existing_client_handler);
    Tcl_CreateChannelHandler((Tcl_Channel)fbs_client_channel(fbsp, i), TCL_READABLE,
	    (Tcl_ChannelProc *)fbs_client_handler(fbsp, i), (ClientData)fbs_client_handler_data(fbsp, i));
#else
    Tcl_CreateFileHandler(fbs_client_fd(fbsp, i), TCL_READABLE,
	    fbs_existing_client_handler, (ClientData)fbs_client_handler_data(fbsp, i));
#endif
}

C_DECL void
tclcad_close_client_handler(struct fbserv_obj *fbsp, int sub)
{
#ifdef USE_TCL_CHAN
    Tcl_Channel chan = (Tcl_Channel)fbs_client_channel(fbsp, sub);
    Tcl_ChannelProc *handler = (Tcl_ChannelProc *)fbs_client_handler(fbsp, sub);
    Tcl_DeleteChannelHandler(chan, handler, (ClientData)(uintptr_t)fbs_client_fd(fbsp, sub));

    Tcl_Close((Tcl_Interp *)fbsp->fbs_interp, chan);
    fbs_set_client_channel(fbsp, sub, NULL);
#else
    Tcl_DeleteFileHandler(fbs_client_fd(fbsp, sub));
#endif
}


#ifdef USE_TCL_CHAN
/**
 * Poll callback for IPC clients on Windows (USE_TCL_CHAN path).
 *
 * Background: on Windows, Tcl_CreateFileHandler is a no-op, and
 * Tcl_MakeFileChannel for pipe handles starts a background reader thread
 * that CONSUMES data from the pipe into Tcl's internal buffer -- which
 * prevents pkg_process() from reading the data via the raw fd.
 *
 * Instead, we install a recurring Tcl timer that uses PeekNamedPipe()
 * (non-destructive) to check whether data is available.  When data is
 * found the same fbs_existing_client_handler() that POSIX file-handler
 * callbacks fire is called directly.  10 ms polling corresponds to the
 * 100 fps maximum framebuffer refresh rate; it adds negligible overhead
 * compared to the actual rendering time.
 */
static void
tcl_ipc_poll_win(ClientData cd)
{
    /* Termination guard: pkg was already closed, stop the timer. */
    if (!fbs_client_data_pkg(cd) || fbs_client_data_pkg(cd) == PKC_NULL) {
	fbs_set_client_data_channel(cd, NULL);
	return;
    }

    /* Peek at the pipe without consuming data. */
    {
	int fd = fbs_client_data_fd(cd);
	int has_data = 0;
	HANDLE h = (fd >= 0) ? (HANDLE)_get_osfhandle(fd) : INVALID_HANDLE_VALUE;
	if (h != INVALID_HANDLE_VALUE) {
	    DWORD avail = 0;
	    if (PeekNamedPipe(h, NULL, 0, NULL, &avail, NULL) && avail > 0)
		has_data = 1;
	}

	/* Reschedule BEFORE calling the handler: the handler may call
	 * drop_client() → tclcad_close_ipc_client_win() which cancels
	 * the new token, so we must store it first.                       */
	{
	    Tcl_TimerToken tok = Tcl_CreateTimerHandler(10, tcl_ipc_poll_win, cd);
	    fbs_set_client_data_channel(cd, (void *)tok);
	}

	if (has_data)
	    fbs_existing_client_handler(cd, TCL_READABLE);
    }
}

/**
 * Open handler for IPC clients on Windows (USE_TCL_CHAN path).
 * Installs the recurring poll timer instead of a file/channel handler.
 */
static void
tclcad_open_ipc_client_win(struct fbserv_obj *fbsp, int i, void *UNUSED(data))
{
    void *client_data = fbs_client_handler_data(fbsp, i);
    Tcl_TimerToken tok =
	Tcl_CreateTimerHandler(10, tcl_ipc_poll_win,
			       (ClientData)client_data);
    fbs_set_client_channel(fbsp, i, (void *)tok);
    fbs_set_client_handler(fbsp, i, NULL);
}

/**
 * Close handler for IPC clients on Windows (USE_TCL_CHAN path).
 * Cancels the poll timer installed by tclcad_open_ipc_client_win().
 */
static void
tclcad_close_ipc_client_win(struct fbserv_obj *fbsp, int sub)
{
    void *chan = fbs_client_channel(fbsp, sub);
    if (chan) {
	Tcl_DeleteTimerHandler((Tcl_TimerToken)chan);
	fbs_set_client_channel(fbsp, sub, NULL);
	fbs_set_client_handler(fbsp, sub, NULL);
    }
}
#endif /* USE_TCL_CHAN */


TCLCAD_EXPORT void
tclcad_fbserv_set_transport(struct fbserv_obj *fbsp)
{
    if (!fbsp)
	return;

#ifdef USE_TCL_CHAN
    static const struct fbserv_transport_ops ops = {
	tclcad_is_listening,
	tclcad_listen_on_port,
	tclcad_open_server_handler,
	tclcad_close_server_handler,
	tclcad_open_client_handler,
	tclcad_close_client_handler,
	tclcad_open_ipc_client_win,
	tclcad_close_ipc_client_win
    };
#else
    static const struct fbserv_transport_ops ops = {
	tclcad_is_listening,
	tclcad_listen_on_port,
	tclcad_open_server_handler,
	tclcad_close_server_handler,
	tclcad_open_client_handler,
	tclcad_close_client_handler,
	tclcad_open_client_handler,
	tclcad_close_client_handler
    };
#endif

    fbs_set_transport(fbsp, &ops);
}


/**
 * Phase 4: IPC listen path for libtclcad.
 *
 * Sets up an IPC-based (pipe/socketpair) framebuffer server on @p fbsp
 * without binding a TCP port.
 *
 * On POSIX, tclcad_open_client_handler / tclcad_close_client_handler use
 * Tcl_CreateFileHandler / Tcl_DeleteFileHandler which work for any fd.
 *
 * On Windows (USE_TCL_CHAN), Tcl_CreateFileHandler is a no-op.  Wrapping
 * a pipe HANDLE in a Tcl channel via Tcl_MakeFileChannel is also not viable
 * because Tcl's Windows pipe channel driver spawns a background reader
 * thread that CONSUMES data from the pipe into its own buffer, preventing
 * pkg_process() from reading via the raw fd.  Instead we install a
 * recurring Tcl timer (tclcad_open_ipc_client_win / tcl_ipc_poll_win) that
 * uses PeekNamedPipe() to check availability without consuming data, then
 * calls fbs_existing_client_handler() when data is ready.
 *
 * If the open succeeds and @p interp is non-NULL, the PKG_ADDR child-end
 * address string is stored in the Tcl variable fbserv(ipc_addr) so that
 * Tcl scripts (rtwizard, MGED's rt.tcl) can pass it to subprocesses.
 *
 * @return BRLCAD_OK on success, BRLCAD_ERROR on failure.
 */
TCLCAD_EXPORT int
tclcad_listen_ipc(struct fbserv_obj *fbsp, Tcl_Interp *interp)
{
    /* Endpoint-backed Obol streams provide a worker transport so synchronous
     * renderer processes do not depend on Tcl servicing their socket.  The
     * Tcl notifier remains the fallback for callers without one. */
    if (!fbs_can_open_ipc(fbsp))
	tclcad_fbserv_set_transport(fbsp);

    if (fbs_open_ipc(fbsp) != BRLCAD_OK) {
	bu_log("tclcad_listen_ipc: fbs_open_ipc failed\n");
	return BRLCAD_ERROR;
    }

    /* Phase 6 support: expose the child-end address as a Tcl variable so
     * rtwizard / MGED Tcl scripts can pass it to spawned subprocesses.     */
    if (interp) {
	const char *addr_env = fbs_ipc_child_addr_env(fbsp);
	if (addr_env) {
	    const char *eq = strchr(addr_env, '=');
	    const char *addr_val = eq ? eq + 1 : addr_env;
	    Tcl_SetVar2(interp, "fbserv", "ipc_addr", addr_val, TCL_GLOBAL_ONLY);
	}
    }

    return BRLCAD_OK;
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
