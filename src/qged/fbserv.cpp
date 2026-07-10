/*                      F B S E R V . C P P
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
/** @file fbserv.cpp
 *
 *  These are the Qt specific callbacks used for I/O between client
 *  and server.
 *
 *  TODO - Look into QLocalSocket, and whether we might be able to
 *  generalize libpkg (or even just use parts of it) to allow us
 *  to communicate using that mechanism...
 *
 *  Initial thought - optional callback functions to replace
 *  select, read, etc - if not set default to current behavior,
 *  if set do the callback instead of those calls...
 */

#include "common.h"

/* bu/ipc.h removed - transport handled by libpkg */
#include "bu/log.h"
#include "imgstream/fbserv.h"
#include "pkg.h"
#include "ged.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "./fbserv.h"
#include "qtcad/QgCanvasBase.h"
#include "qtcad/QgObolWindowHost.h"
#include "qtcad/QgTypes.h"
#include "qtcad/QgView.h"
#include <QSize>
#include <QWidget>

#include <algorithm>
#include <cmath>
#include <cstring>

namespace {
static const qint64 DRAIN_TIME_BUDGET_MS = 4;
static const size_t DRAIN_BYTES_BUDGET = 512 * 1024;
static const int DRAIN_ITERATION_CAP = 256;

static struct fbserv_obj *
qdm_fbserv_obj(void *fbsp)
{
    return static_cast<struct fbserv_obj *>(fbsp);
}

class QDMObolFramebufferBridge {
public:
    explicit QDMObolFramebufferBridge(QgView *v = nullptr)
    {
	setDisplay(v);
    }

    void setDisplay(QgView *v)
    {
	display = v;
	host.setCanvas(display ? display->canvasBase() : nullptr);
    }

    int configure(struct ged *new_gedp, QgView *v)
    {
	if (!new_gedp || !v || !v->canvasBase() || !v->obolViewController())
	    return -1;

	/* Detach image nodes before moving the persistent host to another view. */
	if (gedp && (gedp != new_gedp || display != v))
	    ged_draw_obol_framebuffer_release(gedp);

	gedp = new_gedp;
	setDisplay(v);
	view_ctx = ged_view_active_ctx(gedp);
	if (!view_ctx)
	    return -1;

	if (!ged_draw_obol_controller_attach_opaque_for_view(gedp,
		view_ctx, v->obolViewController(), 1))
	    return -1;

	QSize size = renderSize();
	return ged_draw_obol_framebuffer_backend_install_for_view(gedp,
	    view_ctx, &host, size.width(), size.height(), 1);
    }

    void notifyUpdated()
    {
	if (display) {
	    if (gedp)
		(void)ged_draw_obol_framebuffer_present(gedp);
	    display->need_update(QG_VIEW_REFRESH);
	}
    }

    int matches(const struct fbserv_obj *fbsp) const
    {
	return gedp && gedp->ged_fbs == fbsp;
    }

private:
    QSize renderSize() const
    {
	QgCanvasBase *canvas = display ? display->canvasBase() : nullptr;
	QWidget *widget = canvas ? canvas->canvasWidget() : nullptr;
	if (!widget)
	    return QSize(512, 512);
	qreal dpr = widget->devicePixelRatioF();
	int rw = std::max(1, (int)std::ceil(widget->width() * dpr));
	int rh = std::max(1, (int)std::ceil(widget->height() * dpr));
	return QSize(rw, rh);
    }

    struct ged *gedp = nullptr;
    void *view_ctx = nullptr;
    QgView *display = nullptr;
    QgObolWindowHost host;
};

static QDMObolFramebufferBridge *qdm_active_obol_bridge = nullptr;

static QDMObolFramebufferBridge *
qdm_obol_bridge(struct fbserv_obj *fbsp)
{
    return qdm_active_obol_bridge && qdm_active_obol_bridge->matches(fbsp) ?
	qdm_active_obol_bridge : nullptr;
}

}

static int qdm_is_listening(struct fbserv_obj *fbsp);
static int qdm_listen_on_port(struct fbserv_obj *fbsp, int available_port);
static void qdm_open_server_handler(struct fbserv_obj *fbsp);
static void qdm_close_server_handler(struct fbserv_obj *fbsp);
static void qdm_open_client_handler(struct fbserv_obj *fbsp, int i,
	void *data);
static void qdm_open_ipc_client_handler(struct fbserv_obj *fbsp, int i,
	void *data);
static void qdm_close_client_handler(struct fbserv_obj *fbsp, int i);
static void qdm_close_ipc_client_handler(struct fbserv_obj *fbsp, int i);

static const struct fbserv_transport_ops qdm_transport_ops = {
    qdm_is_listening,
    qdm_listen_on_port,
    qdm_open_server_handler,
    qdm_close_server_handler,
    qdm_open_client_handler,
    qdm_close_client_handler,
    qdm_open_ipc_client_handler,
    qdm_close_ipc_client_handler
};

void
QFBSocket::client_handler()
{
    QTCAD_SLOT("QFBSocket::client_handler", 1);
    bu_log("client_handler\n");

    struct fbserv_obj *fbs = qdm_fbserv_obj(fbsp);

    /* If our client slot has already been torn down (e.g. socket
     * disconnected, drop_client called), there is nothing to do. */
    if (!fbs_client_active(fbs, ind))
	return;
    struct pkg_conn *pkc = fbs_client_pkg(fbs, ind);
    if (!pkc || !s)
	return;

    // Read data.  NOTE:  we're using the Qt read routines rather than
    // pkg_suckin, so the legacy transport helper cannot be used here.
    // Initially tried pkg_suckin, but it didn't seem to work with the socket
    // as set up by Qt.
    QByteArray dbuff = s->read(s->bytesAvailable());

    // We may not have processed all the read data last time, so append
    // this to anything left over from before
    buff.append(dbuff);

    // If we don't have anything, we're done
    if (!buff.length())
	return;

    unsigned char *remaining = nullptr;
    size_t remaining_size = 0;
    struct fbserv_process_result result;
    int process_status = fbs_process_client_bytes(fbs, ind,
	    reinterpret_cast<const unsigned char *>(buff.constData()),
	    static_cast<size_t>(buff.size()), &remaining, &remaining_size,
	    &result);
    buff.clear();
    if (remaining && remaining_size)
	buff.append(reinterpret_cast<const char *>(remaining),
	    static_cast<int>(remaining_size));
    free(remaining);

    if (process_status != BRLCAD_OK)
	bu_log("client_handler framebuffer protocol processing failed\n");
    if (result.messages_processed > 0)
	emit updated();
}

/* Phase D2 (ert reliability): when the remote rt closes the TCP
 * connection, tear down the fbserv client slot cleanly so the
 * pkg_conn / fd / chan are released and a follow-on rt session can
 * occupy the slot.  Without this, a finished rt leaves an orphan
 * QFBSocket and pkg_conn alive until application exit. */
void
QFBSocket::on_disconnected()
{
    QTCAD_SLOT("QFBSocket::on_disconnected", 1);
    struct fbserv_obj *fbs = qdm_fbserv_obj(fbsp);
    if (!fbs_client_active(fbs, ind))
	return;	/* already dropped */
    fbs_drop_client(fbs, ind);
}


QFBServer::QFBServer(void *fp)
{
    fbsp = fp;
}

QFBServer::~QFBServer()
{
}

void
QFBServer::on_Connect()
{
    QTCAD_SLOT("QFBServer::on_Connect", 1);
    // Have a new connection pending, accept it.
    QTcpSocket *tcps = nextPendingConnection();

    bu_log("new connection");

    QFBSocket *fs = new QFBSocket;
    fs->s = tcps;
    fs->fbsp = fbsp;

    int fd = tcps->socketDescriptor();
    bu_log("fd: %d\n", fd);
    struct pkg_conn *pc = pkg_adopt_socket(fd, fbs_pkg_switch(), 0);
    if (pc == PKC_ERROR) {
	bu_log("new connection failed (pkg_adopt_socket)");
	tcps->close();
	delete fs;
	return;
    }

    fs->ind = fbs_new_client(qdm_fbserv_obj(fbsp), pc, (void *)fs);
    if (fs->ind == -1) {
	bu_log("new connection failed");
	pkg_close(pc);
	tcps->close();
	delete fs;
    }
}

/* Check if we're already listening. */
static int
qdm_is_listening(struct fbserv_obj *fbsp)
{
    bu_log("is_listening\n");
    if (fbs_listener_fd(fbsp) >= 0) {
	return 1;
    }
    return 0;
}

static int
qdm_listen_on_port(struct fbserv_obj *fbsp, int available_port)
{
    bu_log("listen on port\n");
    QFBServer *nl = new QFBServer(fbsp);
    nl->port = available_port;
    if (!nl->listen(QHostAddress::LocalHost, available_port)) {
	bu_log("Failed to start listening on %d\n", available_port);
	delete nl;
	return 0;
    }
    fbs_set_listener_channel(fbsp, (void *)nl);
    fbs_set_listener_fd(fbsp, nl->socketDescriptor());
    if (fbs_listener_fd(fbsp) >= 0)
	return 1;
    return 0;
}

static void
qdm_open_server_handler(struct fbserv_obj *fbsp)
{
    bu_log("open_server_handler\n");
    QFBServer *nl = (QFBServer *)fbs_listener_channel(fbsp);
    if (!nl->isListening())
	bu_log("not listening!\n");
    QObject::connect(nl, &QTcpServer::newConnection, nl, &QFBServer::on_Connect, Qt::QueuedConnection);
}

static void
qdm_close_server_handler(struct fbserv_obj *fbsp)
{
    bu_log("close_server_handler\n");
    QFBServer *nl = (QFBServer *)fbs_listener_channel(fbsp);
    delete nl;
    fbs_set_listener_channel(fbsp, NULL);
}

static void
qdm_open_client_handler(struct fbserv_obj *fbsp, int i, void *data)
{
    bu_log("open_client_handler\n");
    fbs_set_client_channel(fbsp, i, data);
    QFBSocket *s = (QFBSocket *)data;
    QObject::connect(s->s, &QTcpSocket::readyRead, s, &QFBSocket::client_handler, Qt::QueuedConnection);
    /* Phase D2: tear down on remote disconnect. */
    QObject::connect(s->s, &QTcpSocket::disconnected, s, &QFBSocket::on_disconnected, Qt::QueuedConnection);

    if (QDMObolFramebufferBridge *bridge = qdm_obol_bridge(fbsp))
	QObject::connect(s, &QFBSocket::updated, s,
			 [bridge]() { bridge->notifyUpdated(); },
			 Qt::QueuedConnection);
}

static void
qdm_close_client_handler(struct fbserv_obj *fbsp, int i)
{
    bu_log("close_client_handler\n");
    QFBSocket *s = (QFBSocket *)fbs_client_channel(fbsp, i);
    if (s) {
	/* Phase D2 (ert reliability): use deleteLater() to be safe when
	 * called from within a slot of `s` itself (e.g. on_disconnected). */
	s->deleteLater();
    }
    fbs_set_client_channel(fbsp, i, NULL);
}

void
qdm_configure_ged_fbserv_handlers(struct ged *gedp, QgView *display)
{
    if (!gedp || !gedp->ged_fbs)
	return;

    static QDMObolFramebufferBridge *bridge = nullptr;
    if (!bridge)
	bridge = new QDMObolFramebufferBridge(display);
    fbs_set_legacy_framebuffer(gedp->ged_fbs, NULL);
    if (bridge->configure(gedp, display) != 0) {
	bu_log("qged fbserv: unable to install shared Obol framebuffer backend\n");
	return;
    }
    qdm_active_obol_bridge = bridge;

    fbs_set_transport(gedp->ged_fbs, &qdm_transport_ops);
    gedp->ged_fbs->fbs_callback = (void (*)(void *))FBS_CALLBACK_NULL;
    gedp->ged_fbs->fbs_clientData = NULL;
}

/* -----------------------------------------------------------------------
 * Phase 5: IPC client handler for qged (pipe/socketpair instead of TCP).
 *
 * QFBIPCSocket registers the pre-connected raw file descriptor from libpkg
 * ipc with Qt's event loop via QSocketNotifier.  When the fd is readable,
 * ipc_handler() reads one message frame from the pkg_conn and dispatches it
 * via pkg_process().
 *
 * This avoids QTcpSocket entirely for the local rt→framebuffer data path.
 * ----------------------------------------------------------------------- */

void
QFBIPCSocket::ipc_handler()
{
    QTCAD_SLOT("QFBIPCSocket::ipc_handler", 1);

    struct fbserv_obj *fbs = qdm_fbserv_obj(fbsp);

    if (!fbs_client_active(fbs, ind))
	return;

    struct pkg_conn *pkc = fbs_client_pkg(fbs, ind);
    if (!pkc) {
	if (notifier) notifier->setEnabled(false);
	return;
    }

    /* Phase H2 (ert reliability): pkc_server_data is set once at slot
     * registration time in fbs_open_ipc(); rewriting it on every event
     * is at best redundant and at worst racy if a handler invoked from
     * pkg_process() recurses through the event loop.  Do not re-stamp
     * it here. */

    /* Phase C3 (ert reliability): drain the fd in a tight loop while
     * notifications are disabled.  Qt's QSocketNotifier is level-
     * triggered: with a slow GUI thread it can re-enqueue activated()
     * events faster than we can service them, and a slot scheduled
     * for an already-drained fd would otherwise call read() and (on a
     * blocking fd) hang the GUI thread inside _pkg_io_read.  We rely
     * on pkg_suckin's EAGAIN handling (Phase C2) to break the loop
     * cleanly when the kernel buffer is empty. */
    if (notifier)
	notifier->setEnabled(false);

    struct fbserv_process_result result;
    int process_status = fbs_drain_client(fbs, ind, DRAIN_BYTES_BUDGET,
	DRAIN_ITERATION_CAP, DRAIN_TIME_BUDGET_MS * 1000, &result);
    if (result.messages_processed > 0)
	emit updated();
    if (process_status != BRLCAD_OK)
	bu_log("QFBIPCSocket::ipc_handler: framebuffer protocol processing failed\n");
    if (result.disconnected)
	return;

    if (notifier)
	notifier->setEnabled(true);
}


static void
qdm_open_ipc_client_handler(struct fbserv_obj *fbsp, int i, void *UNUSED(data))
{
    bu_log("open_ipc_client_handler\n");

    QFBIPCSocket *s = new QFBIPCSocket;
    s->ind  = i;
    s->fbsp = fbsp;
    s->notifier = new QSocketNotifier(fbs_client_fd(fbsp, i),
				      QSocketNotifier::Read, s);
    fbs_set_client_channel(fbsp, i, (void *)s);

    QObject::connect(s->notifier, &QSocketNotifier::activated,
		     s, &QFBIPCSocket::ipc_handler, Qt::QueuedConnection);

    if (QDMObolFramebufferBridge *bridge = qdm_obol_bridge(fbsp))
	QObject::connect(s, &QFBIPCSocket::updated, s,
			 [bridge]() { bridge->notifyUpdated(); },
			 Qt::QueuedConnection);
}

static void
qdm_close_ipc_client_handler(struct fbserv_obj *fbsp, int i)
{
    bu_log("close_ipc_client_handler\n");
    QFBIPCSocket *s = (QFBIPCSocket *)fbs_client_channel(fbsp, i);
    if (s) {
	/* Phase D1 (ert reliability): use deleteLater() so the close
	 * handler is safe to call from within a slot of `s` itself
	 * (e.g. ipc_handler() detecting EOF and asking for a drop).
	 * Disable the notifier first so no further activations fire
	 * after the QObject pointer is otherwise considered dead. */
	if (s->notifier)
	    s->notifier->setEnabled(false);
	s->deleteLater();
    }
    fbs_set_client_channel(fbsp, i, NULL);
}


// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
