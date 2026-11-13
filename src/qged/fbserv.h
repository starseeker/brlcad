/*                        F B S E R V . H
 * BRL-CAD
 *
 * Copyright (c) 2021-2026 United States Government as represented by
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
/** @file fbserv.h
 *
 * Brief description
 *
 */

#ifndef QGED_FBSERV_H
#define QGED_FBSERV_H

#include "common.h"

#include <QByteArray>
#include <QHostAddress>
#include <QObject>
#include <QSocketNotifier>
#include <QTcpServer>
#include <QTcpSocket>
#include <QTimer>
#include <iostream>

struct ged;
class QgView;

// Per client info (TCP path)
class QFBSocket : public QObject
{
    Q_OBJECT

    public:
	QTcpSocket *s;
	int ind;
	void *fbsp;

    signals:
	void updated();

    public slots:
	void client_handler();
	void on_disconnected();

    private:
        QByteArray buff;
};

// Per client info (IPC path — pipe/socketpair instead of TCP)
class QFBIPCSocket : public QObject
{
    Q_OBJECT

    public:
	QSocketNotifier *notifier = nullptr;
	QTimer *timer = nullptr;
	int ind = -1;
	void *fbsp = nullptr;

    signals:
	void updated();

    public slots:
	void ipc_handler();
};

// Overall server that sets up clients
// in response to connection requests
class QFBServer : public QTcpServer
{
	Q_OBJECT

    public:
	QFBServer(void *fp = nullptr);
	~QFBServer();

	int port = -1;
	void *fbsp;

    public slots:
	void on_Connect();
};


extern void
qged_fbserv_configure_ged_handlers(struct ged *gedp, QgView *display);

/** Clear QGED's non-owning Qt transport target before GED tears down fbserv. */
extern void
qged_fbserv_release_ged_handlers(struct ged *gedp);

#endif /* QGED_FBSERV_H */

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
