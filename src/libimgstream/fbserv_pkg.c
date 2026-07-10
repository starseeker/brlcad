/*                    F B S E R V _ P K G . C
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
/** @file libimgstream/fbserv_pkg.c */

#include "common.h"

#include <stdlib.h>

#include "bu/defines.h"
#include "bu/log.h"
#include "imgstream/fbserv.h"
#include "pkg.h"


static void
fbserv_pkg_missing_handler(struct pkg_conn *pcp, char *buf)
{
    if (pcp)
	bu_log("fbserv: no handler installed for message type %d\n",
	    pcp->pkc_type);
    if (buf)
	(void)free(buf);
}


#define FBSERV_HANDLER(_handlers, _name) \
    ((_handlers) && (_handlers)->_name ? (_handlers)->_name : fbserv_pkg_missing_handler)


static void
fbserv_pkg_switch_set(struct pkg_switch *entry,
	unsigned short type,
	pkg_callback handler,
	const char *title)
{
    entry->pks_type = type;
    entry->pks_handler = handler;
    entry->pks_title = title;
    entry->pks_user_data = NULL;
}


int
fbserv_pkg_switch_init(struct pkg_switch *pswitch,
	size_t switch_count,
	const struct fbserv_pkg_handlers *handlers)
{
    if (!pswitch || switch_count < FBSERV_PKG_SWITCH_COUNT)
	return BRLCAD_ERROR;

    fbserv_pkg_switch_set(&pswitch[0], FBSERV_MSG_FBAUTH,
	FBSERV_HANDLER(handlers, auth), "Session Authentication");
    fbserv_pkg_switch_set(&pswitch[1], FBSERV_MSG_FBOPEN,
	FBSERV_HANDLER(handlers, open), "Open Framebuffer");
    fbserv_pkg_switch_set(&pswitch[2], FBSERV_MSG_FBCLOSE,
	FBSERV_HANDLER(handlers, close), "Close Framebuffer");
    fbserv_pkg_switch_set(&pswitch[3], FBSERV_MSG_FBCLEAR,
	FBSERV_HANDLER(handlers, clear), "Clear Framebuffer");
    fbserv_pkg_switch_set(&pswitch[4], FBSERV_MSG_FBREAD,
	FBSERV_HANDLER(handlers, read), "Read Pixels");
    fbserv_pkg_switch_set(&pswitch[5], FBSERV_MSG_FBWRITE,
	FBSERV_HANDLER(handlers, write), "Write Pixels");
    fbserv_pkg_switch_set(&pswitch[6],
	FBSERV_MSG_FBWRITE + FBSERV_MSG_NORETURN,
	FBSERV_HANDLER(handlers, write), "Asynch write");
    fbserv_pkg_switch_set(&pswitch[7], FBSERV_MSG_FBCURSOR,
	FBSERV_HANDLER(handlers, cursor), "Cursor");
    fbserv_pkg_switch_set(&pswitch[8], FBSERV_MSG_FBGETCURSOR,
	FBSERV_HANDLER(handlers, getcursor), "Get Cursor");
    fbserv_pkg_switch_set(&pswitch[9], FBSERV_MSG_FBSCURSOR,
	FBSERV_HANDLER(handlers, scursor), "Screen Cursor");
    fbserv_pkg_switch_set(&pswitch[10], FBSERV_MSG_FBWINDOW,
	FBSERV_HANDLER(handlers, window), "Window");
    fbserv_pkg_switch_set(&pswitch[11], FBSERV_MSG_FBZOOM,
	FBSERV_HANDLER(handlers, zoom), "Zoom");
    fbserv_pkg_switch_set(&pswitch[12], FBSERV_MSG_FBVIEW,
	FBSERV_HANDLER(handlers, view), "View");
    fbserv_pkg_switch_set(&pswitch[13], FBSERV_MSG_FBGETVIEW,
	FBSERV_HANDLER(handlers, getview), "Get View");
    fbserv_pkg_switch_set(&pswitch[14], FBSERV_MSG_FBRMAP,
	FBSERV_HANDLER(handlers, rmap), "R Map");
    fbserv_pkg_switch_set(&pswitch[15], FBSERV_MSG_FBWMAP,
	FBSERV_HANDLER(handlers, wmap), "W Map");
    fbserv_pkg_switch_set(&pswitch[16], FBSERV_MSG_FBHELP,
	FBSERV_HANDLER(handlers, help), "Help Request");
    fbserv_pkg_switch_set(&pswitch[17], FBSERV_MSG_ERROR,
	FBSERV_HANDLER(handlers, unknown), "Error Message");
    fbserv_pkg_switch_set(&pswitch[18], FBSERV_MSG_CLOSE,
	FBSERV_HANDLER(handlers, unknown), "Close Connection");
    fbserv_pkg_switch_set(&pswitch[19], FBSERV_MSG_FBREADRECT,
	FBSERV_HANDLER(handlers, readrect), "Read Rectangle");
    fbserv_pkg_switch_set(&pswitch[20], FBSERV_MSG_FBWRITERECT,
	FBSERV_HANDLER(handlers, writerect), "Write Rectangle");
    fbserv_pkg_switch_set(&pswitch[21],
	FBSERV_MSG_FBWRITERECT + FBSERV_MSG_NORETURN,
	FBSERV_HANDLER(handlers, writerect), "Write Rectangle");
    fbserv_pkg_switch_set(&pswitch[22], FBSERV_MSG_FBBWREADRECT,
	FBSERV_HANDLER(handlers, bwreadrect), "Read BW Rectangle");
    fbserv_pkg_switch_set(&pswitch[23], FBSERV_MSG_FBBWWRITERECT,
	FBSERV_HANDLER(handlers, bwwriterect), "Write BW Rectangle");
    fbserv_pkg_switch_set(&pswitch[24],
	FBSERV_MSG_FBBWWRITERECT + FBSERV_MSG_NORETURN,
	FBSERV_HANDLER(handlers, bwwriterect), "Write BW Rectangle");
    fbserv_pkg_switch_set(&pswitch[25], FBSERV_MSG_FBFLUSH,
	FBSERV_HANDLER(handlers, flush), "Flush Output");
    fbserv_pkg_switch_set(&pswitch[26],
	FBSERV_MSG_FBFLUSH + FBSERV_MSG_NORETURN,
	FBSERV_HANDLER(handlers, flush), "Flush Output");
    fbserv_pkg_switch_set(&pswitch[27], FBSERV_MSG_FBFREE,
	FBSERV_HANDLER(handlers, free), "Free Resources");
    fbserv_pkg_switch_set(&pswitch[28], FBSERV_MSG_FBPOLL,
	FBSERV_HANDLER(handlers, poll), "Handle Events");
    fbserv_pkg_switch_set(&pswitch[29], FBSERV_MSG_FBSETCURSOR,
	FBSERV_HANDLER(handlers, setcursor), "Set Cursor Shape");
    fbserv_pkg_switch_set(&pswitch[30],
	FBSERV_MSG_FBSETCURSOR + FBSERV_MSG_NORETURN,
	FBSERV_HANDLER(handlers, setcursor), "Set Cursor Shape");
    fbserv_pkg_switch_set(&pswitch[31], 0, NULL, NULL);

    return BRLCAD_OK;
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
