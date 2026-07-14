/*                        F B S E R V . C
 * BRL-CAD
 *
 * Copyright (c) 1995-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file mged/fbserv.c
 *
 * MGED host wiring for libimgstream's framebuffer server.  Protocol and
 * connection lifecycle behavior belongs to libimgstream; Tcl notification
 * belongs to libtclcad.
 */

#include "common.h"

#include "bu/log.h"
#include "ged/draw_obol.h"
#include "imgstream/fbserv.h"
#include "tclcad.h"

#include "./mged.h"
#include "./mged_display.h"

void
fbserv_set_port(const struct bu_structparse *UNUSED(sp),
	const char *UNUSED(c1), void *UNUSED(v1),
	const char *UNUSED(c2), void *UNUSED(v2))
{
    struct mged_state *s = MGED_STATE;
    struct fbserv_obj *fbsp;

    if (!s || !s->gedp || !s->mged_curr_display)
	return;
    fbsp = s->gedp->ged_fbs;
    if (!fbsp)
	return;

    if (fbs_can_close(fbsp) && fbs_listener_fd(fbsp) >= 0)
	(void)fbs_close(fbsp);

    if (!mged_variables->mv_listen)
	return;
    if (!mged_variables->mv_fb || !mged_obol_framebuffer_ensure(s) ||
	!fbs_framebuffer_backend_installed(fbsp)) {
	mged_variables->mv_listen = 0;
	return;
    }

    fbsp->fbs_interp = s->interp;
    tclcad_fbserv_set_transport(fbsp);

    (void)fbs_generate_token(fbsp);
    if (fbs_open(fbsp, mged_variables->mv_port) != BRLCAD_OK) {
	bu_log("fbserv_set_port: failed to open a framebuffer listener\n");
	mged_variables->mv_listen = 0;
	return;
    }

    mged_variables->mv_port = fbs_listener_port(fbsp);
}
