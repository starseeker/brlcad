/*                      Q G V I E W C T R L . C P P
 * BRL-CAD
 *
 * Copyright (c) 2022-2026 United States Government as represented by
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
/** @file QgViewCtrl.cpp
 *
 * Qt BRL-CAD View Control button panel implementation
 *
 */

#include "common.h"

#include "BObol/BDisplayEndpoint.h"
#include "bv.h"
#include "bu/env.h"
#include "ged.h"
#include "qtcad/QgSession.h"
#include "qtcad/QgViewCtrl.h"
#include "qtcad/QgSignalFlags.h"

static struct bv_context *
qgviewctrl_active_view(const QgViewCtrl *ctrl)
{
	QgSession *session = ctrl ? ctrl->session() : nullptr;
	return session ? session->activeViewContext() : nullptr;
}

static struct ged *
qgviewctrl_ged(const QgViewCtrl *ctrl)
{
	QgSession *session = ctrl ? ctrl->session() : nullptr;
	return session ? session->ged() : nullptr;
}

QgViewCtrl::QgViewCtrl(QWidget *pparent, QgSession *session) : QToolBar(pparent)
{
	m_session = session;

	this->setStyleSheet("QToolButton{margin:0px;}");

	sca = addAction(QIcon(QPixmap(":images/view/view_scale.png")), "Scale");
	rot = addAction(QIcon(QPixmap(":images/view/view_rotate.png")), "Rotate");
	tra = addAction(QIcon(QPixmap(":images/view/view_translate.png")), "Translate");
	center = addAction(QIcon(QPixmap(":images/view/view_center.png")), "Center");

	addSeparator();

	raytrace = addAction(QIcon(QPixmap(":images/view/raytrace.png")), "Raytrace");
	fb_mode = addAction(QIcon(QPixmap(":images/view/framebuffer_off.png")), "Framebuffer Off");
	fb_clear = addAction(QIcon(QPixmap(":images/view/framebuffer_clear.png")), "Clear Framebuffer");

	// Connect buttons to standard actions
	connect(sca, &QAction::triggered, this, &QgViewCtrl::sca_mode);
	connect(rot, &QAction::triggered, this, &QgViewCtrl::rot_mode);
	connect(tra, &QAction::triggered, this, &QgViewCtrl::tra_mode);
	connect(center, &QAction::triggered, this, &QgViewCtrl::center_mode);
	connect(raytrace, &QAction::triggered, this, &QgViewCtrl::raytrace_cmd);
	connect(fb_clear, &QAction::triggered, this, &QgViewCtrl::fbclear_cmd);
	connect(fb_mode, &QAction::triggered, this, &QgViewCtrl::fb_mode_cmd);
}

QgViewCtrl::~QgViewCtrl()
{
}

QgSession *
QgViewCtrl::session() const
{
	return m_session;
}

void
QgViewCtrl::setSession(QgSession *session)
{
	m_session = session;
}

void
QgViewCtrl::sca_mode()
{
	QTCAD_SLOT("QgViewCtrl::sca_mode", 1);
	emit lmouse_mode(BV_ADJUST_SCALE);
}

void
QgViewCtrl::rot_mode()
{
	QTCAD_SLOT("QgViewCtrl::rot_mode", 1);
	emit lmouse_mode(BV_ADJUST_ROT);
}

void
QgViewCtrl::tra_mode()
{
	QTCAD_SLOT("QgViewCtrl::tra_mode", 1);
	emit lmouse_mode(BV_ADJUST_TRANS);
}

void
QgViewCtrl::center_mode()
{
	QTCAD_SLOT("QgViewCtrl::center_mode", 1);
	emit lmouse_mode(BV_ADJUST_CENTER);
}


void
QgViewCtrl::fbclear_cmd()
{
	QTCAD_SLOT("QgViewCtrl::fbclear_cmd", 1);
	struct ged *gedp = qgviewctrl_ged(this);
	if (!gedp)
		return;
	const char *av[2] = {nullptr};
	av[0] = "fbclear";
	ged_exec_fbclear(gedp, 1, (const char **)av);
	emit view_changed(QG_VIEW_REFRESH);
}

void
QgViewCtrl::fb_mode_cmd()
{
	QTCAD_SLOT("QgViewCtrl::fb_mode_cmd", 1);
	struct bv_context *view_ctx = qgviewctrl_active_view(this);
	if (!view_ctx)
		return;
	struct bobol_endpoint_property_value value =
	    BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	if (ged_view_context_display_property_get(
		ged_view_context_from_bv(view_ctx),
		"composition.framebuffer.mode", &value) !=
	    BOBOL_ENDPOINT_PROPERTY_OK)
		return;
	value.type = BOBOL_ENDPOINT_PROPERTY_ENUM;
	if (BU_STR_EQUAL(value.string_value, "off"))
		value.string_value = "underlay";
	else if (BU_STR_EQUAL(value.string_value, "underlay"))
		value.string_value = "interlay";
	else if (BU_STR_EQUAL(value.string_value, "interlay"))
		value.string_value = "overlay";
	else if (BU_STR_EQUAL(value.string_value, "overlay"))
		value.string_value = "off";
	else
		return;
	if (ged_view_context_display_property_set(
		ged_view_context_from_bv(view_ctx),
		"composition.framebuffer.mode", &value) !=
	    BOBOL_ENDPOINT_PROPERTY_OK)
		return;
	emit view_changed(QG_VIEW_REFRESH);
}

void
QgViewCtrl::do_view_update(QgViewUpdateFlags flags)
{
	QTCAD_SLOT("QgViewCtrl::do_view_update", 1);
	struct bv_context *view_ctx = qgviewctrl_active_view(this);
	if (!view_ctx || !flags)
		return;
	struct bobol_endpoint_property_value value =
	    BOBOL_ENDPOINT_PROPERTY_VALUE_INIT;
	if (ged_view_context_display_property_get(
		ged_view_context_from_bv(view_ctx),
		"composition.framebuffer.mode", &value) !=
	    BOBOL_ENDPOINT_PROPERTY_OK)
		return;
	if (BU_STR_EQUAL(value.string_value, "off")) {
		fb_mode->setIcon(QIcon(QPixmap(":images/view/framebuffer_off.png")));
		fb_mode->setText("Framebuffer Off");
	} else if (BU_STR_EQUAL(value.string_value, "overlay")) {
		fb_mode->setIcon(QIcon(QPixmap(":images/view/framebuffer.png")));
		fb_mode->setText("Framebuffer Overlay");
	} else if (BU_STR_EQUAL(value.string_value, "underlay")) {
		fb_mode->setIcon(QIcon(QPixmap(":images/view/framebuffer_underlay.png")));
		fb_mode->setText("Framebuffer Underlay");
	} else if (BU_STR_EQUAL(value.string_value, "interlay")) {
		/* Reuse the overlay glyph until an interlay-specific Qt asset exists. */
		fb_mode->setIcon(QIcon(QPixmap(":images/view/framebuffer.png")));
		fb_mode->setText("Framebuffer Interlay");
	}
}

int
rt_cmd_start(int UNUSED(ac), const char **UNUSED(av), void *pid_ptr, void *ctx)
{
	int *pidp = (int *)pid_ptr;
	if (!pidp)
		return BRLCAD_OK;
	QgViewCtrl *vctrl = (QgViewCtrl *)ctx;
	vctrl->raytrace_start(*pidp);
	return BRLCAD_OK;
}

int
rt_cmd_done(int UNUSED(ac), const char **UNUSED(av), void *UNUSED(u1), void *ctx)
{
	QgViewCtrl *vctrl = (QgViewCtrl *)ctx;
	vctrl->raytrace_done();
	return BRLCAD_OK;
}

void
QgViewCtrl::raytrace_cmd()
{
	QTCAD_SLOT("QgViewCtrl::raytrace_cmd", 1);
	struct ged *gedp = qgviewctrl_ged(this);
	if (!gedp)
		return;
	const char *av[4] = {nullptr};
	struct bu_vls pid_str = BU_VLS_INIT_ZERO;

	ged_clbk_set(gedp, "ert", BU_CLBK_DURING, &rt_cmd_start, (void *)this);
	ged_clbk_set(gedp, "ert", BU_CLBK_LINGER, &rt_cmd_done, (void *)this);

	if (raytrace_running) {
		if (pid < 0)
			goto cmd_cleanup;
		bu_vls_sprintf(&pid_str, "%d", pid);
		av[0] = "process";
		av[1] = "pabort";
		av[2] = bu_vls_cstr(&pid_str);
		ged_exec_process(gedp, 3, (const char **)av);
		goto cmd_cleanup;
	}

	av[0] = "ert";
	ged_exec_ert(gedp, 1, (const char **)av);
	emit view_changed(QG_VIEW_REFRESH);

cmd_cleanup:
	ged_clbk_set(gedp, "ert", BU_CLBK_DURING, nullptr, nullptr);
	ged_clbk_set(gedp, "ert", BU_CLBK_LINGER, nullptr, nullptr);
	bu_vls_free(&pid_str);
}

void
QgViewCtrl::raytrace_start(int rpid)
{
	QTCAD_SLOT("QgViewCtrl::raytrace_start", 1);
	raytrace->setIcon(QIcon(QPixmap(":images/view/raytrace_abort.png")));
	raytrace_running = true;
	pid = rpid;
}

void
QgViewCtrl::raytrace_done()
{
	QTCAD_SLOT("QgViewCtrl::raytrace_done", 1);
	raytrace->setIcon(QIcon(QPixmap(":images/view/raytrace.png")));
	raytrace_running = false;
	pid = -1;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
