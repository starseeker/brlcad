/*                       Q B O T . C P P
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
/** @file QBot.cpp
 *
 * Edit-preview BOT editor.
 *
 * The transient preview feature owns typed edit-preview callbacks and routes
 * lifecycle events through the shared qged preview helper.
 */

#include "common.h"
#include <QEvent>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QGroupBox>
#include <QVBoxLayout>
#include "ged.h"
#include "rt/db_io.h"
#include "rt/directory.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgSignalFlags.h"
#include "../qged_edit_preview_util.h"
#include "QBot.h"


/* ---- edit-preview callbacks -------------------------------------------- */

static uint64_t
_bot_preview_revision(void *UNUSED(preview_ctx))
{
    /* Placeholder — a real implementation would track the edit epoch. */
    return 0;
}

static int
_bot_preview_update(void *preview_ctx)
{
    QBot *self = (QBot *)preview_ctx;
    if (!self)
	return 0;
    QMetaObject::invokeMethod(self, "update_obj_wireframe", Qt::DirectConnection);
    return 1;
}

/* pick_cb stub: ray-pick for face / vertex selection.
 * Filled in once the BOT editor gains interactive selection logic. */
static int
_bot_preview_pick(void *UNUSED(preview_ctx), int UNUSED(x), int UNUSED(y),
		  void *UNUSED(pick_out))
{
    /* TODO: intersect sample ray against BOT faces/vertices. */
    return 0;
}

/* ---- QBot constructor --------------------------------------------------- */

QBot::QBot()
    : QWidget()
{
    QVBoxLayout *l = new QVBoxLayout;

    QLabel *name_label = new QLabel("Object name:");
    l->addWidget(name_label);
    bot_name = new QLineEdit();
    l->addWidget(bot_name);

    QGroupBox *mode_box = new QGroupBox("Edit mode");
    QVBoxLayout *mbl = new QVBoxLayout;
    edit_mode = new QComboBox();
    edit_mode->addItem("Vertex");
    edit_mode->addItem("Face");
    edit_mode->addItem("Edge");
    mbl->addWidget(edit_mode);
    mode_box->setLayout(mbl);
    l->addWidget(mode_box);

    QGroupBox *ac_box = new QGroupBox("Actions");
    QVBoxLayout *acl = new QVBoxLayout;
    write_edit = new QPushButton("Apply");
    acl->addWidget(write_edit);
    reset_values = new QPushButton("Reset");
    acl->addWidget(reset_values);
    ac_box->setLayout(acl);
    l->addWidget(ac_box);

    l->setAlignment(Qt::AlignTop);
    this->setLayout(l);

    QObject::connect(bot_name, &QLineEdit::textChanged,
		     this, &QBot::update_viewobj_name);
    QObject::connect(write_edit, &QPushButton::clicked,
		     this, &QBot::write_to_db);
    QObject::connect(reset_values, &QPushButton::clicked,
		     this, &QBot::read_from_db);
}

QBot::~QBot()
{
    qged_edit_feature_clear_geometry(p);
    if (!qged_edit_feature_ref_is_null(p) && m_ctx) {
	qged_edit_feature_remove(m_ctx, "_bot_edit");
	p = QGED_EDIT_FEATURE_REF_NULL;
    }
    bu_vls_free(&oname);
}

struct ged *
QBot::getGed() const
{
    if (!m_ctx)
	return nullptr;
    return m_ctx->getGed();
}

void
QBot::read_from_db()
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip || !bu_vls_strlen(&oname))
	return;

    /* BOT loading is deferred until full edit logic is implemented. */
    update_obj_wireframe();
    emit view_updated(QG_VIEW_REFRESH);
}

void
QBot::write_to_db()
{
    qged_edit_preview_publish_event(m_ctx, p, QGED_EDIT_PREVIEW_COMMIT,
	    bu_vls_cstr(&oname));
    emit view_updated(QG_VIEW_DB);
}

void
QBot::update_obj_wireframe()
{
    struct ged *gedp = getGed();
    if (!gedp)
	return;

    struct qged_edit_preview_callbacks callbacks = QGED_EDIT_PREVIEW_CALLBACKS_INIT;
    callbacks.revision_cb = _bot_preview_revision;
    callbacks.update_cb = _bot_preview_update;
    callbacks.pick_cb = _bot_preview_pick;
    p = qged_edit_feature_overlay_ensure(m_ctx, "_bot_edit", this, this,
	    &callbacks, bu_vls_cstr(&oname));
    if (qged_edit_feature_ref_is_null(p))
	return;

    if (!gedp->dbip || !bu_vls_strlen(&oname)) {
	qged_edit_feature_clear_geometry(p);
	qged_edit_feature_set_visible(p, 0);
	return;
    }

    dp = db_lookup(gedp->dbip, bu_vls_cstr(&oname), LOOKUP_QUIET);
    if (!dp || dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	qged_edit_feature_clear_geometry(p);
	qged_edit_feature_set_visible(p, 0);
	return;
    }

    qged_edit_feature_clear_geometry(p);
    qged_edit_feature_set_view(p, m_ctx);
    qged_edit_feature_set_visible(p, 1);

    /* Load and plot the BOT geometry. */
    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    if (rt_db_get_internal(&intern, dp, gedp->dbip, NULL) < 0)
	return;
    if (intern.idb_minor_type != DB5_MINORTYPE_BRLCAD_BOT) {
	rt_db_free_internal(&intern);
	return;
    }

    const struct rt_bot_internal *bot_ip =
	(const struct rt_bot_internal *)intern.idb_ptr;
    if (qged_edit_feature_replace_bot_face_lines(p,
	    QGED_EDIT_FEATURE_TRANSIENT_PREVIEW, bot_ip) &&
	    !qged_edit_feature_ref_is_null(p)) {
	qged_edit_preview_publish_event(m_ctx, p, QGED_EDIT_PREVIEW_UPDATE,
		bu_vls_cstr(&oname));
    }
    rt_db_free_internal(&intern);

    const char *wcolor = "255/255/255";
    const char *av[2] = {wcolor, NULL};
    struct bu_color cval;
    bu_opt_color(NULL, 1, (const char **)&av[0], (void *)&cval);
    unsigned char rgb[3] = {0, 0, 0};
    bu_color_to_rgb_chars(&cval, rgb);
    qged_edit_feature_set_color(p, rgb[0], rgb[1], rgb[2]);

    qged_edit_feature_touch(p);
}

void
QBot::update_viewobj_name(const QString &ostr)
{
    bu_vls_sprintf(&oname, "%s", ostr.toLocal8Bit().data());
    update_obj_wireframe();
    emit view_updated(QG_VIEW_REFRESH);
}

bool
QBot::eventFilter(QObject *, QEvent *)
{
    return false;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
