/*                   Q R E V O L V E . C P P
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
/** @file QRevolve.cpp
 *
 * Edit-preview revolve-solid editor.
 *
 * The transient preview feature owns typed edit-preview callbacks and routes
 * begin/update/commit events through the shared qged preview helper.
 */

#include "common.h"
#include <QEvent>
#include <QLabel>
#include <QLineEdit>
#include <QGroupBox>
#include <QVBoxLayout>
#include "ged.h"
#include "rt/db_io.h"
#include "rt/directory.h"
#include "qtcad/QgGedEventBatch.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgSignalFlags.h"
#include "../qged_edit_preview_util.h"
#include "QRevolve.h"


/* ---- edit-preview callbacks -------------------------------------------- */

static uint64_t
_revolve_preview_revision(void *UNUSED(preview_ctx))
{
    return 0;
}

static int
_revolve_preview_update(void *preview_ctx)
{
    QRevolve *self = (QRevolve *)preview_ctx;
    if (!self)
	return 0;
    QMetaObject::invokeMethod(self, "update_obj_wireframe", Qt::DirectConnection);
    return 1;
}

/* ---- QRevolve constructor ----------------------------------------------- */

QRevolve::QRevolve()
    : QWidget()
{
    rev.magic = RT_REVOLVE_INTERNAL_MAGIC;

    QVBoxLayout *l = new QVBoxLayout;

    QLabel *name_label = new QLabel("Object name:");
    l->addWidget(name_label);
    revolve_name = new QLineEdit();
    l->addWidget(revolve_name);

    QGroupBox *pbox = new QGroupBox("Parameters");
    QVBoxLayout *pbl = new QVBoxLayout;
    pbl->setAlignment(Qt::AlignTop);
    QLabel *a_label = new QLabel("Angle (deg):");
    pbl->addWidget(a_label);
    angle = new QLineEdit();
    pbl->addWidget(angle);
    pbox->setLayout(pbl);
    l->addWidget(pbox);

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

    QObject::connect(revolve_name, &QLineEdit::textChanged,
		     this, &QRevolve::update_viewobj_name);
    QObject::connect(write_edit, &QPushButton::clicked,
		     this, &QRevolve::write_to_db);
    QObject::connect(reset_values, &QPushButton::clicked,
		     this, &QRevolve::read_from_db);
}

QRevolve::~QRevolve()
{
    qged_edit_feature_clear_geometry(p);
    if (!qged_edit_feature_ref_is_null(p) && m_ctx) {
	qged_edit_feature_remove(m_ctx, "_revolve_edit");
	p = QGED_EDIT_FEATURE_REF_NULL;
    }
    bu_vls_free(&oname);
}

struct ged *
QRevolve::getGed() const
{
    if (!m_ctx)
	return nullptr;
    return m_ctx->getGed();
}

void
QRevolve::read_from_db()
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip || !bu_vls_strlen(&oname))
	return;

    struct directory *ldp = db_lookup(gedp->dbip, bu_vls_cstr(&oname), LOOKUP_QUIET);
    if (!ldp || ldp->d_minor_type != DB5_MINORTYPE_BRLCAD_REVOLVE)
	return;

    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    if (rt_db_get_internal(&intern, ldp, gedp->dbip, NULL) < 0)
	return;

    struct rt_revolve_internal *rp = (struct rt_revolve_internal *)intern.idb_ptr;
    RT_REVOLVE_CK_MAGIC(rp);

    rev = *rp;

    rt_db_free_internal(&intern);
    update_obj_wireframe();
    emit view_updated(QG_VIEW_REFRESH);
}

void
QRevolve::write_to_db()
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip || !bu_vls_strlen(&oname))
	return;

    struct directory *ldp = db_lookup(gedp->dbip, bu_vls_cstr(&oname), LOOKUP_QUIET);
    if (!ldp || ldp->d_minor_type != DB5_MINORTYPE_BRLCAD_REVOLVE)
	return;

    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_REVOLVE;
    intern.idb_ptr = &rev;
    intern.idb_meth = &OBJ[ID_REVOLVE];

    {
	QgGedEventBatch event_batch(gedp);
	if (rt_db_put_internal(ldp, gedp->dbip, &intern) < 0)
	    return;
    }

    qged_edit_preview_publish_event(m_ctx, p, QGED_EDIT_PREVIEW_COMMIT,
	    bu_vls_cstr(&oname));
    emit view_updated(QG_VIEW_DB);
}

void
QRevolve::update_obj_wireframe()
{
    struct ged *gedp = getGed();
    if (!gedp)
	return;

    struct qged_edit_preview_callbacks callbacks = QGED_EDIT_PREVIEW_CALLBACKS_INIT;
    callbacks.revision_cb = _revolve_preview_revision;
    callbacks.update_cb = _revolve_preview_update;
    p = qged_edit_feature_overlay_ensure(m_ctx, "_revolve_edit", this, this,
	    &callbacks, bu_vls_cstr(&oname));
    if (qged_edit_feature_ref_is_null(p))
	return;

    if (!gedp->dbip || !bu_vls_strlen(&oname)) {
	qged_edit_feature_clear_geometry(p);
	qged_edit_feature_set_visible(p, 0);
	return;
    }

    dp = db_lookup(gedp->dbip, bu_vls_cstr(&oname), LOOKUP_QUIET);
    if (!dp || dp->d_minor_type != DB5_MINORTYPE_BRLCAD_REVOLVE) {
	qged_edit_feature_clear_geometry(p);
	qged_edit_feature_set_visible(p, 0);
	return;
    }

    qged_edit_feature_clear_geometry(p);
    qged_edit_feature_set_view(p, m_ctx);
    qged_edit_feature_set_visible(p, 1);

    struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
    if (!wdbp)
	return;
    if (qged_edit_feature_replace_revolve_wireframe(p,
	    QGED_EDIT_FEATURE_TRANSIENT_PREVIEW, &rev,
	    &wdbp->wdb_ttol) && !qged_edit_feature_ref_is_null(p))
	qged_edit_preview_publish_event(m_ctx, p, QGED_EDIT_PREVIEW_UPDATE,
		bu_vls_cstr(&oname));

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
QRevolve::update_viewobj_name(const QString &ostr)
{
    bu_vls_sprintf(&oname, "%s", ostr.toLocal8Bit().data());
    update_obj_wireframe();
    emit view_updated(QG_VIEW_REFRESH);
}

bool
QRevolve::eventFilter(QObject *, QEvent *)
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
