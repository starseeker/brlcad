/*                         Q E L L . C P P
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
/** @file QEll.cpp
 *
 */

#include "common.h"
#include <string.h>
#include <QEvent>
#include <QLabel>
#include <QLineEdit>
#include <QButtonGroup>
#include <QGroupBox>
#include "ged.h"
#include "ged/draw.h"
#include "ged/selection_state.h"
#include "rt/db_io.h"
#include "rt/directory.h"
#include "qtcad/QgGedEventBatch.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgSignalFlags.h"
#include "../qged_edit_preview_util.h"
#include "QEll.h"


QEll::QEll()
    : QWidget()
{
    // TODO - in an ideal world the "default" values would be set
    // and updated in response to view changes (if widget is continually
    // visible) or when it becomes visible...
    ell.magic = RT_ELL_INTERNAL_MAGIC;
    VSET(ell.v, 0, 0, 0);
    VSET(ell.a, 100, 0, 0);
    VSET(ell.b, 0, 200, 0);
    VSET(ell.c, 0, 0, 300);


    QVBoxLayout *l = new QVBoxLayout;

    QLabel *ell_name_label = new QLabel("Object name:");
    l->addWidget(ell_name_label);
    ell_name = new QLineEdit();
    l->addWidget(ell_name);

    QGroupBox *abox = new QGroupBox("Elements");
    QVBoxLayout *abl = new QVBoxLayout;
    abl->setAlignment(Qt::AlignTop);
    O_pnt = new QCheckBox("O:");
    abl->addWidget(O_pnt);
    A_axis = new QCheckBox("A:");
    abl->addWidget(A_axis);
    B_axis = new QCheckBox("B:");
    abl->addWidget(B_axis);
    C_axis = new QCheckBox("C:");
    abl->addWidget(C_axis);
    abox->setLayout(abl);
    l->addWidget(abox);

    QGroupBox *ac_box = new QGroupBox("Actions");
    QVBoxLayout *acl = new QVBoxLayout;
    write_edit = new QPushButton("Apply");
    acl->addWidget(write_edit);
    make_sph = new QPushButton("Make sph");
    acl->addWidget(make_sph);
    reset_values = new QPushButton("Reset");
    acl->addWidget(reset_values);
    ac_box->setLayout(acl);
    l->addWidget(ac_box);


    l->setAlignment(Qt::AlignTop);
    this->setLayout(l);
    setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);


    QObject::connect(ell_name, &QLineEdit::textEdited, this, &QEll::update_viewobj_name);
    QObject::connect(write_edit, &QPushButton::clicked,
	this, &QEll::write_to_db);
    QObject::connect(reset_values, &QPushButton::clicked,
	this, &QEll::read_from_db);
}

QEll::~QEll()
{
    if (!qged_edit_feature_ref_is_null(p))
	qged_edit_preview_publish_event(m_ctx, p, "_ell_edit",
	    QGED_EDIT_PREVIEW_CANCEL, bu_vls_cstr(&oname));
    qged_edit_feature_clear_geometry(p);

    // Empty the realized label node before retiring its record.  Controller
    // replacement can make a retained reference stale, while a still-painted
    // frame may continue to show the old node until the next view operation.
    if (!qged_edit_feature_ref_is_null(labels_p)) {
	qged_edit_feature_labels_replace(labels_p, NULL, 0, NULL);
	qged_edit_feature_set_visible(labels_p, 0);
	qged_edit_feature_remove_ref(labels_p);
    }
    if (m_ctx) {
	qged_edit_feature_remove(m_ctx, "_ell_edit");
	qged_edit_feature_remove(m_ctx, "_ell_edit_labels");
    }
    p = QGED_EDIT_FEATURE_REF_NULL;
    labels_p = QGED_EDIT_FEATURE_REF_NULL;
    bu_vls_free(&oname);
}

void
QEll::clear_labels()
{
    if (!m_ctx) {
	labels_p = QGED_EDIT_FEATURE_REF_NULL;
	return;
    }

    // Resolve against the current active-view controller.  The edit widget
    // can outlive an endpoint/controller replacement, in which case its old
    // value reference correctly becomes stale but the current label feature
    // still needs to be cleared.
    labels_p = qged_edit_feature_label_ensure(m_ctx, "_ell_edit_labels", this);
    if (qged_edit_feature_ref_is_null(labels_p))
	return;
    qged_edit_feature_labels_replace(labels_p, NULL, 0, NULL);
    qged_edit_feature_set_visible(labels_p, 0);
}

struct ged *
QEll::getGed() const
{
    return m_ctx ? m_ctx->getGed() : nullptr;
}

void
QEll::read_from_db()
{
    struct ged *gedp = getGed();
    if (!gedp)
	return;
    struct db_i *dbip = gedp->dbip;
    if (!dbip)
	return;

    if (!dp || dp->d_minor_type != DB5_MINORTYPE_BRLCAD_ELL) {
	return;
    }

    if (!qged_edit_feature_ref_is_null(p))
	qged_edit_preview_publish_event(m_ctx, p, "_ell_edit",
	    QGED_EDIT_PREVIEW_CANCEL, bu_vls_cstr(&oname));
    p = QGED_EDIT_FEATURE_REF_NULL;

    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
	return;
    struct rt_ell_internal *ellp = (struct rt_ell_internal *)intern.idb_ptr;
    RT_ELL_CK_MAGIC(ellp);
    VMOVE(ell.v, ellp->v);
    VMOVE(ell.a, ellp->a);
    VMOVE(ell.b, ellp->b);
    VMOVE(ell.c, ellp->c);
    rt_db_free_internal(&intern);

    // We have pulled new data from disk - let the wireframe know
    update_obj_wireframe();
}

void
QEll::write_to_db()
{
    if (!bu_vls_strlen(&oname))
	return;
    struct ged *gedp = getGed();
    if (!gedp)
	return;
    struct db_i *dbip = gedp->dbip;
    if (!dbip)
	return;

    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_ELL;
    intern.idb_ptr = &ell;
    intern.idb_meth = &OBJ[intern.idb_type];

    dp = db_lookup(dbip, bu_vls_cstr(&oname), LOOKUP_QUIET);
    const int adding = (dp == RT_DIR_NULL) ? 1 : 0;

    {
	QgGedEventBatch event_batch(gedp);

	if (dp == RT_DIR_NULL)
	    dp = db_diradd(dbip, bu_vls_cstr(&oname), RT_DIR_PHONY_ADDR, 0, RT_DIR_SOLID, (void *)&intern.idb_type);

	if (dp == RT_DIR_NULL) {
	    rt_db_free_internal(&intern);
	    return;
	}

	if (rt_db_put_internal(dp, dbip, &intern) < 0) {
	    rt_db_free_internal(&intern);
	    return;
	}
	if (adding)
	    (void)ged_event_notify_object_added(gedp, bu_vls_cstr(&oname),
		NULL);
	else
	    (void)ged_event_notify_object_modified(gedp,
		bu_vls_cstr(&oname), 1, NULL);
    }

    rt_db_free_internal(&intern);

    qged_edit_preview_publish_event(m_ctx, p, "_ell_edit",
	QGED_EDIT_PREVIEW_COMMIT, bu_vls_cstr(&oname));
    p = QGED_EDIT_FEATURE_REF_NULL;

    emit view_updated(QG_VIEW_DB);
}

void
QEll::update_obj_wireframe()
{
    struct ged *gedp = getGed();
    if (!gedp)
	return;

    // Resolve the edit object fresh in case it was removed externally
    // (e.g. by a clear/zap command).
    p = qged_edit_feature_overlay_ensure(m_ctx, "_ell_edit", this, NULL);
    if (qged_edit_feature_ref_is_null(p))
	return;

    // No active db or object name means there is nothing to edit - make sure
    // the edit wireframe is hidden.
    if (!gedp->dbip || !bu_vls_strlen(&oname)) {
	qged_edit_feature_clear_geometry(p);
	qged_edit_feature_set_visible(p, 0);
	clear_labels();
	return;
    }

    // Resolve the complete instance path and accumulate its member matrices.
    // db_lookup() also accepts paths, but returns only the leaf; using it here
    // loses the instance transform and plots labels/edit geometry at the raw
    // primitive coordinates instead of on the selected scene occurrence.
    struct db_full_path full_path;
    db_full_path_init(&full_path);
    mat_t path_mat;
    MAT_IDN(path_mat);
    const int path_ok = db_string_to_path(&full_path, gedp->dbip,
	bu_vls_cstr(&oname)) == 0 &&
	db_path_to_mat(gedp->dbip, &full_path, path_mat,
	    (int)full_path.fp_len - 1);
    dp = path_ok ? DB_FULL_PATH_CUR_DIR(&full_path) : RT_DIR_NULL;
    if (!dp || dp->d_minor_type != DB5_MINORTYPE_BRLCAD_ELL) {
	db_free_full_path(&full_path);
	qged_edit_feature_clear_geometry(p);
	qged_edit_feature_set_visible(p, 0);
	clear_labels();
	return;
    }

    struct rt_ell_internal display_ell = ell;
    MAT4X3PNT(display_ell.v, path_mat, ell.v);
    MAT4X3VEC(display_ell.a, path_mat, ell.a);
    MAT4X3VEC(display_ell.b, path_mat, ell.b);
    MAT4X3VEC(display_ell.c, path_mat, ell.c);
    db_free_full_path(&full_path);

    qged_edit_feature_clear_geometry(p);
    qged_edit_feature_set_view(p, m_ctx);

    // Set up the rt_db_internal and trigger the plotting routine with the
    // current ell parameters
    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_ELL;
    intern.idb_ptr = &display_ell;
    intern.idb_meth = &OBJ[intern.idb_type];
    struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
    if (!wdbp)
	return;
    struct bn_tol *tol = &wdbp->wdb_tol;
    if (qged_edit_feature_replace_ell_wireframe(p,
	    QGED_EDIT_FEATURE_TRANSIENT_PREVIEW,
	    (const struct rt_ell_internal *)intern.idb_ptr))
	qged_edit_preview_publish_event(m_ctx, p, "_ell_edit",
	    QGED_EDIT_PREVIEW_UPDATE, bu_vls_cstr(&oname));

    // At least for now, mimic the MGED behavior and make editing wireframes white
    const char *wcolor = "255/255/255";
    const char *av[2] = {wcolor, NULL};
    struct bu_color cval;
    bu_opt_color(NULL, 1, (const char **)&av[0], (void *)&cval);
    unsigned char rgb[3] = {0, 0, 0};
    bu_color_to_rgb_chars(&cval, rgb);
    qged_edit_feature_set_color(p, rgb[0], rgb[1], rgb[2]);

    // When editing, we show the labels (if any)
    struct rt_point_labels pl[8+1];
    int lcnt = 0;
    mat_t idn_mat;
    MAT_IDN(idn_mat);
    if (intern.idb_meth->ft_labels)
	lcnt = intern.idb_meth->ft_labels(pl, 8, idn_mat, &intern, tol);

    labels_p = qged_edit_feature_label_ensure(m_ctx, "_ell_edit_labels", this);
    if (!qged_edit_feature_ref_is_null(labels_p)) {
	const unsigned char label_color[3] = {255, 255, 0};
	qged_edit_feature_labels_replace(labels_p, pl, lcnt, label_color);
	qged_edit_feature_set_visible(labels_p, lcnt > 0 ? 1 : 0);
    }

    qged_edit_feature_set_visible(p, 1);
    // TODO - we should be able to set UP or DOWN on the various labels
    // when their respective controls are enabled/disabled...

    emit view_updated(QG_VIEW_REFRESH);
}

void
QEll::update_viewobj_name(const QString &)
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip)
	return;

    // Resolve/create the edit view feature.  Don't trust cached pointers here
    // since clear/zap may have removed it.
    p = qged_edit_feature_overlay_ensure(m_ctx, "_ell_edit", this, NULL);
    if (qged_edit_feature_ref_is_null(p))
	return;

    // Make sure the view feature names match whatever the dialog says
    // is the current (proposed) name for the written object
    bu_vls_trunc(&oname, 0);
    if (ell_name->placeholderText().length())
	bu_vls_sprintf(&oname, "%s", ell_name->placeholderText().toLocal8Bit().data());
    if (ell_name->text().length())
	bu_vls_sprintf(&oname, "%s", ell_name->text().toLocal8Bit().data());
    if (!bu_vls_strlen(&oname)) {
	dp = RT_DIR_NULL;
	ged_draw_set_highlighted_shape_ref(gedp, GED_DRAW_SHAPE_REF_NULL);
	qged_edit_feature_clear_geometry(p);
	qged_edit_feature_set_visible(p, 0);
	clear_labels();
	emit view_updated(QG_VIEW_REFRESH);
	return;
    }

    struct db_full_path full_path;
    db_full_path_init(&full_path);
    const int valid_path = db_string_to_path(&full_path, gedp->dbip,
	bu_vls_cstr(&oname)) == 0;
    struct directory *ndp = valid_path ?
	DB_FULL_PATH_CUR_DIR(&full_path) : RT_DIR_NULL;
    db_free_full_path(&full_path);

    const bool object_changed = ndp != dp;
    dp = ndp;
    if (dp && dp->d_minor_type == DB5_MINORTYPE_BRLCAD_ELL) {
	ged_draw_highlight_shape_ref_by_name(gedp, bu_vls_cstr(&oname));
	if (object_changed)
	    read_from_db();
	else
	    update_obj_wireframe();
	return;
    }

    ged_draw_set_highlighted_shape_ref(gedp, GED_DRAW_SHAPE_REF_NULL);
    qged_edit_feature_clear_geometry(p);
    qged_edit_feature_set_visible(p, 0);
    clear_labels();
    emit view_updated(QG_VIEW_REFRESH);
}

void
QEll::sync_selection()
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip)
	return;

    QString selected_path;
    if (ged_selection_count(gedp, NULL) == 1) {
	struct bu_vls paths = BU_VLS_INIT_ZERO;
	ged_selection_list_paths(gedp, NULL, &paths);
	selected_path = QString::fromUtf8(bu_vls_cstr(&paths)).trimmed();
	bu_vls_free(&paths);

	struct db_full_path full_path;
	db_full_path_init(&full_path);
	const bool ell_path = !selected_path.isEmpty() &&
	    db_string_to_path(&full_path, gedp->dbip,
		selected_path.toUtf8().constData()) == 0 &&
	    DB_FULL_PATH_CUR_DIR(&full_path) &&
	    DB_FULL_PATH_CUR_DIR(&full_path)->d_minor_type ==
		DB5_MINORTYPE_BRLCAD_ELL;
	db_free_full_path(&full_path);
	if (!ell_path)
	    selected_path.clear();
    }

    if (selected_path.isEmpty()) {
	// Clear only a value that this synchronization installed.  A manual
	// editor target is allowed to outlive an unrelated selection change.
	if (!selection_path.isEmpty() && ell_name->text() == selection_path) {
	    ell_name->clear();
	    update_viewobj_name(QString());
	}
	selection_path.clear();
	return;
    }

    if (ell_name->text() != selected_path) {
	ell_name->setText(selected_path);
	update_viewobj_name(selected_path);
    }
    selection_path = selected_path;
}

void
QEll::reset_for_database()
{
    /* A path and its directory pointer are meaningful only in the database
     * that resolved them.  Do not leave the old target visible while a new
     * database is active, even if it happens to contain the same leaf name. */
    selection_path.clear();
    ell_name->clear();
    bu_vls_trunc(&oname, 0);
    dp = NULL;

    struct ged *gedp = getGed();
    if (gedp)
	ged_draw_set_highlighted_shape_ref(gedp, GED_DRAW_SHAPE_REF_NULL);
    if (!qged_edit_feature_ref_is_null(p)) {
	qged_edit_feature_clear_geometry(p);
	qged_edit_feature_set_visible(p, 0);
    }
    clear_labels();
    emit view_updated(QG_VIEW_REFRESH);
}

bool
QEll::eventFilter(QObject *, QEvent *e)
{
    (void)e;
    return false;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
