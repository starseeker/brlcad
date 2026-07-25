/*                      Q B R E P . C P P
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
 * You should have received a copy of the GNU Lesser General Public License
 * along with this file; see the file named COPYING for more information.
 */
/** @file QBrep.cpp
 *
 * qged B-rep NURBS CV editing panel.
 */

#include "common.h"

#include <QGroupBox>
#include <QHBoxLayout>
#include <QSignalBlocker>
#include <QVBoxLayout>

#include "brep/edit.h"
#include "ged.h"
#include "qtcad/QgGedEventBatch.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgView.h"
#include "raytrace.h"
#include "rt/geom.h"
#include "rt/rt_ecmds.h"

#include "QBrep.h"

QBrep::QBrep()
    : QWidget()
{
    BN_TOL_INIT(&m_tol);

    QVBoxLayout *layout = new QVBoxLayout;

    QGroupBox *object_box = new QGroupBox("B-rep object");
    QHBoxLayout *object_layout = new QHBoxLayout;
    m_name = new QLineEdit;
    QPushButton *load = new QPushButton("Load");
    object_layout->addWidget(m_name);
    object_layout->addWidget(load);
    object_box->setLayout(object_layout);
    layout->addWidget(object_box);

    QGroupBox *selection_box = new QGroupBox("Control vertex");
    QVBoxLayout *selection_layout = new QVBoxLayout;
    m_selection = new QLabel("No CV selected");
    m_selection->setWordWrap(true);
    selection_layout->addWidget(m_selection);
    QHBoxLayout *weight_layout = new QHBoxLayout;
    m_weight = new QDoubleSpinBox;
    m_weight->setDecimals(8);
    m_weight->setRange(1.0e-9, 1.0e9);
    m_weight->setValue(1.0);
    m_set_weight = new QPushButton("Set weight");
    m_set_weight->setEnabled(false);
    weight_layout->addWidget(m_weight);
    weight_layout->addWidget(m_set_weight);
    selection_layout->addLayout(weight_layout);
    selection_box->setLayout(selection_layout);
    layout->addWidget(selection_box);

    QLabel *policy = new QLabel(
	    "Left-drag a control point in the view plane.  Points whose "
	    "basis support reaches a trim are selectable but locked.");
    policy->setWordWrap(true);
    layout->addWidget(policy);

    QGroupBox *actions_box = new QGroupBox("Transaction");
    QHBoxLayout *actions_layout = new QHBoxLayout;
    m_apply = new QPushButton("Apply");
    m_reset = new QPushButton("Reset");
    m_apply->setEnabled(false);
    m_reset->setEnabled(false);
    actions_layout->addWidget(m_apply);
    actions_layout->addWidget(m_reset);
    actions_box->setLayout(actions_layout);
    layout->addWidget(actions_box);

    m_status = new QLabel;
    m_status->setWordWrap(true);
    layout->addWidget(m_status);
    layout->setAlignment(Qt::AlignTop);
    setLayout(layout);

    m_filter = new QgBrepMoveCVFilter(this);
    connect(m_filter, &QgBrepFilter::brep_selection_changed,
	    this, &QBrep::selection_changed);
    connect(m_filter, &QgBrepFilter::brep_changed,
	    this, &QBrep::edit_changed);
    connect(m_filter, &QgBrepFilter::edit_rejected,
	    this, &QBrep::edit_rejected);
    connect(m_filter, &QgViewFilter::view_updated,
	    this, &QBrep::view_updated);

    connect(load, &QPushButton::clicked, this, &QBrep::load_from_db);
    connect(m_name, &QLineEdit::returnPressed, this, &QBrep::load_from_db);
    connect(m_set_weight, &QPushButton::clicked,
	    this, &QBrep::set_selected_weight);
    connect(m_apply, &QPushButton::clicked, this, &QBrep::apply_to_db);
    connect(m_reset, &QPushButton::clicked, this, &QBrep::reset_from_db);
}

QBrep::~QBrep()
{
    if (m_view)
	m_view->clearFilter(m_filter);
    clear_preview(true);
    clear_edit_state();
}

struct ged *
QBrep::getGed() const
{
    return m_ctx ? m_ctx->getGed() : nullptr;
}

void
QBrep::attachToView(QgView *view)
{
    if (!view || view == m_view)
	return;
    if (m_view)
	m_view->clearFilter(m_filter);
    m_view = view;
    m_view->installFilter(m_filter);
    if (m_es) {
	struct rt_edit_view edit_view;
	rt_edit_view_from_context(&edit_view, m_view->viewContext());
	rt_edit_set_view(m_es, &edit_view);
	update_preview();
    }
}

void
QBrep::detachFromView(QgView *view)
{
    if (!view || view != m_view)
	return;
    view->clearFilter(m_filter);
    m_filter->set_view_widget(nullptr);
    m_view = nullptr;
}

void
QBrep::set_status(const QString &message, bool error)
{
    m_status->setText(message);
    m_status->setStyleSheet(error ? "QLabel { color: #b00020; }" : QString());
    if (error && m_ctx)
	m_ctx->logMessage(message);
}

void
QBrep::clear_edit_state()
{
    m_filter->es = nullptr;
    if (m_es)
	rt_edit_destroy(m_es);
    m_es = nullptr;
    m_dp = nullptr;
    m_baseline_valid = false;
    m_dirty = false;
    m_selection->setText("No CV selected");
    m_set_weight->setEnabled(false);
    m_apply->setEnabled(false);
    m_reset->setEnabled(false);
}

void
QBrep::clear_preview(bool cancel)
{
    if (cancel && !qged_edit_feature_ref_is_null(m_preview) && m_ctx) {
	const QByteArray name = m_name->text().toUtf8();
	qged_edit_preview_publish_event(m_ctx, m_preview, "_brep_edit",
		QGED_EDIT_PREVIEW_CANCEL, name.constData());
    }
    if (!qged_edit_feature_ref_is_null(m_preview))
	qged_edit_feature_remove_ref(m_preview);
    if (!qged_edit_feature_ref_is_null(m_cv_overlay))
	qged_edit_feature_remove_ref(m_cv_overlay);
    m_preview = QGED_EDIT_FEATURE_REF_NULL;
    m_cv_overlay = QGED_EDIT_FEATURE_REF_NULL;
}

void
QBrep::load_from_db()
{
    struct ged *gedp = getGed();
    const QByteArray name = m_name->text().trimmed().toUtf8();
    if (!gedp || !gedp->dbip || name.isEmpty()) {
	set_status("Open a database and enter a B-rep object name.", true);
	return;
    }

    struct directory *dp = db_lookup(gedp->dbip, name.constData(), LOOKUP_QUIET);
    if (!dp || dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BREP) {
	set_status(QString("'%1' is not a B-rep primitive.")
		.arg(QString::fromUtf8(name)), true);
	return;
    }

    clear_preview(true);
    clear_edit_state();

    struct db_full_path path;
    db_full_path_init(&path);
    db_add_node_to_full_path(&path, dp);
    const void *view_context = m_view ? m_view->viewContext()
	: (m_ctx ? m_ctx->activeViewContext() : nullptr);
    m_es = rt_edit_create_context(&path, gedp->dbip, &m_tol, view_context);
    db_free_full_path(&path);
    if (!m_es) {
	set_status("Unable to start a B-rep edit transaction.", true);
	return;
    }

    m_dp = dp;
    m_filter->es = m_es;
    const struct rt_brep_internal *bip =
	(const struct rt_brep_internal *)m_es->es_int.idb_ptr;
    m_baseline_valid = bip && bip->brep && brep_is_valid(bip->brep);
    rt_edit_checkpoint(m_es);
    m_reset->setEnabled(true);
    update_preview();
    set_status(m_baseline_valid
	    ? "Loaded. Click or drag a control point."
	    : "Loaded with pre-existing openNURBS validity errors; only "
	      "trim-independent edits will be offered.");
}

void
QBrep::refresh()
{
    if (!m_es || !m_dp)
	return;
    update_preview();
}

void
QBrep::reset_from_db()
{
    if (!m_dp)
	return;
    load_from_db();
    set_status("Edit transaction reset from the database.");
}

void
QBrep::apply_to_db()
{
    struct ged *gedp = getGed();
    if (!gedp || !gedp->dbip || !m_dp || !m_es)
	return;

    const struct rt_brep_internal *bip =
	(const struct rt_brep_internal *)m_es->es_int.idb_ptr;
    if (!bip || !bip->brep) {
	set_status("The edit transaction has no B-rep data.", true);
	return;
    }
    if (m_baseline_valid && !brep_is_valid(bip->brep)) {
	set_status("Apply refused: the edit would make a valid B-rep invalid.",
		true);
	return;
    }

    {
	QgGedEventBatch event_batch(gedp);
	if (rt_db_put_internal(m_dp, gedp->dbip, &m_es->es_int) < 0) {
	    set_status("Failed to write the edited B-rep.", true);
	    return;
	}
	(void)ged_event_notify_object_modified(gedp, m_dp->d_namep, 1, NULL);
    }

    const QByteArray name = m_name->text().toUtf8();
    qged_edit_preview_publish_event(m_ctx, m_preview, "_brep_edit",
	    QGED_EDIT_PREVIEW_COMMIT, name.constData());
    m_preview = QGED_EDIT_FEATURE_REF_NULL;
    if (!qged_edit_feature_ref_is_null(m_cv_overlay))
	qged_edit_feature_remove_ref(m_cv_overlay);
    m_cv_overlay = QGED_EDIT_FEATURE_REF_NULL;
    m_dirty = false;
    emit view_updated(QG_VIEW_DB);
    load_from_db();
    set_status("B-rep edit committed.");
}

void
QBrep::selection_changed(int face, int cv_u, int cv_v, bool topology_safe)
{
    m_selection->setText(QString("Face %1, CV (%2, %3): %4")
	    .arg(face).arg(cv_u).arg(cv_v)
	    .arg(topology_safe ? "editable" : "locked (influences a trim)"));
    m_set_weight->setEnabled(topology_safe);

    if (!m_es)
	return;
    const struct rt_brep_internal *bip =
	(const struct rt_brep_internal *)m_es->es_int.idb_ptr;
    const ON_NurbsSurface *surface = bip && bip->brep
	? dynamic_cast<const ON_NurbsSurface *>(
	    bip->brep->m_F[face].SurfaceOf()) : nullptr;
    ON_4dPoint cv;
    if (surface && surface->GetCV(
		cv_u, cv_v, ON::euclidean_rational, &cv.x)) {
	QSignalBlocker blocker(m_weight);
	m_weight->setValue(cv.w);
    }
    set_status(topology_safe
	    ? "Drag in the view plane or enter a rational weight."
	    : "This CV is locked because moving it can change a trim locus.");
}

void
QBrep::set_selected_weight()
{
    if (!m_es || !m_filter->selected_cv_is_topology_safe())
	return;
    EDOBJ[m_es->es_int.idb_type].ft_set_edit_mode(
	    m_es, ECMD_BREP_SRF_CV_WEIGHT);
    m_es->e_inpara = 1;
    m_es->e_para[0] = m_weight->value();
    rt_edit_process(m_es);
    edit_changed();
}

void
QBrep::edit_changed()
{
    m_dirty = true;
    m_apply->setEnabled(true);
    update_preview();
    set_status("Modified in memory. Apply commits; Reset discards.");
    emit view_updated(QG_VIEW_REFRESH);
}

void
QBrep::edit_rejected(const QString &reason)
{
    set_status(reason, true);
}

void
QBrep::update_preview()
{
    struct ged *gedp = getGed();
    if (!m_ctx || !gedp || !m_es || !m_dp)
	return;

    m_preview = qged_edit_feature_overlay_ensure(
	    m_ctx, "_brep_edit", this, m_dp->d_namep);
    if (qged_edit_feature_ref_is_null(m_preview))
	return;

    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_TOL;
    qged_edit_feature_set_view(m_preview, m_ctx);
    qged_edit_feature_set_visible(m_preview, 1);
    qged_edit_feature_set_color(m_preview, 255, 255, 255);
    if (qged_edit_feature_replace_primitive_wireframe(m_preview,
		gedp->dbip, &m_es->es_int, m_es->e_mat, &ttol, &m_tol)) {
	qged_edit_feature_touch(m_preview);
	qged_edit_preview_publish_event(m_ctx, m_preview, "_brep_edit",
		QGED_EDIT_PREVIEW_UPDATE, m_dp->d_namep);
    }
    update_cv_overlay();
}

void
QBrep::update_cv_overlay()
{
    if (!m_ctx || !m_es || !m_dp)
	return;
    const struct rt_brep_internal *bip =
	(const struct rt_brep_internal *)m_es->es_int.idb_ptr;
    if (!bip || !bip->brep)
	return;

    m_cv_overlay = qged_edit_feature_overlay_ensure(
	    m_ctx, "_brep_edit_cvs", m_filter, m_dp->d_namep);
    if (qged_edit_feature_ref_is_null(m_cv_overlay))
	return;

    struct qged_edit_preview_lines lines;
    qged_edit_preview_lines_init(&lines);
    for (int fi = 0; fi < bip->brep->m_F.Count(); ++fi) {
	const ON_NurbsSurface *surface =
	    dynamic_cast<const ON_NurbsSurface *>(
		bip->brep->m_F[fi].SurfaceOf());
	if (!surface)
	    continue;

	for (int i = 0; i < surface->CVCount(0); ++i) {
	    for (int j = 0; j < surface->CVCount(1); ++j) {
		ON_3dPoint cv;
		if (!surface->GetCV(i, j, cv))
		    continue;
		point_t local = VINIT_ZERO;
		point_t model = VINIT_ZERO;
		VSET(local, cv.x, cv.y, cv.z);
		MAT4X3PNT(model, m_es->e_mat, local);
		qged_edit_preview_lines_append(&lines, model,
			j ? QGED_EDIT_PREVIEW_LINE_DRAW
			  : QGED_EDIT_PREVIEW_LINE_MOVE);
	    }
	}
	for (int j = 0; j < surface->CVCount(1); ++j) {
	    for (int i = 0; i < surface->CVCount(0); ++i) {
		ON_3dPoint cv;
		if (!surface->GetCV(i, j, cv))
		    continue;
		point_t local = VINIT_ZERO;
		point_t model = VINIT_ZERO;
		VSET(local, cv.x, cv.y, cv.z);
		MAT4X3PNT(model, m_es->e_mat, local);
		qged_edit_preview_lines_append(&lines, model,
			i ? QGED_EDIT_PREVIEW_LINE_DRAW
			  : QGED_EDIT_PREVIEW_LINE_MOVE);
		qged_edit_preview_lines_append(&lines, model,
			QGED_EDIT_PREVIEW_POINT_DRAW);
	    }
	}
    }

    qged_edit_feature_set_view(m_cv_overlay, m_ctx);
    qged_edit_feature_set_visible(m_cv_overlay, 1);
    qged_edit_feature_set_color(m_cv_overlay, 255, 176, 0);
    qged_edit_preview_lines_replace(m_cv_overlay,
	    QGED_EDIT_FEATURE_TRANSIENT_PREVIEW, &lines);
    qged_edit_preview_lines_free(&lines);
    qged_edit_feature_touch(m_cv_overlay);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
