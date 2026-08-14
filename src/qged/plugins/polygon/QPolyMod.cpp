/*                     Q P O L Y M O D . C P P
 * BRL-CAD
 *
 * Copyright (c) 2014-2026 United States Government as represented by
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
/** @file QPolyMod.cpp
 *
 */

#include "common.h"

#include "bv.h"
#include <vector>
#include <QLabel>
#include <QLineEdit>
#include <QButtonGroup>
#include <QGroupBox>
#include <QtGlobal>
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgGedEventBatch.h"
#include "qtcad/QgPolyFilter.h"
#include "qtcad/QgView.h"
#include "ged.h"
#include "ged/draw.h"
#include "rt/directory.h"
#include "rt/db_io.h"
#include "QPolyCreate.h"
#include "QPolyMod.h"

/* Collect polygon objects (optionally excluding one). */
struct _qpolymod_poly_collect {
    std::vector<qg_polygon_ref> *polys;
    qg_polygon_ref exclude;
};
extern "C" int
_qpolymod_poly_collect_cb(qg_polygon_ref ref, const qg_polygon_record *, void *data)
{
    struct _qpolymod_poly_collect *s = (struct _qpolymod_poly_collect *)data;
    if (ref.owner != s->exclude.owner || ref.id != s->exclude.id ||
	ref.generation != s->exclude.generation)
	s->polys->push_back(ref);
    return 1;
}

static void *
qpolymod_view(const QgPluginContext *ctx)
{
    return ctx ? ctx->activeViewContext() : nullptr;
}

static struct ged_view_context *
qpolymod_ged_view(void *view)
{
    return ged_view_context_from_bv(static_cast<struct bv_context *>(view));
}

static struct ged_view_context *
qpolymod_current_ged_view(const QgPluginContext *ctx)
{
    return qpolymod_ged_view(qpolymod_view(ctx));
}

static QgView *
qpolymod_view_widget(const QgPluginContext *ctx)
{
    return ctx ? ctx->getViewWidget() : nullptr;
}

QPolyMod::QPolyMod()
    : QWidget()
{
    const QString testPrefix =
	QStringLiteral("org.brlcad.qged.polygon.modify");
    QVBoxLayout *l = new QVBoxLayout;

    QButtonGroup *t_grp = new QButtonGroup();

    select_mode = new QRadioButton("Select");
    select_mode->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".select"));
    t_grp->addButton(select_mode);
    l->addWidget(select_mode);
    QLabel *cs_label = new QLabel("Currently selected polygon:");
    l->addWidget(cs_label);
    mod_names = new QComboBox(this);
    mod_names->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".polygon-name"));
    QObject::connect(mod_names, &QComboBox::currentTextChanged, this, &QPolyMod::select);
    l->addWidget(mod_names);


    QGroupBox *defaultBox = new QGroupBox("Settings");
    QVBoxLayout *default_gl = new QVBoxLayout;
    default_gl->setAlignment(Qt::AlignTop);
    ps = new QPolySettings();
    ps->setTestIdPrefix(testPrefix);
    default_gl->addWidget(ps);
    defaultBox->setLayout(default_gl);
    l->addWidget(defaultBox);
    QObject::connect(ps, &QPolySettings::settings_changed, this, &QPolyMod::polygon_update_props);
    // We'll need to be aware if we specify a colliding polygon object name
    QObject::connect(ps->view_name, &QLineEdit::textEdited, this, &QPolyMod::view_name_edit_str);
    QObject::connect(ps->view_name, &QLineEdit::editingFinished, this, &QPolyMod::view_name_update);
    // The sketch name gets enabled/disabled and in some modes syncs with the
    // view name.
    QObject::connect(ps->sketch_sync, &QCheckBox::toggled, this, &QPolyMod::sketch_sync_bool);
    QObject::connect(ps->sketch_name, &QLineEdit::textEdited, this, &QPolyMod::sketch_name_edit_str);
    QObject::connect(ps->view_name, &QLineEdit::editingFinished, this, &QPolyMod::sketch_name_edit);
    QObject::connect(ps->sketch_name, &QLineEdit::editingFinished, this, &QPolyMod::sketch_name_update);
    QObject::connect(ps, &QPolySettings::grid_snapping_changed, this, &QPolyMod::toggle_grid_snapping);
    QObject::connect(ps, &QPolySettings::line_snapping_changed, this, &QPolyMod::toggle_line_snapping);

    QGroupBox *modpolyBox = new QGroupBox("Modify Polygon");
    QVBoxLayout *mod_poly_gl = new QVBoxLayout;

    move_mode = new QRadioButton("Move");
    move_mode->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".move"));
    QObject::connect(move_mode, &QRadioButton::toggled, this, &QPolyMod::toplevel_config);
    t_grp->addButton(move_mode);
    update_mode = new QRadioButton("Update geometry");
    update_mode->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".update"));
    t_grp->addButton(update_mode);

    close_general_poly = new QCheckBox("Close polygon");
    close_general_poly->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".close"));
    // Disabled if we're not a general polygon
    close_general_poly->setChecked(true);
    close_general_poly->setDisabled(true);
    // Append polygon pnt if above is true (switch to select).
    QObject::connect(close_general_poly, &QCheckBox::toggled, this, &QPolyMod::toggle_closed_poly);

    // Basic shape updating is simple from an interaction standpoint, but the
    // general case is a bit more involved.
    general_mode_opts = new QGroupBox("General Polygon Modes");
    QVBoxLayout *go_l = new QVBoxLayout();
    go_l->setAlignment(Qt::AlignTop);
    //go_l->setSpacing(0);
    //go_l->setContentsMargins(1,1,1,1);

    QButtonGroup *gm_box = new QButtonGroup();
    append_pnt = new QRadioButton("Append polygon pnt");
    append_pnt->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".append-point"));
    append_pnt->setChecked(true);
    gm_box->addButton(append_pnt);
    go_l->addWidget(append_pnt);
    select_pnt = new QRadioButton("Select polygon pnt");
    select_pnt->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".select-point"));
    QObject::connect(select_pnt, &QRadioButton::toggled, this, &QPolyMod::clear_pnt_selection);
    gm_box->addButton(select_pnt);
    go_l->addWidget(select_pnt);
    general_mode_opts->setLayout(go_l);
    QObject::connect(update_mode, &QRadioButton::toggled, this, &QPolyMod::toplevel_config);
    general_mode_opts->setDisabled(true);

    mod_poly_gl->addWidget(move_mode);
    mod_poly_gl->addWidget(update_mode);
    mod_poly_gl->addWidget(close_general_poly);
    mod_poly_gl->addWidget(general_mode_opts);

    modpolyBox->setLayout(mod_poly_gl);
    l->addWidget(modpolyBox);


    QGroupBox *boolBox = new QGroupBox("Apply Boolean Op");
    QVBoxLayout *bool_gl = new QVBoxLayout;
    bool_gl->setAlignment(Qt::AlignTop);
    QLabel *csg_modes_label = new QLabel("Current Operation:");
    bool_gl->addWidget(csg_modes_label);
    csg_modes = new QComboBox();
    csg_modes->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".boolean-mode"));
    csg_modes->addItem("Union");
    csg_modes->addItem("Subtraction");
    csg_modes->addItem("Intersection");
    csg_modes->setCurrentIndex(0);
    bool_gl->addWidget(csg_modes);
    apply_bool = new QPushButton("Apply");
    apply_bool->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".boolean-apply"));
    bool_gl->addWidget(apply_bool);
    QObject::connect(apply_bool, &QPushButton::released, this, &QPolyMod::apply_bool_op);
    boolBox->setLayout(bool_gl);
    l->addWidget(boolBox);


    QGroupBox *viewsnapBox = new QGroupBox("Align View to Polygon");
    QVBoxLayout *viewsnap_gl = new QVBoxLayout;
    viewsnap_gl->setAlignment(Qt::AlignTop);
    viewsnap_poly = new QPushButton("Align");
    viewsnap_poly->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".align"));
    viewsnap_gl->addWidget(viewsnap_poly);
    QObject::connect(viewsnap_poly, &QPushButton::released, this, &QPolyMod::align_to_poly);
    viewsnapBox->setLayout(viewsnap_gl);
    l->addWidget(viewsnapBox);


    QGroupBox *removeBox = new QGroupBox("Remove Polygon");
    QVBoxLayout *remove_gl = new QVBoxLayout;
    remove_gl->setAlignment(Qt::AlignTop);
    remove_poly = new QPushButton("Delete");
    remove_poly->setProperty(
	"qgTestId", testPrefix + QStringLiteral(".delete"));
    remove_poly->setStyleSheet("QPushButton {color: red;}");
    remove_gl->addWidget(remove_poly);
    QObject::connect(remove_poly, &QPushButton::released, this, &QPolyMod::delete_poly);
    removeBox->setLayout(remove_gl);
    l->addWidget(removeBox);


    l->setAlignment(Qt::AlignTop);
    this->setLayout(l);
    setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);

    select_mode->setChecked(true);
    mod_names_reset();
    toplevel_config(true);

    puf = new QgPolyUpdateFilter();
    psf = new QgPolySelectFilter();
    ppf = new QgPolyPointFilter();
    pmf = new QgPolyMoveFilter();
}

QPolyMod::~QPolyMod()
{
}

struct ged *
QPolyMod::getGed() const
{
    return m_ctx ? m_ctx->getGed() : nullptr;
}

void
QPolyMod::app_mod_names_reset(void *)
{
    mod_names_reset();
}

void
QPolyMod::mod_names_reset()
{
    struct ged *gedp = getGed();
    if (!gedp)
	return;
    void *v = qpolymod_view(m_ctx);
    if (!v)
	return;

    /* Geometry-only updates must not call this routine: rebuilding the combo
     * box during every drag perturbs editor selection and adds work to the
     * latency-sensitive mouse path.  Membership/name changes and tool
     * activation call it explicitly.  Preserve the selected scene handle
     * when it still exists; otherwise select the first live record. */
    QString selected_name;
    if (!ged_view_polygon_ref_is_null(p)) {
	qg_polygon_record selected_rec;
	if (ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), p,
		&selected_rec) && selected_rec.name)
	    selected_name = QString::fromLocal8Bit(selected_rec.name);
    }

    mod_names->blockSignals(true);
    mod_names->clear();
    if (gedp) {
	std::vector<qg_polygon_ref> polyvec;
	struct _qpolymod_poly_collect pc;
	pc.polys = &polyvec;
	pc.exclude = GED_VIEW_POLYGON_REF_NULL;
	ged_view_polygon_visit_records(qpolymod_ged_view(v),
		_qpolymod_poly_collect_cb, &pc);
	for (auto s : polyvec) {
	    qg_polygon_record rec;
	    if (ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), s, &rec) && rec.name)
		mod_names->addItem(rec.name);
	}
    }
    int selected_index = selected_name.isEmpty() ? -1 :
	mod_names->findText(selected_name);
    if (selected_index < 0 && mod_names->count() > 0)
	selected_index = 0;
    mod_names->setCurrentIndex(selected_index);
    if (selected_index >= 0) {
	select(mod_names->currentText());
    } else {
	p = GED_VIEW_POLYGON_REF_NULL;
	ps->view_name->clear();
    }
    mod_names->blockSignals(false);
}

void
QPolyMod::poly_type_settings(const struct ged_view_polygon_record *ip)
{
    if (!ip || !ip->contour_count)
	return;
    if (ip->type == QG_POLYGON_GENERAL) {
	general_mode_opts->setEnabled(true);
	close_general_poly->setEnabled(true);
	close_general_poly->blockSignals(true);
	append_pnt->blockSignals(true);
	select_pnt->blockSignals(true);
	if (!ip->first_contour_open) {
	    close_general_poly->setChecked(true);
	    select_pnt->setChecked(true);
	    append_pnt->setChecked(false);
	    append_pnt->setEnabled(false);
	} else {
	    close_general_poly->setChecked(false);
	    append_pnt->setEnabled(true);
	    append_pnt->setChecked(true);
	    select_pnt->setChecked(false);
	}
	close_general_poly->blockSignals(false);
	append_pnt->blockSignals(false);
	select_pnt->blockSignals(false);
    } else {
	close_general_poly->blockSignals(true);
	close_general_poly->setChecked(true);
	close_general_poly->setEnabled(false);
	close_general_poly->blockSignals(false);
	general_mode_opts->setEnabled(false);
    }
}

void
QPolyMod::polygon_update_props()
{
    struct ged *gedp = getGed();
    if (!gedp || ged_view_polygon_ref_is_null(p))
	return;

    ged_view_polygon_set_visual(qpolymod_current_ged_view(m_ctx), p, &ps->edge_color->bc, &ps->fill_color->bc,
	    (fastf_t)(ps->fill_slope_x->text().toDouble()),
	    (fastf_t)(ps->fill_slope_y->text().toDouble()),
	    (fastf_t)(ps->fill_density->text().toDouble()),
	    (fastf_t)(ps->vZ->text().toDouble()),
	    ps->fill_poly->isChecked() ? 1 : 0);
    emit view_updated(QG_VIEW_REFRESH);
}

void
QPolyMod::toplevel_config(bool)
{
    // Initialize
    struct ged *gedp = getGed();
    void *v = qpolymod_view(m_ctx);
    if (!gedp || !v)
	return;
    bool draw_change = false;

    // This function is called when a top level mode change was initiated
    // by a selection button.  Clear any selected points being displayed -
    // when we're switching modes at this level, we always start with a
    // blank slate for points.
    if (gedp)
	draw_change = ged_view_polygon_clear_point_selection(
		qpolymod_ged_view(v)) ? true : false;

    // Make sure the Combo box list is current.
    mod_names_reset();

    // If we have a current p, we know the current type - enable/disable
    // the general settings on that basis

    if (!ged_view_polygon_ref_is_null(p)) {
	qg_polygon_record rec;
	if (ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), p, &rec))
	    poly_type_settings(&rec);
    }

    if (draw_change && gedp)
	emit view_updated(QG_VIEW_REFRESH);
}


void
QPolyMod::clear_pnt_selection(bool checked)
{
    if (checked)
	return;
    int ptype = -1;
    qg_polygon_record rec;
    if (!ged_view_polygon_ref_is_null(p)) {
	if (!ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), p, &rec))
	    return;
	ptype = rec.type;
    }
    if (ptype != QG_POLYGON_GENERAL) {
	return;
    }
    bu_log("got pnt selection clear\n");
    ged_view_polygon_clear_selected_point(qpolymod_current_ged_view(m_ctx), p);

    struct ged *gedp = getGed();
    if (!gedp)
	return;

    emit view_updated(QG_VIEW_REFRESH);
}

void
QPolyMod::select(const QString &poly)
{
    struct ged *gedp = getGed();
    if (!gedp)
	return;
    void *v = qpolymod_view(m_ctx);
    if (!v)
	return;

    p = GED_VIEW_POLYGON_REF_NULL;
    p = ged_view_polygon_find(qpolymod_ged_view(v),
	    poly.toLocal8Bit().data());
    if (!ged_view_polygon_ref_is_null(p)) {
	qg_polygon_record rec;
	if (!ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), p, &rec))
	    return;
	poly_type_settings(&rec);
	ps->settings_sync(&rec);
	ps->view_name->setText(poly);
	if (rec.user_data) {
	    struct directory *dp = (struct directory *)rec.user_data;
	    ps->sketch_sync->blockSignals(true);
	    ps->sketch_sync->setChecked(true);
	    ps->sketch_sync->blockSignals(false);
	    ps->sketch_name->blockSignals(true);
	    ps->sketch_name->setText(dp->d_namep);
	    ps->sketch_name->setEnabled(true);
	    ps->sketch_name->blockSignals(false);
	}

	return;
    }
}

void
QPolyMod::toggle_closed_poly(bool checked)
{
    int ptype = -1;
    qg_polygon_record rec;
    if (!ged_view_polygon_ref_is_null(p)) {
	if (ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), p, &rec))
	    ptype = rec.type;
    }

    if (ptype != QG_POLYGON_GENERAL) {
	if (!checked) {
	    close_general_poly->blockSignals(true);
	    close_general_poly->setChecked(true);
	    close_general_poly->setEnabled(false);
	    close_general_poly->blockSignals(false);
	    general_mode_opts->setEnabled(false);
	}
	return;
    }

    clear_pnt_selection(false);

    // If we don't have a valid polygon, bail
    if (!rec.contour_count)
	return;

    if (checked && ptype == QG_POLYGON_GENERAL) {
	// A contour with less than 3 points can't be closed
	if (rec.point_count < 3)
	    ged_view_polygon_set_all_contours_open(qpolymod_current_ged_view(m_ctx), p, 1);
	else
	    ged_view_polygon_set_all_contours_open(qpolymod_current_ged_view(m_ctx), p, 0);
    } else {
	ged_view_polygon_set_all_contours_open(qpolymod_current_ged_view(m_ctx), p, 1);
    }

    ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), p, &rec);

    close_general_poly->blockSignals(true);
    append_pnt->blockSignals(true);
    select_pnt->blockSignals(true);

    if (!rec.first_contour_open) {
	close_general_poly->setChecked(true);
	select_pnt->setChecked(true);
	append_pnt->setChecked(false);
	append_pnt->setEnabled(false);
    } else {
	close_general_poly->setChecked(false);
	append_pnt->setEnabled(true);
	append_pnt->setChecked(true);
	select_pnt->setChecked(false);
    }

    close_general_poly->blockSignals(false);
    append_pnt->blockSignals(false);
    select_pnt->blockSignals(false);

    struct ged *gedp = getGed();
    if (!gedp)
	return;
    void *v = qpolymod_view(m_ctx);
    if (!v)
	return;

    if (do_bool && rec.type == QG_POLYGON_GENERAL && close_general_poly->isChecked()) {
	std::vector<qg_polygon_ref> targets;
	struct _qpolymod_poly_collect pc;
	pc.polys = &targets;
	pc.exclude = p;
	ged_view_polygon_visit_records(qpolymod_ged_view(v),
		_qpolymod_poly_collect_cb, &pc);

	enum bg_polygon_boolean_op op = BG_POLYGON_BOOLEAN_UNION;
	if (do_bool) {
	    if (csg_modes->currentText() == "Subtraction") {
		op = BG_POLYGON_BOOLEAN_DIFFERENCE;
	    }
	    if (csg_modes->currentText() == "Intersection") {
		op = BG_POLYGON_BOOLEAN_INTERSECTION;
	    }
	}
	// If we're closing a general polygon and we're in boolean op mode,
	// that's our signal to complete the operation
	int pcnt = 0;
	std::vector<qg_polygon_ref> cleanup;
	for (auto target : targets) {
	    pcnt += ged_view_polygon_csg(qpolymod_current_ged_view(m_ctx), target, p, op);
	    qg_polygon_record target_rec;
	    if (!ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), target, &target_rec) || !target_rec.contour_count)
		cleanup.push_back(target);
	}
	for (size_t i = 0; i < cleanup.size(); i++) {
	    ged_view_polygon_remove(qpolymod_current_ged_view(m_ctx), cleanup[i]);
	}
	if (pcnt || op != BG_POLYGON_BOOLEAN_UNION) {
	    ged_view_polygon_remove(qpolymod_current_ged_view(m_ctx), p);
	    p = GED_VIEW_POLYGON_REF_NULL;
	}
	do_bool = false;
    }

    if (!ged_view_polygon_ref_is_null(p))
	ged_view_polygon_update(qpolymod_ged_view(v), p,
		QG_POLYGON_UPDATE_DEFAULT);

    toplevel_config(false);

    emit view_updated(QG_VIEW_REFRESH);
}

void
QPolyMod::apply_bool_op()
{
    struct ged *gedp = getGed();
    if (ged_view_polygon_ref_is_null(p))
	return;
    if (!gedp)
	return;
    void *v = qpolymod_view(m_ctx);
    if (!v)
	return;

    qg_polygon_record rec;
    if (!ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), p, &rec) || !rec.contour_count)
	return;

    if (rec.first_contour_open) {
	return;
    }

    enum bg_polygon_boolean_op op = BG_POLYGON_BOOLEAN_UNION;
    if (csg_modes->currentText() == "Union") {
	op = BG_POLYGON_BOOLEAN_UNION;
    }
    if (csg_modes->currentText() == "Subtraction") {
	op = BG_POLYGON_BOOLEAN_DIFFERENCE;
    }
    if (csg_modes->currentText() == "Intersection") {
	op = BG_POLYGON_BOOLEAN_INTERSECTION;
    }

    std::vector<qg_polygon_ref> targets;
    struct _qpolymod_poly_collect pc;
    pc.polys = &targets;
    pc.exclude = p;
    ged_view_polygon_visit_records(qpolymod_ged_view(v),
	    _qpolymod_poly_collect_cb, &pc);

    std::vector<qg_polygon_ref> cleanup;
    for (auto target : targets) {
	ged_view_polygon_csg(qpolymod_current_ged_view(m_ctx), target, p, op);
	qg_polygon_record target_rec;
	if (!ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), target, &target_rec) || !target_rec.contour_count)
	    cleanup.push_back(target);
    }
    for (size_t i = 0; i < cleanup.size(); i++) {
	ged_view_polygon_remove(qpolymod_current_ged_view(m_ctx), cleanup[i]);
    }

    mod_names_reset();
    emit view_updated(QG_VIEW_REFRESH);
}

void
QPolyMod::align_to_poly()
{
    struct ged *gedp = getGed();
    if (ged_view_polygon_ref_is_null(p))
	return;
    if (!gedp)
	return;

    qg_polygon_record rec;
    if (!ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), p, &rec))
	return;

    point_t center;
    vect_t dir = VINIT_ZERO;
    VSET(dir, rec.vp[0], rec.vp[1], rec.vp[2]);
    bg_plane_pt_at(&center, &rec.vp, 0, 0);
    void *v = qpolymod_view(m_ctx);
    if (!v)
	return;

    bv_center_set(bv_context_view(static_cast<struct bv_context *>(v)), center);
    vect_t aet = VINIT_ZERO;
    bn_ae_vec(&aet[0], &aet[1], dir);
    aet[2] = 0;
    bv_aet_set(bv_context_view(static_cast<struct bv_context *>(v)), aet);
    bv_context_update(static_cast<struct bv_context *>(v), BV_CONTEXT_CHANGED_VIEW);

    emit view_updated(QG_VIEW_REFRESH);
}

void
QPolyMod::delete_poly()
{
    struct ged *gedp = getGed();
    if (ged_view_polygon_ref_is_null(p))
	return;
    if (!gedp)
	return;

    ged_view_polygon_remove(qpolymod_current_ged_view(m_ctx), p);
    p = GED_VIEW_POLYGON_REF_NULL;
    mod_names_reset();

    emit view_updated(QG_VIEW_REFRESH);
}

void
QPolyMod::sketch_sync_bool(bool)
{
    sketch_name_edit();
}

void
QPolyMod::sketch_name_edit_str(const QString &)
{
    sketch_name_edit();
}


void
QPolyMod::sketch_name_edit()
{
    struct ged *gedp = getGed();
    if (!gedp) {
	ps->sketch_name->setPlaceholderText("No .g file open");
	ps->sketch_name->setStyleSheet("color: rgb(200,200,200)");
	ps->sketch_name->setEnabled(false);
	return;
    }

    if (ps->sketch_sync->isChecked()) {
	char *sname = NULL;
	if (!ps->sketch_name->placeholderText().length()) {
	    if (ps->view_name->placeholderText().length()) {
		ps->sketch_name->setPlaceholderText(ps->view_name->placeholderText());
		sname = bu_strdup(ps->sketch_name->placeholderText().toLocal8Bit().data());
	    }
	} else {
	    sname = bu_strdup(ps->sketch_name->placeholderText().toLocal8Bit().data());
	}
	if (!ps->sketch_name->text().length()) {
	    if (ps->view_name->text().length()) {
		ps->sketch_name->setPlaceholderText(ps->view_name->text());
		if (sname)
		    bu_free(sname, "sname");
		sname = bu_strdup(ps->sketch_name->placeholderText().toLocal8Bit().data());
	    }
	} else {
	    if (sname)
		bu_free(sname, "sname");
	    sname = bu_strdup(ps->sketch_name->text().toLocal8Bit().data());
	}

	void *poly_udata = ged_view_polygon_user_data(qpolymod_current_ged_view(m_ctx), p);
	if (!sname && poly_udata) {
	    struct directory *dp = (struct directory *)poly_udata;
	    ps->sketch_name->setPlaceholderText(QString(dp->d_namep));
	    sname = bu_strdup(ps->sketch_name->placeholderText().toLocal8Bit().data());
	}

	if (sname) {
	    // There may be a dp associated with the polygon.  If so, and the name
	    // matches, then a save operation will replace the old copy with an
	    // updated version.
	    bool match_curr = false;
	    poly_udata = ged_view_polygon_user_data(qpolymod_current_ged_view(m_ctx), p);
	    if (poly_udata) {
		struct directory *dp = (struct directory *)poly_udata;
		if (BU_STR_EQUAL(dp->d_namep, sname)) {
		    match_curr = true;
		}
	    }
	    if (!match_curr) {
		// Not a name match to existing dp - fall back on db check.
		if (db_lookup(gedp->dbip, sname, LOOKUP_QUIET) != RT_DIR_NULL) {
		    ps->sketch_name->setStyleSheet("color: rgb(255,0,0)");
		} else {
		    ps->sketch_name->setStyleSheet("");
		}
	    } else {
		ps->sketch_name->setStyleSheet("");
	    }
	    bu_free(sname, "sname");
	} else {
	    ps->sketch_name->setStyleSheet("");
	}
	ps->sketch_name->setEnabled(true);
    } else {
	ps->sketch_name->setPlaceholderText("Enable to save sketch");
	ps->sketch_name->setStyleSheet("color: rgb(200,200,200)");
	ps->sketch_name->setEnabled(false);
    }
}

void
QPolyMod::sketch_name_update()
{
    struct ged *gedp = getGed();
    if (!gedp)
	return;

    if (ged_view_polygon_ref_is_null(p) || !ps->sketch_sync->isChecked()) {
	return;
    }

    char *sk_name = NULL;
    if (!ps->sketch_name->placeholderText().length()) {
	if (ps->view_name->placeholderText().length()) {
	    ps->sketch_name->setPlaceholderText(ps->view_name->placeholderText());
	    sk_name = bu_strdup(ps->sketch_name->placeholderText().toLocal8Bit().data());
	}
    } else {
	sk_name = bu_strdup(ps->sketch_name->placeholderText().toLocal8Bit().data());
    }
    if (!ps->sketch_name->text().length()) {
	if (ps->view_name->text().length()) {
	    ps->sketch_name->setPlaceholderText(ps->view_name->text());
	    bu_free(sk_name, "sk_name");
	    sk_name = bu_strdup(ps->sketch_name->placeholderText().toLocal8Bit().data());
	}
    } else {
	bu_free(sk_name, "sk_name");
	sk_name = bu_strdup(ps->sketch_name->text().toLocal8Bit().data());
    }

    if (!sk_name)
	return;

    struct directory *sk_dp = RT_DIR_NULL;
    {
	QgGedEventBatch event_batch(gedp);
	void *poly_udata = ged_view_polygon_user_data(qpolymod_current_ged_view(m_ctx), p);
	if (poly_udata) {
	    // remove previous dp, if name is different.  If name
	    // matches, we're done
	    struct directory *dp = (struct directory *)poly_udata;
	    if (BU_STR_EQUAL(sk_name, dp->d_namep)) {
		bu_free(sk_name, "name copy");
		return;
	    }

	    // If the proposed new name collides with an existing object,
	    // we're done.
	    if (db_lookup(gedp->dbip, sk_name, LOOKUP_QUIET) != RT_DIR_NULL) {
		bu_free(sk_name, "name copy");
		return;
	    }

	    // Passed the tests - remove old object.
	    int ac = 2;
	    const char *av[2] = {"kill", NULL};
	    av[0] = "kill";
	    av[1] = dp->d_namep;
	    ged_exec_kill(gedp, ac, av);
	}

	sk_dp = ged_view_polygon_export_sketch(qpolymod_current_ged_view(m_ctx), gedp->dbip, sk_name, p);
    }

    ged_view_polygon_user_data_set(qpolymod_current_ged_view(m_ctx), p, (void *)sk_dp);
    emit view_updated(QG_VIEW_DB);

    bu_free(sk_name, "name copy");
}


void
QPolyMod::view_name_edit_str(const QString &)
{
    view_name_edit();
}


void
QPolyMod::view_name_edit()
{
    struct ged *gedp = getGed();
    if (!gedp)
	return;

    if (!ps->uniq_obj_name(NULL, m_ctx)) {
	ps->view_name->setStyleSheet("color: rgb(255,0,0)");
    } else {
	ps->view_name->setStyleSheet("");
    }
}

void
QPolyMod::view_name_update()
{
    struct ged *gedp = getGed();
    if (ged_view_polygon_ref_is_null(p))
	return;
    if (!gedp)
	return;
    void *v = qpolymod_view(m_ctx);
    if (!v)
	return;

    // Make sure the name is unique
    struct bu_vls vname = BU_VLS_INIT_ZERO;
    if (!ps->uniq_obj_name(&vname, m_ctx)) {
	bu_vls_free(&vname);
	return;
    }

    ged_view_polygon_set_name(qpolymod_current_ged_view(m_ctx), p, bu_vls_cstr(&vname));
    bu_vls_free(&vname);
    mod_names_reset();
    emit view_updated(QG_VIEW_REFRESH);
}


void
QPolyMod::toggle_line_snapping(bool s)
{
    void *v = qpolymod_view(m_ctx);
    qg_polygon_ref co = (cf) ? cf->polygon : GED_VIEW_POLYGON_REF_NULL;
    if (!v || ged_view_polygon_ref_is_null(co))
	return;

    struct bv *view = bv_context_view(static_cast<struct bv_context *>(v));
    bv_snap_source_flags_set(view, BV_SNAP_VIEW);
    if (!s) {
	bv_snap_lines_set(view, 0);
    } else {
	ged_view_polygon_snap_exclude_set(qpolymod_ged_view(v), co);
	bv_snap_lines_set(view,
		ged_view_polygon_snap_count(qpolymod_ged_view(v), co) ? 1 : 0);
    }

    emit settings_changed(QG_VIEW_DRAWN);
}

void
QPolyMod::toggle_grid_snapping(bool s)
{
    void *v = qpolymod_view(m_ctx);
    if (!v)
	return;

    struct bv *view = bv_context_view(static_cast<struct bv_context *>(v));
    bv_snap_source_flags_set(view, BV_SNAP_VIEW);
    struct bv_grid_state grid;
    if (!bv_grid_state_get(&grid, view))
	return;
    grid.snap = s ? 1 : 0;
    bv_grid_state_set(view, &grid);

    emit settings_changed(QG_VIEW_DRAWN);
}

void
QPolyMod::checkbox_refresh(unsigned long long)
{
    void *v = qpolymod_view(m_ctx);
    if (!v)
	return;

    ps->grid_snapping->blockSignals(true);
    const struct bv *view = bv_context_view_const(static_cast<const struct bv_context *>(v));
    struct bv_grid_state grid = {};
    (void)bv_grid_state_get(&grid, view);
    if (grid.snap) {
	ps->grid_snapping->setCheckState(Qt::Checked);
    } else {
	ps->grid_snapping->setCheckState(Qt::Unchecked);
    }
    ps->grid_snapping->blockSignals(false);

    ps->line_snapping->blockSignals(true);
    if (bv_snap_lines_get(view)) {
	ps->line_snapping->setCheckState(Qt::Checked);
    } else {
	ps->line_snapping->setCheckState(Qt::Unchecked);
    }
    ps->line_snapping->blockSignals(false);
}

void
QPolyMod::propagate_update(QgViewUpdateFlags)
{
    emit view_updated(QG_VIEW_REFRESH);
}


bool
QPolyMod::eventFilter(QObject *, QEvent *e)
{
    struct ged *gedp = getGed();
    if (!gedp)
	return false;
    QgView *display = qpolymod_view_widget(m_ctx);
    if (!display || !display->viewContext())
	return false;

    // We might be selecting or modifying - if the former, we may
    // not have a current polygon.
    qg_polygon_record rec;
    int have_rec = ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), p, &rec);

    // The mouse filter to use depends on the mode - find out
    cf = puf;
    if (select_mode->isChecked())
	cf = psf;
    if (move_mode->isChecked())
	cf = pmf;
    if (update_mode->isChecked() && have_rec && rec.type == QG_POLYGON_GENERAL) {
	if (append_pnt->isChecked() || select_pnt->isChecked())
	    cf = ppf;
    }

    // Set libqtcad know what the current polygon is
    cf->polygon = p;
    cf->set_view_widget(display);
    cf->ptype = have_rec ? rec.type : QG_POLYGON_GENERAL;
    checkbox_refresh(0);

    // Connect whatever the current filter is to pass on updating signals from
    // the libqtcad logic.
    QObject::connect(cf, &QgPolyFilter::view_updated, this, &QPolyMod::propagate_update);

    // Match the settings
    cf->op = BG_POLYGON_BOOLEAN_NONE;
    cf->fill_poly = (ps->fill_poly->isChecked()) ? true : false;
    cf->fill_slope_x = (fastf_t)(ps->fill_slope_x->text().toDouble());
    cf->fill_slope_y = (fastf_t)(ps->fill_slope_y->text().toDouble());
    cf->fill_density = (fastf_t)(ps->fill_density->text().toDouble());
    BU_COLOR_CPY(&cf->fill_color, &ps->fill_color->bc);
    BU_COLOR_CPY(&cf->edge_color, &ps->edge_color->bc);
    cf->vZ = (fastf_t)(ps->vZ->text().toDouble());

    // Run the guts of the libqtcad filter
    bool ret = cf->eventFilter(NULL, e);

    // Retrieve the polygon ref from the libqtcad data container.
    p = cf->polygon;
    have_rec = ged_view_polygon_record_get(qpolymod_current_ged_view(m_ctx), p, &rec);

    // If we need to, update our selected list entry
    if (select_mode->isChecked() && have_rec) {
	int cind = mod_names->findText(rec.name ? rec.name : "");
	mod_names->blockSignals(true);
	mod_names->setCurrentIndex(cind);
	mod_names->blockSignals(false);
	poly_type_settings(&rec);
	ps->settings_sync(&rec);
    }

    // Because the active filter may change, we only maintain the
    // signal connection for the duration of the event
    QObject::disconnect(cf, &QgPolyFilter::view_updated, this, &QPolyMod::propagate_update);

    return ret;
}


// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
