/*                 Q P O L Y S E T T I N G S . C P P
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
/** @file polygon_settings.cpp
 *
 */

#include <climits>

#include <QLabel>
#include <QString>
#include <QMessageBox>
#include "bu/malloc.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "bu/color.h"
#include "bv/view.h"
#include "ged/draw.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgView.h"
#include "QPolySettings.h"

static void *
qpolysettings_view(const QgPluginContext *ctx)
{
    return ctx ? ctx->activeViewContext() : nullptr;
}

static int
qpolysettings_unique_polygon_name(struct bu_vls *oname,
	const char *seed,
	void *view)
{
    if (!oname || !view)
	return 0;

    struct bu_vls vseed = BU_VLS_INIT_ZERO;
    if (seed && seed[0]) {
	bu_vls_sprintf(&vseed, "%s", seed);
    } else {
	const char *view_name = bv_name_get(bv_context_view_const(static_cast<const struct bv_context *>(view)));
	bu_vls_sprintf(&vseed, "%s:obj_0", view_name ? view_name : "view");
    }

    const char *npattern = "([-_:]*[0-9]+[-_:]*)[^0-9]*$";
    long int loop_guard = 0;
    while (!ged_view_polygon_ref_is_null(
	    ged_view_polygon_find(ged_view_context_from_bv(
		    static_cast<struct bv_context *>(view)),
		    bu_vls_cstr(&vseed))) &&
	    loop_guard < LONG_MAX) {
	(void)bu_vls_incr(&vseed, npattern, NULL, NULL, NULL);
	loop_guard++;
    }

    if (loop_guard >= LONG_MAX) {
	bu_vls_free(&vseed);
	return 0;
    }

    bu_vls_sprintf(oname, "%s", bu_vls_cstr(&vseed));
    bu_vls_free(&vseed);
    return 1;
}

QPolySettings::QPolySettings()
    : QWidget()
{
    QVBoxLayout *l = new QVBoxLayout;
    l->setAlignment(Qt::AlignTop);
    l->setSpacing(1);
    l->setContentsMargins(1,1,1,1);

    QLabel *vn_label = new QLabel("Polygon object name:");
    view_name = new QLineEdit();
    l->addWidget(vn_label);
    l->addWidget(view_name);

    QLabel *sn_label = new QLabel("Polygon sketch (.g) name:");
    l->addWidget(sn_label);
    QWidget *sn = new QWidget();
    QHBoxLayout *snl = new QHBoxLayout();
    snl->setSpacing(0);
    snl->setContentsMargins(1,1,1,1);
    sketch_sync = new QCheckBox();
    snl->addWidget(sketch_sync);
    sketch_name = new QLineEdit();
    sketch_name->setPlaceholderText("Enable to save sketch");
    sketch_name->setStyleSheet("color: rgb(200,200,200)");
    sketch_name->setEnabled(false);
    snl->addWidget(sketch_name);
    QObject::connect(sketch_sync, &QCheckBox::toggled, this, &QPolySettings::sketch_sync_toggled);
    sn->setLayout(snl);
    l->addWidget(sn);

    edge_color = new QgColorRGB(this, "Edge:", QColor(Qt::yellow));
    l->addWidget(edge_color);
    QObject::connect(edge_color, &QgColorRGB::color_changed, this, &QPolySettings::do_settings_changed);
    fill_poly = new QCheckBox("Shade polygon interiors");
    l->addWidget(fill_poly);
    QObject::connect(fill_poly, &QCheckBox::toggled, this, &QPolySettings::do_settings_changed);
    fill_color = new QgColorRGB(this, "Fill:", QColor(Qt::blue));
    l->addWidget(fill_color);
    QObject::connect(fill_color, &QgColorRGB::color_changed, this, &QPolySettings::do_settings_changed);

    QFont f("");
    f.setStyleHint(QFont::Monospace);

    QWidget *xw = new QWidget();
    QHBoxLayout *hx = new QHBoxLayout();
    hx->setSpacing(0);
    hx->setContentsMargins(1,1,1,1);
    QLabel *xlbl = new QLabel(QString("<font face=\"monospace\">Fill Slope X:</font>"));
    hx->addWidget(xlbl);
    fill_slope_x = new QLineEdit(QString("1"));
    fill_slope_x->setFont(f);
    QDoubleValidator *hxv = new QDoubleValidator(0, DBL_MAX, 5);
    fill_slope_x->setValidator(hxv);
    QObject::connect(fill_slope_x, &QLineEdit::editingFinished, this, &QPolySettings::do_settings_changed);
    hx->addWidget(fill_slope_x);
    xw->setLayout(hx);
    l->addWidget(xw);

    QWidget *yw = new QWidget();
    QHBoxLayout *hy = new QHBoxLayout();
    hy->setSpacing(0);
    hy->setContentsMargins(1,1,1,1);
    QLabel *ylbl = new QLabel(QString("<font face=\"monospace\">Fill Slope Y:</font>"));
    hy->addWidget(ylbl);
    fill_slope_y = new QLineEdit(QString("1"));
    fill_slope_y->setFont(f);
    QDoubleValidator *hyv = new QDoubleValidator(0, DBL_MAX, 5);
    fill_slope_y->setValidator(hyv);
    QObject::connect(fill_slope_y, &QLineEdit::editingFinished, this, &QPolySettings::do_settings_changed);
    hy->addWidget(fill_slope_y);
    yw->setLayout(hy);
    l->addWidget(yw);


    QWidget *dw = new QWidget();
    QHBoxLayout *hd = new QHBoxLayout();
    hd->setSpacing(0);
    hd->setContentsMargins(1,1,1,1);
    QLabel *dlbl = new QLabel(QString("<font face=\"monospace\">Fill Spacing:</font>"));
    hd->addWidget(dlbl);
    fill_density = new QLineEdit(QString("10"));
    fill_density->setFont(f);
    QDoubleValidator *hdv = new QDoubleValidator(0, DBL_MAX, 5);
    fill_density->setValidator(hdv);
    QObject::connect(fill_density, &QLineEdit::editingFinished, this, &QPolySettings::do_settings_changed);
    hd->addWidget(fill_density);
    dw->setLayout(hd);
    l->addWidget(dw);

    QWidget *zw = new QWidget();
    QHBoxLayout *zl = new QHBoxLayout();
    zl->setSpacing(0);
    zl->setContentsMargins(1,1,1,1);
    QLabel *zlbl = new QLabel("vZ");
    vZ = new QLineEdit(QString("0"));
    QObject::connect(vZ, &QLineEdit::editingFinished, this, &QPolySettings::do_settings_changed);
    zl->addWidget(zlbl);
    zl->addWidget(vZ);
    zw->setLayout(zl);
    l->addWidget(zw);

    line_snapping = new QCheckBox("Line Snapping");
    l->addWidget(line_snapping);
    QObject::connect(line_snapping, &QCheckBox::toggled, this, &QPolySettings::do_line_snapping_changed);

    grid_snapping = new QCheckBox("Grid Snapping");
    l->addWidget(grid_snapping);
    QObject::connect(grid_snapping, &QCheckBox::toggled, this, &QPolySettings::do_grid_snapping_changed);

    this->setLayout(l);
}

QPolySettings::~QPolySettings()
{
}

void
QPolySettings::setTestIdPrefix(const QString &prefix)
{
    setProperty("qgTestId", prefix + QStringLiteral(".settings"));
    edge_color->setProperty(
	"qgTestId", prefix + QStringLiteral(".edge-color"));
    fill_poly->setProperty(
	"qgTestId", prefix + QStringLiteral(".fill-enabled"));
    fill_color->setProperty(
	"qgTestId", prefix + QStringLiteral(".fill-color"));
    fill_slope_x->setProperty(
	"qgTestId", prefix + QStringLiteral(".fill-slope-x"));
    fill_slope_y->setProperty(
	"qgTestId", prefix + QStringLiteral(".fill-slope-y"));
    fill_density->setProperty(
	"qgTestId", prefix + QStringLiteral(".fill-spacing"));
    view_name->setProperty(
	"qgTestId", prefix + QStringLiteral(".view-name"));
    sketch_name->setProperty(
	"qgTestId", prefix + QStringLiteral(".sketch-name"));
    sketch_sync->setProperty(
	"qgTestId", prefix + QStringLiteral(".sketch-sync"));
    vZ->setProperty("qgTestId", prefix + QStringLiteral(".view-z"));
    line_snapping->setProperty(
	"qgTestId", prefix + QStringLiteral(".line-snapping"));
    grid_snapping->setProperty(
	"qgTestId", prefix + QStringLiteral(".grid-snapping"));
}

bool
QPolySettings::uniq_obj_name(struct bu_vls *oname, const QgPluginContext *ctx)
{
    void *v = qpolysettings_view(ctx);
    if (!v)
	return false;

    char *vname = NULL;
    if (view_name->placeholderText().length()) {
	vname = bu_strdup(view_name->placeholderText().toLocal8Bit().data());
    }
    // If we have a full entry, overwrite the placeholder
    if (view_name->text().length()) {
	bu_free(vname, "vname");
	vname = bu_strdup(view_name->text().toLocal8Bit().data());
    }

    // See if the supplied name will collide.  If it will, then reject.  If we want
    // an output name, fail with a message box
    struct bu_vls ovname = BU_VLS_INIT_ZERO;
    if (!qpolysettings_unique_polygon_name(&ovname, vname, v)) {
	if (vname)
	    bu_free(vname, "vname");
	return false;
    }
    const char *candidate = vname ? vname : bu_vls_cstr(&ovname);
    if (!BU_STR_EQUAL(bu_vls_cstr(&ovname), candidate)) {
	QMessageBox msgBox;
	if (oname) {
	    msgBox.setText("Proposed object name already exists in view.");
	    msgBox.exec();
	}
	bu_vls_free(&ovname);
	if (vname)
	    bu_free(vname, "vname");
	return false;
    }

    // Unique.  If we want it returned, do the printing
    if (oname)
	bu_vls_sprintf(oname, "%s", candidate);

    bu_vls_free(&ovname);
    if (vname)
	bu_free(vname, "vname");
    return true;
}

void
QPolySettings::sketch_sync_toggled(bool)
{
    if (sketch_sync->isChecked()) {
	if (sketch_name->placeholderText() == QString("Enable to save sketch"))
	    sketch_name->setPlaceholderText("");
	sketch_name->setStyleSheet("");
	sketch_name->setEnabled(true);
    } else {
	sketch_name->setPlaceholderText("Enable to save sketch");
	sketch_name->setStyleSheet("color: rgb(200,200,200)");
	sketch_name->setEnabled(false);
    }
}

void
QPolySettings::do_settings_changed()
{
    emit settings_changed();
}


void
QPolySettings::do_line_snapping_changed()
{
    emit line_snapping_changed(line_snapping->isChecked());
}

void
QPolySettings::do_grid_snapping_changed()
{
    emit grid_snapping_changed(grid_snapping->isChecked());
}


void
QPolySettings::settings_sync(const struct ged_view_polygon_record *p)
{
    if (!p)
	return;


    edge_color->blockSignals(true);
    edge_color->rgbtext->setText(QString("%1/%2/%3").arg(p->edge_color[0]).arg(p->edge_color[1]).arg(p->edge_color[2]));
    edge_color->blockSignals(false);

    unsigned char frgb[3];
    bu_color_to_rgb_chars(&p->fill_color, frgb);
    fill_color->blockSignals(true);
    fill_color->rgbtext->setText(QString("%1/%2/%3").arg(frgb[0]).arg(frgb[1]).arg(frgb[2]));
    fill_color->blockSignals(false);

    fill_slope_x->blockSignals(true);
    fill_slope_x->setText(QString("%1").arg(p->fill_dir[0]));
    fill_slope_x->blockSignals(false);

    fill_slope_y->blockSignals(true);
    fill_slope_y->setText(QString("%1").arg(p->fill_dir[1]));
    fill_slope_y->blockSignals(false);

    fill_density->blockSignals(true);
    fill_density->setText(QString("%1").arg(p->fill_delta));
    fill_density->blockSignals(false);

    fill_poly->blockSignals(true);
    if (p->fill_flag) {
	fill_poly->setChecked(true);
    } else {
	fill_poly->setChecked(false);
    }
    fill_poly->blockSignals(false);

    vZ->blockSignals(true);
    vZ->setText(QVariant(p->vZ).toString());
    vZ->blockSignals(false);

    // Values set, now update the button colors
    this->blockSignals(true);
    edge_color->set_color_from_text();
    fill_color->set_color_from_text();
    this->blockSignals(false);
}



// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
