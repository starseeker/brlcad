/*              C A D V I E W M E A S U R E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2023-2026 United States Government as represented by
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
/** @file CADViewMeasure.cpp
 *
 * Brief description
 *
 */

#include "common.h"
#include <QMouseEvent>
#include <QVBoxLayout>
#include <QtGlobal>
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgColorRGB.h"
#include "qtcad/QgMeasureFilter.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

#include "bu/opt.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bu/units.h"
#include "bg/aabb_ray.h"
#include "bg/plane.h"
#include "ged.h"

#include "./CADViewMeasure.h"

static void *
qged_measure_view(const QgPluginContext *ctx)
{
    return ctx ? ctx->activeViewContext() : nullptr;
}

const BObolInputActionLayer *
CADViewMeasure::inputActionLayer()
{
    static const unsigned int modifier_mask = BOBOL_INPUT_MOD_SHIFT |
	BOBOL_INPUT_MOD_CONTROL | BOBOL_INPUT_MOD_ALT |
	BOBOL_INPUT_MOD_META;
    static const BObolInputBinding bindings[] = {
	{BOBOL_INPUT_POINTER_PRESS, BOBOL_INPUT_ANY, 0, 0, modifier_mask,
	 10, QG_MEASURE_INPUT_BEGIN},
	{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY, 0,
	 modifier_mask, 10, QG_MEASURE_INPUT_UPDATE},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, 0, 0, modifier_mask,
	 10, QG_MEASURE_INPUT_COMMIT}
    };
    static const BObolInputActionLayer layer = {
	"qged-measure", bindings, sizeof(bindings) / sizeof(bindings[0]),
	CADViewMeasure::inputActionDispatch
    };
    return &layer;
}

static QgView *
qged_measure_view_from_event_object(QObject *object)
{
    QWidget *widget = qobject_cast<QWidget *>(object);
    while (widget) {
	QgView *view = qobject_cast<QgView *>(widget);
	if (view)
	    return view;
	widget = widget->parentWidget();
    }
    return nullptr;
}

CADViewMeasure::CADViewMeasure(QWidget *)
{
    QVBoxLayout *wl = new QVBoxLayout;
    wl->setAlignment(Qt::AlignTop);

    measure_3d = new QCheckBox("Use 3D hit points");
    wl->addWidget(measure_3d);

    QLabel *ml1_label = new QLabel("Measured Length #1:");
    length1_report = new QLineEdit();
    length1_report->setReadOnly(true);
    wl->addWidget(ml1_label);
    wl->addWidget(length1_report);

    QLabel *ml2_label = new QLabel("Measured Length #2:");
    length2_report = new QLineEdit();
    length2_report->setReadOnly(true);
    wl->addWidget(ml2_label);
    wl->addWidget(length2_report);

    report_radians = new QCheckBox("Report angle in radians");
    wl->addWidget(report_radians);
#if QT_VERSION < QT_VERSION_CHECK(6, 7, 0)
    QObject::connect(report_radians, &QCheckBox::stateChanged, this, &CADViewMeasure::adjust_text);
#else
    QObject::connect(report_radians, &QCheckBox::checkStateChanged, this, &CADViewMeasure::adjust_text);
#endif

    ma_label = new QLabel("Measured Angle (deg):");
    angle_report = new QLineEdit();
    angle_report->setReadOnly(true);
    wl->addWidget(ma_label);
    wl->addWidget(angle_report);

    color_2d = new QgColorRGB(this, "2D:", QColor(Qt::yellow));
    wl->addWidget(color_2d);
    color_3d = new QgColorRGB(this, "3D:", QColor(Qt::green));
    wl->addWidget(color_3d);
    QObject::connect(color_2d, &QgColorRGB::color_changed, this, &CADViewMeasure::update_color);
    QObject::connect(color_3d, &QgColorRGB::color_changed, this, &CADViewMeasure::update_color);

    this->setLayout(wl);

    f2d = new QMeasure2DFilter();
    f3d = new QMeasure3DFilter();
    mf = (measure_3d->isChecked()) ? (QgMeasureFilter *)f3d : (QgMeasureFilter *)f2d;
}

CADViewMeasure::~CADViewMeasure()
{
    detachFromView(nullptr);
    delete f2d;
    delete f3d;
}

void
CADViewMeasure::attachToView(QgView *view)
{
    if (!view)
	return;
    if (m_input_view && m_input_view != view)
	detachFromView(m_input_view);

    m_input_view = view;
    m_input_endpoint = view->displayEndpoint();
    if (m_input_endpoint && bobol_display_endpoint_input_action_layer_set(
	m_input_endpoint, CADViewMeasure::inputActionLayer(), this, this))
	return;

    m_input_endpoint = nullptr;
    view->add_event_filter(this);
    m_qt_filter_installed = true;
}

void
CADViewMeasure::detachFromView(QgView *view)
{
    if (!m_input_view || (view && view != m_input_view))
	return;
    if (m_input_endpoint)
	(void)bobol_display_endpoint_input_action_layer_clear_if(
	    m_input_endpoint, this);
    if (m_qt_filter_installed)
	m_input_view->clear_event_filter(this);
    if (f2d)
	f2d->set_view_widget(nullptr);
    if (f3d)
	f3d->set_view_widget(nullptr);
    m_input_view = nullptr;
    m_input_endpoint = nullptr;
    m_qt_filter_installed = false;
}

void
CADViewMeasure::update_color()
{
    if (!mf)
	return;
    if (measure_3d->isChecked()) {
	mf->update_color(&color_3d->bc);
	emit view_updated(QG_VIEW_REFRESH);
	return;
    }
    mf->update_color(&color_2d->bc);
    emit view_updated(QG_VIEW_REFRESH);
}

void
CADViewMeasure::adjust_text_db(void *)
{
    adjust_text();
}

void
CADViewMeasure::adjust_text()
{
    struct ged *gedp = m_ctx ? m_ctx->getGed() : nullptr;
    if (!gedp || !gedp->dbip || !mf)
	return;


    double angle;
    if (report_radians->isChecked()) {
	ma_label->setText("Measured Angle (rad):");
	angle = mf->angle(true);
    } else {
	ma_label->setText("Measured Angle (deg):");
	angle = mf->angle(false);
    }

    struct bu_vls buffer = BU_VLS_INIT_ZERO;
    bu_vls_sprintf(&buffer, "%.15f %s", mf->length1()*gedp->dbip->dbi_base2local, bu_units_string(gedp->dbip->dbi_local2base));
    length1_report->setText(bu_vls_cstr(&buffer));

    bu_vls_sprintf(&buffer, "%.15f %s", mf->length2()*gedp->dbip->dbi_base2local, bu_units_string(gedp->dbip->dbi_local2base));
    length2_report->setText(bu_vls_cstr(&buffer));

    bu_vls_sprintf(&buffer, "%.15f", angle);
    angle_report->setText(bu_vls_cstr(&buffer));

    bu_vls_free(&buffer);
}

void
CADViewMeasure::do_filter_view_update()
{
    adjust_text();
    emit view_updated(QG_VIEW_REFRESH);
}


bool
CADViewMeasure::eventFilter(QObject *o, QEvent *e)
{
    struct ged *gedp = m_ctx ? m_ctx->getGed() : nullptr;
    if (!gedp)
	return false;
    QgView *display = qged_measure_view_from_event_object(o);
    if (!display && m_ctx)
	display = m_ctx->getViewWidget();
    void *v = display ? display->viewContext() : qged_measure_view(m_ctx);
    if (!v)
	return false;

    mf = (measure_3d->isChecked()) ? (QgMeasureFilter *)f3d : (QgMeasureFilter *)f2d;

    mf->set_view_widget(display);
    update_color();

    // Connect whatever the current filter is to pass on updating signals from
    // the libqtcad logic.
    QObject::connect(mf, &QgMeasureFilter::view_updated, this, &CADViewMeasure::do_filter_view_update);

    bool ret = mf->eventFilter(NULL, e);

    QObject::disconnect(mf, &QgMeasureFilter::view_updated, this, &CADViewMeasure::do_filter_view_update);

    return ret;
}

int
CADViewMeasure::inputActionDispatch(void *user_data,
	BObolInputAction action, const BObolInputEvent *event)
{
    CADViewMeasure *measure = static_cast<CADViewMeasure *>(user_data);
    return measure ? measure->applyInputAction(action, event) :
	BOBOL_INPUT_RESULT_ERROR;
}

int
CADViewMeasure::applyInputAction(BObolInputAction action,
	const BObolInputEvent *event)
{
    struct ged *gedp = m_ctx ? m_ctx->getGed() : nullptr;
    if (!event || !gedp || !m_input_view ||
	(action != QG_MEASURE_INPUT_BEGIN &&
	 action != QG_MEASURE_INPUT_UPDATE &&
	 action != QG_MEASURE_INPUT_COMMIT &&
	 action != QG_MEASURE_INPUT_CANCEL))
	return BOBOL_INPUT_RESULT_UNHANDLED;

    mf = measure_3d->isChecked() ?
	(QgMeasureFilter *)f3d : (QgMeasureFilter *)f2d;
    mf->set_view_widget(m_input_view);
    update_color();

    QObject::connect(mf, &QgMeasureFilter::view_updated, this,
	&CADViewMeasure::do_filter_view_update);
    bool handled = mf->semanticInput(action, event);
    QObject::disconnect(mf, &QgMeasureFilter::view_updated, this,
	&CADViewMeasure::do_filter_view_update);

    return handled ? BOBOL_INPUT_RESULT_HANDLED :
	BOBOL_INPUT_RESULT_UNHANDLED;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
