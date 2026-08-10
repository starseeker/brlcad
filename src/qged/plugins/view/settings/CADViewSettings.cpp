/*             C A D V I E W S E T T I N G S . C P P
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
/** @file CADViewSettings.cpp
 *
 * Widget implementation for viewing and controlling faceplate settings.
 *
 */

#include <QHBoxLayout>
#include <QLabel>
#include <QVBoxLayout>

#include "bu/opt.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "BObol/BDisplayEndpoint.h"
#include "ged/draw.h"
#include "ged/view.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgSignalFlags.h"

#include "CADViewSettings.h"

static struct bv_context *
qged_settings_view(const QgPluginContext *ctx)
{
    return ctx ? ctx->activeViewContext() : nullptr;
}

static int
qged_framebuffer_mode_get(const struct ged_view_context *view_ctx, int *mode)
{
    if (!view_ctx || !mode)
	return 0;

    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    if (ged_view_context_display_property_get(view_ctx,
	    "composition.framebuffer.mode", &value) !=
	BV_DISPLAY_PROPERTY_OK || !value.string_value)
	return 0;

    if (BU_STR_EQUAL(value.string_value, "off"))
	*mode = 0;
    else if (BU_STR_EQUAL(value.string_value, "overlay"))
	*mode = 1;
    else if (BU_STR_EQUAL(value.string_value, "underlay"))
	*mode = 2;
	else if (BU_STR_EQUAL(value.string_value, "interlay"))
	*mode = 3;
    else
	return 0;
    return 1;
}

static int
qged_framebuffer_mode_set(struct ged_view_context *view_ctx, int mode)
{
    static const char *const modes[] = {
	"off", "overlay", "underlay", "interlay"
    };
    if (!view_ctx || mode < 0 ||
	mode >= static_cast<int>(sizeof(modes) / sizeof(modes[0])))
	return 0;

    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_ENUM;
    value.string_value = modes[mode];
    return ged_view_context_display_property_set(view_ctx,
	"composition.framebuffer.mode", &value) ==
	BV_DISPLAY_PROPERTY_OK;
}

static int
qged_faceplate_property_get(const struct ged_view_context *view_ctx, const char *name,
	int *enabled)
{
    if (!view_ctx || !name || !enabled)
	return 0;

    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    if (ged_view_context_display_property_get(view_ctx, name, &value) !=
	BV_DISPLAY_PROPERTY_OK)
	return 0;
    *enabled = value.bool_value ? 1 : 0;
    return 1;
}

static int
qged_faceplate_property_set(struct ged_view_context *view_ctx, const char *name, int enabled)
{
    if (!view_ctx || !name)
	return 0;

    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_BOOL;
    value.bool_value = enabled ? 1 : 0;
    return ged_view_context_display_property_set(view_ctx, name, &value) ==
	BV_DISPLAY_PROPERTY_OK;
}

/* Helper: update a checkbox to reflect an integer flag, blocking signals
 * to prevent triggering a spurious view_refresh round-trip. */
static void
set_ckbx(QCheckBox *cb, int val)
{
    cb->blockSignals(true);
    cb->setCheckState(val ? Qt::Checked : Qt::Unchecked);
    cb->blockSignals(false);
}

/* Helper: return 1 if the checkbox is checked, 0 otherwise. */
static int
ckbx_val(QCheckBox *cb)
{
    return (cb->checkState() == Qt::Checked) ? 1 : 0;
}

static void
qged_faceplate_checkbox_refresh(const struct ged_view_context *view_ctx, QCheckBox *checkbox,
	const char *property)
{
    int enabled = 0;
    const int supported = qged_faceplate_property_get(view_ctx, property,
	&enabled);
    set_ckbx(checkbox, enabled);
    checkbox->setEnabled(supported ? true : false);
}

CADViewSettings::CADViewSettings(QWidget *)
{
    QVBoxLayout *wl = new QVBoxLayout;
    wl->setAlignment(Qt::AlignTop);

    /* Top-level faceplate toggles */
    acsg_ckbx = new QCheckBox("Adaptive Plotting (CSG)");
    amesh_ckbx = new QCheckBox("Adaptive Plotting (Mesh)");
    adc_ckbx = new QCheckBox("Angle/Dist. Cursor");
    cdot_ckbx = new QCheckBox("Center Dot");
    grid_ckbx = new QCheckBox("Grid");
    mdlaxes_ckbx = new QCheckBox("Model Axes");
    scale_ckbx = new QCheckBox("Scale");
    viewaxes_ckbx = new QCheckBox("View Axes");

    /* Framebuffer mode selector: Off / Overlay / Underlay / Interlay */
    QHBoxLayout *fbl = new QHBoxLayout;
    fbl->addWidget(new QLabel("Framebuffer:"));
    fb_mode_combo = new QComboBox;
    fb_mode_combo->addItem("Off");      /* index 0 -> gv_fb_mode = 0 */
    fb_mode_combo->addItem("Overlay");  /* index 1 -> gv_fb_mode = 1 */
    fb_mode_combo->addItem("Underlay"); /* index 2 -> gv_fb_mode = 2 */
    fb_mode_combo->addItem("Interlay"); /* index 3 -> Obol interlay */
    fbl->addWidget(fb_mode_combo);
    fbl->addStretch();

    /* Parameters group: master toggle + per-element sub-flags */
    params_grp = new QGroupBox("Parameters");
    QVBoxLayout *pgl = new QVBoxLayout;
    params_ckbx = new QCheckBox("Enable");
    params_size_ckbx = new QCheckBox("Size");
    params_center_ckbx = new QCheckBox("Center");
    params_az_ckbx = new QCheckBox("Azimuth");
    params_el_ckbx = new QCheckBox("Elevation");
    params_tw_ckbx = new QCheckBox("Twist");
    params_fps_ckbx = new QCheckBox("FPS");
    pgl->addWidget(params_ckbx);
    pgl->addWidget(params_size_ckbx);
    pgl->addWidget(params_center_ckbx);
    pgl->addWidget(params_az_ckbx);
    pgl->addWidget(params_el_ckbx);
    pgl->addWidget(params_tw_ckbx);
    pgl->addWidget(params_fps_ckbx);
    params_grp->setLayout(pgl);

    /* Wire up signals -> view_refresh.
     * QCheckBox::stateChanged(int) works across all Qt5/Qt6 versions.
     * checkStateChanged(Qt::CheckState) was only added in Qt 6.7, so we
     * use the version-specific signal to avoid deprecation warnings when
     * building with a sufficiently new Qt. */
#if QT_VERSION < QT_VERSION_CHECK(6, 7, 0)
    QObject::connect(acsg_ckbx,         &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(amesh_ckbx,        &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(adc_ckbx,          &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(cdot_ckbx,         &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(grid_ckbx,         &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(mdlaxes_ckbx,      &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(scale_ckbx,        &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(viewaxes_ckbx,     &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_ckbx,       &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_size_ckbx,  &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_center_ckbx,&QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_az_ckbx,    &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_el_ckbx,    &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_tw_ckbx,    &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_fps_ckbx,   &QCheckBox::stateChanged, this, &CADViewSettings::view_update_int);
#else
    QObject::connect(acsg_ckbx,         &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(amesh_ckbx,        &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(adc_ckbx,          &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(cdot_ckbx,         &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(grid_ckbx,         &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(mdlaxes_ckbx,      &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(scale_ckbx,        &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(viewaxes_ckbx,     &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_ckbx,       &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_size_ckbx,  &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_center_ckbx,&QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_az_ckbx,    &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_el_ckbx,    &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_tw_ckbx,    &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
    QObject::connect(params_fps_ckbx,   &QCheckBox::checkStateChanged, this, &CADViewSettings::view_update_int);
#endif
    QObject::connect(fb_mode_combo,
		     QOverload<int>::of(&QComboBox::currentIndexChanged),
		     this, &CADViewSettings::view_update_int);

    /* Assemble the top-level layout */
    wl->addWidget(acsg_ckbx);
    wl->addWidget(amesh_ckbx);
    wl->addWidget(adc_ckbx);
    wl->addWidget(cdot_ckbx);
    wl->addWidget(grid_ckbx);
    wl->addWidget(mdlaxes_ckbx);
    wl->addWidget(scale_ckbx);
    wl->addWidget(viewaxes_ckbx);
    wl->addLayout(fbl);
    wl->addWidget(params_grp);

    this->setLayout(wl);
}

CADViewSettings::~CADViewSettings()
{
}

void
CADViewSettings::checkbox_update()
{
    checkbox_refresh(0);
}

void
CADViewSettings::view_update()
{
    view_refresh(0);
}

void
CADViewSettings::view_update_int(int)
{
    view_refresh(0);
}

/* Read current view state and update all widgets to match, without
 * triggering a spurious write-back via the signal connections. */
void
CADViewSettings::checkbox_refresh(unsigned long long)
{
    struct bv_context *bv_ctx = qged_settings_view(m_ctx);
    if (!bv_ctx)
	return;
    struct ged_view_context *view_ctx = ged_view_context_from_bv(bv_ctx);

    ged_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
    (void)ged_view_lod_policy_get(&lod_policy, view_ctx);

    set_ckbx(acsg_ckbx,     lod_policy.csg_enabled);
    set_ckbx(amesh_ckbx,    lod_policy.mesh_enabled);
    qged_faceplate_checkbox_refresh(view_ctx, adc_ckbx,
	"view.faceplate.adc.visible");
    qged_faceplate_checkbox_refresh(view_ctx, cdot_ckbx,
	"view.faceplate.center_dot.visible");
    qged_faceplate_checkbox_refresh(view_ctx, grid_ckbx,
	"view.faceplate.grid.visible");
    qged_faceplate_checkbox_refresh(view_ctx, mdlaxes_ckbx,
	"view.faceplate.model_axes.visible");
    qged_faceplate_checkbox_refresh(view_ctx, scale_ckbx,
	"view.faceplate.scale.visible");
    qged_faceplate_checkbox_refresh(view_ctx, viewaxes_ckbx,
	"view.faceplate.view_axes.visible");

    /* Framebuffer composition is rendered by Obol, while GED retains the
     * requested mode for endpoint-less views awaiting presentation. */
    int fb_mode = 0;
    const int framebuffer_supported =
	qged_framebuffer_mode_get(view_ctx, &fb_mode);
    fb_mode_combo->blockSignals(true);
    fb_mode_combo->setCurrentIndex(fb_mode);
    fb_mode_combo->setEnabled(framebuffer_supported ? true : false);
    fb_mode_combo->blockSignals(false);

    qged_faceplate_checkbox_refresh(view_ctx, params_ckbx,
	"view.faceplate.params.visible");
    qged_faceplate_checkbox_refresh(view_ctx, params_size_ckbx,
	"view.faceplate.params.size");
    qged_faceplate_checkbox_refresh(view_ctx, params_center_ckbx,
	"view.faceplate.params.center");
    qged_faceplate_checkbox_refresh(view_ctx, params_az_ckbx,
	"view.faceplate.params.azimuth");
    qged_faceplate_checkbox_refresh(view_ctx, params_el_ckbx,
	"view.faceplate.params.elevation");
    qged_faceplate_checkbox_refresh(view_ctx, params_tw_ckbx,
	"view.faceplate.params.twist");
    qged_faceplate_checkbox_refresh(view_ctx, params_fps_ckbx,
	"view.faceplate.params.fps");
}

/* Read all widget states and write them back to the view, then signal
 * that the view needs to be redrawn. */
void
CADViewSettings::view_refresh(unsigned long long)
{
    struct bv_context *bv_ctx = qged_settings_view(m_ctx);
    if (!bv_ctx)
	return;
    struct ged_view_context *view_ctx = ged_view_context_from_bv(bv_ctx);

    /* Preserve non-widget LoD policy fields and update only the settings
     * owned by this widget. */
    ged_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
    (void)ged_view_lod_policy_get(&lod_policy, view_ctx);
    lod_policy.csg_enabled = ckbx_val(acsg_ckbx);
    lod_policy.mesh_enabled = ckbx_val(amesh_ckbx);
    lod_policy.zoom_refresh =
	lod_policy.csg_enabled || lod_policy.mesh_enabled;
    (void)ged_view_lod_policy_apply(view_ctx, &lod_policy);
    (void)qged_framebuffer_mode_set(view_ctx, fb_mode_combo->currentIndex());
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.adc.visible", ckbx_val(adc_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.center_dot.visible", ckbx_val(cdot_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.grid.visible", ckbx_val(grid_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.model_axes.visible", ckbx_val(mdlaxes_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.scale.visible", ckbx_val(scale_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.view_axes.visible", ckbx_val(viewaxes_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.params.visible", ckbx_val(params_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.params.size", ckbx_val(params_size_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.params.center", ckbx_val(params_center_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.params.azimuth", ckbx_val(params_az_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.params.elevation", ckbx_val(params_el_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.params.twist", ckbx_val(params_tw_ckbx));
    (void)qged_faceplate_property_set(view_ctx,
	"view.faceplate.params.fps", ckbx_val(params_fps_ckbx));

    emit settings_changed(QG_VIEW_DRAWN);
}


// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
