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
#include "ged/draw.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgSignalFlags.h"
#include "rt/view.h"

#include "CADViewSettings.h"
#include "QgLegacyViewContext.h"

static qg_legacy_view *
qged_settings_view(const QgPluginContext *ctx)
{
    return ctx ? ctx->activeView() : nullptr;
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

    /* Framebuffer mode selector: Off / Overlay / Underlay */
    QHBoxLayout *fbl = new QHBoxLayout;
    fbl->addWidget(new QLabel("Framebuffer:"));
    fb_mode_combo = new QComboBox;
    fb_mode_combo->addItem("Off");      /* index 0 -> gv_fb_mode = 0 */
    fb_mode_combo->addItem("Overlay");  /* index 1 -> gv_fb_mode = 1 */
    fb_mode_combo->addItem("Underlay"); /* index 2 -> gv_fb_mode = 2 */
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
    qg_legacy_view *v = qged_settings_view(m_ctx);
    if (!v)
	return;
    void *view_ctx = qg_legacy_view_to_context(v);
    struct bv *view = qg_legacy_view_bv(v);

    ged_draw_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
    (void)ged_draw_view_context_lod_policy_get(&lod_policy, view_ctx);

    /* Top-level faceplate elements */
    struct bv_adc_state adc = {};
    struct bv_other_state center_dot = {};
    struct bv_grid_state grid = {};
    struct bv_axes_state model_axes = {};
    struct bv_other_state scale_state = {};
    struct bv_axes_state view_axes = {};
    struct bv_params_state params = {};
    (void)bv_adc_state_get(&adc, view);
    (void)bv_center_dot_state_get(&center_dot, view);
    (void)bv_grid_state_get(&grid, view);
    (void)bv_model_axes_state_get(&model_axes, view);
    (void)bv_scale_overlay_state_get(&scale_state, view);
    (void)bv_view_axes_state_get(&view_axes, view);
    (void)bv_params_state_get(&params, view);

    set_ckbx(acsg_ckbx,     lod_policy.csg_enabled);
    set_ckbx(amesh_ckbx,    lod_policy.mesh_enabled);
    set_ckbx(adc_ckbx,      adc.draw);
    set_ckbx(cdot_ckbx,     center_dot.gos_draw);
    set_ckbx(grid_ckbx,     grid.draw);
    set_ckbx(mdlaxes_ckbx,  model_axes.draw);
    set_ckbx(scale_ckbx,    scale_state.gos_draw);
    set_ckbx(viewaxes_ckbx, view_axes.draw);

    /* Framebuffer mode (0=off, 1=overlay, 2=underlay) maps directly to
     * combo index. Clamp to a valid range in case of unexpected values. */
    int fb_mode = bv_framebuffer_mode_get(view);
    if (fb_mode < 0 || fb_mode > 2)
	fb_mode = 0;
    fb_mode_combo->blockSignals(true);
    fb_mode_combo->setCurrentIndex(fb_mode);
    fb_mode_combo->blockSignals(false);

    /* Parameters group: master draw toggle + per-element sub-flags */
    set_ckbx(params_ckbx,        params.draw);
    set_ckbx(params_size_ckbx,   params.draw_size);
    set_ckbx(params_center_ckbx, params.draw_center);
    set_ckbx(params_az_ckbx,     params.draw_az);
    set_ckbx(params_el_ckbx,     params.draw_el);
    set_ckbx(params_tw_ckbx,     params.draw_tw);
    set_ckbx(params_fps_ckbx,    params.draw_fps);
}

/* Read all widget states and write them back to the view, then signal
 * that the view needs to be redrawn. */
void
CADViewSettings::view_refresh(unsigned long long)
{
    qg_legacy_view *v = qged_settings_view(m_ctx);
    if (!v)
	return;
    void *view_ctx = qg_legacy_view_to_context(v);
    struct bv *view = qg_legacy_view_bv(v);

    /* Preserve non-widget LoD policy fields and update only the settings
     * owned by this widget. */
    ged_draw_view_lod_policy lod_policy = BV_LOD_POLICY_INIT;
    (void)ged_draw_view_context_lod_policy_get(&lod_policy, view_ctx);
    lod_policy.csg_enabled = ckbx_val(acsg_ckbx);
    lod_policy.mesh_enabled = ckbx_val(amesh_ckbx);
    lod_policy.zoom_refresh =
	lod_policy.csg_enabled || lod_policy.mesh_enabled;
    (void)ged_draw_view_context_lod_policy_apply(view_ctx, &lod_policy);
    (void)bv_framebuffer_mode_set(view, fb_mode_combo->currentIndex());

    struct bv_adc_state adc = {};
    struct bv_other_state center_dot = {};
    struct bv_grid_state grid = {};
    struct bv_axes_state model_axes = {};
    struct bv_other_state scale_state = {};
    struct bv_axes_state view_axes = {};
    struct bv_params_state params = {};
    (void)bv_adc_state_get(&adc, view);
    (void)bv_center_dot_state_get(&center_dot, view);
    (void)bv_grid_state_get(&grid, view);
    (void)bv_model_axes_state_get(&model_axes, view);
    (void)bv_scale_overlay_state_get(&scale_state, view);
    (void)bv_view_axes_state_get(&view_axes, view);
    (void)bv_params_state_get(&params, view);

    adc.draw = ckbx_val(adc_ckbx);
    center_dot.gos_draw = ckbx_val(cdot_ckbx);
    grid.draw = ckbx_val(grid_ckbx);
    model_axes.draw = ckbx_val(mdlaxes_ckbx);
    scale_state.gos_draw = ckbx_val(scale_ckbx);
    view_axes.draw = ckbx_val(viewaxes_ckbx);

    params.draw        = ckbx_val(params_ckbx);
    params.draw_size   = ckbx_val(params_size_ckbx);
    params.draw_center = ckbx_val(params_center_ckbx);
    params.draw_az     = ckbx_val(params_az_ckbx);
    params.draw_el     = ckbx_val(params_el_ckbx);
    params.draw_tw     = ckbx_val(params_tw_ckbx);
    params.draw_fps    = ckbx_val(params_fps_ckbx);

    bv_adc_state_set(view, &adc);
    bv_center_dot_state_set(view, &center_dot);
    bv_grid_state_set(view, &grid);
    bv_model_axes_state_set(view, &model_axes);
    bv_scale_overlay_state_set(view, &scale_state);
    bv_view_axes_state_set(view, &view_axes);
    bv_params_state_set(view, &params);

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
