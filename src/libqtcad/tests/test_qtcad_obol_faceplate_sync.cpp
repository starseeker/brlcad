/*      T E S T _ Q T C A D _ O B O L _ F A C E P L A T E _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/adc.h"
#include "brlobol/axes.h"
#include "brlobol/grid.h"
#include "brlobol/view_controller.h"
#include "brlobol/vlist_shape.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"
#include "vmath.h"

#include <Inventor/nodes/SoGroup.h>

#include <QApplication>

#include <math.h>
#include <stdio.h>
#include <string.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static SoGroup *
scene_group(BRLObolViewController *controller)
{
    if (!controller)
	return NULL;
    SoNode *root = controller->getSceneRoot();
    if (!root || !root->isOfType(SoGroup::getClassTypeId()))
	return NULL;
    return static_cast<SoGroup *>(root);
}

static SoNode *
find_overlay(SoGroup *group, const char *overlayId)
{
    if (!group || !overlayId)
	return NULL;
    for (int i = 0; i < group->getNumChildren(); i++) {
	SoNode *node = group->getChild(i);
	if (!node)
	    continue;
	if (node->isOfType(SoBRLGrid::getClassTypeId())) {
	    SoBRLGrid *grid = static_cast<SoBRLGrid *>(node);
	    if (strcmp(grid->overlayId.getValue().getString(), overlayId) == 0)
		return node;
	}
	if (node->isOfType(SoBRLAxes::getClassTypeId())) {
	    SoBRLAxes *axes = static_cast<SoBRLAxes *>(node);
	    if (strcmp(axes->overlayId.getValue().getString(), overlayId) == 0)
		return node;
	}
	if (node->isOfType(SoBRLADC::getClassTypeId())) {
	    SoBRLADC *adc = static_cast<SoBRLADC *>(node);
	    if (strcmp(adc->overlayId.getValue().getString(), overlayId) == 0)
		return node;
	}
    }
    return NULL;
}

static int
near_float(float a, float b)
{
    return fabsf(a - b) < 0.001f;
}

int
main(int argc, char **argv)
{
    QApplication app(argc, argv);

    QgView view(NULL, QgView_SW);
    view.resize(160, 120);

    if (!view.view())
	FAIL("QgView should expose transitional view state");
    if (view.legacyBackendInitialized())
	FAIL("test should start before legacy DM initialization");

    BRLObolViewController *controller = view.obolViewController();
    SoGroup *root = scene_group(controller);
    if (!root)
	FAIL("QgView should expose an Obol scene root");

    struct rt_view_grid_state grid = {};
    struct rt_view_axes_state modelAxes = {};
    struct rt_view_axes_state viewAxes = {};
    struct rt_view_adc_state adc = {};
    (void)qg_legacy_view_grid_state_get(view.view(), &grid);
    (void)qg_legacy_view_model_axes_state_get(view.view(), &modelAxes);
    (void)qg_legacy_view_view_axes_state_get(view.view(), &viewAxes);
    (void)qg_legacy_view_adc_state_get(view.view(), &adc);

    grid.draw = 1;
    VSET(grid.anchor, 1.0, 2.0, 3.0);
    grid.res_h = 2.0;
    grid.res_v = 4.0;
    grid.res_major_h = 3;
    grid.res_major_v = 2;
    qg_legacy_view_grid_state_set(view.view(), &grid);

    modelAxes.draw = 1;
    VSET(modelAxes.axes_pos, 4.0, 5.0, 6.0);
    modelAxes.axes_size = 7.0;
    qg_legacy_view_model_axes_state_set(view.view(), &modelAxes);

    viewAxes.draw = 1;
    VSET(viewAxes.axes_pos, -0.8, -0.7, 0.0);
    viewAxes.axes_size = 0.5;
    qg_legacy_view_view_axes_state_set(view.view(), &viewAxes);

    adc.draw = 1;
    VSET(adc.pos_model, 8.0, 9.0, 0.0);
    adc.a1 = 30.0;
    adc.dst = 12.0;
    qg_legacy_view_adc_state_set(view.view(), &adc);

    controller->clearRenderRequest();
    view.need_update(QG_VIEW_DRAWN);
    if (view.legacyBackendInitialized())
	FAIL("Obol faceplate sync should not require legacy DM initialization");
    if (!controller->isRenderRequested() ||
	    strcmp(controller->getRenderReason().getString(), "faceplate") != 0)
	FAIL("faceplate sync should request an Obol render");

    SoBRLGrid *obolGrid = static_cast<SoBRLGrid *>(
	    find_overlay(root, "faceplate::grid"));
    SoBRLAxes *obolModelAxes = static_cast<SoBRLAxes *>(
	    find_overlay(root, "faceplate::model_axes"));
    SoBRLAxes *obolViewAxes = static_cast<SoBRLAxes *>(
	    find_overlay(root, "faceplate::view_axes"));
    SoBRLADC *obolAdc = static_cast<SoBRLADC *>(
	    find_overlay(root, "faceplate::adc"));

    if (!obolGrid || !obolModelAxes || !obolViewAxes || !obolAdc)
	FAIL("faceplate sync should publish grid, axes, and ADC Obol nodes");
    if (!near_float(obolGrid->center.getValue()[0], 1.0f) ||
	    !near_float(obolGrid->center.getValue()[1], 2.0f) ||
	    !near_float(obolGrid->center.getValue()[2], 3.0f) ||
	    !near_float(obolGrid->spacing.getValue(), 2.0f) ||
	    obolGrid->divisions.getValue() != 3 ||
	    !obolGrid->getGeometryShape() ||
	    obolGrid->getGeometryShape()->getSegmentCount() != 14)
	FAIL("grid state should map to Obol grid geometry");
    if (!near_float(obolModelAxes->origin.getValue()[0], 4.0f) ||
	    !near_float(obolModelAxes->origin.getValue()[1], 5.0f) ||
	    !near_float(obolModelAxes->origin.getValue()[2], 6.0f) ||
	    !near_float(obolModelAxes->size.getValue(), 7.0f) ||
	    !obolModelAxes->getGeometryShape() ||
	    obolModelAxes->getGeometryShape()->getSegmentCount() != 3)
	FAIL("model axes state should map to Obol axes geometry");
    if (!near_float(obolViewAxes->origin.getValue()[0], -0.8f) ||
	    !near_float(obolViewAxes->origin.getValue()[1], -0.7f) ||
	    !near_float(obolViewAxes->size.getValue(), 0.5f) ||
	    !obolViewAxes->getGeometryShape() ||
	    obolViewAxes->getGeometryShape()->getSegmentCount() != 3)
	FAIL("view axes state should map to Obol axes geometry");
    if (!near_float(obolAdc->center.getValue()[0], 8.0f) ||
	    !near_float(obolAdc->center.getValue()[1], 9.0f) ||
	    !near_float(obolAdc->angleDegrees.getValue(), 30.0f) ||
	    !near_float(obolAdc->distance.getValue(), 12.0f) ||
	    !obolAdc->getGeometryShape() ||
	    obolAdc->getGeometryShape()->getSegmentCount() != 4)
	FAIL("ADC state should map to Obol ADC geometry");

    grid.draw = 0;
    modelAxes.draw = 0;
    viewAxes.draw = 0;
    adc.draw = 0;
    qg_legacy_view_grid_state_set(view.view(), &grid);
    qg_legacy_view_model_axes_state_set(view.view(), &modelAxes);
    qg_legacy_view_view_axes_state_set(view.view(), &viewAxes);
    qg_legacy_view_adc_state_set(view.view(), &adc);
    view.need_update(QG_VIEW_DRAWN);

    if (find_overlay(root, "faceplate::grid") ||
	    find_overlay(root, "faceplate::model_axes") ||
	    find_overlay(root, "faceplate::view_axes") ||
	    find_overlay(root, "faceplate::adc"))
	FAIL("disabled faceplate settings should remove Obol nodes");

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
