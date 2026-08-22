/*      T E S T _ Q T C A D _ O B O L _ F A C E P L A T E _ S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"

#include "bu/str.h"

#include "bv.h"

#include "BObol/BADC.h"
#include "BObol/BAxes.h"
#include "BObol/BGrid.h"
#include "BObol/BViewController.h"
#include "BObol/BVListShape.h"
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
scene_group(BObolViewController *controller)
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
	    if (bu_strcmp(grid->overlayId.getValue().getString(), overlayId) == 0)
		return node;
	}
	if (node->isOfType(SoBRLAxes::getClassTypeId())) {
	    SoBRLAxes *axes = static_cast<SoBRLAxes *>(node);
	    if (bu_strcmp(axes->overlayId.getValue().getString(), overlayId) == 0)
		return node;
	}
	if (node->isOfType(SoBRLADC::getClassTypeId())) {
	    SoBRLADC *adc = static_cast<SoBRLADC *>(node);
	    if (bu_strcmp(adc->overlayId.getValue().getString(), overlayId) == 0)
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
    bu_setprogname(argv[0]);
    QApplication app(argc, argv);

    QgView view(NULL, QgViewType::SW);
    view.resize(160, 120);

    if (!view.viewContext())
	FAIL("QgView should expose its GED view context");

    BObolViewController *controller = view.obolViewController();
    SoGroup *root = scene_group(controller);
    if (!root)
	FAIL("QgView should expose an Obol scene root");

    struct bv_grid_state grid = {};
    struct bv_axes_state modelAxes = {};
    struct bv_axes_state viewAxes = {};
    struct bv_adc_state adc = {};
    struct bv *bv = bv_context_view(static_cast<struct bv_context *>(view.viewContext()));
    (void)bv_grid_state_get(&grid, bv);
    (void)bv_model_axes_state_get(&modelAxes, bv);
    (void)bv_view_axes_state_get(&viewAxes, bv);
    (void)bv_adc_state_get(&adc, bv);

    grid.draw = 1;
    VSET(grid.anchor, 1.0, 2.0, 3.0);
    grid.res_h = 2.0;
    grid.res_v = 4.0;
    grid.res_major_h = 3;
    grid.res_major_v = 2;
    bv_grid_state_set(bv, &grid);

    modelAxes.draw = 1;
    VSET(modelAxes.axes_pos, 4.0, 5.0, 6.0);
    modelAxes.axes_size = 7.0;
    bv_model_axes_state_set(bv, &modelAxes);

    viewAxes.draw = 1;
    VSET(viewAxes.axes_pos, -0.8, -0.7, 0.0);
    viewAxes.axes_size = 0.5;
    bv_view_axes_state_set(bv, &viewAxes);

    adc.draw = 1;
    VSET(adc.pos_model, 8.0, 9.0, 0.0);
    adc.a1 = 30.0;
    adc.dst = 12.0;
    VSET(adc.line_color, 10, 20, 30);
    VSET(adc.tick_color, 40, 50, 60);
    adc.line_width = 3;
    bv_adc_state_set(bv, &adc);

    controller->clearRenderRequest();
    view.need_update(QG_VIEW_DRAWN);
    if (!controller->isRenderRequested())
	FAIL("faceplate sync should request an Obol render");
    if (bu_strcmp(controller->getRenderReason().getString(),
	    "qt-semantic-refresh") != 0)
	FAIL("draw refresh should retain its semantic Obol render reason");

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
	!near_float(obolGrid->spacingV.getValue(), 4.0f) ||
	obolGrid->majorDivisionsH.getValue() != 3 ||
	obolGrid->majorDivisionsV.getValue() != 2 ||
	obolGrid->getMinorSegmentCount() <= 0 ||
	obolGrid->getMajorSegmentCount() <= 0 ||
	obolGrid->getAxisSegmentCount() != 2 ||
	obolGrid->getTotalSegmentCount() !=
	obolGrid->getMinorSegmentCount() +
	obolGrid->getMajorSegmentCount() +
	obolGrid->getAxisSegmentCount())
	FAIL("grid state should map to adaptive Obol grid hierarchy");
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
	!obolAdc->getTickGeometryShape() ||
	obolAdc->getGeometryShape()->getSegmentCount() != 3 ||
	obolAdc->getTickGeometryShape()->getSegmentCount() != 1 ||
	!near_float(obolAdc->lineColor.getValue()[0], 10.0f / 255.0f) ||
	!near_float(obolAdc->tickColor.getValue()[2], 60.0f / 255.0f) ||
	obolAdc->lineWidth.getValue() != 3)
	FAIL("ADC state should map to separately styled Obol geometry");

    grid.draw = 0;
    modelAxes.draw = 0;
    viewAxes.draw = 0;
    adc.draw = 0;
    bv_grid_state_set(bv, &grid);
    bv_model_axes_state_set(bv, &modelAxes);
    bv_view_axes_state_set(bv, &viewAxes);
    bv_adc_state_set(bv, &adc);
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
