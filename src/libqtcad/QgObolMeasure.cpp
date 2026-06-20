/*               Q G O B O L M E A S U R E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolMeasure.cpp */

#include "common.h"

#include "qtcad/QgObolMeasure.h"

#include "bg/line_layer.h"
#include "brlobol/view_controller.h"
#include "bu/color.h"
#include "qtcad/QgObolPick.h"
#include "qtcad/QgSignalFlags.h"
#include "qtcad/QgView.h"

#include <vector>

int
qg_obol_measure_pick_point(QgView *display,
	int x,
	int y,
	SbVec3f &point,
	std::string *path)
{
    std::vector<QgObolPickRecord> picks;
    if (qg_obol_pick_point(display, x, y, 6.0f, false, picks) <= 0)
	return 0;

    point = picks[0].point;
    if (path)
	*path = picks[0].path;
    return 1;
}

static void
qg_obol_measure_rgb(const struct bu_color *color, int &r, int &g, int &b)
{
    r = 255;
    g = 255;
    b = 0;

    if (!color)
	return;

    unsigned char rgb[3] = {255, 255, 0};
    bu_color_to_rgb_chars(color, rgb);
    r = rgb[0];
    g = rgb[1];
    b = rgb[2];
}

static void
qg_obol_measure_set_point(point_t p, const SbVec3f &point)
{
    VSET(p, point[0], point[1], point[2]);
}

int
qg_obol_measure_update_overlay(QgView *display,
	const char *overlayId,
	const SbVec3f *points,
	int pointCount,
	const struct bu_color *color)
{
    if (!display || !overlayId || !overlayId[0] || !points || pointCount <= 0)
	return 0;

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    if (pointCount > 3)
	pointCount = 3;

    int r = 255;
    int g = 255;
    int b = 0;
    qg_obol_measure_rgb(color, r, g, b);

    struct bg_line_layer_builder *builder = bg_line_layer_builder_create();
    if (!builder)
	return 0;

    int ok = 1;
    point_t p;
    qg_obol_measure_set_point(p, points[0]);
    if (pointCount == 1) {
	ok = bg_line_layer_builder_add(builder, r, g, b, p,
		BG_GEOMETRY_POINT_DRAW);
    } else {
	ok = bg_line_layer_builder_add(builder, r, g, b, p,
		BG_GEOMETRY_LINE_MOVE);
    }

    for (int i = 1; ok && i < pointCount; i++) {
	qg_obol_measure_set_point(p, points[i]);
	ok = bg_line_layer_builder_add(builder, r, g, b, p,
		BG_GEOMETRY_LINE_DRAW);
    }

    int realized = 0;
    if (ok)
	realized = obol->replaceLineLayerOverlay(overlayId, builder, 0,
		FALSE);
    bg_line_layer_builder_free(builder);

    if (realized < 0)
	return 0;

    display->need_update(QG_VIEW_DRAWN);
    return realized > 0 ? realized : 1;
}

int
qg_obol_measure_clear_overlay(QgView *display, const char *overlayId)
{
    if (!display || !overlayId || !overlayId[0])
	return 0;

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    int ret = obol->replaceLineLayerOverlay(overlayId, NULL, 0, FALSE);
    if (ret < 0)
	return 0;

    display->need_update(QG_VIEW_DRAWN);
    return 1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
