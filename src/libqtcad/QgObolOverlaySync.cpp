/*             Q G O B O L O V E R L A Y S Y N C . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolOverlaySync.cpp */

#include "common.h"

#include "QgObolOverlaySyncPrivate.h"

#include "brlobol/view_controller.h"
#include "ged.h"
#include "QgObolViewSyncPrivate.h"
#include "qtcad/QgView.h"

static int
clamp_rgb(int c)
{
    if (c < 0)
	return 0;
    if (c > 255)
	return 255;
    return c;
}

int
qg_obol_sync_line_layer_overlay(struct ged *gedp,
	const char *name,
	const struct bg_line_layer_builder *builder,
	QgView *display)
{
    if (!display || !name || !builder)
	return 0;

    if (!qg_obol_display_accepts_ged_active_view(gedp, display))
	return 0;

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    int realized = obol->replaceLineLayerOverlay(name, builder, 0, TRUE);
    if (realized < 0)
	return 0;

    display->need_update(QG_VIEW_DRAWN);
    return 1;
}

int
qg_obol_sync_hud_label_overlay(struct ged *gedp,
	const struct ged_diagnostic_hud_label *label,
	QgView *display)
{
    if (!display || !label || !label->label_id)
	return 0;

    if (!qg_obol_display_accepts_ged_active_view(gedp, display))
	return 0;

    BRLObolViewController *obol = display->obolViewController();
    if (!obol)
	return 0;

    const SbVec2f position(static_cast<float>(label->position[0]),
	    static_cast<float>(label->position[1]));
    const SbColor color(static_cast<float>(clamp_rgb(label->color[0])) / 255.0f,
	    static_cast<float>(clamp_rgb(label->color[1])) / 255.0f,
	    static_cast<float>(clamp_rgb(label->color[2])) / 255.0f);
    const int realized = obol->replaceHUDLabelOverlay(label->label_id,
	    label->text, position, color, static_cast<float>(label->font_size),
	    label->source_id);
    if (realized < 0)
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
