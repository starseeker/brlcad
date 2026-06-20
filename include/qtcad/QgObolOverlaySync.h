/*               Q G O B O L O V E R L A Y S Y N C . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/QgObolOverlaySync.h */

#ifndef QGOBOLOVERLAYSYNC_H
#define QGOBOLOVERLAYSYNC_H

#include "qtcad/defines.h"

class QgView;
struct bg_line_layer_builder;
struct ged;
struct ged_diagnostic_hud_label;

/**
 * Synchronize one GED diagnostic line-layer publication into a QgView's Obol
 * controller.  Returns non-zero when the overlay was accepted by the Obol
 * scene.
 */
QTCAD_EXPORT int qg_obol_sync_line_layer_overlay(struct ged *gedp,
	const char *name,
	const struct bg_line_layer_builder *builder,
	QgView *display);

/**
 * Synchronize one GED diagnostic HUD label publication into a QgView's Obol
 * controller.  Returns non-zero when the label was accepted by the Obol scene.
 */
QTCAD_EXPORT int qg_obol_sync_hud_label_overlay(struct ged *gedp,
	const struct ged_diagnostic_hud_label *label,
	QgView *display);

#endif /* QGOBOLOVERLAYSYNC_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
