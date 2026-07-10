/*         Q G O B O L O V E R L A Y S Y N C P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolOverlaySyncPrivate.h
 *
 * Private qtcad/qged bridge for GED diagnostic overlay publications.
 */

#ifndef QGOBOLOVERLAYSYNCPRIVATE_H
#define QGOBOLOVERLAYSYNCPRIVATE_H

#include "qtcad/defines.h"

class QgView;
struct bg_line_layer_builder;
struct ged;
struct ged_diagnostic_hud_label;

QTCAD_EXPORT int qg_obol_sync_line_layer_overlay(struct ged *gedp,
	const char *name,
	const struct bg_line_layer_builder *builder,
	QgView *display);

QTCAD_EXPORT int qg_obol_sync_hud_label_overlay(struct ged *gedp,
	const struct ged_diagnostic_hud_label *label,
	QgView *display);

#endif /* QGOBOLOVERLAYSYNCPRIVATE_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
