/*                  Q G S C E N E S Y N C P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgSceneSyncPrivate.h
 *
 * Private qtcad bridge from committed GED scene deltas to view invalidation.
 */

#ifndef QGSCENESYNCPRIVATE_H
#define QGSCENESYNCPRIVATE_H

#include "qtcad/defines.h"

class QgView;
struct ged;
struct ged_scene_delta;

QTCAD_EXPORT int qg_scene_bind(struct ged *gedp, QgView *display);

QTCAD_EXPORT int qg_scene_delta_accepts_view(
	const struct ged_scene_delta *delta, QgView *display);

QTCAD_EXPORT int qg_scene_delta_has_view(
	const struct ged_scene_delta *delta);

QTCAD_EXPORT int qg_scene_delta_notify(
	struct ged *gedp, const struct ged_scene_delta *delta, QgView *display);

#endif /* QGSCENESYNCPRIVATE_H */

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
