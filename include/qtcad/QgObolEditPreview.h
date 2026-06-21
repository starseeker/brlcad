/*              Q G O B O L E D I T P R E V I E W . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file qtcad/QgObolEditPreview.h */

#ifndef QGOBOL_EDIT_PREVIEW_H
#define QGOBOL_EDIT_PREVIEW_H

#include "qtcad/defines.h"

#include <Inventor/SbVec3f.h>

#include <stdint.h>

class QgView;

enum QgObolEditPreviewCommand {
    QG_OBOL_EDIT_PREVIEW_MOVE = 0,
    QG_OBOL_EDIT_PREVIEW_DRAW = 1,
    QG_OBOL_EDIT_PREVIEW_POINT = 2
};

QTCAD_EXPORT int qg_obol_edit_preview_update(QgView *display,
	const char *previewId,
	const char *identity,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	uint32_t sourceRevision,
	uint32_t inputsRevision);

QTCAD_EXPORT int qg_obol_edit_preview_update_with_intent(QgView *display,
	const char *previewId,
	const char *identity,
	const char *editIntentId,
	const char *editIntentRole,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	uint32_t sourceRevision,
	uint32_t inputsRevision);

QTCAD_EXPORT int qg_obol_edit_preview_clear(QgView *display,
	const char *previewId);

#endif /* QGOBOL_EDIT_PREVIEW_H */

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
