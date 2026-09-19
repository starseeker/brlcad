/*                    Q G S K E T C H E D I T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
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
/** @file qtcad/QgSketchEdit.h
 *
 * Session-backed sketch interaction support shared by Qt widgets and retained
 * scene manipulators.  GED owns all mutable edit state; this class holds only
 * a refreshable presentation snapshot and instance-coordinate transforms.
 */

#ifndef QGSKETCHEDIT_H
#define QGSKETCHEDIT_H

#include <stddef.h>

#include <QString>

#include "ged/edit.h"
#include "qtcad/defines.h"
#include "vmath.h"

struct ged;
struct ged_view_context;

class QTCAD_EXPORT QgSketchEdit
{
    public:
	enum FeatureDomain {
	    Vertex = 0,
	    Segment = 1
	};

	QgSketchEdit();
	~QgSketchEdit();

	QgSketchEdit(const QgSketchEdit &) = delete;
	QgSketchEdit &operator=(const QgSketchEdit &) = delete;

	/**
	 * Refresh the retained topology and coordinate mapping for @p session.
	 * @p path is the occurrence path edited by the session, not merely the
	 * terminal sketch name.
	 */
	bool refresh(struct ged *gedp, ged_edit_session_ref session,
	    const QString &path);
	void clear();

	bool isValid() const;
	QString errorString() const;
	ged_edit_session_ref session() const;

	int commandId(const char *name) const;
	int selectCommand(FeatureDomain domain) const;
	int moveCommand(FeatureDomain domain) const;

	size_t vertexCount() const;
	size_t segmentCount() const;
	size_t linePointCount() const;
	bool vertexDisplayPoint(size_t index, point_t point) const;
	bool lineDisplayPoint(size_t index, point_t point, int *command,
	    int *segment) const;
	int segmentType(size_t index) const;
	void pathMatrix(mat_t matrix) const;

	/** Convert a displayed model-space point to sketch UV in base units. */
	bool displayPointToUV(const point_t point, fastf_t *u,
	    fastf_t *v) const;

	/** Convert a display-space delta to sketch UV in base units. */
	bool displayDeltaToUV(const vect_t delta, fastf_t *u,
	    fastf_t *v) const;

	/** Apply an arbitrary descriptor command to the authoritative session. */
	enum ged_edit_status applyCommand(int commandId, const fastf_t *values,
	    size_t valueCount, struct ged_view_context *view = nullptr);

	enum ged_edit_status selectFeature(FeatureDomain domain, int index,
	    struct ged_view_context *view = nullptr);
	enum ged_edit_status moveVertexToDisplayPoint(const point_t point,
	    struct ged_view_context *view = nullptr);
	enum ged_edit_status moveSegmentByDisplayDelta(const vect_t delta,
	    struct ged_view_context *view = nullptr);
	enum ged_edit_status addVertexAtDisplayPoint(const point_t point,
	    struct ged_view_context *view = nullptr);
	enum ged_edit_status setArcRadius(fastf_t radiusBaseUnits,
	    struct ged_view_context *view = nullptr);
	enum ged_edit_status setTangency(int adjacentSegment,
	    fastf_t angleRadians = 0.0,
	    struct ged_view_context *view = nullptr);

    private:
	class Private;
	Private *d;
};

#endif /* QGSKETCHEDIT_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
