/*                  Q G S K E T C H E D I T . C P P
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
/** @file QgSketchEdit.cpp */

#include "common.h"

#include <array>
#include <vector>

#include "bn/mat.h"
#include "bu/str.h"
#include "ged.h"
#include "raytrace.h"
#include "rt/db_fullpath.h"
#include "rt/edit.h"
#include "rt/geom.h"
#include "rt/primitives/sketch.h"
#include "wdb.h"
#include "qtcad/QgSketchEdit.h"


class QgSketchEdit::Private
{
    public:
	struct Point {
	    std::array<fastf_t, 3> coordinates;
	};

	struct ged *gedp = nullptr;
	ged_edit_session_ref session = GED_EDIT_SESSION_REF_NULL;
	const struct rt_edit_prim_desc *descriptor = nullptr;
	QString path;
	QString error;
	mat_t pathMatrix = MAT_INIT_IDN;
	mat_t inversePathMatrix = MAT_INIT_IDN;
	point_t origin = VINIT_ZERO;
	vect_t uAxis = VINIT_ZERO;
	vect_t vAxis = VINIT_ZERO;
	fastf_t base2local = 1.0;
	bool valid = false;
	std::vector<Point> vertices;
	std::vector<Point> linePoints;
	std::vector<int> lineCommands;
	std::vector<int> lineSegments;
	std::vector<int> segmentTypes;
};


static bool
sketch_uv_from_offset(const vect_t offset, const vect_t uAxis,
	const vect_t vAxis, fastf_t *u, fastf_t *v)
{
    const fastf_t uu = VDOT(uAxis, uAxis);
    const fastf_t uv = VDOT(uAxis, vAxis);
    const fastf_t vv = VDOT(vAxis, vAxis);
    const fastf_t determinant = uu * vv - uv * uv;
    if (determinant <= SMALL_FASTF)
	return false;
    const fastf_t offsetU = VDOT(offset, uAxis);
    const fastf_t offsetV = VDOT(offset, vAxis);
    *u = (offsetU * vv - offsetV * uv) / determinant;
    *v = (offsetV * uu - offsetU * uv) / determinant;
    return true;
}


QgSketchEdit::QgSketchEdit() : d(new Private)
{
}


QgSketchEdit::~QgSketchEdit()
{
    delete d;
}


void
QgSketchEdit::clear()
{
    d->gedp = nullptr;
    d->session = GED_EDIT_SESSION_REF_NULL;
    d->descriptor = nullptr;
    d->path.clear();
    d->error.clear();
    MAT_IDN(d->pathMatrix);
    MAT_IDN(d->inversePathMatrix);
    VSETALL(d->origin, 0.0);
    VSETALL(d->uAxis, 0.0);
    VSETALL(d->vAxis, 0.0);
    d->base2local = 1.0;
    d->valid = false;
    d->vertices.clear();
    d->linePoints.clear();
    d->lineCommands.clear();
    d->lineSegments.clear();
    d->segmentTypes.clear();
}


bool
QgSketchEdit::refresh(struct ged *gedp, ged_edit_session_ref session,
	const QString &path)
{
    clear();
    if (!gedp || !gedp->dbip || path.trimmed().isEmpty() ||
	ged_edit_session_ref_is_null(session)) {
	d->error = QStringLiteral("No active sketch edit session");
	return false;
    }

    const struct rt_edit_prim_desc *descriptor = nullptr;
    if (ged_edit_session_descriptor_get(gedp, session, &descriptor) !=
	GED_EDIT_OK || !descriptor || !descriptor->prim_type ||
	!BU_STR_EQUAL(descriptor->prim_type, "sketch")) {
	d->error = QStringLiteral("The active edit session is not a sketch");
	return false;
    }

    struct db_full_path fullPath;
    db_full_path_init(&fullPath);
    const QByteArray pathBytes = path.toUtf8();
    const bool pathOk = db_string_to_path(&fullPath, gedp->dbip,
	pathBytes.constData()) == 0 && fullPath.fp_len > 0 &&
	db_path_to_mat(gedp->dbip, &fullPath, d->pathMatrix,
	    static_cast<int>(fullPath.fp_len) - 1);
    db_free_full_path(&fullPath);
    if (!pathOk) {
	d->error = QStringLiteral("Unable to resolve sketch occurrence path");
	return false;
    }
    bn_mat_inv(d->inversePathMatrix, d->pathMatrix);

    struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
    if (ged_edit_session_internal_copy(gedp, session, &intern) !=
	GED_EDIT_OK || intern.idb_type != ID_SKETCH || !intern.idb_ptr) {
	rt_db_free_internal(&intern);
	d->error = QStringLiteral("Unable to snapshot the active sketch");
	return false;
    }
    const struct rt_sketch_internal *sketch =
	static_cast<const struct rt_sketch_internal *>(intern.idb_ptr);
    RT_SKETCH_CK_MAGIC(sketch);

    struct rt_wdb *wdbp = wdb_dbopen(gedp->dbip, RT_WDB_TYPE_DB_DEFAULT);
    struct rt_sketch_edit_geometry geometry = {};
    if (!wdbp || rt_sketch_edit_geometry_get(&geometry, &intern,
	    &wdbp->wdb_ttol) != BRLCAD_OK) {
	rt_sketch_edit_geometry_free(&geometry);
	rt_db_free_internal(&intern);
	d->error = QStringLiteral("Unable to realize retained sketch topology");
	return false;
    }

    d->vertices.reserve(geometry.vertex_count);
    for (size_t i = 0; i < geometry.vertex_count; i++) {
	Private::Point point;
	for (int axis = 0; axis < 3; axis++)
	    point.coordinates[axis] = geometry.vertices[i][axis];
	d->vertices.push_back(point);
    }
    d->linePoints.reserve(geometry.line_count);
    d->lineCommands.reserve(geometry.line_count);
    d->lineSegments.reserve(geometry.line_count);
    for (size_t i = 0; i < geometry.line_count; i++) {
	Private::Point point;
	for (int axis = 0; axis < 3; axis++)
	    point.coordinates[axis] = geometry.line_points[i][axis];
	d->linePoints.push_back(point);
	d->lineCommands.push_back(geometry.line_commands[i]);
	d->lineSegments.push_back(geometry.line_segments[i]);
    }
    if (geometry.segment_count)
	d->segmentTypes.assign(geometry.segment_types,
	    geometry.segment_types + geometry.segment_count);

    VMOVE(d->origin, sketch->V);
    VMOVE(d->uAxis, sketch->u_vec);
    VMOVE(d->vAxis, sketch->v_vec);
    d->gedp = gedp;
    d->session = session;
    d->descriptor = descriptor;
    d->path = path;
    d->base2local = gedp->dbip->dbi_base2local;
    d->valid = true;
    d->error.clear();

    rt_sketch_edit_geometry_free(&geometry);
    rt_db_free_internal(&intern);
    return true;
}


bool
QgSketchEdit::isValid() const
{
    return d->valid;
}


QString
QgSketchEdit::errorString() const
{
    return d->error;
}


ged_edit_session_ref
QgSketchEdit::session() const
{
    return d->session;
}


int
QgSketchEdit::commandId(const char *name) const
{
    if (!d->descriptor || !name)
	return 0;
    for (int i = 0; i < d->descriptor->ncmd; i++) {
	if (d->descriptor->cmds[i].name &&
	    BU_STR_EQUAL(d->descriptor->cmds[i].name, name))
	    return d->descriptor->cmds[i].cmd_id;
    }
    return 0;
}


int
QgSketchEdit::selectCommand(FeatureDomain domain) const
{
    return commandId(domain == Vertex ? "ECMD_SKETCH_PICK_VERTEX" :
	"ECMD_SKETCH_PICK_SEGMENT");
}


int
QgSketchEdit::moveCommand(FeatureDomain domain) const
{
    return commandId(domain == Vertex ? "ECMD_SKETCH_MOVE_VERTEX" :
	"ECMD_SKETCH_MOVE_SEGMENT");
}


size_t
QgSketchEdit::vertexCount() const
{
    return d->vertices.size();
}


size_t
QgSketchEdit::segmentCount() const
{
    return d->segmentTypes.size();
}


size_t
QgSketchEdit::linePointCount() const
{
    return d->linePoints.size();
}


bool
QgSketchEdit::vertexDisplayPoint(size_t index, point_t point) const
{
    if (!d->valid || index >= d->vertices.size() || !point)
	return false;
    point_t source;
    for (int axis = 0; axis < 3; axis++)
	source[axis] = d->vertices[index].coordinates[axis];
    MAT4X3PNT(point, d->pathMatrix, source);
    return true;
}


bool
QgSketchEdit::lineDisplayPoint(size_t index, point_t point, int *command,
	int *segment) const
{
    if (!d->valid || index >= d->linePoints.size() || !point)
	return false;
    point_t source;
    for (int axis = 0; axis < 3; axis++)
	source[axis] = d->linePoints[index].coordinates[axis];
    MAT4X3PNT(point, d->pathMatrix, source);
    if (command)
	*command = d->lineCommands[index];
    if (segment)
	*segment = d->lineSegments[index];
    return true;
}


int
QgSketchEdit::segmentType(size_t index) const
{
    return index < d->segmentTypes.size() ? d->segmentTypes[index] : -1;
}


void
QgSketchEdit::pathMatrix(mat_t matrix) const
{
    if (matrix)
	MAT_COPY(matrix, d->pathMatrix);
}


bool
QgSketchEdit::displayPointToUV(const point_t point, fastf_t *u,
	fastf_t *v) const
{
    if (!d->valid || !point || !u || !v)
	return false;
    point_t source;
    vect_t offset;
    MAT4X3PNT(source, d->inversePathMatrix, point);
    VSUB2(offset, source, d->origin);
    return sketch_uv_from_offset(offset, d->uAxis, d->vAxis, u, v);
}


bool
QgSketchEdit::displayDeltaToUV(const vect_t delta, fastf_t *u,
	fastf_t *v) const
{
    if (!d->valid || !delta || !u || !v)
	return false;
    vect_t sourceDelta;
    MAT4X3VEC(sourceDelta, d->inversePathMatrix, delta);
    return sketch_uv_from_offset(sourceDelta, d->uAxis, d->vAxis, u, v);
}


enum ged_edit_status
QgSketchEdit::applyCommand(int commandId, const fastf_t *values,
	size_t valueCount, struct ged_view_context *view)
{
    if (!d->valid || !d->gedp || !commandId)
	return GED_EDIT_UNSUPPORTED;
    struct ged_edit_command_input input = {};
    input.command_id = commandId;
    input.values = values;
    input.value_count = valueCount;
    input.view = view;
    return ged_edit_session_apply(d->gedp, d->session, &input);
}


enum ged_edit_status
QgSketchEdit::selectFeature(FeatureDomain domain, int index,
	struct ged_view_context *view)
{
    const fastf_t value = static_cast<fastf_t>(index);
    return applyCommand(selectCommand(domain), &value, 1, view);
}


enum ged_edit_status
QgSketchEdit::moveVertexToDisplayPoint(const point_t point,
	struct ged_view_context *view)
{
    fastf_t values[2];
    if (!displayPointToUV(point, &values[0], &values[1]))
	return GED_EDIT_REJECTED;
    values[0] *= d->base2local;
    values[1] *= d->base2local;
    return applyCommand(moveCommand(Vertex), values, 2, view);
}


enum ged_edit_status
QgSketchEdit::moveSegmentByDisplayDelta(const vect_t delta,
	struct ged_view_context *view)
{
    fastf_t values[2];
    if (!displayDeltaToUV(delta, &values[0], &values[1]))
	return GED_EDIT_REJECTED;
    values[0] *= d->base2local;
    values[1] *= d->base2local;
    return applyCommand(moveCommand(Segment), values, 2, view);
}


enum ged_edit_status
QgSketchEdit::addVertexAtDisplayPoint(const point_t point,
	struct ged_view_context *view)
{
    fastf_t values[2];
    if (!displayPointToUV(point, &values[0], &values[1]))
	return GED_EDIT_REJECTED;
    values[0] *= d->base2local;
    values[1] *= d->base2local;
    return applyCommand(commandId("ECMD_SKETCH_ADD_VERTEX"), values, 2,
	view);
}


enum ged_edit_status
QgSketchEdit::setArcRadius(fastf_t radiusBaseUnits,
	struct ged_view_context *view)
{
    const fastf_t value = radiusBaseUnits * d->base2local;
    return applyCommand(commandId("ECMD_SKETCH_SET_ARC_RADIUS"), &value,
	1, view);
}


enum ged_edit_status
QgSketchEdit::setTangency(int adjacentSegment, fastf_t angleRadians,
	struct ged_view_context *view)
{
    const fastf_t values[2] = {
	static_cast<fastf_t>(adjacentSegment), angleRadians
    };
    return applyCommand(commandId("ECMD_SKETCH_SET_TANGENCY"), values, 2,
	view);
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
