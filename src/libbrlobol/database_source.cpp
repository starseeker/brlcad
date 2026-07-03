/*                D A T A B A S E _ S O U R C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol/database_source.h"
#include "brlobol/lod_mesh_shape.h"
#include "brlobol/lod_realization.h"
#include "brlobol/material_object.h"
#include "brlobol/mesh_shape.h"
#include "brlobol/vlist_shape.h"

#include "bu/color.h"
#include "bu/list.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "nmg.h"
#include "raytrace.h"
#include "rt/func.h"
#include "rt/global.h"
#include "rt/nongeom.h"
#include "rt/db_fullpath.h"
#include "rt/primitives/annot.h"
#include "rt/tree.h"
#include "rt/vlist.h"
#include "rt/view.h"

#include <Inventor/SbName.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/sensors/SoFieldSensor.h>

#include <limits.h>
#include <map>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>

SO_NODE_SOURCE(SoBRLDatabaseSource);

static int
database_source_float_different(float a, float b)
{
    return fabsf(a - b) > 1.0e-6f;
}

static int
database_source_color_equal(const SbColor &a, const SbColor &b)
{
    return !database_source_float_different(a[0], b[0]) &&
	!database_source_float_different(a[1], b[1]) &&
	!database_source_float_different(a[2], b[2]);
}

static int
database_source_vec3f_equal(const SbVec3f &a, const SbVec3f &b)
{
    return !database_source_float_different(a[0], b[0]) &&
	!database_source_float_different(a[1], b[1]) &&
	!database_source_float_different(a[2], b[2]);
}

static int
database_source_string_equal(const SbString &a, const char *b)
{
    return strcmp(a.getString(), b ? b : "") == 0;
}

template <typename FieldT>
static void
database_source_assign_string(FieldT &field, const char *value)
{
    const char *nextValue = value ? value : "";
    if (database_source_string_equal(field.getValue(), nextValue))
	return;

    const std::string stableValue(nextValue);
    field = stableValue.c_str();
}

template <typename FieldT>
static void
database_source_assign_string(FieldT &field, const SbString &value)
{
    database_source_assign_string(field, value.getString());
}

static SbBox3f
database_source_box_from_minmax(const SbVec3f &bmin, const SbVec3f &bmax)
{
    SbBox3f bounds;
    bounds.makeEmpty();
    bounds.extendBy(bmin);
    bounds.extendBy(bmax);
    return bounds;
}

static SbBox3f
database_source_transform_bounds(const SbBox3f &bounds,
	const SbMatrix &matrix)
{
    SbBox3f transformed;
    transformed.makeEmpty();
    if (bounds.isEmpty())
	return transformed;

    const SbVec3f bmin = bounds.getMin();
    const SbVec3f bmax = bounds.getMax();
    for (int xi = 0; xi < 2; xi++) {
	for (int yi = 0; yi < 2; yi++) {
	    for (int zi = 0; zi < 2; zi++) {
		const SbVec3f corner(
			xi ? bmax[0] : bmin[0],
			yi ? bmax[1] : bmin[1],
			zi ? bmax[2] : bmin[2]);
		SbVec3f transformedCorner;
		matrix.multVecMatrix(corner, transformedCorner);
		transformed.extendBy(transformedCorner);
	    }
	}
    }

    return transformed;
}

BRLObolDatabaseSourceSummary::BRLObolDatabaseSourceSummary(void) :
    valid(FALSE),
    path(""),
    instanceKey(""),
    displayName(""),
    hasParent(FALSE),
    drawTreeDepth(0),
    parentGroupPath(""),
    representationKey(""),
    representationMode(-1),
    drawMode(SoBRLDatabaseSource::WIREFRAME),
    sourceRevision(0),
    inputsRevision(0),
    viewRevision(0),
    realizedRevision(0),
    realizedSourceRevision(0),
    realizedInputsRevision(0),
    realizedViewRevision(0),
    realizationStatus(SoBRLDatabaseSource::UNREALIZED),
    realizationDiagnostic(""),
    realizationIdentity(""),
    realizationRoleFlags(SoBRLDatabaseSource::REALIZATION_ROLE_NONE),
    realizationViewDependent(FALSE),
    realizationViewScale(0.0f),
    realizationBotThreshold(0),
    realizationCurveScale(0.0f),
    realizationPointScale(0.0f),
    visible(TRUE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0),
    materialPolicy(SoBRLDatabaseSource::MATERIAL_INHERIT),
    colorOverride(FALSE),
    color(1.0f, 1.0f, 1.0f),
    drawMatrixValid(FALSE),
    drawMatrix(SbMatrix::identity()),
    drawCenterValid(FALSE),
    drawCenter(0.0f, 0.0f, 0.0f),
    drawSizeValid(FALSE),
    drawSize(0.0f),
    sourceBoundsValid(FALSE),
    sourceBounds(),
    stale(TRUE),
    staleReason(SoBRLDatabaseSource::STALE_SOURCE),
    realizedShapeCount(0),
    realizedMeshCount(0),
    realizedMaterialObjectCount(0)
{
    this->sourceBounds.makeEmpty();
}

BRLObolAuxiliaryLineSetDisplayState::BRLObolAuxiliaryLineSetDisplayState(void) :
    valid(FALSE),
    drawMode(BRLOBOL_LOD_DRAW_WIRE),
    visible(TRUE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0f),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0)
{
}

BRLObolExternalLineSet::BRLObolExternalLineSet(void) :
    points(NULL),
    commands(NULL),
    precisePoints(NULL),
    count(0),
    sourceType("line-set"),
    geometryKind("line")
{
}

BRLObolExternalPointSet::BRLObolExternalPointSet(void) :
    points(NULL),
    precisePoints(NULL),
    count(0),
    sourceType("point-set"),
    geometryKind("point")
{
}

BRLObolExternalTriangleMesh::BRLObolExternalTriangleMesh(void) :
    points(NULL),
    pointCount(0),
    indices(NULL),
    indexCount(0),
    sourceType("indexed-face-set"),
    geometryKind("surface")
{
}

BRLObolExternalAnnotationSegment::BRLObolExternalAnnotationSegment(void) :
    kind(SEGMENT_NONE),
    lineStart(0),
    lineEnd(0),
    textRefPoint(0),
    text(NULL)
{
}

BRLObolExternalAnnotation::BRLObolExternalAnnotation(void) :
    basePoint(0.0f, 0.0f, 0.0f),
    linePoints(NULL),
    lineCommands(NULL),
    preciseLinePoints(NULL),
    linePointCount(0),
    annotationPoints(NULL),
    preciseAnnotationPoints(NULL),
    annotationPointCount(0),
    segments(NULL),
    segmentCount(0),
    sourceType("annotation"),
    geometryKind("annotation")
{
}

BRLObolDatabaseSourceDisplayPatch::BRLObolDatabaseSourceDisplayPatch(void) :
    visibleValid(FALSE),
    visible(TRUE),
    highlightedValid(FALSE),
    highlighted(FALSE),
    lineStyleValid(FALSE),
    lineStyle(0),
    lineWidthValid(FALSE),
    lineWidth(0),
    transparencyValid(FALSE),
    transparency(0.0f),
    colorOverrideValid(FALSE),
    colorOverride(FALSE),
    colorValid(FALSE),
    color(1.0f, 1.0f, 1.0f)
{
}

BRLObolRealizedShapeSummary::BRLObolRealizedShapeSummary(void) :
    valid(FALSE),
    shapeKind(SHAPE_UNKNOWN),
    path(""),
    sourceName(""),
    sourceType(""),
    sourceId(0),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    ownerDrawMode(0),
    ownerSourceRevision(0),
    ownerInputsRevision(0),
    ownerViewRevision(0),
    ownerRealizedRevision(0),
    ownerRealizedSourceRevision(0),
    ownerRealizedInputsRevision(0),
    ownerRealizedViewRevision(0),
    ownerRealizationStatus(SoBRLDatabaseSource::UNREALIZED),
    ownerRealizationDiagnostic(""),
    ownerRealizationIdentity(""),
    ownerSourceStale(FALSE),
    ownerStaleReason(SoBRLDatabaseSource::STALE_NONE),
    displayName(""),
    geometryName(""),
    cacheIdentity(""),
    sourceIdentity(""),
    databaseIntent(FALSE),
    overlayIntent(FALSE),
    hudIntent(FALSE),
    localSource(FALSE),
    sharedSource(FALSE),
    nonDatabaseSource(FALSE),
    drawMode(0),
    recordRole(""),
    geometryKind(""),
    regionId(0),
    airCode(0),
    materialId(0),
    los(0),
    materialColorValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialShader(""),
    materialRevision(0),
    visible(FALSE),
    selectable(FALSE),
    selected(FALSE),
    highlighted(FALSE),
    ghosted(FALSE),
    hiddenLine(FALSE),
    editEmphasis(FALSE),
    editIntentId(""),
    editIntentRole(""),
    lodPolicy(0),
    pointCount(0),
    commandCount(0),
    segmentCount(0),
    pointPrimitiveCount(0),
    triangleCount(0),
    indexCount(0),
    boundsValid(FALSE),
    bounds()
{
    this->bounds.makeEmpty();
}

BRLObolRealizedMaterialSummary::BRLObolRealizedMaterialSummary(void) :
    valid(FALSE),
    sourcePath(""),
    sourceName(""),
    sourceType(""),
    sourceId(0),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    ownerDrawMode(0),
    ownerSourceRevision(0),
    ownerInputsRevision(0),
    ownerViewRevision(0),
    ownerRealizedRevision(0),
    ownerRealizedSourceRevision(0),
    ownerRealizedInputsRevision(0),
    ownerRealizedViewRevision(0),
    ownerRealizationStatus(SoBRLDatabaseSource::UNREALIZED),
    ownerRealizationDiagnostic(""),
    ownerRealizationIdentity(""),
    ownerSourceStale(FALSE),
    ownerStaleReason(SoBRLDatabaseSource::STALE_NONE),
    materialName(""),
    parentName(""),
    materialSource(""),
    propertyCount(0)
{
}

BRLObolSceneTreeSummary::BRLObolSceneTreeSummary(void) :
    valid(FALSE),
    nodeKind(NODE_UNKNOWN),
    isGroup(FALSE),
    isShape(FALSE),
    isDatabaseSource(FALSE),
    isMaterialObject(FALSE),
    hasParent(FALSE),
    drawTreeDepth(0),
    childCount(0),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path(""),
    sourceName(""),
    sourceType(""),
    sourceId(0),
    displayName(""),
    geometryName("")
{
}

BRLObolSceneDisplaySummary::BRLObolSceneDisplaySummary(void) :
    valid(FALSE),
    nodeKind(BRLObolSceneTreeSummary::NODE_UNKNOWN),
    isDatabaseSource(FALSE),
    hasDrawIntent(FALSE),
    intentPath(""),
    intentDrawMode(-1),
    visible(TRUE),
    highlighted(FALSE),
    lineStyle(0),
    lineWidth(0),
    transparency(0.0),
    drawMode(0),
    materialValid(FALSE),
    materialColor(1.0f, 1.0f, 1.0f),
    materialRevision(0),
    drawMatrixValid(FALSE),
    drawMatrix(SbMatrix::identity()),
    drawCenterValid(FALSE),
    drawCenter(0.0f, 0.0f, 0.0f),
    drawSizeValid(FALSE),
    drawSize(0.0f),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path("")
{
}

BRLObolSceneMaterialSummary::BRLObolSceneMaterialSummary(void) :
    valid(FALSE),
    nodeKind(BRLObolSceneTreeSummary::NODE_UNKNOWN),
    materialValid(FALSE),
    materialRevision(0),
    materialColor(1.0f, 1.0f, 1.0f),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path("")
{
}

BRLObolSceneBoundsSummary::BRLObolSceneBoundsSummary(void) :
    valid(FALSE),
    nodeKind(BRLObolSceneTreeSummary::NODE_UNKNOWN),
    boundsValid(FALSE),
    bounds(),
    ownerSourceIndex(-1),
    ownerSourcePath(""),
    ownerSourceInstanceKey(""),
    path("")
{
    this->bounds.makeEmpty();
}

static const char *
lookup_name_from_path(const SbString &path)
{
    const char *name = path.getString();
    while (*name == '/')
	name++;
    if (!name[0])
	return NULL;
    return name;
}

static const char *
database_source_skip_leading_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static SbString
record_identity_with_revision(const char *identity, uint32_t revision)
{
    SbString ret = identity ? identity : "";
    char revisionString[64] = {0};
    snprintf(revisionString, sizeof(revisionString), "#%u", revision);
    ret += revisionString;
    return ret;
}

static int
source_record_draw_mode(const SoBRLDatabaseSource *source)
{
    if (!source)
	return BRLOBOL_LOD_DRAW_WIRE;

    switch (source->representationMode.getValue()) {
	case SoBRLDatabaseSource::REPRESENTATION_SHADED_BOTS:
	    return BRLOBOL_LOD_DRAW_SHADED_BOTS;
	case SoBRLDatabaseSource::REPRESENTATION_SHADED:
	    return BRLOBOL_LOD_DRAW_SHADED;
	case SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE:
	    return BRLOBOL_LOD_DRAW_HIDDEN_LINE;
	case SoBRLDatabaseSource::REPRESENTATION_EVAL_POINTS:
	    return BRLOBOL_LOD_DRAW_POINTS;
	case SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE:
	case SoBRLDatabaseSource::REPRESENTATION_WIRE:
	    return BRLOBOL_LOD_DRAW_WIRE;
	default:
	    break;
    }

    if (source->drawMode.getValue() == SoBRLDatabaseSource::SHADED)
	return BRLOBOL_LOD_DRAW_SHADED;
    return BRLOBOL_LOD_DRAW_WIRE;
}

static SbString
source_effective_instance_key(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";
    const SbString key = source->instanceKey.getValue();
    if (key.getLength() > 0)
	return key;
    return source->path.getValue();
}

static SbString
source_effective_representation_key(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";
    const SbString key = source->representationKey.getValue();
    if (key.getLength() > 0)
	return key;
    return source_effective_instance_key(source);
}

static int
source_instance_matches_record_path(const SbString &instanceKey,
	const char *recordPath)
{
    if (!recordPath || !recordPath[0])
	return 0;
    if (instanceKey.getLength() == 0)
	return 0;

    const char *key = instanceKey.getString();
    if (strcmp(key, recordPath) == 0)
	return 1;
    return strcmp(key[0] == '/' ? key + 1 : key,
	    recordPath[0] == '/' ? recordPath + 1 : recordPath) == 0;
}

static SbString
source_record_identity(const SoBRLDatabaseSource *source,
	const char *recordPath)
{
    const char *path = recordPath ? recordPath : "";
    if (!source)
	return path;

    const SbString instanceKey = source_effective_instance_key(source);
    const SbString representationKey =
	source_effective_representation_key(source);
    if (source_instance_matches_record_path(instanceKey, path)) {
	if (representationKey.getLength() == 0 ||
		source_instance_matches_record_path(representationKey, path))
	    return path;
	SbString identity = representationKey;
	if (path[0]) {
	    identity += "::";
	    identity += path;
	}
	return identity;
    }
    if (instanceKey.getLength() == 0 && representationKey.getLength() == 0)
	return path;

    SbString identity = representationKey.getLength() > 0 ?
	representationKey : instanceKey;
    if (path[0]) {
	identity += "::";
	identity += path;
    }
    return identity;
}

static SbString
source_realization_identity(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";

    SbString identity;
    identity.sprintf("dbsource:%s@%s#repr=%s;repr_mode=%d;mode=%d;source=%u;inputs=%u;view=%u;lod=%u;abs=%.9g;rel=%.9g;norm=%.9g",
	    source_effective_instance_key(source).getString(),
	    source->path.getValue().getString(),
	    source_effective_representation_key(source).getString(),
	    source->representationMode.getValue(),
	    source->drawMode.getValue(),
	    source->sourceRevision.getValue(),
	    source->inputsRevision.getValue(),
	    source->viewRevision.getValue(),
	    source->lodBotThreshold.getValue(),
	    static_cast<double>(source->tessellationAbsTol.getValue()),
	    static_cast<double>(source->tessellationRelTol.getValue()),
	    static_cast<double>(source->tessellationNormTol.getValue()));
    return identity;
}

static void
convert_vlist(std::vector<SbVec3f> &points, std::vector<int32_t> &commands, const struct bu_list *vhead)
{
    rt_vlist *vp = NULL;

    BU_LIST_EACH(vhead, vp, rt_vlist) {
	for (size_t i = 0; i < vp->nused; i++) {
	    int cmd = -1;
	    switch (vp->cmd[i]) {
		case RT_VLIST_LINE_MOVE:
		case RT_VLIST_POLY_MOVE:
		case RT_VLIST_TRI_MOVE:
		    cmd = SoBRLVListShape::MOVE;
		    break;
		case RT_VLIST_LINE_DRAW:
		case RT_VLIST_POLY_DRAW:
		case RT_VLIST_POLY_END:
		case RT_VLIST_TRI_DRAW:
		case RT_VLIST_TRI_END:
		    cmd = SoBRLVListShape::DRAW;
		    break;
		case RT_VLIST_POINT_DRAW:
		    cmd = SoBRLVListShape::POINT;
		    break;
		default:
		    break;
	    }
	    if (cmd >= 0) {
		points.push_back(SbVec3f(static_cast<float>(vp->pt[i][0]),
			    static_cast<float>(vp->pt[i][1]),
			    static_cast<float>(vp->pt[i][2])));
		commands.push_back(cmd);
	    }
	}
    }
}

static double
nonnegative_or_default(float value, double defaultValue)
{
    return value >= 0.0f ? static_cast<double>(value) : defaultValue;
}

static struct bg_tess_tol
source_tess_tol(const SoBRLDatabaseSource *source)
{
    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_TOL;

    if (source) {
	ttol.abs = nonnegative_or_default(source->tessellationAbsTol.getValue(), ttol.abs);
	ttol.rel = nonnegative_or_default(source->tessellationRelTol.getValue(), ttol.rel);
	ttol.norm = nonnegative_or_default(source->tessellationNormTol.getValue(), ttol.norm);
    }

    return ttol;
}

static SbBool
node_is_auxiliary_vlist(const SoNode *node)
{
    if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	return FALSE;

    const SoBRLVListShape *shape = static_cast<const SoBRLVListShape *>(node);
    return strcmp(shape->recordRole.getValue().getString(), "auxiliary") == 0 ?
	TRUE : FALSE;
}

static SbBool
node_is_auxiliary_source(const SoNode *node)
{
    if (!node || !node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return FALSE;

    const SoBRLDatabaseSource *source =
	static_cast<const SoBRLDatabaseSource *>(node);
    return source->auxiliarySource.getValue();
}

static const char *
database_source_placement_transform_name(void)
{
    return "__brlobol_source_placement";
}

static SbBool
node_is_source_placement_transform(const SoNode *node)
{
    return node &&
	node->isOfType(SoMatrixTransform::getClassTypeId()) &&
	node->getName() == SbName(database_source_placement_transform_name()) ?
	TRUE : FALSE;
}

static SoMatrixTransform *
source_placement_transform(SoBRLDatabaseSource *source)
{
    if (!source)
	return NULL;

    for (int i = 0; i < source->getNumChildren(); i++) {
	SoNode *child = source->getChild(i);
	if (node_is_source_placement_transform(child))
	    return static_cast<SoMatrixTransform *>(child);
    }

    return NULL;
}

static const SoMatrixTransform *
source_placement_transform(const SoBRLDatabaseSource *source)
{
    return source_placement_transform(const_cast<SoBRLDatabaseSource *>(source));
}

static SoMatrixTransform *
ensure_source_placement_transform(SoBRLDatabaseSource *source)
{
    if (!source)
	return NULL;

    SoMatrixTransform *transform = source_placement_transform(source);
    if (!transform) {
	transform = new SoMatrixTransform;
	transform->setName(SbName(database_source_placement_transform_name()));
	source->insertChild(transform, 0);
	return transform;
    }

    int transformIndex = -1;
    for (int i = 0; i < source->getNumChildren(); i++) {
	if (source->getChild(i) == transform) {
	    transformIndex = i;
	    break;
	}
    }
    if (transformIndex > 0) {
	transform->ref();
	source->removeChild(transformIndex);
	source->insertChild(transform, 0);
	transform->unref();
    }

    return transform;
}

static int
remove_source_placement_transform(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;

    int removed = 0;
    for (int i = source->getNumChildren() - 1; i >= 0; i--) {
	if (node_is_source_placement_transform(source->getChild(i))) {
	    source->removeChild(i);
	    removed = 1;
	}
    }
    return removed;
}

static int
sync_source_placement_transform(SoBRLDatabaseSource *source)
{
    if (!source)
	return 0;

    if (!source->drawMatrixValid.getValue()) {
	return remove_source_placement_transform(source);
    }

    const int hadTransform = source_placement_transform(source) ? 1 : 0;
    SoMatrixTransform *transform = ensure_source_placement_transform(source);
    if (!transform)
	return 0;

    int changed = hadTransform ? 0 : 1;
    const SbMatrix matrix = source->drawMatrix.getValue();
    if (!transform->matrix.getValue().equals(matrix, 0.000001f)) {
	transform->matrix = matrix;
	changed = 1;
    }
    return changed;
}

static void
database_source_add_realized_child(SoBRLDatabaseSource *source,
	SoNode *child)
{
    if (!source || !child)
	return;

    (void)sync_source_placement_transform(source);
    source->addChild(child);
}

static void
remove_non_auxiliary_children(SoGroup *group)
{
    if (!group)
	return;

    for (int i = group->getNumChildren() - 1; i >= 0; i--) {
	SoNode *child = group->getChild(i);
	if (node_is_source_placement_transform(child) ||
		node_is_auxiliary_vlist(child) ||
		node_is_auxiliary_source(child))
	    continue;
	group->removeChild(i);
    }
}

struct realize_walk_data {
    SoBRLDatabaseSource *source;
    uint32_t revision;
    int realized_shapes;
    int failed_shapes;
    SbString diagnostic;
    std::map<std::string, SoBRLVListShape *> sharedWireGeometry;
    std::map<std::string, SoBRLVListShape *> sharedMeshVListGeometry;
    std::map<std::string, SoBRLMeshShape *> sharedMeshGeometry;
};

static SbMatrix
mat_to_sbmatrix(const mat_t mat)
{
    return SbMatrix(
	    static_cast<float>(mat[0]), static_cast<float>(mat[4]),
	    static_cast<float>(mat[8]), static_cast<float>(mat[12]),
	    static_cast<float>(mat[1]), static_cast<float>(mat[5]),
	    static_cast<float>(mat[9]), static_cast<float>(mat[13]),
	    static_cast<float>(mat[2]), static_cast<float>(mat[6]),
	    static_cast<float>(mat[10]), static_cast<float>(mat[14]),
	    static_cast<float>(mat[3]), static_cast<float>(mat[7]),
	    static_cast<float>(mat[11]), static_cast<float>(mat[15]));
}

static SoSeparator *
realize_instance_leaf_separator(const struct db_tree_state *tsp)
{
    SoSeparator *leaf = new SoSeparator;
    if (!tsp || bn_mat_is_identity(tsp->ts_mat))
	return leaf;

    SoMatrixTransform *transform = new SoMatrixTransform;
    transform->matrix = mat_to_sbmatrix(tsp->ts_mat);
    leaf->addChild(transform);
    return leaf;
}

static std::string
realize_geometry_cache_key(const struct directory *dp)
{
    return (dp && dp->d_namep) ? std::string(dp->d_namep) : std::string();
}

static void
set_walk_diagnostic(struct realize_walk_data *data,
	const struct db_full_path *pathp,
	const char *reason)
{
    if (!data || data->diagnostic.getLength() > 0)
	return;

    char *path = pathp ? db_path_to_string(pathp) : NULL;
    SbString msg;
    msg.sprintf("%s: %s", path ? path : "<unknown>", reason ? reason : "realization failed");
    data->diagnostic = msg;
    if (path)
	bu_free(path, "db_path_to_string");
}

static const char *
primitive_type_label(const struct rt_db_internal *intern)
{
    if (!intern || !intern->idb_meth || !intern->idb_meth->ft_label[0])
	return "unknown";
    return intern->idb_meth->ft_label;
}

static bool
primitive_is_annotation(int internalType, const char *typeLabel)
{
    if (internalType == ID_ANNOT)
	return true;
    return typeLabel && (BU_STR_EQUAL(typeLabel, "annot") ||
	    BU_STR_EQUAL(typeLabel, "annotation"));
}

static uint32_t
internal_payload_magic(const struct rt_db_internal *intern)
{
    if (!intern || !intern->idb_ptr)
	return 0;
    return *((const uint32_t *)intern->idb_ptr);
}

static bool
internal_payload_magic_valid(const struct rt_db_internal *intern)
{
    if (!intern || !intern->idb_ptr || !intern->idb_meth)
	return false;
    if (intern->idb_meth->magic != RT_FUNCTAB_MAGIC)
	return false;

    const uint32_t expected = intern->idb_meth->ft_internal_magic;
    if (!expected)
	return true;

    return internal_payload_magic(intern) == expected;
}

struct valid_walk_internal {
    struct rt_db_internal local;
    struct rt_db_internal *intern;
    bool ownsLocal;
    bool refetched;

    valid_walk_internal(void) : intern(NULL), ownsLocal(false), refetched(false)
    {
	RT_DB_INTERNAL_INIT(&local);
    }

    ~valid_walk_internal(void)
    {
	if (ownsLocal)
	    rt_db_free_internal(&local);
    }
};

static struct rt_db_internal *
fetch_local_walk_internal(struct db_tree_state *tsp,
	const struct db_full_path *pathp,
	struct valid_walk_internal *handle)
{
    if (!handle || !tsp || !tsp->ts_dbip || !pathp)
	return NULL;

    struct directory *dp = DB_FULL_PATH_CUR_DIR(pathp);
    if (!dp)
	return NULL;

    handle->refetched = true;
    if (rt_db_get_internal(&handle->local, dp, tsp->ts_dbip, NULL) >= 0 &&
	    internal_payload_magic_valid(&handle->local)) {
	handle->intern = &handle->local;
	handle->ownsLocal = true;
	return handle->intern;
    }

    return NULL;
}

static void
set_invalid_internal_diagnostic(struct realize_walk_data *data,
	const struct db_full_path *pathp,
	const struct rt_db_internal *intern,
	bool refetched)
{
    if (!data)
	return;

    const char *typeLabel = primitive_type_label(intern);
    char reason[256] = {0};
    if (intern && intern->idb_ptr && intern->idb_meth &&
	    intern->idb_meth->ft_internal_magic) {
	snprintf(reason, sizeof(reason),
		"invalid primitive payload for type '%s' (magic 0x%08x, expected 0x%08x)%s",
		typeLabel,
		(unsigned int)internal_payload_magic(intern),
		(unsigned int)intern->idb_meth->ft_internal_magic,
		refetched ? "; refetch failed validation" : "");
    } else {
	snprintf(reason, sizeof(reason),
		"invalid primitive payload for type '%s'%s",
		typeLabel,
		refetched ? "; refetch failed validation" : "");
    }
    set_walk_diagnostic(data, pathp, reason);
}

static void
assign_annotation_record(SoBRLVListShape *shape,
	const struct rt_annot_internal *annot)
{
    if (!shape || !annot)
	return;

    RT_ANNOT_CK_MAGIC(annot);

    shape->sourceType = "annotation";
    shape->geometryKind = "annotation";
    shape->annotationBasePoint = SbVec3f(
	    static_cast<float>(annot->V[X]),
	    static_cast<float>(annot->V[Y]),
	    static_cast<float>(annot->V[Z]));

    const int pointCount = (annot->vert_count > static_cast<size_t>(INT_MAX)) ?
	INT_MAX : static_cast<int>(annot->vert_count);
    if (pointCount > 0 && annot->verts) {
	std::vector<double> points(static_cast<size_t>(pointCount) * 3);
	for (int i = 0; i < pointCount; i++) {
	    const size_t offset = static_cast<size_t>(i) * 3;
	    points[offset + 0] = annot->verts[i][X];
	    points[offset + 1] = annot->verts[i][Y];
	    points[offset + 2] = 0.0;
	}
	shape->setPreciseAnnotationPoints(points.data(), pointCount);
    } else {
	shape->setPreciseAnnotationPoints(NULL, 0);
    }

    const int segmentCount = (annot->ant.count > static_cast<size_t>(INT_MAX)) ?
	INT_MAX : static_cast<int>(annot->ant.count);
    shape->annotationSegmentTextValid.setNum(segmentCount);
    shape->annotationSegmentKind.setNum(segmentCount);
    shape->annotationSegmentStart.setNum(segmentCount);
    shape->annotationSegmentEnd.setNum(segmentCount);
    shape->annotationTextRefPoint.setNum(segmentCount);
    shape->annotationText.setNum(segmentCount);

    for (int i = 0; i < segmentCount; i++) {
	shape->annotationSegmentTextValid.set1Value(i, FALSE);
	shape->annotationSegmentKind.set1Value(i,
		SoBRLVListShape::ANNOTATION_SEGMENT_NONE);
	shape->annotationSegmentStart.set1Value(i, 0);
	shape->annotationSegmentEnd.set1Value(i, 0);
	shape->annotationTextRefPoint.set1Value(i, 0);
	shape->annotationText.set1Value(i, "");

	const uint32_t *magic = annot->ant.segments ?
	    static_cast<const uint32_t *>(annot->ant.segments[i]) : NULL;
	if (!magic)
	    continue;

	switch (*magic) {
	    case CURVE_LSEG_MAGIC: {
		const struct line_seg *lsg =
		    reinterpret_cast<const struct line_seg *>(magic);
		shape->annotationSegmentKind.set1Value(i,
			SoBRLVListShape::ANNOTATION_SEGMENT_LINE);
		shape->annotationSegmentStart.set1Value(i, lsg->start);
		shape->annotationSegmentEnd.set1Value(i, lsg->end);
		break;
	    }
	    case ANN_TSEG_MAGIC: {
		const struct txt_seg *tsg =
		    reinterpret_cast<const struct txt_seg *>(magic);
		const char *label = BU_VLS_IS_INITIALIZED(&tsg->label) ?
		    bu_vls_cstr(&tsg->label) : "";
		shape->annotationSegmentKind.set1Value(i,
			SoBRLVListShape::ANNOTATION_SEGMENT_TEXT);
		shape->annotationTextRefPoint.set1Value(i, tsg->ref_pt);
		shape->annotationText.set1Value(i, label ? label : "");
		shape->annotationSegmentTextValid.set1Value(i,
			(label && label[0]) ? TRUE : FALSE);
		break;
	    }
	    default:
		break;
	}
    }
}

static void
free_annotation_record_copy(struct rt_annot_internal *annot)
{
    if (!annot)
	return;

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_ANNOT;
    intern.idb_meth = &OBJ[ID_ANNOT];
    intern.idb_ptr = annot;
    rt_db_free_internal(&intern);
}

template <typename ShapeT>
static void
assign_realized_identity(ShapeT *shape,
	const struct db_tree_state *tsp,
	const char *path,
	const char *sourceName,
	const char *sourceType,
	uint32_t sourceId,
	const SoBRLDatabaseSource *source)
{
    if (!shape)
	return;

    shape->sourcePath = path ? path : "";
    shape->sourceName = sourceName ? sourceName : "";
    shape->sourceType = sourceType ? sourceType : "";
    shape->sourceId = sourceId;
    if (source && source->displayName.getValue().getLength() > 0)
	shape->displayName = source->displayName.getValue();
    else
	shape->displayName = sourceName ? sourceName : (path ? path : "");
    shape->geometryName = sourceName ? sourceName : "";
    shape->sourceIdentity = source_record_identity(source, path);
    shape->cacheIdentity = record_identity_with_revision(
	    shape->sourceIdentity.getValue().getString(), sourceId);
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = FALSE;
    if (shape->databaseIntent.getValue() ||
	    !shape->nonDatabaseSource.getValue())
	shape->drawMode = source_record_draw_mode(source);
    shape->hiddenLine = (source_record_draw_mode(source) ==
	    BRLOBOL_LOD_DRAW_HIDDEN_LINE) ? TRUE : FALSE;
    shape->recordRole = "database";
    shape->geometryKind = "";
    if (source) {
	shape->visible = source->visible.getValue();
	shape->highlighted = source->highlighted.getValue();
	shape->lineStyle = source->lineStyle.getValue();
	shape->lineWidth = source->lineWidth.getValue();
	shape->transparency = source->transparency.getValue();
	shape->colorOverride = source->colorOverride.getValue();
	shape->color = source->color.getValue();
	shape->materialColorValid = source->materialColorValid.getValue();
	shape->materialColor = source->materialColor.getValue();
	shape->materialRevision = source->materialRevision.getValue();
    }

    if (!tsp) {
	shape->regionId = 0;
	shape->airCode = 0;
	shape->materialId = 0;
	shape->los = 0;
	if (!source || !source->materialColorValid.getValue()) {
	    shape->materialColorValid = FALSE;
	    shape->materialColor = SbColor(1.0f, 1.0f, 1.0f);
	}
	shape->materialShader = "";
	return;
    }

    shape->regionId = tsp->ts_regionid;
    shape->airCode = tsp->ts_aircode;
    shape->materialId = tsp->ts_gmater;
    shape->los = tsp->ts_los;
    shape->materialColorValid = tsp->ts_mater.ma_color_valid ? TRUE : FALSE;
    shape->materialColor = SbColor(
	    static_cast<float>(tsp->ts_mater.ma_color[0]),
	    static_cast<float>(tsp->ts_mater.ma_color[1]),
	    static_cast<float>(tsp->ts_mater.ma_color[2]));
    if (source && source->materialColorValid.getValue()) {
	shape->materialColorValid = TRUE;
	shape->materialColor = source->materialColor.getValue();
    }
    shape->materialShader = tsp->ts_mater.ma_shader ? tsp->ts_mater.ma_shader : "";
}

template <typename ShapeT>
static void
assign_shared_geometry_identity(ShapeT *shape,
	const char *sourceName,
	const char *sourceType,
	uint32_t sourceId,
	const char *geometryKind)
{
    if (!shape)
	return;

    const char *name = sourceName ? sourceName : "";
    shape->sourcePath = name;
    shape->sourceName = name;
    shape->sourceType = sourceType ? sourceType : "";
    shape->sourceId = sourceId;
    shape->displayName = name;
    shape->geometryName = name;
    shape->sourceIdentity = name;
    shape->cacheIdentity = record_identity_with_revision(name, sourceId);
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = TRUE;
    shape->nonDatabaseSource = FALSE;
    shape->recordRole = "shared-geometry";
    if (geometryKind && geometryKind[0])
	shape->geometryKind = geometryKind;
}

static void
retarget_realized_node_source(SoNode *node,
	const SoBRLDatabaseSource *source,
	const char *oldSourcePath,
	const char *newSourcePath,
	uint32_t revision);

template <typename ShapeT>
static void
sync_shape_owner_state(ShapeT *shape, const SoBRLDatabaseSource *source)
{
    if (!shape || !source)
	return;

    database_source_assign_string(shape->ownerSourcePath,
	    source->path.getValue());
    database_source_assign_string(shape->ownerSourceInstanceKey,
	    source_effective_instance_key(source));
    shape->ownerSourceRevision = source->sourceRevision.getValue();
    shape->ownerInputsRevision = source->inputsRevision.getValue();
    shape->ownerViewRevision = source->viewRevision.getValue();
    shape->ownerRealizedRevision = source->realizedRevision.getValue();
    shape->ownerRealizedSourceRevision =
	source->realizedSourceRevision.getValue();
    shape->ownerRealizedInputsRevision =
	source->realizedInputsRevision.getValue();
    shape->ownerRealizedViewRevision = source->realizedViewRevision.getValue();
    shape->ownerRealizationStatus = source->realizationStatus.getValue();
    database_source_assign_string(shape->ownerRealizationDiagnostic,
	    source->realizationDiagnostic.getValue());
    database_source_assign_string(shape->ownerRealizationIdentity,
	    source->realizationIdentity.getValue());
    shape->ownerSourceStale = source->stale.getValue();
    shape->ownerStaleReason = source->staleReason.getValue();
}

template <typename ShapeT>
static void
sync_shape_placement_state(ShapeT *shape, const SoBRLDatabaseSource *source)
{
    if (!shape || !source)
	return;

    shape->drawMatrixValid = source->drawMatrixValid.getValue();
    shape->drawMatrix = source->drawMatrix.getValue();
    shape->drawCenterValid = source->drawCenterValid.getValue();
    shape->drawCenter = source->drawCenter.getValue();
    shape->drawSizeValid = source->drawSizeValid.getValue();
    shape->drawSize = source->drawSize.getValue();
}

template <typename ShapeT>
static void
sync_shape_display_state(ShapeT *shape, const SoBRLDatabaseSource *source)
{
    if (!shape || !source)
	return;

    if (!shape->databaseIntent.getValue() &&
	    shape->nonDatabaseSource.getValue())
	return;

    shape->drawMode = source_record_draw_mode(source);
    shape->hiddenLine = (shape->drawMode.getValue() ==
	    BRLOBOL_LOD_DRAW_HIDDEN_LINE) ? TRUE : FALSE;
    shape->visible = source->visible.getValue();
    shape->highlighted = source->highlighted.getValue();
    shape->lineStyle = source->lineStyle.getValue();
    shape->lineWidth = source->lineWidth.getValue();
    shape->transparency = source->transparency.getValue();
    shape->colorOverride = source->colorOverride.getValue();
    shape->color = source->color.getValue();
    if (source->materialColorValid.getValue()) {
	shape->materialColorValid = TRUE;
	shape->materialColor = source->materialColor.getValue();
	shape->materialRevision = source->materialRevision.getValue();
    }
}

template <typename ShapeT>
static void
sync_shape_display_name(ShapeT *shape, const SoBRLDatabaseSource *source)
{
    if (!shape || !source)
	return;

    if (!shape->databaseIntent.getValue() &&
	    shape->nonDatabaseSource.getValue())
	return;

    std::string nextDisplayName;
    if (source->displayName.getValue().getLength() > 0) {
	nextDisplayName = source->displayName.getValue().getString();
    } else if (shape->sourceName.getValue().getLength() > 0) {
	nextDisplayName = shape->sourceName.getValue().getString();
    } else {
	nextDisplayName = shape->sourcePath.getValue().getString();
    }
    database_source_assign_string(shape->displayName,
	    nextDisplayName.c_str());
}

static void
sync_realized_shape_owner_state_in_node(SoNode *node,
	const SoBRLDatabaseSource *source)
{
    if (!node || !source)
	return;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	sync_shape_owner_state(shape, source);
	sync_shape_placement_state(shape, source);
	sync_shape_display_state(shape, source);
	sync_shape_display_name(shape, source);
	return;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);
	sync_shape_owner_state(shape, source);
	sync_shape_placement_state(shape, source);
	sync_shape_display_state(shape, source);
	sync_shape_display_name(shape, source);
	return;
    }
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    sync_realized_shape_owner_state_in_node(group->getChild(i), source);
    }
}

static union tree *
make_nop_tree(void)
{
    union tree *tp = NULL;
    BU_GET(tp, union tree);
    RT_TREE_INIT(tp);
    tp->tr_op = OP_NOP;
    return tp;
}

static void
material_object_add_properties(SoBRLMaterialObject *object,
	const char *group,
	const struct bu_attribute_value_set *properties)
{
    const struct bu_attribute_value_pair *avpp = NULL;

    if (!object || !properties)
	return;

    for (BU_AVS_FOR(avpp, properties))
	object->addProperty(group, avpp->name, avpp->value);
}

static SoBRLMaterialObject *
material_object_from_internal(struct rt_material_internal *material)
{
    if (!material)
	return NULL;
    RT_CHECK_MATERIAL(material);

    SoBRLMaterialObject *object = new SoBRLMaterialObject;
    object->materialName = bu_vls_cstr(&material->name);
    object->parentName = bu_vls_cstr(&material->parent);
    object->materialSource = bu_vls_cstr(&material->source);
    material_object_add_properties(object, "physical",
	    &material->physicalProperties);
    material_object_add_properties(object, "mechanical",
	    &material->mechanicalProperties);
    material_object_add_properties(object, "optical",
	    &material->opticalProperties);
    material_object_add_properties(object, "thermal",
	    &material->thermalProperties);
    return object;
}

static void
assign_material_identity(SoBRLMaterialObject *object,
	const char *path,
	const char *sourceName,
	const char *sourceType,
	uint32_t sourceId)
{
    if (!object)
	return;

    object->sourcePath = path ? path : "";
    object->sourceName = sourceName ? sourceName : "";
    object->sourceType = sourceType ? sourceType : "";
    object->sourceId = sourceId;
}

static SoBRLVListShape *
vlist_from_plot_internal(struct rt_db_internal *intern,
	const SoBRLDatabaseSource *source)
{
    if (!internal_payload_magic_valid(intern))
	return NULL;

    point_t annotationBase = VINIT_ZERO;
    struct rt_annot_internal *annotation = NULL;
    const char *typeLabel = primitive_type_label(intern);
    const bool addAnnotationBase = primitive_is_annotation(intern->idb_type,
	    typeLabel);
    if (addAnnotationBase) {
	struct rt_annot_internal *annot =
	    static_cast<struct rt_annot_internal *>(intern->idb_ptr);
	RT_ANNOT_CK_MAGIC(annot);
	annotation = rt_copy_annot(annot);
	VMOVE(annotationBase, annot->V);
    }

    struct bu_list vhead;
    BU_LIST_INIT(&vhead);
    struct bg_tess_tol ttol = source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    int ret = rt_obj_plot(&vhead, intern, &ttol, &tol);
    if (ret < 0) {
	free_annotation_record_copy(annotation);
	RT_FREE_VLIST(&rt_vlfree, &vhead);
	return NULL;
    }

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    convert_vlist(points, commands, &vhead);
    RT_FREE_VLIST(&rt_vlfree, &vhead);
    if (addAnnotationBase) {
	for (size_t i = 0; i < points.size(); i++) {
	    points[i].setValue(points[i][0] +
		    static_cast<float>(annotationBase[X]),
		    points[i][1] + static_cast<float>(annotationBase[Y]),
		    points[i][2] + static_cast<float>(annotationBase[Z]));
	}
    }
    if (points.empty() || points.size() != commands.size()) {
	free_annotation_record_copy(annotation);
	return NULL;
    }

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->setLineSet(points.data(), commands.data(),
	    static_cast<int>(points.size()));
    assign_annotation_record(shape, annotation);
    free_annotation_record_copy(annotation);
    return shape;
}

static union tree *
realize_leaf(struct db_tree_state *tsp,
	const struct db_full_path *pathp,
	struct directory *dp,
	void *client_data)
{
    struct realize_walk_data *data = static_cast<struct realize_walk_data *>(client_data);
    if (!data || !data->source || !pathp || !tsp || !tsp->ts_dbip ||
	    !dp)
	return TREE_NULL;

    SoBRLVListShape *sharedShape = NULL;
    const std::string cacheKey = realize_geometry_cache_key(dp);
    std::map<std::string, SoBRLVListShape *>::iterator found =
	data->sharedWireGeometry.find(cacheKey);
    if (found != data->sharedWireGeometry.end())
	sharedShape = found->second;

    const char *typeLabel = sharedShape ?
	sharedShape->sourceType.getValue().getString() : NULL;
    if (!sharedShape) {
	valid_walk_internal validInternal;
	struct rt_db_internal *localIntern =
	    fetch_local_walk_internal(tsp, pathp, &validInternal);
	if (!localIntern) {
	    data->failed_shapes++;
	    set_invalid_internal_diagnostic(data, pathp, NULL,
		    validInternal.refetched);
	    return TREE_NULL;
	}

	typeLabel = primitive_type_label(localIntern);
	if (localIntern->idb_type == ID_MATERIAL) {
	    SoBRLMaterialObject *materialObject =
		material_object_from_internal(static_cast<struct rt_material_internal *>(localIntern->idb_ptr));
	    if (!materialObject) {
		data->failed_shapes++;
		set_walk_diagnostic(data, pathp,
			"material object realization failed");
		return TREE_NULL;
	    }
	    char *path = db_path_to_string(pathp);
	    SoSeparator *leaf = new SoSeparator;
	    assign_material_identity(materialObject, path, dp->d_namep,
		    typeLabel, data->revision);
	    leaf->addChild(materialObject);
	    data->source->addChild(leaf);
	    data->realized_shapes++;
	    if (path)
		bu_free(path, "db_path_to_string");
	    return make_nop_tree();
	}

	sharedShape = vlist_from_plot_internal(localIntern, data->source);
	if (!sharedShape) {
	    char reason[256] = {0};
	    data->failed_shapes++;
	    snprintf(reason, sizeof(reason),
		    "wireframe plot produced no usable geometry for primitive type '%s'",
		    typeLabel);
	    set_walk_diagnostic(data, pathp, reason);
	    return TREE_NULL;
	}
	assign_shared_geometry_identity(sharedShape, dp->d_namep, typeLabel,
		data->revision, "line");
	if (primitive_is_annotation(localIntern->idb_type, typeLabel)) {
	    sharedShape->sourceType = "annotation";
	    sharedShape->geometryKind = "annotation";
	}
	data->sharedWireGeometry[cacheKey] = sharedShape;
	typeLabel = sharedShape->sourceType.getValue().getString();
    }

    char *path = db_path_to_string(pathp);
    SoSeparator *leaf = realize_instance_leaf_separator(tsp);
    SoBRLVListShape *shape = new SoBRLVListShape;
    assign_realized_identity(shape, tsp, path, dp->d_namep, typeLabel,
	    data->revision, data->source);
    shape->setSharedGeometry(sharedShape);
    const char *geometryKind = sharedShape->geometryKind.getValue().getString();
    shape->geometryKind = (geometryKind && geometryKind[0]) ?
	geometryKind : "line";
    if (geometryKind && BU_STR_EQUAL(geometryKind, "annotation"))
	shape->sourceType = "annotation";
    leaf->addChild(shape);
    data->source->addChild(leaf);
    data->realized_shapes++;
    if (path)
	bu_free(path, "db_path_to_string");

    return make_nop_tree();
}

static SoBRLVListShape *
vlist_from_pnts(const struct rt_pnts_internal *pnts)
{
    if (!pnts || !pnts->point || pnts->count == 0)
	return NULL;
    RT_PNTS_CK_MAGIC(pnts);

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    std::vector<int> colorValid;
    std::vector<SbColor> colors;
    std::vector<int> scaleValid;
    std::vector<float> scales;
    std::vector<int> normalValid;
    std::vector<SbVec3f> normals;
    points.reserve(pnts->count);
    commands.reserve(pnts->count);

    colorValid.reserve(pnts->count);
    colors.reserve(pnts->count);
    scaleValid.reserve(pnts->count);
    scales.reserve(pnts->count);
    normalValid.reserve(pnts->count);
    normals.reserve(pnts->count);

    const double defaultScale = pnts->scale;
    auto appendPoint = [&](const fastf_t *v, const struct bu_color *c,
	    const fastf_t *s, const fastf_t *n) {
	points.push_back(SbVec3f(static_cast<float>(v[X]),
		static_cast<float>(v[Y]), static_cast<float>(v[Z])));
	commands.push_back(SoBRLVListShape::POINT);
	if (c) {
	    colorValid.push_back(1);
	    colors.push_back(SbColor(
		    static_cast<float>(c->buc_rgb[RED]),
		    static_cast<float>(c->buc_rgb[GRN]),
		    static_cast<float>(c->buc_rgb[BLU])));
	} else {
	    colorValid.push_back(0);
	    colors.push_back(SbColor(1.0f, 1.0f, 1.0f));
	}
	if (s && *s > 0.0) {
	    scaleValid.push_back(1);
	    scales.push_back(static_cast<float>(*s));
	} else if (!s && defaultScale > 0.0) {
	    scaleValid.push_back(1);
	    scales.push_back(static_cast<float>(defaultScale));
	} else {
	    scaleValid.push_back(0);
	    scales.push_back(0.0f);
	}
	if (n) {
	    normalValid.push_back(1);
	    normals.push_back(SbVec3f(static_cast<float>(n[X]),
		    static_cast<float>(n[Y]), static_cast<float>(n[Z])));
	} else {
	    normalValid.push_back(0);
	    normals.push_back(SbVec3f(0.0f, 0.0f, 1.0f));
	}
    };

    switch (pnts->type) {
	case RT_PNT_TYPE_PNT:
	    {
		const struct pnt *point = NULL;
		for (BU_LIST_FOR(point, pnt, &(((struct pnt *)pnts->point)->l)))
		    appendPoint(point->v, NULL, NULL, NULL);
	    }
	    break;
	case RT_PNT_TYPE_COL:
	    {
		const struct pnt_color *point = NULL;
		for (BU_LIST_FOR(point, pnt_color, &(((struct pnt_color *)pnts->point)->l)))
		    appendPoint(point->v, &point->c, NULL, NULL);
	    }
	    break;
	case RT_PNT_TYPE_SCA:
	    {
		const struct pnt_scale *point = NULL;
		for (BU_LIST_FOR(point, pnt_scale, &(((struct pnt_scale *)pnts->point)->l)))
		    appendPoint(point->v, NULL, &point->s, NULL);
	    }
	    break;
	case RT_PNT_TYPE_NRM:
	    {
		const struct pnt_normal *point = NULL;
		for (BU_LIST_FOR(point, pnt_normal, &(((struct pnt_normal *)pnts->point)->l)))
		    appendPoint(point->v, NULL, NULL, point->n);
	    }
	    break;
	case RT_PNT_TYPE_COL_SCA:
	    {
		const struct pnt_color_scale *point = NULL;
		for (BU_LIST_FOR(point, pnt_color_scale, &(((struct pnt_color_scale *)pnts->point)->l)))
		    appendPoint(point->v, &point->c, &point->s, NULL);
	    }
	    break;
	case RT_PNT_TYPE_COL_NRM:
	    {
		const struct pnt_color_normal *point = NULL;
		for (BU_LIST_FOR(point, pnt_color_normal, &(((struct pnt_color_normal *)pnts->point)->l)))
		    appendPoint(point->v, &point->c, NULL, point->n);
	    }
	    break;
	case RT_PNT_TYPE_SCA_NRM:
	    {
		const struct pnt_scale_normal *point = NULL;
		for (BU_LIST_FOR(point, pnt_scale_normal, &(((struct pnt_scale_normal *)pnts->point)->l)))
		    appendPoint(point->v, NULL, &point->s, point->n);
	    }
	    break;
	case RT_PNT_TYPE_COL_SCA_NRM:
	    {
		const struct pnt_color_scale_normal *point = NULL;
		for (BU_LIST_FOR(point, pnt_color_scale_normal, &(((struct pnt_color_scale_normal *)pnts->point)->l)))
		    appendPoint(point->v, &point->c, &point->s, point->n);
	    }
	    break;
	default:
	    return NULL;
    }

    if (points.empty() || points.size() != commands.size())
	return NULL;

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->setLineSet(points.data(), commands.data(), static_cast<int>(points.size()));
    shape->setPointAttributes(colorValid.data(), colors.data(),
	    scaleValid.data(), scales.data(),
	    normalValid.data(), normals.data(),
	    static_cast<int>(points.size()));
    return shape;
}

static SoBRLMeshShape *
mesh_from_arb(const struct rt_arb_internal *arb)
{
    if (!arb)
	return NULL;
    RT_ARB_CK_MAGIC(arb);

    SbVec3f points[8];
    for (int i = 0; i < 8; i++) {
	points[i] = SbVec3f(static_cast<float>(arb->pt[i][X]),
		static_cast<float>(arb->pt[i][Y]),
		static_cast<float>(arb->pt[i][Z]));
    }

    static const int32_t indices[36] = {
	3, 2, 1, 3, 1, 0,
	4, 5, 6, 4, 6, 7,
	4, 7, 3, 4, 3, 0,
	2, 6, 5, 2, 5, 1,
	1, 5, 4, 1, 4, 0,
	7, 6, 2, 7, 2, 3
    };

    SoBRLMeshShape *shape = new SoBRLMeshShape;
    shape->setIndexedTriangles(points, 8, indices, 36);
    return shape;
}

static SoBRLMeshShape *
mesh_from_bot(const struct rt_bot_internal *bot,
	const SoBRLDatabaseSource *source)
{
    if (!bot || !bot->vertices || !bot->faces ||
	    bot->num_vertices == 0 || bot->num_faces == 0 ||
	    bot->num_vertices > INT_MAX || bot->num_faces > INT_MAX / 3)
	return NULL;
    RT_BOT_CK_MAGIC(bot);

    std::vector<SbVec3f> points;
    points.reserve(bot->num_vertices);
    for (size_t i = 0; i < bot->num_vertices; i++) {
	points.push_back(SbVec3f(static_cast<float>(bot->vertices[i * 3]),
		    static_cast<float>(bot->vertices[i * 3 + 1]),
		    static_cast<float>(bot->vertices[i * 3 + 2])));
    }

    std::vector<int32_t> indices;
    indices.reserve(bot->num_faces * 3);
    for (size_t i = 0; i < bot->num_faces * 3; i++)
	indices.push_back(static_cast<int32_t>(bot->faces[i]));

    uint32_t threshold = source ? source->lodBotThreshold.getValue() : 0;
    SoBRLMeshShape *shape = (threshold > 0 &&
	    bot->num_faces >= static_cast<size_t>(threshold)) ?
	new SoBRLLodMeshShape : new SoBRLMeshShape;
    shape->setIndexedTriangles(points.data(), static_cast<int>(points.size()),
	    indices.data(), static_cast<int>(indices.size()));
    return shape;
}

static void
primitive_indexed_face_set_free(struct rt_primitive_indexed_face_set *faceSet)
{
    if (!faceSet)
	return;
    if (faceSet->points)
	bu_free(faceSet->points, "primitive indexed-face points");
    if (faceSet->normals)
	bu_free(faceSet->normals, "primitive indexed-face normals");
    if (faceSet->indices)
	bu_free(faceSet->indices, "primitive indexed-face indices");
    memset(faceSet, 0, sizeof(*faceSet));
}

static int
indexed_face_finish(std::vector<int32_t> &face,
	size_t pointCount,
	std::vector<int32_t> &triangles,
	size_t *faceCount,
	unsigned int *faceStamp,
	std::vector<unsigned int> &seen)
{
    if (face.size() < 3)
	return 0;

    for (size_t i = 1; i + 1 < face.size(); i++) {
	triangles.push_back(face[0]);
	triangles.push_back(face[i]);
	triangles.push_back(face[i + 1]);
    }

    face.clear();
    if (faceCount)
	(*faceCount)++;
    if (faceStamp && seen.size() == pointCount) {
	if (*faceStamp == UINT_MAX) {
	    for (size_t i = 0; i < seen.size(); i++)
		seen[i] = 0;
	    *faceStamp = 1;
	} else {
	    (*faceStamp)++;
	}
    }
    return 1;
}

static int
indexed_faces_to_triangles(const int *indices,
	size_t indexCount,
	size_t pointCount,
	std::vector<int32_t> &triangles)
{
    if (!indices || !indexCount || !pointCount ||
	    pointCount > static_cast<size_t>(INT_MAX) ||
	    indexCount > static_cast<size_t>(INT_MAX))
	return 0;

    size_t faceCount = 0;
    unsigned int faceStamp = 1;
    std::vector<unsigned int> seen(pointCount, 0);
    std::vector<int32_t> face;

    for (size_t i = 0; i < indexCount; i++) {
	const int idx = indices[i];
	if (idx < 0) {
	    if (idx != -1 || !indexed_face_finish(face, pointCount,
		    triangles, &faceCount, &faceStamp, seen))
		return 0;
	    continue;
	}

	if (static_cast<size_t>(idx) >= pointCount)
	    return 0;
	if (seen[static_cast<size_t>(idx)] == faceStamp)
	    return 0;
	seen[static_cast<size_t>(idx)] = faceStamp;
	face.push_back(static_cast<int32_t>(idx));
    }

    if (!face.empty() && !indexed_face_finish(face, pointCount,
	    triangles, &faceCount, &faceStamp, seen))
	return 0;
    return faceCount > 0 && !triangles.empty();
}

static SoBRLMeshShape *
mesh_from_indexed_face_set(const struct rt_primitive_indexed_face_set *faceSet,
	const SoBRLDatabaseSource *source)
{
    if (!faceSet || !faceSet->points || !faceSet->point_count ||
	    !faceSet->indices || !faceSet->index_count ||
	    faceSet->point_count > static_cast<size_t>(INT_MAX))
	return NULL;

    std::vector<int32_t> triangles;
    if (!indexed_faces_to_triangles(faceSet->indices, faceSet->index_count,
	    faceSet->point_count, triangles))
	return NULL;
    if (triangles.size() > static_cast<size_t>(INT_MAX))
	return NULL;

    std::vector<SbVec3f> points;
    points.reserve(faceSet->point_count);
    for (size_t i = 0; i < faceSet->point_count; i++) {
	points.push_back(SbVec3f(
		static_cast<float>(faceSet->points[i][X]),
		static_cast<float>(faceSet->points[i][Y]),
		static_cast<float>(faceSet->points[i][Z])));
    }

    uint32_t threshold = source ? source->lodBotThreshold.getValue() : 0;
    SoBRLMeshShape *shape = (threshold > 0 &&
	    triangles.size() / 3 >= static_cast<size_t>(threshold)) ?
	new SoBRLLodMeshShape : new SoBRLMeshShape;
    shape->setIndexedTriangles(points.data(), static_cast<int>(points.size()),
	    triangles.data(), static_cast<int>(triangles.size()));
    return shape;
}

static SoBRLMeshShape *
mesh_from_primitive_face_set(struct rt_db_internal *intern,
	const SoBRLDatabaseSource *source)
{
    if (!internal_payload_magic_valid(intern) || !intern->idb_meth->ft_indexed_face_set)
	return NULL;

    struct rt_primitive_indexed_face_set faceSet;
    struct rt_view_info viewInfo = RT_VIEW_INFO_INIT;
    struct bg_tess_tol ttol = source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    memset(&faceSet, 0, sizeof(faceSet));

    int ret = intern->idb_meth->ft_indexed_face_set(&faceSet, intern,
	    &ttol, &tol, &viewInfo);
    if (ret != BRLCAD_OK) {
	primitive_indexed_face_set_free(&faceSet);
	return NULL;
    }

    SoBRLMeshShape *shape = mesh_from_indexed_face_set(&faceSet, source);
    primitive_indexed_face_set_free(&faceSet);
    return shape;
}

static void
publish_lod_result_metadata(SoBRLMeshShape *shape,
	const BRLObolLodResult &result)
{
    if (!shape || result.providerStatus != BRLOBOL_LOD_PROVIDER_READY)
	return;

    (void)shape->applyStagedLodResult(result, &result.request);
}

static void
publish_lod_metadata_if_cached(SoBRLMeshShape *shape, struct db_i *dbip,
	const char *sourceName)
{
    if (!shape || !dbip || !sourceName)
	return;

    struct rt_mesh_lod_cache_status status = RT_MESH_LOD_CACHE_STATUS_INIT;
    (void)db_mesh_lod_status(dbip, sourceName, &status);

    struct rt_mesh_lod *lod = db_mesh_lod_get(dbip, sourceName);
    if (!lod)
	return;

    struct rt_mesh_lod_info info = RT_MESH_LOD_INFO_INIT;
    if (rt_mesh_lod_load_view(lod, NULL, 0) >= 0 &&
	    rt_mesh_lod_info_get(lod, &info)) {
	BRLObolLodRequest request;
	request.databaseId = dbip->dbi_filename ? dbip->dbi_filename : "";
	request.sourceContentHash = status.cache_key;
	request.objectPath = sourceName;
	request.objectName = sourceName;
	request.drawMode = BRLOBOL_LOD_DRAW_SHADED;
	request.lodPolicy = shape->lodPolicy.getValue();
	request.providerId = "rt_mesh_lod";
	request.providerVersion = "rt-cache-v1";
	request.qualityTier = BRLOBOL_LOD_QUALITY_FAST_DISPLAY;
	request.sourceCounts.faceCount = info.face_count;
	request.sourceCounts.pointCount = info.point_count;
	request.sourceCounts.originalPointCount = info.point_orig_count;
	request.sourceCounts.normalCount = info.normal_count;
	request.bounds = SbBox3f(
		SbVec3f(static_cast<float>(info.bmin[X]),
		    static_cast<float>(info.bmin[Y]),
		    static_cast<float>(info.bmin[Z])),
		SbVec3f(static_cast<float>(info.bmax[X]),
		    static_cast<float>(info.bmax[Y]),
		    static_cast<float>(info.bmax[Z])));
	BRLObolLodResult result =
	    brlobol_lod_result_from_rt_mesh_info(request, info, &status);
	publish_lod_result_metadata(shape, result);
    }

    rt_mesh_lod_destroy(lod);
}

static void
append_nmg_triangle_points(std::vector<SbVec3f> &points,
	std::vector<int32_t> &indices,
	const struct vertex *v0,
	const struct vertex *v1,
	const struct vertex *v2)
{
    if (!v0 || !v1 || !v2 || !v0->vg_p || !v1->vg_p || !v2->vg_p)
	return;

    int32_t base = static_cast<int32_t>(points.size());
    points.push_back(SbVec3f(static_cast<float>(v0->vg_p->coord[X]),
		static_cast<float>(v0->vg_p->coord[Y]),
		static_cast<float>(v0->vg_p->coord[Z])));
    points.push_back(SbVec3f(static_cast<float>(v1->vg_p->coord[X]),
		static_cast<float>(v1->vg_p->coord[Y]),
		static_cast<float>(v1->vg_p->coord[Z])));
    points.push_back(SbVec3f(static_cast<float>(v2->vg_p->coord[X]),
		static_cast<float>(v2->vg_p->coord[Y]),
		static_cast<float>(v2->vg_p->coord[Z])));
    indices.push_back(base);
    indices.push_back(base + 1);
    indices.push_back(base + 2);
}

static SoBRLMeshShape *
mesh_from_nmg_region(struct nmgregion *r, struct bu_list *vlfree,
	const struct bn_tol *tol)
{
    std::vector<SbVec3f> points;
    std::vector<int32_t> indices;

    if (!r || !r->m_p || !vlfree || !tol)
	return NULL;

    nmg_triangulate_model(r->m_p, vlfree, tol);

    for (struct shell *s = BU_LIST_FIRST(shell, &r->s_hd);
	    BU_LIST_NOT_HEAD(s, &r->s_hd);
	    s = BU_LIST_PNEXT(shell, &s->l)) {
	NMG_CK_SHELL(s);
	for (struct faceuse *fu = BU_LIST_FIRST(faceuse, &s->fu_hd);
		BU_LIST_NOT_HEAD(fu, &s->fu_hd);
		fu = BU_LIST_PNEXT(faceuse, &fu->l)) {
	    NMG_CK_FACEUSE(fu);
	    if (fu->orientation != OT_SAME)
		continue;
	    for (struct loopuse *lu = BU_LIST_FIRST(loopuse, &fu->lu_hd);
		    BU_LIST_NOT_HEAD(lu, &fu->lu_hd);
		    lu = BU_LIST_PNEXT(loopuse, &lu->l)) {
		struct vertex *verts[3] = {NULL, NULL, NULL};
		int vert_count = 0;

		NMG_CK_LOOPUSE(lu);
		if (BU_LIST_FIRST_MAGIC(&lu->down_hd) != NMG_EDGEUSE_MAGIC)
		    continue;

		for (struct edgeuse *eu = BU_LIST_FIRST(edgeuse, &lu->down_hd);
			BU_LIST_NOT_HEAD(eu, &lu->down_hd);
			eu = BU_LIST_PNEXT(edgeuse, &eu->l)) {
		    NMG_CK_EDGEUSE(eu);
		    if (vert_count < 3)
			verts[vert_count] = eu->vu_p ? eu->vu_p->v_p : NULL;
		    vert_count++;
		}
		if (vert_count == 3)
		    append_nmg_triangle_points(points, indices, verts[0], verts[1], verts[2]);
	    }
	}
    }

    if (points.empty() || indices.empty())
	return NULL;

    SoBRLMeshShape *shape = new SoBRLMeshShape;
    shape->setIndexedTriangles(points.data(), static_cast<int>(points.size()),
	    indices.data(), static_cast<int>(indices.size()));
    return shape;
}

static SoBRLMeshShape *
mesh_from_tessellated_internal(struct rt_db_internal *intern,
	const SoBRLDatabaseSource *source)
{
    if (!internal_payload_magic_valid(intern))
	return NULL;

    struct bg_tess_tol ttol = source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct bu_list vlfree;
    struct model *m = NULL;
    struct nmgregion *r = NULL;
    SoBRLMeshShape *shape = NULL;
    int ret = -1;

    BU_LIST_INIT(&vlfree);
    m = nmg_mm();
    if (!m)
	return NULL;

    if (!BU_SETJUMP) {
	ret = rt_obj_tess(&r, m, intern, &ttol, &tol);
    } else {
	BU_UNSETJUMP;
	nmg_km(m);
	bu_list_free(&vlfree);
	return NULL;
    } BU_UNSETJUMP;

    if (ret == 0 && r)
	shape = mesh_from_nmg_region(r, &vlfree, &tol);

    nmg_km(m);
    bu_list_free(&vlfree);
    return shape;
}

static SoBRLMeshShape *
mesh_from_internal(struct rt_db_internal *intern,
	const SoBRLDatabaseSource *source)
{
    if (!internal_payload_magic_valid(intern))
	return NULL;

    SoBRLMeshShape *faceSetShape = mesh_from_primitive_face_set(intern,
	    source);
    if (faceSetShape)
	return faceSetShape;

    switch (intern->idb_type) {
	case ID_ARB8:
	    return mesh_from_arb(static_cast<const struct rt_arb_internal *>(intern->idb_ptr));
	case ID_BOT:
	    return mesh_from_bot(static_cast<const struct rt_bot_internal *>(intern->idb_ptr), source);
	default:
	    break;
    }

    return mesh_from_tessellated_internal(intern, source);
}

static SoBRLMeshShape *
mesh_instance_for_shared_geometry(const SoBRLMeshShape *sharedShape)
{
    if (sharedShape &&
	    sharedShape->isOfType(SoBRLLodMeshShape::getClassTypeId()))
	return new SoBRLLodMeshShape;
    return new SoBRLMeshShape;
}

static union tree *
realize_mesh_leaf(struct db_tree_state *tsp,
	const struct db_full_path *pathp,
	struct directory *dp,
	void *client_data)
{
    struct realize_walk_data *data = static_cast<struct realize_walk_data *>(client_data);
    if (!data || !data->source || !pathp || !tsp || !tsp->ts_dbip ||
	    !dp)
	return TREE_NULL;

    SoBRLVListShape *sharedVListShape = NULL;
    SoBRLMeshShape *sharedMeshShape = NULL;
    const std::string cacheKey = realize_geometry_cache_key(dp);
    std::map<std::string, SoBRLVListShape *>::iterator foundVList =
	data->sharedMeshVListGeometry.find(cacheKey);
    if (foundVList != data->sharedMeshVListGeometry.end()) {
	sharedVListShape = foundVList->second;
    } else {
	std::map<std::string, SoBRLMeshShape *>::iterator foundMesh =
	    data->sharedMeshGeometry.find(cacheKey);
	if (foundMesh != data->sharedMeshGeometry.end())
	    sharedMeshShape = foundMesh->second;
    }

    const char *typeLabel = sharedVListShape ?
	sharedVListShape->sourceType.getValue().getString() :
	(sharedMeshShape ? sharedMeshShape->sourceType.getValue().getString() : NULL);
    if (!sharedVListShape && !sharedMeshShape) {
	valid_walk_internal validInternal;
	struct rt_db_internal *localIntern =
	    fetch_local_walk_internal(tsp, pathp, &validInternal);
	if (!localIntern) {
	    data->failed_shapes++;
	    set_invalid_internal_diagnostic(data, pathp, NULL,
		    validInternal.refetched);
	    return TREE_NULL;
	}

	typeLabel = primitive_type_label(localIntern);
	const int internalType = localIntern->idb_type;
	if (internalType == ID_MATERIAL) {
	    SoBRLMaterialObject *materialObject =
		material_object_from_internal(static_cast<struct rt_material_internal *>(localIntern->idb_ptr));
	    if (!materialObject) {
		data->failed_shapes++;
		set_walk_diagnostic(data, pathp,
			"material object realization failed");
		return TREE_NULL;
	    }
	    char *path = db_path_to_string(pathp);
	    SoSeparator *leaf = new SoSeparator;
	    assign_material_identity(materialObject, path, dp->d_namep,
		    typeLabel, data->revision);
	    leaf->addChild(materialObject);
	    data->source->addChild(leaf);
	    data->realized_shapes++;
	    if (path)
		bu_free(path, "db_path_to_string");
	    return make_nop_tree();
	}

	if (internalType == ID_PNTS) {
	    sharedVListShape = vlist_from_pnts(
		    static_cast<const struct rt_pnts_internal *>(localIntern->idb_ptr));
	    if (sharedVListShape)
		assign_shared_geometry_identity(sharedVListShape,
			dp->d_namep, typeLabel, data->revision, "point");
	} else if (internalType == ID_SKETCH || internalType == ID_ANNOT) {
	    sharedVListShape = vlist_from_plot_internal(localIntern,
		    data->source);
	    if (sharedVListShape) {
		assign_shared_geometry_identity(sharedVListShape,
			dp->d_namep, typeLabel, data->revision, "line");
		if (primitive_is_annotation(internalType, typeLabel)) {
		    sharedVListShape->sourceType = "annotation";
		    sharedVListShape->geometryKind = "annotation";
		}
	    }
	} else {
	    sharedMeshShape = mesh_from_internal(localIntern, data->source);
	    if (sharedMeshShape)
		assign_shared_geometry_identity(sharedMeshShape,
			dp->d_namep, typeLabel, data->revision, "surface");
	}

	if (!sharedVListShape && !sharedMeshShape) {
	    char reason[256] = {0};
	    data->failed_shapes++;
	    snprintf(reason, sizeof(reason),
		    "unsupported or failed mesh conversion/tessellation for primitive type '%s'",
		    typeLabel);
	    set_walk_diagnostic(data, pathp, reason);
	    return TREE_NULL;
	}

	if (sharedVListShape) {
	    data->sharedMeshVListGeometry[cacheKey] = sharedVListShape;
	    typeLabel = sharedVListShape->sourceType.getValue().getString();
	} else {
	    data->sharedMeshGeometry[cacheKey] = sharedMeshShape;
	    typeLabel = sharedMeshShape->sourceType.getValue().getString();
	}
    }

    char *path = db_path_to_string(pathp);
    SoSeparator *leaf = realize_instance_leaf_separator(tsp);
    if (sharedVListShape) {
	SoBRLVListShape *vlistShape = new SoBRLVListShape;
	assign_realized_identity(vlistShape, tsp, path, dp->d_namep, typeLabel,
		data->revision, data->source);
	vlistShape->setSharedGeometry(sharedVListShape);
	const char *geometryKind =
	    sharedVListShape->geometryKind.getValue().getString();
	vlistShape->geometryKind = (geometryKind && geometryKind[0]) ?
	    geometryKind : "line";
	if (geometryKind && BU_STR_EQUAL(geometryKind, "annotation"))
	    vlistShape->sourceType = "annotation";
	leaf->addChild(vlistShape);
    } else {
	SoBRLMeshShape *shape = mesh_instance_for_shared_geometry(sharedMeshShape);
	assign_realized_identity(shape, tsp, path, dp->d_namep, typeLabel,
		data->revision, data->source);
	shape->setSharedGeometry(sharedMeshShape);
	const char *geometryKind =
	    sharedMeshShape->geometryKind.getValue().getString();
	shape->geometryKind = (geometryKind && geometryKind[0]) ?
	    geometryKind : "surface";
	if (typeLabel && BU_STR_EQUAL(typeLabel, "bot"))
	    publish_lod_metadata_if_cached(shape, tsp->ts_dbip, dp->d_namep);
	leaf->addChild(shape);
    }
    data->source->addChild(leaf);
    data->realized_shapes++;
    if (path)
	bu_free(path, "db_path_to_string");

    return make_nop_tree();
}

SoBRLDatabaseSource::SoBRLDatabaseSource(void) :
    dbip(NULL),
    meshLod(NULL),
    meshLodBoundsValid(FALSE),
    meshLodBoundsMin(0.0f, 0.0f, 0.0f),
    meshLodBoundsMax(0.0f, 0.0f, 0.0f),
    pathSensor(NULL),
    instanceKeySensor(NULL),
    representationKeySensor(NULL),
    representationModeSensor(NULL),
    drawModeSensor(NULL),
    tessellationAbsTolSensor(NULL),
    tessellationRelTolSensor(NULL),
    tessellationNormTolSensor(NULL),
    lodBotThresholdSensor(NULL),
    sourceRevisionSensor(NULL),
    inputsRevisionSensor(NULL),
    viewRevisionSensor(NULL)
{
    SO_NODE_CONSTRUCTOR(SoBRLDatabaseSource);

    SO_NODE_DEFINE_ENUM_VALUE(DrawMode, WIREFRAME);
    SO_NODE_DEFINE_ENUM_VALUE(DrawMode, SHADED);
    SO_NODE_DEFINE_ENUM_VALUE(RealizationStatus, UNREALIZED);
    SO_NODE_DEFINE_ENUM_VALUE(RealizationStatus, REALIZED);
    SO_NODE_DEFINE_ENUM_VALUE(RealizationStatus, FAILED);
    SO_NODE_DEFINE_ENUM_VALUE(MaterialPolicy, MATERIAL_INHERIT);
    SO_NODE_DEFINE_ENUM_VALUE(MaterialPolicy, MATERIAL_DATABASE);

    SO_NODE_ADD_FIELD(instanceKey, (""));
    SO_NODE_ADD_FIELD(path, (""));
    SO_NODE_ADD_FIELD(displayName, (""));
    SO_NODE_ADD_FIELD(representationKey, (""));
    SO_NODE_ADD_FIELD(representationMode, (-1));
    SO_NODE_ADD_FIELD(auxiliarySource, (FALSE));
    SO_NODE_ADD_FIELD(drawMode, (WIREFRAME));
    SO_NODE_SET_SF_ENUM_TYPE(drawMode, DrawMode);
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(highlighted, (FALSE));
    SO_NODE_ADD_FIELD(lineStyle, (0));
    SO_NODE_ADD_FIELD(lineWidth, (0));
    SO_NODE_ADD_FIELD(transparency, (0.0f));
    SO_NODE_ADD_FIELD(materialColorValid, (FALSE));
    SO_NODE_ADD_FIELD(materialColor, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(materialRevision, (0));
    SO_NODE_ADD_FIELD(materialPolicy, (MATERIAL_INHERIT));
    SO_NODE_SET_SF_ENUM_TYPE(materialPolicy, MaterialPolicy);
    SO_NODE_ADD_FIELD(colorOverride, (FALSE));
    SO_NODE_ADD_FIELD(color, (SbColor(1.0f, 1.0f, 1.0f)));
    SO_NODE_ADD_FIELD(drawMatrixValid, (FALSE));
    SO_NODE_ADD_FIELD(drawMatrix, (SbMatrix::identity()));
    SO_NODE_ADD_FIELD(drawCenterValid, (FALSE));
    SO_NODE_ADD_FIELD(drawCenter, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(drawSizeValid, (FALSE));
    SO_NODE_ADD_FIELD(drawSize, (0.0f));
    SO_NODE_ADD_FIELD(sourceBoundsValid, (FALSE));
    SO_NODE_ADD_FIELD(sourceBoundsMin, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(sourceBoundsMax, (SbVec3f(0.0f, 0.0f, 0.0f)));
    SO_NODE_ADD_FIELD(tessellationAbsTol, (0.0f));
    SO_NODE_ADD_FIELD(tessellationRelTol, (0.01f));
    SO_NODE_ADD_FIELD(tessellationNormTol, (0.0f));
    SO_NODE_ADD_FIELD(lodBotThreshold, (0));
    SO_NODE_ADD_FIELD(sourceRevision, (0));
    SO_NODE_ADD_FIELD(inputsRevision, (0));
    SO_NODE_ADD_FIELD(viewRevision, (0));
    SO_NODE_ADD_FIELD(realizedRevision, (0));
    SO_NODE_ADD_FIELD(realizedSourceRevision, (0));
    SO_NODE_ADD_FIELD(realizedInputsRevision, (0));
    SO_NODE_ADD_FIELD(realizedViewRevision, (0));
    SO_NODE_ADD_FIELD(realizationStatus, (UNREALIZED));
    SO_NODE_SET_SF_ENUM_TYPE(realizationStatus, RealizationStatus);
    SO_NODE_ADD_FIELD(realizationDiagnostic, (""));
    SO_NODE_ADD_FIELD(realizationIdentity, (""));
    SO_NODE_ADD_FIELD(realizationRoleFlags, (REALIZATION_ROLE_NONE));
    SO_NODE_ADD_FIELD(realizationViewDependent, (FALSE));
    SO_NODE_ADD_FIELD(realizationViewScale, (0.0f));
    SO_NODE_ADD_FIELD(realizationBotThreshold, (0));
    SO_NODE_ADD_FIELD(realizationCurveScale, (0.0f));
    SO_NODE_ADD_FIELD(realizationPointScale, (0.0f));
    SO_NODE_ADD_FIELD(stale, (TRUE));
    SO_NODE_ADD_FIELD(staleReason, (STALE_SOURCE));

    this->attachFieldSensors();
}

SoBRLDatabaseSource::~SoBRLDatabaseSource(void)
{
    this->detachFieldSensors();
    this->clearMeshLod();
}

void
SoBRLDatabaseSource::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLDatabaseSource, SoSeparator, "Separator");
}

void
SoBRLDatabaseSource::fieldSensorCB(void *data, SoSensor *sensor)
{
    SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(data);
    if (!source)
	return;

    if (sensor == source->inputsRevisionSensor)
	source->markStale(STALE_INPUTS);
    else if (sensor == source->viewRevisionSensor ||
	    sensor == source->lodBotThresholdSensor)
	source->markStale(STALE_VIEW);
    else if (sensor == source->drawModeSensor ||
	    sensor == source->representationModeSensor)
	source->markStale(STALE_DRAW);
    else if (sensor == source->tessellationAbsTolSensor ||
	    sensor == source->tessellationRelTolSensor ||
	    sensor == source->tessellationNormTolSensor)
	source->markStale(STALE_TESSELLATION);
    else
	source->markStale(STALE_SOURCE);
}

void
SoBRLDatabaseSource::attachFieldSensors(void)
{
    this->pathSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->pathSensor->setPriority(0);
    this->pathSensor->attach(&this->path);

    this->instanceKeySensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->instanceKeySensor->setPriority(0);
    this->instanceKeySensor->attach(&this->instanceKey);

    this->representationKeySensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->representationKeySensor->setPriority(0);
    this->representationKeySensor->attach(&this->representationKey);

    this->representationModeSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->representationModeSensor->setPriority(0);
    this->representationModeSensor->attach(&this->representationMode);

    this->drawModeSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->drawModeSensor->setPriority(0);
    this->drawModeSensor->attach(&this->drawMode);

    this->tessellationAbsTolSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->tessellationAbsTolSensor->setPriority(0);
    this->tessellationAbsTolSensor->attach(&this->tessellationAbsTol);

    this->tessellationRelTolSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->tessellationRelTolSensor->setPriority(0);
    this->tessellationRelTolSensor->attach(&this->tessellationRelTol);

    this->tessellationNormTolSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->tessellationNormTolSensor->setPriority(0);
    this->tessellationNormTolSensor->attach(&this->tessellationNormTol);

    this->lodBotThresholdSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->lodBotThresholdSensor->setPriority(0);
    this->lodBotThresholdSensor->attach(&this->lodBotThreshold);

    this->sourceRevisionSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->sourceRevisionSensor->setPriority(0);
    this->sourceRevisionSensor->attach(&this->sourceRevision);

    this->inputsRevisionSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->inputsRevisionSensor->setPriority(0);
    this->inputsRevisionSensor->attach(&this->inputsRevision);

    this->viewRevisionSensor = new SoFieldSensor(SoBRLDatabaseSource::fieldSensorCB, this);
    this->viewRevisionSensor->setPriority(0);
    this->viewRevisionSensor->attach(&this->viewRevision);
}

void
SoBRLDatabaseSource::detachFieldSensors(void)
{
    if (this->pathSensor)
	this->pathSensor->detach();
    delete this->pathSensor;
    this->pathSensor = NULL;
    if (this->instanceKeySensor)
	this->instanceKeySensor->detach();
    delete this->instanceKeySensor;
    this->instanceKeySensor = NULL;
    if (this->representationKeySensor)
	this->representationKeySensor->detach();
    delete this->representationKeySensor;
    this->representationKeySensor = NULL;
    if (this->representationModeSensor)
	this->representationModeSensor->detach();
    delete this->representationModeSensor;
    this->representationModeSensor = NULL;
    if (this->drawModeSensor)
	this->drawModeSensor->detach();
    delete this->drawModeSensor;
    this->drawModeSensor = NULL;
    if (this->tessellationAbsTolSensor)
	this->tessellationAbsTolSensor->detach();
    delete this->tessellationAbsTolSensor;
    this->tessellationAbsTolSensor = NULL;
    if (this->tessellationRelTolSensor)
	this->tessellationRelTolSensor->detach();
    delete this->tessellationRelTolSensor;
    this->tessellationRelTolSensor = NULL;
    if (this->tessellationNormTolSensor)
	this->tessellationNormTolSensor->detach();
    delete this->tessellationNormTolSensor;
    this->tessellationNormTolSensor = NULL;
    if (this->lodBotThresholdSensor)
	this->lodBotThresholdSensor->detach();
    delete this->lodBotThresholdSensor;
    this->lodBotThresholdSensor = NULL;
    if (this->sourceRevisionSensor)
	this->sourceRevisionSensor->detach();
    delete this->sourceRevisionSensor;
    this->sourceRevisionSensor = NULL;
    if (this->inputsRevisionSensor)
	this->inputsRevisionSensor->detach();
    delete this->inputsRevisionSensor;
    this->inputsRevisionSensor = NULL;
    if (this->viewRevisionSensor)
	this->viewRevisionSensor->detach();
    delete this->viewRevisionSensor;
    this->viewRevisionSensor = NULL;
}

void
SoBRLDatabaseSource::markStale(void)
{
    this->markStale(STALE_SOURCE);
}

void
SoBRLDatabaseSource::markStale(uint32_t reason)
{
    if (!reason)
	return;
    this->stale = TRUE;
    this->staleReason = this->staleReason.getValue() | reason;
    this->realizationStatus = UNREALIZED;
    this->realizationDiagnostic = "";
    this->syncRealizedShapeOwnerState();
}

int
SoBRLDatabaseSource::setDrawModeState(int nextDrawMode)
{
    if (nextDrawMode != SHADED)
	nextDrawMode = WIREFRAME;
    if (this->drawMode.getValue() == nextDrawMode)
	return 0;

    const SbBool preserveExternalRealization =
	(this->realizationRoleFlags.getValue() & REALIZATION_ROLE_EXTERNAL) &&
	this->realizationStatus.getValue() == REALIZED &&
	!this->stale.getValue() &&
	((nextDrawMode == SHADED && this->getRealizedMeshCount() > 0) ||
	 (nextDrawMode == WIREFRAME && this->getRealizedShapeCount() > 0));

    if (this->drawModeSensor)
	this->drawModeSensor->detach();
    this->drawMode = nextDrawMode;
    if (this->drawModeSensor)
	this->drawModeSensor->attach(&this->drawMode);
    if (!preserveExternalRealization)
	this->markStale(STALE_DRAW);
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::setDisplayNameState(const char *name)
{
    const char *nextName = name ? name : "";
    if (database_source_string_equal(this->displayName.getValue(), nextName))
	return 0;

    this->displayName = nextName;
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::setRepresentationState(
	const char *sourceRepresentationKey,
	int sourceRepresentationMode)
{
    const char *nextKey = sourceRepresentationKey ?
	sourceRepresentationKey : "";
    const int nextMode = sourceRepresentationMode;
    if (database_source_string_equal(this->representationKey.getValue(),
	    nextKey) && this->representationMode.getValue() == nextMode)
	return 0;

    if (this->representationKeySensor)
	this->representationKeySensor->detach();
    if (this->representationModeSensor)
	this->representationModeSensor->detach();
    this->representationKey = nextKey;
    this->representationMode = nextMode;
    if (this->representationKeySensor)
	this->representationKeySensor->attach(&this->representationKey);
    if (this->representationModeSensor)
	this->representationModeSensor->attach(&this->representationMode);

    this->markStale(STALE_DRAW);
    const char *sourcePath = this->path.getValue().getString();
    for (int i = 0; i < this->getNumChildren(); i++)
	retarget_realized_node_source(this->getChild(i), this,
		sourcePath, sourcePath,
		this->sourceRevision.getValue());
    this->realizationIdentity = source_realization_identity(this);
    return 1;
}

int
SoBRLDatabaseSource::setMaterialPolicyState(int nextMaterialPolicy)
{
    if (nextMaterialPolicy != MATERIAL_DATABASE)
	nextMaterialPolicy = MATERIAL_INHERIT;
    if (this->materialPolicy.getValue() == nextMaterialPolicy)
	return 0;

    this->materialPolicy = nextMaterialPolicy;
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::setRealizationState(int nextStatus,
	uint32_t nextRealizedSourceRevision,
	uint32_t nextRealizedInputsRevision,
	uint32_t nextStaleReason,
	const char *diagnostic)
{
    if (nextStatus != REALIZED && nextStatus != UNREALIZED &&
	    nextStatus != FAILED)
	nextStatus = UNREALIZED;

    const SbBool realized = nextStatus == REALIZED ? TRUE : FALSE;
    if (realized) {
	nextRealizedSourceRevision = nextRealizedSourceRevision ?
	    nextRealizedSourceRevision : this->sourceRevision.getValue();
	nextRealizedInputsRevision = nextRealizedInputsRevision ?
	    nextRealizedInputsRevision : this->inputsRevision.getValue();
	nextStaleReason = STALE_NONE;
    } else if (!nextStaleReason) {
	nextStaleReason = STALE_SOURCE;
    }

    const SbBool nextStale = realized ? FALSE : TRUE;
    SbString nextDiagnostic = diagnostic ? diagnostic : "";

    int changed = 0;
    if (this->realizationStatus.getValue() != nextStatus) {
	this->realizationStatus = nextStatus;
	changed = 1;
    }
    if (this->stale.getValue() != nextStale) {
	this->stale = nextStale;
	changed = 1;
    }
    if (this->staleReason.getValue() != nextStaleReason) {
	this->staleReason = nextStaleReason;
	changed = 1;
    }
    if (strcmp(this->realizationDiagnostic.getValue().getString(),
	    nextDiagnostic.getString()) != 0) {
	this->realizationDiagnostic = nextDiagnostic;
	changed = 1;
    }
    if (realized &&
	    this->realizedSourceRevision.getValue() !=
	    nextRealizedSourceRevision) {
	this->realizedSourceRevision = nextRealizedSourceRevision;
	changed = 1;
    }
    if (realized &&
	    this->realizedInputsRevision.getValue() !=
	    nextRealizedInputsRevision) {
	this->realizedInputsRevision = nextRealizedInputsRevision;
	changed = 1;
    }
    if (realized &&
	    this->realizedViewRevision.getValue() !=
	    this->viewRevision.getValue()) {
	this->realizedViewRevision = this->viewRevision.getValue();
	changed = 1;
    }
    if (realized &&
	    this->realizedRevision.getValue() != nextRealizedSourceRevision) {
	this->realizedRevision = nextRealizedSourceRevision;
	changed = 1;
    }
    if (changed)
	this->syncRealizedShapeOwnerState();
    return changed;
}

int
SoBRLDatabaseSource::setRealizationRoleFlags(int roleFlags)
{
    const int validFlags = REALIZATION_ROLE_CSG | REALIZATION_ROLE_MESH |
	REALIZATION_ROLE_EXTERNAL;
    roleFlags &= validFlags;
    if (this->realizationRoleFlags.getValue() == roleFlags)
	return 0;

    this->realizationRoleFlags = roleFlags;
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::setRealizationViewPolicy(SbBool viewDependent,
	float viewScale,
	uint32_t botThreshold,
	float curveScale,
	float pointScale)
{
    int changed = 0;
    if (this->realizationViewDependent.getValue() != viewDependent) {
	this->realizationViewDependent = viewDependent;
	changed = 1;
    }
    if (database_source_float_different(
	    this->realizationViewScale.getValue(), viewScale)) {
	this->realizationViewScale = viewScale;
	changed = 1;
    }
    if (this->realizationBotThreshold.getValue() != botThreshold) {
	this->realizationBotThreshold = botThreshold;
	changed = 1;
    }
    if (database_source_float_different(
	    this->realizationCurveScale.getValue(), curveScale)) {
	this->realizationCurveScale = curveScale;
	changed = 1;
    }
    if (database_source_float_different(
	    this->realizationPointScale.getValue(), pointScale)) {
	this->realizationPointScale = pointScale;
	changed = 1;
    }
    if (changed)
	this->syncRealizedShapeOwnerState();
    return changed;
}

int
SoBRLDatabaseSource::setDisplayState(SbBool sourceRevisionValid,
	uint32_t nextSourceRevision,
	uint32_t nextInputsRevision,
	SbBool nextVisible,
	SbBool nextHighlighted,
	int nextLineStyle,
	int nextLineWidth,
	float nextTransparency,
	SbBool nextColorOverride,
	const SbColor &nextColor,
	SbBool nextMaterialColorValid,
	const SbColor &nextMaterialColor,
	uint32_t nextMaterialRevision)
{
    int changed = 0;
    const SbBool hadMaterialOverride = this->materialColorValid.getValue();
    if (sourceRevisionValid &&
	    this->sourceRevision.getValue() != nextSourceRevision) {
	this->sourceRevision = nextSourceRevision;
	changed = 1;
    }
    if (this->inputsRevision.getValue() != nextInputsRevision) {
	this->inputsRevision = nextInputsRevision;
	changed = 1;
    }
    if (this->visible.getValue() != nextVisible) {
	this->visible = nextVisible;
	changed = 1;
    }
    if (this->highlighted.getValue() != nextHighlighted) {
	this->highlighted = nextHighlighted;
	changed = 1;
    }
    if (this->lineStyle.getValue() != nextLineStyle) {
	this->lineStyle = nextLineStyle;
	changed = 1;
    }
    if (this->lineWidth.getValue() != nextLineWidth) {
	this->lineWidth = nextLineWidth;
	changed = 1;
    }
    if (database_source_float_different(this->transparency.getValue(),
	    nextTransparency)) {
	this->transparency = nextTransparency;
	changed = 1;
    }
    if (this->colorOverride.getValue() != nextColorOverride) {
	this->colorOverride = nextColorOverride;
	changed = 1;
    }
    if (nextColorOverride &&
	    !database_source_color_equal(this->color.getValue(), nextColor)) {
	this->color = nextColor;
	changed = 1;
    }
    if (this->materialColorValid.getValue() != nextMaterialColorValid) {
	this->materialColorValid = nextMaterialColorValid;
	changed = 1;
    }
    if (nextMaterialColorValid &&
	    !database_source_color_equal(this->materialColor.getValue(),
		nextMaterialColor)) {
	this->materialColor = nextMaterialColor;
	changed = 1;
    }
    if (this->materialRevision.getValue() != nextMaterialRevision) {
	this->materialRevision = nextMaterialRevision;
	changed = 1;
    }
    if (hadMaterialOverride && !nextMaterialColorValid)
	this->markStale(STALE_SOURCE);
    if (changed)
	this->syncRealizedShapeOwnerState();
    return changed;
}

int
SoBRLDatabaseSource::applyDisplayPatch(
	const BRLObolDatabaseSourceDisplayPatch &patch)
{
    int changed = 0;
    if (patch.visibleValid && this->visible.getValue() != patch.visible) {
	this->visible = patch.visible;
	changed = 1;
    }
    if (patch.highlightedValid &&
	    this->highlighted.getValue() != patch.highlighted) {
	this->highlighted = patch.highlighted;
	changed = 1;
    }
    if (patch.lineStyleValid &&
	    this->lineStyle.getValue() != patch.lineStyle) {
	this->lineStyle = patch.lineStyle;
	changed = 1;
    }
    if (patch.lineWidthValid &&
	    this->lineWidth.getValue() != patch.lineWidth) {
	this->lineWidth = patch.lineWidth;
	changed = 1;
    }
    if (patch.transparencyValid &&
	    database_source_float_different(this->transparency.getValue(),
		patch.transparency)) {
	this->transparency = patch.transparency;
	changed = 1;
    }
    if (patch.colorOverrideValid &&
	    this->colorOverride.getValue() != patch.colorOverride) {
	this->colorOverride = patch.colorOverride;
	changed = 1;
    }
    if (patch.colorValid &&
	    !database_source_color_equal(this->color.getValue(), patch.color)) {
	this->color = patch.color;
	changed = 1;
    }
    if (changed)
	this->syncRealizedShapeOwnerState();
    return changed;
}

int
SoBRLDatabaseSource::setPlacementState(SbBool nextDrawMatrixValid,
	const SbMatrix &nextDrawMatrix,
	SbBool nextDrawCenterValid,
	const SbVec3f &nextDrawCenter,
	SbBool nextDrawSizeValid,
	float nextDrawSize)
{
    int changed = 0;
    if (this->drawMatrixValid.getValue() != nextDrawMatrixValid) {
	this->drawMatrixValid = nextDrawMatrixValid;
	changed = 1;
    }
    if (!this->drawMatrix.getValue().equals(nextDrawMatrix, 0.000001f)) {
	this->drawMatrix = nextDrawMatrix;
	changed = 1;
    }
    if (this->drawCenterValid.getValue() != nextDrawCenterValid) {
	this->drawCenterValid = nextDrawCenterValid;
	changed = 1;
    }
    const SbVec3f currentCenter = this->drawCenter.getValue();
    if (database_source_float_different(currentCenter[0],
	    nextDrawCenter[0]) ||
	    database_source_float_different(currentCenter[1],
		nextDrawCenter[1]) ||
	    database_source_float_different(currentCenter[2],
		nextDrawCenter[2])) {
	this->drawCenter = nextDrawCenter;
	changed = 1;
    }
    if (this->drawSizeValid.getValue() != nextDrawSizeValid) {
	this->drawSizeValid = nextDrawSizeValid;
	changed = 1;
    }
    if (database_source_float_different(this->drawSize.getValue(),
	    nextDrawSize)) {
	this->drawSize = nextDrawSize;
	changed = 1;
    }
    if (sync_source_placement_transform(this))
	changed = 1;
    if (changed)
	this->syncRealizedShapeOwnerState();
    return changed;
}

int
SoBRLDatabaseSource::setSourceBoundsState(SbBool nextBoundsValid,
	const SbVec3f &nextBoundsMin,
	const SbVec3f &nextBoundsMax)
{
    SbVec3f sanitizedMin(0.0f, 0.0f, 0.0f);
    SbVec3f sanitizedMax(0.0f, 0.0f, 0.0f);
    if (nextBoundsValid) {
	const SbBox3f bounds =
	    database_source_box_from_minmax(nextBoundsMin, nextBoundsMax);
	sanitizedMin = bounds.getMin();
	sanitizedMax = bounds.getMax();
    }

    int changed = 0;
    if (this->sourceBoundsValid.getValue() != nextBoundsValid) {
	this->sourceBoundsValid = nextBoundsValid;
	changed = 1;
    }
    if (!database_source_vec3f_equal(this->sourceBoundsMin.getValue(),
	    sanitizedMin)) {
	this->sourceBoundsMin = sanitizedMin;
	changed = 1;
    }
    if (!database_source_vec3f_equal(this->sourceBoundsMax.getValue(),
	    sanitizedMax)) {
	this->sourceBoundsMax = sanitizedMax;
	changed = 1;
    }

    return changed;
}

void
SoBRLDatabaseSource::clearSourceBounds(void)
{
    (void)this->setSourceBoundsState(FALSE,
	    SbVec3f(0.0f, 0.0f, 0.0f),
	    SbVec3f(0.0f, 0.0f, 0.0f));
}

SbBool
SoBRLDatabaseSource::getSourceBounds(SbBox3f &bounds) const
{
    bounds.makeEmpty();
    if (!this->sourceBoundsValid.getValue())
	return FALSE;

    bounds = database_source_box_from_minmax(
	    this->sourceBoundsMin.getValue(),
	    this->sourceBoundsMax.getValue());
    return bounds.isEmpty() ? FALSE : TRUE;
}

SbBool
SoBRLDatabaseSource::getEffectiveSourceBounds(SbBox3f &bounds) const
{
    if (!this->getSourceBounds(bounds))
	return FALSE;

    const SoMatrixTransform *transform = source_placement_transform(this);
    if (transform) {
	bounds = database_source_transform_bounds(bounds,
		transform->matrix.getValue());
    } else if (this->drawMatrixValid.getValue()) {
	bounds = database_source_transform_bounds(bounds,
		this->drawMatrix.getValue());
    }

    return bounds.isEmpty() ? FALSE : TRUE;
}

void
SoBRLDatabaseSource::setDatabase(struct db_i *database)
{
    if (this->dbip != database) {
	this->dbip = database;
	this->markStale(STALE_DATABASE);
    }
}

struct db_i *
SoBRLDatabaseSource::getDatabase(void) const
{
    return this->dbip;
}

void
SoBRLDatabaseSource::setMeshLod(struct rt_mesh_lod *lod)
{
    if (this->meshLod == lod) {
	if (!lod)
	    this->clearMeshLodBounds();
	return;
    }

    if (this->meshLod)
	rt_mesh_lod_destroy(this->meshLod);

    this->meshLod = lod;
    this->clearMeshLodBounds();
}

struct rt_mesh_lod *
SoBRLDatabaseSource::getMeshLod(void) const
{
    return this->meshLod;
}

void
SoBRLDatabaseSource::clearMeshLod(void)
{
    if (this->meshLod)
	rt_mesh_lod_destroy(this->meshLod);
    this->meshLod = NULL;
    this->clearMeshLodBounds();
}

int
SoBRLDatabaseSource::setMeshLodBounds(const SbVec3f &bmin,
	const SbVec3f &bmax)
{
    if (bmin[0] > bmax[0] || bmin[1] > bmax[1] || bmin[2] > bmax[2])
	return 0;

    this->meshLodBoundsMin = bmin;
    this->meshLodBoundsMax = bmax;
    this->meshLodBoundsValid = TRUE;
    return 1;
}

SbBool
SoBRLDatabaseSource::getMeshLodBounds(SbVec3f &bmin,
	SbVec3f &bmax) const
{
    if (!this->meshLodBoundsValid)
	return FALSE;

    bmin = this->meshLodBoundsMin;
    bmax = this->meshLodBoundsMax;
    return TRUE;
}

void
SoBRLDatabaseSource::clearMeshLodBounds(void)
{
    this->meshLodBoundsValid = FALSE;
    this->meshLodBoundsMin.setValue(0.0f, 0.0f, 0.0f);
    this->meshLodBoundsMax.setValue(0.0f, 0.0f, 0.0f);
}

void
SoBRLDatabaseSource::configureDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int mode,
	uint32_t revision)
{
    this->configureDatabaseSourceInstance(sourcePath, sourcePath, database,
	    mode, revision);
}

void
SoBRLDatabaseSource::configureDatabaseSourceInstance(
	const char *sourceInstanceKey,
	const char *sourcePath,
	struct db_i *database,
	int mode,
	uint32_t revision)
{
    this->configureDatabaseSourceInstanceRepresentation(sourceInstanceKey,
	    sourcePath, NULL, -1, database, mode, revision);
}

void
SoBRLDatabaseSource::configureDatabaseSourceInstanceRepresentation(
	const char *sourceInstanceKey,
	const char *sourcePath,
	const char *sourceRepresentationKey,
	int sourceRepresentationMode,
	struct db_i *database,
	int mode,
	uint32_t revision)
{
    int sanitizedMode = mode == SHADED ? SHADED : WIREFRAME;
    uint32_t reason = STALE_NONE;
    if (this->dbip != database)
	reason |= STALE_DATABASE;
    if (this->drawMode.getValue() != sanitizedMode)
	reason |= STALE_DRAW;
    const std::string stableSourcePath = sourcePath ? sourcePath : "";
    const std::string stableInstanceKey =
	(sourceInstanceKey && sourceInstanceKey[0]) ?
	sourceInstanceKey : stableSourcePath.c_str();
    const char *effectiveInstanceKey = stableInstanceKey.c_str();
    if (strcmp(this->instanceKey.getValue().getString(),
	    effectiveInstanceKey) != 0)
	reason |= STALE_SOURCE;
    const std::string stableRepresentationKey =
	(sourceRepresentationKey && sourceRepresentationKey[0]) ?
	sourceRepresentationKey : stableInstanceKey.c_str();
    const char *effectiveRepresentationKey = stableRepresentationKey.c_str();
    if (strcmp(this->representationKey.getValue().getString(),
	    effectiveRepresentationKey) != 0)
	reason |= STALE_SOURCE;
    if (this->representationMode.getValue() != sourceRepresentationMode)
	reason |= STALE_DRAW;
    if (strcmp(this->path.getValue().getString(),
	    stableSourcePath.c_str()) != 0)
	reason |= STALE_SOURCE;
    if (this->sourceRevision.getValue() != revision)
	reason |= STALE_SOURCE;

    this->detachFieldSensors();
    this->dbip = database;
    this->instanceKey = effectiveInstanceKey;
    this->path = stableSourcePath.c_str();
    this->representationKey = effectiveRepresentationKey;
    this->representationMode = sourceRepresentationMode;
    this->auxiliarySource = FALSE;
    this->drawMode = sanitizedMode;
    this->sourceRevision = revision;
    this->markStale(reason);
    this->attachFieldSensors();
}

template <typename ShapeT>
static void
retarget_realized_shape_source(ShapeT *shape,
	const SoBRLDatabaseSource *source,
	const char *oldSourcePath,
	const char *newSourcePath,
	uint32_t revision)
{
    if (!shape || !newSourcePath)
	return;

    const char *oldName = oldSourcePath ? lookup_name_from_path(oldSourcePath) :
	NULL;
    const char *newName = lookup_name_from_path(SbString(newSourcePath));
    if (!newName)
	newName = newSourcePath;

    const std::string recordRole = shape->recordRole.getValue().getString();
    const int auxiliary = BU_STR_EQUAL(recordRole.c_str(), "auxiliary");
    const std::string shapeSourceName =
	shape->sourceName.getValue().getString();
    const std::string displayName = shape->displayName.getValue().getString();
    const std::string geometryName =
	shape->geometryName.getValue().getString();
    const std::string sourceDisplayName = source ?
	source->displayName.getValue().getString() : "";

    database_source_assign_string(shape->sourcePath, newSourcePath);
    shape->sourceId = revision;
    if (!auxiliary &&
	    (shapeSourceName.empty() ||
	     (oldName && BU_STR_EQUAL(shapeSourceName.c_str(), oldName)) ||
	     (oldSourcePath && BU_STR_EQUAL(shapeSourceName.c_str(), oldSourcePath))))
	database_source_assign_string(shape->sourceName, newName);
    if (!auxiliary && !sourceDisplayName.empty())
	database_source_assign_string(shape->displayName,
		sourceDisplayName.c_str());
    else if (!auxiliary &&
	    (displayName.empty() ||
	     (oldName && BU_STR_EQUAL(displayName.c_str(), oldName)) ||
	     (oldSourcePath && BU_STR_EQUAL(displayName.c_str(), oldSourcePath))))
	database_source_assign_string(shape->displayName, newName);
    if (!auxiliary &&
	    (geometryName.empty() ||
	     (oldName && BU_STR_EQUAL(geometryName.c_str(), oldName)) ||
	     (oldSourcePath && BU_STR_EQUAL(geometryName.c_str(), oldSourcePath))))
	database_source_assign_string(shape->geometryName, newName);

    SbString identity = source_record_identity(source, newSourcePath);
    if (auxiliary) {
	const char *auxName = shape->geometryName.getValue().getString();
	if (auxName && auxName[0]) {
	    identity += "::";
	    identity += auxName;
	}
    }
    database_source_assign_string(shape->sourceIdentity, identity);
    database_source_assign_string(shape->cacheIdentity,
	    record_identity_with_revision(identity.getString(), revision));
    sync_shape_owner_state(shape, source);
}

static void
retarget_material_object_source(SoBRLMaterialObject *object,
	const char *oldSourcePath,
	const char *newSourcePath,
	uint32_t revision)
{
    if (!object || !newSourcePath)
	return;

    const char *oldName = oldSourcePath ? lookup_name_from_path(oldSourcePath) :
	NULL;
    const char *newName = lookup_name_from_path(SbString(newSourcePath));
    if (!newName)
	newName = newSourcePath;

    const char *sourceName = object->sourceName.getValue().getString();
    object->sourcePath = newSourcePath;
    object->sourceId = revision;
    if (!sourceName || !sourceName[0] ||
	    (oldName && BU_STR_EQUAL(sourceName, oldName)) ||
	    (oldSourcePath && BU_STR_EQUAL(sourceName, oldSourcePath)))
	object->sourceName = newName;
}

static void
retarget_realized_node_source(SoNode *node,
	const SoBRLDatabaseSource *source,
	const char *oldSourcePath,
	const char *newSourcePath,
	uint32_t revision)
{
    if (!node || !newSourcePath)
	return;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	retarget_realized_shape_source(static_cast<SoBRLVListShape *>(node),
		source, oldSourcePath, newSourcePath, revision);
	return;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	retarget_realized_shape_source(static_cast<SoBRLMeshShape *>(node),
		source, oldSourcePath, newSourcePath, revision);
	return;
    }
    if (node->isOfType(SoBRLMaterialObject::getClassTypeId())) {
	retarget_material_object_source(
		static_cast<SoBRLMaterialObject *>(node), oldSourcePath,
		newSourcePath, revision);
	return;
    }
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    retarget_realized_node_source(group->getChild(i), source,
		    oldSourcePath, newSourcePath, revision);
    }
}

int
SoBRLDatabaseSource::retargetDatabaseSource(const char *sourcePath,
	uint32_t revision)
{
    return this->retargetDatabaseSourceInstance(sourcePath, sourcePath,
	    revision);
}

int
SoBRLDatabaseSource::retargetDatabaseSourceInstance(
	const char *sourceInstanceKey,
	const char *sourcePath,
	uint32_t revision)
{
    if (!sourcePath || !sourcePath[0])
	return -1;

    const std::string stableSourcePath = sourcePath;
    const std::string stableInstanceKey =
	(sourceInstanceKey && sourceInstanceKey[0]) ?
	sourceInstanceKey : stableSourcePath.c_str();
    const char *effectiveInstanceKey = stableInstanceKey.c_str();
    const char *oldPath = this->path.getValue().getString();
    const char *oldInstanceKey = this->instanceKey.getValue().getString();
    const int pathChanged = !oldPath ||
	!BU_STR_EQUAL(oldPath, stableSourcePath.c_str());
    const int instanceChanged = !oldInstanceKey ||
	!BU_STR_EQUAL(oldInstanceKey, effectiveInstanceKey);
    if (revision == 0)
	revision = this->sourceRevision.getValue() +
	    ((pathChanged || instanceChanged) ? 1 : 0);
    const int revisionChanged = this->sourceRevision.getValue() != revision;
    if (!pathChanged && !instanceChanged && !revisionChanged)
	return 0;
    const int sourceChanged = pathChanged || revisionChanged;

    std::string oldPathCopy = oldPath ? oldPath : "";

    this->detachFieldSensors();
    this->instanceKey = effectiveInstanceKey;
    this->path = stableSourcePath.c_str();
    this->sourceRevision = revision;
    if (sourceChanged)
	this->markStale(STALE_SOURCE);
    for (int i = 0; i < this->getNumChildren(); i++)
	retarget_realized_node_source(this->getChild(i), this,
		oldPathCopy.c_str(), stableSourcePath.c_str(), revision);
    this->realizationIdentity = source_realization_identity(this);
    this->attachFieldSensors();
    return 1;
}

SbBool
SoBRLDatabaseSource::needsRealization(void) const
{
    return this->stale.getValue() ||
	this->realizedSourceRevision.getValue() != this->sourceRevision.getValue() ||
	this->realizedInputsRevision.getValue() != this->inputsRevision.getValue() ||
	this->realizedViewRevision.getValue() != this->viewRevision.getValue();
}

SbBool
SoBRLDatabaseSource::realizeDatabaseWireframe(void)
{
    this->realizationDiagnostic = "";
    if (!this->dbip) {
	this->realizationDiagnostic = "database source has no database";
	return FALSE;
    }

    const char *treeName = lookup_name_from_path(this->path.getValue());
    if (!treeName) {
	this->realizationDiagnostic = "database source path is empty";
	return FALSE;
    }

    remove_non_auxiliary_children(this);
    (void)remove_source_placement_transform(this);

    struct db_tree_state init_state;
    db_init_db_tree_state(&init_state, this->dbip);
    init_state.ts_stop_at_regions = 0;

    struct realize_walk_data data;
    data.source = this;
    data.revision = this->sourceRevision.getValue();
    data.realized_shapes = 0;
    data.failed_shapes = 0;

    const char *av[1] = { treeName };
    int ret = db_walk_tree_leaf_instances(this->dbip, 1, av, 1, &init_state,
	    NULL, NULL, realize_leaf, &data);
    db_free_db_tree_state(&init_state);

    if (ret < 0 || data.realized_shapes <= 0 || data.failed_shapes > 0) {
	remove_non_auxiliary_children(this);
	this->realizationIdentity = "";
	if (data.diagnostic.getLength() > 0) {
	    this->realizationDiagnostic = data.diagnostic;
	} else if (data.realized_shapes <= 0) {
	    SbString msg;
	    msg.sprintf("%s: no drawable wireframe geometry realized", treeName);
	    this->realizationDiagnostic = msg;
	} else {
	    SbString msg;
	    msg.sprintf("%s: wireframe realization failed", treeName);
	    this->realizationDiagnostic = msg;
	}
	return FALSE;
    }

    this->realizedRevision = this->sourceRevision.getValue();
    this->realizedSourceRevision = this->sourceRevision.getValue();
    this->realizedInputsRevision = this->inputsRevision.getValue();
    this->realizedViewRevision = this->viewRevision.getValue();
    this->realizationStatus = REALIZED;
    this->realizationDiagnostic = "";
    this->realizationIdentity = source_realization_identity(this);
    this->stale = FALSE;
    this->staleReason = STALE_NONE;
    this->syncRealizedShapeOwnerState();
    return TRUE;
}

SbBool
SoBRLDatabaseSource::realizeDatabaseMesh(void)
{
    this->realizationDiagnostic = "";
    if (!this->dbip) {
	this->realizationDiagnostic = "database source has no database";
	return FALSE;
    }

    const char *treeName = lookup_name_from_path(this->path.getValue());
    if (!treeName) {
	this->realizationDiagnostic = "database source path is empty";
	return FALSE;
    }

    remove_non_auxiliary_children(this);
    (void)remove_source_placement_transform(this);

    struct db_tree_state init_state;
    db_init_db_tree_state(&init_state, this->dbip);
    init_state.ts_stop_at_regions = 0;

    struct realize_walk_data data;
    data.source = this;
    data.revision = this->sourceRevision.getValue();
    data.realized_shapes = 0;
    data.failed_shapes = 0;

    const char *av[1] = { treeName };
    int ret = db_walk_tree_leaf_instances(this->dbip, 1, av, 1, &init_state,
	    NULL, NULL, realize_mesh_leaf, &data);
    db_free_db_tree_state(&init_state);

    if (ret < 0 || data.realized_shapes <= 0 || data.failed_shapes > 0) {
	remove_non_auxiliary_children(this);
	this->realizationIdentity = "";
	if (data.diagnostic.getLength() > 0) {
	    this->realizationDiagnostic = data.diagnostic;
	} else if (data.realized_shapes <= 0) {
	    SbString msg;
	    msg.sprintf("%s: no drawable mesh geometry realized", treeName);
	    this->realizationDiagnostic = msg;
	} else {
	    SbString msg;
	    msg.sprintf("%s: mesh realization failed", treeName);
	    this->realizationDiagnostic = msg;
	}
	return FALSE;
    }

    this->realizedRevision = this->sourceRevision.getValue();
    this->realizedSourceRevision = this->sourceRevision.getValue();
    this->realizedInputsRevision = this->inputsRevision.getValue();
    this->realizedViewRevision = this->viewRevision.getValue();
    this->realizationStatus = REALIZED;
    this->realizationDiagnostic = "";
    this->realizationIdentity = source_realization_identity(this);
    this->stale = FALSE;
    this->staleReason = STALE_NONE;
    this->syncRealizedShapeOwnerState();
    return TRUE;
}

SbBool
SoBRLDatabaseSource::realizePrototypeWireframe(void)
{
    this->realizationDiagnostic = "";
    SoBRLVListShape *shape = this->getRealizedShape();
    if (!shape) {
	shape = new SoBRLVListShape;
	database_source_add_realized_child(this, shape);
    }

    const float halfExtent = 1.0f + 0.25f * static_cast<float>(this->sourceRevision.getValue() % 4);
    SbVec3f points[5] = {
	SbVec3f(-halfExtent, -halfExtent, 0.0f),
	SbVec3f( halfExtent, -halfExtent, 0.0f),
	SbVec3f( halfExtent,  halfExtent, 0.0f),
	SbVec3f(-halfExtent,  halfExtent, 0.0f),
	SbVec3f(-halfExtent, -halfExtent, 0.0f)
    };
    int32_t commands[5] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW
    };

    shape->sourcePath = this->path.getValue();
    shape->sourceName = lookup_name_from_path(this->path.getValue()) ?
	lookup_name_from_path(this->path.getValue()) : this->path.getValue().getString();
    shape->sourceType = "prototype";
    shape->sourceId = this->sourceRevision.getValue();
    shape->displayName = this->displayName.getValue().getLength() > 0 ?
	this->displayName.getValue() : shape->sourceName.getValue();
    shape->geometryName = shape->sourceName.getValue();
    shape->sourceIdentity = source_record_identity(this,
	    shape->sourcePath.getValue().getString());
    shape->cacheIdentity = record_identity_with_revision(
	    shape->sourceIdentity.getValue().getString(),
	    shape->sourceId.getValue());
    shape->databaseIntent = FALSE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = TRUE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = TRUE;
    shape->drawMode = BRLOBOL_LOD_DRAW_DIAGNOSTIC;
    shape->recordRole = "prototype";
    shape->geometryKind = "";
    shape->regionId = 0;
    shape->airCode = 0;
    shape->materialId = 0;
    shape->los = 0;
    shape->materialColorValid = FALSE;
    shape->materialColor = SbColor(1.0f, 1.0f, 1.0f);
    shape->materialShader = "";
    shape->setLineSet(points, commands, 5);

    this->realizedRevision = this->sourceRevision.getValue();
    this->realizedSourceRevision = this->sourceRevision.getValue();
    this->realizedInputsRevision = this->inputsRevision.getValue();
    this->realizedViewRevision = this->viewRevision.getValue();
    this->realizationStatus = REALIZED;
    this->realizationDiagnostic = "";
    this->realizationIdentity = source_realization_identity(this);
    this->stale = FALSE;
    this->staleReason = STALE_NONE;
    this->syncRealizedShapeOwnerState();
    return TRUE;
}

static const char *
external_string_or_default(const char *value, const char *fallback)
{
    return value && value[0] ? value : fallback;
}

static const char *
external_source_leaf_name(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";

    const char *sourcePath = source->path.getValue().getString();
    const char *leaf = lookup_name_from_path(source->path.getValue());
    if (leaf && leaf[0])
	return leaf;
    return sourcePath ? sourcePath : "";
}

template <typename ShapeT>
static void
assign_external_primary_identity(ShapeT *shape,
	const SoBRLDatabaseSource *source,
	const char *sourceType,
	const char *geometryKind)
{
    if (!shape || !source)
	return;

    const char *sourcePath = source->path.getValue().getString();
    const char *sourceName = external_source_leaf_name(source);
    assign_realized_identity(shape, NULL, sourcePath, sourceName, sourceType,
	    source->sourceRevision.getValue(), source);
    shape->geometryKind = geometryKind ? geometryKind : "";
    sync_shape_placement_state(shape, source);
}

static int
external_vlist_command_valid(int32_t command)
{
    return command == SoBRLVListShape::MOVE ||
	command == SoBRLVListShape::DRAW ||
	command == SoBRLVListShape::POINT;
}

static SbBool
external_bounds_from_points(const SbVec3f *points,
	int count,
	SbVec3f &boundsMin,
	SbVec3f &boundsMax)
{
    if (!points || count <= 0)
	return FALSE;

    boundsMin = points[0];
    boundsMax = points[0];
    for (int i = 1; i < count; i++) {
	for (int axis = 0; axis < 3; axis++) {
	    if (points[i][axis] < boundsMin[axis])
		boundsMin[axis] = points[i][axis];
	    if (points[i][axis] > boundsMax[axis])
		boundsMax[axis] = points[i][axis];
	}
    }

    return TRUE;
}

static void
set_external_bounds_from_points(SoBRLDatabaseSource *source,
	const SbVec3f *points,
	int count)
{
    if (!source)
	return;

    SbVec3f boundsMin;
    SbVec3f boundsMax;
    if (external_bounds_from_points(points, count, boundsMin, boundsMax))
	(void)source->setSourceBoundsState(TRUE, boundsMin, boundsMax);
    else
	source->clearSourceBounds();
}

static void
mark_external_primary_published_current(SoBRLDatabaseSource *source)
{
    if (!source)
	return;

    const uint32_t sourceRevision = source->sourceRevision.getValue();
    source->realizedRevision = sourceRevision;
    source->realizedSourceRevision = sourceRevision;
    source->realizedInputsRevision = source->inputsRevision.getValue();
    source->realizedViewRevision = source->viewRevision.getValue();
    source->realizationStatus = SoBRLDatabaseSource::REALIZED;
    source->realizationDiagnostic = "";
    source->realizationIdentity = source_realization_identity(source);
    source->realizationRoleFlags = SoBRLDatabaseSource::REALIZATION_ROLE_EXTERNAL;
    source->stale = FALSE;
    source->staleReason = SoBRLDatabaseSource::STALE_NONE;
}

static SoBRLVListShape *
first_direct_primary_vlist_child(SoBRLDatabaseSource *source)
{
    if (!source)
	return NULL;

    for (int i = 0; i < source->getNumChildren(); i++) {
	SoNode *child = source->getChild(i);
	if (!child || !child->isOfType(SoBRLVListShape::getClassTypeId()))
	    continue;
	if (node_is_auxiliary_vlist(child))
	    continue;
	return static_cast<SoBRLVListShape *>(child);
    }

    return NULL;
}

static SoBRLMeshShape *
first_direct_primary_mesh_child(SoBRLDatabaseSource *source)
{
    if (!source)
	return NULL;

    for (int i = 0; i < source->getNumChildren(); i++) {
	SoNode *child = source->getChild(i);
	if (child && child->isOfType(SoBRLMeshShape::getClassTypeId()))
	    return static_cast<SoBRLMeshShape *>(child);
    }

    return NULL;
}

static int
remove_external_primary_children_except(SoBRLDatabaseSource *source,
	SoNode *keepVList,
	SoNode *keepMesh)
{
    if (!source)
	return 0;

    int removed = 0;
    for (int i = source->getNumChildren() - 1; i >= 0; i--) {
	SoNode *child = source->getChild(i);
	if (child == keepVList || child == keepMesh ||
		node_is_source_placement_transform(child) ||
		node_is_auxiliary_vlist(child) ||
		node_is_auxiliary_source(child))
	    continue;
	source->removeChild(i);
	removed++;
    }
    return removed;
}

static void
clear_external_primary_vlist(SoBRLDatabaseSource *source,
	SoBRLVListShape *shape)
{
    if (!source || !shape)
	return;

    const SbString sourceType = shape->sourceType.getValue();
    const SbString geometryKind = shape->geometryKind.getValue();
    assign_external_primary_identity(shape, source,
	    external_string_or_default(sourceType.getString(), "line-set"),
	    external_string_or_default(geometryKind.getString(), "line"));
    shape->setLineSet(NULL, NULL, 0);
    shape->setPrecisePoints(NULL, 0);
}

static void
clear_external_primary_mesh(SoBRLDatabaseSource *source,
	SoBRLMeshShape *shape)
{
    if (!source || !shape)
	return;

    const SbString sourceType = shape->sourceType.getValue();
    const SbString geometryKind = shape->geometryKind.getValue();
    assign_external_primary_identity(shape, source,
	    external_string_or_default(sourceType.getString(), "indexed-face-set"),
	    external_string_or_default(geometryKind.getString(), "surface"));
    shape->setIndexedTriangles(NULL, 0, NULL, 0);
}

int
SoBRLDatabaseSource::clearRealizedGeometry(SbBool preserveAuxiliary)
{
    const int before = this->getNumChildren();
    if (preserveAuxiliary) {
	remove_non_auxiliary_children(this);
    } else {
	for (int i = this->getNumChildren() - 1; i >= 0; i--) {
	    if (node_is_source_placement_transform(this->getChild(i)))
		continue;
	    this->removeChild(i);
	}
    }

    return before != this->getNumChildren() ? 1 : 0;
}

int
SoBRLDatabaseSource::clearExternalPrimaryGeometry(void)
{
    SoBRLVListShape *primaryVList = first_direct_primary_vlist_child(this);
    SoBRLMeshShape *primaryMesh = first_direct_primary_mesh_child(this);

    (void)remove_external_primary_children_except(this, primaryVList,
	    primaryMesh);

    if (!primaryVList && !primaryMesh) {
	primaryVList = new SoBRLVListShape;
	database_source_add_realized_child(this, primaryVList);
    }

    clear_external_primary_vlist(this, primaryVList);
    clear_external_primary_mesh(this, primaryMesh);
    this->clearSourceBounds();
    mark_external_primary_published_current(this);
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::publishExternalLineSet(
	const BRLObolExternalLineSet &lineSet)
{
    if (lineSet.count < 0 || (lineSet.count > 0 && !lineSet.points))
	return 0;

    std::vector<int32_t> fallbackCommands;
    const int32_t *commands = lineSet.commands;
    if (lineSet.count > 0) {
	if (!commands) {
	    fallbackCommands.reserve(lineSet.count);
	    for (int i = 0; i < lineSet.count; i++)
		fallbackCommands.push_back(i == 0 ? SoBRLVListShape::MOVE :
			SoBRLVListShape::DRAW);
	    commands = fallbackCommands.data();
	}
	for (int i = 0; i < lineSet.count; i++) {
	    if (!external_vlist_command_valid(commands[i]))
		return 0;
	}
    }

    if (lineSet.count == 0)
	return this->clearExternalPrimaryGeometry();

    (void)this->clearRealizedGeometry(TRUE);
    SoBRLVListShape *shape = new SoBRLVListShape;
    database_source_add_realized_child(this, shape);

    assign_external_primary_identity(shape, this,
	    external_string_or_default(lineSet.sourceType, "line-set"),
	    external_string_or_default(lineSet.geometryKind, "line"));
    shape->setLineSet(lineSet.points, commands, lineSet.count);
    shape->setPrecisePoints(lineSet.precisePoints, lineSet.count);
    set_external_bounds_from_points(this, lineSet.points, lineSet.count);
    mark_external_primary_published_current(this);
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::publishExternalPointSet(
	const BRLObolExternalPointSet &pointSet)
{
    if (pointSet.count < 0 || (pointSet.count > 0 && !pointSet.points))
	return 0;

    if (pointSet.count == 0)
	return this->clearExternalPrimaryGeometry();

    std::vector<int32_t> commands;
    commands.reserve(pointSet.count);
    for (int i = 0; i < pointSet.count; i++)
	commands.push_back(SoBRLVListShape::POINT);

    BRLObolExternalLineSet lineSet;
    lineSet.points = pointSet.points;
    lineSet.commands = commands.data();
    lineSet.precisePoints = pointSet.precisePoints;
    lineSet.count = pointSet.count;
    lineSet.sourceType =
	external_string_or_default(pointSet.sourceType, "point-set");
    lineSet.geometryKind =
	external_string_or_default(pointSet.geometryKind, "point");
    return this->publishExternalLineSet(lineSet);
}

int
SoBRLDatabaseSource::publishExternalTriangleMesh(
	const BRLObolExternalTriangleMesh &triangleMesh)
{
    if (triangleMesh.pointCount < 0 || triangleMesh.indexCount < 0 ||
	    (triangleMesh.pointCount > 0 && !triangleMesh.points) ||
	    (triangleMesh.indexCount > 0 && !triangleMesh.indices) ||
	    (triangleMesh.indexCount % 3) != 0)
	return 0;

    if (triangleMesh.pointCount == 0 || triangleMesh.indexCount == 0)
	return this->clearExternalPrimaryGeometry();

    for (int i = 0; i < triangleMesh.indexCount; i++) {
	const int32_t index = triangleMesh.indices[i];
	if (index < 0 || index >= triangleMesh.pointCount)
	    return 0;
    }

    (void)this->clearRealizedGeometry(TRUE);
    SoBRLMeshShape *shape = new SoBRLMeshShape;
    database_source_add_realized_child(this, shape);

    assign_external_primary_identity(shape, this,
	    external_string_or_default(triangleMesh.sourceType,
		"indexed-face-set"),
	    external_string_or_default(triangleMesh.geometryKind, "surface"));
    shape->setIndexedTriangles(triangleMesh.points, triangleMesh.pointCount,
	    triangleMesh.indices, triangleMesh.indexCount);
    set_external_bounds_from_points(this, triangleMesh.points,
	    triangleMesh.pointCount);
    mark_external_primary_published_current(this);
    this->syncRealizedShapeOwnerState();
    return 1;
}

int
SoBRLDatabaseSource::publishExternalAnnotation(
	const BRLObolExternalAnnotation &annotation)
{
    if (annotation.linePointCount < 0 || annotation.annotationPointCount < 0 ||
	    annotation.segmentCount < 0 ||
	    (annotation.linePointCount > 0 && !annotation.linePoints) ||
	    (annotation.annotationPointCount > 0 &&
	     !annotation.annotationPoints) ||
	    (annotation.segmentCount > 0 && !annotation.segments))
	return 0;

    if (annotation.lineCommands) {
	for (int i = 0; i < annotation.linePointCount; i++) {
	    if (!external_vlist_command_valid(annotation.lineCommands[i]))
		return 0;
	}
    }

    std::vector<int32_t> fallbackCommands;
    const int32_t *lineCommands = annotation.lineCommands;
    if (annotation.linePointCount > 0 && !lineCommands) {
	fallbackCommands.reserve(annotation.linePointCount);
	for (int i = 0; i < annotation.linePointCount; i++)
	    fallbackCommands.push_back(i == 0 ? SoBRLVListShape::MOVE :
		    SoBRLVListShape::DRAW);
	lineCommands = fallbackCommands.data();
    }

    (void)this->clearRealizedGeometry(TRUE);
    SoBRLVListShape *shape = new SoBRLVListShape;
    database_source_add_realized_child(this, shape);

    assign_external_primary_identity(shape, this,
	    external_string_or_default(annotation.sourceType, "annotation"),
	    external_string_or_default(annotation.geometryKind, "annotation"));
    shape->setLineSet(annotation.linePoints, lineCommands,
	    annotation.linePointCount);
    shape->setPrecisePoints(annotation.preciseLinePoints,
	    annotation.linePointCount);
    shape->annotationBasePoint = annotation.basePoint;
    if (annotation.annotationPointCount > 0)
	shape->annotationPoint.setValues(0, annotation.annotationPointCount,
		annotation.annotationPoints);
    else
	shape->annotationPoint.setNum(0);
    shape->setPreciseAnnotationPoints(annotation.preciseAnnotationPoints,
	    annotation.annotationPointCount);

    shape->annotationSegmentKind.setNum(annotation.segmentCount);
    shape->annotationSegmentStart.setNum(annotation.segmentCount);
    shape->annotationSegmentEnd.setNum(annotation.segmentCount);
    shape->annotationTextRefPoint.setNum(annotation.segmentCount);
    shape->annotationText.setNum(annotation.segmentCount);
    shape->annotationSegmentTextValid.setNum(annotation.segmentCount);
    for (int i = 0; i < annotation.segmentCount; i++) {
	const BRLObolExternalAnnotationSegment &segment =
	    annotation.segments[i];
	int kind = SoBRLVListShape::ANNOTATION_SEGMENT_NONE;
	if (segment.kind == BRLObolExternalAnnotationSegment::SEGMENT_LINE)
	    kind = SoBRLVListShape::ANNOTATION_SEGMENT_LINE;
	else if (segment.kind == BRLObolExternalAnnotationSegment::SEGMENT_TEXT)
	    kind = SoBRLVListShape::ANNOTATION_SEGMENT_TEXT;
	shape->annotationSegmentKind.set1Value(i, kind);
	shape->annotationSegmentStart.set1Value(i, segment.lineStart);
	shape->annotationSegmentEnd.set1Value(i, segment.lineEnd);
	shape->annotationTextRefPoint.set1Value(i, segment.textRefPoint);
	shape->annotationText.set1Value(i,
		(segment.text && segment.text[0]) ? segment.text : "");
	shape->annotationSegmentTextValid.set1Value(i,
		(kind == SoBRLVListShape::ANNOTATION_SEGMENT_TEXT &&
		 segment.text && segment.text[0]) ? TRUE : FALSE);
    }

    if (annotation.linePointCount > 0) {
	set_external_bounds_from_points(this, annotation.linePoints,
		annotation.linePointCount);
    } else {
	this->clearSourceBounds();
    }
    mark_external_primary_published_current(this);
    this->syncRealizedShapeOwnerState();
    return 1;
}

void
SoBRLDatabaseSource::syncRealizedShapeOwnerState(void)
{
    sync_realized_shape_owner_state_in_node(this, this);
}

SoBRLVListShape *
SoBRLDatabaseSource::getRealizedShape(void) const
{
    return this->getRealizedShape(0);
}

static SoBRLVListShape *
find_shape_in_node(SoNode *node, int &index)
{
    if (!node)
	return NULL;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	if (index == 0)
	    return static_cast<SoBRLVListShape *>(node);
	index--;
	return NULL;
    }

    if (node->isOfType(SoGroup::getClassTypeId()) &&
	    !node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SoBRLVListShape *shape = find_shape_in_node(group->getChild(i), index);
	    if (shape)
		return shape;
	}
    }

    return NULL;
}

SoBRLVListShape *
SoBRLDatabaseSource::getRealizedShape(int index) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoBRLVListShape *shape = find_shape_in_node(this->getChild(i), index);
	if (shape)
	    return shape;
    }
    return NULL;
}

static SbBool
vlist_shape_is_auxiliary(const SoBRLVListShape *shape)
{
    return shape &&
	strcmp(shape->recordRole.getValue().getString(), "auxiliary") == 0 ?
	TRUE : FALSE;
}


SoBRLVListShape *
SoBRLDatabaseSource::findAuxiliaryVListShape(const char *name) const
{
    if (!name || !name[0])
	return NULL;

    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (!node || !node->isOfType(SoBRLVListShape::getClassTypeId()))
	    continue;

	SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	if (vlist_shape_is_auxiliary(shape) &&
		strcmp(shape->geometryName.getValue().getString(), name) == 0)
	    return shape;
    }

    return NULL;
}

SoBRLDatabaseSource *
SoBRLDatabaseSource::findAuxiliarySource(const char *sourcePath) const
{
    if (!sourcePath || !sourcePath[0])
	return NULL;

    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (!node || !node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	    continue;

	SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);
	if (!source->auxiliarySource.getValue())
	    continue;

	const char *candidate = source->path.getValue().getString();
	if (strcmp(candidate, sourcePath) == 0 ||
		strcmp(database_source_skip_leading_slash(candidate),
		    database_source_skip_leading_slash(sourcePath)) == 0)
	    return source;
    }

    return NULL;
}


int
SoBRLDatabaseSource::setAuxiliaryLineSet(const char *name,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BRLObolAuxiliaryLineSetDisplayState *displayState)
{
    if (!name || !name[0] || count < 0 || (count > 0 && !points))
	return 0;

    SoBRLVListShape *shape = this->findAuxiliaryVListShape(name);
    if (count == 0) {
	if (!shape)
	    return 0;
	for (int i = 0; i < this->getNumChildren(); i++) {
	    if (this->getChild(i) == shape) {
		this->removeChild(i);
		return 1;
	    }
	}
	return 0;
    }

    if (!shape) {
	shape = new SoBRLVListShape;
	shape->setName(SbName(name));
	database_source_add_realized_child(this, shape);
    }

    std::vector<int32_t> fallbackCommands;
    if (!commands) {
	fallbackCommands.reserve(count);
	for (int i = 0; i < count; i++)
	    fallbackCommands.push_back(i == 0 ? SoBRLVListShape::MOVE :
		    SoBRLVListShape::DRAW);
	commands = fallbackCommands.data();
    }

    const char *sourcePath = this->path.getValue().getString();
    const uint32_t revision = this->sourceRevision.getValue();
    SbString identity = source_record_identity(this, sourcePath);
    if (identity.getLength() > 0)
	identity += "::";
    identity += name;

    shape->sourcePath = sourcePath ? sourcePath : "";
    shape->sourceName = name;
    shape->sourceType = "auxiliary-line-set";
    shape->sourceId = revision;
    shape->displayName = name;
    shape->geometryName = name;
    shape->sourceIdentity = identity;
    shape->cacheIdentity = record_identity_with_revision(identity.getString(),
	    revision);
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = FALSE;
    shape->drawMode = source_record_draw_mode(this);
    shape->recordRole = "auxiliary";
    shape->geometryKind = "line";
    if (displayState && displayState->valid) {
	shape->drawMode = displayState->drawMode;
	shape->visible = displayState->visible;
	shape->highlighted = displayState->highlighted;
	shape->lineStyle = displayState->lineStyle;
	shape->lineWidth = displayState->lineWidth;
	shape->transparency = displayState->transparency;
	shape->materialColorValid = displayState->materialColorValid;
	shape->materialColor = displayState->materialColor;
	shape->materialRevision = displayState->materialRevision;
	shape->colorOverride = FALSE;
	shape->color = displayState->materialColor;
    } else {
	sync_shape_display_state(shape, this);
    }
    sync_shape_placement_state(shape, this);
    shape->setLineSet(points, commands, count);
    sync_shape_owner_state(shape, this);
    return 1;
}


int
SoBRLDatabaseSource::setAuxiliarySourceLineSet(const char *sourcePath,
	const char *auxDisplayName,
	const SbVec3f *points,
	const int32_t *commands,
	int count,
	const BRLObolAuxiliaryLineSetDisplayState *displayState)
{
    if (!sourcePath || !sourcePath[0] || count < 0 ||
	    (count > 0 && !points))
	return 0;

    SoBRLDatabaseSource *source = this->findAuxiliarySource(sourcePath);
    if (count == 0) {
	if (!source)
	    return 0;
	for (int i = 0; i < this->getNumChildren(); i++) {
	    if (this->getChild(i) == source) {
		this->removeChild(i);
		return 1;
	    }
	}
	return 0;
    }

    if (!source) {
	source = new SoBRLDatabaseSource;
	source->setName(SbName(lookup_name_from_path(SbString(sourcePath))));
	database_source_add_realized_child(this, source);
    }

    const uint32_t revision = this->sourceRevision.getValue();
    source->configureDatabaseSourceInstance(sourcePath, sourcePath,
	    this->dbip, this->drawMode.getValue(), revision);
    source->auxiliarySource = TRUE;
    source->displayName = (auxDisplayName && auxDisplayName[0]) ?
	auxDisplayName : lookup_name_from_path(SbString(sourcePath));
    source->visible = this->visible.getValue();
    source->highlighted = this->highlighted.getValue();
    source->lineStyle = this->lineStyle.getValue();
    source->lineWidth = this->lineWidth.getValue();
    source->transparency = this->transparency.getValue();
    source->materialColorValid = this->materialColorValid.getValue();
    source->materialColor = this->materialColor.getValue();
    source->materialRevision = this->materialRevision.getValue();
    source->colorOverride = this->colorOverride.getValue();
    source->color = this->color.getValue();
    source->drawMatrixValid = this->drawMatrixValid.getValue();
    source->drawMatrix = this->drawMatrix.getValue();
    source->drawCenterValid = this->drawCenterValid.getValue();
    source->drawCenter = this->drawCenter.getValue();
    source->drawSizeValid = this->drawSizeValid.getValue();
    source->drawSize = this->drawSize.getValue();
    (void)remove_source_placement_transform(source);

    const char *shapeName = (auxDisplayName && auxDisplayName[0]) ?
	auxDisplayName : lookup_name_from_path(SbString(sourcePath));
    const int changed = source->setAuxiliaryLineSet(shapeName, points,
	    commands, count, displayState);
    if (changed > 0) {
	source->realizedRevision = revision;
	source->realizedSourceRevision = revision;
	source->realizedInputsRevision = source->inputsRevision.getValue();
	source->realizedViewRevision = source->viewRevision.getValue();
	source->realizationStatus = REALIZED;
	source->realizationDiagnostic = "";
	source->realizationIdentity = source_realization_identity(source);
	source->stale = FALSE;
	source->staleReason = STALE_NONE;
	source->syncRealizedShapeOwnerState();
    }
    return changed;
}


int
SoBRLDatabaseSource::clearAuxiliaryShapes(void)
{
    int removed = 0;
    for (int i = this->getNumChildren() - 1; i >= 0; i--) {
	SoNode *node = this->getChild(i);
	if (!node)
	    continue;

	if (node_is_auxiliary_source(node)) {
	    this->removeChild(i);
	    removed++;
	    continue;
	}

	if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	    SoBRLVListShape *shape = static_cast<SoBRLVListShape *>(node);
	    if (!vlist_shape_is_auxiliary(shape))
		continue;

	    this->removeChild(i);
	    removed++;
	}
    }
    return removed;
}


SoBRLMeshShape *
SoBRLDatabaseSource::getRealizedMesh(void) const
{
    return this->getRealizedMesh(0);
}

static SoBRLMeshShape *
find_mesh_in_node(SoNode *node, int &index)
{
    if (!node)
	return NULL;

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	if (index == 0)
	    return static_cast<SoBRLMeshShape *>(node);
	index--;
	return NULL;
    }

    if (node->isOfType(SoGroup::getClassTypeId()) &&
	    !node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SoBRLMeshShape *shape = find_mesh_in_node(group->getChild(i), index);
	    if (shape)
		return shape;
	}
    }

    return NULL;
}

SoBRLMeshShape *
SoBRLDatabaseSource::getRealizedMesh(int index) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoBRLMeshShape *shape = find_mesh_in_node(this->getChild(i), index);
	if (shape)
	    return shape;
    }
    return NULL;
}

SoBRLMaterialObject *
SoBRLDatabaseSource::getRealizedMaterialObject(void) const
{
    return this->getRealizedMaterialObject(0);
}

static SoBRLMaterialObject *
find_material_object_in_node(SoNode *node, int &index)
{
    if (!node)
	return NULL;

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId())) {
	if (index == 0)
	    return static_cast<SoBRLMaterialObject *>(node);
	index--;
	return NULL;
    }

    if (node->isOfType(SoGroup::getClassTypeId()) &&
	    !node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SoBRLMaterialObject *object =
		find_material_object_in_node(group->getChild(i), index);
	    if (object)
		return object;
	}
    }

    return NULL;
}

SoBRLMaterialObject *
SoBRLDatabaseSource::getRealizedMaterialObject(int index) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoBRLMaterialObject *object =
	    find_material_object_in_node(this->getChild(i), index);
	if (object)
	    return object;
    }
    return NULL;
}

static int
count_shapes_in_node(SoNode *node)
{
    if (!node)
	return 0;

    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return 1;

    int ret = 0;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    ret += count_shapes_in_node(group->getChild(i));
    }

    return ret;
}

int
SoBRLDatabaseSource::getRealizedShapeCount(void) const
{
    int ret = 0;
    for (int i = 0; i < this->getNumChildren(); i++)
	ret += count_shapes_in_node(this->getChild(i));
    return ret;
}

static int
count_meshes_in_node(SoNode *node)
{
    if (!node)
	return 0;

    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return 1;

    int ret = 0;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    ret += count_meshes_in_node(group->getChild(i));
    }

    return ret;
}

int
SoBRLDatabaseSource::getRealizedMeshCount(void) const
{
    int ret = 0;
    for (int i = 0; i < this->getNumChildren(); i++)
	ret += count_meshes_in_node(this->getChild(i));
    return ret;
}

static int
count_material_objects_in_node(SoNode *node)
{
    if (!node)
	return 0;

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return 1;

    int ret = 0;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    ret += count_material_objects_in_node(group->getChild(i));
    }

    return ret;
}

int
SoBRLDatabaseSource::getRealizedMaterialObjectCount(void) const
{
    int ret = 0;
    for (int i = 0; i < this->getNumChildren(); i++)
	ret += count_material_objects_in_node(this->getChild(i));
    return ret;
}

static void
realized_material_summary_owner(const SoBRLDatabaseSource *source,
	int sourceIndex, BRLObolRealizedMaterialSummary &summary);

static void
realized_material_summary(const SoBRLMaterialObject *object,
	BRLObolRealizedMaterialSummary &summary)
{
    summary = BRLObolRealizedMaterialSummary();
    if (!object)
	return;

    summary.valid = TRUE;
    summary.sourcePath = object->sourcePath.getValue();
    summary.sourceName = object->sourceName.getValue();
    summary.sourceType = object->sourceType.getValue();
    summary.sourceId = object->sourceId.getValue();
    summary.materialName = object->materialName.getValue();
    summary.parentName = object->parentName.getValue();
    summary.materialSource = object->materialSource.getValue();
    summary.propertyCount = object->getPropertyCount();
}

int
SoBRLDatabaseSource::getRealizedMaterialSummaryCount(void) const
{
    return this->getRealizedMaterialObjectCount();
}

SbBool
SoBRLDatabaseSource::getRealizedMaterialSummary(int index,
	BRLObolRealizedMaterialSummary &summary) const
{
    summary = BRLObolRealizedMaterialSummary();
    if (index < 0)
	return FALSE;

    SoBRLMaterialObject *object = this->getRealizedMaterialObject(index);
    if (!object)
	return FALSE;

    realized_material_summary(object, summary);
    realized_material_summary_owner(this, -1, summary);
    return TRUE;
}

SbBool
SoBRLDatabaseSource::getRealizedMaterialProperty(int materialIndex,
	int propertyIndex, SbString &groupOut, SbString &nameOut,
	SbString &valueOut) const
{
    if (materialIndex < 0 || propertyIndex < 0)
	return FALSE;

    SoBRLMaterialObject *object =
	this->getRealizedMaterialObject(materialIndex);
    if (!object)
	return FALSE;

    return object->getProperty(propertyIndex, groupOut, nameOut, valueOut);
}

static void
realized_summary_owner_instance_from_node(const SoNode *node,
	SbString &ownerSourceInstanceKey)
{
    if (!node || ownerSourceInstanceKey.getLength() > 0)
	return;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	ownerSourceInstanceKey = source_effective_instance_key(source);
	return;
    }
    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	ownerSourceInstanceKey =
	    static_cast<const SoBRLVListShape *>(node)->
	    ownerSourceInstanceKey.getValue();
	return;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	ownerSourceInstanceKey =
	    static_cast<const SoBRLMeshShape *>(node)->
	    ownerSourceInstanceKey.getValue();
	return;
    }
}

static void
realized_tree_summary_fill(const SoNode *node, int depth, SbBool hasParent,
	int ownerSourceIndex, const SbString &ownerSourcePath,
	BRLObolSceneTreeSummary &summary)
{
    summary = BRLObolSceneTreeSummary();
    if (!node)
	return;

    summary.valid = TRUE;
    summary.hasParent = hasParent;
    summary.drawTreeDepth = depth;
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    realized_summary_owner_instance_from_node(node,
	    summary.ownerSourceInstanceKey);
    summary.isGroup =
	node->isOfType(SoGroup::getClassTypeId()) ? TRUE : FALSE;
    summary.isDatabaseSource =
	node->isOfType(SoBRLDatabaseSource::getClassTypeId()) ? TRUE : FALSE;
    summary.isMaterialObject =
	node->isOfType(SoBRLMaterialObject::getClassTypeId()) ? TRUE : FALSE;
    summary.isShape =
	(node->isOfType(SoBRLVListShape::getClassTypeId()) ||
	 node->isOfType(SoBRLMeshShape::getClassTypeId())) ? TRUE : FALSE;

    if (summary.isGroup)
	summary.childCount = static_cast<const SoGroup *>(node)->getNumChildren();

    if (summary.isDatabaseSource) {
	const SoBRLDatabaseSource *sourceNode =
	    static_cast<const SoBRLDatabaseSource *>(node);
	summary.nodeKind = BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE;
	summary.path = sourceNode->path.getValue();
	summary.displayName = sourceNode->displayName.getValue();
	if (summary.ownerSourcePath.getLength() == 0)
	    summary.ownerSourcePath = sourceNode->path.getValue();
	if (summary.ownerSourceInstanceKey.getLength() == 0)
	    summary.ownerSourceInstanceKey =
		source_effective_instance_key(sourceNode);
	return;
    }

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	const SoBRLVListShape *shape =
	    static_cast<const SoBRLVListShape *>(node);
	summary.nodeKind = BRLObolSceneTreeSummary::NODE_VLIST_SHAPE;
	summary.path = shape->sourcePath.getValue();
	summary.sourceName = shape->sourceName.getValue();
	summary.sourceType = shape->sourceType.getValue();
	summary.sourceId = shape->sourceId.getValue();
	summary.displayName = shape->displayName.getValue();
	summary.geometryName = shape->geometryName.getValue();
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	const SoBRLMeshShape *shape =
	    static_cast<const SoBRLMeshShape *>(node);
	summary.nodeKind = BRLObolSceneTreeSummary::NODE_MESH_SHAPE;
	summary.path = shape->sourcePath.getValue();
	summary.sourceName = shape->sourceName.getValue();
	summary.sourceType = shape->sourceType.getValue();
	summary.sourceId = shape->sourceId.getValue();
	summary.displayName = shape->displayName.getValue();
	summary.geometryName = shape->geometryName.getValue();
	return;
    }

    if (summary.isMaterialObject) {
	const SoBRLMaterialObject *object =
	    static_cast<const SoBRLMaterialObject *>(node);
	summary.nodeKind = BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT;
	summary.path = object->sourcePath.getValue();
	summary.sourceName = object->sourceName.getValue();
	summary.sourceType = object->sourceType.getValue();
	summary.sourceId = object->sourceId.getValue();
	summary.displayName = object->materialName.getValue();
	return;
    }

    summary.nodeKind = summary.isGroup ?
	BRLObolSceneTreeSummary::NODE_GROUP :
	BRLObolSceneTreeSummary::NODE_OTHER;
}

static int
realized_tree_summary_node_count(const SoNode *node)
{
    if (!node)
	return 0;
    if (node_is_source_placement_transform(node))
	return 0;

    int count = 1;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    count += realized_tree_summary_node_count(group->getChild(i));
    }
    return count;
}

static SbBool
find_realized_tree_summary_in_node(const SoNode *node, int &index, int depth,
	SbBool hasParent, int ownerSourceIndex,
	const SbString &ownerSourcePath, BRLObolSceneTreeSummary &summary)
{
    if (!node)
	return FALSE;
    if (node_is_source_placement_transform(node))
	return FALSE;

    if (index == 0) {
	realized_tree_summary_fill(node, depth, hasParent, ownerSourceIndex,
		ownerSourcePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_tree_summary_in_node(group->getChild(i),
		    index, depth + 1, TRUE, ownerSourceIndex,
		    ownerSourcePath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedTreeSummaryCount(void) const
{
    return realized_tree_summary_node_count(this);
}

SbBool
SoBRLDatabaseSource::getRealizedTreeSummary(int index,
	BRLObolSceneTreeSummary &summary) const
{
    summary = BRLObolSceneTreeSummary();
    if (index < 0)
	return FALSE;

    const SbString ownerPath = this->path.getValue();
    const SbBool ret = find_realized_tree_summary_in_node(this, index, 0,
	    FALSE, -1, ownerPath, summary);
    if (ret && summary.ownerSourceInstanceKey.getLength() == 0)
	summary.ownerSourceInstanceKey = source_effective_instance_key(this);
    return ret;
}

static void
realized_display_summary_fill_common(BRLObolSceneDisplaySummary &summary,
	int nodeKind, int ownerSourceIndex, const SbString &ownerSourcePath,
	const SbString &nodePath)
{
    summary = BRLObolSceneDisplaySummary();
    summary.valid = TRUE;
    summary.nodeKind = nodeKind;
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    summary.path = nodePath;
}

template <typename ShapeT>
static void
realized_display_summary_fill_shape(const ShapeT *shape,
	int nodeKind, int ownerSourceIndex, const SbString &ownerSourcePath,
	BRLObolSceneDisplaySummary &summary)
{
    realized_display_summary_fill_common(summary, nodeKind, ownerSourceIndex,
	    ownerSourcePath, shape ? shape->sourcePath.getValue() : SbString(""));
    if (!shape)
	return;

    summary.ownerSourceInstanceKey = shape->ownerSourceInstanceKey.getValue();
    summary.hasDrawIntent = shape->sourcePath.getValue().getLength() > 0;
    summary.intentPath = shape->sourcePath.getValue();
    summary.intentDrawMode = shape->drawMode.getValue();
    summary.visible = shape->visible.getValue();
    summary.highlighted = shape->highlighted.getValue();
    summary.lineStyle = shape->lineStyle.getValue();
    summary.lineWidth = shape->lineWidth.getValue();
    summary.transparency = shape->transparency.getValue();
    summary.drawMode = shape->drawMode.getValue();
    summary.materialValid = TRUE;
    summary.materialRevision = shape->materialRevision.getValue();
    if (shape->materialColorValid.getValue())
	summary.materialColor = shape->materialColor.getValue();
    else if (shape->colorOverride.getValue())
	summary.materialColor = shape->color.getValue();
    summary.drawMatrixValid = shape->drawMatrixValid.getValue();
    summary.drawMatrix = shape->drawMatrix.getValue();
    summary.drawCenterValid = shape->drawCenterValid.getValue();
    summary.drawCenter = shape->drawCenter.getValue();
    summary.drawSizeValid = shape->drawSizeValid.getValue();
    summary.drawSize = shape->drawSize.getValue();
}

static void
realized_display_summary_fill(const SoNode *node, int ownerSourceIndex,
	const SbString &ownerSourcePath, BRLObolSceneDisplaySummary &summary)
{
    summary = BRLObolSceneDisplaySummary();
    if (!node)
	return;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	const SoBRLDatabaseSource *source =
	    static_cast<const SoBRLDatabaseSource *>(node);
	realized_display_summary_fill_common(summary,
		BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE,
		ownerSourceIndex, ownerSourcePath, source->path.getValue());
	summary.ownerSourceInstanceKey = source_effective_instance_key(source);
	summary.isDatabaseSource = TRUE;
	summary.hasDrawIntent = source->path.getValue().getLength() > 0;
	summary.intentPath = source->path.getValue();
	summary.intentDrawMode = source->drawMode.getValue();
	summary.visible = source->visible.getValue();
	summary.highlighted = source->highlighted.getValue();
	summary.lineStyle = source->lineStyle.getValue();
	summary.lineWidth = source->lineWidth.getValue();
	summary.transparency = source->transparency.getValue();
	summary.drawMode = source->drawMode.getValue();
	summary.materialValid = TRUE;
	summary.materialRevision = source->materialRevision.getValue();
	if (source->materialColorValid.getValue())
	    summary.materialColor = source->materialColor.getValue();
	else if (source->colorOverride.getValue())
	    summary.materialColor = source->color.getValue();
	return;
    }

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	realized_display_summary_fill_shape(
		static_cast<const SoBRLVListShape *>(node),
		BRLObolSceneTreeSummary::NODE_VLIST_SHAPE,
		ownerSourceIndex, ownerSourcePath, summary);
	return;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	realized_display_summary_fill_shape(
		static_cast<const SoBRLMeshShape *>(node),
		BRLObolSceneTreeSummary::NODE_MESH_SHAPE,
		ownerSourceIndex, ownerSourcePath, summary);
	return;
    }

    if (node->isOfType(SoBRLMaterialObject::getClassTypeId())) {
	const SoBRLMaterialObject *object =
	    static_cast<const SoBRLMaterialObject *>(node);
	realized_display_summary_fill_common(summary,
		BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT,
		ownerSourceIndex, ownerSourcePath,
		object->sourcePath.getValue());
	return;
    }

    const int nodeKind = node->isOfType(SoGroup::getClassTypeId()) ?
	BRLObolSceneTreeSummary::NODE_GROUP :
	BRLObolSceneTreeSummary::NODE_OTHER;
    realized_display_summary_fill_common(summary, nodeKind, ownerSourceIndex,
	    ownerSourcePath, "");
}

static SbBool
find_realized_display_summary_in_node(const SoNode *node, int &index,
	int ownerSourceIndex, const SbString &ownerSourcePath,
	BRLObolSceneDisplaySummary &summary)
{
    if (!node)
	return FALSE;
    if (node_is_source_placement_transform(node))
	return FALSE;

    if (index == 0) {
	realized_display_summary_fill(node, ownerSourceIndex,
		ownerSourcePath, summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_display_summary_in_node(group->getChild(i),
		    index, ownerSourceIndex, ownerSourcePath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedDisplaySummaryCount(void) const
{
    return this->getRealizedTreeSummaryCount();
}

SbBool
SoBRLDatabaseSource::getRealizedDisplaySummary(int index,
	BRLObolSceneDisplaySummary &summary) const
{
    summary = BRLObolSceneDisplaySummary();
    if (index < 0)
	return FALSE;

    const SbString ownerPath = this->path.getValue();
    const SbBool ret = find_realized_display_summary_in_node(this, index,
	    -1, ownerPath, summary);
    if (ret && summary.ownerSourceInstanceKey.getLength() == 0)
	summary.ownerSourceInstanceKey = source_effective_instance_key(this);
    return ret;
}

static void
scene_material_summary_from_display(const BRLObolSceneDisplaySummary &display,
	BRLObolSceneMaterialSummary &summary)
{
    summary = BRLObolSceneMaterialSummary();
    if (!display.valid)
	return;

    summary.valid = TRUE;
    summary.nodeKind = display.nodeKind;
    summary.materialValid =
	(display.nodeKind == BRLObolSceneTreeSummary::NODE_VLIST_SHAPE ||
	 display.nodeKind == BRLObolSceneTreeSummary::NODE_MESH_SHAPE) ?
	display.materialValid : FALSE;
    summary.materialRevision = display.materialRevision;
    summary.materialColor = display.materialColor;
    summary.ownerSourceIndex = display.ownerSourceIndex;
    summary.ownerSourcePath = display.ownerSourcePath;
    summary.ownerSourceInstanceKey = display.ownerSourceInstanceKey;
    summary.path = display.path;
}

int
SoBRLDatabaseSource::getRealizedSceneMaterialSummaryCount(void) const
{
    return this->getRealizedDisplaySummaryCount();
}

SbBool
SoBRLDatabaseSource::getRealizedSceneMaterialSummary(int index,
	BRLObolSceneMaterialSummary &summary) const
{
    summary = BRLObolSceneMaterialSummary();
    BRLObolSceneDisplaySummary display;
    if (!this->getRealizedDisplaySummary(index, display))
	return FALSE;

    scene_material_summary_from_display(display, summary);
    return TRUE;
}

static int
realized_bounds_node_kind(const SoNode *node)
{
    if (!node)
	return BRLObolSceneTreeSummary::NODE_UNKNOWN;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return BRLObolSceneTreeSummary::NODE_DATABASE_SOURCE;
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return BRLObolSceneTreeSummary::NODE_VLIST_SHAPE;
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return BRLObolSceneTreeSummary::NODE_MESH_SHAPE;
    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return BRLObolSceneTreeSummary::NODE_MATERIAL_OBJECT;
    if (node->isOfType(SoGroup::getClassTypeId()))
	return BRLObolSceneTreeSummary::NODE_GROUP;
    return BRLObolSceneTreeSummary::NODE_OTHER;
}

static SbString
realized_bounds_node_path(const SoNode *node)
{
    if (!node)
	return "";
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return static_cast<const SoBRLDatabaseSource *>(node)->path.getValue();
    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return static_cast<const SoBRLVListShape *>(node)->sourcePath.getValue();
    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return static_cast<const SoBRLMeshShape *>(node)->sourcePath.getValue();
    if (node->isOfType(SoBRLMaterialObject::getClassTypeId()))
	return static_cast<const SoBRLMaterialObject *>(node)->sourcePath.getValue();
    return "";
}

static SbBool
realized_bounds_for_vlist_shape(const SoBRLVListShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    const SoBRLVListShape *geom = shape->getGeometrySource();
    for (int i = 0; i < geom->point.getNum(); i++)
	bounds.extendBy(geom->point[i]);
    return geom->point.getNum() > 0;
}

static SbBool
realized_bounds_for_mesh_shape(const SoBRLMeshShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    for (int i = 0; i < geom->point.getNum(); i++)
	bounds.extendBy(geom->point[i]);
    return geom->point.getNum() > 0;
}

static SbBool
realized_bounds_for_node(const SoNode *node, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!node)
	return FALSE;

    SoGetBoundingBoxAction bboxAction(SbViewportRegion(1, 1));
    bboxAction.apply(const_cast<SoNode *>(node));
    bounds = bboxAction.getBoundingBox();
    if (!bounds.isEmpty())
	return TRUE;

    if (node->isOfType(SoBRLVListShape::getClassTypeId()))
	return realized_bounds_for_vlist_shape(
		static_cast<const SoBRLVListShape *>(node), bounds);

    if (node->isOfType(SoBRLMeshShape::getClassTypeId()))
	return realized_bounds_for_mesh_shape(
		static_cast<const SoBRLMeshShape *>(node), bounds);

    SbBool valid = FALSE;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    SbBox3f childBounds;
	    if (realized_bounds_for_node(group->getChild(i),
		    childBounds)) {
		bounds.extendBy(childBounds);
		valid = TRUE;
	    }
	}
    }

    if (valid)
	return TRUE;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()) &&
	    static_cast<const SoBRLDatabaseSource *>(node)->
	    getEffectiveSourceBounds(bounds))
	return TRUE;

    return valid;
}

static void
realized_bounds_summary_fill(const SoNode *node, int ownerSourceIndex,
	const SbString &ownerSourcePath, BRLObolSceneBoundsSummary &summary)
{
    summary = BRLObolSceneBoundsSummary();
    if (!node)
	return;

    summary.valid = TRUE;
    summary.nodeKind = realized_bounds_node_kind(node);
    summary.ownerSourceIndex = ownerSourceIndex;
    summary.ownerSourcePath = ownerSourcePath;
    realized_summary_owner_instance_from_node(node,
	    summary.ownerSourceInstanceKey);
    summary.path = realized_bounds_node_path(node);
    summary.boundsValid = realized_bounds_for_node(node, summary.bounds);
}

static SbBool
find_realized_bounds_summary_in_node(const SoNode *node, int &index,
	int ownerSourceIndex, const SbString &ownerSourcePath,
	BRLObolSceneBoundsSummary &summary)
{
    if (!node)
	return FALSE;
    if (node_is_source_placement_transform(node))
	return FALSE;

    if (index == 0) {
	realized_bounds_summary_fill(node, ownerSourceIndex, ownerSourcePath,
		summary);
	return TRUE;
    }
    index--;

    if (node->isOfType(SoGroup::getClassTypeId())) {
	const SoGroup *group = static_cast<const SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_bounds_summary_in_node(group->getChild(i),
		    index, ownerSourceIndex, ownerSourcePath, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedBoundsSummaryCount(void) const
{
    return this->getRealizedTreeSummaryCount();
}

SbBool
SoBRLDatabaseSource::getRealizedBoundsSummary(int index,
	BRLObolSceneBoundsSummary &summary) const
{
    summary = BRLObolSceneBoundsSummary();
    if (index < 0)
	return FALSE;

    const SbString ownerPath = this->path.getValue();
    const SbBool ret = find_realized_bounds_summary_in_node(this, index,
	    -1, ownerPath, summary);
    if (ret && summary.ownerSourceInstanceKey.getLength() == 0)
	summary.ownerSourceInstanceKey = source_effective_instance_key(this);
    return ret;
}

template <typename ShapeT>
static void
realized_shape_summary_common(const ShapeT *shape,
	BRLObolRealizedShapeSummary &summary)
{
    summary.path = shape->sourcePath.getValue();
    summary.sourceName = shape->sourceName.getValue();
    summary.sourceType = shape->sourceType.getValue();
    summary.sourceId = shape->sourceId.getValue();
    summary.displayName = shape->displayName.getValue();
    summary.geometryName = shape->geometryName.getValue();
    summary.cacheIdentity = shape->cacheIdentity.getValue();
    summary.sourceIdentity = shape->sourceIdentity.getValue();
    summary.databaseIntent = shape->databaseIntent.getValue();
    summary.overlayIntent = shape->overlayIntent.getValue();
    summary.hudIntent = shape->hudIntent.getValue();
    summary.localSource = shape->localSource.getValue();
    summary.sharedSource = shape->sharedSource.getValue();
    summary.nonDatabaseSource = shape->nonDatabaseSource.getValue();
    summary.drawMode = shape->drawMode.getValue();
    summary.recordRole = shape->recordRole.getValue();
    summary.geometryKind = shape->geometryKind.getValue();
    summary.regionId = shape->regionId.getValue();
    summary.airCode = shape->airCode.getValue();
    summary.materialId = shape->materialId.getValue();
    summary.los = shape->los.getValue();
    summary.materialColorValid = shape->materialColorValid.getValue();
    summary.materialColor = shape->materialColor.getValue();
    summary.materialShader = shape->materialShader.getValue();
    summary.materialRevision = shape->materialRevision.getValue();
    summary.visible = shape->visible.getValue();
    summary.selectable = shape->selectable.getValue();
    summary.selected = shape->selected.getValue();
    summary.highlighted = shape->highlighted.getValue();
    summary.ghosted = shape->ghosted.getValue();
    summary.hiddenLine = shape->hiddenLine.getValue();
    summary.editEmphasis = shape->editEmphasis.getValue();
    summary.editIntentId = shape->editIntentId.getValue();
    summary.editIntentRole = shape->editIntentRole.getValue();
    summary.lodPolicy = shape->lodPolicy.getValue();
}

static void
realized_shape_summary_owner(const SoBRLDatabaseSource *source,
	int sourceIndex, BRLObolRealizedShapeSummary &summary)
{
    if (!source)
	return;

    summary.ownerSourceIndex = sourceIndex;
    summary.ownerSourcePath = source->path.getValue();
    summary.ownerSourceInstanceKey = source_effective_instance_key(source);
    summary.ownerDrawMode = source->drawMode.getValue();
    summary.ownerSourceRevision = source->sourceRevision.getValue();
    summary.ownerInputsRevision = source->inputsRevision.getValue();
    summary.ownerViewRevision = source->viewRevision.getValue();
    summary.ownerRealizedRevision = source->realizedRevision.getValue();
    summary.ownerRealizedSourceRevision =
	source->realizedSourceRevision.getValue();
    summary.ownerRealizedInputsRevision =
	source->realizedInputsRevision.getValue();
    summary.ownerRealizedViewRevision =
	source->realizedViewRevision.getValue();
    summary.ownerRealizationStatus = source->realizationStatus.getValue();
    summary.ownerRealizationDiagnostic =
	source->realizationDiagnostic.getValue();
    summary.ownerRealizationIdentity =
	source->realizationIdentity.getValue();
    summary.ownerSourceStale = source->stale.getValue();
    summary.ownerStaleReason = source->staleReason.getValue();
}

static void
realized_material_summary_owner(const SoBRLDatabaseSource *source,
	int sourceIndex, BRLObolRealizedMaterialSummary &summary)
{
    if (!source)
	return;

    summary.ownerSourceIndex = sourceIndex;
    summary.ownerSourcePath = source->path.getValue();
    summary.ownerSourceInstanceKey = source_effective_instance_key(source);
    summary.ownerDrawMode = source->drawMode.getValue();
    summary.ownerSourceRevision = source->sourceRevision.getValue();
    summary.ownerInputsRevision = source->inputsRevision.getValue();
    summary.ownerViewRevision = source->viewRevision.getValue();
    summary.ownerRealizedRevision = source->realizedRevision.getValue();
    summary.ownerRealizedSourceRevision =
	source->realizedSourceRevision.getValue();
    summary.ownerRealizedInputsRevision =
	source->realizedInputsRevision.getValue();
    summary.ownerRealizedViewRevision =
	source->realizedViewRevision.getValue();
    summary.ownerRealizationStatus = source->realizationStatus.getValue();
    summary.ownerRealizationDiagnostic =
	source->realizationDiagnostic.getValue();
    summary.ownerRealizationIdentity =
	source->realizationIdentity.getValue();
    summary.ownerSourceStale = source->stale.getValue();
    summary.ownerStaleReason = source->staleReason.getValue();
}

static void
realized_shape_summary_bounds(const SoMFVec3f &points,
	BRLObolRealizedShapeSummary &summary)
{
    summary.bounds.makeEmpty();
    summary.boundsValid = FALSE;
    for (int i = 0; i < points.getNum(); i++) {
	summary.bounds.extendBy(points[i]);
	summary.boundsValid = TRUE;
    }
}

static void
realized_vlist_shape_summary(const SoBRLVListShape *shape,
	BRLObolRealizedShapeSummary &summary)
{
    summary = BRLObolRealizedShapeSummary();
    if (!shape)
	return;

    summary.valid = TRUE;
    summary.shapeKind = BRLObolRealizedShapeSummary::SHAPE_VLIST;
    realized_shape_summary_common(shape, summary);
    const SoBRLVListShape *geom = shape->getGeometrySource();
    summary.pointCount = geom->point.getNum();
    summary.commandCount = geom->command.getNum();
    summary.segmentCount = shape->getSegmentCount();
    summary.pointPrimitiveCount = shape->getPointPrimitiveCount();
    realized_shape_summary_bounds(geom->point, summary);
}

static void
realized_mesh_shape_summary(const SoBRLMeshShape *shape,
	BRLObolRealizedShapeSummary &summary)
{
    summary = BRLObolRealizedShapeSummary();
    if (!shape)
	return;

    summary.valid = TRUE;
    summary.shapeKind = BRLObolRealizedShapeSummary::SHAPE_MESH;
    realized_shape_summary_common(shape, summary);
    const SoBRLMeshShape *geom = shape->getGeometrySource();
    summary.pointCount = geom->point.getNum();
    summary.indexCount = geom->coordIndex.getNum();
    summary.triangleCount = shape->getTriangleCount();
    realized_shape_summary_bounds(geom->point, summary);
}

static SbBool
find_realized_shape_summary_in_node(SoNode *node, int &index,
	BRLObolRealizedShapeSummary &summary)
{
    if (!node)
	return FALSE;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	if (index == 0) {
	    realized_vlist_shape_summary(
		    static_cast<SoBRLVListShape *>(node), summary);
	    return TRUE;
	}
	index--;
	return FALSE;
    }

    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	if (index == 0) {
	    realized_mesh_shape_summary(
		    static_cast<SoBRLMeshShape *>(node), summary);
	    return TRUE;
	}
	index--;
	return FALSE;
    }

    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++) {
	    if (find_realized_shape_summary_in_node(group->getChild(i),
		    index, summary))
		return TRUE;
	}
    }

    return FALSE;
}

int
SoBRLDatabaseSource::getRealizedShapeSummaryCount(void) const
{
    return this->getRealizedShapeCount() + this->getRealizedMeshCount();
}

SbBool
SoBRLDatabaseSource::getRealizedShapeSummary(int index,
	BRLObolRealizedShapeSummary &summary) const
{
    summary = BRLObolRealizedShapeSummary();
    if (index < 0)
	return FALSE;

    int remaining = index;
    for (int i = 0; i < this->getNumChildren(); i++) {
	if (find_realized_shape_summary_in_node(this->getChild(i),
		remaining, summary)) {
	    realized_shape_summary_owner(this, -1, summary);
	    return TRUE;
	}
    }

    return FALSE;
}

SbBool
SoBRLDatabaseSource::getSummary(BRLObolDatabaseSourceSummary &summary) const
{
    summary = BRLObolDatabaseSourceSummary();
    summary.valid = TRUE;
    summary.path = this->path.getValue();
    summary.instanceKey = source_effective_instance_key(this);
    summary.displayName = this->displayName.getValue();
    summary.representationKey =
	source_effective_representation_key(this);
    summary.representationMode = this->representationMode.getValue();
    summary.drawMode = this->drawMode.getValue();
    summary.sourceRevision = this->sourceRevision.getValue();
    summary.inputsRevision = this->inputsRevision.getValue();
    summary.viewRevision = this->viewRevision.getValue();
    summary.realizedRevision = this->realizedRevision.getValue();
    summary.realizedSourceRevision = this->realizedSourceRevision.getValue();
    summary.realizedInputsRevision = this->realizedInputsRevision.getValue();
    summary.realizedViewRevision = this->realizedViewRevision.getValue();
    summary.realizationStatus = this->realizationStatus.getValue();
    summary.realizationDiagnostic = this->realizationDiagnostic.getValue();
    summary.realizationIdentity = this->realizationIdentity.getValue();
    summary.realizationRoleFlags = this->realizationRoleFlags.getValue();
    summary.realizationViewDependent =
	this->realizationViewDependent.getValue();
    summary.realizationViewScale = this->realizationViewScale.getValue();
    summary.realizationBotThreshold =
	this->realizationBotThreshold.getValue();
    summary.realizationCurveScale = this->realizationCurveScale.getValue();
    summary.realizationPointScale = this->realizationPointScale.getValue();
    summary.visible = this->visible.getValue();
    summary.highlighted = this->highlighted.getValue();
    summary.lineStyle = this->lineStyle.getValue();
    summary.lineWidth = this->lineWidth.getValue();
    summary.transparency = this->transparency.getValue();
    summary.materialColorValid = this->materialColorValid.getValue();
    summary.materialColor = this->materialColor.getValue();
    summary.materialRevision = this->materialRevision.getValue();
    summary.materialPolicy = this->materialPolicy.getValue();
    summary.colorOverride = this->colorOverride.getValue();
    summary.color = this->color.getValue();
    summary.drawMatrixValid = this->drawMatrixValid.getValue();
    summary.drawMatrix = this->drawMatrix.getValue();
    summary.drawCenterValid = this->drawCenterValid.getValue();
    summary.drawCenter = this->drawCenter.getValue();
    summary.drawSizeValid = this->drawSizeValid.getValue();
    summary.drawSize = this->drawSize.getValue();
    summary.sourceBoundsValid = this->getSourceBounds(
	    summary.sourceBounds);
    summary.stale = this->stale.getValue();
    summary.staleReason = this->staleReason.getValue();
    summary.realizedShapeCount = this->getRealizedShapeCount();
    summary.realizedMeshCount = this->getRealizedMeshCount();
    summary.realizedMaterialObjectCount =
	this->getRealizedMaterialObjectCount();
    return TRUE;
}
