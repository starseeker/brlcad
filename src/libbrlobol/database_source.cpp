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
#include "nmg.h"
#include "raytrace.h"
#include "rt/func.h"
#include "rt/global.h"
#include "rt/nongeom.h"
#include "rt/db_fullpath.h"
#include "rt/tree.h"
#include "rt/vlist.h"
#include "rt/view.h"

#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/sensors/SoFieldSensor.h>

#include <vector>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

SO_NODE_SOURCE(SoBRLDatabaseSource);

static int
database_source_float_different(float a, float b)
{
    return fabsf(a - b) > 1.0e-6f;
}

BRLObolDatabaseSourceSummary::BRLObolDatabaseSourceSummary(void) :
    valid(FALSE),
    path(""),
    hasParent(FALSE),
    drawTreeDepth(0),
    parentGroupPath(""),
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
    stale(TRUE),
    staleReason(SoBRLDatabaseSource::STALE_SOURCE),
    realizedShapeCount(0),
    realizedMeshCount(0),
    realizedMaterialObjectCount(0)
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
    if (source && source->drawMode.getValue() == SoBRLDatabaseSource::SHADED)
	return BRLOBOL_LOD_DRAW_SHADED;
    return BRLOBOL_LOD_DRAW_WIRE;
}

static SbString
source_realization_identity(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";

    SbString identity;
    identity.sprintf("dbsource:%s#mode=%d;source=%u;inputs=%u;view=%u;lod=%u;abs=%.9g;rel=%.9g;norm=%.9g",
	    source->path.getValue().getString(),
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

static SbMatrix
mat_to_sbmatrix(const mat_t mat)
{
    return SbMatrix(
	    static_cast<float>(mat[0]), static_cast<float>(mat[4]), static_cast<float>(mat[8]),  static_cast<float>(mat[12]),
	    static_cast<float>(mat[1]), static_cast<float>(mat[5]), static_cast<float>(mat[9]),  static_cast<float>(mat[13]),
	    static_cast<float>(mat[2]), static_cast<float>(mat[6]), static_cast<float>(mat[10]), static_cast<float>(mat[14]),
	    static_cast<float>(mat[3]), static_cast<float>(mat[7]), static_cast<float>(mat[11]), static_cast<float>(mat[15]));
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

struct realize_walk_data {
    SoBRLDatabaseSource *source;
    uint32_t revision;
    int realized_shapes;
    int failed_shapes;
    SbString diagnostic;
};

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
    shape->displayName = sourceName ? sourceName : (path ? path : "");
    shape->geometryName = sourceName ? sourceName : "";
    shape->sourceIdentity = path ? path : "";
    shape->cacheIdentity = record_identity_with_revision(
	    path ? path : (sourceName ? sourceName : ""), sourceId);
    shape->databaseIntent = TRUE;
    shape->overlayIntent = FALSE;
    shape->hudIntent = FALSE;
    shape->localSource = FALSE;
    shape->sharedSource = FALSE;
    shape->nonDatabaseSource = FALSE;
    shape->drawMode = source_record_draw_mode(source);
    shape->recordRole = "database";
    shape->geometryKind = "";

    if (!tsp) {
	shape->regionId = 0;
	shape->airCode = 0;
	shape->materialId = 0;
	shape->los = 0;
	shape->materialColorValid = FALSE;
	shape->materialColor = SbColor(1.0f, 1.0f, 1.0f);
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
    shape->materialShader = tsp->ts_mater.ma_shader ? tsp->ts_mater.ma_shader : "";
}

template <typename ShapeT>
static void
sync_shape_owner_state(ShapeT *shape, const SoBRLDatabaseSource *source)
{
    if (!shape || !source)
	return;

    shape->ownerSourcePath = source->path.getValue();
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
    shape->ownerRealizationDiagnostic =
	source->realizationDiagnostic.getValue();
    shape->ownerRealizationIdentity = source->realizationIdentity.getValue();
    shape->ownerSourceStale = source->stale.getValue();
    shape->ownerStaleReason = source->staleReason.getValue();
}

static void
sync_realized_shape_owner_state_in_node(SoNode *node,
	const SoBRLDatabaseSource *source)
{
    if (!node || !source)
	return;

    if (node->isOfType(SoBRLVListShape::getClassTypeId())) {
	sync_shape_owner_state(static_cast<SoBRLVListShape *>(node), source);
	return;
    }
    if (node->isOfType(SoBRLMeshShape::getClassTypeId())) {
	sync_shape_owner_state(static_cast<SoBRLMeshShape *>(node), source);
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
    if (!intern || !intern->idb_ptr)
	return NULL;

    struct bu_list vhead;
    BU_LIST_INIT(&vhead);
    struct bg_tess_tol ttol = source_tess_tol(source);
    struct bn_tol tol = BN_TOL_INIT_TOL;
    int ret = rt_obj_plot(&vhead, intern, &ttol, &tol);
    if (ret < 0) {
	RT_FREE_VLIST(&rt_vlfree, &vhead);
	return NULL;
    }

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    convert_vlist(points, commands, &vhead);
    RT_FREE_VLIST(&rt_vlfree, &vhead);
    if (points.empty() || points.size() != commands.size())
	return NULL;

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->setLineSet(points.data(), commands.data(),
	    static_cast<int>(points.size()));
    return shape;
}

static union tree *
realize_leaf(struct db_tree_state *tsp,
	const struct db_full_path *pathp,
	struct rt_db_internal *UNUSED(intern),
	void *client_data)
{
    struct realize_walk_data *data = static_cast<struct realize_walk_data *>(client_data);
    if (!data || !data->source || !pathp || !tsp || !tsp->ts_dbip)
	return TREE_NULL;

    struct directory *dp = DB_FULL_PATH_CUR_DIR(pathp);
    if (!dp)
	return TREE_NULL;

    struct rt_db_internal local_intern;
    RT_DB_INTERNAL_INIT(&local_intern);
    if (rt_db_get_internal(&local_intern, dp, tsp->ts_dbip, NULL) < 0) {
	data->failed_shapes++;
	set_walk_diagnostic(data, pathp, "rt_db_get_internal failed");
	return TREE_NULL;
    }

    const char *typeLabel = primitive_type_label(&local_intern);
    if (local_intern.idb_type == ID_MATERIAL) {
	SoBRLMaterialObject *materialObject =
	    material_object_from_internal(static_cast<struct rt_material_internal *>(local_intern.idb_ptr));
	rt_db_free_internal(&local_intern);
	if (!materialObject) {
	    data->failed_shapes++;
	    set_walk_diagnostic(data, pathp,
		    "material object realization failed");
	    return TREE_NULL;
	}
	char *path = db_path_to_string(pathp);
	SoSeparator *leaf = new SoSeparator;
	assign_material_identity(materialObject, path, dp->d_namep, typeLabel,
		data->revision);
	leaf->addChild(materialObject);
	data->source->addChild(leaf);
	data->realized_shapes++;
	if (path)
	    bu_free(path, "db_path_to_string");
	return make_nop_tree();
    }

    struct bu_list vhead;
    BU_LIST_INIT(&vhead);
    struct bg_tess_tol ttol = source_tess_tol(data->source);
    struct bn_tol tol = BN_TOL_INIT_TOL;

    int ret = rt_obj_plot(&vhead, &local_intern, &ttol, &tol);
    if (ret < 0) {
	char reason[256] = {0};
	data->failed_shapes++;
	if (ret == -4) {
	    snprintf(reason, sizeof(reason),
		    "unsupported wireframe primitive type '%s'", typeLabel);
	} else {
	    snprintf(reason, sizeof(reason),
		    "wireframe plot failed for primitive type '%s'", typeLabel);
	}
	set_walk_diagnostic(data, pathp, reason);
	rt_db_free_internal(&local_intern);
	RT_FREE_VLIST(&rt_vlfree, &vhead);
	return TREE_NULL;
    }
    rt_db_free_internal(&local_intern);

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    convert_vlist(points, commands, &vhead);
    RT_FREE_VLIST(&rt_vlfree, &vhead);

    if (points.empty() || points.size() != commands.size()) {
	char reason[256] = {0};
	data->failed_shapes++;
	snprintf(reason, sizeof(reason),
		"wireframe plot produced no usable geometry for primitive type '%s'",
		typeLabel);
	set_walk_diagnostic(data, pathp, reason);
	return TREE_NULL;
    }

    char *path = db_path_to_string(pathp);
    SoSeparator *leaf = new SoSeparator;
    SoMatrixTransform *transform = new SoMatrixTransform;
    SoBRLVListShape *shape = new SoBRLVListShape;
    transform->matrix = mat_to_sbmatrix(tsp->ts_mat);
    assign_realized_identity(shape, tsp, path, dp->d_namep, typeLabel,
	    data->revision, data->source);
    shape->setLineSet(points.data(), commands.data(), static_cast<int>(points.size()));
    leaf->addChild(transform);
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
    if (!intern || !intern->idb_ptr)
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
    if (!intern || !intern->idb_ptr)
	return NULL;

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

static union tree *
realize_mesh_leaf(struct db_tree_state *tsp,
	const struct db_full_path *pathp,
	struct rt_db_internal *UNUSED(intern),
	void *client_data)
{
    struct realize_walk_data *data = static_cast<struct realize_walk_data *>(client_data);
    if (!data || !data->source || !pathp || !tsp || !tsp->ts_dbip)
	return TREE_NULL;

    struct directory *dp = DB_FULL_PATH_CUR_DIR(pathp);
    if (!dp)
	return TREE_NULL;

    struct rt_db_internal local_intern;
    RT_DB_INTERNAL_INIT(&local_intern);
    if (rt_db_get_internal(&local_intern, dp, tsp->ts_dbip, NULL) < 0) {
	data->failed_shapes++;
	set_walk_diagnostic(data, pathp, "rt_db_get_internal failed");
	return TREE_NULL;
    }

    const char *typeLabel = primitive_type_label(&local_intern);
    const int internalType = local_intern.idb_type;
    if (internalType == ID_MATERIAL) {
	SoBRLMaterialObject *materialObject =
	    material_object_from_internal(static_cast<struct rt_material_internal *>(local_intern.idb_ptr));
	rt_db_free_internal(&local_intern);
	if (!materialObject) {
	    data->failed_shapes++;
	    set_walk_diagnostic(data, pathp,
		    "material object realization failed");
	    return TREE_NULL;
	}
	char *path = db_path_to_string(pathp);
	SoSeparator *leaf = new SoSeparator;
	assign_material_identity(materialObject, path, dp->d_namep, typeLabel,
		data->revision);
	leaf->addChild(materialObject);
	data->source->addChild(leaf);
	data->realized_shapes++;
	if (path)
	    bu_free(path, "db_path_to_string");
	return make_nop_tree();
    }
    SoBRLVListShape *vlistShape = NULL;
    SoBRLMeshShape *shape = NULL;
    if (internalType == ID_PNTS)
	vlistShape = vlist_from_pnts(static_cast<const struct rt_pnts_internal *>(local_intern.idb_ptr));
    else if (internalType == ID_SKETCH || internalType == ID_ANNOT)
	vlistShape = vlist_from_plot_internal(&local_intern, data->source);
    else
	shape = mesh_from_internal(&local_intern, data->source);
    rt_db_free_internal(&local_intern);
    if (!vlistShape && !shape) {
	char reason[256] = {0};
	data->failed_shapes++;
	snprintf(reason, sizeof(reason),
		"unsupported or failed mesh conversion/tessellation for primitive type '%s'",
		typeLabel);
	set_walk_diagnostic(data, pathp, reason);
	return TREE_NULL;
    }

    char *path = db_path_to_string(pathp);
    SoSeparator *leaf = new SoSeparator;
    SoMatrixTransform *transform = new SoMatrixTransform;
    transform->matrix = mat_to_sbmatrix(tsp->ts_mat);
    if (vlistShape) {
	assign_realized_identity(vlistShape, tsp, path, dp->d_namep, typeLabel,
		data->revision, data->source);
    } else {
	assign_realized_identity(shape, tsp, path, dp->d_namep, typeLabel,
		data->revision, data->source);
	if (internalType == ID_BOT)
	    publish_lod_metadata_if_cached(shape, tsp->ts_dbip, dp->d_namep);
    }
    leaf->addChild(transform);
    leaf->addChild(vlistShape ? static_cast<SoNode *>(vlistShape) : static_cast<SoNode *>(shape));
    data->source->addChild(leaf);
    data->realized_shapes++;
    if (path)
	bu_free(path, "db_path_to_string");

    return make_nop_tree();
}

SoBRLDatabaseSource::SoBRLDatabaseSource(void) :
    dbip(NULL),
    pathSensor(NULL),
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

    SO_NODE_ADD_FIELD(path, (""));
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
    else if (sensor == source->drawModeSensor)
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
    this->stale = TRUE;
    this->staleReason = this->staleReason.getValue() | reason;
    this->realizationStatus = UNREALIZED;
    this->realizationDiagnostic = "";
    this->syncRealizedShapeOwnerState();
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
    const int validFlags = REALIZATION_ROLE_CSG | REALIZATION_ROLE_MESH;
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
SoBRLDatabaseSource::configureDatabaseSource(const char *sourcePath,
	struct db_i *database,
	int mode,
	uint32_t revision)
{
    int sanitizedMode = mode == SHADED ? SHADED : WIREFRAME;
    uint32_t reason = STALE_SOURCE;
    if (this->dbip != database)
	reason |= STALE_DATABASE;
    if (this->drawMode.getValue() != sanitizedMode)
	reason |= STALE_DRAW;

    this->detachFieldSensors();
    this->dbip = database;
    this->path = sourcePath ? sourcePath : "";
    this->drawMode = sanitizedMode;
    this->sourceRevision = revision;
    this->markStale(reason);
    this->attachFieldSensors();
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

    this->removeAllChildren();

    struct db_tree_state init_state;
    db_init_db_tree_state(&init_state, this->dbip);
    init_state.ts_stop_at_regions = 0;

    struct realize_walk_data data;
    data.source = this;
    data.revision = this->sourceRevision.getValue();
    data.realized_shapes = 0;
    data.failed_shapes = 0;

    const char *av[1] = { treeName };
    int ret = db_walk_tree(this->dbip, 1, av, 1, &init_state,
	    NULL, NULL, realize_leaf, &data);
    db_free_db_tree_state(&init_state);

    if (ret < 0 || data.realized_shapes <= 0 || data.failed_shapes > 0) {
	this->removeAllChildren();
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

    this->removeAllChildren();

    struct db_tree_state init_state;
    db_init_db_tree_state(&init_state, this->dbip);
    init_state.ts_stop_at_regions = 0;

    struct realize_walk_data data;
    data.source = this;
    data.revision = this->sourceRevision.getValue();
    data.realized_shapes = 0;
    data.failed_shapes = 0;

    const char *av[1] = { treeName };
    int ret = db_walk_tree(this->dbip, 1, av, 1, &init_state,
	    NULL, NULL, realize_mesh_leaf, &data);
    db_free_db_tree_state(&init_state);

    if (ret < 0 || data.realized_shapes <= 0 || data.failed_shapes > 0) {
	this->removeAllChildren();
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
	this->addChild(shape);
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
    shape->displayName = shape->sourceName.getValue();
    shape->geometryName = shape->sourceName.getValue();
    shape->sourceIdentity = shape->sourcePath.getValue();
    shape->cacheIdentity = record_identity_with_revision(
	    shape->sourcePath.getValue().getString(),
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

    if (node->isOfType(SoGroup::getClassTypeId())) {
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

    if (node->isOfType(SoGroup::getClassTypeId())) {
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

    if (node->isOfType(SoGroup::getClassTypeId())) {
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
	if (summary.ownerSourcePath.getLength() == 0)
	    summary.ownerSourcePath = sourceNode->path.getValue();
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
    return find_realized_tree_summary_in_node(this, index, 0, FALSE, -1,
	    ownerPath, summary);
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
    return find_realized_display_summary_in_node(this, index, -1,
	    ownerPath, summary);
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
    for (int i = 0; i < shape->point.getNum(); i++)
	bounds.extendBy(shape->point[i]);
    return shape->point.getNum() > 0;
}

static SbBool
realized_bounds_for_mesh_shape(const SoBRLMeshShape *shape, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!shape)
	return FALSE;
    for (int i = 0; i < shape->point.getNum(); i++)
	bounds.extendBy(shape->point[i]);
    return shape->point.getNum() > 0;
}

static SbBool
realized_bounds_for_node(const SoNode *node, SbBox3f &bounds)
{
    bounds.makeEmpty();
    if (!node)
	return FALSE;

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
    return find_realized_bounds_summary_in_node(this, index, -1,
	    ownerPath, summary);
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
    summary.pointCount = shape->point.getNum();
    summary.commandCount = shape->command.getNum();
    summary.segmentCount = shape->getSegmentCount();
    summary.pointPrimitiveCount = shape->getPointPrimitiveCount();
    realized_shape_summary_bounds(shape->point, summary);
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
    summary.pointCount = shape->point.getNum();
    summary.indexCount = shape->coordIndex.getNum();
    summary.triangleCount = shape->getTriangleCount();
    realized_shape_summary_bounds(shape->point, summary);
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
    summary.stale = this->stale.getValue();
    summary.staleReason = this->staleReason.getValue();
    summary.realizedShapeCount = this->getRealizedShapeCount();
    summary.realizedMeshCount = this->getRealizedMeshCount();
    summary.realizedMaterialObjectCount =
	this->getRealizedMaterialObjectCount();
    return TRUE;
}
