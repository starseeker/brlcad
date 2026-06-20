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
#include "brlobol/mesh_shape.h"
#include "brlobol/vlist_shape.h"

#include "bu/list.h"
#include "nmg.h"
#include "raytrace.h"
#include "rt/func.h"
#include "rt/global.h"
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
#include <stdint.h>
#include <stdio.h>
#include <string.h>

SO_NODE_SOURCE(SoBRLDatabaseSource);

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
	uint32_t sourceId)
{
    if (!shape)
	return;

    shape->sourcePath = path ? path : "";
    shape->sourceName = sourceName ? sourceName : "";
    shape->sourceType = sourceType ? sourceType : "";
    shape->sourceId = sourceId;

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

static union tree *
make_nop_tree(void)
{
    union tree *tp = NULL;
    BU_GET(tp, union tree);
    RT_TREE_INIT(tp);
    tp->tr_op = OP_NOP;
    return tp;
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

    struct bu_list vhead;
    BU_LIST_INIT(&vhead);
    struct bg_tess_tol ttol = source_tess_tol(data->source);
    struct bn_tol tol = BN_TOL_INIT_TOL;

    int ret = rt_obj_plot(&vhead, &local_intern, &ttol, &tol);
    const char *typeLabel = primitive_type_label(&local_intern);
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
	    data->revision);
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

    if (pnts->type != RT_PNT_TYPE_PNT)
	return NULL;

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(pnts->count);
    commands.reserve(pnts->count);

    const struct pnt *point = NULL;
    for (BU_LIST_FOR(point, pnt, &(((struct pnt *)pnts->point)->l))) {
	points.push_back(SbVec3f(static_cast<float>(point->v[X]),
		static_cast<float>(point->v[Y]),
		static_cast<float>(point->v[Z])));
	commands.push_back(SoBRLVListShape::POINT);
    }

    if (points.empty() || points.size() != commands.size())
	return NULL;

    SoBRLVListShape *shape = new SoBRLVListShape;
    shape->setLineSet(points.data(), commands.data(), static_cast<int>(points.size()));
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
    SoBRLVListShape *pointShape = NULL;
    SoBRLMeshShape *shape = NULL;
    if (internalType == ID_PNTS)
	pointShape = vlist_from_pnts(static_cast<const struct rt_pnts_internal *>(local_intern.idb_ptr));
    else
	shape = mesh_from_internal(&local_intern, data->source);
    rt_db_free_internal(&local_intern);
    if (!pointShape && !shape) {
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
    if (pointShape) {
	assign_realized_identity(pointShape, tsp, path, dp->d_namep, typeLabel,
		data->revision);
    } else {
	assign_realized_identity(shape, tsp, path, dp->d_namep, typeLabel,
		data->revision);
	if (internalType == ID_BOT)
	    publish_lod_metadata_if_cached(shape, tsp->ts_dbip, dp->d_namep);
    }
    leaf->addChild(transform);
    leaf->addChild(pointShape ? static_cast<SoNode *>(pointShape) : static_cast<SoNode *>(shape));
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

    SO_NODE_ADD_FIELD(path, (""));
    SO_NODE_ADD_FIELD(drawMode, (WIREFRAME));
    SO_NODE_SET_SF_ENUM_TYPE(drawMode, DrawMode);
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
    this->stale = FALSE;
    this->staleReason = STALE_NONE;
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
    this->stale = FALSE;
    this->staleReason = STALE_NONE;
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
    this->stale = FALSE;
    this->staleReason = STALE_NONE;
    return TRUE;
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
