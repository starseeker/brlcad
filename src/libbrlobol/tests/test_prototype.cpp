/*                    T E S T _ P R O T O T Y P E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol.h"

#include "bg/line_layer.h"
#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "raytrace.h"
#include "rt/view.h"
#include "wdb.h"
#include "opennurbs_sphere.h"

#include <Inventor/SbMatrix.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/SoDB.h>
#include <Inventor/SoRenderManager.h>
#include <Inventor/SoViewport.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/annex/HUD/nodes/SoHUDLabel.h>
#include <Inventor/nodes/SoMatrixTransform.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/sensors/SoOneShotSensor.h>
#include <Inventor/sensors/SoSensorManager.h>

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

class PrototypeContextManager : public SoDB::ContextManager {
public:
    virtual void *createOffscreenContext(unsigned int UNUSED(width), unsigned int UNUSED(height))
    {
	return NULL;
    }

    virtual SbBool makeContextCurrent(void *UNUSED(context))
    {
	return FALSE;
    }

    virtual void restorePreviousContext(void *UNUSED(context))
    {
    }

    virtual void destroyContext(void *UNUSED(context))
    {
    }
};

static int
nearly_equal(float a, float b)
{
    return fabsf(a - b) < 0.0001f;
}

static SoBRLVListShape *
shape_with_path(SoBRLDatabaseSource *source, const char *path)
{
    for (int i = 0; i < source->getRealizedShapeCount(); i++) {
	SoBRLVListShape *shape = source->getRealizedShape(i);
	if (shape && BU_STR_EQUAL(shape->sourcePath.getValue().getString(), path))
	    return shape;
    }
    return NULL;
}

static int
total_segment_count(SoBRLDatabaseSource *source)
{
    int ret = 0;
    for (int i = 0; i < source->getRealizedShapeCount(); i++) {
	SoBRLVListShape *shape = source->getRealizedShape(i);
	if (shape)
	    ret += shape->getSegmentCount();
    }
    return ret;
}

static int
total_triangle_count(SoBRLDatabaseSource *source)
{
    int ret = 0;
    for (int i = 0; i < source->getRealizedMeshCount(); i++) {
	SoBRLMeshShape *mesh = source->getRealizedMesh(i);
	if (mesh)
	    ret += mesh->getTriangleCount();
    }
    return ret;
}

static int
export_path_count(const SoBRLExportAction &exportAction, const char *path)
{
    int ret = 0;
    for (int i = 0; i < exportAction.getLineCount(); i++) {
	const SoBRLExportAction::LineRecord &line = exportAction.getLine(i);
	if (BU_STR_EQUAL(line.path.getString(), path))
	    ret++;
    }
    return ret;
}

static const SoBRLExportAction::LineRecord *
export_line_with_path(const SoBRLExportAction &exportAction, const char *path)
{
    for (int i = 0; i < exportAction.getLineCount(); i++) {
	const SoBRLExportAction::LineRecord &line = exportAction.getLine(i);
	if (BU_STR_EQUAL(line.path.getString(), path))
	    return &line;
    }
    return NULL;
}

static int
export_point_path_count(const SoBRLExportAction &exportAction, const char *path)
{
    int ret = 0;
    for (int i = 0; i < exportAction.getPointCount(); i++) {
	const SoBRLExportAction::PointRecord &point = exportAction.getPoint(i);
	if (BU_STR_EQUAL(point.path.getString(), path))
	    ret++;
    }
    return ret;
}

static const SoBRLExportAction::PointRecord *
export_point_with_path(const SoBRLExportAction &exportAction, const char *path)
{
    for (int i = 0; i < exportAction.getPointCount(); i++) {
	const SoBRLExportAction::PointRecord &point = exportAction.getPoint(i);
	if (BU_STR_EQUAL(point.path.getString(), path))
	    return &point;
    }
    return NULL;
}

static int
export_triangle_path_count(const SoBRLExportAction &exportAction, const char *path)
{
    int ret = 0;
    for (int i = 0; i < exportAction.getTriangleCount(); i++) {
	const SoBRLExportAction::TriangleRecord &triangle = exportAction.getTriangle(i);
	if (BU_STR_EQUAL(triangle.path.getString(), path))
	    ret++;
    }
    return ret;
}

static const SoBRLExportAction::TriangleRecord *
export_triangle_with_path(const SoBRLExportAction &exportAction, const char *path)
{
    for (int i = 0; i < exportAction.getTriangleCount(); i++) {
	const SoBRLExportAction::TriangleRecord &triangle = exportAction.getTriangle(i);
	if (BU_STR_EQUAL(triangle.path.getString(), path))
	    return &triangle;
    }
    return NULL;
}

struct source_sensor_update_data {
    SoBRLDatabaseSource *source;
    uint32_t sourceRevision;
    int fired;
};

static void
source_sensor_update_cb(void *client_data, SoSensor *UNUSED(sensor))
{
    struct source_sensor_update_data *data = static_cast<struct source_sensor_update_data *>(client_data);
    if (!data || !data->source)
	return;

    data->fired++;
    data->source->sourceRevision = data->sourceRevision;
}

static int
shape_extents_match(SoBRLVListShape *shape,
	float xmin, float xmax,
	float ymin, float ymax,
	float zmin, float zmax)
{
    if (!shape || shape->point.getNum() <= 0)
	return 0;

    SbBox3f box;
    box.makeEmpty();
    for (int i = 0; i < shape->point.getNum(); i++)
	box.extendBy(shape->point[i]);

    if (box.isEmpty())
	return 0;

    return nearly_equal(box.getMin()[0], xmin) &&
	nearly_equal(box.getMax()[0], xmax) &&
	nearly_equal(box.getMin()[1], ymin) &&
	nearly_equal(box.getMax()[1], ymax) &&
	nearly_equal(box.getMin()[2], zmin) &&
	nearly_equal(box.getMax()[2], zmax);
}

static fastf_t *
ars_make_ring(size_t n_sides, double r, double cx, double cy, double z)
{
    fastf_t *ring = (fastf_t *)bu_malloc((n_sides + 1) * 3 * sizeof(fastf_t), "ars ring");
    for (size_t i = 0; i <= n_sides; i++) {
	double a = M_2PI * (double)i / (double)n_sides;
	ring[i*3+0] = (fastf_t)(cx + r * cos(a));
	ring[i*3+1] = (fastf_t)(cy + r * sin(a));
	ring[i*3+2] = (fastf_t)z;
    }
    return ring;
}

static fastf_t *
ars_make_cap(size_t n_sides, double cx, double cy, double z)
{
    fastf_t *ring = (fastf_t *)bu_malloc((n_sides + 1) * 3 * sizeof(fastf_t), "ars cap");
    for (size_t i = 0; i <= n_sides; i++) {
	ring[i*3+0] = (fastf_t)cx;
	ring[i*3+1] = (fastf_t)cy;
	ring[i*3+2] = (fastf_t)z;
    }
    return ring;
}

static int
make_ars_cylinder(struct rt_wdb *wdbp, const char *name)
{
    const size_t ncurves = 4;
    const size_t n_sides = 8;
    const size_t ppc = n_sides + 1;
    fastf_t **curves = (fastf_t **)bu_calloc(ncurves, sizeof(fastf_t *), "ars curves");

    curves[0] = ars_make_cap(n_sides, 120.0, 20.0, 0.0);
    curves[1] = ars_make_ring(n_sides, 5.0, 120.0, 20.0, 0.0);
    curves[2] = ars_make_ring(n_sides, 5.0, 120.0, 20.0, 10.0);
    curves[3] = ars_make_cap(n_sides, 120.0, 20.0, 10.0);

    return mk_ars(wdbp, name, ncurves, ppc, curves) == 0;
}

static void
free_local_sketch(struct rt_sketch_internal *skt)
{
    if (!skt)
	return;

    if (skt->curve.segment) {
	for (size_t i = 0; i < skt->curve.count; i++) {
	    if (skt->curve.segment[i])
		bu_free(skt->curve.segment[i], "sketch line segment");
	}
	bu_free(skt->curve.segment, "sketch segments");
    }
    if (skt->curve.reverse)
	bu_free(skt->curve.reverse, "sketch reverse flags");
    if (skt->verts)
	bu_free(skt->verts, "sketch verts");

    memset(skt, 0, sizeof(*skt));
}

static int
make_rect_sketch(struct rt_wdb *wdbp, const char *name)
{
    struct rt_sketch_internal skt;
    memset(&skt, 0, sizeof(skt));
    skt.magic = RT_SKETCH_INTERNAL_MAGIC;
    VSET(skt.V, 140.0, 20.0, 0.0);
    VSET(skt.u_vec, 1.0, 0.0, 0.0);
    VSET(skt.v_vec, 0.0, 1.0, 0.0);

    skt.vert_count = 4;
    skt.verts = (point2d_t *)bu_calloc(skt.vert_count, sizeof(point2d_t), "sketch verts");
    V2SET(skt.verts[0], 1.0, 0.0);
    V2SET(skt.verts[1], 4.0, 0.0);
    V2SET(skt.verts[2], 4.0, 3.0);
    V2SET(skt.verts[3], 1.0, 3.0);

    skt.curve.count = 4;
    skt.curve.segment = (void **)bu_calloc(skt.curve.count, sizeof(void *), "sketch segments");
    skt.curve.reverse = (int *)bu_calloc(skt.curve.count, sizeof(int), "sketch reverse flags");
    for (size_t i = 0; i < skt.curve.count; i++) {
	struct line_seg *ls;
	BU_ALLOC(ls, struct line_seg);
	ls->magic = CURVE_LSEG_MAGIC;
	ls->start = (int)i;
	ls->end = (int)((i + 1) % skt.curve.count);
	skt.curve.segment[i] = (void *)ls;
    }

    int ret = (mk_sketch(wdbp, name, &skt) == 0);
    free_local_sketch(&skt);
    return ret;
}

static void
set_poly_face(struct rt_pg_face_internal *face,
	const point_t a,
	const point_t b,
	const point_t c)
{
    vect_t ab;
    vect_t ac;
    vect_t normal;
    const fastf_t *pts[3] = {a, b, c};

    face->npts = 3;
    face->verts = (fastf_t *)bu_calloc(9, sizeof(fastf_t), "poly verts");
    face->norms = (fastf_t *)bu_calloc(9, sizeof(fastf_t), "poly norms");

    VSUB2(ab, b, a);
    VSUB2(ac, c, a);
    VCROSS(normal, ab, ac);
    VUNITIZE(normal);

    for (size_t i = 0; i < 3; i++) {
	VMOVE(&face->verts[i*3], pts[i]);
	VMOVE(&face->norms[i*3], normal);
    }
}

static int
make_poly_tet(struct rt_wdb *wdbp, const char *name)
{
    point_t p0 = {130.0, 40.0, 0.0};
    point_t p1 = {134.0, 40.0, 0.0};
    point_t p2 = {130.0, 44.0, 0.0};
    point_t p3 = {130.0, 40.0, 4.0};
    struct rt_pg_internal *pg;

    BU_ALLOC(pg, struct rt_pg_internal);
    pg->magic = RT_PG_INTERNAL_MAGIC;
    pg->npoly = 4;
    pg->poly = (struct rt_pg_face_internal *)bu_calloc(pg->npoly,
	    sizeof(struct rt_pg_face_internal), "poly faces");
    pg->max_npts = 3;

    set_poly_face(&pg->poly[0], p0, p1, p2);
    set_poly_face(&pg->poly[1], p0, p3, p1);
    set_poly_face(&pg->poly[2], p1, p3, p2);
    set_poly_face(&pg->poly[3], p2, p3, p0);

    return wdb_export(wdbp, name, (void *)pg, ID_POLY, 1.0) == 0;
}

static int
make_nmg_tet(struct rt_wdb *wdbp, const char *name)
{
    struct bn_tol tol = BN_TOL_INIT_TOL;
    struct model *m = nmg_mm();
    struct nmgregion *r = nmg_mrsv(m);
    struct shell *s = BU_LIST_FIRST(shell, &r->s_hd);
    struct vertex *verts[4] = {NULL, NULL, NULL, NULL};
    point_t pts[4] = {
	{150.0, 40.0, 0.0},
	{154.0, 40.0, 0.0},
	{150.0, 44.0, 0.0},
	{150.0, 40.0, 4.0}
    };
    struct vertex **fv0[3] = {&verts[0], &verts[1], &verts[2]};
    struct vertex **fv1[3] = {&verts[0], &verts[1], &verts[3]};
    struct vertex **fv2[3] = {&verts[0], &verts[2], &verts[3]};
    struct vertex **fv3[3] = {&verts[1], &verts[2], &verts[3]};
    struct faceuse *fu = NULL;

    if (!wdbp || !name || !m || !r || !s) {
	if (m)
	    nmg_km(m);
	return 0;
    }

    fu = nmg_cmface(s, fv0, 3);
    nmg_vertex_gv(verts[0], pts[0]);
    nmg_vertex_gv(verts[1], pts[1]);
    nmg_vertex_gv(verts[2], pts[2]);
    if (!fu || nmg_fu_planeeqn(fu, &tol) != 0) {
	nmg_km(m);
	return 0;
    }

    fu = nmg_cmface(s, fv1, 3);
    nmg_vertex_gv(verts[3], pts[3]);
    if (!fu || nmg_fu_planeeqn(fu, &tol) != 0) {
	nmg_km(m);
	return 0;
    }

    fu = nmg_cmface(s, fv2, 3);
    if (!fu || nmg_fu_planeeqn(fu, &tol) != 0) {
	nmg_km(m);
	return 0;
    }

    fu = nmg_cmface(s, fv3, 3);
    if (!fu || nmg_fu_planeeqn(fu, &tol) != 0) {
	nmg_km(m);
	return 0;
    }

    return mk_nmg(wdbp, name, m) == 0;
}

static int
make_revolve_rect(struct rt_wdb *wdbp, const char *name, const char *sketch_name)
{
    struct rt_revolve_internal *rip;

    BU_ALLOC(rip, struct rt_revolve_internal);
    rip->magic = RT_REVOLVE_INTERNAL_MAGIC;
    VSET(rip->v3d, 140.0, 20.0, 0.0);
    VSET(rip->axis3d, 0.0, 0.0, 1.0);
    V2SET(rip->v2d, 0.0, 0.0);
    V2SET(rip->axis2d, 0.0, 1.0);
    VSET(rip->r, 1.0, 0.0, 0.0);
    rip->ang = M_2PI;
    bu_vls_init(&rip->sketch_name);
    bu_vls_strcpy(&rip->sketch_name, sketch_name);
    rip->skt = NULL;

    return wdb_export(wdbp, name, (void *)rip, ID_REVOLVE, 1.0) == 0;
}

static int
make_superell_primitive(struct rt_wdb *wdbp, const char *name)
{
    struct rt_superell_internal *superell;

    BU_ALLOC(superell, struct rt_superell_internal);
    superell->magic = RT_SUPERELL_INTERNAL_MAGIC;
    VSET(superell->v, 160.0, 20.0, 3.0);
    VSET(superell->a, 4.0, 0.0, 0.0);
    VSET(superell->b, 0.0, 3.0, 0.0);
    VSET(superell->c, 0.0, 0.0, 2.0);
    superell->n = 1.0;
    superell->e = 1.0;

    return wdb_export(wdbp, name, (void *)superell, ID_SUPERELL, 1.0) == 0;
}

static int
make_script_object(struct rt_wdb *wdbp, const char *name)
{
    struct rt_script_internal *script;

    BU_ALLOC(script, struct rt_script_internal);
    script->script_magic = RT_SCRIPT_INTERNAL_MAGIC;
    bu_vls_init(&script->s_type);
    bu_vls_strcpy(&script->s_type, "test-obol-script");

    return mk_script(wdbp, name, script) == 0;
}

static int
make_constraint_object(struct rt_wdb *wdbp, const char *name)
{
    struct rt_constraint_internal *constraint;
    struct rt_db_internal intern;

    RT_DB_INTERNAL_INIT(&intern);
    BU_ALLOC(constraint, struct rt_constraint_internal);
    constraint->magic = RT_CONSTRAINT_MAGIC;
    constraint->id = 1432;
    constraint->type = 323;
    bu_vls_init(&constraint->expression);
    bu_vls_strcpy(&constraint->expression, "x=1");

    intern.idb_major_type = DB5_MAJORTYPE_BRLCAD;
    intern.idb_type = ID_CONSTRAINT;
    intern.idb_ptr = (void *)constraint;
    intern.idb_meth = &OBJ[ID_CONSTRAINT];

    return wdb_put_internal(wdbp, name, &intern, mk_conv2mm) == 0;
}

static int
make_line_annot(struct rt_wdb *wdbp, const char *name)
{
    struct rt_annot_internal annot;
    struct line_seg *ls;

    memset(&annot, 0, sizeof(annot));
    annot.magic = RT_ANNOT_INTERNAL_MAGIC;
    VSET(annot.V, 180.0, 20.0, 0.0);
    annot.vert_count = 2;
    annot.verts = (point2d_t *)bu_calloc(annot.vert_count, sizeof(point2d_t), "annot verts");
    V2SET(annot.verts[0], 0.0, 0.0);
    V2SET(annot.verts[1], 5.0, 3.0);

    annot.ant.count = 1;
    annot.ant.reverse = (int *)bu_calloc(annot.ant.count, sizeof(int), "annot reverse");
    annot.ant.segments = (void **)bu_calloc(annot.ant.count, sizeof(void *), "annot segments");
    BU_ALLOC(ls, struct line_seg);
    ls->magic = CURVE_LSEG_MAGIC;
    ls->start = 0;
    ls->end = 1;
    annot.ant.segments[0] = (void *)ls;

    int ret = (mk_annot(wdbp, name, &annot) == 0);
    bu_free(ls, "annot line segment");
    bu_free(annot.ant.segments, "annot segments");
    bu_free(annot.ant.reverse, "annot reverse");
    bu_free(annot.verts, "annot verts");
    return ret;
}

static int
make_joint_object(struct rt_wdb *wdbp, const char *name)
{
    struct rt_joint_internal *joint;

    BU_ALLOC(joint, struct rt_joint_internal);
    joint->magic = RT_JOINT_INTERNAL_MAGIC;
    VSET(joint->location, 200.0, 20.0, 0.0);
    VSET(joint->vector1, 5.0, 0.0, 0.0);
    VSET(joint->vector2, 0.0, 5.0, 0.0);
    joint->value = 1.0;
    bu_vls_init(&joint->reference_path_1);
    bu_vls_init(&joint->reference_path_2);
    bu_vls_strcpy(&joint->reference_path_1, "box.s");
    bu_vls_strcpy(&joint->reference_path_2, "ball.s");

    return wdb_export(wdbp, name, (void *)joint, ID_JOINT, 1.0) == 0;
}

static int
make_datum_object(struct rt_wdb *wdbp, const char *name)
{
    struct rt_datum_internal *datum;

    BU_ALLOC(datum, struct rt_datum_internal);
    datum->magic = RT_DATUM_INTERNAL_MAGIC;
    VSET(datum->pnt, 190.0, 40.0, 0.0);
    VSET(datum->dir, 0.0, 0.0, 10.0);
    datum->w = 0.0;
    datum->next = NULL;

    return wdb_export(wdbp, name, (void *)datum, ID_DATUM, 1.0) == 0;
}

static int
make_hf_object(struct rt_wdb *wdbp, const char *name, const char *data_path)
{
    struct rt_hf_internal *hf;

    BU_ALLOC(hf, struct rt_hf_internal);
    hf->magic = RT_HF_INTERNAL_MAGIC;
    hf->cfile[0] = '\0';
    bu_strlcpy(hf->dfile, data_path, sizeof(hf->dfile));
    bu_strlcpy(hf->fmt, "hus", sizeof(hf->fmt));
    hf->w = 3;
    hf->n = 3;
    hf->shorts = 1;
    hf->file2mm = 1.0;
    VSET(hf->v, 250.0, 20.0, 0.0);
    VSET(hf->x, 1.0, 0.0, 0.0);
    VSET(hf->y, 0.0, 1.0, 0.0);
    hf->xlen = 4.0;
    hf->ylen = 4.0;
    hf->zscale = 0.1;
    hf->mp = NULL;

    return wdb_export(wdbp, name, (void *)hf, ID_HF, 1.0) == 0;
}

static int
make_brep_sphere(struct rt_wdb *wdbp, const char *name)
{
    ON_Sphere sphere(ON_3dPoint(260.0, 20.0, 2.0), 3.0);
    ON_Brep *brep = ON_BrepSphere(sphere);
    int ret = (brep && mk_brep(wdbp, name, (void *)brep) == 0);
    delete brep;
    return ret;
}

static int
make_empty_brep(struct rt_wdb *wdbp, const char *name)
{
    ON_Brep *brep = ON_Brep::New();
    int ret = (brep && mk_brep(wdbp, name, (void *)brep) == 0);
    delete brep;
    return ret;
}

static int
make_bspline_surface(struct rt_wdb *wdbp, const char *name)
{
    struct face_g_snurb *srf = nmg_nurb_new_snurb(
	    2, 2,
	    5, 5,
	    3, 3,
	    RT_NURB_MAKE_PT_TYPE(3, RT_NURB_PT_XYZ, RT_NURB_PT_NONRAT));
    if (!srf)
	return 0;

    srf->u.knots[0] = 0.0; srf->u.knots[1] = 0.0;
    srf->u.knots[2] = 1.0;
    srf->u.knots[3] = 2.0; srf->u.knots[4] = 2.0;

    srf->v.knots[0] = 0.0; srf->v.knots[1] = 0.0;
    srf->v.knots[2] = 1.0;
    srf->v.knots[3] = 2.0; srf->v.knots[4] = 2.0;

    fastf_t *cp = srf->ctl_points;
    VSET(cp +  0, 270.0, 20.0, 0.0);
    VSET(cp +  3, 275.0, 20.0, 0.0);
    VSET(cp +  6, 280.0, 20.0, 0.0);
    VSET(cp +  9, 270.0, 25.0, 0.0);
    VSET(cp + 12, 275.0, 25.0, 2.0);
    VSET(cp + 15, 280.0, 25.0, 0.0);
    VSET(cp + 18, 270.0, 30.0, 0.0);
    VSET(cp + 21, 275.0, 30.0, 0.0);
    VSET(cp + 24, 280.0, 30.0, 0.0);

    struct face_g_snurb **surfs = (struct face_g_snurb **)bu_calloc(2,
	    sizeof(struct face_g_snurb *), "bspline surfaces");
    surfs[0] = srf;
    surfs[1] = NULL;

    return mk_bspline(wdbp, name, surfs) == 0;
}

static int
write_test_db(char *dbpath, size_t dbpath_len)
{
    FILE *fp = bu_temp_file(dbpath, dbpath_len);
    if (!fp)
	return 0;
    fclose(fp);

    struct rt_wdb *wdbp = wdb_fopen(dbpath);
    if (!wdbp)
	return 0;

    point_t min = {-2.0, -3.0, -4.0};
    point_t max = { 3.0,  4.0,  5.0};
    if (mk_rpp(wdbp, "box.s", min, max) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t sph_center = {0.0, 0.0, 0.0};
    if (mk_sph(wdbp, "ball.s", sph_center, 2.0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t ell_center = {8.0, 0.0, 0.0};
    vect_t ell_a = {2.0, 0.0, 0.0};
    vect_t ell_b = {0.0, 1.0, 0.0};
    vect_t ell_c = {0.0, 0.0, 3.0};
    if (mk_ell(wdbp, "ell.s", ell_center, ell_a, ell_b, ell_c) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t tgc_base = {0.0, 8.0, 0.0};
    vect_t tgc_height = {0.0, 0.0, 4.0};
    vect_t tgc_a = {1.5, 0.0, 0.0};
    vect_t tgc_b = {0.0, 1.0, 0.0};
    vect_t tgc_c = {0.75, 0.0, 0.0};
    vect_t tgc_d = {0.0, 0.5, 0.0};
    if (mk_tgc(wdbp, "tgc.s", tgc_base, tgc_height, tgc_a, tgc_b, tgc_c, tgc_d) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t tor_center = {16.0, 0.0, 0.0};
    vect_t tor_norm = {0.0, 0.0, 1.0};
    if (mk_tor(wdbp, "tor.s", tor_center, tor_norm, 3.0, 0.75) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t rpc_V  = {0.0, 20.0, 0.0};
    vect_t  rpc_H  = {0.0, 0.0, 10.0};
    vect_t  rpc_B  = {5.0, 0.0, 0.0};
    if (mk_rpc(wdbp, "rpc.s", rpc_V, rpc_H, rpc_B, 4.0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t rhc_V  = {20.0, 20.0, 0.0};
    vect_t  rhc_H  = {0.0, 0.0, 10.0};
    vect_t  rhc_B  = {5.0, 0.0, 0.0};
    if (mk_rhc(wdbp, "rhc.s", rhc_V, rhc_H, rhc_B, 4.0, 2.0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t epa_V  = {40.0, 20.0, 0.0};
    vect_t  epa_H  = {0.0, 0.0, 10.0};
    vect_t  epa_A  = {1.0, 0.0, 0.0};
    if (mk_epa(wdbp, "epa.s", epa_V, epa_H, epa_A, 5.0, 4.0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t ehy_V  = {60.0, 20.0, 0.0};
    vect_t  ehy_H  = {0.0, 0.0, 10.0};
    vect_t  ehy_A  = {1.0, 0.0, 0.0};
    if (mk_ehy(wdbp, "ehy.s", ehy_V, ehy_H, ehy_A, 4.0, 2.0, 1.0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t eto_V  = {0.0, 40.0, 0.0};
    vect_t  eto_N  = {0.0, 0.0, 1.0};
    vect_t  eto_C  = {8.0, 0.0, 2.0};
    if (mk_eto(wdbp, "eto.s", eto_V, eto_N, eto_C, 12.0, 3.0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t hyp_V  = {20.0, 40.0, 0.0};
    vect_t  hyp_H  = {0.0, 0.0, 10.0};
    vect_t  hyp_A  = {5.0, 0.0, 0.0};
    if (mk_hyp(wdbp, "hyp.s", hyp_V, hyp_H, hyp_A, 4.0, 0.4) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t part_V  = {40.0, 40.0, 0.0};
    vect_t  part_H  = {0.0, 0.0, 8.0};
    if (mk_particle(wdbp, "part.s", part_V, part_H, 5.0, 3.0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t cline_V = {80.0, 40.0, 0.0};
    vect_t cline_H = {0.0, 0.0, 15.0};
    if (mk_cline(wdbp, "cline.s", cline_V, cline_H, 2.0, 0.2) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    point_t hrt_V = {100.0, 20.0, 0.0};
    vect_t hrt_X = {4.0, 0.0, 0.0};
    vect_t hrt_Y = {0.0, 5.0, 0.0};
    vect_t hrt_Z = {0.0, 0.0, 6.0};
    if (mk_hrt(wdbp, "hrt.s", hrt_V, hrt_X, hrt_Y, hrt_Z, 2.0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    struct bu_list pipe_head;
    mk_pipe_init(&pipe_head);
    point_t pp0 = {100.0, 40.0, 0.0};
    point_t pp1 = {100.0, 40.0, 10.0};
    point_t pp2 = {108.0, 40.0, 16.0};
    mk_add_pipe_pnt(&pipe_head, pp0, 4.0, 2.0, 5.0);
    mk_add_pipe_pnt(&pipe_head, pp1, 4.0, 2.0, 5.0);
    mk_add_pipe_pnt(&pipe_head, pp2, 4.0, 2.0, 5.0);
    if (mk_pipe(wdbp, "pipe.s", &pipe_head) != 0) {
	mk_pipe_free(&pipe_head);
	wdb_close(wdbp);
	return 0;
    }
    mk_pipe_free(&pipe_head);

    plane_t arbn_eqn[6];
    VSET(arbn_eqn[0], 1.0, 0.0, 0.0);
    arbn_eqn[0][W] = 125.0;
    VSET(arbn_eqn[1], -1.0, 0.0, 0.0);
    arbn_eqn[1][W] = -115.0;
    VSET(arbn_eqn[2], 0.0, 1.0, 0.0);
    arbn_eqn[2][W] = 45.0;
    VSET(arbn_eqn[3], 0.0, -1.0, 0.0);
    arbn_eqn[3][W] = -35.0;
    VSET(arbn_eqn[4], 0.0, 0.0, 1.0);
    arbn_eqn[4][W] = 10.0;
    VSET(arbn_eqn[5], 0.0, 0.0, -1.0);
    arbn_eqn[5][W] = 0.0;
    if (mk_arbn(wdbp, "arbn.s", 6, (const plane_t *)arbn_eqn) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    if (!make_ars_cylinder(wdbp, "ars.s")) {
	wdb_close(wdbp);
	return 0;
    }

    if (!make_nmg_tet(wdbp, "nmg.s")) {
	wdb_close(wdbp);
	return 0;
    }

    if (!make_rect_sketch(wdbp, "rect.sketch")) {
	wdb_close(wdbp);
	return 0;
    }

    point_t extrude_V = {140.0, 20.0, 0.0};
    vect_t extrude_H = {0.0, 0.0, 6.0};
    vect_t extrude_U = {1.0, 0.0, 0.0};
    vect_t extrude_Vvec = {0.0, 1.0, 0.0};
    if (mk_extrusion(wdbp, "extrude.s", "rect.sketch", extrude_V,
	    extrude_H, extrude_U, extrude_Vvec, 0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    if (!make_revolve_rect(wdbp, "revolve.s", "rect.sketch")) {
	wdb_close(wdbp);
	return 0;
    }

    if (!make_superell_primitive(wdbp, "superell.s")) {
	wdb_close(wdbp);
	return 0;
    }

    fastf_t metaball_p0[5] = {170.0, 20.0, 0.0, 2.0, 1.0};
    fastf_t metaball_p1[5] = {176.0, 20.0, 0.0, 2.0, 1.0};
    const fastf_t *metaball_points[5] = {
	metaball_p0, metaball_p1, NULL, NULL, NULL
    };
    if (mk_metaball(wdbp, "metaball.s", 2, METABALL_ISOPOTENTIAL,
	    1.0, metaball_points) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    unsigned short dsp_data[9] = {
	0, 50, 0,
	50, 100, 50,
	0, 50, 0
    };
    if (mk_binunif(wdbp, "dsp.data", dsp_data,
	    WDB_BINUNIF_UINT16, 9) != 0) {
	wdb_close(wdbp);
	return 0;
    }
    mat_t dsp_mat;
    MAT_IDN(dsp_mat);
    MAT_DELTAS(dsp_mat, 220.0, 20.0, 0.0);
    if (mk_dsp_obj(wdbp, "dsp.s", "dsp.data", 3, 3, dsp_mat) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    unsigned char ebm_data[16] = {
	0, 0, 0, 0,
	0, 1, 1, 0,
	0, 1, 1, 0,
	0, 0, 0, 0
    };
    if (mk_binunif(wdbp, "ebm.data", ebm_data,
	    WDB_BINUNIF_UINT8, 16) != 0) {
	wdb_close(wdbp);
	return 0;
    }
    mat_t ebm_mat;
    MAT_IDN(ebm_mat);
    MAT_DELTAS(ebm_mat, 230.0, 20.0, 0.0);
    if (mk_ebm_obj(wdbp, "ebm.s", "ebm.data", 4, 4, 5.0,
	    ebm_mat) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    unsigned char vol_data[8] = {
	0, 0,
	0, 0,
	0, 1,
	0, 0
    };
    if (mk_binunif(wdbp, "vol.data", vol_data,
	    WDB_BINUNIF_UINT8, 8) != 0) {
	wdb_close(wdbp);
	return 0;
    }
    vect_t vol_cellsize = {2.0, 2.0, 2.0};
    mat_t vol_mat;
    MAT_IDN(vol_mat);
    MAT_DELTAS(vol_mat, 240.0, 20.0, 0.0);
    if (mk_vol(wdbp, "vol.s", RT_VOL_SRC_OBJ, "vol.data",
	    2, 2, 2, 1, 255, vol_cellsize, vol_mat) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    if (!make_brep_sphere(wdbp, "brep.s")) {
	wdb_close(wdbp);
	return 0;
    }
    if (!make_empty_brep(wdbp, "empty.brep")) {
	wdb_close(wdbp);
	return 0;
    }
    if (!make_bspline_surface(wdbp, "bspline.s")) {
	wdb_close(wdbp);
	return 0;
    }

    char binunif_data[4] = {1, 2, 3, 4};
    if (mk_binunif(wdbp, "payload.binunif", binunif_data,
	    WDB_BINUNIF_INT8, 4) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    if (!make_constraint_object(wdbp, "constraint.meta")) {
	wdb_close(wdbp);
	return 0;
    }

    struct bu_attribute_value_set physicalProperties;
    struct bu_attribute_value_set mechanicalProperties;
    struct bu_attribute_value_set opticalProperties;
    struct bu_attribute_value_set thermalProperties;
    bu_avs_init_empty(&physicalProperties);
    bu_avs_init_empty(&mechanicalProperties);
    bu_avs_init_empty(&opticalProperties);
    bu_avs_init_empty(&thermalProperties);
    bu_avs_add(&physicalProperties, "density", "1.0");
    if (mk_material(wdbp, "material.meta", "obol-test-material", "",
	    "test", &physicalProperties, &mechanicalProperties,
	    &opticalProperties, &thermalProperties) != 0) {
	bu_avs_free(&physicalProperties);
	bu_avs_free(&mechanicalProperties);
	bu_avs_free(&opticalProperties);
	bu_avs_free(&thermalProperties);
	wdb_close(wdbp);
	return 0;
    }
    bu_avs_free(&physicalProperties);
    bu_avs_free(&mechanicalProperties);
    bu_avs_free(&opticalProperties);
    bu_avs_free(&thermalProperties);

    if (!make_script_object(wdbp, "script.meta")) {
	wdb_close(wdbp);
	return 0;
    }

    fastf_t pnts_vertices[6] = {
	180.0, 40.0, 0.0,
	185.0, 42.0, 1.0
    };
    if (mk_pnts(wdbp, "pnts.s", RT_PNT_TYPE_PNT, 0.0, 2,
	    pnts_vertices, NULL, NULL, NULL) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    if (!make_line_annot(wdbp, "annot.s")) {
	wdb_close(wdbp);
	return 0;
    }

    if (!make_datum_object(wdbp, "datum.s")) {
	wdb_close(wdbp);
	return 0;
    }

    if (!make_joint_object(wdbp, "joint.s")) {
	wdb_close(wdbp);
	return 0;
    }

    vect_t half_norm = {0.0, 0.0, 1.0};
    if (mk_half(wdbp, "half.s", half_norm, 210.0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    if (mk_submodel(wdbp, "submodel.s", NULL, "box.s", 0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    fastf_t bot_vertices[12] = {
	0.0, 0.0, 0.0,
	2.0, 0.0, 0.0,
	0.0, 2.0, 0.0,
	0.0, 0.0, 2.0
    };
    int bot_faces[12] = {
	0, 1, 2,
	0, 3, 1,
	1, 3, 2,
	2, 3, 0
    };
    if (mk_bot(wdbp, "tet.bot", RT_BOT_SOLID, RT_BOT_UNORIENTED, 0,
	    4, 4, bot_vertices, bot_faces, NULL, NULL) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    struct wmember left;
    BU_LIST_INIT(&left.l);
    mat_t left_mat;
    MAT_IDN(left_mat);
    MAT_DELTAS(left_mat, 10.0, 0.0, 0.0);
    unsigned char left_rgb[3] = {230, 26, 51};
    if (!mk_addmember("box.s", &left.l, left_mat, WMOP_UNION) ||
	    mk_lrcomb(wdbp, "left.c", &left, 1, "plastic", "di=.8 sp=.2",
		left_rgb, 101, 7, 42, 9, 0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    struct wmember right;
    BU_LIST_INIT(&right.l);
    mat_t right_mat;
    MAT_IDN(right_mat);
    MAT_DELTAS(right_mat, 30.0, 0.0, 0.0);
    unsigned char right_rgb[3] = {26, 77, 204};
    if (!mk_addmember("box.s", &right.l, right_mat, WMOP_UNION) ||
	    mk_lrcomb(wdbp, "right.c", &right, 1, "plastic", "di=.5 sp=.1",
		right_rgb, 102, 8, 43, 10, 0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    struct wmember assembly;
    BU_LIST_INIT(&assembly.l);
    if (!mk_addmember("left.c", &assembly.l, NULL, WMOP_UNION) ||
	    !mk_addmember("right.c", &assembly.l, NULL, WMOP_UNION) ||
	    mk_lcomb(wdbp, "assembly.c", &assembly, 0, NULL, NULL, NULL, 0) != 0) {
	wdb_close(wdbp);
	return 0;
    }

    wdb_close(wdbp);
    return 1;
}

static int
write_hf_v4_test_db(char *dbpath,
	size_t dbpath_len,
	char *datapath,
	size_t datapath_len)
{
    FILE *datafp = bu_temp_file(datapath, datapath_len);
    if (!datafp)
	return 0;

    unsigned short hf_data[9] = {
	0, 10, 0,
	10, 30, 10,
	0, 10, 0
    };
    int data_ok = (fwrite(hf_data, sizeof(hf_data[0]), 9, datafp) == 9);
    fclose(datafp);
    if (!data_ok) {
	bu_file_delete(datapath);
	return 0;
    }

    FILE *fp = bu_temp_file(dbpath, dbpath_len);
    if (!fp) {
	bu_file_delete(datapath);
	return 0;
    }
    fclose(fp);

    struct rt_wdb *wdbp = wdb_fopen_v(dbpath, 4);
    if (!wdbp) {
	bu_file_delete(datapath);
	return 0;
    }

    if (!make_hf_object(wdbp, "hf.s", datapath)) {
	wdb_close(wdbp);
	bu_file_delete(datapath);
	return 0;
    }

    wdb_close(wdbp);
    return 1;
}

static int
write_poly_v4_test_db(char *dbpath, size_t dbpath_len)
{
    FILE *fp = bu_temp_file(dbpath, dbpath_len);
    if (!fp)
	return 0;
    fclose(fp);

    struct rt_wdb *wdbp = wdb_fopen_v(dbpath, 4);
    if (!wdbp)
	return 0;

    if (!make_poly_tet(wdbp, "poly.s")) {
	wdb_close(wdbp);
	return 0;
    }

    wdb_close(wdbp);
    return 1;
}

static int
open_database(const char *dbpath, struct db_i **dbipp)
{
    struct db_i *dbip = db_open(dbpath, DB_OPEN_READONLY);
    if (!dbip)
	return 0;
    if (db_dirbuild(dbip) < 0) {
	db_close(dbip);
	return 0;
    }
    *dbipp = dbip;
    return 1;
}

static int
exercise_required_hierarchy_model(const char *model_file,
	const char *root_name,
	int min_wire_shapes,
	int min_wire_segments,
	int min_mesh_shapes,
	int min_mesh_triangles,
	const SbViewportRegion &viewport)
{
    char dbpath[MAXPATHLEN] = {0};
    const char *brlcadRoot = getenv("BRLCAD_ROOT");
    if (brlcadRoot) {
	snprintf(dbpath, MAXPATHLEN, "%s/share/db/%s", brlcadRoot, model_file);
    } else {
	snprintf(dbpath, MAXPATHLEN, "share/db/%s", model_file);
    }

    if (!bu_file_exists(dbpath, NULL)) {
	fprintf(stderr, "required Obol hierarchy model is missing: %s\n", dbpath);
	return 0;
    }

    struct db_i *dbip = NULL;
    if (!open_database(dbpath, &dbip)) {
	fprintf(stderr, "failed to open required Obol hierarchy model: %s\n", dbpath);
	return 0;
    }

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = root_name;
    source->sourceRevision = 101;
    root->addChild(source);

    SoBRLRealizeAction hierarchyRealize;
    hierarchyRealize.apply(root);
    if (hierarchyRealize.getRealizedSourceCount() != 1 ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED) {
	fprintf(stderr, "%s wire hierarchy realization failed: %s\n",
		model_file,
		hierarchyRealize.getDiagnostics().getLength() > 0 ?
		hierarchyRealize.getDiagnostics().getString() :
		source->realizationDiagnostic.getValue().getString());
	root->unref();
	db_close(dbip);
	return 0;
    }

    int wireShapeCount = source->getRealizedShapeCount();
    int wireSegmentCount = total_segment_count(source);
    if (wireShapeCount < min_wire_shapes || wireSegmentCount < min_wire_segments) {
	fprintf(stderr, "%s wire realization was too small: shapes=%d segments=%d\n",
		model_file, wireShapeCount, wireSegmentCount);
	root->unref();
	db_close(dbip);
	return 0;
    }

    SoGetBoundingBoxAction hierarchyBBoxAction(viewport);
    hierarchyBBoxAction.apply(root);
    SbBox3f bbox = hierarchyBBoxAction.getBoundingBox();
    if (bbox.isEmpty()) {
	fprintf(stderr, "%s wire hierarchy produced an empty bounding box\n", model_file);
	root->unref();
	db_close(dbip);
	return 0;
    }

    SoBRLExportAction hierarchyExport;
    hierarchyExport.apply(root);
    if (hierarchyExport.getLineCount() != wireSegmentCount ||
	    hierarchyExport.getBounds().isEmpty()) {
	fprintf(stderr, "%s wire hierarchy export did not match realized line geometry\n",
		model_file);
	root->unref();
	db_close(dbip);
	return 0;
    }

    if (min_mesh_shapes > 0 || min_mesh_triangles > 0) {
	source->drawMode = SoBRLDatabaseSource::SHADED;
	source->sourceRevision = 102;
	source->tessellationRelTol = 0.2f;
	if (!source->needsRealization()) {
	    fprintf(stderr, "%s shaded mode change did not invalidate realization\n",
		    model_file);
	    root->unref();
	    db_close(dbip);
	    return 0;
	}

	SoBRLRealizeAction hierarchyMeshRealize;
	hierarchyMeshRealize.apply(root);
	if (hierarchyMeshRealize.getRealizedSourceCount() != 1 ||
		source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED) {
	    fprintf(stderr, "%s shaded hierarchy realization failed: %s\n",
		    model_file,
		    hierarchyMeshRealize.getDiagnostics().getLength() > 0 ?
		    hierarchyMeshRealize.getDiagnostics().getString() :
		    source->realizationDiagnostic.getValue().getString());
	    root->unref();
	    db_close(dbip);
	    return 0;
	}
	if (source->getRealizedMeshCount() < min_mesh_shapes) {
	    fprintf(stderr, "%s shaded realization produced too few mesh leaves: %d\n",
		    model_file, source->getRealizedMeshCount());
	    root->unref();
	    db_close(dbip);
	    return 0;
	}

	int hierarchyTriangleCount = total_triangle_count(source);
	if (hierarchyTriangleCount < min_mesh_triangles) {
	    fprintf(stderr, "%s shaded realization produced too few triangles: %d\n",
		    model_file, hierarchyTriangleCount);
	    root->unref();
	    db_close(dbip);
	    return 0;
	}

	SoBRLExportAction hierarchyMeshExport;
	hierarchyMeshExport.apply(root);
	if (hierarchyMeshExport.getLineCount() != 0 ||
		hierarchyMeshExport.getTriangleCount() != hierarchyTriangleCount ||
		hierarchyMeshExport.getBounds().isEmpty()) {
	    fprintf(stderr, "%s shaded export did not match realized mesh geometry\n",
		    model_file);
	    root->unref();
	    db_close(dbip);
	    return 0;
	}

	SoBRLMeasureAction hierarchyMeshMeasure;
	hierarchyMeshMeasure.apply(root);
	if (!hierarchyMeshMeasure.hasFaces() ||
		hierarchyMeshMeasure.getShapeCount() != source->getRealizedMeshCount() ||
		hierarchyMeshMeasure.getTriangleCount() != hierarchyTriangleCount ||
		hierarchyMeshMeasure.getSurfaceArea() <= 0.0f) {
	    fprintf(stderr, "%s shaded measure did not report mesh face metrics\n",
		    model_file);
	    root->unref();
	    db_close(dbip);
	    return 0;
	}
    }

    root->unref();
    db_close(dbip);
    return 1;
}

static int
exercise_generated_primitive(struct db_i *dbip,
	const char *name,
	int min_wire_segments,
	int min_mesh_triangles,
	const SbViewportRegion &viewport,
	float shaded_rel_tol = 0.15f)
{
    if (!dbip || !name)
	return 0;

    SbString fullPath;
    fullPath.sprintf("/%s", name);

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = name;
    source->sourceRevision = 201;
    root->addChild(source);

    SoBRLRealizeAction wireRealize;
    wireRealize.apply(root);
    if (wireRealize.getRealizedSourceCount() != 1 ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED) {
	fprintf(stderr, "%s wire realization failed: %s\n",
		name,
		wireRealize.getDiagnostics().getLength() > 0 ?
		wireRealize.getDiagnostics().getString() :
		source->realizationDiagnostic.getValue().getString());
	root->unref();
	return 0;
    }

    if (source->getRealizedShapeCount() != 1) {
	fprintf(stderr, "%s wire realization produced %d shape leaves\n",
		name, source->getRealizedShapeCount());
	root->unref();
	return 0;
    }

    SoBRLVListShape *shape = source->getRealizedShape();
    if (!shape ||
	    shape->getSegmentCount() < min_wire_segments ||
	    strcmp(shape->sourcePath.getValue().getString(), fullPath.getString()) != 0) {
	fprintf(stderr, "%s wire shape did not meet segment/path expectations\n", name);
	root->unref();
	return 0;
    }

    SoGetBoundingBoxAction wireBBoxAction(viewport);
    wireBBoxAction.apply(root);
    if (wireBBoxAction.getBoundingBox().isEmpty()) {
	fprintf(stderr, "%s wire realization produced an empty bounding box\n", name);
	root->unref();
	return 0;
    }

    SoBRLExportAction wireExport;
    wireExport.apply(root);
    if (wireExport.getLineCount() != shape->getSegmentCount() ||
	    export_path_count(wireExport, fullPath.getString()) != shape->getSegmentCount() ||
	    wireExport.getBounds().isEmpty()) {
	fprintf(stderr, "%s wire export did not preserve line geometry/path identity\n", name);
	root->unref();
	return 0;
    }

    SoBRLMeasureAction wireMeasure;
    wireMeasure.apply(root);
    if (!wireMeasure.hasSegments() ||
	    wireMeasure.getShapeCount() != 1 ||
	    wireMeasure.getSegmentCount() != shape->getSegmentCount() ||
	    wireMeasure.getTotalLength() <= 0.0f ||
	    wireMeasure.getBounds().isEmpty()) {
	fprintf(stderr, "%s wire measure did not report line metrics\n", name);
	root->unref();
	return 0;
    }

    if (min_mesh_triangles > 0) {
	source->drawMode = SoBRLDatabaseSource::SHADED;
	source->sourceRevision = 202;
	source->tessellationRelTol = shaded_rel_tol;
	if (!source->needsRealization()) {
	    fprintf(stderr, "%s shaded mode change did not invalidate realization\n", name);
	    root->unref();
	    return 0;
	}

	SoBRLRealizeAction meshRealize;
	meshRealize.apply(root);
	if (meshRealize.getRealizedSourceCount() != 1 ||
		source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED) {
	    fprintf(stderr, "%s shaded realization failed: %s\n",
		    name,
		    meshRealize.getDiagnostics().getLength() > 0 ?
		    meshRealize.getDiagnostics().getString() :
		    source->realizationDiagnostic.getValue().getString());
	    root->unref();
	    return 0;
	}
	if (source->getRealizedShapeCount() != 0 ||
		source->getRealizedMeshCount() != 1) {
	    fprintf(stderr, "%s shaded realization did not replace wire leaves with one mesh\n", name);
	    root->unref();
	    return 0;
	}

	SoBRLMeshShape *mesh = source->getRealizedMesh();
	if (!mesh ||
		mesh->getTriangleCount() < min_mesh_triangles ||
		strcmp(mesh->sourcePath.getValue().getString(), fullPath.getString()) != 0) {
	    fprintf(stderr, "%s mesh did not meet triangle/path expectations\n", name);
	    root->unref();
	    return 0;
	}

	SoBRLExportAction meshExport;
	meshExport.apply(root);
	if (meshExport.getLineCount() != 0 ||
		meshExport.getTriangleCount() != mesh->getTriangleCount() ||
		export_triangle_path_count(meshExport, fullPath.getString()) != mesh->getTriangleCount() ||
		meshExport.getBounds().isEmpty()) {
	    fprintf(stderr, "%s mesh export did not preserve triangle geometry/path identity\n", name);
	    root->unref();
	    return 0;
	}

	SoBRLMeasureAction meshMeasure;
	meshMeasure.apply(root);
	if (!meshMeasure.hasFaces() ||
		meshMeasure.getShapeCount() != 1 ||
		meshMeasure.getTriangleCount() != mesh->getTriangleCount() ||
		meshMeasure.getSurfaceArea() <= 0.0f ||
		meshMeasure.getBounds().isEmpty()) {
	    fprintf(stderr, "%s mesh measure did not report face metrics\n", name);
	    root->unref();
	    return 0;
	}
    }

    root->unref();
    return 1;
}

static int
exercise_generated_primitive_wire_diagnostic(struct db_i *dbip,
	const char *name,
	const char *type_name,
	const SbViewportRegion &viewport)
{
    if (!dbip || !name || !type_name)
	return 0;

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = name;
    source->sourceRevision = 301;
    root->addChild(source);

    SoBRLRealizeAction wireRealize;
    wireRealize.apply(root);
    if (wireRealize.getFailedSourceCount() != 1 ||
	    wireRealize.getRealizedSourceCount() != 0 ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::FAILED ||
	    source->getRealizedMeshCount() != 0 ||
	    source->getRealizedShapeCount() != 0 ||
	    source->realizationDiagnostic.getValue().getLength() == 0 ||
	    !strstr(source->realizationDiagnostic.getValue().getString(), name) ||
	    !strstr(source->realizationDiagnostic.getValue().getString(), type_name) ||
	    wireRealize.getDiagnostics().getLength() == 0 ||
	    !strstr(wireRealize.getDiagnostics().getString(), type_name)) {
	fprintf(stderr, "%s wire diagnostic was not explicit: %s\n",
		name, source->realizationDiagnostic.getValue().getString());
	root->unref();
	return 0;
    }

    SoGetBoundingBoxAction bboxAction(viewport);
    bboxAction.apply(root);
    if (!bboxAction.getBoundingBox().isEmpty()) {
	fprintf(stderr, "%s failed wire realization left drawable scene content\n", name);
	root->unref();
	return 0;
    }

    root->unref();
    return 1;
}

static int
exercise_generated_primitive_mesh_diagnostic(struct db_i *dbip,
	const char *name,
	const char *type_name,
	const SbViewportRegion &viewport)
{
    if (!dbip || !name || !type_name)
	return 0;

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = name;
    source->drawMode = SoBRLDatabaseSource::SHADED;
    source->sourceRevision = 301;
    root->addChild(source);

    SoBRLRealizeAction meshRealize;
    meshRealize.apply(root);
    if (meshRealize.getFailedSourceCount() != 1 ||
	    meshRealize.getRealizedSourceCount() != 0 ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::FAILED ||
	    source->getRealizedMeshCount() != 0 ||
	    source->getRealizedShapeCount() != 0 ||
	    source->realizationDiagnostic.getValue().getLength() == 0 ||
	    !strstr(source->realizationDiagnostic.getValue().getString(), name) ||
	    !strstr(source->realizationDiagnostic.getValue().getString(), type_name) ||
	    meshRealize.getDiagnostics().getLength() == 0 ||
	    !strstr(meshRealize.getDiagnostics().getString(), type_name)) {
	fprintf(stderr, "%s shaded diagnostic was not explicit: %s\n",
		name, source->realizationDiagnostic.getValue().getString());
	root->unref();
	return 0;
    }

    SoGetBoundingBoxAction bboxAction(viewport);
    bboxAction.apply(root);
    if (!bboxAction.getBoundingBox().isEmpty()) {
	fprintf(stderr, "%s failed shaded realization left drawable scene content\n", name);
	root->unref();
	return 0;
    }

    root->unref();
    return 1;
}

static int
exercise_generated_pnts_shaded_points(struct db_i *dbip,
	const char *name,
	const SbViewportRegion &viewport)
{
    if (!dbip || !name)
	return 0;

    SbString fullPath;
    fullPath.sprintf("/%s", name);

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = name;
    source->drawMode = SoBRLDatabaseSource::SHADED;
    source->sourceRevision = 302;
    root->addChild(source);

    SoBRLRealizeAction pointRealize;
    pointRealize.apply(root);
    if (pointRealize.getRealizedSourceCount() != 1 ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED) {
	fprintf(stderr, "%s shaded point realization failed: %s\n",
		name,
		pointRealize.getDiagnostics().getLength() > 0 ?
		pointRealize.getDiagnostics().getString() :
		source->realizationDiagnostic.getValue().getString());
	root->unref();
	return 0;
    }

    if (source->getRealizedShapeCount() != 1 ||
	    source->getRealizedMeshCount() != 0) {
	fprintf(stderr, "%s shaded point realization produced shape/mesh counts %d/%d\n",
		name, source->getRealizedShapeCount(), source->getRealizedMeshCount());
	root->unref();
	return 0;
    }

    SoBRLVListShape *shape = source->getRealizedShape();
    if (!shape ||
	    shape->getSegmentCount() != 0 ||
	    shape->getPointPrimitiveCount() != 2 ||
	    !BU_STR_EQUAL(shape->sourcePath.getValue().getString(), fullPath.getString()) ||
	    !BU_STR_EQUAL(shape->sourceType.getValue().getString(), "pnts")) {
	fprintf(stderr, "%s shaded point shape did not meet point/path/type expectations\n", name);
	root->unref();
	return 0;
    }

    int primitiveIndex = -1;
    SbVec3f point;
    if (!shape->getPointPrimitive(0, primitiveIndex, point) ||
	    primitiveIndex != 0 ||
	    !nearly_equal(point[0], 180.0f) ||
	    !nearly_equal(point[1], 40.0f) ||
	    !nearly_equal(point[2], 0.0f)) {
	fprintf(stderr, "%s first shaded point primitive is wrong\n", name);
	root->unref();
	return 0;
    }
    if (!shape->getPointPrimitive(1, primitiveIndex, point) ||
	    primitiveIndex != 1 ||
	    !nearly_equal(point[0], 185.0f) ||
	    !nearly_equal(point[1], 42.0f) ||
	    !nearly_equal(point[2], 1.0f)) {
	fprintf(stderr, "%s second shaded point primitive is wrong\n", name);
	root->unref();
	return 0;
    }

    SoGetBoundingBoxAction bboxAction(viewport);
    bboxAction.apply(root);
    SbBox3f bbox = bboxAction.getBoundingBox();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 180.0f) ||
	    !nearly_equal(bbox.getMin()[1], 40.0f) ||
	    !nearly_equal(bbox.getMin()[2], 0.0f) ||
	    !nearly_equal(bbox.getMax()[0], 185.0f) ||
	    !nearly_equal(bbox.getMax()[1], 42.0f) ||
	    !nearly_equal(bbox.getMax()[2], 1.0f)) {
	fprintf(stderr, "%s shaded point bounding box is wrong\n", name);
	root->unref();
	return 0;
    }

    SoRayPickAction pickAction(viewport);
    pickAction.setRay(SbVec3f(180.0f, 40.0f, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    pickAction.setRadius(15.0f);
    pickAction.apply(root);
    const SoPickedPoint *pickedPoint = pickAction.getPickedPoint();
    const SoDetail *rawDetail = pickedPoint ? pickedPoint->getDetail(shape) : NULL;
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId())) {
	fprintf(stderr, "%s shaded point pick did not return BRL-CAD detail\n", name);
	root->unref();
	return 0;
    }
    const SoBRLPickDetail *pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (pickDetail->getPrimitiveKind() != SoBRLPickDetail::POINT ||
	    pickDetail->getPrimitiveIndex() != 0 ||
	    !BU_STR_EQUAL(pickDetail->getPath().getString(), fullPath.getString()) ||
	    !BU_STR_EQUAL(pickDetail->getSourceType().getString(), "pnts") ||
	    !nearly_equal(pickDetail->getModelPoint()[0], 180.0f) ||
	    !nearly_equal(pickDetail->getModelPoint()[1], 40.0f) ||
	    !nearly_equal(pickDetail->getModelPoint()[2], 0.0f)) {
	fprintf(stderr, "%s shaded point pick detail is wrong\n", name);
	root->unref();
	return 0;
    }

    SoBRLExportAction pointExport;
    pointExport.apply(root);
    if (pointExport.getLineCount() != 0 ||
	    pointExport.getTriangleCount() != 0 ||
	    pointExport.getPointCount() != 2 ||
	    export_point_path_count(pointExport, fullPath.getString()) != 2 ||
	    pointExport.getBounds().isEmpty()) {
	fprintf(stderr, "%s shaded point export did not preserve point geometry/path identity\n", name);
	root->unref();
	return 0;
    }
    const SoBRLExportAction::PointRecord *firstPoint =
	export_point_with_path(pointExport, fullPath.getString());
    const SoBRLExportAction::PointRecord &secondPoint = pointExport.getPoint(1);
    if (!firstPoint ||
	    firstPoint->primitiveIndex != 0 ||
	    !BU_STR_EQUAL(firstPoint->sourceType.getString(), "pnts") ||
	    !nearly_equal(firstPoint->point[0], 180.0f) ||
	    !nearly_equal(firstPoint->point[1], 40.0f) ||
	    !nearly_equal(firstPoint->point[2], 0.0f) ||
	    secondPoint.primitiveIndex != 1 ||
	    !nearly_equal(secondPoint.point[0], 185.0f) ||
	    !nearly_equal(secondPoint.point[1], 42.0f) ||
	    !nearly_equal(secondPoint.point[2], 1.0f)) {
	fprintf(stderr, "%s shaded point export records are wrong\n", name);
	root->unref();
	return 0;
    }

    SoBRLSnapAction pointSnap;
    pointSnap.setEnabledKinds(SoBRLSnapAction::ENDPOINT);
    pointSnap.setQueryPoint(SbVec3f(180.1f, 40.0f, 0.0f));
    pointSnap.setTolerance(0.25f);
    pointSnap.apply(root);
    if (!pointSnap.hasCandidate() ||
	    pointSnap.getKind() != SoBRLSnapAction::ENDPOINT ||
	    pointSnap.getPrimitiveIndex() != 0 ||
	    !BU_STR_EQUAL(pointSnap.getPath().getString(), fullPath.getString()) ||
	    !nearly_equal(pointSnap.getPoint()[0], 180.0f) ||
	    !nearly_equal(pointSnap.getPoint()[1], 40.0f) ||
	    !nearly_equal(pointSnap.getPoint()[2], 0.0f)) {
	fprintf(stderr, "%s shaded point snap did not preserve endpoint identity\n", name);
	root->unref();
	return 0;
    }

    root->unref();
    return 1;
}

int
main(int UNUSED(argc), const char **UNUSED(argv))
{
    brlobol_init(NULL);
    PrototypeContextManager prototypeManager;
    brlobol_init(&prototypeManager);
    if (SoDB::getContextManager() != &prototypeManager)
	FAIL("brlobol_init should update the context manager after null initialization");

    SoSeparator *root = new SoSeparator;
    root->ref();

    SoBRLDatabaseSource *source = new SoBRLDatabaseSource;
    source->path = "/prototype/arb8";
    source->sourceRevision = 3;
    root->addChild(source);

    if (!source->needsRealization())
	FAIL("new source should start stale");

    SoBRLRealizeAction realize;
    realize.apply(root);
    if (realize.getVisitedSourceCount() != 1)
	FAIL("realize action should visit one database source");
    if (realize.getRealizedSourceCount() != 1)
	FAIL("realize action should realize the stale source");
    if (source->needsRealization())
	FAIL("source should be clean after realization");
    if (source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED)
	FAIL("source realization status should be REALIZED");

    SoBRLVListShape *shape = source->getRealizedShape();
    if (!shape)
	FAIL("source should own a realized Obol vlist shape child");
    if (shape->getSegmentCount() != 4)
	FAIL("prototype wireframe should have four line segments");
    if (!shape->isOfType(SoBRLVListShape::getClassTypeId()))
	FAIL("realized child should have native Obol type identity");

    const float halfExtent = 1.75f;
    SbViewportRegion viewport(100, 100);
    SoGetBoundingBoxAction bboxAction(viewport);
    bboxAction.apply(root);
    SbBox3f bbox = bboxAction.getBoundingBox();
    if (bbox.isEmpty())
	FAIL("realized source should contribute a bounding box");
    if (!nearly_equal(bbox.getMin()[0], -halfExtent) ||
	    !nearly_equal(bbox.getMax()[0], halfExtent) ||
	    !nearly_equal(bbox.getMin()[1], -halfExtent) ||
	    !nearly_equal(bbox.getMax()[1], halfExtent))
	FAIL("bounding box should come from realized Obol geometry");

    SoRayPickAction pickAction(viewport);
    pickAction.setRay(SbVec3f(0.0f, -halfExtent, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    pickAction.apply(root);
    const SoPickedPoint *pickedPoint = pickAction.getPickedPoint();
    if (!pickedPoint)
	FAIL("ray pick should hit realized vlist shape");

    const SoDetail *rawDetail = pickedPoint->getDetail(shape);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("pick should return a BRL-CAD Obol pick detail");

    const SoBRLPickDetail *pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(), "/prototype/arb8") != 0)
	FAIL("pick detail should preserve database path identity");
    if (pickDetail->getPrimitiveKind() != SoBRLPickDetail::LINE_SEGMENT)
	FAIL("pick detail should identify a line segment");

    SoBRLSnapAction snapAction;
    snapAction.setQueryPoint(SbVec3f(halfExtent, 0.2f, 0.0f));
    snapAction.setTolerance(0.3f);
    snapAction.apply(root);
    if (!snapAction.hasCandidate())
	FAIL("snap action should find a vlist line candidate");
    if (snapAction.getKind() != SoBRLSnapAction::LINE_NEAREST)
	FAIL("snap action should prefer nearest point on line for this query");
    if (strcmp(snapAction.getPath().getString(), "/prototype/arb8") != 0)
	FAIL("snap action should preserve database path identity");
    if (!nearly_equal(snapAction.getPoint()[0], halfExtent) ||
	    !nearly_equal(snapAction.getPoint()[1], 0.2f))
	FAIL("snap point should lie on the realized vertical segment");

    SoBRLSnapAction midpointSnap;
    midpointSnap.setEnabledKinds(SoBRLSnapAction::MIDPOINT);
    midpointSnap.setQueryPoint(SbVec3f(0.0f, -halfExtent, 0.0f));
    midpointSnap.setTolerance(0.1f);
    midpointSnap.apply(root);
    if (!midpointSnap.hasCandidate() ||
	    midpointSnap.getKind() != SoBRLSnapAction::MIDPOINT ||
	    strcmp(midpointSnap.getPath().getString(), "/prototype/arb8") != 0 ||
	    !nearly_equal(midpointSnap.getPoint()[0], 0.0f) ||
	    !nearly_equal(midpointSnap.getPoint()[1], -halfExtent))
	FAIL("snap action should support explicit midpoint-only policy");

    SoBRLSnapAction centerSnap;
    centerSnap.setEnabledKinds(SoBRLSnapAction::CENTER);
    centerSnap.setQueryPoint(SbVec3f(0.0f, 0.0f, 0.0f));
    centerSnap.setTolerance(0.1f);
    centerSnap.apply(root);
    if (!centerSnap.hasCandidate() ||
	    centerSnap.getKind() != SoBRLSnapAction::CENTER ||
	    strcmp(centerSnap.getPath().getString(), "/prototype/arb8") != 0 ||
	    !nearly_equal(centerSnap.getPoint()[0], 0.0f) ||
	    !nearly_equal(centerSnap.getPoint()[1], 0.0f))
	FAIL("snap action should support explicit center-only policy");

    SoBRLSnapAction planeSnap;
    planeSnap.setEnabledKinds(SoBRLSnapAction::CONSTRUCTION_PLANE);
    planeSnap.setConstructionPlane(SbVec3f(0.0f, 0.0f, 0.0f),
	    SbVec3f(0.0f, 0.0f, 1.0f),
	    SbString("construction::workplane"));
    planeSnap.setQueryPoint(SbVec3f(0.4f, 0.6f, 2.0f));
    planeSnap.setTolerance(2.1f);
    planeSnap.apply(root);
    if (!planeSnap.hasCandidate() ||
	    planeSnap.getKind() != SoBRLSnapAction::CONSTRUCTION_PLANE ||
	    strcmp(planeSnap.getPath().getString(), "construction::workplane") != 0 ||
	    !nearly_equal(planeSnap.getPoint()[0], 0.4f) ||
	    !nearly_equal(planeSnap.getPoint()[1], 0.6f) ||
	    !nearly_equal(planeSnap.getPoint()[2], 0.0f))
	FAIL("snap action should support opt-in construction-plane policy");

    SoBRLSnapAction prioritySnap;
    prioritySnap.setPriorityPolicy(SoBRLSnapAction::FEATURE_PRIORITY);
    prioritySnap.setQueryPoint(SbVec3f(0.0f, -halfExtent, 0.0f));
    prioritySnap.setTolerance(0.1f);
    prioritySnap.apply(root);
    if (!prioritySnap.hasCandidate() ||
	    prioritySnap.getKind() != SoBRLSnapAction::MIDPOINT)
	FAIL("snap action should use feature priority to break equal-distance candidates");

    SoBRLMeasureAction measureAction;
    measureAction.setQueryPoint(SbVec3f(halfExtent, 0.2f, 0.0f));
    measureAction.apply(root);
    if (!measureAction.hasSegments())
	FAIL("measure action should find realized vlist segments");
    if (measureAction.getShapeCount() != 1 || measureAction.getSegmentCount() != 4)
	FAIL("measure action should count synthetic source geometry");
    if (!nearly_equal(measureAction.getTotalLength(), 14.0f))
	FAIL("measure action should total synthetic wireframe length");
    bbox = measureAction.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], -halfExtent) ||
	    !nearly_equal(bbox.getMax()[0], halfExtent) ||
	    !nearly_equal(bbox.getMin()[1], -halfExtent) ||
	    !nearly_equal(bbox.getMax()[1], halfExtent))
	FAIL("measure action should report synthetic source bounds");
    if (!measureAction.hasNearestSegment() ||
	    strcmp(measureAction.getNearestPath().getString(), "/prototype/arb8") != 0 ||
	    !nearly_equal(measureAction.getNearestPoint()[0], halfExtent) ||
	    !nearly_equal(measureAction.getNearestPoint()[1], 0.2f))
	FAIL("measure action should report nearest synthetic segment identity and point");

    SoBRLMeasureAction angleMeasure;
    angleMeasure.setQueryPoint(SbVec3f(halfExtent, halfExtent, 0.0f));
    angleMeasure.apply(root);
    if (!angleMeasure.hasAngle() ||
	    strcmp(angleMeasure.getAnglePath().getString(), "/prototype/arb8") != 0 ||
	    !nearly_equal(angleMeasure.getAngleDegrees(), 90.0f) ||
	    !nearly_equal(angleMeasure.getAnglePoint()[0], halfExtent) ||
	    !nearly_equal(angleMeasure.getAnglePoint()[1], halfExtent) ||
	    angleMeasure.getAnglePrimitiveIndexA() < 0 ||
	    angleMeasure.getAnglePrimitiveIndexB() < 0)
	FAIL("measure action should report connected-segment corner angles");

    source->sourceRevision = 4;
    if (!source->needsRealization() ||
	    source->realizationStatus.getValue() != SoBRLDatabaseSource::UNREALIZED ||
	    !source->stale.getValue())
	FAIL("database source field changes should mark realized geometry stale");
    SoBRLRealizeAction reRealize;
    reRealize.apply(root);
    if (reRealize.getRealizedSourceCount() != 1 ||
	    source->needsRealization() ||
	    source->realizedRevision.getValue() != 4)
	FAIL("realize action should refresh a source invalidated by field sensors");

    {
	SoPerspectiveCamera *camera = new SoPerspectiveCamera;
	BRLObolViewController viewController(root, camera);
	if (viewController.getSceneRoot() != root ||
		viewController.getRenderRoot() == NULL ||
		viewController.getViewport()->getSceneGraph() != root ||
		viewController.getRenderManager()->getSceneGraph() != viewController.getRenderRoot() ||
		viewController.getRenderManager()->getCamera() != camera ||
		viewController.getCamera() != camera ||
		viewController.getViewport()->getCamera() != camera ||
		viewController.getSceneController()->getSceneRoot() != root)
	    FAIL("view controller should expose the Obol scene root, viewport, and camera");
	viewController.clearRenderRequest();
	if (viewController.isRenderRequested())
	    FAIL("view controller should clear render requests");
	viewController.setViewportSize(320, 240);
	SbString renderReason;
	if (!viewController.isRenderRequested() ||
		strcmp(viewController.getRenderReason().getString(), "viewport-size") != 0 ||
		viewController.getViewportRegion().getWindowSize()[0] != 320 ||
		viewController.getViewportRegion().getWindowSize()[1] != 240 ||
		viewController.getViewport()->getViewportRegion().getWindowSize()[0] != 320 ||
		viewController.getViewport()->getViewportRegion().getWindowSize()[1] != 240 ||
		viewController.getRenderManager()->getViewportRegion().getWindowSize()[0] != 320 ||
		viewController.getRenderManager()->getViewportRegion().getWindowSize()[1] != 240)
	    FAIL("view controller should track viewport-driven render requests");
	if (!viewController.consumeRenderRequest(&renderReason) ||
		strcmp(renderReason.getString(), "viewport-size") != 0 ||
		viewController.isRenderRequested() ||
		viewController.consumeRenderRequest(NULL))
	    FAIL("view controller should consume render requests atomically");
	if (viewController.renderPending(FALSE, FALSE, NULL) ||
		viewController.isRenderRequested())
	    FAIL("view controller renderPending should be a no-op without queued work");
	viewController.clearRenderRequest();
	source->sourceRevision = 8;
	if (!viewController.realizePending() ||
		viewController.getLastVisitedSourceCount() != 1 ||
		viewController.getLastRealizedSourceCount() != 1 ||
		source->realizedSourceRevision.getValue() != 8 ||
		!viewController.isRenderRequested())
	    FAIL("view controller should delegate realization and request redraws");
    }

    root->unref();

    root = new SoSeparator;
    root->ref();

    source = new SoBRLDatabaseSource;
    source->path = "/sensor/source";
    source->sourceRevision = 1;
    root->addChild(source);

    SoBRLSceneController sensorController(root);
    if (!sensorController.realizePending())
	FAIL("scene controller should realize sensor test source");
    if (source->needsRealization() ||
	    source->realizedSourceRevision.getValue() != 1)
	FAIL("sensor test source should start realized");

    struct source_sensor_update_data sensorData;
    sensorData.source = source;
    sensorData.sourceRevision = 6;
    sensorData.fired = 0;

    SoOneShotSensor dynamicUpdate(source_sensor_update_cb, &sensorData);
    dynamicUpdate.schedule();
    SoDB::getSensorManager()->processDelayQueue(TRUE);
    SoDB::getSensorManager()->processDelayQueue(TRUE);

    if (sensorData.fired != 1 ||
	    source->sourceRevision.getValue() != 6)
	FAIL("Obol one-shot sensor should mutate BRL-CAD scene node fields");
    if (!source->needsRealization() ||
	    source->staleReason.getValue() != SoBRLDatabaseSource::STALE_SOURCE)
	FAIL("sensor-driven source field mutation should mark the source stale");
    if (!sensorController.realizePending() ||
	    sensorController.getLastRealizedSourceCount() != 1 ||
	    source->needsRealization() ||
	    source->realizedSourceRevision.getValue() != 6)
	FAIL("scene controller should refresh sensor-invalidated source state");

    SoGetBoundingBoxAction sensorBBoxAction(viewport);
    sensorBBoxAction.apply(root);
    bbox = sensorBBoxAction.getBoundingBox();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], -1.5f) ||
	    !nearly_equal(bbox.getMax()[0], 1.5f) ||
	    !nearly_equal(bbox.getMin()[1], -1.5f) ||
	    !nearly_equal(bbox.getMax()[1], 1.5f))
	FAIL("sensor-driven realization should publish updated Obol geometry");

    root->unref();

    root = new SoSeparator;
    root->ref();

    SoBRLDatabaseSource *leftSource = new SoBRLDatabaseSource;
    leftSource->path = "/controller/left";
    leftSource->sourceRevision = 1;
    root->addChild(leftSource);

    SoBRLDatabaseSource *rightSource = new SoBRLDatabaseSource;
    rightSource->path = "/controller/right";
    rightSource->sourceRevision = 2;
    root->addChild(rightSource);

    {
	SoBRLSceneController controller(root);
	if (!controller.realizePending())
	    FAIL("scene controller should realize a valid scene");
	if (controller.getLastVisitedSourceCount() != 2 ||
		controller.getLastRealizedSourceCount() != 2 ||
		controller.getLastFailedSourceCount() != 0)
	    FAIL("scene controller should report initial source realization diagnostics");
	if (leftSource->needsRealization() || rightSource->needsRealization())
	    FAIL("scene controller should leave realized sources current");
	if (leftSource->staleReason.getValue() != SoBRLDatabaseSource::STALE_NONE ||
		rightSource->staleReason.getValue() != SoBRLDatabaseSource::STALE_NONE)
	    FAIL("scene controller realization should clear stale reasons");

	if (!controller.realizePending())
	    FAIL("scene controller should allow no-op update passes");
	if (controller.getLastVisitedSourceCount() != 2 ||
		controller.getLastRealizedSourceCount() != 0)
	    FAIL("scene controller should not refresh clean sources");

	leftSource->inputsRevision = 5;
	if (!leftSource->needsRealization() ||
		rightSource->needsRealization() ||
		leftSource->staleReason.getValue() != SoBRLDatabaseSource::STALE_INPUTS)
	    FAIL("source input revision change should invalidate only that source");
	if (!controller.realizePending())
	    FAIL("scene controller should refresh field-invalidated source");
	if (controller.getLastVisitedSourceCount() != 2 ||
		controller.getLastRealizedSourceCount() != 1 ||
		controller.getLastFailedSourceCount() != 0)
	    FAIL("scene controller should report partial source refresh");
	if (leftSource->needsRealization() ||
		leftSource->realizedInputsRevision.getValue() != 5 ||
		rightSource->realizedInputsRevision.getValue() != 0)
	    FAIL("partial source refresh should update only the stale source revisions");
    }

    rightSource->viewRevision = 7;
    {
	SoBRLSceneController rebuiltController(root);
	if (!rebuiltController.realizePending())
	    FAIL("rebuilt scene controller should realize from the existing Obol scene");
	if (rebuiltController.getLastVisitedSourceCount() != 2 ||
		rebuiltController.getLastRealizedSourceCount() != 1)
	    FAIL("rebuilt scene controller should use scene-node stale state, not controller side tables");
	if (rightSource->needsRealization() ||
		rightSource->realizedViewRevision.getValue() != 7 ||
		rightSource->staleReason.getValue() != SoBRLDatabaseSource::STALE_NONE)
	    FAIL("rebuilt scene controller should clear view-revision invalidation");
    }

    root->unref();

    root = new SoSeparator;
    root->ref();

    SoBRLEditPreview *preview = new SoBRLEditPreview;
    preview->previewId = "preview::drag-box";
    preview->markSourceRevision(4);
    preview->markInputsRevision(2);
    if (!preview->needsRealization())
	FAIL("edit preview should become stale after revision changes");

    SbVec3f previewPoints[5] = {
	SbVec3f(-1.0f, -1.0f, 0.0f),
	SbVec3f( 1.0f, -1.0f, 0.0f),
	SbVec3f( 1.0f,  1.0f, 0.0f),
	SbVec3f(-1.0f,  1.0f, 0.0f),
	SbVec3f(-1.0f, -1.0f, 0.0f)
    };
    int32_t previewCommands[5] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW
    };
    SbMatrix previewMatrix;
    previewMatrix.setTranslate(SbVec3f(40.0f, 0.0f, 0.0f));
    SoBRLVListShape *previewShape = preview->setTransformedLineSet(
	    "preview::drag-box/outline", previewMatrix, previewPoints, previewCommands, 5);
    if (!previewShape)
	FAIL("edit preview should accept transformed line geometry");
    if (preview->needsRealization())
	FAIL("edit preview should be current after publishing geometry");
    if (preview->previewStatus.getValue() != SoBRLEditPreview::CURRENT)
	FAIL("edit preview status should be CURRENT after publishing geometry");
    if (preview->realizedSourceRevision.getValue() != 4 ||
	    preview->realizedInputsRevision.getValue() != 2)
	FAIL("edit preview should synchronize realized revision fields");
    if (!shape_extents_match(previewShape, -1.0f, 1.0f, -1.0f, 1.0f, 0.0f, 0.0f))
	FAIL("edit preview shape should keep local coordinates under Obol transform nodes");
    root->addChild(preview);

    SoGetBoundingBoxAction previewBBoxAction(viewport);
    previewBBoxAction.apply(root);
    bbox = previewBBoxAction.getBoundingBox();
    if (bbox.isEmpty())
	FAIL("edit preview should contribute an Obol bounding box");
    if (!nearly_equal(bbox.getMin()[0], 39.0f) ||
	    !nearly_equal(bbox.getMax()[0], 41.0f) ||
	    !nearly_equal(bbox.getMin()[1], -1.0f) ||
	    !nearly_equal(bbox.getMax()[1], 1.0f))
	FAIL("edit preview bounding box should reflect its Obol transform");

    SoRayPickAction previewPick(viewport);
    previewPick.setRay(SbVec3f(40.0f, -1.0f, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    previewPick.apply(root);
    pickedPoint = previewPick.getPickedPoint();
    if (!pickedPoint)
	FAIL("edit preview ray pick should hit transformed preview geometry");
    rawDetail = pickedPoint->getDetail(previewShape);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("edit preview pick should return a BRL-CAD Obol pick detail");
    pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(), "preview::drag-box/outline") != 0)
	FAIL("edit preview pick should preserve preview identity");

    SoBRLSnapAction previewSnap;
    previewSnap.setQueryPoint(SbVec3f(41.0f, 0.2f, 0.0f));
    previewSnap.setTolerance(0.3f);
    previewSnap.apply(root);
    if (!previewSnap.hasCandidate())
	FAIL("edit preview snap should find transformed preview geometry");
    if (strcmp(previewSnap.getPath().getString(), "preview::drag-box/outline") != 0)
	FAIL("edit preview snap should preserve preview identity");
    if (!nearly_equal(previewSnap.getPoint()[0], 41.0f) ||
	    !nearly_equal(previewSnap.getPoint()[1], 0.2f))
	FAIL("edit preview snap point should be in transformed model coordinates");

    SoBRLMeasureAction previewMeasure;
    previewMeasure.setQueryPoint(SbVec3f(41.0f, 0.2f, 0.0f));
    previewMeasure.apply(root);
    if (!previewMeasure.hasSegments() ||
	    previewMeasure.getShapeCount() != 1 ||
	    previewMeasure.getSegmentCount() != 4)
	FAIL("edit preview measure should find transformed preview segments");
    if (!nearly_equal(previewMeasure.getTotalLength(), 8.0f))
	FAIL("edit preview measure should total transformed preview wire length");
    bbox = previewMeasure.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 39.0f) ||
	    !nearly_equal(bbox.getMax()[0], 41.0f) ||
	    !nearly_equal(bbox.getMin()[1], -1.0f) ||
	    !nearly_equal(bbox.getMax()[1], 1.0f))
	FAIL("edit preview measure should report transformed bounds");
    if (!previewMeasure.hasNearestSegment() ||
	    strcmp(previewMeasure.getNearestPath().getString(), "preview::drag-box/outline") != 0 ||
	    !nearly_equal(previewMeasure.getNearestPoint()[0], 41.0f) ||
	    !nearly_equal(previewMeasure.getNearestPoint()[1], 0.2f))
	FAIL("edit preview measure should report transformed preview identity");

    SoBRLExportAction previewExport;
    previewExport.apply(root);
    if (previewExport.getLineCount() != 4)
	FAIL("edit preview export should collect transformed line records");
    if (export_path_count(previewExport, "preview::drag-box/outline") != 4)
	FAIL("edit preview export should preserve preview line identity");
    if (!nearly_equal(previewExport.getLine(0).a[0], 39.0f) ||
	    !nearly_equal(previewExport.getLine(0).b[0], 41.0f))
	FAIL("edit preview export should apply Obol transform state");

    SbVec3f replacementPoints[5] = {
	SbVec3f(-2.0f, -0.5f, 0.0f),
	SbVec3f( 2.0f, -0.5f, 0.0f),
	SbVec3f( 2.0f,  0.5f, 0.0f),
	SbVec3f(-2.0f,  0.5f, 0.0f),
	SbVec3f(-2.0f, -0.5f, 0.0f)
    };
    int32_t replacementCommands[5] = {
	SoBRLVListShape::MOVE,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW,
	SoBRLVListShape::DRAW
    };
    SbMatrix replacementMatrix;
    replacementMatrix.setTranslate(SbVec3f(45.0f, 0.0f, 0.0f));
    SoBRLVListShape *replacementShape = preview->setTransformedLineSet(
	    "preview::drag-box/replacement-outline", replacementMatrix,
	    replacementPoints, replacementCommands, 5);
    if (!replacementShape)
	FAIL("edit preview should accept replacement geometry");
    if (preview->getNumChildren() != 1)
	FAIL("edit preview replacement should leave exactly one current child subtree");
    if (preview->needsRealization() ||
	    preview->previewStatus.getValue() != SoBRLEditPreview::CURRENT)
	FAIL("edit preview replacement should publish current scene content");
    if (!shape_extents_match(replacementShape, -2.0f, 2.0f, -0.5f, 0.5f, 0.0f, 0.0f))
	FAIL("replacement preview shape should keep its own local coordinates");

    SoBRLExportAction replacementExport;
    replacementExport.apply(root);
    if (replacementExport.getLineCount() != 4)
	FAIL("edit preview replacement should not leave stale exported line records");
    if (export_path_count(replacementExport, "preview::drag-box/outline") != 0)
	FAIL("edit preview replacement export should not expose old preview identity");
    if (export_path_count(replacementExport, "preview::drag-box/replacement-outline") != 4)
	FAIL("edit preview replacement export should expose only the current identity");
    bbox = replacementExport.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 43.0f) ||
	    !nearly_equal(bbox.getMax()[0], 47.0f) ||
	    !nearly_equal(bbox.getMin()[1], -0.5f) ||
	    !nearly_equal(bbox.getMax()[1], 0.5f))
	FAIL("edit preview replacement export should use current transformed bounds");

    SoRayPickAction stalePreviewPick(viewport);
    stalePreviewPick.setRay(SbVec3f(40.0f, -1.0f, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    stalePreviewPick.apply(root);
    if (stalePreviewPick.getPickedPoint())
	FAIL("edit preview replacement should remove old pickable geometry from the scene");

    SoRayPickAction replacementPick(viewport);
    replacementPick.setRay(SbVec3f(45.0f, -0.5f, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    replacementPick.apply(root);
    pickedPoint = replacementPick.getPickedPoint();
    if (!pickedPoint)
	FAIL("edit preview replacement ray pick should hit current geometry");
    rawDetail = pickedPoint->getDetail(replacementShape);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("edit preview replacement pick should return a current BRL-CAD detail");
    pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(), "preview::drag-box/replacement-outline") != 0)
	FAIL("edit preview replacement pick should preserve current preview identity");

    SoBRLSnapAction replacementSnap;
    replacementSnap.setQueryPoint(SbVec3f(47.0f, 0.25f, 0.0f));
    replacementSnap.setTolerance(0.3f);
    replacementSnap.apply(root);
    if (!replacementSnap.hasCandidate() ||
	    strcmp(replacementSnap.getPath().getString(), "preview::drag-box/replacement-outline") != 0)
	FAIL("edit preview replacement snap should use only current preview identity");
    if (!nearly_equal(replacementSnap.getPoint()[0], 47.0f) ||
	    !nearly_equal(replacementSnap.getPoint()[1], 0.25f))
	FAIL("edit preview replacement snap should use current transformed geometry");

    SoBRLMeasureAction replacementMeasure;
    replacementMeasure.setQueryPoint(SbVec3f(47.0f, 0.25f, 0.0f));
    replacementMeasure.apply(root);
    if (!replacementMeasure.hasSegments() ||
	    replacementMeasure.getShapeCount() != 1 ||
	    replacementMeasure.getSegmentCount() != 4 ||
	    !nearly_equal(replacementMeasure.getTotalLength(), 10.0f))
	FAIL("edit preview replacement measure should not include stale segments");
    if (!replacementMeasure.hasNearestSegment() ||
	    strcmp(replacementMeasure.getNearestPath().getString(), "preview::drag-box/replacement-outline") != 0)
	FAIL("edit preview replacement measure should report current preview identity");

    preview->inputsRevision = 3;
    if (!preview->needsRealization() ||
	    preview->previewStatus.getValue() != SoBRLEditPreview::STALE)
	FAIL("edit preview field changes should mark published geometry stale");
    preview->setLineSet("preview::drag-box/final-outline", previewPoints, previewCommands, 5);
    if (preview->needsRealization() ||
	    preview->realizedInputsRevision.getValue() != 3)
	FAIL("edit preview should become current after republishing field-invalidated geometry");
    if (preview->getNumChildren() != 1)
	FAIL("edit preview field republish should replace previous preview content");
    SoBRLExportAction finalPreviewExport;
    finalPreviewExport.apply(root);
    if (finalPreviewExport.getLineCount() != 4 ||
	    export_path_count(finalPreviewExport, "preview::drag-box/replacement-outline") != 0 ||
	    export_path_count(finalPreviewExport, "preview::drag-box/final-outline") != 4)
	FAIL("edit preview field republish should expose only final preview content");
    preview->previewId = "preview::drag-box-renamed";
    if (!preview->needsRealization() ||
	    preview->previewStatus.getValue() != SoBRLEditPreview::STALE)
	FAIL("edit preview identity field changes should invalidate published geometry");
    preview->clearPreview();
    if (preview->needsRealization() ||
	    preview->previewStatus.getValue() != SoBRLEditPreview::EMPTY ||
	    preview->getNumChildren() != 0)
	FAIL("edit preview clear should remove transient scene content and become current-empty");

    root->unref();

    root = new SoSeparator;
    root->ref();

    SoMatrixTransform *meshTransform = new SoMatrixTransform;
    SbMatrix meshMatrix;
    meshMatrix.setTranslate(SbVec3f(60.0f, 0.0f, 0.0f));
    meshTransform->matrix = meshMatrix;
    root->addChild(meshTransform);

    SoBRLMeshShape *mesh = new SoBRLMeshShape;
    mesh->sourcePath = "mesh::tri-test";
    mesh->sourceId = 77;
    SbVec3f meshPoints[4] = {
	SbVec3f(0.0f, 0.0f, 0.0f),
	SbVec3f(2.0f, 0.0f, 0.0f),
	SbVec3f(0.0f, 2.0f, 0.0f),
	SbVec3f(2.0f, 2.0f, 0.0f)
    };
    int32_t meshIndices[6] = {0, 1, 2, 1, 3, 2};
    mesh->setIndexedTriangles(meshPoints, 4, meshIndices, 6);
    if (mesh->getTriangleCount() != 2)
	FAIL("mesh shape should expose indexed triangle count");
    root->addChild(mesh);

    SoGetBoundingBoxAction meshBBoxAction(viewport);
    meshBBoxAction.apply(root);
    bbox = meshBBoxAction.getBoundingBox();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 60.0f) ||
	    !nearly_equal(bbox.getMax()[0], 62.0f) ||
	    !nearly_equal(bbox.getMin()[1], 0.0f) ||
	    !nearly_equal(bbox.getMax()[1], 2.0f))
	FAIL("mesh shape should contribute transformed Obol bounds");

    SoRayPickAction meshPick(viewport);
    meshPick.setRay(SbVec3f(60.25f, 0.25f, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    meshPick.apply(root);
    pickedPoint = meshPick.getPickedPoint();
    if (!pickedPoint)
	FAIL("mesh ray pick should hit triangle geometry");
    rawDetail = pickedPoint->getDetail(mesh);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("mesh pick should return a BRL-CAD Obol pick detail");
    pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(), "mesh::tri-test") != 0 ||
	    pickDetail->getPrimitiveKind() != SoBRLPickDetail::FACE ||
	    pickDetail->getSourceId() != 77)
	FAIL("mesh pick should preserve face identity");
    if (pickDetail->getPrimitiveIndex() != 0 ||
	    pickDetail->getFaceVertexIndexA() != 0 ||
	    pickDetail->getFaceVertexIndexB() != 1 ||
	    pickDetail->getFaceVertexIndexC() != 2 ||
	    pickDetail->getFaceVertexIndex(3) != -1)
	FAIL("mesh pick should preserve picked triangle vertex identity");
    if (pickDetail->getNearestFaceEdgeSlot() != 0 ||
	    pickDetail->getNearestFaceEdgeVertexIndexA() != 0 ||
	    pickDetail->getNearestFaceEdgeVertexIndexB() != 1 ||
	    pickDetail->getNearestFaceVertexSlot() != 0 ||
	    pickDetail->getNearestFaceVertexIndex() != 0)
	FAIL("mesh pick should preserve nearest edge and vertex identity");
    if (!nearly_equal(pickDetail->getModelPoint()[0], 0.25f) ||
	    !nearly_equal(pickDetail->getModelPoint()[1], 0.25f) ||
	    !nearly_equal(pickDetail->getModelPoint()[2], 0.0f))
	FAIL("mesh pick should preserve object-space hit point");

    SoBRLExportAction meshExport;
    meshExport.apply(root);
    if (meshExport.getLineCount() != 0 ||
	    meshExport.getTriangleCount() != 2)
	FAIL("mesh export should collect triangle records without wire records");
    const SoBRLExportAction::TriangleRecord &tri = meshExport.getTriangle(0);
    if (strcmp(tri.path.getString(), "mesh::tri-test") != 0 ||
	    tri.sourceId != 77 ||
	    tri.primitiveIndex != 0 ||
	    tri.vertexIndexA != 0 ||
	    tri.vertexIndexB != 1 ||
	    tri.vertexIndexC != 2)
	FAIL("mesh export should preserve triangle identity");
    if (tri.lodAvailable ||
	    tri.lodActiveLevel != -1 ||
	    tri.lodFaceCount != 0 ||
	    tri.lodPointCount != 0 ||
	    tri.lodOriginalPointCount != 0 ||
	    tri.lodNormalCount != 0 ||
	    tri.lodHasSnappedPoints ||
	    tri.lodHasNormals)
	FAIL("mesh export should carry unavailable LoD metadata defaults");

    BRLObolLodRequest meshLodRequest;
    meshLodRequest.databaseId = "db://prototype";
    meshLodRequest.sourceRevision = 77;
    meshLodRequest.objectPath = "mesh::tri-test";
    meshLodRequest.objectName = "tri-test";
    meshLodRequest.viewRevision = 1;
    meshLodRequest.policyRevision = 1;
    meshLodRequest.drawMode = BRLOBOL_LOD_DRAW_SHADED;
    meshLodRequest.providerId = "prototype-provider";
    meshLodRequest.providerVersion = "1";
    meshLodRequest.qualityTier = BRLOBOL_LOD_QUALITY_PROXY;
    meshLodRequest.bounds = SbBox3f(SbVec3f(0.0f, 0.0f, 0.0f),
	    SbVec3f(2.0f, 2.0f, 1.0f));

    std::vector<BRLObolLodDependency> meshDependencies;
    BRLObolLodDependency meshDependency;
    meshDependency.objectPath = "mesh::tri-test/child";
    meshDependency.objectName = "child";
    meshDependency.sourceRevision = 123;
    meshDependency.sourceContentHash = 456;
    meshDependency.requiredQualityTier = BRLOBOL_LOD_QUALITY_PROXY;
    meshDependency.optional = TRUE;
    meshDependencies.push_back(meshDependency);
    BRLObolLodResult stagedResult =
	brlobol_lod_directory_result(meshLodRequest, meshDependencies);
    if (!mesh->applyStagedLodResult(stagedResult, &meshLodRequest) ||
	    !mesh->lodStagedAvailable.getValue() ||
	    mesh->lodResultKind.getValue() != BRLOBOL_LOD_RESULT_DIRECTORY ||
	    mesh->lodDependencyPath.getNum() != 1 ||
	    strcmp(mesh->lodDependencyPath[0].getString(),
		"mesh::tri-test/child") != 0 ||
	    strcmp(mesh->lodDependencySourceRevision[0].getString(),
		"123") != 0 ||
	    !mesh->lodDependencyOptional[0])
	FAIL("mesh shape should consume staged LoD dependency results");

    std::vector<BRLObolLodAttribute> meshAttributes;
    BRLObolLodAttribute meshAttribute;
    meshAttribute.name = "display.color";
    meshAttribute.value = "255 0 0";
    meshAttributes.push_back(meshAttribute);
    stagedResult = brlobol_lod_attributes_result(meshLodRequest,
	    meshAttributes);
    if (!mesh->applyStagedLodResult(stagedResult, &meshLodRequest) ||
	    mesh->lodDependencyPath.getNum() != 0 ||
	    mesh->lodAttributeName.getNum() != 1 ||
	    strcmp(mesh->lodAttributeName[0].getString(),
		"display.color") != 0 ||
	    strcmp(mesh->lodAttributeValue[0].getString(), "255 0 0") != 0)
	FAIL("mesh shape should consume staged LoD attribute results");

    BRLObolLodProxy meshProxy;
    meshProxy.kind = BRLOBOL_LOD_PROXY_OBB;
    meshProxy.bounds = meshLodRequest.bounds;
    meshProxy.center.setValue(1.0f, 1.0f, 0.5f);
    meshProxy.axisX.setValue(1.0f, 0.0f, 0.0f);
    meshProxy.axisY.setValue(0.0f, 1.0f, 0.0f);
    meshProxy.axisZ.setValue(0.0f, 0.0f, 1.0f);
    meshProxy.halfExtents.setValue(1.0f, 1.0f, 0.5f);
    stagedResult = brlobol_lod_proxy_result(meshLodRequest, meshProxy, NULL);
    if (!mesh->applyStagedLodResult(stagedResult, &meshLodRequest) ||
	    mesh->lodProxyKind.getValue() != BRLOBOL_LOD_PROXY_OBB ||
	    !nearly_equal(mesh->lodProxyCenter.getValue()[0], 1.0f) ||
	    !nearly_equal(mesh->lodProxyHalfExtents.getValue()[1], 1.0f))
	FAIL("mesh shape should consume staged LoD proxy results");

    BRLObolLodResult staleResult = stagedResult;
    staleResult.request.viewRevision++;
    staleResult.cacheKey = brlobol_lod_cache_key(staleResult.request);
    if (mesh->applyStagedLodResult(staleResult, &meshLodRequest) ||
	    mesh->lodStagedAvailable.getValue() ||
	    mesh->lodProviderStatus.getValue() != BRLOBOL_LOD_PROVIDER_STALE ||
	    mesh->lodAttributeName.getNum() != 0)
	FAIL("mesh shape should reject stale staged LoD results");

    mesh->clearStagedLodResult();
    if (mesh->lodStagedAvailable.getValue() ||
	    mesh->lodProxyKind.getValue() != BRLOBOL_LOD_PROXY_NONE ||
	    mesh->lodAttributeName.getNum() != 0)
	FAIL("mesh shape should clear staged LoD result fields");

    if (!nearly_equal(tri.a[0], 60.0f) ||
	    !nearly_equal(tri.b[0], 62.0f) ||
	    !nearly_equal(tri.c[1], 2.0f))
	FAIL("mesh export should apply Obol transform state to triangle vertices");

    SoBRLMeasureAction meshMeasure;
    meshMeasure.setQueryPoint(SbVec3f(61.0f, 1.0f, 1.0f));
    meshMeasure.apply(root);
    if (!meshMeasure.hasFaces() ||
	    meshMeasure.getShapeCount() != 1 ||
	    meshMeasure.getTriangleCount() != 2 ||
	    !nearly_equal(meshMeasure.getSurfaceArea(), 4.0f))
	FAIL("mesh measure should report transformed mesh face metrics");
    bbox = meshMeasure.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 60.0f) ||
	    !nearly_equal(bbox.getMax()[0], 62.0f) ||
	    !nearly_equal(bbox.getMin()[1], 0.0f) ||
	    !nearly_equal(bbox.getMax()[1], 2.0f))
	FAIL("mesh measure should report transformed mesh bounds");
    if (!meshMeasure.hasNearestPrimitive() ||
	    meshMeasure.getNearestPrimitiveKind() != SoBRLMeasureAction::FACE ||
	    strcmp(meshMeasure.getNearestPath().getString(), "mesh::tri-test") != 0 ||
	    !nearly_equal(meshMeasure.getNearestPoint()[0], 61.0f) ||
	    !nearly_equal(meshMeasure.getNearestPoint()[1], 1.0f) ||
	    !nearly_equal(meshMeasure.getNearestPoint()[2], 0.0f))
	FAIL("mesh measure should report nearest transformed mesh face identity");

    SoBRLSnapAction meshFaceSnap;
    meshFaceSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    meshFaceSnap.setQueryPoint(SbVec3f(61.0f, 1.0f, 1.0f));
    meshFaceSnap.setTolerance(1.1f);
    meshFaceSnap.apply(root);
    if (!meshFaceSnap.hasCandidate() ||
	    meshFaceSnap.getKind() != SoBRLSnapAction::FACE_NEAREST ||
	    strcmp(meshFaceSnap.getPath().getString(), "mesh::tri-test") != 0 ||
	    !nearly_equal(meshFaceSnap.getPoint()[0], 61.0f) ||
	    !nearly_equal(meshFaceSnap.getPoint()[1], 1.0f) ||
	    !nearly_equal(meshFaceSnap.getPoint()[2], 0.0f))
	FAIL("snap action should support mesh face-nearest policy");

    SoBRLMeasureAction localMeshMeasure;
    localMeshMeasure.setCoordinateSpace(SoBRLMeasureAction::PATH_LOCAL_SPACE);
    localMeshMeasure.setQueryPoint(SbVec3f(1.0f, 1.0f, 1.0f));
    localMeshMeasure.apply(root);
    bbox = localMeshMeasure.getBounds();
    if (!localMeshMeasure.hasFaces() ||
	    bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 0.0f) ||
	    !nearly_equal(bbox.getMax()[0], 2.0f) ||
	    !localMeshMeasure.hasNearestPrimitive() ||
	    !nearly_equal(localMeshMeasure.getNearestPoint()[0], 1.0f) ||
	    !nearly_equal(localMeshMeasure.getNearestPoint()[1], 1.0f) ||
	    !nearly_equal(localMeshMeasure.getNearestPoint()[2], 0.0f))
	FAIL("measure action should support path-local coordinate policy on transformed geometry");

    SoBRLSnapAction localMeshSnap;
    localMeshSnap.setCoordinateSpace(SoBRLSnapAction::PATH_LOCAL_SPACE);
    localMeshSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    localMeshSnap.setQueryPoint(SbVec3f(1.0f, 1.0f, 1.0f));
    localMeshSnap.setTolerance(1.1f);
    localMeshSnap.apply(root);
    if (!localMeshSnap.hasCandidate() ||
	    localMeshSnap.getKind() != SoBRLSnapAction::FACE_NEAREST ||
	    strcmp(localMeshSnap.getPath().getString(), "mesh::tri-test") != 0 ||
	    !nearly_equal(localMeshSnap.getPoint()[0], 1.0f) ||
	    !nearly_equal(localMeshSnap.getPoint()[1], 1.0f) ||
	    !nearly_equal(localMeshSnap.getPoint()[2], 0.0f))
	FAIL("snap action should support path-local coordinate policy on transformed geometry");

    SoBRLSnapAction selectedOnlyMiss;
    selectedOnlyMiss.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    selectedOnlyMiss.setSelectionFilter(SoBRLSnapAction::SELECTED_ONLY);
    selectedOnlyMiss.setQueryPoint(SbVec3f(61.0f, 1.0f, 1.0f));
    selectedOnlyMiss.setTolerance(1.1f);
    selectedOnlyMiss.apply(root);
    if (selectedOnlyMiss.hasCandidate())
	FAIL("selected-only snap policy should ignore unselected mesh geometry");

    mesh->selectedPrimitive.set1Value(0, 1);
    mesh->highlightedPrimitive.set1Value(0, 1);
    SoBRLSnapAction selectedOnlySnap;
    selectedOnlySnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    selectedOnlySnap.setSelectionFilter(SoBRLSnapAction::SELECTED_ONLY);
    selectedOnlySnap.setQueryPoint(SbVec3f(61.5f, 1.5f, 1.0f));
    selectedOnlySnap.setTolerance(1.1f);
    selectedOnlySnap.apply(root);
    if (!selectedOnlySnap.hasCandidate() ||
	    selectedOnlySnap.getKind() != SoBRLSnapAction::FACE_NEAREST ||
	    selectedOnlySnap.getPrimitiveIndex() != 1 ||
	    strcmp(selectedOnlySnap.getPath().getString(), "mesh::tri-test") != 0)
	FAIL("selected-only snap policy should accept selected mesh face geometry");

    SoBRLMeasureAction selectedOnlyMeshMeasure;
    selectedOnlyMeshMeasure.setSelectionFilter(SoBRLMeasureAction::SELECTED_ONLY);
    selectedOnlyMeshMeasure.setQueryPoint(SbVec3f(61.5f, 1.5f, 1.0f));
    selectedOnlyMeshMeasure.apply(root);
    if (!selectedOnlyMeshMeasure.hasFaces() ||
	    selectedOnlyMeshMeasure.getTriangleCount() != 1 ||
	    selectedOnlyMeshMeasure.getNearestPrimitiveIndex() != 1 ||
	    !nearly_equal(selectedOnlyMeshMeasure.getSurfaceArea(), 2.0f))
	FAIL("selected-only measure policy should accept selected mesh face geometry");

    SoBRLMeasureAction highlightedOnlyMeshMeasure;
    highlightedOnlyMeshMeasure.setHighlightFilter(SoBRLMeasureAction::HIGHLIGHTED_ONLY);
    highlightedOnlyMeshMeasure.setQueryPoint(SbVec3f(61.5f, 1.5f, 1.0f));
    highlightedOnlyMeshMeasure.apply(root);
    if (!highlightedOnlyMeshMeasure.hasFaces() ||
	    highlightedOnlyMeshMeasure.getTriangleCount() != 1 ||
	    highlightedOnlyMeshMeasure.getNearestPrimitiveIndex() != 1 ||
	    !nearly_equal(highlightedOnlyMeshMeasure.getSurfaceArea(), 2.0f))
	FAIL("highlighted-only measure policy should accept highlighted mesh face geometry");

    SoBRLExportAction meshSubEntityExport;
    meshSubEntityExport.apply(root);
    if (meshSubEntityExport.getTriangleCount() != 2 ||
	    meshSubEntityExport.getTriangle(0).selected ||
	    meshSubEntityExport.getTriangle(0).highlighted ||
	    !meshSubEntityExport.getTriangle(1).selected ||
	    !meshSubEntityExport.getTriangle(1).highlighted)
	FAIL("mesh export should carry per-face selected and highlighted draw intent");

    mesh->selected = TRUE;
    mesh->highlighted = TRUE;
    mesh->ghosted = TRUE;
    mesh->hiddenLine = TRUE;
    mesh->editEmphasis = TRUE;
    mesh->lodPolicy = 3;
    mesh->colorOverride = TRUE;
    mesh->color = SbColor(0.2f, 0.4f, 0.9f);

    SoBRLExportAction meshIntentExport;
    meshIntentExport.apply(root);
    if (meshIntentExport.getTriangleCount() != 2 ||
	    !meshIntentExport.getTriangle(0).selected ||
	    !meshIntentExport.getTriangle(0).highlighted ||
	    !meshIntentExport.getTriangle(0).ghosted ||
	    !meshIntentExport.getTriangle(0).hiddenLine ||
	    !meshIntentExport.getTriangle(0).editEmphasis ||
	    meshIntentExport.getTriangle(0).lodPolicy != 3 ||
	    !meshIntentExport.getTriangle(0).colorOverride ||
	    !nearly_equal(meshIntentExport.getTriangle(0).color[0], 0.2f) ||
	    !nearly_equal(meshIntentExport.getTriangle(0).color[1], 0.4f) ||
	    !nearly_equal(meshIntentExport.getTriangle(0).color[2], 0.9f))
	FAIL("mesh export should carry complete draw intent");

    mesh->visible = FALSE;
    SoGetBoundingBoxAction hiddenMeshBBoxAction(viewport);
    hiddenMeshBBoxAction.apply(root);
    if (!hiddenMeshBBoxAction.getBoundingBox().isEmpty())
	FAIL("hidden mesh should not contribute an Obol bounding box");
    SoBRLExportAction hiddenMeshExport;
    hiddenMeshExport.apply(root);
    if (hiddenMeshExport.getTriangleCount() != 0)
	FAIL("hidden mesh should not be exported");
    SoBRLMeasureAction hiddenMeshMeasure;
    hiddenMeshMeasure.apply(root);
    if (hiddenMeshMeasure.hasFaces())
	FAIL("hidden mesh should not contribute face measurements");
    SoRayPickAction hiddenMeshPick(viewport);
    hiddenMeshPick.setRay(SbVec3f(60.25f, 0.25f, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    hiddenMeshPick.apply(root);
    if (hiddenMeshPick.getPickedPoint())
	FAIL("hidden mesh should not be pickable");
    SoBRLSnapAction hiddenMeshSnap;
    hiddenMeshSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    hiddenMeshSnap.setQueryPoint(SbVec3f(61.0f, 1.0f, 1.0f));
    hiddenMeshSnap.setTolerance(1.1f);
    hiddenMeshSnap.apply(root);
    if (hiddenMeshSnap.hasCandidate())
	FAIL("hidden mesh should not contribute snap candidates");

    root->unref();

    root = new SoSeparator;
    root->ref();

    SoBRLGrid *grid = new SoBRLGrid;
    grid->overlayId = "overlay::work-grid";
    grid->center = SbVec3f(0.0f, 0.0f, 0.0f);
    grid->spacing = 1.0f;
    grid->divisions = 2;
    SoBRLVListShape *gridShape = grid->rebuildGeometry();
    if (!gridShape)
	FAIL("grid overlay should build Obol line geometry");
    if (gridShape->getSegmentCount() != 10)
	FAIL("grid overlay should build expected grid line count");
    root->addChild(grid);

    SoBRLAxes *axes = new SoBRLAxes;
    axes->overlayId = "overlay::model-axes";
    axes->origin = SbVec3f(10.0f, 0.0f, 0.0f);
    axes->size = 3.0f;
    SoBRLVListShape *axesShape = axes->rebuildGeometry();
    if (!axesShape)
	FAIL("axes overlay should build Obol line geometry");
    if (axesShape->getSegmentCount() != 3)
	FAIL("axes overlay should build three axis line segments");
    root->addChild(axes);

    SoBRLADC *adc = new SoBRLADC;
    adc->overlayId = "overlay::adc";
    adc->center = SbVec3f(20.0f, 0.0f, 0.0f);
    adc->angleDegrees = 0.0f;
    adc->distance = 4.0f;
    adc->crosshairSize = 1.0f;
    adc->tickSize = 0.5f;
    SoBRLVListShape *adcShape = adc->rebuildGeometry();
    if (!adcShape)
	FAIL("ADC overlay should build Obol line geometry");
    if (adcShape->getSegmentCount() != 4)
	FAIL("ADC overlay should build crosshair, angle, and distance segments");
    root->addChild(adc);

    SoGetBoundingBoxAction overlayBBoxAction(viewport);
    overlayBBoxAction.apply(root);
    bbox = overlayBBoxAction.getBoundingBox();
    if (bbox.isEmpty())
	FAIL("overlay nodes should contribute Obol bounding boxes");
    if (!nearly_equal(bbox.getMin()[0], -2.0f) ||
	    !nearly_equal(bbox.getMax()[0], 24.0f) ||
	    !nearly_equal(bbox.getMin()[1], -2.0f) ||
	    !nearly_equal(bbox.getMax()[1], 3.0f) ||
	    !nearly_equal(bbox.getMin()[2], 0.0f) ||
	    !nearly_equal(bbox.getMax()[2], 3.0f))
	FAIL("overlay bounding box should reflect grid, axes, and ADC field geometry");

    SoRayPickAction overlayPick(viewport);
    overlayPick.setRay(SbVec3f(11.0f, 0.0f, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    overlayPick.apply(root);
    pickedPoint = overlayPick.getPickedPoint();
    if (!pickedPoint)
	FAIL("overlay ray pick should hit axes geometry");
    rawDetail = pickedPoint->getDetail(axesShape);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("overlay pick should return a BRL-CAD Obol pick detail");
    pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(), "overlay::model-axes") != 0)
	FAIL("overlay pick should preserve axes identity");

    SoRayPickAction adcPick(viewport);
    adcPick.setRay(SbVec3f(22.0f, 0.0f, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    adcPick.apply(root);
    pickedPoint = adcPick.getPickedPoint();
    if (!pickedPoint)
	FAIL("ADC ray pick should hit angle-distance geometry");
    rawDetail = pickedPoint->getDetail(adcShape);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("ADC pick should return a BRL-CAD Obol pick detail");
    pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(), "overlay::adc") != 0)
	FAIL("ADC pick should preserve overlay identity");

    SoBRLSnapAction overlaySnap;
    overlaySnap.setQueryPoint(SbVec3f(2.0f, 0.3f, 0.0f));
    overlaySnap.setTolerance(0.4f);
    overlaySnap.apply(root);
    if (!overlaySnap.hasCandidate())
	FAIL("overlay snap should find grid geometry");
    if (strcmp(overlaySnap.getPath().getString(), "overlay::work-grid") != 0)
	FAIL("overlay snap should preserve grid identity");
    if (!nearly_equal(overlaySnap.getPoint()[0], 2.0f) ||
	    !nearly_equal(overlaySnap.getPoint()[1], 0.3f))
	FAIL("overlay snap point should lie on the grid line");

    SoBRLSnapAction adcSnap;
    adcSnap.setQueryPoint(SbVec3f(22.0f, 0.3f, 0.0f));
    adcSnap.setTolerance(0.4f);
    adcSnap.apply(root);
    if (!adcSnap.hasCandidate())
	FAIL("ADC snap should find angle-distance geometry");
    if (strcmp(adcSnap.getPath().getString(), "overlay::adc") != 0)
	FAIL("ADC snap should preserve overlay identity");
    if (adcSnap.getKind() != SoBRLSnapAction::LINE_NEAREST ||
	    !nearly_equal(adcSnap.getPoint()[0], 22.0f) ||
	    !nearly_equal(adcSnap.getPoint()[1], 0.0f))
	FAIL("ADC snap should return nearest point on the angle-distance line");

    SoBRLMeasureAction overlayMeasure;
    overlayMeasure.setQueryPoint(SbVec3f(11.2f, 0.0f, 0.0f));
    overlayMeasure.apply(root);
    if (!overlayMeasure.hasSegments() ||
	    overlayMeasure.getShapeCount() != 3 ||
	    overlayMeasure.getSegmentCount() != 17)
	FAIL("overlay measure should count grid, axes, and ADC segments");
    if (!nearly_equal(overlayMeasure.getTotalLength(), 58.0f))
	FAIL("overlay measure should total grid, axes, and ADC line length");
    bbox = overlayMeasure.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], -2.0f) ||
	    !nearly_equal(bbox.getMax()[0], 24.0f) ||
	    !nearly_equal(bbox.getMin()[1], -2.0f) ||
	    !nearly_equal(bbox.getMax()[1], 3.0f) ||
	    !nearly_equal(bbox.getMin()[2], 0.0f) ||
	    !nearly_equal(bbox.getMax()[2], 3.0f))
	FAIL("overlay measure should report grid, axes, and ADC bounds");
    if (!overlayMeasure.hasNearestSegment() ||
	    strcmp(overlayMeasure.getNearestPath().getString(), "overlay::model-axes") != 0 ||
	    !nearly_equal(overlayMeasure.getNearestPoint()[0], 11.2f) ||
	    !nearly_equal(overlayMeasure.getNearestPoint()[1], 0.0f))
	FAIL("overlay measure should report nearest overlay segment identity");

    SoBRLMeasureAction adcMeasure;
    adcMeasure.setQueryPoint(SbVec3f(22.0f, 0.3f, 0.0f));
    adcMeasure.apply(adc);
    if (!adcMeasure.hasSegments() ||
	    adcMeasure.getShapeCount() != 1 ||
	    adcMeasure.getSegmentCount() != 4)
	FAIL("ADC measure should count its own field-derived segments");
    if (!adcMeasure.hasNearestSegment() ||
	    strcmp(adcMeasure.getNearestPath().getString(), "overlay::adc") != 0 ||
	    !nearly_equal(adcMeasure.getNearestPoint()[0], 22.0f) ||
	    !nearly_equal(adcMeasure.getNearestPoint()[1], 0.0f))
	FAIL("ADC measure should report nearest angle-distance segment identity");

    SoBRLExportAction overlayExport;
    overlayExport.apply(root);
    if (overlayExport.getLineCount() != 17)
	FAIL("overlay export should collect grid, axes, and ADC line records");
    if (export_path_count(overlayExport, "overlay::work-grid") != 10 ||
	    export_path_count(overlayExport, "overlay::model-axes") != 3 ||
	    export_path_count(overlayExport, "overlay::adc") != 4)
	FAIL("overlay export should preserve per-overlay identities");
    bbox = overlayExport.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], -2.0f) ||
	    !nearly_equal(bbox.getMax()[0], 24.0f))
	FAIL("overlay export should report transformed export bounds");

    grid->divisions = 1;
    gridShape = grid->rebuildGeometry();
    if (!gridShape || gridShape->getSegmentCount() != 6)
	FAIL("grid overlay should rebuild line geometry after field changes");
    SoGetBoundingBoxAction gridBBoxAction(viewport);
    gridBBoxAction.apply(grid);
    bbox = gridBBoxAction.getBoundingBox();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], -1.0f) ||
	    !nearly_equal(bbox.getMax()[0], 1.0f) ||
	    !nearly_equal(bbox.getMin()[1], -1.0f) ||
	    !nearly_equal(bbox.getMax()[1], 1.0f))
	FAIL("grid overlay rebuild should update field-derived bounds");

    root->unref();

    root = new SoSeparator;
    root->ref();

    struct bg_line_layer_builder *diagnosticBuilder = bg_line_layer_builder_create();
    if (!diagnosticBuilder)
	FAIL("failed to allocate diagnostic line-layer builder");

    point_t lp;
    VSET(lp, 100.0, 0.0, 0.0);
    if (!bg_line_layer_builder_add(diagnosticBuilder, 255, 0, 0, lp, BG_GEOMETRY_LINE_MOVE))
	FAIL("failed to add diagnostic red line-layer move");
    VSET(lp, 102.0, 0.0, 0.0);
    if (!bg_line_layer_builder_add(diagnosticBuilder, 255, 0, 0, lp, BG_GEOMETRY_LINE_DRAW))
	FAIL("failed to add diagnostic red line-layer draw");
    VSET(lp, 101.0, -1.0, 0.0);
    if (!bg_line_layer_builder_add(diagnosticBuilder, 255, 0, 0, lp, BG_GEOMETRY_LINE_MOVE))
	FAIL("failed to add diagnostic red vertical move");
    VSET(lp, 101.0, 1.0, 0.0);
    if (!bg_line_layer_builder_add(diagnosticBuilder, 255, 0, 0, lp, BG_GEOMETRY_LINE_DRAW))
	FAIL("failed to add diagnostic red vertical draw");
    VSET(lp, 101.0, 0.0, 0.0);
    if (!bg_line_layer_builder_add(diagnosticBuilder, 255, 0, 0, lp, BG_GEOMETRY_POINT_DRAW))
	FAIL("failed to add diagnostic red point");
    VSET(lp, 104.0, 0.0, 0.0);
    if (!bg_line_layer_builder_add(diagnosticBuilder, 0, 255, 255, lp, BG_GEOMETRY_LINE_MOVE))
	FAIL("failed to add diagnostic cyan line-layer move");
    VSET(lp, 106.0, 0.0, 0.0);
    if (!bg_line_layer_builder_add(diagnosticBuilder, 0, 255, 255, lp, BG_GEOMETRY_LINE_DRAW))
	FAIL("failed to add diagnostic cyan line-layer draw");

    SoBRLLineLayerOverlay *lineOverlay = new SoBRLLineLayerOverlay;
    lineOverlay->overlayId = "overlay::diagnostic";
    lineOverlay->sourceId = 55;
    if (lineOverlay->rebuildGeometry(diagnosticBuilder) != 2 ||
	    lineOverlay->getLayerShapeCount() != 2 ||
	    lineOverlay->getPointCount() != 7)
	FAIL("line-layer overlay should realize one Obol shape per color layer");
    SoBRLVListShape *redDiagnostic = lineOverlay->getLayerShape(0);
    SoBRLVListShape *cyanDiagnostic = lineOverlay->getLayerShape(1);
    if (!redDiagnostic || !cyanDiagnostic)
	FAIL("line-layer overlay should expose realized color-layer shapes");
    if (strcmp(redDiagnostic->sourcePath.getValue().getString(),
		"overlay::diagnostic/rgb_255_000_000") != 0 ||
	    strcmp(redDiagnostic->sourceName.getValue().getString(),
		"rgb_255_000_000") != 0 ||
	    strcmp(redDiagnostic->sourceType.getValue().getString(), "line-layer") != 0 ||
	    redDiagnostic->sourceId.getValue() != 55 ||
	    redDiagnostic->getSegmentCount() != 2 ||
	    !redDiagnostic->colorOverride.getValue() ||
	    !nearly_equal(redDiagnostic->color.getValue()[0], 1.0f) ||
	    !nearly_equal(redDiagnostic->color.getValue()[1], 0.0f) ||
	    !nearly_equal(redDiagnostic->color.getValue()[2], 0.0f))
	FAIL("line-layer overlay should preserve red layer identity, color, and geometry");
    if (strcmp(cyanDiagnostic->sourcePath.getValue().getString(),
		"overlay::diagnostic/rgb_000_255_255") != 0 ||
	    strcmp(cyanDiagnostic->sourceName.getValue().getString(),
		"rgb_000_255_255") != 0 ||
	    strcmp(cyanDiagnostic->sourceType.getValue().getString(), "line-layer") != 0 ||
	    cyanDiagnostic->sourceId.getValue() != 55 ||
	    cyanDiagnostic->getSegmentCount() != 1 ||
	    !cyanDiagnostic->colorOverride.getValue() ||
	    !nearly_equal(cyanDiagnostic->color.getValue()[0], 0.0f) ||
	    !nearly_equal(cyanDiagnostic->color.getValue()[1], 1.0f) ||
	    !nearly_equal(cyanDiagnostic->color.getValue()[2], 1.0f))
	FAIL("line-layer overlay should preserve cyan layer identity, color, and geometry");
    root->addChild(lineOverlay);

    SoGetBoundingBoxAction lineLayerBBoxAction(viewport);
    lineLayerBBoxAction.apply(root);
    bbox = lineLayerBBoxAction.getBoundingBox();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 100.0f) ||
	    !nearly_equal(bbox.getMax()[0], 106.0f) ||
	    !nearly_equal(bbox.getMin()[1], -1.0f) ||
	    !nearly_equal(bbox.getMax()[1], 1.0f))
	FAIL("line-layer overlay should contribute Obol bounds from diagnostic geometry");

    SoRayPickAction lineLayerPick(viewport);
    lineLayerPick.setRay(SbVec3f(101.0f, 0.5f, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    lineLayerPick.apply(root);
    pickedPoint = lineLayerPick.getPickedPoint();
    if (!pickedPoint)
	FAIL("line-layer overlay should be pickable");
    rawDetail = pickedPoint->getDetail(redDiagnostic);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("line-layer overlay pick should return BRL-CAD Obol detail");
    pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(),
		"overlay::diagnostic/rgb_255_000_000") != 0 ||
	    strcmp(pickDetail->getSourceType().getString(), "line-layer") != 0 ||
	    pickDetail->getSourceId() != 55 ||
	    pickDetail->getPrimitiveKind() != SoBRLPickDetail::LINE_SEGMENT)
	FAIL("line-layer overlay pick should preserve layer and primitive identity");

    SoBRLSnapAction lineLayerSnap;
    lineLayerSnap.setQueryPoint(SbVec3f(104.5f, 0.25f, 0.0f));
    lineLayerSnap.setTolerance(0.3f);
    lineLayerSnap.apply(root);
    if (!lineLayerSnap.hasCandidate() ||
	    strcmp(lineLayerSnap.getPath().getString(),
		"overlay::diagnostic/rgb_000_255_255") != 0 ||
	    !nearly_equal(lineLayerSnap.getPoint()[0], 104.5f) ||
	    !nearly_equal(lineLayerSnap.getPoint()[1], 0.0f))
	FAIL("line-layer overlay snap should preserve color-layer identity");

    SoBRLMeasureAction lineLayerMeasure;
    lineLayerMeasure.setQueryPoint(SbVec3f(104.5f, 0.25f, 0.0f));
    lineLayerMeasure.apply(root);
    if (!lineLayerMeasure.hasSegments() ||
	    lineLayerMeasure.getShapeCount() != 2 ||
	    lineLayerMeasure.getSegmentCount() != 3 ||
	    !nearly_equal(lineLayerMeasure.getTotalLength(), 6.0f) ||
	    strcmp(lineLayerMeasure.getNearestPath().getString(),
		"overlay::diagnostic/rgb_000_255_255") != 0)
	FAIL("line-layer overlay measure should traverse realized layer geometry");

    SoBRLExportAction lineLayerExport;
    lineLayerExport.apply(root);
    if (lineLayerExport.getLineCount() != 3 ||
	    export_path_count(lineLayerExport, "overlay::diagnostic/rgb_255_000_000") != 2 ||
	    export_path_count(lineLayerExport, "overlay::diagnostic/rgb_000_255_255") != 1)
	FAIL("line-layer overlay export should preserve color-layer paths");
    const SoBRLExportAction::LineRecord *redLine =
	export_line_with_path(lineLayerExport, "overlay::diagnostic/rgb_255_000_000");
    if (!redLine ||
	    strcmp(redLine->sourceType.getString(), "line-layer") != 0 ||
	    redLine->sourceId != 55 ||
	    !redLine->colorOverride ||
	    !nearly_equal(redLine->color[0], 1.0f) ||
	    !nearly_equal(redLine->color[1], 0.0f) ||
	    !nearly_equal(redLine->color[2], 0.0f))
	FAIL("line-layer overlay export should carry layer color and identity");

    lineOverlay->selectable = FALSE;
    if (lineOverlay->rebuildGeometry(diagnosticBuilder) != 2)
	FAIL("line-layer overlay should rebuild after selectable field changes");
    redDiagnostic = lineOverlay->getLayerShape(0);
    SoRayPickAction unselectableLineLayerPick(viewport);
    unselectableLineLayerPick.setRay(SbVec3f(101.0f, 0.5f, 5.0f),
	    SbVec3f(0.0f, 0.0f, -1.0f));
    unselectableLineLayerPick.apply(root);
    if (unselectableLineLayerPick.getPickedPoint())
	FAIL("unselectable line-layer overlay should not be pickable");
    SoBRLExportAction unselectableLineLayerExport;
    unselectableLineLayerExport.apply(root);
    if (unselectableLineLayerExport.getLineCount() != 3)
	FAIL("unselectable line-layer overlay should remain exportable");

    lineOverlay->visible = FALSE;
    if (lineOverlay->rebuildGeometry(diagnosticBuilder) != 0 ||
	    lineOverlay->getLayerShapeCount() != 0)
	FAIL("hidden line-layer overlay should clear realized geometry");
    SoGetBoundingBoxAction hiddenLineLayerBBoxAction(viewport);
    hiddenLineLayerBBoxAction.apply(root);
    if (!hiddenLineLayerBBoxAction.getBoundingBox().isEmpty())
	FAIL("hidden line-layer overlay should not contribute bounds");

    SoSeparator *controllerOverlayRoot = new SoSeparator;
    controllerOverlayRoot->ref();
    {
	BRLObolViewController overlayController(controllerOverlayRoot, NULL);
	overlayController.clearRenderRequest();
	if (overlayController.replaceLineLayerOverlay("gqa::overlaps",
		    diagnosticBuilder, 88, FALSE) != 2 ||
		controllerOverlayRoot->getNumChildren() != 1 ||
		!overlayController.isRenderRequested() ||
		strcmp(overlayController.getRenderReason().getString(),
		    "line-layer-overlay") != 0)
	    FAIL("view controller should publish a line-layer overlay into the Obol scene");

	SoNode *overlayNode = controllerOverlayRoot->getChild(0);
	if (!overlayNode ||
		!overlayNode->isOfType(SoBRLLineLayerOverlay::getClassTypeId()))
	    FAIL("view controller line-layer publication should create an Obol overlay node");
	SoBRLLineLayerOverlay *controllerOverlay =
	    static_cast<SoBRLLineLayerOverlay *>(overlayNode);
	if (strcmp(controllerOverlay->overlayId.getValue().getString(),
		    "gqa::overlaps") != 0 ||
		controllerOverlay->sourceId.getValue() != 88 ||
		controllerOverlay->getLayerShapeCount() != 2 ||
		controllerOverlay->getPointCount() != 7)
	    FAIL("view controller line-layer publication should preserve overlay identity");
	SoBRLVListShape *controllerLayer = controllerOverlay->getLayerShape(0);
	if (!controllerLayer ||
		controllerLayer->selectable.getValue() ||
		strcmp(controllerLayer->sourcePath.getValue().getString(),
		    "gqa::overlaps/rgb_255_000_000") != 0)
	    FAIL("view controller line-layer publication should apply selectable state and layer paths");

	if (overlayController.replaceLineLayerOverlay("gqa::overlaps",
		    diagnosticBuilder, 89, TRUE) != 2 ||
		controllerOverlayRoot->getNumChildren() != 1)
	    FAIL("view controller line-layer publication should replace an existing overlay");
	controllerOverlay = static_cast<SoBRLLineLayerOverlay *>(
		controllerOverlayRoot->getChild(0));
	controllerLayer = controllerOverlay->getLayerShape(0);
	if (controllerOverlay->sourceId.getValue() != 89 ||
		!controllerLayer ||
		!controllerLayer->selectable.getValue())
	    FAIL("view controller line-layer replacement should refresh overlay fields");

	if (overlayController.replaceLineLayerOverlay("gqa::overlaps",
		    NULL, 0, TRUE) != 0 ||
		controllerOverlayRoot->getNumChildren() != 0)
	    FAIL("view controller line-layer publication should remove named overlays");
    }
    controllerOverlayRoot->unref();

    SoBRLHUDLabelOverlay *labelOverlay = new SoBRLHUDLabelOverlay;
    labelOverlay->ref();
    labelOverlay->labelId = "diagnostic::label";
    labelOverlay->sourceId = 77;
    labelOverlay->text = "overlap candidate";
    labelOverlay->position = SbVec2f(12.0f, 34.0f);
    labelOverlay->color = SbColor(0.25f, 0.5f, 0.75f);
    labelOverlay->fontSize = 13.0f;
    SoHUDLabel *hudLabel = labelOverlay->rebuildGeometry();
    if (!hudLabel ||
	    !labelOverlay->getHUDKit() ||
	    labelOverlay->getHUDLabel() != hudLabel)
	FAIL("HUD label overlay should rebuild direct Obol HUD nodes");
    if (strcmp(labelOverlay->labelId.getValue().getString(),
		"diagnostic::label") != 0 ||
	    labelOverlay->sourceId.getValue() != 77 ||
	    strcmp(hudLabel->string[0].getString(), "overlap candidate") != 0 ||
	    !nearly_equal(hudLabel->position.getValue()[0], 12.0f) ||
	    !nearly_equal(hudLabel->position.getValue()[1], 34.0f) ||
	    !nearly_equal(hudLabel->color.getValue()[0], 0.25f) ||
	    !nearly_equal(hudLabel->color.getValue()[1], 0.5f) ||
	    !nearly_equal(hudLabel->color.getValue()[2], 0.75f) ||
	    !nearly_equal(hudLabel->fontSize.getValue(), 13.0f))
	FAIL("HUD label overlay should preserve label identity and presentation fields");
    labelOverlay->visible = FALSE;
    if (labelOverlay->rebuildGeometry() ||
	    labelOverlay->getHUDKit() ||
	    labelOverlay->getNumChildren() != 0)
	FAIL("hidden HUD label overlay should clear direct Obol HUD nodes");
    labelOverlay->unref();

    SoSeparator *controllerLabelRoot = new SoSeparator;
    controllerLabelRoot->ref();
    {
	BRLObolViewController labelController(controllerLabelRoot, NULL);
	labelController.clearRenderRequest();
	if (labelController.replaceHUDLabelOverlay("gqa::label",
		    "overlap 1", SbVec2f(20.0f, 30.0f),
		    SbColor(1.0f, 0.25f, 0.0f), 10.0f, 101) != 1 ||
		controllerLabelRoot->getNumChildren() != 1 ||
		!labelController.isRenderRequested() ||
		strcmp(labelController.getRenderReason().getString(),
		    "hud-label-overlay") != 0)
	    FAIL("view controller should publish a HUD label overlay into the Obol scene");

	SoNode *labelNode = controllerLabelRoot->getChild(0);
	if (!labelNode ||
		!labelNode->isOfType(SoBRLHUDLabelOverlay::getClassTypeId()))
	    FAIL("view controller HUD label publication should create a BRL-CAD label node");
	SoBRLHUDLabelOverlay *controllerLabel =
	    static_cast<SoBRLHUDLabelOverlay *>(labelNode);
	hudLabel = controllerLabel->getHUDLabel();
	if (!hudLabel ||
		controllerLabel->sourceId.getValue() != 101 ||
		strcmp(controllerLabel->labelId.getValue().getString(),
		    "gqa::label") != 0 ||
		strcmp(hudLabel->string[0].getString(), "overlap 1") != 0 ||
		!nearly_equal(hudLabel->position.getValue()[0], 20.0f) ||
		!nearly_equal(hudLabel->fontSize.getValue(), 10.0f))
	    FAIL("view controller HUD label should preserve identity, text, and placement");

	if (labelController.replaceHUDLabelOverlay("gqa::label",
		    "overlap 2", SbVec2f(22.0f, 32.0f),
		    SbColor(0.0f, 1.0f, 0.0f), 14.0f, 102) != 1 ||
		controllerLabelRoot->getNumChildren() != 1)
	    FAIL("view controller HUD label publication should replace existing labels");
	controllerLabel = static_cast<SoBRLHUDLabelOverlay *>(
	    controllerLabelRoot->getChild(0));
	hudLabel = controllerLabel->getHUDLabel();
	if (!hudLabel ||
		controllerLabel->sourceId.getValue() != 102 ||
		strcmp(hudLabel->string[0].getString(), "overlap 2") != 0 ||
		!nearly_equal(hudLabel->position.getValue()[0], 22.0f) ||
		!nearly_equal(hudLabel->fontSize.getValue(), 14.0f))
	    FAIL("view controller HUD label replacement should refresh fields");

	if (labelController.replaceHUDLabelOverlay("gqa::label",
		    NULL, SbVec2f(0.0f, 0.0f), SbColor(1.0f, 1.0f, 1.0f)) != 1 ||
		controllerLabelRoot->getNumChildren() != 0)
	    FAIL("view controller HUD label publication should remove named labels");
    }
    controllerLabelRoot->unref();

    bg_line_layer_builder_free(diagnosticBuilder);
    root->unref();

    char dbpath[MAXPATHLEN] = {0};
    char lodCacheDir[MAXPATHLEN] = {0};
    bu_dir(lodCacheDir, MAXPATHLEN, BU_DIR_CURR, "brlobol_prototype_lod_cache", NULL);
    bu_dirclear(lodCacheDir);
    bu_mkdir(lodCacheDir);
    bu_setenv("BU_DIR_CACHE", lodCacheDir, 1);

    if (!write_test_db(dbpath, MAXPATHLEN))
	FAIL("failed to create prototype test database");

    struct db_i *dbip = NULL;
    if (!open_database(dbpath, &dbip)) {
	bu_file_delete(dbpath);
	FAIL("failed to open prototype test database");
    }

    root = new SoSeparator;
    root->ref();

    SoBRLDatabaseSource *missingSource = new SoBRLDatabaseSource;
    missingSource->setDatabase(dbip);
    missingSource->path = "missing.s";
    missingSource->drawMode = SoBRLDatabaseSource::SHADED;
    missingSource->sourceRevision = 1;
    root->addChild(missingSource);

    SoBRLRealizeAction failedRealize;
    failedRealize.apply(root);
    if (failedRealize.getFailedSourceCount() != 1 ||
	    failedRealize.getRealizedSourceCount() != 0 ||
	    missingSource->realizationStatus.getValue() != SoBRLDatabaseSource::FAILED)
	FAIL("database realization failure should be counted and mark the source failed");
    if (missingSource->realizationDiagnostic.getValue().getLength() == 0 ||
	    !strstr(missingSource->realizationDiagnostic.getValue().getString(), "missing.s") ||
	    failedRealize.getDiagnostics().getLength() == 0 ||
	    !strstr(failedRealize.getDiagnostics().getString(), "missing.s"))
	FAIL("database realization failure should report source and action diagnostics");

    SoBRLSceneController failedController(root);
    if (failedController.realizePending() ||
	    failedController.getLastFailedSourceCount() != 1 ||
	    failedController.getLastDiagnostics().getLength() == 0 ||
	    !strstr(failedController.getLastDiagnostics().getString(), "missing.s"))
	FAIL("scene controller should expose realization failure diagnostics");

    root->unref();

    root = new SoSeparator;
    root->ref();

    source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "box.s";
    source->sourceRevision = 7;
    root->addChild(source);

    if (source->getDatabase() != dbip)
	FAIL("database source should retain application-provided db_i pointer");

    SoBRLRealizeAction dbRealize;
    dbRealize.apply(root);
    if (dbRealize.getRealizedSourceCount() != 1)
	FAIL("database-backed realization should realize one source");
    if (source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED)
	FAIL("database-backed source realization status should be REALIZED");

    shape = source->getRealizedShape();
    if (!shape)
	FAIL("database-backed source should produce an Obol vlist shape");
    if (source->getRealizedShapeCount() != 1)
	FAIL("direct primitive realization should produce one leaf Obol shape");
    if (shape->getSegmentCount() <= 4)
	FAIL("database-backed source should not fall back to the synthetic square");
    if (strcmp(shape->sourcePath.getValue().getString(), "/box.s") != 0)
	FAIL("direct primitive shape should preserve its full database path");
    if (strcmp(shape->sourceName.getValue().getString(), "box.s") != 0 ||
	    strcmp(shape->sourceType.getValue().getString(), "arb8") != 0 ||
	    shape->sourceId.getValue() != 7 ||
	    shape->materialColorValid.getValue() ||
	    shape->regionId.getValue() != 0)
	FAIL("direct primitive shape should preserve primitive identity fields");

    SoGetBoundingBoxAction dbBBoxAction(viewport);
    dbBBoxAction.apply(root);
    bbox = dbBBoxAction.getBoundingBox();
    if (bbox.isEmpty())
	FAIL("database-backed source should contribute a bounding box");
    if (!nearly_equal(bbox.getMin()[0], -2.0f) ||
	    !nearly_equal(bbox.getMax()[0], 3.0f) ||
	    !nearly_equal(bbox.getMin()[1], -3.0f) ||
	    !nearly_equal(bbox.getMax()[1], 4.0f) ||
	    !nearly_equal(bbox.getMin()[2], -4.0f) ||
	    !nearly_equal(bbox.getMax()[2], 5.0f))
	FAIL("database-backed bounding box should come from real .g wireframe geometry");

    SoSeparator *controllerDbRoot = new SoSeparator;
    controllerDbRoot->ref();
    controllerDbRoot->addChild(new SoSeparator);
    {
	BRLObolViewController dbController(controllerDbRoot, NULL);
	dbController.clearRenderRequest();
	if (dbController.replaceDatabaseSource("box.s", dbip,
		    SoBRLDatabaseSource::WIREFRAME, 71) != 1 ||
		dbController.getDatabaseSourceCount() != 1 ||
		controllerDbRoot->getNumChildren() != 2 ||
		!dbController.isRenderRequested() ||
		strcmp(dbController.getRenderReason().getString(),
		    "database-source") != 0)
	    FAIL("view controller should publish database sources into the Obol scene");

	SoBRLDatabaseSource *controllerSource = dbController.getDatabaseSource(0);
	if (!controllerSource ||
		controllerSource->getDatabase() != dbip ||
		strcmp(controllerSource->path.getValue().getString(), "box.s") != 0 ||
		controllerSource->drawMode.getValue() != SoBRLDatabaseSource::WIREFRAME ||
		controllerSource->sourceRevision.getValue() != 71 ||
		!controllerSource->needsRealization())
	    FAIL("view controller database source should preserve path, database, mode, and revision");
	if (!dbController.realizePending() ||
		controllerSource->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
		controllerSource->getRealizedShapeCount() != 1)
	    FAIL("view controller database source should realize through the scene controller");
	SoDB::getSensorManager()->processDelayQueue(TRUE);
	if (controllerSource->needsRealization() ||
		controllerSource->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED)
	    FAIL("view controller database source replacement should not leave delayed stale sensors after realization");

	if (dbController.replaceDatabaseSource("/box.s", dbip,
		    SoBRLDatabaseSource::WIREFRAME, 72) != 1 ||
		dbController.getDatabaseSourceCount() != 1 ||
		controllerDbRoot->getNumChildren() != 2 ||
		dbController.getDatabaseSource(0)->sourceRevision.getValue() != 72)
	    FAIL("view controller database source replacement should match slash-equivalent paths");

	if (dbController.removeDatabaseSource("box.s") != 1 ||
		dbController.getDatabaseSourceCount() != 0 ||
		controllerDbRoot->getNumChildren() != 1)
	    FAIL("view controller database source removal should leave non-source scene nodes");

	if (dbController.replaceDatabaseSource("box.s", dbip,
		    SoBRLDatabaseSource::WIREFRAME, 73) != 1 ||
		dbController.replaceDatabaseSource("ball.s", dbip,
		    SoBRLDatabaseSource::SHADED, 74) != 1 ||
		dbController.getDatabaseSourceCount() != 2 ||
		dbController.clearDatabaseSources() != 2 ||
		dbController.getDatabaseSourceCount() != 0 ||
		controllerDbRoot->getNumChildren() != 1)
	    FAIL("view controller should clear only database source scene nodes");
    }
    controllerDbRoot->unref();

    source->drawMode = SoBRLDatabaseSource::SHADED;
    if (!source->needsRealization() ||
	    source->staleReason.getValue() != SoBRLDatabaseSource::STALE_DRAW)
	FAIL("database source draw-mode change should invalidate only draw realization state");

    SoBRLRealizeAction dbMeshRealize;
    dbMeshRealize.apply(root);
    if (dbMeshRealize.getRealizedSourceCount() != 1)
	FAIL("database-backed shaded realization should realize one source");
    if (source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED)
	FAIL("database-backed shaded source realization status should be REALIZED");
    if (source->getRealizedShapeCount() != 0)
	FAIL("database-backed shaded realization should replace stale wireframe children");
    if (source->getRealizedMeshCount() != 1)
	FAIL("database-backed shaded realization should produce one Obol mesh shape");

    mesh = source->getRealizedMesh();
    if (!mesh)
	FAIL("database-backed shaded source should expose an Obol mesh shape");
    if (mesh->getTriangleCount() != 12)
	FAIL("database-backed ARB mesh should contain triangulated faces");
    if (strcmp(mesh->sourcePath.getValue().getString(), "/box.s") != 0)
	FAIL("database-backed mesh should preserve its full database path");
    if (strcmp(mesh->sourceName.getValue().getString(), "box.s") != 0 ||
	    strcmp(mesh->sourceType.getValue().getString(), "arb8") != 0 ||
	    mesh->sourceId.getValue() != 7 ||
	    mesh->materialColorValid.getValue() ||
	    mesh->regionId.getValue() != 0)
	FAIL("database-backed mesh should preserve primitive identity fields");

    SoGetBoundingBoxAction dbMeshBBoxAction(viewport);
    dbMeshBBoxAction.apply(root);
    bbox = dbMeshBBoxAction.getBoundingBox();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], -2.0f) ||
	    !nearly_equal(bbox.getMax()[0], 3.0f) ||
	    !nearly_equal(bbox.getMin()[1], -3.0f) ||
	    !nearly_equal(bbox.getMax()[1], 4.0f) ||
	    !nearly_equal(bbox.getMin()[2], -4.0f) ||
	    !nearly_equal(bbox.getMax()[2], 5.0f))
	FAIL("database-backed mesh bounding box should come from Obol mesh fields");

    SoRayPickAction dbMeshPick(viewport);
    dbMeshPick.setRay(SbVec3f(0.0f, 0.0f, 10.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    dbMeshPick.apply(root);
    pickedPoint = dbMeshPick.getPickedPoint();
    if (!pickedPoint)
	FAIL("database-backed mesh ray pick should hit shaded geometry");
    rawDetail = pickedPoint->getDetail(mesh);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("database-backed mesh pick should return a BRL-CAD Obol pick detail");
    pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(), "/box.s") != 0 ||
	    pickDetail->getPrimitiveKind() != SoBRLPickDetail::FACE)
	FAIL("database-backed mesh pick detail should preserve face identity");
    if (strcmp(pickDetail->getSourceName().getString(), "box.s") != 0 ||
	    strcmp(pickDetail->getSourceType().getString(), "arb8") != 0 ||
	    pickDetail->getSourceId() != 7 ||
	    pickDetail->hasMaterialColor() ||
	    pickDetail->getRegionId() != 0)
	FAIL("database-backed mesh pick detail should preserve primitive identity fields");

    SoBRLExportAction dbMeshExport;
    dbMeshExport.apply(root);
    if (dbMeshExport.getLineCount() != 0 ||
	    dbMeshExport.getTriangleCount() != 12)
	FAIL("database-backed shaded export should collect mesh triangles, not wire lines");
    if (export_triangle_path_count(dbMeshExport, "/box.s") != 12)
	FAIL("database-backed shaded export should preserve mesh path identity");
    const SoBRLExportAction::TriangleRecord *boxTriangle =
	export_triangle_with_path(dbMeshExport, "/box.s");
    if (!boxTriangle ||
	    strcmp(boxTriangle->sourceName.getString(), "box.s") != 0 ||
	    strcmp(boxTriangle->sourceType.getString(), "arb8") != 0 ||
	    boxTriangle->sourceId != 7 ||
	    boxTriangle->materialColorValid ||
	    boxTriangle->regionId != 0)
	FAIL("database-backed shaded export should preserve primitive identity");
    bbox = dbMeshExport.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], -2.0f) ||
	    !nearly_equal(bbox.getMax()[0], 3.0f) ||
	    !nearly_equal(bbox.getMin()[2], -4.0f) ||
	    !nearly_equal(bbox.getMax()[2], 5.0f))
	FAIL("database-backed shaded export should report mesh bounds");

    SoBRLMeasureAction dbMeshMeasure;
    dbMeshMeasure.setQueryPoint(SbVec3f(0.0f, 0.0f, 10.0f));
    dbMeshMeasure.apply(root);
    if (!dbMeshMeasure.hasFaces() ||
	    dbMeshMeasure.getShapeCount() != 1 ||
	    dbMeshMeasure.getTriangleCount() != 12 ||
	    !nearly_equal(dbMeshMeasure.getSurfaceArea(), 286.0f))
	FAIL("database-backed shaded measure should report mesh face metrics");
    bbox = dbMeshMeasure.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], -2.0f) ||
	    !nearly_equal(bbox.getMax()[0], 3.0f) ||
	    !nearly_equal(bbox.getMin()[1], -3.0f) ||
	    !nearly_equal(bbox.getMax()[1], 4.0f) ||
	    !nearly_equal(bbox.getMin()[2], -4.0f) ||
	    !nearly_equal(bbox.getMax()[2], 5.0f))
	FAIL("database-backed shaded measure should report mesh field bounds");
    if (!dbMeshMeasure.hasNearestPrimitive() ||
	    dbMeshMeasure.getNearestPrimitiveKind() != SoBRLMeasureAction::FACE ||
	    strcmp(dbMeshMeasure.getNearestPath().getString(), "/box.s") != 0 ||
	    !nearly_equal(dbMeshMeasure.getNearestPoint()[0], 0.0f) ||
	    !nearly_equal(dbMeshMeasure.getNearestPoint()[1], 0.0f) ||
	    !nearly_equal(dbMeshMeasure.getNearestPoint()[2], 5.0f))
	FAIL("database-backed shaded measure should report nearest mesh face identity");

    source->path = "ball.s";
    source->sourceRevision = 17;
    if (!source->needsRealization() ||
	    source->staleReason.getValue() == SoBRLDatabaseSource::STALE_NONE)
	FAIL("database source path change should invalidate shaded realization state");

    SoBRLRealizeAction sphereMeshRealize;
    sphereMeshRealize.apply(root);
    if (sphereMeshRealize.getRealizedSourceCount() != 1)
	FAIL("database-backed sphere shaded realization should realize one source");
    if (source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED)
	FAIL("database-backed sphere shaded source realization status should be REALIZED");
    if (source->getRealizedMeshCount() != 1)
	FAIL("database-backed sphere shaded realization should produce one Obol mesh shape");

    mesh = source->getRealizedMesh();
    if (!mesh)
	FAIL("database-backed sphere source should expose an Obol mesh shape");
    if (mesh->getTriangleCount() <= 12)
	FAIL("database-backed sphere mesh should come from generalized tessellation, not ARB triangulation");
    if (strcmp(mesh->sourcePath.getValue().getString(), "/ball.s") != 0)
	FAIL("database-backed tessellated mesh should preserve its full database path");
    int defaultSphereTriangleCount = mesh->getTriangleCount();

    SoGetBoundingBoxAction sphereBBoxAction(viewport);
    sphereBBoxAction.apply(root);
    bbox = sphereBBoxAction.getBoundingBox();
    if (bbox.isEmpty() ||
	    bbox.getMin()[0] > -1.9f ||
	    bbox.getMax()[0] < 1.9f ||
	    bbox.getMin()[1] > -1.9f ||
	    bbox.getMax()[1] < 1.9f ||
	    bbox.getMin()[2] > -1.9f ||
	    bbox.getMax()[2] < 1.9f)
	FAIL("database-backed tessellated sphere mesh should report spherical bounds");

    SoRayPickAction spherePick(viewport);
    spherePick.setRay(SbVec3f(0.0f, 0.0f, 5.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    spherePick.apply(root);
    pickedPoint = spherePick.getPickedPoint();
    if (!pickedPoint)
	FAIL("database-backed tessellated sphere pick should hit shaded geometry");
    rawDetail = pickedPoint->getDetail(mesh);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("database-backed tessellated sphere pick should return a BRL-CAD Obol pick detail");
    pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(), "/ball.s") != 0 ||
	    pickDetail->getPrimitiveKind() != SoBRLPickDetail::FACE)
	FAIL("database-backed tessellated sphere pick detail should preserve face identity");

    SoBRLExportAction sphereExport;
    sphereExport.apply(root);
    if (sphereExport.getLineCount() != 0 ||
	    sphereExport.getTriangleCount() != mesh->getTriangleCount() ||
	    export_triangle_path_count(sphereExport, "/ball.s") != mesh->getTriangleCount())
	FAIL("database-backed tessellated sphere export should collect mesh triangles with path identity");

    SoBRLMeasureAction sphereMeasure;
    sphereMeasure.setQueryPoint(SbVec3f(0.0f, 0.0f, 3.0f));
    sphereMeasure.apply(root);
    if (!sphereMeasure.hasFaces() ||
	    sphereMeasure.getShapeCount() != 1 ||
	    sphereMeasure.getTriangleCount() != mesh->getTriangleCount() ||
	    sphereMeasure.getSurfaceArea() <= 35.0f ||
	    sphereMeasure.getSurfaceArea() >= 70.0f)
	FAIL("database-backed tessellated sphere measure should report mesh face metrics");
    if (!sphereMeasure.hasNearestPrimitive() ||
	    sphereMeasure.getNearestPrimitiveKind() != SoBRLMeasureAction::FACE ||
	    strcmp(sphereMeasure.getNearestPath().getString(), "/ball.s") != 0 ||
	    sphereMeasure.getNearestPoint()[2] < 1.8f ||
	    sphereMeasure.getNearestPoint()[2] > 2.1f)
	FAIL("database-backed tessellated sphere measure should report nearest mesh face identity");

    SoBRLSnapAction sphereSnap;
    sphereSnap.setEnabledKinds(SoBRLSnapAction::FACE_NEAREST);
    sphereSnap.setQueryPoint(SbVec3f(0.0f, 0.0f, 3.0f));
    sphereSnap.setTolerance(1.25f);
    sphereSnap.apply(root);
    if (!sphereSnap.hasCandidate() ||
	    sphereSnap.getKind() != SoBRLSnapAction::FACE_NEAREST ||
	    strcmp(sphereSnap.getPath().getString(), "/ball.s") != 0 ||
	    sphereSnap.getPoint()[2] < 1.8f ||
	    sphereSnap.getPoint()[2] > 2.1f)
	FAIL("database-backed tessellated sphere snap should report nearest mesh face identity");

    source->tessellationRelTol = 0.25f;
    if (!source->needsRealization() ||
	    !(source->staleReason.getValue() & SoBRLDatabaseSource::STALE_TESSELLATION))
	FAIL("database source tessellation tolerance changes should invalidate shaded realization state");

    SoBRLRealizeAction coarseSphereRealize;
    coarseSphereRealize.apply(root);
    if (coarseSphereRealize.getRealizedSourceCount() != 1)
	FAIL("database-backed coarse sphere realization should realize one source");
    mesh = source->getRealizedMesh();
    if (!mesh || mesh->getTriangleCount() <= 0)
	FAIL("database-backed coarse sphere realization should produce mesh geometry");
    if (mesh->getTriangleCount() >= defaultSphereTriangleCount)
	FAIL("coarser tessellation tolerance should reduce generated sphere triangle count");

    SoBRLExportAction coarseSphereExport;
    coarseSphereExport.apply(root);
    if (coarseSphereExport.getTriangleCount() != mesh->getTriangleCount() ||
	    export_triangle_path_count(coarseSphereExport, "/ball.s") != mesh->getTriangleCount())
	FAIL("database-backed coarse sphere export should use re-realized mesh geometry");

    source->path = "ell.s";
    source->sourceRevision = 18;
    source->tessellationRelTol = 0.05f;
    if (!source->needsRealization() ||
	    source->staleReason.getValue() == SoBRLDatabaseSource::STALE_NONE)
	FAIL("database source curved primitive/tolerance changes should invalidate realization state");

    SoBRLRealizeAction ellipsoidRealize;
    ellipsoidRealize.apply(root);
    if (ellipsoidRealize.getRealizedSourceCount() != 1)
	FAIL("database-backed ellipsoid shaded realization should realize one source");
    mesh = source->getRealizedMesh();
    if (!mesh || mesh->getTriangleCount() <= 12)
	FAIL("database-backed ellipsoid mesh should come from generalized tessellation");
    if (strcmp(mesh->sourcePath.getValue().getString(), "/ell.s") != 0)
	FAIL("database-backed tessellated ellipsoid should preserve its full database path");

    SoGetBoundingBoxAction ellipsoidBBoxAction(viewport);
    ellipsoidBBoxAction.apply(root);
    bbox = ellipsoidBBoxAction.getBoundingBox();
    if (bbox.isEmpty() ||
	    bbox.getMin()[0] > 6.1f ||
	    bbox.getMax()[0] < 9.9f ||
	    bbox.getMin()[1] > -0.9f ||
	    bbox.getMax()[1] < 0.9f ||
	    bbox.getMin()[2] > -2.9f ||
	    bbox.getMax()[2] < 2.9f)
	FAIL("database-backed tessellated ellipsoid mesh should report ellipsoid bounds");

    SoBRLExportAction ellipsoidExport;
    ellipsoidExport.apply(root);
    if (ellipsoidExport.getTriangleCount() != mesh->getTriangleCount() ||
	    export_triangle_path_count(ellipsoidExport, "/ell.s") != mesh->getTriangleCount())
	FAIL("database-backed tessellated ellipsoid export should collect mesh triangles with path identity");

    SoBRLMeasureAction ellipsoidMeasure;
    ellipsoidMeasure.setQueryPoint(SbVec3f(8.0f, 0.0f, 5.0f));
    ellipsoidMeasure.apply(root);
    if (!ellipsoidMeasure.hasFaces() ||
	    ellipsoidMeasure.getTriangleCount() != mesh->getTriangleCount() ||
	    ellipsoidMeasure.getSurfaceArea() <= 20.0f)
	FAIL("database-backed tessellated ellipsoid measure should report mesh face metrics");

    source->path = "tgc.s";
    source->sourceRevision = 19;
    if (!source->needsRealization())
	FAIL("database source TGC path change should invalidate shaded realization state");

    SoBRLRealizeAction tgcRealize;
    tgcRealize.apply(root);
    if (tgcRealize.getRealizedSourceCount() != 1)
	FAIL("database-backed TGC shaded realization should realize one source");
    mesh = source->getRealizedMesh();
    if (!mesh || mesh->getTriangleCount() <= 12)
	FAIL("database-backed TGC mesh should come from generalized tessellation");
    if (strcmp(mesh->sourcePath.getValue().getString(), "/tgc.s") != 0)
	FAIL("database-backed tessellated TGC should preserve its full database path");

    SoBRLExportAction tgcExport;
    tgcExport.apply(root);
    bbox = tgcExport.getBounds();
    if (tgcExport.getTriangleCount() != mesh->getTriangleCount() ||
	    export_triangle_path_count(tgcExport, "/tgc.s") != mesh->getTriangleCount() ||
	    bbox.isEmpty() ||
	    bbox.getMin()[0] > -1.4f ||
	    bbox.getMax()[0] < 1.4f ||
	    bbox.getMin()[1] > 7.1f ||
	    bbox.getMax()[1] < 8.9f ||
	    bbox.getMin()[2] > 0.1f ||
	    bbox.getMax()[2] < 3.9f)
	FAIL("database-backed tessellated TGC export should report mesh identity and bounds");

    source->path = "tor.s";
    source->sourceRevision = 20;
    if (!source->needsRealization())
	FAIL("database source torus path change should invalidate shaded realization state");

    SoBRLRealizeAction torusRealize;
    torusRealize.apply(root);
    if (torusRealize.getRealizedSourceCount() != 1)
	FAIL("database-backed torus shaded realization should realize one source");
    mesh = source->getRealizedMesh();
    if (!mesh || mesh->getTriangleCount() <= 24)
	FAIL("database-backed torus mesh should come from generalized tessellation");
    if (strcmp(mesh->sourcePath.getValue().getString(), "/tor.s") != 0)
	FAIL("database-backed tessellated torus should preserve its full database path");

    SoBRLMeasureAction torusMeasure;
    torusMeasure.setQueryPoint(SbVec3f(16.0f, 0.0f, 2.0f));
    torusMeasure.apply(root);
    bbox = torusMeasure.getBounds();
    if (!torusMeasure.hasFaces() ||
	    torusMeasure.getTriangleCount() != mesh->getTriangleCount() ||
	    torusMeasure.getSurfaceArea() <= 70.0f ||
	    bbox.isEmpty() ||
	    bbox.getMin()[0] > 12.5f ||
	    bbox.getMax()[0] < 19.5f ||
	    bbox.getMin()[1] > -3.6f ||
	    bbox.getMax()[1] < 3.3f ||
	    bbox.getMin()[2] > -0.7f ||
	    bbox.getMax()[2] < 0.7f)
	FAIL("database-backed tessellated torus measure should report mesh metrics and bounds");

    root->unref();

    if (!exercise_generated_primitive(dbip, "rpc.s", 1, 1, viewport))
	FAIL("database-backed RPC primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "rhc.s", 1, 1, viewport))
	FAIL("database-backed RHC primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "epa.s", 1, 1, viewport))
	FAIL("database-backed EPA primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "ehy.s", 1, 1, viewport))
	FAIL("database-backed EHY primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "eto.s", 1, 1, viewport))
	FAIL("database-backed ETO primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "hyp.s", 1, 1, viewport))
	FAIL("database-backed HYP primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "part.s", 1, 1, viewport))
	FAIL("database-backed particle primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "cline.s", 1, 1, viewport))
	FAIL("database-backed CLINE primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "hrt.s", 1, 0, viewport))
	FAIL("database-backed HRT wire primitive coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "hrt.s", "hrt", viewport))
	FAIL("database-backed HRT shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive(dbip, "pipe.s", 1, 1, viewport))
	FAIL("database-backed pipe primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "arbn.s", 1, 1, viewport))
	FAIL("database-backed ARBN primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "ars.s", 1, 1, viewport))
	FAIL("database-backed ARS primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "nmg.s", 1, 1, viewport))
	FAIL("database-backed NMG primitive coverage should pass");

    char polydbpath[MAXPATHLEN] = {0};
    if (!write_poly_v4_test_db(polydbpath, MAXPATHLEN))
	FAIL("failed to create prototype v4 POLY test database");
    struct db_i *polydbip = NULL;
    if (!open_database(polydbpath, &polydbip)) {
	bu_file_delete(polydbpath);
	FAIL("failed to open prototype v4 POLY test database");
    }
    if (!exercise_generated_primitive(polydbip, "poly.s", 1, 1, viewport))
	FAIL("database-backed legacy POLY primitive coverage should pass");
    db_close(polydbip);
    bu_file_delete(polydbpath);

    char hfdbpath[MAXPATHLEN] = {0};
    char hfdatapath[MAXPATHLEN] = {0};
    if (!write_hf_v4_test_db(hfdbpath, MAXPATHLEN, hfdatapath, MAXPATHLEN))
	FAIL("failed to create prototype v4 HF test database");
    struct db_i *hfdbip = NULL;
    if (!open_database(hfdbpath, &hfdbip)) {
	bu_file_delete(hfdbpath);
	bu_file_delete(hfdatapath);
	FAIL("failed to open prototype v4 HF test database");
    }
    if (!exercise_generated_primitive(hfdbip, "hf.s", 1, 0, viewport))
	FAIL("database-backed legacy HF wire primitive coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(hfdbip, "hf.s", "hf", viewport))
	FAIL("database-backed legacy HF shaded diagnostic coverage should pass");
    db_close(hfdbip);
    bu_file_delete(hfdbpath);
    bu_file_delete(hfdatapath);

    if (!exercise_generated_primitive(dbip, "rect.sketch", 1, 0, viewport))
	FAIL("database-backed sketch wire primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "extrude.s", 1, 1, viewport))
	FAIL("database-backed extrude primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "revolve.s", 1, 0, viewport))
	FAIL("database-backed revolve wire primitive coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "revolve.s", "revolve", viewport))
	FAIL("database-backed revolve shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive(dbip, "metaball.s", 1, 1, viewport, 0.01f))
	FAIL("database-backed metaball primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "superell.s", 1, 0, viewport))
	FAIL("database-backed superell wire primitive coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "superell.s", "superell", viewport))
	FAIL("database-backed superell shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive(dbip, "dsp.s", 1, 1, viewport))
	FAIL("database-backed DSP object-data primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "ebm.s", 1, 1, viewport))
	FAIL("database-backed EBM object-data primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "vol.s", 1, 1, viewport))
	FAIL("database-backed VOL object-data primitive coverage should pass");
    if (!exercise_generated_primitive(dbip, "brep.s", 1, 0, viewport))
	FAIL("database-backed BREP wire primitive coverage should pass");
    if (!exercise_generated_primitive_wire_diagnostic(dbip, "empty.brep", "brep", viewport))
	FAIL("database-backed empty BREP wire diagnostic coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "empty.brep", "brep", viewport))
	FAIL("database-backed empty BREP shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive(dbip, "bspline.s", 1, 0, viewport))
	FAIL("database-backed BSpline wire primitive coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "bspline.s", "bspline", viewport))
	FAIL("database-backed BSpline shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive_wire_diagnostic(dbip, "payload.binunif", "binunif", viewport))
	FAIL("database-backed binunif wire diagnostic coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "payload.binunif", "binunif", viewport))
	FAIL("database-backed binunif shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive_wire_diagnostic(dbip, "constraint.meta", "constrnt", viewport))
	FAIL("database-backed constraint wire diagnostic coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "constraint.meta", "constrnt", viewport))
	FAIL("database-backed constraint shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive_wire_diagnostic(dbip, "material.meta", "material", viewport))
	FAIL("database-backed material wire diagnostic coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "material.meta", "material", viewport))
	FAIL("database-backed material shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive_wire_diagnostic(dbip, "script.meta", "script", viewport))
	FAIL("database-backed script wire diagnostic coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "script.meta", "script", viewport))
	FAIL("database-backed script shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive(dbip, "pnts.s", 1, 0, viewport))
	FAIL("database-backed PNTS wire primitive coverage should pass");
    if (!exercise_generated_pnts_shaded_points(dbip, "pnts.s", viewport))
	FAIL("database-backed PNTS shaded point realization should pass");
    if (!exercise_generated_primitive(dbip, "annot.s", 1, 0, viewport))
	FAIL("database-backed annotation wire primitive coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "annot.s", "annot", viewport))
	FAIL("database-backed annotation shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive(dbip, "datum.s", 1, 0, viewport))
	FAIL("database-backed datum wire primitive coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "datum.s", "datum", viewport))
	FAIL("database-backed datum shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive(dbip, "joint.s", 1, 0, viewport))
	FAIL("database-backed joint wire primitive coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "joint.s", "joint", viewport))
	FAIL("database-backed joint shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive(dbip, "half.s", 1, 0, viewport))
	FAIL("database-backed halfspace wire primitive coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "half.s", "half", viewport))
	FAIL("database-backed halfspace shaded diagnostic coverage should pass");
    if (!exercise_generated_primitive(dbip, "submodel.s", 1, 0, viewport))
	FAIL("database-backed submodel wire primitive coverage should pass");
    if (!exercise_generated_primitive_mesh_diagnostic(dbip, "submodel.s", "submodel", viewport))
	FAIL("database-backed submodel shaded diagnostic coverage should pass");

    root = new SoSeparator;
    root->ref();

    source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "tet.bot";
    source->sourceRevision = 8;
    source->drawMode = SoBRLDatabaseSource::SHADED;
    if (source->lodBotThreshold.getValue() != 0)
	FAIL("database-backed BoT LoD threshold should default to the full-mesh path");
    source->lodBotThreshold = 8;
    root->addChild(source);

    SoBRLRealizeAction botBelowThresholdRealize;
    botBelowThresholdRealize.apply(root);
    if (botBelowThresholdRealize.getRealizedSourceCount() != 1)
	FAIL("database-backed BoT below LoD threshold should realize one source");
    mesh = source->getRealizedMesh();
    if (!mesh || source->getRealizedMeshCount() != 1)
	FAIL("database-backed BoT below LoD threshold should produce one Obol mesh shape");
    if (mesh->isLodBackedMesh() ||
	    mesh->isOfType(SoBRLLodMeshShape::getClassTypeId()))
	FAIL("database-backed BoT below LoD threshold should keep the plain mesh shape");

    root->unref();

    root = new SoSeparator;
    root->ref();

    source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "tet.bot";
    source->sourceRevision = 8;
    source->drawMode = SoBRLDatabaseSource::SHADED;
    source->lodBotThreshold = 4;
    if (db_mesh_lod_update(dbip, "tet.bot") != BRLCAD_OK)
	FAIL("database-backed BoT LoD cache setup should pass");
    root->addChild(source);

    SoBRLRealizeAction botMeshRealize;
    botMeshRealize.apply(root);
    if (botMeshRealize.getRealizedSourceCount() != 1)
	FAIL("database-backed BoT mesh realization should realize one source");
    if (source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED)
	FAIL("database-backed BoT mesh realization status should be REALIZED");
    mesh = source->getRealizedMesh();
    if (!mesh || source->getRealizedMeshCount() != 1)
	FAIL("database-backed BoT should produce one Obol mesh shape");
    if (mesh->getTriangleCount() != 4)
	FAIL("database-backed BoT mesh should preserve source triangle count");
    if (strcmp(mesh->sourcePath.getValue().getString(), "/tet.bot") != 0)
	FAIL("database-backed BoT mesh should preserve its database path");
    if (!mesh->isLodBackedMesh() ||
	    !mesh->isOfType(SoBRLLodMeshShape::getClassTypeId()))
	FAIL("database-backed BoT at or above LoD threshold should use the LoD-backed mesh shape");
    if (!mesh->lodAvailable.getValue() ||
	    mesh->lodActiveLevel.getValue() < 0 ||
	    mesh->lodFaceCount.getValue() == 0 ||
	    mesh->lodFaceCount.getValue() > 4 ||
	    mesh->lodPointCount.getValue() == 0 ||
	    mesh->lodPointCount.getValue() > 4 ||
	    mesh->lodOriginalPointCount.getValue() == 0 ||
	    mesh->lodOriginalPointCount.getValue() > 4)
	FAIL("database-backed BoT mesh should publish cached RT LoD metadata");
    if (!nearly_equal(mesh->lodBoundsMin.getValue()[0], 0.0f) ||
	    !nearly_equal(mesh->lodBoundsMin.getValue()[1], 0.0f) ||
	    !nearly_equal(mesh->lodBoundsMin.getValue()[2], 0.0f) ||
	    !nearly_equal(mesh->lodBoundsMax.getValue()[0], 2.0f) ||
	    !nearly_equal(mesh->lodBoundsMax.getValue()[1], 2.0f) ||
	    !nearly_equal(mesh->lodBoundsMax.getValue()[2], 2.0f))
	FAIL("database-backed BoT mesh should publish cached RT LoD bounds");
    if (!mesh->lodStagedAvailable.getValue() ||
	    mesh->lodResultKind.getValue() != BRLOBOL_LOD_RESULT_MESH ||
	    mesh->lodQualityTier.getValue() != BRLOBOL_LOD_QUALITY_FAST_DISPLAY ||
	    strcmp(mesh->lodProviderId.getValue().getString(),
		"rt_mesh_lod") != 0 ||
	    mesh->lodCacheKey.getValue().getLength() == 0)
	FAIL("database-backed BoT mesh should publish staged RT LoD result fields");

    SoBRLExportAction botMeshExport;
    botMeshExport.apply(root);
    if (botMeshExport.getTriangleCount() != 4 ||
	    export_triangle_path_count(botMeshExport, "/tet.bot") != 4)
	FAIL("database-backed BoT export should collect current mesh triangles");
    const SoBRLExportAction::TriangleRecord *botTriangle =
	export_triangle_with_path(botMeshExport, "/tet.bot");
    if (!botTriangle ||
	    !botTriangle->lodAvailable ||
	    botTriangle->lodActiveLevel != mesh->lodActiveLevel.getValue() ||
	    botTriangle->lodFaceCount != mesh->lodFaceCount.getValue() ||
	    botTriangle->lodPointCount != mesh->lodPointCount.getValue() ||
	    botTriangle->lodOriginalPointCount != mesh->lodOriginalPointCount.getValue() ||
	    botTriangle->lodNormalCount != mesh->lodNormalCount.getValue() ||
	    botTriangle->lodHasSnappedPoints != mesh->lodHasSnappedPoints.getValue() ||
	    botTriangle->lodHasNormals != mesh->lodHasNormals.getValue())
	FAIL("database-backed BoT export should carry cached RT LoD metadata");
    if (!nearly_equal(botTriangle->lodBoundsMin[0], mesh->lodBoundsMin.getValue()[0]) ||
	    !nearly_equal(botTriangle->lodBoundsMin[1], mesh->lodBoundsMin.getValue()[1]) ||
	    !nearly_equal(botTriangle->lodBoundsMin[2], mesh->lodBoundsMin.getValue()[2]) ||
	    !nearly_equal(botTriangle->lodBoundsMax[0], mesh->lodBoundsMax.getValue()[0]) ||
	    !nearly_equal(botTriangle->lodBoundsMax[1], mesh->lodBoundsMax.getValue()[1]) ||
	    !nearly_equal(botTriangle->lodBoundsMax[2], mesh->lodBoundsMax.getValue()[2]))
	FAIL("database-backed BoT export should carry cached RT LoD bounds");
    bbox = botMeshExport.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 0.0f) ||
	    !nearly_equal(bbox.getMax()[0], 2.0f) ||
	    !nearly_equal(bbox.getMin()[1], 0.0f) ||
	    !nearly_equal(bbox.getMax()[1], 2.0f) ||
	    !nearly_equal(bbox.getMin()[2], 0.0f) ||
	    !nearly_equal(bbox.getMax()[2], 2.0f))
	FAIL("database-backed BoT export should report mesh field bounds");

    root->unref();

    root = new SoSeparator;
    root->ref();

    source = new SoBRLDatabaseSource;
    source->setDatabase(dbip);
    source->path = "assembly.c";
    source->sourceRevision = 9;
    root->addChild(source);

    SoBRLRealizeAction assemblyRealize;
    assemblyRealize.apply(root);
    if (assemblyRealize.getRealizedSourceCount() != 1)
	FAIL("assembly realization should realize one source");
    if (source->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED)
	FAIL("assembly source realization status should be REALIZED");
    if (source->getRealizedShapeCount() != 2)
	FAIL("assembly realization should preserve one Obol shape per database leaf instance");

    SoBRLVListShape *left_shape = shape_with_path(source, "/assembly.c/left.c/box.s");
    SoBRLVListShape *right_shape = shape_with_path(source, "/assembly.c/right.c/box.s");
    if (!left_shape || !right_shape)
	FAIL("assembly leaf shapes should preserve full BRL-CAD instance paths");
    if (!shape_extents_match(left_shape, -2.0f, 3.0f, -3.0f, 4.0f, -4.0f, 5.0f) ||
	    !shape_extents_match(right_shape, -2.0f, 3.0f, -3.0f, 4.0f, -4.0f, 5.0f))
	FAIL("assembly leaf shapes should keep local primitive coordinates under Obol transform nodes");
    if (strcmp(left_shape->sourceName.getValue().getString(), "box.s") != 0 ||
	    strcmp(left_shape->sourceType.getValue().getString(), "arb8") != 0 ||
	    strcmp(right_shape->sourceName.getValue().getString(), "box.s") != 0 ||
	    strcmp(right_shape->sourceType.getValue().getString(), "arb8") != 0)
	FAIL("assembly leaf shapes should preserve primitive name and type identity");
    if (!left_shape->materialColorValid.getValue() ||
	    left_shape->regionId.getValue() != 101 ||
	    left_shape->airCode.getValue() != 7 ||
	    left_shape->materialId.getValue() != 42 ||
	    left_shape->los.getValue() != 9 ||
	    !nearly_equal(left_shape->materialColor.getValue()[0], 230.0f / 255.0f) ||
	    !nearly_equal(left_shape->materialColor.getValue()[1], 26.0f / 255.0f) ||
	    !nearly_equal(left_shape->materialColor.getValue()[2], 51.0f / 255.0f) ||
	    !strstr(left_shape->materialShader.getValue().getString(), "plastic"))
	FAIL("assembly left region should publish inherited material identity on Obol shape fields");
    if (!right_shape->materialColorValid.getValue() ||
	    right_shape->regionId.getValue() != 102 ||
	    right_shape->airCode.getValue() != 8 ||
	    right_shape->materialId.getValue() != 43 ||
	    right_shape->los.getValue() != 10 ||
	    !nearly_equal(right_shape->materialColor.getValue()[0], 26.0f / 255.0f) ||
	    !nearly_equal(right_shape->materialColor.getValue()[1], 77.0f / 255.0f) ||
	    !nearly_equal(right_shape->materialColor.getValue()[2], 204.0f / 255.0f))
	FAIL("assembly right region should publish inherited material identity on Obol shape fields");

    SoGetBoundingBoxAction assemblyBBoxAction(viewport);
    assemblyBBoxAction.apply(root);
    bbox = assemblyBBoxAction.getBoundingBox();
    if (bbox.isEmpty())
	FAIL("assembly source should contribute a transformed bounding box");
    if (!nearly_equal(bbox.getMin()[0], 8.0f) ||
	    !nearly_equal(bbox.getMax()[0], 33.0f) ||
	    !nearly_equal(bbox.getMin()[1], -3.0f) ||
	    !nearly_equal(bbox.getMax()[1], 4.0f) ||
	    !nearly_equal(bbox.getMin()[2], -4.0f) ||
	    !nearly_equal(bbox.getMax()[2], 5.0f))
	FAIL("assembly bounding box should reflect per-instance member matrices");

    SoRayPickAction assemblyPick(viewport);
    assemblyPick.setRay(SbVec3f(30.5f, -3.0f, 10.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    assemblyPick.apply(root);
    pickedPoint = assemblyPick.getPickedPoint();
    if (!pickedPoint)
	FAIL("assembly ray pick should hit transformed leaf geometry");
    rawDetail = pickedPoint->getDetail(right_shape);
    if (!rawDetail || !rawDetail->isOfType(SoBRLPickDetail::getClassTypeId()))
	FAIL("assembly pick should return a BRL-CAD Obol pick detail for the instance");
    pickDetail = static_cast<const SoBRLPickDetail *>(rawDetail);
    if (strcmp(pickDetail->getPath().getString(), "/assembly.c/right.c/box.s") != 0)
	FAIL("assembly pick detail should preserve the full transformed instance path");
    if (strcmp(pickDetail->getSourceName().getString(), "box.s") != 0 ||
	    strcmp(pickDetail->getSourceType().getString(), "arb8") != 0 ||
	    !pickDetail->hasMaterialColor() ||
	    pickDetail->getRegionId() != 102 ||
	    pickDetail->getAirCode() != 8 ||
	    pickDetail->getMaterialId() != 43 ||
	    pickDetail->getLos() != 10 ||
	    !nearly_equal(pickDetail->getMaterialColor()[2], 204.0f / 255.0f))
	FAIL("assembly pick detail should preserve primitive and material identity");

    SoBRLSnapAction assemblySnap;
    assemblySnap.setQueryPoint(SbVec3f(13.0f, 0.0f, -4.0f));
    assemblySnap.setTolerance(0.1f);
    assemblySnap.apply(root);
    if (!assemblySnap.hasCandidate())
	FAIL("assembly snap should find transformed leaf geometry");
    if (strcmp(assemblySnap.getPath().getString(), "/assembly.c/left.c/box.s") != 0)
	FAIL("assembly snap should preserve transformed instance path identity");
    if (!nearly_equal(assemblySnap.getPoint()[0], 13.0f))
	FAIL("assembly snap point should be in transformed model coordinates");

    SoBRLMeasureAction assemblyMeasure;
    assemblyMeasure.setQueryPoint(SbVec3f(30.5f, -3.0f, 0.0f));
    assemblyMeasure.apply(root);
    if (!assemblyMeasure.hasSegments())
	FAIL("assembly measure should find transformed leaf geometry");
    if (assemblyMeasure.getShapeCount() != 2)
	FAIL("assembly measure should preserve one measured shape per transformed leaf");
    bbox = assemblyMeasure.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 8.0f) ||
	    !nearly_equal(bbox.getMax()[0], 33.0f) ||
	    !nearly_equal(bbox.getMin()[1], -3.0f) ||
	    !nearly_equal(bbox.getMax()[1], 4.0f) ||
	    !nearly_equal(bbox.getMin()[2], -4.0f) ||
	    !nearly_equal(bbox.getMax()[2], 5.0f))
	FAIL("assembly measure should report transformed instance bounds");
    if (!assemblyMeasure.hasNearestSegment() ||
	    strcmp(assemblyMeasure.getNearestPath().getString(), "/assembly.c/right.c/box.s") != 0)
	FAIL("assembly measure should report transformed instance path identity");

    SoBRLExportAction assemblyExport;
    assemblyExport.apply(root);
    if (assemblyExport.getLineCount() != total_segment_count(source))
	FAIL("database export should collect realized hierarchy line records");
    if (export_path_count(assemblyExport, "/assembly.c/left.c/box.s") <= 0 ||
	    export_path_count(assemblyExport, "/assembly.c/right.c/box.s") <= 0)
	FAIL("database export should preserve full BRL-CAD hierarchy paths");
    const SoBRLExportAction::LineRecord *leftLine =
	export_line_with_path(assemblyExport, "/assembly.c/left.c/box.s");
    if (!leftLine ||
	    strcmp(leftLine->sourceName.getString(), "box.s") != 0 ||
	    strcmp(leftLine->sourceType.getString(), "arb8") != 0 ||
	    !leftLine->materialColorValid ||
	    leftLine->regionId != 101 ||
	    leftLine->airCode != 7 ||
	    leftLine->materialId != 42 ||
	    leftLine->los != 9 ||
	    !nearly_equal(leftLine->materialColor[0], 230.0f / 255.0f) ||
	    !nearly_equal(leftLine->materialColor[1], 26.0f / 255.0f) ||
	    !nearly_equal(leftLine->materialColor[2], 51.0f / 255.0f) ||
	    !strstr(leftLine->materialShader.getString(), "plastic"))
	FAIL("database export should carry primitive and material identity");
    bbox = assemblyExport.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 8.0f) ||
	    !nearly_equal(bbox.getMax()[0], 33.0f))
	FAIL("database export should apply combination transform state");

    if (!left_shape->visible.getValue() ||
	    !left_shape->selectable.getValue() ||
	    left_shape->colorOverride.getValue() ||
	    left_shape->selected.getValue() ||
	    left_shape->highlighted.getValue() ||
	    left_shape->hiddenLine.getValue() ||
	    left_shape->editEmphasis.getValue() ||
	    left_shape->lodPolicy.getValue() != 0)
	FAIL("realized database shapes should expose default draw-intent fields");

    left_shape->selectedPrimitive.set1Value(0, 0);
    left_shape->highlightedPrimitive.set1Value(0, 0);
    SoBRLMeasureAction selectedSegmentMeasure;
    selectedSegmentMeasure.setSelectionFilter(SoBRLMeasureAction::SELECTED_ONLY);
    selectedSegmentMeasure.apply(root);
    if (!selectedSegmentMeasure.hasSegments() ||
	    selectedSegmentMeasure.getSegmentCount() != 1)
	FAIL("selected-only measure policy should accept selected wire segment geometry");

    SoBRLMeasureAction highlightedSegmentMeasure;
    highlightedSegmentMeasure.setHighlightFilter(SoBRLMeasureAction::HIGHLIGHTED_ONLY);
    highlightedSegmentMeasure.apply(root);
    if (!highlightedSegmentMeasure.hasSegments() ||
	    highlightedSegmentMeasure.getSegmentCount() != 1)
	FAIL("highlighted-only measure policy should accept highlighted wire segment geometry");

    SoBRLExportAction segmentIntentExport;
    segmentIntentExport.apply(root);
    if (segmentIntentExport.getLineCount() < 2 ||
	    !segmentIntentExport.getLine(0).selected ||
	    !segmentIntentExport.getLine(0).highlighted ||
	    segmentIntentExport.getLine(1).selected ||
	    segmentIntentExport.getLine(1).highlighted)
	FAIL("wire export should carry per-segment selected and highlighted draw intent");

    left_shape->colorOverride = TRUE;
    left_shape->color = SbColor(0.9f, 0.1f, 0.2f);
    left_shape->selected = TRUE;
    left_shape->highlighted = TRUE;
    left_shape->ghosted = TRUE;
    left_shape->hiddenLine = TRUE;
    left_shape->editEmphasis = TRUE;
    left_shape->lodPolicy = 5;
    right_shape->visible = FALSE;

    if (!left_shape->selected.getValue() ||
	    !left_shape->highlighted.getValue() ||
	    !left_shape->ghosted.getValue() ||
	    !left_shape->hiddenLine.getValue() ||
	    !left_shape->editEmphasis.getValue() ||
	    left_shape->lodPolicy.getValue() != 5 ||
	    right_shape->visible.getValue())
	FAIL("per-instance draw-intent fields should be inspectable on Obol shapes");

    SoGetBoundingBoxAction intentBBoxAction(viewport);
    intentBBoxAction.apply(root);
    bbox = intentBBoxAction.getBoundingBox();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 8.0f) ||
	    !nearly_equal(bbox.getMax()[0], 13.0f))
	FAIL("visibility draw intent should affect Obol bounding boxes per instance");

    SoRayPickAction hiddenPick(viewport);
    hiddenPick.setRay(SbVec3f(30.5f, -3.0f, 10.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    hiddenPick.apply(root);
    if (hiddenPick.getPickedPoint())
	FAIL("hidden database instance should not be pickable");

    SoBRLSnapAction hiddenSnap;
    hiddenSnap.setQueryPoint(SbVec3f(33.0f, 0.0f, -4.0f));
    hiddenSnap.setTolerance(0.1f);
    hiddenSnap.apply(root);
    if (hiddenSnap.hasCandidate())
	FAIL("hidden database instance should not contribute snap candidates");

    SoBRLMeasureAction intentMeasure;
    intentMeasure.setQueryPoint(SbVec3f(30.5f, -3.0f, 0.0f));
    intentMeasure.apply(root);
    if (!intentMeasure.hasSegments() ||
	    intentMeasure.getShapeCount() != 1 ||
	    !intentMeasure.hasNearestSegment() ||
	    strcmp(intentMeasure.getNearestPath().getString(), "/assembly.c/left.c/box.s") != 0)
	FAIL("visibility draw intent should filter measure traversal per instance");
    bbox = intentMeasure.getBounds();
    if (bbox.isEmpty() ||
	    !nearly_equal(bbox.getMin()[0], 8.0f) ||
	    !nearly_equal(bbox.getMax()[0], 13.0f))
	FAIL("visibility-filtered measure should report only visible instance bounds");

    SoBRLExportAction intentExport;
    intentExport.apply(root);
    if (export_path_count(intentExport, "/assembly.c/right.c/box.s") != 0)
	FAIL("hidden database instance should not be exported");
    if (export_path_count(intentExport, "/assembly.c/left.c/box.s") != left_shape->getSegmentCount())
	FAIL("visible selected instance should remain exported");
    if (intentExport.getLineCount() <= 0 ||
	    !intentExport.getLine(0).selected ||
	    !intentExport.getLine(0).highlighted ||
	    !intentExport.getLine(0).ghosted ||
	    !intentExport.getLine(0).hiddenLine ||
	    !intentExport.getLine(0).editEmphasis ||
	    intentExport.getLine(0).lodPolicy != 5 ||
	    !intentExport.getLine(0).colorOverride ||
	    !nearly_equal(intentExport.getLine(0).color[0], 0.9f) ||
	    !nearly_equal(intentExport.getLine(0).color[1], 0.1f) ||
	    !nearly_equal(intentExport.getLine(0).color[2], 0.2f))
	FAIL("export should carry complete wire draw intent");

    right_shape->visible = TRUE;
    right_shape->selectable = FALSE;

    SoRayPickAction unselectablePick(viewport);
    unselectablePick.setRay(SbVec3f(30.5f, -3.0f, 10.0f), SbVec3f(0.0f, 0.0f, -1.0f));
    unselectablePick.apply(root);
    if (unselectablePick.getPickedPoint())
	FAIL("unselectable visible database instance should not be pickable");

    SoBRLSnapAction unselectableSnap;
    unselectableSnap.setQueryPoint(SbVec3f(33.0f, 0.0f, -4.0f));
    unselectableSnap.setTolerance(0.1f);
    unselectableSnap.apply(root);
    if (unselectableSnap.hasCandidate())
	FAIL("unselectable visible database instance should not contribute snap candidates");

    SoBRLExportAction unselectableExport;
    unselectableExport.apply(root);
    if (export_path_count(unselectableExport, "/assembly.c/right.c/box.s") <= 0)
	FAIL("unselectable visible instance should remain available to export traversal");

    {
	BRLObolViewController rebuiltController(root, NULL);
	if (rebuiltController.getDatabaseSourceCount() != 1)
	    FAIL("rebuilt view controller should discover existing database source nodes");
	SoBRLDatabaseSource *rebuiltSource = rebuiltController.getDatabaseSource(0);
	if (rebuiltSource != source ||
		strcmp(rebuiltSource->path.getValue().getString(), "assembly.c") != 0 ||
		rebuiltSource->realizationStatus.getValue() != SoBRLDatabaseSource::REALIZED ||
		rebuiltSource->needsRealization() ||
		rebuiltSource->getRealizedShapeCount() != 2)
	    FAIL("rebuilt view controller should preserve realized database source state");

	SoBRLVListShape *rebuiltLeft = shape_with_path(rebuiltSource, "/assembly.c/left.c/box.s");
	SoBRLVListShape *rebuiltRight = shape_with_path(rebuiltSource, "/assembly.c/right.c/box.s");
	if (!rebuiltLeft || !rebuiltRight ||
		!rebuiltLeft->selected.getValue() ||
		!rebuiltLeft->highlighted.getValue() ||
		!rebuiltLeft->ghosted.getValue() ||
		!rebuiltLeft->hiddenLine.getValue() ||
		!rebuiltLeft->editEmphasis.getValue() ||
		rebuiltLeft->lodPolicy.getValue() != 5 ||
		!rebuiltLeft->colorOverride.getValue() ||
		!nearly_equal(rebuiltLeft->color.getValue()[0], 0.9f) ||
		!rebuiltRight->visible.getValue() ||
		rebuiltRight->selectable.getValue())
	    FAIL("rebuilt view controller should preserve Obol-held draw-intent fields");

	rebuiltController.clearRenderRequest();
	if (!rebuiltController.realizePending() ||
		rebuiltController.getLastVisitedSourceCount() != 1 ||
		rebuiltController.getLastRealizedSourceCount() != 0 ||
		rebuiltController.getLastFailedSourceCount() != 0 ||
		rebuiltSource->needsRealization() ||
		!rebuiltController.isRenderRequested())
	    FAIL("rebuilt view controller should inspect existing realized scene without forcing refresh");
	if (strcmp(rebuiltController.getRenderReason().getString(), "realize-failed") == 0)
	    FAIL("rebuilt view controller should not report failed realization for current scene");

	if (rebuiltController.replaceDatabaseSource("/assembly.c", dbip,
		    SoBRLDatabaseSource::WIREFRAME, 10) != 1 ||
		rebuiltController.getDatabaseSourceCount() != 1 ||
		rebuiltController.getDatabaseSource(0) != rebuiltSource ||
		rebuiltSource->sourceRevision.getValue() != 10 ||
		!rebuiltSource->needsRealization())
	    FAIL("rebuilt view controller should partial-refresh existing source by path identity");

	if (!rebuiltController.realizePending() ||
		rebuiltController.getLastVisitedSourceCount() != 1 ||
		rebuiltController.getLastRealizedSourceCount() != 1 ||
		rebuiltController.getLastFailedSourceCount() != 0 ||
		rebuiltSource->getRealizedShapeCount() != 2 ||
		!shape_with_path(rebuiltSource, "/assembly.c/left.c/box.s") ||
		!shape_with_path(rebuiltSource, "/assembly.c/right.c/box.s"))
	    FAIL("rebuilt view controller should re-realize partial refresh from existing scene");
    }

    root->unref();
    db_close(dbip);
    bu_file_delete(dbpath);
    bu_dirclear(lodCacheDir);

    if (!exercise_required_hierarchy_model("pinewood.g", "pinewood",
	    2, 41, 21, 501, viewport))
	FAIL("required pinewood hierarchy coverage should pass");

    if (!exercise_required_hierarchy_model("havoc.g", "havoc",
	    10, 100, 10, 100, viewport))
	FAIL("required havoc hierarchy coverage should pass");

    return 0;
}
