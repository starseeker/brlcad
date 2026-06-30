/*                    T E S T _ V I E W _ S T O R E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "brlobol.h"
#include "bu/app.h"
#include "bu/file.h"
#include "raytrace.h"
#include "wdb.h"

#include <Inventor/SoType.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoText2.h>

#include <stdio.h>
#include <vector>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static int
count_nodes_of_type(SoNode *node, SoType type)
{
    if (!node)
	return 0;

    int count = node->isOfType(type) ? 1 : 0;
    if (node->isOfType(SoGroup::getClassTypeId())) {
	SoGroup *group = static_cast<SoGroup *>(node);
	for (int i = 0; i < group->getNumChildren(); i++)
	    count += count_nodes_of_type(group->getChild(i), type);
    }
    return count;
}

static int
test_feature_nodes(BRLObolViewController &view)
{
    std::vector<BRLObolLabel> labels;
    BRLObolLabel label;
    label.text = "store label";
    label.point = SbVec3f(1.0f, 2.0f, 0.0f);
    label.hasLeader = TRUE;
    label.target = SbVec3f(0.0f, 0.0f, 0.0f);
    labels.push_back(label);

    BRLObolFeatureHandle labelHandle = view.features().publishLabels(
	    "labels",
	    BRLObolFeatureScope::Shared,
	    labels);
    if (!labelHandle.isValid())
	FAIL("publishLabels should return a valid handle");
    if (count_nodes_of_type(view.features().node(labelHandle),
		SoText2::getClassTypeId()) != 1)
	FAIL("label feature should realize a SoText2 node");

    std::vector<SbVec3f> points;
    points.push_back(SbVec3f(0.0f, 0.0f, 0.0f));
    points.push_back(SbVec3f(1.0f, 0.0f, 0.0f));
    points.push_back(SbVec3f(1.0f, 1.0f, 0.0f));
    points.push_back(SbVec3f(0.0f, 1.0f, 0.0f));
    std::vector<SbVec3f> normals;
    std::vector<int32_t> indices;
    indices.push_back(0);
    indices.push_back(1);
    indices.push_back(2);
    indices.push_back(0);
    indices.push_back(2);
    indices.push_back(3);

    BRLObolFeatureHandle meshHandle = view.features().publishIndexedFaceSet(
	    "surface",
	    BRLObolFeatureScope::Shared,
	    points,
	    normals,
	    indices);
    if (!meshHandle.isValid())
	FAIL("publishIndexedFaceSet should return a valid handle");

    SoNode *node = view.features().node(meshHandle);
    if (!node || !node->isOfType(SoBRLMeshShape::getClassTypeId()))
	FAIL("indexed face feature should realize a SoBRLMeshShape");
    SoBRLMeshShape *mesh = static_cast<SoBRLMeshShape *>(node);
    if (mesh->getTriangleCount() != 2)
	FAIL("indexed face feature should triangulate to two triangles");

    return 0;
}

static int
test_polygon_nodes_and_sketch(BRLObolViewController &view)
{
    plane_t viewPlane;
    HSET(viewPlane, 0.0, 0.0, 1.0, 0.0);

    BRLObolPolygonHandle handle = view.polygons().create(
	    "poly",
	    BRLObolFeatureScope::Shared,
	    BRLObolPolygonType::Square,
	    SbVec3f(0.0f, 0.0f, 0.0f),
	    viewPlane,
	    0.0f);
    if (!handle.isValid())
	FAIL("polygon create should return a valid handle");
    if (!view.polygons().updateScreenPoint(handle, 10, 8,
		BRLObolPolygonUpdate::Default))
	FAIL("square polygon update should succeed");
    if (!view.polygons().setFill(handle, TRUE, SbVec2f(1.0f, 0.0f), 2.0f))
	FAIL("polygon setFill should succeed");
    if (!view.polygons().setFillColor(handle, SbColor(0.2f, 0.4f, 0.8f)))
	FAIL("polygon setFillColor should succeed");
    if (!view.polygons().setEdgeColor(handle, SbColor(1.0f, 0.0f, 0.0f)))
	FAIL("polygon setEdgeColor should succeed");

    SoNode *node = view.polygons().node(handle);
    if (!node || !node->isOfType(SoGroup::getClassTypeId()))
	FAIL("filled polygon should realize as a group");
    if (count_nodes_of_type(node, SoBRLMeshShape::getClassTypeId()) != 1)
	FAIL("filled polygon should have one mesh fill child");
    if (count_nodes_of_type(node, SoBRLVListShape::getClassTypeId()) != 1)
	FAIL("filled polygon should have one outline child");

    char dbpath[MAXPATHLEN];
    bu_dir(dbpath, MAXPATHLEN, BU_DIR_CURR,
	    "brlobol_view_store_test.g", NULL);
    bu_file_delete(dbpath);

    struct rt_wdb *wdbp = wdb_fopen_v(dbpath, 5);
    if (!wdbp)
	FAIL("wdb_fopen_v should succeed");
    struct db_i *dbip = wdbp->dbip;

    if (!view.polygons().exportSketch(handle, dbip, "poly.s")) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("polygon exportSketch should succeed");
    }

    struct directory *dp = db_lookup(dbip, "poly.s", LOOKUP_QUIET);
    if (dp == RT_DIR_NULL) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("exported sketch should be in the database directory");
    }

    BRLObolPolygonStore imported;
    BRLObolPolygonHandle importedHandle = imported.importSketch(
	    "imported",
	    BRLObolFeatureScope::Shared,
	    dbip,
	    dp);
    if (!importedHandle.isValid()) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("polygon importSketch should succeed");
    }

    BRLObolPolygonRecord record;
    if (!imported.record(importedHandle, record)) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("imported polygon record lookup should succeed");
    }
    if (record.pointCount != 4 || record.contourCount != 1) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("imported polygon should preserve square contour geometry");
    }
    if (!record.fill || record.fillSpacing != 2.0f ||
	    record.type != BRLObolPolygonType::Square) {
	wdb_close(wdbp);
	bu_file_delete(dbpath);
	FAIL("imported polygon should preserve visual/type attributes");
    }

    wdb_close(wdbp);
    bu_file_delete(dbpath);
    return 0;
}

int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);
    if (argc != 1)
	FAIL("unexpected arguments");

    brlobol_init(NULL);

    SoSeparator *root = new SoSeparator;
    root->ref();
    {
	BRLObolViewController view(root);
	if (test_feature_nodes(view) != 0) {
	    root->unref();
	    return 1;
	}
	if (test_polygon_nodes_and_sketch(view) != 0) {
	    root->unref();
	    return 1;
	}
    }
    root->unref();

    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
