/* T E S T _ S C E N E _ R E P R E S E N T A T I O N _ L I F E C Y C L E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BInit.h"
#include "BObol/BSceneController.h"
#include "bu/app.h"
#include "bu/file.h"
#include "raytrace.h"
#include "rt/wdb.h"
#include "wdb.h"

#include <Inventor/nodes/SoGroup.h>

#include <cstdio>
#include <cstring>
#include <set>
#include <string>
#include <vector>

#define FAIL(_msg) \
    do { \
	std::fprintf(stderr, "FAIL: %s\n", _msg); \
	return 1; \
    } while (0)

static void
collect_sources(SoNode *node, std::vector<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;
    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	sources.push_back(static_cast<SoBRLDatabaseSource *>(node));
    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;
    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	collect_sources(group->getChild(i), sources);
}

static int
verify_scene_index(BObolSceneController &scene, SoNode *root,
	const char *const *keys, size_t key_count)
{
    std::vector<SoBRLDatabaseSource *> tree_sources;
    collect_sources(root, tree_sources);
    if (tree_sources.size() != key_count ||
	scene.getDatabaseSourceCount() != static_cast<int>(key_count))
	return 0;

    std::set<SoBRLDatabaseSource *> tree_set(tree_sources.begin(),
	tree_sources.end());
    std::set<SoBRLDatabaseSource *> index_set;
    for (int i = 0; i < scene.getDatabaseSourceCount(); i++) {
	SoBRLDatabaseSource *source = scene.getDatabaseSource(i);
	if (!source || tree_set.find(source) == tree_set.end() ||
	    !index_set.insert(source).second)
	    return 0;
    }
    if (index_set != tree_set)
	return 0;

    for (size_t i = 0; i < key_count; i++) {
	SoBRLDatabaseSource *indexed = scene.findDatabaseSourceInstance(keys[i]);
	BObolDatabaseSourceSummary summary;
	if (!indexed || tree_set.find(indexed) == tree_set.end() ||
	    !indexed->getSummary(summary) || !summary.valid ||
	    !BU_STR_EQUAL(summary.instanceKey.getString(), keys[i]))
	    return 0;
    }
    return 1;
}

static int
make_database(const char *path)
{
    struct rt_wdb *wdbp = wdb_fopen_v(path, 5);
    if (!wdbp)
	return 0;
    point_t min = {-1.0, -1.0, -1.0};
    point_t max = {1.0, 1.0, 1.0};
    const int ret = mk_rpp(wdbp, "box.s", min, max) == 0 ? 1 : 0;
    wdb_close(wdbp);
    return ret;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    if (argc != 1)
	FAIL("unexpected arguments");
    bobol_init(NULL);

    char dbpath[MAXPATHLEN] = {0};
    FILE *fp = bu_temp_file(dbpath, sizeof(dbpath));
    if (!fp)
	FAIL("could not create temporary database path");
    std::fclose(fp);
    if (!make_database(dbpath)) {
	(void)bu_file_delete(dbpath);
	FAIL("could not create database fixture");
    }

    struct db_i *dbip = db_open(dbpath, DB_OPEN_READONLY);
    if (!dbip || db_dirbuild(dbip) < 0) {
	if (dbip)
	    db_close(dbip);
	(void)bu_file_delete(dbpath);
	FAIL("could not open database fixture");
    }

    SoGroup *root = new SoGroup;
    root->ref();
    BObolSceneController scene(root);
    const char *wire_key = "box.s::wire";
    const char *eval_key = "box.s::eval-wire";
    const char *both_keys[] = {wire_key, eval_key};
    const char *wire_only[] = {wire_key};
    int ret = 0;

    if (scene.replaceDatabaseSourceInstanceRepresentation(wire_key, "box.s",
	    wire_key, SoBRLDatabaseSource::REPRESENTATION_WIRE, dbip,
	    SoBRLDatabaseSource::WIREFRAME, 1) < 0 ||
	scene.replaceDatabaseSourceInstanceRepresentation(eval_key, "box.s",
	    eval_key, SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE, dbip,
	    SoBRLDatabaseSource::WIREFRAME, 1) < 0 ||
	!verify_scene_index(scene, root, both_keys, 2)) {
	std::fprintf(stderr, "FAIL: two representations must have matching scene and index membership\n");
	ret = 1;
	goto cleanup;
    }

    /* The action must traverse the same evaluated source returned by the
     * index.  A valid simple primitive makes this a current-result check,
     * rather than merely an enumeration test. */
    if (!scene.realizePending()) {
	std::fprintf(stderr, "FAIL: simple representation fixture should realize\n");
	ret = 1;
	goto cleanup;
    }
    {
	SoBRLDatabaseSource *evaluated = scene.findDatabaseSourceInstance(eval_key);
	BObolDatabaseSourceSummary summary;
	if (!evaluated || !evaluated->getSummary(summary) || !summary.valid ||
	    summary.representationMode !=
	SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE || summary.stale ||
	    summary.realizationStatus != SoBRLDatabaseSource::REALIZED) {
	    std::fprintf(stderr, "FAIL: evaluated representation should realize as current\n");
	    ret = 1;
	    goto cleanup;
	}
    }

    if (scene.removeDatabaseSourceInstance(eval_key) != 1 ||
	!verify_scene_index(scene, root, wire_only, 1) ||
	scene.findDatabaseSourceInstance(eval_key)) {
	std::fprintf(stderr, "FAIL: erase must remove the exact representation from scene and index\n");
	ret = 1;
	goto cleanup;
    }

    if (scene.replaceDatabaseSourceInstanceRepresentation(eval_key, "box.s",
	    eval_key, SoBRLDatabaseSource::REPRESENTATION_EVAL_WIRE, dbip,
	    SoBRLDatabaseSource::WIREFRAME, 2) < 0 ||
	!verify_scene_index(scene, root, both_keys, 2) || !scene.realizePending()) {
	std::fprintf(stderr, "FAIL: redraw must restore exact representation and index membership\n");
	ret = 1;
    }

cleanup:
    root->unref();
    db_close(dbip);
    (void)bu_file_delete(dbpath);
    return ret;
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
