/*           S C E N E _ N O D E _ H A N D L E S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file scene_node_handles.cpp
 *
 * Focused contract tests for GED-owned scene-node value handles.
 */

#include "common.h"

#include "BObol/BInit.h"
#include "bu/app.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "ged.h"
#include "ged/commands.h"
#include "ged/scene.h"
#include "wdb.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define FAIL(_msg) \
    do { \
	fprintf(stderr, "FAIL: %s\n", (_msg)); \
	return 1; \
    } while (0)

static int
make_handle_db(const char *path)
{
    struct rt_wdb *wdbp = wdb_fopen(path);
    if (!wdbp)
	return 0;

    point_t first = {0.0, 0.0, 0.0};
    point_t second = {4.0, 0.0, 0.0};
    int ret = mk_id(wdbp, "GED scene-node handle test") == 0 &&
	mk_sph(wdbp, "first.s", first, 1.0) == 0 &&
	mk_sph(wdbp, "second.s", second, 1.0) == 0;
    wdb_close(wdbp);
    return ret;
}

static int
node_is_live(struct ged *gedp, ged_scene_node_ref node)
{
    struct ged_draw_scene_tree_summary summary;
    memset(&summary, 0, sizeof(summary));
    return !ged_scene_node_ref_is_null(node) &&
	ged_scene_node_name(gedp, node) &&
	ged_scene_node_tree_summary(gedp, node, &summary) &&
	summary.valid;
}

static int
node_is_rejected(struct ged *gedp, ged_scene_node_ref node)
{
    struct ged_draw_scene_tree_summary summary;
    memset(&summary, 0, sizeof(summary));
    return !ged_scene_node_name(gedp, node) &&
	!ged_scene_node_tree_summary(gedp, node, &summary) &&
	!summary.valid;
}

int
main(int argc, char **argv)
{
    (void)argc;
    bu_setprogname(argv[0]);
    bobol_init(NULL);

    const char *first_db = "ged_scene_node_handles_a.g";
    const char *second_db = "ged_scene_node_handles_b.g";
    bu_file_delete(first_db);
    bu_file_delete(second_db);
    if (!make_handle_db(first_db) || !make_handle_db(second_db))
	FAIL("could not create scene-node handle databases");

    struct ged *gedp = ged_open("db", first_db, 1);
    struct ged *other = ged_open("db", second_db, 1);
    if (!gedp || !other)
	FAIL("could not open scene-node handle databases");

    const char *draw_first[] = {"draw", "first.s"};
    if (ged_exec_draw(gedp, 2, draw_first) != BRLCAD_OK)
	FAIL("initial draw failed");

    ged_scene_node_ref first = ged_scene_first_node(gedp);
    if (!node_is_live(gedp, first) ||
	ged_scene_node_ref_is_null(first) ||
	!ged_scene_node_ref_equal(first, first))
	FAIL("a freshly issued scene-node handle must resolve");

    if (!node_is_rejected(other, first))
	FAIL("a scene-node handle must not resolve through another GED owner");

    ged_scene_node_ref tampered = first;
    tampered.owner = UINT64_MAX;
    if (!node_is_rejected(gedp, tampered))
	FAIL("a scene-node handle with the wrong owner must be rejected");
    tampered = first;
    tampered.id = UINT64_MAX;
    if (!node_is_rejected(gedp, tampered))
	FAIL("a scene-node handle with the wrong id must be rejected");
    tampered = first;
    tampered.generation++;
    if (!node_is_rejected(gedp, tampered))
	FAIL("a scene-node handle with the wrong generation must be rejected");
    tampered = first;
    tampered.scene_revision++;
    if (!node_is_rejected(gedp, tampered))
	FAIL("a scene-node handle with the wrong revision must be rejected");

    const uint64_t revision_before = ged_draw_scene_revision(gedp);
    const char *draw_second[] = {"draw", "second.s"};
    if (ged_exec_draw(gedp, 2, draw_second) != BRLCAD_OK ||
	ged_draw_scene_revision(gedp) <= revision_before)
	FAIL("a structural draw must advance the scene revision");
    if (!node_is_rejected(gedp, first))
	FAIL("a scene-node handle must be stale after structural mutation");

    ged_scene_node_ref current = ged_scene_first_node(gedp);
    if (!node_is_live(gedp, current) ||
	ged_scene_node_ref_equal(first, current) ||
	current.scene_revision != ged_draw_scene_revision(gedp))
	FAIL("reacquiring after mutation must issue a current handle");

    const char *current_name = ged_scene_node_name(gedp, current);
    char *erase_name = current_name ? bu_strdup(current_name) : NULL;
    if (!erase_name)
	FAIL("current scene-node handle must expose an eraseable name");
    const char *erase_current[] = {"erase", erase_name};
    int erase_ret = ged_exec_erase(gedp, 2, erase_current);
    bu_free(erase_name, "scene-node erase name");
    if (erase_ret != BRLCAD_OK)
	FAIL("erase of referenced scene node failed");
    if (!node_is_rejected(gedp, current))
	FAIL("a scene-node handle must be stale after object removal");

    ged_scene_node_ref surviving = ged_scene_first_node(gedp);
    if (!node_is_live(gedp, surviving))
	FAIL("the surviving object must be reacquirable after removal");

    ged_scene_node_ref closed_owner_ref = surviving;
    ged_close(gedp);
    gedp = NULL;

    if (!node_is_rejected(other, closed_owner_ref))
	FAIL("GED teardown must leave handles rejected by a live replacement owner");

    ged_close(other);
    bu_file_delete(first_db);
    bu_file_delete(second_db);
    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
