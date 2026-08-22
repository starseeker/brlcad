/*                 D B _ D I S C O V E R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "raytrace.h"
#include "rt/db_discovery.h"


struct discovery_capture {
    std::vector<std::pair<std::string, std::string>> arcs;
};


static void
capture_arc(struct db_i *UNUSED(dbip), struct directory *parent,
	struct directory *UNUSED(child), const char *child_name,
	db_op_t UNUSED(operation), matp_t UNUSED(matrix), void *data)
{
    discovery_capture *capture = static_cast<discovery_capture *>(data);
    if (!capture || !parent || !parent->d_namep || !child_name)
	return;
    capture->arcs.emplace_back(parent->d_namep, child_name);
}


static bool
capture_counts(struct db_i *dbip, std::map<std::string, size_t> &counts)
{
    if (!dbip)
	return false;
    struct directory *dp = RT_DIR_NULL;
    FOR_ALL_DIRECTORY_START(dp, dbip) {
	if (dp && dp->d_namep)
	    counts[dp->d_namep] = dp->d_nref;
    } FOR_ALL_DIRECTORY_END;
    return true;
}


int
main(int argc, char **argv)
{
    if (argc != 2)
	return 2;

    discovery_capture legacy_arcs;
    std::map<std::string, size_t> legacy_counts;
    struct db_i *legacy = db_open(argv[1], DB_OPEN_READONLY);
    if (!legacy)
	return 3;
    if (db_add_update_nref_clbk(legacy, capture_arc, &legacy_arcs) < 0 ||
	db_dirbuild(legacy) < 0 ||
	!capture_counts(legacy, legacy_counts)) {
	db_close(legacy);
	return 4;
    }
    db_close(legacy);

    discovery_capture parallel_arcs;
    std::map<std::string, size_t> parallel_counts;
    struct db_i *parallel = db_open(argv[1], DB_OPEN_READONLY);
    if (!parallel)
	return 5;
    if (db_add_update_nref_clbk(parallel, capture_arc, &parallel_arcs) < 0) {
	db_close(parallel);
	return 6;
    }
    struct rt_db_discovery_options options;
    rt_db_discovery_options_init(&options);
    options.max_workers = 4;
    struct rt_db_hierarchy *hierarchy = nullptr;
    if (rt_db_discovery_build(parallel, &options, &hierarchy, nullptr) < 0 ||
	!hierarchy || !capture_counts(parallel, parallel_counts)) {
	rt_db_hierarchy_destroy(hierarchy);
	db_close(parallel);
	return 7;
    }
    const size_t hierarchy_arcs = rt_db_hierarchy_arc_count(hierarchy);
    rt_db_hierarchy_destroy(hierarchy);
    db_close(parallel);

    std::sort(legacy_arcs.arcs.begin(), legacy_arcs.arcs.end());
    std::sort(parallel_arcs.arcs.begin(), parallel_arcs.arcs.end());
    if (legacy_counts != parallel_counts ||
	legacy_arcs.arcs != parallel_arcs.arcs ||
	hierarchy_arcs != parallel_arcs.arcs.size())
	return 8;
    return 0;
}

