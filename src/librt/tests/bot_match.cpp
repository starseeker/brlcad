/*                    B O T _ M A T C H . C P P
 * BRL-CAD
 *
 * Copyright (c) 2013-2026 United States Government as represented by
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

#include "common.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "vmath.h"
#include "bu/app.h"
#include "bu/log.h"
#include "bu/time.h"
#include "bu/units.h"
#include "bg/pca.h"
#include "bg/trimesh.h"
#include "raytrace.h"

struct bot_signature {
    struct bg_trimesh_pca_signature pca;
    size_t vertex_count;
    size_t face_count;
};

struct hash_results {
    std::unordered_map<struct directory *, struct bot_signature> signatures;
    std::unordered_map<unsigned long long, std::unordered_set<struct directory *>> bot_groups_hashed;
    size_t rejected = 0;
};

// Test hash method
void
test_hash(struct hash_results *h, struct db_i *dbip, std::vector<struct directory *> &bots)
{
    if (!h || !dbip)
	return;

    int64_t start = bu_gettime();
    int64_t msgtime = bu_gettime();

    // Iterate over all BoTs
    for(size_t i = 0; i < bots.size(); i++) {
        struct directory *wdp = bots[i];

	if (bu_gettime() - msgtime > 5000000.0) {
	    bu_log("Processing %s (%zd of %zd)\n", wdp->d_namep, i, bots.size());
	    msgtime = bu_gettime();
	}

	struct rt_db_internal intern = RT_DB_INTERNAL_INIT_ZERO;
	if (rt_db_get_internal(&intern, wdp, dbip, NULL) < 0)
	    continue;
	struct rt_bot_internal *orig_bot = (struct rt_bot_internal *)(intern.idb_ptr);
	struct bot_signature signature;
	signature.vertex_count = orig_bot->num_vertices;
	signature.face_count = orig_bot->num_faces;
	if (bg_trimesh_pca_get_signature(&signature.pca, orig_bot->faces,
		orig_bot->num_faces, (const point_t *)orig_bot->vertices,
		orig_bot->num_vertices, VUNITIZE_TOL, 1.0e-6) != BRLCAD_OK) {
	    h->rejected++;
	    rt_db_free_internal(&intern);
	    continue;
	}
	rt_db_free_internal(&intern);

	// The signature is a candidate key only; test_diff validates every match.
	h->signatures[wdp] = signature;
	h->bot_groups_hashed[signature.pca.hash].insert(wdp);
    }

    std::unordered_map<unsigned long long, std::unordered_set<struct directory *>>::iterator b_it;
    size_t gcnt = 0;
    size_t gobjs = 0;
    for (b_it = h->bot_groups_hashed.begin(); b_it != h->bot_groups_hashed.end(); ++b_it) {
	// Anything by itself in a set didn't form any groups
	if (b_it->second.size() == 1)
	    continue;

	// OK, found a non-trivial group
	gcnt++;
	gobjs += b_it->second.size();
    }

    int64_t elapsed = bu_gettime() - start;
    fastf_t seconds = elapsed / 1000000.0;

    bu_log("PCA signature complete (%f sec) - %zd candidate groups found, %zd of %zd BoTs are part of a group (%zd PCA-ambiguous BoTs skipped).\n", seconds, gcnt, gobjs, bots.size(), h->rejected);
}

struct diff_results {
    std::unordered_map<struct directory *, std::unordered_set<struct directory *>> bot_groups;
};

// Test diff method
void
test_diff(struct diff_results *d, struct db_i *dbip, const struct hash_results *h)
{

    if (!d || !dbip || !h)
	return;

    int64_t start = bu_gettime();
    size_t comparisons = 0;
    for (const auto &bucket : h->bot_groups_hashed) {
	if (bucket.second.size() < 2)
	    continue;
	std::vector<struct directory *> candidates(bucket.second.begin(), bucket.second.end());
	std::vector<bool> matched(candidates.size(), false);
	for (size_t i = 0; i < candidates.size(); i++) {
	    if (matched[i])
		continue;
	    struct directory *wdp = candidates[i];
	    const struct bot_signature &firstSignature = h->signatures.at(wdp);
	    struct rt_db_internal firstInternal = RT_DB_INTERNAL_INIT_ZERO;
	    if (rt_db_get_internal(&firstInternal, wdp, dbip, NULL) < 0)
		continue;
	    const struct rt_bot_internal *firstBot =
		(const struct rt_bot_internal *)firstInternal.idb_ptr;
	    for (size_t j = i + 1; j < candidates.size(); j++) {
		if (matched[j])
		    continue;
		struct directory *cdp = candidates[j];
		const struct bot_signature &secondSignature = h->signatures.at(cdp);
		if (firstSignature.vertex_count != secondSignature.vertex_count ||
		    firstSignature.face_count != secondSignature.face_count)
		    continue;
		struct rt_db_internal secondInternal = RT_DB_INTERNAL_INIT_ZERO;
		if (rt_db_get_internal(&secondInternal, cdp, dbip, NULL) < 0)
		    continue;
		const struct rt_bot_internal *secondBot =
		    (const struct rt_bot_internal *)secondInternal.idb_ptr;
		comparisons++;
		if (bg_trimesh_pca_equal(&firstSignature.pca, firstBot->faces,
			firstBot->num_faces, (const point_t *)firstBot->vertices,
			firstBot->num_vertices, &secondSignature.pca, secondBot->faces,
			secondBot->num_faces, (const point_t *)secondBot->vertices,
			secondBot->num_vertices, VUNITIZE_TOL) == 0) {
		    d->bot_groups[wdp].insert(cdp);
		    matched[j] = true;
		}
		rt_db_free_internal(&secondInternal);
	    }
	    rt_db_free_internal(&firstInternal);
	}
    }

    int64_t elapsed = bu_gettime() - start;
    fastf_t seconds = elapsed / 1000000.0;

    // Get a count of all objects that would up part of some group
    std::unordered_map<struct directory *, std::unordered_set<struct directory *>>::iterator bg_it;
    size_t gobjs = 0;
    for (bg_it = d->bot_groups.begin(); bg_it != d->bot_groups.end(); ++bg_it) {
	gobjs++;  // Account for the key dp
	gobjs += bg_it->second.size(); // Add the matches
    }

    bu_log("Candidate validation complete (%f sec, %zd comparisons) - %zd match groups found, %zd BoTs are part of a group.\n",seconds, comparisons, d->bot_groups.size(), gobjs);
}


int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);

    if (argc != 2) {
	bu_exit(1, "Usage: %s file.g", argv[0]);
    }

    int64_t start = bu_gettime();

    struct db_i *dbip = db_open(argv[1], DB_OPEN_READONLY);
    if (dbip == DBI_NULL)
	bu_exit(1, "ERROR: Unable to read from %s\n", argv[1]);

    if (db_dirbuild(dbip) < 0)
	bu_exit(1, "ERROR: Unable to read from %s\n", argv[1]);

    db_update_nref(dbip);

    int64_t elapsed = bu_gettime() - start;
    fastf_t seconds = elapsed / 1000000.0;

    bu_log("Setup time: %g seconds\n", seconds);

    // Find all BoT objects in the .g database
    std::vector<struct directory *> bots;
    start = bu_gettime();
#if 0
    // NOTE - several seconds slower than just directly
    // iterating over the hashes on a large model...
    struct bu_ptbl bot_objs = BU_PTBL_INIT_ZERO;
    const char *bot_search = "-type bot";
    (void)db_search(&bot_objs, DB_SEARCH_RETURN_UNIQ_DP|DB_SEARCH_FLAT, bot_search, 0, NULL, dbip, NULL, NULL, NULL);
    for (size_t i = 0; i < BU_PTBL_LEN(&bot_objs); i++)
	bots.push_back((struct directory *)BU_PTBL_GET(&bot_objs, i));
    db_search_free(&bot_objs);
#else
    for (int i = 0; i < RT_DBNHASH; i++) {
	struct directory *dp;
	for (dp = db_dirptr(dbip, i); dp != RT_DIR_NULL; dp = dp->d_forw) {
	    if (dp->d_addr == RT_DIR_PHONY_ADDR)
		continue;
	    if (dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
		continue;
	    bots.push_back(dp);
	}
    }
#endif
    elapsed = bu_gettime() - start;
    seconds = elapsed / 1000000.0;

    bu_log("Initial BoT search time: %g seconds\n", seconds);

    struct hash_results h;
    test_hash(&h, dbip, bots);

    struct diff_results d;
    test_diff(&d, dbip, &h);

#if 0
    // Print any groups found
    for(bg_it = bot_groups.begin(); bg_it != bot_groups.end(); ++bg_it) {
	if (!bg_it->second.size())
	    continue;
	bu_log("%s:\n", bg_it->first->d_namep);
	for (d_it = bg_it->second.begin(); d_it != bg_it->second.end(); ++d_it) {
	    bu_log("\t%s\n", (*d_it)->d_namep);
	}
    }
#endif

    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
