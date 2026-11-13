/*          T E S T _ L O D _ S O U R C E _ E V I D E N C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"
#include "bu/app.h"

#include "../lod_source_evidence_private.h"

#include <cstdio>

using Domain = BObolLodAdmissionRevisionDomain;
using Snapshot = BObolLodSourceSnapshotSet;

static bool
expect(bool condition, const char *message)
{
    if (!condition)
	std::fprintf(stderr, "FAIL: %s\n", message);
    return condition;
}

static bool
expect_observation(BObolLodSourceEvidence &evidence, const Snapshot &snapshot,
    std::optional<Domain> expected, const char *message)
{
    return expect(evidence.observe(snapshot) == expected, message);
}

static int
test_pending_observations(void)
{
    BObolLodSourceEvidence evidence;
    BObolLodSourceSnapshot source;
    source.routingId.set(1);
    source.inventoryRevision.set(1);
    source.visibilityRevision = 1;
    source.path = "root.c";
    const Snapshot baseline({source});
    if (!expect_observation(evidence, Snapshot(), std::nullopt,
	    "empty initial observation is inert") ||
	!expect_observation(evidence, baseline, Domain::INVENTORY,
	    "first source publishes inventory") ||
	!expect(evidence.consume(baseline), "first source is consumed") ||
	!expect(&evidence.submitted() == &baseline.sources(),
	    "consumed source evidence shares immutable storage"))
	return 1;

    source.visibilityRevision++;
    const Snapshot hidden({source});
    if (!expect_observation(evidence, hidden, Domain::VISIBILITY,
	    "first pending visibility change publishes") ||
	!expect_observation(evidence, Snapshot({source}), std::nullopt,
	    "duplicate values in different storage are inert") ||
	!expect(evidence.submitted()[0].visibilityRevision == 1,
	    "observation preserves the sparse submission baseline"))
	return 1;

    source.visibilityRevision++;
    const Snapshot restored({source});
    if (!expect_observation(evidence, restored, Domain::VISIBILITY,
	    "second pending visibility change publishes") ||
	!expect(!evidence.consume(hidden),
	    "superseded submission cannot consume newer evidence"))
	return 1;

    source.inventoryRevision.set(2);
    const Snapshot appended({source});
    if (!expect_observation(evidence, appended, Domain::INVENTORY,
	    "inventory following pending visibility publishes inventory"))
	return 1;
    source.visibilityRevision++;
    const Snapshot afterAppend({source});
    if (!expect_observation(evidence, afterAppend, Domain::VISIBILITY,
	    "visibility following pending inventory publishes visibility") ||
	!expect(evidence.consume(afterAppend), "latest observation is consumed") ||
	!expect_observation(evidence, afterAppend, std::nullopt,
	    "consuming published evidence does not publish it twice") ||
	!expect(&evidence.submitted() == &afterAppend.sources(),
	    "observation and submission share the current source vector"))
	return 1;

    if (!expect_observation(evidence, Snapshot(), Domain::INVENTORY,
	    "source removal publishes inventory") ||
	!expect(evidence.consume(Snapshot()), "empty source set is consumed") ||
	!expect_observation(evidence, Snapshot(), std::nullopt,
	    "duplicate source removal is inert") ||
	!expect_observation(evidence, afterAppend, Domain::INVENTORY,
	    "reappearing source publishes inventory"))
	return 1;

    evidence.reset();
    if (!expect(evidence.submitted().empty(), "retirement releases sources") ||
	!expect(!evidence.consume(afterAppend),
	    "retired evidence cannot be consumed") ||
	!expect_observation(evidence, afterAppend, Domain::INVENTORY,
	    "same source after retirement is a new observation"))
	return 1;
    return 0;
}

static int
test_identity_domains(void)
{
    enum IdentityField {
	ROUTE, DATABASE, DATABASE_ID, PATH, DRAW_MODE, REPRESENTATION,
	VISIBILITY, BOT_THRESHOLD, SOURCE_REVISION, INPUTS_REVISION, FIELD_COUNT
    };
    BObolLodSourceSnapshot initial;
    initial.routingId.set(1);
    initial.inventoryRevision.set(1);
    initial.visibilityRevision = 1;
    initial.path = "root.c";
    const Snapshot baseline({initial});
    int databaseIdentity = 0;
    for (int field = ROUTE; field < FIELD_COUNT; ++field) {
	BObolLodSourceEvidence evidence;
	(void)evidence.observe(baseline);
	BObolLodSourceSnapshot changed = initial;
	switch (field) {
	    case ROUTE: changed.routingId.set(2); break;
	    case DATABASE:
		/* Identity comparison never dereferences the database. */
		changed.database = reinterpret_cast<struct db_i *>(&databaseIdentity);
		break;
	    case DATABASE_ID: changed.databaseId = "replacement.g"; break;
	    case PATH: changed.path = "replacement.c"; break;
	    case DRAW_MODE: changed.drawMode++; break;
	    case REPRESENTATION: changed.representationMode++; break;
	    case VISIBILITY: changed.visible = TRUE; break;
	    case BOT_THRESHOLD: changed.lodBotThreshold++; break;
	    case SOURCE_REVISION: changed.sourceRevision++; break;
	    case INPUTS_REVISION: changed.inputsRevision++; break;
	}
	const Snapshot replacement({changed});
	if (!expect_observation(evidence, replacement, Domain::INVENTORY,
		"a changed source identity publishes inventory while pending") ||
	    !expect(!evidence.consume(baseline),
		"a replacement rejects predecessor consumption") ||
	    !expect(evidence.consume(replacement),
		"replacement evidence is consumable"))
	    return 1;
    }

    BObolLodSourceEvidence evidence;
    (void)evidence.observe(baseline);
    (void)evidence.consume(baseline);
    BObolLodSourceSnapshot second = initial;
    second.routingId.set(2);
    second.path = "second.c";
    const Snapshot twoSources({initial, second});
    if (!expect_observation(evidence, twoSources, Domain::INVENTORY,
	    "adding a source publishes inventory") ||
	!expect(evidence.consume(twoSources), "both sources consumed") ||
	!expect_observation(evidence, Snapshot({second}), Domain::INVENTORY,
	    "removing one source publishes inventory"))
	return 1;
    return 0;
}

int
main(void)
{
    bu_setprogname("test_lod_source_evidence");
    return test_pending_observations() || test_identity_domains();
}
