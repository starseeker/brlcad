/*      T E S T _ C O M P A C T _ O C C U R R E N C E _ S T R E A M . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "bu/app.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <thread>
#include <vector>

static BObolCompactOccurrence
occurrence(size_t index)
{
    BObolCompactOccurrence result;
    result.occurrenceIndex = static_cast<uint32_t>(index);
    return result;
}

static int
test_priority_and_state(void)
{
    BObolCompactOccurrenceStream stream;
    stream.setExpectedCount(100000);
    size_t preparationCompleted = 1;
    size_t preparationTotal = 1;
    stream.getPreparationWorkCount(preparationCompleted, preparationTotal);
    if (preparationCompleted != 0 || preparationTotal != 0) {
	std::fprintf(stderr, "FAIL: initial preparation work rank\n");
	return 1;
    }
    stream.setPreparationWorkCount(3);
    stream.notePreparationWorkCompleted();
    stream.notePreparationWorkCompleted();
    stream.notePreparationWorkCompleted();
    stream.notePreparationWorkCompleted();
    stream.getPreparationWorkCount(preparationCompleted, preparationTotal);
    if (preparationCompleted != 3 || preparationTotal != 3) {
	std::fprintf(stderr, "FAIL: bounded preparation work rank\n");
	return 1;
    }
    stream.setPreparationWorkCount(100000);
    stream.completePreparationWork();
    stream.getPreparationWorkCount(preparationCompleted, preparationTotal);
    if (preparationCompleted != 100000 || preparationTotal != 100000) {
	std::fprintf(stderr, "FAIL: post-walk preparation closure\n");
	return 1;
    }
    BObolCompactSourceProfile profile;
    profile.valid = TRUE;
    profile.occurrenceCount = 100000;
    profile.uniqueAssetCount = 25000;
    profile.encodedSourceBytes = 64ULL * 1024ULL * 1024ULL;
    profile.largestAssetBytes = 4ULL * 1024ULL * 1024ULL;
    profile.reusedOccurrenceCount = 75000;
    stream.setSourceProfile(profile);
    BObolCompactSourceProfile observedProfile;
    if (!stream.getSourceProfile(observedProfile) ||
	observedProfile.occurrenceCount != profile.occurrenceCount ||
	observedProfile.uniqueAssetCount != profile.uniqueAssetCount ||
	observedProfile.encodedSourceBytes != profile.encodedSourceBytes ||
	observedProfile.largestAssetBytes != profile.largestAssetBytes ||
	observedProfile.reusedOccurrenceCount !=
	    profile.reusedOccurrenceCount) {
	std::fprintf(stderr, "FAIL: immutable source profile publication\n");
	return 1;
    }
    BObolCompactSourceProfile conflictingProfile = profile;
    conflictingProfile.uniqueAssetCount++;
    conflictingProfile.reusedOccurrenceCount--;
    stream.setExpectedCount(1);
    stream.setSourceProfile(conflictingProfile);
    if (stream.getExpectedCount() != 100000 ||
	!stream.getSourceProfile(observedProfile) ||
	observedProfile.uniqueAssetCount != profile.uniqueAssetCount) {
	std::fprintf(stderr, "FAIL: stream discovery contract changed after certification\n");
	return 1;
    }
    stream.setWarmCensusComplete(true);
    stream.setWarmCoverageComplete(true);
    const SbBox3f conservativeBounds(SbVec3f(-20.0f, -30.0f, -40.0f),
	SbVec3f(50.0f, 60.0f, 70.0f));
    const SbBox3f exactBounds(SbVec3f(-10.0f, -20.0f, -30.0f),
	SbVec3f(40.0f, 50.0f, 60.0f));
    stream.setCoverageBounds(conservativeBounds);
    SbBox3f publishedBounds;
    if (stream.hasCoverageBoundsComplete() ||
	stream.hasCoverageBoundsDrained() ||
	!stream.getCoverageBounds(publishedBounds) ||
	publishedBounds != conservativeBounds) {
	std::fprintf(stderr,
	    "FAIL: conservative coverage was incorrectly terminal\n");
	return 1;
    }

    stream.push(occurrence(10));
    stream.push(occurrence(11));
    std::vector<BObolCompactOccurrence> pushedBatch;
    pushedBatch.push_back(occurrence(12));
    pushedBatch.push_back(occurrence(13));
    stream.pushBatch(std::move(pushedBatch));
    stream.pushPriority(occurrence(1));
    stream.pushPriority(occurrence(2));
    stream.push(occurrence(14));

    std::vector<BObolCompactOccurrence> overview;
    if (stream.drain(overview, 1) != 1 || overview.size() != 1 ||
	overview[0].occurrenceIndex != 2 ||
	!stream.hasCoverageBoundsDrained() ||
	stream.hasCoverageBoundsComplete()) {
	std::fprintf(stderr,
	    "FAIL: conservative whole-target overview publication\n");
	return 1;
    }

    stream.pushPriority(occurrence(3));
    stream.setCoverageBounds(exactBounds);
    stream.setCoverageBoundsComplete(true);

    if (stream.getExpectedCount() != 100000 ||
	!stream.hasWarmCensusComplete() ||
	!stream.hasWarmCoverageComplete() ||
	!stream.hasCoverageBoundsComplete() ||
	!stream.getCoverageBounds(publishedBounds) ||
	publishedBounds != exactBounds ||
	stream.hasCoverageBoundsDrained() ||
	stream.isCancelled() || stream.size() != 6) {
	std::fprintf(stderr, "FAIL: stream state\n");
	return 1;
    }

    std::vector<BObolCompactOccurrence> first;
    if (stream.drain(first, 3) != 3 || first.size() != 3 ||
	first[0].occurrenceIndex != 3 ||
	first[1].occurrenceIndex != 10 ||
	first[2].occurrenceIndex != 11 ||
	!stream.hasCoverageBoundsDrained() ||
	!stream.getCoverageBounds(publishedBounds) ||
	publishedBounds != exactBounds ||
	stream.size() != 3) {
	std::fprintf(stderr, "FAIL: priority drain order\n");
	return 1;
    }

    std::vector<BObolCompactOccurrence> second;
    if (stream.drain(second, 0) != 3 || second.size() != 3 ||
	second[0].occurrenceIndex != 12 ||
	second[1].occurrenceIndex != 13 ||
	second[2].occurrenceIndex != 14 ||
	stream.size() != 0) {
	std::fprintf(stderr, "FAIL: pending drain order\n");
	return 1;
    }

    stream.requestCancel();
    if (!stream.isCancelled()) {
	std::fprintf(stderr, "FAIL: cancellation publication\n");
	return 1;
    }
    return 0;
}

static int
test_warm_terminal_subset(void)
{
    BObolCompactOccurrenceStream stream;
    stream.setWarmCensusComplete(true);

    BObolCompactManifestOccurrence terminal;
    terminal.path = "all/analytic";
    terminal.sourceName = "analytic";
    terminal.bounds = SbBox3f(SbVec3f(-1.0f, -2.0f, -3.0f),
	SbVec3f(4.0f, 5.0f, 6.0f));
    stream.recordWarmTerminalOccurrence(terminal);

    BObolCompactManifestOccurrence lazyMesh = terminal;
    lazyMesh.path = "all/mesh";
    lazyMesh.sourceName = "mesh";
    lazyMesh.sourceMeshRequestValid = TRUE;
    stream.recordWarmTerminalOccurrence(lazyMesh);

    std::vector<BObolCompactManifestOccurrence> pending;
    if (!stream.hasWarmCensusComplete() ||
	stream.hasWarmCoverageComplete() ||
	!stream.takeWarmTerminalOccurrences(pending) ||
	pending.size() != 1 ||
	pending[0].path != terminal.path ||
	stream.takeWarmTerminalOccurrences(pending)) {
	std::fprintf(stderr, "FAIL: warm terminal subset contract\n");
	return 1;
    }

    stream.setWarmCoverageComplete(true);
    if (!stream.hasWarmCoverageComplete()) {
	std::fprintf(stderr, "FAIL: warm representation completion\n");
	return 1;
    }
    return 0;
}

static int
test_cancelled_publication(void)
{
    BObolCompactOccurrenceStream stream;
    BObolCompactOccurrence leaf = occurrence(1);
    leaf.summary.valid = TRUE;
    leaf.summary.boundsValid = TRUE;
    leaf.summary.path = "root/leaf";
    leaf.summary.sourceName = "leaf";
    leaf.summary.bounds = SbBox3f(SbVec3f(-1, -1, -1), SbVec3f(1, 1, 1));
    stream.push(leaf);
    stream.pushPriority(leaf);
    stream.recordManifestOccurrence(leaf);
    if (!stream.sealManifest(1)) {
	std::fprintf(stderr, "FAIL: cancelled manifest setup\n");
	return 1;
    }
    BObolCompactManifestOccurrence terminal;
    terminal.path = leaf.summary.path;
    terminal.sourceName = leaf.summary.sourceName;
    terminal.bounds = leaf.summary.bounds;
    stream.recordWarmTerminalOccurrence(terminal);
    BObolCompactSourceProfile profile;
    profile.valid = TRUE;
    profile.occurrenceCount = 1;
    profile.uniqueAssetCount = 1;
    profile.encodedSourceBytes = 1024;
    profile.largestAssetBytes = 1024;
    stream.setExpectedCount(1);
    stream.setSourceProfile(profile);
    stream.requestCancel();
    stream.requestCancel();

    /* A producer can already have passed its cancellation check. The stream
     * owns the publication boundary and must reject every late writer. */
    stream.push(leaf);
    stream.push(occurrence(2));
    stream.pushBatch({leaf, leaf});
    stream.pushPriority(leaf);
    stream.recordManifestOccurrence(leaf);
    stream.recordWarmTerminalOccurrence(terminal);
    BObolCompactSourceProfile conflictingProfile = profile;
    conflictingProfile.occurrenceCount = 2;
    conflictingProfile.reusedOccurrenceCount = 1;
    stream.setExpectedCount(2);
    stream.setSourceProfile(conflictingProfile);
    std::vector<BObolCompactManifestOccurrence> manifest;
    std::vector<BObolCompactOccurrence> drained;
    BObolCompactSourceProfile retainedProfile;
    if (!stream.isCancelled() || stream.size() ||
	stream.drain(drained, 0) || stream.sealManifest(1) ||
	stream.takeManifest(manifest) || stream.takeWarmTerminalOccurrences(manifest) ||
	stream.getExpectedCount() != 1 ||
	!stream.getSourceProfile(retainedProfile) ||
	retainedProfile.occurrenceCount != 1) {
	std::fprintf(stderr, "FAIL: cancelled stream retained or accepted publication\n");
	return 1;
    }
    return 0;
}

static int
test_cancelled_owner_release(void)
{
    BObolCompactOccurrenceStream stream;
    bool releasedOutsideLock = false;
    struct Import {
	BObolCompactOccurrenceStream *stream;
	bool *released;
	point_t points[3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
	int faces[3] = {0, 1, 2};
	~Import(void)
	{
	    /* Querying here deadlocks if cancellation destroys owners under the
	     * stream mutex. The storage owner can have independent cleanup. */
	    *released = stream->isCancelled() && stream->size() == 0;
	}
    };
    auto imported = std::make_shared<Import>();
    imported->stream = &stream;
    imported->released = &releasedOutsideLock;
    auto staged = std::make_shared<BObolStagedSourceMesh>();
    staged->owner = imported;
    staged->points = imported->points;
    staged->faces = imported->faces;
    staged->pointCount = 3;
    staged->faceCount = 1;
    staged->assetName = "release-probe";
    staged->byteCount = sizeof(Import);
    if (!stream.retainStagedSource(staged)) {
	std::fprintf(stderr, "FAIL: cancellation release probe setup\n");
	return 1;
    }
    imported.reset();
    staged.reset();
    stream.requestCancel();
    if (!releasedOutsideLock || stream.stagedSourceByteCount()) {
	std::fprintf(stderr, "FAIL: cancellation did not release import ownership\n");
	return 1;
    }
    return 0;
}

static int
test_concurrent_producers(void)
{
    const size_t producerCount = 4;
    const size_t itemsPerProducer = 2000;
    const size_t total = producerCount * itemsPerProducer;

    BObolCompactOccurrenceStream stream;
    stream.setExpectedCount(total);
    std::atomic<size_t> producersDone {0};
    std::vector<std::thread> producers;
    producers.reserve(producerCount);
    for (size_t producer = 0; producer < producerCount; producer++) {
	producers.emplace_back([producer, itemsPerProducer, &stream,
			   &producersDone]() {
	    const size_t base = producer * itemsPerProducer;
	    const size_t batchSize = 31;
	    for (size_t first = 0; first < itemsPerProducer;
		 first += batchSize) {
		std::vector<BObolCompactOccurrence> batch;
		const size_t count = std::min(batchSize,
		    itemsPerProducer - first);
		batch.reserve(count);
		for (size_t i = 0; i < count; ++i)
		    batch.push_back(occurrence(base + first + i));
		stream.pushBatch(std::move(batch));
	    }
	    producersDone.fetch_add(1, std::memory_order_release);
	});
    }

    std::vector<unsigned char> seen(total, 0);
    size_t consumed = 0;
    while (producersDone.load(std::memory_order_acquire) < producerCount ||
	stream.size() > 0) {
	std::vector<BObolCompactOccurrence> batch;
	stream.drain(batch, 127);
	for (const BObolCompactOccurrence &item : batch) {
	    const size_t index = item.occurrenceIndex;
	    if (index >= total || seen[index]) {
		std::fprintf(stderr,
		    "FAIL: duplicate/out-of-range concurrent item\n");
		for (std::thread &producer : producers)
		    producer.join();
		return 1;
	    }
	    seen[index] = 1;
	    consumed++;
	}
	if (batch.empty())
	    std::this_thread::yield();
    }

    for (std::thread &producer : producers)
	producer.join();
    if (consumed != total || stream.size() != 0) {
	std::fprintf(stderr,
	    "FAIL: concurrent stream lost items (%zu/%zu)\n",
	    consumed, total);
	return 1;
    }
    return 0;
}

int
main(int argc, char **argv)
{
    (void)argc;
    bu_setprogname(argv[0]);
    if (test_priority_and_state())
	return 1;
    if (test_warm_terminal_subset())
	return 1;
    if (test_cancelled_publication())
	return 1;
    if (test_cancelled_owner_release())
	return 1;
    if (test_concurrent_producers())
	return 1;
    return 0;
}
