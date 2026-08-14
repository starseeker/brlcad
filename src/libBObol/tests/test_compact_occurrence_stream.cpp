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
    stream.setWarmCoverageComplete(true);
    const SbBox3f exactBounds(SbVec3f(-10.0f, -20.0f, -30.0f),
	SbVec3f(40.0f, 50.0f, 60.0f));
    stream.setCoverageBounds(exactBounds);

    stream.push(occurrence(10));
    stream.push(occurrence(11));
    std::vector<BObolCompactOccurrence> pushedBatch;
    pushedBatch.push_back(occurrence(12));
    pushedBatch.push_back(occurrence(13));
    stream.pushBatch(std::move(pushedBatch));
    stream.pushPriority(occurrence(1));
    stream.pushPriority(occurrence(2));
    stream.push(occurrence(14));
    stream.setCoverageBoundsComplete(true);

    SbBox3f publishedBounds;

    if (stream.getExpectedCount() != 100000 ||
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
	first[0].occurrenceIndex != 2 ||
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
    if (test_concurrent_producers())
	return 1;
    return 0;
}
