/*         C O M P A C T _ O C C U R R E N C E _ S T R E A M . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"
#include "database_source_realization.h"
#include "bu/str.h"

#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <deque>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

static void compact_stream_staged_source_release(size_t bytes);

struct BObolCompactOccurrenceStream::Impl {
    struct Payload {
	~Payload(void)
	{
	    stagedSources.clear();
	    compact_stream_staged_source_release(stagedSourceBytes);
	}

	std::vector<BObolCompactOccurrence> priority;
	size_t priorityOffset = 0;
	/* Producer-local vectors become queue nodes without relocating every
	 * rich record in a growing cold stream. */
	std::deque<std::vector<BObolCompactOccurrence>> pendingBatches;
	size_t pendingBatchOffset = 0;
	size_t pendingCount = 0;
	std::deque<std::shared_ptr<const BObolStagedSourceMesh>> stagedSources;
	size_t stagedSourceBytes = 0;
	/* Persistence needs its own mesh-buffer-free journal: the consumer may
	 * drain the geometry queue before source production completes. */
	std::vector<BObolCompactManifestOccurrence> manifestOccurrences;
	std::unordered_map<std::string, size_t> manifestIndexByPath;
	bool manifestComplete = false;
	std::vector<BObolCompactManifestOccurrence> warmTerminalOccurrences;
    };

    std::mutex mutex;
    /* Cancellation detaches all unconsumed ownership without allocating.
     * Destruction happens outside the stream lock. The small progress facts
     * remain readable for diagnostics while sibling sources finish. */
    std::unique_ptr<Payload> payload{new Payload};
    bool coverageOverviewQueued = false;
    bool coverageOverviewDrained = false;
    std::atomic<bool> warmCensusComplete {false};
    std::atomic<bool> warmCoverageComplete {false};
    std::atomic<bool> coverageBoundsComplete {false};
    SbBox3f coverageBounds;
    std::atomic<bool> cancelled {false};
    std::weak_ptr<BObolSourceRealizationCoordinatorPrivate> sourceCoordinator;
    std::atomic<size_t> expectedCount {0};
    std::atomic<size_t> preparationWorkCompleted {0};
    std::atomic<size_t> preparationWorkTotal {0};
    BObolCompactSourceProfile sourceProfile;
};

static size_t
compact_stream_staged_source_limit(void)
{
    static const size_t limit = []() {
	const size_t mebibyte = 1024ULL * 1024ULL;
	const char *configured = getenv("BOBOL_LOD_STAGED_SOURCE_MB");
	if (configured && configured[0]) {
	    char *end = NULL;
	    const unsigned long long value = strtoull(configured, &end, 10);
	    if (end && end != configured && *end == '\0') {
		if (value == 0)
		    return static_cast<size_t>(0);
		if (value > SIZE_MAX / mebibyte)
		    return SIZE_MAX;
		return static_cast<size_t>(value) * mebibyte;
	    }
	}
	return static_cast<size_t>(512ULL * mebibyte);
    }();
    return limit;
}

/*
 * Lock order:
 *
 *   stream Impl::mutex -> compact_stream_staged_source_budget_mutex
 *
 * The global-budget helpers never acquire a stream mutex, and none of the
 * stream methods invokes a callback while either lock is held.
 */
static std::mutex compact_stream_staged_source_budget_mutex;
static size_t compact_stream_staged_source_budget_bytes = 0;

static bool
compact_stream_staged_source_reserve(size_t bytes)
{
    const size_t limit = compact_stream_staged_source_limit();
    if (!bytes || !limit)
	return false;
    std::lock_guard<std::mutex> guard(
	compact_stream_staged_source_budget_mutex);
    if (bytes > limit) {
	if (compact_stream_staged_source_budget_bytes != 0)
	    return false;
	compact_stream_staged_source_budget_bytes = bytes;
	return true;
    }
    if (compact_stream_staged_source_budget_bytes > limit - bytes)
	return false;
    compact_stream_staged_source_budget_bytes += bytes;
    return true;
}

static void
compact_stream_staged_source_release(size_t bytes)
{
    std::lock_guard<std::mutex> guard(
	compact_stream_staged_source_budget_mutex);
    compact_stream_staged_source_budget_bytes =
	bytes >= compact_stream_staged_source_budget_bytes ?
	0 : compact_stream_staged_source_budget_bytes - bytes;
}

BObolCompactOccurrenceStream::BObolCompactOccurrenceStream(void) :
    d(new Impl)
{
}

BObolCompactOccurrenceStream::~BObolCompactOccurrenceStream(void) = default;

BObolCompactManifestOccurrence::BObolCompactManifestOccurrence(void) :
    orientedBoundsValid(FALSE),
    booleanOperation(SoBRLDatabaseSource::BOOLEAN_UNION), occurrenceIndex(0),
    sourceMeshRequestValid(FALSE), meshAssetContentHash(0),
    meshAssetTessellationAbsTol(0.0), meshAssetTessellationRelTol(0.0),
    meshAssetTessellationNormTol(0.0), sourceFaceCount(0), sourcePointCount(0),
    regionId(0), airCode(0), materialId(0), los(0),
    materialColorValid(FALSE), materialColor(0.0f, 0.0f, 0.0f)
{
    localTransform.makeIdentity();
    bounds.makeEmpty();
    meshAssetBounds.makeEmpty();
    meshAssetTransform.makeIdentity();
}

void
BObolCompactOccurrenceStream::push(
    const BObolCompactOccurrence &occurrence)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return;
    if (payload->pendingBatches.empty() ||
	payload->pendingBatches.back().size() >= 64)
	payload->pendingBatches.emplace_back();
    payload->pendingBatches.back().push_back(occurrence);
    payload->pendingCount++;
}

void
BObolCompactOccurrenceStream::push(BObolCompactOccurrence &&occurrence)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return;
    if (payload->pendingBatches.empty() ||
	payload->pendingBatches.back().size() >= 64)
	payload->pendingBatches.emplace_back();
    payload->pendingBatches.back().push_back(std::move(occurrence));
    payload->pendingCount++;
}

void
BObolCompactOccurrenceStream::pushBatch(
    std::vector<BObolCompactOccurrence> &&occurrences)
{
    if (occurrences.empty())
	return;
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return;
    const size_t count = occurrences.size();
    payload->pendingBatches.push_back(std::move(occurrences));
    payload->pendingCount += count;
}

void
BObolCompactOccurrenceStream::pushPriority(
    const BObolCompactOccurrence &occurrence)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return;

    /*
     * This lane is the current whole-target extent for one realization
     * stream, not an event history.  A bounds worker may publish provisional
     * snapshots while discovery is running; retaining all of them puts the
     * final exact overview behind an increasingly stale priority backlog.
     * The owner would then know the exact source bounds but deliberately
     * defer autoview until it had merged every obsolete box.  Keep only the
     * newest undrained snapshot.  An occurrence already moved into a consumer
     * batch is independent and may finish its merge safely.
     */
    payload->priority.clear();
    payload->priorityOffset = 0;
    payload->priority.push_back(occurrence);
    this->d->coverageOverviewQueued = true;
    this->d->coverageOverviewDrained = false;
}

void
BObolCompactOccurrenceStream::recordManifestOccurrence(
    const BObolCompactOccurrence &occurrence)
{
    if (!occurrence.summary.valid || !occurrence.summary.boundsValid ||
	occurrence.summary.bounds.isEmpty() ||
	occurrence.summary.path.getLength() == 0 ||
	BU_STR_EQUAL(occurrence.summary.recordRole.getString(), "lod-overview"))
	return;

    BObolCompactManifestOccurrence record;
    record.path = occurrence.summary.path;
    record.sourceName = occurrence.summary.sourceName;
    record.localTransform = occurrence.localTransform;
    record.bounds = occurrence.summary.bounds;
    if (occurrence.geometry && occurrence.geometry->aggregateProxyCorners) {
	record.orientedBoundsValid = TRUE;
	record.orientedBounds = *occurrence.geometry->aggregateProxyCorners;
    }
    record.booleanOperation = occurrence.booleanOperation;
    record.occurrenceIndex = occurrence.occurrenceIndex;
    record.regionId = occurrence.summary.regionId;
    record.airCode = occurrence.summary.airCode;
    record.materialId = occurrence.summary.materialId;
    record.los = occurrence.summary.los;
    record.materialColorValid = occurrence.summary.materialColorValid;
    record.materialColor = occurrence.summary.materialColor;
    record.materialShader = occurrence.summary.materialShader;
    record.sourceMeshRequestValid = occurrence.sourceMeshRequestValid;
    if (record.sourceMeshRequestValid) {
	const BObolSourceMeshRequest &request = occurrence.sourceMeshRequest;
	record.sourceType = request.sourceType;
	record.meshAssetPath = request.meshAssetPath.getLength() > 0 ?
	    request.meshAssetPath : occurrence.summary.path;
	record.meshAssetName = request.meshAssetName.getLength() > 0 ?
	    request.meshAssetName : occurrence.summary.sourceName;
	record.meshAssetContentHash = request.meshAssetContentHash;
	record.meshAssetTessellationAbsTol =
	    request.meshAssetTessellationAbsTol;
	record.meshAssetTessellationRelTol =
	    request.meshAssetTessellationRelTol;
	record.meshAssetTessellationNormTol =
	    request.meshAssetTessellationNormTol;
	record.meshAssetBounds = !request.meshAssetBounds.isEmpty() ?
	    request.meshAssetBounds :
	    (!request.bounds.isEmpty() ? request.bounds : record.bounds);
	record.meshAssetTransform = request.meshAssetTransform;
	record.sourceFaceCount = request.faceCount;
	record.sourcePointCount = request.pointCount;
    }

    const std::string key = record.path.getString();
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return;
    if (payload->manifestComplete)
	return;
    const auto found = payload->manifestIndexByPath.find(key);
    if (found == payload->manifestIndexByPath.end()) {
	const size_t index = payload->manifestOccurrences.size();
	payload->manifestOccurrences.push_back(std::move(record));
	payload->manifestIndexByPath.emplace(key, index);
    } else if (found->second < payload->manifestOccurrences.size()) {
	payload->manifestOccurrences[found->second] = std::move(record);
    }
}

bool
BObolCompactOccurrenceStream::sealManifest(size_t expectedCount)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return false;
    bool complete = expectedCount > 0 &&
	payload->manifestOccurrences.size() == expectedCount;
    for (const BObolCompactManifestOccurrence &record :
	 payload->manifestOccurrences) {
	if (!complete)
	    break;
	complete = record.path.getLength() > 0 &&
	    record.sourceName.getLength() > 0 && !record.bounds.isEmpty();
	if (complete && record.sourceMeshRequestValid) {
	    complete = record.meshAssetPath.getLength() > 0 &&
		record.meshAssetName.getLength() > 0 &&
		!record.meshAssetBounds.isEmpty();
	}
    }
    payload->manifestComplete = complete;
    return complete;
}

bool
BObolCompactOccurrenceStream::takeManifest(
    std::vector<BObolCompactManifestOccurrence> &occurrences)
{
    occurrences.clear();
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return false;
    if (!payload->manifestComplete)
	return false;
    occurrences.swap(payload->manifestOccurrences);
    payload->manifestIndexByPath.clear();
    payload->manifestComplete = false;
    return !occurrences.empty();
}

SbBool
BObolCompactOccurrenceStream::retainStagedSource(
    const std::shared_ptr<const BObolStagedSourceMesh> &source)
{
    if (!source || !source->isValid() || !source->byteCount)
	return FALSE;
    const size_t limit = compact_stream_staged_source_limit();
    if (!limit)
	return FALSE;

    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return FALSE;
    /* Keep an exceptional source larger than the ordinary window only while
     * it is the sole lease.  This enables a Lucy/one-huge-part handoff without
     * letting a many-leaf coverage pass retain several exceptional imports. */
    while (!payload->stagedSources.empty() &&
	(source->byteCount > limit ||
	 payload->stagedSourceBytes > limit - source->byteCount)) {
	const std::shared_ptr<const BObolStagedSourceMesh> &oldest =
	    payload->stagedSources.front();
	const size_t bytes = oldest ? oldest->byteCount : 0;
	payload->stagedSourceBytes =
	    bytes >= payload->stagedSourceBytes ?
	    0 : payload->stagedSourceBytes - bytes;
	payload->stagedSources.pop_front();
	compact_stream_staged_source_release(bytes);
    }
    if (!compact_stream_staged_source_reserve(source->byteCount))
	return FALSE;
    if (source->byteCount <= SIZE_MAX - payload->stagedSourceBytes)
	payload->stagedSourceBytes += source->byteCount;
    else
	payload->stagedSourceBytes = SIZE_MAX;
    payload->stagedSources.push_back(source);
    return TRUE;
}

std::shared_ptr<const BObolStagedSourceMesh>
BObolCompactOccurrenceStream::claimStagedSource(
    const std::weak_ptr<const BObolStagedSourceMesh> &source)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return std::shared_ptr<const BObolStagedSourceMesh>();
    const std::shared_ptr<const BObolStagedSourceMesh> requested =
	source.lock();
    if (!requested)
	return std::shared_ptr<const BObolStagedSourceMesh>();

    const auto found = std::find_if(payload->stagedSources.begin(),
	payload->stagedSources.end(),
	[&requested](
	    const std::shared_ptr<const BObolStagedSourceMesh> &candidate) {
	    return candidate.get() == requested.get();
	});
    if (found == payload->stagedSources.end())
	return std::shared_ptr<const BObolStagedSourceMesh>();

    std::shared_ptr<const BObolStagedSourceMesh> claimed = *found;
    const size_t bytes = claimed ? claimed->byteCount : 0;
    payload->stagedSourceBytes =
	bytes >= payload->stagedSourceBytes ?
	0 : payload->stagedSourceBytes - bytes;
    payload->stagedSources.erase(found);
    compact_stream_staged_source_release(bytes);
    return claimed;
}

size_t
BObolCompactOccurrenceStream::stagedSourceByteCount(void)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return 0;
    return payload->stagedSourceBytes;
}

size_t
BObolCompactOccurrenceStream::drain(
    std::vector<BObolCompactOccurrence> &out, size_t cap)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return 0;
    const size_t priorityAvailable =
	payload->priority.size() - payload->priorityOffset;
    const size_t pendingAvailable = payload->pendingCount;
    const size_t available = priorityAvailable + pendingAvailable;
    const size_t count =
	(cap == 0 || cap >= available) ? available : cap;
    if (!count)
	return 0;
    out.reserve(out.size() + count);
    const size_t priorityCount = std::min(count, priorityAvailable);
    for (size_t i = 0; i < priorityCount; i++)
	out.push_back(std::move(
	    payload->priority[payload->priorityOffset + i]));
    payload->priorityOffset += priorityCount;
    if (priorityCount &&
	payload->priorityOffset == payload->priority.size())
	this->d->coverageOverviewDrained = true;
    size_t pendingToDrain = count - priorityCount;
    while (pendingToDrain && !payload->pendingBatches.empty()) {
	std::vector<BObolCompactOccurrence> &batch =
	    payload->pendingBatches.front();
	const size_t batchAvailable = batch.size() -
	    payload->pendingBatchOffset;
	const size_t batchCount = std::min(pendingToDrain, batchAvailable);
	for (size_t i = 0; i < batchCount; ++i)
	    out.push_back(std::move(
		batch[payload->pendingBatchOffset + i]));
	payload->pendingBatchOffset += batchCount;
	pendingToDrain -= batchCount;
	payload->pendingCount = batchCount > payload->pendingCount ?
	    0 : payload->pendingCount - batchCount;
	if (payload->pendingBatchOffset == batch.size()) {
	    payload->pendingBatches.pop_front();
	    payload->pendingBatchOffset = 0;
	}
    }

    /* Priority contains at most the newest aggregate extent in normal use.
     * Keep its cursor logic independent of the queued leaf batches so replacing
     * an overview never disturbs producer-owned occurrence storage. */
    if (payload->priorityOffset == payload->priority.size()) {
	payload->priority.clear();
	payload->priorityOffset = 0;
    } else if (payload->priorityOffset >= 64 &&
	payload->priorityOffset >= payload->priority.size() / 2) {
	payload->priority.erase(payload->priority.begin(),
	    payload->priority.begin() + payload->priorityOffset);
	payload->priorityOffset = 0;
    }
    return count;
}

size_t
BObolCompactOccurrenceStream::size(void)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return 0;
    return
	(payload->priority.size() - payload->priorityOffset) +
	payload->pendingCount;
}

void
BObolCompactOccurrenceStream::setExpectedCount(size_t count)
{
    /* Expected population is progress metadata, not a queue-capacity demand.
     * A cold producer learns this value only after its hierarchy walk, while
     * the owner is already draining completed leaves.  Reserving the total at
     * that point copies the entire live backlog and allocates space for tens
     * of thousands of records which will never coexist.  pushBatch() queues
     * only the actual producer/consumer backlog. */
    if (!count)
	return;
    std::lock_guard<std::mutex> guard(this->d->mutex);
    if (this->d->cancelled.load(std::memory_order_acquire) ||
	(this->d->sourceProfile.valid &&
	 this->d->sourceProfile.occurrenceCount !=
	     static_cast<uint64_t>(count)))
	return;
    size_t unset = 0;
    (void)this->d->expectedCount.compare_exchange_strong(unset, count,
	std::memory_order_release, std::memory_order_acquire);
}

size_t
BObolCompactOccurrenceStream::getExpectedCount(void) const
{
    return this->d->expectedCount.load(std::memory_order_acquire);
}

void
BObolCompactOccurrenceStream::setPreparationWorkCount(size_t count)
{
    this->d->preparationWorkCompleted.store(0, std::memory_order_release);
    this->d->preparationWorkTotal.store(count, std::memory_order_release);
}

void
BObolCompactOccurrenceStream::notePreparationWorkCompleted(void)
{
    const size_t total = this->d->preparationWorkTotal.load(
	std::memory_order_acquire);
    if (!total)
	return;

    size_t completed = this->d->preparationWorkCompleted.load(
	std::memory_order_acquire);
    while (completed < total &&
	!this->d->preparationWorkCompleted.compare_exchange_weak(
	    completed, completed + 1, std::memory_order_release,
	    std::memory_order_acquire)) {
    }
}

void
BObolCompactOccurrenceStream::completePreparationWork(void)
{
    const size_t total = this->d->preparationWorkTotal.load(
	std::memory_order_acquire);
    this->d->preparationWorkCompleted.store(total,
	std::memory_order_release);
}

void
BObolCompactOccurrenceStream::getPreparationWorkCount(
    size_t &completed, size_t &total) const
{
    total = this->d->preparationWorkTotal.load(std::memory_order_acquire);
    completed = std::min(this->d->preparationWorkCompleted.load(
	std::memory_order_acquire), total);
}

void
BObolCompactOccurrenceStream::setSourceProfile(
    const BObolCompactSourceProfile &profile)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    if (this->d->cancelled.load(std::memory_order_acquire))
	return;
    const size_t expectedCount = this->d->expectedCount.load(
	std::memory_order_acquire);
    if (!profile.isValid(expectedCount))
	return;
    /* A profile describes one complete discovery and is immutable for the
     * stream epoch. Repeated publication is idempotent. */
    if (!this->d->sourceProfile.valid)
	this->d->sourceProfile = profile;
}

SbBool
BObolCompactOccurrenceStream::getSourceProfile(
    BObolCompactSourceProfile &profile) const
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    profile = this->d->sourceProfile;
    return profile.valid;
}

void
BObolCompactOccurrenceStream::setWarmCensusComplete(bool complete)
{
    this->d->warmCensusComplete.store(complete, std::memory_order_release);
}

bool
BObolCompactOccurrenceStream::hasWarmCensusComplete(void) const
{
    return this->d->warmCensusComplete.load(std::memory_order_acquire);
}

void
BObolCompactOccurrenceStream::setWarmCoverageComplete(bool complete)
{
    this->d->warmCoverageComplete.store(
	complete, std::memory_order_release);
}

bool
BObolCompactOccurrenceStream::hasWarmCoverageComplete(void) const
{
    return this->d->warmCoverageComplete.load(std::memory_order_acquire);
}

void
BObolCompactOccurrenceStream::recordWarmTerminalOccurrence(
    const BObolCompactManifestOccurrence &occurrence)
{
    if (occurrence.path.getLength() == 0 ||
	occurrence.sourceName.getLength() == 0 || occurrence.bounds.isEmpty() ||
	occurrence.sourceMeshRequestValid)
	return;
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return;
    payload->warmTerminalOccurrences.push_back(occurrence);
}

bool
BObolCompactOccurrenceStream::takeWarmTerminalOccurrences(
    std::vector<BObolCompactManifestOccurrence> &occurrences)
{
    occurrences.clear();
    std::lock_guard<std::mutex> guard(this->d->mutex);
    Impl::Payload *payload = this->d->payload.get();
    if (!payload)
	return false;
    occurrences.swap(payload->warmTerminalOccurrences);
    return !occurrences.empty();
}

void
BObolCompactOccurrenceStream::setCoverageBounds(const SbBox3f &bounds)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    this->d->coverageBounds = bounds;
}

bool
BObolCompactOccurrenceStream::getCoverageBounds(SbBox3f &bounds)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    bounds = this->d->coverageBounds;
    return !bounds.isEmpty();
}

void
BObolCompactOccurrenceStream::setCoverageBoundsComplete(bool complete)
{
    this->d->coverageBoundsComplete.store(
	complete, std::memory_order_release);
}

bool
BObolCompactOccurrenceStream::hasCoverageBoundsComplete(void) const
{
    return this->d->coverageBoundsComplete.load(
	std::memory_order_acquire);
}

bool
BObolCompactOccurrenceStream::hasCoverageBoundsDrained(void)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    return !this->d->coverageBounds.isEmpty() &&
	this->d->coverageOverviewQueued &&
	this->d->coverageOverviewDrained;
}

bool
BObolCompactOccurrenceStream::bindRealizationCoordinator(
    const std::shared_ptr<BObolSourceRealizationCoordinatorPrivate> &coordinator)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    this->d->sourceCoordinator = coordinator;
    return this->d->cancelled.load(std::memory_order_acquire);
}

std::weak_ptr<BObolSourceRealizationCoordinatorPrivate>
BObolCompactOccurrenceStream::cancelPublication(void)
{
    std::unique_ptr<Impl::Payload> retired;
    std::weak_ptr<BObolSourceRealizationCoordinatorPrivate> coordinator;
    {
	std::lock_guard<std::mutex> guard(this->d->mutex);
	if (this->d->cancelled.load(std::memory_order_acquire))
	    return coordinator;
	retired = std::move(this->d->payload);
	this->d->cancelled.store(true, std::memory_order_release);
	coordinator = this->d->sourceCoordinator;
    }
    retired.reset();
    return coordinator;
}

void
BObolCompactOccurrenceStream::requestCancel(void)
{
    bobol_source_realization_cancel_queued(this->cancelPublication());
}

bool
BObolCompactOccurrenceStream::isCancelled(void) const
{
    return this->d->cancelled.load(std::memory_order_acquire);
}
