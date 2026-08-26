/*         C O M P A C T _ O C C U R R E N C E _ S T R E A M . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BDatabaseSource.h"

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

struct BObolCompactOccurrenceStream::Impl {
    std::mutex mutex;
    std::vector<BObolCompactOccurrence> priority;
    size_t priorityOffset = 0;
    bool coverageOverviewQueued = false;
    bool coverageOverviewDrained = false;
    /* Producer-local vectors become queue nodes without moving their rich
     * occurrence values.  The former monolithic vector repeatedly relocated
     * every live record as a 150k-leaf cold stream grew concurrently with
     * GUI drains. */
    std::deque<std::vector<BObolCompactOccurrence>> pendingBatches;
    size_t pendingBatchOffset = 0;
    size_t pendingCount = 0;
    std::deque<std::shared_ptr<const BObolStagedSourceMesh>> stagedSources;
    size_t stagedSourceBytes = 0;
    /* Geometry-free journal used only to create the next warm-start manifest.
     * The GUI may drain pendingBatches before realization completes, so the
     * queue itself cannot be the persistence authority. */
    std::vector<BObolCompactManifestOccurrence> manifestOccurrences;
    std::unordered_map<std::string, size_t> manifestIndexByPath;
    bool manifestComplete = false;
    /* A persisted leaf manifest already supplied every leaf AABB and immutable
     * source-mesh request.  The authoritative semantics walk may therefore
     * skip its otherwise redundant full-BoT coverage import pass. */
    std::atomic<bool> warmCoverageComplete {false};
    std::atomic<bool> coverageBoundsComplete {false};
    SbBox3f coverageBounds;
    std::atomic<bool> cancelled {false};
    std::atomic<size_t> expectedCount {0};
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

BObolCompactOccurrenceStream::~BObolCompactOccurrenceStream(void)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    compact_stream_staged_source_release(this->d->stagedSourceBytes);
    this->d->stagedSourceBytes = 0;
    this->d->stagedSources.clear();
}

BObolCompactManifestOccurrence::BObolCompactManifestOccurrence(void) :
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
    if (this->d->pendingBatches.empty() ||
	this->d->pendingBatches.back().size() >= 64)
	this->d->pendingBatches.emplace_back();
    this->d->pendingBatches.back().push_back(occurrence);
    this->d->pendingCount++;
}

void
BObolCompactOccurrenceStream::push(BObolCompactOccurrence &&occurrence)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    if (this->d->pendingBatches.empty() ||
	this->d->pendingBatches.back().size() >= 64)
	this->d->pendingBatches.emplace_back();
    this->d->pendingBatches.back().push_back(std::move(occurrence));
    this->d->pendingCount++;
}

void
BObolCompactOccurrenceStream::pushBatch(
    std::vector<BObolCompactOccurrence> &&occurrences)
{
    if (occurrences.empty())
	return;
    std::lock_guard<std::mutex> guard(this->d->mutex);
    const size_t count = occurrences.size();
    this->d->pendingCount += count;
    this->d->pendingBatches.push_back(std::move(occurrences));
}

void
BObolCompactOccurrenceStream::pushPriority(
    const BObolCompactOccurrence &occurrence)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);

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
    this->d->priority.clear();
    this->d->priorityOffset = 0;
    this->d->priority.push_back(occurrence);
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
    if (this->d->manifestComplete)
	return;
    const auto found = this->d->manifestIndexByPath.find(key);
    if (found == this->d->manifestIndexByPath.end()) {
	const size_t index = this->d->manifestOccurrences.size();
	this->d->manifestOccurrences.push_back(std::move(record));
	this->d->manifestIndexByPath.emplace(key, index);
    } else if (found->second < this->d->manifestOccurrences.size()) {
	this->d->manifestOccurrences[found->second] = std::move(record);
    }
}

bool
BObolCompactOccurrenceStream::sealManifest(size_t expectedCount)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    bool complete = expectedCount > 0 &&
	this->d->manifestOccurrences.size() == expectedCount;
    for (const BObolCompactManifestOccurrence &record :
	 this->d->manifestOccurrences) {
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
    this->d->manifestComplete = complete;
    return complete;
}

bool
BObolCompactOccurrenceStream::takeManifest(
    std::vector<BObolCompactManifestOccurrence> &occurrences)
{
    occurrences.clear();
    std::lock_guard<std::mutex> guard(this->d->mutex);
    if (!this->d->manifestComplete)
	return false;
    occurrences.swap(this->d->manifestOccurrences);
    this->d->manifestIndexByPath.clear();
    this->d->manifestComplete = false;
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
    /* Keep an exceptional source larger than the ordinary window only while
     * it is the sole lease.  This enables a Lucy/one-huge-part handoff without
     * letting a many-leaf coverage pass retain several exceptional imports. */
    while (!this->d->stagedSources.empty() &&
	(source->byteCount > limit ||
	 this->d->stagedSourceBytes > limit - source->byteCount)) {
	const std::shared_ptr<const BObolStagedSourceMesh> &oldest =
	    this->d->stagedSources.front();
	const size_t bytes = oldest ? oldest->byteCount : 0;
	this->d->stagedSourceBytes =
	    bytes >= this->d->stagedSourceBytes ?
	    0 : this->d->stagedSourceBytes - bytes;
	this->d->stagedSources.pop_front();
	compact_stream_staged_source_release(bytes);
    }
    if (!compact_stream_staged_source_reserve(source->byteCount))
	return FALSE;
    if (source->byteCount <= SIZE_MAX - this->d->stagedSourceBytes)
	this->d->stagedSourceBytes += source->byteCount;
    else
	this->d->stagedSourceBytes = SIZE_MAX;
    this->d->stagedSources.push_back(source);
    return TRUE;
}

std::shared_ptr<const BObolStagedSourceMesh>
BObolCompactOccurrenceStream::claimStagedSource(
    const std::weak_ptr<const BObolStagedSourceMesh> &source)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    const std::shared_ptr<const BObolStagedSourceMesh> requested =
	source.lock();
    if (!requested)
	return std::shared_ptr<const BObolStagedSourceMesh>();

    const auto found = std::find_if(this->d->stagedSources.begin(),
	this->d->stagedSources.end(),
	[&requested](
	    const std::shared_ptr<const BObolStagedSourceMesh> &candidate) {
	    return candidate.get() == requested.get();
	});
    if (found == this->d->stagedSources.end())
	return std::shared_ptr<const BObolStagedSourceMesh>();

    std::shared_ptr<const BObolStagedSourceMesh> claimed = *found;
    const size_t bytes = claimed ? claimed->byteCount : 0;
    this->d->stagedSourceBytes =
	bytes >= this->d->stagedSourceBytes ?
	0 : this->d->stagedSourceBytes - bytes;
    this->d->stagedSources.erase(found);
    compact_stream_staged_source_release(bytes);
    return claimed;
}

size_t
BObolCompactOccurrenceStream::stagedSourceByteCount(void)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    return this->d->stagedSourceBytes;
}

size_t
BObolCompactOccurrenceStream::drain(
    std::vector<BObolCompactOccurrence> &out, size_t cap)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    const size_t priorityAvailable =
	this->d->priority.size() - this->d->priorityOffset;
    const size_t pendingAvailable = this->d->pendingCount;
    const size_t available = priorityAvailable + pendingAvailable;
    const size_t count =
	(cap == 0 || cap >= available) ? available : cap;
    if (!count)
	return 0;
    out.reserve(out.size() + count);
    const size_t priorityCount = std::min(count, priorityAvailable);
    for (size_t i = 0; i < priorityCount; i++)
	out.push_back(std::move(
	    this->d->priority[this->d->priorityOffset + i]));
    this->d->priorityOffset += priorityCount;
    if (priorityCount &&
	this->d->priorityOffset == this->d->priority.size())
	this->d->coverageOverviewDrained = true;
    size_t pendingToDrain = count - priorityCount;
    while (pendingToDrain && !this->d->pendingBatches.empty()) {
	std::vector<BObolCompactOccurrence> &batch =
	    this->d->pendingBatches.front();
	const size_t batchAvailable = batch.size() -
	    this->d->pendingBatchOffset;
	const size_t batchCount = std::min(pendingToDrain, batchAvailable);
	for (size_t i = 0; i < batchCount; ++i)
	    out.push_back(std::move(
		batch[this->d->pendingBatchOffset + i]));
	this->d->pendingBatchOffset += batchCount;
	pendingToDrain -= batchCount;
	this->d->pendingCount = batchCount > this->d->pendingCount ?
	    0 : this->d->pendingCount - batchCount;
	if (this->d->pendingBatchOffset == batch.size()) {
	    this->d->pendingBatches.pop_front();
	    this->d->pendingBatchOffset = 0;
	}
    }

    /* Priority contains at most the newest aggregate extent in normal use.
     * Keep its cursor logic independent of the queued leaf batches so replacing
     * an overview never disturbs producer-owned occurrence storage. */
    if (this->d->priorityOffset == this->d->priority.size()) {
	this->d->priority.clear();
	this->d->priorityOffset = 0;
    } else if (this->d->priorityOffset >= 64 &&
	this->d->priorityOffset >= this->d->priority.size() / 2) {
	this->d->priority.erase(this->d->priority.begin(),
	    this->d->priority.begin() + this->d->priorityOffset);
	this->d->priorityOffset = 0;
    }
    return count;
}

size_t
BObolCompactOccurrenceStream::size(void)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    return
	(this->d->priority.size() - this->d->priorityOffset) +
	this->d->pendingCount;
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
    this->d->expectedCount.store(count, std::memory_order_release);
}

size_t
BObolCompactOccurrenceStream::getExpectedCount(void) const
{
    return this->d->expectedCount.load(std::memory_order_acquire);
}

void
BObolCompactOccurrenceStream::setSourceProfile(
    const BObolCompactSourceProfile &profile)
{
    std::lock_guard<std::mutex> guard(this->d->mutex);
    /* A complete profile is monotonic within one stream epoch. */
    if (!this->d->sourceProfile.valid || profile.valid)
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

void
BObolCompactOccurrenceStream::requestCancel(void)
{
    this->d->cancelled.store(true, std::memory_order_release);
}

bool
BObolCompactOccurrenceStream::isCancelled(void) const
{
    return this->d->cancelled.load(std::memory_order_acquire);
}
