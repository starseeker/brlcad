/*                 L O D _ T E L E M E T R Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_TELEMETRY_PRIVATE_H
#define LIBBOBOL_LOD_TELEMETRY_PRIVATE_H

#include "BObol/BLodRealization.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

class BObolViewController;
struct BObolLodControlTransitionRecord;

/* Privacy-safe summary of one owner-thread result drain.  It records only
 * finite enums, counters, replay causes, and source-local dense indices;
 * object names, paths, cache keys, and diagnostics never enter this record. */
struct BObolLodPublicationTelemetryRecord {
    static constexpr size_t providerStatusCount =
	static_cast<size_t>(BOBOL_LOD_PROVIDER_SUPERSEDED) + 1;
    static constexpr size_t maximumRetryEntryIndices = 64;
    static constexpr size_t maximumResultSamples = 64;

    struct ResultSample {
	uint64_t sourceRoutingId = 0;
	uint32_t sourceEntryIndex = UINT32_MAX;
	int submissionReason = BOBOL_LOD_SUBMISSION_UNSPECIFIED;
	int drawMode = BOBOL_LOD_DRAW_UNKNOWN;
	int normalStyle = BOBOL_LOD_NORMAL_AUTHORED;
	int requestCut = -1;
	int resolvedCut = -1;
	int incomingActiveCut = -1;
	int incomingResidentCut = -1;
	int incomingPresentationAdmissionCut = -1;
	size_t requestChunkCount = 0;
	uint64_t requestChunkHash = 0;
	size_t incomingLayerCount = 0;
	uint64_t incomingMeshRevision = 0;
	uint64_t incomingPreparedRevision = 0;
	bool incomingPreparedGeometry = false;
	uint64_t incomingFaceCount = 0;
	uint64_t incomingPointCount = 0;
	bool incomingTerminal = false;
	bool incomingMemoryLimited = false;
	bool residentBefore = false;
	int beforeActiveCut = -1;
	int beforeResidentCut = -1;
	int beforeRequestedCut = -1;
	int beforeAllocatedCut = -1;
	int beforeNormalStyle = BOBOL_LOD_NORMAL_AUTHORED;
	size_t beforeRequiredChunkCount = 0;
	uint64_t beforeRequiredChunkHash = 0;
	size_t beforePresentedChunkCount = 0;
	uint64_t beforePresentedChunkHash = 0;
	size_t beforeLayerCount = 0;
	uint64_t beforeMeshRevision = 0;
	uint64_t beforePreparedRevision = 0;
	bool beforePreparedGeometry = false;
	uint64_t beforeFaceCount = 0;
	uint64_t beforePointCount = 0;
	bool beforeNormalPresentationMatches = false;
	bool beforeRequestedCutDrawable = false;
	bool beforeAllocatedCutDrawable = false;
	bool sameProgressiveAsset = false;
	bool accepted = false;
	bool unchanged = false;
	bool retryCurrentDemand = false;
	bool residentAfter = false;
	int afterActiveCut = -1;
	int afterResidentCut = -1;
	int afterRequestedCut = -1;
	int afterAllocatedCut = -1;
	int afterNormalStyle = BOBOL_LOD_NORMAL_AUTHORED;
	size_t afterRequiredChunkCount = 0;
	uint64_t afterRequiredChunkHash = 0;
	size_t afterPresentedChunkCount = 0;
	uint64_t afterPresentedChunkHash = 0;
	size_t afterLayerCount = 0;
	uint64_t afterMeshRevision = 0;
	uint64_t afterPreparedRevision = 0;
	bool afterPreparedGeometry = false;
	uint64_t afterFaceCount = 0;
	uint64_t afterPointCount = 0;
	bool afterNormalPresentationMatches = false;
	bool afterRequestedCutDrawable = false;
	bool afterAllocatedCutDrawable = false;
	bool semanticStateChanged = false;
    };

    void noteProviderStatus(int status)
    {
	if (status >= 0 && static_cast<size_t>(status) <
		providerStatusCounts.size())
	    providerStatusCounts[static_cast<size_t>(status)]++;
	else
	    unknownProviderStatusCount++;
    }

    void noteRetryEntry(uint32_t entry)
    {
	retrySourceEntryCount++;
	if (entry != UINT32_MAX &&
	    retrySourceEntryIndices.size() < maximumRetryEntryIndices)
	    retrySourceEntryIndices.push_back(entry);
    }

    void noteResultSample(const ResultSample &sample)
    {
	resultSampleCount++;
	if (resultSamples.size() < maximumResultSamples)
	    resultSamples.push_back(sample);
    }

    size_t processed = 0;
    size_t matched = 0;
    size_t applied = 0;
    size_t rejected = 0;
    size_t unmatched = 0;
    std::array<size_t, providerStatusCount> providerStatusCounts = {};
    size_t unknownProviderStatusCount = 0;
    size_t authenticationPublishCount = 0;
    size_t authenticationTerminalFailureCount = 0;
    size_t authenticationRetryCount = 0;
    size_t authenticationSupersedeCount = 0;
    size_t sourceAcceptedCount = 0;
    size_t sourceUnchangedCount = 0;
    size_t sourceRejectedCount = 0;
    size_t sourceRetryCount = 0;
    size_t updateRetryCount = 0;
    size_t sourceRouteMismatchCount = 0;
    size_t sourcePopulationMismatchCount = 0;
    size_t demandMismatchCount = 0;
    size_t retrySourceEntryCount = 0;
    std::vector<uint32_t> retrySourceEntryIndices;
    size_t resultSampleCount = 0;
    std::vector<ResultSample> resultSamples;
    bool replayFromAuthentication = false;
    bool replayFromSource = false;
    bool replayFromUpdate = false;
    bool replayFromRetainedPublication = false;
    bool replayFromPartialRefinement = false;
    bool actionableQualityDebt = false;
    bool replayRequested = false;
    bool submissionActiveBeforeReplay = false;
    bool submissionActiveAfterReplay = false;
    bool rescanPendingAfterReplay = false;
    size_t submissionSourceIndexBeforeReplay = 0;
    size_t submissionEntryOffsetBeforeReplay = 0;
    size_t submissionSourceIndexAfterReplay = 0;
    size_t submissionEntryOffsetAfterReplay = 0;
};

/* Process-shared JSONL output with controller-local inventory cursors and
 * privacy aliases.  Every public operation is best-effort: diagnostic I/O
 * must never change the drawing control path. */
class BObolLodTelemetry {
public:
    static std::unique_ptr<BObolLodTelemetry> fromEnvironment(void) noexcept;
    ~BObolLodTelemetry(void);

    BObolLodTelemetry(const BObolLodTelemetry &) = delete;
    BObolLodTelemetry &operator=(const BObolLodTelemetry &) = delete;

    void recordTransition(const BObolViewController &controller,
	const BObolLodControlTransitionRecord &transition) noexcept;
    void recordPublication(
	const BObolLodPublicationTelemetryRecord &publication) noexcept;
    void recordGap(void) noexcept;

private:
    class Impl;
    explicit BObolLodTelemetry(std::unique_ptr<Impl> implementation);
    std::unique_ptr<Impl> impl;
};

#endif /* LIBBOBOL_LOD_TELEMETRY_PRIVATE_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8 cino=N-s
 */
