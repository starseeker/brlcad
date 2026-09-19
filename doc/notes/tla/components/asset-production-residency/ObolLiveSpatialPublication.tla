-------------------- MODULE ObolLiveSpatialPublication --------------------
\* Bounded ownership model for live spatial-page publication from a cold
\* serialized mesh.  It deliberately models publication state, not mesh
\* topology, page bytes, GPU timing, or visual-error measurement.
\*
\* The model captures the implementation boundary described in
\* live_spatial_publication.md.  A certified coverage part is the standing
\* whole-object presentation.  Completed private pages are additive layers;
\* they cannot replace coverage individually.  Presentation is not a cache
\* construction prerequisite: only a complete final hierarchy may retire
\* coverage and publish the durable cache marker, even when some completed
\* live pages were never independently presented.  Durable replacement also
\* retains an already-published name mapping until the complete candidate's
\* mapping commits; a denied commit leaves that predecessor discoverable.

EXTENDS FiniteSets, TLC

Pages == {"first", "second", "third"}
CacheMappings == {"none", "prior", "candidate"}

VARIABLES partitioned, readyPages, publishedPages, coverage, finalGeometry,
          cacheMarked, cancelled, constrained, framePending,
          cancelledReadyPages, cancelledPublishedPages, cancelledFinalGeometry,
          cacheMapping, baselineCacheMapping, cacheUnavailable

vars == <<partitioned, readyPages, publishedPages, coverage, finalGeometry,
          cacheMarked, cancelled, constrained, framePending,
          cancelledReadyPages, cancelledPublishedPages, cancelledFinalGeometry,
          cacheMapping, baselineCacheMapping, cacheUnavailable>>

TypeOK ==
    /\ partitioned \in BOOLEAN
    /\ readyPages \subseteq Pages
    /\ publishedPages \subseteq Pages
    /\ coverage \in BOOLEAN
    /\ finalGeometry \in BOOLEAN
    /\ cacheMarked \in BOOLEAN
    /\ cancelled \in BOOLEAN
    /\ constrained \in BOOLEAN
    /\ framePending \in BOOLEAN
    /\ cancelledReadyPages \subseteq Pages
    /\ cancelledPublishedPages \subseteq Pages
    /\ cancelledFinalGeometry \in BOOLEAN
    /\ cacheMapping \in CacheMappings
    /\ baselineCacheMapping \in {"none", "prior"}
    /\ cacheUnavailable \in BOOLEAN

Init ==
    /\ partitioned = FALSE
    /\ readyPages = {}
    /\ publishedPages = {}
    /\ coverage = TRUE
    /\ finalGeometry = FALSE
    /\ cacheMarked = FALSE
    /\ cancelled = FALSE
    /\ constrained = FALSE
    /\ framePending = FALSE
    /\ cancelledReadyPages = {}
    /\ cancelledPublishedPages = {}
    /\ cancelledFinalGeometry = FALSE
    /\ cacheMapping \in {"none", "prior"}
    /\ baselineCacheMapping = cacheMapping
    /\ cacheUnavailable = FALSE

Partition ==
    /\ ~partitioned /\ ~cancelled /\ ~constrained
    /\ partitioned' = TRUE
    /\ UNCHANGED <<readyPages, publishedPages, coverage, finalGeometry,
                    cacheMarked, cancelled, constrained, framePending,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry, cacheMapping,
                    baselineCacheMapping, cacheUnavailable>>

CompletePage(page) ==
    /\ page \in Pages
    /\ partitioned /\ ~cancelled /\ ~constrained
    /\ page \notin readyPages
    /\ readyPages' = readyPages \cup {page}
    /\ UNCHANGED <<partitioned, publishedPages, coverage, finalGeometry,
                    cacheMarked, cancelled, constrained, framePending,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry, cacheMapping,
                    baselineCacheMapping, cacheUnavailable>>

PublishPage(page) ==
    /\ page \in readyPages \ publishedPages
    /\ ~cancelled /\ ~constrained
    /\ publishedPages' = publishedPages \cup {page}
    /\ framePending' = TRUE
    /\ UNCHANGED <<partitioned, readyPages, coverage, finalGeometry,
                    cacheMarked, cancelled, constrained,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry, cacheMapping,
                    baselineCacheMapping, cacheUnavailable>>

Finalize ==
    /\ partitioned /\ readyPages = Pages
    /\ ~cancelled /\ ~constrained
    /\ finalGeometry' = TRUE
    /\ coverage' = FALSE
    /\ framePending' = TRUE
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, cacheMarked,
                    cancelled, constrained, cancelledReadyPages,
                    cancelledPublishedPages, cancelledFinalGeometry,
                    cacheMapping, baselineCacheMapping, cacheUnavailable>>

MarkCache ==
    /\ finalGeometry /\ ~cancelled /\ ~cacheMarked /\ ~cacheUnavailable
    /\ cacheMarked' = TRUE
    /\ cacheMapping' = "candidate"
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, coverage,
                    finalGeometry, cancelled, constrained, framePending,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry, baselineCacheMapping,
                    cacheUnavailable>>

DenyCacheCommit ==
    /\ finalGeometry /\ ~cancelled /\ ~cacheMarked /\ ~cacheUnavailable
    /\ cacheUnavailable' = TRUE
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, coverage,
                    finalGeometry, cacheMarked, cancelled, constrained,
                    framePending, cancelledReadyPages,
                    cancelledPublishedPages, cancelledFinalGeometry,
                    cacheMapping, baselineCacheMapping>>

Constrain ==
    /\ ~finalGeometry /\ ~cancelled /\ ~constrained
    /\ constrained' = TRUE
    /\ framePending' = TRUE
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, coverage,
                    finalGeometry, cacheMarked, cancelled,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry, cacheMapping,
                    baselineCacheMapping, cacheUnavailable>>

Cancel ==
    /\ ~cacheMarked /\ ~cancelled
    /\ cancelled' = TRUE
    /\ cancelledReadyPages' = readyPages
    /\ cancelledPublishedPages' = publishedPages
    /\ cancelledFinalGeometry' = finalGeometry
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, coverage,
                    finalGeometry, cacheMarked, constrained, framePending>>
    /\ UNCHANGED <<cacheMapping, baselineCacheMapping, cacheUnavailable>>

CompleteFrame ==
    /\ framePending
    /\ framePending' = FALSE
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, coverage,
                    finalGeometry, cacheMarked, cancelled, constrained,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry, cacheMapping,
                    baselineCacheMapping, cacheUnavailable>>

Next ==
    \/ Partition
    \/ \E page \in Pages: CompletePage(page)
    \/ \E page \in Pages: PublishPage(page)
    \/ Finalize
    \/ MarkCache
    \/ DenyCacheCommit
    \/ Constrain
    \/ Cancel
    \/ CompleteFrame

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(Partition)
    /\ \A page \in Pages: WF_vars(CompletePage(page))
    /\ \A page \in Pages: WF_vars(PublishPage(page))
    /\ WF_vars(Finalize)
    /\ WF_vars(MarkCache)
    /\ WF_vars(CompleteFrame)

PagesAreCompletedBeforePublication == publishedPages \subseteq readyPages

CancellationFreezesProducerPublication ==
    cancelled => /\ readyPages = cancelledReadyPages
                 /\ publishedPages = cancelledPublishedPages
                 /\ finalGeometry = cancelledFinalGeometry

FinalRequiresCompletePageSet ==
    finalGeometry => readyPages = Pages

CacheMarkerRequiresFinalHierarchy ==
    cacheMarked =>
        /\ finalGeometry /\ readyPages = Pages
        /\ cacheMapping = "candidate"
        /\ ~cacheUnavailable

UnpublishedCacheMappingIsStable ==
    ~cacheMarked => cacheMapping = baselineCacheMapping

DeniedCacheCommitPreservesMapping ==
    cacheUnavailable =>
        /\ ~cacheMarked
        /\ cacheMapping = baselineCacheMapping

CancellationNeverPublishesCacheMarker == cancelled => ~cacheMarked

\* A local completed page is richer information but does not replace the
\* global preview.  The final hierarchy is the only state allowed to retire
\* coverage, so an incomplete source can never collapse to a local fragment.
IncompletePublicationRetainsCoverage ==
    ~finalGeometry => coverage

AlwaysDisplayable == coverage \/ finalGeometry \/ publishedPages # {}

Stable ==
    /\ ~framePending
    /\ (cancelled \/ constrained \/ cacheMarked \/ cacheUnavailable)

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Stable

EventuallyStable == <>Stable

=============================================================================
