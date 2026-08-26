-------------------- MODULE ObolLiveSpatialPublication --------------------
\* Bounded ownership model for live spatial-page publication from a cold
\* serialized mesh.  It deliberately models publication state, not mesh
\* topology, page bytes, GPU timing, or visual-error measurement.
\*
\* The model captures the implementation boundary described in
\* live_spatial_publication.md.  A certified coverage part is the standing
\* whole-object presentation.  Completed private pages are additive layers;
\* they cannot replace coverage individually.  Only a complete final
\* hierarchy may retire coverage and publish the durable cache marker.

EXTENDS FiniteSets, TLC

Pages == {"first", "second", "third"}

VARIABLES partitioned, readyPages, publishedPages, coverage, finalGeometry,
          cacheMarked, cancelled, constrained, framePending,
          cancelledReadyPages, cancelledPublishedPages, cancelledFinalGeometry

vars == <<partitioned, readyPages, publishedPages, coverage, finalGeometry,
          cacheMarked, cancelled, constrained, framePending,
          cancelledReadyPages, cancelledPublishedPages, cancelledFinalGeometry>>

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

Partition ==
    /\ ~partitioned /\ ~cancelled /\ ~constrained
    /\ partitioned' = TRUE
    /\ UNCHANGED <<readyPages, publishedPages, coverage, finalGeometry,
                    cacheMarked, cancelled, constrained, framePending,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry>>

CompletePage(page) ==
    /\ page \in Pages
    /\ partitioned /\ ~cancelled /\ ~constrained
    /\ page \notin readyPages
    /\ readyPages' = readyPages \cup {page}
    /\ UNCHANGED <<partitioned, publishedPages, coverage, finalGeometry,
                    cacheMarked, cancelled, constrained, framePending,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry>>

PublishPage(page) ==
    /\ page \in readyPages \ publishedPages
    /\ ~cancelled /\ ~constrained
    /\ publishedPages' = publishedPages \cup {page}
    /\ framePending' = TRUE
    /\ UNCHANGED <<partitioned, readyPages, coverage, finalGeometry,
                    cacheMarked, cancelled, constrained,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry>>

Finalize ==
    /\ partitioned /\ readyPages = Pages /\ publishedPages = Pages
    /\ ~cancelled /\ ~constrained
    /\ finalGeometry' = TRUE
    /\ coverage' = FALSE
    /\ framePending' = TRUE
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, cacheMarked,
                    cancelled, constrained, cancelledReadyPages,
                    cancelledPublishedPages, cancelledFinalGeometry>>

MarkCache ==
    /\ finalGeometry /\ ~cancelled /\ ~cacheMarked
    /\ cacheMarked' = TRUE
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, coverage,
                    finalGeometry, cancelled, constrained, framePending,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry>>

Constrain ==
    /\ ~finalGeometry /\ ~cancelled /\ ~constrained
    /\ constrained' = TRUE
    /\ framePending' = TRUE
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, coverage,
                    finalGeometry, cacheMarked, cancelled,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry>>

Cancel ==
    /\ ~cacheMarked /\ ~cancelled
    /\ cancelled' = TRUE
    /\ cancelledReadyPages' = readyPages
    /\ cancelledPublishedPages' = publishedPages
    /\ cancelledFinalGeometry' = finalGeometry
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, coverage,
                    finalGeometry, cacheMarked, constrained, framePending>>

CompleteFrame ==
    /\ framePending
    /\ framePending' = FALSE
    /\ UNCHANGED <<partitioned, readyPages, publishedPages, coverage,
                    finalGeometry, cacheMarked, cancelled, constrained,
                    cancelledReadyPages, cancelledPublishedPages,
                    cancelledFinalGeometry>>

Next ==
    \/ Partition
    \/ \E page \in Pages: CompletePage(page)
    \/ \E page \in Pages: PublishPage(page)
    \/ Finalize
    \/ MarkCache
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

FinalRequiresCompletePublishedPageSet ==
    finalGeometry => publishedPages = Pages

CacheMarkerRequiresFinalHierarchy ==
    cacheMarked => finalGeometry /\ publishedPages = Pages

CancellationNeverPublishesCacheMarker == cancelled => ~cacheMarked

\* A local completed page is richer information but does not replace the
\* global preview.  The final hierarchy is the only state allowed to retire
\* coverage, so an incomplete source can never collapse to a local fragment.
IncompletePublicationRetainsCoverage ==
    ~finalGeometry => coverage

AlwaysDisplayable == coverage \/ finalGeometry \/ publishedPages # {}

Stable ==
    /\ ~framePending
    /\ (cancelled \/ constrained \/ cacheMarked)

EventuallyStable == <>Stable

=============================================================================
