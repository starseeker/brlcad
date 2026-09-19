-------------------------- MODULE ObolSpatialProducer -------------------------
\* Bounded control-plane contract for serialized spatial PoP production.
\*
\* A real leaf carries many faces and its data-plane work is deliberately not
\* modeled here.  TLC checks the delicate lifecycle boundary instead: a
\* transient source bootstrap or a first cache page may be presented directly
\* while a cache write is active; durable cache discovery is atomic; and a
\* capacity-constrained live page is useful only to its generating process.

EXTENDS Naturals, Sequences, TLC

CONSTANT LeafCount

Leaves == 1..LeafCount
LeafStates == {"unseen", "durable", "live"}
Phases == {"classifying", "complete", "constrained", "cancelled"}
PreviewStates == {"none", "present", "retired"}

VARIABLES phase,
          nextLeaf,
          leafState,
          writerLock,
          descriptorPublished,
          cacheReadInFlight,
          previewState,
          bootstrapPublished

vars == <<phase, nextLeaf, leafState, writerLock, descriptorPublished,
          cacheReadInFlight, previewState, bootstrapPublished>>

TypeOK == /\ phase \in Phases
          /\ nextLeaf \in 1..(LeafCount + 1)
          /\ leafState \in [Leaves -> LeafStates]
          /\ writerLock \in BOOLEAN
          /\ descriptorPublished \in BOOLEAN
          /\ cacheReadInFlight \in BOOLEAN
          /\ previewState \in PreviewStates
          /\ bootstrapPublished \in BOOLEAN

Init == /\ phase = "classifying"
        /\ nextLeaf = 1
        /\ leafState = [leaf \in Leaves |-> "unseen"]
        /\ writerLock = TRUE
        /\ descriptorPublished = FALSE
        /\ cacheReadInFlight = FALSE
        /\ previewState = "none"
        /\ bootstrapPublished = FALSE

\* A mapped-source bootstrap is bounded and process-local.  It is allowed
\* before the global source identity and first durable cache page exist, but
\* cannot itself make a cache record discoverable or satisfy completion.
PublishBootstrapPreview == /\ phase = "classifying"
                           /\ previewState = "none"
                           /\ ~bootstrapPublished
                           /\ cacheReadInFlight = FALSE
                           /\ previewState' = "present"
                           /\ bootstrapPublished' = TRUE
                           /\ UNCHANGED <<phase, nextLeaf, leafState,
                                           writerLock, descriptorPublished,
                                           cacheReadInFlight>>

\* Each page write either commits one immutable durable record or discovers a
\* capacity edge.  The latter may retain just that page in the producer, but
\* must stop classification before a partial descriptor can be published.
BuildLeaf == /\ phase = "classifying"
             /\ nextLeaf \in Leaves
             /\ \/ /\ leafState' = [leafState EXCEPT ![nextLeaf] = "durable"]
                   /\ nextLeaf' = nextLeaf + 1
                   /\ UNCHANGED <<phase, writerLock, descriptorPublished,
                                   cacheReadInFlight, previewState,
                                   bootstrapPublished>>
                \/ /\ leafState' = [leafState EXCEPT ![nextLeaf] = "live"]
                   /\ phase' = "constrained"
                   /\ writerLock' = FALSE
                   /\ nextLeaf' = nextLeaf + 1
                   /\ UNCHANGED <<descriptorPublished, cacheReadInFlight,
                                   previewState, bootstrapPublished>>

\* The first page is handed directly to the renderer.  It must never reopen
\* its just-written cache record: the producer owns the write transaction and
\* a cache read at this point is a self-deadlock risk.
PublishFirstPreview == /\ previewState = "none"
                       /\ leafState[1] # "unseen"
                       /\ cacheReadInFlight = FALSE
                       /\ previewState' = "present"
                       /\ UNCHANGED <<phase, nextLeaf, leafState, writerLock,
                                       descriptorPublished, cacheReadInFlight,
                                       bootstrapPublished>>

PublishDescriptor == /\ phase = "classifying"
                     /\ nextLeaf = LeafCount + 1
                     /\ \A leaf \in Leaves : leafState[leaf] = "durable"
                     /\ phase' = "complete"
                     /\ writerLock' = FALSE
                     /\ descriptorPublished' = TRUE
                     /\ previewState' = "retired"
                     /\ UNCHANGED <<nextLeaf, leafState, cacheReadInFlight,
                                     bootstrapPublished>>

Cancel == /\ phase = "classifying"
          /\ phase' = "cancelled"
          /\ writerLock' = FALSE
          /\ previewState' = IF previewState = "present" THEN "retired"
                             ELSE previewState
          /\ UNCHANGED <<nextLeaf, leafState, descriptorPublished,
                          cacheReadInFlight, bootstrapPublished>>

\* Only a completed descriptor may be discovered by a later consumer.  Live
\* records are intentionally process-local and are not a cache namespace.
BeginCacheRead == /\ descriptorPublished
                  /\ ~writerLock
                  /\ ~cacheReadInFlight
                  /\ cacheReadInFlight' = TRUE
                  /\ UNCHANGED <<phase, nextLeaf, leafState, writerLock,
                                  descriptorPublished, previewState,
                                  bootstrapPublished>>

EndCacheRead == /\ cacheReadInFlight
                /\ cacheReadInFlight' = FALSE
                /\ UNCHANGED <<phase, nextLeaf, leafState, writerLock,
                                descriptorPublished, previewState,
                                bootstrapPublished>>

Next == \/ PublishBootstrapPreview
        \/ BuildLeaf
        \/ PublishFirstPreview
        \/ PublishDescriptor
        \/ Cancel
        \/ BeginCacheRead
        \/ EndCacheRead

Spec == Init /\ [][Next]_vars
        /\ WF_vars(BuildLeaf)
        /\ WF_vars(PublishDescriptor)

AtomicDescriptor ==
    descriptorPublished =>
        /\ phase = "complete"
        /\ \A leaf \in Leaves : leafState[leaf] = "durable"

NoPartialCacheDiscovery == ~descriptorPublished => ~cacheReadInFlight

NoCacheReentryWhileWriting == ~(writerLock /\ cacheReadInFlight)

PreviewHasMaterializedFirstPage ==
    previewState = "present" =>
        bootstrapPublished \/ leafState[1] # "unseen"

ConstrainedPageIsProcessLocal ==
    (\E leaf \in Leaves : leafState[leaf] = "live") =>
        /\ phase = "constrained"
        /\ ~descriptorPublished
        /\ ~writerLock

TerminalHasNoWriter ==
    phase \in {"complete", "constrained", "cancelled"} => ~writerLock

Terminal == phase \in {"complete", "constrained", "cancelled"}
DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyTerminal == <>(phase \in {"complete", "constrained", "cancelled"})

=============================================================================
