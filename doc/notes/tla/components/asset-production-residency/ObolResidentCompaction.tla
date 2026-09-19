----------------------- MODULE ObolResidentCompaction -----------------------
\* Refinement of the canonical progressive-pipeline rule that background
\* resident-memory work is not presentation authority.  A compaction target
\* is certified by one complete consumer-demand snapshot and semantic epoch.
\* Demand invalidation may race every worker phase; stale work cannot commit.

EXTENDS Naturals, TLC

CONSTANT MaxDemandRevision, MaxSemanticEpoch

Quality == 0..2
WorkStates == {"idle", "queued", "preparing", "complete"}

VARIABLES demandRevision, semanticEpoch,
          snapshotDemandRevision, snapshotSemanticEpoch,
          workState, workDemandRevision, workSemanticEpoch,
          targetQuality, resultDemandRevision,
          residentQuality, committedQuality

vars == <<demandRevision, semanticEpoch,
          snapshotDemandRevision, snapshotSemanticEpoch,
          workState, workDemandRevision, workSemanticEpoch,
          targetQuality, resultDemandRevision,
          residentQuality, committedQuality>>

SnapshotCurrent ==
    snapshotDemandRevision = demandRevision /\
    snapshotSemanticEpoch = semanticEpoch

WorkCurrent ==
    workDemandRevision = demandRevision /\
    workSemanticEpoch = semanticEpoch /\ SnapshotCurrent

Init ==
    /\ demandRevision = 1
    /\ semanticEpoch = 1
    /\ snapshotDemandRevision = 0
    /\ snapshotSemanticEpoch = 0
    /\ workState = "idle"
    /\ workDemandRevision = 0
    /\ workSemanticEpoch = 0
    /\ targetQuality = 0
    /\ resultDemandRevision = 0
    /\ residentQuality = 2
    /\ committedQuality = 1

InstallSnapshot(target) ==
    /\ workState = "idle"
    /\ target \in committedQuality..residentQuality
    /\ snapshotDemandRevision' = demandRevision
    /\ snapshotSemanticEpoch' = semanticEpoch
    /\ workState' = "queued"
    /\ workDemandRevision' = demandRevision
    /\ workSemanticEpoch' = semanticEpoch
    /\ targetQuality' = target
    /\ resultDemandRevision' = 0
    /\ UNCHANGED <<demandRevision, semanticEpoch, residentQuality,
                    committedQuality>>

StartWork ==
    /\ workState = "queued"
    /\ workState' = "preparing"
    /\ UNCHANGED <<demandRevision, semanticEpoch,
                    snapshotDemandRevision, snapshotSemanticEpoch,
                    workDemandRevision, workSemanticEpoch, targetQuality,
                    resultDemandRevision, residentQuality, committedQuality>>

FinishCurrentWork ==
    /\ workState = "preparing"
    /\ WorkCurrent
    /\ workState' = "complete"
    /\ resultDemandRevision' = workDemandRevision
    /\ residentQuality' = targetQuality
    /\ UNCHANGED <<demandRevision, semanticEpoch,
                    snapshotDemandRevision, snapshotSemanticEpoch,
                    workDemandRevision, workSemanticEpoch, targetQuality,
                    committedQuality>>

DiscardStaleWork ==
    /\ workState \in {"queued", "preparing"}
    /\ ~WorkCurrent
    /\ workState' = "idle"
    /\ resultDemandRevision' = 0
    /\ UNCHANGED <<demandRevision, semanticEpoch,
                    snapshotDemandRevision, snapshotSemanticEpoch,
                    workDemandRevision, workSemanticEpoch, targetQuality,
                    residentQuality, committedQuality>>

PublishCurrentResult ==
    /\ workState = "complete"
    /\ resultDemandRevision = demandRevision
    /\ workSemanticEpoch = semanticEpoch
    /\ targetQuality >= committedQuality
    /\ committedQuality' = targetQuality
    /\ workState' = "idle"
    /\ resultDemandRevision' = 0
    /\ UNCHANGED <<demandRevision, semanticEpoch,
                    snapshotDemandRevision, snapshotSemanticEpoch,
                    workDemandRevision, workSemanticEpoch, targetQuality,
                    residentQuality>>

DiscardStaleResult ==
    /\ workState = "complete"
    /\ \/ resultDemandRevision # demandRevision
       \/ workSemanticEpoch # semanticEpoch
    /\ workState' = "idle"
    /\ resultDemandRevision' = 0
    /\ UNCHANGED <<demandRevision, semanticEpoch,
                    snapshotDemandRevision, snapshotSemanticEpoch,
                    workDemandRevision, workSemanticEpoch, targetQuality,
                    residentQuality, committedQuality>>

InvalidateSemanticInput ==
    /\ semanticEpoch < MaxSemanticEpoch
    /\ semanticEpoch' = semanticEpoch + 1
    /\ snapshotDemandRevision' = 0
    /\ snapshotSemanticEpoch' = 0
    /\ UNCHANGED <<demandRevision, workState, workDemandRevision,
                    workSemanticEpoch, targetQuality, resultDemandRevision,
                    residentQuality, committedQuality>>

ChangeDemand ==
    /\ demandRevision < MaxDemandRevision
    /\ semanticEpoch < MaxSemanticEpoch
    /\ demandRevision' = demandRevision + 1
    /\ semanticEpoch' = semanticEpoch + 1
    /\ snapshotDemandRevision' = 0
    /\ snapshotSemanticEpoch' = 0
    /\ UNCHANGED <<workState, workDemandRevision, workSemanticEpoch,
                    targetQuality, resultDemandRevision, residentQuality,
                    committedQuality>>

Next ==
    \/ \E target \in Quality: InstallSnapshot(target)
    \/ StartWork
    \/ FinishCurrentWork
    \/ DiscardStaleWork
    \/ PublishCurrentResult
    \/ DiscardStaleResult
    \/ InvalidateSemanticInput
    \/ ChangeDemand

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(StartWork)
    /\ WF_vars(FinishCurrentWork)
    /\ WF_vars(DiscardStaleWork)
    /\ WF_vars(PublishCurrentResult)
    /\ WF_vars(DiscardStaleResult)

TypeOK ==
    /\ demandRevision \in 1..MaxDemandRevision
    /\ semanticEpoch \in 1..MaxSemanticEpoch
    /\ snapshotDemandRevision \in 0..MaxDemandRevision
    /\ snapshotSemanticEpoch \in 0..MaxSemanticEpoch
    /\ workState \in WorkStates
    /\ workDemandRevision \in 0..MaxDemandRevision
    /\ workSemanticEpoch \in 0..MaxSemanticEpoch
    /\ targetQuality \in Quality
    /\ resultDemandRevision \in 0..MaxDemandRevision
    /\ residentQuality \in Quality
    /\ committedQuality \in Quality

StaleWorkCannotCommit ==
    workState = "preparing" /\ ~WorkCurrent => ~ENABLED FinishCurrentWork

CompletedResultIsCertified ==
    workState = "complete" => resultDemandRevision = workDemandRevision

CompactionDoesNotRegressCommittedQuality == committedQuality >= 1

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => workState = "idle"

EventuallyRetireFinalWork ==
    []((demandRevision = MaxDemandRevision /\
        semanticEpoch = MaxSemanticEpoch /\ workState # "idle") =>
       <>(workState = "idle"))

=============================================================================
