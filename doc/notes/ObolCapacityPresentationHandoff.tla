------------------- MODULE ObolCapacityPresentationHandoff -------------------
\* Ordering contract between one frozen occurrence allocation and the
\* renderer-capacity sample which measures it.
\*
\* The production allocator may change many occurrence-local cuts in one
\* bounded pass while an older renderer-wide ceiling still protects the last
\* coherent framebuffer.  The changed pass is valid application evidence: it
\* must be allowed to hand off to a ceiling-free successor frame.  A sample is
\* consumable only after every assigned cut is active and the ceiling is
\* absent or provably inert.  Frames observed before that point are not samples
\* and cannot spend the bounded sample allowance.
\*
\* A complete allocation may itself exceed its certified budget.  That state
\* belongs to the same handoff: a finite local-representation reduction either
\* makes it fit or publishes an explicit bounded constraint.  It must never
\* fall back into another identical generic allocation.
\*
\* A completed mechanical scan which selects no successor is not a semantic
\* population edge.  It preserves both the population revision and the
\* allocation-certificate revision so the current handoff can consume that
\* certificate.  Reopening the revision at this boundary would permit an
\* infinite apply/scan/handoff loop with no changed population.
\*
\* Certification may select one final occurrence allocation.  Applying that
\* cut and completing its generic mechanical pass do not invalidate the
\* semantic capacity problem.  The terminal certificate remains authoritative
\* through both actions; an unchanged completed pass cannot rearm it.
\*
\* This model proves ownership, safety, and liveness for a frozen allocation.
\* It deliberately does not model frame-time classification, aggregate cost,
\* or visual quality; those remain pure C++ and graphical tests.  It does
\* model the exact-empty boundary: zero renderer work is a valid sample only
\* after an authoritative visibility census proves that no occurrence should
\* have been presented.  It also models the resident-availability retry which
\* precedes this handoff: an allocation selected against an older resident
\* prefix may not consume a newer capacity revision as a no-op.  The retry
\* must either make the assigned cut drawable or publish a bounded availability
\* constraint.

EXTENDS Naturals, TLC

CONSTANTS SampleLimit, ReductionLimit, BudgetLevelCount

ASSUME SampleLimit > 0
ASSUME ReductionLimit > 0
ASSUME BudgetLevelCount > 1

Phases == {"apply", "handoff", "present_exact", "measure",
           "apply_certificate", "complete_certificate", "terminal"}
CeilingStates == {"effective", "inert", "none"}
BudgetStates == {"fits", "over"}
VisibilityStates == {"unknown", "empty", "nonempty"}
AvailabilityRetryStates == {"required", "requested", "resolved", "constrained"}

VARIABLES phase,
          cutsApplied,
          assignedCutsDrawable,
	  availabilityRetryState,
          ceilingState,
          allocationBudgetState,
          reductionAttempts,
          samplesRemaining,
          samplesConsumed,
          invalidFramesObserved,
          noOpScanAvailable,
          noOpScanCompleted,
          populationRevision,
	  allocationRevision,
	  ceilingFreeFrameRequested,
	  certificatePublished,
	  terminalBudgetGuard,
	  certificateApplied,
	  certificatePassCompleted,
	  constraintPublished,
	  reconciliationBudget,
	  activeBudget,
	  actualVisible,
	  visibilityCensus,
	  presentedWork

vars == <<phase, cutsApplied, assignedCutsDrawable, availabilityRetryState,
          ceilingState,
          allocationBudgetState,
          reductionAttempts, samplesRemaining, samplesConsumed,
          invalidFramesObserved, noOpScanAvailable, noOpScanCompleted,
	  populationRevision, allocationRevision,
	  ceilingFreeFrameRequested, certificatePublished,
	  terminalBudgetGuard,
	  certificateApplied, certificatePassCompleted,
	  constraintPublished, reconciliationBudget, activeBudget,
	  actualVisible, visibilityCensus, presentedWork>>

PresentationExact == cutsApplied /\ ceilingState # "effective"

ExactVisibilityPresentation ==
    /\ visibilityCensus # "unknown"
    /\ (presentedWork <=> visibilityCensus = "nonempty")

MeasurablePresentation == PresentationExact /\ ExactVisibilityPresentation

TypeOK ==
    /\ phase \in Phases
    /\ cutsApplied \in BOOLEAN
    /\ assignedCutsDrawable \in BOOLEAN
    /\ availabilityRetryState \in AvailabilityRetryStates
    /\ ceilingState \in CeilingStates
    /\ allocationBudgetState \in BudgetStates
    /\ reductionAttempts \in 0..ReductionLimit
    /\ samplesRemaining \in 0..SampleLimit
    /\ samplesConsumed \in 0..SampleLimit
    /\ invalidFramesObserved \in 0..2
    /\ noOpScanAvailable \in BOOLEAN
    /\ noOpScanCompleted \in BOOLEAN
    /\ populationRevision \in Nat
    /\ allocationRevision \in Nat
    /\ ceilingFreeFrameRequested \in BOOLEAN
    /\ certificatePublished \in BOOLEAN
    /\ terminalBudgetGuard \in BOOLEAN
    /\ certificateApplied \in BOOLEAN
    /\ certificatePassCompleted \in BOOLEAN
    /\ constraintPublished \in BOOLEAN
    /\ reconciliationBudget \in 0..BudgetLevelCount
    /\ activeBudget \in 0..BudgetLevelCount
    /\ actualVisible \in BOOLEAN
    /\ visibilityCensus \in VisibilityStates
    /\ presentedWork \in BOOLEAN

CurrentAllocationCertificate ==
    populationRevision = allocationRevision

SampleAccounting ==
    samplesRemaining + samplesConsumed = SampleLimit

AppliedCutsAreDrawable ==
    cutsApplied => assignedCutsDrawable

AvailabilityRetryIsBounded ==
    /\ (assignedCutsDrawable <=> availabilityRetryState = "resolved")
    /\ (availabilityRetryState \in {"required", "requested"} =>
           phase = "apply" /\ ~cutsApplied)
    /\ (availabilityRetryState = "constrained" =>
           /\ phase = "terminal"
           /\ ~cutsApplied
           /\ reconciliationBudget = 0
           /\ activeBudget = 0)

SampleRequiresExactPresentation ==
    samplesConsumed > 0 =>
        /\ MeasurablePresentation
        /\ allocationBudgetState = "fits"

VisibilityCensusIsAuthoritative ==
    /\ (visibilityCensus = "empty" => ~actualVisible)
    /\ (visibilityCensus = "nonempty" => actualVisible)

ZeroWorkSampleRequiresExactEmpty ==
    samplesConsumed > 0 /\ ~presentedWork => visibilityCensus = "empty"

TerminalCertificateIsExact ==
    certificatePublished =>
        /\ phase \in {"apply_certificate", "complete_certificate", "terminal"}
        /\ MeasurablePresentation
        /\ allocationBudgetState = "fits"
        /\ samplesRemaining = 0
        /\ samplesConsumed = SampleLimit

TerminalCertificateApplication ==
    /\ certificateApplied => certificatePublished
    /\ certificatePassCompleted => certificateApplied
    /\ (phase = "apply_certificate" => ~certificateApplied)
    /\ (phase = "complete_certificate" =>
           certificateApplied /\ ~certificatePassCompleted)
    /\ (phase = "terminal" /\ certificatePublished =>
           certificatePassCompleted)

TerminalCertificateGuardIsAtomic ==
    terminalBudgetGuard = certificatePublished

TerminalCertificateRetiresHandoffDebt ==
    certificatePublished =>
        /\ reconciliationBudget = 0
        /\ activeBudget = 0

CanonicalReconciliationBudget ==
    activeBudget = reconciliationBudget

TerminalConstraintIsBounded ==
    constraintPublished =>
        /\ phase = "terminal"
        /\ cutsApplied
        /\ allocationBudgetState = "over"
        /\ reductionAttempts = ReductionLimit

NoPrematureTerminal ==
    phase = "terminal" =>
        certificatePublished \/ constraintPublished \/
        availabilityRetryState = "constrained"

CeilingFreeFrameOwnership ==
    /\ (phase = "present_exact" =>
           /\ ceilingFreeFrameRequested
           /\ PresentationExact
           /\ CurrentAllocationCertificate)
    /\ (phase \in {"measure", "terminal"} => ~ceilingFreeFrameRequested)

Init ==
    /\ assignedCutsDrawable \in BOOLEAN
    /\ availabilityRetryState =
           IF assignedCutsDrawable THEN "resolved" ELSE "required"
    /\ cutsApplied \in BOOLEAN
    /\ cutsApplied => assignedCutsDrawable
    /\ ceilingState \in {"effective", "inert"}
    /\ allocationBudgetState \in BudgetStates
    /\ reductionAttempts = 0
    /\ phase = IF cutsApplied THEN "handoff" ELSE "apply"
    /\ samplesRemaining = SampleLimit
    /\ samplesConsumed = 0
    /\ invalidFramesObserved = 0
    /\ noOpScanAvailable \in BOOLEAN
    /\ noOpScanCompleted = FALSE
    /\ populationRevision = 1
    /\ allocationRevision = 1
    /\ ceilingFreeFrameRequested = FALSE
    /\ certificatePublished = FALSE
    /\ terminalBudgetGuard = FALSE
    /\ certificateApplied = FALSE
    /\ certificatePassCompleted = FALSE
    /\ constraintPublished = FALSE
    /\ reconciliationBudget = 0
    /\ activeBudget = 0
    /\ actualVisible \in BOOLEAN
    /\ visibilityCensus = "unknown"
    /\ presentedWork = FALSE

RequestAvailabilityRetry ==
    /\ phase = "apply"
    /\ ~cutsApplied
    /\ ~assignedCutsDrawable
    /\ availabilityRetryState = "required"
    /\ availabilityRetryState' = "requested"
    /\ UNCHANGED <<phase, cutsApplied, assignedCutsDrawable, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, ceilingFreeFrameRequested,
		    certificatePublished,
		    terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus, presentedWork>>

CompleteAvailabilityRetry ==
    /\ phase = "apply"
    /\ ~cutsApplied
    /\ ~assignedCutsDrawable
    /\ availabilityRetryState = "requested"
    /\ assignedCutsDrawable' = TRUE
    /\ availabilityRetryState' = "resolved"
    /\ UNCHANGED <<phase, cutsApplied, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, ceilingFreeFrameRequested,
		    certificatePublished,
		    terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus, presentedWork>>

PublishAvailabilityConstraint ==
    /\ phase = "apply"
    /\ ~cutsApplied
    /\ ~assignedCutsDrawable
    /\ availabilityRetryState = "requested"
    /\ phase' = "terminal"
    /\ availabilityRetryState' = "constrained"
    /\ reconciliationBudget' = 0
    /\ activeBudget' = 0
    /\ UNCHANGED <<cutsApplied, assignedCutsDrawable, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, ceilingFreeFrameRequested,
		    certificatePublished,
		    terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, actualVisible, visibilityCensus,
		    presentedWork>>

ApplyAssignedCuts ==
    /\ phase = "apply"
    /\ ~cutsApplied
    /\ assignedCutsDrawable
    /\ cutsApplied' = TRUE
    /\ phase' = "handoff"
    /\ UNCHANGED <<assignedCutsDrawable, availabilityRetryState, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, ceilingFreeFrameRequested,
		    certificatePublished,
		    terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus, presentedWork>>

CompleteNoOpAllocationScan ==
    /\ phase \in {"handoff", "measure"}
    /\ noOpScanAvailable
    /\ ~noOpScanCompleted
    /\ noOpScanCompleted' = TRUE
    /\ UNCHANGED <<phase, cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, noOpScanAvailable,
		    populationRevision, allocationRevision,
		    ceilingFreeFrameRequested, certificatePublished,
		    terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus, presentedWork>>

ReduceOverBudgetAllocation ==
    /\ phase = "handoff"
    /\ cutsApplied
    /\ allocationBudgetState = "over"
    /\ reductionAttempts < ReductionLimit
    /\ reductionAttempts' = reductionAttempts + 1
    /\ allocationBudgetState' \in BudgetStates
    /\ UNCHANGED <<phase, cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
                    samplesRemaining,
                    samplesConsumed, invalidFramesObserved,
                    noOpScanAvailable, noOpScanCompleted,
		    populationRevision, allocationRevision,
		    ceilingFreeFrameRequested, certificatePublished,
		    terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus, presentedWork>>

ReconcileCeiling ==
    /\ phase = "handoff"
    /\ cutsApplied
    /\ allocationBudgetState = "fits"
    /\ ceilingState' = IF ceilingState = "effective"
                           THEN "none"
                           ELSE ceilingState
    /\ phase' = "present_exact"
    /\ ceilingFreeFrameRequested' = TRUE
    /\ UNCHANGED <<cutsApplied, assignedCutsDrawable,
                    availabilityRetryState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, certificatePublished,
		    terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus, presentedWork>>

CertifyVisibilityCensus ==
    /\ phase = "present_exact"
    /\ visibilityCensus = "unknown"
    /\ visibilityCensus' = IF actualVisible THEN "nonempty" ELSE "empty"
    /\ UNCHANGED <<phase, cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, ceilingFreeFrameRequested,
		    certificatePublished, terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, presentedWork>>

PublishPresentedWork ==
    /\ phase = "present_exact"
    /\ visibilityCensus = "nonempty"
    /\ ~presentedWork
    /\ presentedWork' = TRUE
    /\ UNCHANGED <<phase, cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, ceilingFreeFrameRequested,
		    certificatePublished, terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus>>

PresentCeilingFreeCandidate ==
    /\ phase = "present_exact"
    /\ ceilingFreeFrameRequested
    /\ MeasurablePresentation
    /\ CurrentAllocationCertificate
    /\ phase' = "measure"
    /\ ceilingFreeFrameRequested' = FALSE
    /\ UNCHANGED <<cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, certificatePublished,
		    terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus, presentedWork>>

PublishBoundedConstraint ==
    /\ phase = "handoff"
    /\ cutsApplied
    /\ allocationBudgetState = "over"
    /\ reductionAttempts = ReductionLimit
    /\ phase' = "terminal"
    /\ constraintPublished' = TRUE
    /\ UNCHANGED <<cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
                    allocationBudgetState,
                    reductionAttempts, samplesRemaining, samplesConsumed,
                    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, ceilingFreeFrameRequested,
		    certificatePublished, certificateApplied,
		    terminalBudgetGuard,
		    certificatePassCompleted, reconciliationBudget,
		    activeBudget, actualVisible, visibilityCensus,
		    presentedWork>>

ObserveInvalidFrame ==
    /\ phase \in {"apply", "handoff"}
    /\ invalidFramesObserved < 2
    /\ invalidFramesObserved' = invalidFramesObserved + 1
    /\ UNCHANGED <<phase, cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
                    allocationBudgetState, reductionAttempts,
                    samplesRemaining, samplesConsumed,
                    noOpScanAvailable, noOpScanCompleted,
		    populationRevision, allocationRevision,
		    ceilingFreeFrameRequested, certificatePublished,
		    terminalBudgetGuard,
		    certificateApplied, certificatePassCompleted,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus, presentedWork>>

ConsumeSample ==
    /\ phase = "measure"
    /\ MeasurablePresentation
    /\ allocationBudgetState = "fits"
    /\ samplesRemaining > 0
    /\ samplesRemaining' = samplesRemaining - 1
    /\ samplesConsumed' = samplesConsumed + 1
    /\ IF samplesRemaining = 1
	  THEN /\ phase' = "apply_certificate"
	       /\ certificatePublished' = TRUE
	       /\ terminalBudgetGuard' = TRUE
	       /\ reconciliationBudget' = 0
	       /\ activeBudget' = 0
          ELSE /\ phase' = phase
	       /\ certificatePublished' = FALSE
	       /\ terminalBudgetGuard' = FALSE
	       /\ UNCHANGED <<reconciliationBudget, activeBudget>>
    /\ UNCHANGED <<cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
                    allocationBudgetState,
                    reductionAttempts, invalidFramesObserved,
		    noOpScanAvailable, noOpScanCompleted,
		    populationRevision, allocationRevision,
		    ceilingFreeFrameRequested, certificateApplied,
		    certificatePassCompleted, constraintPublished,
		    actualVisible,
		    visibilityCensus, presentedWork>>

ApplyCertifiedAllocation ==
    /\ phase = "apply_certificate"
    /\ certificatePublished
    /\ ~certificateApplied
    /\ phase' = "complete_certificate"
    /\ certificateApplied' = TRUE
    /\ UNCHANGED <<cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
		    allocationBudgetState, reductionAttempts,
		    samplesRemaining, samplesConsumed,
		    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, ceilingFreeFrameRequested,
		    certificatePublished, certificatePassCompleted,
		    terminalBudgetGuard,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus, presentedWork>>

CompleteCertifiedAllocationPass ==
    /\ phase = "complete_certificate"
    /\ certificatePublished
    /\ certificateApplied
    /\ ~certificatePassCompleted
    /\ phase' = "terminal"
    /\ certificatePassCompleted' = TRUE
    /\ UNCHANGED <<cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
		    allocationBudgetState, reductionAttempts,
		    samplesRemaining, samplesConsumed,
		    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, ceilingFreeFrameRequested,
		    certificatePublished, certificateApplied,
		    terminalBudgetGuard,
		    constraintPublished, reconciliationBudget, activeBudget,
		    actualVisible, visibilityCensus, presentedWork>>

\* Repeated callers share one canonical reconciliation request.  A tighter
\* request may reduce it; a weaker reassertion is a semantic no-op and cannot
\* overwrite the active scalar used by the allocation pass.
RequestReconciliation ==
    /\ phase # "terminal"
    /\ ~certificatePublished
    /\ \E budget \in 1..BudgetLevelCount :
         LET canonicalBudget ==
             IF reconciliationBudget = 0 \/ budget < reconciliationBudget
             THEN budget
             ELSE reconciliationBudget
         IN /\ reconciliationBudget' = canonicalBudget
            /\ activeBudget' = canonicalBudget
    /\ UNCHANGED <<phase, cutsApplied, assignedCutsDrawable,
                    availabilityRetryState, ceilingState,
		    allocationBudgetState, reductionAttempts,
		    samplesRemaining, samplesConsumed,
		    invalidFramesObserved, noOpScanAvailable,
		    noOpScanCompleted, populationRevision,
		    allocationRevision, ceilingFreeFrameRequested,
		    certificatePublished, certificateApplied,
		    terminalBudgetGuard,
		    certificatePassCompleted, constraintPublished,
		    actualVisible, visibilityCensus, presentedWork>>

Next ==
    \/ RequestAvailabilityRetry
    \/ CompleteAvailabilityRetry
    \/ PublishAvailabilityConstraint
    \/ ApplyAssignedCuts
    \/ CompleteNoOpAllocationScan
    \/ ReduceOverBudgetAllocation
    \/ ReconcileCeiling
    \/ CertifyVisibilityCensus
    \/ PublishPresentedWork
    \/ PresentCeilingFreeCandidate
    \/ PublishBoundedConstraint
    \/ ObserveInvalidFrame
    \/ ConsumeSample
    \/ ApplyCertifiedAllocation
    \/ CompleteCertifiedAllocationPass
    \/ RequestReconciliation

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(RequestAvailabilityRetry)
    /\ WF_vars(CompleteAvailabilityRetry)
    /\ WF_vars(PublishAvailabilityConstraint)
    /\ WF_vars(ApplyAssignedCuts)
    /\ WF_vars(CompleteNoOpAllocationScan)
    /\ WF_vars(ReduceOverBudgetAllocation)
    /\ WF_vars(ReconcileCeiling)
    /\ WF_vars(CertifyVisibilityCensus)
    /\ WF_vars(PublishPresentedWork)
    /\ WF_vars(PresentCeilingFreeCandidate)
    /\ WF_vars(PublishBoundedConstraint)
    /\ WF_vars(ConsumeSample)
    /\ WF_vars(ApplyCertifiedAllocation)
    /\ WF_vars(CompleteCertifiedAllocationPass)

EventuallyReady == <>(phase = "terminal" /\
    (certificatePassCompleted \/ constraintPublished \/
     availabilityRetryState = "constrained"))

\* Production refinement includes every ordinary planning pump, not only the
\* bounded search.  Once terminal evidence is current, no unmodeled planner
\* may make the same semantic problem active again.
TerminalDoesNotReopen ==
    [](phase = "terminal" => [](phase = "terminal"))

=============================================================================
