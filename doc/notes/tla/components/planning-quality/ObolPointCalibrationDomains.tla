-------------------- MODULE ObolPointCalibrationDomains --------------------
\* Point aggregation has two distinct timing contracts.  Responsive evidence
\* targets the ordinary quiet redraw cadence.  Static evidence targets the
\* longer, interruptible endpoint deadline.  A frame can be unsafe for the
\* former and safe for the latter, so neither evidence nor pressure actions
\* may cross that boundary.
\*
\* Once the finite source coverage epoch closes, terminal refinement owns a
\* bounded sequence of threshold/preload steps.  Discovery-time pressure is
\* disabled throughout that sequence; otherwise it can undo a static step and
\* create the observed responsive -> static -> responsive cycle.

EXTENDS Naturals

CONSTANTS MaxResponsiveSteps, MaxStaticSteps

Owners == {"responsive", "static", "preload", "idle"}

VARIABLES owner, deadlineClass, coverageComplete,
          responsiveEvidenceRevision, staticEvidenceRevision,
          staticStepsRemaining, terminal

vars == <<owner, deadlineClass, coverageComplete, responsiveEvidenceRevision,
          staticEvidenceRevision, staticStepsRemaining, terminal>>

TypeOK ==
    /\ owner \in Owners
    /\ deadlineClass \in {"responsive", "static", "none"}
    /\ coverageComplete \in BOOLEAN
    /\ responsiveEvidenceRevision \in 0..MaxResponsiveSteps
    /\ staticEvidenceRevision \in Nat
    /\ staticStepsRemaining \in 0..MaxStaticSteps
    /\ terminal \in BOOLEAN

Init ==
    /\ owner = "responsive"
    /\ deadlineClass = "responsive"
    /\ coverageComplete = FALSE
    /\ responsiveEvidenceRevision = 0
    /\ staticEvidenceRevision = 0
    /\ staticStepsRemaining = MaxStaticSteps
    /\ terminal = FALSE

ResponsiveObservation ==
    /\ ~coverageComplete
    /\ owner = "responsive"
    /\ ~terminal
    /\ responsiveEvidenceRevision < MaxResponsiveSteps
    /\ responsiveEvidenceRevision' = responsiveEvidenceRevision + 1
    /\ UNCHANGED <<owner, deadlineClass, coverageComplete,
                    staticEvidenceRevision,
                    staticStepsRemaining, terminal>>

CloseCoverage ==
    /\ ~coverageComplete
    /\ owner = "responsive"
    /\ coverageComplete' = TRUE
    /\ owner' = "static"
    /\ deadlineClass' = "static"
    /\ UNCHANGED <<responsiveEvidenceRevision, staticEvidenceRevision,
                    staticStepsRemaining, terminal>>

BeginStaticPreload ==
    /\ coverageComplete
    /\ owner = "static"
    /\ staticStepsRemaining > 0
    /\ ~terminal
    /\ owner' = "preload"
    /\ staticEvidenceRevision' = staticEvidenceRevision + 1
    /\ UNCHANGED <<deadlineClass, coverageComplete,
                    responsiveEvidenceRevision,
                    staticStepsRemaining, terminal>>

CompleteStaticPreload ==
    /\ coverageComplete
    /\ owner = "preload"
    /\ staticStepsRemaining > 0
    /\ ~terminal
    /\ owner' = "static"
    /\ staticStepsRemaining' = staticStepsRemaining - 1
    /\ UNCHANGED <<deadlineClass, coverageComplete,
                    responsiveEvidenceRevision,
                    staticEvidenceRevision, terminal>>

ReachStaticEndpoint ==
    /\ coverageComplete
    /\ owner = "static"
    /\ staticStepsRemaining = 0
    /\ ~terminal
    /\ owner' = "idle"
    /\ deadlineClass' = "none"
    /\ terminal' = TRUE
    /\ UNCHANGED <<coverageComplete, responsiveEvidenceRevision,
                    staticEvidenceRevision, staticStepsRemaining>>

TerminalIdle ==
    /\ terminal
    /\ owner = "idle"
    /\ UNCHANGED vars

Next ==
    \/ ResponsiveObservation
    \/ CloseCoverage
    \/ BeginStaticPreload
    \/ CompleteStaticPreload
    \/ ReachStaticEndpoint
    \/ TerminalIdle

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(BeginStaticPreload)
    /\ WF_vars(CompleteStaticPreload)
    /\ WF_vars(ReachStaticEndpoint)

\* Responsive pressure has no enabled transition after coverage closes.
NoResponsiveOwnerAfterCoverage ==
    coverageComplete /\ ~terminal => owner \in {"static", "preload"}

\* Static observations never overwrite or consume responsive evidence.
EvidenceDomainsIndependent ==
    staticEvidenceRevision <= MaxStaticSteps

\* Consuming a short-lived calibration latch cannot change the timing
\* contract of its evidence domain.  Static/preload work always retains the
\* event-driven endpoint deadline; responsive work always retains the
\* ordinary redraw deadline.
DeadlineBoundToEvidenceDomain ==
    /\ (owner = "responsive" => deadlineClass = "responsive")
    /\ (owner \in {"static", "preload"} => deadlineClass = "static")
    /\ (owner = "idle" => deadlineClass = "none")

EventuallyStaticEndpointAfterCoverage ==
    [](coverageComplete => <>terminal)

=============================================================================
