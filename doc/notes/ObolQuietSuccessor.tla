------------------------- MODULE ObolQuietSuccessor -------------------------
\* Schedule-independent quiet successor below the canonical progressive
\* pipeline.
\*
\* The two transient controls represent different intermediate publication
\* schedules for identical semantic input.  They may require different
\* temporary renderer handoffs, but cannot select different quiet targets.
\* A revision-matched certificate has fixed precedence.  Without one, the
\* stable one-pixel demand is fixed and any motion ceiling is only a bounded
\* presentation-safety obligation.

EXTENDS Naturals, TLC

Phases == {"interacting", "restored", "planning", "terminal"}
TransientControls == {"rich", "coarseA", "coarseB"}
Sources == {"stable", "prior", "proven", "exact"}
Targets == {"stableQuality", "priorQuality", "provenQuality", "exactQuality"}
Handoffs == {"none", "allocation", "presentation"}

SelectedSource(exact, retained, prior, proven) ==
    IF exact THEN "exact"
    ELSE IF ~retained THEN "stable"
    ELSE IF prior THEN "prior"
    ELSE IF proven THEN "proven"
    ELSE "stable"

SelectedTarget(source) ==
    IF source = "exact" THEN "exactQuality"
    ELSE IF source = "prior" THEN "priorQuality"
    ELSE IF source = "proven" THEN "provenQuality"
    ELSE "stableQuality"

SelectedHandoff(source, retained, transient) ==
    IF source = "prior" THEN "presentation"
    ELSE IF retained /\ source \in {"stable", "proven"} /\
            transient # "rich" THEN "allocation"
    ELSE "none"

VARIABLES phase,
          retained,
          priorValid,
          provenValid,
          exactValid,
          transientA,
          transientB,
          sourceA,
          sourceB,
          targetA,
          targetB,
          handoff,
          planCount,
          inputClosed

vars == <<phase, retained, priorValid, provenValid, exactValid,
          transientA, transientB, sourceA, sourceB, targetA, targetB,
          handoff, planCount, inputClosed>>

TypeOK ==
    /\ phase \in Phases
    /\ retained \in BOOLEAN
    /\ priorValid \in BOOLEAN
    /\ provenValid \in BOOLEAN
    /\ exactValid \in BOOLEAN
    /\ transientA \in TransientControls
    /\ transientB \in TransientControls
    /\ sourceA \in Sources \union {"none"}
    /\ sourceB \in Sources \union {"none"}
    /\ targetA \in Targets \union {"none"}
    /\ targetB \in Targets \union {"none"}
    /\ handoff \in Handoffs
    /\ planCount \in 0..1
    /\ inputClosed \in BOOLEAN

Init ==
    /\ phase = "interacting"
    /\ retained \in BOOLEAN
    /\ priorValid \in BOOLEAN
    /\ provenValid \in BOOLEAN
    /\ exactValid \in BOOLEAN
    /\ transientA \in TransientControls
    /\ transientB \in TransientControls
    /\ sourceA = "none"
    /\ sourceB = "none"
    /\ targetA = "none"
    /\ targetB = "none"
    /\ handoff = "none"
    /\ planCount = 0
    /\ inputClosed = FALSE

EndInteraction ==
    /\ phase = "interacting"
    /\ sourceA' = SelectedSource(exactValid, retained,
                                  priorValid, provenValid)
    /\ sourceB' = SelectedSource(exactValid, retained,
                                  priorValid, provenValid)
    /\ targetA' = SelectedTarget(sourceA')
    /\ targetB' = SelectedTarget(sourceB')
    /\ handoff' = SelectedHandoff(sourceA', retained, transientA)
    /\ phase' = "restored"
    /\ UNCHANGED <<retained, priorValid, provenValid, exactValid,
                    transientA, transientB, planCount, inputClosed>>

AuthorizeSuccessor ==
    /\ phase = "restored"
    /\ handoff # "none"
    /\ planCount = 0
    /\ phase' = "planning"
    /\ planCount' = 1
    /\ UNCHANGED <<retained, priorValid, provenValid, exactValid,
                    transientA, transientB, sourceA, sourceB,
                    targetA, targetB, handoff, inputClosed>>

CompleteSuccessor ==
    /\ phase = "planning"
    /\ phase' = "terminal"
    /\ handoff' = "none"
    /\ UNCHANGED <<retained, priorValid, provenValid, exactValid,
                    transientA, transientB, sourceA, sourceB,
                    targetA, targetB, planCount, inputClosed>>

FinishWithoutSuccessor ==
    /\ phase = "restored"
    /\ handoff = "none"
    /\ phase' = "terminal"
    /\ UNCHANGED <<retained, priorValid, provenValid, exactValid,
                    transientA, transientB, sourceA, sourceB,
                    targetA, targetB, handoff, planCount,
                    inputClosed>>

CloseInput ==
    /\ ~inputClosed
    /\ inputClosed' = TRUE
    /\ UNCHANGED <<phase, retained, priorValid, provenValid, exactValid,
                    transientA, transientB, sourceA, sourceB,
                    targetA, targetB, handoff, planCount>>

Next ==
    \/ EndInteraction
    \/ AuthorizeSuccessor
    \/ CompleteSuccessor
    \/ FinishWithoutSuccessor
    \/ CloseInput

Spec ==
    /\ Init
    /\ [][Next]_vars
    /\ WF_vars(EndInteraction)
    /\ WF_vars(AuthorizeSuccessor)
    /\ WF_vars(CompleteSuccessor)
    /\ WF_vars(FinishWithoutSuccessor)

ScheduleIndependentTarget ==
    phase # "interacting" => sourceA = sourceB /\ targetA = targetB

CertificatePrecedence ==
    phase # "interacting" =>
        /\ (exactValid => sourceA = "exact")
        /\ (~exactValid /\ retained /\ priorValid => sourceA = "prior")
        /\ (~exactValid /\ retained /\ ~priorValid /\ provenValid =>
                sourceA = "proven")
        /\ (~exactValid /\ (~retained \/ (~priorValid /\ ~provenValid)) =>
                sourceA = "stable")

AtMostOneSuccessor == planCount <= 1

PriorPoseRequiresCurrentProof ==
    phase = "restored" /\ sourceA = "prior" => handoff = "presentation"

EventuallyTerminal == <> (phase = "terminal")

=============================================================================
