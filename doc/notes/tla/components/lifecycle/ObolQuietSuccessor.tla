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
Targets == {"stableQuality", "priorQuality", "provenQuality",
            "provenPresentationStableDemand", "exactQuality"}
Handoffs == {"none", "allocation", "presentation"}
ChangeKinds == {"pose", "scale"}

SelectedSource(exact, retained, prior, proven, changeKind) ==
    IF exact THEN "exact"
    ELSE IF ~retained THEN "stable"
    \* A completed frame at the new scale is stronger than controls captured
    \* before the zoom.  Without that proof, the prior stable control vector
    \* is still the deterministic starting point for orthographic demand.
    ELSE IF proven /\ changeKind = "scale" THEN "proven"
    ELSE IF prior THEN "prior"
    ELSE "stable"

SelectedTarget(source, prior) ==
    IF source = "exact" THEN "exactQuality"
    ELSE IF source = "prior" THEN "priorQuality"
    ELSE IF source = "proven" /\ prior
         THEN "provenPresentationStableDemand"
    ELSE IF source = "proven" THEN "provenQuality"
    ELSE "stableQuality"

SelectedHandoff(source, retained, transient, changeKind) ==
    IF source = "prior" /\ changeKind = "pose" THEN "presentation"
    ELSE IF source = "prior" THEN "allocation"
    ELSE IF retained /\ source \in {"stable", "proven"} /\
            transient # "rich" THEN "allocation"
    ELSE "none"

VARIABLES phase,
          retained,
          changeKind,
          populationSame,
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

vars == <<phase, retained, changeKind, populationSame, priorValid,
          provenValid, exactValid,
          transientA, transientB, sourceA, sourceB, targetA, targetB,
          handoff, planCount, inputClosed>>

TypeOK ==
    /\ phase \in Phases
    /\ retained \in BOOLEAN
    /\ changeKind \in ChangeKinds
    /\ populationSame \in BOOLEAN
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
    /\ changeKind \in ChangeKinds
    /\ populationSame \in BOOLEAN
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
                                  priorValid, provenValid, changeKind)
    /\ sourceB' = SelectedSource(exactValid, retained,
                                  priorValid, provenValid, changeKind)
    /\ targetA' = SelectedTarget(sourceA', priorValid)
    /\ targetB' = SelectedTarget(sourceB', priorValid)
    /\ handoff' = IF sourceA' = "prior" /\
                       (changeKind = "scale" \/ ~populationSame)
                   THEN "allocation"
                   ELSE SelectedHandoff(sourceA', retained, transientA,
                                        changeKind)
    /\ phase' = "restored"
    /\ UNCHANGED <<retained, changeKind, populationSame, priorValid,
                    provenValid, exactValid,
                    transientA, transientB, planCount, inputClosed>>

AuthorizeSuccessor ==
    /\ phase = "restored"
    /\ handoff # "none"
    /\ planCount = 0
    /\ phase' = "planning"
    /\ planCount' = 1
    /\ UNCHANGED <<retained, changeKind, populationSame, priorValid,
                    provenValid, exactValid,
                    transientA, transientB, sourceA, sourceB,
                    targetA, targetB, handoff, inputClosed>>

CompleteSuccessor ==
    /\ phase = "planning"
    /\ phase' = "terminal"
    /\ handoff' = "none"
    /\ UNCHANGED <<retained, changeKind, populationSame, priorValid,
                    provenValid, exactValid,
                    transientA, transientB, sourceA, sourceB,
                    targetA, targetB, planCount, inputClosed>>

FinishWithoutSuccessor ==
    /\ phase = "restored"
    /\ handoff = "none"
    /\ phase' = "terminal"
    /\ UNCHANGED <<retained, changeKind, populationSame, priorValid,
                    provenValid, exactValid,
                    transientA, transientB, sourceA, sourceB,
                    targetA, targetB, handoff, planCount,
                    inputClosed>>

CloseInput ==
    /\ ~inputClosed
    /\ inputClosed' = TRUE
    /\ UNCHANGED <<phase, retained, changeKind, populationSame, priorValid,
                    provenValid, exactValid,
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
    /\ WF_vars(AuthorizeSuccessor)
    /\ WF_vars(CompleteSuccessor)
    /\ WF_vars(FinishWithoutSuccessor)

ScheduleIndependentTarget ==
    phase # "interacting" => sourceA = sourceB /\ targetA = targetB

CertificatePrecedence ==
    phase # "interacting" =>
        /\ (exactValid => sourceA = "exact")
        /\ (~exactValid /\ retained /\ changeKind = "scale" /\
             provenValid =>
                sourceA = "proven")
        /\ (~exactValid /\ retained /\
             ~(changeKind = "scale" /\ provenValid) /\ priorValid =>
                sourceA = "prior")
        /\ (~exactValid /\ (~retained \/
             (~priorValid /\ ~(changeKind = "scale" /\ provenValid))) =>
                sourceA = "stable")

AtMostOneSuccessor == planCount <= 1

PriorPoseRequiresCurrentProof ==
    phase = "restored" /\ sourceA = "prior" /\ changeKind = "pose" /\
        populationSame =>
        handoff = "presentation"

PriorScaleRequiresAllocation ==
    phase = "restored" /\ sourceA = "prior" /\
        (changeKind = "scale" \/ ~populationSame) =>
        handoff = "allocation"

ProvenScaleKeepsPriorDemand ==
    phase # "interacting" /\ sourceA = "proven" /\ priorValid =>
        targetA = "provenPresentationStableDemand"

Ready == phase = "terminal"
Constrained == FALSE
Failed == FALSE
Terminal == Ready \/ Constrained \/ Failed

DeadlockOnlyAtTerminal == ~ENABLED <<Next>>_vars => Terminal

EventuallyTerminalAfterInteraction ==
    [](phase # "interacting" => <> (phase = "terminal"))

=============================================================================
