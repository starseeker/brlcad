--------------------------- MODULE obol_lod_control ---------------------------
\* A deliberately small TLC model of one CAD occurrence's control plane.
\*
\* This is not a renderer or mesh model.  It captures the contracts that are
\* easy to violate in the real asynchronous code: a source can first publish
\* a fallback part, later supply mesh data, or change the base part identity
\* without changing a numeric geometry revision.  A repair request is legal
\* only while a genuine LoD fallback remains presented.
\*
\* See libbobol_architecture.md for the source-level mapping and intended TLC
\* extensions (two occurrences and a bounded camera-epoch queue).
EXTENDS Naturals, TLC

Parts == {"fallback", "mesh", "authored"}

VARIABLES entryPart,
          activePart,
          registered,
	  registeredRole,
          payloadReady,
          requestOutstanding,
          queued,
          inFlight,
          resultReady,
          structuralBox,
          renderAcknowledged

vars == <<entryPart, activePart, registered, registeredRole, payloadReady,
          requestOutstanding, queued, inFlight, resultReady, structuralBox,
          renderAcknowledged>>

DesiredPart == IF payloadReady THEN "mesh" ELSE entryPart

TypeOK == /\ entryPart \in {"fallback", "authored"}
          /\ activePart \in Parts
          /\ registered \in [Parts -> BOOLEAN]
	  /\ registeredRole \in [Parts -> ({"none"} \union Parts)]
	  /\ \A p \in Parts : registered[p] = (registeredRole[p] # "none")
          /\ payloadReady \in BOOLEAN
          /\ requestOutstanding \in BOOLEAN
          /\ queued \in BOOLEAN
          /\ inFlight \in BOOLEAN
          /\ resultReady \in BOOLEAN
          /\ structuralBox \in BOOLEAN
          /\ renderAcknowledged \in BOOLEAN

Init == /\ entryPart = "fallback"
        /\ activePart = "fallback"
        /\ registered = [p \in Parts |-> p = "fallback"]
	/\ registeredRole = [p \in Parts |->
	       IF p = "fallback" THEN "fallback" ELSE "none"]
        /\ payloadReady = FALSE
        /\ requestOutstanding = TRUE
        /\ queued = TRUE
        /\ inFlight = FALSE
        /\ resultReady = FALSE
        /\ structuralBox = TRUE
        /\ renderAcknowledged = FALSE

StartLoad == /\ queued /\ ~inFlight
             /\ inFlight' = TRUE
             /\ queued' = FALSE
             /\ UNCHANGED <<entryPart, activePart, registered,
			    registeredRole, payloadReady,
                            requestOutstanding, resultReady, structuralBox,
                            renderAcknowledged>>

FinishLoad == /\ inFlight
              /\ inFlight' = FALSE
              /\ resultReady' = TRUE
              /\ registered' = [registered EXCEPT !["mesh"] = TRUE]
	      /\ registeredRole' =
		    [registeredRole EXCEPT !["mesh"] = "mesh"]
              /\ UNCHANGED <<entryPart, activePart, payloadReady,
                             requestOutstanding, queued, structuralBox,
                             renderAcknowledged>>

AcceptResult == /\ resultReady
                /\ payloadReady' = TRUE
                /\ resultReady' = FALSE
                /\ requestOutstanding' = FALSE
                /\ UNCHANGED <<entryPart, activePart, registered,
			       registeredRole, queued,
                               inFlight, structuralBox, renderAcknowledged>>

\* The key identity handoff: a source may replace a temporary fallback part
\* with an authored part while retaining the same numeric geometry revision.
\* The old bug skipped RegisterEntryPart on this transition.
SourceIdentityChange == /\ entryPart = "fallback"
                        /\ entryPart' = "authored"
                        /\ payloadReady' = FALSE
                        /\ requestOutstanding' = FALSE
                        /\ queued' = FALSE
                        /\ inFlight' = FALSE
                        /\ resultReady' = FALSE
                        \* The former box is no longer an LoD obligation as
                        \* soon as its source semantics change.  Publication
                        \* may still visually replace it on the next frame,
                        \* but repair admission must not see it as a reason
                        \* to submit zero useful work.
                        /\ structuralBox' = FALSE
			/\ UNCHANGED <<activePart, registered, registeredRole,
                                       renderAcknowledged>>

RegisterEntryPart == /\ ~registered[DesiredPart]
                     /\ registered' = [registered EXCEPT ![DesiredPart] = TRUE]
		     /\ registeredRole' =
			   [registeredRole EXCEPT ![DesiredPart] = DesiredPart]
                     /\ UNCHANGED <<entryPart, activePart, payloadReady,
                                    requestOutstanding, queued, inFlight,
                                    resultReady, structuralBox,
                                    renderAcknowledged>>

PublishDesiredPart == /\ registered[DesiredPart]
                       /\ activePart # DesiredPart
                       /\ activePart' = DesiredPart
                       /\ structuralBox' = (DesiredPart = "fallback")
                       /\ renderAcknowledged' = FALSE
		       /\ UNCHANGED <<entryPart, registered, registeredRole,
				      payloadReady,
                                      requestOutstanding, queued, inFlight,
                                      resultReady>>

RenderAcknowledge == /\ activePart = DesiredPart
                     /\ renderAcknowledged' = TRUE
		     /\ UNCHANGED <<entryPart, activePart, registered,
				    registeredRole,
                                    payloadReady, requestOutstanding, queued,
                                    inFlight, resultReady, structuralBox>>

\* Repair is intentionally unable to self-arm after a mesh or authored part
\* is presented.  It always starts actual loader work, so TLC can reject a
\* zero-work repair loop rather than treating it as harmless churn.
RequestRepair == /\ ~requestOutstanding
                 /\ ~queued
                 /\ ~inFlight
                 /\ ~resultReady
                 /\ structuralBox
                 /\ activePart = "fallback"
                 /\ entryPart = "fallback"
		 /\ ~payloadReady
                 /\ requestOutstanding' = TRUE
                 /\ queued' = TRUE
                 /\ renderAcknowledged' = FALSE
		 /\ UNCHANGED <<entryPart, activePart, registered,
				registeredRole,
                                payloadReady, inFlight, resultReady,
                                structuralBox>>

Next == \/ StartLoad
        \/ FinishLoad
        \/ AcceptResult
        \/ SourceIdentityChange
        \/ RegisterEntryPart
        \/ PublishDesiredPart
        \/ RenderAcknowledge
        \/ RequestRepair

Spec == Init /\ [][Next]_vars
        /\ WF_vars(StartLoad)
        /\ WF_vars(FinishLoad)
        /\ WF_vars(AcceptResult)
        /\ WF_vars(RegisterEntryPart)
        /\ WF_vars(PublishDesiredPart)
        /\ WF_vars(RenderAcknowledge)

\* Every live part is registered before it can be rendered.
PublishedPartRegistered == registered[activePart]

\* A stable PartId is not sufficient: the immutable geometry registered
\* behind it must still have the semantic role for which that ID was issued.
\* The production counterexample was an expired raw-address memo entry that
\* made an authored line inherit a fallback-box PartId after allocator reuse.
RegisteredPartSemantic ==
    \A p \in Parts : registered[p] => registeredRole[p] = p

\* The only legal visible box is a still-unresolved LoD fallback.
NoFalseStructuralObligation ==
    structuralBox => /\ activePart = "fallback"
                     /\ entryPart = "fallback"

\* Once a non-fallback target is available, fairness of registration and
\* publication guarantees a non-box acknowledged state.  This is the useful
\* liveness property for the real controller's idle contract.
EventuallySettled ==
    []((DesiredPart # "fallback") =>
       <>(activePart = DesiredPart /\ ~structuralBox /\ renderAcknowledged))

=============================================================================
