-------------------- MODULE ObolSparsePlanIdentity --------------------
\* Stable identity ownership across sparse retained-plan replacement.
\*
\* A streamed replacement hides its compiled slot and appends a new active
\* slot with the same InstanceId.  Later cut-bin maintenance may reorder the
\* hidden slot inside its old group.  Such a swap may move an ID-to-index map
\* entry only when that entry still owns the swapped source slot; otherwise a
\* tombstone can steal the active replacement's index.

EXTENDS FiniteSets, Naturals, TLC

CONSTANT MaxSwaps

Identities == {"target", "peer"}
Slots == 0..2
SlotIdentities == Identities \union {"none"}
Phases == {"compiled", "reordering", "patched", "done"}

VARIABLES phase,
          slotIdentity,
          slotActive,
          indexByIdentity,
          swaps

vars == <<phase, slotIdentity, slotActive, indexByIdentity, swaps>>

TypeOK ==
    /\ phase \in Phases
    /\ slotIdentity \in [Slots -> SlotIdentities]
    /\ slotActive \in [Slots -> BOOLEAN]
    /\ indexByIdentity \in [Identities -> Slots]
    /\ swaps \in 0..MaxSwaps

Init ==
    /\ phase = "compiled"
    /\ slotIdentity = [slot \in Slots |->
        IF slot = 0 THEN "peer"
        ELSE IF slot = 1 THEN "target"
        ELSE "none"]
    /\ slotActive = [slot \in Slots |-> slot # 2]
    /\ indexByIdentity = [identity \in Identities |->
        IF identity = "peer" THEN 0 ELSE 1]
    /\ swaps = 0

\* Slot 1 remains as the hidden compiled target record.  Slot 2 is the
\* authoritative active replacement and therefore acquires the mapping.
RebindTarget ==
    /\ phase = "compiled"
    /\ phase' = "reordering"
    /\ slotIdentity' = [slotIdentity EXCEPT ![2] = "target"]
    /\ slotActive' = [slotActive EXCEPT ![1] = FALSE, ![2] = TRUE]
    /\ indexByIdentity' =
        [indexByIdentity EXCEPT !["target"] = 2]
    /\ swaps' = 0

\* The old cut group owns slots 0 and 1.  Update only mappings which still
\* point at one of those source slots.  In particular target->2 is unchanged.
SwapOldGroup ==
    /\ phase = "reordering"
    /\ swaps < MaxSwaps
    /\ slotIdentity' = [slot \in Slots |->
        IF slot = 0 THEN slotIdentity[1]
        ELSE IF slot = 1 THEN slotIdentity[0]
        ELSE slotIdentity[slot]]
    /\ slotActive' = [slot \in Slots |->
        IF slot = 0 THEN slotActive[1]
        ELSE IF slot = 1 THEN slotActive[0]
        ELSE slotActive[slot]]
    /\ indexByIdentity' = [identity \in Identities |->
        LET source == indexByIdentity[identity]
        IN IF source = 0 THEN 1
           ELSE IF source = 1 THEN 0
           ELSE source]
    /\ swaps' = swaps + 1
    /\ UNCHANGED phase

\* A later PoP cut update must resolve target to its active replacement.
PatchTarget ==
    /\ phase = "reordering"
    /\ slotActive[indexByIdentity["target"]]
    /\ slotIdentity[indexByIdentity["target"]] = "target"
    /\ phase' = "patched"
    /\ UNCHANGED <<slotIdentity, slotActive, indexByIdentity, swaps>>

Finish ==
    /\ phase = "patched"
    /\ phase' = "done"
    /\ UNCHANGED <<slotIdentity, slotActive, indexByIdentity, swaps>>

Done ==
    /\ phase = "done"
    /\ UNCHANGED vars

Next ==
    \/ RebindTarget
    \/ SwapOldGroup
    \/ PatchTarget
    \/ Finish
    \/ Done

Spec == Init /\ [][Next]_vars
        /\ WF_vars(RebindTarget)
        /\ WF_vars(PatchTarget)
        /\ WF_vars(Finish)

ActiveIndexOwned ==
    \A identity \in Identities:
        /\ slotActive[indexByIdentity[identity]]
        /\ slotIdentity[indexByIdentity[identity]] = identity

OneActiveSlotPerIdentity ==
    \A identity \in Identities:
        Cardinality({slot \in Slots:
            slotActive[slot] /\ slotIdentity[slot] = identity}) = 1

EventuallyPatched == <>[](phase \in {"patched", "done"})

=======================================================================
