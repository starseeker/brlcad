------------------------------ MODULE LodCommon ------------------------------
EXTENDS Integers, FiniteSets, Sequences

\* Shared finite-domain vocabulary for the BObol LoD specifications.  The
\* strings deliberately match the C++ enum names closely enough that a TLC
\* trace can be translated back into a coordinator event sequence by hand.

Fallback   == "FALLBACK"
Coverage   == "COVERAGE"
Interactive == "INTERACTIVE"
Settling   == "SETTLING"
Stable     == "STABLE"
Compacting == "COMPACTING"

Phases == {Fallback, Coverage, Interactive, Settling, Stable, Compacting}

Initialize          == "INITIALIZE"
FrameCompleted      == "FRAME_COMPLETED"
WorkScheduled       == "WORK_SCHEDULED"
WorkPumped          == "WORK_PUMPED"
ResultPublished     == "RESULT_PUBLISHED"
ServiceChanged      == "SERVICE_CHANGED"
GenerationCancelled == "GENERATION_CANCELLED"
AutoSubmitChanged   == "AUTO_SUBMIT_CHANGED"
ViewInvalidated     == "VIEW_INVALIDATED"
PolicyChanged       == "POLICY_CHANGED"
InteractionStarted  == "INTERACTION_STARTED"
InteractionEnded    == "INTERACTION_ENDED"
ViewObserved        == "VIEW_OBSERVED"

Events == {
    Initialize, FrameCompleted, WorkScheduled, WorkPumped, ResultPublished,
    ServiceChanged, GenerationCancelled, AutoSubmitChanged, ViewInvalidated,
    PolicyChanged, InteractionStarted, InteractionEnded, ViewObserved
}

MinNat(a, b) == IF a <= b THEN a ELSE b
MaxNat(a, b) == IF a >= b THEN a ELSE b

=============================================================================
