BRL-CAD Obol endpoint ownership and lifecycle
==============================================

Public boundary
---------------

A hosted GED view owns or borrows one `bobol_display_endpoint`.  Public code
attaches, replaces, queries, and detaches rendering state only with
`ged_view_context_display_endpoint_ensure`, `_get`, and `_set`.  Direct Obol
controller and scene-controller attachment adapters are private libged
implementation details.

Ownership diagram
-----------------

```
ged
 |- ged_drawable
 |   `- ged_obol_state                         (one shared-scene owner)
 |       `- shared BObolViewController         (owned until GED teardown)
 |           `- shared BObolSceneController    (all shared database sources)
 `- GED view host records[]                    (one record per hosted view)
     `- bobol_display_endpoint                 (owned or borrowed by record)
         `- BObolViewController                (owned according to endpoint flags)
             |- BObolViewAttachment            (also borrowed by host record)
             |- progressive providers          (controller-owned callback data)
             |- managed BObolLodService         (controller-owned, optional)
             `- render root
                 |- shared scene root           (borrowed when view is shared)
                 |- framebuffer interlay
                 `- view-local root
```

The per-GED shared controller is the sole owner of database scene state.  An
endpoint controller owns only its view-local stores, camera, framebuffer
layers, progressive provider, and composed render root.  Endpoint replacement
therefore never transfers or replays database sources: the new endpoint simply
borrows the unchanged shared scene root.  The GED shared-scene state does not
mirror endpoint controllers; the live display endpoint is their sole registry.

Lifecycle
---------

1. `ged_view_context_host_attach` associates the view with one GED.
2. Endpoint ensure/set creates the shared scene on demand, binds the endpoint
   controller attachment to the view, registers controller-owned progressive
   state, and composes its render root.
3. Setting a different endpoint unregisters the old per-view services, detaches
   framebuffer providers, destroys the old endpoint only when the view owned
   it, and binds the replacement to the same shared scene.
4. Setting `NULL` performs the same per-view detach without destroying a
   borrowed endpoint.  Other views and the shared scene remain live.
5. View destruction detaches its endpoint before freeing the host record.  GED
   destruction first detaches all view endpoints, then destroys the shared
   scene owner and transaction observer.

Visibility changes are endpoint-host properties and do not affect ownership.
A hidden endpoint stays attached and receives scene transactions; showing it
again presents the retained shared scene.

Internal binding invariants
---------------------------

`ged_obol_state` contains only the shared controller, its transaction-observer
token, and the shared scene's one-time bootstrap state.  Per-view
synchronization visits live GED host records and derives each controller from
its display endpoint.  Independent/shared composition is determined from
current view policy and the render graph, not a cached bridge flag.

Progressive callback data and managed LoD services are controller-owned and are
released by provider unregistration or controller destruction.  A host record
borrows the controller attachment while the endpoint is installed and restores
a value-only attachment when it is removed.  There are no parallel
primary-controller pointers, per-view attachment vectors, ownership flags, or
source preservation/replay side channels on `ged_drawable`.
