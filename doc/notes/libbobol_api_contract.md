# libBObol API, ownership, and ABI contract

This document is the reviewed contract for libBObol 1.x.  The header lists in
`include/BObol/api_tiers.cmake` are the machine-readable source of truth used
by installation and header-compilation tests.

## API tiers

The **stable source tier** is the host-facing surface included by `BObol.h`:

`BDatabaseSource.h`, `BDefines.h`, `BDisplayEndpoint.h`,
`BDisplaySession.h`, `BExportAction.h`, `BFramebuffer.h`,
`BHeadlessWindowHost.h`, `BHostFactory.h`, `BInit.h`, `BMeasureAction.h`,
`BInput.h`, `BPickDetail.h`, `BRtRender.h`, `BSceneController.h`, `BSceneGroup.h`,
`BSnapAction.h`, `BSourceMeshRequest.h`, `BViewController.h`,
`BViewQuery.h`, `BViewStore.h`, and `BWindowHost.h`.

The **advanced source tier** is installed for BRL-CAD drawing integrations but
is not pulled in by the umbrella header:

`BADC.h`, `BAxes.h`, `BDrawCache.h`, `BEditPreview.h`,
`BEvaluatedPoints.h`, `BGrid.h`, `BHUDLabelOverlay.h`, `BImagePlane.h`,
`BImageSource.h`, `BLineLayerOverlay.h`, `BLodMeshShape.h`,
`BLodRealization.h`, `BLodService.h`, `BLodUpdateAction.h`,
`BMaterialObject.h`, `BMeshLodCache.h`, `BMeshLodSubmitAction.h`,
`BMeshResidencyAction.h`, `BMeshShape.h`, `BPerformance.h`,
`BRealizeAction.h`, `BViewAttachment.h`, `BViewLod.h`,
`BViewportImage.h`, and `BVListShape.h`.

Headers under `src/libBObol`, including realization repositories, performance
implementation helpers, and CAD-assembly build state, are private.  They are
not installed and their symbols must remain hidden.

## ABI policy

The stable C entry points and opaque C handles are the binary compatibility
boundary for libBObol 1.x.  Exported host-facing C++ facades use pImpl where
their Coin runtime contract does not require public fields.  Coin node/action
classes and all advanced C++ APIs are source-supported but their C++ ABI is
experimental: consumers must rebuild them with the matching BRL-CAD release.
An incompatible change to the stable C ABI requires a libBObol SOVERSION bump.
The current `SOVERSION 1` therefore promises the stable C ABI, not a frozen ABI
for every installed C++ class.

The machine-readable `LIBBOBOL_STABLE_C_ABI_HEADERS` list identifies the
headers parsed by the zero-warning documentation gate.  Full reference output
also includes the stable-source and advanced C++ tiers, but undocumented
experimental C++ declarations do not silently expand the SOVERSION promise.

## Owners and lifetime

The display endpoint is the sole attachment boundary.  A host owns or borrows
an endpoint according to the endpoint creation flags; the endpoint owns its
per-view `BObolViewController`.  A GED owns one shared scene controller and
may attach that scene to several endpoints.  Replacing or destroying an
endpoint never transfers scene ownership.

Feature, polygon, compact-instance, and GED scene references are values.  They
contain an owner identity, object identity, and generation or revision.  A
value does not retain its owner.  Every operation resolves all identity fields
before using backing storage, and a removed object, rebuilt store, replaced
controller, or changed scene rejects the stale value.

Pointers returned by `get*`, `find*`, or controller/store accessors are
borrowed.  Unless a declaration explicitly says otherwise, they remain valid
only until the next mutation of that owner and must not be deleted or retained
past owner teardown.

## Thread and callback rules

Scene nodes, controllers, stores, endpoints, and rendering providers are
owner-thread objects.  Their public mutating and query methods run on the host
owner thread.  LoD workers receive database snapshots or plain request values
and return plain result data; they do not mutate Coin nodes, global type state,
or an OpenGL context.  The owner thread drains and publishes worker results.

Frame-request and presentation callbacks may be invoked while the controller
is active, but never after callback removal returns.  A callback must not
destroy the controller from within itself.  It may queue work or request a
subsequent frame.  Publication callbacks receive stack-owned context and must
not retain it.

## Naming and result glossary

`BObol*` names are services and value types.  `SoBRL*` names are Inventor
nodes, actions, details, or elements.  `bobol_*` names are C ABI functions and
opaque C values.  Variables use their role—`endpoint`, `scene`, `source`,
`feature`, and `handle`—rather than a historical renderer layer.

Return values follow these rules:

| Form | Meaning |
|---|---|
| `SbBool` / `bool` | success or a documented predicate; no count is encoded |
| pointer | borrowed/created object as documented; null means unavailable or failure |
| `size_t` getter | a count; it never encodes an error |
| `int` mutator returning `-1/0/1` | error / no change / changed |
| `int` collection operation | negative error, otherwise an explicit count |
| `BRLCAD_OK` / `BRLCAD_ERROR` C adapter | command-style success / failure |

Declarations whose historical integer convention differs must state it at the
declaration.  New stable APIs use a typed status when more than the states in
the table are needed; they do not overload one integer as an undocumented
boolean, count, and error.

## Coordinates and units

Database-source geometry is source-local.  `drawMatrix` and compact occurrence
transforms carry placement separately.  View-controller camera and clipping
methods use BRL-CAD model units; screen coordinates are pixels in the bound
viewport; colors use normalized Inventor channels unless a C declaration says
RGB bytes.  LoD time budgets are microseconds and render timing values are
nanoseconds, as named in their declarations.
