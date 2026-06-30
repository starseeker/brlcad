/*                    D R A W _ O B O L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file ged/draw_obol.h
 *
 * C++ bridge between libged's neutral draw transaction API and an
 * Obol/libbrlobol view controller.
 */

#ifndef GED_DRAW_OBOL_H
#define GED_DRAW_OBOL_H

#include "common.h"

#include <stdint.h>

#include "ged/defines.h"

#ifdef __cplusplus

class BRLObolViewController;
class SoBRLSceneController;
struct ged_draw_transaction;
struct ged_draw_transaction_result;

/**
 * Attach a borrowed Obol scene controller to @p gedp's draw state.
 *
 * The controller remains owned by the caller.  While attached, successful GED
 * draw transactions are mirrored into the controller's database-source scene.
 * When @p sync_current_scene is non-zero, the current retained GED draw set is
 * also synchronized immediately.
 *
 * This is the preferred Capability 1 boundary: libged mutates the
 * libbrlobol/Obol scene model without depending on GUI view policy.
 *
 * Returns 1 on successful attachment, 0 on invalid input or observer
 * registration failure.
 */
GED_EXPORT int
ged_draw_obol_scene_controller_attach(struct ged *gedp,
				      SoBRLSceneController *controller,
				      int sync_current_scene = 1);

/**
 * Stop mirroring GED draw transactions to the currently attached Obol scene
 * controller.
 */
GED_EXPORT void
ged_draw_obol_scene_controller_detach(struct ged *gedp);

/**
 * Return the borrowed scene controller currently attached to @p gedp, or NULL.
 */
GED_EXPORT SoBRLSceneController *
ged_draw_obol_scene_controller(struct ged *gedp);

/**
 * Return the active Obol scene controller, creating a libged-owned controller
 * when none is attached.
 *
 * The owned controller is synchronized from the current GED draw set when
 * @p sync_current_scene is non-zero, receives subsequent draw transactions
 * through the same observer path as borrowed controllers, and is destroyed
 * automatically by ged_draw_obol_scene_controller_detach or GED teardown.
 */
GED_EXPORT SoBRLSceneController *
ged_draw_obol_scene_controller_ensure(struct ged *gedp,
				      int sync_current_scene = 1);

/**
 * Return non-zero when @p gedp's active Obol scene controller is owned by
 * libged rather than borrowed from qtcad/qged or another caller.
 */
GED_EXPORT int
ged_draw_obol_scene_controller_owned(struct ged *gedp);

/**
 * Mirror one completed GED draw transaction into @p controller.
 *
 * If @p controller is NULL, the currently attached scene controller is used.
 * Returns 1 if the Obol scene changed, 0 otherwise.
 */
GED_EXPORT int
ged_draw_obol_scene_sync_transaction(
	struct ged *gedp,
	const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result,
	SoBRLSceneController *controller = NULL);

/**
 * Rebuild @p controller's database-source scene from the current GED draw set.
 *
 * If @p controller is NULL, the currently attached scene controller is used.
 * A zero @p source_revision requests a folded value derived from GED's current
 * draw scene revision.  Returns 1 if the Obol scene changed, 0 otherwise.
 */
GED_EXPORT int
ged_draw_obol_scene_sync_full_scene(struct ged *gedp,
				    void *view_ctx,
				    uint32_t source_revision = 0,
				    SoBRLSceneController *controller = NULL);

/**
 * Attach a borrowed Obol view controller to @p gedp's draw state.
 *
 * This compatibility wrapper attaches the view controller's scene controller.
 * New libged integrations should prefer ged_draw_obol_scene_controller_attach.
 */
GED_EXPORT int
ged_draw_obol_controller_attach(struct ged *gedp,
				BRLObolViewController *controller,
				int sync_current_scene = 1);

/**
 * Stop mirroring GED draw transactions to the currently attached controller.
 */
GED_EXPORT void
ged_draw_obol_controller_detach(struct ged *gedp);

/**
 * Return the borrowed controller currently attached to @p gedp, or NULL.
 */
GED_EXPORT BRLObolViewController *
ged_draw_obol_controller(struct ged *gedp);

/**
 * Mirror one completed GED draw transaction into @p controller.
 *
 * This compatibility wrapper delegates to ged_draw_obol_scene_sync_transaction.
 * If @p controller is NULL, the currently attached scene controller is used.
 * Returns 1 if the Obol scene changed, 0 otherwise.
 */
GED_EXPORT int
ged_draw_obol_sync_transaction(
	struct ged *gedp,
	const struct ged_draw_transaction *txn,
	const struct ged_draw_transaction_result *result,
	BRLObolViewController *controller = NULL);

/**
 * Rebuild @p controller's database-source scene from the current GED draw set.
 *
 * This compatibility wrapper delegates to ged_draw_obol_scene_sync_full_scene.
 * If @p controller is NULL, the currently attached scene controller is used.  A
 * zero @p source_revision requests a folded value derived from GED's current
 * draw scene revision.  Returns 1 if the Obol scene changed, 0 otherwise.
 */
GED_EXPORT int
ged_draw_obol_sync_full_scene(struct ged *gedp,
			      void *view_ctx,
			      uint32_t source_revision = 0,
			      BRLObolViewController *controller = NULL);

#endif /* __cplusplus */

#endif /* GED_DRAW_OBOL_H */
