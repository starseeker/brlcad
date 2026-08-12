/*                         E V E N T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file ged/event.h
 *
 * Renderer-neutral database mutation and reconciliation events.
 *
 * Commands publish post-mutation semantic events instead of repairing draw,
 * selection, edit, or application state independently.  Events are copied by
 * the owning GED context, coalesced at the outermost batch boundary, and then
 * reconciled before post-reconcile observers run.
 */

#ifndef GED_EVENT_H
#define GED_EVENT_H

#include <stddef.h>
#include <stdint.h>

#include "ged/defines.h"

struct ged_event_result;

/** Kind of semantic database mutation represented by an event. */
enum ged_event_kind {
    GED_EVENT_NONE = 0,
    GED_EVENT_OBJECT_ADDED,
    GED_EVENT_OBJECT_REMOVED,
    GED_EVENT_OBJECT_RENAMED,
    GED_EVENT_OBJECT_MODIFIED,
    GED_EVENT_OBJECT_VISIBILITY_CHANGED,
    GED_EVENT_COMB_TREE_CHANGED,
    GED_EVENT_COMB_INSTANCE_REMOVED,
    GED_EVENT_OBJECT_REFERENCES_REMOVED,
    GED_EVENT_ATTRIBUTE_CHANGED,
    GED_EVENT_MATERIAL_CHANGED,
    GED_EVENT_BATCH_REBUILD,
    GED_EVENT_DATABASE_METADATA_CHANGED
};

/** Observer point relative to core GED reconciliation. */
enum ged_event_observer_phase {
    GED_EVENT_OBSERVER_INTERNAL = 0,
    GED_EVENT_OBSERVER_POST_RECONCILE
};

/** Status returned by event mutation operations. */
enum ged_event_status {
    GED_EVENT_OK = 0,
    GED_EVENT_INVALID,
    GED_EVENT_UNAVAILABLE,
    GED_EVENT_ERROR
};

/** Opaque token identifying one observer registration. */
typedef uintptr_t ged_event_observer_token;

/**
 * Borrowed input describing one semantic database mutation.
 *
 * Every string is borrowed only for the duration of the publishing call.
 * NULL strings mean that the corresponding identity is not known or does not
 * apply to this event kind.  Use ged_event_init() before assigning fields.
 */
struct ged_event {
    enum ged_event_kind kind;
    const char *name;
    const char *new_name;
    const char *parent_name;
    const char *child_name;
    const char *path;
    int redraw;
};

/**
 * Callback invoked synchronously for a reconciled event batch.
 *
 * @p events, @p result, and all borrowed values reachable from them remain
 * valid only for the callback.  An observer may remove itself while running;
 * events published reentrantly are deferred until the current dispatch ends.
 */
typedef void (*ged_event_observer_func_t)(
    struct ged *gedp,
    const struct ged_event *events,
    size_t event_count,
    const struct ged_event_result *result,
    void *client_data);

__BEGIN_DECLS

/** Initialize a caller-owned event descriptor to GED_EVENT_NONE. */
GED_EXPORT extern void
ged_event_init(struct ged_event *event);

/** Allocate an empty opaque event result. */
GED_EXPORT extern struct ged_event_result *
ged_event_result_create(void);

/** Clear a reusable event result while retaining its allocation. */
GED_EXPORT extern void
ged_event_result_clear(struct ged_event_result *result);

/** Destroy an event result allocated by ged_event_result_create(). */
GED_EXPORT extern void
ged_event_result_destroy(struct ged_event_result *result);

/** Return the operation status stored in an event result. */
GED_EXPORT extern enum ged_event_status
ged_event_result_status(const struct ged_event_result *result);

/** Return the number of input events represented by a result. */
GED_EXPORT extern size_t
ged_event_result_event_count(const struct ged_event_result *result);

/** Return the number of events remaining after semantic coalescing. */
GED_EXPORT extern size_t
ged_event_result_coalesced_count(const struct ged_event_result *result);

/** Return the number of semantic draw changes made during reconciliation. */
GED_EXPORT extern size_t
ged_event_result_draw_change_count(const struct ged_event_result *result);

/** Return the number of database-index changes made during reconciliation. */
GED_EXPORT extern size_t
ged_event_result_db_index_change_count(
    const struct ged_event_result *result);

/** Return non-zero when selection reconciliation changed semantic state. */
GED_EXPORT extern int
ged_event_result_selection_changed(const struct ged_event_result *result);

/** Return the semantic scene revision before reconciliation. */
GED_EXPORT extern uint64_t
ged_event_result_scene_revision_before(
    const struct ged_event_result *result);

/** Return the semantic scene revision after reconciliation. */
GED_EXPORT extern uint64_t
ged_event_result_scene_revision_after(
    const struct ged_event_result *result);

/** Return the number of structured affected paths in a result. */
GED_EXPORT extern size_t
ged_event_result_path_count(const struct ged_event_result *result);

/** Return a result-owned affected path, or NULL when @p index is invalid. */
GED_EXPORT extern const char *
ged_event_result_path_at(const struct ged_event_result *result, size_t index);

/** Return a borrowed diagnostic string, or an empty string when none exists. */
GED_EXPORT extern const char *
ged_event_result_diagnostic(const struct ged_event_result *result);

/** Return non-zero when semantic event reconciliation is enabled. */
GED_EXPORT extern int
ged_event_available(struct ged *gedp);

/** Suspend semantic event reconciliation and discard queued work. */
GED_EXPORT extern enum ged_event_status
ged_event_disable(struct ged *gedp);

/** Release one event-reconciliation suspension. */
GED_EXPORT extern enum ged_event_status
ged_event_enable(struct ged *gedp);

/** Return non-zero when the event service has a live reconciliation consumer. */
GED_EXPORT extern int
ged_event_has_live_consumers(struct ged *gedp);

/** Begin a bulk database mutation scope. */
GED_EXPORT extern enum ged_event_status
ged_event_bulk_begin(struct ged *gedp);

/** End a bulk mutation scope and publish one rebuild event when warranted. */
GED_EXPORT extern enum ged_event_status
ged_event_bulk_end(struct ged *gedp, struct ged_event_result *result);

/** Return non-zero when a bulk database mutation scope is active. */
GED_EXPORT extern int
ged_event_bulk_active(struct ged *gedp);

/** Enable the fallback librt database-change callback bridge. */
GED_EXPORT extern enum ged_event_status
ged_event_librt_callbacks_enable(struct ged *gedp);

/** Disable the fallback librt database-change callback bridge. */
GED_EXPORT extern enum ged_event_status
ged_event_librt_callbacks_disable(struct ged *gedp);

/** Begin a nestable semantic event batch. */
GED_EXPORT extern enum ged_event_status
ged_event_batch_begin(struct ged *gedp);

/** End an event batch and reconcile at its outermost boundary. */
GED_EXPORT extern enum ged_event_status
ged_event_batch_end(struct ged *gedp, struct ged_event_result *result);

/** Publish and copy a sequence of typed semantic database events. */
GED_EXPORT extern enum ged_event_status
ged_event_publish(struct ged *gedp,
		  const struct ged_event *events,
		  size_t event_count,
		  struct ged_event_result *result);

/** Publish an object-added event. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_object_added(struct ged *gedp, const char *name,
			      struct ged_event_result *result);

/** Publish an object-removed event. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_object_removed(struct ged *gedp, const char *name,
				struct ged_event_result *result);

/** Publish an object-renamed event. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_object_renamed(struct ged *gedp, const char *old_name,
				const char *new_name,
				struct ged_event_result *result);

/** Publish an object-modified event with optional immediate redraw. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_object_modified(struct ged *gedp, const char *name,
				 int redraw,
				 struct ged_event_result *result);

/** Publish an object-visibility event. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_object_visibility_changed(
    struct ged *gedp, const char *name, struct ged_event_result *result);

/** Publish a combination-tree event with optional immediate redraw. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_comb_tree_changed(struct ged *gedp, const char *name,
				   int redraw,
				   struct ged_event_result *result);

/** Publish removal of one combination child instance. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_comb_instance_removed(
    struct ged *gedp, const char *parent_name, const char *child_name,
    const char *path, struct ged_event_result *result);

/** Publish removal of all references to an object. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_object_references_removed(
    struct ged *gedp, const char *name, struct ged_event_result *result);

/** Publish removal of an object reference from one parent/path. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_object_reference_removed_from_parent(
    struct ged *gedp, const char *name, const char *parent_name,
    const char *path, struct ged_event_result *result);

/** Publish an attribute-change event with optional immediate redraw. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_attribute_changed(struct ged *gedp, const char *name,
				   int redraw,
				   struct ged_event_result *result);

/** Publish a database-wide material-change event. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_material_changed(struct ged *gedp,
				  struct ged_event_result *result);

/** Publish a material-change event for one named database object. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_object_material_changed(
    struct ged *gedp, const char *name, struct ged_event_result *result);

/** Publish a complete semantic index/scene rebuild event. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_batch_rebuild(struct ged *gedp,
			       struct ged_event_result *result);

/** Publish a database-metadata-only change event. */
GED_EXPORT extern enum ged_event_status
ged_event_notify_database_metadata_changed(
    struct ged *gedp, struct ged_event_result *result);

/** Register an owner-context event observer. */
GED_EXPORT extern ged_event_observer_token
ged_event_observer_add(struct ged *gedp,
		       enum ged_event_observer_phase phase,
		       ged_event_observer_func_t func,
		       void *client_data);

/** Remove an event observer; removal during dispatch is supported. */
GED_EXPORT extern int
ged_event_observer_remove(struct ged *gedp, ged_event_observer_token token);

__END_DECLS

#endif /* GED_EVENT_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
