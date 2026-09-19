/*                        V I E W . H
 * BRL-CAD
 *
 * Copyright (c) 2008-2026 United States Government as represented by
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
/** @addtogroup ged_view
 *
 * Geometry EDiting Library Database View Related Functions.
 *
 */
/** @{ */
/** @file ged/view.h */

#ifndef GED_VIEW_H
#define GED_VIEW_H

#include "bn/tol.h"
#include "bg/polygon.h"
#include "rt/db_fullpath.h"
#include "rt/db_instance.h"
#include "rt/view.h"
#include "ged/defines.h"
#include "ged/display.h"
#include "ged/scene.h"
#include "ged/view_export.h"

__BEGIN_DECLS

struct bu_ptbl;

typedef void (*ged_view_context_update_callback_t)(struct ged_view_context *view, void *data);

typedef enum ged_view_clear_flags {
    GED_VIEW_CLEAR_DB = 0x01,
    GED_VIEW_CLEAR_VIEW = 0x02,
    GED_VIEW_CLEAR_LOCAL = 0x04
} ged_view_clear_flags;

#define GED_VIEW_REFRESH_VIEW        0x00000001u
#define GED_VIEW_REFRESH_DRAW        0x00000002u
#define GED_VIEW_REFRESH_EDIT        0x00000004u
#define GED_VIEW_REFRESH_FRAMEBUFFER 0x00000008u
#define GED_VIEW_REFRESH_OVERLAY     0x00000010u
#define GED_VIEW_REFRESH_FORCE       0x80000000u
#define GED_VIEW_REFRESH_ALL         0xffffffffu

/** Check if a drawable exists */
#define GED_CHECK_DRAWABLE(_gedp, _flags) \
    if (!ged_scene_available(_gedp)) { \
	int ged_check_drawable_quiet = (_flags) & GED_QUIET; \
	if (!ged_check_drawable_quiet) { \
	    bu_vls_trunc((_gedp)->ged_result_str, 0); \
	    bu_vls_printf((_gedp)->ged_result_str, "A drawable does not exist."); \
	} \
	return (_flags); \
    }

/** Check if a view exists */
#define GED_CHECK_VIEW(_gedp, _flags) \
    if (ged_view_active_ctx(_gedp) == GED_VIEW_NULL) { \
	int ged_check_view_quiet = (_flags) & GED_QUIET; \
	if (!ged_check_view_quiet) { \
	    bu_vls_trunc((_gedp)->ged_result_str, 0); \
	    bu_vls_printf((_gedp)->ged_result_str, "A view does not exist."); \
	} \
	return (_flags); \
    }

/**
 * Rotate angle degrees about the specified axis
 */
GED_EXPORT extern int ged_arot_args(struct ged *gedp, int argc, const char *argv[], mat_t rmat);

/**
 * Rotate the view.
 */
GED_EXPORT extern int ged_rot_args(struct ged *gedp, int argc, const char *argv[], char *coord, mat_t rmat);

/**
 * Scale the view.
 */
GED_EXPORT extern int ged_scale_args(struct ged *gedp, int argc, const char *argv[], fastf_t *sf1, fastf_t *sf2, fastf_t *sf3);

/** Return non-zero if @p view owns an independent semantic scene scope. */
GED_EXPORT extern int ged_view_context_is_independent(const struct ged_view_context *view);

/** Test for an independent scene scope, optionally creating it. */
GED_EXPORT extern int ged_view_context_independent_scope_is_null(struct ged_view_context *view, int create);

/** Destroy the semantic scene scope owned only by @p view. */
GED_EXPORT extern void ged_view_context_independent_scope_destroy(struct ged_view_context *view);

/** Return non-zero if a presentation backend is attached to @p view. */
GED_EXPORT extern int ged_view_context_scene_attached(const struct ged_view_context *view);

/** Clear the view-owned state selected by @p flags and return the change count. */
GED_EXPORT extern size_t ged_view_context_clear(struct ged_view_context *view, int flags);

/** Return the borrowed GED owner of @p view, or NULL for a standalone view. */
GED_EXPORT extern struct ged *ged_view_context_owner(
	const struct ged_view_context *view);

/** Copy the effective view LoD policy into caller-owned storage. */
GED_EXPORT extern int
ged_view_lod_policy_get(ged_view_lod_policy *policy,
	const struct ged_view_context *view);

/** Sanitize and apply a complete LoD policy to a view. */
GED_EXPORT extern int
ged_view_lod_policy_apply(struct ged_view_context *view,
	const ged_view_lod_policy *policy);

/** Apply a LoD policy with an explicit automatic-mode BoT threshold. */
GED_EXPORT extern int
ged_view_lod_policy_apply_bot_threshold(struct ged_view_context *view,
	const ged_view_lod_policy *policy, size_t bot_threshold);

/** Recompute cached LoD bounds after the view geometry changes. */
GED_EXPORT extern int
ged_view_lod_bounds_update(struct ged_view_context *view);

/** Replace the borrowed callback table associated with @p view. */
GED_EXPORT extern int ged_view_context_callbacks_set(struct ged_view_context *view, struct bu_ptbl *callbacks);

/** Allocate an independent GED view context. */
GED_EXPORT extern struct ged_view_context *ged_view_context_create(void);

/** Allocate a GED view context attached to @p set. */
GED_EXPORT extern struct ged_view_context *ged_view_context_create_with_set(struct ged_view_set *set);

/** Allocate a copy of @p source attached to @p set. */
GED_EXPORT extern struct ged_view_context *ged_view_context_create_copy_with_set(const struct ged_view_context *source, struct ged_view_set *set);

/** Release a view context created by the GED view API. */
GED_EXPORT extern void ged_view_context_free(struct ged_view_context *view);

/** Attach @p view to the semantic scene host owned by @p gedp. */
GED_EXPORT extern int ged_view_context_host_attach(struct ged *gedp, struct ged_view_context *view);

/** Add @p view to @p set without transferring view ownership. */
GED_EXPORT extern int ged_view_set_context_add(struct ged_view_set *set, struct ged_view_context *view);

/** Remove @p view from @p set. */
GED_EXPORT extern int ged_view_set_context_remove(struct ged_view_set *set, struct ged_view_context *view);

/** Move @p view to @p set. */
GED_EXPORT extern int ged_view_context_view_set_attach(struct ged_view_context *view, struct ged_view_set *set);

/** Set the update callback and borrowed callback data for @p view. */
GED_EXPORT extern int ged_view_context_update_callback_set(struct ged_view_context *view, ged_view_context_update_callback_t callback, void *data);

/** Notify the registered callback that @p view changed. */
GED_EXPORT extern int ged_view_context_update(struct ged_view_context *view);

/**
 * Translate the view.
 */
GED_EXPORT extern int ged_tra_args(struct ged *gedp, int argc, const char *argv[], char *coord, vect_t tvec);



/**
 * Return ged selections for specified object. Created if it doesn't
 * exist.
 */
GED_EXPORT struct rt_object_selections *ged_get_object_selections(struct ged *gedp,
                                                                  const char *object_name);

/**
 * Return ged selections of specified kind for specified object.
 * Created if it doesn't exist.
 */
GED_EXPORT struct rt_selection_set *ged_get_selection_set(struct ged *gedp,
                                                          const char *object_name,
                                                          const char *selection_name);



typedef void (*ged_drawable_notify_func_t)(int);

/** Set the legacy display-list notification callback for @p gedp. */
GED_EXPORT void
ged_dl_notify_func_set(struct ged *gedp, ged_drawable_notify_func_t f);

/** Return the legacy display-list notification callback for @p gedp. */
GED_EXPORT ged_drawable_notify_func_t
ged_dl_notify_func_get(struct ged *gedp);




/* NMG specific visualizations that (currently) need libged routines.
 *
 * This will almost certainly move elsewhere - its presence here should be
 * considered temporary and not relied on from an API design perspective.
 */
GED_EXPORT extern void nmg_plot_eu(struct ged *gedp, struct edgeuse *es_eu, const struct bn_tol *tol);


__END_DECLS

#endif /* GED_VIEW_H */

/** @} */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
