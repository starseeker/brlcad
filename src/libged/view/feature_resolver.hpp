/*           F E A T U R E _ R E S O L V E R . H P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libged/view/feature_resolver.hpp
 *
 * Command-local typed resolution for objects managed by `view feature`.
 * References borrow their owners and are valid only on the GED/Obol owner
 * thread.  Feature and polygon references validate their value handles;
 * database references validate the scene structural revision captured at
 * resolution time.
 */

#ifndef LIBGED_VIEW_FEATURE_RESOLVER_HPP
#define LIBGED_VIEW_FEATURE_RESOLVER_HPP

#include "BObol/BDatabaseSource.h"
#include "BObol/BViewStore.h"

#include <stdint.h>
#include <string>

struct bu_vls;
struct ged;
struct ged_view_context;
class BObolSceneController;
class BObolViewController;

enum class GedViewManagedFeatureKind {
    None,
    Feature,
    Polygon,
    DatabaseGroup,
    DatabaseSource
};

enum GedViewManagedFeatureDomain {
    GED_VIEW_MANAGED_FEATURES = 1,
    GED_VIEW_MANAGED_DATABASE = 2,
    GED_VIEW_MANAGED_ALL = 3
};

struct GedViewManagedFeatureRef {
    GedViewManagedFeatureKind kind = GedViewManagedFeatureKind::None;
    std::string query;
    BObolViewController *view = nullptr;
    BObolFeatureHandle feature;
    BObolPolygonHandle polygon;
    BObolSceneController *scene = nullptr;
    uint64_t sceneRevision = 0;
    std::string groupPath;
    BObolDatabaseSourceSummary source;

    bool valid(void) const;
    bool current(void) const;
};

/* Return 1 for one current match, 0 for no match, and -1 for ambiguity.
 * When diagnostics is non-NULL, ambiguity is reported with candidate kinds;
 * missing names are left to the calling command so it can choose wording. */
int ged_view_managed_feature_resolve(
    struct ged *gedp,
    struct ged_view_context *view_ctx,
    const char *name,
    unsigned int domains,
    bool local_only,
    GedViewManagedFeatureRef &result,
    struct bu_vls *diagnostics = nullptr);

#endif /* LIBGED_VIEW_FEATURE_RESOLVER_HPP */
