/*           F E A T U R E _ R E S O L V E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libged/view/feature_resolver.cpp */

#include "common.h"

#include "feature_resolver.hpp"

#include "../ged_bobol_private.hpp"
#include "BObol/BSceneController.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BViewController.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "ged/view.h"

#include <Inventor/nodes/SoGroup.h>

#include <algorithm>
#include <vector>

static const char *
managed_kind_name(GedViewManagedFeatureKind kind)
{
    switch (kind) {
	case GedViewManagedFeatureKind::Feature: return "feature";
	case GedViewManagedFeatureKind::Polygon: return "polygon";
	case GedViewManagedFeatureKind::DatabaseGroup: return "database-group";
	case GedViewManagedFeatureKind::DatabaseSource: return "database-source";
	case GedViewManagedFeatureKind::None:
	default: return "none";
    }
}

static const char *
skip_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static bool
path_equal(const char *a, const char *b)
{
    return BU_STR_EQUAL(skip_slash(a), skip_slash(b));
}

bool
GedViewManagedFeatureRef::valid(void) const
{
    switch (this->kind) {
	case GedViewManagedFeatureKind::Feature:
	    return this->view && this->feature.isValid();
	case GedViewManagedFeatureKind::Polygon:
	    return this->view && this->polygon.isValid();
	case GedViewManagedFeatureKind::DatabaseGroup:
	    return this->scene && !this->groupPath.empty();
	case GedViewManagedFeatureKind::DatabaseSource:
	    return this->scene && this->source.valid;
	case GedViewManagedFeatureKind::None:
	default:
	    return false;
    }
}

bool
GedViewManagedFeatureRef::current(void) const
{
    if (!this->valid())
	return false;
    if (this->kind == GedViewManagedFeatureKind::Feature) {
	BObolFeatureRecord record;
	return this->view->features().record(this->feature, record);
    }
    if (this->kind == GedViewManagedFeatureKind::Polygon) {
	BObolPolygonRecord record;
	return this->view->polygons().record(this->polygon, record);
    }
    if (this->scene->getStructuralRevision() != this->sceneRevision)
	return false;
    if (this->kind == GedViewManagedFeatureKind::DatabaseGroup) {
	SoGroup *group = this->scene->findGroup(this->groupPath.c_str());
	return group && group->isOfType(SoBRLSceneGroup::getClassTypeId());
    }
    if (this->source.instanceKey.getLength())
	return this->scene->findDatabaseSourceInstance(
	    this->source.instanceKey.getString()) != nullptr;
    return this->scene->findDatabaseSource(
	this->source.path.getString()) != nullptr;
}

struct FeatureCollectState {
    const char *name = nullptr;
    BObolViewController *view = nullptr;
    std::vector<GedViewManagedFeatureRef> *matches = nullptr;
    unsigned int polygonScopeMask = BOBOL_FEATURE_SCOPE_ALL;
};

static int
collect_feature(const BObolFeatureRecord &record, void *data)
{
    FeatureCollectState *state = static_cast<FeatureCollectState *>(data);
    if (!state || !state->matches || !state->view || !state->name ||
	!BU_STR_EQUAL(record.name.getString(), state->name) ||
	(record.scope == BObolFeatureScope::Local &&
	 !(state->polygonScopeMask & BOBOL_FEATURE_SCOPE_LOCAL)) ||
	(record.scope == BObolFeatureScope::Shared &&
	 !(state->polygonScopeMask & BOBOL_FEATURE_SCOPE_SHARED)))
	return 1;
    GedViewManagedFeatureRef ref;
    ref.kind = GedViewManagedFeatureKind::Feature;
    ref.query = state->name;
    ref.view = state->view;
    ref.feature = record.handle;
    state->matches->push_back(ref);
    return 1;
}

static int
collect_polygon(const BObolPolygonRecord &record, void *data)
{
    FeatureCollectState *state = static_cast<FeatureCollectState *>(data);
    if (!state || !state->matches || !state->view || !state->name ||
	!BU_STR_EQUAL(record.name.getString(), state->name))
	return 1;
    GedViewManagedFeatureRef ref;
    ref.kind = GedViewManagedFeatureKind::Polygon;
    ref.query = state->name;
    ref.view = state->view;
    ref.polygon = record.handle;
    state->matches->push_back(ref);
    return 1;
}

static void
collect_view_matches(std::vector<GedViewManagedFeatureRef> &matches,
    struct ged *gedp, struct ged_view_context *view_ctx, const char *name,
    bool local_only)
{
    BObolViewController *local = ged_bobol_view_controller(view_ctx);
    FeatureCollectState state;
    state.name = name;
    state.matches = &matches;
    if (local) {
	state.view = local;
	const BObolFeatureOwner owner = ged_bobol_view_feature_owner(view_ctx);
	local->features().visitRecords(collect_feature, &state,
	    BOBOL_FEATURE_SCOPE_LOCAL, &owner);
	state.polygonScopeMask = BOBOL_FEATURE_SCOPE_LOCAL;
	local->polygons().visitRecords(collect_polygon, &state);
	if (!local_only) {
	    local->features().visitRecords(collect_feature, &state,
		BOBOL_FEATURE_SCOPE_SHARED, nullptr);
	    state.polygonScopeMask = BOBOL_FEATURE_SCOPE_SHARED;
	    local->polygons().visitRecords(collect_polygon, &state);
	}
    }
    if (local_only)
	return;

    BObolViewController *shared = ged_bobol_shared_view_controller(gedp);
    if (shared && shared != local) {
	state.view = shared;
	shared->features().visitRecords(collect_feature, &state,
	    BOBOL_FEATURE_SCOPE_SHARED, nullptr);
	state.polygonScopeMask = BOBOL_FEATURE_SCOPE_SHARED;
	shared->polygons().visitRecords(collect_polygon, &state);
    }
}

static void
collect_database_matches(std::vector<GedViewManagedFeatureRef> &matches,
    struct ged *gedp, const char *name)
{
    BObolSceneController *scene = ged_bobol_scene(gedp);
    if (!scene)
	return;
    const uint64_t revision = scene->getStructuralRevision();

    SoGroup *group = scene->findGroup(name);
    if (group && group->isOfType(SoBRLSceneGroup::getClassTypeId())) {
	GedViewManagedFeatureRef ref;
	ref.kind = GedViewManagedFeatureKind::DatabaseGroup;
	ref.query = name;
	ref.scene = scene;
	ref.sceneRevision = revision;
	ref.groupPath = name;
	matches.push_back(ref);
    }

    std::vector<BObolDatabaseSourceSummary> exact_instances;
    std::vector<BObolDatabaseSourceSummary> named;
    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
	    continue;
	if (summary.instanceKey.getLength() &&
	    BU_STR_EQUAL(summary.instanceKey.getString(), name)) {
	    exact_instances.push_back(summary);
	    continue;
	}
	if (path_equal(summary.path.getString(), name) ||
	    (summary.displayName.getLength() &&
	     BU_STR_EQUAL(summary.displayName.getString(), name)))
	    named.push_back(summary);
    }
    const std::vector<BObolDatabaseSourceSummary> &sources =
	!exact_instances.empty() ? exact_instances : named;
    for (const BObolDatabaseSourceSummary &summary : sources) {
	GedViewManagedFeatureRef ref;
	ref.kind = GedViewManagedFeatureKind::DatabaseSource;
	ref.query = name;
	ref.scene = scene;
	ref.sceneRevision = revision;
	ref.source = summary;
	matches.push_back(ref);
    }
}

static bool
same_ref(const GedViewManagedFeatureRef &a,
    const GedViewManagedFeatureRef &b)
{
    if (a.kind != b.kind)
	return false;
    if (a.kind == GedViewManagedFeatureKind::Feature)
	return a.view == b.view && a.feature == b.feature;
    if (a.kind == GedViewManagedFeatureKind::Polygon)
	return a.view == b.view && a.polygon == b.polygon;
    if (a.kind == GedViewManagedFeatureKind::DatabaseGroup)
	return a.scene == b.scene && a.groupPath == b.groupPath;
    if (a.kind == GedViewManagedFeatureKind::DatabaseSource) {
	if (a.scene != b.scene)
	    return false;
	if (a.source.instanceKey.getLength() || b.source.instanceKey.getLength())
	    return a.source.instanceKey == b.source.instanceKey;
	return a.source.path == b.source.path;
    }
    return true;
}

int
ged_view_managed_feature_resolve(struct ged *gedp,
    struct ged_view_context *view_ctx, const char *name, unsigned int domains,
    bool local_only, GedViewManagedFeatureRef &result,
    struct bu_vls *diagnostics)
{
    result = GedViewManagedFeatureRef();
    if (!gedp || !view_ctx || !name || !name[0])
	return 0;

    std::vector<GedViewManagedFeatureRef> matches;
    if (domains & GED_VIEW_MANAGED_FEATURES)
	collect_view_matches(matches, gedp, view_ctx, name, local_only);
    if (!local_only && (domains & GED_VIEW_MANAGED_DATABASE))
	collect_database_matches(matches, gedp, name);

    std::vector<GedViewManagedFeatureRef> unique;
    for (const GedViewManagedFeatureRef &candidate : matches) {
	if (!candidate.current())
	    continue;
	if (std::find_if(unique.begin(), unique.end(),
		[&candidate](const GedViewManagedFeatureRef &existing) {
		    return same_ref(candidate, existing);
		}) == unique.end())
	    unique.push_back(candidate);
    }
    if (unique.size() == 1) {
	result = unique.front();
	return 1;
    }
    if (unique.empty())
	return 0;

    if (diagnostics) {
	bu_vls_printf(diagnostics,
	    "View feature name %s is ambiguous (%zu matches:", name,
	    unique.size());
	for (const GedViewManagedFeatureRef &candidate : unique)
	    bu_vls_printf(diagnostics, " %s", managed_kind_name(candidate.kind));
	bu_vls_printf(diagnostics, ")\n");
    }
    return -1;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
