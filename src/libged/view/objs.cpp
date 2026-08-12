/*                        O B J S . C P P
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
/** @file libged/view/objs.c
 *
 * Commands for view features.
 *
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include <ctype.h>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <vector>

extern "C" {
#include "bu/cmd.h"
#include "bu/color.h"
#include "bu/opt.h"
#include "bu/path.h"
#include "bu/vls.h"
#include "bv.h"
#include "raytrace.h"
#include "ged/draw.h"
}
#include "../ged_bobol_private.hpp"
#include "BObol/BSceneController.h"
#include "BObol/BSceneGroup.h"
#include "feature_resolver.hpp"
#include "./ged_view.h"
#include "../ged_private.h"

struct view_obj_list_state {
    std::set<std::string> names;
    const char *glob = nullptr;
};

static int
_view_obj_list_match(const char *name, const char *glob)
{
    return name && name[0] &&
	(!glob || !glob[0] || bu_path_match(glob, name, 0) == 0);
}

static int
_view_obj_list_feature_cb(const BObolFeatureRecord &record, void *ud)
{
    struct view_obj_list_state *ctx =
	static_cast<struct view_obj_list_state *>(ud);
    if (ctx && _view_obj_list_match(record.name.getString(), ctx->glob))
	ctx->names.insert(record.name.getString());
    return 1;
}

static int
_view_obj_list_polygon_cb(const BObolPolygonRecord &record, void *ud)
{
    struct view_obj_list_state *ctx =
	static_cast<struct view_obj_list_state *>(ud);
    if (ctx && _view_obj_list_match(record.name.getString(), ctx->glob))
	ctx->names.insert(record.name.getString());
    return 1;
}

static void
_view_obj_json_string(struct bu_vls *out, const char *value)
{
    bu_vls_putc(out, '"');
    for (const unsigned char *cp =
	 reinterpret_cast<const unsigned char *>(value ? value : "");
	 *cp; cp++) {
	switch (*cp) {
	    case '"': bu_vls_strcat(out, "\\\""); break;
	    case '\\': bu_vls_strcat(out, "\\\\"); break;
	    case '\b': bu_vls_strcat(out, "\\b"); break;
	    case '\f': bu_vls_strcat(out, "\\f"); break;
	    case '\n': bu_vls_strcat(out, "\\n"); break;
	    case '\r': bu_vls_strcat(out, "\\r"); break;
	    case '\t': bu_vls_strcat(out, "\\t"); break;
	    default:
		if (*cp < 0x20)
		    bu_vls_printf(out, "\\u%04x", static_cast<unsigned int>(*cp));
		else
		    bu_vls_putc(out, static_cast<char>(*cp));
		break;
	}
    }
    bu_vls_putc(out, '"');
}

static void
_view_obj_list(struct bu_vls *out, struct ged_view_context *view_ctx,
    int list_view, int list_db, int local_only, const char *glob, int json)
{
    if (!out || !view_ctx)
	return;

    struct view_obj_list_state ctx;
    ctx.glob = glob;
    struct ged *gedp = ged_view_context_owner(view_ctx);

    if (list_view) {
	BObolViewController *local = ged_bobol_view_controller(view_ctx);
	const BObolFeatureOwner owner =
	    ged_bobol_view_feature_owner(view_ctx);
	if (local) {
	    local->features().visitRecords(_view_obj_list_feature_cb, &ctx,
		BOBOL_FEATURE_SCOPE_LOCAL, &owner);
	    local->polygons().visitRecords(_view_obj_list_polygon_cb, &ctx);
	    if (!local_only)
		local->features().visitRecords(_view_obj_list_feature_cb, &ctx,
		    BOBOL_FEATURE_SCOPE_SHARED, nullptr);
	}
	if (!local_only) {
	    BObolViewController *shared =
		ged_bobol_shared_view_controller(gedp);
	    if (shared && shared != local) {
		shared->features().visitRecords(_view_obj_list_feature_cb, &ctx,
		    BOBOL_FEATURE_SCOPE_SHARED, nullptr);
		shared->polygons().visitRecords(_view_obj_list_polygon_cb, &ctx);
	    }
	}
    }

    if (list_db && !local_only) {
	BObolSceneController *scene = ged_bobol_scene(gedp);
	for (int i = 0; scene && i < scene->getDatabaseSourceCount(); i++) {
	    BObolDatabaseSourceSummary summary;
	    if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid)
		continue;
	    const char *name = summary.path.getString();
	    if (_view_obj_list_match(name, glob))
		ctx.names.insert(name);
	}
    }

    if (!json) {
	for (std::set<std::string>::iterator it = ctx.names.begin();
	     it != ctx.names.end(); ++it)
	    bu_vls_printf(out, "%s\n", it->c_str());
	return;
    }
    bu_vls_putc(out, '[');
    bool first = true;
    for (const std::string &name : ctx.names) {
	if (!first)
	    bu_vls_putc(out, ',');
	_view_obj_json_string(out, name.c_str());
	first = false;
    }
    bu_vls_strcat(out, "]\n");
}

static const char *
_view_obj_mode_name(int mode)
{
    switch (mode) {
	case GED_DRAW_MODE_WIRE:
	    return "wireframe";
	case GED_DRAW_MODE_SHADED_BOTS:
	case GED_DRAW_MODE_SHADED:
	    return "shaded";
	case GED_DRAW_MODE_EVAL_WIRE:
	    return "evaluated_wireframe";
	case GED_DRAW_MODE_HIDDEN_LINE:
	    return "hidden_line";
	case GED_DRAW_MODE_EVAL_POINTS:
	    return "evaluated_points";
	default:
	    return "unknown";
    }
}

static void
_view_obj_mode_value_string(struct bu_vls *out, int mode)
{
    if (out)
	bu_vls_strcat(out, _view_obj_mode_name(mode));
}

struct view_obj_record {
    int found;
    GedViewManagedFeatureRef ref;
    std::string name;
    std::string type_name;
    int draw_mode;
    int visible;
    unsigned char color[3];
    size_t vlist_structure_count;
    std::string kind_name;
    std::string scope_name;
    std::string overlay_class_name;
    std::string lifecycle_name;
    std::string owner;
    std::string owner_role;
    uint64_t owner_generation;
    int transient_preview;
    int command_result;
};

static unsigned char _view_obj_color_channel(float value);

static const char *
_view_obj_feature_type_name(BObolFeatureKind kind)
{
    switch (kind) {
	case BObolFeatureKind::Lines: return "line-set";
	case BObolFeatureKind::IndexedLines: return "indexed-line-set";
	case BObolFeatureKind::Points: return "point-set";
	case BObolFeatureKind::Labels: return "labels";
	case BObolFeatureKind::Arrow: return "arrow";
	case BObolFeatureKind::Axes: return "axes";
	case BObolFeatureKind::LineLayer: return "line-layer";
	case BObolFeatureKind::EditPreview: return "edit-preview";
	case BObolFeatureKind::IndexedFaceSet: return "indexed-face-set";
	case BObolFeatureKind::PolygonOverlay: return "polygon-overlay";
	case BObolFeatureKind::HudLabel: return "hud-label";
	case BObolFeatureKind::CustomNode: return "custom-node";
	case BObolFeatureKind::Unknown:
	default: return "view-feature";
    }
}

static const char *
_view_obj_feature_kind_name(BObolFeatureKind kind)
{
    switch (kind) {
	case BObolFeatureKind::Lines: return "lines";
	case BObolFeatureKind::IndexedLines: return "indexed-lines";
	case BObolFeatureKind::Points: return "points";
	case BObolFeatureKind::Labels: return "labels";
	case BObolFeatureKind::Arrow: return "arrow";
	case BObolFeatureKind::Axes: return "axes";
	case BObolFeatureKind::LineLayer: return "line-layer";
	case BObolFeatureKind::EditPreview: return "edit-preview";
	case BObolFeatureKind::IndexedFaceSet: return "indexed-face-set";
	case BObolFeatureKind::PolygonOverlay: return "polygon-overlay";
	case BObolFeatureKind::HudLabel: return "hud-label";
	case BObolFeatureKind::CustomNode: return "custom-node";
	case BObolFeatureKind::Unknown:
	default: return "unknown";
    }
}

static const char *
_view_obj_scope_name(BObolFeatureScope scope)
{
    return scope == BObolFeatureScope::Local ? "local" : "shared";
}

static const char *
_view_obj_overlay_class_name(BObolOverlayClass overlay_class)
{
    switch (overlay_class) {
	case BObolOverlayClass::None: return "none";
	case BObolOverlayClass::Faceplate: return "faceplate";
	case BObolOverlayClass::EditHandle: return "edit-handle";
	case BObolOverlayClass::Measure: return "measure";
	case BObolOverlayClass::SelectionRubberBand: return "selection-rubber-band";
	case BObolOverlayClass::SnapGuide: return "snap-guide";
	case BObolOverlayClass::CommandResult: return "command-result";
	case BObolOverlayClass::Diagnostic: return "diagnostic";
	case BObolOverlayClass::TclOverlay: return "tcl-overlay";
	case BObolOverlayClass::PolygonEdit: return "polygon-edit";
	case BObolOverlayClass::SketchEdit: return "sketch-edit";
	case BObolOverlayClass::UserAnnotation: return "user-annotation";
	default: return "unknown";
    }
}

static const char *
_view_obj_lifecycle_name(BObolOverlayLifecycle lifecycle)
{
    switch (lifecycle) {
	case BObolOverlayLifecycle::None: return "none";
	case BObolOverlayLifecycle::Persistent: return "persistent";
	case BObolOverlayLifecycle::PerFrame: return "per-frame";
	case BObolOverlayLifecycle::PerCommand: return "per-command";
	case BObolOverlayLifecycle::PerTool: return "per-tool";
	case BObolOverlayLifecycle::PerView: return "per-view";
	case BObolOverlayLifecycle::SharedViewSet: return "shared-view-set";
	case BObolOverlayLifecycle::AutoRemoveOnSource: return "auto-remove-on-source";
	default: return "unknown";
    }
}

static size_t
_view_obj_feature_structure_count(const BObolFeatureRecord &record)
{
    switch (record.kind) {
	case BObolFeatureKind::Points:
	    return record.points.size();
	case BObolFeatureKind::Labels:
	case BObolFeatureKind::HudLabel:
	    return record.labels.size();
	case BObolFeatureKind::Axes:
	    return record.axesCenters.size();
	case BObolFeatureKind::LineLayer: {
	    size_t count = 0;
	    for (const BObolLineLayer &layer : record.layers)
		count += layer.points.size() / 2;
	    return count;
	}
	case BObolFeatureKind::IndexedFaceSet:
	    return record.indices.size() / 3;
	default:
	    if (!record.commands.empty()) {
		size_t count = 0;
		for (int32_t command : record.commands)
		    if (command != static_cast<int32_t>(BObolLineCommand::Move))
			count++;
		return count;
	    }
	    return record.points.size() / 2;
    }
}

static int
_view_obj_record_find(struct ged_view_context *view_ctx,
		      const char *name,
		      int list_view,
		      int list_db,
		      int local_only,
		      struct view_obj_record *out)
{
    if (out) {
	*out = view_obj_record();
	out->kind_name = "unknown";
	out->scope_name = "unknown";
	out->overlay_class_name = "unknown";
	out->lifecycle_name = "unknown";
    }
    if (!view_ctx || !name || !name[0] || !out)
	return 0;
    struct ged *gedp = ged_view_context_owner(view_ctx);
    unsigned int domains = 0;
    if (list_view)
	domains |= GED_VIEW_MANAGED_FEATURES;
    if (list_db)
	domains |= GED_VIEW_MANAGED_DATABASE;
    const int resolved = ged_view_managed_feature_resolve(gedp, view_ctx,
	name, domains, local_only != 0, out->ref, gedp->ged_result_str);
    if (resolved != 1)
	return resolved;

    out->found = 1;
    out->name = name;
    if (out->ref.kind == GedViewManagedFeatureKind::Feature) {
	BObolFeatureRecord record;
	if (!out->ref.view->features().record(out->ref.feature, record))
	    return 0;
	out->type_name = _view_obj_feature_type_name(record.kind);
	out->kind_name = _view_obj_feature_kind_name(record.kind);
	out->scope_name = _view_obj_scope_name(record.scope);
	if (record.overlay.isOverlay) {
	    out->overlay_class_name =
		_view_obj_overlay_class_name(record.overlay.overlayClass);
	    out->lifecycle_name =
		_view_obj_lifecycle_name(record.overlay.lifecycle);
	}
	out->owner = record.owner.ownerId.getString();
	out->owner_role = record.owner.ownerRole.getString();
	out->owner_generation = record.owner.generation;
	out->transient_preview = record.kind == BObolFeatureKind::EditPreview;
	out->command_result = record.overlay.isOverlay &&
	    record.overlay.overlayClass == BObolOverlayClass::CommandResult;
	out->draw_mode = record.kind == BObolFeatureKind::IndexedFaceSet ?
	    GED_DRAW_MODE_SHADED : GED_DRAW_MODE_WIRE;
	out->visible = !record.style.hasVisible || record.style.visible;
	const SbColor color = record.style.hasColor ? record.style.color :
	    SbColor(1.0f, 1.0f, 1.0f);
	out->color[0] = _view_obj_color_channel(color[0]);
	out->color[1] = _view_obj_color_channel(color[1]);
	out->color[2] = _view_obj_color_channel(color[2]);
	out->vlist_structure_count = _view_obj_feature_structure_count(record);
	return 1;
    }
    if (out->ref.kind == GedViewManagedFeatureKind::Polygon) {
	BObolPolygonRecord record;
	if (!out->ref.view->polygons().record(out->ref.polygon, record))
	    return 0;
	out->type_name = "polygon";
	out->kind_name = "polygon";
	out->scope_name = _view_obj_scope_name(record.scope);
	out->draw_mode = GED_DRAW_MODE_WIRE;
	out->visible = 1;
	out->color[0] = _view_obj_color_channel(record.edgeColor[0]);
	out->color[1] = _view_obj_color_channel(record.edgeColor[1]);
	out->color[2] = _view_obj_color_channel(record.edgeColor[2]);
	out->vlist_structure_count = record.contourCount;
	return 1;
    }
    if (out->ref.kind == GedViewManagedFeatureKind::DatabaseGroup) {
	SoGroup *group = out->ref.scene->findGroup(out->ref.groupPath.c_str());
	if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	    return 0;
	SoBRLSceneGroup *scene_group = static_cast<SoBRLSceneGroup *>(group);
	out->type_name = "database-group";
	out->kind_name = "database-group";
	out->scope_name = "database";
	out->draw_mode = scene_group->drawMode.getValue();
	out->visible = scene_group->visible.getValue();
	const SbColor color = scene_group->colorOverride.getValue() ?
	    scene_group->color.getValue() :
	    (scene_group->materialColorValid.getValue() ?
	     scene_group->materialColor.getValue() :
	     scene_group->color.getValue());
	out->color[0] = _view_obj_color_channel(color[0]);
	out->color[1] = _view_obj_color_channel(color[1]);
	out->color[2] = _view_obj_color_channel(color[2]);
	out->vlist_structure_count = static_cast<size_t>(
	    out->ref.scene->getGroupDatabaseSourceCount(
		out->ref.groupPath.c_str()));
	return 1;
    }
    if (out->ref.kind == GedViewManagedFeatureKind::DatabaseSource) {
	const BObolDatabaseSourceSummary &summary = out->ref.source;
	out->type_name = "database-source";
	out->kind_name = "database-source";
	out->scope_name = "database";
	out->draw_mode = summary.drawMode;
	out->visible = summary.visible;
	const SbColor color = summary.colorOverride ? summary.color :
	    (summary.materialColorValid ? summary.materialColor : summary.color);
	out->color[0] = _view_obj_color_channel(color[0]);
	out->color[1] = _view_obj_color_channel(color[1]);
	out->color[2] = _view_obj_color_channel(color[2]);
	out->vlist_structure_count = static_cast<size_t>(
	    summary.realizedShapeCount);
	return 1;
    }
    return 0;
}

static unsigned char
_view_obj_color_channel(float value)
{
    if (value <= 0.0f)
	return 0;
    if (value >= 1.0f)
	return 255;
    return static_cast<unsigned char>(value * 255.0f + 0.5f);
}

static void
_view_obj_record_json(struct bu_vls *out, const struct view_obj_record &rec)
{
    bu_vls_putc(out, '{');
    bu_vls_strcat(out, "\"name\":");
    _view_obj_json_string(out, rec.name.c_str());
    bu_vls_strcat(out, ",\"type\":");
    _view_obj_json_string(out, rec.type_name.c_str());
    bu_vls_strcat(out, ",\"kind\":");
    _view_obj_json_string(out, rec.kind_name.c_str());
    bu_vls_strcat(out, ",\"scope\":");
    _view_obj_json_string(out, rec.scope_name.c_str());
    bu_vls_strcat(out, ",\"mode\":");
    _view_obj_json_string(out, _view_obj_mode_name(rec.draw_mode));
    bu_vls_printf(out,
	",\"visible\":%s,\"color\":[%u,%u,%u],\"structure_count\":%zu",
	rec.visible ? "true" : "false",
	static_cast<unsigned int>(rec.color[0]),
	static_cast<unsigned int>(rec.color[1]),
	static_cast<unsigned int>(rec.color[2]), rec.vlist_structure_count);
    bu_vls_strcat(out, ",\"overlay_class\":");
    _view_obj_json_string(out, rec.overlay_class_name.c_str());
    bu_vls_strcat(out, ",\"lifecycle\":");
    _view_obj_json_string(out, rec.lifecycle_name.c_str());
    bu_vls_strcat(out, ",\"owner\":");
    _view_obj_json_string(out, rec.owner.c_str());
    bu_vls_strcat(out, ",\"owner_role\":");
    _view_obj_json_string(out, rec.owner_role.c_str());
    bu_vls_printf(out,
	",\"owner_generation\":%llu,\"transient_preview\":%s,"
	"\"command_result\":%s}\n",
	static_cast<unsigned long long>(rec.owner_generation),
	rec.transient_preview ? "true" : "false",
	rec.command_result ? "true" : "false");
}

static void
_view_obj_style_to_ged(struct ged_view_feature_style *dst,
    const BObolFeatureStyle &src)
{
    struct ged_view_feature_style init = ged_view_feature_style_default();
    *dst = init;
    if (src.hasVisible)
	dst->visible = src.visible ? 1 : 0;
    if (src.hasSelectable)
	dst->selectable = src.selectable ? 1 : 0;
    if (src.hasColor) {
	dst->color_valid = 1;
	dst->color[0] = _view_obj_color_channel(src.color[0]);
	dst->color[1] = _view_obj_color_channel(src.color[1]);
	dst->color[2] = _view_obj_color_channel(src.color[2]);
    }
    if (src.hasLineWidth)
	dst->line_width = src.lineWidth;
    if (src.hasLineStyle)
	dst->line_style = src.lineStyle;
    if (src.hasArrow)
	dst->arrow = src.arrow ? 1 : 0;
    if (src.hasArrowTip) {
	dst->arrow_tip_length = src.arrowTipLength;
	dst->arrow_tip_width = src.arrowTipWidth;
    }
}

static BObolFeatureStyle
_view_obj_style_from_ged(const struct ged_view_feature_style *src)
{
    BObolFeatureStyle dst;
    if (!src)
	return dst;
    if (src->visible >= 0) {
	dst.hasVisible = TRUE;
	dst.visible = src->visible ? TRUE : FALSE;
    }
    if (src->selectable >= 0) {
	dst.hasSelectable = TRUE;
	dst.selectable = src->selectable ? TRUE : FALSE;
    }
    if (src->color_valid) {
	dst.hasColor = TRUE;
	dst.color = SbColor(static_cast<float>(src->color[0]) / 255.0f,
	    static_cast<float>(src->color[1]) / 255.0f,
	    static_cast<float>(src->color[2]) / 255.0f);
    }
    if (src->line_width >= 0) {
	dst.hasLineWidth = TRUE;
	dst.lineWidth = src->line_width;
    }
    if (src->line_style >= 0) {
	dst.hasLineStyle = TRUE;
	dst.lineStyle = src->line_style;
    }
    if (src->arrow >= 0) {
	dst.hasArrow = TRUE;
	dst.arrow = src->arrow ? TRUE : FALSE;
    }
    if (src->arrow_tip_length >= 0.0 || src->arrow_tip_width >= 0.0) {
	dst.hasArrowTip = TRUE;
	dst.arrowTipLength = src->arrow_tip_length >= 0.0 ?
	    static_cast<float>(src->arrow_tip_length) : 0.0f;
	dst.arrowTipWidth = src->arrow_tip_width >= 0.0 ?
	    static_cast<float>(src->arrow_tip_width) : 0.0f;
    }
    return dst;
}

static int
_view_obj_resolve(struct _ged_view_info *gd,
    GedViewManagedFeatureRef &ref, unsigned int domains = GED_VIEW_MANAGED_ALL)
{
    if (!gd || !gd->gedp || !gd->cv || !gd->vobj || !gd->vobj[0])
	return 0;
    const int ret = ged_view_managed_feature_resolve(gd->gedp, gd->cv,
	gd->vobj, domains, gd->local_obj != 0, ref,
	gd->gedp->ged_result_str);
    if (ret == 0)
	bu_vls_printf(gd->gedp->ged_result_str,
	    "No view feature named %s\n", gd->vobj);
    return ret;
}

static const char *
_view_obj_path_skip_slash(const char *path)
{
    if (!path)
	return "";
    while (*path == '/')
	path++;
    return path;
}

static bool
_view_obj_path_contains(const char *parent, const char *child)
{
    const char *p = _view_obj_path_skip_slash(parent);
    const char *c = _view_obj_path_skip_slash(child);
    const size_t plen = strlen(p);
    return BU_STR_EQUAL(p, c) || (plen && strlen(c) > plen &&
	bu_strncmp(p, c, plen) == 0 && c[plen] == '/');
}

static SoBRLSceneGroup *
_view_obj_database_group(const GedViewManagedFeatureRef &ref,
	BObolSceneController **scene_out = nullptr)
{
    if (scene_out)
	*scene_out = nullptr;
    if (ref.kind != GedViewManagedFeatureKind::DatabaseGroup ||
	!ref.current())
	return nullptr;

    BObolSceneController *scene = ref.scene;
    SoGroup *group = scene ? scene->findGroup(ref.groupPath.c_str()) : nullptr;
    if (!group || !group->isOfType(SoBRLSceneGroup::getClassTypeId()))
	return nullptr;
    if (scene_out)
	*scene_out = scene;
    return static_cast<SoBRLSceneGroup *>(group);
}

static int
_view_obj_database_patch(BObolSceneController *scene,
	const BObolDatabaseSourceSummary &summary,
	const BObolDatabaseSourceDisplayPatch &patch)
{
    if (!scene)
	return 0;
    if (summary.instanceKey.getLength())
	return scene->setDatabaseSourceInstanceDisplayPatch(
	    summary.instanceKey.getString(), patch);
    return scene->setDatabaseSourceDisplayPatch(summary.path.getString(),
	patch);
}

static bool
_view_obj_database_style_get(const GedViewManagedFeatureRef &ref,
	struct ged_view_feature_style *style)
{
    SoBRLSceneGroup *group = _view_obj_database_group(ref);
    if (group) {
	const struct ged_view_feature_style init = ged_view_feature_style_default();
	*style = init;
	style->visible = group->visible.getValue() ? 1 : 0;
	style->color_valid = 1;
	const SbColor color = group->colorOverride.getValue() ?
	    group->color.getValue() :
	    (group->materialColorValid.getValue() ?
	     group->materialColor.getValue() : group->color.getValue());
	style->color[0] = _view_obj_color_channel(color[0]);
	style->color[1] = _view_obj_color_channel(color[1]);
	style->color[2] = _view_obj_color_channel(color[2]);
	style->line_style = group->lineStyle.getValue();
	style->line_width = group->lineWidth.getValue();
	style->arrow = 0;
	return true;
    }

    if (!style || ref.kind != GedViewManagedFeatureKind::DatabaseSource ||
	!ref.current())
	return false;
    const BObolDatabaseSourceSummary &summary = ref.source;

    const struct ged_view_feature_style init = ged_view_feature_style_default();
    *style = init;
    style->visible = summary.visible ? 1 : 0;
    style->color_valid = 1;
    const SbColor color = summary.colorOverride ? summary.color :
	(summary.materialColorValid ? summary.materialColor : summary.color);
    style->color[0] = _view_obj_color_channel(color[0]);
    style->color[1] = _view_obj_color_channel(color[1]);
    style->color[2] = _view_obj_color_channel(color[2]);
    style->line_style = summary.lineStyle;
    style->line_width = summary.lineWidth;
    style->arrow = 0;
    return true;
}

static bool
_view_obj_database_style_apply(const GedViewManagedFeatureRef &ref,
	const struct ged_view_feature_style *style, bool recursive)
{
    if (!style)
	return false;
    if (style->arrow >= 0 || style->arrow_tip_length >= 0.0 ||
	style->arrow_tip_width >= 0.0) {
	return false;
    }

    BObolSceneController *group_scene = nullptr;
    SoBRLSceneGroup *group = _view_obj_database_group(ref, &group_scene);
    if (group) {
	const SbBool visible = style->visible >= 0 ?
	    (style->visible ? TRUE : FALSE) : group->visible.getValue();
	const int line_style = style->line_style >= 0 ?
	    style->line_style : group->lineStyle.getValue();
	const int line_width = style->line_width >= 0 ?
	    style->line_width : group->lineWidth.getValue();
	const SbBool color_override = style->color_valid ? TRUE :
	    group->colorOverride.getValue();
	const SbColor color = style->color_valid ?
	    SbColor(static_cast<float>(style->color[0]) / 255.0f,
		static_cast<float>(style->color[1]) / 255.0f,
		static_cast<float>(style->color[2]) / 255.0f) :
	    group->color.getValue();
	const int changed = group_scene->setGroupDisplayState(
	    ref.groupPath.c_str(),
	    visible, group->selected.getValue(), group->highlighted.getValue(),
	    line_style, line_width, group->transparency.getValue(),
	    color_override, color, group->materialColorValid.getValue(),
	    group->materialColor.getValue(), group->materialRevision.getValue());
	(void)recursive;
	return changed >= 0;
    }

    if (ref.kind != GedViewManagedFeatureKind::DatabaseSource ||
	!ref.current())
	return false;
    const BObolDatabaseSourceSummary &target = ref.source;
    BObolSceneController *scene = ref.scene;

    BObolDatabaseSourceDisplayPatch patch;
    if (style->visible >= 0) {
	patch.visibleValid = TRUE;
	patch.visible = style->visible ? TRUE : FALSE;
    }
    if (style->color_valid) {
	patch.colorOverrideValid = TRUE;
	patch.colorOverride = TRUE;
	patch.colorValid = TRUE;
	patch.color = SbColor(static_cast<float>(style->color[0]) / 255.0f,
	    static_cast<float>(style->color[1]) / 255.0f,
	    static_cast<float>(style->color[2]) / 255.0f);
    }
    if (style->line_style >= 0) {
	patch.lineStyleValid = TRUE;
	patch.lineStyle = style->line_style;
    }
    if (style->line_width >= 0) {
	patch.lineWidthValid = TRUE;
	patch.lineWidth = style->line_width;
    }

    if (!recursive)
	return _view_obj_database_patch(scene, target, patch) >= 0;

    bool changed = false;
    for (int i = 0; i < scene->getDatabaseSourceCount(); i++) {
	BObolDatabaseSourceSummary summary;
	if (!scene->getDatabaseSourceSummary(i, summary) || !summary.valid ||
	    !_view_obj_path_contains(target.path.getString(),
		summary.path.getString()))
	    continue;
	changed = _view_obj_database_patch(scene, summary, patch) >= 0 || changed;
    }
    return changed;
}

static bool
_view_obj_database_remove(const GedViewManagedFeatureRef &ref)
{
    BObolSceneController *group_scene = nullptr;
    if (_view_obj_database_group(ref, &group_scene))
	return group_scene->removeGroup(ref.groupPath.c_str()) > 0;

    if (ref.kind != GedViewManagedFeatureKind::DatabaseSource ||
	!ref.current())
	return false;
    const BObolDatabaseSourceSummary &summary = ref.source;
    BObolSceneController *scene = ref.scene;
    if (summary.instanceKey.getLength())
	return scene->removeDatabaseSourceInstance(
	    summary.instanceKey.getString()) > 0;
    return scene->removeDatabaseSource(summary.path.getString()) > 0;
}

static bool
_view_obj_database_realize(const GedViewManagedFeatureRef &ref)
{
    BObolSceneController *group_scene = nullptr;
    if (_view_obj_database_group(ref, &group_scene)) {
	bool found = false;
	bool realized = true;
	for (int i = 0; i < group_scene->getDatabaseSourceCount(); i++) {
	    BObolDatabaseSourceSummary summary;
	    if (!group_scene->getDatabaseSourceSummary(i, summary) ||
		!summary.valid ||
		!_view_obj_path_contains(ref.groupPath.c_str(),
		    summary.parentGroupPath.getString()))
		continue;
	    found = true;
	    SoBRLDatabaseSource *source = summary.instanceKey.getLength() ?
		group_scene->findDatabaseSourceInstance(
		    summary.instanceKey.getString()) :
		group_scene->findDatabaseSource(summary.path.getString());
	    if (!source) {
		realized = false;
		continue;
	    }
	    if (source->needsRealization())
		(void)group_scene->realizePending();
	    BObolDatabaseSourceSummary updated;
	    realized = source->getSummary(updated) && updated.valid &&
		updated.realizationStatus == SoBRLDatabaseSource::REALIZED &&
		!source->needsRealization() && realized;
	}
	return found && realized;
    }

    if (ref.kind != GedViewManagedFeatureKind::DatabaseSource ||
	!ref.current())
	return false;
    const BObolDatabaseSourceSummary &summary = ref.source;
    BObolSceneController *scene = ref.scene;
    SoBRLDatabaseSource *source = summary.instanceKey.getLength() ?
	scene->findDatabaseSourceInstance(summary.instanceKey.getString()) :
	scene->findDatabaseSource(summary.path.getString());
    if (!source)
	return false;
    if (source->needsRealization())
	(void)scene->realizePending();
    BObolDatabaseSourceSummary updated;
    return source->getSummary(updated) && updated.valid &&
	updated.realizationStatus == SoBRLDatabaseSource::REALIZED &&
	!source->needsRealization();
}

static bool
_view_obj_feature_style_get(const GedViewManagedFeatureRef &ref,
    struct ged_view_feature_style *style)
{
    BObolFeatureStyle obol_style;
    if (ref.kind != GedViewManagedFeatureKind::Feature || !ref.current() ||
	!style || !ref.view->features().style(ref.feature, obol_style))
	return false;
    _view_obj_style_to_ged(style, obol_style);
    if (style->arrow < 0)
	style->arrow = 0;
    return true;
}

static bool
_view_obj_feature_style_apply(const GedViewManagedFeatureRef &ref,
    const struct ged_view_feature_style *style, bool recursive)
{
    return ref.kind == GedViewManagedFeatureKind::Feature && ref.current() &&
	style && ref.view->features().applyStyle(ref.feature,
	    _view_obj_style_from_ged(style), recursive ? TRUE : FALSE);
}

int
_objs_cmd_draw(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature style set <name> visible [0|1|UP|DOWN]";
    const char *purpose_string = "toggle view feature visibility";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    GedViewManagedFeatureRef ref;
    if (_view_obj_resolve(gd, ref) != 1)
	return BRLCAD_ERROR;

    if (argc == 0) {
	int visible = 0;
	if (ref.kind == GedViewManagedFeatureKind::Feature) {
	    BObolFeatureStyle style;
	    if (!ref.view->features().style(ref.feature, style))
		return BRLCAD_ERROR;
	    visible = !style.hasVisible || style.visible ? 1 : 0;
	} else if (ref.kind == GedViewManagedFeatureKind::Polygon) {
	    SbBool polygon_visible = FALSE;
	    if (!ref.view->polygons().isVisible(ref.polygon, polygon_visible))
		return BRLCAD_ERROR;
	    visible = polygon_visible ? 1 : 0;
	} else {
	    struct ged_view_feature_style style =
		ged_view_feature_style_default();
	    if (!_view_obj_database_style_get(ref, &style))
		return BRLCAD_ERROR;
	    visible = style.visible > 0 ? 1 : 0;
	}
	bu_vls_printf(gedp->ged_result_str, "%s\n", visible ? "UP" : "DOWN");
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "DOWN")) {
	if (ref.kind == GedViewManagedFeatureKind::Feature)
	    return ref.view->features().setVisible(ref.feature,
		FALSE) ? BRLCAD_OK : BRLCAD_ERROR;
	if (ref.kind == GedViewManagedFeatureKind::Polygon)
	    return ref.view->polygons().setVisible(ref.polygon,
		FALSE) ? BRLCAD_OK : BRLCAD_ERROR;
	struct ged_view_feature_style style = ged_view_feature_style_default();
	style.visible = 0;
	return _view_obj_database_style_apply(ref, &style, false) ?
	    BRLCAD_OK : BRLCAD_ERROR;
    }
    if (BU_STR_EQUAL(argv[0], "UP")) {
	if (ref.kind == GedViewManagedFeatureKind::Feature)
	    return ref.view->features().setVisible(ref.feature,
		TRUE) ? BRLCAD_OK : BRLCAD_ERROR;
	if (ref.kind == GedViewManagedFeatureKind::Polygon)
	    return ref.view->polygons().setVisible(ref.polygon,
		TRUE) ? BRLCAD_OK : BRLCAD_ERROR;
	struct ged_view_feature_style style = ged_view_feature_style_default();
	style.visible = 1;
	return _view_obj_database_style_apply(ref, &style, false) ?
	    BRLCAD_OK : BRLCAD_ERROR;
    }

    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
    return BRLCAD_ERROR;
}

int
_objs_cmd_delete(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature delete <name>";
    const char *purpose_string = "delete view feature";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    GedViewManagedFeatureRef ref;
    if (_view_obj_resolve(gd, ref) != 1)
	return BRLCAD_ERROR;
    if (ref.kind == GedViewManagedFeatureKind::Feature)
	return ref.view->features().remove(ref.feature) ?
	    BRLCAD_OK : BRLCAD_ERROR;
    if (ref.kind == GedViewManagedFeatureKind::Polygon)
	return ref.view->polygons().remove(ref.polygon) ?
	    BRLCAD_OK : BRLCAD_ERROR;
    if (_view_obj_database_remove(ref))
	return BRLCAD_OK;
    return BRLCAD_ERROR;
}

int
_objs_cmd_color(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature style <get|set> <name> color [r/g/b]";
    const char *purpose_string = "show/set feature color";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    int recurse = 0;

    struct bu_opt_desc d[2];
    BU_OPT(d[0], "r", "recursive",       "",  NULL,  &recurse,  "Report (or set) color of all child objects");
    BU_OPT_NULL(d[1]);

    int ac = bu_opt_parse(NULL, argc, argv, d);

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    GedViewManagedFeatureRef ref;
    if (_view_obj_resolve(gd, ref) != 1)
	return BRLCAD_ERROR;

    struct ged_view_feature_style style =
	ged_view_feature_style_default();

    if (ac == 0) {
	if (!_view_obj_feature_style_get(ref, &style) &&
	    !_view_obj_database_style_get(ref, &style))
	    return BRLCAD_ERROR;
	bu_vls_printf(gedp->ged_result_str, "%d/%d/%d\n",
		style.color[0], style.color[1], style.color[2]);
	return BRLCAD_OK;
    }

    struct bu_color val;
    if (bu_opt_color(NULL, 1, (const char **)&argv[0], (void *)&val) != 1) {
	bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	return BRLCAD_ERROR;
    }

    bu_color_to_rgb_chars(&val, style.color);
    style.color_valid = 1;
    if (ref.kind == GedViewManagedFeatureKind::Feature)
	return _view_obj_feature_style_apply(ref, &style, recurse != 0) ?
	    BRLCAD_OK : BRLCAD_ERROR;
    return _view_obj_database_style_apply(ref, &style, recurse != 0) ?
	BRLCAD_OK : BRLCAD_ERROR;
}

int
_objs_cmd_arrow(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature style <get|set> <name> arrow [0|1] [width [#]] [length [#]]";
    const char *purpose_string = "toggle arrow drawing, for features that support it";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    GedViewManagedFeatureRef ref;
    if (_view_obj_resolve(gd, ref) != 1)
	return BRLCAD_ERROR;

    struct ged_view_feature_style style =
	ged_view_feature_style_default();
    if (!_view_obj_feature_style_get(ref, &style) &&
	!_view_obj_database_style_get(ref, &style))
	return BRLCAD_ERROR;

    if (style.arrow < 0 && style.arrow_tip_width < 0.0 &&
	    style.arrow_tip_length < 0.0) {
	bu_vls_printf(gedp->ged_result_str,
		"View feature %s does not support arrow settings\n", gd->vobj);
	return BRLCAD_ERROR;
    }

    if (argc == 0) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", style.arrow > 0 ? 1 : 0);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(argv[0], "0") || BU_STR_EQUAL(argv[0], "1")) {
	style.arrow = BU_STR_EQUAL(argv[0], "1") ? 1 : 0;
	if (ref.kind == GedViewManagedFeatureKind::Feature)
	    return _view_obj_feature_style_apply(ref, &style, false) ?
		BRLCAD_OK : BRLCAD_ERROR;
	return _view_obj_database_style_apply(ref, &style, false) ?
	    BRLCAD_OK : BRLCAD_ERROR;
    }
    if (BU_STR_EQUAL(argv[0], "width"))  {
	if (argc == 2) {
	    fastf_t width = 0.0;
	    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1], (void *)&width) != 1) {
		bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
		return BRLCAD_ERROR;
	    }
	    style.arrow_tip_width = width;
	    if (ref.kind == GedViewManagedFeatureKind::Feature)
		return _view_obj_feature_style_apply(ref, &style, false) ?
		    BRLCAD_OK : BRLCAD_ERROR;
	    return _view_obj_database_style_apply(ref, &style, false) ?
		BRLCAD_OK : BRLCAD_ERROR;
	}
	bu_vls_printf(gedp->ged_result_str, "%f\n", style.arrow_tip_width);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(argv[0], "length"))  {
	if (argc == 2) {
	    fastf_t length = 0.0;
	    if (bu_opt_fastf_t(NULL, 1, (const char **)&argv[1], (void *)&length) != 1) {
		bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
		return BRLCAD_ERROR;
	    }
	    style.arrow_tip_length = length;
	    if (ref.kind == GedViewManagedFeatureKind::Feature)
		return _view_obj_feature_style_apply(ref, &style, false) ?
		    BRLCAD_OK : BRLCAD_ERROR;
	    return _view_obj_database_style_apply(ref, &style, false) ?
		BRLCAD_OK : BRLCAD_ERROR;
	}
	bu_vls_printf(gedp->ged_result_str, "%f\n", style.arrow_tip_length);
	return BRLCAD_OK;
    }

    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
    return BRLCAD_ERROR;
}

int
_objs_cmd_lcnt(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature info <name> lcnt";
    const char *purpose_string = "print the number of vlist entities";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    struct view_obj_record rec;
    const int found = _view_obj_record_find(gd->cv, gd->vobj, 1, 1,
	gd->local_obj, &rec);
    if (found == 1) {
	bu_vls_printf(gedp->ged_result_str, "%zu\n", rec.vlist_structure_count);
	return BRLCAD_OK;
    }
    if (found < 0)
	return BRLCAD_ERROR;
    bu_vls_printf(gedp->ged_result_str, "0\n");
    return BRLCAD_OK;
}

int
_objs_cmd_update(void *bs, int argc, const char **argv)
{
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;
    const char *usage_string = "view feature realize <name> [x y]";
    const char *purpose_string = "update feature";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;

    argc--; argv++;

    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    if (argc && (argc != 2)) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    int have_xy = 0;
    int x = 0;
    int y = 0;
    if (argc) {
	if (bu_opt_int(NULL, 1, (const char **)&argv[0], (void *)&x) != 1 || x < 0) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[0]);
	    return BRLCAD_ERROR;
	}
	if (bu_opt_int(NULL, 1, (const char **)&argv[1], (void *)&y) != 1 || y < 0) {
	    bu_vls_printf(gedp->ged_result_str, "Invalid argument %s\n", argv[1]);
	    return BRLCAD_ERROR;
	}
	have_xy = 1;
	bv_mouse_state_set(bv_context_view(ged_view_context_bv(gd->cv)), x, y);
    }

    GedViewManagedFeatureRef ref;
    if (_view_obj_resolve(gd, ref) != 1)
	return BRLCAD_ERROR;
    if (ref.kind == GedViewManagedFeatureKind::Polygon) {
	if (have_xy) {
	    point_t model_point = VINIT_ZERO;
	    if (!bv_screen_to_model(model_point,
		    bv_context_view_const(ged_view_context_bv_const(gd->cv)),
		    static_cast<fastf_t>(x), static_cast<fastf_t>(y)) ||
		!ref.view->polygons().updateModelPoint(
		    ref.polygon,
		    SbVec3f(static_cast<float>(model_point[X]),
			static_cast<float>(model_point[Y]),
			static_cast<float>(model_point[Z])),
		    BObolPolygonUpdate::Default))
		return BRLCAD_ERROR;
	} else if (!ref.view->polygons().update(ref.polygon,
		BObolPolygonUpdate::Default)) {
	    return BRLCAD_ERROR;
	}
	return BRLCAD_OK;
    }

    if (ref.kind == GedViewManagedFeatureKind::Feature) {
	(void)ref.view->features().realize(ref.feature, TRUE);
	return BRLCAD_OK;
    }

    if (_view_obj_database_realize(ref))
	return BRLCAD_OK;
    return BRLCAD_ERROR;
}

static int
_view_object_scene_ready(struct _ged_view_info *gd)
{
    if (!gd || !gd->gedp)
	return 0;
    struct ged *gedp = gd->gedp;
    if (!gd->cv) {
	bu_vls_printf(gedp->ged_result_str, ": no view current in GED");
	return 0;
    }
    if (!ged_view_context_obol_endpoint_get(gd->cv) &&
	!ged_view_context_display_endpoint_ensure(gd->cv)) {
	bu_vls_printf(gedp->ged_result_str,
	    ": unable to create display endpoint for current view");
	return 0;
    }
    if (!ged_view_context_scene_attached(gd->cv)) {
	struct ged_view_context *cv = ged_view_active_ctx(gedp);
	ged_view_active_ctx_set(gedp, gd->cv);
	ged_draw_ensure_root_attached(gedp);
	ged_view_active_ctx_set(gedp, cv);
    }
    return 1;
}

static int
_view_object_first_pos(int argc, const char **argv)
{
    for (int i = 0; i < argc; i++) {
	if (!argv[i])
	    return -1;
	if (argv[i][0] == '-')
	    continue;
	return i;
    }
    return -1;
}

static void
_view_object_list_defaults(int annotations_requested, int db_requested,
	int all_requested, int *list_view, int *list_db)
{
    if (all_requested || (!annotations_requested && !db_requested)) {
	*list_view = 1;
	*list_db = 1;
	return;
    }
    *list_view = annotations_requested ? 1 : 0;
    *list_db = db_requested ? 1 : 0;
}

static int
_view_object_info(struct _ged_view_info *gd,
	const char *name,
	const char *field,
	int list_view,
	int list_db,
	int json)
{
    struct ged *gedp = gd->gedp;
    struct view_obj_record rec;
    const int found = _view_obj_record_find(gd->cv, name, list_view, list_db,
	gd->local_obj, &rec);
    if (found != 1) {
	if (found == 0)
	bu_vls_printf(gedp->ged_result_str, "No view feature named %s\n", name);
	return BRLCAD_ERROR;
    }

    if (json) {
	_view_obj_record_json(gedp->ged_result_str, rec);
	return BRLCAD_OK;
    }

    if (!field) {
	bu_vls_printf(gedp->ged_result_str, "%s %s\n", name, rec.type_name.c_str());
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(field, "mode")) {
	_view_obj_mode_value_string(gedp->ged_result_str, rec.draw_mode);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "color")) {
	bu_vls_printf(gedp->ged_result_str, "%d/%d/%d\n",
		rec.color[0], rec.color[1], rec.color[2]);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "visible") || BU_STR_EQUAL(field, "draw")) {
	bu_vls_printf(gedp->ged_result_str, "%s\n", rec.visible ? "UP" : "DOWN");
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "lcnt")) {
	bu_vls_printf(gedp->ged_result_str, "%zu\n", rec.vlist_structure_count);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "type") || BU_STR_EQUAL(field, "class")) {
	bu_vls_printf(gedp->ged_result_str, "%s\n", rec.type_name.c_str());
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "kind")) {
	bu_vls_printf(gedp->ged_result_str, "%s\n", rec.kind_name.c_str());
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "scope")) {
	bu_vls_printf(gedp->ged_result_str, "%s\n", rec.scope_name.c_str());
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "overlay_class")) {
	bu_vls_printf(gedp->ged_result_str, "%s\n",
	    rec.overlay_class_name.c_str());
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "lifecycle")) {
	bu_vls_printf(gedp->ged_result_str, "%s\n",
	    rec.lifecycle_name.c_str());
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "owner")) {
	bu_vls_printf(gedp->ged_result_str, "%s\n", rec.owner.c_str());
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "owner_role")) {
	bu_vls_printf(gedp->ged_result_str, "%s\n", rec.owner_role.c_str());
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "owner_generation")) {
	bu_vls_printf(gedp->ged_result_str, "%llu\n",
		(unsigned long long)rec.owner_generation);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "transient_preview")) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", rec.transient_preview);
	return BRLCAD_OK;
    }
    if (BU_STR_EQUAL(field, "command_result")) {
	bu_vls_printf(gedp->ged_result_str, "%d\n", rec.command_result);
	return BRLCAD_OK;
    }

    bu_vls_printf(gedp->ged_result_str, "Unsupported object info field %s", field);
    return BRLCAD_ERROR;
}

static int
_view_object_style_get(struct _ged_view_info *gd, const char *field)
{
    if (!field) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"Usage: view feature style get <name> <color|visible|arrow|arrow_width|arrow_length>");
	return BRLCAD_ERROR;
    }
    if (BU_STR_EQUAL(field, "color")) {
	const char *cargv[2] = {"color", NULL};
	return _objs_cmd_color(gd, 1, cargv);
    }
    if (BU_STR_EQUAL(field, "visible") || BU_STR_EQUAL(field, "draw")) {
	const char *dargv[2] = {"visible", NULL};
	return _objs_cmd_draw(gd, 1, dargv);
    }
    if (BU_STR_EQUAL(field, "arrow")) {
	const char *aargv[2] = {"arrow", NULL};
	return _objs_cmd_arrow(gd, 1, aargv);
    }
    if (BU_STR_EQUAL(field, "arrow_width")) {
	const char *aargv[3] = {"arrow", "width", NULL};
	return _objs_cmd_arrow(gd, 2, aargv);
    }
    if (BU_STR_EQUAL(field, "arrow_length")) {
	const char *aargv[3] = {"arrow", "length", NULL};
	return _objs_cmd_arrow(gd, 2, aargv);
    }
    bu_vls_printf(gd->gedp->ged_result_str, "Unsupported style field %s", field);
    return BRLCAD_ERROR;
}

static int
_view_object_style_set(struct _ged_view_info *gd, int argc, const char **argv)
{
    if (argc < 1) {
	bu_vls_printf(gd->gedp->ged_result_str,
		"Usage: view feature style set <name> <field> <value>");
	return BRLCAD_ERROR;
    }
    if (BU_STR_EQUAL(argv[0], "color")) {
	return _objs_cmd_color(gd, argc, argv);
    }
    if (BU_STR_EQUAL(argv[0], "visible") || BU_STR_EQUAL(argv[0], "draw")) {
	if (argc != 2) {
	    bu_vls_printf(gd->gedp->ged_result_str,
		    "Usage: view feature style set <name> visible <0|1|UP|DOWN>");
	    return BRLCAD_ERROR;
	}
	const char *state = argv[1];
	if (BU_STR_EQUAL(state, "1"))
	    state = "UP";
	else if (BU_STR_EQUAL(state, "0"))
	    state = "DOWN";
	const char *dargv[3] = {"visible", state, NULL};
	return _objs_cmd_draw(gd, 2, dargv);
    }
    if (BU_STR_EQUAL(argv[0], "arrow")) {
	return _objs_cmd_arrow(gd, argc, argv);
    }
    if (BU_STR_EQUAL(argv[0], "arrow_width")) {
	if (argc != 2) {
	    bu_vls_printf(gd->gedp->ged_result_str,
		    "Usage: view feature style set <name> arrow_width <value>");
	    return BRLCAD_ERROR;
	}
	const char *aargv[4] = {"arrow", "width", argv[1], NULL};
	return _objs_cmd_arrow(gd, 3, aargv);
    }
    if (BU_STR_EQUAL(argv[0], "arrow_length")) {
	if (argc != 2) {
	    bu_vls_printf(gd->gedp->ged_result_str,
		    "Usage: view feature style set <name> arrow_length <value>");
	    return BRLCAD_ERROR;
	}
	const char *aargv[4] = {"arrow", "length", argv[1], NULL};
	return _objs_cmd_arrow(gd, 3, aargv);
    }
    bu_vls_printf(gd->gedp->ged_result_str, "Unsupported style field %s", argv[0]);
    return BRLCAD_ERROR;
}

extern "C" int
_view_cmd_feature(void *bs, int argc, const char **argv)
{
    int help = 0;
    int local_only = 0;
    int annotations_requested = 0;
    int db_requested = 0;
    int all_requested = 0;
    int json = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view [options] feature [options] <list|info|show|hide|delete|style|realize> ...";
    const char *purpose_string = "manage typed Obol view features";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;
    if (!_view_object_scene_ready(gd))
	return BRLCAD_ERROR;

    argc--; argv++;

    struct bu_opt_desc d[7];
    BU_OPT(d[0], "h", "help",        "",  NULL,  &help,                  "Print help");
    BU_OPT(d[1], "L", "local",       "",  NULL,  &local_only,            "Restrict feature lookup to current/specified view");
    BU_OPT(d[2], "A", "annotations", "",  NULL,  &annotations_requested, "List or query annotation features");
    BU_OPT(d[3], "D", "database",    "",  NULL,  &db_requested,          "List or query database display features");
    BU_OPT(d[4], "a", "all",         "",  NULL,  &all_requested,         "List or query all feature classes");
    BU_OPT(d[5], "j", "json",        "",  NULL,  &json,                  "Write list/info results as JSON");
    BU_OPT_NULL(d[6]);
    gd->gopts = d;

    int first_pos = _view_object_first_pos(argc, argv);
    int acnt = (first_pos >= 0) ? first_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);
    gd->local_obj = local_only;

    int list_view = 0;
    int list_db = 0;
    _view_object_list_defaults(annotations_requested, db_requested, all_requested, &list_view, &list_db);

    if (help) {
	_ged_cmd_help(gedp, usage_string, d);
	return GED_HELP;
    }

    if (first_pos < 0) {
	_view_obj_list(gedp->ged_result_str, gd->cv, list_view, list_db,
	    gd->local_obj, NULL, json);
	return BRLCAD_OK;
    }

    const char *cmd = argv[first_pos];
    int sub_argc = argc - first_pos;
    const char **sub_argv = argv + first_pos;

    if (BU_STR_EQUAL(cmd, "list")) {
	if (sub_argc > 2) {
	    bu_vls_printf(gedp->ged_result_str,
		    "Usage: view feature [--all|--annotations|--database] [-L] list [glob_pattern]");
	    return BRLCAD_ERROR;
	}
	const char *glob = (sub_argc > 1) ? sub_argv[1] : NULL;
	_view_obj_list(gedp->ged_result_str, gd->cv, list_view, list_db,
	    gd->local_obj, glob, json);
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(cmd, "info")) {
	if (sub_argc < 2 || sub_argc > 3) {
	    bu_vls_printf(gedp->ged_result_str,
		    "Usage: view feature info <name> "
		    "[mode|color|visible|lcnt|type|kind|scope|overlay_class|lifecycle|"
		    "owner|owner_role|owner_generation|transient_preview|command_result]");
	    return BRLCAD_ERROR;
	}
	return _view_object_info(gd, sub_argv[1], (sub_argc == 3) ? sub_argv[2] : NULL,
		list_view, list_db, json);
    }

    if (BU_STR_EQUAL(cmd, "show") || BU_STR_EQUAL(cmd, "hide") ||
	    BU_STR_EQUAL(cmd, "delete") || BU_STR_EQUAL(cmd, "realize")) {
	if (sub_argc < 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view feature %s <name>", cmd);
	    return BRLCAD_ERROR;
	}
	gd->vobj = sub_argv[1];
	if (BU_STR_EQUAL(cmd, "show")) {
	    const char *dargv[3] = {"visible", "UP", NULL};
	    return _objs_cmd_draw(gd, 2, dargv);
	}
	if (BU_STR_EQUAL(cmd, "hide")) {
	    const char *dargv[3] = {"visible", "DOWN", NULL};
	    return _objs_cmd_draw(gd, 2, dargv);
	}
	if (BU_STR_EQUAL(cmd, "delete")) {
	    if (sub_argc != 2) {
		bu_vls_printf(gedp->ged_result_str, "Usage: view feature delete <name>");
		return BRLCAD_ERROR;
	    }
	    const char *rargv[2] = {"delete", NULL};
	    return _objs_cmd_delete(gd, 1, rargv);
	}
	std::vector<const char *> uargv;
	uargv.push_back("realize");
	for (int i = 2; i < sub_argc; i++)
	    uargv.push_back(sub_argv[i]);
	uargv.push_back(NULL);
	return _objs_cmd_update(gd, (int)uargv.size() - 1, uargv.data());
    }

    if (BU_STR_EQUAL(cmd, "style")) {
	if (sub_argc < 4) {
	    bu_vls_printf(gedp->ged_result_str,
		    "Usage: view feature style <get|set> <name> <field> [value...]");
	    return BRLCAD_ERROR;
	}
	const char *style_cmd = sub_argv[1];
	gd->vobj = sub_argv[2];
	if (BU_STR_EQUAL(style_cmd, "get")) {
	    if (sub_argc != 4) {
		bu_vls_printf(gedp->ged_result_str,
			"Usage: view feature style get <name> <field>");
		return BRLCAD_ERROR;
	    }
	    return _view_object_style_get(gd, sub_argv[3]);
	}
	if (BU_STR_EQUAL(style_cmd, "set")) {
	    return _view_object_style_set(gd, sub_argc - 3, sub_argv + 3);
	}
	bu_vls_printf(gedp->ged_result_str, "Unsupported style operation %s", style_cmd);
	return BRLCAD_ERROR;
    }

    bu_vls_printf(gedp->ged_result_str,
	    "Unsupported feature subcommand %s (valid: list, info, show, hide, delete, style, realize)",
	    cmd);
    return BRLCAD_ERROR;
}

extern "C" int
_view_cmd_annotation(void *bs, int argc, const char **argv)
{
    int help = 0;
    int local_only = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view [options] annotation [-L] <line|arrow|label|axes> <verb> <name> [args]";
    const char *purpose_string = "create and edit user annotation objects";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;
    if (!_view_object_scene_ready(gd))
	return BRLCAD_ERROR;

    argc--; argv++;

    struct bu_opt_desc d[3];
    BU_OPT(d[0], "h", "help",  "",  NULL,  &help,       "Print help");
    BU_OPT(d[1], "L", "local", "",  NULL,  &local_only, "Feature is scoped only to current/specified view");
    BU_OPT_NULL(d[2]);
    gd->gopts = d;

    int first_pos = _view_object_first_pos(argc, argv);
    int acnt = (first_pos >= 0) ? first_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);
    gd->local_obj = local_only;

    if (help) {
	_ged_cmd_help(gedp, usage_string, d);
	return GED_HELP;
    }
    if (first_pos < 0 || argc - first_pos < 3) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    const char **sub_argv = argv + first_pos;
    int sub_argc = argc - first_pos;
    const char *type = sub_argv[0];
    const char *verb = sub_argv[1];
    gd->vobj = sub_argv[2];

    std::vector<const char *> cargv;
    if (BU_STR_EQUAL(type, "arrow")) {
	if (!BU_STR_EQUAL(verb, "create")) {
	    bu_vls_printf(gedp->ged_result_str,
		    "Unsupported arrow annotation operation %s", verb);
	    return BRLCAD_ERROR;
	}
	cargv.push_back("line");
	cargv.push_back("create");
	for (int i = 3; i < sub_argc; i++)
	    cargv.push_back(sub_argv[i]);
	cargv.push_back(NULL);
	int ret = _view_cmd_lines(gd, (int)cargv.size() - 1, cargv.data());
	if (ret != BRLCAD_OK)
	    return ret;
	const char *aargv[3] = {"arrow", "1", NULL};
	return _objs_cmd_arrow(gd, 2, aargv);
    }

    cargv.push_back(type);
    cargv.push_back(verb);
    for (int i = 3; i < sub_argc; i++)
	cargv.push_back(sub_argv[i]);
    cargv.push_back(NULL);

    if (BU_STR_EQUAL(type, "line"))
	return _view_cmd_lines(gd, (int)cargv.size() - 1, cargv.data());
    if (BU_STR_EQUAL(type, "label"))
	return _view_cmd_labels(gd, (int)cargv.size() - 1, cargv.data());
    if (BU_STR_EQUAL(type, "axes"))
	return _view_cmd_axes(gd, (int)cargv.size() - 1, cargv.data());

    bu_vls_printf(gedp->ged_result_str,
	    "Unsupported annotation type %s (valid: line, arrow, label, axes)",
	    type);
    return BRLCAD_ERROR;
}

extern "C" int
_view_cmd_polygon(void *bs, int argc, const char **argv)
{
    int help = 0;
    int local_only = 0;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view [options] polygon [-L] <verb> <name> [args]";
    const char *purpose_string = "create and edit polygon annotation objects";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;
    if (!_view_object_scene_ready(gd))
	return BRLCAD_ERROR;

    argc--; argv++;

    struct bu_opt_desc d[3];
    BU_OPT(d[0], "h", "help",  "",  NULL,  &help,       "Print help");
    BU_OPT(d[1], "L", "local", "",  NULL,  &local_only, "Object is scoped only to current/specified view");
    BU_OPT_NULL(d[2]);
    gd->gopts = d;

    int first_pos = _view_object_first_pos(argc, argv);
    int acnt = (first_pos >= 0) ? first_pos : argc;
    (void)bu_opt_parse(NULL, acnt, argv, d);
    gd->local_obj = local_only;

    if (help) {
	_ged_cmd_help(gedp, usage_string, d);
	return GED_HELP;
    }
    if (first_pos < 0 || argc - first_pos < 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s\n", usage_string);
	return BRLCAD_ERROR;
    }

    const char **sub_argv = argv + first_pos;
    int sub_argc = argc - first_pos;
    const char *verb = sub_argv[0];
    gd->vobj = sub_argv[1];

    if (BU_STR_EQUAL(verb, "update")) {
	std::vector<const char *> uargv;
	uargv.push_back("update");
	for (int i = 2; i < sub_argc; i++)
	    uargv.push_back(sub_argv[i]);
	uargv.push_back(NULL);
	return _objs_cmd_update(gd, (int)uargv.size() - 1, uargv.data());
    }

    std::vector<const char *> cargv;
    cargv.push_back("polygon");
    cargv.push_back(verb);
    for (int i = 2; i < sub_argc; i++)
	cargv.push_back(sub_argv[i]);
    cargv.push_back(NULL);
    return _view_cmd_polygons(gd, (int)cargv.size() - 1, cargv.data());
}

extern "C" int
_view_cmd_db_objects(void *bs, int argc, const char **argv)
{
    int help = 0;
    struct bu_vls as_name = BU_VLS_INIT_ZERO;
    struct _ged_view_info *gd = (struct _ged_view_info *)bs;
    struct ged *gedp = gd->gedp;

    const char *usage_string = "view [options] db <add|delete|list> ...";
    const char *purpose_string = "manage Obol database display objects";
    if (_view_cmd_msgs(bs, argc, argv, usage_string, purpose_string))
	return BRLCAD_OK;
    if (!_view_object_scene_ready(gd))
	return BRLCAD_ERROR;

    argc--; argv++;

    if (argc <= 0) {
	_view_obj_list(gedp->ged_result_str, gd->cv, 0, 1, 0, NULL, 0);
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "list")) {
	const char *glob = (argc > 1) ? argv[1] : NULL;
	if (argc > 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view db list [glob_pattern]");
	    return BRLCAD_ERROR;
	}
	_view_obj_list(gedp->ged_result_str, gd->cv, 0, 1, 0, glob, 0);
	return BRLCAD_OK;
    }

    if (BU_STR_EQUAL(argv[0], "add")) {
	if (argc < 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view db add <dbpath> [--as <name>]");
	    return BRLCAD_ERROR;
	}
	struct bu_opt_desc d[3];
	BU_OPT(d[0], "h", "help",  "",      NULL,        &help,       "Print help");
	BU_OPT(d[1], "",  "as",    "name",  &bu_opt_vls, &as_name,    "View feature name");
	BU_OPT_NULL(d[2]);
	gd->gopts = d;
	const char **add_argv = argv + 1;
	int acnt = bu_opt_parse(NULL, argc - 1, add_argv, d);
	if (help) {
	    _ged_cmd_help(gedp, "view db add <dbpath> [--as <name>]", d);
	    bu_vls_free(&as_name);
	    return GED_HELP;
	}
	if (acnt < 1) {
	    bu_vls_free(&as_name);
	    bu_vls_printf(gedp->ged_result_str, "Usage: view db add <dbpath> [--as <name>]");
	    return BRLCAD_ERROR;
	}
	const char *dbpath = add_argv[0];
	const char *name = bu_vls_strlen(&as_name) ? bu_vls_cstr(&as_name) : dbpath;
	int ret = ged_view_feature_gobject_create(gedp, gd->cv, dbpath,
		name, gedp->ged_result_str) ? BRLCAD_OK : BRLCAD_ERROR;
	bu_vls_free(&as_name);
	return ret;
    }

    if (BU_STR_EQUAL(argv[0], "delete")) {
	if (argc != 2) {
	    bu_vls_printf(gedp->ged_result_str, "Usage: view db delete <name>");
	    return BRLCAD_ERROR;
	}
	gd->vobj = argv[1];
	const char *rargv[2] = {"delete", NULL};
	return _objs_cmd_delete(gd, 1, rargv);
    }

    bu_vls_printf(gedp->ged_result_str,
	    "Unsupported db subcommand %s (valid: add, delete, list)", argv[0]);
    return BRLCAD_ERROR;
}


// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
