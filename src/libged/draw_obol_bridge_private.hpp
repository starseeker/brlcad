/*        D R A W _ O B O L _ B R I D G E _ P R I V A T E . H P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file draw_obol_bridge_private.hpp
 *
 * Hidden value-conversion helpers shared by the single-responsibility Obol
 * bridge units.  None of these functions owns, attaches, or preserves a
 * controller, scene, or reference token.
 */

#ifndef LIBGED_DRAW_OBOL_BRIDGE_PRIVATE_HPP
#define LIBGED_DRAW_OBOL_BRIDGE_PRIVATE_HPP

#include "BObol/BViewStore.h"
#include "vmath.h"

#include <cstdint>
#include <string>
#include <vector>

class BObolViewController;
class BObolLodService;
class BObolSceneController;
class SoNode;
class SoBRLDatabaseSource;
class SoBRLSceneGroup;
class SoBRLVListShape;
struct BObolProgressiveOptions;
struct BObolProgressiveStatus;
struct ged;
struct ged_bobol_publication_context;
struct ged_draw_obol_database_source_record;
struct ged_drawable;
struct ged_diagnostic_hud_label;
struct ged_view_context;
struct ged_view_feature_style;

std::string ged_obol_view_scope_name(struct ged_view_context *view_ctx);
struct ged_drawable *ged_obol_gdp(struct ged *gedp);
void ged_obol_state_destroy(struct ged_drawable *gdp);
uint32_t ged_obol_fold_revision(uint64_t revision);
int ged_obol_view_scope_is_independent(struct ged_view_context *view_ctx);
int ged_obol_observer_ensure(struct ged *gedp, struct ged_drawable *gdp);
void ged_obol_configure_compact_realization(struct ged *gedp,
    struct ged_view_context *view_ctx,
    BObolSceneController *scene_controller);
int ged_obol_register_progressive_provider(struct ged *gedp,
    struct ged_view_context *view_ctx, BObolViewController *controller);
int ged_obol_progressive_advance_provider(
    BObolViewController *controller, void *user_data,
    const BObolProgressiveOptions *options,
    BObolProgressiveStatus *status);
BObolLodService *ged_obol_lod_service_ensure(
    BObolViewController *controller);
int ged_obol_bind_view_render_root(struct ged_view_context *view_ctx,
    BObolSceneController *shared_scene,
    BObolViewController *view_controller);
int ged_obol_node_contains(SoNode *root, SoNode *target);
int ged_obol_lod_draw_mode_from_ged(int mode);
int ged_obol_lod_draw_mode_to_ged(int mode);
int ged_obol_lod_draw_mode_from_database_source(
    const SoBRLDatabaseSource *source);
int ged_obol_database_draw_mode_from_ged(int mode);
int ged_obol_database_draw_mode_to_ged(int mode);
int ged_obol_database_representation_mode_from_ged(int mode);
const char *ged_obol_skip_leading_slash(const char *path);
int ged_obol_path_equal(const char *a, const char *b);
int ged_obol_path_has_prefix(const char *path, const char *prefix);
int ged_obol_path_has_component_name(
    const char *path, const char *name, size_t first_idx);
std::string ged_obol_semantic_path_string(const char *path);
int ged_obol_path_has_semantic_prefix(
    const char *path, const char *prefix);
void ged_obol_append_unique_path(
    std::vector<std::string> &paths, const char *path);
std::string ged_obol_group_intent_path(const char *group_path);
std::string ged_obol_group_path_from_record_path(const char *path);
bool ged_obol_group_parent_leaf(
    const std::string &path, std::string &parent, std::string &leaf);
int ged_obol_intent_is_ged_draw_group(const SbString &intent);
const char *ged_obol_group_record_path(
    const SoBRLSceneGroup *scene_group);
int ged_obol_scene_clear_controller(BObolSceneController *scene);
void ged_obol_frontier_visibility_changes_apply(struct ged *gedp);
void ged_obol_frontier_visibility_snapshot_apply(struct ged *gedp);
const char *ged_obol_shape_node_record_role(SoNode *node);
int ged_obol_database_source_summary_ged_mode(
    const BObolDatabaseSourceSummary &summary);
int ged_obol_database_source_summary_matches_mode(
    const BObolDatabaseSourceSummary &summary, int draw_mode);
int ged_obol_database_source_instance_in_scope(
    const BObolDatabaseSourceSummary &summary,
    struct ged_view_context *view_ctx);
std::string ged_obol_database_source_owner_group_path_from_summary(
    const BObolDatabaseSourceSummary &summary);
void ged_obol_append_database_source_instance_key(
    std::vector<std::string> &instance_keys,
    const BObolDatabaseSourceSummary &summary);
std::vector<std::string> ged_obol_database_source_instance_keys_for_path(
    struct ged *gedp, const char *path, int draw_mode,
    int allow_path_prefix,
    const struct ged_bobol_publication_context *publication = NULL);
SoBRLDatabaseSource *ged_obol_database_source_for_instance_key(
    BObolSceneController *scene, const std::string &source_instance_key);
int ged_obol_remove_paths(
    const std::vector<std::string> &paths,
    struct ged_view_context *view_ctx,
    BObolSceneController *scene, int draw_mode = -1);
int ged_obol_remove_instance_keys(
    const std::vector<std::string> &instance_keys,
    BObolSceneController *scene);
int ged_obol_clear_database_sources_in_scope(
    BObolSceneController *scene, struct ged_view_context *view_ctx);
int ged_obol_prune_empty_groups(BObolSceneController *scene);
SbColor ged_obol_summary_material_color(
    const BObolDatabaseSourceSummary &summary);
SoBRLDatabaseSource *ged_obol_owned_database_source_for_path(
    struct ged *gedp, const char *path);
bool ged_obol_database_source_controller_summary_for_path(
    BObolSceneController *scene, const char *path,
    BObolDatabaseSourceSummary &summary);
bool ged_obol_database_source_controller_summary_for_source(
    BObolSceneController *scene, SoBRLDatabaseSource *source,
    BObolDatabaseSourceSummary &summary);
int ged_obol_database_source_record_from_summary(
    struct ged_draw_obol_database_source_record *out,
    const BObolDatabaseSourceSummary &summary);
int ged_obol_database_source_exact_draw_mode_to_ged(
    struct ged *gedp, const BObolDatabaseSourceSummary &summary,
    SoBRLDatabaseSource *source);
SoBRLVListShape *ged_obol_owned_annotation_vlist_shape_for_source(
    SoBRLDatabaseSource *source, const char *fallback_path);
int32_t ged_obol_vlist_command_from_ged(int command, size_t index);
int ged_obol_vlist_shape_is_annotation(SoBRLVListShape *shape);
int ged_obol_vlist_shape_has_annotation_record(SoBRLVListShape *shape);
int ged_obol_database_source_scene_instance_for_path(
    struct ged *gedp, const char *path,
    BObolSceneController **scene_out,
    std::string &source_instance_key,
    const struct ged_bobol_publication_context *publication = NULL);

BObolViewController *ged_obol_view_controller_for_scope(
    struct ged_view_context *view_ctx, int local, int sync_current_scene);

BObolFeatureStyle ged_obol_feature_style_from_ged(
    const struct ged_view_feature_style *style);
int32_t ged_obol_line_command_from_ged(int command, size_t index);
std::vector<SbVec3f> ged_obol_points_from_ged(
    const point_t *points, size_t point_count);
std::vector<int32_t> ged_obol_commands_from_ged(
    const int *commands, size_t point_count);
std::vector<int32_t> ged_obol_indices_from_ged(
    const int *indices, size_t index_count);
std::vector<SbVec3f> ged_obol_vectors_from_ged(
    const vect_t *vectors, size_t vector_count);

BObolOverlayInfo ged_obol_model_overlay_info(
    struct ged_view_context *view_ctx,
    BObolOverlayClass overlay_class,
    BObolOverlayLifecycle lifecycle,
    BObolOverlayOrder order,
    const char *source_path);
int ged_obol_feature_mark_overlay(
    BObolViewController *controller,
    BObolFeatureHandle handle,
    const BObolOverlayInfo &overlay);
BObolLabel ged_obol_label_from_hud(
    const struct ged_diagnostic_hud_label &label);

#endif /* LIBGED_DRAW_OBOL_BRIDGE_PRIVATE_HPP */
