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
struct ged_diagnostic_hud_label;
struct ged_view_context;
struct ged_view_feature_style;

std::string ged_obol_view_scope_name(struct ged_view_context *view_ctx);

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
