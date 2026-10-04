/*                    R E P A I R _ B O U N D A R Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * Distributed under the terms of the GNU Lesser General Public License
 * (LGPL), version 2.1.
 */
/** @file cdt/repair_boundary.h
 *
 * Internal boundary constraints and bounded local mesh reconstruction.
 */

#ifndef LIBBREP_CDT_REPAIR_BOUNDARY_H
#define LIBBREP_CDT_REPAIR_BOUNDARY_H

#include "common.h"

#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "brep/cdt.h"

struct ON_Brep_CDT_State;

typedef std::array<double, 3> repair_point_key;
typedef std::array<repair_point_key, 3> repair_triangle_key;
typedef std::array<repair_point_key, 2> repair_triangle_edge_key;

struct repair_fast_trim_sample {
    fastf_t parameter;
    fastf_t uv[2];
    fastf_t point[3];
    ON_3dPoint *source_point = NULL;
};

struct repair_fast_constraint_store {
    int face_index = -1;
    std::map<int, std::vector<repair_fast_trim_sample>> trims;
    size_t constrained_edges = 0;
    size_t constrained_samples = 0;
};

struct repair_boundary_patch {
    std::vector<int> faces;
    std::vector<fastf_t> vertices;
    std::vector<ON_3dPoint *> source_points;
    std::vector<double> direct_surface_deviations;
    size_t constrained_edges = 0;
    size_t constrained_samples = 0;
    size_t periodic_u_samples = 0;
    size_t periodic_v_samples = 0;
    size_t strip_interior_rows = 0;
    int strip_correspondence_step = 0;
    size_t strip_correspondence_shift = 0;
    bool analytic_torus_strip = false;
    bool analytic_cylinder_strip = false;
};

struct repair_degenerate_neighborhood_stats {
    int components = 0;
    int removed_faces = 0;
    int added_faces = 0;
    int boundary_edges = 0;
    double max_center_offset = 0.0;
    bool inferred_topology = false;
    std::set<int> source_faces;
};

/* Patches retain source point identities until the repair driver assembles
 * them.  Construction alone does not certify a patch for export. */
repair_triangle_edge_key
repair_triangle_edge(const repair_point_key &first,
	const repair_point_key &second);

bool
repair_hybrid_fallback_preflight(size_t geometric_degenerate_faces,
	size_t invalid_vertex_links, size_t face_count,
	size_t unmatched_edges, size_t local_repair_limit);

size_t
repair_adaptive_hole_edge_budget(
	const struct brep_cdt_repair_settings *settings,
	const struct bg_trimesh_repair_report *mesh_report);

fastf_t
repair_adaptive_hole_area_budget(
	const struct brep_cdt_repair_settings *settings,
	const struct bg_trimesh_repair_settings *mesh_settings,
	const struct bg_trimesh_repair_report *mesh_report);

bool
repair_closed_periodic_surface(struct ON_Brep_CDT_State *s_cdt,
	int face_index, size_t max_points, fastf_t allowed_deviation,
	int64_t deadline, repair_boundary_patch &patch);

bool
repair_closed_spherical_surface(struct ON_Brep_CDT_State *s_cdt,
	int face_index, size_t max_points, fastf_t allowed_deviation,
	int64_t deadline, repair_boundary_patch &patch);

bool
repair_constrained_periodic_strip(struct ON_Brep_CDT_State *s_cdt,
	int face_index, const repair_fast_constraint_store &constraints,
	size_t max_points, fastf_t allowed_deviation, int64_t deadline,
	repair_boundary_patch &patch, int required_step = 0,
	size_t required_shift = std::numeric_limits<size_t>::max());

bool
repair_constrained_spherical_cap(struct ON_Brep_CDT_State *s_cdt,
	int face_index, const repair_fast_constraint_store &constraints,
	size_t max_points, fastf_t allowed_deviation, int64_t deadline,
	repair_boundary_patch &patch);

bool
repair_constrained_topological_disk(struct ON_Brep_CDT_State *s_cdt,
	int face_index, const repair_fast_constraint_store &constraints,
	size_t max_points, repair_boundary_patch &patch);

bool
repair_degenerate_neighborhoods(int **faces, int *face_count,
	fastf_t **vertices, int *vertex_count, int rigorous_face_count,
	std::vector<int> &source_faces, double allowed_deviation,
	size_t max_boundary_edges, repair_degenerate_neighborhood_stats *stats);

bool
repair_failed_face_from_rigorous_boundary(
	struct ON_Brep_CDT_State *s_cdt, int **faces, int *face_count,
	fastf_t **vertices, int *vertex_count, int rigorous_face_count,
	std::vector<int> &source_faces, double allowed_deviation,
	size_t max_boundary_edges, repair_degenerate_neighborhood_stats *stats,
	int requested_source_face = -1,
	const std::set<std::pair<int, int>> *requested_boundary = NULL,
	bool require_solid = true, std::vector<int> *vertex_remap = NULL);

repair_fast_constraint_store
repair_fast_face_constraints(struct ON_Brep_CDT_State *s_cdt,
	int face_index);

bool
repair_internal_face_seam(const ON_BrepTrim *trim,
	const ON_BrepFace &face);

bool
repair_pair_near_boundary_cracks(int *faces, int face_count,
	const fastf_t *vertices, int vertex_count, int rigorous_face_count,
	double allowed_deviation, size_t *welded_vertices);

bool
repair_patch_matches_source_boundary(const int *faces, int face_count,
	const std::vector<ON_3dPoint *> &sources,
	const repair_boundary_patch &reference);

bool
repair_topological_disk_stage(int stage);

size_t
repair_weld_roundoff_vertices(int *faces, int face_count,
	const fastf_t *vertices, int vertex_count);

#endif /* LIBBREP_CDT_REPAIR_BOUNDARY_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
