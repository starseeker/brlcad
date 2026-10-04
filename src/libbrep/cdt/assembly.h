/*                     A S S E M B L Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * Distributed under the terms of the GNU Lesser General Public License
 * (LGPL), version 2.1.
 */
/** @file cdt/assembly.h
 *
 * Shared validation and transactional cleanup of assembled CDT meshes.
 */

#ifndef LIBBREP_CDT_ASSEMBLY_H
#define LIBBREP_CDT_ASSEMBLY_H

#include "common.h"
#include "brep/defines.h"

struct assembled_mesh_validation {
    size_t invalid_indices = 0;
    size_t nonfinite_vertices = 0;
    size_t unused_vertices = 0;
    size_t degenerate_faces = 0;
    size_t invalid_vertex_links = 0;
    size_t intersecting_triangle_pairs = 0;
    int first_intersection[2] = {-1, -1};
    bool first_intersection_point_valid = false;
    double first_intersection_point[3] = {0.0, 0.0, 0.0};
};

/* These checks also apply to repair candidates before certification. */
bool
assembled_mesh_triangle_degenerate(const ON_3dPoint points[3]);

bool
assembled_mesh_validate(int vertex_count, int face_count,
	const fastf_t *vertices, const int *faces,
	assembled_mesh_validation *validation, bool check_intersections = true);

int
closed_mesh_orientation_sync(int *faces, int face_count,
	fastf_t *vertices, int vertex_count, int *face_normals);

int
closed_mesh_component_filter(int *faces, int *face_count,
	fastf_t *vertices, int vertex_count, int *face_normals,
	int *source_faces = NULL, double overlap_tolerance = 0.0,
	bool source_is_solid = true, bool emit_diagnostics = true,
	size_t *source_triangles = NULL);

#endif /* LIBBREP_CDT_ASSEMBLY_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
