/*                  A S S E M B L Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * Distributed under the terms of the GNU Lesser General Public License
 * (LGPL), version 2.1.
 */
/** @file cdt/assembly.cpp
 *
 * Shared validation and transactional cleanup of assembled CDT meshes.
 */

#include "common.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include "bg/tri_pt.h"
#include "bg/tri_tri.h"
#include "cdt.h"
#include "assembly.h"

static bool
assembled_mesh_collect_candidate(size_t triangle, void *context)
{
    std::vector<size_t> *candidates =
	(std::vector<size_t> *)context;
    candidates->push_back(triangle);
    return true;
}

bool
assembled_mesh_triangle_degenerate(const ON_3dPoint points[3])
{
    const double ROUNDING_ERROR_FACTOR = 64.0;
    const ON_3dVector ab = points[1] - points[0];
    const ON_3dVector ac = points[2] - points[0];
    const ON_3dVector bc = points[2] - points[1];
    const double longest_sq = std::max(ab.LengthSquared(),
	std::max(ac.LengthSquared(), bc.LengthSquared()));
    const double doubled_area = ON_CrossProduct(ab, ac).Length();
    return !(longest_sq > 0.0) || !std::isfinite(doubled_area) ||
	doubled_area <= ROUNDING_ERROR_FACTOR * std::numeric_limits<double>::epsilon() *
	longest_sq;
}

bool
assembled_mesh_validate(int vertex_count, int face_count,
	const fastf_t *vertices, const int *faces,
	assembled_mesh_validation *validation, bool check_intersections)
{
    if (!validation || vertex_count <= 0 || face_count <= 0 ||
	    !vertices || !faces)
	return false;
    *validation = assembled_mesh_validation();

    typedef std::pair<int, int> mesh_edge;
    std::vector<std::vector<mesh_edge>> vertex_link_edges(
	(size_t)vertex_count);
    std::vector<bool> used_vertices((size_t)vertex_count, false);
    for (int vertex = 0; vertex < vertex_count; ++vertex) {
	for (int coordinate = 0; coordinate < 3; ++coordinate) {
	    if (!std::isfinite(vertices[3 * vertex + coordinate])) {
		validation->nonfinite_vertices++;
		break;
	    }
	}
    }

    for (int face = 0; face < face_count; ++face) {
	ON_3dPoint points[3];
	for (int corner = 0; corner < 3; ++corner) {
	    const int vertex = faces[3 * face + corner];
	    if (vertex < 0 || vertex >= vertex_count) {
		validation->invalid_indices++;
		return false;
	    }
	    used_vertices[(size_t)vertex] = true;
	    points[corner] = ON_3dPoint(vertices[3 * vertex],
		vertices[3 * vertex + 1], vertices[3 * vertex + 2]);
	}
	for (int corner = 0; corner < 3; ++corner) {
	    const int vertex = faces[3 * face + corner];
	    vertex_link_edges[(size_t)vertex].push_back(mesh_edge(
		faces[3 * face + (corner + 1) % 3],
		faces[3 * face + (corner + 2) % 3]));
	}
	if (assembled_mesh_triangle_degenerate(points))
	    validation->degenerate_faces++;
    }

    for (bool used : used_vertices) {
	if (!used)
	    validation->unused_vertices++;
    }

    for (int vertex = 0; vertex < vertex_count; ++vertex) {
	const std::vector<mesh_edge> &links =
	    vertex_link_edges[(size_t)vertex];
	if (links.empty())
	    continue;
	std::map<int, std::vector<int>> adjacency;
	for (const mesh_edge &link : links) {
	    adjacency[link.first].push_back(link.second);
	    adjacency[link.second].push_back(link.first);
	}
	std::set<int> reached;
	std::queue<int> work;
	work.push(adjacency.begin()->first);
	reached.insert(adjacency.begin()->first);
	bool valid_link = true;
	while (!work.empty()) {
	    const int current = work.front();
	    work.pop();
	    const std::vector<int> &neighbors = adjacency[current];
	    valid_link = valid_link && neighbors.size() == 2;
	    for (int neighbor : neighbors) {
		if (reached.insert(neighbor).second)
		    work.push(neighbor);
	    }
	}
	if (!valid_link || reached.size() != adjacency.size())
	    validation->invalid_vertex_links++;
    }

    /* Edge incidence and vertex links do not detect two otherwise valid
     * closed sheets crossing away from shared topology.  Query an R-tree of
     * prior triangle bounds and apply the validated library predicate only to
     * nonadjacent candidates.  One certified counterexample is sufficient
     * to fail closed and bounds work on malformed, heavily overlapping
     * meshes. */
    if (check_intersections && !validation->nonfinite_vertices &&
	    !validation->degenerate_faces) {
	RTree<size_t, double, 3> triangle_index;
	for (int face = 0; face < face_count; ++face) {
	    double minimum[3] = {
		std::numeric_limits<double>::infinity(),
		std::numeric_limits<double>::infinity(),
		std::numeric_limits<double>::infinity()
	    };
	    double maximum[3] = {
		-std::numeric_limits<double>::infinity(),
		-std::numeric_limits<double>::infinity(),
		-std::numeric_limits<double>::infinity()
	    };
	    point_t first_points[3];
	    for (int corner = 0; corner < 3; ++corner) {
		const int vertex = faces[3 * face + corner];
		VSET(first_points[corner], vertices[3 * vertex],
		    vertices[3 * vertex + 1], vertices[3 * vertex + 2]);
		for (int axis = 0; axis < 3; ++axis) {
		    minimum[axis] = std::min(minimum[axis],
			(double)first_points[corner][axis]);
		    maximum[axis] = std::max(maximum[axis],
			(double)first_points[corner][axis]);
		}
	    }
	    std::vector<size_t> candidates;
	    triangle_index.Search(minimum, maximum,
		assembled_mesh_collect_candidate, &candidates);
	    std::sort(candidates.begin(), candidates.end());
	    for (size_t candidate : candidates) {
		const int *first = &faces[3 * face];
		const int *second = &faces[3 * candidate];
		bool adjacent = false;
		for (int i = 0; i < 3 && !adjacent; ++i) {
		    for (int j = 0; j < 3; ++j) {
			if (first[i] == second[j]) {
			    adjacent = true;
			    break;
			}
		    }
		}
		if (adjacent)
		    continue;
		point_t second_points[3];
		for (int corner = 0; corner < 3; ++corner) {
		    const int vertex = second[corner];
		    VSET(second_points[corner], vertices[3 * vertex],
			vertices[3 * vertex + 1],
			vertices[3 * vertex + 2]);
		}
		const bool intersects = cdt_tri_tri_intersection(first_points,
		    second_points);
		if (!intersects)
		    continue;
		validation->intersecting_triangle_pairs = 1;
		validation->first_intersection[0] = (int)candidate;
		validation->first_intersection[1] = face;
		int coplanar = 0;
		point_t start = VINIT_ZERO;
		point_t end = VINIT_ZERO;
		if (bg_tri_tri_isect_with_line(first_points[0],
			first_points[1], first_points[2], second_points[0],
			second_points[1], second_points[2], &coplanar,
			&start, &end) && !coplanar) {
		    for (int axis = 0; axis < 3; ++axis)
			validation->first_intersection_point[axis] =
			    0.5 * (start[axis] + end[axis]);
		    validation->first_intersection_point_valid = true;
		}
		break;
	    }
	    if (validation->intersecting_triangle_pairs)
		break;
	    triangle_index.Insert(minimum, maximum, (size_t)face);
	}
    }

    return !validation->invalid_indices && !validation->nonfinite_vertices &&
	!validation->unused_vertices && !validation->degenerate_faces &&
	!validation->invalid_vertex_links &&
	!validation->intersecting_triangle_pairs;
}

int
cdt_test_assembled_mesh_validation(void)
{
    fastf_t vertices[] = {
	0.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0,
	-1.0, 0.0, 0.0,
	0.0, -1.0, 0.0,
	0.0, 0.0, -1.0
    };
    int tetrahedron[] = {
	0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3
    };
    assembled_mesh_validation validation;
    if (!assembled_mesh_validate(4, 4, vertices, tetrahedron,
	    &validation))
	return 1;

    /* Two otherwise closed tetrahedra meeting only at vertex zero have two
     * disjoint cycles in that vertex's link.  Edge-incidence checks alone do
     * not detect this bow-tie topology. */
    int bow_tie[] = {
	0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3,
	0, 5, 4, 0, 4, 6, 0, 6, 5, 4, 5, 6
    };
    if (assembled_mesh_validate(7, 8, vertices, bow_tie, &validation) ||
	    validation.invalid_vertex_links != 1)
	return 2;

    fastf_t vertices_with_unused[24];
    std::copy(vertices, vertices + 21, vertices_with_unused);
    vertices_with_unused[21] = 2.0;
    vertices_with_unused[22] = 2.0;
    vertices_with_unused[23] = 2.0;
    if (assembled_mesh_validate(8, 4, vertices_with_unused, tetrahedron,
	    &validation) || validation.unused_vertices != 4)
	return 3;

    int degenerate[] = {
	0, 2, 1, 0, 1, 3, 0, 3, 2, 0, 0, 1
    };
    if (assembled_mesh_validate(4, 4, vertices, degenerate, &validation) ||
	    validation.degenerate_faces != 1)
	return 4;

    fastf_t crossing_vertices[] = {
	0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
	0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
	0.25, 0.25, 0.25, 1.25, 0.25, 0.25,
	0.25, 1.25, 0.25, 0.25, 0.25, 1.25
    };
    int crossing_tetrahedra[] = {
	0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3,
	4, 6, 5, 4, 5, 7, 4, 7, 6, 5, 6, 7
    };
    if (assembled_mesh_validate(8, 8, crossing_vertices,
	    crossing_tetrahedra, &validation) ||
	    validation.intersecting_triangle_pairs != 1 ||
	    validation.first_intersection[0] < 0 ||
	    validation.first_intersection[1] < 0 ||
	    !validation.first_intersection_point_valid)
	return 5;

    /* Adjacent face chords can have touching bounding boxes and separated,
     * collinear boundary intervals.  The base predicate's coarse epsilon
     * reports this NIST-derived configuration as an intersection, but its
     * reported endpoint does not lie on both triangles. */
    point_t separated_first[3] = {
	{514.97292083479124, -327.29326229964499, 77.978026208414306},
	{515.54826425965757, -326.98097089301785, 78.005348166208819},
	{515.6018756882213, -326.66867948639072, 78.032670124003332}
    };
    point_t separated_second[3] = {
	{515.45117516473545, -327.29326229964499, 77.978026208414306},
	{515.45520000264924, -327.36064849937765, 77.971900884443286},
	{515.46589947082487, -327.3585443137215, 77.97209910798874}
    };
    if (cdt_tri_tri_intersection(separated_first, separated_second))
	return 6;

    return 0;
}

struct compact_mesh_edge_use {
    int first;
    int second;
    int face;
};

static int
compact_mesh_component_root(std::vector<int> &parents, int item)
{
    int root = item;
    while (parents[(size_t)root] >= 0)
	root = parents[(size_t)root];
    while (item != root) {
	const int next = parents[(size_t)item];
	parents[(size_t)item] = root;
	item = next;
    }
    return root;
}

/* Find face components through exactly two-use edges.  Sorting a packed edge
 * array is substantially smaller than storing millions of edge keys and tiny
 * neighbor vectors in separately allocated tree nodes. */
static void
compact_mesh_face_components(const int *faces, int face_count,
	std::vector<int> &component_id,
	std::vector<std::vector<int>> &components)
{
    component_id.assign((size_t)face_count, -1);
    components.clear();
    if (!faces || face_count <= 0)
	return;

    std::vector<compact_mesh_edge_use> edge_uses;
    edge_uses.reserve((size_t)face_count * 3);
    for (int face = 0; face < face_count; ++face) {
	for (int corner = 0; corner < 3; ++corner) {
	    int first = faces[(size_t)face * 3 + corner];
	    int second = faces[(size_t)face * 3 + (corner + 1) % 3];
	    if (first > second)
		std::swap(first, second);
	    compact_mesh_edge_use use = {first, second, face};
	    edge_uses.push_back(use);
	}
    }
    std::sort(edge_uses.begin(), edge_uses.end(),
	[](const compact_mesh_edge_use &a,
	    const compact_mesh_edge_use &b) {
	    if (a.first != b.first)
		return a.first < b.first;
	    if (a.second != b.second)
		return a.second < b.second;
	    return a.face < b.face;
	});

    std::vector<int> parents((size_t)face_count, -1);
    for (size_t begin = 0; begin < edge_uses.size();) {
	size_t end = begin + 1;
	while (end < edge_uses.size() &&
		edge_uses[end].first == edge_uses[begin].first &&
		edge_uses[end].second == edge_uses[begin].second)
	    ++end;
	if (end - begin == 2) {
	    int first_root = compact_mesh_component_root(parents,
		edge_uses[begin].face);
	    int second_root = compact_mesh_component_root(parents,
		edge_uses[begin + 1].face);
	    if (first_root != second_root) {
		if (parents[(size_t)first_root] >
			parents[(size_t)second_root])
		    std::swap(first_root, second_root);
		parents[(size_t)first_root] += parents[(size_t)second_root];
		parents[(size_t)second_root] = first_root;
	    }
	}
	begin = end;
    }

    std::vector<int> root_component((size_t)face_count, -1);
    for (int face = 0; face < face_count; ++face) {
	const int root = compact_mesh_component_root(parents, face);
	if (root_component[(size_t)root] < 0) {
	    root_component[(size_t)root] = (int)components.size();
	    components.push_back(std::vector<int>());
	}
	const int component = root_component[(size_t)root];
	component_id[(size_t)face] = component;
	components[(size_t)component].push_back(face);
    }
}

/* Correct triangle winding only after the complete mesh proves closed and
 * manifold and misorientation is its sole defect.  bg_trimesh_sync establishes
 * a consistent orientation for each connected component, but its seed may
 * choose either of the two possible component orientations.  Preserve the
 * OpenNURBS face orientation by selecting the result that changes the fewest
 * source triangles in each component.  A tied component has no authoritative
 * majority, so leave the mesh unchanged rather than guessing.
 *
 * This function is transactional: faces and face normals are not modified
 * unless the candidate subsequently passes the complete solid test. */
int
closed_mesh_orientation_sync(int *faces, int face_count,
	fastf_t *vertices, int vertex_count, int *face_normals)
{
    const bool debug_repair_topology =
	getenv("BRLCAD_CDT_DEBUG_REPAIR_TOPOLOGY") != NULL;
    if (!faces || !vertices || face_count < 4 || vertex_count < 4)
	return -1;

    struct bg_trimesh_solid_errors errors = BG_TRIMESH_SOLID_ERRORS_INIT_NULL;
    const int initial_invalid = bg_trimesh_solid2(vertex_count, face_count,
	vertices, faces, &errors);
    if (!initial_invalid) {
	bg_free_trimesh_solid_errors(&errors);
	return 0;
    }

    const bool only_misoriented = errors.misoriented.count > 0 &&
	errors.degenerate.count == 0 && errors.unmatched.count == 0 &&
	errors.excess.count == 0;
    if (debug_repair_topology)
	bu_log("Orientation sync input: degenerate %d, unmatched %d, "
	    "excess %d, misoriented %d\n", errors.degenerate.count,
	    errors.unmatched.count, errors.excess.count,
	    errors.misoriented.count);
    bg_free_trimesh_solid_errors(&errors);
    if (!only_misoriented)
	return -1;

    std::vector<int> original(faces, faces + 3 * face_count);
    std::vector<int> candidate(original);

    if (bg_trimesh_sync(candidate.data(), candidate.data(), face_count) <= 0) {
	if (debug_repair_topology)
	    bu_log("Orientation sync could not establish face adjacency\n");
	return -1;
    }

    /* Determine the connected closed-shell components from shared edges. */
    std::vector<int> component_id;
    std::vector<std::vector<int>> components;
    compact_mesh_face_components(original.data(), face_count,
	component_id, components);

    for (const std::vector<int> &component : components) {
	size_t changed = 0;
	for (int face : component) {
	    const size_t offset = (size_t)face * 3;
	    if (candidate[offset] != original[offset] ||
		    candidate[offset + 1] != original[offset + 1] ||
		    candidate[offset + 2] != original[offset + 2])
		changed++;
	}
	if (changed * 2 == component.size()) {
	    if (debug_repair_topology)
		bu_log("Orientation sync component tied: %zu/%zu changed\n",
		    changed, component.size());
	    return -1;
	}
	if (changed * 2 > component.size()) {
	    for (int face : component)
		std::swap(candidate[(size_t)face * 3],
		    candidate[(size_t)face * 3 + 1]);
	}
    }

    struct bg_trimesh_solid_errors candidate_errors =
	BG_TRIMESH_SOLID_ERRORS_INIT_NULL;
    const int candidate_invalid = bg_trimesh_solid2(vertex_count,
	face_count, vertices, candidate.data(), &candidate_errors);
    if (candidate_invalid != 0) {
	if (debug_repair_topology)
	    bu_log("Orientation sync candidate remained non-solid: "
		"degenerate %d, unmatched %d, excess %d, misoriented %d\n",
		candidate_errors.degenerate.count,
		candidate_errors.unmatched.count,
		candidate_errors.excess.count,
		candidate_errors.misoriented.count);
	bg_free_trimesh_solid_errors(&candidate_errors);
	return -1;
    }
    bg_free_trimesh_solid_errors(&candidate_errors);

    int changed_count = 0;
    for (int face = 0; face < face_count; ++face) {
	const size_t offset = (size_t)face * 3;
	if (candidate[offset] == original[offset] &&
		candidate[offset + 1] == original[offset + 1] &&
		candidate[offset + 2] == original[offset + 2])
	    continue;
	changed_count++;
	if (face_normals)
	    std::swap(face_normals[offset], face_normals[offset + 1]);
    }

    std::copy(candidate.begin(), candidate.end(), faces);
    return changed_count;
}

/* A failed local face remesh can leave a few open triangle islands in addition
 * to the useful tessellation of the source solid.  Those islands are removed
 * only under one of two proofs: every retained component is independently
 * closed and solid, or (for a partial display mesh) every point of an island
 * is already covered within tolerance by retained triangles from the same
 * source BREP face.  Valid disjoint shells are retained, and no change is made
 * if neither proof succeeds.
 *
 * The operation is transactional with respect to the caller's face arrays.
 * Vertex arrays may retain now-unused entries; that is legal for a BoT and
 * keeps the face-normal and point ownership maps stable. */
int
closed_mesh_component_filter(int *faces, int *face_count,
	fastf_t *vertices, int vertex_count, int *face_normals,
	int *source_faces, double overlap_tolerance,
	bool source_is_solid, bool emit_diagnostics,
	size_t *source_triangles)
{
    if (!faces || !face_count || !vertices || *face_count < 2 ||
	    vertex_count < 4)
	return 0;

    /* Validate all indices before making the first transactional edit. */
    for (int face = 0; face < *face_count; ++face) {
	for (int corner = 0; corner < 3; ++corner) {
	    const int vertex = faces[(size_t)face * 3 + corner];
	    if (vertex < 0 || vertex >= vertex_count)
		return 0;
	}
    }

    /* Distinct face-mesh vertices may collapse to the same shared output
     * point.  Cull the resulting zero-area index triangles before component
     * discovery; otherwise their repeated edge can falsely bridge two
     * components in this audit even though rt_bot_split correctly ignores
     * that bridge. */
    const int original_face_count = *face_count;
    std::vector<int> compact_faces;
    std::vector<int> compact_normals;
    std::vector<int> compact_sources;
    std::vector<size_t> compact_triangles;
    compact_faces.reserve((size_t)original_face_count * 3);
    if (face_normals)
	compact_normals.reserve((size_t)original_face_count * 3);
    if (source_faces)
	compact_sources.reserve((size_t)original_face_count);
    if (source_triangles)
	compact_triangles.reserve((size_t)original_face_count);
    for (int face = 0; face < original_face_count; ++face) {
	const size_t offset = (size_t)face * 3;
	if (faces[offset] == faces[offset + 1] ||
		faces[offset + 1] == faces[offset + 2] ||
		faces[offset + 2] == faces[offset])
	    continue;
	for (int corner = 0; corner < 3; ++corner) {
	    compact_faces.push_back(faces[offset + corner]);
	    if (face_normals)
		compact_normals.push_back(face_normals[offset + corner]);
	}
	if (source_faces)
	    compact_sources.push_back(source_faces[face]);
	if (source_triangles)
	    compact_triangles.push_back(source_triangles[face]);
    }
    int removed_count = original_face_count -
	(int)(compact_faces.size() / 3);
    if (removed_count) {
	std::copy(compact_faces.begin(), compact_faces.end(), faces);
	if (face_normals)
	    std::copy(compact_normals.begin(), compact_normals.end(),
		face_normals);
	*face_count = (int)(compact_faces.size() / 3);
	if (source_faces)
	    std::copy(compact_sources.begin(), compact_sources.end(),
		source_faces);
	if (source_triangles)
	    std::copy(compact_triangles.begin(), compact_triangles.end(),
		source_triangles);
    }
    if (*face_count < 2)
	return removed_count;

    std::vector<int> component_id;
    std::vector<std::vector<int>> components;
    compact_mesh_face_components(faces, *face_count, component_id,
	components);
    if (components.size() < 2)
	return removed_count;

    std::vector<bool> retain(components.size(), false);
    size_t retained_face_count = 0;
    for (size_t component = 0;
	    source_is_solid && component < components.size(); ++component) {
	std::vector<int> component_faces;
	component_faces.reserve(components[component].size() * 3);
	for (int face : components[component]) {
	    const size_t offset = (size_t)face * 3;
	    component_faces.push_back(faces[offset]);
	    component_faces.push_back(faces[offset + 1]);
	    component_faces.push_back(faces[offset + 2]);
	}
	if (component_faces.size() >= 12 &&
		bg_trimesh_solid2(vertex_count,
		    (int)components[component].size(), vertices,
		    component_faces.data(), NULL) == 0) {
	    retain[component] = true;
	    retained_face_count += components[component].size();
	}
    }
    if (retained_face_count == (size_t)*face_count)
	return removed_count;

    /* Closed sub-shells are not by themselves proof that every open
     * component is redundant.  In particular, a valid multi-face solid can
     * contain a small independently closed feature alongside a much larger
     * component whose shared-edge assembly is defective.  Preserve the full
     * candidate unless the same-face overlap proof below can account for
     * every discarded triangle. */
    if (retained_face_count && (!source_faces || !(overlap_tolerance > 0.0)))
	return removed_count;

    /* Open local-remesh islands may be removed only when every vertex and
     * centroid lies on retained triangles from the same source face.  This is
     * an explicit redundancy proof within the tessellation tolerance, whether
     * the retained reference is a certified closed shell or, for a partial
     * display mesh, the largest open component. */
    if (source_faces && overlap_tolerance > 0.0) {
	size_t largest = 0;
	if (!retained_face_count) {
	    for (size_t component = 1; component < components.size(); ++component) {
		if (components[component].size() > components[largest].size())
		    largest = component;
	    }
	    retain[largest] = true;
	    retained_face_count = components[largest].size();
	}
	double minimum_failed_distance = DBL_MAX;
	int failed_source_face = -1;
	const auto point_covered =
	    [&](const point_t point, int source_face) {
		double minimum_distance = DBL_MAX;
		for (int candidate = 0; candidate < *face_count; ++candidate) {
		    if (!retain[(size_t)component_id[(size_t)candidate]])
			continue;
		    if (source_faces[candidate] != source_face)
			continue;
		    point_t triangle[3];
		    for (int corner = 0; corner < 3; ++corner) {
			const int vertex =
			    faces[(size_t)candidate * 3 + corner];
			VMOVE(triangle[corner],
			    &vertices[(size_t)vertex * 3]);
		    }
		    const double distance = bg_tri_closest_pt(NULL, point,
			triangle[0], triangle[1], triangle[2]);
		    minimum_distance = std::min(minimum_distance, distance);
		    if (distance <= overlap_tolerance)
			return true;
		}
		minimum_failed_distance = minimum_distance;
		failed_source_face = source_face;
		return false;
	    };
	bool redundant = true;
	for (size_t component = 0;
		component < components.size() && redundant; ++component) {
	    if (retain[component])
		continue;
	    for (int face : components[component]) {
		point_t centroid = VINIT_ZERO;
		for (int corner = 0; corner < 3; ++corner) {
		    const int vertex = faces[(size_t)face * 3 + corner];
		    point_t point;
		    VMOVE(point, &vertices[(size_t)vertex * 3]);
		    if (!point_covered(point, source_faces[face])) {
			redundant = false;
			break;
		    }
		    VADD2(centroid, centroid, point);
		}
		VSCALE(centroid, centroid, 1.0 / 3.0);
		if (redundant &&
			!point_covered(centroid, source_faces[face]))
		    redundant = false;
		if (!redundant)
		    break;
	    }
	}
	if (!redundant) {
	    if (!emit_diagnostics)
		return removed_count;
	    bu_log("CDT retained an open disconnected component because source "
		"face %d was not redundantly covered (distance %.17g, "
		"tolerance %.17g)\n", failed_source_face,
		minimum_failed_distance, overlap_tolerance);
	    return removed_count;
	}
    }
    if (!retained_face_count)
	return removed_count;

    std::vector<int> candidate_faces;
    std::vector<int> candidate_normals;
    std::vector<int> candidate_sources;
    std::vector<size_t> candidate_triangles;
    candidate_faces.reserve(retained_face_count * 3);
    if (face_normals)
	candidate_normals.reserve(retained_face_count * 3);
    if (source_faces)
	candidate_sources.reserve(retained_face_count);
    if (source_triangles)
	candidate_triangles.reserve(retained_face_count);
    for (int face = 0; face < *face_count; ++face) {
	if (!retain[(size_t)component_id[(size_t)face]])
	    continue;
	const size_t offset = (size_t)face * 3;
	for (int corner = 0; corner < 3; ++corner) {
	    candidate_faces.push_back(faces[offset + corner]);
	    if (face_normals)
		candidate_normals.push_back(face_normals[offset + corner]);
	}
	if (source_faces)
	    candidate_sources.push_back(source_faces[face]);
	if (source_triangles)
	    candidate_triangles.push_back(source_triangles[face]);
    }
    if (!source_faces &&
	    bg_trimesh_solid2(vertex_count, (int)retained_face_count,
		vertices, candidate_faces.data(), NULL) != 0)
	return removed_count;

    const int removed = removed_count +
	*face_count - (int)retained_face_count;
    std::copy(candidate_faces.begin(), candidate_faces.end(), faces);
    if (face_normals)
	std::copy(candidate_normals.begin(), candidate_normals.end(),
	    face_normals);
    if (source_faces)
	std::copy(candidate_sources.begin(), candidate_sources.end(),
	    source_faces);
    if (source_triangles)
	std::copy(candidate_triangles.begin(), candidate_triangles.end(),
	    source_triangles);
    *face_count = (int)retained_face_count;
    return removed;
}

// TODO - get rid of all BN_TOL_DIST-only tolerances - if the object is
// very small, that distance is too big (e.g. for linearity testing).

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
