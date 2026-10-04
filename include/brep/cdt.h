/*                      C D T . H
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
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
/** @{ */
/** @file brep/cdt.h */
/** @addtogroup brep_util
 *
 * @brief
 * Constrained Delaunay Triangulation of brep solids.
 *
 */

#ifndef BREP_CDT_H
#define BREP_CDT_H

#include "common.h"

#include "bv/vlist.h"
#include "bn/tol.h"
#include "bg/defines.h"
#include "brep/defines.h"

__BEGIN_DECLS

/* Container that holds the state of a triangulation */
struct ON_Brep_CDT_State;

/* Conversion-quality tessellation result categories.  These supplement the
 * legacy integer return from ON_Brep_CDT_Tessellate with a stable reason for
 * failure. */
#define BREP_CDT_RESULT_UNATTEMPTED 0
#define BREP_CDT_RESULT_SUCCESS 1
#define BREP_CDT_RESULT_PARTIAL 2
#define BREP_CDT_RESULT_INVALID_BREP -1
#define BREP_CDT_RESULT_INVALID_TOLERANCE -2
#define BREP_CDT_RESULT_INITIALIZATION_FAILED -3
#define BREP_CDT_RESULT_FACE_FAILED -4
#define BREP_CDT_RESULT_MESH_EXPORT_FAILED -5
#define BREP_CDT_RESULT_NON_SOLID -6
#define BREP_CDT_RESULT_INVALID_PSLG -7
#define BREP_CDT_RESULT_DETRIA_FAILED -8
#define BREP_CDT_RESULT_CERTIFICATION_FAILED -9
#define BREP_CDT_RESULT_CHART_FAILED -10
#define BREP_CDT_RESULT_REFINEMENT_LIMIT -11
#define BREP_CDT_RESULT_GEOMETRIC_FAILED -12

#define BREP_CDT_STAGE_NONE 0
#define BREP_CDT_STAGE_INPUT 1
#define BREP_CDT_STAGE_TOPOLOGY 2
#define BREP_CDT_STAGE_EDGE_INITIALIZATION 3
#define BREP_CDT_STAGE_FACE_TRIANGULATION 4
#define BREP_CDT_STAGE_MESH_ASSEMBLY 5
#define BREP_CDT_STAGE_SOLID_VALIDATION 6
#define BREP_CDT_STAGE_PSLG_VALIDATION 7
#define BREP_CDT_STAGE_DETRIA 8
#define BREP_CDT_STAGE_CHART_CONSTRUCTION 9
#define BREP_CDT_STAGE_ADAPTIVE_REFINEMENT 10
#define BREP_CDT_STAGE_GEOMETRIC_VALIDATION 11

/* Create and initialize a CDT state with default tolerances.  bv
 * must be a pointer to an ON_Brep object. */
extern BREP_EXPORT struct ON_Brep_CDT_State *
ON_Brep_CDT_Create(void *bv, const char *objname);

/* Destroy a CDT state */
extern BREP_EXPORT void
ON_Brep_CDT_Destroy(struct ON_Brep_CDT_State *s);

extern BREP_EXPORT const char *
ON_Brep_CDT_ObjName(struct ON_Brep_CDT_State *s);

/* Set/get the CDT tolerances. */
extern BREP_EXPORT void
ON_Brep_CDT_Tol_Set(struct ON_Brep_CDT_State *s, const struct bg_tess_tol *t);
extern BREP_EXPORT void
ON_Brep_CDT_Tol_Get(struct bg_tess_tol *t, const struct ON_Brep_CDT_State *s);

/* Return the ON_Brep associated with state s. */
extern BREP_EXPORT void *
ON_Brep_CDT_Brep(struct ON_Brep_CDT_State *s);

/* Given a state, produce a triangulation.  Returns 0 if a solid, valid
 * triangulation was produced, 1 if a triangulation was produced but it
 * isn't solid, and -1 if no triangulation could be produced. If faces is
 * non-null, the triangulation will only attempt to triangulate the
 * specified face(s) and the return code will be the number of successfully
 * triangulated faces.  If the CDT tolerances have been updated since the
 * last Tessellate call, the old tessellation information will be replaced. */
extern BREP_EXPORT int
ON_Brep_CDT_Tessellate(struct ON_Brep_CDT_State *s, int face_cnt, int *faces);

/* Given a state, report the status of its triangulation. -3 indicates a
 * failed attempt to tessellate, -2 indicates a non-solid tessellation is
 * present after an attempt to tessellate all faces, -1 is a state which
 * has had no tessellation attempt made, 0 indicates a solid, valid full
 * brep tessellation is present, and >0 indicates that number of faces has
 * been tessellated but not the full brep. */
extern BREP_EXPORT int
ON_Brep_CDT_Status(struct ON_Brep_CDT_State *s);

/* Construct a vlist plot from the tessellation.  Modes are:
 *
 * 0 - shaded 3D triangles
 * 1 - 3D triangle wireframe
 * 2 - 2D triangle wireframe (from parametric space)
 *
 * Returns 0 if vlist was successfully generated, else -1
 */
extern BREP_EXPORT int
ON_Brep_CDT_VList(
    struct bv_vlblock *vbp,
    struct bu_list *vlfree,
    struct bu_color *c,
    int mode,
    struct ON_Brep_CDT_State *s);

/* Given two or more triangulation states, refine them to clear any face
 * overlaps introduced by the triangulation.  If any of the states are
 * un-tessellated, first perform the tessellation indicated by the state
 * settings and then proceed to resolve after all states have an initial
 * tessellation.  Returns 0 if no changes were needed, the number of
 * updated CDT states if changes were made, and -1 if one or more
 * unresolvable overlaps were encountered.  Individual CDT states may
 * subsequently be queried for other information about their specific
 * states with other function calls - this function returns only the
 * overall result. */
extern BREP_EXPORT int
ON_Brep_CDT_Ovlp_Resolve(struct ON_Brep_CDT_State **s_a, int s_cnt, double lthreshold, int timeout);

#if 0
/* Report the number of other tessellation states which manifest unresolvable
 * overlaps with state s.  If the ovlps argument is non-null, populate with
 * the problematic states.  If no resolve step was performed on s, return -1 */
extern BREP_EXPORT int
ON_Brep_CDT_UnResolvable_Ovlps(std::vector<struct ON_Brep_CDT_State *> *ovlps, struct ON_Brep_CDT_State *s);
#endif

/* Retrieve the face, vertex and normal information from a tessellation state
 * in the form of integer and fastf_t arrays. */
/* TODO - need to allow optional specification of specific faces here -
 * have already hit one scenario where I want triangle information from
 * specific faces. */
extern BREP_EXPORT int
ON_Brep_CDT_Mesh(
    int **faces, int *fcnt,
    fastf_t **vertices, int *vcnt,
    int **face_normals, int *fn_cnt,
    fastf_t **normals, int *ncnt,
    struct ON_Brep_CDT_State *s,
    int exp_face_cnt, int *exp_faces
    );

#ifdef __cplusplus
/* Original (fast but not watertight) routine used for plotting */
extern BREP_EXPORT int
brep_facecdt_plot(struct bu_vls *vls, const char *solid_name,
	const struct bg_tess_tol *ttol, const struct bn_tol *tol,
	const ON_Brep *brep, struct bu_list *p_vhead,
	struct bv_vlblock *vbp, struct bu_list *vlfree,
      	int index, int plottype, int num_points);

/* Routine to capture the triangles from the fast CDT process
 * for caching */
extern BREP_EXPORT int
brep_cdt_fast(int **faces, int *face_cnt, vect_t **pnt_norms, point_t **pnts, int *pntcnt,
	const ON_Brep *brep, int index, const struct bg_tess_tol *ttol, const struct bn_tol *tol);

/* Resource controls and diagnostics for display-quality tessellation.  Zero
 * resource and tolerance values select library defaults; adaptive_quality is
 * a boolean switch.  max_time_ms is checked
 * between faces; the per-face samplers also have fixed progress and recursion
 * guards to prevent non-terminating refinement.  face_status, when non-NULL,
 * is called exactly once per requested face during serial result assembly.
 * face_output, when non-NULL, reports the contiguous output ranges assigned
 * to each completed face with drawable geometry.  The optional trim-sample
 * callbacks replace sampling for trims where trim_sample_count returns at
 * least two points.  Samples must be ordered in trim direction and supply
 * the trim parameter, face UV coordinate, and exact 3D boundary coordinate.
 * trim_sample_source may associate an opaque identity with each supplied
 * sample.  point_source reports that identity for output points derived from
 * supplied samples, or NULL for generated points.  Its point_index is in the
 * final, concatenated output point array.  The identity is never dereferenced
 * by libbrep and need remain valid only until brep_cdt_fast_ex returns.
 * Callbacks may be invoked concurrently when max_workers exceeds one, except
 * face_status, face_diagnostic, and face_output, which run during serial
 * result assembly. */
struct brep_cdt_fast_options {
    size_t max_workers;
    size_t max_result_bytes;
    size_t max_points;
    long max_time_ms;
    int allow_partial;
    void (*face_status)(int face_index, int status, void *data);
    void *face_status_data;
    void (*face_diagnostic)(int face_index, int result, int stage,
	const char *message, void *data);
    void *face_diagnostic_data;
    void (*face_output)(int face_index, size_t first_face,
	size_t face_count, size_t first_point, size_t point_count, void *data);
    void *face_output_data;
    size_t (*trim_sample_count)(int face_index, int trim_index, void *data);
    int (*trim_sample)(int face_index, int trim_index, size_t sample_index,
	fastf_t *trim_parameter, point2d_t uv, point_t point, void *data);
    void *trim_sample_data;
    const void *(*trim_sample_source)(int face_index, int trim_index,
	size_t sample_index, void *data);
    void (*point_source)(int face_index, size_t point_index,
	const void *source, void *data);
    void *point_source_data;
    /* Shared estimated transient-memory allowance for concurrent face jobs.
     * This is independent of max_result_bytes, which bounds retained output.
     * A zero value selects an availability-calibrated library default. */
    size_t max_working_bytes;
    /* Display triangle target.  Adaptive whole-object budgeting may select a
     * lower target according to B-Rep topology.  Authoritative trim
     * boundaries are retained and may make the final count exceed it. */
    size_t max_triangles;
    /* Generate a complete coarse mesh before refining toward the requested
     * tolerance.  Intended for visual display; rigorous callers may disable
     * it to request the specified tolerance directly. */
    int adaptive_quality;
    /* Initial relative display tolerance and per-face unsigned area-change
     * convergence threshold.  Zero values select library defaults. */
    double coarse_relative_tolerance;
    double area_change_tolerance;
    /* Retain validated pcurve-repair samples for cross-face mesh assembly.
     * Display callers normally simplify them within the repair tolerance. */
    int preserve_pullback_samples;
};

#define BREP_CDT_FAST_FACE_COMPLETED 0
#define BREP_CDT_FAST_FACE_FAILED 1
#define BREP_CDT_FAST_FACE_SKIPPED_DEGENERATE 2
#define BREP_CDT_FAST_FACE_NOT_PROCESSED 3
#define BREP_CDT_FAST_FACE_APPROXIMATED 4
#define BREP_CDT_FAST_FACE_SKIPPED_TOLERANCE 5

/* completed_faces includes faces with no drawable output.  Exact zero-area
 * faces and faces which collapse at the requested display tolerance are
 * reported separately. */
struct brep_cdt_fast_report {
    int requested_faces;
    int completed_faces;
    int failed_faces;
    size_t result_bytes;
    int hit_time_limit;
    int hit_memory_limit;
    int hit_point_limit;
    /* Completed faces proven to have exactly zero parametric area. */
    int skipped_degenerate_faces;
    /* Completed faces with no boundary resolvable at the requested
     * tessellation and model tolerances. */
    int skipped_tolerance_faces;
    /* Peak sum of conservative face-work reservations. */
    size_t peak_working_bytes;
    /* Adaptive display-quality provenance. */
    size_t triangle_budget;
    int approximated_faces;
    int area_converged_faces;
    int triangle_budget_limited_faces;
    int refinement_passes;
    int refinement_time_limited;
    /* Faces whose realized triangles do not cover the sampled B-Rep edge
     * envelope within the requested display tolerance. */
    int boundary_envelope_incomplete_faces;
};

#define BREP_CDT_FAST_OK 0
#define BREP_CDT_FAST_PARTIAL 1
#define BREP_CDT_FAST_ERROR -1
#define BREP_CDT_FAST_LIMIT -2

/**
 * Initialize display-quality CDT options with library defaults.
 *
 * A NULL options pointer is ignored.  Callers may override any nonzero
 * resource limit or other option after initialization.
 *
 * @param options options structure to initialize
 */
extern BREP_EXPORT void
brep_cdt_fast_options_default(struct brep_cdt_fast_options *options);

/**
 * Tessellate one B-Rep face, or all faces when @p index is -1, using the
 * bounded display-quality CDT path.
 *
 * Output arrays are allocated with the BRL-CAD allocator and are owned by
 * the caller.  They are initialized to NULL/zero before any validation or
 * processing.  The optional report receives resource and per-face status
 * information.  A partial result is returned only when the options permit
 * partial output.  An all-faces request for a B-Rep with no faces succeeds
 * with empty output; selecting an individual nonexistent face is an error.
 *
 * Failed trim boundaries may be retried after bounded topology healing on
 * an owned copy.  The source is unchanged and no holes are capped.  This
 * recovery is disabled when authoritative trim samples, source callbacks,
 * or preserve_pullback_samples require the original topology identities.
 *
 * @param faces receives triangle vertex indices
 * @param face_cnt receives the number of triangles
 * @param pnt_norms receives three face-oriented corner normals per triangle,
 *        in the same order as faces; smooth normals opposing a realized
 *        nondegenerate facet are replaced by its geometric normal
 * @param pnts receives output vertices
 * @param pntcnt receives the number of output vertices
 * @param brep source B-Rep
 * @param index face index, or -1 for all faces
 * @param ttol tessellation tolerance
 * @param tol geometric tolerance
 * @param options optional resource and callback settings; NULL selects
 *        defaults
 * @param report optional diagnostic report
 * @return BREP_CDT_FAST_OK for complete output, BREP_CDT_FAST_PARTIAL for
 *         permitted partial output, BREP_CDT_FAST_LIMIT for a bounded
 *         resource stop, or BREP_CDT_FAST_ERROR for invalid input or other
 *         failure
 */
extern BREP_EXPORT int
brep_cdt_fast_ex(int **faces, int *face_cnt, vect_t **pnt_norms,
	point_t **pnts, int *pntcnt, const ON_Brep *brep, int index,
	const struct bg_tess_tol *ttol, const struct bn_tol *tol,
	const struct brep_cdt_fast_options *options,
	struct brep_cdt_fast_report *report);
#endif

__END_DECLS

/** @} */

#endif  /* BREP_CDT_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
