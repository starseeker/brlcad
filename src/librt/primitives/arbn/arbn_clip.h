#ifndef RT_PRIMITIVES_ARBN_CLIP_H
#define RT_PRIMITIVES_ARBN_CLIP_H

#include "common.h"
#include "vmath.h"
#include "raytrace.h"

/* Incremental half-space intersection based tessellation for ARBN solids.
 * References (succinct):
 *  - Incremental convex polyhedron clipping: Preparata & Shamos 1985; Seidel 1991.
 *  - Robust plane intersection and tolerance use: Shewchuk 1997.
 *  - Vertex/plane redundancy filtering ideas: Edelsbrunner 1987.
 */

struct rt_arbn_internal;

struct arbn_clip_vertex {
    point_t p;
    int alive;
};

struct arbn_clip_face {
    plane_t plane;
    int *vids;      /* vertex indices forming convex polygon loop (CCW) */
    int vcnt;
};

struct arbn_clip_poly {
    struct arbn_clip_vertex *verts;
    size_t vcnt;
    struct arbn_clip_face  *faces;
    size_t fcnt;
};

/* Build clipped polyhedron from ARBN planes, returns NULL on failure */
RT_EXPORT extern struct arbn_clip_poly *rt_arbn_clip_build(const struct rt_arbn_internal *aip,
                                          const struct bn_tol *tol);

/* Convert clip poly to NMG (replacement for rt_arbn_tess path) */
RT_EXPORT extern int rt_arbn_clip_to_nmg(struct nmgregion **r,
                        struct model *m,
                        const struct arbn_clip_poly *poly,
                        const struct bn_tol *tol);

/* Free clip poly resources */
RT_EXPORT extern void rt_arbn_clip_free(struct arbn_clip_poly *poly);

#endif
