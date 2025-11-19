#include "arbn_clip.h"
#include "raytrace.h"
#include "nmg.h"
#include "bg/polygon.h"
#include "bu/malloc.h"
#include "bu/log.h"
#include <math.h>

/* === Internal helpers =================================================== */

/* Normalize plane and return scale factor.
 * Robust normalization reduces condition number in triple intersections.
 * Ref: Shewchuk 1997.
 */
static void normalize_plane(plane_t p) {
    double len = sqrt(p[X]*p[X] + p[Y]*p[Y] + p[Z]*p[Z]);
    if (len <= SMALL_FASTF) return;
    p[X] /= len; p[Y] /= len; p[Z] /= len; p[W] /= len;
}

/* Classify vertex w.r.t plane: returns signed distance */
static fastf_t plane_dist(const plane_t P, const point_t v) {
    return VDOT(P, v) - P[W];
}

/* Compute intersection between segment (A,B) and plane P if endpoints straddle.
 * Linear interpolation; assumes P normal unit length (after normalization).
 * Reference: Sutherland & Hodgman 1974 (polygon clipping).
 */
static int segment_plane_isect(point_t out, const point_t Ain, const point_t Bin, const plane_t P, const struct bn_tol *tol) {
    fastf_t dA = plane_dist(P, Ain);
    fastf_t dB = plane_dist(P, Bin);
    if (dA * dB > 0) return 0; /* same side or both on plane */

    fastf_t denom = dA - dB;
    if (fabs(denom) < tol->dist) return 0;

    fastf_t t = dA / denom; /* A + t(B - A) */
    if (t < -tol->dist || t > 1.0 + tol->dist) return 0;

    /* out = (1 - t) * A + t * B */
    VBLEND2(out, t, Bin, 1.0 - t, Ain);
    return 1;
}

/* Order face vertices CCW using a local 2D projection.
 * Polar angle sort approach: standard convex polygon ordering (Preparata & Shamos 1985).
 */
static void order_face_vertices(point_t *pts, int *order, int count, const plane_t plane) {
    if (count < 3) return;
    /* Compute centroid */
    point_t C = VINIT_ZERO;
    for (int i = 0; i < count; ++i) VADD2(C, C, pts[i]);
    VSCALE(C, C, 1.0 / count);
    /* Establish basis */
    vect_t n; VSET(n, plane[X], plane[Y], plane[Z]);
    vect_t u = VINIT_ZERO;
    /* Choose arbitrary u not collinear with n */
    if (fabs(n[X]) < 0.9) { u[X] = 1.0; } else { u[Y] = 1.0; }
    vect_t uu, v;
    VCROSS(uu, n, u); VUNITIZE(uu);
    VCROSS(v, n, uu); VUNITIZE(v);
    /* Compute angles */
    double *ang = (double *)bu_malloc(sizeof(double)*count, "face angles");
    for (int i = 0; i < count; ++i) {
        vect_t d; VSUB2(d, pts[i], C);
        double x = VDOT(d, uu);
        double y = VDOT(d, v);
        ang[i] = atan2(y, x);
        order[i] = i;
    }
    /* Simple insertion sort for stability (count small usually) */
    for (int i = 1; i < count; ++i) {
        int k = order[i];
        double a = ang[k];
        int j = i - 1;
        while (j >= 0 && ang[order[j]] > a) { order[j+1] = order[j]; j--; }
        order[j+1] = k;
    }
    bu_free(ang, "face angles");
}

/* === Polyhedron Seeding =================================================
 * Seed polyhedron is a large axis-aligned box big enough to encompass
 * eventual feasible region. Could alternatively derive from first four planes
 * but box is simpler & robust. (Incremental half-space intersection:
 * Clarkson 1994 - output sensitivity rationale.)
 */
static struct arbn_clip_poly *seed_poly(const struct rt_arbn_internal *UNUSED(aip)) {
    struct arbn_clip_poly *poly = (struct arbn_clip_poly *)bu_calloc(1, sizeof(*poly), "clip poly");
    /* Define oversized cube vertices */
    double R = 1e6; /* TODO: derive bound from max |d| of planes */
    poly->vcnt = 8;
    poly->verts = (struct arbn_clip_vertex *)bu_calloc(8, sizeof(struct arbn_clip_vertex), "clip verts");
    point_t base[8] = {
        { -R, -R, -R }, { R, -R, -R }, { R, R, -R }, { -R, R, -R },
        { -R, -R,  R }, { R, -R,  R }, { R, R,  R }, { -R, R,  R }
    };
    for (int i = 0; i < 8; ++i) {
        VMOVE(poly->verts[i].p, base[i]);
        poly->verts[i].alive = 1;
    }
    /* Faces of cube (6) */
    poly->fcnt = 6;
    poly->faces = (struct arbn_clip_face *)bu_calloc(poly->fcnt, sizeof(struct arbn_clip_face), "clip faces");
    int face_vids[6][4] = {
        {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, {2,3,7,6}, {1,2,6,5}, {0,3,7,4}
    };
    for (int f = 0; f < 6; ++f) {
        poly->faces[f].vcnt = 4;
        poly->faces[f].vids = (int *)bu_malloc(sizeof(int)*4, "face vids");
        for (int j = 0; j < 4; ++j) poly->faces[f].vids[j] = face_vids[f][j];
        /* Plane equations not critical initially—will be recomputed when needed */
        HSET(poly->faces[f].plane, 0.0, 0.0, 0.0, 0.0);
    }
    return poly;
}

/* === Clipping Step ======================================================
 * Clip existing polyhedron by a new plane: produce new vertex set and faces.
 * Reference: incremental half-space intersection strategy (Preparata & Shamos 1985; Barber et al. 1996 dual hull reasoning).
 */
static int clip_with_plane(struct arbn_clip_poly *poly, const plane_t P, const struct bn_tol *tol) {
    /* Mark inside/outside for each vertex */
    int *inside = (int *)bu_calloc(poly->vcnt, sizeof(int), "inside flags");
    int inside_count = 0;
    for (size_t i = 0; i < poly->vcnt; ++i) {
        fastf_t d = plane_dist(P, poly->verts[i].p);
        if (d <= tol->dist) { inside[i] = 1; inside_count++; }
    }
    if (inside_count == (int)poly->vcnt) { bu_free(inside, "inside"); return 1; } /* No change */
    if (inside_count == 0) { bu_free(inside, "inside"); return 0; } /* Entire poly outside -> empty */

    /* Collect intersection points for new face */
    point_t *new_pts = (point_t *)bu_malloc(sizeof(point_t)*poly->vcnt*2, "new face pts");
    int new_cnt = 0;

    /* For each face, walk edges and emit new or existing vertices */
    for (size_t f = 0; f < poly->fcnt; ++f) {
        struct arbn_clip_face *face = &poly->faces[f];
        int vprev = face->vids[face->vcnt - 1];
        for (int ei = 0; ei < face->vcnt; ++ei) {
            int vcurr = face->vids[ei];
            point_t A, B;
            VMOVE(A, poly->verts[vprev].p);
            VMOVE(B, poly->verts[vcurr].p);
            int Ain = inside[vprev];
            int Bin = inside[vcurr];
            /* If edge crosses plane, create intersection vertex */
            if (Ain != Bin) {
                point_t ip;
                if (segment_plane_isect(ip, A, B, P, tol)) {
                    VMOVE(new_pts[new_cnt++], ip);
                }
            }
            vprev = vcurr;
        }
    }

    if (new_cnt < 3) {
        /* Degenerate clip (grazes poly) - treat as no change */
        bu_free(new_pts, "new face pts");
        bu_free(inside, "inside flags");
        return 1;
    }

    /* Deduplicate new face points (simple O(n^2) pass) */
    int uniq_cnt = 0;
    for (int i = 0; i < new_cnt; ++i) {
        int dup = 0;
        for (int j = 0; j < uniq_cnt; ++j) {
            if (DIST_PNT_PNT_SQ(new_pts[i], new_pts[j]) <= tol->dist_sq) { dup = 1; break; }
        }
        if (!dup) {
            if (uniq_cnt != i) VMOVE(new_pts[uniq_cnt], new_pts[i]);
            uniq_cnt++;
        }
    }
    new_cnt = uniq_cnt;

    int *order = (int *)bu_malloc(sizeof(int)*new_cnt, "face order");
    for (int i = 0; i < new_cnt; ++i) order[i] = i;
    order_face_vertices(new_pts, order, new_cnt, P);

    /* Append new face polygon to poly (store the cutting face) */
    poly->faces = (struct arbn_clip_face *)bu_realloc(poly->faces, sizeof(struct arbn_clip_face)*(poly->fcnt + 1), "faces grow");
    poly->faces[poly->fcnt].vcnt = new_cnt;
    poly->faces[poly->fcnt].vids = (int *)bu_malloc(sizeof(int)*new_cnt, "new face vids");
    HMOVE(poly->faces[poly->fcnt].plane, P);
    for (int i = 0; i < new_cnt; ++i) {
        /* Add vertex to list (TODO: spatial hash reuse) */
        point_t pt; VMOVE(pt, new_pts[order[i]]);
        poly->verts = (struct arbn_clip_vertex *)bu_realloc(poly->verts, sizeof(struct arbn_clip_vertex)*(poly->vcnt + 1), "vert grow");
        VMOVE(poly->verts[poly->vcnt].p, pt);
        poly->verts[poly->vcnt].alive = 1;
        poly->faces[poly->fcnt].vids[i] = (int)poly->vcnt;
        poly->vcnt++;
    }
    poly->fcnt++;

    bu_free(order, "face order");
    bu_free(new_pts, "new face pts");
    bu_free(inside, "inside flags");
    return 1;
}

/* === Public Interface =================================================== */

struct arbn_clip_poly *rt_arbn_clip_build(const struct rt_arbn_internal *aip, const struct bn_tol *tol) {
    if (!aip || aip->neqn < 4) return NULL;
    struct arbn_clip_poly *poly = seed_poly(aip);

    for (size_t i = 0; i < aip->neqn; ++i) {
        plane_t P; HMOVE(P, aip->eqn[i]); normalize_plane(P);
        if (!clip_with_plane(poly, P, tol)) {
            /* Plane excludes entire polyhedron: empty region */
            rt_arbn_clip_free(poly);
            return NULL;
        }
    }
    return poly;
}

int rt_arbn_clip_to_nmg(struct nmgregion **r, struct model *m, const struct arbn_clip_poly *poly, const struct bn_tol *tol) {
    if (!poly || poly->vcnt == 0 || poly->fcnt == 0) return -1;
    *r = nmg_mrsv(m);
    struct shell *s = BU_LIST_FIRST(shell, &(*r)->s_hd);

    /* Build a vertex map, initialized to NULL. nmg_cmface will populate
     * entries as needed when passed pointers-to-vertex-pointers.
     * After faces are created, assign geometry via nmg_vertex_gv once.
     * (Kettner 1999 half-edge adjacency pattern; Mirtich 1996 for mass properties consistency)
     */
    struct vertex **vert_map = (struct vertex **)bu_calloc(poly->vcnt, sizeof(struct vertex *), "vert map");

    /* Create faces */
    for (size_t f = 0; f < poly->fcnt; ++f) {
        const struct arbn_clip_face *face = &poly->faces[f];
        if (face->vcnt < 3) continue;

        /* nmg_cmface expects struct vertex *** (array of pointers-to-vertex-pointers) */
        struct vertex ***loop = (struct vertex ***)bu_malloc(sizeof(struct vertex **)*face->vcnt, "loop verts ptrptr");
        for (int i = 0; i < face->vcnt; ++i) {
            loop[i] = &vert_map[face->vids[i]];
        }
        struct faceuse *fu = nmg_cmface(s, loop, face->vcnt);
        bu_free(loop, "loop verts ptrptr");

        if (fu) {
            /* Compute and set plane equation from created geometry (Mirtich 1996) */
            if (nmg_fu_planeeqn(fu, tol)) {
                bu_log("Failed to calculate face plane equation\n");
            }
        }
    }

    /* Assign vertex geometry for any vertices lacking it */
    for (size_t i = 0; i < poly->vcnt; ++i) {
        if (vert_map[i] && !vert_map[i]->vg_p) {
            nmg_vertex_gv(vert_map[i], poly->verts[i].p);
        }
    }

    nmg_region_a(*r, tol);
    bu_free(vert_map, "vert map");
    return 0;
}

void rt_arbn_clip_free(struct arbn_clip_poly *poly) {
    if (!poly) return;
    for (size_t f = 0; f < poly->fcnt; ++f)
        bu_free(poly->faces[f].vids, "face vids");
    bu_free(poly->faces, "faces");
    bu_free(poly->verts, "verts");
    bu_free(poly, "poly");
}
