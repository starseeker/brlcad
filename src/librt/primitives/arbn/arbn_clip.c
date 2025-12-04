#include "common.h"

#include "arbn_clip.h"
#include "raytrace.h"
#include "nmg.h"
#include "bg/polygon.h"
#include "bu/malloc.h"
#include "bu/log.h"
#include "bu/env.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* === Configuration via environment variables ============================ */

/* Check if spatial hashing is enabled (default: on - bug #16 fixed) */
static int use_spatial_hash(void) {
    const char *env = getenv("BRLCAD_ARBN_CLIP_SPATIAL_HASH");
    if (!env) return 1; /* default enabled - bug #16 fixed with hash rebuild */
    return !BU_STR_EQUAL(env, "off") && !BU_STR_EQUAL(env, "0");
}

/* Get maximum plane count limit (default: 10000) */
static size_t get_max_planes(void) {
    const char *env = getenv("BRLCAD_ARBN_CLIP_MAX_PLANES");
    if (!env) return 10000;
    long val = atol(env);
    return (val > 0) ? (size_t)val : 10000;
}

/* === Spatial hash for vertex deduplication ==============================
 * Grid-based spatial hashing (Teschner et al. 2003) for O(1) vertex lookup.
 */

#define SPATIAL_HASH_BINS 1024

struct spatial_hash_entry {
    int vertex_id;
    struct spatial_hash_entry *next;
};

struct spatial_hash {
    struct spatial_hash_entry *bins[SPATIAL_HASH_BINS];
    double cell_size;
    const struct bn_tol *tol;
};

static struct spatial_hash *spatial_hash_create(double cell_size, const struct bn_tol *tol) {
    struct spatial_hash *hash = (struct spatial_hash *)bu_calloc(1, sizeof(*hash), "spatial hash");
    hash->cell_size = cell_size;
    hash->tol = tol;
    return hash;
}

static void spatial_hash_free(struct spatial_hash *hash) {
    if (!hash) return;
    for (int i = 0; i < SPATIAL_HASH_BINS; ++i) {
        struct spatial_hash_entry *e = hash->bins[i];
        while (e) {
            struct spatial_hash_entry *next = e->next;
            bu_free(e, "hash entry");
            e = next;
        }
    }
    bu_free(hash, "spatial hash");
}

static size_t spatial_hash_key(const point_t p, double cell_size) {
    /* Quantize to grid cells (compute first, then cast to avoid -Wbad-function-cast) */
    double fx = floor(p[X] / cell_size);
    double fy = floor(p[Y] / cell_size);
    double fz = floor(p[Z] / cell_size);
    long ix = (long)fx;
    long iy = (long)fy;
    long iz = (long)fz;
    size_t h = (size_t)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791));
    return h % SPATIAL_HASH_BINS;
}

/* Find existing vertex within tolerance, or return -1 */
static int spatial_hash_find(struct spatial_hash *hash, const point_t p, const struct arbn_clip_poly *poly) {
    size_t key = spatial_hash_key(p, hash->cell_size);
    struct spatial_hash_entry *e = hash->bins[key];
    while (e) {
        if (e->vertex_id < (int)poly->vcnt && poly->verts[e->vertex_id].alive) {
            if (DIST_PNT_PNT_SQ(p, poly->verts[e->vertex_id].p) <= hash->tol->dist_sq) {
                return e->vertex_id;
            }
        }
        e = e->next;
    }
    return -1;
}

/* Add vertex to hash */
static void spatial_hash_insert(struct spatial_hash *hash, const point_t p, int vertex_id) {
    size_t key = spatial_hash_key(p, hash->cell_size);
    struct spatial_hash_entry *e = (struct spatial_hash_entry *)bu_malloc(sizeof(*e), "hash entry");
    e->vertex_id = vertex_id;
    e->next = hash->bins[key];
    hash->bins[key] = e;
}

/* Clear all entries from the hash (but keep the structure) */
static void spatial_hash_clear(struct spatial_hash *hash) {
    if (!hash) return;
    for (int i = 0; i < SPATIAL_HASH_BINS; ++i) {
        struct spatial_hash_entry *e = hash->bins[i];
        while (e) {
            struct spatial_hash_entry *next = e->next;
            bu_free(e, "hash entry");
            e = next;
        }
        hash->bins[i] = NULL;
    }
}

/* Rebuild hash to only contain alive vertices from poly */
static void spatial_hash_rebuild(struct spatial_hash *hash, const struct arbn_clip_poly *poly) {
    if (!hash) return;
    
    /* Clear existing entries */
    spatial_hash_clear(hash);
    
    /* Re-insert only alive vertices */
    for (size_t i = 0; i < poly->vcnt; ++i) {
        if (poly->verts[i].alive) {
            spatial_hash_insert(hash, poly->verts[i].p, (int)i);
        }
    }
}

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

/* Check if two planes are duplicates (same normal direction and offset within tolerance) */
static int planes_duplicate(const plane_t p1, const plane_t p2, const struct bn_tol *tol) {
    /* Check if normals point in the same direction (not just parallel) */
    vect_t n1, n2;
    VSET(n1, p1[X], p1[Y], p1[Z]);
    VSET(n2, p2[X], p2[Y], p2[Z]);
    double dot = VDOT(n1, n2);
    if (dot < (1.0 - 1e-6)) return 0;  /* Not pointing in same direction */
    
    /* Check if offsets are the same */
    double d_diff = fabs(p1[W] - p2[W]);
    return (d_diff < tol->dist);
}

/* Preprocess planes: normalize, deduplicate, estimate bounds
 * Reference: Plane filtering approaches from Edelsbrunner 1987.
 */
static plane_t *preprocess_planes(const struct rt_arbn_internal *aip,
                                  const struct bn_tol *tol,
                                  size_t *out_count,
                                  struct arbn_clip_stats *stats) {
    size_t max_planes = get_max_planes();
    if (aip->neqn > max_planes) {
        bu_log("ARBN plane count (%zu) exceeds limit (%zu)\n", aip->neqn, max_planes);
        return NULL;
    }

    plane_t *planes = (plane_t *)bu_calloc(aip->neqn, sizeof(plane_t), "normalized planes");
    for (size_t i = 0; i < aip->neqn; ++i) {
        HMOVE(planes[i], aip->eqn[i]);
        normalize_plane(planes[i]);
    }

    int *keep = (int *)bu_calloc(aip->neqn, sizeof(int), "keep flags");
    size_t unique_count = 0;
    for (size_t i = 0; i < aip->neqn; ++i) {
        keep[i] = 1;
        for (size_t j = 0; j < i; ++j) {
            if (keep[j] && planes_duplicate(planes[i], planes[j], tol)) {
                keep[i] = 0;
                if (stats) stats->duplicate_planes++;
                break;
            }
        }
        if (keep[i]) unique_count++;
    }

    plane_t *unique_planes = (plane_t *)bu_calloc(unique_count, sizeof(plane_t), "unique planes");
    size_t idx = 0;
    for (size_t i = 0; i < aip->neqn; ++i) {
        if (keep[i]) {
            HMOVE(unique_planes[idx++], planes[i]);
        }
    }

    bu_free(keep, "keep flags");
    bu_free(planes, "normalized planes");
    *out_count = unique_count;
    return unique_planes;
}

/* Compute bounding radius from plane offsets (tight initial bound)
 * Reference: Preparata & Shamos 1985 - geometric bounds from linear constraints
 */
static double compute_bounding_radius(plane_t *planes, size_t nplanes) {
    if (nplanes == 0) return 1e6;
    double max_offset = 0.0;
    for (size_t i = 0; i < nplanes; ++i) {
        double offset = fabs(planes[i][W]);
        if (offset > max_offset) max_offset = offset;
    }
    return (max_offset > SMALL_FASTF) ? (2.0 * max_offset) : 1e6;
}

/* Classify vertex w.r.t plane: returns signed distance */
static fastf_t plane_dist(const plane_t P, const point_t v) {
    return VDOT(P, v) - P[W];
}

/* Compute intersection between segment (A,B) and plane P if endpoints straddle.
 * Linear interpolation; assumes P normal unit length.
 * Reference: Sutherland & Hodgman 1974.
 */
static int segment_plane_isect(point_t out, const point_t Ain, const point_t Bin, const plane_t P, const struct bn_tol *tol) {
    fastf_t dA = plane_dist(P, Ain);
    fastf_t dB = plane_dist(P, Bin);
    if (dA * dB > 0) return 0;

    fastf_t denom = dA - dB;
    if (fabs(denom) < tol->dist) return 0;

    fastf_t t = dA / denom;
    if (t < -tol->dist || t > 1.0 + tol->dist) return 0;

    VBLEND2(out, t, Bin, 1.0 - t, Ain);
    return 1;
}

/* Order face vertices CCW using a local 2D projection.
 * Reference: Preparata & Shamos 1985 (polar angle sort).
 */
static void order_face_vertices(point_t *pts, int *order, int count, const plane_t plane) {
    if (count < 3) return;
    point_t C = VINIT_ZERO;
    for (int i = 0; i < count; ++i) VADD2(C, C, pts[i]);
    VSCALE(C, C, 1.0 / count);
    vect_t n; VSET(n, plane[X], plane[Y], plane[Z]);
    vect_t u = VINIT_ZERO;
    if (fabs(n[X]) < 0.9) u[X] = 1.0; else u[Y] = 1.0;
    vect_t uu, v;
    VCROSS(uu, n, u); VUNITIZE(uu);
    VCROSS(v, n, uu); VUNITIZE(v);
    double *ang = (double *)bu_malloc(sizeof(double)*count, "face angles");
    for (int i = 0; i < count; ++i) {
        vect_t d; VSUB2(d, pts[i], C);
        double x = VDOT(d, uu);
        double y = VDOT(d, v);
        ang[i] = atan2(y, x);
        order[i] = i;
    }
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
 * Seed polyhedron is an axis-aligned box based on computed bound.
 * Reference: Incremental half-space intersection (Preparata & Shamos 1985;
 * Clarkson 1994).
 */
static struct arbn_clip_poly *seed_poly(double bounding_radius) {
    struct arbn_clip_poly *poly = (struct arbn_clip_poly *)bu_calloc(1, sizeof(*poly), "clip poly");
    double R = bounding_radius;
    poly->vcnt = 8;
    poly->vcap = 8;
    poly->verts = (struct arbn_clip_vertex *)bu_calloc(8, sizeof(struct arbn_clip_vertex), "clip verts");
    point_t base[8] = {
        { -R, -R, -R }, { R, -R, -R }, { R, R, -R }, { -R, R, -R },
        { -R, -R,  R }, { R, -R,  R }, { R, R,  R }, { -R, R,  R }
    };
    for (int i = 0; i < 8; ++i) {
        VMOVE(poly->verts[i].p, base[i]);
        poly->verts[i].alive = 1;
    }
    poly->fcnt = 6;
    poly->fcap = 6;
    poly->faces = (struct arbn_clip_face *)bu_calloc(6, sizeof(struct arbn_clip_face), "clip faces");
    int face_vids[6][4] = {
        {0,1,2,3}, {4,5,6,7}, {0,1,5,4}, {2,3,7,6}, {1,2,6,5}, {0,3,7,4}
    };
    for (int f = 0; f < 6; ++f) {
        poly->faces[f].vcnt = 4;
        poly->faces[f].vids = (int *)bu_malloc(sizeof(int)*4, "face vids");
        for (int j = 0; j < 4; ++j) poly->faces[f].vids[j] = face_vids[f][j];
        poly->faces[f].alive = 1;
        HSET(poly->faces[f].plane, 0.0, 0.0, 0.0, 0.0);
    }
    return poly;
}

/* === Clipping Step ======================================================
 * Clip existing polyhedron by a new plane: update faces and add new cut face.
 * References: Incremental clipping (Preparata & Shamos 1985; Barber et al. 1996).
 */
static int clip_with_plane(struct arbn_clip_poly *poly, const plane_t P,
                           const struct bn_tol *tol, struct spatial_hash *hash) {
    size_t original_vcnt = poly->vcnt;  /* Save original vertex count before adding new ones */
    int *inside = (int *)bu_calloc(poly->vcnt, sizeof(int), "inside flags");
    int inside_count = 0;
    for (size_t i = 0; i < poly->vcnt; ++i) {
        if (!poly->verts[i].alive) continue;
        fastf_t d = plane_dist(P, poly->verts[i].p);
        if (d <= tol->dist) { inside[i] = 1; inside_count++; }
    }

    int total_alive = 0;
    for (size_t i = 0; i < poly->vcnt; ++i) if (poly->verts[i].alive) total_alive++;

    if (inside_count == total_alive) {
        bu_free(inside, "inside"); return 1;
    }
    if (inside_count == 0) {
        bu_free(inside, "inside"); return 0;
    }

    /* Allocate enough space for new intersection points - worst case is each face edge can create an intersection */
    size_t max_new_pts = poly->fcnt * 20;  /* Conservative estimate: 20 intersections per face */
    point_t *new_pts = (point_t *)bu_malloc(sizeof(point_t) * max_new_pts, "new face pts");
    int new_cnt = 0;

    for (size_t f = 0; f < poly->fcnt; ++f) {
        struct arbn_clip_face *face = &poly->faces[f];
        if (!face->alive || face->vcnt < 3) continue;

        int all_in = 1, all_out = 1;
        for (int i = 0; i < face->vcnt; ++i) {
            int vid = face->vids[i];
            if (inside[vid]) all_out = 0;
            else all_in = 0;
        }

        if (all_out) {
            /* Free the vids since this face is being removed */
            bu_free(face->vids, "dead face vids");
            face->vids = NULL;
            face->alive = 0;
            continue;
        }
        if (all_in) continue;

        int vprev = face->vids[face->vcnt - 1];
        int *new_vids = (int *)bu_malloc(sizeof(int) * face->vcnt * 2, "new face vids");
        int new_vid_cnt = 0;

        for (int ei = 0; ei < face->vcnt; ++ei) {
            int vcurr = face->vids[ei];
            int Ain = inside[vprev];
            int Bin = inside[vcurr];

            if (Bin) {
                new_vids[new_vid_cnt++] = vcurr;
            }

            if (Ain != Bin) {
                point_t ip;
                if (segment_plane_isect(ip, poly->verts[vprev].p, poly->verts[vcurr].p, P, tol)) {
                    if ((size_t)new_cnt >= max_new_pts) {
                        bu_log("ERROR: new_cnt (%d) >= max_new_pts (%zu)\n", new_cnt, max_new_pts);
                        bu_free(new_vids, "overflow face vids");
                        bu_free(new_pts, "overflow pts");
                        bu_free(inside, "overflow inside");
                        return 0;
                    }
                    VMOVE(new_pts[new_cnt], ip);
                    /* Use -(new_cnt+1) as placeholder to track which new_pts entry this refers to */
                    new_vids[new_vid_cnt++] = -(new_cnt + 1);
                    new_cnt++;
                }
            }
            vprev = vcurr;
        }

        if (new_vid_cnt >= 3) {
            bu_free(face->vids, "old face vids");
            face->vids = new_vids;
            face->vcnt = new_vid_cnt;
        } else {
            bu_free(new_vids, "new face vids");
            /* Keep the old vids since we're marking face as dead */
            bu_free(face->vids, "dead face vids");
            face->vids = NULL;
            face->alive = 0;
        }
    }

    if (new_cnt < 3) {
        bu_free(new_pts, "new face pts");
        bu_free(inside, "inside flags");
        return 1;
    }

    int *new_vids_map = (int *)bu_malloc(sizeof(int) * new_cnt, "new vids map");
    int original_new_cnt = new_cnt;  /* Save original count before deduplication */
    int uniq_cnt = 0;

    for (int i = 0; i < new_cnt; ++i) {
        int existing_vid = -1;
        if (hash) {
            existing_vid = spatial_hash_find(hash, new_pts[i], poly);
        } else {
            /* First check against existing vertices from previous clipping operations */
            for (size_t j = 0; j < poly->vcnt; ++j) {
                if (poly->verts[j].alive && 
                    DIST_PNT_PNT_SQ(new_pts[i], poly->verts[j].p) <= tol->dist_sq) {
                    existing_vid = (int)j;
                    break;
                }
            }
            /* Then check against new points being added in this operation */
            if (existing_vid < 0) {
                for (int j = 0; j < uniq_cnt; ++j) {
                    if (DIST_PNT_PNT_SQ(new_pts[i], new_pts[j]) <= tol->dist_sq) {
                        existing_vid = new_vids_map[j];
                        break;
                    }
                }
            }
        }

        if (existing_vid >= 0) {
            new_vids_map[i] = existing_vid;
        } else {
            if (poly->vcnt >= poly->vcap) {
                poly->vcap = poly->vcap * 2 + 8;
                poly->verts = (struct arbn_clip_vertex *)bu_realloc(poly->verts,
                    sizeof(struct arbn_clip_vertex) * poly->vcap, "vert grow");
            }
            VMOVE(poly->verts[poly->vcnt].p, new_pts[i]);
            poly->verts[poly->vcnt].alive = 1;
            new_vids_map[i] = (int)poly->vcnt;

            if (hash) spatial_hash_insert(hash, new_pts[i], (int)poly->vcnt);

            if (uniq_cnt != i) {
                VMOVE(new_pts[uniq_cnt], new_pts[i]);
            }
            uniq_cnt++;
            poly->vcnt++;
        }
    }

    new_cnt = uniq_cnt;

    if (new_cnt < 3) {
        bu_free(new_vids_map, "new vids map");
        bu_free(new_pts, "new face pts");
        bu_free(inside, "inside flags");
        return 1;
    }

    int *order = (int *)bu_malloc(sizeof(int)*new_cnt, "face order");
    for (int i = 0; i < new_cnt; ++i) order[i] = i;
    order_face_vertices(new_pts, order, new_cnt, P);

    if (poly->fcnt >= poly->fcap) {
        size_t old_fcap = poly->fcap;
        poly->fcap = poly->fcap * 2 + 6;
        poly->faces = (struct arbn_clip_face *)bu_realloc(poly->faces,
            sizeof(struct arbn_clip_face) * poly->fcap, "faces grow");
        /* Initialize new face structures */
        for (size_t i = old_fcap; i < poly->fcap; i++) {
            poly->faces[i].vids = NULL;
            poly->faces[i].vcnt = 0;
            poly->faces[i].alive = 0;
        }
    }

    poly->faces[poly->fcnt].vcnt = new_cnt;
    poly->faces[poly->fcnt].vids = (int *)bu_malloc(sizeof(int)*new_cnt, "new face vids");
    poly->faces[poly->fcnt].alive = 1;
    HMOVE(poly->faces[poly->fcnt].plane, P);

    for (int i = 0; i < new_cnt; ++i) {
        poly->faces[poly->fcnt].vids[i] = new_vids_map[order[i]];
    }
    poly->fcnt++;

    /* Replace placeholders with actual vertex IDs after deduplication
     * Placeholders were created as -(index+1) where index is in [0, original_new_cnt)
     */
    for (size_t f = 0; f < poly->fcnt - 1; ++f) {
        struct arbn_clip_face *face = &poly->faces[f];
        if (!face->alive) continue;
        for (int i = 0; i < face->vcnt; ++i) {
            if (face->vids[i] < 0) {
                /* Placeholder is -(index+1), so recover the original index */
                int original_idx = -(face->vids[i] + 1);
                if (original_idx >= 0 && original_idx < original_new_cnt) {
                    face->vids[i] = new_vids_map[original_idx];
                } else {
                    bu_log("ERROR: Invalid placeholder index %d (valid range: 0-%d)\n", 
                           original_idx, original_new_cnt-1);
                    face->vids[i] = 0;  /* Fallback to vertex 0 */
                }
            }
        }
    }

    bu_free(new_vids_map, "new vids map");
    bu_free(order, "face order");
    bu_free(new_pts, "new face pts");
    
    /* Mark original vertices outside the clipping plane as dead.
     * Only check vertices that existed before clipping - new intersection 
     * vertices are always on the plane (inside by definition). */
    for (size_t i = 0; i < original_vcnt; ++i) {
        if (poly->verts[i].alive && !inside[i]) {
            poly->verts[i].alive = 0;
        }
    }
    
    /* Rebuild spatial hash to remove dead vertices - this fixes bug #16
     * where orphaned vertices from previous operations caused incorrect matches */
    if (hash) {
        spatial_hash_rebuild(hash, poly);
    }
    
    bu_free(inside, "inside flags");
    return 1;
}

/* === Public Interface =================================================== */

struct arbn_clip_poly *rt_arbn_clip_build_with_stats(const struct rt_arbn_internal *aip,
                                                     const struct bn_tol *tol,
                                                     struct arbn_clip_stats *stats) {
    if (!aip || aip->neqn < 4) return NULL;

    if (stats) {
        memset(stats, 0, sizeof(*stats));
        stats->input_planes = aip->neqn;
        stats->spatial_hash_enabled = use_spatial_hash();
    }

    size_t nplanes = 0;
    plane_t *planes = preprocess_planes(aip, tol, &nplanes, stats);
    if (!planes) return NULL;

    double bounding_radius = compute_bounding_radius(planes, nplanes);
    if (stats) stats->bounding_radius = bounding_radius;

    struct arbn_clip_poly *poly = seed_poly(bounding_radius);

    struct spatial_hash *hash = NULL;
    if (use_spatial_hash()) {
        hash = spatial_hash_create(tol->dist * 10.0, tol);
    }

    for (size_t i = 0; i < nplanes; ++i) {
        if (!clip_with_plane(poly, planes[i], tol, hash)) {
            if (hash) spatial_hash_free(hash);
            rt_arbn_clip_free(poly);
            bu_free(planes, "unique planes");
            return NULL;
        }
    }

    if (hash) spatial_hash_free(hash);
    bu_free(planes, "unique planes");

    /* Mark vertices as dead if they're not referenced by any alive face */
    int *vertex_used = (int *)bu_calloc(poly->vcnt, sizeof(int), "vertex used flags");
    for (size_t f = 0; f < poly->fcnt; ++f) {
        if (!poly->faces[f].alive) continue;
        for (int v = 0; v < poly->faces[f].vcnt; ++v) {
            int vid = poly->faces[f].vids[v];
            if (vid < 0 || (size_t)vid >= poly->vcnt) {
                /* Invalid vertex ID indicates data corruption - should not happen
                 * in correct implementation. Log at debug level for development. */
                if (vid < 0) {
                    bu_log("rt_arbn_clip: face %zu has placeholder vertex ID %d\n", f, vid);
                } else {
                    bu_log("rt_arbn_clip: face %zu has out-of-bounds vertex ID %d (max %zu)\n", 
                           f, vid, poly->vcnt - 1);
                }
                continue;
            }
            vertex_used[vid] = 1;
        }
    }
    for (size_t i = 0; i < poly->vcnt; ++i) {
        if (!vertex_used[i]) {
            poly->verts[i].alive = 0;
        }
    }
    bu_free(vertex_used, "vertex used flags");

    if (stats) {
        stats->active_planes   = nplanes;
        stats->final_vertices  = 0;
        stats->final_faces     = 0;
        for (size_t i = 0; i < poly->vcnt; ++i)
            if (poly->verts[i].alive) stats->final_vertices++;
        for (size_t i = 0; i < poly->fcnt; ++i)
            if (poly->faces[i].alive) stats->final_faces++;
    }

    const char *stats_path = getenv("BRLCAD_ARBN_CLIP_STATS");
    if (stats_path && stats) {
        FILE *fp = fopen(stats_path, "a");
        if (fp) {
            fprintf(fp, "ARBN Clipping Statistics:\n");
            fprintf(fp, "  Input planes: %zu\n", stats->input_planes);
            fprintf(fp, "  Duplicate planes removed: %zu\n", stats->duplicate_planes);
            fprintf(fp, "  Redundant planes removed: %zu\n", stats->redundant_planes);
            fprintf(fp, "  Active planes: %zu\n", stats->active_planes);
            fprintf(fp, "  Final vertices: %zu\n", stats->final_vertices);
            fprintf(fp, "  Final faces: %zu\n", stats->final_faces);
            fprintf(fp, "  Bounding radius: %.6f\n", stats->bounding_radius);
            fprintf(fp, "  Spatial hash: %s\n", stats->spatial_hash_enabled ? "enabled" : "disabled");
            fprintf(fp, "\n");
            fclose(fp);
        }
    }

    return poly;
}

struct arbn_clip_poly *rt_arbn_clip_build(const struct rt_arbn_internal *aip, const struct bn_tol *tol) {
    return rt_arbn_clip_build_with_stats(aip, tol, NULL);
}

int rt_arbn_clip_to_nmg(struct nmgregion **r, struct model *m, const struct arbn_clip_poly *poly, const struct bn_tol *tol) {
    if (!poly || poly->vcnt == 0 || poly->fcnt == 0) return -1;
    *r = nmg_mrsv(m);
    struct shell *s = BU_LIST_FIRST(shell, &(*r)->s_hd);

    struct vertex **vert_map = (struct vertex **)bu_calloc(poly->vcnt, sizeof(struct vertex *), "vert map");

    for (size_t f = 0; f < poly->fcnt; ++f) {
        const struct arbn_clip_face *face = &poly->faces[f];
        if (!face->alive || face->vcnt < 3) continue;

        struct vertex ***loop = (struct vertex ***)bu_malloc(sizeof(struct vertex **)*face->vcnt, "loop verts ptrptr");
        for (int i = 0; i < face->vcnt; ++i) {
            loop[i] = &vert_map[face->vids[i]];
        }
        struct faceuse *fu = nmg_cmface(s, loop, face->vcnt);
        bu_free(loop, "loop verts ptrptr");

        if (fu) {
            if (nmg_fu_planeeqn(fu, tol)) {
                bu_log("Failed to calculate face plane equation\n");
            }
        }
    }

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
    for (size_t f = 0; f < poly->fcnt; ++f) {
        if (poly->faces[f].vids) {
            bu_free(poly->faces[f].vids, "face vids");
        }
    }
    bu_free(poly->faces, "faces");
    bu_free(poly->verts, "verts");
    bu_free(poly, "poly");
}