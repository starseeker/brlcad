/*                    B V / V I E W . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
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
/** @file libbv/view.c
 *
 * Scene-graph-free view mechanics.
 */

#include "common.h"

#include <math.h>
#include <string.h>

#include "vmath.h"
#include "bg/plane.h"
#include "bn/mat.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bv/view.h"

static void
bv_mat_aet(struct bv *v)
{
    mat_t tmat;
    fastf_t twist;
    fastf_t c_twist;
    fastf_t s_twist;

    if (!v)
	return;

    bn_mat_angles(v->rotation, 270.0 + v->aet[Y], 0.0, 270.0 - v->aet[X]);
    twist = -v->aet[Z] * DEG2RAD;
    c_twist = cos(twist);
    s_twist = sin(twist);
    bn_mat_zrot(tmat, s_twist, c_twist);
    bn_mat_mul2(tmat, v->rotation);
}

struct bv *
bv_create(void)
{
    struct bv *v = NULL;
    BU_ALLOC(v, struct bv);
    bv_init(v);
    return v;
}

void
bv_destroy(struct bv *v)
{
    if (!v)
	return;

    bv_free(v);
    bu_free(v, "bv");
}

void
bv_init(struct bv *v)
{
    if (!v)
	return;

    memset(v, 0, sizeof(*v));
    v->magic = BV_STATE_MAGIC;
    bu_vls_init(&v->name);
    v->scale = BV_DEFAULT_SCALE;
    v->initial_scale = v->scale;
    v->absolute_scale = 0.0;
    v->size = 2.0 * v->scale;
    v->inverse_size = 1.0 / v->size;
    v->local2base = 1.0;
    v->base2local = 1.0;
    v->refresh_enabled = 1;
    VSET(v->aet, 35.0, 25.0, 0.0);
    VSET(v->eye_pos, 0.0, 0.0, 1.0);
    v->coord_mode = 'v';
    v->rotate_about = 'v';
    v->min_mouse_delta = -20.0;
    v->max_mouse_delta = 20.0;
    v->rotate_scale = 0.4;
    v->scale_scale = 2.0;
    MAT_IDN(v->rotation);
    MAT_IDN(v->center);
    MAT_IDN(v->model2view);
    MAT_IDN(v->pmodel2view);
    MAT_IDN(v->view2model);
    MAT_IDN(v->pmat);
    v->radius = 1.0;
    bv_knobs_reset(&v->knobs, BV_KNOBS_ALL);
    bv_mat_aet(v);
    bv_update(v);
}

void
bv_free(struct bv *v)
{
    if (!v)
	return;

    bu_vls_free(&v->name);
    v->magic = 0;
}

int
bv_is_valid(const struct bv *v)
{
    return (v && v->magic == BV_STATE_MAGIC);
}

int
bv_copy(struct bv *dst, const struct bv *src)
{
    if (!dst || !bv_is_valid(src))
	return 0;

    if (dst == src)
	return 1;

    if (bv_is_valid(dst))
	bv_free(dst);

    bv_init(dst);
    bu_vls_sprintf(&dst->name, "%s", bu_vls_cstr(&src->name));
    dst->user_data = src->user_data;
    dst->width = src->width;
    dst->height = src->height;
    dst->scale = src->scale;
    dst->initial_scale = src->initial_scale;
    dst->absolute_scale = src->absolute_scale;
    dst->size = src->size;
    dst->inverse_size = src->inverse_size;
    dst->perspective = src->perspective;
    dst->local2base = src->local2base;
    dst->base2local = src->base2local;
    dst->frame_revision = src->frame_revision;
    dst->refresh_dirty = src->refresh_dirty;
    dst->refresh_enabled = src->refresh_enabled;
    dst->refresh_suppressed = src->refresh_suppressed;
    dst->refresh_drawn_count = src->refresh_drawn_count;
    dst->frametime = src->frametime;
    dst->zclip = src->zclip;
    dst->framebuffer_mode = src->framebuffer_mode;
    dst->cleared = src->cleared;
    VMOVE(dst->aet, src->aet);
    VMOVE(dst->eye_pos, src->eye_pos);
    VMOVE(dst->keypoint, src->keypoint);
    VMOVE(dst->current_point, src->current_point);
    dst->mouse_x = src->mouse_x;
    dst->mouse_y = src->mouse_y;
    dst->previous_mouse_x = src->previous_mouse_x;
    dst->previous_mouse_y = src->previous_mouse_y;
    dst->coord_mode = src->coord_mode;
    dst->rotate_about = src->rotate_about;
    dst->min_mouse_delta = src->min_mouse_delta;
    dst->max_mouse_delta = src->max_mouse_delta;
    dst->rotate_scale = src->rotate_scale;
    dst->scale_scale = src->scale_scale;
    MAT_COPY(dst->rotation, src->rotation);
    MAT_COPY(dst->center, src->center);
    MAT_COPY(dst->model2view, src->model2view);
    MAT_COPY(dst->pmodel2view, src->pmodel2view);
    MAT_COPY(dst->view2model, src->view2model);
    MAT_COPY(dst->pmat, src->pmat);
    VMOVE(dst->obb_center, src->obb_center);
    VMOVE(dst->obb_extent1, src->obb_extent1);
    VMOVE(dst->obb_extent2, src->obb_extent2);
    VMOVE(dst->obb_extent3, src->obb_extent3);
    dst->radius = src->radius;
    dst->knobs = src->knobs;
    return 1;
}

int
bv_name_set(struct bv *v, const char *name)
{
    if (!bv_is_valid(v))
	return 0;

    bu_vls_sprintf(&v->name, "%s", name ? name : "");
    return 1;
}

const char *
bv_name_get(const struct bv *v)
{
    return (bv_is_valid(v) && bu_vls_strlen(&v->name)) ?
	bu_vls_cstr(&v->name) : NULL;
}

int
bv_user_data_set(struct bv *v, void *user_data)
{
    if (!bv_is_valid(v))
	return 0;

    v->user_data = user_data;
    return 1;
}

void *
bv_user_data_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->user_data : NULL;
}

int
bv_dimensions_set(struct bv *v, int width, int height)
{
    if (!bv_is_valid(v))
	return 0;

    v->width = width;
    v->height = height;
    return 1;
}

int
bv_refresh_request(struct bv *v, uint32_t flags)
{
    if (!bv_is_valid(v))
	return 0;

    if (v->refresh_enabled && !v->refresh_suppressed)
	v->refresh_dirty |= flags;
    return 1;
}

int
bv_refresh_dirty_get(const struct bv *v)
{
    return (bv_is_valid(v) && v->refresh_dirty) ? 1 : 0;
}

uint32_t
bv_refresh_consume(struct bv *v)
{
    uint32_t dirty;

    if (!bv_is_valid(v))
	return 0;

    dirty = v->refresh_dirty;
    v->refresh_dirty = 0;
    return dirty;
}

int
bv_refresh_complete(struct bv *v)
{
    if (!bv_is_valid(v))
	return 0;

    v->refresh_dirty = 0;
    return 1;
}

int
bv_refresh_enabled_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->refresh_enabled : 0;
}

int
bv_refresh_enabled_set(struct bv *v, int enabled)
{
    if (!bv_is_valid(v))
	return 0;

    v->refresh_enabled = enabled ? 1 : 0;
    return 1;
}

int
bv_refresh_suppressed_get(const struct bv *v)
{
    return (bv_is_valid(v) && v->refresh_suppressed > 0) ? 1 : 0;
}

int
bv_refresh_suppress_begin(struct bv *v)
{
    if (!bv_is_valid(v))
	return 0;

    v->refresh_suppressed++;
    return 1;
}

int
bv_refresh_suppress_end(struct bv *v)
{
    if (!bv_is_valid(v))
	return 0;

    if (v->refresh_suppressed > 0)
	v->refresh_suppressed--;
    return 1;
}

int
bv_refresh_drawn_count_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->refresh_drawn_count : 0;
}

int
bv_refresh_drawn_count_set(struct bv *v, int count)
{
    if (!bv_is_valid(v))
	return 0;

    v->refresh_drawn_count = count;
    return 1;
}

int
bv_frametime_set(struct bv *v, uint64_t frametime)
{
    if (!bv_is_valid(v))
	return 0;

    v->frametime = frametime;
    return 1;
}

uint64_t
bv_frametime_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->frametime : 0;
}

int
bv_zclip_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->zclip : 0;
}

int
bv_zclip_set(struct bv *v, int zclip)
{
    if (!bv_is_valid(v))
	return 0;

    v->zclip = zclip ? 1 : 0;
    return 1;
}

int
bv_framebuffer_mode_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->framebuffer_mode : 0;
}

int
bv_framebuffer_mode_set(struct bv *v, int mode)
{
    if (!bv_is_valid(v))
	return 0;

    v->framebuffer_mode = mode;
    return 1;
}

int
bv_cleared_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->cleared : 0;
}

int
bv_cleared_set(struct bv *v, int cleared)
{
    if (!bv_is_valid(v))
	return 0;

    v->cleared = cleared ? 1 : 0;
    return 1;
}

int
bv_update(struct bv *v)
{
    vect_t work;
    vect_t work1;
    vect_t temp;
    vect_t temp1;

    if (!bv_is_valid(v))
	return 0;

    bn_mat_mul(v->model2view, v->rotation, v->center);
    v->model2view[15] = v->scale;
    bn_mat_inv(v->view2model, v->model2view);

    VSET(work, 0.0, 0.0, 1.0);
    MAT4X3VEC(temp, v->view2model, work);
    VSET(work1, 1.0, 0.0, 0.0);
    MAT4X3VEC(temp1, v->view2model, work1);
    bn_aet_vec(&v->aet[X], &v->aet[Y], &v->aet[Z], temp, temp1,
	    (fastf_t)0.005);

    if ((NEAR_EQUAL(v->aet[Y], 90.0, (fastf_t)0.005) ||
		NEAR_EQUAL(v->aet[Y], -90.0, (fastf_t)0.005)) &&
	    v->aet[X] < 0.0 && !NEAR_ZERO(v->aet[X], (fastf_t)0.005))
	v->aet[X] += 360.0;
    else if (NEAR_ZERO(v->aet[X], (fastf_t)0.005))
	v->aet[X] = 0.0;

    bn_mat_mul(v->pmodel2view, v->pmat, v->model2view);
    v->frame_revision++;
    return 1;
}

int
bv_model2view_get(mat_t model2view, const struct bv *v)
{
    if (!model2view)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(model2view);
	return 0;
    }
    MAT_COPY(model2view, v->model2view);
    return 1;
}

int
bv_view2model_get(mat_t view2model, const struct bv *v)
{
    if (!view2model)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(view2model);
	return 0;
    }
    MAT_COPY(view2model, v->view2model);
    return 1;
}

int
bv_pmodel2view_get(mat_t pmodel2view, const struct bv *v)
{
    if (!pmodel2view)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(pmodel2view);
	return 0;
    }
    MAT_COPY(pmodel2view, v->pmodel2view);
    return 1;
}

int
bv_pmat_get(mat_t pmat, const struct bv *v)
{
    if (!pmat)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(pmat);
	return 0;
    }
    MAT_COPY(pmat, v->pmat);
    return 1;
}

int
bv_pmat_set(struct bv *v, const mat_t pmat)
{
    if (!bv_is_valid(v) || !pmat)
	return 0;

    MAT_COPY(v->pmat, pmat);
    bn_mat_mul(v->pmodel2view, v->pmat, v->model2view);
    v->frame_revision++;
    return 1;
}

int
bv_rotation_get(mat_t rotation, const struct bv *v)
{
    if (!rotation)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(rotation);
	return 0;
    }
    MAT_COPY(rotation, v->rotation);
    return 1;
}

int
bv_rotation_set(struct bv *v, const mat_t rotation)
{
    if (!bv_is_valid(v) || !rotation)
	return 0;

    MAT_COPY(v->rotation, rotation);
    return bv_update(v);
}

int
bv_center_get(point_t center, const struct bv *v)
{
    if (!center)
	return 0;
    if (!bv_is_valid(v)) {
	VSETALL(center, 0.0);
	return 0;
    }
    MAT_DELTAS_GET_NEG(center, v->center);
    return 1;
}

int
bv_center_set(struct bv *v, const point_t center)
{
    if (!bv_is_valid(v) || !center)
	return 0;

    MAT_DELTAS_VEC_NEG(v->center, center);
    return bv_update(v);
}

int
bv_scale_set(struct bv *v, fastf_t scale)
{
    if (!bv_is_valid(v))
	return 0;

    v->scale = (scale < BV_MIN_SCALE) ? BV_MIN_SCALE : scale;
    v->size = 2.0 * v->scale;
    v->inverse_size = (v->size > 0.0) ? 1.0 / v->size : 0.0;
    return bv_update(v);
}

int
bv_size_set(struct bv *v, fastf_t size)
{
    if (!bv_is_valid(v))
	return 0;

    v->size = (size < BV_MIN_SIZE) ? BV_MIN_SIZE : size;
    v->scale = 0.5 * v->size;
    if (v->scale < BV_MIN_SCALE)
	v->scale = BV_MIN_SCALE;
    v->inverse_size = (v->size > 0.0) ? 1.0 / v->size : 0.0;
    return bv_update(v);
}

int
bv_aet_get(vect_t aet, const struct bv *v)
{
    if (!aet)
	return 0;
    if (!bv_is_valid(v)) {
	VSETALL(aet, 0.0);
	return 0;
    }
    VMOVE(aet, v->aet);
    return 1;
}

int
bv_aet_set(struct bv *v, const vect_t aet)
{
    if (!bv_is_valid(v) || !aet)
	return 0;

    VMOVE(v->aet, aet);
    bv_mat_aet(v);
    return bv_update(v);
}

int
bv_aet_state_set(struct bv *v, const vect_t aet)
{
    if (!bv_is_valid(v) || !aet)
	return 0;

    VMOVE(v->aet, aet);
    v->frame_revision++;
    return 1;
}

int
bv_eye_pos_get(point_t eye_pos, const struct bv *v)
{
    if (!eye_pos)
	return 0;
    if (!bv_is_valid(v)) {
	VSET(eye_pos, 0.0, 0.0, 1.0);
	return 0;
    }
    VMOVE(eye_pos, v->eye_pos);
    return 1;
}

int
bv_eye_pos_set(struct bv *v, const point_t eye_pos)
{
    if (!bv_is_valid(v) || !eye_pos)
	return 0;

    VMOVE(v->eye_pos, eye_pos);
    v->frame_revision++;
    return 1;
}

int
bv_keypoint_get(point_t keypoint, const struct bv *v)
{
    if (!keypoint)
	return 0;
    if (!bv_is_valid(v)) {
	VSETALL(keypoint, 0.0);
	return 0;
    }
    VMOVE(keypoint, v->keypoint);
    return 1;
}

int
bv_keypoint_set(struct bv *v, const point_t keypoint)
{
    if (!bv_is_valid(v) || !keypoint)
	return 0;

    VMOVE(v->keypoint, keypoint);
    v->frame_revision++;
    return 1;
}

int
bv_autoview_bounds(struct bv *v, fastf_t scale, const point_t min, const point_t max)
{
    vect_t center = VINIT_ZERO;
    vect_t radial;
    vect_t sqrt_small;

    if (!bv_is_valid(v) || !min || !max)
	return 0;

    if (scale < SQRT_SMALL_FASTF)
	scale = 2.0;

    VSETALL(sqrt_small, SQRT_SMALL_FASTF);
    VADD2SCALE(center, max, min, 0.5);
    VSUB2(radial, max, center);
    VMAX(radial, sqrt_small);
    if (VNEAR_ZERO(radial, SQRT_SMALL_FASTF))
	VSETALL(radial, 1.0);

    MAT_IDN(v->center);
    MAT_DELTAS_VEC_NEG(v->center, center);
    v->scale = radial[X];
    V_MAX(v->scale, radial[Y]);
    V_MAX(v->scale, radial[Z]);
    v->size = scale * v->scale;
    v->inverse_size = (v->size > 0.0) ? 1.0 / v->size : 0.0;
    return bv_update(v);
}

int
bv_screen_to_view(fastf_t *vx, fastf_t *vy, const struct bv *v, fastf_t sx, fastf_t sy)
{
    fastf_t aspect;

    if (vx)
	*vx = 0.0;
    if (vy)
	*vy = 0.0;
    if (!vx || !vy || !bv_is_valid(v) || !v->width || !v->height)
	return 0;

    *vx = sx / (fastf_t)v->width * 2.0 - 1.0;
    aspect = (fastf_t)v->width / (fastf_t)v->height;
    *vy = (sy / (fastf_t)v->height * -2.0 + 1.0) / aspect;
    return 1;
}

int
bv_screen_to_model(point_t model_point, const struct bv *v, fastf_t sx, fastf_t sy)
{
    fastf_t vx = 0.0;
    fastf_t vy = 0.0;
    point_t vpt = VINIT_ZERO;

    if (model_point)
	VSETALL(model_point, 0.0);
    if (!model_point || !bv_screen_to_view(&vx, &vy, v, sx, sy))
	return 0;

    VSET(vpt, vx, vy, 0.0);
    MAT4X3PNT(model_point, v->view2model, vpt);
    return 1;
}

int
bv_plane_get(plane_t *plane, const struct bv *v)
{
    point_t center;
    vect_t normal;

    if (!plane || !bv_is_valid(v))
	return -1;

    MAT_DELTAS_GET_NEG(center, v->center);
    VMOVEN(normal, v->rotation + 8, 3);
    VUNITIZE(normal);
    VSCALE(normal, normal, -1.0);
    return bg_plane_pt_nrml(plane, center, normal);
}

int
bv_mouse_state_set(struct bv *v, int x, int y)
{
    if (!bv_is_valid(v))
	return 0;

    v->previous_mouse_x = (fastf_t)v->mouse_x;
    v->previous_mouse_y = (fastf_t)v->mouse_y;
    v->mouse_x = x;
    v->mouse_y = y;
    return bv_screen_to_model(v->current_point, v, (fastf_t)x, (fastf_t)y);
}

int
bv_adjust(struct bv *v, int dx, int dy, const point_t keypoint, int UNUSED(mode), unsigned long long flags)
{
    if (!bv_is_valid(v) || flags == BV_ADJUST_IDLE)
	return 0;

    if (flags & BV_ADJUST_ROT) {
	point_t rot_pt;
	point_t new_origin;
	mat_t viewchg;
	mat_t viewchginv;
	point_t new_cent_view;
	point_t new_cent_model;
	mat_t newrot;

	if (!keypoint)
	    return 0;
	bn_mat_angles(newrot, (fastf_t)dy * 0.25, (fastf_t)dx * 0.25, 0.0);
	MAT4X3PNT(rot_pt, v->model2view, keypoint);
	bn_mat_xform_about_pnt(viewchg, newrot, rot_pt);
	bn_mat_inv(viewchginv, viewchg);
	VSET(new_origin, 0.0, 0.0, 0.0);
	MAT4X3PNT(new_cent_view, viewchginv, new_origin);
	MAT4X3PNT(new_cent_model, v->view2model, new_cent_view);
	MAT_DELTAS_VEC_NEG(v->center, new_cent_model);
	bn_mat_mul2(newrot, v->rotation);
	return bv_update(v);
    }

    if (flags & BV_ADJUST_TRANS) {
	fastf_t aspect;
	fastf_t fx;
	fastf_t fy;
	vect_t tt;
	point_t delta;
	point_t work;
	point_t vc;
	point_t nvc;

	if (!v->width || !v->height)
	    return 0;
	aspect = (fastf_t)v->width / (fastf_t)v->height;
	fx = (fastf_t)dx / (fastf_t)v->width * 2.0;
	fy = -(fastf_t)dy / (fastf_t)v->height / aspect * 2.0;
	VSET(tt, fx, fy, 0.0);
	MAT4X3PNT(work, v->view2model, tt);
	MAT_DELTAS_GET_NEG(vc, v->center);
	VSUB2(delta, work, vc);
	VSUB2(nvc, vc, delta);
	MAT_DELTAS_VEC_NEG(v->center, nvc);
	return bv_update(v);
    }

    if (flags & BV_ADJUST_SCALE) {
	double f;

	if (!dx || !dy)
	    return 0;
	f = (double)dy / (double)dx;
	if (NEAR_ZERO(f, SQRT_SMALL_FASTF))
	    return 0;
	v->scale /= f;
	if (v->scale < BV_MIN_SCALE)
	    v->scale = BV_MIN_SCALE;
	v->size = 2.0 * v->scale;
	v->inverse_size = (v->size > 0.0) ? 1.0 / v->size : 0.0;
	return bv_update(v);
    }

    if (flags & BV_ADJUST_CENTER) {
	point_t center;
	if (!bv_screen_to_model(center, v, (fastf_t)dx, (fastf_t)dy))
	    return 0;
	MAT_DELTAS_VEC_NEG(v->center, center);
	return bv_update(v);
    }

    return 0;
}

int
bv_knobs_reset(struct bv_knobs *knobs, int category)
{
    if (!knobs)
	return 0;

    if (!category || category == BV_KNOBS_RATE) {
	knobs->rot_model_active = 0;
	VSETALL(knobs->rot_model, 0.0);
	knobs->rot_object_active = 0;
	VSETALL(knobs->rot_object, 0.0);
	knobs->rot_view_active = 0;
	VSETALL(knobs->rot_view, 0.0);
	knobs->trans_model_active = 0;
	VSETALL(knobs->trans_model, 0.0);
	knobs->trans_view_active = 0;
	VSETALL(knobs->trans_view, 0.0);
	knobs->scale_rate_active = 0;
	knobs->scale_rate = 0.0;
    }

    if (!category || category == BV_KNOBS_ABS) {
	VSETALL(knobs->abs_rot_model, 0.0);
	VSETALL(knobs->abs_rot_model_last, 0.0);
	VSETALL(knobs->abs_rot_object, 0.0);
	VSETALL(knobs->abs_rot_object_last, 0.0);
	VSETALL(knobs->abs_rot_view, 0.0);
	VSETALL(knobs->abs_rot_view_last, 0.0);
	VSETALL(knobs->abs_trans_model, 0.0);
	VSETALL(knobs->abs_trans_model_last, 0.0);
	VSETALL(knobs->abs_trans_view, 0.0);
	VSETALL(knobs->abs_trans_view_last, 0.0);
	knobs->abs_scale = 0.0;
    }

    return 1;
}

unsigned long long
bv_knobs_hash(const struct bv_knobs *knobs, struct bu_data_hash_state *state)
{
    int own_state = 0;
    unsigned long long hv;

    if (!knobs)
	return 0ULL;

    if (!state) {
	state = bu_data_hash_create();
	if (!state)
	    return 0ULL;
	own_state = 1;
    }

    bu_data_hash_update(state, &knobs->rot_model, sizeof(knobs->rot_model));
    bu_data_hash_update(state, &knobs->rot_model_active, sizeof(knobs->rot_model_active));
    bu_data_hash_update(state, &knobs->rot_model_origin, sizeof(knobs->rot_model_origin));
    bu_data_hash_update(state, &knobs->rot_object, sizeof(knobs->rot_object));
    bu_data_hash_update(state, &knobs->rot_object_active, sizeof(knobs->rot_object_active));
    bu_data_hash_update(state, &knobs->rot_object_origin, sizeof(knobs->rot_object_origin));
    bu_data_hash_update(state, &knobs->rot_view, sizeof(knobs->rot_view));
    bu_data_hash_update(state, &knobs->rot_view_active, sizeof(knobs->rot_view_active));
    bu_data_hash_update(state, &knobs->rot_view_origin, sizeof(knobs->rot_view_origin));
    bu_data_hash_update(state, &knobs->scale_rate, sizeof(knobs->scale_rate));
    bu_data_hash_update(state, &knobs->scale_rate_active, sizeof(knobs->scale_rate_active));
    bu_data_hash_update(state, &knobs->trans_model, sizeof(knobs->trans_model));
    bu_data_hash_update(state, &knobs->trans_model_active, sizeof(knobs->trans_model_active));
    bu_data_hash_update(state, &knobs->trans_view, sizeof(knobs->trans_view));
    bu_data_hash_update(state, &knobs->trans_view_active, sizeof(knobs->trans_view_active));
    bu_data_hash_update(state, &knobs->abs_rot_model, sizeof(knobs->abs_rot_model));
    bu_data_hash_update(state, &knobs->abs_rot_model_last, sizeof(knobs->abs_rot_model_last));
    bu_data_hash_update(state, &knobs->abs_rot_object, sizeof(knobs->abs_rot_object));
    bu_data_hash_update(state, &knobs->abs_rot_object_last, sizeof(knobs->abs_rot_object_last));
    bu_data_hash_update(state, &knobs->abs_rot_view, sizeof(knobs->abs_rot_view));
    bu_data_hash_update(state, &knobs->abs_rot_view_last, sizeof(knobs->abs_rot_view_last));
    bu_data_hash_update(state, &knobs->abs_scale, sizeof(knobs->abs_scale));
    bu_data_hash_update(state, &knobs->abs_trans_model, sizeof(knobs->abs_trans_model));
    bu_data_hash_update(state, &knobs->abs_trans_model_last, sizeof(knobs->abs_trans_model_last));
    bu_data_hash_update(state, &knobs->abs_trans_view, sizeof(knobs->abs_trans_view));
    bu_data_hash_update(state, &knobs->abs_trans_view_last, sizeof(knobs->abs_trans_view_last));

    if (!own_state)
	return 0ULL;

    hv = bu_data_hash_val(state);
    bu_data_hash_destroy(state);
    return hv;
}

int
bv_knobs_cmd_process(vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, struct bv *v, const char *cmd, fastf_t factor, char origin, int model_flag, int incr_flag)
{
    struct bv_knobs *k;
    int ind = -1;
    char c;

    if (!bv_is_valid(v) || !cmd || !rvec || !do_rot || !tvec || !do_tran)
	return BRLCAD_ERROR;

    k = &v->knobs;
    c = (cmd[1] == '\0') ? cmd[0] : cmd[1];
    switch (c) {
	case 'x':
	case 'X':
	    ind = X;
	    break;
	case 'y':
	case 'Y':
	    ind = Y;
	    break;
	case 'z':
	case 'Z':
	    ind = Z;
	    break;
	case 'S':
	    ind = 0;
	    break;
	default:
	    return BRLCAD_ERROR;
    }

    if (cmd[1] == '\0') {
	if (cmd[0] == 'x' || cmd[0] == 'y' || cmd[0] == 'z') {
	    fastf_t *coord = model_flag ? &k->rot_model[ind] : &k->rot_view[ind];
	    if (incr_flag)
		*coord += factor;
	    else
		*coord = factor;
	    if (model_flag)
		k->rot_model_origin = origin;
	    else
		k->rot_view_origin = origin;
	    return BRLCAD_OK;
	}
	if (cmd[0] == 'X' || cmd[0] == 'Y' || cmd[0] == 'Z') {
	    fastf_t *coord = model_flag ? &k->trans_model[ind] : &k->trans_view[ind];
	    if (incr_flag)
		*coord += factor;
	    else
		*coord = factor;
	    return BRLCAD_OK;
	}
	if (cmd[0] == 'S') {
	    if (incr_flag)
		k->scale_rate += factor;
	    else
		k->scale_rate = factor;
	    return BRLCAD_OK;
	}
	return BRLCAD_ERROR;
    }

    if (cmd[0] == 'a' && strlen(cmd) == 2) {
	if (cmd[1] == 'x' || cmd[1] == 'y' || cmd[1] == 'z') {
	    fastf_t *abs_rot = model_flag ? k->abs_rot_model : k->abs_rot_view;
	    fastf_t *last = model_flag ? k->abs_rot_model_last : k->abs_rot_view_last;
	    if (incr_flag) {
		abs_rot[ind] += factor;
		(*rvec)[ind] = factor;
	    } else {
		(*rvec)[ind] = factor - last[ind];
		abs_rot[ind] = factor;
	    }
	    if (abs_rot[ind] < -180.0)
		abs_rot[ind] += 360.0;
	    else if (abs_rot[ind] > 180.0)
		abs_rot[ind] -= 360.0;
	    last[ind] = abs_rot[ind];
	    *do_rot = 1;
	    return BRLCAD_OK;
	}
	if (cmd[1] == 'X' || cmd[1] == 'Y' || cmd[1] == 'Z') {
	    fastf_t *abs_tran = model_flag ? k->abs_trans_model : k->abs_trans_view;
	    fastf_t *last = model_flag ? k->abs_trans_model_last : k->abs_trans_view_last;
	    fastf_t scale = (v->scale > SMALL_FASTF) ? v->scale : 1.0;
	    fastf_t sf = factor * v->local2base / scale;
	    if (incr_flag) {
		abs_tran[ind] += sf;
		last[ind] = abs_tran[ind];
		(*tvec)[ind] = factor;
	    } else {
		(*tvec)[ind] = factor - last[ind] * v->scale * v->base2local;
		abs_tran[ind] = sf;
		last[ind] = abs_tran[ind];
	    }
	    *do_tran = 1;
	    return BRLCAD_OK;
	}
	if (cmd[1] == 'S') {
	    if (incr_flag)
		v->absolute_scale += factor;
	    else
		v->absolute_scale = factor;
	    if (-SMALL_FASTF < v->absolute_scale && v->absolute_scale < SMALL_FASTF)
		v->scale = v->initial_scale;
	    else if (v->absolute_scale > 0.0)
		v->scale = v->initial_scale * (1.0 - v->absolute_scale);
	    else
		v->scale = v->initial_scale * (1.0 + (v->absolute_scale * -9.0));
	    if (v->scale < BV_MIN_SCALE)
		v->scale = BV_MIN_SCALE;
	    v->size = 2.0 * v->scale;
	    v->inverse_size = (v->size > 0.0) ? 1.0 / v->size : 0.0;
	    return bv_update(v) ? BRLCAD_OK : BRLCAD_ERROR;
	}
    }

    return BRLCAD_ERROR;
}

int
bv_knobs_translate(struct bv *v, const vect_t tvec, int model_flag)
{
    point_t vc;
    point_t nvc;

    if (!bv_is_valid(v) || !tvec)
	return 0;

    if (model_flag || v->coord_mode == 'o') {
	point_t delta;
	VSCALE(delta, tvec, v->local2base);
	MAT_DELTAS_GET_NEG(vc, v->center);
	VSUB2(nvc, vc, delta);
    } else {
	vect_t tt;
	point_t delta;
	point_t work;
	fastf_t scale = (v->scale > SMALL_FASTF) ? v->scale : 1.0;
	VSCALE(tt, tvec, v->local2base / scale);
	MAT4X3PNT(work, v->view2model, tt);
	MAT_DELTAS_GET_NEG(vc, v->center);
	VSUB2(delta, work, vc);
	VSUB2(nvc, vc, delta);
    }
    MAT_DELTAS_VEC_NEG(v->center, nvc);
    return bv_update(v);
}

int
bv_knobs_rotate(struct bv *v, const vect_t rvec, char origin, char coords, const matp_t obj_rot, const pointp_t pvt_pt)
{
    mat_t base;
    mat_t view_rot;
    point_t pivot_view;
    int recenter = 0;

    if (!bv_is_valid(v) || !rvec)
	return 0;

    MAT_IDN(base);
    bn_mat_angles(base, rvec[X], rvec[Y], rvec[Z]);
    MAT_IDN(view_rot);

    if (coords == 'm') {
	mat_t t1;
	mat_t rvinv;
	bn_mat_inv(rvinv, v->rotation);
	bn_mat_mul(t1, v->rotation, base);
	bn_mat_mul(view_rot, t1, rvinv);
    } else if (coords == 'o' && obj_rot) {
	mat_t obj_inv;
	mat_t model_rot;
	mat_t t1;
	mat_t t2;
	mat_t rvinv;
	bn_mat_inv(obj_inv, obj_rot);
	bn_mat_mul(t2, obj_rot, base);
	bn_mat_mul(model_rot, t2, obj_inv);
	bn_mat_inv(rvinv, v->rotation);
	bn_mat_mul(t1, v->rotation, model_rot);
	bn_mat_mul(view_rot, t1, rvinv);
    } else {
	MAT_COPY(view_rot, base);
    }

    switch (origin) {
	case 'e':
	    VSET(pivot_view, 0.0, 0.0, 1.0);
	    recenter = 1;
	    break;
	case 'm': {
	    point_t mzero = VINIT_ZERO;
	    MAT4X3PNT(pivot_view, v->model2view, mzero);
	    recenter = 1;
	    break;
	}
	case 'k': {
	    point_t mp;
	    if (pvt_pt)
		VMOVE(mp, pvt_pt);
	    else
		VSETALL(mp, 0.0);
	    MAT4X3PNT(pivot_view, v->model2view, mp);
	    recenter = 1;
	    break;
	}
	case 'v':
	default:
	    VSETALL(pivot_view, 0.0);
	    break;
    }

    if (recenter) {
	mat_t about;
	mat_t about_inv;
	point_t new_origin = VINIT_ZERO;
	point_t new_cent_view;
	point_t new_cent_model;
	bn_mat_xform_about_pnt(about, view_rot, pivot_view);
	bn_mat_inv(about_inv, about);
	MAT4X3PNT(new_cent_view, about_inv, new_origin);
	MAT4X3PNT(new_cent_model, v->view2model, new_cent_view);
	MAT_DELTAS_VEC_NEG(v->center, new_cent_model);
    }

    bn_mat_mul2(view_rot, v->rotation);
    return bv_update(v);
}

int
bv_knobs_update_rate_flags(struct bv *v)
{
    struct bv_knobs *k;

    if (!bv_is_valid(v))
	return 0;

    k = &v->knobs;
    k->rot_model_active = (!ZERO(k->rot_model[X]) || !ZERO(k->rot_model[Y]) || !ZERO(k->rot_model[Z]));
    k->trans_model_active = (!ZERO(k->trans_model[X]) || !ZERO(k->trans_model[Y]) || !ZERO(k->trans_model[Z]));
    k->rot_view_active = (!ZERO(k->rot_view[X]) || !ZERO(k->rot_view[Y]) || !ZERO(k->rot_view[Z]));
    k->trans_view_active = (!ZERO(k->trans_view[X]) || !ZERO(k->trans_view[Y]) || !ZERO(k->trans_view[Z]));
    k->scale_rate_active = !ZERO(k->scale_rate);
    return 1;
}

int
bv_knob_values_get(struct bv_knob_values *values, const struct bv *v)
{
    if (!values)
	return 0;

    memset(values, 0, sizeof(*values));
    if (!bv_is_valid(v))
	return 0;

    VMOVE(values->rate_rotation, v->knobs.rot_view);
    VMOVE(values->rate_translation, v->knobs.trans_view);
    values->rate_scale = v->knobs.scale_rate;
    VMOVE(values->absolute_rotation, v->knobs.abs_rot_view);
    VMOVE(values->absolute_translation, v->knobs.abs_trans_view);
    values->absolute_scale = v->absolute_scale;
    return 1;
}

int
bv_knobs_calibrate(struct bv *v)
{
    if (!bv_is_valid(v))
	return 0;

    VSETALL(v->knobs.abs_trans_view, 0.0);
    VSETALL(v->knobs.abs_trans_view_last, 0.0);
    VSETALL(v->knobs.abs_trans_model, 0.0);
    VSETALL(v->knobs.abs_trans_model_last, 0.0);
    return 1;
}

struct bv_set *
bv_set_create(void)
{
    struct bv_set *s = NULL;
    BU_ALLOC(s, struct bv_set);
    bv_set_init(s);
    return s;
}

void
bv_set_destroy(struct bv_set *s)
{
    if (!s)
	return;

    bv_set_free(s);
    bu_free(s, "bv_set");
}

void
bv_set_init(struct bv_set *s)
{
    if (!s)
	return;

    memset(s, 0, sizeof(*s));
    s->magic = BV_SET_MAGIC;
    BU_PTBL_INIT(&s->views);
}

void
bv_set_free(struct bv_set *s)
{
    if (!s || s->magic != BV_SET_MAGIC)
	return;

    bu_ptbl_free(&s->views);
    s->magic = 0;
}

int
bv_set_add(struct bv_set *s, struct bv *v)
{
    if (!s || s->magic != BV_SET_MAGIC || !bv_is_valid(v))
	return 0;

    if (v->view_set && v->view_set != s)
	bv_set_remove(v->view_set, v);

    bu_ptbl_ins_unique(&s->views, (long *)v);
    v->view_set = s;
    return 1;
}

int
bv_set_remove(struct bv_set *s, struct bv *v)
{
    if (!s || s->magic != BV_SET_MAGIC || !v)
	return 0;

    if (v->view_set == s)
	v->view_set = NULL;
    (void)bu_ptbl_rm(&s->views, (long *)v);
    return 1;
}

struct bu_ptbl *
bv_set_views(struct bv_set *s)
{
    return (s && s->magic == BV_SET_MAGIC) ? &s->views : NULL;
}

struct bv *
bv_set_find(struct bv_set *s, const char *name)
{
    if (!s || s->magic != BV_SET_MAGIC || !name)
	return NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&s->views); i++) {
	struct bv *v = (struct bv *)BU_PTBL_GET(&s->views, i);
	const char *vname = bv_name_get(v);
	if (vname && BU_STR_EQUAL(vname, name))
	    return v;
    }

    return NULL;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
