/*                   E D I T _ T E S T _ V I E W . H
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
/** @file edit_test_view.h
 *
 * Test helpers for constructing RT-owned edit view snapshots.
 */

#ifndef LIBRT_TESTS_EDIT_TEST_VIEW_H
#define LIBRT_TESTS_EDIT_TEST_VIEW_H

#include <math.h>

#include "vmath.h"
#include "bn/mat.h"
#include "rt/edit.h"
#include "rt/view.h"

static inline void
rt_edit_test_view_update(struct rt_edit_view *v)
{
    if (!v)
	return;

    bn_mat_mul(v->gv_model2view, v->gv_rotation, v->gv_center);
    v->gv_model2view[15] = v->gv_scale;
    bn_mat_inv(v->gv_view2model, v->gv_model2view);
}

static inline void
rt_edit_test_view_set_aet(struct rt_edit_view *v,
	fastf_t az, fastf_t el, fastf_t twist, fastf_t size)
{
    mat_t tmat;
    fastf_t twist_rad;

    if (!v)
	return;

    bn_mat_angles(v->gv_rotation, 270.0 + el, 0.0, 270.0 - az);

    twist_rad = -twist * DEG2RAD;
    bn_mat_zrot(tmat, sin(twist_rad), cos(twist_rad));
    bn_mat_mul2(tmat, v->gv_rotation);

    v->gv_scale = 0.5 * size;
    rt_edit_test_view_update(v);
}

static inline void
rt_edit_test_view_set_size(struct rt_edit_view *v, fastf_t size)
{
    if (!v)
	return;

    v->gv_scale = 0.5 * size;
    rt_edit_test_view_update(v);
}

static inline void
rt_edit_test_view_init_size(struct rt_edit_view *v, fastf_t size)
{
    if (!v)
	return;

    rt_edit_view_init(v);
    rt_edit_test_view_set_aet(v, 45.0, 35.0, 0.0, size);
}

static inline void
rt_edit_test_view_init_identity_size(struct rt_edit_view *v, fastf_t size)
{
    if (!v)
	return;

    rt_edit_view_init(v);
    rt_edit_test_view_set_size(v, size);
}

static inline void
rt_edit_test_view_init(struct rt_edit_view *v)
{
    if (!v)
	return;

    rt_edit_test_view_init_size(v, 73.3197);
}

#endif /* LIBRT_TESTS_EDIT_TEST_VIEW_H */
