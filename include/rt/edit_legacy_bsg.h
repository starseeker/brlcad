/*                    E D I T _ L E G A C Y _ B S G . H
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
/** @file edit_legacy_bsg.h
 *
 * Transitional adapters for callers that still own BSG view state while the
 * librt edit API migrates to rt_edit_view.
 */

#ifndef RT_EDIT_LEGACY_BSG_H
#define RT_EDIT_LEGACY_BSG_H

#include "rt/edit.h"

__BEGIN_DECLS

RT_EXPORT extern void
rt_edit_view_from_bsg(struct rt_edit_view *ev, const void *view_ctx);

RT_EXPORT extern struct rt_edit *
rt_edit_create_bsg(struct db_full_path *dfp, struct db_i *dbip,
		   struct bn_tol *, const void *view_ctx);

RT_EXPORT extern int
rt_edit_reinit_bsg(struct rt_edit *s, struct db_full_path *dfp,
		   struct db_i *dbip, struct bn_tol *tol,
		   const void *view_ctx);

RT_EXPORT extern int
rt_edit_knob_cmd_process_bsg(
	struct rt_edit *s,
	vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, int *do_sca,
	const void *view_ctx, const char *cmd, fastf_t f,
	char origin, int incr_flag, void *u_data
	);

__END_DECLS

#endif /* RT_EDIT_LEGACY_BSG_H */
