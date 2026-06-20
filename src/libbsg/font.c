/*                          F O N T . C
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

#include "common.h"

#include "vmath.h"
#include "bg/vlist.h"
#include "bsg/vlist.h"

void
bsg_vlist_3string(struct bu_list *vhead,
		 struct bu_list *free_hd, /* source of free vlists */
		 const char *string,    /* string of chars to be plotted */
		 const vect_t origin,	/* lower left corner of 1st char */
		 const mat_t rot,	/* Transform matrix (WARNING: may xlate) */
		 double scale)    	/* scale factor to change 1x1 char sz */
{
    bg_vlist_3string(vhead, free_hd, string, origin, rot, scale);
}


void
bsg_vlist_2string(struct bu_list *vhead,
		 struct bu_list *free_hd,
		 const char *string,
		 double x,
		 double y,
		 double scale,
		 double theta)
{
    bg_vlist_2string(vhead, free_hd, string, x, y, scale, theta);
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
