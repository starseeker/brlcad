/*                       D R A W . C
 * BRL-CAD
 *
 * Copyright (c) 2000-2026 United States Government as represented by
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
/** @addtogroup libtclcad */
/** @{ */
/** @file libtclcad/view/draw.c
 *
 */
/** @} */

#include "common.h"
#include "dm/view.h"
#include "ged.h"
#include "ged/bsg_ged_draw.h"
#include "rt/view_legacy_bsg.h"
#include "tclcad.h"

/* Private headers */
#include "../tclcad_private.h"

#include "../view/view.h"



void
go_draw(struct bsg_view *gdvp)
{
    struct dm *dmp = (struct dm *)rt_view_display_manager_from_bsg(gdvp);
    if (!dmp)
	return;

    mat_t model2view;
    mat_t pmat;
    fastf_t perspective = rt_view_perspective_from_bsg(gdvp);

    rt_view_model2view_from_bsg(model2view, gdvp);
    rt_view_pmat_from_bsg(pmat, gdvp);
    (void)dm_loadmatrix(dmp, model2view, 0);

    if (SMALL_FASTF < perspective)
	(void)dm_loadpmatrix(dmp, pmat);
    else
	(void)dm_loadpmatrix(dmp, (fastf_t *)NULL);

    dm_draw_objs(gdvp);
}


struct to_edit_redraw_data {
    struct ged *gedp;
    const char **argv;
    const struct db_full_path *subpath;
};


static int
to_edit_redraw_record_cb(const struct ged_draw_view_db_object_record *rec,
			 void *userdata)
{
    struct to_edit_redraw_data *d = (struct to_edit_redraw_data *)userdata;
    if (!rec || !rec->path || !d || !d->gedp || !d->subpath)
	return 1;

    struct db_full_path rec_path;
    db_full_path_init(&rec_path);
    if (db_string_to_path(&rec_path, d->gedp->dbip, rec->path) < 0) {
	db_free_full_path(&rec_path);
	return 1;
    }

    size_t pi;
    for (pi = 0; pi < d->subpath->fp_len; pi++) {
	if (!db_full_path_search(&rec_path, d->subpath->fp_names[pi]))
	    continue;

	/* Match found -- re-execute draw for this path */
	struct bu_vls mflag = BU_VLS_INIT_ZERO;
	struct bu_vls xflag = BU_VLS_INIT_ZERO;
	char *av[5] = {0};
	int arg = 0;

	av[arg++] = (char *)d->argv[0];
	if (rec->draw_mode == 4) {
	    av[arg++] = "-h";
	} else {
	    bu_vls_printf(&mflag, "-m%d", rec->draw_mode);
	    bu_vls_printf(&xflag, "-x%f", rec->transparency);
	    av[arg++] = bu_vls_addr(&mflag);
	    av[arg++] = bu_vls_addr(&xflag);
	}
	av[arg] = bu_strdup(rec->path);

	ged_exec(d->gedp, arg + 1, (const char **)av);

	bu_free(av[arg], "to_edit_redraw");
	bu_vls_free(&mflag);
	bu_vls_free(&xflag);
	break;
    }

    db_free_full_path(&rec_path);
    return 1;
}


int
to_edit_redraw(struct ged *gedp,
	       int argc,
	       const char *argv[])
{
    if (argc != 2)
	return BRLCAD_ERROR;

    struct db_full_path subpath;
    if (db_string_to_path(&subpath, gedp->dbip, argv[1]) != 0)
	return BRLCAD_OK;  /* path not found — nothing to do */

    struct to_edit_redraw_data d = {gedp, argv, &subpath};

    /* Iterate view-scoped DB object export records through the GED draw
     * boundary instead of walking local TclCAD draw state. */
    struct bu_ptbl *views = ged_view_set_views_ctx(gedp);
    size_t vi;
    for (vi = 0; vi < BU_PTBL_LEN(views); vi++) {
	struct bsg_view *v = (struct bsg_view *)BU_PTBL_GET(views, vi);
	ged_draw_foreach_view_db_object_record(v, to_edit_redraw_record_cb, &d);
    }

    db_free_full_path(&subpath);
    to_refresh_all_views(current_top);
    return BRLCAD_OK;
}

int
to_redraw(struct ged *gedp,
	  int argc,
	  const char *argv[],
	  ged_func_ptr UNUSED(func),
	  const char *usage,
	  int UNUSED(maxargs))
{
    /* initialize result */
    bu_vls_trunc(gedp->ged_result_str, 0);

    /* must be wanting help */
    if (argc == 1) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return GED_HELP;
    }

    if (argc != 2) {
	bu_vls_printf(gedp->ged_result_str, "Usage: %s %s", argv[0], usage);
	return BRLCAD_ERROR;
    }

    return to_edit_redraw(gedp, argc, argv);
}

int
to_blast(struct ged *gedp,
	 int argc,
	 const char *argv[],
	 ged_func_ptr UNUSED(func),
	 const char *UNUSED(usage),
	 int UNUSED(maxargs))
{
    int ret;

    ret = ged_exec(gedp, argc, argv);

    if (ret != BRLCAD_OK)
	return ret;

    to_autoview_all_views(current_top);

    return ret;
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
