/*                     T E S T _ G Q A . C
 * BRL-CAD
 *
 * Copyright (c) 2018-2026 United States Government as represented by
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
/** @file test_gqa.c
 *
 * Extract plotting data from a gqa run
 *
 */

#include "common.h"

#include <stdio.h>
#include <string.h>
#include <bu.h>
#include "bg/plot3.h"
#include <ged.h>
#include <ged/bsg_ged_draw.h>

struct gqa_segment_writer {
    FILE *fp;
    size_t count;
};

struct gqa_label_capture {
    size_t count;
    char text[256];
};

struct gqa_record_capture {
    struct gqa_segment_writer *writer;
    const char *plot_fname;
};

static int
capture_gqa_label(struct ged *UNUSED(gedp),
	const struct ged_diagnostic_hud_label *label,
	void *data)
{
    struct gqa_label_capture *capture = (struct gqa_label_capture *)data;
    if (!capture || !label || !label->label_id || !label->text)
	return 0;
    if (!strstr(label->label_id, "gqa::overlaps::summary"))
	return 0;
    capture->count++;
    bu_strlcpy(capture->text, label->text, sizeof(capture->text));
    return 1;
}

static int
write_gqa_segment(const point_t a, const point_t b, void *data)
{
    struct gqa_segment_writer *writer = (struct gqa_segment_writer *)data;
    pdv_3line(writer->fp, a, b);
    writer->count++;
    return 1;
}

static int
capture_gqa_record(const struct ged_draw_view_db_object_record *rec, void *data)
{
    struct gqa_record_capture *capture = (struct gqa_record_capture *)data;
    struct gqa_segment_writer *writer = capture ? capture->writer : NULL;

    if (!writer || !rec || !rec->path || !strstr(rec->path, "gqa::overlaps"))
	return 1;
    if (!rec->vlist_point_count)
	return 1;

    if (!writer->fp) {
	writer->fp = fopen(capture->plot_fname, "wb");
	if (!writer->fp)
	    bu_exit(EXIT_FAILURE, "Could not open %s for writing\n",
		    capture->plot_fname);
	pl_color(writer->fp, 255, 255, 0);
    }
    printf("found %s;\n", rec->path);
    (void)ged_draw_view_db_object_record_foreach_segment(rec,
	    write_gqa_segment, writer);

    return 1;
}

int
main(int ac, char *av[]) {
    struct ged *gedp;
    const char *gqa_plot_fname = "gqa_ovlps.plot3";
    const char *gqa[4] = {"gqa", "-Aop", "ovlp", NULL};
    struct gqa_segment_writer writer = {NULL, 0};
    struct gqa_label_capture label_capture = {0, {0}};

    bu_setprogname(av[0]);

    if (ac != 2) {
	printf("Usage: %s file.g\n", av[0]);
	return 1;
    }
    if (!bu_file_exists(av[1], NULL)) {
	printf("ERROR: [%s] does not exist, expecting .g file\n", av[1]);
	return 2;
    }

    gedp = ged_open("db", av[1], 1);
    ged_diagnostic_hud_label_handler_set(gedp, capture_gqa_label,
	    &label_capture);
    ged_exec_gqa(gedp, 3, gqa);
    ged_diagnostic_hud_label_handler_set(gedp, NULL, NULL);
    printf("%s\n", bu_vls_cstr(gedp->ged_result_str));
    if (!label_capture.count ||
	    !strstr(label_capture.text, "gqa overlaps:"))
	bu_exit(EXIT_FAILURE, "No GQA overlap summary HUD label published.\n");

    void *view_ctx = ged_view_active_ctx(gedp);
    if (!view_ctx)
	bu_exit(EXIT_FAILURE, "No active GED view context available.\n");

    struct ged_draw_view_record_query query;
    memset(&query, 0, sizeof(query));
    query.flags = GED_DRAW_VIEW_RECORD_QUERY_VIEW_OBJECTS;
    query.draw_mode = -1;
    struct gqa_record_capture record_capture = {&writer, gqa_plot_fname};
    ged_draw_foreach_view_record_query(view_ctx, &query, capture_gqa_record,
	    &record_capture);

    if (writer.count) {
	printf("Writing plot data to %s for inspection with overlay command\n", gqa_plot_fname);
	fclose(writer.fp);
    } else {
	if (writer.fp)
	    fclose(writer.fp);
	bu_exit(EXIT_FAILURE, "No GQA plotting data found.\n");
    }

    ged_close(gedp);

    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
