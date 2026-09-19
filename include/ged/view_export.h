/*                         V I E W _ E X P O R T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged/view_export.h */

#ifndef GED_VIEW_EXPORT_H
#define GED_VIEW_EXPORT_H

#include "ged/defines.h"
#include "ged/view_export_types.h"

__BEGIN_DECLS

/** Summarize visible database presentation without exposing render records. */
GED_EXPORT extern int
ged_view_database_export_summary_get(
    struct ged_view_context *view_ctx,
    struct ged_view_database_export_summary *summary);

/** Visit model-space line segments of all visible database presentation. */
GED_EXPORT extern size_t
ged_view_database_export_segments_visit(
    struct ged_view_context *view_ctx,
    ged_view_export_segment_func_t callback,
    void *client_data);

__END_DECLS

#endif /* GED_VIEW_EXPORT_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
