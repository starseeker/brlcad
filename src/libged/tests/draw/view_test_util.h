/*                 V I E W _ T E S T _ U T I L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file view_test_util.h
 *
 * Local helpers for draw tests that need passive view mechanics.  Tests still
 * get contexts through GED so owner/host policy is exercised, but passive view
 * reads and writes go directly through libbv rather than the transitional RT
 * view-context facade.
 */

#ifndef GED_TEST_DRAW_VIEW_TEST_UTIL_H
#define GED_TEST_DRAW_VIEW_TEST_UTIL_H

#include "bv.h"
#include "ged/view_types.h"

struct ged;

__BEGIN_DECLS

int draw_test_obol_view_init(struct ged *gedp, struct ged_view_context *view_ctx,
	int width, int height);

/* Advance the endpoint controller until no progressive work remains. */
int draw_test_obol_progressive_drain(struct ged *gedp, struct ged_view_context *view_ctx,
	unsigned int max_attempts, unsigned int sleep_milliseconds);

/* Emit the retained Obol scene diagnostics when GED_TEST_DRAW_OBOL_DEBUG is
 * enabled. */
void draw_test_obol_debug_scene(struct ged *gedp, int id,
	struct ged_view_context *view_ctx);

__END_DECLS

static inline struct bv *
draw_test_bv(struct ged_view_context *view_ctx)
{
    return bv_context_view(ged_view_context_bv(view_ctx));
}

static inline const struct bv *
draw_test_bv_const(const struct ged_view_context *view_ctx)
{
    return bv_context_view_const(ged_view_context_bv_const(view_ctx));
}

#define DRAW_TEST_BV(_ctx) draw_test_bv((_ctx))
#define DRAW_TEST_BV_CONST(_ctx) draw_test_bv_const((_ctx))

#endif /* GED_TEST_DRAW_VIEW_TEST_UTIL_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
