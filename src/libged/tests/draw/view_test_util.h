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

static inline struct bv *
draw_test_bv(void *view_ctx)
{
    return bv_context_view((struct bv_context *)view_ctx);
}

static inline const struct bv *
draw_test_bv_const(const void *view_ctx)
{
    return bv_context_view_const((const struct bv_context *)view_ctx);
}

#define DRAW_TEST_BV(_ctx) draw_test_bv((void *)(_ctx))
#define DRAW_TEST_BV_CONST(_ctx) draw_test_bv_const((const void *)(_ctx))

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
