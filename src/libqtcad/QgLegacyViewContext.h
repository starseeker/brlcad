/*              Q G L E G A C Y V I E W C O N T E X T . H
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
/** @file QgLegacyViewContext.h
 *
 * Private adapter between qtcad's opaque legacy view handle and opaque
 * retained view-context pointers.  This is intentionally not an installed
 * qtcad header; public qtcad code should use QgLegacyView.h helpers instead.
 */

#ifndef QGLEGACYVIEWCONTEXT_H
#define QGLEGACYVIEWCONTEXT_H

#include "bv.h"
#include "qtcad/QgLegacyView.h"

static inline qg_legacy_view *
qg_legacy_view_from_context(void *view_ctx)
{
    return reinterpret_cast<qg_legacy_view *>(view_ctx);
}

static inline const qg_legacy_view *
qg_legacy_view_from_context(const void *view_ctx)
{
    return reinterpret_cast<const qg_legacy_view *>(view_ctx);
}

static inline void *
qg_legacy_view_to_context(qg_legacy_view *view)
{
    return reinterpret_cast<void *>(view);
}

static inline const void *
qg_legacy_view_to_context(const qg_legacy_view *view)
{
    return reinterpret_cast<const void *>(view);
}

static inline struct bv_context *
qg_legacy_view_context(qg_legacy_view *view)
{
    return reinterpret_cast<struct bv_context *>(view);
}

static inline const struct bv_context *
qg_legacy_view_context_const(const qg_legacy_view *view)
{
    return reinterpret_cast<const struct bv_context *>(view);
}

static inline struct bv *
qg_legacy_view_bv(qg_legacy_view *view)
{
    return bv_context_view(qg_legacy_view_context(view));
}

static inline const struct bv *
qg_legacy_view_bv_const(const qg_legacy_view *view)
{
    return bv_context_view_const(qg_legacy_view_context_const(view));
}

static inline struct bv_context *
qg_legacy_context_bv_context(void *view_ctx)
{
    return reinterpret_cast<struct bv_context *>(view_ctx);
}

static inline const struct bv_context *
qg_legacy_context_bv_context_const(const void *view_ctx)
{
    return reinterpret_cast<const struct bv_context *>(view_ctx);
}

static inline struct bv *
qg_legacy_context_bv(void *view_ctx)
{
    return bv_context_view(qg_legacy_context_bv_context(view_ctx));
}

static inline const struct bv *
qg_legacy_context_bv_const(const void *view_ctx)
{
    return bv_context_view_const(qg_legacy_context_bv_context_const(view_ctx));
}

#endif /* QGLEGACYVIEWCONTEXT_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
