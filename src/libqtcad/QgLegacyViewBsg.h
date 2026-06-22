/*              Q G L E G A C Y V I E W B S G . H
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
/** @file QgLegacyViewBsg.h
 *
 * Private adapter between qtcad's opaque legacy view handle and the retained
 * staged BSG view pointer.  This is intentionally not an installed qtcad
 * header; public qtcad code should use QgLegacyView.h helpers instead.
 */

#ifndef QGLEGACYVIEWBSG_H
#define QGLEGACYVIEWBSG_H

#include "qtcad/QgLegacyView.h"

struct bsg_view;

static inline qg_legacy_view *
qg_legacy_view_from_bsg(struct bsg_view *view)
{
    return reinterpret_cast<qg_legacy_view *>(view);
}

static inline const qg_legacy_view *
qg_legacy_view_from_bsg(const struct bsg_view *view)
{
    return reinterpret_cast<const qg_legacy_view *>(view);
}

static inline struct bsg_view *
qg_legacy_view_to_bsg(qg_legacy_view *view)
{
    return reinterpret_cast<struct bsg_view *>(view);
}

static inline const struct bsg_view *
qg_legacy_view_to_bsg(const qg_legacy_view *view)
{
    return reinterpret_cast<const struct bsg_view *>(view);
}

#endif /* QGLEGACYVIEWBSG_H */

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
