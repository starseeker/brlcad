/*                       D A T U M . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @addtogroup rt_datum */
/** @{ */
/** @file rt/primitives/datum.h */

#ifndef RT_PRIMITIVES_DATUM_H
#define RT_PRIMITIVES_DATUM_H

#include "common.h"
#include "bu/vls.h"
#include "rt/defines.h"
#include "rt/geom.h"

__BEGIN_DECLS

/** Resolve RT_DATUM_AUTO using the historical pnt/dir/w convention. */
RT_EXPORT extern rt_datum_type rt_datum_resolved_type(
    const struct rt_datum_internal *datum);

/** Validate a NULL-terminated datum chain.  Returns zero when valid.  If
 * messages is non-NULL, validation failures are appended to it. */
RT_EXPORT extern int rt_datum_validate(
    const struct rt_datum_internal *datum,
    struct bu_vls *messages);

struct rt_db_internal;
struct rt_primitive_lod_realization;

RT_EXPORT extern int rt_datum_wireframe_line_set(struct rt_primitive_lod_realization *realization, struct rt_db_internal *ip);

/** @} */

__END_DECLS

#endif /* RT_PRIMITIVES_DATUM_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
