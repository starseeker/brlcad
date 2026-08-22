/*                   D B _ D I S C O V E R Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */
/** @file rt/db_discovery.h
 *
 * Bounded, read-only database hierarchy discovery.
 *
 * This API is intentionally separate from the traditional database I/O
 * entry points while its large-model behavior is evaluated.  Discovery may
 * decode independent referencing objects in parallel, but publishes normal
 * librt reference counts and callbacks serially.  The returned inventory is
 * immutable and may also seed a client's derived hierarchy without reading
 * every combination a second time.
 */

#ifndef RT_DB_DISCOVERY_H
#define RT_DB_DISCOVERY_H

#include "common.h"

#include <stddef.h>
#include <stdint.h>

#include "vmath.h"
#include "rt/defines.h"
#include "rt/op.h"

__BEGIN_DECLS

struct db_i;
struct directory;
struct rt_db_hierarchy;

struct rt_db_discovery_options {
    /** Zero selects a bounded automatic worker count. */
    size_t max_workers;

    /**
     * Bound the sum of serialized referencing records being decoded.  Zero
     * selects a conservative automatic limit.  A single record larger than
     * the limit is still admitted, alone.
     */
    size_t max_inflight_bytes;
};

struct rt_db_discovery_stats {
    uint64_t directory_microseconds;
    uint64_t decode_microseconds;
    uint64_t publish_microseconds;
    uint64_t object_count;
    uint64_t referencing_object_count;
    uint64_t reference_count;
    uint32_t worker_count;
    uint32_t used_mapped_view;
};

struct rt_db_hierarchy_arc {
    struct directory *parent;
    struct directory *child;
    const char *child_name;
    db_op_t operation;
    int matrix_valid;
    mat_t matrix;
};

RT_EXPORT extern void
rt_db_discovery_options_init(struct rt_db_discovery_options *options);

/**
 * Build a database directory and discover its reference hierarchy.
 *
 * The database must not be edited while this call is active.  The directory
 * scan remains ordered; independent referencing records are decoded in a
 * bounded scatter/gather pass.  On success d_nref and registered
 * dbi_update_nref_t callbacks have the same ordinary state produced by
 * db_update_nref().
 */
RT_EXPORT extern int
rt_db_discovery_build(struct db_i *dbip,
		      const struct rt_db_discovery_options *options,
		      struct rt_db_hierarchy **hierarchy,
		      struct rt_db_discovery_stats *stats);

RT_EXPORT extern void
rt_db_hierarchy_destroy(struct rt_db_hierarchy *hierarchy);

RT_EXPORT extern size_t
rt_db_hierarchy_arc_count(const struct rt_db_hierarchy *hierarchy);

/**
 * Copy arc @p index into caller-owned storage.  child_name remains borrowed
 * from @p hierarchy and is valid until rt_db_hierarchy_destroy().
 */
RT_EXPORT extern int
rt_db_hierarchy_arc_get(const struct rt_db_hierarchy *hierarchy,
			size_t index,
			struct rt_db_hierarchy_arc *arc);

__END_DECLS

#endif /* RT_DB_DISCOVERY_H */

