/*                   P E R F O R M A N C E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/performance.h */

#ifndef BRLOBOL_PERFORMANCE_H
#define BRLOBOL_PERFORMANCE_H

#include "brlobol/defines.h"

#include <stdint.h>

__BEGIN_DECLS

struct BRLOBOL_EXPORT BRLObolPerformanceCounters {
    uint64_t realize_calls;
    uint64_t realize_total_us;
    uint64_t realize_seed_us;
    uint64_t realize_walk_us;
    uint64_t sources_visited;
    uint64_t sources_realized;
    uint64_t sources_failed;

    uint64_t wire_realize_calls;
    uint64_t wire_realize_us;
    uint64_t mesh_realize_calls;
    uint64_t mesh_realize_us;

    uint64_t direct_leaf_calls;
    uint64_t direct_leaf_us;
    uint64_t direct_leaf_realized;
    uint64_t direct_leaf_failed;
    uint64_t direct_leaf_fallback;

    uint64_t wire_cache_hits;
    uint64_t wire_cache_misses;
    uint64_t mesh_cache_hits;
    uint64_t mesh_cache_misses;
    uint64_t plot_calls;
    uint64_t plot_us;
    uint64_t vlist_convert_calls;
    uint64_t vlist_convert_us;
    uint64_t vlist_points;
    uint64_t realized_instance_nodes;
    uint64_t realized_instance_node_us;
    uint64_t cad_compact_attempts;
    uint64_t cad_compact_sources;
    uint64_t cad_compact_instances;
    uint64_t cad_compact_us;

    uint64_t source_replace_calls;
    uint64_t source_replace_us;
    uint64_t source_state_calls;
    uint64_t source_state_us;
    uint64_t source_move_calls;
    uint64_t source_move_us;
    uint64_t source_index_rebuild_calls;
    uint64_t source_index_rebuild_us;
};

BRLOBOL_EXPORT void
brlobol_performance_counters_init(struct BRLObolPerformanceCounters *counters);

BRLOBOL_EXPORT void
brlobol_performance_counters_reset(void);

BRLOBOL_EXPORT void
brlobol_performance_counters_set_enabled(int enabled);

BRLOBOL_EXPORT int
brlobol_performance_counters_enabled(void);

BRLOBOL_EXPORT void
brlobol_performance_counters_get(
    struct BRLObolPerformanceCounters *counters);

__END_DECLS

#endif /* BRLOBOL_PERFORMANCE_H */
