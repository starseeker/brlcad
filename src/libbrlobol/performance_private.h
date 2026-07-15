/*              P E R F O R M A N C E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBRLOBOL_PERFORMANCE_PRIVATE_H
#define LIBBRLOBOL_PERFORMANCE_PRIVATE_H

#include "brlobol/performance.h"

#include <stdint.h>

enum BRLObolPerformanceCounterId {
    BRLOBOL_PERF_REALIZE_CALLS,
    BRLOBOL_PERF_REALIZE_TOTAL_US,
    BRLOBOL_PERF_REALIZE_SEED_US,
    BRLOBOL_PERF_REALIZE_WALK_US,
    BRLOBOL_PERF_SOURCES_VISITED,
    BRLOBOL_PERF_SOURCES_REALIZED,
    BRLOBOL_PERF_SOURCES_FAILED,
    BRLOBOL_PERF_WIRE_REALIZE_CALLS,
    BRLOBOL_PERF_WIRE_REALIZE_US,
    BRLOBOL_PERF_MESH_REALIZE_CALLS,
    BRLOBOL_PERF_MESH_REALIZE_US,
    BRLOBOL_PERF_DIRECT_LEAF_CALLS,
    BRLOBOL_PERF_DIRECT_LEAF_US,
    BRLOBOL_PERF_DIRECT_LEAF_REALIZED,
    BRLOBOL_PERF_DIRECT_LEAF_FAILED,
    BRLOBOL_PERF_DIRECT_LEAF_FALLBACK,
    BRLOBOL_PERF_WIRE_CACHE_HITS,
    BRLOBOL_PERF_WIRE_CACHE_MISSES,
    BRLOBOL_PERF_MESH_CACHE_HITS,
    BRLOBOL_PERF_MESH_CACHE_MISSES,
    BRLOBOL_PERF_PLOT_CALLS,
    BRLOBOL_PERF_PLOT_US,
    BRLOBOL_PERF_VLIST_CONVERT_CALLS,
    BRLOBOL_PERF_VLIST_CONVERT_US,
    BRLOBOL_PERF_VLIST_POINTS,
    BRLOBOL_PERF_REALIZED_INSTANCE_NODES,
    BRLOBOL_PERF_REALIZED_INSTANCE_NODE_US,
    BRLOBOL_PERF_CAD_COMPACT_ATTEMPTS,
    BRLOBOL_PERF_CAD_COMPACT_SOURCES,
    BRLOBOL_PERF_CAD_COMPACT_INSTANCES,
    BRLOBOL_PERF_CAD_COMPACT_US,
    BRLOBOL_PERF_SOURCE_REPLACE_CALLS,
    BRLOBOL_PERF_SOURCE_REPLACE_US,
    BRLOBOL_PERF_SOURCE_STATE_CALLS,
    BRLOBOL_PERF_SOURCE_STATE_US,
    BRLOBOL_PERF_SOURCE_MOVE_CALLS,
    BRLOBOL_PERF_SOURCE_MOVE_US,
    BRLOBOL_PERF_SOURCE_INDEX_REBUILD_CALLS,
    BRLOBOL_PERF_SOURCE_INDEX_REBUILD_US,
    BRLOBOL_PERF_COUNTER_COUNT
};

void brlobol_performance_counter_add(enum BRLObolPerformanceCounterId id,
    uint64_t value);

int64_t brlobol_performance_time_now(void);

class BRLObolPerformanceTimer {
public:
    explicit BRLObolPerformanceTimer(enum BRLObolPerformanceCounterId id) :
	counter(id),
	start(brlobol_performance_time_now())
    {
    }

    ~BRLObolPerformanceTimer(void)
    {
	this->stop();
    }

    int active(void) const
    {
	return this->start > 0;
    }

    uint64_t stop(void)
    {
	if (this->start <= 0)
	    return 0;
	const int64_t elapsed = brlobol_performance_time_now() - this->start;
	this->start = 0;
	if (elapsed > 0)
	    brlobol_performance_counter_add(this->counter,
		static_cast<uint64_t>(elapsed));
	return elapsed > 0 ? static_cast<uint64_t>(elapsed) : 0;
    }

private:
    enum BRLObolPerformanceCounterId counter;
    int64_t start;
};

#endif /* LIBBRLOBOL_PERFORMANCE_PRIVATE_H */
