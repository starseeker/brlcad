/*              P E R F O R M A N C E _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_PERFORMANCE_PRIVATE_H
#define LIBBOBOL_PERFORMANCE_PRIVATE_H

#include "BObol/BPerformance.h"

#include <stdint.h>

enum BObolPerformanceCounterId {
    BOBOL_PERF_REALIZE_CALLS,
    BOBOL_PERF_REALIZE_TOTAL_US,
    BOBOL_PERF_REALIZE_SEED_US,
    BOBOL_PERF_REALIZE_WALK_US,
    BOBOL_PERF_SOURCES_VISITED,
    BOBOL_PERF_SOURCES_REALIZED,
    BOBOL_PERF_SOURCES_FAILED,
    BOBOL_PERF_WIRE_REALIZE_CALLS,
    BOBOL_PERF_WIRE_REALIZE_US,
    BOBOL_PERF_MESH_REALIZE_CALLS,
    BOBOL_PERF_MESH_REALIZE_US,
    BOBOL_PERF_DIRECT_LEAF_CALLS,
    BOBOL_PERF_DIRECT_LEAF_US,
    BOBOL_PERF_DIRECT_LEAF_REALIZED,
    BOBOL_PERF_DIRECT_LEAF_FAILED,
    BOBOL_PERF_DIRECT_LEAF_FALLBACK,
    BOBOL_PERF_WIRE_CACHE_HITS,
    BOBOL_PERF_WIRE_CACHE_MISSES,
    BOBOL_PERF_MESH_CACHE_HITS,
    BOBOL_PERF_MESH_CACHE_MISSES,
    BOBOL_PERF_PLOT_CALLS,
    BOBOL_PERF_PLOT_US,
    BOBOL_PERF_VLIST_CONVERT_CALLS,
    BOBOL_PERF_VLIST_CONVERT_US,
    BOBOL_PERF_VLIST_POINTS,
    BOBOL_PERF_REALIZED_INSTANCE_NODES,
    BOBOL_PERF_REALIZED_INSTANCE_NODE_US,
    BOBOL_PERF_CAD_COMPACT_ATTEMPTS,
    BOBOL_PERF_CAD_COMPACT_SOURCES,
    BOBOL_PERF_CAD_COMPACT_INSTANCES,
    BOBOL_PERF_CAD_COMPACT_US,
    BOBOL_PERF_SOURCE_REPLACE_CALLS,
    BOBOL_PERF_SOURCE_REPLACE_US,
    BOBOL_PERF_SOURCE_STATE_CALLS,
    BOBOL_PERF_SOURCE_STATE_US,
    BOBOL_PERF_SOURCE_MOVE_CALLS,
    BOBOL_PERF_SOURCE_MOVE_US,
    BOBOL_PERF_SOURCE_INDEX_REBUILD_CALLS,
    BOBOL_PERF_SOURCE_INDEX_REBUILD_US,
    BOBOL_PERF_COUNTER_COUNT
};

void bobol_performance_counter_add(enum BObolPerformanceCounterId id,
    uint64_t value);

int64_t bobol_performance_time_now(void);

class BObolPerformanceTimer {
public:
    explicit BObolPerformanceTimer(enum BObolPerformanceCounterId id) :
	counter(id),
	start(bobol_performance_time_now())
    {
    }

    ~BObolPerformanceTimer(void)
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
	const int64_t elapsed = bobol_performance_time_now() - this->start;
	this->start = 0;
	if (elapsed > 0)
	    bobol_performance_counter_add(this->counter,
		static_cast<uint64_t>(elapsed));
	return elapsed > 0 ? static_cast<uint64_t>(elapsed) : 0;
    }

private:
    enum BObolPerformanceCounterId counter;
    int64_t start;
};

#endif /* LIBBOBOL_PERFORMANCE_PRIVATE_H */
