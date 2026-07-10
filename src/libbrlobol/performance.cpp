/*                 P E R F O R M A N C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "performance_private.h"

#include <atomic>
#include <string.h>

#include "bu/time.h"

static std::atomic<int> brlobol_perf_enabled(0);
static std::atomic<uint64_t> brlobol_perf_counters[BRLOBOL_PERF_COUNTER_COUNT];

void
brlobol_performance_counters_init(struct BRLObolPerformanceCounters *counters)
{
    if (counters)
	memset(counters, 0, sizeof(*counters));
}

void
brlobol_performance_counters_reset(void)
{
    for (size_t i = 0; i < BRLOBOL_PERF_COUNTER_COUNT; i++)
	brlobol_perf_counters[i].store(0, std::memory_order_relaxed);
}

void
brlobol_performance_counters_set_enabled(int enabled)
{
    brlobol_perf_enabled.store(enabled ? 1 : 0, std::memory_order_relaxed);
}

int
brlobol_performance_counters_enabled(void)
{
    return brlobol_perf_enabled.load(std::memory_order_relaxed);
}

void
brlobol_performance_counter_add(enum BRLObolPerformanceCounterId id,
				uint64_t value)
{
    if (id < 0 || id >= BRLOBOL_PERF_COUNTER_COUNT || value == 0 ||
	!brlobol_performance_counters_enabled())
	return;

    brlobol_perf_counters[id].fetch_add(value, std::memory_order_relaxed);
}

int64_t
brlobol_performance_time_now(void)
{
    return brlobol_performance_counters_enabled() ? bu_gettime() : 0;
}

void
brlobol_performance_counters_get(
    struct BRLObolPerformanceCounters *counters)
{
    if (!counters)
	return;

    brlobol_performance_counters_init(counters);
    uint64_t values[BRLOBOL_PERF_COUNTER_COUNT] = {0};
    for (size_t i = 0; i < BRLOBOL_PERF_COUNTER_COUNT; i++)
	values[i] = brlobol_perf_counters[i].load(std::memory_order_relaxed);

    counters->realize_calls = values[BRLOBOL_PERF_REALIZE_CALLS];
    counters->realize_total_us = values[BRLOBOL_PERF_REALIZE_TOTAL_US];
    counters->realize_seed_us = values[BRLOBOL_PERF_REALIZE_SEED_US];
    counters->realize_walk_us = values[BRLOBOL_PERF_REALIZE_WALK_US];
    counters->sources_visited = values[BRLOBOL_PERF_SOURCES_VISITED];
    counters->sources_realized = values[BRLOBOL_PERF_SOURCES_REALIZED];
    counters->sources_failed = values[BRLOBOL_PERF_SOURCES_FAILED];
    counters->wire_realize_calls = values[BRLOBOL_PERF_WIRE_REALIZE_CALLS];
    counters->wire_realize_us = values[BRLOBOL_PERF_WIRE_REALIZE_US];
    counters->mesh_realize_calls = values[BRLOBOL_PERF_MESH_REALIZE_CALLS];
    counters->mesh_realize_us = values[BRLOBOL_PERF_MESH_REALIZE_US];
    counters->direct_leaf_calls = values[BRLOBOL_PERF_DIRECT_LEAF_CALLS];
    counters->direct_leaf_us = values[BRLOBOL_PERF_DIRECT_LEAF_US];
    counters->direct_leaf_realized = values[BRLOBOL_PERF_DIRECT_LEAF_REALIZED];
    counters->direct_leaf_failed = values[BRLOBOL_PERF_DIRECT_LEAF_FAILED];
    counters->direct_leaf_fallback = values[BRLOBOL_PERF_DIRECT_LEAF_FALLBACK];
    counters->wire_cache_hits = values[BRLOBOL_PERF_WIRE_CACHE_HITS];
    counters->wire_cache_misses = values[BRLOBOL_PERF_WIRE_CACHE_MISSES];
    counters->plot_calls = values[BRLOBOL_PERF_PLOT_CALLS];
    counters->plot_us = values[BRLOBOL_PERF_PLOT_US];
    counters->vlist_convert_calls = values[BRLOBOL_PERF_VLIST_CONVERT_CALLS];
    counters->vlist_convert_us = values[BRLOBOL_PERF_VLIST_CONVERT_US];
    counters->vlist_points = values[BRLOBOL_PERF_VLIST_POINTS];
    counters->realized_instance_nodes =
	values[BRLOBOL_PERF_REALIZED_INSTANCE_NODES];
    counters->realized_instance_node_us =
	values[BRLOBOL_PERF_REALIZED_INSTANCE_NODE_US];
    counters->cad_compact_attempts =
	values[BRLOBOL_PERF_CAD_COMPACT_ATTEMPTS];
    counters->cad_compact_sources =
	values[BRLOBOL_PERF_CAD_COMPACT_SOURCES];
    counters->cad_compact_instances =
	values[BRLOBOL_PERF_CAD_COMPACT_INSTANCES];
    counters->cad_compact_us = values[BRLOBOL_PERF_CAD_COMPACT_US];
    counters->source_replace_calls =
	values[BRLOBOL_PERF_SOURCE_REPLACE_CALLS];
    counters->source_replace_us = values[BRLOBOL_PERF_SOURCE_REPLACE_US];
    counters->source_state_calls = values[BRLOBOL_PERF_SOURCE_STATE_CALLS];
    counters->source_state_us = values[BRLOBOL_PERF_SOURCE_STATE_US];
    counters->source_move_calls = values[BRLOBOL_PERF_SOURCE_MOVE_CALLS];
    counters->source_move_us = values[BRLOBOL_PERF_SOURCE_MOVE_US];
    counters->source_index_rebuild_calls =
	values[BRLOBOL_PERF_SOURCE_INDEX_REBUILD_CALLS];
    counters->source_index_rebuild_us =
	values[BRLOBOL_PERF_SOURCE_INDEX_REBUILD_US];
}
