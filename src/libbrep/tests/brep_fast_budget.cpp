/*                   B R E P _ F A S T _ B U D G E T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * Distributed under the terms of the GNU Lesser General Public License
 * (LGPL), version 2.1.
 */

#include "common.h"

#include <limits>

#include "bu/app.h"
#include "vmath.h"
#include "../cdt/fast_budget.h"

int
main(int UNUSED(argc), const char **argv)
{
    bu_setprogname(argv[0]);
    /* A replacement with fewer triangles can still need more vertices.
     * Rejecting its point growth must not release bytes to another face. */
    const size_t triangle_bytes = 3 * sizeof(int) + 9 * sizeof(fastf_t);
    const size_t point_bytes = 3 * sizeof(fastf_t);
    const fast_output_size original = {10 * triangle_bytes + 10 * point_bytes, 10};
    const fast_output_size fewer_triangles = {8 * triangle_bytes + 12 * point_bytes, 12};
    const fast_output_size larger = {11 * triangle_bytes + 11 * point_bytes, 11};
    fast_output_budget budget({2 * original.bytes, 21});
    if (budget.replace({}, original) != fast_output_budget::accepted ||
	budget.replace({}, original) != fast_output_budget::accepted ||
	budget.replace(original, fewer_triangles) != fast_output_budget::point_limit ||
	budget.replace(original, larger) != fast_output_budget::byte_limit ||
	budget.replace(original, {}) != fast_output_budget::accepted ||
	budget.replace(original, larger) != fast_output_budget::accepted)
	return 1;

    /* Capacity checks must not overflow when the caller selects SIZE_MAX. */
    const size_t maximum = std::numeric_limits<size_t>::max();
    fast_output_budget large({maximum, maximum});
    if (large.replace({}, {maximum, maximum}) != fast_output_budget::accepted ||
	large.replace({}, {1, 0}) != fast_output_budget::byte_limit ||
	large.replace({}, {0, 1}) != fast_output_budget::point_limit ||
	large.replace({maximum, maximum}, {1, 1}) != fast_output_budget::accepted)
	return 1;
    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
