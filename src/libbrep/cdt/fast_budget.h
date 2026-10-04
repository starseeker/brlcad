/*                      F A S T _ B U D G E T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * Distributed under the terms of the GNU Lesser General Public License
 * (LGPL), version 2.1.
 */

#ifndef LIBBREP_CDT_FAST_BUDGET_H
#define LIBBREP_CDT_FAST_BUDGET_H

#include "common.h"

#include <cstddef>
#include <mutex>

#include "bu/assert.h"

struct fast_output_size {
    size_t bytes = 0;
    size_t points = 0;
};

class fast_output_budget {
public:
    enum admission { accepted, byte_limit, point_limit };

    explicit fast_output_budget(fast_output_size limit) : m_limit(limit) {}

    /* Replacing a mesh can reduce bytes while increasing points.  Publish
     * both changes together so another worker cannot consume storage needed
     * to retain the original mesh when either limit rejects its replacement. */
    admission replace(fast_output_size previous, fast_output_size next)
    {
	std::lock_guard<std::mutex> lock(m_mutex);
	BU_ASSERT(previous.bytes <= m_retained.bytes);
	BU_ASSERT(previous.points <= m_retained.points);
	const fast_output_size remaining = {
	    m_retained.bytes - previous.bytes,
	    m_retained.points - previous.points
	};
	if (next.bytes > m_limit.bytes - remaining.bytes)
	    return byte_limit;
	if (next.points > m_limit.points - remaining.points)
	    return point_limit;
	m_retained = {remaining.bytes + next.bytes, remaining.points + next.points};
	return accepted;
    }

private:
    const fast_output_size m_limit;
    fast_output_size m_retained;
    std::mutex m_mutex;
};

#endif

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
