/*                      O W N E D . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * Distributed under the terms of the GNU Lesser General Public License
 * (LGPL), version 2.1.
 */

#ifndef LIBBREP_CDT_OWNED_H
#define LIBBREP_CDT_OWNED_H

#include "common.h"

#include <memory>

#include "bu/malloc.h"

struct cdt_free_array {
    void operator()(void *memory) const
    {
	bu_free(memory, "CDT array");
    }
};

template <typename T>
using cdt_owned_array = std::unique_ptr<T[], cdt_free_array>;

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
