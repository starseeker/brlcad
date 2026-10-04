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

/* C APIs may free and replace an existing array.  Transfer ownership for
 * the duration of that call and reclaim the returned pointer on every exit.
 */
template <typename T>
class cdt_array_edit {
public:
    explicit cdt_array_edit(cdt_owned_array<T> &array)
	: owner(&array), pointer(array.release()) {}
    cdt_array_edit(cdt_array_edit &&other) noexcept
	: owner(other.owner), pointer(other.pointer)
    {
	other.owner = NULL;
    }
    ~cdt_array_edit()
    {
	if (owner)
	    owner->reset(pointer);
    }
    operator T **() { return &pointer; }

private:
    cdt_array_edit(const cdt_array_edit &) = delete;
    cdt_array_edit &operator=(const cdt_array_edit &) = delete;
    cdt_owned_array<T> *owner;
    T *pointer;
};

template <typename T>
cdt_array_edit<T>
cdt_array_address(cdt_owned_array<T> &array)
{
    return cdt_array_edit<T>(array);
}

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
