/*                         E R T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2008-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libged/dm/ert.cpp
 *
 * Render the active view with the external rt executable into the endpoint
 * framebuffer stream.  ERT owns the command-level policy and callbacks; the
 * process and stream setup is shared with other external image commands.
 */

#include "common.h"

#include "../ged_private.h"


extern "C" int
ged_ert_core(struct ged *gedp, int argc, const char *argv[])
{
    return _ged_external_rt_to_endpoint(gedp, argc, argv, "rt", "ert");
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
