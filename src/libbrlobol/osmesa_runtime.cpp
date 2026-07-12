/*                 O S M E S A _ R U N T I M E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

/* Obol's dual-GL shared package historically recorded libosmesa as a private
 * link dependency.  Keep the renderer runtime explicit at libbrlobol's owning
 * boundary so --as-needed consumers do not depend on an accidental RUNPATH. */
extern "C" void *OSMesaGetProcAddress(const char *);

void
brlobol_osmesa_runtime_link(void)
{
    (void)OSMesaGetProcAddress("glGetString");
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
