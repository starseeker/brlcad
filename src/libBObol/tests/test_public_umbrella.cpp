/*              T E S T _ P U B L I C _ U M B R E L L A . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "BObol.h"

int
main(void)
{
    /* This is deliberately a link test as well as an include test.  The
     * property inventory is process-local and does not initialize a display
     * backend, native host, or Coin scene graph. */
    return bobol_display_endpoint_property_count() > 0 ? 0 : 1;
}
