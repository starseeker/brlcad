/*                  M A I N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "BObol.h"
#include "bu/app.h"

int
main(void)
{
    bu_setprogname("libbobol_installed_consumer");
    return bobol_display_endpoint_property_count() > 0 ? 0 : 1;
}
