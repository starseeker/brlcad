/*                 F B S E R V _ L E G A C Y . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file dm/fbserv_legacy.h
 *
 * Narrow adapter from the legacy libdm framebuffer API to libimgstream's
 * framebuffer server backend contract.  New display hosts should install a
 * native fbserv backend instead.
 */

#ifndef DM_FBSERV_LEGACY_H
#define DM_FBSERV_LEGACY_H

#include "common.h"
#include "dm/defines.h"
#include "imgstream/fbserv.h"

__BEGIN_DECLS

DM_EXPORT int dm_fbserv_set_framebuffer(struct fbserv_obj *fbsp,
	struct fb *framebuffer);

__END_DECLS

#endif /* DM_FBSERV_LEGACY_H */
