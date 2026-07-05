/*                         O B O L . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file dm/obol.h
 *
 * Public libdm hooks for the Obol display-manager backend.
 */

#ifndef DM_OBOL_H
#define DM_OBOL_H

#include "common.h"

#include "./defines.h"

struct dm;

__BEGIN_DECLS

/**
 * Return the backend-owned BRLObolViewController for an Obol DM as an opaque
 * pointer, or NULL when @p dmp is not an Obol display manager.
 *
 * The controller remains owned by libdm.  This narrow accessor lets higher
 * layers attach their Obol scene producers without making libdm depend on
 * libged or exposing C++ Obol types through the C libdm API.
 */
DM_EXPORT extern void *dm_obol_controller(struct dm *dmp);

__END_DECLS

#endif /* DM_OBOL_H */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
