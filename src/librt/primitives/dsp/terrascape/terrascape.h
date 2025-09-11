/*                    T E R R A S C A P E . H
 * BRL-CAD
 *
 * Copyright (c) 2025 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file terrascape.h
 *
 * C interface for TerraScape volumetric tessellation library for DSP primitives.
 * Wraps the C++ TerraScape library to provide tessellation functionality.
 *
 */

#ifndef LIBRT_PRIMITIVES_DSP_TERRASCAPE_TERRASCAPE_H
#define LIBRT_PRIMITIVES_DSP_TERRASCAPE_TERRASCAPE_H

#include "common.h"
#include "vmath.h"
#include "rt/geom.h"

__BEGIN_DECLS

/* Structure to hold volumetric tessellation parameters */
struct terrascape_params {
    uint32_t magic;
    double tolerance;       /* Tessellation tolerance/error threshold */
    int point_limit;        /* Maximum number of points in tessellation */
    double z_base;          /* Base level for volumetric extrusion */
    int ensure_watertight;  /* Guarantee watertight output */
};
#define TERRASCAPE_PARAMS_MAGIC 0x54455253

/* Structure to hold generated BOT mesh data */
struct terrascape_bot_data {
    uint32_t magic;
    size_t num_vertices;
    size_t num_faces;
    point_t *vertices;
    int *faces;             /* Array of triangle vertex indices */
};
#define TERRASCAPE_BOT_DATA_MAGIC 0x54424F54

/**
 * Initialize Terrascape tessellation parameters with default values
 */
void terrascape_params_init(struct terrascape_params *params);

/**
 * Perform volumetric tessellation of DSP data using TerraScape library
 *
 * @param dsp_ip Pointer to DSP internal structure
 * @param params Tessellation parameters
 * @param bot_data Output BOT mesh data
 * @return 0 on success, -1 on failure
 */
int terrascape_tessellate_volumetric(struct rt_dsp_internal *dsp_ip,
                                     struct terrascape_params *params,
                                     struct terrascape_bot_data *bot_data);

/**
 * Free BOT data allocated by Terrascape tessellation
 */
void terrascape_bot_data_free(struct terrascape_bot_data *bot_data);

__END_DECLS

#endif /* LIBRT_PRIMITIVES_DSP_TERRASCAPE_TERRASCAPE_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */