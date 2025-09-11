/*                    T E R R A S C A P E . C
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
/** @file terrascape.c
 *
 * C wrapper for TerraScape volumetric tessellation library.
 * Provides a C interface to the C++ TerraScape library for DSP primitives.
 *
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bu/malloc.h"
#include "bu/log.h"
#include "vmath.h"
#include "rt/geom.h"

/* Include the DSP primitive header for DSP access macro */
#include "../dsp.h"
#include "terrascape.h"

/* We need to include the C++ TerraScape header and provide C++ wrapper functions */
#ifdef __cplusplus
extern "C" {
#endif

/* Forward declarations for C++ wrapper functions */
struct terrascape_cpp_result {
    float *vertices;    /* x,y,z coordinates */
    int *triangles;     /* triangle vertex indices */
    size_t num_vertices;
    size_t num_triangles;
};

int terrascape_cpp_tessellate(int width, int height, const unsigned short *elevations,
                             float z_base, float error_threshold, int point_limit,
                             struct terrascape_cpp_result *result);
void terrascape_cpp_free_result(struct terrascape_cpp_result *result);

#ifdef __cplusplus
}
#endif


void
terrascape_params_init(struct terrascape_params *params)
{
    if (!params) return;
    
    params->magic = TERRASCAPE_PARAMS_MAGIC;
    params->tolerance = 1.0;        /* Default error threshold */
    params->point_limit = 10000;    /* Default point limit */
    params->z_base = 0.0;           /* Default base level */
    params->ensure_watertight = 1;
}


void
terrascape_bot_data_free(struct terrascape_bot_data *bot_data)
{
    if (!bot_data) return;
    
    if (bot_data->vertices) {
        bu_free(bot_data->vertices, "terrascape vertices");
        bot_data->vertices = NULL;
    }
    
    if (bot_data->faces) {
        bu_free(bot_data->faces, "terrascape faces");
        bot_data->faces = NULL;
    }
    
    bot_data->num_vertices = 0;
    bot_data->num_faces = 0;
    bot_data->magic = 0;
}


/**
 * Perform volumetric tessellation of DSP data using the TerraScape library.
 * This function serves as a C wrapper around the C++ TerraScape API.
 */
int
terrascape_tessellate_volumetric(struct rt_dsp_internal *dsp_ip,
                                 struct terrascape_params *params,
                                 struct terrascape_bot_data *bot_data)
{
    struct terrascape_cpp_result cpp_result;
    int ret;
    size_t i;

    if (!dsp_ip || !params || !bot_data) {
        bu_log("terrascape_tessellate_volumetric: NULL pointer(s)\n");
        return -1;
    }
    
    if (params->magic != TERRASCAPE_PARAMS_MAGIC) {
        bu_log("terrascape_tessellate_volumetric: bad params magic\n");
        return -1;
    }
    
    if (!dsp_ip->dsp_buf) {
        bu_log("terrascape_tessellate_volumetric: no DSP data buffer\n");
        return -1;
    }
    
    /* Initialize output structure */
    bot_data->magic = TERRASCAPE_BOT_DATA_MAGIC;
    bot_data->num_vertices = 0;
    bot_data->num_faces = 0;
    bot_data->vertices = NULL;
    bot_data->faces = NULL;
    
    int xcnt = dsp_ip->dsp_xcnt;
    int ycnt = dsp_ip->dsp_ycnt;
    
    if (xcnt < 2 || ycnt < 2) {
        bu_log("terrascape_tessellate_volumetric: DSP too small (%dx%d)\n", xcnt, ycnt);
        return -1;
    }
    
    /* Call the C++ TerraScape wrapper function */
    ret = terrascape_cpp_tessellate(xcnt, ycnt, dsp_ip->dsp_buf,
                                   (float)params->z_base, 
                                   (float)params->tolerance,
                                   params->point_limit,
                                   &cpp_result);
    
    if (ret != 0) {
        bu_log("terrascape_tessellate_volumetric: TerraScape tessellation failed\n");
        return -1;
    }
    
    /* Convert the result from C++ to our C structure */
    bot_data->num_vertices = cpp_result.num_vertices;
    bot_data->num_faces = cpp_result.num_triangles;
    
    /* Allocate memory for vertices and faces */
    if (bot_data->num_vertices > 0) {
        bot_data->vertices = (point_t *)bu_calloc(bot_data->num_vertices, sizeof(point_t), "terrascape vertices");
        if (!bot_data->vertices) {
            bu_log("terrascape_tessellate_volumetric: vertex allocation failed\n");
            terrascape_cpp_free_result(&cpp_result);
            return -1;
        }
        
        /* Copy vertices from float array to point_t array */
        for (i = 0; i < bot_data->num_vertices; i++) {
            bot_data->vertices[i][X] = cpp_result.vertices[i*3 + 0];
            bot_data->vertices[i][Y] = cpp_result.vertices[i*3 + 1];
            bot_data->vertices[i][Z] = cpp_result.vertices[i*3 + 2];
        }
    }
    
    if (bot_data->num_faces > 0) {
        bot_data->faces = (int *)bu_calloc(bot_data->num_faces * 3, sizeof(int), "terrascape faces");
        if (!bot_data->faces) {
            bu_log("terrascape_tessellate_volumetric: face allocation failed\n");
            terrascape_bot_data_free(bot_data);
            terrascape_cpp_free_result(&cpp_result);
            return -1;
        }
        
        /* Copy triangle indices */
        for (i = 0; i < bot_data->num_faces * 3; i++) {
            bot_data->faces[i] = cpp_result.triangles[i];
        }
    }
    
    /* Clean up the C++ result */
    terrascape_cpp_free_result(&cpp_result);
    
    return 0;
}


/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
    params->tolerance = 0.1;
    params->ensure_watertight = 1;
}


void
terrascape_bot_data_free(struct terrascape_bot_data *bot_data)
{
    if (!bot_data) return;
    
    if (bot_data->vertices) {
        bu_free(bot_data->vertices, "terrascape vertices");
        bot_data->vertices = NULL;
    }
    
    if (bot_data->faces) {
        bu_free(bot_data->faces, "terrascape faces");
        bot_data->faces = NULL;
    }
    
    bot_data->num_vertices = 0;
    bot_data->num_faces = 0;
    bot_data->magic = 0;
}


/**
 * Generate a simple volumetric tessellation of DSP data.
 * This creates a basic volumetric mesh by extruding the DSP surface
 * down to a base level, ensuring watertight topology.
 */
int
terrascape_tessellate_volumetric(struct rt_dsp_internal *dsp_ip,
                                 struct terrascape_params *params,
                                 struct terrascape_bot_data *bot_data)
{
    if (!dsp_ip || !params || !bot_data) {
        bu_log("terrascape_tessellate_volumetric: NULL pointer(s)\n");
        return -1;
    }
    
    if (params->magic != TERRASCAPE_PARAMS_MAGIC) {
        bu_log("terrascape_tessellate_volumetric: bad params magic\n");
        return -1;
    }
    
    if (!dsp_ip->dsp_buf) {
        bu_log("terrascape_tessellate_volumetric: no DSP data buffer\n");
        return -1;
    }
    
    /* Initialize output structure */
    bot_data->magic = TERRASCAPE_BOT_DATA_MAGIC;
    bot_data->num_vertices = 0;
    bot_data->num_faces = 0;
    bot_data->vertices = NULL;
    bot_data->faces = NULL;
    
    int xcnt = dsp_ip->dsp_xcnt;
    int ycnt = dsp_ip->dsp_ycnt;
    
    if (xcnt < 2 || ycnt < 2) {
        bu_log("terrascape_tessellate_volumetric: DSP too small (%dx%d)\n", xcnt, ycnt);
        return -1;
    }
    
    /* Calculate number of vertices and faces for volumetric mesh
     * We create two layers: top surface matching DSP heights, bottom at z=0
     * Plus side walls to make it watertight
     */
    size_t top_vertices = xcnt * ycnt;
    size_t bottom_vertices = xcnt * ycnt;
    size_t total_vertices = top_vertices + bottom_vertices;
    
    /* Faces: 
     * - Top surface: 2 triangles per cell = 2 * (xcnt-1) * (ycnt-1)
     * - Bottom surface: 2 triangles per cell = 2 * (xcnt-1) * (ycnt-1)  
     * - Side walls: 2 triangles per edge = 2 * (2 * (xcnt-1) + 2 * (ycnt-1))
     */
    size_t top_faces = 2 * (xcnt - 1) * (ycnt - 1);
    size_t bottom_faces = 2 * (xcnt - 1) * (ycnt - 1);
    size_t side_faces = 2 * (2 * (xcnt - 1) + 2 * (ycnt - 1));
    size_t total_faces = top_faces + bottom_faces + side_faces;
    
    /* Allocate memory */
    bot_data->vertices = (point_t *)bu_calloc(total_vertices, sizeof(point_t), "terrascape vertices");
    bot_data->faces = (int *)bu_calloc(total_faces * 3, sizeof(int), "terrascape faces");
    
    if (!bot_data->vertices || !bot_data->faces) {
        bu_log("terrascape_tessellate_volumetric: memory allocation failed\n");
        terrascape_bot_data_free(bot_data);
        return -1;
    }
    
    bot_data->num_vertices = total_vertices;
    bot_data->num_faces = total_faces;
    
    /* Generate vertices */
    size_t vert_idx = 0;
    
    /* Top surface vertices (matching DSP heights) */
    for (int y = 0; y < ycnt; y++) {
        for (int x = 0; x < xcnt; x++) {
            unsigned short height = DSP(dsp_ip, x, y);
            VSET(bot_data->vertices[vert_idx], x, y, height);
            vert_idx++;
        }
    }
    
    /* Bottom surface vertices (at z=0) */
    for (int y = 0; y < ycnt; y++) {
        for (int x = 0; x < xcnt; x++) {
            VSET(bot_data->vertices[vert_idx], x, y, 0);
            vert_idx++;
        }
    }
    
    /* Generate faces */
    size_t face_idx = 0;
    
    /* Top surface faces */
    for (int y = 0; y < ycnt - 1; y++) {
        for (int x = 0; x < xcnt - 1; x++) {
            int v0 = y * xcnt + x;
            int v1 = y * xcnt + (x + 1);
            int v2 = (y + 1) * xcnt + x;
            int v3 = (y + 1) * xcnt + (x + 1);
            
            /* First triangle */
            bot_data->faces[face_idx * 3 + 0] = v0;
            bot_data->faces[face_idx * 3 + 1] = v1;
            bot_data->faces[face_idx * 3 + 2] = v2;
            face_idx++;
            
            /* Second triangle */
            bot_data->faces[face_idx * 3 + 0] = v1;
            bot_data->faces[face_idx * 3 + 1] = v3;
            bot_data->faces[face_idx * 3 + 2] = v2;
            face_idx++;
        }
    }
    
    /* Bottom surface faces (opposite orientation) */
    for (int y = 0; y < ycnt - 1; y++) {
        for (int x = 0; x < xcnt - 1; x++) {
            int v0 = top_vertices + y * xcnt + x;
            int v1 = top_vertices + y * xcnt + (x + 1);
            int v2 = top_vertices + (y + 1) * xcnt + x;
            int v3 = top_vertices + (y + 1) * xcnt + (x + 1);
            
            /* First triangle (reversed winding) */
            bot_data->faces[face_idx * 3 + 0] = v0;
            bot_data->faces[face_idx * 3 + 1] = v2;
            bot_data->faces[face_idx * 3 + 2] = v1;
            face_idx++;
            
            /* Second triangle (reversed winding) */
            bot_data->faces[face_idx * 3 + 0] = v1;
            bot_data->faces[face_idx * 3 + 1] = v2;
            bot_data->faces[face_idx * 3 + 2] = v3;
            face_idx++;
        }
    }
    
    /* Side faces - connect top and bottom perimeter */
    
    /* X-min edge */
    for (int y = 0; y < ycnt - 1; y++) {
        int top0 = y * xcnt;
        int top1 = (y + 1) * xcnt;
        int bot0 = top_vertices + y * xcnt;
        int bot1 = top_vertices + (y + 1) * xcnt;
        
        /* Triangle 1 */
        bot_data->faces[face_idx * 3 + 0] = top0;
        bot_data->faces[face_idx * 3 + 1] = bot0;
        bot_data->faces[face_idx * 3 + 2] = top1;
        face_idx++;
        
        /* Triangle 2 */
        bot_data->faces[face_idx * 3 + 0] = top1;
        bot_data->faces[face_idx * 3 + 1] = bot0;
        bot_data->faces[face_idx * 3 + 2] = bot1;
        face_idx++;
    }
    
    /* X-max edge */
    for (int y = 0; y < ycnt - 1; y++) {
        int top0 = y * xcnt + (xcnt - 1);
        int top1 = (y + 1) * xcnt + (xcnt - 1);
        int bot0 = top_vertices + y * xcnt + (xcnt - 1);
        int bot1 = top_vertices + (y + 1) * xcnt + (xcnt - 1);
        
        /* Triangle 1 */
        bot_data->faces[face_idx * 3 + 0] = top0;
        bot_data->faces[face_idx * 3 + 1] = top1;
        bot_data->faces[face_idx * 3 + 2] = bot0;
        face_idx++;
        
        /* Triangle 2 */
        bot_data->faces[face_idx * 3 + 0] = top1;
        bot_data->faces[face_idx * 3 + 1] = bot1;
        bot_data->faces[face_idx * 3 + 2] = bot0;
        face_idx++;
    }
    
    /* Y-min edge */
    for (int x = 0; x < xcnt - 1; x++) {
        int top0 = x;
        int top1 = x + 1;
        int bot0 = top_vertices + x;
        int bot1 = top_vertices + x + 1;
        
        /* Triangle 1 */
        bot_data->faces[face_idx * 3 + 0] = top0;
        bot_data->faces[face_idx * 3 + 1] = top1;
        bot_data->faces[face_idx * 3 + 2] = bot0;
        face_idx++;
        
        /* Triangle 2 */
        bot_data->faces[face_idx * 3 + 0] = top1;
        bot_data->faces[face_idx * 3 + 1] = bot1;
        bot_data->faces[face_idx * 3 + 2] = bot0;
        face_idx++;
    }
    
    /* Y-max edge */
    for (int x = 0; x < xcnt - 1; x++) {
        int top0 = (ycnt - 1) * xcnt + x;
        int top1 = (ycnt - 1) * xcnt + x + 1;
        int bot0 = top_vertices + (ycnt - 1) * xcnt + x;
        int bot1 = top_vertices + (ycnt - 1) * xcnt + x + 1;
        
        /* Triangle 1 */
        bot_data->faces[face_idx * 3 + 0] = top0;
        bot_data->faces[face_idx * 3 + 1] = bot0;
        bot_data->faces[face_idx * 3 + 2] = top1;
        face_idx++;
        
        /* Triangle 2 */
        bot_data->faces[face_idx * 3 + 0] = top1;
        bot_data->faces[face_idx * 3 + 1] = bot0;
        bot_data->faces[face_idx * 3 + 2] = bot1;
        face_idx++;
    }
    
    if (face_idx != total_faces) {
        bu_log("terrascape_tessellate_volumetric: face count mismatch (%zu vs %zu)\n", 
               face_idx, total_faces);
        terrascape_bot_data_free(bot_data);
        return -1;
    }
    
    return 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */