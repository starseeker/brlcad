/*                 T E R R A S C A P E _ C P P . C P P
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
/** @file terrascape_cpp.cpp
 *
 * C++ implementation wrapper for TerraScape library.
 * Provides C-compatible interface to the C++ TerraScape library.
 *
 */

#include "common.h"

extern "C" {
#include "bu/malloc.h"
#include "bu/log.h"
}

#include "TerraScape.hpp"
#include <vector>
#include <stdexcept>

extern "C" {

struct terrascape_cpp_result {
    float *vertices;    /* x,y,z coordinates */
    int *triangles;     /* triangle vertex indices */
    size_t num_vertices;
    size_t num_triangles;
};

/**
 * C++ wrapper function that calls TerraScape library to perform volumetric tessellation
 */
int terrascape_cpp_tessellate(int width, int height, const unsigned short *elevations,
                             float z_base, float error_threshold, int point_limit,
                             struct terrascape_cpp_result *result)
{
    if (!elevations || !result) {
        return -1;
    }
    
    /* Initialize result structure */
    result->vertices = nullptr;
    result->triangles = nullptr;
    result->num_vertices = 0;
    result->num_triangles = 0;
    
    try {
        /* Call TerraScape volumetric tessellation function */
        TerraScape::MeshResult mesh = TerraScape::grid_to_mesh_volumetric<unsigned short>(
            width, height, elevations,
            z_base, error_threshold, point_limit,
            TerraScape::MeshRefineStrategy::AUTO
        );
        
        /* Copy results to C-compatible format */
        result->num_vertices = mesh.vertices.size();
        result->num_triangles = mesh.triangles.size();
        
        if (result->num_vertices > 0) {
            /* Allocate vertex array (3 floats per vertex) */
            result->vertices = (float*)bu_malloc(result->num_vertices * 3 * sizeof(float), "terrascape cpp vertices");
            if (!result->vertices) {
                return -1;
            }
            
            /* Copy vertex data */
            for (size_t i = 0; i < result->num_vertices; i++) {
                result->vertices[i*3 + 0] = mesh.vertices[i].x;
                result->vertices[i*3 + 1] = mesh.vertices[i].y;
                result->vertices[i*3 + 2] = mesh.vertices[i].z;
            }
        }
        
        if (result->num_triangles > 0) {
            /* Allocate triangle array (3 ints per triangle) */
            result->triangles = (int*)bu_malloc(result->num_triangles * 3 * sizeof(int), "terrascape cpp triangles");
            if (!result->triangles) {
                if (result->vertices) {
                    bu_free(result->vertices, "terrascape cpp vertices");
                    result->vertices = nullptr;
                }
                return -1;
            }
            
            /* Copy triangle data */
            for (size_t i = 0; i < result->num_triangles; i++) {
                result->triangles[i*3 + 0] = mesh.triangles[i].v0;
                result->triangles[i*3 + 1] = mesh.triangles[i].v1;
                result->triangles[i*3 + 2] = mesh.triangles[i].v2;
            }
        }
        
        return 0;
        
    } catch (const std::exception& e) {
        bu_log("TerraScape tessellation failed: %s\n", e.what());
        terrascape_cpp_free_result(result);
        return -1;
    } catch (...) {
        bu_log("TerraScape tessellation failed: unknown exception\n");
        terrascape_cpp_free_result(result);
        return -1;
    }
}

/**
 * Free memory allocated by terrascape_cpp_tessellate
 */
void terrascape_cpp_free_result(struct terrascape_cpp_result *result)
{
    if (!result) return;
    
    if (result->vertices) {
        bu_free(result->vertices, "terrascape cpp vertices");
        result->vertices = nullptr;
    }
    
    if (result->triangles) {
        bu_free(result->triangles, "terrascape cpp triangles");
        result->triangles = nullptr;
    }
    
    result->num_vertices = 0;
    result->num_triangles = 0;
}

} /* extern "C" */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */