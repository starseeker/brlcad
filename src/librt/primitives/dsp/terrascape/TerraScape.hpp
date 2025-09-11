#pragma once

/*
 * TerraScape - Advanced Terrain Mesh Generation Library
 * 
 * This library provides efficient grid-to-mesh conversion with advanced
 * Delaunay triangulation and greedy refinement algorithms.
 * 
 * Based on the Terra and Scape terrain simplification algorithms with
 * modern C++ implementation using robust geometric predicates.
 */

#include <vector>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <limits>
#include <queue>
#include <memory>
#include <iostream>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__linux__)
#include <unistd.h>
#elif defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

namespace TerraScape {

// --- Strategy selection for mesh refinement ---
enum class MeshRefineStrategy {
    AUTO,      // Automatically select method based on grid size and RAM
    HEAP,      // Use heap + affected update (optimal, high memory) - uses detria
    HYBRID,    // Hybrid: sparse sampling, then local refinement - uses detria  
    SPARSE     // Sparse sampling only (fastest, lowest quality) - uses detria
};

// --- Helper: RAM detection (platform-specific) ---
inline size_t get_available_ram_bytes() {
#if defined(_WIN32)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullAvailPhys;
#elif defined(__linux__)
    long pages = sysconf(_SC_AVPHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return size_t(pages) * size_t(page_size);
#elif defined(__APPLE__)
    int64_t mem = 0;
    size_t len = sizeof(mem);
    sysctlbyname("hw.memsize", &mem, &len, NULL, 0);
    return size_t(mem);
#else
    // Fallback: Assume 1GB available
    return size_t(1) << 30;
#endif
}

// --- Core data structures ---
struct Vertex { 
    float x, y, z; 
};

struct Triangle { 
    int v0, v1, v2; 
};

struct MeshResult {
    std::vector<Vertex> vertices;
    std::vector<Triangle> triangles;
    bool is_volumetric = false;  // Indicates if this is a volumetric (closed manifold) mesh
};

// Result structure for volumetric mesh generation that separates positive and negative volumes
struct VolumetricMeshResult {
    MeshResult positive_volume;  // Mesh for areas where surface height > z_base
    MeshResult negative_volume;  // Mesh for areas where surface height < z_base (with reversed normals)
    bool has_positive_volume = false;
    bool has_negative_volume = false;
};

// --- Forward declarations for implementation ---
class DetriaTriangulationManager;
class GreedyMeshRefiner;

// --- Helper: Enhanced grid accessor ---
template<typename T>
class Grid {
public:
    int width, height;
    const T* data;
    
    Grid(int w, int h, const T* d) : width(w), height(h), data(d) {}
    
    T operator()(int x, int y) const { 
        return data[y * width + x]; 
    }
    
    // Bilinear interpolation
    float interp(float x, float y) const {
        int ix = std::max(0, std::min(int(x), width - 2));
        int iy = std::max(0, std::min(int(y), height - 2));
        float fx = x - ix, fy = y - iy;
        float h00 = float((*this)(ix, iy)), h10 = float((*this)(ix + 1, iy));
        float h01 = float((*this)(ix, iy + 1)), h11 = float((*this)(ix + 1, iy + 1));
        return h00 * (1 - fx) * (1 - fy) + h10 * fx * (1 - fy) +
               h01 * (1 - fx) * fy + h11 * fx * fy;
    }
};

// --- Utility: barycentric triangle routines ---
inline bool point_in_triangle(float px, float py,
                              const Vertex& v0, const Vertex& v1, const Vertex& v2)
{
    float dX = px - v2.x;
    float dY = py - v2.y;
    float dX21 = v2.x - v1.x;
    float dY12 = v1.y - v2.y;
    float D = (v0.x - v2.x) * dY12 + (v0.y - v2.y) * dX21;
    float s = ((v0.x - v2.x) * dY + (v0.y - v2.y) * dX) / D;
    float t = (dY12 * dX + dX21 * dY) / D;
    return (s >= 0) && (t >= 0) && (s + t <= 1);
}

inline float barycentric_interp(float px, float py,
                                const Vertex& v0, const Vertex& v1, const Vertex& v2)
{
    float denom = ((v1.y - v2.y)*(v0.x - v2.x) + (v2.x - v1.x)*(v0.y - v2.y));
    if (std::abs(denom) < 1e-10f) return v0.z; // degenerate triangle
    
    float a = ((v1.y - v2.y)*(px - v2.x) + (v2.x - v1.x)*(py - v2.y)) / denom;
    float b = ((v2.y - v0.y)*(px - v2.x) + (v0.x - v2.x)*(py - v2.y)) / denom;
    float c = 1 - a - b;
    return a*v0.z + b*v1.z + c*v2.z;
}

// --- Forward declaration of implementation function ---
template<typename T>
MeshResult grid_to_mesh_detria(
    int width, int height, const T* elevations,
    float error_threshold, int point_limit,
    MeshRefineStrategy strategy);

// --- Helper functions for volumetric mesh generation ---
struct Edge {
    int v0, v1;
    
    Edge(int a, int b) : v0(std::min(a, b)), v1(std::max(a, b)) {}
    
    bool operator<(const Edge& other) const {
        return std::tie(v0, v1) < std::tie(other.v0, other.v1);
    }
    
    bool operator==(const Edge& other) const {
        return v0 == other.v0 && v1 == other.v1;
    }
};

// Find boundary edges of a surface mesh (edges that belong to only one triangle)
std::vector<Edge> find_boundary_edges(const std::vector<Triangle>& triangles);

// Convert surface mesh to volumetric mesh by adding base and side faces
MeshResult make_volumetric_mesh(const MeshResult& surface_mesh, float z_base);

// Convert surface mesh to volumetric mesh with separated positive/negative volumes
VolumetricMeshResult make_volumetric_mesh_separated(const MeshResult& surface_mesh, float z_base);

// --- Forward declaration of volumetric mesh function ---
template<typename T>
MeshResult grid_to_mesh_volumetric(
    int width, int height, const T* elevations,
    float z_base = 0.0f,
    float error_threshold = 1.0f, int point_limit = 10000,
    MeshRefineStrategy strategy = MeshRefineStrategy::AUTO);

// --- Forward declaration of separated volumetric mesh function ---
template<typename T>
VolumetricMeshResult grid_to_mesh_volumetric_separated(
    int width, int height, const T* elevations,
    float z_base = 0.0f,
    float error_threshold = 1.0f, int point_limit = 10000,
    MeshRefineStrategy strategy = MeshRefineStrategy::AUTO);

// --- Main API: Core mesh generation with detria integration ---
// This function provides the primary interface for converting elevation grids
// to triangle meshes using advanced Delaunay triangulation and greedy refinement.
template<typename T>
MeshResult grid_to_mesh(
    int width, int height, const T* elevations,
    float error_threshold = 1.0f, int point_limit = 10000,
    MeshRefineStrategy strategy = MeshRefineStrategy::AUTO)
{
    // --- Strategy selection ---
    if (strategy == MeshRefineStrategy::AUTO) {
        size_t grid_size = size_t(width) * size_t(height);
        size_t avail_ram = get_available_ram_bytes();
        // More conservative memory estimation for stability
        size_t heap_mem_needed = grid_size * 32;
        
        // Use extremely conservative thresholds to ensure robustness with Hawaii-scale datasets  
        // HYBRID/HEAP strategies have assertion failures with larger datasets, so prioritize SPARSE
        if (grid_size <= 500 && heap_mem_needed < avail_ram / 8) {  // Only very tiny grids
            strategy = MeshRefineStrategy::HEAP;
        } else if (grid_size <= 1000 && heap_mem_needed < avail_ram / 4) {  // Small grids only
            strategy = MeshRefineStrategy::HYBRID;
        } else {
            strategy = MeshRefineStrategy::SPARSE;  // Use SPARSE for virtually all real datasets
        }
        
        std::cout << "AUTO strategy selected: ";
        switch (strategy) {
            case MeshRefineStrategy::HEAP: std::cout << "HEAP"; break;
            case MeshRefineStrategy::HYBRID: std::cout << "HYBRID"; break;
            case MeshRefineStrategy::SPARSE: std::cout << "SPARSE"; break;
            default: break;
        }
        std::cout << " for " << grid_size << " points\n";
    }
    
    // Delegate to detria-based implementation with robust error handling
    try {
        return grid_to_mesh_detria(width, height, elevations, error_threshold, point_limit, strategy);
    } catch (...) {
        // If the selected strategy fails (e.g., assertion in triangulation), 
        // fall back to SPARSE strategy which is more robust for large datasets
        if (strategy != MeshRefineStrategy::SPARSE) {
            std::cout << "Strategy failed, falling back to SPARSE for robustness\n";
            return grid_to_mesh_detria(width, height, elevations, error_threshold, point_limit, MeshRefineStrategy::SPARSE);
        } else {
            // If even SPARSE fails, re-throw the exception
            throw;
        }
    }
}

} // namespace TerraScape

// --- Implementation section ---
// Implementation is in TerraScapeImpl.cpp