/*        D R A W _ O B O L _ O V E R L A Y _ P R I V A T E . H P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file draw_obol_overlay_private.hpp
 *
 * Hidden GED/Inventor value-boundary conversions shared by the retained
 * database bridge and overlay publication unit.
 */

#ifndef LIBGED_DRAW_OBOL_OVERLAY_PRIVATE_HPP
#define LIBGED_DRAW_OBOL_OVERLAY_PRIVATE_HPP

#include "vmath.h"

#include <Inventor/SbColor.h>
#include <Inventor/SbMatrix.h>

void ged_obol_rgb_from_color(const SbColor &color, unsigned char rgb[3]);
void ged_obol_mat_from_sbmatrix(const SbMatrix &matrix, mat_t mat);
SbMatrix ged_obol_sbmatrix_from_mat(const mat_t mat);
SbColor ged_obol_color_from_rgb(const unsigned char rgb[3]);
float ged_obol_transparency_from_appearance_opacity(fastf_t opacity);
fastf_t ged_obol_appearance_opacity_from_transparency(float transparency);
fastf_t ged_obol_reported_transparency(float transparency);

#endif /* LIBGED_DRAW_OBOL_OVERLAY_PRIVATE_HPP */
