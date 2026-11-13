/*                 G E D _ O B O L _ O U T P U T _ P R I V A T E . H P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
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

#ifndef GED_OBOL_OUTPUT_PRIVATE_HPP
#define GED_OBOL_OUTPUT_PRIVATE_HPP

#include "common.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <Inventor/SbString.h>

constexpr size_t GED_OBOL_LINE_PATTERN_BITS = 16u;
constexpr double GED_OBOL_MINIMUM_OUTPUT_UNIT = 0.001;

struct ged_obol_dash_pattern {
    bool visible = true;
    std::vector<double> lengths;
    double phase = 0.0;
};

inline bool
ged_obol_output_is_annotation(const SbString &geometryKind)
{
    return geometryKind == "annotation";
}

inline float
ged_obol_output_composite_channel(float foreground, float background,
	float transparency)
{
    const float alpha = 1.0f - (std::max)(0.0f,
	(std::min)(1.0f, transparency));
    return foreground * alpha + background * (1.0f - alpha);
}

inline bool
ged_obol_output_pattern_draws(uint16_t pattern, uint16_t factor,
	size_t rasterStep)
{
    if (!pattern)
	return false;
    const size_t repeat = (std::max)(size_t(1), static_cast<size_t>(factor));
    const size_t bit = (rasterStep / repeat) % GED_OBOL_LINE_PATTERN_BITS;
    return (pattern & (uint16_t(1u) << bit)) != 0;
}

inline ged_obol_dash_pattern
ged_obol_output_dash_pattern(uint16_t pattern, uint16_t factor,
	double outputUnit)
{
    ged_obol_dash_pattern result;
    if (!pattern) {
	result.visible = false;
	return result;
    }
    if (pattern == 0xffffu)
	return result;

    const double unit = (std::max)(outputUnit, GED_OBOL_MINIMUM_OUTPUT_UNIT) *
	(static_cast<double>((std::max)(uint16_t(1), factor)));
    size_t firstOn = 0;
    for (; firstOn < GED_OBOL_LINE_PATTERN_BITS; ++firstOn) {
	const size_t previous =
	    (firstOn + GED_OBOL_LINE_PATTERN_BITS - 1u) %
	    GED_OBOL_LINE_PATTERN_BITS;
	if ((pattern & (uint16_t(1u) << firstOn)) &&
		!(pattern & (uint16_t(1u) << previous)))
	    break;
    }
    if (firstOn == GED_OBOL_LINE_PATTERN_BITS)
	return result;

    bool on = true;
    size_t run = 0;
    for (size_t offset = 0; offset < GED_OBOL_LINE_PATTERN_BITS; ++offset) {
	const size_t bit = (firstOn + offset) % GED_OBOL_LINE_PATTERN_BITS;
	const bool value = (pattern & (uint16_t(1u) << bit)) != 0;
	if (value == on) {
	    ++run;
	    continue;
	}
	result.lengths.push_back(static_cast<double>(run) * unit);
	run = 1;
	on = value;
    }
    result.lengths.push_back(static_cast<double>(run) * unit);
    result.phase = static_cast<double>(
	(GED_OBOL_LINE_PATTERN_BITS - firstOn) % GED_OBOL_LINE_PATTERN_BITS) * unit;
    return result;
}

inline const char *
ged_obol_plot_line_mode(uint16_t pattern)
{
    switch (pattern) {
	case 0xffffu: return "solid";
	case 0x1111u: return "dotted";
	case 0x00ffu: return "shortdashed";
	case 0x18ffu: return "dotdashed";
	case 0x28ffu: return "longdashed";
	default: return pattern ? "dotdashed" : "invisible";
    }
}

inline unsigned int
ged_obol_output_pixel_width(float width)
{
    if (!std::isfinite(width) || width <= 1.0f)
	return 1u;
    return static_cast<unsigned int>(std::lround(width));
}

#endif /* GED_OBOL_OUTPUT_PRIVATE_HPP */
