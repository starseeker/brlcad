/*                 I M A G E _ S I Z E . C
 * BRL-CAD
 *
 * Copyright (c) 2004-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libimgstream/image_size.c
 *
 * Infer raw image dimensions from conventional file names and pixel counts.
 */

#include "common.h"

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

#include "imgstream/fb_compat.h"

struct image_size_entry {
    size_t width;
    size_t height;
};

static const struct image_size_entry common_image_sizes[] = {
    {160, 120}, {192, 128}, {320, 200}, {320, 240}, {384, 256},
    {640, 512}, {640, 400}, {640, 480}, {640, 485}, {640, 486},
    {720, 485}, {720, 486}, {768, 512}, {1024, 768}, {1152, 900},
    {1280, 960}, {1280, 1024}, {1440, 972}, {1536, 1024},
    {1600, 1200}, {1600, 1280}, {3072, 2048}, {3200, 4000},
    {3400, 4400}, {4700, 3300}, {0, 0}
};

int
imgstream_image_name_size(size_t *width, size_t *height, const char *name)
{
    if (!width || !height || !name)
	return 0;

    const char *cp = name;
    while (*cp) {
	cp = strchr(cp, '-');
	if (!cp)
	    break;
	unsigned long parsed_width = 0;
	unsigned long parsed_height = 0;
	if (sscanf(cp, "-w%lu-n%lu", &parsed_width, &parsed_height) == 2 &&
	    parsed_width <= SIZE_MAX && parsed_height <= SIZE_MAX) {
	    *width = (size_t)parsed_width;
	    *height = (size_t)parsed_height;
	    return 1;
	}
	cp++;
    }
    return 0;
}

int
imgstream_image_size(size_t *width, size_t *height, size_t pixel_count)
{
    if (!width || !height || !pixel_count)
	return 0;

    for (const struct image_size_entry *entry = common_image_sizes;
	 entry->width; entry++) {
	if (pixel_count == entry->width * entry->height) {
	    *width = entry->width;
	    *height = entry->height;
	    return 1;
	}
    }

    size_t low = 1;
    size_t high = pixel_count < SIZE_MAX / 2 ? pixel_count + 1 : SIZE_MAX;
    while (low < high) {
	size_t midpoint = low + (high - low) / 2;
	if (midpoint <= pixel_count / midpoint)
	    low = midpoint + 1;
	else
	    high = midpoint;
    }
    size_t root = low - 1;
    if (root == pixel_count / root && root * root == pixel_count) {
	*width = root;
	*height = root;
	return 1;
    }
    return 0;
}

int
imgstream_image_file_size(size_t *width, size_t *height,
	const char *filename, int bytes_per_pixel)
{
    if (!width || !height)
	return 0;
    *width = 0;
    *height = 0;
    if (!filename || !filename[0] || bytes_per_pixel <= 0)
	return 0;

    const char *basename = strrchr(filename, '/');
    basename = basename ? basename + 1 : filename;
    if (imgstream_image_name_size(width, height, basename))
	return 1;

    struct stat file_info;
    if (stat(filename, &file_info) != 0 || file_info.st_size < 0)
	return 0;
    return imgstream_image_size(width, height,
	(size_t)file_info.st_size / (size_t)bytes_per_pixel);
}
