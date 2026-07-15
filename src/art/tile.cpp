/*                        T I L E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2022-2026 United States Government as represented by
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

/* interface header */
#include "tile.h"

#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic push
#endif
#if defined(__clang__)
#  pragma clang diagnostic push
#endif
#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic ignored "-Wall"
#  pragma GCC diagnostic ignored "-Wextra"
#endif
#if defined(__clang__)
#  pragma clang diagnostic ignored "-Wdocumentation"
#  pragma clang diagnostic ignored "-Wfloat-equal"
#  pragma clang diagnostic ignored "-Wunused-parameter"
#  pragma clang diagnostic ignored "-Wpedantic"
#  pragma clang diagnostic ignored "-Wignored-qualifiers"
#endif

#include "renderer/modeling/frame/frame.h"
#include "foundation/image/pixel.h"
#include "foundation/image/image.h"
#include "foundation/image/tile.h"

#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic pop
#endif
#if defined(__clang__)
#  pragma clang diagnostic pop
#endif

#include <vector>

#include "bu/parallel.h"
#include "imgstream/fb_compat.h"


using namespace foundation;

extern imgstream_fb_t *fbp;	/* Framebuffer handle */


void
ArtTileCallback::on_tile_end(const renderer::Frame* frame, const size_t tile_x, const size_t tile_y)
{
    foundation::Tile& t = frame->image().tile(tile_x, tile_y);
    const foundation::Tile rgb(t, PixelFormatUInt8);
    const foundation::CanvasProperties& props = frame->image().properties();
    const size_t x_coord = tile_x * props.m_tile_width;
    const size_t top = tile_y * props.m_tile_height;
    const size_t y_coord = props.m_canvas_height - top - rgb.get_height();
    std::vector<unsigned char> pixels(rgb.get_width() * rgb.get_height() * 3);

    for (size_t y = 0; y < rgb.get_height(); y++) {
	const size_t src_y = rgb.get_height() - 1 - y;
	for (size_t x = 0; x < rgb.get_width(); x++) {
	    const unsigned char *src = rgb.get_storage() +
		(src_y * rgb.get_width() + x) * 4;
	    unsigned char *dst = pixels.data() + (y * rgb.get_width() + x) * 3;
	    dst[0] = src[0];
	    dst[1] = src[1];
	    dst[2] = src[2];
	}
    }

    bu_semaphore_acquire(BU_SEM_SYSCALL);
    (void)imgstream_fb_writerect(fbp, (int)x_coord, (int)y_coord,
	    (int)rgb.get_width(), (int)rgb.get_height(), pixels.data());
    bu_semaphore_release(BU_SEM_SYSCALL);
}


// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
