/*                           A N I M . H
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
/** @addtogroup libicv
 *
 * Functions for creating and manipulating animations (APNG, MJPG)
 * provided by the LIBICV image processing library.
 *
 */
/** @{ */
/** @file include/icv/anim.h */

#ifndef ICV_ANIM_H
#define ICV_ANIM_H

#include "common.h"
#include "icv/defines.h"

__BEGIN_DECLS

/** Supported animation formats */
typedef enum {
    ICV_ANIM_UNKNOWN = 0,
    ICV_ANIM_APNG    = 1,
    ICV_ANIM_MJPG    = 2
} icv_anim_format_t;

/** Opaque animation object */
typedef struct icv_anim icv_anim_t;

/** Opaque bounded-memory animation writer. */
typedef struct icv_anim_writer icv_anim_writer_t;

/**
 * Create a new empty animation object.
 *
 * @param format The output format (ICV_ANIM_APNG or ICV_ANIM_MJPG).
 * @param width The canvas width (can be 0 to auto-detect based on frames).
 * @param height The canvas height (can be 0 to auto-detect based on frames).
 * @param fps The default frames per second for the animation.
 */
ICV_EXPORT extern icv_anim_t *icv_anim_create(icv_anim_format_t format, uint32_t width, uint32_t height, int fps);

/**
 * Read an animation from file into an animation object.
 *
 * @param filename The path to the animation file.
 */
ICV_EXPORT extern icv_anim_t *icv_anim_read(const char *filename);

/**
 * Write the animation object to a file.
 *
 * @param anim The animation object.
 * @param filename The output file path.
 */
ICV_EXPORT extern int icv_anim_write(const icv_anim_t *anim, const char *filename);

/**
 * Free an animation object.
 *
 * @param anim The animation object.
 * @return 0 on success.
 */
ICV_EXPORT extern int icv_anim_destroy(icv_anim_t *anim);

/**
 * Append a frame to the end of the animation.
 *
 * @param anim The animation object.
 * @param img The frame image data.
 */
ICV_EXPORT extern int icv_anim_add_frame(icv_anim_t *anim, const icv_image_t *img);

/**
 * Insert a frame into the animation at the specified 0-based index.
 *
 * @param anim The animation object.
 * @param index The 0-based insertion index.
 * @param img The frame image data.
 */
ICV_EXPORT extern int icv_anim_insert_frame(icv_anim_t *anim, size_t index, const icv_image_t *img);

/**
 * Replace a frame in the animation at the specified 0-based index.
 *
 * @param anim The animation object.
 * @param index The 0-based replacement index.
 * @param img The frame image data.
 */
ICV_EXPORT extern int icv_anim_replace_frame(icv_anim_t *anim, size_t index, const icv_image_t *img);

/**
 * Remove a frame from the animation at the specified 0-based index.
 *
 * @param anim The animation object.
 * @param index The 0-based removal index.
 */
ICV_EXPORT extern int icv_anim_remove_frame(icv_anim_t *anim, size_t index);

/**
 * Extract a specific frame from the animation.
 *
 * @param anim The animation object.
 * @param index The 0-based index of the frame to extract.
 * @return A newly allocated icv_image_t containing the frame. Must be freed by caller.
 */
ICV_EXPORT extern icv_image_t *icv_anim_get_frame(const icv_anim_t *anim, size_t index);

/**
 * Get the number of frames currently in the animation.
 *
 * @param anim The animation object.
 * @return The number of frames.
 */
ICV_EXPORT extern size_t icv_anim_num_frames(const icv_anim_t *anim);

/**
 * Set the global frames per second for the animation.
 *
 * @param anim The animation object.
 * @param fps The frames per second.
 * @return 0 on success, -1 on invalid input.
 */
ICV_EXPORT extern int icv_anim_set_fps(icv_anim_t *anim, int fps);

/**
 * Set a specific frame's delay in microseconds. This overrides the global FPS
 * for this specific frame. (Note: Only fully supported by APNG format. MJPG
 * will use the global FPS and generally ignore per-frame delays).
 *
 * @param anim The animation object.
 * @param index The 0-based index of the frame.
 * @param delay_usec The delay in microseconds (1,000,000 = 1 second).
 */
ICV_EXPORT extern int icv_anim_set_frame_delay(icv_anim_t *anim, size_t index, uint32_t delay_usec);

/**
 * Open a streaming APNG writer.  Unlike icv_anim_t, this interface encodes
 * each frame as it arrives and does not retain the complete animation in
 * memory.  The declared frame count must match the number of frames added
 * before close.
 *
 * @param filename Output APNG path.
 * @param width Fixed canvas width.
 * @param height Fixed canvas height.
 * @param frame_count Number of frames that will be added.
 * @param fps Default frame rate used when a frame delay is zero.
 */
ICV_EXPORT extern icv_anim_writer_t *icv_anim_writer_create(
    const char *filename, uint32_t width, uint32_t height,
    size_t frame_count, int fps);

/**
 * Encode one image into an open streaming APNG.
 *
 * The image may be smaller than the fixed canvas and is placed at the given
 * upper-left APNG canvas offset.  It uses source replacement semantics; a
 * caller that needs a complete opaque frame should supply a canvas-sized
 * image.
 *
 * @param writer Open writer.
 * @param img Frame image.
 * @param delay_usec Presentation duration; zero selects the writer FPS.
 * @param x_offset Horizontal canvas offset.
 * @param y_offset Vertical canvas offset.
 * @return 0 on success, -1 on invalid input or encoding failure.
 */
ICV_EXPORT extern int icv_anim_writer_add_frame(
    icv_anim_writer_t *writer, const icv_image_t *img,
    uint32_t delay_usec, uint32_t x_offset, uint32_t y_offset);

/**
 * Encode one top-down, tightly packed RGBA8 frame without converting through
 * libicv's double-precision image storage.  This is the bounded-memory path
 * for native framebuffer captures.
 *
 * @param writer Open writer.
 * @param rgba Top-down width*height*4 byte array.
 * @param width Frame width.
 * @param height Frame height.
 * @param delay_usec Presentation duration; zero selects the writer FPS.
 * @param x_offset Horizontal canvas offset.
 * @param y_offset Vertical canvas offset.
 * @return 0 on success, -1 on invalid input or encoding failure.
 */
ICV_EXPORT extern int icv_anim_writer_add_rgba8(
    icv_anim_writer_t *writer, const unsigned char *rgba,
    uint32_t width, uint32_t height, uint32_t delay_usec,
    uint32_t x_offset, uint32_t y_offset);

/**
 * Finalize and close a streaming APNG writer.  The writer is destroyed even
 * when finalization fails.
 *
 * @return 0 on success, -1 on output failure or a frame-count mismatch.
 */
ICV_EXPORT extern int icv_anim_writer_close(icv_anim_writer_t *writer);

__END_DECLS

#endif /* ICV_ANIM_H */

/** @} */
/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
