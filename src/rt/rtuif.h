/*                         R T U I F . H
 * BRL-CAD
 *
 * Copyright (c) 1985-2026 United States Government as represented by
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
/** @file rt/rtuif.h
 *
 * Private API header for developing applications that utilize rt's
 * common Ray Trace User Interface Framework (RTUIF).
 *
 * While not (yet) public API, it's a usable interface for rapidly
 * developing applications involving grids of rays.  It provides the
 * same command-line interface as 'rt' and is the foundation for
 * numerous BRL-CAD ray tracing based applications such as rtweight,
 * rtarea, and others.  See viewdummy.c for an example template.
 *
 */

#ifndef RT_RTUIF_H
#define RT_RTUIF_H

#include "common.h"

#include <stdio.h>

#include "raytrace.h"
#include "optical.h"

typedef int (*rt_frame_execute_callback)(int framenumber, void *data);
typedef int (*rt_frame_runner_callback)(rt_frame_execute_callback execute,
	int framenumber, void *data);
typedef void (*rt_framebuffer_flush_callback)(void *data);

/* Optional application-owned frame runner.  rt uses this only when a visible
 * display session needs the caller thread to pump its native event loop. */
extern void rt_frame_runner_set(rt_frame_runner_callback callback, void *data);
/* Invoke one frame through the installed runner, or directly when no display
 * session owns the frame.  This is internal RTUIF plumbing shared by do.c and
 * the application setup code. */
extern int rt_frame_runner_run(rt_frame_execute_callback execute,
	int framenumber, void *data);

/* The common view code requests a progressive flush after framebuffer writes.
 * Headless consumers retain the default no-op; a visible host registers its
 * owner-thread-safe presentation callback. */
extern void rt_framebuffer_flush_callback_set(
	rt_framebuffer_flush_callback callback, void *data);
extern void rt_fb_progressive_flush(void);


/**
 * Called by main() at the start of a run.  Should set up rayhit() and
 * raymiss() in the application structure.  The rayhit() callback
 * function is called via the a_hit linkage from rt_shootray() when a
 * ray hits.  The raymiss() callback function is called via the a_miss
 * linkage from rt_shootray() when a ray misses.
 *
 * Returns 1 if framebuffer should be opened, else 0.
 */
extern int view_init(struct application *ap, char *file, char *obj, int minus_o, int minus_F);

/**
 * Called by do_prep(), just before rt_prep() is called, in
 * “do.c”. This allows the lighting model to get set up for this
 * frame, e.g., generate lights, associate materials routines, etc.
 */
extern void view_setup(struct rt_i *rtip);

/**
 * Called at the beginning of a frame. Called by do_frame() just
 * before raytracing starts.
 */
extern void view_2init(struct application *ap, char *framename);

/**
 * Called by worker() after the end of processing for each pixel.
 */
extern void view_pixel(struct application *ap);

/**
 * Called after the end of each ray trace scanline.
 */
extern void view_eol(struct application *ap);

/**
 * Called in do_frame() at the end of a frame, just after raytracing
 * completes.
 */
extern void view_end(struct application *ap);

/**
 * Called before rt_clean() in do.c, for releasing resources allocated
 * and opened during ray tracing.
 */
extern void view_cleanup(struct rt_i *rtip);

#endif  /* RT_RTUIF_H */
/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
