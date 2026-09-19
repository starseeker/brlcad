/*                    P R O G R E S S I V E . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "rtuif.h"

static rt_framebuffer_flush_callback flush_callback = NULL;
static void *flush_callback_data = NULL;
static rt_frame_runner_callback frame_runner_callback = NULL;
static void *frame_runner_data = NULL;

void
rt_frame_runner_set(rt_frame_runner_callback callback, void *data)
{
    frame_runner_callback = callback;
    frame_runner_data = callback ? data : NULL;
}

int
rt_frame_runner_run(rt_frame_execute_callback execute, int framenumber,
	void *data)
{
    if (!execute)
	return -1;
    if (frame_runner_callback)
	return frame_runner_callback(execute, framenumber, frame_runner_data);
    return execute(framenumber, data);
}

void
rt_framebuffer_flush_callback_set(rt_framebuffer_flush_callback callback,
	void *data)
{
    flush_callback = callback;
    flush_callback_data = callback ? data : NULL;
}

void
rt_fb_progressive_flush(void)
{
    if (flush_callback)
	flush_callback(flush_callback_data);
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
