/*                         U I _ F B . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file ui_fb.c
 *
 * Windowed front end for the BRL-CAD launcher.  Menu pixels use imgstream;
 * a display session owns presentation and normalized application input.
 */

#include "common.h"

#include <stdio.h>
#include <string.h>

#include "bu.h"
#include "bu/color.h"
#include "bu/mime.h"
#include "bu/snooze.h"
#include "icv.h"
#include "imgstream/fb_compat.h"
#include "BObol/BDisplaySession.h"
#if defined(HAVE_QTCAD_OBOL_DISPLAY_PROVIDER)
#  include "qtcad/display_provider.h"
#elif defined(HAVE_TCLCAD_OBOL_DISPLAY_PROVIDER)
#  include "tclcad/setup.h"
#endif

#include "launcher.h"
#include "fbtext.h"

#define QUIT_ACTION (-1)

/* Palette (R,G,B). */
static const unsigned char col_bg[3]        = { 26, 28, 34 };
static const unsigned char col_panel[3]     = { 52, 58, 70 };
static const unsigned char col_panel_off[3] = { 38, 40, 46 };
static const unsigned char col_hover[3]     = { 178, 66, 42 };
static const unsigned char col_border[3]    = { 84, 92, 108 };
static const unsigned char col_accent[3]    = { 200, 80, 50 };
static const unsigned char col_text[3]      = { 235, 238, 242 };
static const unsigned char col_text_dim[3]  = { 150, 156, 168 };
static const unsigned char col_text_off[3]  = { 96, 100, 110 };

struct button {
    int x, yb, w, h;	/* framebuffer coords (origin lower-left) */
    int action;		/* registry index, or QUIT_ACTION */
};


static void
fill_rect(imgstream_fb_t *fbp, int x, int yb, int w, int h, const unsigned char c[3])
{
    int rW = imgstream_fb_width(fbp);
    int rH = imgstream_fb_height(fbp);
    unsigned char *row;
    int i, yy;

    if (x < 0) { w += x; x = 0; }
    if (yb < 0) { h += yb; yb = 0; }
    if (x + w > rW) w = rW - x;
    if (yb + h > rH) h = rH - yb;
    if (w <= 0 || h <= 0)
	return;

    row = (unsigned char *)bu_malloc(w * 3, "fill_rect row");
    for (i = 0; i < w; i++) {
	row[i * 3 + RED] = c[RED];
	row[i * 3 + GRN] = c[GRN];
	row[i * 3 + BLU] = c[BLU];
    }
    for (yy = yb; yy < yb + h; yy++)
	(void)imgstream_fb_write(fbp, x, yy, row, w);
    bu_free(row, "fill_rect row");
}


static void
draw_border(imgstream_fb_t *fbp, int x, int yb, int w, int h, const unsigned char c[3])
{
    fill_rect(fbp, x, yb, w, 1, c);
    fill_rect(fbp, x, yb + h - 1, w, 1, c);
    fill_rect(fbp, x, yb, 1, h, c);
    fill_rect(fbp, x + w - 1, yb, 1, h, c);
}


/*
 * Blit an icv image upright with its lower-left corner at (sx, sy) in
 * framebuffer coordinates.  libicv's uchar row buffer is bottom-up (row 0 is
 * the bottom of the image), matching imgstream's bottom-up rows, so the copy is a
 * straight 1:1 mapping.  Doing this directly (rather than via imgstream_fb_read_icv)
 * keeps the orientation consistent with the rest of the imgstream_fb_write drawing below.
 */
static void
blit_icv(imgstream_fb_t *fbp, icv_image_t *img, int sx, int sy)
{
    unsigned char *rgb;
    int iw, ih, i;
    int rH = imgstream_fb_height(fbp);

    if (!img)
	return;
    rgb = icv_data2uchar(img);
    if (!rgb)
	return;
    iw = (int)img->width;
    ih = (int)img->height;
    for (i = 0; i < ih; i++) {
	int fy = sy + i;
	if (fy < 0 || fy >= rH)
	    continue;
	(void)imgstream_fb_write(fbp, sx, fy, rgb + (size_t)i * iw * 3, iw);
    }
    bu_free(rgb, "blit_icv rgb");
}


/*
 * Load and scale the splash image.  Prefers a rendered splash asset, falling
 * back to the bundled BRL-CAD logo.  Returns NULL if neither is present.
 */
static icv_image_t *
load_splash(int maxw, int maxh)
{
    char path[MAXPATHLEN] = {0};
    icv_image_t *img = NULL;
    double sc;

    bu_dir(path, sizeof(path), BU_DIR_DATA, "launcher", "splash.png", NULL);
    if (bu_file_exists(path, NULL))
	img = icv_read(path, BU_MIME_IMAGE_AUTO, 0, 0);

    if (!img) {
	bu_dir(path, sizeof(path), BU_DIR_DATA, "images", "brlLogo-nobg.png", NULL);
	if (bu_file_exists(path, NULL))
	    img = icv_read(path, BU_MIME_IMAGE_AUTO, 0, 0);
    }

    if (!img)
	return NULL;

    sc = 1.0;
    if ((int)img->width > maxw)
	sc = (double)maxw / (double)img->width;
    if ((double)img->height * sc > (double)maxh)
	sc = (double)maxh / (double)img->height;

    if (sc < 1.0) {
	size_t nw = (size_t)((double)img->width * sc);
	size_t nh = (size_t)((double)img->height * sc);
	if (nw > 0 && nh > 0)
	    (void)icv_resize(img, ICV_RESIZE_BINTERP, nw, nh, 0);
    }

    return img;
}


static void
redraw(imgstream_fb_t *fbp, struct fbtext *ft, struct app_registry *r,
       icv_image_t *splash, struct button *btns, int nbtns, int hover)
{
    int rW = imgstream_fb_width(fbp);
    int rH = imgstream_fb_height(fbp);
    int i;
    int lh = fbtext_line_height(ft);
    const int text_padding = 18;
    int description_offset = fbtext_string_width(ft, "Quit") + 2 * text_padding;
    for (i = 0; i < r->count; i++) {
	const int name_extent = fbtext_string_width(ft, r->apps[i].name) + 2 * text_padding;
	if (name_extent > description_offset)
	    description_offset = name_extent;
    }

    fill_rect(fbp, 0, 0, rW, rH, col_bg);

    /* Splash across the top, centered. */
    if (splash) {
	int iw = (int)splash->width;
	int ih = (int)splash->height;
	int sx = (rW - iw) / 2;
	int sy = rH - 24 - ih;		/* 24px top margin, image top-aligned */
	if (sx < 0) sx = 0;
	if (sy < 0) sy = 0;
	blit_icv(fbp, splash, sx, sy);
	/* Accent rule under the splash. */
	fill_rect(fbp, 40, sy - 12, rW - 80, 2, col_accent);
    }

    /* Menu buttons. */
    for (i = 0; i < nbtns; i++) {
	struct button *b = &btns[i];
	int avail = 1;
	const char *name = "Quit";
	const char *desc = "Exit the launcher";
	const unsigned char *tcol = col_text;
	const unsigned char *dcol = col_text_dim;

	if (b->action != QUIT_ACTION) {
	    struct app_entry *e = &r->apps[b->action];
	    name = e->name;
	    desc = e->description;
	    avail = e->available;
	}

	if (i == hover)
	    fill_rect(fbp, b->x, b->yb, b->w, b->h, col_hover);
	else
	    fill_rect(fbp, b->x, b->yb, b->w, b->h, avail ? col_panel : col_panel_off);
	draw_border(fbp, b->x, b->yb, b->w, b->h, col_border);

	if (!avail) {
	    tcol = col_text_off;
	    dcol = col_text_off;
	}

	/* Name at the left, description further right, vertically centered. */
	fbtext_draw(fbp, ft, b->x + 18, b->yb + (b->h - lh) / 2, name, tcol);
	if (desc && desc[0]) {
	    struct bu_vls label = BU_VLS_INIT_ZERO;
	    if (b->action != QUIT_ACTION && !avail)
		bu_vls_sprintf(&label, "%s  (not installed)", desc);
	    else
		bu_vls_sprintf(&label, "%s", desc);
	    fbtext_draw(fbp, ft, b->x + description_offset, b->yb + (b->h - lh) / 2, bu_vls_cstr(&label), dcol);
	    bu_vls_free(&label);
	}
    }

    /* Footer hint. */
    fbtext_draw(fbp, ft, 44, 18,
		"Click an entry.  q / Esc: Quit",
		col_text_dim);

    imgstream_fb_flush(fbp);
}


static int
hit_test(struct button *btns, int nbtns, int x, int y)
{
    int i;
    for (i = 0; i < nbtns; i++) {
	struct button *b = &btns[i];
	if (x >= b->x && x < b->x + b->w && y >= b->yb && y < b->yb + b->h)
	    return i;
    }
    return -1;
}


/* Perform the action for a button; sets *running to 0 to exit. */
static void
activate(struct app_registry *r, struct button *b, int *running)
{
    if (b->action == QUIT_ACTION) {
	*running = 0;
	return;
    }
    {
	struct app_entry *e = &r->apps[b->action];
	if (e->available)
	    (void)app_launch(e);
    }
}


enum { LAUNCHER_MENU_INPUT = 1, LAUNCHER_PRIMARY_BUTTON = 0 };

struct menu_input {
    struct app_registry *registry;
    struct button *buttons;
    int count;
    int hover;
    int running;
    int redraw;
    int viewport_changed;
    int image_width;
    int image_height;
    unsigned int viewport_width;
    unsigned int viewport_height;
};

static int
menu_input_event(void *data, BObolInputAction UNUSED(action), const BObolInputEvent *event)
{
    struct menu_input *menu = (struct menu_input *)data;
    int hit;
    switch (event->type) {
	case BOBOL_INPUT_CLOSE:
	    menu->running = 0;
	    break;
	case BOBOL_INPUT_POINTER_MOTION:
	case BOBOL_INPUT_POINTER_RELEASE:
	    if (!menu->viewport_width || !menu->viewport_height)
		return BOBOL_INPUT_RESULT_UNHANDLED;
	    hit = hit_test(menu->buttons, menu->count,
		(int)((double)event->x * menu->image_width / menu->viewport_width),
		menu->image_height - 1 -
		(int)((double)event->y * menu->image_height / menu->viewport_height));
	    if (event->type == BOBOL_INPUT_POINTER_MOTION && hit != menu->hover) {
		menu->hover = hit;
		menu->redraw = 1;
	    }
	    if (event->type == BOBOL_INPUT_POINTER_RELEASE && event->button == LAUNCHER_PRIMARY_BUTTON && hit >= 0) {
		activate(menu->registry, &menu->buttons[hit], &menu->running);
		menu->redraw = 1;
	    }
	    break;
	case BOBOL_INPUT_KEY_PRESS:
	    if (event->key == 'Q' || event->key == 'q' || event->key == 27) {
		menu->running = 0;
	    } else if (event->key >= '1' && event->key <= '9' &&
		event->key - '1' < menu->count) {
		activate(menu->registry, &menu->buttons[event->key - '1'], &menu->running);
		menu->redraw = 1;
	    } else if ((event->key == '\r' || event->key == '\n') && menu->hover >= 0) {
		activate(menu->registry, &menu->buttons[menu->hover], &menu->running);
		menu->redraw = 1;
	    }
	    break;
	case BOBOL_INPUT_RESIZE:
	    menu->viewport_width = event->width;
	    menu->viewport_height = event->height;
	    menu->viewport_changed = 1;
	    menu->redraw = 1;
	    break;
	case BOBOL_INPUT_EXPOSE:
	    menu->redraw = 1;
	    break;
	default:
	    return BOBOL_INPUT_RESULT_UNHANDLED;
    }
    return BOBOL_INPUT_RESULT_HANDLED;
}

int
ui_fb_run(struct app_registry *r)
{
    imgstream_fb_t *fbp;
    struct fbtext ft;
    icv_image_t *splash;
    struct button *btns;
    int nbtns;
    int rW, rH;
    int i;
    bobol_display_session_t *session;
    long rate;

    const int win_w = 820;
    const int win_h = 780;
    int bw;
    const int bh = 48;
    const int gap = 10;

    if (!r || r->count == 0)
	return -1;

#if defined(HAVE_QTCAD_OBOL_DISPLAY_PROVIDER)
    if (!qtcad_obol_display_provider_register())
	return -1;
#elif defined(HAVE_TCLCAD_OBOL_DISPLAY_PROVIDER)
    if (!tclcad_obol_display_provider_register())
	return -1;
#else
    return -1;
#endif
    session = bobol_display_session_open("/dev/swrast", win_w, win_h, "BRL-CAD");
    if (!session)
	return -1;
    fbp = bobol_display_session_framebuffer(session);

    rW = imgstream_fb_width(fbp);
    rH = imgstream_fb_height(fbp);

    if (fbtext_open(&ft, NULL) != 0) {
	/* No font available -- not usable as a graphical menu. */
	bobol_display_session_close(session);
	return -1;
    }


    splash = load_splash(rW - 40, 300);

    /* Build the button list: one per registry entry, plus Quit. */
    bw = rW - 60;
    nbtns = r->count + 1;
    btns = (struct button *)bu_calloc(nbtns, sizeof(struct button), "buttons");
    {
	int splash_h = splash ? (int)splash->height : 0;
	int menu_top = 24 + splash_h + 28;	/* top-down y of first button */
	int bx = (rW - bw) / 2;
	for (i = 0; i < nbtns; i++) {
	    int ty = menu_top + i * (bh + gap);
	    btns[i].x = bx;
	    btns[i].w = bw;
	    btns[i].h = bh;
	    btns[i].yb = rH - (ty + bh);	/* convert top-down to fb bottom */
	    btns[i].action = (i < r->count) ? i : QUIT_ACTION;
	}
    }

    struct menu_input menu = {r, btns, nbtns, -1, 1, 1, 1, rW, rH,
	(unsigned int)win_w, (unsigned int)win_h};
    const BObolInputBinding bindings[] = {
	{BOBOL_INPUT_POINTER_MOTION, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY, 0, 0, 0, LAUNCHER_MENU_INPUT},
	{BOBOL_INPUT_POINTER_RELEASE, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY, 0, 0, 0, LAUNCHER_MENU_INPUT},
	{BOBOL_INPUT_KEY_PRESS, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY, 0, 0, 0, LAUNCHER_MENU_INPUT},
	{BOBOL_INPUT_RESIZE, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY, 0, 0, 0, LAUNCHER_MENU_INPUT},
	{BOBOL_INPUT_CLOSE, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY, 0, 0, 0, LAUNCHER_MENU_INPUT},
	{BOBOL_INPUT_EXPOSE, BOBOL_INPUT_ANY, BOBOL_INPUT_ANY, 0, 0, 0, LAUNCHER_MENU_INPUT}
    };
    const BObolInputActionLayer layer = {"launcher", bindings,
	sizeof(bindings) / sizeof(bindings[0]), menu_input_event};
    bobol_display_endpoint_t *endpoint = bobol_display_session_endpoint(session);
    if (!bobol_display_endpoint_input_action_layer_set(endpoint, &layer, &menu, &menu)) {
	if (splash)
	    icv_destroy(splash);
	bu_free(btns, "buttons");
	fbtext_close(&ft);
	bobol_display_session_close(session);
	return -1;
    }
    rate = bobol_display_session_poll_rate(session);
    int status = 0;
    while (menu.running) {
	if (menu.viewport_changed && menu.viewport_width && menu.viewport_height) {
	    if (imgstream_fb_viewport(fbp, 0, 0, (int)menu.viewport_width,
		    (int)menu.viewport_height) != 0) {
		bu_log("brlcad-launcher: unable to resize menu presentation\n");
		status = -1;
		break;
	    }
	    menu.viewport_changed = 0;
	}
	if (menu.redraw) {
	    redraw(fbp, &ft, r, splash, btns, nbtns, menu.hover);
	    menu.redraw = 0;
	}
	const int poll_status = bobol_display_session_poll(session);
	if (poll_status != 0) {
	    if (poll_status < 0) {
		bu_log("brlcad-launcher: display event processing failed\n");
		status = -1;
	    }
	    break;
	}
	if (rate > 0)
	    bu_snooze(rate);
    }
    (void)bobol_display_endpoint_input_action_layer_clear_if(endpoint, &menu);

    if (splash)
	icv_destroy(splash);
    bu_free(btns, "buttons");
    fbtext_close(&ft);
    bobol_display_session_close(session);
    return status;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
