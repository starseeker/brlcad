/*                         F B _ I O . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libimgstream/fb_io.c
 *
 * Neutral file/image transfer helpers for imgstream framebuffers.
 */

#include "common.h"

#include "bu/str.h"

#include <limits.h>
#include <string.h>

#include "bio.h"
#include "bu/malloc.h"
#include "icv.h"
#include "imgstream/fb_compat.h"

#define RGB_CHANNELS 3

static int
fb_dimensions(const imgstream_fb_t *fb, int *width, int *height)
{
    if (!fb || !width || !height)
	return -1;
    size_t w = imgstream_fb_width(fb);
    size_t h = imgstream_fb_height(fb);
    if (!w || !h || w > INT_MAX || h > INT_MAX)
	return -1;
    *width = (int)w;
    *height = (int)h;
    return 0;
}

static int
fb_set_import_view(imgstream_fb_t *fb, int width, int height,
	int xoff, int yoff, int image_width, int image_height,
	int inverse, int zoom)
{
    if (!zoom)
	return imgstream_fb_view(fb, width / 2, height / 2, 1, 1);
    if (image_width <= 0 || image_height <= 0)
	return -1;

    int factor = width / image_width;
    int yfactor = height / image_height;
    if (yfactor < factor)
	factor = yfactor;
    if (factor < 1)
	factor = 1;
    int ycenter = inverse ? height - 1 - (yoff + image_height / 2) :
	yoff + image_height / 2;
    return imgstream_fb_view(fb, xoff + image_width / 2, ycenter,
	factor, factor);
}

int
imgstream_fb_export_pix_fp(const imgstream_fb_t *fb, FILE *fp,
	size_t requested_width, size_t requested_height, int crunch, int inverse)
{
    if (!fb || !fp)
	return -1;

    size_t width = imgstream_fb_width(fb);
    size_t height = imgstream_fb_height(fb);
    if (requested_width && requested_width < width)
	width = requested_width;
    if (requested_height && requested_height < height)
	height = requested_height;
    if (!width || !height || width > ((size_t)-1) / RGB_CHANNELS)
	return -1;

    struct imgstream_fb_colormap cmap;
    if (crunch && (imgstream_fb_rmap(fb, &cmap) != 0 ||
	imgstream_fb_colormap_is_linear(&cmap)))
	crunch = 0;

    size_t row_bytes = width * RGB_CHANNELS;
    unsigned char *row = (unsigned char *)bu_malloc(row_bytes,
	"imgstream PIX export row");
    int ret = 0;
    for (size_t line = 0; line < height; line++) {
	size_t y = inverse ? height - 1 - line : line;
	if (imgstream_fb_read(fb, 0, (int)y, row, width) != (ssize_t)width) {
	    ret = -1;
	    break;
	}
	if (crunch) {
	    for (size_t x = 0; x < width; x++) {
		row[x * RGB_CHANNELS + 0] =
		    (unsigned char)(cmap.red[row[x * RGB_CHANNELS + 0]] >> 8);
		row[x * RGB_CHANNELS + 1] =
		    (unsigned char)(cmap.green[row[x * RGB_CHANNELS + 1]] >> 8);
		row[x * RGB_CHANNELS + 2] =
		    (unsigned char)(cmap.blue[row[x * RGB_CHANNELS + 2]] >> 8);
	    }
	}
	if (fwrite(row, 1, row_bytes, fp) != row_bytes) {
	    ret = -1;
	    break;
	}
    }
    bu_free(row, "imgstream PIX export row");
    return ret;
}

static int
fb_skip_fd(int fd, size_t bytes, unsigned char *scratch, size_t scratch_size)
{
    while (bytes) {
	size_t chunk = bytes < scratch_size ? bytes : scratch_size;
	ssize_t got = read(fd, scratch, chunk);
	if (got <= 0)
	    return -1;
	bytes -= (size_t)got;
    }
    return 0;
}

int
imgstream_fb_import_pix_fd(imgstream_fb_t *fb, int fd,
	const char *filename, size_t file_width, size_t file_height, int autosize,
	const struct imgstream_fb_import_options *options)
{
    struct imgstream_fb_import_options defaults =
	IMGSTREAM_FB_IMPORT_OPTIONS_INIT;
    const struct imgstream_fb_import_options *opts = options ? options :
	&defaults;
    int screen_width = 0;
    int screen_height = 0;
    if (!fb || fd < 0 || fb_dimensions(fb, &screen_width, &screen_height) != 0)
	return -1;

    if (autosize && filename && filename[0] && bu_strcmp(filename, "-") != 0) {
	size_t width = 0;
	size_t height = 0;
	if (imgstream_image_file_size(&width, &height, filename,
		RGB_CHANNELS)) {
	    file_width = width;
	    file_height = height;
	}
    }
    if (!file_width || !file_height || file_width > INT_MAX ||
	file_height > INT_MAX || file_width > ((size_t)-1) / RGB_CHANNELS)
	return -1;

    const int file_xoff = opts->file_xoff > 0 ? opts->file_xoff : 0;
    const int file_yoff = opts->file_yoff > 0 ? opts->file_yoff : 0;
    int source_x = file_xoff;
    int destination_x = opts->screen_xoff;
    if (destination_x < 0) {
	source_x -= destination_x;
	destination_x = 0;
    }
    if ((size_t)source_x >= file_width || destination_x >= screen_width ||
	(size_t)file_yoff >= file_height)
	return 0;

    size_t output_width = file_width - (size_t)source_x;
    if (output_width > (size_t)(screen_width - destination_x))
	output_width = (size_t)(screen_width - destination_x);
    size_t output_height = file_height - (size_t)file_yoff;
    size_t row_bytes = file_width * RGB_CHANNELS;
    unsigned char *row = (unsigned char *)bu_malloc(row_bytes,
	"imgstream PIX import row");

    if (opts->clear) {
	const unsigned char black[3] = {0, 0, 0};
	if (imgstream_fb_clear(fb, black) != 0) {
	    bu_free(row, "imgstream PIX import row");
	    return -1;
	}
    }
    if (file_yoff && fb_skip_fd(fd, (size_t)file_yoff * row_bytes, row,
	row_bytes) != 0) {
	bu_free(row, "imgstream PIX import row");
	return -1;
    }

    int ret = 0;
    for (size_t source_y = 0; source_y < output_height; source_y++) {
	ssize_t got = bu_mread(fd, (char *)row, row_bytes);
	if (got == 0)
	    break;
	if (got < 0 || (size_t)got < row_bytes) {
	    ret = -1;
	    break;
	}
	int logical_y = opts->screen_yoff + (int)source_y;
	int destination_y = opts->inverse ? screen_height - 1 - logical_y :
	    logical_y;
	if (destination_y < 0 || destination_y >= screen_height)
	    continue;
	if (imgstream_fb_write(fb, destination_x, destination_y,
		row + (size_t)source_x * RGB_CHANNELS,
		output_width) != (ssize_t)output_width) {
	    ret = -1;
	    break;
	}
    }
    bu_free(row, "imgstream PIX import row");
    if (ret != 0)
	return ret;
    return fb_set_import_view(fb, screen_width, screen_height,
	destination_x, opts->screen_yoff, (int)output_width,
	(int)output_height, opts->inverse, opts->zoom);
}

int
imgstream_fb_import_icv(imgstream_fb_t *fb, const struct icv_image *image,
	const struct imgstream_fb_import_options *options)
{
    struct imgstream_fb_import_options defaults =
	IMGSTREAM_FB_IMPORT_OPTIONS_INIT;
    const struct imgstream_fb_import_options *opts = options ? options :
	&defaults;
    int screen_width = 0;
    int screen_height = 0;
    if (!fb || !image || fb_dimensions(fb, &screen_width, &screen_height) != 0)
	return -1;
    if (!ICV_IMAGE_IS_INITIALIZED(image) || !image->data ||
	!image->width || !image->height || !image->channels ||
	image->width > INT_MAX || image->height > INT_MAX)
	return -1;

    int source_x = opts->file_xoff > 0 ? opts->file_xoff : 0;
    int source_y = opts->file_yoff > 0 ? opts->file_yoff : 0;
    int destination_x = opts->screen_xoff;
    if (destination_x < 0) {
	source_x -= destination_x;
	destination_x = 0;
    }
    if ((size_t)source_x >= image->width || (size_t)source_y >= image->height ||
	destination_x >= screen_width)
	return 0;

    size_t output_width = image->width - (size_t)source_x;
    if (output_width > (size_t)(screen_width - destination_x))
	output_width = (size_t)(screen_width - destination_x);
    size_t output_height = image->height - (size_t)source_y;
    unsigned char *pixels = icv_data2uchar(image);
    if (!pixels)
	return -1;
    unsigned char *row = (unsigned char *)bu_malloc(output_width *
	RGB_CHANNELS, "imgstream ICV import row");

    if (opts->clear) {
	const unsigned char black[3] = {0, 0, 0};
	if (imgstream_fb_clear(fb, black) != 0) {
	    bu_free(row, "imgstream ICV import row");
	    bu_free(pixels, "imgstream ICV pixels");
	    return -1;
	}
    }

    int ret = 0;
    for (size_t line = 0; line < output_height; line++) {
	int logical_y = opts->screen_yoff + (int)line;
	int destination_y = opts->inverse ? screen_height - 1 - logical_y :
	    logical_y;
	if (destination_y < 0 || destination_y >= screen_height)
	    continue;
	const unsigned char *source = pixels +
	    (((size_t)source_y + line) * image->width + (size_t)source_x) *
	    image->channels;
	for (size_t x = 0; x < output_width; x++) {
	    if (image->channels == 1) {
		row[x * RGB_CHANNELS + 0] = source[x];
		row[x * RGB_CHANNELS + 1] = source[x];
		row[x * RGB_CHANNELS + 2] = source[x];
	    } else {
		row[x * RGB_CHANNELS + 0] = source[x * image->channels + 0];
		row[x * RGB_CHANNELS + 1] = source[x * image->channels + 1];
		row[x * RGB_CHANNELS + 2] = source[x * image->channels + 2];
	    }
	}
	if (imgstream_fb_write(fb, destination_x, destination_y, row,
		output_width) != (ssize_t)output_width) {
	    ret = -1;
	    break;
	}
    }

    bu_free(row, "imgstream ICV import row");
    bu_free(pixels, "imgstream ICV pixels");
    if (ret != 0)
	return ret;
    return fb_set_import_view(fb, screen_width, screen_height,
	destination_x, opts->screen_yoff, (int)output_width,
	(int)output_height, opts->inverse, opts->zoom);
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
