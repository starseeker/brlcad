/*                  W I N D O W _ H O S T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/window_host.h */

#ifndef BRLOBOL_WINDOW_HOST_H
#define BRLOBOL_WINDOW_HOST_H

#include "brlobol/defines.h"
#include "brlobol/input.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbString.h>

#include <stddef.h>

class BRLObolViewController;
class SoBRLImageSource;
class SoBRLViewportImage;

struct imgstream_fb;
typedef struct imgstream_fb imgstream_fb_t;
struct imgstream_fb_spec_info_s;
typedef struct imgstream_fb_spec_info_s imgstream_fb_spec_info_t;
struct imgstream_fb_view_s;
typedef struct imgstream_fb_view_s imgstream_fb_view_t;
struct imgstream_fb_cursor_s;
typedef struct imgstream_fb_cursor_s imgstream_fb_cursor_t;

enum BRLObolWindowMode {
    BRLOBOL_WINDOW_TOPLEVEL = 0,
    BRLOBOL_WINDOW_EMBEDDED = 1,
    BRLOBOL_WINDOW_HEADLESS = 2,
    BRLOBOL_WINDOW_DIAGNOSTIC = 3
};

enum BRLObolWindowBackend {
    BRLOBOL_WINDOW_BACKEND_AUTO = 0,
    BRLOBOL_WINDOW_BACKEND_QT = 1,
    BRLOBOL_WINDOW_BACKEND_TK = 2,
    BRLOBOL_WINDOW_BACKEND_OPENGL = 3,
    BRLOBOL_WINDOW_BACKEND_OFFSCREEN = 4,
    BRLOBOL_WINDOW_BACKEND_DIAGNOSTIC = 5
};

struct BRLOBOL_EXPORT BRLObolWindowDesc {
    BRLObolWindowMode mode;
    BRLObolWindowBackend backend;
    unsigned int width;
    unsigned int height;
    SbString title;
    SbString display;
    SbString nativeIdHint;
    SbBool visible;
};

struct BRLObolWindowHostPrivate;

class BRLOBOL_EXPORT BRLObolWindowHost {
public:
    BRLObolWindowHost(void);
    virtual ~BRLObolWindowHost(void);

    virtual int open(const BRLObolWindowDesc *desc = NULL);
    virtual void close(void);
    SbBool isOpen(void) const;
    const BRLObolWindowDesc &getDesc(void) const;

    void attachController(BRLObolViewController *controller, SbBool takeOwnership = FALSE);
    BRLObolViewController *getController(void) const;

    virtual int handleInputEvent(const BRLObolInputEvent *event,
	const BRLObolInputProfile *profile = NULL);
    virtual int applyInputAction(BRLObolInputAction action,
	const BRLObolInputEvent *event = NULL);
    virtual int poll(const BRLObolInputProfile *profile = NULL);
    virtual long pollRate(void) const;

    virtual int openFramebuffer(imgstream_fb_t *fb,
	const imgstream_fb_spec_info_t *info);
    virtual void closeFramebuffer(imgstream_fb_t *fb);
    virtual int flushFramebuffer(imgstream_fb_t *fb);
    virtual int resetFramebuffer(imgstream_fb_t *fb);
    virtual int setFramebufferViewport(imgstream_fb_t *fb,
	int left, int top, int right, int bottom);
    virtual int setFramebufferView(imgstream_fb_t *fb,
	const imgstream_fb_view_t *view);
    virtual int setFramebufferCursor(imgstream_fb_t *fb,
	const imgstream_fb_cursor_t *cursor);
    virtual int setFramebufferScreenCursor(imgstream_fb_t *fb,
	int mode, int x, int y);
    virtual int setFramebufferCursorShape(imgstream_fb_t *fb,
	const unsigned char *bits, int xbits, int ybits, int xorig, int yorig);

    int getFramebufferCount(void) const;
    SoBRLImageSource *getFramebufferImageSource(imgstream_fb_t *fb) const;
    SoBRLViewportImage *getFramebufferViewportImage(imgstream_fb_t *fb) const;

private:
    BRLObolWindowHostPrivate *p;
};

BRLOBOL_EXPORT imgstream_fb_t *brlobol_window_host_open_display_framebuffer(
    BRLObolWindowHost *host, const char *spec, size_t width, size_t height);

#endif /* BRLOBOL_WINDOW_HOST_H */
