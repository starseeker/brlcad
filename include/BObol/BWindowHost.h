/*                 B W I N D O W H O S T . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BWindowHost.h */

#ifndef BOBOL_BWINDOWHOST_H
#define BOBOL_BWINDOWHOST_H

#include "BObol/BDefines.h"
#include "BObol/BFramebuffer.h"
#include "BObol/BInput.h"

#include <Inventor/SbBasic.h>
#include <Inventor/SbString.h>

#include <stddef.h>

class BObolViewController;
class BObolFramebufferStream;
class SoBRLImageSource;
class SoBRLViewportImage;
struct BObolFramebufferBridge;

struct imgstream_fb;
typedef struct imgstream_fb imgstream_fb_t;
struct imgstream_fb_spec_info_s;
typedef struct imgstream_fb_spec_info_s imgstream_fb_spec_info_t;
struct imgstream_fb_view_s;
typedef struct imgstream_fb_view_s imgstream_fb_view_t;
struct imgstream_fb_cursor_s;
typedef struct imgstream_fb_cursor_s imgstream_fb_cursor_t;

enum BObolWindowMode {
    BOBOL_WINDOW_TOPLEVEL = 0,
    BOBOL_WINDOW_EMBEDDED = 1,
    BOBOL_WINDOW_HEADLESS = 2,
    BOBOL_WINDOW_DIAGNOSTIC = 3
};

enum BObolWindowBackend {
    BOBOL_WINDOW_BACKEND_AUTO = 0,
    BOBOL_WINDOW_BACKEND_QT = 1,
    BOBOL_WINDOW_BACKEND_TK = 2,
    BOBOL_WINDOW_BACKEND_OPENGL = 3,
    BOBOL_WINDOW_BACKEND_OFFSCREEN = 4,
    BOBOL_WINDOW_BACKEND_DIAGNOSTIC = 5
};

struct BOBOL_EXPORT BObolWindowDesc {
    BObolWindowMode mode;
    BObolWindowBackend backend;
    unsigned int width;
    unsigned int height;
    SbString title;
    SbString display;
    SbString nativeIdHint;
    SbBool visible;
};

struct BObolWindowHostPrivate;

class BOBOL_EXPORT BObolWindowHost {
public:
    BObolWindowHost(void);
    virtual ~BObolWindowHost(void);

    virtual int open(const BObolWindowDesc *desc = NULL);
    virtual void close(void);
    SbBool isOpen(void) const;
    const BObolWindowDesc &getDesc(void) const;

    void attachController(BObolViewController *controller, SbBool takeOwnership = FALSE);
    BObolViewController *getController(void) const;

    virtual int handleInputEvent(const BObolInputEvent *event,
	const BObolInputProfile *profile = NULL);
    virtual int applyInputAction(BObolInputAction action,
	const BObolInputEvent *event = NULL);
    virtual int poll(const BObolInputProfile *profile = NULL);
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
    virtual int setFramebufferComposition(imgstream_fb_t *fb,
	BObolFramebufferComposition composition);

    int getFramebufferCount(void) const;
    SoBRLImageSource *getFramebufferImageSource(imgstream_fb_t *fb) const;
    SoBRLViewportImage *getFramebufferViewportImage(imgstream_fb_t *fb) const;

private:
    friend class BObolFramebufferStream;
    friend struct BObolFramebufferBridge;

    void registerFramebufferStream(BObolFramebufferStream *stream);
    void unregisterFramebufferStream(BObolFramebufferStream *stream);
    void registerDisplayFramebuffer(imgstream_fb_t *fb);
    void unregisterDisplayFramebuffer(imgstream_fb_t *fb);
    void detachFramebufferStreams(void);
    void detachDisplayFramebuffers(void);

    BObolWindowHostPrivate *p;
};

BOBOL_EXPORT imgstream_fb_t *bobol_window_host_open_display_framebuffer(
    BObolWindowHost *host, const char *spec, size_t width, size_t height);

#endif /* BOBOL_BWINDOWHOST_H */
