/*
 * Togl - a Tk OpenGL widget
 *
 * Copyright 1996-2002  Brian Paul and Ben Bederson
 * Copyright 2005-2009  Greg Couch
 * See the LICENSE file for copyright details.
 *
 * This is an enhanced version of Togl called tcl3dTogl.
 * It adds functionality, so that Tcl-wrapped OpenGL functions
 * can be evaluated by a Togl widget.
 * It is part of the Tcl3D package.
 * The modifications are based on Togl Version 1.7.
 * Changes and additions are marked with TCL3D and are
 * Copyright 2005-2025 Paul Obermeier
 */

/*
 * Currently we support X11, Win32 and Mac OS X only
 */

#include <tcl.h>
#include <tk.h>

#ifdef TOGL_WGL
    #if TK_MAJOR_VERSION >= 9 || (TK_MAJOR_VERSION == 8 && TK_MINOR_VERSION >= 6)
        #ifndef UNICODE
            #define UNICODE
        #endif
        #ifndef _UNICODE
            #define _UNICODE
        #endif
        #define TOGL_CLASS_NAME TEXT("Togl Class")
        #ifndef TK_WIN_CHILD_CLASS_NAME
            #define TK_WIN_CHILD_CLASS_NAME TEXT("TkChild")
        #endif
    #else
        #define TOGL_CLASS_NAME "Togl Class"
        #ifndef TK_WIN_CHILD_CLASS_NAME
            #define TK_WIN_CHILD_CLASS_NAME "TkChild"
        #endif
    #endif
#endif

#include "togl.h"
#include "toglFont.h"

#if !defined(TOGL_AGL)
/* TkpMakeWindow is the one Tk-private entry point Togl still requires.  Keep
 * Tk's internal structures opaque rather than importing version-coupled
 * private headers. */
typedef struct TkWindow TkWindow;
extern Window TkpMakeWindow(TkWindow *winPtr, Window parent);
#endif

/* Check, if Tcl version supports Tcl_Size,
   which was introduced in Tcl 8.7 and 9.
*/
#ifndef TCL_SIZE_MAX
    #include <limits.h>
    #define TCL_SIZE_MAX INT_MAX

    #ifndef Tcl_Size
        typedef int Tcl_Size;
    #endif

    #define TCL_SIZE_MODIFIER ""
    #define Tcl_GetSizeIntFromObj Tcl_GetIntFromObj
#endif

/* Use WIDGREC to cast widgRec arguments */
#define WIDGREC (char *)

/*** Windows headers ***/
#if defined(TOGL_WGL)
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  undef WIN32_LEAN_AND_MEAN
#  include <winnt.h>

/* WGL extension constants used by dynamically loaded entry points.  Keeping
 * these local avoids a dependency on GLEW solely for declarations. */
#  ifndef WGL_DRAW_TO_WINDOW_ARB
#    define WGL_DRAW_TO_WINDOW_ARB 0x2001
#    define WGL_ACCELERATION_ARB 0x2003
#    define WGL_DOUBLE_BUFFER_ARB 0x2011
#    define WGL_RED_BITS_ARB 0x2015
#    define WGL_GREEN_BITS_ARB 0x2017
#    define WGL_BLUE_BITS_ARB 0x2019
#    define WGL_ALPHA_BITS_ARB 0x201B
#    define WGL_ACCUM_RED_BITS_ARB 0x201E
#    define WGL_ACCUM_GREEN_BITS_ARB 0x201F
#    define WGL_ACCUM_BLUE_BITS_ARB 0x2020
#    define WGL_ACCUM_ALPHA_BITS_ARB 0x2021
#    define WGL_DEPTH_BITS_ARB 0x2022
#    define WGL_STENCIL_BITS_ARB 0x2023
#    define WGL_AUX_BUFFERS_ARB 0x2024
#    define WGL_FULL_ACCELERATION_ARB 0x2027
#    define WGL_SAMPLE_BUFFERS_ARB 0x2041
#    define WGL_SAMPLES_ARB 0x2042
#    define WGL_CONTEXT_MAJOR_VERSION_ARB 0x2091
#    define WGL_CONTEXT_MINOR_VERSION_ARB 0x2092
#    define WGL_CONTEXT_FLAGS_ARB 0x2094
#    define WGL_CONTEXT_DEBUG_BIT_ARB 0x0001
#    define WGL_CONTEXT_FORWARD_COMPATIBLE_BIT_ARB 0x0002
#    define WGL_CONTEXT_PROFILE_MASK_ARB 0x9126
#    define WGL_CONTEXT_CORE_PROFILE_BIT_ARB 0x00000001
#    define WGL_CONTEXT_COMPATIBILITY_PROFILE_BIT_ARB 0x00000002
#  endif

#  undef near
#  undef far

/*** X Window System headers ***/
#elif defined(TOGL_X11)
#  include <X11/Xlib.h>
#  include <X11/Xutil.h>
#  include <X11/Xatom.h>        /* for XA_RGB_DEFAULT_MAP atom */
#  include <GL/glx.h>

/*** Mac headers ***/
#elif defined(TOGL_AGL)
#  define Cursor QDCursor
#  include <AGL/agl.h>
#  undef Cursor
#  include <tkInt.h>
#  include <tkMacOSX.h>
#  include <tkMacOSXInt.h>      /* usa MacDrawable */
#  include <ApplicationServices/ApplicationServices.h>
#  define Togl_MacOSXGetDrawablePort(togl) TkMacOSXGetDrawablePort((Drawable) ((TkWindow *) togl->TkWin)->privatePtr)


#else /* make sure only one platform defined */
#  error Unsupported platform, or confused platform defines...
#endif

/*** Standard C headers ***/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef TOGL_WGL
#  include <tkPlatDecls.h>
#endif

#if TK_MAJOR_VERSION < 8
#  error Sorry Togl requires Tcl/Tk ver 8.0 or higher.
#endif

#if defined(TOGL_AGL)
#  if TK_MAJOR_VERSION < 8 || (TK_MAJOR_VERSION == 8 && TK_MINOR_VERSION < 4)
#    error Sorry Mac Aqua version requires Tcl/Tk ver 8.4.0 or higher.
#  endif
#endif /* TOGL_AGL */

#if TK_MAJOR_VERSION == 8 && TK_MINOR_VERSION < 6
static void (*SetClassProcsPtr) (Tk_Window, Tk_ClassProcs *, ClientData);
#else
static void (*SetClassProcsPtr) (Tk_Window, const Tk_ClassProcs *, ClientData);
#endif

/* Defaults */
#define DEFAULT_WIDTH           "400"
#define DEFAULT_HEIGHT          "400"
#define DEFAULT_IDENT           ""
#define DEFAULT_TIME            "1"


#ifdef TOGL_WGL
/* Maximum size of a logical palette corresponding to a colormap in color index
 * mode. */
#  define MAX_CI_COLORMAP_SIZE 4096
#  define MAX_CI_COLORMAP_BITS 12

/*
 * OPA TODO. Possible incompatibility because of use of Tk internal structure.
 * Copy of TkWinColormap from tkWinInt.h
 */

typedef struct
{
    HPALETTE palette;           /* Palette handle used when drawing. */
    UINT    size;               /* Number of entries in the palette. */
    int     stale;              /* 1 if palette needs to be realized, otherwise
                                 * 0.  If the palette is stale, then an idle
                                 * handler is scheduled to realize the palette. */
    Tcl_HashTable refCounts;    /* Hash table of palette entry reference counts
                                 * indexed by pixel value. */
} TkWinColormap;

static  LRESULT(CALLBACK *tkWinChildProc) (HWND hwnd, UINT message,
        WPARAM wParam, LPARAM lParam) = NULL;

#endif /* TOGL_WGL */

/* The constant DUMMY_WINDOW is used to signal window creation failure from
 * Togl_MakeWindow() */
#define DUMMY_WINDOW ((Window) NULL)

#define ALL_EVENTS_MASK         \
   (KeyPressMask                \
   |KeyReleaseMask              \
   |ButtonPressMask             \
   |ButtonReleaseMask           \
   |EnterWindowMask             \
   |LeaveWindowMask             \
   |PointerMotionMask           \
   |ExposureMask                \
   |VisibilityChangeMask        \
   |FocusChangeMask             \
   |PropertyChangeMask          \
   |ColormapChangeMask)

/* TCL3D Begin add Debug */
static Tcl_Channel sDebugOutChan;   /* Output channel for debug messages */
static char        sDebugBuf[512];  /* Buffer for debug messages. */
#ifdef DEBUG
#  define DEBUG_OUT(buf) if (sDebugOutChan) { Tcl_WriteChars (sDebugOutChan, buf, -1); }
#else
#  define DEBUG_OUT(buf)
#endif
/* TCL3D End add Debug */

/*
 * Stuff we initialize on a per package (Togl_Init) basis.
 * Since Tcl uses one interpreter per thread, any per-thread
 * data goes here.
 */
struct Togl_PackageGlobals
{
    Tk_OptionTable optionTable; /* Used to parse options */
    Togl   *toglHead;           /* Head of linked list of all Togl widgets */
    int     nextContextTag;     /* Used to assign similar context tags */
};
typedef struct Togl_PackageGlobals Togl_PackageGlobals;

struct Togl
{
    Togl   *Next;               /* next in linked list */

#if defined(TOGL_WGL)
    HGLRC   Ctx;                /* OpenGL rendering context to be made current */
    HDC     tglGLHdc;           /* Device context of device that OpenGL calls
                                 * will be drawn on */
    int     CiColormapSize;     /* (Maximum) size of colormap in color index
                                 * mode */
#elif defined(TOGL_X11)
    GLXContext Ctx;             /* Normal planes GLX context */
#elif defined(TOGL_AGL)
    AGLContext Ctx;
#endif
    int     contextTag;         /* all contexts with same tag share display
                                 * lists */

    XVisualInfo *VisInfo;       /* Visual info of the current */

    Display *display;           /* X's token for the window's display. */
    Tk_Window TkWin;            /* Tk window structure */
    Tcl_Interp *Interp;         /* Tcl interpreter */
    Tcl_Command widgetCmd;      /* Token for togl's widget command */
    Togl_PackageGlobals *tpg;   /* Used to access globals */
#ifndef NO_TK_CURSOR
    Tk_Cursor Cursor;           /* The widget's cursor */
#endif
    int     Width, Height;      /* Dimensions of window */
    int     SetGrid;            /* positive is grid size for window manager */
    int     TimerInterval;      /* Time interval for timer in milliseconds */
    Tcl_TimerToken timerHandler; /* Token for togl's timer handler */
    Bool    RgbaFlag;           /* configuration flags (ala GLX parameters) */
    int     RgbaRed;
    int     RgbaGreen;
    int     RgbaBlue;
    Bool    DoubleFlag;
    Bool    DepthFlag;
    int     DepthSize;
    Bool    AccumFlag;
    int     AccumRed;
    int     AccumGreen;
    int     AccumBlue;
    int     AccumAlpha;
    Bool    AlphaFlag;
    int     AlphaSize;
    Bool    StencilFlag;
    int     StencilSize;
    Bool    PrivateCmapFlag;
    Bool    OverlayFlag;
    int     AuxNumber;
    Bool    Indirect;
    int     PixelFormat;
    /* TCL3D Begin add ContextAttribs */
    int     MajorVersion;
    int     MinorVersion;
    Bool    ForwardCompatibleContext;
    Bool    DebugContext;
    Bool    CoreProfile;
    /* TCL3D End add ContextAttribs */
    /* TCL3D Begin add SwapInterval */
    int     SwapInterval;
    /* TCL3D End add SwapInterval */
    /* TCL3D Begin add Multisample */
    int     MultisampleBuffers;
    int     MultisampleSamples;
    /* TCL3D End add Multisample */
    Bool    FullscreenFlag;
    const char *ShareList;      /* name (ident) of Togl to share dlists with */
    const char *ShareContext;   /* name (ident) to share OpenGL context with */

    const char *Ident;          /* User's identification string */
    ClientData Client_Data;     /* Pointer to user data */

    Bool    UpdatePending;      /* Should normal planes be redrawn? */

    Tcl_Obj *CreateProc;        /* Callback when widget is realized */
    Tcl_Obj *DisplayProc;       /* Callback when widget is redrawn */
    Tcl_Obj *ReshapeProc;       /* Callback when window size changes */
    Tcl_Obj *DestroyProc;       /* Callback when widget is destroyed */
    Tcl_Obj *TimerProc;         /* Callback when widget is idle */

    /* Overlay stuff */
#if defined(TOGL_X11)
    GLXContext OverlayCtx;      /* Overlay planes OpenGL context */
#elif defined(TOGL_WGL)
    HGLRC   tglGLOverlayHglrc;
#endif

    Window  OverlayWindow;              /* The overlay window, or 0 */
    Tcl_Obj *OverlayDisplayProc;        /* Overlay redraw proc */
    Bool    OverlayUpdatePending;       /* Should overlay be redrawn? */
    Colormap OverlayCmap;               /* colormap for overlay is created */
    int     OverlayTransparentPixel;    /* transparent pixel */
    Bool    OverlayIsMapped;

    GLfloat *RedMap;            /* Index2RGB Maps for Color index modes */
    GLfloat *GreenMap;
    GLfloat *BlueMap;
    GLint   MapSize;            /* = Number of indices in our Togl */
    int     badWindow;          /* true when Togl_MakeWindow fails or should
                                 * create a dummy window */
};



/*
 * Prototypes for functions local to this file
 */
static int Togl_ObjCmd(ClientData clientData, Tcl_Interp *interp,
                       int objc, Tcl_Obj *const *objv);
static void Togl_ObjCmdDelete(ClientData clientData);
static void Togl_EventProc(ClientData clientData, XEvent *eventPtr);
static Window Togl_MakeWindow(Tk_Window, Window, ClientData);
static void Togl_WorldChanged(ClientData);

#ifdef MESA_COLOR_HACK
static int get_free_color_cells(Display *display, int screen,
        Colormap colormap);
static void free_default_color_cells(Display *display, Colormap colormap);
#endif
static void ToglCmdDeletedProc(ClientData);

#if defined(TOGL_AGL)
static void SetMacBufRect(Togl *togl);
#endif

static void Togl_PostRedisplay(Togl *togl);
static void Togl_SwapBuffers(const Togl *togl);
static GLuint Togl_LoadBitmapFont (const Togl *togl, const char *fontname,
                                   const char *fontsize, const char *weight,
                                   const char *slant);
static void Togl_UnloadBitmapFont(const Togl *togl, GLuint fontbase);
static void Togl_UseLayer(Togl *togl, int layer);
static void Togl_ShowOverlay(Togl *togl);
static void Togl_HideOverlay(Togl *togl);
static void Togl_PostOverlayRedisplay(Togl *togl);
static int Togl_ExistsOverlay(const Togl *togl);
static int Togl_IsMappedOverlay(const Togl *togl);
static int Togl_GetOverlayTransparentValue(const Togl *togl);
static int Togl_ContextTag(const Togl *togl);
static int Togl_Width(const Togl *togl);
static int Togl_Height(const Togl *togl);

/*
 * Setup Togl widget configuration options:
 */

#define GEOMETRY_MASK 0x1       /* widget geometry */
#define FORMAT_MASK 0x2         /* pixel format */
#define CURSOR_MASK 0x4
#define TIMER_MASK 0x8
#define OVERLAY_MASK 0x10
#define SWAP_MASK 0x20
#define STEREO_MASK 0x40
#define STEREO_FORMAT_MASK 0x80
#define MULTISAMPLE_MASK 0x100
#define PROFILE_MASK 0x200

static Tk_OptionSpec optionSpecs[] = {
    {TK_OPTION_PIXELS, "-height", "height", "Height",
        DEFAULT_HEIGHT, -1, offsetof(Togl, Height), 0, NULL, GEOMETRY_MASK},
    {TK_OPTION_PIXELS, "-width", "width", "Width",
        DEFAULT_WIDTH, -1, offsetof(Togl, Width), 0, NULL, GEOMETRY_MASK},
    {TK_OPTION_BOOLEAN, "-rgba", "rgba", "Rgba",
        "true", -1, offsetof(Togl, RgbaFlag), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-redsize", "redsize", "RedSize",
        "1", -1, offsetof(Togl, RgbaRed), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-greensize", "greensize", "GreenSize",
        "1", -1, offsetof(Togl, RgbaGreen), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-bluesize", "bluesize", "BlueSize",
        "1", -1, offsetof(Togl, RgbaBlue), 0, NULL, FORMAT_MASK},
    {TK_OPTION_BOOLEAN, "-double", "double", "Double",
        "false", -1, offsetof(Togl, DoubleFlag), 0, NULL, FORMAT_MASK},
    {TK_OPTION_BOOLEAN, "-depth", "depth", "Depth",
        "false", -1, offsetof(Togl, DepthFlag), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-depthsize", "depthsize", "DepthSize",
        "1", -1, offsetof(Togl, DepthSize), 0, NULL, FORMAT_MASK},
    {TK_OPTION_BOOLEAN, "-accum", "accum", "Accum",
        "false", -1, offsetof(Togl, AccumFlag), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-accumredsize", "accumredsize", "AccumRedSize",
        "1", -1, offsetof(Togl, AccumRed), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-accumgreensize", "accumgreensize",
        "AccumGreenSize", "1", -1, offsetof(Togl, AccumGreen), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-accumbluesize", "accumbluesize",
        "AccumBlueSize", "1", -1, offsetof(Togl, AccumBlue), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-accumalphasize", "accumalphasize",
        "AccumAlphaSize", "1", -1, offsetof(Togl, AccumAlpha), 0, NULL, FORMAT_MASK},
    {TK_OPTION_BOOLEAN, "-alpha", "alpha", "Alpha",
        "false", -1, offsetof(Togl, AlphaFlag), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-alphasize", "alphasize", "AlphaSize",
        "1", -1, offsetof(Togl, AlphaSize), 0, NULL, FORMAT_MASK},
    {TK_OPTION_BOOLEAN, "-stencil", "stencil", "Stencil",
        "false", -1, offsetof(Togl, StencilFlag), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-stencilsize", "stencilsize", "StencilSize",
        "1", -1, offsetof(Togl, StencilSize), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-auxbuffers", "auxbuffers", "AuxBuffers",
        "0", -1, offsetof(Togl, AuxNumber), 0, NULL, FORMAT_MASK},
    {TK_OPTION_BOOLEAN, "-privatecmap", "privateCmap", "PrivateCmap",
        "false", -1, offsetof(Togl, PrivateCmapFlag), 0, NULL, FORMAT_MASK},
    {TK_OPTION_BOOLEAN, "-overlay", "overlay", "Overlay",
        "false", -1, offsetof(Togl, OverlayFlag), 0, NULL, OVERLAY_MASK},
    {TK_OPTION_CURSOR, "-cursor", "cursor", "Cursor",
        "", -1, offsetof(Togl, Cursor), TK_OPTION_NULL_OK, NULL, CURSOR_MASK},
    {TK_OPTION_INT, "-setgrid", "setGrid", "SetGrid",
        "0", -1, offsetof(Togl, SetGrid), 0, NULL, GEOMETRY_MASK},
    {TK_OPTION_INT, "-time", "time", "Time",
        DEFAULT_TIME, -1, offsetof(Togl, TimerInterval), 0, NULL, TIMER_MASK},
    {TK_OPTION_STRING, "-sharelist", "sharelist", "ShareList",
        NULL, -1, offsetof(Togl, ShareList), 0, NULL, FORMAT_MASK},
    {TK_OPTION_STRING, "-sharecontext", "sharecontext",
        "ShareContext", NULL, -1, offsetof(Togl, ShareContext), 0, NULL, FORMAT_MASK},
    {TK_OPTION_STRING, "-ident", "ident", "Ident",
        DEFAULT_IDENT, -1, offsetof(Togl, Ident), 0, NULL, 0},
    {TK_OPTION_BOOLEAN, "-indirect", "indirect", "Indirect",
        "false", -1, offsetof(Togl, Indirect), 0, NULL, FORMAT_MASK},
    {TK_OPTION_INT, "-pixelformat", "pixelFormat", "PixelFormat",
        "0", -1, offsetof(Togl, PixelFormat), 0, NULL, FORMAT_MASK},
    /* TCL3D Begin add ContextAttribs */
    {TK_OPTION_BOOLEAN, "-coreprofile", "coreProfile", "CoreProfile",
        "false", -1, offsetof(Togl, CoreProfile), 0, NULL, PROFILE_MASK},
    {TK_OPTION_BOOLEAN, "-debug", "debug", "Debug",
        "false", -1, offsetof(Togl, DebugContext), 0, NULL, PROFILE_MASK},
    {TK_OPTION_BOOLEAN, "-forwardcompatible", "forwardcompatible", "Forwardcompatible",
        "false", -1, offsetof(Togl, ForwardCompatibleContext), 0, NULL, PROFILE_MASK},
    {TK_OPTION_INT, "-major", "major", "Major",
        "1", -1, offsetof(Togl, MajorVersion), 0, NULL, PROFILE_MASK},
    {TK_OPTION_INT, "-minor", "minor", "Minor",
        "0", -1, offsetof(Togl, MinorVersion), 0, NULL, PROFILE_MASK},
    /* TCL3D End add ContextAttribs */
    /* TCL3D Begin add SwapInterval */
    {TK_OPTION_INT, "-swapinterval", "swapInterval", "SwapInterval",
        "0", -1, offsetof(Togl, SwapInterval), 0, NULL, SWAP_MASK},
    /* TCL3D End add SwapInterval */
    /* TCL3D Begin add Multisample */
    {TK_OPTION_INT, "-multisamplebuffers", "multisampleBuffers",
        "MultisampleBuffers", "0", -1, offsetof(Togl, MultisampleBuffers),
        0, NULL, MULTISAMPLE_MASK},
    {TK_OPTION_INT, "-multisamplesamples", "multisampleSamples",
        "MultisampleSamples", "2", -1, offsetof(Togl, MultisampleSamples),
        0, NULL, MULTISAMPLE_MASK},
    /* TCL3D End add Multisample */
    {TK_OPTION_BOOLEAN, "-fullscreen", "fullscreen", "Fullscreen",
        "false", -1, offsetof(Togl, FullscreenFlag), 0, NULL, GEOMETRY_MASK|FORMAT_MASK},
    {TK_OPTION_STRING, "-createcommand", "createCommand",
        "CallbackCommand", NULL, offsetof(Togl, CreateProc), -1, TK_OPTION_NULL_OK, NULL, 0},
    {TK_OPTION_SYNONYM, "-create", NULL, NULL,
        NULL, -1, -1, 0, (ClientData) "-createcommand", 0},
    {TK_OPTION_STRING, "-displaycommand", "displayCommand",
        "CallbackCommand", NULL, offsetof(Togl, DisplayProc), -1, TK_OPTION_NULL_OK, NULL, 0},
    {TK_OPTION_SYNONYM, "-display", NULL, NULL,
        NULL, -1, -1, 0, (ClientData) "-displaycommand", 0},
    {TK_OPTION_STRING, "-reshapecommand", "reshapeCommand",
        "CallbackCommand", NULL, offsetof(Togl, ReshapeProc), -1, TK_OPTION_NULL_OK, NULL, 0},
    {TK_OPTION_SYNONYM, "-reshape", NULL, NULL,
        NULL, -1, -1, 0, (ClientData) "-reshapecommand", 0},
    {TK_OPTION_STRING, "-destroycommand", "destroyCommand",
        "CallbackCommand", NULL, offsetof(Togl, DestroyProc), -1, TK_OPTION_NULL_OK, NULL, 0},
    {TK_OPTION_SYNONYM, "-destroy", NULL, NULL,
        NULL, -1, -1, 0, (ClientData) "-destroycommand", 0},
    {TK_OPTION_STRING, "-timercommand", "timerCommand",
        "CallbackCommand", NULL, offsetof(Togl, TimerProc), -1, TK_OPTION_NULL_OK, NULL, 0},
    {TK_OPTION_SYNONYM, "-timer", NULL, NULL,
        NULL, -1, -1, 0, (ClientData) "-timercommand", 0},
    {TK_OPTION_STRING, "-overlaydisplaycommand",
        "overlaydisplayCommand", "CallbackCommand", NULL,
        offsetof(Togl, OverlayDisplayProc), -1,
        TK_OPTION_NULL_OK, NULL, OVERLAY_MASK},
    {TK_OPTION_SYNONYM, "-overlaydisplay", NULL, NULL,
        NULL, -1, -1, 0, (ClientData) "-overlaydisplaycommand", 0},
    /* Tcl3D backwards compatibility */
    {TK_OPTION_SYNONYM, "-createproc", NULL, NULL,
        NULL, -1, -1, 0, (ClientData) "-createcommand", 0},
    {TK_OPTION_SYNONYM, "-displayproc", NULL, NULL,
        NULL, -1, -1, 0, (ClientData) "-displaycommand", 0},
    {TK_OPTION_SYNONYM, "-reshapeproc", NULL, NULL,
        NULL, -1, -1, 0, (ClientData) "-reshapecommand", 0},
    /* end Tcl3D compatibility */
    {TK_OPTION_END, NULL, NULL, NULL, NULL, -1, -1, 0, NULL, 0}
};

/*
 * Add given togl widget to linked list.
 */
static void
AddToList(Togl *t)
{
    t->Next = t->tpg->toglHead;
    t->tpg->toglHead = t;
}
/*
 * Remove given togl widget from linked list.
 */
static void
RemoveFromList(Togl *t)
{
    Togl *prev;
    Togl *cur;

    for (cur = t->tpg->toglHead, prev = NULL; cur; prev = cur, cur = cur->Next) {
        if (t != cur)
            continue;
        if (prev) {
            prev->Next = cur->Next;
        } else {
            t->tpg->toglHead = cur->Next;
        }
        break;
    }
    if (cur)
        cur->Next = NULL;
}

/*
 * Return pointer to togl widget given a user identifier string.
 */
static Togl *
FindTogl(Togl *togl, const char *ident)
{
    Togl *t;

    if (ident[0] != '.') {
        for (t = togl->tpg->toglHead; t; t = t->Next) {
            if (strcmp(t->Ident, ident) == 0)
                break;
        }
    } else {
        for (t = togl->tpg->toglHead; t; t = t->Next) {
            const char *pathname = Tk_PathName(t->TkWin);

            if (strcmp(pathname, ident) == 0)
                break;
        }
    }
    return t;
}


/*
 * Return pointer to another togl widget with same OpenGL context.
 */
static Togl *
FindToglWithSameContext(const Togl *togl)
{
    Togl *t;

    for (t = togl->tpg->toglHead; t != NULL; t = t->Next) {
        if (t == togl)
            continue;

        if (t->Ctx == togl->Ctx) {
            return t;
        }
    }
    return NULL;
}

#if TOGL_USE_OVERLAY
/*
 * Return pointer to another togl widget with same OpenGL overlay context.
 */
static Togl *
FindToglWithSameOverlayContext(const Togl *togl)
{
    Togl *t;

    for (t = togl->tpg->toglHead; t != NULL; t = t->Next) {
        if (t == togl)
            continue;
#if defined(TOGL_X11)
        if (t->OverlayCtx == togl->OverlayCtx)
#elif defined(TOGL_WGL)
        if (t->tglGLOverlayHglrc == togl->tglGLOverlayHglrc)
#endif
        {
            return t;
        }
    }
    return NULL;
}
#endif

#if defined(TOGL_X11)
/*
 * Return an X colormap to use for OpenGL RGB-mode rendering.
 * Input:  dpy - the X display
 *         scrnum - the X screen number
 *         visinfo - the XVisualInfo as returned by glXChooseVisual()
 * Return:  an X Colormap or 0 if there's a _serious_ error.
 */
static Colormap
get_rgb_colormap(Display *dpy,
        int scrnum, const XVisualInfo *visinfo, Tk_Window tkwin)
{
    Atom    hp_cr_maps;
    Status  status;
    int     numCmaps;
    int     i;
    XStandardColormap *standardCmaps;
    Window  root = XRootWindow(dpy, scrnum);
    Bool    using_mesa;

    /*
     * First check if visinfo's visual matches the default/root visual.
     */
    if (visinfo->visual == Tk_Visual(tkwin)) {
        /* use the default/root colormap */
        Colormap cmap;

        cmap = Tk_Colormap(tkwin);
#  ifdef MESA_COLOR_HACK
        (void) get_free_color_cells(dpy, scrnum, cmap);
#  endif
        return cmap;
    }

    /*
     * Check if we're using Mesa.
     */
    if (strstr(glXQueryServerString(dpy, scrnum, GLX_VERSION), "Mesa")) {
        using_mesa = True;
    } else {
        using_mesa = False;
    }

    /*
     * Next, if we're using Mesa and displaying on an HP with the "Color
     * Recovery" feature and the visual is 8-bit TrueColor, search for a
     * special colormap initialized for dithering.  Mesa will know how to
     * dither using this colormap.
     */
    if (using_mesa) {
        hp_cr_maps = XInternAtom(dpy, "_HP_RGB_SMOOTH_MAP_LIST", True);
        if (hp_cr_maps
#  ifdef __cplusplus
                && visinfo->visual->c_class == TrueColor
#  else
                && visinfo->visual->class == TrueColor
#  endif
                && visinfo->depth == 8) {
            status = XGetRGBColormaps(dpy, root, &standardCmaps,
                    &numCmaps, hp_cr_maps);
            if (status) {
                for (i = 0; i < numCmaps; i++) {
                    if (standardCmaps[i].visualid == visinfo->visual->visualid) {
                        Colormap cmap = standardCmaps[i].colormap;

                        (void) XFree(standardCmaps);
                        return cmap;
                    }
                }
                (void) XFree(standardCmaps);
            }
        }
    }

    /*
     * Next, try to find a standard X colormap.
     */
#  if !defined(HP) && !defined(SUN)
#    ifndef SOLARIS_BUG
    status = XGetRGBColormaps(dpy, root, &standardCmaps, &numCmaps, XA_RGB_DEFAULT_MAP);
    if (status == 1) {
        for (i = 0; i < numCmaps; i++) {
            if (standardCmaps[i].visualid == visinfo->visualid) {
                Colormap cmap = standardCmaps[i].colormap;
                (void) XFree(standardCmaps);
                return cmap;
            }
        }
        (void) XFree(standardCmaps);
    }
#    endif
#  endif

    /*
     * If we get here, give up and just allocate a new colormap.
     */
    return XCreateColormap(dpy, root, visinfo->visual, AllocNone);
}
#elif defined(TOGL_WGL)

/* Code to create RGB palette is taken from the GENGL sample program of Win32
 * SDK */

static const unsigned char threeto8[8] = {
    0, 0111 >> 1, 0222 >> 1, 0333 >> 1, 0444 >> 1, 0555 >> 1, 0666 >> 1, 0377
};

static const unsigned char twoto8[4] = {
    0, 0x55, 0xaa, 0xff
};

static const unsigned char oneto8[2] = {
    0, 255
};

static const int defaultOverride[13] = {
    0, 3, 24, 27, 64, 67, 88, 173, 181, 236, 247, 164, 91
};

static const PALETTEENTRY defaultPalEntry[20] = {
    {0, 0, 0, 0},
    {0x80, 0, 0, 0},
    {0, 0x80, 0, 0},
    {0x80, 0x80, 0, 0},
    {0, 0, 0x80, 0},
    {0x80, 0, 0x80, 0},
    {0, 0x80, 0x80, 0},
    {0xC0, 0xC0, 0xC0, 0},

    {192, 220, 192, 0},
    {166, 202, 240, 0},
    {255, 251, 240, 0},
    {160, 160, 164, 0},

    {0x80, 0x80, 0x80, 0},
    {0xFF, 0, 0, 0},
    {0, 0xFF, 0, 0},
    {0xFF, 0xFF, 0, 0},
    {0, 0, 0xFF, 0},
    {0xFF, 0, 0xFF, 0},
    {0, 0xFF, 0xFF, 0},
    {0xFF, 0xFF, 0xFF, 0}
};

static unsigned char
ComponentFromIndex(int i, UINT nbits, UINT shift)
{
    unsigned char val;

    val = (unsigned char) (i >> shift);
    switch (nbits) {

      case 1:
          val &= 0x1;
          return oneto8[val];

      case 2:
          val &= 0x3;
          return twoto8[val];

      case 3:
          val &= 0x7;
          return threeto8[val];

      default:
          return 0;
    }
}

static Colormap
Win32CreateRgbColormap(PIXELFORMATDESCRIPTOR pfd)
{
    TkWinColormap *cmap = (TkWinColormap *) ckalloc(sizeof (TkWinColormap));
    LOGPALETTE *pPal;
    int     n, i;

    n = 1 << pfd.cColorBits;
    pPal = (PLOGPALETTE) LocalAlloc(LMEM_FIXED, sizeof (LOGPALETTE)
            + n * sizeof (PALETTEENTRY));
    pPal->palVersion = 0x300;
    pPal->palNumEntries = n;
    for (i = 0; i < n; i++) {
        pPal->palPalEntry[i].peRed =
                ComponentFromIndex(i, pfd.cRedBits, pfd.cRedShift);
        pPal->palPalEntry[i].peGreen =
                ComponentFromIndex(i, pfd.cGreenBits, pfd.cGreenShift);
        pPal->palPalEntry[i].peBlue =
                ComponentFromIndex(i, pfd.cBlueBits, pfd.cBlueShift);
        pPal->palPalEntry[i].peFlags = 0;
    }

    /* fix up the palette to include the default GDI palette */
    if ((pfd.cColorBits == 8)
            && (pfd.cRedBits == 3) && (pfd.cRedShift == 0)
            && (pfd.cGreenBits == 3) && (pfd.cGreenShift == 3)
            && (pfd.cBlueBits == 2) && (pfd.cBlueShift == 6)) {
        for (i = 1; i <= 12; i++)
            pPal->palPalEntry[defaultOverride[i]] = defaultPalEntry[i];
    }

    cmap->palette = CreatePalette(pPal);
    LocalFree(pPal);
    cmap->size = n;
    cmap->stale = 0;

    /* Since this is a private colormap of a fix size, we do not need a valid
     * hash table, but a dummy one */

    Tcl_InitHashTable(&cmap->refCounts, TCL_ONE_WORD_KEYS);
    return (Colormap) cmap;
}

static Colormap
Win32CreateCiColormap(Togl *togl)
{
    /* Create a colormap with size of togl->CiColormapSize and set all entries
     * to black */

    LOGPALETTE logPalette;
    TkWinColormap *cmap = (TkWinColormap *) ckalloc(sizeof (TkWinColormap));

    logPalette.palVersion = 0x300;
    logPalette.palNumEntries = 1;
    logPalette.palPalEntry[0].peRed = 0;
    logPalette.palPalEntry[0].peGreen = 0;
    logPalette.palPalEntry[0].peBlue = 0;
    logPalette.palPalEntry[0].peFlags = 0;

    cmap->palette = CreatePalette(&logPalette);
    cmap->size = togl->CiColormapSize;
    ResizePalette(cmap->palette, cmap->size);   /* sets new entries to black */
    cmap->stale = 0;

    /* Since this is a private colormap of a fix size, we do not need a valid
     * hash table, but a dummy one */

    Tcl_InitHashTable(&cmap->refCounts, TCL_ONE_WORD_KEYS);
    return (Colormap) cmap;
}
#endif

#ifndef STRINGIFY
#  define STRINGIFY(x) STRINGIFY1(x)
#  define STRINGIFY1(x) #x
#endif

/*
 * Togl_Init
 *
 *   Called upon system startup to create Togl command.
 */
int
BrlcadTkObolHost_Init(Tcl_Interp *interp)
{
#ifdef USE_TCL_STUBS
    if (Tcl_InitStubs(interp, "8.6-", 0) == NULL) {
        return TCL_ERROR;
    }
#endif
#ifdef USE_TK_STUBS
    if (Tk_InitStubs(interp, "8.6-", 0) == NULL) {
        return TCL_ERROR;
    }
#endif

#ifdef USE_TK_STUBS
    SetClassProcsPtr = tkStubsPtr->tk_SetClassProcs;
#else
    SetClassProcsPtr = Tk_SetClassProcs;
#endif

    if (Tcl_FindCommand(interp, "::brlcad::tkobol::host", NULL,
            TCL_GLOBAL_ONLY)) {
        return TCL_OK;
    }

    if (!Tcl_FindNamespace(interp, "::brlcad", NULL, TCL_GLOBAL_ONLY) &&
            !Tcl_CreateNamespace(interp, "::brlcad", NULL, NULL)) {
        return TCL_ERROR;
    }
    if (!Tcl_FindNamespace(interp, "::brlcad::tkobol", NULL,
            TCL_GLOBAL_ONLY) && !Tcl_CreateNamespace(interp,
            "::brlcad::tkobol", NULL, NULL)) {
        return TCL_ERROR;
    }

    if (Tcl_CreateObjCommand(interp, "::brlcad::tkobol::host", Togl_ObjCmd, NULL,
                    Togl_ObjCmdDelete) == NULL) {
        return TCL_ERROR;
    }

    if (Tcl_PkgProvide(interp, "BrlcadTkObolHost", TOGL_VERSION) != TCL_OK) {
        return TCL_ERROR;
    }

    return TCL_OK;
}

/*
 * Togl_CallCallback
 *
 * Call command with togl widget as only argument
 */

int
Togl_CallCallback(Togl *togl, Tcl_Obj *cmd)
{
    int     result;
    Tcl_Obj *objv[3];

    if (cmd == NULL || togl->widgetCmd == NULL) {
        return TCL_OK;
    }

    objv[0] = cmd;
    Tcl_IncrRefCount(objv[0]);
    objv[1] =
            Tcl_NewStringObj(Tcl_GetCommandName(togl->Interp, togl->widgetCmd),
            -1);
    Tcl_IncrRefCount(objv[1]);
    objv[2] = NULL;
    result = Tcl_EvalObjv(togl->Interp, 2, objv, TCL_EVAL_GLOBAL);
    Tcl_DecrRefCount(objv[1]);
    Tcl_DecrRefCount(objv[0]);
    if (result != TCL_OK) {
        /* Use the debug buffer to write error message to stdout. This is needed,
         * if the callback is called before the event loop is entered.
         */
        sprintf (sDebugBuf, "Error in %s: %s\n", Tcl_GetString(cmd), Tcl_GetStringResult(togl->Interp));
        if (sDebugOutChan) { Tcl_WriteChars (sDebugOutChan, sDebugBuf, -1); }
        Tcl_BackgroundException(togl->Interp, result);
    }
    return result;
}

/*
 * Togl_Timer
 *
 * Gets called from Tk_CreateTimerHandler.
 */
static void
Togl_Timer(ClientData clientData)
{
    Togl *togl = (Togl *) clientData;

    if (togl->TimerProc) {
        if (Togl_CallCallback(togl, togl->TimerProc) != TCL_OK) {
            togl->timerHandler = NULL;
            Tcl_BackgroundError (togl->Interp);
            return;
        }
        /*
         * Re-register this callback since Tcl/Tk timers are "one-shot".
         * That is, after the timer callback is called it not normally
         * called again.  That's not the behavior we want for Togl.
         */
        togl->timerHandler =
                Tcl_CreateTimerHandler(togl->TimerInterval, Togl_Timer,
                (ClientData) togl);
    }
}

/*
 * Togl_MakeCurrent
 *
 *   Bind the OpenGL rendering context to the specified
 *   Togl widget.  If given a NULL argument, then the
 *   OpenGL context is released without assigning a new one.
 */
void
Togl_MakeCurrent(const Togl *togl)
{
#if defined(TOGL_WGL)
    if (togl == NULL) {
        HDC hdc = wglGetCurrentDC();
        if (hdc != NULL) {
            (void) wglMakeCurrent(hdc, NULL);
        }
    } else {
        (void) wglMakeCurrent(togl->tglGLHdc, togl->Ctx);
    }

#elif defined(TOGL_X11)
    Display *display = togl ? togl->display : glXGetCurrentDisplay();

    if (display) {
        GLXDrawable drawable;

        if (!togl)
            drawable = None;
        else if (togl->TkWin)
            drawable = Tk_WindowId(togl->TkWin);
        else
            drawable = None;
        (void) glXMakeCurrent(display, drawable, drawable ? togl->Ctx : NULL);
    }

#elif defined(TOGL_AGL)
   if (togl == NULL || togl->Ctx == NULL) {
        (void) aglSetCurrentContext(NULL);
    } else {
        (void) aglSetCurrentContext(togl->Ctx);
    }
#endif
}

/* TCL3D Begin add SwapInterval */
Bool
Togl_SwapInterval(const Togl *togl, int interval)
{
#ifdef TOGL_AGL
    GLint swapInterval = interval;

    return aglSetInteger(togl->Ctx, AGL_SWAP_INTERVAL, &swapInterval);
#endif

#ifdef TOGL_WGL
    typedef BOOL (WINAPI *BOOLFuncInt)(int);
    typedef const char *(WINAPI *StrFuncHDC)(HDC);
    static BOOLFuncInt swapInterval = NULL;
    static BOOL initialized = False;

    if (!initialized) {
        const char *extensions;
        StrFuncHDC getExtensionsString;

        getExtensionsString = (StrFuncHDC) wglGetProcAddress("wglGetExtensionsStringARB");
        if (getExtensionsString == NULL) {
            getExtensionsString = (StrFuncHDC) wglGetProcAddress("wglGetExtensionsStringEXT");
        }
        if (getExtensionsString) {
            extensions = getExtensionsString(togl->tglGLHdc);
            if (strstr(extensions, "WGL_EXT_swap_control") != NULL) {
                swapInterval = (BOOLFuncInt) wglGetProcAddress("wglSwapIntervalEXT");
            }
        }
        initialized = True;
    }
    if (swapInterval) {
        return swapInterval(interval);
    }
    return False;
#endif

#ifdef TOGL_X11
    typedef int (*IntFuncInt)(int);
    static IntFuncInt swapInterval = NULL;
    static int initialized = False;
    typedef void (*ExtFunc)(Display *, GLXDrawable, int);
    static ExtFunc swapIntervalExt = NULL;

    if (!initialized) {
        const char *extensions = glXQueryExtensionsString(togl->display, Tk_ScreenNumber(togl->TkWin));

        if (strstr(extensions, "GLX_EXT_swap_control") != NULL) {
            swapIntervalExt = (ExtFunc) glXGetProcAddressARB((const GLubyte *)"glXSwapIntervalEXT");
        } else if (strstr(extensions, "GLX_MESA_swap_control") != NULL) {
            swapInterval = (IntFuncInt) glXGetProcAddressARB((const GLubyte *)"glXSwapIntervalMESA");
        } else if (strstr(extensions, "GLX_SGI_swap_control") != NULL) {
            swapInterval = (IntFuncInt) glXGetProcAddressARB((const GLubyte *)"glXSwapIntervalSGI");
        }
        initialized = True;
    }
    if (swapIntervalExt) {
        swapIntervalExt(togl->display, Tk_WindowId(togl->TkWin), interval);
        return True;
    }
    if (swapInterval) {
        return swapInterval(interval) == 0;
    }
    return False;
#endif
}
/* TCL3D End add SwapInterval */

#if defined(TOGL_AGL)
/* tell OpenGL which part of the Mac window to render to */
static void
SetMacBufRect(Togl *togl)
{
    GLint   wrect[4];
    Rect    r;
    MacDrawable *d = ((TkWindow *) togl->TkWin)->privatePtr;

    /* set wrect[0,1] to lower left corner of widget */
    wrect[2] = Tk_Width(togl->TkWin);
    wrect[3] = Tk_Height(togl->TkWin);
    wrect[0] = d->xOff;

    GetPortBounds(Togl_MacOSXGetDrawablePort(togl), &r);
    wrect[1] = r.bottom - wrect[3] - d->yOff;

    if (togl->FullscreenFlag) {
        aglEnable(togl->Ctx, AGL_FS_CAPTURE_SINGLE);
        aglSetFullScreen(togl->Ctx, 0, 0, 0, 0);
    } else {
        aglUpdateContext(togl->Ctx);
    }
    aglSetInteger(togl->Ctx, AGL_BUFFER_RECT, wrect);
    aglEnable(togl->Ctx, AGL_BUFFER_RECT);
}

static void
ReconfigureCB(CGDirectDisplayID display, CGDisplayChangeSummaryFlags flags,
        void *closure)
{
    /* Display reconfiguration callback. Documented as needed by Apple QA1209.
     * Updated for 10.3 (and later) to use
     * CGDisplayRegisterReconfigurationCallback. */
    Togl   *togl = (Togl *) closure;

    if (0 != (flags & kCGDisplayBeginConfigurationFlag))
        return;                 /* wait until display is reconfigured */

    SetMacBufRect(togl);
    Togl_MakeCurrent(togl);
    if (togl->Ctx) {
        if (togl->ReshapeProc) {
            Togl_CallCallback(togl, togl->ReshapeProc);
        } else {
            glViewport(0, 0, togl->Width, togl->Height);
        }
    }
}
#endif

/*
 * Called when the widget's contents must be redrawn.  Basically, we
 * just call the user's render callback function.
 *
 * Note that the parameter type is ClientData so this function can be
 * passed to Tcl_DoWhenIdle().
 */
static void
Togl_Render(ClientData clientData)
{
    Togl   *togl = (Togl *) clientData;

    if (togl->DisplayProc) {
        Togl_MakeCurrent(togl);
        if (Togl_CallCallback(togl, togl->DisplayProc) != TCL_OK) {
            Tcl_BackgroundError (togl->Interp);
            return;
        }
    }
    togl->UpdatePending = False;
}


static void
Togl_RenderOverlay(ClientData clientData)
{
    Togl *togl = (Togl *) clientData;

    if (togl->OverlayFlag && togl->OverlayDisplayProc) {
#if defined(TOGL_WGL)
        int res = wglMakeCurrent(togl->tglGLHdc, togl->tglGLOverlayHglrc);
        if (!res) {
            Tcl_Panic("wglMakeCurrent overlay");
        }

#elif defined(TOGL_X11)
        (void) glXMakeCurrent(Tk_Display(togl->TkWin),
                togl->OverlayWindow, togl->OverlayCtx);

#endif /* TOGL_WGL */

        if (Togl_CallCallback(togl, togl->OverlayDisplayProc) != TCL_OK) {
            Tcl_BackgroundError (togl->Interp);
            return;
        }
    }
    togl->OverlayUpdatePending = False;
}

/*
 * See domentation about what can't be changed
 */
static int
Togl_ObjConfigure(Tcl_Interp *interp, Togl *togl,
        int objc, Tcl_Obj *const *objv)
{
    Tk_SavedOptions savedOptions;
    int     error, mask;
    int     undoMask = 0;
    Tcl_Obj *errorResult = NULL;

    for (error = 0; error <= 1; ++error, mask = undoMask) {
        if (error == 0) {
            /*
             * Tk_SetOptions parses the command arguments
             * and looks for defaults in the resource database.
             */
            if (Tk_SetOptions(interp, WIDGREC togl, togl->tpg->optionTable,
                            objc, objv, togl->TkWin, &savedOptions, &mask)
                    != TCL_OK) {
                /* previous values are restored, so nothing to do */
                return TCL_ERROR;
            }
        } else {
            /*
             * Restore options from saved values
             */
            errorResult = Tcl_GetObjResult(interp);
            Tcl_IncrRefCount(errorResult);
            Tk_RestoreSavedOptions(&savedOptions);
        }

        if (togl->Ident && togl->Ident[0] == '.') {
            Tcl_AppendResult(interp, "Can not set ident to a window path name", NULL);
            continue;
        }

        if (togl->FullscreenFlag) {
            /* override width and height */
            togl->Width = WidthOfScreen(Tk_Screen(togl->TkWin));
            togl->Height = HeightOfScreen(Tk_Screen(togl->TkWin));
            undoMask |= GEOMETRY_MASK;
        }

        if (mask & GEOMETRY_MASK) {
            Togl_WorldChanged((ClientData) togl);
            undoMask |= GEOMETRY_MASK;
        }

        if (mask & OVERLAY_MASK) {
#if !TOGL_USE_OVERLAY
            if (togl->OverlayFlag) {
                Tcl_AppendResult(interp, "Sorry, overlay was disabled", NULL);
                continue;
            }
#else
#  if defined(TOGL_X11)
            if (togl->OverlayCtx)
#  elif defined(TOGL_WGL)
            if (togl->tglGLOverlayHglrc)
#  endif
            {
                /*
                 * Trying to change existing pixel format/graphics context
                 */
                Tcl_AppendResult(interp, "Unable to change overlay pixel format", NULL);
                continue;
            }
#endif
        }

        if (mask & SWAP_MASK) {
            if (togl->Ctx) {
                /*
                 * Change existing swap interval
                 */
                Togl_MakeCurrent(togl); /* TODO: needed? */
                Togl_SwapInterval(togl, togl->SwapInterval);
                undoMask |= SWAP_MASK;
            }
        }

        if (mask & FORMAT_MASK) {
            if (togl->Ctx) {
                /*
                 * Trying to change existing pixel format/graphics context
                 * TODO: (re)create graphics context
                 *
                 * save old graphics context
                 * try to create new one and share display lists
                 * if failure, then restore old one
                 */
                Tcl_AppendResult(interp, "Unable to change pixel format", NULL);
                continue;
            }
            if (togl->ShareContext && togl->ShareList) {
                Tcl_AppendResult(interp,
                        "only one of -sharelist and -sharecontext allowed",
                        NULL);
                continue;
            }
            if (togl->FullscreenFlag) {
#ifndef TOGL_AGL
#  if TK_MAJOR_VERSION < 8 || (TK_MAJOR_VERSION == 8 && TK_MINOR_VERSION < 5)
                Tcl_AppendResult(interp, "Need Tk 8.5 or later for fullscreen support", NULL);
                continue;
#  endif
#endif
            }

            /* Whether or not the format is okay is figured out when togl tries
             * to create the window. */
#ifdef MESA_COLOR_HACK
            free_default_color_cells(Tk_Display(togl->TkWin),
                    Tk_Colormap(togl->TkWin));
#endif
            undoMask |= FORMAT_MASK;
        }

        if (mask & TIMER_MASK) {
            if (togl->timerHandler != NULL) {
                Tcl_DeleteTimerHandler(togl->timerHandler);
            }
            if (togl->TimerProc) {
                togl->timerHandler =
                        Tcl_CreateTimerHandler(togl->TimerInterval, Togl_Timer,
                        (ClientData) togl);
            }
            undoMask |= TIMER_MASK;
        }
        break;
    }

    if (error == 0) {
        Tk_FreeSavedOptions(&savedOptions);
        return TCL_OK;
    } else {
        Tcl_SetObjResult(interp, errorResult);
        Tcl_DecrRefCount(errorResult);
        return TCL_ERROR;
    }
}

static int
Togl_ObjWidget(ClientData clientData, Tcl_Interp *interp, int objc,
               Tcl_Obj *const *objv)
{
    Togl *togl = (Togl *) clientData;
    static const char *commands[] = {
        "cget", "configure", "extensions",
        "postredisplay", "render",
        "swapbuffers", "makecurrent",
        "loadbitmapfont", "unloadbitmapfont",
        "uselayer", "showoverlay", "hideoverlay",
        "postredisplayoverlay", "renderoverlay",
        "existsoverlay", "ismappedoverlay",
        "getoverlaytransparentvalue",
        "contexttag", "width", "height",
        NULL
    };
    enum command
    {
        TOGL_CGET, TOGL_CONFIGURE, TOGL_GL_EXTENSIONS,
        TOGL_POSTREDISPLAY, TOGL_RENDER,
        TOGL_SWAPBUFFERS, TOGL_MAKECURRENT,
        TOGL_LOADBITMAPFONT, TOGL_UNLOADBITMAPFONT,
        TOGL_USELAYER, TOGL_SHOWOVERLAY, TOGL_HIDEOVERLAY,
        TOGL_POSTREDISPLAYOVERLAY, TOGL_RENDEROVERLAY,
        TOGL_EXISTSOVERLAY, TOGL_ISMAPPEDOVERLAY,
        TOGL_GETOVERLAYTRANSPARENTVALUE,
        TOGL_CONTEXTTAG, TOGL_WIDTH, TOGL_HEIGHT
    };
    int     result = TCL_OK;
    Tcl_Obj *objPtr;
    int     index;

    if (objc < 2) {
        Tcl_WrongNumArgs(interp, 1, objv, "command ?arg arg ...?");
        return TCL_ERROR;
    }

    Tcl_Preserve((ClientData) togl);

    result = Tcl_GetIndexFromObj(interp, objv[1], commands, "option", 0, &index);

    if (result == TCL_OK) {
        switch (index) {
            case TOGL_CGET:
                if (objc != 3) {
                    Tcl_WrongNumArgs(interp, 2, objv, "option");
                    result = TCL_ERROR;
                    break;
                }
                objPtr = Tk_GetOptionValue(interp, WIDGREC togl,
                         togl->tpg->optionTable, (objc == 3) ? objv[2] : NULL,
                         togl->TkWin);
                if (objPtr == NULL) {
                    result = TCL_ERROR;
                    break;
                }
                Tcl_SetObjResult(interp, objPtr);
                break;

            case TOGL_CONFIGURE:
                if (objc <= 3) {
                    /*
                     * Return one item if the option is given,
                     * or return all configuration information
                     */
                    objPtr = Tk_GetOptionInfo(interp, WIDGREC togl,
                             togl->tpg->optionTable, (objc == 3) ? objv[2] : NULL,
                             togl->TkWin);
                    if (objPtr == NULL) {
                        result = TCL_ERROR;
                    } else {
                        Tcl_SetObjResult(interp, objPtr);
                    }
                } else {
                    /* Execute a configuration change */
                    result = Togl_ObjConfigure(interp, togl, objc - 2, objv + 2);
                }
                break;

            case TOGL_GL_EXTENSIONS: {
                /* Return a list of available OpenGL extensions. */
                const char *extensions = NULL;
                Tcl_Size length = -1;

                if (objc == 2 || objc == 3) {
                    if (objc == 2) {
                        extensions = (const char *) glGetString(GL_EXTENSIONS);
                    } else if (objc == 3) {
                        const char *option;
                        option = Tcl_GetString(objv[2]);
                        if (strcmp (option, "-gl") == 0) {
                            extensions = (const char *) glGetString(GL_EXTENSIONS);
                        } else if (strcmp (option, "-glu") == 0) {
                            extensions = (const char *) gluGetString(GLU_EXTENSIONS);
                        } else if (strcmp (option, "-platform") == 0) {
#ifdef TOGL_WGL
                            typedef const char *(WINAPI *StrFuncHDC)(HDC);
                            StrFuncHDC getExtensionsString;

                            getExtensionsString = (StrFuncHDC) wglGetProcAddress("wglGetExtensionsStringARB");
                            if (getExtensionsString == NULL) {
                                getExtensionsString = (StrFuncHDC) wglGetProcAddress("wglGetExtensionsStringEXT");
                            }
                            if (getExtensionsString) {
                                extensions = getExtensionsString(togl->tglGLHdc);
                            }
#endif
#ifdef TOGL_X11
                            extensions = glXQueryExtensionsString(togl->display, Tk_ScreenNumber(togl->TkWin));
#endif
                        } else {
                            Tcl_ResetResult(interp);
                            Tcl_AppendResult(interp, "Invalid extensions option \"", option, "\"", (char *)NULL);
                            result = TCL_ERROR;
                        }
                    }
                    if (extensions) {
                        objPtr = Tcl_NewStringObj(extensions, -1);
                        /* convert to list by asking for its length */
                        (void) Tcl_ListObjLength(interp, objPtr, &length);
                        Tcl_SetObjResult(interp, objPtr);
                    } else {
                        Tcl_SetResult(interp, "Extensions string not available.", TCL_STATIC);
                        result = TCL_ERROR;
                    }
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;
            }

            case TOGL_POSTREDISPLAY:
                /* schedule the widget to be redrawn */
                if (objc == 2) {
                    Togl_PostRedisplay(togl);
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_RENDER:
                /* force the widget to be redrawn */
                if (objc == 2) {
                    Togl_Render((ClientData) togl);
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_SWAPBUFFERS:
                /* force the widget to be redrawn */
                if (objc == 2) {
                    Togl_SwapBuffers(togl);
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_MAKECURRENT:
                /* force the widget to be redrawn */
                if (objc == 2) {
                    Togl_MakeCurrent(togl);
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            /* TCL3D Begin change Font */
            case TOGL_LOADBITMAPFONT:
            {
                GLuint fontbase;
                Tcl_Obj *fontbaseAsTclObject;
                const char *name    = NULL;
                const char *names[] = { "fixed", "courier", "helvetica" };
                const char *size    = "12";
                const char *weight  = "*";
                const char *slant   = "*";

                int i = 2;
                const char *optName;
                while (i < objc) {
                    optName  = Tcl_GetString(objv[i]);

                    if (strncmp (optName, "-family", 7) == 0) {
                        if (i < objc-1) {
                            name = Tcl_GetString(objv[++i]);
                        }
                    } else if (strncmp (optName, "-slant", 6) == 0) {
                        if (i < objc-1) {
                            slant = Tcl_GetString(objv[++i]);
                        }
                    } else if (strncmp (optName, "-weight", 7) == 0) {
                        if (i < objc-1) {
                            weight = Tcl_GetString(objv[++i]);
                        }
                    } else if (strncmp (optName, "-size", 5) == 0) {
                        if (i < objc-1) {
                            size = Tcl_GetString(objv[++i]);
                        }
                    } else {
                        /* No option given. The supplied parameter should be a XLFD
                         * specification. We pass it over to Togl_LoadBitmapFont, where
                         * it will be checked for correctness. */
                        name = Tcl_GetString(objv[i]);
                        break;
                    }
                    i++;
                }
                if (name == NULL) {
                    /* No font name specified. Try the default ones. */
                    for (i = 0; i < (int)(sizeof(names) / sizeof(char *)); i++) {
                        name = names[i];
                        fontbase = Togl_LoadBitmapFont(togl, name, size, weight, slant);
                        if (fontbase) {
                            break;
                        }
                    }
                } else {
                    fontbase = Togl_LoadBitmapFont(togl, name, size, weight, slant);
                }
                if (fontbase) {
                    fontbaseAsTclObject = Tcl_NewIntObj(fontbase);
                    Tcl_SetObjResult(interp, fontbaseAsTclObject);
                    result = TCL_OK;
                } else {
                    Tcl_AppendResult(interp, "Could not allocate font \"", name, "\"", NULL);
                    result = TCL_ERROR;
                }
                break;
            }
            /* TCL3D End change Font */

            case TOGL_UNLOADBITMAPFONT:
                if (objc == 3) {
                    int fontbase;
                    result = Tcl_GetIntFromObj(interp, objv[2], &fontbase);
                    if (result == TCL_ERROR)
                        break;
                    Togl_UnloadBitmapFont(togl, fontbase);
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, "fontbase");
                    result = TCL_ERROR;
                }
                break;

            case TOGL_USELAYER:
                if (objc == 3) {
                    int     layer;

                    result = Tcl_GetIntFromObj(interp, objv[2], &layer);
                    if (result == TCL_OK) {
                        Togl_UseLayer(togl, layer);
                    }
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, "layer");
                    result = TCL_ERROR;
                }
                break;

            case TOGL_SHOWOVERLAY:
                if (objc == 2) {
                    Togl_ShowOverlay(togl);
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_HIDEOVERLAY:
                if (objc == 2) {
                    Togl_HideOverlay(togl);
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_POSTREDISPLAYOVERLAY:
                if (objc == 2) {
                    Togl_PostOverlayRedisplay(togl);
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_RENDEROVERLAY:
                /* force the overlay to be redrawn */
                if (objc == 2) {
                    Togl_RenderOverlay((ClientData) togl);
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_EXISTSOVERLAY:
                if (objc == 2) {
                    Tcl_SetObjResult(interp, Tcl_NewIntObj(Togl_ExistsOverlay(togl)));
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_ISMAPPEDOVERLAY:
                if (objc == 2) {
                    Tcl_SetObjResult(interp,
                        Tcl_NewIntObj(Togl_IsMappedOverlay(togl)));
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_GETOVERLAYTRANSPARENTVALUE:
                if (objc == 2) {
                    Tcl_SetObjResult(interp,
                            Tcl_NewIntObj(Togl_GetOverlayTransparentValue(togl)));
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_CONTEXTTAG:
                if (objc == 2) {
                    Tcl_SetObjResult(interp, Tcl_NewIntObj(Togl_ContextTag(togl)));
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            /* TCL3D Begin add Configure */
            case TOGL_WIDTH:
                if (objc == 2) {
                    Tcl_SetObjResult(interp, Tcl_NewIntObj(Togl_Width(togl)));
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;

            case TOGL_HEIGHT:
                if (objc == 2) {
                    Tcl_SetObjResult(interp, Tcl_NewIntObj(Togl_Height(togl)));
                } else {
                    Tcl_WrongNumArgs(interp, 2, objv, NULL);
                    result = TCL_ERROR;
                }
                break;
            /* TCL3D End add Configure */
        }
    }

    Tcl_Release((ClientData) togl);
    return result;
}

/*
 * Togl_ObjCmdDelete
 *
 * Called when togl command is removed from interpreter.
 */

static void
Togl_ObjCmdDelete(ClientData clientData)
{
    if (clientData != NULL) {
        Togl_PackageGlobals *tpg = (Togl_PackageGlobals *) clientData;

        Tk_DeleteOptionTable(tpg->optionTable);
        ckfree((char *) clientData);
    }
}


/*
 * Togl_ObjCmd
 *
 *   Called when Togl is executed - creation of a Togl widget.
 *     * Creates a new window
 *     * Creates an 'Togl' data structure
 *     * Creates an event handler for this window
 *     * Creates a command that handles this object
 *     * Configures this Togl for the given arguments
 */
int
Togl_ObjCmd(ClientData clientData, Tcl_Interp *interp, int objc,
        Tcl_Obj *const *objv)
{
    Togl_PackageGlobals *tpg;
    Togl   *togl;
    Tk_Window tkwin;
    Tcl_InterpState state;

    if (objc <= 1) {
        Tcl_WrongNumArgs(interp, 1, objv, "pathName ?options?");
        return TCL_ERROR;
    }
    tpg = (Togl_PackageGlobals *) clientData;
    if (tpg == NULL) {
        Tcl_CmdInfo info;
        const char *name;

        /*
         * Initialize the Togl_PackageGlobals for this widget the
         * first time a Togl widget is created.  The globals are
         * saved as our client data.
         */

        tpg = (Togl_PackageGlobals *) ckalloc(sizeof (Togl_PackageGlobals));
        if (tpg == NULL) {
            return TCL_ERROR;
        }
        tpg->nextContextTag = 0;
        tpg->optionTable = Tk_CreateOptionTable(interp, optionSpecs);
        tpg->toglHead = NULL;

        name = Tcl_GetString(objv[0]);
        Tcl_GetCommandInfo(interp, name, &info);
        info.objClientData = (ClientData) tpg;
        Tcl_SetCommandInfo(interp, name, &info);
    }

    /* Create the window. */
    tkwin = Tk_CreateWindowFromPath(interp, Tk_MainWindow(interp),
            Tcl_GetString(objv[1]), NULL);
    if (tkwin == NULL) {
        return TCL_ERROR;
    }

    Tk_SetClass(tkwin, "Togl");

    /* Create Togl data structure */
    togl = (Togl *) ckalloc(sizeof (Togl));
    if (togl == NULL) {
        return TCL_ERROR;
    }

    /* initialize Togl data structures values */
    togl->Next = NULL;
    togl->Ctx = NULL;
#if defined(TOGL_WGL)
    togl->tglGLHdc = NULL;
    togl->tglGLOverlayHglrc = NULL;
#elif defined(TOGL_X11)
    togl->OverlayCtx = NULL;
#endif
    togl->contextTag = 0;
    togl->display = Tk_Display(tkwin);
    togl->TkWin = tkwin;
    togl->Interp = interp;
    togl->VisInfo = NULL;
    togl->OverlayWindow = None;
    togl->OverlayCmap = None;
    togl->OverlayTransparentPixel = 0;
    togl->OverlayIsMapped = False;

    togl->UpdatePending = False;
    togl->OverlayUpdatePending = False;
    togl->tpg = tpg;
    togl->Client_Data = NULL;

    /* for color index mode photos */
    togl->RedMap = togl->GreenMap = togl->BlueMap = NULL;
    togl->MapSize = 0;

#ifndef NO_TK_CURSOR
    togl->Cursor = NULL;
#endif
    togl->Width = 0;
    togl->Height = 0;
    togl->SetGrid = 0;
    togl->TimerInterval = 0;
    togl->RgbaFlag = True;
    togl->RgbaRed = 1;
    togl->RgbaGreen = 1;
    togl->RgbaBlue = 1;
    togl->DoubleFlag = False;
    togl->DepthFlag = False;
    togl->DepthSize = 1;
    togl->AccumFlag = False;
    togl->AccumRed = 1;
    togl->AccumGreen = 1;
    togl->AccumBlue = 1;
    togl->AccumAlpha = 1;
    togl->AlphaFlag = False;
    togl->AlphaSize = 1;
    togl->StencilFlag = False;
    togl->StencilSize = 1;
    togl->OverlayFlag = False;
    togl->AuxNumber = 0;
    togl->Indirect = False;
    togl->PixelFormat = 0;
    /* TCL3D Begin add ContextAttribs */
    togl->MajorVersion = 1;
    togl->MinorVersion = 0;
    togl->DebugContext = False;
    togl->ForwardCompatibleContext = False;
    togl->CoreProfile  = False;
    /* TCL3D End add ContextAttribs */
    /* TCL3D Begin add SwapInterval */
    togl->SwapInterval = 0;
    /* TCL3D End add SwapInterval */
    /* TCL3D Begin add Multisample */
    togl->MultisampleBuffers = 0;
    togl->MultisampleSamples = 2;
    /* TCL3D End add Multisample */
    togl->FullscreenFlag = False;

    togl->CreateProc = NULL;
    togl->DisplayProc = NULL;
    togl->ReshapeProc = NULL;
    togl->DestroyProc = NULL;
    togl->TimerProc = NULL;
    togl->timerHandler = NULL;
    togl->OverlayDisplayProc = NULL;
    togl->ShareList = NULL;
    togl->ShareContext = NULL;
    togl->Ident = NULL;
    togl->PrivateCmapFlag = False;
#ifdef HAVE_AUTOSTEREO
    togl->as_initialized = False;
    togl->ash = -1;
#endif
    togl->badWindow = False;

    sDebugOutChan = Tcl_GetStdChannel (TCL_STDOUT);

    /* Create command event handler */
    togl->widgetCmd = Tcl_CreateObjCommand(interp, Tk_PathName(tkwin),
            Togl_ObjWidget, (ClientData) togl, ToglCmdDeletedProc);

    /*
     * Setup the Tk_ClassProcs callbacks to point at our own window creation
     * function.
     */
    if (SetClassProcsPtr != NULL) {     /* use public API (Tk 8.4+) */
        Tk_ClassProcs *procsPtr;

        procsPtr = (Tk_ClassProcs *) ckalloc(sizeof (Tk_ClassProcs));
        procsPtr->size = sizeof (Tk_ClassProcs);
        procsPtr->createProc = Togl_MakeWindow;
        procsPtr->worldChangedProc = Togl_WorldChanged;
        procsPtr->modalProc = NULL;
        /* Tk_SetClassProcs(togl->TkWin, procsPtr, (ClientData) togl); */
        (SetClassProcsPtr) (togl->TkWin, procsPtr, (ClientData) togl);
    }

    Tk_CreateEventHandler(tkwin, ExposureMask | StructureNotifyMask,
            Togl_EventProc, (ClientData) togl);

    /* Configure Togl widget */
    if (Tk_InitOptions(interp, WIDGREC togl, tpg->optionTable, tkwin) != TCL_OK
            || Togl_ObjConfigure(interp, togl, objc - 2, objv + 2) != TCL_OK) {
        goto error;
    }

    /*
     * If OpenGL window wasn't already created by Togl_ObjConfigure() we
     * create it now.  We can tell by checking if the OpenGL context has
     * been initialized.
     */
    if (!togl->Ctx) {
        Tk_MakeWindowExist(togl->TkWin);
        if (togl->badWindow) {
            goto error;
        }
    }
    if (!togl->Ctx) {
        goto error;
    }
    Togl_MakeCurrent(togl);
    if (togl->contextTag == 0) {
        togl->contextTag = ++tpg->nextContextTag;
    }

    /* TCL3D Begin add SwapInterval */
    (void) Togl_SwapInterval(togl, togl->SwapInterval);
    /* TCL3D End add SwapInterval */

    /* If defined, call create callback */
    if (togl->CreateProc) {
        if (Togl_CallCallback(togl, togl->CreateProc) != TCL_OK) {
            goto error;
        }
    }

#ifdef TOGL_AGL
    SetMacBufRect(togl);
#endif
    /* If defined, call reshape proc */
    if (togl->ReshapeProc) {
        if (Togl_CallCallback(togl, togl->ReshapeProc) != TCL_OK) {
            goto error;
        }
    } else {
        glViewport(0, 0, togl->Width, togl->Height);
#if defined(TOGL_X11)
        if (togl->OverlayFlag) {
            Togl_UseLayer(togl, TOGL_OVERLAY);
            glViewport(0, 0, togl->Width, togl->Height);
            Togl_UseLayer(togl, TOGL_NORMAL);
        }
#endif
    }

    Tcl_AppendResult(interp, Tk_PathName(tkwin), NULL);

    /* Add to linked list */
    AddToList(togl);

    return TCL_OK;

error:
    sprintf (sDebugBuf, "ObjCmd Error: %s\n", Tcl_GetStringResult(interp)); DEBUG_OUT (sDebugBuf);
    state = Tcl_SaveInterpState(interp, TCL_OK);
    togl->badWindow = True;
    (void) Tcl_DeleteCommandFromToken(interp, togl->widgetCmd);
    Tcl_RestoreInterpState(interp, state);
    return TCL_ERROR;
}

#if TOGL_USE_OVERLAY

/*
 * Do all the setup for overlay planes
 * Return:   TCL_OK or TCL_ERROR
 */
static int
SetupOverlay(Togl *togl)
{
#  if defined(TOGL_X11)

#    ifdef GLX_TRANSPARENT_TYPE_EXT
    static int ovAttributeList[] = {
        GLX_BUFFER_SIZE, 2,
        GLX_LEVEL, 1,
        GLX_TRANSPARENT_TYPE_EXT, GLX_TRANSPARENT_INDEX_EXT,
        None
    };
#    else
    static int ovAttributeList[] = {
        GLX_BUFFER_SIZE, 2,
        GLX_LEVEL, 1,
        None
    };
#    endif

    XVisualInfo *visinfo;
    TkWindow *winPtr = (TkWindow *) togl->TkWin;

    XSetWindowAttributes swa;
    Tcl_HashEntry *hPtr;
    int     new_flag;

    visinfo = glXChooseVisual(togl->display, Tk_ScreenNumber(winPtr), ovAttributeList);
    if (!visinfo) {
        Tcl_AppendResult(togl->Interp, Tk_PathName(winPtr),
                ": No suitable overlay index visual available", NULL);
        togl->OverlayCtx = NULL;
        togl->OverlayWindow = 0;
        togl->OverlayCmap = 0;
        return TCL_ERROR;
    }
#    ifdef GLX_TRANSPARENT_INDEX_EXT
    {
        int     fail = glXGetConfig(togl->display, visinfo,
                GLX_TRANSPARENT_INDEX_VALUE_EXT,
                &togl->OverlayTransparentPixel);

        if (fail)
            togl->OverlayTransparentPixel = 0;  /* maybe, maybe ... */
    }
#    else
    togl->OverlayTransparentPixel = 0;  /* maybe, maybe ... */
#    endif

    /* share display lists with normal layer context */
    togl->OverlayCtx = glXCreateContext(togl->display, visinfo, togl->Ctx,
                       !togl->Indirect);

    swa.colormap = XCreateColormap(togl->display,
                   XRootWindow(togl->display, visinfo->screen),
                   visinfo->visual, AllocNone);
    togl->OverlayCmap = swa.colormap;

    swa.border_pixel = 0;
    swa.event_mask = ALL_EVENTS_MASK;
    togl->OverlayWindow = XCreateWindow(togl->display,
            Tk_WindowId(togl->TkWin), 0, 0,
            togl->Width, togl->Height, 0,
            visinfo->depth, InputOutput,
            visinfo->visual, CWBorderPixel | CWColormap | CWEventMask, &swa);

    hPtr = Tcl_CreateHashEntry(&winPtr->dispPtr->winTable,
            (const char *) togl->OverlayWindow, &new_flag);
    Tcl_SetHashValue(hPtr, winPtr);

    togl->OverlayIsMapped = False;

    /* Make sure window manager installs our colormap */
    XSetWMColormapWindows(togl->display, togl->OverlayWindow,
                          &togl->OverlayWindow, 1);

    return TCL_OK;

#  elif defined(TOGL_WGL) || defined(TOGL_AGL)
    /* not yet implemented on these */
    return TCL_ERROR;
#  endif
}

#endif /* TOGL_USE_OVERLAY */



#ifdef TOGL_WGL
static Bool ToglClassInitialized = False;

static LRESULT CALLBACK
Win32WinProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    LONG    result;
    Togl   *togl = (Togl *) GetWindowLongPtr(hwnd, 0);
    WNDCLASS childClass;

    switch (message) {
      case WM_WINDOWPOSCHANGED:
          /* Should be processed by DefWindowProc, otherwise a double buffered
           * context is not properly resized when the corresponding window is
           * resized. */
          break;
      case WM_DESTROY:
          if (togl && togl->TkWin != NULL) {
              if (togl->SetGrid > 0) {
                  Tk_UnsetGrid(togl->TkWin);
              }
              (void) Tcl_DeleteCommandFromToken(togl->Interp, togl->widgetCmd);
          }
          break;

      case WM_ERASEBKGND:
          /* We clear our own window */
          return 1;

      default:
#  if USE_STATIC_LIB
          return TkWinChildProc(hwnd, message, wParam, lParam);
#  else
          /*
           * OK, since TkWinChildProc is not explicitly exported in the
           * dynamic libraries, we have to retrieve it from the class info
           * registered with windows.
           *
           */
          if (tkWinChildProc == NULL) {
              GetClassInfo(Tk_GetHINSTANCE(), TK_WIN_CHILD_CLASS_NAME,
                      &childClass);
              tkWinChildProc = childClass.lpfnWndProc;
          }
          return tkWinChildProc(hwnd, message, wParam, lParam);
#  endif
    }
    result = DefWindowProc(hwnd, message, wParam, lParam);
    Tcl_ServiceAll();
    return result;
}
#endif /* TOGL_WGL */


/*
 * Togl_MakeWindow
 *
 *   Window creation function, invoked as a callback from Tk_MakeWindowExist.
 *   This is called instead of Tk_MakeWindow and must always succeed.
 */
static Window
Togl_MakeWindow(Tk_Window tkwin, Window parent, ClientData instanceData)
{
    Togl   *togl = (Togl *) instanceData;
    XVisualInfo *visinfo = NULL;
    Display *dpy;
    Colormap cmap;
    int     scrnum;
    Window  window = None;
    Bool contextAttribSuccess = False;

#if defined(TOGL_X11)
    Bool    directCtx = True;
    int     attrib_list[1000];
    int     attrib_count;
    int     dummy;
    XSetWindowAttributes swa;

#  define MAX_ATTEMPTS 12
    static int ci_depths[MAX_ATTEMPTS] = {
        8, 4, 2, 1, 12, 16, 8, 4, 2, 1, 12, 16
    };
    static int dbl_flags[MAX_ATTEMPTS] = {
        0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1
    };
    typedef GLXContext (*CREATECONTEXTATTRIBSARB)(Display *dpy,
            GLXFBConfig config, GLXContext share_context, Bool direct,
            const int *attrib_list);
    static CREATECONTEXTATTRIBSARB createContextAttribs;
    GLXFBConfig* fbc = None;
    int fbcCount = 0;
#elif defined(TOGL_WGL)
    typedef BOOL (WINAPI *ChooseFunc)(HDC, const int *, const FLOAT *, UINT, int *, UINT *);
    typedef const char *(WINAPI *StrFuncHDC)(HDC);
    typedef HGLRC (WINAPI * CREATECONTEXTATTRIBSARB) (HDC, HGLRC, const int *);
    static ChooseFunc choosePixelFormat = NULL;
    static CREATECONTEXTATTRIBSARB createContextAttribs;
    static BOOL initialized = False;
    HWND    hwnd, parentWin;
    int     pixelformat;
    HINSTANCE hInstance;
    WNDCLASS ToglClass;
    PIXELFORMATDESCRIPTOR pfd;
    UINT    matching;
    float   fAttribs[1] = { 0 };
    int     iAttribs[100];
    int     attrib_count;
    HGLRC hglrcAttrib = NULL;
#elif defined(TOGL_AGL)
    GLint   attribs[20];
    int     na;
    AGLPixelFormat fmt;
#endif /* TOGL_X11 */

    if (togl->badWindow) {
#if TK_MAJOR_VERSION >= 9
        Tk_Window winPtr = (Tk_Window ) tkwin;
        return Tk_MakeWindow(winPtr, parent);
#else
        TkWindow *winPtr = (TkWindow *) tkwin;
        return TkpMakeWindow(winPtr, parent);
#endif
    }

    dpy = Tk_Display(tkwin);

#if defined(TOGL_X11)
    /* Make sure OpenGL's GLX extension supported */
    if (!glXQueryExtension(dpy, &dummy, &dummy)) {
        Tcl_SetResult(togl->Interp,
                "X server has no OpenGL GLX extension", TCL_STATIC);
        return DUMMY_WINDOW;
    }

    if (togl->ShareContext && FindTogl(togl, togl->ShareContext)) {
        /* share OpenGL context with existing Togl widget */
        Togl   *shareWith = FindTogl(togl, togl->ShareContext);

        assert(shareWith != NULL);
        assert(shareWith->Ctx != NULL);
        togl->Ctx = shareWith->Ctx;
        togl->contextTag = shareWith->contextTag;
        togl->VisInfo = shareWith->VisInfo;
        visinfo = togl->VisInfo;
    } else {
        if (togl->PixelFormat) {
            XVisualInfo template;
            int     count = 1;

            template.visualid = togl->PixelFormat;
            visinfo = XGetVisualInfo(dpy, VisualIDMask, &template, &count);
            if (visinfo == NULL) {
                Tcl_SetResult(togl->Interp,
                        "Could not choose pixel format", TCL_STATIC);
                return DUMMY_WINDOW;
            }
            /* fill in flags normally passed in that affect behavior */
            (void) glXGetConfig(dpy, visinfo, GLX_RGBA, &togl->RgbaFlag);
            (void) glXGetConfig(dpy, visinfo, GLX_DOUBLEBUFFER,
                    &togl->DoubleFlag);
        } else {
            int attempt;

            /* It may take a few tries to get a visual */
            for (attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
                attrib_count = 0;
                if (togl->RgbaFlag) {
                    /* RGB[A] mode */
                    attrib_list[attrib_count++] = GLX_RED_SIZE;
                    attrib_list[attrib_count++] = togl->RgbaRed;
                    attrib_list[attrib_count++] = GLX_GREEN_SIZE;
                    attrib_list[attrib_count++] = togl->RgbaGreen;
                    attrib_list[attrib_count++] = GLX_BLUE_SIZE;
                    attrib_list[attrib_count++] = togl->RgbaBlue;
                    if (togl->AlphaFlag) {
                        attrib_list[attrib_count++] = GLX_ALPHA_SIZE;
                        attrib_list[attrib_count++] = togl->AlphaSize;
                    }

                    if (togl->RedMap)
                        free(togl->RedMap);
                    if (togl->GreenMap)
                        free(togl->GreenMap);
                    if (togl->BlueMap)
                        free(togl->BlueMap);
                    togl->RedMap = togl->GreenMap = togl->BlueMap =
                            NULL;
                    togl->MapSize = 0;
                } else {
                    /* Color index mode */
                    int depth;

                    attrib_list[attrib_count++] = GLX_BUFFER_SIZE;
                    depth = ci_depths[attempt];
                    attrib_list[attrib_count++] = depth;
                }
                if (togl->DepthFlag) {
                    attrib_list[attrib_count++] = GLX_DEPTH_SIZE;
                    attrib_list[attrib_count++] = togl->DepthSize;
                }
                if (togl->DoubleFlag || dbl_flags[attempt]) {
                    attrib_list[attrib_count++] = GLX_DOUBLEBUFFER;
                    attrib_list[attrib_count++] = togl->DoubleFlag;
                }
                if (togl->StencilFlag) {
                    attrib_list[attrib_count++] = GLX_STENCIL_SIZE;
                    attrib_list[attrib_count++] = togl->StencilSize;
                }
                if (togl->AccumFlag) {
                    attrib_list[attrib_count++] = GLX_ACCUM_RED_SIZE;
                    attrib_list[attrib_count++] = togl->AccumRed;
                    attrib_list[attrib_count++] = GLX_ACCUM_GREEN_SIZE;
                    attrib_list[attrib_count++] = togl->AccumGreen;
                    attrib_list[attrib_count++] = GLX_ACCUM_BLUE_SIZE;
                    attrib_list[attrib_count++] = togl->AccumBlue;
                    if (togl->AlphaFlag) {
                        attrib_list[attrib_count++] = GLX_ACCUM_ALPHA_SIZE;
                        attrib_list[attrib_count++] = togl->AccumAlpha;
                    }
                }

                /* TCL3D Begin add Multisample */
                if (togl->MultisampleBuffers) {
                    attrib_list[attrib_count++] = GLX_SAMPLE_BUFFERS_ARB;
                    attrib_list[attrib_count++] = togl->MultisampleBuffers;
                    attrib_list[attrib_count++] = GLX_SAMPLES_ARB;
                    attrib_list[attrib_count++] = togl->MultisampleSamples;
                }
                /* TCL3D End add Multisample */

                if (togl->AuxNumber != 0) {
                    attrib_list[attrib_count++] = GLX_AUX_BUFFERS;
                    attrib_list[attrib_count++] = togl->AuxNumber;
                }
                if (togl->Indirect) {
                    directCtx = False;
                }
                attrib_list[attrib_count++] = None;

                if (!(fbc = glXChooseFBConfig( dpy, Tk_ScreenNumber(tkwin), attrib_list, &fbcCount ))) {
                    continue;
                }
                visinfo = glXGetVisualFromFBConfig( dpy, fbc[ 0 ] );
                if (visinfo) {
                    /* found a GLX visual! */
                    break;
                }
            }

            togl->VisInfo = visinfo;

            if (visinfo == NULL) {
                Tcl_SetResult(togl->Interp,
                        "Could not choose pixel format", TCL_STATIC);
                return DUMMY_WINDOW;
            }

            /* TCL3D Begin add ContextAttribs */
            createContextAttribs = (CREATECONTEXTATTRIBSARB)
                                   glXGetProcAddressARB((const GLubyte *)"glXCreateContextAttribsARB");
            contextAttribSuccess = createContextAttribs? True: False;
            if (createContextAttribs) {
                int attribList[] = {
                    GLX_CONTEXT_MAJOR_VERSION_ARB, 1,
                    GLX_CONTEXT_MINOR_VERSION_ARB, 0,
                    GLX_CONTEXT_FLAGS_ARB, 0,
                    GLX_CONTEXT_PROFILE_MASK_ARB, 0,
                    0
                };
                int contextFlags = 0;
                if (togl->DebugContext) {
                    contextFlags |= GLX_CONTEXT_DEBUG_BIT_ARB;
                }
                if (togl->ForwardCompatibleContext) {
                    contextFlags |= GLX_CONTEXT_FORWARD_COMPATIBLE_BIT_ARB;
                }

                attribList[1] = togl->MajorVersion;
                attribList[3] = togl->MinorVersion;
                attribList[5] = contextFlags;
                attribList[7] = togl->CoreProfile? GLX_CONTEXT_CORE_PROFILE_BIT_ARB:
                                                   GLX_CONTEXT_COMPATIBILITY_PROFILE_BIT_ARB;

                if (!(togl->Ctx = createContextAttribs(dpy, fbc[0], None, directCtx, attribList))) {
                    contextAttribSuccess = False;
                }
            }
            /* TCL3D End add ContextAttribs */
        }
        if (!contextAttribSuccess) {
            togl->Ctx = glXCreateContext(dpy, visinfo, None, directCtx);
        }

        /*
         * Create a new OpenGL rendering context.
         */
        if (togl->ShareList) {
            /* share display lists with existing togl widget */
            Togl   *shareWith = FindTogl(togl, togl->ShareList);
            GLXContext shareCtx;

            if (shareWith) {
                shareCtx = shareWith->Ctx;
                togl->contextTag = shareWith->contextTag;
            } else {
                shareCtx = None;
            }
            togl->Ctx = glXCreateContext(dpy, visinfo, shareCtx, directCtx);
        }

        if (togl->Ctx == NULL) {
            Tcl_SetResult(togl->Interp,
                    "Could not create rendering context",
                    TCL_STATIC);
            return DUMMY_WINDOW;
        }

    }
#endif /* TOGL_X11 */

#ifdef TOGL_WGL
    parentWin = Tk_GetHWND(parent);
    hInstance = Tk_GetHINSTANCE();
    if (!ToglClassInitialized) {
        ToglClassInitialized = True;
        ToglClass.style = CS_HREDRAW | CS_VREDRAW;
        ToglClass.cbClsExtra = 0;
        ToglClass.cbWndExtra = sizeof (LONG_PTR);       /* to save Togl* */
        ToglClass.hInstance = hInstance;
        ToglClass.hbrBackground = NULL;
        ToglClass.lpszMenuName = NULL;
        ToglClass.lpszClassName = TOGL_CLASS_NAME;
        ToglClass.lpfnWndProc = Win32WinProc;
        ToglClass.hIcon = NULL;
        ToglClass.hCursor = NULL;
        if (!RegisterClass(&ToglClass)) {
            Tcl_SetResult(togl->Interp,
                    "unable register Togl window class", TCL_STATIC);
            return DUMMY_WINDOW;
        }
    }

    hwnd = CreateWindow(TOGL_CLASS_NAME, NULL,
            WS_CHILD | WS_CLIPCHILDREN | WS_CLIPSIBLINGS, 0, 0,
            togl->Width, togl->Height, parentWin, NULL, hInstance, NULL);
    SetWindowPos(hwnd, HWND_TOP, 0, 0, 0, 0,
            SWP_NOACTIVATE | SWP_NOMOVE | SWP_NOSIZE);

    togl->tglGLHdc = GetDC(hwnd);

    pfd.nSize = sizeof (PIXELFORMATDESCRIPTOR);
    pfd.nVersion = 1;
    pfd.dwFlags =
            PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_SUPPORT_COMPOSITION;
    if (togl->DoubleFlag) {
        pfd.dwFlags |= PFD_DOUBLEBUFFER;
    }

    if (togl->PixelFormat) {
        pixelformat = togl->PixelFormat;
    } else {
        pfd.cColorBits = togl->RgbaRed + togl->RgbaGreen + togl->RgbaBlue;
        pfd.iPixelType = togl->RgbaFlag ? PFD_TYPE_RGBA : PFD_TYPE_COLORINDEX;
        /* Alpha bitplanes are not supported in the current generic OpenGL
         * implementation, but may be supported by specific hardware devices. */
        pfd.cAlphaBits = togl->AlphaFlag ? togl->AlphaSize : 0;
        pfd.cAccumBits = togl->AccumFlag ? (togl->AccumRed + togl->AccumGreen +
                togl->AccumBlue + togl->AccumAlpha) : 0;
        pfd.cDepthBits = togl->DepthFlag ? togl->DepthSize : 0;
        pfd.cStencilBits = togl->StencilFlag ? togl->StencilSize : 0;
        /* Auxiliary buffers are not supported in the current generic OpenGL
         * implementation, but may be supported by specific hardware devices. */
        pfd.cAuxBuffers = togl->AuxNumber;
        pfd.iLayerType = PFD_MAIN_PLANE;

        /* TCL3D Begin add Multisample */
        attrib_count = 0;
        iAttribs[attrib_count++] = WGL_DRAW_TO_WINDOW_ARB;
        iAttribs[attrib_count++] = GL_TRUE;
        iAttribs[attrib_count++] = WGL_ACCELERATION_ARB;
        iAttribs[attrib_count++] = WGL_FULL_ACCELERATION_ARB;
        if (togl->RgbaFlag) {
            /* RGB[A] mode */
            iAttribs[attrib_count++] = WGL_RED_BITS_ARB;
            iAttribs[attrib_count++] = togl->RgbaRed;
            iAttribs[attrib_count++] = WGL_GREEN_BITS_ARB;
            iAttribs[attrib_count++] = togl->RgbaGreen;
            iAttribs[attrib_count++] = WGL_BLUE_BITS_ARB;
            iAttribs[attrib_count++] = togl->RgbaBlue;
            if (togl->AlphaFlag) {
                iAttribs[attrib_count++] = WGL_ALPHA_BITS_ARB;
                iAttribs[attrib_count++] = togl->AlphaSize;
            }
        } else {
            /* Color index mode */
            /* OPA TODO */
        }
        if (togl->DepthFlag) {
            iAttribs[attrib_count++] = WGL_DEPTH_BITS_ARB;
            iAttribs[attrib_count++] = togl->DepthSize;
        }
        if (togl->DoubleFlag) {
            iAttribs[attrib_count++] = WGL_DOUBLE_BUFFER_ARB;
            iAttribs[attrib_count++] = GL_TRUE;
        }
        if (togl->StencilFlag) {
            iAttribs[attrib_count++] = WGL_STENCIL_BITS_ARB;
            iAttribs[attrib_count++] = togl->StencilSize;
        }
        if (togl->AccumFlag) {
            iAttribs[attrib_count++] = WGL_ACCUM_RED_BITS_ARB;
            iAttribs[attrib_count++] = togl->AccumRed;
            iAttribs[attrib_count++] = WGL_ACCUM_GREEN_BITS_ARB;
            iAttribs[attrib_count++] = togl->AccumGreen;
            iAttribs[attrib_count++] = WGL_ACCUM_BLUE_BITS_ARB;
            iAttribs[attrib_count++] = togl->AccumBlue;
            if (togl->AlphaFlag) {
                iAttribs[attrib_count++] = WGL_ACCUM_ALPHA_BITS_ARB;
                iAttribs[attrib_count++] = togl->AccumAlpha;
            }
        }
        if (togl->MultisampleBuffers) {
            iAttribs[attrib_count++] = WGL_SAMPLE_BUFFERS_ARB;
            iAttribs[attrib_count++] = togl->MultisampleBuffers;
            iAttribs[attrib_count++] = WGL_SAMPLES_ARB;
            iAttribs[attrib_count++] = togl->MultisampleSamples;
        }
        if (togl->AuxNumber != 0) {
            iAttribs[attrib_count++] = WGL_AUX_BUFFERS_ARB;
            iAttribs[attrib_count++] = togl->AuxNumber;
        }
        iAttribs[attrib_count] = 0;

        if (!initialized) {
            const char *extensions;
            StrFuncHDC getExtensionsString;
            HWND hwnd1;
            HDC hdc;
            HGLRC hglrc;
            int pformat;

            hwnd1 = CreateWindow(TOGL_CLASS_NAME, TEXT(""), WS_POPUP | WS_DISABLED,
                                0, 0, 10, 10,
                                NULL, NULL, hInstance, NULL);
            hdc = GetDC(hwnd1);
            pformat = ChoosePixelFormat(hdc, &pfd);
            SetPixelFormat(hdc, pformat, &pfd);

            hglrc = wglCreateContext(hdc);
            if ( hglrc ) {
                wglMakeCurrent(hdc, hglrc);
            }

            getExtensionsString = (StrFuncHDC)
                                wglGetProcAddress("wglGetExtensionsStringARB");
            if (getExtensionsString == NULL) {
                getExtensionsString = (StrFuncHDC)
                                wglGetProcAddress("wglGetExtensionsStringEXT");
            }
            if (getExtensionsString) {
                extensions = getExtensionsString(togl->tglGLHdc);
                if (strstr(extensions, "WGL_ARB_pixel_format") != NULL) {
                    choosePixelFormat = (ChooseFunc) wglGetProcAddress("wglChoosePixelFormatARB");
                }
            }

            createContextAttribs = (CREATECONTEXTATTRIBSARB)
                            wglGetProcAddress("wglCreateContextAttribsARB");

            initialized = True;
            if ( hglrc ) {
                wglMakeCurrent(NULL, NULL);
                wglDeleteContext(hglrc);
            }
            ReleaseDC(hwnd1, hdc);
            /* OPA TODO
            DestroyWindow(hwnd1);
            */
        }

        /* Choose and set the closest available pixel format */
        if (!choosePixelFormat ||
            !choosePixelFormat(togl->tglGLHdc, iAttribs, fAttribs, 1, &pixelformat, &matching) ||
            !matching) {
            pixelformat = ChoosePixelFormat(togl->tglGLHdc, &pfd);
        }
        if (!pixelformat) {
            Tcl_SetResult(togl->Interp,
                    "Togl: Could not choose pixel format",
                    TCL_STATIC);
            return DUMMY_WINDOW;
        }
        /* TCL3D End add Multisample */
    }
    if (SetPixelFormat(togl->tglGLHdc, pixelformat, &pfd) == FALSE) {
        Tcl_SetResult(togl->Interp,
                "Could not choose pixel format", TCL_STATIC);
        ReleaseDC(hwnd, togl->tglGLHdc);
        togl->tglGLHdc = NULL;
        DestroyWindow(hwnd);
        return DUMMY_WINDOW;
    }

    /* Get the actual pixel format */
    DescribePixelFormat(togl->tglGLHdc, pixelformat, sizeof (pfd), &pfd);
    if (togl->PixelFormat) {
        /* fill in flags normally passed in that affect behavior */
        togl->RgbaFlag = pfd.iPixelType == PFD_TYPE_RGBA;
        togl->DoubleFlag = pfd.cDepthBits > 0;
        /* TODO: set depth flag, and more */
    }

    if (togl->ShareContext && FindTogl(togl, togl->ShareContext)) {
        /* share OpenGL context with existing Togl widget */
        Togl   *shareWith = FindTogl(togl, togl->ShareContext);

        assert(shareWith);
        assert(shareWith->Ctx);
        togl->Ctx = shareWith->Ctx;
        togl->contextTag = shareWith->contextTag;
        togl->VisInfo = shareWith->VisInfo;
        visinfo = togl->VisInfo;
    } else {
        /*
         * Create a new OpenGL rendering context. And check to share lists.
         */

        /* TCL3D Begin add ContextAttribs */
        contextAttribSuccess = createContextAttribs? True: False;
        if (createContextAttribs) {
            int attribList[] = {
                WGL_CONTEXT_MAJOR_VERSION_ARB, 1,
                WGL_CONTEXT_MINOR_VERSION_ARB, 0,
                WGL_CONTEXT_FLAGS_ARB, 0,
                WGL_CONTEXT_PROFILE_MASK_ARB, 0,
                0
            };
            int contextFlags = 0;
            if (togl->DebugContext) {
                contextFlags |= WGL_CONTEXT_DEBUG_BIT_ARB;
            }
            if (togl->ForwardCompatibleContext) {
                contextFlags |= WGL_CONTEXT_FORWARD_COMPATIBLE_BIT_ARB;
            }

            attribList[1] = togl->MajorVersion;
            attribList[3] = togl->MinorVersion;
            attribList[5] = contextFlags;
            attribList[7] = togl->CoreProfile? WGL_CONTEXT_CORE_PROFILE_BIT_ARB:
                                               WGL_CONTEXT_COMPATIBILITY_PROFILE_BIT_ARB;

            if (!(hglrcAttrib = createContextAttribs(togl->tglGLHdc, 0, attribList))) {
                contextAttribSuccess = False;
            }
            wglMakeCurrent(togl->tglGLHdc, hglrcAttrib);
            togl->Ctx = hglrcAttrib;
            /* TCL3D End add ContextAttribs */
        }
        if (!contextAttribSuccess) {
            togl->Ctx = wglCreateContext(togl->tglGLHdc);
        }

        if (togl->ShareList) {
            /* share display lists with existing togl widget */
            Togl   *shareWith = FindTogl(togl, togl->ShareList);

            if (shareWith) {
                if (!wglShareLists(shareWith->Ctx, togl->Ctx)) {
                    Tcl_SetResult(togl->Interp,
                            "unable to share display lists",
                            TCL_STATIC);
                    ReleaseDC(hwnd, togl->tglGLHdc);
                    togl->tglGLHdc = NULL;
                    DestroyWindow(hwnd);
                    return DUMMY_WINDOW;
                }
                togl->contextTag = shareWith->contextTag;
            }
        }

        if (!togl->Ctx) {
            Tcl_SetResult(togl->Interp,
                    "Could not create rendering context",
                    TCL_STATIC);
            ReleaseDC(hwnd, togl->tglGLHdc);
            togl->tglGLHdc = NULL;
            DestroyWindow(hwnd);
            return DUMMY_WINDOW;
        }

        /* Just for portability, define the simplest visinfo */
        visinfo = (XVisualInfo *) malloc(sizeof (XVisualInfo));
        visinfo->visual = DefaultVisual(dpy, DefaultScreen(dpy));
        visinfo->depth = visinfo->visual->bits_per_rgb;
        togl->VisInfo = visinfo;
    }
    SetWindowLongPtr(hwnd, 0, (LONG_PTR) togl);
#endif /* TOGL_WGL */

    /* for TOGL_AGL, we create the window first, then choose the pixel format */

    /*
     * find a colormap
     */
    scrnum = Tk_ScreenNumber(tkwin);
    if (togl->RgbaFlag) {
        /* Colormap for RGB mode */
#if defined(TOGL_X11)
        cmap = get_rgb_colormap(dpy, scrnum, visinfo, tkwin);

#elif defined(TOGL_WGL)
        if (pfd.dwFlags & PFD_NEED_PALETTE) {
            cmap = Win32CreateRgbColormap(pfd);
        } else {
            cmap = DefaultColormap(dpy, scrnum);
        }

        if (togl->RedMap)
            free(togl->RedMap);
        if (togl->GreenMap)
            free(togl->GreenMap);
        if (togl->BlueMap)
            free(togl->BlueMap);
        togl->RedMap = togl->GreenMap = togl->BlueMap = NULL;
        togl->MapSize = 0;

#elif defined(TOGL_AGL)
        cmap = DefaultColormap(dpy, scrnum);

        if (togl->RedMap)
            free(togl->RedMap);
        if (togl->GreenMap)
            free(togl->GreenMap);
        if (togl->BlueMap)
            free(togl->BlueMap);
        togl->RedMap = togl->GreenMap = togl->BlueMap = NULL;
        togl->MapSize = 0;
#endif
    } else {
        /* Colormap for CI mode */
#ifdef TOGL_WGL
        /* this logic is to overcome a combination driver/compiler bug: 1.
         * cColorBits may be unusally large (e.g., 32 instead of 8 or 12) 2. 1
         * << 32 might be 1 instead of zero (gcc for ia32) */
        if (pfd.cColorBits >= MAX_CI_COLORMAP_BITS) {
            togl->CiColormapSize = MAX_CI_COLORMAP_SIZE;
        } else {
            togl->CiColormapSize = 1 << pfd.cColorBits;
            if (togl->CiColormapSize >= MAX_CI_COLORMAP_SIZE)
                togl->CiColormapSize = MAX_CI_COLORMAP_SIZE;
        }

#endif
        if (togl->PrivateCmapFlag) {
            /* need read/write colormap so user can store own color entries */
#if defined(TOGL_X11)
            cmap = XCreateColormap(dpy, XRootWindow(dpy, visinfo->screen),
                    visinfo->visual, AllocAll);
#elif defined(TOGL_WGL)
            cmap = Win32CreateCiColormap(togl);
#elif defined(TOGL_AGL)
            /* need to figure out how to do this correctly on Mac... */
            cmap = DefaultColormap(dpy, scrnum);
#endif
        } else {
            if (visinfo->visual == DefaultVisual(dpy, scrnum)) {
                /* share default/root colormap */
                cmap = Tk_Colormap(tkwin);
            } else {
                /* make a new read-only colormap */
                cmap = XCreateColormap(dpy, XRootWindow(dpy, visinfo->screen),
                        visinfo->visual, AllocNone);
            }
        }
    }

#if !defined(TOGL_AGL)
    /* Make sure Tk knows to switch to the new colormap when the cursor is over
     * this window when running in color index mode. */
    (void) Tk_SetWindowVisual(tkwin, visinfo->visual, visinfo->depth, cmap);
#endif

#ifdef TOGL_WGL
    /* Install the colormap */
    SelectPalette(togl->tglGLHdc, ((TkWinColormap *) cmap)->palette, TRUE);
    RealizePalette(togl->tglGLHdc);
#endif

#if defined(TOGL_X11)
    swa.background_pixmap = None;
    swa.border_pixel = 0;
    swa.colormap = cmap;
    swa.event_mask = ALL_EVENTS_MASK;
    window = XCreateWindow(dpy, parent,
            0, 0, togl->Width, togl->Height,
            0, visinfo->depth,
            InputOutput, visinfo->visual,
            CWBackPixmap | CWBorderPixel | CWColormap | CWEventMask, &swa);
    /* Make sure window manager installs our colormap */
    (void) XSetWMColormapWindows(dpy, window, &window, 1);

   if (!togl->DoubleFlag) {
        int     dbl_flag;

        /* See if we requested single buffering but had to accept a double
         * buffered visual.  If so, set the GL draw buffer to be the front
         * buffer to simulate single buffering. */
        if (glXGetConfig(dpy, togl->VisInfo, GLX_DOUBLEBUFFER, &dbl_flag)) {
            if (dbl_flag) {
                glXMakeCurrent(dpy, window, togl->Ctx);
                glDrawBuffer(GL_FRONT);
                glReadBuffer(GL_FRONT);
            }
        }
    }
#elif defined(TOGL_WGL)
    window = Tk_AttachHWND(tkwin, hwnd);

#elif defined(TOGL_AGL)
    {

#if TK_MAJOR_VERSION >= 9
        Tk_Window winPtr = (Tk_Window ) tkwin;
        window = Tk_MakeWindow(winPtr, parent);
#else
        TkWindow *winPtr = (TkWindow *) tkwin;
        window = TkpMakeWindow(winPtr, parent);
#endif
    }
#endif

#if TOGL_USE_OVERLAY
    if (togl->OverlayFlag) {
        if (SetupOverlay(togl) == TCL_ERROR) {
            fprintf(stderr, "Warning: Could not setup overlay.\n");
            togl->OverlayFlag = False;
        }
    }
#endif

    /* Request the X window to be displayed */
    (void) XMapWindow(dpy, window);

#if defined(TOGL_AGL)
    if (togl->ShareContext && FindTogl(togl, togl->ShareContext)) {
        /* share OpenGL context with existing Togl widget */
        Togl   *shareWith = FindTogl(togl, togl->ShareContext);

        assert(shareWith);
        assert(shareWith->Ctx);
        togl->Ctx = shareWith->Ctx;
        togl->contextTag = shareWith->contextTag;
        togl->VisInfo = shareWith->VisInfo;
        visinfo = togl->VisInfo;

    } else {
        AGLContext shareCtx = NULL;

        if (togl->PixelFormat) {
            /* fill in RgbaFlag, DoubleFlag */
            fmt = (AGLPixelFormat) togl->PixelFormat;
            GLint   has_rgba, has_doublebuf;

            if (aglDescribePixelFormat(fmt, AGL_RGBA, &has_rgba)
                    && aglDescribePixelFormat(fmt, AGL_DOUBLEBUFFER,
                            &has_doublebuf)) {
                togl->RgbaFlag = (has_rgba ? True : False);
                togl->DoubleFlag = (has_doublebuf ? True : False);
            } else {
                Tcl_SetResult(togl->Interp,
                        "failed querying pixel format attributes",
                        TCL_STATIC);
                return DUMMY_WINDOW;
            }
        } else {

            /* Need to do this after mapping window, so MacDrawable structure
             * is more completely filled in */
            na = 0;
            attribs[na++] = AGL_MINIMUM_POLICY;
            /* ask for hardware-accelerated onscreen */
            attribs[na++] = AGL_ACCELERATED;
            attribs[na++] = AGL_NO_RECOVERY;
            if (togl->RgbaFlag) {
                /* RGB[A] mode */
                attribs[na++] = AGL_RGBA;
                attribs[na++] = AGL_RED_SIZE;
                attribs[na++] = togl->RgbaRed;
                attribs[na++] = AGL_GREEN_SIZE;
                attribs[na++] = togl->RgbaGreen;
                attribs[na++] = AGL_BLUE_SIZE;
                attribs[na++] = togl->RgbaBlue;
                if (togl->AlphaFlag) {
                    attribs[na++] = AGL_ALPHA_SIZE;
                    attribs[na++] = togl->AlphaSize;
                }
            } else {
                /* Color index mode */
                attribs[na++] = AGL_BUFFER_SIZE;
                attribs[na++] = 8;
            }
            if (togl->DepthFlag) {
                attribs[na++] = AGL_DEPTH_SIZE;
                attribs[na++] = togl->DepthSize;
            }
            if (togl->DoubleFlag) {
                attribs[na++] = AGL_DOUBLEBUFFER;
            }
            if (togl->StencilFlag) {
                attribs[na++] = AGL_STENCIL_SIZE;
                attribs[na++] = togl->StencilSize;
            }
            if (togl->AccumFlag) {
                attribs[na++] = AGL_ACCUM_RED_SIZE;
                attribs[na++] = togl->AccumRed;
                attribs[na++] = AGL_ACCUM_GREEN_SIZE;
                attribs[na++] = togl->AccumGreen;
                attribs[na++] = AGL_ACCUM_BLUE_SIZE;
                attribs[na++] = togl->AccumBlue;
                if (togl->AlphaFlag) {
                    attribs[na++] = AGL_ACCUM_ALPHA_SIZE;
                    attribs[na++] = togl->AccumAlpha;
                }
            }
            if (togl->AuxNumber != 0) {
                attribs[na++] = AGL_AUX_BUFFERS;
                attribs[na++] = togl->AuxNumber;
            }
            attribs[na++] = AGL_NONE;

            if ((fmt = aglChoosePixelFormat(NULL, 0, attribs)) == NULL) {
                Tcl_SetResult(togl->Interp,
                        "Could not choose pixel format", TCL_STATIC);
                return DUMMY_WINDOW;
            }
        }

        /*
         * Check whether to share lists.
         */
        if (togl->ShareList) {
            /* share display lists with existing togl widget */
            Togl   *shareWith = FindTogl(togl, togl->ShareList);

            if (shareWith) {
                shareCtx = shareWith->Ctx;
                togl->contextTag = shareWith->contextTag;
            }
        }
        if ((togl->Ctx = aglCreateContext(fmt, shareCtx)) == NULL) {
            GLenum  err = aglGetError();

            aglDestroyPixelFormat(fmt);
            if (err == AGL_BAD_MATCH)
                Tcl_SetResult(togl->Interp,
                        "unable to share display lists"
                        ": shared context doesn't match", TCL_STATIC);
            else if (err == AGL_BAD_CONTEXT)
                Tcl_SetResult(togl->Interp,
                        "unable to share display lists"
                        ": bad shared context", TCL_STATIC);
            else if (err == AGL_BAD_PIXELFMT)
                Tcl_SetResult(togl->Interp,
                        "Could not create rendering context"
                        ": bad pixel format", TCL_STATIC);
            else
                Tcl_SetResult(togl->Interp,
                        "Could not create rendering context"
                        ": unknown reason", TCL_STATIC);
            return DUMMY_WINDOW;
        }

        aglDestroyPixelFormat(fmt);
        if (!aglSetDrawable(togl->Ctx,
                        ((MacDrawable *) (window))->toplevel->grafPtr)) {
            aglDestroyContext(togl->Ctx);
            Tcl_SetResult(togl->Interp,
                    "Could not set drawable", TCL_STATIC);
            return DUMMY_WINDOW;
        }

        /* Just for portability, define the simplest visinfo */
        visinfo = (XVisualInfo *) malloc(sizeof (XVisualInfo));
        visinfo->visual = DefaultVisual(dpy, DefaultScreen(dpy));
        visinfo->depth = visinfo->visual->bits_per_rgb;

        Tk_SetWindowVisual(tkwin, visinfo->visual, visinfo->depth, cmap);
    }
#endif /* TOGL_AGL */

#if defined(TOGL_X11)
    /* Check for a single/double buffering snafu */
    {
        int     dbl_flag;

        if (glXGetConfig(dpy, visinfo, GLX_DOUBLEBUFFER, &dbl_flag)) {
            if (!togl->DoubleFlag && dbl_flag) {
                /* We requested single buffering but had to accept a */
                /* double buffered visual.  Set the GL draw buffer to */
                /* be the front buffer to simulate single buffering. */
                glDrawBuffer(GL_FRONT);
            }
        }
    }
#endif

    if (!togl->RgbaFlag) {
        int     index_size;

#if defined(TOGL_X11) || defined(TOGL_AGL)
        GLint   index_bits;

        glGetIntegerv(GL_INDEX_BITS, &index_bits);
        index_size = 1 << index_bits;
#elif defined(TOGL_WGL)
        index_size = togl->CiColormapSize;
#endif
        if (togl->MapSize != index_size) {
            if (togl->RedMap)
                free(togl->RedMap);
            if (togl->GreenMap)
                free(togl->GreenMap);
            if (togl->BlueMap)
                free(togl->BlueMap);
            togl->MapSize = index_size;
            togl->RedMap   = (GLfloat *) calloc(index_size, sizeof (GLfloat));
            togl->GreenMap = (GLfloat *) calloc(index_size, sizeof (GLfloat));
            togl->BlueMap  = (GLfloat *) calloc(index_size, sizeof (GLfloat));
        }
    }
    return window;
}

/*
 * Togl_WorldChanged
 *
 *    Add support for setgrid option.
 */
static void
Togl_WorldChanged(ClientData instanceData)
{
    Togl *togl = (Togl *) instanceData;

    Tk_GeometryRequest(togl->TkWin, togl->Width, togl->Height);
    Tk_SetInternalBorder(togl->TkWin, 0);
    if (togl->SetGrid > 0) {
        Tk_SetGrid(togl->TkWin, togl->Width / togl->SetGrid,
                togl->Height / togl->SetGrid, togl->SetGrid, togl->SetGrid);
    } else {
        Tk_UnsetGrid(togl->TkWin);
    }
}

/*
 * ToglFree
 *
 *      Wrap the ckfree macro.
 */
#if TCL_MAJOR_VERSION >= 9 || (TCL_MAJOR_VERSION == 8 && TCL_MINOR_VERSION > 7)
static void ToglFree(void *clientData)
#else
static void ToglFree(char *clientData)
#endif
{
    ckfree(clientData);
}

/*
 * ToglCmdDeletedProc
 *
 *      This procedure is invoked when a widget command is deleted.  If
 *      the widget isn't already in the process of being destroyed,
 *      this command destroys it.
 *
 * Results:
 *      None.
 *
 * Side effects:
 *      The widget is destroyed.
 *
 *----------------------------------------------------------------------
 */
static void
ToglCmdDeletedProc(ClientData clientData)
{
    Togl   *togl = (Togl *) clientData;
    Tk_Window tkwin = togl->TkWin;

    /*
     * This procedure could be invoked either because the window was
     * destroyed and the command was then deleted (in which case tkwin
     * is NULL) or because the command was deleted, and then this procedure
     * destroys the widget.
     */

    if (tkwin) {
        Tk_DeleteEventHandler(tkwin, ExposureMask | StructureNotifyMask,
                Togl_EventProc, (ClientData) togl);
    }

    Tcl_Preserve((ClientData) togl);
    Tcl_EventuallyFree((ClientData) togl, ToglFree);

    if (togl->DestroyProc) {
        /* call user's cleanup code */
        if (Togl_CallCallback(togl, togl->DestroyProc) != TCL_OK) {
            Tcl_BackgroundError (togl->Interp);
            return;
        }
    }

    if (togl->TimerProc != NULL) {
        Tcl_DeleteTimerHandler(togl->timerHandler);
        togl->timerHandler = NULL;
    }
    if (togl->UpdatePending) {
        Tcl_CancelIdleCall(Togl_Render, (ClientData) togl);
        togl->UpdatePending = False;
    }
#ifndef NO_TK_CURSOR
    if (togl->Cursor != NULL) {
        Tk_FreeCursor(togl->display, togl->Cursor);
        togl->Cursor = NULL;
    }
#endif

    /* remove from linked list */
    RemoveFromList(togl);

    togl->TkWin = NULL;
    if (tkwin != NULL && Tk_WindowId(tkwin) != DUMMY_WINDOW) {

#if defined(TOGL_AGL)
        if (togl->Ctx) {
            aglDestroyContext(togl->Ctx);
            togl->Ctx = NULL;
        }
#endif

#if defined(TOGL_X11)
        if (togl->Ctx) {
            if (FindToglWithSameContext(togl) == NULL)
                glXDestroyContext(togl->display, togl->Ctx);
            togl->Ctx = NULL;
        }
#  if TOGL_USE_OVERLAY
        if (togl->OverlayCtx) {
            Tcl_HashEntry *entryPtr;
            TkWindow *winPtr = (TkWindow *) tkwin;

            if (winPtr) {
                entryPtr = Tcl_FindHashEntry(&winPtr->dispPtr->winTable,
                        (const char *) togl->OverlayWindow);
                Tcl_DeleteHashEntry(entryPtr);
            }
            if (FindToglWithSameOverlayContext(togl) == NULL)
                glXDestroyContext(togl->display, togl->OverlayCtx);
            togl->OverlayCtx = NULL;
        }
#  endif /* TOGL_USE_OVERLAY */
#endif
#if defined(TOGL_WGL)
        if (togl->Ctx) {
            if (FindToglWithSameContext(togl) == NULL)
                wglDeleteContext(togl->Ctx);
            togl->Ctx = NULL;
        }
        if (tkwin && togl->tglGLHdc) {
            HWND    hwnd = Tk_GetHWND(Tk_WindowId(tkwin));

            ReleaseDC(hwnd, togl->tglGLHdc);
            togl->tglGLHdc = NULL;
        }
#endif
        /* TODO: delete contexts on other platforms */

        if (togl->SetGrid > 0) {
            Tk_UnsetGrid(tkwin);
        }
        Tk_DestroyWindow(tkwin);
    }

    Tcl_Release((ClientData) togl);
}

/*
 * This gets called to handle Togl window configuration events
 */
static void
Togl_EventProc(ClientData clientData, XEvent *eventPtr)
{
    Togl   *togl = (Togl *) clientData;
    switch (eventPtr->type) {
      case Expose:
          if (eventPtr->xexpose.count == 0) {
              if (!togl->UpdatePending
                      && eventPtr->xexpose.window == Tk_WindowId(togl->TkWin)) {
                  Togl_PostRedisplay(togl);
              }
#if defined(TOGL_X11)
              if (!togl->OverlayUpdatePending && togl->OverlayFlag
                      && togl->OverlayIsMapped
                      && eventPtr->xexpose.window == togl->OverlayWindow) {
                  Togl_PostOverlayRedisplay(togl);
              }
#endif
          }
          break;
      case ConfigureNotify:
         if (togl->Width == Tk_Width(togl->TkWin)
                          && togl->Height == Tk_Height(togl->TkWin)) {
#ifdef TOGL_AGL
              // Even though the size hasn't changed,
              // it's position on the screen may have.
              SetMacBufRect(togl);
#endif
              break;
          }
          togl->Width = Tk_Width(togl->TkWin);
          togl->Height = Tk_Height(togl->TkWin);
          (void) XResizeWindow(Tk_Display(togl->TkWin),
                  Tk_WindowId(togl->TkWin), togl->Width, togl->Height);
#if defined(TOGL_X11)
          if (togl->OverlayFlag) {
              (void) XResizeWindow(Tk_Display(togl->TkWin),
                      togl->OverlayWindow, togl->Width, togl->Height);
              (void) XRaiseWindow(Tk_Display(togl->TkWin), togl->OverlayWindow);
          }
#endif
#ifdef TOGL_AGL
          SetMacBufRect(togl);
#endif
          Togl_MakeCurrent(togl);
          if (togl->ReshapeProc) {
              if (Togl_CallCallback(togl, togl->ReshapeProc) != TCL_OK) {
                  Tcl_BackgroundError (togl->Interp);
                  return;
              }
          } else {
              glViewport(0, 0, togl->Width, togl->Height);
#if defined(TOGL_X11)
              if (togl->OverlayFlag) {
                      Togl_UseLayer(togl, TOGL_OVERLAY);
                      glViewport(0, 0, togl->Width, togl->Height);
                      Togl_UseLayer(togl, TOGL_NORMAL);
                  }
#endif
          }
          break;
      case MapNotify:
#if defined(TOGL_AGL)
      {
          /*
           * See comment for the UnmapNotify case below.
           */
          AGLDrawable d = Togl_MacOSXGetDrawablePort(togl);

          /* aglSetDrawable is deprecated in OS X 10.5 */
          aglSetDrawable(togl->Ctx, d);
          SetMacBufRect(togl);
      }
#endif
          break;
      case UnmapNotify:
#if defined(TOGL_AGL)
      {
          /*
           * For Mac OS X Aqua, Tk subwindows are not implemented as
           * separate Aqua windows.  They are just different regions of
           * a single Aqua window.  To unmap them they are just not drawn.
           * Have to disconnect the AGL context otherwise they will continue
           * to be displayed directly by Aqua.
           */
          aglSetDrawable(togl->Ctx, NULL);
      }
#endif
          break;
      case DestroyNotify:
          if (togl->TkWin != NULL) {
#ifdef TOGL_WGL
              HWND    hwnd = Tk_GetHWND(Tk_WindowId(togl->TkWin));

              /* Prevent Win32WinProc from calling Tcl_DeleteCommandFromToken
               * a second time */
              SetWindowLongPtr(hwnd, 0, (LONG_PTR) 0);
#endif
              if (togl->SetGrid > 0) {
                  Tk_UnsetGrid(togl->TkWin);
              }
              (void) Tcl_DeleteCommandFromToken(togl->Interp, togl->widgetCmd);
          }
          break;
      default:
          /* nothing */
          ;
    }
}

static void
Togl_PostRedisplay(Togl *togl)
{
    if (!togl->UpdatePending) {
        togl->UpdatePending = True;
        Tcl_DoWhenIdle(Togl_Render, (ClientData) togl);
    }
}

Bool
Togl_UpdatePending(const Togl *togl)
{
    return togl->UpdatePending;
}

static void
Togl_SwapBuffers(const Togl *togl)
{
    if (togl->DoubleFlag) {
#if defined(TOGL_WGL)
        int res = SwapBuffers(togl->tglGLHdc);
        if (!res) {
            Tcl_Panic("SwapBuffers");
        }
#elif defined(TOGL_X11)
        glXSwapBuffers(Tk_Display(togl->TkWin), Tk_WindowId(togl->TkWin));
#elif defined(TOGL_AGL)
        aglSwapBuffers(togl->Ctx);
#endif
    } else {
        glFlush();
    }
}

static int
Togl_Width(const Togl *togl)
{
    return togl->Width;
}

static int
Togl_Height(const Togl *togl)
{
    return togl->Height;
}

Tcl_Interp *
Togl_Interp(const Togl *togl)
{
    return togl->Interp;
}


Tk_Window
Togl_TkWin(const Togl *togl)
{
    return togl->TkWin;
}

const char *
Togl_CommandName(const Togl *togl)
{
    return Tcl_GetCommandName(togl->Interp, togl->widgetCmd);
}

static int
Togl_ContextTag(const Togl *togl)
{
    return togl->contextTag;
}


#if defined(TOGL_X11)
/*
 * A replacement for XAllocColor.  This function should never
 * fail to allocate a color.  When XAllocColor fails, we return
 * the nearest matching color.  If we have to allocate many colors
 * this function isn't too efficient; the XQueryColors() could be
 * done just once.
 * Written by Michael Pichler, Brian Paul, Mark Kilgard
 * Input:  dpy - X display
 *         cmap - X colormap
 *         cmapSize - size of colormap
 * In/Out: color - the XColor struct
 * Output:  exact - 1=exact color match, 0=closest match
 */
static void
noFaultXAllocColor(Display *dpy, Colormap cmap, int cmapSize,
        XColor *color, int *exact)
{
    XColor *ctable, subColor;
    int     i, bestmatch;
    double  mindist;            /* 3*2^16^2 exceeds long int precision. */

    /* First try just using XAllocColor. */
    if (XAllocColor(dpy, cmap, color)) {
        *exact = 1;
        return;
    }

    /* Retrieve color table entries. */
    /* XXX alloca candidate. */
    ctable = (XColor *) ckalloc(cmapSize * sizeof (XColor));
    for (i = 0; i < cmapSize; i++) {
        ctable[i].pixel = i;
    }
    (void) XQueryColors(dpy, cmap, ctable, cmapSize);

    /* Find best match. */
    bestmatch = -1;
    mindist = 0;
    for (i = 0; i < cmapSize; i++) {
        double  dr = (double) color->red - (double) ctable[i].red;
        double  dg = (double) color->green - (double) ctable[i].green;
        double  db = (double) color->blue - (double) ctable[i].blue;
        double  dist = dr * dr + dg * dg + db * db;

        if (bestmatch < 0 || dist < mindist) {
            bestmatch = i;
            mindist = dist;
        }
    }

    /* Return result. */
    subColor.red = ctable[bestmatch].red;
    subColor.green = ctable[bestmatch].green;
    subColor.blue = ctable[bestmatch].blue;
    ckfree((char *)ctable);
    /* Try to allocate the closest match color.  This should only fail if the
     * cell is read/write.  Otherwise, we're incrementing the cell's reference
     * count. */
    if (!XAllocColor(dpy, cmap, &subColor)) {
        /* do this to work around a problem reported by Frank Ortega */
        subColor.pixel = (unsigned long) bestmatch;
        subColor.red = ctable[bestmatch].red;
        subColor.green = ctable[bestmatch].green;
        subColor.blue = ctable[bestmatch].blue;
        subColor.flags = DoRed | DoGreen | DoBlue;
    }
    *color = subColor;
}

#elif defined(TOGL_WGL)

static UINT
Win32AllocColor(const Togl *togl, float red, float green, float blue)
{
    /* Modified version of XAllocColor emulation of Tk. - returns index,
     * instead of color itself - allocates logical palette entry even for
     * non-palette devices */

    TkWinColormap *cmap = (TkWinColormap *) Tk_Colormap(togl->TkWin);
    UINT    index;
    COLORREF newColor, closeColor;
    PALETTEENTRY entry, closeEntry;
    int    isNew;
    size_t refCount;
    Tcl_HashEntry *entryPtr;

    entry.peRed = (unsigned char) (red * 255 + .5);
    entry.peGreen = (unsigned char) (green * 255 + .5);
    entry.peBlue = (unsigned char) (blue * 255 + .5);
    entry.peFlags = 0;

    /*
     * Find the nearest existing palette entry.
     */

    newColor = RGB(entry.peRed, entry.peGreen, entry.peBlue);
    index = GetNearestPaletteIndex(cmap->palette, newColor);
    GetPaletteEntries(cmap->palette, index, 1, &closeEntry);
    closeColor = RGB(closeEntry.peRed, closeEntry.peGreen, closeEntry.peBlue);

    /*
     * If this is not a duplicate and colormap is not full, allocate a new entry.
     */

    if (newColor != closeColor) {
        if (cmap->size == (unsigned int) togl->CiColormapSize) {
            entry = closeEntry;
        } else {
            cmap->size++;
            ResizePalette(cmap->palette, cmap->size);
            index = cmap->size - 1;
            SetPaletteEntries(cmap->palette, index, 1, &entry);
            SelectPalette(togl->tglGLHdc, cmap->palette, TRUE);
            RealizePalette(togl->tglGLHdc);
        }
    }
    newColor = PALETTERGB(entry.peRed, entry.peGreen, entry.peBlue);
    entryPtr = Tcl_CreateHashEntry(&cmap->refCounts,
            (const char *) &newColor, &isNew);
    if (isNew) {
        refCount = 1;
    } else {
        refCount = ((size_t)Tcl_GetHashValue(entryPtr)) + 1;
    }
    Tcl_SetHashValue(entryPtr, (void *)refCount);

    togl->RedMap[index] = (GLfloat) (entry.peRed / 255.0);
    togl->GreenMap[index] = (GLfloat) (entry.peGreen / 255.0);
    togl->BlueMap[index] = (GLfloat) (entry.peBlue / 255.0);
    return index;
}

static void
Win32FreeColor(const Togl *togl, unsigned long index)
{
    TkWinColormap *cmap = (TkWinColormap *) Tk_Colormap(togl->TkWin);
    COLORREF cref;
    UINT   count;
    size_t refCount;
    PALETTEENTRY entry, *entries;
    Tcl_HashEntry *entryPtr;

    if (index >= cmap->size) {
        Tcl_Panic("Tried to free a color that isn't allocated.");
    }
    GetPaletteEntries(cmap->palette, index, 1, &entry);
    cref = PALETTERGB(entry.peRed, entry.peGreen, entry.peBlue);
    entryPtr = Tcl_FindHashEntry(&cmap->refCounts, (const char *) &cref);
    if (!entryPtr) {
        Tcl_Panic("Tried to free a color that isn't allocated.");
    }
    refCount = (size_t)Tcl_GetHashValue(entryPtr) - 1;
    if (refCount == 0) {
        count = cmap->size - index;
        entries = (PALETTEENTRY *) ckalloc(sizeof (PALETTEENTRY) * count);
        GetPaletteEntries(cmap->palette, index + 1, count, entries);
        SetPaletteEntries(cmap->palette, index, count, entries);
        SelectPalette(togl->tglGLHdc, cmap->palette, TRUE);
        RealizePalette(togl->tglGLHdc);
        ckfree((char *) entries);
        cmap->size--;
        Tcl_DeleteHashEntry(entryPtr);
    } else {
        Tcl_SetHashValue(entryPtr, (void *)refCount);
    }
}

static void
Win32SetColor(const Togl *togl,
        unsigned long index, float red, float green, float blue)
{
    TkWinColormap *cmap = (TkWinColormap *) Tk_Colormap(togl->TkWin);
    PALETTEENTRY entry;

    entry.peRed = (unsigned char) (red * 255 + .5);
    entry.peGreen = (unsigned char) (green * 255 + .5);
    entry.peBlue = (unsigned char) (blue * 255 + .5);
    entry.peFlags = 0;
    SetPaletteEntries(cmap->palette, index, 1, &entry);
    SelectPalette(togl->tglGLHdc, cmap->palette, TRUE);
    RealizePalette(togl->tglGLHdc);

    togl->RedMap[index] = (GLfloat) (entry.peRed / 255.0);
    togl->GreenMap[index] = (GLfloat) (entry.peGreen / 255.0);
    togl->BlueMap[index] = (GLfloat) (entry.peBlue / 255.0);
}
#endif /* TOGL_X11 */


unsigned long
Togl_AllocColor(const Togl *togl, float red, float green, float blue)
{
    if (togl->RgbaFlag) {
        (void) fprintf(stderr,
                "Error: Togl_AllocColor illegal in RGBA mode.\n");
        return 0;
    }
    /* TODO: maybe not... */
    if (togl->PrivateCmapFlag) {
        (void) fprintf(stderr,
                "Error: Togl_AllocColor illegal with private colormap\n");
        return 0;
    }
#if defined(TOGL_X11)
    {
        XColor  xcol;
        int     exact;

        xcol.red = (short) (red * 65535.0);
        xcol.green = (short) (green * 65535.0);
        xcol.blue = (short) (blue * 65535.0);

        noFaultXAllocColor(Tk_Display(togl->TkWin), Tk_Colormap(togl->TkWin),
                Tk_Visual(togl->TkWin)->map_entries, &xcol, &exact);
        togl->RedMap[xcol.pixel] = (float) xcol.red / 65535.0;
        togl->GreenMap[xcol.pixel] = (float) xcol.green / 65535.0;
        togl->BlueMap[xcol.pixel] = (float) xcol.blue / 65535.0;

        return xcol.pixel;
    }

#elif defined(TOGL_WGL)
    return Win32AllocColor(togl, red, green, blue);

#elif defined(TOGL_AGL)
    /* still need to implement this on Mac... */
    return 0;

#endif
}



void
Togl_FreeColor(const Togl *togl, unsigned long pixel)
{
    if (togl->RgbaFlag) {
        (void) fprintf(stderr, "Error: Togl_FreeColor illegal in RGBA mode.\n");
        return;
    }
    /* TODO: maybe not... */
    if (togl->PrivateCmapFlag) {
        (void) fprintf(stderr,
                "Error: Togl_FreeColor illegal with private colormap\n");
        return;
    }
#if defined(TOGL_X11)
    (void) XFreeColors(Tk_Display(togl->TkWin), Tk_Colormap(togl->TkWin),
            &pixel, 1, 0);
#elif defined(TOGL_WGL)
    Win32FreeColor(togl, pixel);
#endif
}



void
Togl_SetColor(const Togl *togl,
        unsigned long index, float red, float green, float blue)
{

    if (togl->RgbaFlag) {
        (void) fprintf(stderr, "Error: Togl_SetColor illegal in RGBA mode.\n");
        return;
    }
    if (!togl->PrivateCmapFlag) {
        (void) fprintf(stderr,
                "Error: Togl_SetColor requires a private colormap\n");
        return;
    }
#if defined(TOGL_X11)
    {
        XColor  xcol;

        xcol.pixel = index;
        xcol.red = (short) (red * 65535.0);
        xcol.green = (short) (green * 65535.0);
        xcol.blue = (short) (blue * 65535.0);
        xcol.flags = DoRed | DoGreen | DoBlue;

        (void) XStoreColor(Tk_Display(togl->TkWin), Tk_Colormap(togl->TkWin),
                &xcol);

        togl->RedMap[xcol.pixel] = (float) xcol.red / 65535.0;
        togl->GreenMap[xcol.pixel] = (float) xcol.green / 65535.0;
        togl->BlueMap[xcol.pixel] = (float) xcol.blue / 65535.0;
    }
#elif defined(TOGL_WGL)
    Win32SetColor(togl, index, red, green, blue);
#endif
}


#  define MAX_FONTS 1000
static GLuint ListBase[MAX_FONTS];
static GLuint ListCount[MAX_FONTS];

/*
 * Load the named bitmap font as a sequence of bitmaps in a display list.
 * fontname may be one of the predefined fonts like TOGL_FONT_MONOSPACE
 * or an X font name.
 */
/* TCL3D Begin change Font */
static GLuint
Togl_LoadBitmapFont (const Togl *togl, const char *fontname,
                     const char *fontsize, const char *weight,
                     const char *slant)
{
    static Bool FirstTime = True;

#  if defined(TOGL_X11)
    XFontStruct *fontinfo;
#  elif defined(TOGL_WGL)
    FontAttributes fa;
    XLFDAttributes xa;
    HFONT   font;
    HFONT   oldFont;
    TEXTMETRIC tm;
#  endif
    int     first, last, count;
    GLuint  fontbase;
    const char *name;
    char simpleFont[100];

    /* Initialize the ListBase and ListCount arrays */
    if (FirstTime) {
        int i;

        for (i = 0; i < MAX_FONTS; i++) {
            ListBase[i] = ListCount[i] = 0;
        }
        FirstTime = False;
    }

    #define FONT_TEMPLATE "-*-%s-%s-%c-*-*-%s-*-*-*-*-*-*-*"

    if (!strchr (fontname, '-')) {
        /* Font specification does not seem to be in XLFD notation. */
        sprintf (simpleFont, FONT_TEMPLATE, fontname, weight, slant[0], fontsize);
        name = simpleFont;
    } else {
        name = (const char *) fontname;
    }
    /* TCL3D End change Font */

    sprintf (sDebugBuf, "FontName: %s\n", name); DEBUG_OUT (sDebugBuf);

    assert(name);

#  if defined(TOGL_X11)
    fontinfo = (XFontStruct *) XLoadQueryFont(Tk_Display(togl->TkWin), name);
    if (!fontinfo) {
        return 0;
    }
    first = fontinfo->min_char_or_byte2;
    last = fontinfo->max_char_or_byte2;

    sprintf (sDebugBuf, "Font: (%d, %d)\n", first, last); DEBUG_OUT (sDebugBuf);

#  elif defined(TOGL_WGL)
    if (TCL_OK != FontParseXLFD(name, &fa, &xa)) {
        return 0;
    }

    font = CreateFontA (fa.size,                   /* Height Of Font */
                       0,                         /* Width Of Font */
                       0,                         /* Angle Of Escapement */
                       0,                         /* Orientation Angle */
                       fa.weight,                 /* Font Weight */
                       fa.slant,                  /* Italic */
                       fa.underline,              /* Underline */
                       fa.overstrike,             /* Strikeout */
                       ANSI_CHARSET,              /* Character Set Identifier */
                       OUT_TT_PRECIS,             /* Output Precision */
                       CLIP_DEFAULT_PRECIS,       /* Clipping Precision */
                       ANTIALIASED_QUALITY,       /* Output Quality */
                       FF_DONTCARE|DEFAULT_PITCH, /* Family And Pitch */
                       fa.family);                /* Font Name */
    if (!font) {
        return 0;
    }
    oldFont = SelectObject(togl->tglGLHdc, font);
    GetTextMetrics(togl->tglGLHdc, &tm);
    first =   0;  /* Was tm.tmFirstChar */
    last  = 255;  /* Was tm.tmLastChar  */

    sprintf (sDebugBuf, "Font: %s Size: %d (%d, %d)\n", fa.family, fa.size, first, last); DEBUG_OUT (sDebugBuf);

#  elif defined(TOGL_AGL)
    first =   0;
    last  = 255;
#  endif

    count = last - first + 1;
    fontbase = glGenLists((GLuint) (last + 1));
    if (fontbase == 0) {
#  ifdef TOGL_WGL
        SelectObject(togl->tglGLHdc, oldFont);
        DeleteObject(font);
#  endif
        return 0;
    }
#  if defined(TOGL_WGL)
    wglUseFontBitmaps(togl->tglGLHdc, first, count, (int) fontbase + first);
    SelectObject(togl->tglGLHdc, oldFont);
    DeleteObject(font);
#  elif defined(TOGL_X11)
    glXUseXFont(fontinfo->fid, first, count, (int) fontbase + first);
#  elif defined(TOGL_AGL)
    aglUseFont(togl->Ctx, 1, 0, atoi(fontsize), first, count, fontbase + first);
#  endif

    /* Record the list base and number of display lists for
     * Togl_UnloadBitmapFont(). */
    {
        int i;
        for (i = 0; i < MAX_FONTS; i++) {
            if (ListBase[i] == 0) {
                ListBase[i] = fontbase;
                ListCount[i] = last + 1;
                break;
            }
        }
    }

    return fontbase;
}


/*
 * Release the display lists which were generated by Togl_LoadBitmapFont().
 */
static void
Togl_UnloadBitmapFont(const Togl *togl, GLuint fontbase)
{
    int i;

    (void) togl;
    for (i = 0; i < MAX_FONTS; i++) {
        if (ListBase[i] == fontbase) {
            glDeleteLists(ListBase[i], ListCount[i]);
            ListBase[i] = ListCount[i] = 0;
            return;
        }
    }
}


/*
 * Overlay functions
 */

static void
Togl_UseLayer(Togl *togl, int layer)
{
    if (layer == TOGL_NORMAL) {
#if defined(TOGL_WGL)
        int res = wglMakeCurrent(togl->tglGLHdc, togl->Ctx);
        if (!res) {
            Tcl_Panic("wglMakeCurrent");
        }
#elif defined(TOGL_X11)
        (void) glXMakeCurrent(Tk_Display(togl->TkWin),
                Tk_WindowId(togl->TkWin), togl->Ctx);
#elif defined(TOGL_AGL)
        (void) aglSetCurrentContext(togl->Ctx);
#endif
    } else if (layer == TOGL_OVERLAY && togl->OverlayWindow) {
#if defined(TOGL_WGL)
        int res = wglMakeCurrent(togl->tglGLHdc, togl->tglGLOverlayHglrc);
        if (!res) {
            Tcl_Panic("wglMakeCurrent overlay");
        }
#elif defined(TOGL_X11)
        (void) glXMakeCurrent(Tk_Display(togl->TkWin),
                togl->OverlayWindow, togl->OverlayCtx);
#elif defined(TOGL_AGL)
#endif
    } else {
        /* error */
    }
}


static void
Togl_ShowOverlay(Togl *togl)
{
#if defined(TOGL_X11)           /* not yet implemented on Windows */
    if (togl->OverlayWindow) {
        (void) XMapWindow(Tk_Display(togl->TkWin), togl->OverlayWindow);
        (void) XInstallColormap(Tk_Display(togl->TkWin), togl->OverlayCmap);
        togl->OverlayIsMapped = True;
    }
#endif
}


static void
Togl_HideOverlay(Togl *togl)
{
    if (togl->OverlayWindow && togl->OverlayIsMapped) {
        (void) XUnmapWindow(Tk_Display(togl->TkWin), togl->OverlayWindow);
        togl->OverlayIsMapped = False;
    }
}


static void
Togl_PostOverlayRedisplay(Togl *togl)
{
    if (!togl->OverlayUpdatePending
            && togl->OverlayWindow && togl->OverlayDisplayProc) {
        Tcl_DoWhenIdle(Togl_RenderOverlay, (ClientData) togl);
        togl->OverlayUpdatePending = True;
    }
}


static int
Togl_ExistsOverlay(const Togl *togl)
{
    return togl->OverlayFlag;
}


static int
Togl_GetOverlayTransparentValue(const Togl *togl)
{
    return togl->OverlayTransparentPixel;
}


static int
Togl_IsMappedOverlay(const Togl *togl)
{
    return togl->OverlayFlag && togl->OverlayIsMapped;
}


unsigned long
Togl_AllocColorOverlay(const Togl *togl, float red, float green, float blue)
{
#if defined(TOGL_X11)           /* not yet implemented on Windows */
    if (togl->OverlayFlag && togl->OverlayCmap) {
        XColor  xcol;

        xcol.red = (short) (red * 65535.0);
        xcol.green = (short) (green * 65535.0);
        xcol.blue = (short) (blue * 65535.0);
        if (!XAllocColor(Tk_Display(togl->TkWin), togl->OverlayCmap, &xcol))
            return (unsigned long) -1;
        return xcol.pixel;
    }
#endif /* TOGL_X11 */
    return (unsigned long) -1;
}


void
Togl_FreeColorOverlay(const Togl *togl, unsigned long pixel)
{
#if defined(TOGL_X11)           /* not yet implemented on Windows */
    if (togl->OverlayFlag && togl->OverlayCmap) {
        (void) XFreeColors(Tk_Display(togl->TkWin), togl->OverlayCmap, &pixel,
                1, 0);
    }
#endif /* TOGL_X11 */
}


/*
 * User client data
 */


ClientData
Togl_GetClientData(const Togl *togl)
{
    return togl->Client_Data;
}


void
Togl_SetClientData(Togl *togl, ClientData clientData)
{
    togl->Client_Data = clientData;
}


#ifdef MESA_COLOR_HACK
/*
 * Let's know how many free colors do we have
 */

#  define RLEVELS     5
#  define GLEVELS     9
#  define BLEVELS     5

/* to free dithered_rgb_colormap pixels allocated by Mesa */
static unsigned long *ToglMesaUsedPixelCells = NULL;
static int ToglMesaUsedFreeCells = 0;

static int
get_free_color_cells(Display *display, int screen, Colormap colormap)
{
    if (!ToglMesaUsedPixelCells) {
        XColor  xcol;
        int     i;
        int     colorsfailed, ncolors = XDisplayCells(display, screen);

        long    r, g, b;

        ToglMesaUsedPixelCells =
                (unsigned long *) ckalloc(ncolors * sizeof (unsigned long));

        /* Allocate X colors and initialize color_table[], red_table[], etc */
        /* de Mesa 2.1: xmesa1.c setup_dithered_(...) */
        i = colorsfailed = 0;
        for (r = 0; r < RLEVELS; r++)
            for (g = 0; g < GLEVELS; g++)
                for (b = 0; b < BLEVELS; b++) {
                    int     exact;

                    xcol.red = (r * 65535) / (RLEVELS - 1);
                    xcol.green = (g * 65535) / (GLEVELS - 1);
                    xcol.blue = (b * 65535) / (BLEVELS - 1);
                    noFaultXAllocColor(display, colormap, ncolors,
                            &xcol, &exact);
                    ToglMesaUsedPixelCells[i++] = xcol.pixel;
                    if (!exact) {
                        colorsfailed++;
                    }
                }
        ToglMesaUsedFreeCells = i;

        XFreeColors(display, colormap, ToglMesaUsedPixelCells,
                ToglMesaUsedFreeCells, 0x00000000);
    }
    return ToglMesaUsedFreeCells;
}


static void
free_default_color_cells(Display *display, Colormap colormap)
{
    if (ToglMesaUsedPixelCells) {
        XFreeColors(display, colormap, ToglMesaUsedPixelCells,
                ToglMesaUsedFreeCells, 0x00000000);
        ckfree((char *)ToglMesaUsedPixelCells);
        ToglMesaUsedPixelCells = NULL;
        ToglMesaUsedFreeCells = 0;
    }
}
#endif

int
Togl_GetToglFromObj(Tcl_Interp *interp, Tcl_Obj *obj, Togl **toglPtr)
{
    Tcl_Command toglCmd;
    Tcl_CmdInfo info;

    toglCmd = Tcl_GetCommandFromObj(interp, obj);
    if (Tcl_GetCommandInfoFromToken(toglCmd, &info) == 0
            || info.objProc != Togl_ObjWidget) {
        Tcl_AppendResult(interp, "Expected togl command argument", NULL);
        return TCL_ERROR;
    }
    *toglPtr = (Togl *) info.objClientData;
    return TCL_OK;
}

int
Togl_GetToglFromName(Tcl_Interp *interp, const char *cmdName, Togl **toglPtr)
{
    Tcl_CmdInfo info;

    if (Tcl_GetCommandInfo(interp, cmdName, &info) == 0
            || info.objProc != Togl_ObjWidget) {
        Tcl_AppendResult(interp, "Expected togl command argument", NULL);
        return TCL_ERROR;
    }
    *toglPtr = (Togl *) info.objClientData;
    return TCL_OK;
}

int
BrlcadTkObolHost_SetSwapInterval(Tcl_Interp *interp,
        const char *widget_command, int interval)
{
    Togl *togl = NULL;

    if (!interp || !widget_command ||
            Togl_GetToglFromName(interp, widget_command, &togl) != TCL_OK ||
            !togl || !togl->Ctx)
        return 0;
    Togl_MakeCurrent(togl);
    if (!Togl_SwapInterval(togl, interval ? 1 : 0))
        return 0;
    togl->SwapInterval = interval ? 1 : 0;
    return 1;
}
