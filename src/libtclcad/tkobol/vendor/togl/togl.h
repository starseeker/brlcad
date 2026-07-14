/* 
 * Togl - a Tk OpenGL widget
 *
 * Copyright 1996-1998  Brian Paul and Ben Bederson
 * Copyright 2005-2008  Greg Couch
 * See the LICENSE file for copyright details.
 *
 * This is an enhanced version of Togl called tcl3dTogl. 
 * It adds functionality, so that Tcl-wrapped OpenGL functions
 * can be evaluated by a Togl widget. 
 * It is part of the Tcl3D package.
 * The modifications are based on Togl Version 2.0.
 * Changes and additions are marked with TCL3D and are
 * Copyright 2005-2025 Paul Obermeier
 */


#ifndef TOGL_H
#  define TOGL_H

#  ifdef TOGL_WGL
#    define WIN32_LEAN_AND_MEAN
#    include <windows.h>
#    undef WIN32_LEAN_AND_MEAN
#    if defined(_MSC_VER)
#       define DllEntryPoint DllMain
#    endif
#  endif

#  define TOGL_EXTERN extern

#  ifdef TOGL_AGL
#    ifndef MAC_OSX_TCL
#      define MAC_OSX_TCL 1
#    endif
#    ifndef MAC_OSX_TK
#      define MAC_OSX_TK 1
#    endif
#  endif

#  include <tcl.h>
#  include <tk.h>
#  include "tclcad/defines.h"

#  ifdef TOGL_AGL
#    include <OpenGL/gl.h>
#    include <OpenGL/glu.h>
#  else
#    include <GL/gl.h>
#    include <GL/glu.h>
#  endif

#  ifdef __sgi
#    include <GL/glx.h>
#  endif

#  ifndef NULL
#    define NULL 0
#  endif

#  ifdef __cplusplus
   extern "C" {
#  endif

#  define TOGL_VERSION "1.0.2"
#  define TOGL_MAJOR_VERSION 1
#  define TOGL_MINOR_VERSION 0

/* 
 * Normal and overlay plane constants
 */
#  define TOGL_NORMAL   1
#  define TOGL_OVERLAY  2

struct Togl;
typedef struct Togl Togl;

typedef void (Togl_Callback) (Togl *togl);
typedef int (Togl_CmdProc) (Togl *togl, int argc, const char *argv[]);

TCLCAD_EXPORT int BrlcadTkObolHost_Init(Tcl_Interp *interp);
TCLCAD_EXPORT int BrlcadTkObolHost_SetSwapInterval(
    Tcl_Interp *interp, const char *widget_command, int interval);

#  ifdef __cplusplus
   }
#  endif

#endif
