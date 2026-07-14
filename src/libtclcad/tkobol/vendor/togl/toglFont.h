/* 
 * Togl - a Tk OpenGL widget
 *
 * Copyright 1996-2002  Brian Paul and Ben Bederson
 * See the LICENSE file for copyright details.
 *
 * This is an enhanced version of Togl called tcl3dTogl. 
 * It adds functionality, so that Tcl-wrapped OpenGL functions
 * can be evaluated by a Togl widget. 
 * It is part of the Tcl3D package.
 * The modifications are based on Togl Version 1.7.
 * Changes and additions are marked with TCL3D and are
 * Copyright 2005-2025 Paul Obermeier
 *
 * This file implements parsing of font names specified in XLFD notation.
 * It is needed for the Windows port only.
 *
 * Code taken from Tk's file tkFont.h.
 *
 * This is the copyright information from tkFont.h
 * Copyright (c) 1996-1997 Sun Microsystems, Inc.
 *
 * See the file "Tcl3D_License.txt" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

#ifndef TOGLFONT_H
#  define TOGLFONT_H

/*
 * The following defines specify the meaning of the fields in a fully
 * qualified XLFD.
 */

#define XLFD_FOUNDRY        0
#define XLFD_FAMILY         1
#define XLFD_WEIGHT         2
#define XLFD_SLANT          3
#define XLFD_SETWIDTH       4
#define XLFD_ADD_STYLE      5
#define XLFD_PIXEL_SIZE     6
#define XLFD_POINT_SIZE     7
#define XLFD_RESOLUTION_X   8
#define XLFD_RESOLUTION_Y   9
#define XLFD_SPACING        10
#define XLFD_AVERAGE_WIDTH  11
#define XLFD_CHARSET        12
#define XLFD_NUMFIELDS      13  /* Number of fields in XLFD. */

#define MAX_FONT_NAME 100       /* The max. length of a font specification */

/*
 * The following structure is used to return attributes when parsing an
 * XLFD.  The extra information is of interest to the Unix-specific code
 * when attempting to find the closest matching font.
 */

typedef struct XLFDAttributes {
    char foundry[MAX_FONT_NAME]; /* The foundry of the font. */
    int slant;                   /* The tristate value for the slant, which
                                  * is significant under X. */
    int setwidth;                /* The proportionate width, see below for
                                  * definition. */
    char charset[MAX_FONT_NAME]; /* The actual charset string. */
} XLFDAttributes;

/*
 * The following structure keeps track of the attributes of a font.  It can
 * be used to keep track of either the desired attributes or the actual
 * attributes gotten when the font was instantiated.
 */

typedef struct FontAttributes {
    char family[MAX_FONT_NAME]; /* Font family, or NULL to represent
                                 * plaform-specific default system font. */
    int size;                   /* Pointsize of font, 0 for default size, or
                                 * negative number meaning pixel size. */
    int weight;                 /* Weight flag; see below for def'n. */
    int slant;                  /* Slant flag; see below for def'n. */
    int underline;              /* Non-zero for underline font. */
    int overstrike;             /* Non-zero for overstrike font. */
} FontAttributes;

/*
 * The following structure is used as a two way map between integers
 * and strings, usually to map between an internal C representation
 * and the strings used in Tcl.
 */

typedef struct StateMap {
    int numKey;                 /* Integer representation of a value. */
    char *strKey;               /* String representation of a value. */
} StateMap;

extern int
FontParseXLFD (const char *string, FontAttributes *faPtr, XLFDAttributes *xaPtr);

#endif /* TOGLFONT_H */
