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
 * The modifications are based on Togl Version 2.0.
 * Changes and additions are marked with TCL3D and are
 * Copyright 2005-2025 Paul Obermeier
 *
 * This file implements parsing of font names specified in XLFD notation.
 * It is needed for the Windows port only.
 *
 * Code taken from Tk's file tkFont.c.
 * This is the copyright information from tkFont.c
 * Copyright (c) 1990-1994 The Regents of the University of California.
 * Copyright (c) 1994-1998 Sun Microsystems, Inc.
 *
 * See the file "Tcl3D_License.txt" for information on usage and redistribution
 * of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 */

/* This file only gets compiled under Windows */

#if defined(TOGL_WGL)

#include <stdlib.h>
#include <tcl.h>
#include <tk.h>

#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  undef WIN32_LEAN_AND_MEAN

#include "toglFont.h"

/*
 *---------------------------------------------------------------------------
 *
 * FindStateNum --
 *
 *      Given a lookup table, map a string to a number in the table.
 *
 * Results:
 *      If strKey was equal to the string keys of one of the elements
 *      in the table, returns the numeric key of that element.
 *      Returns the numKey associated with the last element (the NULL
 *      string one) in the table if strKey was not equal to any of the
 *      string keys in the table.  In that case, an error message is
 *      also left in the interp's result (if interp is not NULL).
 *
 * Side effects.
 *      None.
 *
 *---------------------------------------------------------------------------
 */

static int
FindStateNum(interp, option, mapPtr, strKey)
    Tcl_Interp *interp;         /* Interp for error reporting. */
    const char *option;         /* String to use when constructing error. */
    const StateMap *mapPtr;     /* Lookup table. */
    const char *strKey;         /* String to try to find in lookup table. */
{
    const StateMap *mPtr;

    for (mPtr = mapPtr; mPtr->strKey != NULL; mPtr++) {
        if (strcmp(strKey, mPtr->strKey) == 0) {
            return mPtr->numKey;
        }
    }
    if (interp != NULL) {
        mPtr = mapPtr;
        Tcl_AppendResult(interp, "bad ", option, " value \"", strKey,
                "\": must be ", mPtr->strKey, (char *) NULL);
        for (mPtr++; mPtr->strKey != NULL; mPtr++) {
            Tcl_AppendResult(interp,
                    ((mPtr[1].strKey != NULL) ? ", " : ", or "),
                    mPtr->strKey, (char *) NULL);
        }
    }
    return mPtr->numKey;
}

/*
 * The following structures are used as two-way maps between the values for
 * the fields in the TkFontAttributes structure and the strings used in
 * Tcl, when parsing both option-value format and style-list format font
 * name strings.
 */

static StateMap weightMap[] = {
    {FW_NORMAL, "normal"},
    {FW_BOLD,   "bold"},
    {FW_DONTCARE,       NULL}
};

static StateMap slantMap[] = {
    {0, "roman"},
    {1, "italic"},
    {0, NULL}
};

static StateMap underlineMap[] = {
    {1,                 "underline"},
    {0,                 NULL}
};

static StateMap overstrikeMap[] = {
    {1,                 "overstrike"},
    {0,                 NULL}
};

/*
 * The following structures are used when parsing XLFD's into a set of
 * TkFontAttributes.
 */

static StateMap xlfdWeightMap[] = {
    {FW_NORMAL, "normal"},
    {FW_NORMAL, "medium"},
    {FW_NORMAL, "book"},
    {FW_NORMAL, "light"},
    {FW_BOLD,   "bold"},
    {FW_BOLD,   "demi"},
    {FW_BOLD,   "demibold"},
    {FW_NORMAL, NULL}           /* Assume anything else is "normal". */
};

#define TOGL_FS_ROMAN   0
#define TOGL_FS_ITALIC  1
#define TOGL_FS_OBLIQUE 0
#define TOGL_FS_ROMAN   0

static StateMap xlfdSlantMap[] = {
    {TOGL_FS_ROMAN,     "r"},
    {TOGL_FS_ITALIC,    "i"},
    {TOGL_FS_OBLIQUE,   "o"},
    {TOGL_FS_ROMAN,     NULL}           /* Assume anything else is "roman". */
};

/*
 * Possible values for the "setwidth" field in a TkXLFDAttributes structure.
 * The setwidth is whether characters are considered wider or narrower than
 * normal.
 */

#define TOGL_SW_NORMAL   0
#define TOGL_SW_CONDENSE 1
#define TOGL_SW_EXPAND   2
#define TOGL_SW_UNKNOWN  3      /* Unknown setwidth.  This value may be
                                 * stored in the setwidth field. */

static StateMap xlfdSetwidthMap[] = {
    {TOGL_SW_NORMAL,    "normal"},
    {TOGL_SW_CONDENSE,  "narrow"},
    {TOGL_SW_CONDENSE,  "semicondensed"},
    {TOGL_SW_CONDENSE,  "condensed"},
    {TOGL_SW_UNKNOWN,   NULL}
};

/*
 *---------------------------------------------------------------------------
 *
 * FieldSpecified --
 *
 *      Helper function for TkParseXLFD().  Determines if a field in the
 *      XLFD was set to a non-null, non-don't-care value.
 *
 * Results:
 *      The return value is 0 if the field in the XLFD was not set and
 *      should be ignored, non-zero otherwise.
 *
 * Side effects:
 *      None.
 *
 *---------------------------------------------------------------------------
 */

static int
FieldSpecified(field)
    const char *field;  /* The field of the XLFD to check.  Strictly
                         * speaking, only when the string is "*" does it mean
                         * don't-care.  However, an unspecified or question
                         * mark is also interpreted as don't-care. */
{
    char ch;

    if (field == NULL) {
        return 0;
    }
    ch = field[0];
    return (ch != '*' && ch != '?');
}

/*
 *---------------------------------------------------------------------------
 *
 * TkFontParseXLFD --
 *
 *      Break up a fully specified XLFD into a set of font attributes.
 *
 * Results:
 *      Return value is TCL_ERROR if string was not a fully specified XLFD.
 *      Otherwise, fills font attribute buffer with the values parsed
 *      from the XLFD and returns TCL_OK.
 *
 * Side effects:
 *      None.
 *
 *---------------------------------------------------------------------------
 */

#define InitFontAttributes(fa)   memset((fa), 0, sizeof(FontAttributes));
#define InitXLFDAttributes(xa)   memset((xa), 0, sizeof(XLFDAttributes));
#define UCHAR(c) ((unsigned char) (c))

int
FontParseXLFD(string, faPtr, xaPtr)
    const char *string;         /* Parseable font description string. */
    FontAttributes *faPtr;      /* Filled with attributes parsed from font
                                 * name.  Any attributes that were not
                                 * specified in font name are filled with
                                 * default values. */
    XLFDAttributes *xaPtr;      /* Filled with X-specific attributes parsed
                                 * from font name.  Any attributes that were
                                 * not specified in font name are filled with
                                 * default values.  May be NULL if such
                                 * information is not desired. */
{
    char *src;
    const char *str;
    int i, j;
    char *field[XLFD_NUMFIELDS + 2];
    Tcl_DString ds;
    XLFDAttributes xa;
    char msgStr[512];

    if (xaPtr == NULL) {
        xaPtr = &xa;
    }
    InitFontAttributes(faPtr);
    InitXLFDAttributes(xaPtr);

    memset(field, '\0', sizeof(field));

    str = string;
    if (*str == '-') {
        str++;
    }

    Tcl_DStringInit(&ds);
    Tcl_DStringAppend(&ds, (char *) str, -1);
    src = Tcl_DStringValue(&ds);

    field[0] = src;
    for (i = 0; *src != '\0'; src++) {
        if (!(*src & 0x80)
                && Tcl_UniCharIsUpper(UCHAR(*src))) {
            *src = (char) Tcl_UniCharToLower(UCHAR(*src));
        }
        if (*src == '-') {
            i++;
            if (i == XLFD_NUMFIELDS) {
                continue;
            }
            *src = '\0';
            field[i] = src + 1;
            if (i > XLFD_NUMFIELDS) {
                break;
            }
        }
    }

    /*
     * An XLFD of the form -adobe-times-medium-r-*-12-*-* is pretty common,
     * but it is (strictly) malformed, because the first * is eliding both
     * the Setwidth and the Addstyle fields.  If the Addstyle field is a
     * number, then assume the above incorrect form was used and shift all
     * the rest of the fields right by one, so the number gets interpreted
     * as a pixelsize.  This fix is so that we don't get a million reports
     * that "it works under X (as a native font name), but gives a syntax
     * error under Windows (as a parsed set of attributes)".
     */

    if ((i > XLFD_ADD_STYLE) && (FieldSpecified(field[XLFD_ADD_STYLE]))) {
        if (atoi(field[XLFD_ADD_STYLE]) != 0) {
            for (j = XLFD_NUMFIELDS - 1; j >= XLFD_ADD_STYLE; j--) {
                field[j + 1] = field[j];
            }
            field[XLFD_ADD_STYLE] = NULL;
            i++;
        }
    }

    /*
     * Bail if we don't have enough of the fields (up to pointsize).
     */

    if (i < XLFD_FAMILY) {
        Tcl_DStringFree(&ds);
        return TCL_ERROR;
    }

    if (FieldSpecified(field[XLFD_FOUNDRY])) {
        strncpy (xaPtr->foundry, field[XLFD_FOUNDRY], MAX_FONT_NAME);
    }

    if (FieldSpecified(field[XLFD_FAMILY])) {
        strncpy (faPtr->family, field[XLFD_FAMILY], MAX_FONT_NAME);
    }
    if (FieldSpecified(field[XLFD_WEIGHT])) {
        faPtr->weight = FindStateNum(NULL, NULL, xlfdWeightMap,
                field[XLFD_WEIGHT]);
    }
    if (FieldSpecified(field[XLFD_SLANT])) {
        xaPtr->slant = FindStateNum(NULL, NULL, xlfdSlantMap,
                field[XLFD_SLANT]);
        if (xaPtr->slant == TOGL_FS_ROMAN) {
            faPtr->slant = TOGL_FS_ROMAN;
        } else {
            faPtr->slant = TOGL_FS_ITALIC;
        }
    }
    if (FieldSpecified(field[XLFD_SETWIDTH])) {
        xaPtr->setwidth = FindStateNum(NULL, NULL, xlfdSetwidthMap,
                field[XLFD_SETWIDTH]);
    }

    /* XLFD_ADD_STYLE ignored. */

    /*
     * Pointsize in tenths of a point, but treat it as tenths of a pixel
     * for historical compatibility.
     */

    faPtr->size = 12;

    if (FieldSpecified(field[XLFD_POINT_SIZE])) {
        if (field[XLFD_POINT_SIZE][0] == '[') {
            /*
             * Some X fonts have the point size specified as follows:
             *
             *      [ N1 N2 N3 N4 ]
             *
             * where N1 is the point size (in points, not decipoints!), and
             * N2, N3, and N4 are some additional numbers that I don't know
             * the purpose of, so I ignore them.
             */

            faPtr->size = atoi(field[XLFD_POINT_SIZE] + 1);
        } else if (Tcl_GetInt(NULL, field[XLFD_POINT_SIZE],
                &faPtr->size) == TCL_OK) {
            faPtr->size /= 10;
        } else {
            return TCL_ERROR;
        }
    }

    /*
     * Pixel height of font.  If specified, overrides pointsize.
     */

    if (FieldSpecified(field[XLFD_PIXEL_SIZE])) {
        if (field[XLFD_PIXEL_SIZE][0] == '[') {
            /*
             * Some X fonts have the pixel size specified as follows:
             *
             *      [ N1 N2 N3 N4 ]
             *
             * where N1 is the pixel size, and where N2, N3, and N4
             * are some additional numbers that I don't know
             * the purpose of, so I ignore them.
             */

            faPtr->size = atoi(field[XLFD_PIXEL_SIZE] + 1);
        } else if (Tcl_GetInt(NULL, field[XLFD_PIXEL_SIZE],
                &faPtr->size) != TCL_OK) {
            return TCL_ERROR;
        }
    }

    faPtr->size = -faPtr->size;

    /* XLFD_RESOLUTION_X ignored. */

    /* XLFD_RESOLUTION_Y ignored. */

    /* XLFD_SPACING ignored. */

    /* XLFD_AVERAGE_WIDTH ignored. */

    if (FieldSpecified(field[XLFD_CHARSET])) {
        strncpy (xaPtr->charset, field[XLFD_CHARSET], MAX_FONT_NAME);
    } else {
        strncpy (xaPtr->charset, "iso8859-1", MAX_FONT_NAME);
    }
    Tcl_DStringFree(&ds);
    return TCL_OK;
}

#endif /* TOGL_WGL */
