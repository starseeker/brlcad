/*                  D M - T K S W R A S T . C P P
 * BRL-CAD
 */

#include "common.h"

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "tcl.h"
#include "tk.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "OSMesa/gl.h"
#include "OSMesa/osmesa.h"

extern "C" {
#include "dm.h"
#include "../dm-gl.h"
#include "../include/private.h"
#include "../swrast/dm-swrast.h"
}

struct tkswrast_vars {
    Tk_Window top;
    Tk_Window xtkwin;
    Display *dpy;
    Window win;
    int (*orig_configureWin)(struct dm *, int);
    int (*orig_close)(struct dm *);
    unsigned char *x_b;
    size_t x_b_size;
};

static int tkswrast_debug(void);
static void tkswrast_log(const char *fmt, ...);

static int
tkswrast_debug(void)
{
    const char *d = getenv("TKSWRAST_DEBUG");
    return (d && BU_STR_EQUAL(d, "1"));
}

static void
tkswrast_log(const char *fmt, ...)
{
    static FILE *lfp = NULL;
    if (!tkswrast_debug())
        return;

    if (!lfp) {
        const char *lfile = getenv("TKSWRAST_LOG_FILE");
        if (!lfile || !strlen(lfile))
            lfile = "swrast_log.txt";
        lfp = fopen(lfile, "a");
        if (!lfp)
            return;
    }

    va_list ap;
    va_start(ap, fmt);
    vfprintf(lfp, fmt, ap);
    va_end(ap);
    fflush(lfp);
}

static void
tkswrast_log_tk_hierarchy(struct dm *dmp)
{
    if (!tkswrast_debug() || !dmp || !dmp->i || !dmp->i->dm_interp)
        return;

    Tcl_Interp *ti = (Tcl_Interp *)dmp->i->dm_interp;
    struct bu_vls cmd = BU_VLS_INIT_ZERO;
    struct bu_vls out = BU_VLS_INIT_ZERO;

    bu_vls_printf(&cmd,
        "set _w %s; set _p [winfo parent $_w]; set _mgr [winfo manager $_w]; "
        "set _g [winfo geometry $_w]; set _pg [winfo geometry $_p]; "
        "set _ch [winfo children $_p]; list $_w $_p $_mgr $_g $_pg $_ch",
        bu_vls_addr(&dmp->i->dm_pathName));
    if (Tcl_Eval(ti, bu_vls_addr(&cmd)) == TCL_OK) {
        bu_vls_printf(&out, "%s", Tcl_GetStringResult(ti));
        tkswrast_log("tkswrast_tk[%s]: %s\n", bu_vls_addr(&dmp->i->dm_pathName), bu_vls_addr(&out));
    } else {
        tkswrast_log("tkswrast_tk[%s]: hierarchy query failed: %s\n",
            bu_vls_addr(&dmp->i->dm_pathName), Tcl_GetStringResult(ti));
    }

    bu_vls_free(&out);
    bu_vls_free(&cmd);
}

static int
tkswrast_configureWin(struct dm *dmp, int force)
{
    struct tkswrast_vars *tv = (struct tkswrast_vars *)dmp->i->dm_udata;
    struct swrast_vars *sv = (struct swrast_vars *)dmp->i->dm_vars.priv_vars;
    if (!tv || !sv || !sv->v)
        return BRLCAD_ERROR;

    int width = Tk_Width(tv->xtkwin);
    int height = Tk_Height(tv->xtkwin);
    if (tkswrast_debug()) {
        tkswrast_log("tkswrast_configureWin[%s]: force=%d tk=%dx%d dm(before)=%dx%d top=%d\n",
            bu_vls_addr(&dmp->i->dm_pathName), force, width, height,
            dmp->i->dm_width, dmp->i->dm_height, dmp->i->dm_top);
        tkswrast_log_tk_hierarchy(dmp);
    }
    if (width < 1 || height < 1)
        return BRLCAD_OK;

    if (!force && dmp->i->dm_width == width && dmp->i->dm_height == height)
        return BRLCAD_OK;

    sv->v->gv_width = width;
    sv->v->gv_height = height;

    dmp->i->dm_width = width;
    dmp->i->dm_height = height;
    dmp->i->dm_aspect = (fastf_t)width / (fastf_t)height;

    if (tkswrast_debug()) {
        tkswrast_log("tkswrast_configureWin[%s]: dm(after)=%dx%d aspect=%g\n",
            bu_vls_addr(&dmp->i->dm_pathName), dmp->i->dm_width, dmp->i->dm_height, dmp->i->dm_aspect);
    }

    return tv->orig_configureWin(dmp, 1);
}

static int
tkswrast_SwapBuffers(struct dm *dmp)
{
    struct tkswrast_vars *tv = (struct tkswrast_vars *)dmp->i->dm_udata;
    struct swrast_vars *sv = (struct swrast_vars *)dmp->i->dm_vars.priv_vars;
    if (!tv || !sv || !sv->ctx || !tv->dpy || !tv->win)
        return BRLCAD_OK;

    int ww = Tk_Width(tv->xtkwin);
    int wh = Tk_Height(tv->xtkwin);
    if (ww > 0 && wh > 0 && (ww != dmp->i->dm_width || wh != dmp->i->dm_height)) {
        (void)dmp->i->dm_configureWin(dmp, 1);
    }

    int src_w = 0, src_h = 0, format = 0;
    void *buffer = NULL;
    if (!OSMesaGetColorBuffer(sv->ctx, &src_w, &src_h, &format, &buffer) || !buffer)
        return BRLCAD_ERROR;

    int width = src_w;
    int height = src_h;
    int dw = Tk_Width(tv->xtkwin);
    int dh = Tk_Height(tv->xtkwin);
    if (dw > 0 && dh > 0) {
        if (dw < width)
            width = dw;
        if (dh < height)
            height = dh;
    }

    if (tkswrast_debug()) {
        tkswrast_log("tkswrast_swap[%s]: tk=%dx%d osmesa=%dx%d blit=%dx%d dm=%dx%d fmt=%d\n",
            bu_vls_addr(&dmp->i->dm_pathName), ww, wh, src_w, src_h, width, height,
            dmp->i->dm_width, dmp->i->dm_height, format);
    }

    size_t need = (size_t)width * (size_t)height * 4;
    if (!tv->x_b || tv->x_b_size < need) {
        tv->x_b = (unsigned char *)bu_realloc(tv->x_b, need, "tkswrast ximage buffer");
        tv->x_b_size = need;
    }

    unsigned char *src = (unsigned char *)buffer;
    for (int y = 0; y < height; y++) {
        int sy = src_h - 1 - y;
        for (int x = 0; x < width; x++) {
            size_t si = ((size_t)sy * (size_t)src_w + (size_t)x) * 4;
            size_t di = ((size_t)y * (size_t)width + (size_t)x) * 4;
            tv->x_b[di + 0] = src[si + 2];
            tv->x_b[di + 1] = src[si + 1];
            tv->x_b[di + 2] = src[si + 0];
            /* Keep destination alpha opaque for visuals/compositors that
             * honor a 32-bit channel, to avoid dark edge artifacts. */
            tv->x_b[di + 3] = 255;
        }
    }

    XImage *img = XCreateImage(tv->dpy,
            DefaultVisual(tv->dpy, DefaultScreen(tv->dpy)),
            DefaultDepth(tv->dpy, DefaultScreen(tv->dpy)),
            ZPixmap,
            0,
            (char *)tv->x_b,
            width,
            height,
            32,
            0);
    if (!img)
        return BRLCAD_ERROR;

    GC gc = DefaultGC(tv->dpy, DefaultScreen(tv->dpy));
    XPutImage(tv->dpy, tv->win, gc, img, 0, 0, 0, 0, width, height);
    XFlush(tv->dpy);

    img->data = NULL;
    XDestroyImage(img);
    return BRLCAD_OK;
}

static int
tkswrast_doevent(struct dm *dmp, void *UNUSED(vclientData), void *veventPtr)
{
    XEvent *eventPtr = (XEvent *)veventPtr;
    if (!eventPtr)
        return TCL_OK;

    if (eventPtr->type == Expose && eventPtr->xexpose.count == 0) {
        dm_set_dirty(dmp, 1);
        return TCL_RETURN;
    }
    return TCL_OK;
}

static int
tkswrast_close(struct dm *dmp)
{
    struct tkswrast_vars *tv = (struct tkswrast_vars *)dmp->i->dm_udata;
    int (*orig_close)(struct dm *) = NULL;
    if (tv) {
        orig_close = tv->orig_close;
        if (tv->xtkwin)
            Tk_DestroyWindow(tv->xtkwin);
        if (tv->x_b)
            bu_free(tv->x_b, "tkswrast ximage buffer");
        bu_free(tv, "tkswrast vars");
        dmp->i->dm_udata = NULL;
    }
    if (orig_close)
        return orig_close(dmp);
    return BRLCAD_OK;
}

static int
tkswrast_viable(const char *UNUSED(dpy_string))
{
    return 1;
}

static struct dm *
tkswrast_open(void *ctx, void *vinterp, int argc, const char **argv)
{
    Tcl_Interp *interp = (Tcl_Interp *)vinterp;
    Tk_Window tkwin = Tk_MainWindow(interp);
    if (!tkwin)
        return DM_NULL;

    struct dm *dmp = dm_open(ctx, interp, "swrast", argc, argv);
    if (!dmp)
        return DM_NULL;

    dmp->i->dm_interp = (void *)interp;

    struct tkswrast_vars *tv = NULL;
    BU_ALLOC(tv, struct tkswrast_vars);
    tv->orig_configureWin = dmp->i->dm_configureWin;
    tv->orig_close = dmp->i->dm_close;

    char *cp = strrchr(bu_vls_addr(&dmp->i->dm_pathName), '.');
    if (dmp->i->dm_top) {
        tv->xtkwin = Tk_CreateWindowFromPath(interp, tkwin,
            bu_vls_addr(&dmp->i->dm_pathName), bu_vls_addr(&dmp->i->dm_dName));
        tv->top = tv->xtkwin;
    } else {
        if (!cp || cp == bu_vls_addr(&dmp->i->dm_pathName)) {
            tv->top = tkwin;
            cp = (char *)bu_vls_addr(&dmp->i->dm_pathName);
            if (*cp == '.')
                cp++;
        } else {
            struct bu_vls top_vls = BU_VLS_INIT_ZERO;
            bu_vls_strncpy(&top_vls, bu_vls_addr(&dmp->i->dm_pathName), cp - bu_vls_addr(&dmp->i->dm_pathName));
            tv->top = Tk_NameToWindow(interp, bu_vls_addr(&top_vls), tkwin);
            bu_vls_free(&top_vls);
            cp++;
        }
        tv->xtkwin = Tk_CreateWindow(interp, tv->top, cp, (char *)NULL);
    }

    if (tkswrast_debug()) {
        tkswrast_log("tkswrast_open: path=%s top=%d\n", bu_vls_addr(&dmp->i->dm_pathName), dmp->i->dm_top);
    }

    if (!tv->xtkwin) {
        bu_free(tv, "tkswrast vars");
        dm_close(dmp);
        return DM_NULL;
    }

    if (dmp->i->dm_top) {
        Tk_GeometryRequest(tv->xtkwin, dmp->i->dm_width, dmp->i->dm_height);
    } else {
        int pw = Tk_Width(tv->top);
        int ph = Tk_Height(tv->top);
        if (tkswrast_debug()) {
            tkswrast_log("tkswrast_open: embedded parent initial size=%dx%d\n", pw, ph);
        }
        if (pw > 1 && ph > 1) {
            struct swrast_vars *sv = (struct swrast_vars *)dmp->i->dm_vars.priv_vars;
            dmp->i->dm_width = pw;
            dmp->i->dm_height = ph;
            if (sv && sv->v) {
                sv->v->gv_width = pw;
                sv->v->gv_height = ph;
            }
        }
    }
    Tk_SetWindowBackground(tv->xtkwin, BlackPixelOfScreen(Tk_Screen(tv->xtkwin)));
    Tk_MakeWindowExist(tv->xtkwin);
    tv->dpy = Tk_Display(tv->top);
    tv->win = Tk_WindowId(tv->xtkwin);

    dmp->i->dm_udata = (void *)tv;
    dmp->i->dm_id = tv->win;
    bu_vls_sprintf(&dmp->i->dm_tkName, "%s", Tk_Name(tv->xtkwin));

    dmp->i->dm_close = tkswrast_close;
    dmp->i->dm_viable = tkswrast_viable;
    dmp->i->dm_configureWin = tkswrast_configureWin;
    dmp->i->dm_SwapBuffers = tkswrast_SwapBuffers;
    dmp->i->dm_doevent = tkswrast_doevent;
    dmp->i->dm_graphical = 1;
    dmp->i->graphics_system = "Tk";
    dmp->i->dm_name = "tkswrast";
    dmp->i->dm_lname = "Tk OSMesa swrast graphics";

    dm_set_zbuffer(dmp, 1);
    {
        struct gl_vars *mvars = (struct gl_vars *)dmp->i->m_vars;
        if (mvars) {
            mvars->adaptive_zclip = 1;
            mvars->adaptive_zclip_factor = 2.0;
            mvars->adaptive_zclip_min = 100.0;
            mvars->adaptive_zclip_max = 5000.0;
        }
    }

    Tk_MapWindow(tv->xtkwin);
    if (tkswrast_debug()) {
        tkswrast_log("tkswrast_open: child mapped size=%dx%d parent size=%dx%d\n",
            Tk_Width(tv->xtkwin), Tk_Height(tv->xtkwin),
            Tk_Width(tv->top), Tk_Height(tv->top));
        tkswrast_log_tk_hierarchy(dmp);
    }
    (void)dmp->i->dm_configureWin(dmp, 1);

    /* Ensure each embedded Tk DM gets configure callbacks, matching the
     * expected Tk-driven resize behavior used by other windowed DMs. */
    {
        Tcl_Interp *ti = (Tcl_Interp *)dmp->i->dm_interp;
        if (ti) {
            struct bu_vls tcmd = BU_VLS_INIT_ZERO;
            bu_vls_printf(&tcmd, "bind %s <Configure> {catch {%%W configure}}", bu_vls_addr(&dmp->i->dm_pathName));
            (void)Tcl_Eval(ti, bu_vls_addr(&tcmd));
            bu_vls_trunc(&tcmd, 0);
            bu_vls_printf(&tcmd, "event generate %s <Configure>", bu_vls_addr(&dmp->i->dm_pathName));
            (void)Tcl_Eval(ti, bu_vls_addr(&tcmd));
            bu_vls_free(&tcmd);
        }
    }

    dm_set_dirty(dmp, 1);

    return dmp;
}

static struct dm_impl dm_tkswrast_impl;

struct dm dm_tkswrast = { DM_MAGIC, &dm_tkswrast_impl, 0 };

#ifdef DM_PLUGIN
static const struct dm_plugin pinfo = { DM_API, &dm_tkswrast };
extern "C" {
COMPILER_DLLEXPORT const struct dm_plugin *dm_plugin_info(void)
{
    dm_tkswrast_impl.dm_open = tkswrast_open;
    dm_tkswrast_impl.dm_close = tkswrast_close;
    dm_tkswrast_impl.dm_viable = tkswrast_viable;
    dm_tkswrast_impl.dm_configureWin = tkswrast_configureWin;
    dm_tkswrast_impl.dm_SwapBuffers = tkswrast_SwapBuffers;
    dm_tkswrast_impl.dm_doevent = tkswrast_doevent;
    dm_tkswrast_impl.dm_graphical = 1;
    dm_tkswrast_impl.graphics_system = "Tk";
    dm_tkswrast_impl.dm_name = "tkswrast";
    dm_tkswrast_impl.dm_lname = "Tk OSMesa swrast graphics";
    return &pinfo;
}
}
#endif
