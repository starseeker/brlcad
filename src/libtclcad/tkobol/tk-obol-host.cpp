/*                  T K - O B O L - H O S T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file libtclcad/tkobol/tk-obol-host.cpp
 *
 * DM-independent Tk display hosts for Obol endpoints.
 */

#include "common.h"

#include "tcl.h"
#include "tk.h"

#ifdef TKOBOL_X11
#  include <X11/Xlib.h>
#  include <GL/gl.h>
#  include <GL/glx.h>
#elif defined(TKOBOL_WGL)
#  include <windows.h>
#  include <GL/gl.h>
#endif

#include "brlobol/host_factory.h"
#include "brlobol/init.h"
#include "brlobol/view_controller.h"
#include "bu/log.h"
#include "bu/malloc.h"
#include "tclcad/setup.h"

#include <Inventor/SoDB.h>
#include <Inventor/SoRenderManager.h>

#include "vendor/togl/togl.h"

#include <algorithm>
#include <atomic>
#include <cstring>
#include <initializer_list>
#include <new>
#include <string>
#include <vector>

static thread_local unsigned int tk_endpoint_offscreen_depth = 0;

class TkEndpointContextManager : public SoDB::ContextManager {
public:
    TkEndpointContextManager(void) : offscreen(SoDB::createOSMesaContextManager()) {}
    ~TkEndpointContextManager(void) override { delete this->offscreen; }

    void *createOffscreenContext(unsigned int width, unsigned int height) override
    {
	return this->offscreen ?
	    this->offscreen->createOffscreenContext(width, height) : NULL;
    }

    SbBool makeContextCurrent(void *context) override
    {
	if (!this->offscreen || !this->offscreen->makeContextCurrent(context))
	    return FALSE;
	tk_endpoint_offscreen_depth++;
	return TRUE;
    }

    void restorePreviousContext(void *context) override
    {
	if (!this->offscreen)
	    return;
	this->offscreen->restorePreviousContext(context);
	if (tk_endpoint_offscreen_depth)
	    tk_endpoint_offscreen_depth--;
    }

    void destroyContext(void *context) override
    {
	if (this->offscreen)
	    this->offscreen->destroyContext(context);
    }

    SbBool isOSMesaContext(void *context) override
    {
	return this->offscreen ? this->offscreen->isOSMesaContext(context) : FALSE;
    }

    void maxOffscreenDimensions(unsigned int &width,
	unsigned int &height) const override
    {
	if (this->offscreen)
	    this->offscreen->maxOffscreenDimensions(width, height);
	else
	    width = height = 0;
    }

    void getActualSurfaceSize(void *context, unsigned int &width,
	unsigned int &height) const override
    {
	if (this->offscreen)
	    this->offscreen->getActualSurfaceSize(context, width, height);
	else
	    width = height = 0;
    }

    void *getProcAddress(const char *name) override
    {
	if (tk_endpoint_offscreen_depth)
	    return this->offscreen ? this->offscreen->getProcAddress(name) : NULL;
#ifdef TKOBOL_X11
	if (name && glXGetCurrentContext())
	    return reinterpret_cast<void *>(glXGetProcAddressARB(
		reinterpret_cast<const GLubyte *>(name)));
#elif defined(TKOBOL_WGL)
	if (name && wglGetCurrentContext()) {
	    PROC proc = wglGetProcAddress(name);
	    if (proc)
		return reinterpret_cast<void *>(proc);
	}
#endif
	return this->offscreen ? this->offscreen->getProcAddress(name) : NULL;
    }

private:
    SoDB::ContextManager *offscreen;
};

static SoDB::ContextManager *
tk_endpoint_context_manager(void)
{
    static TkEndpointContextManager manager;
    return &manager;
}

class TkObolEndpointHost;

struct TkEndpointEvent {
    Tcl_Event header;
    TkObolEndpointHost *host;
};

class TkObolEndpointHost {
public:
    TkObolEndpointHost(Tcl_Interp *new_interp, bool use_software) :
	interp(new_interp), owner_thread(Tcl_GetCurrentThread()), software(use_software),
	controller(NULL), tkwin(NULL), container(NULL), photo(NULL), opened(false),
	closing(false), event_queued(false), width(1), height(1)
    {
    }

    ~TkObolEndpointHost(void)
    {
	this->close();
    }

    int bind(BRLObolViewController *new_controller)
    {
	this->controller = new_controller;
	return 1;
    }

    int open(const struct brlobol_host_desc *desc)
    {
	if (this->opened || !this->interp || !Tk_MainWindow(this->interp) ||
	    !this->controller || !desc)
	    return 0;
	if (BrlcadTkObolHost_Init(this->interp) != TCL_OK)
	    return 0;

	brlobol_init(tk_endpoint_context_manager());
	this->width = desc->width ? desc->width : 1;
	this->height = desc->height ? desc->height : 1;
	this->controller->setViewportSize(this->width, this->height);

	static std::atomic<unsigned long> next_id(0);
	const unsigned long id = ++next_id;
	this->container_path = desc->native_id_hint && desc->native_id_hint[0] ?
	    desc->native_id_hint : ".brlobol_tk" + std::to_string(id);
	this->container_command_name = "::brlcad::tkobol::endpoint_container_" +
	    std::to_string(id);
	const bool toplevel = desc->mode == BRLOBOL_HOST_MODE_TOPLEVEL;
	if (!toplevel && (!desc->native_id_hint || !desc->native_id_hint[0]))
	    return 0;
	this->widget_path = toplevel ? this->container_path + ".__obol" :
	    this->container_path;
	this->display_command = "::brlcad::tkobol::endpoint_display_" +
	    std::to_string(id);
	this->reshape_command = "::brlcad::tkobol::endpoint_reshape_" +
	    std::to_string(id);
	this->widget_command_name = "::brlcad::tkobol::endpoint_widget_" +
	    std::to_string(id);
	this->photo_name = "::brlcad::tkobol::endpoint_photo_" +
	    std::to_string(id);

	Tcl_CreateObjCommand(this->interp, this->display_command.c_str(),
	    display_callback, this, NULL);
	Tcl_CreateObjCommand(this->interp, this->reshape_command.c_str(),
	    reshape_callback, this, NULL);

	if (toplevel) {
	    if (!this->eval({"toplevel", this->container_path.c_str()}))
		return this->open_failed();
	    this->container = Tk_NameToWindow(this->interp,
		this->container_path.c_str(), Tk_MainWindow(this->interp));
	    if (!this->container)
		return this->open_failed();
	    if (!this->eval({"rename", this->container_path.c_str(),
		this->container_command_name.c_str()}))
		return this->open_failed();
	}

	const std::string swidth = std::to_string(this->width);
	const std::string sheight = std::to_string(this->height);
	if (this->software) {
	    if (!this->eval({"image", "create", "photo", this->photo_name.c_str(),
		"-width", swidth.c_str(), "-height", sheight.c_str()}) ||
		!this->eval({"label", this->widget_path.c_str(), "-image",
		    this->photo_name.c_str(), "-borderwidth", "0",
		    "-highlightthickness", "0"}))
		return this->open_failed();
	    this->photo = Tk_FindPhoto(this->interp, this->photo_name.c_str());
	    if (!this->photo)
		return this->open_failed();
	} else {
	    if (!this->eval({"::brlcad::tkobol::host", this->widget_path.c_str(),
		"-width", swidth.c_str(), "-height", sheight.c_str(),
		"-rgba", "true", "-double", "true", "-depth", "true",
		"-depthsize", "24", "-major", "2", "-minor", "1",
		"-coreprofile", "false", "-swapinterval", "1",
		"-displayproc", this->display_command.c_str(),
		"-reshapeproc", this->reshape_command.c_str()}))
		return this->open_failed();
	    /* TclCAD may use the view path as its own command name.  Preserve
	     * Togl's widget command before that later registration replaces it. */
	    if (!this->eval({"rename", this->widget_path.c_str(),
		this->widget_command_name.c_str()}))
		return this->open_failed();
	}

	this->tkwin = Tk_NameToWindow(this->interp, this->widget_path.c_str(),
	    Tk_MainWindow(this->interp));
	if (!this->tkwin)
	    return this->open_failed();
	if (!this->container)
	    this->container = this->tkwin;
	if (toplevel && !this->eval({"pack", this->widget_path.c_str(),
	    "-expand", "true", "-fill", "both"}))
	    return this->open_failed();

	Tk_MapWindow(this->tkwin);
	if (toplevel) {
	    Tk_MakeWindowExist(this->container);
	    Tk_MapWindow(this->container);
	    if (desc->title && desc->title[0])
		(void)this->eval({"wm", "title", this->container_path.c_str(),
		    desc->title});
	    if (desc->visible)
		(void)this->eval({"wm", "deiconify", this->container_path.c_str()});
	}

	this->opened = true;
	this->closing = false;
	return this->request("Tk endpoint open");
    }

    void close(void)
    {
	if (this->closing)
	    return;
	this->closing = true;
	if (Tcl_GetCurrentThread() == this->owner_thread)
	    Tcl_DeleteEvents(delete_event, this);
	this->event_queued.store(false);

	if (this->container)
	    Tk_DestroyWindow(this->container);
	this->container = NULL;
	this->tkwin = NULL;
	if (this->interp && this->software && !this->photo_name.empty())
	    (void)this->eval({"image", "delete", this->photo_name.c_str()});
	if (this->interp && !this->display_command.empty())
	    Tcl_DeleteCommand(this->interp, this->display_command.c_str());
	if (this->interp && !this->reshape_command.empty())
	    Tcl_DeleteCommand(this->interp, this->reshape_command.c_str());
	this->photo = NULL;
	this->opened = false;
	this->controller = NULL;
    }

    int request(const char *reason)
    {
	if (!this->controller || this->closing)
	    return 0;
	this->controller->requestRender(reason);
	bool expected = false;
	if (!this->event_queued.compare_exchange_strong(expected, true))
	    return 1;
	TkEndpointEvent *event = reinterpret_cast<TkEndpointEvent *>(
	    Tcl_Alloc(sizeof(TkEndpointEvent)));
	event->header.proc = event_proc;
	event->host = this;
	Tcl_ThreadQueueEvent(this->owner_thread, &event->header, TCL_QUEUE_TAIL);
	Tcl_ThreadAlert(this->owner_thread);
	return 1;
    }

    int resize(unsigned int new_width, unsigned int new_height,
	double device_pixel_ratio)
    {
	if (!this->opened || !this->tkwin || !new_width || !new_height ||
	    device_pixel_ratio <= 0.0)
	    return 0;
	this->width = new_width;
	this->height = new_height;
	Tk_GeometryRequest(this->tkwin,
	    std::max(1, (int)(new_width / device_pixel_ratio)),
	    std::max(1, (int)(new_height / device_pixel_ratio)));
	this->controller->setViewportSize(new_width, new_height);
	return this->request("Tk endpoint resize");
    }

    int capture(unsigned char **pixels, size_t *size, unsigned int *out_width,
	unsigned int *out_height, unsigned int *components)
    {
	if (!this->opened || !this->controller || !pixels || !size ||
	    !out_width || !out_height || !components)
	    return 0;
	*pixels = NULL;
	if (this->software) {
	    if (this->controller->renderToImage(pixels, 0, 1, NULL,
		brlobol_headless_context_manager(), NULL) != BRLCAD_OK || !*pixels)
		return 0;
	} else {
	    if (!this->render_hardware(false) ||
		!this->widget_command("makecurrent"))
		return 0;
	    const size_t byte_count = (size_t)this->width * this->height * 4;
	    *pixels = static_cast<unsigned char *>(bu_malloc(byte_count,
		"Tk endpoint capture"));
	    GLboolean double_buffered = GL_FALSE;
	    GLint alignment = 4;
	    GLint read_buffer = GL_BACK;
	    glGetBooleanv(GL_DOUBLEBUFFER, &double_buffered);
	    glGetIntegerv(GL_PACK_ALIGNMENT, &alignment);
	    glGetIntegerv(GL_READ_BUFFER, &read_buffer);
	    glPixelStorei(GL_PACK_ALIGNMENT, 1);
	    glReadBuffer(double_buffered ? GL_BACK : GL_FRONT);
	    glReadPixels(0, 0, (GLsizei)this->width, (GLsizei)this->height,
		GL_RGBA, GL_UNSIGNED_BYTE, *pixels);
	    glFinish();
	    glReadBuffer(read_buffer);
	    glPixelStorei(GL_PACK_ALIGNMENT, alignment);
	}
	*out_width = this->width;
	*out_height = this->height;
	*components = 4;
	*size = (size_t)this->width * this->height * 4;
	return 1;
    }

    int dimensions(unsigned int *out_width, unsigned int *out_height,
	double *device_pixel_ratio)
    {
	if (!this->opened || !this->tkwin || !out_width || !out_height ||
	    !device_pixel_ratio)
	    return 0;
	this->sync_size();
	*out_width = this->width;
	*out_height = this->height;
	*device_pixel_ratio = 1.0;
	return 1;
    }

private:
    bool eval(std::initializer_list<const char *> args)
    {
	if (!this->interp)
	    return false;
	std::vector<Tcl_Obj *> objv;
	objv.reserve(args.size());
	for (const char *arg : args) {
	    Tcl_Obj *obj = Tcl_NewStringObj(arg ? arg : "", -1);
	    Tcl_IncrRefCount(obj);
	    objv.push_back(obj);
	}
	Tcl_InterpState state = Tcl_SaveInterpState(this->interp, TCL_OK);
	const int ret = Tcl_EvalObjv(this->interp, (int)objv.size(), objv.data(),
	    TCL_EVAL_GLOBAL);
	if (ret != TCL_OK)
	    bu_log("Tk Obol endpoint command failed: %s\n",
		Tcl_GetStringResult(this->interp));
	Tcl_RestoreInterpState(this->interp, state);
	for (std::vector<Tcl_Obj *>::reverse_iterator it = objv.rbegin();
	    it != objv.rend(); ++it)
	    Tcl_DecrRefCount(*it);
	return ret == TCL_OK;
    }

    bool widget_command(const char *command)
    {
	const char *widget_command = this->software ? this->widget_path.c_str() :
	    this->widget_command_name.c_str();
	return this->eval({widget_command, command});
    }

    void sync_size(void)
    {
	if (!this->tkwin || !this->controller)
	    return;
	const int tk_width = Tk_Width(this->tkwin);
	const int tk_height = Tk_Height(this->tkwin);
	if (tk_width > 1)
	    this->width = (unsigned int)tk_width;
	if (tk_height > 1)
	    this->height = (unsigned int)tk_height;
	this->controller->setViewportSize(this->width, this->height);
    }

    bool present_software(void)
    {
	this->sync_size();
	unsigned char *pixels = NULL;
	BRLObolProgressiveStatus progressive;
	if (this->controller->renderToImage(&pixels, 1, 0, NULL,
		brlobol_headless_context_manager(), &progressive) != BRLCAD_OK ||
	    !pixels)
	    return false;

	Tk_PhotoImageBlock block;
	block.pixelPtr = pixels;
	block.width = (int)this->width;
	block.height = (int)this->height;
	block.pitch = (int)this->width * 3;
	block.pixelSize = 3;
	block.offset[0] = 0;
	block.offset[1] = 1;
	block.offset[2] = 2;
	block.offset[3] = 0;
#if TK_MAJOR_VERSION < 9
	Tk_PhotoSetSize(this->interp, this->photo, block.width, block.height);
	Tk_PhotoPutBlock(this->interp, this->photo, &block, 0, 0,
	    block.width, block.height, TK_PHOTO_COMPOSITE_SET);
#else
	Tk_PhotoSetSize(this->photo, block.width, block.height);
	Tk_PhotoPutBlock(this->photo, &block, 0, 0, block.width, block.height,
	    TK_PHOTO_COMPOSITE_SET);
#endif
	bu_free(pixels, "Tk endpoint software frame");
	if (progressive.hasMore || this->controller->hasProgressiveWorkPending())
	    (void)this->request("Tk endpoint progressive");
	return true;
    }

    bool render_hardware(bool present)
    {
	if (!this->widget_command("makecurrent"))
	    return false;
	this->sync_size();
	GLboolean double_buffered = GL_FALSE;
	glGetBooleanv(GL_DOUBLEBUFFER, &double_buffered);
	this->controller->getRenderManager()->setDoubleBuffer(
	    double_buffered ? TRUE : FALSE);
	glDrawBuffer(double_buffered ? GL_BACK : GL_FRONT);
	(void)this->controller->realizePending();
	(void)this->controller->advanceProgressiveWork();
	this->controller->renderBackground();
	const uint64_t started = this->controller->beginRenderTiming();
	if (this->controller->getCamera() && this->controller->getRenderRoot())
	    this->controller->getRenderManager()->render(FALSE, FALSE);
	this->controller->completeRenderTiming(started);
	this->controller->clearRenderRequest();
	if (present && double_buffered && !this->widget_command("swapbuffers"))
	    return false;
	return true;
    }

    bool render(void)
    {
	if (!this->opened || this->closing || !this->controller)
	    return false;
	return this->software ? this->present_software() :
	    this->render_hardware(true);
    }

    int open_failed(void)
    {
	this->close();
	return 0;
    }

    static int display_callback(ClientData data, Tcl_Interp *interp,
	int objc, Tcl_Obj *const *objv)
    {
	TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(data);
	if (!host || objc != 2) {
	    Tcl_WrongNumArgs(interp, 1, objv, "widget");
	    return TCL_ERROR;
	}
	return host->render() ? TCL_OK : TCL_ERROR;
    }

    static int reshape_callback(ClientData data, Tcl_Interp *interp,
	int objc, Tcl_Obj *const *objv)
    {
	TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(data);
	if (!host || objc != 2) {
	    Tcl_WrongNumArgs(interp, 1, objv, "widget");
	    return TCL_ERROR;
	}
	host->sync_size();
	return host->request("Tk endpoint reshape") ? TCL_OK : TCL_ERROR;
    }

    static int event_proc(Tcl_Event *event_ptr, int UNUSED(flags))
    {
	TkEndpointEvent *event = reinterpret_cast<TkEndpointEvent *>(event_ptr);
	TkObolEndpointHost *host = event ? event->host : NULL;
	if (!host)
	    return 1;
	host->event_queued.store(false);
	if (!host->software && host->opened && !host->closing)
	    (void)host->widget_command("postredisplay");
	else if (host->software)
	    (void)host->render();
	return 1;
    }

    static int delete_event(Tcl_Event *event_ptr, ClientData data)
    {
	TkEndpointEvent *event = reinterpret_cast<TkEndpointEvent *>(event_ptr);
	return event && event->host == data;
    }

    Tcl_Interp *interp;
    Tcl_ThreadId owner_thread;
    bool software;
    BRLObolViewController *controller;
    Tk_Window tkwin;
    Tk_Window container;
    Tk_PhotoHandle photo;
    bool opened;
    bool closing;
    std::atomic<bool> event_queued;
    unsigned int width;
    unsigned int height;
    std::string container_path;
    std::string container_command_name;
    std::string widget_path;
    std::string display_command;
    std::string reshape_command;
    std::string widget_command_name;
    std::string photo_name;
};

struct TkFactoryKind {
    bool software;
};

static int
tk_factory_probe(const struct brlobol_host_desc *desc, void *UNUSED(data))
{
    Tcl_Interp *interp = desc ?
	static_cast<Tcl_Interp *>(desc->application_context) : NULL;
    if (!interp || !Tk_MainWindow(interp) ||
	(desc->mode != BRLOBOL_HOST_MODE_TOPLEVEL &&
	 desc->mode != BRLOBOL_HOST_MODE_EMBEDDED))
	return 0;
    return desc->mode != BRLOBOL_HOST_MODE_EMBEDDED ||
	(desc->native_id_hint && desc->native_id_hint[0]);
}

static void *
tk_factory_create(const struct brlobol_host_desc *desc, void *data)
{
    TkFactoryKind *kind = static_cast<TkFactoryKind *>(data);
    Tcl_Interp *interp = desc ?
	static_cast<Tcl_Interp *>(desc->application_context) : NULL;
    return kind && interp ?
	new (std::nothrow) TkObolEndpointHost(interp, kind->software) : NULL;
}

static void
tk_factory_destroy(void *instance, void *UNUSED(data))
{
    delete static_cast<TkObolEndpointHost *>(instance);
}

static int
tk_factory_bind(void *instance, void *controller, void *UNUSED(data))
{
    TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(instance);
    return host ? host->bind(static_cast<BRLObolViewController *>(controller)) : 0;
}

static int
tk_factory_open(void *instance, const struct brlobol_host_desc *desc,
	void *UNUSED(data))
{
    TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(instance);
    return host ? host->open(desc) : 0;
}

static void
tk_factory_close(void *instance, void *UNUSED(data))
{
    TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(instance);
    if (host)
	host->close();
}

static int
tk_factory_request(void *instance, const char *reason, void *UNUSED(data))
{
    TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(instance);
    return host ? host->request(reason) : 0;
}

static int
tk_factory_resize(void *instance, unsigned int width, unsigned int height,
	double device_pixel_ratio, void *UNUSED(data))
{
    TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(instance);
    return host ? host->resize(width, height, device_pixel_ratio) : 0;
}

static int
tk_factory_capture(void *instance, unsigned char **pixels, size_t *size,
	unsigned int *width, unsigned int *height, unsigned int *components,
	void *UNUSED(data))
{
    TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(instance);
    return host ? host->capture(pixels, size, width, height, components) : 0;
}

static int
tk_factory_dimensions(void *instance, unsigned int *width,
	unsigned int *height, double *device_pixel_ratio, void *UNUSED(data))
{
    TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(instance);
    return host ? host->dimensions(width, height, device_pixel_ratio) : 0;
}

static brlobol_host_factory_token_t *
tk_factory_register(const char *name, int priority, uint64_t capabilities,
	TkFactoryKind *kind)
{
    struct brlobol_host_factory factory;
    std::memset(&factory, 0, sizeof(factory));
    factory.abi_version = BRLOBOL_HOST_FACTORY_ABI_VERSION;
    factory.struct_size = sizeof(factory);
    factory.name = name;
    factory.priority = priority;
    factory.capabilities = capabilities;
    factory.user_data = kind;
    factory.probe = tk_factory_probe;
    factory.create = tk_factory_create;
    factory.destroy = tk_factory_destroy;
    factory.bind_controller = tk_factory_bind;
    factory.open = tk_factory_open;
    factory.close = tk_factory_close;
    factory.request_frame = tk_factory_request;
    factory.resize = tk_factory_resize;
    factory.capture = tk_factory_capture;
    factory.dimensions = tk_factory_dimensions;
    return brlobol_host_factory_register(&factory);
}

extern "C" TCLCAD_EXPORT int
tclcad_obol_host_factories_register(void)
{
    static TkFactoryKind gl_kind = {false};
    static TkFactoryKind photo_kind = {true};
    const uint64_t common = BRLOBOL_HOST_CAP_TOPLEVEL |
	BRLOBOL_HOST_CAP_EMBEDDED | BRLOBOL_HOST_CAP_PROGRESSIVE_PRESENT |
	BRLOBOL_HOST_CAP_INPUT | BRLOBOL_HOST_CAP_READBACK |
	BRLOBOL_HOST_CAP_FRAMEBUFFER_PRESENT | BRLOBOL_HOST_CAP_MULTI_VIEW |
	BRLOBOL_HOST_CAP_THREAD_AFFINE;
    static brlobol_host_factory_token_t *gl_token = tk_factory_register(
	"tk-gl", 40, common | BRLOBOL_HOST_CAP_SYSTEM_GL, &gl_kind);
    static brlobol_host_factory_token_t *photo_token = tk_factory_register(
	"tk-photo", 30, common | BRLOBOL_HOST_CAP_PIXEL_PRESENT, &photo_kind);
    return gl_token && photo_token ? 1 : 0;
}
