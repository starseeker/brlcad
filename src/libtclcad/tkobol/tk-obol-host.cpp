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
#  include <errno.h>
#  include <fcntl.h>
#  include <unistd.h>
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
#include <Inventor/actions/SoGLRenderAction.h>

#include "vendor/togl/togl.h"

#include <algorithm>
#include <atomic>
#include <cstring>
#include <initializer_list>
#include <new>
#include <string>
#include <vector>

static thread_local unsigned int tk_endpoint_offscreen_depth = 0;

#if defined(TKOBOL_WGL) && defined(TCL_THREADS)
static Tcl_Mutex tk_endpoint_event_mutex;
#endif

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
	/* A null handle classifies the host's on-screen action, which is Togl
	 * system GL.  Only contexts explicitly created through the delegated
	 * offscreen manager are OSMesa contexts. */
	return context && this->offscreen ?
	    this->offscreen->isOSMesaContext(context) : FALSE;
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

#ifdef TKOBOL_X11
static bool
tk_endpoint_set_nonblocking(int fd)
{
    const int flags = fcntl(fd, F_GETFL, 0);
    return flags >= 0 && fcntl(fd, F_SETFL, flags | O_NONBLOCK) == 0;
}
#endif

class TkObolEndpointHost {
public:
    TkObolEndpointHost(Tcl_Interp *new_interp, bool use_software) :
	interp(new_interp), owner_thread(Tcl_GetCurrentThread()), software(use_software),
	controller(NULL), tkwin(NULL), container(NULL), photo(NULL), opened(false),
	closing(false), toplevel(false), visible(false), event_queued(false),
	input_handler_registered(false), container_key_handler_registered(false),
	input_dispatch(NULL),
	input_dispatch_data(NULL), last_motion_x(0), last_motion_y(0),
	last_motion_valid(false), width(1), height(1)
#ifdef TKOBOL_X11
	, request_pipe_read(-1), request_pipe_write(-1)
#endif
    {
    }

    ~TkObolEndpointHost(void)
    {
	this->close();
    }

    int bind(BRLObolViewController *new_controller)
    {
	this->controller = new_controller;
	SoRenderManager *manager = this->controller ?
	    this->controller->getRenderManager() : NULL;
	SoGLRenderAction *action = manager ? manager->getGLRenderAction() : NULL;
	if (action)
	    action->setContextManager(tk_endpoint_context_manager());
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
	this->input_dispatch = desc->input_dispatch;
	this->input_dispatch_data = desc->input_dispatch_data;
	this->controller->setViewportSize(this->width, this->height);

	static std::atomic<unsigned long> next_id(0);
	const unsigned long id = ++next_id;
	this->container_path = desc->native_id_hint && desc->native_id_hint[0] ?
	    desc->native_id_hint : ".brlobol_tk" + std::to_string(id);
	this->container_command_name = "::brlcad::tkobol::endpoint_container_" +
	    std::to_string(id);
	this->toplevel = desc->mode == BRLOBOL_HOST_MODE_TOPLEVEL;
	const bool is_toplevel = this->toplevel;
	if (!is_toplevel && (!desc->native_id_hint || !desc->native_id_hint[0]))
	    return 0;
	if (!this->setup_request_notifier())
	    return this->open_failed();
	this->widget_path = is_toplevel ? this->container_path + ".__obol" :
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

	if (is_toplevel) {
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
	const char *swap_interval = desc->vsync == BRLOBOL_HOST_VSYNC_OFF ?
	    "0" : "1";
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
		"-coreprofile", "false", "-swapinterval", swap_interval,
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
	#ifdef TKOBOL_X11
	if (this->input_dispatch) {
	    Tk_CreateEventHandler(this->tkwin,
		KeyPressMask | KeyReleaseMask | ButtonPressMask |
		ButtonReleaseMask | PointerMotionMask | ExposureMask |
		StructureNotifyMask | FocusChangeMask,
		input_event_handler, this);
	    /* Tk directs keyboard focus to a toplevel while pointer events target
	     * its child widget.  A generic handler observes the former without
	     * attaching an event record to Tk's renamed toplevel command. */
	    if (this->container != this->tkwin)
		Tk_CreateGenericHandler(container_key_event_handler, this);
	    this->container_key_handler_registered =
		this->container != this->tkwin;
	    this->input_handler_registered = true;
	}
	#endif
	if (!this->software && desc->vsync != BRLOBOL_HOST_VSYNC_AUTO &&
	    !BrlcadTkObolHost_SetSwapInterval(this->interp,
		this->widget_command_name.c_str(),
		desc->vsync == BRLOBOL_HOST_VSYNC_OFF ? 0 : 1))
	    return this->open_failed();
	if (!this->container)
	    this->container = this->tkwin;
	if (is_toplevel && !this->eval({"pack", this->widget_path.c_str(),
	    "-expand", "true", "-fill", "both"}))
	    return this->open_failed();

	Tk_MapWindow(this->tkwin);
	if (is_toplevel) {
	    Tk_MakeWindowExist(this->container);
	    Tk_MapWindow(this->container);
	    if (desc->title && desc->title[0])
		(void)this->eval({"wm", "title", this->container_path.c_str(),
		    desc->title});
	    (void)this->eval({"wm", desc->visible ? "deiconify" : "withdraw",
		this->container_path.c_str()});
	}

	this->opened = true;
	this->visible = !is_toplevel || desc->visible != 0;
	this->closing = false;
	return this->visible ? this->request("Tk endpoint open") : 1;
    }

    void close(void)
    {
	if (this->closing)
	    return;
	this->closing = true;
	if (Tcl_GetCurrentThread() == this->owner_thread)
	    Tcl_DeleteEvents(delete_event, this);
	this->event_queued.store(false);
	this->teardown_request_notifier();
	#ifdef TKOBOL_X11
	if (this->input_handler_registered) {
	    if (this->tkwin)
		Tk_DeleteEventHandler(this->tkwin,
		    KeyPressMask | KeyReleaseMask | ButtonPressMask |
		    ButtonReleaseMask | PointerMotionMask | ExposureMask |
		    StructureNotifyMask | FocusChangeMask,
		    input_event_handler, this);
	}
	if (this->container_key_handler_registered)
	    Tk_DeleteGenericHandler(container_key_event_handler, this);
	#endif
	this->input_handler_registered = false;
	this->container_key_handler_registered = false;
	this->input_dispatch = NULL;
	this->input_dispatch_data = NULL;
	this->last_motion_valid = false;

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
	this->toplevel = false;
	this->visible = false;
	this->controller = NULL;
    }

    int request(const char *reason)
    {
	if (!this->controller || this->closing)
	    return 0;
	if (!this->visible)
	    return 1;
	this->controller->requestRender(reason);
	bool expected = false;
	if (!this->event_queued.compare_exchange_strong(expected, true))
	    return 1;
	return this->queue_frame_event();
    }

    int queue_frame_event(void)
    {
#ifdef TKOBOL_X11
	if (this->request_pipe_write < 0) {
	    this->event_queued.store(false);
	    return 0;
	}
	char marker = 0;
	ssize_t written = -1;
	do {
	    written = write(this->request_pipe_write, &marker, 1);
	} while (written < 0 && errno == EINTR);
	if (written == 1 || (written < 0 &&
		(errno == EAGAIN || errno == EWOULDBLOCK)))
	    return 1;
	this->event_queued.store(false);
	return 0;
#else
	TkEndpointEvent *event = reinterpret_cast<TkEndpointEvent *>(
	    Tcl_Alloc(sizeof(TkEndpointEvent)));
	event->header.proc = event_proc;
	event->host = this;
	if (Tcl_GetCurrentThread() == this->owner_thread) {
	    Tcl_QueueEvent(&event->header, TCL_QUEUE_TAIL);
	    return 1;
	}
#  ifdef TCL_THREADS
	Tcl_MutexLock(&tk_endpoint_event_mutex);
	Tcl_ThreadQueueEvent(this->owner_thread, &event->header, TCL_QUEUE_TAIL);
	Tcl_ThreadAlert(this->owner_thread);
	Tcl_MutexUnlock(&tk_endpoint_event_mutex);
	return 1;
#  else
	Tcl_Free(reinterpret_cast<char *>(event));
	this->event_queued.store(false);
	return 0;
#  endif
#endif
	}

    void dispatch_frame_event(void)
    {
	this->event_queued.store(false);
	if (!this->opened || this->closing)
	    return;
	if (!this->software)
	    (void)this->widget_command("postredisplay");
	else
	    (void)this->render();
	}

    bool setup_request_notifier(void)
    {
#ifdef TKOBOL_X11
	if (this->request_pipe_read >= 0 || this->request_pipe_write >= 0)
	    return true;
	int fds[2] = {-1, -1};
	if (pipe(fds) != 0 || !tk_endpoint_set_nonblocking(fds[0]) ||
	    !tk_endpoint_set_nonblocking(fds[1])) {
	    if (fds[0] >= 0)
		::close(fds[0]);
	    if (fds[1] >= 0)
		::close(fds[1]);
	    return false;
	}
	this->request_pipe_read = fds[0];
	this->request_pipe_write = fds[1];
	Tcl_CreateFileHandler(this->request_pipe_read, TCL_READABLE,
	    request_pipe_handler, this);
#endif
	return true;
	}

    void teardown_request_notifier(void)
    {
#ifdef TKOBOL_X11
	if (this->request_pipe_read >= 0) {
	    Tcl_DeleteFileHandler(this->request_pipe_read);
	    ::close(this->request_pipe_read);
	}
	if (this->request_pipe_write >= 0)
	    ::close(this->request_pipe_write);
	this->request_pipe_read = -1;
	this->request_pipe_write = -1;
#endif
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

    int set_title(const char *title)
    {
	if (!this->opened || !this->toplevel || !title ||
	    Tcl_GetCurrentThread() != this->owner_thread)
	    return 0;
	return this->eval({"wm", "title", this->container_path.c_str(), title}) ?
	    1 : 0;
    }

    int set_visible(int mapped)
    {
	if (!this->opened || !this->toplevel ||
	    Tcl_GetCurrentThread() != this->owner_thread)
	    return 0;
	if (!mapped) {
	    /* Togl may have a toolkit redisplay idle callback in addition to the
	     * endpoint's Tcl_Event.  Drain it while the GL context is still mapped
	     * before suppressing future frame requests and withdrawing the window. */
	    (void)this->eval({"update", "idletasks"});
	    this->visible = false;
	    Tcl_DeleteEvents(delete_event, this);
	    this->event_queued.store(false);
	    if (this->controller)
		this->controller->clearRenderRequest();
	    return this->eval({"wm", "withdraw",
		this->container_path.c_str()}) ? 1 : 0;
	}
	if (!this->eval({"wm", "deiconify", this->container_path.c_str()}))
	    return 0;
	this->visible = true;
	return this->request("Tk endpoint visible");
    }

    int set_vsync(int enabled)
    {
	if (!this->opened || this->software ||
	    Tcl_GetCurrentThread() != this->owner_thread)
	    return 0;
	return BrlcadTkObolHost_SetSwapInterval(this->interp,
	    this->widget_command_name.c_str(), enabled ? 1 : 0);
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
	if (!this->opened || this->closing || !this->controller ||
	    !this->visible)
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
	if (!host->opened || host->closing || !host->visible)
	    return TCL_OK;
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

#ifdef TKOBOL_X11
    static unsigned int input_modifiers(unsigned int state)
    {
	unsigned int modifiers = BRLOBOL_INPUT_MOD_NONE;
	if (state & ShiftMask)
	    modifiers |= BRLOBOL_INPUT_MOD_SHIFT;
	if (state & ControlMask)
	    modifiers |= BRLOBOL_INPUT_MOD_CONTROL;
	if (state & Mod1Mask)
	    modifiers |= BRLOBOL_INPUT_MOD_ALT;
	if (state & Mod4Mask)
	    modifiers |= BRLOBOL_INPUT_MOD_META;
	return modifiers;
    }

    static unsigned int input_buttons(unsigned int state)
    {
	unsigned int buttons = 0;
	if (state & Button1Mask)
	    buttons |= 1u << 0;
	if (state & Button2Mask)
	    buttons |= 1u << 1;
	if (state & Button3Mask)
	    buttons |= 1u << 2;
	return buttons;
    }

    static int input_button(unsigned int state)
    {
	if (state & Button1Mask)
	    return 0;
	if (state & Button2Mask)
	    return 1;
	if (state & Button3Mask)
	    return 2;
	return BRLOBOL_INPUT_ANY;
    }

    static int input_button_number(unsigned int button)
    {
	return button >= Button1 && button <= Button3 ?
	    static_cast<int>(button - Button1) : BRLOBOL_INPUT_ANY;
    }

    void dispatch_input(const BRLObolInputEvent &event)
    {
	if (!this->input_dispatch || this->closing)
	    return;
	if (this->input_dispatch(this->input_dispatch_data, &event) > 0)
	    (void)this->request("Tk input action");
    }

    static void input_event_handler(ClientData data, XEvent *xevent)
    {
	TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(data);
	if (!host || !xevent || !host->input_dispatch || host->closing)
	    return;

	BRLObolInputEvent event;
	switch (xevent->type) {
	    case KeyPress:
	    case KeyRelease: {
		event.type = xevent->type == KeyPress ? BRLOBOL_INPUT_KEY_PRESS :
	    BRLOBOL_INPUT_KEY_RELEASE;
		event.timestamp = xevent->xkey.time;
		event.x = xevent->xkey.x;
		event.y = xevent->xkey.y;
		event.modifiers = input_modifiers(xevent->xkey.state);
		KeySym key = XLookupKeysym(&xevent->xkey, 0);
		if (key >= XK_a && key <= XK_z)
		    key -= XK_a - XK_A;
		event.key = static_cast<int>(key);
		host->dispatch_input(event);
		return;
	    }
	    case ButtonPress:
	    case ButtonRelease: {
		if (xevent->type == ButtonPress &&
		    (xevent->xbutton.button == Button4 ||
		     xevent->xbutton.button == Button5)) {
		    event.type = BRLOBOL_INPUT_WHEEL;
		    event.timestamp = xevent->xbutton.time;
		    event.x = xevent->xbutton.x;
		    event.y = xevent->xbutton.y;
		    event.modifiers = input_modifiers(xevent->xbutton.state);
		    /* Match Qt's angleDelta()/8 normalization. */
		    event.wheelDelta = xevent->xbutton.button == Button4 ? -15 : 15;
		    host->dispatch_input(event);
		    return;
		}
		event.type = xevent->type == ButtonPress ? BRLOBOL_INPUT_POINTER_PRESS :
	    BRLOBOL_INPUT_POINTER_RELEASE;
		event.timestamp = xevent->xbutton.time;
		event.x = xevent->xbutton.x;
		event.y = xevent->xbutton.y;
		event.button = input_button_number(xevent->xbutton.button);
		event.buttons = input_buttons(xevent->xbutton.state);
		if (xevent->type == ButtonPress && event.button != BRLOBOL_INPUT_ANY)
		    event.buttons |= 1u << event.button;
	if (xevent->type == ButtonRelease && event.button != BRLOBOL_INPUT_ANY)
	    event.buttons &= ~(1u << event.button);
	event.modifiers = input_modifiers(xevent->xbutton.state);
	if (xevent->type == ButtonPress) {
	    host->last_motion_x = event.x;
	    host->last_motion_y = event.y;
	    host->last_motion_valid = true;
	} else {
	    host->last_motion_valid = false;
	}
	host->dispatch_input(event);
		return;
	    }
	    case MotionNotify:
		event.type = BRLOBOL_INPUT_POINTER_MOTION;
		event.timestamp = xevent->xmotion.time;
		event.x = xevent->xmotion.x;
		event.y = xevent->xmotion.y;
		event.buttons = input_buttons(xevent->xmotion.state);
		event.button = input_button(xevent->xmotion.state);
		event.modifiers = input_modifiers(xevent->xmotion.state);
		if (host->last_motion_valid) {
		    event.dx = event.x - host->last_motion_x;
		    event.dy = event.y - host->last_motion_y;
		}
		host->last_motion_x = event.x;
		host->last_motion_y = event.y;
		host->last_motion_valid = true;
		host->dispatch_input(event);
		return;
	    case ConfigureNotify:
		event.type = BRLOBOL_INPUT_RESIZE;
		event.width = static_cast<unsigned int>(std::max(1,
		    xevent->xconfigure.width));
		event.height = static_cast<unsigned int>(std::max(1,
		    xevent->xconfigure.height));
		host->dispatch_input(event);
		return;
	    case FocusIn:
	    case FocusOut:
		event.type = BRLOBOL_INPUT_FOCUS;
		event.key = xevent->type == FocusIn ? 1 : 0;
		host->dispatch_input(event);
		return;
	    case Expose:
		event.type = BRLOBOL_INPUT_EXPOSE;
		host->dispatch_input(event);
		return;
	default:
	    return;
	}
    }

    static int container_key_event_handler(ClientData data, XEvent *xevent)
    {
	TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(data);
	if (!host || !xevent || !host->container || host->closing ||
	    (xevent->type != KeyPress && xevent->type != KeyRelease) ||
	    xevent->xany.window != Tk_WindowId(host->container))
	    return 0;
	input_event_handler(data, xevent);
	return 0;
    }
#endif

    static int event_proc(Tcl_Event *event_ptr, int UNUSED(flags))
    {
	TkEndpointEvent *event = reinterpret_cast<TkEndpointEvent *>(event_ptr);
	TkObolEndpointHost *host = event ? event->host : NULL;
	if (!host)
	    return 1;
	host->dispatch_frame_event();
	return 1;
    }

#ifdef TKOBOL_X11
    static void request_pipe_handler(ClientData data, int UNUSED(mask))
    {
	TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(data);
	if (!host || host->request_pipe_read < 0)
	    return;
	char buffer[64];
	while (read(host->request_pipe_read, buffer, sizeof(buffer)) > 0) {
	    /* Drain all coalesced worker requests before rendering once. */
	}
	host->dispatch_frame_event();
    }
#endif

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
    bool toplevel;
    bool visible;
    std::atomic<bool> event_queued;
    bool input_handler_registered;
    bool container_key_handler_registered;
    BRLObolInputEventHandler input_dispatch;
    void *input_dispatch_data;
    int last_motion_x;
    int last_motion_y;
    bool last_motion_valid;
    unsigned int width;
    unsigned int height;
#ifdef TKOBOL_X11
    int request_pipe_read;
    int request_pipe_write;
#endif
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

static int
tk_factory_set_title(void *instance, const char *title, void *UNUSED(data))
{
    TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(instance);
    return host ? host->set_title(title) : 0;
}

static int
tk_factory_set_visible(void *instance, int visible, void *UNUSED(data))
{
    TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(instance);
    return host ? host->set_visible(visible) : 0;
}

static int
tk_factory_set_vsync(void *instance, int enabled, void *UNUSED(data))
{
    TkObolEndpointHost *host = static_cast<TkObolEndpointHost *>(instance);
    return host ? host->set_vsync(enabled) : 0;
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
    factory.set_title = tk_factory_set_title;
    factory.set_visible = tk_factory_set_visible;
    factory.set_vsync = tk_factory_set_vsync;
    return brlobol_host_factory_register(&factory);
}

extern "C" TCLCAD_EXPORT int
tclcad_obol_host_factories_register(void)
{
    static TkFactoryKind gl_kind = {false};
    static TkFactoryKind photo_kind = {true};
    uint64_t common = BRLOBOL_HOST_CAP_TOPLEVEL |
	BRLOBOL_HOST_CAP_EMBEDDED | BRLOBOL_HOST_CAP_PROGRESSIVE_PRESENT |
	BRLOBOL_HOST_CAP_READBACK |
	BRLOBOL_HOST_CAP_FRAMEBUFFER_PRESENT | BRLOBOL_HOST_CAP_MULTI_VIEW;
#if defined(TKOBOL_X11) || defined(TCL_THREADS)
    common |= BRLOBOL_HOST_CAP_THREAD_AFFINE;
#endif
#ifdef TKOBOL_X11
    common |= BRLOBOL_HOST_CAP_INPUT;
#endif
    static brlobol_host_factory_token_t *gl_token = tk_factory_register(
	"tk-gl", 40, common | BRLOBOL_HOST_CAP_SYSTEM_GL |
	BRLOBOL_HOST_CAP_PRESENT_VSYNC, &gl_kind);
    static brlobol_host_factory_token_t *photo_token = tk_factory_register(
	"tk-photo", 30, common | BRLOBOL_HOST_CAP_PIXEL_PRESENT, &photo_kind);
    return gl_token && photo_token ? 1 : 0;
}
