/*            T K O B O L _ T H R E A D _ A F F I N I T Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"

#include "brlobol/display_endpoint.h"
#include "brlobol/view_controller.h"
#include "bu/log.h"
#include "bu/vls.h"
#include "tclcad/setup.h"

#include <Inventor/SbColor.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>

#include <tcl.h>
#include <tk.h>

#include <chrono>
#include <thread>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\\n", _msg); \
	return 1; \
    } \
} while (0)

static int
update_tk_events(Tcl_Interp *interp)
{
    const int ret = Tcl_EvalEx(interp, "update", -1, TCL_EVAL_DIRECT);
    if (ret != TCL_OK)
	Tcl_ResetResult(interp);
    return ret == TCL_OK ? 1 : 0;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    Tcl_FindExecutable(argc > 0 ? argv[0] : "tkobol_thread_affinity");
    Tcl_Interp *interp = Tcl_CreateInterp();
    CHECK(interp != NULL, "Tk affinity test creates an interpreter");
    struct bu_vls tcl_log = BU_VLS_INIT_ZERO;
    if (tclcad_init(interp, 1, &tcl_log) != TCL_OK) {
	bu_log("FAIL: Tk affinity test initializes Tcl/Tk: %s\\n",
	    bu_vls_addr(&tcl_log));
	bu_vls_free(&tcl_log);
	return 1;
    }
    bu_vls_free(&tcl_log);
    CHECK(tclcad_obol_host_factories_register(),
	"Tk affinity test registers Obol host factories");

    struct brlobol_host_desc desc = {};
    desc.struct_size = sizeof(desc);
    desc.mode = BRLOBOL_HOST_MODE_TOPLEVEL;
    desc.width = 64;
    desc.height = 48;
    desc.device_pixel_ratio = 1.0;
    desc.visible = 1;
    desc.required_capabilities = BRLOBOL_HOST_CAP_PIXEL_PRESENT |
	BRLOBOL_HOST_CAP_THREAD_AFFINE;
    desc.application_context = interp;

    brlobol_display_endpoint_t *endpoint =
	brlobol_display_endpoint_create(NULL, 0);
    CHECK(endpoint != NULL, "Tk affinity test creates endpoint");
    CHECK(brlobol_display_endpoint_host_open(endpoint, "tk-photo", &desc),
	"Tk affinity test opens TkPhoto host");

    BRLObolViewController *controller =
	static_cast<BRLObolViewController *>(
	    brlobol_display_endpoint_controller(endpoint));
    CHECK(controller != NULL, "Tk affinity test obtains endpoint controller");
    SoPerspectiveCamera *camera = new SoPerspectiveCamera;
    camera->position = SbVec3f(0.0f, 0.0f, 1.0f);
    controller->setCamera(camera);

    CHECK(update_tk_events(interp), "Tk affinity test drains the host-open frame");
    controller->clearProgressiveWorkPending();
    controller->setBackgroundColors(SbColor(0.8f, 0.0f, 0.0f),
	SbColor(0.8f, 0.0f, 0.0f));
    controller->clearRenderRequest();
    const uint64_t render_time_before =
	controller->getLastRenderTimeNanoseconds();

    int request_ok = 0;
    std::thread worker([endpoint, &request_ok]() {
	request_ok = brlobol_display_endpoint_request_frame(endpoint,
	    "worker-thread-frame");
    });
    worker.join();

    const std::chrono::steady_clock::time_point deadline =
	std::chrono::steady_clock::now() + std::chrono::seconds(2);
    int dispatched = 0;
    while (std::chrono::steady_clock::now() < deadline) {
	CHECK(update_tk_events(interp), "Tk affinity test services the Tcl event loop");
	if (!controller->isRenderRequested() &&
	    controller->getLastRenderTimeNanoseconds() != render_time_before) {
	    dispatched = 1;
	    break;
	}
	std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }

    CHECK(request_ok && dispatched,
	"TkPhoto host renders a worker-requested frame on the Tcl event loop");

    controller->clearProgressiveWorkPending();
    brlobol_display_endpoint_destroy(endpoint);
    (void)update_tk_events(interp);
    Tcl_DeleteInterp(interp);
    return 0;
}
