/*        T E S T _ Q G E D _ O B O L _ S E T T I N G S . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "bu/app.h"

#include "bu/str.h"

#include "BObol/BDisplayEndpoint.h"
#include "bu/log.h"
#include "ged.h"
#include "ged/view.h"
#include "qtcad/QgPluginContext.h"
#include "qtcad/QgView.h"

#include "CADViewSettings.h"

#include <QApplication>

#include <cstring>

#define CHECK(_expr, _msg) do { \
    if (!(_expr)) { \
	bu_log("FAIL: %s\n", _msg); \
	fail = 1; \
	goto cleanup; \
    } \
} while (0)

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    QApplication app(argc, argv);
    int fail = 0;
    struct ged *gedp = NULL;

    {
	QgView view(NULL, QgViewType::SW);
	QgPluginContext context;
	CADViewSettings settings;
	struct bv_display_property_value value =
	    BV_DISPLAY_PROPERTY_VALUE_INIT;
	struct ged_view_context *view_ctx = NULL;
	view.resize(160, 120);
	CHECK(view.viewContext() && view.displayEndpoint(),
	    "settings test creates an endpoint-backed QgView");
	view_ctx = ged_view_context_from_bv(view.viewContext());

	gedp = ged_open("inmem", "0x0", 0);
	CHECK(gedp, "settings test creates an in-memory GED context");
	ged_view_active_ctx_set(gedp, view_ctx);
	CHECK(ged_view_context_host_attach(gedp, view_ctx),
	    "settings test attaches the QgView context to GED");
	CHECK(ged_view_context_obol_endpoint_set(view_ctx,
	    view.displayEndpoint(), 0),
	    "settings test binds the QgView endpoint to its GED context");

	context.viewWidgetAccessor = [&view]() { return &view; };
	settings.setContext(&context);
	settings.checkbox_refresh(0);
	CHECK(settings.fb_mode_combo->isEnabled() &&
	    settings.fb_mode_combo->currentIndex() == 0,
	    "settings reads the endpoint-owned initial framebuffer mode");
	CHECK(settings.params_ckbx->isEnabled() &&
	    settings.params_ckbx->checkState() == Qt::Checked &&
	    settings.params_fps_ckbx->isEnabled() &&
	    settings.params_fps_ckbx->checkState() == Qt::Checked,
	    "parameter telemetry and FPS are effectively visible by default");

	settings.fb_mode_combo->setCurrentIndex(2);
	value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	CHECK(ged_view_context_display_property_get(view_ctx,
	    "composition.framebuffer.mode", &value) ==
	    BV_DISPLAY_PROPERTY_OK && value.string_value &&
	    bu_strcmp(value.string_value, "underlay") == 0,
	    "settings writes framebuffer composition through the endpoint");

	value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	value.type = BV_DISPLAY_PROPERTY_ENUM;
	value.string_value = "overlay";
	CHECK(ged_view_context_display_property_set(view_ctx,
	    "composition.framebuffer.mode", &value) ==
	    BV_DISPLAY_PROPERTY_OK,
	    "endpoint accepts an external framebuffer composition update");
	settings.checkbox_refresh(0);
	CHECK(settings.fb_mode_combo->isEnabled() &&
	    settings.fb_mode_combo->currentIndex() == 1,
	    "settings reflects the endpoint-owned framebuffer mode");

	settings.fb_mode_combo->setCurrentIndex(3);
	value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	CHECK(ged_view_context_display_property_get(view_ctx,
	    "composition.framebuffer.mode", &value) ==
	    BV_DISPLAY_PROPERTY_OK && value.string_value &&
	    bu_strcmp(value.string_value, "interlay") == 0,
	    "settings writes the Obol interlay framebuffer composition");

	settings.grid_ckbx->setCheckState(Qt::Checked);
	value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	CHECK(ged_view_context_display_property_get(view_ctx,
	    "view.faceplate.grid.visible", &value) ==
	    BV_DISPLAY_PROPERTY_OK && value.bool_value,
	    "settings writes faceplate visibility through the endpoint");

	value = BV_DISPLAY_PROPERTY_VALUE_INIT;
	value.type = BV_DISPLAY_PROPERTY_BOOL;
	value.bool_value = 1;
	CHECK(ged_view_context_display_property_set(view_ctx,
	    "view.faceplate.params.fps", &value) ==
	    BV_DISPLAY_PROPERTY_OK,
	    "endpoint accepts an external faceplate update");
	settings.checkbox_refresh(0);
	CHECK(settings.params_fps_ckbx->isEnabled() &&
	    settings.params_fps_ckbx->checkState() == Qt::Checked,
	    "settings reflects endpoint-owned faceplate visibility");

	CHECK(ged_view_context_obol_endpoint_set(view_ctx, NULL, 0),
	    "settings test detaches the endpoint property owner");
	settings.checkbox_refresh(0);
	CHECK(settings.fb_mode_combo->isEnabled() &&
	    settings.fb_mode_combo->currentIndex() == 3,
	    "settings retains framebuffer composition before endpoint attachment");

cleanup:
	if (gedp)
	    ged_close(gedp);
    }

    return fail;
}

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
