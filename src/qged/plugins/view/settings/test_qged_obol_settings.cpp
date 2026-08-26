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

#include <cmath>
#include <cstdio>

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
	const char *cutting_get[] = {"view", "cutting"};
	const char *cutting_external[] = {
	    "view", "cutting", "origin", "-4", "8", "2"
	};
	int cuttingEnabled = 0;
	double cuttingOrigin[3] = {};
	double cuttingNormal[3] = {};
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

	context.gedAccessor = [&gedp]() { return gedp; };
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
	CHECK(settings.cutting_grp->isEnabled() &&
	    !settings.cutting_enabled_ckbx->isChecked() &&
	    std::fabs(settings.cutting_normal[2]->value() - 1.0) < 1.0e-12,
	    "settings reads the disabled default cutting plane");

	settings.cutting_origin[0]->setValue(12.5);
	settings.cutting_origin[1]->setValue(-3.0);
	settings.cutting_origin[2]->setValue(7.25);
	settings.cutting_normal[0]->setValue(0.0);
	settings.cutting_normal[1]->setValue(1.0);
	settings.cutting_normal[2]->setValue(0.0);
	settings.cutting_enabled_ckbx->setChecked(true);
	CHECK(ged_exec(gedp, 2, cutting_get) == BRLCAD_OK &&
	    std::sscanf(bu_vls_cstr(gedp->ged_result_str),
	    "enable %d\norigin %lf %lf %lf\nnormal %lf %lf %lf",
	    &cuttingEnabled, &cuttingOrigin[0], &cuttingOrigin[1],
	    &cuttingOrigin[2], &cuttingNormal[0], &cuttingNormal[1],
	    &cuttingNormal[2]) == 7 && cuttingEnabled == 1 &&
	    /* A plane canonically retains normal * distance, so tangential
	     * components of the input point are intentionally absent on readback. */
	    std::fabs(cuttingOrigin[0]) < 1.0e-12 &&
	    std::fabs(cuttingOrigin[1] + 3.0) < 1.0e-12 &&
	    std::fabs(cuttingOrigin[2]) < 1.0e-12 &&
	    std::fabs(cuttingNormal[0]) < 1.0e-12 &&
	    std::fabs(cuttingNormal[1] - 1.0) < 1.0e-12 &&
	    std::fabs(cuttingNormal[2]) < 1.0e-12,
	    "settings writes the world-space cutting plane through GED");

	CHECK(ged_exec(gedp, 6, cutting_external) == BRLCAD_OK,
	    "external cutting origin command succeeds");
	settings.checkbox_refresh(0);
	CHECK(std::fabs(settings.cutting_origin[0]->value()) < 1.0e-12 &&
	    std::fabs(settings.cutting_origin[1]->value() - 8.0) < 1.0e-12 &&
	    std::fabs(settings.cutting_origin[2]->value()) < 1.0e-12,
	    "settings reflects command-owned cutting-plane updates");
	settings.grid_ckbx->setCheckState(Qt::Checked);
	CHECK(ged_exec(gedp, 2, cutting_get) == BRLCAD_OK &&
	    std::sscanf(bu_vls_cstr(gedp->ged_result_str),
	    "enable %d\norigin %lf %lf %lf\nnormal %lf %lf %lf",
	    &cuttingEnabled, &cuttingOrigin[0], &cuttingOrigin[1],
	    &cuttingOrigin[2], &cuttingNormal[0], &cuttingNormal[1],
	    &cuttingNormal[2]) == 7 && cuttingEnabled == 1 &&
	    std::fabs(cuttingOrigin[1] - 8.0) < 1.0e-12,
	    "unrelated faceplate settings do not overwrite the cutting plane");
	settings.cutting_normal[1]->setValue(0.0);
	CHECK(ged_exec(gedp, 2, cutting_get) == BRLCAD_OK &&
	    std::sscanf(bu_vls_cstr(gedp->ged_result_str),
	    "enable %d\norigin %lf %lf %lf\nnormal %lf %lf %lf",
	    &cuttingEnabled, &cuttingOrigin[0], &cuttingOrigin[1],
	    &cuttingOrigin[2], &cuttingNormal[0], &cuttingNormal[1],
	    &cuttingNormal[2]) == 7 && cuttingEnabled == 1 &&
	    std::fabs(cuttingNormal[1] - 1.0) < 1.0e-12 &&
	    std::fabs(settings.cutting_normal[1]->value() - 1.0) < 1.0e-12,
	    "invalid cutting normal restores authoritative plane state");

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
