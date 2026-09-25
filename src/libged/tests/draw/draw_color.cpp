/*                      D R A W _ C O L O R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 */

#include "common.h"

#include <algorithm>
#include <array>
#include <cstdio>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "bu/app.h"
#include "bu/file.h"
#include "ged.h"
#include "rt/geom.h"
#include "wdb.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BSceneController.h"
#include "BObol/BViewController.h"
#include "ged/display_obol_private.h"
#include "ged/view.h"
#include "view_test_util.h"


using Color = std::array<unsigned char, 3>;
using ExpectedColors = std::map<std::string, Color>;

static bool
write_group(struct rt_wdb *wdbp, const char *name,
    const std::vector<std::string> &children, const unsigned char *color, int inherit)
{
    struct wmember members;
    BU_LIST_INIT(&members.l);
    for (const auto &child : children) {
	if (!mk_addmember(child.c_str(), &members.l, NULL, WMOP_UNION)) {
	    mk_freemembers(&members.l);
	    return false;
	}
    }
    return mk_comb(wdbp, name, &members.l, 0, NULL, NULL, color,
	0, 0, 0, 0, inherit, 0, 0) == 0;
}

static bool
create_fixture(const char *path, ExpectedColors &expected)
{
    std::unique_ptr<struct rt_wdb, decltype(&wdb_close)> database(wdb_fopen(path), wdb_close);
    if (!database)
	return false;

    ON_PlaneSurface surface(ON_xy_plane);
    surface.SetExtents(0, ON_Interval(0.0, 1.0));
    surface.SetExtents(1, ON_Interval(0.0, 1.0));
    ON_Brep brep;
    if (!brep.NewFace(surface) || !brep.IsValid())
	return false;
    fastf_t vertices[] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    int faces[] = {0, 1, 2};

    const Color parent_color = {{20, 40, 60}};
    const struct {
	const char *name;
	const char *attribute;
	bool colored_parent;
	int inherit;
	Color expected;
	const char *rgb_alias = NULL;
    } cases[] = {
	{"default", NULL, false, 0, {{255, 0, 0}}},
	{"primitive", "126/137/141", false, 0, {{126, 137, 141}}},
	{"inherited", NULL, true, 1, parent_color},
	{"primitive_parent", "126/137/141", true, 0, {{126, 137, 141}}},
	{"primitive_inherit", "126/137/141", true, 1, {{126, 137, 141}}},
	{"invalid", "invalid", true, 0, parent_color},
	{"negative", "-1/20/30", true, 0, parent_color},
	{"clamped", "999/10/20", false, 0, {{255, 10, 20}}},
	{"alias", NULL, true, 1, {{80, 120, 160}}, "80/120/160"},
	{"canonical_alias", "126/137/141", true, 1, {{126, 137, 141}}, "80/120/160"}
    };

    std::vector<std::string> groups;
    for (const auto &test : cases) {
	const std::string brep_name = std::string(test.name) + ".brep";
	const std::string bot_name = std::string(test.name) + ".bot";
	if (mk_brep(database.get(), brep_name.c_str(), &brep) < 0 ||
	    mk_bot(database.get(), bot_name.c_str(), RT_BOT_SURFACE, RT_BOT_CCW, 0,
		3, 1, vertices, faces, NULL, NULL) < 0)
	    return false;
	const std::vector<std::string> children = {brep_name, bot_name};
	for (const auto &child : children) {
	    expected.emplace(child, test.expected);
	    if (test.attribute && db5_update_attribute(child.c_str(), db5_standard_attribute(ATTR_COLOR),
		test.attribute, database->dbip) < 0)
		return false;
	    if (test.rgb_alias && db5_update_attribute(child.c_str(), "rgb",
		test.rgb_alias, database->dbip) < 0)
		return false;
	}
	if (!write_group(database.get(), test.name, children,
	    test.colored_parent ? parent_color.data() : NULL, test.inherit))
	    return false;
	groups.emplace_back(test.name);
    }
    return write_group(database.get(), "all", groups, NULL, 0);
}

static bool
check_drawing(struct ged *gedp, const char *path, const char *mode,
    const ExpectedColors &expected, bool override_color)
{
    const char *zap[] = {"zap", NULL};
    if (ged_exec_zap(gedp, 1, zap) != BRLCAD_OK)
	return false;

    const Color override_rgb = {{9, 80, 150}};
    const char *draw[] = {"draw", "-m", mode, path, NULL};
    const char *draw_override[] = {"draw", "-m", mode, "-C", "9/80/150", path, NULL};
    if ((override_color ? ged_exec_draw(gedp, 6, draw_override) : ged_exec_draw(gedp, 4, draw)) != BRLCAD_OK) {
	bu_log("draw mode %s failed: %s\n", mode, bu_vls_cstr(gedp->ged_result_str));
	return false;
    }

    const char *autoview[] = {"autoview", NULL};
    if (ged_exec_autoview(gedp, 1, autoview) != BRLCAD_OK)
	return false;

    ExpectedColors remaining = expected;
    bool passed = true;
    struct ged_view_context *view = ged_view_active_ctx(gedp);
    if (!draw_test_obol_progressive_drain(gedp, view, 2000, 1)) {
	bu_log("draw mode %s did not settle\n", mode);
	return false;
    }
    BObolViewController *controller = static_cast<BObolViewController *>(
	bobol_display_endpoint_controller(ged_view_context_obol_endpoint_get(view)));
    if (!controller || !controller->getRenderSceneRoot())
	return false;
    BObolSceneController render_scene(controller->getRenderSceneRoot());
    BObolSceneController *scene = &render_scene;
    for (int i = 0; i < scene->getRealizedShapeSummaryCount(); ++i) {
	BObolRealizedShapeSummary object;
	if (!scene->getRealizedShapeSummary(i, object) || !object.valid)
	    return false;
	if (object.recordRole == "lod-overview")
	    continue;
	const std::string object_path = object.path.getString();
	const size_t separator = object_path.find_last_of('/');
	const std::string name = object_path.substr(separator == std::string::npos ? 0 : separator + 1);
	const auto found = remaining.find(name);
	if (found == remaining.end()) {
	    bu_log("unexpected or duplicate drawn object: %s, path %s, kind %s, role %s\n", name.c_str(), object.path.getString(), object.geometryKind.getString(), object.recordRole.getString());
	    return false;
	}
	const Color &color = override_color ? override_rgb : found->second;
	const SbColor actual = object.colorOverride ? object.color : object.materialColor;
	const SbColor expected_color(color[0] / 255.0f, color[1] / 255.0f, color[2] / 255.0f);
	if ((!object.colorOverride && !object.materialColorValid) ||
	    object.geometryKind.getLength() == 0 ||
	    !actual.equals(expected_color, SMALL_FASTF)) {
	    bu_log("%s, mode %s, override %d: expected %d/%d/%d, got %.0f/%.0f/%.0f\n",
		name.c_str(), mode, override_color, color[0], color[1], color[2],
		actual[0] * 255.0f, actual[1] * 255.0f, actual[2] * 255.0f);
	    passed = false;
	}
	remaining.erase(found);
    }
    if (!remaining.empty()) {
	for (const auto &missing : remaining)
	    bu_log("draw %s, mode %s, override %d: omitted %s\n",
		path, mode, override_color, missing.first.c_str());
	draw_test_obol_debug_scene(gedp, 0, view);
    }
    return passed && remaining.empty();
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    if (argc != 1)
	return 1;
    char path[MAXPATHLEN] = {0};
    FILE *temporary = bu_temp_file(path, sizeof(path));
    if (!temporary)
	return 1;
    if (std::fclose(temporary)) {
	bu_file_delete(path);
	return 1;
    }

    ON::Begin();
    ExpectedColors expected;
    bool passed = create_fixture(path, expected);
    std::unique_ptr<struct ged, decltype(&ged_close)> context(
	passed ? ged_open("db", path, 1) : NULL, ged_close);
    if (context) {
	// Differently colored siblings detect material state leaking between paths.
	if (draw_test_obol_view_init(context.get(), ged_view_active_ctx(context.get()),
	    512, 512) != BRLCAD_OK)
	    return EXIT_FAILURE;
	for (const char *mode : {"0", "1", "2"})
	    for (bool override_color : {false, true})
		passed = check_drawing(context.get(), "all", mode, expected, override_color) && passed;
	for (const char *primitive_path : {"alias.bot", "alias.brep"}) {
	    const ExpectedColors primitive = {{primitive_path, expected.at(primitive_path)}};
	    for (const char *mode : {"0", "1", "2"})
		passed = check_drawing(context.get(), primitive_path, mode, primitive, false) && passed;
	}
    } else {
	bu_log("could not create or open drawing color fixture\n");
	passed = false;
    }
    context.reset();
    ON::End();
    bu_file_delete(path);
    return passed ? 0 : 1;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
