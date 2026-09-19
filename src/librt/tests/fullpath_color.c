/*                   F U L L P A T H _ C O L O R . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"

#include "bu/color.h"
#include "raytrace.h"
#include "wdb.h"


static int
add_sphere(struct rt_wdb *wdbp, const char *name)
{
    point_t center = VINIT_ZERO;

    return mk_sph(wdbp, name, center, 1.0);
}


static int
add_comb(struct rt_wdb *wdbp,
	const char *name,
	const char *member,
	int region,
	const unsigned char *rgb,
	int region_id,
	int inherit)
{
    struct wmember wm;

    BU_LIST_INIT(&wm.l);
    if (!mk_addmember(member, &wm.l, NULL, WMOP_UNION))
	return 1;

    return mk_comb(wdbp, name, &wm.l, region, NULL, NULL, rgb,
	    region_id, 0, 0, 0, inherit, 0, 0);
}


static int
path_color_rgb(struct db_i *dbip, const char *path, unsigned char rgb[3])
{
    struct db_full_path fp;
    struct bu_color color = BU_COLOR_INIT_ZERO;
    int ret;

    db_full_path_init(&fp);
    if (db_string_to_path(&fp, dbip, path) != 0)
	return 0;

    db_full_path_color(&color, &fp, dbip);
    ret = bu_color_to_rgb_chars(&color, rgb);
    db_free_full_path(&fp);
    return ret;
}


static int
expect_rgb(struct db_i *dbip,
	const char *label,
	const char *path,
	unsigned char r,
	unsigned char g,
	unsigned char b)
{
    unsigned char rgb[3] = {0, 0, 0};

    if (!path_color_rgb(dbip, path, rgb)) {
	bu_log("%s: failed to resolve path color for %s\n", label, path);
	return 1;
    }

    if (rgb[0] != r || rgb[1] != g || rgb[2] != b) {
	bu_log("%s: expected %s color %u/%u/%u, got %u/%u/%u\n",
		label, path, r, g, b, rgb[0], rgb[1], rgb[2]);
	return 1;
    }

    return 0;
}


int
main(int UNUSED(argc), const char **UNUSED(argv))
{
    bu_setprogname("fullpath_color");
    struct db_i *dbip = db_create_inmem();
    struct rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_INMEM);
    const unsigned char red[3] = {200, 0, 0};
    const unsigned char blue[3] = {0, 0, 200};
    int failures = 0;

    if (!dbip || !wdbp)
	return 1;

    failures += add_sphere(wdbp, "leaf.s");
    failures += db5_update_attribute("leaf.s", "color", "250/250/0", dbip);
    db_mater_add(dbip, 77, 77, 10, 20, 255, MATER_NO_ADDR);

    failures += add_comb(wdbp, "region_table.r", "leaf.s", 1, NULL, 77, 0);
    failures += add_comb(wdbp, "region_blue.r", "leaf.s", 1, blue, 77, 0);
    failures += add_comb(wdbp, "parent_inherit.c", "region_blue.r", 0, red, 0, 1);
    failures += add_comb(wdbp, "parent_lower.c", "region_blue.r", 0, red, 0, 0);
    failures += add_comb(wdbp, "parent_table.c", "region_table.r", 0, NULL, 0, 0);

    if (failures) {
	db_close(dbip);
	return 1;
    }

    failures += expect_rgb(dbip, "inherited parent explicit color",
	    "parent_inherit.c/region_blue.r/leaf.s", 200, 0, 0);
    failures += expect_rgb(dbip, "lower region explicit color",
	    "parent_lower.c/region_blue.r/leaf.s", 0, 0, 200);
    failures += expect_rgb(dbip, "region-id table fallback",
	    "parent_table.c/region_table.r/leaf.s", 10, 20, 255);

    db_close(dbip);
    return failures ? 1 : 0;
}
