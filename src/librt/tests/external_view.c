/*                    E X T E R N A L _ V I E W . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include <string.h>

#include "bu/app.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "raytrace.h"
#include "wdb.h"


int
main(int argc, char **argv)
{
    char path[MAXPATHLEN] = {0};
    struct rt_wdb *wdbp;
    struct db_i *dbip;
    struct directory *dp;
    struct bu_external copied;
    const unsigned char *view;
    size_t view_bytes = 0;
    point_t center = VINIT_ZERO;
    FILE *fp;
    int ret = 0;

    (void)argc;
    bu_setprogname(argv[0]);
    fp = bu_temp_file(path, sizeof(path));
    if (!fp)
	return 1;
    fclose(fp);
    bu_file_delete(path);

    wdbp = wdb_fopen_v(path, 5);
    if (!wdbp)
	return 1;
    if (mk_id(wdbp, "borrowed serialized object test") ||
	mk_sph(wdbp, "sample.s", center, 7.0)) {
	db_close(wdbp->dbip);
	bu_file_delete(path);
	return 1;
    }
    db_close(wdbp->dbip);

    dbip = db_open(path, DB_OPEN_READONLY);
    if (!dbip || db_dirbuild(dbip) < 0) {
	if (dbip)
	    db_close(dbip);
	bu_file_delete(path);
	return 1;
    }
    dp = db_lookup(dbip, "sample.s", LOOKUP_QUIET);
    BU_EXTERNAL_INIT(&copied);
    view = db_external_view(dbip, dp, &view_bytes);
    if (!dp || !view || view_bytes != dp->d_len ||
	db_get_external(&copied, dp, dbip) < 0 ||
	copied.ext_nbytes != view_bytes ||
	memcmp(view, copied.ext_buf, view_bytes) != 0) {
	bu_log("db_external_view did not match db_get_external\n");
	ret = 1;
    }
    bu_free_external(&copied);

    view_bytes = 99;
    if (db_external_view(dbip, NULL, &view_bytes) != NULL || view_bytes != 0) {
	bu_log("db_external_view did not reject a NULL directory\n");
	ret = 1;
    }

    db_close(dbip);
    bu_file_delete(path);
    return ret;
}
