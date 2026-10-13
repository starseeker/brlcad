/* BRL-CAD
 *
 * Copyright (c) 1986-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */

#include "common.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "raytrace.h"
#include "rt/edit.h"
#include "rt/geom.h"
#include "rt/misc.h"
#include "rt/primitives/arb8.h"

static int
facedef_error(struct rt_edit *edit, const char *message)
{
    bu_vls_printf(edit->log_str, "facedef: %s\n", message);
    return BRLCAD_ERROR;
}

int
rt_arb_facedef(struct rt_edit *edit, int argc, const char *argv[])
{
    if (!edit || !edit->log_str)
        return BRLCAD_ERROR;
    if (edit->es_int.idb_type != ID_ARB8)
        return facedef_error(edit, "solid type must be ARB");
    if (argc < 3)
        return facedef_error(edit, "expected face and definition form");
    int type = rt_arb_std_type(&edit->es_int, edit->tol);
    int product = 1;
    size_t digits = strlen(argv[1]);
    if (digits < 3 || digits > 4)
        return facedef_error(edit, "expected three or four face vertices");
    for (size_t i = 0; i < digits; i++) {
        if (argv[1][i] < '1' || argv[1][i] > '8')
            return facedef_error(edit, "invalid face vertex");
        product *= argv[1][i] - '0';
    }
    int plane;
    /* These products identify the conventional compact ARB face labels. */
    switch (product) {
        case 6: case 24: plane = type == ARB4 && product == 24 ? 2 : 0; break;
        case 8: case 180: case 210: case 1680: plane = 1; break;
        case 30: case 120: case 20: case 160: plane = type == ARB5 ? 4 : 2; break;
        case 12: case 10: case 252: plane = type == ARB5 ? 1 : 3; break;
        case 72: case 60: plane = type == ARB5 ? 3 : 4; break;
        case 420: case 672: plane = 5; break;
        default: return facedef_error(edit, "unknown face");
    }
    const char form = argv[2][0];
    if (!argv[2][0] || argv[2][1])
        return facedef_error(edit, "invalid definition form");
    if (form == 'q')
        return BRLCAD_OK;
    const bool constrained = type == ARB7 && plane != 0 && plane != 3;
    if (constrained && form != 'c')
        return facedef_error(edit, "cannot redefine that ARB7 plane");
    int count = form == 'a' ? 4 : form == 'b' ? 9 : form == 'd' ? 3 :
        form == 'c' ? (constrained ? 2 : argc > 5 && argv[5][0] == 'v' ? 3 : 5) : -1;
    if (count < 0) {
        bu_vls_printf(edit->log_str, "facedef: %s is not an option\n", argv[2]);
        return BRLCAD_ERROR;
    }
    if (argc != count + 3)
        return facedef_error(edit, "incorrect number of definition parameters");
    fastf_t values[9] = {0};
    for (int i = 0; i < count; i++) {
        if (form == 'c' && i == 2 && argv[5][0] == 'v')
            continue;
        char *end = NULL;
        values[i] = strtod(argv[i + 3], &end);
        if (!end || end == argv[i + 3] || *end || !std::isfinite(values[i]))
            return facedef_error(edit, "expected finite numeric definition parameters");
    }
    struct rt_db_internal transformed;
    RT_DB_INTERNAL_INIT(&transformed);
    if (rt_generic_xform(&transformed, edit->e_mat, &edit->es_int, 0, edit->dbip) < 0)
        return facedef_error(edit, "unable to transform working geometry");
    struct rt_arb_internal *arb = static_cast<struct rt_arb_internal *>(transformed.idb_ptr);
    plane_t planes[6];
    int ret = rt_arb_calc_planes(edit->log_str, arb, type, planes, edit->tol);
    if (ret) {
        rt_db_free_internal(&transformed);
        return BRLCAD_ERROR;
    }
    fastf_t *definition = planes[plane];
    point_t point;
    if (form == 'a') {
        VMOVE(definition, values);
        if (MAGNITUDE(definition) <= SMALL_FASTF) {
            rt_db_free_internal(&transformed);
            return facedef_error(edit, "plane normal must be nonzero");
        }
        VUNITIZE(definition);
        definition[W] = values[3] * edit->local2base;
    } else if (form == 'b') {
        point_t points[3];
        for (int i = 0; i < 3; i++)
            VSCALE(points[i], values + i * 3, edit->local2base);
        ret = bg_make_plane_3pnts(definition, points[0], points[1], points[2], edit->tol);
    } else if (form == 'd') {
        VSCALE(point, values, edit->local2base);
        definition[W] = VDOT(definition, point);
    } else {
        fastf_t rotation = values[0] * DEG2RAD;
        fastf_t fallback = values[1] * DEG2RAD;
        VSET(definition, cos(fallback) * cos(rotation),
            cos(fallback) * sin(rotation), sin(fallback));
        if (constrained || argv[5][0] == 'v') {
            char *end = NULL;
            long vertex = constrained ? 5 : strtol(argv[5] + 1, &end, 10);
            if (!constrained && (!end || end == argv[5] + 1 || *end || vertex < 1 || vertex > 8)) {
                rt_db_free_internal(&transformed);
                return facedef_error(edit, "invalid fixed vertex");
            }
            VMOVE(point, arb->pt[vertex - 1]);
        } else {
            VSCALE(point, values + 2, edit->local2base);
        }
        definition[W] = VDOT(definition, point);
    }
    if (!ret)
        ret = rt_arb_calc_points(arb, type, (const plane_t *)planes, edit->tol);
    if (ret) {
        rt_db_free_internal(&transformed);
        return facedef_error(edit, "unable to construct the redefined solid");
    }
    struct rt_arb_internal *original = static_cast<struct rt_arb_internal *>(edit->es_int.idb_ptr);
    for (int i = 0; i < 8; i++)
        MAT4X3PNT(original->pt[i], edit->e_invmat, arb->pt[i]);
    rt_db_free_internal(&transformed);
    return BRLCAD_OK;
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
