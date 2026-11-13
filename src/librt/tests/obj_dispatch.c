/*                    O B J _ D I S P A T C H . C
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/* Verify generic geometry dispatch honors the internal's runtime type. */

#include "common.h"

#include "bu/app.h"
#include "raytrace.h"

static int
wireframe_providers_registered(void)
{
    static const int provider_ids[] = {
	ID_TOR,
	ID_ELL,
	ID_SPH,
	ID_BSPLINE,
	ID_EBM,
	ID_VOL,
	ID_PIPE,
	ID_RPC,
	ID_RHC,
	ID_EPA,
	ID_EHY,
	ID_ETO,
	ID_HF,
	ID_DSP,
	ID_SKETCH,
	ID_EXTRUDE,
	ID_SUPERELL,
	ID_BREP,
	ID_REVOLVE,
	ID_HRT,
	ID_DATUM
    };

    for (size_t i = 0; i < sizeof(provider_ids) / sizeof(provider_ids[0]); i++) {
	const struct rt_functab *ft = &OBJ[provider_ids[i]];

	if (!ft->ft_wireframe_line_set) {
	    bu_log("missing canonical wireframe provider for %s\n", ft->ft_name);
	    return 0;
	}
    }

    return 1;
}

int
main(int UNUSED(argc), const char **argv)
{
    struct rt_db_internal intern;
    struct nmgregion *region = NULL;
    struct bu_list vhead;

    bu_setprogname(argv[0]);

    if (!wireframe_providers_registered())
	return 1;
    if (!OBJ[ID_VOL].ft_indexed_face_set) {
	bu_log("missing canonical indexed-face provider for vol\n");
	return 1;
    }

    RT_DB_INTERNAL_INIT(&intern);
    intern.idb_major_type = DB5_MAJORTYPE_BINARY_UNIF;
    intern.idb_minor_type = DB5_MINORTYPE_BINU_8BITINT_U;
    intern.idb_meth = &OBJ[ID_BINUNIF];
    BU_LIST_INIT(&vhead);

    /* This binunif subtype's numeric value overlaps a geometry primitive.
     * Both generic operations must reject the unsupported runtime method,
     * not reinterpret the payload according to that numeric coincidence. */
    if (rt_obj_tess(&region, NULL, &intern, NULL, NULL) != -4 || region ||
	rt_obj_plot(&vhead, &intern, NULL, NULL) != -4 ||
	!BU_LIST_IS_EMPTY(&vhead))
	return 2;

    bu_avs_free(&intern.idb_avs);
    return 0;
}
