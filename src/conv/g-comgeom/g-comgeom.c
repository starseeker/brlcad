/*                     G - C O M G E O M . C
 * BRL-CAD
 *
 * Export BRL-CAD database geometry to COMGEOM format (versions 0,1,4,5).
 * Successor to vdeck: non‑interactive, testable, scriptable.
 *
 * References (succinct):
 *  - Primitive classification (tgc family, ellipsoid reduction): Requicha & Voelcker 1982.
 *  - Boolean tree serialization (future): Foley & Van Dam 1996 (scene graph).
 *  - Data exchange design: Hoffmann 1989.
 */

#include "common.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "bu/app.h"
#include "bu/getopt.h"
#include "bu/units.h"
#include "bu/log.h"
#include "bu/vls.h"
#include "vmath.h"
#include "raytrace.h"
#include "rt/geom.h"
#include "wdb.h"

/* Configuration structure */
struct config {
    int version;              /* COMGEOM version: 0, 1, 4, or 5 */
    char output_prefix[256];  /* Output file prefix */
    int solid_num;            /* Starting solid number */
    int region_num;           /* Starting region number */
    int verbose;              /* Verbosity level */
    int json_report;          /* Generate JSON report */
    char json_path[256];      /* JSON report output path */
    struct bu_vls filter;     /* Object filter pattern (future) */
};

/* Statistics for JSON report */
struct export_stats {
    size_t solid_count;
    size_t region_count;
    size_t unsupported_count;
    size_t error_count;
    size_t warning_count;

    /* Primitive breakdown */
    size_t tor_count;
    size_t tgc_count;
    size_t ell_count;
    size_t arb_count;
    size_t arbn_count;
    size_t ars_count;
    size_t half_count;
    size_t bot_count;
    size_t brep_count;
    size_t other_count;
};

/* File handles for output */
static FILE *solid_fp = NULL;
static FILE *region_fp = NULL;
static FILE *ident_fp  = NULL;

/* Forward declarations */
static int export_solid(struct rt_db_internal *ip, const char *name, struct config *cfg, struct export_stats *stats);
static int export_region(struct directory *dp, struct db_i *dbip, struct config *cfg, struct export_stats *stats);
static void write_json_report(struct config *cfg, struct export_stats *stats);

/* ---- Utility printing helpers ----------------------------------------- */
static void print_vec(FILE *fp, const vect_t v) {
    fprintf(fp, "%10.4f%10.4f%10.4f", v[X], v[Y], v[Z]);
}

static void print_point(FILE *fp, const point_t p) {
    fprintf(fp, "%10.4f%10.4f%10.4f", p[X], p[Y], p[Z]);
}

/* ---- Usage ------------------------------------------------------------- */
static void
usage(const char *progname)
{
    bu_log("Usage: %s [options] input.g [objects...]\n", progname);
    bu_log("Options:\n");
    bu_log("  -V version    COMGEOM version (0, 1, 4, or 5) [default: 5]\n");
    bu_log("  -o prefix     Output file prefix [default: output]\n");
    bu_log("  -s num        Starting solid number [default: 1]\n");
    bu_log("  -r num        Starting region number [default: 1]\n");
    bu_log("  -j path       Generate JSON report at path\n");
    bu_log("  -v            Verbose output\n");
    bu_log("  -h            This help message\n");
    bu_log("\nSupported primitives: TOR, TGC (RCC/TRC/REC/TEC/TGC), ELL (SPH/ELL1/ELLG), ARB8, ARBN, ARS, HALF\n");
}

/* ---- Argument parsing -------------------------------------------------- */
static int
parse_args(int argc, char **argv, struct config *cfg)
{
    int c;

    cfg->version = 5;
    bu_strlcpy(cfg->output_prefix, "output", sizeof(cfg->output_prefix));
    cfg->solid_num = 1;
    cfg->region_num = 1;
    cfg->verbose = 0;
    cfg->json_report = 0;
    cfg->json_path[0] = '\0';
    bu_vls_init(&cfg->filter);

    while ((c = bu_getopt(argc, argv, "V:o:s:r:j:vh?")) != -1) {
        switch (c) {
        case 'V':
            cfg->version = atoi(bu_optarg);
            if (cfg->version != 0 && cfg->version != 1 && cfg->version != 4 && cfg->version != 5) {
                bu_log("ERROR: Invalid version %d (must be 0, 1, 4, or 5)\n", cfg->version);
                return 0;
            }
            break;
        case 'o':
            bu_strlcpy(cfg->output_prefix, bu_optarg, sizeof(cfg->output_prefix));
            break;
        case 's':
            cfg->solid_num = atoi(bu_optarg);
            if (cfg->solid_num < 1) {
                bu_log("ERROR: Invalid starting solid number\n");
                return 0;
            }
            break;
        case 'r':
            cfg->region_num = atoi(bu_optarg);
            if (cfg->region_num < 1) {
                bu_log("ERROR: Invalid starting region number\n");
                return 0;
            }
            break;
        case 'j':
            cfg->json_report = 1;
            bu_strlcpy(cfg->json_path, bu_optarg, sizeof(cfg->json_path));
            break;
        case 'v':
            cfg->verbose++;
            break;
        case 'h':
        case '?':
        default:
            usage(argv[0]);
            return 0;
        }
    }

    return 1;
}

/* ---- Output file management ------------------------------------------- */
static int
open_output_files(struct config *cfg)
{
    char fname[512];

    snprintf(fname, sizeof(fname), "%s.sol", cfg->output_prefix);
    solid_fp = fopen(fname, "w");
    if (!solid_fp) { bu_log("ERROR: Cannot open %s\n", fname); return 0; }

    snprintf(fname, sizeof(fname), "%s.reg", cfg->output_prefix);
    region_fp = fopen(fname, "w");
    if (!region_fp) {
        bu_log("ERROR: Cannot open %s\n", fname);
        fclose(solid_fp); solid_fp = NULL;
        return 0;
    }

    snprintf(fname, sizeof(fname), "%s.id", cfg->output_prefix);
    ident_fp = fopen(fname, "w");
    if (!ident_fp) {
        bu_log("ERROR: Cannot open %s\n", fname);
        fclose(solid_fp); fclose(region_fp);
        solid_fp = region_fp = NULL;
        return 0;
    }

    return 1;
}

static void
close_output_files(void)
{
    if (solid_fp) { fclose(solid_fp); solid_fp = NULL; }
    if (region_fp){ fclose(region_fp); region_fp = NULL; }
    if (ident_fp) { fclose(ident_fp); ident_fp = NULL; }
}

/* ---- Primitive exporters ---------------------------------------------- */

static int
export_tor(struct rt_tor_internal *tor, const char *name, struct config *cfg, struct export_stats *stats)
{
    (void)name;
    if (cfg->verbose) bu_log("Export TOR %s\n", name);
    /* Card 1: center + axis */
    fprintf(solid_fp, "%5d TOR  ", cfg->solid_num);
    print_point(solid_fp, tor->v);
    print_vec(solid_fp, tor->h);
    fprintf(solid_fp, "\n");
    /* Card 2: radii */
    fprintf(solid_fp, "%5d      %10.4f%10.4f\n", cfg->solid_num, tor->r_a, tor->r_h);
    stats->tor_count++; stats->solid_count++; cfg->solid_num++;
    return 0;
}

/* Classify TGC into RCC, TRC, REC, TEC, or remain TGC */
static const char *
classify_tgc(const struct rt_tgc_internal *tgc)
{
    double a = MAGNITUDE(tgc->a);
    double b = MAGNITUDE(tgc->b);
    double c = MAGNITUDE(tgc->c);
    double d = MAGNITUDE(tgc->d);

    int a_eq_b = NEAR_EQUAL(a, b, 1e-6);
    int c_eq_d = NEAR_EQUAL(c, d, 1e-6);
    int a_eq_c = NEAR_EQUAL(a, c, 1e-6);

    if (a_eq_b && c_eq_d) {
        if (a_eq_c) return "RCC";
        return "TRC";
    }
    if (c_eq_d && !a_eq_b) return "REC";

    if (!c_eq_d && !a_eq_b) {
        if (!NEAR_ZERO(c, 1e-6) && !NEAR_ZERO(a, 1e-6)) {
            double r1 = b/a;
            double r2 = d/c;
            if (NEAR_EQUAL(r1, r2, 1e-6)) return "TEC";
        }
    }
    return "TGC";
}

static int
export_tgc(struct rt_tgc_internal *tgc, const char *name, struct config *cfg, struct export_stats *stats)
{
    (void)name;
    if (cfg->verbose) bu_log("Export TGC %s\n", name);

    const char *type = classify_tgc(tgc);

    fprintf(solid_fp, "%5d %-5s", cfg->solid_num, type);
    print_point(solid_fp, tgc->v);
    print_vec(solid_fp, tgc->h);
    fprintf(solid_fp, "\n");

    if (!strcmp(type, "RCC")) {
        fprintf(solid_fp, "%5d      %10.4f\n", cfg->solid_num, MAGNITUDE(tgc->a));
    } else if (!strcmp(type, "TRC")) {
        fprintf(solid_fp, "%5d      %10.4f%10.4f\n", cfg->solid_num, MAGNITUDE(tgc->a), MAGNITUDE(tgc->c));
    } else if (!strcmp(type, "REC")) {
        fprintf(solid_fp, "%5d      ", cfg->solid_num);
        print_vec(solid_fp, tgc->a);
        print_vec(solid_fp, tgc->b);
        fprintf(solid_fp, "\n");
    } else if (!strcmp(type, "TEC")) {
        fprintf(solid_fp, "%5d      ", cfg->solid_num);
        print_vec(solid_fp, tgc->a);
        print_vec(solid_fp, tgc->b);
        fprintf(solid_fp, "\n");
        double ratio = (MAGNITUDE(tgc->a) > SMALL_FASTF && MAGNITUDE(tgc->c) > SMALL_FASTF)
                       ? (MAGNITUDE(tgc->c)/MAGNITUDE(tgc->a)) : 0.0;
        fprintf(solid_fp, "%5d      %10.4f\n", cfg->solid_num, ratio);
    } else {
        fprintf(solid_fp, "%5d      ", cfg->solid_num);
        print_vec(solid_fp, tgc->a);
        print_vec(solid_fp, tgc->b);
        fprintf(solid_fp, "\n");
        fprintf(solid_fp, "%5d      ", cfg->solid_num);
        print_vec(solid_fp, tgc->c);
        print_vec(solid_fp, tgc->d);
        fprintf(solid_fp, "\n");
    }

    stats->tgc_count++; stats->solid_count++; cfg->solid_num++;
    return 0;
}

static int
export_ell(struct rt_ell_internal *ell, const char *name, struct config *cfg, struct export_stats *stats)
{
    (void)name;
    if (cfg->verbose) bu_log("Export ELL %s\n", name);

    double a = MAGNITUDE(ell->a);
    double b = MAGNITUDE(ell->b);
    double c = MAGNITUDE(ell->c);

    const char *type = "ELLG";
    if (NEAR_EQUAL(a, b, 1e-6) && NEAR_EQUAL(b, c, 1e-6)) type = "SPH";
    else if (NEAR_EQUAL(a, b, 1e-6) || NEAR_EQUAL(b, c, 1e-6) || NEAR_EQUAL(a, c, 1e-6))
        type = "ELL1";

    fprintf(solid_fp, "%5d %-5s", cfg->solid_num, type);
    print_point(solid_fp, ell->v);
    print_vec(solid_fp, ell->a);
    fprintf(solid_fp, "\n");

    if (strcmp(type, "SPH") != 0) {
        fprintf(solid_fp, "%5d      ", cfg->solid_num);
        print_vec(solid_fp, ell->b);
        print_vec(solid_fp, ell->c);
        fprintf(solid_fp, "\n");
    }

    stats->ell_count++; stats->solid_count++; cfg->solid_num++;
    return 0;
}

static int
export_arb(struct rt_arb_internal *arb, const char *name, struct config *cfg, struct export_stats *stats)
{
    (void)name;
    if (cfg->verbose) bu_log("Export ARB8 %s\n", name);

    fprintf(solid_fp, "%5d ARB8 ", cfg->solid_num);
    print_point(solid_fp, arb->pt[0]); print_point(solid_fp, arb->pt[1]); fprintf(solid_fp, "\n");
    fprintf(solid_fp, "%5d      ", cfg->solid_num);
    print_point(solid_fp, arb->pt[2]); print_point(solid_fp, arb->pt[3]); fprintf(solid_fp, "\n");
    fprintf(solid_fp, "%5d      ", cfg->solid_num);
    print_point(solid_fp, arb->pt[4]); print_point(solid_fp, arb->pt[5]); fprintf(solid_fp, "\n");
    fprintf(solid_fp, "%5d      ", cfg->solid_num);
    print_point(solid_fp, arb->pt[6]); print_point(solid_fp, arb->pt[7]); fprintf(solid_fp, "\n");

    stats->arb_count++; stats->solid_count++; cfg->solid_num++;
    return 0;
}

static int
export_arbn(struct rt_arbn_internal *arbn, const char *name, struct config *cfg, struct export_stats *stats)
{
    (void)name;
    if (cfg->verbose) bu_log("Export ARBN %s (%zu planes)\n", name, arbn->neqn);

    fprintf(solid_fp, "%5d ARBN %5zu\n", cfg->solid_num, arbn->neqn);
    for (size_t i = 0; i < arbn->neqn; ++i) {
        fprintf(solid_fp, "%5d      %10.4f%10.4f%10.4f%10.4f\n",
                cfg->solid_num,
                arbn->eqn[i][X], arbn->eqn[i][Y],
                arbn->eqn[i][Z], arbn->eqn[i][W]);
    }

    stats->arbn_count++; stats->solid_count++; cfg->solid_num++;
    return 0;
}

static int
export_ars(struct rt_ars_internal *ars, const char *name, struct config *cfg, struct export_stats *stats)
{
    (void)name;
    if (cfg->verbose) bu_log("Export ARS %s (%zu curves)\n", name, ars->ncurves);

    fprintf(solid_fp, "%5d ARS  %5zu%5zu\n", cfg->solid_num, ars->ncurves, ars->pts_per_curve);

    for (size_t i = 0; i < ars->ncurves; ++i) {
        fastf_t *curve = ars->curves[i];
        for (size_t j = 0; j < ars->pts_per_curve; ++j) {
            const fastf_t *p = &curve[3*j];
            fprintf(solid_fp, "%5d      %10.4f%10.4f%10.4f\n",
                    cfg->solid_num, p[0], p[1], p[2]);
        }
    }

    stats->ars_count++; stats->solid_count++; cfg->solid_num++;
    return 0;
}

static int
export_half(struct rt_half_internal *half, const char *name, struct config *cfg, struct export_stats *stats)
{
    (void)name;
    if (cfg->verbose) bu_log("Export HALF %s\n", name);

    fprintf(solid_fp, "%5d HALF ", cfg->solid_num);
    fprintf(solid_fp, "%10.4f%10.4f%10.4f%10.4f\n",
            half->eqn[X], half->eqn[Y], half->eqn[Z], half->eqn[W]);

    stats->half_count++; stats->solid_count++; cfg->solid_num++;
    return 0;
}

/* ---- Dispatch for primitive export ------------------------------------ */
static int
export_solid(struct rt_db_internal *ip, const char *name, struct config *cfg, struct export_stats *stats)
{
    int ret = 0;

    switch (ip->idb_type) {
    case ID_TOR:
        ret = export_tor((struct rt_tor_internal *)ip->idb_ptr, name, cfg, stats);
        break;
    case ID_TGC:
        ret = export_tgc((struct rt_tgc_internal *)ip->idb_ptr, name, cfg, stats);
        break;
    case ID_ELL:
        ret = export_ell((struct rt_ell_internal *)ip->idb_ptr, name, cfg, stats);
        break;
    case ID_ARB8:
        ret = export_arb((struct rt_arb_internal *)ip->idb_ptr, name, cfg, stats);
        break;
    case ID_ARBN:
        ret = export_arbn((struct rt_arbn_internal *)ip->idb_ptr, name, cfg, stats);
        break;
    case ID_ARS:
        ret = export_ars((struct rt_ars_internal *)ip->idb_ptr, name, cfg, stats);
        break;
    case ID_HALF:
        ret = export_half((struct rt_half_internal *)ip->idb_ptr, name, cfg, stats);
        break;
    case ID_BOT:
        bu_log("WARNING: BOT '%s' not supported, skipping\n", name);
        stats->bot_count++; stats->unsupported_count++; stats->warning_count++; ret = -1;
        break;
    case ID_BREP:
        bu_log("WARNING: BREP '%s' not supported, skipping\n", name);
        stats->brep_count++; stats->unsupported_count++; stats->warning_count++; ret = -1;
        break;
    default:
        bu_log("WARNING: Unsupported primitive id=%d '%s'\n", ip->idb_type, name);
        stats->other_count++; stats->unsupported_count++; stats->warning_count++; ret = -1;
        break;
    }

    if (ret != 0) stats->error_count += (ret < 0);
    else stats->solid_count++;
    return ret;
}

/* ---- Region export (placeholder) -------------------------------------- */
static int
export_region(struct directory *dp, struct db_i *dbip, struct config *cfg, struct export_stats *stats)
{
    (void)dbip;
    if (cfg->verbose) bu_log("Export region (stub): %s\n", dp->d_namep);

    fprintf(region_fp, "%5d  + %5d\n", cfg->region_num, cfg->solid_num - 1);

    fprintf(ident_fp, "%5d%5d%5d%5d%5d\n",
            cfg->region_num, cfg->region_num, 0, 1, 100);

    stats->region_count++; cfg->region_num++;
    return 0;
}

/* ---- JSON report ------------------------------------------------------- */
static void
write_json_report(struct config *cfg, struct export_stats *stats)
{
    FILE *fp = fopen(cfg->json_path, "w");
    if (!fp) {
        bu_log("ERROR: Cannot open JSON report file %s\n", cfg->json_path);
        return;
    }

    fprintf(fp, "{\n");
    fprintf(fp, "  \"version\": %d,\n", cfg->version);
    fprintf(fp, "  \"output_prefix\": \"%s\",\n", cfg->output_prefix);
    fprintf(fp, "  \"statistics\": {\n");
    fprintf(fp, "    \"solids\": %zu,\n", stats->solid_count);
    fprintf(fp, "    \"regions\": %zu,\n", stats->region_count);
    fprintf(fp, "    \"unsupported\": %zu,\n", stats->unsupported_count);
    fprintf(fp, "    \"errors\": %zu,\n", stats->error_count);
    fprintf(fp, "    \"warnings\": %zu\n", stats->warning_count);
    fprintf(fp, "  },\n");
    fprintf(fp, "  \"primitives\": {\n");
    fprintf(fp, "    \"TOR\": %zu,\n", stats->tor_count);
    fprintf(fp, "    \"TGC\": %zu,\n", stats->tgc_count);
    fprintf(fp, "    \"ELL\": %zu,\n", stats->ell_count);
    fprintf(fp, "    \"ARB\": %zu,\n", stats->arb_count);
    fprintf(fp, "    \"ARBN\": %zu,\n", stats->arbn_count);
    fprintf(fp, "    \"ARS\": %zu,\n", stats->ars_count);
    fprintf(fp, "    \"HALF\": %zu,\n", stats->half_count);
    fprintf(fp, "    \"BOT\": %zu,\n", stats->bot_count);
    fprintf(fp, "    \"BREP\": %zu,\n", stats->brep_count);
    fprintf(fp, "    \"Other\": %zu\n", stats->other_count);
    fprintf(fp, "  }\n");
    fprintf(fp, "}\n");

    fclose(fp);
    if (cfg->verbose) bu_log("JSON report written: %s\n", cfg->json_path);
}

/* ---- Main -------------------------------------------------------------- */
int
main(int argc, char **argv)
{
    struct config cfg;
    struct export_stats stats;
    struct db_i *dbip;
    struct directory *dp;
    struct rt_db_internal intern;
    int i;

    bu_setprogname(argv[0]);
    memset(&stats, 0, sizeof(stats));

    if (!parse_args(argc, argv, &cfg))
        return 1;

    if (bu_optind >= argc) {
        bu_log("ERROR: No input database specified\n");
        usage(argv[0]);
        return 1;
    }

    const char *db_name = argv[bu_optind++];
    dbip = db_open(db_name, DB_OPEN_READONLY);
    if (dbip == DBI_NULL) {
        bu_log("ERROR: Cannot open database %s\n", db_name);
        return 1;
    }

    if (db_dirbuild(dbip) < 0) {
        bu_log("ERROR: Cannot build directory for %s\n", db_name);
        db_close(dbip);
        return 1;
    }

    if (!open_output_files(&cfg)) {
        db_close(dbip);
        return 1;
    }

    fprintf(solid_fp, "BRL-CAD COMGEOM v%d export\n", cfg.version);
    fprintf(solid_fp, "%5d%5d\n", 0, 0);

    if (bu_optind < argc) {
        for (i = bu_optind; i < argc; ++i) {
            dp = db_lookup(dbip, argv[i], LOOKUP_QUIET);
            if (dp == RT_DIR_NULL) {
                bu_log("WARNING: Object '%s' not found\n", argv[i]);
                stats.warning_count++;
                continue;
            }
            if (dp->d_flags & RT_DIR_SOLID) {
                if (rt_db_get_internal(&intern, dp, dbip, NULL, &rt_uniresource) >= 0) {
                    export_solid(&intern, dp->d_namep, &cfg, &stats);
                    rt_db_free_internal(&intern);
                }
            } else if (dp->d_flags & RT_DIR_REGION) {
                export_region(dp, dbip, &cfg, &stats);
            }
        }
    } else {
        for (i = 0; i < RT_DBNHASH; i++) {
            for (dp = dbip->dbi_Head[i]; dp != RT_DIR_NULL; dp = dp->d_forw) {
                if (dp->d_flags & RT_DIR_SOLID) {
                    if (rt_db_get_internal(&intern, dp, dbip, NULL, &rt_uniresource) >= 0) {
                        export_solid(&intern, dp->d_namep, &cfg, &stats);
                        rt_db_free_internal(&intern);
                    }
                }
            }
        }
    }

    rewind(solid_fp);
    fprintf(solid_fp, "BRL-CAD COMGEOM v%d export\n", cfg.version);
    fprintf(solid_fp, "%5zu%5zu\n", stats.solid_count, stats.region_count);

    close_output_files();
    db_close(dbip);
    bu_vls_free(&cfg.filter);

    if (cfg.json_report)
        write_json_report(&cfg, &stats);

    if (cfg.verbose) {
        bu_log("Summary:\n");
        bu_log("  Solids: %zu\n", stats.solid_count);
        bu_log("  Regions: %zu\n", stats.region_count);
        bu_log("  Unsupported: %zu\n", stats.unsupported_count);
        bu_log("  Errors: %zu\n", stats.error_count);
    }

    return (stats.error_count > 0) ? 1 : 0;
}