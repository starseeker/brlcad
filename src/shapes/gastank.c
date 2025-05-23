/*                       G A S T A N K . C
 * BRL-CAD
 *
 * Copyright (c) 2004-2025 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file shapes/gastank.c
 *
 * Program to create a gas tank using libwdb.  All dimensions are in
 * mm.  The gas tank is composed of 3 arb8s, 8 spheres, and 12
 * cylinders.  The gas tank is solid and centered at (0, 0, 0).
 * Note that it is not approparte as-is for any application which
 * expects the internal volume of the tank to be modeled as well.
 *
 * Introduced in BRL-CAD release 4.4.
 */

#include "common.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rt/db4.h"
#include "vmath.h"
#include "bu/app.h"
#include "raytrace.h"
#include "wdb.h"

void
explain(void)
{
	fprintf(stderr,"This program constructs a solid gas tank with all\n");
	fprintf(stderr,"edges and corners rounded.  If not used interactively:\n");
	fprintf(stderr,"Usage: gastank [-f mged_file_name] [-n #_of_gastanks] [-H gas_tank_height]\n");
	fprintf(stderr,"       [-w gas_tank_width] [-d gas_tank_depth] [-r radius_of_corners]\n");
	fprintf(stderr,"       (units of mm)\n");
}

int
main(int argc, char **argv)
{
    /* START # 1 */
    struct rt_wdb *fpw;		/* File to be written to. */
#define NAME_LEN 256
    char filemged[NAME_LEN+1] = {0};	/* Mged file create. */
    double hgt=0;       	/* Height, width, & depth of gas tank. */
    double wid=0;
    double dpt=0;
    double rds=0;		/* Radius of the corner of gas tank. */
    point_t pts[8];		/* Points for arb8. */
    point_t bs;			/* Base of cylinder. */
    vect_t ht;			/* Height of cylinder. */
    fastf_t rad;		/* Radius of cylinder & sphere. */
    point_t cent;		/* Center of sphere. */

    /* point_t and vect_t are set using typedef of type fastf_t. */
    /* fastf_t is a type that is machine dependent. */

    char *temp;			/* Temporary character string. */
    char temp1[NAME_LEN+1];		/* Temporary character string. */

    char solnam[9];		/* Solid name. */
    char regnam[9];		/* Region name. */
    char grpnam[5];		/* Group name. */
    int numtnk=0;		/* Number of gas tanks to be created */
				/* (<=maxnumtnk). */
    int maxnumtnk = NAME_LEN;

    struct wmember comb;	/* Used to make regions. */
    struct wmember comb1;	/* Used to make groups. */

    int i, j, k;		/* Loop counters. */
    int ret;

    bu_setprogname(argv[0]);

    /* Set up solid, region, and group names. */
    solnam[0] = 's';
    solnam[1] = '.';
    solnam[2] = 't';
    solnam[3] = 'n';
    solnam[4] = 'k';
    solnam[5] = ' ';
    solnam[6] = '#';
    solnam[7] = '#';
    solnam[8] = '\0';
    regnam[0] = 'r';
    regnam[1] = '.';
    regnam[2] = 't';
    regnam[3] = 'n';
    regnam[4] = 'k';
    regnam[5] = ' ';
    regnam[6] = '#';
    regnam[7] = '#';
    regnam[8] = '\0';
    grpnam[0] = 't';
    grpnam[1] = 'n';
    grpnam[2] = 'k';
    grpnam[3] = ' ';
    grpnam[4] = '\0';

    /* If there are no arguments ask questions. */
    if (argc == 1) {
	/* START # 3 */

	/* Print info about the window. */
	explain();

	/* Find name of mged file to be created. */
	printf("Enter the mged file to be created (%d char max).\n\t", NAME_LEN);
	(void)fflush(stdout);
	ret = scanf(CPP_SCAN(NAME_LEN), filemged);
	if (ret == 0)
	    perror("scanf");
	if (BU_STR_EQUAL(filemged, ""))
	    bu_strlcpy(filemged, "gastank.g", sizeof(filemged));

	/* Find the number of gas tanks to create. */
	printf("Enter the number of gas tanks to create (%d max).\n\t", maxnumtnk);
	(void)fflush(stdout);
	ret = scanf("%d", &numtnk);
	if (ret == 0) {
	    perror("scanf");
	    numtnk = 1;
	}
	else if (numtnk < 1)
	    numtnk = 1;
	else if (numtnk > maxnumtnk)
	    numtnk = maxnumtnk;

	/* Find the dimensions of the gas tanks. */
	printf("Enter the height, width, and depth of the gas tank (units mm).\n\t");
	(void)fflush(stdout);
	ret = scanf("%lf %lf %lf", &hgt, &wid, &dpt);
	if (ret == 0) {
	    perror("scanf");
	    hgt = 1000.0;
	    wid = 1000.0;
	    dpt = 1000.0;
	}
	if (hgt < SMALL_FASTF)
	    hgt = SMALL_FASTF;
	if (wid < SMALL_FASTF)
	    wid = SMALL_FASTF;
	if (dpt < SMALL_FASTF)
	    dpt = SMALL_FASTF;

	printf("Enter the radius of the corners (units mm).\n\t");
	(void)fflush(stdout);
	ret = scanf("%lf", &rds);
	if (ret == 0) {
	    perror("scanf");
	    rds = 10.0;
	}
	if (rds < SMALL_FASTF)
	    rds = SMALL_FASTF;
    }							/* END # 3 */

    /* If there are arguments get answers from arguments. */
    else {
	/* START # 4 */
	/* List options. */
	/* -h or -? help page */
	/* -fname - name = mged file name. */
	/* -n# - # = number of gas tanks. */
	/* -H# - # = height of gas tank in mm. */
	/* -w# - # = width of gas tank in mm. */
	/* -d# - # = depth of gas tank in mm. */
	/* -r# - # = radius of corners in mm. */

	for (i=1; i<argc; i++) {
	    /* START # 5 */
	    /* Put argument in temporary character string. */
	    temp = argv[i];

	    if (temp[1] == 'h' || temp[1] == '?') {
	    	explain();
		bu_exit(2,NULL);
	    }

	    /* -f - mged file. */
	    if (temp[1] == 'f') {
		/* START # 6 */
		j = 2;
		k = 0;
		while ((temp[j] != '\0') && (k < NAME_LEN)) {
		    /* START # 7 */
		    filemged[k] = temp[j];
		    j++;
		    k++;
		}					/* END # 7 */
		filemged[k] = '\0';
	    }						/* END # 6 */

	    /* All other options. */
	    else {
		/* START # 8 */
		/* Set up temporary character string. */
		j = 2;
		k = 0;
		while ((temp[j] != '\0') && (k < NAME_LEN)) {
		    /* START # 9 */
		    temp1[k] = temp[j];
		    j++;
		    k++;
		}					/* END # 9 */
		temp1[k] = '\0';
		if (temp[1] == 'n') {
		    sscanf(temp1, "%d", &numtnk);
		    if (numtnk < 1)
			numtnk = 1;
		    else if (numtnk > maxnumtnk)
			numtnk = maxnumtnk;
		} else if (temp[1] == 'H') {
		    sscanf(temp1, "%lf", &hgt);
		} else if (temp[1] == 'w') {
		    sscanf(temp1, "%lf", &wid);
		} else if (temp[1] == 'd') {
		    sscanf(temp1, "%lf", &dpt);
		} else if (temp[1] == 'r') {
		    sscanf(temp1, "%lf", &rds);
		}
	    }						/* END # 8 */
	}						/* END # 5 */
    }							/* END # 4 */

    /* Print out all info. */
    printf("\nmged file:  %s\n", filemged);
    printf("height of gas tank:  %f mm\n", hgt);
    printf("width of gas tank:  %f mm\n", wid);
    printf("depth of gas tank:  %f mm\n", dpt);
    printf("radius of corner:  %f mm\n", rds);
    printf("number of gas tanks:  %d\n\n", numtnk);
    (void)fflush(stdout);

    /* Open mged file. */
    fpw = wdb_fopen(filemged);

    /* Write ident record. */
    mk_id(fpw, "windows");

    for (i=0; i<numtnk; i++) {
	/* START # 2 */

	/* Create all solids. */

	/* Create the 3 arb8s. */

	pts[0][0] = (fastf_t)(dpt / 2.0);
	pts[0][1] = (fastf_t)(wid / 2.0 - rds);
	pts[0][2] = (fastf_t)(hgt / 2.0 - rds);
	pts[1][0] = pts[0][0];
	pts[1][1] = pts[0][1];
	pts[1][2] = (-pts[0][2]);
	pts[2][0] = pts[0][0];
	pts[2][1] = (-pts[0][1]);
	pts[2][2] = pts[1][2];
	pts[3][0] = pts[0][0];
	pts[3][1] = pts[2][1];
	pts[3][2] = pts[0][2];
	pts[4][0] = (-pts[0][0]);
	pts[4][1] = pts[0][1];
	pts[4][2] = pts[0][2];
	pts[5][0] = pts[4][0];
	pts[5][1] = pts[0][1];
	pts[5][2] = (-pts[0][2]);
	pts[6][0] = pts[4][0];
	pts[6][1] = (-pts[0][1]);
	pts[6][2] = pts[1][2];
	pts[7][0] = pts[4][0];
	pts[7][1] = pts[2][1];
	pts[7][2] = pts[0][2];
	solnam[5] = 97 + i;
	solnam[6] = '0';
	solnam[7] = '1';
	mk_arb8(fpw, solnam, &pts[0][X]);

	pts[0][0] = (fastf_t)(dpt / 2.0 - rds);
	pts[0][1] = (fastf_t)(wid / 2.0);
	pts[1][0] = pts[0][0];
	pts[1][1] = pts[0][1];
	pts[2][0] = pts[0][0];
	pts[2][1] = (-pts[0][1]);
	pts[3][0] = pts[0][0];
	pts[3][1] = pts[2][1];
	pts[4][0] = (-pts[0][0]);
	pts[4][1] = pts[0][1];
	pts[5][0] = pts[4][0];
	pts[5][1] = pts[0][1];
	pts[6][0] = pts[4][0];
	pts[6][1] = pts[2][1];
	pts[7][0] = pts[4][0];
	pts[7][1] = pts[2][1];
	solnam[7] = '2';
	mk_arb8(fpw, solnam, &pts[0][X]);

	pts[0][1] = (fastf_t)(wid / 2.0 - rds);
	pts[0][2] = (fastf_t)(hgt / 2.0);
	pts[1][1] = pts[0][1];
	pts[1][2] = (-pts[0][2]);
	pts[2][1] = (-pts[0][1]);
	pts[2][2] = pts[1][2];
	pts[3][1] = pts[2][1];
	pts[3][2] = pts[0][2];
	pts[4][1] = pts[0][1];
	pts[4][2] = pts[0][2];
	pts[5][1] = pts[0][1];
	pts[5][2] = (-pts[0][2]);
	pts[6][1] = (-pts[0][1]);
	pts[6][2] = pts[1][2];
	pts[7][1] = pts[2][1];
	pts[7][2] = pts[0][2];
	solnam[7] = '3';
	mk_arb8(fpw, solnam, &pts[0][X]);

	/* Make 8 spheres. */

	cent[0] = (fastf_t)(dpt / 2.0 - rds);
	cent[1] = (fastf_t)(wid / 2.0 - rds);
	cent[2] = (fastf_t)(hgt / 2.0 - rds);
	rad = (fastf_t)(rds);
	solnam[7] = '4';
	mk_sph(fpw, solnam, cent, rad);

	cent[2] = (-cent[2]);
	solnam[7] = '5';
	mk_sph(fpw, solnam, cent, rad);

	cent[1] = (-cent[1]);
	solnam[7] = '6';
	mk_sph(fpw, solnam, cent, rad);

	cent[2] = (-cent[2]);
	solnam[7] = '7';
	mk_sph(fpw, solnam, cent, rad);

	cent[0] = (-cent[0]);
	cent[1] = (-cent[1]);
	solnam[7] = '8';
	mk_sph(fpw, solnam, cent, rad);

	cent[2] = (-cent[2]);
	solnam[7] = '9';
	mk_sph(fpw, solnam, cent, rad);

	cent[1] = (-cent[1]);
	solnam[6] = '1';
	solnam[7] = '0';
	mk_sph(fpw, solnam, cent, rad);

	cent[2] = (-cent[2]);
	solnam[7] = '1';
	mk_sph(fpw, solnam, cent, rad);

	/* Make 12 cylinders. */

	bs[0] = (fastf_t)(dpt / 2.0 - rds);
	bs[1] = (fastf_t)(wid / 2.0 - rds);
	bs[2] = (fastf_t)(hgt / 2.0 - rds);
	ht[0] = (fastf_t)(0.0);
	ht[1] = (fastf_t)(-wid + 2 * rds);
	ht[2] = (fastf_t)(0.0);
	solnam[7] = '2';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[2] = (-bs[2]);
	solnam[7] = '3';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[0] = (-bs[0]);
	solnam[7] = '4';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[2] = (-bs[2]);
	solnam[7] = '5';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[0] = (fastf_t)(dpt / 2.0 - rds);
	bs[1] = (fastf_t)(wid / 2.0 - rds);
	bs[2] = (fastf_t)(hgt / 2.0 - rds);
	ht[0] = (fastf_t)(0.0);
	ht[1] = (fastf_t)(0.0);
	ht[2] = (fastf_t)(-hgt + 2 * rds);
	solnam[7] = '6';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[1] = (-bs[1]);
	solnam[7] = '7';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[0] = (-bs[0]);
	solnam[7] = '8';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[1] = (-bs[1]);
	solnam[7] = '9';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[0] = (fastf_t)(dpt / 2.0 - rds);
	bs[1] = (fastf_t)(wid / 2.0 - rds);
	bs[2] = (fastf_t)(hgt / 2.0 - rds);
	ht[0] = (fastf_t)(-dpt + 2 * rds);
	ht[1] = (fastf_t)(0.0);
	ht[2] = (fastf_t)(0.0);
	solnam[6] = '2';
	solnam[7] = '0';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[2] = (-bs[2]);
	solnam[7] = '1';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[1] = (-bs[1]);
	solnam[7] = '2';
	mk_rcc(fpw, solnam, bs, ht, rad);

	bs[2] = (-bs[2]);
	solnam[7] = '3';
	mk_rcc(fpw, solnam, bs, ht, rad);

	/* Make all regions. */

	/* Initialize list. */
	BU_LIST_INIT(&comb.l);

	/* Region 1. */
	solnam[5] = 97 + i;
	solnam[6] = '0';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[5] = 97 + i;
	regnam[6] = '0';
	regnam[7] = '1';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 2. */
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[7] = '2';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 3. */
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	regnam[7] = '3';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 4. */
	solnam[7] = '4';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '1';
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '6';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[6] = '2';
	solnam[7] = '0';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[7] = '4';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 5. */
	solnam[6] = '0';
	solnam[7] = '5';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '1';
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '6';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[6] = '2';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[7] = '5';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 6. */
	solnam[6] = '0';
	solnam[7] = '6';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '1';
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '7';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[6] = '2';
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[7] = '6';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 7. */
	solnam[6] = '0';
	solnam[7] = '7';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '1';
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '7';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[6] = '2';
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[7] = '7';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 8. */
	solnam[6] = '0';
	solnam[7] = '8';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '1';
	solnam[7] = '5';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '9';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[6] = '2';
	solnam[7] = '0';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[7] = '8';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 9. */
	solnam[6] = '0';
	solnam[7] = '9';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '1';
	solnam[7] = '4';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '9';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[6] = '2';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[7] = '9';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 10. */
	solnam[6] = '1';
	solnam[7] = '0';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[7] = '4';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '8';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[6] = '2';
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '1';
	regnam[7] = '0';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 11. */
	solnam[6] = '1';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[7] = '5';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '8';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[6] = '2';
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '1';
	regnam[7] = '1';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 12. */
	solnam[6] = '1';
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '1';
	regnam[7] = '2';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 13. */
	solnam[6] = '1';
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '1';
	regnam[7] = '3';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 14. */
	solnam[6] = '1';
	solnam[7] = '4';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '1';
	regnam[7] = '4';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 15. */
	solnam[6] = '1';
	solnam[7] = '5';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '1';
	regnam[7] = '5';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 16. */
	solnam[6] = '1';
	solnam[7] = '6';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '1';
	regnam[7] = '6';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 17. */
	solnam[6] = '1';
	solnam[7] = '7';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '1';
	regnam[7] = '7';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 18. */
	solnam[6] = '1';
	solnam[7] = '8';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '1';
	regnam[7] = '8';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 19. */
	solnam[6] = '1';
	solnam[7] = '9';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '1';
	regnam[7] = '9';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 20. */
	solnam[6] = '2';
	solnam[7] = '0';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '2';
	regnam[7] = '0';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 21. */
	solnam[6] = '2';
	solnam[7] = '1';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '2';
	regnam[7] = '1';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 22. */
	solnam[6] = '2';
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '2';
	regnam[7] = '2';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Region 23. */
	solnam[6] = '2';
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_INTERSECT);
	solnam[6] = '0';
	solnam[7] = '2';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	solnam[7] = '3';
	(void)mk_addmember(solnam, &comb.l, NULL, WMOP_SUBTRACT);
	regnam[6] = '2';
	regnam[7] = '3';
	mk_lfcomb(fpw, regnam, &comb, 1);

	/* Create group. */

	/* Initialize list. */
	BU_LIST_INIT(&comb1.l);

	regnam[6] = '0';
	regnam[7] = '1';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '2';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '3';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '4';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '5';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '6';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '7';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '8';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '9';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[6] = '1';
	regnam[7] = '0';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '1';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '2';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '3';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '4';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '5';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '6';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '7';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '8';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '9';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[6] = '2';
	regnam[7] = '0';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '1';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '2';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	regnam[7] = '3';
	(void)mk_addmember(regnam, &comb1.l, NULL, WMOP_UNION);
	grpnam[3] = 97 + i;
	mk_lfcomb(fpw, grpnam, &comb1, 0);

    }							/* START # 2 */

    /* Close file. */
    db_close(fpw->dbip);
    return 0;
}							/* END # 1 */

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
