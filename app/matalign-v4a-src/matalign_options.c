/****************************************************************
 * MatrixAligner (version 4a)					*
 * Copyright (C) 2003, 2004, 2005, 2006                         *
 * Washington University, St. Louis, MO, USA.                   *
 * All Rights Reserved.						*
 * Author:							*
 *  Ting Wang                                                   *
 *  Laboratory of Dr. Gary D. Stormo	                        *
 *  Department of Genetics					*
 *  Washington University in St. Louis				*
 *  Campus Box 8232						*
 *  St. Louis MO 63110						*
 *								*
 *  twang@ural.wustl.edu					*
 ****************************************************************/


/* External global variables corresponding to various options,
 * the matrix of sequences, and the list of the summary matrices. */

#include <stdio.h>
#include "matalign_options.h"



int   PID;

int   Global      = 0;
int   Comp_status = 1;
int   EZ          = 0;
int   List_num    = 0;

char  *file1 = (char *)NULL;		/* Name of matrix description file 1 */
char  *file2 = (char *)NULL; 		/* Name of matrix description file 2 */
char  *path1 = (char *)NULL;		/* Directory for matrices in list 1 */
char  *path2 = (char *)NULL;		/* Directory for matrices in list 2 */
char **Mat_names_1 = (char **)NULL;	/* Matrix names of list 1 */
char **Mat_names_2 = (char **)NULL;  	/* Matrix names of list 2 */
int    Mat_num_1 = 0;			/* Number of matrices in list 1 */
int    Mat_num_2 = 0;			/* Number of matrices in list 2 */


PROFILE **Profiles_1;	/* Matrix list 1 */
PROFILE **Profiles_2;	/* Matrix list 2 */
HSP     **Hsp_Q;	/* List of HSPs from comparisons */

/* Information for ALLR statistic */
char   Dcode[15];		/* Degenerate code array */
double Dcode_f[15];		/* Degenerate letter frequency */
int    SM1[15][15];		/* ALLR scoring matrix 1 */
int    SM2[15][15];		/* ALLR scoring matrix 2 */

/* Information for Karlin/Altschul statistic */
double M;
double N;
double K;			/* Parameter K */
double Lambda;			/* Parameter Lambda */
double H;			/* Parameter H */


