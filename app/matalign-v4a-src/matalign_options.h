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


#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "tw_alloc.h"
#include "tw_error.h"
#include "tw_file.h"
#include "matalign_struct.h"
#include "alpha.h"


#define NO 0
#define YES 1
#define LN2 0.6931471805599453094172321214581765680755   /* ln(2) */
#define MAXLINE 4096
#define CHUNK_SIZE 32
#define BIG_CHUNK_SIZE 512
#define PSEUDOCNT 1.0

extern int       PID;
extern int       Global;
extern int       Comp_status;
extern int       EZ;
extern int       List_num;
extern char     *file1;
extern char     *file2;
extern char     *path1;
extern char     *path2;
extern char    **Mat_names_1;
extern char    **Mat_names_2;
extern int       Mat_num_1;
extern int       Mat_num_2;
extern PROFILE **Profiles_1;	/* Matrix list 1 */
extern PROFILE **Profiles_2;	/* Matrix list 2 */
extern HSP     **Hsp_Q;		/* List of HSPs from comparisons */

/* Information for ALLR statistic */
extern char   Dcode[15];
extern double Dcode_f[15];
extern int    SM1[15][15];
extern int    SM2[15][15];

/* Information for Karlin/Altschul statistic */
extern double M;
extern double N;
extern double K;
extern double Lambda;
extern double H;

