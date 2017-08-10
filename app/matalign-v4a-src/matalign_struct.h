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


/* Structure for holding information about a profile.
 */
typedef struct struct_PROFILE
{
  char           *mat_name;	/* Name of the matrix file. */
  int           **int_matrix;	/* Integer count matrix, first dimension width,
  				 * second dimension alphabet. */
  double        **f_matrix; 	/* Frequency matrix */
  int      	 *dcode;	/* Degenerate code of the profile */
  int      	  width;	/* The width of the matrix */
  int     	 *num_seq;	/* Number of sequences that form the profile */
  struct PROFILE *mirror;	/* Reverse complement profile */
} PROFILE;


/* Structure for holding information about an HSP
 */
typedef struct struct_HSP
{
  PROFILE      *profile_1;  /* Parent profile that gives rise to the HSP */
  PROFILE      *profile_2;  /* Parent profile that gives rise to the HSP */
  double        hsp_score;  /* hsp score of the HSP */
  int           m;	    /* Length of profile_1 */
  int           n;	    /* Length of profile_2 */
  int           i_start;    /* Start index of HSP in profile_1 */
  int           j_start;    /* Start index of HSP in profile_2 */
  int           i_end;	    /* End index of HSP in profile_1 */
  int           j_end;	    /* End index of HSP in profile_2 */
  int           length;     /* Length of the HSP */
} HSP;

/* Structure for holding coordinates of a point
 */
typedef struct struct_COORD
{
  int   i;
  int   j;
} COORD;

typedef struct struct_MATRIX
{
  int **matrix;
  int   i;
  int   j;
} MATRIX;


/* Structure for holding a column of a matrices. */
typedef struct struct_COLUMN
{
  int    *column;      /* The column of integers. */
  double *frequency;   /* The column of frequencies */
  int     num_seqs;
  double  info;
} COLUMN;

