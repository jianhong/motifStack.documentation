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


#include "matalign_options.h"

int Hsp_num;

void     compare_profiles();
HSP     *compare_2_profiles( PROFILE *, PROFILE * );
double **score_function( PROFILE *, PROFILE * );
void     calculate_distance( double *, double *, HSP * );
PROFILE *merge_HSP( HSP * );
void     print_HSP_info();
void     print_consensus( PROFILE * );
double   get_hsp_E_value ( HSP * );
double   e_value_to_p_value( double );
void     free_HSP( HSP * );
void     free_scores( double **, int m, int n );
void     free_PROFILE( PROFILE * );

/* Main function for comparing matrices */
void compare_profiles()
{
  int  i, j;
  int  hsp_cnt;
  HSP *hsp_f;
  HSP *hsp_r;
  
  HSP *compare_2_profiles();
  void print_HSP_info();
  
  /* Determine types of comparison */
  if ( List_num == 1 )
  {
    Hsp_num = Mat_num_1 * (Mat_num_1-1) /2;
  }
  else if ( List_num == 2 )
  {
    Hsp_num = Mat_num_1 * Mat_num_2;
  }
  
  Hsp_Q = (HSP **)tw_calloc( Hsp_num,
  			     sizeof(HSP *),
  			     "Hsp_Q",
  			     "compare_profiles()" );
  
  hsp_cnt = 0;
  hsp_f   = (HSP *)NULL;
  hsp_r   = (HSP *)NULL;
  
  /* One list, pair-wise compare all matrices against each other */
  if ( List_num == 1 )
  {
    for ( i=0; i<Mat_num_1-1; ++i )
    {
      for ( j=i+1; j<Mat_num_1; ++j )
      {
      	if ( Comp_status )
      	{
      	  hsp_f = compare_2_profiles( Profiles_1[i], Profiles_1[j] );
      	  hsp_r = compare_2_profiles( Profiles_1[i], Profiles_1[j]->mirror );
      	  if ( hsp_f->hsp_score >= hsp_r->hsp_score )
      	  {
      	    Hsp_Q[hsp_cnt] = hsp_f;
      	    free_HSP( hsp_r );
      	  }
      	  else
      	  {
      	    Hsp_Q[hsp_cnt] = hsp_r;
      	    free_HSP( hsp_f );
      	  }
      	}
      	else
      	{
      	  Hsp_Q[hsp_cnt] = compare_2_profiles( Profiles_1[i], Profiles_1[j] );
      	}
      	hsp_cnt++;
      }
    }
  }
  /* Two lists, compare all matrices in list 1 to all matrices in list 2 */
  else if ( List_num == 2 )
  {
    for ( i=0; i<Mat_num_1; ++i )
    {
      for ( j=0; j<Mat_num_2; ++j )
      {
        if ( Comp_status )
      	{
      	  hsp_f = compare_2_profiles( Profiles_1[i], Profiles_2[j] );
      	  hsp_r = compare_2_profiles( Profiles_1[i], Profiles_2[j]->mirror );
      	  if ( hsp_f->hsp_score >= hsp_r->hsp_score )
      	  {
      	    Hsp_Q[hsp_cnt] = hsp_f;
      	    free_HSP( hsp_r );
      	  }
      	  else
      	  {
      	    Hsp_Q[hsp_cnt] = hsp_r;
      	    free_HSP( hsp_f );
      	  }
      	}
      	else
      	{
      	  Hsp_Q[hsp_cnt] = compare_2_profiles( Profiles_1[i], Profiles_2[j] );
      	}
      	hsp_cnt++;
      }
    }
  }
  
  /* Check number of comparisons */
  if ( hsp_cnt != Hsp_num )
  {
    printf ( "Warning: HSP numbers not consistent -- %d vs %d\n", hsp_cnt, Hsp_num );
  }
  
  /* Generate output table */
  print_HSP_info();

}


/* Compare two profiles, generate one HSP */
HSP *compare_2_profiles( PROFILE *profile_1, PROFILE *profile_2 )
{
  int       width_1;
  int       width_2;
  int       i;
  
  /* Variables to hold scores 
   */
  double  **score;
  double  **dp_score;
  
  /* Variables to hold information about HSPs 
   */
  HSP     *hsp;
  
  double **score_function();
  double **DP_function();
  double **DP_function_g();
  HSP     *trace_back();
  HSP     *trace_back_g();
  void     free_scores();
  
  width_1    = profile_1->width;
  width_2    = profile_2->width;

  /* Calculate pairwise comparison scores */
  score = score_function( profile_1, profile_2 );
  
  /* Dynamic programming function */
  if ( Global == 0 )
  {
    dp_score = DP_function( score, width_1, width_2 );
  }
  else
  {
    dp_score = DP_function_g( score, width_1, width_2 );
  }
  
  /* HSP generated from this comparison */
  if ( Global == 0 )
  {
    hsp = trace_back( dp_score, score, width_1, width_2 );
  }
  else
  {
    hsp = trace_back_g( dp_score, score, width_1, width_2 );
  }
  hsp->profile_1 = profile_1;
  hsp->profile_2 = profile_2;

  /* Free the space used for comparison */
  free_scores( score, width_1, width_2 );
  free_scores( dp_score, width_1, width_2 );
  
  return hsp;
  
}  


/* Calculate pair-wise position score */
double **score_function( PROFILE *profile_1, PROFILE *profile_2 )
{
  int      i, j;
  double **score;
  int    **int_matrix_1;
  int    **int_matrix_2;
  double **f_matrix_1;
  double **f_matrix_2;
  int     *num_seq_1;
  int     *num_seq_2;
  int      m, n;
  
  double   ALLR_function();
  double   ALLR_score();
  double   calculate_ALLR_double();
  
  int_matrix_1 = profile_1->int_matrix;
  int_matrix_2 = profile_2->int_matrix;
  f_matrix_1   = profile_1->f_matrix;
  f_matrix_2   = profile_2->f_matrix;
  num_seq_1    = profile_1->num_seq;
  num_seq_2    = profile_2->num_seq;
  m            = profile_1->width;
  n            = profile_2->width;
  
  
  score = (double **)tw_calloc( m, 
                                sizeof(double *),
				"STscore",
				"score_function()" );

  for ( i=0; i<m; ++i )
  {
    score[i] = (double *)tw_calloc( n,
    				    sizeof(double),
				    "STscore[]",
				    "score_function()" );
  }
    
  for ( i=0; i<m; ++i )
  {
    for ( j=0; j<n; ++j )
    { 
      score[i][j] = ALLR_function( int_matrix_1[i], f_matrix_1[i], num_seq_1[i],
      				   int_matrix_2[j], f_matrix_2[j], num_seq_2[j] );
      //score[i][j] = ALLR_score( int_matrix_1[i], int_matrix_2[j] );
      //score[i][j] = calculate_ALLR_double( f_matrix_1[i], num_seq_1[i], 
      //						f_matrix_2[j], num_seq_2[j] );
    }
  }
  return ( score );
}


/* Dynamic programming function, local version */
double **DP_function( double **score, int m, int n )
{
  int        i, j;
  double     s = 0;
  double   **DP;
  
  DP = (double **)tw_calloc( m, sizeof(double *), "DP", "DP_function" );
  
  for ( i=0; i<m; ++i )
  {
    DP[i] = (double *)tw_calloc( n,
    				 sizeof(double),
				 "DP[]",
				 "DP_function()" );
  }     
  
  for ( i=0; i<m; ++i )
  {
    DP[i][0] = (score[i][0]>=0) ? score[i][0] : 0;
  }
  for ( j=0; j<n; ++j )
  {
    DP[0][j] = (score[0][j]>=0) ? score[0][j] : 0;
  }
  
  for ( i=1; i<m; ++i )
  {
    for ( j=1; j<n; ++j )
    {
      s = DP[i-1][j-1] + score[i][j];
      DP[i][j] = (s>=0) ? s : 0;
    }
  }
  return ( DP );
}

/* Global alignment version */
double **DP_function_g( double **score, int m, int n )
{
  int        i, j;
  double     s = 0;
  double   **DP;
  
  DP = (double **)tw_calloc( m, sizeof(double *), "DP", "DP_function" );
  
  for ( i=0; i<m; ++i )
  {
    DP[i] = (double *)tw_calloc( n,
    				 sizeof(double),
				 "DP[]",
				 "DP_function_g()" );
  }     
  
  for ( i=0; i<m; ++i )
  {
    DP[i][0] = score[i][0]; //(score[i][0]>=0) ? score[i][0] : 0;
  }
  for ( j=0; j<n; ++j )
  {
    DP[0][j] = score[0][j]; //(score[0][j]>=0) ? score[0][j] : 0;
  }
  
  for ( i=1; i<m; ++i )
  {
    for ( j=1; j<n; ++j )
    {
      s = DP[i-1][j-1] + score[i][j];
      DP[i][j] = s;
    }
  }
  return ( DP );
}


/* Trace back */
HSP *trace_back( double **DP, double **score, int m, int n )
{
  int       i, j;
  double    max;
  int       i_start, j_start;
  int       i_end, j_end;
  HSP      *hsp;
  
  hsp = (HSP *)tw_calloc( 1, sizeof(HSP), "hsp", "trace_back()" );
  
  max = 0;
  
  for ( i=0; i<m; ++i )
  {
    for ( j=0; j<n; ++j )
    {
      if ( DP[i][j]>=max )
      {
        max = DP[i][j];
	i_end = i;
	j_end = j;
      }
    }
  }

  i = i_end;
  j = j_end;
  
  while ( (i>0) && (j>0) )
  {
    if ( DP[i][j]>0 )
    {
      --i;
      --j;
    }
    else
    {
      ++i;
      ++j;
      break;
    }
  }
  
  i_start = i;
  j_start = j;
  
  hsp->profile_1 = NULL;
  hsp->profile_2 = NULL;
  hsp->hsp_score = max;
  hsp->m = m;
  hsp->n = n;
  hsp->i_start = i_start;
  hsp->j_start = j_start;
  hsp->i_end = i_end;
  hsp->j_end = j_end;
  hsp->length = i_end - i_start + 1;
  
  /*
  hsp->match = (int *)tw_calloc( hsp->length,
  				 sizeof(int),
  				 "hsp->match",
  				 "trace_back()" );
  for ( i=0; i<hsp->length; ++i )
  {
    if ( score[ i_start+i ][ j_start+i ] >= 0 )
    {
      hsp->match[i] = 1;
    }
    else
    {
      hsp->match[i] = 0;
    }
  }
  */
  
  return hsp;
}

/* Trace back, global version */
HSP *trace_back_g( double **DP, double **score, int m, int n )
{
  int       i, j;
  double    max;
  int       i_start, j_start;
  int       i_end, j_end;
  HSP      *hsp;
  
  hsp = (HSP *)tw_calloc( 1, sizeof(HSP), "hsp", "trace_back()" );
  
  max = DP[m-1][n-1];
  i_end = m-1;
  j_end = n-1;
  for ( i=0; i<m; ++i )
  {
    if ( DP[i][n-1] >= max )
    {
      max = DP[i][n-1];
      i_end = i;
      j_end = n-1;
    }
  }
  for ( j=0; j<n; ++j )
  {
    if ( DP[m-1][j] >= max )
    {
      max = DP[m-1][j];
      i_end = m-1;
      j_end = j;
    }
  }
  
  i = i_end;
  j = j_end;
  
  while ( (i>0) && (j>0) )
  {
    --i;
    --j;
  }
  
  i_start = i;
  j_start = j;
  
  
  hsp->profile_1 = NULL;
  hsp->profile_2 = NULL;
  hsp->hsp_score = max;
  hsp->m = m;
  hsp->n = n;
  hsp->i_start = i_start;
  hsp->j_start = j_start;
  hsp->i_end = i_end;
  hsp->j_end = j_end;
  hsp->length = i_end - i_start + 1;
  
  /*
  hsp->match = (int *)tw_calloc( hsp->length,
  				 sizeof(int),
  				 "hsp->match",
  				 "trace_back_g()" );
  for ( i=0; i<hsp->length; ++i )
  {
    if ( score[ i_start+i ][ j_start+i ] > 0 )
    {
      hsp->match[i] = 1;
    }
    else
    {
      hsp->match[i] = 0;
    }
  }
  */
  
  return hsp;
}


/* Generate output data table */
void print_HSP_info()
{
  int i;
  double   dist1;
  double   dist2;
  PROFILE *new_profile;
  double   e_v;
  double   p_v;
  
  void print_consensus();
  double get_ALLR_E_value();
  double e_value_to_p_value();
  
  printf ( "Number of pair-wise comparisons: %d\n", Hsp_num );
  
  printf ( "%s\t", "Matrix_1" );
  if ( !EZ )
  {
    printf ( "%s\t", "Consensus" );
  }
  printf ( "%s\t", "Matrix_2" );
  if ( !EZ )
  {
    printf ( "%s\t", "Consensus" );
  }
  printf ( "%s\t", "ALLR" );
  printf ( "%s\t", "Overlap" );
  
  if ( !EZ )
  {
    printf ( "%s\t", "Distance" );
    printf ( "%s\t", "Aligned_Dist" );
    printf ( "%s\t", "Shared_Consensus" );
    printf ( "%s\t", "E_value" );
    printf ( "%s\t", "P_value" );
  }
  printf ( "\n" );

  
  for ( i=0; i<Hsp_num; i++ )
  {
    printf ( "%s\t", Hsp_Q[i]->profile_1->mat_name );
    if ( !EZ )
    {
      print_consensus( Hsp_Q[i]->profile_1 );
      printf ( "\t" );
    }
    printf ( "%s\t", Hsp_Q[i]->profile_2->mat_name );
    if ( !EZ )
    {
      print_consensus( Hsp_Q[i]->profile_2 );
      printf ( "\t" );
    }
    printf ( "%.3f\t", Hsp_Q[i]->hsp_score );
    printf ( "%d\t", Hsp_Q[i]->length );
    
    if ( !EZ )
    {
      /* calculate distances */
      calculate_distance( &dist1, &dist2, Hsp_Q[i] );
      printf ( "%.3f\t", dist1 );
      printf ( "%.3f\t", dist2 );
      
      /* generate a new matrix, print consensus sequence */
      new_profile = merge_HSP ( Hsp_Q[i] );
      print_consensus( new_profile );
      printf ( "\t" );
      free_PROFILE( new_profile );
      
      /* calculate e_value and p_value */
      e_v = get_hsp_E_value( Hsp_Q[i] );
      p_v = e_value_to_p_value( e_v );
      printf ( "%.4g\t", e_v );
      printf ( "%.4g\t", p_v );
      
      //printf ( "%d\t%d", Hsp_Q[i]->i_start, Hsp_Q[i]->j_start );
    }
    
    printf ( "\n" );
  }
}

/* print the consensus pattern of a profile */
void print_consensus( PROFILE *profile )
{
  int i;
  
  for ( i=0; i<profile->width; ++i )
  {
    printf ( "%c", Dcode[ profile->dcode[i] ] );
  }
}

/* Calculate distances between two profiles */
void calculate_distance( double *dist1, double *dist2, HSP *hsp )
{
  PROFILE *profile_1;
  PROFILE *profile_2;
  double   full_1;
  double   full_2;
  double   aligned_1;
  double   aligned_2;
  int      i, j;
  
  double   ALLR_function();
  
  profile_1 = hsp->profile_1;
  profile_2 = hsp->profile_2;
  full_1    = 0;
  full_2    = 0;
  aligned_1 = 0;
  aligned_2 = 0;
  
  for ( i=0; i<profile_1->width; ++i )
  {
    full_1 += ALLR_function( profile_1->int_matrix[i],
    			     profile_1->f_matrix[i],
    			     profile_1->num_seq[i],
    			     profile_1->int_matrix[i],
    			     profile_1->f_matrix[i],
    			     profile_1->num_seq[i]
    			   );
  }
  for ( i=hsp->i_start; i<=hsp->i_end; i++ )
  {
    aligned_1 += ALLR_function( profile_1->int_matrix[i],
    			     	profile_1->f_matrix[i],
    			     	profile_1->num_seq[i],
    			     	profile_1->int_matrix[i],
    			     	profile_1->f_matrix[i],
    			     	profile_1->num_seq[i]
    			      );
  }
  for ( i=0; i<profile_2->width; ++i )
  {
    full_2 += ALLR_function( profile_2->int_matrix[i],
    			     profile_2->f_matrix[i],
    			     profile_2->num_seq[i],
    			     profile_2->int_matrix[i],
    			     profile_2->f_matrix[i],
    			     profile_2->num_seq[i]
    			   );
  }
  for ( i=hsp->j_start; i<=hsp->j_end; i++ )
  {
    aligned_2 += ALLR_function( profile_2->int_matrix[i],
    			     	profile_2->f_matrix[i],
    			     	profile_2->num_seq[i],
    			     	profile_2->int_matrix[i],
    			     	profile_2->f_matrix[i],
    			     	profile_2->num_seq[i]
    			      );
  }
  *dist1 = full_1 + full_2 - 2*hsp->hsp_score;
  *dist2 = aligned_1 + aligned_2 - 2*hsp->hsp_score;
  
}


/* Create a profile based on a HSP */     
PROFILE *merge_HSP( HSP * hsp )
{
  int        i, i1, i2, j;
  PROFILE   *new_profile;
  PROFILE   *profile_1;
  PROFILE   *profile_2;
  int      **int_matrix;
  double   **f_matrix;
  int       *dcode;
  int        width;
  int       *num_seq;

  int map_column_dcode();

  new_profile = (PROFILE *)tw_malloc( sizeof(PROFILE),
                                      "new_profile", 
				      "merge_HSP()" );
  width      = hsp->length;
  profile_1  = hsp->profile_1;
  profile_2  = hsp->profile_2;
  int_matrix = (int **)tw_calloc( width,
  				  sizeof(int *),
  				  "int_matrix",
  				  "merge_HSP()" );
  f_matrix   = (double **)tw_calloc( width,
  				     sizeof(double *),
  				     "f_matrix",
  				     "merge_HSP()" );
  dcode      = (int *)tw_calloc( width,
  				 sizeof(int),
  				 "dcode",
  				 "merge_HSP()" );
  num_seq    = (int *)tw_calloc( width,
  				 sizeof(int),
  				 "num_seq",
  				 "merge_HSP()" );
  for ( i=0; i<width; ++i )
  {
    int_matrix[i] = (int *)tw_calloc( A_size,
    				      sizeof(int),
    				      "int_matrix[]",
    				      "merge_HSP()" );
    f_matrix[i] = (double *)tw_calloc( A_size,
    				       sizeof(double),
    				       "f_matrix[]",
    				       "merge_HSP()" );
  }
  
  for ( i=0, i1=hsp->i_start, i2=hsp->j_start;
        i<width; 
	++i, ++i1, ++i2 )
  {
    num_seq[i] = profile_1->num_seq[i1]+profile_2->num_seq[i2];
    for ( j=0; j<A_size; ++j )
    {
      int_matrix[i][j] = profile_1->int_matrix[i1][j]+profile_2->int_matrix[i2][j];
      f_matrix[i][j] = ( int_matrix[i][j]+P[j]*(double)PSEUDOCNT )/(num_seq[i]+(double)PSEUDOCNT );
    }
    dcode[i] = map_column_dcode( f_matrix[i] );
  }
  
  new_profile->mat_name   = "New";
  new_profile->int_matrix = int_matrix;
  new_profile->f_matrix   = f_matrix;
  new_profile->dcode      = dcode;
  new_profile->num_seq    = num_seq;
  new_profile->width      = width;
  new_profile->mirror     = (PROFILE *)NULL;
  
  return new_profile;
}


/* Free score matrix space */
void free_scores( double **score, int m, int n )
{
  int i, j;
  for ( i=0; i<m; ++i )
  {
    tw_free ( score[i], "score[]", "free_scores()" );
  }
  tw_free( score, "score", "free_scores()" );
}    

/* Free HSP structure */
void free_HSP( HSP *hsp )
{
  hsp->profile_1 = (PROFILE *)NULL;
  hsp->profile_2 = (PROFILE *)NULL;
  tw_free( hsp, "hsp", "free_HSP()" );
}

/* Free PROFILE structure */
void free_PROFILE( PROFILE *profile )
{
  int i;
  PROFILE *mirror;
  
  mirror = profile->mirror;
  //if ( profile->mat_name )
  //{
  //  tw_free( profile->mat_name, "mat_name", "free_PROFILE()" );
  //}
  for ( i=0; i<profile->width; ++i )
  {
    tw_free( profile->int_matrix[i], "int_matrix[]", "free_PROFILE()" );
    tw_free( profile->f_matrix[i], "f_matrix[]", "free_PROFILE()" );
  }
  tw_free( profile->int_matrix, "int_matrix", "free_PROFILE()" );
  tw_free( profile->f_matrix, "f_matrix", "free_PROFILE()" );
  tw_free( profile->dcode, "dcode", "free_PROFILE()" );
  tw_free( profile->num_seq, "dcode", "free_PROFILE()" );
  
  tw_free( profile, "profile", "free_PROFILE()" );
  
  if ( mirror )
  {
    free_PROFILE( mirror );
  }
}

