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

char Dcode[15] = 
	{ 'A', 'T', 'C', 'G', 'M', 'K', 'R', 'Y', 'W', 'S', 'a', 't', 'c', 'g', 'n' };

double Dcode_f[15] = 
	{ 0.0877,0.0850,0.0418,0.0418,0.0146,0.0142,0.0311,0.0313,0.0221,0.0108,0.1328,0.1312,0.0883,0.0879,0.1795 };

int SM1[15][15] = 
	{
        {  144, -399, -400, -400,  -35, -396,  -35, -396,  -34, -397,   60, -309, -319, -288, -167 },
        { -399,  144, -400, -400, -396,  -35, -396,  -35,  -34, -397, -309,   59, -288, -320, -168 },
        { -400, -400,  195, -401,    3, -398, -398,    3, -398,    2, -324, -287,  105, -327, -163 },
        { -400, -400, -401,  195, -398,    3,    3, -398, -397,    2, -288, -326, -329,  104, -165 },
        {  -35, -396,    3, -398,   84, -394, -169, -142, -168, -143,  -55, -277,  -20, -289, -113 },
        { -396,  -35, -398,    3, -394,   84, -142, -169, -168, -143, -278,  -56, -289,  -21, -114 },
        {  -35, -396, -398,    3, -169, -142,   84, -394, -168, -143,  -34, -300, -307,   -2, -114 },
        { -396,  -35,    3, -398, -142, -169, -394,   84, -168, -143, -299,  -33,   -2, -307, -113 },
        {  -34,  -34, -398, -397, -168, -168, -168, -168,   59, -394,  -46,  -46, -283, -284, -117 },
        { -397, -397,    2,    2, -143, -143, -143, -143, -394,  108, -287, -287,  -26,  -25, -112 },
        {   60, -309, -324, -288,  -55, -278,  -34, -299,  -46, -287,   18, -233, -250, -202, -120 },
        { -309,   59, -287, -326, -277,  -56, -300,  -33,  -46, -287, -233,   19, -202, -252, -121 },
        { -319, -288,  105, -329,  -20, -289, -307,   -2, -283,  -26, -250, -202,   56, -258, -118 },
        { -288, -320, -327,  104, -289,  -21,   -2, -307, -284,  -25, -202, -252, -258,   56, -119 },
        { -167, -168, -163, -165, -113, -114, -114, -113, -117, -112, -120, -121, -118, -119,  -74 },
	};

int SM2[15][15] =
	{
        {   73, -260, -282, -250,  -48, -221,  -38, -236,  -43, -252,   13, -148, -173, -150,  -77 },
        { -260,   73, -250, -283, -220,  -48, -235,  -37,  -43, -252, -147,   13, -149, -174,  -76 },
        { -282, -250,  115, -293,  -14, -241, -254,   -5, -204,  -16, -168, -148,   49, -188, -109 },
        { -250, -283, -293,  115, -242,  -14,   -5, -257, -204,  -16, -148, -170, -188,   49, -110 },
        {  -48, -220,  -14, -242,   30, -196, -104,  -68,  -77,  -85,  -26, -107,   -2, -132,  -46 },
        { -221,  -48, -241,  -14, -196,   30,  -68, -103,  -77,  -84, -107,  -26, -131,   -2,  -46 },
        {  -38, -235, -254,   -5, -144,  -68,   33, -209,  -80,  -87,  -19, -119, -144,    5,  -48 },
        { -236,  -37,   -5, -257,  -68, -103, -209,   33,  -80,  -87, -119,  -18,    5, -145,  -48 },
        {  -43,  -43, -204, -204,  -77,  -77,  -80,  -80,    8, -186,  -13,  -13, -102, -103,  -20 },
        { -252, -252,  -16,  -16,  -85,  -84,  -87,  -87, -186,   56, -134, -133,  -11,  -11,  -75 },
        {   13, -147, -168, -148,  -26, -107,  -19, -119,  -13, -134,   16,  -61,  -81,  -66,  -18 },
        { -148,   13, -148, -170, -107,  -26, -119,  -18,  -13, -133,  -61,   16,  -66,  -82,  -18 },
        { -173, -149,   49, -188,   -2, -131, -144,    5, -102,  -11,  -81,  -66,   38,  -98,  -39 },
        { -150, -174, -188,   49, -132,   -2,    5, -145, -103,  -11,  -66,  -82,  -98,   39,  -39 },
        {  -77,  -76, -109, -110,  -46,  -46,  -48,  -48,  -20,  -75,  -18,  -18,  -39,  -39,   -4 },
	};

int    Table_size = 101;
static double ***C_LN_C;
static double ***C_LN_N;
static double  **C_LN_P;

double ALLR_score( int *, int * );
double ALLR_function( int *, double *, int, int *, double *, int );
double calculate_ALLR_double( double *, int, double *, int );
int    map_column_dcode( double * );
void   init_allr_lookup_table();


/* ALLR function */
double ALLR_function( int *col_1, double *f_1, int num_seq_1, 
		      int *col_2, double *f_2, int num_seq_2 )
{
  int    i;
  double score     = 0; 
  //double s1, s2;
  double calculate_ALLR_double();
  
  /***********************************************************
   * ALLR function elements:
    
     ALLR = [ n1*f1*log(f2/p) + n2*f2*log(f1/p) ]/[ n1+n2 ]
          = [ c1*log(c2/(n2*p)) + c2*log(c1/(n1*p)) ]/[ n1+n2 ]
          = [ c1*log(c2) - c1*log(n2) - c1*log(p) 
             +c2*log(c1) - c2*log(n1) - c2*log(p) ]/[ n1+n2 ]
     
     look up:
     	C_LN_C, [A_size][count_1][count_2]
     	C_LN_N, [A_size][count_1][count_2]
     	C_LN_P; [A_size][count_1]
   *
   ***********************************************************/

  if ( 1 ) //num_seq_1 >= Table_size || num_seq_2 >= Table_size )
  {
    score = calculate_ALLR_double( f_1, num_seq_1, f_2, num_seq_2 );
    //score = ALLR_score( col_1, col_2 );
    return score;
  }
  else
  {
    for ( i=0; i<A_size; ++i )
    {
      score += C_LN_C[i][ col_2[i] ][ col_1[i] ] +
               C_LN_C[i][ col_1[i] ][ col_2[i] ] -
	       C_LN_N[i][ col_2[i] ][ num_seq_1+1 ] -
	       C_LN_N[i][ col_1[i] ][ num_seq_2+1 ] -
	       C_LN_P[i][ col_1[i] ] -
	       C_LN_P[i][ col_2[i] ];
    }
    score /= num_seq_1 + num_seq_2 + 2;
    
    //s1 = score;
    //s2 = calculate_ALLR_double( f_1, num_seq_1, f_2, num_seq_2 );
    //printf ( "%d:%d\t%.3f\t%.3f\n", num_seq_1, num_seq_2, s1, s2 );
  }
  return ( score );
}


/* ALLR function */
double calculate_ALLR_double( double *f_1, int num_seq_1,
		              double *f_2, int num_seq_2 )
{
  double score;
  double I_1;
  double I_2;

  double cal_I();

  score = 0;
  I_1 = cal_I( f_1, f_2 );
  I_2 = cal_I( f_2, f_1 );
  score = ( I_1*num_seq_2+I_2*num_seq_1 )/( num_seq_1+num_seq_2 );
  
  /*score = ( I_1*(num_seq_2+1)+I_2*(num_seq_1+1) )/( num_seq_1+num_seq_2+2 );*/
   /*printf ( "%d-%d: %.2f  %.2f  %.3f\n", num_seq_1, num_seq_2, I_1, I_2, score );
   
   printf ( "%.2f-%.2f-%.2f-%.2f :  %.2f-%.2f-%.2f-%.2f %.3f\n",
            f_1[0], f_1[1], f_1[2], f_1[3], f_2[0], f_2[1], f_2[2], f_2[3], score );
   */
  return ( score );
}

/* Calculate I' */
double cal_I( double *f_1, double *f_2 )
{
  int       i;
  double    I = 0;
  
  for ( i=0; i<A_size; ++i )
  {
    I += f_2[i] * log(f_1[i]/P[i])/LN2;
  }
  return ( I );
}



/* Function to initialize ALLR lookup table */
void init_allr_lookup_table()
{
  int i, j, k;
  int D1 = A_size;
  int D2 = Table_size;
  int D3 = Table_size;
  
  /***********************************************************
   * ALLR function elements:
    
     ALLR = [ n1*f1*log(f2/p) + n2*f2*log(f1/p) ]/[ n1+n2 ]
          = [ c1*log(c2/(n2*p)) + c2*log(c1/(n1*p)) ]/[ n1+n2 ]
          = [ c1*log(c2) - c1*log(n2) - c1*log(p) 
             +c2*log(c1) - c2*log(n1) - c2*log(p) ]/[ n1+n2 ]
     
     look up:
     	C_LN_C, [A_size][count_1][count_2]
     	C_LN_N, [A_size][count_1][count_2]
     	C_LN_P; [A_size][count_1]
   *
   ***********************************************************/
  
  C_LN_C = (double ***)tw_calloc( D1, 
    				  sizeof(double **), 
				  "C_LN_C",
				  "init_allr_lookup_table()" );
  for ( i=0; i<D1; ++i )
  {
    C_LN_C[i] = (double **)tw_calloc( D2,
    				      sizeof(double *),
				      "C_LN_C[]",
				      "init_allr_lookup_table()" );
    for ( j=0; j<D2; ++j )
    {
      C_LN_C[i][j] = (double *)tw_calloc( D3,
      					  sizeof(double),
					  "C_LN_C[][]",
					  "init_allr_lookup_table()" );
    }
  }
  
  C_LN_N = (double ***)tw_calloc( A_size, 
    				  sizeof(double **), 
				  "C_LN_N",
				  "init_allr_lookup_table()" );  
  for ( i=0; i<D1; ++i )
  {
    C_LN_N[i] = (double **)tw_calloc( D2,
    				      sizeof(double *),
				      "C_LN_N[]",
				      "init_allr_lookup_table()" );
    for ( j=0; j<D2; ++j )
    {
      C_LN_N[i][j] = (double *)tw_calloc( D3,
      					  sizeof(double),
					  "C_LN_N[][]",
					  "init_allr_lookup_table()" );
    }
  }
				       
  C_LN_P = (double **)tw_calloc( A_size, 
    				 sizeof(double *), 
				 "C_LN_P",
				 "init_allr_lookup_table()" );
  for ( i=0; i<D1; ++i )
  {
    C_LN_P[i] = (double *)tw_calloc( D2,
    				     sizeof(double),
				     "C_LN_P[]",
				     "init_allr_lookup_table()" );
  }
  
  for ( i=0; i<D1; ++i )
  {
    for ( j=0; j<D2; ++j )
    {
      for ( k=0; k<D3; ++k )
      {
        
        C_LN_C[i][j][k] = ( j /*+ PSEUDOCNT*P[i]*/ ) * log( k + P[i]*PSEUDOCNT ) / LN2;
	C_LN_N[i][j][k] = j * log( k + PSEUDOCNT ) / LN2;
      }
      C_LN_P[i][j] = (j /*+ PSEUDOCNT*P[i]*/) * log( P[i] ) / LN2;
    }
  }

}


/* Map a frequency vector to its degenerate dcode */
int map_column_dcode( double *base )
{
  int i;
  int top1;
  int top2;
  int d = 14;
  
  top1 = 0;
  for ( i=1; i<A_size; ++i )
  {
    if ( base[i]>base[top1] )
    {
      top1 = i;
    }
  }
  top2 = ( top1==0 )? 1 : 0;
  for ( i=0; i<A_size; ++i )
  {
    if ( i!=top1 )
    {
      if ( base[i]>base[top2] )
      {
        top2 = i;
      }
    }
  }
  
  if ( base[top1] >= 0.7 )
  {
    return top1;
  }
  
  if ( base[top1]>=0.5 )
  {
    if ( base[top2]<0.35 )
    {
      if ( top1==0 )
      {
      	d=10; // a
      }
      else if ( top1==1 )
      {
      	d=11; // t
      }
      else if ( top1==2 )
      {
      	d=12; // c
      }
      else 
      {
      	d=13; // g
      }
      return d;
    }
    else
    {
      if ( ( top1==0 && top2==1 ) || (top1==1 && top2==0) )
      {
      	d=8; // c='W';
      }
      else if ( ( top1==0 && top2==2 ) || (top1==2 && top2==0) )
      {
      	d=4; // c='M';
      }
      else if ( ( top1==0 && top2==3 ) || (top1==3 && top2==0) )
      {
      	d=6; // c='R';
      }
      else if ( ( top1==1 && top2==2 ) || (top1==2 && top2==1) )
      {
      	d=7; // c='Y';
      }
      else if ( ( top1==1 && top2==3 ) || (top1==3 && top2==1) )
      {
      	d=5; // c='K';
      }
      else if ( ( top1==2 && top2==3 ) || (top1==3 && top2==2) )
      {
      	d=9; // c='S';
      }
      return d;
    }
  }
  
  if ( base[top1]>=0.35 && base[top2]>=0.35 )
  {
    if ( ( top1==0 && top2==1 ) || (top1==1 && top2==0) )
    {
      	d=8; // c='W';
    }
    else if ( ( top1==0 && top2==2 ) || (top1==2 && top2==0) )
    {
      	d=4; // c='M';
    }
    else if ( ( top1==0 && top2==3 ) || (top1==3 && top2==0) )
    {
      	d=6; // c='R';
    }
    else if ( ( top1==1 && top2==2 ) || (top1==2 && top2==1) )
    {
      	d=7; // c='Y';
    }
    else if ( ( top1==1 && top2==3 ) || (top1==3 && top2==1) )
    {
      	d=5; // c='K';
    }
    else if ( ( top1==2 && top2==3 ) || (top1==3 && top2==2) )
    {
      	d=9; // c='S';
    }
    return d;
  }
  
  return 14;
} 


/* ALLR function */
double ALLR_score( int *col_1, int *col_2 )
{
  int     i;
  double  score = 0;
  double *f_1;
  double *f_2;  
  double  I_1;
  double  I_2;
  int     num_seq_1=0;
  int     num_seq_2=0;
  
  double   *cal_f();
  double    cal_I();
  
  f_1 = cal_f( col_1 );
  f_2 = cal_f( col_2 );
  
  I_1 = cal_I( f_1, f_2 );
  I_2 = cal_I( f_2, f_1 );
  
  for ( i=0; i<A_size; ++i )
  {
    num_seq_1 += col_1[i];
    num_seq_2 += col_2[i];
  }
  
  score = ( I_1*num_seq_2+I_2*num_seq_1 )/( num_seq_1+num_seq_2 );
  //score = ( I_1*num_seq_2+I_2*num_seq_1 );
  /*score = ( I_1*(num_seq_2+1)+I_2*(num_seq_1+1) )/( num_seq_1+num_seq_2+2 );*/
  
  tw_free( f_1, "f_1", "ALLR_score()" );
  tw_free( f_2, "f_2", "ALLR_score()" );
  
  return ( score );
}


/* Calcualte frequency from count vector */
double * cal_f( int *col )
{
  int     i;
  double *f;
  int     num_seq = 0;
  double  pseudo = 1.0; /* --> PSEUDOCOUNT */
  
  f = (double *)tw_calloc( A_size, sizeof(double), "f", "cal_f()" );
  
  for ( i=0; i<A_size; ++i )
  {
    num_seq += col[i];
  }
  
  /*pseudo = sqrt( (double)num_seq );*/
  for ( i=0; i<A_size; ++i )
  {
    f[i] = ( col[i] + P[i]*(double)pseudo ) / (double)(num_seq+pseudo) ;
  }
  
  return ( f );
}





