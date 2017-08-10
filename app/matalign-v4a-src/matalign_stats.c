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
 /* Implementation of Karlin Altschul statistics is based on 
  * Warren Gish and Sean Eddy, personal communication. 
  */

#include "matalign_options.h"

#define MAXIT     200
#define SUMLIMIT  0.00001
#define TOLERANCE 1e-6      /* how close to 0 root-finders must get */
#define ZERO 5.0E-324

static double newtonraphson( int **, double *f, int );
static double val  (int **s, double *f, int K, double lambda);
static double deriv(int **s, double *f, int K, double lambda);

void karlin_altschul();
void find_score_spectrum();
//void find_query_frequency();
void karlin_lambda_K_H();
int  greatest_common_denominator( int, int );

double get_ALLR_E_value( double );
double get_hsp_E_value ( HSP * );
double get_hsp_P_valie ( HSP * );
double e_value_to_p_value( double );
double log_poisson_unit( double, int );


void karlin_altschul()
{
  int     i, j;
  int   **sm;
  double *sm_f;
  double *query_f;
  int     dim;
  int     min;
  int     max;
  double *score_spectrum;
  
  void find_score_spectrum( double **, int **, double *, double *, int, int, int );
  //void find_query_frequency( double **, int );
  void karlin_lambda_K_H( int, int, double *, double *, double *, double * );
  
  /* Point to the scoring matrix */
  dim = 15;
  sm = (int **)tw_calloc( dim, sizeof(int *), "sm", "karlin_altschul()" );
  for ( i=0; i<dim; ++i )
  {
    sm[i] = (int *)tw_calloc( dim, sizeof(int), "sm[]", "karlin_altschul()" );
    for ( j=0; j<dim; ++j )
    {
      sm[i][j] = SM1[i][j];
    }
  }
  
  sm_f = Dcode_f;
  //find_query_frequency( &query_f, dim );
  
  /* Found the highest and lowest value in scoring matrix */
  min = max = sm[0][0];
  for ( i=0; i<dim; ++i )
  {
    for ( j=0; j<dim; ++j )
    {
      min = min < sm[i][j] ? min : sm[i][j];
      max = max > sm[i][j] ? max : sm[i][j];
    }
  }
  
  /* Found the range of scores and their probabilities:
   * Note: second sm_f is query_f in WRG version 
   * The sm_f and query_f should result in using a smaller lambda,
   * based on conversation with Warren Gish.
   */
  find_score_spectrum( &score_spectrum, sm, sm_f, sm_f, min, max, dim );
  
  /* Calculate lambda, K and H */
  karlin_lambda_K_H( min, max, score_spectrum, &Lambda, &K, &H );
  
  /* Calculate lambda */
  if ( 0 )
  {
    Lambda = newtonraphson( sm, sm_f, dim );
  }
  
  //printf ( "\n" );
  //printf ( "Parameters for Karlin-Altschul statistics:\n" );
  //printf ( "\tLambda: %.5f\n", Lambda );
  //printf ( "\tK: %.5f\n", K );
  //printf ( "\tH: %.5f bits\n", H );
  
}



void karlin_lambda_K_H( int low, 
			int high,
			double *prob,
			double *lambda_r,
			double *K_r,
			double *H_r )
{
  int range;
  int i, j;
  int lo, hi, first, last;
  double sum;
  double up, new, Sum, av;
  double *p0 = NULL, *P0 = NULL;
  double *p, *P, *ptrP, *ptr1, *ptr2, *ptr1e;
  
  int greatest_common_denominator( int, int );
  
  
  if ( low >= 0 )
  {
    printf ( "Error: The lowest score in the substitution matrix must be negative.\n" );
  }
  if ( high <= 0 )
  {
    printf ( "Error: The highest score in the substitution matrix must be positive.\n" );
  }
  
  range = high - low;
  for ( i=range; i>-low && prob[i]==0; --i )
  {
    ;
  }
  if ( i<=-low )
  {
    printf ( "Error: A positive score is impossible with the given substitution matrix " );
    printf ( "       and residue composition.\n" );
  }
  
  for ( sum=0, i=0; i<=range; sum+=prob[i++] )
  {
    if ( prob[i]<0 )
    {
      printf ( "Error: Negative probabilities are disallowed.\n" );
    }
  }
  
  p0 = p = (double *)tw_malloc( sizeof(*p)*(range+1), 
  				"p0", 
  				"karlin_lambda_K_H()" );
  P0 = P = (double *)tw_malloc( MAXIT*sizeof(*P)*(range+1),
  				"P0",
  				"karlin_lambda_K_H()" );
  
  for ( i=0; i<=range; ++i )
  {
    prob[i] /= sum;
  }
  for ( sum=0, i=0; i<=range; sum+=prob[i++] )
  {
    if ( prob[i]<0 )
    {
      printf ( "Error: Negative probabilities are disallowed.\n" );
    }
  }
  
  if ( sum<0.99995 || sum>1.00005 )
  {
    printf ( "Error: Probabilities sum to %.5f and will be normalized to 1.\n",
    	     sum );
  }
  for ( Sum=low, i=0; i<=range; ++i )
  {
    Sum += i * (p[i] = prob[i]/sum);
  }
  
  if ( Sum>=0 )
  {
    printf ( "Error: Invalid (non-negative) expected score: %#0.21g", Sum );
  }
  
  /* Calculate the parameter Lambda */
  up = 0.5;
  do
  {
    up *= 2;
    ptr1 = p;
    for( sum=0, i=low; i<=high; ++i )
    {
      sum+= *ptr1++ * exp(up*i);
    }
  }
  while ( sum<1.0 );
  
  for ( *lambda_r=0, j=0; j<15; ++j )
  {
    new = ( *lambda_r+up )/2.0;
    ptr1 = p;
    for ( sum=0, i=low; i<=high; ++i )
    {
      sum += *ptr1++ * exp(new*i);
    }
    if ( sum>1.0 )
    {
      up = new;
    }
    else
    {
      *lambda_r = new;
    }
  }
  
  /* Calculate the relative entropy of the p's and q's, the parameter H */
  ptr1 = p;
  for ( av=0, i=low; i<=high; ++i )
  {
    av += *ptr1++ *i*exp( *lambda_r*i );
  }
  *H_r = *lambda_r * av/LN2;
  
  /* Calculate the parameter K */
  if ( low == -1 || high == 1 )
  {
    *K_r = ( high==1 ? av : Sum*Sum/av );
    *K_r *= 1.0 - exp( - *lambda_r );
    goto OKExit;
  }
  Sum = 0;
  lo = hi = 0;
  for ( *P=sum=1, j=1; j<=MAXIT && sum>SUMLIMIT; Sum+=sum/=j++ )
  {
    first = last = range;
    for ( ptrP=P+(hi+=high)-(lo+=low); ptrP>=P; *ptrP-- = sum )
    {
      ptr1 = ptrP-first;
      ptr1e = ptrP-last;
      ptr2 = p+first;
      for ( sum=0; ptr1>=ptr1e; )
      {
        sum += *ptr1-- * *ptr2++;
      }
      if ( first )
      {
      	--first;
      }
      if ( ptrP-P <=range )
      {
      	--last;
      }
    }
    for ( sum=0, i=lo; i; ++i )
    {
      sum += *++ptrP * exp(*lambda_r*i);
    }
    for ( ; i<=hi; ++i )
    {
      sum += *++ptrP;
    }
  }
  
  if ( j>MAXIT )
  {
    printf ( "Error: Value for K may be too large due to insufficient iteration.\n" );
  }
  
  for ( i=low; p[i-low]==0; ++i )
  {
    ;
  }
  for ( j=-i; i<high&&j>1; )
  {
    if ( p[++i-low] )
    {
      j = greatest_common_denominator( j, i );
    }
  }
  
  *K_r = (j*exp(-2.0*Sum))/(av*(1.0-exp(- *lambda_r*j )));
  
OKExit:
  if ( p0 != NULL )
  {
    tw_free( p0, "p0", "karlin_lambda_K_H()" );
  }
  if ( P0 != NULL )
  {
    tw_free( P0, "P0", "karlin_lambda_K_H()" );
  }
  
}

  
/* Calcualte possible score frequency */
void find_score_spectrum( double **ret_spectrum, 
			  int    **sm,
			  double  *sm_f,
			  double  *query_f,
			  int      min,
			  int      max,
			  int      dim )
{
  int     i, j;
  int     range;
  double *spectrum;
  
  range = max - min;
  
  spectrum = (double *)tw_calloc( range+1,
  				  sizeof(double),
  				  "spectrum",
  				  "find_score_spectrum()" );
  for ( i=0; i<=range; ++i )
  {
    spectrum[i] = 0.;
  }
  
  for ( i=0; i<dim; ++i )
  {
    for ( j=0; j<dim; ++j )
    {
      spectrum[ sm[i][j]-min ] += query_f[i]*sm_f[j];
    }
  }
  
  *ret_spectrum = spectrum;
}
  

/* Calculate the frequency of dcode in query sequence */
/*void find_query_frequency( double **ret_f, int dim )
{
  int     i, j;
  double *query_f;
  
  query_f = (double *)tw_calloc( dim,
  				 sizeof(double),
  				 "query_f",
  				 "find_query_frequency()" );
  for ( i=0; i<dim; ++i )
  {
    query_f[i] = 0;
  }
  
  if ( EZ )
  {
    for ( i=0; i<Query_sequence_len; ++i )
    {
      if ( Core_promoter->base_seq[i]->num_ele > 0 )
      {
        query_f[ Core_promoter->base_seq[i]->d_code ]++;
      }
      if ( Core_promoter_comp->base_seq[i]->num_ele > 0 )
      {
        query_f[ Core_promoter_comp->base_seq[i]->d_code ]++;
      }
    }
    for ( i=0; i<dim; ++i )
    {
      query_f[i] /= 2*Query_effective_len;
    }
  }
  else
  {
    for ( i=0; i<Core_profile_num; ++i )
    {
      for ( j=0; j<Core_profile[i]->width; ++j )
      {
        query_f[ (Core_profile[i]->dcode)[j] ]++;
        query_f[ (Core_profile_comp[i]->dcode)[j] ]++;
      }
    }
    for ( i=0; i<dim; ++i )
    {
      query_f[i] /= 2*Query_profile_len;
    }
  }
  
  *ret_f = query_f;

}*/


static double newtonraphson( int **s, double *f, int K )
{
  double lambda, newlambda;
  int    iter;
  double fx, dfx;
  
  /* There's two zeros, not one: lambda=0 is a solution.
   * To find the nonzero zero (ahem), we've got to make sure 
   * we start to the right of it: where the function is positive.
   * Start with anything, then move right until fx > 0.
   */
  lambda = 0.1;
  for ( iter = 0; iter < 100; iter++ ) 
  {
    if ( val(s, f, K, lambda) > 0. )
    {
      break;
    }
    lambda *= 2.;
  }
  
  if ( iter == 100 )
  { 
    fprintf(stderr, "failed to find a start pt for newton/raphson\n");
    exit(1);
  }

  /* Now, starting from there, come back in with the Newton/Raphson
   * algorithm. At least in theory, because we know what f(x) looks
   * like, this should be well-behaved, smoothly converging to the
   * solution.
   */
  for ( iter = 0; iter < 100; iter++ ) 
  {
    fx = val( s, f, K, lambda );
    if ( fabs(fx) < TOLERANCE )
    {
      break;   /* success! */
    }
    dfx = deriv( s, f, K, lambda );
    newlambda = lambda - fx / dfx; /* This update defines Newton/Raphson */
    if ( newlambda <= 0 )
    {
      newlambda = 0.000001; /* this shouldn't happen */
    }
    lambda = newlambda;
  }
  if ( iter == 100 ) 
  { 
    fprintf(stderr, "Newton/Raphson search failed\n"); 
    exit(1); 
  }
  
  return ( lambda );
}


/* Calculate and return f(lambda):  \sum_ab f_a f_b e^{lambda s_ab} - 1
 */                                     
static double val( int **s, double *f, int K, double lambda )
{
  int   a,b;
  double total = -1.;

  for (a = 0; a < K; a++)
    for (b = 0; b < K; b++)
      total += f[a] * f[b] * exp(lambda * s[a][b]);
  return total;
}


/* First derivative is:  \sum_ab f_a f_b s_ab e^{lambda s_ab}
 */
static double deriv( int **s, double *f, int K, double lambda )
{
  int   a,b;
  double deriv = 0.;

  for (a = 0; a < K; a++)
    for (b = 0; b < K; b++)
      deriv += f[a] * f[b] * s[a][b] * exp(lambda * s[a][b]);
  return deriv;
}

double get_hsp_E_value( HSP *hsp )
{
  double e;
  
  e = K * (hsp->profile_1->width) * (hsp->profile_2->width) * 
  		exp( -Lambda*(hsp->hsp_score*100) );
  
  return ( e );
}

double get_ALLR_E_value( double allr )
{
  double e;
  
  e = K*M*N*exp(-Lambda*allr);
  
  //printf ( "%d, %d\n", profile[0]->width, profile[1]->width );
  //printf ( "K=%.4f\n", K );
  //printf ( "M=%.2f, N=%.2f\n", M, N );
  //printf ( "Lambda=%.4f\n", Lambda );
  //printf ( "allr=%.4f\n", allr );
  //printf ( "E=%G\n", e );
  
  //printf ( "K=%.4f, M=%d, N=%d, Lambda=%.4f, allr=%.4f, e=%.4f\n", K, M, N, Lambda, allr, e );
  
  return ( e );
}


double e_value_to_p_value( double e )
{
  double p;
  
  /* P( S >= x )   <=   1 - exp [ - KN exp ( - lambda * x ) ] */
  
  p = 1 - exp( -e );
  
  return ( p );
}


double e_to_p_helper_1( double y, int m )
{
  int    i;
  double p = 0;
  
  
  double pp = 0;
  double log_pp = 0;
  
  log_pp = m*log(y) - y;
  for ( i=1; i<=m; ++i )
  {
    log_pp -= log(i);
  }
  p = pp = exp(log_pp);
  
  for ( i=m+1; i<2*m+100; ++i )
  {
    log_pp += log( y );
    log_pp -= log( i );
    pp = exp( log_pp );
    p += pp;
    if ( pp <= ZERO )
    {
      break;
    }
  }
  
  
  /*
  for ( i=m; i<2*m+100; ++i )
  {
    pp = exp( log_poisson_unit( e, i ) );
    p += pp;
    if ( pp<=ZERO )
    {
      break;
    }
  }*/
  //printf ( "use helper 1, e=%.2f, m=%d\n", e, m );
  return ( p ); 
}

double e_to_p_helper_2( double e, int m )
{
  int    i;
  double p, pp;
  
  p  = 1;
  pp = 1;
  for ( i=1; i<m; i++ )
  {
    pp *= e;
    pp /= i;
    p += pp;
  }
  
  p *= exp(-e);
  p = 1-p;
  
  //printf ( "use helper 2, e=%.2f, m=%d\n", e, m );
  return ( p );
}


double log_poisson_unit( double y, int m )
{
  double u;
  int    i;
  
  u = m*log(y) - y;
  for ( i=1; i<=m; i++ )
  {
    u -= log(i);
  }
  return ( u );
}


/* Find the greatest common denominator of 2 integers */
int greatest_common_denominator( int a, int b )
{
  if ( b==0 )
  {
    return ( a );
  }
  else
  {
    return greatest_common_denominator( b, a%b );
  }
}

/***************************
 * End of matalign_stats.c *
 ***************************/

