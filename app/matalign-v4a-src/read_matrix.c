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

static char STbuffer[MAXLINE];

void      read_matrices();
void      read_list_file( char *, char **, char ***, int * );
char     *read_first_string();
PROFILE **read_matrix_files( char *, char **, int * );
PROFILE  *read_matrix( char *, char * );
PROFILE  *make_revcomp_profile( PROFILE * );


/* Read matrices */
void read_matrices()
{
  int       i;
  
  void      read_list_file();
  PROFILE **read_matrix_files();
  PROFILE  *make_revcomp_profile();
  
  /* Read file containing list of matrix names */
  read_list_file( file1, &path1, &Mat_names_1, &Mat_num_1 );
  if ( List_num == 2 )
  {
    read_list_file( file2, &path2, &Mat_names_2, &Mat_num_2 );
  }
  
  /* Read matrix files */
  Profiles_1 = read_matrix_files( path1, Mat_names_1, &Mat_num_1 );
  if ( List_num == 2 )
  {
    Profiles_2 = read_matrix_files( path2, Mat_names_2, &Mat_num_2 );
    for ( i=0; i<Mat_num_2; ++i )
    {
      Profiles_2[i]->mirror = make_revcomp_profile( Profiles_2[i] );
    }
    
  }
  else
  {
    for ( i=0; i<Mat_num_1; ++i )
    {
      Profiles_1[i]->mirror = make_revcomp_profile( Profiles_1[i] );
    }
  }
}

/* Read matrices into profiles */
PROFILE **read_matrix_files( char *path, char **names, int *ret_num )
{
  int       i, j;
  int       mat_cnt;
  char     *mat_file;
  PROFILE  *profile;
  PROFILE **profiles;
  
  PROFILE *read_matrix();
  
  mat_cnt  = 0;
  mat_file = (char *)NULL;
  profile  = (PROFILE *)NULL;
  
  profiles = (PROFILE **)tw_calloc( *ret_num,
  				    sizeof(PROFILE *),
  				    "profiles",
  				    "read_matrix_files()" );
  for ( i=0; i<*ret_num; ++i )
  {
    mat_file = File_put_path( path, names[i] );
    if ( File_exists( mat_file ) )
    {
      profile = read_matrix( path, names[i] );
      if ( profile )
      {
        profiles[mat_cnt] = profile;
        mat_cnt++;
      }
    }
    else
    {
      printf ( "Warning: Matrix file %s doesn't exist.\n", mat_file );
    }
  }
  profiles = (PROFILE **)tw_recalloc( profiles,
  				      mat_cnt,
  				      sizeof(PROFILE *),
  				      "profiles",
  				      "read_matrix_files()" );
  *ret_num = mat_cnt;
  return ( profiles );
}


/* Read list of matrix names. 
 * This function opens a file and returns the first array in the 
 * first string of the line as the path, then an array that contains 
 * the first string in each line of the rest of the file, and the size of the
 * array.
 */
void read_list_file( char *list_f, char **ret_path, char ***ret_name, int *ret_num )
{
  FILE  *fp;
  char **first_strings;
  char  *string1;
  int    str_len    = 0;
  int    num        = 0;
  int    array_size = CHUNK_SIZE;
  
  char *read_1st_string();
  
  first_strings = (char **)tw_calloc( array_size, 
                		      sizeof(char *),
                		      "first_strings",
                		      "read_list_file()" );
  fp = fopen( list_f, "r" );
  if ( fp == NULL )
  {
    tw_free ( first_strings, 
              "first_strings",
              "read_list_file()" );
    fprintf ( stderr, "Cannot open file %s. Exit.\n", list_f );
    exit(0);
  }
  
  /* Read the first line, which contains directory information */
  if ( fgets( STbuffer, BIG_CHUNK_SIZE-1, fp ) )
  {
    string1 = read_1st_string( STbuffer );
    str_len = strlen( string1 );
    if ( string1[str_len-1] != '\/' )
    {
      string1    = (char *)tw_recalloc( string1,
                  			str_len+2,
                  			sizeof(char),
                  			"string1",
                  			"read_list_file()" );
      string1[str_len] = '\/';
      string1[str_len+1]   = '\0';
    } 
    *ret_path = string1;
  }
  
  /* Read the rest of the file, save the first string of each line */
  while ( fgets( STbuffer, BIG_CHUNK_SIZE-1, fp ) )
  {
    num++;
    string1 = read_1st_string( STbuffer );
    if ( num >= array_size-2 )
    {
      array_size += CHUNK_SIZE;
      first_strings = (char **)tw_recalloc( first_strings,
                    			    array_size,
                    			    sizeof(char *),
                    			    "first_strings",
                    			    "read_list_file()" );
    }
    first_strings[num-1] = string1;
  }
  fclose( fp );
  *ret_name = first_strings;
  *ret_num  = num;
}


/* Read the first string from the char array */
char *read_1st_string( char *buffer )
{
  int   i, j;
  char *string1;
  
  if ( strlen(buffer) == 0 )
  {
    return NULL;
  }
  string1 = (char *)tw_calloc( strlen(buffer)+1,
               		       sizeof(char),
               		       "string1", 
               		       "read_1st_string()" );
  for ( i=0, j=0; i<strlen(buffer); ++i )
  {
    if ( isspace( buffer[i] ) )
    {
      if ( strlen( string1 ) == 0 ) {;}
      else if ( strlen( string1 ) > 0 )
      {
        break;
      }
    }
    else
    {
      string1[j] = buffer[i];
      ++j;
    }
  }
  string1[j] = '\0';
  string1    = (char *)tw_recalloc( string1,
              			    j+1,
              			    sizeof(char),
              			    "string1",
              			    "read_1st_string()" );
  return string1;
}


/* Read matrix */
PROFILE *read_matrix( char *path, char *name )
{
  int       i, j;
  char     *mat_file;
  PROFILE  *profile;
  int     **int_matrix;
  double  **f_matrix;
  int      *dcode;
  int      *num_seq;
  int       length;
  
  int     **read_int_matrix();
  char     *File_put_path();
  int       mat_column_dcode();
  
  mat_file   = File_put_path( path, name );
  int_matrix = read_int_matrix( mat_file, &length );
  f_matrix   = (double **)tw_calloc( length,
  				     sizeof(double *),
  				     "f_matrix",
  				     "read_matrix()" );
  dcode = (int *)tw_calloc( length,
            		    sizeof(int),
            		    "dcode",
            		    "read_matrix()" );
  
  num_seq = (int *)tw_calloc( length,
  			      sizeof (int),
  			      "num_seq",
  			      "read_matrix()" );
  
  for ( i=0; i<length; ++i )
  {
    num_seq[i] = 0;
    f_matrix[i] = (double *)tw_calloc( A_size,
    				       sizeof(double),
    				       "f_matrix[]",
    				       "read_matrix()" );
    for ( j=0; j<A_size; ++j )
    {
      num_seq[i] += int_matrix[i][j];
    }
    for ( j=0; j<A_size; ++j )
    {
      f_matrix[i][j] = ( int_matrix[i][j] + P[j]*(double)PSEUDOCNT )
    		       / (double)( num_seq[i]+(double)PSEUDOCNT );
    }
    dcode[i] = map_column_dcode( f_matrix[i] );
  }
  
  /* Construct profile */
  profile = (PROFILE *)tw_calloc( 1, 
            			  sizeof(PROFILE),
          			  "profile",
          			  "read_matrix()" );
  profile->mat_name   = name;
  profile->int_matrix = int_matrix;
  profile->f_matrix   = f_matrix;
  profile->dcode      = dcode;
  profile->num_seq    = num_seq;
  profile->width      = length;
  profile->mirror     = (PROFILE *)NULL;

  return ( profile );

}


/* Read an integer matrix */
int **read_int_matrix( char *file, int *ret_length )
{
  int   **matrix;
  FILE   *fp;
  char    buffer[MAXLINE];
  int     firstline[MAXLINE/10];
  char   *s;
  int     length;
  int     i;
  int     let_idx;
  
  /* Open matrix file */
  if ( ( fp=fopen( file, "r" ) ) == NULL )
  {
    fprintf( stderr, "Error in opening file %s\n", file );
    exit (1);
  }

  /* Get the length of the matrix and first letter index */
  length = -1;    /* format: A | 10 ... */
  fgets ( buffer, MAXLINE, fp );
  s = strtok( buffer, " \t\n" );
  /* Find the letter of the first line */
  for ( i=0; i<A_size; ++i )
  {
    if ( toascii((int)s[0]) == A[i] )
    {
      let_idx = i;
      break;
    }
  }
  while ( ( s = strtok( NULL, " \t\n" ) ) != NULL )
  {
    length++;
    /* Record the first line, starting from '|'. */
    if ( length>0 )
    {
      firstline[length-1] = atoi( s );
    }
  }

  /* Allocate space for the counts */
  matrix = (int **)tw_calloc( length,
              		      sizeof(int *),
            		      "matrix",
            		      "read_int_matrix()" );
  for ( i=0; i<length; ++i )
  {
    matrix[i] = (int *)tw_calloc( A_size,
              			  sizeof(int),
          			  "matrix[]",
          			  "read_int_matrix()" );
  }
  
  /* Read one line */
  for ( i=0; i<length; ++i )
  {
    matrix[i][let_idx] = firstline[i];
  }

  while ( (fgets( buffer, MAXLINE, fp ) != NULL) && 
        ( (s = strtok( buffer, " \t\n" )) != NULL) )
  {
    for ( i=0; i<A_size; ++i )
    {
      if ( toascii((int)s[0]) == A[i] )
      {
        let_idx = i;
        break;
      }
    }

    for ( i=-1; i<length; ++i )
    {
      s = strtok( NULL, " \t\n" );
      
      if ( s!=NULL && s[0]!='|' )
      {
        matrix[i][let_idx] = atoi( s );
      }
    }
  }

  fclose (fp);
  *ret_length = length;
  return ( matrix );
    
}  


/* Make reverse complement profile */
PROFILE *make_revcomp_profile( PROFILE *profile )
{
  PROFILE  *mirror;
  int     **int_matrix;
  double  **f_matrix;
  int     **pi_matrix;
  double  **pf_matrix;
  int      *dcode;
  int       width;
  int      *num_seq;
  int       i, j;
  int       comp_l;
  int       comp_a;
  
  int map_column_dcode();
  
  pi_matrix = profile->int_matrix;
  pf_matrix = profile->f_matrix;
  width     = profile->width;
  
  mirror = (PROFILE *)tw_malloc( sizeof(PROFILE), 
           			 "mirror", 
           			 "make_revcomp_profile()" );
  int_matrix = (int **)tw_calloc( width,
  			      	  sizeof(int *),
  			      	  "int_matrix",
  			          "make_revcomp_profile()" );
  f_matrix = (double **)tw_calloc( width,
  				   sizeof(double *),
  				   "f_matrix",
  				   "make_revcomp_profile()" );
  for ( i=0; i<width; ++i )
  {
    int_matrix[i] = (int *)tw_calloc( A_size,
    				      sizeof(int),
    				      "int_matrix[]",
    				      "make_revcomp_profile()" );
    f_matrix[i] = (double *)tw_calloc( A_size,
    				       sizeof(double),
    				       "f_matrix[]",
    				       "make_revcomp_profile()" );
  } 				
  dcode  = (int *)tw_calloc( width, 
            		     sizeof(int),
            		     "dcode",
            		     "make_revcomp_profile()" );
  num_seq = (int *)tw_calloc( width,
  			      sizeof(int),
  			      "num_seq",
  			      "make_revcomp_profile()" );
  
  for ( i=0; i<width; ++i )
  {
    comp_l = width - i - 1;
    for ( j=0; j<A_size; ++j )
    {
      comp_a = A_comp[j];
      int_matrix[i][j] = pi_matrix[comp_l][comp_a];
      f_matrix[i][j]   = pf_matrix[comp_l][comp_a];
    }
    num_seq[i] = (profile->num_seq)[comp_l];
    dcode[i] = map_column_dcode( f_matrix[i] );
  }
  
  mirror->mat_name   = profile->mat_name;
  mirror->int_matrix = int_matrix;
  mirror->f_matrix   = f_matrix;
  mirror->width      = width;
  mirror->num_seq    = num_seq;
  mirror->dcode      = dcode;
  mirror->mirror     = profile;
  
  return mirror;
}

