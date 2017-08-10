/****************************************************************
 * Ting Wang's C Utility Collection				*
 * Copyright (C) 2003, 2004, 2005                               *
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

/* Dynamic Memory Allocation
 * Based on Gerald Hertz(Copyright1990-2001)'s consensus package
 */


#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include "tw_alloc.h"


/* The following functions allocate dynamic memory with markers preceding and
 * following the memory to make sure the memory has not been written over
 * when freed or reallocated.  The int containing ISBEGIN is followed by
 * an int indicating the size of the memory allocated.  The concept behind
 * these functions is motivated by a GNU debugging "malloc()". */
 
#define ISBEGIN (11111)    /* Magic int preceding dynamic memory. */
#define ISEND (123)        /* '{': Magic char following dynamic memory. */
#define ISFREED (-11111)   /* Magic int replacing ISBEGIN to indicate
			    * that the memory has been freed. */


static void *malloc_check( char   *alloc_function,
			   size_t  size,
			   char   *variable,
			   char   *function );
static void *realloc_check( char   *alloc_function,
			    void   *array,
			    size_t  size,
			    char   *variable,
			    char   *function );
static size_t determine_gross_size( size_t size );
static void *determine_gross_array( void *array, size_t size );


/* Malloc with an error message and an exit if unsuccessful. 
 */
void *tw_malloc(
      size_t size,      /* Size of the array. */
      char *variable,   /* Name of the variable. */
      char *function    /* Name of the calling function. */
      )
{
  void *array;
  array = malloc_check( "tw_malloc",
  			size,
  			variable,
  			function );
  return( array );
}


/* Calloc (actually "malloc") with an error message and exit if unsuccessful.
 */
void *tw_calloc(
      int elt_count,     /* Number of elements in the array. */
      size_t elt_size,   /* Size of each element. */
      char *variable,    /* Name of the variable. */
      char *function     /* Name of the calling function. */
      )
{
  void *array;
  array = malloc_check( "tw_calloc",
  			elt_count * elt_size,
  			variable,
  			function );
  return( array );
}


/* Realloc with an error message and an exit if unsuccessful.  
 */
void *tw_realloc(
      void *array,     /* The array whose size is to be modified. */
      size_t size,     /* Size of the array. */
      char *variable,  /* Name of the variable. */
      char *function   /* Name of the calling function. */
      )
{
  if ( array == (void *)NULL )
  {
    return( malloc_check( "tw_realloc", 
    			  size, 
    			  variable, 
    			  function ) );
  }
  else
  {
    return( realloc_check( "tw_realloc", 
    			   array, 
    			   size, 
    			   variable, 
    			   function ) );
  }
}


/* Realloc an array with an error message and an exit if unsuccessful.
 */
void *tw_recalloc(
      void *array,      /* The array whose size is to be modified. */
      int elt_count,    /* Number of elements in the array. */
      size_t elt_size,  /* Size of each element. */
      char *variable,   /* Name of the variable. */
      char *function    /* Name of the calling function. */
      )
{

  if ( array == (void *)NULL )
  {
    return( malloc_check( "tw_realloc", 
    			  elt_count * elt_size, 
    			  variable, 
    			  function ) );
  }
  else
  {
    return( realloc_check( "tw_realloc", 
    			   array, 
    			   elt_count * elt_size, 
    			   variable, 
    			   function ) );
  }
}


/* Free dynamic memory, but do nothing if handed a NULL pointer. */
void tw_free(
      void *array,      /* The array being freed. */
      char *variable,   /* Name of the variable. */
      char *function)   /* Name of the calling function. */
{
  void *gross_array;     /* Array, including the magic ints. */

  if ( array != (void *)NULL )
  {
    gross_array = memory_integrity_check( "free", 
    					  array, 
    					  variable, 
    					  function );
    ((int *)gross_array)[0] = ISFREED;
    free( gross_array );
  }
  else
  {
    fprintf( stderr, "Warning: " );
    fprintf( stderr, "\"%s\" is NULL -- address %p -- in ", variable, array );
    fprintf( stderr, "the \"%s\" function.\n", function );
  }
}


/* Determine the amount of space to allocate so that there is room for
 * the magic bytes. */
static size_t determine_gross_size( size_t size )
 			/* The number of bytes originally called for. */
{
  return( size + 2 * sizeof(int) + sizeof(char) );
}


/* Set the magic ints and char.  Return a pointer to the
 * beginning of the usable memory. */
static void *determine_gross_array(
     	     void *array,  /* The allocated memory, 
     	     		    * including the magic ints and char. */
     	     size_t size   /* The amount of usable memory. */
     	     )
{
  int  *int_array  = (int *)array;
  char *char_array = (char *)(int_array + 2);

  int_array[0] = ISBEGIN;
  int_array[1] = size;
  char_array[size] = ISEND;

  return( (void *)char_array );
}


/* Malloc with an error message tuned to the type of alloc,
 * and an exit if unsuccessful.  
 * Returns a pointer to the usable memory. */
static void *malloc_check(
     	     char *alloc_function,  /* Name of allocating function*/
     	     size_t size,           /* Size of the array. */
     	     char *variable,        /* Name of the variable. */
     	     char *function         /* Name of calling function. */
     	     )
{
  void *gross_array;  /* Array, including the magic ints and char. */
  size_t gross_size;  /* The gross size of the memory allocated. */

  gross_size  = determine_gross_size( size );
  gross_array = malloc( gross_size );
  if ( gross_array == (void *)NULL )
  {
    fprintf( stderr, "Cannot %s space for \"%s\" ", alloc_function, variable );
    fprintf( stderr, "in the \"%s\" function.\n", function );
    exit( 1 );
  }
  else 
  {
    return( determine_gross_array( gross_array, size ) );
  }
}


/* Realloc with an error message tuned to the type of alloc,
 * and an exit if unsuccessful. */
static void *realloc_check(
     	     char *alloc_function,  /*Name of allocating function*/
     	     void *array,           /* Array whose size is to be modified.*/
     	     size_t size,           /* Size of the array. */
     	     char *variable,        /* Name of the variable. */
     	     char *function         /* Name of calling function. */
     	     )
{
  void *gross_array;   /* Array, including the magic ints and char. */
  size_t gross_size;   /* The gross size of the memory allocated. */

  gross_array = memory_integrity_check( alloc_function, 
  					array, 
  					variable, 
  					function );
  gross_size  = determine_gross_size( size );
  gross_array = realloc( gross_array, gross_size );

  if ( gross_array == (void *)NULL )
  {
    fprintf( stderr, "Cannot %s space for \"%s\" ", alloc_function, variable );
    fprintf( stderr, "in the \"%s\" function.\n", function );
    exit( 1 );
  }
  else
  {
    return( determine_gross_array( gross_array, size ) );
  }
}


/* Check the array to make sure the magic ints and char are all consistent
 * with unmolested memory.  Returns a void pointer to the gross memory.
 * When called outside of "tw_alloc_helper.c",  alloc_function = "check". */
void *memory_integrity_check(
      char *alloc_function,  /* Name of allocating function.*/
      void *array,           /* Array being checked. */
      char *variable,        /* Name of the variable. */
      char *function         /* Name of the calling function. */
      )
{
  char error = 0;           /* Flag indicating whether an error was detected.
			     * 0: no error; 1: an error was detected. */
  int *int_array = ((int *)array) - 2;

  if ( int_array[0] == ISFREED )
  {
    fprintf( stderr, "The memory being %s", alloc_function );
    fprintf( stderr, "ed for\nthe \"%s\" variable in the\n", variable );
    fprintf( stderr, "\"%s\" function has already been freed.\n\n", function );
    error = 1;
  }

  else if ( int_array[0] != ISBEGIN )
  {
    fprintf( stderr, "The beginning of the memory being %s", alloc_function );
    fprintf( stderr, "ed for\nthe \"%s\" variable in the\n", variable );
    fprintf( stderr, "\"%s\" function has been overwritten.\n\n", function );
    error = 1;
  }

  if ( ((char *)array)[int_array[1]] != ISEND)
  {
    fprintf( stderr, "The end of the memory being %s", alloc_function );
    fprintf( stderr, "ed for\nthe \"%s\" variable in the\n", variable );
    fprintf( stderr, "\"%s\" function has been overwritten.\n\n", function );
    error = 1;
  }

  if ( error == 1 ) 
  {
    exit(1);
  }
  else 
  {
    return( (void *)int_array );
  }
}
