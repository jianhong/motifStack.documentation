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
 * Based on Gerald Hertz(Copyright1990-1997)'s consensus package
 */

extern void *tw_malloc( 
	     size_t  size,     /* Size of the element. */
       	     char   *variable,  /* Name of the variable. */
       	     char   *function   /* Name of the calling function. */
       	     );


/* Calloc with an error message and an exit if unsuccessful.
 */
extern void *tw_calloc( 
	     int     elt_count, /* Number of elements in the array. */
       	     size_t  elt_size,  /* Size of each element. */
             char   *variable,  /* Name of the variable. */
       	     char   *function   /* Name of the calling function. */
             );


/* Realloc with an error message and an exit if unsuccessful. 
 */
extern void *tw_realloc( 
	     void   *array,    /* The array to be modified. */
       	     size_t  size,     /* Size of the array. */
       	     char   *variable, /* Name of the variable. */
       	     char   *function  /* Name of the calling function. */
             );
                

/* Realloc an array with an error message and an exit if unsuccessful.
 */
extern void *tw_recalloc( 
	     void   *array,     /* The array to be modified. */
       	     int     elt_count, /* Number of elements in the array. */
       	     size_t  elt_size,  /* Size of each element. */
       	     char   *variable,  /* Name of the variable. */
             char   *function   /* Name of the calling function. */
       	     );


/* Free dynamic memory, but do nothing if handed a NULL pointer. */
extern void tw_free( 
	    void *array,     /* Array to be freed */
       	    char *variable,  /* Name of the variable. */
       	    char *function   /* Name of the calling function. */
       	    );


/* NULL function.
 * Defined in "tw_alloc_helper.c"
 */
extern void *memory_integrity_check(
       	     char *alloc_function,  /* Name of the allocating function. */
       	     void *array,           /* Array being checked. */
       	     char *variable,        /* Name of the variable. */
       	     char *function         /* Name of the calling function. */
             );
