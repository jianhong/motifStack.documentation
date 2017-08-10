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
 
#include <stdio.h>
#include <stdlib.h>
#include "tw_error.h"

/* A general function for reporting program bugs. */
void bug_report( char *function )
{
  fprintf( stderr, "PROGRAM BUG: %s\n", function );
  exit( 1 );
}



