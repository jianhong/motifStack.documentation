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

/* File System Management
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tw_file.h"
#include "tw_alloc.h"

#define YES 1
#define NO 0
#define DIRSLASH '/'        /* UNIX directory paths have /foo/bar */

/* Source for the following functions:
 * File_get_path()
 * File_put_path()
 * File_exists()
 */


/* Function: File_get_path(char *file)
 * Purpose:  Returns the path from a filename:
 *             "/foo/bar/baz"  -> "/foo/bar"
 *             "foo/bar"       -> "foo" 
 *             "foo"           -> "."
 *             "/"             -> "/" 
 *           i.e. the string will be non-NULL; it will
 *           contain the string up to but not including the
 *           last '/' character; returns "." if
 *           there are no '/' characters, and returns "/"
 *           if the last slash is the first character.
 *           Modeled on Tcl's "file dirname" command.
 * Args:     file - name of file "/foo/bar/baz".           
 * Return:   ptr to malloc'ed string "/foo/bar".          
 */
char * File_get_path(char *file)
{
  char *dirname;
  char *lastslash;
  int   len;
  
  lastslash = strrchr(file, DIRSLASH);
  len =  (lastslash == NULL) ? 0 : (int) (lastslash - file);
  dirname = (char *)tw_malloc( sizeof(char) * (len+2),
  				    "dirname",
  				    "File_get_path()" );
  if (len > 0)
  {
    strncpy(dirname, file, len);
  }
  else if (*file != DIRSLASH) 
  { 
    *dirname = '.';      
    len = 1; 
  }
  else
  { 
    *dirname = DIRSLASH; 
    len = 1; 
  }
  dirname[len] = '\0';
  return dirname;
}

/* Function: File_put_path()
 * Purpose:  Concatenate a directory path and a file name,
 *           returning a pointer to a malloc'ed string with the
 *           full filename.
 */
char * File_put_path( char *dir, char *file )
{
  char *full;

  full = (char *)tw_calloc( strlen(dir)+strlen(file)+2,
  			       sizeof(char),
  			       "full",
  			       "File_put_path()" );
  
  if ( *file == DIRSLASH ) 
  {
    strcpy(full, file); /* file = "/foo", ignore directory. */
  }
  else   
  {
    sprintf( full, "%s%c%s", dir, DIRSLASH, file );
  }
  return full;
}


/* Function: File_exists()
 * Purpose:  Return YES if filename exists.
 */
int File_exists(char *filename)
{
  FILE *fp;
  if ( (fp = fopen(filename, "r")) ) 
  { 
    fclose(fp); 
    return YES; 
  }
  return NO;
}


