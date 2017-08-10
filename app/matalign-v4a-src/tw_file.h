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
extern char *File_get_path( char *file );


/* Function: File_put_path()
 * Purpose:  Concatenate a directory path and a file name,
 *           returning a pointer to a malloc'ed string with the
 *           full filename.
 */
extern char *File_put_path( char *dir, char *file );


/* Function: File_exists()
 * Purpose:  Return YES if filename exists.
 */
extern int File_exists( char *filename );
  
