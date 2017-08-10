/****************************************************************
 * MatrixAligner (version 4a)					*
 * Copyright (C) 2006                              		*
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

/* main() function for the MatrixAligner (version 4a) program.
 */

main( int argc, char ** argv )
{
  void command_line();  
  void read_matrices();
  void init_allr_lookup_table();
  void karlin_altschul();	     
  void compare_profiles();   
  
  /* Read the command line options.
   * Print to standard output the PID and the options chosen. 
   */
  command_line(argc, argv);

  /* Read matrices */
  read_matrices();
  
  /* Initialize statistic calculation for MatAlign */
  init_allr_lookup_table();
  karlin_altschul();
  
  /* Compare matrices */
  compare_profiles();

  exit(0);
}
