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
#include "parse-line.h"

#define MESSAGE "Use the \"-h\" option for detailed directions"

/* Functions specific to this particular program and
 * needed in the OPTION vector. */
extern int pl_Help();    /* Call the function "print_directions". */

/* Vector holding the possible command line options. */
static OPTION ST_options[] =
{
  "h",   OPT, pl_Help,      NULL,              "Print directions",
  "f1",  REQ, pl_String,    &file1,            "Name of matrix list file 1",
  "f2",  OPT, pl_String,    &file2,            "Name of matrix list file 2",
  "c0",  OPT, pl_NBool,     &Comp_status,      "Compare only the forward orientation",
  "c1",  OPT, pl_Bool,      &Comp_status,      "Compare both orientations (default)",
  "g",   OPT, pl_Bool,      &Global,           "Global alignment, (default: 0, for local alignment)",
  "ez",  OPT, pl_Bool,      &EZ,	       "For simple output (only ALLR and Overlap)", 
  "a",   OPT, pl_Alpha_af,  NULL,              "Name of ascii alphabet file",
  "i",   OPT, pl_Alpha_if,  NULL,              "Name of integer alphabet file",
  "A",   OPT, pl_Alpha_ac,  NULL,              "Ascii alphabet information on the command line",
  "CS",  OPT, pl_NBool,     &Case_sensitive,   "Ascii alphabet is case sensitive (default: ascii alphabets are case insensitive)",
  END
};


/* Vector indicating mutually exclusive options. */
static char *ST_exclusive[] =
{
  "c0", "c1", END,
  "a", "i", "A", END,
  "CS", "i", END,
  END
};


/* Read the command line options.  Print to standard output the PID
 * and the options chosen. */
void command_line( int argc, char ** argv)
{
  void check_options();
  void print_options();
  int  getpid();
  void usage(OPTION *, char *, char **); /* defined in parse-line.c */

  /* Read the command line. */
  if ( parse_line( ST_options, ST_exclusive, argc, argv ) == NO )
  {
    usage( ST_options, MESSAGE, argv );
  }

  /* Determine the process ID. */
  PID = (int)getpid();

  /* Read alphabet from indicated file, if not already read from command
   * line.  Make sure the alphabet is complementary if both strands are
   * being used.  Convert to uppercase letters if (Case_sensitive != 0).
   * Make sure none of the characters in the alphabet occur more than
   * once.  Convert the alphabet normalizations to frequencies. */
  adjust_alphabet();

  /* Check for options not yet implemented. */
  check_options();

  /* Print the options chosen. */
  print_options( argc, argv );

}


/* Print the options chosen. */
void print_options( int argc, char *argv[] )
{
  int i;
  char *indent = "                   ";

  /* Echo the command line. */
  printf( "COMMAND LINE:" );
  for ( i = 0; i < argc; ++i ) 
  {
    printf( " %s", argv[i] );
  }
  printf( "\n\n" );

  /* Display general information */
  printf( "***** PID: %d *****\n", PID );
  
  if ( List_num == 1 )
  {
    printf ( "Pair-wise comparison among matrices specified in %s.\n", file1 );
  }
  else if ( List_num == 2 )
  {
    printf ( "Compare all matrices specified in %s to all matrices specified in %s.\n", file1, file2 );
  }
  else
  {
    printf ( "Number of matrix list is %d. What's going on?\n", List_num );
    exit(0);
  }
  
  printf ( "Algorithm options:\n" );
  if ( Comp_status )
  {
    printf ( "  Compare both orientations.\n" );
  }
  else
  {
    printf ( "  Compare only the forward orientation.\n" );
  }
  if ( Global )
  {
    printf ( "  Global alignment.\n" );
  }
  else
  {
    printf ( "  Local alignment.\n" );
  }
  
  printf ( "Output option:\n" );
  if ( EZ )
  {
    printf ( "  Simple output, only names of matrices, ALLR and Overlap.\n" );
  }
  else
  {
    printf ( "  Complete output, including all statistics.\n" );
  }

}

void check_options()
{
  if ( file2 )
  {
    List_num = 2;
  }
  else
  {
    List_num = 1;
  }
}
