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


#ifndef OPTIONS
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "tw_alloc.h"
#else
#include "matalign_options.h"
#endif
#include "parse-line.h"

/*
 * Functions used for processing the command line.
 *    int parse_line: read the command line and
 *                    check for presence of required variables;
 *    void usage: a usage description;
 */

/* Command specific functions.  
 * Returns the index of the next command line variable. 
 * Return values <= 0 indicate an error.
 * The following is its format: /
int function(
      void *variable,     / Address of the variable to be updated. /
      int argc,           / Number of variables on the command line. /
      char *argv[],       / The table of command line strings. /
      int i,              / The index to the current command line string. /
      int k              / Index to current position of the current string. /
      )

 * The general command specific functions not defined in "parse-line.c".
 *    pl_Help: page directions from the function "text_directions(FILE *fp)";
 *    pl_Alpha_af: read ascii alphabet information from a file;
 *    pl_Alpha_if: read integer alphabet information from a file;
 *    pl_Alpha_ac: read alphabet information from command line;

 * The general command specific functions defined in "parse-line.c".
 *    pl_NBool: assign the integer 0 to the named char;
 *    pl_Bool: assign the integer 1 to the named char;
 *    pl_Int_2: assign the integer 2 to the named char;
 *    pl_Int_3: assign the integer 3 to named char;
 *    pl_Int_4: assign the integer 4 to named char;
 *    pl_IChar: assign the following int to the named char;
 *    pl_Int: assign the following int to the named int;
 *    pl_P_Int: assign the following positive int to the named int;
 *    pl_Nn_Int: assign the following non-negative int to the named int;
 *    pl_Double: assign the following double to the named double;
 *    pl_P_Double: assign the following positive double to the named double;
 *    pl_Nn_Double: assign the following non-negative double to named double;
 *    pl_String: assign the following string to the named (char *).
 *    pl_Char: assign the last char of the option's name to the named char.
 */


/* "pl_Help()" is defined in the file "parse-line-help.c".
 * If pl_Help is called, a function "text_directions()" for printing
 * the text of the directions needs to be defined in another file.
   void text_directions(
	  FILE *fp      / Output stream for printing the directions. /
          )
   {
     fprintf(fp, "Text of the directions.\n");
   }
 */

/* The functions "pl_Alpha_af()", "pl_Alpha_if()", and "pl_Alpha_ac()"
 * for reading alphabet information are defined in "alpha.c". */


/*
 * Read the command line.  Returns 1 if everything is OK; otherwise,
 * returns 0.  The options are strings; however, no option can be a
 * prefix of another option.  The options are described in an array of
 * OPTION structures.  An OPTION structure has the following entries:
 *   (1) the command line string (the '-' is assumed); 
 *   (2) macro determining whether option is required (REQ), optional  
 *       (OPT), or hidden (HID); 
 *   (3) pointer to the function that will process the command line option,
 *       returns the index to the next command line variable to be 
 *       processed by parse_line (functions are described and declared 
 *       in "parse_line.h"); 
 *   (4) the address of the variable to be changed (cast to (void *)); 
 *   (5) string describing the particular option.  
 * Mutually exclusive options are described in an array of (char *)'s.
 * Each set of mutually exclusive options are separated by the macro END; 
 * in addition, an extra END is placed at the end of the list. 
 * A NULL pointer indicates that none of the options are mutually exclusive.
 */


int parse_line(
      OPTION *options,      /* Vector describing the options. */
      char **exclusive,     /* Vector indicating mutually exclusive options. */
      int argc,             /* Number of variables on the command line. */
      char *argv[]          /* The table of command line strings. */
      )
{
  char *flag;         /* Flag for whether an option has been used. */
  int *length;        /* The length of each option string. */
  int status = 1;     /* The return status of the program. */
  int num_opt;        /* The number of options listed in "options". */
  int flag_sum;       /* The sum of "flag[]" when testing for
		       * mutually exclusive options. */
  char opt_flag = ' ';      /* "opt_flag" entry for set of mutually exclusive 
                       * options. */
  int i, j, k, l;


  /* Determine "num_opt". The last item in the OPTION vector is
   * always "END". The following for-loop increment num_opt until
   * the "END" is reached. */
  for ( num_opt = 0; options[num_opt].opt_string[0] != '\0'; ++num_opt );

  /* Allocate space for and initialize "flag[]" and "length[]". */
  flag = (char *)tw_calloc( num_opt, sizeof(char), "flag", "parse_line()" );
  length = (int *)tw_calloc( num_opt, sizeof(int), "length", "parse_line()" );
  for ( i = 0; i < num_opt; ++i )
  {
    flag[i] = 0;    /* flag this option to be "not used yet" */
    length[i] = strlen( options[i].opt_string );
  }

  /* Read the command line. 
   * Process each option read and advance to the next option, until
   * all the command line argument is processed */
  for ( i = 1; i < argc; )
  {
    /* Make sure the command begins with a '-'. */
    if ( argv[i][0] != '-' )
    {
      fprintf( stderr, "\"%s\" (item %d on the command line) ", argv[i], i );
      fprintf( stderr, "does not begin with a \"-\".\n\n" );
      return( NO );
    }

    /* Check the command line against the different options. */
    for ( j = 0; j < num_opt; ++j )
    {
      if ( !strncmp( &argv[i][1], options[j].opt_string, length[j] ) )
      {
	/* Mark the "flag[]" array. */
	if ( flag[j] == 1 )
	{
	  fprintf( stderr,"Option \"-%s\" occu", options[j].opt_string );
	  fprintf( stderr,"rs more than once on the command line.\n\n" );
	  return( NO );
	}
	else
	{
	  flag[j] = 1;
	}

	/* Determine where the command specific function should look
	 * for its arguments. */
	k = 1 + length[j];
	if ( argv[i][k] == '\0' )
	{
	  ++i;
	  k = 0;
        }
	/* Command specific function.  Returns the index
	 * of the next command line variable. */
	i = (*options[j].opt_func)( options[j].opt_var, argc, argv, i, k );
	if ( i <= 0 ) 
	{
	  return( NO );
	}
	else 
	{
	  break;
	}
      }
    }

    /* Make sure the command line matched an option. */
    if ( j == num_opt )
    {
      fprintf( stderr, "\"%s\" (item %d on the command line) ", argv[i], i );
      fprintf( stderr, "does not match any of the legal options.\n\n" );
      return( NO );
    }
  }

  /* Check for mutually exclusive options. */
  if ( exclusive != (char **)NULL )
  {
    for ( i=0; exclusive[i][0]!='\0'; i=i+k+1 )
    {
      for ( k=0, flag_sum=0; exclusive[i+k][0]!='\0'; ++k )
      {
	for ( j=0; ( j<num_opt ) &&
		   ( strcmp(exclusive[i+k], options[j].opt_string) != 0 ); 
		   ++j )
        {
	  ;
	}

	/* Make sure the options in exclusive[] occur in options[]. */
	if ( j == num_opt )
	{
	  fprintf( stderr, "There is an error in the " );
	  fprintf( stderr, "\"parse_line\" function.\n" );
	  fprintf( stderr, "An option in the \"exclusive\" array does" );
	  fprintf( stderr, " not appear in the \"options\" array.\n" );
	  exit( 1 );
	}

	/* Make sure if one mutually exclusive option is required,
	 * so are the others. */
	if ( k == 0 )
	{
	  opt_flag = options[j].opt_flag;
	}
	else if ( (opt_flag != options[j].opt_flag) &&
		  ((opt_flag == REQ) || (options[j].opt_flag == REQ)) )
	{
	  fprintf( stderr, "There is an error in the " );
	  fprintf( stderr, "\"parse_line\" function.\n" );
	  fprintf( stderr, "A required option is mutually \n" );
	  fprintf( stderr, "exclusive with a non-required option.\n" );
	  exit( 1 );
	}

	/* Count the number of mutually exclusive options used. */
	flag_sum += (int)flag[j];

	/* If the options are required, mark them all as being used,
	 * but make sure "flag_sum" equals 1 below. */
	if ( opt_flag == REQ )
	{
	  flag[j] = 1;
	}
      }

      /* Make sure only 1 mutually exclusive option is used. */
      if ( flag_sum > 1 )
      {
	 status = 0;
	 fprintf( stderr, "The " );
	 for ( l=0; l<k-1; ++l )
	 {
	   fprintf( stderr, "\"-%s\", ", exclusive[i + l] );
	 }
	 fprintf( stderr,"and \"-%s\" options are mutually exclusive.\n\n",
		  exclusive[i + l] );
      }

      /* Make sure one representative of a required option is used. */
      if ( (opt_flag == REQ) && (flag_sum == 0) )
      {
	status = 0;
	fprintf( stderr, "One of the following mutually exclusive " );
	fprintf( stderr, "options is required on the command\nline: " );
	for ( l=0; l<k-1; ++l )
	{
	  fprintf( stderr, "\"-%s\", ", exclusive[i + l] );
	}
	fprintf( stderr,"and \"-%s\".\n\n", exclusive[i + l] );
      }
    }
  }

  /* Make sure all the required options have been given. */
  for ( j = 0; j < num_opt; ++j )
  {
    if ( (options[j].opt_flag == REQ) && (flag[j] == 0) )
    {
      status = 0;
      fprintf( stderr, "The required option \"-%s\" is missing\n",
	       options[j].opt_string );
    }
  }

  /* Free dynamic memory. */
  tw_free(flag, "flag", "parse_line()");
  tw_free(length, "length", "parse_line()");

  return(status);
  
}


/* Usage message based on the "opt_descrpt" entries of the OPTION array.
 * Also allows for an additional message.  If "message" is handed a NULL
 * pointer or a '\0', the message will be ignored. */
void usage(
       OPTION *options,      /* Vector describing the options. */
       char *message,        /* An additional help message. */
       char *argv[]         /* The table of command line strings. */
       )
{
  int i;

  fprintf( stderr, "Usage: %s\n", argv[0] );
  if ( (message != (char *)NULL) && (message[0] != '\0') )
  {
    fprintf(stderr, "%s\n", message);
  }

  /* Print the REQuired options. */
  for ( i = 0; options[i].opt_string[0] != '\0'; ++i )
  {
    if ( options[i].opt_flag == REQ )
    {
      if ( (options[i].opt_descrpt == (char *)NULL) ||
	   (options[i].opt_descrpt[0] == '\0') )
      {
	fprintf(stderr, "-%s\n", options[i].opt_string);
      }
      else
      {
        fprintf( stderr, "-%s <%s>\n",
		 options[i].opt_string, options[i].opt_descrpt );
      }
    }
  }

  /* Print the OPTional options. */
  for (i = 0; options[i].opt_string[0] != '\0'; ++i)
    {
      if (options[i].opt_flag == OPT)
	{
	  if ((options[i].opt_descrpt == (char *)NULL) ||
	      (options[i].opt_descrpt[0] == '\0'))
	    fprintf(stderr, "[-%s]\n", options[i].opt_string);
	  else
	    fprintf(stderr, "[-%s <%s>]\n",
		   options[i].opt_string, options[i].opt_descrpt);
	}
    }

  exit(1);

}


/* Command specific functions.  Returns the index of the next
 * command line variable.  Return values <= 0 indicate an error.
 * The following is its format: /
int function(
     void *variable,     / Address of the variable to be updated. /
     int argc,           / Number of variables on the command line. /
     char *argv[],       / The table of command line strings. /
     int i,              / The index to the current command line string. /
     int k)              / Index to current position of the current string. /
 */


/**************************
 * The general functions. *
 **************************/

/* Assign the integer 0 to the named char. */
int pl_NBool(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  if ( k != 0 )
  {
    fprintf( stderr, "\"%s\" (item %d on the command line) ", argv[i], i );
    fprintf( stderr, "does not match any of the legal options.\n\n" );
    return( NO );
  }
  else
  {
    *(char *)variable = 0;
    return( i );
  }
}


/* Assign the integer 1 to the named char. */
int pl_Bool(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k              /* Index to current position of the current string. */
      )
{
  if ( k != 0 )
  {
    fprintf( stderr, "\"%s\" (item %d on the command line) ", argv[i], i );
    fprintf( stderr, "does not match any of the legal options.\n\n" );
    return( NO );
  }
  else
  {
    *(char *)variable = 1;
    return( i );
  }
}


/* Assign the integer 2 to the named char. */
int pl_Int_2(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  if ( k != 0 )
  {
    fprintf( stderr, "\"%s\" (item %d on the command line) ", argv[i], i );
    fprintf( stderr, "does not match any of the legal options.\n\n" );
    return( NO );
  }
  else
  {
    *(char *)variable = 2;
    return( i );
  }
}


/* Assign the integer 3 to the named char. */
int pl_Int_3(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  if ( k != 0 )
  {
    fprintf( stderr, "\"%s\" (item %d on the command line) ", argv[i], i );
    fprintf( stderr, "does not match any of the legal options.\n\n" );
    return( NO );
  }
  else
  {
    *(char *)variable = 3;
    return(i);
  }
}


/* Assign the integer 4 to the named char. */
int pl_Int_4(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  if ( k != 0 )
  {
    fprintf( stderr, "\"%s\" (item %d on the command line) ", argv[i], i );
    fprintf( stderr, "does not match any of the legal options.\n\n" );
    return(NO);
  }
  else
  {
    *(char *)variable = 4;
    return(i);
  }
}


/* Assign the following int to the named char. */
int pl_IChar(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  int int_variable;
  char c = '\0';

  if ( i >= argc )
  {
    fprintf( stderr,"The command line is missing an integer at its end.\n\n" );
    return( NO );
  }

  else if ( (sscanf( &argv[i][k], "%d%c", &int_variable, &c) != 1)
	    || (c!='\0') )
  {
    fprintf( stderr, "\"%s\" (from item %d on the command line) ",
	     &argv[i][k], i );
    fprintf( stderr, "is not an integer.\n\n" );
    return( NO );
  }

  else
  {
    *(char *)variable = (char)int_variable;
    return(i + 1);
  }
  
}


/* Assign the following int to the named int. */
int pl_Int(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  char c = '\0';

  if ( i >= argc )
  {
    fprintf( stderr,"The command line is missing an integer at its end.\n\n" );
    return( NO );
  }

  else if ( (sscanf(&argv[i][k], "%d%c", (int *)variable, &c) != 1)
	    || (c != '\0') )
  {
    fprintf( stderr, "\"%s\" (from item %d on the command line) ",
	     &argv[i][k], i );
    fprintf( stderr, "is not an integer.\n\n" );
    return( NO );
  }

  else
  {
    return(i + 1);
  }
  
}


/* Assign the following positive int to the named int. */
int pl_P_Int(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  char c = '\0';

  if ( i >= argc )
  {
    fprintf( stderr, "The command line is missing a positive " );
    fprintf( stderr, "integer at its end.\n\n" );
    return( NO );
  }

  else if ( (sscanf(&argv[i][k], "%d%c", (int *)variable, &c) != 1)
	    || (c != '\0') || (*((int *)variable) <= 0) )
  {
    fprintf( stderr, "\"%s\" (from item %d on the command line) ",
	     &argv[i][k], i );
    fprintf( stderr, "is not a positive integer.\n\n" );
    return( NO );
  }

  else
  {
    return(i + 1);
  }
  
}


/* Assign the following non-negative int to the named int. */
int pl_Nn_Int(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  char c = '\0';

  if ( i >= argc )
  {
    fprintf( stderr, "The command line is missing a non-negative " );
    fprintf( stderr, "integer at its end.\n\n" );
    return( NO );
  }

  else if ( (sscanf(&argv[i][k], "%d%c", (int *)variable, &c) != 1)
	    || (c != '\0') || (*((int *)variable) < 0) )
  {
    fprintf( stderr, "\"%s\" (from item %d on the command line) ",
	     &argv[i][k], i );
    fprintf( stderr, "is not a non-negative integer.\n\n" );
    return( NO );
  }

  else
  {
    return(i + 1);
  }
  
}


/* Assign the following double to the named double. */
int pl_Double(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  char c = '\0';

  if ( i >= argc )
  {
    fprintf( stderr, "The command line is missing a floating " );
    fprintf( stderr, "point number at its end.\n\n" );
    return( NO );
  }

  else if ( (sscanf(&argv[i][k], "%lf%c", (double *)variable, &c) != 1)
	    || (c != '\0') )
  {
    fprintf( stderr, "\"%s\" (from item %d on the command line) ",
	     &argv[i][k], i );
    fprintf( stderr, "is not a floating point number.\n\n" );
    return( NO );
  }

  else
  {
    return(i + 1);
  }
  
}


/* Assign the following positive double to the named double. */
int pl_P_Double(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k              /* Index to current position of the current string. */
      )
{
  char c = '\0';

  if ( i >= argc )
  {
    fprintf( stderr, "The command line is missing a positive floating " );
    fprintf( stderr, "point number at its end.\n\n" );
    return( NO );
  }

  else if ( (sscanf(&argv[i][k], "%lf%c", (double *)variable, &c) != 1)
	    || (c != '\0') || (*((double *)variable) <= 0) )
  {
    fprintf( stderr, "\"%s\" (from item %d on the command line) ",
	     &argv[i][k], i );
    fprintf( stderr, "is not a positive floating point number.\n\n" );
    return( NO );
  }

  else
  {
    return(i + 1);
  }
  
}


/* Assign the following non-negative double to the named double. */
int pl_Nn_Double(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k              /* Index to current position of the current string. */
      )
{
  char c = '\0';

  if ( i >= argc )
  {
    fprintf( stderr, "The command line is missing a non-negative floating " );
    fprintf( stderr, "point number at its end.\n\n" );
    return( NO );
  }

  else if ( (sscanf(&argv[i][k], "%lf%c", (double *)variable, &c) != 1 )
	    || (c != '\0') || (*((double *)variable) < 0) )
  {
    fprintf( stderr, "\"%s\" (from item %d on the command line) ",
	     &argv[i][k], i );
    fprintf( stderr, "is not a non-negative floating point number.\n\n" );
    return( NO );
  }

  else
  {
    return(i + 1);
  }
}


/* Assign the following string to the named (char *). */
int pl_String(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  if ( i < argc )
  {
    *(char **)variable = &argv[i][k];
    return(i + 1);
  }

  else
  {
    fprintf( stderr, "The command line is missing a " );
    fprintf( stderr, "string of characters at its end\n\n" );
    return( NO );
  }
  
}


/* Assign the last char of the option's name to the named char. */
int pl_Char(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k              /* Index to current position of the current string. */
      )
{
  int length;

  if (k == 0) 
  {
    --i;
  }

  length = strlen(argv[i]);
  *(char *)variable = argv[i][length - 1];
  return(i + 1);
  
}


/* Call the function "print_directions". */
int pl_Help(
      void *variable,     /* Address of the variable to be updated. */
      int argc,           /* Number of variables on the command line. */
      char *argv[],       /* The table of command line strings. */
      int i,              /* The index to the current command line string. */
      int k               /* Index to current position of the current string. */
      )
{
  void print_directions(void);

  if ( k != 0 )
  {
    fprintf( stderr, "\"%s\" (item %d on the command line) ", argv[i], i );
    fprintf( stderr, "does not match any of the legal options.\n\n" );
    return( NO );
  }
  else
  {
    print_directions();
    exit( 0 );
  }
}

/* Page the directions. */
void print_directions( void )
{
  FILE *fp;                    /* Pointer to the paging format. */
  char *pager;                 /* Value of the environment variable PAGER. */
  void text_directions( FILE *fp );      /* The text of the directions. */

#ifndef unix
  text_directions( stdout );
#else
  /* Determine the format for printing the directions. */
  pager = getenv( "PAGER" );

  /* Determine the pointer to the printing format. */
  if ( pager == (char *)NULL )
  {
    fp = popen( "more", "w" );
  }
  else
  {
    fp = popen( pager, "w" );
  }

  /* Print the directions and close. */
  if ( fp == (FILE *)NULL )
  {
    fp = stdout;
    text_directions(fp);
  }
  else
  {
    text_directions(fp);
    pclose(fp);
  }
#endif
}

