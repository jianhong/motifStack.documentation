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
 /* Based on Gerald Hertz(Copyright1990-1997)'s consensus package */

#ifndef OPTIONS
#include <stdio.h>
#include <ctype.h>
#include "tw_alloc.h"
#include "tw_error.h"
#include "alpha.h"
#else
#include "matalign_options.h"
#endif
#include "parse-line.h"

#define BIG_ALPHA_SIZE 512
#define ALPHA_SIZE 16


/* Functions for constructing alphabet when using the function "parse_line()".
 * 1) The file "alpha.h" must be "#included" in the file containing the
 *    description of the "parse_line" vector.
 * 2) After the command line is read, the function "adjust_alphabet()"
 *    must be called.
 * 3) The function "print_alpha()" can be used to print the alphabet
 *    information.
 * 4) In addition, the function "find_line()" is useful for other functions.
 */

/* Functions used by the function "parse_line()": declared in "parse-line.h" /
extern int pl_Alpha_af();  / Read ascii alphabet information from a file. /
extern int pl_Alpha_if();  / Read integer alphabet information from a file. /
extern int pl_Alpha_ac();  / Read alphabet information from command line. /
 */

/* The following external variables are required and
 * can be accessed by the outside. */
char *Alpha_file = "alphabet";    /* The name of the alphabet file. */
//char Comp_status = 0;    
			 /* Flag for whether to scan complementary strand.
                          * 0: ignore the complementary strand;
			  * 1: include both strands as separate sequences;
			  * 2: include both strands as a single sequence
			  *    (i.e., orientation unknown);
			  * 3: assume pattern is symmetrical. */
char Ascii = YES;        /* Flag for whether alphabet characters are
                          * ASCII (YES) or integer (NO). */
char Case_sensitive = 1; /* Flag for whether letter alphabet is case sensitive.
                          *     0: case sensitive; 1: not case sensitive;
                          *     2: not case sensitive, but indicate the
                          *        locations of lowercase letters. */
char Comp_flag = 0;      /* Flag indicating whether alphabet is complementary
                          * and whether any letters are their own complement.
                          *     0: no complements; 1: no own; 2: some own */
int A_size = 0;          /* Alphabet size---i.e., number of letters. */
int *A;                  /* Translation between an integral index and the
                          * corresponding letter/integer.  Convention for
                          * nucleic acid alphabets: if "n" is the index of a
                          * letter, its complement will have the index (n+1),
                          * if the letters are not the same. */
int *A_comp;             /* Translation between the index of a letter and
                          * the index of its complement. */
double *P;               /* The prior frequencies of the letters: obtained
                          * from the normalization data */


/* The following are the descriptions of the functions accessed by
 * the outside. */

/* Read alphabet from indicated file, if not already read from command line.
 * Make sure the alphabet is complementary if both strands are being used.
 * Convert to uppercase letters if (Case_sensitive != 0).
 * Make sure none of the characters in the alphabet occur more than once.
 * Convert the alphabet normalizations to frequencies. /
void adjust_alphabet(); */

/* Print the alphabet information. /
void print_alpha(); */

/* Find beginning of data line.  Scans over any comments (see below)
 * or whitespace.  Return values:
 *     position within the file
 *     -1L: EOF reached. /
long find_line(fp, file_name)
     FILE *fp;
     char *file_name; */

/* The following three functions are declared in "parse-line.h" and
 * used by the function "parse_line()" /
int pl_Alpha_af();      / Read ascii alphabet information from a file. /
int pl_Alpha_if();      / Read integer alphabet information from a file. /
int pl_Alpha_ac();      / Read alphabet information from the command line. /
*/



/*
The alphabet format is based on the following three mutually exclusive
alphabet options.  In addition, the portion of a line following a ';',
'%', or '#' is considered a comment and is ignored.  Comments can
begin anywhere in a line and always end at the end of the line.

     -a filename: file containing the alphabet and normalization information.
        Each line contains a letter (a symbol in the alphabet) followed by
        an optional normalization number (floating point; default is 1.0).
        In nucleic acid alphabets, a letter and its complement appear on
        the same line, separated by a colon (a letter can be its own
        complement, e.g.  when using a dimer alphabet).  The two letters
        may or may not have the same normalization number.  Only the
        standard 26 letters are permissible; however, the alphabet is case
        sensitive so that a total of 52 different characters are possible.

        POSSIBLE LINE FORMATS WITHOUT COMPLEMENTARY LETTERS:
        letter
        letter normalization

        POSSIBLE LINE FORMATS WITH COMPLEMENTARY LETTERS:
        letter:complement
        letter:complement normalization
        letter:complement normalization:complement_normalization

     -i filename: same as the "-a" option, except that the symbols of
        the alphabet are represented by integers rather than by letters.
        The number of symbols is currently limited to 127, although any
        integer permitted by the machine is a permissible symbol.

     -A alphabet_and_normalization_information: same as "-a" option, except
        information appears on the command line (e.g., -A a:t 30 c:g 20).
*/



static char STalpha_flag = 1;   /* Flag indicating which option will be
				 * used to describe the alphabet.
				 * 1: ASCII alphabet read from a file;
				 * 2: integer alphabet read from a file;
				 * 3: ASCII alphabet read from command line. */

/* Read the next character on the command line.  Returns values:
 *     an alpha letter;
 *     EOF: no alpha letter was found.
 */
static int ALPHA_read_letter(
     int argc,           /* Number of variables on the command line. */
     char *argv[],       /* The table of command line strings. */
     int *i,             /* The index to the current command line string. */
     int *k)             /* Index to current position of the current string. */
{
  while (*i < argc)
    {
      /* Skip over any whitespace. */
      for ( ; isspace((int)argv[*i][*k]) != NO; ++(*k));
      if (argv[*i][*k] == '\0')
	{
	  ++(*i);
	  (*k) = 0;
	}
      else
	break;
    }

  /* Make sure the end of the alphabet information has not been reached
   * and that the character is a letter. */
  if ((*i >= argc) || (isalpha((int)argv[*i][*k]) == NO))
    return(EOF);

  /* Return the letter. */
  else
    {
      ++(*k);
      return(argv[*i][(*k) - 1]);
    }
}

/* Read a positive floating point number from the command line.
 * Returns values:
 *     the positive number;
 *     0.0: a positive number was not found (or it was
 *          not followed by whitespace or a '\0').
 */
static double ALPHA_read_number(
     int argc,           /* Number of variables on the command line. */
     char *argv[],       /* The table of command line strings. */
     int *i,             /* The index to the current command line string. */
     int *k)             /* Index to current position of the current string. */
{
  double number = 0.0;
  char c = '\0';
  int status;

  while (*i < argc)
    {
      /* Skip over any whitespace. */
      for ( ; isspace((int)argv[*i][*k]) != NO; ++(*k));
      if (argv[*i][*k] == '\0')
	{
	  ++(*i);
	  (*k) = 0;
	}
      else
	break;
    }
  /* Make sure the end of the alphabet information has not been reached. */
  if ((*i >= argc) || (argv[*i][*k] == '-'))
    return(0.0);

  /* Make sure a number followed either by whitespace,
   * a colon, or '\0' is found. */
  status = sscanf(&argv[*i][*k], "%lf%c", &number, &c);
  if (number <= 0.0)
    return(0.0);
  else if (status == 1)
    {
      ++(*i);
      (*k) = 0;
    }
  else if (isspace((int)c) != NO)   /* Indicates that 'c' is whitespace. */
    for (++(*k) ; isspace((int)argv[*i][*k]) == NO; ++(*k));
  else if (c == ':')           /* Indicates that 'c' is a colon. */
    for (++(*k) ; argv[*i][*k] != ':'; ++(*k));
  else
    return(0.0);

  /* Return the number. */
  return(number);
}

/* Reads a colon from the command line.
 * Returns values:
 *     YES: found a colon;
 *     EOF: colon NOT found.
 */
static int ALPHA_read_colon(
     int argc,           /* Number of variables on the command line. */
     char *argv[],       /* The table of command line strings. */
     int *i,             /* The index to the current command line string. */
     int *k)             /* Index to current position of the current string. */
{
  while (*i < argc)
    {
      /* Skip over any whitespace. */
      for ( ; isspace((int)argv[*i][*k]) != NO; ++(*k));
      if (argv[*i][*k] == '\0')
	{
	  ++(*i);
	  (*k) = 0;
	}
      else
	break;
    }

  /* Make sure the end of the alphabet information has not been reached
   * and that the character is a colon. */
  if ((*i >= argc) || (argv[*i][*k] != ':'))
    return(EOF);

  /* Return a YES if the letter is a colon. */
  else
    {
      ++(*k);
      return(YES);
    }
}

/* Allocate space for "A[]", "A_comp[]", and "P[]".
 * Returns the amount of space allocated. */
static int ALPHA_alloc(void)
{
  int array_size = ALPHA_SIZE;   /* Determine the array_size. */


  A = (int *)tw_calloc(array_size, sizeof(int), "A", "ALPHA_alloc()");

  A_comp = (int *)tw_calloc(array_size, sizeof(int),
			       "A_comp", "ALPHA_alloc()");

  P = (double *)tw_calloc(array_size, sizeof(double), "P", "ALPHA_alloc()");

  return(array_size);
}

/* Reallocate space for "A[]", "A_comp[]", and "P[]".
 * Returns the amount of space reallocated. */
static int ALPHA_realloc(
     int array_size)             /* The amount of space currently allocated. */
{
  /* Determine the new array_size. */
  array_size = array_size + ALPHA_SIZE;

  A = (int *)tw_recalloc(A, array_size, sizeof(int),
			    "A", "ALPHA_realloc()");

  A_comp = (int *)tw_recalloc(A_comp, array_size, sizeof(int),
				 "A_comp", "ALPHA_realloc()");

  P = (double *)tw_recalloc(P, array_size, sizeof(double),
			       "P", "ALPHA_realloc()");

  return(array_size);
}

/* Truncate "A[]", "A_comp[]", "P[]" to proper size.*/
static void ALPHA_truncate(void)
{
  A = (int *)tw_recalloc(A, A_size, sizeof(int),
			    "A", "ALPHA_truncate()");

  if (Comp_flag == 0)
    tw_free((char *)A_comp, "A_comp", "ALPHA_truncate()");
  else
    A_comp = (int *)tw_recalloc(A_comp, A_size, sizeof(int),
				   "A_comp", "ALPHA_truncate()");

  P = (double *)tw_recalloc(P, A_size, sizeof(double),
			       "P", "ALPHA_truncate()");
}

/* Find beginning of data line in a file (cannot be used with stdin).
 * Scans over any comments or whitespace.  Return values:
 *     position within the file
 *     -1L: EOF reached.
 */
long find_line(
     FILE *fp,         /* Pointer to the file's input stream. */
     char *file)       /* The name of the input file. */
{
  long line_beginning;   /* The position within a file of the current line. */
  int letter;            /* The letter just read from the file. */

  for (line_beginning = ftell(fp), letter = getc(fp); ;
       line_beginning = ftell(fp), letter = getc(fp))
    {
      if (line_beginning == -1L) /* Make sure there was not an "ftell" error.*/
	{
	  fprintf(stderr, "\"ftell\" error in file \"%s\".\n", file);
	  exit(1);
	}
      else if (letter == EOF)      /* Double check of EOF. */
	{
	  if (feof(fp) == NO)
	    {
	      fprintf(stderr, "An error occurred while reading file \"%s\".\n",
		      file);
	      exit(1);
	    }
	  else
	    return(-1L);
	}
      else if (isspace(letter) != NO);
      else if ((letter == ';') || (letter == '%') || (letter == '#'))
	{
	  while ((letter != '\n') && (letter != EOF))
	    letter = getc(fp);
	}
      else if (fseek(fp, line_beginning, 0) != 0)
	{
	  fprintf(stderr, "\"fseek\" error in file \"%s\".\n\n", file);
	  exit(1);
	}
      else
	return(line_beginning);
    }
}

/* Read the next character in the file.  Returns values:
 *     an alpha letter;
 *     EOF: not an alpha letter.
 */
static int ALPHA_get_letter(
     FILE *fp)
{
  int letter = EOF;     /* The letter just read from the file. */

  /* Ignore whitespace unless it is a newline. */
  for (letter = getc(fp); isspace(letter) != NO; letter = getc(fp))
    {
      if (letter == '\n')
	{
	  if (ungetc(letter, fp) == EOF)
	    {
	      fprintf(stderr, "Cannot \"ungetc\" %c to file \"%s\".\n",
		      letter, Alpha_file);
	      exit(1);
	    }
	  return(EOF);
	}
    }

  /* Double check of EOF. */
  if (letter == EOF)
    {
      if (feof(fp) == NO)
	{
	  fprintf(stderr, "An error occurred while reading file \"%s\".\n",
		  Alpha_file);
	  exit(1);
	}
    }

  /* Make sure the character is a letter. */
  else if (isalpha(letter) == NO)
    {
      if (ungetc(letter, fp) == EOF)
	{
	  fprintf(stderr, "Cannot \"ungetc\" %c to file \"%s\".\n",
		  letter, Alpha_file);
	  exit(1);
	}
      letter = EOF;
    }

  return(letter);
}

/* Read a positive floating point number from the file.
 * Returns values:
 *     the positive number;
 *     0.0: positive number was not found.
 */
static double ALPHA_get_double(
     FILE *fp)
{
  double number = 0.0;   /* The number just read from the file. */
  int letter;            /* The letter just read from the file. */

  /* Ignore whitespace unless it is a newline. */
  for (letter = getc(fp); isspace(letter) != NO; letter = getc(fp))
    {
      if (letter == '\n')
	{
	  if (ungetc(letter, fp) == EOF)
	    {
	      fprintf(stderr, "Cannot \"ungetc\" '\\n' to file \"%s\".\n",
		      Alpha_file);
	      exit(1);
	    }
	  return(number);
	}
    }

  /* Double check of EOF. */
  if (letter == EOF)
    {
      if (feof(fp) == NO)
	{
	  fprintf(stderr, "An error occurred while reading file \"%s\".\n",
		  Alpha_file);
	  exit(1);
	}
    }

  /* Make sure the character is not a minus sign. */
  else if (letter == '-')
    {
      if (ungetc(letter, fp) == EOF)
	{
	  fprintf(stderr,
		  "Cannot \"ungetc\" '-' to file \"%s\".\n", Alpha_file);
	  exit(1);
	}
    }

  /* Try to read the number. */
  else
    {
      if (ungetc(letter, fp) == EOF)
	{
	  fprintf(stderr, "Cannot \"ungetc\" %c to file \"%s\".\n",
		  letter, Alpha_file);
	  exit(1);
	}
      else if ((fscanf(fp, "%lf", &number) == 1) && (number <= 0.0))
	{
	  fprintf(stderr, "The normalization value \"%G\" is ", number);
	  fprintf(stderr, "less than or equal to zero.\n");
	  exit(1);
	}
    }

  return(number);
}

/* Reads a colon from the file.
 * Returns values:
 *     YES: found a colon;
 *     EOF: colon NOT found.
 */
static int ALPHA_get_colon(
     FILE *fp)
{
  int letter;            /* The letter just read from the file. */

  /* Ignore whitespace unless it is a newline. */
  for (letter = getc(fp); isspace(letter) != NO; letter = getc(fp))
    {
      if (letter == '\n')
	{
	  if (ungetc(letter, fp) == EOF)
	    {
	      fprintf(stderr, "Cannot \"ungetc\" %c to file \"%s\".\n",
		      letter, Alpha_file);
	      exit(1);
	    }
	  return(EOF);
	}
    }

  /* Double check of EOF. */
  if (letter == EOF)
    {
      if (feof(fp) == NO)
	{
	  fprintf(stderr, "An error occurred while reading file \"%s\".\n",
		  Alpha_file);
	  exit(1);
	}
      else
	return(EOF);
    }

  else if (letter == ':')
    return(YES);

  else
    {
      if (ungetc(letter, fp) == EOF)
	{
	  fprintf(stderr, "Cannot \"ungetc\" %c to file \"%s\".\n",
		  letter, Alpha_file);
	  exit(1);
	}
      return(EOF);
    }
}

/* Read an integer from the file.  Returns values:
 *     YES: an integer was found;
 *     EOF: an integer was not found.
 */
static int ALPHA_get_int(
     FILE *fp,
     int *number)        /* The integer to be read. */
{
  int letter;            /* The letter just read from the file. */

  /* Ignore whitespace unless it is a newline. */
  for (letter = getc(fp); isspace(letter) != NO; letter = getc(fp))
    {
      if (letter == '\n')
	{
	  if (ungetc(letter, fp) == EOF)
	    {
	      fprintf(stderr, "Cannot \"ungetc\" '\\n' to file \"%s\".\n",
		      Alpha_file);
	      exit(1);
	    }
	  return(EOF);
	}
    }

  /* Double check of EOF. */
  if (letter == EOF)
    {
      if (feof(fp) == NO)
	{
	  fprintf(stderr, "An error occurred while reading file \"%s\".\n",
		  Alpha_file);
	  exit(1);
	}
      else
	return(EOF);
    }

  /* Try to read the number. */
  else
    {
      if (ungetc(letter, fp) == EOF)
	{
	  fprintf(stderr, "Cannot \"ungetc\" '%c' to file \"%s\".\n",
		  (char)letter, Alpha_file);
	  exit(1);
	}
      if (fscanf(fp, "%d", number) == 1)
	return(YES);
      else
	return(EOF);
    }
}


/* Make sure none of the characters in the alphabet occur more than once. */
static void ALPHA_repeat(void)
{
  int repeat;      /* Flag for whether any letters in alphabet are repeated. */
  int i, j;

  /* Make sure none of the characters in the alphabet occur more than once. */
  for (i = 0, repeat = NO; i < A_size; ++i)
    {
      for (j = i + 1; j < A_size; ++j)
	{
	  if (A[i] == A[j])
	    {
	      repeat = YES;
	      if (Ascii == YES)
		{
		  fprintf(stderr, "The letter \"%c\" is repe", (char)A[i]);
		  fprintf(stderr, "ated more than once in the alphabet.\n");
		}
	      else if (Ascii == NO)
		{
		  fprintf(stderr, "The integer \"%d\" is repeated ", A[i]);
		  fprintf(stderr, "more than once in the alphabet.\n");
		}
	      else bug_report("ALPHA_repeat()");
	    }
	}
    }
  if (repeat == YES)
    exit(1);
}


/* The name of the file from which to read the ascii alphabet. */
int pl_Alpha_af(
     void *variable,     /* Address of the variable to be updated. */
     int argc,           /* Number of variables on the command line. */
     char *argv[],       /* The table of command line strings. */
     int i,              /* The index to the current command line string. */
     int k)              /* Index to current position of the current string. */
{
  STalpha_flag = 1;
  Ascii = 1;
  return(pl_String((void *)&Alpha_file, argc, argv, i, k));
}

/* The name of the file from which to read the integer alphabet. */
int pl_Alpha_if(
     void *variable,     /* Address of the variable to be updated. */
     int argc,           /* Number of variables on the command line. */
     char *argv[],       /* The table of command line strings. */
     int i,              /* The index to the current command line string. */
     int k)              /* Index to current position of the current string. */
{
  STalpha_flag = 2;
  Ascii = 0;
  return(pl_String((void *)&Alpha_file, argc, argv, i, k));
}

/* Read the alphabet information from the command line. */
int pl_Alpha_ac(
     void *variable,     /* Address of the variable to be updated. */
     int argc,           /* Number of variables on the command line. */
     char *argv[],       /* The table of command line strings. */
     int i,              /* The index to the current command line string. */
     int k)              /* Index to current position of the current string. */
{
  char comp_flag = 0;    /* Same as "Comp_flag", except it refers only
			  * to the letter currently being read. */
  int letter;            /* The letter currently being read. */
  int array_size;        /* Current size of "A[]","A_comp[]",and"P[]" arrays.*/

  /* Allocate space for "A[]", "A_comp[]", and "P[]". */
  int ALPHA_alloc(void);

  /* Reallocate space for "A[]","A_comp[]",and "P[]". */
  int ALPHA_realloc(int array_size);

  /* Truncate "A[]", "A_comp[]", "P[]" to proper size.*/
  void ALPHA_truncate(void);

  /* Read the next character on the command line.  Returns values:
   *     an alpha letter;
   *     EOF: no alpha letter was found.
   */
  int ALPHA_read_letter(int argc, char *argv[], int *i, int *k);

  /* Read a positive floating point number from the command line.
   * Returns values:
   *     the positive number;
   *     0.0: positive number was not found.
   */
  double ALPHA_read_number(int argc, char *argv[], int *i, int *k);

  /* Reads a colon from the command line.
   * Returns values:
   *     YES: found a colon;
   *     EOF: colon NOT found.
   */
  int ALPHA_read_colon(int argc, char *argv[], int *i, int *k);


  STalpha_flag = 3;
  Ascii = 1;

  /* Allocate space for "A[]", "A_comp[]", and "P[]" arrays. */
  array_size = ALPHA_alloc();


  /* Read the first set of letters and normalization information. */

  /* Read the first letter. */
  if ((A[0] = ALPHA_read_letter(argc, argv, &i, &k)) == EOF)
    {
      fprintf(stderr, "The \"-A\" command line option must be followed by ");
      fprintf(stderr, "a string describing the alphabet.\n\n");
      return(NO);
    }
  else
    ++A_size;

  /* Determine whether alphabet is complementary and
   * read the complementary letter. */
  if (ALPHA_read_colon(argc, argv, &i, &k) == YES)
    {
      Comp_flag = 1;

      if ((letter = ALPHA_read_letter(argc, argv, &i, &k)) == EOF)
	{
	  fprintf(stderr, "\"%s\" (item %d on the command line) ", argv[i], i);
	  fprintf(stderr, "has the wrong format for the alphabet.\n\n");
	  return(NO);
	}
      if (A[0] == letter)
	{
	  Comp_flag = 2;
	  comp_flag = 2;
	  A_comp[0] = 0;
	}
      else
	{
	  comp_flag = 1;
	  ++A_size;
	  A[1] = letter;
	  A_comp[0] = 1;
	  A_comp[1] = 0;
	}
    }
  /* Read the normalization number. */
  if ((P[0] = ALPHA_read_number(argc, argv, &i, &k)) == 0.0)
    {
      P[0] = 1.0;
      if (comp_flag == 1)
	P[1] = 1.0;
    }
  else if (comp_flag == 1)
    {
      if (ALPHA_read_colon(argc, argv, &i, &k) == EOF)
	P[1] = P[0];
      else if ((P[1] = ALPHA_read_number(argc, argv, &i, &k)) == 0.0)
	{
	  fprintf(stderr, "\"%s\" (item %d on the command line) ", argv[i], i);
	  fprintf(stderr, "has the wrong format for the alphabet.\n\n");
	  return(NO);
	}
    }

  /* Read the rest of the letters and normalization information. */

  for ( ; ; )
    {
      /* Allocate space for "A[]", "A_comp[]", and "P[]" arrays. */
      if (array_size < A_size + 2)
	array_size = ALPHA_realloc(array_size);

      /* Read the next letter. */
      if ((A[A_size] = ALPHA_read_letter(argc, argv, &i, &k)) == EOF)
	{
	  ALPHA_truncate();

	  /* Make sure scan is at the beginning of a command line variable. */
	  if (k == 0)
	    return(i);
	  else
	    {
	      fprintf(stderr, "\"%s\" (item %d on the command lin", argv[i],i);
	      fprintf(stderr, "e) has the wrong format for the alphabet.\n\n");
	      return(NO);
	    }
	}
      else
	++A_size;

      /* Read the complement if the alphabet is complementary. */
      if (Comp_flag != 0)
	{
	  if (ALPHA_read_colon(argc, argv, &i, &k) == EOF)
	    {
	      fprintf(stderr, "\"%s\" (item %d on the command lin", argv[i],i);
	      fprintf(stderr, "e) has the wrong format for the alphabet.\n\n");
	      return(NO);
	    }
	  else if ((letter = ALPHA_read_letter(argc, argv, &i, &k)) == EOF)
	    {
	      fprintf(stderr, "\"%s\" (item %d on the command lin", argv[i],i);
	      fprintf(stderr, "e) has the wrong format for the alphabet.\n\n");
	      return(NO);
	    }
	  else if (A[A_size - 1] == letter)
	    {
	      Comp_flag = 2;
	      comp_flag = 2;
	      A_comp[A_size - 1] = A_size - 1;
	    }
	  else
	    {
	      comp_flag = 1;
	      A_comp[A_size - 1] = A_size;
	      A_comp[A_size] = A_size - 1;
	      A[A_size] = letter;
	    }
	}

      /* Read the normalization number. */
      if ((P[A_size - 1] = ALPHA_read_number(argc, argv, &i, &k)) == 0.0)
	{
	  P[A_size - 1] = 1.0;
	  if (comp_flag == 1)
	    {
	      P[A_size] = 1.0;
	      ++A_size;
	    }
	}
      else if (comp_flag == 1)
	{
	  if (ALPHA_read_colon(argc, argv, &i, &k) == EOF)
	    P[A_size] = P[A_size - 1];
	  else if ((P[A_size] = ALPHA_read_number(argc, argv, &i, &k)) == 0.0)
	    {
	      fprintf(stderr, "\"%s\" (item %d on the command lin", argv[i],i);
	      fprintf(stderr, "e) has the wrong format for the alphabet.\n\n");
	      return(NO);
	    }
	  ++A_size;
	}
    }
}

/* Read an integer alphabet from the file named in "Alpha_file". */
static void ALPHA_read_int(void)
{
  FILE *fp;
  char comp_flag = 0;    /* Same as "Comp_flag", except it refers only
			  * to the letter currently being read. */
  int letter;            /* The letter currently being read. */
  int array_size;        /* Current size of "A[]","A_comp[]",and"P[]" arrays.*/
  char error_line[BIG_ALPHA_SIZE];     /* An erroneous line. */
  long line_beginning;   /* Location---from "find-line()"---of a data line. */

  /* Allocate space for "A[]", "A_comp[]", and "P[]". */
  int ALPHA_alloc(void);

  /* Reallocate space for "A[]","A_comp[]",and "P[]". */
  int ALPHA_realloc(int array_size);

  /* Truncate "A[]", "A_comp[]", "P[]" to proper size.*/
  void ALPHA_truncate(void);

  /* Find beginning of data line.  Return values:
   *     position within the file
   *     EOF: EOF reached.
   */
  long find_line(FILE *fp, char *file);

  /* Read the next integer character in the file.  Returns values:
   *     YES: an integer was found;
   *     EOF: an integer was not found.
   */
  int ALPHA_get_int(FILE *fp, int *number);

  /* Read a positive floating point number from the file.
   * Returns values:
   *     the positive number;
   *     0.0: positive number was not found.
   */
  double ALPHA_get_double(FILE *fp);

  /* Reads a colon from the file.
   * Returns values:
   *     YES: found a colon;
   *     EOF: colon NOT found.
   */
  int ALPHA_get_colon(FILE *fp);


  /* Open the "Alpha_file". */
  if ((fp = fopen(Alpha_file, "r")) == (FILE *)NULL)
    {
      fprintf(stderr, "Cannot open the alphabet file \"%s\".\n", Alpha_file);
      exit(1);
    }
  /* Allocate space for "A[]", "A_comp[]", and "P[]" arrays. */
  array_size = ALPHA_alloc();


  /* Read the first set of letters and normalization information. */

  /* Read the first letter. */
  if ((line_beginning = find_line(fp, Alpha_file)) == EOF)
    {
      fprintf(stderr, "File \"%s\" contains no alphabet information.\n",
	      Alpha_file);
      exit(1);
    }
  if (ALPHA_get_int(fp, &A[0]) == EOF)
    {
      if (fseek(fp, line_beginning, 0) != 0)
	fprintf(stderr, "\"fseek\" error in file \"%s\".\n\n", Alpha_file);
      fprintf(stderr, "The following line in the alphabet file ");
      fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
	      Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
      exit(1);
    }
  else
    ++A_size;

  /* Determine whether alphabet is complementary and
   * read the complementary letter. */
  if (ALPHA_get_colon(fp) == YES)
    {
      Comp_flag = 1;

      if (ALPHA_get_int(fp, &letter) == EOF)
	{
	  if (fseek(fp, line_beginning, 0) != 0)
	    fprintf(stderr, "\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	  fprintf(stderr, "The following line in the alphabet file ");
	  fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		  Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	  exit(1);
	}
      else if (A[0] == letter)
	{
	  Comp_flag = 2;
	  comp_flag = 2;
	  A_comp[0] = 0;
	}
      else
	{
	  comp_flag = 1;
	  ++A_size;
	  A[1] = letter;
	  A_comp[0] = 1;
	  A_comp[1] = 0;
	}
    }

  /* Read the normalization number. */
  if ((P[0] = ALPHA_get_double(fp)) == 0.0)
    {
      P[0] = 1.0;
      if (comp_flag == 1)
	P[1] = 1.0;
    }
  else if (comp_flag == 1)
    {
      if (ALPHA_get_colon(fp) == EOF)
	P[1] = P[0];
      else if ((P[1] = ALPHA_get_double(fp)) == 0.0)
	{
	  if (fseek(fp, line_beginning, 0) != 0)
	    fprintf(stderr, "\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	  fprintf(stderr, "The following line in the alphabet file ");
	  fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		  Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	  exit(1);
	}
    }


  /* Read the rest of the letters and normalization information. */

  while ((line_beginning = find_line(fp, Alpha_file)) != EOF)
    {
      /* Allocate space for "A[]", "A_comp[]", and "P[]" arrays. */
      if (array_size < A_size + 2)
	array_size = ALPHA_realloc(array_size);

      /* Read the next letter. */
      if (ALPHA_get_int(fp, &A[A_size]) == EOF)
	{
	  if (find_line(fp, Alpha_file) == EOF)
	    {
	      ALPHA_truncate();
	      fclose(fp);
	      return;
	    }
	  else
	    {
	      if (fseek(fp, line_beginning, 0) != 0)
		fprintf(stderr,
			"\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	      fprintf(stderr, "The following line in the alphabet file ");
	      fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		      Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	      exit(1);
	    }
	}
      else
	++A_size;

      /* Read the complement if the alphabet is complementary. */
      if (Comp_flag != 0)
	{
	  if (ALPHA_get_colon(fp) == EOF)
	    {
	      if (fseek(fp, line_beginning, 0) != 0)
		fprintf(stderr,
			"\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	      fprintf(stderr, "The following line in the alphabet file ");
	      fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		      Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	      exit(1);
	    }
	  else if (ALPHA_get_int(fp, &letter) == EOF)
	    {
	      if (fseek(fp, line_beginning, 0) != 0)
		fprintf(stderr,
			"\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	      fprintf(stderr, "The following line in the alphabet file ");
	      fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		      Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	      exit(1);
	    }
	  else if (A[A_size - 1] == letter)
	    {
	      Comp_flag = 2;
	      comp_flag = 2;
	      A_comp[A_size - 1] = A_size - 1;
	    }
	  else
	    {
	      comp_flag = 1;
	      A_comp[A_size - 1] = A_size;
	      A_comp[A_size] = A_size - 1;
	      A[A_size] = letter;
	    }
	}

      /* Read the normalization number. */
      if ((P[A_size - 1] = ALPHA_get_double(fp)) == 0.0)
	{
	  P[A_size - 1] = 1.0;
	  if (comp_flag == 1)
	    {
	      P[A_size] = 1.0;
	      ++A_size;
	    }
	}
      else if (comp_flag == 1)
	{
	  if (ALPHA_get_colon(fp) == EOF)
	    P[A_size] = P[A_size - 1];
	  else if ((P[A_size] = ALPHA_get_double(fp)) == 0.0)
	    {
	      if (fseek(fp, line_beginning, 0) != 0)
		fprintf(stderr,
			"\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	      fprintf(stderr, "The following line in the alphabet file ");
	      fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		      Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	      exit(1);
	    }
	  ++A_size;
	}
    }
  ALPHA_truncate();
  fclose(fp);
  return;
}


/* Read an ascii alphabet from the file named in "Alpha_file". */
static void ALPHA_read_ascii(void)
{
  FILE *fp;
  char comp_flag = 0;    /* Same as "Comp_flag", except it refers only
			  * to the letter currently being read. */
  int letter;            /* The letter currently being read. */
  int array_size;        /* Current size of "A[]","A_comp[]",and"P[]" arrays.*/
  char error_line[BIG_ALPHA_SIZE];     /* An erroneous line. */
  long line_beginning;   /* Location---from "find-line()"---of a data line. */

  /* Allocate space for "A[]", "A_comp[]", and "P[]". */
  int ALPHA_alloc(void);

  /* Reallocate space for "A[]","A_comp[]",and "P[]". */
  int ALPHA_realloc(int array_size);

  /* Truncate "A[]", "A_comp[]", "P[]" to proper size.*/
  void ALPHA_truncate(void);

  /* Find beginning of data line.  Return values:
   *     position within the file
   *     EOF: EOF reached.
   */
  long find_line(FILE *fp, char *file);

  /* Read the next character in the file.  Returns values:
   *     an alpha letter;
   *     EOF: not an alpha letter.
   */
  int ALPHA_get_letter(FILE *fp);

  /* Read a positive floating point number from the file.
   * Returns values:
   *     the positive number;
   *     0.0: positive number was not found.
   */
  double ALPHA_get_double(FILE *fp);

  /* Reads a colon from the file.
   * Returns values:
   *     YES: found a colon;
   *     EOF: colon NOT found.
   */
  int ALPHA_get_colon(FILE *fp);


  /* Open the "Alpha_file". */
  if ((fp = fopen(Alpha_file, "r")) == (FILE *)NULL)
    {
      fprintf(stderr, "Cannot open the alphabet file \"%s\".\n", Alpha_file);
      exit(1);
    }
  /* Allocate space for "A[]", "A_comp[]", and "P[]" arrays. */
  array_size = ALPHA_alloc();


  /* Read the first set of letters and normalization information. */

  /* Read the first letter. */
  if ((line_beginning = find_line(fp, Alpha_file)) == EOF)
    {
      fprintf(stderr, "File \"%s\" contains no alphabet information.\n",
	      Alpha_file);
      exit(1);
    }
  if ((A[0] = ALPHA_get_letter(fp)) == EOF)
    {
      if (fseek(fp, line_beginning, 0) != 0)
	fprintf(stderr, "\"fseek\" error in file \"%s\".\n\n", Alpha_file);
      fprintf(stderr, "The following line in the alphabet file ");
      fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
	      Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
      exit(1);
    }
  else
    ++A_size;

  /* Determine whether alphabet is complementary and
   * read the complementary letter. */
  if (ALPHA_get_colon(fp) == YES)
    {
      Comp_flag = 1;

      if ((letter = ALPHA_get_letter(fp)) == EOF)
	{
	  if (fseek(fp, line_beginning, 0) != 0)
	    fprintf(stderr, "\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	  fprintf(stderr, "The following line in the alphabet file ");
	  fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		  Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	  exit(1);
	}
      else if (A[0] == letter)
	{
	  Comp_flag = 2;
	  comp_flag = 2;
	  A_comp[0] = 0;
	}
      else
	{
	  comp_flag = 1;
	  ++A_size;
	  A[1] = letter;
	  A_comp[0] = 1;
	  A_comp[1] = 0;
	}
    }

  /* Read the normalization number. */
  if ((P[0] = ALPHA_get_double(fp)) == 0.0)
    {
      P[0] = 1.0;
      if (comp_flag == 1)
	P[1] = 1.0;
    }
  else if (comp_flag == 1)
    {
      if (ALPHA_get_colon(fp) == EOF)
	P[1] = P[0];
      else if ((P[1] = ALPHA_get_double(fp)) == 0.0)
	{
	  if (fseek(fp, line_beginning, 0) != 0)
	    fprintf(stderr, "\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	  fprintf(stderr, "The following line in the alphabet file ");
	  fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		  Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	  exit(1);
	}
    }


  /* Read the rest of the letters and normalization information. */

  while ((line_beginning = find_line(fp, Alpha_file)) != EOF)
    {
      /* Allocate space for "A[]", "A_comp[]", and "P[]" arrays. */
      if (array_size < A_size + 2)
	array_size = ALPHA_realloc(array_size);

      /* Read the next letter. */
      if ((A[A_size] = ALPHA_get_letter(fp)) == EOF)
	{
	  if (find_line(fp, Alpha_file) == EOF)
	    {
	      ALPHA_truncate();
	      fclose(fp);
	      return;
	    }
	  else
	    {
	      if (fseek(fp, line_beginning, 0) != 0)
		fprintf(stderr,
			"\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	      fprintf(stderr, "The following line in the alphabet file ");
	      fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		      Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	      exit(1);
	    }
	}
      else
	++A_size;

      /* Read the complement if the alphabet is complementary. */
      if (Comp_flag != 0)
	{
	  if (ALPHA_get_colon(fp) == EOF)
	    {
	      if (fseek(fp, line_beginning, 0) != 0)
		fprintf(stderr,
			"\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	      fprintf(stderr, "The following line in the alphabet file ");
	      fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		      Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	      exit(1);
	    }
	  else if ((letter = ALPHA_get_letter(fp)) == EOF)
	    {
	      if (fseek(fp, line_beginning, 0) != 0)
		fprintf(stderr,
			"\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	      fprintf(stderr, "The following line in the alphabet file ");
	      fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		      Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	      exit(1);
	    }
	  else if (A[A_size - 1] == letter)
	    {
	      Comp_flag = 2;
	      comp_flag = 2;
	      A_comp[A_size - 1] = A_size - 1;
	    }
	  else
	    {
	      comp_flag = 1;
	      A_comp[A_size - 1] = A_size;
	      A_comp[A_size] = A_size - 1;
	      A[A_size] = letter;
	    }
	}

      /* Read the normalization number. */
      if ((P[A_size - 1] = ALPHA_get_double(fp)) == 0.0)
	{
	  P[A_size - 1] = 1.0;
	  if (comp_flag == 1)
	    {
	      P[A_size] = 1.0;
	      ++A_size;
	    }
	}
      else if (comp_flag == 1)
	{
	  if (ALPHA_get_colon(fp) == EOF)
	    P[A_size] = P[A_size - 1];
	  else if ((P[A_size] = ALPHA_get_double(fp)) == 0.0)
	    {
	      if (fseek(fp, line_beginning, 0) != 0)
		fprintf(stderr,
			"\"fseek\" error in file \"%s\".\n\n", Alpha_file);
	      fprintf(stderr, "The following line in the alphabet file ");
	      fprintf(stderr, "\"%s\" has the wrong format:\n   \"%s\"\n",
		      Alpha_file, fgets(error_line, BIG_ALPHA_SIZE, fp));
	      exit(1);
	    }
	  ++A_size;
	}
    }
  ALPHA_truncate();
  fclose(fp);
  return;
}


/* Read alphabet from indicated file, if not already read from command line.
 * Make sure the alphabet is complementary if both strands are being used.
 * Convert to uppercase letters if (Case_sensitive != 0).
 * Make sure none of the characters in the alphabet occur more than once.
 * Convert complementary letters to the same a priori probabilities
 * if both strands are being used.
 * Convert the alphabet normalizations to frequencies. */
void adjust_alphabet(void)
{
  int i, i_comp;
  double norm_sum;                /* Sum of the normalizations. */
  void ALPHA_read_ascii(void);    /* Read an ascii alphabet from a file. */
  void ALPHA_read_int(void);      /* Read an integer alphabet from a file. */
  void ALPHA_repeat(void);        /* Make sure no letters are repeated. */


  /* Read the alphabet from the indicated file, unless it has
   * already been read from the command line. */
  if (STalpha_flag == 1) ALPHA_read_ascii();
  else if (STalpha_flag == 2) ALPHA_read_int();
  else if (STalpha_flag != 3) bug_report("adjust_alpha()");

  /* Make sure the alphabet is complementary if both strands are being used. */
  if((Comp_status != 0) && (Comp_flag == 0))
    {
      fprintf(stderr,"Complemetary alphabet is mandatory when a complement ");
      fprintf(stderr,"option is chosen\n\n");
      exit(1);
    }

  /* Convert to uppercase letters if (Case_sensitive != 0). */
  if ((Ascii == YES) && (Case_sensitive != 0))
    {
      for (i = 0; i < A_size; ++i)
	if (islower(A[i]) != NO) A[i] = toupper(A[i]);
    }

  /* Make sure none of the characters in the alphabet occur more than once. */
  ALPHA_repeat();

  /* Convert complementary letters to the same a priori probabilities
   * if both strands are being used. */
  if (Comp_status != 0)
    {
      for (i = 0; i < A_size; ++i)
	{
	  i_comp = A_comp[i];
	  P[i] = P[i_comp] = (P[i] + P[i_comp]) / 2.0;
	  i = i_comp;
	}
    }

  /* Convert the alphabet normalizations to frequencies. */
  for (i = 0, norm_sum = 0.0; i < A_size; ++i) norm_sum = norm_sum + P[i];
  for (i = 0; i < A_size; ++i) P[i] = P[i] / norm_sum;
}



/* Print the alphabet informations. */
void print_alpha(void)
{
  int i, k;

  if (STalpha_flag == 1)
    {
      printf("#***** Prior frequencies for the alphabet from ");
      printf("file \"%s\". *****#\n", Alpha_file);
    }
  else if (STalpha_flag == 2)
    {
      printf("***** Information for the integral alphabet from ");
      printf("file \"%s\". *****\n", Alpha_file);
    }
  else if (STalpha_flag == 3)
    {
      printf("***** Information for the alphabet from ");
      printf("the command line. *****\n");
    }
  else bug_report("print_alpha()");

  /* Print an ascii alphabet. */
  if (Ascii == YES)
    {
      if (Comp_flag == 0)
	{
	  for(i = 0; i < A_size; ++i)
	    printf("letter %3d: %c  prior frequency = %#G\n",
		   i + 1, A[i], P[i]);
	}
      else if (Comp_flag == 1)
	{
	  for(i = 0, k = 1; i < A_size; i = i + 2, ++k)
	    {
	      printf("letter %3d: %c (complement: ", k, A[i]);
	      printf("%c)  prior frequency = %#G\n", A[A_comp[i]], P[i]);
	    }
	  for(i = A_size - 1; i > 0; i = i - 2, ++k)
	    {
	      printf("letter %3d: %c (complement: ", k, A[i]);
	      printf("%c)  prior frequency = %#G\n", A[A_comp[i]], P[i]);
	    }
	}
      else if (Comp_flag == 2)
	{
	  for(i = 0; i < A_size; ++i)
	    {
	      printf("letter %3d: %c (complement: ", i + 1, A[i]);
	      printf("%c)  prior frequency = %#G\n", A[A_comp[i]], P[i]);
	    }
	}
      else
	{
	  fprintf(stderr, "The variable \"Comp_flag\" has an illegal value ");
	  fprintf(stderr, "in the \"print_alpha()\" function.\n");
	}
    }

  /* Print an integer alphabet. */
  else if (Ascii == NO)
    {
      if (Comp_flag == 0)
	{
	  for(i = 0; i < A_size; ++i)
	    printf("letter %3d: %3d  prior frequency = %#G\n",
		   i + 1, A[i], P[i]);
	}
      else if (Comp_flag == 1)
	{
	  for(i = 0, k = 1; i < A_size; i = i + 2, ++k)
	    {
	      printf("letter %3d: %3d (complement: ", k, A[i]);
	      printf("%3d)  prior frequency = %#G\n", A[A_comp[i]], P[i]);
	    }
	  for(i = A_size - 1; i > 0; i = i - 2, ++k)
	    {
	      printf("letter %3d: %3d (complement: ", k, A[i]);
	      printf("%3d)  prior frequency = %#G\n", A[A_comp[i]], P[i]);
	    }
	}
      else if (Comp_flag == 2)
	{
	  for(i = 0; i < A_size; ++i)
	    {
	      printf("letter %3d: %3d (complement: ", i + 1, A[i]);
	      printf("%3d)  prior frequency = %#G\n", A[A_comp[i]], P[i]);
	    }
	}
      else
	{
	  fprintf(stderr, "The variable \"Comp_flag\" has an illegal value ");
	  fprintf(stderr, "in the \"print_alpha()\" function.\n");
	}
    }

  else
    {
      fprintf(stderr, "The variable \"Ascii\" has an illegal value in ");
      fprintf(stderr, "the \"print_alpha\" function.\n");
      exit(1);
    }
}

/******************
 * End of alpha.c *
 ******************/

