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


/* Declaration of functions defined in "alpha.c" and
 * used to read alphabet information. */

/* Read alphabet from indicated file.
 * Make sure the alphabet is complementary if both strands are being used.
 * Make sure none of the characters in the alphabet occur more than once.
 * Convert complementary letters to the same a priori probabilities
 * if both strands are being used.
 * Convert the alphabet normalizations to frequencies. */
extern void adjust_alphabet(void);

/* Print the alphabet information. */
extern void print_alpha(void);

/* Find beginning of data line in a file (cannot be used with stdin).
 * Scans over any comments or whitespace.  Return values:
 *     position within the file
 *     -1L: EOF reached.
 */
extern long find_line(
     FILE *fp,         /* Pointer to the file's input stream. */
     char *file);      /* The name of the input file. */


/* Functions used by the function "parse_line()": declared in "parse-line.h" /
extern int pl_Alpha_af();  / Read ascii alphabet information from a file. /
extern int pl_Alpha_if();  / Read integer alphabet information from a file. /
extern int pl_Alpha_ac();  / Read alphabet information from command line. /
 */


/* The following external variables are required and
 * can be accessed by the outside. */
extern char *Alpha_file; /* The name of the alphabet file. */

extern char Ascii;       /* Flag for whether alphabet characters are
                          * ASCII (YES) or integer (NO). */
extern char Case_sensitive;/* Flag whether letter alphabet is case sensitive.
                          *     0: case sensitive; 1: not case sensitive;
                          *     2: not case sensitive, but indicate the
                          *        locations of lowercase letters. */
extern char Comp_flag;   /* Flag indicating whether alphabet is complementary
                          * and whether any letters are their own complement.
                          *     0: no complements; 1: no own; 2: some own */
extern int A_size;       /* Alphabet size---i.e., number of letters. */
extern int *A;           /* Translation between an integral index and the
                          * corresponding letter/integer.  Convention for
                          * nucleic acid alphabets: if "n" is the index of a
                          * letter, its complement will have the index (n+1),
                          * if the letters are not the same. */
extern int *A_comp;      /* Translation between the index of a letter and
                          * the index of its complement. */
extern double *P;        /* The prior frequencies of the letters: obtained
                          * from the normalization data */
