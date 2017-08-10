# Description file for the "MatrixAligner" (version 4a) program.

###################################################
# Explicitly define the compiling rules assumed
COMPILE.c=$(CC) $(CFLAGS) $(CPPFLAGS) -c
.c.o:
	$(COMPILE.c) $(OUTPUT_OPTION) $<
###################################################

CFLAGS		= -O2 -w	# The basic optimized version, no warning
#CFLAGS		= -O2 -pg	# Version for profiling
#CFLAGS		= -O2 -pg -Wall # Version for profiling
#CFLAGS 	= -g         	# Version for debugging (limited optimization)
CPPFLAGS 	= -D "OPTIONS=1"
LIBS		= -lm
#CC 		= cc
CC 		= gcc
PROGRAM 	= matalign-v4a

SRC 		= main.c \
		  tw_alloc_helper.c \
		  tw_error.c \
		  tw_file.c \
		  matalign_options.c \
		  matalign_command.c \
		  alpha.c \
		  parse-line.c \
		  matalign_direction.c \
		  read_matrix.c \
		  compare_profile.c \
		  matalign_allr.c \
		  matalign_stats.c \

SRC_H		= alpha.h  \
		  matalign_options.h  \
		  matalign_struct.h  \
		  parse-line.h  \
		  tw_alloc.h  \
		  tw_error.h  \
		  tw_file.h


OBJS		= $(SRC:.c=.o)

$(PROGRAM) : $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $(OBJS) $(LIBS)

tw_alloc_helper.o      : tw_alloc.h

tw_error.o             : tw_error.h

tw_file.o	       : tw_file.h tw_alloc.h

matalign_options.o     : tw_alloc.h \
			 tw_error.h \
			 alpha.h \
			 matalign_options.h \
			 matalign_struct.h \
			 
matalign_command.o     : tw_alloc.h \
			 tw_error.h alpha.h \
			 matalign_options.h \
			 matalign_struct.h \
			 parse-line.h \
			 
parse-line.o           : tw_alloc.h \
			 tw_error.h \
			 alpha.h \
			 matalign_options.h \
			 matalign_struct.h \
			 parse-line.h \
			 
matalign_direction.o   : tw_alloc.h \
			 tw_error.h \
			 alpha.h \
			 matalign_options.h \
			 
alpha.o                : tw_alloc.h \
			 tw_error.h \
			 alpha.h \
			 matalign_options.h \
			 matalign_struct.h \
			 parse-line.h \

read_matrix.o          : tw_alloc.h \
			 tw_error.h \
			 alpha.h \
			 matalign_options.h \
			 matalign_struct.h \
			 parse-line.h \

compare_profile.o      : tw_alloc.h \
			 tw_error.h \
			 alpha.h \
			 matalign_options.h \
			 matalign_struct.h \
			 parse-line.h \
			 			 
matalign_allr.o        : tw_alloc.h \
			 tw_error.h \
			 alpha.h \
			 matalign_options.h \
			 matalign_struct.h \
			 
matalign_stats.o       : tw_alloc.h \
			 tw_error.h \
			 alpha.h \
			 matalign_options.h \
			 matalign_struct.h \
			 


clean:
	-rm -f $(PROGRAM) $(OBJS)

lint :
	lint $(CPPFLAGS) $(SRC)

tar :
	tar -chvf matalign-v4a.tar README_v4a makefile $(SRC_H) $(SRC)

