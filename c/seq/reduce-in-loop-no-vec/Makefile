SRC		= heat_dissipation.c util.c
HDR		= util.h
PROGNAME	= heat_dissipation

INCLUDE_DIRS	= -I.
LIB_DIRS	= -L.
LIB		= 

CC		= gcc
CFLAGS		= -O3 -Wall -fno-tree-vectorize -fopt-info-vec-optimized $(INCLUDE_DIRS)

LD		= gcc
LDFLAGS		= $(LIB_DIRS) $(LIB)

OBJ		= $(SRC:.c=.o)

.PHONY:		clean


$(PROGNAME): $(OBJ)
	$(LD) -o $@ $(OBJ) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	rm -f $(PROGNAME) *.o $(DEPFILE)

