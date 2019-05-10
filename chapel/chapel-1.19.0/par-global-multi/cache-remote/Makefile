# The source files for Chapel
CHAPEL_SRC	= heat_dissipation.chpl util.chpl

# The name of the program that we are making
PROGNAME	= heat_dissipation

# The compiler and its flags.
CC		= chpl
CFLAGS		= --fast --cache-remote

# The runtime flags.  With this flag we set environment variables that help the
# compiler to generate code for the right runtime.
RT_FLAGS	= CHPL_LAUNCHER=slurm-gasnetrun_ibv CHPL_COMM=gasnet \
		  CHPL_COMM_SUBSTRATE=ibv

# The following targets do not represent files and are therefore phony targets
.PHONY:		clean cleanslurm cleanoutput cleanall

# The program depends on the Chapel and C source
$(PROGNAME): $(CHAPEL_SRC)
	$(RT_FLAGS) $(CC) -o $@ $(CFLAGS) $<

# Clean the generated binaries
clean: 
	rm -f $(PROGNAME) $(PROGNAME)_real

# Clean the Slurm output files
cleanslurm:
	rm -f *.out

# Clean everyting
cleanall: clean cleanslurm


