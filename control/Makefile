#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.2.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      May 2012
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

ZSRC := \
	zpanel_to_q.cpp	\
	zprint.cpp	\

ifeq ($(PGIFORTRAN), 0)
ZSRC := $(ZSRC) 		\
	magma_zf77.cpp		\
	magma_param.f90		\
	magma_zfortran.f90	\
	magma.f90
else
ZSRC := $(ZSRC) 		\
	magma_zf77pgi.cpp
endif

-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC)	\
	auxiliary.cpp	\
	xerbla.cpp 	\
	get_nb_tesla.cpp	\
	get_nb_fermi.cpp	\
	timer.cpp	\
	trace.cpp	\
	strlcpy.cpp	\
	bulge_auxiliary.cpp

ALLOBJ := $(ALLSRC:.cpp=.o)
ALLOBJ := $(ALLOBJ:.f90=.o)

all: $(LIBMAGMA) $(PTREXEC)

$(LIBMAGMA): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@
	cp *.mod ../include

$(PTREXEC): $(PTRFILE)
	$(CC) $(OPTS) -o $@ $<

%.o: %.cpp
	$(CC) $(OPTS) $(INC) -c $< -o $@

%.o: %.f90 $(PTREXEC)
	$(FORT) $(FOPTS) $(INC) -c $< -o $@

clean:
	rm -f *.cu_o *.o *~ *.linkinfo *.mod

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMA)

# Extra dependencies
magma.o: magma_zfortran.o magma_cfortran.o magma_dfortran.o magma_sfortran.o 
