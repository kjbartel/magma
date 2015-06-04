#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.4.0-beta2) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      June 2013
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

ZSRC  = blas_z.cpp	\
	zset_get.cpp	\

-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC) \
	alloc.cpp	\
	error.cpp	\
	interface.cpp	\
	set_get.cpp	\

ALLOBJ := $(ALLSRC:.cpp=.o)
ALLOBJ := $(ALLOBJ:.f=.o)

.PHONY: all clean cleanall

all: $(LIBMAGMA)

$(LIBMAGMA): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~ *.linkinfo

cleanall: clean
	rm -f $(LIBMAGMA)
