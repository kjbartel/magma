###
#
# @file Makefile
#
#  MAGMA Makefile
#  MAGMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver,
#  and INRIA Bordeaux Sud-Ouest
#
#  @version 2.3.1
#  @author Mathieu Faverge
#  @date 2011-06-01
#
###

# Overwritten in make.inc
MAGMA_DIR = ..
include $(MAGMA_DIR)/Makefile.internal

all: 
	(cd include && $(MAKE) )

ifeq (${MORSE_USE_QUARK}, 1)
all: libmorse_quark
endif

ifeq (${MORSE_USE_STARPU}, 1)
all: libmorse_starpu
endif

libmorse_quark:
	(cd quark && $(MAKE))

libmorse_starpu:
	(cd starpu && $(MAKE))

clean:
	(cd quark  && $(MAKE) clean )
	(cd starpu && $(MAKE) clean )

cleanall: 
	(cd include  && $(MAKE) cleanall )
	(cd quark    && $(MAKE) cleanall )
	(cd starpu   && $(MAKE) cleanall )

generate:
	(cd include  && $(MAKE) generate )
	(cd quark    && $(MAKE) generate )
	(cd starpu   && $(MAKE) generate )
