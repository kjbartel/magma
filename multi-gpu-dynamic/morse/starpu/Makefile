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
MAGMA_DIR = ../..
include $(MAGMA_DIR)/Makefile.internal

all: libmorse_starpu

libmorse_starpu: generate
	(cd codelets && $(MAKE) )
	(cd common   && $(MAKE) )

clean:
	(cd include  && $(MAKE) clean )
	(cd codelets && $(MAKE) clean )
	(cd common   && $(MAKE) clean )

cleanall: 
	(cd include  && $(MAKE) cleanall )
	(cd codelets && $(MAKE) cleanall )
	(cd common   && $(MAKE) cleanall )

generate:
	(cd include  && $(MAKE) generate )
	(cd codelets && $(MAKE) generate )
	(cd common   && $(MAKE) generate )
