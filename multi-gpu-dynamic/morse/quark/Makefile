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

all: libmorse_quark

libmorse_quark:
	(cd common   && $(MAKE) )
	(cd codelets && $(MAKE) )

clean:
	(cd common   && $(MAKE) clean )
	(cd codelets && $(MAKE) clean )

cleanall: 
	(cd common   && $(MAKE) cleanall )
	(cd codelets && $(MAKE) cleanall )

generate:
	(cd common   && $(MAKE) generate )
	(cd codelets && $(MAKE) generate )
