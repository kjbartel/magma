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
MAGMA_DIR = .
include ./Makefile.internal

all: lib timings
#test example timings

lib: generation libmorse libmagma_mgpu

libmorse: 
	( cd morse && $(MAKE) )

libmagma_mgpu:
	(cd common  && $(MAKE) )
	(cd compute && $(MAKE) )

test: testplasma testlapack

testplasma: lib
	(cd testing && $(MAKE))

testlapack: lib
	(cd testing/lin && $(MAKE))

example: lib
	(cd examples && $(MAKE) )

timings: lib
	(cd timing && $(MAKE) )

clean:
	(cd include     && $(MAKE) clean )
	(cd compute     && $(MAKE) clean )
	(cd common      && $(MAKE) clean )
	(cd morse       && $(MAKE) clean )
#	(cd testing     && $(MAKE) clean )
#	(cd testing/lin && $(MAKE) clean )
#	(cd examples    && $(MAKE) clean )
	(cd timing      && $(MAKE) clean )

cleanall: 
	(cd include     && $(MAKE) cleanall )
	(cd compute     && $(MAKE) cleanall )
	(cd common      && $(MAKE) cleanall )
	(cd morse       && $(MAKE) cleanall )
#	(cd testing     && $(MAKE) cleanall )
#	(cd testing/lin && $(MAKE) cleanall )
#	(cd examples    && $(MAKE) cleanall )
	(cd timing      && $(MAKE) cleanall )
	(cd lib         && rm -f *.a )

##############################################################
#            Installation

dir:
	mkdir -p $(prefix)
	mkdir -p $(prefix)/include
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/lib/pkgconfig

install: dir
#       MAGMA_MGPU
	cp $(MAGMA_DIR)/include/*.h   $(prefix)/include
	cp $(LIBCODELETS)             $(prefix)/lib
	cp $(LIBMORSE)                $(prefix)/lib
	cp $(LIBMAGMA_MGPU)           $(prefix)/lib
#       pkgconfig
	cat $(MAGMA_DIR)/lib/pkgconfig/magma_mgpu.pc.in | 	\
	    sed -e s:\__PREFIX:"$(prefix)":          	| 	\
	    sed -e s:\__LIBEXT:"$(LIBEXT)":          	| 	\
	    sed -e s:\__REQUIRE:"$(require)":          		\
	    > $(prefix)/lib/pkgconfig/magma_mgpu.pc

