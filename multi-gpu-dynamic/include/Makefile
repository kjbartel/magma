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

MAGMA_DIR = ..
include ../Makefile.internal

ZHDR =  magma_morse_z.h
CHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p c)
DHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p "d ds")
SHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p s)

ALLHDR= magma_morse.h	\
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)


all: $(ALLHDR)

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a
