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

MAGMA_DIR = ../../..
include $(MAGMA_DIR)/Makefile.internal

# Add the headers from src directory
INC := -I../include -I$(MAGMA_DIR)/common $(INC)

ZSRC =  codelet_zlacpy.c                                                   \
    codelet_zplrnt.c codelet_zplghe.c codelet_zplgsy.c codelet_zunmlq.c    \
    codelet_zgemm.c  codelet_ztrsm.c  codelet_zherk.c  codelet_ztsmlq.c    \
    codelet_zpotrf.c codelet_zgelqt.c codelet_zlaset.c codelet_ztslqt.c    \
    codelet_zgetrl.c codelet_ztstrf.c codelet_zgessm.c codelet_zssssm.c    \
    codelet_zgeqrt.c codelet_ztsqrt.c codelet_zunmqr.c codelet_ztsmqr.c    \
    codelet_zaxpy.c  codelet_zbrdalg.c codelet_zttqrt.c codelet_zttlqt.c   \
    codelet_zttmlq.c codelet_zsymm.c codelet_zhemm.c       \
    codelet_zher2k.c codelet_zsyr2k.c

CSRC =     $(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p c)
DSRC =     $(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p "d ds")
SSRC =     $(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p s)

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

all: $(LIBMORSE_QUARK)

$(LIBMORSE_QUARK): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.PHONY: all clean cleanall
