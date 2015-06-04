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
include ../Makefile.internal

INC := -I$(MAGMA_DIR)/common/ $(INC)

ZHDR =  
CHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p c)
DHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p "d ds")
SHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p s)

ZSRC =  zgels.c  zgelqf.c zgelqs.c zgeqrf.c zgeqrs.c 	\
        zgesv.c  zgetrf.c zgetrs.c             		\
        zposv.c  zpotrf.c zpotrs.c             		\
	ztrsm.c  zgemm.c  zsymm.c zsyrk.c ztrsmpl.c  	\
	zpotri.c ztrtri.c zlauum.c ztrmm.c zsyr2k.c	\
	zlange.c zlansy.c				\
	zunmqr.c zunmlq.c zungqr.c zunglq.c 		\
	zhemm.c  zherk.c  zlanhe.c zher2k.c   		\
	zcposv.c zcgesv.c zcungesv.c			\
        zgecfi2.c zgecfi.c zgetmi.c                     \
	zplrnt.c zplghe.c zplgsy.c			\
	\
	pzgelqf.c pzgeqrf.c pzpotrf.c pztrtri.c pzlauum.c  	\
        pzgelqfrh.c pzgeqrfrh.c                                 \
        pztrsm.c  pzgemm.c  pzsymm.c  pzsyrk.c  pztrsmpl.c 	\
	pzlange.c pzlansy.c pzsyr2k.c				\
	pzlacpy.c pzaxpy.c  pzgetrf.c pztrmm.c			\
	pzunmqr.c pzunmlq.c pzungqr.c pzunglq.c 		\
	pzunmqrrh.c pzunmlqrh.c pzungqrrh.c pzunglqrh.c         \
        pzlanhe.c pzhemm.c pzherk.c pzher2k.c			\
	pzlag2c.c                                               \
        pzpack.c  pzshift.c pzgetmi2.c				\

ZSRC =  zplrnt.c zplghe.c zplgsy.c ztile.c 		\
	zgemm.c ztrsm.c ztrsmpl.c 			\
	zpotrf.c zpotrs.c zposv.c     			\
	zgetrf_incpiv.c zgetrs_incpiv.c zgesv_incpiv.c	\
	zgeqrf.c					\
	zlauum.c ztrtri.c zpotri.c                      \
	\
	pzplrnt.c pzplghe.c pzplgsy.c pztile.c 		\
	pzgemm.c pztrsm.c pztrsmpl.c			\
	pzpotrf.c pzgetrf_incpiv.c pzgeqrf.c		\
	pztrtri.c pzlauum.c				\

# 	\
#	zpotrf_extra.c					\
#	zpotri_extra.c                                	\
#	ztrtri_extra.c                                	\


CSRC = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p c)
DSRC = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p "d ds")
SSRC = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p s)

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)

all: $(LIBMAGMAMGPU)

$(LIBMAGMAMGPU): $(HEADERS) $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

zgecfi2.c zgecfi.c : $(ZHDR)
cgecfi2.c cgecfi.c : $(CHDR)
dgecfi2.c dgecfi.c : $(DHDR)
sgecfi2.c sgecfi.c : $(SHDR)

.PHONY: clean cleanall
