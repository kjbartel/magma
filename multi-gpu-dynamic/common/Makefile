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

ZHDR =  compute_z.h
CHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p c)
DHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p "d ds")
SHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p s)

ALLHDR= common.h	\
	context.h	\
	async.h		\
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)

# ZSRC = 
# CSRC = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p c)
# DSRC = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p "d ds")
# SSRC = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p s)

ALLSRC= context.c control.c 	\
	descriptor.c		\
	async.c			\
	$(ZSRC) $(CSRC) $(DSRC) $(SSRC)

ALLOBJ = $(ALLSRC:.c=.o)

all: $(LIBMAGMAMGPU)

$(LIBMAGMAMGPU): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~

cleanall: clean
	rm -f *.a

.c.o: $(ALLHDR)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

.f.o:
	$(FC) $(FFLAGS) $(INC) -c $< -o $@

# Headers dependencies
$(ALLOBJ) : $(ALLHDR)

.PHONY: clean cleanall
