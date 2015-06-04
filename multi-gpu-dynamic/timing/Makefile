###
#
# @file Makefile
#
#  MAGMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.3.1
# @author Mathieu Faverge
# @date 2010-11-15
#
###

MAGMA_DIR = ..
include ../Makefile.internal

CFLAGS  += -I../control/ $(INC)
LDFLAGS = $(LIB)

ZHDR =  zauxiliary.h zauxiliary.c
CHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p c)
DHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p "d ds")
SHDR = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZHDR)" -o -p s)

#timing.c is in headers because it's a template for every timing routine
ALLHDR= auxiliary.h auxiliary.c timing.h timing.c \
	$(ZHDR) $(CHDR) $(DHDR) $(SHDR)

ZSRC  = time_zgemm.c 		\
	time_zgesv.c 		\
	time_zgetrf.c 		\
	time_zgels.c 		\
	time_zgeqrf.c 		\
	time_zgemm_tile.c 	\
	time_zpotri_tile.c	\
	time_zgesv_tile.c 	\
	time_zgetrf_tile.c 	\
	time_zgels_tile.c 	\
	time_zgeqrf_tile.c      \
        time_zgecfi.c           \
        time_zlapack2tile.c

ZSRC  = time_zgemm.c			\
	time_zgemm_tile.c 		\
	time_zgeqrf_tile.c		\
	time_zgesv_incpiv.c		\
	time_zgesv_incpiv_tile.c	\
	time_zgetrf_incpiv.c		\
	time_zgetrf_incpiv_tile.c	\
	time_zposv_tile.c 		\
	time_zpotrf_tile.c		\
	time_zpotri_tile.c		\
	time_ztrsm.c

#	time_zgels.c			\
#	time_zgels_tile.c		\
#	time_zgeqrf.c			\

CSRC = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p c)
DSRC = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p "d ds")
SSRC = 	$(shell $(PYTHON) $(CODEGEN) -f "$(ZSRC)" -o -p s)

ALLSRC = $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ = $(ALLSRC:.c=.o)
EXE = $(ALLOBJ:.o=)

all: $(EXE)

$(EXE) : $(LIBMAGMAMGPU) $(LIBMORSE)

clean:
	rm -f *.o *~

cleanall: clean
	rm -f $(EXE)

time_%.o : time_%.c $(ALLHDR)
	$(CC) $(CFLAGS) -c $< -o $@

%auxiliary.o : %auxiliary.c $(ALLHDR)
	$(CC) $(CFLAGS) -c $< -o $@

time_s% : time_s%.o sauxiliary.o
	  $(LOADER) $^ -o $@ $(LDFLAGS)
time_d% : time_d%.o dauxiliary.o
	$(LOADER) $^ -o $@ $(LDFLAGS)
time_c% : time_c%.o cauxiliary.o
	  $(LOADER) $^ -o $@ $(LDFLAGS)
time_z% : time_z%.o zauxiliary.o
	  $(LOADER) $^ -o $@ $(LDFLAGS)


time_zlapack2tile time_clapack2tile time_dlapack2tile time_slapack2tile : auxiliary.o
time_zgecfi time_cgecfi time_dgecfi time_sgecfi : auxiliary.o

png: time_dgesv time_zgesv
	./time_dgesv --gnuplot-png | tee timeplot.gpt ; gnuplot timeplot.gpt

.PHONY: all clean
