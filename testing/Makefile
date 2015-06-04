#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.3.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      November 2012
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

# ----------
# BLAS
ZSRC += \
	testing_zgemm.cpp	\
	testing_zgemv.cpp	\
	testing_zhemv.cpp	\
	testing_zsymv.cpp	\
	\
	testing_auxiliary.cpp	\
	testing_zblas.cpp	\
	testing_zgeadd.cpp	\
	testing_zlacpy.cpp	\
	testing_zsymmetrize.cpp	\
	testing_zswap.cpp	\
	testing_zutil.cpp	\
	testing_zhemv_mgpu.cpp	\
#	\
#	testing_zhemm_mgpu.cpp	\
#	testing_zher2k_mgpu.cpp	\

# ----------
# Cholesky, GPU interface
ZSRC += \
	testing_zcposv_gpu.cpp	\
	\
	testing_zposv_gpu.cpp	\
	testing_zpotrf_gpu.cpp	\
	testing_zpotri_gpu.cpp	\
	\
	testing_zpotrf_mgpu.cpp	\

# Cholesky, CPU interface
ZSRC += \
	testing_zposv.cpp	\
	testing_zpotrf.cpp	\
	testing_zpotri.cpp	\

# ----------
# LU, GPU interface
ZSRC += \
	testing_zcgesv_gpu.cpp	\
	\
	testing_zgesv_gpu.cpp	\
	testing_zgetrf_gpu.cpp	\
	testing_zgetri_gpu.cpp	\
	\
	testing_zgetrf_mgpu.cpp	\

# LU, CPU interface
ZSRC += \
	testing_zgesv.cpp	\
	testing_zgetrf.cpp	\

# ----------
# QR and least squares, GPU interface
ZSRC += \
	testing_zcgeqrsv_gpu.cpp	\
	\
	testing_zgelqf_gpu.cpp	\
	testing_zgeqrf_gpu.cpp	\
	testing_zgeqr2_gpu.cpp  \
	testing_zgeqrs3_gpu.cpp	\
	testing_zgeqrs_gpu.cpp	\
	testing_zlarfb_gpu.cpp	\
	testing_zungqr_gpu.cpp	\
	testing_zunmqr_gpu.cpp	\
	\
	testing_zgeqrf_mgpu.cpp	\

# QR, CPU interface
ZSRC += \
	testing_zgelqf.cpp	\
	testing_zgeqlf.cpp	\
	testing_zgeqrf.cpp	\
	testing_zungqr.cpp	\
	testing_zunmqr.cpp	\
	\
	testing_zgeqp3.cpp	\

# ----------
# symmetric eigenvalues, GPU interface
ZSRC += \
	testing_dsyevd_gpu.cpp	\
	testing_zheevd_gpu.cpp	\
	testing_zhetrd_gpu.cpp	\

# symmetric eigenvalues, CPU interface
ZSRC += \
	testing_dsyevd.cpp	\
	testing_zheevd.cpp	\
	testing_zhetrd.cpp	\

# symmetric eigenvalues, 2-stage
ZSRC += \
	testing_zhetrd_he2hb.cpp	\
	testing_zhegvdx_2stage.cpp	\
#	testing_zhetrd_he2hb_mgpu.cpp	\

# generalized symmetric eigenvalues
ZSRC += \
	testing_dsygvd.cpp	\
	testing_zhegvd.cpp	\
	testing_zhegvdx.cpp	\

# ----------
# non-symmetric eigenvalues
ZSRC += \
	testing_dgeev.cpp	\
	testing_zgeev.cpp	\
	testing_zgehrd.cpp	\

# ----------
# SVD
ZSRC += \
	testing_zgesvd.cpp	\
	testing_zgebrd.cpp	\

# ----------
ifeq ($(FORT), pgfortran)
	FOBJ  = fortran_thunking.o
	ZSRCF = testing_zgetrf_gpu_f.cuf
else
	FOBJ  = fortran.o
	ZSRCF = testing_zgetrf_gpu_f.F90	\
		testing_zgetrf_f.f90
endif

ZSRC += $(ZSRCF)


# ------------------------------------------------------------
-include Makefile.local
-include Makefile.src

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ := $(ALLSRC:.f90=.o)
ALLOBJ := $(ALLOBJ:.F90=.o)
ALLOBJ := $(ALLOBJ:.cuf=.o)
ALLOBJ := $(ALLOBJ:.cpp=.o)
EXE    := $(ALLOBJ:.o=)

LIBTEST := libtest.a
LIBSRC  := testing_zutil.o testing_cutil.o testing_dutil.o testing_sutil.o testing_util.o
LIBOBJ  := $(LIBSRC:.cpp=.o)
LIBEXE  := $(LIBOBJ:.o=)

# strip out fake exe files generated from LIBSRC
EXE := $(filter-out $(LIBEXE), $(EXE))

# only testing_auxiliary needs control
INC += -I../control

all: $(EXE)

lib: $(LIBMAGMA) $(LIBMAGMABLAS)

$(LIBTEST): $(LIBOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $^
	$(RANLIB) $@

LIBLAPACKTEST = lin/liblapacktest.a

$(LIBLAPACKTEST):
	cd lin && make liblapacktest.a

$(EXE): $(LIBMAGMA) $(LIBMAGMABLAS) $(LIBLAPACKTEST) $(LIBTEST)

clean:
	rm -f *.o *~ $(LIBTEST)

cleanall: clean
	rm -f $(EXE)

ifeq ($(FORT), pgfortran)
testing_%.o : testing_%.cuf
	$(FORT) $(FOPTS) $(INC) -c $< -o $@ 
else
testing_%.o : testing_%.f90
	$(FORT) $(FOPTS) $(INC) -c $< -o $@ 
endif

testing_%_f : testing_%_f.o $(FOBJ)
	$(FORT) $(LDOPTS) $< $(FOBJ) -o $@ \
	$(LIBTEST) $(LIBLAPACKTEST) -L../lib -lmagma -lmagmablas -lmagma \
	$(LIBDIR) \
	$(LIB) 

testing_% : testing_%.o
	$(CC) $(LDOPTS) $(NOMAIN) $< $(OBJ) -o $@ \
	$(LIBTEST) $(LIBLAPACKTEST) -L../lib -lmagma -lmagmablas -lmagma \
	$(LIBDIR) \
	$(LIB) 

#
# These two files are provided by cuda
#
fortran.o : $(CUDADIR)/src/fortran.c
	$(CC) $(OPTS) $(INC) -c $< -o $@

fortran_thunking.o : $(CUDADIR)/src/fortran_thunking.c
	$(CC) $(OPTS) $(INC) -c $< -o $@

# keep intermediate object files, to speed up re-compiling when library changes
.SECONDARY: $(ALLOBJ)
