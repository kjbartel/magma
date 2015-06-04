#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      November 2011
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

INC += -I$(MAGMA_DIR)/control

ZSRC  = zpotrf.cpp	\
        zpotrf_gpu.cpp	\
        zpotrf_ooc.cpp	\
        zpotrf_mgpu.cpp \
        zpotrf2_ooc.cpp \
        zpotrf2_mgpu.cpp \
        zpotrs_gpu.cpp	\
        zposv.cpp	\
        zposv_gpu.cpp	\
        zpotri.cpp	\
        zpotri_gpu.cpp	\
	\
        zgetrf.cpp	\
        zgetrf-v2.cpp   \
        zgetrf-v3.cpp   \
        zgetrf_gpu.cpp 	\
        zgetrf_mgpu.cpp 	\
        zgetrf_nopiv.cpp\
        zgetrf2_ooc.cpp \
	zgetf2_nopiv.cpp\
        zgetrf_nopiv_gpu.cpp\
        zgetrs_gpu.cpp 	\
        zgesv.cpp       \
        zgesv_gpu.cpp 	\
        zgetri_gpu.cpp 	\
	\
	zgeqrf.cpp	\
	zgeqrf-v4.cpp   \
        zgeqrf_ooc.cpp  \
	zgeqlf.cpp	\
	zgelqf.cpp	\
        zgelqf_gpu.cpp  \
        zungqr.cpp	\
        zungtr.cpp      \
	zgeqrf_gpu.cpp  \
        zgeqrf_mgpu-trace.cpp \
	zgeqrf2_gpu.cpp	\
	zgeqrf3_gpu.cpp \
	zgeqrs_gpu.cpp	\
	zgeqrs3_gpu.cpp \
	zgels_gpu.cpp	\
	zgels3_gpu.cpp  \
	zlarfb_gpu.cpp	\
        zunmqr_gpu.cpp	\
        zungqr_gpu.cpp	\
	\
	zgebrd.cpp 	\
	zlabrd_gpu.cpp	\
	zgehrd.cpp 	\
        zgehrd2.cpp     \
        zlahru.cpp 	\
        zlahr2.cpp 	\
        zunghr.cpp      \
        zhetrd.cpp 	\
        zhetrd_gpu.cpp  \
        zhetrd2_gpu.cpp \
        zlatrd.cpp 	\
        zlatrd2.cpp     \
        zunmtr.cpp 	\
        zunmtr_gpu.cpp  \
        zunmqr.cpp 	\
        zunmqr2_gpu.cpp \
	zunmql.cpp      \
        zunmql2_gpu.cpp \
        \
        zgeev.cpp       \
        dgeev.cpp       \
        zheevd.cpp      \
        zheevd_gpu.cpp  \
	zheevdx_gpu.cpp \
	zheevr.cpp      \
	zheevr_gpu.cpp  \
	zheevx.cpp      \
	zheevx_gpu.cpp  \
        dsyevd.cpp      \
        dsyevd_gpu.cpp  \
        zhegvd.cpp      \
	zhegvdx.cpp     \
        zhegvx.cpp      \
        zhegvr.cpp      \
        dsygvd.cpp      \
        zhegst.cpp      \
        zhegst_gpu.cpp  \
        zgesvd.cpp	\
	dgesvd.cpp	\
	\
	zcposv_gpu.cpp	\
	zcgesv_gpu.cpp	\
	zcgetrs_gpu.cpp	\
	zcgeqrsv_gpu.cpp\
	\
	zgessm_gpu.cpp	\
	zssssm_gpu.cpp  \
	zgetrl_gpu.cpp	\
	ztstrf_gpu.cpp  \
	\
	ztrtri.cpp	\
	ztrtri_gpu.cpp	\
	zlauum.cpp	\
	zlauum_gpu.cpp


-include Makefile.src

ALLSRC= $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ  = $(ALLSRC:.cpp=.o)

all: $(LIBMAGMA)

$(LIBMAGMA): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

.cpp.o:
	$(CC) $(OPTS) $(INC) -c $< -o $@

clean:
	rm -f *.o *~ *.linkinfo

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMA)
