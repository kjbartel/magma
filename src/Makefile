#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.2.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      June 2012
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

INC += -I$(MAGMA_DIR)/control

# ------------------------------------------------------------
# Loosely, the order is alphabetic in each section, but mixed precision,
# driver, and multi-GPU routines are separated where it made sense.
# Double precision "d" files that generate "s" files must also be listed
# in ZSRC, despite the name.

# Cholesky, GPU interface
ZSRC += \
	zcposv_gpu.cpp		\
	\
	zposv_gpu.cpp		\
	zpotrf_gpu.cpp		\
	zpotri_gpu.cpp		\
	zpotrs_gpu.cpp		\
	zlauum_gpu.cpp		\
	ztrtri_gpu.cpp		\
	\
	zpotrf_mgpu.cpp		\

# Cholesky, CPU interface
ZSRC += \
	zposv.cpp		\
	zpotrf.cpp		\
	zpotrf2_ooc.cpp		\
	zpotrf3_mgpu.cpp	\
	zpotri.cpp		\
	zlauum.cpp		\
	ztrtri.cpp		\

# LU, GPU interface
ZSRC += \
	zcgesv_gpu.cpp		\
	zcgetrs_gpu.cpp		\
	\
	zgesv_gpu.cpp		\
	zgetrf_gpu.cpp		\
	zgetrf_nopiv_gpu.cpp	\
	zgetri_gpu.cpp		\
	zgetrs_gpu.cpp		\
	\
	zgetrf_mgpu.cpp		\

# LU, auxiliary for multi-gpu-dynamic
ZSRC += \
	zgetrl_gpu.cpp		\
	zgessm_gpu.cpp		\
	zssssm_gpu.cpp		\
	ztstrf_gpu.cpp		\

# LU, CPU interface
ZSRC += \
	zgesv.cpp		\
	zgetrf.cpp		\
	zgetrf3_ooc.cpp		\
	zgetrf1_mgpu.cpp	\
	zgetf2_nopiv.cpp	\
	zgetrf_nopiv.cpp	\

# QR and least squares, GPU interface
ZSRC += \
	zcgeqrsv_gpu.cpp	\
	\
	zgelqf_gpu.cpp		\
	zgels3_gpu.cpp		\
	zgels_gpu.cpp		\
	zgeqrf2_gpu.cpp		\
	zgeqrf3_gpu.cpp		\
	zgeqrf_gpu.cpp		\
	zgeqrs3_gpu.cpp		\
	zgeqrs_gpu.cpp		\
	zlarfb_gpu.cpp		\
	zungqr_gpu.cpp		\
	zunmql2_gpu.cpp		\
	zunmqr2_gpu.cpp		\
	zunmqr_gpu.cpp		\
	\
	zgeqp3.cpp		\
	dgeqp3.cpp		\
	zlaqps.cpp		\
	\
	zgeqrf_mgpu.cpp		\

# QR, CPU interface
ZSRC += \
	zgelqf.cpp		\
	zgeqlf.cpp		\
	zgeqrf-v4.cpp		\
	zgeqrf.cpp		\
	zgeqrf_ooc.cpp		\
	zungqr.cpp		\
	zunmql.cpp		\
	zunmqr.cpp		\
	zunmqr_m.cpp		\

# symmetric eigenvalue, GPU interface
ZSRC += \
	dsyevd_gpu.cpp		\
	zheevd_gpu.cpp		\
	zheevdx_gpu.cpp		\
	zheevr_gpu.cpp		\
	zheevx_gpu.cpp		\
	zhetrd2_gpu.cpp		\
	zhetrd_gpu.cpp		\
	zunmtr_gpu.cpp		\

# symmetric eigenvalues, CPU interface
ZSRC += \
	dsyevd.cpp		\
	dsyevdx.cpp		\
	zheevd.cpp		\
	zheevdx.cpp		\
	zheevr.cpp		\
	zheevx.cpp		\
	\
	dlaex0.cpp		\
	dlaex1.cpp		\
	dlaex3.cpp		\
	dstedx.cpp		\
	zhetrd.cpp		\
	zlatrd.cpp		\
	zlatrd2.cpp		\
	zstedx.cpp		\
	zungtr.cpp		\
	zunmtr.cpp		\
	\
	zheevd_m.cpp		\
	dlaex0_m.cpp		\
	dlaex1_m.cpp		\
	dlaex3_m.cpp		\
	dstedx_m.cpp		\
	zstedx_m.cpp		\
	zunmtr_m.cpp		\

# symmetric eigenvalues 2-stage
ZSRC += \
	dbulge_aux.cpp		\
	zbulge_applyQ.cpp	\
	zbulge_aux.cpp		\
	zbulge_kernel.cpp	\
	zhebbd.cpp		\
	zhetrd_bhe2trc_v3.cpp	\
	zungqr_2stage_gpu.cpp	\

# generalized symmetric eigenvalues
ZSRC += \
	dsygvd.cpp		\
	zhegst.cpp		\
	zhegvd.cpp		\
	zhegvdx.cpp		\
	zhegvr.cpp		\
	zhegvx.cpp		\
	zhegst_gpu.cpp		\
	\
	zhegst_m.cpp		\
	zhegvd_m.cpp		\
	ztrsm_m.cpp		\

# non-symmetric eigenvalues
ZSRC += \
	dgeev.cpp		\
	zgeev.cpp		\
	zgehrd.cpp		\
	zgehrd2.cpp		\
	zlahr2.cpp		\
	zlahru.cpp		\
	zunghr.cpp		\

# SVD
ZSRC += \
	dgesvd.cpp		\
	zgebrd.cpp		\
	zgesvd.cpp		\
	zlabrd_gpu.cpp		\


# ------------------------------------------------------------
-include Makefile.local
-include Makefile.src

ALLSRC_F = dlaed2.f \
           slaed2.f

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC)
ALLOBJ := $(ALLSRC:.cpp=.o) $(ALLSRC_F:.f=.o)

all: $(LIBMAGMA)

$(LIBMAGMA): $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

clean:
	rm -f *.o *~ *.linkinfo

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMA)
