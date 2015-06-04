#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.2.0) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of Colorado, Denver
#      May 2012
#//////////////////////////////////////////////////////////////////////////////

MAGMA_DIR = ..
include ../Makefile.internal

INC += -I$(MAGMA_DIR)/control

ZHDR  = commonblas_z.h

ZSRC  = \
	zauxiliary.cu	\
	zpermute.cu	\
	zpermute-v2.cu	\
	zpermute-v3.cu	\
	zswap.cu	\
	zswapblk.cu	\
	zswapdblk.cu	\
	zlacpy.cu	\
	zlange.cu	\
	zlanhe.cu	\
	zlascl.cu	\
	ztranspose.cu	\
	ztranspose-v2.cu\
	zinplace_transpose.cu	\
	zgetmatrix_transpose.cu	\
	zgetmatrix_transpose2.cu\
	zsetmatrix_transpose.cu	\
	zsetmatrix_transpose2.cu\
	zbcyclic.cu	\
	zclaswp.cu	\
	zcaxpycp.cu	\
	zlag2c.cu	\
	clag2z.cu	\
	zlat2c.cu	\
	strsm.cu	\
	dtrsm.cu	\
	ssyr2k.cu	\
	csymv_tesla.cu	\
	zhemv_fermi.cu	\
	zsymv_fermi.cu	\

# Card-specific versions that need precisions generated must be above in ZSRC,
# so they are all put into Makefile.src in the release. The code is #ifdef'd
# for the correct card.

-include Makefile.src

FERMISRC = \
	sgemv_fermi.cu	\
	dgemv_fermi.cu	\
	cgemv_fermi.cu	\
	zgemv_fermi.cu	\
	zgemvt_fermi.cu	\
	\
	sgemm_fermi.cu	\
	sgemm_fermi80.cu\
	sgemm_fermi64.cu\
	dgemm_fermi.cu	\
	zgemm_fermi.cu	\

TESLASRC = \
	sgemv_tesla.cu	\
	dgemv_tesla.cu	\
	cgemv_tesla.cu	\
	zgemv_tesla.cu	\
	gemv32_tesla.cu	\
	zsymv_tesla.cu	\
	zhemv_tesla.cu	\
	chemv_tesla.cu	\
	\
	sgemm_tesla.cu                            \
	sgemm_tesla_a_0.cu                        \
	sgemm_tesla_ab_0.cu                       \
	sgemm_tesla_N_N_64_16_16_16_4.cu          \
	sgemm_tesla_N_N_64_16_16_16_4_special.cu  \
	sgemm_tesla_N_T_64_16_4_16_4.cu           \
	sgemm_tesla_T_N_32_32_8_8_8.cu            \
	sgemm_tesla_T_T_64_16_16_16_4.cu          \
	\
	dgemm_tesla.cu                            \
	dgemm_tesla_a_0.cu                        \
	dgemm_tesla_ab_0.cu                       \
	dgemm_tesla_N_N_64_16_16_16_4.cu          \
	dgemm_tesla_N_N_64_16_16_16_4_special.cu  \
	dgemm_tesla_N_T_64_16_4_16_4.cu           \
	dgemm_tesla_T_N_32_32_8_8_8.cu            \
	dgemm_tesla_T_T_64_16_16_16_4.cu          \
	dgemm_tesla_T_T_64_16_16_16_4_v2.cu       \

ALLHDR := $(ZHDR) $(CHDR) $(DHDR) $(SHDR) \
	commonblas.h

ALLSRC := $(ZSRC) $(CSRC) $(DSRC) $(SSRC) \
	dgemv_MLU.cu	\
	stream.cu	\

ifeq (${GPU_TARGET}, Tesla)
    ALLSRC += $(TESLASRC)
    GPUOPT  = ${TESLAOPT}
else ifeq (${GPU_TARGET}, Fermi)
    ALLSRC += $(FERMISRC)
    GPUOPT  = ${FERMIOPT}
endif

ALLOBJ := $(ALLSRC:.cu=.cu_o) 

all: $(LIBMAGMABLAS) 

$(LIBMAGMABLAS): $(ALLHDR) $(ALLOBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(ALLOBJ)
	$(RANLIB) $@

%.o: %.cpp 
	$(CC) $(OPTS) $(INC) -c $< -o $@

%.cu_o: %.cu
	$(NVCC) $(NVOPTS) $(GPUOPT) $(INC) -c $< -o $@

clean:
	rm -f *.cu_o *.o *~ *.linkinfo

cleanall: clean
	rm -f *.a
	rm -f $(LIBMAGMABLAS)
