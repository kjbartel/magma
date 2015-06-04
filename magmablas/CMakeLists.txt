###
#
# @file      : CMakeLists.txt
#
# @description   : Project MAGMA
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 20-01-2012
# @last modified : mer. 23 mai 2012 14:37:32 CEST
#
###


CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(RulesPrecisions)

# Generate the magma sources for all possible precisions
# ------------------------------------------------------
SET(MAGMA_SRCS_GENERATED "")
SET(ZHDR commonblas_z.h)

# Alphabetical order
SET(ZSRC
        clag2z.cu
        csymv_tesla.cu
        dtrsm.cu
        ssyr2k.cu
        strsm.cu
        dtrsm_work.cu
        strsm_work.cu
        zauxiliary.cu
        zcaxpycp.cu
        zclaswp.cu
        zgeadd.cu
        zgetmatrix_transpose.cu
        zhemv_fermi.cu
        zinplace_transpose.cu
        zlacpy.cu
        zlag2c.cu
        zlange.cu
        zlanhe.cu
        zlascl.cu
        zlaswp.cu
        zlat2c.cu
        zsetmatrix_transpose.cu
        zswap.cu
        zswapblk.cu
        zswapdblk.cu
        zsymmetrize.cu
        zsymmetrize_tiles.cu
        zsymv_fermi.cu
        ztranspose-v2.cu
        ztranspose.cu
        dznrm2.cu
        zlarf.cu
        zlarfg-v2.cu
# Multi-gpus
        zbcyclic.cu
        zhemm_1gpu.cpp
        zhemm_1gpu_old.cpp
        zhemm_mgpu.cpp
        zhemm_mgpu_com.cpp
        zhemm_mgpu_spec.cpp
        zher2k_mgpu.cpp
        zher2k_mgpu_spec.cpp
        zgetmatrix_transpose_mgpu.cu
        zsetmatrix_transpose_mgpu.cu
        zhemv_fermi_mgpu_32.cu
        zhemv_fermi_mgpu.cu
         )

SET(FERMISRC
  sgemv_fermi.cu
  dgemv_fermi.cu
  cgemv_fermi.cu
  zgemv_fermi.cu
  zgemvt_fermi.cu
  sgemm_fermi.cu
  sgemm_fermi80.cu
  sgemm_fermi64.cu
  dgemm_fermi.cu
  zgemm_fermi.cu
  )

SET(TESLASRC
  sgemv_tesla.cu
  dgemv_tesla.cu
  cgemv_tesla.cu
  zgemv_tesla.cu
  gemv32_tesla.cu
  zsymv_tesla.cu
  zhemv_tesla.cu
  chemv_tesla.cu
  ###
  sgemm_tesla.cu
  sgemm_tesla_a_0.cu
  sgemm_tesla_ab_0.cu
  sgemm_tesla_N_N_64_16_16_16_4.cu
  sgemm_tesla_N_N_64_16_16_16_4_special.cu
  sgemm_tesla_N_T_64_16_4_16_4.cu
  sgemm_tesla_T_N_32_32_8_8_8.cu
  sgemm_tesla_T_T_64_16_16_16_4.cu
  ###
  dgemm_tesla.cu
  dgemm_tesla_a_0.cu
  dgemm_tesla_ab_0.cu
  dgemm_tesla_N_N_64_16_16_16_4.cu
  dgemm_tesla_N_N_64_16_16_16_4_special.cu
  dgemm_tesla_N_T_64_16_4_16_4.cu
  dgemm_tesla_T_N_32_32_8_8_8.cu
  dgemm_tesla_T_T_64_16_16_16_4.cu
  dgemm_tesla_T_T_64_16_16_16_4_v2.cu
  )

precisions_rules_py(MAGMABLAS_HDRS_GENERATED "${ZHDR}"
                    PRECISIONS "${MAGMA_PRECISION}"
                    DICTIONARY "MAGMA")

precisions_rules_py(MAGMABLAS_SRCS_GENERATED "${ZSRC}"
                    PRECISIONS "${MAGMA_PRECISION}"
                    DICTIONARY "MAGMA")

# Define the list of headers
# --------------------------
SET(MAGMABLAS_HDRS
    ${MAGMABLAS_HDRS_GENERATED}
    commonblas.h
   )

# Define the list of sources
# --------------------------
SET(MAGMABLAS_SRCS
    ${MAGMABLAS_SRCS_GENERATED}
    dgemv_MLU.cu
    stream.cu
   )

IF(MAGMA_USE_FERMI)
    SET(MAGMABLAS_SRCS ${MAGMABLAS_SRCS} ${FERMISRC})
ELSE()
    SET(MAGMABLAS_SRCS ${MAGMABLAS_SRCS} ${TESLASRC})
ENDIF()


# Force generation of headers
# ---------------------------
ADD_CUSTOM_TARGET(magmablas_headers ALL SOURCES ${MAGMABLAS_HDRS})
ADD_CUSTOM_TARGET(magmablas_sources ALL SOURCES ${MAGMABLAS_SRCS})
ADD_DEPENDENCIES(magmablas_headers magma_include)
ADD_DEPENDENCIES(magmablas_sources magmablas_headers)

# Compile step
# ------------
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

CUDA_ADD_LIBRARY(magmablas ${MAGMABLAS_SRCS})
ADD_DEPENDENCIES(magmablas magmablas_sources)
ADD_DEPENDENCIES(magmablas magma)

# Install step
# ------------
INSTALL(TARGETS magmablas
        DESTINATION lib)

###
### END CMakeLists.txt
###
