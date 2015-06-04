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
SET(ZSRC zauxiliary.cu
         zpermute.cu
         zpermute-v2.cu
         zpermute-v3.cu
         zswap.cu
         zswapblk.cu
         zswapdblk.cu
         zlacpy.cu
         zlange.cu
         zlanhe.cu
         zlascl.cu
         ztranspose.cu
         ztranspose-v2.cu
         zinplace_transpose.cu
         zgetmatrix_transpose.cu
         zgetmatrix_transpose2.cu
         zsetmatrix_transpose.cu
         zsetmatrix_transpose2.cu
         zbcyclic.cu
         zclaswp.cu
         zcaxpycp.cu
         zlag2c.cu
         clag2z.cu
         zlat2c.cu
         strsm.cu
         dtrsm.cu
         ssyr2k.cu
         csymv_tesla.cu
         zhemv_fermi.cu
         zsymv_fermi.cu
         zsymmetrize.cu
         )

SET(FERMISRC sgemv_fermi.cu
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

SET(TESLASRC sgemv_tesla.cu
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

precisions_rules_py(MAGMABLAS_HDRS_GENERATED
            "${ZHDR}"
            PRECISIONS ${MAGMA_PRECISION})

precisions_rules_py(MAGMABLAS_SRCS_GENERATED
            "${ZSRC}"
            PRECISIONS ${MAGMA_PRECISION})

# Define the list of headers
# --------------------------
SET(MAGMABLAS_HDRS ${MAGMABLAS_HDRS_GENERATED}
           commonblas.h
           )

# Define the list of sources
# --------------------------
SET(MAGMABLAS_SRCS ${MAGMABLAS_SRCS_GENERATED}
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
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/control)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/control)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/magmablas)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/magmablas)

CUDA_ADD_LIBRARY(magmablas ${MAGMABLAS_SRCS})

ADD_DEPENDENCIES(magmablas magmablas_sources)
ADD_DEPENDENCIES(magmablas magma)

FOREACH(_external_package "blas" "lapack" "plasma")
    STRING(TOUPPER "${_external_package}" _NAMEVAR)
    IF(DEFINED ${_NAMEVAR}_BUILD_MODE)
        ADD_DEPENDENCIES(magmablas ${_external_package}_install)
    ENDIF()
ENDFOREACH()

# Install step
# ------------
INSTALL(TARGETS magmablas
    DESTINATION lib)

###
### END CMakeLists.txt
###
