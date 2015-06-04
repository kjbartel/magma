###
#
# @file      : CMakeLists.txt
#
# @description   : Project MAGMA
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 20-01-2012
# @last modified : mer. 16 mai 2012 11:37:18 CEST
#
###


CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(RulesPrecisions)

# Generate the magma sources for all possible precisions
# ------------------------------------------------------
SET(MAGMA_CTRL_GENERATED "")
SET(MAGMA_CTRLF_GENERATED "")
SET(ZSRCF "")
SET(ZSRC 
         ../control/zprint.cpp
         ../control/zpanel_to_q.cpp
         ../interface_cuda/blas_z.cpp
         ../interface_cuda/zset_get.cpp
         )
IF(HAVE_PGI)
    SET(ZSRC ${ZSRC}
             ../control/magma_zf77pgi.cpp
             )
ELSE() 
    SET(ZSRCF ${ZSRCF}
              ../control/magma_zfortran.F90
              )
    SET(ZSRC ${ZSRC}
             ../control/magma_zf77.cpp
             )

    SET(MAGMA_CTRLF 
                    ../control/magma_param.f90
                    ../control/magma.f90
                    )
ENDIF()

precisions_rules_py(MAGMA_CTRL_GENERATED
                    "${ZSRC}"
                    PRECISIONS ${MAGMA_PRECISION})
precisions_rules_py(MAGMA_CTRLF_GENERATED
                    "${ZSRCF}"
                    PRECISIONS ${MAGMA_PRECISION})

# Define the list of sources
# --------------------------
SET(MAGMA_CTRLF ${MAGMA_CTRLF_GENERATED}
                ${MAGMA_CTRLF}
                )

SET(MAGMA_CTRL ${MAGMA_CTRLF}
               ${MAGMA_CTRL_GENERATED}
               ../control/auxiliary.cpp
               ../control/xerbla.cpp
               ../control/get_nb_tesla.cpp
               ../control/get_nb_fermi.cpp
               ../control/timer.cpp
               ../control/trace.cpp
               ../control/strlcpy.cpp
               ../control/bulge_auxiliary.cpp
               ../interface_cuda/alloc.cpp
               ../interface_cuda/error.cpp
               ../interface_cuda/set_get.cpp
               ../interface_cuda/interface.cpp
               )

# Generate the magma sources for all possible precisions
# ------------------------------------------------------
SET(MAGMA_SRCS "")
SET(MAGMA_SRCS_GENERATED "")
SET(ZSRC 
###
     zpotrf.cpp
     zpotrf_gpu.cpp
     zpotrf_mgpu.cpp
     zpotrf2_ooc.cpp
     zpotrf3_mgpu.cpp
     zpotrs_gpu.cpp
     zposv.cpp
     zposv_gpu.cpp
     zpotri.cpp
     zpotri_gpu.cpp
###
     zgetrf.cpp
     zgetrf_gpu.cpp
     zgetrf_mgpu.cpp
     zgetrf3_ooc.cpp
     zgetrf1_mgpu.cpp
     zgetf2_nopiv.cpp
     zgetrf_nopiv.cpp
     zgetrf_nopiv_gpu.cpp
     zgetrs_gpu.cpp
     zgesv.cpp
     zgesv_gpu.cpp
     zgetri_gpu.cpp
###
     zgeqrf.cpp
     zgeqrf-v4.cpp
     zgeqrf_ooc.cpp
     zgeqlf.cpp
     zgelqf.cpp
     zgelqf_gpu.cpp
     zungqr.cpp
     zungtr.cpp
     zgeqrf_gpu.cpp
     zgeqrf_mgpu-trace.cpp
     zgeqrf2_gpu.cpp
     zgeqrf3_gpu.cpp
     zgeqrs_gpu.cpp
     zgeqrs3_gpu.cpp
     zgels_gpu.cpp
     zgels3_gpu.cpp
     zlarfb_gpu.cpp
     zunmqr_gpu.cpp
     zungqr_gpu.cpp
###
     zgebrd.cpp
     zlabrd_gpu.cpp
     zgehrd.cpp
     zgehrd2.cpp
     zlahru.cpp
     zlahr2.cpp
     zunghr.cpp
     zhetrd.cpp
     zhetrd_gpu.cpp
     zhetrd2_gpu.cpp
     zhebbd.cpp
     zlatrd.cpp
     zlatrd2.cpp
     zunmtr.cpp
     zunmtr_gpu.cpp
     zunmqr.cpp
     zunmqr2_gpu.cpp
     zunmql.cpp
     zunmql2_gpu.cpp
     zunmtr_m.cpp
     zunmqr_m.cpp
###
     zgeev.cpp
     dgeev.cpp
     zheevd.cpp
     zheevdx.cpp
     zheevd_gpu.cpp
     zheevdx_gpu.cpp
     zheevr.cpp
     zheevr_gpu.cpp
     zheevx.cpp
     zheevx_gpu.cpp
     dsyevd.cpp
     dsyevdx.cpp
     dsyevd_gpu.cpp
     zhegvd.cpp
     zhegvdx.cpp
     zhegvx.cpp
     zhegvr.cpp
     dsygvd.cpp
     zhegst.cpp
     zhegst_gpu.cpp
     zgesvd.cpp
     dgesvd.cpp
     dstedx.cpp
     zstedx.cpp
     dlaex0.cpp
     dlaex1.cpp
     dlaex3.cpp
     zheevd_m.cpp
     zhegvd_m.cpp
     zhegst_m.cpp
     ztrsm_m.cpp
     dstedx_m.cpp
     zstedx_m.cpp
     dlaex0_m.cpp
     dlaex1_m.cpp
     dlaex3_m.cpp
###
     zcposv_gpu.cpp
     zcgesv_gpu.cpp
     zcgetrs_gpu.cpp
     zcgeqrsv_gpu.cpp
###
     ztrtri.cpp
     ztrtri_gpu.cpp
     zlauum.cpp
     zlauum_gpu.cpp
     zungqr_2stage_gpu.cpp
     zbulge_kernel.cpp
     zbulge_applyQ.cpp
     dbulge_aux.cpp
     zbulge_aux.cpp
     zhetrd_bhe2trc_v3.cpp
     zher2k_mgpu.cpp
     zhemm_mgpu.cpp
     )

IF(MAGMA_USE_PLASMA)
    SET(ZSRC ${ZSRC}
             zgessm_gpu.cpp
             zssssm_gpu.cpp
             zgetrl_gpu.cpp
             ztstrf_gpu.cpp
             )
ENDIF()
precisions_rules_py(MAGMA_SRCS_GENERATED
                    "${ZSRC}"
                    PRECISIONS ${MAGMA_PRECISION})

# Define the list of sources
# --------------------------
SET(MAGMA_SRCS 
               ${MAGMA_CTRL}
               ${MAGMA_SRCS_GENERATED}
               dlaed2.f
               slaed2.f
               )

# Add dependencies on the header files
# ------------------------------------
SET_SOURCE_FILES_PROPERTIES(${MAGMA_SRCS} 
    PROPERTIES OBJECTS_DEPENDS magma_include
    )

# Compile step
# ------------
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/control)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/interface_cuda)
ADD_LIBRARY(magma ${MAGMA_SRCS})
SET_PROPERTY(TARGET magma PROPERTY LINKER_LANGUAGE Fortran)

FOREACH(_external_package "blas" "lapack" "plasma")
    STRING(TOUPPER "${_external_package}" _NAMEVAR)
    IF(DEFINED ${_NAMEVAR}_BUILD_MODE)
        ADD_DEPENDENCIES(magma ${_external_package}_install)
    ENDIF()
ENDFOREACH()

# Install step
# ------------
INSTALL(TARGETS magma
        DESTINATION lib)

FOREACH(_src ${MAGMA_CTRLF})
    GET_FILENAME_COMPONENT(_namesrc ${_src} NAME)
    STRING(REGEX REPLACE ".[Ff]90" ".mod"  _namemod "${_namesrc}")
    INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${_namemod}
            DESTINATION include)
ENDFOREACH()

###
### END CMakeLists.txt
###
