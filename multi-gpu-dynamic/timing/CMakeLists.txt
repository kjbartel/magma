###
#
# @file      : CMakeLists.txt
#
# @description   : Project MORSE
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 23-01-2012
# @last modified : mar. 22 mai 2012 09:39:39 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(RulesPrecisions)

# Generate the morse auxiliary headers for all possible precisions
# ----------------------------------------------------------------
SET(MORSE_AUX_HDRS_GENERATED "")
SET(ZHDR zauxiliary.h)
precisions_rules_py(MORSE_AUX_HDRS_GENERATED
            "${ZHDR}"
            PRECISIONS ${MAGMA_PRECISION})
SET(MORSE_AUX_HDRS auxiliary.h
                   timing.h
                   timing.c
                   ${MORSE_AUX_HDRS_GENERATED}
                   )

# Force generation of headers
# ---------------------------
ADD_CUSTOM_TARGET(morse_auxiliary_headers ALL SOURCES ${MORSE_AUX_HDRS})

# Generate the morse auxiliary sources for all possible precisions
# ----------------------------------------------------------------
SET(MORSE_AUX_SRCS_GENERATED "")
SET(ZSRC zauxiliary.c)
precisions_rules_py(MORSE_AUX_SRCS_GENERATED
            "${ZSRC}"
            PRECISIONS ${MAGMA_PRECISION})
SET(MORSE_SRCS auxiliary.c
               ${MORSE_AUX_SRCS_GENERATED}
               )

# Force generation of headers
# ---------------------------
ADD_CUSTOM_TARGET(morse_auxiliary_sources ALL SOURCES ${MORSE_AUX_SRCS_GENERATED})

# Generate the morse testing sources for all possible precisions
# --------------------------------------------------------------
SET(MORSE_TESTS_GENERATED "")
SET(ZTEST 
      time_zgemm_tile.c
      time_zgeqrf_tile.c
      time_zgesv_incpiv.c
      time_zgesv_incpiv_tile.c
      time_zgetrf_incpiv.c
      time_zgetrf_incpiv_tile.c
      time_zposv_tile.c
      time_zpotrf_tile.c
      time_zpotri_tile.c
      ####
      #time_zgemm.c
      #time_ztrsm.c
      )

precisions_rules_py(MORSE_TESTS_GENERATED
            "${ZTEST}"
            PRECISIONS ${MAGMA_PRECISION})

SET(MORSE_TESTS ${MORSE_TESTS_GENERATED})

# Add dependencies on the header files
# ------------------------------------
#SET_SOURCE_FILES_PROPERTIES(${MORSE_TESTS} 
#    PROPERTIES OBJECTS_DEPENDS morse_auxiliary_headers
#    )
#SET_SOURCE_FILES_PROPERTIES(${MORSE_TESTS} 
#    PROPERTIES OBJECTS_DEPENDS morse_auxiliary_sources
#    )

# Add include and link directories
# --------------------------------
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)

# Create libauxilary.a
# --------------------
ADD_LIBRARY(auxilary STATIC ${MORSE_SRCS})
SET_PROPERTY(TARGET auxilary PROPERTY LINKER_LANGUAGE Fortran)
ADD_DEPENDENCIES(auxilary morse_auxiliary_headers)
ADD_DEPENDENCIES(auxilary morse_auxiliary_sources)

FOREACH(_external_package "blas" "lapack" "starpu" "plasma")
    STRING(TOUPPER "${_external_package}" _NAMEVAR)
    IF(DEFINED ${_NAMEVAR}_BUILD_MODE)
        ADD_DEPENDENCIES(auxilary ${_external_package}_install)
    ENDIF()
ENDFOREACH()

# Add tests (C/CPP)
# -----------------
SET(libs_for_tests "${CMAKE_EXTRA_LDFLAGS_F}")

IF(MORSE_SCHED_STARPU)
    SET(libs_for_tests ${libs_for_tests} morse_starpu ${STARPU_LIBRARIES})
ELSEIF( MORSE_SCHED_QUARK )
    SET(libs_for_tests ${libs_for_tests} morse_quark ${QUARK_LIBRARIES})
ENDIF(MORSE_SCHED_STARPU)

IF(MORSE_USE_CUDA)
    SET(libs_for_tests ${libs_for_tests} magma magmablas magma cublas)
ENDIF(MORSE_USE_CUDA)

IF(MORSE_USE_MPI)
    SET(libs_for_tests ${libs_for_tests} ${MPI_LIBRARIES})
ENDIF(MORSE_USE_MPI)

SET(libs_for_tests ${libs_for_tests} ${PLASMA_LIBRARIES})
SET(libs_for_tests auxilary magma_mgpu ${libs_for_tests} ${EXTRA_LIBS})

FOREACH(_test_src ${MORSE_TESTS_GENERATED})
    GET_FILENAME_COMPONENT(_name_exe ${_test_src} NAME_WE)
    ADD_EXECUTABLE(${_name_exe} ${_test_src})
    SET_PROPERTY(TARGET ${_name_exe} PROPERTY LINKER_LANGUAGE Fortran)
    TARGET_LINK_LIBRARIES(${_name_exe} ${libs_for_tests})
    ADD_TEST(NAME ${_name_exe} 
             COMMAND ${_name_exe} --n_range=1000:2000:1000 --check)
ENDFOREACH()
