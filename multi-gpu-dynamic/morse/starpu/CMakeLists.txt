###
#
# @file          : CMakeLists.txt
#
# @description   : Project MORSE
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 25-01-2012
# @last modified : mar. 22 mai 2012 10:31:00 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(RulesPrecisions)

# Add StarPU headers
#-------------------
INCLUDE_DIRECTORIES(${STARPU_INCLUDE_PATH})

# Generate the magma headers for all possible precisions
# ------------------------------------------------------
SET(MORSE_HEADERS_GENERATED "")
SET(ZHDR include/codelet_z.h)

precisions_rules_py(MORSE_HEADERS_GENERATED
            "${ZHDR}"
            PRECISIONS "z;c;d;s;zc;ds")

# Define the list of headers
# --------------------------
SET(MORSE_HEADERS 
    include/morse_starpu.h
    include/codelets.h
    ${MORSE_HEADERS_GENERATED}
    )

# Force generation of headers
# ---------------------------
ADD_CUSTOM_TARGET(morse_starpu_codelets ALL SOURCES ${MORSE_HEADERS_GENERATED})

# Generate the morse common for all possible precisions
# -----------------------------------------------------
SET(MORSE_COMMON_GENERATED "")
SET(ZSRC 
    common/zprofiling.c
    common/zlocality.c
    )

precisions_rules_py(MORSE_COMMON_GENERATED
            "${ZSRC}"
            PRECISIONS ${MAGMA_PRECISION})

SET(MORSE_COMMON
    common/control.c
    common/context.c
    common/descriptor.c
    common/async.c
    common/options.c
    common/profiling.c
    common/workspace.c
    ${MORSE_COMMON_GENERATED}
    )

# Generate the morse sources for all possible precisions
# ------------------------------------------------------
SET(MORSE_SRCS_GENERATED "")
SET(ZSRC 
    codelets/codelet_zlacpy.c
    codelets/codelet_zplrnt.c 
    codelets/codelet_zplghe.c 
    codelets/codelet_zplgsy.c
    codelets/codelet_zgemm.c  
    codelets/codelet_ztrsm.c  
    codelets/codelet_zherk.c
    codelets/codelet_zpotrf.c
    codelets/codelet_zgetrl.c 
    codelets/codelet_ztstrf.c 
    codelets/codelet_zgessm.c 
    codelets/codelet_zssssm.c
    codelets/codelet_zgeqrt.c 
    codelets/codelet_ztsqrt.c 
    codelets/codelet_zunmqr.c 
    codelets/codelet_ztsmqr.c
    codelets/codelet_ztrtri.c
    codelets/codelet_zlauum.c
    codelets/codelet_ztrmm.c
    codelets/codelet_zcallback.c
    )

precisions_rules_py(MORSE_SRCS_GENERATED
            "${ZSRC}"
            PRECISIONS ${MAGMA_PRECISION})

SET(MORSE_SRCS ${MORSE_COMMON}
               ${MORSE_SRCS_GENERATED}
               )

# Include directories
# -------------------
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# Add library
# -----------
ADD_LIBRARY(morse_starpu ${MORSE_SRCS})
ADD_DEPENDENCIES(morse_starpu magma_mgpu)
SET_PROPERTY(TARGET morse_starpu PROPERTY LINKER_LANGUAGE Fortran)

ADD_DEPENDENCIES(morse_starpu morse_include)
ADD_DEPENDENCIES(morse_starpu morse_common_include)
ADD_DEPENDENCIES(morse_starpu morse_starpu_codelets)

FOREACH(_external_package "blas" "lapack" "starpu" "plasma")
    STRING(TOUPPER "${_external_package}" _NAMEVAR)
    IF(DEFINED ${_NAMEVAR}_BUILD_MODE)
        ADD_DEPENDENCIES(morse_starpu ${_external_package}_install)
    ENDIF()
ENDFOREACH()

# installation
# ------------
INSTALL(TARGETS morse_starpu 
        DESTINATION lib)

#INSTALL(FILES ${MORSE_HEADERS} 
#        DESTINATION include)

###
### END CMakeLists.txt
###
