###
#
# @file          : CMakeLists.txt
#
# @description   : Project MORSE
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 24-01-2012
# @last modified : mar. 22 mai 2012 10:30:47 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(RulesPrecisions)

# Generate the morse common for all possible precisions
# -----------------------------------------------------
SET(MORSE_COMMON_GENERATED "")
SET(ZSRC ../common/zprofiling.c
         ../common/zlocality.c
         )

precisions_rules_py(MORSE_COMMON_GENERATED
            "${ZSRC}"
            PRECISIONS ${MAGMA_PRECISION})

SET(MORSE_COMMON
    ../common/control.c
    ../common/context.c
    ../common/descriptor.c
    ../common/async.c
    ../common/options.c
    ../common/profiling.c
    ${MORSE_COMMON_GENERATED}
    )

# Generate the morse sources for all possible precisions
# ------------------------------------------------------
SET(MORSE_SRCS_GENERATED "")
SET(ZSRC codelet_zlacpy.c
         codelet_zplrnt.c codelet_zplghe.c codelet_zplgsy.c codelet_zunmlq.c
         codelet_zgemm.c  codelet_ztrsm.c  codelet_zherk.c  codelet_ztsmlq.c
         codelet_zpotrf.c codelet_zgelqt.c codelet_zlaset.c codelet_ztslqt.c
         codelet_zgetrl.c codelet_ztstrf.c codelet_zgessm.c codelet_zssssm.c
         codelet_zgeqrt.c codelet_ztsqrt.c codelet_zunmqr.c codelet_ztsmqr.c
         codelet_zaxpy.c  codelet_zttqrt.c codelet_zttlqt.c
         codelet_zttmlq.c codelet_zsymm.c codelet_zhemm.c
         codelet_zher2k.c codelet_zsyr2k.c
         codelet_ztrtri.c
         codelet_zlauum.c
         codelet_ztrmm.c
#        codelet_zbrdalg.c 
         )

precisions_rules_py(MORSE_SRCS_GENERATED
            "${ZSRC}"
            PRECISIONS ${MAGMA_PRECISION})

SET(MORSE_SRCS ${MORSE_COMMON}
               ${MORSE_SRCS_GENERATED}
               )

# Add dependencies on the header files
# ------------------------------------
SET_SOURCE_FILES_PROPERTIES(${MORSE_SRCS} 
    PROPERTIES OBJECTS_DEPENDS mgpu_include
    )
SET_SOURCE_FILES_PROPERTIES(${MORSE_SRCS} 
    PROPERTIES OBJECTS_DEPENDS common_include
    )
SET_SOURCE_FILES_PROPERTIES(${MORSE_SRCS} 
    PROPERTIES OBJECTS_DEPENDS morse_include
    )

# Compile step
# ------------
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../common)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../common)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../include)

ADD_LIBRARY(morse_quark ${MORSE_SRCS})
ADD_DEPENDENCIES(morse_quark magma_mgpu)
SET_PROPERTY(TARGET morse_quark PROPERTY LINKER_LANGUAGE Fortran)

FOREACH(_external_package "blas" "lapack" "plasma")
    STRING(TOUPPER "${_external_package}" _NAMEVAR)
    IF(DEFINED ${_NAMEVAR}_BUILD_MODE)
        ADD_DEPENDENCIES(morse_quark ${_external_package}_install)
    ENDIF()
ENDFOREACH()

# installation
# ------------
INSTALL(TARGETS morse_quark
        DESTINATION lib)

###
### END CMakeLists.txt
###
