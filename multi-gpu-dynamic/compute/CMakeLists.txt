###
#
# @file          : CMakeLists.txt
#
# @description   : Project MORSE
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 24-01-2012
# @last modified : mar. 22 mai 2012 10:30:14 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(RulesPrecisions)

# Define the list of sources
# --------------------------
SET(MORSE_COMMON ../common/context.c
                 ../common/control.c
                 ../common/descriptor.c
                 ../common/async.c
                 ../common/workspace.c
                 )

# Generate the morse sources for all possible precisions
# ------------------------------------------------------
SET(MORSE_SRCS_GENERATED "")
SET(ZSRC 
    zplrnt.c zplghe.c zplgsy.c ztile.c
    zgemm.c ztrsm.c ztrsmpl.c
    zpotrf.c zpotrs.c zposv.c
    zgetrf_incpiv.c zgetrs_incpiv.c zgesv_incpiv.c
    zgeqrf.c
    zlauum.c ztrtri.c zpotri.c
    ###
    pzplrnt.c pzplghe.c pzplgsy.c pztile.c
    pzgemm.c pztrsm.c pztrsmpl.c
    pzpotrf.c pzgetrf_incpiv.c pzgeqrf.c
    pztrtri.c pzlauum.c
    ###
    ../common/workspace_z.c
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
    PROPERTIES OBJECTS_DEPENDS morse_include
    )
SET_SOURCE_FILES_PROPERTIES(${MORSE_SRCS} 
    PROPERTIES OBJECTS_DEPENDS morse_common_include
    )

# Compile step
# ------------
ADD_LIBRARY(magma_mgpu ${MORSE_SRCS})
IF(MAGMA_1GPU)
	ADD_DEPENDENCIES(magma_mgpu magma)
	ADD_DEPENDENCIES(magma_mgpu magmablas)
ENDIF(MAGMA_1GPU)
SET_PROPERTY(TARGET magma_mgpu PROPERTY LINKER_LANGUAGE Fortran)

FOREACH(_external_package "blas" "lapack" "starpu" "plasma")
    STRING(TOUPPER "${_external_package}" _NAMEVAR)
    IF(DEFINED ${_NAMEVAR}_BUILD_MODE)
        ADD_DEPENDENCIES(magma_mgpu ${_external_package}_install)
    ENDIF()
ENDFOREACH()

# installation
# ------------
INSTALL(TARGETS magma_mgpu
        DESTINATION lib)

###
### END CMakeLists.txt
###
