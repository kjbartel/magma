###
#
# @file      : installPLASMA.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 21-01-2012
# @last modified : mar. 05 juin 2012 15:37:33 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(definePACKAGE)
INCLUDE(downloadPACKAGE)
INCLUDE(installPACKAGE)
INCLUDE(infoPLASMA)

MACRO(INSTALL_PLASMA _MODE)

    # Get info for this package
    # -------------------------
    PLASMA_INFO_INSTALL()

    # Search for dependencies
    # -----------------------
    DEFINE_PACKAGE("LAPACK" "depends")
    DEFINE_PACKAGE("CBLAS" "depends")
    DEFINE_PACKAGE("LAPACKE" "depends")

    # Define prefix paths
    # -------------------
    IF(MORSE_SEPARATE_PROJECTS)
        SET(PLASMA_PATH ${CMAKE_INSTALL_PREFIX}/plasma)
    ELSE(MORSE_SEPARATE_PROJECTS)
        SET(PLASMA_PATH ${CMAKE_INSTALL_PREFIX})
    ENDIF(MORSE_SEPARATE_PROJECTS)

    # Create make.inc
    # ---------------
    IF(EXISTS ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc)
        FILE(REMOVE ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc)
    ENDIF()
    STRING(TOUPPER "${CMAKE_BUILD_TYPE}" TOUPPER_BUILD_TYPE)
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "prefix     = ${PLASMA_PATH}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "CC         = ${CMAKE_C_COMPILER}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "FC         = ${CMAKE_Fortran_COMPILER}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "LOADER     = ${CMAKE_Fortran_COMPILER}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "ARCH       = ${CMAKE_AR}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "ARCHFLAGS  = cr\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "RANLIB     = ${CMAKE_RANLIB}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "CFLAGS     = ${CMAKE_C_COMPILER_${TOUPPER_BUILD_TYPE}} ${CMAKE_EXTRA_CFLAGS} ${FORTRAN_MANGLING_DETECTED}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "FFLAGS     = ${CMAKE_Fortran_COMPILER_${TOUPPER_BUILD_TYPE}} ${CMAKE_EXTRA_FFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "LDFLAGS    = ${CMAKE_EXTRA_LDFLAGS_F} ${FORTRAN_MANGLING_DETECTED}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "LIBBLAS    = ${BLAS_LDFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "LIBCBLAS   = ${CBLAS_LDFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "LIBLAPACK  = ${LAPACK_LDFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "INCCLAPACK = -I${LAPACKE_INCLUDE_PATH}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc "LIBCLAPACK = ${LAPACKE_LDFLAGS}\n")

    # Define steps of installation
    # ----------------------------
    SET(PLASMA_CONFIG_CMD ${CMAKE_COMMAND} -E copy
                        ${CMAKE_SOURCE_DIR}/externals/plasma_make.inc
                        ${CMAKE_BINARY_DIR}/externals/plasma/make.inc)
    SET(PLASMA_MAKE_CMD ${CMAKE_MAKE_PROGRAM} lib)
    SET(PLASMA_MAKEINSTALL_CMD ${CMAKE_MAKE_PROGRAM} install)

    # Define options
    # --------------
    SET(PLASMA_OPTIONS "") 

    # Install the external package
    # ----------------------------
    DEFINE_DOWNLOAD_PACKAGE("plasma" "${_MODE}")
    INSTALL_PACKAGE("plasma" "${PLASMA_BUILD_MODE}")

    # Set linker flags
    # ----------------
    SET(PLASMA_LIBRARY_PATH ${PLASMA_PATH}/lib)
    SET(PLASMA_INCLUDE_PATH ${PLASMA_PATH}/include)
    SET(PLASMA_LDFLAGS "-L${PLASMA_LIBRARY_PATH} -lplasma -lcoreblas -lquark -l${LAPACKE_LIBRARIES} -l${CBLAS_LIBRARIES} ${LAPACK_LDFLAGS} ${BLAS_LDFLAGS}")
    SET(PLASMA_LIBRARIES "plasma;coreblas;quark;${LAPACKE_LIBRARIES};${CBLAS_LIBRARIES};${LAPACK_LIBRARIES};${BLAS_LIBRARIES}")

ENDMACRO(INSTALL_PLASMA)

###
### END installPLASMA.cmake
###
