###
#
# @file      : installREFBLAS.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 21-01-2012
# @last modified : mar. 12 juin 2012 16:44:31 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(installPACKAGE)
INCLUDE(downloadPACKAGE)

MACRO(INSTALL_REFBLAS _MODE)

    # Define prefix paths
    # -------------------
    IF(MORSE_SEPARATE_PROJECTS)
        SET(BLAS_PATH ${CMAKE_INSTALL_PREFIX}/blas)
    ELSE(MORSE_SEPARATE_PROJECTS)
        SET(BLAS_PATH ${CMAKE_INSTALL_PREFIX})
    ENDIF(MORSE_SEPARATE_PROJECTS)

    # Create make.inc
    # ---------------
    IF(EXISTS ${CMAKE_SOURCE_DIR}/externals/blas_make.inc)
        FILE(REMOVE ${CMAKE_SOURCE_DIR}/externals/blas_make.inc)
    ENDIF()
    STRING(TOUPPER "${CMAKE_BUILD_TYPE}" TOUPPER_BUILD_TYPE)
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/blas_make.inc "FORTRAN = ${CMAKE_Fortran_COMPILER}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/blas_make.inc "OPTS = ${CMAKE_Fortran_COMPILER_${TOUPPER_BUILD_TYPE}} ${CMAKE_EXTRA_FFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/blas_make.inc "DRVOPTS = $(OPTS)\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/blas_make.inc "LOADER = ${CMAKE_Fortran_COMPILER} ${CMAKE_EXTRA_LDFLAGS_F}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/blas_make.inc "ARCH = ${CMAKE_AR}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/blas_make.inc "ARCHFLAGS = cr\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/blas_make.inc "RANLIB = ${CMAKE_RANLIB}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/blas_make.inc "BLASLIB = libblas.a\n")

    # Define steps of installation
    # ----------------------------
    SET(BLAS_CONFIG_CMD ${CMAKE_COMMAND} -E copy
                        ${CMAKE_SOURCE_DIR}/externals/blas_make.inc
                        ${CMAKE_BINARY_DIR}/externals/blas/make.inc)
    SET(BLAS_MAKE_CMD ${CMAKE_MAKE_PROGRAM})
    SET(BLAS_MAKEINSTALL_CMD ${CMAKE_COMMAND} -E copy
                             ${CMAKE_BINARY_DIR}/externals/blas/libblas.a
                             ${BLAS_PATH}/lib/libblas.a)

    # Define additional step
    # ----------------------
    SET(BLAS_ADD_STEP blas_create_prefix)
    SET(blas_create_prefix_CMD ${CMAKE_COMMAND} -E make_directory ${BLAS_PATH}/lib)
    SET(blas_create_prefix_DIR ${CMAKE_INSTALL_PREFIX})
    SET(blas_create_prefix_DEP_BEFORE build)
    SET(blas_create_prefix_DEP_AFTER install)

    # Define options
    # --------------
    SET(BLAS_OPTIONS "")

    # Install the external package
    # ----------------------------
    DEFINE_DOWNLOAD_PACKAGE("blas" "${_MODE}")
    INSTALL_PACKAGE("blas" "${BLAS_BUILD_MODE}")
    MESSAGE(STATUS "Installing BLAS - refblas version")
    MESSAGE(STATUS "Installing BLAS - do not expect high performance from this reference library!")

    # Set linker flags
    # ----------------
    SET(BLAS_LIBRARY_PATH ${BLAS_PATH}/lib)
    SET(BLAS_LDFLAGS "-L${BLAS_LIBRARY_PATH} -lblas")
    SET(BLAS_LIBRARIES "blas")

ENDMACRO(INSTALL_REFBLAS)

###
### END installREFBLAS.cmake
###
