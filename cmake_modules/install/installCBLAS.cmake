###
#
# @file      : installCBLAS.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 21-01-2012
# @last modified : mer. 16 mai 2012 10:16:06 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(installPACKAGE)
INCLUDE(downloadPACKAGE)
INCLUDE(infoCBLAS)

MACRO(INSTALL_CBLAS _MODE)

    # Get info for this package
    # -------------------------
    CBLAS_INFO_INSTALL()

    # Define prefix paths
    # -------------------
    IF(MORSE_SEPARATE_PROJECTS)
        SET(CBLAS_PATH ${CMAKE_INSTALL_PREFIX}/cblas)
    ELSE(MORSE_SEPARATE_PROJECTS)
        SET(CBLAS_PATH ${CMAKE_INSTALL_PREFIX})
    ENDIF(MORSE_SEPARATE_PROJECTS)

    # Create make.inc
    # ---------------
    IF(EXISTS ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc)
        FILE(REMOVE ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc)
    ENDIF()
    STRING(TOUPPER "${CMAKE_BUILD_TYPE}" TOUPPER_BUILD_TYPE)
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "SHELL  = /bin/sh\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "CC     = ${CMAKE_C_COMPILER}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "FC     = ${CMAKE_Fortran_COMPILER}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "LOADER = ${CMAKE_Fortran_COMPILER} ${CMAKE_EXTRA_LDFLAGS_F}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "CFLAGS = ${CMAKE_C_COMPILER_${TOUPPER_BUILD_TYPE}} ${CMAKE_EXTRA_CFLAGS} ${FORTRAN_MANGLING_DETECTED}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "FFLAGS = ${CMAKE_Fortran_COMPILER_${TOUPPER_BUILD_TYPE}} ${CMAKE_EXTRA_FFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "ARCH   = ${CMAKE_AR}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "ARCHFLAGS = cr\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "RANLIB = ${CMAKE_RANLIB}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "BLLIB  = ${BLAS_LDFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc "CBLIB  = libcblas.a\n")

    # Define steps of installation
    # ----------------------------
    SET(CBLAS_CONFIG_CMD ${CMAKE_COMMAND} -E copy
                        ${CMAKE_SOURCE_DIR}/externals/cblas_make.inc
                        ${CMAKE_BINARY_DIR}/externals/cblas/Makefile.in)
    SET(CBLAS_MAKE_CMD ${CMAKE_MAKE_PROGRAM} alllib)
    SET(CBLAS_MAKEINSTALL_CMD ${CMAKE_COMMAND} -E copy
                             ${CMAKE_BINARY_DIR}/externals/cblas/src/libcblas.a
                             ${CBLAS_PATH}/lib/libcblas.a)

    # Define additional step
    # ----------------------
    SET(CBLAS_ADD_STEP cblas_create_prefix_lib cblas_create_prefix_include cblas_copy_include)
    SET(cblas_create_prefix_lib_CMD ${CMAKE_COMMAND} -E make_directory ${CBLAS_PATH}/lib)
    SET(cblas_create_prefix_lib_DIR ${CMAKE_INSTALL_PREFIX})
    SET(cblas_create_prefix_lib_DEP_BEFORE build)
    SET(cblas_create_prefix_lib_DEP_AFTER install)
    SET(cblas_create_prefix_include_CMD ${CMAKE_COMMAND} -E make_directory ${CBLAS_PATH}/include)
    SET(cblas_create_prefix_include_DIR ${CMAKE_INSTALL_PREFIX})
    SET(cblas_create_prefix_include_DEP_BEFORE build)
    SET(cblas_create_prefix_include_DEP_AFTER install)
    SET(cblas_copy_include_CMD ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/externals/cblas/include/cblas.h .)
    SET(cblas_copy_include_DIR ${CBLAS_PATH}/include)
    SET(cblas_copy_include_DEP_BEFORE cblas_create_prefix_include)
    SET(cblas_copy_include_DEP_AFTER install)

    # Define options
    # --------------
    SET(CBLAS_OPTIONS "")

    # Install the external package
    # ----------------------------
    DEFINE_DOWNLOAD_PACKAGE("cblas" "${_MODE}")
    INSTALL_PACKAGE("cblas" "${CBLAS_BUILD_MODE}")

    # Set linker flags
    # ----------------
    SET(CBLAS_LIBRARY_PATH ${CBLAS_PATH}/lib)
    SET(CBLAS_INCLUDE_PATH ${CBLAS_PATH}/include)
    SET(CBLAS_LDFLAGS "-L${CBLAS_LIBRARY_PATH} -lcblas")
    SET(CBLAS_LIBRARIES "cblas")

ENDMACRO(INSTALL_CBLAS)

###
### END installCBLAS.cmake
###
