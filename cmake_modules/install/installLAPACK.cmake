###
#
# @file      : installLAPACK.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 21-01-2012
# @last modified : mar. 05 juin 2012 16:12:15 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(definePACKAGE)
INCLUDE(downloadPACKAGE)
INCLUDE(installPACKAGE)
INCLUDE(infoLAPACK)

MACRO(INSTALL_LAPACK _MODE)

    # Get info for this package
    # -------------------------
    LAPACK_INFO_INSTALL()

    # Search for dependencies
    # -----------------------
    DEFINE_PACKAGE("BLAS" "depends")

    # Define prefix paths
    # -------------------
    IF(MORSE_SEPARATE_PROJECTS)
        SET(LAPACK_PATH ${CMAKE_INSTALL_PREFIX}/lapack)
    ELSE(MORSE_SEPARATE_PROJECTS)
        SET(LAPACK_PATH ${CMAKE_INSTALL_PREFIX})
    ENDIF(MORSE_SEPARATE_PROJECTS)

    # Create make.inc
    # ---------------
    IF(EXISTS ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc)
        FILE(REMOVE ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc)
    ENDIF()
    STRING(TOUPPER "${CMAKE_BUILD_TYPE}" TOUPPER_BUILD_TYPE)
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "SHELL    = /bin/sh\n")

    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "FORTRAN  = ${CMAKE_Fortran_COMPILER}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "OPTS     = ${CMAKE_Fortran_COMPILER_${TOUPPER_BUILD_TYPE}} ${CMAKE_EXTRA_FFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "DRVOPTS  = ${CMAKE_Fortran_COMPILER_${TOUPPER_BUILD_TYPE}} ${CMAKE_EXTRA_FFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "NOOPT    = ${CMAKE_EXTRA_NOOPT}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "LOADER   = ${CMAKE_Fortran_COMPILER}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "LOADOPTS = ${CMAKE_Fortran_COMPILER_${TOUPPER_BUILD_TYPE}} ${CMAKE_EXTRA_LDFLAGS_F}\n")

    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "TIMER = INT_CPU_TIME\n")

    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "ARCH      = ${CMAKE_AR}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "ARCHFLAGS = cr\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "RANLIB    = ${CMAKE_RANLIB}\n")

    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "BLASLIB   = ${BLAS_LDFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "XBLASLIB  = \n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "LAPACKLIB = liblapack.a\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "TMGLIB    = libtmg.a\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "EIGSRCLIB = libeigsrc.a\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc "LINSRCLIB = liblinsrc.a\n")

    # Define steps of installation
    # ----------------------------
    SET(LAPACK_CONFIG_CMD ${CMAKE_COMMAND} -E copy
                        ${CMAKE_SOURCE_DIR}/externals/lapack_make.inc
                        ${CMAKE_BINARY_DIR}/externals/lapack/make.inc)
    SET(LAPACK_MAKE_CMD ${CMAKE_MAKE_PROGRAM} lapacklib)
    SET(LAPACK_MAKEINSTALL_CMD ${CMAKE_COMMAND} -E copy
                             ${CMAKE_BINARY_DIR}/externals/lapack/liblapack.a
                             ${LAPACK_PATH}/lib/liblapack.a)

    # Define additional step
    # ----------------------
    SET(LAPACK_ADD_STEP lapack_create_prefix)
    SET(lapack_create_prefix_CMD ${CMAKE_COMMAND} -E make_directory ${LAPACK_PATH}/lib)
    SET(lapack_create_prefix_DIR ${CMAKE_INSTALL_PREFIX})
    SET(lapack_create_prefix_DEP_BEFORE build)
    SET(lapack_create_prefix_DEP_AFTER install)

    # Define options
    # --------------
    SET(LAPACK_OPTIONS "")

    # Install the external package
    # ----------------------------
    DEFINE_DOWNLOAD_PACKAGE("lapack" "${_MODE}")
    INSTALL_PACKAGE("lapack" "${LAPACK_BUILD_MODE}")

    # Set linker flags
    # ----------------
    STRING(TOUPPER "${MORSE_USE_BLAS}" VALUE_MORSE_USE_BLAS)
    IF("${VALUE_MORSE_USE_BLAS}" MATCHES "EIGEN")
        SET(LAPACK_LIBRARY_PATH ${EIGENLAPACK_LIBRARY_PATH} ${LAPACK_PATH}/lib)
        LIST(REMOVE_DUPLICATES LAPACK_LIBRARY_PATH)
        SET(LAPACK_LDFLAGS "${EIGENLAPACK_LDFLAGS} -L${LAPACK_PATH}/lib -llapack")
        SET(LAPACK_LIBRARIES "${EIGENLAPACK_LIBRARIES};lapack")
    ELSE()
        SET(LAPACK_LIBRARY_PATH ${LAPACK_PATH}/lib)
        SET(LAPACK_LDFLAGS "-L${LAPACK_LIBRARY_PATH} -llapack")
        SET(LAPACK_LIBRARIES "lapack")
    ENDIF()

ENDMACRO(INSTALL_LAPACK)

###
### END installLAPACK.cmake
###
