###
#
# @file      : installLAPACKE.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 21-01-2012
# @last modified : mer. 16 mai 2012 10:18:29 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(installPACKAGE)
INCLUDE(downloadPACKAGE)
INCLUDE(infoLAPACKE)

MACRO(INSTALL_LAPACKE _MODE)

    # Get info for this package
    # -------------------------
    LAPACKE_INFO_INSTALL()

    # Define prefix paths
    # -------------------
    IF(MORSE_SEPARATE_PROJECTS)
        SET(LAPACKE_PATH ${CMAKE_INSTALL_PREFIX}/lapacke)
    ELSE(MORSE_SEPARATE_PROJECTS)
        SET(LAPACKE_PATH ${CMAKE_INSTALL_PREFIX})
    ENDIF(MORSE_SEPARATE_PROJECTS)

    # Create make.inc
    # ---------------
    IF(EXISTS ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc)
        FILE(REMOVE ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc)
    ENDIF()
    STRING(TOUPPER "${CMAKE_BUILD_TYPE}" TOUPPER_BUILD_TYPE)
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc "SHELL   = /bin/sh\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc "CC      = ${CMAKE_C_COMPILER} ${CMAKE_C_COMPILER_${TOUPPER_BUILD_TYPE}} ${CMAKE_EXTRA_CFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc "FORTRAN = ${CMAKE_Fortran_COMPILER} ${CMAKE_Fortran_COMPILER_${TOUPPER_BUILD_TYPE}} ${CMAKE_EXTRA_FFLAGS}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc "LOADER  = ${CMAKE_Fortran_COMPILER} ${CMAKE_EXTRA_LDFLAGS_F}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc "OPTS    = \n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc "NOOPT   = ${CMAKE_EXTRA_NOOPT}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc "ARCH    = ${CMAKE_AR}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc "ARCHFLAGS = cr\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc "RANLIB  = ${CMAKE_RANLIB}\n")
    FILE(APPEND ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc "LAPACKE = liblapacke.a\n")

    # Define steps of installation
    # ----------------------------
    SET(LAPACKE_CONFIG_CMD ${CMAKE_COMMAND} -E copy
                        ${CMAKE_SOURCE_DIR}/externals/lapacke_make.inc
                        ${CMAKE_BINARY_DIR}/externals/lapacke/make.inc)
    SET(LAPACKE_MAKE_CMD ${CMAKE_MAKE_PROGRAM} lapacke)
    SET(LAPACKE_MAKEINSTALL_CMD ${CMAKE_COMMAND} -E copy
                             ${CMAKE_BINARY_DIR}/externals/lapacke/liblapacke.a
                             ${LAPACKE_PATH}/lib/liblapacke.a)

    # Define additional step
    # ----------------------
    SET(LAPACKE_ADD_STEP lapacke_create_prefix_lib lapacke_create_prefix_include lapacke_copy_include)
    SET(lapacke_create_prefix_lib_CMD ${CMAKE_COMMAND} -E make_directory ${LAPACKE_PATH}/lib)
    SET(lapacke_create_prefix_lib_DIR ${CMAKE_INSTALL_PREFIX})
    SET(lapacke_create_prefix_lib_DEP_BEFORE build)
    SET(lapacke_create_prefix_lib_DEP_AFTER install)
    SET(lapacke_create_prefix_include_CMD ${CMAKE_COMMAND} -E make_directory ${LAPACKE_PATH}/include)
    SET(lapacke_create_prefix_include_DIR ${CMAKE_INSTALL_PREFIX})
    SET(lapacke_create_prefix_include_DEP_BEFORE build)
    SET(lapacke_create_prefix_include_DEP_AFTER install)
    SET(lapacke_copy_include_CMD ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/externals/lapacke/include/lapacke.h .)
    SET(lapacke_copy_include_DIR ${LAPACKE_PATH}/include)
    SET(lapacke_copy_include_DEP_BEFORE lapacke_create_prefix_include)
    SET(lapacke_copy_include_DEP_AFTER install)

    # Define options
    # --------------
    SET(LAPACKE_OPTIONS "")

    # Install the external package
    # ----------------------------
    DEFINE_DOWNLOAD_PACKAGE("lapacke" "${_MODE}")
    INSTALL_PACKAGE("lapacke" "${LAPACKE_BUILD_MODE}")

    # Set linker flags
    # ----------------
    SET(LAPACKE_LIBRARY_PATH "${LAPACKE_PATH}/lib")
    SET(LAPACKE_INCLUDE_PATH "${LAPACKE_PATH}/include")
    SET(LAPACKE_LDFLAGS "-L${LAPACKE_LIBRARY_PATH} -llapacke")
    SET(LAPACKE_LIBRARIES "lapacke")

ENDMACRO(INSTALL_LAPACKE)

###
### END installLAPACKE.cmake
###
