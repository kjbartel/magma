###
#
#  @file infoPACKAGE.cmake
#
#  @project MAGMA
#  MAGMA is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.1.0
#  @author Cedric Castagnede
#  @date 13-07-2012
#
###

###
#
# INFO_INSTALL_PACKAGE: Macro to define test functions
#
###
MACRO(INFO_DEPS_PACKAGE _NAME)

    # Define input variable to uppercase
    # ----------------------------------
    STRING(TOUPPER "${_NAME}" _UPIDP)

    # Looking for info for MAGMA_ALL_PACKAGES
    # ---------------------------------------
    IF("^BLAS$" STREQUAL "^${_UPIDP}$")
        INCLUDE(infoBLAS)
        BLAS_INFO_DEPS()

    ELSEIF("^TMG$" STREQUAL "^${_UPIDP}$")
        INCLUDE(infoTMG)
        TMG_INFO_DEPS()

    ELSEIF("^LAPACK$" STREQUAL "^${_UPIDP}$")
        INCLUDE(infoLAPACK)
        LAPACK_INFO_DEPS()

    ELSEIF("^CBLAS$" STREQUAL "^${_UPIDP}$")
        INCLUDE(infoCBLAS)
        CBLAS_INFO_DEPS()

    ELSEIF("^PLASMA$" STREQUAL "^${_UPIDP}$")
        INCLUDE(infoPLASMA)
        PLASMA_INFO_DEPS()

    ELSEIF("^CUDA$" STREQUAL "^${_UPIDP}$")
        INCLUDE(infoCUDA)
        CUDA_INFO_DEPS()

    ELSE()
        MESSAGE(FATAL_ERROR "DEV_ERROR - Macro INFO_DEPS_PACKAGE: ${_UPIDP}")

    ENDIF()

ENDMACRO(INFO_DEPS_PACKAGE)

###
#
# INFO_FIND_PACKAGE: Macro to define infos
#
###
MACRO(INFO_FIND_PACKAGE _NAME)

    # Define input variable to uppercase
    # ----------------------------------
    STRING(TOUPPER "${_NAME}" _UPIFP)

   # Looking for info for MAGMA_ALL_PACKAGES
    # --------------------------------------
    IF("^BLAS$" STREQUAL "^${_UPIFP}$")
        INCLUDE(infoBLAS)
        BLAS_INFO_FIND()

    ELSEIF("^TMG$" STREQUAL "^${_UPIFP}$")
        INCLUDE(infoTMG)
        TMG_INFO_FIND()

    ELSEIF("^LAPACK$" STREQUAL "^${_UPIFP}$")
        INCLUDE(infoLAPACK)
        LAPACK_INFO_FIND()

    ELSEIF("^CBLAS$" STREQUAL "^${_UPIFP}$")
        INCLUDE(infoCBLAS)
        CBLAS_INFO_FIND()

    ELSEIF("^PLASMA$" STREQUAL "^${_UPIFP}$")
        INCLUDE(infoPLASMA)
        PLASMA_INFO_FIND()

    ELSE()
        MESSAGE(FATAL_ERROR "DEV_ERROR - Macro INFO_FIND_PACKAGE: ${_UPIFP}")

    ENDIF()

ENDMACRO(INFO_FIND_PACKAGE)

##
## @end file infoPACKAGE.cmake
##
