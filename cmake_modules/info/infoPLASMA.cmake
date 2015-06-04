###
#
#  @file infoPLASMA.cmake
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

MACRO(PLASMA_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(PLASMA_type_library        "C;Fortran"                          )
    SET(PLASMA_name_library        "plasma;coreblas;quark"              )
    SET(PLASMA_name_pkgconfig      "plasma"                             )
    SET(PLASMA_name_include        "plasma.h;core_blas.h;quark.h"       )
    SET(PLASMA_name_include_suffix "PLASMA_name_include_suffix-NOTFOUND")
    SET(PLASMA_name_binary         "PLASMA_name_binary-NOTFOUND"        )
    IF(BUILD_SINGLE)
        LIST(APPEND PLASMA_name_fct_test "PLASMA_sgetrf"                )
    ENDIF()
    IF(BUILD_DOUBLE)
        LIST(APPEND PLASMA_name_fct_test "PLASMA_dgetrf"                )
    ENDIF()
    IF(BUILD_COMPLEX)
        LIST(APPEND PLASMA_name_fct_test "PLASMA_cgetrf"                )
    ENDIF()
    IF(BUILD_COMPLEX16)
        LIST(APPEND PLASMA_name_fct_test "PLASMA_zgetrf"                )
    ENDIF()

ENDMACRO(PLASMA_INFO_FIND)

##
## @end file infoPLASMA.cmake
##
