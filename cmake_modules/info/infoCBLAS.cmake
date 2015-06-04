###
#
#  @file infoCBLAS.cmake
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

MACRO(CBLAS_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    UNSET(CBLAS_name_fct_test)
    SET(CBLAS_type_library        "C"                                 )
    SET(CBLAS_name_library        "cblas"                             )
    SET(CBLAS_name_include        "CBLAS_name_include-NOTFOUND"       )
    SET(CBLAS_name_pkgconfig      "CBLAS_name_pkgconfig-NOTFOUND"     )
    SET(CBLAS_name_include_suffix "CBLAS_name_include_suffix-NOTFOUND")
    SET(CBLAS_name_binary         "CBLAS_name_binary-NOTFOUND"        )
    IF(BUILD_SINGLE)
        LIST(APPEND CBLAS_name_fct_test "cblas_sgemm"                 )
    ENDIF()
    IF(BUILD_DOUBLE)
        LIST(APPEND CBLAS_name_fct_test "cblas_dgemm"                 )
    ENDIF()
    IF(BUILD_COMPLEX)
        LIST(APPEND CBLAS_name_fct_test "cblas_cgemm"                 )
    ENDIF()
    IF(BUILD_COMPLEX16)
        LIST(APPEND CBLAS_name_fct_test "cblas_zgemm"                 )
    ENDIF()

ENDMACRO(CBLAS_INFO_FIND)

##
## @end file infoCBLAS.cmake
##
