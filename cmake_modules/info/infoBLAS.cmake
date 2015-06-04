###
#
#  @file infoBLAS.cmake
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

MACRO(BLAS_INFO_FIND)

    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(BLAS_type_library        "Fortran"                          )
    SET(BLAS_name_library        "BLAS_name_library-NOTFOUND"       )
    SET(BLAS_name_pkgconfig      "BLAS_name_pkgconfig-NOTFOUND"     )
    SET(BLAS_name_include        "BLAS_name_include-NOTFOUND"       )
    SET(BLAS_name_include_suffix "BLAS_name_include_suffix-NOTFOUND")
    SET(BLAS_name_binary         "BLAS_name_binary-NOTFOUND"        )

    UNSET(LAS_name_fct_test)
    IF(BUILD_SINGLE)
        LIST(APPEND BLAS_name_fct_test "sgemm")
    ENDIF()
    IF(BUILD_DOUBLE)
        LIST(APPEND BLAS_name_fct_test "dgemm")
    ENDIF()
    IF(BUILD_COMPLEX)
        LIST(APPEND BLAS_name_fct_test "cgemm")
    ENDIF()
    IF(BUILD_COMPLEX16)
        LIST(APPEND BLAS_name_fct_test "zgemm")
    ENDIF()

ENDMACRO(BLAS_INFO_FIND)

##
## @end file infoBLAS.cmake
##
