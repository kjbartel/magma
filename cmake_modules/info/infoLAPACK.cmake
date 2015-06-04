###
#
#  @file infoLAPACK.cmake
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

MACRO(LAPACK_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(LAPACK_type_library        "C;Fortran"                          )
    SET(LAPACK_name_library        "LAPACK_name_library-NOTFOUND"       )
    SET(LAPACK_name_pkgconfig      "LAPACK_name_pkgconfig-NOTFOUND"     )
    SET(LAPACK_name_include        "LAPACK_name_include-NOTFOUND"       )
    SET(LAPACK_name_include_suffix "LAPACK_name_include_suffix-NOTFOUND")
    SET(LAPACK_name_binary         "LAPACK_name_binary-NOTFOUND"        )
    IF(BUILD_SINGLE)
        LIST(APPEND LAPACK_name_fct_test "ssyev"                        )
    ENDIF()
    IF(BUILD_DOUBLE)
        LIST(APPEND LAPACK_name_fct_test "dsyev"                        )
    ENDIF()
    IF(BUILD_COMPLEX)
        LIST(APPEND LAPACK_name_fct_test "cheev"                        )
    ENDIF()
    IF(BUILD_COMPLEX16)
        LIST(APPEND LAPACK_name_fct_test "zheev"                        )
    ENDIF()

ENDMACRO(LAPACK_INFO_FIND)

##
## @end file infoLAPACK.cmake
##
