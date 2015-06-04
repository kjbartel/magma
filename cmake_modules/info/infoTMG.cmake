###
#
#  @file infoTMG.cmake
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

MACRO(TMG_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    UNSET(TMG_name_fct_test)
    SET(TMG_type_library        "C;Fortran"                       )
    SET(TMG_name_library        "tmg"                             )
    SET(TMG_name_include        "TMG_name_include-NOTFOUND"       )
    SET(TMG_name_pkgconfig      "TMG_name_pkgconfig-NOTFOUND"     )
    SET(TMG_name_include_suffix "TMG_name_include_suffix-NOTFOUND")
    SET(TMG_name_binary         "TMG_name_binary-NOTFOUND"        )
    IF(BUILD_SINGLE)
        LIST(APPEND TMG_name_fct_test "slatms"                    )
    ENDIF()
    IF(BUILD_DOUBLE)
        LIST(APPEND TMG_name_fct_test "dlatms"                    )
    ENDIF()
    IF(BUILD_COMPLEX)
        LIST(APPEND TMG_name_fct_test "clatms"                    )
    ENDIF()
    IF(BUILD_COMPLEX16)
        LIST(APPEND TMG_name_fct_test "zlatms"                    )
    ENDIF()

ENDMACRO(TMG_INFO_FIND)

##
## @end file infoTMG.cmake
##
