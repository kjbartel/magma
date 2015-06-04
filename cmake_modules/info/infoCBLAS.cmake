###
#
# @file          : infoCBLAS.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mer. 16 mai 2012 10:14:06 CEST
#
###

MACRO(CBLAS_INFO_INSTALL)
    # Define web link of cblas
    # ------------------------
    IF(NOT DEFINED CBLAS_URL)
        SET(CBLAS_URL     "http://www.netlib.org/blas/blast-forum/cblas.tgz")
        SET(CBLAS_MD5SUM  "1e8830f622d2112239a4a8a83b84209a"                )
    ENDIF()

    # Define tarball of cblas
    # -----------------------
    IF(NOT DEFINED CBLAS_TARBALL)
        SET(CBLAS_TARBALL "cblas.tgz"                       )
        SET(CBLAS_MD5SUM  "1e8830f622d2112239a4a8a83b84209a")
    ENDIF()

    # Define repository of cblas
    # --------------------------
    IF(NOT DEFINED CBLAS_SVN_REP)
        SET(CBLAS_REPO_MODE "SVN")
        SET(CBLAS_SVN_REP   ""   )
        SET(CBLAS_SVN_ID    ""   )
        SET(CBLAS_SVN_PWD   ""   )
    ENDIF()

   # Define dependencies
   # -------------------
   SET(CBLAS_DEPENDENCIES "blas")

ENDMACRO(CBLAS_INFO_INSTALL)

MACRO(CBLAS_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(CBLAS_name_library        "cblas"                             )
    SET(CBLAS_name_pkgconfig      "CBLAS_name_pkgconfig-NOTFOUND"     )
    SET(CBLAS_name_include        "cblas.h"                           )
    SET(CBLAS_name_include_suffix "CBLAS_name_include_suffix-NOTFOUND")
    SET(CBLAS_name_fct_test       "cblas_dgemm"                       )
    SET(CBLAS_name_binary         "CBLAS_name_binary-NOTFOUND"        )
ENDMACRO(CBLAS_INFO_FIND)

###
### END infoCBLAS.cmake
###
