###
#
# @file          : infoREFBLAS.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mar. 12 juin 2012 16:31:23 CEST
#
###

MACRO(REFBLAS_INFO_INSTALL)
    # Define web link of blas
    # -----------------------
    IF(NOT DEFINED REFBLAS_URL)
        SET(REFBLAS_URL     "http://netlib.org/blas/blas.tgz")
    ENDIF()

    # Define tarball of blas
    # ----------------------
    IF(NOT DEFINED REFBLAS_TARBALL)
        SET(REFBLAS_TARBALL "blas.tgz")
    ENDIF()

    # Define md5sum of blas
    # ---------------------
    IF(DEFINED REFBLAS_URL OR DEFINED REFBLAS_TARBALL)
        SET(REFBLAS_MD5SUM  "5e99e975f7a1e3ea6abcad7c6e7e42e6")
    ENDIF()

    # Define repository of blas
    # -------------------------
    IF(NOT DEFINED REFBLAS_SVN_REP)
        SET(REFBLAS_REPO_MODE "SVN")
        SET(REFBLAS_SVN_REP   ""   )
        SET(REFBLAS_SVN_ID    ""   )
        SET(REFBLAS_SVN_PWD   ""   )
    ENDIF()

   # Define dependencies
   # -------------------
   SET(REFBLAS_DEPENDENCIES "")

ENDMACRO(REFBLAS_INFO_INSTALL)

MACRO(REFBLAS_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(REFBLAS_name_library        "BLAS_name_library-NOTFOUND"       )
    SET(REFBLAS_name_pkgconfig      "BLAS_name_pkgconfig-NOTFOUND"     )
    SET(REFBLAS_name_include        "BLAS_name_include-NOTFOUND"       )
    SET(REFBLAS_name_include_suffix "BLAS_name_include_suffix-NOTFOUND")
    SET(REFBLAS_name_fct_test       "dgemm"                            )
    SET(REFBLAS_name_binary         "BLAS_name_binary-NOTFOUND"        )

ENDMACRO(REFBLAS_INFO_FIND)

###
### END infoREFBLAS.cmake
###
