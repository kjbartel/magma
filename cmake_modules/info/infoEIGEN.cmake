###
#
# @file          : infoEIGEN.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mar. 12 juin 2012 16:30:50 CEST
#
###

MACRO(EIGEN_INFO_INSTALL)
    # Define web link of eigen
    # ------------------------
    IF(NOT DEFINED EIGEN_URL)
        SET(EIGEN_URL     "http://bitbucket.org/eigen/eigen/get/3.1.0-alpha2.tar.gz")
    ENDIF()

    # Define tarball of eigen
    # -----------------------
    IF(NOT DEFINED EIGEN_TARBALL)
        SET(EIGEN_TARBALL "3.1.0-alpha2.tar.gz")
    ENDIF()

    # Define md5sum of eigen
    # ----------------------
    IF(DEFINED EIGEN_URL OR DEFINED EIGEN_TARBALL)
        SET(EIGEN_MD5SUM  "48ee0624e75a054587e15367ec221907")
    ENDIF()

    # Define repository of eigen
    # --------------------------
    IF(NOT DEFINED EIGEN_SVN_REP)
        SET(EIGEN_REPO_MODE "")
        SET(EIGEN_SVN_REP   "")
        SET(EIGEN_SVN_ID    "")
        SET(EIGEN_SVN_PWD   "")
    ENDIF()

   # Define dependencies
   # -------------------
   SET(EIGEN_DEPENDENCIES "")

ENDMACRO(EIGEN_INFO_INSTALL)

MACRO(EIGEN_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(EIGEN_name_library        "BLAS_name_library-NOTFOUND"       )
    SET(EIGEN_name_pkgconfig      "BLAS_name_pkgconfig-NOTFOUND"     )
    SET(EIGEN_name_include        "BLAS_name_include-NOTFOUND"       )
    SET(EIGEN_name_include_suffix "BLAS_name_include_suffix-NOTFOUND")
    SET(EIGEN_name_fct_test       "dgemm"                            )
    SET(EIGEN_name_binary         "BLAS_name_binary-NOTFOUND"        )

ENDMACRO(EIGEN_INFO_FIND)

###
### END infoEIGEN.cmake
###
