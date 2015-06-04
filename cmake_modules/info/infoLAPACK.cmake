###
#
# @file          : infoLAPACK.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mer. 16 mai 2012 10:17:28 CEST
#
###

MACRO(LAPACK_INFO_INSTALL)
    # Define web link of lapack
    # -------------------------
    IF(NOT DEFINED LAPACK_URL)
        SET(LAPACK_URL     "http://www.netlib.org/lapack/lapack.tgz")
        SET(LAPACK_MD5SUM  "44c3869c38c8335c2b9c2a8bb276eb55"       )
    ENDIF()

    # Define tarball of lapack
    # ------------------------
    IF(NOT DEFINED LAPACK_TARBALL)
        SET(LAPACK_TARBALL "lapack.tgz"                      )
        SET(LAPACK_MD5SUM  "44c3869c38c8335c2b9c2a8bb276eb55")
    ENDIF()

    # Define repository of lapack
    # ---------------------------
    IF(NOT DEFINED LAPACK_SVN_REP)
        SET(LAPACK_REPO_MODE "SVN"                                               )
        SET(LAPACK_SVN_REP   "https://icl.cs.utk.edu/svn/lapack-dev/lapack/trunk")
        SET(LAPACK_SVN_ID    ""                                                  )
        SET(LAPACK_SVN_PWD   ""                                                  )
    ENDIF()

   # Define dependencies
   # -------------------
   SET(LAPACK_DEPENDENCIES "blas")

ENDMACRO(LAPACK_INFO_INSTALL)

MACRO(LAPACK_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(LAPACK_name_library        "LAPACK_name_library-NOTFOUND"       )
    SET(LAPACK_name_pkgconfig      "LAPACK_name_pkgconfig-NOTFOUND"     )
    SET(LAPACK_name_include        "LAPACK_name_include-NOTFOUND"       )
    SET(LAPACK_name_include_suffix "LAPACK_name_include_suffix-NOTFOUND")
    SET(LAPACK_name_fct_test       "cheev"                              )
    SET(LAPACK_name_binary         "LAPACK_name_binary-NOTFOUND"        )
ENDMACRO(LAPACK_INFO_FIND)

###
### END infoLAPACK.cmake
###
