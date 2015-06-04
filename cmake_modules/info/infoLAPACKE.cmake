###
#
# @file          : infoLAPACKE.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mer. 16 mai 2012 10:18:07 CEST
#
###

MACRO(LAPACKE_INFO_INSTALL)
    # Define web link of lapacke
    # --------------------------
    IF(NOT DEFINED LAPACKE_URL)
        SET(LAPACKE_URL     "http://icl.cs.utk.edu/projectsfiles/plasma/pubs/lapacke-3.3.0.tgz")
        SET(LAPACKE_MD5SUM  "1e2c8764f41f75bf14919b105513fd25"                                 )
    ENDIF()

    # Define tarball of lapacke
    # -------------------------
    IF(NOT DEFINED LAPACKE_TARBALL)
        SET(LAPACKE_TARBALL "lapacke-3.3.0.tgz"               )
        SET(LAPACKE_MD5SUM  "1e2c8764f41f75bf14919b105513fd25")
    ENDIF()

    # Define repository of lapacke
    # ----------------------------
    IF(NOT DEFINED LAPACKE_SVN_REP)
        SET(LAPACKE_REPO_MODE "SVN")
        SET(LAPACKE_SVN_REP   ""   )
        SET(LAPACKE_SVN_ID    ""   )
        SET(LAPACKE_SVN_PWD   ""   )
    ENDIF()

   # Define dependencies
   # -------------------
   SET(LAPACKE_DEPENDENCIES "blas;lapack")

ENDMACRO(LAPACKE_INFO_INSTALL)

MACRO(LAPACKE_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(LAPACKE_name_library        "lapacke"                             )
    SET(LAPACKE_name_pkgconfig      "LAPACKE_name_pkgconfig-NOTFOUND"     )
    SET(LAPACKE_name_include        "lapacke.h"                           )
    SET(LAPACKE_name_include_suffix "LAPACKE_name_include_suffix-NOTFOUND")
    SET(LAPACKE_name_fct_test       "LAPACKE_cheev"                       )
    SET(LAPACKE_name_binary         "LAPACKE_name_binary-NOTFOUND"        )
ENDMACRO(LAPACKE_INFO_FIND)

###
### END infoLAPACKE.cmake
###
