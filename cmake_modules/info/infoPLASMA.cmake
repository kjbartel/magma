###
#
# @file          : infoPLASMA.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mer. 16 mai 2012 10:20:16 CEST
#
###

MACRO(PLASMA_INFO_INSTALL)
    # Define web link of plasma
    # -------------------------
    IF(NOT DEFINED PLASMA_URL)
        SET(PLASMA_URL     "http://icl.cs.utk.edu/projectsfiles/plasma/pubs/plasma_2.4.5.tar.gz")
        SET(PLASMA_MD5SUM  "077ce8df88d48e3c2be115b459bb0231"                                   )
    ENDIF()

    # Define tarball of plasma
    # ------------------------
    IF(NOT DEFINED PLASMA_TARBALL)
        SET(PLASMA_TARBALL "plasma_2.4.5.tar.gz"             )
        SET(PLASMA_MD5SUM  "077ce8df88d48e3c2be115b459bb0231")
    ENDIF()

    # Define repository of plasma
    # ---------------------------
    IF(NOT DEFINED PLASMA_SVN_REP)
        SET(PLASMA_REPO_MODE "SVN")
        SET(PLASMA_SVN_REP   ""   )
        SET(PLASMA_SVN_ID    ""   )
        SET(PLASMA_SVN_PWD   ""   )
    ENDIF()

   # Define dependencies
   # -------------------
   SET(PLASMA_DEPENDENCIES "hwloc;blas;lapack;cblas;lapacke")

ENDMACRO(PLASMA_INFO_INSTALL)

MACRO(PLASMA_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(PLASMA_name_library        "plasma;coreblas;quark"              )
    SET(PLASMA_name_pkgconfig      "plasma"                             )
    SET(PLASMA_name_include        "plasma.h;core_blas.h;quark.h"       )
    SET(PLASMA_name_include_suffix "PLASMA_name_include_suffix-NOTFOUND")
    SET(PLASMA_name_fct_test       "PLASMA_dgetrf"                      )
    SET(PLASMA_name_binary         "PLASMA_name_binary-NOTFOUND"        )
ENDMACRO(PLASMA_INFO_FIND)

###
### END infoPLASMA.cmake
###
