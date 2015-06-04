###
#
# @file          : infoHWLOC.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mer. 16 mai 2012 10:17:01 CEST
#
###

MACRO(HWLOC_INFO_INSTALL)
    # Define web link of hwloc
    # ------------------------
    IF(NOT DEFINED HWLOC_URL)
        SET(HWLOC_URL     "http://www.open-mpi.org/software/hwloc/v1.4/downloads/hwloc-1.4.1.tar.gz")
        SET(HWLOC_MD5SUM  "d8b28f45d6a841087a6094a2183f4673")
    ENDIF()

    # Define tarball of hwloc
    # -----------------------
    IF(NOT DEFINED HWLOC_TARBALL)
        SET(HWLOC_TARBALL "hwloc-1.4.1.tar.gz")
        SET(HWLOC_MD5SUM  "d8b28f45d6a841087a6094a2183f4673")
    ENDIF()

    # Define repository of hwloc
    # --------------------------
    IF(NOT DEFINED HWLOC_SVN_REP)
        SET(HWLOC_REPO_MODE "SVN"                                     )
        SET(HWLOC_SVN_REP   "https://svn.open-mpi.org/svn/hwloc/trunk")
        SET(HWLOC_SVN_ID    ""                                        )
        SET(HWLOC_SVN_PWD   ""                                        )
    ENDIF()

   # Define dependencies
   # -------------------
   SET(HWLOC_DEPENDENCIES "")

ENDMACRO(HWLOC_INFO_INSTALL)

MACRO(HWLOC_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(HWLOC_name_library        "hwloc"                             )
    SET(HWLOC_name_pkgconfig      "hwloc"                             )
    SET(HWLOC_name_include        "hwloc.h"                           )
    SET(HWLOC_name_include_suffix "HWLOC_name_include_suffix-NOTFOUND")
    SET(HWLOC_name_fct_test       "hwloc_set_cpubind"                 )
    SET(HWLOC_name_binary         "hwloc-ls"                          )
ENDMACRO(HWLOC_INFO_FIND)

###
### END infoHWLOC.cmake
###
