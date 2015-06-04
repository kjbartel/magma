###
#
# @file          : infoFXT.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mer. 16 mai 2012 10:16:32 CEST
#
###

MACRO(FXT_INFO_INSTALL)
    # Define web link of fxt
    # ----------------------
    IF(NOT DEFINED FXT_URL)
        SET(FXT_URL     "http://download.savannah.gnu.org/releases/fkt/fxt-0.2.5.tar.gz")
        SET(FXT_MD5SUM  "0ea5a93ff69dbb94bb0b99714df4a004")
    ENDIF()

    # Define tarball of fxt
    # ---------------------
    IF(NOT DEFINED FXT_TARBALL)
        SET(FXT_TARBALL "fxt-0.2.5.tar.gz")
        SET(FXT_MD5SUM  "0ea5a93ff69dbb94bb0b99714df4a004")
    ENDIF()

    # Define repository of fxt
    # ------------------------
    IF(NOT DEFINED FXT_SVN_REP)
        SET(FXT_REPO_MODE "CVS"                                 )
        SET(FXT_CVS_REP   "cvs.savannah.nongnu.org:/sources/fkt")
        SET(FXT_CVS_MOD   "fkt"                                 )
    ENDIF()

   # Define dependencies
   # -------------------
   SET(FXT_DEPENDENCIES "")

ENDMACRO(FXT_INFO_INSTALL)

MACRO(FXT_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(FXT_name_library        "fxt"                             )
    SET(FXT_name_pkgconfig      "fxt"                             )
    SET(FXT_name_include        "fxt.h"                           )
    SET(FXT_name_include_suffix "FXT_name_include_suffix-NOTFOUND")
    SET(FXT_name_fct_test       "fxt_infos"                       )
    SET(FXT_name_binary         "fxt_print"                       )
ENDMACRO(FXT_INFO_FIND)

###
### END infoFXT.cmake
###
