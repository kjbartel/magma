###
#
# @file          : infoOPENMPI.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mer. 16 mai 2012 10:18:53 CEST
#
###

MACRO(MPI_INFO_INSTALL)
    # Define web link of openmpi
    # --------------------------
    IF(NOT DEFINED MPI_URL)
        SET(MPI_URL     "http://www.open-mpi.org/software/ompi/v1.5/downloads/openmpi-1.5.4.tar.gz")
        SET(MPI_MD5SUM  "cb11927986419374c54aa4c878209913"                                         )
    ENDIF()

    # Define tarball of openmpi
    # -------------------------
    IF(NOT DEFINED MPI_TARBALL)
        SET(MPI_TARBALL "openmpi-1.5.4.tar.gz"            )
        SET(MPI_MD5SUM  "cb11927986419374c54aa4c878209913")
    ENDIF()

    # Define repository of openmpi
    # ----------------------------
    IF(NOT DEFINED MPI_SVN_REP)
        SET(MPI_REPO_MODE "SVN"                                   )
        SET(MPI_SVN_REP   "http://svn.open-mpi.org/svn/ompi/trunk")
        SET(MPI_SVN_ID    ""                                      )
        SET(MPI_SVN_PWD   ""                                      )
    ENDIF()

   # Define dependencies
   # -------------------
   SET(MPI_DEPENDENCIES "")

ENDMACRO(MPI_INFO_INSTALL)

MACRO(MPI_INFO_FIND)
    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(MPI_name_library        "MPI_name_library-NOTFOUND"       )
    SET(MPI_name_pkgconfig      "MPI_name_pkgconfig-NOTFOUND"     )
    SET(MPI_name_include        "MPI_name_include-NOTFOUND"       )
    SET(MPI_name_include_suffix "MPI_name_include_suffix-NOTFOUND")
    SET(MPI_name_fct_test       "MPI_INIT"                        )
    SET(MPI_name_binary         "MPI_name_binary-NOTFOUND"        )
ENDMACRO(MPI_INFO_FIND)

###
### END infoOPENMPI.cmake
###
