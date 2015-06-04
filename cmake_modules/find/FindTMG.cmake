###
#
#  @file FindTMG.cmake
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
#
# This module finds an installed TMG library.
# TMG is distributed at : 
#
# This module can be tuned for searching in a specific folder:
#  - TMG_DIR      : path of the presumed top directory of TMG
#               (the library will be searched in this in the first time).
#
#
# This module sets the following variables:
#  - TMG_FOUND        : set to true if a library is found.
#
#  - TMG_PATH         : path of the top directory of tmg installation.
#  - TMG_VERSION      : version of tmg.
#
#  - TMG_LDFLAGS      : string of all required linker flags.
#  - TMG_LIBRARY      : list of all required library  (absolute path).
#  - TMG_LIBRARIES    : list of required linker flags (without -l and -L).
#  - TMG_LIBRARY_PATH : path of the library directory of tmg installation.
#
###

# Early exit if already searched
IF(TMG_FOUND)
    MESSAGE(STATUS "Looking for TMG - already found")
    RETURN()
ENDIF(TMG_FOUND)

# Load required modules
INCLUDE(populatePACKAGE)
INCLUDE(findPACKAGE)
INCLUDE(infoTMG)

# Begin section - Looking for TMG
MESSAGE(STATUS "Looking for TMG")

# Define parameters for FIND_MY_PACKAGE
TMG_INFO_FIND()
SET(BACKUP_TMG_name_library "${TMG_name_library}")

# Search for the library
MESSAGE(STATUS "Looking for TMG - check only with given ${BLAS_VENDOR} BLAS library")
SET(TMG_name_library "${BLAS_LIBRARIES}")
FIND_MY_PACKAGE("TMG")

# If the BLAS detected implementation do not contains TMG,
# we try to find the Netlib version
IF(TMG_FOUND)
    SET(TMG_VENDOR "${BLAS_VENDOR}")

ELSE(TMG_FOUND)
    # Looking for dependencies - LAPACK
    FIND_PACKAGE(LAPACK QUIET)
    IF(LAPACK_FOUND)
        POPULATE_COMPILE_SYSTEM("LAPACK")
        SET(TMG_FIND_DEPS "LAPACK")
    ENDIF(LAPACK_FOUND)

    # Looking for dependencies - BLAS
    FIND_PACKAGE(BLAS QUIET)
    IF(BLAS_FOUND)
        POPULATE_COMPILE_SYSTEM("BLAS")
        SET(TMG_FIND_DEPS "BLAS")
    ENDIF(BLAS_FOUND)

    # Search for the library
    MESSAGE(STATUS "Looking for TMG - check Netlib implementation")
    SET(TMG_name_library "${BACKUP_TMG_name_library}")
    FIND_MY_PACKAGE("TMG")

ENDIF(TMG_FOUND)

# Begin section - Looking for TMG
IF(TMG_FOUND)
    MESSAGE(STATUS "Looking for TMG - found")
ELSE(TMG_FOUND)
    MESSAGE(STATUS "Looking for TMG - not found")
ENDIF(TMG_FOUND)

##
## @end file FindTMG.cmake
##
