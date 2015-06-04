###
#
#  @file FindCBLAS.cmake
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
# This module finds an installed CBLAS library.
# CBLAS is distributed at : 
#
# This module can be tuned for searching in a specific folder:
#  - CBLAS_DIR      : path of the presumed top directory of CBLAS
#               (the library will be searched in this in the first time).
#
#
# This module sets the following variables:
#  - CBLAS_FOUND        : set to true if a library is found.
#
#  - CBLAS_PATH         : path of the top directory of cblas installation.
#  - CBLAS_VERSION      : version of cblas.
#
#  - CBLAS_LDFLAGS      : string of all required linker flags.
#  - CBLAS_LIBRARY      : list of all required library  (absolute path).
#  - CBLAS_LIBRARIES    : list of required linker flags (without -l and -L).
#  - CBLAS_LIBRARY_PATH : path of the library directory of cblas installation.
#
#  - CBLAS_INLCUDE_PATH : path of the include directory of cblas installation.
#
###

# Early exit if already searched
IF(CBLAS_FOUND)
    MESSAGE(STATUS "Looking for CBLAS - already found")
    RETURN()
ENDIF(CBLAS_FOUND)

# Load required modules
INCLUDE(populatePACKAGE)
INCLUDE(findPACKAGE)
INCLUDE(infoCBLAS)

# Begin section - Looking for CBLAS
MESSAGE(STATUS "Looking for CBLAS")

# Define parameters for FIND_MY_PACKAGE
CBLAS_INFO_FIND()
SET(BACKUP_CBLAS_name_library "${CBLAS_name_library}")

# Search for the library
MESSAGE(STATUS "Looking for CBLAS - check only with given ${BLAS_VENDOR} BLAS library")
SET(CBLAS_name_library "${BLAS_LIBRARIES}")
FIND_MY_PACKAGE("CBLAS")

# If the BLAS detected implementation do not contains CBLAS,
# we try to find the Netlib version
IF(CBLAS_FOUND)
    SET(CBLAS_VENDOR "${BLAS_VENDOR}")

ELSE(CBLAS_FOUND)
    # Looking for dependencies
    FIND_PACKAGE(BLAS QUIET)
    IF(BLAS_FOUND)
        POPULATE_COMPILE_SYSTEM("BLAS")
        SET(CBLAS_FIND_DEPS "BLAS")
    ENDIF(BLAS_FOUND)

    # Search for the library
    MESSAGE(STATUS "Looking for CBLAS - check Netlib implementation")
    SET(CBLAS_name_library "${BACKUP_CBLAS_name_library}")
    FIND_MY_PACKAGE("CBLAS")

ENDIF(CBLAS_FOUND)

# Begin section - Looking for CBLAS
IF(CBLAS_FOUND)
    MESSAGE(STATUS "Looking for CBLAS - found")
ELSE(CBLAS_FOUND)
    MESSAGE(STATUS "Looking for CBLAS - not found")
ENDIF(CBLAS_FOUND)

##
## @end file FindCBLAS.cmake
##
