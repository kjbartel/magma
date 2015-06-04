###
#
# @file      : FindCBLAS.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : mer. 04 avril 2012 10:55:10 CEST
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

INCLUDE(findPACKAGE)
INCLUDE(infoCBLAS)

# Begin section - Looking for CBLAS
MESSAGE(STATUS "Looking for CBLAS")

# Define extra directory to look for
IF(MORSE_SEPARATE_PROJECTS)
    SET(CBLAS_DIR ${CMAKE_INSTALL_PREFIX}/cblas)
ELSE(MORSE_SEPARATE_PROJECTS)
    SET(CBLAS_DIR ${CMAKE_INSTALL_PREFIX})
ENDIF(MORSE_SEPARATE_PROJECTS)

# Define parameters for FIND_MY_PACKAGE
CBLAS_INFO_FIND()

# Search for the library
FIND_MY_PACKAGE("CBLAS"
                FALSE TRUE TRUE FALSE
                TRUE FALSE)

# Begin section - Looking for CBLAS
IF(CBLAS_FOUND)
    MESSAGE(STATUS "Looking for CBLAS - found")
ELSE(CBLAS_FOUND)
    MESSAGE(STATUS "Looking for CBLAS - not found")
ENDIF(CBLAS_FOUND)

###
### END FindCBLAS.cmake
###
