###
#
# @file      : FindLAPACKE.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : mer. 04 avril 2012 10:56:12 CEST
#
###
#
# This module finds an installed LAPACKE library.
# LAPACKE is distributed at : 
#
# This module can be tuned for searching in a specific folder:
#  - LAPACKE_DIR      : path of the presumed top directory of LAPACKE
#               (the library will be searched in this in the first time).
#
#
# This module sets the following variables:
#  - LAPACKE_FOUND        : set to true if a library is found.
#
#  - LAPACKE_PATH         : path of the top directory of lapacke installation.
#  - LAPACKE_VERSION      : version of lapacke.
#
#  - LAPACKE_LDFLAGS      : string of all required linker flags.
#  - LAPACKE_LIBRARY      : list of all required library  (absolute path).
#  - LAPACKE_LIBRARIES    : list of required linker flags (without -l and -L).
#  - LAPACKE_LIBRARY_PATH : path of the library directory of lapacke installation.
#
#  - LAPACKE_INLCUDE_PATH : path of the include directory of lapacke installation.
#
###

INCLUDE(findPACKAGE)
INCLUDE(infoLAPACKE)

# Begin section - Looking for LAPACKE
MESSAGE(STATUS "Looking for LAPACKE")

# Define extra directory to look for
IF(MORSE_SEPARATE_PROJECTS)
    SET(LAPACKE_DIR ${CMAKE_INSTALL_PREFIX}/lapacke)
ELSE(MORSE_SEPARATE_PROJECTS)
    SET(LAPACKE_DIR ${CMAKE_INSTALL_PREFIX})
ENDIF(MORSE_SEPARATE_PROJECTS)

# Define parameters for FIND_MY_PACKAGE
LAPACKE_INFO_FIND()

# Search for the library
FIND_MY_PACKAGE("LAPACKE"
                FALSE TRUE TRUE FALSE
                TRUE FALSE)

# Begin section - Looking for LAPACKE
IF(LAPACKE_FOUND)
    MESSAGE(STATUS "Looking for LAPACKE - found")
ELSE(LAPACKE_FOUND)
    MESSAGE(STATUS "Looking for LAPACKE - not found")
ENDIF(LAPACKE_FOUND)

###
### END FindLAPACKE.cmake
###
