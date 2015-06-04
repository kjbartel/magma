###
#
# @file      : FindFXT.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : mer. 04 avril 2012 10:55:29 CEST
#
###
#
# This module finds an installed FXT library.
# FXT is distributed at : 
#
# This module can be tuned for searching in a specific folder:
#  - FXT_DIR      : path of the presumed top directory of FXT
#               (the library will be searched in this in the first time).
#
#
# This module sets the following variables:
#  - FXT_FOUND        : set to true if a library is found.
#
#  - FXT_PATH         : path of the top directory of fxt installation.
#  - FXT_VERSION      : version of fxt.
#
#  - FXT_LDFLAGS      : string of all required linker flags (with -l and -L).
#  - FXT_LIBRARY      : list of all required library  (absolute path).
#  - FXT_LIBRARIES    : list of required linker flags (without -l and -L).
#  - FXT_LIBRARY_PATH : path of the library directory of fxt installation.
#
#  - FXT_INLCUDE_PATH : path of the include directory of fxt installation.
#
###

INCLUDE(findPACKAGE)
INCLUDE(infoFXT)

# Begin section - Looking for FXT
MESSAGE(STATUS "Looking for FXT")
    
# Define extra directory to look for
IF(MORSE_SEPARATE_PROJECTS)
    SET(FXT_DIR ${CMAKE_INSTALL_PREFIX}/fxt)
ELSE(MORSE_SEPARATE_PROJECTS)
    SET(FXT_DIR ${CMAKE_INSTALL_PREFIX})
ENDIF(MORSE_SEPARATE_PROJECTS)

# Define parameters for FIND_MY_PACKAGE
FXT_INFO_FIND()

# Search for the library
FIND_MY_PACKAGE("FXT"
                TRUE TRUE TRUE TRUE
                TRUE FALSE)

# Begin section - Looking for FXT
IF(FXT_FOUND)
    MESSAGE(STATUS "Looking for FXT - found")
ELSE(FXT_FOUND)
    MESSAGE(STATUS "Looking for FXT - not found")
ENDIF(FXT_FOUND)

###
### END FindFXT.cmake
###
