###
#
# @file      : FindSTARPU.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 20-01-2012
# @last modified : mer. 04 avril 2012 10:56:41 CEST
#
###
#
# This module finds an installed STARPU library.
# STARPU is distributed at : 
#
# This module can be tuned for searching in a specific folder:
#  - STARPU_DIR      : path of the presumed top directory of STARPU
#               (the library will be searched in this in the first time).
#
#
# This module sets the following variables:
#  - STARPU_FOUND        : set to true if a library is found.
#
#  - STARPU_PATH         : path of the top directory of starpu installation.
#  - STARPU_VERSION      : version of starpu.
#
#  - STARPU_LDFLAGS      : string of all required linker flags (with -l and -L).
#  - STARPU_LIBRARY      : list of all required library  (absolute path).
#  - STARPU_LIBRARIES    : list of required linker flags (without -l and -L).
#  - STARPU_LIBRARY_PATH : path of the library directory of starpu installation.
#
#  - STARPU_INLCUDE_PATH : path of the include directory of starpu installation.
#
#  - STARPU_BINARY_PATH  : path of the binary directory of starpu installation.
#
###

INCLUDE(findPACKAGE)
INCLUDE(infoSTARPU)

# Begin section - Looking for STARPU
MESSAGE(STATUS "Looking for STARPU")

# Define extra directory to look for
IF(MORSE_SEPARATE_PROJECTS)
    SET(STARPU_DIR ${CMAKE_INSTALL_PREFIX}/starpu)
ELSE(MORSE_SEPARATE_PROJECTS)
    SET(STARPU_DIR ${CMAKE_INSTALL_PREFIX})
ENDIF(MORSE_SEPARATE_PROJECTS)

# Define parameters for FIND_MY_PACKAGE
STARPU_INFO_FIND()

# Search for the library
FIND_MY_PACKAGE("STARPU"
                TRUE TRUE TRUE TRUE
                TRUE FALSE)

# Begin section - Looking for STARPU
IF(STARPU_FOUND)
    MESSAGE(STATUS "Looking for STARPU - found")
ELSE(STARPU_FOUND)
    MESSAGE(STATUS "Looking for STARPU - not found")
ENDIF(STARPU_FOUND)

# Define if we have to use deprecated API
#IF("0.9.9" LESS ${STARPU_VERSION})
#    SET(STARPU_NEED_DEPRECATED_API ON)
#ELSE()
#    SET(STARPU_NEED_DEPRECATED_API OFF)
#ENDIF()

###
### END FindSTARPU.cmake
###
