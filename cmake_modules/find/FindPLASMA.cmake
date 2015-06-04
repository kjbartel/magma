###
#
# @file      : FindPLASMA.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : mer. 04 avril 2012 12:06:39 CEST
#
###
#
# This module finds an installed PLASMA library.
# PLASMA is distributed at : 
#
# This module can be tuned for searching in a specific folder:
#  - PLASMA_DIR      : path of the presumed top directory of PLASMA
#               (the library will be searched in this in the first time).
#
#
# This module sets the following variables:
#  - PLASMA_FOUND        : set to true if a library is found.
#
#  - PLASMA_PATH         : path of the top directory of plasma installation.
#  - PLASMA_VERSION      : version of plasma.
#
#  - PLASMA_LDFLAGS      : string of all required linker flags.
#  - PLASMA_LIBRARY      : list of all required library  (absolute path).
#  - PLASMA_LIBRARIES    : list of required linker flags (without -l and -L).
#  - PLASMA_LIBRARY_PATH : path of the library directory of plasma installation.
#
#  - PLASMA_INLCUDE_PATH : path of the include directory of plasma installation.
#
###

INCLUDE(findPACKAGE)
INCLUDE(infoPLASMA)

# Begin section - Looking for PLASMA
MESSAGE(STATUS "Looking for PLASMA")

# Define extra directory to look for
IF(MORSE_SEPARATE_PROJECTS)
    SET(PLASMA_DIR ${CMAKE_INSTALL_PREFIX}/plasma)
ELSE(MORSE_SEPARATE_PROJECTS)
    SET(PLASMA_DIR ${CMAKE_INSTALL_PREFIX})
ENDIF(MORSE_SEPARATE_PROJECTS)

# Define parameters for FIND_MY_PACKAGE
PLASMA_INFO_FIND()

# Search for the library
FIND_MY_PACKAGE("PLASMA"
                FALSE TRUE TRUE TRUE
                TRUE TRUE)

IF(PLASMA_FOUND)
    # Fix about link : plasma;coreblas => plasma;coreblas;plasma
    LIST(FIND PLASMA_LIBRARIES "coreblas" _index)
    IF(NOT ${_index} EQUAL "-1")
        MATH(EXPR _index "${_index}+1")
        LIST(INSERT PLASMA_LIBRARIES "${_index}" "plasma")
    ENDIF()
    
    # If PLASMA is found, we have to set BLAS, CBLAS, LAPACK and LAPACKE variables correctly
    SET(_clean_libs ${PLASMA_LIBRARIES})
    LIST(REMOVE_ITEM _clean_libs plasma coreblas quark hwloc pthread m)
    FOREACH(_package BLAS CBLAS LAPACK LAPACKE)
        # Remove the fact that we have to build ${_package}
        UNSET(${_package}_BUILD_MODE)
        #STRING(TOLOWER "${_package}" _package_target)
        #SET_TARGET_PROPERTIES(${_package_target}_install PROPERTIES EXCLUDE_FROM_ALL "TRUE")

        # set ${_package} variables
        SET(${_package}_FOUND        TRUE                      )
        SET(${_package}_PATH         "${PLASMA_PATH}"          ) 
        SET(${_package}_VERSION      "${PLASMA_VERSION}"       )
        SET(${_package}_LDFLAGS      "-L${PLASMA_LIBRARY_PATH}")
        FOREACH(_libs ${_clean_libs})
            SET(${_package}_LDFLAGS  " -l${_libs}"             )
        ENDFOREACH()
        SET(${_package}_LIBRARY      ""                        )
        SET(${_package}_LIBRARIES    "${_clean_libs}"          )
        SET(${_package}_LIBRARY_PATH "${PLASMA_LIBRARY_PATH}"  )
        SET(${_package}_INLCUDE_PATH "${PLASMA_INCLUDE_PATH}"  )
    ENDFOREACH()

ENDIF()

# Begin section - Looking for PLASMA
IF(PLASMA_FOUND)
    MESSAGE(STATUS "Looking for PLASMA - found")
ELSE(PLASMA_FOUND)
    MESSAGE(STATUS "Looking for PLASMA - not found")
ENDIF(PLASMA_FOUND)

###
### END FindPLASMA.cmake
###
