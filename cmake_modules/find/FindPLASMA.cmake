###
#
#  @file FindPLASMA.cmake
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

# Early exit if already searched
IF(PLASMA_FOUND)
    MESSAGE(STATUS "Looking for PLASMA - already found")
    RETURN()
ENDIF(PLASMA_FOUND)

# Load required modules
INCLUDE(BuildSystemTools)
INCLUDE(populatePACKAGE)
INCLUDE(findPACKAGE)
INCLUDE(infoPLASMA)

# Begin section - Looking for PLASMA
MESSAGE(STATUS "Looking for PLASMA")

# Define parameters for FIND_MY_PACKAGE
PLASMA_INFO_FIND()

# Search for the library
FIND_MY_PACKAGE("PLASMA")

IF(PLASMA_FOUND)
    # Fix about link : plasma;coreblas => plasma;coreblas;plasma
    LIST(FIND PLASMA_LIBRARIES "coreblas" _index)
    IF(NOT ${_index} EQUAL "-1")
        MATH(EXPR _index "${_index}+1")
        LIST(INSERT PLASMA_LIBRARIES "${_index}" "plasma")
    ENDIF()

    SET(_idx_cb "-1")
    SET(_idx_cb "-1")
    LIST(LENGTH PLASMA_LIBRARY size_list)
    MATH(EXPR size_list "${size_list}-1")
    FOREACH(i RANGE 0 ${size_list})
        LIST(GET PLASMA_LIBRARY ${i} el_list)
        IF("${el_list}" MATCHES ".*/libplasma\\.([a-z]*)") 
            SET(_idx_pl ${i})
            SET(_add_pl ${el_list})
        ELSEIF("${el_list}" MATCHES ".*/libcoreblas\\.([a-z]*)")
            SET(_idx_cb ${i})
       ENDIF() 
    ENDFOREACH()
    IF(NOT ${_idx_pl} EQUAL "-1" AND NOT ${_idx_cb} EQUAL "-1")
        MATH(EXPR _index "${_idx_cb}+1")
        LIST(INSERT PLASMA_LIBRARY "${_index}" "${_add_pl}")
    ENDIF()

    # If PLASMA is found, we have to set BLAS, CBLAS, LAPACK and LAPACKE variables correctly
    SET(_clean_libs ${PLASMA_LIBRARIES})
    LIST(REMOVE_ITEM _clean_libs    plasma coreblas quark hwloc tmg pthread)
    SET(_clean_library ${PLASMA_LIBRARY})
    FOREACH(_library ${PLASMA_LIBRARY})
        IF("^${_library}$" MATCHES "^(.*)libplasma(.*)$"   OR
           "^${_library}$" MATCHES "^(.*)libcoreblas(.*)$" OR
           "^${_library}$" MATCHES "^(.*)libquark(.*)$"    OR
           "^${_library}$" MATCHES "^(.*)libhwloc(.*)$"    OR
           "^${_library}$" MATCHES "^(.*)libtmg(.*)$"      OR
           "^${_library}$" MATCHES "^(.*)libpthread(.*)$")
            LIST(REMOVE_ITEM _clean_library ${_library})
        ENDIF()
    ENDFOREACH()
    FOREACH(_package BLAS CBLAS LAPACK LAPACKE)
        # Message to prevent users
        MESSAGE(STATUS "Looking for PLASMA - FindPLASMA will set ${_package} flags")
        IF("^${MAGMA_USE_${_package}}$" STREQUAL "^ON$" AND NOT "^${MAGMA_USE_${_package}}$" STREQUAL "^OFF$" AND NOT"${MAGMA_USE_${_package}}" STREQUAL "^$")
            MESSAGE(STATUS "Looking for PLASMA - MAGMA_USE_${_package} was defined as '${MAGMA_USE_${_package}}' but not used")
        ENDIF()

        # Remove the fact that we have to build ${_package}
        SET(${_package}_USED_MODE "FIND")
        SET(HAVE_${_package}      ON    )

        # set ${_package} variables
        SET(${_package}_FOUND        TRUE                      )
        SET(${_package}_PATH         "${PLASMA_PATH}"          ) 
        SET(${_package}_VERSION      "${PLASMA_VERSION}"       )
        SET(${_package}_LDFLAGS      ""                        )
        FOREACH(_path ${PLASMA_LIBRARY_PATH})
            ADD_FLAGS(${_package}_LDFLAGS "-L${_path}"         ) 
        ENDFOREACH()
        FOREACH(_lib ${_clean_libs})
            ADD_FLAGS(${_package}_LDFLAGS  "-l${_lib}"         )
        ENDFOREACH()
        SET(${_package}_LIBRARY      "${_clean_library}"       )
        SET(${_package}_LIBRARIES    "${_clean_libs}"          )
        SET(${_package}_LIBRARY_PATH "${PLASMA_LIBRARY_PATH}"  )
    ENDFOREACH()
    SET(BLAS_VENDOR   "provided by PLASMA")
    SET(LAPACK_VENDOR "provided by PLASMA")
    SET(CBLAS_INCLUDE_PATH   "${PLASMA_INCLUDE_PATH}"  )
    SET(LAPACKE_INCLUDE_PATH "${PLASMA_INCLUDE_PATH}"  )

    # If PLASMA is found, we have to set QUARK variables correctly
    MESSAGE(STATUS "Looking for PLASMA - FindPLASMA will set QUARK flags")
    IF("^${MAGMA_USE_QUARK}$" STREQUAL "^ON$" AND NOT "^${MAGMA_USE_QUARK}$" STREQUAL "^OFF$" AND NOT"${MAGMA_USE_${_package}}" STREQUAL "^$")
        MESSAGE(STATUS "Looking for PLASMA - MAGMA_USE_QUARK was defined as '${MAGMA_USE_QUARK}' but not used")
    ENDIF()
    SET(QUARK_USED_MODE    "FIND"                           )
    SET(HAVE_QUARK         ON                               )
    SET(QUARK_FOUND        TRUE                             )
    SET(QUARK_PATH         "${PLASMA_PATH}"                 ) 
    SET(QUARK_VERSION      "${PLASMA_VERSION}"              )
    SET(QUARK_LIBRARY_PATH "${PLASMA_LIBRARY_PATH}"         )
    SET(QUARK_INCLUDE_PATH "${PLASMA_INCLUDE_PATH}"         )
    SET(QUARK_LDFLAGS      "-L${QUARK_LIBRARY_PATH} -lquark")
    SET(QUARK_LIBRARY      ""                               )
    SET(QUARK_LIBRARIES    "quark"                          )
    FOREACH(_library ${PLASMA_LIBRARY})
        IF("^${_library}$" MATCHES "^(.*)libquark(.*)$")
            LIST(APPEND    QUARK_LIBRARY ${_library}        )
        ENDIF()
    ENDFOREACH()

    # Populate the dependencies
    POPULATE_COMPILE_SYSTEM("BLAS")
    POPULATE_COMPILE_SYSTEM("CBLAS")
    POPULATE_COMPILE_SYSTEM("LAPACK")
    POPULATE_COMPILE_SYSTEM("LAPACKE")
    POPULATE_COMPILE_SYSTEM("QUARK")
    
ENDIF()

# Begin section - Looking for PLASMA
IF(PLASMA_FOUND)
    MESSAGE(STATUS "Looking for PLASMA - found")
ELSE(PLASMA_FOUND)
    MESSAGE(STATUS "Looking for PLASMA - not found")
ENDIF(PLASMA_FOUND)

##
## @end file FindPLASMA.cmake
##
