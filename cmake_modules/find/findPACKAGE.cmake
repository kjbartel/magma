###
#
#  @file findPACKAGE.cmake
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
# This module finds an installed ${_NAMEVAR} library.
#
# This module sets the following variables:
#  - <PACKAGE>_FOUND        : set to true if a library is found.
#
#  - <PACKAGE>_PATH         : path of the top directory of <PACKAGE> installation.
#  - <PACKAGE>_VERSION      : version of the library.
#
#  - <PACKAGE>_LDFLAGS      : string of all required linker flags (with -l and -L). 
#  - <PACKAGE>_LIBRARY      : list of all required library  (absolute path).
#  - <PACKAGE>_LIBRARIES    : list of required linker flags (without -l and -L).
#  - <PACKAGE>_LIBRARY_PATH : path of the library directory of <PACKAGE> installation.
#
#  - <PACKAGE>_INLCUDE_PATH : path of the include directory of <PACKAGE> installation.
#
#  - <PACKAGE>_BINARY_PATH  : path of the binary directory of <PACKAGE> installation.
#
###
#
# You can set the variable DEBUG_MORSE_FINDPACKAGE to ON for printing all debug informations
#   eg: cmake <path> -DDEBUG_MORSE_FINDPACKAGE=ON <options>
#
###

INCLUDE(BuildSystemTools)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckFortranFunctionExists)
INCLUDE(FindPackageHandleStandardArgs)

# Find package
# ------------
MACRO(FIND_MY_PACKAGE _NAMEVAR)

    # Set ${_NAMEVAR}_FOUND at FALSE by default
    # -----------------------------------------
    SET(${_NAMEVAR}_FOUND FALSE)
    SET(${_NAMEVAR}_ERROR_OCCURRED FALSE)

    # Optionally use pkg-config to detect include/library paths (if pkg-config is available)
    # --------------------------------------------------------------------------------------
    FIND_PACKAGE(PkgConfig QUIET)
    IF(PKG_CONFIG_EXECUTABLE AND ${_NAMEVAR}_name_pkgconfig)
        STRING(REPLACE ":" ";" PATH_PKGCONFIGPATH "$ENV{PKG_CONFIG_PATH}")
        FIND_FILE(${_NAMEVAR}_PKG_FILE_FOUND
                  NAME  ${${_NAMEVAR}_name_pkgconfig}.pc
                  PATHS ${PATH_PKGCONFIGPATH})
        MARK_AS_ADVANCED(${_NAMEVAR}_PKG_FILE_FOUND)
        IF(${_NAMEVAR}_PKG_FILE_FOUND)
            pkg_search_module(PC_${_NAMEVAR} QUIET ${${_NAMEVAR}_name_pkgconfig})
            DEBUG_PKG_CONFIG_OUTPUT(${_NAMEVAR})

        ELSE(${_NAMEVAR}_PKG_FILE_FOUND)
            MESSAGE(STATUS "Looking for ${_NAMEVAR} - pkgconfig not used")

        ENDIF(${_NAMEVAR}_PKG_FILE_FOUND)

    ELSE(PKG_CONFIG_EXECUTABLE AND ${_NAMEVAR}_name_pkgconfig)
        MESSAGE(STATUS "Looking for ${_NAMEVAR} - pkgconfig not used")

    ENDIF(PKG_CONFIG_EXECUTABLE AND ${_NAMEVAR}_name_pkgconfig)

    # Looking for version
    # -------------------
    IF(PC_${_NAMEVAR}_FOUND AND ${_NAMEVAR}_name_pkgconfig)
        SET(${_NAMEVAR}_VERSION ${PC_${_NAMEVAR}_VERSION})
    ENDIF(PC_${_NAMEVAR}_FOUND AND ${_NAMEVAR}_name_pkgconfig)

    # Looking for top directory
    # -------------------------
    IF(PC_${_NAMEVAR}_FOUND AND ${_NAMEVAR}_name_pkgconfig)
        SET(${_NAMEVAR}_PATH ${PC_${_NAMEVAR}_PREFIX})
    ENDIF(PC_${_NAMEVAR}_FOUND AND ${_NAMEVAR}_name_pkgconfig)

    # Clean old values
    # ----------------
    UNSET(${_NAMEVAR}_LIBRARY)
    UNSET(${_NAMEVAR}_LDFLAGS)
    UNSET(${_NAMEVAR}_LIBRARIES)
    UNSET(${_NAMEVAR}_LIBRARY_PATH)
    UNSET(${_NAMEVAR}_INCLUDE_PATH)

    # Add some path to help
    # ---------------------
    UNSET(_lib_env)
    IF(WIN32)
        STRING(REPLACE ":" ";" _lib_env "$ENV{LIB}")
    ELSE(WIN32)
        IF(APPLE)
            STRING(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
        ELSE()
            STRING(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
        ENDIF()
        LIST(APPEND _lib_env "/usr/local/lib64")
        LIST(APPEND _lib_env "/usr/lib64")
        LIST(APPEND _lib_env "/usr/local/lib")
        LIST(APPEND _lib_env "/usr/lib")
    ENDIF()

    # Define path we have to looking for first
    # ----------------------------------------
    UNSET(CMAKE_PREFIX_PATH)
    LIST(APPEND CMAKE_PREFIX_PATH ${PC_${_NAMEVAR}_LIBDIR})
    LIST(APPEND CMAKE_PREFIX_PATH ${PC_${_NAMEVAR}_LIBRARY_DIRS})
    LIST(APPEND CMAKE_PREFIX_PATH ${${_NAMEVAR}_DIR})
    LIST(APPEND CMAKE_PREFIX_PATH ${_lib_env})

    # Get and preporcess of collected informations
    # ---------------------------------------------
    SET(_all_libs "")
    LIST(APPEND _all_libs ${${_NAMEVAR}_name_library})
    LIST(APPEND _all_libs ${PC_${_NAMEVAR}_LIBRARIES})
    LIST(REMOVE_DUPLICATES _all_libs)
    LIST(FIND _all_libs "stdc++" is_extra_flags)
    IF("-1" MATCHES ${is_extra_flags})
        SET(is_extra_flags FALSE)
    ELSE()
        SET(is_extra_flags TRUE)
        LIST(REMOVE_ITEM _all_libs "stdc++")
    ENDIF()
    IF(${_NAMEVAR}_name_library)
        FOREACH(_lib_file ${_all_libs})
            FIND_LIBRARY(${_lib_file}_path
                         NAMES ${_lib_file}
                         PATH_SUFFIXES lib lib64 lib32
                        )
            MARK_AS_ADVANCED(${_lib_file}_path)
            IF(${_lib_file}_path)
                GET_FILENAME_COMPONENT(${_lib_file}_filename ${${_lib_file}_path} NAME)
                GET_FILENAME_COMPONENT(${_lib_file}_lib_path ${${_lib_file}_path} PATH)
                LIST(FIND ${_NAMEVAR}_LIBRARIES "${_lib_file}" IS_${_lib_file})
                IF("^${IS_${_lib_file}}$" STREQUAL "^-1$")
                    LIST(APPEND ${_NAMEVAR}_LIBRARIES    "${_lib_file}")
                ENDIF()
                LIST(FIND ${_NAMEVAR}_LIBRARY "${${_lib_file}_path}" IS_${_lib_file}_path)
                IF("^${IS_${_lib_file}_path}$" STREQUAL "^-1$")
                    LIST(APPEND ${_NAMEVAR}_LIBRARY      "${${_lib_file}_path}")
                ENDIF()
                LIST(APPEND ${_NAMEVAR}_LIBRARY_PATH "${${_lib_file}_lib_path}")
                LIST(REMOVE_DUPLICATES ${_NAMEVAR}_LIBRARY_PATH)

            ELSE(${_lib_file}_path)
                SET(${_NAMEVAR}_ERROR_OCCURED TRUE)
                MESSAGE(STATUS "Looking for ${_lib_file} - not found")

            ENDIF(${_lib_file}_path)
        ENDFOREACH()

        SET(${_NAMEVAR}_LDFLAGS "")
        FOREACH(_lib_dir ${${_NAMEVAR}_LIBRARY_PATH})
            SET(${_NAMEVAR}_LDFLAGS "${${_NAMEVAR}_LDFLAGS}-L${_lib_dir} ")
        ENDFOREACH()
        FOREACH(_lib_flags ${${_NAMEVAR}_LIBRARIES})
            SET(${_NAMEVAR}_LDFLAGS "${${_NAMEVAR}_LDFLAGS}-l${_lib_flags} ")
        ENDFOREACH()

        # Fix case of stdc++
        IF(${is_extra_flags})
            LIST(APPEND ${_NAMEVAR}_LIBRARIES "stdc++")
            SET(${_NAMEVAR}_LDFLAGS "${${_NAMEVAR}_LDFLAGS}-lstdc++")

        ENDIF(${is_extra_flags})
    ENDIF(${_NAMEVAR}_name_library)

    # Add some path to help
    # ---------------------
    UNSET(_lib_env)
    UNSET(_inc_env)
    IF(WIN32)
        STRING(REPLACE ":" ";" _lib_env "$ENV{INCLUDE}")
    ELSE(WIN32)
        STRING(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
        LIST(APPEND _inc_env "${_path_env}")
        STRING(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
        LIST(APPEND _inc_env "${_path_env}")
        STRING(REPLACE ":" ";" _path_env "$ENV{CPATH}")
        LIST(APPEND _inc_env "${_path_env}")
        STRING(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
        LIST(APPEND _inc_env "${_path_env}")
    ENDIF()

    # Define path we have to looking for first
    # ----------------------------------------
    UNSET(CMAKE_PREFIX_PATH)
    LIST(APPEND CMAKE_PREFIX_PATH ${PC_${_NAMEVAR}_INCLUDEDIR})
    LIST(APPEND CMAKE_PREFIX_PATH ${PC_${_NAMEVAR}_INCLUDE_DIRS})
    LIST(APPEND CMAKE_PREFIX_PATH ${${_NAMEVAR}_DIR})
    LIST(APPEND CMAKE_PREFIX_PATH ${_inc_env})

    # Looking for include
    # -------------------
    SET(_all_headers "")
    LIST(APPEND _all_headers ${${_NAMEVAR}_name_include})
    IF(${_NAMEVAR}_name_include)
        FOREACH(_header_file ${_all_headers})
            FIND_PATH(${_header_file}_path
                         NAMES ${_header_file}
                         PATH_SUFFIXES include ${${_NAMEVAR}_name_include_suffix} include/${${_NAMEVAR}_name_include_suffix}
                        )
            MARK_AS_ADVANCED(${_header_file}_path)
            IF(${_header_file}_path)
                LIST(APPEND ${_NAMEVAR}_INCLUDE_PATH "${${_header_file}_path}")
                LIST(REMOVE_DUPLICATES ${_NAMEVAR}_INCLUDE_PATH)

            ELSE(${_header_file}_path)
                SET(${_NAMEVAR}_ERROR_OCCURED TRUE)
                MESSAGE(STATUS "Looking for ${_header_file} - not found")

            ENDIF(${_header_file}_path)
        ENDFOREACH()

    ENDIF(${_NAMEVAR}_name_include)

    # Define path we have to looking for first
    # ----------------------------------------
    UNSET(CMAKE_PREFIX_PATH)

    # Looking for binary
    # ------------------
    SET(_all_binary "")
    LIST(APPEND _all_binary ${${_NAMEVAR}_name_binary})
    IF(${_NAMEVAR}_name_binary)

        FOREACH(_bin_file ${_all_binary})
            FIND_PATH(${_bin_file}_path
                         NAMES ${_bin_file}
                         PATHS ${${_NAMEVAR}_PATH} ${${_NAMEVAR}_DIR}
                         PATH_SUFFIXES bin bin64 bin32
                        )
            MARK_AS_ADVANCED(${_bin_file}_path)
            IF(${_bin_file}_path)
                LIST(APPEND ${_NAMEVAR}_BINARY_PATH "${${_bin_file}_path}")
                LIST(REMOVE_DUPLICATES ${_NAMEVAR}_BINARY_PATH)

            ELSE(${_bin_file}_path)
                SET(${_NAMEVAR}_ERROR_OCCURED TRUE)
                MESSAGE(STATUS "Looking for ${_bin_file} - not found")

            ENDIF(${_bin_file}_path)
        ENDFOREACH()

    ENDIF(${_NAMEVAR}_name_binary)

    # Print the results obtained
    # --------------------------
    DEBUG_FIND_OUTPUT(${_NAMEVAR})

    # Check if the library works
    # --------------------------
    INCLUDE(checkPACKAGE)
    CHECK_PACKAGE(${_NAMEVAR})
    IF(${_NAMEVAR}_ERROR_OCCURRED)
        PACKAGE_CLEAN_VALUE(${_NAMEVAR})
        MESSAGE(STATUS "Looking for ${_NAMEVAR} - not working")
        SET(${_NAMEVAR}_FOUND FALSE)

    ELSE(${_NAMEVAR}_ERROR_OCCURRED)
        SET(${_NAMEVAR}_USED_MODE "FIND")
        SET(${_NAMEVAR}_FOUND     TRUE  )

    ENDIF(${_NAMEVAR}_ERROR_OCCURRED)

ENDMACRO(FIND_MY_PACKAGE)

###
#
#
#
###
MACRO(PACKAGE_CLEAN_VALUE _NAMEVAR)

    SET(${_NAMEVAR}_VERSION          "${_NAMEVAR}_VERSION-NOTFOUND"     )
    SET(${_NAMEVAR}_PATH             "${_NAMEVAR}_PATH-NOTFOUND"        )
    IF(${_NAMEVAR}_name_library)
        SET(${_NAMEVAR}_LDFLAGS      "${_NAMEVAR}_LDFLAGS-NOTFOUND"     )
        SET(${_NAMEVAR}_LIBRARY      "${_NAMEVAR}_LIBRARY-NOTFOUND"     )
        SET(${_NAMEVAR}_LIBRARIES    "${_NAMEVAR}_LIBRARIES-NOTFOUND"   )
        SET(${_NAMEVAR}_LIBRARY_PATH "${_NAMEVAR}_LIBRARY_PATH-NOTFOUND")
    ENDIF(${_NAMEVAR}_name_library)
    IF(${_NAMEVAR}_name_binary)
        SET(${_NAMEVAR}_BINARY_PATH  "${_NAMEVAR}_BINARY_PATH-NOTFOUND" )
    ENDIF(${_NAMEVAR}_name_binary)
    IF(${_NAMEVAR}_name_include)
        SET(${_NAMEVAR}_INCLUDE_PATH "${_NAMEVAR}_INCLUDE_PATH-NOTFOUND")
    ENDIF(${_NAMEVAR}_name_include)

ENDMACRO(PACKAGE_CLEAN_VALUE)

###
#
# FIND_MY_PACAKGE_DEBUG: print all the results
#
###
MACRO(DEBUG_FIND_OUTPUT _NAMEVAR)

    IF(MAGMA_DEBUG_CMAKE)

        # Load infoPACKAGE
        # ----------------
        #INCLUDE(infoPACKAGE)
        #INFO_FIND_PACKAGE(${_NAMEVAR})

        # Debug post-treatment
        # --------------------
        MESSAGE(STATUS "  * debug:")
        MESSAGE(STATUS "  * debug: Find${_NAMEVAR} - obtained values")
        MESSAGE(STATUS "  * debug:")
        MESSAGE(STATUS "  * debug: ${_NAMEVAR}_ERROR_OCCURRED = ${${_NAMEVAR}_ERROR_OCCURRED}")
        MESSAGE(STATUS "  * debug: ${_NAMEVAR}_VERSION        = ${${_NAMEVAR}_VERSION}"       )
        MESSAGE(STATUS "  * debug: ${_NAMEVAR}_PATH           = ${${_NAMEVAR}_PATH}"          )
        IF(${_NAMEVAR}_name_library)
            MESSAGE(STATUS "  * debug: ${_NAMEVAR}_LIBRARY_PATH   = ${${_NAMEVAR}_LIBRARY_PATH}")
            MESSAGE(STATUS "  * debug: ${_NAMEVAR}_LDFLAGS        = ${${_NAMEVAR}_LDFLAGS}"     )
            MESSAGE(STATUS "  * debug: ${_NAMEVAR}_LIBRARY        = ${${_NAMEVAR}_LIBRARY}"     )
            MESSAGE(STATUS "  * debug: ${_NAMEVAR}_LIBRARIES      = ${${_NAMEVAR}_LIBRARIES}"   )
        ENDIF(${_NAMEVAR}_name_library)
        IF(${_NAMEVAR}_name_include)
            MESSAGE(STATUS "  * debug: ${_NAMEVAR}_INCLUDE_PATH   = ${${_NAMEVAR}_INCLUDE_PATH}")
        ENDIF(${_NAMEVAR}_name_include)
        IF(${_NAMEVAR}_name_binary)
            MESSAGE(STATUS "  * debug: ${_NAMEVAR}_BINARY_PATH    = ${${_NAMEVAR}_BINARY_PATH}" )
        ENDIF(${_NAMEVAR}_name_binary)
        MESSAGE(STATUS "  * debug:")

    ENDIF(MAGMA_DEBUG_CMAKE)

ENDMACRO(DEBUG_FIND_OUTPUT)

###
#
# DEBUG_PKG_CONFIG_OUTPUT: print all the informations collect by `pkg_search_module'
#
###
MACRO(DEBUG_PKG_CONFIG_OUTPUT _NAMEVAR)

    IF(MAGMA_DEBUG_CMAKE)

        # Load infoPACKAGE
        # ----------------
        #INCLUDE(infoPACKAGE)
        #INFO_FIND_PACKAGE(${_NAMEVAR})

        # Debug PKG-CONFIG
        # ----------------
        MESSAGE(STATUS "  * debug:")
        MESSAGE(STATUS "  * debug: Find${_NAMEVAR} - output form pkg-config for ${${_NAMEVAR}_name_pkgconfig}")
        MESSAGE(STATUS "  * debug:")
        MESSAGE(STATUS "  * debug: PC_${_NAMEVAR}_FOUND        = ${PC_${_NAMEVAR}_FOUND}"       )
        MESSAGE(STATUS "  * debug: PC_${_NAMEVAR}_VERSION      = ${PC_${_NAMEVAR}_VERSION}"     )
        MESSAGE(STATUS "  * debug: PC_${_NAMEVAR}_PREFIX       = ${PC_${_NAMEVAR}_PREFIX}"      )
        IF(${_NAMEVAR}_name_library)
            MESSAGE(STATUS "  * debug: PC_${_NAMEVAR}_LIBRARY_DIRS = ${PC_${_NAMEVAR}_LIBRARY_DIRS}")
            MESSAGE(STATUS "  * debug: PC_${_NAMEVAR}_LIBDIR       = ${PC_${_NAMEVAR}_LIBDIR}"      )
        ENDIF(${_NAMEVAR}_name_library)
        IF(${_NAMEVAR}_name_include)
            MESSAGE(STATUS "  * debug: PC_${_NAMEVAR}_INCLUDEDIR   = ${PC_${_NAMEVAR}_INCLUDEDIR}"  )
            MESSAGE(STATUS "  * debug: PC_${_NAMEVAR}_INCLUDE_DIRS = ${PC_${_NAMEVAR}_INCLUDE_DIRS}")
        ENDIF(${_NAMEVAR}_name_include)
        MESSAGE(STATUS "  * debug:")

    ENDIF(MAGMA_DEBUG_CMAKE)

ENDMACRO(DEBUG_PKG_CONFIG_OUTPUT)

##
## @end file findPACKAGE.cmake
##
