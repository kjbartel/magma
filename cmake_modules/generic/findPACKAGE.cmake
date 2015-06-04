###
#
# @file      : findPACKAGE.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 20-01-2012
# @last modified : mer. 06 juin 2012 17:27:19 CEST
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

INCLUDE(CheckFunctionExists)
INCLUDE(CheckFortranFunctionExists)
INCLUDE(FindPackageHandleStandardArgs)

# Find package
# ------------
MACRO(FIND_MY_PACKAGE _NAMEVAR 
                      _have_binary _have_library _have_include _have_pkgconfig
                      _check_c _check_fortran)

    # Found env var to check
    IF(WIN32)
        SET(_libdir LIB)
    ELSE(WIN32)
        IF(APPLE)
            SET(_libdir DYLD_LIBRARY_PATH)
        ELSE(APPLE)
            SET(_libdir LD_LIBRARY_PATH)
        ENDIF(APPLE)
    ENDIF(WIN32)

    # Set ${_NAMEVAR}_FOUND at FALSE by default
    SET(${_NAMEVAR}_FOUND FALSE)
    SET(${_NAMEVAR}_ERROR_OCCURED FALSE)

    # Optionally use pkg-config to detect include/library paths (if pkg-config is available)
    FIND_PACKAGE(PkgConfig QUIET)
    IF(PKG_CONFIG_EXECUTABLE AND ${_have_pkgconfig})
        pkg_search_module(PC_${_NAMEVAR} QUIET ${${_NAMEVAR}_name_pkgconfig})
    ELSE(PKG_CONFIG_EXECUTABLE AND ${_have_pkgconfig})
        MESSAGE(STATUS "Looking for ${_NAMEVAR} - pkgconfig not used")
    ENDIF(PKG_CONFIG_EXECUTABLE AND ${_have_pkgconfig})

    # looking for version
    IF(PC_${_NAMEVAR}_FOUND)
        SET(${_NAMEVAR}_VERSION ${PC_${_NAMEVAR}_VERSION})
    ENDIF(PC_${_NAMEVAR}_FOUND)

    # looking for top directory
    IF(PC_${_NAMEVAR}_FOUND)
        SET(${_NAMEVAR}_PATH ${PC_${_NAMEVAR}_PREFIX})
    ENDIF(PC_${_NAMEVAR}_FOUND)

    # looking for library
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
    IF(${_have_library})
        FOREACH(_lib_file ${_all_libs})
            FIND_LIBRARY(${_lib_file}_path
                         NAMES ${_lib_file}
                         PATHS ${PC_${_NAMEVAR}_LIBDIR} ${PC_${_NAMEVAR}_LIBRARY_DIRS} ${${_NAMEVAR}_DIR} ENV ${_libdir} 
                         PATH_SUFFIXES lib lib64 lib32
                        )
            IF(${_lib_file}_path)
                GET_FILENAME_COMPONENT(${_lib_file}_filename    ${${_lib_file}_path} NAME   )
                STRING(REGEX REPLACE "(.*)/${${_lib_file}_filename}" "\\1" ${_lib_file}_lib_path "${${_lib_file}_path}")
                LIST(APPEND ${_NAMEVAR}_LIBRARIES    "${_lib_file}")
                LIST(APPEND ${_NAMEVAR}_LIBRARY      "${${_lib_file}_path}")
                LIST(APPEND ${_NAMEVAR}_LIBRARY_PATH "${${_lib_file}_lib_path}")
                LIST(REMOVE_DUPLICATES ${_NAMEVAR}_LIBRARY_PATH)

            ELSE(${_lib_file}_path)
                SET(${_NAMEVAR}_ERROR_OCCURED TRUE)
                MESSAGE(STATUS "Looking for ${_NAMEVAR} : library ${_lib_file} not found")

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
    ENDIF(${_have_library})

    # looking for include
    SET(_all_headers "")
    LIST(APPEND _all_headers ${${_NAMEVAR}_name_include})
    IF(${_have_include})
        FOREACH(_header_file ${_all_headers})
            FIND_PATH(${_header_file}_path
                         NAMES ${_header_file}
                         PATHS ${PC_${_NAMEVAR}_INCLUDEDIR} ${PC_${_NAMEVAR}_INCLUDE_DIRS} ${${_NAMEVAR}_DIR}
                         PATH_SUFFIXES include ${${_NAMEVAR}_name_include_suffix} include/${${_NAMEVAR}_name_include_suffix}
                        )
            IF(${_header_file}_path)
                STRING(REGEX REPLACE "(.*)/${_header_file}" "\\1" ${_header_file}_header_path "${${_header_file}_path}")
                LIST(APPEND ${_NAMEVAR}_INCLUDE_PATH "${${_header_file}_header_path}")
                LIST(REMOVE_DUPLICATES ${_NAMEVAR}_INCLUDE_PATH)

            ELSE(${_header_file}_path)
                SET(${_NAMEVAR}_ERROR_OCCURED TRUE)
                MESSAGE(STATUS "Looking for ${_NAMEVAR} : header ${_header_file} not found")

            ENDIF(${_header_file}_path)
        ENDFOREACH()
    ENDIF(${_have_include})

    # looking for binary
    SET(_all_binary "")
    LIST(APPEND _all_binary ${${_NAMEVAR}_name_binary})
    IF(${_have_binary})

        FOREACH(_bin_file ${_all_binary})
            FIND_PATH(${_bin_file}_path
                         NAMES ${_bin_file}
                         PATHS ${${_NAMEVAR}_PATH} ${${_NAMEVAR}_DIR}
                         PATH_SUFFIXES bin bin64 bin32
                        )
            IF(${_bin_file}_path)
                STRING(REGEX REPLACE "(.*)/${_bin_file}" "\\1" ${_bin_file}_bin_path "${${_bin_file}_path}")
                LIST(APPEND ${_NAMEVAR}_BINARY_PATH "${${_bin_file}_bin_path}")
                LIST(REMOVE_DUPLICATES ${_NAMEVAR}_BINARY_PATH)

            ELSE(${_bin_file}_path)
                SET(${_NAMEVAR}_ERROR_OCCURED TRUE)
                MESSAGE(STATUS "Looking for ${_NAMEVAR} : binary ${_bin_file} not found")

            ENDIF(${_bin_file}_path)
        ENDFOREACH()
    ENDIF(${_have_binary})

    # Check if the library works
    CHECK_MY_PACKAGE(${_NAMEVAR} ${_have_include} ${_check_c} ${_check_fortran})
    IF(NOT ${_NAMEVAR}_C_FOUND AND NOT ${_NAMEVAR}_F_FOUND)
        SET(${_NAMEVAR}_ERROR_OCCURED TRUE)

    ENDIF()

    # Debug
    #FIND_MY_PACAKGE_DEBUG(${_have_binary} ${_have_library} ${_have_include})

    # Warning users on screen
    PACKAGE_HANDLE(${_have_binary} ${_have_library} ${_have_include})

ENDMACRO(FIND_MY_PACKAGE)

MACRO(PACKAGE_CLEAN_VALUE _have_binary _have_library _have_include)

    SET(${_NAMEVAR}_VERSION          "${_NAMEVAR}_VERSION-NOTFOUND"     )
    SET(${_NAMEVAR}_PATH             "${_NAMEVAR}_PATH-NOTFOUND"        )
    IF(${_have_library})
        SET(${_NAMEVAR}_LDFLAGS      "${_NAMEVAR}_LDFLAGS-NOTFOUND"     )
        SET(${_NAMEVAR}_LIBRARY      "${_NAMEVAR}_LIBRARY-NOTFOUND"     )
        SET(${_NAMEVAR}_LIBRARIES    "${_NAMEVAR}_LIBRARIES-NOTFOUND"   )
        SET(${_NAMEVAR}_LIBRARY_PATH "${_NAMEVAR}_LIBRARY_PATH-NOTFOUND")
    ENDIF(${_have_library})
    IF(${_have_binary})
        SET(${_NAMEVAR}_BINARY_PATH  "${_NAMEVAR}_BINARY_PATH-NOTFOUND" )
    ENDIF(${_have_binary})
    IF(${_have_include})
        SET(${_NAMEVAR}_INCLUDE_PATH "${_NAMEVAR}_INCLUDE_PATH-NOTFOUND")
    ENDIF(${_have_include})

ENDMACRO(PACKAGE_CLEAN_VALUE)

MACRO(PACKAGE_HANDLE _have_binary _have_library _have_include)

    # If one error occurred, we clean all variables
    # ---------------------------------------------
    IF(${${_NAMEVAR}_ERROR_OCCURED})
        PACKAGE_CLEAN_VALUE(${_have_binary} ${_have_library} ${_have_include})
    ENDIF(${${_NAMEVAR}_ERROR_OCCURED})

    # handle the QUIETLY and REQUIRED arguments and set ${_NAMEVAR}_FOUND to TRUE if all listed variables are TRUE
    # ------------------------------------------------------------------------------------------------------------
    find_package_handle_standard_args(${_NAMEVAR}
                                      DEFAULT_MSG
                                      ${_NAMEVAR}_PATH 
                                      ${_NAMEVAR}_VERSION
                                     )

    IF(${_have_library})
        find_package_handle_standard_args(${_NAMEVAR}
                                          DEFAULT_MSG
                                          ${_NAMEVAR}_LDFLAGS
                                          ${_NAMEVAR}_LIBRARY
                                          ${_NAMEVAR}_LIBRARIES
                                          ${_NAMEVAR}_LIBRARY_PATH
                                         )
    ENDIF(${_have_library})

    IF(${_have_include})
        find_package_handle_standard_args(${_NAMEVAR}
                                          DEFAULT_MSG
                                          ${_NAMEVAR}_INCLUDE_PATH
                                         )
    ENDIF(${_have_include})

    IF(${_have_binary})
        find_package_handle_standard_args(${_NAMEVAR}
                                          DEFAULT_MSG
                                          ${_NAMEVAR}_BINARY_PATH
                                         )

    ENDIF(${_have_binary})

    # if ${_NAMEVAR}_FOUND is not true, we clean all variables
    # --------------------------------------------------------
    IF(NOT ${_NAMEVAR}_FOUND)
        PACKAGE_CLEAN_VALUE(${_have_binary} ${_have_library} ${_have_include})
    ENDIF()

ENDMACRO(PACKAGE_HANDLE)


# Check library
# -------------
MACRO(CHECK_MY_PACKAGE _NAMEVAR _have_include _check_c _check_f)

    SET(${_NAMEVAR}_C_FOUND FALSE)
    SET(${_NAMEVAR}_F_FOUND FALSE)

    SET(CMAKE_REQUIRED_FLAGS ${${_NAMEVAR}_LDFLAGS})

    IF(${_have_include})
        SET(CMAKE_REQUIRED_INCLUDES ${${_NAMEVAR}_INCLUDE_PATH})
    ENDIF(${_have_include})

    IF(${_check_c})
        check_function_exists(${${_NAMEVAR}_name_fct_test} _status_C)
        IF(_status_C)
            SET(${_NAMEVAR}_C_FOUND TRUE)
        ENDIF(_status_C)
    ENDIF(${_check_c})

    IF(${_check_f})
        check_fortran_function_exists(${${_NAMEVAR}_name_fct_test} _status_F)
        IF(_status_F)
            SET(${_NAMEVAR}_F_FOUND TRUE)
        ENDIF(_status_F)
    ENDIF(${_check_f})

ENDMACRO(CHECK_MY_PACKAGE)

MACRO(FIND_MY_PACAKGE_DEBUG _have_binary _have_library _have_include)

    # Info PACKAGE
    # ------------
    MESSAGE(STATUS "DEBUG: Detection of ${_NAMEVAR}")
    MESSAGE(STATUS "DEBUG: ${_NAMEVAR} - have_library = ${_have_library}")
    MESSAGE(STATUS "DEBUG: ${_NAMEVAR} - have_include = ${_have_include}")
    MESSAGE(STATUS "DEBUG: ${_NAMEVAR} - have_binary  = ${_have_binary}")

    # Debug PKG-CONFIG
    # ----------------
    MESSAGE(STATUS "DEBUG: Output form PKG-CONFIG")
    MESSAGE(STATUS "DEBUG: PC_${_NAMEVAR}_FOUND        = ${PC_${_NAMEVAR}_FOUND}")
    MESSAGE(STATUS "DEBUG: PC_${_NAMEVAR}_VERSION      = ${PC_${_NAMEVAR}_VERSION}")
    MESSAGE(STATUS "DEBUG: PC_${_NAMEVAR}_PREFIX       = ${PC_${_NAMEVAR}_PREFIX}")
    IF(${_have_library})
        MESSAGE(STATUS "DEBUG: PC_${_NAMEVAR}_LIBRARY_DIRS = ${PC_${_NAMEVAR}_LIBRARY_DIRS}")
        MESSAGE(STATUS "DEBUG: PC_${_NAMEVAR}_LIBDIR       = ${PC_${_NAMEVAR}_LIBDIR}")
    ENDIF(${_have_library})
    IF(${_have_include})
        MESSAGE(STATUS "DEBUG: PC_${_NAMEVAR}_INCLUDEDIR   = ${PC_${_NAMEVAR}_INCLUDEDIR}")
        MESSAGE(STATUS "DEBUG: PC_${_NAMEVAR}_INCLUDE_DIRS = ${PC_${_NAMEVAR}_INCLUDE_DIRS}")
    ENDIF(${_have_include})

    # Debug post-treatment
    # --------------------
    MESSAGE(STATUS "DEBUG: Value for ${_NAMEVAR}")
    MESSAGE(STATUS "DEBUG: ${_NAMEVAR}_VERSION        = ${${_NAMEVAR}_VERSION}")
    MESSAGE(STATUS "DEBUG: ${_NAMEVAR}_PATH           = ${${_NAMEVAR}_PATH}")
    IF(${_have_library})
        MESSAGE(STATUS "DEBUG: ${_NAMEVAR}_LIBRARY_PATH   = ${${_NAMEVAR}_LIBRARY_PATH}")
        MESSAGE(STATUS "DEBUG: ${_NAMEVAR}_LDFLAGS        = ${${_NAMEVAR}_LDFLAGS}"     )
        MESSAGE(STATUS "DEBUG: ${_NAMEVAR}_LIBRARY        = ${${_NAMEVAR}_LIBRARY}"     )
        MESSAGE(STATUS "DEBUG: ${_NAMEVAR}_LIBRARIES      = ${${_NAMEVAR}_LIBRARIES}"   )
    ENDIF(${_have_library})
    IF(${_have_include})
        MESSAGE(STATUS "DEBUG: ${_NAMEVAR}_INCLUDE_PATH   = ${${_NAMEVAR}_INCLUDE_PATH}")
    ENDIF(${_have_include})
    IF(${_have_binary})
        MESSAGE(STATUS "DEBUG: ${_NAMEVAR}_BINARY_PATH    = ${${_NAMEVAR}_BINARY_PATH}")
    ENDIF(${_have_binary})

    # Debug tests
    # -----------
    MESSAGE(STATUS "${_NAMEVAR}_C_FOUND               = ${${_NAMEVAR}_C_FOUND}")
    MESSAGE(STATUS "${_NAMEVAR}_F_FOUND               = ${${_NAMEVAR}_F_FOUND}")

ENDMACRO(FIND_MY_PACAKGE_DEBUG)


###
### END findPACKAGE.cmake
###
