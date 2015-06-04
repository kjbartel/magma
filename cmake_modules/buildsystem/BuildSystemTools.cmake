###
#
#  @file BuildSystemTools.cmake 
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

###
#
# GET_VERSION: Get the version of the software by parsing a file
#
###

MACRO(GET_VERSION _PACKAGE _filepath)

    FILE(READ "${_filepath}" _file)
    STRING(REGEX REPLACE "(.*)define([ \t]*)${_PACKAGE}_VERSION_MAJOR([ \t]*)([0-9]+)(.*)" "\\4" ${_PACKAGE}_VERSION_MAJOR "${_file}")
    STRING(REGEX REPLACE "(.*)define([ \t]*)${_PACKAGE}_VERSION_MINOR([ \t]*)([0-9]+)(.*)" "\\4" ${_PACKAGE}_VERSION_MINOR "${_file}")
    STRING(REGEX REPLACE "(.*)define([ \t]*)${_PACKAGE}_VERSION_MICRO([ \t]*)([0-9]+)(.*)" "\\4" ${_PACKAGE}_VERSION_PATCH "${_file}")
    SET(${_PACKAGE}_VERSION_NUMBER "${${_PACKAGE}_VERSION_MAJOR}.${${_PACKAGE}_VERSION_MINOR}.${${_PACKAGE}_VERSION_PATCH}")
    #MESSAGE(STATUS "${_PACKAGE}_VERSION_MAJOR = -${${_PACKAGE}_VERSION_MAJOR}-")
    #MESSAGE(STATUS "${_PACKAGE}_VERSION_MINOR = -${${_PACKAGE}_VERSION_MINOR}-")
    #MESSAGE(STATUS "${_PACKAGE}_VERSION_MICRO = -${${_PACKAGE}_VERSION_MICRO}-")

ENDMACRO(GET_VERSION)


###
#
# ADD_FLAGS: Add _FLAGS in _VARIABLE with correct whitespace.
#
###
MACRO(ADD_FLAGS _VARIABLE _FLAGS)

    IF("${${_VARIABLE}}" STREQUAL "")
        SET(${_VARIABLE} "${_FLAGS}")
    ELSE()
        SET(${_VARIABLE} "${${_VARIABLE}} ${_FLAGS}")
    ENDIF()

ENDMACRO(ADD_FLAGS _VARIABLE _FLAGS)

###
#
# ADD_FLAGS: Add _FLAGS in _VARIABLE with correct whitespace.
#
###
MACRO(DEFINE_LIBRARY _NAME)

    MESSAGE(STATUS "Looking for lib${_NAME}")
    STRING(TOUPPER "${_NAME}" _UPPERNAME)
    STRING(REPLACE "+" "P" _UPPERNAME "${_UPPERNAME}")
    FIND_LIBRARY(${_UPPERNAME}_LIBRARY
                 NAME ${_NAME}
                 PATHS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
    MARK_AS_ADVANCED(${_UPPERNAME}_LIBRARY)

    IF(${_UPPERNAME}_LIBRARY)
        MESSAGE(STATUS "Looking for lib${_NAME} - found")
        GET_FILENAME_COMPONENT(_libpath ${${_UPPERNAME}_LIBRARY} PATH)
        SET(${_UPPERNAME}_LDFLAGS      "-L${_libpath} -l${_NAME}")
        SET(${_UPPERNAME}_LIST_LDFLAGS "-L${_libpath};-l${_NAME}")
        SET(${_UPPERNAME}_LIBRARY_PATH "${_libpath}")
        SET(${_UPPERNAME}_LIBRARIES    "${_NAME}")
        SET(${_UPPERNAME}_FOUND        TRUE)
        SET(HAVE_${_UPPERNAME}         ON)
        LIST(APPEND MAGMA_EXTRA_LIBS "-l${_NAME}")

    ELSE(${_UPPERNAME}_LIBRARY)
        MESSAGE(STATUS "Looking for lib${_NAME} - not found")

    ENDIF(${_UPPERNAME}_LIBRARY)

ENDMACRO(DEFINE_LIBRARY)

###
#
# GENERATE_PKGCONFIG_FILE: generate a file .pc according to the options 
#
###
MACRO(GENERATE_PKGCONFIG_FILE _PACKAGE _file) 

    # Clear entry
    # -----------
    SET(${_PACKAGE}_REQUIRED "")

    # Define required package
    # -----------------------
    IF(MAGMA)
        IF(HAVE_PLASMA)
            SET(MAGMA_REQUIRED "${MAGMA_REQUIRED} plasma")
        ENDIF()
    ENDIF()

    # Create .pc file
    # ---------------
    GET_FILENAME_COMPONENT(_output_file ${_file} NAME)
    STRING(REPLACE ".in" "" _output_file "${_output_file}")
    SET(_output_file "${CMAKE_BINARY_DIR}/${_output_file}")
    CONFIGURE_FILE("${_file}" "${_output_file}" @ONLY)

    # installation
    # ------------
    INSTALL(FILES ${_output_file} DESTINATION lib/pkgconfig)

ENDMACRO(GENERATE_PKGCONFIG_FILE)

##
## @end file BuildSystemTools.cmake
##
