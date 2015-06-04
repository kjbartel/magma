###
#
#  @file SetDependencies.cmake
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

# Macro to set manually howto link
# --------------------------------
MACRO(SET_PACKAGE _NAME)

    # Set internal variable
    # ---------------------
    STRING(TOUPPER "${_NAME}" _NAMEVAR)

    # Clean internal variables
    # ------------------------
    UNSET(${_NAMEVAR}_LDFLAGS)
    UNSET(${_NAMEVAR}_LIBRARY_PATH)
    UNSET(${_NAMEVAR}_LIBRARIES)
    UNSET(${_NAMEVAR}_INCLUDE_PATH)

    # Load infoPACKAGE
    # ----------------
    INCLUDE(infoPACKAGE)
    INFO_FIND_PACKAGE(${_NAMEVAR})

    # Check info<PACKAGE> to know if there is library to link
    # -------------------------------------------------------
    # Check if <PACKAGE>_LIB was defined
    IF(DEFINED ${_NAMEVAR}_LIB)
        # Get dir and lib
        STRING(REPLACE " " ";" ${_NAMEVAR}_LIST_LIB "${${_NAMEVAR}_LIB}")
        FOREACH(_param ${${_NAMEVAR}_LIST_LIB})

            # Library directory processing
            IF(_param MATCHES "^-L")
                STRING(REGEX REPLACE "^-L(.*)$" "\\1" _dir "${_param}")
                IF(IS_DIRECTORY ${_dir})
                    SET(${_NAMEVAR}_LDFLAGS "${${_NAMEVAR}_LDFLAGS} -L${_dir}")
                    LIST(APPEND ${_NAMEVAR}_LIBRARY_PATH ${_dir})
                    LINK_DIRECTORIES(${_dir})
                ELSE()
                    MESSAGE(FATAL_ERROR "Setting ${_NAMEVAR} - ${_dir} is not a directory")
                ENDIF()

            # Flag processing
            ELSEIF(_param MATCHES "^-l")
                STRING(REGEX REPLACE "^-l(.*)$" "\\1" _lib "${_param}")
                SET(${_NAMEVAR}_LDFLAGS "${${_NAMEVAR}_LDFLAGS} -l${_lib}")
                LIST(APPEND ${_NAMEVAR}_LIBRARIES ${_lib})

            # Include directory processing
            ELSEIF(_param MATCHES "^-I")
                STRING(REGEX REPLACE "^-I(.*)$" "\\1" _dir "${_param}")
                IF(IS_DIRECTORY ${_dir})
                    LIST(APPEND ${_NAMEVAR}_INCLUDE_PATH ${_dir})
                    INCLUDE_DIRECTORIES(${_dir})
                ELSE()
                    MESSAGE(FATAL_ERROR "Setting ${_NAMEVAR} - ${_dir} is not a directory")
                ENDIF()
            ENDIF()
        ENDFOREACH()

        ELSE(DEFINED ${_NAMEVAR}_LIB)
            MESSAGE(FATAL_ERROR "Setting ${_NAMEVAR} - ${_NAMEVAR}_LIB was not defined")

    ENDIF(DEFINED ${_NAMEVAR}_LIB)

    # Check info<PACKAGE> to know if there is an include to find
    # ----------------------------------------------------------
    # Check if <PACKAGE>_INC was defined
    IF(DEFINED ${_NAMEVAR}_INC)
        IF(IS_DIRECTORY ${${_NAMEVAR}_INC})
            LIST(APPEND ${_NAMEVAR}_INCLUDE_PATH ${${_NAMEVAR}_INC})
            INCLUDE_DIRECTORIES( ${${_NAMEVAR}_INC})
        ELSE()
            MESSAGE(FATAL_ERROR "Setting ${_NAMEVAR} - ${${_NAMEVAR}_INC} is not a directory")
        ENDIF()

    ELSE(DEFINED ${_NAMEVAR}_INC)
        MESSAGE(FATAL_ERROR "Setting ${_NAMEVAR} - ${_NAMEVAR}_INC was not defined")

    ENDIF(DEFINED ${_NAMEVAR}_INC)

    # Test library
    # ------------
    INCLUDE(checkPACKAGE)
    CHECK_PACKAGE(${_NAMEVAR})
    IF(${_NAMEVAR}_ERROR_OCCURRED)
        MESSAGE(FATAL_ERROR "Setting ${_NAMEVAR} - The value that you give is not working")

    ELSE(${_NAMEVAR}_ERROR_OCCURRED)
        MESSAGE(STATUS "Setting ${_NAMEVAR} - worked")
        SET(HAVE_${_NAMEVAR} ON)
        SET(${_NAMEVAR}_SET  ON)
        SET(${_NAMEVAR}_USED_MODE "SET")
        SET(${_NAMEVAR}_VENDOR "user's desire")

    ENDIF(${_NAMEVAR}_ERROR_OCCURRED)

ENDMACRO(SET_PACKAGE)

##
## @end file SetDependencies.cmake
##
