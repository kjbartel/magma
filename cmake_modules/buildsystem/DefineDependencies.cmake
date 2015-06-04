###
#
#  @file DefineDependencies.cmake
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
# * Users variables to help MAGMA dependencies system
#   - ${_UPPERWANTED}_LIB     :
#   - ${_UPPERWANTED}_INC     :
#   - ${_UPPERWANTED}_DIR     :
#   - ${_UPPERWANTED}_URL     :
#   - ${_UPPERWANTED}_TARBALL :
#
# * Variables of control the policy about MAGMA dependencies
#   - ${_UPPERWANTED}_USE_LIB     :
#   - ${_UPPERWANTED}_USE_SYSTEM  : 
#   - ${_UPPERWANTED}_USE_TARBALL : 
#   - ${_UPPERWANTED}_USE_WEB     :
#   - ${_UPPERWANTED}_USE_AUTO    :
#
###

INCLUDE(populatePACKAGE)
INCLUDE(SetDependencies)

# Macro to define if you have to download or use a package
# --------------------------------------------------------
MACRO(DEFINE_PACKAGE _WANTED)

    # Define input variable to uppercase
    # ----------------------------------
    STRING(TOUPPER "${_WANTED}" _UPPERWANTED)
    STRING(TOLOWER "${_WANTED}" _LOWERWANTED)


    # Process to define if you have to download or use a package
    # ----------------------------------------------------------
    UNSET(HAVE_${_UPPERWANTED})

    # Define manually processing
    IF(DEFINED ${_NAMEVAR}_LIB OR DEFINED ${_NAMEVAR}_INC)
        MESSAGE(STATUS "Setting ${_UPPERWANTED}")
        SET_PACKAGE(${_UPPERWANTED})

    # Find in the system processing
    ELSE()
        MESSAGE(STATUS "Looking for ${_UPPERWANTED}")
        FIND_PACKAGE(${_UPPERWANTED} REQUIRED)
        IF(${_UPPERWANTED}_FOUND)
            POPULATE_COMPILE_SYSTEM(${_UPPERWANTED})
        ELSE()
            MESSAGE(FATAL_ERROR "Looking for ${_UPPERWANTED} - not found")
        ENDIF()

    ENDIF()

ENDMACRO(DEFINE_PACKAGE)

##
## @end file DefineDependencies.cmake
##
