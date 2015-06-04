###
#
# @file      : update_env.cmake 
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : sam. 21 janv. 2012 19:05:24 CET
#
###

#
#
# @file search_library.cmake
#
#  Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#  (experimental version)
#
# @version 1.2.0
# @author Cedric Castagnede
# @date 2011-12-22
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Macro to update user's environnement
# ------------------------------------
MACRO(UPDATE_ENV _VARNAME _ADD)
    STRING(REGEX MATCH "${_ADD}" _STATUS "$ENV{${_VARNAME}}")
    IF(DEFINED _STATUS)
        SET(ENV{${_VARNAME}} "${_ADD}:$ENV{${_VARNAME}}")
    ENDIF()
ENDMACRO(UPDATE_ENV)

