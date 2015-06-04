###
#
# @file      : downloadPACKAGE.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : mer. 04 avril 2012 12:03:14 CEST

#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Macro to download tarball
# -------------------------
MACRO(DOWNLOAD_PACKAGE _NAME)

    STRING(TOUPPER "${_NAME}" _NAMEVAR)
    MESSAGE(STATUS "Download: ${${_NAMEVAR}_URL}")
    IF(DEFINED ${_NAMEVAR}_MD5SUM)
        FILE(DOWNLOAD ${${_NAMEVAR}_URL}
             ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_USED_FILE}
             EXPECTED_MD5 ${${_NAMEVAR}_MD5SUM}
             STATUS IS_GOT
             SHOW_PROGRESS
             )
    ELSE()
        MESSAGE(STATUS "Warning: ${_NAMEVAR} checksum is not performed")
        FILE(DOWNLOAD ${${_NAMEVAR}_URL}
             ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_USED_FILE}
             STATUS IS_GOT
             SHOW_PROGRESS
             )
    ENDIF()
    LIST(GET IS_GOT 0 IS_GOT_CODE)
    LIST(GET IS_GOT 1 IS_GOT_MSG)
    IF(IS_GOT_CODE)
        MESSAGE(FATAL_ERROR "ERROR: ${_NAMEVAR} tarball -- not downloaded")
    ENDIF(IS_GOT_CODE)
    MESSAGE(STATUS "Download: ${${_NAMEVAR}_URL} - Done")

ENDMACRO(DOWNLOAD_PACKAGE)

# Macro to download tarball (alternative methode)
# -----------------------------------------------
MACRO(DOWNLOAD_PACKAGE_ALTERNATIVE _NAME)

    STRING(TOUPPER "${_NAME}" _NAMEVAR)

    IF(WIN32)
        DOWNLOAD_PACKAGE("${_NAMEVAR}")
    ENDIF(WIN32)

    IF(APPLE)
        MESSAGE(STATUS "Download: ${${_NAMEVAR}_URL}")
        EXECUTE_PROCESS(
            COMMAND curl ${${_NAMEVAR}_URL} -o ${${_NAMEVAR}_USED_FILE}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/externals/
            )
        MESSAGE(STATUS "Download: ${${_NAMEVAR}_URL} - Done")

    ELSE(APPLE)
        MESSAGE(STATUS "Download: ${${_NAMEVAR}_URL}")
        EXECUTE_PROCESS(
            COMMAND wget ${${_NAMEVAR}_URL} -O ${${_NAMEVAR}_USED_FILE}
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/externals/
            )
        MESSAGE(STATUS "Download: ${${_NAMEVAR}_URL} - Done")

    ENDIF(APPLE)

ENDMACRO(DOWNLOAD_PACKAGE_ALTERNATIVE)

# Macro to know if you have to download tarball
# ---------------------------------------------
MACRO(DEFINE_DOWNLOAD_PACKAGE _NAME _MODE)
    
    # DEBUG
    # -----
    #MESSAGE(STATUS "MACRO DEFINE_DOWNLOAD_PACKAGE : ${_NAME} - ${_MODE}")

    STRING(TOUPPER "${_NAME}" _NAMEVAR)
    UNSET(${_NAMEVAR}_BUILD_MODE)

    IF(${_MODE} MATCHES "TARBALL")
         IF(NOT EXISTS ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_TARBALL})
            MESSAGE(FATAL_ERROR "Looking for ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_TARBALL} - not found")
        ENDIF()
        IF(DEFINED ${_NAMEVAR}_MD5SUM)
            EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E md5sum ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_TARBALL}
                            OUTPUT_VARIABLE ${_NAMEVAR}_TARBALL_MD5_FOUND
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                           )
            IF(NOT "${${_NAMEVAR}_TARBALL_MD5_FOUND}" MATCHES "${${_NAMEVAR}_MD5SUM}")
                MESSAGE(FATAL_ERROR "Looking for ${${_NAMEVAR}_TARBALL} - checksum failed")
            ELSE()
                MESSAGE(STATUS "Looking for ${${_NAMEVAR}_TARBALL} - checksum succeed")
            ENDIF()
        ELSE()
            MESSAGE(STATUS "Warning: ${${_NAMEVAR}_TARBALL} checksum is not performed")
        ENDIF()
        SET(${_NAMEVAR}_BUILD_MODE "TARBALL")
        SET(${_NAMEVAR}_USED_FILE "${${_NAMEVAR}_TARBALL}")

    ELSEIF(${_MODE} MATCHES "WEB")
        STRING(REGEX REPLACE ".*/(.*)$" "\\1" ${_NAMEVAR}_USED_FILE "${${_NAMEVAR}_URL}")
        IF(EXISTS ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_TARBALL})
            FILE(REMOVE ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_TARBALL})
        ENDIF()
        DOWNLOAD_PACKAGE(${_NAMEVAR})
        SET(${_NAMEVAR}_BUILD_MODE "WEB")

    ELSEIF(${_MODE} MATCHES "AUTO")
        IF(NOT EXISTS ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_TARBALL})
            STRING(REGEX REPLACE ".*/(.*)$" "\\1" ${_NAMEVAR}_USED_FILE "${${_NAMEVAR}_URL}")
            DOWNLOAD_PACKAGE(${_NAMEVAR})
            SET(${_NAMEVAR}_BUILD_MODE "WEB")
        ELSE()
            IF(DEFINED ${_NAMEVAR}_MD5SUM)
                EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E md5sum ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_TARBALL}
                                OUTPUT_VARIABLE ${_NAMEVAR}_TARBALL_MD5_FOUND
                                OUTPUT_STRIP_TRAILING_WHITESPACE
                               )
                IF(NOT "${${_NAMEVAR}_TARBALL_MD5_FOUND}" MATCHES "${${_NAMEVAR}_MD5SUM}")
                    MESSAGE(FATAL_ERROR "Looking for ${${_NAMEVAR}_TARBALL} - checksum failed")
                ELSE()
                    MESSAGE(STATUS "Looking for ${${_NAMEVAR}_TARBALL} - checksum succeed")
                ENDIF()
            ELSE()
                MESSAGE(STATUS "Warning: ${${_NAMEVAR}_TARBALL} checksum is not performed")
            ENDIF()
            SET(${_NAMEVAR}_BUILD_MODE "TARBALL")
            SET(${_NAMEVAR}_USED_FILE "${${_NAMEVAR}_TARBALL}")
        ENDIF()

    ELSEIF(${_MODE} MATCHES "REPO")
       SET(${_NAMEVAR}_BUILD_MODE "${${_NAMEVAR}_REPO_MODE}")

    ELSE()
        MESSAGE(FATAL_ERROR "DEV_ERROR - Macro DEFINE_DOWNLOAD_PACKAGE")

    ENDIF()
    
    # DEBUG
    #MESSAGE(STATUS "MACRO DEFINE_DOWNLOAD_PACKAGE : ${_NAME} - ${${_NAMEVAR}_BUILD_MODE}")
    #MESSAGE(STATUS "MACRO DEFINE_DOWNLOAD_PACKAGE : ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_TARBALL}")

ENDMACRO(DEFINE_DOWNLOAD_PACKAGE)

###
### END downloadPACKAGE.cmake
###
