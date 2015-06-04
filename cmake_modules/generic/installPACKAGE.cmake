###
#
# @file      : installPACKAGE.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : mer. 16 mai 2012 10:12:38 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(ExternalProject)

MACRO(INSTALL_PACKAGE _NAME _MODE)

    STRING(TOUPPER "${_NAME}" _NAMEVAR)
    MESSAGE(STATUS "Installing ${_NAMEVAR}")

    # Define dependencies
    # -------------------
    SET(${_NAMEVAR}_MAKEDEP "")
    FOREACH(_dep ${${_NAMEVAR}_DEPENDENCIES})
        STRING(TOUPPER "${_dep}" _DEPVAR)
        IF(DEFINED ${_DEPVAR}_BUILD_MODE)
            SET(${_NAMEVAR}_MAKEDEP ${${_NAMEVAR}_MAKEDEP} ${_dep}_install)
        ENDIF(DEFINED ${_DEPVAR}_BUILD_MODE)
    ENDFOREACH()

    # Looking for installation mode
    # -----------------------------
    IF(${_MODE} MATCHES "TARBALL" OR ${_MODE} MATCHES "WEB")
        SET(${_NAMEVAR}_GETSOURCE ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_USED_FILE})
        IF(DEFINED ${_NAMEVAR}_MD5SUM)
            ExternalProject_Add(${_NAME}_install
                    DEPENDS      ${${_NAMEVAR}_MAKEDEP}
                    DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/externals
                    SOURCE_DIR   ${CMAKE_BINARY_DIR}/externals/${_NAME}
                    BINARY_DIR   ${CMAKE_BINARY_DIR}/externals/${_NAME}
                    INSTALL_DIR  ${${_NAMEVAR}_PATH}
                    URL          ${${_NAMEVAR}_GETSOURCE}
                    URL_MD5      ${${_NAMEVAR}_MD5SUM}
                    CONFIGURE_COMMAND ${${_NAMEVAR}_CONFIG_CMD} ${${_NAMEVAR}_OPTIONS} 
                    BUILD_COMMAND ${${_NAMEVAR}_MAKE_CMD} 
                    INSTALL_COMMAND ${${_NAMEVAR}_MAKEINSTALL_CMD})
        ELSE()
            ExternalProject_Add(${_NAME}_install
                    DEPENDS      ${${_NAMEVAR}_MAKEDEP}
                    DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/externals
                    SOURCE_DIR   ${CMAKE_BINARY_DIR}/externals/${_NAME}
                    BINARY_DIR   ${CMAKE_BINARY_DIR}/externals/${_NAME}
                    INSTALL_DIR  ${${_NAMEVAR}_PATH}
                    URL          ${${_NAMEVAR}_GETSOURCE}
                    CONFIGURE_COMMAND ${${_NAMEVAR}_CONFIG_CMD} ${${_NAMEVAR}_OPTIONS} 
                    BUILD_COMMAND ${${_NAMEVAR}_MAKE_CMD} 
                    INSTALL_COMMAND ${${_NAMEVAR}_MAKEINSTALL_CMD})
        ENDIF()
        MESSAGE(STATUS "Installing ${_NAMEVAR} - ${_MODE}")

    ELSEIF(${_MODE} MATCHES "SVN")
        SET(${_NAMEVAR}_GETSOURCE ${${_NAMEVAR}_SVN_REP})
        ExternalProject_Add(${_NAME}_install
                DEPENDS        ${${_NAMEVAR}_MAKEDEP}
                DOWNLOAD_DIR   ${CMAKE_SOURCE_DIR}/externals
                SOURCE_DIR     ${CMAKE_BINARY_DIR}/externals/${_NAME} 
                BINARY_DIR     ${CMAKE_BINARY_DIR}/externals/${_NAME}
                INSTALL_DIR    ${${_NAMEVAR}_PATH}
                SVN_REPOSITORY ${${_NAMEVAR}_GETSOURCE}
                SVN_USERNAME   ${${_NAMEVAR}_SVN_ID}
                SVN_PASSWORD   ${${_NAMEVAR}_SVN_PWD}
                CONFIGURE_COMMAND ${${_NAMEVAR}_CONFIG_CMD} ${${_NAMEVAR}_OPTIONS} 
                BUILD_COMMAND ${${_NAMEVAR}_MAKE_CMD}
                INSTALL_COMMAND ${${_NAMEVAR}_MAKEINSTALL_CMD})
        MESSAGE(STATUS "Installing ${_NAMEVAR} - ${_MODE}")

    ELSEIF(${_MODE} MATCHES "CVS")
        SET(${_NAMEVAR}_GETSOURCE ${${_NAMEVAR}_CVS_REP})
        ExternalProject_Add(${_NAME}_install
                DEPENDS        ${${_NAMEVAR}_MAKEDEP}
                DOWNLOAD_DIR   ${CMAKE_SOURCE_DIR}/externals
                SOURCE_DIR     ${CMAKE_BINARY_DIR}/externals/${_NAME}
                BINARY_DIR     ${CMAKE_BINARY_DIR}/externals/${_NAME}
                INSTALL_DIR    ${${_NAMEVAR}_PATH}}
                CVS_REPOSITORY ${${_NAMEVAR}_GETSOURCE}
                CVS_MODULE     ${${_NAMEVAR}_CVS_MOD}
                CONFIGURE_COMMAND ${${_NAMEVAR}_CONFIG_CMD} ${${_NAMEVAR}_OPTIONS} 
                BUILD_COMMAND ${${_NAMEVAR}_MAKE_CMD}
                INSTALL_COMMAND ${${_NAMEVAR}_MAKEINSTALL_CMD})
        MESSAGE(STATUS "Installing ${_NAMEVAR} - ${_MODE}")

    ELSE()
        MESSAGE(FATA_ERROR "DEV_ERROR - Macro INSTALL_PACKAGE")

    ENDIF()

    # Define additional step of the installation
    # ------------------------------------------
    IF(DEFINED ${_NAMEVAR}_ADD_STEP)
        FOREACH(_STEP ${${_NAMEVAR}_ADD_STEP})
            ExternalProject_Add_Step(${_NAME}_install ${_STEP}
                     COMMAND ${${_STEP}_CMD}
                     WORKING_DIRECTORY ${${_STEP}_DIR}
                     DEPENDEES ${${_STEP}_DEP_BEFORE}
                     DEPENDERS ${${_STEP}_DEP_AFTER})
        ENDFOREACH()
    ENDIF(DEFINED ${_NAMEVAR}_ADD_STEP)

    # Message to advise installation mode
    # -----------------------------------
    #MESSAGE(STATUS "${_NAMEVAR} install mode           : ${${_NAMEVAR}_GETMODE}")
    #MESSAGE(STATUS "${_NAMEVAR} will be installed from : ${${_NAMEVAR}_GETSOURCE}")
    #MESSAGE(STATUS "${_NAMEVAR} will be installed in   : ${CMAKE_INSTALL_PREFIX}/${_NAME}")
    #MESSAGE(STATUS "${_NAMEVAR} configuration options  : ${${_NAMEVAR}_OPTIONS}")

ENDMACRO(INSTALL_PACKAGE)

###
### END installPACKAGE.cmake
###
