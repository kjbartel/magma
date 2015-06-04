###
#
# @file          : definePACKAGE.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 09-03-2012
# @last modified : mar. 12 juin 2012 16:43:45 CEST
#
###

# Macro to define test functions 
# ------------------------------
MACRO(INFO_PACKAGE)

    INCLUDE(infoBLAS)
    BLAS_INFO_INSTALL()
    INCLUDE(infoLAPACK)
    LAPACK_INFO_INSTALL()
    INCLUDE(infoCBLAS)
    CBLAS_INFO_INSTALL()
    INCLUDE(infoLAPACKE)
    LAPACKE_INFO_INSTALL()
    INCLUDE(infoPLASMA)
    PLASMA_INFO_INSTALL()
    INCLUDE(infoOPENMPI)
    MPI_INFO_INSTALL()
    INCLUDE(infoFXT)
    FXT_INFO_INSTALL()
    INCLUDE(infoHWLOC)
    HWLOC_INFO_INSTALL()
    INCLUDE(infoSTARPU)
    STARPU_INFO_INSTALL()

ENDMACRO(INFO_PACKAGE)

# Macro to call installation
# --------------------------
MACRO(MORSE_INSTALL_PACAKGE _NAME _MODE)

    STRING(TOUPPER "${_NAME}" _NAMEVAR)
    IF("BLAS" MATCHES "${_NAMEVAR}")
        INCLUDE(installBLAS)
        INSTALL_BLAS("${_MODE}")

    ELSEIF("LAPACK" MATCHES "${_NAMEVAR}")
        INCLUDE(installLAPACK)
        INSTALL_LAPACK("${_MODE}")

    ELSEIF("CBLAS" MATCHES "${_NAMEVAR}")
        INCLUDE(installCBLAS)
        INSTALL_CBLAS("${_MODE}")

    ELSEIF("LAPACKE" MATCHES "${_NAMEVAR}")
        INCLUDE(installLAPACKE)
        INSTALL_LAPACKE("${_MODE}")

    ELSEIF("PLASMA" MATCHES "${_NAMEVAR}")
        INCLUDE(installPLASMA)
        INSTALL_PLASMA("${_MODE}")

    ELSEIF("MPI" MATCHES "${_NAMEVAR}")
        INCLUDE(installOPENMPI)
        INSTALL_OPENMPI("${_MODE}")

    ELSEIF("HWLOC" MATCHES "${_NAMEVAR}")
        INCLUDE(installHWLOC)
        INSTALL_HWLOC("${_MODE}")

    ELSEIF("FXT" MATCHES "${_NAMEVAR}")
        INCLUDE(installFXT)
        INSTALL_FXT("${_MODE}")

    ELSEIF("STARPU" MATCHES "${_NAMEVAR}")
        INCLUDE(installSTARPU)
        INSTALL_STARPU("${_MODE}")

    ELSE()
        MESSAGE(FATAL_ERROR "DEV_ERROR - Macro MORSE_INSTALL_PACAKGE")

    ENDIF()

ENDMACRO(MORSE_INSTALL_PACAKGE)

# Macro to define if you have to download or use a package
# --------------------------------------------------------
MACRO(DEFINE_PACKAGE _NAME _TYPE_DEP)

    STRING(TOUPPER "${_NAME}" _NAMEVAR)

    # Default options
    # ---------------
    IF(${_TYPE_DEP} MATCHES "suggests")
        SET(${_NAMEVAR}_TRYUSE OFF)
        SET(${_NAMEVAR}_REQUIRED OFF)

    ELSEIF(${_TYPE_DEP} MATCHES "recommands")
        SET(${_NAMEVAR}_TRYUSE ON)
        SET(${_NAMEVAR}_REQUIRED OFF)

    ELSEIF(${_TYPE_DEP} MATCHES "depends")
        SET(${_NAMEVAR}_TRYUSE ON)
        SET(${_NAMEVAR}_REQUIRED ON)

    ELSE()
        MESSAGE(FATAL_ERROR "DEV_ERROR - Macro DEFINE_PACKAGE")

    ENDIF()
    MARK_AS_ADVANCED(${_NAMEVAR}_TRYUSE ${_NAMEVAR}_REQUIRED)

    # Switch options according to user's desires
    # ------------------------------------------
    SET(MORSE_USE_${_NAMEVAR} ""
        CACHE STRING "Enable/Disable ${_NAMEVAR} dependency (ON/OFF/<not-defined>)")

    IF("${MORSE_USE_${_NAMEVAR}}" MATCHES "OFF")
        SET(${_NAMEVAR}_TRYUSE OFF)
        SET(${_NAMEVAR}_REQUIRED OFF)
    ELSEIF(NOT "${MORSE_USE_${_NAMEVAR}}" MATCHES "")
        SET(${_NAMEVAR}_TRYUSE ON)
        SET(${_NAMEVAR}_REQUIRED ON)
    ENDIF()

    # Defaut options for a package
    # ----------------------------
    IF(DEFINED ${_NAMEVAR}_LIB OR ${_NAMEVAR}_USE_LIB)
        OPTION(${_NAMEVAR}_USE_LIB "Enable/Disable to link with ${_NAMEVAR}" ON)
        OPTION(${_NAMEVAR}_USE_AUTO "Enable/Disable to install ${_NAMEVAR} if it was not found" OFF)
        SET(${_NAMEVAR}_TRYUSE ON)
        SET(${_NAMEVAR}_REQUIRED ON)
    ELSE()
        OPTION(${_NAMEVAR}_USE_LIB "Enable/Disable to link with ${_NAMEVAR}" OFF)
        OPTION(${_NAMEVAR}_USE_AUTO "Enable/Disable to install ${_NAMEVAR} if it was not found" ON)
    ENDIF()

    IF(DEFINED ${_NAMEVAR}_URL OR ${_NAMEVAR}_USE_WEB)
        OPTION(${_NAMEVAR}_USE_WEB "Enable/Disable to install ${_NAMEVAR} form ${${_NAMEVAR}_URL}" ON)
        OPTION(${_NAMEVAR}_USE_AUTO "Enable/Disable to install ${_NAMEVAR} if it was not found" OFF)
        SET(${_NAMEVAR}_TRYUSE ON)
        SET(${_NAMEVAR}_REQUIRED ON)
    ELSE()
        OPTION(${_NAMEVAR}_USE_WEB "Enable/Disable to install ${_NAMEVAR} form ${${_NAMEVAR}_URL}" OFF)
        OPTION(${_NAMEVAR}_USE_AUTO "Enable/Disable to install ${_NAMEVAR} if it was not found" ON)
    ENDIF()

    IF(DEFINED ${_NAMEVAR}_TARBALL OR ${_NAMEVAR}_USE_TARBALL)
        OPTION(${_NAMEVAR}_USE_TARBALL "Enable/Disable to install ${_NAMEVAR} form ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_TARBALL}" ON)
        OPTION(${_NAMEVAR}_USE_AUTO "Enable/Disable to install ${_NAMEVAR} if it was not found" OFF)
        SET(${_NAMEVAR}_TRYUSE ON)
        SET(${_NAMEVAR}_REQUIRED ON)
    ELSE()
        OPTION(${_NAMEVAR}_USE_TARBALL "Enable/Disable to install ${_NAMEVAR} form ${CMAKE_SOURCE_DIR}/externals/${${_NAMEVAR}_TARBALL}" OFF)
        OPTION(${_NAMEVAR}_USE_AUTO "Enable/Disable to install ${_NAMEVAR} if it was not found" ON)
    ENDIF()

    OPTION(${_NAMEVAR}_USE_SYSTEM "Enable/Disable to look for ${_NAMEVAR} installation in environment" OFF)
    OPTION(${_NAMEVAR}_USE_SVN "Enable/Disable to install ${_NAMEVAR} form the repository" OFF)

    # Force installation for eigen
    # ----------------------------
    IF("${_NAMEVAR}" MATCHES "BLAS")
        STRING(TOUPPER "${MORSE_USE_BLAS}" VALUE_MORSE_USE_BLAS)
        IF("${VALUE_MORSE_USE_BLAS}" MATCHES "EIGEN" OR
           "${VALUE_MORSE_USE_BLAS}" MATCHES "REFBLAS")
            SET(BLAS_TRYUSE      ON )
            SET(BLAS_REQUIRED    ON )
            SET(BLAS_USE_LIB     OFF)
            SET(BLAS_USE_SYSTEM  OFF)
            SET(BLAS_USE_TARBALL ON )
            SET(BLAS_USE_WEB     ON )
            SET(BLAS_USE_AUTO    OFF)
        ENDIF()
    ENDIF()

    # Status of internal variables
    # ----------------------------
    #MESSAGE(STATUS "Status of ${_NAMEVAR}:")
    #MESSAGE(STATUS "  --> ${_NAMEVAR}_TRYUSE  : ${${_NAMEVAR}_TRYUSE}")
    #MESSAGE(STATUS "  --> ${_NAMEVAR}_REQUIRED: ${${_NAMEVAR}_REQUIRED}")

    # Process to define if you have to download or use a package
    # ----------------------------------------------------------
    UNSET(HAVE_${_NAMEVAR})
    UNSET(${_NAMEVAR}_BUILD_MODE)
    IF(${${_NAMEVAR}_TRYUSE})
        IF(${_NAMEVAR}_USE_LIB)
            MESSAGE(STATUS "Looking for ${_NAMEVAR}")
            IF(DEFINED ${_NAMEVAR}_LIB)
                # Get dir and lib
                UNSET(${_NAMEVAR}_LDFLAGS)
                UNSET(${_NAMEVAR}_LIBRARY_PATH)
                UNSET(${_NAMEVAR}_LIBRARIES)
                UNSET(${_NAMEVAR}_INCLUDE_PATH)
                SET(${_NAMEVAR}_LDFLAGS "${${_NAMEVAR}_LIB}")
                STRING(REPLACE " " ";" ${_NAMEVAR}_LIST_LIB "${${_NAMEVAR}_LIB}")
                FOREACH(_param ${${_NAMEVAR}_LIST_LIB})
                    IF(_param MATCHES "^-L")
                        STRING(REGEX REPLACE "^-L(.*)$" "\\1" _dir "${_param}")
                        IF(IS_DIRECTORY ${_dir})
                            LIST(APPEND ${_NAMEVAR}_LIBRARY_PATH ${_dir})
                            LINK_DIRECTORIES(${_dir})
                        ELSE()
                            MESSAGE(FATAL_ERROR "Looking for ${_NAMEVAR} - ${_dir} is not a directory")
                        ENDIF()
                    ELSEIF(_param MATCHES "^-l")
                        STRING(REGEX REPLACE "^-l(.*)$" "\\1" _lib "${_param}")
                        LIST(APPEND ${_NAMEVAR}_LIBRARIES ${_lib})
                    ELSEIF(_param MATCHES "^-I")
                        STRING(REGEX REPLACE "^-I(.*)$" "\\1" _dir "${_param}")
                        IF(IS_DIRECTORY ${_dir})
                            LIST(APPEND ${_NAMEVAR}_INCLUDE_PATH ${_dir})
                            INCLUDE_DIRECTORIES(${_dir})
                        ELSE()
                            MESSAGE(FATAL_ERROR "Looking for ${_NAMEVAR} - ${_dir} is not a directory")
                        ENDIF()
                    ELSE()
                        MESSAGE(FATAL_ERROR "Looking for ${_NAMEVAR} - ${_param} unrecognized")
                    ENDIF()    
                ENDFOREACH()
                IF(DEFINED ${_NAMEVAR}_INC)
                    IF(IS_DIRECTORY ${${_NAMEVAR}_INC})
                        LIST(APPEND ${_NAMEVAR}_INCLUDE_PATH ${${_NAMEVAR}_INC})
                        INCLUDE_DIRECTORIES(${${_NAMEVAR}_INCLUDE_PATH})
                    ELSE()
                        MESSAGE(FATAL_ERROR "Looking for ${_NAMEVAR} - ${${_NAMEVAR}_INC} is not a directory")
                    ENDIF()
                ENDIF(DEFINED ${_NAMEVAR}_INC)

                # Test library
                INCLUDE(findPACKAGE)
                INFO_PACKAGE()
                IF(DEFINED ${_NAMEVAR}_INCLUDE_PATH)
                    CHECK_MY_PACKAGE(${_NAMEVAR} TRUE TRUE TRUE)
                ELSE()
                    CHECK_MY_PACKAGE(${_NAMEVAR} FALSE TRUE TRUE)
                ENDIF()
                IF(${_NAMEVAR}_C_FOUND OR ${_NAMEVAR}_F_FOUND)
                    SET(${_NAMEVAR}_FOUND TRUE)
                    SET(HAVE_${_NAMEVAR} ON)
                    MESSAGE(STATUS "Looking for ${_NAMEVAR} - found")
                ELSE()
                    MESSAGE(FATAL_ERROR "Looking for ${_NAMEVAR} - not working")
                ENDIF()
            ELSE()
                MESSAGE(FATAL_ERROR "Looking for ${_NAMEVAR} - not found")
            ENDIF()

        ELSEIF(${_NAMEVAR}_USE_SYSTEM)
            #
            # TODO: SET REQUIRED and search only in system path only
            #
            #IF(${${_NAMEVAR}_REQUIRED})
            #    FIND_PACKAGE(${_NAMEVAR} REQUIRED)
            #ELSE(${${_NAMEVAR}_REQUIRED})
            #    FIND_PACKAGE(${_NAMEVAR} QUIET)
            #ENDIF(${${_NAMEVAR}_REQUIRED})
            FIND_PACKAGE(${_NAMEVAR})
            IF(${_NAMEVAR}_FOUND)
                SET(HAVE_${_NAMEVAR} ON)
                LINK_DIRECTORIES(${${_NAMEVAR}_LIBRARY_PATH})
                INCLUDE_DIRECTORIES(${${_NAMEVAR}_INCLUDE_PATH})
                LIST(INSERT EXTRA_LIBS 0 ${${_NAMEVAR}_LIBRARIES})
            ELSE(${_NAMEVAR}_FOUND)
                MESSAGE(FATAL_ERROR "Looking for ${_NAMEVAR} - not found")
            ENDIF()

        ELSEIF(${_NAMEVAR}_USE_TARBALL)
            MORSE_INSTALL_PACAKGE(${_NAMEVAR} "TARBALL")
            SET(HAVE_${_NAMEVAR} ON)
            LINK_DIRECTORIES(${${_NAMEVAR}_LIBRARY_PATH})
            INCLUDE_DIRECTORIES(${${_NAMEVAR}_INCLUDE_PATH})
            LIST(INSERT EXTRA_LIBS 0 ${${_NAMEVAR}_LIBRARIES})

        ELSEIF(${_NAMEVAR}_USE_WEB)
            MORSE_INSTALL_PACAKGE(${_NAMEVAR} "WEB")
            SET(HAVE_${_NAMEVAR} ON)
            LINK_DIRECTORIES(${${_NAMEVAR}_LIBRARY_PATH})
            INCLUDE_DIRECTORIES(${${_NAMEVAR}_INCLUDE_PATH})
            LIST(INSERT EXTRA_LIBS 0 ${${_NAMEVAR}_LIBRARIES})

        ELSEIF(${_NAMEVAR}_USE_SVN)
            MORSE_INSTALL_PACAKGE(${_NAMEVAR} "REPO")
            SET(HAVE_${_NAMEVAR} ON)
            LINK_DIRECTORIES(${${_NAMEVAR}_LIBRARY_PATH})
            INCLUDE_DIRECTORIES(${${_NAMEVAR}_INCLUDE_PATH})
            LIST(INSERT EXTRA_LIBS 0 ${${_NAMEVAR}_LIBRARIES})

        ELSEIF(${_NAMEVAR}_USE_AUTO)
            #
            # TODO: SET REQUIRED and search only in system path only
            #
            #IF(${${_NAMEVAR}_REQUIRED})
            #    FIND_PACKAGE(${_NAMEVAR} REQUIRED)
            #ELSE(${${_NAMEVAR}_REQUIRED})
            #    FIND_PACKAGE(${_NAMEVAR} QUIET)
            #ENDIF(${${_NAMEVAR}_REQUIRED})
            FIND_PACKAGE(${_NAMEVAR})
            IF(${_NAMEVAR}_FOUND)
                SET(HAVE_${_NAMEVAR} ON)
                LINK_DIRECTORIES(${${_NAMEVAR}_LIBRARY_PATH})
                INCLUDE_DIRECTORIES(${${_NAMEVAR}_INCLUDE_PATH})
                LIST(INSERT EXTRA_LIBS 0 ${${_NAMEVAR}_LIBRARIES})

            ELSE(${_NAMEVAR}_FOUND)
                IF(${${_NAMEVAR}_REQUIRED})
                    MORSE_INSTALL_PACAKGE(${_NAMEVAR} "AUTO")
                    SET(HAVE_${_NAMEVAR} ON)
                    LINK_DIRECTORIES(${${_NAMEVAR}_LIBRARY_PATH})
                    INCLUDE_DIRECTORIES(${${_NAMEVAR}_INCLUDE_PATH})
                    LIST(INSERT EXTRA_LIBS 0 ${${_NAMEVAR}_LIBRARIES})

                ELSE(${_NAMEVAR}_REQUIRED)
                    MESSAGE(STATUS "Looking for ${_NAMEVAR} - not found")
                    MESSAGE(STATUS "Installing ${_NAMEVAR} - not necessary")

                ENDIF(${${_NAMEVAR}_REQUIRED})
            ENDIF()

        ELSE()
            MESSAGE(STATUS "Looking for ${_NAMEVAR} - not found")
            MESSAGE(STATUS "Installing ${_NAMEVAR} - no")

        ENDIF()

    ELSE(${${_NAMEVAR}_TRYUSE})
        MESSAGE(STATUS "Looking for ${_NAMEVAR} - no")
        MESSAGE(STATUS "Installing ${_NAMEVAR} - no")

    ENDIF(${${_NAMEVAR}_TRYUSE})

ENDMACRO(DEFINE_PACKAGE)

###
### END definePACKAGE.cmake
###
