###
#
#  @file checkPACKAGE.cmake
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

MACRO(CHECK_PACKAGE _NAME)

    # Import cmake modules
    # --------------------
    INCLUDE(CheckIncludeFiles)
    INCLUDE(CheckFunctionExists)
    INCLUDE(CheckFortranFunctionExists)

    # Set internal variable
    # ---------------------
    STRING(TOUPPER "${_NAME}" _NAMEVAR)

    # Set EXTRA_DEPENDENCIES
    # ----------------------
    SET(_EXTRA_DEPENDENCIES ${${_NAMEVAR}_FIND_DEPS})

    # Load infoPACKAGE
    # ----------------
    #INCLUDE(infoPACKAGE)
    #INFO_FIND_PACKAGE(${_NAMEVAR})

    # Set ${_NAMEVAR}_ERROR_OCCURRED to TRUE if an error is encountered
    # -----------------------------------------------------------------
    IF(NOT DEFINED ${_NAMEVAR}_ERROR_OCCURRED)
        SET(${_NAMEVAR}_ERROR_OCCURRED FALSE)
    ENDIF(NOT DEFINED ${_NAMEVAR}_ERROR_OCCURRED)

    # Check if headers are provided
    # -----------------------------
    IF(${_NAMEVAR}_name_include)

        # Add include path
        SET(${_NAMEVAR}_TEST_INCS ${${_NAMEVAR}_INCLUDE_PATH})
        IF(${_NAMEVAR}_name_include_suffix)
            LIST(APPEND "${${_NAMEVAR}_INCLUDE_PATH}/${${_NAMEVAR}_name_include_suffix}")
        ENDIF(${_NAMEVAR}_name_include_suffix)
        FOREACH(_package ${_EXTRA_DEPENDENCIES})
            IF(NOT "^${${_package}_INCLUDE_PATH}$" STREQUAL "^$")
                LIST(APPEND ${_NAMEVAR}_TEST_INCS "${${_package}_INCLUDE_PATH}")
            ENDIF()
        ENDFOREACH()
        SET(CMAKE_REQUIRED_INCLUDES ${${_NAMEVAR}_TEST_INCS})

        # Add linker flags
        SET(${_NAMEVAR}_TEST_FLAGS ${${_NAMEVAR}_LDFLAGS})
        FOREACH(_package ${_EXTRA_DEPENDENCIES})
            IF(NOT "^${${_package}_LDFLAGS}$" STREQUAL "^$")
                SET(${_NAMEVAR}_TEST_FLAGS "${${_NAMEVAR}_TEST_FLAGS} ${${_package}_LDFLAGS}")
            ENDIF()
        ENDFOREACH()
        SET(CMAKE_REQUIRED_FLAGS ${${_NAMEVAR}_TEST_FLAGS})

        # Launch tests
        STRING(TOUPPER "${_header}" _FILE)
        STRING(REPLACE "." "_" _FILE "${_FILE}")
        CHECK_INCLUDE_FILES("${${_NAMEVAR}_name_include}" HAVE_ALL_${_NAMEVAR}_HEADER)
        IF(NOT HAVE_ALL_${_NAMEVAR}_HEADER)
            SET(${_NAMEVAR}_ERROR_OCCURRED TRUE)
        ENDIF(NOT HAVE_ALL_${_NAMEVAR}_HEADER)

    ENDIF(${_NAMEVAR}_name_include)

    # Check if library worked
    # -----------------------
    IF(${_NAMEVAR}_name_library)
        # Check C funtion
        IF("${${_NAMEVAR}_type_library}" MATCHES "C")
            FOREACH(_fct ${${_NAMEVAR}_name_fct_test})
                # Define a hash for return value
                STRING(REPLACE ";" "_" hres_prefix "${_EXTRA_DEPENDENCIES}")

                # Add include path
                SET(${_NAMEVAR}_TEST_INCS ${${_NAMEVAR}_INCLUDE_PATH})
                FOREACH(_package ${_EXTRA_DEPENDENCIES})
                    IF(NOT "^${${_package}_INCLUDE_PATH}$" STREQUAL "^$")
                        LIST(APPEND ${_NAMEVAR}_TEST_INCS "${${_package}_INCLUDE_PATH}")
                    ENDIF()
                ENDFOREACH()
                SET(CMAKE_REQUIRED_INCLUDES ${${_NAMEVAR}_TEST_INCS})

                # Add linker flags
                SET(${_NAMEVAR}_TEST_FLAGS ${${_NAMEVAR}_LDFLAGS})
                FOREACH(_package ${_EXTRA_DEPENDENCIES})
                    IF(NOT "^${${_package}_LDFLAGS}$" STREQUAL "^$")
                        SET(${_NAMEVAR}_TEST_FLAGS "${${_NAMEVAR}_TEST_FLAGS} ${${_package}_LDFLAGS}")
                    ENDIF()
                ENDFOREACH()
                SET(CMAKE_REQUIRED_FLAGS ${${_NAMEVAR}_TEST_FLAGS})

                # Add library
                SET(${_NAMEVAR}_TEST_LIBS ${${_NAMEVAR}_LIBRARY})
                FOREACH(_package ${_EXTRA_DEPENDENCIES};GFORTRAN;IFCORE)
                    IF(NOT "^${${_package}_LIBRARY}$" STREQUAL "^$")
                        LIST(APPEND ${_NAMEVAR}_TEST_LIBS "${${_package}_LIBRARY}")
                    ENDIF()
                ENDFOREACH()
                SET(CMAKE_REQUIRED_LIBRARIES ${${_NAMEVAR}_TEST_LIBS})

                # Check library
                CHECK_FUNCTION_EXISTS(${_fct} ${hres_prefix}_${_fct}_C_FOUND)

                # Debug message
                IF(MAGMA_DEBUG_CMAKE)
                    MESSAGE(STATUS "  * debug:")
                    MESSAGE(STATUS "  * debug: CHECK ${_PACKAGE}:")
                    MESSAGE(STATUS "  * debug:")
                    MESSAGE(STATUS "  * debug: --> extra libraries:")
                    MESSAGE(STATUS "  * debug:     - _EXTRA_DEPENDENCIES    : ${_EXTRA_DEPENDENCIES}")
                    MESSAGE(STATUS "  * debug: --> used parameters:")
                    MESSAGE(STATUS "  * debug:     - ${_NAMEVAR}_TEST_INCS  : ${${_NAMEVAR}_TEST_INCS}")
                    MESSAGE(STATUS "  * debug:     - ${_NAMEVAR}_TEST_FLAGS : ${${_NAMEVAR}_TEST_FLAGS}")
                    MESSAGE(STATUS "  * debug:     - ${_NAMEVAR}_TEST_LIBS  : ${${_NAMEVAR}_TEST_LIBS}")
                    MESSAGE(STATUS "  * debug: --> return code:")
                    MESSAGE(STATUS "  * debug:     - ${_fct}_C_FOUND        : ${${hres_prefix}_${_fct}_C_FOUND}")
                    MESSAGE(STATUS "  * debug:")
                ENDIF(MAGMA_DEBUG_CMAKE)

                # Handle callback
                IF(${${hres_prefix}_${_fct}_C_FOUND})
                    SET(${_NAMEVAR}_C_FOUND TRUE)
                ELSE(${${hres_prefix}_${_fct}_C_FOUND})
                    SET(${_NAMEVAR}_C_FOUND FALSE)
                    BREAK()
                ENDIF(${${hres_prefix}_${_fct}_C_FOUND})
            ENDFOREACH()
        ENDIF()

        # Check Fortran funtion
        IF("${${_NAMEVAR}_type_library}" MATCHES "Fortran")
            FOREACH(_fct ${${_NAMEVAR}_name_fct_test})
                # Define a hash for return value
                STRING(REPLACE ";" "_" hres_prefix "${_EXTRA_DEPENDENCIES}")

                # Add library
                SET(${_NAMEVAR}_TEST_LIBS ${${_NAMEVAR}_LIBRARY})
                FOREACH(_package ${_EXTRA_DEPENDENCIES};GFORTRAN;IFCORE)
                    IF(NOT "^${${_package}_LIBRARY}$" STREQUAL "^$")
                        LIST(APPEND ${_NAMEVAR}_TEST_LIBS "${${_package}_LIBRARY}")
                    ENDIF()
                ENDFOREACH()
                SET(CMAKE_REQUIRED_LIBRARIES ${${_NAMEVAR}_TEST_LIBS})

                # Check library
                CHECK_FORTRAN_FUNCTION_EXISTS(${_fct} ${hres_prefix}_${_fct}_F_FOUND)

                # Debug message
                IF(MAGMA_DEBUG_CMAKE)
                    MESSAGE(STATUS "  * debug:")
                    MESSAGE(STATUS "  * debug: CHECK ${_PACKAGE}:")
                    MESSAGE(STATUS "  * debug:")
                    MESSAGE(STATUS "  * debug: --> extra libraries:")
                    MESSAGE(STATUS "  * debug:     - _EXTRA_DEPENDENCIES    : ${_EXTRA_DEPENDENCIES}")
                    MESSAGE(STATUS "  * debug: --> used parameters:")
                    MESSAGE(STATUS "  * debug:     - ${_NAMEVAR}_TEST_LIBS  : ${${_NAMEVAR}_TEST_LIBS}")
                    MESSAGE(STATUS "  * debug: --> return code:")
                    MESSAGE(STATUS "  * debug:     - ${_fct}_F_FOUND        : ${${hres_prefix}_${_fct}_F_FOUND}")
                    MESSAGE(STATUS "  * debug:")
                ENDIF(MAGMA_DEBUG_CMAKE)

                # Handle callback
                IF(${${hres_prefix}_${_fct}_F_FOUND})
                    SET(${_NAMEVAR}_F_FOUND TRUE)
                ELSE(${${hres_prefix}_${_fct}_F_FOUND})
                    SET(${_NAMEVAR}_F_FOUND FALSE)
                    BREAK()
                ENDIF(${${hres_prefix}_${_fct}_F_FOUND})
            ENDFOREACH()
        ENDIF()

        # Provide an error if necessary
        IF(${_NAMEVAR}_C_FOUND OR ${_NAMEVAR}_F_FOUND)
            MESSAGE(STATUS "Looking for ${_NAMEVAR} - working")
        ELSE()
            SET(${_NAMEVAR}_ERROR_OCCURRED TRUE)
        ENDIF()

    ENDIF(${_NAMEVAR}_name_library)

ENDMACRO(CHECK_PACKAGE)

##
## @end file checkPACKAGE.cmake
##
