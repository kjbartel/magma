###     
#
#  @file Summmary.cmake
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

MACRO(SUMMARY)
########################
#
# MAGMA Summary
#
########################
    MESSAGE(STATUS "-----------------------------------------------------------------------------")
    MESSAGE(STATUS "Status of your MAGMA project configuration:")
    MESSAGE(STATUS "-----------------------------------------------------------------------------")
    IF(MAGMA)
        MESSAGE(STATUS "")
        MESSAGE(STATUS "Configuration of the MAGMA project:")
        MESSAGE(STATUS "  - MAGMA_USE_FERMI         : ${MAGMA_USE_FERMI}")
        MESSAGE(STATUS "  - MAGMA_USE_PLASMA        : ${MAGMA_USE_PLASMA}")
        MESSAGE(STATUS "  - MAGMA_ENABLE_TESTING : ${MAGMA_ENABLE_TESTING}")
    ENDIF(MAGMA)

########################
#
# IMPLEMENTATION Summary
#
########################

    MESSAGE(STATUS "")
    MESSAGE(STATUS "Configuration of the external libraries: implementation")
    FOREACH(_lib CUDA OPENCL)
        IF(HAVE_${_lib})
            MESSAGE(STATUS "  - ${_lib} :")
            MESSAGE(STATUS "    --> ${_lib}_MODE         : ${${_lib}_USED_MODE}")
            IF(NOT "^${_lib}$" STREQUAL "^MPI$")
            MESSAGE(STATUS "    --> ${_lib}_LIBRARY_PATH : ${${_lib}_LIBRARY_PATH}")
            ENDIF()
            MESSAGE(STATUS "    --> ${_lib}_INCLUDE_PATH : ${${_lib}_INCLUDE_PATH}")
            MESSAGE(STATUS "    --> ${_lib}_LIBRARIES    : ${${_lib}_LIBRARIES}") 
            IF(MAGMA_DEBUG_CMAKE)
            MESSAGE(STATUS "    --> ${_lib}_LIBRARY      : ${${_lib}_LIBRARY}") 
            ENDIF(MAGMA_DEBUG_CMAKE)
            MESSAGE(STATUS "")
        ENDIF(HAVE_${_lib})
    ENDFOREACH()

########################
#
# LINEAR ALGEBRA Summary
#
########################

    MESSAGE(STATUS "Configuration of the external libraries: linear algebra")
    FOREACH(_lib BLAS LAPACK TMG CBLAS)
        IF(HAVE_${_lib})
            MESSAGE(STATUS "  - ${_lib} :")
            MESSAGE(STATUS "    --> ${_lib}_MODE         : ${${_lib}_USED_MODE}")
            IF("^${_lib}$" STREQUAL "^TMG$"   OR
               "^${_lib}$" STREQUAL "^BLAS$"  OR
               "^${_lib}$" STREQUAL "^CBLAS$" OR
               "^${_lib}$" STREQUAL "^LAPACK$"  )
            MESSAGE(STATUS "    --> ${_lib}_VENDOR       : ${${_lib}_VENDOR}")
            ENDIF()
            MESSAGE(STATUS "    --> ${_lib}_LIBRARY_PATH : ${${_lib}_LIBRARY_PATH}")
            IF(NOT "^${_lib}$" STREQUAL "^BLAS$" OR
               NOT "^${_lib}$" STREQUAL "^LAPACK$")
            MESSAGE(STATUS "    --> ${_lib}_INCLUDE_PATH : ${${_lib}_INCLUDE_PATH}")
            ENDIF()
            MESSAGE(STATUS "    --> ${_lib}_LIBRARIES    : ${${_lib}_LIBRARIES}")
            IF(MAGMA_DEBUG_CMAKE)
            MESSAGE(STATUS "    --> ${_lib}_LIBRARY      : ${${_lib}_LIBRARY}")
            ENDIF(MAGMA_DEBUG_CMAKE)
            MESSAGE(STATUS "")
        ENDIF(HAVE_${_lib})
    ENDFOREACH()

#######################
#
# Help 
#
#######################
    MESSAGE(STATUS "")
    MESSAGE(STATUS "--------------+--------------------------------------------------------------")
    MESSAGE(STATUS "Command       | Description")
    MESSAGE(STATUS "--------------+--------------------------------------------------------------")
    MESSAGE(STATUS "make          | Compile MAGMA project in:")
    MESSAGE(STATUS "              |     ${CMAKE_BINARY_DIR}")
    MESSAGE(STATUS "              |")
    MESSAGE(STATUS "make test     | Build and run the unit-tests.")
    MESSAGE(STATUS "              |")
    MESSAGE(STATUS "make install  | Install MAGMA projecs in:")
    MESSAGE(STATUS "              |     ${CMAKE_INSTALL_PREFIX}")
    MESSAGE(STATUS "              | To change that:")
    MESSAGE(STATUS "              |     cmake ../ -DCMAKE_INSTALL_PREFIX=yourpath")
    MESSAGE(STATUS "              |")
    MESSAGE(STATUS "--------------+--------------------------------------------------------------")

ENDMACRO(SUMMARY)

##
## @end file Summary.cmake
##
