###
#
#  @file populatePACKAGE.cmake
#
#  @project MAGMA
#  MAGMA is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver. 
# 
#  @version 1.0.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @date 13-07-2012
#   
###

###
#
# POPULATE_COMPILE_SYSTEM: Call INCLUDE_DIRECTORIES, LINK_DIRECTORIES, etc...
#
###
MACRO(POPULATE_COMPILE_SYSTEM _PCS_PACKAGE)

    # Define input variable to uppercase
    # ----------------------------------
    STRING(TOUPPER "${_PCS_PACKAGE}" _UPPERPCS_PACKAGE)
    STRING(TOLOWER "${_PCS_PACKAGE}" _LOWERPCS_PACKAGE)
    IF("^${_UPPERPCS_PACKAGE}$" STREQUAL "^CUDA$")
        MESSAGE(STATUS "Looking for CUDA - found")

        # Define all flags needed
        IF(CUDA_VERSION VERSION_LESS "4.0")
            SET(CUDA_HOST_COMPILATION_CPP OFF)
        ENDIF(CUDA_VERSION VERSION_LESS "4.0")
        SET(CUDA_BUILD_EMULATION OFF)
        IF(${ARCH_X86_64})
            LINK_DIRECTORIES(${CUDA_TOOLKIT_ROOT_DIR}/lib64)
            SET(CUDA_LIBRARY_PATH "${CUDA_TOOLKIT_ROOT_DIR}/lib64")
        ELSE()
            LINK_DIRECTORIES(${CUDA_TOOLKIT_ROOT_DIR}/lib)
            SET(CUDA_LIBRARY_PATH "${CUDA_TOOLKIT_ROOT_DIR}/lib")
        ENDIF()
        SET(CUDA_LIBRARIES "cudart;cublas;cuda")
        SET(CUDA_LDFLAGS "-L${CUDA_LIBRARY_PATH} -lcudart -lcuda")
        SET(CUDA_LIBRARY "${CUDA_CUDART_LIBRARY};${CUDA_CUDA_LIBRARY}")
        SET(CUDA_INCLUDE_PATH "${CUDA_INCLUDE_DIRS}")
        INCLUDE_DIRECTORIES(${CUDA_INCLUDE_PATH})
        IF("^${MAGMA_USE_CUDA}$" STREQUAL "^$")
            SET(MAGMA_USE_CUDA "ON" CACHE STRING
                "Enable/Disable CUDA dependency (ON/OFF/<not-defined>)" FORCE)
        ENDIF()
        SET(CUDA_USED_MODE "FIND")
        SET(HAVE_CUDA ON)
        SET(HAVE_CUBLAS ON)
        
        # Should move them to a magma_config.h.in 
        ADD_DEFINITIONS(-DHAVE_CUDA)
        ADD_DEFINITIONS(-DHAVE_CUBLAS)

        # Fortran mangling for cuda >= 5.0
        SET(MORSE_CUDA_WITH_FORTRAN ON)
        IF(CUDA_VERSION VERSION_GREATER "4.9")
            IF(FORTRAN_MANGLING_DETECTED)
                IF("${FORTRAN_MANGLING_DETECTED}" STREQUAL "-DADD_")
                    ADD_DEFINITIONS(-DCUBLAS_GFORTRAN)
                ELSEIF("${FORTRAN_MANGLING_DETECTED}" STREQUAL "-DUPCASE")
                    ADD_DEFINITIONS(-DCUBLAS_INTEL_FORTRAN)
                ELSE()
                    SET(MORSE_CUDA_WITH_FORTRAN OFF)
                ENDIF()
            ENDIF()
        ENDIF()

    ELSE()
        LINK_DIRECTORIES(${${_UPPERPCS_PACKAGE}_LIBRARY_PATH})
        INCLUDE_DIRECTORIES(${${_UPPERPCS_PACKAGE}_INCLUDE_PATH})
        IF("^${MORSE_USE_${_UPPERPCS_PACKAGE}}$" STREQUAL "^$")
            SET(MORSE_USE_${_UPPERPCS_PACKAGE} "ON" CACHE STRING
                "Enable/Disable ${_UPPERPCS_PACKAGE} dependency (ON/OFF/<not-defined>)" FORCE)
        ENDIF()
        SET(HAVE_${_UPPERPCS_PACKAGE} ON)

    ENDIF()

    # Fill the following variables
    #   - ${_UPPERPCS_PACKAGE}_LDFLAGS
    #   - ${_UPPERPCS_PACKAGE}_LIST_LDFLAGS
    # -------------------------------------
    #UNSET(${_UPPERPCS_PACKAGE}_LDFLAGS)
    UNSET(${_UPPERPCS_PACKAGE}_LIST_LDFLAGS)
    FOREACH(_path ${${_UPPERPCS_PACKAGE}_LIBRARY_PATH})
        LIST(APPEND ${_UPPERPCS_PACKAGE}_LIST_LDFLAGS "-L${_path}")
    ENDFOREACH()
    FOREACH(_lib ${${_UPPERPCS_PACKAGE}_LIBRARIES})
        LIST(APPEND ${_UPPERPCS_PACKAGE}_LIST_LDFLAGS "-l${_lib}")
    ENDFOREACH()

ENDMACRO(POPULATE_COMPILE_SYSTEM)

##
## @end file populatePACKAGE.cmake 
##
