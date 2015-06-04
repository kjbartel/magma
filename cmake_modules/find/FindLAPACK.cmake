###
#
#  @file FindLAPACK.cmake
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
# This module sets the following variables:
#  LAPACK_VENDOR       : if set checks only the specified vendor, if not set checks all the possibilities
#  LAPACK_FOUND        : set to true if a library implementing the LAPACK interface is found
#
#  LAPACK_LDFLAGS      : uncached string of all required linker flags.
#  LAPACK_LIBRARY      : uncached list of library
#  LAPACK_LIBRARIES    : uncached list of required linker flags (with -l and -L).
#  LAPACK_LIBRARY_PATH : uncached path of the library directory of LAPACK installation.
#
###
#
# You can choose your prefered LAPACK with this option : -DMAGMA_USE_LAPACK="VENDOR"
#    NB: "VENDOR" can be in uppercase or lowercase
#
# You can defined LAPACK_DIR to help the search of LAPACK
#
###


# Early exit if already defined
# -----------------------------
IF(DEFINED LAPACK_USED_MODE)
    MESSAGE(STATUS "Looking for LAPACK - mode ${${_DEPVAR}_USED_MODE} already defined")
    RETURN()
ENDIF(DEFINED LAPACK_USED_MODE)

# Early exit if already searched
# ------------------------------
IF(FINDLAPACK_TESTED)
    IF(LAPACK_FOUND)
        MESSAGE(STATUS "Looking for LAPACK - already found")
    ELSE(LAPACK_FOUND)
        MESSAGE(STATUS "Looking for LAPACK - already NOT found")
    ENDIF(LAPACK_FOUND)
    RETURN()
ENDIF(FINDLAPACK_TESTED)
SET(FINDLAPACK_TESTED TRUE)

# Message
# -------
MESSAGE(STATUS "Looking for LAPACK")

# Include
# -------
INCLUDE(CheckLibraryExists)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckFortranTypeSizes)
INCLUDE(CheckFortranFunctionExists)

# Check the language being used
# -----------------------------
GET_PROPERTY(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
IF(_LANGUAGES_ MATCHES Fortran)
    SET(_CHECK_FORTRAN TRUE)

ELSEIF(_LANGUAGES_ MATCHES C OR _LANGUAGES_ MATCHES CXX)
    SET(_CHECK_FORTRAN FALSE)

ELSE()
    IF(LAPACK_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Looking for LAPACK - requires Fortran, C, or C++ to be enabled.")
    ELSE(LAPACK_FIND_REQUIRED)
        MESSAGE(STATUS "Looking for LAPACK - not found (unsupported languages)")
        RETURN()
    ENDIF(LAPACK_FIND_REQUIRED)

ENDIF()

# Add some path to help
# ---------------------
IF(WIN32)
    STRING(REPLACE ":" ";" _lib_env "$ENV{LIB}")
ELSE(WIN32)
    IF(APPLE)
        STRING(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
    ELSE()
        STRING(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
    ENDIF()
    LIST(APPEND _lib_env "/usr/local/lib64")
    LIST(APPEND _lib_env "/usr/lib64")
    LIST(APPEND _lib_env "/usr/local/lib")
    LIST(APPEND _lib_env "/usr/lib")
ENDIF()

UNSET(LAPACK_EXTRA_PATH)
IF(NOT "^${BLAS_DIR}$" STREQUAL "^$")
    LIST(APPEND LAPACK_EXTRA_PATH ${BLAS_DIR})
ENDIF()
IF(NOT "^${LAPACK_DIR}$" STREQUAL "^$")
    LIST(APPEND LAPACK_EXTRA_PATH ${LAPACK_DIR})
ENDIF()
IF(NOT "^${_lib_env}$" STREQUAL "^$")
    LIST(APPEND LAPACK_EXTRA_PATH ${_lib_env})
ENDIF()

# Determine the default integer size
# ----------------------------------
CHECK_FORTRAN_TYPE_SIZES()
IF(NOT SIZEOF_INTEGER)
    MESSAGE(WARNING "Looking for LAPACK - Unable to determine default integer size")
    MESSAGE(WARNING "Looking for LAPACK - Assuming integer*4")
    SET(SIZEOF_INTEGER 4)

ENDIF(NOT SIZEOF_INTEGER)

# Macro to locate a library and check for a specified symbol
# ----------------------------------------------------------
MACRO(LAPACK_LOCATE_AND_TEST _prefix _lapack_libnames _flags _blas_ldflags)

    # Status of internal variables
    IF(MAGMA_DEBUG_CMAKE)
        MESSAGE( STATUS "Looking for LAPACK - try ${_prefix}")
        MESSAGE(STATUS "  * debug:")
        MESSAGE(STATUS "  * debug: Looking for LAPACK - status of ${_prefix}")
        MESSAGE(STATUS "  * debug:   - libraries             : ${_lapack_libnames}")
        MESSAGE(STATUS "  * debug:   - additional flags      : ${_flags}")
    ENDIF(MAGMA_DEBUG_CMAKE)

    # Initialize values
    SET(LAPACK_${_prefix}_FOUND        FALSE)
    SET(LAPACK_${_prefix}_LDFLAGS      ""   )
    SET(LAPACK_${_prefix}_LIBRARY      ""   )
    SET(LAPACK_${_prefix}_LIBRARIES    ""   )
    SET(LAPACK_${_prefix}_LIBRARY_PATH ""   )
    MARK_AS_ADVANCED(LAPACK_${_prefix}_LDFLAGS)
    MARK_AS_ADVANCED(LAPACK_${_prefix}_LIBRARY)
    MARK_AS_ADVANCED(LAPACK_${_prefix}_LIBRARIES)
    MARK_AS_ADVANCED(LAPACK_${_prefix}_LIBRARY_PATH)

    # Set the library suffix to look for
    IF(WIN32)
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".dll;.lib")
    ELSE(WIN32)
        IF(APPLE)
            SET(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib;.a")
        ELSE(APPLE)
            # for ubuntu's libblas3gf and liblapack3gf packages
            SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a;.so.3gf")
        ENDIF(APPLE)
    ENDIF(WIN32)

    # Set paths where we are looking for
    UNSET(LAPACK_SEARCH_PATH)
    IF(NOT "^${LAPCK_${_prefix}_DIR}$" STREQUAL "^$")
        LIST(APPEND LAPACK_SEARCH_PATH ${LAPACK_${_prefix}_DIR})
    ENDIF()
    IF(NOT "^${LAPACK_EXTRA_PATH}$" STREQUAL "^$")
        LIST(APPEND LAPACK_SEARCH_PATH ${LAPACK_EXTRA_PATH})
    ENDIF()
    IF(MAGMA_DEBUG_CMAKE)
        MESSAGE(STATUS "  * debug:   - searching in this paths : ${LAPACK_SEARCH_PATH}")
    ENDIF(MAGMA_DEBUG_CMAKE)

    # Define path we have to looking for first
    SET(CMAKE_PREFIX_PATH ${LAPACK_SEARCH_PATH})

    # Looking for path of libraries
    SET({LAPACK_${_prefix}_LDFLAGS ${_lapack_libnames})
    FOREACH(_library ${_lapack_libnames})
        FIND_LIBRARY(LAPACK_${_prefix}_${_library}_LIBRARY
                     NAMES ${_library}
                     PATHS ${LAPACK_SEARCH_PATH}
                    )
        MARK_AS_ADVANCED(LAPACK_${_prefix}_${_library}_LIBRARY)

        IF(MAGMA_DEBUG_CMAKE)
            MESSAGE(STATUS "  * debug:   - obtained library path : ${LAPACK_${_prefix}_${_library}_LIBRARY}")
        ENDIF(MAGMA_DEBUG_CMAKE)

        IF(LAPACK_${_prefix}_${_library}_LIBRARY)
            LIST(APPEND LAPACK_${_prefix}_LIBRARY   "${LAPACK_${_prefix}_${_library}_LIBRARY}")
            LIST(APPEND LAPACK_${_prefix}_LIBRARIES "${_library}"                           )
            GET_FILENAME_COMPONENT(LAPACK_${_prefix}_${_library}_LIBRARY_FILENAME
                                   ${LAPACK_${_prefix}_${_library}_LIBRARY}
                                   NAME)
            STRING(REGEX REPLACE "(.*)/${LAPACK_${_prefix}_${_library}_LIBRARY_FILENAME}" "\\1"
                   LAPACK_${_prefix}_${_library}_LIBRARY_PATH "${LAPACK_${_prefix}_${_library}_LIBRARY}")
            LIST(APPEND LAPACK_${_prefix}_LIBRARY_PATH   "${LAPACK_${_prefix}_${_library}_LIBRARY_PATH}")
            LIST(REMOVE_DUPLICATES LAPACK_${_prefix}_LIBRARY_PATH)

        ELSE()
            SET(LAPACK_${_prefix}_LIBRARY      "LAPACK_${_prefix}_LIBRARY-NOTFOUND"     )
            SET(LAPACK_${_prefix}_LIBRARIES    "LAPACK_${_prefix}_LIBRARIES-NOTFOUND"   )
            SET(LAPACK_${_prefix}_LIBRARY_PATH "LAPACK_${_prefix}_LIBRARY_PATH-NOTFOUND")
            BREAK()
        ENDIF()

    ENDFOREACH()

    # Define path we have to looking for first
    UNSET(CMAKE_PREFIX_PATH)

    # Test this combination of libraries.
    IF(LAPACK_${_prefix}_LIBRARIES)
        # Set CMAKE_REQUIRED_LIBRARIES to test the function(s)
        UNSET(CMAKE_REQUIRED_LIBRARIES)
        SET(LAPACK_${_prefix}_EXTRA_FLAGS ${_flags})
        SET(LAPACK_${_prefix}_LDFLAGS      ""      )
        FOREACH(_path ${LAPACK_${_prefix}_LIBRARY_PATH})
            SET(LAPACK_${_prefix}_LDFLAGS "${LAPACK_${_prefix}_LDFLAGS} -L${_path}")
        ENDFOREACH()
        FOREACH(_flag ${LAPACK_${_prefix}_LIBRARIES})
            SET(LAPACK_${_prefix}_LDFLAGS "${LAPACK_${_prefix}_LDFLAGS} -l${_flag}")
        ENDFOREACH()
        FOREACH(_flag ${LAPACK_${_prefix}_EXTRA_FLAGS})
            SET(LAPACK_${_prefix}_LDFLAGS "${LAPACK_${_prefix}_LDFLAGS} ${_flag}")
        ENDFOREACH()
        STRING(REGEX REPLACE "^ " "" LAPACK_${_prefix}_LDFLAGS "${LAPACK_${_prefix}_LDFLAGS}")
        SET(CMAKE_REQUIRED_LIBRARIES "${LAPACK_${_prefix}_LDFLAGS} ${_blas_ldflags}")
        STRING(REPLACE " " ";" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

        # Load infoPACKAGE
        INCLUDE(infoLAPACK)
        LAPACK_INFO_FIND()

        # Test [X]heev
        FOREACH(_fct ${LAPACK_name_fct_test})
            IF(_CHECK_FORTRAN)
                CHECK_FORTRAN_FUNCTION_EXISTS("${_fct}" LAPACK_${_prefix}_${_fct}_WORKS)

            ELSE(_CHECK_FORTRAN)
                IF(${FORTRAN_MANGLING_DETECTED} MATCHES "-DADD_")
                    SET(_fct_C "${_fct}_")

                ELSEIF(${FORTRAN_MANGLING_DETECTED} MATCHES "-DNOCHANGE")
                    SET(_fct_C "${_fct}")

                ELSEIF(${FORTRAN_MANGLING_DETECTED} MATCHES "-DfcIsF2C")
                    SET(_fct_C "${_fct}__")

                ELSEIF(${FORTRAN_MANGLING_DETECTED} MATCHES "-DUPCASE")
                    STRING(TOUPPER "${_fct}" _fct_C)

                ELSE()
                    SET(_fct_C "${_fct}")

                ENDIF()
                CHECK_FUNCTION_EXISTS("${_fxt_C}" BLAS_${_prefix}_${_fct}_WORKS)

            ENDIF(_CHECK_FORTRAN)

            # break at first error occurred
            IF(LAPACK_${_prefix}_${_fct}_WORKS)
                SET(LAPACK_${_prefix}_WORKS TRUE)
            ELSE(LAPACK_${_prefix}_${_fct}_WORKS)
                SET(LAPACK_${_prefix}_WORKS FALSE)
                BREAK()
            ENDIF()
        ENDFOREACH()

    ENDIF(LAPACK_${_prefix}_LIBRARIES)

    # Record the fact that LAPACK_${prefix} works
    IF(LAPACK_${_prefix}_WORKS)
        # Add to the list of working vendor
        LIST(APPEND LAPACK_VENDOR_WORKING "${_prefix}")

        # We suppposed that the order we test all LAPACK implementation is our preference
        IF(NOT LAPACK_FOUND)
            SET(LAPACK_FOUND  TRUE      )
            SET(LAPACK_VENDOR ${_prefix})

        ENDIF(NOT LAPACK_FOUND)

    ELSE(LAPACK_${_prefix}_WORKS)
        SET(LAPACK_${_prefix}_LIBRARY      "LAPACK_${_prefix}_LIBRARY-NOTFOUND"     )
        SET(LAPACK_${_prefix}_LIBRARIES    "LAPACK_${_prefix}_LIBRARIES-NOTFOUND"   )
        SET(LAPACK_${_prefix}_LIBRARY_PATH "LAPACK_${_prefix}_LIBRARY_PATH-NOTFOUND")

    ENDIF(LAPACK_${_prefix}_WORKS)

    # Status of internal variables
    IF(MAGMA_DEBUG_CMAKE)
        MESSAGE(STATUS "  * debug:   - LAPACK_${_prefix}_LIBRARY_PATH : ${LAPACK_${_prefix}_LIBRARY_PATH}")
        MESSAGE(STATUS "  * debug:   - LAPACK_${_prefix}_LDFLAGS      : ${LAPACK_${_prefix}_LDFLAGS}"     )
        MESSAGE(STATUS "  * debug:   - LAPACK_${_prefix}_CHEEV        : ${LAPACK_${_prefix}_WORKS}")
        MESSAGE(STATUS "  * debug:")
    ENDIF(MAGMA_DEBUG_CMAKE)

ENDMACRO(LAPACK_LOCATE_AND_TEST)

# LAPACK requires BLAS
# --------------------
MESSAGE(STATUS "Looking for LAPACK - LAPACK requires BLAS")
IF(NOT BLAS_SET)
    FIND_PACKAGE(BLAS)
    IF(NOT ${BLAS_FOUND})
        MESSAGE(STATUS "Looking for LAPACK - LAPACK requires BLAS but BLAS not found")
        MESSAGE(STATUS "Looking for LAPACK - not found")
        RETURN()
    ENDIF(NOT ${BLAS_FOUND})
ENDIF()


# Define LAPACK_TESTED_VENDOR
# -------------------------
IF(${BLAS_FOUND})
    IF(NOT "${MAGMA_USE_LAPACK}" STREQUAL  "^$")
        MESSAGE(STATUS "Looking for LAPACK - LAPACK tested vendor: ${LAPACK_VENDOR}")
        STRING(TOUPPER "${MAGMA_USE_LAPACK}" LAPACK_TESTED_VENDOR)

    ELSE()
        IF("^${BLAS_VENDOR}$" STREQUAL "MKL"        OR
           "^${BLAS_VENDOR}$" STREQUAL "MKL_SEQ"    OR
           "^${BLAS_VENDOR}$" STREQUAL "MKL_RT"     OR
           "^${BLAS_VENDOR}$" STREQUAL "ACML"       OR
           "^${BLAS_VENDOR}$" STREQUAL "GOTO"       OR
           "^${BLAS_VENDOR}$" STREQUAL "ACCELERATE" OR
           "^${BLAS_VENDOR}$" STREQUAL "VECLIB"       )
            MESSAGE(STATUS "Looking for LAPACK - LAPACK tested vendor: ${BLAS_VENDOR}")
            STRING(TOUPPER "${BLAS_VENDOR}" LAPACK_TESTED_VENDOR)
        ENDIF()

    ENDIF()
ENDIF(${BLAS_FOUND})

# Put BLAS_LDFLAGS in a the list LAPACK_EXTRA_FLAGS
# -------------------------------------------------
IF(${BLAS_FOUND})
    STRING(REPLACE " " ";" LAPACK_EXTRA_FLAGS "${BLAS_LDFLAGS}")
    MARK_AS_ADVANCED(LAPACK_EXTRA_FLAGS)
ENDIF(${BLAS_FOUND})

# Initialized values
# ------------------
SET(LAPACK_FOUND FALSE)
SET(LAPACK_VENDOR_WORKING "")
IF(NOT LAPACK_TESTED_VENDOR)
    SET(TEST_ALL_LAPACK_VENDOR ON)
ENDIF(NOT LAPACK_TESTED_VENDOR)

# Try to find LAPACK in BLAS given by the user
# --------------------------------------------
IF(BLAS_SET)
    # Locate and test this implementation
    LAPACK_LOCATE_AND_TEST("user's desire" "${BLAS_LIBRARIES}" "" "")
ENDIF(BLAS_SET)


# Intel MKL (sequential)
# ----------------------
IF(LAPACK_TESTED_VENDOR STREQUAL "MKL_SEQ" OR LAPACK_TESTED_VENDOR STREQUAL "MKL" OR TEST_ALL_LAPACK_VENDOR)
    # Determine the architecture
    IF(SIZEOF_INTEGER EQUAL 4)
        SET(_LAPACK_MKL_PRECISION "_lp64" )
    ELSEIF(SIZEOF_INTEGER EQUAL 8)
        SET(_LAPACK_MKL_PRECISION "_ilp64" )
    ENDIF()

    # Locate and test this implementation
    LAPACK_LOCATE_AND_TEST(
        "MKL_SEQ"
        "mkl_intel${_LAPACK_MKL_PRECISION};mkl_sequential;mkl_core"
        "-lm"
        "${LAPACK_EXTRA_FLAGS}")

ENDIF(LAPACK_TESTED_VENDOR STREQUAL "MKL_SEQ" OR LAPACK_TESTED_VENDOR STREQUAL "MKL" OR TEST_ALL_LAPACK_VENDOR)


# Intel MKL (multi-threaded)
# --------------------------
IF(LAPACK_TESTED_VENDOR STREQUAL "MKL_MT" OR TEST_ALL_LAPACK_VENDOR)
    # Determine the architecture
    IF(SIZEOF_INTEGER EQUAL 4)
        SET(_LAPACK_MKL_PRECISION "_lp64" )
    ELSEIF(SIZEOF_INTEGER EQUAL 8)
        SET(_LAPACK_MKL_PRECISION "_ilp64" )
    ENDIF()

    # Locate and test this implementation depends on threading layer
    # --------------------------------------------------------------
    IF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        LAPACK_LOCATE_AND_TEST(
            "MKL_MT"
            "mkl_intel${_LAPACK_MKL_PRECISION};mkl_intel_thread;mkl_core"
            "-openmp;-lm"
            "${LAPACK_EXTRA_FLAGS}")

    ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        LAPACK_LOCATE_AND_TEST(
            "MKL_MT"
            "mkl_intel${_LAPACK_MKL_PRECISION};mkl_gnu_thread;mkl_core"
            "-fopenmp;-lm"
            "${LAPACK_EXTRA_FLAGS}")

    ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
        LAPACK_LOCATE_AND_TEST(
            "MKL_MT"
            "mkl_intel${_LAPACK_MKL_PRECISION};mkl_pgi_thread;mkl_core"
            "-mp;-lm"
            "${LAPACK_EXTRA_FLAGS}")

    ENDIF()

ENDIF(LAPACK_TESTED_VENDOR STREQUAL "MKL_MT" OR TEST_ALL_LAPACK_VENDOR)

# AMD ACML 
# --------
IF(LAPACK_TESTED_VENDOR STREQUAL "ACML.*" OR TEST_ALL_LAPACK_VENDOR)
    # Locate and test this implementation
    IF(LAPACK_TESTED_VENDOR STREQUAL "ACML_MP" )
        SET(LAPACK_ACML_MP_DIR ${BLAS_DIR})
        LAPACK_LOCATE_AND_TEST("ACML_MP" "acml_mp;acml_mv" "" "${LAPACK_EXTRA_FLAGS}")

    ELSE() #IF(LAPACK_TESTED_VENDOR STREQUAL "ACML")
        SET(LAPACK_ACML_DIR ${BLAS_DIR})
        LAPACK_LOCATE_AND_TEST("ACML" "acml;acml_mv" "" "${LAPACK_EXTRA_FLAGS}")

    ENDIF()


ENDIF(LAPACK_TESTED_VENDOR STREQUAL "ACML.*" OR TEST_ALL_LAPACK_VENDOR)


# GotoBLAS2 (http://www.tacc.utexas.edu/tacc-projects/gotoblas2)
# --------------------------------------------------------------
IF(LAPACK_TESTED_VENDOR STREQUAL "GOTO" OR TEST_ALL_LAPACK_VENDOR)
    # Locate and test this implementation
    LAPACK_LOCATE_AND_TEST("GOTO" "goto2" "" "${LAPACK_EXTRA_FLAGS}")

ENDIF(LAPACK_TESTED_VENDOR STREQUAL "GOTO" OR TEST_ALL_LAPACK_VENDOR)


# Eigen (http://eigen.tuxfamily.org)
# ----------------------------------
IF(LAPACK_TESTED_VENDOR STREQUAL "EIGEN" OR TEST_ALL_LAPACK_VENDOR)
    # Locate and test this implementation
    LAPACK_LOCATE_AND_TEST("EIGEN" "eigen_lapack;lapack" "" "${LAPACK_EXTRA_FLAGS}")

ENDIF(LAPACK_TESTED_VENDOR STREQUAL "EIGEN" OR TEST_ALL_LAPACK_VENDOR)


# APPLE Accelerate
# ----------------
IF(LAPACK_TESTED_VENDOR STREQUAL "ACCELERATE" OR TEST_ALL_LAPACK_VENDOR)
    # Locate and test this implementation
    IF(APPLE)
        LAPACK_LOCATE_AND_TEST("ACCELERATE" "Accelerate" "-framework Accelerate" "${LAPACK_EXTRA_FLAGS}")

    ENDIF(APPLE)

ENDIF(LAPACK_TESTED_VENDOR STREQUAL "ACCELERATE" OR TEST_ALL_LAPACK_VENDOR)

# VECLIB
# ------
IF(LAPACK_TESTED_VENDOR STREQUAL "VECLIB" OR TEST_ALL_LAPACK_VENDOR)
    # Locate and test this implementation
    IF(NOT APPLE)
        IF(SIZEOF_INTEGER EQUAL 4)
            LAPACK_LOCATE_AND_TEST("VECLIB" "veclib" "" "${LAPACK_EXTRA_FLAGS}")

         ELSEIF(SIZEOF_INTEGER EQUAL 8)
            LAPACK_LOCATE_AND_TEST("VECLIB" "veclib8" "" "${LAPACK_EXTRA_FLAGS}")

         ENDIF()
    ENDIF(NOT APPLE)

ENDIF(LAPACK_TESTED_VENDOR STREQUAL "VECLIB" OR TEST_ALL_LAPACK_VENDOR)


# GENERIC LAPACK
# ------------
IF(LAPACK_TESTED_VENDOR STREQUAL "GENERIC" OR TEST_ALL_LAPACK_VENDOR)
    # Locate and test this implementation
    LAPACK_LOCATE_AND_TEST("GENERIC" "lapack" "" "${LAPACK_EXTRA_FLAGS}")

ENDIF(LAPACK_TESTED_VENDOR STREQUAL "GENERIC" OR TEST_ALL_LAPACK_VENDOR)


# Information about we find out
# -----------------------------
IF(${LAPACK_VENDOR})
    MESSAGE(STATUS "Looking for LAPACK - not found")
ELSE(${LAPACK_VENDOR})
    MESSAGE(STATUS "Looking for LAPACK - found vendors: ${LAPACK_VENDOR_WORKING}")
    MESSAGE(STATUS "Looking for LAPACK - picked vendor: ${LAPACK_VENDOR}")
ENDIF(${LAPACK_VENDOR})

# Set variable for LAPACK
# ---------------------
IF(LAPACK_FOUND)
    # Set status of LAPACK
    SET(LAPACK_USED_MODE "FIND")
    SET(HAVE_LAPACK ON         )

    # Set value for LAPACK
    SET(LAPACK_LDFLAGS      "${LAPACK_${LAPACK_VENDOR}_LDFLAGS}"     )
    SET(LAPACK_LIBRARY      "${LAPACK_${LAPACK_VENDOR}_LIBRARY}"     )
    SET(LAPACK_LIBRARIES    "${LAPACK_${LAPACK_VENDOR}_LIBRARIES}"   )
    SET(LAPACK_LIBRARY_PATH "${LAPACK_${LAPACK_VENDOR}_LIBRARY_PATH}")

    # Status of internal variables
    IF(MAGMA_DEBUG_CMAKE)
        MESSAGE(STATUS "  * debug:")
        MESSAGE(STATUS "  * debug: Looking for LAPACK - picked vendor")
        MESSAGE(STATUS "  * debug:   - LAPACK_LIBRARY_PATH     : ${LAPACK_LIBRARY_PATH}")
        MESSAGE(STATUS "  * debug:   - LAPACK_LIBRARIES        : ${LAPACK_LIBRARIES}"   )
        MESSAGE(STATUS "  * debug:   - LAPACK_LDFLAGS          : ${LAPACK_LDFLAGS}"     )
        MESSAGE(STATUS "  * debug:")
    ENDIF(MAGMA_DEBUG_CMAKE)
ENDIF(LAPACK_FOUND)

##
## @end file FindLAPACK.cmake
##
