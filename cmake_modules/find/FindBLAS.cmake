###
#
#  @file FindBLAS.cmake
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
#  BLAS_VENDOR       : if set checks only the specified vendor, if not set checks all the possibilities
#  BLAS_FOUND        : set to true if a library implementing the BLAS interface is found
#
#  BLAS_LDFLAGS      : uncached string of all required linker flags.
#  BLAS_LIBRARY      : uncached list of library
#  BLAS_LIBRARIES    : uncached list of required linker flags (with -l and -L).
#  BLAS_LIBRARY_PATH : uncached path of the library directory of BLAS installation.
#
###
#
# You can choose your prefered BLAS with this option : -DMAGMA_USE_BLAS="VENDOR"
#    NB: "VENDOR" can be in uppercase or lowercase
#
# You can defined BLAS_DIR to help the search of BLAS
#
###
#
# The preferential order of the chosen BLAS is:
#     * MKL_SEQ OR MKL - Intel Math Kernel Library (sequential interface)
#     * MKL_MT         - Intel Math Kernel Library (multi-threaded interface)
#                         --> See http://software.intel.com/en-us/intel-mkl
#
#     * ACML           - Single threaded version of the AMD Core Math Librar
#     * ACML_MP        - Multithreaded version of the AMD Core Math Library using OpenMP
#                         --> See http://developer.amd.com/cpu/Libraries/acml
#
#     * EIGEN          - Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms
#                         --> See http://eigen.tuxfamily.org
#
#     * GOTO           - Goto BLAS v2
#                         --> See http://www.tacc.utexas.edu/tacc-projects/gotoblas2
#
#     * ATLAS          - Automatically Tuned Linear Algebra Software
#                         --> See http://math-atlas.sourceforge.net/
#
#     * SCS            - SGI's Scientific Computing Software Library
#     * SCSL           - SGI's Scientific Computing Software Library
#                         --> See http://www.sgi.com/products/software/irix/scsl.html
#
#     * SGIMATH        - SGIMATH library
#                         --> See ???
#
#     * SUNPERF        - Oracle Performance Library (formerly Sun Performance Library)
#                         --> See http://www.oracle.com/technetwork/server-storage/solarisstudio
#
#     * ESSLSMP        - IBM's Engineering and Scientific Subroutine Library (smp)
#     * ESSLIBM        - IBM's Engineering and Scientific Subroutine Library
#                         --> See http://www-03.ibm.com/systems/software/essl/
#
#     * ACCELERATE     - Apple's Accelerate library
#                         --> http://developer.apple.com/performance/accelerateframework
#
#     * VECLIB         - HP's Math Library: VECLIB
#                         --> See ???
#
#     * PHIPACK        - PhiPACK libraries
#                         --> See ???
#
#     * CXML           - Alpha CXML library
#                         --> See ???
#
#     * DXML           - Alpha DXML library
#                         --> See ???
#
#     * REFBLAS        - Search for refblas
#
#     * GENERIC        - Search for a generic libblas
#
###

# Early exit if already found or already tested
# ---------------------------------------------
IF(DEFINED BLAS_USED_MODE)
    MESSAGE(STATUS "Looking for BLAS - mode ${${_DEPVAR}_USED_MODE} already defined")
    RETURN()
ENDIF(DEFINED BLAS_USED_MODE)

# Early exit if already searched
# ------------------------------
IF(FINDBLAS_TESTED)
    IF(BLAS_FOUND)
        MESSAGE(STATUS "Looking for BLAS - already found")
    ELSE(BLAS_FOUND)
        MESSAGE(STATUS "Looking for BLAS - already NOT found")
    ENDIF(BLAS_FOUND)
    RETURN()
ENDIF(FINDBLAS_TESTED)
SET(FINDBLAS_TESTED TRUE)

# Message
# -------
MESSAGE(STATUS "Looking for BLAS")

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
    IF(BLAS_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Looking for BLAS - requires Fortran, C, or C++ to be enabled.")
    ELSE(BLAS_FIND_REQUIRED)
        MESSAGE(STATUS "Looking for BLAS - not found (unsupported languages)")
        RETURN()
    ENDIF(BLAS_FIND_REQUIRED)

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

UNSET(BLAS_EXTRA_PATH)
IF(NOT "^${BLAS_DIR}$" STREQUAL "^$")
    LIST(APPEND BLAS_EXTRA_PATH ${BLAS_DIR})
ENDIF()
IF(NOT "^${_lib_env}$" STREQUAL "^$")
    LIST(APPEND BLAS_EXTRA_PATH ${_lib_env})
ENDIF()

# Determine the default integer size
# ----------------------------------
CHECK_FORTRAN_TYPE_SIZES()
IF(NOT SIZEOF_INTEGER)
    MESSAGE(WARNING "Looking for BLAS - Unable to determine default integer size")
    MESSAGE(WARNING "Looking for BLAS - Assuming integer*4")
    SET(SIZEOF_INTEGER 4)

ENDIF(NOT SIZEOF_INTEGER)

# Macro to locate a library and check for a specified symbol
# ----------------------------------------------------------
MACRO(BLAS_LOCATE_AND_TEST _prefix _blas_libnames _flags)

    # Status of internal variables
    IF(MAGMA_DEBUG_CMAKE)
        MESSAGE( STATUS "Looking for BLAS - try ${_prefix}")
        MESSAGE(STATUS "  * debug:")
        MESSAGE(STATUS "  * debug: Looking for BLAS - status of ${_prefix}")
        MESSAGE(STATUS "  * debug:   - libraries               : ${_blas_libnames}")
        MESSAGE(STATUS "  * debug:   - additional flags        : ${_flags}")
    ENDIF(MAGMA_DEBUG_CMAKE)

    # Initialize values
    SET(BLAS_${_prefix}_FOUND        FALSE)
    SET(BLAS_${_prefix}_LDFLAGS      ""   )
    SET(BLAS_${_prefix}_LIBRARY      ""   )
    SET(BLAS_${_prefix}_LIBRARIES    ""   )
    SET(BLAS_${_prefix}_LIBRARY_PATH ""   )
    MARK_AS_ADVANCED(BLAS_${_prefix}_LDFLAGS)
    MARK_AS_ADVANCED(BLAS_${_prefix}_LIBRARY)
    MARK_AS_ADVANCED(BLAS_${_prefix}_LIBRARIES)
    MARK_AS_ADVANCED(BLAS_${_prefix}_LIBRARY_PATH)

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
    UNSET(BLAS_SEARCH_PATH)
    IF(NOT "^${BLAS_${_prefix}_DIR}$" STREQUAL "^$")
        LIST(APPEND BLAS_SEARCH_PATH ${BLAS_${_prefix}_DIR})
    ENDIF()
    IF(NOT "^${BLAS_EXTRA_PATH}$" STREQUAL "^$")
        LIST(APPEND BLAS_SEARCH_PATH ${BLAS_EXTRA_PATH})
    ENDIF()
    IF(MAGMA_DEBUG_CMAKE)
        MESSAGE(STATUS "  * debug:   - searching in this paths : ${BLAS_SEARCH_PATH}")
    ENDIF(MAGMA_DEBUG_CMAKE)

    # Define path we have to looking for first
    SET(CMAKE_PREFIX_PATH ${BLAS_SEARCH_PATH})

    # Looking for path of libraries
    SET({BLAS_${_prefix}_LDFLAGS ${_blas_libnames})
    FOREACH(_library ${_blas_libnames})
        FIND_LIBRARY(BLAS_${_prefix}_${_library}_LIBRARY
                     NAMES ${_library}
                     PATHS ${BLAS_SEARCH_PATH}
                    )
        MARK_AS_ADVANCED(BLAS_${_prefix}_${_library}_LIBRARY)

        IF(MAGMA_DEBUG_CMAKE)
            MESSAGE(STATUS "  * debug:   - obtained library path   : ${BLAS_${_prefix}_${_library}_LIBRARY}")
        ENDIF(MAGMA_DEBUG_CMAKE)

        IF(BLAS_${_prefix}_${_library}_LIBRARY)
            LIST(APPEND BLAS_${_prefix}_LIBRARY   "${BLAS_${_prefix}_${_library}_LIBRARY}")
            LIST(APPEND BLAS_${_prefix}_LIBRARIES "${_library}"                           )
            GET_FILENAME_COMPONENT(BLAS_${_prefix}_${_library}_LIBRARY_FILENAME
                                   ${BLAS_${_prefix}_${_library}_LIBRARY}
                                   NAME)
            STRING(REGEX REPLACE "(.*)/${BLAS_${_prefix}_${_library}_LIBRARY_FILENAME}" "\\1"
                   BLAS_${_prefix}_${_library}_LIBRARY_PATH "${BLAS_${_prefix}_${_library}_LIBRARY}")
            LIST(APPEND BLAS_${_prefix}_LIBRARY_PATH   "${BLAS_${_prefix}_${_library}_LIBRARY_PATH}")
            LIST(REMOVE_DUPLICATES BLAS_${_prefix}_LIBRARY_PATH)

        ELSE()
            SET(BLAS_${_prefix}_LIBRARY      "BLAS_${_prefix}_LIBRARY-NOTFOUND"     )
            SET(BLAS_${_prefix}_LIBRARIES    "BLAS_${_prefix}_LIBRARIES-NOTFOUND"   )
            SET(BLAS_${_prefix}_LIBRARY_PATH "BLAS_${_prefix}_LIBRARY_PATH-NOTFOUND")
            BREAK()
        ENDIF()

    ENDFOREACH()

    # Define path we have to looking for first
    UNSET(CMAKE_PREFIX_PATH)

    # Test this combination of libraries.
    IF(BLAS_${_prefix}_LIBRARIES)
        # Set CMAKE_REQUIRED_LIBRARIES to test the function(s)
        UNSET(CMAKE_REQUIRED_LIBRARIES)
        SET(BLAS_${_prefix}_EXTRA_FLAGS ${_flags})
        SET(BLAS_${_prefix}_LDFLAGS      ""      )
        FOREACH(_path ${BLAS_${_prefix}_LIBRARY_PATH})
            SET(BLAS_${_prefix}_LDFLAGS "${BLAS_${_prefix}_LDFLAGS} -L${_path}")
        ENDFOREACH()
        FOREACH(_flag ${BLAS_${_prefix}_LIBRARIES})
            SET(BLAS_${_prefix}_LDFLAGS "${BLAS_${_prefix}_LDFLAGS} -l${_flag}")
        ENDFOREACH()
        FOREACH(_flag ${BLAS_${_prefix}_EXTRA_FLAGS})
            SET(BLAS_${_prefix}_LDFLAGS "${BLAS_${_prefix}_LDFLAGS} ${_flag}")
        ENDFOREACH()
        STRING(REGEX REPLACE "^ " "" BLAS_${_prefix}_LDFLAGS "${BLAS_${_prefix}_LDFLAGS}")
        SET(CMAKE_REQUIRED_LIBRARIES "${BLAS_${_prefix}_LDFLAGS}")
        STRING(REPLACE " " ";" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

        # Load infoPACKAGE
        INCLUDE(infoBLAS)
        BLAS_INFO_FIND()

        # Test [X]gemm
        FOREACH(_fct ${BLAS_name_fct_test})
            # test gemm at a given precision
            IF(_CHECK_FORTRAN)
                CHECK_FORTRAN_FUNCTION_EXISTS("${_fct}" BLAS_${_prefix}_${_fct}_WORKS)

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
                CHECK_FUNCTION_EXISTS("${_fxt_C}" BLAS_${_prefix}_${fct}_WORKS)

            ENDIF(_CHECK_FORTRAN)

            # break at first error occurred
            IF(BLAS_${_prefix}_${_fct}_WORKS)
                SET(BLAS_${_prefix}_WORKS TRUE)
            ELSE(BLAS_${_prefix}_${_fct}_WORKS)
                SET(BLAS_${_prefix}_WORKS FALSE)
                BREAK()
            ENDIF()
        ENDFOREACH()

    ENDIF(BLAS_${_prefix}_LIBRARIES)

    # Record the fact that BLAS_${prefix} works
    IF(BLAS_${_prefix}_WORKS)
        # Add to the list of working vendor
        LIST(APPEND BLAS_VENDOR_WORKING "${_prefix}")

        # We suppposed that the order we test all BLAS implementation is our preference
        IF(NOT BLAS_FOUND)
            SET(BLAS_FOUND  TRUE      )
            SET(BLAS_VENDOR ${_prefix})

        ENDIF(NOT BLAS_FOUND)

    ELSE(BLAS_${_prefix}_WORKS)
        SET(BLAS_${_prefix}_LIBRARY      "BLAS_${_prefix}_LIBRARY-NOTFOUND"     )
        SET(BLAS_${_prefix}_LIBRARIES    "BLAS_${_prefix}_LIBRARIES-NOTFOUND"   )
        SET(BLAS_${_prefix}_LIBRARY_PATH "BLAS_${_prefix}_LIBRARY_PATH-NOTFOUND")

    ENDIF(BLAS_${_prefix}_WORKS)

    # Status of internal variables
    IF(MAGMA_DEBUG_CMAKE)
        MESSAGE(STATUS "  * debug:   - BLAS_${_prefix}_LIBRARY_PATH : ${BLAS_${_prefix}_LIBRARY_PATH}")
        MESSAGE(STATUS "  * debug:   - BLAS_${_prefix}_LDFLAGS      : ${BLAS_${_prefix}_LDFLAGS}"     )
        MESSAGE(STATUS "  * debug:   - BLAS_${_prefix}_DGEMM        : ${BLAS_${_prefix}_WORKS}")
        MESSAGE(STATUS "  * debug:")
    ENDIF(MAGMA_DEBUG_CMAKE)

ENDMACRO(BLAS_LOCATE_AND_TEST)

# Define BLAS_TESTED_VENDOR
# -------------------------
IF(NOT "^${MAGMA_USE_BLAS}$" STREQUAL "^$")
    STRING(TOUPPER "${MAGMA_USE_BLAS}" BLAS_TESTED_VENDOR)
    MESSAGE(STATUS "Looking for BLAS - BLAS tested vendor: ${BLAS_TESTED_VENDOR}")
ENDIF()

# Initialized values
# ------------------
SET(BLAS_FOUND FALSE)
SET(BLAS_VENDOR_WORKING "")
IF(NOT BLAS_TESTED_VENDOR)
    SET(TEST_ALL_BLAS_VENDOR ON)
ENDIF(NOT BLAS_TESTED_VENDOR)

# Intel MKL (sequential)
# ----------------------
IF(BLAS_TESTED_VENDOR STREQUAL "MKL_SEQ" OR BLAS_TESTED_VENDOR STREQUAL "MKL" OR TEST_ALL_BLAS_VENDOR)
    # Determine the architecture
    IF(SIZEOF_INTEGER EQUAL 4)
        SET(_BLAS_MKL_PRECISION "_lp64" )
    ELSEIF(SIZEOF_INTEGER EQUAL 8)
        SET(_BLAS_MKL_PRECISION "_ilp64" )
    ENDIF()

    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST(
        "MKL_SEQ"
        "mkl_intel${_BLAS_MKL_PRECISION};mkl_sequential;mkl_core"
        "-lpthread;-lm")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "MKL_SEQ" OR BLAS_TESTED_VENDOR STREQUAL "MKL" OR TEST_ALL_BLAS_VENDOR)


# Intel MKL (multi-threaded)
# --------------------------
IF(BLAS_TESTED_VENDOR STREQUAL "MKL_MT" OR TEST_ALL_BLAS_VENDOR)
    # Determine the architecture
    IF(SIZEOF_INTEGER EQUAL 4)
        SET(_BLAS_MKL_PRECISION "_lp64" )
    ELSEIF(SIZEOF_INTEGER EQUAL 8)
        SET(_BLAS_MKL_PRECISION "_ilp64" )
    ENDIF()

    # Locate and test this implementation depends on threading layer
    # --------------------------------------------------------------
    IF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        BLAS_LOCATE_AND_TEST(
            "MKL_MT"
            "mkl_intel${_BLAS_MKL_PRECISION};mkl_intel_thread;mkl_core"
            "-openmp;-lpthread;-lm")

    ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        BLAS_LOCATE_AND_TEST(
            "MKL_MT"
            "mkl_gf${_BLAS_MKL_PRECISION};mkl_gnu_thread;mkl_core"
            "-fopenmp;-lpthread;-lm")

    ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
        BLAS_LOCATE_AND_TEST(
            "MKL_MT"
            "mkl_intel${_BLAS_MKL_PRECISION};mkl_pgi_thread;mkl_core"
            "-mp;-lpthread;-lm")

    ENDIF()

ENDIF(BLAS_TESTED_VENDOR STREQUAL "MKL_MT" OR TEST_ALL_BLAS_VENDOR)

# AMD ACML
# --------
IF(BLAS_TESTED_VENDOR STREQUAL "ACML" OR BLAS_TESTED_VENDOR STREQUAL "ACML_MP" OR TEST_ALL_BLAS_VENDOR)
    # Looking for EULA
    IF(WIN32)
        FILE(GLOB _ACML_ROOT "C:/AMD/acml*/ACML-EULA.txt")
    ELSE(WIN32)
        FILE( GLOB _ACML_ROOT "/opt/acml*/ACML-EULA.txt")
    ENDIF(WIN32)

    # Define parameters
    IF(_ACML_ROOT)
        # Looking for ACML prefix
        GET_FILENAME_COMPONENT(_ACML_ROOT ${_ACML_ROOT} PATH)
        IF( SIZEOF_INTEGER EQUAL 8)
            SET(_ACML_PATH_SUFFIX "_int64")
        ELSE()
            SET(_ACML_PATH_SUFFIX "")
        ENDIF()
        IF(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
            SET(_ACML_COMPILER32 "IFort32")
            SET(_ACML_COMPILER64 "IFort64")
        ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "SunPro")
            SET( _ACML_COMPILER32 "sun32")
            SET( _ACML_COMPILER64 "sun64")
        ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
            SET(_ACML_COMPILER32 "pgi32")
            IF(WIN32)
                SET(_ACML_COMPILER64 "win64")
            ELSE()
                SET(_ACML_COMPILER64 "pgi64")
            ENDIF()
        ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "Open64")
            # 32 bit builds not supported on Open64 but for code simplicity
            # We'll just use the same directory twice
            SET(_ACML_COMPILER32 "open64_64")
            SET(_ACML_COMPILER64 "open64_64")
        ELSEIF(CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
            SET(_ACML_COMPILER32 "nag32")
            SET(_ACML_COMPILER64 "nag64")
        ELSE() #IF(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
            SET(_ACML_COMPILER32 "gfortran32")
            SET(_ACML_COMPILER64 "gfortran64")
        ENDIF()

        IF(BLAS_TESTED_VENDOR STREQUAL "ACML_MP")
            SET(BLAS_ACML_MP_DIR
                "${_ACML_ROOT}/${_ACML_COMPILER32}_mp${_ACML_PATH_SUFFIX}/lib"
                "${_ACML_ROOT}/${_ACML_COMPILER64}_mp${_ACML_PATH_SUFFIX}/lib")
        ELSE() #IF(BLAS_TESTED_VENDOR STREQUAL "ACML")
            SET(BLAS_ACML_DIR
                "${_ACML_ROOT}/${_ACML_COMPILER32}${_ACML_PATH_SUFFIX}/lib"
                "${_ACML_ROOT}/${_ACML_COMPILER64}${_ACML_PATH_SUFFIX}/lib")
        ENDIF()
    ENDIF()

    # Locate and test this implementation
    IF(BLAS_TESTED_VENDOR STREQUAL "ACML_MP" )
        BLAS_LOCATE_AND_TEST("ACML_MP" "acml_mp;acml_mv" "")

    ELSE() #IF(_BLAS_TESTED_VENDOR STREQUAL "ACML")
        BLAS_LOCATE_AND_TEST("ACML" "acml;acml_mv" "")

    ENDIF()

ENDIF(BLAS_TESTED_VENDOR STREQUAL "ACML" OR BLAS_TESTED_VENDOR STREQUAL "ACML_MP" OR TEST_ALL_BLAS_VENDOR)


# Eigen (http://eigen.tuxfamily.org)
# ----------------------------------
IF(BLAS_TESTED_VENDOR STREQUAL "EIGEN" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST("EIGEN" "eigen_blas" "")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "EIGEN" OR TEST_ALL_BLAS_VENDOR)


# GotoBLAS2 (http://www.tacc.utexas.edu/tacc-projects/gotoblas2)
# --------------------------------------------------------------
IF(BLAS_TESTED_VENDOR STREQUAL "GOTO" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST("GOTO" "goto2" "")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "GOTO" OR TEST_ALL_BLAS_VENDOR)


# ATLAS (http://math-atlas.sourceforge.net/)
# ------------------------------------------
IF(BLAS_TESTED_VENDOR STREQUAL "ATLAS" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST("ATLAS" "f77blas;atlas" "")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "ATLAS" OR TEST_ALL_BLAS_VENDOR)


# SGI
# ---
IF(BLAS_TESTED_VENDOR STREQUAL "SCS" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    IF(SIZEOF_INTEGER EQUAL 4)
        BLAS_LOCATE_AND_TEST("SCS" "scs" "")

    ELSEIF(SIZEOF_INTEGER EQUAL 8)
        BLAS_LOCATE_AND_TEST("SCS" "scs_i8" "")

    ENDIF()

ENDIF(BLAS_TESTED_VENDOR STREQUAL "SCS" OR TEST_ALL_BLAS_VENDOR)

# SGI/Cray Scientific Library
# ---------------------------
IF(BLAS_TESTED_VENDOR STREQUAL "SCSL" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST("SCSL" "scsl" "")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "SCSL" OR TEST_ALL_BLAS_VENDOR)


# SGIMATH library
# ---------------
IF(BLAS_TESTED_VENDOR STREQUAL "SGIMATH" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST("SGIMATH" "complib.sgimath" "")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "SGIMATH" OR TEST_ALL_BLAS_VENDOR)


# Sun / Oracle PerfLib
# --------------------
IF(BLAS_TESTED_VENDOR STREQUAL "SUNPERF" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    IF(CMAKE_Fortran_COMPILER_ID MATCHES "SunPro")
        BLAS_LOCATE_AND_TEST("SUNPERF" "sunperf;sunmath" "-xlic_lib=sunperf")

    ELSE(CMAKE_Fortran_COMPILER_ID MATCHES "SunPro")
        BLAS_LOCATE_AND_TEST("SUNPERF" "sunperf;mtsk" "")

    ENDIF(CMAKE_Fortran_COMPILER_ID MATCHES "SunPro")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "SUNPERF" OR TEST_ALL_BLAS_VENDOR)


# IBM ESSL (SMP Version)
# ----------------------
IF(BLAS_TESTED_VENDOR STREQUAL "ESSLSMP" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    IF(SIZEOF_INTEGER EQUAL 4)
        BLAS_LOCATE_AND_TEST("ESSLSMP" "esslsmp;blas" "")

    ELSEIF(SIZEOF_INTEGER EQUAL 8)
        BLAS_LOCATE_AND_TEST("ESSLSMP" "esslsmp6464;blas" "")

    ENDIF()

ENDIF(BLAS_TESTED_VENDOR STREQUAL "ESSLSMP" OR TEST_ALL_BLAS_VENDOR)


# IBM ESSL
# --------
IF(BLAS_TESTED_VENDOR STREQUAL "ESSLIBM" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    IF(SIZEOF_INTEGER EQUAL 4)
        BLAS_LOCATE_AND_TEST("ESSLIBM" "essl;blas" "")

    ELSEIF(SIZEOF_INTEGER EQUAL 8)
        BLAS_LOCATE_AND_TEST("ESSLIBM" "essl6464;blas" "")

    ENDIF()

ENDIF(BLAS_TESTED_VENDOR STREQUAL "ESSLIBM" OR TEST_ALL_BLAS_VENDOR)


# APPLE Accelerate
# ----------------
IF(BLAS_TESTED_VENDOR STREQUAL "ACCELERATE" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    IF(APPLE)
        BLAS_LOCATE_AND_TEST("ACCELERATE" "Accelerate" "-framework Accelerate")

    ENDIF(APPLE)

ENDIF(BLAS_TESTED_VENDOR STREQUAL "ACCELERATE" OR TEST_ALL_BLAS_VENDOR)


# VECLIB
# ------
IF(BLAS_TESTED_VENDOR STREQUAL "VECLIB" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    IF(NOT APPLE)
        IF(SIZEOF_INTEGER EQUAL 4)
            BLAS_LOCATE_AND_TEST("VECLIB" "veclib" "")

         ELSEIF(SIZEOF_INTEGER EQUAL 8)
            BLAS_LOCATE_AND_TEST("VECLIB" "veclib8" "")

         ENDIF()
    ENDIF(NOT APPLE)

ENDIF(BLAS_TESTED_VENDOR STREQUAL "VECLIB" OR TEST_ALL_BLAS_VENDOR)


# PhiPACK libraries
# -----------------
IF(BLAS_TESTED_VENDOR STREQUAL "PHIPACK" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST("PHIPACK" "sgemm;dgemm;blas" "")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "PHIPACK" OR TEST_ALL_BLAS_VENDOR)


# Alpha CXML library
# ------------------
IF(BLAS_TESTED_VENDOR STREQUAL "CXML" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST("CXML" "cxml" "")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "CXML" OR TEST_ALL_BLAS_VENDOR)


# Alpha DXML library
# ------------------
IF(BLAS_TESTED_VENDOR STREQUAL "DXML" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST("DXML" "dxml" "")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "DXML" OR TEST_ALL_BLAS_VENDOR)


# REFBLAS library
# ---------------
IF(BLAS_TESTED_VENDOR STREQUAL "REFBLAS" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST("REFBLAS" "refblas" "")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "REFBLAS" OR TEST_ALL_BLAS_VENDOR)


# GENERIC BLAS
# ------------
IF(BLAS_TESTED_VENDOR STREQUAL "GENERIC" OR TEST_ALL_BLAS_VENDOR)
    # Locate and test this implementation
    BLAS_LOCATE_AND_TEST("GENERIC" "blas" "")

ENDIF(BLAS_TESTED_VENDOR STREQUAL "GENERIC" OR TEST_ALL_BLAS_VENDOR)


# Set variable for BLAS
# ---------------------
IF(BLAS_FOUND)
    # Set status of BLAS
    MESSAGE(STATUS "Looking for BLAS - found vendors: ${BLAS_VENDOR_WORKING}")
    MESSAGE(STATUS "Looking for BLAS - picked vendor: ${BLAS_VENDOR}")
    IF("^${MAGMA_USE_CUDA}$" STREQUAL "^$")
        SET(MAGMA_USE_CUDA "${BLAS_VENDOR}" CACHE STRING
            "Enable/Disable CUDA dependency (ON/OFF/<not-defined>)" FORCE)
     ENDIF()    
    SET(BLAS_USED_MODE "FIND")
    SET(HAVE_BLAS ON         )

    # Set value for BLAS
    SET(BLAS_LDFLAGS      "${BLAS_${BLAS_VENDOR}_LDFLAGS}"     )
    SET(BLAS_LIBRARY      "${BLAS_${BLAS_VENDOR}_LIBRARY}"     )
    SET(BLAS_LIBRARIES    "${BLAS_${BLAS_VENDOR}_LIBRARIES}"   )
    SET(BLAS_LIBRARY_PATH "${BLAS_${BLAS_VENDOR}_LIBRARY_PATH}")

    # Status of internal variables
    IF(MAGMA_DEBUG_CMAKE)
        MESSAGE(STATUS "  * debug:")
        MESSAGE(STATUS "  * debug: Looking for BLAS - picked vendor"                )
        MESSAGE(STATUS "  * debug:   - BLAS_LIBRARY_PATH     : ${BLAS_LIBRARY_PATH}")
        MESSAGE(STATUS "  * debug:   - BLAS_LIBRARIES        : ${BLAS_LIBRARIES}"   )
        MESSAGE(STATUS "  * debug:   - BLAS_LDFLAGS          : ${BLAS_LDFLAGS}"     )
        MESSAGE(STATUS "  * debug:")
    ENDIF(MAGMA_DEBUG_CMAKE)

ELSE(BLAS_FOUND)
    MESSAGE(STATUS "Looking for BLAS - not found")

ENDIF(BLAS_FOUND)

##
## @end file FindBLAS.cmake
##
