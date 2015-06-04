###
#
# @file      : CMakeLists.txt
#
# @description   : Project MAGMA
#
# @version       : 
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : lun. 11 juin 2012 16:08:54 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#############################################
#                                           #
#      Force compilation out-of-source      #
#                                           #
#############################################
STRING(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" MAGMA_COMPILE_INPLACE)
IF(MAGMA_COMPILE_INPLACE)
    MESSAGE(FATAL_ERROR "${PROJECT_NAME} requires an out of source build.
  1/ remove CMakeCache.txt and CMakeFiles in ${CMAKE_SOURCE_DIR}
  2/ do mkdir build && cd build and tape cmake ../")
ENDIF(MAGMA_COMPILE_INPLACE)

#############################################
#                                           #
#        Definition of the project          #
#                                           #
#############################################
PROJECT(MAGMA Fortran C CXX)

#############################################
#                                           #
#      Systems checks and parameters        #
#                                           #
#############################################

# CMake can find FindFoo.cmake in this directory
# ----------------------------------------------
SET(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake_modules
    ${CMAKE_SOURCE_DIR}/cmake_modules/info
    ${CMAKE_SOURCE_DIR}/cmake_modules/find
    ${CMAKE_SOURCE_DIR}/cmake_modules/install
    ${CMAKE_SOURCE_DIR}/cmake_modules/generic)

# Set the version of this project
# -------------------------------
INCLUDE(VersionMAGMA)

# Active cmake options
# --------------------
INCLUDE(FortranMangling)
INCLUDE(CMakeDetermineSystem)
INCLUDE(CheckCCompilerFlag)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckIncludeFiles)
INCLUDE(definePACKAGE)
INCLUDE(Summary)

# Active CTEST
# ------------
INCLUDE(CTest)
ENABLE_TESTING()

# Define the default install directory
# ------------------------------------
IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install directory used by install" FORCE)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# Organize output files in a tree
# -------------------------------
#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${MAGMA_BINARY_DIR}/bin)
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${MAGMA_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${MAGMA_BINARY_DIR}/lib)

# Set the build type
# ------------------
IF(NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release 
        CACHE STRING "Choose the type of build, options are None, Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#############################################
#                                           #
#      Definition of the MAGMA options      #
#                                           #
#############################################

# Misc options
# ------------
OPTION(BUILD_SHARED_LIBS
       "Build shared libraries" OFF)
OPTION(BUILD_64bits
       "Build 64 bits mode" ON)

# Versions of MAGMA which need to be compiled
# -------------------------------------------
INCLUDE(ConfigMAGMA)
#INCLUDE(ConfigMORSE)

# MAGMA options
# -------------
OPTION(MAGMA_USE_PLASMA 
       "Enable/Disable kernels using plasma-like kernels in MAGMA" OFF )
OPTION(MAGMA_USE_FERMI   
       "Switch between Tesla (OFF) and Fermi (ON) cards" ON)

# MORSE options
# -------------
OPTION(MORSE_USE_MPI
       "Enable/Disable distributed version of MORSE" OFF)
#OPTION(MORSE_USE_MULTICORE
#  "Enable/Disable multi-threaded kernels in MORSE" OFF)
#MARK_AS_ADVANCED(MORSE_USE_MULTICORE)

OPTION(MORSE_SCHED_STARPU
      "Enable/Disable StarPU scheduler in MORSE" ON)
OPTION(MORSE_SCHED_QUARK
      "Enable/Disable Quark scheduler in MORSE" OFF)

# MORSE options
# -------------
OPTION(MORSE_SEPARATE_PROJECTS
       "Enable/disable the separation of different projects in the installation tree" OFF)

#############################################
#                                           #
#       Check options compatibility         #
#                                           #
#############################################

# Force check for externals libraries
# -----------------------------------
IF(MAGMA_MORSE)

    # If MAGMA_MORSE is true, MAGMA_MORSE need to take the priority over MAGMA_1GPU
    # So MAGMA_1GPU will be activated only if CUDA_FOUND is true
    SET(MAGMA_1GPU       OFF CACHE BOOL "Enable Build of magma 1gpu lib"      FORCE)
    SET(MAGMA_USE_PLASMA OFF CACHE BOOL "Enable PLASMA-like kernels in MAGMA" FORCE)
  
    # Require PLASMA
    #SET(REQUIRE_BLAS    ON)
    #SET(REQUIRE_CBLAS   ON)
    SET(REQUIRE_LAPACK   ON)
    #SET(REQUIRE_LAPACKE ON)
    SET(REQUIRE_HWLOC    ON)
    #SET(REQUIRE_QUARK   ON)
    SET(REQUIRE_PLASMA   ON)

  IF(MORSE_USE_CUDA)
    SET(MAGMA_1GPU       ON CACHE BOOL "Enable Build of magma 1gpu lib"      FORCE)
    SET(MAGMA_USE_PLASMA ON CACHE BOOL "Enable PLASMA-like kernels in MAGMA" FORCE)
    SET(REQUIRE_CUDA     ON)
    ADD_DEFINITIONS(-DMORSE_USE_CUDA)
  ENDIF(MORSE_USE_CUDA)

  IF(MORSE_USE_MPI)
    SET(REQUIRE_MPI      ON)
    ADD_DEFINITIONS(-DMORSE_USE_MPI)
  ENDIF(MORSE_USE_MPI)

  IF(MORSE_SCHED_STARPU)
     SET(REQUIRE_HWLOC   ON) 
     SET(REQUIRE_STARPU  ON) 
  ENDIF(MORSE_SCHED_STARPU)  

  IF(MORSE_SCHED_QUARK)
     SET(REQUIRE_HWLOC   ON) 
     #SET(REQUIRE_QUARK  ON) 
  ENDIF(MORSE_SCHED_QUARK)  

ENDIF()

IF(MAGMA_1GPU)
    SET(REQUIRE_CUDA     ON)
    #SET(REQUIRE_BLAS    ON)
    SET(REQUIRE_LAPACK   ON)
ENDIF(MAGMA_1GPU)


# PLASMA-like kernels
# -------------------
IF(MAGMA_USE_PLASMA)
    #SET(REQUIRE_BLAS    ON)
    #SET(REQUIRE_CBLAS   ON)
    SET(REQUIRE_LAPACK   ON)
    #SET(REQUIRE_LAPACKE ON)
    SET(REQUIRE_HWLOC    ON)
    #SET(REQUIRE_QUARK   ON)
    SET(REQUIRE_PLASMA   ON)

    ADD_DEFINITIONS(-DMAGMA_WITH_PLASMA)
ENDIF()

# Define the type of the GPU card
# -------------------------------
IF(MAGMA_USE_FERMI)
	ADD_DEFINITIONS(-DGPUSHMEM=200)
	SET(MAGMA_GPU_OPTS "-arch=sm_20")
ELSE(MAGMA_USE_FERMI)
	ADD_DEFINITIONS(-DGPUSHMEM=130)
	SET(MAGMA_GPU_OPTS "-arch=sm_13;--generate-code;arch=compute_13,code=compute_13;--generate-code;arch=compute_10,code=compute_10")
ENDIF(MAGMA_USE_FERMI)

# Detect the system and set specific options
# ------------------------------------------
INCLUDE(SystemDetection)

# Define Fortran mangling
# -----------------------
FortranMangling()

#############################################
#                                           #
#    Look for the required dependencies     #
#                                           #
#############################################

# PKG_CONFIG
# ----------
FIND_PACKAGE(PkgConfig QUIET)
IF(PKG_CONFIG_EXECUTABLE)
    MESSAGE(STATUS "Found pkg-config: ${PKG_CONFIG_EXECUTABLE}")
ELSE(PKG_CONFIG_EXECUTABLE)
    MESSAGE(STATUS "Looking for pkg-config - not found")
ENDIF(PKG_CONFIG_EXECUTABLE)

# PYTHON
# ------
FIND_PACKAGE(PythonInterp REQUIRED)

# CUDA
# ----
IF(REQUIRE_CUDA)
    MESSAGE(STATUS "Looking for CUDA - required dependency")
    FIND_PACKAGE(CUDA REQUIRED)
ELSE()
    IF("${MORSE_USE_CUDA}" MATCHES "OFF")
        MESSAGE(STATUS "Looking for CUDA - not asked")
    ELSE()
        MESSAGE(STATUS "Looking for CUDA - optional dependency")
        FIND_PACKAGE(CUDA QUIET)
        IF(CUDA_FOUND)
            SET(MORSE_USE_CUDA "ON"
                CACHE STRING "Enable/Disable CUDA dependency (ON/OFF/<not-defined>)")
        ELSE(CUDA_FOUND)
            SET(MORSE_USE_CUDA "OFF"
                CACHE STRING "Enable/Disable CUDA dependency (ON/OFF/<not-defined>)")
        ENDIF(CUDA_FOUND)
    ENDIF()
ENDIF()

IF(CUDA_FOUND)
    MESSAGE(STATUS "Looking for CUDA - found")
    IF(CUDA_VERSION VERSION_LESS "4.0")
        SET(CUDA_HOST_COMPILATION_CPP OFF)
    ENDIF(CUDA_VERSION VERSION_LESS "4.0")
    SET(CUDA_BUILD_EMULATION OFF)

    # Define all flags needed
    IF(${ARCH_X86_64})
        LINK_DIRECTORIES(${CUDA_TOOLKIT_ROOT_DIR}/lib64)
        SET(CUDA_LIBRARY_PATH "${CUDA_TOOLKIT_ROOT_DIR}/lib64")
    ELSE()
        LINK_DIRECTORIES(${CUDA_TOOLKIT_ROOT_DIR}/lib)
        SET(CUDA_LIBRARY_PATH "${CUDA_TOOLKIT_ROOT_DIR}/lib")
    ENDIF()
    SET(CUDA_LDFLAGS "-L${CUDA_LIBRARY_PATH} -lcudart -lcuda")
    SET(CUDA_LIBRARY "${CUDA_CUDART_LIBRARY};${CUDA_CUDA_LIBRARY}")
    SET(CUDA_INCLUDE_PATH "${CUDA_INCLUDE_DIRS}")
    INCLUDE_DIRECTORIES(${CUDA_INCLUDE_PATH})
    LIST(INSERT EXTRA_LIBS 0 ${CUDA_CUDA_LIBRARY})
    LIST(INSERT EXTRA_LIBS 0 ${CUDA_CUDART_LIBRARY})
    SET(HAVE_CUDA ON)
    SET(HAVE_CUBLAS ON)

    # Should move them to a magma_config.h.in 
    ADD_DEFINITIONS(-DHAVE_CUDA)
    ADD_DEFINITIONS(-DHAVE_CUBLAS)

    # FIX for MAGMA_MORSE priority
    IF(MAGMA_MORSE)
        SET(MAGMA_1GPU       ON CACHE BOOL "Enable Build of magma 1gpu lib"             FORCE)
        SET(MAGMA_USE_PLASMA ON CACHE BOOL "Enable PLASMA-like kernels in MAGMA"        FORCE)
        SET(MORSE_USE_CUDA   ON CACHE BOOL "Enable/Disable GPU kernel support in MORSE" FORCE)
        SET(REQUIRE_CUDA     ON)
        ADD_DEFINITIONS(-DMORSE_USE_CUDA)
    
    ENDIF(MAGMA_MORSE)

ELSE(CUDA_FOUND)
    MESSAGE(STATUS "Looking for CUDA - not found")
    SET(HAVE_CUDA OFF)

ENDIF(CUDA_FOUND)

# MPI
# ---
IF(REQUIRE_MPI)
    DEFINE_PACKAGE("MPI" "depends")
ELSE()
    DEFINE_PACKAGE("MPI" "recommands")
ENDIF()

# HWLOC
# -----
IF(REQUIRE_HWLOC)
    DEFINE_PACKAGE("HWLOC" "depends")
ENDIF(REQUIRE_HWLOC)

# FXT
# ---
DEFINE_PACKAGE("FXT" "suggests")

# StarPU
# ------
IF(REQUIRE_STARPU)
    DEFINE_PACKAGE("STARPU" "depends")
ELSE()
    DEFINE_PACKAGE("STARPU" "recommands")
ENDIF()

# PLASMA
# ------
IF(REQUIRE_PLASMA)
    DEFINE_PACKAGE("PLASMA" "depends")
ENDIF(REQUIRE_PLASMA)

# LAPACK / BLAS
#--------------
IF(REQUIRE_LAPACK AND NOT HAVE_LAPACK AND NOT HAVE_PLASMA)
    DEFINE_PACKAGE("LAPACK" "depends")
    SET(MAGMA_EXTRA_LDFLAGS "${LAPACK_LDFLAGS} ${BLAS_LDFLAGS}")
ENDIF()

#############################################
#                                           #
#          Define parameters for :          #
#           * magma.pc                      #
#           * morse-starpu.pc               #
#           * magma-quark.pc                #
#                                           #
#############################################
SET(MAGMA_REQUIRED "")
SET(MORSE_REQUIRED "")
IF(MAGMA_1GPU)
    SET(MORSE_REQUIRED "${MORSE_REQUIRED} magma")
ENDIF()
IF(REQUIRE_PLASMA)
    SET(MAGMA_REQUIRED "${MAGMA_REQUIRED} plasma")
    SET(MORSE_REQUIRED "${MORSE_REQUIRED} plasma")
ENDIF()
IF(MORSE_SCHED_STARPU)
    IF(MPI_FOUND OR DEFINED MPI_BUILD_MODE)
        SET(MORSE_REQUIRED "${MORSE_REQUIRED} libstarpumpi")
    ELSE()
        SET(MORSE_REQUIRED "${MORSE_REQUIRED} libstarpu")
    ENDIF()
ENDIF()
IF(REQUIRE_HWLOC)
    SET(MAGMA_REQUIRED "${MAGMA_REQUIRED} hwloc")
    SET(MORSE_REQUIRED "${MORSE_REQUIRED} hwloc")
ENDIF()
IF(FXT_FOUND OR DEFINED FXT_BUILD_MODE)
    SET(MORSE_REQUIRED "${MORSE_REQUIRED} fxt")
ENDIF()

#############################################
#                                           #
#        Compilation of the MAGMA           #
#                                           #
#############################################

# Define the pointer size for fortran compilation
# -----------------------------------------------
ADD_DEFINITIONS(\"-Dmagma_devptr_t=integer\(kind=${SIZEOF_VOID_PTR}\)\")

# Define precision supported by MAGMA
# -----------------------------------
SET(MAGMA_PRECISION "z;c;d;s;zc;ds")

# Define compilation flags according to the compiler
# --------------------------------------------------
IF(CUDA_FOUND)
    SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};${MAGMA_GPU_OPTS};--compiler-options;-fno-strict-aliasing;-DUNIX")
ENDIF()

# Add MAGMA includes to everything
# --------------------------------
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/include)

# Check for the subdirectories 
# ----------------------------
IF(MAGMA_1GPU)
    ADD_SUBDIRECTORY(include)
    ADD_SUBDIRECTORY(src)
    ADD_SUBDIRECTORY(magmablas)
    ADD_SUBDIRECTORY(testing)

    # Create file magma.pc
    # --------------------
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/magma.pc.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/lib/pkgconfig/magma.pc"
                   @ONLY)

    # installation
    # ------------
    INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib/pkgconfig/magma.pc
            DESTINATION lib/pkgconfig)

ENDIF()

IF(MAGMA_MORSE)
    #IF(MORSE_SCHED_QUARK AND NOT HAVE_QUARK)
    #    ADD_SUBDIRECTORY(quark)
    #ENDIF()
    ADD_SUBDIRECTORY(multi-gpu-dynamic)
ENDIF()

#############################################
#                                           #
#    Print the summary of configuration     #
#                                           #
#############################################
SUMMARY()

###
### END CMakeLists.txt
###
