###
#
#  @file CMakeLists.txt
#
#  @project MAGMA
#  MAGMA is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.1.0
#  @author Mathieu Faverge 
#  @author Cedric Castagnede
#  @date 13-07-2012
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE(CheckCCompilerFlag)
INCLUDE(CheckCSourceCompiles)
INCLUDE(CheckFunctionExists)
INCLUDE(CheckLibraryExists)
INCLUDE(CheckIncludeFiles)
INCLUDE(CheckStructHasMember)
INCLUDE(CheckTypeSize)
INCLUDE(BuildSystemTools)


#########################################
#                                       #
#     Check the capabilities of the     #
#      system we are building for       #
#                                       #
#########################################

# check for the CPU we build for
# ------------------------------
SITE_NAME(CMAKE_HOSTNAME)
MESSAGE(STATUS "Building on ${CMAKE_HOSTNAME}")
MESSAGE(STATUS "Building for target ${CMAKE_SYSTEM_NAME}")
STRING(REGEX MATCH "(i.86-*)|(athlon-*)|(pentium-*)" _mach_x86 ${CMAKE_SYSTEM_PROCESSOR})
IF (_mach_x86)
    MESSAGE(STATUS "Build for X86")
    SET(ARCH_X86 1)
ENDIF (_mach_x86)

STRING(REGEX MATCH "(x86_64-*)|(X86_64-*)|(AMD64-*)|(amd64-*)" _mach_x86_64 ${CMAKE_SYSTEM_PROCESSOR})
IF (_mach_x86_64)
    MESSAGE(STATUS "Build for X86_64")
#    OPTION(BUILD_64bits "Build 64 bits mode" ON)
    SET(ARCH_X86_64 1)
ENDIF (_mach_x86_64)

STRING(REGEX MATCH "(ppc-*)|(powerpc-*)" _mach_ppc ${CMAKE_SYSTEM_PROCESSOR})
IF (_mach_ppc)
    MESSAGE(STATUS "Build for PPC")
    SET(ARCH_PPC 1)
ENDIF (_mach_ppc)

# Looking for architecture parameters
# -----------------------------------
CHECK_TYPE_SIZE(void* SIZEOF_VOID_PTR)
IF("^${SIZEOF_VOID_PTR}$" STREQUAL "^8$")
    MESSAGE(STATUS "Build for 64 bits ")
    OPTION(BUILD_64bits "Build 64 bits mode" ON)
ELSE()
    MESSAGE(STATUS "Build for 32 bits")
    OPTION(BUILD_64bits "Build 64 bits mode" OFF)
ENDIF()

# Define extension for compiled library without CMake
# ---------------------------------------------------
IF(WIN32)
    IF(BUILD_SHARED_LIBS)
        SET(MAGMA_LIBRARY_EXTENSION "dll")
    ELSE()
        SET(MAGMA_LIBRARY_EXTENSION "lib")
    ENDIF()
ENDIF(WIN32)
IF(APPLE)
    IF(BUILD_SHARED_LIBS)
        SET(MAGMA_LIBRARY_EXTENSION "dylib")
    ELSE()
        SET(MAGMA_LIBRARY_EXTENSION "a")
    ENDIF()
ELSE(APPLE)
    IF(BUILD_SHARED_LIBS)
        SET(MAGMA_LIBRARY_EXTENSION "so")
    ELSE()
        SET(MAGMA_LIBRARY_EXTENSION "a")
    ENDIF()
ENDIF(APPLE)

#########################################
#                                       #
#    Check the capabilities of the      #
#    compiler we are building with      #
#                                       #
#########################################

# Looking for the type of build
# -----------------------------
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" TYPE)

# Check the you have GNU Compiler
# -------------------------------
IF(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    SET(HAVE_GNU ON)
ENDIF()

# Check the you have Intel Compiler
# -------------------------------
IF(CMAKE_C_COMPILER_ID STREQUAL "Intel")
    SET(HAVE_INTEL ON)
ENDIF()

# Check the you have PGI Compiler
# -------------------------------
IF(CMAKE_C_COMPILER_ID STREQUAL "PGI")
    SET(HAVE_PGI ON)
ENDIF()

# Intel tricks
# ------------
STRING(REGEX MATCH ".*icc$" _match_icc ${CMAKE_C_COMPILER})
IF(_match_icc)
    ADD_FLAGS(CMAKE_C_FLAGS "-diag-disable vec")
ENDIF(_match_icc)

STRING(REGEX MATCH ".*ifort$" _match_ifort ${CMAKE_Fortran_COMPILER})
IF(_match_ifort)
    MESSAGE(STATUS "Add -nofor_main to the Fortran linker")
    ADD_FLAGS(CMAKE_Fortran_FLAGS   "-diag-disable vec")
    ADD_FLAGS(CMAKE_Fortran_FLAGS   "-fltconsistency")
    ADD_FLAGS(CMAKE_Fortran_FLAGS   "-fp_port -mp")
    ADD_FLAGS(CMAKE_Fortran_FLAGS   "-nofor_main")
    ADD_FLAGS(CMAKE_Fortran_LDFLAGS "-nofor_main")
ENDIF(_match_ifort)

# Free Fortran Compiler tricks
# ----------------------------
STRING(REGEX MATCH ".*ftn$" _match_ftn ${CMAKE_Fortran_COMPILER})
IF(_match_ftn)
    MESSAGE(STATUS "Add -Mnomain to the Fortran linker")
    ADD_FLAGS(CMAKE_Fortran_FLAGS   "-Mnomain")
    ADD_FLAGS(CMAKE_Fortran_LDFLAGS "-Mnomain")
    ADD_FLAGS(CMAKE_Fortran_LDFLAGS "-Bstatic")
ENDIF(_match_ftn)

# IBM tricks
# ----------
STRING(REGEX MATCH ".*xlc$" _match_xlc ${CMAKE_C_COMPILER})
IF(_match_xlc)
    MESSAGE(ERROR "Please use the thread-safe version of the xlc compiler (xlc_r)")
ENDIF (_match_xlc)

STRING(REGEX MATCH ".*xlc_r$" _match_xlc ${CMAKE_C_COMPILER})
IF(_match_xlc)
    ADD_FLAGS(CMAKE_C_FLAGS "-qstrict")
    ADD_FLAGS(CMAKE_C_FLAGS "-qthreaded")
ENDIF(_match_xlc)

STRING(REGEX MATCH ".*xlf$" _match_xlf ${CMAKE_Fortran_COMPILER})
IF(_match_xlf)
    MESSAGE(ERROR "Please use the thread-safe version of the xlf compiler (xlf_r)")
ENDIF(_match_xlf)

STRING(REGEX MATCH ".*xlf_r$" _match_xlf ${CMAKE_Fortran_COMPILER})
IF(_match_xlf)
    MESSAGE(STATUS "Add -nofor_main to the Fortran linker")
    ADD_FLAGS(CMAKE_Fortran_FLAGS   "-qstrict")
    ADD_FLAGS(CMAKE_Fortran_FLAGS   "-qthreaded")
    ADD_FLAGS(CMAKE_Fortran_FLAGS   "-nofor_main")
    ADD_FLAGS(CMAKE_Fortran_LDFLAGS "-nofor_main")
ENDIF(_match_xlf)

# Fix the building system for 32 or 64 bits.
# ------------------------------------------
#
# On MAC OS X there is a easy solution, by setting the 
# CMAKE_OSX_ARCHITECTURES to a subset of the following values:
# ppc;ppc64;i386;x86_64.
# On Linux this is a little bit tricky. We have to check that the
# compiler supports the -m32/-m64 flags as well as the linker.
# Once this issue resolved the CMAKE_C_FLAGS and CMAKE_C_LDFLAGS
# have to be updated accordingly.
#
# TODO: Same trick for the Fortran compiler...
#       no idea how to correctly detect if the required/optional
#      libraries are in the correct format.
#
IF(BUILD_64bits)
    IF(_match_xlc OR _match_xlf)
        SET(ARCH_BUILD "-q64")
    ELSE(_match_xlc OR _match_xlf)
        SET(ARCH_BUILD "-m64")
    ENDIF(_match_xlc OR _match_xlf)
ELSE(BUILD_64bits)
    IF(_match_xlc OR _match_xlf)
        SET(ARCH_BUILD "-q32")
    ELSE(_match_xlc OR _match_xlf)
        SET(ARCH_BUILD "-m32")
    ENDIF(_match_xlc OR _match_xlf)
ENDIF(BUILD_64bits)

SET(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${ARCH_BUILD}")
CHECK_C_COMPILER_FLAG(${ARCH_BUILD} C_M32or64)

IF(C_M32or64)
    ADD_FLAGS(CMAKE_C_FLAGS         "${ARCH_BUILD}")
    ADD_FLAGS(CMAKE_C_LDFLAGS       "${ARCH_BUILD}")
    ADD_FLAGS(CMAKE_Fortran_FLAGS   "${ARCH_BUILD}")
    ADD_FLAGS(CMAKE_Fortran_LDFLAGS "${ARCH_BUILD}")
    MESSAGE(STATUS "Add ${ARCH_BUILD} to the C/Fortran compiler/linker")
ENDIF(C_M32or64)

# Set warnings for debug builds 
# -----------------------------
CHECK_C_COMPILER_FLAG("-Wall" HAVE_WALL)
IF(HAVE_WALL)
    ADD_FLAGS(C_WFLAGS "-Wall")
ENDIF(HAVE_WALL)
CHECK_C_COMPILER_FLAG("-Wextra" HAVE_WEXTRA)
IF(HAVE_WEXTRA)
    ADD_FLAGS(C_WFLAGS "-Wextra")
ENDIF(HAVE_WEXTRA)

# flags for the overly verbose icc
# --------------------------------
CHECK_C_COMPILER_FLAG("-wd424" HAVE_WD)
IF(HAVE_WD) 
    # 424: checks for duplicate ";"
    # 981: every volatile triggers a "unspecified evaluation order", obnoxious
    #      but might be useful for some debugging sessions. 
    # 1419: warning about extern functions being declared in .c
    #       files
    # 1572: cuda compares floats with 0.0f. 
    ADD_FLAGS(C_WFLAGS "-wd424")
    ADD_FLAGS(C_WFLAGS "-wd981")
    ADD_FLAGS(C_WFLAGS "-wd1419")
    ADD_FLAGS(C_WFLAGS "-wd1572")
ENDIF(HAVE_WD)
ADD_FLAGS(CMAKE_C_FLAGS_DEBUG "${C_WFLAGS}")
ADD_FLAGS(CMAKE_C_FLAGS_RELWITHDEBINFO "${C_WFLAGS}")

########################################
#                                      #
#     Looking for others libraries     #
#      stored in MAGMA_EXTRA_LIBS      #
#       thanks to DEFINE_LIBRARY       #
#                                      #
########################################

# Check if libgfortran is available
# ---------------------------------
IF(HAVE_GNU)
    DEFINE_LIBRARY("gfortran")
ENDIF(HAVE_GNU)

# Check if libifcore is available
# ---------------------------------
IF(HAVE_INTEL)
    DEFINE_LIBRARY("ifcore")
ENDIF(HAVE_INTEL)

# Check if libstdc++ is available
# -------------------------------
CHECK_C_COMPILER_FLAG("-lstdc++" HAVE_STDCPP)
IF(HAVE_STDCPP)
    DEFINE_LIBRARY("stdc++")
ENDIF(HAVE_STDCPP)

# Check if libdl is available 
# ---------------------------
CHECK_C_COMPILER_FLAG("-ldl" HAVE_DL)
IF(HAVE_DL)
    DEFINE_LIBRARY("dl")
ENDIF(HAVE_DL)

# Check if libm is available 
# ---------------------------
CHECK_C_COMPILER_FLAG("-lm" HAVE_M)
IF(HAVE_M)
    DEFINE_LIBRARY("m")
ENDIF(HAVE_M)

# Check if pthread is available
# -----------------------------
FIND_PACKAGE(Threads)
IF(Threads_FOUND)
    SET(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES};${CMAKE_THREAD_LIBS_INIT}")
    CHECK_FUNCTION_EXISTS(pthread_create HAVE_PTHREAD)
    IF(HAVE_PTHREAD)
        LIST(APPEND MAGMA_EXTRA_LIBS "${CMAKE_THREAD_LIBS_INIT}")
    ENDIF(HAVE_PTHREAD)
ENDIF(Threads_FOUND)

CHECK_FUNCTION_EXISTS(sched_setaffinity HAVE_SCHED_SETAFFINITY)
if( NOT HAVE_SCHED_SETAFFINITY )
  CHECK_LIBRARY_EXISTS(rt sched_setaffinity "" HAVE_SCHED_SETAFFINITY)
endif( NOT HAVE_SCHED_SETAFFINITY )

# Timeval, timespec, realtime clocks, etc
# ---------------------------------------
CHECK_STRUCT_HAS_MEMBER("struct timespec" tv_nsec time.h HAVE_TIMESPEC_TV_NSEC)
IF(NOT HAVE_TIMESPEC_TV_NSEC)
    ADD_DEFINITIONS(-D_GNU_SOURCE)
    CHECK_STRUCT_HAS_MEMBER("struct timespec" tv_nsec time.h HAVE_TIMESPEC_TV_NSEC)
ENDIF(NOT HAVE_TIMESPEC_TV_NSEC)

CHECK_LIBRARY_EXISTS(rt clock_gettime "" HAVE_CLOCK_GETTIME)
IF(HAVE_CLOCK_GETTIME)
    LIST(APPEND MAGMA_EXTRA_LIBS "rt")
ENDIF(HAVE_CLOCK_GETTIME)

########################################
#                                      #
#       Looking for some headers       #
#                                      #
########################################

# stdlib, stdio, string, getopt, etc
# ----------------------------------
CHECK_INCLUDE_FILES(stdarg.h HAVE_STDARG_H)

# va_copy is special as it is not required to be a function.
# ----------------------------------------------------------
IF(HAVE_STDARG_H)
    CHECK_C_SOURCE_COMPILES("
      #include <stdarg.h>
      int main(void) {
        va_list a, b;
        va_copy(a, b);
        return 0;
      }"
      HAVE_VA_COPY
      )
 
    IF(NOT HAVE_VA_COPY)
        CHECK_C_SOURCE_COMPILES("
      #include <stdarg.h>
      int main(void) {
        va_list a, b;
        __va_copy(a, b);
        return 0;
      }"
        HAVE_UNDERSCORE_VA_COPY
        )
    ENDIF(NOT HAVE_VA_COPY)
ENDIF(HAVE_STDARG_H)

CHECK_FUNCTION_EXISTS(asprintf HAVE_ASPRINTF)
CHECK_FUNCTION_EXISTS(vasprintf HAVE_VASPRINTF)
CHECK_FUNCTION_EXISTS(getopt_long HAVE_GETOPT_LONG)
CHECK_FUNCTION_EXISTS(getrusage HAVE_GETRUSAGE)
CHECK_INCLUDE_FILES(errno.h HAVE_ERRNO_H)
CHECK_INCLUDE_FILES(stddef.h HAVE_STDDEF_H)
CHECK_INCLUDE_FILES(getopt.h HAVE_GETOPT_H)
CHECK_INCLUDE_FILES(unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES(limits.h HAVE_LIMITS_H)
CHECK_INCLUDE_FILES(string.h HAVE_STRING_H)

##
## @end file SystemDetection.cmake
##
