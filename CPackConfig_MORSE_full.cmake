###
#
# @file          : CPackConfig.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 22-03-2012
# @last modified : mar. 12 juin 2012 15:40:21 CEST
#
###
# This file will be configured to contain variables for CPack. These variables
# should be set in the CMake list file of the project before CPack module is
# included. Example variables are:
#   CPACK_GENERATOR                     - Generator used to create package
#   CPACK_INSTALL_CMAKE_PROJECTS        - For each project (path, name, component)
#   CPACK_CMAKE_GENERATOR               - CMake Generator used for the projects
#   CPACK_INSTALL_COMMANDS              - Extra commands to install components
#   CPACK_INSTALL_DIRECTORIES           - Extra directories to install
#   CPACK_PACKAGE_DESCRIPTION_FILE      - Description file for the package
#   CPACK_PACKAGE_DESCRIPTION_SUMMARY   - Summary of the package
#   CPACK_PACKAGE_EXECUTABLES           - List of pairs of executables and labels
#   CPACK_PACKAGE_FILE_NAME             - Name of the package generated
#   CPACK_PACKAGE_ICON                  - Icon used for the package
#   CPACK_PACKAGE_INSTALL_DIRECTORY     - Name of directory for the installer
#   CPACK_PACKAGE_NAME                  - Package project name
#   CPACK_PACKAGE_VENDOR                - Package project vendor
#   CPACK_PACKAGE_VERSION               - Package project version
#   CPACK_PACKAGE_VERSION_MAJOR         - Package project version (major)
#   CPACK_PACKAGE_VERSION_MINOR         - Package project version (minor)
#   CPACK_PACKAGE_VERSION_PATCH         - Package project version (patch)
###
#
# cpack --config CPackConfig_MORSE_full.cmake
#
###

# Allow to add a suffix frr the package version 
SET(SUFFIX_VERSION "")

# Get the top directory of the project
EXECUTE_PROCESS(COMMAND pwd
                OUTPUT_VARIABLE SOURCE_DIRECTORY
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )

# The following lines can be modify by hydra to generate the package
#SET(REFBLAS_URL    "")
#SET(CBLAS_URL      "")
#SET(LAPACK_URL     "")
#SET(LAPACKE_URL    "")
#SET(PLASMA_URL     "")
#SET(MPI_URL        "")
#SET(FXT_URL        "")
#SET(HWLOC_URL      "")
#SET(STARPU_URL     "")
#SET(EIGEN_URL      "")

# Include info of package dependencies
SET(CMAKE_MODULE_PATH ${SOURCE_DIRECTORY}/cmake_modules/info)
include(infoREFBLAS)
REFBLAS_INFO_INSTALL()
include(infoCBLAS)
CBLAS_INFO_INSTALL()
include(infoLAPACK)
LAPACK_INFO_INSTALL()
include(infoLAPACKE)
LAPACKE_INFO_INSTALL()
include(infoPLASMA)
PLASMA_INFO_INSTALL()
include(infoOPENMPI)
MPI_INFO_INSTALL()
include(infoFXT)
FXT_INFO_INSTALL()
include(infoHWLOC)
HWLOC_INFO_INSTALL()
include(infoSTARPU)
STARPU_INFO_INSTALL()
include(infoEIGEN)
EIGEN_INFO_INSTALL()

# Remove all stuff in externals
EXECUTE_PROCESS(COMMAND ls ${SOURCE_DIRECTORY}/externals OUTPUT_VARIABLE LIST_FILE_TO_REMOVE OUTPUT_STRIP_TRAILING_WHITESPACE)
STRING(REPLACE "\n" ";" LIST_FILE_TO_REMOVE "${LIST_FILE_TO_REMOVE}")
FOREACH(_file ${LIST_FILE_TO_REMOVE})
    IF(IS_DIRECTORY ${SOURCE_DIRECTORY}/externals/${_file})
        EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove_directory ${SOURCE_DIRECTORY}/externals/${_file})
    ELSE()
        EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove -f ${SOURCE_DIRECTORY}/externals/${_file})
    ENDIF()
ENDFOREACH()

# Write a file in externals to force cpack to add this folder
FILE(WRITE ${SOURCE_DIRECTORY}/externals/info "Put all tarballs in this directory.\n")

# Get and put all packages in externals
FOREACH(_prefix REFBLAS CBLAS LAPACK LAPACKE PLASMA MPI FXT HWLOC STARPU EIGEN)
    IF(EXISTS ${${_prefix}_URL})
        EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${${_prefix}_URL} ${SOURCE_DIRECTORY}/externals/${${_prefix}_TARBALL})
        EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E md5sum ${SOURCE_DIRECTORY}/externals/${${_prefix}_TARBALL}
                        OUTPUT_VARIABLE ${_prefix}_MD5SUM_FOUND
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                       )
        IF(NOT "${${_prefix}_MD5SUM_FOUND}" MATCHES "${${_prefix}_MD5SUM}")
            MESSAGE(FATAL_ERROR "${_prefix} tarball -- not copied")
        ENDIF()
    ELSE(EXISTS ${${_prefix}_URL})
        FILE(DOWNLOAD ${${_prefix}_URL}
                      ${SOURCE_DIRECTORY}/externals/${${_prefix}_TARBALL}
                      EXPECTED_MD5 ${${_prefix}_MD5SUM}
                      STATUS IS_GOT
                      SHOW_PROGRESS
            )
        LIST(GET IS_GOT 0 IS_GOT_CODE)
        LIST(GET IS_GOT 1 IS_GOT_MSG)
        IF(IS_GOT_CODE)
            MESSAGE(FATAL_ERROR "${_prefix} tarball -- not downloaded")
        ENDIF(IS_GOT_CODE) 
    ENDIF(EXISTS ${${_prefix}_URL})
ENDFOREACH()

# Load MORSE default options
EXECUTE_PROCESS(COMMAND sed -i "s/INCLUDE(ConfigMAGMA)/INCLUDE(ConfigMORSE)/" ${SOURCE_DIRECTORY}/CMakeLists.txt)
MESSAGE(WARNING "
Run this following command to revert to original CMakeLists.txt:
  > sed -i \"s/INCLUDE(ConfigMORSE)/INCLUDE(ConfigMAGMA)/\" ${SOURCE_DIRECTORY}/CMakeLists.txt")

# Generate PACKAGE
SET(CPACK_GENERATOR                   "TGZ")
SET(CPACK_PACKAGE_NAME                "MORSE")
SET(CPACK_PACKAGE_VERSION_MAJOR       "1")
SET(CPACK_PACKAGE_VERSION_MINOR       "2")
SET(CPACK_PACKAGE_VERSION_PATCH       "0")
SET(MORSE_VERSION                     "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_PACKAGE_DESCRIPTION_FILE    "${SOURCE_DIRECTORY}/COPYRIGHT")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MORSE - Matrices Over Runtime Systems @ Exascale")
SET(CPACK_INSTALLED_DIRECTORIES       "${SOURCE_DIRECTORY};/")
SET(CPACK_INSTALL_CMAKE_PROJECTS      "")
SET(CPACK_PACKAGE_FILE_NAME           "magma-morse-full-${MORSE_VERSION}${SUFFIX_VERSION}")
LIST(APPEND CPACK_IGNORE_FILES        "/\\.svn/;\\.swp$;\\.#;/#;~$;.DS_Store;\\.old")
LIST(APPEND CPACK_IGNORE_FILES        "/exp/;/build/")
LIST(APPEND CPACK_IGNORE_FILES        "CPackConfig_MAGMA.cmake;CPackConfig_MORSE_full.cmake;CPackConfig_MORSE_light.cmake")
LIST(APPEND CPACK_IGNORE_FILES        "magma-1gpu-${MORSE_VERSION}${SUFFIX_VERSION}.tar.gz")
LIST(APPEND CPACK_IGNORE_FILES        "magma-morse-light-${MORSE_VERSION}${SUFFIX_VERSION}.tar.gz")
SET(CPACK_PACKAGE_VENDOR
"University of Tennessee,
Inria Bordeaux - Sud-Ouest,
University Paris-Sud,
King Abdullah University of Science and Technology,
University of California,
Berkeley, University of Colorado Denver")

###
### END CPackConfig.cmake
###
