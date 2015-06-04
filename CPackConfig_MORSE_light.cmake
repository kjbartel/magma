###
#
# @file          : CPackConfig.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 22-03-2012
# @last modified : mer. 09 mai 2012 16:07:04 CEST
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

SET(SUFFIX_VERSION "")

EXECUTE_PROCESS(COMMAND pwd
                OUTPUT_VARIABLE SOURCE_DIRECTORY
                OUTPUT_STRIP_TRAILING_WHITESPACE
                )
EXECUTE_PROCESS(COMMAND ls ${SOURCE_DIRECTORY}/externals OUTPUT_VARIABLE LIST_FILE_TO_REMOVE OUTPUT_STRIP_TRAILING_WHITESPACE)
STRING(REPLACE "\n" ";" LIST_FILE_TO_REMOVE "${LIST_FILE_TO_REMOVE}")
FOREACH(_file ${LIST_FILE_TO_REMOVE})
    EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove -f ${SOURCE_DIRECTORY}/externals/${_file})
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
SET(CPACK_PACKAGE_FILE_NAME           "magma-morse-light-${MORSE_VERSION}${SUFFIX_VERSION}")
LIST(APPEND CPACK_IGNORE_FILES        "/\\.svn/;\\.swp$;\\.#;/#;~$;.DS_Store;\\.old;\\.tar$;\\.tar.gz$;\\.tgz$")
LIST(APPEND CPACK_IGNORE_FILES        "/exp/;/build/")
LIST(APPEND CPACK_IGNORE_FILES        "CPackConfig_MAGMA.cmake;CPackConfig_MORSE_full.cmake;CPackConfig_MORSE_light.cmake")
LIST(APPEND CPACK_IGNORE_FILES        "magma-1gpu-${MORSE_VERSION}${SUFFIX_VERSION}.tar.gz")
LIST(APPEND CPACK_IGNORE_FILES        "magma-morse-full-${MORSE_VERSION}${SUFFIX_VERSION}.tar.gz")
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
