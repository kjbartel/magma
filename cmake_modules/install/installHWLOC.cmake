###
#
# @file      : installHWLOC.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : mer. 16 mai 2012 10:17:10 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(installPACKAGE)
INCLUDE(downloadPACKAGE)
INCLUDE(infoHWLOC)

MACRO(INSTALL_HWLOC _MODE)

    # Get info for this package
    # -------------------------
    HWLOC_INFO_INSTALL()

    # Define prefix paths
    # -------------------
    IF(MORSE_SEPARATE_PROJECTS)
        SET(HWLOC_PATH ${CMAKE_INSTALL_PREFIX}/hwloc)
    ELSE(MORSE_SEPARATE_PROJECTS)
        SET(HWLOC_PATH ${CMAKE_INSTALL_PREFIX})
    ENDIF(MORSE_SEPARATE_PROJECTS)

    # Define steps of installation
    # ----------------------------
    SET(HWLOC_CONFIG_CMD ./configure)
    SET(HWLOC_MAKE_CMD ${CMAKE_MAKE_PROGRAM})
    SET(HWLOC_MAKEINSTALL_CMD ${CMAKE_MAKE_PROGRAM} install)

    # Define options
    # --------------
    SET(HWLOC_OPTIONS --prefix=${HWLOC_PATH})

    # Install the external package
    # ----------------------------
    DEFINE_DOWNLOAD_PACKAGE("hwloc" "${_MODE}")
    INSTALL_PACKAGE("hwloc" "${HWLOC_BUILD_MODE}")

    # Set linker flags
    # ----------------
    SET(HWLOC_BINARY_PATH ${HWLOC_PATH}/bin)
    SET(HWLOC_LIBRARY_PATH ${HWLOC_PATH}/lib)
    SET(HWLOC_INCLUDE_PATH ${HWLOC_PATH}/include)
    SET(HWLOC_LDFLAGS "-L${HWLOC_LIBRARY_PATH} -lhwloc")
    SET(HWLOC_LIBRARIES "hwloc")

ENDMACRO(INSTALL_HWLOC)

###
### END installHWLOC.cmake
###

