###
#
# @file      : installFXT.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : mer. 16 mai 2012 10:16:27 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(installPACKAGE)
INCLUDE(downloadPACKAGE)
INCLUDE(infoFXT)

MACRO(INSTALL_FXT _MODE)

    # Get info for this package
    # -------------------------
    FXT_INFO_INSTALL()

    # Define prefix paths
    # -------------------
    IF(MORSE_SEPARATE_PROJECTS)
        SET(FXT_PATH ${CMAKE_INSTALL_PREFIX}/fxt)
    ELSE(MORSE_SEPARATE_PROJECTS)
        SET(FXT_PATH ${CMAKE_INSTALL_PREFIX})
    ENDIF(MORSE_SEPARATE_PROJECTS)

    # Define steps of installation
    # ----------------------------
    SET(FXT_CONFIG_CMD ./configure)
    SET(FXT_MAKE_CMD ${CMAKE_MAKE_PROGRAM})
    SET(FXT_MAKEINSTALL_CMD ${CMAKE_MAKE_PROGRAM} install)

    # Define options
    # --------------
    SET(FXT_OPTIONS --prefix=${FXT_PATH})

    # Install the external package
    # -----------------------------
    DEFINE_DOWNLOAD_PACKAGE("fxt" "${_MODE}")
    INSTALL_PACKAGE("fxt" "${FXT_BUILD_MODE}")

    # Set linker flags
    # ----------------
    SET(FXT_BINARY_PATH ${FXT_PATH}/bin)
    SET(FXT_LIBRARY_PATH ${FXT_PATH}/lib)
    SET(FXT_INCLUDE_PATH ${FXT_PATH}/include/fxt)
    SET(FXT_LDFLAGS "-L${FXT_LIBRARY_PATH} -lfxt")
    SET(FXT_LIBRARIES "fxt")

ENDMACRO(INSTALL_FXT)

###
### END installFXT.cmake
###

