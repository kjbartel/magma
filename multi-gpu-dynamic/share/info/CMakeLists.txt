###
#
# @file          : CMakeLists.txt
#
# @description   : Project MORSE
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 22-05-2012
# @last modified : mar. 29 mai 2012 10:35:54 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Create file version.texi
# ------------------------
INCLUDE(VersionMORSE)
CONFIGURE_FILE("version.texi.in"
               "version.texi"
               @ONLY)

# Looking for makeinfo
# --------------------
FIND_PROGRAM(MAKEINFO_COMPILER makeinfo)
IF(MAKEINFO_COMPILER)
    # Add target
    # ----------
    ADD_CUSTOM_TARGET(info ALL 
                      COMMAND ${MAKEINFO_COMPILER} morse.texi)
    ADD_CUSTOM_TARGET(html ALL
                      COMMAND ${MAKEINFO_COMPILER} --html --no-split morse.texi --css-include=morse.css)

    # Installation
    # ------------
    INSTALL(FILES morse.info
            DESTINATION share/info)
    INSTALL(FILES morse.html
            DESTINATION share/html)
ELSE(MAKEINFO_COMPILER)
    MESSAGE(STATUS "Looking for makeinfo - not found")
ENDIF(MAKEINFO_COMPILER)

# Looking for texi2dvi
# --------------------
FIND_PROGRAM(TEXI2DVI_COMPILER texi2dvi)
IF(TEXI2DVI_COMPILER)
    # Add target
    # ----------
    ADD_CUSTOM_TARGET(pdf ALL COMMAND ${TEXI2DVI_COMPILER} --pdf --batch morse.texi)

    # Installation
    # ------------
    INSTALL(FILES morse.pdf
            DESTINATION share/pdf)
ELSE(TEXI2DVI_COMPILER)
    MESSAGE(STATUS "Looking for texi2dvi - not found")
ENDIF(TEXI2DVI_COMPILER)

###
### END CMakeLists.txt
###
