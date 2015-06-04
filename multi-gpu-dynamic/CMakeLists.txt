###
#
# @file          : CMakeLists.txt
#
# @description   : Project MORSE
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 24-01-2012
# @last modified : mar. 29 mai 2012 10:25:55 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(VersionMORSE)

#############################################
#                                           #
#        Compilation of the MORSE           #
#                                           #
#############################################

# Include headers
# ---------------
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/common)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/common)

# Check for the subdirectories 
# ----------------------------
ADD_SUBDIRECTORY(include)
ADD_SUBDIRECTORY(common)
ADD_SUBDIRECTORY(compute)
ADD_SUBDIRECTORY(morse)
ADD_SUBDIRECTORY(timing)

# Create file morse_starpu.pc
# ---------------------------
IF(MORSE_SCHED_STARPU)
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/morse-starpu.pc.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/lib/pkgconfig/morse-starpu.pc"
                   @ONLY)
    INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib/pkgconfig/morse-starpu.pc
            DESTINATION lib/pkgconfig)
ENDIF()

# Create file morse_quark.pc
# --------------------------
IF(MORSE_SCHED_QUARK)
    CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/lib/pkgconfig/morse-quark.pc.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/lib/pkgconfig/morse-quark.pc"
                   @ONLY)
    INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib/pkgconfig/morse-quark.pc
            DESTINATION lib/pkgconfig)
ENDIF(MORSE_SCHED_QUARK)

# Generate documentation
# ----------------------
INCLUDE(share/info)

###
### END CMakeLists.txt
###
