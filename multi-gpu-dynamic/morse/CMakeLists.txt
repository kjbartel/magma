###
#
# @file          : CMakeLists.txt
#
# @description   : Project MORSE
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 25-01-2012
# @last modified : mer. 25 janv. 2012 21:15:36 CET
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# Check for the subdirectories 
# ----------------------------
IF(MORSE_SCHED_QUARK)
    ADD_SUBDIRECTORY(quark)
ENDIF()

IF(MORSE_SCHED_STARPU)
    ADD_SUBDIRECTORY(starpu)
ENDIF()

###
### END CMakeLists.txt
###
