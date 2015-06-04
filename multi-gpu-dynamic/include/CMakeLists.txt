###
#
# @file          : CMakeLists.txt
#
# @description   : Project MORSE
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 24-01-2012
# @last modified : lun. 19 mars 2012 17:30:53 CET
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(RulesPrecisions)

# Generate the morse headers for all possible precisions
# ------------------------------------------------------
SET(MORSE_HEADERS_GENERATED "")
SET(ZHDR 
    magma_morse_z.h
    morse_z.h
)

precisions_rules_py(MORSE_HEADERS_GENERATED
            "${ZHDR}"
            PRECISIONS "z;c;d;s;zc;ds")

# Define the list of headers
# --------------------------
SET(MORSE_HEADERS 
    magma_morse.h
    morse.h
    morse_kernels.h
    ${MORSE_HEADERS_GENERATED}
    )

# Force generation of headers
# ---------------------------
ADD_CUSTOM_TARGET(morse_include ALL SOURCES ${MORSE_HEADERS})

# installation
# ------------
INSTALL(FILES ${MORSE_HEADERS} 
        DESTINATION include)

###
### END CMakeLists.txt
###
