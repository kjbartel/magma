###
#
# @file          : CMakeLists.txt
#
# @description   : Project MORSE
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 24-01-2012
# @last modified : mer. 25 janv. 2012 21:53:49 CET
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(RulesPrecisions)

# Generate the morse headers for all possible precisions
# ------------------------------------------------------
SET(MORSE_HEADERS_GENERATED "")
SET(ZHDR 
    compute_z.h
)

precisions_rules_py(MORSE_HEADERS_GENERATED
            "${ZHDR}"
            PRECISIONS ${MAGMA_PRECISION})

# Define the list of headers
# --------------------------
SET(MORSE_HEADERS 
    common.h
    context.h
    async.h
    ${MORSE_HEADERS_GENERATED}
    )

# Force generation of headers
# ---------------------------
ADD_CUSTOM_TARGET(morse_common_include ALL SOURCES ${MORSE_HEADERS})

# installation
# ------------
INSTALL(FILES ${MORSE_HEADERS}
        DESTINATION include)

###
### END CMakeLists.txt
###
