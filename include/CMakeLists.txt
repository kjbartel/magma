###
#
# @file      : CMakeLists.txt
#
# @description   : Project MAGMA
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 19-01-2012
# @last modified : mer. 25 janv. 2012 17:36:23 CET
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(RulesPrecisions)

# Generate the magma headers for all possible precisions
# ------------------------------------------------------
SET(MAGMA_HEADERS_GENERATED "")
SET(ZHDR magma_z.h
         magma_zc.h
         magmablas_z.h
         magmablas_zc.h
         magma_zlapack.h
         magma_zbulge.h
         magma_zbulgeinc.h
     )

precisions_rules_py(MAGMA_HEADERS_GENERATED "${ZHDR}"
                    PRECISIONS "${MAGMA_PRECISION}"
                    DICTIONARY "MAGMA")

# Define the list of headers
# --------------------------
SET(MAGMA_HEADERS 
    cblas.h
    magma.h
    magmablas.h
    magma_lapack.h
    ${MAGMA_HEADERS_GENERATED}
    )

# Force generation of headers
# ---------------------------
ADD_CUSTOM_TARGET(magma_include ALL SOURCES ${MAGMA_HEADERS})

# installation
# ------------
INSTALL(FILES ${MAGMA_HEADERS}
        DESTINATION include)

###
### END CMakeLists.txt
### 
