###
#
# @file          : installBLAS.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mar. 12 juin 2012 17:13:36 CEST
#
###

MACRO(INSTALL_BLAS _MODE)

    # Get the information of BLAS
    # ---------------------------
    INCLUDE(infoBLAS)
    BLAS_INFO_INSTALL()

    # Start the installation of BLAS
    # ------------------------------
    STRING(TOUPPER "${MORSE_USE_BLAS}" VALUE_MORSE_USE_BLAS)
    IF("${VALUE_MORSE_USE_BLAS}"MATCHES "EIGEN")
        INCLUDE(installEIGEN)
        INSTALL_EIGEN("${_MODE}")

    ELSEIF("${VALUE_MORSE_USE_BLAS}" MATCHES "REFBLAS")
        INCLUDE(installREFBLAS)
        INSTALL_REFBLAS("${_MODE}")

    ELSE()
        SET(VALUE_MORSE_USE_BLAS "REFBLAS")
        INCLUDE(installREFBLAS)
        INSTALL_REFBLAS("${_MODE}")

    ENDIF()

ENDMACRO(INSTALL_BLAS)

###
### END installBLAS.cmake
###
