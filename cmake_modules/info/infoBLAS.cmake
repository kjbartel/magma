###
#
# @file          : infoBLAS.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-04-2012
# @last modified : mar. 12 juin 2012 16:31:54 CEST
#
###

MACRO(BLAS_INFO_INSTALL)

    # Get the implementation of BLAS
    # ------------------------------
    STRING(TOUPPER "${MORSE_USE_BLAS}" VALUE_MORSE_USE_BLAS)
    IF("${VALUE_MORSE_USE_BLAS}"MATCHES "EIGEN")
        INCLUDE(infoEIGEN)
        EIGEN_INFO_INSTALL()

    ELSEIF("${VALUE_MORSE_USE_BLAS}" MATCHES "REFBLAS")
        INCLUDE(infoREFBLAS)
        REFBLAS_INFO_INSTALL()

    ELSE()
        SET(VALUE_MORSE_USE_BLAS "REFBLAS")
        INCLUDE(infoREFBLAS)
        REFBLAS_INFO_INSTALL()

    ENDIF()


    # Define web link of blas
    # -----------------------
    IF(NOT DEFINED BLAS_URL)
        SET(BLAS_URL    ${${VALUE_MORSE_USE_BLAS}_URL})
    ENDIF()

    # Define tarball of blas
    # ----------------------
    IF(NOT DEFINED BLAS_TARBALL)
        SET(BLAS_TARBALL ${${VALUE_MORSE_USE_BLAS}_TARBALL})
    ENDIF()

    # Define md5sum of blas
    # ---------------------
    IF(DEFINED BLAS_URL OR DEFINED BLAS_TARBALL)
        SET(BLAS_MD5SUM  ${${VALUE_MORSE_USE_BLAS}_MD5SUM})
    ENDIF()

    # Define repository of blas
    # -------------------------
    IF(NOT DEFINED BLAS_SVN_REP)
        SET(BLAS_REPO_MODE ${${VALUE_MORSE_USE_BLAS}_REPO_MODE})
        SET(BLAS_SVN_REP   ${${VALUE_MORSE_USE_BLAS}_SVN_REP}  )
        SET(BLAS_SVN_ID    ${${VALUE_MORSE_USE_BLAS}_SVN_ID}   )
        SET(BLAS_SVN_PWD   ${${VALUE_MORSE_USE_BLAS}_SVN_PWD}  )
    ENDIF()

   # Define dependencies
   # -------------------
    SET(BLAS_DEPENDENCIES ${${VALUE_MORSE_USE_BLAS}_DEPENDENCIES})

ENDMACRO(BLAS_INFO_INSTALL)

MACRO(BLAS_INFO_FIND)

    # Get the implementation of BLAS
    # ------------------------------
    STRING(TOUPPER "${MORSE_USE_BLAS}" VALUE_MORSE_USE_BLAS)
    IF("${VALUE_MORSE_USE_BLAS}"MATCHES "EIGEN")
        INCLUDE(infoEIGEN)
        EIGEN_INFO_FIND()

    ELSEIF("${VALUE_MORSE_USE_BLAS}" MATCHES "REFBLAS")
        INCLUDE(infoREFBLAS)
        REFBLAS_INFO_FIND()

    ELSE()
        SET(VALUE_MORSE_USE_BLAS "REFBLAS")
        INCLUDE(infoREFBLAS)
        REFBLAS_INFO_FIND()

    ENDIF()

    # Define parameters for FIND_MY_PACKAGE
    # -------------------------------------
    SET(BLAS_name_library        "${${VALUE_MORSE_USE_BLAS}_name_library}"       )
    SET(BLAS_name_pkgconfig      "${${VALUE_MORSE_USE_BLAS}_name_pkgconfig}"     )
    SET(BLAS_name_include        "${${VALUE_MORSE_USE_BLAS}_name_include}"       )
    SET(BLAS_name_include_suffix "${${VALUE_MORSE_USE_BLAS}_name_include_suffix}")
    SET(BLAS_name_fct_test       "${${VALUE_MORSE_USE_BLAS}_name_fct_test}"      )
    SET(BLAS_name_binary         "${${VALUE_MORSE_USE_BLAS}_name_binary}"        )

ENDMACRO(BLAS_INFO_FIND)

###
### END infoBLAS.cmake
###
