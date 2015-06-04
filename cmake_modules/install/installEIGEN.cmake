###
#
# @file      : installEIGEN.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 21-01-2012
# @last modified : mar. 12 juin 2012 16:44:46 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(installPACKAGE)
INCLUDE(downloadPACKAGE)

MACRO(INSTALL_EIGEN _MODE)

    # Define extension of the library
    # -------------------------------
    IF(APPLE)
        SET(EIGEN_EXTENSION dylib)
    ELSE(APPLE)
        SET(EIGEN_EXTENSION so)
    ENDIF(APPLE)

    # Define prefix paths
    # -------------------
    IF(MORSE_SEPARATE_PROJECTS)
        SET(BLAS_PATH ${CMAKE_INSTALL_PREFIX}/eigen)
    ELSE(MORSE_SEPARATE_PROJECTS)
        SET(BLAS_PATH ${CMAKE_INSTALL_PREFIX})
    ENDIF(MORSE_SEPARATE_PROJECTS)

    # Define steps of installation
    # ----------------------------
    SET(BLAS_CONFIG_CMD cd build && ${CMAKE_COMMAND} ../)
    SET(BLAS_MAKE_CMD cd build && ${CMAKE_MAKE_PROGRAM} blas lapack)
    SET(BLAS_MAKEINSTALL_CMD ${CMAKE_COMMAND} -E echo)

    # Define additional step
    # ----------------------
    SET(BLAS_ADD_STEP eigen_create_build_out_source eigen_create_prefix eigen_install_blas eigen_install_lapack)
    SET(eigen_create_build_out_source_CMD ${CMAKE_COMMAND} -E make_directory build)
    SET(eigen_create_build_out_source_DIR ${CMAKE_BINARY_DIR}/externals/blas)
    SET(eigen_create_build_out_source_DEP_BEFORE download)
    SET(eigen_create_build_out_source_DEP_AFTER configure)
    SET(eigen_create_prefix_CMD ${CMAKE_COMMAND} -E make_directory ${BLAS_PATH}/lib)
    SET(eigen_create_prefix_DIR ${CMAKE_INSTALL_PREFIX})
    SET(eigen_create_prefix_DEP_BEFORE build)
    SET(eigen_create_prefix_DEP_AFTER install)
    SET(eigen_install_blas_CMD ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/externals/blas/build/blas/libeigen_blas.${EIGEN_EXTENSION} ${BLAS_PATH}/lib/libeigen_blas.${EIGEN_EXTENSION})
    SET(eigen_install_blas_DIR ${CMAKE_INSTALL_PREFIX})
    SET(eigen_install_blas_DEP_BEFORE eigen_create_prefix)
    SET(eigen_install_blas_DEP_AFTER install)
    SET(eigen_install_lapack_CMD ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/externals/blas/build/lapack/libeigen_lapack.${EIGEN_EXTENSION} ${BLAS_PATH}/lib/libeigen_lapack.${EIGEN_EXTENSION})
    SET(eigen_install_lapack_DIR ${CMAKE_INSTALL_PREFIX})
    SET(eigen_install_lapack_DEP_BEFORE eigen_create_prefix)
    SET(eigen_install_lapack_DEP_AFTER install)

    # Define options
    # --------------
    SET(BLAS_OPTIONS -DCMAKE_INSTALL_PREFIX=${BLAS_PATH})
    SET(BLAS_OPTIONS ${BLAS_OPTIONS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
    SET(BLAS_OPTIONS ${BLAS_OPTIONS} -DEIGEN_BUILD_PKGCONFIG=OFF)
    
    # Install the external package
    # ----------------------------
    DEFINE_DOWNLOAD_PACKAGE("blas" "${_MODE}")
    INSTALL_PACKAGE("blas" "${BLAS_BUILD_MODE}")
    MESSAGE(STATUS "Installing BLAS - eigen version")

    # Set linker flags
    # ----------------
    SET(BLAS_LIBRARY_PATH ${BLAS_PATH}/lib)
    SET(BLAS_LDFLAGS "-L${BLAS_LIBRARY_PATH} -leigen_blas")
    SET(BLAS_LIBRARIES "eigen_blas")
    SET(EIGENLAPACK_LIBRARY_PATH ${BLAS_PATH}/lib)
    SET(EIGENLAPACK_LDFLAGS "-L${BLAS_LIBRARY_PATH} -leigen_lapack")
    SET(EIGENLAPACK_LIBRARIES "eigen_lapack")


ENDMACRO(INSTALL_EIGEN)

###
### END installEIGEN.cmake
###
