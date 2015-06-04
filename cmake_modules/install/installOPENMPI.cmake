###
#
# @file      : installOPENMPI.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 21-01-2012
# @last modified : mer. 16 mai 2012 10:19:06 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
INCLUDE(installPACKAGE)
INCLUDE(downloadPACKAGE)
INCLUDE(infoOPENMPI)

MACRO(INSTALL_OPENMPI _MODE)

    # Get info for this package
    # -------------------------
    MPI_INFO_INSTALL()

    # Define prefix paths
    # -------------------
    IF(MORSE_SEPARATE_PROJECTS)
        SET(MPI_PATH ${CMAKE_INSTALL_PREFIX}/mpi)
    ELSE(MORSE_SEPARATE_PROJECTS)
        SET(MPI_PATH ${CMAKE_INSTALL_PREFIX})
    ENDIF(MORSE_SEPARATE_PROJECTS)

    # Define steps of installation
    # ----------------------------
    SET(MPI_CONFIG_CMD ./configure)
    SET(MPI_MAKE_CMD ${CMAKE_MAKE_PROGRAM})
    SET(MPI_MAKEINSTALL_CMD ${CMAKE_MAKE_PROGRAM} install)

    # Define options
    # --------------
    SET(MPI_OPTIONS --prefix=${MPI_PATH} --enable-mpi-thread-multiple)

    # Install the external package
    # ----------------------------
    # ----------------------------
    DEFINE_DOWNLOAD_PACKAGE("mpi" "${_MODE}")
    INSTALL_PACKAGE("mpi" "${MPI_BUILD_MODE}")

    # Set linker flags
    # ----------------
    SET(MPIEXEC               "${MPI_PATH}/bin/mpiexec")
    SET(MPIEXEC_MAX_NUMPROCS  "2"                      )
    SET(MPIEXEC_NUMPROC_FLAGS "-np"                    )
    SET(MPIEXEC_POSTFLAGS     ""                       )
    SET(MPIEXEC_PREFLAGS      ""                       )

    SET(MPI_CXX_COMPILER      "${MPI_PATH}/bin/mpicixx"                                )
    SET(MPI_CXX_COMPILE_FLAGS ""                                                       )
    SET(MPI_CXX_INCLUDE_PATH  "${MPI_PATH}/include"                                    )
    SET(MPI_CXX_LIBRARIES     "${MPI_PATH}/lib/libmpi_cxx.so;${MPI_PATH}/lib/libmpi.so")
    SET(MPI_CXX_LINK_FLAGS    "-Wl,--export-dynamic"                                   )

    SET(MPI_C_COMPILER      "${MPI_PATH}/bin/mpicc"    )
    SET(MPI_C_COMPILE_FLAGS ""                         )
    SET(MPI_C_INCLUDE_PATH  "${MPI_PATH}/include"      )
    SET(MPI_C_LIBRARIES     "${MPI_PATH}/lib/libmpi.so")
    SET(MPI_C_LINK_FLAGS    "-Wl,--export-dynamic"     )

    SET(MPI_Fortran_COMPILER      "${MPI_PATH}/bin/mpif90"                                                               )
    SET(MPI_Fortran_COMPILE_FLAGS ""                                                                                     )
    SET(MPI_Fortran_INCLUDE_PATH  "${MPI_PATH}/include"                                                                  )
    SET(MPI_Fortran_LIBRARIES     "${MPI_PATH}/lib/libmpi_f90.so;${MPI_PATH}/lib/libmpi_f77.so;${MPI_PATH}/lib/libmpi.so")
    SET(MPI_Fortran_LINK_FLAGS    "-Wl,--export-dynamic"                                                                 )

    SET(MPI_LIBRARY       "${MPI_PATH}lib/libmpi_cxx.so")
    SET(MPI_EXTRA_LIBRARY "${MPI_PATH}/lib/libmpi.so"   )

    SET(MPI_LDFLAGS      "-L${MPI_PATH}/lib -lmpi_cxx -lmpi -lmpi_f90")
    SET(MPI_LIBRARIES    "mpi_cxx;mpi;mpi_f90")
    SET(MPI_BINARY_PATH  ${MPI_PATH}/bin)
    SET(MPI_LIBRARY_PATH ${MPI_PATH}/lib)
    SET(MPI_INCLUDE_PATH ${MPI_PATH}/include)

ENDMACRO(INSTALL_OPENMPI)

###
### END installOPENMPI.cmake
###

