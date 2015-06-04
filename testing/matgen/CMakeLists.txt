###
#
# @file      : CMakeLists.txt
#
# @description   : Project MAGMA
#
# @version       :
# @created by    : 
# @creation date : 23-01-2012
# @last modIFied : lun. 23 janv. 2012 19:44:56 CET
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# List of sources files
# ---------------------

SET(SCATGEN  slatm1.f slaran.f slarnd.f)

SET(SMATGEN  slatms.f slatme.f slatmr.f slatmt.f 
   slagge.f slagsy.f slakf2.f slarge.f slaror.f slarot.f slatm2.f 
   slatm3.f slatm5.f slatm6.f slatm7.f slahilb.f)

SET(CMATGEN  clatms.f clatme.f clatmr.f clatmt.f 
   clagge.f claghe.f clagsy.f clakf2.f clarge.f claror.f clarot.f 
   clatm1.f clarnd.f clatm2.f clatm3.f clatm5.f clatm6.f clahilb.f)

SET(DZATGEN  dlatm1.f dlaran.f dlarnd.f)

SET(DMATGEN  dlatms.f dlatme.f dlatmr.f dlatmt.f 
   dlagge.f dlagsy.f dlakf2.f dlarge.f dlaror.f dlarot.f dlatm2.f 
   dlatm3.f dlatm5.f dlatm6.f dlatm7.f dlahilb.f)

SET(ZMATGEN  zlatms.f zlatme.f zlatmr.f zlatmt.f 
  zlagge.f zlaghe.f zlagsy.f zlakf2.f zlarge.f zlaror.f zlarot.f 
  zlatm1.f zlarnd.f zlatm2.f zlatm3.f zlatm5.f zlatm6.f zlahilb.f)

SET(ALLOBJ ${SMATGEN}
           ${CMATGEN}
           ${SCATGEN}
           ${DMATGEN}
           ${ZMATGEN}
           ${DZATGEN}
           )

IF(BUILD_SINGLE)
  SET(ALLOBJ ${SMATGEN} ${SCATGEN})
ENDIF()

IF(BUILD_DOUBLE)
  SET(ALLOBJ ${DMATGEN} ${DZATGEN})
ENDIF()

IF(BUILD_COMPLEX)
  SET(ALLOBJ  ${CMATGEN} ${SCATGEN})
ENDIF()

IF(BUILD_COMPLEX16)
  SET(ALLOBJ ${ZMATGEN} ${DZATGEN})
ENDIF()

# Create libtmglib.a
# ------------------
ADD_LIBRARY(tmglib ${ALLOBJ})
SET_PROPERTY(TARGET tmglib PROPERTY LINKER_LANGUAGE Fortran)
TARGET_LINK_LIBRARIES(tmglib ${EXTRA_LIBS})
FOREACH(_external_package "blas" "lapack" "plasma")
    STRING(TOUPPER "${_external_package}" _NAMEVAR)
    IF(DEFINED ${_NAMEVAR}_BUILD_MODE)
        ADD_DEPENDENCIES(tmglib ${_external_package}_install)
    ENDIF()
ENDFOREACH()

# Install step of libtmglib.a
# ---------------------------
INSTALL(TARGETS tmglib
    DESTINATION lib)

###
### END CMakeLists.txt
###
