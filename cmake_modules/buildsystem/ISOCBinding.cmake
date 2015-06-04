###
#
#  @file ISOCBinding.cmake
#
#  @project MAGMA
#  MAGMA is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.2.0
#  @author Cedric Castagnede
#  @date 13-07-2012
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

MACRO(ISOCBinding)

    MESSAGE(STATUS "Looking Fortran ISO C Binding")
    FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test_iso_c_binding)
    FILE(WRITE ${CMAKE_BINARY_DIR}/test_iso_c_binding/test.f "
       program binding
       use iso_c_binding
       end
       \n")
    EXECUTE_PROCESS(COMMAND ${CMAKE_Fortran_COMPILER} test.f -o test_binding
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test_iso_c_binding
                    RESULT_VARIABLE FORTRAN_BINDING_MISSING
                    OUTPUT_QUIET
                    ERROR_QUIET
                   )
    IF(FORTRAN_BINDING_MISSING)
        SET(HAVE_ISO_C_BINDING FALSE)
        MESSAGE(STATUS "Looking Fortran ISO C Binding - not found")
    ELSE()
        SET(HAVE_ISO_C_BINDING TRUE)
        MESSAGE(STATUS "Looking Fortran ISO C Binding - working")
    ENDIF()

ENDMACRO(ISOCBinding)

###
## @end file ISOCBinding.cmake
###
