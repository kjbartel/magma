###
#
# @file          : FortranMangling.cmake
#
# @description   : Project MORSE (http://hiepacs.bordeaux.inria.fr/eagullo/morse)
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 21-03-2012
# @last modified : mer. 21 mars 2012 15:26:07 CET
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

MACRO(FortranMangling)

    MESSAGE(STATUS "Looking Fortran mangling")
    FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/test_mangle)
    FILE(WRITE ${CMAKE_BINARY_DIR}/test_mangle/test.f "
       program intface
       external c_intface
       integer i
       call c_intface(i)
       stop
       end
       \n")
    FILE(WRITE ${CMAKE_BINARY_DIR}/test_mangle/test.c "
       #include <stdio.h>
       void c_intface_(int *i){fprintf(stdout, \"-DADD_\");fflush(stdout);}
       void c_intface(int *i){fprintf(stdout, \"-DNOCHANGE\");fflush(stdout);}
       void c_intface__(int *i){fprintf(stdout, \"-DfcIsF2C\");fflush(stdout);}
       void C_INTFACE(int *i){fprintf(stdout, \"-DUPCASE\");fflush(stdout);}
       \n")
    EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -c test.c
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test_mangle
                    OUTPUT_QUIET
                    ERROR_QUIET
                   )
    EXECUTE_PROCESS(COMMAND ${CMAKE_Fortran_COMPILER} test.f test.o -o test_mangle
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test_mangle
                    OUTPUT_QUIET
                    ERROR_QUIET
                   )
    EXECUTE_PROCESS(COMMAND ./test_mangle
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test_mangle
                    RESULT_VARIABLE FORTRAN_MANGLING_RESULT
                    OUTPUT_VARIABLE FORTRAN_MANGLING_DETECTED
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                   )
    IF(NOT FORTRAN_MANGLING_RESULT)
        ADD_DEFINITIONS(${FORTRAN_MANGLING_DETECTED})
        MESSAGE(STATUS "Looking Fortran mangling - Done")
    ELSE()
        MESSAGE(FATAL_ERROR "Looking Fortran mangling - error")
    ENDIF()

ENDMACRO(FortranMangling)

###
### END FortranMangling.cmake
###
