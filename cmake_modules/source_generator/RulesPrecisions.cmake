###
#
#  @file RulesPrecisions.cmake 
# 
#  @project MAGMA
#  MAGMA is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
# 
#  @version 0.1.0
#  @author MAGMA developers
#  @date 13-07-2012
#   
###
# 
# MORSE Internal: generation of various floating point precision files from a template.
#
###


INCLUDE(ParseArguments)
FIND_PACKAGE(PythonInterp REQUIRED)
#
# Generates a rule for every SOURCES file, to create the precisions in PRECISIONS. If TARGETDIR
# is not empty then all generated files will be prepended with the $TARGETDIR/.
# A new file is created, from a copy by default
# If the first precision is "/", all occurences of the basename in the file are remplaced by 
# "pbasename" where p is the selected precision. 
# the target receives a -DPRECISION_p in its cflags. 
#
MACRO(precisions_rules_py)
  PARSE_ARGUMENTS(PREC_RULE "TARGETDIR;PRECISIONS;DICTIONARY" "" ${ARGN})

  IF("^${PREC_RULE_DICTIONARY}$" STREQUAL "^MAGMA$")
    SET(GENDEPENDENCIES  ${CMAKE_SOURCE_DIR}/tools/genDependencies.py)
    SET(PRECISIONPP      ${CMAKE_SOURCE_DIR}/tools/codegen.py)
    SET(PRECISIONPP_subs ${CMAKE_SOURCE_DIR}/tools/magmasubs.py)
    MESSAGE(STATUS "Generate precision dependencies with magmasubs in ${CMAKE_CURRENT_SOURCE_DIR}")

  ELSEIF("^${PREC_RULE_DICTIONARY}$" STREQUAL "^MORSE$")
    SET(GENDEPENDENCIES  ${CMAKE_SOURCE_DIR}/tools/precision-generator/genDependencies.py)
    SET(PRECISIONPP      ${CMAKE_SOURCE_DIR}/tools/precision-generator/codegen.py)
    SET(PRECISIONPP_subs ${CMAKE_SOURCE_DIR}/tools/precision-generator/morsesubs.py)
    MESSAGE(STATUS "Generate precision dependencies with morsesubs in ${CMAKE_CURRENT_SOURCE_DIR}")

  ELSE()
    SET(GENDEPENDENCIES  ${CMAKE_SOURCE_DIR}/tools/genDependencies.py)
    SET(PRECISIONPP      ${CMAKE_SOURCE_DIR}/tools/codegen.py)
    SET(PRECISIONPP_subs ${CMAKE_SOURCE_DIR}/tools/magmasubs.py)
    MESSAGE(STATUS "Generate precision dependencies with magmasubs in ${CMAKE_CURRENT_SOURCE_DIR}")

  ENDIF()


  # The first is the output variable list
  CAR(OUTPUTLIST ${PREC_RULE_DEFAULT_ARGS})
  # Everything else should be source files.
  CDR(SOURCES ${PREC_RULE_DEFAULT_ARGS})
  # By default the TARGETDIR is the current binary directory
  IF( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
    SET(PREC_RULE_TARGETDIR "./")
    SET(PRECISIONPP_prefix "./")
    SET(PRECISIONPP_arg "-P")
  ELSE( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
    IF(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    ELSE(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
      FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    ENDIF(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    SET(PRECISIONPP_arg "-P")
    SET(PRECISIONPP_prefix "${PREC_RULE_TARGETDIR}")
  ENDIF( "${PREC_RULE_TARGETDIR}" STREQUAL "" )

  #set(PRECISIONPP_prefix "${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR}")

  SET(options_list "")
  FOREACH(prec_rules_PREC ${PREC_RULE_PRECISIONS})
    SET(options_list "${options_list} ${prec_rules_PREC}")
  ENDFOREACH()

  SET(sources_list "")
  FOREACH(_src ${SOURCES})
    SET(sources_list "${sources_list} ${_src}")
  ENDFOREACH()

  #MESSAGE(STATUS ${prec_options})
  #MESSAGE(STATUS ${sources_list})

  SET(gencmd ${PYTHON_EXECUTABLE} ${GENDEPENDENCIES} -f "${sources_list}" -p "${options_list}" -s "${CMAKE_CURRENT_SOURCE_DIR}" ${PRECISIONPP_arg} ${PRECISIONPP_prefix})
  EXECUTE_PROCESS(COMMAND ${gencmd} OUTPUT_VARIABLE dependencies_list)  

  #MESSAGE(STATUS "List: ${dependencies_list}")
  
  FOREACH(_dependency ${dependencies_list})
    
    STRING(STRIP "${_dependency}" _dependency)
    STRING(COMPARE NOTEQUAL "${_dependency}" "" not_empty)
    IF( not_empty )
      #MESSAGE(STATUS "dependency: ${_dependency}")

      STRING(REGEX REPLACE "^(.*),(.*),(.*)$" "\\1" _dependency_INPUT "${_dependency}")
      SET(_dependency_PREC   "${CMAKE_MATCH_2}")
      SET(_dependency_OUTPUT "${CMAKE_MATCH_3}")
      
      #MESSAGE(STATUS "input :${_dependency_INPUT}")
      #MESSAGE(STATUS "prec  :${_dependency_PREC}")
      #MESSAGE(STATUS "output:${_dependency_OUTPUT}")
      
      SET(pythoncmd ${PYTHON_EXECUTABLE} ${PRECISIONPP} -f ${CMAKE_CURRENT_SOURCE_DIR}/${_dependency_INPUT} -p ${_dependency_PREC} ${PRECISIONPP_arg} ${PRECISIONPP_prefix})
      
      STRING(STRIP "${_dependency_OUTPUT}" _dependency_OUTPUT)
      STRING(COMPARE NOTEQUAL "${_dependency_OUTPUT}" "" got_file)
      STRING(COMPARE NOTEQUAL "${_dependency_OUTPUT}" "${_dependency_INPUT}" generate_out )
      
      # Force the copy of the original files in the binary_dir
      # for VPATH compilation
      IF( NOT MORSE_COMPILE_INPLACE )
        SET(generate_out 1)
      ENDIF()
      
      # We generate a dependency only if a file will be generated
      IF( got_file )
        IF( generate_out )
          # the custom command is executed in CMAKE_CURRENT_BINARY_DIR
          ADD_CUSTOM_COMMAND(
            OUTPUT ${_dependency_OUTPUT}
            COMMAND ${CMAKE_COMMAND} -E remove -f ${_dependency_OUTPUT} && ${pythoncmd} && chmod a-w ${_dependency_OUTPUT}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_dependency_INPUT} ${PRECISIONPP} ${PRECISIONPP_subs})
          
          SET(_dependency_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${_dependency_OUTPUT}")
          SET_SOURCE_FILES_PROPERTIES(${_dependency_OUTPUT} PROPERTIES GENERATED 1 IS_IN_BINARY_DIR 1 ) 
          #"COMPILE_FLAGS -DPRECISION_${_dependency_PREC}" 
          
        ELSE( generate_out )
          SET_SOURCE_FILES_PROPERTIES(${_dependency_OUTPUT} PROPERTIES GENERATED 0 )
          #"COMPILE_FLAGS -DPRECISION_${_dependency_PREC}" 
        ENDIF( generate_out )
        
        #MESSAGE(STATUS "OUTPUT: ${_dependency_OUTPUT}")
        LIST(APPEND ${OUTPUTLIST} ${_dependency_OUTPUT})
      ENDIF( got_file )
    ENDIF()
  ENDFOREACH()

  IF("^${PREC_RULE_DICTIONARY}$" STREQUAL "^MAGMA$")
    MESSAGE(STATUS "Generate precision dependencies with magmasubs in ${CMAKE_CURRENT_SOURCE_DIR} - Done")
  ELSEIF("^${PREC_RULE_DICTIONARY}$" STREQUAL "^MORSE$")
    MESSAGE(STATUS "Generate precision dependencies with morsesubs in ${CMAKE_CURRENT_SOURCE_DIR} - Done")
  ELSE()
    MESSAGE(STATUS "Generate precision dependencies with magmasubs in ${CMAKE_CURRENT_SOURCE_DIR} - Done")
  ENDIF()

ENDMACRO(precisions_rules_py)


#macro(precisions_rules_py_backup)
#  PARSE_ARGUMENTS(PREC_RULE
#    "TARGETDIR;PRECISIONS"
#    ""
#    ${ARGN})
#
#  MESSAGE(STATUS "Generate precision dependencies in ${CMAKE_CURRENT_SOURCE_DIR}")
#
#  # The first is the output variable list
#  CAR(OUTPUTLIST ${PREC_RULE_DEFAULT_ARGS})
#  # Everything else should be source files.
#  CDR(SOURCES ${PREC_RULE_DEFAULT_ARGS})
#  # By default the TARGETDIR is the current binary directory
#  if( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
#    set(PREC_RULE_TARGETDIR "./")
#    set(PRECISIONPP_prefix "./")
#    set(PRECISIONPP_arg "-P")
#  else( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
#    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
#    else(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
#      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
#    endif(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
#    set(PRECISIONPP_arg "-P")
#    set(PRECISIONPP_prefix "${PREC_RULE_TARGETDIR}")
#  endif( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
#
#  #set(PRECISIONPP_prefix "${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR}")
#
#  foreach(prec_rules_SOURCE ${SOURCES})
#    foreach(prec_rules_PREC ${PREC_RULE_PRECISIONS})
#      
#      set(pythoncmd ${PYTHON_EXECUTABLE} ${PRECISIONPP} -f ${CMAKE_CURRENT_SOURCE_DIR}/${prec_rules_SOURCE} -p ${prec_rules_PREC} ${PRECISIONPP_arg} ${PRECISIONPP_prefix})
#      execute_process(COMMAND ${pythoncmd} --out
#        OUTPUT_VARIABLE prec_rules_OSRC)
#
#      #message(STATUS "${pythoncmd}")
#
#      string(STRIP "${prec_rules_OSRC}" prec_rules_OSRC)
#      string(COMPARE NOTEQUAL "${prec_rules_OSRC}" "" got_file)
#      string(COMPARE NOTEQUAL "${prec_rules_OSRC}" "${prec_rules_SOURCE}" generate_out )
#
#      # Force the copy of the original files in the binary_dir
#      # for VPATH compilation
#      if( NOT MAGMA_COMPILE_INPLACE )
#        set(generate_out 1)
#      endif()
#
#      #message(STATUS "OUTPUT: ${prec_rules_OSRC}")
#
#      # We generate a dependency only if a file will be generated
#      if( got_file )
#
#        if( generate_out )
#          # the custom command is executed in CMAKE_CURRENT_BINARY_DIR
#          add_custom_command(
#            OUTPUT ${prec_rules_OSRC}
#            COMMAND ${CMAKE_COMMAND} -E remove -f ${prec_rules_OSRC} && ${pythoncmd} && chmod a-w ${prec_rules_OSRC}
#            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${prec_rules_SOURCE} ${PRECISIONPP} ${PRECISIONPP_subs})
#
#          set(prec_rules_OSRC "${CMAKE_CURRENT_BINARY_DIR}/${prec_rules_OSRC}")
#          set_source_files_properties(${prec_rules_OSRC} PROPERTIES GENERATED 1 IS_IN_BINARY_DIR 1 ) 
#          #"COMPILE_FLAGS -DPRECISION_${prec_rules_PREC}" 
#
#        else( generate_out )
#          set_source_files_properties(${prec_rules_OSRC} PROPERTIES GENERATED 0 )
#          #"COMPILE_FLAGS -DPRECISION_${prec_rules_PREC}" 
#        endif( generate_out )
#
#        #message(STATUS "OUTPUT: ${prec_rules_OSRC}")
#        list(APPEND ${OUTPUTLIST} ${prec_rules_OSRC})
#      endif( got_file )
#    endforeach()
#  endforeach()
#
#  MESSAGE(STATUS "Generate precision dependencies in ${CMAKE_CURRENT_SOURCE_DIR} - Done")
#
#endmacro(precisions_rules_py_backup)

##
## @end file RulesPrecisions.cmake
##
