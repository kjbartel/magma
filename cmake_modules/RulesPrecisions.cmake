# 
# MAGMA Internal: generation of various floating point precision files from a template.
#

set(GENDEPENDENCIES  ${CMAKE_SOURCE_DIR}/tools/genDependencies.py)
set(PRECISIONPP      ${CMAKE_SOURCE_DIR}/tools/codegen.py)
set(PRECISIONPP_subs ${CMAKE_SOURCE_DIR}/tools/magmasubs.py)

include(ParseArguments)
FIND_PACKAGE(PythonInterp REQUIRED)
#
# Generates a rule for every SOURCES file, to create the precisions in PRECISIONS. If TARGETDIR
# is not empty then all generated files will be prepended with the $TARGETDIR/.
# A new file is created, from a copy by default
# If the first precision is "/", all occurences of the basename in the file are remplaced by 
# "pbasename" where p is the selected precision. 
# the target receives a -DPRECISION_p in its cflags. 
#
macro(precisions_rules_py)
  PARSE_ARGUMENTS(PREC_RULE
    "TARGETDIR;PRECISIONS"
    ""
    ${ARGN})

  MESSAGE(STATUS "Generate precision dependencies in ${CMAKE_CURRENT_SOURCE_DIR}")

  # The first is the output variable list
  CAR(OUTPUTLIST ${PREC_RULE_DEFAULT_ARGS})
  # Everything else should be source files.
  CDR(SOURCES ${PREC_RULE_DEFAULT_ARGS})
  # By default the TARGETDIR is the current binary directory
  if( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
    set(PREC_RULE_TARGETDIR "./")
    set(PRECISIONPP_prefix "./")
    set(PRECISIONPP_arg "-P")
  else( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    else(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    endif(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    set(PRECISIONPP_arg "-P")
    set(PRECISIONPP_prefix "${PREC_RULE_TARGETDIR}")
  endif( "${PREC_RULE_TARGETDIR}" STREQUAL "" )

  #set(PRECISIONPP_prefix "${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR}")

  set(options_list "")
  foreach(prec_rules_PREC ${PREC_RULE_PRECISIONS})
    set(options_list "${options_list} ${prec_rules_PREC}")
  endforeach()

  set(sources_list "")
  foreach(_src ${SOURCES})
    set(sources_list "${sources_list} ${_src}")
  endforeach()

  #message(STATUS ${prec_options})
  #message(STATUS ${sources_list})

  set(gencmd ${PYTHON_EXECUTABLE} ${GENDEPENDENCIES} -f "${sources_list}" -p "${options_list}" -s "${CMAKE_CURRENT_SOURCE_DIR}" ${PRECISIONPP_arg} ${PRECISIONPP_prefix})
  execute_process(COMMAND ${gencmd} OUTPUT_VARIABLE dependencies_list)  

  #message(STATUS "List: ${dependencies_list}")
  
  foreach(_dependency ${dependencies_list})
    
    string(STRIP "${_dependency}" _dependency)
    string(COMPARE NOTEQUAL "${_dependency}" "" not_empty)
    if( not_empty )
      #message(STATUS "dependency: ${_dependency}")

      string(REGEX REPLACE "^(.*),(.*),(.*)$" "\\1" _dependency_INPUT "${_dependency}")
      set(_dependency_PREC   "${CMAKE_MATCH_2}")
      set(_dependency_OUTPUT "${CMAKE_MATCH_3}")
      
      #message(STATUS "input :${_dependency_INPUT}")
      #message(STATUS "prec  :${_dependency_PREC}")
      #message(STATUS "output:${_dependency_OUTPUT}")
      
      set(pythoncmd ${PYTHON_EXECUTABLE} ${PRECISIONPP} -f ${CMAKE_CURRENT_SOURCE_DIR}/${_dependency_INPUT} -p ${_dependency_PREC} ${PRECISIONPP_arg} ${PRECISIONPP_prefix})
      
      string(STRIP "${_dependency_OUTPUT}" _dependency_OUTPUT)
      string(COMPARE NOTEQUAL "${_dependency_OUTPUT}" "" got_file)
      string(COMPARE NOTEQUAL "${_dependency_OUTPUT}" "${_dependency_INPUT}" generate_out )
      
      # Force the copy of the original files in the binary_dir
      # for VPATH compilation
      if( NOT MAGMA_COMPILE_INPLACE )
        set(generate_out 1)
      endif()
      
      # We generate a dependency only if a file will be generated
      if( got_file )
        if( generate_out )
          # the custom command is executed in CMAKE_CURRENT_BINARY_DIR
          add_custom_command(
            OUTPUT ${_dependency_OUTPUT}
            COMMAND ${CMAKE_COMMAND} -E remove -f ${_dependency_OUTPUT} && ${pythoncmd} && chmod a-w ${_dependency_OUTPUT}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_dependency_INPUT} ${PRECISIONPP} ${PRECISIONPP_subs})
          
          set(_dependency_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${_dependency_OUTPUT}")
          set_source_files_properties(${_dependency_OUTPUT} PROPERTIES GENERATED 1 IS_IN_BINARY_DIR 1 ) 
          #"COMPILE_FLAGS -DPRECISION_${_dependency_PREC}" 
          
        else( generate_out )
          set_source_files_properties(${_dependency_OUTPUT} PROPERTIES GENERATED 0 )
          #"COMPILE_FLAGS -DPRECISION_${_dependency_PREC}" 
        endif( generate_out )
        
        #message(STATUS "OUTPUT: ${_dependency_OUTPUT}")
        list(APPEND ${OUTPUTLIST} ${_dependency_OUTPUT})
      endif( got_file )
    endif()
  endforeach()

  MESSAGE(STATUS "Generate precision dependencies in ${CMAKE_CURRENT_SOURCE_DIR} - Done")

endmacro(precisions_rules_py)


macro(precisions_rules_py_backup)
  PARSE_ARGUMENTS(PREC_RULE
    "TARGETDIR;PRECISIONS"
    ""
    ${ARGN})

  MESSAGE(STATUS "Generate precision dependencies in ${CMAKE_CURRENT_SOURCE_DIR}")

  # The first is the output variable list
  CAR(OUTPUTLIST ${PREC_RULE_DEFAULT_ARGS})
  # Everything else should be source files.
  CDR(SOURCES ${PREC_RULE_DEFAULT_ARGS})
  # By default the TARGETDIR is the current binary directory
  if( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
    set(PREC_RULE_TARGETDIR "./")
    set(PRECISIONPP_prefix "./")
    set(PRECISIONPP_arg "-P")
  else( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    else(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    endif(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    set(PRECISIONPP_arg "-P")
    set(PRECISIONPP_prefix "${PREC_RULE_TARGETDIR}")
  endif( "${PREC_RULE_TARGETDIR}" STREQUAL "" )

  #set(PRECISIONPP_prefix "${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR}")

  foreach(prec_rules_SOURCE ${SOURCES})
    foreach(prec_rules_PREC ${PREC_RULE_PRECISIONS})
      
      set(pythoncmd ${PYTHON_EXECUTABLE} ${PRECISIONPP} -f ${CMAKE_CURRENT_SOURCE_DIR}/${prec_rules_SOURCE} -p ${prec_rules_PREC} ${PRECISIONPP_arg} ${PRECISIONPP_prefix})
      execute_process(COMMAND ${pythoncmd} --out
        OUTPUT_VARIABLE prec_rules_OSRC)

      #message(STATUS "${pythoncmd}")

      string(STRIP "${prec_rules_OSRC}" prec_rules_OSRC)
      string(COMPARE NOTEQUAL "${prec_rules_OSRC}" "" got_file)
      string(COMPARE NOTEQUAL "${prec_rules_OSRC}" "${prec_rules_SOURCE}" generate_out )

      # Force the copy of the original files in the binary_dir
      # for VPATH compilation
      if( NOT MAGMA_COMPILE_INPLACE )
        set(generate_out 1)
      endif()

      #message(STATUS "OUTPUT: ${prec_rules_OSRC}")

      # We generate a dependency only if a file will be generated
      if( got_file )

        if( generate_out )
          # the custom command is executed in CMAKE_CURRENT_BINARY_DIR
          add_custom_command(
            OUTPUT ${prec_rules_OSRC}
            COMMAND ${CMAKE_COMMAND} -E remove -f ${prec_rules_OSRC} && ${pythoncmd} && chmod a-w ${prec_rules_OSRC}
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${prec_rules_SOURCE} ${PRECISIONPP} ${PRECISIONPP_subs})

          set(prec_rules_OSRC "${CMAKE_CURRENT_BINARY_DIR}/${prec_rules_OSRC}")
          set_source_files_properties(${prec_rules_OSRC} PROPERTIES GENERATED 1 IS_IN_BINARY_DIR 1 ) 
          #"COMPILE_FLAGS -DPRECISION_${prec_rules_PREC}" 

        else( generate_out )
          set_source_files_properties(${prec_rules_OSRC} PROPERTIES GENERATED 0 )
          #"COMPILE_FLAGS -DPRECISION_${prec_rules_PREC}" 
        endif( generate_out )

        #message(STATUS "OUTPUT: ${prec_rules_OSRC}")
        list(APPEND ${OUTPUTLIST} ${prec_rules_OSRC})
      endif( got_file )
    endforeach()
  endforeach()

  MESSAGE(STATUS "Generate precision dependencies in ${CMAKE_CURRENT_SOURCE_DIR} - Done")

endmacro(precisions_rules_py_backup)


