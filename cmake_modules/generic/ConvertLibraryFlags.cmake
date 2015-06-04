###
#
# @file          : ConvertLibraryFlags.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 09-03-2012
# @last modified : ven. 09 mars 2012 09:52:26 CET
#
###
#- convert from -L and -l flags to full paths of libraries
#  convert_library_flags(<variable> [flags ...])
#
# Given some compiler flags, replace the -L... and -l... flags with full
# paths to libraries, and store the result into <variable>. This is
# useful for converting the output from scripts like pkg-config.
#
###

function(convert_library_flags variable)

  # grab libdirs from the -L flags
  set(libdirs)
  foreach(flag ${ARGN})
    if(flag MATCHES "^-L")
      # chop -L and append to libdirs
      string(REGEX REPLACE "^-L(.*)$" "\\1" dir ${flag})
      list(APPEND libdirs ${dir})
    endif()
  endforeach(flag)

  # now convert flags to result
  set(result)
  foreach(flag ${ARGN})
    if(flag MATCHES "^-L")
      # do nothing, removes -L flags from result
    elseif(flag MATCHES "^-l")
      # chop -l
      string(REGEX REPLACE "^-l(.*)$" "\\1" lib ${flag})

      # We cannot use find_library, because we do not want a cache
      # variable. So do the search manually. This uses three nested
      # foreach loops (for dir, prefix, suffix).
      #
      # FIXME - Where does the compiler look for libraries, when
      # there is no -L flag? Assuming /usr/lib and /usr/local/lib
      # but this is wrong and nonportable.
      #
      # FIXME - This fails to find shared libraries in OpenBSD,
      # because of no "*.so" symlinks without a version number.
      #
      set(go TRUE)
      foreach(dir ${libdirs} /usr/lib /usr/local/lib)
        foreach(prefix ${CMAKE_FIND_LIBRARY_PREFIXES})
          foreach(suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
            if(go)
              set(file ${dir}/${prefix}${lib}${suffix})
              if(EXISTS ${file})
                # found it! append to result
                list(APPEND result ${file})
                set(go FALSE) # break from nested loops
              endif()
            endif()
          endforeach(suffix)
        endforeach(prefix)
      endforeach(dir)

      if(go)
        message(SEND_ERROR "library for flag ${flag}: not found")
      endif(go)

    else()
      # Flag is not -L or -l, might be something like -pthread, so
      # just preserve it.
      list(APPEND result ${flag})
    endif()
  endforeach(flag)

  # return the result
  set("${variable}" ${result} PARENT_SCOPE)

endfunction(convert_library_flags)

###
### END ConvertLibraryFlags.cmake
###
