###
#
# @file          : ConfigMORSE.cmake
#
# @description   :
#
# @version       :
# @created by    : Cedric Castagnede
# @creation date : 04-05-2012
# @last modified : lun. 11 juin 2012 15:51:34 CEST
#
###

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

OPTION(MAGMA_MORSE   
  "Enable/Disable MORSE version of MAGMA" ON)
OPTION(MAGMA_1GPU
  "Enable/Disable compilation of MAGMA mono-GPU" OFF)
#OPTION(MAGMA_MGPUS_STATIC
#  "Enable/Disable compilation of MAGMA multi-GPUs static" OFF)

SET(MORSE_USE_CUDA ""
   CACHE STRING "Enable/Disable CUDA dependency (ON/OFF/<not-defined>)")

###
### END ConfigMORSE.cmake
###
