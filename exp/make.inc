#//////////////////////////////////////////////////////////////////////////////
#   -- MAGMA (version 1.1) --
#      Univ. of Tennessee, Knoxville
#      Univ. of California, Berkeley
#      Univ. of	Colorado, Denver
#      November 2011
#//////////////////////////////////////////////////////////////////////////////

#
# GPU_TARGET specifies for which GPU you want to compile MAGMA
#      0: Tesla family
#      1: Fermi Family
#
GPU_TARGET = 1

CC        = gcc
NVCC      = nvcc
FORT      = gfortran

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

OPTS      = -O3 -DADD_
FOPTS     = -O3 -DADD_ -x f95-cpp-input
NVOPTS    = --compiler-options -fno-strict-aliasing -DUNIX -O3 -DADD_
LDOPTS    = -fPIC -Xlinker -zmuldefs

LIB       = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -lpthread -lcublas -lm 
#LIB       = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core  -liomp5 -lpthread -lcublas -lm

CUDADIR   = /sw/keeneland/cuda/4.0/linux_binary

LIBDIR    = -L/opt/intel/composerxe-2011.1.107/mkl/lib/intel64 \
            -L/opt/intel/Compiler/11.1/038/lib/intel64 \
            -L$(CUDADIR)/lib64 
INC       = -I$(CUDADIR)/include 

LIBMAGMA     = ../lib/libmagma.a
LIBMAGMABLAS = ../../lib/libmagmablas.a




