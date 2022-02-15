####################################################################
#  LAPACKE TEST-SUITE make include file.                           #
#  LAPACKE, Version 3.10.0                                         #
#  DECEMBER 2021                                                   #
####################################################################

SHELL = /bin/sh

#  GTEST Related build variables
GTESTDIR=
GTESTLIBDIR=$(GTESTDIR)/lib
GTESTINC=$(GTESTDIR)/include  
GTESTLIB=$(GTESTLIBDIR)/libgtest.a $(GTESTLIBDIR)/libgtest_main.a

#  The location of the libraries to which you will link.  (The
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
LIBFLAME_INSTALL_DIR =
BLIS_INSTALL_DIR =

#Netlib reference libraries location 
NETLIB_INSTALL_DIR=

BLASLIB      =$(BLIS_INSTALL_DIR)/lib/libblis-mt.a
LAPACKLIB    =$(LIBFLAME_INSTALL_DIR)/lib/libflame.a
LAPACKELIB   =$(LIBFLAME_INSTALL_DIR)/lib/liblapacke.a

NETLIB_LAPACKE_LIB="$(NETLIB_INSTALL_DIR)/liblapacke.so.3.10.0"
NETLIB_LAPACK_LIB="$(NETLIB_INSTALL_DIR)/liblapack.so.3.10.0"
NETLIB_BLAS_LIB="$(NETLIB_INSTALL_DIR)/libblas.so.3.10.0"

#  CC is the C compiler, normally invoked with options CFLAGS.
#
CC = gcc
#CC = clang

CFLAGS = -O3 -DREF_LPKE_LIB=$(NETLIB_LAPACKE_LIB) -DREF_BLAS_LIB=$(NETLIB_BLAS_LIB)

CXX = g++
#CXX = clang++

EXTFLAGS = -DREF_LPKE_LIB=$(NETLIB_LAPACKE_LIB) -DREF_BLAS_LIB=$(NETLIB_BLAS_LIB) 

CXXFLAGS = -O3  -DREF_LPKE_LIB=$(NETLIB_LAPACKE_LIB) -DREF_BLAS_LIB=$(NETLIB_BLAS_LIB)  -Wno-comment -fopenmp -I$(GTESTINC) -I$(LIBFLAME_INSTALL_DIR)/include/

#  Modify the FC and FFLAGS definitions to the desired compiler
#  and desired compiler options for your machine.  NOOPT refers to
#  the compiler options desired when NO OPTIMIZATION is selected.
#
#  Note: During a regular execution, LAPACK might create NaN and Inf
#  and handle these quantities appropriately. As a consequence, one
#  should not compile LAPACK with flags such as -ffpe-trap=overflow.
#
FC = gfortran
#FC = flang

FCLAGS = -O3  -frecursive
FFLAGS = -O3  -frecursive
FFLAGS_DRV = $(FFLAGS)
FFLAGS_NOOPT = -O0 -frecursive

#  Define LDFLAGS to the desired linker options for your machine.
#
LDFLAGS =  -fopenmp -ldl

#  The archiver and the flag(s) to use when building an archive
#  (library).  If your system has no ranlib, set RANLIB = echo.
#
AR = ar
ARFLAGS = cr
RANLIB = ranlib

#  Uncomment the following line to include deprecated routines in
#  the LAPACK library.
#
#BUILD_DEPRECATED = Yes

#  LAPACKE has the interface to some routines from tmglib.
#  If LAPACKE_WITH_TMG is defined, add those routines to LAPACKE.
#
#LAPACKE_WITH_TMG = Yes

#  Location of the extended-precision BLAS (XBLAS) Fortran library
#  used for building and testing extended-precision routines.  The
#  relevant routines will be compiled and XBLAS will be linked only
#  if USEXBLAS is defined.
#
#USEXBLAS = Yes
#XBLASLIB = -lxblas
