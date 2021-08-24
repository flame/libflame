####################################################################
#  LAPACK make include file.                                       #
#  LAPACK, Version 3.8.0                                           #
#  November 2017                                                   #
####################################################################

SHELL = /bin/sh

#  CC is the C compiler, normally invoked with options CFLAGS.
#
CC = gcc

#CFLAGS = -O0 -g
CFLAGS = -O3

#  GTEST Related build variables
GTESTDIR=/home/amd/nagendra/softwares/GoogleTest/googletest-master/Gtest_Install
GTESTLIBDIR=$(GTESTDIR)/lib64
GTESTINC=$(GTESTDIR)/include  
GTESTLIB=$(GTESTLIBDIR)/libgtest.a $(GTESTLIBDIR)/libgtest_main.a

#  GBENCH Related build variables
GBENCHDIR=/home/amd/nagendra/softwares/Google_Bench/benchmark-main/Gbench_Install
GBENCHLIBDIR=$(GBENCHDIR)/lib64
GBENCHINC=$(GBENCHDIR)/include
GBENCHLIB=$(GBENCHLIBDIR)/libbenchmark.a $(GBENCHLIBDIR)/libbenchmark_main.a

CXX = g++
#CXXFLAGS = -O3  -Wall -Wno-comment -fopenmp -I$(GTESTINC) -I$(GBENCHINC) -I../include/
CXXFLAGS = -O0 -g  -Wall -Wno-comment -fopenmp -I$(GTESTINC) -I$(GBENCHINC) -I../include/

#  Modify the FC and FFLAGS definitions to the desired compiler
#  and desired compiler options for your machine.  NOOPT refers to
#  the compiler options desired when NO OPTIMIZATION is selected.
#
#  Note: During a regular execution, LAPACK might create NaN and Inf
#  and handle these quantities appropriately. As a consequence, one
#  should not compile LAPACK with flags such as -ffpe-trap=overflow.
#
FC = gfortran
#FCLAGS = -O0 -g  -frecursive
#FFLAGS = -O0 -g  -frecursive
FCLAGS = -O3  -frecursive
FFLAGS = -O3  -frecursive
FFLAGS_DRV = $(FFLAGS)
FFLAGS_NOOPT = -O0 -frecursive

#  Define LDFLAGS to the desired linker options for your machine.
#
LDFLAGS =  -fopenmp

#  The archiver and the flag(s) to use when building an archive
#  (library).  If your system has no ranlib, set RANLIB = echo.
#
AR = ar
ARFLAGS = cr
RANLIB = ranlib

#  Timer for the SECOND and DSECND routines
#
#  Default:  SECOND and DSECND will use a call to the
#  EXTERNAL FUNCTION ETIME
#TIMER = EXT_ETIME
#  For RS6K:  SECOND and DSECND will use a call to the
#  EXTERNAL FUNCTION ETIME_
#TIMER = EXT_ETIME_
#  For gfortran compiler:  SECOND and DSECND will use a call to the
#  INTERNAL FUNCTION ETIME
TIMER = INT_ETIME
#  If your Fortran compiler does not provide etime (like Nag Fortran
#  Compiler, etc...) SECOND and DSECND will use a call to the
#  INTERNAL FUNCTION CPU_TIME
#TIMER = INT_CPU_TIME
#  If none of these work, you can use the NONE value.
#  In that case, SECOND and DSECND will always return 0.
#TIMER = NONE

#  Uncomment the following line to include deprecated routines in
#  the LAPACK library.
#
#BUILD_DEPRECATED = Yes

#  Location of the extended-precision BLAS (XBLAS) Fortran library
#  used for building and testing extended-precision routines.  The
#  relevant routines will be compiled and XBLAS will be linked only
#  if USEXBLAS is defined.
#
#USEXBLAS = Yes
#XBLASLIB = -lxblas

#  The location of the libraries to which you will link.  (The
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#

#BLASLIB      = $(TOPSRCDIR)/librefblas.a
BLASLIB      = /home/amd/nagendra/work/LibFlame/LibFlame_perf_test/blis/Blis_GCC_Install/lib/libblis.a

#CBLASLIB     = $(TOPSRCDIR)/libcblas.a
LAPACKLIB    =/home/amd/nagendra/work/LibFlame/LibFlame_perf_test/LF_3p1/libflame/libflame_gcc_install/lib/libflame.a


#TMGLIB       = $(TOPSRCDIR)/libtmglib.a
#TMGLIB       = /home/nagendra/Netlib_ref_libs/lapack-3.9.0/libtmglib.a


