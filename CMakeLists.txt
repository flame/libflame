##Copyright (C) 2020, Advanced Micro Devices, Inc.##

cmake_minimum_required(VERSION 3.0.0)

project(AOCL-LibFLAME-Win)
#change the dependency path to blis library path
SET(CMAKE_EXT_BLIS_LIBRARY_DEPENDENCY_PATH "${CMAKE_SOURCE_DIR}/bin/Release" CACHE STRING "AOCL-BLIS library path name")
SET(EXT_BLIS_LIBNAME "AOCL-LibBlis-Win-MT-dll.lib" CACHE STRING "AOCL-BLIS library name")

#gccopenmp
SET(OpenMP_libomp_LIBRARY "C:/Program Files/LLVM/lib/libomp.lib" CACHE STRING "openmp library path")
#Logic to create AOCL-libFLAME library name according to the linking AOCL-BLIS library name.
string(REPLACE "LibBlis" "LibFlame" PROJECT_NAME ${EXT_BLIS_LIBNAME})
string(REPLACE ".lib" "" PROJECT_NAME ${PROJECT_NAME})

#set the output directory path
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")


#enabling the configuration options
option (ENABLE_AMD_OPT "Enable AMD optimized path" OFF)
option (ENABLE_AMD_FLAGS "Enable AMD FLAGS internally enables ENABLE_BLAS_EXT_GEMMT, ENABLE_AMD_OPT, ENABLE_BUILTIN_LAPACK2FLAME, ENABLE_EXT_LAPACK_INTERFACE, ENABLE_F2C_DOTC, ENABLE_VOID_RETURN_COMPLEX_FUNCTION, ENABLE_MULTITHREADING" ON)
option (ENABLE_ILP64 "Enable ILP64" OFF)
option (ENABLE_UPPERCASE "Enable uppercase API's" OFF)
option (ENABLE_WRAPPER "Enable wrapper code" ON)
option (INCLUDE_LAPACKE "INCLUDE LAPACKE Code" ON)
option (BUILD_SHARED_LIBS "Building shared Library" ON)
option (BUILD_TEST "Building test directory" OFF)
option (BUILD_NETLIB_TEST "Building test directory" OFF)
option (ENABLE_NON_CRITICAL_CODE "Enable non-critical" ON)
option (ENABLE_BUILTIN_LAPACK2FLAME "lapack2flame library" OFF)
option (ENABLE_BUILTIN_BLAS "Built-in blas library" OFF)
option (ENABLE_BLIS1_USE_OF_FLA_MALLOC "Enable use of FLA-MALLOC" ON)
option (ENABLE_EXT_LAPACK_SUBPROBLEMS "Enable external lapack for subproblems" OFF)
option (ENABLE_EXT_LAPACK_INTERFACE "Enable external lapack interface" OFF)
option (ENABLE_BLAS3_FRNTEND_CNTL_TREES "Uses control trees" OFF)
option (ENABLE_PORTABLE_TIMER "Uses portable timer" ON)
option (ENABLE_INTERNAL_ERROR_CHECKING "Enable internal checking levels" ON)
set (ERROR_CHECKING_LEVEL 2) #1-Minimal error cheking 2-full error cheking 3-no error cheking
option (ENABLE_MEMLK_CNTR "Enable memory leak counter" OFF)
option (ENABLE_MULTITHREADING "Enables multithreading" ON)
#option (ENABLE_OPENMP "Enables Threading model as Openmp" OFF)
option (ENABLE_SUPER_MATRIX "Enable super matrix" OFF)
option (ENABLE_CBLAS_INTERFACES "Enable cblas interfaces" OFF)
option (ENABLE_GOTO_INTERFACES "Enable goto interfaces" OFF)
option (ENABLE_DEFAULT_BLKSZ "Enable default blocksize" OFF)
if (ENABLE_DEFAULT_BLKSZ)
    set (FLA_DEFAULT_M_BLKSZ 128)
    set (FLA_DEFAULT_N_BLKSZ 128)
    set (FLA_DEFAULT_K_BLKSZ 128)
endif ()
option (ENABLE_MEMORY_ALIGNMENT "Align Dynamic memory boundary" OFF)
option (ENABLE_LIDM_ALIGNMENT "Enable leading dimension alignment" OFF)
option (ENABLE_AUTODETECT_F77_UNDERSCORING "Enable autodetect F77 Underscore" OFF)
option (ENABLE_UPPERCASE_BLAS "Enable invocation of Blas with Flame" OFF)
option (ENABLE_UPPERCASE_LAPACK "Enable invocation of Lapack with Flame" OFF)
option (ENABLE_GPU "Enable gpu advantage" OFF)
option (ENABLE_XBLAS "Enable xblas API's" OFF)
option (ENABLE_F2C_DOTC "Enable invoking Complex datatype BLAS functions DOTC and DOTU using Fortran To C calling convention with void return instead of complex return type" OFF)
option (ENABLE_VOID_RETURN_COMPLEX_FUNCTION "Enable 'void' return instead of complex return type for libFLAME functions such as cladiv/zladiv that by default have complex return type" OFF)
option (ENABLE_SET_LIB_VERSION "Set library version" OFF)
option (ENABLE_BLAS_EXT_GEMMT "Enable BLAS GEMMT Usage" ON)

set (ENABLE_WINDOWS_BUILD "Enables windows build" ON)
SET(ENABLE_SIMD_FLAGS "AVX2" CACHE STRING "Set compiler SIMD flags")
SET_PROPERTY(CACHE ENABLE_SIMD_FLAGS PROPERTY STRINGS none SSE2 AVX AVX2)

if(${ENABLE_SIMD_FLAGS} MATCHES "AVX2")
   add_definitions(/arch:AVX2)
elseif(${ENABLE_SIMD_FLAGS} MATCHES "AVX")
   add_definitions(/arch:AVX)
elseif(${ENABLE_SIMD_FLAGS} MATCHES "SSE2")
   add_definitions(/arch:SSE2)
endif()


if(ENABLE_AMD_FLAGS)
    if(ENABLE_BLAS_EXT_GEMMT EQUAL ON)
    	set (ENABLE_BLAS_EXT_GEMMT TRUE)
    endif ()
    set (ENABLE_AMD_OPT TRUE)
    set (ENABLE_BUILTIN_LAPACK2FLAME TRUE)
    set (ENABLE_EXT_LAPACK_INTERFACE  TRUE)
    set (ENABLE_F2C_DOTC TRUE)
    set (ENABLE_VOID_RETURN_COMPLEX_FUNCTION TRUE)
    if(ENABLE_MULTITHREADING EQUAL ON)
    	set (ENABLE_MULTITHREADING TRUE)
    	set (MULTITHREADING_MODEL 1)
    	set (FLA_MULTITHREADING_MODEL TRUE)
    endif ()
    message(STATUS "Enabled ENABLE_AMD_FLAGS : ${ENABLE_AMD_FLAGS} ")
endif()

if(ENABLE_MULTITHREADING)
    set (ENABLE_OPENMP TRUE)
    message(STATUS "Enabled ENABLE_MULTITHREADING : ${ENABLE_MULTITHREADING} ")
endif()


if(ENABLE_UPPERCASE)
add_definitions(-DBLIS1_ENABLE_UPPERCASE_F77)
add_definitions(-DUPPER)
elseif(NOT ENABLE_UPPERCASE)
add_definitions(-DNOCHANGE)
endif()

if (ENABLE_WRAPPER)
add_definitions(-DFLA_ENABLE_EXTRA_SYMBOLS)
else()
   if(ENABLE_UPPERCASE)
      add_definitions(-DBLIS1_ENABLE_UPPERCASE_F77)
      add_definitions(-DUPPER)
   elseif(NOT ENABLE_UPPERCASE)
      add_definitions(-DNOCHANGE)
   endif()
endif()

if(ENABLE_XBLAS)
add_definitions(-DFLA_ENABLE_XBLAS)
endif()
if (ENABLE_NON_CRITICAL_CODE)
    set (FLA_ENABLE_NON_CRITICAL_CODE TRUE)
endif ()

if(NOT BUILD_SHARED_LIBS)
add_definitions(-DLIBFLAME_IS_BUILDING_LIBRARY)
endif()

if (ENABLE_BUILTIN_BLAS)
    set (FLA_ENABLE_BUILTIN_BLAS TRUE)
endif ()

if (ENABLE_BLIS1_USE_OF_FLA_MALLOC)
    set (FLA_ENABLE_BLIS1_USE_OF_FLA_MALLOC TRUE)
endif ()
if (ENABLE_BUILTIN_LAPACK2FLAME)
    set (FLA_ENABLE_LAPACK2FLAME TRUE)
    message(STATUS "Enabled FLA_ENABLE_LAPACK2FLAME : ${ENABLE_BUILTIN_LAPACK2FLAME} ")
endif ()

if (ENABLE_EXT_LAPACK_SUBPROBLEMS)
    if (ENABLE_BUILTIN_LAPACK2)
        message (FATAL_ERROR "Configuring libflame to enable \
            external-lapack-for-subproblems with lapack2flame \
            is not allowed")
    elseif (ENABLE_EXT_LAPACK_INTERFACE)
        set (FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES  TRUE)
        set (FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS TRUE)
    else ()
        message (FATAL_ERROR "Configuring libflame to enable \
            external-lapack-for-subproblems without \
            external-lapack-interfaces is not allowed")
    endif ()
endif ()

if (ENABLE_EXT_LAPACK_INTERFACE)
    set (FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES  TRUE)
    message(STATUS "Enabled FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES : ${ENABLE_EXT_LAPACK_INTERFACE} ")
endif ()


if (ENABLE_BLAS3_FRNTEND_CNTL_TREES)
    set (FLA_ENABLE_BLAS3_FRONT_END_CNTL_TREES TRUE)
endif ()

if (ENABLE_WINDOWS_BUILD)
	set(FLA_ENABLE_WINDOWS_BUILD TRUE)
endif ()

if (ENABLE_PORTABLE_TIMER)
    set (FLA_ENABLE_PORTABLE_TIMER TRUE)
endif ()


if (ENABLE_INTERNAL_ERROR_CHECKING)
    set (FLA_ENABLE_INTERNAL_ERROR_CHECKING TRUE)
    set (FLA_INTERNAL_ERROR_CHECKING_LEVEL TRUE)
endif ()

if (ENABLE_MEMLK_CNTR)
    set (FLA_ENABLE_MEMORY_LEAK_COUNTER TRUE)
endif ()

if (ENABLE_MULTITHREADING)
        if (ENABLE_AMD_FLAGS)
		set (FLA_OPENMP_MULTITHREADING TRUE)
		add_definitions(-DFLA_OPENMP_MULTITHREADING)
		message(STATUS "Enabled FLA_OPENMP_MULTITHREADING : ${FLA_OPENMP_MULTITHREADING} ")
	else ()
    		set (FLA_ENABLE_MULTITHREADING TRUE)
		add_definitions(-DFLA_ENABLE_MULTITHREADING)
		message(STATUS "Enabled FLA_ENABLE_MULTITHREADING : ${ENABLE_MULTITHREADING} ")
	endif ()

        if (ENABLE_OPENMP)
            find_package(OpenMP)
            if (OPENMP_FOUND)
                set (MULTITHREADING_MODEL 1)
                set (FLA_MULTITHREADING_MODEL TRUE)
		add_definitions(-DFLA_MULTITHREADING_MODEL)
				set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
				set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
				set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
				message(STATUS "Enabled ENABLE_OPENMP : ${ENABLE_OPENMP} ")
            else ()
                message (FATAL_ERROR "Openmp Not Found")
            endif ()
        endif ()
endif ()

if (ENABLE_SUPER_MATRIX)
    set (FLA_ENABLE_SUPERMATRIX TRUE)
endif ()

if (ENABLE_GOTO_INTERFACES)
    set (FLA_ENABLE_GOTO_INTERFACES TRUE)
endif ()

if (ENABLE_CBLAS_INTERFACES)
    set(FLA_ENABLE_CBLAS_INTERFACES TRUE)
endif ()

if (ENABLE_DEFAULT_BLKSZ)
    if (ENABLE_GOTO_INTERFACES)
        message (FATAL_ERROR "Configuring libflame to enable \
            defaultblocksize with goto interfaces is not allowed")
    endif ()
    set (FLA_DEFAULT_M_BLOCKSIZE TRUE)
    set (FLA_DEFAULT_N_BLOCKSIZE TRUE)
    set (FLA_DEFAULT_K_BLOCKSIZE TRUE)
endif ()

if (ENABLE_MEMORY_ALIGNMENT)
    set (FLA_ENABLE_MEMORY_ALIGNMENT TRUE)
    set (FLA_MEMORY_ALIGNMENT_BOUNDARY TRUE)
endif ()

if (ENABLE_LIDM_ALIGNMENT)
    set (FLA_ENABLE_LDIM_ALIGNMENT TRUE)
endif ()

if (ENABLE_AUTODETECT_F77_UNDERSCORING)
    set (FLA_ENABLE_AUTODETECT_F77_UNDERSCORING TRUE)
    set (FLA_F77_UNDERSCORE TRUE)
    set (FLA_F77_EXTRA_UNDERSCORE TRUE)
endif ()

if (ENABLE_UPPERCASE_BLAS)
    set (FLA_ENABLE_UPPERCASE_BLAS TRUE)
endif ()

if (ENABLE_UPPERCASE_LAPACK)
    set (FLA_ENABLE_UPPERCASE_LAPACK TRUE)
endif ()

if (ENABLE_GPU)
        set (FLA_ENABLE_GPU  TRUE)
endif ()

if(BUILD_SHARED_LIBS AND ENABLE_WINDOWS_BUILD)
  set (WINDOWS_FLA_SHARED_BUILD TRUE)
  add_compile_definitions(WINDOWS_FLA_SHARED_BUILD=1)
endif()

if(BUILD_TEST)
  add_compile_definitions(WINDOWS_FLA_TEST=1)
endif()

add_definitions(-DLAPACK_COMPLEX_STRUCTURE)
add_definitions(-DBLIS1_ENABLE_WINDOWS_BUILD)
add_definitions(-DHAVE_LAPACK_CONFIG_H)
add_definitions(-DFLA_ENABLE_LAPACK2FLAME)

if(ENABLE_BLAS_EXT_GEMMT)
set(FLA_ENABLE_BLAS_EXT_GEMMT TRUE)
add_definitions(-DFLA_ENABLE_BLAS_EXT_GEMMT)
message(STATUS "Enabled FLA_ENABLE_BLAS_EXT_GEMMT : ${ENABLE_BLAS_EXT_GEMMT} ")
endif()

if(ENABLE_F2C_DOTC)
set(FLA_ENABLE_F2C_DOTC TRUE)
add_definitions(-DFLA_ENABLE_F2C_DOTC)
message(STATUS "Enabled FLA_ENABLE_F2C_DOTC : ${ENABLE_F2C_DOTC} ")
endif()

if(ENABLE_VOID_RETURN_COMPLEX_FUNCTION)
set(FLA_ENABLE_VOID_RETURN_COMPLEX_FUNCTION TRUE)
add_definitions(-DFLA_ENABLE_VOID_RETURN_COMPLEX_FUNCTION)
message(STATUS "Enabled FLA_ENABLE_VOID_RETURN_COMPLEX_FUNCTION : ${ENABLE_VOID_RETURN_COMPLEX_FUNCTION} ")
endif()

if(ENABLE_ILP64)
set(FLA_ENABLE_ILP64 TRUE)
add_definitions(-DFLA_ENABLE_ILP64)
message(STATUS "Enabled FLA_ENABLE_ILP64 : ${ENABLE_ILP64} ")
endif()

if(ENABLE_AMD_OPT)
set(FLA_ENABLE_AMD_OPT TRUE)
add_definitions(-DFLA_ENABLE_AMD_OPT)
message(STATUS "Enabled FLA_ENABLE_AMD_OPT : ${ENABLE_AMD_OPT} ")
endif()

if(ENABLE_SET_LIB_VERSION)
	add_definitions(-DFLA_LIBFLAME_VERSION=${ENABLE_SET_LIB_VERSION})
else ()
	string(TIMESTAMP TODAY "%Y%m%d")
	file (STRINGS "so_version" VERSION_TEXT)
	list(GET VERSION_TEXT 0 MAJOR_VERSION)
	list(GET VERSION_TEXT 1 MINOR_VERSION)
	set( LIBRARY_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}\tBuild\t${TODAY}")
	add_definitions(-DFLA_LIBFLAME_VERSION=${LIBRARY_VERSION})
endif()

#generates the FLA_config.h header file
configure_file (cmake.config.h.in ${CMAKE_SOURCE_DIR}/windows/build/FLA_config.h @ONLY)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0 ")

#genreating the monolithic headers
find_package(PythonLibs 3 REQUIRED)

set (MK_HEADER_DIR_PATHS "${CMAKE_SOURCE_DIR}/windows/build")
string (APPEND MK_HEADER_DIR_PATHS
	" ${CMAKE_SOURCE_DIR}/src/aocl_dtl/"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/include/"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/blis/include/"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/control/"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/control/init"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/blis/include/"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/control/"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/control/init/"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/hierarchy/include/"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/include/"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/supermatrix/gpu/include/"
    " ${CMAKE_SOURCE_DIR}/src/base/flamec/wrappers/base/wrapper/"
	" ${CMAKE_SOURCE_DIR}/src/base/flamec/supermatrix/include/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/axpy/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/axpy/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/axpyt/c/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/axpyt/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/axpyt/h/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/axpyt/n/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/axpyt/t/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/copy/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/copy/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/copyr/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/copyr/l/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/copyr/u/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/copyt/c/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/copyt/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/copyt/h/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/copyt/n/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/copyt/t/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/scal/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/scal/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/scalr/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/scalr/l/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/1/scalr/u/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/gemv/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/gemv/h/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/gemv/n/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/gemv/t/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/trsv/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/trsv/lc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/trsv/ln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/trsv/lt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/trsv/uc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/trsv/un/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/2/trsv/ut/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/cc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/ch/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/cn/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/ct/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/hc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/hh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/hn/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/ht/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/nc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/nh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/nn/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/nt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/tc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/th/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/tn/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/gemm/tt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/hemm/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/hemm/ll/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/hemm/lu/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/hemm/rl/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/hemm/ru/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/her2k/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/her2k/lh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/her2k/ln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/her2k/uh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/her2k/un/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/herk/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/herk/lh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/herk/ln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/herk/uh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/herk/un/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/symm/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/symm/ll/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/symm/lu/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/symm/rl/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/symm/ru/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/syr2k/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/syr2k/ln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/syr2k/lt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/syr2k/un/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/syr2k/ut/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/syrk/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/syrk/ln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/syrk/lt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/syrk/un/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/syrk/ut/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/llc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/llh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/lln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/llt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/luc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/luh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/lun/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/lut/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/rlc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/rlh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/rln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/rlt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/ruc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/ruh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/run/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trmm/rut/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/llc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/llh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/lln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/llt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/luc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/luh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/lun/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/lut/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/rlc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/rlh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/rln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/rlt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/ruc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/ruh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/run/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/blas/3/trsm/rut/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/bsvd/ext/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/bsvd/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/bsvd/n/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/bsvd/v/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/chol/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/chol/l/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/chol/u/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/hevd/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/hevd/ln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/hevd/lv/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/lu/incpiv/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/lu/incpiv/zaux/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/lu/nopiv/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/lu/nopiv/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/lu/piv/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/lu/piv/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/caqr2ut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/caqr2ut/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/caqrutinc/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/lqut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/lqut/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/qr2ut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/qr2ut/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/qrut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/qrut/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/qrutinc/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/qrutpiv/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/q/qrutpiv/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/svd/ext/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/svd/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/svd/uv/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/tevd/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/tevd/n/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/dec/tevd/v/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/inv/spd/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/inv/tri/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/inv/tri/ln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/inv/tri/lu/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/inv/tri/un/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/inv/tri/uu/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/misc/ttmm/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/misc/ttmm/l/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/misc/ttmm/u/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/misc/uddateut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/misc/uddateut/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/misc/uddateutinc/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/bidiagut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/bidiagut/u/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/eig/gest/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/eig/gest/il/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/eig/gest/iu/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/eig/gest/nl/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/eig/gest/nu/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/hessut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/hessut/vars/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/tridiagut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/red/tridiagut/l/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/sol/lyap/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/sol/lyap/h/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/sol/lyap/n/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/sol/sylv/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/sol/sylv/hh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/sol/sylv/hn/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/sol/sylv/nh/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/sol/sylv/nn/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/accum/tut/fc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/accum/tut/fr/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/accum/tut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/caq2ut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/caq2ut/lhfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/caqutinc/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/caqutinc/lhfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/givens/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/givens/lb/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/givens/lf/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/givens/rb/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/givens/rf/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/h2ut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/h2ut/l/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/h2ut/r/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/hudut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/hudut/l/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/piv/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/piv/ln/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/piv/lt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/piv/rn/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/piv/rt/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/q2ut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/q2ut/lhfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/q2ut/lnfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qudut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qudut/lhfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qudutinc/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qudutinc/lhfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/lhbc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/lhbr/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/lhfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/lhfr/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/lnbc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/lnbr/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/lnfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/lnfr/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/rhbc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/rhbr/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/rhfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/rhfr/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/rnbc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/rnbr/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/rnfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qut/rnfr/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qutinc/front/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qutinc/lhfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/util/app/qutinc/lnfc/flamec/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/x86/"
	" ${CMAKE_SOURCE_DIR}/src/lapacke/LAPACKE/example/"
	" ${CMAKE_SOURCE_DIR}/src/lapacke/LAPACKE/include/"
	" ${CMAKE_SOURCE_DIR}/src/lapacke/LAPACKE/src/"
	" ${CMAKE_SOURCE_DIR}/src/lapacke/LAPACKE/utils/"
	" ${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/"
    " ${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/check"
	" ${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c"
	" ${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/c"
	" ${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/flamec"
	" ${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/install/util/io/"
    " ${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/install/static/"
	)
set (C_COMMENT "-c")
set (VERBOSE "-v1")
set (BASE_INC_PATH "include")
file(MAKE_DIRECTORY ${BASE_INC_PATH})

# Run python script to generate monolithic header at configuration time
# this will generate  blis1.h file in windows/include
set (BLIS1_TARGET_PATH "${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/blis1.h")
set (BLIS1__H_SRC_PATH "${CMAKE_SOURCE_DIR}/src/base/flamec/blis/include/blis1.h" )
set (BLIS_FLATTEN_HEADER "${CMAKE_SOURCE_DIR}/build/flatten-headers.py")
set (PYTHON_EXE "python")

execute_process(
    COMMAND ${PYTHON_EXE} ${BLIS_FLATTEN_HEADER} "${C_COMMENT}" "${VERBOSE}"
	"${BLIS1__H_SRC_PATH}" "${BLIS1_TARGET_PATH}" "${BASE_INC_PATH}" "${MK_HEADER_DIR_PATHS}"
	RESULT_VARIABLE CMD_RESULT
    OUTPUT_VARIABLE CMD_OUTPUT
	)

# Run python script to generate monolithic header at configuration time
# this will generate FLA_f2c.h in windows/include
set (FLAME_TARGET_PATH "${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/FLA_f2c.h")
set (FLAME_H_SRC_PATH  "${CMAKE_SOURCE_DIR}/src/base/flamec/include/FLA_f2c.h" )
execute_process(
    COMMAND ${PYTHON_EXE} ${BLIS_FLATTEN_HEADER} "${C_COMMENT}" "${VERBOSE}"
	"${FLAME_H_SRC_PATH}" "${FLAME_TARGET_PATH}" "${BASE_INC_PATH}" "${MK_HEADER_DIR_PATHS}"
	RESULT_VARIABLE CMD_RESULT
    OUTPUT_VARIABLE CMD_OUTPUT
	)

# Run python script to generate monolithic header at configuration time
# this will generate FLAME.h in windows/include
set (FLAME_TARGET_PATH "${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/FLAME.h")
set (FLAME_H_SRC_PATH  "${CMAKE_SOURCE_DIR}/src/base/flamec/include/FLAME.h" )
execute_process(
    COMMAND ${PYTHON_EXE} ${BLIS_FLATTEN_HEADER} "${C_COMMENT}" "${VERBOSE}"
	"${FLAME_H_SRC_PATH}" "${FLAME_TARGET_PATH}" "${BASE_INC_PATH}" "${MK_HEADER_DIR_PATHS}"
	RESULT_VARIABLE CMD_RESULT
    OUTPUT_VARIABLE CMD_OUTPUT
	)

include_directories(${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/.)
include_directories(${CMAKE_SOURCE_DIR}/src/aocl_dtl)
include_directories(${CMAKE_SOURCE_DIR}/src/base/flamec/include)
include_directories(${CMAKE_SOURCE_DIR}/src/base/flamec/supermatrix/include)
include_directories(${CMAKE_SOURCE_DIR}/src/base/flamec/supermatrix/gpu/include)
include_directories(${CMAKE_SOURCE_DIR}/src/base/flamec/blis/include)
include_directories(${CMAKE_SOURCE_DIR}/src/base/flamec/hierarchy/include)
include_directories(${CMAKE_SOURCE_DIR}/src/base/flamec/control)
include_directories(${CMAKE_SOURCE_DIR}/src/base/flamec/wrappers/base/wrapper)
include_directories(${CMAKE_SOURCE_DIR}/src/map)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/check)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/c)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/flamec)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/flamec/front)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/flamec/gelq)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/flamec/geqr)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/flamec/hetd)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/flamec/spffrt)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/install)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/install/static)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec/f2c/install/util)
include_directories(${CMAKE_SOURCE_DIR}/src/lapacke/LAPACKE/include)
include_directories(${CMAKE_SOURCE_DIR}/src/lapacke/LAPACKE/src)
include_directories(${CMAKE_SOURCE_DIR}/src/lapacke/LAPACKE/utils)
include_directories(${CMAKE_SOURCE_DIR}/src/lapack/x86/front)
include_directories(${CMAKE_SOURCE_DIR}/src/lapack/x86/avx2)
file(GLOB headers ${CMAKE_SOURCE_DIR}/*.h)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MT ")

if (BUILD_SHARED_LIBS)
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
	add_library("${PROJECT_NAME}" SHARED ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/blis1.h
                 ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/FLA_f2c.h
                 ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/FLAME.h
                 )
	if(ENABLE_OPENMP)
		target_link_libraries(${PROJECT_NAME}  OpenMP::OpenMP_CXX)
	endif ()
	target_compile_definitions("${PROJECT_NAME}" PUBLIC -DLIBFLAME_IS_BUILDING_LIBRARY)
	target_compile_definitions("${PROJECT_NAME}" PUBLIC -DLIBFLAME_ENABLE_SHARED)
else ()
	add_library("${PROJECT_NAME}" STATIC ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/blis1.h
                 ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/FLA_f2c.h
                 ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/FLAME.h
                 )
	if(ENABLE_OPENMP)
        	set_target_properties("${PROJECT_NAME}" PROPERTIES STATIC_LIBRARY_OPTIONS "${OpenMP_libomp_LIBRARY}")
	endif ()
endif ()
# TODO : enable this flag to build source supplied blas library
# link externally build blis library
if (NOT ENABLE_BUILTIN_BLAS)
target_link_libraries(${PROJECT_NAME} ${CMAKE_EXT_BLIS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLIS_LIBNAME})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_EXT_BLIS_LIBRARY_DEPENDENCY_PATH})
endif ()

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
add_definitions(-DEXPMODULE)
add_subdirectory(src)

if (BUILD_TEST)
add_subdirectory(test/legacyflame)
add_subdirectory(test/aoclflaprogress)
#add_subdirectory(testcpp)
endif()

if (BUILD_NETLIB_TEST)
add_subdirectory(netlib-test)
elseif(NOT BUILD_NETLIB_TEST)
add_subdirectory(test/main)
endif()

