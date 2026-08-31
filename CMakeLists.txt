# Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.

cmake_minimum_required(VERSION 3.26.0)

# Reading program version from file
file(STRINGS "so_version" VERSION_DATA)
list(JOIN VERSION_DATA "." PROJ_VERSION)

string(APPEND LIBS_LINK_FLAGS "-lpthread -lm ")

if(WIN32)
	project(AOCL-LibFLAME-Win VERSION ${PROJ_VERSION})
elseif(UNIX)
	project(flame VERSION ${PROJ_VERSION})
else()
	message( FATAL_ERROR "OS Undetected")
endif()

# --- Target-based build support -------------------------------------------
# libflame uses the top-level source-dir builtin throughout to refer to its OWN
# source tree, which only holds when libflame is the top-level project. Under
# add_subdirectory/FetchContent that builtin points at the parent (aocl) project
# and every such path breaks. Capture libflame's root ONCE into FLAME_ROOT and
# use ${FLAME_ROOT} for all path references -- uniform with the other unified
# AOCL components' <PKG>_ROOT pattern (a cache variable is visible in all scopes
# and unaffected by nesting).
set(FLAME_ROOT "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "libflame source root")
# --------------------------------------------------------------------------

if(WIN32)
	SET(OpenMP_libomp_LIB "C:/Program Files/LLVM/lib/libomp.lib" CACHE STRING "openmp library path")
	SET(CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH "${FLAME_ROOT}/bin/Release" CACHE STRING "BLAS library path")
	SET(EXT_BLAS_LIBNAME "" CACHE STRING "BLAS library name")
	SET(LIBAOCLUTILS_LIBRARY_PATH "" CACHE STRING "AOCLUtils library name")
	SET(LIBAOCLUTILS_INCLUDE_PATH "" CACHE STRING "AOCLUtils library header files path")
endif()

message(STATUS "External blas library Path : ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}")
message(STATUS "External blas library name : ${EXT_BLAS_LIBNAME}")
SET(LF_ISA_CONFIG "avx2" CACHE STRING "LibFlame ISA configuration: none, avx, avx2, avx2-strict, avx512, avx512-strict, auto")

if(WIN32)
	SET(OpenMP_libomp_LIBRARY "${OpenMP_libomp_LIB}" CACHE STRING "openmp library path")
endif()

message(STATUS "CMAKE_C_COMPILER_ID: ${CMAKE_C_COMPILER_ID}")
message(STATUS "C compiler version: ${CMAKE_C_COMPILER_VERSION}")


#set the output directory path for windows and Linux
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${FLAME_ROOT}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${FLAME_ROOT}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${FLAME_ROOT}/bin")
# Added for debugging purpose
set(CMAKE_VERBOSE_MAKEFILE OFF)

# getting python3 interpreter
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Note: ISA validation and configuration is now handled by cmake/CompilerFlags.cmake
# via the configure_compiler_flags() macro called below

#set AOCL Utilities library path values
set(LIBAOCLUTILS_PROJECT "libaoclutils")
set(LIBAOCLUTILS_INSTALL_DIR_NAME "install-aoclutils")
set(LIBAOCLUTILS_GIT_URL "https://github.com/amd/aocl-utils.git" CACHE STRING "libaoclutils repository URL to be used")
set(LIBAOCLUTILS_GIT_TAG "5.0" CACHE STRING "libaoclutils repository tag to be used")

if(UNIX)
	set(LIBAOCLUTILS_SHAREDLIB "libaoclutils.so")
	set(LIBAOCLUTILS_STATICLIB "libaoclutils.a")
else()
	set(LIBAOCLUTILS_SHAREDLIB "libaoclutils.dll")
	set(LIBAOCLUTILS_STATICLIB "libaoclutils_static.lib")
endif ()

#enabling the configuration options
option (ENABLE_AMD_OPT "Enable AMD optimized path" OFF)
option (ENABLE_AMD_COMMON_FLAGS "Enable AMD FLAGS internally enables ENABLE_BLAS_EXT_GEMMT, ENABLE_AMD_OPT, ENABLE_BUILTIN_LAPACK2FLAME, ENABLE_EXT_LAPACK_INTERFACE, ENABLE_MULTITHREADING" OFF)
# hiding this option from cmake gui
mark_as_advanced(ENABLE_AMD_COMMON_FLAGS)
option (ENABLE_AMD_FLAGS "Enable AMD FLAGS internally enables ENABLE_AMD_COMMON_FLAGS" OFF)
option (ENABLE_AMD_AOCC_FLAGS "Enable AMD FLAGS internally enables Enable AMD FLAGS internally enables ENABLE_AMD_COMMON_FLAGS, ENABLE_F2C_DOTC, ENABLE_VOID_RETURN_COMPLEX_FUNCTION" OFF)
# Since AOCC isn't used in windows hiding this option from cmake gui
if(WIN32)
	mark_as_advanced(ENABLE_AMD_AOCC_FLAGS)
endif()
option (ENABLE_AOCL_BLAS "Enables tight coupling with AOCL-BLAS library in order to use AOCL-BLAS internal routines" OFF)
option (ENABLE_ILP64 "Enable ILP64" OFF)
option (ENABLE_UPPERCASE "Enable uppercase API's" OFF)
option (ENABLE_WRAPPER "Enable wrapper code" ON)
option (INCLUDE_LAPACKE "INCLUDE LAPACKE Code" ON)
option (BUILD_SHARED_LIBS "Building shared Library" ON)
option (BUILD_TEST "Building main test directory" OFF)
option (BUILD_LEGACY_TEST "Building legacy test directory" OFF)
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
option (BUILD_DOC "Build documentation for AOCL-LAPACK (Disabled by default)" OFF)
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
set(ENABLE_AOCL_DTL OFF CACHE STRING "Enable DTL Usage. USE with Option ALL to Enable both Trace and Log, TRACE to enable Trace, and LOG to ENABLE only Input Logging ")
set_property(CACHE ENABLE_AOCL_DTL PROPERTY STRINGS ALL TRACE_AND_LOG TRACE LOG OFF )
option (ENABLE_EMBED_AOCLUTILS "Enables embedding aocl-utils library in AOCL-LAPACK" OFF)
option (ENABLE_LAPACKE_LSAME_FORTRAN_STRLEN " For LAPACKE interface enables lsame arguments length based on Fortran interface " OFF)
option (ENABLE_ASAN "Enable Address Sanitizer (Debug build)" OFF)
# By default build LAPACKE Index-32 API library
option(BUILD_INDEX64 "Build LAPACKE Index-64 API libraries" OFF)
if(BUILD_INDEX64)
  set(LAPACKELIB "lapacke64")
else()
  set(LAPACKELIB "lapacke")
endif()

# By default build LAPACKE Index-64 API as extended API for supported compilers only.
#TODO : Require lapack_64.h to build LAPACKE Index-64 API.
option(BUILD_INDEX64_EXT_API
  "Build LAPACKE Index-64 API as extended API with _64 suffix" OFF)
message(STATUS "Build LAPACKE Index-64 API as extended API with _64 suffix: ${BUILD_INDEX64_EXT_API}")

message(STATUS "CMAKE_BUILD_TYPE : ${CMAKE_BUILD_TYPE}")
if(WIN32)
    set (ENABLE_WINDOWS_BUILD "Enables windows build " ON)
    message(STATUS "ENABLE_WINDOWS_BUILD : ${ENABLE_WINDOWS_BUILD}")
    cmake_policy(SET CMP0091 NEW)
    if(BUILD_SHARED_LIBS)
       set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
    else()
       set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
endif()
set(AOCL_ROOT "" CACHE STRING "AOCL-BLAS installation path")

# ============================================================================
# Include CMake modules for better organization
# ============================================================================
include(cmake/CompilerFlags.cmake)

# Configure all compiler flags using the modular approach
configure_compiler_flags()

message(STATUS "Compiler flags configured")

# Adding missing Macros
if(UNIX)
	set(FLA_PORTABLE_TIMER_IS_CLOCK_GETTIME TRUE)
	set(FLA_VECTOR_INTRINSIC_TYPE TRUE)
	set(HAVE_ASSERT_H TRUE)
	set(HAVE_FCNTL_H TRUE)
	set(HAVE_INTTYPES_H TRUE)
	set(HAVE_LIBM TRUE)
	set(HAVE_MATH_H TRUE)
	set(HAVE_MEMORY_H TRUE)
	set(HAVE_SIGNAL_H TRUE)
	set(HAVE_STDINT_H TRUE)
	set(HAVE_STDLIB_H TRUE)
	set(HAVE_STRINGS_H TRUE)
	set(HAVE_STRING_H TRUE)
	set(HAVE_SYS_STAT_H TRUE)
	set(HAVE_SYS_TIME_H TRUE)
	set(HAVE_SYS_TYPES_H TRUE)
	set(HAVE_UNISTD_H TRUE)
	set(PROTOTYPES TRUE)
	set(STDC_HEADERS TRUE)
	set(TIME_WITH_SYS_TIME TRUE)
	set(_GNU_SOURCE TRUE)
	set(__PROTOTYPES TRUE)
endif()

# ENABLE AOCL Trace , Log and Dump form CMAKE
if(ENABLE_AOCL_DTL)
	message(STATUS "ENABLE_AOCL_DTL : ${ENABLE_AOCL_DTL}")
	if(ENABLE_AOCL_DTL STREQUAL "ALL" OR ENABLE_AOCL_DTL STREQUAL "TRACE_AND_LOG" )
		add_definitions(-DLF_AOCL_DTL_TRACE_ENABLE=1 -DLF_AOCL_DTL_LOG_ENABLE=1)
	elseif(ENABLE_AOCL_DTL STREQUAL "TRACE")
		add_definitions(-DLF_AOCL_DTL_TRACE_ENABLE=1)
	elseif(ENABLE_AOCL_DTL STREQUAL "LOG")
		add_definitions(-DLF_AOCL_DTL_LOG_ENABLE=1)
	else()
	message(FATAL_ERROR "INVALID PARAM for -DENABLE_AOCL_DTL=${ENABLE_AOCL_DTL} \n PLEASE PASS \n \
	 1. \"ALL\" to ENABLE TRACE and LOG \n \
	 2. \"TRACE\" to ENABLE TRACE and \n \
	 3. \"LOG\" to ENABLE LOG \n \
	 4. \"OFF\" to Disable trace and log ")
	endif()
endif()

if(ENABLE_AMD_FLAGS)
	set(ENABLE_AMD_COMMON_FLAGS TRUE)
	if(WIN32)
	# Since windows uses only enable AMD flags
		set (ENABLE_F2C_DOTC TRUE)
		set (ENABLE_VOID_RETURN_COMPLEX_FUNCTION TRUE)
	endif()
	message(STATUS "Enabled ENABLE_AMD_FLAGS : ${ENABLE_AMD_FLAGS}")
endif()

if(ENABLE_AMD_AOCC_FLAGS)
	# time.h requires std=c11 when compiling with AOCC otherwise we get compilation issues
	if(ENABLE_AOCL_DTL)
		set	(CMAKE_C_STANDARD 11)
		set	(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11")
	endif()
    set (ENABLE_AMD_COMMON_FLAGS TRUE)
    set (ENABLE_F2C_DOTC TRUE)
    set (ENABLE_VOID_RETURN_COMPLEX_FUNCTION TRUE)
	message(STATUS "Enabled ENABLE_AMD_AOCC_FLAGS : ${ENABLE_AMD_AOCC_FLAGS}")
endif()

if(ENABLE_AMD_COMMON_FLAGS)
    if(ENABLE_BLAS_EXT_GEMMT)
    	set (ENABLE_BLAS_EXT_GEMMT TRUE)
    endif ()
    set (ENABLE_AMD_OPT TRUE)
    set (ENABLE_BUILTIN_LAPACK2FLAME TRUE)
    set (ENABLE_EXT_LAPACK_INTERFACE TRUE)
    if(ENABLE_MULTITHREADING)
    	set (ENABLE_MULTITHREADING TRUE)
    	set (MULTITHREADING_MODEL 1)
		set(FLA_MULTITHREADING_MODEL 1)
    endif ()
	message(STATUS "Enabled ENABLE_AMD_COMMON_FLAGS : ${ENABLE_AMD_COMMON_FLAGS}")
endif()

# STATUS of AMD/AOCC FLAGS
message(STATUS "Enabled ENABLE_BLAS_EXT_GEMMT : ${ENABLE_BLAS_EXT_GEMMT} ")
message(STATUS "Enabled ENABLE_AMD_OPT : ${ENABLE_AMD_OPT} ")
message(STATUS "Enabled ENABLE_BUILTIN_LAPACK2FLAME : ${ENABLE_BUILTIN_LAPACK2FLAME} ")
message(STATUS "Enabled ENABLE_EXT_LAPACK_INTERFACE : ${ENABLE_EXT_LAPACK_INTERFACE} ")
message(STATUS "Enabled ENABLE_F2C_DOTC : ${ENABLE_F2C_DOTC} ")
message(STATUS "Enabled ENABLE_VOID_RETURN_COMPLEX_FUNCTION	: ${ENABLE_VOID_RETURN_COMPLEX_FUNCTION} ")
message(STATUS "Enabled ENABLE_MULTITHREADING : ${ENABLE_MULTITHREADING} ")
message(STATUS "Enabled MULTITHREADING_MODEL : ${MULTITHREADING_MODEL} ")

if(ENABLE_MULTITHREADING)
    set (ENABLE_OPENMP TRUE)
endif()
message(STATUS "Enabled ENABLE_OPENMP : ${ENABLE_OPENMP} ")

if(ENABLE_UPPERCASE)
add_definitions(-DBLIS1_ENABLE_UPPERCASE_F77)
add_definitions(-DUPPER)
elseif(NOT ENABLE_UPPERCASE)
add_definitions(-DNOCHANGE)
endif()

if (ENABLE_WRAPPER)
set(FLA_ENABLE_EXTRA_SYMBOLS TRUE)
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
endif ()

if (ENABLE_EXT_LAPACK_SUBPROBLEMS)
    if (ENABLE_BUILTIN_LAPACK2)
        message (FATAL_ERROR "Configuring libflame to enable \
            external-lapack-for-subproblems with lapack2flame \
            is not allowed")
    elseif (ENABLE_EXT_LAPACK_INTERFACE)
	set(FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES TRUE)
	set(FLA_ENABLE_EXTERNAL_LAPACK_FOR_SUBPROBLEMS TRUE)
    else ()
        message (FATAL_ERROR "Configuring libflame to enable \
            external-lapack-for-subproblems without \
            external-lapack-interfaces is not allowed")
    endif ()
endif ()

if (ENABLE_EXT_LAPACK_INTERFACE)
	set(FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES TRUE)
endif ()


if (ENABLE_BLAS3_FRNTEND_CNTL_TREES)
	set(FLA_ENABLE_BLAS3_FRONT_END_CNTL_TREES TRUE)
endif ()

if (ENABLE_WINDOWS_BUILD)
	add_definitions(-DFLA_ENABLE_WINDOWS_BUILD)
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
    if (ENABLE_AMD_FLAGS OR ENABLE_AMD_AOCC_FLAGS)
	set (FLA_OPENMP_MULTITHREADING TRUE)
	add_definitions(-DFLA_OPENMP_MULTITHREADING)
	message(STATUS "Enabled FLA_OPENMP_MULTITHREADING : ${FLA_OPENMP_MULTITHREADING} ")
    else ()
    	set (FLA_ENABLE_MULTITHREADING TRUE)
	add_definitions(-DFLA_ENABLE_MULTITHREADING=)
	message(STATUS "Enabled FLA_ENABLE_MULTITHREADING : ${FLA_ENABLE_MULTITHREADING} ")
    endif ()

    if (ENABLE_OPENMP)
		set (MULTITHREADING_MODEL 1)
		set (FLA_MULTITHREADING_MODEL TRUE)
		add_definitions(-DFLA_MULTITHREADING_MODEL)
		if(NOT EXT_OPENMP_PATH)
			find_package(OpenMP)
			if (OPENMP_FOUND)
					set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
					set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
					set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
					message(STATUS "Enabled ENABLE_OPENMP : ${ENABLE_OPENMP} ")
			else ()
					message (FATAL_ERROR "Openmp Not Found")
			endif ()
		endif()
    endif ()
endif ()

if (ENABLE_SUPER_MATRIX)
	set(FLA_ENABLE_SUPERMATRIX TRUE)
endif ()

if (ENABLE_GOTO_INTERFACES)
	set(FLA_ENABLE_GOTO_INTERFACES TRUE)
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

if(BUILD_LEGACY_TEST)
	if(WIN32)
		add_compile_definitions(WINDOWS_FLA_TEST=1)
	endif()
endif()

if (NOT ENABLE_AOCL_BLAS)
add_definitions(-DBLIS1_FROM_LIBFLAME)
endif()

add_definitions(-DLAPACK_COMPLEX_STRUCTURE)


add_definitions(-DHAVE_LAPACK_CONFIG_H)
add_definitions(-DFLA_ENABLE_LAPACK2FLAME)

if(ENABLE_BLAS_EXT_GEMMT)
	set(FLA_ENABLE_BLAS_EXT_GEMMT TRUE)
endif()

if(ENABLE_F2C_DOTC)
	set(FLA_ENABLE_F2C_DOTC TRUE)
	add_definitions(-DFLA_ENABLE_F2C_DOTC=1)
endif()

if(ENABLE_VOID_RETURN_COMPLEX_FUNCTION)
	set(FLA_ENABLE_VOID_RETURN_COMPLEX_FUNCTION TRUE)
	add_definitions(-DFLA_ENABLE_VOID_RETURN_COMPLEX_FUNCTION)
	message(STATUS "Enabled FLA_ENABLE_VOID_RETURN_COMPLEX_FUNCTION : ${ENABLE_VOID_RETURN_COMPLEX_FUNCTION} ")
endif()

if(ENABLE_ILP64)
	set(FLA_ENABLE_ILP64 TRUE)
	#BUILD_INDEX64_EXT_API is disabled until _64 suffix APIs is supported from dependent libraries.
	set(BUILD_INDEX64_EXT_API OFF)
	add_definitions(-DFLA_ENABLE_ILP64=1)
	message(STATUS "Enabled FLA_ENABLE_ILP64 : ${ENABLE_ILP64} ")
endif()

if(ENABLE_LAPACKE_LSAME_FORTRAN_STRLEN)
    add_definitions(-DLAPACK_LSAME_FORTRAN_STRLEN_END)
endif()

set(REQ_LIBAOCLUTILS_PKGNAME "")
# Configuration for pkg-config
if(NOT ENABLE_EMBED_AOCLUTILS)
	# set the variable if aocl-utils library is not embedded
	set(REQ_LIBAOCLUTILS_PKGNAME "aocl-utils")
endif()
# Get the required AOCL-BLAS library pkg name
if (ENABLE_MULTITHREADING)
    set(REQ_BLAS_PKGNAME "blis-mt")
else()
	set(REQ_BLAS_PKGNAME "blis")
endif()

# Include the AOCL-BLAS path check cmake directory.
# --- Target-based build support -------------------------------------------
# When libflame is pulled into a parent build (FetchContent/add_subdirectory)
# that already provides AOCL-BLAS as an in-tree CMake target, the BLAS library
# is NOT installed yet at configure time, so the install-path discovery in
# find_dep_libs.cmake cannot find it. In that case, consume the BLAS target
# directly: target-linking it propagates its headers to our compilation, and
# the parent supplies the BLAS include dir via AOCL_TB_BLAS_INCLUDE_DIR.
if(TARGET AOCL::BLAS)
	set(BLAS_LIBRARY AOCL::BLAS)
	# Prefer the include dirs propagated by the target; fall back to the
	# parent-provided cache variable (if any) for unified builds.
	get_target_property(_blas_incdirs AOCL::BLAS INTERFACE_INCLUDE_DIRECTORIES)
	if(_blas_incdirs)
		set(BLAS_INCLUDE_DIR "${_blas_incdirs}")
	else()
		set(BLAS_INCLUDE_DIR "${AOCL_TB_BLAS_INCLUDE_DIR}")
	endif()
	message(STATUS "[libflame] target-based build: using in-tree AOCL::BLAS target")
	message(STATUS "[libflame] BLAS include dir: ${BLAS_INCLUDE_DIR}")
else()
	include(cmake/find_dep_libs.cmake)
endif()
# --------------------------------------------------------------------------

# AOCL-BLAS library linking option. If ENABLE_AOCL_BLAS is set, then the AOCL-BLAS library will be linked.
# Note: This block has to come after the FLA_ENABLE_ILP64 block, as that decides the type of AOCL-BLAS library thats linked.
if (ENABLE_AOCL_BLAS)
    message(STATUS "AOCL-BLAS library linking option selected")

	if(NOT BLAS_LIBRARY)
		message (FATAL_ERROR "Error: Either AOCL_ROOT environment variable is not set or could not find a suitable installation of AOCL-BLAS Library in \$AOCL_ROOT=${AOCL_ROOT}")
	endif()
endif()

# If the external BLAS library path is not set, then check if the BLAS_LIBRARY is set from AOCL_ROOT environment. 
# If it is set, then use the BLAS_LIBRARY path as the AOCL-BLAS path.
if ((ENABLE_AOCL_BLAS OR NOT (CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH AND EXT_BLAS_LIBNAME AND BLAS_HEADER_PATH))
	AND (WIN32 OR BUILD_TEST OR BUILD_LEGACY_TEST OR BUILD_NETLIB_TEST) 
	AND BLAS_LIBRARY)

	if (ENABLE_AOCL_BLAS AND (CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH OR EXT_BLAS_LIBNAME OR BLAS_HEADER_PATH))
		message(WARNING "CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH and/or EXT_BLAS_LIBNAME are set. These will be ignored since BLAS library will be chosen from AOCL_ROOT path when ENABLE_AOCL_BLAS is set")
	endif()

	set(CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH ${BLAS_LIB_PATH})
	set(EXT_BLAS_LIBNAME ${BLAS_LIB_NAME})
	set(BLAS_HEADER_PATH ${BLAS_INCLUDE_DIR})
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
	if(UNIX)
	set(FLA_LIBFLAME_VERSION TRUE)
	set( LIBRARY_VERSION "\"${MAJOR_VERSION}.${MINOR_VERSION}\tBuild ${TODAY}\"")
	elseif(WIN32)
	set( LIBRARY_VERSION "\"${MAJOR_VERSION}.${MINOR_VERSION}\tBuild ${TODAY}\"")
	endif()
	add_definitions(-DFLA_LIBFLAME_VERSION=${LIBRARY_VERSION})
endif()
#generates the FLA_config.h header file
if(WIN32)
	configure_file (cmake.config.h.in ${FLAME_ROOT}/windows/build/FLA_config.h @ONLY)
elseif(UNIX)
# As per the current build system
	configure_file(cmake.config.h.in ${FLAME_ROOT}/build/FLA_config.h @ONLY)
else()
	message(FATAL_ERROR " OS not found")
endif()

if(WIN32)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
endif()

#genreating the monolithic headers
if(WIN32)
	set (MK_HEADER_DIR_PATHS "${FLAME_ROOT}/windows/build")
elseif(UNIX)
	set(MK_HEADER_DIR_PATHS "${FLAME_ROOT}/build")
else()
	message(FATAL_ERROR "OS NOT DETECTED")
endif()
string (APPEND MK_HEADER_DIR_PATHS
	" ${FLAME_ROOT}/src/aocl_dtl/"
	" ${FLAME_ROOT}/src/base/flamec/include/"
	" ${FLAME_ROOT}/src/base/flamec/blis/include/"
	" ${FLAME_ROOT}/src/base/flamec/control/"
	" ${FLAME_ROOT}/src/base/flamec/control/init"
	" ${FLAME_ROOT}/src/base/flamec/blis/include/"
	" ${FLAME_ROOT}/src/base/flamec/control/"
	" ${FLAME_ROOT}/src/base/flamec/control/init/"
	" ${FLAME_ROOT}/src/base/flamec/hierarchy/include/"
	" ${FLAME_ROOT}/src/base/flamec/include/"
	" ${FLAME_ROOT}/src/base/flamec/supermatrix/gpu/include/"
    " ${FLAME_ROOT}/src/base/flamec/wrappers/base/wrapper/"
	" ${FLAME_ROOT}/src/base/flamec/supermatrix/include/"
	" ${FLAME_ROOT}/src/blas/1/axpy/front/flamec/"
	" ${FLAME_ROOT}/src/blas/1/axpy/vars/flamec/"
	" ${FLAME_ROOT}/src/blas/1/axpyt/c/flamec/"
	" ${FLAME_ROOT}/src/blas/1/axpyt/front/flamec/"
	" ${FLAME_ROOT}/src/blas/1/axpyt/h/flamec/"
	" ${FLAME_ROOT}/src/blas/1/axpyt/n/flamec/"
	" ${FLAME_ROOT}/src/blas/1/axpyt/t/flamec/"
	" ${FLAME_ROOT}/src/blas/1/copy/front/flamec/"
	" ${FLAME_ROOT}/src/blas/1/copy/vars/flamec/"
	" ${FLAME_ROOT}/src/blas/1/copyr/front/flamec/"
	" ${FLAME_ROOT}/src/blas/1/copyr/l/flamec/"
	" ${FLAME_ROOT}/src/blas/1/copyr/u/flamec/"
	" ${FLAME_ROOT}/src/blas/1/copyt/c/flamec/"
	" ${FLAME_ROOT}/src/blas/1/copyt/front/flamec/"
	" ${FLAME_ROOT}/src/blas/1/copyt/h/flamec/"
	" ${FLAME_ROOT}/src/blas/1/copyt/n/flamec/"
	" ${FLAME_ROOT}/src/blas/1/copyt/t/flamec/"
	" ${FLAME_ROOT}/src/blas/1/scal/front/flamec/"
	" ${FLAME_ROOT}/src/blas/1/scal/vars/flamec/"
	" ${FLAME_ROOT}/src/blas/1/scalr/front/flamec/"
	" ${FLAME_ROOT}/src/blas/1/scalr/l/flamec/"
	" ${FLAME_ROOT}/src/blas/1/scalr/u/flamec/"
	" ${FLAME_ROOT}/src/blas/2/gemv/front/flamec/"
	" ${FLAME_ROOT}/src/blas/2/gemv/h/flamec/"
	" ${FLAME_ROOT}/src/blas/2/gemv/n/flamec/"
	" ${FLAME_ROOT}/src/blas/2/gemv/t/flamec/"
	" ${FLAME_ROOT}/src/blas/2/trsv/front/flamec/"
	" ${FLAME_ROOT}/src/blas/2/trsv/lc/flamec/"
	" ${FLAME_ROOT}/src/blas/2/trsv/ln/flamec/"
	" ${FLAME_ROOT}/src/blas/2/trsv/lt/flamec/"
	" ${FLAME_ROOT}/src/blas/2/trsv/uc/flamec/"
	" ${FLAME_ROOT}/src/blas/2/trsv/un/flamec/"
	" ${FLAME_ROOT}/src/blas/2/trsv/ut/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/cc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/ch/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/cn/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/ct/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/front/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/hc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/hh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/hn/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/ht/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/nc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/nh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/nn/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/nt/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/tc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/th/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/tn/flamec/"
	" ${FLAME_ROOT}/src/blas/3/gemm/tt/flamec/"
	" ${FLAME_ROOT}/src/blas/3/hemm/front/flamec/"
	" ${FLAME_ROOT}/src/blas/3/hemm/ll/flamec/"
	" ${FLAME_ROOT}/src/blas/3/hemm/lu/flamec/"
	" ${FLAME_ROOT}/src/blas/3/hemm/rl/flamec/"
	" ${FLAME_ROOT}/src/blas/3/hemm/ru/flamec/"
	" ${FLAME_ROOT}/src/blas/3/her2k/front/flamec/"
	" ${FLAME_ROOT}/src/blas/3/her2k/lh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/her2k/ln/flamec/"
	" ${FLAME_ROOT}/src/blas/3/her2k/uh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/her2k/un/flamec/"
	" ${FLAME_ROOT}/src/blas/3/herk/front/flamec/"
	" ${FLAME_ROOT}/src/blas/3/herk/lh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/herk/ln/flamec/"
	" ${FLAME_ROOT}/src/blas/3/herk/uh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/herk/un/flamec/"
	" ${FLAME_ROOT}/src/blas/3/symm/front/flamec/"
	" ${FLAME_ROOT}/src/blas/3/symm/ll/flamec/"
	" ${FLAME_ROOT}/src/blas/3/symm/lu/flamec/"
	" ${FLAME_ROOT}/src/blas/3/symm/rl/flamec/"
	" ${FLAME_ROOT}/src/blas/3/symm/ru/flamec/"
	" ${FLAME_ROOT}/src/blas/3/syr2k/front/flamec/"
	" ${FLAME_ROOT}/src/blas/3/syr2k/ln/flamec/"
	" ${FLAME_ROOT}/src/blas/3/syr2k/lt/flamec/"
	" ${FLAME_ROOT}/src/blas/3/syr2k/un/flamec/"
	" ${FLAME_ROOT}/src/blas/3/syr2k/ut/flamec/"
	" ${FLAME_ROOT}/src/blas/3/syrk/front/flamec/"
	" ${FLAME_ROOT}/src/blas/3/syrk/ln/flamec/"
	" ${FLAME_ROOT}/src/blas/3/syrk/lt/flamec/"
	" ${FLAME_ROOT}/src/blas/3/syrk/un/flamec/"
	" ${FLAME_ROOT}/src/blas/3/syrk/ut/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/front/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/llc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/llh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/lln/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/llt/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/luc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/luh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/lun/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/lut/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/rlc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/rlh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/rln/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/rlt/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/ruc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/ruh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/run/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trmm/rut/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/front/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/llc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/llh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/lln/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/llt/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/luc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/luh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/lun/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/lut/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/rlc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/rlh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/rln/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/rlt/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/ruc/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/ruh/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/run/flamec/"
	" ${FLAME_ROOT}/src/blas/3/trsm/rut/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/bsvd/ext/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/bsvd/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/bsvd/n/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/bsvd/v/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/chol/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/chol/l/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/chol/u/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/hevd/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/hevd/ln/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/hevd/lv/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/lu/incpiv/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/lu/incpiv/zaux/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/lu/nopiv/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/lu/nopiv/vars/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/lu/piv/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/lu/piv/vars/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/caqr2ut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/caqr2ut/vars/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/caqrutinc/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/lqut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/lqut/vars/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/qr2ut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/qr2ut/vars/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/qrut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/qrut/vars/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/qrutinc/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/qrutpiv/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/q/qrutpiv/vars/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/svd/ext/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/svd/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/svd/uv/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/tevd/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/tevd/n/flamec/"
	" ${FLAME_ROOT}/src/lapack/dec/tevd/v/flamec/"
	" ${FLAME_ROOT}/src/lapack/fblas/"
	" ${FLAME_ROOT}/src/lapack/fblas/front"
	" ${FLAME_ROOT}/src/lapack/fblas/src"
	" ${FLAME_ROOT}/src/lapack/inv/spd/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/inv/tri/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/inv/tri/ln/flamec/"
	" ${FLAME_ROOT}/src/lapack/inv/tri/lu/flamec/"
	" ${FLAME_ROOT}/src/lapack/inv/tri/un/flamec/"
	" ${FLAME_ROOT}/src/lapack/inv/tri/uu/flamec/"
	" ${FLAME_ROOT}/src/lapack/misc/ttmm/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/misc/ttmm/l/flamec/"
	" ${FLAME_ROOT}/src/lapack/misc/ttmm/u/flamec/"
	" ${FLAME_ROOT}/src/lapack/misc/uddateut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/misc/uddateut/vars/flamec/"
	" ${FLAME_ROOT}/src/lapack/misc/uddateutinc/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/bidiagut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/bidiagut/u/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/eig/gest/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/eig/gest/il/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/eig/gest/iu/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/eig/gest/nl/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/eig/gest/nu/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/hessut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/hessut/vars/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/tridiagut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/red/tridiagut/l/flamec/"
	" ${FLAME_ROOT}/src/lapack/sol/lyap/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/sol/lyap/h/flamec/"
	" ${FLAME_ROOT}/src/lapack/sol/lyap/n/flamec/"
	" ${FLAME_ROOT}/src/lapack/sol/sylv/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/sol/sylv/hh/flamec/"
	" ${FLAME_ROOT}/src/lapack/sol/sylv/hn/flamec/"
	" ${FLAME_ROOT}/src/lapack/sol/sylv/nh/flamec/"
	" ${FLAME_ROOT}/src/lapack/sol/sylv/nn/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/accum/tut/fc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/accum/tut/fr/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/accum/tut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/caq2ut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/caq2ut/lhfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/caqutinc/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/caqutinc/lhfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/givens/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/givens/lb/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/givens/lf/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/givens/rb/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/givens/rf/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/h2ut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/h2ut/l/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/h2ut/r/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/hudut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/hudut/l/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/piv/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/piv/ln/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/piv/lt/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/piv/rn/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/piv/rt/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/q2ut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/q2ut/lhfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/q2ut/lnfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qudut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qudut/lhfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qudutinc/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qudutinc/lhfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/lhbc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/lhbr/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/lhfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/lhfr/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/lnbc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/lnbr/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/lnfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/lnfr/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/rhbc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/rhbr/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/rhfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/rhfr/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/rnbc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/rnbr/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/rnfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qut/rnfr/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qutinc/front/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qutinc/lhfc/flamec/"
	" ${FLAME_ROOT}/src/lapack/util/app/qutinc/lnfc/flamec/"
	" ${FLAME_ROOT}/src/lapacke/LAPACKE/example/"
	" ${FLAME_ROOT}/src/lapacke/LAPACKE/include/"
	" ${FLAME_ROOT}/src/lapacke/LAPACKE/src/"
	" ${FLAME_ROOT}/src/lapacke/LAPACKE/utils/"
	" ${FLAME_ROOT}/src/map/lapack2flamec/"
    " ${FLAME_ROOT}/src/map/lapack2flamec/check"
	" ${FLAME_ROOT}/src/map/lapack2flamec/f2c"
	" ${FLAME_ROOT}/src/map/lapack2flamec/f2c/c"
	" ${FLAME_ROOT}/src/map/lapack2flamec/f2c/flamec"
	" ${FLAME_ROOT}/src/map/lapack2flamec/f2c/install/util/io/"
    " ${FLAME_ROOT}/src/map/lapack2flamec/f2c/install/static/"
	)

set (C_COMMENT "-c")
set (VERBOSE "-v1")
set (BASE_INC_PATH "include")
file(MAKE_DIRECTORY ${BASE_INC_PATH})
# Run python script to generate monolithic header at configuration time
# this will generate  blis1.h file in windows/include
set (BLIS1_TARGET_PATH "${FLAME_ROOT}/${BASE_INC_PATH}/blis1.h")
set (BLIS1__H_SRC_PATH "${FLAME_ROOT}/src/base/flamec/blis/include/blis1.h" )
set (FLA_FLATTEN_HEADER "${FLAME_ROOT}/build/flatten-headers.py")

execute_process(
    COMMAND ${Python3_EXECUTABLE} ${FLA_FLATTEN_HEADER} "${C_COMMENT}" "${VERBOSE}"
	"${BLIS1__H_SRC_PATH}" "${BLIS1_TARGET_PATH}" "${BASE_INC_PATH}" "${MK_HEADER_DIR_PATHS}"
	RESULT_VARIABLE CMD_RESULT
    OUTPUT_VARIABLE CMD_OUTPUT
	)

# Run python script to generate monolithic header at configuration time
# this will generate FLA_f2c.h in windows/include
set (FLAME_TARGET_PATH "${FLAME_ROOT}/${BASE_INC_PATH}/FLA_f2c.h")
set (FLAME_H_SRC_PATH  "${FLAME_ROOT}/src/base/flamec/include/FLA_f2c.h" )

execute_process(
    COMMAND ${Python3_EXECUTABLE} ${FLA_FLATTEN_HEADER} "${C_COMMENT}" "${VERBOSE}"
	"${FLAME_H_SRC_PATH}" "${FLAME_TARGET_PATH}" "${BASE_INC_PATH}" "${MK_HEADER_DIR_PATHS}"
	RESULT_VARIABLE CMD_RESULT
    OUTPUT_VARIABLE CMD_OUTPUT
	)

# Run python script to generate monolithic header at configuration time
# this will generate FLAME.h in windows/include
set (FLAME_TARGET_PATH "${FLAME_ROOT}/${BASE_INC_PATH}/FLAME.h")
set (FLAME_H_SRC_PATH  "${FLAME_ROOT}/src/base/flamec/include/FLAME.h" )
execute_process(
    COMMAND ${Python3_EXECUTABLE} ${FLA_FLATTEN_HEADER} "${C_COMMENT}" "${VERBOSE}"
	"${FLAME_H_SRC_PATH}" "${FLAME_TARGET_PATH}" "${BASE_INC_PATH}" "${MK_HEADER_DIR_PATHS}"
	RESULT_VARIABLE CMD_RESULT
    OUTPUT_VARIABLE CMD_OUTPUT
	)

# Run python script to generate libflame_interface header at configuration time
set (INTERFACE_HEADER_DIR_PATHS "${FLAME_ROOT}/src/src_cpp")
set (INTERFACE_TARGET_PATH "${FLAME_ROOT}/${BASE_INC_PATH}/libflame_interface.hh")
execute_process(
	COMMAND ${Python3_EXECUTABLE} ${FLA_FLATTEN_HEADER}  "${VERBOSE}"
	"${INTERFACE_HEADER_DIR_PATHS}/libflame_interface.hh" "${INTERFACE_TARGET_PATH}"
	"${BASE_INC_PATH}" "${INTERFACE_HEADER_DIR_PATHS}"
	RESULT_VARIABLE CMD_RESULT
    OUTPUT_VARIABLE CMD_OUTPUT
	)

# Run python script to generate libflame_interface header at configuration time
set (LAPACKE_HEADER_SRC_PATH "${FLAME_ROOT}/src/lapacke/LAPACKE/include/")
set (LAPACKE_HEADER_TARGET_PATH "${FLAME_ROOT}/${BASE_INC_PATH}/lapacke.h")
execute_process(
	COMMAND ${Python3_EXECUTABLE} ${FLA_FLATTEN_HEADER}  "${VERBOSE}"
	"${LAPACKE_HEADER_SRC_PATH}/lapacke.h" "${LAPACKE_HEADER_TARGET_PATH}"
	"${BASE_INC_PATH}" "${LAPACKE_HEADER_SRC_PATH}"
	RESULT_VARIABLE CMD_RESULT
    OUTPUT_VARIABLE CMD_OUTPUT
	)

include_directories(${FLAME_ROOT}/${BASE_INC_PATH}/.)
if(ENABLE_AOCL_DTL)
	#only enabling aocl dtl when ENABLE_AOCL_DTL is set to a valid value
	include_directories(${FLAME_ROOT}/src/aocl_dtl)
endif()
include_directories(${FLAME_ROOT}/src/base/flamec/include)
include_directories(${FLAME_ROOT}/src/base/flamec/supermatrix/include)
include_directories(${FLAME_ROOT}/src/base/flamec/supermatrix/gpu/include)
include_directories(${FLAME_ROOT}/src/base/flamec/blis/include)
include_directories(${FLAME_ROOT}/src/base/flamec/hierarchy/include)
include_directories(${FLAME_ROOT}/src/base/flamec/control)
include_directories(${FLAME_ROOT}/src/base/flamec/wrappers/base/wrapper)
include_directories(${FLAME_ROOT}/src/map)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/check)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c/c)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c/flamec)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c/flamec/front)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c/flamec/gelq)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c/flamec/geqr)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c/flamec/hetd)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c/flamec/spffrt)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c/install)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c/install/static)
include_directories(${FLAME_ROOT}/src/map/lapack2flamec/f2c/install/util)
include_directories(${FLAME_ROOT}/src/lapacke/LAPACKE/include)
include_directories(${FLAME_ROOT}/src/lapacke/LAPACKE/src)
include_directories(${FLAME_ROOT}/src/lapacke/LAPACKE/utils)
include_directories(${FLAME_ROOT}/src/lapack/x86/front)
include_directories(${FLAME_ROOT}/src/lapack/x86/avx2)
include_directories(${FLAME_ROOT}/src/lapack/x86/avx512)
include_directories(${FLAME_ROOT}/src/lapack/fblas/front)
include_directories(${FLAME_ROOT}/src/lapack/fblas/src)
if(ENABLE_AOCL_BLAS)
	include_directories(${BLAS_INCLUDE_DIR})
endif()

file(GLOB headers ${FLAME_ROOT}/*.h)

if(WIN32)
	cmake_policy(SET CMP0091 NEW)
	if(BUILD_SHARED_LIBS)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
	else()
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()
endif()

set(OBJECT_LIBRARIES
    $<TARGET_OBJECTS:FLA_LAPACK_AVX512>
    $<TARGET_OBJECTS:FLA_LAPACK_AVX2>
    $<TARGET_OBJECTS:${LAPACKELIB}_obj>
    $<$<BOOL:${BUILD_INDEX64_EXT_API}>: $<TARGET_OBJECTS:${LAPACKELIB}_64_obj>>
)

if (BUILD_SHARED_LIBS)
	if(WIN32)
		set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
	endif()
	add_library("${PROJECT_NAME}" SHARED ${FLAME_ROOT}/${BASE_INC_PATH}/blis1.h
                 ${FLAME_ROOT}/${BASE_INC_PATH}/FLA_f2c.h
                 ${FLAME_ROOT}/${BASE_INC_PATH}/FLAME.h
                 ${OBJECT_LIBRARIES}
                 )
	message(STATUS "Shared library PROJECT_NAME: ${PROJECT_NAME}")
	target_compile_definitions("${PROJECT_NAME}" PUBLIC -DLIBFLAME_IS_BUILDING_LIBRARY)
	target_compile_definitions("${PROJECT_NAME}" PUBLIC -DLIBFLAME_ENABLE_SHARED)
else ()

	add_library("${PROJECT_NAME}" STATIC ${FLAME_ROOT}/${BASE_INC_PATH}/blis1.h
                 ${FLAME_ROOT}/${BASE_INC_PATH}/FLA_f2c.h
               	 ${FLAME_ROOT}/${BASE_INC_PATH}/FLAME.h
                 ${OBJECT_LIBRARIES}
		)

	message(STATUS "Static library PROJECT_NAME: ${PROJECT_NAME}")
endif ()

# AOCL Build-It-Yourself unified build: contribute libflame's object files to the
# single libaocl. OBJECT_LIBRARIES is a list of $<TARGET_OBJECTS:...> genexes.
if(AOCL_TB_UNIFIED_BUILD AND COMMAND aocl_tb_add_objects)
	aocl_tb_add_objects(${OBJECT_LIBRARIES})
endif()

if(WIN32)
	set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "AOCL-LAPACK lib/ include/ installation path")
endif()
message(STATUS "CMAKE_INSTALL_PREFIX : ${CMAKE_INSTALL_PREFIX}")

set(LIBFLAME_PUBLIC_HEADERS	"include/FLAME.h"
				"${FLAME_ROOT}/${BASE_INC_PATH}/lapacke.h"
				"${FLAME_ROOT}/${BASE_INC_PATH}/libflame_interface.hh")

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${LIBFLAME_PUBLIC_HEADERS}")

#Strip DWARF debug info from the release shared library.
if(BUILD_SHARED_LIBS AND UNIX AND CMAKE_STRIP
   AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug"
   AND NOT CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
	add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
		COMMAND ${CMAKE_STRIP} --strip-debug $<TARGET_FILE:${PROJECT_NAME}>)
endif()

install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
								ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
								RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
								PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

#Logic to create AOCL-libFLAME library name according to the linking AOCL-BLIS library name.
if(WIN32)
	if(EXT_BLAS_LIBNAME)
		string(REPLACE "LibBlis" "LibFlame" TEMP_VAR ${EXT_BLAS_LIBNAME})
		string(REPLACE ".lib" "" TEMP_VAR ${TEMP_VAR})
		set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${TEMP_VAR}")
		message(STATUS "Library Name : ${TEMP_VAR}${CMAKE_SHARED_LIBRARY_SUFFIX}")
	else(EXT_BLAS_LIBNAME)
		message(FATAL_ERROR "BLAS Library name not provided! Set BLAS library name using flag EXT_BLAS_LIBNAME")
	endif()
endif()

if(ENABLE_OPENMP)
	if(EXT_OPENMP_PATH)
		if(NOT EXISTS ${EXT_OPENMP_PATH}/${EXT_OPENMP_LIB})
			message(FATAL_ERROR "\n Invalid path to openmp Library \n \
		${EXT_OPENMP_PATH}/${EXT_OPENMP_LIB} \
		Does not exist.\n Please check the path again")
		endif()

		find_library(LIBIOMP ${EXT_OPENMP_LIB} ${EXT_OPENMP_PATH} NO_DEFAULT_PATH)
        	message("User provided OpenMP library path ${LIBIOMP}")
		
		target_link_libraries(${PROJECT_NAME}  ${LIBIOMP})

		string(APPEND LIBS_LINK_FLAGS "-L${EXT_OPENMP_PATH} -l${EXT_OPENMP_LIB} ")

		if(UNIX)
				set	(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
	        	target_link_libraries(${PROJECT_NAME} -Wl,--as-needed -ldl)
		endif()
	else()
		target_link_libraries(${PROJECT_NAME}  OpenMP::OpenMP_CXX)
		string(APPEND LIBS_LINK_FLAGS "-fopenmp ")
	endif()
endif ()

#Set variable for static library path string
if(NOT BUILD_SHARED_LIBS)
	if(UNIX)
		set(LIBFLAME_STATICLIB_NAME lib${PROJECT_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
		set(LIBFLAME_STATICLIB_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${LIBFLAME_STATICLIB_NAME})
		message(STATUS "libflame static library path ${LIBFLAME_STATICLIB_PATH}")
	else ()
		set(LIBFLAME_STATICLIB_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/$<CONFIG>/${TEMP_VAR}${CMAKE_STATIC_LIBRARY_SUFFIX})
		message(STATUS "libflame static library path ${LIBFLAME_STATICLIB_PATH}")
	endif ()
endif()

if (ENABLE_AOCL_BLAS)
  	add_definitions(-DFLA_ENABLE_AOCL_BLAS)
  	set(FLA_ENABLE_AOCL_BLAS 1)
	message(STATUS "Enabled FLA_ENABLE_AOCL_BLAS flag")
  	target_link_libraries(${PROJECT_NAME} ${BLAS_LIBRARY})
endif ()

	target_include_directories(${PROJECT_NAME} PUBLIC src/lapack/x86/avx2)
	target_include_directories(${PROJECT_NAME} PUBLIC src/lapack/x86/avx512)
	target_include_directories(${PROJECT_NAME} PUBLIC src/lapack/x86/front)
	target_include_directories(${PROJECT_NAME} PUBLIC src/lapack/x86/generic)

# IF AOCL-UTILS paths is not provided by user, try to get it from AOCL_ROOT environment vriable
if (NOT LIBAOCLUTILS_INCLUDE_PATH OR NOT LIBAOCLUTILS_LIBRARY_PATH)
	if (NOT LIBAOCLUTILS_INCLUDE_PATH)
		if(AOCL_UTILS_LIBRARY)
			message("LIBAOCLUTILS_INCLUDE_PATH not set. AOCL-UTILS library path provided by AOCL_ROOT environment variable will be used")
			set(LIBAOCLUTILS_INCLUDE_PATH ${AOCL_UTILS_INCLUDE_DIR})
			set(LIBAOCLUTILS_LIBRARY_PATH ${AOCL_UTILS_LIBRARY})
		else()
			message(FATAL_ERROR "AOCL-UTILS library header path not provided! Set the same using flag LIBAOCLUTILS_INCLUDE_PATH or AOCL_ROOT environment variable")
		endif()
	elseif(NOT LIBAOCLUTILS_LIBRARY_PATH AND (BUILD_SHARED_LIBS OR BUILD_TEST OR BUILD_LEGACY_TEST OR BUILD_NETLIB_TEST))
		if (AOCL_UTILS_LIBRARY)
			message("LIBAOCLUTILS_LIBRARY_PATH not set. AOCL-UTILS library path provided by AOCL_ROOT environment variable will be used")
			set(LIBAOCLUTILS_LIBRARY_PATH ${AOCL_UTILS_LIBRARY})
			set(LIBAOCLUTILS_INCLUDE_PATH ${AOCL_UTILS_INCLUDE_DIR})
		else()
			message(FATAL_ERROR "AOCL-UTILS library path not provided! Set the same using flag LIBAOCLUTILS_LIBRARY_PATH or AOCL_ROOT environment variable")
		endif()
	endif()
endif()

message(STATUS "AOCL-UTILS library path : ${LIBAOCLUTILS_LIBRARY_PATH}")
message(STATUS "AOCL-UTILS library headers path : ${LIBAOCLUTILS_INCLUDE_PATH}")

#If embedding AOCL-UTILS library in AOCL-LAPACK option is selected
if(ENABLE_EMBED_AOCLUTILS)
  #AOCLUtils library path handling
  #If user has provided the library path, use it.
  if(LIBAOCLUTILS_LIBRARY_PATH AND LIBAOCLUTILS_INCLUDE_PATH)

	message(STATUS "User provided AOCLUtils library path : " ${LIBAOCLUTILS_LIBRARY_PATH})
	message(STATUS "User provided AOCLUtils library headers path : " ${LIBAOCLUTILS_INCLUDE_PATH})

	target_include_directories(${PROJECT_NAME} PUBLIC ${LIBAOCLUTILS_INCLUDE_PATH})

  else(LIBAOCLUTILS_LIBRARY_PATH AND LIBAOCLUTILS_INCLUDE_PATH)
	#Else build from source if EMBED_AOCLUTILS is selected
	message(STATUS "LIBAOCLUTILS_LIBRARY_PATH and LIBAOCLUTILS_INCLUDE_PATH flags not provded")

	set(LIBAOCLUTILS_SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/${LIBAOCLUTILS_PROJECT}-prefix/src/${LIBAOCLUTILS_PROJECT}")
	set(LIBAOCLUTILS_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/${LIBAOCLUTILS_PROJECT}-prefix/src/${LIBAOCLUTILS_PROJECT}-build/${LIBAOCLUTILS_INSTALL_DIR_NAME}")
        if(UNIX)
		set(CHECK_LIB_EXISTS ${LIBAOCLUTILS_INSTALL_DIR})
	elseif(WIN32)
		set(CHECK_LIB_EXISTS "${LIBAOCLUTILS_INSTALL_DIR}/lib/*.lib")
	endif ()
	message(STATUS "Check AOCL Utils path ${LIBAOCLUTILS_SOURCE_DIR} exists..")
	if (EXISTS ${CHECK_LIB_EXISTS})
		message(STATUS "Found AOCL Utils source directory. Won't rebuild it")
	else()
		message(STATUS "AOCLUtils library neither provided by user nor earlier build exists. Downloading from GitHub repository")
		message(STATUS "libaoclutils  will be cloned from URL ${LIBAOCLUTILS_GIT_URL} and tag ${LIBAOCLUTILS_GIT_TAG}")
		#Download AOCLUtils library from Github
		include(ExternalProject)
		ExternalProject_Add(
			${LIBAOCLUTILS_PROJECT}
			GIT_REPOSITORY ${LIBAOCLUTILS_GIT_URL}
			GIT_TAG ${LIBAOCLUTILS_GIT_TAG}
			CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${LIBAOCLUTILS_INSTALL_DIR_NAME}
			UPDATE_DISCONNECTED TRUE
			LOG_DOWNLOAD TRUE
			LOG_CONFIGURE TRUE
			LOG_BUILD TRUE
		)

		ExternalProject_Get_Property(${LIBAOCLUTILS_PROJECT} source_dir)
		ExternalProject_Get_Property(${LIBAOCLUTILS_PROJECT} binary_dir)

		set(LIBAOCLUTILS_SOURCE_DIR ${source_dir})
		set(LIBAOCLUTILS_INSTALL_DIR ${binary_dir}/${LIBAOCLUTILS_INSTALL_DIR_NAME})

		message(STATUS "libaoclutils source directory ${LIBAOCLUTILS_SOURCE_DIR}")
		message(STATUS "libaoclutils install directory ${LIBAOCLUTILS_INSTALL_DIR}")

		add_dependencies(${PROJECT_NAME} ${LIBAOCLUTILS_PROJECT})
	endif()

	set(LIBAOCLUTILS_LIBRARY_PATH ${LIBAOCLUTILS_INSTALL_DIR}/lib/${LIBAOCLUTILS_STATICLIB})
	target_include_directories(${PROJECT_NAME} PUBLIC ${LIBAOCLUTILS_SOURCE_DIR}/include)

  endif()

  #For libflame shared library, Link directly with libaoclutils library
  #For libflame static library, merge libaoclutils and libflame binaries into a single static library
  if(BUILD_SHARED_LIBS)
	target_link_libraries(${PROJECT_NAME} ${LIBAOCLUTILS_LIBRARY_PATH})
  else ()
	if(UNIX)
		#Merge libflame and libaoclutils libraries in post build step using archive command
		add_custom_command(TARGET ${PROJECT_NAME}
				POST_BUILD
				COMMAND ar -x ${LIBAOCLUTILS_LIBRARY_PATH}
				COMMAND ar -x ${LIBFLAME_STATICLIB_PATH}
				COMMAND ar -crs ${LIBFLAME_STATICLIB_NAME} *.o
				COMMAND rm *.o
				COMMAND mv ${LIBFLAME_STATICLIB_NAME} ${LIBFLAME_STATICLIB_PATH}
				COMMENT "Merging libflame and libaoclutils library object files into single binary"
				)
	elseif(WIN32)
		#Merge libaoclutils static library with libflame static library in post build step
		add_custom_command(TARGET ${PROJECT_NAME}
				POST_BUILD
				COMMAND lib.exe /OUT:${LIBFLAME_STATICLIB_PATH} ${LIBFLAME_STATICLIB_PATH} ${LIBAOCLUTILS_LIBRARY_PATH}
				COMMENT "Merging libflame and libaoclutils libraries"
				)
	endif ()
  endif()

  if(UNIX)
	#Linking with libstdc++ library is needed for libaoclutils library dependency
	target_link_libraries(${PROJECT_NAME} "libstdc++.so")
  endif ()

elseif(BUILD_SHARED_LIBS)
	if(LIBAOCLUTILS_LIBRARY_PATH AND LIBAOCLUTILS_INCLUDE_PATH)
		target_include_directories(${PROJECT_NAME} PUBLIC ${LIBAOCLUTILS_INCLUDE_PATH})
		target_link_libraries(${PROJECT_NAME} ${LIBAOCLUTILS_LIBRARY_PATH})
	else()
		message(FATAL_ERROR "Shared Library build requires AOCL-Utils library and header path to be set! Please set the same using LIBAOCLUTILS_LIBRARY_PATH and LIBAOCLUTILS_INCLUDE_PATH options respectively")

	endif()

elseif(LIBAOCLUTILS_INCLUDE_PATH)

	message(STATUS "User provided AOCLUtils library headers path : " ${LIBAOCLUTILS_INCLUDE_PATH})

	target_include_directories(${PROJECT_NAME} PUBLIC ${LIBAOCLUTILS_INCLUDE_PATH})

else()
	message(FATAL_ERROR "Header file path of aocl-utils library not set! Please set the same using LIBAOCLUTILS_INCLUDE_PATH option")
endif()

# TODO : enable this flag to build source supplied blas library
# link externally built blas library
if(WIN32)
	if (NOT ENABLE_BUILTIN_BLAS)
		target_link_libraries(${PROJECT_NAME} ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLAS_LIBNAME})
		target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH})
	endif ()
endif()
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
add_definitions(-DEXPMODULE)
add_subdirectory(src)

# AOCL Build-It-Yourself unified build: the f2c LAPACK translation units (which
# define ilaver, sbbcsd, aocl_lapack_dlaset, ...) are attached to the main
# libflame target via target_sources() inside add_subdirectory(src) above, so
# they compile straight into the STATIC/SHARED libflame target rather than into a
# standalone OBJECT library. The OBJECT_LIBRARIES list registered earlier only
# covers the arch kernels (FLA_LAPACK_AVX2/AVX512) and LAPACKE; recompile the
# target's own sources into an OBJECT twin so those f2c symbols are contributed to
# the single libaocl too. Must run AFTER add_subdirectory(src) so SOURCES is
# complete, and in this (target-defining) directory scope so the twin inherits the
# same include directories and compile definitions.
if(AOCL_TB_UNIFIED_BUILD AND COMMAND aocl_tb_objectify)
	aocl_tb_objectify(${PROJECT_NAME})
endif()

if (BUILD_TEST OR BUILD_LEGACY_TEST OR BUILD_NETLIB_TEST)

	if (NOT CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH OR NOT EXT_BLAS_LIBNAME)
		message(FATAL_ERROR "BLAS library path not set! Please set the same using CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH option or AOCL_ROOT environment variable or make sure blis pkg-config file is in PKG_CONFIG_PATH")
	endif()

	if(BUILD_TEST AND NOT ENABLE_AOCL_BLAS AND NOT BLAS_HEADER_PATH)
		message(FATAL_ERROR "BLAS header path not set! Please set the same using BLAS_HEADER_PATH option")
	endif()

	# print paths of external libraries
	message(STATUS "External BLAS library path for test suite: ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}")
	message(STATUS "External BLAS library name for test suite: ${EXT_BLAS_LIBNAME}")
	message(STATUS "External BLAS header path for test suite: ${BLAS_HEADER_PATH}")
	
	# if user wants to run tests on an external lapack library
	# needs to be modified for external lapack libraries same way as done for blas
	if((NOT EXT_LAPACK_LIBRARY_PATH) OR (NOT EXT_LAPACK_LIBNAME))

		if (UNIX)
			if (BUILD_SHARED_LIBS)
				set(LIBFLAME_LIB_NAME lib${PROJECT_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX})
				set(LIBFLAME_LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
			else()
				set(LIBFLAME_LIB_NAME lib${PROJECT_NAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
				set(LIBFLAME_LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
			endif()
		else ()
			if (BUILD_SHARED_LIBS)
				set(LIBFLAME_LIB_NAME ${TEMP_VAR}${CMAKE_SHARED_LIBRARY_SUFFIX})
				set(LIBFLAME_LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/$<CONFIG>/})
			else()
				set(LIBFLAME_LIB_NAME ${TEMP_VAR}${CMAKE_STATIC_LIBRARY_SUFFIX})
				set(LIBFLAME_LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/$<CONFIG>/)
			endif()
		endif()

		set(EXT_LAPACK_LIBRARY_PATH ${LIBFLAME_LIBRARY_OUTPUT_PATH})
		set(EXT_LAPACK_LIBNAME ${LIBFLAME_LIB_NAME})
		message(STATUS "Missing External lapack library. Using ${LIBFLAME_LIBRARY_OUTPUT_PATH}/${LIBFLAME_LIB_NAME} ")
	elseif(NOT EXISTS ${EXT_LAPACK_LIBRARY_PATH}/${EXT_LAPACK_LIBNAME})
		message(FATAL_ERROR "\n Invalid path to lapack Library \n \
		${EXT_LAPACK_LIBRARY_PATH}/${EXT_LAPACK_LIBNAME} \
		Does not exist.\n Please check the path again, or avoid passing \n \
		-DEXT_LAPACK_LIBRARY_PATH and -DEXT_LAPACK_LIBNAME to use libflame")
	endif()

endif()

if (BUILD_TEST OR BUILD_LEGACY_TEST)
	# enable ctest
	enable_testing()
	message(STATUS "CTEST ENABLED")

	add_library(blas SHARED IMPORTED)
	add_library(aoclutils SHARED IMPORTED)
	if(NOT ENABLE_EMBED_AOCLUTILS AND NOT EXISTS ${LIBAOCLUTILS_LIBRARY_PATH})
		message(FATAL_ERROR "aocl-utils library path not set! Please set the same using LIBAOCLUTILS_LIBRARY_PATH \
		option or make sure aocl-utils pkg-config file is in PKG_CONFIG_PATH")
	endif()
		
	if(UNIX)
		# CHECKS if blas library path is explicitly set
		if (CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH AND EXT_BLAS_LIBNAME AND (NOT ENABLE_AOCL_BLAS))
			# CHECKS if blas library exists
			set(BLAS_LIB_LIST ${EXT_BLAS_LIBNAME})

			# If there are multiple libraries that are part of given BLAS, each will be separated by space
			separate_arguments(BLAS_LIB_LIST)

			foreach(BLAS_LIB ${BLAS_LIB_LIST})
				if(NOT EXISTS ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${BLAS_LIB})
					message(FATAL_ERROR "\n Invalid path to BLAS Library \n \
					${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLAS_LIBNAME} \
					Does not exist.\n Please check the path again")
				endif()
			endforeach()

			foreach(BLAS_LIB ${BLAS_LIB_LIST})
				list(APPEND BLAS_LINK_LIST ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${BLAS_LIB})
			endforeach()
	
		elseif(ENABLE_AOCL_BLAS)
			# If BLAS library is provided from environment, then use it
			list(APPEND BLAS_LINK_LIST ${BLAS_LIBRARY})
			set(BLAS_HEADER_PATH ${BLAS_INCLUDE_DIR})			
			message(STATUS "BLAS library from environment")
		endif()

		set_target_properties(aoclutils PROPERTIES
		IMPORTED_LOCATION ${LIBAOCLUTILS_LIBRARY_PATH}
		POSITION_INDEPENDENT_CODE ON
		)

	elseif(WIN32)
		# CHECKS if blas library is explicitly set
		if (CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH AND EXT_BLAS_LIBNAME AND (NOT ENABLE_AOCL_BLAS))
			if(NOT EXISTS ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLAS_LIBNAME})
				message(FATAL_ERROR "\n Invalid path to BLAS Library \n \
					${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLAS_LIBNAME} \
					Does not exist.\n Please check the path again")
			endif()

			set_target_properties(blas PROPERTIES
			IMPORTED_IMPLIB ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLAS_LIBNAME}
			POSITION_INDEPENDENT_CODE ON
			)
		elseif (ENABLE_AOCL_BLAS)
			# If BLAS library is provided from environment, then use it
			set_target_properties(blas PROPERTIES
			IMPORTED_IMPLIB ${BLAS_LIBRARY}
			POSITION_INDEPENDENT_CODE ON
			)
		endif()

		set_target_properties(aoclutils PROPERTIES
		IMPORTED_IMPLIB ${LIBAOCLUTILS_LIBRARY_PATH}
		POSITION_INDEPENDENT_CODE ON
		)
	endif()
endif()

if(BUILD_TEST)
	add_subdirectory(test/main)
endif()
if(BUILD_LEGACY_TEST)
	add_subdirectory(test/legacyflame)
endif()

if (BUILD_NETLIB_TEST)

	if (BUILD_SHARED_LIBS)
		message(FATAL_ERROR "Netlib test is not supported for shared library build")
	endif()

	if(WIN32)
		add_subdirectory(netlib-test)
	else()
		include(${CMAKE_CURRENT_SOURCE_DIR}/netlib-test/netlib_ctest.cmake)
	endif()
endif()

if(BUILD_LAPACKE_TEST)
	add_subdirectory(test/lapacke)
endif()


set(REQ_PKG_LIB ${PROJECT_NAME})
configure_file(aocl-lapack.pc.in flame.pc @ONLY)
# configure_file writes flame.pc into THIS project's binary dir
# (CMAKE_CURRENT_BINARY_DIR). Under the AOCL-BIY in-tree build CMAKE_BINARY_DIR is
# the parent (aocl) tree, so install must read from CMAKE_CURRENT_BINARY_DIR; the
# two coincide for the standalone build.
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/flame.pc DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)

# To View Added Definitions and C Flags
get_directory_property( DirDefs COMPILE_DEFINITIONS )
message(STATUS "COMPILE_DEFINITIONS : ${DirDefs} " )
message(STATUS "CMAKE_C_FLAGS : ${CMAKE_C_FLAGS}")
# make clean for additional header files any further files can be added in the list
set_property(
        TARGET ${PROJECT_NAME}
        APPEND
        PROPERTY ADDITIONAL_CLEAN_FILES
		${FLAME_ROOT}/build/FLA_config.h ${FLAME_ROOT}/include/
)
if (BUILD_DOC)
    add_subdirectory(docs/libflame)
endif()
