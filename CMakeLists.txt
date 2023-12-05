##Copyright (C) 2023, Advanced Micro Devices, Inc.##

cmake_minimum_required(VERSION 3.15.0)

if(WIN32)
	project(AOCL-LibFLAME-Win)
elseif(UNIX)
	project(flame)
else()
	message( FATAL_ERROR "OS Undetected")
endif()

if(WIN32)
	SET(CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH "${CMAKE_SOURCE_DIR}/bin/Release" CACHE STRING "AOCL-BLIS library path name")
	SET(EXT_BLAS_LIBNAME "AOCL-LibBlis-Win-MT-dll.lib" CACHE STRING "AOCL-BLIS library name")
	SET(OpenMP_libomp_LIB "C:/Program Files/LLVM/lib/libomp.lib" CACHE STRING "openmp library path")
	SET(LIBAOCLUTILS_LIBRARY_PATH "" CACHE STRING "AOCLUtils library name")
	SET(LIBAOCLUTILS_INCLUDE_PATH "" CACHE STRING "AOCLUtils library header files path")
endif()

message(STATUS "External blas library Path : ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}")
message(STATUS "External blas library name : ${EXT_BLAS_LIBNAME}")
SET(LF_ISA_CONFIG "AVX2" CACHE STRING "LibFlame Instruction set architecture Configuration")

#gccopenmp
if(WIN32)
	SET(OpenMP_libomp_LIBRARY "${OpenMP_libomp_LIB}" CACHE STRING "openmp library path")
endif()

message(STATUS "CMAKE_C_COMPILER_ID: ${CMAKE_C_COMPILER_ID}")
message(STATUS "C compiler version: ${CMAKE_C_COMPILER_VERSION}")


#set the output directory path for windows and Linux
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
# Added for debugging purpose
set(CMAKE_VERBOSE_MAKEFILE OFF)

if(WIN32)
	set (PYTHON_EXE "python")
else()
	set (PYTHON_EXE "python3")
endif()
### Making LF_ISA_CONFIG case insensitive ###
string(TOLOWER "${LF_ISA_CONFIG}" lf_isa_config_lower)

if ("${lf_isa_config_lower}" STREQUAL "none" OR "${lf_isa_config_lower}" STREQUAL "avx" OR "${lf_isa_config_lower}" STREQUAL "avx2" OR "${lf_isa_config_lower}" STREQUAL "avx512" OR "${lf_isa_config_lower}" STREQUAL "auto")
   message("Machine ISA configuration set to ${LF_ISA_CONFIG}")
else()
   message(FATAL_ERROR "Machine ISA configuration error, valid values are AVX, AVX2, AVX512 and auto")
endif()

if (${lf_isa_config_lower} STREQUAL "auto")
    set(AUTO_CONFIG_PY "${CMAKE_SOURCE_DIR}/build/auto_config.py")
    # Run python script to find the architecture family name
    execute_process(
          COMMAND ${PYTHON_EXE} ${AUTO_CONFIG_PY}
          RESULT_VARIABLE CMD_RESULT
          OUTPUT_VARIABLE CMD_OUTPUT
          OUTPUT_STRIP_TRAILING_WHITESPACE)
    message( STATUS "Auto configuring the family :" ${CMD_OUTPUT})
    if(${CMD_OUTPUT} STREQUAL "zen" OR ${CMD_OUTPUT} STREQUAL "zen2" OR ${CMD_OUTPUT} STREQUAL "zen3")
      set(lf_isa_config_lower "avx2") 
    elseif(${CMD_OUTPUT} STREQUAL "zen4")
      set(lf_isa_config_lower "avx512")
    else()
      set(lf_isa_config_lower "none")
    endif()
endif ()
message("LF_ISA_CONFIG selected:${LF_ISA_CONFIG}")
#set AOCL Utilities library path values
set(LIBAOCLUTILS_PROJECT "libaoclutils")
set(LIBAOCLUTILS_INSTALL_DIR_NAME "install-aoclutils")
set(LIBAOCLUTILS_GIT_URL "https://github.com/amd/aocl-utils.git" CACHE STRING "libaoclutils repository URL to be used")
set(LIBAOCLUTILS_GIT_TAG "4.1" CACHE STRING "libaoclutils repository tag to be used")

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
# hiding this option form cmake gui
mark_as_advanced(ENABLE_AMD_COMMON_FLAGS)
option (ENABLE_AMD_FLAGS "Enable AMD FLAGS internally enables ENABLE_AMD_COMMON_FLAGS" OFF)
option (ENABLE_AMD_AOCC_FLAGS "Enable AMD FLAGS internally enables Enable AMD FLAGS internally enables ENABLE_AMD_COMMON_FLAGS, ENABLE_F2C_DOTC, ENABLE_VOID_RETURN_COMPLEX_FUNCTION" OFF)
# Since AOCC isn't used in windows hiding this option form cmake gui
if(WIN32)
	mark_as_advanced(ENABLE_AMD_AOCC_FLAGS)
endif()
option (ENABLE_AOCL_BLAS "Enables tight coupling with AOCL-BLAS library in order to use AOCL-BLAS internal routines" OFF)
option (ENABLE_ILP64 "Enable ILP64" OFF)
option (ENABLE_UPPERCASE "Enable uppercase API's" OFF)
option (ENABLE_WRAPPER "Enable wrapper code" ON)
option (INCLUDE_LAPACKE "INCLUDE LAPACKE Code" ON)
option (BUILD_SHARED_LIBS "Building shared Library" ON)
option (BUILD_TEST "Building aoclflaprogress and main test directory" OFF)
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
set_property(CACHE ENABLE_AOCL_DTL PROPERTY STRINGS TRACE_AND_LOG TRACE LOG OFF )
option (ENABLE_EMBED_AOCLUTILS "Enables embedding aocl-utils library in AOCL-LAPACK" OFF)

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

if (ENABLE_AOCL_BLAS)
       message(STATUS "AOCL-BLAS library linking option selected")
       include(cmake/find_aocl-blas.cmake)
endif()


message(STATUS "checking SIMD  ")
if(WIN32)
	if(${lf_isa_config_lower} STREQUAL "avx512")
	add_definitions(/arch:AVX512)
	elseif(${lf_isa_config_lower} STREQUAL "avx2" OR ${lf_isa_config_lower} STREQUAL "avx")
	add_definitions(/arch:AVX2)
	else()
	message(STATUS "No ISA flag set")
	endif()
elseif(UNIX)
	set(COMPILER_OPTIMIZATION_FLAGS	"-mtune=native -O3")
	if(${lf_isa_config_lower} STREQUAL "avx512")
		set(COMPILER_OPTIMIZATION_FLAGS	"${COMPILER_OPTIMIZATION_FLAGS} -mavx512f -mfma")
	elseif(${lf_isa_config_lower} STREQUAL "avx2" OR ${lf_isa_config_lower} STREQUAL "avx")
		set(COMPILER_OPTIMIZATION_FLAGS	"${COMPILER_OPTIMIZATION_FLAGS} -mavx2 -mfma")
	else()
		message(STATUS "No ISA flag set")
	endif()
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMPILER_OPTIMIZATION_FLAGS}")
else()
	message(STATUS "OS UNKNOWN CANNOT SET SIMD")
endif()

if(ENABLE_GCOV)
	set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE ON)
	add_definitions(--coverage)
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
endif()

#Adding PIC flag and other flags
if(UNIX)
	SET(COMPILER_DEBUG_FLAG "-g0")
	set(GCC_WARNING_FLAGS "-Wall -Wno-comment")
	set(COMPILER_LANGUAGE_FLAGS " -std=c11 -D_GNU_SOURCE -Wno-unused-function -Wno-parentheses -Wfatal-errors ")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMPILER_DEBUG_FLAG} ${GCC_WARNING_FLAGS} ${COMPILER_LANGUAGE_FLAGS}")
	add_compile_options(-fPIC)
endif()

#Temporary hack for Clang. Enable -Wno-error=implicit-function-declaration
#and -Wno-error=incompatible-function-pointer-types
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error=implicit-function-declaration -Wno-error=incompatible-function-pointer-types")
endif ()

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
		set(LF_AOCL_DTL_TRACE_ENABLE 1)
		set(LF_AOCL_DTL_LOG_ENABLE 1)
	elseif(ENABLE_AOCL_DTL STREQUAL "TRACE")
		set(LF_AOCL_DTL_TRACE_ENABLE 1)
	elseif(ENABLE_AOCL_DTL STREQUAL "LOG")
		set(LF_AOCL_DTL_LOG_ENABLE 1)
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
    if (ENABLE_AMD_FLAGS OR ENABLE_AMD_AOCC_FLAGS)
	set (FLA_OPENMP_MULTITHREADING TRUE)
	add_definitions(-DFLA_OPENMP_MULTITHREADING)
	message(STATUS "Enabled FLA_OPENMP_MULTITHREADING : ${FLA_OPENMP_MULTITHREADING} ")
    else ()
    	set (FLA_ENABLE_MULTITHREADING TRUE)
	add_definitions(-DFLA_ENABLE_MULTITHREADING)
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

if(WIN32)
	add_definitions(-DBLIS1_ENABLE_WINDOWS_BUILD)
endif()

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
	add_definitions(-DFLA_ENABLE_ILP64=1)
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
	if(UNIX)
	set(FLA_LIBFLAME_VERSION TRUE)
	set( LIBRARY_VERSION "\"${MAJOR_VERSION}.${MINOR_VERSION}\tBuild\t${TODAY}\"")
	elseif(WIN32)
		set( LIBRARY_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}\tBuild\t${TODAY}")
	endif()
	add_definitions(-DFLA_LIBFLAME_VERSION=${LIBRARY_VERSION})
endif()
#generates the FLA_config.h header file
if(WIN32)
	configure_file (cmake.config.h.in ${CMAKE_SOURCE_DIR}/windows/build/FLA_config.h @ONLY)
elseif(UNIX)
# As per the current build system
	configure_file(cmake.config.h.in ${CMAKE_SOURCE_DIR}/build/FLA_config.h @ONLY)
else()
	message(FATAL_ERROR " OS not found")
endif()

if(WIN32)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
endif()

#genreating the monolithic headers
find_package(PythonInterp REQUIRED)
if(WIN32)
	set (MK_HEADER_DIR_PATHS "${CMAKE_SOURCE_DIR}/windows/build")
elseif(UNIX)
	set(MK_HEADER_DIR_PATHS "${CMAKE_SOURCE_DIR}/build")
else()
	message(FATAL_ERROR "OS NOT DETECTED")
endif()
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
	" ${CMAKE_SOURCE_DIR}/src/lapack/fblas/"
	" ${CMAKE_SOURCE_DIR}/src/lapack/fblas/front"
	" ${CMAKE_SOURCE_DIR}/src/lapack/fblas/src"
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
if(ENABLE_AOCL_DTL)
	#only enabling aocl dtl when ENABLE_AOCL_DTL is set to a valid value
	include_directories(${CMAKE_SOURCE_DIR}/src/aocl_dtl)
endif()
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
include_directories(${CMAKE_SOURCE_DIR}/src/lapack/x86/avx512)
include_directories(${CMAKE_SOURCE_DIR}/src/lapack/fblas/front)
include_directories(${CMAKE_SOURCE_DIR}/src/lapack/fblas/src)
if(ENABLE_AOCL_BLAS)
	include_directories(${BLAS_INCLUDE_DIR}) 
endif()

file(GLOB headers ${CMAKE_SOURCE_DIR}/*.h)

if(WIN32)
	cmake_policy(SET CMP0091 NEW)
	if(BUILD_SHARED_LIBS)
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
	else()
		set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
	endif()
endif()

set(OBJECT_LIBRARIES 
    $<TARGET_OBJECTS:AVX512_LU>
    $<TARGET_OBJECTS:AVX2_LU>
)

if (BUILD_SHARED_LIBS)
	if(WIN32)
		set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
	endif()
	add_library("${PROJECT_NAME}" SHARED ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/blis1.h
                 ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/FLA_f2c.h
                 ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/FLAME.h
                 ${OBJECT_LIBRARIES}
                 )
	message(STATUS "Shared library PROJECT_NAME: ${PROJECT_NAME}")
	target_compile_definitions("${PROJECT_NAME}" PUBLIC -DLIBFLAME_IS_BUILDING_LIBRARY)
	target_compile_definitions("${PROJECT_NAME}" PUBLIC -DLIBFLAME_ENABLE_SHARED)
else ()

	add_library("${PROJECT_NAME}" STATIC ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/blis1.h
                 ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/FLA_f2c.h
               	 ${CMAKE_SOURCE_DIR}/${BASE_INC_PATH}/FLAME.h
                 ${OBJECT_LIBRARIES}
		)

	message(STATUS "Static library PROJECT_NAME: ${PROJECT_NAME}")
endif ()

if(WIN32)
	set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "AOCL-LAPACK lib/ include/ installation path")
endif()
message(STATUS "CMAKE_INSTALL_PREFIX : ${CMAKE_INSTALL_PREFIX}")

set(LIBFLAME_PUBLIC_HEADERS	"include/FLAME.h"
				"src/lapacke/LAPACKE/include/lapacke_mangling.h"
				"src/lapacke/LAPACKE/include/lapacke.h"
				"src/lapacke/LAPACKE/include/lapack.h"
				"src/src_cpp/libflame_interface.hh"
				"src/src_cpp/libflame.hh")

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${LIBFLAME_PUBLIC_HEADERS}")

install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
								ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
								RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
								PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

#Logic to create AOCL-libFLAME library name according to the linking AOCL-BLIS library name.
if(WIN32)
	string(REPLACE "LibBlis" "LibFlame" TEMP_VAR ${EXT_BLAS_LIBNAME})
	string(REPLACE ".lib" "" TEMP_VAR ${TEMP_VAR})
	set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${TEMP_VAR}")
	message(STATUS "Library Name : ${TEMP_VAR}${CMAKE_SHARED_LIBRARY_SUFFIX}")
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
		if(UNIX)
	        	target_link_libraries(${PROJECT_NAME} -Wl,--as-needed -ldl)
		endif()
	else()
		target_link_libraries(${PROJECT_NAME}  OpenMP::OpenMP_CXX)
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

elseif(WIN32 AND BUILD_SHARED_LIBS) 
	if(LIBAOCLUTILS_LIBRARY_PATH AND LIBAOCLUTILS_INCLUDE_PATH)
		target_include_directories(${PROJECT_NAME} PUBLIC ${LIBAOCLUTILS_INCLUDE_PATH})	
		target_link_libraries(${PROJECT_NAME} ${LIBAOCLUTILS_LIBRARY_PATH})
	else()
		message(FATAL_ERROR "Windows Shared Library build requires AOCL-Utils library and header path to be set! Please set the same using LIBAOCLUTILS_LIBRARY_PATH and LIBAOCLUTILS_INCLUDE_PATH options respectively")

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

if (BUILD_TEST OR BUILD_LEGACY_TEST)
# enable ctest
enable_testing()
message(STATUS "CTEST ENABLED")
# CHECKS if blas library exists 
	if(NOT EXISTS ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLAS_LIBNAME})
		message(FATAL_ERROR "\n Invalid path to blas Library \n \
${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLAS_LIBNAME} \
Does not exist.\n Please check the path again")
	else()
		add_library(blas SHARED IMPORTED)
		add_library(aoclutils SHARED IMPORTED)
		if(NOT ENABLE_EMBED_AOCLUTILS AND NOT EXISTS ${LIBAOCLUTILS_LIBRARY_PATH})
			message(FATAL_ERROR "aocl-utils library path not set! Please set the same using LIBAOCLUTILS_LIBRARY_PATH option")
		endif()
		if(UNIX)
			set_target_properties(blas PROPERTIES
			IMPORTED_LOCATION ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLAS_LIBNAME}
			POSITION_INDEPENDENT_CODE ON
			)
			set_target_properties(aoclutils PROPERTIES
			IMPORTED_LOCATION ${LIBAOCLUTILS_LIBRARY_PATH}
			POSITION_INDEPENDENT_CODE ON
			)
		elseif(WIN32)
			set_target_properties(blas PROPERTIES
			IMPORTED_IMPLIB ${CMAKE_EXT_BLAS_LIBRARY_DEPENDENCY_PATH}/${EXT_BLAS_LIBNAME}
			POSITION_INDEPENDENT_CODE ON
			)
			set_target_properties(aoclutils PROPERTIES
			IMPORTED_IMPLIB ${LIBAOCLUTILS_LIBRARY_PATH}
			POSITION_INDEPENDENT_CODE ON
			)
		endif()
	endif()

	# if user wants to run tests on an external lapack library
	# needs to be modified for external lapack libraries same way as done for blas
	if((NOT EXT_LAPACK_LIBRARY_PATH) OR (NOT EXT_LAPACK_LIBNAME))
		set(EXT_LAPACK_LIBRARY_PATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
		message(STATUS "Missing External lapack library. Using ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/${EXT_LAPACK_LIBNAME} ")
	elseif(NOT EXISTS ${EXT_LAPACK_LIBRARY_PATH}/${EXT_LAPACK_LIBNAME})
		message(FATAL_ERROR "\n Invalid path to lapack Library \n \
	${EXT_LAPACK_LIBRARY_PATH}/${EXT_LAPACK_LIBNAME} \
	Does not exist.\n Please check the path again, or avoid passing \n \
	-DEXT_LAPACK_LIBRARY_PATH and -DEXT_LAPACK_LIBNAME to use libflame")
	endif()
	#add_subdirectory(testcpp)
endif()

if(BUILD_TEST)
	add_subdirectory(test/main)
	add_subdirectory(test/aoclflaprogress)
endif()
if(BUILD_LEGACY_TEST)
	add_subdirectory(test/legacyflame)
endif()

if (BUILD_NETLIB_TEST)
	if(WIN32)
		add_subdirectory(netlib-test)
	else()
		include(${CMAKE_CURRENT_SOURCE_DIR}/netlib-test/netlib_ctest.cmake)
	endif()
endif()

if(BUILD_LAPACKE_TEST)
	add_subdirectory(test/lapacke)
endif()

# To View Added Definitions and C Flags 
get_directory_property( DirDefs COMPILE_DEFINITIONS ) 
message(STATUS "COMPILE_DEFINITIONS : ${DirDefs} " )
message(STATUS "CMAKE_C_FLAGS : ${CMAKE_C_FLAGS}")
# make clean for additional header files any further files can be added in the list
set_property(
        TARGET ${PROJECT_NAME}
        APPEND
        PROPERTY ADDITIONAL_CLEAN_FILES 
		${CMAKE_SOURCE_DIR}/build/FLA_config.h ${CMAKE_SOURCE_DIR}/include/
)
