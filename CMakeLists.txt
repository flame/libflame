##Copyright (C) 2020, Advanced Micro Devices, Inc.##

cmake_minimum_required(VERSION 3.0.0)

project(AOCL-LibFLAME-Win)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
	option(BUILD_SHARED_LIBS "Build shared library" ON)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MP")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W0 ")
	add_definitions(-DFLA_ENABLE_WINDOWS_BUILD)	
	
include_directories(${CMAKE_SOURCE_DIR}/windows/include/.)
include_directories(${CMAKE_SOURCE_DIR}/src/map/lapack2flamec)
include_directories(${CMAKE_SOURCE_DIR}/src/base/flamec/blis/include)
file(GLOB headers ${CMAKE_SOURCE_DIR}/*.h)
if(BUILD_SHARED_LIBS)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MD ")
add_library("${PROJECT_NAME}" SHARED ${CMAKE_SOURCE_DIR}/windows/include/blis1.h
                                     ${CMAKE_SOURCE_DIR}/windows/include/FLA_f2c.h
						             ${CMAKE_SOURCE_DIR}/windows/include/FLAME.h
                                 )
endif()								 
if(NOT BUILD_SHARED_LIBS)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /MT ")
add_library("${PROJECT_NAME}" STATIC ${CMAKE_SOURCE_DIR}/windows/include/blis1.h
                                     ${CMAKE_SOURCE_DIR}/windows/include/FLA_f2c.h
						             ${CMAKE_SOURCE_DIR}/windows/include/FLAME.h
                                 )
endif()


link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
add_definitions(-DEXPMODULE)
add_subdirectory(src)
add_subdirectory(test)
#add_subdirectory(testcpp)