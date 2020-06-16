##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Cntl_init.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Cntl_init_flamec.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Cntl_init_flamec.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Cntl_init_flash.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Cntl_init_flash.h
	)
	
add_subdirectory(base)
add_subdirectory(blas)
add_subdirectory(lapack)