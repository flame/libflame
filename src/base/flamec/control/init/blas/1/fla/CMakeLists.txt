##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Axpy_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Axpyt_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Copy_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Copyr_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Copyt_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Scal_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Scalr_cntl_init.c
	)