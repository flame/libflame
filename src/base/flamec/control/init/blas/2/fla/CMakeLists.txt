##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemv_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsv_cntl_init.c
	)