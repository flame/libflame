##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemv_task.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsv_task.c
	)