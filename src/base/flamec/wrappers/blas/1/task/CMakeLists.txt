##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Axpy_task.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Axpyt_task.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Copy_task.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Copyr_task.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Copyt_task.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Scal_task.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Scalr_task.c
	)