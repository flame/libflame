##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Obj_create_buffer_task.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Obj_free_buffer_task.c
	)