##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/bl1_abort.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_param_map.c
	)