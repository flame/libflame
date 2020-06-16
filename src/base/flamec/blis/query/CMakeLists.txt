##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/bl1_check.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_does.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_is.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_proj.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_vector.c
	)