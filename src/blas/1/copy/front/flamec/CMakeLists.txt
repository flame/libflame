##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Copy.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Copy.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Copy_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Copy.c
	)