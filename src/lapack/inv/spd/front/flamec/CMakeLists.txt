##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_SPDinv.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_SPDinv.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_SPDinv_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_SPDinv.c
	)