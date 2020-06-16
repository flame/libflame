##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Scal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Scal.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Scal_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Scal.c
	)