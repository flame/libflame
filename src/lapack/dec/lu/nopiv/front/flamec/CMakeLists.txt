##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_nopiv.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_nopiv.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_nopiv_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_nopiv_solve.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_LU_nopiv.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_LU_nopiv_solve.c
	)