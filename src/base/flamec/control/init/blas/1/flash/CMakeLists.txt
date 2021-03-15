##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Axpy_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Axpyt_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Copy_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Copyr_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Copyt_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Scal_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Scalr_cntl_init.c
	)