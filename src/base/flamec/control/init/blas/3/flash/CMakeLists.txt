##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Gemm_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Hemm_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Her2k_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Herk_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Symm_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Syr2k_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Syrk_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Trmm_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Trsm_cntl_init.c
	)