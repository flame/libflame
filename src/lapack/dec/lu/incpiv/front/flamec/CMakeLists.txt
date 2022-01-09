##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_incpiv.h
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_FS_incpiv.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_LU_incpiv.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_LU_incpiv_create_hier_matrices.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_LU_incpiv_noopt.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_LU_incpiv_opt1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_LU_incpiv_solve.c
	)