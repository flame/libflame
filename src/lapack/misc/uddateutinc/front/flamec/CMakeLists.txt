##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_UDdate_UT_inc.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_UDdate_UT_inc_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_UDdate_UT_inc.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_UDdate_UT_inc_create_hier_matrices.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_UDdate_UT_inc_solve.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_UDdate_UT_inc_update_rhs.c
	)