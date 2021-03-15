##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_blk_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_blk_var4.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_blk_var5.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_opt_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_opt_var4.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_opt_var5.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_unb_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_unb_var3b.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_unb_var4.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_unb_var5.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_vars.h
	)