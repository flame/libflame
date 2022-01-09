##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_ruh.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_ruh_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_ruh_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_ruh_blk_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_ruh_blk_var4.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_ruh_unb_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_ruh_unb_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_ruh_unb_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_ruh_unb_var4.c
	)