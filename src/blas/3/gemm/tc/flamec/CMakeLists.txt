##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_blk_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_blk_var4.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_blk_var5.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_blk_var6.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_unb_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_unb_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_unb_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_unb_var4.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_unb_var5.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_tc_unb_var6.c
	)