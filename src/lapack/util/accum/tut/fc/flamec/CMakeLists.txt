##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Accum_T_UT_fc.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Accum_T_UT_fc_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Accum_T_UT_fc_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Accum_T_UT_fc_unb_var1.c
	)