##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Accum_T_UT_fr.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Accum_T_UT_fr_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Accum_T_UT_fr_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Accum_T_UT_fr_unb_var1.c
	)