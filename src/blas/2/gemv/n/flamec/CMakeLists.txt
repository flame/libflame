##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemv_n.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemv_n_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemv_n_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemv_n_blk_var5.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemv_n_blk_var6.c
	)