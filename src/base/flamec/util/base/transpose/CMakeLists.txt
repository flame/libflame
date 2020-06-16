##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Swap_t_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Swap_t_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Transpose.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Transpose_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Transpose_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Transpose_unb_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Transpose_unb_var2.c
	)