##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_Q_UT_rnbc.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_Q_UT_rnbc_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_Q_UT_rnbc_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_Q_UT_rnbc_blk_var3.c
	)