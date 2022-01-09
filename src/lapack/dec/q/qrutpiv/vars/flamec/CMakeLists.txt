##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_H2_UT_piv_row.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR_UT_piv_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR_UT_piv_blk_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR_UT_piv_unb_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR_UT_piv_unb_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR_UT_piv_vars.h
	)