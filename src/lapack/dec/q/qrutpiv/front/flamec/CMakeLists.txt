##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR_UT_piv.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR_UT_piv.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR_UT_piv_colnorm.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_QR_UT_piv_internal.c
	)