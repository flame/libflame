##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Hess_UT.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Hess_UT.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Hess_UT_create_T.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Hess_UT_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Hess_UT_recover_tau.c
	)