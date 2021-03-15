##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_Q_UT_inc.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_Q_UT_inc_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_Q_UT_inc_lhfc.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_Q_UT_inc_lnfc.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Apply_Q_UT_inc.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Apply_Q_UT_inc_create_workspace.c
	)