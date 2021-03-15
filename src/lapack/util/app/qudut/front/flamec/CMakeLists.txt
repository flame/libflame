##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_QUD_UT.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_QUD_UT.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_QUD_UT_create_workspace.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_QUD_UT_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_QUD_UT_lhfc.c
	)