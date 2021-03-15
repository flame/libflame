##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_CAQ2_UT.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_CAQ2_UT_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_CAQ2_UT_lhfc.c
	)