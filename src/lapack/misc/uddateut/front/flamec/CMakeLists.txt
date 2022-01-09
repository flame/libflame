##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_UDdate_UT.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_UDdate_UT.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_UDdate_UT_create_T.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_UDdate_UT_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_UDdate_UT_solve.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_UDdate_UT_update_rhs.c
	)