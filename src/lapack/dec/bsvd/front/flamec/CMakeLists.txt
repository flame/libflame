##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd.c
	${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_create_workspace.c
	${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_ext.c
	)