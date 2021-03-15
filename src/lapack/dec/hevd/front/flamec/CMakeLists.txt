##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/FLA_Hevd.c
    ${CMAKE_CURRENT_SOURCE_DIR}/FLA_Hevd.h
	${CMAKE_CURRENT_SOURCE_DIR}/FLA_Hevd_compute_scaling.c
	)