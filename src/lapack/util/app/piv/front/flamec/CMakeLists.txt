##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_pivots.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_pivots.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_pivots_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_pivots_ln.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_pivots_lt.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_pivots_rn.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_pivots_rt.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Apply_pivots.c
	)