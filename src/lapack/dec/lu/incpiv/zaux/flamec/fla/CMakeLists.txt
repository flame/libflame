##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_SA_Apply_pivots.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_SA_FS_blk.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_SA_LU_blk.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_SA_LU_unb.c
	)