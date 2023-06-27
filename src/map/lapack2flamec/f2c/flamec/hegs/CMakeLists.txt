##Copyright (C) 2023, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/chegs2_fla.c
${CMAKE_CURRENT_SOURCE_DIR}/chegst_fla.c
${CMAKE_CURRENT_SOURCE_DIR}/zhegs2_fla.c
${CMAKE_CURRENT_SOURCE_DIR}/zhegst_fla.c
	)