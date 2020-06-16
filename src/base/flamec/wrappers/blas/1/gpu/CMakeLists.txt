##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Axpy_external_gpu.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Copy_external_gpu.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Scal_external_gpu.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Scalr_external_gpu.c
	)