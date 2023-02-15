##Copyright (C) 2023, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_fblas_kernels.h
${CMAKE_CURRENT_SOURCE_DIR}/fla_idamax_small.c
	)
