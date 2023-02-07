##Copyright (C) 2022, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_avx2_kernels.h
${CMAKE_CURRENT_SOURCE_DIR}/fla_dhrot3_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_drot_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_zscal_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_zgetrf_avx2.c
	)
