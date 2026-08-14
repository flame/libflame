##Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.##
	target_sources("${PROJECT_NAME}" PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_lu_small_kernels_common.h
	${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_lu_small_kernels_d.h
	${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_lu_small_kernels_s.h
	${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_lu_small_kernels_c.h
	${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_lu_getri_small_d.c
	${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_cholesky_small_kernels.h
	)
