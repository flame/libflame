##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_cvar13.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_cvar15.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_cvar31.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_cvar35.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_cvar51.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_cvar53.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_var2.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_var4.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_var5.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_nn_omp_var6.c
	)