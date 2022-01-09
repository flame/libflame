##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Hemm_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Her2k_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Herk_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Symm_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Syr2k_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Syrk_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trmm_cntl_init.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_cntl_init.c
	)