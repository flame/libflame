##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Gemm_external.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Hemm_external.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Her2k_external.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Herk_external.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Symm_external.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Syr2k_external.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Syrk_external.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trmm_external.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trmmsx_external.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsm_external.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Trsmsx_external.c
	)