##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_CAQR_UT_inc.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_CAQR_UT_inc_blk_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_CAQR_UT_inc_compute_p_length.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_CAQR_UT_inc_copy_triangles.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_CAQR_UT_inc_factorize_panels.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_CAQR_UT_inc_init_structure.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_CAQR_UT_inc.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_CAQR_UT_inc_create_hier_matrices.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_CAQR_UT_inc_noopt.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_CAQR_UT_inc_solve.c
	)