##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_create_T.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_extract_diagonals.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_extract_real_diagonals.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_form_Q.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_l.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_realify.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_realify_subdiagonal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_recover_tau.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_scale_diagonals.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_shift_U.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tridiag_UT_u.c
	)