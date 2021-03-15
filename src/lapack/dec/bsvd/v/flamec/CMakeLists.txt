##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_compute_shift.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_compute_tol_thresh.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_find_converged.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_find_max_min.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_find_submatrix.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_francis_v.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_francis_v_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_iteracc_v.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_iteracc_v_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_sinval_v.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_sinval_v_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_v.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_v_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Bsvd_v_opt_var2.c
	)