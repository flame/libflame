##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_compute_scaling.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_eigval_v.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_eigval_v_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_eigval_v_opt_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_find_perfshift.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_find_submatrix.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_francis_v.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_francis_v_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_iteracc_v.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_iteracc_v_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_iteracc_v_opt_var3.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_v.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_v_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_v_opt_var2.c
	)