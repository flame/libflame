##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_eigval_n.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_eigval_n_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_francis_n.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_francis_n_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_iteracc_n.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_iteracc_n_opt_var1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_n.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Tevd_n_opt_var1.c
	)