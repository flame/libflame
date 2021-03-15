##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Accum_T_UT_check.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_H2_UT_check.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Apply_HUD_UT_check.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Mach_params_check.c
	)