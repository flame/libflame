##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Hermitianize.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_LU_find_zero_on_diagonal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Max_elemwise_diff.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Norm1.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Obj_create_diag_panel.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Random_matrix.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Random_spd_matrix.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Set.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Shift_diag.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Triangularize.c
	)