##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_il.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_iu.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nl.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_Eig_gest_nu.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Eig_gest.c
	)