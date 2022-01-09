##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv.h
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_internal.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_solve.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_LU_piv.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_LU_piv_solve.c
${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_piv_small_d.c
${CMAKE_CURRENT_SOURCE_DIR}/lapack_sgetrf.c
${CMAKE_CURRENT_SOURCE_DIR}/lapack_sgetf2.c
${CMAKE_CURRENT_SOURCE_DIR}/lapack_cgetrf.c
${CMAKE_CURRENT_SOURCE_DIR}/lapack_cgetf2.c
${CMAKE_CURRENT_SOURCE_DIR}/lapack_zgetrf.c
${CMAKE_CURRENT_SOURCE_DIR}/lapack_zgetf2.c


	)
