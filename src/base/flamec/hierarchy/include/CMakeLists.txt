##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLASH.h
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_blas1_prototypes.h
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_blas2_prototypes.h
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_blas3_prototypes.h
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_lapack_prototypes.h
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_macro_defs.h
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_main_prototypes.h
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_util_prototypes.h
	)