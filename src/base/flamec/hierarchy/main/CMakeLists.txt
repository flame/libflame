##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Axpy_other.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Copy_other.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_Obj.c
${CMAKE_CURRENT_SOURCE_DIR}/FLASH_View.c
	)