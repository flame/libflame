##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/FLA_LU_incpiv_aux.h
	)
	
add_subdirectory(fla)
add_subdirectory(flash)