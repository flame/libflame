##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/fla_dlamch.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_lsame.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_slamch.c
	)