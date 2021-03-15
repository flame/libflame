##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/bl1_axmyv2.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_axpyv2b.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_axpyv2bdotaxpy.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_axpyv3b.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_dotaxmyv2.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_dotaxpy.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_dotsv2.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_dotsv3.c
${CMAKE_CURRENT_SOURCE_DIR}/bl1_dotv2axpyv2b.c
	)