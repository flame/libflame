##Copyright (C) 2020, Advanced Micro Devices, Inc.## 
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/cspffrt2_fla.c
${CMAKE_CURRENT_SOURCE_DIR}/cspffrtx_fla.c
${CMAKE_CURRENT_SOURCE_DIR}/dspffrt2_fla.c
${CMAKE_CURRENT_SOURCE_DIR}/dspffrtx_fla.c
${CMAKE_CURRENT_SOURCE_DIR}/sspffrt2_fla.c
${CMAKE_CURRENT_SOURCE_DIR}/sspffrtx_fla.c
${CMAKE_CURRENT_SOURCE_DIR}/zspffrt2_fla.c
${CMAKE_CURRENT_SOURCE_DIR}/zspffrtx_fla.c
	)