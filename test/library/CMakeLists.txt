###############################################################################
# Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################

target_sources(test_libFLAME_main
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/test_common.c
${CMAKE_CURRENT_SOURCE_DIR}/test_common.h
${CMAKE_CURRENT_SOURCE_DIR}/test_linear_solvers.c
${CMAKE_CURRENT_SOURCE_DIR}/test_linear_solvers.h
${CMAKE_CURRENT_SOURCE_DIR}/test_prototype.h
)
