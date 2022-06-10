###############################################################################
# Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
###############################################################################

target_sources(test_libFLAME_main
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/test_syevd.c
${CMAKE_CURRENT_SOURCE_DIR}/test_orgqr.c
${CMAKE_CURRENT_SOURCE_DIR}/test_potrs.c
${CMAKE_CURRENT_SOURCE_DIR}/test_getrf.c
${CMAKE_CURRENT_SOURCE_DIR}/test_getri.c
${CMAKE_CURRENT_SOURCE_DIR}/test_getrs.c
${CMAKE_CURRENT_SOURCE_DIR}/test_geqrf.c
${CMAKE_CURRENT_SOURCE_DIR}/test_gerq2.c
${CMAKE_CURRENT_SOURCE_DIR}/test_gerqf.c
${CMAKE_CURRENT_SOURCE_DIR}/test_geevx.c
${CMAKE_CURRENT_SOURCE_DIR}/test_gesdd.c
${CMAKE_CURRENT_SOURCE_DIR}/test_potrf.c
${CMAKE_CURRENT_SOURCE_DIR}/test_lapack.c
${CMAKE_CURRENT_SOURCE_DIR}/test_lapack.h
${CMAKE_CURRENT_SOURCE_DIR}/test_routines.h
)
