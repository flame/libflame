##Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.##

add_library(FLA_LAPACK_AVX512
     OBJECT
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_avx512_kernels.h
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_dgetrf_small_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_sgetrf_small_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_cgetrf_small_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_zgetrf_small_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_dscal_ix1_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_lange_kernels_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_drot_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_zrot_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_dhrot3_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_sscal_ix1_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_zscal_ix1_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_dlarf_kernels_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_dnrm2_blas_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_zlarf_kernels_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_dpotrf_small_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_dspr_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_dcopy_scal_avx512.c
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_dgtsv_kernel_avx512.c
)

if(WIN32)
target_compile_options(FLA_LAPACK_AVX512 PRIVATE /arch:AVX512)
else(UNIX)
target_compile_options(FLA_LAPACK_AVX512 PRIVATE -mavx512f -mavx512dq -mfma)
endif()
