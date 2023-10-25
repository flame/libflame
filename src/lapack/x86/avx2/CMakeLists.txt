##Copyright (C) 2023, Advanced Micro Devices, Inc.## 
add_library(AVX2_LU
     OBJECT
${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_avx2_kernels.h
${CMAKE_CURRENT_SOURCE_DIR}/fla_dgeqrf_small_avx2.h
${CMAKE_CURRENT_SOURCE_DIR}/fla_dgeqrf_small_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_dgetrf_small_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_dhrot3_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_drot_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_dscal_ix1_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_sger_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_sscal_ix1_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_zgetrf_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_zrot_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_zscal_ix1_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_lu_piv_small_d_update_tr_matrix_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_dgesvd_nn_small10_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_dgesvd_small6T_avx2.c
${CMAKE_CURRENT_SOURCE_DIR}/fla_dgesvd_nn_small1T_avx2.c
)

if(WIN32)
target_compile_options(AVX2_LU PRIVATE /arch:AVX2)
else(UNIX)
target_compile_options(AVX2_LU PRIVATE -mavx2 -mfma)
endif()
