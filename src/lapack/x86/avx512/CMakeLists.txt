##Copyright (C) 2023, Advanced Micro Devices, Inc.## 
add_library(AVX512_LU
     OBJECT
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_lapack_avx512_kernels.h
     ${CMAKE_CURRENT_SOURCE_DIR}/fla_dgetrf_small_avx512.c
)

if(WIN32)
target_compile_options(AVX512_LU PRIVATE /arch:AVX512)
else(UNIX)
target_compile_options(AVX512_LU PRIVATE -mavx512f -mfma)
endif()
