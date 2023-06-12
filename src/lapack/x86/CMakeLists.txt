##Copyright (C) 2023, Advanced Micro Devices, Inc.## 
add_subdirectory(front)
remove_definitions(/arch:AVX)
add_subdirectory(avx2)
remove_definitions(/arch:AVX2)
add_subdirectory(avx512)
