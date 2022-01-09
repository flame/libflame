##Copyright (C) 2020, Advanced Micro Devices, Inc.##
target_sources("${PROJECT_NAME}" 
     PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/blis_f77_name_mangling.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_macro_defs.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_prototypes_blas.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_prototypes_cblas.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_prototypes_fused1.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_prototypes_level1.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_prototypes_level2.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_prototypes_level3.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_prototypes_misc.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_prototypes_query.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_prototypes_util.h
${CMAKE_CURRENT_SOURCE_DIR}/blis_type_defs.h
${CMAKE_CURRENT_SOURCE_DIR}/blis1.h
	)