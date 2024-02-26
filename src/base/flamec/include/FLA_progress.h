/*
    Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
*/
/* typedef long integer integer; */
#ifdef __cplusplus
  // For C++, include stdint.h.
#include <stdint.h> // skipped
#elif __STDC_VERSION__ >= 199901L
  // For C99 (or later), include stdint.h.
#include <stdint.h> // skipped
#else
  // When stdint.h is not available, manually typedef the types we will use.
#ifdef _WIN32
typedef          __int32  int32_t;
typedef unsigned __int32 uint32_t;
typedef          __int64  int64_t;
typedef unsigned __int64 uint64_t;
#else
#error "Attempting to compile on pre-C99 system without stdint.h."
#endif
#endif
/* typedef long integer integer; */
#ifdef FLA_ENABLE_ILP64
typedef int64_t integer;
typedef uint64_t uinteger;
#else
typedef int integer;
typedef unsigned long int uinteger;
#endif

#define AOCL_FLA_PROGRESS_H 1
typedef int (*aocl_fla_progress_callback)(
const char* const api,
const integer lenapi,
const integer* const progress,
const integer* const current_thread,
const integer* const total_threads
);

void aocl_fla_set_progress(aocl_fla_progress_callback func);
extern volatile aocl_fla_progress_callback aocl_fla_progress_glb_ptr;
#ifndef FLA_ENABLE_WINDOWS_BUILD  
__attribute__((weak))
int aocl_fla_progress(
const char* const api,
const integer lenapi,
const integer* const progress,
const integer* const current_thread,
const integer* const total_threads
);
#endif
// Macro to send update using api name
#define AOCL_FLA_PROGRESS_FUNC_PTR(api,lenapi,progress,tid,nt) \
         if((*aocl_fla_progress_ptr) (api,lenapi, progress, tid, nt)){\
            printf("stop computation \n");\
			exit(0);\
         }\

#if FLA_OPENMP_MULTITHREADING

#define AOCL_FLA_PROGRESS_VAR \
        aocl_fla_progress_callback aocl_fla_progress_ptr = aocl_fla_progress_glb_ptr;\
        static TLS_CLASS_SPEC integer progress_step_count = 0;\
        static TLS_CLASS_SPEC integer progress_thread_id = 0;\
        static TLS_CLASS_SPEC integer progress_total_threads = 1;\
        progress_thread_id = omp_get_thread_num();\
        progress_total_threads = omp_get_num_threads();\

#else

#define AOCL_FLA_PROGRESS_VAR \
        aocl_fla_progress_callback aocl_fla_progress_ptr = aocl_fla_progress_glb_ptr;\
        static TLS_CLASS_SPEC integer progress_step_count = 0;\
        static TLS_CLASS_SPEC integer progress_thread_id = 0;\
        static TLS_CLASS_SPEC integer progress_total_threads = 1;\
        progress_thread_id = 0;\
        progress_total_threads = 1;\

#endif
