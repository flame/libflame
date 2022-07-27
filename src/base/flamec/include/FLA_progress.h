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
char* api,
integer lenapi,
integer *progress,
integer *current_thread,
integer *total_threads
);

void aocl_fla_set_progress(aocl_fla_progress_callback func);
extern aocl_fla_progress_callback aocl_fla_progress_ptr;
#ifndef FLA_ENABLE_WINDOWS_BUILD  
__attribute__((weak))
int aocl_fla_progress(
char* api,
integer lenapi,
integer *progress,
integer *current_thread,
integer *total_threads
);
#endif
// Macro to send update using api name
#define AOCL_FLA_PROGRESS_FUNC_PTR(api,lenapi,progress,tid,nt) \
         if((*aocl_fla_progress_ptr) (api,lenapi, progress, tid, nt)){\
            printf("stop computation \n");\
			exit(0);\
         }\

#define AOCL_FLA_PROGRESS_VAR \
        static TLS_CLASS_SPEC integer step_count=0;\
        static TLS_CLASS_SPEC integer size=0;\
        static TLS_CLASS_SPEC integer thread_id = 0;\
        static TLS_CLASS_SPEC integer total_threads = 1;\
        if(aocl_fla_progress_ptr)\
        {\
        /* Current implementation returns threadid as 0 and total_threads as 1*/ \
        /* even if invoked from multithreaded application. */ \
        /* Support for actual thread number will be added in future */ \
            thread_id = 0;\
            total_threads =  1;\
        }\
