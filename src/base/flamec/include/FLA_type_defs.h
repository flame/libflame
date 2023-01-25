/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef FLA_TYPE_DEFS_H
#define FLA_TYPE_DEFS_H

#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
#ifdef FLA_ENABLE_TIDSP
#include <ti/omp/omp.h>
#else
#include <omp.h>
#endif
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
#include <pthread.h>
#endif


// --- Complex type definitions -----------------------------------------------

#ifndef _DEFINED_SCOMPLEX
#define _DEFINED_SCOMPLEX
typedef struct scomplex
{
  float real, imag;
} scomplex;
#endif

#ifndef _DEFINED_DCOMPLEX
#define _DEFINED_DCOMPLEX
typedef struct dcomplex
{
  double real, imag;
} dcomplex;
#endif


// --- Parameter and return type definitions ----------------------------------

typedef int FLA_Bool;
typedef int FLA_Error;
typedef int FLA_Quadrant;
typedef int FLA_Datatype;
typedef int FLA_Elemtype;
typedef int FLA_Side;
typedef int FLA_Uplo;
typedef int FLA_Trans;
typedef int FLA_Conj;
typedef int FLA_Diag;
typedef int FLA_Dimension;
typedef int FLA_Pivot_type;
typedef int FLA_Direct;
typedef int FLA_Store;
typedef int FLA_Matrix_type;
typedef int FLA_Precision;
typedef int FLA_Domain;
typedef int FLA_Inv;
typedef int FLA_Evd_type;
typedef int FLA_Svd_type;
typedef int FLA_Machval;
typedef int FLA_Diag_off;

#ifndef _DEFINED_DIM_T
#define _DEFINED_DIM_T
typedef unsigned long dim_t;
#endif

// --- Intrinsic/assembly definitions ----------------------------------------

#if FLA_VECTOR_INTRINSIC_TYPE == FLA_SSE_INTRINSICS

#include "pmmintrin.h"

//typedef double v2df __attribute__ ((vector_size (16)));

typedef union
{
    __m128  v; 
    float   f[4];
} v4sf_t;

typedef union
{
    __m128d v; 
    double  d[2];
} v2df_t;

#endif

// --- FLAME object definitions -----------------------------------------------

typedef struct FLA_Lock_s     FLA_Lock;
typedef struct FLA_RWLock_s   FLA_RWLock;

//#ifdef FLA_ENABLE_MULTITHREADING
struct FLA_Lock_s
{
  // Implementation-specific lock object
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_lock_t       lock;
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_mutex_t  lock;
#endif
};
struct FLA_RWLock_s
{
  // Implementation-specific lock object
#if   FLA_MULTITHREADING_MODEL == FLA_OPENMP
  omp_lock_t       lock;
#elif FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  pthread_rwlock_t lock;
#endif
};
//#endif

#ifdef FLA_ENABLE_SUPERMATRIX
typedef int                   FLASH_Verbose;
typedef int                   FLASH_Data_aff;

typedef struct FLASH_Queue_s  FLASH_Queue;
typedef struct FLASH_Task_s   FLASH_Task;
typedef struct FLASH_Dep_s    FLASH_Dep;
#endif
typedef struct FLASH_Thread_s FLASH_Thread;

typedef struct FLA_Obj_struct
{
  // Basic object description fields
  FLA_Datatype  datatype;
  FLA_Elemtype  elemtype;
  dim_t         m;
  dim_t         n;
  dim_t         rs;
  dim_t         cs;
  dim_t         m_inner;
  dim_t         n_inner;
  unsigned long id;
  dim_t         m_index;
  dim_t         n_index;

  dim_t         n_elem_alloc;
  void*         buffer;
  int           buffer_info;

  FLA_Uplo      uplo;

#ifdef FLA_ENABLE_SUPERMATRIX
  // Fields for supermatrix
  int           n_read_blocks;
  int           n_write_blocks;

  // All the tasks that previously read this block, anti-dependency
  int           n_read_tasks;
  FLASH_Dep*    read_task_head;
  FLASH_Dep*    read_task_tail;

  // Task that last overwrote this block, flow dependency
  FLASH_Task*   write_task;
#endif
} FLA_Base_obj;

typedef struct FLA_Obj_view
{
  // Basic object view description fields
  dim_t         offm;
  dim_t         offn;
  dim_t         m;
  dim_t         n;
  dim_t         m_inner;
  dim_t         n_inner;

  FLA_Base_obj* base;

} FLA_Obj;

#ifdef FLA_ENABLE_SUPERMATRIX
struct FLASH_Queue_s
{
  // Number of tasks currently in queue
  unsigned int  n_tasks;

  // Pointers to head (front) and tail (back) of queue
  FLASH_Task*   head;
  FLASH_Task*   tail;
};

struct FLASH_Task_s
{
  // Execution information
  int           n_ready;

  // Labels
  int           order;
  int           queue;
  int           height;
  int           thread;
  int           cache;
  FLA_Bool      hit;
      
  // Function pointer
  void*         func;

  // Control tree pointer
  void*         cntl;

  // Name of task
  char*         name;

  // GPU enabled task
  FLA_Bool      enabled_gpu;

  // HIP enabled task
  FLA_Bool      enabled_hip;

  // Integer arguments
  int           n_int_args;
  int*          int_arg;

  // Constant FLA_Obj arguments
  int           n_fla_args;
  FLA_Obj*      fla_arg;

  // Input FLA_Obj arguments
  int           n_input_args;
  FLA_Obj*      input_arg;

  // Output FLA_Obj argument
  int           n_output_args;
  FLA_Obj*      output_arg;

  // Number of blocks within all macroblocks
  int           n_macro_args;

  // Number of write after read dependencies
  int           n_war_args;

  // Dependence information
  int           n_dep_args;
  FLASH_Dep*    dep_arg_head;
  FLASH_Dep*    dep_arg_tail;
  
  // Support for a doubly linked list of tasks
  FLASH_Task*   prev_task;
  FLASH_Task*   next_task;

  // Support for a doubly linked list for wait queue
  FLASH_Task*   prev_wait;
  FLASH_Task*   next_wait;
};

struct FLASH_Dep_s
{
  // Task yielding dependency
  FLASH_Task*   task;

  // Support for linked list of FLASH_Deps
  FLASH_Dep*    next_dep;
};
#endif // FLA_ENABLE_SUPERMATRIX

struct FLASH_Thread_s
{
  // The thread's unique identifier
  int       id;

  // Pointer to variables needed to execute SuperMatrix mechanism
  void*     args;

#if FLA_MULTITHREADING_MODEL == FLA_PTHREADS
  // The thread object. Only needed for the POSIX threads implementation.
  pthread_t pthread_obj;
#endif
};

#endif // FLA_TYPE_DEFS_H
