/* ************************************************************************
 * Copyright (c) 2022-2023 Advanced Micro Devices, Inc. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * ************************************************************************ */

#pragma once
#ifndef FLA_CONTEXT_H
#define FLA_CONTEXT_H

// -- Type and macro definitions -----------------------------------------------

#if defined(FLA_NO_CONTEXT)

// This branch defines a pthread-like API, fla_pthread_*(), and implements it
// in terms of "dummy" code that doesn't depend on POSIX threads or any other
// threading mechanism.
// NOTE: THIS CODE DOES NOT IMPLEMENT THREADING AND IS NOT THREAD-SAFE!

// -- pthread types --

typedef int fla_pthread_mutex_t;
typedef int fla_pthread_once_t;

// -- pthreads macros --

#define FLA_PTHREAD_MUTEX_INITIALIZER 0
#define FLA_PTHREAD_ONCE_INIT 0

#elif defined(_MSC_VER) // !defined(FLA_NO_CONTEXT)

#include <windows.h>
// This branch defines a pthread-like API, fla_pthread_*(), and implements it
// in terms of Windows API calls.

// -- pthread types --
typedef SRWLOCK   fla_pthread_mutex_t;
typedef INIT_ONCE fla_pthread_once_t;

// -- pthreads macros --

#define FLA_PTHREAD_MUTEX_INITIALIZER SRWLOCK_INIT
#define FLA_PTHREAD_ONCE_INIT INIT_ONCE_STATIC_INIT

#else // !defined(FLA_NO_CONTEXT) && !defined(_MSC_VER)

#include <pthread.h>

// This branch defines a pthreads-like API, fla_pthreads_*(), and implements it
// in terms of the corresponding pthreads_*() types, macros, and function calls.

// -- pthread types --

typedef pthread_mutex_t fla_pthread_mutex_t;
typedef pthread_once_t  fla_pthread_once_t;

// -- pthreads macros --

#define FLA_PTHREAD_MUTEX_INITIALIZER PTHREAD_MUTEX_INITIALIZER
#define FLA_PTHREAD_ONCE_INIT PTHREAD_ONCE_INIT

#endif

// -- Function definitions -----------------------------------------------------

// -- pthread_mutex_*() --

int fla_pthread_mutex_lock(fla_pthread_mutex_t *mutex);

int fla_pthread_mutex_unlock(fla_pthread_mutex_t *mutex);

// -- pthread_once() --

void fla_pthread_once(fla_pthread_once_t *once, void (*init)(void));

/******************************************************************************************
 * \brief fla_context is a structure holding the number of threads, ISA information
 * It gets initialised by fla_init_once().
 *****************************************************************************************/
typedef struct _fla_context
{
    // num of threads
    int num_threads;
    FLA_Bool    is_fma;
    FLA_Bool    is_avx2;
    FLA_Bool    is_avx512;
    FLA_Bool    libflame_mt; // num_threads is set using libFLAME environment variable or using OpenMP.
} fla_context;

#define FLA_CONTEXT_INITIALIZER \
    { \
      .num_threads = -1, \
      .is_fma      = FALSE, \
      .is_avx2     = FALSE, \
      .is_avx512   = FALSE, \
      .libflame_mt = FALSE, \
    }

extern fla_context global_context;

typedef struct _fla_tl_context
{
    // num of threads
    int num_threads;
    FLA_Bool    libflame_mt; // num_threads is set using libFLAME environment variable or using OpenMP.
} fla_tl_context;

#define FLA_TL_CONTEXT_INITIALIZER \
    { \
      .num_threads = -1, \
      .libflame_mt = FALSE, \
    }

extern TLS_CLASS_SPEC fla_tl_context tl_context;

/*! \ingroup aux_module
 *  \brief Initialise various framework variables including
 *  1.context
 *  2.number of threads from environment
 *
 *  \retval none.
 */
void aocl_fla_init();

/*! \ingroup aux_module
 *  \brief Deallocate and clean all initalized buffers
 */
void aocl_fla_finalize();

#endif // FLA_CONTEXT_H
