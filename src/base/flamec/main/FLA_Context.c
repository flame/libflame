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

#include "FLAME.h"
#include "alci/arch.h"

#if defined(FLA_NO_CONTEXT)

// This branch defines a pthread-like API, fla_pthread_*(), and implements it
// in terms of "dummy" code that doesn't depend on POSIX threads or any other
// threading mechanism.
// NOTE: THIS CODE DOES NOT IMPLEMENT THREADING AND IS NOT THREAD-SAFE!

int fla_pthread_mutex_lock(fla_pthread_mutex_t *mutex)
{
    //return pthread_mutex_lock( mutex );
    return 0;
}

int fla_pthread_mutex_unlock(fla_pthread_mutex_t *mutex)
{
    //return pthread_mutex_unlock( mutex );
    return 0;
}

// -- pthread_once() --

void fla_pthread_once(fla_pthread_once_t *once, void (*init)(void))
{
    //pthread_once( once, init );
    return;
}

#elif defined(_MSC_VER) // !defined(FLA_DISABLE_SYSTEM)

#include <errno.h>

// This branch defines a pthread-like API, fla_pthread_*(), and implements it
// in terms of Windows API calls.

// -- pthread_mutex_*() --

int fla_pthread_mutex_lock(fla_pthread_mutex_t *mutex)
{
    AcquireSRWLockExclusive(mutex);
    return 0;
}

int fla_pthread_mutex_unlock(fla_pthread_mutex_t *mutex)
{
    ReleaseSRWLockExclusive(mutex);
    return 0;
}

// -- pthread_once() --

static FLA_Bool
    fla_init_once_wrapper(fla_pthread_once_t *once, void *param, void **context)
{
    (void)once;
    (void)context;
    typedef void (*callback)(void);
    ((callback)param)();
    return TRUE;
}

void fla_pthread_once(fla_pthread_once_t *once, void (*init)(void))
{
    InitOnceExecuteOnce(once, fla_init_once_wrapper, init, NULL);
}

#else // !defined(FLA_NO_CONTEXT) && !defined(_MSC_VER)

// This branch defines a pthreads-like API, fla_pthreads_*(), and implements it
// in terms of the corresponding pthreads_*() types, macros, and function calls.
// This branch is compiled for Linux and other non-Windows environments where
// we assume that *some* implementation of pthreads is provided (although it
// may lack barriers--see below).

// -- pthread_mutex_*() --

int fla_pthread_mutex_lock(fla_pthread_mutex_t *mutex)
{
    return pthread_mutex_lock(mutex);
}

int fla_pthread_mutex_unlock(fla_pthread_mutex_t *mutex)
{
    return pthread_mutex_unlock(mutex);
}

// -- pthread_once() --

void fla_pthread_once(fla_pthread_once_t *once, void (*init)(void))
{
    pthread_once(once, init);
}

#endif // !defined(FLA_NO_CONTEXT) && !defined(_MSC_VER)

// The global fla_context structure, which holds the global thread count
// and ISA settings
fla_context global_context = FLA_CONTEXT_INITIALIZER;

// The global fla_context structure, which holds the updated thread-local
// thread count
TLS_CLASS_SPEC fla_tl_context tl_context = FLA_TL_CONTEXT_INITIALIZER;
TLS_CLASS_SPEC FLA_Bool tl_context_init = FALSE;

// A mutex to allow synchronous access to global_thread.
fla_pthread_mutex_t global_thread_mutex = FLA_PTHREAD_MUTEX_INITIALIZER;

/********************************************************************************
 * \brief fla_env_get_var is a function used to query the environment
 * variable and convert the string into integer and return the same
 ********************************************************************************/
int fla_env_get_var(const char *env, int fallback)
{
    int r_val;
    char *str;

    // Query the environment variable and store the result in str.
    str = getenv(env);

    // Set the return value based on the string obtained from getenv().
    if(str != NULL)
    {
        // If there was no error, convert the string to an integer and
        // prepare to return that integer.
        r_val = (int)strtol(str, NULL, 10);
    }
    else
    {
        // If there was an error, use the "fallback" as the return value.
        r_val = fallback;
    }

    return r_val;
}

// This updates global_context
void fla_thread_init_rntm_from_env(fla_context *context)
{
    int nt;
    FLA_Bool libflame_mt;

#ifdef FLA_OPENMP_MULTITHREADING
    // Try to read FLA_NUM_THREADS first.
    nt = fla_env_get_var("FLA_NUM_THREADS", -1);

    // If FLA_NUM_THREADS was not set, set OpenMP threading in a
    // subsequent call to fla_thread_update_rntm_from_env().
    if(nt == -1)
    {
        libflame_mt = FALSE;
    }
    else
    {
        libflame_mt = TRUE;
    }
#else
    // If multi-thread mode not configured, set maximum threads as 1
    nt = 1;
    libflame_mt = FALSE;
#endif

    context->num_threads = nt;
    context->libflame_mt = libflame_mt;
}

// This updates tl_context
void fla_thread_update_rntm_from_env(fla_tl_context *context)
{

#ifdef FLA_OPENMP_MULTITHREADING

    if( !tl_context_init )
    {
        // On first call for each thread, need to check settings from
        // BLIS environment variables in global_context. First, set
        // tl_context_init to TRUE for subsequent calls.
        tl_context_init = TRUE;

        // Acquire the mutex protecting global_thread.
        fla_pthread_mutex_lock(&global_thread_mutex);

        // Copy values from global_context.
        context->num_threads = global_context.num_threads;
        context->libflame_mt = global_context.libflame_mt;

        // Release the mutex protecting global_thread.
        fla_pthread_mutex_unlock(&global_thread_mutex);
    }

    // If FLA_NUM_THREADS was not set, read OpenMP's omp_get_max_threads()
    // to get maximum number of threads that library can use. We also
    // need to consider the number of active OpenMP levels and which
    // level we are at.
    if( !context->libflame_mt )
    {
        int active_level = omp_get_active_level();
        int max_levels = omp_get_max_active_levels();
        if ( active_level < max_levels )
        {
            context->num_threads = omp_get_max_threads();
        }
        else
        {
            context->num_threads = 1;
        }
    }

#else

    if( !tl_context_init )
    {
        // First, set tl_context_init to TRUE for subsequent calls.
        tl_context_init = TRUE;

        // Always set maximum threads as 1. These should never be
        // changed so only set on first call.
        context->num_threads = 1;
        context->libflame_mt = FALSE;
    }

#endif

}

void fla_isa_init(fla_context *context)
{
    // Check if the target supports AVX2
    if(alcpu_flag_is_available(ALC_E_FLAG_AVX2))
    {
        context->is_avx2 = TRUE;
    }
    if (alcpu_flag_is_available(ALC_E_FLAG_AVX512F))
    {
        context->is_avx512 = TRUE;
    }
}

// -----------------------------------------------------------------------------
void fla_context_init(void)
{
    // Read the environment variables and use them to initialize the
    // global runtime object.
    fla_thread_init_rntm_from_env(&global_context);

    // Read target ISA if it supports avx512
    fla_isa_init(&global_context);
}

// -----------------------------------------------------------------------------

void fla_context_finalize(void) {}

// -----------------------------------------------------------------------------

// A pthread_once_t variable is a pthread structure used in pthread_once().
// pthread_once() is guaranteed to execute exactly once among all threads that
// pass in this control object. Thus, we need one for initialization and a
// separate one for finalization.
static fla_pthread_once_t once_init     = FLA_PTHREAD_ONCE_INIT;
static fla_pthread_once_t once_finalize = FLA_PTHREAD_ONCE_INIT;

void aocl_fla_init(void)
{
    fla_pthread_once(&once_init, fla_context_init);
}

void aocl_fla_finalize(void)
{
    fla_pthread_once(&once_finalize, fla_context_finalize);
}

int fla_thread_get_num_threads(void)
{
    // We must ensure that global_context and tl_context have been initialized.
    aocl_fla_init();

    // Update the OpenMP information from the runtime, unless FLA_NUM_THREADS
    // was set or fla_thread_set_num_threads() was called.
    fla_thread_update_rntm_from_env(&tl_context);

    return tl_context.num_threads;
}

void fla_thread_set_num_threads(int n_threads)
{

#ifdef FLA_OPENMP_MULTITHREADING

    // We must ensure that global_thread has been initialized.
    aocl_fla_init();

    // Update values in tl_context for future reference
    tl_context.num_threads = n_threads;
    tl_context.libflame_mt = TRUE;

#endif

}
