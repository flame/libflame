/* ************************************************************************
 * Copyright (c) 2022 Advanced Micro Devices, Inc.
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

// The global fla_context structure, which holds the global thread,ISA settings
fla_context global_context;

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

void fla_thread_init_rntm_from_env(fla_context *context)
{
    int nt;

#ifdef FLA_OPENMP_MULTITHREADING
    // Try to read FLA_NUM_THREADS first.
    nt = fla_env_get_var("FLA_NUM_THREADS", -1);

    // If FLA_NUM_THREADS was not set, read OpenMP's omp_get_max_threads() to get maximum number
    // of threads that library can use
    if(nt == -1)
        nt = omp_get_max_threads();
#else
    // If multi-thread mode not configured, set maximum threads as 1
    nt = 1;
#endif

    context->num_threads = nt;
}

void fla_isa_init(fla_context *context)
{
    // Check if the target supports AVX512
    if(alc_cpu_has_avx512f())
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
    // We must ensure that global_rntm has been initialized.
    aocl_fla_init();

    return global_context.num_threads;
}

void fla_thread_set_num_threads(int n_threads)
{
    // We must ensure that global_thread has been initialized.
    aocl_fla_init();

    // Acquire the mutex protecting global_thread.
    fla_pthread_mutex_lock(&global_thread_mutex);

    global_context.num_threads = n_threads;

    // Release the mutex protecting global_thread.
    fla_pthread_mutex_unlock(&global_thread_mutex);
}
