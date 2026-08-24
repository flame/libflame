/* ************************************************************************
 * Copyright (C) 2022-2026, Advanced Micro Devices, Inc. All rights reserved.
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
#include "Capi/au/cpuid/cpuid.h"

#if defined(FLA_NO_CONTEXT)

// This branch defines a pthread-like API, fla_pthread_*(), and implements it
// in terms of "dummy" code that doesn't depend on POSIX threads or any other
// threading mechanism.
// NOTE: THIS CODE DOES NOT IMPLEMENT THREADING AND IS NOT THREAD-SAFE!

int fla_pthread_mutex_lock(fla_pthread_mutex_t *mutex)
{
    // return pthread_mutex_lock( mutex );
    return 0;
}

int fla_pthread_mutex_unlock(fla_pthread_mutex_t *mutex)
{
    // return pthread_mutex_unlock( mutex );
    return 0;
}

// -- pthread_once() --

void fla_pthread_once(fla_pthread_once_t *once, void (*init)(void))
{
    // pthread_once( once, init );
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

static FLA_Bool fla_init_once_wrapper(fla_pthread_once_t *once, void *param, void **context)
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
fla_context fla_global_context = FLA_CONTEXT_INITIALIZER;

// The global fla_context structure, which holds the updated thread-local
// thread count
TLS_CLASS_SPEC fla_tl_context_t fla_tl_context = FLA_TL_CONTEXT_INITIALIZER;
TLS_CLASS_SPEC FLA_Bool fla_tl_context_init = FALSE;

// variable to check if cpu archtecture is explicitly set
TLS_CLASS_SPEC fla_dim_t __attribute__((unused)) fla_req_id = -1;

// Variable to get the ISA architecture to use
TLS_CLASS_SPEC fla_arch_t fla_arch_id = -1;

// Keep track if AOCL_ENABLE_INSTRUCTIONS environment variable was set.
TLS_CLASS_SPEC integer __attribute__((unused)) fla_aocl_e_i = FALSE;

// A mutex to allow synchronous access to global_thread.
TLS_CLASS_SPEC fla_pthread_mutex_t fla_global_thread_mutex = FLA_PTHREAD_MUTEX_INITIALIZER;

/********************************************************************************
 * \brief fla_env_get_var is a function used to query the environment
 * variable and convert the string into fla_dim_t and return the same
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
        // If there was no error, convert the string to an fla_dim_t and
        // prepare to return that fla_dim_t.
        r_val = (int)strtol(str, NULL, 10);
    }
    else
    {
        // If there was an error, use the "fallback" as the return value.
        r_val = fallback;
    }

    return r_val;
}

/* Get value of AOCL_ENABLE_INSTRUCTIONS environment variable if set*/
int fla_env_get_var_arch_type(const char *env, int fallback)
{
    int r_val, size, i;
    char *str;
    str = getenv(env);

    if(str != NULL)
    {
        // If there was no error, convert the string to an fla_dim_t
        r_val = (int)strtol(str, NULL, 10);

        if(r_val == 0)
        {
            // Look for non-numeric values

            // convert string to lowercase
            size = strlen(str);
            for(i = 0; i <= size; i++)
            {
                str[i] = tolower(str[i]);
            }

            if(strcmp(str, "avx512") == 0)
            {
                r_val = FLA_ARCH_AVX512;
            }
            else if(strcmp(str, "avx2") == 0)
            {
                r_val = FLA_ARCH_AVX2;
            }
            else if(strcmp(str, "avx") == 0)
            {
                r_val = FLA_ARCH_AVX;
            }
            else if(strcmp(str, "sse2") == 0)
            {
                r_val = FLA_ARCH_SSE2;
            }
            else if(strcmp(str, "generic") == 0)
            {
                r_val = FLA_ARCH_GENERIC;
            }
            else
            {
                // Invalid value was set for env variable
                // Print error message and abort
                fprintf(stderr,
                    "Invalid architecture id value set for env var AOCL_ENABLE_INSTRUCTIONS.\n");
                FLA_Abort();
            }
        }
        else
        {
            // Invalid value was set for env variable
            // Print error message and abort
            fprintf(stderr,
                    "Invalid architecture id value set for env var AOCL_ENABLE_INSTRUCTIONS.\n");
            FLA_Abort();
        }
    }
    else
    {
        r_val = fallback;
    }

    return r_val;
}

/* Return the CPU ISA architecture to use */
fla_arch_t fla_arch_query_id(void)
{
    return fla_arch_id;
}

/* Return boolean that indicates if AOCL_ENABLE_INSTRUCTIONS environment variable has been set */
integer fla_aocl_enable_instruction_query(void)
{
    return fla_aocl_e_i;
}

/* Get value of AOCL_ENABLE_INSTRUCTIONS environment variable
 * Sets boolean fla_aocl_e_i if valid value returned*/
void fla_get_arch_info_from_env(fla_context *context)
{
    fla_req_id = fla_env_get_var_arch_type("AOCL_ENABLE_INSTRUCTIONS", -1);
    if(fla_req_id != -1)
    {
        fla_aocl_e_i = TRUE;
    }
}

/* Update global_context if FLA_NUM_THREADS is set
 * OpenMP threading is set in subsequent call to
 * fla_thread_update_rntm_from_env() */
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

/* Update global context and thread local context
 * with appropriate threads to use*/
void fla_thread_update_rntm_from_env(fla_tl_context_t *context)
{

#ifdef FLA_OPENMP_MULTITHREADING

    if( !fla_tl_context_init )
    {
        // On first call for each thread, need to check settings from
        // BLIS environment variables in global_context. First, set
        // fla_tl_context_init to TRUE for subsequent calls.
        fla_tl_context_init = TRUE;

        // Acquire the mutex protecting global_thread.
        fla_pthread_mutex_lock(&fla_global_thread_mutex);

        // Copy values from global_context.
        context->num_threads = fla_global_context.num_threads;
        context->libflame_mt = fla_global_context.libflame_mt;

        // Release the mutex protecting global_thread.
        fla_pthread_mutex_unlock(&fla_global_thread_mutex);
    }

    // If FLA_NUM_THREADS was not set, read OpenMP's omp_get_max_threads()
    // to get maximum number of threads that library can use. We also
    // need to consider the number of active OpenMP levels and which
    // level we are at.
    if(!context->libflame_mt)
    {
        int active_level = omp_get_active_level();
        int max_levels = omp_get_max_active_levels();
        if(active_level < max_levels)
        {
            context->num_threads = omp_get_max_threads();
        }
        else
        {
            context->num_threads = 1;
        }
    }

#else

    if( !fla_tl_context_init )
    {
        // First, set fla_tl_context_init to TRUE for subsequent calls.
        fla_tl_context_init = TRUE;

        // Always set maximum threads as 1. These should never be
        // changed so only set on first call.
        context->num_threads = 1;
        context->libflame_mt = FALSE;
    }

#endif
}

/* Find the ISA to use for code execution paths
 * based on target CPU architecture and preference
 * set by user using AOCL_ENABLE_INSTRUCTIONS env
 * variable */
void fla_isa_init(fla_context *context)
{
    fla_arch_id = FLA_ARCH_GENERIC;
    const char* const flags_array[] = {"avx2", "avx512f", "avx512dq"};
    au_cpu_num_t cpu_num = AU_CURRENT_CPU_NUM;

    integer cpu_supports_isa = au_cpuid_has_flags(cpu_num, flags_array, 1);

    // Check if the target supports AVX2/AVX512
    if(cpu_supports_isa)
    {
        context->is_avx2 = TRUE;
        fla_arch_id = FLA_ARCH_AVX2;
    }

    /* AVX512DQ is superset of AVX512F
     * if AVX512DQ is supported, then AVX512F is also supported */
    cpu_supports_isa = au_cpuid_has_flags(cpu_num, &flags_array[2], 1);
    if(cpu_supports_isa)
    {
        context->is_avx512 = TRUE;
        fla_arch_id = FLA_ARCH_AVX512DQ;
    }
    else
    {
        cpu_supports_isa = au_cpuid_has_flags(cpu_num, &flags_array[1], 1);
        if(cpu_supports_isa)
        {
            context->is_avx512 = TRUE;
            fla_arch_id = FLA_ARCH_AVX512;
        }
    }

#ifdef FLA_STRICT_ARCH
    /* Compile-time override: force architecture ID, but cap to CPU capability.
       Usage: add -DFLA_STRICT_ARCH=FLA_ARCH_AVX2 (or other FLA_ARCH_* value) in CMake
       This overrides AOCL_ENABLE_INSTRUCTIONS behavior. */
    {
        fla_arch_t strict_arch_id = (fla_arch_t) FLA_STRICT_ARCH;

        fla_arch_id = fla_arch_id < strict_arch_id ? fla_arch_id : strict_arch_id;
    }

#else
    // Check user has set AOCL_ENABLE_INSTRUCTIONS env variable
    fla_get_arch_info_from_env(context);

    if(fla_aocl_e_i)
    {
        // We always select the best suitable ISA. If user has chosen
        // higher level ISA than supported by target CPU, we choose
        // best supported architecture on target CPU.
        // If user has chosen a lower level ISA, then same will
        // be used

        if(fla_req_id == FLA_ARCH_AVX2 && fla_arch_id == FLA_ARCH_AVX512)
        {
            fla_arch_id = FLA_ARCH_AVX2;
        }
        else if(fla_req_id < FLA_ARCH_AVX2)
        {
            // For any ISA path less than AVX2, we set generic(c reference) code path
            fla_arch_id = FLA_ARCH_GENERIC;
        }
    }
#endif

    // set the ARCH ID to use
    context->arch_id = fla_arch_id;
}

// -----------------------------------------------------------------------------
void fla_context_init(void)
{
    // Read the environment variables and use them to initialize the
    // global runtime object.
    fla_thread_init_rntm_from_env(&fla_global_context);

    // Read target ISA path to run
    fla_isa_init(&fla_global_context);
}

// -----------------------------------------------------------------------------

void fla_context_finalize(void) {}

// -----------------------------------------------------------------------------

// A pthread_once_t variable is a pthread structure used in pthread_once().
// pthread_once() is guaranteed to execute exactly once among all threads that
// pass in this control object. Thus, we need one for initialization and a
// separate one for finalization.
static fla_pthread_once_t once_init = FLA_PTHREAD_ONCE_INIT;
static fla_pthread_once_t once_finalize = FLA_PTHREAD_ONCE_INIT;
static FLA_Bool fla_context_initialized = FALSE;

void aocl_fla_init(void)
{    
    // Fast path: if already initialized, return immediately
    if (fla_context_initialized == TRUE)
        return;
    
    // Thread-safe initialization using pthread_once
    fla_pthread_once(&once_init, fla_context_init);

    // Mark as initialized
    fla_context_initialized = TRUE;
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
    fla_thread_update_rntm_from_env(&fla_tl_context);

    return fla_tl_context.num_threads;
}

void fla_thread_set_num_threads(int n_threads)
{

#ifdef FLA_OPENMP_MULTITHREADING

    // We must ensure that global_thread has been initialized.
    aocl_fla_init();

    // Update values in tl_context for future reference
    fla_tl_context.num_threads = n_threads;
    fla_tl_context.libflame_mt = TRUE;

#endif
}
