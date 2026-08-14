/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
    Modifications Copyright (c) 2021-2026 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_prototypes.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_util_defs.h"
#include "fla_lapack_avx2_kernels.h"
#include "fla_lapack_avx512_kernels.h"
#include "fla_lapack_lu_small_kernels_d.h"
#include "fla_lapack_lu_small_kernels_s.h"
#include "fla_lapack_lu_small_kernels_c.h"
#include "fla_lapack_x86_common.h"

/*
  GETRF computes an LU factorization of a general M-by-N matrix A
  using partial pivoting with row interchanges.

  INFO
  = 0: successful exit
  < 0: if INFO = -i, the i-th argument had an illegal value - LAPACK_getrf_op_check
  > 0: if INFO = i, U(i,i) is exactly zero. The factorization - FLA_LU_piv
  has been completed, but the factor U is exactly
  singular, and division by zero will occur if it is used
  to solve a system of equations.
*/

extern void DTL_Trace(uint8 ui8LogLevel, uint8 ui8LogType, const int8 *pi8FileName,
                      const int8 *pi8FunctionName, uint32 ui32LineNumber, const int8 *pi8Message);
int fla_thread_get_num_threads(void);
/** Generated wrapper function */
void sgetrf_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgetrf(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgetrf(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgetrf_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgetrf(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgetrf(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void cgetrf_(aocl_int_t *m, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgetrf(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgetrf(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zgetrf_(aocl_int_t *m, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zgetrf(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zgetrf(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void sgetf2_(aocl_int_t *m, aocl_int_t *n, real *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_sgetf2(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_sgetf2(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void dgetf2_(aocl_int_t *m, aocl_int_t *n, doublereal *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_dgetf2(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_dgetf2(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void cgetf2_(aocl_int_t *m, aocl_int_t *n, scomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_cgetf2(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_cgetf2(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

/** Generated wrapper function */
void zgetf2_(aocl_int_t *m, aocl_int_t *n, dcomplex *buff_A, aocl_int_t *ldim_A, aocl_int_t *buff_p, aocl_int_t *info)
{
#if FLA_ENABLE_ILP64
    aocl_lapack_zgetf2(m, n, buff_A, ldim_A, buff_p, info);
#else
    aocl_int64_t m_64 = *m;
    aocl_int64_t n_64 = *n;
    aocl_int64_t ldim_A_64 = *ldim_A;
    aocl_int64_t info_64 = *info;

    aocl_lapack_zgetf2(&m_64, &n_64, buff_A, &ldim_A_64, buff_p, &info_64);

    *info = (aocl_int_t)info_64;
#endif
}

#define LAPACK_getrf(prefix)                                                                 \
    void aocl_lapack_##prefix##getrf(aocl_int64_t *m, aocl_int64_t *n, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, \
                             aocl_int64_t * ldim_A, aocl_int_t * buff_p, aocl_int64_t * info)

#ifndef FLA_ENABLE_SUPERMATRIX

#if FLA_ENABLE_AMD_OPT /* FLA_ENABLE_AMD_OPT */
/* FLA_ENABLE_AMD_OPT enables the code which selects algorithm variants based on size */
#define LAPACK_getrf_body_d(prefix)                                                         \
    extern fla_context fla_global_context;                                                  \
    aocl_int64_t i = 0;                                                                     \
    if(*m == 2 && *n == 2)                                                                  \
    {                                                                                       \
        FLA_LU_PIV_SMALL_D_2x2(i, *n, buff_A, ldim_A, buff_p, *info);                       \
    }                                                                                       \
    else if(*m == 3 && *n == 3)                                                             \
    {                                                                                       \
        FLA_LU_PIV_SMALL_D_3x3(i, *n, buff_A, ldim_A, buff_p, *info);                       \
    }                                                                                       \
    else if(*m == 4 && *n == 4)                                                             \
    {                                                                                       \
        FLA_LU_PIV_SMALL_D_4x4(i, *n, buff_A, ldim_A, buff_p, *info);                       \
    }                                                                                       \
    else if(*m <= FLA_DGETRF_SMALL_THRESH0 && *n <= FLA_DGETRF_SMALL_THRESH0)               \
    {                                                                                       \
        FLA_LU_piv_small_d_var0(m, n, buff_A, ldim_A, buff_p, info);                        \
    }                                                                                       \
    else                                                                                    \
    {                                                                                       \
        /* Initialize global context data */                                                \
        aocl_fla_init();                                                                    \
        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2) && *m < FLA_DGETRF_SMALL_AVX2_THRESH0          \
           && *n < FLA_DGETRF_SMALL_AVX2_THRESH0)                                           \
        {                                                                                   \
            /* Calling vectorized code when avx2 supported architecture detected */         \
            fla_dgetrf_small_avx2(m, n, buff_A, ldim_A, buff_p, info);                      \
        }                                                                                   \
        else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512) && *m < FLA_DGETRF_SMALL_AVX512_THRESH0 \
                && *n < FLA_DGETRF_SMALL_AVX512_THRESH0)                                    \
        {                                                                                   \
            if(fla_thread_get_num_threads() != 1 && *m < FLA_DGETRF_SMALL_AVX512_THRESH1    \
               && *n < FLA_DGETRF_SMALL_AVX512_THRESH1)                                     \
            {                                                                               \
                fla_dgetrf_small_avx512(m, n, buff_A, ldim_A, buff_p, info);                \
            }                                                                               \
            else                                                                            \
            {                                                                               \
                aocl_lapack_dgetrf2(m, n, buff_A, ldim_A, buff_p, info);                    \
            }                                                                               \
        }                                                                                   \
        else                                                                                \
        {                                                                                   \
            FLA_LU_piv_d_parallel(m, n, buff_A, ldim_A, buff_p, info);                      \
        }                                                                                   \
    }

#ifdef FLA_OPENMP_MULTITHREADING

#define LAPACK_getrf_body_z(prefix)                                    \
    if(*m <= FLA_ZGETRF_SMALL_THRESH && *n <= FLA_ZGETRF_SMALL_THRESH) \
    {                                                                  \
        FLA_LU_piv_z_var0(m, n, buff_A, ldim_A, buff_p, info);         \
    }                                                                  \
    else                                                               \
    {                                                                  \
        FLA_LU_piv_z_parallel(m, n, buff_A, ldim_A, buff_p, info);     \
    }

#define LAPACK_getrf_body_c(prefix)                                                         \
    aocl_int64_t i = 0;                                                                          \
    if(*m == 1 && *n == 1)                                                                  \
    {                                                                                       \
        FLA_LU_PIV_SMALL_C_1x1(i, *n, buff_A, ldim_A, buff_p, *info);                       \
    }                                                                                       \
    else if(*m == 2 && *n == 2)                                                             \
    {                                                                                       \
        FLA_LU_PIV_SMALL_C_2x2(i, *n, buff_A, ldim_A, buff_p, *info);                        \
    }                                                                                       \
    else if(*m == 3 && *n == 3)                                                             \
    {                                                                                       \
        FLA_LU_PIV_SMALL_C_3x3(i, *n, buff_A, ldim_A, buff_p, *info);                        \
    }                                                                                       \
    else if(*m == 4 && *n == 4)                                                             \
    {                                                                                       \
        FLA_LU_PIV_SMALL_C_4x4(i, *n, buff_A, ldim_A, buff_p, *info);                        \
    }                                                                                       \
    else if((*m <= FLA_CGETRF_SMALL_THRESH && *n <= FLA_CGETRF_SMALL_THRESH) || *m <= FLA_CGETRF_SMALL_THRESH0 \
|| *n <= FLA_CGETRF_SMALL_THRESH0)                  \
    {                                                                                       \
        FLA_LU_piv_c_var0(m, n, buff_A, ldim_A, buff_p, info);                              \
    }                                                                                       \
    else                                                                                    \
    {                                                                                       \
        FLA_LU_piv_c_parallel(m, n, buff_A, ldim_A, buff_p, info);                          \
    }

#else

#define LAPACK_getrf_body_z(prefix) FLA_LU_piv_z_var0(m, n, buff_A, ldim_A, buff_p, info);

#define LAPACK_getrf_body_c(prefix) FLA_LU_piv_c_var0(m, n, buff_A, ldim_A, buff_p, info);

#endif

#define LAPACK_getrf_body_s(prefix)                                                         \
    extern fla_context fla_global_context;                                                  \
    aocl_int64_t i = 0;                                                                          \
    if(*m == 2 && *n == 2)                                                                  \
    {                                                                                       \
        FLA_LU_PIV_SMALL_S_2x2(i, *n, buff_A, ldim_A, buff_p, *info);                       \
    }                                                                                       \
    else if(*m == 3 && *n == 3)                                                             \
    {                                                                                       \
        FLA_LU_PIV_SMALL_S_3x3(i, *n, buff_A, ldim_A, buff_p, *info);                       \
    }                                                                                       \
    else if(*m == 4 && *n == 4)                                                             \
    {                                                                                       \
        aocl_fla_init();                                                                    \
        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))                                               \
        {                                                                                   \
            fla_sgetrf_4x4_avx2(m, n, buff_A, ldim_A, buff_p, info);                        \
        }                                                                                   \
        else                                                                                \
        {                                                                                   \
            FLA_LU_PIV_SMALL_S_4x4(i, *n, buff_A, ldim_A, buff_p, *info);                   \
        }                                                                                   \
    }                                                                                       \
    else if(*m <= FLA_SGETRF_SMALL_THRESH0 && *n <= FLA_SGETRF_SMALL_THRESH0)               \
    {                                                                                       \
        FLA_LU_piv_small_s_var0(m, n, buff_A, ldim_A, buff_p, info);                        \
    }                                                                                       \
    else                                                                                    \
    {                                                                                       \
        /* Initialize global context data */                                                \
        aocl_fla_init();                                                                    \
        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2) && *m < FLA_SGETRF_SMALL_AVX2_THRESH0          \
           && *n < FLA_SGETRF_SMALL_AVX2_THRESH0)                                           \
        {                                                                                   \
            /* Calling vectorized code when avx2 supported architecture detected */         \
            fla_sgetrf_small_avx2(m, n, buff_A, ldim_A, buff_p, info);                      \
        }                                                                                   \
        else if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX512) && *m < FLA_SGETRF_SMALL_AVX512_THRESH0 \
                && *n < FLA_SGETRF_SMALL_AVX512_THRESH0)                                    \
        {                                                                                   \
            /* Calling vectorized code when avx512 supported architecture detected */       \
            fla_sgetrf_small_avx512(m, n, buff_A, ldim_A, buff_p, info);                    \
        }                                                                                   \
        else                                                                                \
        {                                                                                   \
            FLA_LU_piv_s_parallel(m, n, buff_A, ldim_A, buff_p, info);                      \
        }                                                                                   \
    }

#else /* FLA_ENABLE_AMD_OPT */

#define LAPACK_getrf_body_z LAPACK_getrf_body

#define LAPACK_getrf_body_c LAPACK_getrf_body

#define LAPACK_getrf_body_s LAPACK_getrf_body

/* Original FLA path */
/* The FLA_INT pivot object stores indices as fla_dim_t (64-bit). buff_p is the   */
/* LAPACK pivot array (aocl_int_t), which is 32-bit in the LP64 ABI. Attaching     */
/* buff_p directly would make FLA_LU_piv write 64-bit pivots into a 32-bit buffer  */
/* and overrun it. Use a fla_dim_t scratch buffer and narrow back to buff_p.       */
#define LAPACK_getrf_body_d(prefix)                                          \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                   \
    FLA_Obj A, p;                                                            \
    aocl_int64_t min_m_n = fla_min(*m, *n);                                  \
    FLA_Error e_val = FLA_SUCCESS;                                           \
    FLA_Error init_result;                                                   \
    fla_dim_t *buff_p64 = (fla_dim_t *)FLA_malloc(sizeof(fla_dim_t) * min_m_n); \
    /* FLA_malloc(0) returns NULL by design; only a nonzero pivot count    */ \
    /* is a genuine allocation that must succeed.                          */ \
    if(min_m_n > 0)                                                          \
        FLA_Check_error_code( FLA_Check_malloc_pointer( buff_p64 ) );        \
                                                                            \
    FLA_Init_safe(&init_result);                                            \
                                                                            \
    FLA_Obj_create_without_buffer(datatype, *m, *n, &A);                    \
    FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);                          \
                                                                            \
    FLA_Obj_create_without_buffer(FLA_INT, min_m_n, 1, &p);                 \
    FLA_Obj_attach_buffer(buff_p64, 1, min_m_n, &p);                        \
                                                                            \
    e_val = FLA_LU_piv(A, p);                                               \
    FLA_Shift_pivots_to(FLA_LAPACK_PIVOTS, p);                              \
                                                                            \
    FLA_Obj_free_without_buffer(&A);                                        \
    FLA_Obj_free_without_buffer(&p);                                        \
                                                                            \
    FLA_Finalize_safe(init_result);                                        \
    for(aocl_int64_t ip = 0; ip < min_m_n; ip++)                            \
        buff_p[ip] = (aocl_int_t)buff_p64[ip];                             \
    FLA_free(buff_p64);                                                     \
                                                                            \
    if(e_val != FLA_SUCCESS)                                                \
        *info = e_val + 1;

#endif /* FLA_ENABLE_AMD_OPT */

// Note that p should be set zero.
#define LAPACK_getrf_body(prefix)                                                            \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                                   \
    FLA_Obj A, p;                                                                            \
    aocl_int64_t min_m_n = fla_min(*m, *n);                                                       \
    FLA_Error e_val = FLA_SUCCESS;                                                           \
    FLA_Error init_result;                                                                   \
    FLA_Bool skip = FALSE;                                                                   \
                                                                                             \
    if(*m < FLA_GETRF_SMALL && *n < FLA_GETRF_SMALL) /* Small sizes- lapack path */          \
    {                                                                                        \
        switch(datatype)                                                                     \
        {                                                                                    \
            case FLA_FLOAT:                                                                  \
            {                                                                                \
                lapack_sgetrf(m, n, (float *)buff_A, ldim_A, buff_p, info);                  \
                break;                                                                       \
            }                                                                                \
            case FLA_COMPLEX:                                                                \
            {                                                                                \
                lapack_cgetrf(m, n, (scomplex *)buff_A, ldim_A, buff_p, info);               \
                break;                                                                       \
            }                                                                                \
            case FLA_DOUBLE_COMPLEX:                                                         \
            {                                                                                \
                lapack_zgetrf(m, n, (dcomplex *)buff_A, ldim_A, buff_p, info);               \
                break;                                                                       \
            }                                                                                \
        }                                                                                    \
        if(*info != 0)                                                                       \
            skip = TRUE;                                                                     \
    }                                                                                        \
    else if((datatype == FLA_FLOAT && *m < FLA_GETRF_FLOAT && *n < FLA_GETRF_FLOAT)          \
            || (datatype == FLA_COMPLEX && *m < FLA_GETRF_COMPLEX && *n < FLA_GETRF_COMPLEX) \
            || (datatype == FLA_DOUBLE_COMPLEX && *m < FLA_GETRF_DOUBLE_COMPLEX              \
                && *n < FLA_GETRF_DOUBLE_COMPLEX))                                           \
    {                                                                                        \
        fla_dim_t *buff_p64 = (fla_dim_t *) FLA_malloc(sizeof(fla_dim_t) * min_m_n);         \
        /* FLA_malloc(0) returns NULL by design when min_m_n == 0; skip the check then. */   \
        if(min_m_n > 0)                                                                       \
            FLA_Check_error_code( FLA_Check_malloc_pointer( buff_p64 ) );                     \
        FLA_Init_safe(&init_result);                                                         \
                                                                                             \
        FLA_Obj_create_without_buffer(datatype, *m, *n, &A);                                 \
        FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);                                       \
                                                                                             \
        FLA_Obj_create_without_buffer(FLA_INT, min_m_n, 1, &p);                              \
        FLA_Obj_attach_buffer(buff_p64, 1, min_m_n, &p);                                     \
                                                                                             \
        e_val = FLA_LU_piv(A, p);                                                            \
        FLA_Shift_pivots_to(FLA_LAPACK_PIVOTS, p);                                           \
                                                                                             \
        FLA_Obj_free_without_buffer(&A);                                                     \
        FLA_Obj_free_without_buffer(&p);                                                     \
                                                                                             \
        FLA_Finalize_safe(init_result);                                                      \
        for(aocl_int64_t ip = 0; ip < min_m_n; ip++)                                         \
            buff_p[ip] = (aocl_int_t)buff_p64[ip];                                           \
        FLA_free(buff_p64);                                                                  \
    }                                                                                        \
    else                                                                                     \
    {                                                                                        \
        switch(datatype)                                                                     \
        {                                                                                    \
            case FLA_FLOAT:                                                                  \
            {                                                                                \
                aocl_lapack_sgetrf2(m, n, (float *)buff_A, ldim_A, buff_p, info);                       \
                break;                                                                       \
            }                                                                                \
            case FLA_COMPLEX:                                                                \
            {                                                                                \
                aocl_lapack_cgetrf2(m, n, (scomplex *)buff_A, ldim_A, buff_p, info);                    \
                break;                                                                       \
            }                                                                                \
            case FLA_DOUBLE_COMPLEX:                                                         \
            {                                                                                \
                aocl_lapack_zgetrf2(m, n, (dcomplex *)buff_A, ldim_A, buff_p, info);                    \
                break;                                                                       \
            }                                                                                \
        }                                                                                    \
        if(*info != 0)                                                                       \
            skip = TRUE;                                                                     \
    }                                                                                        \
                                                                                             \
    if(e_val != FLA_SUCCESS)                                                                 \
        *info = e_val + 1;                                                                   \
    else if(skip != TRUE)                                                                    \
        *info = 0;

#else /* FLA_ENABLE_SUPERMATRIX */

#define LAPACK_getrf_body_s LAPACK_getrf_body
#define LAPACK_getrf_body_d LAPACK_getrf_body
#define LAPACK_getrf_body_c LAPACK_getrf_body
#define LAPACK_getrf_body_z LAPACK_getrf_body

// Note that p should be set zero.
/* The FLA_INT pivot object stores indices as fla_dim_t (64-bit). buff_p is the   */
/* LAPACK pivot array (aocl_int_t), which is 32-bit in the LP64 ABI. Use a         */
/* fla_dim_t scratch buffer for the FLAME/FLASH pivots and narrow back to buff_p.  */
#define LAPACK_getrf_body(prefix)                            \
    FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);   \
    FLA_Obj A, p, AH, ph;                                    \
    aocl_int64_t min_m_n = fla_min(*m, *n);                  \
    fla_dim_t nth, b_flash;                                  \
    FLA_Error e_val;                                         \
    FLA_Error init_result;                                   \
    fla_dim_t *buff_p64 = (fla_dim_t *)FLA_malloc(sizeof(fla_dim_t) * min_m_n); \
    /* FLA_malloc(0) returns NULL by design when min_m_n == 0; skip check. */ \
    if(min_m_n > 0)                                          \
        FLA_Check_error_code( FLA_Check_malloc_pointer( buff_p64 ) ); \
                                                             \
    nth = FLASH_get_num_threads(1);                          \
    b_flash = FLA_EXT_HIER_BLOCKSIZE;                        \
                                                             \
    FLA_Init_safe(&init_result);                             \
    FLASH_Queue_set_num_threads(nth);                        \
                                                             \
    FLA_Obj_create_without_buffer(datatype, *m, *n, &A);     \
    FLA_Obj_attach_buffer(buff_A, 1, *ldim_A, &A);           \
                                                             \
    FLA_Obj_create_without_buffer(FLA_INT, min_m_n, 1, &p);  \
    FLA_Obj_attach_buffer(buff_p64, 1, min_m_n, &p);         \
                                                             \
    FLA_Set(FLA_ZERO, p);                                    \
                                                             \
    FLASH_Obj_create_hier_copy_of_flat(A, 1, &b_flash, &AH); \
    FLASH_Obj_create_hier_copy_of_flat(p, 1, &b_flash, &ph); \
                                                             \
    e_val = FLASH_LU_piv(AH, ph);                            \
                                                             \
    FLASH_Obj_flatten(AH, A);                                \
    FLASH_Obj_flatten(ph, p);                                \
    FLA_Shift_pivots_to(FLA_LAPACK_PIVOTS, p);               \
                                                             \
    FLA_Obj_free(&AH);                                       \
    FLA_Obj_free(&ph);                                       \
                                                             \
    FLA_Obj_free_without_buffer(&A);                         \
    FLA_Obj_free_without_buffer(&p);                         \
                                                             \
    FLA_Finalize_safe(init_result);                          \
    for(aocl_int64_t ip = 0; ip < min_m_n; ip++)             \
        buff_p[ip] = (aocl_int_t)buff_p64[ip];               \
    FLA_free(buff_p64);                                      \
                                                             \
    if(e_val != FLA_SUCCESS)                                 \
        *info = e_val + 1;                                   \
    else                                                     \
        *info = 0;

#endif /* FLA_ENABLE_SUPERMATRIX */

LAPACK_getrf(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgetrf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);

    {
        LAPACK_RETURN_CHECK_VAR1(sgetrf_check(m, n, buff_A, ldim_A, buff_p, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_getrf_body_s(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_getrf(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgetrf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dgetrf_check(m, n, buff_A, ldim_A, buff_p, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_getrf_body_d(d)
            /* fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_getrf(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgetrf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cgetrf_check(m, n, buff_A, ldim_A, buff_p, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        extern fla_context fla_global_context;
        aocl_fla_init();
        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
        {
            LAPACK_getrf_body_c(c)
        }
        else
        {
            LAPACK_getrf_body(c)
        }

        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }

    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_getrf(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgetrf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zgetrf_check(m, n, buff_A, ldim_A, buff_p, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        extern fla_context fla_global_context;
        aocl_fla_init();
        if(FLA_IS_MIN_ARCH_ID(FLA_ARCH_AVX2))
        {
            LAPACK_getrf_body_z(z)
        }
        else
        {
            LAPACK_getrf_body(z)
        }

        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#define LAPACK_getf2(prefix)                                                                 \
    void aocl_lapack_##prefix##getf2(aocl_int64_t *m, aocl_int64_t *n, PREFIX2LAPACK_TYPEDEF(prefix) * buff_A, \
                             aocl_int64_t * ldim_A, aocl_int_t * buff_p, aocl_int64_t * info)

LAPACK_getf2(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgetf2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(sgetf2_check(m, n, buff_A, ldim_A, buff_p, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_getrf_body(s)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_getf2(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgetf2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dgetf2_check(m, n, buff_A, ldim_A, buff_p, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_getrf_body_d(d)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_getf2(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgetf2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cgetf2_check(m, n, buff_A, ldim_A, buff_p, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_getrf_body(c)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}
LAPACK_getf2(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgetf2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n,
                      *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zgetf2_check(m, n, buff_A, ldim_A, buff_p, info), fla_error)
    }
    if(fla_error == LAPACK_SUCCESS)
    {
        LAPACK_getrf_body(z)
            /** fla_error set to 0 on LAPACK_SUCCESS */
            fla_error
            = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return;
}

#endif
